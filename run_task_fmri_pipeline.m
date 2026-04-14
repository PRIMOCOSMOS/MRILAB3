function run_task_fmri_pipeline()
% RUN_TASK_FMRI_PIPELINE
% -------------------------------------------------------------------------
% 基于 DPABI/SPM 思路实现任务态 fMRI-BOLD 完整流程：
% 1) 原始数据(DICOM/NIfTI)准备
% 2) 预处理：去前若干时间点、切片时序、重对齐、配准、分割+DARTEL、标准化、平滑
% 3) 一级统计：模型指定、估计、对比
% 4) 结果：阈值化统计图 + MNI 坐标标注 + 立体视图导出
%
% 运行前请先修改配置文件：
%   cfg = task_fmri_pipeline_config();
%
% 依赖：
%   - MATLAB
%   - SPM12 (已加入 MATLAB path)
%
% 数据组织参考 DPABI 风格（可在配置中微调）：
%   DataRaw/
%     Sub001/
%       anat/ 或 T1Img/
%       func/ 或 FunImg/
% -------------------------------------------------------------------------

    cfg = task_fmri_pipeline_config();

    % --------------------------- 环境检查 -------------------------------
    assert(exist('spm', 'file') == 2, 'SPM 未在 MATLAB 路径中，请先 addpath(spm12).');
    spm('defaults', 'fmri');
    spm_jobman('initcfg');

    subjects = discover_subjects(cfg.paths.dataRawDir);
    assert(~isempty(subjects), '未在 %s 中发现受试者目录。', cfg.paths.dataRawDir);

    % --------------------------- 逐被试处理 -----------------------------
    preprocInfo = struct([]);
    for s = 1:numel(subjects)
        subjID = subjects{s};
        fprintf('\n[%s] ===== 处理开始: %s =====\n', datestr(now, 31), subjID);
        preprocInfo(s) = run_preprocessing_for_subject(subjID, cfg); %#ok<AGROW>
    end

    % ------------------------- DARTEL 群体模板 ---------------------------
    run_group_dartel(preprocInfo, cfg);

    % ------------------ 标准化+平滑 + 一级统计 + 结果可视化 ---------------
    for s = 1:numel(preprocInfo)
        subjID = preprocInfo(s).subjID;
        fprintf('\n[%s] ===== 一级分析开始: %s =====\n', datestr(now, 31), subjID);
        normalize_and_smooth_subject(preprocInfo(s), cfg);
        run_first_level_glm(preprocInfo(s), cfg);
        render_subject_results(preprocInfo(s), cfg);
    end

    fprintf('\n[%s] 全流程完成。\n', datestr(now, 31));
end

% ============================== 功能块 1 ================================
% 被试发现
function subjects = discover_subjects(dataRawDir)
    d = dir(dataRawDir);
    d = d([d.isdir]);
    names = {d.name};
    names = names(~startsWith(names, '.'));
    subjects = sort(names);
end

% ============================== 功能块 2 ================================
% 单被试预处理（至分割输出）
function info = run_preprocessing_for_subject(subjID, cfg)
    subjRawDir = fullfile(cfg.paths.dataRawDir, subjID);
    subjWorkDir = fullfile(cfg.paths.derivativeDir, subjID);
    ensure_dir(subjWorkDir);

    % ---------- 参数局部化：目录候选 ----------
    anatDir = locate_first_existing(subjRawDir, cfg.paths.anatDirCandidates);
    funcDir = locate_first_existing(subjRawDir, cfg.paths.funcDirCandidates);
    assert(~isempty(anatDir), '受试者 %s 未找到解剖目录。', subjID);
    assert(~isempty(funcDir), '受试者 %s 未找到功能目录。', subjID);

    % ---------- 原始数据转换：DICOM -> NIfTI ----------
    anatNii = ensure_single_nifti(anatDir, fullfile(subjWorkDir, 'anat_nifti'));
    runNii4D = ensure_run_niftis(funcDir, fullfile(subjWorkDir, 'func_nifti'));
    assert(~isempty(runNii4D), '受试者 %s 未找到功能 NIfTI。', subjID);

    % ---------- 去前若干时间点 ----------
    trimmed4D = cell(size(runNii4D));
    for r = 1:numel(runNii4D)
        outDir = fullfile(subjWorkDir, sprintf('run%02d_trim', r));
        ensure_dir(outDir);
        trimmed4D{r} = remove_first_volumes(runNii4D{r}, cfg.preproc.removeFirstN, outDir);
    end

    % ---------- 分离4D为3D(供SPM批处理) ----------
    runScans = cell(size(trimmed4D));
    for r = 1:numel(trimmed4D)
        runDir = fullfile(subjWorkDir, sprintf('run%02d_3d', r));
        ensure_dir(runDir);
        copyfile(trimmed4D{r}, fullfile(runDir, 'trimmed.nii'));
        spm_file_split(fullfile(runDir, 'trimmed.nii'));
        runScans{r} = cellstr(spm_select('ExtFPList', runDir, '^trimmed_.*\.nii$', Inf));
        assert(~isempty(runScans{r}), '受试者 %s run %d 分离3D失败。', subjID, r);
    end

    % ---------- 切片时序 ----------
    runScans = run_slice_timing(runScans, cfg);

    % ---------- 重对齐（含重采样） ----------
    [runScans, meanFuncPath, motionParamFiles] = run_realign_estwrite(runScans, cfg);

    % ---------- T1 到功能均值图配准 ----------
    coregAnat = run_coregister_estwrite(anatNii, meanFuncPath, cfg);

    % ---------- New Segment ----------
    [rc1, rc2, flowField] = run_new_segment(coregAnat, cfg);

    info.subjID = subjID;
    info.subjWorkDir = subjWorkDir;
    info.runScansRealigned = runScans;
    info.motionParamFiles = motionParamFiles;
    info.coregAnat = coregAnat;
    info.rc1 = rc1;
    info.rc2 = rc2;
    info.flowField = flowField;
end

% ============================== 功能块 3 ================================
% DARTEL 群体模板构建
function run_group_dartel(preprocInfo, cfg)
    rc1_all = cell(numel(preprocInfo), 1);
    rc2_all = cell(numel(preprocInfo), 1);
    for i = 1:numel(preprocInfo)
        rc1_all{i} = preprocInfo(i).rc1;
        rc2_all{i} = preprocInfo(i).rc2;
    end

    ensure_dir(cfg.paths.dartelDir);
    prev = pwd;
    cleanupObj = onCleanup(@() cd(prev)); %#ok<NASGU>
    cd(cfg.paths.dartelDir);

    matlabbatch = [];
    matlabbatch{1}.spm.tools.dartel.warp.images = {rc1_all, rc2_all};
    matlabbatch{1}.spm.tools.dartel.warp.settings.template = cfg.preproc.dartelTemplateName;
    matlabbatch{1}.spm.tools.dartel.warp.settings.rform = 0;
    matlabbatch{1}.spm.tools.dartel.warp.settings.param = cfg.preproc.dartelParam;
    matlabbatch{1}.spm.tools.dartel.warp.settings.optim = cfg.preproc.dartelOptim;
    spm_jobman('run', matlabbatch);
end

% ============================== 功能块 4 ================================
% DARTEL 标准化 + 平滑
function normalize_and_smooth_subject(info, cfg)
    subjDir = info.subjWorkDir;
    normDir = fullfile(subjDir, 'normalized');
    smoothDir = fullfile(subjDir, 'smoothed');
    ensure_dir(normDir);
    ensure_dir(smoothDir);

    allRealigned = vertcat(info.runScansRealigned{:});
    allRealigned = to_scan_refs(allRealigned);

    matlabbatch = [];

    % ---------- DARTEL 到 MNI ----------
    matlabbatch{1}.spm.tools.dartel.mni_norm.template = {fullfile(cfg.paths.dartelDir, cfg.preproc.dartelLastTemplate)};
    matlabbatch{1}.spm.tools.dartel.mni_norm.data.subj.flowfield = {info.flowField};
    matlabbatch{1}.spm.tools.dartel.mni_norm.data.subj.images = {allRealigned};
    matlabbatch{1}.spm.tools.dartel.mni_norm.vox = cfg.preproc.normalizeVoxelSize;
    matlabbatch{1}.spm.tools.dartel.mni_norm.bb = cfg.preproc.normalizeBoundingBox;
    matlabbatch{1}.spm.tools.dartel.mni_norm.preserve = 0;
    matlabbatch{1}.spm.tools.dartel.mni_norm.fwhm = [0 0 0];

    spm_jobman('run', matlabbatch);

    % 收集标准化输出
    wImgs = cellstr(spm_select('FPListRec', subjDir, '^wra.*\.nii$'));
    assert(~isempty(wImgs), '受试者 %s 未找到标准化图像。', info.subjID);
    for i = 1:numel(wImgs)
        copyfile(strtrim(wImgs{i}), fullfile(normDir, [num2str(i, '%05d') '.nii']));
    end

    % ---------- 空间平滑 ----------
    normImgs = cellstr(spm_select('FPList', normDir, '^\d+\.nii$'));
    normImgs = to_scan_refs(normImgs);

    matlabbatch = [];
    matlabbatch{1}.spm.spatial.smooth.data = normImgs;
    matlabbatch{1}.spm.spatial.smooth.fwhm = cfg.preproc.smoothFWHM;
    matlabbatch{1}.spm.spatial.smooth.dtype = 0;
    matlabbatch{1}.spm.spatial.smooth.im = 0;
    matlabbatch{1}.spm.spatial.smooth.prefix = 's';
    spm_jobman('run', matlabbatch);
end

% ============================== 功能块 5 ================================
% 一级 GLM：Specify -> Estimate -> Contrast
function run_first_level_glm(info, cfg)
    firstLevelDir = fullfile(info.subjWorkDir, 'first_level');
    ensure_dir(firstLevelDir);

    condFile = fullfile(cfg.paths.onsetDir, info.subjID, cfg.firstLevel.conditionFileName);
    assert(exist(condFile, 'file') == 2, ...
        '未找到条件文件: %s (请参照 README 提供 names/onsets/durations)', condFile);
    S = load(condFile);
    assert(isfield(S, 'names') && isfield(S, 'onsets') && isfield(S, 'durations'), ...
        '条件文件必须包含 names/onsets/durations.');

    scans = cellstr(spm_select('FPListRec', fullfile(info.subjWorkDir, 'smoothed'), '^s\d+\.nii$'));
    scans = to_scan_refs(scans);
    assert(~isempty(scans), '受试者 %s 未找到平滑后功能图像。', info.subjID);

    % -------- 参数局部化：默认单 session；如多 run 可拓展为逐 run 指定 -------
    matlabbatch = [];
    matlabbatch{1}.spm.stats.fmri_spec.dir = {firstLevelDir};
    matlabbatch{1}.spm.stats.fmri_spec.timing.units = cfg.firstLevel.units;
    matlabbatch{1}.spm.stats.fmri_spec.timing.RT = cfg.firstLevel.TR;
    matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t = cfg.firstLevel.fmri_t;
    matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t0 = cfg.firstLevel.fmri_t0;
    matlabbatch{1}.spm.stats.fmri_spec.sess.scans = scans;

    for c = 1:numel(S.names)
        matlabbatch{1}.spm.stats.fmri_spec.sess.cond(c).name = S.names{c};
        matlabbatch{1}.spm.stats.fmri_spec.sess.cond(c).onset = S.onsets{c};
        matlabbatch{1}.spm.stats.fmri_spec.sess.cond(c).duration = S.durations{c};
        matlabbatch{1}.spm.stats.fmri_spec.sess.cond(c).tmod = 0;
        matlabbatch{1}.spm.stats.fmri_spec.sess.cond(c).pmod = struct('name', {}, 'param', {}, 'poly', {});
        matlabbatch{1}.spm.stats.fmri_spec.sess.cond(c).orth = 1;
    end

    if cfg.firstLevel.includeMotionRegressors
        rp = info.motionParamFiles;
        rp = rp(cellfun(@(x) exist(x, 'file') == 2, rp));
        if ~isempty(rp)
            % 多 run 时建议按 session 分配；此处合并附加
            allRP = fullfile(firstLevelDir, 'motion_regressors.txt');
            concat_motion_params(rp, allRP);
            matlabbatch{1}.spm.stats.fmri_spec.sess.multi_reg = {allRP};
        else
            matlabbatch{1}.spm.stats.fmri_spec.sess.multi_reg = {''};
        end
    else
        matlabbatch{1}.spm.stats.fmri_spec.sess.multi_reg = {''};
    end

    matlabbatch{1}.spm.stats.fmri_spec.sess.hpf = cfg.firstLevel.hpf;
    matlabbatch{1}.spm.stats.fmri_spec.fact = struct('name', {}, 'levels', {});
    matlabbatch{1}.spm.stats.fmri_spec.bases.hrf.derivs = cfg.firstLevel.hrfDerivs;
    matlabbatch{1}.spm.stats.fmri_spec.volt = 1;
    matlabbatch{1}.spm.stats.fmri_spec.global = 'None';
    matlabbatch{1}.spm.stats.fmri_spec.mthresh = 0.8;
    matlabbatch{1}.spm.stats.fmri_spec.mask = {''};
    matlabbatch{1}.spm.stats.fmri_spec.cvi = cfg.firstLevel.cvi;
    spm_jobman('run', matlabbatch);

    matlabbatch = [];
    matlabbatch{1}.spm.stats.fmri_est.spmmat = {fullfile(firstLevelDir, 'SPM.mat')};
    matlabbatch{1}.spm.stats.fmri_est.method.Classical = 1;
    spm_jobman('run', matlabbatch);

    matlabbatch = [];
    matlabbatch{1}.spm.stats.con.spmmat = {fullfile(firstLevelDir, 'SPM.mat')};
    matlabbatch{1}.spm.stats.con.delete = 1;
    for i = 1:numel(cfg.firstLevel.contrasts)
        matlabbatch{1}.spm.stats.con.consess{i}.tcon.name = cfg.firstLevel.contrasts(i).name;
        matlabbatch{1}.spm.stats.con.consess{i}.tcon.weights = cfg.firstLevel.contrasts(i).weights;
        matlabbatch{1}.spm.stats.con.consess{i}.tcon.sessrep = 'none';
    end
    spm_jobman('run', matlabbatch);
end

% ============================== 功能块 6 ================================
% 统计阈值与立体视图标注激活区域
function render_subject_results(info, cfg)
    firstLevelDir = fullfile(info.subjWorkDir, 'first_level');
    resultDir = fullfile(firstLevelDir, 'results');
    ensure_dir(resultDir);

    xSPM = struct();
    xSPM.swd = firstLevelDir;
    xSPM.Ic = cfg.results.contrastIndex;
    xSPM.u = cfg.results.voxelPThreshold;
    xSPM.k = cfg.results.clusterExtent;
    xSPM.thresDesc = cfg.results.thresDesc;
    xSPM.title = cfg.results.title;
    xSPM.Im = [];
    xSPM.pm = [];
    xSPM.Ex = [];
    xSPM.n = 1;

    [xSPM, SPM] = spm_getSPM(xSPM); %#ok<ASGLU>
    hReg = spm_results_ui('Setup', xSPM);
    spm_list('List', xSPM, hReg);

    % 导出阈值结果表（含显著峰 MNI 坐标）
    [TabDat, xSPM, hReg] = spm_list('Table', xSPM, hReg); %#ok<NASGU,ASGLU>
    save(fullfile(resultDir, 'xSPM_thresholded.mat'), 'xSPM', 'TabDat');

    % 立体视图（Orthviews）并保存快照
    spm_orthviews('Reset');
    bg = cfg.results.backgroundImage;
    if exist(bg, 'file') ~= 2
        bg = fullfile(spm('Dir'), 'canonical', 'single_subj_T1.nii');
    end
    spm_orthviews('Image', bg, [0.05 0.05 0.9 0.9]);
    spm_orthviews('AddColouredBlobs', 1, xSPM.XYZ, xSPM.Z, xSPM.M, [1 0 0]);
    spm_orthviews('Redraw');

    fig = spm_figure('FindWin', 'Graphics');
    if ~isempty(fig)
        exportgraphics(fig, fullfile(resultDir, 'orthview_activation.png'), 'Resolution', 300);
    end

    % 3D 渲染图（可选）
    try
        spm_render(xSPM.Z, xSPM.XYZ, fullfile(spm('Dir'), 'rend', 'render_single_subj.mat'));
        fig = spm_figure('FindWin', 'Render');
        if ~isempty(fig)
            exportgraphics(fig, fullfile(resultDir, 'render_activation.png'), 'Resolution', 300);
        end
    catch
        warning('3D 渲染失败，已保留正交视图与坐标表。');
    end
end

% ============================== 工具函数 ================================
function outPath = remove_first_volumes(inNii, nDrop, outDir)
    V = spm_vol(inNii);
    assert(numel(V) > nDrop, '文件 %s 的体积数不足以去除前 %d 个时间点。', inNii, nDrop);

    Y = spm_read_vols(V);
    Y = Y(:, :, :, (nDrop + 1):end);

    outPath = fullfile(outDir, 'trimmed.nii');
    Vo = V(1);
    Vo.fname = outPath;
    Vo.n = [1 1];
    spm_write_vol(Vo, Y(:, :, :, 1));
    for t = 2:size(Y, 4)
        Vo.n = [t 1];
        spm_write_vol(Vo, Y(:, :, :, t));
    end
end

function outScans = run_slice_timing(runScans, cfg)
    outScans = runScans;
    for r = 1:numel(runScans)
        matlabbatch = [];
        matlabbatch{1}.spm.temporal.st.scans = {runScans{r}};
        matlabbatch{1}.spm.temporal.st.nslices = cfg.preproc.nslices;
        matlabbatch{1}.spm.temporal.st.tr = cfg.preproc.TR;
        matlabbatch{1}.spm.temporal.st.ta = cfg.preproc.TR - cfg.preproc.TR / cfg.preproc.nslices;
        matlabbatch{1}.spm.temporal.st.so = cfg.preproc.sliceOrder;
        matlabbatch{1}.spm.temporal.st.refslice = cfg.preproc.refSlice;
        matlabbatch{1}.spm.temporal.st.prefix = cfg.preproc.sliceTimingPrefix;
        spm_jobman('run', matlabbatch);

        p = fileparts(runScans{r}{1});
        outScans{r} = cellstr(spm_select('ExtFPList', p, ...
            ['^' regexptranslate('escape', cfg.preproc.sliceTimingPrefix) 'trimmed_.*\.nii$'], Inf));
    end
end

function [outScans, meanFuncPath, motionParamFiles] = run_realign_estwrite(runScans, cfg)
    allScans = vertcat(runScans{:});
    allScans = to_scan_refs(allScans);

    matlabbatch = [];
    matlabbatch{1}.spm.spatial.realign.estwrite.data = {allScans};
    matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.quality = cfg.preproc.realign.quality;
    matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.sep = cfg.preproc.realign.sep;
    matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.fwhm = cfg.preproc.realign.fwhm;
    matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.rtm = cfg.preproc.realign.rtm;
    matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.interp = cfg.preproc.realign.interp;
    matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.wrap = cfg.preproc.realign.wrap;
    matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.weight = '';
    matlabbatch{1}.spm.spatial.realign.estwrite.roptions.which = cfg.preproc.realign.which;
    matlabbatch{1}.spm.spatial.realign.estwrite.roptions.interp = cfg.preproc.realign.writeInterp;
    matlabbatch{1}.spm.spatial.realign.estwrite.roptions.wrap = cfg.preproc.realign.wrap;
    matlabbatch{1}.spm.spatial.realign.estwrite.roptions.mask = 1;
    matlabbatch{1}.spm.spatial.realign.estwrite.roptions.prefix = cfg.preproc.realignPrefix;
    spm_jobman('run', matlabbatch);

    % 输出收集
    outScans = runScans;
    motionParamFiles = cell(numel(runScans), 1);
    for r = 1:numel(runScans)
        p = fileparts(runScans{r}{1});
        outScans{r} = cellstr(spm_select('ExtFPList', p, ...
            ['^' regexptranslate('escape', cfg.preproc.realignPrefix) ...
            regexptranslate('escape', cfg.preproc.sliceTimingPrefix) 'trimmed_.*\.nii$'], Inf));
        rp = cellstr(spm_select('FPList', p, '^rp_.*\.txt$'));
        if ~isempty(rp)
            motionParamFiles{r} = strtrim(rp{1});
        else
            motionParamFiles{r} = '';
        end
    end

    firstRunDir = fileparts(runScans{1}{1});
    meanImg = cellstr(spm_select('FPList', firstRunDir, '^mean.*\.nii$'));
    assert(~isempty(meanImg), '未找到重对齐均值图像 mean*.nii。');
    meanFuncPath = strtrim(meanImg{1});
end

function coregAnat = run_coregister_estwrite(anatNii, meanFuncPath, cfg)
    matlabbatch = [];
    matlabbatch{1}.spm.spatial.coreg.estwrite.ref = {to_scan_ref(meanFuncPath)};
    matlabbatch{1}.spm.spatial.coreg.estwrite.source = {to_scan_ref(anatNii)};
    matlabbatch{1}.spm.spatial.coreg.estwrite.other = {''};
    matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.cost_fun = cfg.preproc.coreg.cost_fun;
    matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.sep = cfg.preproc.coreg.sep;
    matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.tol = cfg.preproc.coreg.tol;
    matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.fwhm = cfg.preproc.coreg.fwhm;
    matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.interp = cfg.preproc.coreg.interp;
    matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.wrap = [0 0 0];
    matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.mask = 0;
    matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.prefix = cfg.preproc.coregPrefix;
    spm_jobman('run', matlabbatch);

    [p, n, e] = fileparts(anatNii);
    coregAnat = fullfile(p, [cfg.preproc.coregPrefix n e]);
    assert(exist(coregAnat, 'file') == 2, '配准输出不存在：%s', coregAnat);
end

function [rc1, rc2, flowField] = run_new_segment(coregAnat, cfg)
    tpm = fullfile(spm('Dir'), 'tpm', 'TPM.nii');
    assert(exist(tpm, 'file') == 2, '未找到 TPM 文件：%s', tpm);

    matlabbatch = [];
    matlabbatch{1}.spm.spatial.preproc.channel.vols = {to_scan_ref(coregAnat)};
    matlabbatch{1}.spm.spatial.preproc.channel.biasreg = cfg.preproc.segment.biasreg;
    matlabbatch{1}.spm.spatial.preproc.channel.biasfwhm = cfg.preproc.segment.biasfwhm;
    matlabbatch{1}.spm.spatial.preproc.channel.write = [0 1];

    for k = 1:6
        matlabbatch{1}.spm.spatial.preproc.tissue(k).tpm = {[tpm ',' num2str(k)]};
        matlabbatch{1}.spm.spatial.preproc.tissue(k).ngaus = cfg.preproc.segment.ngaus(k);
        matlabbatch{1}.spm.spatial.preproc.tissue(k).native = cfg.preproc.segment.native{k};
        matlabbatch{1}.spm.spatial.preproc.tissue(k).warped = [0 0];
    end

    matlabbatch{1}.spm.spatial.preproc.warp.mrf = cfg.preproc.segment.warp.mrf;
    matlabbatch{1}.spm.spatial.preproc.warp.cleanup = cfg.preproc.segment.warp.cleanup;
    matlabbatch{1}.spm.spatial.preproc.warp.reg = cfg.preproc.segment.warp.reg;
    matlabbatch{1}.spm.spatial.preproc.warp.affreg = cfg.preproc.segment.warp.affreg;
    matlabbatch{1}.spm.spatial.preproc.warp.fwhm = cfg.preproc.segment.warp.fwhm;
    matlabbatch{1}.spm.spatial.preproc.warp.samp = cfg.preproc.segment.warp.samp;
    matlabbatch{1}.spm.spatial.preproc.warp.write = [1 1];
    spm_jobman('run', matlabbatch);

    [p, n, e] = fileparts(coregAnat);
    rc1 = fullfile(p, ['rc1' n e]);
    rc2 = fullfile(p, ['rc2' n e]);
    flowField = fullfile(p, ['u_rc1' n e]);
    assert(exist(rc1, 'file') == 2 && exist(rc2, 'file') == 2 && exist(flowField, 'file') == 2, ...
        '分割输出不完整，请检查 New Segment 结果。');
end

function concat_motion_params(inFiles, outFile)
    allM = [];
    for i = 1:numel(inFiles)
        if exist(inFiles{i}, 'file') ~= 2
            continue;
        end
        m = load(inFiles{i});
        allM = [allM; m]; %#ok<AGROW>
    end
    dlmwrite(outFile, allM, 'delimiter', '\t');
end

function pathOut = locate_first_existing(baseDir, candidates)
    pathOut = '';
    for i = 1:numel(candidates)
        p = fullfile(baseDir, candidates{i});
        if exist(p, 'dir') == 7
            pathOut = p;
            return;
        end
    end
end

function nii = ensure_single_nifti(inDir, outDir)
    ensure_dir(outDir);
    niiList = cellstr(spm_select('FPList', inDir, '.*\.nii$'));
    if isempty(niiList)
        nii = convert_dicom_dir_to_nifti(inDir, outDir, true);
    else
        nii = strtrim(niiList{1});
    end
end

function runs = ensure_run_niftis(funcDir, outDir)
    ensure_dir(outDir);
    runs = {};

    runFolders = dir(funcDir);
    runFolders = runFolders([runFolders.isdir]);
    runFolders = runFolders(~startsWith({runFolders.name}, '.'));

    % 情况 A：func 目录内直接是 NIfTI（单 run）
    topNii = cellstr(spm_select('FPList', funcDir, '.*\.nii$'));
    if ~isempty(topNii)
        runs{1} = strtrim(topNii{1});
        return;
    end

    % 情况 B：每个子目录一个 run（DICOM 或 NIfTI）
    for r = 1:numel(runFolders)
        rDir = fullfile(funcDir, runFolders(r).name);
        rNii = cellstr(spm_select('FPList', rDir, '.*\.nii$'));
        if ~isempty(rNii)
            runs{end + 1} = strtrim(rNii{1}); %#ok<AGROW>
        else
            runOut = fullfile(outDir, sprintf('run%02d', r));
            ensure_dir(runOut);
            runs{end + 1} = convert_dicom_dir_to_nifti(rDir, runOut, false); %#ok<AGROW>
        end
    end
end

function nii = convert_dicom_dir_to_nifti(dicomDir, outDir, singleOnly)
    dcm = cellstr(spm_select('FPList', dicomDir, '.*\.(dcm|IMA)$'));
    if isempty(dcm)
        dcm = cellstr(spm_select('FPList', dicomDir, '.*$'));
        dcm = dcm(cellfun(@(x) exist(strtrim(x), 'file') == 2, dcm));
    end
    assert(~isempty(dcm), '目录 %s 中未发现可转换 DICOM。', dicomDir);

    hdr = spm_dicom_headers(char(dcm));
    prev = pwd;
    cleanupObj = onCleanup(@() cd(prev)); %#ok<NASGU>
    cd(outDir);
    spm_dicom_convert(hdr, 'all', 'flat', 'nii');
    niiFiles = cellstr(spm_select('FPList', outDir, '.*\.nii$'));
    assert(~isempty(niiFiles), 'DICOM 转 NIfTI 失败：%s', dicomDir);

    if singleOnly
        nii = strtrim(niiFiles{1});
    else
        % 如果输出多个3D文件，则合并为4D
        if numel(niiFiles) == 1
            nii = strtrim(niiFiles{1});
        else
            nii = fullfile(outDir, 'merged_run.nii');
            spm_file_merge(char(niiFiles), nii);
        end
    end
end

function ensure_dir(d)
    if exist(d, 'dir') ~= 7
        mkdir(d);
    end
end

function out = to_scan_refs(in)
    in = cellfun(@strtrim, in, 'UniformOutput', false);
    out = cellfun(@to_scan_ref, in, 'UniformOutput', false);
end

function out = to_scan_ref(in)
    if contains(in, ',')
        out = in;
    else
        out = [in ',1'];
    end
end
