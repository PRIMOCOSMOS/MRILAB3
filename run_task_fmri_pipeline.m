function run_task_fmri_pipeline()
% RUN_TASK_FMRI_PIPELINE
% -------------------------------------------------------------------------
% Standalone task-fMRI pipeline (MATLAB 2025a, SPM25-compatible environment)
% 说明：
% 1) 本脚本不调用 DPABI/SPM 预处理与统计函数；
% 2) 每一步显式实现核心计算逻辑，便于学习和复现；
% 3) 使用 MATLAB 内置能力（niftiread/niftiwrite/imreg*/imwarp/volshow 等）。
% -------------------------------------------------------------------------

    cfg = task_fmri_pipeline_config();
    ensure_dir(cfg.paths.derivativeDir);
    tplCtx = resolve_template_context(cfg);

    subjects = discover_subjects(cfg.paths);
    assert(~isempty(subjects), '未发现被试目录（已检查 FunRaw/T1Raw 与 DataRaw 回退路径）。');

    fprintf('\n=== [1/4] 预处理与配准 ===\n');
    subjData = struct([]);
    for i = 1:numel(subjects)
        subjData(i) = preprocess_subject(subjects{i}, cfg, tplCtx); %#ok<AGROW>
    end

    fprintf('\n=== [2/4] 模板构建与标准化 ===\n');
    tpl = build_group_template(subjData, cfg, tplCtx);
    for i = 1:numel(subjData)
        subjData(i) = normalize_subject_to_template(subjData(i), tpl, cfg);
    end

    fprintf('\n=== [3/4] 一级 GLM 分析 ===\n');
    for i = 1:numel(subjData)
        subjData(i).glm = first_level_glm(subjData(i), cfg);
    end

    fprintf('\n=== [4/4] 现代化可视化 ===\n');
    for i = 1:numel(subjData)
        visualize_subject_result(subjData(i), tpl, cfg);
    end

    fprintf('\nPipeline finished.\n');
end

% ================================ 被试发现 ================================
function subjects = discover_subjects(paths)
    if isfolder(paths.funcRawDir) && isfolder(paths.t1RawDir)
        funSubs = list_subdirs(paths.funcRawDir);
        t1Subs = list_subdirs(paths.t1RawDir);
        if ~isempty(funSubs) && ~isempty(t1Subs)
            subjects = sort(intersect(funSubs, t1Subs));
            return;
        end
    end
    subjects = sort(list_subdirs(paths.dataRawDir));
end

% ============================ 单被试预处理块 ==============================
function out = preprocess_subject(subjID, cfg, tplCtx)
    fprintf('[%s] preprocess...\n', subjID);

    if isfolder(cfg.paths.funcRawDir) && isfolder(cfg.paths.t1RawDir)
        anatDir = resolve_subject_input_dir(cfg.paths.t1RawDir, subjID);
        funcDir = resolve_subject_input_dir(cfg.paths.funcRawDir, subjID);
    else
        subjRaw = fullfile(cfg.paths.dataRawDir, subjID);
        anatDir = locate_first_existing(subjRaw, cfg.paths.anatDirCandidates);
        funcDir = locate_first_existing(subjRaw, cfg.paths.funcDirCandidates);
    end
    subjDer = fullfile(cfg.paths.derivativeDir, subjID);
    ensure_dir(subjDer);
    assert(~isempty(anatDir) && ~isempty(funcDir), '被试 %s 的原始结构像或功能像目录不完整。', subjID);

    % --- 读取并保存结构像 ---
    [anatVol, anatInfo] = read_single_volume(anatDir, fullfile(subjDer, 'anat'));

    % --- 读取功能像 run ---
    [funcRuns, funcInfo] = read_functional_runs(funcDir, fullfile(subjDer, 'func_raw'));
    assert(~isempty(funcRuns), '被试 %s 无功能像。', subjID);

    % --- run 级预处理 ---
    preprocRuns = cell(numel(funcRuns), 1);
    motionRuns = cell(numel(funcRuns), 1);
    for r = 1:numel(funcRuns)
        run4d = single(funcRuns{r});
        run4d = module_remove_first_timepoints(run4d, cfg.preproc);
        run4d = module_slice_timing(run4d, cfg.preproc);
        [run4d, motion] = module_realign(run4d, cfg.preproc);
        preprocRuns{r} = run4d;
        motionRuns{r} = motion;
    end

    % --- 合并 run ---
    all4d = cat(4, preprocRuns{:});
    motionAll = vertcat(motionRuns{:});

    % --- 结构像到功能均值配准 ---
    meanFunc = mean(all4d, 4, 'omitnan');
    [anatCoreg, tformAF] = module_coregister(anatVol, meanFunc, cfg.preproc);

    % --- 结构像分割（GM/WM/CSF） ---
    seg = module_segment(anatCoreg, cfg.preproc, tplCtx);

    % --- 写中间结果 ---
    write_nifti_like(anatCoreg, anatInfo, fullfile(subjDer, 'anat', 'anat_coreg.nii'));
    write_nifti_like(meanFunc, funcInfo{1}, fullfile(subjDer, 'func_preproc', 'mean_func.nii'));
    save(fullfile(subjDer, 'func_preproc', 'motion_6dof.mat'), 'motionAll');
    write_nifti_like(seg.gm, anatInfo, fullfile(subjDer, 'anat', 'gm_prob.nii'));
    write_nifti_like(seg.wm, anatInfo, fullfile(subjDer, 'anat', 'wm_prob.nii'));
    write_nifti_like(seg.csf, anatInfo, fullfile(subjDer, 'anat', 'csf_prob.nii'));

    out.subjID = subjID;
    out.subjDer = subjDer;
    out.anat = anatCoreg;
    out.anatInfo = anatInfo;
    out.func4d = all4d;
    out.funcInfo = funcInfo{1};
    out.motion = motionAll;
    out.seg = seg;
    out.tformAF = tformAF;
end

function V = module_remove_first_timepoints(V, P)
    if isfield(P, 'pipeline') && isfield(P.pipeline, 'removeFirstTimePoints') ...
            && isfield(P.pipeline.removeFirstTimePoints, 'enabled') ...
            && P.pipeline.removeFirstTimePoints.enabled
        V = drop_initial_volumes(V, P.pipeline.removeFirstTimePoints.nVolumes);
    end
end

function V = module_slice_timing(V, P)
    if isfield(P, 'pipeline') && isfield(P.pipeline, 'sliceTiming') ...
            && isfield(P.pipeline.sliceTiming, 'enabled') ...
            && P.pipeline.sliceTiming.enabled
        st = P.pipeline.sliceTiming;
        st.TR = P.TR;
        st.sliceTimingMode = st.mode;
        st.sliceTimingMs = st.timingMs;
        V = slice_timing_correction(V, st);
    end
end

function [V, motion] = module_realign(V, P)
    if isfield(P, 'pipeline') && isfield(P.pipeline, 'realign') ...
            && isfield(P.pipeline.realign, 'enabled') ...
            && P.pipeline.realign.enabled
        rp = P.pipeline.realign;
        rp.realignPyramidLevels = rp.pyramidLevels;
        [V, motion] = realign_4d(V, rp);
    else
        motion = zeros(size(V, 4), 6, 'single');
    end
end

function [anatReg, tform] = module_coregister(anat, meanFunc, P)
    if isfield(P, 'pipeline') && isfield(P.pipeline, 'coregister') ...
            && isfield(P.pipeline.coregister, 'enabled') ...
            && P.pipeline.coregister.enabled
        cp = P.pipeline.coregister;
        cp.coregPyramidLevels = cp.pyramidLevels;
        [anatReg, tform] = coregister_anat_to_func(anat, meanFunc, cp);
    else
        anatReg = anat;
        tform = affine3d(eye(4));
    end
end

function seg = module_segment(anat, P, tplCtx)
    if isfield(P, 'pipeline') && isfield(P.pipeline, 'segment') ...
            && isfield(P.pipeline.segment, 'enabled') ...
            && P.pipeline.segment.enabled
        sp = P.pipeline.segment;
        seg = segment_t1(anat, sp, tplCtx.segmentation);
    else
        seg.csf = zeros(size(anat), 'single');
        seg.gm = zeros(size(anat), 'single');
        seg.wm = zeros(size(anat), 'single');
    end
end

% ============================ 群体模板构建块 ==============================
function tpl = build_group_template(subjData, cfg, tplCtx)
    fprintf('[group] building template...\n');
    ensure_dir(cfg.paths.templateDir);

    % 初始模板：所有受试者结构像平均
    allAnat = cat(4, subjData.anat);
    template = mean(allAnat, 4, 'omitnan');

    if ~cfg.preproc.pipeline.normalize.enabled
        tpl.volume = template;
        tpl.info = subjData(1).anatInfo;
        tpl.groupTemplate = template;
        tpl.groupInfo = subjData(1).anatInfo;
        tpl.groupToTargetField = [];
        tpl.hasTargetTemplate = false;
        tpl.path = fullfile(cfg.paths.templateDir, 'group_template.nii');
        write_nifti_like(template, tpl.info, tpl.path);
        return;
    end

    % 迭代：每轮将各受试者 anat 非线性配准到当前模板，再更新模板
    NP = cfg.preproc.pipeline.normalize;
    for it = 1:NP.templateIters
        warpedAll = zeros([size(template), numel(subjData)], 'single');
        for i = 1:numel(subjData)
            moving = subjData(i).anat;
            [D, movingReg] = imregdemons(moving, template, ...
                NP.demonsIters, ...
                'AccumulatedFieldSmoothing', NP.demonsSmoothing);
            warpedAll(:, :, :, i) = movingReg;
            subjData(i).normField = D; %#ok<AGROW>
        end
        template = mean(warpedAll, 4, 'omitnan');
        fprintf('[group] template iteration %d done.\n', it);
    end

    tpl.volume = template;
    tpl.info = subjData(1).anatInfo;
    tpl.groupTemplate = template;
    tpl.groupInfo = subjData(1).anatInfo;
    tpl.groupToTargetField = [];
    tpl.hasTargetTemplate = false;
    tpl.path = fullfile(cfg.paths.templateDir, 'group_template.nii');
    write_nifti_like(template, tpl.info, tpl.path);

    if tplCtx.normalize.hasTargetTemplate
        fprintf('[group] mapping group template to target template: %s\n', tplCtx.normalize.targetPath);
        [Dgt, groupInTarget] = imregdemons(template, tplCtx.normalize.targetVolume, ...
            tplCtx.normalize.groupToTargetDemonsIters, ...
            'AccumulatedFieldSmoothing', tplCtx.normalize.groupToTargetDemonsSmoothing);
        tpl.groupToTargetField = Dgt;
        tpl.hasTargetTemplate = true;
        tpl.targetTemplate = tplCtx.normalize.targetVolume;
        tpl.targetInfo = tplCtx.normalize.targetInfo;
        tpl.targetPath = tplCtx.normalize.targetPath;
        tpl.volume = groupInTarget;
        tpl.info = tpl.targetInfo;
        tpl.path = fullfile(cfg.paths.templateDir, 'group_template_in_target_space.nii');
        write_nifti_like(groupInTarget, tpl.targetInfo, tpl.path);
        write_nifti_like(tplCtx.normalize.targetVolume, tpl.targetInfo, ...
            fullfile(cfg.paths.templateDir, 'target_template_used.nii'));
    end
end

% ============================== 标准化和平滑 ==============================
function out = normalize_subject_to_template(in, tpl, cfg)
    fprintf('[%s] normalize + smooth...\n', in.subjID);
    NP = cfg.preproc.pipeline.normalize;
    V = size(in.func4d, 4);
    D = [];
    if NP.enabled
        [D, ~] = imregdemons(in.anat, tpl.groupTemplate, ...
            NP.demonsIters, ...
            'AccumulatedFieldSmoothing', NP.demonsSmoothing);
        funcNorm = zeros(size(in.func4d), 'single');
        for t = 1:V
            funcNorm(:, :, :, t) = imwarp(in.func4d(:, :, :, t), D, 'cubic', ...
                'OutputView', imref3d(size(tpl.groupTemplate)));
        end
        if tpl.hasTargetTemplate
            funcNormTarget = zeros([size(tpl.targetTemplate), V], 'single');
            for t = 1:V
                funcNormTarget(:, :, :, t) = imwarp(funcNorm(:, :, :, t), tpl.groupToTargetField, ...
                    'cubic', 'OutputView', imref3d(size(tpl.targetTemplate)));
            end
            funcNorm = funcNormTarget;
        end
    else
        funcNorm = in.func4d;
    end

    SP = cfg.preproc.pipeline.smooth;
    if SP.enabled
        sigma = fwhm_to_sigma(SP.fwhm, SP.voxelSize);
        funcSmooth = zeros(size(funcNorm), 'single');
        for t = 1:V
            funcSmooth(:, :, :, t) = imgaussfilt3(funcNorm(:, :, :, t), sigma);
        end
    else
        funcSmooth = funcNorm;
    end

    out = in;
    out.normField = D;
    out.funcNorm = funcNorm;
    out.funcSmooth = funcSmooth;

    outDir = fullfile(in.subjDer, 'func_norm');
    ensure_dir(outDir);
    if tpl.hasTargetTemplate
        outRefInfo = tpl.targetInfo;
    else
        outRefInfo = tpl.info;
    end
    write_nifti_4d(funcSmooth, outRefInfo, fullfile(outDir, 'func_smooth_norm_4d.nii'));
end

% ============================== 一级统计分析 ==============================
function glm = first_level_glm(subj, cfg)
    fprintf('[%s] first-level GLM...\n', subj.subjID);

    C = resolve_first_level_design(subj.subjID, cfg.firstLevel);

    Y4d = subj.funcSmooth;
    T = size(Y4d, 4);
    Y = reshape(Y4d, [], T)';         % T x V

    % --- 构建设计矩阵 ---
    Xtask = build_task_regressors(C.names, C.onsetsSec, C.durationsSec, T, cfg.firstLevel);
    Xnuis = build_nuisance_regressors(subj.motion, T, cfg.firstLevel);
    X = [Xtask, Xnuis];
    X = zscore_cols(X);
    X = [X, ones(T, 1)];              % 截距

    % --- 高通滤波 ---
    [Yhp, Xhp] = highpass_dct(Y, X, cfg.firstLevel.hpf, cfg.preproc.TR);

    % --- AR(1) 预白化 ---
    rho = estimate_global_ar1(Yhp);
    [Yw, Xw] = ar1_prewhiten(Yhp, Xhp, rho);

    % --- GLS(whitened OLS) ---
    XtX = Xw' * Xw;
    XtY = Xw' * Yw;
    beta = XtX \ XtY;                  % P x V
    res = Yw - Xw * beta;              % T x V
    dof = size(Xw, 1) - rank(Xw);
    sigma2 = sum(res.^2, 1) / max(dof, 1);   % 1 x V

    % --- Contrast/T-map ---
    c = cfg.firstLevel.contrastWeights(:);
    assert(numel(c) == size(Xtask, 2), 'contrastWeights 与任务回归列数不一致。');
    cFull = [c; zeros(size(Xnuis, 2), 1); 0];
    cbeta = cFull' * beta;                             % 1 x V
    varC = (cFull' * (XtX \ cFull));                   % scalar
    tmap = cbeta ./ sqrt(max(eps, sigma2 * varC));     % 1 x V
    tmap3 = reshape(single(tmap), size(Y4d, 1), size(Y4d, 2), size(Y4d, 3));

    % --- 双阈值：voxel p + cluster extent ---
    mask = threshold_tmap(tmap3, dof, cfg.results.voxelPThreshold, cfg.results.clusterExtent);

    outDir = fullfile(subj.subjDer, 'first_level');
    ensure_dir(outDir);
    write_nifti_like(tmap3, subj.funcInfo, fullfile(outDir, 'tmap.nii'));
    write_nifti_like(single(mask), subj.funcInfo, fullfile(outDir, 'mask_thresholded.nii'));

    glm.X = X;
    glm.beta = beta;
    glm.dof = dof;
    glm.rho = rho;
    glm.tmap = tmap3;
    glm.mask = mask;
    glm.outDir = outDir;
end

function C = resolve_first_level_design(subjID, P)
    assert(isfield(P, 'design') && ~isempty(P.design), 'firstLevel.design 未配置。');
    D = P.design;
    assert(isfield(D, 'names') && isfield(D, 'onsets') && isfield(D, 'durations'), ...
        'firstLevel.design 需要包含 names/onsets/durations。');

    names = D.names;
    onsets = D.onsets;
    durations = D.durations;
    assert(iscell(names) && iscell(onsets) && iscell(durations), ...
        'firstLevel.design 中 names/onsets/durations 必须为 cell。');
    assert(numel(names) == numel(onsets) && numel(names) == numel(durations), ...
        'firstLevel.design 条件数不一致。');

    units = 'secs';
    if isfield(P, 'timingUnits') && ~isempty(P.timingUnits)
        units = lower(char(P.timingUnits));
    end
    assert(any(strcmp(units, {'secs','seconds','scans'})), ...
        'timingUnits 必须为 ''secs''、''seconds'' 或 ''scans''。');

    onsetsSec = cell(size(onsets));
    durationsSec = cell(size(durations));
    for i = 1:numel(names)
        condName = char(names{i});
        on = double(onsets{i}(:)');
        du = double(durations{i}(:)');
        assert(~isempty(on), '条件 %s 的 onsets 不能为空。', condName);
        assert(~isempty(du), '条件 %s 的 durations 不能为空。', condName);
        if isscalar(du)
            du = repmat(du, size(on));
        end
        assert(numel(du) == numel(on), ...
            '条件 %s 的 durations 若非标量，长度需与 onsets 一致。', condName);

        if strcmp(units, 'scans')
            onsetBase = 0;
            if isfield(P, 'scanOnsetIndexBase') && ~isempty(P.scanOnsetIndexBase)
                onsetBase = double(P.scanOnsetIndexBase);
            end
            assert(any(onsetBase == [0 1]), ...
                'scanOnsetIndexBase 必须为 0 或 1。');
            assert(all(on >= onsetBase), ...
                '存在小于 scanOnsetIndexBase 的 onset。');
            on = (on - onsetBase) * P.TR;
            du = du * P.TR;
        end

        onsetsSec{i} = on;
        durationsSec{i} = du;
    end

    C.subjID = subjID;
    C.names = names;
    C.onsetsSec = onsetsSec;
    C.durationsSec = durationsSec;
end

% ============================== 现代可视化块 ==============================
function visualize_subject_result(subj, tpl, cfg)
    fprintf('[%s] visualization...\n', subj.subjID);
    tmap = subj.glm.tmap;
    mask = subj.glm.mask;
    anat = tpl.volume;

    outDir = fullfile(subj.subjDer, 'visualization');
    ensure_dir(outDir);

    % 峰值坐标表
    peaks = extract_peaks(tmap, mask, cfg.results.numPeaks);
    writetable(peaks, fullfile(outDir, 'activation_peaks.csv'));

    % 图1：多切面叠加热图
    f1 = figure('Visible', 'off', 'Color', 'w', 'Position', [80 80 1400 900]);
    tiledlayout(2, 2, 'Padding', 'compact', 'TileSpacing', 'compact');
    show_overlay_slices(anat, tmap, mask, cfg.results);
    exportgraphics(f1, fullfile(outDir, 'modern_slice_overlay.png'), 'Resolution', 300);
    close(f1);

    % 图2：3D 表面 + 激活团块
    f2 = figure('Visible', 'off', 'Color', 'w', 'Position', [80 80 1400 900]);
    show_3d_surfaces(anat, tmap, mask);
    exportgraphics(f2, fullfile(outDir, 'modern_3d_surface.png'), 'Resolution', 300);
    close(f2);

    % 图3：体绘制（若可用）
    if cfg.results.enableVolshow && exist('volshow', 'file') == 2
        try
            f3 = figure('Visible', 'off', 'Color', 'w', 'Position', [80 80 1400 900]);
            volshow(abs(tmap) .* single(mask), 'BackgroundColor', [0 0 0], ...
                'Colormap', turbo(256), 'Renderer', 'VolumeRendering');
            exportgraphics(f3, fullfile(outDir, 'modern_volshow.png'), 'Resolution', 300);
            close(f3);
        catch
            warning('volshow 导出失败，已跳过。');
        end
    end
end

% ============================== 模板解析与管理 ==============================
function tplCtx = resolve_template_context(cfg)
    tplCtx = struct();
    tplCtx.segmentation = struct('enabled', false, 'hasTpm', false, ...
        'tpmPath', '', 'tpmVolumeIndices', [1 2 3], 'priorWeight', 0);
    tplCtx.normalize = struct('hasTargetTemplate', false, 'targetPath', '', ...
        'targetInfo', struct(), 'targetVolume', [], ...
        'groupToTargetDemonsIters', [80 40 20], ...
        'groupToTargetDemonsSmoothing', 1.0);

    if ~isfield(cfg, 'templates') || ~cfg.templates.enabled
        fprintf('[template] disabled.\n');
        return;
    end

    Tc = cfg.templates;
    candidateFiles = discover_nifti_templates(Tc.searchDirs);

    % ---- segmentation priors (TPM) ----
    if isfield(Tc, 'segmentation') && Tc.segmentation.enabled
        tplCtx.segmentation.enabled = true;
        tplCtx.segmentation.tpmVolumeIndices = Tc.segmentation.tpmVolumeIndices;
        tplCtx.segmentation.priorWeight = Tc.segmentation.priorWeight;

        tpmPath = '';
        if isfield(Tc.segmentation, 'tpmPath') && ~isempty(Tc.segmentation.tpmPath) ...
                && isfile(Tc.segmentation.tpmPath)
            tpmPath = Tc.segmentation.tpmPath;
        else
            tpmPath = select_best_template(candidateFiles, {'tpm'}, ...
                logical_if_exists(Tc.segmentation, 'preferEastAsian'), ...
                {'east', 'asian', 'china', 'chinese', 'eastern'});
        end

        if ~isempty(tpmPath)
            tplCtx.segmentation.hasTpm = true;
            tplCtx.segmentation.tpmPath = tpmPath;
            fprintf('[template] segmentation TPM: %s\n', tpmPath);
        else
            fprintf('[template] segmentation TPM not found, fallback to intensity-only segmentation.\n');
        end
    end

    % ---- normalization target template (e.g., MNI152) ----
    if isfield(Tc, 'normalize')
        tplCtx.normalize.groupToTargetDemonsIters = Tc.normalize.groupToTargetDemonsIters;
        tplCtx.normalize.groupToTargetDemonsSmoothing = Tc.normalize.groupToTargetDemonsSmoothing;
        tgtPath = '';
        if isfield(Tc.normalize, 'targetTemplatePath') && ~isempty(Tc.normalize.targetTemplatePath) ...
                && isfile(Tc.normalize.targetTemplatePath)
            tgtPath = Tc.normalize.targetTemplatePath;
        else
            tgtPath = select_best_template(candidateFiles, {'mni', 'icbm'}, ...
                logical_if_exists(Tc.normalize, 'preferMNI'), ...
                {'mni152', 'mni', 'icbm', 'template'});
            if isempty(tgtPath)
                tgtPath = select_best_template(candidateFiles, {'template'}, ...
                    logical_if_exists(Tc.normalize, 'preferMNI'), ...
                    {'mni152', 'mni', 'icbm', 'template'});
            end
        end

        if ~isempty(tgtPath)
            [tgtVol, tgtInfo] = read_template_volume(tgtPath);
            tplCtx.normalize.hasTargetTemplate = true;
            tplCtx.normalize.targetPath = tgtPath;
            tplCtx.normalize.targetVolume = tgtVol;
            tplCtx.normalize.targetInfo = tgtInfo;
            fprintf('[template] normalization target: %s\n', tgtPath);
        else
            fprintf('[template] normalization target not found, keep group-template space.\n');
        end
    end
end

function tf = logical_if_exists(S, fieldName)
    tf = false;
    if isfield(S, fieldName) && ~isempty(S.(fieldName))
        tf = logical(S.(fieldName));
    end
end

function files = discover_nifti_templates(searchDirs)
    files = {};
    if isempty(searchDirs)
        return;
    end
    for i = 1:numel(searchDirs)
        d = searchDirs{i};
        if ~ischar(d) && ~isstring(d)
            continue;
        end
        d = char(d);
        if ~isfolder(d)
            continue;
        end
        nii = dir(fullfile(d, '**', '*.nii'));
        niigz = dir(fullfile(d, '**', '*.nii.gz'));
        hit = [nii; niigz];
        for k = 1:numel(hit)
            files{end + 1} = fullfile(hit(k).folder, hit(k).name); %#ok<AGROW>
        end
    end
    files = unique(files);
end

function bestPath = select_best_template(files, mustContain, preferSpecial, specialKeywords)
    specialBoostWhenPreferred = 5;     % 开启偏好时，显著提升特定关键词（如 East/MNI）候选得分
    specialBoostWhenNotPreferred = 1;  % 未开启偏好时，仍保留轻微加分，避免同分情况下完全忽略
    bestPath = '';
    bestScore = -inf;
    for i = 1:numel(files)
        p = lower(files{i});
        hitAny = false;
        for j = 1:numel(mustContain)
            if contains(p, lower(mustContain{j}))
                hitAny = true;
            end
        end
        if ~hitAny
            continue;
        end

        score = 0;
        for j = 1:numel(mustContain)
            if contains(p, lower(mustContain{j})), score = score + 3; end
        end
        if contains(p, 'template'), score = score + 2; end
        if contains(p, '152'), score = score + 1; end

        specialHit = false;
        for j = 1:numel(specialKeywords)
            if contains(p, lower(specialKeywords{j}))
                specialHit = true;
                break;
            end
        end
        if specialHit
            if preferSpecial
                score = score + specialBoostWhenPreferred;
            else
                score = score + specialBoostWhenNotPreferred;
            end
        end
        if score > bestScore
            bestScore = score;
            bestPath = files{i};
        end
    end
end

function [V, info] = read_template_volume(pathIn)
    info = niftiinfo(pathIn);
    V = single(niftiread(info));
    if ndims(V) == 4
        % 目标标准空间模板通常为 3D 参考体；若是 4D 文件，默认取首个体积作为空间参考。
        V = V(:, :, :, 1);
    end
end

function pri = load_segmentation_priors(anatSize, segTpl)
    pri = struct('available', false, 'csf', [], 'gm', [], 'wm', []);
    if ~isfield(segTpl, 'enabled') || ~segTpl.enabled || ~segTpl.hasTpm
        return;
    end
    if ~isfile(segTpl.tpmPath)
        return;
    end
    try
        info = niftiinfo(segTpl.tpmPath);
        tpm = single(niftiread(info));
        if ndims(tpm) < 4 || size(tpm, 4) < 3
            warning('TPM 文件不是有效4D组织概率模板：%s', segTpl.tpmPath);
            return;
        end
        % TPM 的组织顺序由配置显式定义（默认: 1->GM, 2->WM, 3->CSF）；
        % 若所用 TPM 顺序不同，请在 cfg.templates.segmentation.tpmVolumeIndices 中调整。
        idx = double(segTpl.tpmVolumeIndices(:)');
        assert(numel(idx) == 3 && all(idx >= 1) && all(mod(idx, 1) == 0) && all(idx <= size(tpm, 4)), ...
            'cfg.templates.segmentation.tpmVolumeIndices 非法：当前值=%s；要求=长度为3的正整数，且每个索引在 [1, %d] 范围内。', ...
            mat2str(idx), size(tpm, 4));
        gm = imresize3(tpm(:, :, :, idx(1)), anatSize, 'linear');
        wm = imresize3(tpm(:, :, :, idx(2)), anatSize, 'linear');
        csf = imresize3(tpm(:, :, :, idx(3)), anatSize, 'linear');
        gm = max(gm, 0); wm = max(wm, 0); csf = max(csf, 0);
        ps = gm + wm + csf + eps('single');
        pri.gm = gm ./ ps;
        pri.wm = wm ./ ps;
        pri.csf = csf ./ ps;
        pri.available = true;
    catch ME
        warning('加载 TPM 先验失败（文件: %s，已回退强度分割）：%s', segTpl.tpmPath, ME.message);
    end
end

% ============================== I/O 与工具函数 ============================
function [V, info] = read_single_volume(inputDir, outDir)
    ensure_dir(outDir);
    outPath = fullfile(outDir, 'anat_raw.nii');
    wroteFromDicom = false;
    nii = dir(fullfile(inputDir, '*.nii*'));
    if ~isempty(nii)
        p = fullfile(nii(1).folder, nii(1).name);
        info = niftiinfo(p);
        V = single(niftiread(info));
        if ndims(V) == 4, V = V(:, :, :, 1); end
    else
        dcm = dir(fullfile(inputDir, '**', '*.dcm'));
        assert(~isempty(dcm), '结构像目录中无 NIfTI 或 DICOM。');
        [V, info] = convert_dicom_to_nifti_3d(dcm, outPath);
        wroteFromDicom = true;
    end
    if ~wroteFromDicom
        write_nifti_like(V, info, outPath);
    end
end

function V = drop_initial_volumes(V, n)
    assert(size(V, 4) > n, '时间点不足以去除前 %d 个。', n);
    V = V(:, :, :, n+1:end);
end

function Vout = slice_timing_correction(Vin, P)
    nx = size(Vin,1); ny = size(Vin,2); nz = size(Vin,3); nt = size(Vin,4);
    assert(nz == P.nslices, '功能像切片数(%d)与配置 nslices(%d)不一致。', nz, P.nslices);
    t = (0:nt-1) * P.TR;
    shifts = resolve_slice_timing_shifts(P, nz);

    Vout = zeros(size(Vin), 'single');
    for z = 1:nz
        tq = t + shifts(z);
        tq = min(max(tq, t(1)), t(end));
        for x = 1:nx
            for y = 1:ny
                s = squeeze(Vin(x, y, z, :));
                Vout(x, y, z, :) = interp1(t, s, tq, 'pchip', 'extrap');
            end
        end
    end
end

function shiftsSec = resolve_slice_timing_shifts(P, nz)
    mode = 'order'; % 回退默认：当配置未显式提供 sliceTimingMode 时使用传统顺序模式
    if isfield(P, 'sliceTimingMode') && ~isempty(P.sliceTimingMode)
        mode = lower(string(P.sliceTimingMode));
    elseif isfield(P, 'mode') && ~isempty(P.mode)
        mode = lower(string(P.mode));
    end

    if mode == "timing_ms"
        if isfield(P, 'sliceTimingMs') && ~isempty(P.sliceTimingMs)
            stMs = double(P.sliceTimingMs(:)');
        elseif isfield(P, 'timingMs') && ~isempty(P.timingMs)
            stMs = double(P.timingMs(:)');
        else
            error('sliceTimingMode=timing_ms 时必须提供 sliceTimingMs/timingMs。');
        end
        assert(numel(stMs) == nz, 'sliceTimingMs 长度(%d)必须等于切片数(%d)。', numel(stMs), nz);

        if isfield(P, 'refTimingMs') && ~isempty(P.refTimingMs)
            refMs = double(P.refTimingMs);
        else
            refIdx = P.refSlice;
            assert(refIdx >= 1 && refIdx <= nz, ...
                'refSlice 越界：%d (有效范围 1-%d)。', refIdx, nz);
            refMs = stMs(refIdx);
        end
        shiftsSec = (stMs - refMs) / 1000;
        return;
    end

    assert(isfield(P, 'sliceOrder') && ~isempty(P.sliceOrder), ...
        'sliceTimingMode=order 时必须提供 sliceOrder。');
    so = double(P.sliceOrder(:)');
    assert(numel(so) == nz, 'sliceOrder 长度(%d)必须等于切片数(%d)。', numel(so), nz);
    counts = histcounts(so, 0.5:1:(nz+0.5));
    missingSlices = find(counts == 0);
    duplicateSlices = find(counts > 1);
    assert(isempty(missingSlices) && isempty(duplicateSlices), ...
        'sliceOrder 非法。缺失切片: [%s]；重复切片: [%s]。', ...
        num2str(missingSlices), num2str(duplicateSlices));
    assert(isfield(P, 'refSlice') && ~isempty(P.refSlice), ...
        'sliceTimingMode=order 时必须提供 refSlice。');
    refOrd = find(so == P.refSlice, 1, 'first');
    dt = P.TR / nz;
    shiftsSec = zeros(1, nz);
    for z = 1:nz
        ord = find(so == z, 1, 'first');
        shiftsSec(z) = (ord - refOrd) * dt;
    end
end

function [Vout, motion] = realign_4d(Vin, P)
    nt = size(Vin, 4);
    ref = Vin(:, :, :, 1);
    Rfixed = imref3d(size(ref));
    optimizer = registration.optimizer.OnePlusOneEvolutionary();
    optimizer.InitialRadius = 6.25e-3;
    metric = registration.metric.MeanSquares();

    Vout = zeros(size(Vin), 'single');
    motion = zeros(nt, 6, 'single'); % tx ty tz rx ry rz (approx)
    Vout(:, :, :, 1) = ref;

    for t = 2:nt
        moving = Vin(:, :, :, t);
        tform = imregtform(moving, ref, 'rigid', optimizer, metric, ...
            'PyramidLevels', P.realignPyramidLevels);
        Vout(:, :, :, t) = imwarp(moving, tform, 'OutputView', Rfixed, ...
            'InterpolationMethod', 'cubic');
        motion(t, :) = affine_to_6dof(tform.A);
    end
end

function [anatReg, tform] = coregister_anat_to_func(anat, meanFunc, P)
    Rfixed = imref3d(size(meanFunc));
    optimizer = registration.optimizer.OnePlusOneEvolutionary();
    optimizer.InitialRadius = 6.25e-3;
    metric = registration.metric.MattesMutualInformation();
    tform = imregtform(anat, meanFunc, 'rigid', optimizer, metric, ...
        'PyramidLevels', P.coregPyramidLevels);
    anatReg = imwarp(anat, tform, 'OutputView', Rfixed, 'InterpolationMethod', 'cubic');
end

function seg = segment_t1(anat, P, segTpl)
    anat = single(anat);
    anat = anat - min(anat(:));
    anat = anat / max(eps, max(anat(:)));

    % 低频 bias 场校正（近似）
    bias = imgaussfilt3(anat, P.biasFieldSigma);
    anatCorr = anat ./ max(bias, eps('single'));
    anatCorr = anatCorr / max(eps, max(anatCorr(:)));

    % 3类 GMM 近似（用 kmeans 初始化）
    x = anatCorr(:);
    x = x(isfinite(x));
    [idx, c] = kmeans(x, 3, 'Replicates', 5, 'MaxIter', 200);
    [~, clustersByIntensity] = sort(c, 'ascend');
    mu = zeros(3,1); sd = zeros(3,1);
    for k = 1:3
        v = x(idx==clustersByIntensity(k));
        mu(k) = mean(v); sd(k) = std(v) + 1e-4;
    end

    p1 = normpdf(anatCorr, mu(1), sd(1));
    p2 = normpdf(anatCorr, mu(2), sd(2));
    p3 = normpdf(anatCorr, mu(3), sd(3));
    ps = p1 + p2 + p3 + eps('single');
    p1 = p1 ./ ps; p2 = p2 ./ ps; p3 = p3 ./ ps;

    % 经验映射：低灰度=CSF，中灰度=GM，高灰度=WM
    seg.csf = p1;
    seg.gm = p2;
    seg.wm = p3;

    % 可选：融合模板先验（TPM），提升组织概率图稳定性
    pri = load_segmentation_priors(size(anatCorr), segTpl);
    if pri.available
        % 将 priorWeight 限制在 [0,1]，确保“强度似然/模板先验”线性融合始终稳定。
        w = min(max(double(segTpl.priorWeight), 0), 1);
        seg.csf = (1 - w) * seg.csf + w * pri.csf;
        seg.gm = (1 - w) * seg.gm + w * pri.gm;
        seg.wm = (1 - w) * seg.wm + w * pri.wm;
        ps2 = seg.csf + seg.gm + seg.wm + eps('single');
        seg.csf = seg.csf ./ ps2;
        seg.gm = seg.gm ./ ps2;
        seg.wm = seg.wm ./ ps2;
    end
end

function Xtask = build_task_regressors(names, onsetsSec, durationsSec, T, P)
    assert(P.microtimeResolution >= 1 && mod(P.microtimeResolution, 1) == 0, ...
        'microtimeResolution 必须为正整数。');
    assert(P.microtimeOnset >= 1 && P.microtimeOnset <= P.microtimeResolution ...
        && mod(P.microtimeOnset, 1) == 0, ...
        'microtimeOnset 必须是 1..microtimeResolution 的整数。');
    dt = P.TR / P.microtimeResolution;
    nMicro = T * P.microtimeResolution;
    tMicro = (0:nMicro-1)' * dt;
    h = canonical_hrf(dt, P.hrf);
    % tMicro 的索引1对应时间0, 因此 fmri_t0 直接映射为每TR内的采样索引（1-based）
    sampleIdx = (0:T-1) * P.microtimeResolution + P.microtimeOnset;
    badIdx = sampleIdx(sampleIdx < 1 | sampleIdx > nMicro);
    assert(isempty(badIdx), ...
        'microtimeOnset 越界: fmri_t=%d, fmri_t0=%d, nMicro=%d, badIdx=[%s]', ...
        P.microtimeResolution, P.microtimeOnset, nMicro, num2str(badIdx));

    Xtask = zeros(T, numel(names));
    for i = 1:numel(names)
        uMicro = zeros(nMicro, 1);
        on = onsetsSec{i};
        du = durationsSec{i};
        for k = 1:numel(on)
            uMicro = uMicro + double(tMicro >= on(k) & tMicro < (on(k) + du(k)));
        end
        xc = conv(uMicro, h);
        xc = xc(1:nMicro);
        Xtask(:, i) = xc(sampleIdx(:));
    end
end

function Xn = build_nuisance_regressors(motion, T, P)
    m = motion;
    if size(m,1) ~= T
        if size(m,1) > T
            m = m(1:T, :);
        else
            m = [m; zeros(T-size(m,1), size(m,2), 'like', m)];
        end
    end
    d = [zeros(1, size(m,2)); diff(m,1,1)];
    Xn = [m, d];

    if P.addQuadraticMotion
        Xn = [Xn, m.^2, d.^2];
    end

    if P.addLinearTrend
        Xn = [Xn, linspace(-1,1,T)'];
    end
end

function h = canonical_hrf(dt, H)
    % 双Gamma HRF:
    % p1/d1: 主峰 shape/scale; p2/d2: undershoot shape/scale; ratio: 相对幅值比
    t = (0:dt:H.length)';
    g1 = gampdf(t, H.p1, H.d1);
    g2 = gampdf(t, H.p2, H.d2);
    h = g1 - g2 / H.ratio;
    h = h / max(eps, sum(h));
end

function [Yf, Xf] = highpass_dct(Y, X, hpfSec, TR)
    T = size(X,1);
    n = floor(2 * (T * TR) / hpfSec + 1);
    C = zeros(T, n);
    for k = 1:n
        C(:, k) = cos((0:T-1)' * k * pi / T);
    end
    R = eye(T) - C * ((C' * C) \ C');
    Yf = R * Y;
    Xf = R * X;
end

function rho = estimate_global_ar1(Y)
    y = mean(Y, 2, 'omitnan');
    y1 = y(2:end); y0 = y(1:end-1);
    rho = (y0' * y1) / max(eps, (y0' * y0));
    rho = max(min(rho, 0.9), -0.9);
end

function [Yw, Xw] = ar1_prewhiten(Y, X, rho)
    Yw = Y(2:end, :) - rho * Y(1:end-1, :);
    Xw = X(2:end, :) - rho * X(1:end-1, :);
end

function mask = threshold_tmap(tmap, dof, pVoxel, kExtent)
    tThr = tinv(1 - pVoxel / 2, dof);
    mask = abs(tmap) > tThr;
    CC = bwconncomp(mask, 26);
    keep = false(size(mask));
    for i = 1:CC.NumObjects
        idx = CC.PixelIdxList{i};
        if numel(idx) >= kExtent
            keep(idx) = true;
        end
    end
    mask = keep;
end

function tbl = extract_peaks(tmap, mask, nPeaks)
    z = tmap;
    z(~mask) = -inf;
    [vals, idx] = maxk(z(:), nPeaks);
    [x, y, zc] = ind2sub(size(tmap), idx);
    tbl = table((1:numel(vals))', vals, x, y, zc, ...
        'VariableNames', {'Rank','TValue','I','J','K'});
end

function show_overlay_slices(anat, tmap, mask, R)
    idxX = round(size(anat,1) * [0.35 0.50]);
    idxY = round(size(anat,2) * [0.45 0.60]);
    idxZ = round(size(anat,3) * [0.35 0.50]);
    cm = turbo(256);

    nexttile(1); overlay_slice(squeeze(anat(idxX(1),:,:))', squeeze(tmap(idxX(1),:,:))', squeeze(mask(idxX(1),:,:))', cm, R);
    title('Sagittal');
    nexttile(2); overlay_slice(squeeze(anat(:,idxY(1),:))', squeeze(tmap(:,idxY(1),:))', squeeze(mask(:,idxY(1),:))', cm, R);
    title('Coronal');
    nexttile(3); overlay_slice(squeeze(anat(:,:,idxZ(1)))', squeeze(tmap(:,:,idxZ(1)))', squeeze(mask(:,:,idxZ(1)))', cm, R);
    title('Axial');
    nexttile(4); overlay_slice(squeeze(anat(:,:,idxZ(2)))', squeeze(tmap(:,:,idxZ(2)))', squeeze(mask(:,:,idxZ(2)))', cm, R);
    title('Axial+');
end

function overlay_slice(bg, ov, m, cm, R)
    imagesc(bg); axis image off; colormap(gca, gray); hold on;
    alphaMap = zeros(size(m), 'single');
    alphaMap(m>0) = R.overlayAlpha;
    ovn = rescale(ov);
    h = imagesc(ovn);
    set(h, 'AlphaData', alphaMap);
    hold off;
    c = colorbar; c.Label.String = 'Activation intensity';
    colormap(gca, cm);
end

function show_3d_surfaces(anat, tmap, mask)
    p1 = patch(isosurface(smooth3(anat), 0.35));
    p1.FaceColor = [0.75 0.75 0.8];
    p1.EdgeColor = 'none';
    p1.FaceAlpha = 0.18;
    hold on;

    act = abs(tmap) .* single(mask);
    if any(act(:) > 0)
        lv = prctile(act(act>0), 70);
        p2 = patch(isosurface(smooth3(act), lv));
        p2.FaceColor = [1 0.2 0.1];
        p2.EdgeColor = 'none';
        p2.FaceAlpha = 0.85;
    end
    camlight; camlight headlight; lighting gouraud;
    axis image off; view(3);
    title('3D cortical-like envelope + activation clusters');
end

function A6 = affine_to_6dof(A)
    tx = A(1,4); ty = A(2,4); tz = A(3,4);
    R = A(1:3,1:3);
    ry = asin(-R(3,1));
    rx = atan2(R(3,2), R(3,3));
    rz = atan2(R(2,1), R(1,1));
    A6 = single([tx ty tz rx ry rz]);
end

function X = zscore_cols(X)
    mu = mean(X, 1, 'omitnan');
    sd = std(X, 0, 1, 'omitnan');
    sd(sd < eps) = 1;
    X = (X - mu) ./ sd;
end

function sigma = fwhm_to_sigma(fwhmMM, voxMM)
    % FWHM = 2*sqrt(2*ln(2))*sigma
    fwhmToSigmaFactor = 2 * sqrt(2 * log(2));
    sigma = (fwhmMM ./ max(voxMM, eps)) / fwhmToSigmaFactor;
end

function p = locate_first_existing(baseDir, names)
    p = '';
    for i = 1:numel(names)
        c = fullfile(baseDir, names{i});
        if isfolder(c)
            p = c;
            return;
        end
    end
end

% 列出 baseDir 下所有非隐藏子目录名称。
% 输入：
%   baseDir - 父目录路径。
% 输出：
%   names   - 非隐藏子目录名称 cell 数组。
function names = list_subdirs(baseDir)
    names = {};
    if ~isfolder(baseDir)
        return;
    end
    d = dir(baseDir);
    d = d([d.isdir]);
    names = {d.name};
    names = names(~startsWith(names, '.'));
end

% 在指定模态根目录下解析被试输入目录。
% 输入：
%   rootDir - 模态根目录（例如 FunRaw 或 T1Raw）。
%   subjID  - 被试目录名。
% 输出：
%   p       - 解析后的目录路径；优先 rootDir/subjID；
%             若不存在且 rootDir 下直接有数据文件，则回退为 rootDir。
function p = resolve_subject_input_dir(rootDir, subjID)
    p = '';
    candidateSubj = fullfile(rootDir, subjID);
    if isfolder(candidateSubj)
        p = candidateSubj;
        return;
    end
    filesNii = dir(fullfile(rootDir, '*.nii*'));
    if ~isempty(filesNii)
        p = rootDir;
        return;
    end
    filesDcm = dir(fullfile(rootDir, '**', '*.dcm'));
    if ~isempty(filesDcm)
        p = rootDir;
    end
end

function ensure_dir(p)
    if ~isfolder(p)
        mkdir(p);
    end
end

function write_nifti_like(vol, refInfo, pathOut)
    ensure_dir(fileparts(pathOut));
    info = sanitize_nifti_info_for_write(refInfo);
    info.ImageSize = size(vol);
    if numel(info.ImageSize)==2, info.ImageSize(3)=1; end
    info.Datatype = 'single';
    safe_niftiwrite(single(vol), pathOut, info);
end

function write_nifti_4d(vol4d, refInfo, pathOut)
    ensure_dir(fileparts(pathOut));
    info = sanitize_nifti_info_for_write(refInfo);
    info.ImageSize = size(vol4d);
    info.Datatype = 'single';
    safe_niftiwrite(single(vol4d), pathOut, info);
end

function infoOut = sanitize_nifti_info_for_write(infoIn)
    infoOut = struct();
    if ~isstruct(infoIn)
        return;
    end
    % 仅保留 niftiwrite 可识别的通用字段，避免版本差异导致的字段报错（如 Description）。
    keep = {'ImageSize','PixelDimensions','Datatype','SpaceUnits','TimeUnits', ...
            'AdditiveOffset','MultiplicativeScaling','Transform','Qfactor'};
    for i = 1:numel(keep)
        k = keep{i};
        if isfield(infoIn, k)
            infoOut.(k) = infoIn.(k);
        end
    end
end

function safe_niftiwrite(vol, pathOut, info)
    infoTry = info;
    % Keep at least a few retries even for tiny headers, and allow a small buffer
    % beyond current field count for version-specific parser behavior.
    % 4: enough attempts for very small headers; +2: tolerate parser-side extra checks.
    minStripAttempts = 4;
    extraStripAttempts = 2;
    maxStripAttempts = max(minStripAttempts, numel(fieldnames(infoTry)) + extraStripAttempts);
    for k = 1:maxStripAttempts
        try
            niftiwrite(vol, pathOut, infoTry, 'Compressed', false);
            return;
        catch ME
            badField = extract_unknown_info_field(ME.message);
            if ~isempty(badField) && isfield(infoTry, badField)
                infoTry = rmfield(infoTry, badField);
                continue;
            end
            break;
        end
    end
    warning('niftiwrite:InfoFallback', ...
        'Failed to keep compatible NIfTI header for %s; fallback to write without Info.', pathOut);
    niftiwrite(vol, pathOut, 'Compressed', false);
end

function fieldName = extract_unknown_info_field(msg)
    fieldName = '';
    % Match localized messages (CN/EN); quoteClass includes: ", “, ”, and '.
    quoteClass = '"“”''';
    cnUnknownFieldPattern = ['无法识别的字段名称[:：\s]*[', quoteClass, ']?([^', quoteClass, '\s]+)'];
    t = regexp(msg, cnUnknownFieldPattern, 'tokens', 'once');
    if ~isempty(t)
        fieldName = t{1};
        return;
    end
    enUnknownFieldPattern = 'Unrecognized field name\s*''([^'']+)''';
    t = regexp(msg, enUnknownFieldPattern, 'tokens', 'once');
    if ~isempty(t)
        fieldName = t{1};
    end
end

function V = dicom_series_to_volume(dcmList)
    files = fullfile({dcmList.folder}, {dcmList.name});
    z = zeros(numel(files), 1);
    for i = 1:numel(files)
        h = dicominfo(files{i});
        if isfield(h, 'InstanceNumber'), z(i) = h.InstanceNumber; else, z(i)=i; end
    end
    [~, ord] = sort(z);
    files = files(ord);
    s = dicomread(files{1});
    V = zeros([size(s), numel(files)], 'single');
    for i = 1:numel(files), V(:,:,i) = single(dicomread(files{i})); end
end

function [V, info] = convert_dicom_to_nifti_3d(dcmList, outPath)
    ensure_dir(fileparts(outPath));
    V = dicom_series_to_volume(dcmList);
    niftiwrite(single(V), outPath, 'Compressed', false);
    info = niftiinfo(outPath);
    V = single(niftiread(info));
    assert(ndims(V) == 3, '结构像 DICOM 转 NIfTI 后应为 3D。');
end
