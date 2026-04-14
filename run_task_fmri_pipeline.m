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

    subjects = discover_subjects(cfg.paths.dataRawDir);
    assert(~isempty(subjects), '未发现被试目录: %s', cfg.paths.dataRawDir);

    fprintf('\n=== [1/4] 预处理与配准 ===\n');
    subjData = struct([]);
    for i = 1:numel(subjects)
        subjData(i) = preprocess_subject(subjects{i}, cfg); %#ok<AGROW>
    end

    fprintf('\n=== [2/4] 模板构建与标准化 ===\n');
    tpl = build_group_template(subjData, cfg);
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
function subjects = discover_subjects(dataRawDir)
    d = dir(dataRawDir);
    d = d([d.isdir]);
    names = {d.name};
    names = names(~startsWith(names, '.'));
    subjects = sort(names);
end

% ============================ 单被试预处理块 ==============================
function out = preprocess_subject(subjID, cfg)
    fprintf('[%s] preprocess...\n', subjID);

    subjRaw = fullfile(cfg.paths.dataRawDir, subjID);
    subjDer = fullfile(cfg.paths.derivativeDir, subjID);
    ensure_dir(subjDer);

    anatDir = locate_first_existing(subjRaw, cfg.paths.anatDirCandidates);
    funcDir = locate_first_existing(subjRaw, cfg.paths.funcDirCandidates);
    assert(~isempty(anatDir) && ~isempty(funcDir), '被试 %s 的 anat/func 目录不完整。', subjID);

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
        run4d = drop_initial_volumes(run4d, cfg.preproc.removeFirstN);
        run4d = slice_timing_correction(run4d, cfg.preproc);
        [run4d, motion] = realign_4d(run4d, cfg.preproc);
        preprocRuns{r} = run4d;
        motionRuns{r} = motion;
    end

    % --- 合并 run ---
    all4d = cat(4, preprocRuns{:});
    motionAll = vertcat(motionRuns{:});

    % --- 结构像到功能均值配准 ---
    meanFunc = mean(all4d, 4, 'omitnan');
    [anatCoreg, tformAF] = coregister_anat_to_func(anatVol, meanFunc, cfg.preproc);

    % --- 结构像分割（GM/WM/CSF） ---
    seg = segment_t1(anatCoreg, cfg.preproc);

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

% ============================ 群体模板构建块 ==============================
function tpl = build_group_template(subjData, cfg)
    fprintf('[group] building template...\n');
    ensure_dir(cfg.paths.templateDir);

    % 初始模板：所有受试者结构像平均
    allAnat = cat(4, subjData.anat);
    template = mean(allAnat, 4, 'omitnan');

    % 迭代：每轮将各受试者 anat 非线性配准到当前模板，再更新模板
    for it = 1:cfg.normalization.templateIters
        warpedAll = zeros([size(template), numel(subjData)], 'single');
        for i = 1:numel(subjData)
            moving = subjData(i).anat;
            [D, movingReg] = imregdemons(moving, template, ...
                cfg.normalization.demonsIters, ...
                'AccumulatedFieldSmoothing', cfg.normalization.demonsSmoothing);
            warpedAll(:, :, :, i) = movingReg;
            subjData(i).normField = D; %#ok<AGROW>
        end
        template = mean(warpedAll, 4, 'omitnan');
        fprintf('[group] template iteration %d done.\n', it);
    end

    tpl.volume = template;
    tpl.info = subjData(1).anatInfo;
    tpl.path = fullfile(cfg.paths.templateDir, 'group_template.nii');
    write_nifti_like(template, tpl.info, tpl.path);
end

% ============================== 标准化和平滑 ==============================
function out = normalize_subject_to_template(in, tpl, cfg)
    fprintf('[%s] normalize + smooth...\n', in.subjID);

    [D, ~] = imregdemons(in.anat, tpl.volume, ...
        cfg.normalization.demonsIters, ...
        'AccumulatedFieldSmoothing', cfg.normalization.demonsSmoothing);

    V = size(in.func4d, 4);
    funcNorm = zeros(size(in.func4d), 'single');
    for t = 1:V
        funcNorm(:, :, :, t) = imwarp(in.func4d(:, :, :, t), D, 'cubic', ...
            'OutputView', imref3d(size(tpl.volume)));
    end

    sigma = fwhm_to_sigma(cfg.preproc.smoothFWHM, cfg.preproc.voxelSize);
    funcSmooth = zeros(size(funcNorm), 'single');
    for t = 1:V
        funcSmooth(:, :, :, t) = imgaussfilt3(funcNorm(:, :, :, t), sigma);
    end

    out = in;
    out.normField = D;
    out.funcNorm = funcNorm;
    out.funcSmooth = funcSmooth;

    outDir = fullfile(in.subjDer, 'func_norm');
    ensure_dir(outDir);
    write_nifti_4d(funcSmooth, tpl.info, fullfile(outDir, 'func_smooth_norm_4d.nii'));
end

% ============================== 一级统计分析 ==============================
function glm = first_level_glm(subj, cfg)
    fprintf('[%s] first-level GLM...\n', subj.subjID);

    condPath = fullfile(cfg.paths.onsetDir, subj.subjID, cfg.firstLevel.conditionFileName);
    assert(exist(condPath, 'file') == 2, '缺少条件文件: %s', condPath);
    C = load(condPath);
    assert(isfield(C, 'names') && isfield(C, 'onsets') && isfield(C, 'durations'), ...
        'conditions.mat 必须包含 names/onsets/durations');

    Y4d = subj.funcSmooth;
    T = size(Y4d, 4);
    Y = reshape(Y4d, [], T)';         % T x V

    % --- 构建设计矩阵 ---
    Xtask = build_task_regressors(C.names, C.onsets, C.durations, T, cfg.firstLevel);
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

% ============================== I/O 与工具函数 ============================
function [V, info] = read_single_volume(inputDir, outDir)
    ensure_dir(outDir);
    nii = dir(fullfile(inputDir, '*.nii*'));
    if ~isempty(nii)
        p = fullfile(nii(1).folder, nii(1).name);
        info = niftiinfo(p);
        V = single(niftiread(info));
        if ndims(V) == 4, V = V(:, :, :, 1); end
    else
        dcm = dir(fullfile(inputDir, '**', '*.dcm'));
        assert(~isempty(dcm), '结构像目录中无 NIfTI 或 DICOM。');
        [V, info] = dicom_series_to_volume(dcm);
    end
    write_nifti_like(V, info, fullfile(outDir, 'anat_raw.nii'));
end

function [runs, infos] = read_functional_runs(funcDir, outDir)
    ensure_dir(outDir);
    runs = {};
    infos = {};

    niiTop = dir(fullfile(funcDir, '*.nii*'));
    if ~isempty(niiTop)
        p = fullfile(niiTop(1).folder, niiTop(1).name);
        info = niftiinfo(p);
        V = single(niftiread(info));
        assert(ndims(V) == 4, '功能像 NIfTI 必须为 4D。');
        runs{1} = V; infos{1} = info; %#ok<AGROW>
        write_nifti_4d(V, info, fullfile(outDir, 'run01_raw.nii'));
        return;
    end

    d = dir(funcDir);
    d = d([d.isdir]);
    d = d(~startsWith({d.name}, '.'));
    for i = 1:numel(d)
        runDir = fullfile(funcDir, d(i).name);
        nii = dir(fullfile(runDir, '*.nii*'));
        if ~isempty(nii)
            p = fullfile(nii(1).folder, nii(1).name);
            info = niftiinfo(p);
            V = single(niftiread(info));
            assert(ndims(V) == 4, 'run %s 需为4D。', d(i).name);
        else
            dcm = dir(fullfile(runDir, '**', '*.dcm'));
            [V, info] = dicom_series_to_4d(dcm);
        end
        runs{end + 1} = V; %#ok<AGROW>
        infos{end + 1} = info; %#ok<AGROW>
        write_nifti_4d(V, info, fullfile(outDir, sprintf('run%02d_raw.nii', i)));
    end
end

function V = drop_initial_volumes(V, n)
    assert(size(V, 4) > n, '时间点不足以去除前 %d 个。', n);
    V = V(:, :, :, n+1:end);
end

function Vout = slice_timing_correction(Vin, P)
    nx = size(Vin,1); ny = size(Vin,2); nz = size(Vin,3); nt = size(Vin,4);
    t = (0:nt-1) * P.TR;
    dt = P.TR / P.nslices;
    shifts = zeros(1, nz);
    for z = 1:nz
        ord = find(P.sliceOrder == z, 1, 'first');
        shifts(z) = (ord - find(P.sliceOrder == P.refSlice, 1, 'first')) * dt;
    end

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

function seg = segment_t1(anat, P)
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
    [~, order] = sort(c, 'ascend');
    mu = zeros(3,1); sd = zeros(3,1);
    for k = 1:3
        v = x(idx==order(k));
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
end

function Xtask = build_task_regressors(names, onsets, durations, T, P)
    t = (0:T-1)' * P.TR;
    h = canonical_hrf(P.TR, P.hrf);
    Xtask = zeros(T, numel(names));
    for i = 1:numel(names)
        u = zeros(T, 1);
        on = onsets{i};
        du = durations{i};
        if isscalar(du), du = repmat(du, size(on)); end
        for k = 1:numel(on)
            u = u + double(t >= on(k) & t < (on(k) + du(k)));
        end
        xc = conv(u, h);
        Xtask(:, i) = xc(1:T);
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

function h = canonical_hrf(TR, H)
    % 双Gamma HRF:
    % p1/d1: 主峰 shape/scale; p2/d2: undershoot shape/scale; ratio: 相对幅值比
    t = (0:TR:H.length)';
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
    tThr = tinv(1 - pVoxel, dof);
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
    tx = A(4,1); ty = A(4,2); tz = A(4,3);
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
    % FWHM = 2*sqrt(2*ln(2))*sigma ≈ 2.3548*sigma
    sigma = (fwhmMM ./ max(voxMM, eps)) / 2.3548;
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

function ensure_dir(p)
    if ~isfolder(p)
        mkdir(p);
    end
end

function write_nifti_like(vol, refInfo, pathOut)
    ensure_dir(fileparts(pathOut));
    info = refInfo;
    if isfield(info, 'Filename'), info = rmfield(info, 'Filename'); end
    if isfield(info, 'Filemoddate'), info = rmfield(info, 'Filemoddate'); end
    if isfield(info, 'Filesize'), info = rmfield(info, 'Filesize'); end
    info.ImageSize = size(vol);
    if numel(info.ImageSize)==2, info.ImageSize(3)=1; end
    niftiwrite(single(vol), pathOut, info, 'Compressed', false);
end

function write_nifti_4d(vol4d, refInfo, pathOut)
    ensure_dir(fileparts(pathOut));
    info = refInfo;
    if isfield(info, 'Filename'), info = rmfield(info, 'Filename'); end
    if isfield(info, 'Filemoddate'), info = rmfield(info, 'Filemoddate'); end
    if isfield(info, 'Filesize'), info = rmfield(info, 'Filesize'); end
    info.ImageSize = size(vol4d);
    niftiwrite(single(vol4d), pathOut, info, 'Compressed', false);
end

function [V, info] = dicom_series_to_volume(dcmList)
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
    info = niftiinfo_from_dicom(files{1}, size(V));
end

function [V4d, info] = dicom_series_to_4d(dcmList)
    assert(~isempty(dcmList), 'DICOM 列表为空。');
    files = fullfile({dcmList.folder}, {dcmList.name});
    hdr = cell(numel(files), 1);
    for i = 1:numel(files), hdr{i} = dicominfo(files{i}); end

    inst = zeros(numel(files),1); time = zeros(numel(files),1);
    for i = 1:numel(files)
        h = hdr{i};
        if isfield(h, 'InstanceNumber'), inst(i)=h.InstanceNumber; else, inst(i)=i; end
        if isfield(h, 'TemporalPositionIdentifier'), time(i)=h.TemporalPositionIdentifier; else, time(i)=1; end
    end

    tVals = unique(time);
    nT = numel(tVals);
    idxT1 = find(time == tVals(1));
    [~, ordZ] = sort(inst(idxT1));
    nZ = numel(idxT1);
    s = dicomread(files{idxT1(ordZ(1))});
    V4d = zeros([size(s), nZ, nT], 'single');

    for t = 1:nT
        idx = find(time == tVals(t));
        [~, o] = sort(inst(idx));
        idx = idx(o);
        for z = 1:nZ
            V4d(:,:,z,t) = single(dicomread(files{idx(z)}));
        end
    end
    info = niftiinfo_from_dicom(files{idxT1(ordZ(1))}, size(V4d));
end

function info = niftiinfo_from_dicom(oneFile, imgSize)
    h = dicominfo(oneFile);
    info = struct();
    info.ImageSize = imgSize;
    info.Datatype = 'single';
    info.SpaceUnits = 'Millimeter';
    info.TimeUnits = 'Second';
    if isfield(h,'PixelSpacing')
        px = double(h.PixelSpacing(:)');
    else
        px = [1 1];
    end
    if isfield(h,'SliceThickness'), st = double(h.SliceThickness); else, st = 1; end
    info.PixelDimensions = [px st 1];
end
