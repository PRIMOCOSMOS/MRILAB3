function cfg = task_fmri_pipeline_config()
% TASK_FMRI_PIPELINE_CONFIG
% -------------------------------------------------------------------------
% Standalone pipeline 参数文件（MATLAB 2025a）
% 所有参数按功能分块集中，便于学习与定位。
% -------------------------------------------------------------------------

    repoRoot = fileparts(mfilename('fullpath'));

    % ================================ 路径块 ================================
    cfg.paths.dataRawDir = fullfile(repoRoot, 'DataRaw');
    cfg.paths.derivativeDir = fullfile(repoRoot, 'Derivatives');
    cfg.paths.onsetDir = fullfile(repoRoot, 'Onsets');
    cfg.paths.templateDir = fullfile(cfg.paths.derivativeDir, 'Template');

    % 输入目录候选（兼容 DPABI 习惯）
    cfg.paths.anatDirCandidates = {'anat', 'T1Img', 'T1', 'Anat'};
    cfg.paths.funcDirCandidates = {'func', 'FunImg', 'fun', 'Func'};

    % ============================== 预处理参数块 =============================
    cfg.preproc.TR = 2.0;
    cfg.preproc.removeFirstN = 6;
    cfg.preproc.nslices = 36;
    cfg.preproc.sliceOrder = [1:2:35, 2:2:36];
    cfg.preproc.refSlice = 18;
    assert(numel(cfg.preproc.sliceOrder) == cfg.preproc.nslices, ...
        'sliceOrder 长度必须与 nslices 一致。');

    % Realign/Coreg 参数
    cfg.preproc.realignPyramidLevels = 3;
    cfg.preproc.coregPyramidLevels = 3;
    cfg.preproc.biasFieldSigma = 18;           % 用于近似 bias field 校正

    % 标准化与平滑
    cfg.preproc.voxelSize = [3 3 3];           % mm
    cfg.preproc.smoothFWHM = [6 6 6];          % mm

    % ============================ 模板/标准化参数块 ==========================
    cfg.normalization.templateIters = 2;       % 群体模板迭代次数
    cfg.normalization.demonsIters = [120 80 40];
    cfg.normalization.demonsSmoothing = 1.2;

    % ============================ 一级统计参数块 =============================
    cfg.firstLevel.conditionFileName = 'conditions.mat';
    cfg.firstLevel.TR = cfg.preproc.TR;
    cfg.firstLevel.hpf = 128;                  % sec
    cfg.firstLevel.addQuadraticMotion = true;
    cfg.firstLevel.addLinearTrend = true;

    % HRF 双Gamma参数
    cfg.firstLevel.hrf.p1 = 6;
    cfg.firstLevel.hrf.d1 = 1;
    cfg.firstLevel.hrf.p2 = 16;
    cfg.firstLevel.hrf.d2 = 1;
    cfg.firstLevel.hrf.ratio = 6;
    cfg.firstLevel.hrf.length = 32;

    % 一级对比（示例：第1个条件 > 第2个条件）
    cfg.firstLevel.contrastWeights = [1 -1];

    % ============================= 可视化参数块 ==============================
    cfg.results.voxelPThreshold = 0.001;
    cfg.results.clusterExtent = 20;
    cfg.results.numPeaks = 30;
    cfg.results.overlayAlpha = 0.50;
    cfg.results.enableVolshow = true;
end
