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
    cfg.paths.templateDir = fullfile(cfg.paths.derivativeDir, 'Template');

    % 输入目录候选（兼容 DPABI 习惯）
    cfg.paths.anatDirCandidates = {'anat', 'T1Img', 'T1', 'Anat'};
    cfg.paths.funcDirCandidates = {'func', 'FunImg', 'fun', 'Func'};

    % ============================== 预处理参数块 =============================
    cfg.preproc.TR = 2.0;

    % DPABI 风格预处理模块白盒配置（步骤开关 + 各步骤参数）
    cfg.preproc.pipeline.removeFirstTimePoints.enabled = true;
    cfg.preproc.pipeline.removeFirstTimePoints.nVolumes = 6;

    cfg.preproc.pipeline.sliceTiming.enabled = true;
    cfg.preproc.pipeline.sliceTiming.nslices = 36;
    cfg.preproc.pipeline.sliceTiming.sliceOrder = [1:2:35, 2:2:36];
    cfg.preproc.pipeline.sliceTiming.refSlice = 18;
    % 1) 'order'     : 使用 sliceOrder + refSlice（传统方式）
    % 2) 'timing_ms' : 使用每层采集时间（毫秒），可兼容多带/同时间采集切片
    cfg.preproc.pipeline.sliceTiming.mode = 'order';
    cfg.preproc.pipeline.sliceTiming.timingMs = [];    % 长度需等于 nslices
    cfg.preproc.pipeline.sliceTiming.refTimingMs = []; % 留空则自动使用 timingMs(refSlice)
    assert(numel(cfg.preproc.pipeline.sliceTiming.sliceOrder) == cfg.preproc.pipeline.sliceTiming.nslices, ...
        'sliceOrder 长度必须与 nslices 一致。');

    % timing_ms 示例（36 层）
    % 说明：重复毫秒值表示这些切片在同一时刻采集（常见于多带采集）
    % cfg.preproc.pipeline.sliceTiming.mode = 'timing_ms';
    % cfg.preproc.pipeline.sliceTiming.timingMs = [ ...
    %     0 1430 880 330 1760 1210 660 110 1540 990 440 1870 1320 770 220 ...
    %     1650 1100 550 0 1430 880 330 1760 1210 660 110 1540 990 440 1870 ...
    %     1320 770 220 1650 1100 550];

    % Realign/Coreg 参数
    cfg.preproc.pipeline.realign.enabled = true;
    cfg.preproc.pipeline.realign.pyramidLevels = 3;

    cfg.preproc.pipeline.coregister.enabled = true;
    cfg.preproc.pipeline.coregister.pyramidLevels = 3;

    cfg.preproc.pipeline.segment.enabled = true;
    cfg.preproc.pipeline.segment.biasFieldSigma = 18;   % 用于近似 bias field 校正

    cfg.preproc.pipeline.normalize.enabled = true;
    cfg.preproc.pipeline.normalize.templateIters = 2;    % 群体模板迭代次数
    cfg.preproc.pipeline.normalize.demonsIters = [120 80 40];
    cfg.preproc.pipeline.normalize.demonsSmoothing = 1.2;

    cfg.preproc.pipeline.smooth.enabled = true;
    cfg.preproc.pipeline.smooth.voxelSize = [3 3 3];     % mm
    cfg.preproc.pipeline.smooth.fwhm = [6 6 6];          % mm

    % ============================ 一级统计参数块 =============================
    cfg.firstLevel.TR = cfg.preproc.TR;
    cfg.firstLevel.hpf = 128;                  % sec
    cfg.firstLevel.addQuadraticMotion = true;
    cfg.firstLevel.addLinearTrend = true;
    cfg.firstLevel.timingUnits = 'scans';      % 与 SPM "Units for design" 对齐: 'scans' 或 'secs'
    cfg.firstLevel.scanOnsetIndexBase = 0;     % 仅 timingUnits='scans' 时生效：0(首体积=0) 或 1(首体积=1)
    cfg.firstLevel.microtimeResolution = 16;   % SPM fmri_t
    cfg.firstLevel.microtimeOnset = 8;         % SPM fmri_t0

    % 设计矩阵条件：全部在代码中定义（不依赖外部 conditions.mat）
    % 下述为示例默认值：1个条件 righthand，5个 onset，duration=15
    cfg.firstLevel.design.names = {'righthand'};
    cfg.firstLevel.design.onsets = {[0 30 60 90 120]};   % 若 timingUnits='scans'，按 scanOnsetIndexBase 解释
    cfg.firstLevel.design.durations = {15};              % 标量表示该条件所有 trial 共享同一时长

    % HRF 双Gamma参数
    cfg.firstLevel.hrf.p1 = 6;
    cfg.firstLevel.hrf.d1 = 1;
    cfg.firstLevel.hrf.p2 = 16;
    cfg.firstLevel.hrf.d2 = 1;
    cfg.firstLevel.hrf.ratio = 6;
    cfg.firstLevel.hrf.length = 32;

    % 一级对比（默认：righthand > implicit baseline；因仅1个任务条件，使用标量1）
    cfg.firstLevel.contrastWeights = 1;

    % ============================= 可视化参数块 ==============================
    cfg.results.voxelPThreshold = 0.001;
    cfg.results.clusterExtent = 20;
    cfg.results.numPeaks = 30;
    cfg.results.overlayAlpha = 0.50;
    cfg.results.enableVolshow = true;
end
