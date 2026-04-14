function cfg = task_fmri_pipeline_config()
% TASK_FMRI_PIPELINE_CONFIG
% -------------------------------------------------------------------------
% Standalone pipeline 参数文件（MATLAB 2025a）
% 所有参数按功能分块集中，便于学习与定位。
% 一参数一注释：每个参数均说明“含义/用途/影响”。
% -------------------------------------------------------------------------

    repoRoot = fileparts(mfilename('fullpath')); % 仓库根目录（用于拼接默认输入输出路径）；用途：避免硬编码绝对路径；影响：路径错误会导致全流程无法读写数据。

    % ================================ 路径块 ================================
    cfg.paths.boldDataDir = fullfile(repoRoot, 'BOLDDATA'); % 原始数据总根目录；用途：统一管理原始输入；影响：改动后会同时改变结构像与功能像的搜索基准。
    cfg.paths.funcRawDir = fullfile(cfg.paths.boldDataDir, 'FunRaw'); % 功能像根目录（SPI/BOLD）；用途：读取 task-fMRI 时间序列；影响：目录错误会导致功能像缺失、无法进行预处理与GLM。
    cfg.paths.t1RawDir = fullfile(cfg.paths.boldDataDir, 'T1Raw'); % 结构像根目录（T1）；用途：读取解剖像用于配准/分割/归一化；影响：目录错误会导致配准和标准化失败。
    cfg.paths.dataRawDir = fullfile(repoRoot, 'DataRaw'); % 旧版兼容目录（仅当 FunRaw/T1Raw 不存在时回退）；用途：兼容历史组织；影响：保证旧数据仍可运行但优先使用新结构。
    cfg.paths.derivativeDir = fullfile(repoRoot, 'Derivatives'); % 输出根目录；用途：保存全部中间结果与最终结果；影响：决定输出位置与覆盖范围。
    cfg.paths.templateDir = fullfile(cfg.paths.derivativeDir, 'Template'); % 群体模板输出目录；用途：保存归一化模板；影响：影响模板复用与后续标准化读取路径。
    cfg.paths.anatDirCandidates = {'anat', 'T1Img', 'T1', 'Anat'}; % 旧版结构像子目录候选名；用途：兼容历史 DataRaw/SubXX 下的命名差异；影响：命中顺序不同会改变实际读取源。
    cfg.paths.funcDirCandidates = {'func', 'FunImg', 'fun', 'Func'}; % 旧版功能像子目录候选名；用途：兼容历史 DataRaw/SubXX 下的命名差异；影响：命中顺序不同会改变实际读取源。

    % ============================== 预处理参数块 =============================
    cfg.preproc.TR = 2.0; % 重复时间 TR（秒）；用途：时间轴换算、切片时序与GLM建模；影响：TR错误会导致时序校正与统计模型时间信息偏差。

    % DPABI 风格预处理模块白盒配置（步骤开关 + 各步骤参数）
    cfg.preproc.pipeline.removeFirstTimePoints.enabled = true; % 去前若干时间点开关；用途：丢弃磁化未稳态体积；影响：开启可提升稳定性，但会减少可用时间点。
    cfg.preproc.pipeline.removeFirstTimePoints.nVolumes = 6; % 去除体积数；用途：控制丢弃长度；影响：过小可能残留不稳态，过大降低统计效能。

    cfg.preproc.pipeline.sliceTiming.enabled = true; % 切片时序校正开关；用途：校正同一TR内不同切片采集时差；影响：开启可改善时序对齐，不匹配参数会引入伪差。
    cfg.preproc.pipeline.sliceTiming.nslices = 36; % 每体积切片数；用途：定义切片时序长度；影响：与真实切片数不一致会导致校正错误。
    cfg.preproc.pipeline.sliceTiming.sliceOrder = [1:2:35, 2:2:36]; % 切片采集顺序（order模式）；用途：计算每层时间偏移；影响：顺序错误会系统性扭曲时间相位。
    cfg.preproc.pipeline.sliceTiming.refSlice = 18; % 参考切片索引（order模式）；用途：统一插值对齐目标；影响：改变会轻微影响事件对齐与统计敏感性。
    % 1) 'order'     : 使用 sliceOrder + refSlice（传统方式）
    % 2) 'timing_ms' : 使用每层采集时间（毫秒），可兼容多带/同时间采集切片
    cfg.preproc.pipeline.sliceTiming.mode = 'order'; % 切片时序模式（order/timing_ms）；用途：选择顺序驱动或毫秒驱动校正；影响：应与采集协议一致，否则校正失真。
    cfg.preproc.pipeline.sliceTiming.timingMs = []; % 每层采集时间毫秒数组（timing_ms模式）；用途：支持多带等非等间隔采集；影响：数值不准会直接影响插值时序。
    cfg.preproc.pipeline.sliceTiming.refTimingMs = []; % 参考层毫秒时间（timing_ms模式）；用途：指定对齐参考时刻；影响：留空将用refSlice对应时间，改变会影响时间零点。
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
    cfg.preproc.pipeline.realign.enabled = true; % 头动校正开关；用途：对齐时间序列各体积；影响：关闭会保留头动伪影并显著影响统计。
    cfg.preproc.pipeline.realign.pyramidLevels = 3; % 头动配准金字塔层级；用途：控制多分辨率优化鲁棒性；影响：更高通常更稳但更慢。

    cfg.preproc.pipeline.coregister.enabled = true; % 结构像到功能像配准开关；用途：建立解剖-功能空间对应；影响：关闭会降低定位准确性。
    cfg.preproc.pipeline.coregister.pyramidLevels = 3; % 结构-功能配准金字塔层级；用途：提升跨模态配准收敛稳定性；影响：层级过低可能陷入局部最优。

    cfg.preproc.pipeline.segment.enabled = true; % 组织分割开关；用途：生成GM/WM/CSF概率图；影响：关闭会削弱后续标准化和可解释性。
    cfg.preproc.pipeline.segment.biasFieldSigma = 18; % 偏置场平滑尺度；用途：近似校正强度不均匀；影响：过小会过拟合噪声，过大可能欠校正。

    cfg.preproc.pipeline.normalize.enabled = true; % 标准化开关；用途：将个体数据映射到统一模板空间；影响：关闭无法进行跨被试空间对齐。
    cfg.preproc.pipeline.normalize.templateIters = 2; % 群体模板迭代次数；用途：迭代优化模板质量；影响：更多迭代通常提升模板稳定性但耗时增加。
    cfg.preproc.pipeline.normalize.demonsIters = [120 80 40]; % demons每层迭代次数；用途：控制非线性配准优化深度；影响：迭代不足配准欠佳，过多可能过拟合。
    cfg.preproc.pipeline.normalize.demonsSmoothing = 1.2; % demons位移场平滑系数；用途：约束形变平滑性；影响：更大更平滑但细节对齐减弱。

    cfg.preproc.pipeline.smooth.enabled = true; % 空间平滑开关；用途：提高信噪比并满足随机场近似；影响：关闭保细节但噪声更高。
    cfg.preproc.pipeline.smooth.voxelSize = [3 3 3]; % 体素尺寸(mm)；用途：将FWHM换算为高斯sigma；影响：与真实分辨率不一致会导致实际平滑核偏差。
    cfg.preproc.pipeline.smooth.fwhm = [6 6 6]; % 平滑核半高全宽(mm)；用途：控制空间平滑强度；影响：更大提升SNR但会降低空间精度。

    % ============================ 一级统计参数块 =============================
    cfg.firstLevel.TR = cfg.preproc.TR; % 一级模型TR（秒）；用途：设计矩阵采样与HRF离散化；影响：与预处理TR不一致会造成模型错配。
    cfg.firstLevel.hpf = 128; % 高通滤波截止周期（秒）；用途：去除低频漂移；影响：过短会滤掉真实慢效应，过长会保留漂移噪声。
    cfg.firstLevel.addQuadraticMotion = true; % 是否加入头动二次项；用途：增强头动混杂回归；影响：可降残余头动影响，但自由度减少。
    cfg.firstLevel.addLinearTrend = true; % 是否加入线性趋势项；用途：建模缓慢漂移；影响：开启可提高稳健性但轻微占用自由度。
    cfg.firstLevel.timingUnits = 'scans'; % 设计输入时间单位（scans/secs/seconds）；用途：解释onset/duration单位；影响：设置错误会导致事件时序整体偏移。
    cfg.firstLevel.scanOnsetIndexBase = 0; % scan单位下起始索引基准（0或1）；用途：兼容不同教材/软件记法；影响：0/1切换会整体平移一个TR。
    cfg.firstLevel.microtimeResolution = 16; % 微时间分辨率fmri_t；用途：HRF卷积的亚TR采样；影响：更高更精细但计算增加。
    cfg.firstLevel.microtimeOnset = 8; % 微时间参考点fmri_t0；用途：定义TR内参考采样位置；影响：改变会微调回归量相位。

    % 设计矩阵条件：全部在代码中定义（不依赖外部 conditions.mat）
    % 下述为示例默认值：1个条件 righthand，5个 onset，duration=15
    cfg.firstLevel.design.names = {'righthand'}; % 条件名称列表；用途：定义任务回归量标签；影响：名称顺序对应对比向量位置。
    cfg.firstLevel.design.onsets = {[0 30 60 90 120]}; % 条件起始时间列表；用途：生成任务设计矩阵；影响：onset偏差会直接降低激活检测准确性。
    cfg.firstLevel.design.durations = {15}; % 条件持续时间列表；用途：确定事件/区块持续建模；影响：持续时间不准会改变效应宽度与统计值。

    % HRF 双Gamma参数
    cfg.firstLevel.hrf.p1 = 6; % 主峰位置参数；用途：控制正峰出现时间；影响：改变会移动主激活响应峰值。
    cfg.firstLevel.hrf.d1 = 1; % 主峰尺度参数；用途：控制主峰宽度；影响：更大可使主峰更宽缓。
    cfg.firstLevel.hrf.p2 = 16; % 负峰位置参数；用途：控制负相位出现时间；影响：改变会影响HRF尾部形状。
    cfg.firstLevel.hrf.d2 = 1; % 负峰尺度参数；用途：控制负峰宽度；影响：影响去氧响应尾部拟合。
    cfg.firstLevel.hrf.ratio = 6; % 主峰/负峰比例；用途：设定双Gamma组合权重；影响：比例变化会改变整体响应对比度。
    cfg.firstLevel.hrf.length = 32; % HRF长度（秒）；用途：截断HRF卷积核；影响：过短可能截断尾部，过长增加计算但收益有限。

    % 一级对比（默认：righthand > implicit baseline；因仅1个任务条件，使用标量1）
    cfg.firstLevel.contrastWeights = 1; % 对比权重；用途：定义待检验线性假设；影响：符号/大小决定激活方向与统计解释。

    % ============================= 可视化参数块 ==============================
    cfg.results.voxelPThreshold = 0.001; % 体素水平阈值（未校正p）；用途：初筛显著体素；影响：阈值更严假阳性少但敏感性下降。
    cfg.results.clusterExtent = 20; % 团簇最小体素数；用途：团簇尺度二次筛选；影响：更大可抑制噪声小团簇但可能漏检小激活区。
    cfg.results.numPeaks = 30; % 峰值报告数量上限；用途：控制导出峰点表长度；影响：仅影响报告完整度，不改变统计结果本身。
    cfg.results.overlayAlpha = 0.50; % 叠加透明度；用途：控制激活图与背景解剖图可视平衡；影响：仅影响显示效果不影响统计。
    cfg.results.enableVolshow = true; % 体绘制开关；用途：启用/禁用3D体渲染输出；影响：仅影响可视化耗时与输出文件。
end
