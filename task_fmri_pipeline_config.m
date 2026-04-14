function cfg = task_fmri_pipeline_config()
% TASK_FMRI_PIPELINE_CONFIG
% -------------------------------------------------------------------------
% 在此处集中维护全流程配置。
% 建议仅在本文件修改参数，主脚本无需改动。
% -------------------------------------------------------------------------

    repoRoot = fileparts(mfilename('fullpath'));

    % ============================ 路径配置块 ============================
    cfg.paths.dataRawDir = fullfile(repoRoot, 'DataRaw');
    cfg.paths.derivativeDir = fullfile(repoRoot, 'Derivatives');
    cfg.paths.dartelDir = fullfile(cfg.paths.derivativeDir, 'DARTEL_Group');
    cfg.paths.onsetDir = fullfile(repoRoot, 'Onsets');

    % DPABI 常见目录命名候选（按优先级）
    cfg.paths.anatDirCandidates = {'anat', 'T1Img', 'T1', 'Anat'};
    cfg.paths.funcDirCandidates = {'func', 'FunImg', 'fun', 'Func'};

    % ========================== 预处理参数块 ============================
    % 基础采集参数
    cfg.preproc.TR = 2.0;
    cfg.preproc.removeFirstN = 6;
    cfg.preproc.nslices = 36;
    cfg.preproc.sliceOrder = [1:2:35, 2:2:36];
    cfg.preproc.refSlice = 18;

    % 文件前缀
    cfg.preproc.sliceTimingPrefix = 'a';
    cfg.preproc.realignPrefix = 'r';
    cfg.preproc.coregPrefix = 'r';

    % Realign 参数
    cfg.preproc.realign.quality = 0.9;
    cfg.preproc.realign.sep = 4;
    cfg.preproc.realign.fwhm = 5;
    cfg.preproc.realign.rtm = 1;
    cfg.preproc.realign.interp = 2;
    cfg.preproc.realign.wrap = [0 0 0];
    cfg.preproc.realign.which = [2 1];
    cfg.preproc.realign.writeInterp = 4;

    % Coreg 参数
    cfg.preproc.coreg.cost_fun = 'nmi';
    cfg.preproc.coreg.sep = [4 2];
    cfg.preproc.coreg.tol = ...
        [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
    cfg.preproc.coreg.fwhm = [7 7];
    cfg.preproc.coreg.interp = 4;

    % New Segment 参数（参考 SPM/DPABI 常用配置）
    cfg.preproc.segment.biasreg = 0.001;
    cfg.preproc.segment.biasfwhm = 60;
    cfg.preproc.segment.ngaus = [1 1 2 3 4 2];
    cfg.preproc.segment.native = {
        [1 1]   % GM: native + dartel imported
        [1 1]   % WM: native + dartel imported
        [1 0]   % CSF
        [0 0]
        [0 0]
        [0 0]
    };
    cfg.preproc.segment.warp.mrf = 1;
    cfg.preproc.segment.warp.cleanup = 1;
    cfg.preproc.segment.warp.reg = [0 0.001 0.5 0.05 0.2];
    cfg.preproc.segment.warp.affreg = 'eastern';  % East Asian
    cfg.preproc.segment.warp.fwhm = 0;
    cfg.preproc.segment.warp.samp = 3;

    % DARTEL 参数
    cfg.preproc.dartelTemplateName = 'Template';
    cfg.preproc.dartelParam(1).its = 3;
    cfg.preproc.dartelParam(1).rparam = [4 2 1e-06];
    cfg.preproc.dartelParam(1).K = 0;
    cfg.preproc.dartelParam(1).slam = 16;
    cfg.preproc.dartelParam(2).its = 3;
    cfg.preproc.dartelParam(2).rparam = [2 1 1e-06];
    cfg.preproc.dartelParam(2).K = 0;
    cfg.preproc.dartelParam(2).slam = 8;
    cfg.preproc.dartelParam(3).its = 3;
    cfg.preproc.dartelParam(3).rparam = [1 0.5 1e-06];
    cfg.preproc.dartelParam(3).K = 1;
    cfg.preproc.dartelParam(3).slam = 4;
    cfg.preproc.dartelParam(4).its = 3;
    cfg.preproc.dartelParam(4).rparam = [0.5 0.25 1e-06];
    cfg.preproc.dartelParam(4).K = 2;
    cfg.preproc.dartelParam(4).slam = 2;
    cfg.preproc.dartelParam(5).its = 3;
    cfg.preproc.dartelParam(5).rparam = [0.25 0.125 1e-06];
    cfg.preproc.dartelParam(5).K = 4;
    cfg.preproc.dartelParam(5).slam = 1;
    cfg.preproc.dartelParam(6).its = 3;
    cfg.preproc.dartelParam(6).rparam = [0.25 0.125 1e-06];
    cfg.preproc.dartelParam(6).K = 6;
    cfg.preproc.dartelParam(6).slam = 0.5;
    cfg.preproc.dartelOptim.lmreg = 0.01;
    cfg.preproc.dartelOptim.cyc = 3;
    cfg.preproc.dartelOptim.its = 3;
    cfg.preproc.dartelLastTemplate = 'Template_6.nii';

    % 标准化/平滑参数（对应 PDF 所述常用设置）
    cfg.preproc.normalizeBoundingBox = [-90 -126 -72; 90 90 108];
    cfg.preproc.normalizeVoxelSize = [3 3 3];
    cfg.preproc.smoothFWHM = [6 6 6];

    % ========================== 一级统计参数块 ==========================
    cfg.firstLevel.conditionFileName = 'conditions.mat';
    cfg.firstLevel.units = 'secs';
    cfg.firstLevel.TR = cfg.preproc.TR;
    cfg.firstLevel.fmri_t = 16;
    cfg.firstLevel.fmri_t0 = 8;
    cfg.firstLevel.hpf = 128;
    cfg.firstLevel.hrfDerivs = [0 0];
    cfg.firstLevel.cvi = 'AR(1)';
    cfg.firstLevel.includeMotionRegressors = true;

    % 对比设置（示例：条件1 > 条件2）
    cfg.firstLevel.contrasts(1).name = 'ConditionA > ConditionB';
    cfg.firstLevel.contrasts(1).weights = [1 -1];
    cfg.firstLevel.contrasts(2).name = 'ConditionB > ConditionA';
    cfg.firstLevel.contrasts(2).weights = [-1 1];

    % ========================== 结果显示参数块 ==========================
    cfg.results.contrastIndex = 1;
    cfg.results.voxelPThreshold = 0.001;
    cfg.results.clusterExtent = 20;
    cfg.results.thresDesc = 'none';     % 可选: 'none'/'FWE'/'FDR'
    cfg.results.title = 'Task Activation (MNI)';
    cfg.results.backgroundImage = '';
end
