# MRILAB3

## Standalone 任务态 fMRI-BOLD 1:1 逻辑复现（MATLAB 2025a）

本仓库当前实现目标：

- **不依赖 DPABI/SPM 的预处理与统计函数调用**
- 以代码显式展示每一步计算逻辑（原始数据→预处理→GLM→阈值→可视化）
- 兼容你使用的 **SPM25 + MATLAB 2025a** 环境（即使安装了 SPM25，本流程也不调用其处理API）
- 输出更现代化图形：切面叠加、3D 表面渲染、可选体绘制（`volshow`）

核心文件：

- `./run_task_fmri_pipeline.m`
- `./task_fmri_pipeline_config.m`
- `./read_ctex_pdf_text.m`

兼容性补充（开箱即用）：
- DICOM 输入会先被转换为标准 NIfTI，再通过 `niftiinfo/niftiread` 回读进入后续流程，避免手工拼装头信息带来的不一致；
- 脚本写 NIfTI 时会自动过滤不同 MATLAB 版本中 `niftiwrite` 不识别的元数据字段（如 `Description`），避免因头信息字段差异直接报错；
- 写出前会强制对齐 `Info.Datatype` 与实际写入数组类型（`single`），并在必要时自动回退为无 `Info` 写出，避免因头信息不一致中断流程；
- 如需完全复现 README 中默认参数与图形表现，建议使用 MATLAB 2025a；
- 若使用其他版本 MATLAB，无需额外手工配置即可运行本流程，但可视化与个别数值细节可能存在轻微版本差异。

---

## 数据目录（当前默认读取方式）

```text
BOLDDATA/
  FunRaw/
    Sub001/
      *.nii(4D) 或 DICOM
      run1/ (可选) -> *.nii(4D) 或 DICOM
      run2/ (可选) -> *.nii(4D) 或 DICOM
    Sub002/
      ...
  T1Raw/
    Sub001/
      *.nii(3D) 或 DICOM
    Sub002/
      ...
```

兼容说明：
- 若 `BOLDDATA/FunRaw` 与 `BOLDDATA/T1Raw` 同时存在，流程优先按上述新结构读取；
- 若新结构不存在，则自动回退到旧版 `DataRaw/SubXXX/anat|func` 兼容模式。

---

## 流程模块（代码内有对应分块注释）

1. 原始数据读取（NIfTI/DICOM）
2. 去前若干时间点
3. Slice Timing（显式插值）
4. Motion Realign（3D rigid 配准）
5. T1→功能均值配准（MI）
6. T1 分割（bias 校正 + 3类概率映射）
7. 群体模板迭代构建 + 非线性标准化（demons）
8. 高斯平滑
9. 一级 GLM（HRF卷积、HPF、AR(1) 预白化、GLS/t-contrast）
10. 双阈值激活图（voxel + cluster）
11. 现代可视化导出（2D/3D/体绘制）

模板逻辑（standalone，不调用 SPM/DPABI API）：
- 分割先验模板：GM/WM/CSF 的 TPM（默认体积索引 `[1 2 3]`，即 `1→GM, 2→WM, 3→CSF`，并与对应组织概率直接融合）
  - 若你的 TPM 组织顺序不同，请调整 `cfg.templates.segmentation.tpmVolumeIndices`
- 群体模板：由本数据集 T1 在仓库内迭代构建（`Derivatives/Template/group_template.nii`）
- 标准空间目标模板（可选）：自动/手动指定 MNI 类模板后，执行“群体模板 -> 目标模板”二段映射
- East Asian 优先：自动发现多个 TPM 候选时，优先匹配 East Asian/Chinese 相关命名

---

## 运行方式

在 MATLAB 2025a 中：

```matlab
run_task_fmri_pipeline
```

通过 Python `matlab.engine` 测试（推荐先 smoke）：

```bash
python /absolute/path/to/MRILAB3/test_pipeline_with_matlab_engine.py --mode smoke
python /absolute/path/to/MRILAB3/test_pipeline_with_matlab_engine.py --mode full
```

说明：
- `smoke`：启动 MATLAB、加载仓库路径、读取 `task_fmri_pipeline_config` 并打印模板配置；
- `full`：直接执行 `run_task_fmri_pipeline` 全流程；
- 若缺少 `matlab.engine`，请在 MATLAB 安装目录执行 Python 引擎安装（`<matlabroot>/extern/engines/python`）。

参数修改位置：

```matlab
./task_fmri_pipeline_config.m
```

一级模型（SPM 风格）关键参数已全部内置在配置中（不依赖外部条件文件）：

- `cfg.firstLevel.timingUnits`（`'scans'` 或 `'secs'`）
- `cfg.firstLevel.scanOnsetIndexBase`（仅 `timingUnits='scans'`：0-based 或 1-based）
- `cfg.firstLevel.TR`
- `cfg.firstLevel.microtimeResolution`（SPM `fmri_t`）
- `cfg.firstLevel.microtimeOnset`（SPM `fmri_t0`）
- `cfg.firstLevel.design.names/onsets/durations`

DPABI 预处理步骤（白盒映射）也全部内置在配置中，并可逐项开关：

- `cfg.preproc.pipeline.removeFirstTimePoints.enabled`
- `cfg.preproc.pipeline.sliceTiming.enabled`
- `cfg.preproc.pipeline.realign.enabled`
- `cfg.preproc.pipeline.coregister.enabled`
- `cfg.preproc.pipeline.segment.enabled`
- `cfg.preproc.pipeline.normalize.enabled`
- `cfg.preproc.pipeline.smooth.enabled`

模板管理关键参数（`task_fmri_pipeline_config.m`）：
- `cfg.templates.searchDirs`
- `cfg.templates.installRoots`
- `cfg.templates.installRootEnvVars`
- `cfg.templates.segmentation.tpmPath`
- `cfg.templates.segmentation.tpmVolumeIndices`
- `cfg.templates.segmentation.preferEastAsian`
- `cfg.templates.segmentation.priorWeight`
- `cfg.templates.normalize.targetTemplatePath`
- `cfg.templates.normalize.preferMNI`

模板路径推断（适配 SPM/DPABI 常见组织）：
- 默认外部安装根目录：`D:\spm`、`D:\DPABI_V9.0_250415\DPARSF`
- 会自动扩展搜索常见子目录：`tpm`、`canonical`、`Templates` 以及上级目录中的 `Templates`
- 支持环境变量覆盖：`SPM_DIR`、`DPABI_DIR`、`DPARSF_DIR`

中文 ctex PDF（如 `MRILAB3.pdf`）建议先提取文本再对照实现：

```matlab
out = read_ctex_pdf_text('./MRILAB3.pdf', './MRILAB3_extracted.txt');
disp(out.method)
```

输出目录：

```text
./Derivatives/
```

每个被试可见：

- `func_preproc/motion_6dof.mat`
- `func_norm/func_smooth_norm_4d.nii`
- `first_level/tmap.nii`
- `first_level/mask_thresholded.nii`
- `visualization/modern_slice_overlay.png`
- `visualization/modern_3d_surface.png`
- `visualization/modern_volshow.png`（若支持）
- `visualization/activation_peaks.csv`

---

## 说明

- 本实现强调“**可学习、可追踪、可独立复现**”；
- 复杂步骤（例如 DARTEL/GRF）在此以独立可读的 MATLAB 近似实现替代；
- 若你希望下一步进一步逼近你课件中的每个数学细节（例如更严格 ReML/GRF/FWE），可以在当前脚本上继续逐块替换实现。

---

## Slice Timing：支持逐层毫秒输入（timing_ms）

在 `task_fmri_pipeline_config.m` 中可切换：

```matlab
cfg.preproc.pipeline.sliceTiming.mode = 'timing_ms';
cfg.preproc.pipeline.sliceTiming.timingMs = [ ...
    0 1430 880 330 1760 1210 660 110 1540 990 440 1870 1320 770 220 ...
    1650 1100 550 0 1430 880 330 1760 1210 660 110 1540 990 440 1870 ...
    1320 770 220 1650 1100 550];
```

说明：
- `timingMs` 长度需等于 `nslices`
- 如不设置 `refTimingMs`，默认使用 `timingMs(refSlice)` 作为参考层时间
- `timing_ms` 模式可兼容同一时间采集多个切片（如多带）
