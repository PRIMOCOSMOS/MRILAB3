# MRILAB3

## 任务态 fMRI-BOLD 全流程 MATLAB 复现

本仓库已提供基于 **SPM12（DPABI任务态处理核心思路）** 的完整脚本：

- `/home/runner/work/MRILAB3/MRILAB3/run_task_fmri_pipeline.m`
- `/home/runner/work/MRILAB3/MRILAB3/task_fmri_pipeline_config.m`

流程覆盖：

1. 原始数据读取（DICOM/NIfTI）
2. 预处理（去前时间点、Slice Timing、Realign、Coreg、New Segment + DARTEL、Normalize、Smooth）
3. 一级统计（Specify/Estimate/Contrast，AR(1)）
4. 结果输出（阈值化结果、MNI坐标表、正交立体视图与3D渲染图）

---

## 目录结构（参考 DPABI 习惯）

在仓库根目录下准备：

```text
DataRaw/
  Sub001/
    anat/ 或 T1Img/
      *.nii 或 DICOM
    func/ 或 FunImg/
      run1/ (可选) -> *.nii 或 DICOM
      run2/ (可选) -> *.nii 或 DICOM
  Sub002/
    ...

Onsets/
  Sub001/
    conditions.mat
  Sub002/
    conditions.mat
```

`conditions.mat` 至少包含 3 个变量（SPM 一级分析标准格式）：

- `names`：条件名 cell，例如 `{'ConditionA','ConditionB'}`
- `onsets`：每个条件 onset（秒）cell，例如 `{[10 40 70],[25 55 85]}`
- `durations`：每个条件持续时间（秒）cell，例如 `{[2 2 2],[2 2 2]}`

---

## 运行方法

1. 打开 MATLAB，加入 SPM12 路径。
2. 按需修改 `/home/runner/work/MRILAB3/MRILAB3/task_fmri_pipeline_config.m` 参数（TR、切片顺序、对比等）。
3. 在 MATLAB 执行：

```matlab
run('/home/runner/work/MRILAB3/MRILAB3/run_task_fmri_pipeline.m')
```

输出目录默认在：

- `/home/runner/work/MRILAB3/MRILAB3/Derivatives/`

其中每个被试会生成预处理结果、一级统计目录及结果图像（`orthview_activation.png`、`render_activation.png`）与阈值结果表（`xSPM_thresholded.mat`）。
