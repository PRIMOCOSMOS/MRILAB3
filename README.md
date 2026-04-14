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

---

## 数据目录（参考 DPABI 约定）

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

`conditions.mat` 至少包含：

- `names`：条件名 cell，例如 `{'ConditionA','ConditionB'}`
- `onsets`：每个条件 onset（秒）cell
- `durations`：每个条件持续时间（秒）cell

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

---

## 运行方式

在 MATLAB 2025a 中：

```matlab
run_task_fmri_pipeline
```

参数修改位置：

```matlab
./task_fmri_pipeline_config.m
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
