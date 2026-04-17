# Copilot Chat Conversation Export: fMRI BOLD data analysis pipeline development

**User:** @PRIMOCOSMOS
**Thread URL:** https://github.com/copilot/c/6b2d5eee-e355-4fd9-842a-f0d4a63fd639

## @PRIMOCOSMOS

我目前有一个库，目标功能是仿照DPABI中的DPARSF功能模块和SPM中的功能，自行实现一个standalone的任务态fMRI-BOLD的数据分析pipeline。我在上下文中已经勾选了我自己的库MRILAB3和两个提及的开源库。但是我自己的库由于未能完全解析处理的逻辑和具体方式，尚不能实现完整的功能。现在你的任务是根据我自己的库和这两个开源库的相关信息，撰写一份prompt给copilot agent，使其能理解开源库的详细的逻辑，并参照完成代码设计，注意并非调用。另外，prompt要以markdown文件的形式给出，我要直接导出。我会提供对pipeline原理的描述。如果开源2项目中有特别重要的参考源代码，请你也直接指出，我可以直接提供给copilot agent。使用工具箱的实验已经成功。我们也可以尝试使用代码自行实现，但前提是明确整个数据处理与分析流程的pipeline。

\paragraph{格式转换层：DICOM to NIFTI (EPI ， T1)}

SPM和基于它的DPABI所接受的文件结构就是结构像和功能像的按被试分类的文件夹。\\

通过MATLAB提供的内置的函数严格地解析DICOM的文件头，将其转换为NIFTI文件，通过仿射变换（Affine Transformation）将散乱的 2D 切片重组成一个四维矩阵。注意，根据预先的解析，EPI功能像是MOSAIC类型的，因此是将3D分层拼接为2D的数据结构，例如实际上fMRI是36层64×64的功能像，但实际上是384*384的2D数据结构。这一点是非常重要的。对于结构像，也是按照类似的方式处理。

\paragraph{去除起始不稳定信号时间点}
在 MRI 脉冲序列刚开始激发时，组织内部的自旋系统（Spin System）由于连续的射频脉冲（RF）激发，纵向磁化矢量（Longitudinal Magnetization, $M_z$)还没有达到动态平衡（Steady State）。前几个 TR 的图像信号通常偏高且不稳定（被称为 T1 弛豫的饱和过程）。
\textbf{底层代码就是一个简单的矩阵截取操作}，截取发生在时间维度。

\paragraph{时间校正：Slice Timing}

这个校正的数学原理和物理上的必要性已在上文中介绍，此处不再赘述，简单地来说就是根据傅里叶变换的时移性质来机型校正。我们设计的预期是给出每一层的相对采集时间（以ms为单位的数组，长度为层数），代码就能完成这一步。

\paragraph{头动校正：Realign}

人在扫描仪里不可避免地会有呼吸和轻微头动。在数学建模上，我们定义：头部被假设为一个刚体（Rigid Body），拥有 6 个自由度（DOF）。\\
fMRI 的 BOLD 信号极其微弱（变化幅度在 1 \% - 5 \% 左右），如果不校正头动，体素边缘的脑脊液或白质移位到灰质区域，会产生巨大的伪影信号，淹没真实的神经活动。因此，我们
使用\textbf{高斯-牛顿优化算法（Gauss-Newton Optimization）}与\textbf{泰勒一阶展开}来迭代求解这6 个参数——3个平移$T_x,T_y,T_z$、3个旋转角$\theta_x,\theta_y,\theta_z$，并最终输出 rp\_*.txt 头动参数文件，就像SPM做的那样。\\

\paragraph{重定位与配准准备：Reorient (Fun* ，T1*) + T1 Coreg to Fun}

将图像的坐标原点手动或自动平移到前连合（Anterior Commissure, AC）。仅仅是修改 NIfTI 的仿射矩阵参数，不改变原始像素值;将高分辨率的 T1 结构像空间对齐配准到低分辨率的 fMRI (Fun) 功能像空间。

\paragraph{分割与构建个体模板：New Segment + DARTEL}

为了进行群体级别的统计分析，我们必须把所有人的大脑拉伸、压缩到一个标准的公共空间（如 MNI 空间），这也就是所谓的\textbf{模板配准}。
DARTEL (Diffeomorphic Anatomical Registration Through Exponentiated Lie Algebra，李代数指数化的微分同胚解剖配准) 是 SPM 里面最高阶、最精确的非线性配准算法；East Asian 意味着算法在贝叶斯估计中使用“东亚人脑概率图”作为先验知识（Prior），这对提高中国受试者的分割精度很有帮助。


\paragraph{空间标准化：Normalize by DARTEL}
这是上一步的延续。将个体经过头动校正和切片时间校正后的 4D 功能数据，全部“扭曲”变形到 MNI (Montreal Neurological Institute) 标准空间。这样，配准才算完成。

\paragraph{空间平滑：Smooth}

原理简单，不再赘述，用于增加信噪比 、解剖差异容错、满足高斯随机场（GRF）理论，有利于后续分析。至此，我们应当得到一个SPM.mat文件，作为处理的结果。

\paragraph{SPM Specify 1st-level and estimate}

这就是一般线性模型（General Linear Model, GLM）的大致分析流程，也就是具体地，构造线性矩阵，估计参数$\bm \beta$，评估方差，并利用多重比较校正去除假阳性区域，得到结果。我们利用的是上文的统计量$T$，故称为T-contrast。之后，我们还能使用renderer进行更为形象的展示。

## Copilot

这份Prompt被设计为可以直接发给Copilot Agent，用以帮助你理解开源库逻辑并在 `PRIMOCOSMOS/MRILAB3` 库中实现这套fMRI-BOLD的预处理和分析pipeline。

在这份Prompt中，我已经提取了你在描述中强调的关键物理/数学原理（例如MOSAIC格式处理、高斯-牛顿优化、傅里叶时移性质、DARTEL等），并使用GitHub代码搜索定位了 `spm/spm` 和 `Chaogan-Yan/DPABI` 中**最核心的参考源码文件**，直接喂给Copilot以提高它生成精确API调用/复现代码的准确度。

请直接复制以下Markdown格式的内容导出：

```markdown
# fMRI BOLD Pipeline Development Prompt for Copilot Agent

## Context and Goal
You are an expert MATLAB developer and neuroimaging methods engineer. My goal is to build a standalone, task-state fMRI-BOLD data analysis pipeline from scratch in my pure MATLAB repository (`PRIMOCOSMOS/MRILAB3`). This pipeline should emulate the data processing and analysis logic found in two major neuroimaging toolboxes:
- **SPM (Statistical Parametric Mapping)** (`spm/spm`)
- **DPABI / DPARSF** (`Chaogan-Yan/DPABI`)

I have successfully used the GUI toolboxes, but my own repository currently lacks the complete logic to automate the workflow programmatically. Your task is to refer to the source code and logic of the open-source projects (`spm/spm` and `Chaogan-Yan/DPABI`) to help me write the pipeline code in `MRILAB3`. **Do not just blindly instruct me to open the GUI.** We want to implement the pipeline script that invokes the underlying algorithms or properly constructs the batch processing structures correctly, mimicking the automated batch logic of SPM/DPARSF.

## Pipeline Steps and Reference Logic

Please review the following steps, read the corresponding open-source codes, and write the appropriate MATLAB pipeline code for my repository.

### 1. Format Conversion: DICOM to NIFTI (EPI, T1)
**Description:** Parse DICOM headers to convert data to NIFTI files using Affine Transformations. Crucially, the fMRI EPI sequences are MOSAIC types (3D slices spliced into 2D data structures, e.g., 64x64 functional slices saved as 384x384 2D images).
**Reference Code to Read:** 
- In `spm/spm`: Look at `spm_dicom_convert.m`, `spm_dicom_headers.m`, and `spm_dicom_metadata.m`.
- In `Chaogan-Yan/DPABI`: Look at the logic in `y_Call_dcm2nii.m` or `y_dicom2nii.m`.

### 2. Remove Initial Unstable Signal Time Points
**Description:** Due to T1 relaxation saturation, the first few TRs are unstable. We need to discard these initial volumes.
**Logic:** Simple matrix truncation along the 4th dimension (time). Read the 4D NIFTI, slice `data(:,:,:, N+1:end)`, and save it back.

### 3. Slice Timing Correction
**Description:** Correct for slice acquisition delays using the Fourier transform's time-shift property. The function needs to accept an array of relative acquisition times (in ms) for each slice.
**Reference Code to Read:**
- In `spm/spm`: Look at `spm_slice_timing.m`.

### 4. Head Motion Correction (Realign)
**Description:** Correct respiratory and minor head movements. Assume the head is a rigid body (6 DOF). Use Gauss-Newton Optimization and 1st-order Taylor expansion to iteratively estimate 3 translations ($T_x,T_y,T_z$) and 3 rotations ($\theta_x,\theta_y,\theta_z$). The output must produce an `rp_*.txt` file identical to SPM's output.
**Reference Code to Read:**
- In `spm/spm`: Look at `spm_realign.m` and `spm_reslice.m`.

### 5. Reorient & T1 Coregistration to Functional Images (Coreg)
**Description:** Manually/Automatically shift the origin to the Anterior Commissure (AC) by modifying the NIfTI affine matrix without changing pixel values. Then, coregister the high-resolution T1 structural image to the low-resolution fMRI (Fun) space.
**Reference Code to Read:**
- In `spm/spm`: Look at `spm_coreg.m` (mutual information based coregistration).

### 6. Segmentation & Individual Template Construction (New Segment + DARTEL)
**Description:** Perform tissue segmentation and warp individual brains to a standard MNI space. DARTEL (Diffeomorphic Anatomical Registration Through Exponentiated Lie Algebra) is the core algorithm here. For our dataset, we use the "East Asian" probability map as the Bayesian prior.
**Reference Code to Read:**
- In `spm/spm`: Look at `spm_preproc_run.m` (New Segment) and `spm_dartel_warp.m`.

### 7. Spatial Normalization (Normalize by DARTEL)
**Description:** Apply the DARTEL flow fields to warp the 4D functional data (already motion and slice-timing corrected) to the MNI standard space.
**Reference Code to Read:**
- In `spm/spm`: Look at `spm_dartel_norm_fun.m`.

### 8. Spatial Smoothing (Smooth)
**Description:** Apply a 3D Gaussian kernel filter to increase SNR, accommodate anatomical variability, and meet Gaussian Random Field (GRF) theory assumptions.
**Reference Code to Read:**
- In `spm/spm`: Look at `spm_smooth.m`.

### 9. First-Level GLM Specification and Estimation
**Description:** Build the General Linear Model (GLM) matrix, estimate the $\bm \beta$ parameters, evaluate variance, and perform multiple comparison corrections. Use T-contrasts to compute the $T$ statistics to find activated regions, yielding a complete `SPM.mat` workspace.
**Reference Code to Read:**
- In `spm/spm`: Look at `spm_fMRI_design.m`, `spm_spm.m` (Estimation), and `spm_contrasts.m`.

## Your Task as Copilot
1. Cross-reference the logic in the provided `.m` files from the SPM and DPABI repositories.
2. Design a clean, modular MATLAB script/function structure for `PRIMOCOSMOS/MRILAB3` that loops through subjects and executes this exact pipeline headless (without GUI).
3. Write the MATLAB code implementing each step. Provide detailed comments referencing the mathematical principles provided in this prompt (e.g., Gauss-Newton, Fourier time-shift, DARTEL Lie Algebra).
4. Formulate the `matlabbatch` struct arrays correctly to automate SPM underlying algorithms programmatically, emulating the `DPARSFA_run.m` script logic.
```
