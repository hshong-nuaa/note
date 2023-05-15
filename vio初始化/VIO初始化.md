# VIO初始化

## 一、VINS

### 1. 滑窗下的纯视觉SfM

检查最新帧与之前帧之间特征的对应关系，如果能够找到稳定的特征跟踪与足够的视差，使用五点算法在两个frame之间恢复相对旋转与最大尺度平移。然后设置尺度并三角化两帧之间观测到的所有特征。基于这些三角化的点，采用PnP的方法来估计滑窗中所有帧。最后采用全局BA的方法最小化全局重投影误差。由于没有全局坐标系，采用相机第一帧作为SfM的参考。所有的帧的位姿$(\overline{p}^{c_0}_{c_k},\overline{R}^{c_0}_{c_k})$和特征位置都表示为$(\cdot)^{c_0}$，给定外参$(p^b_c,R^b_c)$，可以得到相机到IMU的转换
$$
\begin{aligned}
\mathbf{R}_{b_k}^{c_0} & =\left(\mathbf{R}_c^b\right)^{-1} \mathbf{R}_{c_k}^{c_0} \cdot\\
s \overline{\mathbf{p}}_{b_k}^{c_0} & =s \overline{\mathbf{p}}_{c_k}^{c_0}-\mathbf{R}_{b_k}^{c_0} \mathbf{p}_c^b
\end{aligned}
$$
其中，$s$是未知的缩放参数，将在下一节求解

### 2. 视觉惯性对齐

#### (1)陀螺零偏估计

考虑滑窗的两个连续帧$b_k,b_{k+1}$，可以通过SfM得到两帧的旋转$R^{C_0}_{b_k},R^{C_0}_{b_{k+1}}$，同时可以得到两帧之间的预积分约束$\vartriangle R^{b_k}_{b_{k+1}}$，将预积分项与陀螺零偏线性化，并使得下面残差最小：
$$
\begin{gathered}
\min _{\delta b_w} \sum_{k \in \mathcal{B}}\left\|\mathbf{R}_{b_{k+1}}^{c_0}{ }^{-1} \cdot \mathbf{R}_{b_k}^{c_0} \cdot \boldsymbol{\vartriangle R}_{b_{k+1}}^{b_k}\right\|^2 \\
\boldsymbol{\vartriangle R}_{b_{k+1}}^{b_k} \approx \vartriangle \hat{\boldsymbol{ R}}_{b_{k+1}}^{b_k} \cdot 
Exp(\mathbf{J}_{b_w} \delta \mathbf{b}_w)
\end{gathered}
$$

#### (2)速度，重力与尺度初始化

当完成陀螺零偏初始化后，将对其他量进行初始化。
$$
\mathcal{X}_I=\left[\mathbf{v}_{b_0}^{b_0}, \mathbf{v}_{b_1}^{b_1}, \ldots \mathbf{v}_{b_n}^{b_n}, \mathbf{g}^{c_0}, s\right]
$$
考虑两帧之间，可以获得以下等式：
$$
\begin{aligned}
\boldsymbol{\Delta p}_{b_{k+1}}^{b_k} & =\mathbf{R}_{c_0}^{b_k}\left(s\left(\overline{\mathbf{p}}_{b_{k+1}}^{c_0}-\overline{\mathbf{p}}_{b_k}^{c_0}\right)+\frac{1}{2} \mathbf{g}^{c_0} \Delta t_k^2-\mathbf{R}_{b_k}^{c_0} \mathbf{v}_{b_k}^{b_k} \Delta t_k\right) \\
\boldsymbol{\Delta v}_{b_{k+1}}^{b_k} & =\mathbf{R}_{c_0}^{b_k}\left(\mathbf{R}_{b_{k+1}}^{c_0} \mathbf{v}_{b_{k+1}}^{b_{k+1}}+\mathbf{g}^{c_0} \Delta t_k-\mathbf{R}_{b_k}^{c_0} \mathbf{v}_{b_k}^{b_k}\right) .
\end{aligned}
$$
联合(4)与(2)可以得到线性测量模型：
$$
\hat{\mathbf{z}}_{b_{k+1}}^{b_k}=\left[\begin{array}{c}
\Delta\hat{\boldsymbol{ p}}_{b_{k+1}}^{b_k}-\mathbf{p}_c^b+\mathbf{R}_{c_0}^{b_k} \mathbf{R}_{b_{k+1}}^{c_0} \mathbf{p}_c^b \\
\Delta\hat{\boldsymbol{ v}}_{b_{k+1}}^{b_k}
\end{array}\right]=\mathbf{H}_{b_{k+1}}^{b_k} \mathcal{X}_I+\mathbf{n}_{b_{k+1}}^{b_k}
$$
其中
$$
\mathbf{H}_{b_{k+1}}^{b_k}=\left[\begin{array}{cccc}
-\mathbf{I} \Delta t_k & \mathbf{0} & \frac{1}{2} \mathbf{R}_{c_0}^{b_k} \Delta t_k^2 \mathbf{R}_{c_0}^{b_k}\left(\overline{\mathbf{p}}_{c_{k+1}}^{c_0}-\overline{\mathbf{p}}_{c_k}^{c_0}\right) \\
-\mathbf{I} & \mathbf{R}_{c_0}^{b_k} \mathbf{R}_{b_{k+1}}^{c_0} & \mathbf{R}_{c_0}^{b_k} \Delta t_k & \mathbf{0}
\end{array}\right] .
$$
通过求解最小二乘问题：
$$
\min _{\mathcal{X}_I} \sum_{k \in \mathcal{B}}\left\|\hat{\mathbf{z}}_{b_{k+1}}^{b_k}-\mathbf{H}_{b_{k+1}}^{b_k} \mathcal{X}_I\right\|^2
$$
可以获得速度，第一帧坐标系下的重力与缩放因子。

#### (3)重力