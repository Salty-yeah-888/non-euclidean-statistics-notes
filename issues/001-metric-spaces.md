# 第 001 期：从度量空间到稳健 Fréchet 回归

**Non-Euclidean Statistics — Learning & Research Notes**

- 日期：2026-11-16
- 关键词：Metric Space、Fréchet Mean、Nonparametric Regression、Robust Estimation
- 今日精读：*Robust Local Fréchet Regression: Nonasymptotic Theory and Minimax Optimality*

> 本系列并行记录非欧统计学的基础知识和前沿研究。基础部分循序渐进，研究部分保留重要的数学公式、理论结论和方法创新。

---

## A. 基础知识：为什么统计学需要非欧空间？

### 1. 从方向数据开始

假设两个观测方向分别为 $10^\circ$ 和 $350^\circ$。

普通算术平均给出：

$$
\frac{10^\circ+350^\circ}{2}=180^\circ.
$$

但是两个方向都接近正北方，合理的平均方向应该是 $0^\circ$。

问题在于：方向数据属于圆周 $S^1$，而不是普通实数轴。角度 $0^\circ$ 与 $360^\circ$ 表示同一个点。

这说明，在分析数据前，我们需要先理解数据所处的空间。

### 2. 度量空间（Metric Space）

一个度量空间表示为：

$$
(\Omega,d),
$$

其中 $\Omega$ 是数据对象的集合，$d$ 是距离函数。

度量满足：

1. $d(x,y)\geq 0$，且 $d(x,y)=0$ 当且仅当 $x=y$；
2. $d(x,y)=d(y,x)$；
3. $d(x,z)\leq d(x,y)+d(y,z)$。

常见的非欧数据包括：

| 数据对象 | 典型空间 |
|---|---|
| 方向数据 | 圆周 $S^1$ |
| 地球表面坐标 | 球面 $S^2$ |
| 协方差矩阵 | SPD 矩阵空间 |
| 概率分布 | Wasserstein 空间 |
| 树和网络 | 树空间、图空间 |

注意：这些对象有时可以嵌入欧氏空间，但它们的自然统计距离不一定是欧氏距离。

### 3. Fréchet Mean

在欧氏空间中，样本均值可以通过最小化平方距离得到：

$$
\bar x
=
\arg\min_{\mu\in\mathbb R^p}
\sum_{i=1}^n\|x_i-\mu\|^2.
$$

这个形式不直接依赖于求和平均，可以推广到度量空间。

**定义：样本 Fréchet mean**

$$
\boxed{
\hat\mu_F\in
\arg\min_{\omega\in\Omega}
\frac1n\sum_{i=1}^n d^2(Y_i,\omega)
}
$$

它寻找一个与所有观测对象的平方距离总和最小的代表对象。

当 $\Omega=\mathbb R^p$ 且 $d$ 为欧氏距离时，Fréchet mean 就是普通样本均值。

对于更一般的空间，其存在性和唯一性需要额外研究。

**今日核心思想：**

非欧统计学并不是简单地把欧氏公式换成复杂的几何公式，而是重新思考均值、方差和回归等统计概念在一般空间中应该如何定义。

---

## B. 前沿研究雷达

以下六项为本期关注的 2026 年研究成果，均附有原始 arXiv 链接。它们属于预印本研究，不应直接视为已经完成期刊同行评审。

### 1. Adaptive Local Fréchet Regression with Polyspherical Predictors

- 作者：Chang Jun Im
- 日期：2026-11-09
- 方向：流形回归、变量选择
- 原文：https://arxiv.org/abs/2611.13240

**研究问题：** 当解释变量分布在多个球面组成的乘积空间中，如何进行局部 Fréchet 回归？

**核心贡献：** 提出针对不同球面因子的局部常数和局部线性估计器，允许各因子使用不同带宽，并通过交叉拟合预测风险进行变量选择。

**值得关注：** 建立统一收敛速率及变量选择一致性结果，为具有复杂几何结构的解释变量提供非参数回归方法。

### 2. Theoretical Foundations of Principal Manifold Estimation with Non-Euclidean Templates

- 作者：Kun Meng, Christopher Perez
- 首次提交：2026-04-05
- 方向：主流形估计、流形学习
- 原文：https://arxiv.org/abs/2604.04272

**研究问题：** 如何从高维观测点云中恢复潜在的低维非线性流形？

**核心贡献：** 为允许非欧拓扑结构的主流形估计建立数学框架，研究估计对象的良定义性、计算算法收敛性与有限样本一致性。

**值得关注：** 将流形学习从单纯的降维问题推进到具有明确统计估计目标和理论保证的框架。

### 3. Robust Local Fréchet Regression: Nonasymptotic Theory and Minimax Optimality

- 作者：Chang Jun Im, Daisuke Kurisu
- 日期：2026-10-16
- 方向：稳健统计、非参数回归、Minimax Theory
- 原文：https://arxiv.org/abs/2610.23757

**研究问题：** 当响应变量位于一般度量空间、并且具有重尾分布时，如何获得可靠的局部回归估计？

**核心贡献：** 构造基于距离的 Median-of-Means 稳健估计器，建立非渐近高概率误差界及对应的 minimax 下界。

**值得关注：** 在相应几何、光滑性和有限二阶矩条件下，获得经典非参数回归的逐点收敛速率，并且不要求响应空间满足特定熵条件。

**本期重点精读。**

### 4. Online Covariance Estimation for Averaged Stochastic Gradient Descent on Riemannian Manifolds

- 作者：Jingrong Guan, Liwei Jiang, Sen Na
- 日期：2026-10-24
- 方向：黎曼优化、在线统计推断
- 原文：https://arxiv.org/abs/2610.33419

**研究问题：** 流形上的随机梯度下降算法如何进行协方差估计与统计推断？

**核心贡献：** 利用向量传输整合不同切空间的信息，构造在线协方差估计器。

**值得关注：** 将欧氏随机优化的协方差推断思想扩展到黎曼流形。

### 5. Estimating the Wasserstein Barycenter of One-Dimensional Distributions Under Sparse Sampling

- 作者：James Peng, Florian Stijven, Linbo Wang, Peter B. Gilbert
- 日期：2026-06-08
- 方向：最优传输、分布值数据
- 原文：https://arxiv.org/abs/2606.10096

**研究问题：** 当每个个体的概率分布只能通过少量样本观测时，如何估计总体 Wasserstein barycenter？

**核心贡献：** 提出 marginal-constructed barycenter（MCB）估计器，避免直接完整估计每个个体分布，并研究其逐点一致性和渐近正态性。

**值得关注：** 将稀疏采样、分布数据和 Wasserstein 几何联系起来。

### 6. Inference for Fréchet Regression

- 作者：Wookyeong Song, Paromita Dubey, Hans-Georg Müller, Alexander Petersen
- 日期：2026-05-19
- 方向：Fréchet 回归、假设检验
- 原文：https://arxiv.org/abs/2605.19519

**研究问题：** 对度量空间值响应进行回归以后，如何检验解释变量是否具有显著作用？

**核心贡献：** 提出整体显著性检验及部分解释变量效应检验，结合随机乘子与 Cauchy combination 等工具进行统计推断。

**值得关注：** 将 Fréchet 回归从估计与预测扩展到正式的统计假设检验。

---

## C. 今日核心论文精读

### Robust Local Fréchet Regression: Nonasymptotic Theory and Minimax Optimality

**作者：** Chang Jun Im, Daisuke Kurisu（2026）

原文：https://arxiv.org/abs/2610.23757

### C.1 研究动机

考虑：

$$
X\in\mathbb R^d,\qquad Y\in(\mathcal Y,d_{\mathcal Y}),
$$

其中 $Y$ 属于一般度量空间。

Fréchet 回归的目标是条件 Fréchet mean：

$$
m_\oplus(x)
\in
\arg\min_{\omega\in\mathcal Y}
\mathbb E[d_{\mathcal Y}^2(Y,\omega)\mid X=x].
$$

经典局部 Fréchet 回归通过最小化核加权平方距离来估计这一函数。

然而，在重尾分布下，少数距离特别大的观测可能对局部估计产生较大影响。

**论文试图解决的问题：**

能否在仅有局部有限二阶矩、且不依赖响应空间熵条件的情况下，获得非欧局部回归的非渐近高概率误差保证？

### C.2 核心方法：Distance-Based Median-of-Means

算法思想可以分为三步。

**第一步：将观测分为 $B$ 个独立分组。**

设分组索引为：

$$
I_1,\ldots,I_B.
$$

**第二步：每个分组计算局部 Fréchet 估计器。**

$$
\hat m_{h,b}(x)
\in
\arg\min_{\omega\in\mathcal Y}
\sum_{i\in I_b}
K_h(X_i,x)d_{\mathcal Y}^2(Y_i,\omega).
$$

其中 $K_h$ 是带宽为 $h$ 的核权重。

**第三步：通过距离进行多数邻居聚合。**

对每个分组估计量，寻找包含严格超过一半分组估计结果的最小球半径：

$$
r_b=
\inf\left\{
r\geq0:
\sum_{\ell=1}^B
\mathbf 1\left[
d_{\mathcal Y}(\hat m_{h,b},\hat m_{h,\ell})\leq r
\right]>\frac B2
\right\}.
$$

最终选择：

$$
\hat b\in\arg\min_b r_b,
\qquad
\boxed{
\hat m^{RLFR}_{h,B}(x)=\hat m_{h,\hat b}(x).
}
$$

与普通均值聚合不同，这个方法只依赖于分组估计量之间的距离。

因此不需要响应变量能够进行向量加法，也不需要再对整个响应空间求一个新的 Fréchet median。

### C.3 主要理论结果

论文在两组不同的局部光滑性条件下建立误差界。

**情形一：局部二阶光滑性**

在二次增长、四点不等式、局部有限二阶矩、核函数和设计分布等相应条件下，论文证明了以下形式的高概率界：

$$
\boxed{
d_{\mathcal Y}
\left(
\hat m_{h,B}^{RLFR}(x_0),
m_\oplus(x_0)
\right)
\lesssim
h^2+
\sigma
\sqrt{
\frac{\log(1/\delta)}
{nh^df_X(x_0)}
}
}
$$

该界以至少 $1-\delta$ 的概率成立，并需要满足论文规定的带宽和局部样本量条件。

其中：

- $h^2$：局部平滑偏差；
- $nh^df_X(x_0)$：局部有效样本量的典型尺度；
- $\delta$：失败概率；
- $\sigma$：局部响应波动规模。

平衡偏差与随机误差后，对于固定置信水平，可得到逐点速率：

$$
\boxed{
n^{-2/(d+4)}.
}
$$

**情形二：条件分布满足 $W_1$-Hölder 连续性**

设：

$$
P_x=\mathcal L(Y\mid X=x),
$$

并假设：

$$
W_1(P_x,P_{x_0})
\leq
L_W\|x-x_0\|^\beta,
\qquad 0<\beta\leq1.
$$

相应误差界具有：

$$
h^\beta+
\sqrt{\frac{\log(1/\delta)}{nh^d}}
$$

的主要结构，其中省略了固定模型常数。

平衡两项后得到逐点速率：

$$
\boxed{
n^{-\beta/(2\beta+d)}.
}
$$

论文还通过嵌入高斯回归子类构造 minimax lower bounds，证明相应速率在所研究的模型类中达到最优阶。

### C.4 为什么能够获得高概率保证？

关键是将两个理论工具结合：

1. **Fréchet 估计稳定性**：借助二次增长条件、四点不等式和样本替换分析，控制单组加权 Fréchet 估计的误差。
2. **Median-of-Means 聚合**：即使每个分组估计只有一个严格大于 $1/2$ 的成功概率，也可以通过足够多的独立分组将整体失败概率显著降低。

由此构成：

$$
\text{有限二阶矩}
\longrightarrow
\text{分组误差控制}
\longrightarrow
\text{稳健聚合}
\longrightarrow
\text{高概率误差界}.
$$

### C.5 方法价值与局限

**主要创新：**

- 构造适用于一般度量空间的距离型稳健聚合方法；
- 处理随机核权重与异质条件响应分布；
- 建立有限样本高概率界；
- 得到与经典非参数回归匹配的逐点 minimax 速率；
- 避免依赖特定的响应空间熵条件。

**值得注意的限制：**

- 仍需要一定的几何结构，例如二次增长和四点不等式；
- 理论重点是逐点估计，不直接提供整个解释变量区域上的统一置信带；
- 分组数增加可以提升高概率稳健性，但也会减少每组样本量，带来效率损失。

### C.6 我的研究思考

以下属于进一步研究的可能方向，而不是论文已经证明的结论。

**方向一：统一收敛理论**

能否进一步控制：

$$
\sup_{x\in\mathcal X}
d_{\mathcal Y}
\left(
\hat m^{RLFR}(x),m_\oplus(x)
\right)?
$$

**方向二：自适应参数选择**

能否同时数据驱动地选择带宽 $h$ 和分组数 $B$，并保持理论误差保证？

**方向三：削弱几何假设**

如果二次增长条件不成立，或 Fréchet mean 不唯一，是否能够建立新的稳健估计理论？

### C.7 与今日基础知识的联系

今天学习了：

$$
\hat\mu_F
\in\arg\min_{\omega}
\sum_i d^2(Y_i,\omega).
$$

而论文研究的是它的条件化、局部化、加权化和稳健化：

$$
\text{Fréchet Mean}
\longrightarrow
\text{Local Fréchet Regression}
\longrightarrow
\text{Robust Local Fréchet Regression}.
$$

因此，这篇论文展示了一个基础数学概念如何发展成现代统计研究方法。

---

## D. 今日小结

### 三个核心收获

1. **数据空间决定统计方法。** 对于方向、矩阵、概率分布等对象，应尊重其天然几何结构。
2. **Fréchet mean 是推广统计平均的重要工具。** 它通过距离最小化定义，不需要向量加法。
3. **现代非欧统计研究不仅涉及估计方法，还研究稳健性、收敛速度和统计最优性。**

### 今日练习

设三个普通数值样本为：

$$
1,\quad3,\quad5.
$$

求其 Fréchet mean。

**答案：**

$$
\hat\mu_F
=
\arg\min_\mu
\left[
(1-\mu)^2+(3-\mu)^2+(5-\mu)^2
\right]
=3.
$$

### 下一期预告

**第 002 期：不同距离函数为什么会改变统计结果？**

基础课将介绍：

- 欧氏距离；
- 球面测地距离；
- 距离选择与 Fréchet mean 的关系。

前沿部分继续追踪新的非欧统计方法和理论研究。

---

## 参考文献

1. Im, C. J., & Kurisu, D. (2026). *Robust Local Fréchet Regression: Nonasymptotic Theory and Minimax Optimality*. [arXiv:2610.23757](https://arxiv.org/abs/2610.23757).
2. Im, C. J. (2026). *Adaptive Local Fréchet Regression with Polyspherical Predictors*. [arXiv:2611.13240](https://arxiv.org/abs/2611.13240).
3. Meng, K., & Perez, C. (2026). *Theoretical Foundations of Principal Manifold Estimation with Non-Euclidean Templates*. [arXiv:2604.04272](https://arxiv.org/abs/2604.04272).
4. Guan, J., Jiang, L., & Na, S. (2026). *Online Covariance Estimation for Averaged Stochastic Gradient Descent on Riemannian Manifolds*. [arXiv:2610.33419](https://arxiv.org/abs/2610.33419).
5. Peng, J., Stijven, F., Wang, L., & Gilbert, P. B. (2026). *Estimating the Wasserstein Barycenter of One-Dimensional Distributions Under Sparse Sampling*. [arXiv:2606.10096](https://arxiv.org/abs/2606.10096).
6. Song, W., Dubey, P., Müller, H.-G., & Petersen, A. (2026). *Inference for Fréchet Regression*. [arXiv:2605.19519](https://arxiv.org/abs/2605.19519).

---

*本期为中文学习与科研阅读笔记。数学结论和方法描述基于相关论文；对研究价值和潜在研究方向的评论为笔记整理者的理解，不代表原作者结论。*
