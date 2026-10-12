# IGA 非线性固体求解器 —— 第三阶段规格:动力学、Generalized-$\alpha$ 时间积分与实时可视化

> 本文档接续《实现规格说明》(单元层)与《第二阶段规格》(装配、牛顿迭代)。
>
> **前置假设(必须全部满足,否则不要进入本阶段):**
> - 单元级有限差分校验通过
> - 全局有限差分校验通过($\mathbf{K}_T$ 与 $\mathbf{R}$ 的一致性)
> - 静力测试 A(轴向自重杆)相对误差 $<10^{-5}$
> - 静力测试 C 观察到**二次收敛**
>
> 理由:动力学在静力残量外只多加了一个惯性项和一层时间递推。若静力层有 bug,动力学里会表现为"能量缓慢漂移""频率偏了几个百分点"这类症状,与时间积分格式本身的误差**完全无法区分**。
>
> **符号约定变更(重要):** 前两阶段用 $C$ 表示 Voigt 本构矩阵($6\times6$)。本阶段中 $\mathbf{C}$ 是时间积分格式里的**全局阻尼矩阵**,是完全不同的对象。**本文档起,Voigt 本构矩阵一律记作 $\mathbf{D}$**,代码里也建议同步改名(`Cmat` → `Dmat`)。这个改名成本极低,但不改的话第二部分的公式没法读。

------

## 第零部分:与静力求解器的关系

**这是本阶段最重要的一张表。** 动力学不是重写求解器,而是在静力求解器外面套一层时间循环 + 在残量里加一项。

| 组件 | 静力 | 动力学 | 复用程度 |
|---|---|---|---|
| 形函数 / Jacobian / $G_A$ / $\tilde w$ | L1 | L1 | **完全复用,一行不改** |
| $H,F,E,S$ | L2 | L2 | **完全复用** |
| $B_A$, $\tilde B_A$, $\mathbf{D}$ | L1/L2 | L1/L2 | **完全复用** |
| 单元刚度 $K_e$ | 每次迭代 | 每次迭代 | **完全复用** |
| 单元内力 $f_e^{\text{int}}$ | 每次迭代 | 每次迭代 | **完全复用** |
| 体力向量 $\mathbf{f}^{\text{ext}}$ | 一次 | 一次(死载荷) + 时间缩放 | 复用,加一个 $t$ 的标量因子 |
| **质量矩阵 $\mathbf{M}$** | 无 | **全程一次** | **新增(第一部分)** |
| 求解的矩阵 | $\mathbf{K}_T$ | $\mathbf{M}^*=\alpha_m\mathbf{M}+\alpha_f\beta\Delta t^2\mathbf{K}_T$ | 新增组合(第二部分) |
| 牛顿未知量 | $\Delta\mathbf{d}$ | $\Delta\mathbf{a}_{n+1}$ | 变更(第二部分) |
| 求值构形 | $\mathbf{d}$ | $\mathbf{d}_{n+\alpha_f}$(中间构形) | 变更 |
| 初始条件 | 无 | $\mathbf{d}_0,\mathbf{v}_0,\mathbf{a}_0$ | **新增(第三部分)** |

**关键认识:$\mathbf{d}_{n+\alpha_f}$ 只是一个普通的位移向量。** 把它喂给现有的 `assemble(d, K, f_int)`,拿回来的 $\mathbf{K}_T$ 和 $\mathbf{f}^{\text{int}}$ 就是在该构形上求值的结果 —— 装配函数**完全不需要知道自己身处动力学计算中**。保持这个接口不变。

**Total Lagrangian 红利在动力学下更大:** L1 层的量($G_A$、$\tilde w$、$B_A$)不但在牛顿迭代中不变,在**整个时间历程中都不变**,因为它们只依赖参考构形。$\mathbf{M}$ 同理 —— 这是下一节的全部内容。

------

## 第一部分:质量矩阵

### 1.1 公式

由半离散弱形式,一致质量矩阵的分块为

$$\mathbf{M}_{AB} = \left(\int_{\Omega_0}\rho_0\,R_A R_B\,\mathrm{d}\Omega_0\right)\mathbf{I}_3$$

单元层实现是**整个项目里最简单的一段代码**:

```cpp
// nen×nen 标量矩阵,每单元
Me_scalar.setZero();
for (q in gauss points) {
    // R, w_tilde 都是 L1 层已有的量,零额外计算
    Me_scalar.selfadjointView<Eigen::Lower>()
             .rankUpdate(R, rho0 * w_tilde);      // Me += rho0*w̃ * R Rᵀ
}
Me_scalar = Me_scalar.selfadjointView<Eigen::Lower>();   // ⚠️ 同样别忘了展开
```

$R$ 是 $n_{en}\times1$ 的形函数值向量。注意这里又是一次 syrk,**与 $K_e$ 一样有"忘记展开对称视图"的陷阱**。

### 1.2 【必做】$\mathbf{M}$ 全程只组装一次

$\rho_0$、$R_A$、$\tilde w$ 全部只依赖参考构形 ⟹ $\mathbf{M}$ 是**常量**。在时间循环开始前组装一次,此后:

- 牛顿迭代中不重算
- 时间步之间不重算
- 载荷分级中不重算

⚠️ 若在牛顿循环内看到 `assemble_mass()`,那是纯粹的浪费,而且会掩盖"$\mathbf{M}$ 其实依赖了 $\mathbf{d}$"这种错误(比如误用了当前构形的 $\det J$)。

### 1.3 【必做】不要做集中质量(lumping)

隐式格式下 $\mathbf{M}$ 无论如何都要进入线性求解,对角化**不带来任何求解速度收益**,却会:

- 降低精度(NURBS 的一致质量矩阵频谱精度远优于 $C^0$ 有限元,这是 IGA 的主要卖点之一,集中化把它扔掉了)
- 引入需要专门验证的新代码路径

(顺带一提:NURBS 基函数满足 $\sum_A R_A=1$ 且 $R_A\ge0$,所以行和集中法是**可用**的、不会产生负质量 —— 这一点优于高阶 Lagrange 单元。但本阶段没有理由用它。将来若做显式动力学再说。)

### 1.4 存储与 $\mathbf{M}^*$ 的组装策略

$\mathbf{M}$ 的稀疏模式与 $\mathbf{K}_T$ **完全相同**(同一套 IEN 散射)。有两条路:

**路线 A(验证阶段推荐):** $\mathbf{M}$ 存成独立的 $3N_{\text{cp}}\times3N_{\text{cp}}$ 稀疏矩阵,每次迭代做
```cpp
Mstar = alpha_m * M + alpha_f * beta * dt*dt * K_T;
```
清晰、易于单独测试 $\mathbf{M}$。稀疏加法有开销但不影响正确性。

**路线 B(功能跑通后再切):** 在单元循环内直接组装有效切线,一次散射搞定:
```cpp
Ke_eff = alpha_f*beta*dt*dt * Ke + alpha_m * Me_expanded;
```
其中 `Me_expanded` 是把 `Me_scalar` 加到每个 $3\times3$ 块对角线上(与几何刚度的处理**完全一样**,复用同一段代码)。省掉一次稀疏加法和一份矩阵的存储。

⚠️ **不要一上来就走路线 B。** 路线 A 下你可以单独断言 $\mathbf{M}$ 的性质(1.5 节),路线 B 里 $\mathbf{M}$ 从来不单独存在,查不了。

### 1.5 $\mathbf{M}$ 的自检(零成本,建议常驻断言)

| 检验 | 内容 | 抓什么 |
|---|---|---|
| **总质量** | $\mathbf{1}^T\mathbf{M}\mathbf{1}/3 = \sum_{A,B}\int\rho_0R_AR_B = \rho_0 V_0$ | $\tilde w$(含 $\det J_{\text{parent}}$)、装配 |
| **对称性** | $\|\mathbf{M}-\mathbf{M}^T\|/\|\mathbf{M}\|<10^{-14}$ | syrk 未展开 |
| **正定性** | `SimplicialLLT` 分解 $\mathbf{M}$ 应成功 | 一致质量矩阵必然 SPD;失败即装配错 |
| **行和** | 第 $A$ 块的行和 $=\int\rho_0R_A\,\mathrm{d}\Omega$,全体求和 $=\rho_0V_0$ | 同第一项,更细粒度 |

第一项特别有用:长方体 $L\times b\times h=10\times1\times1$、$\rho_0=1$ 时应精确得到 $10.0$(几何是仿射的,Gauss 积分精确)。**这一项过不了就不要往下走**,它和静力阶段的"$\sum_A\mathbf{f}^{\text{ext}}_A=\rho_0\mathbf{g}V_0$"是同一类锐利测试。

------

## 第二部分:Generalized-$\alpha$ 时间积分格式

### 2.1 待求解的方程

组装后的半离散运动方程:

$$\mathbf{M}\mathbf{a} + \mathbf{f}^{\text{int}}(\mathbf{d}) = \mathbf{f}^{\text{ext}}(t)$$

其中 $\mathbf{v}=\dot{\mathbf{d}}$、$\mathbf{a}=\ddot{\mathbf{d}}$。这是一个二阶非线性常微分方程组,用 Generalized-$\alpha$ 离散。

### 2.2 Newmark 更新式(以 $\mathbf{a}_{n+1}$ 为未知量)

$$\mathbf{d}_{n+1} = \mathbf{d}_n + \Delta t\,\mathbf{v}_n + (\Delta t)^2\left[\left(\tfrac12-\beta\right)\mathbf{a}_n + \beta\,\mathbf{a}_{n+1}\right]$$
$$\mathbf{v}_{n+1} = \mathbf{v}_n + \Delta t\left[(1-\gamma)\mathbf{a}_n + \gamma\,\mathbf{a}_{n+1}\right]$$

⚠️ **牛顿迭代的未知量是 $\mathbf{a}_{n+1}$,不是 $\mathbf{d}_{n+1}$。** 也可以改写成以 $\mathbf{d}_{n+1}$ 为未知量(两者等价,有效切线差一个 $\beta\Delta t^2$ 的缩放),但**必须选定一个并全程贯彻**。以 $\mathbf{a}$ 为未知量是 Hughes 书中 predictor-multicorrector 的写法,与文献对照方便,本文采用它。

### 2.3 中间时刻量

$$\mathbf{d}_{n+\alpha_f} = \mathbf{d}_n + \alpha_f(\mathbf{d}_{n+1}-\mathbf{d}_n)$$
$$\mathbf{v}_{n+\alpha_f} = \mathbf{v}_n + \alpha_f(\mathbf{v}_{n+1}-\mathbf{v}_n)$$
$$\mathbf{a}_{n+\alpha_m} = \mathbf{a}_n + \alpha_m(\mathbf{a}_{n+1}-\mathbf{a}_n)$$

⚠️ **位移/速度用 $\alpha_f$,加速度用 $\alpha_m$,两者不同。** 把它们写反是本节最典型的 bug,症状是数值阻尼行为不对(高频不衰减或过度衰减)、且 $\Delta t$ 收敛阶掉到一阶。

本问题中 $\mathbf{f}^{\text{int}}$ 不依赖速度,$\mathbf{v}_{n+\alpha_f}$ **实际不参与残量计算**,代码里可以不算。列在这里是因为 2.5 节的链式法则要用到它(将来加 Rayleigh 阻尼时就需要了)。

### 2.4 残量

$$\boxed{\mathbf{R}_{n+1} = \mathbf{M}\,\mathbf{a}_{n+\alpha_m} + \mathbf{f}^{\text{int}}\big(\mathbf{d}_{n+\alpha_f}\big) - \mathbf{f}^{\text{ext}}\big(t_{n+\alpha_f}\big)}$$

其中 $t_{n+\alpha_f} = t_n + \alpha_f\Delta t$。

**与静力残量 $\mathbf{R}=\mathbf{f}^{\text{int}}-\mathbf{f}^{\text{ext}}$ 的差别只有两处:** 多了 $\mathbf{M}\mathbf{a}_{n+\alpha_m}$;$\mathbf{f}^{\text{int}}$ 的求值点从 $\mathbf{d}$ 换成 $\mathbf{d}_{n+\alpha_f}$。

⚠️ **外力也要取在 $t_{n+\alpha_f}$,不是 $t_{n+1}$。** 常值重力下二者相同,所以这个错误在测试算例里**完全不会暴露**,直到某天加了时变载荷。一开始就写对。

### 2.5 有效切线

先给三个偏导的定义:

$$\mathbf{M}\equiv\frac{\partial\mathbf{R}}{\partial\mathbf{a}_{n+\alpha_m}},\qquad
\mathbf{C}\equiv\frac{\partial\mathbf{R}}{\partial\mathbf{v}_{n+\alpha_f}}=\mathbf{0},\qquad
\mathbf{K}\equiv\frac{\partial\mathbf{R}}{\partial\mathbf{d}_{n+\alpha_f}}=\mathbf{K}_T\big(\mathbf{d}_{n+\alpha_f}\big)$$

$\mathbf{C}=\mathbf{0}$ 是因为超弹性无粘性项。$\mathbf{K}$ **就是第二阶段装配出来的那个切线刚度矩阵**,只是求值构形换了。

由于未知量是 $\mathbf{a}_{n+1}$,把上述偏导通过 2.2、2.3 传递过去:

$$\frac{\partial\mathbf{a}_{n+\alpha_m}}{\partial\mathbf{a}_{n+1}}=\alpha_m,\qquad
\frac{\partial\mathbf{v}_{n+\alpha_f}}{\partial\mathbf{a}_{n+1}}=\alpha_f\gamma\Delta t,\qquad
\frac{\partial\mathbf{d}_{n+\alpha_f}}{\partial\mathbf{a}_{n+1}}=\alpha_f\beta(\Delta t)^2$$

链式法则给出

$$\mathbf{M}^*\equiv\frac{\partial\mathbf{R}_{n+1}}{\partial\mathbf{a}_{n+1}}
=\alpha_m\mathbf{M}+\alpha_f\gamma\Delta t\,\mathbf{C}+\alpha_f\beta(\Delta t)^2\mathbf{K}$$

本问题 $\mathbf{C}=\mathbf{0}$,故

$$\boxed{\mathbf{M}^* = \alpha_m\mathbf{M} + \alpha_f\beta(\Delta t)^2\,\mathbf{K}_T\big(\mathbf{d}_{n+\alpha_f}\big)}$$

**性质:** $\mathbf{M}$ 对称正定、$\mathbf{K}_T$ 对称 ⟹ $\mathbf{M}^*$ **对称**。第一阶段的 syrk 优化、`SimplicialLDLT` 全部继续有效。

⚠️ 更进一步:$\Delta t$ 较小时 $\alpha_m\mathbf{M}$ 主导,$\mathbf{M}^*$ 往往**比 $\mathbf{K}_T$ 条件数更好**,即动力学的线性系统比静力的更好解。这是好消息,但也意味着"求解器不报错"在动力学里比在静力里更没有诊断价值。

### 2.6 参数取法

二阶精度要求:
$$\gamma = \tfrac12 - \alpha_f + \alpha_m,\qquad \beta = \tfrac14(1-\alpha_f+\alpha_m)^2$$

无条件稳定要求:
$$\alpha_m \ge \alpha_f \ge \tfrac12$$

最优高频阻尼下,由谱半径 $\rho_\infty\in[0,1]$ 唯一确定:
$$\alpha_m = \frac{2-\rho_\infty}{1+\rho_\infty},\qquad \alpha_f = \frac{1}{1+\rho_\infty}$$

⚠️ **$\alpha$ 的约定在文献中有两套。** 本文(与 Hughes 的 IGA 书、Jansen–Collis–Whiting 一致)中 $\alpha=1$ 对应全隐式;另一套(Chung–Hulbert 原文)用 $1-\alpha$,取值范围完全不同。**照抄公式前先确认是哪一套**,否则会得到一个稳定但精度一阶的格式 —— 症状是"能跑、能收敛,但 $\Delta t$ 减半误差只降一半"。

**推荐取值表:**

| 用途 | $\rho_\infty$ | $\alpha_m$ | $\alpha_f$ | $\gamma$ | $\beta$ |
|---|---|---|---|---|---|
| 能量守恒测试(线性区) | $1.0$ | $1/2$ | $1/2$ | $1/2$ | $1/4$ |
| **常规计算(默认)** | $0.8$ | $2/3$ | $5/9$ | $11/18$ | $25/81$ |
| 强非线性 / 需要压制高频 | $0.5$ | $1$ | $2/3$ | $5/6$ | $25/36$ |
| 高频一步湮灭 | $0.0$ | $2$ | $1$ | $3/2$ | $9/4$ |

$\rho_\infty=1$ 退化为中点法则(无数值阻尼),**只用于能量测试**,非线性问题下鲁棒性不足。

⚠️ **这些参数必须写成由 $\rho_\infty$ 计算的函数,不要硬编码数值。** 上表的值应作为单元测试的期望值。

------

## 第三部分:初始条件与边界条件

### 3.1 【最容易错】初始加速度必须解出来,不能置零

给定 $\mathbf{d}_0$、$\mathbf{v}_0$ 后,$\mathbf{a}_0$ **不是自由的**,由 $t=0$ 时刻的运动方程确定:

$$\mathbf{M}\,\mathbf{a}_0 = \mathbf{f}^{\text{ext}}(0) - \mathbf{f}^{\text{int}}(\mathbf{d}_0)$$

这是一个**线性**方程组(系数矩阵是常量 $\mathbf{M}$),在自由自由度上求解,受约束自由度上 $\mathbf{a}_0=\mathbf{0}$。

⚠️ **直接令 $\mathbf{a}_0=\mathbf{0}$ 是错的**,除非初始时刻恰好静力平衡。症状:第一步残量异常大、前几步出现虚假高频振荡、能量曲线开头有一个台阶。这个错误很隐蔽,因为格式是稳定的,虚假振荡会被数值阻尼慢慢吃掉,看起来"最后还是对的"。

**建议:** $\mathbf{M}$ 的 `SimplicialLLT` 分解在 1.5 节的正定性自检中已经做过一次,直接复用它解 $\mathbf{a}_0$。

### 3.2 边界条件

齐次 Dirichlet 在动力学下的处理是静力的直接推广:

- 初始条件在受约束自由度上必须满足 BC:$\mathbf{d}_0=\mathbf{v}_0=\mathbf{a}_0=\mathbf{0}$
- 牛顿迭代中 $\Delta\mathbf{a}$ 在受约束自由度上恒为 $\mathbf{0}$
- 对 $\mathbf{M}^*$ **行、列都置零 + 对角置 1**,右端项对应位置置零(与静力完全相同的 `apply_bc`)

**由 2.2 的更新式,$\Delta\mathbf{a}=\mathbf{0}$ 自动保证 $\mathbf{d}_{n+1}$、$\mathbf{v}_{n+1}$ 在受约束自由度上保持为零。** 不需要额外处理。

⚠️ 装配时跳过受约束 triplet 的做法(第二阶段 3.3 节)在这里要**同时应用于 $\mathbf{M}$ 和 $\mathbf{K}_T$**。若只对 $\mathbf{K}_T$ 跳过而 $\mathbf{M}$ 全量组装,$\mathbf{M}^*$ 在受约束行上会残留质量项,再被 `apply_bc` 覆盖掉 —— 结果碰巧是对的,但如果哪天改成路线 B(1.4 节)就会出错。建议 $\mathbf{M}$ 也全量组装、统一由 `apply_bc` 处理,即路线 A 保持简单。

### 3.3 收敛判据(与静力不同)

静力用 $\|\mathbf{R}\|_{\text{free}}/\|\mathbf{f}^{\text{ext}}\|_{\text{free}}$。动力学下这个分母可能为零(自由振动、无外载),必须换:

$$\frac{\|\mathbf{R}^{(i)}\|_{\text{free}}}{\|\mathbf{R}^{(0)}\|_{\text{free}}} < \texttt{tol}_{\text{rel}}
\qquad\text{或}\qquad
\|\mathbf{R}^{(i)}\|_{\text{free}} < \texttt{tol}_{\text{abs}}$$

其中 $\mathbf{R}^{(0)}$ 是**本时间步**预测子处的残量。典型 $\texttt{tol}_{\text{rel}}=10^{-8}$。两个判据取"或",避免 $\mathbf{R}^{(0)}$ 本身已经很小时(接近平衡态)白白多迭代。

⚠️ **范数仍然只在自由自由度上取。** 受约束自由度上的残量现在是"支反力 + 惯性反力",量级更大,算进去同样导致永不收敛。

### 3.4 预测子

$$\mathbf{a}_{n+1}^{(0)} = \frac{\gamma-1}{\gamma}\,\mathbf{a}_n$$

代入 2.2 的速度更新式可验证这恰好给出 $\mathbf{v}_{n+1}^{(0)}=\mathbf{v}_n$(**同速度预测子**),$\mathbf{d}_{n+1}^{(0)}$ 随之由位移更新式确定。这是标准选择,通常比零阶预测子 $\mathbf{a}_{n+1}^{(0)}=\mathbf{a}_n$ 少一次迭代。

⚠️ 预测子只影响迭代次数,**不影响收敛后的解**。若换预测子后答案变了,说明牛顿没有真正收敛(tol 太松)。这是一个有用的自检。

------

## 第四部分:算法流程

### 4.1 完整伪代码

```
// ===== 全程一次(L0) =====
build space, span arrays, univariate basis tables
D = SVK Voigt matrix;  L_chol = cholesky(D)          // 注意改名 C→D
constrained_dofs = collect_constrained_dofs()
M = assemble_mass()                                   // 第一部分,只此一次
assert(mass self-checks pass)                         // 1.5
M_llt.compute(M_free)                                 // 供 3.1 用
f_ext_total = assemble_body_force()                   // 死载荷

// ===== 时间积分参数 =====
(alpha_m, alpha_f) = from_rho_inf(rho_inf)            // 2.6,不要硬编码
gamma = 0.5 - alpha_f + alpha_m
beta  = 0.25 * (1 - alpha_f + alpha_m)^2

// ===== 初始条件 =====
d = d0;  v = v0                                       // 满足 BC
assemble(d) -> _, f_int                               // 只要内力
a = M_llt.solve(f_ext(0) - f_int)  on free dofs; 0 on constrained

// ===== 时间循环 =====
for n in 0 .. N_steps-1:
    t_n     = n * dt
    t_alpha = t_n + alpha_f * dt

    // ---- 预测 ----
    a_new = (gamma - 1) / gamma * a                   // 3.4
    d_new = d + dt*v + dt^2*((0.5-beta)*a + beta*a_new)
    v_new = v + dt*((1-gamma)*a + gamma*a_new)

    for it in 1 .. max_iter:
        // ---- 中间时刻量 ----
        d_af = d + alpha_f*(d_new - d)
        a_am = a + alpha_m*(a_new - a)

        // ---- 残量(复用静力装配) ----
        assemble(d_af) -> K_T, f_int                  // ⚠️ 在中间构形上求值
        R = M*a_am + f_int - f_ext(t_alpha)

        res = norm(R restricted to free dofs)
        if it == 1: res0 = res
        print(n, it, res, res/res0)                   // ⚠️ 必须打印
        if res/res0 < tol_rel or res < tol_abs: break

        // ---- 有效切线 ----
        Mstar = alpha_m*M + alpha_f*beta*dt*dt*K_T    // 2.5
        apply_bc(Mstar, -R)
        solve Mstar * da = -R
        assert(solver.info() == Success)

        // ---- 更新(顺序不能乱) ----
        a_new += da
        d_new += alpha_f_free_note: beta*dt*dt * da   // 见 4.2
        v_new += gamma*dt * da
        assert(all det F > 0 at d_af)

    assert(converged)
    d = d_new;  v = v_new;  a = a_new
    record_state(t_n + dt, d)                         // 供可视化,第七部分
```

### 4.2 ⚠️ 迭代内更新 $\mathbf{d}$、$\mathbf{v}$ 的增量系数

这是本节最容易写错的三行。由 2.2 的更新式对 $\mathbf{a}_{n+1}$ 求导:

$$\Delta\mathbf{d}_{n+1} = \beta(\Delta t)^2\,\Delta\mathbf{a},\qquad
\Delta\mathbf{v}_{n+1} = \gamma\Delta t\,\Delta\mathbf{a}$$

**注意这里没有 $\alpha_f$。** $\alpha_f$ 只出现在从 $\mathbf{d}_{n+1}$ 到 $\mathbf{d}_{n+\alpha_f}$ 的那一步(以及因此出现在 $\mathbf{M}^*$ 里)。把 $\alpha_f$ 错误地放进增量更新,会得到一个**仍然收敛但收敛到错误解**的格式 —— 残量下降正常、二次收敛可能仍然观察得到,但时间精度阶数不对。

**唯一可靠的检查手段是 4.3 的有限差分校验。**

### 4.3 【最关键】有效切线的有限差分校验

这是本阶段的 4.7(单元 FD)/5.4(全局 FD)的对应物,**优先级最高**。

固定一个时间步的 $(\mathbf{d}_n,\mathbf{v}_n,\mathbf{a}_n)$ 和 $\Delta t$,把残量看成 $\mathbf{a}_{n+1}$ 的函数 $\mathbf{R}(\mathbf{a}_{n+1})$,做中心差分:

$$\mathbf{M}^*_{IJ}\approx\frac{R_I(\mathbf{a}_{n+1}+\varepsilon\mathbf{e}_J)-R_I(\mathbf{a}_{n+1}-\varepsilon\mathbf{e}_J)}{2\varepsilon}$$

与解析的 $\mathbf{M}^*$ 逐项比对,**收敛阶应看到二阶**。用 $2\times1\times1$ 的小网格。

**这是唯一能同时抓住以下所有错误的检查:**
- $\alpha_m$ / $\alpha_f$ 用反(2.3)
- $\beta\Delta t^2$ / $\gamma\Delta t$ 系数错(2.5)
- 增量更新里多了或少了 $\alpha_f$(4.2)
- $\alpha$ 约定用错了那一套(2.6)
- $\mathbf{M}$ 装配错

**不通过就不要往下跑任何物理算例。** 上述错误中有好几个都会给出"稳定、收敛、看起来合理"的结果。

⚠️ 校验时 $\Delta t$ 不要取得太小,否则 $\alpha_m\mathbf{M}$ 完全主导,$\mathbf{K}_T$ 那一项的错误会被淹没在舍入误差里。取 $\Delta t$ 使 $\alpha_m\|\mathbf{M}\|\sim\alpha_f\beta\Delta t^2\|\mathbf{K}_T\|$ 量级相当,即 $\Delta t\sim\sqrt{\|\mathbf{M}\|/\|\mathbf{K}\|}$,大致是最高阶模态的周期。**分别用大、小两个 $\Delta t$ 各跑一次。**

------

## 第五部分:时间步长与算例设计

### 5.1 $\Delta t$ 不受稳定性限制,但受精度限制

Generalized-$\alpha$ 无条件稳定 ⟹ 没有 CFL 条件。$\Delta t$ 完全由**你想解析的模态**决定:

$$\Delta t \lesssim \frac{T_{\text{目标}}}{20}$$

其中 $T_{\text{目标}}$ 是最高的、你关心其精度的模态周期。更高的模态会被数值阻尼(由 $\rho_\infty$ 控制)压制掉,这正是 Generalized-$\alpha$ 的设计意图。

⚠️ **不要用最高频模态定 $\Delta t$。** 空间离散的最高频模态($\sim c/h$)是纯粹的离散产物、没有物理意义,按它定步长会得到一个能跑一辈子也算不完的程序。

⚠️ 反过来也别过头:$\Delta t$ 大到把**你关心的**模态也压掉了,结果会平滑地收敛到一个错误答案,不报任何错。

### 5.2 阶跃载荷会激发所有模态

从静止状态突然施加重力,是一个理想阶跃载荷。后果:

- 理论上激发全部模态,含最高频的非物理模态 ⟹ **必须开数值阻尼**($\rho_\infty<1$)
- 线性系统下响应是"绕静力解的等幅振荡",峰值 $=2\times$ 静力解(**动力放大系数 2**)

第二点是一个极好的测试判据(见 6.4)。

若要避免阶跃,把载荷按 $\mathbf{f}^{\text{ext}}(t)=s(t)\,\mathbf{f}^{\text{ext}}_{\text{total}}$ 缓升,$s(t)$ 用光滑函数(如 $s=\tfrac12(1-\cos\pi t/T_r)$,$T_r\gtrsim3T_1$)。**但验证阶段建议就用阶跃**,因为它有解析判据。

### 5.3 与静力解的一致性

**动力学求解器必须能复现静力解。** 两条路:

1. 缓升载荷 + 强阻尼($\rho_\infty=0$)+ 大 $\Delta t$,跑到定常 ⟹ 应收敛到静力解
2. 阶跃载荷 + 长时间平均 ⟹ 平均值应等于静力解

第 1 条是最直接的回归测试,**建议作为动力学的第一个物理算例**(见 6.3)。

------

## 第六部分:验证测试(按发现 bug 的效率排序)

### 6.1 $\mathbf{M}$ 的自检

见 1.5。**先做这个**,零成本。

### 6.2 有效切线的 FD 校验

见 4.3。**不通过不要继续。**

### 6.3 测试 D:自由平动(最锐利)

**设置:** 无 BC(自由漂浮)、无外力、$\mathbf{d}_0=\mathbf{0}$、$\mathbf{v}_0=$ 均匀速度场(所有控制点同一个 $\mathbf{v}^*$)。

**精确解:** 刚体匀速平动。$\mathbf{d}(t)$ 处处等于 $\mathbf{v}^*t$,$\mathbf{a}\equiv\mathbf{0}$,$E\equiv\mathbf{0}$,$\mathbf{f}^{\text{int}}\equiv\mathbf{0}$。

**判据(极强):**
- 每个时间步的残量在**第一次迭代**就应满足收敛判据(预测子已经是精确解)
- $\max|E| < 10^{-14}$,全时程
- 总线动量 $\mathbf{p}=\mathbf{M}\mathbf{v}$ 的每个分量守恒到机器精度
- $\mathbf{d}(t)$ 与 $\mathbf{v}^*t$ 逐分量相差 $<10^{-12}$

**这个测试能抓:** $\mathbf{M}$ 装配、初始条件设置、Newmark 更新式、自由度排序。而且**离散误差与时间离散误差都严格为零**,任何偏差都是 bug。地位等同于静力阶段的测试 A。

⚠️ 注意此时 $\mathbf{K}_T$ 奇异(6 个刚体模态),但 $\mathbf{M}^*=\alpha_m\mathbf{M}+\dots$ **仍然正定**(质量矩阵没有零模)。这是动力学相对静力的一个便利:自由体问题是可解的。

### 6.4 测试 E:悬臂梁阶跃重力(动力放大系数)

**设置:** 静力测试 B 的完全相同配置($L=10$、$b=h=1$、$p=3$、$16\times2\times2$、$E=10^7$、$\nu=0$、$\rho_0g_z=1$),加上 $\rho_0=1$,$t=0$ 时刻突然施加重力,$\mathbf{v}_0=\mathbf{0}$。

**参考量(先用静力解和梁理论算出来):**

| 量 | 值 |
|---|---|
| 静力挠度 $w_{\text{static}}$(EB) | $1.5\times10^{-3}$ |
| 一阶弯曲频率 $\omega_1=(1.875104/L)^2\sqrt{EI/(\rho_0A)}$ | $32.10$ rad/s |
| 一阶周期 $T_1$ | $0.1957$ s |
| 建议 $\Delta t = T_1/40$ | $4.9\times10^{-3}$ s |
| 建议总时长 | $3T_1\approx0.6$ s |

($I=bh^3/12=1/12$,$A=1$,$EI=8.333\times10^5$,$\sqrt{EI/\rho_0A}=912.87$)

**判据:**
- 梁端挠度在 $0$ 与 $\approx2w_{\text{static}}=3.0\times10^{-3}$ 之间振荡
- 振荡周期与 $T_1$ 相差**几个百分点内**(三维实体解含剪切与转动惯量,比 EB 略低,$L/h=10$ 时约低 $2$–$3\%$)
- 时间平均值 $\approx w_{\text{static}}$
- $\rho_\infty=1$ 时振幅在若干周期内基本不衰减;$\rho_\infty=0.8$ 时缓慢衰减

⚠️ **梁端挠度的后处理必须在 $\boldsymbol\xi^*=(1,0.5,0.5)$ 处求值**(第二阶段 5.2 的陷阱),不能取 $\max_A|\mathbf{d}_A|$。

⚠️ 若观察到的周期偏差大于 $10\%$,先怀疑 $\rho_0$ 的量纲/取值,再怀疑 $\mathbf{M}$,最后才怀疑积分格式。

### 6.5 测试 F:纵向自由振动(频率判据最锐利)

**设置:** 静力测试 A 的杆($L=10$、$\nu=0$、$E=10^7$、$\rho_0=1$),$\xi=0$ 面固定,无重力,初始位移取静力测试 A 的解(或线性一阶模态形状),$\mathbf{v}_0=\mathbf{0}$。

**精确解:** 固定–自由杆的一维波动方程,$c=\sqrt{E/\rho_0}=3162.28$,

$$\omega_1 = \frac{\pi c}{2L} = 496.73\ \text{rad/s},\qquad T_1 = 1.2649\times10^{-2}\ \text{s}$$

**判据:** 数值周期与 $T_1$ 相差 $<1\%$($\nu=0$ 使三维解严格退化为一维,没有 6.4 那样的模型误差,所以判据可以卡得更紧)。

⚠️ 这个算例对 $\mathbf{M}$ 的错误极其敏感:$\mathbf{M}$ 整体差一个常数因子 $\kappa$,频率就差 $\sqrt\kappa$。若频率恰好差 $\sqrt{2}$、$\sqrt{3}$ 之类的因子,先去查 $\rho_0$ 和 $\tilde w$。

### 6.6 $\Delta t$ 收敛阶(验证二阶精度)

取 6.5 的算例,用 $\Delta t$、$\Delta t/2$、$\Delta t/4$ 各跑到同一个 $t_{\text{end}}$,以最细的解为参考,误差应按 $4$ 倍下降。

**这是检验 $\gamma$、$\beta$、$\alpha_m$、$\alpha_f$ 之间一致性的最终手段。** 若只看到一阶收敛,说明 2.6 的二阶精度条件没被满足 —— 最常见原因是 $\alpha$ 约定用错了那一套。

⚠️ 用 $\rho_\infty=1$ 做这个测试,数值阻尼会污染收敛阶的观测。

### 6.7 能量守恒

取 $\rho_\infty=1$(中点法则)、无外力、小初始扰动(线性区),监测

$$\mathcal{E}(t) = \underbrace{\tfrac12\mathbf{v}^T\mathbf{M}\mathbf{v}}_{\text{动能}} + \underbrace{\int_{\Omega_0}\tfrac12 S:E\,\mathrm{d}\Omega_0}_{\text{应变能(SVK)}}$$

应变能可以在单元循环里顺手累加,几乎零成本。

**判据(注意分两种情况):**
- **线性/小应变区:** 中点法则精确守恒能量,$\mathcal{E}$ 的相对漂移应 $<10^{-10}$
- **大变形非线性区:** 中点法则**不**严格守恒能量(需要 Simo–Tarnow 类能量–动量守恒格式才行),会看到 $O(\Delta t^2)$ 的缓慢漂移。**这是格式的已知性质,不是 bug。**

⚠️ 不知道第二条的话,会花很多时间去查一个不存在的 bug。

### 6.8 补充检验

| 检验 | 内容 | 抓什么 |
|---|---|---|
| **静力回归** | 缓升载荷 + $\rho_\infty=0$ + 大 $\Delta t$ ⟹ 收敛到静力解 | 整条管路的一致性 |
| **零初速零载荷** | $\mathbf{d}_0=\mathbf{v}_0=\mathbf{0}$、$\mathbf{f}^{\text{ext}}=\mathbf{0}$ ⟹ 全时程 $\mathbf{d}\equiv\mathbf{0}$ | $\mathbf{a}_0$ 求解、预测子 |
| **二次收敛** | 每个时间步内 $\mathbf{R}$ 的下降序列 | $\mathbf{M}^*$ 的正确性 |
| **$\mathbf{M}^*$ 对称性** | $\|\mathbf{M}^*-\mathbf{M}^{*T}\|/\|\mathbf{M}^*\|<10^{-12}$ | syrk 未展开、BC 只置零了行 |
| **角动量** | 自由体、无外力 ⟹ 总角动量守恒 | $\mathbf{M}$、内力自平衡性 |
| **$\rho_\infty$ 效果** | 给一个高频初始扰动,$\rho_\infty=0$ 时应几步内消失 | $\alpha_m,\alpha_f$ 是否真的起作用 |

------

## 第七部分:Qt6 + OpenGL 实时可视化

### 7.1 线程模型

⚠️ **求解绝不能跑在 GUI 线程上。** 一个时间步可能要几百毫秒,期间界面完全冻结,连"停止"按钮都点不了。

```
GUI 线程                          求解线程 (QThread)
--------                          ------------------
QOpenGLWidget::paintGL()          for n in steps:
  读 front buffer                     Newton solve
  更新 VBO                            写 back buffer
  绘制                                原子交换 buffer
                                      emit stepFinished(n, t)
```

**状态传递用三缓冲 + 原子指针交换**,不要用互斥锁保护整个位移向量 —— 渲染会被求解阻塞。$\mathbf{d}$ 只有 $3N_{\text{cp}}$ 个 double(测试算例约 $1425$ 个,即 $11$ KB),三份拷贝毫无压力。

Qt 侧用 `QObject::moveToThread` + 信号槽(`Qt::QueuedConnection`)传递"第 $n$ 步完成"的通知,**不要用信号槽传位移数据本身**(会拷贝且有队列延迟)。

### 7.2 【必做】可视化采样矩阵只算一次

这是本节最重要的优化,量级级别。

在参考构形上固定一组可视化采样点 $\{\boldsymbol\xi^*_k\}$(通常是表面上的一张结构化网格)。由于**参考构形永不改变**,每个采样点处的 $R_A(\boldsymbol\xi^*_k)$ 是常量。把它们组装成一个稀疏矩阵

$$\boldsymbol\Phi \in \mathbb{R}^{n_{\text{vis}}\times N_{\text{cp}}},\qquad \Phi_{kA}=R_A(\boldsymbol\xi^*_k)$$

每行只有 $d^3$ 个非零元。于是每帧的顶点位置只是**一次稀疏矩阵乘法**:

$$\mathbf{x}_{\text{vis}} = \boldsymbol\Phi\,(\mathbf{P}+\mathbf{d})$$

⚠️ **绝对不要每帧调 `FindSpan` + NURBS 求值。** 那是把 L1 层的工作搬到了渲染循环里,帧率会被彻底毁掉。$\boldsymbol\Phi$ 在程序启动时构造一次,这正是 `FindSpan` 的正当用途之一(与第二阶段 5.2 的后处理求值同源,可以复用同一段代码)。

法向同理:预存 $\boldsymbol\Phi_\xi$、$\boldsymbol\Phi_\eta$(采样点处的参数导数),每帧算
$$\mathbf{a}_1=\boldsymbol\Phi_\xi(\mathbf{P}+\mathbf{d}),\quad \mathbf{a}_2=\boldsymbol\Phi_\eta(\mathbf{P}+\mathbf{d}),\quad \mathbf{n}=\mathbf{a}_1\times\mathbf{a}_2$$
逐顶点归一化。**变形后法向必须重算**,否则光照在大变形下明显错误。

### 7.3 时间轴解耦

⚠️ **仿真 $\Delta t$ 与渲染帧间隔是两回事。** 直接"算一步画一帧"会导致:$\Delta t=5$ ms 而求解耗时 $200$ ms ⟹ 动画比真实慢 40 倍;反之 $\Delta t$ 很小时又快得看不清。

**做法:** 求解线程把每步的 $(t_n,\mathbf{d}_n)$ 存进一个历史缓冲;渲染线程按**墙钟时间** $\times$ 播放倍率去查这个缓冲,在相邻两步之间**线性插值**:

$$\mathbf{d}(t) = (1-s)\,\mathbf{d}_n + s\,\mathbf{d}_{n+1},\qquad s=\frac{t-t_n}{\Delta t}$$

位移的线性插值在视觉上完全够用(不需要用 Newmark 公式做三次插值)。

**存储代价可以忽略:** $1425$ 个 double $\times$ $8$ 字节 $\times$ $1000$ 步 $\approx11$ MB。**建议全程保存**,这样可以做进度条拖拽、回放、逐帧检查 —— 对调试非线性动力学极其有用。

### 7.4 ⚠️ 不要渲染中间构形

$\mathbf{d}_{n+\alpha_f}$ 和牛顿迭代中的 $\mathbf{d}^{(i)}_{n+1}$ **都不是物理状态**。只在时间步收敛后推送 $\mathbf{d}_{n+1}$。

调试时若确实想看迭代过程,单独做一个"逐迭代模式"并在界面上明确标注,不要混进正常播放。

### 7.5 ⚠️ 位移放大系数

静力测试 A 的位移量级是 $10^{-6}$,几何尺寸是 $10$ —— 直接画出来**完全看不出变形**。必须提供放大系数:

$$\mathbf{x}_{\text{vis}} = \boldsymbol\Phi\,(\mathbf{P}+s_{\text{amp}}\mathbf{d})$$

界面上给一个滑块,并**始终显示当前放大倍数**,否则很容易把放大 1000 倍的图误当作真实变形去判断结果对错。

### 7.6 OpenGL 实现要点

- `QOpenGLWidget` + `QOpenGLFunctions_3_3_Core`,在 `initializeGL()` 里建 VAO/VBO/着色器
- 顶点缓冲用 `QOpenGLBuffer::DynamicDraw`,每帧 `write()`(内部即 `glBufferSubData`)。索引缓冲(三角形拓扑)**只上传一次** —— 拓扑不随变形改变
- 位置和法向分开放两个 VBO,或交错放一个;都可以,但**索引缓冲一定要独立且静态**
- 用 `update()` / `requestUpdate()` 驱动重绘,配 `QTimer` 或 `QElapsedTimer` 控制播放速度。不要在 `paintGL()` 里 `sleep`
- 深度测试开、背面剔除**先别开**(变形可能翻转朝向,剔除掉会让人误以为几何破了)
- 颜色映射建议用 von Mises 应力或位移模长。**应力应在求解线程算好**(每个采样点一个标量,随 $\mathbf{d}$ 一起推送),不要在渲染线程重算 $F,E,S$

### 7.7 建议的诊断界面

动力学调试非常依赖时间序列曲线,QtCharts 画三条就够:

1. **监测点位移 vs 时间** —— 直接对照 6.4、6.5 的周期判据
2. **总能量 vs 时间**(动能/应变能分开画)—— 6.7 的漂移检查
3. **每步牛顿迭代次数 + 残量下降序列** —— 二次收敛的可视化

第 3 条尤其值得做:静力阶段是靠打印看收敛阶,动力学有几百上千个时间步,打印会淹没在滚屏里。**"每步迭代次数"这一条曲线突然抬高,通常就是物理上出问题的地方**(接近失稳、单元将要翻转)。

------

## 第八部分:bug 症状对照表

| 症状 | 最可能的原因 |
|---|---|
| 第一步残量巨大,之后正常 | $\mathbf{a}_0$ 被置零而非解出来(3.1) |
| 频率整体差一个 $\sqrt{2}$、$\sqrt{3}$ 之类的因子 | $\rho_0$ 取值/量纲错,或 $\mathbf{M}$ 漏了 $\det J_{\text{parent}}$ |
| 频率偏高、且加密网格不改善 | $\mathbf{M}$ 做了集中化,或 $\mathbf{M}$ 只填了对角块 |
| 振幅无缘无故衰减 | $\rho_\infty<1$ 的数值阻尼(**正常现象**),不是 bug |
| $\rho_\infty=1$ 下能量仍缓慢漂移 | 大变形区中点法则本就不守恒能量(6.7),先确认是否在线性区 |
| $\Delta t$ 减半误差只降一半(一阶) | $\alpha$ 约定用错那一套(2.6),或 $\gamma,\beta$ 与 $\alpha_m,\alpha_f$ 不匹配 |
| 稳定收敛但答案与静力回归不符 | 增量更新里的 $\beta\Delta t^2$/$\gamma\Delta t$ 系数错(4.2)——**跑 4.3** |
| 数值阻尼行为反了(改 $\rho_\infty$ 没效果) | $\alpha_m$ 与 $\alpha_f$ 用反了(2.3) |
| 每步都要迭代很多次,残量线性下降 | $\mathbf{M}^*$ 有错(继承自静力的 $\mathbf{K}_T$ 问题,或 $\alpha$ 系数错) |
| 前几步剧烈高频振荡后平静 | 阶跃载荷 + $\rho_\infty=1$(5.2);或 $\mathbf{a}_0$ 错 |
| 求解器报 `NumericalIssue` | 自由体问题下正常吗?$\mathbf{M}^*$ 应正定 —— 若报错,查 $\mathbf{M}$ 装配 |
| 时变载荷下结果有 $O(\Delta t)$ 偏差 | $\mathbf{f}^{\text{ext}}$ 取在了 $t_{n+1}$ 而非 $t_{n+\alpha_f}$(2.4) |
| 界面卡死 | 求解跑在 GUI 线程(7.1) |
| 帧率极低但求解很快 | 每帧重新做 NURBS 求值(7.2) |
| 动画速度与 $\Delta t$ 绑定 | 时间轴未解耦(7.3) |
| 画面看不出变形 | 缺放大系数(7.5) |

------

## 附录:实施顺序

1. **符号改名** $C\to\mathbf{D}$(Voigt 本构矩阵),腾出 $\mathbf{C}$ —— 5 分钟,先做
2. **质量矩阵组装**(第一部分)—— 跑 1.5 的全部自检,**总质量必须精确**
3. **时间积分参数函数** `from_rho_inf()` —— 用 2.6 的表做单元测试
4. **初始加速度求解**(3.1)
5. **残量 + 有效切线**(2.4、2.5)—— **跑 4.3 的 FD 校验,大小两个 $\Delta t$ 各一次**
6. **时间循环**(4.1)—— 跑 6.8 的"零初速零载荷"
7. **测试 D:自由平动**(6.3)—— 判据 $10^{-12}$,最锐利
8. **测试 F:纵向自由振动**(6.5)—— 频率 $<1\%$
9. **$\Delta t$ 收敛阶**(6.6)—— **确认二阶**
10. **测试 E:悬臂梁阶跃重力**(6.4)—— 动力放大系数 $\approx2$
11. **能量守恒**(6.7)—— 先在线性区
12. **静力回归**(6.8)
13. **可视化:采样矩阵 $\boldsymbol\Phi$ + 静态渲染**(7.2)—— 先渲染静力解,确认几何正确
14. **可视化:线程模型 + 时间轴解耦**(7.1、7.3)
15. **诊断界面**(7.7)
16. 后续:Rayleigh 阻尼(公式已在 2.5 中备好)、面力/随动载荷、路线 B 的 $\mathbf{M}^*$ 直接组装、CSR 索引优化

## 附录:参考文献

- Cottrell, Hughes & Bazilevs, *Isogeometric Analysis: Toward Integration of CAD and FEA*, ch. 7 —— Generalized-$\alpha$ 格式、NURBS 离散的频谱性质、质量矩阵
- Chung & Hulbert (1993), *A time integration algorithm for structural dynamics with improved numerical dissipation: the generalized-$\alpha$ method* —— 原始格式(注意 $\alpha$ 约定与本文相反)
- Jansen, Whiting & Hulbert (2000), *A generalized-$\alpha$ method for integrating the filtered Navier–Stokes equations* —— 本文采用的 $\alpha$ 约定、predictor-multicorrector 流程
- Hughes, *The Finite Element Method: Linear Static and Dynamic Finite Element Analysis*, ch. 9 —— Newmark 族、稳定性与精度分析
- Belytschko et al., *Nonlinear Finite Elements for Continua and Structures*, ch. 6 —— 非线性动力学的隐式求解、能量守恒讨论
- Simo & Tarnow (1992), *The discrete energy-momentum method* —— 若将来需要严格能量守恒格式
