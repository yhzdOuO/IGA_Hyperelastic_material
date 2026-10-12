# IGA 非线性固体求解器 —— 实现规格说明

> 目标:三维超弹性 Total Lagrangian 非线性 IGA 求解器,单 patch。
> 当前阶段的首要目的是**验证公式推导的正确性**,而非构造高性能生产级求解器。
> 因此:**代码清晰、与公式逐项对应、可断言验证** 优先于极致优化。
> 但下文标注为「必做」的优化都是量级级别或几乎零成本的,应当直接采用。

---

## 第一部分:已经完成的部分

### 1.1 已有的数据结构

**`KnotVector`**(`include/nurbs/KnotVector.hpp`)

封装结点向量 $U = \{U_0, \dots, U_m\}$,内部为 `std::vector<double>`。

已提供接口:
- `operator[](int i)` —— 访问 $U_i$
- `size()` —— 结点个数,即 $m+1$
- `lastIndex()` —— $m$
- `inferLastCtrlPtIndex(int p)` —— 返回 $n = m - p - 1$(该方向控制点个数为 $n+1$)
- `frontValue()` / `backValue()` / `raw()`

构造时校验:至少 2 个结点、序列非递减。

**结论:这个类不需要改动。** 它是一维标量序列,`std::vector<double>` 是恰当的选择,无需迁移到 Eigen。

**`ControlPoint<Dim>` / `ControlPointArray<Dim>`**(`include/nurbs/ControlPoints.hpp`)

当前实现:`ControlPoint<Dim>` 含 `std::array<double,Dim> coord` 与 `double weight`;`ControlPointArray<Dim>` 是其 `std::vector` 封装。

**结论:这个需要改动,见 2.1 节。**

另有自由函数 `isKnotControlConsistent(numKnots, p, numCtrlPts)`,校验 $m = n + p + 1$。保留。

### 1.2 已有的算法

三个算法已实现于 `include/nurbs/NurbsAlgorithms.hpp` + `src/NurbsAlgorithms.cpp`,已通过独立的递归 Cox–de Boor 参考实现交叉验证(含 $p=2$ 非均匀重结点、$p=3$ 两组场景,partition of unity、导数和为零、超阶导数为零等检查全部通过)。

**`int FindSpan(int n, int p, double u, const KnotVector& U)`**
The NURBS Book 算法 A2.1。返回跨度下标 $i$ 使 $u \in [U_i, U_{i+1})$。已加入定义域检查(越界抛 `std::domain_error`)。

**`void BasisFuns(int i, double u, int p, const KnotVector& U, std::vector<double>& N)`**
算法 A2.2。返回 $p+1$ 个非零 B 样条基函数值,`N[j]` $= N_{i-p+j,p}(u)$。

**`void DersBasisFuns(int i, double u, int p, int nDer, const KnotVector& U, std::vector<std::vector<double>>& ders)`**
算法 A2.3。返回 `ders[k][j]` $= \mathrm{d}^k N_{i-p+j,p}/\mathrm{d}u^k$。

注意:书中 A2.1 与 A2.3 都用符号 $n$ 但语义完全不同(前者是最后一个控制点下标,后者是最高求导阶数),本实现已将后者重命名为 `nDer`。`nDer > p` 时高阶导数按定义返回 0(已 clamp,不会越界)。

### 1.3 当前能力边界

**目前只能做到:** 给定单一参数方向、给定跨度下标,计算该跨度内非零的 **B 样条**基函数及其对**参数域**的导数。

**还不能做到:** 有理化(NURBS)、张量积、物理域导数、以及任何力学量。

---

## 第二部分:需要完成的部分

按实现顺序排列。**每一层验证通过后再进入下一层**,不要全部写完再一起调试。

### 2.1 【必做】数据结构迁移到 Eigen

采用 Eigen。理由不是"矩阵运算方便",而是固定尺寸类型的实际收益:`Matrix3d` 是栈上 9 个 double,`.inverse()` / `.determinant()` 对 $3\times3$ 有闭式特化并完全展开,无堆分配、无循环。

**类型选择表:**

| 对象 | 类型 |
|---|---|
| $J$, $J^{-1}$, $F$, $H$, $E$, $S$, $C_{\text{右柯西}}$ | `Eigen::Matrix3d` |
| $G_A = \nabla_X R_A$ | `Eigen::Vector3d` |
| $B_A$, $\tilde B_A$ | `Eigen::Matrix<double,6,3>` |
| 本构矩阵 $C$(Voigt $6\times6$) | `Eigen::Matrix<double,6,6>` |
| $\nabla_\xi N$, $\nabla_X N$, $P$(控制点坐标), $u$(控制点位移) | `Eigen::Matrix<double,Eigen::Dynamic,3>`,行数 $=n_{en}$ |
| $\tilde B$(拼接) | `Eigen::Matrix<double,6,Eigen::Dynamic>`,列数 $=3n_{en}$ |
| $K^e$, $f^e$ | `Eigen::MatrixXd` / `Eigen::VectorXd` |
| 全局 $K$ | `Eigen::SparseMatrix<double>` + `SimplicialLDLT` 或 `SparseLU` |

**要改的两处:**

1. `ControlPointArray<Dim>` → 改为 `Eigen::Matrix<double,Dynamic,3>`(记作 `P`,第 $A$ 行是控制点 $A$ 的坐标)+ 独立的 `Eigen::VectorXd weights`。这样 $J = P^T \nabla_\xi N$ 就是一个表达式而不是手写循环。
2. `DersBasisFuns` 的输出类型 `std::vector<std::vector<double>>` **必须换掉**。嵌套 vector 每行一次堆分配、内存散布、cache 极差,而这个函数在预计算阶段被反复调用。改为扁平数组(`std::vector<double>` + 手工 stride)或 `Eigen::Matrix<double,Dynamic,Dynamic>`。

**Eigen 使用禁忌(务必遵守):**

- **绝对不要对 Eigen 表达式用 `auto`。** `auto x = A * B` 存的是惰性表达式模板,操作数一旦失效即悬垂引用。这是 Eigen 头号陷阱。显式写出类型或用 `.eval()`。
- 别名:`A = A * B` Eigen 会自动开临时量;确定无别名时用 `.noalias()` 跳过。
- Debug 构建慢十倍以上,**验证阶段也要开 `-O2`**,否则会误判性能。
- `Matrix3d` 无对齐问题;但若类成员含 `Vector4d`/`Matrix4d`,需 `EIGEN_MAKE_ALIGNED_OPERATOR_NEW`。

### 2.2 【必做】张量积空间与 IEN

**关键概念区分(不要混淆):**
- **参数维度 $d_p$** —— $\xi,\eta,\zeta$ 有几个方向,决定基函数如何做张量积。本项目 $d_p = 3$。
- **物理嵌入维度** —— 控制点是 2D 还是 3D 坐标。本项目为 3。

这两者独立(NURBS 曲面 $d_p=2$ 但控制点可在 3D 空间中)。IEN 只关心 $d_p$。

**建议新建类 `TensorProductNurbsSpace`**,持有:
- 3 个 `KnotVector`
- 3 个 degree $p_1,p_2,p_3$
- 控制点矩阵 `P` 与权重 `weights`

**单元的定义:** 每个参数方向上,**非零长度**的结点跨度(即满足 $U_i \neq U_{i+1}$ 且 $p \le i \le n$ 的下标 $i$),三个方向做张量积。

**必须预先建立的:** 每个方向一个数组
```
span[d][e_d] = 该方向第 e_d 个非空跨度对应的 FindSpan 意义下的跨度下标
```
总共只有约 $3n$ 个整数。构造方式:扫一遍结点向量,把满足上述条件的 $i$ 依次记录。

**编号约定(必须全局统一,否则一定出 bug):**

局部基函数编号,$\xi$ 方向最快变化:
$$a = a_1 + a_2 d_1 + a_3 d_1 d_2, \qquad a_\alpha \in [0, p_\alpha],\ d_\alpha = p_\alpha + 1$$

全局基函数编号(0-based),同样 $\xi$ 最快:
$$A = i_1 + i_2 n_1 + i_3 n_1 n_2, \qquad i_\alpha = s_\alpha - p_\alpha + a_\alpha$$
其中 $n_\alpha$ 是该方向控制点个数。

**IEN 不要物化成数组。** 理由:$\mathrm{IEN}$ 有 $n^3 d^3$ 个条目,$n=100, p=3$ 时约 256 MB 内存,而它可以 O(1) 现算。实现为**惰性函数/视图**:

```
给定单元 e → 解出 (e1,e2,e3) → 查 span[] 得 (s1,s2,s3)
A_base = (s1-p1) + (s2-p2)*n1 + (s3-p3)*n1*n2
然后三重循环 a3,a2,a1:  A = A_base + a1 + a2*n1 + a3*n1*n2
```

用**强度削减**:$a_1$ 递增则 $A$ 加 1,$a_2$ 递增则加 $n_1$,$a_3$ 递增则加 $n_1n_2$ —— 内层全是加法,零乘法。局部编号 $a$ 就是顺序写入位置,不需计算。

保留一个 `IEN(a, e)` 的**函数**接口供调试与对照书上表格。

**装配循环内不要调 `FindSpan`。** 跨度下标已由 `span[]` 给出,二分搜索是重复劳动。`FindSpan` 的用途在任意点求值(后处理、collocation、点反演、点载荷)。

### 2.3 【必做】单变量基函数预计算表

$\xi$ 方向的单变量基函数值只依赖(该方向单元号 $e_1$, 该方向 Gauss 点序号 $q_1$),**与另两个方向完全无关**。因此 $(e_1,\cdot,\cdot)$ 这 $n^2$ 个单元共享同一组值。

每个方向建表:
```
table[d] : [该方向单元数] × [d_α 个 Gauss 点] × [2 × d_α 个 double]   (值 + 一阶导数)
```
总内存约 $6nd^2$ 个 double(几十 KB)。之后 `DersBasisFuns` 一次都不用再调。

**注意:即使不建全 patch 表,也必须做到"单元内按方向缓存"** —— 即每个单元只做 $3d$ 次单变量求值,而**不是**每个积分点做 $3d^3$ 次。这是白捡的,只需把求值放在单元层而非积分点层。

### 2.4 NURBS 形函数与参数域梯度

见第三部分 3.1。这是当前的核心待实现项。

### 2.5 两层 Jacobian 与物理域梯度

见 3.2。**注意有两个 Jacobian**,只有一个单元的验证算例里第二个恰好可能为 1 而被漏掉却不报错。

### 2.6 力学量:$H, F, E, S$

见 3.3。**注意 $E$ 不要用 $\frac12(F^TF - I)$**,有严重的浮点抵消。

### 2.7 $B_A$、$\tilde B_A$ 与切线刚度

见 3.4、3.5。

### 2.8 装配与牛顿迭代

本次不实现,但代码结构应为之留出接口。

---

## 第三部分:算法流程

### 计算分层(核心组织原则)

设 $d = p+1$(各方向非零基函数个数),$n_{en} = d^3$。

| 层级 | 内容 | 执行次数 | 牛顿迭代中是否重算 |
|---|---|---|---|
| **L0** | 本构矩阵 $C$(仅 SVK)、母元 Gauss 点与权重、$C$ 的 Cholesky 因子、单变量基函数表 | 全程 1 次 | 否 |
| **L1** | $R_A$、$\nabla_\xi R_A$、$J$、$J^{-1}$、$G_A = \nabla_X R_A$、$B_A$、积分权重 $\tilde w$ | 每积分点 1 次 | **否(TL 红利)** |
| **L2** | $H, F, E, S$;逐 $A$:$\tilde B_A$ | 每积分点每次迭代 | 是 |
| **L3** | $K^e$ 的双重循环 / 大矩阵乘 | 每积分点每次迭代 | 是 |

**L1 全程不变是 Total Lagrangian 格式的直接红利**:所有积分都在参考构形 $\Omega_0$ 上、所有导数都对参考坐标 $X$ 取,而参考构形永不改变。(Updated Lagrangian 下 $\nabla_x N$ 每次迭代都变,存不了。)

**内存决策:** 不要缓存 $\nabla_X N$ 到全部积分点($n=100,p=3$ 时约 98 GB)。值得缓存的是 $\tilde w = \det J \cdot \det J_{\text{parent}} \cdot w_{\text{gauss}}$(每积分点 1 个标量)和单变量表。

---

### 3.1 NURBS 形函数值与参数域梯度

**⚠️ 符号冲突警告:** 本节中权重和 $\sum W$ 与力学中的第二类 PK 应力 $S$ 是不同的量。本文档中前者统一记作 `Wsum`,后者记作 $S$。实现时务必用不同变量名。

**输入:** 单元 $e$、积分点三重下标 $(q_1,q_2,q_3)$
**输出:** $R_{abc}$ 及 $\partial R_{abc}/\partial\xi$、$\partial R_{abc}/\partial\eta$、$\partial R_{abc}/\partial\zeta$,各 $d^3$ 个

**数学依据:** NURBS 基函数
$$R_{abc} = \frac{\omega_{abc} N_a M_b L_c}{\sum_{a'b'c'} \omega_{a'b'c'} N_{a'} M_{b'} L_{c'}}$$
导数由商法则给出。

**Step 1** —— 由单元 $e$ 解出 $(e_1,e_2,e_3)$,查 `span[]` 得 $(s_1,s_2,s_3)$。

**Step 2** —— 查单变量表(2.3 节)取得该积分点上的:
- $\xi$ 向:$N_a, N'_a$,$a \in [0,p_1]$
- $\eta$ 向:$M_b, M'_b$,$b \in [0,p_2]$
- $\zeta$ 向:$L_c, L'_c$,$c \in [0,p_3]$

**Step 3【优化:因子共享】** —— 先建三张 $d^2$ 小表:
```
Pbc[b][c] = M_b  * L_c
Qbc[b][c] = M'_b * L_c
Tbc[b][c] = M_b  * L'_c
```

**Step 4** —— 对每个 $(a,b,c)$,按 2.2 节的局部编号约定算出扁平下标 $a_{\text{loc}}$,然后:
```
z  = ω[a][b][c] * N_a
z' = ω[a][b][c] * N'_a

W  [a_loc] = z  * Pbc[b][c]     // 加权基函数
Wxi[a_loc] = z' * Pbc[b][c]     // ∂W/∂ξ
Wet[a_loc] = z  * Qbc[b][c]     // ∂W/∂η
Wze[a_loc] = z  * Tbc[b][c]     // ∂W/∂ζ
```
每项 6 次乘法(朴素的三因子外积 + 乘权重需 12 次)。

**Step 5** —— 求四个标量和:
```
Wsum    = Σ W
Wsum_xi = Σ Wxi
Wsum_et = Σ Wet
Wsum_ze = Σ Wze
```

**Step 6【优化:消除除法】** —— 算**一次**倒数,并复用 $R$:
```
inv = 1.0 / Wsum

R[A]      = inv * W[A]
dR_dxi[A] = inv * (Wxi[A] - R[A] * Wsum_xi)
dR_det[A] = inv * (Wet[A] - R[A] * Wsum_et)
dR_dze[A] = inv * (Wze[A] - R[A] * Wsum_ze)
```

这与商法则原式 $\dfrac{\text{Wsum}\cdot W_\xi - W\cdot \text{Wsum}_\xi}{\text{Wsum}^2}$ 数学等价:提出 $1/\text{Wsum}$ 后用已算好的 $R$ 替换 $W/\text{Wsum}$。

**优化收益说明:** 朴素实现约 $21d^3$ 乘 $+\ 4d^3$ 除;Step 3 + Step 6 后约 $13d^3 + 3d^2$ 乘 $+$ **1 次除**。$d=4$ 时除法从 256 次降到 1 次(除法延迟是乘法的 10–20 倍)。

**数值稳定性说明:** $\text{Wsum} \in [\min\omega, \max\omega]$(因非零 $N$ 之和恒为 1、权重恒正),不会接近 0,无需额外保护。

**注:** 本例程在总时间中约占 1%(单元刚度形成是 $O(d^9)$,本例程是 $O(d^6)$)。上述优化几乎零成本故采纳,但不必为此牺牲代码可读性。

---

### 3.2 两层 Jacobian 与物理域梯度

记 $\nabla_\xi N$ 为 $n_{en}\times3$ 矩阵,第 $A$ 行为 $[\partial R_A/\partial\xi,\ \partial R_A/\partial\eta,\ \partial R_A/\partial\zeta]$。

**几何 Jacobian:**
$$J = P^T \nabla_\xi N, \qquad J_{ij} = \frac{\partial X_i}{\partial \xi_j}$$
($P$ 为 $n_{en}\times3$ 控制点坐标矩阵。Eigen:`Matrix3d J = P.transpose() * dN_dxi;`)

**必须断言 $\det J > 0$** 于所有积分点。参数化退化时这是第一道防线。

**物理域梯度(注意方向,不要转置):**
$$\nabla_X N = \nabla_\xi N \cdot J^{-1}$$
即 $n_{en}\times3$ 右乘 $3\times3$。分量形式 $\dfrac{\partial R_A}{\partial X_i} = \sum_k (J^{-1})_{ki}\dfrac{\partial R_A}{\partial \xi_k}$。

记 $G_A = \nabla_X R_A$ 为 $\nabla_X N$ 的第 $A$ 行(转置为 `Vector3d`)。

**母元 Jacobian(容易漏!):** Gauss 点定义在母元 $[-1,1]^3$ 上,须先仿射映到参数域:
$$J_{\text{parent}} = \mathrm{diag}\!\left(\frac{\xi_{s_1+1}-\xi_{s_1}}{2},\ \frac{\eta_{s_2+1}-\eta_{s_2}}{2},\ \frac{\zeta_{s_3+1}-\zeta_{s_3}}{2}\right)$$
是对角阵,$\det J_{\text{parent}}$ 为三项之积。

**积分权重:**
$$\tilde w = \det J \cdot \det J_{\text{parent}} \cdot w_{\text{gauss}}$$

⚠️ 只有一个单元的验证算例中 $\det J_{\text{parent}}$ 可能恰好为 1,漏掉不会报错,加密网格后才暴露。**一开始就显式写出来。**

---

### 3.3 位移梯度、变形梯度、应变、应力

记 $u$ 为 $n_{en}\times3$ 控制点位移矩阵。

**位移梯度:**
$$H = u^T \nabla_X N, \qquad H_{iJ} = \frac{\partial u_i}{\partial X_J}$$

**变形梯度 —— 用 $I + H$,不要用 $(P+u)^T\nabla_X N$:**
$$F = I + H$$

理由:$P^T\nabla_X N$ 在精确算术下恒等于 $I$,但浮点下只是近似 $I$,那点误差会污染 $F$。直接写 $I$ 更准。

**Green–Lagrange 应变 —— 必须用 $H$ 形式:**
$$E = \tfrac12\left(H + H^T + H^TH\right)$$

**不要用 $E = \frac12(F^TF - I)$。** 两式数学等价,但后者在小位移时是两个几乎相等的量相减,有效数字大量抵消。实测相对误差对比(有理数精确解为基准):

| $\|H\|$ 量级 | $\frac12(F^TF-I)$ | $\frac12(H+H^T+H^TH)$ |
|---|---|---|
| $10^{-2}$ | 3.6e-15 | 6.8e-17 |
| $10^{-4}$ | 4.2e-13 | 4.7e-17 |
| $10^{-6}$ | 1.0e-10 | 0 |
| $10^{-8}$ | 5.3e-09 | 1.6e-16 |
| $10^{-10}$ | 3.8e-07 | 0 |

朴素式误差随位移变小而线性恶化。**这对验证阶段尤其要紧** —— 会先跑小位移算例与线性解析解对比,恰好落在抵消最严重的区间,否则会看到"精度上不去"却查不出原因。

**第二类 PK 应力(SVK) —— 直接用张量形式:**
$$S = \lambda\,\mathrm{tr}(E)\,I + 2\mu E$$

$3\times3$ 矩阵的标量乘加,**不需要 Voigt、不需要逐分量写**。$E$ 与 $S$ 全程以 `Matrix3d` 存在,Voigt 只在做 $\tilde B_A^T C \tilde B_B$ 时才临时用到。

(换 neo-Hookean 也一样是张量运算:$S = \mu(I - C^{-1}) + \lambda\ln(\det F)\,C^{-1}$,其中 $C = F^TF$。)

**SVK 适用范围提醒:** SVK 在大压缩下失去椭圆性(非物理软化)。做验证、中等应变没问题,但不要用它跑大变形压缩算例然后困惑于牛顿不收敛。

---

### 3.4 应变–位移矩阵

**⚠️ 必须严格区分两个矩阵,不要混淆:**

$$\delta E^v = \tilde B_A\,\delta d_A, \qquad \tilde B_A = B_A F^T$$

- $B_A$($6\times3$)**只依赖** $G_A$,与 $F$、与位移无关 ⟹ 属于 **L1**,全程不变
- $\tilde B_A$($6\times3$)含 $F$ ⟹ 属于 **L2**,每次迭代重算

**Voigt 约定(必须固定,这是本类代码最高发的 bug):**
$$\{E\} = [E_{11},\ E_{22},\ E_{33},\ 2E_{12},\ 2E_{23},\ 2E_{13}]^T \quad(\text{工程剪应变,带因子 2})$$
$$\{S\} = [S_{11},\ S_{22},\ S_{33},\ S_{12},\ S_{23},\ S_{13}]^T \quad(\text{不带因子 2})$$

这样 $\{S\}^T\{\delta E\} = S:\delta E$ 才成立。若两者都带因子 2,剪切能量会被双重计数 —— 症状是"结果差一点点但收敛正常",极难查。

**$B_A$ 的显式形式:**
$$B_A=\begin{bmatrix} G_{A1} & 0 & 0\\ 0 & G_{A2} & 0\\ 0 & 0 & G_{A3}\\ G_{A2} & G_{A1} & 0\\ 0 & G_{A3} & G_{A2}\\ G_{A3} & 0 & G_{A1} \end{bmatrix}$$

这**就是小应变里最标准的 B 矩阵**。$F = I$ 时 $\tilde B_A = B_A$ 自动退化 —— 大变形的全部复杂性被隔离在 $F^T$ 这一个因子里。

**推导依据(供核对):** $\delta F = \delta d_A \otimes G_A$,代入 $\delta E = \frac12(F^T\delta F + \delta F^T F)$ 得
$$\delta E_{IJ} = \tfrac12\big(G_{AJ}F_{kI} + G_{AI}F_{kJ}\big)\delta d_{Ak}$$
按上述 Voigt 序展开,$F$ 的部分正好凑成 $F^T$ 右乘。

**本构矩阵 $C$(SVK,各向同性,$6\times6$):**
```
C[0:3, 0:3] : 对角 λ+2μ, 非对角 λ
C[3,3] = C[4,4] = C[5,5] = μ
其余为 0
```
$C$ 与 $x$ 和待求变量 $d$ 都无关 ⟹ **仅 SVK 下属于 L0,全程算一次**。

⚠️ **代码结构应按"$C$ 属于 L2"来组织**(即允许它每积分点每迭代变化),这样以后换 neo-Hookean(其 $C$ 依赖 $F$)不用动骨架。但仍在双重循环**之外**。

**$C$ 的特征值**为 $3\lambda+2\mu$、$2\mu$(二重)、$\mu$(三重),只要 $\mu>0$ 且体积模量为正即**对称正定**,这一点在 3.5 节会用到。

---

### 3.5 单元切线刚度矩阵

$$\boxed{K_{AB} = \int_\Omega \left[(G_A^T S\, G_B)\,I + F B_A^T C B_B F^T\right]\mathrm{d}\Omega}$$

**第一步:利用结合律把 $F$ 提出双重循环。** 这是关键恒等式:
$$F B_A^T C B_B F^T = (B_AF^T)^T C\,(B_BF^T) = \tilde B_A^T C\, \tilde B_B$$

原式中 $F$ 夹在双重循环里,每对 $(A,B)$ 要做两次 $3\times3$ 乘法;提出后每个 $A$ 只算一次 $\tilde B_A$。$\tilde B_A$ 必须在 L2 中**一次性算好存成局部数组**,否则会被重算 $d^3$ 次。

**第二步:大矩阵路线 + 对称性(两者不冲突)。**

把所有 $\tilde B_A$ 横向拼成 $\tilde B$($6 \times 3n_{en}$)。由于 $C$ 对称正定,做 Cholesky $C = LL^T$:
$$K_{\text{mat}} = \tilde B^TC\tilde B = \tilde B^TLL^T\tilde B = (L^T\tilde B)^T(L^T\tilde B) = Y^TY$$

$Y^TY$ 是标准的对称秩-$k$ 更新(syrk)形式,**同时获得"一半 flop"与"向量化稠密内核"**:

```
L0:  Eigen::LLT<Matrix<double,6,6>> llt(C);   L = llt.matrixL();   // 全程一次
L2:  Matrix<double,6,Dynamic> Y = L.transpose() * Btilde;          // 21n 乘(L^T 上三角)
L3:  Ke.selfadjointView<Eigen::Lower>().rankUpdate(Y.transpose(), w_tilde);
```

代价从 $6n^2$ 降到 $3n^2$($n = 3d^3$),已数值验证两条路线一致(误差 7e-15)。

⚠️ **若以后换到有软化的本构,$C$ 可能失去正定性,Cholesky 会失败。** 加 `if (llt.info() != Eigen::Success)` 回退到稠密 `Btilde.transpose() * C * Btilde`。

**第三步:几何项。** 记 $G$ 为 $n_{en}\times3$ 的 $\nabla_X N$,则几何刚度的标量部分是 $GSG^T$($n_{en}\times n_{en}$ 对称)。

$$K_{\text{geo}} = (G S G^T)\otimes I_3$$

**但不要显式构造 Kronecker 积再相加。** $n_{en}\times n_{en}$ 扩成 $3n_{en}\times3n_{en}$ 后,$9d^6$ 个位置里有 $6d^6$ 是零。正确做法:算出 $n_{en}\times n_{en}$ 标量矩阵,然后**只往每个 $3\times3$ 块的对角线上加** —— $d^6$ 次加法而非 $9d^6$。

几何项只占材料项的约 5.6%,**建议直接 `G * S * G.transpose()` 稠密算完事**,不必为它再做对称优化。

(若确实要优化:$S$ 是应力,受压时有负特征值,**不能 Cholesky**。可走特征分解 $S = \sum_k\lambda_k v_kv_k^T$,则 $GSG^T = \sum_k\lambda_k(Gv_k)(Gv_k)^T$,三次 rank-1 更新,已验证正确。收益 $3d^6 \to 1.5d^6$,即总量的 3%,不值得。)

**第四步:展开对称视图。** syrk 只填下三角,最后必须
```
Ke = Ke.selfadjointView<Eigen::Lower>();
```
⚠️ **忘记这一步是极隐蔽的 bug** —— 症状是求解器报错或收敛异常,但看不出是哪一步错的。

**内力向量(用张量形式,更便宜且不经过 Voigt):**
$$f_A^{\text{int}} = F S\, G_A \qquad (3\text{ 维向量})$$

**代价总结:** 朴素约 $216d^6$ 乘/单元 → 上述路线约 $28d^6$,约 7–8 倍。$d^6$ 这个量级本身动不了(需要 sum factorization,把 $d^9$ 降到 $d^7$,改动大,现阶段不划算)。

---

### 3.6 完整单元循环伪代码

```
// ===== L0:全程一次 =====
build univariate basis tables for 3 directions      // 2.3
build span[3][] arrays                              // 2.2
C = SVK constitutive matrix                         // 3.4
L = cholesky(C).matrixL()                           // 3.5
gauss points & weights on [-1,1]^3

// ===== 每个单元 =====
for e in elements:
    (e1,e2,e3) = unflatten(e)
    (s1,s2,s3) = (span[0][e1], span[1][e2], span[2][e3])
    A_base     = (s1-p1) + (s2-p2)*n1 + (s3-p3)*n1*n2
    局部→全局映射:按 2.2 强度削减法生成 IEN 视图
    gather P_e (nen×3), weights_e, u_e (nen×3)

    Ke.setZero();  fe.setZero();

    for q in gauss points (q1,q2,q3):
        // ---- L1(牛顿迭代中可缓存/不变) ----
        R, dN_dxi = NurbsShapeFunctions(...)          // 3.1
        J        = P_e.transpose() * dN_dxi           // 3.2
        assert(J.determinant() > 0)
        Jinv     = J.inverse()
        dN_dX    = dN_dxi * Jinv                      // = G,  nen×3
        w_tilde  = J.determinant() * detJparent * w_gauss
        build B_A for each A from dN_dX rows          // 3.4

        // ---- L2 ----
        H = u_e.transpose() * dN_dX                   // 3.3
        F = Matrix3d::Identity() + H
        E = 0.5*(H + H.transpose() + H.transpose()*H)
        S = lambda*E.trace()*Matrix3d::Identity() + 2*mu*E
        for A: Btilde.block(0,3A,6,3) = B_A * F.transpose()
        Y = L.transpose() * Btilde

        // ---- L3 ----
        Ke.selfadjointView<Lower>().rankUpdate(Y.transpose(), w_tilde)
        Kg_scalar = dN_dX * S * dN_dX.transpose()     // nen×nen
        add w_tilde*Kg_scalar onto the 3 diagonal entries of each 3×3 block
        for A: fe.segment(3A,3) += w_tilde * (F * S * dN_dX.row(A).transpose())

    Ke = Ke.selfadjointView<Lower>()                  // ⚠️ 不要漏
    scatter Ke, fe into global K, f via IEN
```

---

## 第四部分:验证测试(按发现 bug 的效率排序)

**强烈建议每写完一层就跑对应的测试,不要全部写完再一起调试。**

### 4.1 $P^T\nabla_X N = I$(打通第一层时最该先跑的)

因为参考坐标 $X$ 本身就是用同一组基函数插值出来的(等参概念):
$$X_i(\xi) = \sum_A P_{Ai}R_A(\xi) \ \Longrightarrow\ \sum_A P_{Ai}\frac{\partial R_A}{\partial X_J} = \frac{\partial X_i}{\partial X_J} = \delta_{iJ}$$

数值上要成立,必须 $J = P^T\nabla_\xi N$ 算对、$J^{-1}$ 求对、$\nabla_XN = \nabla_\xi N\cdot J^{-1}$ 的乘法顺序与转置都没错 —— 三步任一处出错结果就不是 $I$。且对**任意**积分点、任意几何都成立,不需构造特殊算例。

配套:$\sum_A R_A = 1$(检验基函数值)、$\sum_A \nabla_X R_A = 0$(检验梯度)。

### 4.2 形函数退化测试

所有权重 $=1$ 时,$R$ 必须退化成 B 样条基函数($\omega$ 全为 1 时 $R = W/\text{Wsum}$ 恰为 $N_aM_bL_c$),且 $\sum R = 1$、$\sum \partial R/\partial\xi = 0$。

⚠️ **有理导数路径只有在权重非 1 的几何上才被真正走到**(如厚壁圆筒),必须专门构造这样的算例测试,否则商法则的 bug 会一直潜伏。

### 4.3 $\det J > 0$ 断言

在所有积分点上断言。参数化退化时的第一道防线。

### 4.4 刚体转动 ⟹ $E \equiv 0$

施加纯刚体转动位移场,$E$ 必须为机器精度零。抓 $E$ 和 $F$ 的错。

### 4.5 $F = I$ 时 $\tilde B_A$ 退化成 $B_A$

即小应变 B 矩阵。零成本的自检。

### 4.6 两条内力路径必须一致

$$B_A^T\{S\}\ \text{(经 } \tilde B_A\text{)} \quad == \quad F S G_A\ \text{(张量形式)}$$

这是检验 $B_A$ 与 Voigt 约定(尤其那个因子 2)是否自洽的**最好办法**。建议写成断言常驻。

### 4.7 切线矩阵的有限差分校验(最关键)

对残量向量做中心差分,与 $K_T$ 逐项比对,**收敛阶应看到二阶**。

这是唯一能同时抓住 $B_A$、几何刚度、Voigt 约定、$C$ 矩阵所有错误的检查。**不通过就不要往下走。**

### 4.8 syrk 版本与朴素版本逐项比对

验证阶段先用最朴素的 `Ke = Btilde.transpose() * C * Btilde` 跑通并过 4.7,**再**换 syrk 版本,两版结果逐项比对。

### 4.9 后续物理验证

- 单轴拉伸对比解析解
- 纯剪切(专门抓因子 2 的错)
- $h$-refinement 收敛率
- 牛顿二次收敛(切线刚度正确性最灵敏的指标)

---

## 附录:实现优先级建议

1. **数据结构迁移**(2.1)—— 先把 `DersBasisFuns` 输出改成扁平/Eigen 类型,`ControlPointArray` 改成 Eigen 矩阵
2. **张量积空间 + span 数组 + 惰性 IEN**(2.2)—— 对照小例子手工核对连接表
3. **单变量基函数表**(2.3)
4. **NURBS 形函数 + 参数域梯度**(3.1)—— 跑 4.2
5. **两层 Jacobian + 物理域梯度**(3.2)—— 跑 4.1、4.3
6. **$H, F, E, S$**(3.3)—— 跑 4.4
7. **$B_A$、$\tilde B_A$、$C$**(3.4)—— 跑 4.5、4.6
8. **切线刚度(先朴素版)**(3.5)—— 跑 4.7
9. **切线刚度(syrk 版)**—— 跑 4.8
10. 装配与牛顿迭代(本次不实现)

## 附录:参考文献

- Piegl & Tiller, *The NURBS Book*, 2nd ed. —— 算法 A2.1/A2.2/A2.3
- Hughes, Cottrell & Bazilevs, *Isogeometric Analysis: Toward Integration of CAD and FEA* —— Appendix 3.A 形函数例程、IEN/INC 数组
- Bonet & Wood, *Nonlinear Continuum Mechanics for Finite Element Analysis*, ch. 5 & 9 —— $\tilde B_A = B_AF^T$ 分解、几何刚度
- Belytschko et al., *Nonlinear Finite Elements for Continua and Structures* —— Total Lagrangian 章节



---

# IGA 非线性固体求解器 —— 第二阶段规格:外力、装配、牛顿迭代与测试

> 本文档接续《IGA 非线性固体求解器 —— 实现规格说明》。 **前置假设:** 单元层已完成并通过验证 —— 能对任一单元、任一积分点算出 $R_A$、$G_A=\nabla_XR_A$、$J$、$\tilde w$、$H$、$F$、$E$、$S$、$B_A$、$\tilde B_A$, 并组装出局部切线刚度 $K_e$($3d^3\times3d^3$)与局部内力 $f_e^{\text{int}}$($3d^3$)。 **尤其:单元级有限差分校验必须已通过。** 若未通过,不要进入本阶段 —— 全局层的 bug 会与单元层的 bug 互相掩盖,极难分离。

------

## 第零部分:自由度编号约定(先固定,后面全部依赖它)

$$\text{全局自由度编号} \quad I = 3A + i,\qquad A \in [0, N_{\text{cp}}), \ i \in {0,1,2}$$

其中 $A$ 是全局基函数(控制点)编号,$i$ 是空间分量。局部同理 $\ 3a+i$。

**必须与单元层的 $\tilde B$ 布局一致。** 前一阶段中 $\tilde B$ 是把每个 $A$ 的 $\tilde B_A$($6\times3$)横向拼接,所以局部自由度就是 `3a+i`,全局必然是 `3*IEN(a,e)+i`。

⚠️ **若误用 $I = A + i,N_{\text{cp}}$(按分量分块)的排法**,矩阵仍对称、求解器仍收敛、残量仍下降 —— **但物理答案是错的**。这是最难发现的一类 bug,务必在装配处写死并加注释。

总自由度 $N_{\text{dof}} = 3 n_1 n_2 n_3$。

------

## 第一部分:外力项

### 1.1 概念前提:控制点不是物理节点

$R_A$ **不满足** Kronecker delta 性质($R_A(\boldsymbol\xi_B)\neq\delta_{AB}$)。因此**不能把集中力"加到控制点上"**,必须走一致载荷向量:

$$\mathbf{f}^{\text{ext}}*A = \underbrace{\int*{\Omega_0} R_A,\rho_0\mathbf{g},\mathrm{d}\Omega_0}*{\text{体力}} + \underbrace{\int*{\Gamma_0} R_A,\mathbf{t}*0,\mathrm{d}\Gamma_0}*{\text{面力}}$$

**本阶段只实现体力项。** 面力需要独立的面积分机制(2D Gauss 积分 + 面 Jacobian $|\mathbf{a}_1\times\mathbf{a}_2|$,不是 $\det J$),与牛顿迭代的 bug 混在一起很难查,推迟到下一阶段。

### 1.2 体力向量的计算

复用**现有的单元循环和积分点**,几乎零额外成本:

```
for e in elements:
    fe_ext.setZero()                          // 3*nen
    for q in gauss points:
        R, w_tilde = (单元层已有的量)
        for a in 0..nen-1:
            fe_ext.segment(3*a, 3) += w_tilde * R[a] * (rho0 * g)
    scatter fe_ext into global f_ext via IEN
```

其中 $\mathbf{g}$ 是 `Vector3d` 重力加速度向量,$\rho_0$ 是**参考构形**密度。

### 1.3 【关键】重力是死载荷 ⟹ $\mathbf{f}^{\text{ext}}$ 只算一次

重力在 TL 格式下作用于参考构形、方向固定、大小固定,**不依赖变形**。因此:

- $\mathbf{f}^{\text{ext}}$ 在**牛顿迭代开始之前算一次**,全程不变
- **没有载荷刚度项**,$\mathbf{K}*T = \partial\mathbf{f}^{\text{int}}/\partial\mathbf{d}$ 即为你已有的 $K*{AB}$ 装配结果
- $\mathbf{K}_T$ 保持对称,前一阶段的 syrk 优化继续有效

⚠️ 对比:**压力载荷/随动载荷**(压强沿当前法向、作用于当前面积)依赖变形,需要额外载荷刚度 $\mathbf{K}^{\text{ext}} = -\partial\mathbf{f}^{\text{ext}}/\partial\mathbf{d}$,**且该矩阵非对称** —— 会破坏对称性、使 syrk 与 `SimplicialLDLT` 都失效。**验证阶段绝不引入。**

------

## 第二部分:装配

### 2.1 使用 triplet 列表(验证阶段的正确选择)

```cpp
std::vector<Eigen::Triplet<double>> triplets;
triplets.reserve(n_elem * (3*nen) * (3*nen));   // 先 reserve,避免反复扩容

for e in elements:
    compute Ke, fe_int
    for a in 0..nen-1:
      for i in 0..2:
        I = 3*IEN(a,e) + i;
        f_int(I) += fe_int(3*a+i);
        for b in 0..nen-1:
          for j in 0..2:
            J = 3*IEN(b,e) + j;
            triplets.emplace_back(I, J, Ke(3*a+i, 3*b+j));

K.setFromTriplets(triplets.begin(), triplets.end());   // 自动累加重复项
```

`setFromTriplets` 会自动把同一 $(I,J)$ 的多个贡献相加,这正是装配需要的语义。

**不要在此阶段做 CSR 索引预计算优化。** 那个优化(稀疏模式全程固定 ⟹ 索引表算一次、后续所有牛顿迭代直接按索引累加)在牛顿迭代下确实划算,能把散射开销砍 3–10 倍,但**等功能跑通、测试全绿之后再做**。

### 2.2 每次迭代都要重装

$\mathbf{K}_T$ 和 $\mathbf{f}^{\text{int}}$ 都依赖当前 $\mathbf{d}$,每次牛顿迭代都要重新装配。$\mathbf{f}^{\text{ext}}$ 不用(第 1.3 节)。

建议把装配写成一个函数 `assemble(d, /*out*/ K, /*out*/ f_int)`,牛顿循环里调用。这样也方便做全局有限差分校验。

### 2.3 $\mathbf{f}^{\text{int}}$ 用同一套 IEN 散射

容易漏。$\mathbf{f}^{\text{int}}$ 和 $\mathbf{K}_T$ 的行索引必须完全一致。

------

## 第三部分:边界条件

### 3.1 IGA 中齐次 Dirichlet 可以精确施加

由于 $R_A$ 不是节点插值函数,一般的 Dirichlet 条件在 IGA 中不能想当然地施加。但**开(clamped)结点向量**下情况可控:在 $\xi=0$ 面上,只有 $i_1=0$ 那一层控制点的基函数非零,且它们在该面上构成单位分解($\sum R_A|_{\xi=0}=1$)。因此:

| 边界条件                            | 能否精确施加 | 做法                                      |
| ----------------------------------- | ------------ | ----------------------------------------- |
| 齐次 $\mathbf{u}=\mathbf{0}$        | **可以**     | 该层全部控制变量置 0                      |
| 常值 $\mathbf{u}=\bar{\mathbf{u}}$  | **可以**     | 该层全部置 $\bar{\mathbf{u}}$(靠单位分解) |
| 一般函数 $\mathbf{u}=\mathbf{g}(X)$ | **不能**     | 需 $L^2$ 投影 / 最小二乘                  |

**本阶段的测试只用第一种。**

### 3.2 受约束自由度的确定

固定 $\xi=0$ 面 ⟹ 所有 $i_1 = 0$ 的控制点,即 $$A = 0 + i_2 n_1 + i_3 n_1 n_2, \qquad i_2\in[0,n_2),\ i_3\in[0,n_3)$$ 每个 $A$ 的三个分量全部约束。共 $3n_2n_3$ 个受约束自由度。

**这足以消除全部 6 个刚体模态**(整个面在三个方向都被固定)。若不足,$\mathbf{K}_T$ 奇异,`SimplicialLDLT` 会返回 `NumericalIssue` 或给出垃圾结果 —— **必须检查求解器返回值并断言**。

### 3.3 施加到线性系统:让初始猜测满足 BC

**推荐做法**:令 $\mathbf{d}^0$ 在受约束自由度上已等于规定值(齐次情形即 0)。于是**所有迭代中 $\Delta\mathbf{d}$ 在受约束自由度上恒为 0**,处理极简:

```
对每个受约束自由度 I:
    K 的第 I 行、第 I 列全部置零        // 行和列都要,否则破坏对称性
    K(I,I) = 1.0
    rhs(I)  = 0.0
```

⚠️ **必须行、列都置零。** 只置行会破坏对称性,`SimplicialLDLT` 直接失效。

Eigen 稀疏矩阵按行/列置零的实用做法:装配时**直接跳过**任一索引受约束的 triplet,然后单独为每个受约束自由度加一个 $(I,I,1)$ 的 triplet。这比事后修改稀疏矩阵干净得多,也更快。

------

## 第四部分:牛顿迭代

### 4.1 残量与迭代

$$\mathbf{R}(\mathbf{d}) = \mathbf{f}^{\text{int}}(\mathbf{d}) - \mathbf{f}^{\text{ext}}$$

每次迭代解 $$\mathbf{K}_T,\Delta\mathbf{d} = -\mathbf{R},\qquad \mathbf{d} \mathrel{+}= \Delta\mathbf{d}$$

### 4.2 收敛判据

$$\frac{|\mathbf{R}|*{\text{free}}}{|\mathbf{f}^{\text{ext}}|*{\text{free}}} < \texttt{tol}\qquad(\texttt{tol} = 10^{-9}\ \text{典型})$$

⚠️ **范数只在自由自由度上取。** 受约束自由度上的残量是支反力,不为零且量级很大,算进去会导致**永远不收敛** —— 这是经典 bug,而且症状(残量停在某个不小的值)很容易被误诊为切线矩阵有错。

若 $|\mathbf{f}^{\text{ext}}|=0$(纯位移加载),回退到绝对判据或用 $|\Delta\mathbf{d}|/|\mathbf{d}|$。

### 4.3 【最重要】必须观察到二次收敛

典型序列:

```
iter 1:  rel_res = 1.0e+00
iter 2:  rel_res = 3.7e-03
iter 3:  rel_res = 1.2e-07
iter 4:  rel_res = 8.4e-15     <- 收敛
```

残量指数**大致翻倍**($10^{-3}\to10^{-7}\to10^{-14}$)即为二次收敛。

**这是检验 $\mathbf{K}_T$ 正确性最灵敏的手段。** 若只看到线性收敛(每次降一个固定倍数,如 $10^{-1}\to10^{-2}\to10^{-3}$),说明**切线矩阵有错而残量向量对**。最常见原因:

- 几何刚度项漏了或系数错
- Voigt 因子 2 用错(${E}$ 带 2、${S}$ 不带)
- $\tilde B_A$ 里 $F$ 的转置反了(写成 $B_AF$ 而非 $B_AF^T$)
- syrk 后忘记 `Ke = Ke.selfadjointView<Lower>()`

⚠️ **这种情况下解最终仍会收敛到正确答案** —— 所以不看收敛阶就发现不了。**必须把每次迭代的残量打印出来。**

### 4.4 载荷分级(robustness)

大变形下从零直接牛顿往往不收敛。把总载荷分成 $N_{\text{step}}$ 步:

```
for n in 1..N_step:
    f_ext_current = (n / N_step) * f_ext_total
    Newton solve, 初始猜测 = 上一步的解
```

线性/小变形算例用 $N_{\text{step}}=1$;强非线性用 10–20 步。

### 4.5 每次迭代后检查 $\det F > 0$

牛顿的步长可能大到把单元翻转过来($\det F \le 0$),之后一切都是垃圾。检测到就把 $\Delta\mathbf{d}$ 减半重试(简易线搜索),重试超过若干次则报错退出。

### 4.6 求解器选择

用 `Eigen::SimplicialLDLT`,**不要用 `SimplicialLLT`** —— 过了临界点 $\mathbf{K}_T$ 可能不定,Cholesky 会失败。调试期可用 `SparseLU` 交叉验证同一系统的解。

**每次求解后必须检查 `solver.info() == Eigen::Success`。**

### 4.7 完整流程伪代码

```
// ---- 一次性 ----
build space, span arrays, univariate basis tables
C = SVK matrix;  L = cholesky(C)
f_ext_total = assemble_body_force()          // 死载荷,只算一次
constrained_dofs = collect_xi0_face_dofs()
d = 0                                        // 已满足齐次 BC

// ---- 载荷分级 ----
for n in 1..N_step:
    f_ext = (n/N_step) * f_ext_total
    for it in 1..max_iter:
        assemble(d) -> K_T, f_int
        R = f_int - f_ext
        rel_res = norm(R restricted to free dofs) / norm(f_ext restricted to free dofs)
        print(it, rel_res)                   // ⚠️ 必须打印,用于判断收敛阶
        if rel_res < tol: break
        apply_bc(K_T, -R)                    // 行列置零 + 对角 1
        solve K_T * delta = -R
        assert(solver.info() == Success)
        d += delta
        assert(all det F > 0)                // 否则 delta 减半重试
    assert(converged)
```

------

## 第五部分:测试算例

### 5.0 对"悬臂梁自重弯曲"方案的评价

**方案本身很合理**,理由:

- 一端固定 ⟹ 只用**齐次** Dirichlet,IGA 中可精确施加(3.1 节)
- 只有体力 ⟹ 不引入面积分机制,复用现有积分点
- 重力是**死载荷** ⟹ $\mathbf{f}^{\text{ext}}$ 常量、无载荷刚度、$\mathbf{K}_T$ 保持对称
- 有 Euler–Bernoulli 解析解可比
- 弯曲状态会真正激活 $B_A$ 的全部行(含剪切行),比纯拉伸更充分

**但有三个盲点必须知道:**

| 盲点                                                         | 后果                                                         | 对策                                                         |
| ------------------------------------------------------------ | ------------------------------------------------------------ | ------------------------------------------------------------ |
| 直梁是长方体,**所有权重 $\omega=1$**                         | 此时 $\text{Wsum}\equiv1$、$\text{Wsum}*\xi\equiv0$,**有理导数的商法则路径形同虚设** —— $\text{Wsum}*\xi$ 项若有 bug 完全看不出来 | 后续必须补一个**非单位权重**几何(如四分之一圆环/厚壁圆筒)的算例 |
| 直梁的几何映射是**仿射的**,$J = \mathrm{diag}(L/n_1', b, h)$ 为常量 | $J^{-1}$ 若只在特定情形下正确、或 $\nabla_XN = \nabla_\xi N\cdot J^{-1}$ 的乘法顺序错但恰好对角情形下无碍,**都不会暴露** | 同上,需曲边几何补测                                          |
| 细长实体梁有**剪切锁定**风险                                 | 位移偏小,会被误诊为公式错                                    | 用 $p=3$、$L/h=10$(不要更细长)、沿长度方向足够单元数         |

**因此建议先加一个更简单、判据更锐利的算例(5.1),再做弯曲算例(5.2)。** 5.1 能把"装配 / BC / 求解 / 牛顿"这条管路的 bug 与"弯曲力学"的 bug 彻底分离。

------

### 5.1 【首先做这个】测试 A:轴向自重杆

**唯一差别是重力方向沿轴向**,代码零改动,但获得**解析精确解**。

**几何与参数:**

```
长方体:  L = 10 (沿 X/ξ),  b = 1 (沿 Y/η),  h = 1 (沿 Z/ζ)
degree:  p1 = p2 = p3 = 3
单元:    8 × 1 × 1   (N_dof = 3·11·4·4 = 528)
材料:    E = 1.0e7,  nu = 0     <- 必须用 nu=0,见下
重力:    g = (g_x, 0, 0),  rho0*g_x = 1.0
BC:      xi=0 面全约束
```

**为什么必须 $\nu=0$:** $\nu=0$ ⟹ $\lambda=0$,于是横向应力 $S_{yy}=S_{zz}=0$ 自动满足自由表面条件,**三维解退化为严格的一维解**,横向位移 $u_y=u_z\equiv0$。有了 $\nu\neq0$ 就有横向收缩,不再有闭式解。

**精确线性解:** $$u_x(X) = \frac{\rho_0g}{E}\left(LX - \frac{X^2}{2}\right),\qquad u_y=u_z=0$$

**这个解是 $X$ 的二次多项式** ⟹ $p\ge2$ 的 NURBS 空间能**精确表示**它 ⟹ **离散误差为零**。

**判据:** $\rho_0g=1$ 时最大应变 $\varepsilon_{\max}=\rho_0gL/E=10^{-6}$,几何非线性带来的模型偏差为 $O(\varepsilon)=10^{-6}$。因此

> **计算解与上述线性精确解的相对误差应 $< 10^{-5}$;超出即有 bug。**

$u_x(L) = 5.0\times10^{-6}$。同时断言 $\max|u_y|,\max|u_z| < 10^{-12}\cdot L$。

**这是一个极其锐利的测试** —— 离散误差为零,任何装配错误、自由度排序错误、BC 施加错误都会立刻表现为明显偏差。

**测试 A2:同一算例,大载荷,对比非线性精确解**

一维 SVK 下有闭式(半解析)解:第一类 PK 应力 $P=FS$ 满足 $\mathrm{d}P/\mathrm{d}X=-\rho_0g$、$P(L)=0$,故 $$P(X)=\rho_0g(L-X),\qquad \frac{E}{2}F(F^2-1)=\rho_0g(L-X)$$ 对每个 $X$ 解这个**三次方程**得 $F(X)$,再积分 $u(X)=\int_0^X(F-1),\mathrm{d}X'$。

取 $\rho_0g = 10^4$($\varepsilon_{\max}=10^{-2}$,与线性解偏差 $\approx1%$)。参考值(已用高精度求根+积分算出,可直接硬编码):

| $X$  | $u_x^{\text{exact}}$   |
| ---- | ---------------------- |
| 0.0  | 0.0                    |
| 2.5  | 2.159257958214e-02     |
| 5.0  | 3.707162832131e-02     |
| 7.5  | 4.639251893304e-02     |
| 10.0 | **4.950974524113e-02** |

(线性解会给出 $u_x(L)=5.0\times10^{-2}$,差 0.99% —— 足以区分"非线性项是否真的起作用")

此时 $F(X)$ 含三次方根,**不是多项式,故离散误差不为零**。判据改为:$8\times1\times1$ 网格相对误差应在 $10^{-5}$ 量级以内,且**加密网格时误差应下降**。

**这个算例是验证 $E$、$S$、$\tilde B_A$、几何刚度的最好载体** —— 有精确解、又真正激活了非线性项。

------

### 5.2 测试 B:悬臂梁自重弯曲(线性区)

**几何与参数:** 同 5.1,但

```
单元:    16 × 2 × 2   (N_dof = 3·19·5·5 = 1425)
材料:    E = 1.0e7,  nu = 0      (先用 nu=0 便于与梁理论对比,避免反曲效应)
重力:    g = (0, 0, -g_z),  rho0*g_z = 1.0
BC:      xi=0 面全约束
```

**Euler–Bernoulli 参考解:** $$w_{\text{tip}} = \frac{qL^4}{8EI},\quad q=\rho_0g,bh,\ I=\frac{bh^3}{12} \ \Longrightarrow\ w_{\text{tip}}=\frac{1.5,\rho_0g,L^4}{E,h^2}$$

$\rho_0g=1$ 时 $w_{\text{tip}}^{\text{EB}} = 1.5\times10^{-3}$,$w/L=1.5\times10^{-4}$(纯线性区),$\varepsilon_{\max}=3\times10^{-5}$。

**判据:** 与 EB 解**相差几个百分点内即可**。不要期望机器精度 —— EB 忽略剪切变形、且是小挠度理论,而三维实体解含剪切与端部效应。$L/h=10$ 时差异通常 1–3%。

⚠️ **若偏差达到 10% 以上且偏小,首先怀疑剪切锁定**(而非公式错)。验证方法:加密网格或提高 $p$,若位移显著增大则是锁定。

**⚠️ 后处理陷阱:** "梁端位移"不能取 $\max_A|\mathbf{d}_A|$。$\mathbf{d}_A$ 是**系数**不是节点位移,可以超出真实位移范围(同控制点可位于几何之外)。必须求值: $$\mathbf{u}(\boldsymbol\xi^*) = \sum_A R_A(\boldsymbol\xi^*),\mathbf{d}_A$$ 在 $\boldsymbol\xi^*=(1, 0.5, 0.5)$(梁端中心)处求值。**这是 `FindSpan` 终于派上用场的地方。**

**载荷步:** $N_{\text{step}}=1$。牛顿应在 2 次迭代内收敛(几乎线性)。

------

### 5.3 测试 C:悬臂梁自重弯曲(非线性区)

同 5.2,但 $\rho_0g = 500$:

| 量                           | 值                 |
| ---------------------------- | ------------------ |
| EB 线性预测 $w_{\text{tip}}$ | $7.5\times10^{-1}$ |
| $w/L$                        | $7.5\times10^{-2}$ |
| $\varepsilon_{\max}$         | $1.5\times10^{-2}$ |

**预期行为:**

- 计算得到的 $w_{\text{tip}}$ 应**明显小于** EB 线性预测(几何刚化效应)。若计算值 $\approx$ 线性预测,说明**非线性项没起作用**(几何刚度漏了、或 $E$ 里的 $H^TH$ 项漏了)。
- 牛顿需 4–5 次迭代,**必须观察到二次收敛**(4.3 节)。这是本测试的主要目的。
- $N_{\text{step}}=1$ 应该能收敛;若不收敛,试 $N_{\text{step}}=5$。

**注:** $\varepsilon_{\max}=1.5%$ 对 SVK 是安全的。**不要把载荷加到应变超过 ~5%** —— SVK 在大压缩下失去椭圆性(非物理软化),会出现"牛顿怎么都不收敛"而其实是本构模型的问题,不是代码 bug。

有用的关系式(便于设计其他载荷水平): $$\frac{w}{L} = \frac{1}{2}\cdot\frac{L}{h}\cdot\varepsilon_{\max}$$ 即细长比换取"低应变下的大转动"。但 $L/h$ 越大剪切锁定越严重,$L/h=10\sim20$ 是折中。

------

### 5.4 补充检验(低成本,建议全部加上)

| 检验                      | 内容                                                         | 抓什么                                 |
| ------------------------- | ------------------------------------------------------------ | -------------------------------------- |
| **全局 FD 校验**          | 对 $\mathbf{R}(\mathbf{d})$ 做中心差分,与装配出的 $\mathbf{K}_T$ 逐项比(小规模网格,如 $2\times1\times1$) | 装配、自由度映射、BC —— **最值得先做** |
| **$\mathbf{K}_T$ 对称性** | $|\mathbf{K}-\mathbf{K}^T|/|\mathbf{K}| < 10^{-12}$          | syrk 未展开、只置零了行没置零列        |
| **零载荷 ⟹ 零位移**       | $\mathbf{f}^{\text{ext}}=0$ 时 $\mathbf{d}=0$,残量首次即为 0 | 内力向量在 $\mathbf{d}=0$ 时应恒为零   |
| **支反力平衡**            | 受约束自由度上的残量之和 $= -\int\rho_0\mathbf{g}$(总重力),三个分量分别核对 | 体力向量、装配、BC                     |
| **刚体转动叠加**          | 给收敛解叠加一个大刚体转动,重算残量应仍为零(客观性)          | $E$、$F$ 的正确性                      |
| **$h$-refinement**        | 加密网格,位移应单调收敛                                      | 整体一致性;也用于区分锁定              |
| **纯剪切算例**            | 后续单独构造                                                 | 专抓 Voigt 因子 2                      |

------

## 第六部分:bug 症状对照表

| 症状                                   | 最可能的原因                                                 |
| -------------------------------------- | ------------------------------------------------------------ |
| 残量停在某个不小的值不再下降           | 收敛判据把**受约束自由度**算进范数了(4.2)                    |
| 只有线性收敛,但最终答案正确            | $\mathbf{K}_T$ 有错、$\mathbf{f}^{\text{int}}$ 对:几何刚度漏项 / Voigt 因子 2 / $B_AF$ 写成 $B_AF^T$ 之外的形式 / syrk 未展开 |
| 求解器报 `NumericalIssue`              | 刚体模态未完全约束;或 BC 只置零了行未置零列;或用了 `LLT` 而非 `LDLT` |
| 答案量级完全不对(差几个数量级)         | 自由度排序不一致(第零部分);或漏了 $\det J_{\text{parent}}$   |
| 答案对称性错(梁向侧向弯)               | 自由度排序 `3A+i` 与 `A+i*N` 混用                            |
| 位移系统性偏小 10%+                    | 剪切锁定(加密/提 $p$ 验证),或积分点数不足                    |
| 位移 $\approx$ 线性预测,大载荷下也如此 | 非线性项失效:$E$ 里漏 $H^TH$,或几何刚度整块漏了              |
| 小载荷下精度上不去($10^{-6}$ 量级)     | $E$ 用了 $\frac12(F^TF-I)$ 而非 $\frac12(H+H^T+H^TH)$        |
| 有理几何上错、直梁上对                 | 商法则 $\text{Wsum}_\xi$ 项有 bug(直梁 $\omega\equiv1$ 掩盖了它) |
| 后处理位移比解析解大很多               | 误把 $\max                                                   |

------

## 附录:实施顺序

1. **体力向量**(第一部分)—— 单独测:$\sum_A \mathbf{f}^{\text{ext}}_A$ 应等于 $\rho_0\mathbf{g}\times$ 总体积
2. **装配**(第二部分)—— 测 $\mathbf{K}_T$ 对称性
3. **BC 施加**(第三部分)—— 测 BC 后仍对称、非奇异
4. **全局 FD 校验**(5.4)—— 小网格,**通过后才继续**
5. **牛顿迭代**(第四部分)—— 打印每次迭代残量
6. **测试 A**(5.1)—— 判据 $<10^{-5}$,最锐利
7. **测试 A2**(5.1)—— 对比非线性精确解
8. **测试 B**(5.2)—— 对比 EB,几个百分点内
9. **测试 C**(5.3)—— **确认二次收敛**
10. **补充检验**(5.4)全部跑一遍
11. 后续:面力积分、非单位权重几何(补 5.0 的盲点)、CSR 索引优化

## 附录:参考文献

- Bonet & Wood, *Nonlinear Continuum Mechanics for Finite Element Analysis* —— 第 9 章牛顿求解、载荷刚度
- Belytschko et al., *Nonlinear Finite Elements for Continua and Structures* —— Total Lagrangian、收敛判据
- Hughes, Cottrell & Bazilevs, *Isogeometric Analysis* —— 第 3 章 Dirichlet 条件在 IGA 中的施加、开结点向量性质
- Timoshenko, *Strength of Materials* —— Euler–Bernoulli 悬臂梁解
