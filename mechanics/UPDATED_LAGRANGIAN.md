# Updated Lagrangian 单元积分

`UpdatedLagrangianElement.*` 与 `TotalLagrangianElement.*` 并列，分别在当前构形和
参考构形积分。两者返回相同的 `ElementResult`：交错自由度布局下的内力和一致切线。
全局散射、约束、静力 Newton 和 generalized-alpha 算法共用。

单元入口分别为 `computeTotalLagrangianElementResponse()` 和
`computeUpdatedLagrangianElementResponse()`。参考梯度统一命名为 `referenceGradients`，
当前梯度为 `spatialGradients`；`referenceJacobian` 和 `currentJacobian` 分别表示
`dX/dxi` 和 `dx/dxi`。TL 的 `materialB` 和 UL 的 `spatialB` 分别对应材料应变和
空间对称梯度算子。

推导对应 Bonet & Wood 第 9 章：9.2 节的运动学、式 (9.15b) 的内力、
式 (9.41) 的材料刚度和式 (9.44) 的初应力刚度。

## 四步积分

1. `evaluateReferenceGradients()` 接收 NURBS 形函数及参数梯度，计算参考梯度。
   不构造 TL 的参考 B 矩阵，也不保存参考逆 Jacobian。
2. `evaluateSpatialPoint()` 根据当前控制点计算当前梯度、积分权重和总运动学量。
3. `HyperelasticMaterial::evaluateSpatial()` 返回 Cauchy 应力和空间切线。
4. `computeUpdatedLagrangianElementResponse()` 在当前构形积分内力及两部分刚度。

当前控制点 `x_A = X_A + u_A` 在积分循环外构造。这里的 `u_A` 始终是相对初始
构形的总位移，不是本次 Newton 增量，也不在迭代后重置材料参考构形。

## 三个映射

用 `xi` 表示 NURBS 参数域坐标，用 `eta` 表示 Gauss 母单元坐标。代码区分：

```text
eta -> xi -> X
          -> x

A0 = dX/dxi
At = dx/dxi
F  = dx/dX = At * inverse(A0)
J  = det(F) = det(At) / det(A0)
```

参数域到参考构形的映射在 UL 中仍用于求参考梯度和总变形梯度，不能完全省略。
UL 不需要保存完整参考映射，但保留其行列式，为体积比和后续平均体积算法提供数据。

项目中形函数梯度按行存储，因此：

```text
G_X = G_xi * inverse(A0)
G_x = G_xi * inverse(At) = G_X * inverse(F)
F   = I + u_e^T * G_X
```

UL 的完整积分权重为：

```text
w_current = det(At) * det(dxi/deta) * w_gauss
          = J * w_reference
```

使用 `det(At)` 后不能再乘一次 `J`。`det(dxi/deta)` 是独立的母单元到参数域因子，
也只能乘一次。参考几何和当前几何均要求正 Jacobian。

## 材料及 Voigt 约定

项目统一使用 `[11,22,33,12,23,13]`，应变/变形率的剪切分量乘 2，应力不乘 2。
这与书中的最后两项顺序不同。

`SpatialMaterialResponse::spatialTangent` 表示
`c = push_forward(C) / J` 的矩阵形式，搭配独立的初应力刚度使用。
它不是普通 Cauchy 应力对应变的导数，也不是带应力修正的 Jaumann 切线。

`evaluate()` 和 `evaluateSpatial()` 均为纯虚函数，材料必须显式实现两个入口。
可以只写一套构形公式，再通过 `pushForwardMaterialResponse()` 或
`pullBackMaterialResponse()` 推导另一套响应。SVK 在 TL 入口使用公式，空间入口调用推前。
两个入口的 `strainEnergyDensity` 都保持为单位参考体积上的能量。

```text
f_a      += w_current * sigma * grad_x(N_a)
K_mat    += w_current * B_x^T * D_spatial * B_x
K_geo_ab += w_current * (grad_x(N_a)^T * sigma * grad_x(N_b)) * I3
```

本次实现只包含局部本构对应的这两部分切线；mean-dilatation 的平均压力、内力修改
及额外 dilatational tangent 需要在之后的单元层扩展中一起加入。

## 选择积分形式

直接装配时传入第四个参数：

```cpp
const auto system = iga::assembleInternalSystem(
    space, displacements, material, iga::ElementFormulation::UpdatedLagrangian);
```

静力和动力学分别通过 `NewtonSettings::elementFormulation` 和
`DynamicsSettings::elementFormulation` 选择，均默认 TL。动力学选择会传到初始加速度
和每一步有效切线；质量矩阵及死载荷仍可在参考域一次装配。

## 验证

`tests/iga_updated_lagrangian_tests.pro` 覆盖：

- 各向异性四阶张量的独立逐指标推前，检查所有剪切分量。
- 非单位参数区间、多 span、弯曲有理 NURBS 几何上的三个映射和积分权重。
- 非均匀变形下 SVK 与测试用 Neo-Hookean 的 TL/UL 单元内力和切线对照。
- UL 内力的独立中心有限差分切线、刚体运动和无效构形拒绝。
- 全局装配及静力、动力学积分形式选择和解的对照。
