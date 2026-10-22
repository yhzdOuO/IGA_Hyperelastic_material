# 动力学调用逻辑：从程序入口到时间步迭代

程序默认以 `--mode static` 运行静力学；本文描述传入 `--mode dynamic`
后实际走到的动力学路径。静力模式直接调用
`solveCantileverVisualization()`，不进入以下后台时间积分链。

## 1. 总调用链

```text
main(argc, argv)                                      main.cpp
└─ QApplication + --mode dynamic 命令行解析
   └─ MainWindow(samplesPerSpan)                     mainwindow.cpp
      ├─ makeCantileverSpace()                       构造分析 patch
      ├─ BoundarySamplingCache(...)                  预计算可视化采样算子
      ├─ IgaVisualizationWidget(...)                 GUI 线程创建 OpenGL 控件
      ├─ DynamicCantileverRunner(...)                创建求解工作对象
      ├─ runner->moveToThread(solverThread)           移到后台线程
      └─ solverThread->start()
         └─ QThread::started -> runner.run()          DynamicCantileverRunner.cpp
            ├─ 创建材料、体力、边界条件和 DynamicsSettings
            └─ iga::solveDynamics(...)               GeneralizedAlphaDynamics.cpp
               ├─ assembleConsistentMass()           全时程一次
               ├─ diagnoseMassMatrix()               全时程一次
               ├─ generalizedAlphaFromSpectralRadius()
               ├─ computeInitialAcceleration()       解 a_0
               ├─ 提交 t=0 状态
               └─ for step=1..numberOfSteps          时间步循环
                  ├─ 同速度预测 d/v/a_(n+1)^(0)
                  ├─ for iteration=1..maxIterations  Newton 循环
                  │  ├─ evaluateDynamicResidual()
                  │  │  ├─ 重建 d_(n+1)
                  │  │  ├─ 形成 d_(n+alpha_f)、a_(n+alpha_m)
                  │  │  ├─ assembleInternalSystem(d_(n+alpha_f))
                  │  │  │  └─ computeTotalLagrangianElementResponse() 对所有单元装配
                  │  │  ├─ R=M*a_am+f_int-f_ext(t_af)
                  │  │  └─ M*=alpha_m*M+alpha_f*beta*dt^2*K_T
                  │  ├─ freeDofNorm() 判断收敛
                  │  ├─ applyHomogeneousDirichlet(M*,-R)
                  │  ├─ SimplicialLDLT 求解 M* Delta_a=-R
                  │  ├─ 同步校正 a、d、v
                  │  └─ hasPositiveDeformationJacobians(d_af)
                  ├─ 提交收敛的 DynamicState(t_(n+1))
                  └─ stepCallback(state)
                     ├─ 计算可视化应力标量
                     ├─ 写入 PlaybackHistory
                     └─ emit stepFinished(step,time)
```

## 2. 程序入口和线程切换

`main()` 只负责 Qt 生命周期和 `--samples-per-span` 参数，然后构造
`MainWindow`。`MainWindow` 在 GUI 线程中建立参考网格和 OpenGL 控件，
但不直接调用求解器。

`DynamicCantileverRunner` 通过 `moveToThread()` 归属后台 `QThread`。
线程启动信号以队列方式调用 `DynamicCantileverRunner::run()`，所以
`solveDynamics()` 的质量装配、单元积分和 Newton 线性求解都不会阻塞 GUI。

`stepFinished` 信号只发送步号和时间；真实位移与应力通过
`PlaybackHistory` 的不可变快照传递，避免在 Qt 事件队列中复制大向量。

## 3. `solveDynamics()` 进入时间循环前

进入 `for (step...)` 之前按以下顺序执行：

1. 检查 `d_0`、`v_0` 的维度、有限性及齐次边界条件。
2. `assembleConsistentMass()` 组装常量一致质量矩阵
   `M_AB=integral(rho_0 R_A R_B)dOmega_0 tensor I_3`。
3. `diagnoseMassMatrix()` 检查总质量、对称性和正定性。
4. 由 `rho_inf` 计算 `alpha_m`、`alpha_f`、`gamma`、`beta`。
5. `computeInitialAcceleration()` 在自由自由度上求解
   `M a_0=f_ext(0)-f_int(d_0)`；不能默认令 `a_0=0`。
6. 把 `(d_0,v_0,a_0)` 作为第一个物理状态写入历史并调用回调。

## 4. 一个时间步内部

### 4.1 预测

Newton 未知量选为 `a_(n+1)`：

```text
a_(n+1)^(0) = (gamma-1)/gamma * a_n
d_(n+1)^(0) = d_n + dt*v_n
                + dt^2*((1/2-beta)*a_n + beta*a_(n+1)^(0))
v_(n+1)^(0) = v_n + dt*((1-gamma)*a_n + gamma*a_(n+1)^(0))
```

该预测使 `v_(n+1)^(0)=v_n`。

### 4.2 Newton 残量和有效切线

每次迭代先形成中间时刻量：

```text
d_af = d_n + alpha_f*(d_(n+1)-d_n)
a_am = a_n + alpha_m*(a_(n+1)-a_n)
t_af = t_n + alpha_f*dt
```

然后调用静力层的 `assembleInternalSystem(d_af)` 得到 `K_T(d_af)` 和
`f_int(d_af)`。动力学层不重写单元公式，只在静力残量外增加惯性项：

```text
R  = M*a_am + f_int(d_af) - f_ext(t_af)
M* = alpha_m*M + alpha_f*beta*dt^2*K_T(d_af)
```

### 4.3 校正与提交

在施加齐次约束后求解：

```text
M* Delta_a = -R
```

并同步校正：

```text
a_(n+1) += Delta_a
d_(n+1) += beta*dt^2*Delta_a
v_(n+1) += gamma*dt*Delta_a
```

这里的后两行不能乘 `alpha_f`。只有残量求值点从 `d_(n+1)` 插值到
`d_(n+alpha_f)` 时才出现 `alpha_f`。

当自由自由度残量满足相对或绝对阈值后，当前 `n+1` 状态才写入
`DynamicsResult::states` 并传给可视化回调。Newton 中间状态和
`d_(n+alpha_f)` 都不是物理帧。

## 5. 主要文件职责

| 文件 | 职责 |
|---|---|
| `ElementCommon.*` | TL/UL 共用的运动学、Voigt 工具和单元返回类型 |
| `HyperelasticMaterial.*` | 参考/空间材料响应，以及默认的应力与切线推前 |
| `TotalLagrangianElement.*` | 在参考构形 Gauss 积分 `K_e` 和 `f_int,e` |
| `UpdatedLagrangianElement.*` | 在当前构形 Gauss 积分相同的 `K_e` 和 `f_int,e` |
| `GlobalNonlinearSystem.*` | 全局装配、边界条件、自由 DOF 范数及静力 Newton |
| `GeneralizedAlphaDynamics.*` | 一致质量、初始加速度、动态残量、时间步和 Newton 循环 |
| `DynamicCantileverRunner.*` | 设置示例参数，在后台线程调用动力学入口并发布收敛状态 |
| `mainwindow.*` | 创建线程、连接信号、读取历史并驱动可视化播放 |

`DynamicsSettings::elementFormulation` 选择 TL 或 UL，默认 TL；该选择同时传入初始
加速度和每一步的动态残量。全局散射与时间积分算法共用，详见 `UPDATED_LAGRANGIAN.md`。
