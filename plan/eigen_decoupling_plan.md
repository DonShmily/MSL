# MSL Eigen 解耦计划

## 一、Eigen 依赖全景

### 1.1 直接依赖 Eigen 的库文件（5个）

| 文件 | 使用的 Eigen 功能 | 依赖强度 |
|---|---|---|
| `msl/matrix/eigen_interface.hpp` | `Eigen::Map`, `Eigen::MatrixXd`, `Eigen::VectorXd`, `PartialPivLU`(solve) | **完全依赖**（此文件本身就是 Eigen 桥） |
| `msl/matrix/martrix_decompose.hpp` | `JacobiSVD`, `EigenSolver`, `GeneralizedEigenSolver`, `PartialPivLU`, `HouseholderQR` | **核心依赖** |
| `msl/matrix/matrix_operation.hpp` | `MatrixXd::inverse()`, `::adjoint()`, `::determinant()` | **核心依赖** |
| `msl/signal/fft.hpp` | `Eigen::FFT<double>`（12处实例化） | **核心依赖** |
| `msl/signal/filtfilt.hpp` | `MatrixXd::inverse()` + 矩阵-向量乘法（求解初始条件） | 核心依赖 |

### 1.2 传递依赖链

```
eigen_interface.hpp ──────────────────────────────────────────────────┐
  ├── martrix_decompose.hpp (SVD/eig/LU/QR)                          │
  │     ├── polynomial/polynomial.hpp  [QR → 最小二乘拟合]             │
  │     │     └── signal/detrend.hpp                                  │
  │     └── matrix_operation.hpp (inverse/adjoint/determinant)        │
  │                                                                   │
  ├── signal/fft.hpp (Eigen::FFT)                                    │
  │     ├── signal/fourier_domain_filter.hpp                          │
  │     └── signal/power_spectral_density.hpp                         │
  │                                                                   │
  └── signal/filtfilt.hpp (矩阵求逆 → 初始条件计算)                    │
        └── signal/butterworth_filter.hpp                             │
```

**结论**：Eigen 禁用后影响的模块 — 多项式拟合、全部矩阵分解、矩阵求逆/伴随/行列式、FFT及其下游（频域滤波、PSD）、零相位滤波。

### 1.3 不受影响的模块

`interp/`, `integral/`, `difference/`, `signal/window.hpp`, `signal/filter.hpp`（因果滤波），`signal/filter_design.hpp`，`polynomial/polynomial.hpp` 的纯系数构造路径。

---

## 二、替代算法复杂度评估

### 2.1 低复杂度（可直接手写，~200行以内）

| 算法 | 文件 | 替代方案 | 估计行数 |
|---|---|---|---|
| **LU 分解** | `martrix_decompose.hpp` | 列主元高斯消元（Doolittle算法） | ~150行 |
| **矩阵求逆** | `matrix_operation.hpp` | 基于 LU 分解求解 AX=I | ~80行 |
| **行列式** | `matrix_operation.hpp` | LU 分解对角元连乘 + 置换奇偶性 | ~30行 |
| **伴随矩阵** | `matrix_operation.hpp` | `adj(A) = det(A) · inv(A)`，已用行列式+逆实现 | ~20行 |
| **filtfilt初始条件** | `filtfilt.hpp` | 用 `lu()` + 前代/回代 替代 `Eigen::inverse()`；矩阵规模=滤波器阶数（通常 ≤ 10），性能无影响 | ~50行调整 |

### 2.2 中等复杂度（~300-600行）

| 算法 | 文件 | 替代方案 | 估计行数 |
|---|---|---|---|
| **QR 分解** | `martrix_decompose.hpp` | Householder 反射（标准算法，Numerical Recipes Chap 11） | ~300行 |
| **FFT** | `fft.hpp` | 基-2 Cooley-Tukey，仅支持 2^k 长度；非 2^k 可补零 | ~400行 |

> **QR 是关键路径**：`polynomial.hpp` 的多项式拟合直接依赖 `matrix::qr()`。如果 QR 不实现，多项式拟合将不可用。

### 2.3 高复杂度（>800行，建议仅提供基础版本或暂不实现）

| 算法 | 文件 | 替代方案 | 估计行数 |
|---|---|---|---|
| **SVD** | `martrix_decompose.hpp` | Golub-Reinsch 或 Jacobi SVD | ~800行 |
| **特征值分解** | `martrix_decompose.hpp` | QR 迭代算法 | ~600行 |
| **广义特征值** | `martrix_decompose.hpp` | QZ 算法 | ~1000行 |

> SVD 和特征值分解目前**无内部使用者**（仅测试用例调用）。可暂时降级为"需 Eigen 开启"功能，后续迭代实现。

---

## 三、方案设计

### 3.1 宏定义方案

采用单一宏 `MSL_USE_EIGEN` 控制：

- **未定义**（默认）：使用纯 C++20 原生实现
- **定义 `MSL_USE_EIGEN`**（用户选择）：启用 Eigen 加速后端

```cpp
// 用户使用方式：
// 默认无 Eigen：
#include "msl/matrix.hpp"       // 使用原生 LU/QR/inverse

// 开启 Eigen：
#define MSL_USE_EIGEN
#include "msl/matrix.hpp"       // 使用 Eigen 加速
```

在 `xmake.lua` 中可通过 `add_defines("MSL_USE_EIGEN")` 开启。

### 3.2 架构模式

每个需要条件编译的文件采用以下模式：

```cpp
#ifdef MSL_USE_EIGEN
    // === Eigen 实现（现有代码） ===
    #include <eigen3/Eigen/Core>
    // ... 现有 Eigen 代码 ...
#else
    // === 原生 C++20 实现 ===
    // ... 纯 C++ 实现，使用 msl::matrix 自身类型 ...
#endif
```

**关键原则**：
1. 公共 API 签名**不变**（函数名、参数类型、返回值类型完全一致）
2. 条件编译在**函数体内部**或**文件级 #ifdef 块**处理
3. `eigen_interface.hpp` 的 `as_eigen()`/`from_eigen()` 在无 Eigen 模式下不提供（用户代码中如果有调用则编译错误，这是预期行为）

### 3.3 文件改动清单

| 文件 | 改动类型 | 改动说明 |
|---|---|---|
| `msl/matrix/eigen_interface.hpp` | **重度修改** | 整个文件用 `#ifdef MSL_USE_EIGEN` 包裹；无Eigen时文件为空（仅含namespace声明） |
| `msl/matrix/martrix_decompose.hpp` | **重度修改** | LU/QR 提供原生实现；SVD/eig 在无Eigen时标记为不可用或抛出编译期错误 |
| `msl/matrix/matrix_operation.hpp` | **中度修改** | inverse/determinant/adjoint 提供基于自身LU的原生实现；transpose/trace 无影响 |
| `msl/signal/fft.hpp` | **中度修改** | 提供原生基-2 FFT 实现 |
| `msl/signal/filtfilt.hpp` | **轻度修改** | 用 `msl::matrix::lu()` 替代 `Eigen::inverse()` |
| `msl/polynomial/polynomial.hpp` | **无需修改** | 已通过 `matrix::qr()` 间接使用，只要 QR 有原生实现即可 |
| `msl/signal/fourier_domain_filter.hpp` | **无需修改** | 已通过 `signal::fft()` 间接使用 |
| `msl/signal/power_spectral_density.hpp` | **无需修改** | 通过 `signal::fft()` 间接使用 |
| `xmake.lua` | **轻度修改** | Eigen include 路径包裹在 `if has_config("use_eigen")` 中 |

### 3.4 原生实现质量要求

所有原生实现需满足：
- **数值精度**：与 Eigen 输出误差 < 1e-10（通过交叉验证测试保证）
- **列主序兼容**：与 MSL 现有矩阵存储格式一致
- **原地操作支持**：尽量不引入额外内存分配
- **异常安全**：输入检查与现有 Eigen 版本行为一致

---

## 四、实施计划（优先级排序）

### Phase 1: 基础设施（低风险，高收益）

| 序号 | 任务 | 文件 | 估算工作量 |
|---|---|---|---|
| 1.1 | 实现原生 LU 分解（列主元高斯消元） | `martrix_decompose.hpp` | 2h |
| 1.2 | 基于 LU 实现 inverse/determinant/adjoint | `matrix_operation.hpp` | 1h |
| 1.3 | 用 `msl::matrix::lu()` 替代 `filtfilt.hpp` 中的 Eigen 直接调用 | `filtfilt.hpp` | 0.5h |
| 1.4 | 添加 `#ifdef MSL_USE_EIGEN` 条件编译框架 | 上述3个文件 + `eigen_interface.hpp` | 1h |
| 1.5 | 更新 `xmake.lua` 支持 optional Eigen | `xmake.lua` | 0.5h |

**里程碑 A**：无 Eigen 环境下，LU、求逆、行列式、伴随矩阵、filtfilt 均可用。

### Phase 2: QR 分解（关键路径）

| 序号 | 任务 | 文件 | 估算工作量 |
|---|---|---|---|
| 2.1 | 实现原生 QR 分解（Householder 反射） | `martrix_decompose.hpp` | 3h |
| 2.2 | 添加条件编译 | `martrix_decompose.hpp` | 0.5h |

**里程碑 B**：无 Eigen 环境下，多项式最小二乘拟合完全可用（`polynomial.hpp` 通过 `matrix::qr()` 路径）。

### Phase 3: FFT（可选但影响面大）

| 序号 | 任务 | 文件 | 估算工作量 |
|---|---|---|---|
| 3.1 | 实现原生基-2 FFT/IFFT（实数+复数） | `fft.hpp` | 4h |
| 3.2 | 实现 2D FFT（行列分离，调用1D FFT） | `fft.hpp` | 2h |
| 3.3 | 添加条件编译 | `fft.hpp` | 0.5h |

**里程碑 C**：无 Eigen 环境下，FFT、频域滤波、PSD 均可用。

### Phase 4: 高级分解（低优先级）

| 序号 | 任务 | 文件 | 估算工作量 |
|---|---|---|---|
| 4.1 | 实现原生 Jacobi SVD | `martrix_decompose.hpp` | 6h |
| 4.2 | 实现原生 QR 迭代特征值 | `martrix_decompose.hpp` | 6h |

**里程碑 D**：完全去除 Eigen 依赖（可选），SVD 和特征值分解可用原生实现。

### Phase 5: 验证与文档

| 序号 | 任务 |
|---|---|
| 5.1 | 编写交叉验证脚本：Eigen vs 原生 输出对比 |
| 5.2 | 更新所有测试，确保两套实现均通过 |
| 5.3 | 在 CI 中增加无 Eigen 构建配置 |

---

## 五、性能评估

### 5.1 预期性能差距

| 算法 | Eigen (开启) | 原生实现 (关闭) | 预期差距 |
|---|---|---|---|
| LU 分解 | 高度优化的 Block LU | Doolittle 逐列消元 | **3-10x**（中小矩阵差距小） |
| QR 分解 | Block Householder | 逐列 Householder | **3-8x** |
| 逆矩阵 | 基于 LU + BLAS | 基于自家 LU + 前代回代 | **2-5x** |
| 行列式 | LU 对角积 | 同上 | **无差距** |
| FFT | Eigen FFT (FFTW后端) | 基-2 Cooley-Tukey | **2-5x**（仅 2^k 长度，需补零） |
| SVD | 高度优化的分治/BiDiag | Jacobi 迭代 | **10-50x** |

### 5.2 性能瓶颈分析

MSL 的典型使用场景（信号处理、数据拟合）中：
- 矩阵规模通常 ≤ 1000×1000
- LU/QR 在此规模下，原生实现与 Eigen 的绝对耗时差距在毫秒级，通常**不是瓶颈**
- FFT 长度通常 ≤ 65536，原生实现差距也在可接受范围

**结论**：对于 MSL 的典型工作负载，原生实现性能完全可接受。仅在需要大规模矩阵运算时，建议开启 `MSL_USE_EIGEN`。

---

## 六、方案评估与建议

### 6.1 宏方案优缺点

**优点**：
- 用户可自主选择，零成本切换
- 默认无外部依赖，降低构建复杂度
- 保留 Eigen 作为可选加速后端

**缺点**：
- 每个文件需维护两套实现，增加代码量
- 两套实现需保持 API 和行为一致，增加测试负担
- `eigen_interface.hpp` 的 `as_eigen()`/`from_eigen()` 用户代码在无Eigen时无法编译

### 6.2 改进建议

除了原始的单一宏方案，可考虑增加**细粒度控制宏**供高级用户使用：

| 宏 | 作用 | 默认值 |
|---|---|---|
| `MSL_USE_EIGEN` | 主开关，启用所有 Eigen 加速 | 未定义 |
| `MSL_USE_EIGEN_FFT` | 仅 FFT 使用 Eigen（需 `MSL_USE_EIGEN` 已定义才有效） | 跟随主开关 |
| `MSL_USE_EIGEN_DECOMPOSE` | 仅分解使用 Eigen | 跟随主开关 |

细粒度宏的实现成本低（在 `#ifdef MSL_USE_EIGEN` 内部再加一层判断），建议在 Phase 1 就预留接口。

### 6.3 最终建议

1. **采用单一宏 `MSL_USE_EIGEN` 方案**，简洁且满足 90% 场景
2. **Phase 1 (LU + inverse + determinant) 最优先**，这是其他所有功能的基础
3. **Phase 2 (QR) 必须完成**，否则 `polynomial` 模块不可用
4. **Phase 3 (FFT) 强烈建议**，覆盖面广（频域滤波、PSD）
5. **Phase 4 (SVD/eig) 可延迟**，无内部依赖，当前也无用户强需求
