# MSL FFT 模块使用指南

## 概述

MSL FFT 模块基于 Eigen 的 unsupported FFT 功能，提供了类型安全、易用的 FFT 接口。支持：
- 1D FFT（实数/复数输入）
- 2D FFT（按行/列/全2D）（TODO：全2D）
- 逆变换
- 频谱分析工具
- 窗函数

## 依赖

```cpp
#include "fft.hpp"
#include <unsupported/Eigen/FFT>  // Eigen FFT (自动包含)
```

**注意**：需要实现 `complex_matrix_owned.hpp` (matrixc)

---

## 1. 基本 1D FFT

### 实数输入 → 复数输出

```cpp
// 生成信号
std::vector<double> signal(128);
for (size_t i = 0; i < 128; ++i) {
    signal[i] = std::sin(2.0 * M_PI * 5.0 * i / 100.0);
}

// 计算 FFT
auto fft_result = msl::signal::fft(signal);

// fft_result.size() = 65 (N/2 + 1，利用共轭对称性)
```

**关键点**：
- 输入 N 个实数样本
- 输出 N/2+1 个复数值（负频率通过共轭对称性隐含）
- 节省内存和计算

### 复数输入 → 复数输出

```cpp
std::vector<std::complex<double>> complex_signal(64);
// ... 初始化 ...

auto fft_result = msl::signal::fft(complex_signal);
// fft_result.size() = 64 (完整的复数 FFT)
```

### 使用 span（零拷贝）

```cpp
double data[256];
// ... 填充数据 ...

// 使用 span 避免拷贝
auto result = msl::signal::fft(std::span<const double>(data, 256));
```

---

## 2. 逆 FFT (IFFT)

### 复数 → 复数

```cpp
auto fft_result = msl::signal::fft(signal);
auto reconstructed = msl::signal::ifft(fft_result);

// reconstructed 应该接近原始信号
```

### 复数 → 实数（用于实数 FFT 的逆变换）

```cpp
std::vector<double> signal(128);
// ... 初始化 ...

auto fft_result = msl::signal::fft(signal);  // 128 → 65 复数值

// 逆变换回实数
auto reconstructed = msl::signal::ifft_real(fft_result, 128);
// 必须指定原始长度！
```

**注意**：必须提供原始信号长度，因为 FFT 结果不包含此信息。

---

## 3. 2D FFT（矩阵）

### 按列 FFT（推荐用于列主序矩阵）

```cpp
// 创建矩阵，每列是一个信号
matrixd signals(256, 10);  // 10个信号，每个256点
// ... 填充数据 ...

// 对每一列独立做 FFT
auto fft_result = msl::signal::fft_columns(signals);
// 结果: (129 x 10) 复数矩阵

// 逆变换
auto reconstructed = msl::signal::ifft_columns_real(fft_result, 256);
```

**用途**：批量处理多个信号（如音频多声道、传感器阵列）

### 按行 FFT

```cpp
matrixd signals(10, 256);  // 10个信号，每个256点（按行排列）

auto fft_result = msl::signal::fft_rows(signals);
// 结果: (10 x 129) 复数矩阵

auto reconstructed = msl::signal::ifft_rows_real(fft_result, 256);
```

**注意**：行访问在列主序矩阵中效率较低，优先使用 `fft_columns`。

### 完整 2D FFT（图像处理）

```cpp
// 创建"图像"（2D 数据）
matrixd image(512, 512);
// ... 填充图像数据 ...

// 完整 2D FFT（在两个维度上都做变换）
auto fft_2d = msl::signal::fft2(image);

// 逆变换
auto reconstructed = msl::signal::ifft2_real(fft_2d, 512, 512);
```

**用途**：图像处理、2D 信号分析

---

## 4. 频谱分析工具

### 计算频率轴

```cpp
size_t n = 256;
double sampling_rate = 1000.0;  // Hz

auto freqs = msl::signal::fft_frequencies(n, sampling_rate);
// freqs[i] = i * (sampling_rate / n)
// 从 0 到 Nyquist 频率
```

### 功率谱

```cpp
auto fft_result = msl::signal::fft(signal);

// 功率谱 = |FFT|²
auto power = msl::signal::power_spectrum(fft_result);

// 找到最大功率的频率
size_t peak_idx = std::max_element(power.begin() + 1, power.end()) 
                  - power.begin();
double peak_freq = freqs[peak_idx];
```

### 幅度谱和相位谱

```cpp
auto magnitude = msl::signal::magnitude_spectrum(fft_result);  // |FFT|
auto phase = msl::signal::phase_spectrum(fft_result);          // arg(FFT)
```

---

## 5. 窗函数

### 为什么需要窗函数？

FFT 假设信号是周期性的。如果信号不是整数个周期，会产生**频谱泄漏**。

### Hamming 窗

```cpp
std::vector<double> signal(256);
// ... 初始化信号 ...

// 应用 Hamming 窗（修改 signal）
msl::signal::apply_hamming_window(signal);

auto fft_result = msl::signal::fft(signal);
```

### Hann 窗

```cpp
msl::signal::apply_hann_window(signal);
```

### 选择窗函数

| 窗函数       | 主瓣宽度 | 旁瓣抑制      | 用途                 |
| ------------ | -------- | ------------- | -------------------- |
| 无窗（矩形） | 窄       | 差 (-13 dB)   | 信号正好是整数周期   |
| Hamming      | 中等     | 好 (-43 dB)   | 通用，频率分辨率重要 |
| Hann         | 中等     | 较好 (-32 dB) | 通用                 |

---

## 6. 完整示例：频率检测

```cpp
#include "fft.hpp"
#include <iostream>
#include <vector>
#include <cmath>

int main() {
    // 参数
    const size_t n = 1024;
    const double fs = 8000.0;  // 8 kHz 采样率
    const double f0 = 440.0;   // A4 音符
    
    // 生成信号
    std::vector<double> signal(n);
    for (size_t i = 0; i < n; ++i) {
        double t = i / fs;
        signal[i] = std::sin(2.0 * M_PI * f0 * t);
    }
    
    // 应用窗函数
    msl::signal::apply_hamming_window(signal);
    
    // 计算 FFT
    auto fft_result = msl::signal::fft(signal);
    auto magnitude = msl::signal::magnitude_spectrum(fft_result);
    auto freqs = msl::signal::fft_frequencies(n, fs);
    
    // 找峰值
    size_t peak_idx = 0;
    double peak_val = 0.0;
    for (size_t i = 1; i < magnitude.size(); ++i) {
        if (magnitude[i] > peak_val) {
            peak_val = magnitude[i];
            peak_idx = i;
        }
    }
    
    std::cout << "检测到的频率: " << freqs[peak_idx] << " Hz\n";
    std::cout << "预期频率: " << f0 << " Hz\n";
    
    return 0;
}
```

---

## 7. 性能考虑

### FFT 大小

- **2的幂次**：最快（如 128, 256, 512, 1024）
- **小素数因子**：较快（如 120 = 2³×3×5）
- **大素数**：慢

```cpp
// 好
fft(signal_256);   // 2⁸
fft(signal_1024);  // 2¹⁰

// 可以接受
fft(signal_240);   // 2⁴×3×5

// 避免
fft(signal_997);   // 素数
```

### 批量处理

使用矩阵版本批量处理多个信号：

```cpp
// 慢：循环调用
for (int i = 0; i < 100; ++i) {
    auto result = fft(signals[i]);
}

// 快：矩阵批处理
matrixd all_signals(256, 100);
// ... 填充 ...
auto all_results = fft_columns(all_signals);
```

### 内存使用

- `fft(real)` → 输出大小 = (N/2+1)
- `fft(complex)` → 输出大小 = N
- 2D FFT 可能需要大量内存

---

## 8. 常见陷阱

### ❌ 错误：忘记指定 IFFT 长度

```cpp
auto fft_result = fft(signal);  // 128 → 65
auto wrong = ifft_real(fft_result);  // ❌ 缺少长度参数！
```

✅ **正确**：
```cpp
auto correct = ifft_real(fft_result, 128);
```

### ❌ 错误：不使用窗函数

```cpp
// 信号不是整数周期
auto fft_result = fft(signal);  // ❌ 频谱泄漏！
```

✅ **正确**：
```cpp
apply_hamming_window(signal);
auto fft_result = fft(signal);
```

### ❌ 错误：混淆采样率和频率分辨率

```cpp
auto freqs = fft_frequencies(n, 1.0);  // ❌ 忘记提供采样率
```

✅ **正确**：
```cpp
double fs = 1000.0;  // 1 kHz
auto freqs = fft_frequencies(n, fs);
// 频率分辨率 = fs / n
```

### ❌ 错误：在列主序矩阵上用 fft_rows

```cpp
matrixd signals(256, 10);  // 列主序
auto result = fft_rows(signals);  // ❌ 效率低！
```

✅ **正确**：
```cpp
auto result = fft_columns(signals);  // 更快
```

---

## 9. 与其他库的比较

| 特性     | MSL FFT     | FFTW   | NumPy FFT |
| -------- | ----------- | ------ | --------- |
| 易用性   | ⭐⭐⭐⭐⭐       | ⭐⭐⭐    | ⭐⭐⭐⭐⭐     |
| 性能     | ⭐⭐⭐⭐        | ⭐⭐⭐⭐⭐  | ⭐⭐⭐⭐      |
| 安装     | Header-only | 需编译 | Python    |
| 类型安全 | ✓           | ✗      | ✗         |
| 矩阵支持 | ✓           | 需手动 | ✓         |

---

## 10. TODO 和未来改进

- [ ] 实现 `complex_matrix_owned` (matrixc)
- [ ] 添加更多窗函数（Blackman, Kaiser 等）
- [ ] 支持原地 FFT（节省内存）
- [ ] 添加 FFTW 后端选项（更快）
- [ ] 支持多线程批处理
- [ ] 添加 STFT（短时傅里叶变换）
- [ ] 添加频谱图生成

---

## 编译示例

```bash
# 需要 Eigen
g++ -std=c++20 -I/path/to/eigen3 example.cpp -o example

# 运行
./example
```

**注意**：Eigen 的 FFT 是 header-only，无需链接额外库。