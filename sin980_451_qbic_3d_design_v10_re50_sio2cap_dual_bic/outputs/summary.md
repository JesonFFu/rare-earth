# v10: 50 nm 稀土膜 + 12.5 nm SiO2 cap 的 980/451 nm 双高 Q BIC 搜索总结

本目录为独立 v10 工作目录：`E:\xitugrating\sin980_451_qbic_3d_design_v10_re50_sio2cap_dual_bic`。没有覆盖 v1-v9、原始 `.mph`、`bic0413.mph` 或材料文件。

## 固定层结构

- 有效 Si / 1500 nm SiO2 / 50 nm 稀土纳米晶自组装膜 / 12.5 nm SiO2 保护层 / 300 nm 图形化 SiN / air。
- 只在 300 nm SiN 层内刻蚀或图形化；50 nm 稀土膜和 12.5 nm SiO2 cap 保持完整连续。
- 980 nm: `n_sin=2.069946018`, `n_re=1.355`, `n_sio2=1.45`。
- 451 nm: `n_sin=2.128247126`, `n_re=1.365`, `n_sio2=1.466`。
- Q 按复本征频率计算：`Q = Re(f)/(2*abs(Im(f)))`。

## 计算策略与耗时

- 先尝试 3D 多结构族粗扫；3D 模型在本机耗时过高，3 小时只稳定落盘了 1 个候选。
- 因此切换为 2D 快速筛选，使用同一 v10 层堆栈，扫描 1D primitive、451/980 hybrid supercell、dimerized grating、trimer supercell 等 75 个粗扫点，并对最优 980 分支做 70 个局部细扫点。
- 2D 筛选用于判断双高 Q 存在性；最终仍建议用 3D 复核任何后续真正满足双阈值的候选。

## 3D 粗扫已得到的点

3D 只完成了一个可用粗扫点：

| model family | geometry | lambda980 | Q980 | lambda451 | Q451 |
|---|---|---:|---:|---:|---:|
| single air-hole slab | P=540 nm, R=74 nm | 978.995452 nm | 2.23376e10 | 451.638257 nm | 1771.10 |

这个 3D 点证明新 stack 下 980 nm 高 Q 可以保留，但 451 nm 仍没有达到 `Q451>=1e5`。

## 2D 最佳折中结构

推荐保存为当前 v10 的 980 优先折中模型：

- 结构：1D 周期 SiN 光栅，`P=562 nm`, `fill=0.72`
- SiN ridge 宽度：`404.64 nm`
- etched gap：`157.36 nm`
- 最小主结构尺寸：`157.36 nm`
- 980 nm: `lambda=980.115156 nm`, `Q=3.54213e9`
- 451 nm: `lambda=453.291596 nm`, `Q=1300.06`

该结构满足 980 nm 高 Q BIC-like 目标，但不满足 451 nm 高 Q 目标。

## BIC 证据

980 nm 分支的 Gamma 附近 Q 迅速下降，符合 BIC-like 行为：

| kx | lambda980 | Q980 |
|---:|---:|---:|
| 0 | 980.115156 nm | 3.54213e9 |
| 0.002 pi/L | 980.119995 nm | 1.36112e6 |
| 0.006 pi/L | 980.158697 nm | 1.51462e5 |
| 0.020 pi/L | 980.597852 nm | 1.38112e4 |

451 nm 分支没有高 Q BIC 行为：

| kx | lambda451 | Q451 |
|---:|---:|---:|
| 0 | 453.291596 nm | 1300.06 |
| 0.002 pi/L | 453.290919 nm | 1295.59 |
| 0.006 pi/L | 453.285502 nm | 1260.93 |
| 0.020 pi/L | 453.223872 nm | 966.922 |

## 451 nm 高 Q 对照

为了确认 451 nm 高 Q 支路本身不是不存在，另存了一个短周期对照：

- 结构：1D 周期 SiN 光栅，`P=294 nm`, `fill=0.42`
- etched gap：`170.52 nm`
- 451 nm: `lambda=445.026698 nm`, `Q=3.35753e9`
- 同一结构 980 nm: `lambda=990.135049 nm`, `Q=7.055`

这说明 451 nm 高 Q 可以在短周期支路中出现，但与 980 nm 高 Q 支路在同一简单 SiN 单层结构中没有自然重合。

## 结论

在当前约束下，本轮没有找到同时满足 `Q980>=1e8` 与 `Q451>=1e5` 的双高 Q 结构。最好的物理解释是：单层 300 nm SiN、连续低折射率稀土膜和 cap 的体系中，980 nm BIC 分支和 451 nm BIC 分支更像两个不同尺度的导模/折叠模；简单 1D/超胞调制可以把其中一支调高 Q，但会让另一支落到低 Q 辐射支路。

下一步如果继续追求真正双高 Q，建议优先增加结构自由度，而不是在单周期附近继续细扫：

- 允许 SiN 内做两个可制造刻蚀深度，例如主光栅 + 宽浅槽，但槽宽/肋宽保持 `>=100 nm`。
- 用 2D PhC slab 的双孔/多孔超胞继续做 3D 复核，但需要更长计算时间或更粗的初筛代理。
- 考虑 2x2 或 3x3 supercell 的多孔/椭圆孔结构，让 451 支路和 980 支路分别由不同折叠模承担。

## 输出文件

- `overall_best_dual_bic.mph`: 当前 980 优先折中结构的 2D 可执行 COMSOL 模型，材料设置为 980 nm。
- `overall_best_dual_bic_451check.mph`: 同一几何的 451 nm 材料设置模型。
- `fast2d/fast2d_tradeoff_980.mph`, `fast2d/fast2d_tradeoff_451.mph`: 折中结构的两套本征模模型。
- `fast2d/fast2d_high451_control_980.mph`, `fast2d/fast2d_high451_control_451.mph`: 451 高 Q 对照模型。
- `dual_best_candidates.csv`: 已保存模型的关键结果。
- `band_980.csv`, `band_451.csv`, `q_near_gamma_980.csv`, `q_near_gamma_451.csv`: 能带与 Gamma 附近 Q 数据。
- `band_980.svg`, `band_451.svg`: 简化能带/Q 图。
- `fast2d/coarse_scan_2d.csv`, `fast2d/fine_scan_2d.csv`: 2D 粗扫和细扫完整数据。
- `coarse_dual_scan.csv`: 3D 粗扫已完成数据。
