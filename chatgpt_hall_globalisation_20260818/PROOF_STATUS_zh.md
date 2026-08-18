# 证明状态与边界审计

## A. 本文直接证明

### A1. 任意秩 mutable Hall 全局化判据

对主系数量子化参数 `S`，证明以下条件等价：

```text
存在统一的 contractible-radical term-lattice 形式
<=> Delta_R^T S Delta_R = G_R^T S G_R 对所有 R 成立
<=> L_R^T S L_R = S 对所有 R 成立。
```

并证明全局 term 形式唯一为 `D^T S D`。

### A2. 全局化参数格

证明全部参数组成 transfer group 的整数不变格

```text
Q_global = Skew_n(Z)^{Hol_T}.
```

有限 atlas 下可由整数线性系统精确计算。

### A3. 规范主系数任意秩定理

证明 `S=0` 对任意 skew-symmetric `B_0` 都全局化；对应 cocycle 为

```text
Omega_0 = - Lambda_Mor.
```

在每个 rigid chamber 中，修正 Hall 乘法与 canonical principal quantum form 的 mutable 乘法一致。

### A4. 秩一墙与辛 transvection

证明相邻转移

```text
T = I + v ell^T,   ell^T v = 0,
```

以及

```text
T^T S T = S  <=>  S v = c ell, c in Z.
```

满秩兼容形式下，全部非零墙是整数辛 transvection。

### A5. 完整 fixed-frozen 判据

证明完整主系数形式在固定 frozen basis 上全局化，当且仅当同时满足

```text
Delta_R^T S Delta_R = G_R^T S G_R,
(Delta_R^T + G_R^T)(d I + S B_0) = 0.
```

### A6. 规范 frozen no-go 与满秩恢复

证明：

- `S=0` 时完整 fixed-frozen 全局化要求 `Delta_R=-G_R`；
- 满秩 `B_0^T Lambda_0=dI` 时，`S=Lambda_0` 使 `dI+Lambda_0B_0=0`，从而整个 principal form 全局化，frozen 方向中心化。

### A7. universal face-colimit 核判据

构造 rigid face algebra 的普适余极限 `U_T`，证明延拓所有局部 Hall–theta 同构的全局同态存在，当且仅当

```text
ker(pi_H) subset ker(pi_theta).
```

### A8. critical assignment 识别原理

证明任何满足初始归一化和一步变异协变的 integrated critical assignment，在全部可达 rigid 态射上必等于

```text
theta_{Gamma_T(I(f))}.
```

## B. 从前序统一稿保留的直接证明

- 态射同调正规形与 Krull–Schmidt 余项消去；
- 完整同调重建 index、两两 `Ext^1`、rigidity 和经典簇角色；
- presentation orbit 余维等于自扩张维数；
- rigid open orbit 唯一性和 finite-field rigid faithfulness；
- fan–orbit equivalence；
- split-index fan 到 cluster `g`-fan 的分片整幺模运输；
- exact Hall term lattice、Euler/skew form 和 differential blindness；
- full-rank cocycle-corrected Hall–theta chamber isomorphism。

本稿对这些内容作了足够的自洽回顾，但没有把已有前序稿的每个证明逐字重复。

## C. 由文献定理严格传递

- 李伟宇：态射范畴 maximal rigid 与 support `tau`-tilting 对应、补全、交换图；
- AIR：support `tau`-tilting/two-term silting 和 Bongartz completion；
- DIJ：`g`-扇、整幺模性、`tau`-tilting finite 完备性；
- Derksen–Fei：general presentations 与 generic `E`-invariant；
- Palu、Plamondon：cluster character 与 generic character；
- Nakanishi–Zelevinsky：tropical duality；
- Berenstein–Zelevinsky：量子种子变异；
- Davison–Mandel：quantum theta basis 及 cluster chamber 中的量子簇单项式识别；
- Zhang：固定长度投射复形/正合 Hall integration 的一般基础；
- Keller–Yang、Bridgeland、Davison：3-CY 变异、Hall scattering 与 cohomological DT 背景。

## D. 有限计算证书

### D1. 完整非遗传 A4

程序穷举：

- 1008 个有标号 paired seeds；
- 42 个无标号 rigid chambers；
- 4032 条有向变异；
- 所有 transfer 不变方程；
- 所有秩一墙分解与辛条件；
- 完整 fixed-frozen 仿射参数空间；
- `C_{S2}` 的显式 canonical frozen no-go。

### D2. 跨秩回归

固定确定性变异词，在满秩 2、4、6 及奇异 3、5 阶矩阵上验证：

- `S=0` 的任意秩 mutable 全局化；
- full-rank compatible form 的辛保持；
- full-rank 完整 fixed-frozen 全局化；
- canonical `S=0` 的完整 frozen 失败样本。

## E. 仍然开放

1. 证明或否定 `ker(pi_H) subset ker(pi_theta)`；
2. 构造不相容 Hall 扩张的自然过滤、商或 completion；
3. 构造具体 rigid morphism 的 canonical framed critical Hall/DT class；
4. 比较 exact morphism orbit strata 与 vanishing-cycle/cohomological Hall strata；
5. 计算更大、无限型家族的 transfer holonomy 与 globalizable skew 格；
6. 非完备 `tau`-tilting fan 补集中的有理/整数点问题；
7. 从一般 bound quiver 自动生成全部输入数据的通用软件。

## F. 不作的过度宣称

- 不声称 rigid Hall basis 对不相容 Hall 乘法闭合；
- 不声称已构造整个 exact Hall algebra 到量子簇代数的代数同态；
- 不声称存在性已由 critical recognition theorem 证明；
- 不声称有限证书替代一般证明；
- 不声称当前检索构成绝对优先权证明。
