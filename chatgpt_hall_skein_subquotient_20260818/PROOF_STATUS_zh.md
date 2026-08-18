# 证明状态与边界审计

## A. 用户提供文献直接支持的输入

李伟宇《态射范畴中的极大 rigid 对象》支持：

- `Mor(proj A)` 中不可分解对象的 `C_M, Z_P, K_P` 框架；
- support τ-tilting pairs 与 maximal rigid morphism objects 的双射；
- rigid 对象补全；
- 基本 maximal rigid 对象的不可分解直和项计数；
- maximal-rigid exchange graph 与 support-τ-tilting graph 同构。

该文不包含本文的同调正规形、Hall cocycle、skein quotient 或 non-rigid theta straightening。

## B. 文献定理传递

1. **AIR / BMRRT**：support τ-tilting、two-term silting 与 type-A cluster category 多边形组合。
2. **Canakci–Schroll**：无穿孔曲面 gentle Jacobian algebra 的 extension–crossing 对应与 middle-term smoothing。
3. **Muller**：marked surface skein basis，以及边界局部化后的 skein = quantum cluster = quantum upper cluster。
4. **Ishibashi–Kano–Yuasa**：walled skein coefficient systems。
5. **Cao–Huang–Wang (2026)**：polygon stated skein 与 frozen-localized quantum cluster algebra 的等同及 theta-basis refinement。
6. **Davison–Mandel**：cluster chamber 中 theta function 等于 normalized quantum cluster monomial。
7. **Caldero–Keller / Palu / Chen–Xiao–Xu / Xiao–Xu–Yang**：2-CY cluster multiplication 背景。
8. **Ding–Xu–Zhang / Fu–Peng–Zhang**：无环 morphism Hall subquotient realization 的比较基准。

## C. 本文直接证明

### C1. 整个 reduced Hall algebra 的生成

对 representation-finite Krull–Schmidt exact category，ordered indecomposable Hall products 对 direct-sum basis 具有可逆的 degeneration-triangular transition，因此 indecomposable classes 生成整个 Hall algebra。应用于 cluster-tilted type A 后，所有 reduced indecomposables 都由多边形对角线索引。

### C2. Crossing-leading triangularity

利用 exact-to-homotopy Ext 识别与 Canakci–Schroll smoothing：

- split extension 给出 crossing direct sum；
- non-split middle term 至少 smooth 一个 crossing；
- 因而其交叉数严格下降。

### C3. 显式 Hall–skein ideal

对每个 ordered crossing pair 定义

```tex
\mathcal R_{\gamma,\delta}
=\varepsilon_\gamma\star\varepsilon_\delta
-\omega^a\varepsilon_{\gamma\#_0\delta}
-\omega^b\varepsilon_{\gamma\#_\infty\delta}.
```

`J_sk` 是这些元素生成的双边理想，不是事后定义的未知 kernel。

### C4. Filtered Diamond lemma

Hall–skein defect 的 leading Hall basis term 是 split crossing object，系数可逆。定向 reduction 严格降低 crossing number。对 disjoint 和 overlapping critical ambiguities：

- skein locality 处理几何部分；
- Hall associativity 使两条 Hall-corrected reduction 的差落在同一理想且最高项抵消；
- 对较低 crossing number 归纳，得到合流。

因此 quotient 的 basis 是 non-crossing multidiagonals。

### C5. 全局 subquotient theorem

由 Hall quotient 与 polygon skein algebra 具有相同 generators、relations 和 non-crossing basis，证明

```tex
\overline{\mathcal H}^{\star}_T/J_{\rm sk}
\cong \mathscr S^{\rm fr}_\omega(P_{n+3})
\cong \mathcal A^{\rm fr}_\omega(P_{n+3}).
```

从而得到整个 centrally reduced exact Hall algebra 到整个 frozen-localized quantum cluster algebra 的满代数同态。

### C6. Non-rigid theta straightening

对任意 Hall basis object，选择 crossing pair，展开真实 Hall product，解出 split leading term并递归。Diamond lemma 保证结果与选择无关，得到唯一有限 theta expansion。

## D. 计算证书

代码在 `omega=1` 的经典 Ptolemy 特化下验证：

- type A4 七边形的 14 条对角线、35 个 crossing pairs、42 个 triangulations；
- 14^3 个 ordered triple associativity checks；
- degree-four multidiagram 的全部 first-crossing confluence；
- type A2–A6 的确定性跨秩回归。

裸弦图不包含 endpoint elevation/state 或 wall data，因此代码不伪造一般量子权。量子合流性由 stated/walled skein 文献定理传递。代码也不自动计算 arbitrary bound-quiver Hall numbers。

## E. 未证明 / 未宣称

1. 任意 Jacobi-finite cluster-tilted algebra的相同 subquotient theorem。
2. punctured surfaces、band objects、essential loops 的完整 Hall–skein quotient。
3. Hall-corrected theta coefficients的正性或 bar-invariance。
4. exact morphism Hall algebra 与 critical/cohomological Hall algebra 的直接函子或同构。
5. arbitrary bound quiver 的自动 Hall-number engine。
6. 当前定向文献检索构成绝对优先权证明。

## F. 投稿前需要独立核查的关键点

- direct-sum extensions 的 simultaneous smoothing 归纳是否覆盖全部 exact middle terms；
- endpoint-state / coefficient normalization 与所选 Hall cocycle 的逐项一致；
- filtered Diamond lemma 对所有 boundary-state critical pairs 的详细局部分类；
- “整个 Hall algebra”表述中的 localization、central specialization 与 base change 的精确定义域。

这些点在正文中给出了证明链，但仍应由 Hall algebra、gentle surface algebra 和 stated skein 方向专家独立审阅。
