# 文献检索与研究定位（截至 2026-08-18）

## 1. 态射范畴与 tau-tilting

用户提供的李伟宇论文证明：投射模态射范畴中的 maximal rigid objects 与 support `tau`-tilting pairs 一一对应；rigid 对象可补全；交换图同构。本文把这些结果作为范畴输入。

2026 年 Schroll–Tattar–Treffinger–Williams 在 *Mathematische Zeitschrift* 发表的工作从 wall-and-chamber 结构定义 `tau`-cluster morphism category。它支持以墙和腔室组织 rigid 对象的几何视角，但没有讨论本文的 Hall term-lattice cocycle、transfer 不变 skew 格或 fixed-frozen 障碍。

## 2. 量子化空间

Gellert–Lampe 研究 cluster algebra 的 quantisation spaces，证明满秩初始交换矩阵存在量子化并给出全部量子化的生成结构。

本文的 `Q_global` 不是替代该理论，而是回答另一个比较问题：

> 一个 principal-compatible skew 参数 `S` 何时能通过全部 presentation-to-cluster fan charts 拉回为同一个 exact Hall term-lattice 形式？

答案是 transfer group 的不变格：

```text
Q_global = Skew_n(Z)^{Hol_T}.
```

## 3. tropical duality 与 fan transport

Nakanishi–Zelevinsky 建立 cluster `g`-vectors 与 `c`-vectors 的 tropical duality。本文的 transfer

```text
L_R = G_R Delta_R^{-1}
```

比较的是两个对象标号的格：split presentation index 与 cluster `g`-vector。满秩时它保持 `Lambda_0`，并在相邻墙上成为秩一辛 transvection。

本轮定向检索没有找到把这一 piecewise transfer 明确解释为“presentation-fan symplectic holonomy”并用于分类 Hall-globalizable quantisations 的现有定理。

## 4. quantum theta basis

Davison–Mandel 构造 quantum theta bases，证明它们扩张量子簇单项式并具有正的结构常数。在 cluster chamber 中，theta label 直接给出 normalized quantum cluster monomial。

本文使用这一已有结果，把 rigid morphism 的量子对象写成

```text
theta_{Gamma_T(I(f))}.
```

新内容不在 theta basis 的构造，而在：

- 从 morphism homology 和 presentation fan 计算 label；
- 构造与 cluster chamber 乘法精确相符的 Hall cocycle；
- 分类哪些 principal gauges 可由一个全局 Hall term form 实现。

## 5. Hall 实现的适用范围

Ding–Xu–Zhang 和 Fu–Peng–Zhang 的直接 morphism-Hall 实现从有限无环箭图的遗传路径代数出发。

Zhang 的固定长度投射复形 Hall 理论给出正合 integration 的一般基础，因此 `Mor(proj A)` 的标准 exact integration 在非遗传有限维代数上仍存在；问题是它对微分盲。

Chen–Hu 2026 年关于 two-term complexes 的 quantum Weyl relations 仍假设 Dynkin quiver 和 `A=kQ`。Contu 2026 年的 semi-derived Hall quantum cluster structure 也以 representation-finite path algebra 为背景。

本文针对非遗传 cluster-tilted algebra 的结论是 rigid chamber 上的 cocycle-corrected Hall–theta atlas，不是整个 Hall basis 上的 global quantum cluster algebra structure。

## 6. 3-Calabi–Yau / critical Hall 比较

Keller–Yang、Bridgeland、Davison 等已经提供 quiver-with-potential mutation、Hall scattering 和 cohomological DT/量子正性的广泛框架。因此“是否存在某种 3-CY Hall wall-crossing”并不是准确的开放问题。

真正缺少的是：

> 对一个具体 projective morphism `f:P^- -> P^0`，规范构造 framed critical Hall/DT class，并把它与 exact morphism object、其同调和 theta label 直接比较。

本文的 recognition theorem 证明：若该类的积分具有正确初始值并满足一步 mutation covariance，则它在全部可达 rigid 对象上必等于 `theta_{Gamma_T(I(f))}`。这缩小了验证任务，但没有构造该类。

## 7. 本轮研究定位

当前检索没有找到将以下四点合为同一主定理的完全相同工作：

1. `L_R^T S L_R=S` 的 Hall-globalization 充要判据；
2. 相邻墙 `T=I+v ell^T` 的整数辛 transvection 判据；
3. `(Delta_R^T+G_R^T)(dI+SB_0)=0` 的 fixed-frozen 主系数障碍；
4. rigid face colimit 的 `ker(pi_H) subset ker(pi_theta)` 延拓判据。

这只是限定关键词、出版社页面、arXiv 和引用链范围内的研究定位，不是绝对优先权证明。正式投稿前应继续检索 MathSciNet、zbMATH、作者主页和后续引用。
