# 文献检索与研究定位（2026-08-18）

## 1. 直接 morphism Hall realization

- Ding–Xu–Zhang, *Acyclic quantum cluster algebras via Hall algebras of morphisms*, Math. Z. 296 (2020), 945–968：对有限无环箭图路径代数，从投射态射范畴 Hall algebra 构造 principal-coefficient quantum cluster algebra 的 subquotient。
- Fu–Peng–Zhang, *Quantum cluster characters of Hall algebras revisited*, Selecta Math. 29 (2023), Paper 4：在 finite acyclic valued quiver 情形构造 bialgebra/integration map，并恢复 DXZ 满同态与 quantum Caldero–Chapoton formula。
- Berenstein–Rupel 以及 Chen–Ding–Xu 的 Hall-to-quantum-cluster homomorphisms 同样以 hereditary/acyclic 输入为主。

这些结果不能直接应用于 nonhereditary cluster-tilted algebra，因为普通 Euler twist、derived reduction 与 extension middle terms 的行为不同。

## 2. Finite-type / 2-CY multiplication

- Caldero–Keller 证明 finite-type cluster algebra 可实现为 cluster category 的 exceptional Hall algebra。
- Palu 给出一般 Hom-finite 2-CY category 的 cluster-character multiplication formula。
- Chen–Xiao–Xu 给出 weighted quantum cluster functions multiplication formulas。
- Xiao–Xu–Yang 给出 motivic refinements，并明确指出一般 2-CY category 中 quantum cluster character 的显式 Laurent formula 仍未知。

这些工作解释两方向 Ext/triangle varieties，但不是从 nonhereditary `Mor(proj A)` 的 exact Hall algebra出发的显式 quotient。

## 3. Surface extensions 与 smoothing

- Canakci–Schroll 把 unpunctured surface gentle Jacobian algebra 的 indecomposable module extensions 与 generalized arcs crossings 对应，并计算 cluster category 中 smoothing middle terms。
- 该结果为正文的 crossing-leading Hall triangularity 提供核心几何输入。

## 4. Skein = quantum cluster

- Muller 证明 marked surface skein algebra 的 multicurve basis，并在适当边界局部化后得到
  `A_q ⊂ Sk_q^o ⊂ U_q`，在每个分支至少有两个 marked points 时三者相等。多边形满足条件。
- Ishibashi–Kano–Yuasa (IMRN 2025) 用 walled surfaces 实现带任意几何系数与量子化的 surface cluster algebras。
- Cao–Huang–Wang, arXiv:2605.12114 (2026) 对 polygon stated SL_n skein algebra证明 frozen-localized skein = quantum cluster = quantum upper cluster，并给出 theta-basis/rotation-invariant basis 结果。
- Cooper–Samuelson 在 partially wrapped Fukaya surface derived Hall algebra 中证明 graded skein relation，显示 Hall 与 skein 的独立拓扑联系。

## 5. 本稿相对于文献的精确增量

本稿组合并新增：

1. source 是 cluster-tilted type A 的 **exact projective-morphism Hall algebra**；
2. 经过前稿的 Hall-compatible cocycle、contractible localization 与 central coefficient specialization；
3. 对所有 ordered crossing arcs 定义显式 Hall–skein defect ideal；
4. 利用 exact Hall triangularity + extension smoothing + filtered Diamond lemma证明 quotient basis 是 non-crossing direct sums；
5. 由 polygon skein theorem 得到整个 Hall subquotient 到整个 frozen-localized quantum cluster algebra的同构；
6. 每个 non-rigid Hall middle term得到唯一 Hall-corrected theta expansion。

截至本轮定向检索，未发现把上述六部分作为同一主定理的相同论文。这不是绝对优先权声明。

## 6. 关键限制

该证明目前严格依赖 polygon/type A：

- indecomposables = arcs；
- non-split extensions = smoothing；
- crossing number提供终止良序；
- skein local relations具有 non-crossing basis；
- finite-type theta basis全部由 cluster monomials覆盖。

一般 Jacobi-finite cluster-tilted algebra尚缺这些统一输入。
