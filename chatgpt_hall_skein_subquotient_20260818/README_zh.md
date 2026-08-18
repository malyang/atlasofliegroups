# Hall–Skein 子商研究包

## 论文题目

**Hall–Skein Subquotients for Cluster-Tilted Type A: Incompatible Hall Products, Non-Rigid Theta Straightening, and a Global Quantum Cluster Realization**

## 本轮解决的问题

此前版本已经证明：

- 完整态射同调重建 reduced 二项对象、index、刚性与量子种子；
- 经过 term-lattice bicharacter 修正，Hall 乘法在每个刚性兼容锥上等于量子簇单项式乘法；
- 但不相容 rigid 对象的 Hall 乘积包含非分裂扩张和 non-rigid middle terms，因此局部 Hall–theta 图表尚未形成整个代数上的同态。

本稿在 **cluster-tilted Dynkin type A** 中补上这一层。利用 `(n+3)` 边形模型，对每一对相交对角线加入一个显式 Hall–skein defect：

```tex
\mathcal R_{\gamma,\delta}
=\varepsilon_\gamma\star\varepsilon_\delta
-\omega^{a(\gamma,\delta)}\varepsilon_{\gamma\#_0\delta}
-\omega^{b(\gamma,\delta)}\varepsilon_{\gamma\#_\infty\delta}.
```

由这些 defect 生成双边理想 `J_sk`。主定理为

```tex
\overline{\mathcal H}^{\star}_T/J_{\rm sk}
\cong \mathscr S^{\rm fr}_\omega(P_{n+3})
\cong \mathcal A^{\rm fr}_\omega(P_{n+3}).
```

这里源代数来自 **整个** 投射态射范畴 exact Hall algebra，经：

1. Hall-compatible bicharacter twist；
2. contractible classes 局部化；
3. 中心系数特化；
4. Hall–skein 理想取商。

因此这是 Ding–Xu–Zhang / Fu–Peng–Zhang 风格的全局 subquotient realization，但当前严格范围是 cluster-tilted type A。

## non-rigid middle terms 的解释

每个 non-rigid Hall basis object 对应一个有交叉的 multidiagonal。正文证明：

- split middle term 是 Hall 乘积中唯一最高交叉项；
- 每个 non-split middle term 都由 smoothing 获得，交叉数严格下降；
- Hall-corrected reduction 终止；
- filtered Diamond lemma 给出合流性；
- normal forms 正好是非交叉 multidiagonals，即 rigid direct sums / quantum cluster monomials / chamber theta functions。

于是任意 Hall basis object `M` 有唯一展开

```tex
\Phi_{\rm sk}(\varepsilon_M)
=\sum_R c_{M,R}(\omega)\,\vartheta_{g(R)}.
```

对于任意、包括不相容的 Hall 元 `x,y`：

```tex
\Phi_{\rm sk}(x\star y)
=\Phi_{\rm sk}(x)\Phi_{\rm sk}(y).
```

## 文件

- `hall_skein_subquotient.tex`：主 TeX 文件
- `sections/*.tex`：论文分节
- `hall_skein.py`：多边形 arcs、crossings、smoothings、递归 straightening 与合流检查
- `verify_A4_hall_skein.py`：type A4 / 七边形完整证书
- `verify_typeA_hall_skein.py`：rank 2–6 回归
- `generated/*.json`：机器可读证书
- `PROOF_STATUS_zh.md`：证明边界审计
- `LITERATURE_NOTES_zh.md`：文献定位
- `reproduce.sh`：一键复现

## 复现

```bash
python -m pip install sympy
./reproduce.sh
```

脚本将重建 JSON 证书、生成 TeX 宏并编译 PDF。

## 严格边界

本稿没有声称对任意 Jacobi-finite / 任意 cluster-tilted 代数都完成全局 Hall 同态。一般情形尚缺：

- 所有 non-rigid objects 的局部 skein presentation；
- extension middle terms 严格降低的全局良序；
- 与 theta basis 匹配的合流 normal-form system；
- exact morphism Hall strata 到 critical/cohomological Hall classes 的直接比较。

代码验证的是多边形 straightening，不代替一般 Hall 结构常数的证明或计算。
