"""
La3Ni2O7 (327) 双层两轨道(e_g) + O 2p 模型：哈密顿量构造

位点/轨道索引（0 起）:
  0  Ni0 d_{x2-y2}          1  Ni0 d_{z2}
  2  Ni1 d_{x2-y2}          3  Ni1 d_{z2}
  4,5  下层平面 O p_x       6,7  下层平面 O p_y
  8,9  上层平面 O p_x       10,11 上层平面 O p_y
  12   顶位 (apical) O p_z

自旋：s=0 上自旋, s=1 下自旋；mode = 2*orb + s

注：本模型为 e_g-only 近似（论文亦称 "limiting the basis set to the e_g
orbitals as an approximation"）。d-d 相互作用取 Racah 参数 A,B,C 在 e_g²
子空间的有效形式 (U, U', J)，对应 d⁸ 谱项：
  ³A₂(S=1) = U' - J = A - 8B（Hund 基态）
  ¹Eg(S=0) = U' + J = A + 2C
  ¹A₁g     = U + J  = A + 8B + 4C
  同轨道双占库仑 U_aa = A + 4B + 3C（= F0 + 4F2 + 36F4）
注意 ¹Eg 不是自由离子的 A - 3B + 2C：¹D 含 t2g 成分，e_g-only 基组取不到该值。
"""
import math
import os
from itertools import combinations

import numpy as np
# qiskit 只在需要导出 FermionicOp / SparsePauliOp 时才导入（见 build_hamiltonian 末尾）：
# 精确对角化走 return_terms=True，可以完全跳过这段 import，每次运行省约 3 s。

################################################################
tpdx = 1.58    # Ni 3dx2-y2 - intra O 2p
tpdz = 1.58/math.sqrt(3)    # Ni 3dz2 - intra O 2p (论文 T_z2,npd = t_pd/√3)
tpd_z = 1.66    # Ni 3dz2 - inter O 2pz (论文 T_dO)
tpp = 0.562      # O 2px - O 2py
tpp_z = 0.487    # intra O 2p - inter O 2pz (论文 T_pO)
A = 6.0       # Racah A
B = 0.15      # Racah B
C = 0.58      # Racah C
up = 4.0      # O 2p 同轨道 Hubbard U（论文 U_p = 4）
eps_dx2 = 0.0  # Ni 3d_x2-y2 在位能（参考零点）
eps_d = 0.095  # Ni 3d_z2   在位能（晶场分裂）
eps_p = 2.90   # 面内 O 2p 在位能
eps_p_z = 3.24 # 顶位 O pz 在位能

# Ni 3d e_g 相互作用：Racah 算符投影到 e_g² 子空间的严格结果，无自由参数。
#   U  = A + 4B + 3C = 8.34   （同轨道双占）
#   U' = A - 4B + C  = 5.98   （异轨道）
#   J  = 4B + C      = 1.18   （Hund 交换 / 配对跳跃）
DD_RACAH = (A + 4.0*B + 3.0*C, A - 4.0*B + C, 4.0*B + C)

################################################################
NI = [0, 1, 2, 3]
O_DN = [4, 5, 6, 7]
O_UP = [8, 9, 10, 11]
O_AP = [12]
O_ALL = O_DN + O_UP + O_AP
N_ORB = len(NI) + len(O_ALL)
N_MODES = 2 * N_ORB


def mode(orb: int, spin: int) -> int:
    return 2 * orb + spin


def add_h_function(terms, modes, f):
    """把只依赖 modes 上空穴数 h 的函数 f(h)（h = 0..len(modes)）写成密度乘积。

    n_i 是本征值 0/1 的占据数，故 C(h, k) 就是"任选 k 个模式同时被占据"的投影
    之和；Newton 展开 f(h) = Σ_k c_k C(h,k) 给出各阶密度乘积的系数。
    """
    from math import comb
    n = len(modes)
    c = [sum((-1) ** (k - j) * comb(k, j) * f[j] for j in range(k + 1))
         for k in range(n + 1)]
    for k in range(n + 1):
        if abs(c[k]) < 1e-12:
            continue
        if k == 0:
            terms[""] = terms.get("", 0.0) + c[0]
            continue
        for subset in combinations(modes, k):
            key = " ".join(f"+_{m} -_{m}" for m in subset)
            terms[key] = terms.get(key, 0.0) + c[k]
    return terms


def cut_interaction_above_d8(terms, modes, G):
    """参考的 d-d 相互作用只在 d⁸（h=2）上生效。

    空穴数 h=3/4 时，上面这套 Kanamori 给出的是与自旋排布无关的常数：
    h=3 必有一个轨道双占，能量 U + 2U' - J = G；h=4 则为 2G。
    因此减去 G·δ_{h,3} + 2G·δ_{h,4} 即可把 h≥3 的相互作用清零。
    """
    return add_h_function(terms, modes, [0.0, 0.0, 0.0, -G, -2.0 * G])


def build_hamiltonian(u_p=up,
                      tpdx=tpdx, tpdz=tpdz, tpd_z=tpd_z, tpp=tpp, tpp_z=tpp_z,
                      eps_dx2=eps_dx2, eps_d=eps_d, eps_p=eps_p, eps_p_z=eps_p_z,
                      reference_zero=None,
                      d_level: bool = True, d8_only: bool = True,
                      return_terms: bool = False):
    """构造哈密顿量（返回 FermionicOp / SparsePauliOp）。

    模型定义在**空穴语言**下：'+_i' = h†_i（产生空穴），'-_i' = h_i，
    n_i = '+_i -_i' 是该自旋轨道上的空穴数。本模型固定 5 个空穴（21 个电子），
    站点上的双占、Hund 交换、配对都由空穴承担，所以密度必须按空穴数写。
    若改用电子数 eps·c†c、U·n↑n↓，换到空穴语言后会变成 eps·(1-n)、
    U·(1-n↑)(1-n↓)，那是另一个算符，而不是差一个常数。

    与 5hole 参考程序一致的约定：
      * 在位能按空穴数计：每个空穴 +eps；
      * d-d 相互作用写成空穴密度（U、U'、U'-J、交换、配对）；
      * d 能级：每个 Ni 加 A + 1.5A·|h-2|，h 为该 Ni 的 e_g 空穴数
        （f(2)=0，f(1)=f(3)=1.5A，f(0)=f(4)=3A）；
      * 相互作用只在 d⁸（h=2）子空间上生效；h=3/4 上 Kanamori 给出的是常数
        G 和 2G（G = U + 2U' - J），按参考截掉。

    返回的项字典就是最终形式（空穴语言），不再做任何语言转换。
    """
    terms = {}

    def add(key, val):
        terms[key] = terms.get(key, 0.0) + val

    def hop(o1, o2, t):
        """-t(c†₁c₂ + h.c.)。

        换到空穴语言 c†₁ → h₁、c₂ → h†₂，即 -t(h₁h†₂ + h.c.)；i≠j 时
        h₁h†₂ = -h†₂h₁，所以每项化成 +t(h†₂h₁ + h.c.)——系数变号、脚标对调。
        相位表照旧，因为跳积分本身没动。
        """
        for s in (0, 1):
            add(f"+_{mode(o2, s)} -_{mode(o1, s)}", t)
            add(f"+_{mode(o1, s)} -_{mode(o2, s)}", t)

    def onsite(o, eps):
        """每个空穴 +eps。"""
        for s in (0, 1):
            add(f"+_{mode(o, s)} -_{mode(o, s)}", eps)

    dphase = {4: 1.0, 5: -1.0, 6: -1.0, 7: 1.0,
              8: 1.0, 9: -1.0, 10: -1.0, 11: 1.0} ## d_x2-y2 相位
    dzphase = {4: -1.0, 5: 1.0, 6: -1.0, 7: 1.0,
               8: -1.0, 9: 1.0, 10: -1.0, 11: 1.0}

    for p in O_DN:
        hop(0, p, tpdx * dphase[p])
        hop(1, p, tpdz * dzphase[p])
    for p in O_UP:
        hop(2, p, tpdx * dphase[p])
        hop(3, p, tpdz * dzphase[p])

    hop(1, O_AP[0], tpd_z)
    hop(3, O_AP[0], -tpd_z)

    # 面内 O p-p：符号随 dphase 乘积交替
    for plane in ((4, 5, 6, 7), (8, 9, 10, 11)):
        for px in plane[:2]:
            for py in plane[2:]:
                hop(px, py, -tpp * dphase[px] * dphase[py])

    # 顶位 pz - 面内 p
    for i, o in enumerate(O_DN + O_UP):
        s = 1.0 if i % 2 == 0 else -1.0
        hop(O_AP[0], o, (s if i < 4 else -s) * tpp_z)

    for d in (0, 2):               # d_x2-y2
        onsite(d, eps_dx2)
    for d in (1, 3):               # d_z2
        onsite(d, eps_d)
    for p in O_DN + O_UP:          # 面内 O
        onsite(p, eps_p)
    onsite(O_AP[0], eps_p_z)       # 顶点 O

    #===== Ni 3d e_g 的 d-d 相互作用：Racah 在 e_g² 子空间的有效值 =====
    #  ³A₂g = U' - J = A - 8B，¹Eg = U' + J = A + 2C，¹A₁g = U + J = A + 8B + 4C
    U_ni, U_prime, J_H = DD_RACAH

    for o in NI + O_ALL:           # 同轨道双占
        u = U_ni if o in NI else u_p
        add(f"+_{mode(o, 0)} -_{mode(o, 0)} "
            f"+_{mode(o, 1)} -_{mode(o, 1)}", u)

    for a, b in ((0, 1), (2, 3)):
        modes = [mode(a, 0), mode(a, 1), mode(b, 0), mode(b, 1)]
        local = {}

        def add_local(key, val):
            local[key] = local.get(key, 0.0) + val

        for s1 in (0, 1):          # 异轨道
            for s2 in (0, 1):
                u = (U_prime - J_H) if s1 == s2 else U_prime
                add_local(f"+_{mode(a, s1)} -_{mode(a, s1)} "
                          f"+_{mode(b, s2)} -_{mode(b, s2)}", u)
        # 自旋翻转（Sz 守恒）：a↓b↑ <-> a↑b↓
   
        add_local(f"+_{mode(a, 0)} -_{mode(a, 1)} "
                  f"+_{mode(b, 1)} -_{mode(b, 0)}", -J_H)
        add_local(f"+_{mode(a, 1)} -_{mode(a, 0)} "
                  f"+_{mode(b, 0)} -_{mode(b, 1)}", -J_H)
        # 配对跳跃：a↑a↓ <-> b↑b↓
        add_local(f"+_{mode(a, 0)} +_{mode(a, 1)} "
                  f"-_{mode(b, 1)} -_{mode(b, 0)}", J_H)
        add_local(f"+_{mode(b, 0)} +_{mode(b, 1)} "
                  f"-_{mode(a, 1)} -_{mode(a, 0)}", J_H)

        for key, val in local.items():
            add(key, val)

    if d8_only:                    # 只保留 d⁸(h=2) 的相互作用
        G = U_ni + 2.0 * U_prime - J_H
        for a, b in ((0, 1), (2, 3)):
            cut_interaction_above_d8(
                terms, [mode(a, 0), mode(a, 1), mode(b, 0), mode(b, 1)], G)

    if d_level:                    # d^n 能级 A + 1.5A|h-2|
        for a, b in ((0, 1), (2, 3)):
            modes = [mode(a, 0), mode(a, 1), mode(b, 0), mode(b, 1)]
            add_h_function(terms, modes,
                           [abs(h - 2) * 1.5 * A for h in range(5)])

    ###################################################################3
    # 能量零点平移：H -> H - E_ref
    if reference_zero is not None:
        terms[""] = terms.get("", 0.0) - reference_zero

    if return_terms:
        # 直接给出项字典（精确对角化用；跳过 JW 映射，省时间）
        return terms

    from qiskit_nature.second_q.mappers import JordanWignerMapper
    from qiskit_nature.second_q.operators import FermionicOp

    ferm_op = FermionicOp(terms, num_spin_orbitals=N_MODES).normal_order()
    qubit_op = JordanWignerMapper().map(ferm_op)
    return ferm_op, qubit_op


# ---------------- JW 结果展示 ----------------
_ORB_LABEL = {
    0: "Ni0.dx2", 1: "Ni0.dz2", 2: "Ni1.dx2", 3: "Ni1.dz2",
    4: "Odn.px1", 5: "Odn.px2", 6: "Odn.py1", 7: "Odn.py2",
    8: "Oup.px1", 9: "Oup.px2", 10: "Oup.py1", 11: "Oup.py2",
    12: "Oap.pz",
}


def mode_label(m: int) -> str:
    """自旋轨道编号 -> '轨道名+自旋'。"""
    orb, spin = divmod(m, 2)
    return f"{_ORB_LABEL.get(orb, f'orb{orb}')}{'↑' if spin == 0 else '↓'}"


def format_jw_hamiltonian(qubit_op, sort_by_magnitude: bool = True) -> str:
    """把 SparsePauliOp 渲染成表格（qubit 0 在最左）。"""
    coeffs = np.asarray(qubit_op.coeffs)
    labels_all = [str(p) for p in qubit_op.paulis]
    order = np.argsort(-np.abs(coeffs)) if sort_by_magnitude else np.arange(len(coeffs))

    n = qubit_op.num_qubits
    head = (f"JW 变换后的哈密顿量（{len(coeffs)} 项，{n} qubit，"
            f"串按 qubit 0 → {n - 1} 从左到右排列）")
    lines = [head, "=" * len(head), ""]
    lines.append(f"{'#':>4}  {'coeff':>13}  {'Pauli string (q0..q%d)' % (n - 1)}"
                 f"  {'作用位置':<}")
    lines.append("-" * 108)

    for rank, k in enumerate(order, start=1):
        c = coeffs[k]
        s = labels_all[k][::-1]
        cs = f"{c.real:+.6f}" if abs(c.imag) < 1e-12 else f"{c:.6f}"
        xy = [m for m, ch in enumerate(s) if ch in "XY"]
        zs = [m for m, ch in enumerate(s) if ch == "Z"]
        parts = []
        if xy:
            parts.append("X/Y: " + ", ".join(f"{m}({mode_label(m)})" for m in xy))
        if zs:
            if len(zs) == 1:
                parts.append(f"Z: {zs[0]}({mode_label(zs[0])})")
            else:
                parts.append(f"Z链: {zs[0]}({mode_label(zs[0])})"
                             f"→{zs[-1]}({mode_label(zs[-1])}) 共{len(zs)}个")
        annot = "    ".join(parts) if parts else "常数项"
        lines.append(f"{rank:>4}  {cs:>13}  {s}  {annot}")
    return "\n".join(lines)


if __name__ == "__main__":
    ferm_op, qubit_op = build_hamiltonian()
    print(f"空间轨道: {N_ORB}  自旋轨道(qubit): {N_MODES}")
    print(f"FermionicOp 项数: {len(ferm_op)}")
    print(f"SparsePauliOp 项数: {len(qubit_op)}")
    print(f"qubit 数: {qubit_op.num_qubits}")
    text = format_jw_hamiltonian(qubit_op)
    print()
    print(text)
    os.makedirs("results", exist_ok=True)
    with open(os.path.join("results", "jw_hamiltonian.txt"), "w",
              encoding="utf-8") as f:
        f.write(text + "\n")
    print("\n已保存: results/jw_hamiltonian.txt")
