"""La3Ni2O7 双层模型（Hamitonian.py，13 轨道 / 26 qubit）精确对角化。

用法:
  py exact_solve.py                 5 空穴：全 Sz 扫描 + 构型/自旋能隙对照（results/exact_5h.txt）
  py exact_solve.py --holes 4       4 空穴：同上（results/exact_4h.txt）
  py exact_solve.py --holes 6       6 空穴：同上（results/exact_6h.txt）
  py exact_solve.py --compare       4+6 与 5+5 空穴序的能量差（不写文件）
  py exact_solve.py --config        只解默认扇区，打印构型表（不写文件）
  py exact_solve.py --config 2 3    同上，但指定 (N↑,N↓)
"""
import importlib.util
import math
import os
import sys
import time
import unicodedata
from dataclasses import dataclass
from itertools import combinations

import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import eigsh

N_HOLES = 6                      # 默认空穴数（可用 --holes 覆盖）
HAM_FILE = "Hamitonian.py"
RESULTS_DIR = "results"


def out_file(n_holes):
    """报表输出路径：results/exact_{N_h}h.txt（按空穴数分文件）。"""
    return os.path.join(RESULTS_DIR, f"exact_{n_holes}h.txt")


OUT_FILE = out_file(N_HOLES)


def default_sector(n_holes):
    """默认基态扇区：多余的一个空穴给 ↑ 自旋。"""
    n_dn = n_holes // 2
    return n_holes - n_dn, n_dn

N_PRINT = 4
N_LOW_STATES = 12
DEG_TOL = 1e-8
COEF_EPS = 1e-12
DENSE_MAX_DIM = 2000      # ≤ 此维度用稠密 eigh，更大一律走 ARPACK（稠密是 O(dim³)）

ORBIT_FAMILIES = [
    ("dx2/下", [0]), ("dz2/下", [1]), ("dx2/上", [2]), ("dz2/上", [3]),
    ("px/下", [4, 5]), ("py/下", [6, 7]),
    ("px/上", [8, 9]), ("py/上", [10, 11]),
    ("pz", [12]),
]


def popcount(x):
    """位串中 1 的个数（SWAR，逐位并行）。

    numpy 数组和 Python 整数都可用，向量化后是 build_sector 的主要热路径之一。
    这里数的是"有几个空穴"：位串第 m 位为 1 表示 mode m 上有空穴。
    """
    x = x - ((x >> 1) & 0x55555555)
    x = (x & 0x33333333) + ((x >> 2) & 0x33333333)
    x = (x + (x >> 4)) & 0x0F0F0F0F
    return (x * 0x01010101) >> 24 & 0xFF


def comb_masks(modes, k):
    """枚举 modes 中取 k 个的全部组合，每组合返回一个位串。

    用来生成一个 Sz 扇区的基矢：上自旋取 n_up 个、下自旋取 n_dn 个，
    两者按位 OR 就是一个 Fock 态。
    """
    for picked in combinations(modes, k):
        bits = 0
        for mode in picked:
            bits |= 1 << mode
        yield bits


def parity_sign(states, mode):
    """费米子交换符号：在 mode 上作用算符时，该 mode 之前有几个已占据。

    a_m |n₀n₁…⟩ = (−1)^{Σ_{j<m} n_j} |…⟩，所以用低位 1 的个数定符号。
    返回 int8 数组 ±1，与 states 同长度。
    """
    lower_bits = states & ((1 << mode) - 1)
    return np.where(popcount(lower_bits) & 1, -1, 1).astype(np.int8)


def parse_operator_string(text):
    """把 "+_3 -_7" 解析成 [("+", 3), ("-", 7)]。

    terms 的键是字符串，作用到基矢上时转成 (符号, mode) 列表；符号只能是
    "+"（产生空穴）或 "-"（湮灭空穴）。
    """
    return [(token.split("_")[0], int(token.split("_")[1]))
            for token in text.split()]


class Model:
    """一个已经构造好的哈密顿量：好用的索引 + 项列表。

    由 load_model() 返回，供 build_sector / 统计层消费。
    """

    def __init__(self, module, terms, params=None):
        """从 Hamitonian 模块和项字典建立模型。

        module 提供 N_ORB（13 个空间轨道），terms 是 {算符串: 系数}，
        params 记录调用时的参数快照（用于复现，不含 reference_zero）。
        项列表这里就从字符串解析成 (算符, 系数)，后面不再解析。
        """
        self.n_orb = module.N_ORB
        self.n_modes = 2 * self.n_orb          # 26 个自旋轨道 = qubit 数
        self.terms = [(parse_operator_string(k), v) for k, v in terms.items()]
        self.params = dict(params or {})
        self.up_modes = list(range(0, self.n_modes, 2))    # 偶数 mode = ↑
        self.dn_modes = list(range(1, self.n_modes, 2))    # 奇数 mode = ↓


def load_module(path=HAM_FILE):
    """按文件路径把 Hamitonian.py 作为模块加载（不要求它在 sys.path 里）。

    用 importlib 而不是 import，避免包名冲突；返回值提供 N_ORB、N_MODES
    和 build_hamiltonian。
    """
    spec = importlib.util.spec_from_file_location("ham", path)
    module = importlib.util.module_from_spec(spec)
    sys.modules["ham"] = module
    spec.loader.exec_module(module)
    return module


def reference_energy(n_holes=4, path=HAM_FILE, **overrides):
    """用同一套参数自己算出参考能量：n_holes 空穴在所有 Sz 扇区里的最低本征值。

    默认 n_holes=4，即 29.5 GPa 下的 d⁸-d⁸ 基态（t_pd = 1.58）。
    """
    model = load_model(path=path, **overrides)
    best = math.inf
    for n_up, n_dn in iter_sectors(model, n_holes):
        H, _, _ = build_sector(model, n_up, n_dn)
        e0 = solve_weights(H)[0][0]
        if e0 < best:
            best = e0
    return best


def load_model(path=HAM_FILE, reference_zero=None, **overrides):
    """构建模型；overrides 直接透传给 Hamitonian.build_hamiltonian。

    Hamitonian.build_hamiltonian 直接用 return_terms=True 给出项字典
    （空穴语言），这里只做零点平移和系数截断。

    reference_zero 除了给具体数值，也可以是字符串 "auto"：此时用同一套参数
    现算 4 空穴基态（见 reference_energy）作为零点。
    """
    module = load_module(path)
    if reference_zero == "auto":
        reference_zero = reference_energy(path=path, **overrides)
    terms = dict(module.build_hamiltonian(return_terms=True, **overrides))
    if reference_zero is not None:
        terms[""] = terms.get("", 0.0) - reference_zero
    terms = {k: v for k, v in terms.items() if abs(v) > COEF_EPS}
    return Model(module, terms, overrides)


def apply_operator_string(states, operators):
    """把一项算符串作用到整批 Fock 态上（全向量化，无 Python 逐态循环）。

    参数
    ----
    states : int64 数组，每个元素是一个位串（第 m 位 = mode m 上有空穴）
    operators : [(符号, mode)]，符号 "+" 产生、"−" 湮灭，按从右到左作用

    返回
    ----
    (new_states, signs, illegal)

    signs 是费米子交换累积出的 ±1；illegal 标记该态上这一步做不下去的情形
    （产生已占的轨道，或湮灭空的轨道），这些分量在拼矩阵时会被丢弃。
    这是"算符串 → 矩阵元素"的唯一入口。
    """
    new_states = states.copy()
    signs = np.ones(len(states), dtype=np.int8)
    illegal = np.zeros(len(states), dtype=bool)

    for symbol, mode in reversed(operators):
        bit = 1 << mode
        occupied = (new_states & bit) != 0
        sign = parity_sign(new_states, mode)

        if symbol == "+":
            illegal |= occupied
            new_states = new_states | bit
        else:
            illegal |= ~occupied
            new_states = new_states & ~bit

        signs *= sign

    return new_states, signs, illegal


def build_sector(model, n_up, n_dn):
    """在固定 (N↑, N↓) 扇区里把哈密顿量拼成稀疏矩阵。

    步骤
    ----
    1. 枚举基矢：上自旋从 13 个 ↑ mode 里选 n_up 个，下自旋选 n_dn 个，按位 OR
       得到 Fock 位串。维度 = C(13,n_up) × C(13,n_dn)。
    2. 对每一项算符串，用 apply_operator_string 批量算出目标位串和费米子符号。
    3. 用 searchsorted 在排好序的基矢里定位行号，只保留仍在本扇区内的跃迁
       （跨扇区的分量被丢弃，这正是 H 按 N↑、N↓ 严格块对角的原因）。

    返回 (H_csr, masks, index_of)。masks[i] 是第 i 个基矢的位串，
    index_of 是位串到行号的字典（⟨S²⟩ 计算要用）。
    """
    masks = np.array([up_bits | dn_bits
                      for up_bits in comb_masks(model.up_modes, n_up)
                      for dn_bits in comb_masks(model.dn_modes, n_dn)],
                     dtype=np.int64)
    index_of = {int(mask): i for i, mask in enumerate(masks)}
    dim = len(masks)

    # 用 searchsorted 批量查行号，代替逐态 dict.get，把 O(项数×维数) 的 Python
    # 循环压成 O(项数) 次向量化操作
    order = np.argsort(masks, kind="stable")
    sorted_masks = masks[order]

    rows, cols, values = [], [], []
    for operators, coefficient in model.terms:
        new_states, signs, illegal = apply_operator_string(masks, operators)
        valid = ~illegal
        if not valid.any():
            continue
        col = np.nonzero(valid)[0]
        target = new_states[col]
        pos = np.searchsorted(sorted_masks, target)
        np.clip(pos, 0, dim - 1, out=pos)
        hit = sorted_masks[pos] == target
        if not hit.any():
            continue
        rows.append(order[pos[hit]])
        cols.append(col[hit])
        values.append(signs[col][hit] * coefficient)

    if values:
        H = coo_matrix((np.concatenate(values),
                        (np.concatenate(rows), np.concatenate(cols))),
                       shape=(dim, dim)).tocsr()
    else:
        H = coo_matrix((dim, dim)).tocsr()
    return H, masks, index_of


def iter_sectors(model, n_holes=N_HOLES):
    """枚举 N_h 个空穴的全部自旋扇区 (N↑, N↓)，只取 N↑ ≥ N↓（Sz ≥ 0）。

    N 个空穴守恒，但自旋可以按任意方式分配；Sz < 0 的扇区与对应正扇区能量
    相同（SU(2) 对称），不必重复求解。5 空穴给出 (5,0)、(4,1)、(3,2) 三个扇区。
    """
    max_up = min(model.n_orb, n_holes)
    min_up = max(0, n_holes - model.n_orb)
    for n_up in range(max_up, min_up - 1, -1):
        n_dn = n_holes - n_up
        if n_up >= n_dn:
            yield n_up, n_dn


@dataclass
class Sector:
    """一个自旋扇区的全部结果：矩阵、基矢、本征值、简并子空间。

    只解一次，后面算 ⟨S²⟩ 和构型权重都复用这里的缓存，避免重复对角化。
    """

    n_up: int                      # 上自旋空穴数
    n_dn: int                      # 下自旋空穴数
    H: object                      # 该扇区的 CSR 稀疏矩阵
    masks: np.ndarray              # 基矢位串
    index_of: dict                 # 位串 -> 行号
    evals: np.ndarray              # 最低的若干本征值
    weight: np.ndarray = None      # 简并子空间平均后的 |v|²（每个基矢一份）
    degeneracy: int = 1            # 基态简并度
    vectors: np.ndarray = None     # 简并子空间的本征矢（dim × g）

    @property
    def dim(self):
        """扇区维数（基矢个数）。"""
        return self.H.shape[0]

    @property
    def sz(self):
        """总自旋 z 分量 Sz = (N↑ − N↓)/2。"""
        return (self.n_up - self.n_dn) / 2


def krylov_size(dim, k):
    """ARPACK 的 Krylov 子空间维数 ncv。

    取 max(2k+8, 60)，并保证 ncv ≤ dim−1（ARPACK 的硬性要求），
    避免小扇区上 ncv 超过维数导致收敛失败。
    """
    return min(dim - 1, max(2 * k + 8, 60))


def solve_subspace(H, k=N_LOW_STATES):
    """最低 k 个本征态；返回 (本征值, 简并子空间平均权重, 简并度, 简并子空间本征矢)。"""
    dim = H.shape[0]
    if dim <= DENSE_MAX_DIM:
        eigenvalues, eigenvectors = np.linalg.eigh(H.toarray())
    else:
        kk = max(1, min(k, dim - 2))
        eigenvalues, eigenvectors = eigsh(H, k=kk, which="SA",
                                          ncv=krylov_size(dim, kk))
        order = np.argsort(eigenvalues)
        eigenvalues, eigenvectors = eigenvalues[order], eigenvectors[:, order]
    degeneracy = int(np.sum(np.abs(eigenvalues - eigenvalues[0]) < DEG_TOL))
    subspace = eigenvectors[:, :degeneracy]
    weight = np.mean(np.abs(subspace) ** 2, axis=1)
    return eigenvalues, weight, degeneracy, subspace


def solve_weights(H, k=N_LOW_STATES):
    """最低若干本征态；返回 (本征值, 简并子空间上的平均权重, 简并度)。"""
    eigenvalues, weight, degeneracy, _ = solve_subspace(H, k)
    return eigenvalues, weight, degeneracy


def solve_spectra(model, out, n_holes=N_HOLES):
    """每个 Sz 扇区只解一次，并缓存简并子空间本征矢供后续统计。"""
    spectra = []
    for n_up, n_dn in iter_sectors(model, n_holes):
        H, masks, index_of = build_sector(model, n_up, n_dn)
        evals, weight, degeneracy, subspace = solve_subspace(H)
        sector = Sector(n_up, n_dn, H, masks, index_of, evals,
                        weight, degeneracy, subspace)
        spectra.append(sector)
        out(f"  (N↑,N↓)=({n_up},{n_dn})  Sz={sector.sz:+.1f}  dim={sector.dim:>7,}  "
            f"E = {', '.join(f'{e:10.6f}' for e in sector.evals[:N_PRINT])}")
    return spectra


def total_spin_square(model, sector):
    """⟨S²⟩（简并子空间平均）。

    S⁺ = Σ_i h†_{i↑} h_{i↓}（hole 语言），S² = Sz(Sz+1) + S⁻S⁺，
    故 ⟨S²⟩ = Sz(Sz+1) + ‖S⁺|v⟩‖²。
    """
    sz = sector.sz
    if sector.vectors is None:
        return sz * (sz + 1)

    n_up, n_dn = sector.n_up + 1, sector.n_dn - 1
    if n_up > model.n_orb or n_dn < 0:
        return sz * (sz + 1)

    _, _, index_up = build_sector(model, n_up, n_dn)
    degeneracy = sector.vectors.shape[1]
    raised = np.zeros((len(index_up), degeneracy))

    for i in range(model.n_orb):
        operators = [("+", 2 * i), ("-", 2 * i + 1)]
        new_states, signs, illegal = apply_operator_string(sector.masks, operators)
        for col in np.nonzero(~illegal)[0]:
            row = index_up.get(int(new_states[col]))
            if row is not None:
                raised[row] += signs[col] * sector.vectors[col, :]

    return sz * (sz + 1) + float(np.sum(raised ** 2)) / degeneracy


def mode_groups(families=ORBIT_FAMILIES):
    """把 9 个轨道族按自旋拆成 18 个统计组：(组名, 该组的 mode 列表)。

    分组是"物理记账"的最小单位：每个 Ni 的 dx2、dz2 各自一组，面内 O 的
    px、py 各自一组，顶位 pz 一组，每组再分 ↑、↓。
    """
    groups = []
    for name, orbitals in families:
        for spin, spin_tag in ((0, "↑"), (1, "↓")):
            groups.append((f"{name}{spin_tag}", [2 * o + spin for o in orbitals]))
    return groups


def block_modes(groups):
    """把 18 个组归并成 5 个物理块，供构型标签使用。

    五个块是：下层 Ni、上层 Ni、下层面内 O、上层面内 O、顶位 pz。
    返回 {块名: [组下标]}，下标对应 groups 的顺序。
    """
    blocks = {
        "Ni下": ["dx2/下", "dz2/下"],
        "Ni上": ["dx2/上", "dz2/上"],
        "O面内下": ["px/下", "py/下"],
        "O面内上": ["px/上", "py/上"],
        "pz": ["pz"],
    }
    position = {name: i for i, (name, _) in enumerate(groups)}
    return {block: [position[f"{family}{spin}"] for family in families
                    for spin in ("↑", "↓")]
            for block, families in blocks.items()}


def count_in(key, indexes):
    """把 key 里指定几个位置的值加起来（key 是各组的空穴数元组）。

    例如 count_in(key, blocks["Ni上"]) 就是上层 Ni 两个轨道上的 e_g 空穴数 h。
    """
    return sum(key[i] for i in indexes)


def type_weights(model, masks, basis_weight, groups):
    """把基矢按"各组空穴数"归并，累加权重与组态数。

    masks 是扇区基矢，basis_weight 是每个基矢的 |v|²（简并已平均）。
    每条基矢数出 18 个组的空穴数，得到一个 18 元组作为 key，把同 key 的
    权重相加。返回 (按权重降序的 [(key, weight)], {key: 组态个数})。
    """
    bits = ((masks[:, None] >> np.arange(model.n_modes)) & 1).astype(bool)
    hole_counts = np.stack([bits[:, modes].sum(axis=1) for _, modes in groups],
                           axis=1)

    weight_sum, basis_count = {}, {}
    for counts, weight in zip(hole_counts.tolist(), basis_weight.tolist()):
        key = tuple(counts)
        weight_sum[key] = weight_sum.get(key, 0.0) + weight
        basis_count[key] = basis_count.get(key, 0) + 1

    ordered = sorted(weight_sum.items(), key=lambda item: -item[1])
    return ordered, basis_count


def layered_weights(ordered, basis_count, blocks=None):
    """按 (h↑, h↓, L↑, L↓, O) 分层标签合并权重与组态数。"""
    if blocks is None:
        blocks = block_modes(mode_groups())
    merged, counts = {}, {}
    for key, weight in ordered:
        label = (count_in(key, blocks["Ni上"]), count_in(key, blocks["Ni下"]),
                 count_in(key, blocks["O面内上"]), count_in(key, blocks["O面内下"]),
                 count_in(key, blocks["pz"]))
        merged[label] = merged.get(label, 0.0) + weight
        counts[label] = counts.get(label, 0) + basis_count[key]
    return merged, counts


class Report:
    """报表输出器：一边打印到屏幕，一边把同样的内容攒起来写文件。

    所有 report_* 函数都接收一个 out 对象并调用 out("...")，因此既能直接
    print，也能在测试里传一个假对象截获文本。
    """

    def __init__(self):
        """建一个空缓冲；每行文本都会追加到 self.lines。"""
        self.lines = []

    def __call__(self, text=""):
        """输出一行：立即打印，同时缓存，供 save() 落盘。"""
        print(text)
        self.lines.append(text)

    def save(self, path=OUT_FILE):
        """把缓存的所有行写入文件（自动建 results/ 目录）。"""
        directory = os.path.dirname(path)
        if directory:
            os.makedirs(directory, exist_ok=True)
        with open(path, "w", encoding="utf-8") as f:
            f.write("\n".join(self.lines) + "\n")
        print(f"\n已保存: {path}")


def report_header(out, model, n_holes=N_HOLES):
    """报表抬头：轨道数、qubit 数、空穴数与整个希尔伯特空间的维数 C(26, N_h)。"""
    out("=" * 78)
    out(f"La3Ni2O7 双层模型精确对角化   {model.n_orb} 轨道 / "
        f"{model.n_modes} 自旋轨道(qubit)")
    out(f"空穴数 N_h = {n_holes}     扇区维度 "
        f"C({model.n_modes},{n_holes}) = {math.comb(model.n_modes, n_holes):,}")
    out("=" * 78)


def report_ground(out, sector, E0, degeneracy):
    """全局基态一行：所在扇区、Sz、能量、维数与简并度。"""
    out(f"\n全局基态: (N↑,N↓)={(sector.n_up, sector.n_dn)}  Sz={sector.sz:+.1f}   "
        f"E0 = {E0:.9f}   dim = {sector.dim:,}   简并度 g = {degeneracy}")


def spin_from_square(s2):
    """由 ⟨S²⟩ 反解 S：S(S+1) = s2 的物理根，S = (−1 + √(1+4s2))/2。

    用来把数值上算出的 ⟨S²⟩ 换成整数或半整数的自旋量子数，
    例如 s2 = 0.75 → S = 0.5，s2 = 3.75 → S = 1.5。
    """
    return (-1.0 + math.sqrt(1.0 + 4.0 * s2)) / 2.0


def report_spin_states(out, spectra, model, n_holes=N_HOLES):
    """自旋态表：每个 Sz 扇区的最低能量、⟨S²⟩、相对基态的激发能。

    有参考读数时并排列出"参考/差"，差全为 0 说明自旋能隙与参考程序一致。
    ⟨S²⟩ 必须是 S(S+1)，这是检查 Kanamori 符号最快的判据。
    """
    ground = min(spectra, key=lambda sector: sector.evals[0])
    E0 = ground.evals[0]
    ref_gaps = REF_SPIN_GAPS.get(n_holes, {})

    out("\n自旋态：各 Sz 扇区的最低能量（S 由 ⟨S²⟩ 判定）")
    if ref_gaps:
        out("  参考程序（t_pd = 1.58，29.5 GPa）相对基态的激发能")
    header = f"  {'S':>6}{'⟨S²⟩':>9}{'Sz':>7}{'E':>15}{'ΔE (eV)':>11}"
    out(header + f"{'参考':>9}{'差':>10}" if ref_gaps else header)

    for sector in sorted(spectra, key=lambda s: s.sz):
        s2 = total_spin_square(model, sector)
        spin = spin_from_square(s2)
        gap = sector.evals[0] - E0
        line = (f"  {spin:>6.2f}{s2:>9.4f}{sector.sz:>7.1f}"
                f"{sector.evals[0]:>15.6f}{gap:>11.4f}")
        if ref_gaps:
            value = ref_gaps.get(round(spin, 1))
            line += (f"{value:>9.4f}{gap - value:>+10.4f}" if value is not None
                     else f"{'-':>9}{'-':>10}")
        out(line)

    if ref_gaps:
        out("  注：'-' 表示该扇区就是基态本身（ΔE = 0），或参考程序无对应读数。")


# 参考程序（5hole 的代码，t_pd = 1.58、29.5 GPa、t_dO = 1.66）算出的基态构型权重。
# 4/5/6 空穴分别在各自基态扇区上统计（Sz = 0 / 0.5 / 0），简并子空间对 |v|² 平均。
# 由 ref_program.py 生成；改了参数就重跑 `python ref_program.py 4 5 6` 更新。
REF_WEIGHTS = {
    4: {"d8-d8": 0.58486, "d8-d9L": 0.26689, "d8-O-d9": 0.07059,
        "d9L-d9L": 0.03045, "d9-O-d9L": 0.01673},
    5: {"d8-d8L": 0.40761, "d7-d8": 0.18565, "d8L-d9L": 0.09697,
        "d8-d9L2": 0.08213, "d8-O-d9L": 0.04882, "d7-d9L": 0.04625,
        "d8-O-d8": 0.03920},
    6: {"d8L-d8L": 0.27565, "d7-d8L": 0.24598, "d8L-d9L2": 0.11195,
        "d8L-O-d9L": 0.05605, "d7-d7": 0.05481, "d8-O-d8L": 0.05629,
        "d7-O-d9L": 0.02428},
}

# 分层标签 = (d空穴_上, d空穴_下, 面内L_上, 面内L_下, 顶点O空穴)；上下层镜像同时列出
CONFIG_LABELS = {
    4: {
        "d8-d8": ((2, 2, 0, 0, 0),),
        "d8-d9L": ((2, 1, 0, 1, 0), (1, 2, 1, 0, 0)),
        "d8-O-d9": ((2, 1, 0, 0, 1), (1, 2, 0, 0, 1)),
        "d9L-d9L": ((1, 1, 1, 1, 0),),
        "d9-O-d9L": ((1, 1, 0, 1, 1), (1, 1, 1, 0, 1)),
    },
    5: {
        "d8-d8L": ((2, 2, 0, 1, 0), (2, 2, 1, 0, 0)),
        "d7-d8": ((3, 2, 0, 0, 0), (2, 3, 0, 0, 0)),
        "d8L-d9L": ((2, 1, 1, 1, 0), (1, 2, 1, 1, 0)),
        "d8-d9L2": ((2, 1, 0, 2, 0), (1, 2, 2, 0, 0)),
        "d8-O-d9L": ((2, 1, 0, 1, 1), (1, 2, 1, 0, 1)),
        "d7-d9L": ((3, 1, 0, 1, 0), (1, 3, 1, 0, 0)),
        "d8-O-d8": ((2, 2, 0, 0, 1),),
    },
    6: {
        "d8L-d8L": ((2, 2, 1, 1, 0),),
        "d7-d8L": ((3, 2, 0, 1, 0), (2, 3, 1, 0, 0)),
        "d8L-d9L2": ((2, 1, 1, 2, 0), (1, 2, 2, 1, 0)),
        "d8L-O-d9L": ((2, 1, 1, 1, 1), (1, 2, 1, 1, 1)),
        "d7-d7": ((3, 3, 0, 0, 0),),
        "d8-O-d8L": ((2, 2, 0, 1, 1), (2, 2, 1, 0, 1)),
        "d7-O-d9L": ((3, 1, 0, 1, 1), (1, 3, 1, 0, 1)),
    },
}

# 参考程序给出的自旋激发能（各 Sz 扇区相对基态，eV），键为总自旋 S。
# 同样由 ref_program.py 生成；S 由 |Sz| 定（各扇区最低态）。
REF_SPIN_GAPS = {
    4: {1: 0.032274, 2: 0.104956},
    5: {1.5: 0.188533, 2.5: 4.493914},
    6: {1: 0.313579, 2: 4.611133, 3: 8.916512},
}


def pad(text, width, align="<"):
    """按终端显示宽度对齐（中文字符占 2 列）。"""
    display = sum(2 if unicodedata.east_asian_width(ch) in "WF" else 1 for ch in text)
    blank = " " * max(0, width - display)
    return blank + text if align == ">" else text + blank


def config_label(h_up, h_dn, l_up, l_dn, apical=0):
    """把分层标签 (h↑, h↓, L↑, L↓, O) 渲染成人类可读的构型名。

    该层的化学式记作 d^{10−h}L^L：h 是 e_g 空穴数，所以空穴 2 个就是 d⁸。
    中间若顶点氧上有空穴，插入一个 O。例如 (2,2,0,1,0) → "d8 - d8L"，
    表示一层是 d⁸、另一层是 d⁸ 带一个面内配体空穴。
    """

    def layer(h, ligand):
        """单层的 d^{10-h}L^ligand 写法。"""
        tag = f"d{10 - h}"
        if ligand:
            tag += "L" if ligand == 1 else f"L^{ligand}"
        return tag

    parts = [layer(h_up, l_up)]
    if apical:
        parts.append("O" if apical == 1 else f"O^{apical}")
    parts.append(layer(h_dn, l_dn))
    return " - ".join(parts)


def report_configurations(out, ordered, basis_count, blocks, n_holes=N_HOLES):
    """构型占比表：按 (h↑, h↓, L↑, L↓, O) 合并权重，与参考程序并排对比。

    有参考读数时，把 CONFIG_LABELS 里同名的若干镜像标签加在一起再比，
    所以 d8-d8L 的 0.408 是上下层两种镜像之和。没有参考读数（其他空穴数）
    时只打印权重最高的 10 个构型。
    """
    merged, counts = layered_weights(ordered, basis_count, blocks)
    ref = REF_WEIGHTS.get(n_holes)
    labels_of = CONFIG_LABELS.get(n_holes, {})

    out("\n构型占比（分层标签，不合并上下层）")
    out("  标签 = 上层 d^{10-h}L - 顶点 O - 下层 d^{10-h}L；"
        "h 为该层 e_g 空穴数，L = 面内 O 空穴，O = 顶点 O 空穴")
    out("  分层标签 (h↑,h↓,L↑,L↓,O) = 上/下层 e_g 空穴数、上/下层面内 O 空穴数、顶点 O 空穴数")

    if ref is None:
        out(f"\n  前 10 个构型（N_h = {n_holes}，无参考程序读数可比）")
        out("  " + pad("构型", 12) + pad("权重", 8, ">")
            + pad("组态数", 9, ">") + "   标签")
        for label, weight in sorted(merged.items(), key=lambda item: -item[1])[:10]:
            out(f"  {config_label(*label):<12}{weight:>8.3f}"
                f"{counts[label]:>9,}   ({','.join(str(v) for v in label)})")
        return

    out(f"\n  参考程序（t_pd = 1.58，29.5 GPa）算出的构型")
    out("  " + pad("构型", 12) + pad("参考", 9, ">") + pad("本文", 9, ">")
        + pad("差", 10, ">") + pad("组态数", 9, ">") + "   标签")
    ref_total = calc_total = 0.0
    for name, value in ref.items():
        labels = labels_of[name]
        weight = sum(merged.get(label, 0.0) for label in labels)
        count = sum(counts.get(label, 0) for label in labels)
        ref_total += value
        calc_total += weight
        tags = "  ".join(f"({','.join(str(v) for v in label)})" for label in labels)
        out(f"  {name:<12}{value:>9.5f}{weight:>9.5f}{weight - value:>+10.5f}"
            f"{count:>9,}   {tags}")
    out("  " + pad("合计", 12) + f"{ref_total:>9.5f}{calc_total:>9.5f}"
        f"{calc_total - ref_total:>+10.5f}")
    out("  注：'合计' 只统计上表列出的构型，未列出的构型不在表内。")


def report_charge_order(out, model):
    """4+6 与 5+5 空穴序的能量差 ΔE = E0(6) + E0(4) − 2·E0(5)。

    论文正文用这个量讨论棋盘型电荷序与均匀 5+5 序的竞争：
    ΔE 在所有压力下为负，但量级 ≤ 10 meV。三个基态都落在默认扇区。
    """
    out("\n电荷序：4+6 与 5+5 空穴团簇的能量比较")
    energies = {}
    for n_holes in (4, 5, 6):
        n_up, n_dn = default_sector(n_holes)
        H, _, _ = build_sector(model, n_up, n_dn)
        evals, _, _ = solve_weights(H, k=2)
        energies[n_holes] = float(evals[0])
        out(f"  E0(N_h={n_holes}) = {energies[n_holes]:.9f}   "
            f"(N↑,N↓)=({n_up},{n_dn})   dim={H.shape[0]:,}")
    delta = energies[6] + energies[4] - 2.0 * energies[5]
    out(f"  ΔE = E0(6) + E0(4) − 2·E0(5) = {delta * 1000:+.2f} meV")
    out("  论文：ΔE < 0 且 |ΔE| ≲ 10 meV，即 4+6 棋盘序只以约 10 meV 的尺度占优；")
    out("        该量级不足以压过团簇间跳跃等修正，原文据此不下定论。")


def quick_config_report(out, model, n_up=None, n_dn=None, n_holes=N_HOLES):
    """只解一个扇区，打印分层构型表（快速核对参数用）。"""
    if n_up is None or n_dn is None:
        n_up, n_dn = default_sector(n_holes)

    H, masks, _ = build_sector(model, n_up, n_dn)
    eigenvalues, basis_weight, degeneracy = solve_weights(H)

    out(f"快速构型对照: N_h = {n_up + n_dn}  (N↑,N↓)=({n_up},{n_dn})  "
        f"Sz={(n_up - n_dn) / 2:+.1f}  dim={H.shape[0]:,}  "
        f"E0 = {eigenvalues[0]:.9f}  g = {degeneracy}")

    groups = mode_groups()
    ordered, basis_count = type_weights(model, masks, basis_weight, groups)
    report_configurations(out, ordered, basis_count, block_modes(groups),
                          n_up + n_dn)


def parse_args(args):
    """解析 --holes N、--config [n↑ n↓]、--compare。

    返回 (空穴数, 是否快速模式, 扇区, 是否做电荷序比较)。
    """
    n_holes, quick, sector, compare, index = N_HOLES, False, None, False, 0
    while index < len(args):
        arg = args[index]
        if arg == "--holes":
            if index + 1 >= len(args):
                raise ValueError("--holes 需要一个整数")
            n_holes = int(args[index + 1])
            index += 2
        elif arg in ("-c", "--config"):
            quick = True
            rest = args[index + 1:]
            if len(rest) >= 2 and rest[0].isdigit() and rest[1].isdigit():
                sector = (int(rest[0]), int(rest[1]))
                index += 3
            else:
                index += 1
        elif arg == "--compare":
            compare = True
            index += 1
        else:
            raise ValueError(f"未知参数 {arg}")
    return n_holes, quick, sector, compare


def main(argv=None):
    """命令行入口：解析参数、建模型、按模式分派。

    四种走法：
      --config [n↑ n↓]  只解一个扇区，打印构型表（不写文件）
      --compare         算 4+6 与 5+5 的能量差，看电荷序竞争
      默认              扫该空穴数的全部 Sz 扇区 + 自旋态表 + 构型表，
                        结果写入 results/exact_{N_h}h.txt
    argv 默认取 sys.argv[1:]，方便在测试里直接传参。
    """
    args = list(sys.argv[1:] if argv is None else argv)
    out = Report()

    if any(arg in ("-h", "--help") for arg in args):
        out(__doc__)
        return
    try:
        n_holes, quick, sector, compare = parse_args(args)
    except (ValueError, IndexError) as exc:
        out(f"参数错误: {exc}\n")
        out(__doc__)
        return

    model = load_model()
    if quick:
        n_up, n_dn = sector or default_sector(n_holes)
        quick_config_report(out, model, n_up, n_dn, n_up + n_dn)
        return

    if compare:
        out("=" * 78)
        out("La3Ni2O7 双层模型精确对角化   电荷序能量比较")
        out("=" * 78)
        report_charge_order(out, model)
        return

    report_header(out, model, n_holes)

    t0 = time.time()
    spectra = solve_spectra(model, out, n_holes)
    out(f"\n求解耗时 {time.time() - t0:.1f} s")

    ground = min(spectra, key=lambda sector: sector.evals[0])
    E0 = float(ground.evals[0])

    report_ground(out, ground, E0, ground.degeneracy)
    report_spin_states(out, spectra, model, n_holes)

    groups = mode_groups()
    ordered, basis_count = type_weights(model, ground.masks, ground.weight, groups)
    report_configurations(out, ordered, basis_count, block_modes(groups), n_holes)

    out(f"\n总计耗时 {time.time() - t0:.1f} s")
    out.save(out_file(n_holes))


if __name__ == "__main__":
    main()
