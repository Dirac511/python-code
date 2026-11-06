"""2x2 格点图（双层）：上层两个格点放 d(x^2-y^2) 与 d(z^2)，下层同样。

虚线沿格点连线：每层内连接两个格点，并在上下层对应格点之间竖直相连。
着色用双光源 Lambert + 弱镜面（按解析法向计算），不做投影阴影。
用法：py plot_eg_orbitals.py
输出：results/eg_orbitals_2x2.png
"""

import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

plt.rcParams["axes.unicode_minus"] = False

COL_POS = np.array([0.28, 0.55, 0.84])     # 相位 +（明亮钢蓝）
COL_NEG = np.array([0.88, 0.45, 0.40])     # 相位 −（明亮砖红）
N_GRID = 241
FATTEN = 0.68          # <1：把 |Y| 的细长瓣压得饱满一点

A = 2.30               # 层内半格距（沿 x）
L = 1.85               # 层间距的一半（沿 z）
R_ORB = 1.00           # 轨道半径

ELEV, AZIM = 18.0, -90.0
# 主光 + 补光
KEY_AZ, KEY_ALT, KEY_W = 215.0, 55.0, 0.85
FILL_AZ, FILL_ALT, FILL_W = 25.0, 12.0, 0.30
AMBIENT = 0.52


def orbital_surface(kind, n=N_GRID):
    """实球谐函数 |Y| 的等值面。kind: 'x2-y2' | 'z2'。"""
    theta = np.linspace(0.0, np.pi, n)
    phi = np.linspace(0.0, 2.0 * np.pi, n)
    th, ph = np.meshgrid(theta, phi)

    if kind == "x2-y2":
        ang = np.cos(2.0 * ph)                     # 四瓣，沿 ±x ±y
        rad = np.abs(ang) * np.sin(th) ** 2
    else:
        ang = 3.0 * np.cos(th) ** 2 - 1.0          # 两瓣沿 ±z + 腰环
        rad = np.abs(ang)

    rad = R_ORB * (rad / rad.max()) ** FATTEN
    x = rad * np.sin(th) * np.cos(ph)
    y = rad * np.sin(th) * np.sin(ph)
    z = rad * np.cos(th)
    return x, y, z, ang


def _unit(v):
    return v / (np.linalg.norm(v) + 1e-12)


def _normals(x, y, z):
    """参数曲面的单位外法向。"""
    d0 = np.stack([np.gradient(x, axis=0), np.gradient(y, axis=0),
                   np.gradient(z, axis=0)], axis=-1)
    d1 = np.stack([np.gradient(x, axis=1), np.gradient(y, axis=1),
                   np.gradient(z, axis=1)], axis=-1)
    n = np.cross(d0, d1)
    n /= np.linalg.norm(n, axis=-1, keepdims=True) + 1e-12

    outward = np.stack([x, y, z], axis=-1)         # 星形曲面：与径矢同向即朝外
    flip = np.where((n * outward).sum(-1) < 0.0, -1.0, 1.0)
    return n * flip[..., None]


def _light(az, alt):
    az, alt = np.deg2rad(az), np.deg2rad(alt)
    return np.array([np.cos(alt) * np.cos(az),
                     np.cos(alt) * np.sin(az),
                     np.sin(alt)])


def shade(x, y, z, rgb):
    """双光源 Lambert 漫反射 + 弱而宽的镜面高光 + 环境光。"""
    n = _normals(x, y, z)
    key = _light(KEY_AZ, KEY_ALT)
    fill = _light(FILL_AZ, FILL_ALT)
    view = _light(AZIM, ELEV)
    half = _unit(key + view)

    diff = (KEY_W * np.clip(n @ key, 0, 1)
            + FILL_W * np.clip(n @ fill, 0, 1))
    spec = 0.13 * np.clip(n @ half, 0, 1) ** 36
    rim = 0.16 * (1.0 - np.clip(np.abs(n @ view), 0, 1)) ** 3   # 边缘光

    shade = AMBIENT + 0.62 * diff
    out = rgb * shade[..., None] + spec[..., None] + rim[..., None]
    rgba = np.empty(x.shape + (4,))
    rgba[..., :3] = np.clip(out, 0.0, 1.0)
    rgba[..., 3] = 1.0
    return rgba


def draw_orbital(ax, kind, cx, cz):
    x, y, z, ang = orbital_surface(kind)

    rgb = np.empty(x.shape + (3,))
    pos = ang >= 0.0
    rgb[pos] = COL_POS
    rgb[~pos] = COL_NEG

    ax.plot_surface(
        x + cx, y, z + cz, facecolors=shade(x, y, z, rgb),
        rstride=1, cstride=1, linewidth=0, antialiased=True,
        shade=False, zorder=4,
    )


def draw_layers(ax):
    """两层各一张极淡的平面，标出上/下层高度。"""
    xx = np.array([[-A - 1.0, A + 1.0], [-A - 1.0, A + 1.0]])
    yy = np.array([[-1.55, -1.55], [1.55, 1.55]])
    for cz in (-L, L):
        ax.plot_surface(
            xx, yy, np.full_like(xx, cz), color="0.55", alpha=0.10,
            linewidth=0, antialiased=True, shade=False, zorder=1,
        )


def draw_bonds(ax):
    """层内相邻格点 + 上下层对应格点之间的虚线。"""
    kw = dict(ls=(0, (6, 4)), lw=1.4, color="0.35", zorder=3)
    for cz in (-L, L):                    # 每层内沿 x 连接两个格点
        ax.plot([-A, A], [0, 0], [cz, cz], **kw)
    for cx in (-A, A):                    # 上下层之间竖直相连
        ax.plot([cx, cx], [0, 0], [-L, L], **kw)


def main():
    fig = plt.figure(figsize=(8.6, 7.4))
    ax = fig.add_axes([0.0, 0.0, 1.0, 1.0], projection="3d")
    ax.computed_zorder = False          # 层平面 < 虚线 < 轨道

    draw_layers(ax)
    draw_bonds(ax)
    for cz in (-L, L):                  # 下层、上层；每层左 x2-y2、右 z2
        for cx, kind in ((-A, "x2-y2"), (A, "z2")):
            draw_orbital(ax, kind, cx, cz)

    pad_x, pad_z = 0.85, 0.55
    ax.set_xlim(-A - pad_x, A + pad_x)
    ax.set_ylim(-R_ORB - 0.6, R_ORB + 0.6)
    ax.set_zlim(-L - R_ORB - pad_z, L + R_ORB + pad_z)
    ax.set_box_aspect(
        (2 * (A + pad_x), 2 * (R_ORB + 0.6), 2 * (L + R_ORB + pad_z)),
        zoom=1.30,
    )
    ax.set_proj_type("persp", focal_length=0.85)
    ax.view_init(elev=ELEV, azim=AZIM)
    ax.set_axis_off()

    os.makedirs("results", exist_ok=True)
    out = os.path.join("results", "eg_orbitals_2x2.png")
    fig.savefig(out, dpi=300)
    print("saved:", os.path.abspath(out))


if __name__ == "__main__":
    main()
