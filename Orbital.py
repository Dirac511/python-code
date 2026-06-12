import numpy as np
from mayavi import mlab

def dx2y2_orbital(x, y, z, center=(0, 0, 0)):
    x_shift = x - center[0]
    y_shift = y - center[1]
    z_shift = z - center[2]
    r = 1.2 * np.sqrt(0.3 * x_shift ** 2 + 0.3 * y_shift ** 2 + 0.7 * z_shift ** 2)
    r_safe = np.where(r == 0, 1e-10, r)
    return ((2 * x_shift ** 2 - 2 * y_shift ** 2) / r_safe ** 2) * np.exp(-r / 1)


def dz2_orbital(x, y, z, center=(0, 0, 0)):
    x_shift = x - center[0]
    y_shift = y - center[1]
    z_shift = z - center[2]
    r = 1.1 * np.sqrt(x_shift ** 2 + y_shift ** 2 + 0.3 * z_shift ** 2)
    r_safe = np.where(r == 0, 1e-10, r)
    return ((3 * z_shift ** 2 - r ** 2) / r_safe ** 2) * np.exp(-r / 1)

def px_orbital(x, y, z, center=(0, 0, 0)):
    x_shift = x - center[0]
    y_shift = y - center[1]
    z_shift = z - center[2]
    r = 1.1 * np.sqrt(x_shift**2 + y_shift**2 + z_shift**2)
    r_safe = np.where(r == 0, 1e-10, r)
    return (1.5*x_shift / r_safe) * np.exp(-r / 1.5)


def py_orbital(x, y, z, center=(0, 0, 0)):
    x_shift = x - center[0]
    y_shift = y - center[1]
    z_shift = z - center[2]
    r = 1.1 * np.sqrt(x_shift**2 + y_shift**2 + z_shift**2)
    r_safe = np.where(r == 0, 1e-10, r)
    return (1.5*y_shift / r_safe) * np.exp(-r / 1.5)


def pz_orbital(x, y, z, center=(0, 0, 0)):
    x_shift = x - center[0]
    y_shift = y - center[1]
    z_shift = z - center[2]
    r = 1.1 * np.sqrt(x_shift**2 + y_shift**2 + z_shift**2)
    r_safe = np.where(r == 0, 1e-10, r)
    return (1.5*z_shift / r_safe) * np.exp(-r / 1.5)


def plot_dx2y2_orbitals():
    x = np.linspace(-10, 10, 100)
    y = np.linspace(-10, 10, 100)
    z = np.linspace(-10, 10, 100)
    X, Y, Z = np.meshgrid(x, y, z, indexing='ij')
    psi = dx2y2_orbital(X, Y, Z, 1)

    mlab.figure(size=(1200, 900), bgcolor=(0.98, 0.98, 0.98))

    color_pos = (0.85, 0.45, 0.20)  # 正相位 - 暖橙色
    color_neg = (0.20, 0.60, 0.60)  # 负相位 - 青绿色

    contour_level = 0.08

    mlab.contour3d(X, Y, Z, psi, contours=[contour_level],
                   color=color_pos, opacity=0.95,
                   transparent=True, name='Positive Phase')

    mlab.contour3d(X, Y, Z, -psi, contours=[contour_level],
                   color=color_neg, opacity=0.95,
                   transparent=True, name='Negative Phase')

    mlab.show()


def plot_dz2_orbitals():
    x = np.linspace(-8, 8, 100)
    y = np.linspace(-8, 8, 100)
    z = np.linspace(-12, 12, 100)
    X, Y, Z = np.meshgrid(x, y, z, indexing='ij')
    psi = dz2_orbital(X, Y, Z)

    mlab.figure(size=(1200, 900), bgcolor=(0.98, 0.98, 0.98))

    color_pos = (0.85, 0.45, 0.20)  # 正相位 - 暖橙色
    color_neg = (0.20, 0.60, 0.60)  # 负相位 - 青绿色

    contour_level = 0.1

    mlab.contour3d(X, Y, Z, psi, contours=[contour_level],
                   color=color_pos, opacity=0.95,
                   transparent=True, name='Positive Phase')

    mlab.contour3d(X, Y, Z, -psi, contours=[contour_level],
                   color=color_neg, opacity=0.95,
                   transparent=True, name='Negative Phase')

    mlab.axes(xlabel='X', ylabel='Y', zlabel='Z',
              color=(0.5, 0.5, 0.5), nb_labels=5,
              ranges=[-8, 16, -8, 8, -8, 8])
    mlab.outline(color=(0.5, 0.5, 0.5))

    mlab.show()


def plot_orbitals():
    x = np.linspace(-25, 10, 100)
    y = np.linspace(-25, 10, 100)
    z = np.linspace(-16, 35, 100)
    X, Y, Z = np.meshgrid(x, y, z, indexing='ij')

    #grid_centers = [(-8, -8, 0),
    #                (-8, 8, 0),
    #                (8, -8, 0),
    #                (8, 8, 0),
    #                (8, 24, 0),
    #                (-8, 24, 0),
    #                (-8, -8, 20),
    #                (-8, 8, 20),
    #                (8, -8, 20),
    #                (8, 8, 20),
    #                (8, 24, 20),
    #                (-8, 24, 20)]

    grid_centerds = [(-8, -8, -4),
                     (-8, -8, 20)]

    grid_centerpxs = [(-20, -8, -4),
                      (4, -8, -4),
                      (-20, -8, 20),
                       (4, -8, 20)]

    grid_centerpys = [(-8, -20, -4),
                      (-8, 4, -4),
                      (-8, -20, 20),
                      (-8, 4, 20)
                      ]

    grid_centerpzs = [(-8, -8, 8)]

    mlab.figure(size=(1200, 900), bgcolor=(1, 1, 1))

    color_pos = (0.90, 0.62, 0.42)  # 香槟橙
    color_neg = (0.42, 0.65, 0.82)  # 静谧蓝

    contour_level = 0.1

    for i, center in enumerate(grid_centerds):
        psi = dz2_orbital(X, Y, Z, center)
        psi1 = dx2y2_orbital(X, Y, Z, center)

        # dz2 轨道
        mlab.contour3d(X, Y, Z, psi, contours=[contour_level],
                       color=color_pos, opacity=1,
                       transparent=True, name=f'dz2_pos_{i}')

        mlab.contour3d(X, Y, Z, -psi, contours=[contour_level],
                       color=color_neg, opacity=1,
                       transparent=True, name=f'dz2_neg_{i}')

        # dx2y2 轨道
        mlab.contour3d(X, Y, Z, psi1, contours=[contour_level],
                       color=color_pos, opacity=1,
                       transparent=True, name=f'dx2y2_pos_{i}')

        mlab.contour3d(X, Y, Z, -psi1, contours=[contour_level],
                       color=color_neg, opacity=1,
                       transparent=True, name=f'dx2y2_neg_{i}')

    for i, center in enumerate(grid_centerpxs):
        psi = px_orbital(X, Y, Z, center)

        # px 轨道
        mlab.contour3d(X, Y, Z, psi, contours=[contour_level],
                       color=color_pos, opacity=1,
                       transparent=True, name=f'px_pos_{i}')

        mlab.contour3d(X, Y, Z, -psi, contours=[contour_level],
                       color=color_neg, opacity=1,
                       transparent=True, name=f'px_neg_{i}')

    for i, center in enumerate(grid_centerpys):
        psi = py_orbital(X, Y, Z, center)

        # py 轨道
        mlab.contour3d(X, Y, Z, psi, contours=[contour_level],
                       color=color_pos, opacity=1,
                       transparent=True, name=f'py_pos_{i}')

        mlab.contour3d(X, Y, Z, -psi, contours=[contour_level],
                       color=color_neg, opacity=1,
                       transparent=True, name=f'py_neg_{i}')

    for i, center in enumerate(grid_centerpzs):
        psi = pz_orbital(X, Y, Z, center)

        # pz 轨道
        mlab.contour3d(X, Y, Z, psi, contours=[contour_level],
                       color=color_pos, opacity=1,
                       transparent=True, name=f'pz_pos_{i}')

        mlab.contour3d(X, Y, Z, -psi, contours=[contour_level],
                       color=color_neg, opacity=1,
                       transparent=True, name=f'pz_neg_{i}')

    # 设置视角
    mlab.view(azimuth=30, elevation=60, distance=100, focalpoint=(10, 10, 10))

    mlab.savefig('nature_style_orbitals.png', size=(2000, 2000), magnification=4)
    mlab.show()


if __name__ == "__main__":
    # plot_dx2y2_orbitals()
    # plot_dz2_orbitals()

    plot_orbitals()