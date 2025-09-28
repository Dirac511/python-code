import numpy as np
from mayavi import mlab


def dx2y2_orbital(x, y, z, center=(0,0,0)):
    """计算d_x²-y²轨道波函数"""
    x_shift = x - center[0]
    y_shift = y - center[1]
    z_shift = z - center[2]
    r = np.sqrt(0.3*x_shift ** 2 + 0.3*y_shift ** 2 + z_shift ** 2)
    r_safe = np.where(r == 0, 1e-10, r)
    return ((1.4*x_shift ** 2 - 1.4*y_shift ** 2) / r_safe ** 2) * np.exp(-r / 1)  # 调整衰减参数

def dz2_orbital(x, y, z, center=(0,0,0)):
    """计算d_z²轨道波函数"""
    x_shift = x - center[0]
    y_shift = y - center[1]
    z_shift = z - center[2]
    r = np.sqrt(x_shift**2 + y_shift**2 + 0.3*z_shift**2)
    r_safe = np.where(r == 0, 1e-10, r)
    return ((3*z_shift**2 - r**2) / r_safe**2) * np.exp(-r/1)

def plot_dx2y2_orbitals():
    """绘制d_x²-y²轨道"""
    # 扩大网格范围
    x = np.linspace(-10, 10, 100)
    y = np.linspace(-10, 10, 100)
    z = np.linspace(-10, 10, 100)
    X, Y, Z = np.meshgrid(x, y, z, indexing='ij')
    # 计算波函数
    psi = dx2y2_orbital(X, Y, Z, 1)



    # 创建图形窗口
    mlab.figure(size=(1200, 900), bgcolor=(1, 1, 1))

    # 使用更低的等值面级别来显示完整轨道
    contour_level = 0.02  # 降低等值面级别


    # 绘制正相位（沿x轴方向，红色）
    mlab.contour3d(X, Y, Z, psi, contours=[contour_level],
                   color=(1, 0, 0), opacity=1.0,
                   transparent=True, name='Positive Phase')

    # 绘制负相位（沿y轴方向，蓝色）
    mlab.contour3d(X, Y, Z, -psi, contours=[contour_level],
                   color=(0, 0, 1), opacity=1.0,
                   transparent=True, name='Negative Phase')

    # 添加坐标轴
    #mlab.axes(xlabel='X', ylabel='Y', zlabel='Z',
    #          color=(0, 0, 0), nb_labels=5,
    #          ranges=[-8, 8, -8, 8, -8, 8])  # 设置坐标轴范围
    #mlab.outline(color=(0, 0, 0))

    # 添加标题
    #mlab.title(r'$d_{x^2-y^2}$ 电子轨道 - 完整显示', size=0.5, height=0.95, color=(0, 0, 0))

    # 设置视角
    #mlab.view(azimuth=45, elevation=60, distance=30)

    # 显示图形
    mlab.show()

def plot_dz2_orbitals():
    """绘制d_x²-y²轨道"""
    # 扩大网格范围
    x = np.linspace(-8, 8, 100)
    y = np.linspace(-8, 8, 100)
    z = np.linspace(-12, 12, 100)
    X, Y, Z = np.meshgrid(x, y, z, indexing='ij')
    # 计算波函数值
    psi = dz2_orbital(X, Y, Z)
    # 创建图形窗口
    mlab.figure(size=(1200, 900), bgcolor=(1, 1, 1))
    # 使用更低的等值面级别来显示完整轨道
    contour_level = 0.02  # 降低等值面级别
    # 绘制正相位（沿x轴方向，红色）
    mlab.contour3d(X, Y, Z, psi, contours=[contour_level],
                   color=(1, 0, 0), opacity=1.0,
                   transparent=True, name='Positive Phase')
    # 绘制负相位（沿y轴方向，蓝色）
    mlab.contour3d(X, Y, Z, -psi, contours=[contour_level],
                   color=(0, 0, 1), opacity=1.0,
                   transparent=True, name='Negative Phase')

    # 添加坐标轴
    mlab.axes(xlabel='X', ylabel='Y', zlabel='Z',
               color=(1, 1, 0), nb_labels=5,
               ranges=[-8, 8, -8, 8, -8, 8])  # 设置坐标轴范围
    mlab.outline(color=(1, 1, 0))

    # 添加标题
    # mlab.title(r'$d_{x^2-y^2}$ 电子轨道 - 完整显示', size=0.5, height=0.95, color=(0, 0, 0))

    # 设置视角
    # mlab.view(azimuth=45, elevation=60, distance=30)

    # 显示图形
    mlab.show()

def plot_orbitals():
    """绘制d_x²-y²轨道"""
    # 扩大网格范围
    x = np.linspace(-20, 20, 100)
    y = np.linspace(-20, 20, 100)
    z = np.linspace(-30, 30, 100)
    X, Y, Z = np.meshgrid(x, y, z, indexing='ij')
    # 计算波函数值
    grid_centers = [(-8, -8, 0),
                    (-8, 8, 0),
                    (8, -8, 0),
                    (8, 8, 0),
                    (-8, -8, 20),
                    (-8, 8, 20),
                    (8, -8, 20),
                    (8, 8, 20)
                    ]

    # 创建图形窗口
    mlab.figure(size=(1200, 900), bgcolor=(1, 1, 1))
    # 使用更低的等值面级别来显示完整轨道
    contour_level = 0.2  # 降低等值面级别
    # 绘制正相位（沿x轴方向，红色）
    for i, center in enumerate(grid_centers):
        psi = dz2_orbital(X, Y, Z,center)
        psi1 = dx2y2_orbital(X, Y, Z,center)
        mlab.contour3d(X, Y, Z, psi, contours=[contour_level],
                       color=(1, 0, 0), opacity=1.0,
                       transparent=True, name='Positive Phase')
        # 绘制负相位（沿y轴方向，蓝色）
        mlab.contour3d(X, Y, Z, -psi, contours=[contour_level],
                       color=(0, 0, 1), opacity=1.0,
                       transparent=True, name='Negative Phase')

        mlab.contour3d(X, Y, Z, psi1, contours=[contour_level],
                       color=(1, 0, 0), opacity=1.0,
                       transparent=True, name='Positive Phase')
        # 绘制负相位（沿y轴方向，蓝色）
        mlab.contour3d(X, Y, Z, -psi1, contours=[contour_level],
                       color=(0, 0, 1), opacity=1.0,
                       transparent=True, name='Negative Phase')
        # 添加坐标轴
    #mlab.axes(xlabel='X', ylabel='Y', zlabel='Z',
    #         color=(0,0,0), nb_labels=5,
    #         ranges=[-50, 100, -50, 100, -50, 100])  # 设置坐标轴范围
    #mlab.outline(color=(0, 0, 0))
    #mlab.plot3d([-20, 20], [20, -20], [0, 0], color=(0.7, 0.7, 0.7), tube_radius=0.01)
    #mlab.plot3d([-20, 20], [20, -20], [20, 20], color=(0.7, 0.7, 0.7), tube_radius=0.01)
    plane_size = 20
    xx, yy = np.mgrid[-plane_size:plane_size:50j, -plane_size:plane_size:50j]
    zz = np.zeros_like(xx)
    zz1 = np.ones_like(xx)*20
    mlab.mesh(xx, yy, zz, color=(0.7, 0.9, 0.9), opacity=0.4)
    mlab.mesh(xx, yy, zz1, color=(0.7, 0.9, 0.9), opacity=0.4)
    # 添加标题
    # mlab.title(r'$d_{x^2-y^2}$ 电子轨道 - 完整显示', size=0.5, height=0.95, color=(0, 0, 0))

    # 设
    mlab.view(azimuth=30, elevation=70, distance=100)

    # 显示图形
    mlab.savefig('1.png', size=(2000, 2000),magnification=4)
    mlab.show()

if __name__ == "__main__":
    #plot_dx2y2_orbitals()
    #plot_dz2_orbitals()
    plot_orbitals()
    plot_orbitals()

