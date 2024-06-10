"""
@Author: Hu yiyue
@Detail: COVID19 streamLine_category.py
@Time: 2024-06-04 11:54
@E-mail: yiyuehuu@gmail.com
@Description: 按照不同的物理量分开画图
"""

import numpy as np
import matplotlib.pyplot as plt

# 打开文件，读取所有行
with open('./dataset/101-51.txt', 'r') as f:
    lines = f.readlines()

# 找到数据开始的行数（跳过头信息）
start_line = 0
for i, line in enumerate(lines):
    if "DATAPACKING=POINT" in line:
        start_line = i + 1
        break

# 读取数据，跳过头信息行
data_lines = lines[start_line:]

# 去除每行的换行符，并将其拆分成列表
data = []
for line in data_lines:
    parts = line.strip().split()
    if len(parts) == 8:  # 确保每行有8个数据点
        data.append(parts)
    else:
        print(line)

# 转换为浮点数数组
data = np.array(data, dtype=float)
print(data.shape)

x = data[:, 0]
y = data[:, 1]
rho = data[:, 2]
u = data[:, 3]
v = data[:, 4]
p = data[:, 5]
m = data[:, 6]

# 计算速度大小
speed = np.sqrt(u**2 + v**2)

# 将数据整理成网格
X = x.reshape((51, 101))  # 根据你的数据维度
Y = y.reshape((51, 101))
R = rho.reshape((51, 101))
U = u.reshape((51, 101))
V = v.reshape((51, 101))
P = p.reshape((51, 101))
M = m.reshape((51, 101))
Speed = speed.reshape((51, 101))

# 创建一个图形窗口和三个子图
fig, axs = plt.subplots(1, 3, figsize=(30, 5))  # 1行3列

# 绘制压力流场图
cax1 = axs[0].contourf(X, Y, P, levels=100, cmap='jet')
fig.colorbar(cax1, ax=axs[0])
# axs[0].streamplot(X, Y, U, V, color='k')  # 注意这里改为 axs[0]
axs[0].set_title('101_51_Pressure')
axs[0].set_xlabel('x')
axs[0].set_ylabel('y')

# 绘制密度流场图
cax2 = axs[1].contourf(X, Y, R, levels=100, cmap='jet')
fig.colorbar(cax2, ax=axs[1])
# axs[1].streamplot(X, Y, U, V, color='k')  # 注意这里改为 axs[1]
axs[1].set_title('101_51_Rho')
axs[1].set_xlabel('x')
axs[1].set_ylabel('y')

# 绘制速度大小流场图
cax3 = axs[2].contourf(X, Y, Speed, levels=100, cmap='jet')  # 使用 Speed
fig.colorbar(cax3, ax=axs[2])
axs[2].streamplot(X, Y, U, V, color='k')  # 注意这里改为 axs[2]
axs[2].set_title('101_51_Velocity')
axs[2].set_xlabel('x')
axs[2].set_ylabel('y')

# 隐藏所有子图的边框
for ax in axs.flat:
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.spines['bottom'].set_visible(False)

# 调整子图布局以减少空隙
plt.tight_layout()

# 如果需要更精细的控制，可以使用subplots_adjust
# hspace 和 wspace 分别控制子图间的垂直和水平间距
plt.subplots_adjust(hspace=0.05, wspace=0.05)

# 保存图像
plt.savefig('./dataset/streamplot_category.png', dpi=300)

# 显示图形
plt.show()