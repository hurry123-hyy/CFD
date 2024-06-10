"""
@Author: Hu yiyue
@Detail: COVID19 drawnInitValue.py
@Time: 2024-06-04 21:45
@E-mail: yiyuehuu@gmail.com
@Description: 绘制初值赋值示意图
"""
import matplotlib.pyplot as plt
from matplotlib import rcParams
#
# # 设置字体以显示中文
# rcParams['font.sans-serif'] = ['SimHei']  # 使用黑体
# rcParams['axes.unicode_minus'] = False  # 解决负号显示问题
#
# # 创建图形和轴
# fig, ax = plt.subplots(figsize=(8, 6))
#

# 创建图形和轴,2个子图
fig, ax = plt.subplots(1, 2, figsize=(16, 12))

# 定义网格维度
grid_lines = 6
grid_columns = 6

# 定义红色和蓝色点的位置和标签
red_points = [(1, 1), (1, 2), (2, 3), (3, 4), (4, 4), (4, 4)]
blue_points = [(3, 1), (4, 1), (3, 2), (4, 2)]
red_labels = ['②', '④', '①', '④', '③']
blue_labels = ['②', '①', '④', '③']

# 画网格线
for i in range(grid_lines):
    for j in range(grid_columns):
        ax[0].plot([i, i], [0, grid_columns - 1], 'k-', linewidth=1)
        ax[0].plot([0, grid_lines - 1], [j, j], 'k-', linewidth=1)

# 画红色点和标签
for (x, y), label in zip(red_points, red_labels):
    ax[0].plot(x, y, 'ro', markersize=6)
    ax[0].text(x, y, label, fontsize=18, ha='right', va='bottom', color='black')

# 画蓝色点和标签
for (x, y), label in zip(blue_points, blue_labels):
    ax[0].plot(x, y, 'bs', markersize=6)
    ax[0].text(x, y, label, fontsize=18, ha='right', va='bottom', color='black')

# 画红色粗线
ax[0].plot([2, 2], [0, 3], 'r-', linewidth=3)
ax[0].plot([2, 5], [3, 3], 'r-', linewidth=3)

# 添加文本注释
annotations = {
    (3, 2): '算术\n平均',
    (4, 2): '镜面\n反射',
    (3, 1): '镜面\n反射',
    (4, 1): '矢量取反\n标量取值',
}

for (x, y), text in annotations.items():
    ax[0].text(x, y, text, fontsize=20, ha='left', va='top', color='black')

ax[0].set_aspect('equal')
ax[0].set_xticks([])
ax[0].set_yticks([])
ax[0].set_xticklabels([])
ax[0].set_yticklabels([])
plt.grid(False)
#
# plt.savefig('./dataset/initValue.png', dpi=300)
#
# plt.show()

# 画左下角点

# 设置字体以显示中文
rcParams['font.sans-serif'] = ['SimHei']  # 使用黑体
rcParams['axes.unicode_minus'] = False  # 解决负号显示问题



# 定义网格维度
grid_lines = 6
grid_columns = 6

# 定义红色和蓝色点的位置和标签
red_points = [(1, 3), (1, 2), (2, 1), (3, 0), (4, 0), (4, 4)]
blue_points = [(3, 2), (4, 2), (3, 3), (4, 3)]
red_labels = ['②', '④', '①', '④', '③']
blue_labels = ['④', '③', '②', '①']

# 画网格线
for i in range(grid_lines):
    for j in range(grid_columns):
        ax[1].plot([i, i], [0, grid_columns - 1], 'k-', linewidth=1)
        ax[1].plot([0, grid_lines - 1], [j, j], 'k-', linewidth=1)

# 画红色点和标签
for (x, y), label in zip(red_points, red_labels):
    ax[1].plot(x, y, 'ro', markersize=6)
    ax[1].text(x, y, label, fontsize=18, ha='right', va='bottom', color='black')

# 画蓝色点和标签
for (x, y), label in zip(blue_points, blue_labels):
    ax[1].plot(x, y, 'bs', markersize=6)
    ax[1].text(x, y, label, fontsize=18, ha='right', va='bottom', color='black')

# 画红色粗线
ax[1].plot([2, 2], [1, 4], 'r-', linewidth=3)
ax[1].plot([2, 5], [1, 1], 'r-', linewidth=3)

# 添加文本注释
annotations = {
    (3, 2): '算术\n平均',
    (4, 2): '镜面\n反射',
    (3, 3): '镜面\n反射',
    (4, 3): '矢量取反\n标量取值',
}

for (x, y), text in annotations.items():
    ax[1].text(x, y, text, fontsize=20, ha='left', va='top', color='black')

ax[1].set_aspect('equal')
ax[1].set_xticks([])
ax[1].set_yticks([])
ax[1].set_xticklabels([])
ax[1].set_yticklabels([])
plt.grid(False)
plt.subplots_adjust(hspace=0.05, wspace=0.05)

plt.savefig('./dataset/initValue.png', dpi=300)

plt.show()
