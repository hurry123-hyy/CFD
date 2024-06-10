"""
@Author: Hu yiyue
@Detail: COVID19 streamLine.py
@Time: 2024-06-04 11:18
@E-mail: yiyuehuu@gmail.com
@Description: 画流场图
"""

import numpy as np
import matplotlib.pyplot as plt

# 打开文件，读取所有行
with open('./dataset/flow.txt', 'r') as f:
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

# 转换为浮点数数组
data = np.array(data, dtype=float)

x = data[:, 0]
y = data[:, 1]
m = data[:, 6]

# 将数据整理成网格
X = x.reshape((201, 401))  # 根据你的数据维度
Y = y.reshape((201, 401))
M = m.reshape((201, 401))

# 绘制流场图
plt.figure(figsize=(8, 6))
contour = plt.contourf(X, Y, M, levels=100, cmap='jet')
plt.colorbar(contour)
plt.title('Roe (MUSCL)')
plt.xlabel('x')
plt.ylabel('y')

# 保存图像
plt.savefig('./dataset/flow_field.png')

plt.show()
