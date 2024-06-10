
"""
@Author: Hu yiyue
@Detail: COVID19 drawnInitValue.py
@Time: 2024-06-04 17:28
@E-mail: yiyuehuu@gmail.com
@Description: 绘制网格示意图
"""

import matplotlib.pyplot as plt

# Create the figure and axis
fig, ax = plt.subplots(figsize=(14, 7))

# Define the grid dimensions and rectangle parameters
ist, jst = 0, 0
ied, jed = 22, 11
lw, jw1, jw2 = 3, 4, 6

# Draw the grid points and lines
for x in range(ied + 1):
    for y in range(jed + 1):
        if x % 2 == 0:
            ax.plot(x, y, 'bs', markersize=6)  # Blue squares
        else:
            ax.plot(x, y, 'ro', markersize=4)  # Red circles

        if x < ied:
            ax.plot([x, x + 1], [y, y], 'k-', linewidth=0.5)
        if y < jed:
            ax.plot([x, x], [y, y + 1], 'k-', linewidth=0.5)

# Draw the inner rectangle
rect = plt.Rectangle((6, 2), 10, 7, linewidth=2, edgecolor='k', facecolor='none')
ax.add_patch(rect)

# Annotate the points and make the corresponding points larger
annotations = {
    "(ist, jst)": (ist+4, jst),
    "(ist, jed)": (ist+4, jed),
    "(ied, jst)": (ied-4, jst),
    "(ied, jed)": (ied-4, jed),
    "(ist+lw, jst+jw1)": (6, 2),
    "(ist+lw, jst+jw2)": (6, 9),
    # "(ied-1, jst)": (ied - 5, jst),
    # "(ied-1, jed)": (ied - 5, jed)
}

for label, (x, y) in annotations.items():
    ax.plot(x, y, 'bo', markersize=10)  # Make the annotated points larger
    ax.text(x, y + 0.2, label, fontsize=12, ha='right')  # Move the text slightly up

ax.set_aspect('equal')
ax.set_xticks([])
ax.set_yticks([])
ax.set_xticklabels([])
ax.set_yticklabels([])
plt.grid(False)

plt.savefig('./dataset/streamplot_category.png', dpi=300)

plt.show()
