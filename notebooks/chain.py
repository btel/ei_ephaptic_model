from matplotlib import patches
import matplotlib.pyplot as plt
import numpy as np
n_nodes = 7

# %%
fig = plt.figure(figsize=(5,5))
ax = plt.subplot(111)
ax.axis('equal')
ax.axis('off')
small_radius = 1./(n_nodes + 1) / 2
y = 0.5
for i in range(n_nodes):
    x = small_radius + i * small_radius * 2
    circle = patches.Circle((x, y),
                            radius=small_radius,
                            edgecolor='k',
                            facecolor='none',
                            transform=ax.transAxes,
                            clip_on=False)
    ax.text(x, y, i, va='center', ha='center', size=12,
            transform=ax.transAxes)
    ax.add_patch(circle)
ax.autoscale_view()
fig.savefig('chain.svg')
