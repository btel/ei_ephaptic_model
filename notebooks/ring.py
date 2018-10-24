from matplotlib import patches
import matplotlib.pyplot as plt
import numpy as np
n_nodes = 7

# %%
fig = plt.figure(figsize=(3,3))
ax = plt.subplot(111)
ax.axis('equal')
ax.axis('off')
big_radius = 20
delta = np.pi * 2 / n_nodes
small_radius = np.sin(delta / 2) * big_radius
for i in range(n_nodes):
    x, y = np.sin(delta * i) * big_radius, np.cos(delta * i) * big_radius
    circle = patches.Circle((x, y), radius=small_radius, edgecolor='k',
                            facecolor='none')
    ax.text(x, y, i, va='center', ha='center', size=15)
    ax.add_patch(circle)
ax.autoscale_view()
fig
fig.savefig('ring.svg', transparent=True)
