# %%
import matplotlib.pyplot as plt
import math
from matplotlib.patches import RegularPolygon
from matplotlib.collections import RegularPolyCollection
from matplotlib import cm
from mpl_toolkits.axes_grid1 import make_axes_locatable

# %%
def plot_map(grid,
             d_matrix):
    """
    Plot hexagon map where each neuron is represented by a hexagon.
    """
    n_centers = grid['centers']
    r = grid['r']
    ax = plt.gca()
    ax.axis('equal')
    # Discover difference between centers
    for xy in n_centers:
        patch = RegularPolygon(xy, numVertices=6, radius=r, color='none',
                               ec='k')
        ax.add_patch(patch)
        patch.set_transform(ax.transData)
    ax.autoscale_view()
    plt.axis('off')

    return ax

# %%
import numpy as np
r = 20.
dx = (math.cos(math.pi / 6)) * r
a = math.sin(math.pi/6) * 2 * r
dy = r + a / 2.
centers = [(-dx, 0),
           (0, dy),
           (2 * dx, dy),
           (3 * dx, 0),
           (2 * dx, -dy),
           (0, -dy),
           (dx, 0)
             ]
x0 = np.array([-dx, 0, 0, dx, 2 * dx, 2 * dx, 3 * dx])
y0 = np.array([dy, 0, 2 * dy, dy, 2 * dy, 0, dy])
#centers = np.vstack((x0, y0)).T
grid = {'centers': centers,
        'r' : r}


# %%
fig = plt.figure()
ax = fig.add_subplot(111)
plot_map(grid, np.random.rand(len(grid['centers'])))
for i, xy in enumerate(centers):
    plt.text(*xy, str(i), ha='center', va='center',
             size=15)
plt.savefig('hexagonmap.svg')
