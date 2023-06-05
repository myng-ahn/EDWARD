import sys
import numpy as np
import matplotlib.pyplot as plt

USAGE = 'python plot_plink.py <plink_eigenvec> <PC_x> <PC_y>'
if len(sys.argv) != 4:
    print(USAGE)
    exit(1)

__, fname, pc_x, pc_y = sys.argv
pc_x = int(pc_x)
pc_y = int(pc_y)
pca_transformed = np.loadtxt(fname, dtype='str', delimiter=' ')[:,2:]

if (pc_x > pca_transformed.shape[1]) or (pc_y > pca_transformed.shape[1]):
    raise ValueError("PC exceeds number of dimensions in output")
pca_transformed = pca_transformed.astype(dtype='float64')

fig, ax = plt.subplots(figsize=(16,9))
ax.scatter(pca_transformed[:,pc_x-1],
           pca_transformed[:,pc_y-1],
           s=2, alpha=0.4)
ax.set_xlabel(f'PC {pc_x}')
ax.set_ylabel(f'PC {pc_y}')
fig.savefig(f'{fname}_pc{pc_x}_vs_pc{pc_y}.png')