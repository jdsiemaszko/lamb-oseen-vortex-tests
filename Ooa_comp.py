from blobs_solver2.blobs.base.merging import remesh
import numpy as np 
import timeit
import matplotlib.pyplot as plt
from pathlib import Path
import blobs_solver2 as pHyFlow
# import pHyFlow
# import VorticityFoamPy as foam
import solvers.particle as particle
from blobs_solver2.blobs.base.induced import vorticity
from matplotlib.lines import Line2D

import os
import sys
import yaml
import re
import csv
import pandas

case = 'remesh'
case_dir = os.getcwd()

methods = ['induced', 'M4prime', 'Stock', 'velocity', 'multigrid_medTH', 'multigrid_lowTH', 'multigrid_highTH']
subfolders = ['order_of_accuracy_' + method for method in methods]
source_dirs = [os.path.join(case_dir, dir) for dir in subfolders]

dest_dir = os.path.join(case_dir, 'OoA_Comparison')
Path(dest_dir).mkdir(parents=True, exist_ok=True)

fig1, ax1 = plt.subplots(1, 1)
fig2, ax2 = plt.subplots(1, 1)
fig3, ax3 = plt.subplots(1, 1)

iterator = iter(['*', 'x', 'o', 'p', 'v', 's', '8'])
for subfolder, source_dir in zip(subfolders,source_dirs):
    arr = np.loadtxt(os.path.join(source_dir, "results_.csv"), delimiter = ' ').T
    compression_list = arr[-5, :]
    h_list = arr[-4, :]
    L2_list = arr[-3, :]
    Linf_list = arr[-2, :]
    circ_list = np.abs(arr[-1, :])

    i = next(iterator)
    ax1.plot(h_list, L2_list, label = subfolder[18:], marker=i)
    ax2.plot(h_list, Linf_list, label= subfolder[18:], marker=i)
    ax3.plot(h_list, circ_list, label= subfolder[18:], marker=i)
    




ticks = [1, ]
ax1.set_xscale('log')
ax1.set_yscale('log')
ax1.set_xlabel('normalized spacing: 1/sqrt(CR)')
ax1.set_ylabel('L2 Error')
ax1.grid(which='both')
ax1.legend()
fig1.savefig("{}/L2_comparison_{}.png".format(dest_dir,case), dpi=300, bbox_inches="tight")
plt.close(fig1)

ax2.set_xscale('log')
ax2.set_yscale('log')
ax2.set_xlabel('normalized spacing: 1/sqrt(CR)')
ax2.set_ylabel('Linf Error')
ax2.grid(which='both')
ax2.legend()
fig2.savefig("{}/Linf_comparison{}.png".format(dest_dir,case), dpi=300, bbox_inches="tight")
plt.close(fig2)

ax3.set_xscale('log')
ax3.set_yscale('log')
ax3.set_xlabel('normalized spacing: 1/sqrt(CR)')
ax3.set_ylabel('Circulation Error')
ax3.grid(which='both')
ax3.legend()
fig3.savefig("{}/circ_comparison{}.png".format(dest_dir,case), dpi=300, bbox_inches="tight")
plt.close(fig3)