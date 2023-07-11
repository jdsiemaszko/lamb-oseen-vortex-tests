import os
import sys
import yaml
import pandas
import numpy as np
import csv
import re
import timeit
from pathlib import Path
import blobs_solver2 as pHyFlow
from blobs_solver2.blobs.base.induced import vorticity
import matplotlib.pyplot as plt

case_dir = os.getcwd()
data_dir = os.path.join(case_dir,'data_avrm_compression_7')
ref_dir = os.path.join(case_dir,'data_avrm')
plots_dir =  os.path.join(case_dir, 'comparison_WORSE')

Path(plots_dir).mkdir(parents=True, exist_ok=True)

case = "gridfill_avrm_rk4"
case_ref = "avrm_rk4"
gammaC = 0.
nTimeSteps = 1000
writeInterval_plots = 100
coreSize = 'variable'
deltaTc = 0.02

# uxNorm = np.array([])
# uyNorm = np.array([])
# omegaNorm = np.array([])
# t_norm = np.array([])
#Line plots
times_file = os.path.join(data_dir,"times_{}.csv".format(case))
times_ref_file = os.path.join(ref_dir,"times_{}.csv".format(case_ref))
times_data = pandas.read_csv(times_file)
times_ref_data = pandas.read_csv(times_ref_file)

time1 = times_data['Time']
noBlobs1 = times_data['NoBlobs']
evolution_time1 = times_data['Evolution_time']
circulation1 = times_data['Circulation']

time2 = times_ref_data['Time']
noBlobs2 = times_ref_data['NoBlobs']
evolution_time2 = times_ref_data['Evolution_time']
circulation2 = times_ref_data['Circulation']

fig, ax = plt.subplots(1,1,figsize=(6,6))
ax.plot(time1,noBlobs1, label='Run w/ Compression')
ax.plot(time2,noBlobs2, label='Run w/o Compression')

plt.grid(color = '#666666', which='major', linestyle = '--', linewidth = 0.5)
plt.minorticks_on()
plt.grid(visible=True, which='minor', color='#999999', linestyle='-', alpha=0.2)
plt.title('Total number of particles')
plt.ylabel('Particles')
plt.xlabel('time $(sec)$')
plt.legend()
plt.savefig("{}/number_of_particles_{}.png".format(plots_dir,case), dpi=300, bbox_inches="tight")

fig, ax = plt.subplots(1,1,figsize=(6,6))
ax.plot(time1,circulation1- gammaC, label='Run w/ Compression')
ax.plot(time2,circulation2- gammaC, label='Run w/o Compression')

plt.grid(color = '#666666', which='major', linestyle = '--', linewidth = 0.5)
plt.minorticks_on()
plt.grid(visible=True, which='minor', color='#999999', linestyle='-', alpha=0.2)
plt.title('Circulation')
plt.ylabel('circulation')
plt.xlabel('time $(sec)$')
plt.legend()
plt.savefig("{}/circulation_{}.png".format(plots_dir,case), dpi=300, bbox_inches="tight")

fig = plt.subplots(figsize=(6,6))
index = np.arange(len(evolution_time1))
width = deltaTc
lagrangian = plt.bar(index[1:nTimeSteps]*deltaTc, evolution_time1[1:nTimeSteps] - evolution_time2[1:nTimeSteps], width)
plt.ylabel('Time (s)')
plt.xlabel('Simulation time(s)')
plt.title('Evolution Time Difference')
plt.savefig("{}/times_dif_{}.png".format(plots_dir,case), dpi=300, bbox_inches="tight")

