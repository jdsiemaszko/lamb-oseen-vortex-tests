from blobs_solver2.blobs.base.merging import remesh, multiscale_remesh
import numpy as np 
import timeit
import matplotlib.pyplot as plt
from pathlib import Path
import blobs_solver2 as pHyFlow
# import pHyFlow
# import VorticityFoamPy as foam
import solvers.particle as particle
from blobs_solver2.blobs.base.induced import vorticity

import os
import sys
import yaml
import re
import csv
import pandas

case = 'remesh'
case_dir = os.getcwd()
methods = ['induced', 'M4prime', 'Stock', 'velocity', 'multigrid_medTH', 'multigrid_lowTH', 'multigrid_highTH']
# methods = ['induced']

compression_list = np.array([0.001, 0.01, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 0.99, 0.999])
# compression_list = np.array([0.95, 0.99])

 # ___________ Params ____________
remesh_params = {'Csigma' : 1., 'margin' : 0.05, 'sigma_max' : 0.1,
                'radius_check' : 1.2, 'ignore':False,
                    'refinement_factor':2.}

for method in methods:
    print(f'______________________ RUNNING METHOD {method} _________________________-')
    folder_name = 'order_of_accuracy_' + method
    plots_dir = os.path.join(case_dir,folder_name)

    Path(plots_dir).mkdir(parents=True, exist_ok=True)

   
    if method == 'multigrid_medTH':
        remesh_params['type'] = 'induced'
        remesh_params['threshold'] = 1e-6
    elif method == 'multigrid_highTH':
        remesh_params['type'] = 'induced'
        remesh_params['threshold'] = 1e-4
    elif method == 'multigrid_lowTH':
        remesh_params['type'] = 'induced'
        remesh_params['threshold'] = 1e-8
    else:
        remesh_params['type'] = method


    gammaC = 1.0
    nu = 0.001
    vInfx = 0.
    vInfy = 0.
    overlap = 1.0
    sigma0 = 0.01
    hBlob = sigma0*overlap
    hMesh = 0.005

    analytical_solver = particle.LambOseenVortexParticle(gammaC,4.0,nu,0.,0.,vInfx,vInfy)

    xMesh, yMesh = np.meshgrid(np.arange(-1,1,hMesh),np.arange(-1,1,hMesh))
    xMeshFlat = xMesh.flatten()
    yMeshFlat = yMesh.flatten()

    xBlobs, yBlobs = np.meshgrid(np.arange(-1,1,hBlob),np.arange(-1,1,hBlob))
    xBlobFlat = xBlobs.flatten()
    yBlobFlat = yBlobs.flatten()

    xyBlobs = np.column_stack((xBlobFlat,yBlobFlat))

    initial_vorticity = analytical_solver.vorticity_initial_blobs(xyBlobs,hBlob)
    g_init = initial_vorticity*hBlob*hBlob
    sigma_init = np.full(g_init.shape, sigma0)

    def loop(compression):
        remesh_params['compression'] = compression
        if method in ['multigrid_medTH', 'multigrid_lowTH', 'multigrid_highTH']:
            x, y, g, sigma = multiscale_remesh(xBlobFlat, yBlobFlat, g_init, sigma_init, remesh_params)
        else:
            x, y, g, sigma = remesh(xBlobFlat, yBlobFlat, g_init, sigma_init, remesh_params)
        # 


        omega_actual = analytical_solver.vorticity(np.vstack((xMeshFlat, yMeshFlat)).T)
        omega_sim = vorticity(np.array(x), np.array(y), np.array(g), np.array(sigma), xEval = np.array(xMeshFlat), yEval = np.array(yMeshFlat))

        # rel_error = abs_error/max_omega*100
        comp_actual = len(x)/len(xBlobFlat)
        L2error = np.linalg.norm(omega_sim-omega_actual)/np.linalg.norm(omega_actual)
        Linferror = np.max(np.abs(omega_sim - omega_actual)/np.max(np.abs(omega_actual)))
        circ_error = np.sum(g) - gammaC
        return L2error, Linferror, comp_actual, circ_error
    
    compression_actual_list = np.zeros(compression_list.shape)

    L2_list = np.zeros(compression_list.shape)
    Linf_list = np.zeros(compression_list.shape)
    circ_list = np.zeros(compression_list.shape)

    for ind,comp in enumerate(compression_list):
        print(f'Running case for compression: {comp}')
        L2, Linf, comp_actual, circ_error = loop(comp)
        L2_list[ind] = L2
        Linf_list[ind] = Linf
        circ_list[ind] = circ_error
        compression_actual_list[ind] = comp_actual
        print(f'Saved error results of compression: {comp}')


    h_list = 1/np.sqrt(compression_actual_list)
    less_than_bound  = np.where(h_list<3)[0]
    h_list2=h_list[less_than_bound]
    compression_actual_list2 = compression_actual_list[less_than_bound]
    L2_list2 = L2_list[less_than_bound]
    Linf_list2 = Linf_list[less_than_bound]
    A = np.vstack((np.ones(len(compression_actual_list2)), np.log(h_list2))).T
    b1 = np.log(L2_list2)
    b2 = np.log(Linf_list2)
    sol1 = np.linalg.lstsq(A, b1)[0]
    sol2 = np.linalg.lstsq(A, b2)[0]
    print(f'Least squares fit of L2 error: {sol1}')
    print(f'Least squares fit of Linf error: {sol2}')

    np.savetxt(os.path.join(plots_dir,"results_.csv"), np.c_[compression_list, compression_actual_list, h_list, L2_list, Linf_list, circ_list], delimiter=' ', header = "Desired CR, Actual CR, SR, L2, Linf, circ_error")
    np.savetxt(os.path.join(plots_dir,"log_.csv"), np.c_[sol1, sol2], delimiter=' ')

    fig, ax = plt.subplots(1,1)

    ax.plot(h_list, L2_list, label = 'L2 error', marker='*')
    ax.plot(h_list, Linf_list, label='Linf error', marker='o')
    


    # list = np.linspace(0., 1., 1000)
    # xval = 10**list
    # yval = 10**(2*list)
    # ax.plot(xval, yval, 'r--', label='Expected Accuracy O(h^2)', alpha=0.5)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel('normalized spacing: 1/sqrt(CR)')
    ax.set_ylabel('Error')
    ax.grid(which='both')
    plt.legend()
    plt.savefig("{}/order_of_accuracy_{}.png".format(plots_dir,case), dpi=300, bbox_inches="tight")
    plt.close(fig)

    fig, ax = plt.subplots(1,1)
    ax.plot(h_list, circ_list, label='Circulation Error', marker='x')
    ax.set_xlabel('normalized spacing: 1/sqrt(CR)')
    ax.set_ylabel('Error')
    ax.grid(which='both')
    plt.legend()
    plt.savefig("{}/circulation_error_{}.png".format(plots_dir,case), dpi=300, bbox_inches="tight")
    plt.close(fig)

    fig, ax = plt.subplots(1,1)

