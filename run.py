import numpy as np 
import timeit
import matplotlib.pyplot as plt
from pathlib import Path
import blobs_solver2 as pHyFlow
# import pHyFlow
# import VorticityFoamPy as foam
import solvers.particle as particle

import os
import sys
import yaml
import re
import csv
import pandas


#---------------Current directory and paths---------------------F---------------
arg = sys.argv
if len(arg) > 2:
    raise Exception("More than two arguments inserted!")
if len(arg) <= 1:
    raise Exception("No config file specificed!")
configFile = arg[1]

#-----------------------Config the yaml file ------------------
loader = yaml.SafeLoader
loader.add_implicit_resolver(
    u'tag:yaml.org,2002:float',
    re.compile(u'''^(?:
     [-+]?(?:[0-9][0-9_]*)\\.[0-9_]*(?:[eE][-+]?[0-9]+)?
    |[-+]?(?:[0-9][0-9_]*)(?:[eE][-+]?[0-9]+)
    |\\.[0-9_]+(?:[eE][-+][0-9]+)?
    |[-+]?[0-9][0-9_]*(?::[0-5]?[0-9])+\\.[0-9_]*
    |[-+]?\\.(?:inf|Inf|INF)
    |\\.(?:nan|NaN|NAN))$''', re.X),
    list(u'-+0123456789.'))
config = yaml.load(open(os.path.join(configFile)),Loader=loader)

case_dir = os.getcwd()
data_dir = os.path.join(case_dir,config["data_folder"])
plots_dir = os.path.join(case_dir,config["plots_folder"])

Path(data_dir).mkdir(parents=True, exist_ok=True)
Path(plots_dir).mkdir(parents=True, exist_ok=True)

start_index = configFile.index('_') + 1
end_index = configFile.index('.')
# case = configFile[start_index:end_index]
case = config['case']

# #--------------Copy values from the config file--------------------

nPlotPoints = config["nPlotPoints"]
xMinPlot = config["xMinPlot"]
xMaxPlot = config["xMaxPlot"]
yMinPlot = config["yMinPlot"]
yMaxPlot = config["yMaxPlot"]

writeInterval_plots = config["writeInterval_plots"]
run_analytical_flag = config["run_analytical_flag"]
plot_flag = config["plot_flag"]

vInfx = config["vInfx"]
vInfy = config["vInfy"]
vInf = np.array([vInfx, vInfy])

nu = config["nu"]
gammaC = config["gammaC"]

coreSize = config["coreSize"]

overlap = config["overlap"]

deltaTc = config["deltaTc"]
nTimeSteps = config["nTimeSteps"]

if coreSize == "fixed":
    sigma = float(2.0*(np.sqrt(6.0*nu*deltaTc)))
    hBlob = sigma*overlap
else :
    sigma = np.array([2.0*(np.sqrt(6.0*nu*deltaTc))])
    hBlob = sigma[0]*overlap

compressionFlag = config['compressionFlag']
compression_method = config["compression_method"]
support_method = config["support_method"]
compression_params = config["compression_params"]
support_params = config["support_params"]

compression_stride = config['compression_stride']
def create_linear(x0, y0, x1, y1):
    if x0 == x1:
        raise ValueError('two reference x values cannot be the same!')
    
    a = (y1 - y0) / (x1 - x0)
    b = y0 - a*x0
    return lambda index, x, y, g, sigma: max(a*x[index].mean()+b, min(y1, y0))

# compression_params['compression'] = create_linear(**config['compression_func_values'])

# step, start, end = config['support_func_values'].values()
# support_params['r'] = np.arange(start, end+step, step)

#--------------------Plot parameters--------------------------------
xplot,yplot = np.meshgrid(np.linspace(xMinPlot,xMaxPlot,nPlotPoints),np.linspace(yMinPlot,yMaxPlot,nPlotPoints))
xplotflat = xplot.flatten()
yplotflat = yplot.flatten()
xyPlot = np.column_stack((xplotflat, yplotflat))

#------------------Parameters for blobs-----------------------------
computationParams = {'hardware':config['hardware'], 'method':config['method']}

blobControlParams = {'methodPopulationControl':config['method_popControl'],'typeOfThresholds':'relative', 'stepRedistribution':config['stepRedistribution'],\
                     'stepPopulationControl':config['stepPopulationControl'], 'gThresholdLocal': config['gThresholdLocal'],\
                     'gThresholdGlobal':config['gThresholdGlobal']}

blobDiffusionParams = {'method' : config['method_diffusion']}

timeIntegrationParams = {'method':config['time_integration_method']}

kernelParams = {'kernel' : config['kernel'], 'coreSize' : config['coreSize']}

avrmParams = {'useRelativeThresholds':True, 'ignoreThreshold' : 1e-12, 'adaptThreshold': 1e-8, 'Clapse' : 0.01,\
                       'merge_flag':True, 'stepMerge':1, 'mergeThreshold':0.0001}

compressionParams = {'method':compression_method, 'support':support_method, 'methodParams':compression_params,\
                     'supportParams':support_params}

xShift = config['xShift'] 
yShift = config['yShift']
#-------------------------------------------------------------------
x0 = np.array([0.0])
y0 = np.array([0.0])

analytical_solver = particle.LambOseenVortexParticle(gammaC,4.0,nu,x0[0],y0[0],vInfx,vInfy)

xBlobs, yBlobs = np.meshgrid(np.arange(-1,1,hBlob),np.arange(-1,1,hBlob))
xBlobFlat = xBlobs.flatten()
yBlobFlat = yBlobs.flatten()

xyBlobs = np.column_stack((xBlobFlat,yBlobFlat))

initial_vorticity = analytical_solver.vorticity_initial_blobs(xyBlobs,hBlob)
g = initial_vorticity*hBlob*hBlob

# g = np.array([gammaC])
if coreSize == 'variable':
    sigma = np.full(len(g),hBlob)

wField = (xBlobFlat, yBlobFlat, g)

# generate the blobs
blobs = pHyFlow.blobs.Blobs(wField,vInf,nu,deltaTc,sigma,overlap,xShift,yShift,
                            kernelParams=kernelParams,
                            diffusionParams=blobDiffusionParams,
                            velocityComputationParams=computationParams,
                            timeIntegrationParams=timeIntegrationParams,
                            blobControlParams=blobControlParams,
                            avrmParams=avrmParams,
                            mergingParams=compressionParams)


header = ['Time', 'NoBlobs', 'Evolution_time', 'Circulation']
with open('{}/times_{}.csv'.format(data_dir,case), 'w', encoding='UTF8') as f:
    writer = csv.writer(f)
    # write the header
    writer.writerow(header)

ux, uy = blobs.evaluateVelocity(xplotflat,yplotflat)
omega = blobs.evaluateVorticity(xplotflat,yplotflat)
np.savetxt(os.path.join(data_dir,"results_{}_{n:06d}.csv".format(case,n=0)), np.c_[xplotflat,yplotflat,ux,uy,omega], delimiter=' ')
if coreSize == 'variable':
    np.savetxt(os.path.join(data_dir, "blobs_{}_{n:06d}.csv".format(case,n=0)), np.c_[blobs.x, blobs.y, blobs.g, blobs.sigma], delimiter= ' ')
else :
    np.savetxt(os.path.join(data_dir, "blobs_{}_{n:06d}.csv".format(case,n=0)), np.c_[blobs.x, blobs.y, blobs.g], delimiter= ' ')
if run_analytical_flag:
    analytical_ux, analytical_uy = analytical_solver.velocity(xyPlot)
    analytical_vorticity = analytical_solver.vorticity(xyPlot)
    np.savetxt(os.path.join(data_dir,"results_analytical_{n:06d}.csv".format(n=0)), np.c_[xplotflat,yplotflat,analytical_ux,analytical_uy,analytical_vorticity], delimiter=' ')

for timeStep in range(1,nTimeSteps+1):
    time_start = timeit.default_timer()
    blobs.evolve()
    if compressionFlag and (timeStep%compression_stride == 0 or timeStep == 1):
        print('----------------Performing Compression--------------')
        nbefore = blobs.numBlobs
        blobs._compress()
        blobs.populationControl()
        nafter = blobs.numBlobs
        print(f'removed {nbefore-nafter} particles')
        print(f'current number of particles: {nafter}')
    time_end = timeit.default_timer()
    print("Time to evolve in timeStep {} is {}".format(timeStep,time_end - time_start))

    evolution_time = time_end - time_start
    T = timeStep*deltaTc

    data = [T,blobs.numBlobs,evolution_time, blobs.g.sum()]

    with open('{}/times_{}.csv'.format(data_dir,case), 'a', encoding='UTF8') as f:
        writer = csv.writer(f)
        writer.writerow(data)

    if run_analytical_flag:
        start_time_analytical = timeit.default_timer()
        analytical_solver.evolve(deltaTc)
        end_time_analytical = timeit.default_timer()

        print("Analytical run in {} seconds".format(end_time_analytical - start_time_analytical))

    if timeStep%writeInterval_plots==0 :
        ux, uy = blobs.evaluateVelocity(xplotflat,yplotflat)
        omega = blobs.evaluateVorticity(xplotflat,yplotflat)
        print(np.max(np.abs(omega)))

        np.savetxt(os.path.join(data_dir,"results_{}_{n:06d}.csv".format(case,n=timeStep)), np.c_[xplotflat,yplotflat,ux,uy,omega], delimiter=' ')
        if coreSize == 'variable':
            np.savetxt(os.path.join(data_dir, "blobs_{}_{n:06d}.csv".format(case,n=timeStep)), np.c_[blobs.x, blobs.y, blobs.g, blobs.sigma], delimiter= ' ')
        else :
            np.savetxt(os.path.join(data_dir, "blobs_{}_{n:06d}.csv".format(case,n=timeStep)), np.c_[blobs.x, blobs.y, blobs.g], delimiter= ' ')

        if run_analytical_flag:
            analytical_ux, analytical_uy = analytical_solver.velocity(xyPlot)
            analytical_vorticity = analytical_solver.vorticity(xyPlot)
            np.savetxt(os.path.join(data_dir,"results_analytical_{n:06d}.csv".format(n=timeStep)), np.c_[xplotflat,yplotflat,analytical_ux,analytical_uy,analytical_vorticity], delimiter=' ')
    


if plot_flag == True:
    uxNorm = np.array([])
    uyNorm = np.array([])
    omegaNorm = np.array([])
    t_norm = np.array([])
    #Line plots
    times_file = os.path.join(data_dir,"times_{}.csv".format(case))
    times_data = pandas.read_csv(times_file)

    time = times_data['Time']
    noBlobs = times_data['NoBlobs']
    evolution_time = times_data['Evolution_time']
    circulation = times_data['Circulation']

    fig, ax = plt.subplots(1,1,figsize=(6,6))
    ax.plot(time,noBlobs, label='No of Particles')
    plt.grid(color = '#666666', which='major', linestyle = '--', linewidth = 0.5)
    plt.minorticks_on()
    plt.grid(visible=True, which='minor', color='#999999', linestyle='-', alpha=0.2)
    plt.title('Total number of particles')
    plt.ylabel('Particles')
    plt.xlabel('time $(sec)$')
    plt.legend()
    plt.savefig("{}/number_of_particles_{}.png".format(plots_dir,case), dpi=300, bbox_inches="tight")

    fig, ax = plt.subplots(1,1,figsize=(6,6))
    ax.plot(time,circulation- gammaC, label='Circulation deficit')
    plt.grid(color = '#666666', which='major', linestyle = '--', linewidth = 0.5)
    plt.minorticks_on()
    plt.grid(visible=True, which='minor', color='#999999', linestyle='-', alpha=0.2)
    plt.title('absolute error in circulation')
    plt.ylabel('circulation')
    plt.xlabel('time $(sec)$')
    plt.legend()
    plt.savefig("{}/circulation_error_{}.png".format(plots_dir,case), dpi=300, bbox_inches="tight")

    fig = plt.subplots(figsize=(6,6))
    index = np.arange(len(evolution_time))
    width = deltaTc
    lagrangian = plt.bar(index[1:]*deltaTc, evolution_time[1:], width)
    plt.ylabel('Time (s)')
    plt.xlabel('Simulation time (s)')
    plt.title('Evolution time')
    plt.savefig("{}/times_{}.png".format(plots_dir,case), dpi=300, bbox_inches="tight")

    for timeStep in range(nTimeSteps+1):
        if timeStep%writeInterval_plots == 0:
            ####Fields
            lagrangian_file = os.path.join(data_dir,'results_{}_{n:06d}.csv'.format(case,n=timeStep))
            lagrangian_data = np.genfromtxt(lagrangian_file)

            xplot = lagrangian_data[:,0]
            yplot = lagrangian_data[:,1]
            length = int(np.sqrt(len(xplot)))
            xPlotMesh = xplot.reshape(length,length)
            yPlotMesh = yplot.reshape(length,length)

            lagrangian_ux = lagrangian_data[:,2]
            lagrangian_uy = lagrangian_data[:,3]
            lagrangian_omega = lagrangian_data[:,4]

            analytical_file = os.path.join(data_dir,'results_analytical_{n:06d}.csv'.format(case,n=timeStep))
            analytical_data = np.genfromtxt(analytical_file)

            analytical_ux = analytical_data[:,2]
            analytical_uy = analytical_data[:,3]
            analytical_omega = analytical_data[:,4]

            xTicks = np.linspace(-2,2,5)
            yTicks = np.linspace(-2,2,5)

            fig, ax = plt.subplots(1,1,figsize=(6,6))
            ax.set_aspect("equal")
            ax.set_xticks(xTicks)
            ax.set_yticks(yTicks)
            plt.grid(color = '#666666', which='major', linestyle = '--', linewidth = 0.5)
            plt.minorticks_on()
            plt.grid(visible=True, which='minor', color='#999999', linestyle='-', alpha=0.2)
            ax.set_xlabel("x")
            ax.set_ylabel("y")
            cax = ax.contourf(xPlotMesh,yPlotMesh,lagrangian_omega.reshape(length,length),levels=100,cmap='RdBu',extend="both")
            cbar = fig.colorbar(cax,format="%.4f")
            cbar.set_label("Vorticity (1/s)")
            plt.tight_layout()
            plt.savefig("{}/vorticity_{}_{}.png".format(plots_dir,case,timeStep), dpi=300, bbox_inches="tight")
            plt.close(fig)
            
            fig, ax = plt.subplots(1,1,figsize=(6,6))
            ax.set_aspect("equal")
            ax.set_xticks(xTicks)
            ax.set_yticks(yTicks)
            plt.grid(color = '#666666', which='major', linestyle = '--', linewidth = 0.5)
            plt.minorticks_on()
            plt.grid(visible=True, which='minor', color='#999999', linestyle='-', alpha=0.2)
            ax.set_xlabel("x")
            ax.set_ylabel("y")
            cax = ax.contourf(xPlotMesh,yPlotMesh,lagrangian_ux.reshape(length,length),levels=100,cmap='RdBu',extend="both")
            cbar = fig.colorbar(cax,format="%.4f")
            cbar.set_label("Velocity (1/s)")
            plt.tight_layout()
            plt.savefig("{}/velocity_{}_{}.png".format(plots_dir,case,timeStep), dpi=300, bbox_inches="tight")
            plt.close(fig)

            if run_analytical_flag==True:
                fig, ax = plt.subplots(1,1,figsize=(6,6))
                ax.set_aspect("equal")
                ax.set_xticks(xTicks)
                ax.set_yticks(yTicks)
                plt.grid(color = '#666666', which='major', linestyle = '--', linewidth = 0.5)
                plt.minorticks_on()
                plt.grid(visible=True, which='minor', color='#999999', linestyle='-', alpha=0.2)
                ax.set_xlabel("x")
                ax.set_ylabel("y")
                cax = ax.contourf(xPlotMesh,yPlotMesh,analytical_omega.reshape(length,length),levels=100,cmap='RdBu',extend="both")
                cbar = fig.colorbar(cax,format="%.4f")
                cbar.set_label("Vorticity (1/s)")
                plt.tight_layout()
                plt.savefig("{}/vorticity_analytical_{}.png".format(plots_dir,timeStep), dpi=300, bbox_inches="tight")
                plt.close(fig)

                fig, ax = plt.subplots(1,1,figsize=(6,6))
                ax.set_aspect("equal")
                ax.set_xticks(xTicks)
                ax.set_yticks(yTicks)
                plt.grid(color = '#666666', which='major', linestyle = '--', linewidth = 0.5)
                plt.minorticks_on()
                plt.grid(visible=True, which='minor', color='#999999', linestyle='-', alpha=0.2)
                ax.set_xlabel("x")
                ax.set_ylabel("y")
                cax = ax.contourf(xPlotMesh,yPlotMesh,analytical_ux.reshape(length,length),levels=100,cmap='RdBu',extend="both")
                cbar = fig.colorbar(cax,format="%.4f")
                cbar.set_label("Velocity (1/s)")
                plt.tight_layout()
                plt.savefig("{}/velocity_analytical_{}.png".format(plots_dir,timeStep), dpi=300, bbox_inches="tight")
                plt.close(fig)

#### Errors
            omegaScale = np.max(np.abs(analytical_omega))

            omega_error = ((lagrangian_omega - analytical_omega)/(omegaScale))*100
            fig, ax = plt.subplots(1,1,figsize=(6,6))
            ax.set_aspect("equal")
            ax.set_xticks(xTicks)
            ax.set_yticks(yTicks)
            plt.grid(color = '#666666', which='major', linestyle = '--', linewidth = 0.5)
            plt.minorticks_on()
            plt.grid(visible=True, which='minor', color='#999999', linestyle='-', alpha=0.2)
            ax.set_xlabel("x")
            ax.set_ylabel("y")
            cax = ax.contourf(xPlotMesh,yPlotMesh,omega_error.reshape(length,length),levels=100,cmap='jet',extend="both")
            cbar = fig.colorbar(cax,format="%.4f")
            cbar.set_label("Vorticity error (%) (1/s)")
            plt.tight_layout()
            plt.savefig("{}/vorticity_error_{}_{}.png".format(plots_dir,case,timeStep), dpi=300, bbox_inches="tight")
            plt.close(fig)


#### Blobs distribution

            blobs_file = os.path.join(data_dir,'blobs_{}_{n:06d}.csv'.format(case,n=timeStep))
            blobs_data = np.genfromtxt(blobs_file)

            blobs_x = blobs_data[:,0]
            blobs_y = blobs_data[:,1]
            blobs_g = blobs_data[:,2]

            if coreSize == 'variable':
                blobs_sigma = blobs_data[:,3]

                fig, ax = plt.subplots(1,1,figsize=(6,6))
                ax.scatter(blobs_x,blobs_y,c=blobs_g, s= blobs_sigma*30)
                plt.savefig("{}/blobs_{}_{}.png".format(plots_dir,case,timeStep), dpi=300, bbox_inches="tight")
                plt.close(fig)
            else:
                fig, ax = plt.subplots(1,1,figsize=(6,6))
                ax.scatter(blobs_x,blobs_y,c=blobs_g, s=0.2)
                plt.savefig("{}/blobs_{}_{}.png".format(plots_dir,case,timeStep), dpi=300, bbox_inches="tight")
                plt.close(fig)


#### L2 errors in vorticity and velocity


            uxNorm = np.append(uxNorm,np.linalg.norm(lagrangian_ux-analytical_ux)/np.linalg.norm(analytical_ux))
            uyNorm = np.append(uyNorm,np.linalg.norm(lagrangian_uy-analytical_uy)/np.linalg.norm(analytical_uy))
            omegaNorm = np.append(omegaNorm,np.linalg.norm(lagrangian_omega-analytical_omega)/np.linalg.norm(analytical_omega))
            t_norm = np.append(t_norm,deltaTc*timeStep)
    fig, ax = plt.subplots(1,1,figsize=(6,6))
    plt.plot(t_norm,uxNorm, label='ux-L2 error')
    plt.plot(t_norm,uyNorm, label='uy-L2 error')
    plt.plot(t_norm,omegaNorm, label='omega-L2 error')
    plt.grid(color = '#666666', which='major', linestyle = '--', linewidth = 0.5)
    plt.minorticks_on()
    plt.grid(visible=True, which='minor', color='#999999', linestyle='-', alpha=0.2)
    plt.ylabel('L2-error')
    plt.xlabel('time')
    plt.legend()
    plt.savefig("{}/L2_error_{}.png".format(plots_dir,case), dpi=300, bbox_inches="tight")
    plt.close(fig)