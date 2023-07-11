import os
import sys
import shutil as sh

import numpy as np
import matplotlib.pyplot as plt

import IPython as IP
import yaml

from particle import LambOseenVortexParticle

def main():
    #-----------#
    # CONFIGURE #
    #-----------#

    arg = sys.argv
    caseDir = arg[1]
    configFile = arg[2]

    # Load config files
    config = yaml.load(open(configFile),Loader=yaml.FullLoader)
    cellCenters = np.genfromtxt("{}/cellCenterCoordinates.dat".format(caseDir))
    if "useFaceCells" in config:
        useFaceCells = config["useFaceCells"]
    else:
        useFaceCells = False
    if useFaceCells:
        faceCenters = np.genfromtxt("{}/faceCellCenterCoordinates.dat".format(caseDir))
    else:
        faceCenters = np.genfromtxt("{}/faceCenterCoordinates.dat".format(caseDir))

    # Parse particle parameters
    Gamma = float(config["particleStrength"])
    x0 = float(config["initialPositionX"])
    y0 = float(config["initialPositionY"])
    uInf = float(config["freestreamVelocityX"])
    vInf = float(config["freestreamVelocityY"])
    t0 = float(config["timeConstant"])
    nu = float(config["kinematicViscosity"])

    # Parse time and space parameters
    timeStep = float(config["timeStep"])
    endTime = float(config["endTime"])
    numCellsX, numCellsY = list(map(int, config["gridSize"].split("x")))

    # Parse output parameters
    writeInterval = int(config["writeInterval"])

    # Define derived parameters
    numSteps = int(endTime/timeStep)
    dataRoot = "{}/data/".format(caseDir)

    # Define particle
    particle = LambOseenVortexParticle(Gamma,t0, nu, x0, y0, uInf, vInf)

    #----------------#
    # INITIALISATION #
    #----------------#

    if not os.path.isdir(dataRoot):
        os.makedirs(dataRoot)

    # Compute velocity, pressure and pressure gradient at cell centers and at
    # face centers along the numerical boundary at t = 0
    omegaInternal = particle.vorticity(cellCenters)
    uInternal, vInternal = particle.velocity(cellCenters)
    # pInternal = particle.pressure(uInternal, vInternal)
    pInternal = particle.pressure(cellCenters)
    dpdxInternal, dpdyInternal = particle.pressure_gradient(cellCenters, uInternal, vInternal)
    uBoundary, vBoundary = particle.velocity(faceCenters)
    # pBoundary = particle.pressure(uBoundary, vBoundary)
    pBoundary = particle.pressure(faceCenters)
    dpdxBoundary, dpdyBoundary = particle.pressure_gradient(faceCenters, uBoundary, vBoundary)
    dudxBoundary, dudyBoundary, dvdxBoundary, dvdyBoundary = particle.velocity_gradient(faceCenters)

    # Export velocity pressure and pressure gradient at t = 0
    np.savetxt(os.path.join(dataRoot, "lagrangian_internal_{n:06d}.dat".format(n=0)), np.c_[uInternal, vInternal, pInternal, dpdxInternal, dpdyInternal, omegaInternal], delimiter=' ')
    np.savetxt(os.path.join(dataRoot, "lagrangian_boundary_{n:06d}.dat".format(n=0)), np.c_[uBoundary, vBoundary, pBoundary, dpdxBoundary, dpdyBoundary, dudxBoundary, dudyBoundary, dvdxBoundary, dvdyBoundary], delimiter=' ')

    #------------#
    # SIMULATION #
    #------------#

    for step in range(1,numSteps+1):
        # Time
        # t = step*timeStep
        Dt = timeStep

        # Evolve particle by one time step
        particle.evolve(Dt)

        # Compute velocity, pressure and pressure gradient at cell centers and
        # at face centers along the numerical boundary
        omegaInternal = particle.vorticity(cellCenters)
        uInternal, vInternal = particle.velocity(cellCenters)
        # pInternal = particle.pressure(uInternal, vInternal)
        pInternal = particle.pressure(cellCenters)
        dpdxInternal, dpdyInternal = particle.pressure_gradient(cellCenters, uInternal, vInternal)
        uBoundary, vBoundary = particle.velocity(faceCenters)
        # pBoundary = particle.pressure(uBoundary, vBoundary)
        pBoundary = particle.pressure(faceCenters)
        dpdxBoundary, dpdyBoundary = particle.pressure_gradient(faceCenters, uBoundary, vBoundary)
        dudxBoundary, dudyBoundary, dvdxBoundary, dvdyBoundary = particle.velocity_gradient(faceCenters)

        # Export velocity pressure and pressure gradient
        np.savetxt(os.path.join(dataRoot, "lagrangian_internal_{n:06d}.dat".format(n=step)), np.c_[uInternal, vInternal, pInternal, dpdxInternal, dpdyInternal, omegaInternal], delimiter=' ')
        np.savetxt(os.path.join(dataRoot, "lagrangian_boundary_{n:06d}.dat".format(n=step)), np.c_[uBoundary, vBoundary, pBoundary, dpdxBoundary, dpdyBoundary, dudxBoundary, dudyBoundary, dvdxBoundary, dvdyBoundary], delimiter=' ')

main()
