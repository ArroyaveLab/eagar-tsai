# Import the required modules
import os
import sys
import numpy as np
import pandas as pd
from scipy.integrate import *
import scipy.version
import ctypes
import time
from itertools import compress
import os.path as path
import concurrent.futures

#import resource
#import gc
#from pympler import muppy,summary

def _run_chunk(params):
    chunk_num, data, domain_um, spatial_res_um, out_dir = params
    data = data.reset_index(drop=True)

    # ========== SIMULATION ====================================================
    domain = np.array(domain_um, dtype='f8')   # microns, expanded if needed
    spatialRes = float(spatial_res_um)         # microns
    # ==========================================================================

    # ========== BEAM and MATERIAL =============================================
    beam1, mat1 = beamFromCSV(data)
    # ==========================================================================

    integrationWarning()

    runSize = np.size(beam1.v, axis=0)
    for i in np.arange(runSize):
        sim1 = simParam(domain, spatialRes)
        try:
            (data.at[i, 'melt_length'],
             data.at[i, 'melt_width'],
             data.at[i, 'melt_depth'],
             data.at[i, 'peakT'],
             data.at[i, 'minT']) = eagarTsaiParam(beam1, mat1, sim1, i)
        except Exception as exc:
            print(f'Error in chunk {chunk_num}: {exc}')

    # Add micron columns for convenience (meters -> microns)
    data['melt_length_um'] = data['melt_length'] * 1.0e6
    data['melt_width_um'] = data['melt_width'] * 1.0e6
    data['melt_depth_um'] = data['melt_depth'] * 1.0e6

    if out_dir is not None:
        os.makedirs(out_dir, exist_ok=True)
        data.to_csv(path.join(out_dir, f'ET_v3_OUT_{chunk_num}.csv'), index=False)

    return data


def compute_melt_pool(
    data,
    domain_um=(1200.0, 1200.0, 1000.0),
    spatial_res_um=1.0,
    chunk_size=1,
    workers=None,
    out_dir=None,
):
    """
    Compute melt pool dimensions for each row in a DataFrame.

    Required columns:
      Velocity_m/s, Power, Beam_diameter_m, Absorptivity,
      T_liquidus, thermal_cond_liq, Density_kg/m3, Cp_J/kg
    """
    if not isinstance(data, pd.DataFrame):
        raise TypeError("data must be a pandas DataFrame")

    chunks = [data[i:i + chunk_size] for i in range(0, data.shape[0], chunk_size)]
    params = [(i, chunk, domain_um, spatial_res_um, out_dir) for i, chunk in enumerate(chunks)]

    if workers is None or workers <= 1:
        results = [_run_chunk(p) for p in params]
    else:
        with concurrent.futures.ProcessPoolExecutor(workers) as executor:
            results = list(executor.map(_run_chunk, params))

    return pd.concat(results, ignore_index=True)


def beamFromCSV(data):
    """ Import beam parameters from csv file """
    #data = pd.read_csv(inputfile,header=0)

	# These we will iterate through
#    v = data['Velocity_m/s'] # m/s
#    P = data['P [W]'] # Watts
#    twoSigma = data["Beam_diameter_m"]  #m # meters
#    A = data['Absorptivity'] # >0, <1
#    tMelt = data['PROP LT (K)']
#    k = data['Prop Liquidus Thermal conductivity (W/(mK))']
#    rho = data["Prop RT Density (kg/m3)"]
#    cp = data['eff_Cp_JkgK']
    
    v = data['Velocity_m/s'] # m/s
    P = data['Power'] # Watts
    twoSigma = data["Beam_diameter_m"]  #m # meters
    A = data['Absorptivity'] # >0, <1
    tMelt = data['T_liquidus']
    k = data['thermal_cond_liq']
    rho = data["Density_kg/m3"]
    cp = data['Cp_J/kg']
    
    # Create the beam object
    return beam(twoSigma,P,v,A), material(tMelt,k,rho,cp)

def eagarTsaiParam(beam,material,simParam,i):
    """ Function that runs each iteration of the EagarTsai Simulation """
    # This version of the code uses a re-formulation of EagarTsai
    # By Sasha Rubenchik - LLNL 2015
    
    # Unpack a few variables
    # Material
    #print(material.tMelt[:])
    tMelt = material.tMelt[i]
    k = material.k[i]
    rho = material.rho[i]
    cp = material.cp[i]
    alpha = k/(rho*cp)
    
    # Beam (this varies each run)
    P = beam.P[i]
    A = beam.A[i]
    v = beam.v[i]
    sigma = beam.sigma[i]
    
    # Simulation params
    delta = simParam.spatialRes
    
    # Now lets define the domain we are going to evaluate temperature over
    # Define the minimums and maximums
    xMin = round(-1.0 * beam.twoSigma[i] * 1.5, 5)
    xMax = simParam.domain[0]
    yMin = 0.0
    yMax = simParam.domain[1]
    zMin = -1.0 * simParam.domain[2]
    zMax = 0.0
    # Find the number of intervals
    nx = int(np.round(abs(xMax-xMin)/delta))+1
    ny = int(np.round(abs(yMax-yMin)/delta))+1
    nz = int(np.round(abs(zMax-zMin)/delta))+1
    # Create the range arrays
    nxrange = np.linspace(xMin,xMax,nx)
    nyrange = np.linspace(yMin,yMax,ny)
    nzrange = np.linspace(zMin,zMax,nz)
    
    # We need to do some checking about which method to use
    # We need to see which version of scipy we are running
    # The compiled integral will only work on scipy 0.15.1 and later
    sciVer = scipy.version.version.split(".")
    # Find the platform type
    osType = sys.platform
    
    if int(sciVer[1]) >= 1:
        # Then we can use compiled code to run the integration
        if osType == 'darwin':   # macOS
            libpath = 'libeagar_tsai_integrand.dylib'
    
        elif osType.startswith('linux'):   # Linux
            libpath = './libeagar_tsai_integrand.so'
        
        elif (osType == 'win32') or (osType == 'cygwin') or (osType == 'msys'):
            # Windows: use DLL if present, otherwise fall back to interpreted code
            libpath = 'libeagar_tsai_integrand.dll'

    else:
        # Old version of SciPy and possibly python
        libpath = 0
        
    # Run integral
    tplanexy, tplanexz = runIntegrate(nxrange,nyrange,nzrange,nx,ny,nz,alpha,sigma,k,v,A,P,libpath)

    # Find the peak temperature
    peakT = np.amax(tplanexy)
    minT = np.amin(tplanexy)
    # Now check to see if the peak temperature is hotter than Tmelt
    if peakT > tMelt:
        # Then there is a melt pool, find the size
        # Section to actually extract the metrics regarding melt pool
        # From the XY Plane
        # Want to find the length
        meltXInd = np.squeeze(np.where(tplanexy[0,:] > tMelt))
        melt_length = np.amax(nxrange[meltXInd])-np.amin(nxrange[meltXInd])
        melt_trail_length = abs(np.amin(nxrange[meltXInd]))
        
        # Now want to find the width + depth (can do it in same outer loop)
        yLength = 0
        zLength = 0
        for i1 in np.arange(np.size(meltXInd,axis=0)):
            meltYInd = np.squeeze(np.where(tplanexy[:,meltXInd[i1]] > tMelt))
            tmpYLength = np.amax(nyrange[meltYInd]) - np.amin(nyrange[meltYInd])
            if tmpYLength > yLength:
                yLength = tmpYLength
            
            meltZInd = np.squeeze(np.where(tplanexz[:,meltXInd[i1]] > tMelt))
            tmpZLength = np.amax(nzrange[meltZInd]) - np.amin(nzrange[meltZInd])
            if tmpZLength > zLength:
                zLength = tmpZLength
		
        # Test to see if the domain is the correct size
        if np.isclose(np.amax(nxrange[meltXInd]),xMax):
            # Then the x domain is not long enough
            print ("The x domain (length) is not large enough, increasing size and re-running")
            simParam.domain[0] += sigma
            # Re-run the analysis
            # melt_width, melt_depth, melt_length = eagarTsaiParam(beam,material,simParam,i)
            melt_length, melt_width, melt_depth = eagarTsaiParam(beam,material,simParam,i)
        elif np.isclose(yLength,abs(yMax - yMin)):
            # Then the y domain is not large enough
            print ("The y domain (width) is not large enough, increasing size and re-running")
            simParam.domain[1] += sigma
            # Re-run the analysis
            melt_length, melt_width, melt_depth = eagarTsaiParam(beam,material,simParam,i)
        elif np.isclose(zLength,abs(zMax-zMin)):
            print ("The z domain (depth) is not large enough, increasing size and re-running")
            simParam.domain[2] += sigma
            # Re-run the analysis
            melt_length, melt_width, melt_depth = eagarTsaiParam(beam,material,simParam,i)
        else:
            # All is good, the melt pool wasnt clipped in the domain
            # Return the values
            melt_width = yLength * 2
            melt_depth = zLength
    else:
        # Then there is no melt pool, return 0 lengths
        melt_width = 0.0
        melt_depth = 0.0
        melt_length = 0.0

    del tplanexy, tplanexz
    return melt_length, melt_width, melt_depth, peakT, minT

def runIntegrate(nxrange,nyrange,nzrange,nx,ny,nz,alpha,sigma,k,v,A,P,libpath):
    """ Function to run the integration """
    # Check if using compiled or interpreted code
    if libpath == 0 or (isinstance(libpath, str) and not path.exists(libpath)):
        # using interpreted code
        func = eagar_tsai_integrand
    else:
        # using compiled code
        lib = ctypes.CDLL(libpath) # Use absolute path to shared library
        func = lib.eagar_tsai_integrand # Assign specific function to name func (for simplicity)
        func.restype = ctypes.c_double
        func.argtypes = (ctypes.c_int, ctypes.c_double)

    # Set up invariant parameters
    # Define the starting temperature
    t0 = 300.0 # Kelvin
    Ts = (A*P)/(np.pi*(k/alpha)*np.sqrt(np.pi*alpha*v*(sigma**3)))
    p = alpha/(v*sigma)
    
    # Create the two temperature planes
    tplanexy = np.zeros(nx*ny,dtype='f8').reshape(ny,nx)
    tplanexz = np.zeros(nx*nz,dtype='f8').reshape(nz,nx)
    
    # Run the integration
    for i1 in np.arange(nx):
        x = nxrange[i1]/sigma # make dimensional
        for i2 in np.arange(ny):
            y = nyrange[i2]/sigma # make dimensional
            tmpTemp = quad(func,0.,np.inf,args=(x, y, 0.0, p))
            tplanexy[i2,i1] = t0 + Ts*tmpTemp[0]
        for i3 in np.arange(nz):
            z = nzrange[i3]/np.sqrt((alpha*sigma)/v) # make dimensional
            tmpTemp = quad(func,0.,np.inf,args=(x, 0.0, z, p))
            tplanexz[i3,i1] = t0 + Ts*tmpTemp[0]

    # Return the temperature planes
    return tplanexy, tplanexz

def eagar_tsai_integrand(t, x, y, z, p):
    # This is Sasha's formulation
    intpre = 1.0/((4*p*t + 1)*np.sqrt(t))
    intexp = (-(z**2)/(4*t))-(((y**2)+(x-t)**2)/(4*p*t + 1))
    return intpre * np.exp(intexp)

class simParam():
    def __init__(self,domain,spatialRes):
        """ Define a simulation parameter object """
        self.domain = domain / 1.e6 # Convert into meters
        self.spatialRes = spatialRes / 1.e6 # Convert into meters

class beam():
    def __init__(self,twoSigma,P,v,A):
        """ Define a beam object """
        self.twoSigma = twoSigma
        # The next sigma is an altered version for Sasha's ET re-interp
        self.sigma = np.sqrt(2.0) * (self.twoSigma / 2.0)
        self.P = P
        self.v = v
        self.A = A

class material():
    def __init__(self,tMelt,k,rho,cp):
        """ Define a material object """
        self.tMelt = tMelt
        self.k = k
        self.rho = rho
        self.cp = cp

def integrationWarning():
    """ Function that warns the user about the integration type """
    # We need to do some checking about which method to use
    # We need to see which version of scipy we are running
    # The compiled integral will only work on scipy 0.15.1 and later
    sciVer = scipy.version.version.split(".")
    
    # Find the platform type
    osType = sys.platform
    if int(sciVer[1]) >= 1:
        # Then we can use compiled code to run the integration
        if osType == 'darwin': # Then you're cool because that's a Mac!
            print ("Using compiled integration code, should be faster")
        elif osType in 'linux2': # It's linux, I guess that's OK
            print ("Using compiled integration code, should be faster")
        elif (osType=='win32') or (osType=='cygwin'):
            # Sorry, windows sucks, run the integration using normal interpreted code
            print ("Using interpreted integration code, will be slow")
            print ("Consider switching to a Mac or Linux")
    else:
        print ("Using old version of Python and SciPy")
        print ("Please consider switching to Python 2.7.10 and SciPy 0.15.1")
        print ("Using interpreted integration code, will be slow")

if __name__ == "__main__":

    ##########################################################################

    results_df = pd.read_excel('et_input_data_example.xlsx')
#    results_df['Velocity_m/s'] = results_df['v']
#    results_df['Power'] = results_df['P']
    results_df['T_liquidus'] = results_df['PROP LT (K)']
    results_df['thermal_cond_liq'] = results_df['PROP LT THCD (W/(mK))']
    results_df['Density_kg/m3'] = results_df['PROP RT Density (kg/m3)']
    results_df['Cp_J/kg'] =  results_df['PROP LT C (J/(kg K))']
    results_df['Beam_diameter_m'] = results_df['Beam Diam (m)']
    #
    elements =	['W',	'Re',	'Nb'	,'Ta',	'Mo',	'Hf',	'V']
    savename = 'prop_out'
    ##########################################################################

    compute_melt_pool(
        results_df,
        chunk_size=1,
        workers=40,
        out_dir='CalcFiles',
    )


    #main()
    #sum1 = summary.summarize(muppy.get_objects())
    #summary.print_(sum1)
