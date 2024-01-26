# -*- coding: utf-8 -*-
"""
Utilities and simulation code underlying the SpinFlux figures.
    - simulate_thetaI(roiparams, pinholeparams, sigmas, theta, thetaI_range, scenario) simulates the CRLB and underlying photon counts for a range of signal photon budget values thetaI.
    - simulate_thetaI_thetab(roiparams, pinholeparams, sigmas, theta, thetab_list, scenario) simulates the CRLB and underlying photon counts for a range of signal photon budget values thetaI and background counts thetab.
    - simulate_xp(roiparams, pinholeparams, sigmas, theta, xp_range, scenario) simulates the CRLB and underlying photon counts for a range of pinhole x-positions xp.
    - simulate_yp(roiparams, pinholeparams, sigmas, theta, yp_range, scenario) simulates the CRLB and underlying photon counts for a range of pinhole y-positions yp.
    - simulate_xp_rp(roiparams, pinholeparams, sigmas, theta, r_p_list, scenario) simulates the CRLB and underlying photon counts for a range of pinhole x-positions x_p and pinhole radii r_p.
    - simulate_yp_rp(roiparams, pinholeparams, sigmas, theta, r_p_list, scenario) simulates the CRLB and underlying photon counts for a range of pinhole y-positions y_p and pinhole radii r_p.
    - simulate_xp_2pattern(roiparams, pinholeparams, sigmas, theta, xp_range, scenario) simulates the CRLB and underlying photon counts for the separated 2-pattern configuration, for a range of pinhole x-positions xp.
    - simulate_yp_2pattern(roiparams, pinholeparams, sigmas, theta, yp_range, scenario) simulates the CRLB and underlying photon counts for the separated 2-pattern configuration, for a range of pinhole y-positions yp.
    - simulate_xp_s_patterns(roiparams, pinholeparams, sigmas, theta, s_list, scenario) simulates the CRLB and underlying photon counts in the x-separated 2-pattern configuration, for a range of focus x-positions x_p and separations s.
    - simulate_yp_s_patterns(roiparams, pinholeparams, sigmas, theta, s_list, scenario) simulates the CRLB and underlying photon counts in the y-separated 2-pattern configuration, for a range of focus y-positions y_p and separations s.
    - simulate_yp_xs_patterns(roiparams, pinholeparams, sigmas, theta, s_list, scenario) simulates the CRLB and underlying photon counts in the x-separated 2-pattern configuration, for a range of focus y-positions y_p and separations s.
    - coordinates_triangulation(center, r, rotation) returns pinhole center coordinates of the triangle configuration
    - simulate_xp_r_triangulation(roiparams, pinholeparams, sigmas, theta, r_list, scenario, rotation) simulates the CRLB and underlying photon counts in the equilateral triangle pattern configuration, for a range of focus x-positions x_p and spacings r.
    - simulate_thetaI_Nm_triangulation(roiparams, pinholeparams, sigmas, theta, thetaI_list, scenario) simulates the CRLB and underlying photon counts in the single- or two-pattern configuration, for a range of mesh pixel amounts Nm and signal photon counts thetaI
    - simulate_thetaI_Nm_triangulation(roiparams, pinholeparams, sigmas, theta, thetaI_list, r, scenario, rotation) simulates the CRLB and underlying photon counts in the equilateral triangle pattern configuration, for a range of mesh pixel amounts Nm and signal photon counts thetaI.
    - simulate_SMLM(roiparams, sigmapsf, theta) simulates the CRLB and underlying photon counts of the SMLM benchmark.
    - simulate(plotnumber) simulates the SpinFlux Cramér-Rao lower bound (CRLB) data underlying the figure specified in plotnumber.
"""

import numpy as np
from spinflux import SpinFlux
from smlm import SMLM

def simulate_thetaI(roiparams, pinholeparams, sigmas, theta, thetaI_range, scenario, illumination_doughnut=False, N_m = 100):
    '''
    simulate_thetaI(roiparams, pinholeparams, sigmas, theta, thetaI_range, scenario) simulates
    the CRLB and underlying photon counts for a range of signal photon budget values thetaI.
    
    Inputs:
        roiparams: [dx_i, Npixels]
            dx_i: Camera pixel size (m)
            Npixels: Amount of camera pixels per direction
        pinholeparams: (K x 3)-array [x_p, y_p, r_p]
            K: Amount of pinholes and patterns
            x_p: (K x 1) array of pinhole x-coordinates (m)
            y_p: (K x 1) array of pinhole y-coordinates (m)
            r_p: (K x 1) array of pinhole radii (m)
        sigmas: [sigmaillum, sigmapsf]
            sigmaillum: Standard deviation of Gaussian illumination PSF (m)
            sigmapsf: Standard deviation of Gaussian emission PSF (m)
        theta: [thetax, thetay, thetaI, thetab]
            thetax: Emitter x-position (m)
            thetay: Emitter y-position (m)
            thetaI: Expected signal photon budget under maximum illumination (photons) (overwritten)
            thetab: Expected background photon budget under maximum illumination (photons/(m*m))
        thetaI_range: range of thetaI values, overwrites thetaI.
        scenario: Simulation scenario under consideration, which models the normalization
            constant (i.e. to benchmark against SMLM) and pattern-dependent or 
            pattern-independent background.
            Options:
                'fixed_budget': Scenario for which the signal photon budget is exhausted
                'fixed_energy': Scenario for which the illumination power and time are constant per pattern
                'fixed_budget_adaptedbg': Scenario for which the signal photon budget is exhausted, with pattern-independent background
       
    Optional input:
        N_m: Amount of mesh pixels per direction, default set to 100.
        illumination_doughnut: If True, doughnut-shaped intensity patterns with a
            zero-intensity minimum in the center will be used for illumination.
            If False, Gaussian illumination will be used.            
    Outputs:
        CRLBx: Cramér-Rao lower bound (CRLB)
        musum: Illumination-intensity normalized total amount of photons in region of interest (ROI)
        signalsum: Illumination-intensity normalized total amount of signal photons in ROI
        Bisum: Sum of pattern-(in)dependent background constants over ROI
        EfBgsum: Illumination-intensity normalized total amount of background photons in ROI
    '''
    CRLBx = np.zeros(np.size(thetaI_range))
    musum = np.zeros(np.size(thetaI_range))
    signalsum = np.zeros(np.size(thetaI_range))
    Bisum = np.zeros(np.size(thetaI_range))
    EfBgsum = np.zeros(np.size(thetaI_range))
    
    for idx, thetaI in enumerate(thetaI_range):
        theta[2] = thetaI
        sf = SpinFlux(theta, sigmas, roiparams, pinholeparams, N_m = N_m, scenario=scenario, illumination_doughnut=illumination_doughnut)
        I, musum[idx], signalsum[idx], Bisum[idx], EfBgsum[idx] = sf.FisherMatrix(approximation='middle')
        I_inv = np.linalg.inv(I)
        CRLBx[idx] = np.sqrt(I_inv[0,0])
    return CRLBx, musum, signalsum, Bisum, EfBgsum

def simulate_thetaI_thetab(roiparams, pinholeparams, sigmas, theta, thetab_list, scenario, illumination_doughnut=False, filenameextension = 'thetaI_thetab', savepath = './SimData/', N_m = 100, start = 2, stop = 5, num = 20):
    '''
    simulate_thetaI_thetab(roiparams, pinholeparams, sigmas, theta, thetab_list, scenario) simulates
    the CRLB and underlying photon counts for a range of signal photon budget values thetaI
    and background counts thetab.
    
    Inputs:
        roiparams: [dx_i, Npixels]
            dx_i: Camera pixel size (m)
            Npixels: Amount of camera pixels per direction
        pinholeparams: (K x 3)-array [x_p, y_p, r_p]
            K: Amount of pinholes and patterns
            x_p: (K x 1) array of pinhole x-coordinates (m)
            y_p: (K x 1) array of pinhole y-coordinates (m)
            r_p: (K x 1) array of pinhole radii (m)
        sigmas: [sigmaillum, sigmapsf]
            sigmaillum: Standard deviation of Gaussian illumination PSF (m)
            sigmapsf: Standard deviation of Gaussian emission PSF (m)
        theta: [thetax, thetay, thetaI, thetab]
            thetax: Emitter x-position (m)
            thetay: Emitter y-position (m)
            thetaI: Expected signal photon budget under maximum illumination (photons) (overwritten)
            thetab: Expected background photon budget under maximum illumination (photons/(m*m)) (overwritten)
        thetab_list: range of thetab values, overwrites thetab.
        scenario: Simulation scenario under consideration, which models the normalization
            constant (i.e. to benchmark against SMLM) and pattern-dependent or 
            pattern-independent background.
            Options:
                'fixed_budget': Scenario for which the signal photon budget is exhausted
                'fixed_energy': Scenario for which the illumination power and time are constant per pattern
                'fixed_budget_adaptedbg': Scenario for which the signal photon budget is exhausted, with pattern-independent background
       
    Optional inputs:
        illumination_doughnut: If True, doughnut-shaped intensity patterns with a
            zero-intensity minimum in the center will be used for illumination.
            If False, Gaussian illumination will be used.
        filenameextension: string, describing simulation file to be saved to.
        savepath: string, describing simulation path to save to.
        N_m: Amount of mesh pixels per direction, default set to 100.
        start: starting order of logspace for signal photon count range thetaI
        stop: stopping order of logspace for signal photon count range thetaI
        num: amount of signal photon count values thetaI to sample
        
    Outputs:
        CRLBx: Cramér-Rao lower bound (CRLB)
        musum: Illumination-intensity normalized total amount of photons in region of interest (ROI)
        signalsum: Illumination-intensity normalized total amount of signal photons in ROI
        Bisum: Sum of pattern-(in)dependent background constants over ROI
        EfBgsum: Illumination-intensity normalized total amount of background photons in ROI
    '''
    print('Starting simulation thetaI_thetab')
    saverplist = 'thetablist_' + filenameextension
    np.save(savepath+saverplist, thetab_list)
    
    thetaI_range = np.logspace(start, stop, num)
    savexrange = 'thetaIrange_' + filenameextension
    np.save(savepath+savexrange,thetaI_range)
    
    CRLBx = np.zeros((np.size(thetab_list),np.size(thetaI_range)))
    musum = np.zeros((np.size(thetab_list),np.size(thetaI_range)))
    signalsum = np.zeros((np.size(thetab_list),np.size(thetaI_range)))
    Bisum = np.zeros((np.size(thetab_list),np.size(thetaI_range)))
    EfBgsum = np.zeros((np.size(thetab_list),np.size(thetaI_range)))
    
    for i in range(np.size(thetab_list)):
        print('Simulation thetaI_thetab: Evaluating CRLB for r_p with index ' + str(int(i)))
        thetab = thetab_list[i]
        theta[3]=thetab
        CRLBx[i,:], musum[i,:], signalsum[i,:], Bisum[i,:], EfBgsum[i,:] = simulate_thetaI(roiparams, pinholeparams, sigmas, theta, thetaI_range, scenario=scenario, illumination_doughnut=illumination_doughnut, N_m = 100)
        
    saveCRLBx = 'CRLBx_' + filenameextension
    savemusum = 'musum_' + filenameextension
    savesignalsum = 'signalsum_' + filenameextension
    saveBisum = 'Bisum_' + filenameextension
    saveEfBgsum = 'EfBgsum_' + filenameextension
    
    np.save(savepath+saveCRLBx,CRLBx)
    np.save(savepath+savemusum,musum)
    np.save(savepath+savesignalsum,signalsum)
    np.save(savepath+saveBisum,Bisum)
    np.save(savepath+saveEfBgsum,EfBgsum)
    
    return CRLBx, musum, signalsum, Bisum, EfBgsum

def simulate_xp(roiparams, pinholeparams, sigmas, theta, xp_range, scenario, N_m, illumination_doughnut=False):
    '''
    simulate_xp(roiparams, pinholeparams, sigmas, theta, xp_range, scenario) simulates
    the CRLB and underlying photon counts for a range of pinhole x-positions xp.
    
    Inputs:
        roiparams: [dx_i, Npixels]
            dx_i: Camera pixel size (m)
            Npixels: Amount of camera pixels per direction
        pinholeparams: (K x 3)-array [x_p, y_p, r_p]
            K: Amount of pinholes and patterns
            x_p: (K x 1) array of pinhole x-coordinates (m) (overwritten)
            y_p: (K x 1) array of pinhole y-coordinates (m)
            r_p: (K x 1) array of pinhole radii (m)
        sigmas: [sigmaillum, sigmapsf]
            sigmaillum: Standard deviation of Gaussian illumination PSF (m)
            sigmapsf: Standard deviation of Gaussian emission PSF (m)
        theta: [thetax, thetay, thetaI, thetab]
            thetax: Emitter x-position (m)
            thetay: Emitter y-position (m)
            thetaI: Expected signal photon budget under maximum illumination (photons)
            thetab: Expected background photon budget under maximum illumination (photons/(m*m))
        xp_range: range of xp values, overwrites x_p.
        scenario: Simulation scenario under consideration, which models the normalization
            constant (i.e. to benchmark against SMLM) and pattern-dependent or 
            pattern-independent background.
            Options:
                'fixed_budget': Scenario for which the signal photon budget is exhausted
                'fixed_energy': Scenario for which the illumination power and time are constant per pattern
                'fixed_budget_adaptedbg': Scenario for which the signal photon budget is exhausted, with pattern-independent background
        N_m: Amount of mesh pixels per direction, default set to 100.

    Optional input:
        illumination_doughnut: If True, doughnut-shaped intensity patterns with a
            zero-intensity minimum in the center will be used for illumination.
            If False, Gaussian illumination will be used.
        
    Outputs:
        CRLBx: Cramér-Rao lower bound (CRLB)
        musum: Illumination-intensity normalized total amount of photons in region of interest (ROI)
        signalsum: Illumination-intensity normalized total amount of signal photons in ROI
        Bisum: Sum of pattern-(in)dependent background constants over ROI
        EfBgsum: Illumination-intensity normalized total amount of background photons in ROI
    '''
    CRLBx = np.zeros(np.size(xp_range))
    musum = np.zeros(np.size(xp_range))
    signalsum = np.zeros(np.size(xp_range))
    Bisum = np.zeros(np.size(xp_range))
    EfBgsum = np.zeros(np.size(xp_range))
    
    for idx, x_p in enumerate(xp_range):
        pinholeparams[:,0] = x_p
        sf = SpinFlux(theta, sigmas, roiparams, pinholeparams, illumination_doughnut=illumination_doughnut, N_m = N_m, scenario=scenario)
        I, musum[idx], signalsum[idx], Bisum[idx], EfBgsum[idx] = sf.FisherMatrix(approximation='middle')
        I_inv = np.linalg.inv(I)
        CRLBx[idx] = np.sqrt(I_inv[0,0])
    return CRLBx, musum, signalsum, Bisum, EfBgsum

def simulate_yp(roiparams, pinholeparams, sigmas, theta, yp_range, scenario, N_m, illumination_doughnut=False):
    '''
    simulate_yp(roiparams, pinholeparams, sigmas, theta, yp_range, scenario) simulates
    the CRLB and underlying photon counts for a range of pinhole y-positions yp.
    
    Inputs:
        roiparams: [dx_i, Npixels]
            dx_i: Camera pixel size (m)
            Npixels: Amount of camera pixels per direction
        pinholeparams: (K x 3)-array [x_p, y_p, r_p]
            K: Amount of pinholes and patterns
            x_p: (K x 1) array of pinhole x-coordinates (m)
            y_p: (K x 1) array of pinhole y-coordinates (m) (overwritten)
            r_p: (K x 1) array of pinhole radii (m)
        sigmas: [sigmaillum, sigmapsf]
            sigmaillum: Standard deviation of Gaussian illumination PSF (m)
            sigmapsf: Standard deviation of Gaussian emission PSF (m)
        theta: [thetax, thetay, thetaI, thetab]
            thetax: Emitter x-position (m)
            thetay: Emitter y-position (m)
            thetaI: Expected signal photon budget under maximum illumination (photons)
            thetab: Expected background photon budget under maximum illumination (photons/(m*m))
        yp_range: range of yp values, overwrites y_p.
        scenario: Simulation scenario under consideration, which models the normalization
            constant (i.e. to benchmark against SMLM) and pattern-dependent or 
            pattern-independent background.
            Options:
                'fixed_budget': Scenario for which the signal photon budget is exhausted
                'fixed_energy': Scenario for which the illumination power and time are constant per pattern
                'fixed_budget_adaptedbg': Scenario for which the signal photon budget is exhausted, with pattern-independent background
        N_m: Amount of mesh pixels per direction, default set to 100.

    Optional input:
        illumination_doughnut: If True, doughnut-shaped intensity patterns with a
            zero-intensity minimum in the center will be used for illumination.
            If False, Gaussian illumination will be used.
        
    Outputs:
        CRLBx: Cramér-Rao lower bound (CRLB)
        musum: Illumination-intensity normalized total amount of photons in region of interest (ROI)
        signalsum: Illumination-intensity normalized total amount of signal photons in ROI
        Bisum: Sum of pattern-(in)dependent background constants over ROI
        EfBgsum: Illumination-intensity normalized total amount of background photons in ROI
    '''
    CRLBx = np.zeros(np.size(yp_range))
    musum = np.zeros(np.size(yp_range))
    signalsum = np.zeros(np.size(yp_range))
    Bisum = np.zeros(np.size(yp_range))
    EfBgsum = np.zeros(np.size(yp_range))
    
    for idx, y_p in enumerate(yp_range):
        pinholeparams[:,1] = y_p
        sf = SpinFlux(theta, sigmas, roiparams, pinholeparams, illumination_doughnut=illumination_doughnut, N_m = N_m, scenario=scenario)
        I, musum[idx], signalsum[idx], Bisum[idx], EfBgsum[idx] = sf.FisherMatrix(approximation='middle')
        I_inv = np.linalg.inv(I)
        CRLBx[idx] = np.sqrt(I_inv[0,0])
    return CRLBx, musum, signalsum, Bisum, EfBgsum

def simulate_xp_rp(roiparams, pinholeparams, sigmas, theta, r_p_list, scenario, illumination_doughnut=False, filenameextension = 'xp_rp', savepath = './SimData/', N_m = 100, start = 0, stop = 11*65*10**-9, num = 50):
    '''
    simulate_xp_rp(roiparams, pinholeparams, sigmas, theta, r_p_list, scenario) simulates
    the CRLB and underlying photon counts for a range of pinhole x-positions x_p
    and pinhole radii r_p.
    
    Inputs:
        roiparams: [dx_i, Npixels]
            dx_i: Camera pixel size (m)
            Npixels: Amount of camera pixels per direction
        pinholeparams: (K x 3)-array [x_p, y_p, r_p]
            K: Amount of pinholes and patterns
            x_p: (K x 1) array of pinhole x-coordinates (m)
            y_p: (K x 1) array of pinhole y-coordinates (m)
            r_p: (K x 1) array of pinhole radii (m) (overwritten)
        sigmas: [sigmaillum, sigmapsf]
            sigmaillum: Standard deviation of Gaussian illumination PSF (m)
            sigmapsf: Standard deviation of Gaussian emission PSF (m)
        theta: [thetax, thetay, thetaI, thetab]
            thetax: Emitter x-position (m)
            thetay: Emitter y-position (m)
            thetaI: Expected signal photon budget under maximum illumination (photons)
            thetab: Expected background photon budget under maximum illumination (photons/(m*m))
        r_p_list: range of r_p values, overwrites r_p.
        scenario: Simulation scenario under consideration, which models the normalization
            constant (i.e. to benchmark against SMLM) and pattern-dependent or 
            pattern-independent background.
            Options:
                'fixed_budget': Scenario for which the signal photon budget is exhausted
                'fixed_energy': Scenario for which the illumination power and time are constant per pattern
                'fixed_budget_adaptedbg': Scenario for which the signal photon budget is exhausted, with pattern-independent background
       
    Optional inputs:
        illumination_doughnut: If True, doughnut-shaped intensity patterns with a
            zero-intensity minimum in the center will be used for illumination.
            If False, Gaussian illumination will be used.
        filenameextension: string, describing simulation file to be saved to.
        savepath: string, describing simulation path to save to.
        N_m: Amount of mesh pixels per direction, default set to 100.
        start: starting order of linspace for pinhole x-positions x_p
        stop: stopping order of linspace for pinhole x-positions x_p
        num: amount of pinhole x-positions x_p to sample
        
    Outputs:
        CRLBx: Cramér-Rao lower bound (CRLB)
        musum: Illumination-intensity normalized total amount of photons in region of interest (ROI)
        signalsum: Illumination-intensity normalized total amount of signal photons in ROI
        Bisum: Sum of pattern-(in)dependent background constants over ROI
        EfBgsum: Illumination-intensity normalized total amount of background photons in ROI
    '''
    print('Starting simulation xp_rp')
    saverplist = 'rplist_' + filenameextension
    np.save(savepath+saverplist, r_p_list)
    
    xrange = np.linspace(start, stop, num)
    savexrange = 'xrange_' + filenameextension
    np.save(savepath+savexrange,xrange)
    
    CRLBx = np.zeros((np.size(r_p_list),np.size(xrange)))
    musum = np.zeros((np.size(r_p_list),np.size(xrange)))
    signalsum = np.zeros((np.size(r_p_list),np.size(xrange)))
    Bisum = np.zeros((np.size(r_p_list),np.size(xrange)))
    EfBgsum = np.zeros((np.size(r_p_list),np.size(xrange)))
    
    for i in range(np.size(r_p_list)):
        print('Simulation xp_rp: Evaluating CRLB for r_p with index ' + str(int(i)))
        r_p = r_p_list[i]
        pinholeparams[:,2] = r_p
        CRLBx[i,:], musum[i,:], signalsum[i,:], Bisum[i,:], EfBgsum[i,:] = simulate_xp(roiparams, pinholeparams, sigmas, theta, xrange, scenario=scenario, illumination_doughnut=illumination_doughnut, N_m = 100)
        
    saveCRLBx = 'CRLBx_' + filenameextension
    savemusum = 'musum_' + filenameextension
    savesignalsum = 'signalsum_' + filenameextension
    saveBisum = 'Bisum_' + filenameextension
    saveEfBgsum = 'EfBgsum_' + filenameextension
    
    np.save(savepath+saveCRLBx,CRLBx)
    np.save(savepath+savemusum,musum)
    np.save(savepath+savesignalsum,signalsum)
    np.save(savepath+saveBisum,Bisum)
    np.save(savepath+saveEfBgsum,EfBgsum)
    
    return CRLBx, musum, signalsum, Bisum, EfBgsum

def simulate_yp_rp(roiparams, pinholeparams, sigmas, theta, r_p_list, scenario, illumination_doughnut=False, filenameextension = 'yp_rp', savepath = './SimData/', N_m = 100, start = 0, stop = 11*65*10**-9, num = 50):
    '''
    simulate_yp_rp(roiparams, pinholeparams, sigmas, theta, r_p_list, scenario) simulates
    the CRLB and underlying photon counts for a range of pinhole y-positions y_p
    and pinhole radii r_p.
    
    Inputs:
        roiparams: [dx_i, Npixels]
            dx_i: Camera pixel size (m)
            Npixels: Amount of camera pixels per direction
        pinholeparams: (K x 3)-array [x_p, y_p, r_p]
            K: Amount of pinholes and patterns
            x_p: (K x 1) array of pinhole x-coordinates (m)
            y_p: (K x 1) array of pinhole y-coordinates (m)
            r_p: (K x 1) array of pinhole radii (m) (overwritten)
        sigmas: [sigmaillum, sigmapsf]
            sigmaillum: Standard deviation of Gaussian illumination PSF (m)
            sigmapsf: Standard deviation of Gaussian emission PSF (m)
        theta: [thetax, thetay, thetaI, thetab]
            thetax: Emitter x-position (m)
            thetay: Emitter y-position (m)
            thetaI: Expected signal photon budget under maximum illumination (photons)
            thetab: Expected background photon budget under maximum illumination (photons/(m*m))
        r_p_list: range of r_p values, overwrites r_p.
        scenario: Simulation scenario under consideration, which models the normalization
            constant (i.e. to benchmark against SMLM) and pattern-dependent or 
            pattern-independent background.
            Options:
                'fixed_budget': Scenario for which the signal photon budget is exhausted
                'fixed_energy': Scenario for which the illumination power and time are constant per pattern
                'fixed_budget_adaptedbg': Scenario for which the signal photon budget is exhausted, with pattern-independent background
       
    Optional inputs:
        illumination_doughnut: If True, doughnut-shaped intensity patterns with a
            zero-intensity minimum in the center will be used for illumination.
            If False, Gaussian illumination will be used.
        filenameextension: string, describing simulation file to be saved to.
        savepath: string, describing simulation path to save to.
        N_m: Amount of mesh pixels per direction, default set to 100.
        start: starting order of linspace for pinhole y-positions y_p
        stop: stopping order of linspace for pinhole y-positions y_p
        num: amount of pinhole y-positions y_p to sample
        
    Outputs:
        CRLBx: Cramér-Rao lower bound (CRLB)
        musum: Illumination-intensity normalized total amount of photons in region of interest (ROI)
        signalsum: Illumination-intensity normalized total amount of signal photons in ROI
        Bisum: Sum of pattern-(in)dependent background constants over ROI
        EfBgsum: Illumination-intensity normalized total amount of background photons in ROI
    '''
    print('Starting simulation yp_rp')
    saverplist = 'rplist_' + filenameextension
    np.save(savepath+saverplist, r_p_list)
    
    yrange = np.linspace(start, stop, num)
    saveyrange = 'yrange_' + filenameextension
    np.save(savepath+saveyrange,yrange)
    
    CRLBx = np.zeros((np.size(r_p_list),np.size(yrange)))
    musum = np.zeros((np.size(r_p_list),np.size(yrange)))
    signalsum = np.zeros((np.size(r_p_list),np.size(yrange)))
    Bisum = np.zeros((np.size(r_p_list),np.size(yrange)))
    EfBgsum = np.zeros((np.size(r_p_list),np.size(yrange)))
    
    for i in range(np.size(r_p_list)):
        print('Simulation yp_rp: Evaluating CRLB for r_p with index ' + str(int(i)))
        r_p = r_p_list[i]
        pinholeparams[:,2] = r_p
        CRLBx[i,:], musum[i,:], signalsum[i,:], Bisum[i,:], EfBgsum[i,:] = simulate_yp(roiparams, pinholeparams, sigmas, theta, yrange, scenario=scenario, illumination_doughnut=illumination_doughnut, N_m = 100)
        
    saveCRLBx = 'CRLBx_' + filenameextension
    savemusum = 'musum_' + filenameextension
    savesignalsum = 'signalsum_' + filenameextension
    saveBisum = 'Bisum_' + filenameextension
    saveEfBgsum = 'EfBgsum_' + filenameextension
    
    np.save(savepath+saveCRLBx,CRLBx)
    np.save(savepath+savemusum,musum)
    np.save(savepath+savesignalsum,signalsum)
    np.save(savepath+saveBisum,Bisum)
    np.save(savepath+saveEfBgsum,EfBgsum)
    
    return CRLBx, musum, signalsum, Bisum, EfBgsum

def simulate_xp_2pattern(roiparams, pinholeparams, sigmas, theta, xp_range, scenario, N_m, illumination_doughnut=False):
    '''
    simulate_xp_2pattern(roiparams, pinholeparams, sigmas, theta, xp_range, scenario) simulates
    the CRLB and underlying photon counts for the separated 2-pattern configuration,
    for a range of pinhole x-positions xp.
    
    Inputs:
        roiparams: [dx_i, Npixels]
            dx_i: Camera pixel size (m)
            Npixels: Amount of camera pixels per direction
        pinholeparams: (K x 3)-array [x_shift, y_p, r_p]
            K: Amount of pinholes and patterns
            x_shift: (K x 1) array of separated pinhole x-coordinates around x_p = 0 (m)
            y_p: (K x 1) array of pinhole y-coordinates (m)
            r_p: (K x 1) array of pinhole radii (m)
        sigmas: [sigmaillum, sigmapsf]
            sigmaillum: Standard deviation of Gaussian illumination PSF (m)
            sigmapsf: Standard deviation of Gaussian emission PSF (m)
        theta: [thetax, thetay, thetaI, thetab]
            thetax: Emitter x-position (m)
            thetay: Emitter y-position (m)
            thetaI: Expected signal photon budget under maximum illumination (photons)
            thetab: Expected background photon budget under maximum illumination (photons/(m*m))
        xp_range: range of xp values, added to x_p.
        scenario: Simulation scenario under consideration, which models the normalization
            constant (i.e. to benchmark against SMLM) and pattern-dependent or 
            pattern-independent background.
            Options:
                'fixed_budget': Scenario for which the signal photon budget is exhausted
                'fixed_energy': Scenario for which the illumination power and time are constant per pattern
                'fixed_budget_adaptedbg': Scenario for which the signal photon budget is exhausted, with pattern-independent background
        N_m: Amount of mesh pixels per direction, default set to 100.

    Optional input:
        illumination_doughnut: If True, doughnut-shaped intensity patterns with a
            zero-intensity minimum in the center will be used for illumination.
            If False, Gaussian illumination will be used.
        
    Outputs:
        CRLBx: Cramér-Rao lower bound (CRLB)
        musum: Illumination-intensity normalized total amount of photons in region of interest (ROI)
        signalsum: Illumination-intensity normalized total amount of signal photons in ROI
        Bisum: Sum of pattern-(in)dependent background constants over ROI
        EfBgsum: Illumination-intensity normalized total amount of background photons in ROI
    '''
    CRLBx = np.zeros(np.size(xp_range))
    musum = np.zeros(np.size(xp_range))
    signalsum = np.zeros(np.size(xp_range))
    Bisum = np.zeros(np.size(xp_range))
    EfBgsum = np.zeros(np.size(xp_range))
    
    for idx, x_p in enumerate(xp_range):
        pinholeparams[:,0] += x_p
        sf = SpinFlux(theta, sigmas, roiparams, pinholeparams, N_m = N_m, scenario=scenario, illumination_doughnut=illumination_doughnut)
        I, musum[idx], signalsum[idx], Bisum[idx], EfBgsum[idx] = sf.FisherMatrix(approximation='middle')
        I_inv = np.linalg.inv(I)
        CRLBx[idx] = np.sqrt(I_inv[0,0])
        pinholeparams[:,0] -= x_p
    return CRLBx, musum, signalsum, Bisum, EfBgsum

def simulate_yp_2pattern(roiparams, pinholeparams, sigmas, theta, xp_range, scenario, N_m, illumination_doughnut=False):
    '''
    simulate_yp_2pattern(roiparams, pinholeparams, sigmas, theta, yp_range, scenario) simulates
    the CRLB and underlying photon counts for the separated 2-pattern configuration,
    for a range of pinhole y-positions yp.
    
    Inputs:
        roiparams: [dx_i, Npixels]
            dx_i: Camera pixel size (m)
            Npixels: Amount of camera pixels per direction
        pinholeparams: (K x 3)-array [x_p, y_p, r_p]
            K: Amount of pinholes and patterns
            x_p: (K x 1) array of pinhole x-coordinates (m)
            y_shift: (K x 1) array of separated pinhole y-coordinates around y_p = 0 (m)
            r_p: (K x 1) array of pinhole radii (m)
        sigmas: [sigmaillum, sigmapsf]
            sigmaillum: Standard deviation of Gaussian illumination PSF (m)
            sigmapsf: Standard deviation of Gaussian emission PSF (m)
        theta: [thetax, thetay, thetaI, thetab]
            thetax: Emitter x-position (m)
            thetay: Emitter y-position (m)
            thetaI: Expected signal photon budget under maximum illumination (photons)
            thetab: Expected background photon budget under maximum illumination (photons/(m*m))
        yp_range: range of yp values, added to y_p.
        scenario: Simulation scenario under consideration, which models the normalization
            constant (i.e. to benchmark against SMLM) and pattern-dependent or 
            pattern-independent background.
            Options:
                'fixed_budget': Scenario for which the signal photon budget is exhausted
                'fixed_energy': Scenario for which the illumination power and time are constant per pattern
                'fixed_budget_adaptedbg': Scenario for which the signal photon budget is exhausted, with pattern-independent background
        N_m: Amount of mesh pixels per direction, default set to 100.

    Optional input:
        illumination_doughnut: If True, doughnut-shaped intensity patterns with a
            zero-intensity minimum in the center will be used for illumination.
            If False, Gaussian illumination will be used.
        
    Outputs:
        CRLBx: Cramér-Rao lower bound (CRLB)
        musum: Illumination-intensity normalized total amount of photons in region of interest (ROI)
        signalsum: Illumination-intensity normalized total amount of signal photons in ROI
        Bisum: Sum of pattern-(in)dependent background constants over ROI
        EfBgsum: Illumination-intensity normalized total amount of background photons in ROI
    '''
    CRLBx = np.zeros(np.size(xp_range))
    musum = np.zeros(np.size(xp_range))
    signalsum = np.zeros(np.size(xp_range))
    Bisum = np.zeros(np.size(xp_range))
    EfBgsum = np.zeros(np.size(xp_range))
    
    for idx, x_p in enumerate(xp_range):
        pinholeparams[:,1] += x_p
        sf = SpinFlux(theta, sigmas, roiparams, pinholeparams, scenario, N_m = N_m, scenario=scenario, illumination_doughnut=illumination_doughnut)
        I, musum[idx], signalsum[idx], Bisum[idx], EfBgsum[idx] = sf.FisherMatrix(approximation='middle')
        I_inv = np.linalg.inv(I)
        CRLBx[idx] = np.sqrt(I_inv[0,0])
        pinholeparams[:,1] -= x_p
    return CRLBx, musum, signalsum, Bisum, EfBgsum

def simulate_xp_s_2pattern(roiparams, pinholeparams, sigmas, theta, s_list, scenario, illumination_doughnut=False, filenameextension = 'xp_s_2pattern', savepath = './SimData/', N_m = 100, start = 0, stop = 11*65*10**-9, num = 50):
    '''
    simulate_xp_s_patterns(roiparams, pinholeparams, sigmas, theta, s_list, scenario) simulates
    the CRLB and underlying photon counts in the x-separated 2-pattern configuration,
    for a range of focus x-positions x_p and separations s.
    
    Inputs:
        roiparams: [dx_i, Npixels]
            dx_i: Camera pixel size (m)
            Npixels: Amount of camera pixels per direction
        pinholeparams: (K x 3)-array [x_p, y_p, r_p]
            K: Amount of pinholes and patterns
            x_p: (K x 1) array of pinhole x-coordinates (m) (overwritten)
            y_p: (K x 1) array of pinhole y-coordinates (m)
            r_p: (K x 1) array of pinhole radii (m)
        sigmas: [sigmaillum, sigmapsf]
            sigmaillum: Standard deviation of Gaussian illumination PSF (m)
            sigmapsf: Standard deviation of Gaussian emission PSF (m)
        theta: [thetax, thetay, thetaI, thetab]
            thetax: Emitter x-position (m)
            thetay: Emitter y-position (m)
            thetaI: Expected signal photon budget under maximum illumination (photons)
            thetab: Expected background photon budget under maximum illumination (photons/(m*m))
        s_list: range of s values, overwrites x_p.
        scenario: Simulation scenario under consideration, which models the normalization
            constant (i.e. to benchmark against SMLM) and pattern-dependent or 
            pattern-independent background.
            Options:
                'fixed_budget': Scenario for which the signal photon budget is exhausted
                'fixed_energy': Scenario for which the illumination power and time are constant per pattern
                'fixed_budget_adaptedbg': Scenario for which the signal photon budget is exhausted, with pattern-independent background
       
    Optional inputs:
        illumination_doughnut: If True, doughnut-shaped intensity patterns with a
            zero-intensity minimum in the center will be used for illumination.
            If False, Gaussian illumination will be used.
        filenameextension: string, describing simulation file to be saved to.
        savepath: string, describing simulation path to save to.
        N_m: Amount of mesh pixels per direction, default set to 100.
        start: starting order of linspace for pinhole x-positions x_p
        stop: stopping order of linspace for pinhole x-positions x_p
        num: amount of focus x-positions x_p to sample
        
    Outputs:
        CRLBx: Cramér-Rao lower bound (CRLB)
        musum: Illumination-intensity normalized total amount of photons in region of interest (ROI)
        signalsum: Illumination-intensity normalized total amount of signal photons in ROI
        Bisum: Sum of pattern-(in)dependent background constants over ROI
        EfBgsum: Illumination-intensity normalized total amount of background photons in ROI
    '''
    print('Starting simulation xp_s_2pattern')
    saveslist = 'slist_' + filenameextension
    np.save(savepath+saveslist, s_list)
    
    xrange = np.linspace(start, stop, num)
    savexrange = 'xrange_' + filenameextension
    np.save(savepath+savexrange,xrange)
    
    CRLBx = np.zeros((np.size(s_list),np.size(xrange)))
    musum = np.zeros((np.size(s_list),np.size(xrange)))
    signalsum = np.zeros((np.size(s_list),np.size(xrange)))
    Bisum = np.zeros((np.size(s_list),np.size(xrange)))
    EfBgsum = np.zeros((np.size(s_list),np.size(xrange)))
    
    for i in range(np.size(s_list)):
        print('Simulation xp_s_2pattern: Evaluating CRLB for s with index ' + str(int(i)))
        s = s_list[i]
        pinholeparams[0,0] = -s/2
        pinholeparams[1,0] = s/2
        CRLBx[i,:], musum[i,:], signalsum[i,:], Bisum[i,:], EfBgsum[i,:] = simulate_xp_2pattern(roiparams, pinholeparams, sigmas, theta, xrange, scenario=scenario, illumination_doughnut=illumination_doughnut, N_m = 100)
        
    saveCRLBx = 'CRLBx_' + filenameextension
    savemusum = 'musum_' + filenameextension
    savesignalsum = 'signalsum_' + filenameextension
    saveBisum = 'Bisum_' + filenameextension
    saveEfBgsum = 'EfBgsum_' + filenameextension
    
    np.save(savepath+saveCRLBx,CRLBx)
    np.save(savepath+savemusum,musum)
    np.save(savepath+savesignalsum,signalsum)
    np.save(savepath+saveBisum,Bisum)
    np.save(savepath+saveEfBgsum,EfBgsum)
    
    return CRLBx, musum, signalsum, Bisum, EfBgsum

def simulate_yp_s_2pattern(roiparams, pinholeparams, sigmas, theta, s_list, scenario, illumination_doughnut=False, filenameextension = 'yp_s_2pattern', savepath = './SimData/', N_m = 100, start = 0, stop = 11*65*10**-9, num = 50):
    '''
    simulate_yp_s_patterns(roiparams, pinholeparams, sigmas, theta, s_list, scenario) simulates
    the CRLB and underlying photon counts in the y-separated 2-pattern configuration,
    for a range of focus y-positions y_p and separations s.
    
    Inputs:
        roiparams: [dx_i, Npixels]
            dx_i: Camera pixel size (m)
            Npixels: Amount of camera pixels per direction
        pinholeparams: (K x 3)-array [x_p, y_p, r_p]
            K: Amount of pinholes and patterns
            x_p: (K x 1) array of pinhole x-coordinates (m)
            y_p: (K x 1) array of pinhole y-coordinates (m) (overwritten)
            r_p: (K x 1) array of pinhole radii (m)
        sigmas: [sigmaillum, sigmapsf]
            sigmaillum: Standard deviation of Gaussian illumination PSF (m)
            sigmapsf: Standard deviation of Gaussian emission PSF (m)
        theta: [thetax, thetay, thetaI, thetab]
            thetax: Emitter x-position (m)
            thetay: Emitter y-position (m)
            thetaI: Expected signal photon budget under maximum illumination (photons)
            thetab: Expected background photon budget under maximum illumination (photons/(m*m))
        s_list: range of s values, overwrites y_p.
        scenario: Simulation scenario under consideration, which models the normalization
            constant (i.e. to benchmark against SMLM) and pattern-dependent or 
            pattern-independent background.
            Options:
                'fixed_budget': Scenario for which the signal photon budget is exhausted
                'fixed_energy': Scenario for which the illumination power and time are constant per pattern
                'fixed_budget_adaptedbg': Scenario for which the signal photon budget is exhausted, with pattern-independent background
       
    Optional inputs:
        illumination_doughnut: If True, doughnut-shaped intensity patterns with a
            zero-intensity minimum in the center will be used for illumination.
            If False, Gaussian illumination will be used.
        filenameextension: string, describing simulation file to be saved to.
        savepath: string, describing simulation path to save to.
        N_m: Amount of mesh pixels per direction, default set to 100.
        start: starting order of linspace for focus y-positions y_p
        stop: stopping order of linspace for focus y-positions y_p
        num: amount of focus y-positions y_p to sample
        
    Outputs:
        CRLBx: Cramér-Rao lower bound (CRLB)
        musum: Illumination-intensity normalized total amount of photons in region of interest (ROI)
        signalsum: Illumination-intensity normalized total amount of signal photons in ROI
        Bisum: Sum of pattern-(in)dependent background constants over ROI
        EfBgsum: Illumination-intensity normalized total amount of background photons in ROI
    '''
    print('Starting simulation yp_s_2pattern')
    saveslist = 'slist_' + filenameextension
    np.save(savepath+saveslist, s_list)
    
    xrange = np.linspace(start, stop, num)
    savexrange = 'xrange_' + filenameextension
    np.save(savepath+savexrange,xrange)
    
    CRLBx = np.zeros((np.size(s_list),np.size(xrange)))
    musum = np.zeros((np.size(s_list),np.size(xrange)))
    signalsum = np.zeros((np.size(s_list),np.size(xrange)))
    Bisum = np.zeros((np.size(s_list),np.size(xrange)))
    EfBgsum = np.zeros((np.size(s_list),np.size(xrange)))
    
    for i in range(np.size(s_list)):
        print('Simulation yp_s_2pattern: Evaluating CRLB for s with index ' + str(int(i)))
        s = s_list[i]
        pinholeparams[0,1] = -s/2
        pinholeparams[1,1] = s/2
        CRLBx[i,:], musum[i,:], signalsum[i,:], Bisum[i,:], EfBgsum[i,:] = simulate_yp_2pattern(roiparams, pinholeparams, sigmas, theta, xrange, scenario=scenario, illumination_doughnut=illumination_doughnut, N_m = 100)
        
    saveCRLBx = 'CRLBx_' + filenameextension
    savemusum = 'musum_' + filenameextension
    savesignalsum = 'signalsum_' + filenameextension
    saveBisum = 'Bisum_' + filenameextension
    saveEfBgsum = 'EfBgsum_' + filenameextension
    
    np.save(savepath+saveCRLBx,CRLBx)
    np.save(savepath+savemusum,musum)
    np.save(savepath+savesignalsum,signalsum)
    np.save(savepath+saveBisum,Bisum)
    np.save(savepath+saveEfBgsum,EfBgsum)
    
    return CRLBx, musum, signalsum, Bisum, EfBgsum

def simulate_yp_xs_2pattern(roiparams, pinholeparams, sigmas, theta, s_list, scenario, illumination_doughnut=False, filenameextension = 'yp_xs_2pattern', savepath = './SimData/', N_m = 100, start = 0, stop = 11*65*10**-9, num = 50):
    '''
    simulate_yp_xs_patterns(roiparams, pinholeparams, sigmas, theta, s_list, scenario) simulates
    the CRLB and underlying photon counts in the x-separated 2-pattern configuration,
    for a range of focus y-positions y_p and separations s.
    
    Inputs:
        roiparams: [dx_i, Npixels]
            dx_i: Camera pixel size (m)
            Npixels: Amount of camera pixels per direction
        pinholeparams: (K x 3)-array [x_p, y_p, r_p]
            K: Amount of pinholes and patterns
            x_p: (K x 1) array of pinhole x-coordinates (m)
            y_p: (K x 1) array of pinhole y-coordinates (m)
            r_p: (K x 1) array of pinhole radii (m)
        sigmas: [sigmaillum, sigmapsf]
            sigmaillum: Standard deviation of Gaussian illumination PSF (m)
            sigmapsf: Standard deviation of Gaussian emission PSF (m)
        theta: [thetax, thetay, thetaI, thetab]
            thetax: Emitter x-position (m)
            thetay: Emitter y-position (m)
            thetaI: Expected signal photon budget under maximum illumination (photons)
            thetab: Expected background photon budget under maximum illumination (photons/(m*m))
        s_list: range of s values, added to x_p.
        scenario: Simulation scenario under consideration, which models the normalization
            constant (i.e. to benchmark against SMLM) and pattern-dependent or 
            pattern-independent background.
            Options:
                'fixed_budget': Scenario for which the signal photon budget is exhausted
                'fixed_energy': Scenario for which the illumination power and time are constant per pattern
                'fixed_budget_adaptedbg': Scenario for which the signal photon budget is exhausted, with pattern-independent background
       
    Optional inputs:
        illumination_doughnut: If True, doughnut-shaped intensity patterns with a
            zero-intensity minimum in the center will be used for illumination.
            If False, Gaussian illumination will be used.
        filenameextension: string, describing simulation file to be saved to.
        savepath: string, describing simulation path to save to.
        N_m: Amount of mesh pixels per direction, default set to 100.
        start: starting order of linspace for focus y-positions y_p
        stop: stopping order of linspace for focus y-positions y_p
        num: amount of focus y-positions y_p to sample
        
    Outputs:
        CRLBx: Cramér-Rao lower bound (CRLB)
        musum: Illumination-intensity normalized total amount of photons in region of interest (ROI)
        signalsum: Illumination-intensity normalized total amount of signal photons in ROI
        Bisum: Sum of pattern-(in)dependent background constants over ROI
        EfBgsum: Illumination-intensity normalized total amount of background photons in ROI
    '''
    print('Starting simulation yp_xs_2pattern')
    saveslist = 'slist_' + filenameextension
    np.save(savepath+saveslist, s_list)
    
    xrange = np.linspace(start, stop, num)
    savexrange = 'xrange_' + filenameextension
    np.save(savepath+savexrange,xrange)
    
    CRLBx = np.zeros((np.size(s_list),np.size(xrange)))
    musum = np.zeros((np.size(s_list),np.size(xrange)))
    signalsum = np.zeros((np.size(s_list),np.size(xrange)))
    Bisum = np.zeros((np.size(s_list),np.size(xrange)))
    EfBgsum = np.zeros((np.size(s_list),np.size(xrange)))
    
    for i in range(np.size(s_list)):
        print('Simulation yp_xs_2pattern: Evaluating CRLB for s with index ' + str(int(i)))
        s = s_list[i]
        pinholeparams[0,0] -= s/2
        pinholeparams[1,0] += s/2
        CRLBx[i,:], musum[i,:], signalsum[i,:], Bisum[i,:], EfBgsum[i,:] = simulate_yp_2pattern(roiparams, pinholeparams, sigmas, theta, xrange, scenario=scenario, illumination_doughnut=illumination_doughnut, N_m = 100)
        pinholeparams[0,0] += s/2
        pinholeparams[1,0] -= s/2
        
    saveCRLBx = 'CRLBx_' + filenameextension
    savemusum = 'musum_' + filenameextension
    savesignalsum = 'signalsum_' + filenameextension
    saveBisum = 'Bisum_' + filenameextension
    saveEfBgsum = 'EfBgsum_' + filenameextension
    
    np.save(savepath+saveCRLBx,CRLBx)
    np.save(savepath+savemusum,musum)
    np.save(savepath+savesignalsum,signalsum)
    np.save(savepath+saveBisum,Bisum)
    np.save(savepath+saveEfBgsum,EfBgsum)
    
    return CRLBx, musum, signalsum, Bisum, EfBgsum

def coordinates_triangulation(center, r, rotation, return_centercoordinates = True):
    '''
    coordinates_triangulation(center, r, rotation) returns pinhole center coordinates
    of the triangle configuration
    
    Inputs:
        center: [x_f, y_f]
            x_f: Pattern focus x-coordinate
            y_f: Pattern focus y-coordinate
        r: Pattern spacing (distance center of triangle to corner of triangle)
        rotation: Counterclockwise rotation angle of triangle (in rad)
    
    Optional input:
        return_centercoordinates: If True, center pinhole coordinates are appended to the output
    
    Output:
        (3 x 2) (or (4 x 2) if return_centercoordinates is set to True) array of
        pinhole (x,y) coordinates.
    '''
    coordinates = np.array([[r,0], [-r*np.cos(np.pi/3), r*np.sin(np.pi/3)], [-r*np.cos(np.pi/3), -r*np.sin(np.pi/3)]])
    coordinates = coordinates @ np.array([[np.cos(rotation), np.sin(rotation)], [-np.sin(rotation), np.cos(rotation)]])
    coordinates += center
    if return_centercoordinates:
        coordinates = np.concatenate((coordinates, np.array([center])))
    return coordinates  

def simulate_xp_r_triangulation(roiparams, pinholeparams, sigmas, theta, r_list, y_center, scenario, rotation = 0, pinhole_in_center = True, illumination_doughnut=False, filenameextension = 'xp_r_triangulation', savepath = './SimData/', N_m = 100, start = 0, stop = 11*65*10**-9, num = 50):
    '''
    simulate_xp_r_triangulation(roiparams, pinholeparams, sigmas, theta, r_list, scenario, rotation) simulates
    the CRLB and underlying photon counts in the equilateral triangle pattern configuration,
    for a range of focus x-positions x_p and spacings r.
    
    Inputs:
        roiparams: [dx_i, Npixels]
            dx_i: Camera pixel size (m)
            Npixels: Amount of camera pixels per direction
        pinholeparams: (K x 3)-array [x_p, y_p, r_p]
            K: Amount of pinholes and patterns
            x_p: (K x 1) array of pinhole x-coordinates (m) (overwritten)
            y_p: (K x 1) array of pinhole y-coordinates (m) (overwritten)
            r_p: (K x 1) array of pinhole radii (m)
        sigmas: [sigmaillum, sigmapsf]
            sigmaillum: Standard deviation of Gaussian illumination PSF (m)
            sigmapsf: Standard deviation of Gaussian emission PSF (m)
        theta: [thetax, thetay, thetaI, thetab]
            thetax: Emitter x-position (m)
            thetay: Emitter y-position (m)
            thetaI: Expected signal photon budget under maximum illumination (photons)
            thetab: Expected background photon budget under maximum illumination (photons/(m*m))
        r_list: range of r values
        y_center: y-coordinate of pattern focus
        scenario: Simulation scenario under consideration, which models the normalization
            constant (i.e. to benchmark against SMLM) and pattern-dependent or 
            pattern-independent background.
            Options:
                'fixed_budget': Scenario for which the signal photon budget is exhausted
                'fixed_energy': Scenario for which the illumination power and time are constant per pattern
                'fixed_budget_adaptedbg': Scenario for which the signal photon budget is exhausted, with pattern-independent background
        rotation: Counterclockwise rotation angle of triangle (in rad)
       
    Optional inputs:
        pinhole_in_center: boolean, set to True to add pinhole in center of triangle configuration
        illumination_doughnut: If True, doughnut-shaped intensity patterns with a
            zero-intensity minimum in the center will be used for illumination.
            If False, Gaussian illumination will be used.
        filenameextension: string, describing simulation file to be saved to.
        savepath: string, describing simulation path to save to.
        N_m: Amount of mesh pixels per direction, default set to 100.
        start: starting order of linspace for focus x-positions x_p
        stop: stopping order of linspace for focus x-positions x_p
        num: amount of focus x-positions x_p to sample
        
    Outputs:
        CRLBx: Cramér-Rao lower bound (CRLB)
        musum: Illumination-intensity normalized total amount of photons in region of interest (ROI)
        signalsum: Illumination-intensity normalized total amount of signal photons in ROI
        Bisum: Sum of pattern-(in)dependent background constants over ROI
        EfBgsum: Illumination-intensity normalized total amount of background photons in ROI
    '''
    
    print('Starting simulation xp_r_triangulation')
    saverlist = 'rlist_' + filenameextension
    np.save(savepath+saverlist, r_list)
    
    xrange = np.linspace(start, stop, num)
    savexrange = 'xrange_' + filenameextension
    np.save(savepath+savexrange,xrange)
    
    CRLBx = np.zeros((np.size(r_list),np.size(xrange)))
    musum = np.zeros((np.size(r_list),np.size(xrange)))
    signalsum = np.zeros((np.size(r_list),np.size(xrange)))
    Bisum = np.zeros((np.size(r_list),np.size(xrange)))
    EfBgsum = np.zeros((np.size(r_list),np.size(xrange)))
    
    for i in range(np.size(r_list)):
        print('Simulation xp_r_triangulation: Evaluating CRLB for r with index ' + str(int(i)))
        r = r_list[i]
        for idx, x in enumerate(xrange):
            pinholeparams[:,[0,1]] = coordinates_triangulation([x,y_center], r, rotation, return_centercoordinates = pinhole_in_center)
            sf = SpinFlux(theta, sigmas, roiparams, pinholeparams, illumination_doughnut=illumination_doughnut, N_m = N_m, scenario=scenario)
            I, musum[i, idx], signalsum[i,idx], Bisum[i,idx], EfBgsum[i,idx] = sf.FisherMatrix(approximation='middle')
            I_inv = np.linalg.inv(I)
            CRLBx[i,idx] = np.sqrt(I_inv[0,0])
        
    saveCRLBx = 'CRLBx_' + filenameextension
    savemusum = 'musum_' + filenameextension
    savesignalsum = 'signalsum_' + filenameextension
    saveBisum = 'Bisum_' + filenameextension
    saveEfBgsum = 'EfBgsum_' + filenameextension
    
    np.save(savepath+saveCRLBx,CRLBx)
    np.save(savepath+savemusum,musum)
    np.save(savepath+savesignalsum,signalsum)
    np.save(savepath+saveBisum,Bisum)
    np.save(savepath+saveEfBgsum,EfBgsum)
    
    return CRLBx, musum, signalsum, Bisum, EfBgsum

def simulate_thetaI_Nm_12pattern(roiparams, pinholeparams, sigmas, theta, thetaI_list, scenario, illumination_doughnut=False, filenameextension = 'thetaI_Nm_12pattern', savepath = './SimData/', start = 10, stop = 200, num = 20):
    '''
    simulate_thetaI_Nm_triangulation(roiparams, pinholeparams, sigmas, theta, thetaI_list, scenario) simulates
    the CRLB and underlying photon counts in the single- or two-pattern configuration,
    for a range of mesh pixel amounts Nm and signal photon counts thetaI
    
    Inputs:
        roiparams: [dx_i, Npixels]
            dx_i: Camera pixel size (m)
            Npixels: Amount of camera pixels per direction
        pinholeparams: (K x 3)-array [x_p, y_p, r_p]
            K: Amount of pinholes and patterns
            x_p: (K x 1) array of pinhole x-coordinates (m)
            y_p: (K x 1) array of pinhole y-coordinates (m)
            r_p: (K x 1) array of pinhole radii (m)
        sigmas: [sigmaillum, sigmapsf]
            sigmaillum: Standard deviation of Gaussian illumination PSF (m)
            sigmapsf: Standard deviation of Gaussian emission PSF (m)
        theta: [thetax, thetay, thetaI, thetab]
            thetax: Emitter x-position (m)
            thetay: Emitter y-position (m)
            thetaI: Expected signal photon budget under maximum illumination (photons) (overwritten)
            thetab: Expected background photon budget under maximum illumination (photons/(m*m))
        thetaI_list: range of thetaI values
        scenario: Simulation scenario under consideration, which models the normalization
            constant (i.e. to benchmark against SMLM) and pattern-dependent or 
            pattern-independent background.
            Options:
                'fixed_budget': Scenario for which the signal photon budget is exhausted
                'fixed_energy': Scenario for which the illumination power and time are constant per pattern
                'fixed_budget_adaptedbg': Scenario for which the signal photon budget is exhausted, with pattern-independent background
       
    Optional inputs:
        pinhole_in_center: boolean, set to True to add pinhole in center of triangle configuration
        illumination_doughnut: If True, doughnut-shaped intensity patterns with a
            zero-intensity minimum in the center will be used for illumination.
            If False, Gaussian illumination will be used.
        filenameextension: string, describing simulation file to be saved to.
        savepath: string, describing simulation path to save to.
        start: starting order of linspace for mesh pixel amount Nm
        stop: stopping order of linspace for mesh pixel amount Nm
        num: amount of mesh pixel amounts Nm to sample
        
    Outputs:
        CRLBx: Cramér-Rao lower bound (CRLB)
        musum: Illumination-intensity normalized total amount of photons in region of interest (ROI)
        signalsum: Illumination-intensity normalized total amount of signal photons in ROI
        Bisum: Sum of pattern-(in)dependent background constants over ROI
        EfBgsum: Illumination-intensity normalized total amount of background photons in ROI
    '''
    print('Starting simulation thetaI_Nm_12pattern')
    saverlist = 'thetaIlist_' + filenameextension
    np.save(savepath+saverlist, thetaI_list)
    
    Nmrange = np.linspace(start, stop, num, dtype=int)
    savexrange = 'Nmrange_' + filenameextension
    np.save(savepath+savexrange,Nmrange)
    
    CRLBx = np.zeros((np.size(thetaI_list),np.size(Nmrange)))
    musum = np.zeros((np.size(thetaI_list),np.size(Nmrange)))
    signalsum = np.zeros((np.size(thetaI_list),np.size(Nmrange)))
    Bisum = np.zeros((np.size(thetaI_list),np.size(Nmrange)))
    EfBgsum = np.zeros((np.size(thetaI_list),np.size(Nmrange)))
    
    for i in range(np.size(thetaI_list)):
        print('Simulation thetaI_Nm_12pattern: Evaluating CRLB for theta_I with index ' + str(int(i)))
        theta[2] = thetaI_list[i]
        for idx, Nm in enumerate(Nmrange):
            sf = SpinFlux(theta, sigmas, roiparams, pinholeparams, illumination_doughnut=illumination_doughnut, N_m = Nm, scenario=scenario)
            I, musum[i, idx], signalsum[i,idx], Bisum[i,idx], EfBgsum[i,idx] = sf.FisherMatrix(approximation='middle')
            I_inv = np.linalg.inv(I)
            CRLBx[i,idx] = np.sqrt(I_inv[0,0])
        
    saveCRLBx = 'CRLBx_' + filenameextension
    savemusum = 'musum_' + filenameextension
    savesignalsum = 'signalsum_' + filenameextension
    saveBisum = 'Bisum_' + filenameextension
    saveEfBgsum = 'EfBgsum_' + filenameextension
    
    np.save(savepath+saveCRLBx,CRLBx)
    np.save(savepath+savemusum,musum)
    np.save(savepath+savesignalsum,signalsum)
    np.save(savepath+saveBisum,Bisum)
    np.save(savepath+saveEfBgsum,EfBgsum)
    
    return CRLBx, musum, signalsum, Bisum, EfBgsum

def simulate_thetaI_Nm_triangulation(roiparams, pinholeparams, sigmas, theta, thetaI_list, r, scenario, rotation, pinhole_in_center = True, illumination_doughnut=False, filenameextension = 'thetaI_Nm_triangulation', savepath = './SimData/', start = 10, stop = 200, num = 20):
    '''
    simulate_thetaI_Nm_triangulation(roiparams, pinholeparams, sigmas, theta, thetaI_list, r, scenario, rotation) simulates
    the CRLB and underlying photon counts in the equilateral triangle pattern configuration,
    for a range of mesh pixel amounts Nm and signal photon counts thetaI
    
    Inputs:
        roiparams: [dx_i, Npixels]
            dx_i: Camera pixel size (m)
            Npixels: Amount of camera pixels per direction
        pinholeparams: (K x 3)-array [x_p, y_p, r_p]
            K: Amount of pinholes and patterns
            x_p: (K x 1) array of pinhole x-coordinates (m) (overwritten)
            y_p: (K x 1) array of pinhole y-coordinates (m) (overwritten)
            r_p: (K x 1) array of pinhole radii (m)
        sigmas: [sigmaillum, sigmapsf]
            sigmaillum: Standard deviation of Gaussian illumination PSF (m)
            sigmapsf: Standard deviation of Gaussian emission PSF (m)
        theta: [thetax, thetay, thetaI, thetab]
            thetax: Emitter x-position (m)
            thetay: Emitter y-position (m)
            thetaI: Expected signal photon budget under maximum illumination (photons) (overwritten)
            thetab: Expected background photon budget under maximum illumination (photons/(m*m))
        thetaI_list: range of thetaI values
        r: Pattern spacing (distance center of triangle to corner of triangle)
        scenario: Simulation scenario under consideration, which models the normalization
            constant (i.e. to benchmark against SMLM) and pattern-dependent or 
            pattern-independent background.
            Options:
                'fixed_budget': Scenario for which the signal photon budget is exhausted
                'fixed_energy': Scenario for which the illumination power and time are constant per pattern
                'fixed_budget_adaptedbg': Scenario for which the signal photon budget is exhausted, with pattern-independent background
        rotation: Counterclockwise rotation angle of triangle (in rad)
       
    Optional inputs:
        pinhole_in_center: boolean, set to True to add pinhole in center of triangle configuration
        illumination_doughnut: If True, doughnut-shaped intensity patterns with a
            zero-intensity minimum in the center will be used for illumination.
            If False, Gaussian illumination will be used.
        filenameextension: string, describing simulation file to be saved to.
        savepath: string, describing simulation path to save to.
        start: starting order of linspace for mesh pixel amount Nm
        stop: stopping order of linspace for mesh pixel amount Nm
        num: amount of mesh pixel amounts Nm to sample
        
    Outputs:
        CRLBx: Cramér-Rao lower bound (CRLB)
        musum: Illumination-intensity normalized total amount of photons in region of interest (ROI)
        signalsum: Illumination-intensity normalized total amount of signal photons in ROI
        Bisum: Sum of pattern-(in)dependent background constants over ROI
        EfBgsum: Illumination-intensity normalized total amount of background photons in ROI
    '''
    print('Starting simulation thetaI_Nm_triangulation')
    saverlist = 'thetaIlist_' + filenameextension
    np.save(savepath+saverlist, thetaI_list)
    
    Nmrange = np.linspace(start, stop, num, dtype=int)
    savexrange = 'Nmrange_' + filenameextension
    np.save(savepath+savexrange,Nmrange)
    
    CRLBx = np.zeros((np.size(thetaI_list),np.size(Nmrange)))
    musum = np.zeros((np.size(thetaI_list),np.size(Nmrange)))
    signalsum = np.zeros((np.size(thetaI_list),np.size(Nmrange)))
    Bisum = np.zeros((np.size(thetaI_list),np.size(Nmrange)))
    EfBgsum = np.zeros((np.size(thetaI_list),np.size(Nmrange)))
    
    for i in range(np.size(thetaI_list)):
        print('Simulation thetaI_Nm_triangulation: Evaluating CRLB for theta_I with index ' + str(int(i)))
        theta[2] = thetaI_list[i]
        for idx, Nm in enumerate(Nmrange):
            pinholeparams[:,[0,1]] = coordinates_triangulation([theta[0],theta[1]], r, rotation, return_centercoordinates = pinhole_in_center)
            sf = SpinFlux(theta, sigmas, roiparams, pinholeparams, illumination_doughnut=illumination_doughnut, N_m = Nm, scenario=scenario)
            I, musum[i, idx], signalsum[i,idx], Bisum[i,idx], EfBgsum[i,idx] = sf.FisherMatrix(approximation='middle')
            I_inv = np.linalg.inv(I)
            CRLBx[i,idx] = np.sqrt(I_inv[0,0])
        
    saveCRLBx = 'CRLBx_' + filenameextension
    savemusum = 'musum_' + filenameextension
    savesignalsum = 'signalsum_' + filenameextension
    saveBisum = 'Bisum_' + filenameextension
    saveEfBgsum = 'EfBgsum_' + filenameextension
    
    np.save(savepath+saveCRLBx,CRLBx)
    np.save(savepath+savemusum,musum)
    np.save(savepath+savesignalsum,signalsum)
    np.save(savepath+saveBisum,Bisum)
    np.save(savepath+saveEfBgsum,EfBgsum)
    
    return CRLBx, musum, signalsum, Bisum, EfBgsum

def simulate_SMLM(roiparams, sigmapsf, theta, filenameextension = 'SMLM', savepath = './SimData/'):
    '''
    simulate_SMLM(roiparams, sigmapsf, theta) simulates the CRLB and underlying 
    photon counts of the SMLM benchmark.
    
    Inputs:
        roiparams: [dx_i, Npixels]
            dx_i: Camera pixel size (m)
            Npixels: Amount of camera pixels per direction
        sigmapsf: Standard deviation of Gaussian emission PSF (m)
        theta: [thetax, thetay, thetaI, thetab]
            thetax: Emitter x-position (m)
            thetay: Emitter y-position (m)
            thetaI: Expected signal photon budget under maximum illumination (photons)
            thetab: Expected background photon budget under maximum illumination (photons/(m*m))
       
    Optional inputs:
        filenameextension: string, describing simulation file to be saved to.
        savepath: string, describing simulation path to save to.
        
    Outputs:
        CRLBx: Cramér-Rao lower bound (CRLB)
        musum: Illumination-intensity normalized total amount of photons in region of interest (ROI)
        signalsum: Illumination-intensity normalized total amount of signal photons in ROI
        Bisum: Sum of pattern-(in)dependent background constants over ROI
        EfBgsum: Illumination-intensity normalized total amount of background photons in ROI
    '''
    print('Starting simulation SMLM')
    smlm = SMLM(theta, sigmapsf, roiparams)
    I, musum, signalsum, Bisum, EfBgsum = smlm.FisherMatrix()
    I_inv = np.linalg.inv(I)
    CRLBx = np.sqrt(I_inv[0,0]) 
    
    saveCRLBx = 'CRLBx_' + filenameextension
    savemusum = 'musum_' + filenameextension
    savesignalsum = 'signalsum_' + filenameextension
    saveBisum = 'Bisum_' + filenameextension
    saveEfBgsum = 'EfBgsum_' + filenameextension
    
    np.save(savepath+saveCRLBx,CRLBx)
    np.save(savepath+savemusum,musum)
    np.save(savepath+savesignalsum,signalsum)
    np.save(savepath+saveBisum,Bisum)
    np.save(savepath+saveEfBgsum,EfBgsum)
    
    return CRLBx, musum, signalsum, Bisum, EfBgsum

def simulate(plotnumber):
    '''
    simulate(plotnumber) simulates the SpinFlux Cramér-Rao lower bound (CRLB) data 
    underlying the figure specified in plotnumber.
    
    Input:
        plotnumber: string, describing the figure of which CRLB is to be simulated.
            Options:
                'SMLM': SMLM Benchmark
                '3-S4', '3-heatmap': Figure 3 & S4: 1 x-offset pattern, fixed signal photon budget
                '4de-S6', '4bc-heatmap': Figure 4b-e & S6: 2 x-separated patterns, fixed signal photon budget
                '4ij-S10', '4gh-heatmap': Figure 4g-j & S10: 3 patterns in equilateral triangle configuration without center pinhole, fixed signal photon budget
                '5-S16', '5-heatmap': Figure 5 & S16: 4 doughnut-patterns in equilateral triangle configuration with center pinhole, fixed signal photon budget
                'S2': Figure S2: Effect of pinhole discretization on CRLB.
                'S5': Figure S5: 1 y-offset pattern, fixed signal photon budget
                'S7': Figure S7: 2 x-separated patterns with y-offset, fixed signal photon budget
                'S8': Figure S8: 2 x-separated patterns without pinholes, fixed signal photon budget
                'S9': Figure S9: 2 y-separated patterns, fixed signal photon budget
                'S11': Figure S11: 3 patterns in 90 degree rotated equilateral triangle configuration without center pinhole, fixed signal photon budget
                'S12': Figure S12: 4 patterns in equilateral triangle configuration with center pinhole, fixed signal photon budget
                'S13': Figure S13: 4 patterns in 90 degree rotated equilateral triangle configuration with center pinhole, fixed signal photon budget
                'S14': Figure S14: 3 doughnut-patterns in equilateral triangle configuration without center pinhole, fixed signal photon budget
                'S15': Figure S15: 3 doughnut-patterns in 90 degree rotated equilateral triangle configuration without center pinhole, fixed signal photon budget
                'S17': Figure S17: 4 doughnut-patterns in 90 degree rotated equilateral triangle configuration with center pinhole, fixed signal photon budget                
                'S18': Figure S18: 1 x-offset pattern, fixed illumination power and time
                'S19': Figure S19: 1 y-offset pattern, fixed illumination power and time
                'S20': Figure S20: 2 x-separated patterns, fixed illumination power and time
                'S21': Figure S21: 2 x-separated patterns with y-offset, fixed illumination power and time
                'S22': Figure S22: 2 x-separated patterns without pinholes, fixed illumination power and time
                'S23': Figure S23: 2 y-separated patterns, fixed illumination power and time
                'S24': Figure S24: 3 patterns in equilateral triangle configuration without center pinhole, fixed illumination power and time
                'S25': Figure S25: 3 patterns in 90 degree rotated equilateral triangle configuration without center pinhole, fixed illumination power and time
                'S26': Figure S26: 4 patterns in equilateral triangle configuration with center pinhole, fixed illumination power and time
                'S27': Figure S27: 4 patterns in 90 degree rotated equilateral triangle configuration with center pinhole, fixed illumination power and time
                'S28': Figure S28: 3 doughnut-patterns in equilateral triangle configuration without center pinhole, fixed illumination power and time
                'S29': Figure S29: 3 doughnut-patterns in 90 degree rotated equilateral triangle configuration without center pinhole, fixed illumination power and time
                'S30': Figure S30: 4 doughnut-patterns in equilateral triangle configuration with center pinhole, fixed fixed illumination power and time
                'S31': Figure S31: 4 doughnut-patterns in 90 degree rotated equilateral triangle configuration with center pinhole, fixed illumination power and time
                'S32': Figure S32: 1 x-offset pattern, fixed signal photon budget with pattern-independent background
                'S33': Figure S33: 1 y-offset pattern, fixed signal photon budget with pattern-independent background
                'S34': Figure S34: 2 x-separated patterns, fixed signal photon budget with pattern-independent background
                'S35': Figure S35: 2 x-separated patterns with y-offset, fixed signal photon budget with pattern-independent background
                'S36': Figure S36: 2 x-separated patterns without pinholes, fixed signal photon budget with pattern-independent background
                'S37': Figure S37: 2 y-separated patterns, fixed signal photon budget with pattern-independent background
                'S38': Figure S38: 3 patterns in equilateral triangle configuration without center pinhole, fixed signal photon budget with pattern-independent background
                'S39': Figure S39: 3 patterns in 90 degree rotated equilateral triangle configuration without center pinhole, fixed signal photon budget with pattern-independent background
                'S40': Figure S40: 4 patterns in equilateral triangle configuration with center pinhole, fixed signal photon budget with pattern-independent background
                'S41': Figure S41: 4 patterns in 90 degree rotated equilateral triangle configuration with center pinhole, fixed signal photon budget with pattern-independent background
                'S42': Figure S42: 3 doughnut-patterns in equilateral triangle configuration without center pinhole, fixed signal photon budget with pattern-independent background
                'S43': Figure S43: 3 doughnut-patterns in 90 degree rotated equilateral triangle configuration without center pinhole, fixed signal photon budget with pattern-independent background
                'S44': Figure S44: 4 doughnut-patterns in equilateral triangle configuration with center pinhole, fixed signal photon budget with pattern-independent background
                'S45': Figure S45: 4 doughnut-patterns in 90 degree rotated equilateral triangle configuration with center pinhole, fixed signal photon budget with pattern-independent background  
    '''
    #Parameters
    L_illum = 546*10**-9        #Illumination wavelength (m) (cf. Sirinakis et al. (2022))
    L_emission = 600*10**(-9)   #Emission wavelength (m) (cf. Sirinakis et al. (2022))
    NA = 1.35                   #Numerical aperture (cf. Sirinakis et al. (2022))
    
    dx_i = 65*10**-9    #Camera pixel size (m)
    Npixels = 10        #Amount of camera pixels in each direction of the region of interest
    
    x_p = Npixels*dx_i/2            #Pinhole x-coordinate (m)
    y_p = Npixels*dx_i/2            #Pinhole y-coordinate (m)
    r_p = 0.21*L_emission/NA        #Radius pinhole (m)
    
    roiparams=[dx_i, Npixels]
    sigmas = [0.21*L_illum/NA, 0.21*L_emission/NA]              #Illumination and emission PSF standard deviation (m) (cf. Huang and Olivo-Marin (2014))
    theta = [Npixels*dx_i/2, Npixels*dx_i/2, 2000, 8/(dx_i**2)] #Emitter x- and y-position (m), expected signal photon budget (photons) and expected background (photons/m**2)
    
    N_m = 100           #Amount of mesh pixels per direction
    nsteps = 50         #Amount of intermediate steps when simulating spatial effect on CRLB
    
    if plotnumber=='SMLM':
        #SMLM Benchmark
        simulate_SMLM(roiparams, sigmas[1], theta, filenameextension = 'SMLM')
    
    elif plotnumber=='3-S4':
        #Figures 3 & S4: 1 x-offset pattern, fixed signal photon budget
        r_p_list = np.array([1, 2, 3, 4, 100000])*r_p
        pinholeparams=np.repeat([[x_p, y_p, r_p]], 1, axis=0)
        simulate_xp_rp(roiparams, pinholeparams, sigmas, theta, r_p_list, scenario='fixed_budget', filenameextension = 'xp_rp_1pattern_fb', N_m = N_m, start = 3*dx_i, stop = 7*dx_i, num = nsteps)
    
        thetab_list = np.array([0, 1, 4, 8, 16])/(dx_i**2)
        pinholeparams=np.repeat([[x_p, y_p, 3*r_p]], 1, axis=0)
        simulate_thetaI_thetab(roiparams, pinholeparams, sigmas, theta, thetab_list, scenario='fixed_budget', filenameextension = 'thetaI_thetab_1pattern_fb', N_m = N_m, start = 2, stop = 4, num = 20)        

    elif plotnumber=='3-heatmap':
        #Figure 3 (heatmap): 1 x-offset pattern shifted through y, fixed signal photon budget
        r_p_list = np.array([3])*r_p
        xnum = 51
        ynum = 51
        ypos = np.linspace(3*dx_i, 7*dx_i, ynum)
        for yindex in range(ynum):
            pinholeparams=np.repeat([[x_p, ypos[yindex], r_p]], 1, axis=0)
            simulate_xp_rp(roiparams, pinholeparams, sigmas, theta, r_p_list, scenario='fixed_budget', filenameextension = ('xp_rp_1pattern_heatmap_' + str(int(yindex)) + '_fb'), N_m = N_m, start = 3*dx_i, stop = 7*dx_i, num = xnum)
        
    elif plotnumber=='4de-S6':
        #Figures 4de & S6: 2 x-separated patterns, fixed signal photon budget
        s_list = np.array([0, 1, 2, 3, 4, 5])*r_p
        pinholeparams=np.repeat([[x_p, y_p, 3*r_p]], 2, axis=0)
        simulate_xp_s_2pattern(roiparams, pinholeparams, sigmas, theta, s_list, scenario='fixed_budget', filenameextension = 'xp_s_2pattern_fb', N_m = N_m, start = 3*dx_i, stop = 7*dx_i, num = nsteps)

        thetab_list = np.array([0, 1, 4, 8, 16])/(dx_i**2)
        pinholeparams=np.array([[x_p-4*r_p/2, y_p, 3*r_p], [x_p+4*r_p/2, y_p, 3*r_p]])
        simulate_thetaI_thetab(roiparams, pinholeparams, sigmas, theta, thetab_list, scenario='fixed_budget', filenameextension = 'thetaI_thetab_2pattern_fb', N_m = N_m, start = 2, stop = 4, num = 20)

    elif plotnumber=='4bc-heatmap':
        #Figure 4bc (heatmap): 2 x-separated patterns shifted through y, fixed signal photon budget
        s_list = np.array([4])*r_p
        xnum = 51
        ynum = 51
        ypos = np.linspace(3*dx_i, 7*dx_i, ynum)
        for yindex in range(ynum):
            pinholeparams=np.repeat([[x_p, ypos[yindex], 3*r_p]], 2, axis=0)   
            simulate_xp_s_2pattern(roiparams, pinholeparams, sigmas, theta, s_list, scenario='fixed_budget', filenameextension = ('xp_s_2pattern_heatmap_' + str(int(yindex)) + '_fb'), N_m = N_m, start = 3*dx_i, stop = 7*dx_i, num = xnum)

    elif plotnumber=='4ij-S10':
        #Figures 4ij & S10: 3 patterns in equilateral triangle configuration without center pinhole, fixed signal photon budget
        r_list = np.array([1, 2, 3, 4, 5])*r_p
        pinholeparams=np.repeat([[0, 0, 3*r_p]], 3, axis=0)
        simulate_xp_r_triangulation(roiparams, pinholeparams, sigmas, theta, r_list, y_p, scenario='fixed_budget', rotation = np.pi/2, pinhole_in_center = False, filenameextension = 'xp_r_triangulation_rot_nocenter_ext_fb', N_m = N_m, start = 3*dx_i, stop = 7*dx_i, num = nsteps)
        
        thetab_list = np.array([0, 1, 4, 8, 16])/(dx_i**2)
        pinholeparams=np.repeat([[0, 0, 3*r_p]], 3, axis=0)
        pinholeparams[:,[0,1]]= coordinates_triangulation([x_p, y_p], 2*r_p, np.pi/2, return_centercoordinates = False)
        simulate_thetaI_thetab(roiparams, pinholeparams, sigmas, theta, thetab_list, scenario='fixed_budget', filenameextension = 'thetaI_thetab_triangulation_rot_nocenter_fb', N_m = N_m, start = 2, stop = 5, num = 20)

    elif plotnumber=='4gh-heatmap':
        #Figure 4gh (heatmap): 3 patterns in equilateral triangle configuration without center pinhole, shifted through y, fixed signal photon budget
        r_list = np.array([2])*r_p
        xnum = 41
        ynum = 41
        ypos = np.linspace(3*dx_i, 7*dx_i, ynum)
        pinholeparams=np.repeat([[0, 0, 3*r_p]], 3, axis=0)
        for yindex in range(ynum):
            simulate_xp_r_triangulation(roiparams, pinholeparams, sigmas, theta, r_list, ypos[yindex], scenario='fixed_budget', rotation = np.pi/2, pinhole_in_center = False, filenameextension = ('xp_r_triangulation_rot_nocenter_ext_heatmap_' + str(int(yindex)) + '_fb'), N_m = N_m, start = 3*dx_i, stop = 7*dx_i, num = xnum)

    elif plotnumber=='5-S16':
        #Figures 5 & S16: 4 doughnut-patterns in equilateral triangle configuration with center pinhole, fixed signal photon budget
        r_list = np.array([0.5, 1, 2, 3, 4])*r_p
        pinholeparams=np.repeat([[0, 0, 3*r_p]], 4, axis=0)
        simulate_xp_r_triangulation(roiparams, pinholeparams, sigmas, theta, r_list, y_p, scenario='fixed_budget', rotation = np.pi/2, pinhole_in_center = True, illumination_doughnut = True, filenameextension = 'xp_r_triangulation_rot_doughnut_ext_fb', N_m = N_m, start = 3*dx_i, stop = 7*dx_i, num = nsteps)
        
        thetab_list = np.array([0, 1, 4, 8, 16])/(dx_i**2)
        pinholeparams=np.repeat([[0, 0, 3*r_p]], 4, axis=0)
        pinholeparams[:,[0,1]]= coordinates_triangulation([x_p, y_p], 2*r_p, np.pi/2, return_centercoordinates = True)
        simulate_thetaI_thetab(roiparams, pinholeparams, sigmas, theta, thetab_list, scenario='fixed_budget', illumination_doughnut = True, filenameextension = 'thetaI_thetab_triangulation_rot_doughnut_fb', N_m = N_m, start = 2, stop = 5, num = 20)


    elif plotnumber=='5-heatmap':
        #Figure 5 (heatmap): 4 doughnut-patterns in equilateral triangle configuration with center pinhole, fixed signal photon budget
        r_list = np.array([3])*r_p
        xnum = 35
        ynum = 35
        ypos = np.linspace(3*dx_i, 7*dx_i, ynum)
        pinholeparams=np.repeat([[0, 0, 3*r_p]], 4, axis=0)
        for yindex in range(ynum):
            simulate_xp_r_triangulation(roiparams, pinholeparams, sigmas, theta, r_list, ypos[yindex], scenario='fixed_budget', rotation = np.pi/2, pinhole_in_center = True, illumination_doughnut = True, filenameextension = ('xp_r_triangulation_rot_doughnut_ext_heatmap_' + str(int(yindex)) + '_fb'), N_m = N_m, start = 3*dx_i, stop = 7*dx_i, num = xnum)
            
    elif plotnumber=='S2':
        #Supplementary Figure 2: Effect of pinhole discretization on CRLB.
        
        #Fixed signal photon budget
        thetaI_list = np.array([500, 1000, 2000, 5000])
        pinholeparams=np.repeat([[x_p, y_p, 3*r_p]], 1, axis=0)
        simulate_thetaI_Nm_12pattern(roiparams, pinholeparams, sigmas, theta, thetaI_list,scenario='fixed_budget',  filenameextension = 'thetaI_Nm_1pattern_fb', savepath = './SimData/', start = 10, stop = 200, num = 20)

        thetaI_list = np.array([500, 1000, 2000, 5000])
        pinholeparams=np.array([[x_p-4*r_p/2, y_p, 3*r_p], [x_p+4*r_p/2, y_p, 3*r_p]])
        simulate_thetaI_Nm_12pattern(roiparams, pinholeparams, sigmas, theta, thetaI_list, scenario='fixed_budget', filenameextension = 'thetaI_Nm_2pattern_fb', savepath = './SimData/', start = 10, stop = 200, num = 20)
        
        thetaI_list = np.array([500, 1000, 2000, 5000])
        pinholeparams=np.repeat([[0, 0, 3*r_p]], 3, axis=0)
        simulate_thetaI_Nm_triangulation(roiparams, pinholeparams, sigmas, theta, thetaI_list, 2*r_p, scenario='fixed_budget', rotation = np.pi/2, pinhole_in_center = False, filenameextension = 'thetaI_Nm_triangulation_rot_nocenter_fb', savepath = './SimData/', start = 10, stop = 200, num = 20)
        
        thetaI_list = np.array([500, 1000, 2000, 5000])
        pinholeparams=np.repeat([[x_p, y_p, 3*r_p]], 1, axis=0)
        simulate_thetaI_Nm_12pattern(roiparams, pinholeparams, sigmas, theta, thetaI_list, scenario='fixed_budget', filenameextension = 'thetaI_Nm_1pattern_gt_fb', savepath = './SimData/', start = 1000, stop = 1000, num = 1)
        
        thetaI_list = np.array([500, 1000, 2000, 5000])
        pinholeparams=np.array([[x_p-4*r_p/2, y_p, 3*r_p], [x_p+4*r_p/2, y_p, 3*r_p]])
        simulate_thetaI_Nm_12pattern(roiparams, pinholeparams, sigmas, theta, thetaI_list, scenario='fixed_budget', filenameextension = 'thetaI_Nm_2pattern_gt_fb', savepath = './SimData/', start = 1000, stop = 1000, num = 1)
        
        thetaI_list = np.array([500, 1000, 2000, 5000])
        pinholeparams=np.repeat([[0, 0, 3*r_p]], 3, axis=0)
        simulate_thetaI_Nm_triangulation(roiparams, pinholeparams, sigmas, theta, thetaI_list, 2*r_p, scenario='fixed_budget', rotation = np.pi/2, pinhole_in_center = False, filenameextension = 'thetaI_Nm_triangulation_rot_nocenter_gt_fb', savepath = './SimData/', start = 1000, stop = 1000, num = 1)
        
        #Fixed illumination energy
        thetaI_list = np.array([500, 1000, 2000, 5000])
        pinholeparams=np.repeat([[x_p, y_p, 3*r_p]], 1, axis=0)
        simulate_thetaI_Nm_12pattern(roiparams, pinholeparams, sigmas, theta, thetaI_list, scenario='fixed_energy', filenameextension = 'thetaI_Nm_1pattern', savepath = './SimData/', start = 10, stop = 200, num = 20)
        
        thetaI_list = np.array([500, 1000, 2000, 5000])
        pinholeparams=np.array([[x_p-4*r_p/2, y_p, 3*r_p], [x_p+4*r_p/2, y_p, 3*r_p]])
        simulate_thetaI_Nm_12pattern(roiparams, pinholeparams, sigmas, theta, thetaI_list, scenario='fixed_energy', filenameextension = 'thetaI_Nm_2pattern', savepath = './SimData/', start = 10, stop = 200, num = 20)
        
        thetaI_list = np.array([500, 1000, 2000, 5000])
        pinholeparams=np.repeat([[0, 0, 3*r_p]], 3, axis=0)
        simulate_thetaI_Nm_triangulation(roiparams, pinholeparams, sigmas, theta, thetaI_list, 2*r_p, scenario='fixed_energy', rotation = np.pi/2, pinhole_in_center = False, filenameextension = 'thetaI_Nm_triangulation_rot_nocenter', savepath = './SimData/', start = 10, stop = 200, num = 20)
        
        thetaI_list = np.array([500, 1000, 2000, 5000])
        pinholeparams=np.repeat([[x_p, y_p, 3*r_p]], 1, axis=0)
        simulate_thetaI_Nm_12pattern(roiparams, pinholeparams, sigmas, theta, thetaI_list, scenario='fixed_energy', filenameextension = 'thetaI_Nm_1pattern_gt', savepath = './SimData/', start = 1000, stop = 1000, num = 1)
        
        thetaI_list = np.array([500, 1000, 2000, 5000])
        pinholeparams=np.array([[x_p-4*r_p/2, y_p, 3*r_p], [x_p+4*r_p/2, y_p, 3*r_p]])
        simulate_thetaI_Nm_12pattern(roiparams, pinholeparams, sigmas, theta, thetaI_list, scenario='fixed_energy', filenameextension = 'thetaI_Nm_2pattern_gt', savepath = './SimData/', start = 1000, stop = 1000, num = 1)
        
        thetaI_list = np.array([500, 1000, 2000, 5000])
        pinholeparams=np.repeat([[0, 0, 3*r_p]], 3, axis=0)
        simulate_thetaI_Nm_triangulation(roiparams, pinholeparams, sigmas, theta, thetaI_list, 2*r_p, rotation = np.pi/2, scenario='fixed_energy', pinhole_in_center = False, filenameextension = 'thetaI_Nm_triangulation_rot_nocenter_gt', savepath = './SimData/', start = 1000, stop = 1000, num = 1)       

        #Fixed signal photon budget, pattern-independent background
        thetaI_list = np.array([500, 1000, 2000, 5000])
        pinholeparams=np.repeat([[x_p, y_p, 3*r_p]], 1, axis=0)
        simulate_thetaI_Nm_12pattern(roiparams, pinholeparams, sigmas, theta, thetaI_list, scenario='fixed_budget_adaptedbg', filenameextension = 'thetaI_Nm_1pattern_adaptedbg', savepath = './SimData/', start = 10, stop = 200, num = 20)

        thetaI_list = np.array([500, 1000, 2000, 5000])
        pinholeparams=np.array([[x_p-4*r_p/2, y_p, 3*r_p], [x_p+4*r_p/2, y_p, 3*r_p]])
        simulate_thetaI_Nm_12pattern(roiparams, pinholeparams, sigmas, theta, thetaI_list, scenario='fixed_budget_adaptedbg', filenameextension = 'thetaI_Nm_2pattern_adaptedbg', savepath = './SimData/', start = 10, stop = 200, num = 20)

        thetaI_list = np.array([500, 1000, 2000, 5000])
        pinholeparams=np.repeat([[0, 0, 3*r_p]], 3, axis=0)
        simulate_thetaI_Nm_triangulation(roiparams, pinholeparams, sigmas, theta, thetaI_list, 2*r_p, scenario='fixed_budget_adaptedbg', rotation = np.pi/2, pinhole_in_center = False, filenameextension = 'thetaI_Nm_triangulation_rot_nocenter_adaptedbg', savepath = './SimData/', start = 10, stop = 200, num = 20)

        thetaI_list = np.array([500, 1000, 2000, 5000])
        pinholeparams=np.repeat([[x_p, y_p, 3*r_p]], 1, axis=0)
        simulate_thetaI_Nm_12pattern(roiparams, pinholeparams, sigmas, theta, thetaI_list, scenario='fixed_budget_adaptedbg', filenameextension = 'thetaI_Nm_1pattern_gt_adaptedbg', savepath = './SimData/', start = 1000, stop = 1000, num = 1)

        thetaI_list = np.array([500, 1000, 2000, 5000])
        pinholeparams=np.array([[x_p-4*r_p/2, y_p, 3*r_p], [x_p+4*r_p/2, y_p, 3*r_p]])
        simulate_thetaI_Nm_12pattern(roiparams, pinholeparams, sigmas, theta, thetaI_list, scenario='fixed_budget_adaptedbg', filenameextension = 'thetaI_Nm_2pattern_gt_adaptedbg', savepath = './SimData/', start = 1000, stop = 1000, num = 1)

        thetaI_list = np.array([500, 1000, 2000, 5000])
        pinholeparams=np.repeat([[0, 0, 3*r_p]], 3, axis=0)
        simulate_thetaI_Nm_triangulation(roiparams, pinholeparams, sigmas, theta, thetaI_list, 2*r_p, scenario='fixed_budget_adaptedbg', rotation = np.pi/2, pinhole_in_center = False, filenameextension = 'thetaI_Nm_triangulation_rot_nocenter_gt_adaptedbg', savepath = './SimData/', start = 1000, stop = 1000, num = 1)

    elif plotnumber=='S5':
        #Supplementary Figure 5: 1 y-offset pattern, fixed signal photon budget
        r_p_list = np.array([1, 2, 3, 4, 100000])*r_p
        pinholeparams=np.repeat([[x_p, y_p, r_p]], 1, axis=0)
        simulate_yp_rp(roiparams, pinholeparams, sigmas, theta, r_p_list, scenario='fixed_budget', filenameextension = 'yp_rp_1pattern_fb', N_m = N_m, start = 3*dx_i, stop = 7*dx_i, num = nsteps)

    elif plotnumber=='S7':
        #Supplementary Figure 7: 2 x-separated patterns with y-offset, fixed signal photon budget
        s_list = np.array([0, 1, 2, 3, 4, 5])*r_p
        pinholeparams=np.repeat([[x_p, 0, 3*r_p]], 2, axis=0)
        simulate_yp_xs_2pattern(roiparams, pinholeparams, sigmas, theta, s_list, scenario='fixed_budget', filenameextension = 'yp_xs_2pattern_fb', N_m = N_m, start = 3*dx_i, stop = 7*dx_i, num = nsteps)        

    elif plotnumber=='S8':
        #Supplementary Figure 8: 2 x-separated patterns without pinholes, fixed signal photon 
        s_list = np.array([0, 1, 2, 3, 4, 5])*r_p
        pinholeparams=np.repeat([[x_p, y_p, 10000*r_p]], 2, axis=0)
        simulate_xp_s_2pattern(roiparams, pinholeparams, sigmas, theta, s_list, scenario='fixed_budget', filenameextension = 'xp_s_2pattern_nopinhole_fb', N_m = N_m, start = 3*dx_i, stop = 7*dx_i, num = nsteps)

        thetab_list = np.array([0, 1, 4, 8, 16])/(dx_i**2)
        pinholeparams=np.array([[x_p-4*r_p/2, y_p, 1000*r_p], [x_p+4*r_p/2, y_p, 1000*r_p]])
        simulate_thetaI_thetab(roiparams, pinholeparams, sigmas, theta, thetab_list, scenario='fixed_budget', filenameextension = 'thetaI_thetab_2pattern_nopinhole_fb', N_m = N_m, start = 2, stop = 4, num = 20)
                        
    elif plotnumber=='S9':
        #Supplementary Figure 9: 2 y-separated patterns, fixed signal photon budget
        s_list = np.array([0, 1, 2, 3, 4, 5])*r_p
        pinholeparams=np.repeat([[x_p, y_p, 3*r_p]], 2, axis=0)
        simulate_yp_s_2pattern(roiparams, pinholeparams, sigmas, theta, s_list, scenario='fixed_budget', filenameextension = 'yp_s_2pattern_fb', N_m = N_m, start = 3*dx_i, stop = 7*dx_i, num = nsteps)
        
        thetab_list = np.array([0, 1, 4, 8, 16])/(dx_i**2)
        pinholeparams=np.array([[x_p, y_p-4*r_p/2, 3*r_p], [x_p, y_p+4*r_p/2, 3*r_p]])
        simulate_thetaI_thetab(roiparams, pinholeparams, sigmas, theta, thetab_list, scenario='fixed_budget', filenameextension = 'thetaI_thetab_2pattern_rot_fb', N_m = N_m, start = 2, stop = 4, num = 20)            

    elif plotnumber=='S11':
        #Supplementary Figure 11: 3 patterns in 90 degree rotated equilateral triangle configuration without center pinhole, fixed signal photon budget
        r_list = np.array([1, 2, 3, 4, 5])*r_p
        pinholeparams=np.repeat([[0, 0, 3*r_p]], 3, axis=0)
        simulate_xp_r_triangulation(roiparams, pinholeparams, sigmas, theta, r_list, y_p, scenario='fixed_budget', rotation = 0, pinhole_in_center = False, filenameextension = 'xp_r_triangulation_nocenter_ext_fb', N_m = N_m, start = 3*dx_i, stop = 7*dx_i, num = nsteps)        
        
        thetab_list = np.array([0, 1, 4, 8, 16])/(dx_i**2)
        pinholeparams=np.repeat([[0, 0, 3*r_p]], 3, axis=0)
        pinholeparams[:,[0,1]]= coordinates_triangulation([x_p, y_p], 2*r_p, 0, return_centercoordinates = False)
        simulate_thetaI_thetab(roiparams, pinholeparams, sigmas, theta, thetab_list, scenario='fixed_budget', filenameextension = 'thetaI_thetab_triangulation_nocenter_fb', N_m = N_m, start = 2, stop = 4, num = 20)    

    elif plotnumber=='S12':
        #Supplementary Figure 12: 4 patterns in equilateral triangle configuration with center pinhole, fixed signal photon budget
        r_list = np.array([0.5, 1, 1.5, 2])*r_p
        pinholeparams=np.repeat([[0, 0, 3*r_p]], 4, axis=0)
        simulate_xp_r_triangulation(roiparams, pinholeparams, sigmas, theta, r_list, y_p, scenario='fixed_budget', rotation = np.pi/2, pinhole_in_center = True, filenameextension = 'xp_r_triangulation_rot_fb', N_m = N_m, start = 3*dx_i, stop = 7*dx_i, num = nsteps)
        
        thetab_list = np.array([0, 1, 4, 8, 16])/(dx_i**2)
        pinholeparams=np.repeat([[0, 0, 3*r_p]], 4, axis=0)
        pinholeparams[:,[0,1]]= coordinates_triangulation([x_p, y_p], 2*r_p, np.pi/2, return_centercoordinates = True)
        simulate_thetaI_thetab(roiparams, pinholeparams, sigmas, theta, thetab_list, scenario='fixed_budget', filenameextension = 'thetaI_thetab_triangulation_rot_fb', N_m = N_m, start = 2, stop = 5, num = 20)

    elif plotnumber=='S13':
        #Supplementary Figure 13: 4 patterns in 90 degree rotated equilateral triangle configuration with center pinhole, fixed signal photon budget
        r_list = np.array([0.5, 1, 1.5, 2])*r_p
        pinholeparams=np.repeat([[0, 0, 3*r_p]], 4, axis=0)
        simulate_xp_r_triangulation(roiparams, pinholeparams, sigmas, theta, r_list, y_p, scenario='fixed_budget', rotation = 0, pinhole_in_center = True, filenameextension = 'xp_r_triangulation_fb', N_m = N_m, start = 3*dx_i, stop = 7*dx_i, num = nsteps)

        thetab_list = np.array([0, 1, 4, 8, 16])/(dx_i**2)
        pinholeparams=np.repeat([[0, 0, 3*r_p]], 4, axis=0)
        pinholeparams[:,[0,1]]= coordinates_triangulation([x_p, y_p], 2*r_p, 0, return_centercoordinates = True)
        simulate_thetaI_thetab(roiparams, pinholeparams, sigmas, theta, thetab_list, scenario='fixed_budget', filenameextension = 'thetaI_thetab_triangulation_fb', N_m = N_m, start = 2, stop = 4, num = 20)

    elif plotnumber=='S14':
        #Supplementary Figure 14: 3 doughnut-patterns in equilateral triangle configuration without center pinhole, fixed signal photon budget
        r_list = np.array([1, 2, 3, 4, 5])*r_p
        pinholeparams=np.repeat([[0, 0, 3*r_p]], 3, axis=0)
        simulate_xp_r_triangulation(roiparams, pinholeparams, sigmas, theta, r_list, y_p, scenario='fixed_budget', rotation = np.pi/2, pinhole_in_center = False, illumination_doughnut = True, filenameextension = 'xp_r_triangulation_rot_nocenter_doughnut_fb', N_m = N_m, start = 3*dx_i, stop = 7*dx_i, num = nsteps)
        
        thetab_list = np.array([0, 1, 4, 8, 16])/(dx_i**2)
        pinholeparams=np.repeat([[0, 0, 3*r_p]], 3, axis=0)
        pinholeparams[:,[0,1]]= coordinates_triangulation([x_p, y_p], 2*r_p, np.pi/2, return_centercoordinates = False)
        simulate_thetaI_thetab(roiparams, pinholeparams, sigmas, theta, thetab_list, scenario='fixed_budget', illumination_doughnut = True, filenameextension = 'thetaI_thetab_triangulation_rot_nocenter_doughnut_fb', N_m = N_m, start = 2, stop = 5, num = 20)

    elif plotnumber=='S15':
        #Supplementary Figure 15: 3 doughnut-patterns in 90 degree rotated equilateral triangle configuration without center pinhole, fixed signal photon budget
        r_list = np.array([1, 2, 3, 4, 5])*r_p
        pinholeparams=np.repeat([[0, 0, 3*r_p]], 3, axis=0)
        simulate_xp_r_triangulation(roiparams, pinholeparams, sigmas, theta, r_list, y_p, scenario='fixed_budget', rotation = 0, pinhole_in_center = False, illumination_doughnut = True, filenameextension = 'xp_r_triangulation_nocenter_doughnut_fb', N_m = N_m, start = 3*dx_i, stop = 7*dx_i, num = nsteps)        
        
        thetab_list = np.array([0, 1, 4, 8, 16])/(dx_i**2)
        pinholeparams=np.repeat([[0, 0, 3*r_p]], 3, axis=0)
        pinholeparams[:,[0,1]]= coordinates_triangulation([x_p, y_p], 2*r_p, 0, return_centercoordinates = False)
        simulate_thetaI_thetab(roiparams, pinholeparams, sigmas, theta, thetab_list, scenario='fixed_budget', illumination_doughnut = True, filenameextension = 'thetaI_thetab_triangulation_nocenter_doughnut_fb', N_m = N_m, start = 2, stop = 4, num = 20)    


    elif plotnumber=='S17':
        #Supplementary Figure 17: 4 doughnut-patterns in 90 degree rotated equilateral triangle configuration with center pinhole, fixed signal photon budget
        r_list = np.array([0.5, 1, 2, 3, 4])*r_p
        pinholeparams=np.repeat([[0, 0, 3*r_p]], 4, axis=0)
        simulate_xp_r_triangulation(roiparams, pinholeparams, sigmas, theta, r_list, y_p, scenario='fixed_budget', rotation = 0, pinhole_in_center = True, illumination_doughnut = True, filenameextension = 'xp_r_triangulation_doughnut_ext_fb', N_m = N_m, start = 3*dx_i, stop = 7*dx_i, num = nsteps)

        thetab_list = np.array([0, 1, 4, 8, 16])/(dx_i**2)
        pinholeparams=np.repeat([[0, 0, 3*r_p]], 4, axis=0)
        pinholeparams[:,[0,1]]= coordinates_triangulation([x_p, y_p], 2*r_p, 0, return_centercoordinates = True)
        simulate_thetaI_thetab(roiparams, pinholeparams, sigmas, theta, thetab_list, scenario='fixed_budget', illumination_doughnut = True, filenameextension = 'thetaI_thetab_triangulation_doughnut_fb', N_m = N_m, start = 2, stop = 4, num = 20)

    elif plotnumber=='S18':
        #Supplementary Figure 18: 1 x-offset pattern, fixed illumination power and time
        r_p_list = np.array([1, 2, 3, 4, 100000])*r_p
        pinholeparams=np.repeat([[x_p, y_p, r_p]], 1, axis=0)
        simulate_xp_rp(roiparams, pinholeparams, sigmas, theta, r_p_list, scenario='fixed_energy', filenameextension = 'xp_rp_1pattern', N_m = N_m, start = 3*dx_i, stop = 7*dx_i, num = nsteps)
        
        thetab_list = np.array([0, 1, 4, 8, 16])/(dx_i**2)
        pinholeparams=np.repeat([[x_p, y_p, 3*r_p]], 1, axis=0)
        simulate_thetaI_thetab(roiparams, pinholeparams, sigmas, theta, thetab_list, scenario='fixed_energy', filenameextension = 'thetaI_thetab_1pattern', N_m = N_m, start = 2, stop = 4, num = 20)

    elif plotnumber=='S19':
        #Supplementary Figure 19: 1 y-offset pattern, fixed illumination power and time
        r_p_list = np.array([1, 2, 3, 4, 100000])*r_p
        pinholeparams=np.repeat([[x_p, y_p, r_p]], 1, axis=0)
        simulate_yp_rp(roiparams, pinholeparams, sigmas, theta, r_p_list, scenario='fixed_energy', filenameextension = 'yp_rp_1pattern', N_m = N_m, start = 3*dx_i, stop = 7*dx_i, num = nsteps)

    elif plotnumber=='S20':
        #Supplementary Figure 20: 2 x-separated patterns, fixed illumination power and time
        s_list = np.array([0, 1, 2, 3, 4, 5])*r_p
        pinholeparams=np.repeat([[x_p, y_p, 3*r_p]], 2, axis=0)
        simulate_xp_s_2pattern(roiparams, pinholeparams, sigmas, theta, s_list, scenario='fixed_energy', filenameextension = 'xp_s_2pattern', N_m = N_m, start = 3*dx_i, stop = 7*dx_i, num = nsteps)

        thetab_list = np.array([0, 1, 4, 8, 16])/(dx_i**2)
        pinholeparams=np.array([[x_p-4*r_p/2, y_p, 3*r_p], [x_p+4*r_p/2, y_p, 3*r_p]])
        simulate_thetaI_thetab(roiparams, pinholeparams, sigmas, theta, thetab_list, scenario='fixed_energy', filenameextension = 'thetaI_thetab_2pattern', N_m = N_m, start = 2, stop = 4, num = 20)

    elif plotnumber=='S21':
        #Supplementary Figure 21: 2 x-separated patterns with y-offset, fixed illumination power and time
        s_list = np.array([0, 1, 2, 3, 4, 5])*r_p
        pinholeparams=np.repeat([[x_p, 0, 3*r_p]], 2, axis=0)
        simulate_yp_xs_2pattern(roiparams, pinholeparams, sigmas, theta, s_list, scenario='fixed_energy', filenameextension = 'yp_xs_2pattern', N_m = N_m, start = 3*dx_i, stop = 7*dx_i, num = nsteps)

    elif plotnumber=='S22':
        #Supplementary Figure 22: 2 x-separated patterns without pinholes, fixed illumination power and time
        s_list = np.array([0, 1, 2, 3, 4, 5])*r_p
        pinholeparams=np.repeat([[x_p, y_p, 10000*r_p]], 2, axis=0)
        simulate_xp_s_2pattern(roiparams, pinholeparams, sigmas, theta, s_list, scenario='fixed_energy', filenameextension = 'xp_s_2pattern_nopinhole', N_m = N_m, start = 3*dx_i, stop = 7*dx_i, num = nsteps)

        thetab_list = np.array([0, 1, 4, 8, 16])/(dx_i**2)
        pinholeparams=np.array([[x_p-4*r_p/2, y_p, 1000*r_p], [x_p+4*r_p/2, y_p, 1000*r_p]])
        simulate_thetaI_thetab(roiparams, pinholeparams, sigmas, theta, thetab_list, scenario='fixed_energy', filenameextension = 'thetaI_thetab_2pattern_nopinhole', N_m = N_m, start = 2, stop = 4, num = 20)

    elif plotnumber=='S23':
        #Supplementary Figure 23: 2 y-separated patterns, fixed illumination power and time
        s_list = np.array([0, 1, 2, 3, 4, 5])*r_p
        pinholeparams=np.repeat([[x_p, y_p, 3*r_p]], 2, axis=0)
        simulate_yp_s_2pattern(roiparams, pinholeparams, sigmas, theta, s_list, scenario='fixed_energy', filenameextension = 'yp_s_2pattern', N_m = N_m, start = 3*dx_i, stop = 7*dx_i, num = nsteps)

        thetab_list = np.array([0, 1, 4, 8, 16])/(dx_i**2)
        pinholeparams=np.array([[x_p, y_p-4*r_p/2, 3*r_p], [x_p, y_p+4*r_p/2, 3*r_p]])
        simulate_thetaI_thetab(roiparams, pinholeparams, sigmas, theta, thetab_list, scenario='fixed_energy', filenameextension = 'thetaI_thetab_2pattern_rot', N_m = N_m, start = 2, stop = 4, num = 20)

    elif plotnumber=='S24':
        #Supplementary Figure 24: 3 patterns in equilateral triangle configuration without center pinhole, fixed illumination power and time
        r_list = np.array([1, 2, 3, 4, 5])*r_p
        pinholeparams=np.repeat([[0, 0, 3*r_p]], 3, axis=0)
        simulate_xp_r_triangulation(roiparams, pinholeparams, sigmas, theta, r_list, y_p, scenario='fixed_energy', rotation = np.pi/2, pinhole_in_center = False, filenameextension = 'xp_r_triangulation_rot_nocenter_ext', N_m = N_m, start = 3*dx_i, stop = 7*dx_i, num = nsteps)
        
        thetab_list = np.array([0, 1, 4, 8, 16])/(dx_i**2)
        pinholeparams=np.repeat([[0, 0, 3*r_p]], 3, axis=0)
        pinholeparams[:,[0,1]]= coordinates_triangulation([x_p, y_p], 2*r_p, np.pi/2, return_centercoordinates = False)
        simulate_thetaI_thetab(roiparams, pinholeparams, sigmas, theta, thetab_list, scenario='fixed_energy', filenameextension = 'thetaI_thetab_triangulation_rot_nocenter', N_m = N_m, start = 2, stop = 5, num = 20)

    elif plotnumber=='S25':
        #Supplementary Figure 25: 3 patterns in 90 degree rotated equilateral triangle configuration without center pinhole, fixed illumination power and time
        r_list = np.array([1, 2, 3, 4, 5])*r_p
        pinholeparams=np.repeat([[0, 0, 3*r_p]], 3, axis=0)
        simulate_xp_r_triangulation(roiparams, pinholeparams, sigmas, theta, r_list, y_p, scenario='fixed_energy', rotation = 0, pinhole_in_center = False, filenameextension = 'xp_r_triangulation_nocenter_ext', N_m = N_m, start = 3*dx_i, stop = 7*dx_i, num = nsteps)
        
        thetab_list = np.array([0, 1, 4, 8, 16])/(dx_i**2)
        pinholeparams=np.repeat([[0, 0, 3*r_p]], 3, axis=0)
        pinholeparams[:,[0,1]]= coordinates_triangulation([x_p, y_p], 2*r_p, 0, return_centercoordinates = False)
        simulate_thetaI_thetab(roiparams, pinholeparams, sigmas, theta, thetab_list, scenario='fixed_energy', filenameextension = 'thetaI_thetab_triangulation_nocenter', N_m = N_m, start = 2, stop = 4, num = 20)

    elif plotnumber=='S26':
        #Supplementary Figure 26: 4 patterns in equilateral triangle configuration with center pinhole, fixed illumination power and time
        r_list = np.array([0.5, 1, 1.5, 2])*r_p
        pinholeparams=np.repeat([[0, 0, 3*r_p]], 4, axis=0)
        simulate_xp_r_triangulation(roiparams, pinholeparams, sigmas, theta, r_list, y_p, scenario='fixed_energy', rotation = np.pi/2, pinhole_in_center = True, filenameextension = 'xp_r_triangulation_rot', N_m = N_m, start = 3*dx_i, stop = 7*dx_i, num = nsteps)

        thetab_list = np.array([0, 1, 4, 8, 16])/(dx_i**2)
        pinholeparams=np.repeat([[0, 0, 3*r_p]], 4, axis=0)
        pinholeparams[:,[0,1]]= coordinates_triangulation([x_p, y_p], 2*r_p, np.pi/2, return_centercoordinates = True)
        simulate_thetaI_thetab(roiparams, pinholeparams, sigmas, theta, thetab_list, scenario='fixed_energy', filenameextension = 'thetaI_thetab_triangulation_rot', N_m = N_m, start = 2, stop = 4, num = 20)

    elif plotnumber=='S27':
        #Supplementary Figure 27: 4 patterns in 90 degree rotated equilateral triangle configuration with center pinhole, fixed illumination power and time
        r_list = np.array([0.5, 1, 1.5, 2])*r_p
        pinholeparams=np.repeat([[0, 0, 3*r_p]], 4, axis=0)
        simulate_xp_r_triangulation(roiparams, pinholeparams, sigmas, theta, r_list, y_p, scenario='fixed_energy', rotation = 0, pinhole_in_center = True, filenameextension = 'xp_r_triangulation', N_m = N_m, start = 3*dx_i, stop = 7*dx_i, num = nsteps)

        thetab_list = np.array([0, 1, 4, 8, 16])/(dx_i**2)
        pinholeparams=np.repeat([[0, 0, 3*r_p]], 4, axis=0)
        pinholeparams[:,[0,1]]= coordinates_triangulation([x_p, y_p], 2*r_p, 0, return_centercoordinates = True)
        simulate_thetaI_thetab(roiparams, pinholeparams, sigmas, theta, thetab_list, scenario='fixed_energy', filenameextension = 'thetaI_thetab_triangulation', N_m = N_m, start = 2, stop = 4, num = 20)

    elif plotnumber=='S28':
        #Supplementary Figure 28: 3 doughnut-patterns in equilateral triangle configuration without center pinhole, fixed illumination power and time
        r_list = np.array([1, 2, 3, 4, 5])*r_p
        pinholeparams=np.repeat([[0, 0, 3*r_p]], 3, axis=0)
        simulate_xp_r_triangulation(roiparams, pinholeparams, sigmas, theta, r_list, y_p, scenario='fixed_energy', rotation = np.pi/2, pinhole_in_center = False, illumination_doughnut = True, filenameextension = 'xp_r_triangulation_rot_nocenter_doughnut', N_m = N_m, start = 3*dx_i, stop = 7*dx_i, num = nsteps)
        
        thetab_list = np.array([0, 1, 4, 8, 16])/(dx_i**2)
        pinholeparams=np.repeat([[0, 0, 3*r_p]], 3, axis=0)
        pinholeparams[:,[0,1]]= coordinates_triangulation([x_p, y_p], 2*r_p, np.pi/2, return_centercoordinates = False)
        simulate_thetaI_thetab(roiparams, pinholeparams, sigmas, theta, thetab_list, scenario='fixed_energy', illumination_doughnut = True, filenameextension = 'thetaI_thetab_triangulation_rot_nocenter_doughnut', N_m = N_m, start = 2, stop = 5, num = 20)

    elif plotnumber=='S29':
        #Supplementary Figure 29: 3 doughnut-patterns in 90 degree rotated equilateral triangle configuration without center pinhole, fixed illumination power and time
        r_list = np.array([1, 2, 3, 4, 5])*r_p
        pinholeparams=np.repeat([[0, 0, 3*r_p]], 3, axis=0)
        simulate_xp_r_triangulation(roiparams, pinholeparams, sigmas, theta, r_list, y_p, scenario='fixed_energy', rotation = 0, pinhole_in_center = False, illumination_doughnut = True, filenameextension = 'xp_r_triangulation_nocenter_doughnut', N_m = N_m, start = 3*dx_i, stop = 7*dx_i, num = nsteps)        
        
        thetab_list = np.array([0, 1, 4, 8, 16])/(dx_i**2)
        pinholeparams=np.repeat([[0, 0, 3*r_p]], 3, axis=0)
        pinholeparams[:,[0,1]]= coordinates_triangulation([x_p, y_p], 2*r_p, 0, return_centercoordinates = False)
        simulate_thetaI_thetab(roiparams, pinholeparams, sigmas, theta, thetab_list, scenario='fixed_energy', illumination_doughnut = True, filenameextension = 'thetaI_thetab_triangulation_nocenter_doughnut', N_m = N_m, start = 2, stop = 4, num = 20)    

    elif plotnumber=='S30':
        #Supplementary Figure 30: 4 doughnut-patterns in equilateral triangle configuration with center pinhole, fixed illumination power and time
        r_list = np.array([0.5, 1, 2, 3, 4])*r_p
        pinholeparams=np.repeat([[0, 0, 3*r_p]], 4, axis=0)
        simulate_xp_r_triangulation(roiparams, pinholeparams, sigmas, theta, r_list, y_p, scenario='fixed_energy', rotation = np.pi/2, pinhole_in_center = True, illumination_doughnut = True, filenameextension = 'xp_r_triangulation_rot_doughnut_ext', N_m = N_m, start = 3*dx_i, stop = 7*dx_i, num = nsteps)
        
        thetab_list = np.array([0, 1, 4, 8, 16])/(dx_i**2)
        pinholeparams=np.repeat([[0, 0, 3*r_p]], 4, axis=0)
        pinholeparams[:,[0,1]]= coordinates_triangulation([x_p, y_p], 2*r_p, np.pi/2, return_centercoordinates = True)
        simulate_thetaI_thetab(roiparams, pinholeparams, sigmas, theta, thetab_list, scenario='fixed_energy', illumination_doughnut = True, filenameextension = 'thetaI_thetab_triangulation_rot_doughnut', N_m = N_m, start = 2, stop = 5, num = 20)

    elif plotnumber=='S31':
        #Supplementary Figure 31: 4 doughnut-patterns in 90 degree rotated equilateral triangle configuration with center pinhole, fixed illumination power and time
        r_list = np.array([0.5, 1, 2, 3, 4])*r_p
        pinholeparams=np.repeat([[0, 0, 3*r_p]], 4, axis=0)
        simulate_xp_r_triangulation(roiparams, pinholeparams, sigmas, theta, r_list, y_p, scenario='fixed_energy', rotation = 0, pinhole_in_center = True, illumination_doughnut = True, filenameextension = 'xp_r_triangulation_doughnut_ext', N_m = N_m, start = 3*dx_i, stop = 7*dx_i, num = nsteps)

        thetab_list = np.array([0, 1, 4, 8, 16])/(dx_i**2)
        pinholeparams=np.repeat([[0, 0, 3*r_p]], 4, axis=0)
        pinholeparams[:,[0,1]]= coordinates_triangulation([x_p, y_p], 2*r_p, 0, return_centercoordinates = True)
        simulate_thetaI_thetab(roiparams, pinholeparams, sigmas, theta, thetab_list, scenario='fixed_energy', illumination_doughnut = True, filenameextension = 'thetaI_thetab_triangulation_doughnut', N_m = N_m, start = 2, stop = 4, num = 20)

    elif plotnumber=='S32':
        #Supplementary Figure 32: 1 x-offset pattern, fixed signal photon budget with pattern-independent background
        r_p_list = np.array([1, 2, 3, 4, 100000])*r_p
        pinholeparams=np.repeat([[x_p, y_p, r_p]], 1, axis=0)
        simulate_xp_rp(roiparams, pinholeparams, sigmas, theta, r_p_list, scenario='fixed_budget_adaptedbg', filenameextension = 'xp_rp_1pattern_adaptedbg', N_m = N_m, start = 3*dx_i, stop = 7*dx_i, num = nsteps)

        thetab_list = np.array([0, 1, 4, 8, 16])/(dx_i**2)
        pinholeparams=np.repeat([[x_p, y_p, 3*r_p]], 1, axis=0)
        simulate_thetaI_thetab(roiparams, pinholeparams, sigmas, theta, thetab_list, scenario='fixed_budget_adaptedbg', filenameextension = 'thetaI_thetab_1pattern_adaptedbg', N_m = N_m, start = 2, stop = 4, num = 20)

    elif plotnumber=='S33':      
        #Supplementary Figure 33: 1 y-offset pattern, fixed signal photon budget with pattern-independent background
        r_p_list = np.array([1, 2, 3, 4, 100000])*r_p
        pinholeparams=np.repeat([[x_p, y_p, r_p]], 1, axis=0)
        simulate_yp_rp(roiparams, pinholeparams, sigmas, theta, r_p_list, scenario='fixed_budget_adaptedbg', filenameextension = 'yp_rp_1pattern_adaptedbg', N_m = N_m, start = 3*dx_i, stop = 7*dx_i, num = nsteps)

    elif plotnumber=='S34':   
        #Supplementary Figure 34: 2 x-separated patterns, fixed signal photon budget with pattern-independent background
        s_list = np.array([0, 1, 2, 3, 4, 5])*r_p
        pinholeparams=np.repeat([[x_p, y_p, 3*r_p]], 2, axis=0)
        simulate_xp_s_2pattern(roiparams, pinholeparams, sigmas, theta, s_list, scenario='fixed_budget_adaptedbg', filenameextension = 'xp_s_2pattern_adaptedbg', N_m = N_m, start = 3*dx_i, stop = 7*dx_i, num = nsteps)

        thetab_list = np.array([0, 1, 4, 8, 16])/(dx_i**2)
        pinholeparams=np.array([[x_p-4*r_p/2, y_p, 3*r_p], [x_p+4*r_p/2, y_p, 3*r_p]])
        simulate_thetaI_thetab(roiparams, pinholeparams, sigmas, theta, thetab_list, scenario='fixed_budget_adaptedbg', filenameextension = 'thetaI_thetab_2pattern_adaptedbg', N_m = N_m, start = 2, stop = 4, num = 20)

    elif plotnumber=='S35':   
        #Supplementary Figure 35: 2 x-separated patterns with y-offset, fixed signal photon budget with pattern-independent background
        s_list = np.array([0, 1, 2, 3, 4, 5])*r_p
        pinholeparams=np.repeat([[x_p, 0, 3*r_p]], 2, axis=0)
        simulate_yp_xs_2pattern(roiparams, pinholeparams, sigmas, theta, s_list, scenario='fixed_budget_adaptedbg', filenameextension = 'yp_xs_2pattern_adaptedbg', N_m = N_m, start = 3*dx_i, stop = 7*dx_i, num = nsteps)

    elif plotnumber=='S36':   
        #Supplementary Figure 36: 2 x-separated patterns without pinholes, fixed signal photon budget with pattern-independent background
        s_list = np.array([0, 1, 2, 3, 4, 5])*r_p
        pinholeparams=np.repeat([[x_p, y_p, 10000*r_p]], 2, axis=0)
        simulate_xp_s_2pattern(roiparams, pinholeparams, sigmas, theta, s_list, scenario='fixed_budget_adaptedbg', filenameextension = 'xp_s_2pattern_nopinhole_adaptedbg', N_m = N_m, start = 3*dx_i, stop = 7*dx_i, num = nsteps)

        thetab_list = np.array([0, 1, 4, 8, 16])/(dx_i**2)
        pinholeparams=np.array([[x_p-4*r_p/2, y_p, 1000*r_p], [x_p+4*r_p/2, y_p, 1000*r_p]])
        simulate_thetaI_thetab(roiparams, pinholeparams, sigmas, theta, thetab_list, scenario='fixed_budget_adaptedbg', filenameextension = 'thetaI_thetab_2pattern_nopinhole_adaptedbg', N_m = N_m, start = 2, stop = 4, num = 20)

    elif plotnumber=='S37':   
        #Supplementary Figure 37: 2 y-separated patterns, fixed signal photon budget with pattern-independent background
        s_list = np.array([0, 1, 2, 3, 4, 5])*r_p
        pinholeparams=np.repeat([[x_p, y_p, 3*r_p]], 2, axis=0)
        simulate_yp_s_2pattern(roiparams, pinholeparams, sigmas, theta, s_list, scenario='fixed_budget_adaptedbg', filenameextension = 'yp_s_2pattern_adaptedbg', N_m = N_m, start = 3*dx_i, stop = 7*dx_i, num = nsteps)

        thetab_list = np.array([0, 1, 4, 8, 16])/(dx_i**2)
        pinholeparams=np.array([[x_p, y_p-4*r_p/2, 3*r_p], [x_p, y_p+4*r_p/2, 3*r_p]])
        simulate_thetaI_thetab(roiparams, pinholeparams, sigmas, theta, thetab_list, scenario='fixed_budget_adaptedbg', filenameextension = 'thetaI_thetab_2pattern_rot_adaptedbg', N_m = N_m, start = 2, stop = 4, num = 20)

    elif plotnumber=='S38':   
        #Supplementary Figure 38: 3 patterns in equilateral triangle configuration without center pinhole, fixed signal photon budget with pattern-independent background
        r_list = np.array([1, 2, 3, 4, 5])*r_p
        pinholeparams=np.repeat([[0, 0, 3*r_p]], 3, axis=0)
        simulate_xp_r_triangulation(roiparams, pinholeparams, sigmas, theta, r_list, y_p, rotation=np.pi/2, pinhole_in_center = False, scenario='fixed_budget_adaptedbg', filenameextension = 'xp_r_triangulation_rot_nocenter_ext_adaptedbg', N_m = N_m, start = 3*dx_i, stop = 7*dx_i, num = nsteps)

        thetab_list = np.array([0, 1, 4, 8, 16])/(dx_i**2)
        pinholeparams=np.repeat([[0, 0, 3*r_p]], 3, axis=0)
        pinholeparams[:,[0,1]]= coordinates_triangulation([x_p, y_p], 2*r_p, np.pi/2, return_centercoordinates = False)
        simulate_thetaI_thetab(roiparams, pinholeparams, sigmas, theta, thetab_list, scenario='fixed_budget_adaptedbg', filenameextension = 'thetaI_thetab_triangulation_rot_nocenter_adaptedbg', N_m = N_m, start = 2, stop = 5, num = 20)

    elif plotnumber=='S39':  
        #Supplementary Figure 39: 3 patterns in 90 degree rotated equilateral triangle configuration without center pinhole, fixed signal photon budget with pattern-independent background
        r_list = np.array([1, 2, 3, 4, 5])*r_p
        pinholeparams=np.repeat([[0, 0, 3*r_p]], 3, axis=0)
        simulate_xp_r_triangulation(roiparams, pinholeparams, sigmas, theta, r_list, y_p, rotation=0, pinhole_in_center = False, scenario='fixed_budget_adaptedbg', filenameextension = 'xp_r_triangulation_nocenter_ext_adaptedbg', N_m = N_m, start = 3*dx_i, stop = 7*dx_i, num = nsteps)

        thetab_list = np.array([0, 1, 4, 8, 16])/(dx_i**2)
        pinholeparams=np.repeat([[0, 0, 3*r_p]], 3, axis=0)
        pinholeparams[:,[0,1]]= coordinates_triangulation([x_p, y_p], 2*r_p, 0, return_centercoordinates = False)
        simulate_thetaI_thetab(roiparams, pinholeparams, sigmas, theta, thetab_list, scenario='fixed_budget_adaptedbg', filenameextension = 'thetaI_thetab_triangulation_nocenter_adaptedbg', N_m = N_m, start = 2, stop = 5, num = 20)

    elif plotnumber=='S40':   
        #Supplementary Figure 40: 4 patterns in equilateral triangle configuration with center pinhole, fixed signal photon budget with pattern-independent background
        r_list = np.array([0.5, 1, 1.5, 2])*r_p
        pinholeparams=np.repeat([[0, 0, 3*r_p]], 4, axis=0)
        simulate_xp_r_triangulation(roiparams, pinholeparams, sigmas, theta, r_list, y_p, rotation=np.pi/2, pinhole_in_center = True, scenario='fixed_budget_adaptedbg', filenameextension = 'xp_r_triangulation_rot_adaptedbg', N_m = N_m, start = 3*dx_i, stop = 7*dx_i, num = nsteps)

        thetab_list = np.array([0, 1, 4, 8, 16])/(dx_i**2)
        pinholeparams=np.repeat([[0, 0, 3*r_p]], 4, axis=0)
        pinholeparams[:,[0,1]]= coordinates_triangulation([x_p, y_p], 2*r_p, np.pi/2, return_centercoordinates = True)
        simulate_thetaI_thetab(roiparams, pinholeparams, sigmas, theta, thetab_list, scenario='fixed_budget_adaptedbg', filenameextension = 'thetaI_thetab_triangulation_rot_adaptedbg', N_m = N_m, start = 2, stop = 5, num = 20)

    elif plotnumber=='S41':  
        #Supplementary Figure 41: 4 patterns in 90 degree rotated equilateral triangle configuration with center pinhole, fixed signal photon budget with pattern-independent background
        r_list = np.array([0.5, 1, 1.5, 2])*r_p
        pinholeparams=np.repeat([[0, 0, 3*r_p]], 4, axis=0)
        simulate_xp_r_triangulation(roiparams, pinholeparams, sigmas, theta, r_list, y_p, rotation=0, pinhole_in_center = True, scenario='fixed_budget_adaptedbg', filenameextension = 'xp_r_triangulation_adaptedbg', N_m = N_m, start = 3*dx_i, stop = 7*dx_i, num = nsteps)
        
        thetab_list = np.array([0, 1, 4, 8, 16])/(dx_i**2)
        pinholeparams=np.repeat([[0, 0, 3*r_p]], 4, axis=0)
        pinholeparams[:,[0,1]]= coordinates_triangulation([x_p, y_p], 2*r_p, 0, return_centercoordinates = True)
        simulate_thetaI_thetab(roiparams, pinholeparams, sigmas, theta, thetab_list, scenario='fixed_budget_adaptedbg', filenameextension = 'thetaI_thetab_triangulation_adaptedbg', N_m = N_m, start = 2, stop = 5, num = 20)
        
    elif plotnumber=='S42':
        #Supplementary Figure 42: 3 doughnut-patterns in equilateral triangle configuration without center pinhole, fixed signal photon budget with pattern-independent background
        r_list = np.array([1, 2, 3, 4, 5])*r_p
        pinholeparams=np.repeat([[0, 0, 3*r_p]], 3, axis=0)
        simulate_xp_r_triangulation(roiparams, pinholeparams, sigmas, theta, r_list, y_p, scenario='fixed_budget_adaptedbg', rotation = np.pi/2, pinhole_in_center = False, illumination_doughnut = True, filenameextension = 'xp_r_triangulation_rot_nocenter_doughnut_adaptedbg', N_m = N_m, start = 3*dx_i, stop = 7*dx_i, num = nsteps)
        
        thetab_list = np.array([0, 1, 4, 8, 16])/(dx_i**2)
        pinholeparams=np.repeat([[0, 0, 3*r_p]], 3, axis=0)
        pinholeparams[:,[0,1]]= coordinates_triangulation([x_p, y_p], 2*r_p, np.pi/2, return_centercoordinates = False)
        simulate_thetaI_thetab(roiparams, pinholeparams, sigmas, theta, thetab_list, scenario='fixed_budget_adaptedbg', illumination_doughnut = True, filenameextension = 'thetaI_thetab_triangulation_rot_nocenter_doughnut_adaptedbg', N_m = N_m, start = 2, stop = 5, num = 20)

    elif plotnumber=='S43':
        #Supplementary Figure 43: 3 doughnut-patterns in 90 degree rotated equilateral triangle configuration without center pinhole, fixed signal photon budget with pattern-independent background
        r_list = np.array([1, 2, 3, 4, 5])*r_p
        pinholeparams=np.repeat([[0, 0, 3*r_p]], 3, axis=0)
        simulate_xp_r_triangulation(roiparams, pinholeparams, sigmas, theta, r_list, y_p, scenario='fixed_budget_adaptedbg', rotation = 0, pinhole_in_center = False, illumination_doughnut = True, filenameextension = 'xp_r_triangulation_nocenter_doughnut_adaptedbg', N_m = N_m, start = 3*dx_i, stop = 7*dx_i, num = nsteps)        
        
        thetab_list = np.array([0, 1, 4, 8, 16])/(dx_i**2)
        pinholeparams=np.repeat([[0, 0, 3*r_p]], 3, axis=0)
        pinholeparams[:,[0,1]]= coordinates_triangulation([x_p, y_p], 2*r_p, 0, return_centercoordinates = False)
        simulate_thetaI_thetab(roiparams, pinholeparams, sigmas, theta, thetab_list, scenario='fixed_budget_adaptedbg', illumination_doughnut = True, filenameextension = 'thetaI_thetab_triangulation_nocenter_doughnut_adaptedbg', N_m = N_m, start = 2, stop = 4, num = 20)    

    elif plotnumber=='S44':
        #Supplementary Figure 44: 4 doughnut-patterns in equilateral triangle configuration with center pinhole, fixed signal photon budget with pattern-independent background
        r_list = np.array([0.5, 1, 2, 3, 4])*r_p
        pinholeparams=np.repeat([[0, 0, 3*r_p]], 4, axis=0)
        simulate_xp_r_triangulation(roiparams, pinholeparams, sigmas, theta, r_list, y_p, scenario='fixed_budget_adaptedbg', rotation = np.pi/2, pinhole_in_center = True, illumination_doughnut = True, filenameextension = 'xp_r_triangulation_rot_doughnut_ext_adaptedbg', N_m = N_m, start = 3*dx_i, stop = 7*dx_i, num = nsteps)
        
        thetab_list = np.array([0, 1, 4, 8, 16])/(dx_i**2)
        pinholeparams=np.repeat([[0, 0, 3*r_p]], 4, axis=0)
        pinholeparams[:,[0,1]]= coordinates_triangulation([x_p, y_p], 2*r_p, np.pi/2, return_centercoordinates = True)
        simulate_thetaI_thetab(roiparams, pinholeparams, sigmas, theta, thetab_list, scenario='fixed_budget_adaptedbg', illumination_doughnut = True, filenameextension = 'thetaI_thetab_triangulation_rot_doughnut_adaptedbg', N_m = N_m, start = 2, stop = 5, num = 20)

    elif plotnumber=='S45':
        #Supplementary Figure 45: 4 doughnut-patterns in 90 degree rotated equilateral triangle configuration with center pinhole, fixed signal photon budget with pattern-independent background
        r_list = np.array([0.5, 1, 2, 3, 4])*r_p
        pinholeparams=np.repeat([[0, 0, 3*r_p]], 4, axis=0)
        simulate_xp_r_triangulation(roiparams, pinholeparams, sigmas, theta, r_list, y_p, scenario='fixed_budget_adaptedbg', rotation = 0, pinhole_in_center = True, illumination_doughnut = True, filenameextension = 'xp_r_triangulation_doughnut_ext_adaptedbg', N_m = N_m, start = 3*dx_i, stop = 7*dx_i, num = nsteps)

        thetab_list = np.array([0, 1, 4, 8, 16])/(dx_i**2)
        pinholeparams=np.repeat([[0, 0, 3*r_p]], 4, axis=0)
        pinholeparams[:,[0,1]]= coordinates_triangulation([x_p, y_p], 2*r_p, 0, return_centercoordinates = True)
        simulate_thetaI_thetab(roiparams, pinholeparams, sigmas, theta, thetab_list, scenario='fixed_budget_adaptedbg', illumination_doughnut = True, filenameextension = 'thetaI_thetab_triangulation_doughnut_adaptedbg', N_m = N_m, start = 2, stop = 4, num = 20)