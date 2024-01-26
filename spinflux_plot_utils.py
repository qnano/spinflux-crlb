# -*- coding: utf-8 -*-
"""
Utilities and plotting code underlying the SpinFlux figures.
    - load_thetaI_Nm(fileextension) loads simulation data for the simulation file described in fileextension.
    - load_thetaI_Nm(fileextension) loads simulation conditions for any configuration, as a function of the signal photon budget thetaI for different background counts thetab.
    - load_xp_rp(fileextension) loads simulation conditions for the single-pattern configuration, as a function of the pattern focus x-position for different pinhole radii rp.
    - load_yp_rp(fileextension) loads simulation conditions for the single-pattern configuration, as a function of the pattern focus y-position for different pinhole radii rp.
    - load_xp_s(fileextension) loads simulation conditions for the two-pattern configuration, as a function of the pattern focus x-position for different separations s.
    - load_xp_r(fileextension) loads simulation conditions for the triangle configuration, as a function of the pattern focus x-position for different spacings r.
    - load_thetaI_Nm(fileextension) loads simulation conditions about the effect of discretization on the CRLB, as a function of the amount of mesh pixels Nm for different values of the signal photon count thetaI.
    - plot(plotnumber, panel) plots the SpinFlux Cramér-Rao lower bound (CRLB) figures.
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import matplotlib as mpl
import matplotlib.colors as mcolors

new_rc_params = {
"font.size": 13,
"axes.labelsize": 13}
mpl.rcParams.update(new_rc_params)

def load_4(filenameextension, loadpath = './SimData/'):
    '''
    load_thetaI_Nm(fileextension) loads simulation data for the simulation file
    described in fileextension.
    
    Input:
        filenameextension: string, describing simulation file to be loaded in.
        
    Outputs:
        CRLBx: Cramér-Rao lower bound (CRLB)
        musum: Illumination-intensity normalized total amount of photons in region of interest (ROI)
        signalsum: Illumination-intensity normalized total amount of signal photons in ROI
        EfBgsum: Illumination-intensity normalized total amount of background photons in ROI
    '''
    saveCRLBx = 'CRLBx_' + filenameextension + '.npy'
    savemusum = 'musum_' + filenameextension + '.npy'
    savesignalsum = 'signalsum_' + filenameextension + '.npy'
    saveEfBgsum = 'EfBgsum_' + filenameextension + '.npy'
    
    CRLBx = np.load(loadpath+saveCRLBx)
    musum = np.load(loadpath+savemusum)
    signalsum = np.load(loadpath+savesignalsum)
    EfBgsum = np.load(loadpath+saveEfBgsum)
    return CRLBx, musum, signalsum, EfBgsum

def load_thetaI_thetab(filenameextension, loadpath = './SimData/'):
    '''
    load_thetaI_Nm(fileextension) loads simulation conditions for any configuration,
    as a function of the signal photon budget thetaI for different background counts thetab.
    
    Input:
        filenameextension: string, describing simulation file to be loaded in.

    Outputs:
        thetab_list: list of background photon count values thetab
        thetaIrange: range of signal photon budget values thetaI
    '''
    savethetablist = 'thetablist_' + filenameextension + '.npy'
    savethetaIrange = 'thetaIrange_' + filenameextension + '.npy'

    thetab_list = np.load(loadpath+savethetablist)
    thetaIrange = np.load(loadpath+savethetaIrange)
    return thetab_list, thetaIrange

def load_xp_rp(filenameextension, loadpath = './SimData/'):
    '''
    load_xp_rp(fileextension) loads simulation conditions for the single-pattern configuration, 
    as a function of the pattern focus x-position for different pinhole radii rp.
    
    Input:
        filenameextension: string, describing simulation file to be loaded in.

    Outputs:
        r_p_list: list of pinhole radii r_p
        xrange: range of pattern focus x-positions
    '''
    saverplist = 'rplist_' + filenameextension + '.npy'
    savexrange = 'xrange_' + filenameextension + '.npy'

    r_p_list = np.load(loadpath+saverplist)
    xrange = np.load(loadpath+savexrange)
    return r_p_list, xrange

def load_yp_rp(filenameextension, loadpath = './SimData/'):
    '''
    load_yp_rp(fileextension) loads simulation conditions for the single-pattern configuration, 
    as a function of the pattern focus y-position for different pinhole radii rp.
    
    Input:
        filenameextension: string, describing simulation file to be loaded in.

    Outputs:
        r_p_list: list of pinhole radii r_p
        yrange: range of pattern focus y-positions
    '''
    saverplist = 'rplist_' + filenameextension + '.npy'
    saveyrange = 'yrange_' + filenameextension + '.npy'

    r_p_list = np.load(loadpath+saverplist)
    yrange = np.load(loadpath+saveyrange)
    return r_p_list, yrange

def load_xp_s(filenameextension, loadpath = './SimData/'):
    '''
    load_xp_s(fileextension) loads simulation conditions for the two-pattern configuration, 
    as a function of the pattern focus x-position for different separations s.
    
    Input:
        filenameextension: string, describing simulation file to be loaded in.

    Outputs:
        s_list: list of pinhole separations s
        xrange: range of pattern focus x-positions
    '''
    saveslist = 'slist_' + filenameextension + '.npy'
    savexrange = 'xrange_' + filenameextension + '.npy'

    s_list = np.load(loadpath+saveslist)
    xrange = np.load(loadpath+savexrange)
    return s_list, xrange

def load_xp_r(filenameextension, loadpath = './SimData/'):
    '''
    load_xp_r(fileextension) loads simulation conditions for the triangle configuration, 
    as a function of the pattern focus x-position for different spacings r.
    
    Input:
        filenameextension: string, describing simulation file to be loaded in.
        
    Outputs:
        r_list: list of pinhole spacings r
        xrange: range of pattern focus x-positions
    '''
    saverlist = 'rlist_' + filenameextension + '.npy'
    savexrange = 'xrange_' + filenameextension + '.npy'

    r_list = np.load(loadpath+saverlist)
    xrange = np.load(loadpath+savexrange)
    return r_list, xrange

def load_thetaI_Nm(filenameextension, loadpath = './SimData/'):
    '''
    load_thetaI_Nm(fileextension) loads simulation conditions about the effect of 
    discretization on the CRLB, as a function of the amount of mesh pixels Nm
    for different values of the signal photon count thetaI.
    
    Input:
        filenameextension: string, describing simulation file to be loaded in.
        
    Outputs:
        thetaI_list: list of signal photon budget values thetaI
        Nmrange: range of mesh pixel values Nm
    '''
    savethetaIlist = 'thetaIlist_' + filenameextension + '.npy'
    saveNmrange = 'Nmrange_' + filenameextension + '.npy'

    thetaI_list = np.load(loadpath+savethetaIlist)
    Nmrange = np.load(loadpath+saveNmrange)
    return thetaI_list, Nmrange

def plot(plotnumber, panel, loadpath = './SimData/'):
    '''
    plot(plotnumber, panel) plots the SpinFlux Cramér-Rao lower bound (CRLB) figures.
    
    Input:
        plotnumber: string, describing the figure of which CRLB is to be simulated.
            Options:
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
        panel: int, describing the panel of the figure to be created.
            Options:
                0: Legend
                1: CRLB as a function of spatial difference between pattern focus and emitter position
                2: Improvement of CRLB over SMLM as a function of spatial difference between pattern focus and emitter position
                3: Illumination-intensity normalized total amount of photons as a function of spatial difference between pattern focus and emitter position
                4: Illumination-intensity normalized total amount of signal photons as a function of spatial difference between pattern focus and emitter position
                5: Illumination-intensity normalized average amount of background photons per pixel as a function of spatial difference between pattern focus and emitter position
                6: CRLB as a function of signal photon budget for different background counts
    
    Optional arguments:
        loadpath: Path to load simulation data from. Defaults to ./SimData/.
    '''
    #Parameters
    L_illum = 546*10**-9        #Illumination wavelength (m) (cf. Sirinakis et al. (2022))
    L_emission = 600*10**(-9)   #Emission wavelength (m) (cf. Sirinakis et al. (2022))
    NA = 1.35                   #Numerical aperture (cf. Sirinakis et al. (2022))
    
    dx_i = 65*10**-9    #Camera pixel size (m)
    Npixels = 10        #Amount of camera pixels in each direction of the region of interest
    
    x_p = Npixels*dx_i/2     #Pinhole x-coordinate (m)
    y_p = Npixels*dx_i/2     #Pinhole y-coordinate (m)
    r_p = 0.21*L_emission/NA        #Radius pinhole (m)
    
    sigmas = [0.21*L_illum/NA, 0.21*L_emission/NA]    #Illumination and emission PSF standard deviation (m) (cf. Huang and Olivo-Marin (2014))
    theta = [Npixels*dx_i/2, Npixels*dx_i/2, 2000, 8/(dx_i**2)] #Emitter x- and y-position (m), expected signal photon budget (photons) and expected background (photons/m**2)
    
    CRLBx_ISM = 1.7736218191848463e-09
    imp_ISM = 1.4765639096640992
    
    CRLBx_FR = 1.246677240431214e-09
    imp_FR = 2.100676809256009

    if plotnumber == '3-S4':
        CRLBx, musum, signalsum, Bisum = load_4('xp_rp_1pattern_fb')
        CRLBx_SMLM, musum_SMLM, signalsum_SMLM, Bisum_SMLM = load_4('SMLM')
        r_p_list, xrange = load_xp_rp('xp_rp_1pattern_fb')
        
        fig, ax = plt.subplots(1,1)
        fig.set_size_inches(5, 5, forward=True)
        
        if panel == 1:
            ax.semilogy((xrange-theta[0])/(10**-9),CRLBx.T[:,0:4]*10**9)
            ax.semilogy((xrange-theta[0])/(10**-9),CRLBx.T[:,-1]*10**9)
            ax.semilogy((xrange-theta[0])/(10**-9),CRLBx_SMLM*10**9*np.ones(np.size(xrange)), color='k')
            ax.semilogy((xrange-theta[0])/(10**-9),CRLBx_ISM*10**9*np.ones(np.size(xrange)), color='k', linestyle='--')
            ax.semilogy((xrange-theta[0])/(10**-9),CRLBx_FR*10**9*np.ones(np.size(xrange)), color='k', linestyle=':')
            ax.set_xlabel('Emitter-pinhole distance ' + r'$(x_p - \theta_x)$' + ' (nm)')
            ax.set_ylabel(r'$x$'+'-Cramér-Rao lower bound (nm)')
            ax.grid()
            ax.set_xlim(left=-130, right = 130)
            ax.set_ylim(bottom = 0.6, top = 60)
            ax.tick_params(axis='both', which='major')
            ax.tick_params(axis='both', which='minor')
            
        if panel == 2:
            ax.plot((xrange-theta[0])/(10**-9),CRLBx_SMLM/CRLBx.T[:,0:4])
            ax.plot((xrange-theta[0])/(10**-9),CRLBx_SMLM/CRLBx.T[:,-1])
            ax.plot((xrange-theta[0])/(10**-9),imp_ISM*np.ones(np.size(xrange)), color='k', linestyle='--')
            ax.plot((xrange-theta[0])/(10**-9),imp_FR*np.ones(np.size(xrange)), color='k', linestyle=':')
            ax.set_xlabel('Emitter-pinhole distance ' + r'$(x_p - \theta_x)$' + ' (nm)')
            ax.set_ylabel('Improvement over SMLM in ' + r'$x$'+ '-direction')
            ax.grid()
            ax.set_xlim(left=-130, right = 130)
            ax.set_ylim(bottom = 0, top = 3.5)
            ax.tick_params(axis='both', which='major')
            ax.tick_params(axis='both', which='minor')                     

        if panel == 3:
            ax.plot((xrange-theta[0])/(10**-9),musum.T[:,0:4])
            ax.plot((xrange-theta[0])/(10**-9),musum.T[:,-1])
            ax.plot((xrange-theta[0])/(10**-9),musum_SMLM*np.ones(np.size(xrange)), color='k')
            ax.set_xlabel('Emitter-pinhole distance ' + r'$(x_p - \theta_x)$' + ' (nm)')
            ax.set_ylabel('Expected amount of photons')
            ax.grid()
            ax.set_xlim(left=-130, right = 130)
            ax.set_ylim(bottom = 0, top = 2900)
            ax.tick_params(axis='both', which='major')
            ax.tick_params(axis='both', which='minor')                  

        if panel == 4:
            ax.plot((xrange-theta[0])/(10**-9),signalsum.T[:,0:4])
            ax.plot((xrange-theta[0])/(10**-9),signalsum.T[:,-1])
            ax.plot((xrange-theta[0])/(10**-9),signalsum_SMLM*np.ones(np.size(xrange)), color='k')
            ax.set_xlabel('Emitter-pinhole distance ' + r'$(x_p - \theta_x)$' + ' (nm)')
            ax.set_ylabel('Illumination-intensity normalized' + '\n' + 'amount of signal photons')
            ax.grid()
            ax.set_xlim(left=-130, right = 130)
            ax.set_ylim(bottom = 0, top = 2100)
            ax.tick_params(axis='both', which='major')
            ax.tick_params(axis='both', which='minor')   

        if panel == 5:
            ax.plot((xrange-theta[0])/(10**-9),Bisum.T[:,0:4]/100)
            ax.plot((xrange-theta[0])/(10**-9),Bisum.T[:,-1]/100)
            ax.plot((xrange-theta[0])/(10**-9),Bisum_SMLM/100*np.ones(np.size(xrange)), color='k')
            ax.set_xlabel('Emitter-pinhole distance ' + r'$(x_p - \theta_x)$' + ' (nm)')
            ax.set_ylabel('Illumination-intensity normalized average' + '\n' + 'amount of background photons per pixel')
            ax.grid()
            ax.set_xlim(left=-130, right = 130)
            ax.set_ylim(bottom = 0, top = 9)
            ax.tick_params(axis='both', which='major')
            ax.tick_params(axis='both', which='minor')   

        if panel == 6:
            CRLBx, musum, signalsum, Bisum = load_4('thetaI_thetab_1pattern_fb')
            thetab_list, thetaIrange = load_thetaI_thetab('thetaI_thetab_1pattern_fb')
            
            ax.loglog(thetaIrange,CRLBx.T*10**9)
            ax.set_xlabel('Expected signal photon budget ' + r'$\theta_I$' + ' (photons)')
            ax.set_ylabel(r'$x$'+'-Cramér-Rao lower bound (nm)')
            ax.grid()
            ax.set_xlim(left=100, right = 10000)
            ax.set_ylim(bottom = 0.1, top = 40)
            ax.tick_params(axis='both', which='major')
            ax.tick_params(axis='both', which='minor')  
            
            #Legend
            C0_line = mlines.Line2D([0], [0],color='C0', label=(r'$\theta_b = 0$'))
            C1_line = mlines.Line2D([0], [0],color='C1', label=(r'$\theta_b = 1$'))
            C2_line = mlines.Line2D([0], [0],color='C2', label=(r'$\theta_b = 4$'))
            C3_line = mlines.Line2D([0], [0],color='C3', label=(r'$\theta_b = 8$'))
            C8_line = mlines.Line2D([0], [0],color='C4', label=(r'$\theta_b = 16$'))
            C9_line = mlines.Line2D([0], [0],color='None', label=('(in phot/px.)'))
    
            ax.legend(handles=[C0_line, C1_line, C2_line, C3_line, C8_line, C9_line], ncol = 2, prop={'size':13}, labelspacing = 1)    
            
        if panel == 0:
            #Legend
            C0_line = mlines.Line2D([0], [0],color='C0', label=r'$r_p = 1\sigma_{\mathrm{PSF}}$')
            C1_line = mlines.Line2D([0], [0],color='C1', label=r'$r_p = 2\sigma_{\mathrm{PSF}}$')
            C2_line = mlines.Line2D([0], [0],color='C2', label=r'$r_p = 3\sigma_{\mathrm{PSF}}$')
            C3_line = mlines.Line2D([0], [0],color='C3', label=r'$r_p = 4\sigma_{\mathrm{PSF}}$')
            C8_line = mlines.Line2D([0], [0],color='C4', label=r'$r_p = \infty$')
            C9_line = mlines.Line2D([0], [0],color='k', label=r'SMLM (Uniform' + '\n' + 'illumination intensity)')
            C10_line = mlines.Line2D([0], [0],color='k', linestyle='--', label=r'Localization on ISM' + '\n' + 'reconstruction data')
            C11_line = mlines.Line2D([0], [0],color='k', linestyle=':', label=r'Localization on' + '\n' + 'Fourier reweighted ISM' + '\n' + 'reconstruction data')
    
            ax.legend(handles=[C0_line, C1_line, C2_line, C3_line, C8_line, C9_line, C10_line, C11_line], prop={'size':13}, loc='upper left', bbox_to_anchor=(1.04,1), labelspacing = 1)
        
        plt.savefig(('./Figures/3-S4/'+str(int(panel))) + '.svg')

    if plotnumber == '3-heatmap':
        fig, ax = plt.subplots(1,1)
        fig.set_size_inches(5, 5, forward=True)
        
        xnum=51
        ynum=51
        CRLBx = np.zeros((xnum, ynum))
        for yindex in range(ynum):
            CRLBx[yindex,:], _, _, _ = load_4(('xp_rp_1pattern_heatmap_' + str(int(yindex)) + '_fb'))
        CRLBx_SMLM, _, _, _ = load_4('SMLM')

        if panel == 1:
            divnorm = mcolors.LogNorm(vmin=2, vmax=np.max(CRLBx*10**9))
            img = ax.imshow(CRLBx*10**9, cmap='viridis', interpolation='catrom', aspect='equal', extent = (-2*dx_i*10**9, 2*dx_i*10**9, -2*dx_i*10**9, 2*dx_i*10**9), norm = divnorm)
            ax.set_xlabel('Emitter-pinhole distance ' + r'$(x_p - \theta_x)$' + ' (nm)')
            ax.set_ylabel('Emitter-pinhole distance ' + r'$(y_p - \theta_y)$' + ' (nm)')
            ax.tick_params(axis='both', which='major')
            ax.tick_params(axis='both', which='minor') 
            
            plt.savefig(('./Figures/3-heatmap/heatmap_'+str(int(panel))) + '.svg')
            
            fig, ax = plt.subplots(1,1)
            fig.set_size_inches(5, 5, forward=True)
            cbar = plt.colorbar(img)
            cbar.set_label(r'$x$'+'-Cramér-Rao lower bound (nm)')
            
            plt.savefig(('./Figures/3-heatmap/colorbar_heatmap_'+str(int(panel))) + '.svg')
            
        if panel == 2: 
            divnorm = mcolors.TwoSlopeNorm(vcenter=1)
            img = ax.imshow(CRLBx_SMLM/CRLBx, cmap='RdYlGn', interpolation='catrom', aspect='equal', extent = (-2*dx_i*10**9, 2*dx_i*10**9, -2*dx_i*10**9, 2*dx_i*10**9), norm=divnorm)
            ax.set_xlabel('Emitter-pinhole distance ' + r'$(x_p - \theta_x)$' + ' (nm)')
            ax.set_ylabel('Emitter-pinhole distance ' + r'$(y_p - \theta_y)$' + ' (nm)')
            ax.tick_params(axis='both', which='major')
            ax.tick_params(axis='both', which='minor') 
            
            plt.savefig(('./Figures/3-heatmap/heatmap_'+str(int(panel))) + '.svg')
            
            fig, ax = plt.subplots(1,1)
            fig.set_size_inches(5, 5, forward=True)
            cbar = plt.colorbar(img)
            cbar.set_label('Improvement over SMLM in ' + r'$x$'+ '-direction')
            
            plt.savefig(('./Figures/3-heatmap/colorbar_heatmap_'+str(int(panel))) + '.svg')            

    if plotnumber == '4de-S6':
        CRLBx, musum, signalsum, Bisum = load_4('xp_s_2pattern_fb')
        CRLBx_SMLM, musum_SMLM, signalsum_SMLM, Bisum_SMLM = load_4('SMLM')
        s_list, xrange = load_xp_s('xp_s_2pattern_fb')
        
        fig, ax = plt.subplots(1,1)
        fig.set_size_inches(5, 5, forward=True)
        
        if panel == 1:
            ax.semilogy((xrange-theta[0])/(10**-9),CRLBx.T*10**9)
            ax.semilogy((xrange-theta[0])/(10**-9),CRLBx_SMLM*10**9*np.ones(np.size(xrange)), color='k')
            ax.semilogy((xrange-theta[0])/(10**-9),CRLBx_ISM*10**9*np.ones(np.size(xrange)), color='k', linestyle='--')
            ax.semilogy((xrange-theta[0])/(10**-9),CRLBx_FR*10**9*np.ones(np.size(xrange)), color='k', linestyle=':')
            ax.set_xlabel('Emitter-focus distance ' + r'$(x_f - \theta_x)$' + ' (nm)')
            ax.set_ylabel(r'$x$'+'-Cramér-Rao lower bound (nm)')
            ax.grid()
            ax.set_xlim(left=-130, right = 130)
            ax.set_ylim(bottom = 0.6, top = 60)
            ax.tick_params(axis='both', which='major')
            ax.tick_params(axis='both', which='minor')
            
        if panel == 2:
            ax.plot((xrange-theta[0])/(10**-9),CRLBx_SMLM/CRLBx.T)
            ax.plot((xrange-theta[0])/(10**-9),imp_ISM*np.ones(np.size(xrange)), color='k', linestyle='--')
            ax.plot((xrange-theta[0])/(10**-9),imp_FR*np.ones(np.size(xrange)), color='k', linestyle=':')
            ax.set_xlabel('Emitter-focus distance ' + r'$(x_f - \theta_x)$' + ' (nm)')
            ax.set_ylabel('Improvement over SMLM in ' + r'$x$'+ '-direction')
            ax.grid()
            ax.set_xlim(left=-130, right = 130)
            ax.set_ylim(bottom = 0, top = 3.5)
            ax.tick_params(axis='both', which='major')
            ax.tick_params(axis='both', which='minor')                     

        if panel == 3:
            ax.plot((xrange-theta[0])/(10**-9),musum.T)
            ax.plot((xrange-theta[0])/(10**-9),musum_SMLM*np.ones(np.size(xrange)), color='k')
            ax.set_xlabel('Emitter-focus distance ' + r'$(x_f - \theta_x)$' + ' (nm)')
            ax.set_ylabel('Expected amount of photons')
            ax.grid()
            ax.set_xlim(left=-130, right = 130)
            ax.set_ylim(bottom = 0, top = 5000)
            ax.tick_params(axis='both', which='major')
            ax.tick_params(axis='both', which='minor')                  

        if panel == 4:
            ax.plot((xrange-theta[0])/(10**-9),signalsum.T)
            ax.plot((xrange-theta[0])/(10**-9),signalsum_SMLM*np.ones(np.size(xrange)), color='k')
            ax.set_xlabel('Emitter-focus distance ' + r'$(x_f - \theta_x)$' + ' (nm)')
            ax.set_ylabel('Illumination-intensity normalized' + '\n' + 'amount of signal photons')
            ax.grid()
            ax.set_xlim(left=-130, right = 130)
            ax.set_ylim(bottom = 0, top = 2100)
            ax.tick_params(axis='both', which='major')
            ax.tick_params(axis='both', which='minor')   

        if panel == 5:
            ax.plot((xrange-theta[0])/(10**-9),Bisum.T/100)
            ax.plot((xrange-theta[0])/(10**-9),Bisum_SMLM/100*np.ones(np.size(xrange)), color='k')
            ax.set_xlabel('Emitter-focus distance ' + r'$(x_f - \theta_x)$' + ' (nm)')
            ax.set_ylabel('Illumination-intensity normalized average' + '\n' + 'amount of background photons per pixel')
            ax.grid()
            ax.set_xlim(left=-130, right = 130)
            ax.set_ylim(bottom = 0, top = 30)
            ax.tick_params(axis='both', which='major')
            ax.tick_params(axis='both', which='minor')   

        if panel == 6:
            CRLBx, musum, signalsum, Bisum = load_4('thetaI_thetab_2pattern_fb')
            thetab_list, thetaIrange = load_thetaI_thetab('thetaI_thetab_2pattern_fb')
            
            ax.loglog(thetaIrange,CRLBx.T*10**9)
            ax.set_xlabel('Expected signal photon budget ' + r'$\theta_I$' + ' (photons)')
            ax.set_ylabel(r'$x$'+'-Cramér-Rao lower bound (nm)')
            ax.grid()
            ax.set_xlim(left=100, right = 10000)
            ax.set_ylim(bottom = 0.1, top = 40)
            ax.tick_params(axis='both', which='major')
            ax.tick_params(axis='both', which='minor')  
            
            #Legend
            C0_line = mlines.Line2D([0], [0],color='C0', label=(r'$\theta_b = 0$'))
            C1_line = mlines.Line2D([0], [0],color='C1', label=(r'$\theta_b = 1$'))
            C2_line = mlines.Line2D([0], [0],color='C2', label=(r'$\theta_b = 4$'))
            C3_line = mlines.Line2D([0], [0],color='C3', label=(r'$\theta_b = 8$'))
            C8_line = mlines.Line2D([0], [0],color='C4', label=(r'$\theta_b = 16$'))
            C9_line = mlines.Line2D([0], [0],color='None', label=('(in phot/px.)'))
    
            ax.legend(handles=[C0_line, C1_line, C2_line, C3_line, C8_line, C9_line], ncol = 2, prop={'size':13}, labelspacing = 1)     
            
        if panel == 0:
            #Legend
            C0_line = mlines.Line2D([0], [0],color='C0', label=r'$s = 0\sigma_{\mathrm{PSF}}$')
            C1_line = mlines.Line2D([0], [0],color='C1', label=r'$s= 1\sigma_{\mathrm{PSF}}$')
            C2_line = mlines.Line2D([0], [0],color='C2', label=r'$s= 2\sigma_{\mathrm{PSF}}$')
            C3_line = mlines.Line2D([0], [0],color='C3', label=r'$s= 3\sigma_{\mathrm{PSF}}$')
            C4_line = mlines.Line2D([0], [0],color='C4', label=r'$s= 4\sigma_{\mathrm{PSF}}$')
            C8_line = mlines.Line2D([0], [0],color='C5', label=r'$s = 5\sigma_{\mathrm{PSF}}$')
            C9_line = mlines.Line2D([0], [0],color='k', label=r'SMLM (Uniform' + '\n' + 'illumination intensity)')
            C10_line = mlines.Line2D([0], [0],color='k', linestyle='--', label=r'Localization on ISM' + '\n' + 'reconstruction data')
            C11_line = mlines.Line2D([0], [0],color='k', linestyle=':', label=r'Localization on' + '\n' + 'Fourier reweighted ISM' + '\n' + 'reconstruction data')
    
            ax.legend(handles=[C0_line, C1_line, C2_line, C3_line, C4_line, C8_line, C9_line, C10_line, C11_line], prop={'size':13}, loc='upper left', bbox_to_anchor=(1.04,1), labelspacing = 1)

        plt.savefig(('./Figures/4de-S6/'+str(int(panel))) + '.svg')

    if plotnumber == '4bc-heatmap':
        fig, ax = plt.subplots(1,1)
        fig.set_size_inches(5, 5, forward=True)
        
        xnum=51
        ynum=51
        CRLBx = np.zeros((xnum, ynum))
        for yindex in range(ynum):
            CRLBx[yindex,:], _, _, _ = load_4(('xp_s_2pattern_heatmap_' + str(int(yindex)) + '_fb'))
        CRLBx_SMLM, _, _, _ = load_4('SMLM')

        if panel == 1:
            divnorm = mcolors.LogNorm(vmin=np.min(CRLBx*10**9), vmax=np.max(CRLBx*10**9))
            img = ax.imshow(CRLBx*10**9, cmap='viridis', interpolation='catrom', aspect='equal', extent = (-2*dx_i*10**9, 2*dx_i*10**9, -2*dx_i*10**9, 2*dx_i*10**9), norm = divnorm)
            ax.set_xlabel('Emitter-focus distance ' + r'$(x_f - \theta_x)$' + ' (nm)')
            ax.set_ylabel('Emitter-focus distance ' + r'$(y_f - \theta_y)$' + ' (nm)')
            ax.tick_params(axis='both', which='major')
            ax.tick_params(axis='both', which='minor') 
            
            plt.savefig(('./Figures/4bc-heatmap/heatmap_'+str(int(panel))) + '.svg')
            
            fig, ax = plt.subplots(1,1)
            fig.set_size_inches(5, 5, forward=True)
            cbar = plt.colorbar(img)
            cbar.set_label(r'$x$'+'-Cramér-Rao lower bound (nm)')
            
            plt.savefig(('./Figures/4bc-heatmap/colorbar_heatmap_'+str(int(panel))) + '.svg')
            
        if panel == 2: 
            divnorm = mcolors.TwoSlopeNorm(vcenter=1)
            img = ax.imshow(CRLBx_SMLM/CRLBx, cmap='RdYlGn', interpolation='catrom', aspect='equal', extent = (-2*dx_i*10**9, 2*dx_i*10**9, -2*dx_i*10**9, 2*dx_i*10**9), norm=divnorm)
            ax.set_xlabel('Emitter-focus distance ' + r'$(x_f - \theta_x)$' + ' (nm)')
            ax.set_ylabel('Emitter-focus distance ' + r'$(y_f - \theta_y)$' + ' (nm)')
            ax.tick_params(axis='both', which='major')
            ax.tick_params(axis='both', which='minor') 
            
            plt.savefig(('./Figures/4bc-heatmap/heatmap_'+str(int(panel))) + '.svg')
            
            fig, ax = plt.subplots(1,1)
            fig.set_size_inches(5, 5, forward=True)
            cbar = plt.colorbar(img)
            cbar.set_label('Improvement over SMLM in ' + r'$x$'+ '-direction')
            
            plt.savefig(('./Figures/4bc-heatmap/colorbar_heatmap_'+str(int(panel))) + '.svg')    

    if plotnumber == '4ij-S10':
        CRLBx, musum, signalsum, Bisum = load_4('xp_r_triangulation_rot_nocenter_ext_fb')
        CRLBx_SMLM, musum_SMLM, signalsum_SMLM, Bisum_SMLM = load_4('SMLM')
        s_list, xrange = load_xp_r('xp_r_triangulation_rot_nocenter_ext_fb')
        
        fig, ax = plt.subplots(1,1)
        fig.set_size_inches(5, 5, forward=True)
        
        if panel == 1:
            ax.semilogy((xrange-theta[0])/(10**-9),CRLBx.T[:,0:4]*10**9)
            ax.semilogy((xrange-theta[0])/(10**-9),CRLBx_SMLM*10**9*np.ones(np.size(xrange)), color='k')
            ax.semilogy((xrange-theta[0])/(10**-9),CRLBx_ISM*10**9*np.ones(np.size(xrange)), color='k', linestyle='--')
            ax.semilogy((xrange-theta[0])/(10**-9),CRLBx_FR*10**9*np.ones(np.size(xrange)), color='k', linestyle=':')
            ax.set_xlabel('Emitter-focus distance ' + r'$(x_f - \theta_x)$' + ' (nm)')
            ax.set_ylabel(r'$x$'+'-Cramér-Rao lower bound (nm)')
            ax.grid()
            ax.set_xlim(left=-130, right = 130)
            ax.set_ylim(bottom = 0.6, top = 60)
            ax.tick_params(axis='both', which='major')
            ax.tick_params(axis='both', which='minor')
            
        if panel == 2:
            ax.plot((xrange-theta[0])/(10**-9),CRLBx_SMLM/CRLBx.T[:,0:4])
            ax.plot((xrange-theta[0])/(10**-9),imp_ISM*np.ones(np.size(xrange)), color='k', linestyle='--')
            ax.plot((xrange-theta[0])/(10**-9),imp_FR*np.ones(np.size(xrange)), color='k', linestyle=':')
            ax.set_xlabel('Emitter-focus distance ' + r'$(x_f - \theta_x)$' + ' (nm)')
            ax.set_ylabel('Improvement over SMLM in ' + r'$x$'+ '-direction')
            ax.grid()
            ax.set_xlim(left=-130, right = 130)
            ax.set_ylim(bottom = 0, top = 3.5)
            ax.tick_params(axis='both', which='major')
            ax.tick_params(axis='both', which='minor')                     

        if panel == 3:
            ax.plot((xrange-theta[0])/(10**-9),musum.T[:,0:4])
            ax.plot((xrange-theta[0])/(10**-9),musum_SMLM*np.ones(np.size(xrange)), color='k')
            ax.set_xlabel('Emitter-focus distance ' + r'$(x_f - \theta_x)$' + ' (nm)')
            ax.set_ylabel('Expected amount of photons')
            ax.grid()
            ax.set_xlim(left=-130, right = 130)
            ax.set_ylim(bottom = 0, top = 7100)
            ax.tick_params(axis='both', which='major')
            ax.tick_params(axis='both', which='minor')                  

        if panel == 4:
            ax.plot((xrange-theta[0])/(10**-9),signalsum.T[:,0:4])
            ax.plot((xrange-theta[0])/(10**-9),signalsum_SMLM*np.ones(np.size(xrange)), color='k')
            ax.set_xlabel('Emitter-focus distance ' + r'$(x_f - \theta_x)$' + ' (nm)')
            ax.set_ylabel('Illumination-intensity normalized' + '\n' + 'amount of signal photons')
            ax.grid()
            ax.set_xlim(left=-130, right = 130)
            ax.set_ylim(bottom = 0, top = 2100)
            ax.tick_params(axis='both', which='major')
            ax.tick_params(axis='both', which='minor')   

        if panel == 5:
            ax.plot((xrange-theta[0])/(10**-9),Bisum.T[:,0:4]/100)
            ax.plot((xrange-theta[0])/(10**-9),Bisum_SMLM/100*np.ones(np.size(xrange)), color='k')
            ax.set_xlabel('Emitter-focus distance ' + r'$(x_f - \theta_x)$' + ' (nm)')
            ax.set_ylabel('Illumination-intensity normalized average' + '\n' + 'amount of background photons per pixel')
            ax.grid()
            ax.set_xlim(left=-130, right = 130)
            ax.set_ylim(bottom = 0, top = 5100)
            ax.tick_params(axis='both', which='major')
            ax.tick_params(axis='both', which='minor')

        if panel == 6:
            CRLBx, musum, signalsum, Bisum = load_4('thetaI_thetab_triangulation_rot_nocenter_fb')
            thetab_list, thetaIrange = load_thetaI_thetab('thetaI_thetab_triangulation_rot_nocenter_fb')
            
            ax.loglog(thetaIrange,CRLBx.T*10**9)
            ax.set_xlabel('Expected signal photon budget ' + r'$\theta_I$' + ' (photons)')
            ax.set_ylabel(r'$x$'+'-Cramér-Rao lower bound (nm)')
            ax.grid()
            ax.set_xlim(left=100, right = 10000)
            ax.set_ylim(bottom = 0.1, top = 40)
            ax.tick_params(axis='both', which='major')
            ax.tick_params(axis='both', which='minor')  
            
            #Legend
            C0_line = mlines.Line2D([0], [0],color='C0', label=(r'$\theta_b = 0$'))
            C1_line = mlines.Line2D([0], [0],color='C1', label=(r'$\theta_b = 1$'))
            C2_line = mlines.Line2D([0], [0],color='C2', label=(r'$\theta_b = 4$'))
            C3_line = mlines.Line2D([0], [0],color='C3', label=(r'$\theta_b = 8$'))
            C8_line = mlines.Line2D([0], [0],color='C4', label=(r'$\theta_b = 16$'))
            C9_line = mlines.Line2D([0], [0],color='None', label=('(in phot/px.)'))
    
            ax.legend(handles=[C0_line, C1_line, C2_line, C3_line, C8_line, C9_line], ncol = 2, prop={'size':13}, labelspacing = 1)  

        if panel == 0:
            #Legend
            C0_line = mlines.Line2D([0], [0],color='C0', label=r'$r = 1\sigma_{\mathrm{PSF}}$')
            C1_line = mlines.Line2D([0], [0],color='C1', label=r'$r = 2\sigma_{\mathrm{PSF}}$')
            C2_line = mlines.Line2D([0], [0],color='C2', label=r'$r = 3\sigma_{\mathrm{PSF}}$')
            C3_line = mlines.Line2D([0], [0],color='C3', label=r'$r = 4\sigma_{\mathrm{PSF}}$')
            C9_line = mlines.Line2D([0], [0],color='k', label=r'SMLM (Uniform' + '\n' + 'illumination intensity)')
            C10_line = mlines.Line2D([0], [0],color='k', linestyle='--', label=r'Localization on ISM' + '\n' + 'reconstruction data')
            C11_line = mlines.Line2D([0], [0],color='k', linestyle=':', label=r'Localization on' + '\n' + 'Fourier reweighted ISM' + '\n' + 'reconstruction data')
    
            ax.legend(handles=[C0_line, C1_line, C2_line, C3_line, C9_line, C10_line, C11_line], prop={'size':13}, loc='upper left', bbox_to_anchor=(1.04,1), labelspacing = 1)

        plt.savefig(('./Figures/4ij-S10/'+str(int(panel))) + '.svg')

    if plotnumber == '4gh-heatmap':
        fig, ax = plt.subplots(1,1)
        fig.set_size_inches(5, 5, forward=True)
        
        xnum=41
        ynum=41
        CRLBx = np.zeros((xnum, ynum))
        for yindex in range(ynum):
            CRLBx[yindex,:], _, _, _ = load_4(('xp_r_triangulation_rot_nocenter_ext_heatmap_' + str(int(yindex)) + '_fb'))
        CRLBx_SMLM, _, _, _ = load_4('SMLM')

        if panel == 1:
            divnorm = mcolors.LogNorm(vmin=np.min(CRLBx*10**9), vmax=np.max(CRLBx*10**9))
            img = ax.imshow(CRLBx*10**9, cmap='viridis', interpolation='catrom', aspect='equal', extent = (-2*dx_i*10**9, 2*dx_i*10**9, -2*dx_i*10**9, 2*dx_i*10**9), norm = divnorm)
            ax.set_xlabel('Emitter-focus distance ' + r'$(x_f - \theta_x)$' + ' (nm)')
            ax.set_ylabel('Emitter-focus distance ' + r'$(y_f - \theta_y)$' + ' (nm)')
            ax.tick_params(axis='both', which='major')
            ax.tick_params(axis='both', which='minor') 
            
            plt.savefig(('./Figures/4gh-heatmap/heatmap_'+str(int(panel))) + '.svg')
            
            fig, ax = plt.subplots(1,1)
            fig.set_size_inches(5, 5, forward=True)
            cbar = plt.colorbar(img)
            cbar.set_label(r'$x$'+'-Cramér-Rao lower bound (nm)')
            
            plt.savefig(('./Figures/4gh-heatmap/colorbar_heatmap_'+str(int(panel))) + '.svg')
            
        if panel == 2: 
            divnorm = mcolors.TwoSlopeNorm(vcenter=1)
            img = ax.imshow(CRLBx_SMLM/CRLBx, cmap='RdYlGn', interpolation='catrom', aspect='equal', extent = (-2*dx_i*10**9, 2*dx_i*10**9, -2*dx_i*10**9, 2*dx_i*10**9), norm=divnorm)
            ax.set_xlabel('Emitter-focus distance ' + r'$(x_f - \theta_x)$' + ' (nm)')
            ax.set_ylabel('Emitter-focus distance ' + r'$(y_f - \theta_y)$' + ' (nm)')
            ax.tick_params(axis='both', which='major')
            ax.tick_params(axis='both', which='minor') 
            
            plt.savefig(('./Figures/4gh-heatmap/heatmap_'+str(int(panel))) + '.svg')
            
            fig, ax = plt.subplots(1,1)
            fig.set_size_inches(5, 5, forward=True)
            cbar = plt.colorbar(img)
            cbar.set_label('Improvement over SMLM in ' + r'$x$'+ '-direction')
            
            plt.savefig(('./Figures/4gh-heatmap/colorbar_heatmap_'+str(int(panel))) + '.svg')  

    if plotnumber == '5-S31':
        CRLBx, musum, signalsum, Bisum = load_4('xp_r_triangulation_rot_doughnut_ext_fb')
        CRLBx_SMLM, musum_SMLM, signalsum_SMLM, Bisum_SMLM = load_4('SMLM')
        s_list, xrange = load_xp_r('xp_r_triangulation_rot_doughnut_ext_fb')
        
        fig, ax = plt.subplots(1,1)
        fig.set_size_inches(5, 5, forward=True)
        
        if panel == 1:
            ax.semilogy((xrange-theta[0])/(10**-9),CRLBx.T*10**9)
            ax.semilogy((xrange-theta[0])/(10**-9),CRLBx_SMLM*10**9*np.ones(np.size(xrange)), color='k')
            ax.semilogy((xrange-theta[0])/(10**-9),CRLBx_ISM*10**9*np.ones(np.size(xrange)), color='k', linestyle='--')
            ax.semilogy((xrange-theta[0])/(10**-9),CRLBx_FR*10**9*np.ones(np.size(xrange)), color='k', linestyle=':')
            ax.set_xlabel('Emitter-focus distance ' + r'$(x_f - \theta_x)$' + ' (nm)')
            ax.set_ylabel(r'$x$'+'-Cramér-Rao lower bound (nm)')
            ax.grid()
            ax.set_xlim(left=-130, right = 130)
            ax.set_ylim(bottom = 0.4, top = 60)
            ax.tick_params(axis='both', which='major')
            ax.tick_params(axis='both', which='minor')
            
        if panel == 2:
            ax.plot((xrange-theta[0])/(10**-9),CRLBx_SMLM/CRLBx.T)
            ax.plot((xrange-theta[0])/(10**-9),imp_ISM*np.ones(np.size(xrange)), color='k', linestyle='--')
            ax.plot((xrange-theta[0])/(10**-9),imp_FR*np.ones(np.size(xrange)), color='k', linestyle=':')
            ax.set_xlabel('Emitter-focus distance ' + r'$(x_f - \theta_x)$' + ' (nm)')
            ax.set_ylabel('Improvement over SMLM in ' + r'$x$'+ '-direction')
            ax.grid()
            ax.set_xlim(left=-130, right = 130)
            ax.set_ylim(bottom = 0, top = 4)
            ax.tick_params(axis='both', which='major')
            ax.tick_params(axis='both', which='minor')                     

        if panel == 3:
            ax.plot((xrange-theta[0])/(10**-9),musum.T)
            ax.plot((xrange-theta[0])/(10**-9),musum_SMLM*np.ones(np.size(xrange)), color='k')
            ax.set_xlabel('Emitter-focus distance ' + r'$(x_f - \theta_x)$' + ' (nm)')
            ax.set_ylabel('Expected amount of photons')
            ax.grid()
            ax.set_xlim(left=-130, right = 130)
            ax.set_ylim(bottom = 0, top = 2900)
            ax.tick_params(axis='both', which='major')
            ax.tick_params(axis='both', which='minor')                  

        if panel == 4:
            ax.plot((xrange-theta[0])/(10**-9),signalsum.T)
            ax.plot((xrange-theta[0])/(10**-9),signalsum_SMLM*np.ones(np.size(xrange)), color='k')
            ax.set_xlabel('Emitter-focus distance ' + r'$(x_f - \theta_x)$' + ' (nm)')
            ax.set_ylabel('Illumination-intensity normalized' + '\n' + 'amount of signal photons')
            ax.grid()
            ax.set_xlim(left=-130, right = 130)
            ax.set_ylim(bottom = 0, top = 2100)
            ax.tick_params(axis='both', which='major')
            ax.tick_params(axis='both', which='minor')   

        if panel == 5:
            ax.plot((xrange-theta[0])/(10**-9),Bisum.T/100)
            ax.plot((xrange-theta[0])/(10**-9),Bisum_SMLM/100*np.ones(np.size(xrange)), color='k')
            ax.set_xlabel('Emitter-focus distance ' + r'$(x_f - \theta_x)$' + ' (nm)')
            ax.set_ylabel('Illumination-intensity normalized average' + '\n' + 'amount of background photons per pixel')
            ax.grid()
            ax.set_xlim(left=-130, right = 130)
            ax.set_ylim(bottom = 0, top = 700)
            ax.tick_params(axis='both', which='major')
            ax.tick_params(axis='both', which='minor')   

        if panel == 6:
            CRLBx, musum, signalsum, Bisum = load_4('thetaI_thetab_triangulation_rot_doughnut_fb')
            thetab_list, thetaIrange = load_thetaI_thetab('thetaI_thetab_triangulation_rot_doughnut_fb')
            
            ax.loglog(thetaIrange,CRLBx.T*10**9)
            ax.set_xlabel('Expected signal photon budget ' + r'$\theta_I$' + ' (photons)')
            ax.set_ylabel(r'$x$'+'-Cramér-Rao lower bound (nm)')
            ax.grid()
            ax.set_xlim(left=100, right = 10000)
            ax.set_ylim(bottom = 0.1, top = 40)
            ax.tick_params(axis='both', which='major')
            ax.tick_params(axis='both', which='minor')  
            
            #Legend
            C0_line = mlines.Line2D([0], [0],color='C0', label=(r'$\theta_b = 0$'))
            C1_line = mlines.Line2D([0], [0],color='C1', label=(r'$\theta_b = 1$'))
            C2_line = mlines.Line2D([0], [0],color='C2', label=(r'$\theta_b = 4$'))
            C3_line = mlines.Line2D([0], [0],color='C3', label=(r'$\theta_b = 8$'))
            C8_line = mlines.Line2D([0], [0],color='C4', label=(r'$\theta_b = 16$'))
            C9_line = mlines.Line2D([0], [0],color='None', label=('(in phot/px.)'))
    
            ax.legend(handles=[C0_line, C1_line, C2_line, C3_line, C8_line, C9_line], ncol = 2, prop={'size':13}, labelspacing = 1) 

        if panel == 0:
            #Legend
            C0_line = mlines.Line2D([0], [0],color='C0', label=r'$r = 0.5\sigma_{\mathrm{PSF}}$')
            C1_line = mlines.Line2D([0], [0],color='C1', label=r'$r = 1\sigma_{\mathrm{PSF}}$')
            C2_line = mlines.Line2D([0], [0],color='C2', label=r'$r = 2\sigma_{\mathrm{PSF}}$')
            C3_line = mlines.Line2D([0], [0],color='C3', label=r'$r = 3\sigma_{\mathrm{PSF}}$')
            C4_line = mlines.Line2D([0], [0],color='C4', label=r'$r = 4\sigma_{\mathrm{PSF}}$')
            C9_line = mlines.Line2D([0], [0],color='k', label=r'SMLM (Uniform' + '\n' + 'illumination intensity)')
            C10_line = mlines.Line2D([0], [0],color='k', linestyle='--', label=r'Localization on ISM' + '\n' + 'reconstruction data')
            C11_line = mlines.Line2D([0], [0],color='k', linestyle=':', label=r'Localization on' + '\n' + 'Fourier reweighted ISM' + '\n' + 'reconstruction data')
    
            ax.legend(handles=[C0_line, C1_line, C2_line, C3_line, C4_line, C9_line, C10_line, C11_line], prop={'size':13}, loc='upper left', bbox_to_anchor=(1.04,1), labelspacing = 1)

        plt.savefig(('./Figures/5-S31/'+str(int(panel))) + '.svg')
            
    if plotnumber == '5-heatmap':
        fig, ax = plt.subplots(1,1)
        fig.set_size_inches(5, 5, forward=True)
        
        xnum=35
        ynum=35
        CRLBx = np.zeros((xnum, ynum))
        for yindex in range(ynum):
            CRLBx[yindex,:], _, _, _ = load_4(('xp_r_triangulation_rot_doughnut_ext_heatmap_' + str(int(yindex)) + '_fb'))
        CRLBx_SMLM, _, _, _ = load_4('SMLM')

        if panel == 1:
            divnorm = mcolors.LogNorm(vmin=np.min(CRLBx*10**9), vmax=np.max(CRLBx*10**9))
            img = ax.imshow(CRLBx*10**9, cmap='viridis', interpolation='catrom', aspect='equal', extent = (-2*dx_i*10**9, 2*dx_i*10**9, -2*dx_i*10**9, 2*dx_i*10**9), norm = divnorm)
            ax.set_xlabel('Emitter-focus distance ' + r'$(x_f - \theta_x)$' + ' (nm)')
            ax.set_ylabel('Emitter-focus distance ' + r'$(y_f - \theta_y)$' + ' (nm)')
            ax.tick_params(axis='both', which='major')
            ax.tick_params(axis='both', which='minor') 
            
            plt.savefig(('./Figures/5-heatmap/heatmap_'+str(int(panel))) + '.svg')
            
            fig, ax = plt.subplots(1,1)
            fig.set_size_inches(5, 5, forward=True)
            cbar = plt.colorbar(img)
            cbar.set_label(r'$x$'+'-Cramér-Rao lower bound (nm)')
            
            plt.savefig(('./Figures/5-heatmap/colorbar_heatmap_'+str(int(panel))) + '.svg')
            
        if panel == 2: 
            divnorm = mcolors.TwoSlopeNorm(vcenter=1)
            img = ax.imshow(CRLBx_SMLM/CRLBx, cmap='RdYlGn', interpolation='catrom', aspect='equal', extent = (-2*dx_i*10**9, 2*dx_i*10**9, -2*dx_i*10**9, 2*dx_i*10**9), norm=divnorm)
            ax.set_xlabel('Emitter-focus distance ' + r'$(x_f - \theta_x)$' + ' (nm)')
            ax.set_ylabel('Emitter-focus distance ' + r'$(y_f - \theta_y)$' + ' (nm)')
            ax.tick_params(axis='both', which='major')
            ax.tick_params(axis='both', which='minor') 
            
            plt.savefig(('./Figures/5-heatmap/heatmap_'+str(int(panel))) + '.svg')
            
            fig, ax = plt.subplots(1,1)
            fig.set_size_inches(5, 5, forward=True)
            cbar = plt.colorbar(img)
            cbar.set_label('Improvement over SMLM in ' + r'$x$'+ '-direction')
            
            plt.savefig(('./Figures/5-heatmap/colorbar_heatmap_'+str(int(panel))) + '.svg')  

    if plotnumber == 'S2':       
        fig, ax = plt.subplots(1,1)
        fig.set_size_inches(5, 5, forward=True)

        if panel == 1:
            CRLBx, musum, signalsum, Bisum = load_4('thetaI_Nm_1pattern_fb')
            CRLBx_SMLM, musum_SMLM, signalsum_SMLM, Bisum_SMLM = load_4('SMLM')
            s_list, xrange = load_thetaI_Nm('thetaI_Nm_1pattern_fb')
            CRLBx_gt, musum_gt, signalsum_gt, Bisum_gt = load_4('thetaI_Nm_1pattern_gt_fb')
            gt = np.repeat([CRLBx_gt[:,-1]],20,axis=0)
            
            ax.plot(xrange,np.abs((CRLBx.T-gt))/np.abs(gt)*100)
            ax.set_xlabel('Amount of mesh pixels per direction')
            ax.set_ylabel('Relative error of ' + r'$x$' +'-CRLB (%)')
            ax.grid()
            ax.set_xlim(left=10, right = 200)
            ax.set_ylim(bottom = 0, top = 1.5)
            ax.tick_params(axis='both', which='major')
            ax.tick_params(axis='both', which='minor')
            
        if panel == 2:
            CRLBx, musum, signalsum, Bisum = load_4('thetaI_Nm_2pattern_fb')
            CRLBx_SMLM, musum_SMLM, signalsum_SMLM, Bisum_SMLM = load_4('SMLM')
            s_list, xrange = load_thetaI_Nm('thetaI_Nm_2pattern_fb')
            CRLBx_gt, musum_gt, signalsum_gt, Bisum_gt = load_4('thetaI_Nm_2pattern_gt_fb')
            gt = np.repeat([CRLBx_gt[:,-1]],20,axis=0)
            
            ax.plot(xrange,np.abs((CRLBx.T-gt))/np.abs(gt)*100)
            ax.set_xlabel('Amount of mesh pixels per direction')
            ax.set_ylabel('Relative error of ' + r'$x$' +'-CRLB (%)')
            ax.grid()
            ax.set_xlim(left=10, right = 200)
            ax.set_ylim(bottom = 0, top = 1.5)
            ax.tick_params(axis='both', which='major')
            ax.tick_params(axis='both', which='minor')            

        if panel == 3:
            CRLBx, musum, signalsum, Bisum = load_4('thetaI_Nm_triangulation_rot_nocenter_fb')
            CRLBx_SMLM, musum_SMLM, signalsum_SMLM, Bisum_SMLM = load_4('SMLM')
            s_list, xrange = load_thetaI_Nm('thetaI_Nm_triangulation_rot_nocenter_fb')
            CRLBx_gt, musum_gt, signalsum_gt, Bisum_gt = load_4('thetaI_Nm_triangulation_rot_nocenter_gt_fb')
            gt = np.repeat([CRLBx_gt[:,-1]],20,axis=0)
            
            ax.plot(xrange,np.abs((CRLBx.T-gt))/np.abs(gt)*100)
            ax.set_xlabel('Amount of mesh pixels per direction')
            ax.set_ylabel('Relative error of ' + r'$x$' +'-CRLB (%)')
            ax.grid()
            ax.set_xlim(left=10, right = 200)
            ax.set_ylim(bottom = 0, top = 1.5)
            ax.tick_params(axis='both', which='major')
            ax.tick_params(axis='both', which='minor')   
        
        if panel == 4:
            CRLBx, musum, signalsum, Bisum = load_4('thetaI_Nm_1pattern')
            CRLBx_SMLM, musum_SMLM, signalsum_SMLM, Bisum_SMLM = load_4('SMLM')
            s_list, xrange = load_thetaI_Nm('thetaI_Nm_1pattern')
            CRLBx_gt, musum_gt, signalsum_gt, Bisum_gt = load_4('thetaI_Nm_1pattern_gt')
            gt = np.repeat([CRLBx_gt[:,-1]],20,axis=0)
            
            ax.plot(xrange,np.abs((CRLBx.T-gt))/np.abs(gt)*100)
            ax.set_xlabel('Amount of mesh pixels per direction')
            ax.set_ylabel('Relative error of ' + r'$x$' +'-CRLB (%)')
            ax.grid()
            ax.set_xlim(left=10, right = 200)
            ax.set_ylim(bottom = 0, top = 1.5)
            ax.tick_params(axis='both', which='major')
            ax.tick_params(axis='both', which='minor')
            
        if panel == 5:
            CRLBx, musum, signalsum, Bisum = load_4('thetaI_Nm_2pattern')
            CRLBx_SMLM, musum_SMLM, signalsum_SMLM, Bisum_SMLM = load_4('SMLM')
            s_list, xrange = load_thetaI_Nm('thetaI_Nm_2pattern')
            CRLBx_gt, musum_gt, signalsum_gt, Bisum_gt = load_4('thetaI_Nm_2pattern_gt')
            gt = np.repeat([CRLBx_gt[:,-1]],20,axis=0)
            
            ax.plot(xrange,np.abs((CRLBx.T-gt))/np.abs(gt)*100)
            ax.set_xlabel('Amount of mesh pixels per direction')
            ax.set_ylabel('Relative error of ' + r'$x$' +'-CRLB (%)')
            ax.grid()
            ax.set_xlim(left=10, right = 200)
            ax.set_ylim(bottom = 0, top = 1.5)
            ax.tick_params(axis='both', which='major')
            ax.tick_params(axis='both', which='minor')            

        if panel == 6:
            CRLBx, musum, signalsum, Bisum = load_4('thetaI_Nm_triangulation_rot_nocenter')
            CRLBx_SMLM, musum_SMLM, signalsum_SMLM, Bisum_SMLM = load_4('SMLM')
            s_list, xrange = load_thetaI_Nm('thetaI_Nm_triangulation_rot_nocenter')
            CRLBx_gt, musum_gt, signalsum_gt, Bisum_gt = load_4('thetaI_Nm_triangulation_rot_nocenter_gt')
            gt = np.repeat([CRLBx_gt[:,-1]],20,axis=0)
            
            ax.plot(xrange,np.abs((CRLBx.T-gt))/np.abs(gt)*100)
            ax.set_xlabel('Amount of mesh pixels per direction')
            ax.set_ylabel('Relative error of ' + r'$x$' +'-CRLB (%)')
            ax.grid()
            ax.set_xlim(left=10, right = 200)
            ax.set_ylim(bottom = 0, top = 1.5)
            ax.tick_params(axis='both', which='major')
            ax.tick_params(axis='both', which='minor')                

        if panel == 7:
            CRLBx, musum, signalsum, Bisum = load_4('thetaI_Nm_1pattern_adaptedbg')
            CRLBx_SMLM, musum_SMLM, signalsum_SMLM, Bisum_SMLM = load_4('SMLM')
            s_list, xrange = load_thetaI_Nm('thetaI_Nm_1pattern_adaptedbg')
            CRLBx_gt, musum_gt, signalsum_gt, Bisum_gt = load_4('thetaI_Nm_1pattern_gt_adaptedbg')
            gt = np.repeat([CRLBx_gt[:,-1]],20,axis=0)
            
            ax.plot(xrange,np.abs((CRLBx.T-gt))/np.abs(gt)*100)
            ax.set_xlabel('Amount of mesh pixels per direction')
            ax.set_ylabel('Relative error of ' + r'$x$' +'-CRLB (%)')
            ax.grid()
            ax.set_xlim(left=10, right = 200)
            ax.set_ylim(bottom = 0, top = 1.5)
            ax.tick_params(axis='both', which='major')
            ax.tick_params(axis='both', which='minor')
            
        if panel == 8:
            CRLBx, musum, signalsum, Bisum = load_4('thetaI_Nm_2pattern_adaptedbg')
            CRLBx_SMLM, musum_SMLM, signalsum_SMLM, Bisum_SMLM = load_4('SMLM')
            s_list, xrange = load_thetaI_Nm('thetaI_Nm_2pattern_adaptedbg')
            CRLBx_gt, musum_gt, signalsum_gt, Bisum_gt = load_4('thetaI_Nm_2pattern_gt_adaptedbg')
            gt = np.repeat([CRLBx_gt[:,-1]],20,axis=0)
            
            ax.plot(xrange,np.abs((CRLBx.T-gt))/np.abs(gt)*100)
            ax.set_xlabel('Amount of mesh pixels per direction')
            ax.set_ylabel('Relative error of ' + r'$x$' +'-CRLB (%)')
            ax.grid()
            ax.set_xlim(left=10, right = 200)
            ax.set_ylim(bottom = 0, top = 1.5)
            ax.tick_params(axis='both', which='major')
            ax.tick_params(axis='both', which='minor')            

        if panel == 9:
            CRLBx, musum, signalsum, Bisum = load_4('thetaI_Nm_triangulation_rot_nocenter_adaptedbg')
            CRLBx_SMLM, musum_SMLM, signalsum_SMLM, Bisum_SMLM = load_4('SMLM')
            s_list, xrange = load_thetaI_Nm('thetaI_Nm_triangulation_rot_nocenter_adaptedbg')
            CRLBx_gt, musum_gt, signalsum_gt, Bisum_gt = load_4('thetaI_Nm_triangulation_rot_nocenter_gt_adaptedbg')
            gt = np.repeat([CRLBx_gt[:,-1]],20,axis=0)
            
            ax.plot(xrange,np.abs((CRLBx.T-gt))/np.abs(gt)*100)
            ax.set_xlabel('Amount of mesh pixels per direction')
            ax.set_ylabel('Relative error of ' + r'$x$' +'-CRLB (%)')
            ax.grid()
            ax.set_xlim(left=10, right = 200)
            ax.set_ylim(bottom = 0, top = 1.5)
            ax.tick_params(axis='both', which='major')
            ax.tick_params(axis='both', which='minor')   

        if panel == 0:
            #Legend
            C0_line = mlines.Line2D([0], [0],color='C0', label=r'$\theta_I = 500$' + ' photons')
            C1_line = mlines.Line2D([0], [0],color='C1', label=r'$\theta_I = 1000$' + ' photons')
            C2_line = mlines.Line2D([0], [0],color='C2', label=r'$\theta_I = 2000$' + ' photons')
            C3_line = mlines.Line2D([0], [0],color='C3', label=r'$\theta_I = 5000$' + ' photons')
            ax.legend(handles=[C0_line, C1_line, C2_line, C3_line], prop={'size':13}, loc='upper left', bbox_to_anchor=(1.04,1), labelspacing = 1)

        plt.savefig(('./Figures/S2/'+str(int(panel))) + '.svg')

    if plotnumber == 'S5':
        CRLBx, musum, signalsum, Bisum = load_4('yp_rp_1pattern_fb')
        CRLBx_SMLM, musum_SMLM, signalsum_SMLM, Bisum_SMLM = load_4('SMLM')
        r_p_list, yrange = load_yp_rp('yp_rp_1pattern_fb')
        
        fig, ax = plt.subplots(1,1)
        fig.set_size_inches(5, 5, forward=True)
        
        if panel == 1:
            ax.semilogy((yrange-theta[1])/(10**-9),CRLBx.T[:,0:4]*10**9)
            ax.semilogy((yrange-theta[1])/(10**-9),CRLBx.T[:,-1]*10**9)
            ax.semilogy((yrange-theta[1])/(10**-9),CRLBx_SMLM*10**9*np.ones(np.size(yrange)), color='k')
            ax.semilogy((yrange-theta[1])/(10**-9),CRLBx_ISM*10**9*np.ones(np.size(yrange)), color='k', linestyle='--')
            ax.semilogy((yrange-theta[1])/(10**-9),CRLBx_FR*10**9*np.ones(np.size(yrange)), color='k', linestyle=':')
            ax.set_xlabel('Emitter-pinhole distance ' + r'$(y_p - \theta_y)$' + ' (nm)')
            ax.set_ylabel(r'$x$'+'-Cramér-Rao lower bound (nm)')
            ax.grid()
            ax.set_xlim(left=-130, right = 130)
            ax.set_ylim(bottom = 0.6, top = 60)
            ax.tick_params(axis='both', which='major')
            ax.tick_params(axis='both', which='minor')
            
        if panel == 2:
            ax.plot((yrange-theta[1])/(10**-9),CRLBx_SMLM/CRLBx.T[:,0:4])
            ax.plot((yrange-theta[1])/(10**-9),CRLBx_SMLM/CRLBx.T[:,-1])
            ax.plot((yrange-theta[1])/(10**-9),imp_ISM*np.ones(np.size(yrange)), color='k', linestyle='--')
            ax.plot((yrange-theta[1])/(10**-9),imp_FR*np.ones(np.size(yrange)), color='k', linestyle=':')
            ax.set_xlabel('Emitter-pinhole distance ' + r'$(y_p - \theta_y)$' + ' (nm)')
            ax.set_ylabel('Improvement over SMLM in ' + r'$x$'+ '-direction')
            ax.grid()
            ax.set_xlim(left=-130, right = 130)
            ax.set_ylim(bottom = 0, top = 3.5)
            ax.tick_params(axis='both', which='major')
            ax.tick_params(axis='both', which='minor')                     

        if panel == 3:
            ax.plot((yrange-theta[1])/(10**-9),musum.T[:,0:4])
            ax.plot((yrange-theta[1])/(10**-9),musum.T[:,-1])
            ax.plot((yrange-theta[1])/(10**-9),musum_SMLM*np.ones(np.size(yrange)), color='k')
            ax.set_xlabel('Emitter-pinhole distance ' + r'$(y_p - \theta_y)$' + ' (nm)')
            ax.set_ylabel('Expected amount of photons')
            ax.grid()
            ax.set_xlim(left=-130, right = 130)
            ax.set_ylim(bottom = 0, top = 2900)
            ax.tick_params(axis='both', which='major')
            ax.tick_params(axis='both', which='minor')                  

        if panel == 4:
            ax.plot((yrange-theta[1])/(10**-9),signalsum.T[:,0:4])
            ax.plot((yrange-theta[1])/(10**-9),signalsum.T[:,-1])
            ax.plot((yrange-theta[1])/(10**-9),signalsum_SMLM*np.ones(np.size(yrange)), color='k')
            ax.set_xlabel('Emitter-pinhole distance ' + r'$(y_p - \theta_y)$' + ' (nm)')
            ax.set_ylabel('Illumination-intensity normalized' + '\n' + 'amount of signal photons')
            ax.grid()
            ax.set_xlim(left=-130, right = 130)
            ax.set_ylim(bottom = 0, top = 2100)
            ax.tick_params(axis='both', which='major')
            ax.tick_params(axis='both', which='minor')   

        if panel == 5:
            ax.plot((yrange-theta[1])/(10**-9),Bisum.T[:,0:4]/100)
            ax.plot((yrange-theta[1])/(10**-9),Bisum.T[:,-1]/100)
            ax.plot((yrange-theta[1])/(10**-9),Bisum_SMLM*np.ones(np.size(yrange)), color='k')
            ax.set_xlabel('Emitter-pinhole distance ' + r'$(y_p - \theta_y)$' + ' (nm)')
            ax.set_ylabel('Illumination-intensity normalized average' + '\n' + 'amount of background photons per pixel')
            ax.grid()
            ax.set_xlim(left=-130, right = 130)
            ax.set_ylim(bottom = 0, top = 9)
            ax.tick_params(axis='both', which='major')
            ax.tick_params(axis='both', which='minor')  
        
        if panel == 0:
            #Legend
            C0_line = mlines.Line2D([0], [0],color='C0', label=r'$r_p = 1\sigma_{\mathrm{PSF}}$')
            C1_line = mlines.Line2D([0], [0],color='C1', label=r'$r_p = 2\sigma_{\mathrm{PSF}}$')
            C2_line = mlines.Line2D([0], [0],color='C2', label=r'$r_p = 3\sigma_{\mathrm{PSF}}$')
            C3_line = mlines.Line2D([0], [0],color='C3', label=r'$r_p = 4\sigma_{\mathrm{PSF}}$')
            C8_line = mlines.Line2D([0], [0],color='C4', label=r'$r_p = \infty$')
            C9_line = mlines.Line2D([0], [0],color='k', label=r'SMLM (Uniform' + '\n' + 'illumination intensity)')
            C10_line = mlines.Line2D([0], [0],color='k', linestyle='--', label=r'Localization on ISM' + '\n' + 'reconstruction data')
            C11_line = mlines.Line2D([0], [0],color='k', linestyle=':', label=r'Localization on' + '\n' + 'Fourier reweighted ISM' + '\n' + 'reconstruction data')
    
            ax.legend(handles=[C0_line, C1_line, C2_line, C3_line, C8_line, C9_line, C10_line, C11_line], prop={'size':13}, loc='upper left', bbox_to_anchor=(1.04,1), labelspacing = 1)

        plt.savefig(('./Figures/S5/'+str(int(panel))) + '.svg')

    if plotnumber == 'S7':
        CRLBx, musum, signalsum, Bisum = load_4('yp_xs_2pattern_fb')
        CRLBx_SMLM, musum_SMLM, signalsum_SMLM, Bisum_SMLM = load_4('SMLM')
        s_list, xrange = load_xp_s('yp_xs_2pattern_fb')
        
        fig, ax = plt.subplots(1,1)
        fig.set_size_inches(5, 5, forward=True)
        
        if panel == 1:
            ax.semilogy((xrange-theta[0])/(10**-9),CRLBx.T*10**9)
            ax.semilogy((xrange-theta[0])/(10**-9),CRLBx_SMLM*10**9*np.ones(np.size(xrange)), color='k')
            ax.semilogy((xrange-theta[0])/(10**-9),CRLBx_ISM*10**9*np.ones(np.size(xrange)), color='k', linestyle='--')
            ax.semilogy((xrange-theta[0])/(10**-9),CRLBx_FR*10**9*np.ones(np.size(xrange)), color='k', linestyle=':')
            ax.set_xlabel('Emitter-focus distance ' + r'$(y_f - \theta_y)$' + ' (nm)')
            ax.set_ylabel(r'$x$'+'-Cramér-Rao lower bound (nm)')
            ax.grid()
            ax.set_xlim(left=-130, right = 130)
            ax.set_ylim(bottom = 0.6, top = 60)
            ax.tick_params(axis='both', which='major')
            ax.tick_params(axis='both', which='minor')
            
        if panel == 2:
            ax.plot((xrange-theta[0])/(10**-9),CRLBx_SMLM/CRLBx.T)
            ax.plot((xrange-theta[0])/(10**-9),imp_ISM*np.ones(np.size(xrange)), color='k', linestyle='--')
            ax.plot((xrange-theta[0])/(10**-9),imp_FR*np.ones(np.size(xrange)), color='k', linestyle=':')
            ax.set_xlabel('Emitter-focus distance ' + r'$(y_f - \theta_y)$' + ' (nm)')
            ax.set_ylabel('Improvement over SMLM in ' + r'$x$'+ '-direction')
            ax.grid()
            ax.set_xlim(left=-130, right = 130)
            ax.set_ylim(bottom = 0, top = 3.5)
            ax.tick_params(axis='both', which='major')
            ax.tick_params(axis='both', which='minor')                     

        if panel == 3:
            ax.plot((xrange-theta[0])/(10**-9),musum.T)
            ax.plot((xrange-theta[0])/(10**-9),musum_SMLM*np.ones(np.size(xrange)), color='k')
            ax.set_xlabel('Emitter-focus distance ' + r'$(y_f - \theta_y)$' + ' (nm)')
            ax.set_ylabel('Expected amount of photons')
            ax.grid()
            ax.set_xlim(left=-130, right = 130)
            ax.set_ylim(bottom = 0, top = 9500)
            ax.tick_params(axis='both', which='major')
            ax.tick_params(axis='both', which='minor')                  

        if panel == 4:
            ax.plot((xrange-theta[0])/(10**-9),signalsum.T)
            ax.plot((xrange-theta[0])/(10**-9),signalsum_SMLM*np.ones(np.size(xrange)), color='k')
            ax.set_xlabel('Emitter-focus distance ' + r'$(y_f - \theta_y)$' + ' (nm)')
            ax.set_ylabel('Illumination-intensity normalized' + '\n' + 'amount of signal photons')
            ax.grid()
            ax.set_xlim(left=-130, right = 130)
            ax.set_ylim(bottom = 0, top = 2100)
            ax.tick_params(axis='both', which='major')
            ax.tick_params(axis='both', which='minor')   

        if panel == 5:
            ax.plot((xrange-theta[0])/(10**-9),Bisum.T/100)
            ax.plot((xrange-theta[0])/(10**-9),Bisum_SMLM/100*np.ones(np.size(xrange)), color='k')
            ax.set_xlabel('Emitter-focus distance ' + r'$(y_f - \theta_y)$' + ' (nm)')
            ax.set_ylabel('Illumination-intensity normalized average' + '\n' + 'amount of background photons per pixel')
            ax.grid()
            ax.set_xlim(left=-130, right = 130)
            ax.set_ylim(bottom = 0, top = 85)
            ax.tick_params(axis='both', which='major')
            ax.tick_params(axis='both', which='minor')   

        if panel == 0:
            #Legend
            C0_line = mlines.Line2D([0], [0],color='C0', label=r'$s = 0\sigma_{\mathrm{PSF}}$')
            C1_line = mlines.Line2D([0], [0],color='C1', label=r'$s= 1\sigma_{\mathrm{PSF}}$')
            C2_line = mlines.Line2D([0], [0],color='C2', label=r'$s= 2\sigma_{\mathrm{PSF}}$')
            C3_line = mlines.Line2D([0], [0],color='C3', label=r'$s= 3\sigma_{\mathrm{PSF}}$')
            C4_line = mlines.Line2D([0], [0],color='C4', label=r'$s= 4\sigma_{\mathrm{PSF}}$')
            C8_line = mlines.Line2D([0], [0],color='C5', label=r'$s = 5\sigma_{\mathrm{PSF}}$')
            C9_line = mlines.Line2D([0], [0],color='k', label=r'SMLM (Uniform' + '\n' + 'illumination intensity)')
            C9_line = mlines.Line2D([0], [0],color='k', label=r'SMLM (Uniform' + '\n' + 'illumination intensity)')
            C10_line = mlines.Line2D([0], [0],color='k', linestyle='--', label=r'Localization on ISM' + '\n' + 'reconstruction data')
            C11_line = mlines.Line2D([0], [0],color='k', linestyle=':', label=r'Localization on' + '\n' + 'Fourier reweighted ISM' + '\n' + 'reconstruction data')
    
            ax.legend(handles=[C0_line, C1_line, C2_line, C3_line, C4_line, C8_line, C9_line, C10_line, C11_line], prop={'size':13}, loc='upper left', bbox_to_anchor=(1.04,1), labelspacing = 1)

        plt.savefig(('./Figures/S7/'+str(int(panel))) + '.svg')

    if plotnumber == 'S8':
        CRLBx, musum, signalsum, Bisum = load_4('xp_s_2pattern_nopinhole_fb')
        CRLBx_SMLM, musum_SMLM, signalsum_SMLM, Bisum_SMLM = load_4('SMLM')
        s_list, xrange = load_xp_s('xp_s_2pattern_nopinhole_fb')
        
        fig, ax = plt.subplots(1,1)
        fig.set_size_inches(5, 5, forward=True)
        
        if panel == 1:
            ax.semilogy((xrange-theta[0])/(10**-9),CRLBx.T*10**9)
            ax.semilogy((xrange-theta[0])/(10**-9),CRLBx_SMLM*10**9*np.ones(np.size(xrange)), color='k')
            ax.semilogy((xrange-theta[0])/(10**-9),CRLBx_ISM*10**9*np.ones(np.size(xrange)), color='k', linestyle='--')
            ax.semilogy((xrange-theta[0])/(10**-9),CRLBx_FR*10**9*np.ones(np.size(xrange)), color='k', linestyle=':')
            ax.set_xlabel('Emitter-focus distance ' + r'$(x_f - \theta_x)$' + ' (nm)')
            ax.set_ylabel(r'$x$'+'-Cramér-Rao lower bound (nm)')
            ax.grid()
            ax.set_xlim(left=-130, right = 130)
            ax.set_ylim(bottom = 0.6, top = 60)
            ax.tick_params(axis='both', which='major')
            ax.tick_params(axis='both', which='minor')
            
        if panel == 2:
            ax.plot((xrange-theta[0])/(10**-9),CRLBx_SMLM/CRLBx.T)
            ax.plot((xrange-theta[0])/(10**-9),imp_ISM*np.ones(np.size(xrange)), color='k', linestyle='--')
            ax.plot((xrange-theta[0])/(10**-9),imp_FR*np.ones(np.size(xrange)), color='k', linestyle=':')
            ax.set_xlabel('Emitter-focus distance ' + r'$(x_f - \theta_x)$' + ' (nm)')
            ax.set_ylabel('Improvement over SMLM in ' + r'$x$'+ '-direction')
            ax.grid()
            ax.set_xlim(left=-130, right = 130)
            ax.set_ylim(bottom = 0, top = 3.5)
            ax.tick_params(axis='both', which='major')
            ax.tick_params(axis='both', which='minor')                     

        if panel == 3:
            ax.plot((xrange-theta[0])/(10**-9),musum.T)
            ax.plot((xrange-theta[0])/(10**-9),musum_SMLM*np.ones(np.size(xrange)), color='k')
            ax.set_xlabel('Emitter-focus distance ' + r'$(x_f - \theta_x)$' + ' (nm)')
            ax.set_ylabel('Expected amount of photons')
            ax.grid()
            ax.set_xlim(left=-130, right = 130)
            ax.set_ylim(bottom = 0, top = 5000)
            ax.tick_params(axis='both', which='major')
            ax.tick_params(axis='both', which='minor')                  

        if panel == 4:
            ax.plot((xrange-theta[0])/(10**-9),signalsum.T)
            ax.plot((xrange-theta[0])/(10**-9),signalsum_SMLM*np.ones(np.size(xrange)), color='k')
            ax.set_xlabel('Emitter-focus distance ' + r'$(x_f - \theta_x)$' + ' (nm)')
            ax.set_ylabel('Illumination-intensity normalized' + '\n' + 'amount of signal photons')
            ax.grid()
            ax.set_xlim(left=-130, right = 130)
            ax.set_ylim(bottom = 0, top = 2100)
            ax.tick_params(axis='both', which='major')
            ax.tick_params(axis='both', which='minor')   

        if panel == 5:
            ax.plot((xrange-theta[0])/(10**-9),Bisum.T/100)
            ax.plot((xrange-theta[0])/(10**-9),Bisum_SMLM/100*np.ones(np.size(xrange)), color='k')
            ax.set_xlabel('Emitter-focus distance ' + r'$(x_f - \theta_x)$' + ' (nm)')
            ax.set_ylabel('Illumination-intensity normalized average' + '\n' + 'amount of background photons per pixel')
            ax.grid()
            ax.set_xlim(left=-130, right = 130)
            ax.set_ylim(bottom = 0, top = 30)
            ax.tick_params(axis='both', which='major')
            ax.tick_params(axis='both', which='minor')   

        if panel == 6:
            CRLBx, musum, signalsum, Bisum = load_4('thetaI_thetab_2pattern_nopinhole_fb')
            thetab_list, thetaIrange = load_thetaI_thetab('thetaI_thetab_2pattern_nopinhole_fb')
            
            ax.loglog(thetaIrange,CRLBx.T*10**9)
            ax.set_xlabel('Expected signal photon budget ' + r'$\theta_I$' + ' (photons)')
            ax.set_ylabel(r'$x$'+'-Cramér-Rao lower bound (nm)')
            ax.grid()
            ax.set_xlim(left=100, right = 10000)
            ax.set_ylim(bottom = 0.1, top = 40)
            ax.tick_params(axis='both', which='major')
            ax.tick_params(axis='both', which='minor')  
            
            #Legend
            C0_line = mlines.Line2D([0], [0],color='C0', label=(r'$\theta_b = 0$'))
            C1_line = mlines.Line2D([0], [0],color='C1', label=(r'$\theta_b = 1$'))
            C2_line = mlines.Line2D([0], [0],color='C2', label=(r'$\theta_b = 4$'))
            C3_line = mlines.Line2D([0], [0],color='C3', label=(r'$\theta_b = 8$'))
            C8_line = mlines.Line2D([0], [0],color='C4', label=(r'$\theta_b = 16$'))
            C9_line = mlines.Line2D([0], [0],color='None', label=('(in phot/px.)'))
    
            ax.legend(handles=[C0_line, C1_line, C2_line, C3_line, C8_line, C9_line], ncol = 2, prop={'size':13}, labelspacing = 1)  

        if panel == 0:
            #Legend
            C0_line = mlines.Line2D([0], [0],color='C0', label=r'$s = 0\sigma_{\mathrm{PSF}}$')
            C1_line = mlines.Line2D([0], [0],color='C1', label=r'$s= 1\sigma_{\mathrm{PSF}}$')
            C2_line = mlines.Line2D([0], [0],color='C2', label=r'$s= 2\sigma_{\mathrm{PSF}}$')
            C3_line = mlines.Line2D([0], [0],color='C3', label=r'$s= 3\sigma_{\mathrm{PSF}}$')
            C4_line = mlines.Line2D([0], [0],color='C4', label=r'$s= 4\sigma_{\mathrm{PSF}}$')
            C8_line = mlines.Line2D([0], [0],color='C5', label=r'$s = 5\sigma_{\mathrm{PSF}}$')
            C9_line = mlines.Line2D([0], [0],color='k', label=r'SMLM (Uniform' + '\n' + 'illumination intensity)')
            C10_line = mlines.Line2D([0], [0],color='k', linestyle='--', label=r'Localization on ISM' + '\n' + 'reconstruction data')
            C11_line = mlines.Line2D([0], [0],color='k', linestyle=':', label=r'Localization on' + '\n' + 'Fourier reweighted ISM' + '\n' + 'reconstruction data')
    
            ax.legend(handles=[C0_line, C1_line, C2_line, C3_line, C4_line, C8_line, C9_line, C10_line, C11_line], prop={'size':13}, loc='upper left', bbox_to_anchor=(1.04,1), labelspacing = 1)

        plt.savefig(('./Figures/S8/'+str(int(panel))) + '.svg')

    if plotnumber == 'S9':
        CRLBx, musum, signalsum, Bisum = load_4('yp_s_2pattern_fb')
        CRLBx_SMLM, musum_SMLM, signalsum_SMLM, Bisum_SMLM = load_4('SMLM')
        s_list, xrange = load_xp_s('yp_s_2pattern_fb')
        
        fig, ax = plt.subplots(1,1)
        fig.set_size_inches(5, 5, forward=True)
        
        if panel == 1:
            ax.semilogy((xrange-theta[0])/(10**-9),CRLBx.T*10**9)
            ax.semilogy((xrange-theta[0])/(10**-9),CRLBx_SMLM*10**9*np.ones(np.size(xrange)), color='k')
            ax.semilogy((xrange-theta[0])/(10**-9),CRLBx_ISM*10**9*np.ones(np.size(xrange)), color='k', linestyle='--')
            ax.semilogy((xrange-theta[0])/(10**-9),CRLBx_FR*10**9*np.ones(np.size(xrange)), color='k', linestyle=':')
            ax.set_xlabel('Emitter-focus distance ' + r'$(y_f - \theta_y)$' + ' (nm)')
            ax.set_ylabel(r'$x$'+'-Cramér-Rao lower bound (nm)')
            ax.grid()
            ax.set_xlim(left=-130, right = 130)
            ax.set_ylim(bottom = 0.6, top = 60)
            ax.tick_params(axis='both', which='major')
            ax.tick_params(axis='both', which='minor')
            
        if panel == 2:
            ax.plot((xrange-theta[0])/(10**-9),CRLBx_SMLM/CRLBx.T)
            ax.plot((xrange-theta[0])/(10**-9),imp_ISM*np.ones(np.size(xrange)), color='k', linestyle='--')
            ax.plot((xrange-theta[0])/(10**-9),imp_FR*np.ones(np.size(xrange)), color='k', linestyle=':')
            ax.set_xlabel('Emitter-focus distance ' + r'$(y_f - \theta_y)$' + ' (nm)')
            ax.set_ylabel('Improvement over SMLM in ' + r'$x$'+ '-direction')
            ax.grid()
            ax.set_xlim(left=-130, right = 130)
            ax.set_ylim(bottom = 0, top = 3.5)
            ax.tick_params(axis='both', which='major')
            ax.tick_params(axis='both', which='minor')                     

        if panel == 3:
            ax.plot((xrange-theta[0])/(10**-9),musum.T)
            ax.plot((xrange-theta[0])/(10**-9),musum_SMLM*np.ones(np.size(xrange)), color='k')
            ax.set_xlabel('Emitter-focus distance ' + r'$(y_f - \theta_y)$' + ' (nm)')
            ax.set_ylabel('Expected amount of photons')
            ax.grid()
            ax.set_xlim(left=-130, right = 130)
            ax.set_ylim(bottom = 0, top = 5000)
            ax.tick_params(axis='both', which='major')
            ax.tick_params(axis='both', which='minor')                  

        if panel == 4:
            ax.plot((xrange-theta[0])/(10**-9),signalsum.T)
            ax.plot((xrange-theta[0])/(10**-9),signalsum_SMLM*np.ones(np.size(xrange)), color='k')
            ax.set_xlabel('Emitter-focus distance ' + r'$(y_f - \theta_y)$' + ' (nm)')
            ax.set_ylabel('Illumination-intensity normalized' + '\n' + 'amount of signal photons')
            ax.grid()
            ax.set_xlim(left=-130, right = 130)
            ax.set_ylim(bottom = 0, top = 2100)
            ax.tick_params(axis='both', which='major')
            ax.tick_params(axis='both', which='minor')   

        if panel == 5:
            ax.plot((xrange-theta[0])/(10**-9),Bisum.T/100)
            ax.plot((xrange-theta[0])/(10**-9),Bisum_SMLM/100*np.ones(np.size(xrange)), color='k')
            ax.set_xlabel('Emitter-focus distance ' + r'$(y_f - \theta_y)$' + ' (nm)')
            ax.set_ylabel('Illumination-intensity normalized average' + '\n' + 'amount of background photons per pixel')
            ax.grid()
            ax.set_xlim(left=-130, right = 130)
            ax.set_ylim(bottom = 0, top = 30)
            ax.tick_params(axis='both', which='major')
            ax.tick_params(axis='both', which='minor')   

        if panel == 6:
            CRLBx, musum, signalsum, Bisum = load_4('thetaI_thetab_2pattern_rot_fb')
            thetab_list, thetaIrange = load_thetaI_thetab('thetaI_thetab_2pattern_rot_fb')
            
            ax.loglog(thetaIrange,CRLBx.T*10**9)
            ax.set_xlabel('Expected signal photon budget ' + r'$\theta_I$' + ' (photons)')
            ax.set_ylabel(r'$x$'+'-Cramér-Rao lower bound (nm)')
            ax.grid()
            ax.set_xlim(left=100, right = 10000)
            ax.set_ylim(bottom = 0.1, top = 40)
            ax.tick_params(axis='both', which='major')
            ax.tick_params(axis='both', which='minor')  
            
            #Legend
            C0_line = mlines.Line2D([0], [0],color='C0', label=(r'$\theta_b = 0$'))
            C1_line = mlines.Line2D([0], [0],color='C1', label=(r'$\theta_b = 1$'))
            C2_line = mlines.Line2D([0], [0],color='C2', label=(r'$\theta_b = 4$'))
            C3_line = mlines.Line2D([0], [0],color='C3', label=(r'$\theta_b = 8$'))
            C8_line = mlines.Line2D([0], [0],color='C4', label=(r'$\theta_b = 16$'))
            C9_line = mlines.Line2D([0], [0],color='None', label=('(in phot/px.)'))
    
            ax.legend(handles=[C0_line, C1_line, C2_line, C3_line, C8_line, C9_line], ncol = 2, prop={'size':13}, labelspacing = 1)  
            
        if panel == 0:
            #Legend
            C0_line = mlines.Line2D([0], [0],color='C0', label=r'$s = 0\sigma_{\mathrm{PSF}}$')
            C1_line = mlines.Line2D([0], [0],color='C1', label=r'$s= 1\sigma_{\mathrm{PSF}}$')
            C2_line = mlines.Line2D([0], [0],color='C2', label=r'$s= 2\sigma_{\mathrm{PSF}}$')
            C3_line = mlines.Line2D([0], [0],color='C3', label=r'$s= 3\sigma_{\mathrm{PSF}}$')
            C4_line = mlines.Line2D([0], [0],color='C4', label=r'$s= 4\sigma_{\mathrm{PSF}}$')
            C8_line = mlines.Line2D([0], [0],color='C5', label=r'$s = 5\sigma_{\mathrm{PSF}}$')
            C9_line = mlines.Line2D([0], [0],color='k', label=r'SMLM (Uniform' + '\n' + 'illumination intensity)')
            C10_line = mlines.Line2D([0], [0],color='k', linestyle='--', label=r'Localization on ISM' + '\n' + 'reconstruction data')
            C11_line = mlines.Line2D([0], [0],color='k', linestyle=':', label=r'Localization on' + '\n' + 'Fourier reweighted ISM' + '\n' + 'reconstruction data')
    
            ax.legend(handles=[C0_line, C1_line, C2_line, C3_line, C4_line, C8_line, C9_line, C10_line, C11_line], prop={'size':13}, loc='upper left', bbox_to_anchor=(1.04,1), labelspacing = 1)

        plt.savefig(('./Figures/S9/'+str(int(panel))) + '.svg')
    
    if plotnumber == 'S11':
        CRLBx, musum, signalsum, Bisum = load_4('xp_r_triangulation_nocenter_ext_fb')
        CRLBx_SMLM, musum_SMLM, signalsum_SMLM, Bisum_SMLM = load_4('SMLM')
        s_list, xrange = load_xp_r('xp_r_triangulation_nocenter_ext_fb')
        
        fig, ax = plt.subplots(1,1)
        fig.set_size_inches(5, 5, forward=True)
        
        if panel == 1:
            ax.semilogy((xrange-theta[0])/(10**-9),CRLBx.T[:,0:4]*10**9)
            ax.semilogy((xrange-theta[0])/(10**-9),CRLBx_SMLM*10**9*np.ones(np.size(xrange)), color='k')
            ax.semilogy((xrange-theta[0])/(10**-9),CRLBx_ISM*10**9*np.ones(np.size(xrange)), color='k', linestyle='--')
            ax.semilogy((xrange-theta[0])/(10**-9),CRLBx_FR*10**9*np.ones(np.size(xrange)), color='k', linestyle=':')
            ax.set_xlabel('Emitter-focus distance ' + r'$(x_f - \theta_x)$' + ' (nm)')
            ax.set_ylabel(r'$x$'+'-Cramér-Rao lower bound (nm)')
            ax.grid()
            ax.set_xlim(left=-130, right = 130)
            ax.set_ylim(bottom = 0.6, top = 100)
            ax.tick_params(axis='both', which='major')
            ax.tick_params(axis='both', which='minor')
            
        if panel == 2:
            ax.plot((xrange-theta[0])/(10**-9),CRLBx_SMLM/CRLBx.T[:,0:4])
            ax.plot((xrange-theta[0])/(10**-9),imp_ISM*np.ones(np.size(xrange)), color='k', linestyle='--')
            ax.plot((xrange-theta[0])/(10**-9),imp_FR*np.ones(np.size(xrange)), color='k', linestyle=':')
            ax.set_xlabel('Emitter-focus distance ' + r'$(x_f - \theta_x)$' + ' (nm)')
            ax.set_ylabel('Improvement over SMLM in ' + r'$x$'+ '-direction')
            ax.grid()
            ax.set_xlim(left=-130, right = 130)
            ax.set_ylim(bottom = 0, top = 3.5)
            ax.tick_params(axis='both', which='major')
            ax.tick_params(axis='both', which='minor')                     

        if panel == 3:
            ax.plot((xrange-theta[0])/(10**-9),musum.T[:,0:4])
            ax.plot((xrange-theta[0])/(10**-9),musum_SMLM*np.ones(np.size(xrange)), color='k')
            ax.set_xlabel('Emitter-focus distance ' + r'$(x_f - \theta_x)$' + ' (nm)')
            ax.set_ylabel('Expected amount of photons')
            ax.grid()
            ax.set_xlim(left=-130, right = 130)
            ax.tick_params(axis='both', which='major')
            ax.tick_params(axis='both', which='minor')                  

        if panel == 4:
            ax.plot((xrange-theta[0])/(10**-9),signalsum.T[:,0:4])
            ax.plot((xrange-theta[0])/(10**-9),signalsum_SMLM*np.ones(np.size(xrange)), color='k')
            ax.set_xlabel('Emitter-focus distance ' + r'$(x_f - \theta_x)$' + ' (nm)')
            ax.set_ylabel('Illumination-intensity normalized' + '\n' + 'amount of signal photons')
            ax.grid()
            ax.set_xlim(left=-130, right = 130)
            ax.set_ylim(bottom = 0, top = 2100)
            ax.tick_params(axis='both', which='major')
            ax.tick_params(axis='both', which='minor')   

        if panel == 5:
            ax.plot((xrange-theta[0])/(10**-9),Bisum.T[:,0:4]/100)
            ax.plot((xrange-theta[0])/(10**-9),Bisum_SMLM/100*np.ones(np.size(xrange)), color='k')
            ax.set_xlabel('Emitter-focus distance ' + r'$(x_f - \theta_x)$' + ' (nm)')
            ax.set_ylabel('Illumination-intensity normalized average' + '\n' + 'amount of background photons per pixel')
            ax.grid()
            ax.set_xlim(left=-130, right = 130)
            ax.set_ylim(bottom = 0, top = 5100)
            ax.tick_params(axis='both', which='major')
            ax.tick_params(axis='both', which='minor')  

        if panel == 6:
            CRLBx, musum, signalsum, Bisum = load_4('thetaI_thetab_triangulation_nocenter_fb')
            thetab_list, thetaIrange = load_thetaI_thetab('thetaI_thetab_triangulation_nocenter_fb')
            
            ax.loglog(thetaIrange,CRLBx.T*10**9)
            ax.set_xlabel('Expected signal photon budget ' + r'$\theta_I$' + ' (photons)')
            ax.set_ylabel(r'$x$'+'-Cramér-Rao lower bound (nm)')
            ax.grid()
            ax.set_xlim(left=100, right = 10000)
            ax.set_ylim(bottom = 0.1, top = 40)
            ax.tick_params(axis='both', which='major')
            ax.tick_params(axis='both', which='minor')  
            
            #Legend
            C0_line = mlines.Line2D([0], [0],color='C0', label=(r'$\theta_b = 0$'))
            C1_line = mlines.Line2D([0], [0],color='C1', label=(r'$\theta_b = 1$'))
            C2_line = mlines.Line2D([0], [0],color='C2', label=(r'$\theta_b = 4$'))
            C3_line = mlines.Line2D([0], [0],color='C3', label=(r'$\theta_b = 8$'))
            C8_line = mlines.Line2D([0], [0],color='C4', label=(r'$\theta_b = 16$'))
            C9_line = mlines.Line2D([0], [0],color='None', label=('(in phot/px.)'))
    
            ax.legend(handles=[C0_line, C1_line, C2_line, C3_line, C8_line, C9_line], ncol = 2, prop={'size':13}, labelspacing = 1)    

        if panel == 0:
            #Legend
            C0_line = mlines.Line2D([0], [0],color='C0', label=r'$r = 1\sigma_{\mathrm{PSF}}$')
            C1_line = mlines.Line2D([0], [0],color='C1', label=r'$r = 2\sigma_{\mathrm{PSF}}$')
            C2_line = mlines.Line2D([0], [0],color='C2', label=r'$r = 3\sigma_{\mathrm{PSF}}$')
            C3_line = mlines.Line2D([0], [0],color='C3', label=r'$r = 4\sigma_{\mathrm{PSF}}$')
            C9_line = mlines.Line2D([0], [0],color='k', label=r'SMLM (Uniform' + '\n' + 'illumination intensity)')
            C10_line = mlines.Line2D([0], [0],color='k', linestyle='--', label=r'Localization on ISM' + '\n' + 'reconstruction data')
            C11_line = mlines.Line2D([0], [0],color='k', linestyle=':', label=r'Localization on' + '\n' + 'Fourier reweighted ISM' + '\n' + 'reconstruction data')
    
            ax.legend(handles=[C0_line, C1_line, C2_line, C3_line, C9_line, C10_line, C11_line], prop={'size':13}, loc='upper left', bbox_to_anchor=(1.04,1), labelspacing = 1)

        plt.savefig(('./Figures/S11/'+str(int(panel))) + '.svg')

    if plotnumber == 'S12':
        CRLBx, musum, signalsum, Bisum = load_4('xp_r_triangulation_rot_fb')
        CRLBx_SMLM, musum_SMLM, signalsum_SMLM, Bisum_SMLM = load_4('SMLM')
        s_list, xrange = load_xp_r('xp_r_triangulation_rot_fb')
        
        fig, ax = plt.subplots(1,1)
        fig.set_size_inches(5, 5, forward=True)
        
        if panel == 1:
            ax.semilogy((xrange-theta[0])/(10**-9),CRLBx.T*10**9)
            ax.semilogy((xrange-theta[0])/(10**-9),CRLBx_SMLM*10**9*np.ones(np.size(xrange)), color='k')
            ax.semilogy((xrange-theta[0])/(10**-9),CRLBx_ISM*10**9*np.ones(np.size(xrange)), color='k', linestyle='--')
            ax.semilogy((xrange-theta[0])/(10**-9),CRLBx_FR*10**9*np.ones(np.size(xrange)), color='k', linestyle=':')
            ax.set_xlabel('Emitter-focus distance ' + r'$(x_f - \theta_x)$' + ' (nm)')
            ax.set_ylabel(r'$x$'+'-Cramér-Rao lower bound (nm)')
            ax.grid()
            ax.set_xlim(left=-130, right = 130)
            ax.set_ylim(bottom = 0.6, top = 60)
            ax.tick_params(axis='both', which='major')
            ax.tick_params(axis='both', which='minor')
            
        if panel == 2:
            ax.plot((xrange-theta[0])/(10**-9),CRLBx_SMLM/CRLBx.T)
            ax.plot((xrange-theta[0])/(10**-9),imp_ISM*np.ones(np.size(xrange)), color='k', linestyle='--')
            ax.plot((xrange-theta[0])/(10**-9),imp_FR*np.ones(np.size(xrange)), color='k', linestyle=':')
            ax.set_xlabel('Emitter-focus distance ' + r'$(x_f - \theta_x)$' + ' (nm)')
            ax.set_ylabel('Improvement over SMLM in ' + r'$x$'+ '-direction')
            ax.grid()
            ax.set_xlim(left=-130, right = 130)
            ax.set_ylim(bottom = 0, top = 3.5)
            ax.tick_params(axis='both', which='major')
            ax.tick_params(axis='both', which='minor')                     

        if panel == 3:
            ax.plot((xrange-theta[0])/(10**-9),musum.T)
            ax.plot((xrange-theta[0])/(10**-9),musum_SMLM*np.ones(np.size(xrange)), color='k')
            ax.set_xlabel('Emitter-focus distance ' + r'$(x_f - \theta_x)$' + ' (nm)')
            ax.set_ylabel('Expected amount of photons')
            ax.grid()
            ax.set_xlim(left=-130, right = 130)
            ax.set_ylim(bottom = 0, top = 2900)
            ax.tick_params(axis='both', which='major')
            ax.tick_params(axis='both', which='minor')                  

        if panel == 4:
            ax.plot((xrange-theta[0])/(10**-9),signalsum.T)
            ax.plot((xrange-theta[0])/(10**-9),signalsum_SMLM*np.ones(np.size(xrange)), color='k')
            ax.set_xlabel('Emitter-focus distance ' + r'$(x_f - \theta_x)$' + ' (nm)')
            ax.set_ylabel('Illumination-intensity normalized' + '\n' + 'amount of signal photons')
            ax.grid()
            ax.set_xlim(left=-130, right = 130)
            ax.set_ylim(bottom = 0, top = 2100)
            ax.tick_params(axis='both', which='major')
            ax.tick_params(axis='both', which='minor')   

        if panel == 5:
            ax.plot((xrange-theta[0])/(10**-9),Bisum.T/100)
            ax.plot((xrange-theta[0])/(10**-9),Bisum_SMLM/100*np.ones(np.size(xrange)), color='k')
            ax.set_xlabel('Emitter-focus distance ' + r'$(x_f - \theta_x)$' + ' (nm)')
            ax.set_ylabel('Illumination-intensity normalized average' + '\n' + 'amount of background photons per pixel')
            ax.grid()
            ax.set_xlim(left=-130, right = 130)
            ax.set_ylim(bottom = 0, top = 9)
            ax.tick_params(axis='both', which='major')
            ax.tick_params(axis='both', which='minor')   

        if panel == 6:
            CRLBx, musum, signalsum, Bisum = load_4('thetaI_thetab_triangulation_rot_fb')
            thetab_list, thetaIrange = load_thetaI_thetab('thetaI_thetab_triangulation_rot_fb')
            
            ax.loglog(thetaIrange,CRLBx.T*10**9)
            ax.set_xlabel('Expected signal photon budget ' + r'$\theta_I$' + ' (photons)')
            ax.set_ylabel(r'$x$'+'-Cramér-Rao lower bound (nm)')
            ax.grid()
            ax.set_xlim(left=100, right = 10000)
            ax.set_ylim(bottom = 0.1, top = 40)
            ax.tick_params(axis='both', which='major')
            ax.tick_params(axis='both', which='minor')  
            
            #Legend
            C0_line = mlines.Line2D([0], [0],color='C0', label=(r'$\theta_b = 0$'))
            C1_line = mlines.Line2D([0], [0],color='C1', label=(r'$\theta_b = 1$'))
            C2_line = mlines.Line2D([0], [0],color='C2', label=(r'$\theta_b = 4$'))
            C3_line = mlines.Line2D([0], [0],color='C3', label=(r'$\theta_b = 8$'))
            C8_line = mlines.Line2D([0], [0],color='C4', label=(r'$\theta_b = 16$'))
            C9_line = mlines.Line2D([0], [0],color='None', label=('(in phot/px.)'))
    
            ax.legend(handles=[C0_line, C1_line, C2_line, C3_line, C8_line, C9_line], ncol = 2, prop={'size':13}, labelspacing = 1) 

        if panel == 0:
            #Legend
            C0_line = mlines.Line2D([0], [0],color='C0', label=r'$r = 0.5\sigma_{\mathrm{PSF}}$')
            C1_line = mlines.Line2D([0], [0],color='C1', label=r'$r = 1\sigma_{\mathrm{PSF}}$')
            C2_line = mlines.Line2D([0], [0],color='C2', label=r'$r = 1.5\sigma_{\mathrm{PSF}}$')
            C3_line = mlines.Line2D([0], [0],color='C3', label=r'$r = 2\sigma_{\mathrm{PSF}}$')
            C9_line = mlines.Line2D([0], [0],color='k', label=r'SMLM (Uniform' + '\n' + 'illumination intensity)')
            C10_line = mlines.Line2D([0], [0],color='k', linestyle='--', label=r'Localization on ISM' + '\n' + 'reconstruction data')
            C11_line = mlines.Line2D([0], [0],color='k', linestyle=':', label=r'Localization on' + '\n' + 'Fourier reweighted ISM' + '\n' + 'reconstruction data')
    
            ax.legend(handles=[C0_line, C1_line, C2_line, C3_line, C9_line, C10_line, C11_line], prop={'size':13}, loc='upper left', bbox_to_anchor=(1.04,1), labelspacing = 1)

        plt.savefig(('./Figures/S12/'+str(int(panel))) + '.svg')

    if plotnumber == 'S13':
        CRLBx, musum, signalsum, Bisum = load_4('xp_r_triangulation_fb')
        CRLBx_SMLM, musum_SMLM, signalsum_SMLM, Bisum_SMLM = load_4('SMLM')
        s_list, xrange = load_xp_r('xp_r_triangulation_fb')
        
        fig, ax = plt.subplots(1,1)
        fig.set_size_inches(5, 5, forward=True)
        
        if panel == 1:
            ax.semilogy((xrange-theta[0])/(10**-9),CRLBx.T*10**9)
            ax.semilogy((xrange-theta[0])/(10**-9),CRLBx_SMLM*10**9*np.ones(np.size(xrange)), color='k')
            ax.semilogy((xrange-theta[0])/(10**-9),CRLBx_ISM*10**9*np.ones(np.size(xrange)), color='k', linestyle='--')
            ax.semilogy((xrange-theta[0])/(10**-9),CRLBx_FR*10**9*np.ones(np.size(xrange)), color='k', linestyle=':')
            ax.set_xlabel('Emitter-focus distance ' + r'$(x_f - \theta_x)$' + ' (nm)')
            ax.set_ylabel(r'$x$'+'-Cramér-Rao lower bound (nm)')
            ax.grid()
            ax.set_xlim(left=-130, right = 130)
            ax.set_ylim(bottom = 0.6, top = 60)
            ax.tick_params(axis='both', which='major')
            ax.tick_params(axis='both', which='minor')
            
        if panel == 2:
            ax.plot((xrange-theta[0])/(10**-9),CRLBx_SMLM/CRLBx.T)
            ax.plot((xrange-theta[0])/(10**-9),imp_ISM*np.ones(np.size(xrange)), color='k', linestyle='--')
            ax.plot((xrange-theta[0])/(10**-9),imp_FR*np.ones(np.size(xrange)), color='k', linestyle=':')
            ax.set_xlabel('Emitter-focus distance ' + r'$(x_f - \theta_x)$' + ' (nm)')
            ax.set_ylabel('Improvement over SMLM in ' + r'$x$'+ '-direction')
            ax.grid()
            ax.set_xlim(left=-130, right = 130)
            ax.set_ylim(bottom = 0, top = 3.5)
            ax.tick_params(axis='both', which='major')
            ax.tick_params(axis='both', which='minor')                     

        if panel == 3:
            ax.plot((xrange-theta[0])/(10**-9),musum.T)
            ax.plot((xrange-theta[0])/(10**-9),musum_SMLM*np.ones(np.size(xrange)), color='k')
            ax.set_xlabel('Emitter-focus distance ' + r'$(x_f - \theta_x)$' + ' (nm)')
            ax.set_ylabel('Expected amount of photons')
            ax.grid()
            ax.set_xlim(left=-130, right = 130)
            ax.set_ylim(bottom = 0, top = 2900)
            ax.tick_params(axis='both', which='major')
            ax.tick_params(axis='both', which='minor')                  

        if panel == 4:
            ax.plot((xrange-theta[0])/(10**-9),signalsum.T)
            ax.plot((xrange-theta[0])/(10**-9),signalsum_SMLM*np.ones(np.size(xrange)), color='k')
            ax.set_xlabel('Emitter-focus distance ' + r'$(x_f - \theta_x)$' + ' (nm)')
            ax.set_ylabel('Illumination-intensity normalized' + '\n' + 'amount of signal photons')
            ax.grid()
            ax.set_xlim(left=-130, right = 130)
            ax.set_ylim(bottom = 0, top = 2100)
            ax.tick_params(axis='both', which='major')
            ax.tick_params(axis='both', which='minor')   

        if panel == 5:
            ax.plot((xrange-theta[0])/(10**-9),Bisum.T/100)
            ax.plot((xrange-theta[0])/(10**-9),Bisum_SMLM/100*np.ones(np.size(xrange)), color='k')
            ax.set_xlabel('Emitter-focus distance ' + r'$(x_f - \theta_x)$' + ' (nm)')
            ax.set_ylabel('Illumination-intensity normalized average' + '\n' + 'amount of background photons per pixel')
            ax.grid()
            ax.set_xlim(left=-130, right = 130)
            ax.set_ylim(bottom = 0, top = 9)
            ax.tick_params(axis='both', which='major')
            ax.tick_params(axis='both', which='minor')   

        if panel == 6:
            CRLBx, musum, signalsum, Bisum = load_4('thetaI_thetab_triangulation_fb')
            thetab_list, thetaIrange = load_thetaI_thetab('thetaI_thetab_triangulation_fb')
            
            ax.loglog(thetaIrange,CRLBx.T*10**9)
            ax.set_xlabel('Expected signal photon budget ' + r'$\theta_I$' + ' (photons)')
            ax.set_ylabel(r'$x$'+'-Cramér-Rao lower bound (nm)')
            ax.grid()
            ax.set_xlim(left=100, right = 10000)
            ax.set_ylim(bottom = 0.1, top = 40)
            ax.tick_params(axis='both', which='major')
            ax.tick_params(axis='both', which='minor')  
            
            #Legend
            C0_line = mlines.Line2D([0], [0],color='C0', label=(r'$\theta_b = 0$'))
            C1_line = mlines.Line2D([0], [0],color='C1', label=(r'$\theta_b = 1$'))
            C2_line = mlines.Line2D([0], [0],color='C2', label=(r'$\theta_b = 4$'))
            C3_line = mlines.Line2D([0], [0],color='C3', label=(r'$\theta_b = 8$'))
            C8_line = mlines.Line2D([0], [0],color='C4', label=(r'$\theta_b = 16$'))
            C9_line = mlines.Line2D([0], [0],color='None', label=('(in phot/px.)'))
    
            ax.legend(handles=[C0_line, C1_line, C2_line, C3_line, C8_line, C9_line], ncol = 2, prop={'size':13}, labelspacing = 1)    
            
        if panel == 0:
            #Legend
            C0_line = mlines.Line2D([0], [0],color='C0', label=r'$r = 0.5\sigma_{\mathrm{PSF}}$')
            C1_line = mlines.Line2D([0], [0],color='C1', label=r'$r = 1\sigma_{\mathrm{PSF}}$')
            C2_line = mlines.Line2D([0], [0],color='C2', label=r'$r = 1.5\sigma_{\mathrm{PSF}}$')
            C3_line = mlines.Line2D([0], [0],color='C3', label=r'$r = 2\sigma_{\mathrm{PSF}}$')
            C9_line = mlines.Line2D([0], [0],color='k', label=r'SMLM (Uniform' + '\n' + 'illumination intensity)')
            C10_line = mlines.Line2D([0], [0],color='k', linestyle='--', label=r'Localization on ISM' + '\n' + 'reconstruction data')
            C11_line = mlines.Line2D([0], [0],color='k', linestyle=':', label=r'Localization on' + '\n' + 'Fourier reweighted ISM' + '\n' + 'reconstruction data')
    
            ax.legend(handles=[C0_line, C1_line, C2_line, C3_line, C9_line, C10_line, C11_line], prop={'size':13}, loc='upper left', bbox_to_anchor=(1.04,1), labelspacing = 1)

        plt.savefig(('./Figures/S13/'+str(int(panel))) + '.svg')

    if plotnumber == 'S14':
        CRLBx, musum, signalsum, Bisum = load_4('xp_r_triangulation_rot_nocenter_doughnut_fb')
        CRLBx_SMLM, musum_SMLM, signalsum_SMLM, Bisum_SMLM = load_4('SMLM')
        s_list, xrange = load_xp_r('xp_r_triangulation_rot_nocenter_doughnut_fb')
        
        fig, ax = plt.subplots(1,1)
        fig.set_size_inches(5, 5, forward=True)
        
        if panel == 1:
            ax.semilogy((xrange-theta[0])/(10**-9),CRLBx.T[:,0:4]*10**9)
            ax.semilogy((xrange-theta[0])/(10**-9),CRLBx_SMLM*10**9*np.ones(np.size(xrange)), color='k')
            ax.semilogy((xrange-theta[0])/(10**-9),CRLBx_ISM*10**9*np.ones(np.size(xrange)), color='k', linestyle='--')
            ax.semilogy((xrange-theta[0])/(10**-9),CRLBx_FR*10**9*np.ones(np.size(xrange)), color='k', linestyle=':')
            ax.set_xlabel('Emitter-focus distance ' + r'$(x_f - \theta_x)$' + ' (nm)')
            ax.set_ylabel(r'$x$'+'-Cramér-Rao lower bound (nm)')
            ax.grid()
            ax.set_xlim(left=-130, right = 130)
            ax.set_ylim(bottom = 0.6, top = 60)
            ax.tick_params(axis='both', which='major')
            ax.tick_params(axis='both', which='minor')
            
        if panel == 2:
            ax.plot((xrange-theta[0])/(10**-9),CRLBx_SMLM/CRLBx.T[:,0:4])
            ax.plot((xrange-theta[0])/(10**-9),imp_ISM*np.ones(np.size(xrange)), color='k', linestyle='--')
            ax.plot((xrange-theta[0])/(10**-9),imp_FR*np.ones(np.size(xrange)), color='k', linestyle=':')
            ax.set_xlabel('Emitter-focus distance ' + r'$(x_f - \theta_x)$' + ' (nm)')
            ax.set_ylabel('Improvement over SMLM in ' + r'$x$'+ '-direction')
            ax.grid()
            ax.set_xlim(left=-130, right = 130)
            ax.set_ylim(bottom = 0, top = 3.5)
            ax.tick_params(axis='both', which='major')
            ax.tick_params(axis='both', which='minor')                     

        if panel == 3:
            ax.plot((xrange-theta[0])/(10**-9),musum.T[:,0:4])
            ax.plot((xrange-theta[0])/(10**-9),musum_SMLM*np.ones(np.size(xrange)), color='k')
            ax.set_xlabel('Emitter-focus distance ' + r'$(x_f - \theta_x)$' + ' (nm)')
            ax.set_ylabel('Expected amount of photons')
            ax.grid()
            ax.set_xlim(left=-130, right = 130)
            ax.set_ylim(bottom = 0, top = 7100)
            ax.tick_params(axis='both', which='major')
            ax.tick_params(axis='both', which='minor')                  

        if panel == 4:
            ax.plot((xrange-theta[0])/(10**-9),signalsum.T[:,0:4])
            ax.plot((xrange-theta[0])/(10**-9),signalsum_SMLM*np.ones(np.size(xrange)), color='k')
            ax.set_xlabel('Emitter-focus distance ' + r'$(x_f - \theta_x)$' + ' (nm)')
            ax.set_ylabel('Illumination-intensity normalized' + '\n' + 'amount of signal photons')
            ax.grid()
            ax.set_xlim(left=-130, right = 130)
            ax.set_ylim(bottom = 0, top = 2100)
            ax.tick_params(axis='both', which='major')
            ax.tick_params(axis='both', which='minor')   

        if panel == 5:
            ax.plot((xrange-theta[0])/(10**-9),Bisum.T[:,0:4]/100)
            ax.plot((xrange-theta[0])/(10**-9),Bisum_SMLM/100*np.ones(np.size(xrange)), color='k')
            ax.set_xlabel('Emitter-focus distance ' + r'$(x_f - \theta_x)$' + ' (nm)')
            ax.set_ylabel('Illumination-intensity normalized average' + '\n' + 'amount of background photons per pixel')
            ax.grid()
            ax.set_xlim(left=-130, right = 130)
            ax.set_ylim(bottom = 0, top = 500)
            ax.tick_params(axis='both', which='major')
            ax.tick_params(axis='both', which='minor')

        if panel == 6:
            CRLBx, musum, signalsum, Bisum = load_4('thetaI_thetab_triangulation_rot_nocenter_doughnut_fb')
            thetab_list, thetaIrange = load_thetaI_thetab('thetaI_thetab_triangulation_rot_nocenter_doughnut_fb')
            
            ax.loglog(thetaIrange,CRLBx.T*10**9)
            ax.set_xlabel('Expected signal photon budget ' + r'$\theta_I$' + ' (photons)')
            ax.set_ylabel(r'$x$'+'-Cramér-Rao lower bound (nm)')
            ax.grid()
            ax.set_xlim(left=100, right = 10000)
            ax.set_ylim(bottom = 0.1, top = 40)
            ax.tick_params(axis='both', which='major')
            ax.tick_params(axis='both', which='minor')  
            
            #Legend
            C0_line = mlines.Line2D([0], [0],color='C0', label=(r'$\theta_b = 0$'))
            C1_line = mlines.Line2D([0], [0],color='C1', label=(r'$\theta_b = 1$'))
            C2_line = mlines.Line2D([0], [0],color='C2', label=(r'$\theta_b = 4$'))
            C3_line = mlines.Line2D([0], [0],color='C3', label=(r'$\theta_b = 8$'))
            C8_line = mlines.Line2D([0], [0],color='C4', label=(r'$\theta_b = 16$'))
            C9_line = mlines.Line2D([0], [0],color='None', label=('(in phot/px.)'))
    
            ax.legend(handles=[C0_line, C1_line, C2_line, C3_line, C8_line, C9_line], ncol = 2, prop={'size':13}, labelspacing = 1)  

        if panel == 0:
            #Legend
            C0_line = mlines.Line2D([0], [0],color='C0', label=r'$r = 1\sigma_{\mathrm{PSF}}$')
            C1_line = mlines.Line2D([0], [0],color='C1', label=r'$r = 2\sigma_{\mathrm{PSF}}$')
            C2_line = mlines.Line2D([0], [0],color='C2', label=r'$r = 3\sigma_{\mathrm{PSF}}$')
            C3_line = mlines.Line2D([0], [0],color='C3', label=r'$r = 4\sigma_{\mathrm{PSF}}$')
            C9_line = mlines.Line2D([0], [0],color='k', label=r'SMLM (Uniform' + '\n' + 'illumination intensity)')
            C10_line = mlines.Line2D([0], [0],color='k', linestyle='--', label=r'Localization on ISM' + '\n' + 'reconstruction data')
            C11_line = mlines.Line2D([0], [0],color='k', linestyle=':', label=r'Localization on' + '\n' + 'Fourier reweighted ISM' + '\n' + 'reconstruction data')
    
            ax.legend(handles=[C0_line, C1_line, C2_line, C3_line, C9_line, C10_line, C11_line], prop={'size':13}, loc='upper left', bbox_to_anchor=(1.04,1), labelspacing = 1)

        plt.savefig(('./Figures/S14/'+str(int(panel))) + '.svg')
        
    if plotnumber == 'S15':
        CRLBx, musum, signalsum, Bisum = load_4('xp_r_triangulation_nocenter_doughnut_fb')
        CRLBx_SMLM, musum_SMLM, signalsum_SMLM, Bisum_SMLM = load_4('SMLM')
        s_list, xrange = load_xp_r('xp_r_triangulation_nocenter_doughnut_fb')
        
        fig, ax = plt.subplots(1,1)
        fig.set_size_inches(5, 5, forward=True)
        
        if panel == 1:
            ax.semilogy((xrange-theta[0])/(10**-9),CRLBx.T[:,0:4]*10**9)
            ax.semilogy((xrange-theta[0])/(10**-9),CRLBx_SMLM*10**9*np.ones(np.size(xrange)), color='k')
            ax.semilogy((xrange-theta[0])/(10**-9),CRLBx_ISM*10**9*np.ones(np.size(xrange)), color='k', linestyle='--')
            ax.semilogy((xrange-theta[0])/(10**-9),CRLBx_FR*10**9*np.ones(np.size(xrange)), color='k', linestyle=':')
            ax.set_xlabel('Emitter-focus distance ' + r'$(x_f - \theta_x)$' + ' (nm)')
            ax.set_ylabel(r'$x$'+'-Cramér-Rao lower bound (nm)')
            ax.grid()
            ax.set_xlim(left=-130, right = 130)
            ax.set_ylim(bottom = 0.6, top = 100)
            ax.tick_params(axis='both', which='major')
            ax.tick_params(axis='both', which='minor')
            
        if panel == 2:
            ax.plot((xrange-theta[0])/(10**-9),CRLBx_SMLM/CRLBx.T[:,0:4])
            ax.plot((xrange-theta[0])/(10**-9),imp_ISM*np.ones(np.size(xrange)), color='k', linestyle='--')
            ax.plot((xrange-theta[0])/(10**-9),imp_FR*np.ones(np.size(xrange)), color='k', linestyle=':')
            ax.set_xlabel('Emitter-focus distance ' + r'$(x_f - \theta_x)$' + ' (nm)')
            ax.set_ylabel('Improvement over SMLM in ' + r'$x$'+ '-direction')
            ax.grid()
            ax.set_xlim(left=-130, right = 130)
            ax.set_ylim(bottom = 0, top = 3.5)
            ax.tick_params(axis='both', which='major')
            ax.tick_params(axis='both', which='minor')                     

        if panel == 3:
            ax.plot((xrange-theta[0])/(10**-9),musum.T[:,0:4])
            ax.plot((xrange-theta[0])/(10**-9),musum_SMLM*np.ones(np.size(xrange)), color='k')
            ax.set_xlabel('Emitter-focus distance ' + r'$(x_f - \theta_x)$' + ' (nm)')
            ax.set_ylabel('Expected amount of photons')
            ax.grid()
            ax.set_xlim(left=-130, right = 130)
            ax.tick_params(axis='both', which='major')
            ax.tick_params(axis='both', which='minor')                  

        if panel == 4:
            ax.plot((xrange-theta[0])/(10**-9),signalsum.T[:,0:4])
            ax.plot((xrange-theta[0])/(10**-9),signalsum_SMLM*np.ones(np.size(xrange)), color='k')
            ax.set_xlabel('Emitter-focus distance ' + r'$(x_f - \theta_x)$' + ' (nm)')
            ax.set_ylabel('Illumination-intensity normalized' + '\n' + 'amount of signal photons')
            ax.grid()
            ax.set_xlim(left=-130, right = 130)
            ax.set_ylim(bottom = 0, top = 2100)
            ax.tick_params(axis='both', which='major')
            ax.tick_params(axis='both', which='minor')   

        if panel == 5:
            ax.plot((xrange-theta[0])/(10**-9),Bisum.T[:,0:4]/100)
            ax.plot((xrange-theta[0])/(10**-9),Bisum_SMLM/100*np.ones(np.size(xrange)), color='k')
            ax.set_xlabel('Emitter-focus distance ' + r'$(x_f - \theta_x)$' + ' (nm)')
            ax.set_ylabel('Illumination-intensity normalized average' + '\n' + 'amount of background photons per pixel')
            ax.grid()
            ax.set_xlim(left=-130, right = 130)
            ax.set_ylim(bottom = 0, top = 500)
            ax.tick_params(axis='both', which='major')
            ax.tick_params(axis='both', which='minor')  

        if panel == 6:
            CRLBx, musum, signalsum, Bisum = load_4('thetaI_thetab_triangulation_nocenter_doughnut_fb')
            thetab_list, thetaIrange = load_thetaI_thetab('thetaI_thetab_triangulation_nocenter_doughnut_fb')
            
            ax.loglog(thetaIrange,CRLBx.T*10**9)
            ax.set_xlabel('Expected signal photon budget ' + r'$\theta_I$' + ' (photons)')
            ax.set_ylabel(r'$x$'+'-Cramér-Rao lower bound (nm)')
            ax.grid()
            ax.set_xlim(left=100, right = 10000)
            ax.set_ylim(bottom = 0.1, top = 40)
            ax.tick_params(axis='both', which='major')
            ax.tick_params(axis='both', which='minor')  
            
            #Legend
            C0_line = mlines.Line2D([0], [0],color='C0', label=(r'$\theta_b = 0$'))
            C1_line = mlines.Line2D([0], [0],color='C1', label=(r'$\theta_b = 1$'))
            C2_line = mlines.Line2D([0], [0],color='C2', label=(r'$\theta_b = 4$'))
            C3_line = mlines.Line2D([0], [0],color='C3', label=(r'$\theta_b = 8$'))
            C8_line = mlines.Line2D([0], [0],color='C4', label=(r'$\theta_b = 16$'))
            C9_line = mlines.Line2D([0], [0],color='None', label=('(in phot/px.)'))
    
            ax.legend(handles=[C0_line, C1_line, C2_line, C3_line, C8_line, C9_line], ncol = 2, prop={'size':13}, labelspacing = 1)    

        if panel == 0:
            #Legend
            C0_line = mlines.Line2D([0], [0],color='C0', label=r'$r = 1\sigma_{\mathrm{PSF}}$')
            C1_line = mlines.Line2D([0], [0],color='C1', label=r'$r = 2\sigma_{\mathrm{PSF}}$')
            C2_line = mlines.Line2D([0], [0],color='C2', label=r'$r = 3\sigma_{\mathrm{PSF}}$')
            C3_line = mlines.Line2D([0], [0],color='C3', label=r'$r = 4\sigma_{\mathrm{PSF}}$')
            C9_line = mlines.Line2D([0], [0],color='k', label=r'SMLM (Uniform' + '\n' + 'illumination intensity)')
            C10_line = mlines.Line2D([0], [0],color='k', linestyle='--', label=r'Localization on ISM' + '\n' + 'reconstruction data')
            C11_line = mlines.Line2D([0], [0],color='k', linestyle=':', label=r'Localization on' + '\n' + 'Fourier reweighted ISM' + '\n' + 'reconstruction data')
    
            ax.legend(handles=[C0_line, C1_line, C2_line, C3_line, C9_line, C10_line, C11_line], prop={'size':13}, loc='upper left', bbox_to_anchor=(1.04,1), labelspacing = 1)

        plt.savefig(('./Figures/S15/'+str(int(panel))) + '.svg')

    if plotnumber == 'S17':
        CRLBx, musum, signalsum, Bisum = load_4('xp_r_triangulation_doughnut_ext_fb')
        CRLBx_SMLM, musum_SMLM, signalsum_SMLM, Bisum_SMLM = load_4('SMLM')
        s_list, xrange = load_xp_r('xp_r_triangulation_doughnut_ext_fb')
        
        fig, ax = plt.subplots(1,1)
        fig.set_size_inches(5, 5, forward=True)
        
        if panel == 1:
            ax.semilogy((xrange-theta[0])/(10**-9),CRLBx.T*10**9)
            ax.semilogy((xrange-theta[0])/(10**-9),CRLBx_SMLM*10**9*np.ones(np.size(xrange)), color='k')
            ax.semilogy((xrange-theta[0])/(10**-9),CRLBx_ISM*10**9*np.ones(np.size(xrange)), color='k', linestyle='--')
            ax.semilogy((xrange-theta[0])/(10**-9),CRLBx_FR*10**9*np.ones(np.size(xrange)), color='k', linestyle=':')
            ax.set_xlabel('Emitter-focus distance ' + r'$(x_f - \theta_x)$' + ' (nm)')
            ax.set_ylabel(r'$x$'+'-Cramér-Rao lower bound (nm)')
            ax.grid()
            ax.set_xlim(left=-130, right = 130)
            ax.set_ylim(bottom = 0.4, top = 60)
            ax.tick_params(axis='both', which='major')
            ax.tick_params(axis='both', which='minor')
            
        if panel == 2:
            ax.plot((xrange-theta[0])/(10**-9),CRLBx_SMLM/CRLBx.T)
            ax.plot((xrange-theta[0])/(10**-9),imp_ISM*np.ones(np.size(xrange)), color='k', linestyle='--')
            ax.plot((xrange-theta[0])/(10**-9),imp_FR*np.ones(np.size(xrange)), color='k', linestyle=':')
            ax.set_xlabel('Emitter-focus distance ' + r'$(x_f - \theta_x)$' + ' (nm)')
            ax.set_ylabel('Improvement over SMLM in ' + r'$x$'+ '-direction')
            ax.grid()
            ax.set_xlim(left=-130, right = 130)
            ax.set_ylim(bottom = 0, top = 4)
            ax.tick_params(axis='both', which='major')
            ax.tick_params(axis='both', which='minor')                     

        if panel == 3:
            ax.plot((xrange-theta[0])/(10**-9),musum.T)
            ax.plot((xrange-theta[0])/(10**-9),musum_SMLM*np.ones(np.size(xrange)), color='k')
            ax.set_xlabel('Emitter-focus distance ' + r'$(x_f - \theta_x)$' + ' (nm)')
            ax.set_ylabel('Expected amount of photons')
            ax.grid()
            ax.set_xlim(left=-130, right = 130)
            ax.set_ylim(bottom = 0, top = 2900)
            ax.tick_params(axis='both', which='major')
            ax.tick_params(axis='both', which='minor')                  

        if panel == 4:
            ax.plot((xrange-theta[0])/(10**-9),signalsum.T)
            ax.plot((xrange-theta[0])/(10**-9),signalsum_SMLM*np.ones(np.size(xrange)), color='k')
            ax.set_xlabel('Emitter-focus distance ' + r'$(x_f - \theta_x)$' + ' (nm)')
            ax.set_ylabel('Illumination-intensity normalized' + '\n' + 'amount of signal photons')
            ax.grid()
            ax.set_xlim(left=-130, right = 130)
            ax.set_ylim(bottom = 0, top = 2100)
            ax.tick_params(axis='both', which='major')
            ax.tick_params(axis='both', which='minor')   

        if panel == 5:
            ax.plot((xrange-theta[0])/(10**-9),Bisum.T/100)
            ax.plot((xrange-theta[0])/(10**-9),Bisum_SMLM/100*np.ones(np.size(xrange)), color='k')
            ax.set_xlabel('Emitter-focus distance ' + r'$(x_f - \theta_x)$' + ' (nm)')
            ax.set_ylabel('Illumination-intensity normalized average' + '\n' + 'amount of background photons per pixel')
            ax.grid()
            ax.set_xlim(left=-130, right = 130)
            ax.set_ylim(bottom = 0, top = 700)
            ax.tick_params(axis='both', which='major')
            ax.tick_params(axis='both', which='minor')   

        if panel == 6:
            CRLBx, musum, signalsum, Bisum = load_4('thetaI_thetab_triangulation_doughnut_fb')
            thetab_list, thetaIrange = load_thetaI_thetab('thetaI_thetab_triangulation_doughnut_fb')
            
            ax.loglog(thetaIrange,CRLBx.T*10**9)
            ax.set_xlabel('Expected signal photon budget ' + r'$\theta_I$' + ' (photons)')
            ax.set_ylabel(r'$x$'+'-Cramér-Rao lower bound (nm)')
            ax.grid()
            ax.set_xlim(left=100, right = 10000)
            ax.set_ylim(bottom = 0.1, top = 40)
            ax.tick_params(axis='both', which='major')
            ax.tick_params(axis='both', which='minor')  
            
            #Legend
            C0_line = mlines.Line2D([0], [0],color='C0', label=(r'$\theta_b = 0$'))
            C1_line = mlines.Line2D([0], [0],color='C1', label=(r'$\theta_b = 1$'))
            C2_line = mlines.Line2D([0], [0],color='C2', label=(r'$\theta_b = 4$'))
            C3_line = mlines.Line2D([0], [0],color='C3', label=(r'$\theta_b = 8$'))
            C8_line = mlines.Line2D([0], [0],color='C4', label=(r'$\theta_b = 16$'))
            C9_line = mlines.Line2D([0], [0],color='None', label=('(in phot/px.)'))
    
            ax.legend(handles=[C0_line, C1_line, C2_line, C3_line, C8_line, C9_line], ncol = 2, prop={'size':13}, labelspacing = 1)    
            
        if panel == 0:
            #Legend
            C0_line = mlines.Line2D([0], [0],color='C0', label=r'$r = 0.5\sigma_{\mathrm{PSF}}$')
            C1_line = mlines.Line2D([0], [0],color='C1', label=r'$r = 1\sigma_{\mathrm{PSF}}$')
            C2_line = mlines.Line2D([0], [0],color='C2', label=r'$r = 2\sigma_{\mathrm{PSF}}$')
            C3_line = mlines.Line2D([0], [0],color='C3', label=r'$r = 3\sigma_{\mathrm{PSF}}$')
            C4_line = mlines.Line2D([0], [0],color='C4', label=r'$r = 4\sigma_{\mathrm{PSF}}$')
            C9_line = mlines.Line2D([0], [0],color='k', label=r'SMLM (Uniform' + '\n' + 'illumination intensity)')
            C10_line = mlines.Line2D([0], [0],color='k', linestyle='--', label=r'Localization on ISM' + '\n' + 'reconstruction data')
            C11_line = mlines.Line2D([0], [0],color='k', linestyle=':', label=r'Localization on' + '\n' + 'Fourier reweighted ISM' + '\n' + 'reconstruction data')
    
            ax.legend(handles=[C0_line, C1_line, C2_line, C3_line, C4_line, C9_line, C10_line, C11_line], prop={'size':13}, loc='upper left', bbox_to_anchor=(1.04,1), labelspacing = 1)

        plt.savefig(('./Figures/S17/'+str(int(panel))) + '.svg')

    if plotnumber == 'S18':
        CRLBx, musum, signalsum, Bisum = load_4('xp_rp_1pattern')
        CRLBx_SMLM, musum_SMLM, signalsum_SMLM, Bisum_SMLM = load_4('SMLM')
        r_p_list, xrange = load_xp_rp('xp_rp_1pattern')
        
        fig, ax = plt.subplots(1,1)
        fig.set_size_inches(5, 5, forward=True)
        
        if panel == 1:
            ax.semilogy((xrange-theta[0])/(10**-9),CRLBx.T[:,0:4]*10**9)
            ax.semilogy((xrange-theta[0])/(10**-9),CRLBx.T[:,-1]*10**9)
            ax.semilogy((xrange-theta[0])/(10**-9),CRLBx_SMLM*10**9*np.ones(np.size(xrange)), color='k')
            ax.semilogy((xrange-theta[0])/(10**-9),CRLBx_ISM*10**9*np.ones(np.size(xrange)), color='k', linestyle='--')
            ax.semilogy((xrange-theta[0])/(10**-9),CRLBx_FR*10**9*np.ones(np.size(xrange)), color='k', linestyle=':')
            ax.set_xlabel('Emitter-pinhole distance ' + r'$(x_p - \theta_x)$' + ' (nm)')
            ax.set_ylabel(r'$x$'+'-Cramér-Rao lower bound (nm)')
            ax.grid()
            ax.set_xlim(left=-130, right = 130)
            ax.set_ylim(bottom = 0.6, top = 60)
            ax.tick_params(axis='both', which='major')
            ax.tick_params(axis='both', which='minor')
            
        if panel == 2:
            ax.plot((xrange-theta[0])/(10**-9),CRLBx_SMLM/CRLBx.T[:,0:4])
            ax.plot((xrange-theta[0])/(10**-9),CRLBx_SMLM/CRLBx.T[:,-1])
            ax.plot((xrange-theta[0])/(10**-9),imp_ISM*np.ones(np.size(xrange)), color='k', linestyle='--')
            ax.plot((xrange-theta[0])/(10**-9),imp_FR*np.ones(np.size(xrange)), color='k', linestyle=':')
            ax.set_xlabel('Emitter-pinhole distance ' + r'$(x_p - \theta_x)$' + ' (nm)')
            ax.set_ylabel('Improvement over SMLM in ' + r'$x$'+ '-direction')
            ax.grid()
            ax.set_xlim(left=-130, right = 130)
            ax.set_ylim(bottom = 0, top = 2.2)
            ax.tick_params(axis='both', which='major')
            ax.tick_params(axis='both', which='minor')                     

        if panel == 3:
            ax.plot((xrange-theta[0])/(10**-9),musum.T[:,0:4])
            ax.plot((xrange-theta[0])/(10**-9),musum.T[:,-1])
            ax.plot((xrange-theta[0])/(10**-9),musum_SMLM*np.ones(np.size(xrange)), color='k')
            ax.set_xlabel('Emitter-pinhole distance ' + r'$(x_p - \theta_x)$' + ' (nm)')
            ax.set_ylabel('Expected amount of photons')
            ax.grid()
            ax.set_xlim(left=-130, right = 130)
            ax.set_ylim(bottom = 0, top = 2900)
            ax.tick_params(axis='both', which='major')
            ax.tick_params(axis='both', which='minor')                  

        if panel == 4:
            ax.plot((xrange-theta[0])/(10**-9),signalsum.T[:,0:4])
            ax.plot((xrange-theta[0])/(10**-9),signalsum.T[:,-1])
            ax.plot((xrange-theta[0])/(10**-9),signalsum_SMLM*np.ones(np.size(xrange)), color='k')
            ax.set_xlabel('Emitter-pinhole distance ' + r'$(x_p - \theta_x)$' + ' (nm)')
            ax.set_ylabel('Illumination-intensity normalized' + '\n' + 'amount of signal photons')
            ax.grid()
            ax.set_xlim(left=-130, right = 130)
            ax.set_ylim(bottom = 0, top = 2100)
            ax.tick_params(axis='both', which='major')
            ax.tick_params(axis='both', which='minor')   

        if panel == 5:
            ax.plot((xrange-theta[0])/(10**-9),Bisum.T[:,0:4]/100)
            ax.plot((xrange-theta[0])/(10**-9),Bisum.T[:,-1]/100)
            ax.plot((xrange-theta[0])/(10**-9),Bisum_SMLM/100*np.ones(np.size(xrange)), color='k')
            ax.set_xlabel('Emitter-pinhole distance ' + r'$(x_p - \theta_x)$' + ' (nm)')
            ax.set_ylabel('Illumination-intensity normalized average' + '\n' + 'amount of background photons per pixel')
            ax.grid()
            ax.set_xlim(left=-130, right = 130)
            ax.set_ylim(bottom = 0, top = 9)
            ax.tick_params(axis='both', which='major')
            ax.tick_params(axis='both', which='minor')   
            
        if panel == 6:
            CRLBx, musum, signalsum, Bisum = load_4('thetaI_thetab_1pattern')
            thetab_list, thetaIrange = load_thetaI_thetab('thetaI_thetab_1pattern')
            
            ax.loglog(thetaIrange,CRLBx.T*10**9)
            ax.set_xlabel('Expected signal photon budget ' + r'$\theta_I$' + ' (photons)')
            ax.set_ylabel(r'$x$'+'-Cramér-Rao lower bound (nm)')
            ax.grid()
            ax.set_xlim(left=100, right = 10000)
            ax.set_ylim(bottom = 0.6, top = 60)
            ax.tick_params(axis='both', which='major')
            ax.tick_params(axis='both', which='minor') 

            #Legend
            C0_line = mlines.Line2D([0], [0],color='C0', label=(r'$\theta_b = 0$'))
            C1_line = mlines.Line2D([0], [0],color='C1', label=(r'$\theta_b = 1$'))
            C2_line = mlines.Line2D([0], [0],color='C2', label=(r'$\theta_b = 4$'))
            C3_line = mlines.Line2D([0], [0],color='C3', label=(r'$\theta_b = 8$'))
            C8_line = mlines.Line2D([0], [0],color='C4', label=(r'$\theta_b = 16$'))
            C9_line = mlines.Line2D([0], [0],color='None', label=('(in phot/px.)'))
    
            ax.legend(handles=[C0_line, C1_line, C2_line, C3_line, C8_line, C9_line], ncol = 2, prop={'size':13}, labelspacing = 1)          

        if panel == 0:
            #Legend
            C0_line = mlines.Line2D([0], [0],color='C0', label=r'$r_p = 1\sigma_{\mathrm{PSF}}$')
            C1_line = mlines.Line2D([0], [0],color='C1', label=r'$r_p = 2\sigma_{\mathrm{PSF}}$')
            C2_line = mlines.Line2D([0], [0],color='C2', label=r'$r_p = 3\sigma_{\mathrm{PSF}}$')
            C3_line = mlines.Line2D([0], [0],color='C3', label=r'$r_p = 4\sigma_{\mathrm{PSF}}$')
            C8_line = mlines.Line2D([0], [0],color='C4', label=r'$r_p = \infty$')
            C9_line = mlines.Line2D([0], [0],color='k', label=r'SMLM (Uniform' + '\n' + 'illumination intensity)')
            C10_line = mlines.Line2D([0], [0],color='k', linestyle='--', label=r'Localization on ISM' + '\n' + 'reconstruction data')
            C11_line = mlines.Line2D([0], [0],color='k', linestyle=':', label=r'Localization on' + '\n' + 'Fourier reweighted ISM' + '\n' + 'reconstruction data')
    
            ax.legend(handles=[C0_line, C1_line, C2_line, C3_line, C8_line, C9_line, C10_line, C11_line], prop={'size':13}, loc='upper left', bbox_to_anchor=(1.04,1), labelspacing = 1)
        
        plt.savefig(('./Figures/S18/'+str(int(panel))) + '.svg')

    if plotnumber == 'S19':
        CRLBx, musum, signalsum, Bisum = load_4('yp_rp_1pattern')
        CRLBx_SMLM, musum_SMLM, signalsum_SMLM, Bisum_SMLM = load_4('SMLM')
        r_p_list, yrange = load_yp_rp('yp_rp_1pattern')
        
        fig, ax = plt.subplots(1,1)
        fig.set_size_inches(5, 5, forward=True)
        
        if panel == 1:
            ax.semilogy((yrange-theta[1])/(10**-9),CRLBx.T[:,0:4]*10**9)
            ax.semilogy((yrange-theta[1])/(10**-9),CRLBx.T[:,-1]*10**9)
            ax.semilogy((yrange-theta[1])/(10**-9),CRLBx_SMLM*10**9*np.ones(np.size(yrange)), color='k')
            ax.semilogy((yrange-theta[1])/(10**-9),CRLBx_ISM*10**9*np.ones(np.size(yrange)), color='k', linestyle='--')
            ax.semilogy((yrange-theta[1])/(10**-9),CRLBx_FR*10**9*np.ones(np.size(yrange)), color='k', linestyle=':')
            ax.set_xlabel('Emitter-pinhole distance ' + r'$(y_p - \theta_y)$' + ' (nm)')
            ax.set_ylabel(r'$x$'+'-Cramér-Rao lower bound (nm)')
            ax.grid()
            ax.set_xlim(left=-130, right = 130)
            ax.set_ylim(bottom = 0.6, top = 60)
            ax.tick_params(axis='both', which='major')
            ax.tick_params(axis='both', which='minor')
            
        if panel == 2:
            ax.plot((yrange-theta[1])/(10**-9),CRLBx_SMLM/CRLBx.T[:,0:4])
            ax.plot((yrange-theta[1])/(10**-9),CRLBx_SMLM/CRLBx.T[:,-1])
            ax.plot((yrange-theta[1])/(10**-9),imp_ISM*np.ones(np.size(yrange)), color='k', linestyle='--')
            ax.plot((yrange-theta[1])/(10**-9),imp_FR*np.ones(np.size(yrange)), color='k', linestyle=':')
            ax.set_xlabel('Emitter-pinhole distance ' + r'$(y_p - \theta_y)$' + ' (nm)')
            ax.set_ylabel('Improvement over SMLM in ' + r'$x$'+ '-direction')
            ax.grid()
            ax.set_xlim(left=-130, right = 130)
            ax.set_ylim(bottom = 0, top = 2.2)
            ax.tick_params(axis='both', which='major')
            ax.tick_params(axis='both', which='minor')                     

        if panel == 3:
            ax.plot((yrange-theta[1])/(10**-9),musum.T[:,0:4])
            ax.plot((yrange-theta[1])/(10**-9),musum.T[:,-1])
            ax.plot((yrange-theta[1])/(10**-9),musum_SMLM*np.ones(np.size(yrange)), color='k')
            ax.set_xlabel('Emitter-pinhole distance ' + r'$(y_p - \theta_y)$' + ' (nm)')
            ax.set_ylabel('Expected amount of photons')
            ax.grid()
            ax.set_xlim(left=-130, right = 130)
            ax.set_ylim(bottom = 0, top = 2900)
            ax.tick_params(axis='both', which='major')
            ax.tick_params(axis='both', which='minor')                  

        if panel == 4:
            ax.plot((yrange-theta[1])/(10**-9),signalsum.T[:,0:4])
            ax.plot((yrange-theta[1])/(10**-9),signalsum.T[:,-1])
            ax.plot((yrange-theta[1])/(10**-9),signalsum_SMLM*np.ones(np.size(yrange)), color='k')
            ax.set_xlabel('Emitter-pinhole distance ' + r'$(y_p - \theta_y)$' + ' (nm)')
            ax.set_ylabel('Illumination-intensity normalized' + '\n' + 'amount of signal photons')
            ax.grid()
            ax.set_xlim(left=-130, right = 130)
            ax.set_ylim(bottom = 0, top = 2100)
            ax.tick_params(axis='both', which='major')
            ax.tick_params(axis='both', which='minor')   

        if panel == 5:
            ax.plot((yrange-theta[1])/(10**-9),Bisum.T[:,0:4]/100)
            ax.plot((yrange-theta[1])/(10**-9),Bisum.T[:,-1]/100)
            ax.plot((yrange-theta[1])/(10**-9),Bisum_SMLM*np.ones(np.size(yrange)), color='k')
            ax.set_xlabel('Emitter-pinhole distance ' + r'$(y_p - \theta_y)$' + ' (nm)')
            ax.set_ylabel('Illumination-intensity normalized average' + '\n' + 'amount of background photons per pixel')
            ax.grid()
            ax.set_xlim(left=-130, right = 130)
            ax.set_ylim(bottom = 0, top = 9)
            ax.tick_params(axis='both', which='major')
            ax.tick_params(axis='both', which='minor')  
        
        if panel == 0:
            #Legend
            C0_line = mlines.Line2D([0], [0],color='C0', label=r'$r_p = 1\sigma_{\mathrm{PSF}}$')
            C1_line = mlines.Line2D([0], [0],color='C1', label=r'$r_p = 2\sigma_{\mathrm{PSF}}$')
            C2_line = mlines.Line2D([0], [0],color='C2', label=r'$r_p = 3\sigma_{\mathrm{PSF}}$')
            C3_line = mlines.Line2D([0], [0],color='C3', label=r'$r_p = 4\sigma_{\mathrm{PSF}}$')
            C8_line = mlines.Line2D([0], [0],color='C4', label=r'$r_p = \infty$')
            C9_line = mlines.Line2D([0], [0],color='k', label=r'SMLM (Uniform' + '\n' + 'illumination intensity)')
            C10_line = mlines.Line2D([0], [0],color='k', linestyle='--', label=r'Localization on ISM' + '\n' + 'reconstruction data')
            C11_line = mlines.Line2D([0], [0],color='k', linestyle=':', label=r'Localization on' + '\n' + 'Fourier reweighted ISM' + '\n' + 'reconstruction data')
    
            ax.legend(handles=[C0_line, C1_line, C2_line, C3_line, C8_line, C9_line, C10_line, C11_line], prop={'size':13}, loc='upper left', bbox_to_anchor=(1.04,1), labelspacing = 1)

        plt.savefig(('./Figures/S19/'+str(int(panel))) + '.svg')
    
    if plotnumber == 'S20':
        CRLBx, musum, signalsum, Bisum = load_4('xp_s_2pattern')
        CRLBx_SMLM, musum_SMLM, signalsum_SMLM, Bisum_SMLM = load_4('SMLM')
        s_list, xrange = load_xp_s('xp_s_2pattern')
        
        fig, ax = plt.subplots(1,1)
        fig.set_size_inches(5, 5, forward=True)
        
        if panel == 1:
            ax.semilogy((xrange-theta[0])/(10**-9),CRLBx.T*10**9)
            ax.semilogy((xrange-theta[0])/(10**-9),CRLBx_SMLM*10**9*np.ones(np.size(xrange)), color='k')
            ax.semilogy((xrange-theta[0])/(10**-9),CRLBx_ISM*10**9*np.ones(np.size(xrange)), color='k', linestyle='--')
            ax.semilogy((xrange-theta[0])/(10**-9),CRLBx_FR*10**9*np.ones(np.size(xrange)), color='k', linestyle=':')
            ax.set_xlabel('Emitter-focus distance ' + r'$(x_f - \theta_x)$' + ' (nm)')
            ax.set_ylabel(r'$x$'+'-Cramér-Rao lower bound (nm)')
            ax.grid()
            ax.set_xlim(left=-130, right = 130)
            ax.set_ylim(bottom = 0.6, top = 60)
            ax.tick_params(axis='both', which='major')
            ax.tick_params(axis='both', which='minor')
            
        if panel == 2:
            ax.plot((xrange-theta[0])/(10**-9),CRLBx_SMLM/CRLBx.T)
            ax.plot((xrange-theta[0])/(10**-9),imp_ISM*np.ones(np.size(xrange)), color='k', linestyle='--')
            ax.plot((xrange-theta[0])/(10**-9),imp_FR*np.ones(np.size(xrange)), color='k', linestyle=':')
            ax.set_xlabel('Emitter-focus distance ' + r'$(x_f - \theta_x)$' + ' (nm)')
            ax.set_ylabel('Improvement over SMLM in ' + r'$x$'+ '-direction')
            ax.grid()
            ax.set_xlim(left=-130, right = 130)
            ax.set_ylim(bottom = 0, top = 2.2)
            ax.tick_params(axis='both', which='major')
            ax.tick_params(axis='both', which='minor')                     

        if panel == 3:
            ax.plot((xrange-theta[0])/(10**-9),musum.T)
            ax.plot((xrange-theta[0])/(10**-9),musum_SMLM*np.ones(np.size(xrange)), color='k')
            ax.set_xlabel('Emitter-focus distance ' + r'$(x_f - \theta_x)$' + ' (nm)')
            ax.set_ylabel('Expected amount of photons')
            ax.grid()
            ax.set_xlim(left=-130, right = 130)
            ax.set_ylim(bottom = 0, top = 2900)
            ax.tick_params(axis='both', which='major')
            ax.tick_params(axis='both', which='minor')                  

        if panel == 4:
            ax.plot((xrange-theta[0])/(10**-9),signalsum.T)
            ax.plot((xrange-theta[0])/(10**-9),signalsum_SMLM*np.ones(np.size(xrange)), color='k')
            ax.set_xlabel('Emitter-focus distance ' + r'$(x_f - \theta_x)$' + ' (nm)')
            ax.set_ylabel('Illumination-intensity normalized' + '\n' + 'amount of signal photons')
            ax.grid()
            ax.set_xlim(left=-130, right = 130)
            ax.set_ylim(bottom = 0, top = 2100)
            ax.tick_params(axis='both', which='major')
            ax.tick_params(axis='both', which='minor')   

        if panel == 5:
            ax.plot((xrange-theta[0])/(10**-9),Bisum.T/100)
            ax.plot((xrange-theta[0])/(10**-9),Bisum_SMLM/100*np.ones(np.size(xrange)), color='k')
            ax.set_xlabel('Emitter-focus distance ' + r'$(x_f - \theta_x)$' + ' (nm)')
            ax.set_ylabel('Illumination-intensity normalized average' + '\n' + 'amount of background photons per pixel')
            ax.grid()
            ax.set_xlim(left=-130, right = 130)
            ax.set_ylim(bottom = 0, top = 9)
            ax.tick_params(axis='both', which='major')
            ax.tick_params(axis='both', which='minor')   

        if panel == 6:
            CRLBx, musum, signalsum, Bisum = load_4('thetaI_thetab_2pattern')
            thetab_list, thetaIrange = load_thetaI_thetab('thetaI_thetab_2pattern')
            
            ax.loglog(thetaIrange,CRLBx.T*10**9)
            ax.set_xlabel('Expected signal photon budget ' + r'$\theta_I$' + ' (photons)')
            ax.set_ylabel(r'$x$'+'-Cramér-Rao lower bound (nm)')
            ax.grid()
            ax.set_xlim(left=100, right = 10000)
            ax.set_ylim(bottom = 0.6, top = 60)
            ax.tick_params(axis='both', which='major')
            ax.tick_params(axis='both', which='minor')  
            
            #Legend
            C0_line = mlines.Line2D([0], [0],color='C0', label=(r'$\theta_b = 0$'))
            C1_line = mlines.Line2D([0], [0],color='C1', label=(r'$\theta_b = 1$'))
            C2_line = mlines.Line2D([0], [0],color='C2', label=(r'$\theta_b = 4$'))
            C3_line = mlines.Line2D([0], [0],color='C3', label=(r'$\theta_b = 8$'))
            C8_line = mlines.Line2D([0], [0],color='C4', label=(r'$\theta_b = 16$'))
            C9_line = mlines.Line2D([0], [0],color='None', label=('(in phot/px.)'))
    
            ax.legend(handles=[C0_line, C1_line, C2_line, C3_line, C8_line, C9_line], ncol = 2, prop={'size':13}, labelspacing = 1)    
            
        if panel == 0:
            #Legend
            C0_line = mlines.Line2D([0], [0],color='C0', label=r'$s = 0\sigma_{\mathrm{PSF}}$')
            C1_line = mlines.Line2D([0], [0],color='C1', label=r'$s= 1\sigma_{\mathrm{PSF}}$')
            C2_line = mlines.Line2D([0], [0],color='C2', label=r'$s= 2\sigma_{\mathrm{PSF}}$')
            C3_line = mlines.Line2D([0], [0],color='C3', label=r'$s= 3\sigma_{\mathrm{PSF}}$')
            C4_line = mlines.Line2D([0], [0],color='C4', label=r'$s= 4\sigma_{\mathrm{PSF}}$')
            C8_line = mlines.Line2D([0], [0],color='C5', label=r'$s = 5\sigma_{\mathrm{PSF}}$')
            C9_line = mlines.Line2D([0], [0],color='k', label=r'SMLM (Uniform' + '\n' + 'illumination intensity)')
            C10_line = mlines.Line2D([0], [0],color='k', linestyle='--', label=r'Localization on ISM' + '\n' + 'reconstruction data')
            C11_line = mlines.Line2D([0], [0],color='k', linestyle=':', label=r'Localization on' + '\n' + 'Fourier reweighted ISM' + '\n' + 'reconstruction data')
    
            ax.legend(handles=[C0_line, C1_line, C2_line, C3_line, C4_line, C8_line, C9_line, C10_line, C11_line], prop={'size':13}, loc='upper left', bbox_to_anchor=(1.04,1), labelspacing = 1)

        plt.savefig(('./Figures/S20/'+str(int(panel))) + '.svg')

    if plotnumber == 'S21':
        CRLBx, musum, signalsum, Bisum = load_4('yp_xs_2pattern')
        CRLBx_SMLM, musum_SMLM, signalsum_SMLM, Bisum_SMLM = load_4('SMLM')
        s_list, xrange = load_xp_s('yp_xs_2pattern')
        
        fig, ax = plt.subplots(1,1)
        fig.set_size_inches(5, 5, forward=True)
        
        if panel == 1:
            ax.semilogy((xrange-theta[0])/(10**-9),CRLBx.T*10**9)
            ax.semilogy((xrange-theta[0])/(10**-9),CRLBx_SMLM*10**9*np.ones(np.size(xrange)), color='k')
            ax.semilogy((xrange-theta[0])/(10**-9),CRLBx_ISM*10**9*np.ones(np.size(xrange)), color='k', linestyle='--')
            ax.semilogy((xrange-theta[0])/(10**-9),CRLBx_FR*10**9*np.ones(np.size(xrange)), color='k', linestyle=':')
            ax.set_xlabel('Emitter-focus distance ' + r'$(y_f - \theta_y)$' + ' (nm)')
            ax.set_ylabel(r'$x$'+'-Cramér-Rao lower bound (nm)')
            ax.grid()
            ax.set_xlim(left=-130, right = 130)
            ax.set_ylim(bottom = 0.6, top = 60)
            ax.tick_params(axis='both', which='major')
            ax.tick_params(axis='both', which='minor')
            
        if panel == 2:
            ax.plot((xrange-theta[0])/(10**-9),CRLBx_SMLM/CRLBx.T)
            ax.plot((xrange-theta[0])/(10**-9),imp_ISM*np.ones(np.size(xrange)), color='k', linestyle='--')
            ax.plot((xrange-theta[0])/(10**-9),imp_FR*np.ones(np.size(xrange)), color='k', linestyle=':')
            ax.set_xlabel('Emitter-focus distance ' + r'$(y_f - \theta_y)$' + ' (nm)')
            ax.set_ylabel('Improvement over SMLM in ' + r'$x$'+ '-direction')
            ax.grid()
            ax.set_xlim(left=-130, right = 130)
            ax.set_ylim(bottom = 0, top = 2.2)
            ax.tick_params(axis='both', which='major')
            ax.tick_params(axis='both', which='minor')                     

        if panel == 3:
            ax.plot((xrange-theta[0])/(10**-9),musum.T)
            ax.plot((xrange-theta[0])/(10**-9),musum_SMLM*np.ones(np.size(xrange)), color='k')
            ax.set_xlabel('Emitter-focus distance ' + r'$(y_f - \theta_y)$' + ' (nm)')
            ax.set_ylabel('Expected amount of photons')
            ax.grid()
            ax.set_xlim(left=-130, right = 130)
            ax.set_ylim(bottom = 0, top = 2900)
            ax.tick_params(axis='both', which='major')
            ax.tick_params(axis='both', which='minor')                  

        if panel == 4:
            ax.plot((xrange-theta[0])/(10**-9),signalsum.T)
            ax.plot((xrange-theta[0])/(10**-9),signalsum_SMLM*np.ones(np.size(xrange)), color='k')
            ax.set_xlabel('Emitter-focus distance ' + r'$(y_f - \theta_y)$' + ' (nm)')
            ax.set_ylabel('Illumination-intensity normalized' + '\n' + 'amount of signal photons')
            ax.grid()
            ax.set_xlim(left=-130, right = 130)
            ax.set_ylim(bottom = 0, top = 2100)
            ax.tick_params(axis='both', which='major')
            ax.tick_params(axis='both', which='minor')   

        if panel == 5:
            ax.plot((xrange-theta[0])/(10**-9),Bisum.T/100)
            ax.plot((xrange-theta[0])/(10**-9),Bisum_SMLM/100*np.ones(np.size(xrange)), color='k')
            ax.set_xlabel('Emitter-focus distance ' + r'$(y_f - \theta_y)$' + ' (nm)')
            ax.set_ylabel('Illumination-intensity normalized average' + '\n' + 'amount of background photons per pixel')
            ax.grid()
            ax.set_xlim(left=-130, right = 130)
            ax.set_ylim(bottom = 0, top = 9)
            ax.tick_params(axis='both', which='major')
            ax.tick_params(axis='both', which='minor')   

        if panel == 0:
            #Legend
            C0_line = mlines.Line2D([0], [0],color='C0', label=r'$s = 0\sigma_{\mathrm{PSF}}$')
            C1_line = mlines.Line2D([0], [0],color='C1', label=r'$s= 1\sigma_{\mathrm{PSF}}$')
            C2_line = mlines.Line2D([0], [0],color='C2', label=r'$s= 2\sigma_{\mathrm{PSF}}$')
            C3_line = mlines.Line2D([0], [0],color='C3', label=r'$s= 3\sigma_{\mathrm{PSF}}$')
            C4_line = mlines.Line2D([0], [0],color='C4', label=r'$s= 4\sigma_{\mathrm{PSF}}$')
            C8_line = mlines.Line2D([0], [0],color='C5', label=r'$s = 5\sigma_{\mathrm{PSF}}$')
            C9_line = mlines.Line2D([0], [0],color='k', label=r'SMLM (Uniform' + '\n' + 'illumination intensity)')
            C10_line = mlines.Line2D([0], [0],color='k', linestyle='--', label=r'Localization on ISM' + '\n' + 'reconstruction data')
            C11_line = mlines.Line2D([0], [0],color='k', linestyle=':', label=r'Localization on' + '\n' + 'Fourier reweighted ISM' + '\n' + 'reconstruction data')
    
            ax.legend(handles=[C0_line, C1_line, C2_line, C3_line, C4_line, C8_line, C9_line, C10_line, C11_line], prop={'size':13}, loc='upper left', bbox_to_anchor=(1.04,1), labelspacing = 1)

        plt.savefig(('./Figures/S21/'+str(int(panel))) + '.svg')

    if plotnumber == 'S22':
        CRLBx, musum, signalsum, Bisum = load_4('xp_s_2pattern_nopinhole')
        CRLBx_SMLM, musum_SMLM, signalsum_SMLM, Bisum_SMLM = load_4('SMLM')
        s_list, xrange = load_xp_s('xp_s_2pattern_nopinhole')
        
        fig, ax = plt.subplots(1,1)
        fig.set_size_inches(5, 5, forward=True)
        
        if panel == 1:
            ax.semilogy((xrange-theta[0])/(10**-9),CRLBx.T*10**9)
            ax.semilogy((xrange-theta[0])/(10**-9),CRLBx_SMLM*10**9*np.ones(np.size(xrange)), color='k')
            ax.semilogy((xrange-theta[0])/(10**-9),CRLBx_ISM*10**9*np.ones(np.size(xrange)), color='k', linestyle='--')
            ax.semilogy((xrange-theta[0])/(10**-9),CRLBx_FR*10**9*np.ones(np.size(xrange)), color='k', linestyle=':')
            ax.set_xlabel('Emitter-focus distance ' + r'$(x_f - \theta_x)$' + ' (nm)')
            ax.set_ylabel(r'$x$'+'-Cramér-Rao lower bound (nm)')
            ax.grid()
            ax.set_xlim(left=-130, right = 130)
            ax.set_ylim(bottom = 0.6, top = 60)
            ax.tick_params(axis='both', which='major')
            ax.tick_params(axis='both', which='minor')
            
        if panel == 2:
            ax.plot((xrange-theta[0])/(10**-9),CRLBx_SMLM/CRLBx.T)
            ax.plot((xrange-theta[0])/(10**-9),imp_ISM*np.ones(np.size(xrange)), color='k', linestyle='--')
            ax.plot((xrange-theta[0])/(10**-9),imp_FR*np.ones(np.size(xrange)), color='k', linestyle=':')
            ax.set_xlabel('Emitter-focus distance ' + r'$(x_f - \theta_x)$' + ' (nm)')
            ax.set_ylabel('Improvement over SMLM in ' + r'$x$'+ '-direction')
            ax.grid()
            ax.set_xlim(left=-130, right = 130)
            ax.set_ylim(bottom = 0, top = 2.2)
            ax.tick_params(axis='both', which='major')
            ax.tick_params(axis='both', which='minor')                     

        if panel == 3:
            ax.plot((xrange-theta[0])/(10**-9),musum.T)
            ax.plot((xrange-theta[0])/(10**-9),musum_SMLM*np.ones(np.size(xrange)), color='k')
            ax.set_xlabel('Emitter-focus distance ' + r'$(x_f - \theta_x)$' + ' (nm)')
            ax.set_ylabel('Expected amount of photons')
            ax.grid()
            ax.set_xlim(left=-130, right = 130)
            ax.set_ylim(bottom = 0, top = 2900)
            ax.tick_params(axis='both', which='major')
            ax.tick_params(axis='both', which='minor')                  

        if panel == 4:
            ax.plot((xrange-theta[0])/(10**-9),signalsum.T)
            ax.plot((xrange-theta[0])/(10**-9),signalsum_SMLM*np.ones(np.size(xrange)), color='k')
            ax.set_xlabel('Emitter-focus distance ' + r'$(x_f - \theta_x)$' + ' (nm)')
            ax.set_ylabel('Illumination-intensity normalized' + '\n' + 'amount of signal photons')
            ax.grid()
            ax.set_xlim(left=-130, right = 130)
            ax.set_ylim(bottom = 0, top = 2100)
            ax.tick_params(axis='both', which='major')
            ax.tick_params(axis='both', which='minor')   

        if panel == 5:
            ax.plot((xrange-theta[0])/(10**-9),Bisum.T/100)
            ax.plot((xrange-theta[0])/(10**-9),Bisum_SMLM/100*np.ones(np.size(xrange)), color='k')
            ax.set_xlabel('Emitter-focus distance ' + r'$(x_f - \theta_x)$' + ' (nm)')
            ax.set_ylabel('Illumination-intensity normalized average' + '\n' + 'amount of background photons per pixel')
            ax.grid()
            ax.set_xlim(left=-130, right = 130)
            ax.set_ylim(bottom = 0, top = 9)
            ax.tick_params(axis='both', which='major')
            ax.tick_params(axis='both', which='minor')   

        if panel == 6:
            CRLBx, musum, signalsum, Bisum = load_4('thetaI_thetab_2pattern_nopinhole')
            thetab_list, thetaIrange = load_thetaI_thetab('thetaI_thetab_2pattern_nopinhole')
            
            ax.loglog(thetaIrange,CRLBx.T*10**9)
            ax.set_xlabel('Expected signal photon budget ' + r'$\theta_I$' + ' (photons)')
            ax.set_ylabel(r'$x$'+'-Cramér-Rao lower bound (nm)')
            ax.grid()
            ax.set_xlim(left=100, right = 10000)
            ax.set_ylim(bottom = 0.6, top = 60)
            ax.tick_params(axis='both', which='major')
            ax.tick_params(axis='both', which='minor')  
            
            #Legend
            C0_line = mlines.Line2D([0], [0],color='C0', label=(r'$\theta_b = 0$'))
            C1_line = mlines.Line2D([0], [0],color='C1', label=(r'$\theta_b = 1$'))
            C2_line = mlines.Line2D([0], [0],color='C2', label=(r'$\theta_b = 4$'))
            C3_line = mlines.Line2D([0], [0],color='C3', label=(r'$\theta_b = 8$'))
            C8_line = mlines.Line2D([0], [0],color='C4', label=(r'$\theta_b = 16$'))
            C9_line = mlines.Line2D([0], [0],color='None', label=('(in phot/px.)'))
    
            ax.legend(handles=[C0_line, C1_line, C2_line, C3_line, C8_line, C9_line], ncol = 2, prop={'size':13}, labelspacing = 1)    

        if panel == 0:
            #Legend
            C0_line = mlines.Line2D([0], [0],color='C0', label=r'$s = 0\sigma_{\mathrm{PSF}}$')
            C1_line = mlines.Line2D([0], [0],color='C1', label=r'$s= 1\sigma_{\mathrm{PSF}}$')
            C2_line = mlines.Line2D([0], [0],color='C2', label=r'$s= 2\sigma_{\mathrm{PSF}}$')
            C3_line = mlines.Line2D([0], [0],color='C3', label=r'$s= 3\sigma_{\mathrm{PSF}}$')
            C4_line = mlines.Line2D([0], [0],color='C4', label=r'$s= 4\sigma_{\mathrm{PSF}}$')
            C8_line = mlines.Line2D([0], [0],color='C5', label=r'$s = 5\sigma_{\mathrm{PSF}}$')
            C9_line = mlines.Line2D([0], [0],color='k', label=r'SMLM (Uniform' + '\n' + 'illumination intensity)')
            C10_line = mlines.Line2D([0], [0],color='k', linestyle='--', label=r'Localization on ISM' + '\n' + 'reconstruction data')
            C11_line = mlines.Line2D([0], [0],color='k', linestyle=':', label=r'Localization on' + '\n' + 'Fourier reweighted ISM' + '\n' + 'reconstruction data')
    
            ax.legend(handles=[C0_line, C1_line, C2_line, C3_line, C4_line, C8_line, C9_line, C10_line, C11_line], prop={'size':13}, loc='upper left', bbox_to_anchor=(1.04,1), labelspacing = 1)

        plt.savefig(('./Figures/S22/'+str(int(panel))) + '.svg')

    if plotnumber == 'S23':
        CRLBx, musum, signalsum, Bisum = load_4('yp_s_2pattern')
        CRLBx_SMLM, musum_SMLM, signalsum_SMLM, Bisum_SMLM = load_4('SMLM')
        s_list, xrange = load_xp_s('yp_s_2pattern')
        
        fig, ax = plt.subplots(1,1)
        fig.set_size_inches(5, 5, forward=True)
        
        if panel == 1:
            ax.semilogy((xrange-theta[0])/(10**-9),CRLBx.T*10**9)
            ax.semilogy((xrange-theta[0])/(10**-9),CRLBx_SMLM*10**9*np.ones(np.size(xrange)), color='k')
            ax.semilogy((xrange-theta[0])/(10**-9),CRLBx_ISM*10**9*np.ones(np.size(xrange)), color='k', linestyle='--')
            ax.semilogy((xrange-theta[0])/(10**-9),CRLBx_FR*10**9*np.ones(np.size(xrange)), color='k', linestyle=':')
            ax.set_xlabel('Emitter-focus distance ' + r'$(y_f - \theta_y)$' + ' (nm)')
            ax.set_ylabel(r'$x$'+'-Cramér-Rao lower bound (nm)')
            ax.grid()
            ax.set_xlim(left=-130, right = 130)
            ax.set_ylim(bottom = 0.6, top = 60)
            ax.tick_params(axis='both', which='major')
            ax.tick_params(axis='both', which='minor')
            
        if panel == 2:
            ax.plot((xrange-theta[0])/(10**-9),CRLBx_SMLM/CRLBx.T)
            ax.plot((xrange-theta[0])/(10**-9),imp_ISM*np.ones(np.size(xrange)), color='k', linestyle='--')
            ax.plot((xrange-theta[0])/(10**-9),imp_FR*np.ones(np.size(xrange)), color='k', linestyle=':')
            ax.set_xlabel('Emitter-focus distance ' + r'$(y_f - \theta_y)$' + ' (nm)')
            ax.set_ylabel('Improvement over SMLM in ' + r'$x$'+ '-direction')
            ax.grid()
            ax.set_xlim(left=-130, right = 130)
            ax.set_ylim(bottom = 0, top = 2.2)
            ax.tick_params(axis='both', which='major')
            ax.tick_params(axis='both', which='minor')                     

        if panel == 3:
            ax.plot((xrange-theta[0])/(10**-9),musum.T)
            ax.plot((xrange-theta[0])/(10**-9),musum_SMLM*np.ones(np.size(xrange)), color='k')
            ax.set_xlabel('Emitter-focus distance ' + r'$(y_f - \theta_y)$' + ' (nm)')
            ax.set_ylabel('Expected amount of photons')
            ax.grid()
            ax.set_xlim(left=-130, right = 130)
            ax.set_ylim(bottom = 0, top = 2900)
            ax.tick_params(axis='both', which='major')
            ax.tick_params(axis='both', which='minor')                  

        if panel == 4:
            ax.plot((xrange-theta[0])/(10**-9),signalsum.T)
            ax.plot((xrange-theta[0])/(10**-9),signalsum_SMLM*np.ones(np.size(xrange)), color='k')
            ax.set_xlabel('Emitter-focus distance ' + r'$(y_f - \theta_y)$' + ' (nm)')
            ax.set_ylabel('Illumination-intensity normalized' + '\n' + 'amount of signal photons')
            ax.grid()
            ax.set_xlim(left=-130, right = 130)
            ax.set_ylim(bottom = 0, top = 2100)
            ax.tick_params(axis='both', which='major')
            ax.tick_params(axis='both', which='minor')   

        if panel == 5:
            ax.plot((xrange-theta[0])/(10**-9),Bisum.T/100)
            ax.plot((xrange-theta[0])/(10**-9),Bisum_SMLM/100*np.ones(np.size(xrange)), color='k')
            ax.set_xlabel('Emitter-focus distance ' + r'$(y_f - \theta_y)$' + ' (nm)')
            ax.set_ylabel('Illumination-intensity normalized average' + '\n' + 'amount of background photons per pixel')
            ax.grid()
            ax.set_xlim(left=-130, right = 130)
            ax.set_ylim(bottom = 0, top = 9)
            ax.tick_params(axis='both', which='major')
            ax.tick_params(axis='both', which='minor')   

        if panel == 6:
            CRLBx, musum, signalsum, Bisum = load_4('thetaI_thetab_2pattern_rot')
            thetab_list, thetaIrange = load_thetaI_thetab('thetaI_thetab_2pattern_rot')
            
            ax.loglog(thetaIrange,CRLBx.T*10**9)
            ax.set_xlabel('Expected signal photon budget ' + r'$\theta_I$' + ' (photons)')
            ax.set_ylabel(r'$x$'+'-Cramér-Rao lower bound (nm)')
            ax.grid()
            ax.set_xlim(left=100, right = 10000)
            ax.set_ylim(bottom = 0.6, top = 200)
            ax.tick_params(axis='both', which='major')
            ax.tick_params(axis='both', which='minor')  
            
            #Legend
            C0_line = mlines.Line2D([0], [0],color='C0', label=(r'$\theta_b = 0$'))
            C1_line = mlines.Line2D([0], [0],color='C1', label=(r'$\theta_b = 1$'))
            C2_line = mlines.Line2D([0], [0],color='C2', label=(r'$\theta_b = 4$'))
            C3_line = mlines.Line2D([0], [0],color='C3', label=(r'$\theta_b = 8$'))
            C8_line = mlines.Line2D([0], [0],color='C4', label=(r'$\theta_b = 16$'))
            C9_line = mlines.Line2D([0], [0],color='None', label=('(in phot/px.)'))
    
            ax.legend(handles=[C0_line, C1_line, C2_line, C3_line, C8_line, C9_line], ncol = 2, prop={'size':13}, labelspacing = 1)    
            
        if panel == 0:
            #Legend
            C0_line = mlines.Line2D([0], [0],color='C0', label=r'$s = 0\sigma_{\mathrm{PSF}}$')
            C1_line = mlines.Line2D([0], [0],color='C1', label=r'$s= 1\sigma_{\mathrm{PSF}}$')
            C2_line = mlines.Line2D([0], [0],color='C2', label=r'$s= 2\sigma_{\mathrm{PSF}}$')
            C3_line = mlines.Line2D([0], [0],color='C3', label=r'$s= 3\sigma_{\mathrm{PSF}}$')
            C4_line = mlines.Line2D([0], [0],color='C4', label=r'$s= 4\sigma_{\mathrm{PSF}}$')
            C8_line = mlines.Line2D([0], [0],color='C5', label=r'$s = 5\sigma_{\mathrm{PSF}}$')
            C9_line = mlines.Line2D([0], [0],color='k', label=r'SMLM (Uniform' + '\n' + 'illumination intensity)')
            C10_line = mlines.Line2D([0], [0],color='k', linestyle='--', label=r'Localization on ISM' + '\n' + 'reconstruction data')
            C11_line = mlines.Line2D([0], [0],color='k', linestyle=':', label=r'Localization on' + '\n' + 'Fourier reweighted ISM' + '\n' + 'reconstruction data')
    
            ax.legend(handles=[C0_line, C1_line, C2_line, C3_line, C4_line, C8_line, C9_line, C10_line, C11_line], prop={'size':13}, loc='upper left', bbox_to_anchor=(1.04,1), labelspacing = 1)

        plt.savefig(('./Figures/S23/'+str(int(panel))) + '.svg')

    if plotnumber == 'S24':
        CRLBx, musum, signalsum, Bisum = load_4('xp_r_triangulation_rot_nocenter_ext')
        CRLBx_SMLM, musum_SMLM, signalsum_SMLM, Bisum_SMLM = load_4('SMLM')
        s_list, xrange = load_xp_r('xp_r_triangulation_rot_nocenter_ext')
        
        fig, ax = plt.subplots(1,1)
        fig.set_size_inches(5, 5, forward=True)
        
        if panel == 1:
            ax.semilogy((xrange-theta[0])/(10**-9),CRLBx.T[:,0:4]*10**9)
            ax.semilogy((xrange-theta[0])/(10**-9),CRLBx_SMLM*10**9*np.ones(np.size(xrange)), color='k')
            ax.semilogy((xrange-theta[0])/(10**-9),CRLBx_ISM*10**9*np.ones(np.size(xrange)), color='k', linestyle='--')
            ax.semilogy((xrange-theta[0])/(10**-9),CRLBx_FR*10**9*np.ones(np.size(xrange)), color='k', linestyle=':')
            ax.set_xlabel('Emitter-focus distance ' + r'$(x_f - \theta_x)$' + ' (nm)')
            ax.set_ylabel(r'$x$'+'-Cramér-Rao lower bound (nm)')
            ax.grid()
            ax.set_xlim(left=-130, right = 130)
            ax.set_ylim(bottom = 0.6, top = 10000)
            ax.tick_params(axis='both', which='major')
            ax.tick_params(axis='both', which='minor')
            
        if panel == 2:
            ax.plot((xrange-theta[0])/(10**-9),CRLBx_SMLM/CRLBx.T[:,0:4])
            ax.plot((xrange-theta[0])/(10**-9),imp_ISM*np.ones(np.size(xrange)), color='k', linestyle='--')
            ax.plot((xrange-theta[0])/(10**-9),imp_FR*np.ones(np.size(xrange)), color='k', linestyle=':')
            ax.set_xlabel('Emitter-focus distance ' + r'$(x_f - \theta_x)$' + ' (nm)')
            ax.set_ylabel('Improvement over SMLM in ' + r'$x$'+ '-direction')
            ax.grid()
            ax.set_xlim(left=-130, right = 130)
            ax.set_ylim(bottom = 0, top = 2.2)
            ax.tick_params(axis='both', which='major')
            ax.tick_params(axis='both', which='minor')                     

        if panel == 3:
            ax.plot((xrange-theta[0])/(10**-9),musum.T[:,0:4])
            ax.plot((xrange-theta[0])/(10**-9),musum_SMLM*np.ones(np.size(xrange)), color='k')
            ax.set_xlabel('Emitter-focus distance ' + r'$(x_f - \theta_x)$' + ' (nm)')
            ax.set_ylabel('Expected amount of photons')
            ax.grid()
            ax.set_xlim(left=-130, right = 130)
            ax.set_ylim(bottom = 0, top = 2900)
            ax.tick_params(axis='both', which='major')
            ax.tick_params(axis='both', which='minor')                  

        if panel == 4:
            ax.plot((xrange-theta[0])/(10**-9),signalsum.T[:,0:4])
            ax.plot((xrange-theta[0])/(10**-9),signalsum_SMLM*np.ones(np.size(xrange)), color='k')
            ax.set_xlabel('Emitter-focus distance ' + r'$(x_f - \theta_x)$' + ' (nm)')
            ax.set_ylabel('Illumination-intensity normalized' + '\n' + 'amount of signal photons')
            ax.grid()
            ax.set_xlim(left=-130, right = 130)
            ax.set_ylim(bottom = 0, top = 2100)
            ax.tick_params(axis='both', which='major')
            ax.tick_params(axis='both', which='minor')   

        if panel == 5:
            ax.plot((xrange-theta[0])/(10**-9),Bisum.T[:,0:4]/100)
            ax.plot((xrange-theta[0])/(10**-9),Bisum_SMLM/100*np.ones(np.size(xrange)), color='k')
            ax.set_xlabel('Emitter-focus distance ' + r'$(x_f - \theta_x)$' + ' (nm)')
            ax.set_ylabel('Illumination-intensity normalized average' + '\n' + 'amount of background photons per pixel')
            ax.grid()
            ax.set_xlim(left=-130, right = 130)
            ax.set_ylim(bottom = 0, top = 9)
            ax.tick_params(axis='both', which='major')
            ax.tick_params(axis='both', which='minor')   

        if panel == 6:
            CRLBx, musum, signalsum, Bisum = load_4('thetaI_thetab_triangulation_rot_nocenter')
            thetab_list, thetaIrange = load_thetaI_thetab('thetaI_thetab_triangulation_rot_nocenter')
            
            ax.loglog(thetaIrange,CRLBx.T*10**9)
            ax.set_xlabel('Expected signal photon budget ' + r'$\theta_I$' + ' (photons)')
            ax.set_ylabel(r'$x$'+'-Cramér-Rao lower bound (nm)')
            ax.grid()
            ax.set_xlim(left=100, right = 10000)
            ax.set_ylim(bottom = 0.1, top = 60)
            ax.tick_params(axis='both', which='major')
            ax.tick_params(axis='both', which='minor')  
            
            #Legend
            C0_line = mlines.Line2D([0], [0],color='C0', label=(r'$\theta_b = 0$'))
            C1_line = mlines.Line2D([0], [0],color='C1', label=(r'$\theta_b = 1$'))
            C2_line = mlines.Line2D([0], [0],color='C2', label=(r'$\theta_b = 4$'))
            C3_line = mlines.Line2D([0], [0],color='C3', label=(r'$\theta_b = 8$'))
            C8_line = mlines.Line2D([0], [0],color='C4', label=(r'$\theta_b = 16$'))
            C9_line = mlines.Line2D([0], [0],color='None', label=('(in phot/px.)'))
    
            ax.legend(handles=[C0_line, C1_line, C2_line, C3_line, C8_line, C9_line], ncol = 2, prop={'size':13}, labelspacing = 1) 

        if panel == 0:
            #Legend
            C0_line = mlines.Line2D([0], [0],color='C0', label=r'$r = 1\sigma_{\mathrm{PSF}}$')
            C1_line = mlines.Line2D([0], [0],color='C1', label=r'$r = 2\sigma_{\mathrm{PSF}}$')
            C2_line = mlines.Line2D([0], [0],color='C2', label=r'$r = 3\sigma_{\mathrm{PSF}}$')
            C3_line = mlines.Line2D([0], [0],color='C3', label=r'$r = 4\sigma_{\mathrm{PSF}}$')
            C9_line = mlines.Line2D([0], [0],color='k', label=r'SMLM (Uniform' + '\n' + 'illumination intensity)')
            C10_line = mlines.Line2D([0], [0],color='k', linestyle='--', label=r'Localization on ISM' + '\n' + 'reconstruction data')
            C11_line = mlines.Line2D([0], [0],color='k', linestyle=':', label=r'Localization on' + '\n' + 'Fourier reweighted ISM' + '\n' + 'reconstruction data')
    
            ax.legend(handles=[C0_line, C1_line, C2_line, C3_line, C9_line, C10_line, C11_line], prop={'size':13}, loc='upper left', bbox_to_anchor=(1.04,1), labelspacing = 1)

        plt.savefig(('./Figures/S24/'+str(int(panel))) + '.svg')
    
    if plotnumber == 'S25':
        CRLBx, musum, signalsum, Bisum = load_4('xp_r_triangulation_nocenter_ext')
        CRLBx_SMLM, musum_SMLM, signalsum_SMLM, Bisum_SMLM = load_4('SMLM')
        s_list, xrange = load_xp_r('xp_r_triangulation_nocenter_ext')
        
        fig, ax = plt.subplots(1,1)
        fig.set_size_inches(5, 5, forward=True)
        
        if panel == 1:
            ax.semilogy((xrange-theta[0])/(10**-9),CRLBx.T[:,0:4]*10**9)
            ax.semilogy((xrange-theta[0])/(10**-9),CRLBx_SMLM*10**9*np.ones(np.size(xrange)), color='k')
            ax.semilogy((xrange-theta[0])/(10**-9),CRLBx_ISM*10**9*np.ones(np.size(xrange)), color='k', linestyle='--')
            ax.semilogy((xrange-theta[0])/(10**-9),CRLBx_FR*10**9*np.ones(np.size(xrange)), color='k', linestyle=':')
            ax.set_xlabel('Emitter-focus distance ' + r'$(x_f - \theta_x)$' + ' (nm)')
            ax.set_ylabel(r'$x$'+'-Cramér-Rao lower bound (nm)')
            ax.grid()
            ax.set_xlim(left=-130, right = 130)
            ax.set_ylim(bottom = 0.6, top = 10000)
            ax.tick_params(axis='both', which='major')
            ax.tick_params(axis='both', which='minor')
            
        if panel == 2:
            ax.plot((xrange-theta[0])/(10**-9),CRLBx_SMLM/CRLBx.T[:,0:4])
            ax.plot((xrange-theta[0])/(10**-9),imp_ISM*np.ones(np.size(xrange)), color='k', linestyle='--')
            ax.plot((xrange-theta[0])/(10**-9),imp_FR*np.ones(np.size(xrange)), color='k', linestyle=':')
            ax.set_xlabel('Emitter-focus distance ' + r'$(x_f - \theta_x)$' + ' (nm)')
            ax.set_ylabel('Improvement over SMLM in ' + r'$x$'+ '-direction')
            ax.grid()
            ax.set_xlim(left=-130, right = 130)
            ax.set_ylim(bottom = 0, top = 2.2)
            ax.tick_params(axis='both', which='major')
            ax.tick_params(axis='both', which='minor')                     
    
        if panel == 3:
            ax.plot((xrange-theta[0])/(10**-9),musum.T[:,0:4])
            ax.plot((xrange-theta[0])/(10**-9),musum_SMLM*np.ones(np.size(xrange)), color='k')
            ax.set_xlabel('Emitter-focus distance ' + r'$(x_f - \theta_x)$' + ' (nm)')
            ax.set_ylabel('Expected amount of photons')
            ax.grid()
            ax.set_xlim(left=-130, right = 130)
            ax.set_ylim(bottom = 0, top = 2900)
            ax.tick_params(axis='both', which='major')
            ax.tick_params(axis='both', which='minor')                  
    
        if panel == 4:
            ax.plot((xrange-theta[0])/(10**-9),signalsum.T[:,0:4])
            ax.plot((xrange-theta[0])/(10**-9),signalsum_SMLM*np.ones(np.size(xrange)), color='k')
            ax.set_xlabel('Emitter-focus distance ' + r'$(x_f - \theta_x)$' + ' (nm)')
            ax.set_ylabel('Illumination-intensity normalized' + '\n' + 'amount of signal photons')
            ax.grid()
            ax.set_xlim(left=-130, right = 130)
            ax.set_ylim(bottom = 0, top = 2100)
            ax.tick_params(axis='both', which='major')
            ax.tick_params(axis='both', which='minor')   
    
        if panel == 5:
            ax.plot((xrange-theta[0])/(10**-9),Bisum.T[:,0:4]/100)
            ax.plot((xrange-theta[0])/(10**-9),Bisum_SMLM/100*np.ones(np.size(xrange)), color='k')
            ax.set_xlabel('Emitter-focus distance ' + r'$(x_f - \theta_x)$' + ' (nm)')
            ax.set_ylabel('Illumination-intensity normalized average' + '\n' + 'amount of background photons per pixel')
            ax.grid()
            ax.set_xlim(left=-130, right = 130)
            ax.set_ylim(bottom = 0, top = 9)
            ax.tick_params(axis='both', which='major')
            ax.tick_params(axis='both', which='minor')   
    
        if panel == 6:
            CRLBx, musum, signalsum, Bisum = load_4('thetaI_thetab_triangulation_nocenter')
            thetab_list, thetaIrange = load_thetaI_thetab('thetaI_thetab_triangulation_nocenter')
            
            ax.loglog(thetaIrange,CRLBx.T*10**9)
            ax.set_xlabel('Expected signal photon budget ' + r'$\theta_I$' + ' (photons)')
            ax.set_ylabel(r'$x$'+'-Cramér-Rao lower bound (nm)')
            ax.grid()
            ax.set_xlim(left=100, right = 10000)
            ax.set_ylim(bottom = 0.1, top = 60)
            ax.tick_params(axis='both', which='major')
            ax.tick_params(axis='both', which='minor')  
            
            #Legend
            C0_line = mlines.Line2D([0], [0],color='C0', label=(r'$\theta_b = 0$'))
            C1_line = mlines.Line2D([0], [0],color='C1', label=(r'$\theta_b = 1$'))
            C2_line = mlines.Line2D([0], [0],color='C2', label=(r'$\theta_b = 4$'))
            C3_line = mlines.Line2D([0], [0],color='C3', label=(r'$\theta_b = 8$'))
            C8_line = mlines.Line2D([0], [0],color='C4', label=(r'$\theta_b = 16$'))
            C9_line = mlines.Line2D([0], [0],color='None', label=('(in phot/px.)'))

            ax.legend(handles=[C0_line, C1_line, C2_line, C3_line, C8_line, C9_line], ncol = 2, prop={'size':13}, labelspacing = 1) 

        if panel == 0:
            #Legend
            C0_line = mlines.Line2D([0], [0],color='C0', label=r'$r = 1\sigma_{\mathrm{PSF}}$')
            C1_line = mlines.Line2D([0], [0],color='C1', label=r'$r = 2\sigma_{\mathrm{PSF}}$')
            C2_line = mlines.Line2D([0], [0],color='C2', label=r'$r = 3\sigma_{\mathrm{PSF}}$')
            C3_line = mlines.Line2D([0], [0],color='C3', label=r'$r = 4\sigma_{\mathrm{PSF}}$')
            C9_line = mlines.Line2D([0], [0],color='k', label=r'SMLM (Uniform' + '\n' + 'illumination intensity)')
            C10_line = mlines.Line2D([0], [0],color='k', linestyle='--', label=r'Localization on ISM' + '\n' + 'reconstruction data')
            C11_line = mlines.Line2D([0], [0],color='k', linestyle=':', label=r'Localization on' + '\n' + 'Fourier reweighted ISM' + '\n' + 'reconstruction data')
    
            ax.legend(handles=[C0_line, C1_line, C2_line, C3_line, C9_line, C10_line, C11_line], prop={'size':13}, loc='upper left', bbox_to_anchor=(1.04,1), labelspacing = 1)
    
        plt.savefig(('./Figures/S25/'+str(int(panel))) + '.svg')

    if plotnumber == 'S26':
        CRLBx, musum, signalsum, Bisum = load_4('xp_r_triangulation_rot')
        CRLBx_SMLM, musum_SMLM, signalsum_SMLM, Bisum_SMLM = load_4('SMLM')
        s_list, xrange = load_xp_r('xp_r_triangulation_rot')
        
        fig, ax = plt.subplots(1,1)
        fig.set_size_inches(5, 5, forward=True)
        
        if panel == 1:
            ax.semilogy((xrange-theta[0])/(10**-9),CRLBx.T*10**9)
            ax.semilogy((xrange-theta[0])/(10**-9),CRLBx_SMLM*10**9*np.ones(np.size(xrange)), color='k')
            ax.semilogy((xrange-theta[0])/(10**-9),CRLBx_ISM*10**9*np.ones(np.size(xrange)), color='k', linestyle='--')
            ax.semilogy((xrange-theta[0])/(10**-9),CRLBx_FR*10**9*np.ones(np.size(xrange)), color='k', linestyle=':')
            ax.set_xlabel('Emitter-focus distance ' + r'$(x_f - \theta_x)$' + ' (nm)')
            ax.set_ylabel(r'$x$'+'-Cramér-Rao lower bound (nm)')
            ax.grid()
            ax.set_xlim(left=-130, right = 130)
            ax.set_ylim(bottom = 0.6, top = 60)
            ax.tick_params(axis='both', which='major')
            ax.tick_params(axis='both', which='minor')
            
        if panel == 2:
            ax.plot((xrange-theta[0])/(10**-9),CRLBx_SMLM/CRLBx.T)
            ax.plot((xrange-theta[0])/(10**-9),imp_ISM*np.ones(np.size(xrange)), color='k', linestyle='--')
            ax.plot((xrange-theta[0])/(10**-9),imp_FR*np.ones(np.size(xrange)), color='k', linestyle=':')
            ax.set_xlabel('Emitter-focus distance ' + r'$(x_f - \theta_x)$' + ' (nm)')
            ax.set_ylabel('Improvement over SMLM in ' + r'$x$'+ '-direction')
            ax.grid()
            ax.set_xlim(left=-130, right = 130)
            ax.set_ylim(bottom = 0, top = 2.2)
            ax.tick_params(axis='both', which='major')
            ax.tick_params(axis='both', which='minor')                     

        if panel == 3:
            ax.plot((xrange-theta[0])/(10**-9),musum.T)
            ax.plot((xrange-theta[0])/(10**-9),musum_SMLM*np.ones(np.size(xrange)), color='k')
            ax.set_xlabel('Emitter-focus distance ' + r'$(x_f - \theta_x)$' + ' (nm)')
            ax.set_ylabel('Expected amount of photons')
            ax.grid()
            ax.set_xlim(left=-130, right = 130)
            ax.set_ylim(bottom = 0, top = 2900)
            ax.tick_params(axis='both', which='major')
            ax.tick_params(axis='both', which='minor')                  

        if panel == 4:
            ax.plot((xrange-theta[0])/(10**-9),signalsum.T)
            ax.plot((xrange-theta[0])/(10**-9),signalsum_SMLM*np.ones(np.size(xrange)), color='k')
            ax.set_xlabel('Emitter-focus distance ' + r'$(x_f - \theta_x)$' + ' (nm)')
            ax.set_ylabel('Illumination-intensity normalized' + '\n' + 'amount of signal photons')
            ax.grid()
            ax.set_xlim(left=-130, right = 130)
            ax.set_ylim(bottom = 0, top = 2100)
            ax.tick_params(axis='both', which='major')
            ax.tick_params(axis='both', which='minor')   

        if panel == 5:
            ax.plot((xrange-theta[0])/(10**-9),Bisum.T/100)
            ax.plot((xrange-theta[0])/(10**-9),Bisum_SMLM/100*np.ones(np.size(xrange)), color='k')
            ax.set_xlabel('Emitter-focus distance ' + r'$(x_f - \theta_x)$' + ' (nm)')
            ax.set_ylabel('Illumination-intensity normalized average' + '\n' + 'amount of background photons per pixel')
            ax.grid()
            ax.set_xlim(left=-130, right = 130)
            ax.set_ylim(bottom = 0, top = 9)
            ax.tick_params(axis='both', which='major')
            ax.tick_params(axis='both', which='minor')   

        if panel == 6:
            CRLBx, musum, signalsum, Bisum = load_4('thetaI_thetab_triangulation_rot')
            thetab_list, thetaIrange = load_thetaI_thetab('thetaI_thetab_triangulation_rot')
            
            ax.loglog(thetaIrange,CRLBx.T*10**9)
            ax.set_xlabel('Expected signal photon budget ' + r'$\theta_I$' + ' (photons)')
            ax.set_ylabel(r'$x$'+'-Cramér-Rao lower bound (nm)')
            ax.grid()
            ax.set_xlim(left=100, right = 10000)
            ax.set_ylim(bottom = 0.1, top = 40)
            ax.tick_params(axis='both', which='major')
            ax.tick_params(axis='both', which='minor')  
            
            #Legend
            C0_line = mlines.Line2D([0], [0],color='C0', label=(r'$\theta_b = 0$'))
            C1_line = mlines.Line2D([0], [0],color='C1', label=(r'$\theta_b = 1$'))
            C2_line = mlines.Line2D([0], [0],color='C2', label=(r'$\theta_b = 4$'))
            C3_line = mlines.Line2D([0], [0],color='C3', label=(r'$\theta_b = 8$'))
            C8_line = mlines.Line2D([0], [0],color='C4', label=(r'$\theta_b = 16$'))
            C9_line = mlines.Line2D([0], [0],color='None', label=('(in phot/px.)'))
    
            ax.legend(handles=[C0_line, C1_line, C2_line, C3_line, C8_line, C9_line], ncol = 2, prop={'size':13}, labelspacing = 1)  

        if panel == 0:
            #Legend
            C0_line = mlines.Line2D([0], [0],color='C0', label=r'$r = 0.5\sigma_{\mathrm{PSF}}$')
            C1_line = mlines.Line2D([0], [0],color='C1', label=r'$r = 1\sigma_{\mathrm{PSF}}$')
            C2_line = mlines.Line2D([0], [0],color='C2', label=r'$r = 1.5\sigma_{\mathrm{PSF}}$')
            C3_line = mlines.Line2D([0], [0],color='C3', label=r'$r = 2\sigma_{\mathrm{PSF}}$')
            C9_line = mlines.Line2D([0], [0],color='k', label=r'SMLM (Uniform' + '\n' + 'illumination intensity)')
            C10_line = mlines.Line2D([0], [0],color='k', linestyle='--', label=r'Localization on ISM' + '\n' + 'reconstruction data')
            C11_line = mlines.Line2D([0], [0],color='k', linestyle=':', label=r'Localization on' + '\n' + 'Fourier reweighted ISM' + '\n' + 'reconstruction data')
    
            ax.legend(handles=[C0_line, C1_line, C2_line, C3_line, C9_line, C10_line, C11_line], prop={'size':13}, loc='upper left', bbox_to_anchor=(1.04,1), labelspacing = 1)

        plt.savefig(('./Figures/S26/'+str(int(panel))) + '.svg')
        
    if plotnumber == 'S27':
        CRLBx, musum, signalsum, Bisum = load_4('xp_r_triangulation')
        CRLBx_SMLM, musum_SMLM, signalsum_SMLM, Bisum_SMLM = load_4('SMLM')
        s_list, xrange = load_xp_r('xp_r_triangulation')
        
        fig, ax = plt.subplots(1,1)
        fig.set_size_inches(5, 5, forward=True)
        
        if panel == 1:
            ax.semilogy((xrange-theta[0])/(10**-9),CRLBx.T*10**9)
            ax.semilogy((xrange-theta[0])/(10**-9),CRLBx_SMLM*10**9*np.ones(np.size(xrange)), color='k')
            ax.semilogy((xrange-theta[0])/(10**-9),CRLBx_ISM*10**9*np.ones(np.size(xrange)), color='k', linestyle='--')
            ax.semilogy((xrange-theta[0])/(10**-9),CRLBx_FR*10**9*np.ones(np.size(xrange)), color='k', linestyle=':')
            ax.set_xlabel('Emitter-focus distance ' + r'$(x_f - \theta_x)$' + ' (nm)')
            ax.set_ylabel(r'$x$'+'-Cramér-Rao lower bound (nm)')
            ax.grid()
            ax.set_xlim(left=-130, right = 130)
            ax.set_ylim(bottom = 0.6, top = 60)
            ax.tick_params(axis='both', which='major')
            ax.tick_params(axis='both', which='minor')
            
        if panel == 2:
            ax.plot((xrange-theta[0])/(10**-9),CRLBx_SMLM/CRLBx.T)
            ax.plot((xrange-theta[0])/(10**-9),imp_ISM*np.ones(np.size(xrange)), color='k', linestyle='--')
            ax.plot((xrange-theta[0])/(10**-9),imp_FR*np.ones(np.size(xrange)), color='k', linestyle=':')
            ax.set_xlabel('Emitter-focus distance ' + r'$(x_f - \theta_x)$' + ' (nm)')
            ax.set_ylabel('Improvement over SMLM in ' + r'$x$'+ '-direction')
            ax.grid()
            ax.set_xlim(left=-130, right = 130)
            ax.set_ylim(bottom = 0, top = 2.2)
            ax.tick_params(axis='both', which='major')
            ax.tick_params(axis='both', which='minor')                     

        if panel == 3:
            ax.plot((xrange-theta[0])/(10**-9),musum.T)
            ax.plot((xrange-theta[0])/(10**-9),musum_SMLM*np.ones(np.size(xrange)), color='k')
            ax.set_xlabel('Emitter-focus distance ' + r'$(x_f - \theta_x)$' + ' (nm)')
            ax.set_ylabel('Expected amount of photons')
            ax.grid()
            ax.set_xlim(left=-130, right = 130)
            ax.set_ylim(bottom = 0, top = 2900)
            ax.tick_params(axis='both', which='major')
            ax.tick_params(axis='both', which='minor')                  

        if panel == 4:
            ax.plot((xrange-theta[0])/(10**-9),signalsum.T)
            ax.plot((xrange-theta[0])/(10**-9),signalsum_SMLM*np.ones(np.size(xrange)), color='k')
            ax.set_xlabel('Emitter-focus distance ' + r'$(x_f - \theta_x)$' + ' (nm)')
            ax.set_ylabel('Illumination-intensity normalized' + '\n' + 'amount of signal photons')
            ax.grid()
            ax.set_xlim(left=-130, right = 130)
            ax.set_ylim(bottom = 0, top = 2100)
            ax.tick_params(axis='both', which='major')
            ax.tick_params(axis='both', which='minor')   

        if panel == 5:
            ax.plot((xrange-theta[0])/(10**-9),Bisum.T/100)
            ax.plot((xrange-theta[0])/(10**-9),Bisum_SMLM/100*np.ones(np.size(xrange)), color='k')
            ax.set_xlabel('Emitter-focus distance ' + r'$(x_f - \theta_x)$' + ' (nm)')
            ax.set_ylabel('Illumination-intensity normalized average' + '\n' + 'amount of background photons per pixel')
            ax.grid()
            ax.set_xlim(left=-130, right = 130)
            ax.set_ylim(bottom = 0, top = 9)
            ax.tick_params(axis='both', which='major')
            ax.tick_params(axis='both', which='minor')   

        if panel == 6:
            CRLBx, musum, signalsum, Bisum = load_4('thetaI_thetab_triangulation')
            thetab_list, thetaIrange = load_thetaI_thetab('thetaI_thetab_triangulation')
            
            ax.loglog(thetaIrange,CRLBx.T*10**9)
            ax.set_xlabel('Expected signal photon budget ' + r'$\theta_I$' + ' (photons)')
            ax.set_ylabel(r'$x$'+'-Cramér-Rao lower bound (nm)')
            ax.grid()
            ax.set_xlim(left=100, right = 10000)
            ax.set_ylim(bottom = 0.6, top = 60)
            ax.tick_params(axis='both', which='major')
            ax.tick_params(axis='both', which='minor')  

            #Legend
            C0_line = mlines.Line2D([0], [0],color='C0', label=(r'$\theta_b = 0$'))
            C1_line = mlines.Line2D([0], [0],color='C1', label=(r'$\theta_b = 1$'))
            C2_line = mlines.Line2D([0], [0],color='C2', label=(r'$\theta_b = 4$'))
            C3_line = mlines.Line2D([0], [0],color='C3', label=(r'$\theta_b = 8$'))
            C8_line = mlines.Line2D([0], [0],color='C4', label=(r'$\theta_b = 16$'))
            C9_line = mlines.Line2D([0], [0],color='None', label=('(in phot/px.)'))
    
            ax.legend(handles=[C0_line, C1_line, C2_line, C3_line, C8_line, C9_line], ncol = 2, prop={'size':13}, labelspacing = 1)    

        if panel == 0:
            #Legend
            C0_line = mlines.Line2D([0], [0],color='C0', label=r'$r = 0.5\sigma_{\mathrm{PSF}}$')
            C1_line = mlines.Line2D([0], [0],color='C1', label=r'$r = 1\sigma_{\mathrm{PSF}}$')
            C2_line = mlines.Line2D([0], [0],color='C2', label=r'$r = 1.5\sigma_{\mathrm{PSF}}$')
            C3_line = mlines.Line2D([0], [0],color='C3', label=r'$r = 2\sigma_{\mathrm{PSF}}$')
            C9_line = mlines.Line2D([0], [0],color='k', label=r'SMLM (Uniform' + '\n' + 'illumination intensity)')
            C10_line = mlines.Line2D([0], [0],color='k', linestyle='--', label=r'Localization on ISM' + '\n' + 'reconstruction data')
            C11_line = mlines.Line2D([0], [0],color='k', linestyle=':', label=r'Localization on' + '\n' + 'Fourier reweighted ISM' + '\n' + 'reconstruction data')
    
            ax.legend(handles=[C0_line, C1_line, C2_line, C3_line, C9_line, C10_line, C11_line], prop={'size':13}, loc='upper left', bbox_to_anchor=(1.04,1), labelspacing = 1)

        plt.savefig(('./Figures/S27/'+str(int(panel))) + '.svg')

    if plotnumber == 'S28':
        CRLBx, musum, signalsum, Bisum = load_4('xp_r_triangulation_rot_nocenter_doughnut')
        CRLBx_SMLM, musum_SMLM, signalsum_SMLM, Bisum_SMLM = load_4('SMLM')
        s_list, xrange = load_xp_r('xp_r_triangulation_rot_nocenter_doughnut')
        
        fig, ax = plt.subplots(1,1)
        fig.set_size_inches(5, 5, forward=True)
        
        if panel == 1:
            ax.semilogy((xrange-theta[0])/(10**-9),CRLBx.T[:,0:4]*10**9)
            ax.semilogy((xrange-theta[0])/(10**-9),CRLBx_SMLM*10**9*np.ones(np.size(xrange)), color='k')
            ax.semilogy((xrange-theta[0])/(10**-9),CRLBx_ISM*10**9*np.ones(np.size(xrange)), color='k', linestyle='--')
            ax.semilogy((xrange-theta[0])/(10**-9),CRLBx_FR*10**9*np.ones(np.size(xrange)), color='k', linestyle=':')
            ax.set_xlabel('Emitter-focus distance ' + r'$(x_f - \theta_x)$' + ' (nm)')
            ax.set_ylabel(r'$x$'+'-Cramér-Rao lower bound (nm)')
            ax.grid()
            ax.set_xlim(left=-130, right = 130)
            ax.set_ylim(bottom = 0.6, top = 1000)
            ax.tick_params(axis='both', which='major')
            ax.tick_params(axis='both', which='minor')
            
        if panel == 2:
            ax.plot((xrange-theta[0])/(10**-9),CRLBx_SMLM/CRLBx.T[:,0:4])
            ax.plot((xrange-theta[0])/(10**-9),imp_ISM*np.ones(np.size(xrange)), color='k', linestyle='--')
            ax.plot((xrange-theta[0])/(10**-9),imp_FR*np.ones(np.size(xrange)), color='k', linestyle=':') 
            ax.set_xlabel('Emitter-focus distance ' + r'$(x_f - \theta_x)$' + ' (nm)')
            ax.set_ylabel('Improvement over SMLM in ' + r'$x$'+ '-direction')
            ax.grid()
            ax.set_xlim(left=-130, right = 130)
            ax.set_ylim(bottom = 0, top = 3.5)
            ax.tick_params(axis='both', which='major')
            ax.tick_params(axis='both', which='minor')                     

        if panel == 3:
            ax.plot((xrange-theta[0])/(10**-9),musum.T[:,0:4])
            ax.plot((xrange-theta[0])/(10**-9),musum_SMLM*np.ones(np.size(xrange)), color='k')
            ax.set_xlabel('Emitter-focus distance ' + r'$(x_f - \theta_x)$' + ' (nm)')
            ax.set_ylabel('Expected amount of photons')
            ax.grid()
            ax.set_xlim(left=-130, right = 130)
            ax.set_ylim(bottom = 0, top = 7100)
            ax.tick_params(axis='both', which='major')
            ax.tick_params(axis='both', which='minor')                  

        if panel == 4:
            ax.plot((xrange-theta[0])/(10**-9),signalsum.T[:,0:4])
            ax.plot((xrange-theta[0])/(10**-9),signalsum_SMLM*np.ones(np.size(xrange)), color='k')
            ax.set_xlabel('Emitter-focus distance ' + r'$(x_f - \theta_x)$' + ' (nm)')
            ax.set_ylabel('Illumination-intensity normalized' + '\n' + 'amount of signal photons')
            ax.grid()
            ax.set_xlim(left=-130, right = 130)
            ax.set_ylim(bottom = 0, top = 2100)
            ax.tick_params(axis='both', which='major')
            ax.tick_params(axis='both', which='minor')   

        if panel == 5:
            ax.plot((xrange-theta[0])/(10**-9),Bisum.T[:,0:4]/100)
            ax.plot((xrange-theta[0])/(10**-9),Bisum_SMLM/100*np.ones(np.size(xrange)), color='k')
            ax.set_xlabel('Emitter-focus distance ' + r'$(x_f - \theta_x)$' + ' (nm)')
            ax.set_ylabel('Illumination-intensity normalized average' + '\n' + 'amount of background photons per pixel')
            ax.grid()
            ax.set_xlim(left=-130, right = 130)
            ax.set_ylim(bottom = 0, top = 9)
            ax.tick_params(axis='both', which='major')
            ax.tick_params(axis='both', which='minor')

        if panel == 6:
            CRLBx, musum, signalsum, Bisum = load_4('thetaI_thetab_triangulation_rot_nocenter_doughnut')
            thetab_list, thetaIrange = load_thetaI_thetab('thetaI_thetab_triangulation_rot_nocenter_doughnut')
            
            ax.loglog(thetaIrange,CRLBx.T*10**9)
            ax.set_xlabel('Expected signal photon budget ' + r'$\theta_I$' + ' (photons)')
            ax.set_ylabel(r'$x$'+'-Cramér-Rao lower bound (nm)')
            ax.grid()
            ax.set_xlim(left=100, right = 10000)
            ax.set_ylim(bottom = 0.1, top = 40)
            ax.tick_params(axis='both', which='major')
            ax.tick_params(axis='both', which='minor')  
            
            #Legend
            C0_line = mlines.Line2D([0], [0],color='C0', label=(r'$\theta_b = 0$'))
            C1_line = mlines.Line2D([0], [0],color='C1', label=(r'$\theta_b = 1$'))
            C2_line = mlines.Line2D([0], [0],color='C2', label=(r'$\theta_b = 4$'))
            C3_line = mlines.Line2D([0], [0],color='C3', label=(r'$\theta_b = 8$'))
            C8_line = mlines.Line2D([0], [0],color='C4', label=(r'$\theta_b = 16$'))
            C9_line = mlines.Line2D([0], [0],color='None', label=('(in phot/px.)'))
    
            ax.legend(handles=[C0_line, C1_line, C2_line, C3_line, C8_line, C9_line], ncol = 2, prop={'size':13}, labelspacing = 1)  

        if panel == 0:
            #Legend
            C0_line = mlines.Line2D([0], [0],color='C0', label=r'$r = 1\sigma_{\mathrm{PSF}}$')
            C1_line = mlines.Line2D([0], [0],color='C1', label=r'$r = 2\sigma_{\mathrm{PSF}}$')
            C2_line = mlines.Line2D([0], [0],color='C2', label=r'$r = 3\sigma_{\mathrm{PSF}}$')
            C3_line = mlines.Line2D([0], [0],color='C3', label=r'$r = 4\sigma_{\mathrm{PSF}}$')
            C9_line = mlines.Line2D([0], [0],color='k', label=r'SMLM (Uniform' + '\n' + 'illumination intensity)')
            C10_line = mlines.Line2D([0], [0],color='k', linestyle='--', label=r'Localization on ISM' + '\n' + 'reconstruction data')
            C11_line = mlines.Line2D([0], [0],color='k', linestyle=':', label=r'Localization on' + '\n' + 'Fourier reweighted ISM' + '\n' + 'reconstruction data')
    
            ax.legend(handles=[C0_line, C1_line, C2_line, C3_line, C9_line, C10_line, C11_line], prop={'size':13}, loc='upper left', bbox_to_anchor=(1.04,1), labelspacing = 1)

        plt.savefig(('./Figures/S28/'+str(int(panel))) + '.svg')
        
    if plotnumber == 'S29':
        CRLBx, musum, signalsum, Bisum = load_4('xp_r_triangulation_nocenter_doughnut')
        CRLBx_SMLM, musum_SMLM, signalsum_SMLM, Bisum_SMLM = load_4('SMLM')
        s_list, xrange = load_xp_r('xp_r_triangulation_nocenter_doughnut')
        
        fig, ax = plt.subplots(1,1)
        fig.set_size_inches(5, 5, forward=True)
        
        if panel == 1:
            ax.semilogy((xrange-theta[0])/(10**-9),CRLBx.T[:,0:4]*10**9)
            ax.semilogy((xrange-theta[0])/(10**-9),CRLBx_SMLM*10**9*np.ones(np.size(xrange)), color='k')
            ax.semilogy((xrange-theta[0])/(10**-9),CRLBx_ISM*10**9*np.ones(np.size(xrange)), color='k', linestyle='--')
            ax.semilogy((xrange-theta[0])/(10**-9),CRLBx_FR*10**9*np.ones(np.size(xrange)), color='k', linestyle=':')
            ax.set_xlabel('Emitter-focus distance ' + r'$(x_f - \theta_x)$' + ' (nm)')
            ax.set_ylabel(r'$x$'+'-Cramér-Rao lower bound (nm)')
            ax.grid()
            ax.set_xlim(left=-130, right = 130)
            ax.set_ylim(bottom = 0.6, top = 1000)
            ax.tick_params(axis='both', which='major')
            ax.tick_params(axis='both', which='minor')
            
        if panel == 2:
            ax.plot((xrange-theta[0])/(10**-9),CRLBx_SMLM/CRLBx.T[:,0:4])
            ax.plot((xrange-theta[0])/(10**-9),imp_ISM*np.ones(np.size(xrange)), color='k', linestyle='--')
            ax.plot((xrange-theta[0])/(10**-9),imp_FR*np.ones(np.size(xrange)), color='k', linestyle=':')
            ax.set_xlabel('Emitter-focus distance ' + r'$(x_f - \theta_x)$' + ' (nm)')
            ax.set_ylabel('Improvement over SMLM in ' + r'$x$'+ '-direction')
            ax.grid()
            ax.set_xlim(left=-130, right = 130)
            ax.set_ylim(bottom = 0, top = 3.5)
            ax.tick_params(axis='both', which='major')
            ax.tick_params(axis='both', which='minor')                     

        if panel == 3:
            ax.plot((xrange-theta[0])/(10**-9),musum.T[:,0:4])
            ax.plot((xrange-theta[0])/(10**-9),musum_SMLM*np.ones(np.size(xrange)), color='k')
            ax.set_xlabel('Emitter-focus distance ' + r'$(x_f - \theta_x)$' + ' (nm)')
            ax.set_ylabel('Expected amount of photons')
            ax.grid()
            ax.set_xlim(left=-130, right = 130)
            ax.tick_params(axis='both', which='major')
            ax.tick_params(axis='both', which='minor')                  

        if panel == 4:
            ax.plot((xrange-theta[0])/(10**-9),signalsum.T[:,0:4])
            ax.plot((xrange-theta[0])/(10**-9),signalsum_SMLM*np.ones(np.size(xrange)), color='k')
            ax.set_xlabel('Emitter-focus distance ' + r'$(x_f - \theta_x)$' + ' (nm)')
            ax.set_ylabel('Illumination-intensity normalized' + '\n' + 'amount of signal photons')
            ax.grid()
            ax.set_xlim(left=-130, right = 130)
            ax.set_ylim(bottom = 0, top = 2100)
            ax.tick_params(axis='both', which='major')
            ax.tick_params(axis='both', which='minor')   

        if panel == 5:
            ax.plot((xrange-theta[0])/(10**-9),Bisum.T[:,0:4]/100)
            ax.plot((xrange-theta[0])/(10**-9),Bisum_SMLM/100*np.ones(np.size(xrange)), color='k')
            ax.set_xlabel('Emitter-focus distance ' + r'$(x_f - \theta_x)$' + ' (nm)')
            ax.set_ylabel('Illumination-intensity normalized average' + '\n' + 'amount of background photons per pixel')
            ax.grid()
            ax.set_xlim(left=-130, right = 130)
            ax.set_ylim(bottom = 0, top = 9)
            ax.tick_params(axis='both', which='major')
            ax.tick_params(axis='both', which='minor')  

        if panel == 6:
            CRLBx, musum, signalsum, Bisum = load_4('thetaI_thetab_triangulation_nocenter_doughnut')
            thetab_list, thetaIrange = load_thetaI_thetab('thetaI_thetab_triangulation_nocenter_doughnut')
            
            ax.loglog(thetaIrange,CRLBx.T*10**9)
            ax.set_xlabel('Expected signal photon budget ' + r'$\theta_I$' + ' (photons)')
            ax.set_ylabel(r'$x$'+'-Cramér-Rao lower bound (nm)')
            ax.grid()
            ax.set_xlim(left=100, right = 10000)
            ax.set_ylim(bottom = 0.1, top = 40)
            ax.tick_params(axis='both', which='major')
            ax.tick_params(axis='both', which='minor')  
            
            #Legend
            C0_line = mlines.Line2D([0], [0],color='C0', label=(r'$\theta_b = 0$'))
            C1_line = mlines.Line2D([0], [0],color='C1', label=(r'$\theta_b = 1$'))
            C2_line = mlines.Line2D([0], [0],color='C2', label=(r'$\theta_b = 4$'))
            C3_line = mlines.Line2D([0], [0],color='C3', label=(r'$\theta_b = 8$'))
            C8_line = mlines.Line2D([0], [0],color='C4', label=(r'$\theta_b = 16$'))
            C9_line = mlines.Line2D([0], [0],color='None', label=('(in phot/px.)'))
    
            ax.legend(handles=[C0_line, C1_line, C2_line, C3_line, C8_line, C9_line], ncol = 2, prop={'size':13}, labelspacing = 1)    

        if panel == 0:
            #Legend
            C0_line = mlines.Line2D([0], [0],color='C0', label=r'$r = 1\sigma_{\mathrm{PSF}}$')
            C1_line = mlines.Line2D([0], [0],color='C1', label=r'$r = 2\sigma_{\mathrm{PSF}}$')
            C2_line = mlines.Line2D([0], [0],color='C2', label=r'$r = 3\sigma_{\mathrm{PSF}}$')
            C3_line = mlines.Line2D([0], [0],color='C3', label=r'$r = 4\sigma_{\mathrm{PSF}}$')
            C9_line = mlines.Line2D([0], [0],color='k', label=r'SMLM (Uniform' + '\n' + 'illumination intensity)')
            C10_line = mlines.Line2D([0], [0],color='k', linestyle='--', label=r'Localization on ISM' + '\n' + 'reconstruction data')
            C11_line = mlines.Line2D([0], [0],color='k', linestyle=':', label=r'Localization on' + '\n' + 'Fourier reweighted ISM' + '\n' + 'reconstruction data')
    
            ax.legend(handles=[C0_line, C1_line, C2_line, C3_line, C9_line, C10_line, C11_line], prop={'size':13}, loc='upper left', bbox_to_anchor=(1.04,1), labelspacing = 1)

        plt.savefig(('./Figures/S29/'+str(int(panel))) + '.svg')

    if plotnumber == 'S30':
        CRLBx, musum, signalsum, Bisum = load_4('xp_r_triangulation_rot_doughnut_ext')
        CRLBx_SMLM, musum_SMLM, signalsum_SMLM, Bisum_SMLM = load_4('SMLM')
        s_list, xrange = load_xp_r('xp_r_triangulation_rot_doughnut_ext')
        
        fig, ax = plt.subplots(1,1)
        fig.set_size_inches(5, 5, forward=True)
        
        if panel == 1:
            ax.semilogy((xrange-theta[0])/(10**-9),CRLBx.T*10**9)
            ax.semilogy((xrange-theta[0])/(10**-9),CRLBx_SMLM*10**9*np.ones(np.size(xrange)), color='k')
            ax.semilogy((xrange-theta[0])/(10**-9),CRLBx_ISM*10**9*np.ones(np.size(xrange)), color='k', linestyle='--')
            ax.semilogy((xrange-theta[0])/(10**-9),CRLBx_FR*10**9*np.ones(np.size(xrange)), color='k', linestyle=':')
            ax.set_xlabel('Emitter-focus distance ' + r'$(x_f - \theta_x)$' + ' (nm)')
            ax.set_ylabel(r'$x$'+'-Cramér-Rao lower bound (nm)')
            ax.grid()
            ax.set_xlim(left=-130, right = 130)
            ax.set_ylim(bottom = 0.6, top = 60)
            ax.tick_params(axis='both', which='major')
            ax.tick_params(axis='both', which='minor')
            
        if panel == 2:
            ax.plot((xrange-theta[0])/(10**-9),CRLBx_SMLM/CRLBx.T)
            ax.plot((xrange-theta[0])/(10**-9),imp_ISM*np.ones(np.size(xrange)), color='k', linestyle='--')
            ax.plot((xrange-theta[0])/(10**-9),imp_FR*np.ones(np.size(xrange)), color='k', linestyle=':')
            ax.set_xlabel('Emitter-focus distance ' + r'$(x_f - \theta_x)$' + ' (nm)')
            ax.set_ylabel('Improvement over SMLM in ' + r'$x$'+ '-direction')
            ax.grid()
            ax.set_xlim(left=-130, right = 130)
            ax.set_ylim(bottom = 0, top = 3.5)
            ax.tick_params(axis='both', which='major')
            ax.tick_params(axis='both', which='minor')                     

        if panel == 3:
            ax.plot((xrange-theta[0])/(10**-9),musum.T)
            ax.plot((xrange-theta[0])/(10**-9),musum_SMLM*np.ones(np.size(xrange)), color='k')
            ax.set_xlabel('Emitter-focus distance ' + r'$(x_f - \theta_x)$' + ' (nm)')
            ax.set_ylabel('Expected amount of photons')
            ax.grid()
            ax.set_xlim(left=-130, right = 130)
            ax.set_ylim(bottom = 0, top = 2900)
            ax.tick_params(axis='both', which='major')
            ax.tick_params(axis='both', which='minor')                  

        if panel == 4:
            ax.plot((xrange-theta[0])/(10**-9),signalsum.T)
            ax.plot((xrange-theta[0])/(10**-9),signalsum_SMLM*np.ones(np.size(xrange)), color='k')
            ax.set_xlabel('Emitter-focus distance ' + r'$(x_f - \theta_x)$' + ' (nm)')
            ax.set_ylabel('Illumination-intensity normalized' + '\n' + 'amount of signal photons')
            ax.grid()
            ax.set_xlim(left=-130, right = 130)
            ax.set_ylim(bottom = 0, top = 2100)
            ax.tick_params(axis='both', which='major')
            ax.tick_params(axis='both', which='minor')   

        if panel == 5:
            ax.plot((xrange-theta[0])/(10**-9),Bisum.T/100)
            ax.plot((xrange-theta[0])/(10**-9),Bisum_SMLM/100*np.ones(np.size(xrange)), color='k')
            ax.set_xlabel('Emitter-focus distance ' + r'$(x_f - \theta_x)$' + ' (nm)')
            ax.set_ylabel('Illumination-intensity normalized average' + '\n' + 'amount of background photons per pixel')
            ax.grid()
            ax.set_xlim(left=-130, right = 130)
            ax.set_ylim(bottom = 0, top = 9)
            ax.tick_params(axis='both', which='major')
            ax.tick_params(axis='both', which='minor')   

        if panel == 6:
            CRLBx, musum, signalsum, Bisum = load_4('thetaI_thetab_triangulation_rot_doughnut')
            thetab_list, thetaIrange = load_thetaI_thetab('thetaI_thetab_triangulation_rot_doughnut')
            
            ax.loglog(thetaIrange,CRLBx.T*10**9)
            ax.set_xlabel('Expected signal photon budget ' + r'$\theta_I$' + ' (photons)')
            ax.set_ylabel(r'$x$'+'-Cramér-Rao lower bound (nm)')
            ax.grid()
            ax.set_xlim(left=100, right = 10000)
            ax.set_ylim(bottom = 0.1, top = 40)
            ax.tick_params(axis='both', which='major')
            ax.tick_params(axis='both', which='minor')  
            
            #Legend
            C0_line = mlines.Line2D([0], [0],color='C0', label=(r'$\theta_b = 0$'))
            C1_line = mlines.Line2D([0], [0],color='C1', label=(r'$\theta_b = 1$'))
            C2_line = mlines.Line2D([0], [0],color='C2', label=(r'$\theta_b = 4$'))
            C3_line = mlines.Line2D([0], [0],color='C3', label=(r'$\theta_b = 8$'))
            C8_line = mlines.Line2D([0], [0],color='C4', label=(r'$\theta_b = 16$'))
            C9_line = mlines.Line2D([0], [0],color='None', label=('(in phot/px.)'))
    
            ax.legend(handles=[C0_line, C1_line, C2_line, C3_line, C8_line, C9_line], ncol = 2, prop={'size':13}, labelspacing = 1) 

        if panel == 0:
            #Legend
            C0_line = mlines.Line2D([0], [0],color='C0', label=r'$r = 0.5\sigma_{\mathrm{PSF}}$')
            C1_line = mlines.Line2D([0], [0],color='C1', label=r'$r = 1\sigma_{\mathrm{PSF}}$')
            C2_line = mlines.Line2D([0], [0],color='C2', label=r'$r = 2\sigma_{\mathrm{PSF}}$')
            C3_line = mlines.Line2D([0], [0],color='C3', label=r'$r = 3\sigma_{\mathrm{PSF}}$')
            C4_line = mlines.Line2D([0], [0],color='C4', label=r'$r = 4\sigma_{\mathrm{PSF}}$')
            C9_line = mlines.Line2D([0], [0],color='k', label=r'SMLM (Uniform' + '\n' + 'illumination intensity)')
            C10_line = mlines.Line2D([0], [0],color='k', linestyle='--', label=r'Localization on ISM' + '\n' + 'reconstruction data')
            C11_line = mlines.Line2D([0], [0],color='k', linestyle=':', label=r'Localization on' + '\n' + 'Fourier reweighted ISM' + '\n' + 'reconstruction data')
    
            ax.legend(handles=[C0_line, C1_line, C2_line, C3_line, C4_line, C9_line, C10_line, C11_line], prop={'size':13}, loc='upper left', bbox_to_anchor=(1.04,1), labelspacing = 1)

        plt.savefig(('./Figures/S30/'+str(int(panel))) + '.svg')

    if plotnumber == 'S31':
        CRLBx, musum, signalsum, Bisum = load_4('xp_r_triangulation_doughnut_ext')
        CRLBx_SMLM, musum_SMLM, signalsum_SMLM, Bisum_SMLM = load_4('SMLM')
        s_list, xrange = load_xp_r('xp_r_triangulation_doughnut_ext')
        
        fig, ax = plt.subplots(1,1)
        fig.set_size_inches(5, 5, forward=True)
        
        if panel == 1:
            ax.semilogy((xrange-theta[0])/(10**-9),CRLBx.T*10**9)
            ax.semilogy((xrange-theta[0])/(10**-9),CRLBx_SMLM*10**9*np.ones(np.size(xrange)), color='k')
            ax.semilogy((xrange-theta[0])/(10**-9),CRLBx_ISM*10**9*np.ones(np.size(xrange)), color='k', linestyle='--')
            ax.semilogy((xrange-theta[0])/(10**-9),CRLBx_FR*10**9*np.ones(np.size(xrange)), color='k', linestyle=':')
            ax.set_xlabel('Emitter-focus distance ' + r'$(x_f - \theta_x)$' + ' (nm)')
            ax.set_ylabel(r'$x$'+'-Cramér-Rao lower bound (nm)')
            ax.grid()
            ax.set_xlim(left=-130, right = 130)
            ax.set_ylim(bottom = 0.6, top = 60)
            ax.tick_params(axis='both', which='major')
            ax.tick_params(axis='both', which='minor')
            
        if panel == 2:
            ax.plot((xrange-theta[0])/(10**-9),CRLBx_SMLM/CRLBx.T)
            ax.plot((xrange-theta[0])/(10**-9),imp_ISM*np.ones(np.size(xrange)), color='k', linestyle='--')
            ax.plot((xrange-theta[0])/(10**-9),imp_FR*np.ones(np.size(xrange)), color='k', linestyle=':')
            ax.set_xlabel('Emitter-focus distance ' + r'$(x_f - \theta_x)$' + ' (nm)')
            ax.set_ylabel('Improvement over SMLM in ' + r'$x$'+ '-direction')
            ax.grid()
            ax.set_xlim(left=-130, right = 130)
            ax.set_ylim(bottom = 0, top = 3.5)
            ax.tick_params(axis='both', which='major')
            ax.tick_params(axis='both', which='minor')                     

        if panel == 3:
            ax.plot((xrange-theta[0])/(10**-9),musum.T)
            ax.plot((xrange-theta[0])/(10**-9),musum_SMLM*np.ones(np.size(xrange)), color='k')
            ax.set_xlabel('Emitter-focus distance ' + r'$(x_f - \theta_x)$' + ' (nm)')
            ax.set_ylabel('Expected amount of photons')
            ax.grid()
            ax.set_xlim(left=-130, right = 130)
            ax.set_ylim(bottom = 0, top = 2900)
            ax.tick_params(axis='both', which='major')
            ax.tick_params(axis='both', which='minor')                  

        if panel == 4:
            ax.plot((xrange-theta[0])/(10**-9),signalsum.T)
            ax.plot((xrange-theta[0])/(10**-9),signalsum_SMLM*np.ones(np.size(xrange)), color='k')
            ax.set_xlabel('Emitter-focus distance ' + r'$(x_f - \theta_x)$' + ' (nm)')
            ax.set_ylabel('Illumination-intensity normalized' + '\n' + 'amount of signal photons')
            ax.grid()
            ax.set_xlim(left=-130, right = 130)
            ax.set_ylim(bottom = 0, top = 2100)
            ax.tick_params(axis='both', which='major')
            ax.tick_params(axis='both', which='minor')   

        if panel == 5:
            ax.plot((xrange-theta[0])/(10**-9),Bisum.T/100)
            ax.plot((xrange-theta[0])/(10**-9),Bisum_SMLM/100*np.ones(np.size(xrange)), color='k')
            ax.set_xlabel('Emitter-focus distance ' + r'$(x_f - \theta_x)$' + ' (nm)')
            ax.set_ylabel('Illumination-intensity normalized average' + '\n' + 'amount of background photons per pixel')
            ax.grid()
            ax.set_xlim(left=-130, right = 130)
            ax.set_ylim(bottom = 0, top = 9)
            ax.tick_params(axis='both', which='major')
            ax.tick_params(axis='both', which='minor')   

        if panel == 6:
            CRLBx, musum, signalsum, Bisum = load_4('thetaI_thetab_triangulation_doughnut')
            thetab_list, thetaIrange = load_thetaI_thetab('thetaI_thetab_triangulation_doughnut')
            
            ax.loglog(thetaIrange,CRLBx.T*10**9)
            ax.set_xlabel('Expected signal photon budget ' + r'$\theta_I$' + ' (photons)')
            ax.set_ylabel(r'$x$'+'-Cramér-Rao lower bound (nm)')
            ax.grid()
            ax.set_xlim(left=100, right = 10000)
            ax.set_ylim(bottom = 0.1, top = 40)
            ax.tick_params(axis='both', which='major')
            ax.tick_params(axis='both', which='minor')  
            
            #Legend
            C0_line = mlines.Line2D([0], [0],color='C0', label=(r'$\theta_b = 0$'))
            C1_line = mlines.Line2D([0], [0],color='C1', label=(r'$\theta_b = 1$'))
            C2_line = mlines.Line2D([0], [0],color='C2', label=(r'$\theta_b = 4$'))
            C3_line = mlines.Line2D([0], [0],color='C3', label=(r'$\theta_b = 8$'))
            C8_line = mlines.Line2D([0], [0],color='C4', label=(r'$\theta_b = 16$'))
            C9_line = mlines.Line2D([0], [0],color='None', label=('(in phot/px.)'))
    
            ax.legend(handles=[C0_line, C1_line, C2_line, C3_line, C8_line, C9_line], ncol = 2, prop={'size':13}, labelspacing = 1)    
            
        if panel == 0:
            #Legend
            C0_line = mlines.Line2D([0], [0],color='C0', label=r'$r = 0.5\sigma_{\mathrm{PSF}}$')
            C1_line = mlines.Line2D([0], [0],color='C1', label=r'$r = 1\sigma_{\mathrm{PSF}}$')
            C2_line = mlines.Line2D([0], [0],color='C2', label=r'$r = 2\sigma_{\mathrm{PSF}}$')
            C3_line = mlines.Line2D([0], [0],color='C3', label=r'$r = 3\sigma_{\mathrm{PSF}}$')
            C4_line = mlines.Line2D([0], [0],color='C4', label=r'$r = 4\sigma_{\mathrm{PSF}}$')
            C9_line = mlines.Line2D([0], [0],color='k', label=r'SMLM (Uniform' + '\n' + 'illumination intensity)')
            C10_line = mlines.Line2D([0], [0],color='k', linestyle='--', label=r'Localization on ISM' + '\n' + 'reconstruction data')
            C11_line = mlines.Line2D([0], [0],color='k', linestyle=':', label=r'Localization on' + '\n' + 'Fourier reweighted ISM' + '\n' + 'reconstruction data')
    
            ax.legend(handles=[C0_line, C1_line, C2_line, C3_line, C4_line, C9_line, C10_line, C11_line], prop={'size':13}, loc='upper left', bbox_to_anchor=(1.04,1), labelspacing = 1)

        plt.savefig(('./Figures/S31/'+str(int(panel))) + '.svg')

    if plotnumber == 'S32':
        CRLBx, musum, signalsum, Bisum = load_4('xp_rp_1pattern_adaptedbg')
        CRLBx_SMLM, musum_SMLM, signalsum_SMLM, Bisum_SMLM = load_4('SMLM')
        r_p_list, xrange = load_xp_rp('xp_rp_1pattern_adaptedbg')
        
        fig, ax = plt.subplots(1,1)
        fig.set_size_inches(5, 5, forward=True)
        
        if panel == 1:
            ax.semilogy((xrange-theta[0])/(10**-9),CRLBx.T[:,0:4]*10**9)
            ax.semilogy((xrange-theta[0])/(10**-9),CRLBx.T[:,-1]*10**9)
            ax.semilogy((xrange-theta[0])/(10**-9),CRLBx_SMLM*10**9*np.ones(np.size(xrange)), color='k')
            ax.semilogy((xrange-theta[0])/(10**-9),CRLBx_ISM*10**9*np.ones(np.size(xrange)), color='k', linestyle='--')
            ax.semilogy((xrange-theta[0])/(10**-9),CRLBx_FR*10**9*np.ones(np.size(xrange)), color='k', linestyle=':')
            ax.set_xlabel('Emitter-pinhole distance ' + r'$(x_p - \theta_x)$' + ' (nm)')
            ax.set_ylabel(r'$x$'+'-Cramér-Rao lower bound (nm)')
            ax.grid()
            ax.set_xlim(left=-130, right = 130)
            ax.set_ylim(bottom = 0.6, top = 60)
            ax.tick_params(axis='both', which='major')
            ax.tick_params(axis='both', which='minor')
            
        if panel == 2:
            ax.plot((xrange-theta[0])/(10**-9),CRLBx_SMLM/CRLBx.T[:,0:4])
            ax.plot((xrange-theta[0])/(10**-9),CRLBx_SMLM/CRLBx.T[:,-1])
            ax.plot((xrange-theta[0])/(10**-9),imp_ISM*np.ones(np.size(xrange)), color='k', linestyle='--')
            ax.plot((xrange-theta[0])/(10**-9),imp_FR*np.ones(np.size(xrange)), color='k', linestyle=':')
            ax.set_xlabel('Emitter-pinhole distance ' + r'$(x_p - \theta_x)$' + ' (nm)')
            ax.set_ylabel('Improvement over SMLM in ' + r'$x$'+ '-direction')
            ax.grid()
            ax.set_xlim(left=-130, right = 130)
            ax.set_ylim(bottom = 0, top = 3.5)
            ax.tick_params(axis='both', which='major')
            ax.tick_params(axis='both', which='minor')                     

        if panel == 3:
            ax.plot((xrange-theta[0])/(10**-9),musum.T[:,0:4])
            ax.plot((xrange-theta[0])/(10**-9),musum.T[:,-1])
            ax.plot((xrange-theta[0])/(10**-9),musum_SMLM*np.ones(np.size(xrange)), color='k')
            ax.set_xlabel('Emitter-pinhole distance ' + r'$(x_p - \theta_x)$' + ' (nm)')
            ax.set_ylabel('Expected amount of photons')
            ax.grid()
            ax.set_xlim(left=-130, right = 130)
            ax.set_ylim(bottom = 0, top = 2900)
            ax.tick_params(axis='both', which='major')
            ax.tick_params(axis='both', which='minor')                  

        if panel == 4:
            ax.plot((xrange-theta[0])/(10**-9),signalsum.T[:,0:4])
            ax.plot((xrange-theta[0])/(10**-9),signalsum.T[:,-1])
            ax.plot((xrange-theta[0])/(10**-9),signalsum_SMLM*np.ones(np.size(xrange)), color='k')
            ax.set_xlabel('Emitter-pinhole distance ' + r'$(x_p - \theta_x)$' + ' (nm)')
            ax.set_ylabel('Illumination-intensity normalized' + '\n' + 'amount of signal photons')
            ax.grid()
            ax.set_xlim(left=-130, right = 130)
            ax.set_ylim(bottom = 0, top = 2100)
            ax.tick_params(axis='both', which='major')
            ax.tick_params(axis='both', which='minor')   

        if panel == 5:
            ax.plot((xrange-theta[0])/(10**-9),Bisum.T[:,0:4]/100)
            ax.plot((xrange-theta[0])/(10**-9),Bisum.T[:,-1]/100)
            ax.plot((xrange-theta[0])/(10**-9),Bisum_SMLM/100*np.ones(np.size(xrange)), color='k')
            ax.set_xlabel('Emitter-pinhole distance ' + r'$(x_p - \theta_x)$' + ' (nm)')
            ax.set_ylabel('Average amount of background' + '\n' + 'photons per pixel')
            ax.grid()
            ax.set_xlim(left=-130, right = 130)
            ax.set_ylim(bottom = 0, top = 9)
            ax.tick_params(axis='both', which='major')
            ax.tick_params(axis='both', which='minor')   

        if panel == 6:
            CRLBx, musum, signalsum, Bisum = load_4('thetaI_thetab_1pattern_adaptedbg')
            thetab_list, thetaIrange = load_thetaI_thetab('thetaI_thetab_1pattern_adaptedbg')
            
            ax.loglog(thetaIrange,CRLBx.T*10**9)
            ax.set_xlabel('Expected signal photon budget ' + r'$\theta_I$' + ' (photons)')
            ax.set_ylabel(r'$x$'+'-Cramér-Rao lower bound (nm)')
            ax.grid()
            ax.set_xlim(left=100, right = 10000)
            ax.set_ylim(bottom = 0.3, top = 40)
            ax.tick_params(axis='both', which='major')
            ax.tick_params(axis='both', which='minor')  
            
            #Legend
            C0_line = mlines.Line2D([0], [0],color='C0', label=(r'$\theta_b = 0$'))
            C1_line = mlines.Line2D([0], [0],color='C1', label=(r'$\theta_b = 1$'))
            C2_line = mlines.Line2D([0], [0],color='C2', label=(r'$\theta_b = 4$'))
            C3_line = mlines.Line2D([0], [0],color='C3', label=(r'$\theta_b = 8$'))
            C8_line = mlines.Line2D([0], [0],color='C4', label=(r'$\theta_b = 16$'))
            C9_line = mlines.Line2D([0], [0],color='None', label=('(in phot/px.)'))
    
            ax.legend(handles=[C0_line, C1_line, C2_line, C3_line, C8_line, C9_line], ncol = 2, prop={'size':13}, labelspacing = 1)    
            
        if panel == 0:
            #Legend
            C0_line = mlines.Line2D([0], [0],color='C0', label=r'$r_p = 1\sigma_{\mathrm{PSF}}$')
            C1_line = mlines.Line2D([0], [0],color='C1', label=r'$r_p = 2\sigma_{\mathrm{PSF}}$')
            C2_line = mlines.Line2D([0], [0],color='C2', label=r'$r_p = 3\sigma_{\mathrm{PSF}}$')
            C3_line = mlines.Line2D([0], [0],color='C3', label=r'$r_p = 4\sigma_{\mathrm{PSF}}$')
            C8_line = mlines.Line2D([0], [0],color='C4', label=r'$r_p = \infty$')
            C9_line = mlines.Line2D([0], [0],color='k', label=r'SMLM (Uniform' + '\n' + 'illumination intensity)')
            C10_line = mlines.Line2D([0], [0],color='k', linestyle='--', label=r'Localization on ISM' + '\n' + 'reconstruction data')
            C11_line = mlines.Line2D([0], [0],color='k', linestyle=':', label=r'Localization on' + '\n' + 'Fourier reweighted ISM' + '\n' + 'reconstruction data')
    
            ax.legend(handles=[C0_line, C1_line, C2_line, C3_line, C8_line, C9_line, C10_line, C11_line], prop={'size':13}, loc='upper left', bbox_to_anchor=(1.04,1), labelspacing = 1)
        
        plt.savefig(('./Figures/S32/'+str(int(panel))) + '.svg')  

    if plotnumber == 'S33':
        CRLBx, musum, signalsum, Bisum = load_4('yp_rp_1pattern_adaptedbg')
        CRLBx_SMLM, musum_SMLM, signalsum_SMLM, Bisum_SMLM = load_4('SMLM')
        r_p_list, yrange = load_yp_rp('yp_rp_1pattern_adaptedbg')
        
        fig, ax = plt.subplots(1,1)
        fig.set_size_inches(5, 5, forward=True)
        
        if panel == 1:
            ax.semilogy((yrange-theta[1])/(10**-9),CRLBx.T[:,0:4]*10**9)
            ax.semilogy((yrange-theta[1])/(10**-9),CRLBx.T[:,-1]*10**9)
            ax.semilogy((yrange-theta[1])/(10**-9),CRLBx_SMLM*10**9*np.ones(np.size(yrange)), color='k')
            ax.semilogy((yrange-theta[1])/(10**-9),CRLBx_ISM*10**9*np.ones(np.size(yrange)), color='k', linestyle='--')
            ax.semilogy((yrange-theta[1])/(10**-9),CRLBx_FR*10**9*np.ones(np.size(yrange)), color='k', linestyle=':')
            ax.set_xlabel('Emitter-pinhole distance ' + r'$(y_p - \theta_y)$' + ' (nm)')
            ax.set_ylabel(r'$x$'+'-Cramér-Rao lower bound (nm)')
            ax.grid()
            ax.set_xlim(left=-130, right = 130)
            ax.set_ylim(bottom = 0.6, top = 60)
            ax.tick_params(axis='both', which='major')
            ax.tick_params(axis='both', which='minor')
            
        if panel == 2:
            ax.plot((yrange-theta[1])/(10**-9),CRLBx_SMLM/CRLBx.T[:,0:4])
            ax.plot((yrange-theta[1])/(10**-9),CRLBx_SMLM/CRLBx.T[:,-1])
            ax.plot((yrange-theta[1])/(10**-9),imp_ISM*np.ones(np.size(yrange)), color='k', linestyle='--')
            ax.plot((yrange-theta[1])/(10**-9),imp_FR*np.ones(np.size(yrange)), color='k', linestyle=':')
            ax.set_xlabel('Emitter-pinhole distance ' + r'$(y_p - \theta_y)$' + ' (nm)')
            ax.set_ylabel('Improvement over SMLM in ' + r'$x$'+ '-direction')
            ax.grid()
            ax.set_xlim(left=-130, right = 130)
            ax.set_ylim(bottom = 0, top = 3.5)
            ax.tick_params(axis='both', which='major')
            ax.tick_params(axis='both', which='minor')                     

        if panel == 3:
            ax.plot((yrange-theta[1])/(10**-9),musum.T[:,0:4])
            ax.plot((yrange-theta[1])/(10**-9),musum.T[:,-1])
            ax.plot((yrange-theta[1])/(10**-9),musum_SMLM*np.ones(np.size(yrange)), color='k')
            ax.set_xlabel('Emitter-pinhole distance ' + r'$(y_p - \theta_y)$' + ' (nm)')
            ax.set_ylabel('Expected amount of photons')
            ax.grid()
            ax.set_xlim(left=-130, right = 130)
            ax.set_ylim(bottom = 0, top = 2900)
            ax.tick_params(axis='both', which='major')
            ax.tick_params(axis='both', which='minor')                  

        if panel == 4:
            ax.plot((yrange-theta[1])/(10**-9),signalsum.T[:,0:4])
            ax.plot((yrange-theta[1])/(10**-9),signalsum.T[:,-1])
            ax.plot((yrange-theta[1])/(10**-9),signalsum_SMLM*np.ones(np.size(yrange)), color='k')
            ax.set_xlabel('Emitter-pinhole distance ' + r'$(y_p - \theta_y)$' + ' (nm)')
            ax.set_ylabel('Illumination-intensity normalized' + '\n' + 'amount of signal photons')
            ax.grid()
            ax.set_xlim(left=-130, right = 130)
            ax.set_ylim(bottom = 0, top = 2100)
            ax.tick_params(axis='both', which='major')
            ax.tick_params(axis='both', which='minor')   

        if panel == 5:
            ax.plot((yrange-theta[1])/(10**-9),Bisum.T[:,0:4]/100)
            ax.plot((yrange-theta[1])/(10**-9),Bisum.T[:,-1]/100)
            ax.plot((yrange-theta[1])/(10**-9),Bisum_SMLM*np.ones(np.size(yrange)), color='k')
            ax.set_xlabel('Emitter-pinhole distance ' + r'$(y_p - \theta_y)$' + ' (nm)')
            ax.set_ylabel('Average amount of background' + '\n' + 'photons per pixel')
            ax.grid()
            ax.set_xlim(left=-130, right = 130)
            ax.set_ylim(bottom = 0, top = 9)
            ax.tick_params(axis='both', which='major')
            ax.tick_params(axis='both', which='minor')  
        
        if panel == 0:
            #Legend
            C0_line = mlines.Line2D([0], [0],color='C0', label=r'$r_p = 1\sigma_{\mathrm{PSF}}$')
            C1_line = mlines.Line2D([0], [0],color='C1', label=r'$r_p = 2\sigma_{\mathrm{PSF}}$')
            C2_line = mlines.Line2D([0], [0],color='C2', label=r'$r_p = 3\sigma_{\mathrm{PSF}}$')
            C3_line = mlines.Line2D([0], [0],color='C3', label=r'$r_p = 4\sigma_{\mathrm{PSF}}$')
            C8_line = mlines.Line2D([0], [0],color='C4', label=r'$r_p = \infty$')
            C9_line = mlines.Line2D([0], [0],color='k', label=r'SMLM (Uniform' + '\n' + 'illumination intensity)')
            C10_line = mlines.Line2D([0], [0],color='k', linestyle='--', label=r'Localization on ISM' + '\n' + 'reconstruction data')
            C11_line = mlines.Line2D([0], [0],color='k', linestyle=':', label=r'Localization on' + '\n' + 'Fourier reweighted ISM' + '\n' + 'reconstruction data')
    
            ax.legend(handles=[C0_line, C1_line, C2_line, C3_line, C8_line, C9_line, C10_line, C11_line], prop={'size':13}, loc='upper left', bbox_to_anchor=(1.04,1), labelspacing = 1)

        plt.savefig(('./Figures/S33/'+str(int(panel))) + '.svg')
        
    if plotnumber == 'S34':
        CRLBx, musum, signalsum, Bisum = load_4('xp_s_2pattern_adaptedbg')
        CRLBx_SMLM, musum_SMLM, signalsum_SMLM, Bisum_SMLM = load_4('SMLM')
        s_list, xrange = load_xp_s('xp_s_2pattern_adaptedbg')
        
        fig, ax = plt.subplots(1,1)
        fig.set_size_inches(5, 5, forward=True)
        
        if panel == 1:
            ax.semilogy((xrange-theta[0])/(10**-9),CRLBx.T*10**9)
            ax.semilogy((xrange-theta[0])/(10**-9),CRLBx_SMLM*10**9*np.ones(np.size(xrange)), color='k')
            ax.semilogy((xrange-theta[0])/(10**-9),CRLBx_ISM*10**9*np.ones(np.size(xrange)), color='k', linestyle='--')
            ax.semilogy((xrange-theta[0])/(10**-9),CRLBx_FR*10**9*np.ones(np.size(xrange)), color='k', linestyle=':')
            ax.set_xlabel('Emitter-focus distance ' + r'$(x_f - \theta_x)$' + ' (nm)')
            ax.set_ylabel(r'$x$'+'-Cramér-Rao lower bound (nm)')
            ax.grid()
            ax.set_xlim(left=-130, right = 130)
            ax.set_ylim(bottom = 0.6, top = 60)
            ax.tick_params(axis='both', which='major')
            ax.tick_params(axis='both', which='minor')
            
        if panel == 2:
            ax.plot((xrange-theta[0])/(10**-9),CRLBx_SMLM/CRLBx.T)
            ax.plot((xrange-theta[0])/(10**-9),imp_ISM*np.ones(np.size(xrange)), color='k', linestyle='--')
            ax.plot((xrange-theta[0])/(10**-9),imp_FR*np.ones(np.size(xrange)), color='k', linestyle=':')
            ax.set_xlabel('Emitter-focus distance ' + r'$(x_f - \theta_x)$' + ' (nm)')
            ax.set_ylabel('Improvement over SMLM in ' + r'$x$'+ '-direction')
            ax.grid()
            ax.set_xlim(left=-130, right = 130)
            ax.set_ylim(bottom = 0, top = 3.5)
            ax.tick_params(axis='both', which='major')
            ax.tick_params(axis='both', which='minor')                     

        if panel == 3:
            ax.plot((xrange-theta[0])/(10**-9),musum.T)
            ax.plot((xrange-theta[0])/(10**-9),musum_SMLM*np.ones(np.size(xrange)), color='k')
            ax.set_xlabel('Emitter-focus distance ' + r'$(x_f - \theta_x)$' + ' (nm)')
            ax.set_ylabel('Expected amount of photons')
            ax.grid()
            ax.set_xlim(left=-130, right = 130)
            ax.set_ylim(bottom = 0, top = 2900)
            ax.tick_params(axis='both', which='major')
            ax.tick_params(axis='both', which='minor')                  

        if panel == 4:
            ax.plot((xrange-theta[0])/(10**-9),signalsum.T)
            ax.plot((xrange-theta[0])/(10**-9),signalsum_SMLM*np.ones(np.size(xrange)), color='k')
            ax.set_xlabel('Emitter-focus distance ' + r'$(x_f - \theta_x)$' + ' (nm)')
            ax.set_ylabel('Illumination-intensity normalized' + '\n' + 'amount of signal photons')
            ax.grid()
            ax.set_xlim(left=-130, right = 130)
            ax.set_ylim(bottom = 0, top = 2100)
            ax.tick_params(axis='both', which='major')
            ax.tick_params(axis='both', which='minor')   

        if panel == 5:
            ax.plot((xrange-theta[0])/(10**-9),Bisum.T/100)
            ax.plot((xrange-theta[0])/(10**-9),Bisum_SMLM/100*np.ones(np.size(xrange)), color='k')
            ax.set_xlabel('Emitter-focus distance ' + r'$(x_f - \theta_x)$' + ' (nm)')
            ax.set_ylabel('Average amount of background' + '\n' + 'photons per pixel')
            ax.grid()
            ax.set_xlim(left=-130, right = 130)
            ax.set_ylim(bottom = 0, top = 9)
            ax.tick_params(axis='both', which='major')
            ax.tick_params(axis='both', which='minor')   

        if panel == 6:
            CRLBx, musum, signalsum, Bisum = load_4('thetaI_thetab_2pattern_adaptedbg')

            thetab_list, thetaIrange = load_thetaI_thetab('thetaI_thetab_2pattern_adaptedbg')
            
            ax.loglog(thetaIrange,CRLBx.T*10**9)
            ax.set_xlabel('Expected signal photon budget ' + r'$\theta_I$' + ' (photons)')
            ax.set_ylabel(r'$x$'+'-Cramér-Rao lower bound (nm)')
            ax.grid()
            ax.set_xlim(left=100, right = 10000)
            ax.set_ylim(bottom = 0.3, top = 40)
            ax.tick_params(axis='both', which='major')
            ax.tick_params(axis='both', which='minor')  
            
            #Legend
            C0_line = mlines.Line2D([0], [0],color='C0', label=(r'$\theta_b = 0$'))
            C1_line = mlines.Line2D([0], [0],color='C1', label=(r'$\theta_b = 1$'))
            C2_line = mlines.Line2D([0], [0],color='C2', label=(r'$\theta_b = 4$'))
            C3_line = mlines.Line2D([0], [0],color='C3', label=(r'$\theta_b = 8$'))
            C8_line = mlines.Line2D([0], [0],color='C4', label=(r'$\theta_b = 16$'))
            C9_line = mlines.Line2D([0], [0],color='None', label=('(in phot/px.)'))
    
            ax.legend(handles=[C0_line, C1_line, C2_line, C3_line, C8_line, C9_line], ncol = 2, prop={'size':13}, labelspacing = 1)    
            
        if panel == 0:
            #Legend
            C0_line = mlines.Line2D([0], [0],color='C0', label=r'$s = 0\sigma_{\mathrm{PSF}}$')
            C1_line = mlines.Line2D([0], [0],color='C1', label=r'$s= 1\sigma_{\mathrm{PSF}}$')
            C2_line = mlines.Line2D([0], [0],color='C2', label=r'$s= 2\sigma_{\mathrm{PSF}}$')
            C3_line = mlines.Line2D([0], [0],color='C3', label=r'$s= 3\sigma_{\mathrm{PSF}}$')
            C4_line = mlines.Line2D([0], [0],color='C4', label=r'$s= 4\sigma_{\mathrm{PSF}}$')
            C8_line = mlines.Line2D([0], [0],color='C5', label=r'$s = 5\sigma_{\mathrm{PSF}}$')
            C9_line = mlines.Line2D([0], [0],color='k', label=r'SMLM (Uniform' + '\n' + 'illumination intensity)')
            C10_line = mlines.Line2D([0], [0],color='k', linestyle='--', label=r'Localization on ISM' + '\n' + 'reconstruction data')
            C11_line = mlines.Line2D([0], [0],color='k', linestyle=':', label=r'Localization on' + '\n' + 'Fourier reweighted ISM' + '\n' + 'reconstruction data')
    
            ax.legend(handles=[C0_line, C1_line, C2_line, C3_line, C4_line, C8_line, C9_line, C10_line, C11_line], prop={'size':13}, loc='upper left', bbox_to_anchor=(1.04,1), labelspacing = 1)

        plt.savefig(('./Figures/S34/'+str(int(panel))) + '.svg')

    if plotnumber == 'S35':
        CRLBx, musum, signalsum, Bisum = load_4('yp_xs_2pattern_adaptedbg')
        CRLBx_SMLM, musum_SMLM, signalsum_SMLM, Bisum_SMLM = load_4('SMLM')
        s_list, xrange = load_xp_s('yp_xs_2pattern_adaptedbg')
        
        fig, ax = plt.subplots(1,1)
        fig.set_size_inches(5, 5, forward=True)
        
        if panel == 1:
            ax.semilogy((xrange-theta[0])/(10**-9),CRLBx.T*10**9)
            ax.semilogy((xrange-theta[0])/(10**-9),CRLBx_SMLM*10**9*np.ones(np.size(xrange)), color='k')
            ax.semilogy((xrange-theta[0])/(10**-9),CRLBx_ISM*10**9*np.ones(np.size(xrange)), color='k', linestyle='--')
            ax.semilogy((xrange-theta[0])/(10**-9),CRLBx_FR*10**9*np.ones(np.size(xrange)), color='k', linestyle=':')
            ax.set_xlabel('Emitter-focus distance ' + r'$(y_f - \theta_y)$' + ' (nm)')
            ax.set_ylabel(r'$x$'+'-Cramér-Rao lower bound (nm)')
            ax.grid()
            ax.set_xlim(left=-130, right = 130)
            ax.set_ylim(bottom = 0.6, top = 60)
            ax.tick_params(axis='both', which='major')
            ax.tick_params(axis='both', which='minor')
            
        if panel == 2:
            ax.plot((xrange-theta[0])/(10**-9),CRLBx_SMLM/CRLBx.T)
            ax.plot((xrange-theta[0])/(10**-9),imp_ISM*np.ones(np.size(xrange)), color='k', linestyle='--')
            ax.plot((xrange-theta[0])/(10**-9),imp_FR*np.ones(np.size(xrange)), color='k', linestyle=':')
            ax.set_xlabel('Emitter-focus distance ' + r'$(y_f - \theta_y)$' + ' (nm)')
            ax.set_ylabel('Improvement over SMLM in ' + r'$x$'+ '-direction')
            ax.grid()
            ax.set_xlim(left=-130, right = 130)
            ax.set_ylim(bottom = 0, top = 3.5)
            ax.tick_params(axis='both', which='major')
            ax.tick_params(axis='both', which='minor')                     

        if panel == 3:
            ax.plot((xrange-theta[0])/(10**-9),musum.T)
            ax.plot((xrange-theta[0])/(10**-9),musum_SMLM*np.ones(np.size(xrange)), color='k')
            ax.set_xlabel('Emitter-focus distance ' + r'$(y_f - \theta_y)$' + ' (nm)')
            ax.set_ylabel('Expected amount of photons')
            ax.grid()
            ax.set_xlim(left=-130, right = 130)
            ax.set_ylim(bottom = 0, top = 2900)
            ax.tick_params(axis='both', which='major')
            ax.tick_params(axis='both', which='minor')                  

        if panel == 4:
            ax.plot((xrange-theta[0])/(10**-9),signalsum.T)
            ax.plot((xrange-theta[0])/(10**-9),signalsum_SMLM*np.ones(np.size(xrange)), color='k')
            ax.set_xlabel('Emitter-focus distance ' + r'$(y_f - \theta_y)$' + ' (nm)')
            ax.set_ylabel('Illumination-intensity normalized' + '\n' + 'amount of signal photons')
            ax.grid()
            ax.set_xlim(left=-130, right = 130)
            ax.set_ylim(bottom = 0, top = 2100)
            ax.tick_params(axis='both', which='major')
            ax.tick_params(axis='both', which='minor')   

        if panel == 5:
            ax.plot((xrange-theta[0])/(10**-9),Bisum.T/100)
            ax.plot((xrange-theta[0])/(10**-9),Bisum_SMLM/100*np.ones(np.size(xrange)), color='k')
            ax.set_xlabel('Emitter-focus distance ' + r'$(y_f - \theta_y)$' + ' (nm)')
            ax.set_ylabel('Average amount of background' + '\n' + 'photons per pixel')
            ax.grid()
            ax.set_xlim(left=-130, right = 130)
            ax.set_ylim(bottom = 0, top = 9)
            ax.tick_params(axis='both', which='major')
            ax.tick_params(axis='both', which='minor')   

        if panel == 0:
            #Legend
            C0_line = mlines.Line2D([0], [0],color='C0', label=r'$s = 0\sigma_{\mathrm{PSF}}$')
            C1_line = mlines.Line2D([0], [0],color='C1', label=r'$s= 1\sigma_{\mathrm{PSF}}$')
            C2_line = mlines.Line2D([0], [0],color='C2', label=r'$s= 2\sigma_{\mathrm{PSF}}$')
            C3_line = mlines.Line2D([0], [0],color='C3', label=r'$s= 3\sigma_{\mathrm{PSF}}$')
            C4_line = mlines.Line2D([0], [0],color='C4', label=r'$s= 4\sigma_{\mathrm{PSF}}$')
            C8_line = mlines.Line2D([0], [0],color='C5', label=r'$s = 5\sigma_{\mathrm{PSF}}$')
            C9_line = mlines.Line2D([0], [0],color='k', label=r'SMLM (Uniform' + '\n' + 'illumination intensity)')
            C10_line = mlines.Line2D([0], [0],color='k', linestyle='--', label=r'Localization on ISM' + '\n' + 'reconstruction data')
            C11_line = mlines.Line2D([0], [0],color='k', linestyle=':', label=r'Localization on' + '\n' + 'Fourier reweighted ISM' + '\n' + 'reconstruction data')
    
            ax.legend(handles=[C0_line, C1_line, C2_line, C3_line, C4_line, C8_line, C9_line, C10_line, C11_line], prop={'size':13}, loc='upper left', bbox_to_anchor=(1.04,1), labelspacing = 1)

        plt.savefig(('./Figures/S35/'+str(int(panel))) + '.svg')

    if plotnumber == 'S36':
        CRLBx, musum, signalsum, Bisum = load_4('xp_s_2pattern_nopinhole_adaptedbg')
        CRLBx_SMLM, musum_SMLM, signalsum_SMLM, Bisum_SMLM = load_4('SMLM')
        s_list, xrange = load_xp_s('xp_s_2pattern_nopinhole_adaptedbg')
        
        fig, ax = plt.subplots(1,1)
        fig.set_size_inches(5, 5, forward=True)
        
        if panel == 1:
            ax.semilogy((xrange-theta[0])/(10**-9),CRLBx.T*10**9)
            ax.semilogy((xrange-theta[0])/(10**-9),CRLBx_SMLM*10**9*np.ones(np.size(xrange)), color='k')
            ax.semilogy((xrange-theta[0])/(10**-9),CRLBx_ISM*10**9*np.ones(np.size(xrange)), color='k', linestyle='--')
            ax.semilogy((xrange-theta[0])/(10**-9),CRLBx_FR*10**9*np.ones(np.size(xrange)), color='k', linestyle=':')
            ax.set_xlabel('Emitter-focus distance ' + r'$(x_f - \theta_x)$' + ' (nm)')
            ax.set_ylabel(r'$x$'+'-Cramér-Rao lower bound (nm)')
            ax.grid()
            ax.set_xlim(left=-130, right = 130)
            ax.set_ylim(bottom = 0.6, top = 60)
            ax.tick_params(axis='both', which='major')
            ax.tick_params(axis='both', which='minor')
            
        if panel == 2:
            ax.plot((xrange-theta[0])/(10**-9),CRLBx_SMLM/CRLBx.T)
            ax.plot((xrange-theta[0])/(10**-9),imp_ISM*np.ones(np.size(xrange)), color='k', linestyle='--')
            ax.plot((xrange-theta[0])/(10**-9),imp_FR*np.ones(np.size(xrange)), color='k', linestyle=':')
            ax.set_xlabel('Emitter-focus distance ' + r'$(x_f - \theta_x)$' + ' (nm)')
            ax.set_ylabel('Improvement over SMLM in ' + r'$x$'+ '-direction')
            ax.grid()
            ax.set_xlim(left=-130, right = 130)
            ax.set_ylim(bottom = 0, top = 3.6)
            ax.tick_params(axis='both', which='major')
            ax.tick_params(axis='both', which='minor')                     

        if panel == 3:
            ax.plot((xrange-theta[0])/(10**-9),musum.T)
            ax.plot((xrange-theta[0])/(10**-9),musum_SMLM*np.ones(np.size(xrange)), color='k')
            ax.set_xlabel('Emitter-focus distance ' + r'$(x_f - \theta_x)$' + ' (nm)')
            ax.set_ylabel('Expected amount of photons')
            ax.grid()
            ax.set_xlim(left=-130, right = 130)
            ax.set_ylim(bottom = 0, top = 2900)
            ax.tick_params(axis='both', which='major')
            ax.tick_params(axis='both', which='minor')                  

        if panel == 4:
            ax.plot((xrange-theta[0])/(10**-9),signalsum.T)
            ax.plot((xrange-theta[0])/(10**-9),signalsum_SMLM*np.ones(np.size(xrange)), color='k')
            ax.set_xlabel('Emitter-focus distance ' + r'$(x_f - \theta_x)$' + ' (nm)')
            ax.set_ylabel('Illumination-intensity normalized' + '\n' + 'amount of signal photons')
            ax.grid()
            ax.set_xlim(left=-130, right = 130)
            ax.set_ylim(bottom = 0, top = 2100)
            ax.tick_params(axis='both', which='major')
            ax.tick_params(axis='both', which='minor')   

        if panel == 5:
            ax.plot((xrange-theta[0])/(10**-9),Bisum.T/100)
            ax.plot((xrange-theta[0])/(10**-9),Bisum_SMLM/100*np.ones(np.size(xrange)), color='k')
            ax.set_xlabel('Emitter-focus distance ' + r'$(x_f - \theta_x)$' + ' (nm)')
            ax.set_ylabel('Average amount of background' + '\n' + 'photons per pixel')
            ax.grid()
            ax.set_xlim(left=-130, right = 130)
            ax.set_ylim(bottom = 0, top = 9)
            ax.tick_params(axis='both', which='major')
            ax.tick_params(axis='both', which='minor')   
            
        if panel == 6:
            CRLBx, musum, signalsum, Bisum = load_4('thetaI_thetab_2pattern_nopinhole_adaptedbg')
            thetab_list, thetaIrange = load_thetaI_thetab('thetaI_thetab_2pattern_nopinhole_adaptedbg')
            
            ax.loglog(thetaIrange,CRLBx.T*10**9)
            ax.set_xlabel('Expected signal photon budget ' + r'$\theta_I$' + ' (photons)')
            ax.set_ylabel(r'$x$'+'-Cramér-Rao lower bound (nm)')
            ax.grid()
            ax.set_xlim(left=100, right = 10000)
            ax.set_ylim(bottom = 0.3, top = 40)
            ax.tick_params(axis='both', which='major')
            ax.tick_params(axis='both', which='minor')  
            
            #Legend
            C0_line = mlines.Line2D([0], [0],color='C0', label=(r'$\theta_b = 0$'))
            C1_line = mlines.Line2D([0], [0],color='C1', label=(r'$\theta_b = 1$'))
            C2_line = mlines.Line2D([0], [0],color='C2', label=(r'$\theta_b = 4$'))
            C3_line = mlines.Line2D([0], [0],color='C3', label=(r'$\theta_b = 8$'))
            C8_line = mlines.Line2D([0], [0],color='C4', label=(r'$\theta_b = 16$'))
            C9_line = mlines.Line2D([0], [0],color='None', label=('(in phot/px.)'))
    
            ax.legend(handles=[C0_line, C1_line, C2_line, C3_line, C8_line, C9_line], ncol = 2, prop={'size':13}, labelspacing = 1)    
        
        if panel == 0:
            #Legend
            C0_line = mlines.Line2D([0], [0],color='C0', label=r'$s = 0\sigma_{\mathrm{PSF}}$')
            C1_line = mlines.Line2D([0], [0],color='C1', label=r'$s= 1\sigma_{\mathrm{PSF}}$')
            C2_line = mlines.Line2D([0], [0],color='C2', label=r'$s= 2\sigma_{\mathrm{PSF}}$')
            C3_line = mlines.Line2D([0], [0],color='C3', label=r'$s= 3\sigma_{\mathrm{PSF}}$')
            C4_line = mlines.Line2D([0], [0],color='C4', label=r'$s= 4\sigma_{\mathrm{PSF}}$')
            C8_line = mlines.Line2D([0], [0],color='C5', label=r'$s = 5\sigma_{\mathrm{PSF}}$')
            C9_line = mlines.Line2D([0], [0],color='k', label=r'SMLM (Uniform' + '\n' + 'illumination intensity)')
            C10_line = mlines.Line2D([0], [0],color='k', linestyle='--', label=r'Localization on ISM' + '\n' + 'reconstruction data')
            C11_line = mlines.Line2D([0], [0],color='k', linestyle=':', label=r'Localization on' + '\n' + 'Fourier reweighted ISM' + '\n' + 'reconstruction data')
    
            ax.legend(handles=[C0_line, C1_line, C2_line, C3_line, C4_line, C8_line, C9_line, C10_line, C11_line], prop={'size':13}, loc='upper left', bbox_to_anchor=(1.04,1), labelspacing = 1)

        plt.savefig(('./Figures/S36/'+str(int(panel))) + '.svg')
    
    if plotnumber == 'S37':
        CRLBx, musum, signalsum, Bisum = load_4('yp_s_2pattern_adaptedbg')
        CRLBx_SMLM, musum_SMLM, signalsum_SMLM, Bisum_SMLM = load_4('SMLM')
        s_list, xrange = load_xp_s('yp_s_2pattern_adaptedbg')
        
        fig, ax = plt.subplots(1,1)
        fig.set_size_inches(5, 5, forward=True)
        
        if panel == 1:
            ax.semilogy((xrange-theta[0])/(10**-9),CRLBx.T*10**9)
            ax.semilogy((xrange-theta[0])/(10**-9),CRLBx_SMLM*10**9*np.ones(np.size(xrange)), color='k')
            ax.semilogy((xrange-theta[0])/(10**-9),CRLBx_ISM*10**9*np.ones(np.size(xrange)), color='k', linestyle='--')
            ax.semilogy((xrange-theta[0])/(10**-9),CRLBx_FR*10**9*np.ones(np.size(xrange)), color='k', linestyle=':')
            ax.set_xlabel('Emitter-focus distance ' + r'$(y_f - \theta_y)$' + ' (nm)')
            ax.set_ylabel(r'$x$'+'-Cramér-Rao lower bound (nm)')
            ax.grid()
            ax.set_xlim(left=-130, right = 130)
            ax.set_ylim(bottom = 0.6, top = 60)
            ax.tick_params(axis='both', which='major')
            ax.tick_params(axis='both', which='minor')
            
        if panel == 2:
            ax.plot((xrange-theta[0])/(10**-9),CRLBx_SMLM/CRLBx.T)
            ax.plot((xrange-theta[0])/(10**-9),imp_ISM*np.ones(np.size(xrange)), color='k', linestyle='--')
            ax.plot((xrange-theta[0])/(10**-9),imp_FR*np.ones(np.size(xrange)), color='k', linestyle=':')
            ax.set_xlabel('Emitter-focus distance ' + r'$(y_f - \theta_y)$' + ' (nm)')
            ax.set_ylabel('Improvement over SMLM in ' + r'$x$'+ '-direction')
            ax.grid()
            ax.set_xlim(left=-130, right = 130)
            ax.set_ylim(bottom = 0, top = 2.2)
            ax.tick_params(axis='both', which='major')
            ax.tick_params(axis='both', which='minor')                     

        if panel == 3:
            ax.plot((xrange-theta[0])/(10**-9),musum.T)
            ax.plot((xrange-theta[0])/(10**-9),musum_SMLM*np.ones(np.size(xrange)), color='k')
            ax.set_xlabel('Emitter-focus distance ' + r'$(y_f - \theta_y)$' + ' (nm)')
            ax.set_ylabel('Expected amount of photons')
            ax.grid()
            ax.set_xlim(left=-130, right = 130)
            ax.set_ylim(bottom = 0, top = 2900)
            ax.tick_params(axis='both', which='major')
            ax.tick_params(axis='both', which='minor')                  

        if panel == 4:
            ax.plot((xrange-theta[0])/(10**-9),signalsum.T)
            ax.plot((xrange-theta[0])/(10**-9),signalsum_SMLM*np.ones(np.size(xrange)), color='k')
            ax.set_xlabel('Emitter-focus distance ' + r'$(y_f - \theta_y)$' + ' (nm)')
            ax.set_ylabel('Illumination-intensity normalized' + '\n' + 'amount of signal photons')
            ax.grid()
            ax.set_xlim(left=-130, right = 130)
            ax.set_ylim(bottom = 0, top = 2100)
            ax.tick_params(axis='both', which='major')
            ax.tick_params(axis='both', which='minor')   

        if panel == 5:
            ax.plot((xrange-theta[0])/(10**-9),Bisum.T/100)
            ax.plot((xrange-theta[0])/(10**-9),Bisum_SMLM/100*np.ones(np.size(xrange)), color='k')
            ax.set_xlabel('Emitter-focus distance ' + r'$(y_f - \theta_y)$' + ' (nm)')
            ax.set_ylabel('Average amount of background' + '\n' + 'photons per pixel')
            ax.grid()
            ax.set_xlim(left=-130, right = 130)
            ax.set_ylim(bottom = 0, top = 9)
            ax.tick_params(axis='both', which='major')
            ax.tick_params(axis='both', which='minor')
            
        if panel == 6:
            CRLBx, musum, signalsum, Bisum = load_4('thetaI_thetab_2pattern_rot_adaptedbg')
            thetab_list, thetaIrange = load_thetaI_thetab('thetaI_thetab_2pattern_rot_adaptedbg')
            
            ax.loglog(thetaIrange,CRLBx.T*10**9)
            ax.set_xlabel('Expected signal photon budget ' + r'$\theta_I$' + ' (photons)')
            ax.set_ylabel(r'$x$'+'-Cramér-Rao lower bound (nm)')
            ax.grid()
            ax.set_xlim(left=100, right = 10000)
            ax.set_ylim(bottom = 0.3, top = 40)
            ax.tick_params(axis='both', which='major')
            ax.tick_params(axis='both', which='minor')  
            
            #Legend
            C0_line = mlines.Line2D([0], [0],color='C0', label=(r'$\theta_b = 0$'))
            C1_line = mlines.Line2D([0], [0],color='C1', label=(r'$\theta_b = 1$'))
            C2_line = mlines.Line2D([0], [0],color='C2', label=(r'$\theta_b = 4$'))
            C3_line = mlines.Line2D([0], [0],color='C3', label=(r'$\theta_b = 8$'))
            C8_line = mlines.Line2D([0], [0],color='C4', label=(r'$\theta_b = 16$'))
            C9_line = mlines.Line2D([0], [0],color='None', label=('(in phot/px.)'))
    
            ax.legend(handles=[C0_line, C1_line, C2_line, C3_line, C8_line, C9_line], ncol = 2, prop={'size':13}, labelspacing = 1)    
        
        if panel == 0:
            #Legend
            C0_line = mlines.Line2D([0], [0],color='C0', label=r'$s = 0\sigma_{\mathrm{PSF}}$')
            C1_line = mlines.Line2D([0], [0],color='C1', label=r'$s= 1\sigma_{\mathrm{PSF}}$')
            C2_line = mlines.Line2D([0], [0],color='C2', label=r'$s= 2\sigma_{\mathrm{PSF}}$')
            C3_line = mlines.Line2D([0], [0],color='C3', label=r'$s= 3\sigma_{\mathrm{PSF}}$')
            C4_line = mlines.Line2D([0], [0],color='C4', label=r'$s= 4\sigma_{\mathrm{PSF}}$')
            C8_line = mlines.Line2D([0], [0],color='C5', label=r'$s = 5\sigma_{\mathrm{PSF}}$')
            C9_line = mlines.Line2D([0], [0],color='k', label=r'SMLM (Uniform' + '\n' + 'illumination intensity)')
            C10_line = mlines.Line2D([0], [0],color='k', linestyle='--', label=r'Localization on ISM' + '\n' + 'reconstruction data')
            C11_line = mlines.Line2D([0], [0],color='k', linestyle=':', label=r'Localization on' + '\n' + 'Fourier reweighted ISM' + '\n' + 'reconstruction data')
    
            ax.legend(handles=[C0_line, C1_line, C2_line, C3_line, C4_line, C8_line, C9_line, C10_line, C11_line], prop={'size':13}, loc='upper left', bbox_to_anchor=(1.04,1), labelspacing = 1)

        plt.savefig(('./Figures/S37/'+str(int(panel))) + '.svg')

    if plotnumber == 'S38':
        CRLBx, musum, signalsum, Bisum = load_4('xp_r_triangulation_rot_nocenter_ext_adaptedbg')
        CRLBx_SMLM, musum_SMLM, signalsum_SMLM, Bisum_SMLM = load_4('SMLM')
        s_list, xrange = load_xp_r('xp_r_triangulation_rot_nocenter_ext_adaptedbg')
        
        fig, ax = plt.subplots(1,1)
        fig.set_size_inches(5, 5, forward=True)
        
        if panel == 1:
            ax.semilogy((xrange-theta[0])/(10**-9),CRLBx.T[:,0:4]*10**9)
            ax.semilogy((xrange-theta[0])/(10**-9),CRLBx_SMLM*10**9*np.ones(np.size(xrange)), color='k')
            ax.semilogy((xrange-theta[0])/(10**-9),CRLBx_ISM*10**9*np.ones(np.size(xrange)), color='k', linestyle='--')
            ax.semilogy((xrange-theta[0])/(10**-9),CRLBx_FR*10**9*np.ones(np.size(xrange)), color='k', linestyle=':')
            ax.set_xlabel('Emitter-focus distance ' + r'$(x_f - \theta_x)$' + ' (nm)')
            ax.set_ylabel(r'$x$'+'-Cramér-Rao lower bound (nm)')
            
            ax.grid()
            ax.set_xlim(left=-130, right = 130)
            ax.set_ylim(bottom = 0.6, top = 60)
            ax.tick_params(axis='both', which='major')
            ax.tick_params(axis='both', which='minor')
            
        if panel == 2:
            ax.plot((xrange-theta[0])/(10**-9),CRLBx_SMLM/CRLBx.T[:,0:4])
            ax.plot((xrange-theta[0])/(10**-9),imp_ISM*np.ones(np.size(xrange)), color='k', linestyle='--')
            ax.plot((xrange-theta[0])/(10**-9),imp_FR*np.ones(np.size(xrange)), color='k', linestyle=':')
            ax.set_xlabel('Emitter-focus distance ' + r'$(x_f - \theta_x)$' + ' (nm)')
            ax.set_ylabel('Improvement over SMLM in ' + r'$x$'+ '-direction')
            ax.grid()
            ax.set_xlim(left=-130, right = 130)
            ax.set_ylim(bottom = 0, top = 3)
            ax.tick_params(axis='both', which='major')
            ax.tick_params(axis='both', which='minor')                     

        if panel == 3:
            ax.plot((xrange-theta[0])/(10**-9),musum.T[:,0:4])
            ax.plot((xrange-theta[0])/(10**-9),musum_SMLM*np.ones(np.size(xrange)), color='k')
            ax.set_xlabel('Emitter-focus distance ' + r'$(x_f - \theta_x)$' + ' (nm)')
            ax.set_ylabel('Expected amount of photons')
            ax.grid()
            ax.set_xlim(left=-130, right = 130)
            ax.set_ylim(bottom = 0, top = 2900)
            ax.tick_params(axis='both', which='major')
            ax.tick_params(axis='both', which='minor')                  

        if panel == 4:
            ax.plot((xrange-theta[0])/(10**-9),signalsum.T[:,0:4])
            ax.plot((xrange-theta[0])/(10**-9),signalsum_SMLM*np.ones(np.size(xrange)), color='k')
            ax.set_xlabel('Emitter-focus distance ' + r'$(x_f - \theta_x)$' + ' (nm)')
            ax.set_ylabel('Illumination-intensity normalized' + '\n' + 'amount of signal photons')
            ax.grid()
            ax.set_xlim(left=-130, right = 130)
            ax.set_ylim(bottom = 0, top = 2100)
            ax.tick_params(axis='both', which='major')
            ax.tick_params(axis='both', which='minor')   

        if panel == 5:
            ax.plot((xrange-theta[0])/(10**-9),Bisum.T[:,0:4]/100)
            ax.plot((xrange-theta[0])/(10**-9),Bisum_SMLM/100*np.ones(np.size(xrange)), color='k')
            ax.set_xlabel('Emitter-focus distance ' + r'$(x_f - \theta_x)$' + ' (nm)')
            ax.set_ylabel('Average amount of background' + '\n' + 'photons per pixel')
            ax.grid()
            ax.set_xlim(left=-130, right = 130)
            ax.set_ylim(bottom = 0, top = 9)
            ax.tick_params(axis='both', which='major')
            ax.tick_params(axis='both', which='minor')   

        if panel == 6:
            CRLBx, musum, signalsum, Bisum = load_4('thetaI_thetab_triangulation_rot_nocenter_adaptedbg')
            thetab_list, thetaIrange = load_thetaI_thetab('thetaI_thetab_triangulation_rot_nocenter_adaptedbg')
            
            ax.loglog(thetaIrange,CRLBx.T*10**9)
            ax.set_xlabel('Expected signal photon budget ' + r'$\theta_I$' + ' (photons)')
            ax.set_ylabel(r'$x$'+'-Cramér-Rao lower bound (nm)')
            ax.grid()
            ax.set_xlim(left=100, right = 10000)
            ax.set_ylim(bottom = 0.1, top = 40)
            ax.tick_params(axis='both', which='major')
            ax.tick_params(axis='both', which='minor')  
            
            #Legend
            C0_line = mlines.Line2D([0], [0],color='C0', label=(r'$\theta_b = 0$'))
            C1_line = mlines.Line2D([0], [0],color='C1', label=(r'$\theta_b = 1$'))
            C2_line = mlines.Line2D([0], [0],color='C2', label=(r'$\theta_b = 4$'))
            C3_line = mlines.Line2D([0], [0],color='C3', label=(r'$\theta_b = 8$'))
            C8_line = mlines.Line2D([0], [0],color='C4', label=(r'$\theta_b = 16$'))
            C9_line = mlines.Line2D([0], [0],color='None', label=('(in phot/px.)'))
    
            ax.legend(handles=[C0_line, C1_line, C2_line, C3_line, C8_line, C9_line], ncol = 2, prop={'size':13}, labelspacing = 1)  

        if panel == 0:
            #Legend
            C0_line = mlines.Line2D([0], [0],color='C0', label=r'$r = 1\sigma_{\mathrm{PSF}}$')
            C1_line = mlines.Line2D([0], [0],color='C1', label=r'$r = 2\sigma_{\mathrm{PSF}}$')
            C2_line = mlines.Line2D([0], [0],color='C2', label=r'$r = 3\sigma_{\mathrm{PSF}}$')
            C3_line = mlines.Line2D([0], [0],color='C3', label=r'$r = 4\sigma_{\mathrm{PSF}}$')
            C9_line = mlines.Line2D([0], [0],color='k', label=r'SMLM (Uniform' + '\n' + 'illumination intensity)')
            C10_line = mlines.Line2D([0], [0],color='k', linestyle='--', label=r'Localization on ISM' + '\n' + 'reconstruction data')
            C11_line = mlines.Line2D([0], [0],color='k', linestyle=':', label=r'Localization on' + '\n' + 'Fourier reweighted ISM' + '\n' + 'reconstruction data')
    
            ax.legend(handles=[C0_line, C1_line, C2_line, C3_line, C9_line, C10_line, C11_line], prop={'size':13}, loc='upper left', bbox_to_anchor=(1.04,1), labelspacing = 1)

        plt.savefig(('./Figures/S38/'+str(int(panel))) + '.svg')

    if plotnumber == 'S39':
        CRLBx, musum, signalsum, Bisum = load_4('xp_r_triangulation_nocenter_ext_adaptedbg')
        CRLBx_SMLM, musum_SMLM, signalsum_SMLM, Bisum_SMLM = load_4('SMLM')
        s_list, xrange = load_xp_r('xp_r_triangulation_nocenter_ext_adaptedbg')
        
        fig, ax = plt.subplots(1,1)
        fig.set_size_inches(5, 5, forward=True)
        
        if panel == 1:
            ax.semilogy((xrange-theta[0])/(10**-9),CRLBx.T[:,0:4]*10**9)
            ax.semilogy((xrange-theta[0])/(10**-9),CRLBx_SMLM*10**9*np.ones(np.size(xrange)), color='k')
            ax.semilogy((xrange-theta[0])/(10**-9),CRLBx_ISM*10**9*np.ones(np.size(xrange)), color='k', linestyle='--')
            ax.semilogy((xrange-theta[0])/(10**-9),CRLBx_FR*10**9*np.ones(np.size(xrange)), color='k', linestyle=':')
            ax.set_xlabel('Emitter-focus distance ' + r'$(x_f - \theta_x)$' + ' (nm)')
            ax.set_ylabel(r'$x$'+'-Cramér-Rao lower bound (nm)')
            ax.grid()
            ax.set_xlim(left=-130, right = 130)
            ax.set_ylim(bottom = 0.6, top = 60)
            ax.tick_params(axis='both', which='major')
            ax.tick_params(axis='both', which='minor')
            
        if panel == 2:
            ax.plot((xrange-theta[0])/(10**-9),CRLBx_SMLM/CRLBx.T[:,0:4])
            ax.plot((xrange-theta[0])/(10**-9),imp_ISM*np.ones(np.size(xrange)), color='k', linestyle='--')
            ax.plot((xrange-theta[0])/(10**-9),imp_FR*np.ones(np.size(xrange)), color='k', linestyle=':')
            ax.set_xlabel('Emitter-focus distance ' + r'$(x_f - \theta_x)$' + ' (nm)')
            ax.set_ylabel('Improvement over SMLM in ' + r'$x$'+ '-direction')
            ax.grid()
            ax.set_xlim(left=-130, right = 130)
            ax.set_ylim(bottom = 0, top = 3)
            ax.tick_params(axis='both', which='major')
            ax.tick_params(axis='both', which='minor')                     

        if panel == 3:
            ax.plot((xrange-theta[0])/(10**-9),musum.T[:,0:4])
            ax.plot((xrange-theta[0])/(10**-9),musum_SMLM*np.ones(np.size(xrange)), color='k')
            ax.set_xlabel('Emitter-focus distance ' + r'$(x_f - \theta_x)$' + ' (nm)')
            ax.set_ylabel('Expected amount of photons')
            ax.grid()
            ax.set_xlim(left=-130, right = 130)
            ax.set_ylim(bottom = 0, top = 2900)
            ax.tick_params(axis='both', which='major')
            ax.tick_params(axis='both', which='minor')                  

        if panel == 4:
            ax.plot((xrange-theta[0])/(10**-9),signalsum.T[:,0:4])
            ax.plot((xrange-theta[0])/(10**-9),signalsum_SMLM*np.ones(np.size(xrange)), color='k')
            ax.set_xlabel('Emitter-focus distance ' + r'$(x_f - \theta_x)$' + ' (nm)')
            ax.set_ylabel('Illumination-intensity normalized' + '\n' + 'amount of signal photons')
            ax.grid()
            ax.set_xlim(left=-130, right = 130)
            ax.set_ylim(bottom = 0, top = 2100)
            ax.tick_params(axis='both', which='major')
            ax.tick_params(axis='both', which='minor')   

        if panel == 5:
            ax.plot((xrange-theta[0])/(10**-9),Bisum.T[:,0:4]/100)
            ax.plot((xrange-theta[0])/(10**-9),Bisum_SMLM/100*np.ones(np.size(xrange)), color='k')
            ax.set_xlabel('Emitter-focus distance ' + r'$(x_f - \theta_x)$' + ' (nm)')
            ax.set_ylabel('Average amount of background' + '\n' + 'photons per pixel')
            ax.grid()
            ax.set_xlim(left=-130, right = 130)
            ax.set_ylim(bottom = 0, top = 9)
            ax.tick_params(axis='both', which='major')
            ax.tick_params(axis='both', which='minor')   

        if panel == 6:
            CRLBx, musum, signalsum, Bisum = load_4('thetaI_thetab_triangulation_nocenter_adaptedbg')
            thetab_list, thetaIrange = load_thetaI_thetab('thetaI_thetab_triangulation_nocenter_adaptedbg')
            
            ax.loglog(thetaIrange,CRLBx.T*10**9)
            ax.set_xlabel('Expected signal photon budget ' + r'$\theta_I$' + ' (photons)')
            ax.set_ylabel(r'$x$'+'-Cramér-Rao lower bound (nm)')
            ax.grid()
            ax.set_xlim(left=100, right = 10000)
            ax.set_ylim(bottom = 0.1, top = 40)
            ax.tick_params(axis='both', which='major')
            ax.tick_params(axis='both', which='minor')  
            
            #Legend
            C0_line = mlines.Line2D([0], [0],color='C0', label=(r'$\theta_b = 0$'))
            C1_line = mlines.Line2D([0], [0],color='C1', label=(r'$\theta_b = 1$'))
            C2_line = mlines.Line2D([0], [0],color='C2', label=(r'$\theta_b = 4$'))
            C3_line = mlines.Line2D([0], [0],color='C3', label=(r'$\theta_b = 8$'))
            C8_line = mlines.Line2D([0], [0],color='C4', label=(r'$\theta_b = 16$'))
            C9_line = mlines.Line2D([0], [0],color='None', label=('(in phot/px.)'))
    
            ax.legend(handles=[C0_line, C1_line, C2_line, C3_line, C8_line, C9_line], ncol = 2, prop={'size':13}, labelspacing = 1) 

        if panel == 0:
            #Legend
            C0_line = mlines.Line2D([0], [0],color='C0', label=r'$r = 1\sigma_{\mathrm{PSF}}$')
            C1_line = mlines.Line2D([0], [0],color='C1', label=r'$r = 2\sigma_{\mathrm{PSF}}$')
            C2_line = mlines.Line2D([0], [0],color='C2', label=r'$r = 3\sigma_{\mathrm{PSF}}$')
            C3_line = mlines.Line2D([0], [0],color='C3', label=r'$r = 4\sigma_{\mathrm{PSF}}$')
            C9_line = mlines.Line2D([0], [0],color='k', label=r'SMLM (Uniform' + '\n' + 'illumination intensity)')
            C10_line = mlines.Line2D([0], [0],color='k', linestyle='--', label=r'Localization on ISM' + '\n' + 'reconstruction data')
            C11_line = mlines.Line2D([0], [0],color='k', linestyle=':', label=r'Localization on' + '\n' + 'Fourier reweighted ISM' + '\n' + 'reconstruction data')
    
            ax.legend(handles=[C0_line, C1_line, C2_line, C3_line, C9_line, C10_line, C11_line], prop={'size':13}, loc='upper left', bbox_to_anchor=(1.04,1), labelspacing = 1)

        plt.savefig(('./Figures/S39/'+str(int(panel))) + '.svg')

    if plotnumber == 'S40':
        CRLBx, musum, signalsum, Bisum = load_4('xp_r_triangulation_rot_adaptedbg')
        CRLBx_SMLM, musum_SMLM, signalsum_SMLM, Bisum_SMLM = load_4('SMLM')
        s_list, xrange = load_xp_r('xp_r_triangulation_rot_adaptedbg')
        
        fig, ax = plt.subplots(1,1)
        fig.set_size_inches(5, 5, forward=True)
        
        if panel == 1:
            ax.semilogy((xrange-theta[0])/(10**-9),CRLBx.T*10**9)
            ax.semilogy((xrange-theta[0])/(10**-9),CRLBx_SMLM*10**9*np.ones(np.size(xrange)), color='k')
            ax.semilogy((xrange-theta[0])/(10**-9),CRLBx_ISM*10**9*np.ones(np.size(xrange)), color='k', linestyle='--')
            ax.semilogy((xrange-theta[0])/(10**-9),CRLBx_FR*10**9*np.ones(np.size(xrange)), color='k', linestyle=':')
            ax.set_xlabel('Emitter-focus distance ' + r'$(x_f - \theta_x)$' + ' (nm)')
            ax.set_ylabel(r'$x$'+'-Cramér-Rao lower bound (nm)')
            ax.grid()
            ax.set_xlim(left=-130, right = 130)
            ax.set_ylim(bottom = 0.6, top = 60)
            ax.tick_params(axis='both', which='major')
            ax.tick_params(axis='both', which='minor')
            
        if panel == 2:
            ax.plot((xrange-theta[0])/(10**-9),CRLBx_SMLM/CRLBx.T)
            ax.plot((xrange-theta[0])/(10**-9),imp_ISM*np.ones(np.size(xrange)), color='k', linestyle='--')
            ax.plot((xrange-theta[0])/(10**-9),imp_FR*np.ones(np.size(xrange)), color='k', linestyle=':')
            ax.set_xlabel('Emitter-focus distance ' + r'$(x_f - \theta_x)$' + ' (nm)')
            ax.set_ylabel('Improvement over SMLM in ' + r'$x$'+ '-direction')
            ax.grid()
            ax.set_xlim(left=-130, right = 130)
            ax.set_ylim(bottom = 0, top = 3.5)
            ax.tick_params(axis='both', which='major')
            ax.tick_params(axis='both', which='minor')                     

        if panel == 3:
            ax.plot((xrange-theta[0])/(10**-9),musum.T)
            ax.plot((xrange-theta[0])/(10**-9),musum_SMLM*np.ones(np.size(xrange)), color='k')
            ax.set_xlabel('Emitter-focus distance ' + r'$(x_f - \theta_x)$' + ' (nm)')
            ax.set_ylabel('Expected amount of photons')
            ax.grid()
            ax.set_xlim(left=-130, right = 130)
            ax.set_ylim(bottom = 0, top = 2900)
            ax.tick_params(axis='both', which='major')
            ax.tick_params(axis='both', which='minor')                  

        if panel == 4:
            ax.plot((xrange-theta[0])/(10**-9),signalsum.T)
            ax.plot((xrange-theta[0])/(10**-9),signalsum_SMLM*np.ones(np.size(xrange)), color='k')
            ax.set_xlabel('Emitter-focus distance ' + r'$(x_f - \theta_x)$' + ' (nm)')
            ax.set_ylabel('Illumination-intensity normalized' + '\n' + 'amount of signal photons')
            ax.grid()
            ax.set_xlim(left=-130, right = 130)
            ax.set_ylim(bottom = 0, top = 2100)
            ax.tick_params(axis='both', which='major')
            ax.tick_params(axis='both', which='minor')   

        if panel == 5:
            ax.plot((xrange-theta[0])/(10**-9),Bisum.T/100)
            ax.plot((xrange-theta[0])/(10**-9),Bisum_SMLM/100*np.ones(np.size(xrange)), color='k')
            ax.set_xlabel('Emitter-focus distance ' + r'$(x_f - \theta_x)$' + ' (nm)')
            ax.set_ylabel('Average amount of background' + '\n' + 'photons per pixel')
            ax.grid()
            ax.set_xlim(left=-130, right = 130)
            ax.set_ylim(bottom = 0, top = 9)
            ax.tick_params(axis='both', which='major')
            ax.tick_params(axis='both', which='minor')   

        if panel == 6:
            CRLBx, musum, signalsum, Bisum = load_4('thetaI_thetab_triangulation_rot_adaptedbg')
            thetab_list, thetaIrange = load_thetaI_thetab('thetaI_thetab_triangulation_rot_adaptedbg')
            
            ax.loglog(thetaIrange,CRLBx.T*10**9)
            ax.set_xlabel('Expected signal photon budget ' + r'$\theta_I$' + ' (photons)')
            ax.set_ylabel(r'$x$'+'-Cramér-Rao lower bound (nm)')
            ax.grid()
            ax.set_xlim(left=100, right = 10000)
            ax.set_ylim(bottom = 0.1, top = 40)
            ax.tick_params(axis='both', which='major')
            ax.tick_params(axis='both', which='minor')  
            
            #Legend
            C0_line = mlines.Line2D([0], [0],color='C0', label=(r'$\theta_b = 0$'))
            C1_line = mlines.Line2D([0], [0],color='C1', label=(r'$\theta_b = 1$'))
            C2_line = mlines.Line2D([0], [0],color='C2', label=(r'$\theta_b = 4$'))
            C3_line = mlines.Line2D([0], [0],color='C3', label=(r'$\theta_b = 8$'))
            C8_line = mlines.Line2D([0], [0],color='C4', label=(r'$\theta_b = 16$'))
            C9_line = mlines.Line2D([0], [0],color='None', label=('(in phot/px.)'))
    
            ax.legend(handles=[C0_line, C1_line, C2_line, C3_line, C8_line, C9_line], ncol = 2, prop={'size':13}, labelspacing = 1)  

        if panel == 0:
            #Legend
            C0_line = mlines.Line2D([0], [0],color='C0', label=r'$r = 0.5\sigma_{\mathrm{PSF}}$')
            C1_line = mlines.Line2D([0], [0],color='C1', label=r'$r = 1\sigma_{\mathrm{PSF}}$')
            C2_line = mlines.Line2D([0], [0],color='C2', label=r'$r = 1.5\sigma_{\mathrm{PSF}}$')
            C3_line = mlines.Line2D([0], [0],color='C3', label=r'$r = 2\sigma_{\mathrm{PSF}}$')
            C9_line = mlines.Line2D([0], [0],color='k', label=r'SMLM (Uniform' + '\n' + 'illumination intensity)')
            C10_line = mlines.Line2D([0], [0],color='k', linestyle='--', label=r'Localization on ISM' + '\n' + 'reconstruction data')
            C11_line = mlines.Line2D([0], [0],color='k', linestyle=':', label=r'Localization on' + '\n' + 'Fourier reweighted ISM' + '\n' + 'reconstruction data')
    
            ax.legend(handles=[C0_line, C1_line, C2_line, C3_line, C9_line, C10_line, C11_line], prop={'size':13}, loc='upper left', bbox_to_anchor=(1.04,1), labelspacing = 1)

        plt.savefig(('./Figures/S40/'+str(int(panel))) + '.svg')

    if plotnumber == 'S41':
        CRLBx, musum, signalsum, Bisum = load_4('xp_r_triangulation_adaptedbg')
        CRLBx_SMLM, musum_SMLM, signalsum_SMLM, Bisum_SMLM = load_4('SMLM')
        s_list, xrange = load_xp_r('xp_r_triangulation_adaptedbg')
        
        fig, ax = plt.subplots(1,1)
        fig.set_size_inches(5, 5, forward=True)
        
        if panel == 1:
            ax.semilogy((xrange-theta[0])/(10**-9),CRLBx.T*10**9)
            ax.semilogy((xrange-theta[0])/(10**-9),CRLBx_SMLM*10**9*np.ones(np.size(xrange)), color='k')
            ax.semilogy((xrange-theta[0])/(10**-9),CRLBx_ISM*10**9*np.ones(np.size(xrange)), color='k', linestyle='--')
            ax.semilogy((xrange-theta[0])/(10**-9),CRLBx_FR*10**9*np.ones(np.size(xrange)), color='k', linestyle=':')
            ax.set_xlabel('Emitter-focus distance ' + r'$(x_f - \theta_x)$' + ' (nm)')
            ax.set_ylabel(r'$x$'+'-Cramér-Rao lower bound (nm)')
            ax.grid()
            ax.set_xlim(left=-130, right = 130)
            ax.set_ylim(bottom = 0.6, top = 60)
            ax.tick_params(axis='both', which='major')
            ax.tick_params(axis='both', which='minor')
            
        if panel == 2:
            ax.plot((xrange-theta[0])/(10**-9),CRLBx_SMLM/CRLBx.T)
            ax.plot((xrange-theta[0])/(10**-9),imp_ISM*np.ones(np.size(xrange)), color='k', linestyle='--')
            ax.plot((xrange-theta[0])/(10**-9),imp_FR*np.ones(np.size(xrange)), color='k', linestyle=':')
            ax.set_xlabel('Emitter-focus distance ' + r'$(x_f - \theta_x)$' + ' (nm)')
            ax.set_ylabel('Improvement over SMLM in ' + r'$x$'+ '-direction')
            ax.grid()
            ax.set_xlim(left=-130, right = 130)
            ax.set_ylim(bottom = 0, top = 3.5)
            ax.tick_params(axis='both', which='major')
            ax.tick_params(axis='both', which='minor')                     

        if panel == 3:
            ax.plot((xrange-theta[0])/(10**-9),musum.T)
            ax.plot((xrange-theta[0])/(10**-9),musum_SMLM*np.ones(np.size(xrange)), color='k')
            ax.set_xlabel('Emitter-focus distance ' + r'$(x_f - \theta_x)$' + ' (nm)')
            ax.set_ylabel('Expected amount of photons')
            ax.grid()
            ax.set_xlim(left=-130, right = 130)
            ax.set_ylim(bottom = 0, top = 2900)
            ax.tick_params(axis='both', which='major')
            ax.tick_params(axis='both', which='minor')                  

        if panel == 4:
            ax.plot((xrange-theta[0])/(10**-9),signalsum.T)
            ax.plot((xrange-theta[0])/(10**-9),signalsum_SMLM*np.ones(np.size(xrange)), color='k')
            ax.set_xlabel('Emitter-focus distance ' + r'$(x_f - \theta_x)$' + ' (nm)')
            ax.set_ylabel('Illumination-intensity normalized' + '\n' + 'amount of signal photons')
            ax.grid()
            ax.set_xlim(left=-130, right = 130)
            ax.set_ylim(bottom = 0, top = 2100)
            ax.tick_params(axis='both', which='major')
            ax.tick_params(axis='both', which='minor')   

        if panel == 5:
            ax.plot((xrange-theta[0])/(10**-9),Bisum.T/100)
            ax.plot((xrange-theta[0])/(10**-9),Bisum_SMLM/100*np.ones(np.size(xrange)), color='k')
            ax.set_xlabel('Emitter-focus distance ' + r'$(x_f - \theta_x)$' + ' (nm)')
            ax.set_ylabel('Average amount of background' + '\n' + 'photons per pixel')
            ax.grid()
            ax.set_xlim(left=-130, right = 130)
            ax.set_ylim(bottom = 0, top = 9)
            ax.tick_params(axis='both', which='major')
            ax.tick_params(axis='both', which='minor')   

        if panel == 6:
            CRLBx, musum, signalsum, Bisum = load_4('thetaI_thetab_triangulation_adaptedbg')
            thetab_list, thetaIrange = load_thetaI_thetab('thetaI_thetab_triangulation_adaptedbg')
            
            ax.loglog(thetaIrange,CRLBx.T*10**9)
            ax.set_xlabel('Expected signal photon budget ' + r'$\theta_I$' + ' (photons)')
            ax.set_ylabel(r'$x$'+'-Cramér-Rao lower bound (nm)')
            ax.grid()
            ax.set_xlim(left=100, right = 10000)
            ax.set_ylim(bottom = 0.3, top = 40)
            ax.tick_params(axis='both', which='major')
            ax.tick_params(axis='both', which='minor')  
            
            #Legend
            C0_line = mlines.Line2D([0], [0],color='C0', label=(r'$\theta_b = 0$'))
            C1_line = mlines.Line2D([0], [0],color='C1', label=(r'$\theta_b = 1$'))
            C2_line = mlines.Line2D([0], [0],color='C2', label=(r'$\theta_b = 4$'))
            C3_line = mlines.Line2D([0], [0],color='C3', label=(r'$\theta_b = 8$'))
            C8_line = mlines.Line2D([0], [0],color='C4', label=(r'$\theta_b = 16$'))
            C9_line = mlines.Line2D([0], [0],color='None', label=('(in phot/px.)'))
    
            ax.legend(handles=[C0_line, C1_line, C2_line, C3_line, C8_line, C9_line], ncol = 2, prop={'size':13}, labelspacing = 1)    
            
        if panel == 0:
            #Legend
            C0_line = mlines.Line2D([0], [0],color='C0', label=r'$r = 0.5\sigma_{\mathrm{PSF}}$')
            C1_line = mlines.Line2D([0], [0],color='C1', label=r'$r = 1\sigma_{\mathrm{PSF}}$')
            C2_line = mlines.Line2D([0], [0],color='C2', label=r'$r = 1.5\sigma_{\mathrm{PSF}}$')
            C3_line = mlines.Line2D([0], [0],color='C3', label=r'$r = 2\sigma_{\mathrm{PSF}}$')
            C9_line = mlines.Line2D([0], [0],color='k', label=r'SMLM (Uniform' + '\n' + 'illumination intensity)')
            C10_line = mlines.Line2D([0], [0],color='k', linestyle='--', label=r'Localization on ISM' + '\n' + 'reconstruction data')
            C11_line = mlines.Line2D([0], [0],color='k', linestyle=':', label=r'Localization on' + '\n' + 'Fourier reweighted ISM' + '\n' + 'reconstruction data')
    
            ax.legend(handles=[C0_line, C1_line, C2_line, C3_line, C9_line, C10_line, C11_line], prop={'size':13}, loc='upper left', bbox_to_anchor=(1.04,1), labelspacing = 1)

        plt.savefig(('./Figures/S41/'+str(int(panel))) + '.svg')
               
    if plotnumber == 'S42':
        CRLBx, musum, signalsum, Bisum = load_4('xp_r_triangulation_rot_nocenter_doughnut_adaptedbg')
        CRLBx_SMLM, musum_SMLM, signalsum_SMLM, Bisum_SMLM = load_4('SMLM')
        s_list, xrange = load_xp_r('xp_r_triangulation_rot_nocenter_doughnut_adaptedbg')
        
        fig, ax = plt.subplots(1,1)
        fig.set_size_inches(5, 5, forward=True)
        
        if panel == 1:
            ax.semilogy((xrange-theta[0])/(10**-9),CRLBx.T[:,0:4]*10**9)
            ax.semilogy((xrange-theta[0])/(10**-9),CRLBx_SMLM*10**9*np.ones(np.size(xrange)), color='k')
            ax.semilogy((xrange-theta[0])/(10**-9),CRLBx_ISM*10**9*np.ones(np.size(xrange)), color='k', linestyle='--')
            ax.semilogy((xrange-theta[0])/(10**-9),CRLBx_FR*10**9*np.ones(np.size(xrange)), color='k', linestyle=':')
            ax.set_xlabel('Emitter-focus distance ' + r'$(x_f - \theta_x)$' + ' (nm)')
            ax.set_ylabel(r'$x$'+'-Cramér-Rao lower bound (nm)')
            ax.grid()
            ax.set_xlim(left=-130, right = 130)
            ax.set_ylim(bottom = 0.6, top = 60)
            ax.tick_params(axis='both', which='major')
            ax.tick_params(axis='both', which='minor')
            
        if panel == 2:
            ax.plot((xrange-theta[0])/(10**-9),CRLBx_SMLM/CRLBx.T[:,0:4])
            ax.plot((xrange-theta[0])/(10**-9),imp_ISM*np.ones(np.size(xrange)), color='k', linestyle='--')
            ax.plot((xrange-theta[0])/(10**-9),imp_FR*np.ones(np.size(xrange)), color='k', linestyle=':')
            ax.set_xlabel('Emitter-focus distance ' + r'$(x_f - \theta_x)$' + ' (nm)')
            ax.set_ylabel('Improvement over SMLM in ' + r'$x$'+ '-direction')
            ax.grid()
            ax.set_xlim(left=-130, right = 130)
            ax.set_ylim(bottom = 0, top = 3.5)
            ax.tick_params(axis='both', which='major')
            ax.tick_params(axis='both', which='minor')                     

        if panel == 3:
            ax.plot((xrange-theta[0])/(10**-9),musum.T[:,0:4])
            ax.plot((xrange-theta[0])/(10**-9),musum_SMLM*np.ones(np.size(xrange)), color='k')
            ax.set_xlabel('Emitter-focus distance ' + r'$(x_f - \theta_x)$' + ' (nm)')
            ax.set_ylabel('Expected amount of photons')
            ax.grid()
            ax.set_xlim(left=-130, right = 130)
            ax.set_ylim(bottom = 0, top = 7100)
            ax.tick_params(axis='both', which='major')
            ax.tick_params(axis='both', which='minor')                  

        if panel == 4:
            ax.plot((xrange-theta[0])/(10**-9),signalsum.T[:,0:4])
            ax.plot((xrange-theta[0])/(10**-9),signalsum_SMLM*np.ones(np.size(xrange)), color='k')
            ax.set_xlabel('Emitter-focus distance ' + r'$(x_f - \theta_x)$' + ' (nm)')
            ax.set_ylabel('Illumination-intensity normalized' + '\n' + 'amount of signal photons')
            ax.grid()
            ax.set_xlim(left=-130, right = 130)
            ax.set_ylim(bottom = 0, top = 2100)
            ax.tick_params(axis='both', which='major')
            ax.tick_params(axis='both', which='minor')   

        if panel == 5:
            ax.plot((xrange-theta[0])/(10**-9),Bisum.T[:,0:4]/100)
            ax.plot((xrange-theta[0])/(10**-9),Bisum_SMLM/100*np.ones(np.size(xrange)), color='k')
            ax.set_xlabel('Emitter-focus distance ' + r'$(x_f - \theta_x)$' + ' (nm)')
            ax.set_ylabel('Illumination-intensity normalized average' + '\n' + 'amount of background photons per pixel')
            ax.grid()
            ax.set_xlim(left=-130, right = 130)
            ax.set_ylim(bottom = 0, top = 9)
            ax.tick_params(axis='both', which='major')
            ax.tick_params(axis='both', which='minor')

        if panel == 6:
            CRLBx, musum, signalsum, Bisum = load_4('thetaI_thetab_triangulation_rot_nocenter_doughnut_adaptedbg')
            thetab_list, thetaIrange = load_thetaI_thetab('thetaI_thetab_triangulation_rot_nocenter_doughnut_adaptedbg')
            
            ax.loglog(thetaIrange,CRLBx.T*10**9)
            ax.set_xlabel('Expected signal photon budget ' + r'$\theta_I$' + ' (photons)')
            ax.set_ylabel(r'$x$'+'-Cramér-Rao lower bound (nm)')
            ax.grid()
            ax.set_xlim(left=100, right = 10000)
            ax.set_ylim(bottom = 0.1, top = 40)
            ax.tick_params(axis='both', which='major')
            ax.tick_params(axis='both', which='minor')  
            
            #Legend
            C0_line = mlines.Line2D([0], [0],color='C0', label=(r'$\theta_b = 0$'))
            C1_line = mlines.Line2D([0], [0],color='C1', label=(r'$\theta_b = 1$'))
            C2_line = mlines.Line2D([0], [0],color='C2', label=(r'$\theta_b = 4$'))
            C3_line = mlines.Line2D([0], [0],color='C3', label=(r'$\theta_b = 8$'))
            C8_line = mlines.Line2D([0], [0],color='C4', label=(r'$\theta_b = 16$'))
            C9_line = mlines.Line2D([0], [0],color='None', label=('(in phot/px.)'))
    
            ax.legend(handles=[C0_line, C1_line, C2_line, C3_line, C8_line, C9_line], ncol = 2, prop={'size':13}, labelspacing = 1)  

        if panel == 0:
            #Legend
            C0_line = mlines.Line2D([0], [0],color='C0', label=r'$r = 1\sigma_{\mathrm{PSF}}$')
            C1_line = mlines.Line2D([0], [0],color='C1', label=r'$r = 2\sigma_{\mathrm{PSF}}$')
            C2_line = mlines.Line2D([0], [0],color='C2', label=r'$r = 3\sigma_{\mathrm{PSF}}$')
            C3_line = mlines.Line2D([0], [0],color='C3', label=r'$r = 4\sigma_{\mathrm{PSF}}$')
            C9_line = mlines.Line2D([0], [0],color='k', label=r'SMLM (Uniform' + '\n' + 'illumination intensity)')
            C10_line = mlines.Line2D([0], [0],color='k', linestyle='--', label=r'Localization on ISM' + '\n' + 'reconstruction data')
            C11_line = mlines.Line2D([0], [0],color='k', linestyle=':', label=r'Localization on' + '\n' + 'Fourier reweighted ISM' + '\n' + 'reconstruction data')
    
            ax.legend(handles=[C0_line, C1_line, C2_line, C3_line, C9_line, C10_line, C11_line], prop={'size':13}, loc='upper left', bbox_to_anchor=(1.04,1), labelspacing = 1)

        plt.savefig(('./Figures/S42/'+str(int(panel))) + '.svg')
        
    if plotnumber == 'S43':
        CRLBx, musum, signalsum, Bisum = load_4('xp_r_triangulation_nocenter_doughnut_adaptedbg')
        CRLBx_SMLM, musum_SMLM, signalsum_SMLM, Bisum_SMLM = load_4('SMLM')
        s_list, xrange = load_xp_r('xp_r_triangulation_nocenter_doughnut_adaptedbg')
        
        fig, ax = plt.subplots(1,1)
        fig.set_size_inches(5, 5, forward=True)
        
        if panel == 1:
            ax.semilogy((xrange-theta[0])/(10**-9),CRLBx.T[:,0:4]*10**9)
            ax.semilogy((xrange-theta[0])/(10**-9),CRLBx_SMLM*10**9*np.ones(np.size(xrange)), color='k')
            ax.semilogy((xrange-theta[0])/(10**-9),CRLBx_ISM*10**9*np.ones(np.size(xrange)), color='k', linestyle='--')
            ax.semilogy((xrange-theta[0])/(10**-9),CRLBx_FR*10**9*np.ones(np.size(xrange)), color='k', linestyle=':')
            ax.set_xlabel('Emitter-focus distance ' + r'$(x_f - \theta_x)$' + ' (nm)')
            ax.set_ylabel(r'$x$'+'-Cramér-Rao lower bound (nm)')
            ax.grid()
            ax.set_xlim(left=-130, right = 130)
            ax.set_ylim(bottom = 0.6, top = 100)
            ax.tick_params(axis='both', which='major')
            ax.tick_params(axis='both', which='minor')
            
        if panel == 2:
            ax.plot((xrange-theta[0])/(10**-9),CRLBx_SMLM/CRLBx.T[:,0:4])
            ax.plot((xrange-theta[0])/(10**-9),imp_ISM*np.ones(np.size(xrange)), color='k', linestyle='--')
            ax.plot((xrange-theta[0])/(10**-9),imp_FR*np.ones(np.size(xrange)), color='k', linestyle=':')
            ax.set_xlabel('Emitter-focus distance ' + r'$(x_f - \theta_x)$' + ' (nm)')
            ax.set_ylabel('Improvement over SMLM in ' + r'$x$'+ '-direction')
            ax.grid()
            ax.set_xlim(left=-130, right = 130)
            ax.set_ylim(bottom = 0, top = 3.5)
            ax.tick_params(axis='both', which='major')
            ax.tick_params(axis='both', which='minor')                     

        if panel == 3:
            ax.plot((xrange-theta[0])/(10**-9),musum.T[:,0:4])
            ax.plot((xrange-theta[0])/(10**-9),musum_SMLM*np.ones(np.size(xrange)), color='k')
            ax.set_xlabel('Emitter-focus distance ' + r'$(x_f - \theta_x)$' + ' (nm)')
            ax.set_ylabel('Expected amount of photons')
            ax.grid()
            ax.set_xlim(left=-130, right = 130)
            ax.tick_params(axis='both', which='major')
            ax.tick_params(axis='both', which='minor')                  

        if panel == 4:
            ax.plot((xrange-theta[0])/(10**-9),signalsum.T[:,0:4])
            ax.plot((xrange-theta[0])/(10**-9),signalsum_SMLM*np.ones(np.size(xrange)), color='k')
            ax.set_xlabel('Emitter-focus distance ' + r'$(x_f - \theta_x)$' + ' (nm)')
            ax.set_ylabel('Illumination-intensity normalized' + '\n' + 'amount of signal photons')
            ax.grid()
            ax.set_xlim(left=-130, right = 130)
            ax.set_ylim(bottom = 0, top = 2100)
            ax.tick_params(axis='both', which='major')
            ax.tick_params(axis='both', which='minor')   

        if panel == 5:
            ax.plot((xrange-theta[0])/(10**-9),Bisum.T[:,0:4]/100)
            ax.plot((xrange-theta[0])/(10**-9),Bisum_SMLM/100*np.ones(np.size(xrange)), color='k')
            ax.set_xlabel('Emitter-focus distance ' + r'$(x_f - \theta_x)$' + ' (nm)')
            ax.set_ylabel('Illumination-intensity normalized average' + '\n' + 'amount of background photons per pixel')
            ax.grid()
            ax.set_xlim(left=-130, right = 130)
            ax.set_ylim(bottom = 0, top = 9)
            ax.tick_params(axis='both', which='major')
            ax.tick_params(axis='both', which='minor')  

        if panel == 6:
            CRLBx, musum, signalsum, Bisum = load_4('thetaI_thetab_triangulation_nocenter_doughnut_adaptedbg')
            thetab_list, thetaIrange = load_thetaI_thetab('thetaI_thetab_triangulation_nocenter_doughnut_adaptedbg')
            
            ax.loglog(thetaIrange,CRLBx.T*10**9)
            ax.set_xlabel('Expected signal photon budget ' + r'$\theta_I$' + ' (photons)')
            ax.set_ylabel(r'$x$'+'-Cramér-Rao lower bound (nm)')
            ax.grid()
            ax.set_xlim(left=100, right = 10000)
            ax.set_ylim(bottom = 0.1, top = 40)
            ax.tick_params(axis='both', which='major')
            ax.tick_params(axis='both', which='minor')  
            
            #Legend
            C0_line = mlines.Line2D([0], [0],color='C0', label=(r'$\theta_b = 0$'))
            C1_line = mlines.Line2D([0], [0],color='C1', label=(r'$\theta_b = 1$'))
            C2_line = mlines.Line2D([0], [0],color='C2', label=(r'$\theta_b = 4$'))
            C3_line = mlines.Line2D([0], [0],color='C3', label=(r'$\theta_b = 8$'))
            C8_line = mlines.Line2D([0], [0],color='C4', label=(r'$\theta_b = 16$'))
            C9_line = mlines.Line2D([0], [0],color='None', label=('(in phot/px.)'))
    
            ax.legend(handles=[C0_line, C1_line, C2_line, C3_line, C8_line, C9_line], ncol = 2, prop={'size':13}, labelspacing = 1)    

        if panel == 0:
            #Legend
            C0_line = mlines.Line2D([0], [0],color='C0', label=r'$r = 1\sigma_{\mathrm{PSF}}$')
            C1_line = mlines.Line2D([0], [0],color='C1', label=r'$r = 2\sigma_{\mathrm{PSF}}$')
            C2_line = mlines.Line2D([0], [0],color='C2', label=r'$r = 3\sigma_{\mathrm{PSF}}$')
            C3_line = mlines.Line2D([0], [0],color='C3', label=r'$r = 4\sigma_{\mathrm{PSF}}$')
            C9_line = mlines.Line2D([0], [0],color='k', label=r'SMLM (Uniform' + '\n' + 'illumination intensity)')
            C10_line = mlines.Line2D([0], [0],color='k', linestyle='--', label=r'Localization on ISM' + '\n' + 'reconstruction data')
            C11_line = mlines.Line2D([0], [0],color='k', linestyle=':', label=r'Localization on' + '\n' + 'Fourier reweighted ISM' + '\n' + 'reconstruction data')
    
            ax.legend(handles=[C0_line, C1_line, C2_line, C3_line, C9_line, C10_line, C11_line], prop={'size':13}, loc='upper left', bbox_to_anchor=(1.04,1), labelspacing = 1)

        plt.savefig(('./Figures/S43/'+str(int(panel))) + '.svg')

    if plotnumber == 'S44':
        CRLBx, musum, signalsum, Bisum = load_4('xp_r_triangulation_rot_doughnut_ext_adaptedbg')
        CRLBx_SMLM, musum_SMLM, signalsum_SMLM, Bisum_SMLM = load_4('SMLM')
        s_list, xrange = load_xp_r('xp_r_triangulation_rot_doughnut_ext_adaptedbg')
        
        fig, ax = plt.subplots(1,1)
        fig.set_size_inches(5, 5, forward=True)
        
        if panel == 1:
            ax.semilogy((xrange-theta[0])/(10**-9),CRLBx.T*10**9)
            ax.semilogy((xrange-theta[0])/(10**-9),CRLBx_SMLM*10**9*np.ones(np.size(xrange)), color='k')
            ax.semilogy((xrange-theta[0])/(10**-9),CRLBx_ISM*10**9*np.ones(np.size(xrange)), color='k', linestyle='--')
            ax.semilogy((xrange-theta[0])/(10**-9),CRLBx_FR*10**9*np.ones(np.size(xrange)), color='k', linestyle=':')
            ax.set_xlabel('Emitter-focus distance ' + r'$(x_f - \theta_x)$' + ' (nm)')
            ax.set_ylabel(r'$x$'+'-Cramér-Rao lower bound (nm)')
            ax.grid()
            ax.set_xlim(left=-130, right = 130)
            ax.set_ylim(bottom = 0.1, top = 60)
            ax.tick_params(axis='both', which='major')
            ax.tick_params(axis='both', which='minor')
            
        if panel == 2:
            ax.plot((xrange-theta[0])/(10**-9),CRLBx_SMLM/CRLBx.T)
            ax.plot((xrange-theta[0])/(10**-9),imp_ISM*np.ones(np.size(xrange)), color='k', linestyle='--')
            ax.plot((xrange-theta[0])/(10**-9),imp_FR*np.ones(np.size(xrange)), color='k', linestyle=':')
            ax.set_xlabel('Emitter-focus distance ' + r'$(x_f - \theta_x)$' + ' (nm)')
            ax.set_ylabel('Improvement over SMLM in ' + r'$x$'+ '-direction')
            ax.grid()
            ax.set_xlim(left=-130, right = 130)
            ax.set_ylim(bottom = 0, top = 20)
            ax.tick_params(axis='both', which='major')
            ax.tick_params(axis='both', which='minor')                     

        if panel == 3:
            ax.plot((xrange-theta[0])/(10**-9),musum.T)
            ax.plot((xrange-theta[0])/(10**-9),musum_SMLM*np.ones(np.size(xrange)), color='k')
            ax.set_xlabel('Emitter-focus distance ' + r'$(x_f - \theta_x)$' + ' (nm)')
            ax.set_ylabel('Expected amount of photons')
            ax.grid()
            ax.set_xlim(left=-130, right = 130)
            ax.set_ylim(bottom = 0, top = 2900)
            ax.tick_params(axis='both', which='major')
            ax.tick_params(axis='both', which='minor')                  

        if panel == 4:
            ax.plot((xrange-theta[0])/(10**-9),signalsum.T)
            ax.plot((xrange-theta[0])/(10**-9),signalsum_SMLM*np.ones(np.size(xrange)), color='k')
            ax.set_xlabel('Emitter-focus distance ' + r'$(x_f - \theta_x)$' + ' (nm)')
            ax.set_ylabel('Illumination-intensity normalized' + '\n' + 'amount of signal photons')
            ax.grid()
            ax.set_xlim(left=-130, right = 130)
            ax.set_ylim(bottom = 0, top = 2100)
            ax.tick_params(axis='both', which='major')
            ax.tick_params(axis='both', which='minor')   

        if panel == 5:
            ax.plot((xrange-theta[0])/(10**-9),Bisum.T/100)
            ax.plot((xrange-theta[0])/(10**-9),Bisum_SMLM/100*np.ones(np.size(xrange)), color='k')
            ax.set_xlabel('Emitter-focus distance ' + r'$(x_f - \theta_x)$' + ' (nm)')
            ax.set_ylabel('Illumination-intensity normalized average' + '\n' + 'amount of background photons per pixel')
            ax.grid()
            ax.set_xlim(left=-130, right = 130)
            ax.set_ylim(bottom = 0, top = 9)
            ax.tick_params(axis='both', which='major')
            ax.tick_params(axis='both', which='minor')   

        if panel == 6:
            CRLBx, musum, signalsum, Bisum = load_4('thetaI_thetab_triangulation_rot_doughnut_adaptedbg')
            thetab_list, thetaIrange = load_thetaI_thetab('thetaI_thetab_triangulation_rot_doughnut_adaptedbg')
            
            ax.loglog(thetaIrange,CRLBx.T*10**9)
            ax.set_xlabel('Expected signal photon budget ' + r'$\theta_I$' + ' (photons)')
            ax.set_ylabel(r'$x$'+'-Cramér-Rao lower bound (nm)')
            ax.grid()
            ax.set_xlim(left=100, right = 10000)
            ax.set_ylim(bottom = 0.1, top = 40)
            ax.tick_params(axis='both', which='major')
            ax.tick_params(axis='both', which='minor')  
            
            #Legend
            C0_line = mlines.Line2D([0], [0],color='C0', label=(r'$\theta_b = 0$'))
            C1_line = mlines.Line2D([0], [0],color='C1', label=(r'$\theta_b = 1$'))
            C2_line = mlines.Line2D([0], [0],color='C2', label=(r'$\theta_b = 4$'))
            C3_line = mlines.Line2D([0], [0],color='C3', label=(r'$\theta_b = 8$'))
            C8_line = mlines.Line2D([0], [0],color='C4', label=(r'$\theta_b = 16$'))
            C9_line = mlines.Line2D([0], [0],color='None', label=('(in phot/px.)'))
    
            ax.legend(handles=[C0_line, C1_line, C2_line, C3_line, C8_line, C9_line], ncol = 2, prop={'size':13}, labelspacing = 1) 

        if panel == 0:
            #Legend
            C0_line = mlines.Line2D([0], [0],color='C0', label=r'$r = 0.5\sigma_{\mathrm{PSF}}$')
            C1_line = mlines.Line2D([0], [0],color='C1', label=r'$r = 1\sigma_{\mathrm{PSF}}$')
            C2_line = mlines.Line2D([0], [0],color='C2', label=r'$r = 2\sigma_{\mathrm{PSF}}$')
            C3_line = mlines.Line2D([0], [0],color='C3', label=r'$r = 3\sigma_{\mathrm{PSF}}$')
            C4_line = mlines.Line2D([0], [0],color='C4', label=r'$r = 4\sigma_{\mathrm{PSF}}$')
            C9_line = mlines.Line2D([0], [0],color='k', label=r'SMLM (Uniform' + '\n' + 'illumination intensity)')
            C10_line = mlines.Line2D([0], [0],color='k', linestyle='--', label=r'Localization on ISM' + '\n' + 'reconstruction data')
            C11_line = mlines.Line2D([0], [0],color='k', linestyle=':', label=r'Localization on' + '\n' + 'Fourier reweighted ISM' + '\n' + 'reconstruction data')
    
            ax.legend(handles=[C0_line, C1_line, C2_line, C3_line, C4_line, C9_line, C10_line, C11_line], prop={'size':13}, loc='upper left', bbox_to_anchor=(1.04,1), labelspacing = 1)

        plt.savefig(('./Figures/S44/'+str(int(panel))) + '.svg')

    if plotnumber == 'S45':
        CRLBx, musum, signalsum, Bisum = load_4('xp_r_triangulation_doughnut_ext_adaptedbg')
        CRLBx_SMLM, musum_SMLM, signalsum_SMLM, Bisum_SMLM = load_4('SMLM')
        s_list, xrange = load_xp_r('xp_r_triangulation_doughnut_ext_adaptedbg')
        
        fig, ax = plt.subplots(1,1)
        fig.set_size_inches(5, 5, forward=True)
        
        if panel == 1:
            ax.semilogy((xrange-theta[0])/(10**-9),CRLBx.T*10**9)
            ax.semilogy((xrange-theta[0])/(10**-9),CRLBx_SMLM*10**9*np.ones(np.size(xrange)), color='k')
            ax.semilogy((xrange-theta[0])/(10**-9),CRLBx_ISM*10**9*np.ones(np.size(xrange)), color='k', linestyle='--')
            ax.semilogy((xrange-theta[0])/(10**-9),CRLBx_FR*10**9*np.ones(np.size(xrange)), color='k', linestyle=':')
            ax.set_xlabel('Emitter-focus distance ' + r'$(x_f - \theta_x)$' + ' (nm)')
            ax.set_ylabel(r'$x$'+'-Cramér-Rao lower bound (nm)')
            ax.grid()
            ax.set_xlim(left=-130, right = 130)
            ax.set_ylim(bottom = 0.1, top = 60)
            ax.tick_params(axis='both', which='major')
            ax.tick_params(axis='both', which='minor')
            
        if panel == 2:
            ax.plot((xrange-theta[0])/(10**-9),CRLBx_SMLM/CRLBx.T)
            ax.plot((xrange-theta[0])/(10**-9),imp_ISM*np.ones(np.size(xrange)), color='k', linestyle='--')
            ax.plot((xrange-theta[0])/(10**-9),imp_FR*np.ones(np.size(xrange)), color='k', linestyle=':')
            ax.set_xlabel('Emitter-focus distance ' + r'$(x_f - \theta_x)$' + ' (nm)')
            ax.set_ylabel('Improvement over SMLM in ' + r'$x$'+ '-direction')
            ax.grid()
            ax.set_xlim(left=-130, right = 130)
            ax.set_ylim(bottom = 0, top = 20)
            ax.tick_params(axis='both', which='major')
            ax.tick_params(axis='both', which='minor')                     

        if panel == 3:
            ax.plot((xrange-theta[0])/(10**-9),musum.T)
            ax.plot((xrange-theta[0])/(10**-9),musum_SMLM*np.ones(np.size(xrange)), color='k')
            ax.set_xlabel('Emitter-focus distance ' + r'$(x_f - \theta_x)$' + ' (nm)')
            ax.set_ylabel('Expected amount of photons')
            ax.grid()
            ax.set_xlim(left=-130, right = 130)
            ax.set_ylim(bottom = 0, top = 2900)
            ax.tick_params(axis='both', which='major')
            ax.tick_params(axis='both', which='minor')                  

        if panel == 4:
            ax.plot((xrange-theta[0])/(10**-9),signalsum.T)
            ax.plot((xrange-theta[0])/(10**-9),signalsum_SMLM*np.ones(np.size(xrange)), color='k')
            ax.set_xlabel('Emitter-focus distance ' + r'$(x_f - \theta_x)$' + ' (nm)')
            ax.set_ylabel('Illumination-intensity normalized' + '\n' + 'amount of signal photons')
            ax.grid()
            ax.set_xlim(left=-130, right = 130)
            ax.set_ylim(bottom = 0, top = 2100)
            ax.tick_params(axis='both', which='major')
            ax.tick_params(axis='both', which='minor')   

        if panel == 5:
            ax.plot((xrange-theta[0])/(10**-9),Bisum.T/100)
            ax.plot((xrange-theta[0])/(10**-9),Bisum_SMLM/100*np.ones(np.size(xrange)), color='k')
            ax.set_xlabel('Emitter-focus distance ' + r'$(x_f - \theta_x)$' + ' (nm)')
            ax.set_ylabel('Illumination-intensity normalized average' + '\n' + 'amount of background photons per pixel')
            ax.grid()
            ax.set_xlim(left=-130, right = 130)
            ax.set_ylim(bottom = 0, top = 9)
            ax.tick_params(axis='both', which='major')
            ax.tick_params(axis='both', which='minor')   

        if panel == 6:
            CRLBx, musum, signalsum, Bisum = load_4('thetaI_thetab_triangulation_doughnut_adaptedbg')
            thetab_list, thetaIrange = load_thetaI_thetab('thetaI_thetab_triangulation_doughnut_adaptedbg')
            
            ax.loglog(thetaIrange,CRLBx.T*10**9)
            ax.set_xlabel('Expected signal photon budget ' + r'$\theta_I$' + ' (photons)')
            ax.set_ylabel(r'$x$'+'-Cramér-Rao lower bound (nm)')
            ax.grid()
            ax.set_xlim(left=100, right = 10000)
            ax.set_ylim(bottom = 0.1, top = 40)
            ax.tick_params(axis='both', which='major')
            ax.tick_params(axis='both', which='minor')  
            
            #Legend
            C0_line = mlines.Line2D([0], [0],color='C0', label=(r'$\theta_b = 0$'))
            C1_line = mlines.Line2D([0], [0],color='C1', label=(r'$\theta_b = 1$'))
            C2_line = mlines.Line2D([0], [0],color='C2', label=(r'$\theta_b = 4$'))
            C3_line = mlines.Line2D([0], [0],color='C3', label=(r'$\theta_b = 8$'))
            C8_line = mlines.Line2D([0], [0],color='C4', label=(r'$\theta_b = 16$'))
            C9_line = mlines.Line2D([0], [0],color='None', label=('(in phot/px.)'))
    
            ax.legend(handles=[C0_line, C1_line, C2_line, C3_line, C8_line, C9_line], ncol = 2, prop={'size':13}, labelspacing = 1)    
            
        if panel == 0:
            #Legend
            C0_line = mlines.Line2D([0], [0],color='C0', label=r'$r = 0.5\sigma_{\mathrm{PSF}}$')
            C1_line = mlines.Line2D([0], [0],color='C1', label=r'$r = 1\sigma_{\mathrm{PSF}}$')
            C2_line = mlines.Line2D([0], [0],color='C2', label=r'$r = 2\sigma_{\mathrm{PSF}}$')
            C3_line = mlines.Line2D([0], [0],color='C3', label=r'$r = 3\sigma_{\mathrm{PSF}}$')
            C4_line = mlines.Line2D([0], [0],color='C4', label=r'$r = 4\sigma_{\mathrm{PSF}}$')
            C9_line = mlines.Line2D([0], [0],color='k', label=r'SMLM (Uniform' + '\n' + 'illumination intensity)')
            C10_line = mlines.Line2D([0], [0],color='k', linestyle='--', label=r'Localization on ISM' + '\n' + 'reconstruction data')
            C11_line = mlines.Line2D([0], [0],color='k', linestyle=':', label=r'Localization on' + '\n' + 'Fourier reweighted ISM' + '\n' + 'reconstruction data')
    
            ax.legend(handles=[C0_line, C1_line, C2_line, C3_line, C4_line, C9_line, C10_line, C11_line], prop={'size':13}, loc='upper left', bbox_to_anchor=(1.04,1), labelspacing = 1)

        plt.savefig(('./Figures/S45/'+str(int(panel))) + '.svg')