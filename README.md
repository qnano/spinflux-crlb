General information
---
This software belongs to the article:
Quantifying the minimum localization uncertainty of image scanning localization microscopy

DOI: https://doi.org/10.1016/j.bpr.2024.100143

Data availability
---
The simulation data underlying the figures in the article can be found at:

DOI: https://doi.org/10.4121/21313230

Windows installation
---
The SpinFlux code was developed and tested in Python v3.7.10, running on Windows 10.

Steps:

1. Clone repository or download software from https://github.com/qnano/spinflux-crlb. If necessary, extract the code in the directory '/spinflux-crlb/'.

2. Install the Anaconda distribution: https://www.anaconda.com/products/individual. 

3. We recommend creating and activating a new environment for spinflux-crlb, as follows:

    - Open the Anaconda prompt.

    - Run the following commands in the terminal:
    ```
        conda create -n vti_env python=3.7.10 anaconda
        conda activate spinflux_env
    ```

SpinFlux class
---
The SpinFlux image formation model and Cramér-Rao lower bound (CRLB), as derived in the publication, are implemented as methods of the class ```SpinFlux(theta, sigmas, roiparams, pinholeparams, scenario, N_m, illumination_doughnut)``` in ```spinflux.py```, with the inputs representing:	
   - theta: [thetax, thetay, thetaI, thetab]
		* thetax: Emitter x-position (m)
		* thetay: Emitter y-position (m)
		* thetaI: Expected signal photon budget under maximum illumination (photons)
		* thetab: Expected background photon budget under maximum illumination (photons/(m*m))
   - sigmas: [sigmaillum, sigmapsf]
		* sigmaillum: Standard deviation of Gaussian illumination PSF (m)
		* sigmapsf: Standard deviation of Gaussian emission PSF (m)
   - roiparams: [dx_i, Npixels]
		* dx_i: Camera pixel size (m)
		* Npixels: Amount of camera pixels per direction
   - pinholeparams: (K x 3)-array [x_p, y_p, r_p]
		* K: Amount of pinholes and patterns
		* x_p: (K x 1) array of pinhole x-coordinates (m)
		* y_p: (K x 1) array of pinhole y-coordinates (m)
		* r_p: (K x 1) array of pinhole radii (m)
   - scenario: Simulation scenario under consideration, which models the normalization constant (i.e. to benchmark against SMLM) and pattern-dependent or pattern-independent background.
		* 'fixed_budget': Scenario for which the signal photon budget is exhausted
		* 'fixed_energy': Scenario for which the illumination power and time are constant per pattern
		* 'fixed_budget_adaptedbg': Scenario for which the signal photon budget is exhausted, with pattern-independent background
   - N_m: (Optional) Amount of mesh pixels per direction, default set to 100.	
   - illumination_doughnut: If True, doughnut-shaped intensity patterns with a zero-intensity minimum in the center will be used for illumination. If False, Gaussian illumination will be used.	

For example, the SpinFlux class for two x-separated patterns can be initiated as follows:
```
#Parameters
L_illum = 546*10**-9        #Illumination wavelength (m)
L_emission = 600*10**(-9)   #Emission wavelength (m)
NA = 1.35                   #Numerical aperture
    
dx_i = 65*10**-9    #Camera pixel size (m)
Npixels = 10        #Amount of camera pixels in each direction of the region of interest
    
x_p = Npixels*dx_i/2            #Pinhole x-coordinate (m)
y_p = Npixels*dx_i/2            #Pinhole y-coordinate (m)
r_p = 0.21*L_emission/NA        #Radius pinhole (m)
pinholeparams=np.array([[x_p-4*r_p/2, y_p, 3*r_p], [x_p+4*r_p/2, y_p, 3*r_p]])  
    
roiparams=[dx_i, Npixels]
sigmas = [0.21*L_illum/NA, 0.21*L_emission/NA]              #Illumination and emission PSF standard deviation (m)
theta = [Npixels*dx_i/2, Npixels*dx_i/2, 2000, 8/(dx_i**2)] #Emitter x- and y-position (m), expected signal photon budget (photons) and expected background (photons/m**2)
    
N_m = 100           #Amount of mesh pixels per direction

sf = SpinFlux(theta, sigmas, roiparams, pinholeparams, N_m = N_m, scenario=scenario, illumination_doughnut=False)
```

To draw sample Poisson-distributed measurements in the region of interest from the image formation model, the method ```Sample(amount)``` is used. The argument amount specifies the amount of sample regions of interest to be generated.
For example, 100 sample regions of interest can be drawn for the above example code as follows:
```
sample = sf.Sample(100)
```

To calculate the Fisher information matrix, the method ```FisherMatrix()``` is used. ```FisherMatrix()``` calculates and returns the Fisher information and underlying photon counts.
For example, the CRLB can be calculated for the above example code as follows:
```
I, musum, signalsum, Bisum, EfBgsum = sf.FisherMatrix()
I_inv = np.linalg.inv(I)
CRLBx = np.sqrt(I_inv[0,0])
```

To generate CRLB data for single-molecule localization microscopy (SMLM) to benchmark SpinFlux against, a class ```SMLM(theta, sigmapsf, roiparams)``` is provided in ```smlm.py``` with similar syntax as the SpinFlux class. For further details, please refer to the final section of this README. 

Data processing
---
We describe the simulation procedure that is needed to reproduce the figures from the article. For convenience, we also provide the existing simulation data on https://doi.org/10.4121/21313230. We therefore also describe how to reproduce the figures from the article from the existing simulation data.

### Simulation procedure

This section describes the simulation procedure to reproduce the simulation data underlying the publication. If you intend to (re)produce results from existing simulation data, please refer to the next subsection.
For the simulations underlying the publication, all necessary simulation code and parameters are already supplied in the file ```spinflux_simulation_utils.py```. To generate results that are not part of the simulation conditions in the publication, please refer to the previous section.

1. Verify that code is contained in the directory "/spinflux-crlb/".

2. Create an empty folder "/spinflux-crlb/SimData/". This folder will be used to store the data.

3. (If changing the simulation parameters is desired, open spinflux_simulation_utils.py. Navigate to the function simulate(plotnumber), in which parameters can be adjusted.)

4. Open ```spinflux_simulation_main.py```.

5. To run simulations:
   - To run all simulations underlying the publication, run spinflux_simulation_main.py. Please note that this requires significant calculation time.
   - To run specific simulations underlying the publication, call the function simulate(plotnumber). The argument plotnumber specifies the simulation type:
		* 'SMLM': SMLM Benchmark
        * '3-S4', '3-heatmap': Figure 3 & S4: 1 x-offset pattern, fixed signal photon budget
        * '4de-S6', '4bc-heatmap': Figure 4b-e & S6: 2 x-separated patterns, fixed signal photon budget
        * '4ij-S10', '4gh-heatmap': Figure 4g-j & S10: 3 patterns in equilateral triangle configuration without center pinhole, fixed signal photon budget
        * '5-S16', '5-heatmap': Figure 5 & S16: 4 doughnut-patterns in equilateral triangle configuration with center pinhole, fixed signal photon budget
        * 'S2': Figure S2: Effect of pinhole discretization on CRLB.
        * 'S5': Figure S5: 1 y-offset pattern, fixed signal photon budget
        * 'S7': Figure S7: 2 x-separated patterns with y-offset, fixed signal photon budget
        * 'S8': Figure S8: 2 x-separated patterns without pinholes, fixed signal photon budget
        * 'S9': Figure S9: 2 y-separated patterns, fixed signal photon budget
        * 'S11': Figure S11: 3 patterns in 90 degree rotated equilateral triangle configuration without center pinhole, fixed signal photon budget
        * 'S12': Figure S12: 4 patterns in equilateral triangle configuration with center pinhole, fixed signal photon budget
        * 'S13': Figure S13: 4 patterns in 90 degree rotated equilateral triangle configuration with center pinhole, fixed signal photon budget
        * 'S14': Figure S14: 3 doughnut-patterns in equilateral triangle configuration without center pinhole, fixed signal photon budget
        * 'S15': Figure S15: 3 doughnut-patterns in 90 degree rotated equilateral triangle configuration without center pinhole, fixed signal photon budget
        * 'S17': Figure S17: 4 doughnut-patterns in 90 degree rotated equilateral triangle configuration with center pinhole, fixed signal photon budget                
        * 'S18': Figure S18: 1 x-offset pattern, fixed illumination power and time
        * 'S19': Figure S19: 1 y-offset pattern, fixed illumination power and time
        * 'S20': Figure S20: 2 x-separated patterns, fixed illumination power and time
        * 'S21': Figure S21: 2 x-separated patterns with y-offset, fixed illumination power and time
        * 'S22': Figure S22: 2 x-separated patterns without pinholes, fixed illumination power and time
        * 'S23': Figure S23: 2 y-separated patterns, fixed illumination power and time
        * 'S24': Figure S24: 3 patterns in equilateral triangle configuration without center pinhole, fixed illumination power and time
        * 'S25': Figure S25: 3 patterns in 90 degree rotated equilateral triangle configuration without center pinhole, fixed illumination power and time
        * 'S26': Figure S26: 4 patterns in equilateral triangle configuration with center pinhole, fixed illumination power and time
        * 'S27': Figure S27: 4 patterns in 90 degree rotated equilateral triangle configuration with center pinhole, fixed illumination power and time
        * 'S28': Figure S28: 3 doughnut-patterns in equilateral triangle configuration without center pinhole, fixed illumination power and time
        * 'S29': Figure S29: 3 doughnut-patterns in 90 degree rotated equilateral triangle configuration without center pinhole, fixed illumination power and time
        * 'S30': Figure S30: 4 doughnut-patterns in equilateral triangle configuration with center pinhole, fixed fixed illumination power and time
        * 'S31': Figure S31: 4 doughnut-patterns in 90 degree rotated equilateral triangle configuration with center pinhole, fixed illumination power and time
        * 'S32': Figure S32: 1 x-offset pattern, fixed signal photon budget with pattern-independent background
        * 'S33': Figure S33: 1 y-offset pattern, fixed signal photon budget with pattern-independent background
        * 'S34': Figure S34: 2 x-separated patterns, fixed signal photon budget with pattern-independent background
        * 'S35': Figure S35: 2 x-separated patterns with y-offset, fixed signal photon budget with pattern-independent background
        * 'S36': Figure S36: 2 x-separated patterns without pinholes, fixed signal photon budget with pattern-independent background
        * 'S37': Figure S37: 2 y-separated patterns, fixed signal photon budget with pattern-independent background
        * 'S38': Figure S38: 3 patterns in equilateral triangle configuration without center pinhole, fixed signal photon budget with pattern-independent background
        * 'S39': Figure S39: 3 patterns in 90 degree rotated equilateral triangle configuration without center pinhole, fixed signal photon budget with pattern-independent background
        * 'S40': Figure S40: 4 patterns in equilateral triangle configuration with center pinhole, fixed signal photon budget with pattern-independent background
        * 'S41': Figure S41: 4 patterns in 90 degree rotated equilateral triangle configuration with center pinhole, fixed signal photon budget with pattern-independent background
        * 'S42': Figure S42: 3 doughnut-patterns in equilateral triangle configuration without center pinhole, fixed signal photon budget with pattern-independent background
        * 'S43': Figure S43: 3 doughnut-patterns in 90 degree rotated equilateral triangle configuration without center pinhole, fixed signal photon budget with pattern-independent background
        * 'S44': Figure S44: 4 doughnut-patterns in equilateral triangle configuration with center pinhole, fixed signal photon budget with pattern-independent background
        * 'S45': Figure S45: 4 doughnut-patterns in 90 degree rotated equilateral triangle configuration with center pinhole, fixed signal photon budget with pattern-independent background  
6. The simulation data will be saved in the directory "/spinflux-crlb/SimData/".

7. To do additional simulations, repeat step 7 to run the simulations.

### Using existing data

This section describes how to use the existing simulation data from https://doi.org/10.4121/21313230. If you intend to (re)produce results by running the simulation yourself, please refer to the previous subsection.

1. Create a directory "/spinflux-crlb/SimData/".

2. Download SimData.zip from https://doi.org/10.4121/21313230. SimData.zip contains all simulation data underlying the publication.

3. Extract SimData.zip to the directory "/spinflux-crlb/SimData/".

### Reproducing results

1. Verify that the directory "/spinflux-crlb/SimData/" exists and that it contains the necessary simulation data. Additionally, verify that code is contained in the directory "/spinflux-crlb/".

2. Create a directory "/spinflux-crlb/Figures/". Within this directory, create subdirectories corresponding to the plotnumbers (see step 4) of the figures to be created. E.g. to create Supplementary Figure 5, make sure the directory "/spinflux-crlb/Figures/S5/" exists.

3. Open ```spinflux_plot_main.py```.

4. To produce the desired figures from the article:
   - To produce all figures underlying the publication, run spinflux_plot_main.py.
   - To produce specific figures underlying the publication, call the function plot(plotnumber, panel). 
   - The argument plotnumber specifies the simulation type:
		* 'SMLM': SMLM Benchmark
        * '3-S4', '3-heatmap': Figure 3 & S4: 1 x-offset pattern, fixed signal photon budget
        * '4de-S6', '4bc-heatmap': Figure 4b-e & S6: 2 x-separated patterns, fixed signal photon budget
        * '4ij-S10', '4gh-heatmap': Figure 4g-j & S10: 3 patterns in equilateral triangle configuration without center pinhole, fixed signal photon budget
        * '5-S16', '5-heatmap': Figure 5 & S16: 4 doughnut-patterns in equilateral triangle configuration with center pinhole, fixed signal photon budget
        * 'S2': Figure S2: Effect of pinhole discretization on CRLB.
        * 'S5': Figure S5: 1 y-offset pattern, fixed signal photon budget
        * 'S7': Figure S7: 2 x-separated patterns with y-offset, fixed signal photon budget
        * 'S8': Figure S8: 2 x-separated patterns without pinholes, fixed signal photon budget
        * 'S9': Figure S9: 2 y-separated patterns, fixed signal photon budget
        * 'S11': Figure S11: 3 patterns in 90 degree rotated equilateral triangle configuration without center pinhole, fixed signal photon budget
        * 'S12': Figure S12: 4 patterns in equilateral triangle configuration with center pinhole, fixed signal photon budget
        * 'S13': Figure S13: 4 patterns in 90 degree rotated equilateral triangle configuration with center pinhole, fixed signal photon budget
        * 'S14': Figure S14: 3 doughnut-patterns in equilateral triangle configuration without center pinhole, fixed signal photon budget
        * 'S15': Figure S15: 3 doughnut-patterns in 90 degree rotated equilateral triangle configuration without center pinhole, fixed signal photon budget
        * 'S17': Figure S17: 4 doughnut-patterns in 90 degree rotated equilateral triangle configuration with center pinhole, fixed signal photon budget                
        * 'S18': Figure S18: 1 x-offset pattern, fixed illumination power and time
        * 'S19': Figure S19: 1 y-offset pattern, fixed illumination power and time
        * 'S20': Figure S20: 2 x-separated patterns, fixed illumination power and time
        * 'S21': Figure S21: 2 x-separated patterns with y-offset, fixed illumination power and time
        * 'S22': Figure S22: 2 x-separated patterns without pinholes, fixed illumination power and time
        * 'S23': Figure S23: 2 y-separated patterns, fixed illumination power and time
        * 'S24': Figure S24: 3 patterns in equilateral triangle configuration without center pinhole, fixed illumination power and time
        * 'S25': Figure S25: 3 patterns in 90 degree rotated equilateral triangle configuration without center pinhole, fixed illumination power and time
        * 'S26': Figure S26: 4 patterns in equilateral triangle configuration with center pinhole, fixed illumination power and time
        * 'S27': Figure S27: 4 patterns in 90 degree rotated equilateral triangle configuration with center pinhole, fixed illumination power and time
        * 'S28': Figure S28: 3 doughnut-patterns in equilateral triangle configuration without center pinhole, fixed illumination power and time
        * 'S29': Figure S29: 3 doughnut-patterns in 90 degree rotated equilateral triangle configuration without center pinhole, fixed illumination power and time
        * 'S30': Figure S30: 4 doughnut-patterns in equilateral triangle configuration with center pinhole, fixed fixed illumination power and time
        * 'S31': Figure S31: 4 doughnut-patterns in 90 degree rotated equilateral triangle configuration with center pinhole, fixed illumination power and time
        * 'S32': Figure S32: 1 x-offset pattern, fixed signal photon budget with pattern-independent background
        * 'S33': Figure S33: 1 y-offset pattern, fixed signal photon budget with pattern-independent background
        * 'S34': Figure S34: 2 x-separated patterns, fixed signal photon budget with pattern-independent background
        * 'S35': Figure S35: 2 x-separated patterns with y-offset, fixed signal photon budget with pattern-independent background
        * 'S36': Figure S36: 2 x-separated patterns without pinholes, fixed signal photon budget with pattern-independent background
        * 'S37': Figure S37: 2 y-separated patterns, fixed signal photon budget with pattern-independent background
        * 'S38': Figure S38: 3 patterns in equilateral triangle configuration without center pinhole, fixed signal photon budget with pattern-independent background
        * 'S39': Figure S39: 3 patterns in 90 degree rotated equilateral triangle configuration without center pinhole, fixed signal photon budget with pattern-independent background
        * 'S40': Figure S40: 4 patterns in equilateral triangle configuration with center pinhole, fixed signal photon budget with pattern-independent background
        * 'S41': Figure S41: 4 patterns in 90 degree rotated equilateral triangle configuration with center pinhole, fixed signal photon budget with pattern-independent background
        * 'S42': Figure S42: 3 doughnut-patterns in equilateral triangle configuration without center pinhole, fixed signal photon budget with pattern-independent background
        * 'S43': Figure S43: 3 doughnut-patterns in 90 degree rotated equilateral triangle configuration without center pinhole, fixed signal photon budget with pattern-independent background
        * 'S44': Figure S44: 4 doughnut-patterns in equilateral triangle configuration with center pinhole, fixed signal photon budget with pattern-independent background
        * 'S45': Figure S45: 4 doughnut-patterns in 90 degree rotated equilateral triangle configuration with center pinhole, fixed signal photon budget with pattern-independent background  
   - The argument panel specifies the panel to create:
		* 0: Legend
		* 1: CRLB as a function of spatial difference between pattern focus and emitter position
		* 2: Improvement of CRLB over SMLM as a function of spatial difference between pattern focus and emitter position
		* 3: Illumination-intensity normalized total amount of photons as a function of spatial difference between pattern focus and emitter position
		* 4: Illumination-intensity normalized total amount of signal photons as a function of spatial difference between pattern focus and emitter position
		* 5: Illumination-intensity normalized average amount of background photons per pixel as a function of spatial difference between pattern focus and emitter position
		* 6: CRLB as a function of signal photon budget for different background counts

5. To produce additional figures, repeat steps 2-4.

Description of individual Python files
---
### smlm.py
SMLM(theta, sigmapsf, roiparams): object containing the SMLM image formation 
model and methods to calculate the Cramér-Rao lower bound.

Inputs:
- theta: [thetax, thetay, thetaI, thetab]
    * thetax: Emitter x-position (m)
    * thetay: Emitter y-position (m)
    * thetaI: Expected signal photon budget under maximum illumination (photons)
    * thetab: Expected background photon budget under maximum illumination (photons/(m*m))
- sigmapsf: Standard deviation of Gaussian emission PSF (m)
- roiparams: [dx_i, Npixels]
    * dx_i: Camera pixel size (m)
    * Npixels: Amount of camera pixels per direction
                
Methods:
- ErfDiff(x, dx, sigma) implements and returns Eq. (S16), i.e. E(x, dx, sigma**2).
- DerivErfDiff(x, dx, sigma) implements and returns the derivative of Eq. (S16), i.e. E(x, dx, sigma**2).
- ErfErf(coords) implements and returns the error function difference over the camera pixel with coordinates coords. Additionally, it returns the derivatives of this term with respect to thetax and thetay.
- ExpectedValue(erferf) returns the Poisson mean mu_{i}. 
- Deriv(sumerferf, derivsumerferfx, derivsumerferfy) returns the derivatives of the Poisson mean mu_{i,k} with respect to the entries of theta.
- FisherMatrix() calculates and returns the Fisher information and underlying photon counts.

### spinflux.py
SpinFlux(theta, sigmas, roiparams, pinholeparams): object containing the SpinFlux 
image formation model and methods to calculate the Cramér-Rao lower bound.

Inputs:
- theta: [thetax, thetay, thetaI, thetab]
	* thetax: Emitter x-position (m)
	* thetay: Emitter y-position (m)
	* thetaI: Expected signal photon budget under maximum illumination (photons)
	* thetab: Expected background photon budget under maximum illumination (photons/(m*m))
- sigmas: [sigmaillum, sigmapsf]
	* sigmaillum: Standard deviation of Gaussian illumination PSF (m)
	* sigmapsf: Standard deviation of Gaussian emission PSF (m)
- roiparams: [dx_i, Npixels]
	* dx_i: Camera pixel size (m)
	* Npixels: Amount of camera pixels per direction
- pinholeparams: (K x 3)-array [x_p, y_p, r_p]
	* K: Amount of pinholes and patterns
	* x_p: (K x 1) array of pinhole x-coordinates (m)
	* y_p: (K x 1) array of pinhole y-coordinates (m)
	* r_p: (K x 1) array of pinhole radii (m)
- scenario: Simulation scenario under consideration, which models the normalization constant (i.e. to benchmark against SMLM) and pattern-dependent or pattern-independent background.
	* 'fixed_budget': Scenario for which the signal photon budget is exhausted
	* 'fixed_energy': Scenario for which the illumination power and time are constant per pattern
	* 'fixed_budget_adaptedbg': Scenario for which the signal photon budget is exhausted, with pattern-independent background
            
Optional input:
- N_m: Amount of mesh pixels per direction, default set to 100.
- illumination_doughnut: If True, doughnut-shaped intensity patterns with a zero-intensity minimum in the center will be used for illumination. If False, Gaussian illumination will be used.

Methods:
- in_Sp(idx, idy, approximation) returns True if the specified coordinate of mesh pixel (idx, idy) lies in the area S_p of the aperture with radius r_p, else it returns False.
- in_Ai(idx, idy, pidx, pidy) returns True if the center coordinate of mesh pixel (idx, idy) lies in the area Ai of camera pixel (pidx, pidy), else it returns False.
- pinhole_mesh() returns a 2-dimensional (N_m x N_m) array, of which the entries are 1 if the associated mesh-coordinates lie in the pinhole aperture, and 0 otherwise.
- pixel_mesh() returns a 2-dimensional (N_m x N_m) array, of which the entries are 1 if the associated mesh-coordinates lie in the pixel with index (pidx, pidy), and 0 otherwise.
- intersection_mesh() returns a 2-dimensional (N_m x N_m) array intersection, of which the entries are 1 if the associated mesh-coordinates lie in the overlap of the pinhole and the pixel with index (pidx, pidy), and 0 otherwise. 
- ErfDiff(x, dx, sigma) implements and returns Eq. (S16), i.e. E(x, dx, sigma**2).
- DoughnutAntiderivative(x,y,k) implements and returns Eq. (S35), i.e. F(x,y).
- DerivErfDiff(x, dx, sigma) implements and returns the derivative of Eq. (S16), i.e. E(x, dx, sigma**2).
- SumErfErf(intcoords) implements and returns the summed products of error function difference over the mesh-pixels in the overlap between pixel and pinhole, as shown in Eq. (S19).
- EffectiveBg(intcoords, k) returns the pattern-dependent background constant (i.e. B_{i,k} as in Eq. (S40)) over all mesh-coordinates contained in intcoords, for pattern k.
- EffectiveBg(intcoords, k) returns the pattern-independent background constant (i.e. B_{i} as in Eq. (S43)) over all mesh-coordinates contained in intcoords.
- ExpectedValue(sumerferf, Bi) returns the Poisson mean mu_{i,k} for all patterns k (as in Eq. (S39)). 
- Deriv(sumerferf, derivsumerferfx, derivsumerferfy, Bi) returns the derivatives of the Poisson mean mu_{i,k} with respect to the entries of theta, for all patterns k.
- Sample(amount) returns amount samples of Poisson-distributed measurements in the region of interest.
- FisherMatrix() calculates and returns the Fisher information and underlying photon counts.

### spinflux_plot_main.py
Plots panels of all the SpinFlux figures.
- plot(plotnumber, panel) plots the SpinFlux Cramér-Rao lower bound (CRLB) figures.

### spinflux_plot_utils.py
Utilities and plotting code underlying the SpinFlux figures.
- load_thetaI_Nm(fileextension) loads simulation data for the simulation file described in fileextension.
- load_thetaI_Nm(fileextension) loads simulation conditions for any configuration, as a function of the signal photon budget thetaI for different background counts thetab.
- load_xp_rp(fileextension) loads simulation conditions for the single-pattern configuration, as a function of the pattern focus x-position for different pinhole radii rp.
- load_yp_rp(fileextension) loads simulation conditions for the single-pattern configuration, as a function of the pattern focus y-position for different pinhole radii rp.
- load_xp_s(fileextension) loads simulation conditions for the two-pattern configuration, as a function of the pattern focus x-position for different separations s.
- load_xp_r(fileextension) loads simulation conditions for the triangle configuration, as a function of the pattern focus x-position for different spacings r.
- load_thetaI_Nm(fileextension) loads simulation conditions about the effect of discretization on the CRLB, as a function of the amount of mesh pixels Nm for different values of the signal photon count thetaI.
- plot(plotnumber, panel) plots the SpinFlux Cramér-Rao lower bound (CRLB) figures.

### spinflux_simulation_main.py
Simulates all SpinFlux Cramér-Rao lower bound data underlying figures.
- simulate(plotnumber) simulates the SpinFlux Cramér-Rao lower bound (CRLB) 
    data underlying the figure specified in plotnumber.

### spinflux_simulation_utils.py
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