# -*- coding: utf-8 -*-
"""
SpinFlux(theta, sigmas, roiparams, pinholeparams): object containing the SpinFlux 
image formation model and methods to calculate the Cramér-Rao lower bound.

Inputs:
    theta: [thetax, thetay, thetaI, thetab]
        thetax: Emitter x-position (m)
        thetay: Emitter y-position (m)
        thetaI: Expected signal photon budget under maximum illumination (photons)
        thetab: Expected background photon budget under maximum illumination (photons/(m*m))
    sigmas: [sigmaillum, sigmapsf]
        sigmaillum: Standard deviation of Gaussian illumination PSF (m)
        sigmapsf: Standard deviation of Gaussian emission PSF (m)
    roiparams: [dx_i, Npixels]
        dx_i: Camera pixel size (m)
        Npixels: Amount of camera pixels per direction
    pinholeparams: (K x 3)-array [x_p, y_p, r_p]
        K: Amount of pinholes and patterns
        x_p: (K x 1) array of pinhole x-coordinates (m)
        y_p: (K x 1) array of pinhole y-coordinates (m)
        r_p: (K x 1) array of pinhole radii (m)
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
"""
import numpy as np
from scipy.special import erf
    
class SpinFlux():
    '''
    SpinFlux(theta, sigmas, roiparams, pinholeparams, scenario): object containing the SpinFlux 
    image formation model and methods to calculate the Cramér-Rao lower bound.

    Inputs:
        theta: [thetax, thetay, thetaI, thetab]
            thetax: Emitter x-position (m)
            thetay: Emitter y-position (m)
            thetaI: Expected signal photon budget under maximum illumination (photons)
            thetab: Expected background photon budget under maximum illumination (photons/(m*m))
        sigmas: [sigmaillum, sigmapsf]
            sigmaillum: Standard deviation of Gaussian illumination PSF (m)
            sigmapsf: Standard deviation of Gaussian emission PSF (m)
        roiparams: [dx_i, Npixels]
            dx_i: Camera pixel size (m)
            Npixels: Amount of camera pixels per direction
        pinholeparams: (K x 3)-array [x_p, y_p, r_p]
            K: Amount of pinholes and patterns
            x_p: (K x 1) array of pinhole x-coordinates (m)
            y_p: (K x 1) array of pinhole y-coordinates (m)
            r_p: (K x 1) array of pinhole radii (m)
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
    '''
    def __init__(self, theta, sigmas, roiparams, pinholeparams, scenario, N_m = 100, illumination_doughnut=False):
        self.theta = theta      #[thetax, thetay, thetaI, thetab]
        self.thetax = theta[0]  #Emitter x-position (m)
        self.thetay = theta[1]  #Emitter y-position (m)
        self.thetaI = theta[2]  #Expected signal photon budget under maximum illumination (photons)
        self.thetab = theta[3]  #Expected background photon budget under maximum illumination (photons/(m*m))

        self.sigmaillum = sigmas[0]     #Standard deviation of Gaussian illumination PSF (m)
        self.sigmapsf = sigmas[1]       #Standard deviation of Gaussian emission PSF (m)

        self.dx_i = roiparams[0]          #Camera pixel size (m)
        self.Npixels = roiparams[1]       #Amount of camera pixels per direction
        self.x_i = np.linspace(self.dx_i/2, self.Npixels*self.dx_i-self.dx_i/2, self.Npixels) #Camera pixel x-coordinates (m)
        self.y_i = np.linspace(self.dx_i/2, self.Npixels*self.dx_i-self.dx_i/2, self.Npixels) #Camera pixel y-coordinates (m)

        self.x_p = pinholeparams[:,0]     #Pinhole x-coordinates (m)
        self.y_p = pinholeparams[:,1]     #Pinhole y-coordinates (m)
        self.r_p = pinholeparams[:,2]     #Radius pinhole (m)
        self.K = np.size(pinholeparams, axis=0) #Amount of patterns
        
        self.N_m = N_m                    #Amount of mesh pixels per direction
        self.dx_m = self.Npixels/self.N_m*self.dx_i                                           #Mesh pixel size (m)
        self.x_m = np.linspace(self.dx_m/2, self.N_m*self.dx_m-self.dx_m/2, self.N_m)         #Mesh pixel x-coordinates (m)
        self.y_m = np.linspace(self.dx_m/2, self.N_m*self.dx_m-self.dx_m/2, self.N_m)         #Mesh pixel x-coordinates (m)
        
        self.illumination_doughnut = illumination_doughnut #If True, doughnut-shaped intensity patterns with a zero-intensity minimum in the center will be used for illumination.
        
        if self.illumination_doughnut:
            self.Illum = np.exp(1)*((self.thetax-self.x_p)**2+(self.thetay-self.y_p)**2)/(2*self.sigmaillum**2)*np.exp((-(self.thetax-self.x_p)**2-(self.thetay-self.y_p)**2)/(2*self.sigmaillum**2))
            self.DerivIllumx = np.exp(1)*(self.thetax-self.x_p)/(self.sigmaillum**2)*np.exp((-(self.thetax-self.x_p)**2-(self.thetay-self.y_p)**2)/(2*self.sigmaillum**2)) + (self.x_p-self.thetax)/(self.sigmaillum**2)*self.Illum
            self.DerivIllumy = np.exp(1)*(self.thetay-self.y_p)/(self.sigmaillum**2)*np.exp((-(self.thetax-self.x_p)**2-(self.thetay-self.y_p)**2)/(2*self.sigmaillum**2)) + (self.y_p-self.thetay)/(self.sigmaillum**2)*self.Illum
        else:
            self.Illum = np.exp((-(self.thetax-self.x_p)**2-(self.thetay-self.y_p)**2)/(2*self.sigmaillum**2)) #Illumination intensity on emitter position
            self.DerivIllumx = (self.x_p-self.thetax)/(self.sigmaillum**2)*self.Illum #Derivative of illumination w.r.t. the emitter x-position
            self.DerivIllumy = (self.y_p-self.thetay)/(self.sigmaillum**2)*self.Illum #Derivative of illumination w.r.t. the emitter y-position
            
        self.scenario = scenario
        #Normalization constant for individual scenarios
        if self.scenario=='fixed_energy':
            self.a = 1/self.K #Normalization constant for the scenario where the illumination power and time are constant per pattern
        elif self.scenario == 'fixed_budget' or scenario == 'fixed_budget_adaptedbg':
            self.a = 1/np.sum(self.Illum) #Normalization constant for the scenario where the signal photon budget is exhausted
        
    def in_Sp(self, idx, idy, approximation):
        '''
        in_Sp(idx, idy, approximation) returns True if the specified coordinate of 
        mesh pixel (idx, idy) lies in the area S_p of the aperture with radius r_p,
        else it returns False.
        
        Inputs:
            idx: mesh pixel x-index
            idy: mesh pixel y-index
        
        The relevant coordinates are specified by the input approximation:
            approximation='left': bottom left coordinates of mesh pixel are used
            approximation='middle': center coordinates of mesh pixel are used
            approximation='right': upper right coordinates of mesh pixel are used
            
        Output: 
            True if the specified coordinate of mesh pixel (idx, idy) lies in the 
            area S_p of the aperture with radius r_p, else False.
        '''
        if approximation=='left':
            return (self.x_m[idx]-0.5*self.dx_m-self.x_p)**2 + (self.y_m[idy]-0.5*self.dx_m-self.y_p)**2 <= self.r_p**2
        elif approximation=='middle':
            return (self.x_m[idx]-self.x_p)**2 + (self.y_m[idy]-self.y_p)**2 <= self.r_p**2
        elif approximation=='right':
            return (self.x_m[idx]+0.5*self.dx_m-self.x_p)**2 + (self.y_m[idy]+0.5*self.dx_m-self.y_p)**2 <= self.r_p**2
        
    def in_Ai(self, idx, idy, pidx, pidy):
        '''
        in_Ai(idx, idy, pidx, pidy) returns True if the center coordinate of 
        mesh pixel (idx, idy) lies in the area Ai of camera pixel (pidx, pidy),
        else it returns False.
        
        Inputs:
            idx: mesh pixel x-index
            idy: mesh pixel y-index
            pidx: camera pixel x-index
            pidy: camera pixel y-index
            
        Output:
            True if the center coordinate of mesh pixel (idx, idy) lies in the 
            area Ai of camera pixel (pidx, pidy), else False.
        '''
        return (self.x_i[pidx] - self.dx_i/2 <= self.x_m[idx] and 
                self.x_i[pidx] + self.dx_i/2 >= self.x_m[idx] and
                self.y_i[pidy] - self.dx_i/2 <= self.y_m[idy] and
                self.y_i[pidy] + self.dx_i/2 >= self.y_m[idy] )
        
    def pinhole_mesh(self, approximation='middle'):
        '''
        pinhole_mesh() returns a 2-dimensional (N_m x N_m) array, of which the entries are 1 if
        the associated mesh-coordinates lie in the pinhole aperture, and 0 otherwise.
        
        The relevant coordinates are specified by the input approximation:
            approximation='left': bottom left coordinates of mesh pixel are used
            approximation='middle': center coordinates of mesh pixel are used
            approximation='right': upper right coordinates of mesh pixel are used
            
        Output:
            (N_m x N_m) array, of which the entries are 1 if the associated 
            mesh-coordinates lie in the pinhole aperture, and 0 otherwise.
        '''
        pinhole = np.zeros((self.N_m, self.N_m, self.K))
        for (idx, idy, k), _ in np.ndenumerate(pinhole):
            if (self.in_Sp(idx, idy, approximation))[k]:
                pinhole[idx, idy, k] = 1
        return pinhole
    
    def pixel_mesh(self, pidx, pidy):
        '''
        pixel_mesh() returns a 2-dimensional (N_m x N_m) array, of which the entries are 1 if
        the associated mesh-coordinates lie in the pixel with index (pidx, pidy),
        and 0 otherwise.
        
        Inputs:
            pidx: camera pixel x-index
            pidy: camera pixel y-index
            
        Output:
            (N_m x N_m) array, of which the entries are 1 if the associated 
            mesh-coordinates lie in the pixel with index (pidx, pidy),
            and 0 otherwise.
        '''
        pixel = np.zeros((self.N_m, self.N_m, self.K))
        for (idx, idy), _ in np.ndenumerate(pixel[:,:,0]):
            if self.in_Ai(idx, idy, pidx, pidy):
                pixel[idx, idy, :] = 1
        return pixel
    
    def intersection_mesh(self, pidx, pidy, approximation='middle'):
        '''
        intersection_mesh() returns a 2-dimensional (N_m x N_m) array intersection, of which the 
        entries are 1 if the associated mesh-coordinates lie in the overlap of the pinhole and the
        pixel with index (pidx, pidy), and 0 otherwise. 
        
        Additionally, it returns a 2-dimensional (N_overlap x 2) array intcoords, containing
        the mesh-coordinates of all N_overlap mesh-pixels in the overlapping area.

        Inputs:
            pidx: camera pixel x-index
            pidy: camera pixel y-index
        
        The relevant coordinates are specified by the input approximation:
            approximation='left': bottom left coordinates of mesh pixel are used
            approximation='middle': center coordinates of mesh pixel are used
            approximation='right': upper right coordinates of mesh pixel are used
            
        Outputs:
            intersection: (N_m x N_m) array, of which the entries are 1 if the 
                associated mesh-coordinates lie in the overlap of the pinhole and
                the pixel with index (pidx, pidy), and 0 otherwise. 
            intcoords:(N_overlap x 2) array, containing the mesh-coordinates of 
                all N_overlap mesh-pixels in the overlapping area.
                
        '''
        intersection = self.pixel_mesh(pidx, pidy)*(self.pinhole_mesh(approximation))
        intcoords=[[] for x in range(self.K)]
        for (idx, idy, k), _ in np.ndenumerate(intersection):
            if intersection[idx,idy,k]:
                intcoords[k].append([self.x_m[idx], self.y_m[idy]])
        return intersection, intcoords
    
    def ErfDiff(self, x, dx, sigma):
        '''
        ErfDiff(x, dx, sigma) implements and returns Eq. (S16), i.e. E(x, dx, sigma**2).
        It represents the difference of two error functions over the boundaries of the pixel,
        to evaluate the definite integral of a Gaussian function over the pixel area.
        
        Inputs:
            x: x-coordinate in E(x, dx, sigma**2)
            dx: pixel size in E(x, dx, sigma**2)
            sigma: variance deviation, i.e. sigma in E(x, dx, sigma**2)
            
        Outputs:
            E(x, dx, sigma**2)
        '''
        return 0.5*erf((x+dx/2)/(sigma*np.sqrt(2)))-0.5*erf((x-dx/2)/(sigma*np.sqrt(2)))

    def DoughnutAntiderivative(self, x, y, k):
        '''
        DoughnutAntiderivative implements and returns Eq. (S35), i.e. F(x,y).
        It represents the antiderivative w.r.t. x and y of the PSF convolved with the
        illumination pattern, needed to compute the pattern-dependent background
        under doughnut-shaped illumination.
        
        Inputs:
            x: x-coordinate
            y: y-coordinate
            k: pattern index
            
        Outputs:
            F(x,y)
        '''
        erfx = erf( (x-self.x_p[k])/(np.sqrt(2*self.sigmapsf**2 + 2*self.sigmaillum**2)) )
        erfy = erf( (y-self.y_p[k])/(np.sqrt(2*self.sigmapsf**2 + 2*self.sigmaillum**2)) )     
        expx = np.exp(-((x-self.x_p[k])**2)/(2*self.sigmapsf**2 + 2*self.sigmaillum**2))
        expy = np.exp(-((y-self.y_p[k])**2)/(2*self.sigmapsf**2 + 2*self.sigmaillum**2))
        factor1 = np.exp(1)*np.pi*self.sigmaillum**2/2
        factor2 = -(np.exp(1)*np.sqrt(np.pi)*self.sigmaillum**4)/(np.sqrt((2*self.sigmapsf**2 + 2*self.sigmaillum**2)**3))
        return factor1*erfx*erfy + factor2*( (x-self.x_p[k])*expx*erfy + (y-self.y_p[k])*expy*erfx)
    
    def DerivErfDiff(self, x, dx, sigma):
        '''
        DerivErfDiff(x, dx, sigma) implements and returns the derivative of Eq. (S16), i.e. E(x, dx, sigma**2).
        It represents the derivative with respect to x of the difference of two error 
        functions over the boundaries of the pixel, to evaluate the definite integral 
        of a Gaussian function over the pixel area.
        
        Inputs:
            x: x-coordinate in E(x, dx, sigma**2)
            dx: pixel size in E(x, dx, sigma**2)
            sigma: variance deviation, i.e. sigma in E(x, dx, sigma**2)
            
        Output:
            Derivative of E(x, dx, sigma**2) with respect to x
        '''        
        return 1/(np.sqrt(2*np.pi)*sigma)*(np.exp(-(x+dx/2)**2/(2*sigma**2)) - np.exp(-(x-dx/2)**2/(2*sigma**2)))
    
    def SumErfErf(self, intcoords):
        '''
        SumErfErf(intcoords) implements and returns the summed products of error function difference
        over the mesh-pixels in the overlap between pixel and pinhole, as shown in Eq. (S19).
        
        Additionally, it returns the derivatives of this sum-term with respect to thetax and thetay.
        
        This function is required to compute the Poisson mean mu as well as its derivatives.
        To avoid computational overhead, this function should be called explicitly instead of
        from within the individual calculations of mu and its derivatives.
        
        Inputs:
            intcoords: (N_overlap x 2) array, containing the mesh-coordinates of 
                all N_overlap mesh-pixels in the overlapping area.
            
        Outputs:
            sumerferf: sum of error function differences of Eq. (S19)
            derivsumerferfx: derivative of sumerferf with respect to thetax
            derivsumerferfy: derivative of sumerferf with respect to thetay
        '''     
        sumerferf = 0
        derivsumerferfx = 0
        derivsumerferfy = 0
        for coords in intcoords:
            sumerferf += self.ErfDiff(coords[0]-self.thetax, self.dx_m, self.sigmapsf)*self.ErfDiff(coords[1]-self.thetay, self.dx_m, self.sigmapsf)
            derivsumerferfx += -self.DerivErfDiff(coords[0]-self.thetax, self.dx_m, self.sigmapsf)*self.ErfDiff(coords[1]-self.thetay, self.dx_m, self.sigmapsf)
            derivsumerferfy += -self.ErfDiff(coords[0]-self.thetax, self.dx_m, self.sigmapsf)*self.DerivErfDiff(coords[1]-self.thetay, self.dx_m, self.sigmapsf)
        return sumerferf, derivsumerferfx, derivsumerferfy
        
    def EffectiveBg(self, intcoords, k):
        '''
        EffectiveBg(intcoords, k) returns the pattern-dependent background constant
        (i.e. B_{i,k} as in Eq. (S40)) over all mesh-coordinates contained in intcoords,
        for pattern k, or the equivalent for doughnut-shaped illumination (Eq. (S34)).
        
        This function is required to compute the Poisson mean mu as well as its derivatives.
        To avoid computational overhead, this function should be called explicitly instead of
        from within the individual calculations of mu and its derivatives.

        Inputs:
            intcoords: (N_overlap x 2) array, containing the mesh-coordinates of 
                all N_overlap mesh-pixels in the overlapping area.
            k: pattern index
            
        Output:
            Pattern-dependent background constant (i.e. B_{i,k} as in Eq. (S40))
        '''
        if self.illumination_doughnut:
            sumbg = 0
            for coords in intcoords:
                sumbg += (self.DoughnutAntiderivative(coords[0]+self.dx_m/2, coords[1]+self.dx_m/2, k) 
                          - self.DoughnutAntiderivative(coords[0]+self.dx_m/2, coords[1]-self.dx_m/2, k) 
                          - self.DoughnutAntiderivative(coords[0]-self.dx_m/2, coords[1]+self.dx_m/2, k) 
                          + self.DoughnutAntiderivative(coords[0]-self.dx_m/2, coords[1]-self.dx_m/2, k))
            return sumbg        
        else:
            sumerferfbg = 0
            sigmabg = np.sqrt(self.sigmapsf**2 + self.sigmaillum**2)
            for coords in intcoords:
                sumerferfbg += self.ErfDiff(coords[0]-self.x_p[k], self.dx_m, sigmabg)*self.ErfDiff(coords[1]-self.y_p[k], self.dx_m, sigmabg)
            return 2*np.pi*self.sigmaillum**2*sumerferfbg

    def EffectiveBg_adaptedbg(self, intcoords):
        '''
        EffectiveBg(intcoords, k) returns the pattern-independent background constant
        (i.e. B_{i} as in Eq. (S43)) over all mesh-coordinates contained in intcoords.
        
        This function is required to compute the Poisson mean mu as well as its derivatives.
        To avoid computational overhead, this function should be called explicitly instead of
        from within the individual calculations of mu and its derivatives.
        
        Inputs:
            intcoords: (N_overlap x 2) array, containing the mesh-coordinates of 
                all N_overlap mesh-pixels in the overlapping area.
            
        Output:
            Pattern-independent background constant (i.e. B_{i} as in Eq. (S43))
        '''
        Bi = 0
        for coords in intcoords:
            Bi += (self.dx_m)**2
        return Bi
    
    def ExpectedValue(self, sumerferf, Bi):
        '''
        ExpectedValue(sumerferf, Bi) returns the Poisson mean mu_{i,k} for all 
        patterns k (as in Eq. (S39)). 
        
        It requires explicit calculation of SumErfErf(intcoords) and 
        EffectiveBg(intcoords,k) for all patterns k, which respectively give 
        sumerferf and Bi.
        
        Inputs:
            sumerferf: sum of error function differences of Eq. (S19)
            Bi: Pattern-(in)dependent background constant B_{i,k} (or B_{i})
            
        Output:
            Poisson mean mu_{i,k} (as in Eq. (S39))
        '''
        if self.scenario=='fixed_budget_adaptedbg':
            return self.a*self.thetaI*self.Illum*sumerferf+1/self.K*self.thetab*Bi        
        else:
            return self.a*self.thetaI*self.Illum*sumerferf+self.a*self.thetab*Bi
    
    def Deriv(self, sumerferf, derivsumerferfx, derivsumerferfy, Bi):
        '''
        Deriv(sumerferf, derivsumerferfx, derivsumerferfy, Bi) returns the derivatives
        of the Poisson mean mu_{i,k} with respect to the entries of theta, for all patterns k.
        
        It requires explicit calculation of SumErfErf(intcoords) and 
        EffectiveBg(intcoords,k) for all patterns k, which respectively give 
        sumerferf, derivsumerferfx, derivsumerferfy and Bi.
        
        Inputs:
            sumerferf: sum of error function differences of Eq. (S19)
            derivsumerferfx: derivative of sumerferf with respect to thetax
            derivsumerferfy: derivative of sumerferf with respect to thetay
            Bi: Pattern-(in)dependent background constant B_{i,k} (or B_{i})
            
        Output:
            deriv: derivatives of the Poisson mean mu_{i,k} with respect to the 
                entries of theta
        '''
        deriv = np.zeros((4,self.K))
        deriv[0,:] = self.a*self.thetaI*self.DerivIllumx*sumerferf + self.a*self.thetaI*self.Illum*derivsumerferfx
        deriv[1,:] = self.a*self.thetaI*self.DerivIllumy*sumerferf + self.a*self.thetaI*self.Illum*derivsumerferfy
        deriv[2,:] = self.a*self.Illum*sumerferf
        if self.scenario=='fixed_budget_adaptedbg':
            deriv[3,:] = 1/self.K*Bi
        else:
            deriv[3,:] = self.a*Bi
        return deriv

    def Sample(self, amount, approximation='middle'):
        '''
        Sample(amount) returns amount samples of Poisson-distributed measurements
        in the region of interest.

        The relevant mesh coordinates are specified by the input approximation:
            approximation='left': bottom left coordinates of mesh pixel are used
            approximation='middle': center coordinates of mesh pixel are used
            approximation='right': upper right coordinates of mesh pixel are used
            
        Input:
            amount: amount of Poisson-distributed measurements in the region of interest.
            
        Output:
            sample: (amount x K x Npixels x Npixels) array, containing Poisson-distributed 
            measurements for each pattern k in the region of interest.
        '''
        sample = np.zeros((amount, self.K, self.Npixels, self.Npixels))
        rnd = np.random.default_rng()
        for i in range(self.Npixels):
            for j in range(self.Npixels):
                _, intcoords = self.intersection_mesh(i,j,approximation)
                sumerferf = np.zeros(self.K)
                derivsumerferfx = np.zeros(self.K)
                derivsumerferfy = np.zeros(self.K)
                Bi = np.zeros(self.K)
                for k in range(self.K):
                    sumerferf[k], derivsumerferfx[k], derivsumerferfy[k] = self.SumErfErf(intcoords[k])
                    if self.scenario=='fixed_budget_adaptedbg':
                        Bi[k] = self.EffectiveBg_adaptedbg(intcoords[k])
                    else:
                        Bi[k] = self.EffectiveBg(intcoords[k],k)
                mu = self.ExpectedValue(sumerferf, Bi)
                sample[:,:,i,j]=rnd.poisson(mu, size=(amount,self.K))
        return sample
        
    def FisherMatrix(self, approximation='middle'):
        '''
        FisherMatrix() calculates and returns the Fisher information and underlying photon counts.
        
        The relevant mesh coordinates are specified by the input approximation:
            approximation='left': bottom left coordinates of mesh pixel are used
            approximation='middle': center coordinates of mesh pixel are used
            approximation='right': upper right coordinates of mesh pixel are used        
        
        Outputs:
            I: the Fisher information matrix
            musum: the total photon count
            signalsum: the total signal photon count
            Bisum: the sum of effective background coefficients B_i
            EfBgsum: the total background photon count in the region of interest
        '''
        I = 0
        musum = 0
        signalsum = 0
        Bisum = 0
        EfBgsum = 0
        for i in range(self.Npixels):
            for j in range(self.Npixels):
                _, intcoords = self.intersection_mesh(i,j,approximation)
                sumerferf = np.zeros(self.K)
                derivsumerferfx = np.zeros(self.K)
                derivsumerferfy = np.zeros(self.K)
                Bi = np.zeros(self.K)
                for k in range(self.K):
                    sumerferf[k], derivsumerferfx[k], derivsumerferfy[k] = self.SumErfErf(intcoords[k])
                    if self.scenario=='fixed_budget_adaptedbg':
                        Bi[k] = self.EffectiveBg_adaptedbg(intcoords[k])
                    else:
                        Bi[k] = self.EffectiveBg(intcoords[k],k)
                mu = self.ExpectedValue(sumerferf, Bi)
                musum += np.sum(mu)
                Bisum += np.sum(Bi)
                if self.scenario=='fixed_budget_adaptedbg':
                    EfBgsum += 1/self.K*self.thetab*np.sum(Bi)
                    signalsum += np.sum(mu) - 1/self.K*self.thetab*np.sum(Bi)
                else:
                    EfBgsum += self.a*self.thetab*np.sum(Bi)
                    signalsum += np.sum(mu) - self.a*self.thetab*np.sum(Bi)
                deriv = self.Deriv(sumerferf, derivsumerferfx, derivsumerferfy, Bi)
                mu[mu<1e-9] = 1e-9
                for k in range(self.K):
                    derivk = deriv[:,k].reshape((4,1))
                    I += 1/mu[k]*derivk@derivk.T
        return I, musum, signalsum, Bisum, EfBgsum