# -*- coding: utf-8 -*-
"""
SMLM(theta, sigmapsf, roiparams): object containing the SMLM image formation 
model and methods to calculate the Cramér-Rao lower bound.

Inputs:
    theta: [thetax, thetay, thetaI, thetab]
        thetax: Emitter x-position (m)
        thetay: Emitter y-position (m)
        thetaI: Expected signal photon budget under maximum illumination (photons)
        thetab: Expected background photon budget under maximum illumination (photons/(m*m))
    sigmapsf: Standard deviation of Gaussian emission PSF (m)
    roiparams: [dx_i, Npixels]
        dx_i: Camera pixel size (m)
        Npixels: Amount of camera pixels per direction
                
Methods:
    - ErfDiff(x, dx, sigma) implements and returns Eq. (S16), i.e. E(x, dx, sigma**2).
    - DerivErfDiff(x, dx, sigma) implements and returns the derivative of Eq. (S16), i.e. E(x, dx, sigma**2).
    - ErfErf(coords) implements and returns the error function difference over the camera pixel with coordinates coords. Additionally, it returns the derivatives of this term with respect to thetax and thetay.
    - ExpectedValue(erferf) returns the Poisson mean mu_{i}. 
    - Deriv(sumerferf, derivsumerferfx, derivsumerferfy) returns the derivatives of the Poisson mean mu_{i,k} with respect to the entries of theta.
    - FisherMatrix() calculates and returns the Fisher information and underlying photon counts.
"""
import numpy as np
from scipy.special import erf

class SMLM():
    '''
    SMLM(theta, sigmapsf, roiparams): object containing the SMLM image formation 
    model and methods to calculate the Cramér-Rao lower bound.

    Inputs:
        theta: [thetax, thetay, thetaI, thetab]
            thetax: Emitter x-position (m)
            thetay: Emitter y-position (m)
            thetaI: Expected signal photon budget under maximum illumination (photons)
            thetab: Expected background photon budget under maximum illumination (photons/(m*m))
        sigmapsf: Standard deviation of Gaussian emission PSF (m)
        roiparams: [dx_i, Npixels]
            dx_i: Camera pixel size (m)
            Npixels: Amount of camera pixels per direction
                
    Methods:
        - ErfDiff(x, dx, sigma) implements and returns Eq. (S16), i.e. E(x, dx, sigma**2).
        - DerivErfDiff(x, dx, sigma) implements and returns the derivative of Eq. (S16), i.e. E(x, dx, sigma**2).
        - ErfErf(coords) implements and returns the error function difference over the camera pixel with coordinates coords. Additionally, it returns the derivatives of this term with respect to thetax and thetay.
        - ExpectedValue(erferf) returns the Poisson mean mu_{i}. 
        - Deriv(sumerferf, derivsumerferfx, derivsumerferfy) returns the derivatives of the Poisson mean mu_{i,k} with respect to the entries of theta.
        - FisherMatrix() calculates and returns the Fisher information and underlying photon counts.
    '''
    def __init__(self, theta, sigmapsf, roiparams):
        self.theta = theta      #[thetax, thetay, thetaI, thetab]
        self.thetax = theta[0]  #Emitter x-position (m)
        self.thetay = theta[1]  #Emitter y-position (m)
        self.thetaI = theta[2]  #Expected signal photon budget under maximum illumination (photons)
        self.thetab = theta[3]  #Expected background photon budget under maximum illumination (photons/(m*m))
        
        self.sigmapsf = sigmapsf        #Standard deviation of Gaussian emission PSF (m)

        self.dx_i = roiparams[0]          #Camera pixel size (m)
        self.Npixels = roiparams[1]       #Amount of camera pixels per direction
        self.x_i = np.linspace(self.dx_i/2, self.Npixels*self.dx_i-self.dx_i/2, self.Npixels) #Camera pixel x-coordinates (m)
        self.y_i = np.linspace(self.dx_i/2, self.Npixels*self.dx_i-self.dx_i/2, self.Npixels) #Camera pixel y-coordinates (m)
    
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
        return 1/(np.sqrt(2*np.pi)*sigma)*(np.exp(-((x+dx/2)**2)/(2*sigma**2)) - np.exp(-((x-dx/2)**2)/(2*sigma**2)))
    
    def ErfErf(self, coords):
        '''
        ErfErf(coords) implements and returns the error function difference
        over the camera pixel with coordinates coords.
        
        Additionally, it returns the derivatives of this term with respect to 
        thetax and thetay.
        
        This function is required to compute the Poisson mean mu as well as its derivatives.
        To avoid computational overhead, this function should be called explicitly instead of
        from within the individual calculations of mu and its derivatives.
        
        Inputs:
            coords: (1 x 2) array, containing the camera -coordinates of 
                the camera pixel
            
        Outputs:
            erferf: error function difference
            deriverferfx: derivative of erferf with respect to thetax
            deriverferfy: derivative of erferf with respect to thetay
        '''  
        erferf = self.ErfDiff(coords[0]-self.thetax, self.dx_i, self.sigmapsf)*self.ErfDiff(coords[1]-self.thetay, self.dx_i, self.sigmapsf)
        deriverferfx = -self.DerivErfDiff(coords[0]-self.thetax, self.dx_i, self.sigmapsf)*self.ErfDiff(coords[1]-self.thetay, self.dx_i, self.sigmapsf)
        deriverferfy = -self.ErfDiff(coords[0]-self.thetax, self.dx_i, self.sigmapsf)*self.DerivErfDiff(coords[1]-self.thetay, self.dx_i, self.sigmapsf)
        return erferf, deriverferfx, deriverferfy
    
    def ExpectedValue(self, erferf):
        '''
        ExpectedValue(erferf) returns the Poisson mean mu_{i}. 
        
        It requires explicit calculation of ErfErf(coords), which gives erferf.
        
        Input:
            erferf: error function difference
            
        Output:
            Poisson mean mu_{i}
        '''
        return self.thetaI*erferf+self.thetab*self.dx_i**2
    
    def Deriv(self, erferf, deriverferfx, deriverferfy):
        '''
        Deriv(sumerferf, derivsumerferfx, derivsumerferfy) returns the derivatives
        of the Poisson mean mu_{i,k} with respect to the entries of theta.
        
        It requires explicit calculation of ErfErf(coords), which gives
        erferf, deriverferfx and deriverferfy.

        Inputs:
            erferf: error function difference
            deriverferfx: derivative of erferf with respect to thetax
            deriverferfy: derivative of erferf with respect to thetay
            
        Output:
            deriv: derivatives of the Poisson mean mu_{i} with respect to the 
                entries of theta
        '''
        deriv = np.zeros((4,1))
        deriv[0,0] = self.thetaI*deriverferfx
        deriv[1,0] = self.thetaI*deriverferfy
        deriv[2,0] = erferf
        deriv[3,0] = self.dx_i**2
        return deriv
        
    def FisherMatrix(self):
        '''
        FisherMatrix() calculates and returns the Fisher information and underlying photon counts.

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
                coords = [self.x_i[i], self.y_i[j]]
                erferf, deriverferfx, deriverferfy = self.ErfErf(coords)
                mu = self.ExpectedValue(erferf)
                musum += mu
                Bisum += self.dx_i**2
                EfBgsum += self.thetab*self.dx_i**2
                signalsum += mu - self.thetab*self.dx_i**2
                deriv = self.Deriv(erferf, deriverferfx, deriverferfy)
                if mu < 1e-9:
                    mu = 1e-9
                I += 1/mu*deriv@deriv.T
        return I, musum, signalsum, Bisum, EfBgsum