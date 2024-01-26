# -*- coding: utf-8 -*-
"""
Plots panels of all the SpinFlux figures.
    - plot(plotnumber, panel) plots the SpinFlux Cramér-Rao lower bound (CRLB) figures.
    
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
"""

from spinflux_plot_utils import plot

#Figure 3 & S4
plot('3-S4', 1)
plot('3-S4', 2)
plot('3-S4', 3)
plot('3-S4', 4)
plot('3-S4', 5)
plot('3-S4', 6)
plot('3-S4', 0)

#Figure 3 (heatmap)
plot('3-heatmap', 1)
plot('3-heatmap', 2)

#Figure 4de & S6
plot('4de-S6', 1)
plot('4de-S6', 2)
plot('4de-S6', 3)
plot('4de-S6', 4)
plot('4de-S6', 5)
plot('4de-S6', 6)
plot('4de-S6', 0)

#Figure 4bc (heatmap)
plot('4bc-heatmap', 1)
plot('4bc-heatmap', 2)

#Figure 4ij & S10
plot('4ij-S10', 1)
plot('4ij-S10', 2)
plot('4ij-S10', 3)
plot('4ij-S10', 4)
plot('4ij-S10', 5)
plot('4ij-S10', 6)
plot('4ij-S10', 0)

#Figure 4gh (heatmap)
plot('4gh-heatmap', 1)
plot('4gh-heatmap', 2)

#Figure 5 & S16
plot('5-S16', 1)
plot('5-S16', 2)
plot('5-S16', 3)
plot('5-S16', 4)
plot('5-S16', 5)
plot('5-S16', 6)
plot('5-S16', 0)

#Figure 5 (heatmap)
plot('5-heatmap', 1)
plot('5-heatmap', 2)

#Figure S2
plot('S2', 1)
plot('S2', 2)
plot('S2', 3)
plot('S2', 4)
plot('S2', 5)
plot('S2', 6)
plot('S2', 7)
plot('S2', 8)
plot('S2', 9)
plot('S2', 0)

#Figure S5
plot('S5', 1)
plot('S5', 2)
plot('S5', 3)
plot('S5', 4)
plot('S5', 5)
plot('S5', 0)

#Figure S7
plot('S7', 1)
plot('S7', 2)
plot('S7', 3)
plot('S7', 4)
plot('S7', 5)
plot('S7', 0)

#Figure S8
plot('S8', 1)
plot('S8', 2)
plot('S8', 3)
plot('S8', 4)
plot('S8', 5)
plot('S8', 6)
plot('S8', 0)

#Figure S9
plot('S9', 1)
plot('S9', 2)
plot('S9', 3)
plot('S9', 4)
plot('S9', 5)
plot('S9', 6)
plot('S9', 0)

#Figure S11
plot('S11', 1)
plot('S11', 2)
plot('S11', 3)
plot('S11', 4)
plot('S11', 5)
plot('S11', 6)
plot('S11', 0)

#Figure S12
plot('S12', 1)
plot('S12', 2)
plot('S12', 3)
plot('S12', 4)
plot('S12', 5)
plot('S12', 6)
plot('S12', 0)

#Figure S13
plot('S13', 1)
plot('S13', 2)
plot('S13', 3)
plot('S13', 4)
plot('S13', 5)
plot('S13', 6)
plot('S13', 0)

#Figure S14
plot('S14', 1)
plot('S14', 2)
plot('S14', 3)
plot('S14', 4)
plot('S14', 5)
plot('S14', 6)
plot('S14', 0)

#Figure S15
plot('S15', 1)
plot('S15', 2)
plot('S15', 3)
plot('S15', 4)
plot('S15', 5)
plot('S15', 6)
plot('S15', 0)

#Figure S17
plot('S17', 1)
plot('S17', 2)
plot('S17', 3)
plot('S17', 4)
plot('S17', 5)
plot('S17', 6)
plot('S17', 0)

#Figure S18
plot('S18', 1)
plot('S18', 2)
plot('S18', 3)
plot('S18', 4)
plot('S18', 5)
plot('S18', 6)
plot('S18', 0)

#Figure S19
plot('S19', 1)
plot('S19', 2)
plot('S19', 3)
plot('S19', 4)
plot('S19', 5)
plot('S19', 0)

#Figure S20
plot('S20', 1)
plot('S20', 2)
plot('S20', 3)
plot('S20', 4)
plot('S20', 5)
plot('S20', 6)
plot('S20', 0)

#Figure S21
plot('S21', 1)
plot('S21', 2)
plot('S21', 3)
plot('S21', 4)
plot('S21', 5)
plot('S21', 0)

#Figure S22
plot('S22', 1)
plot('S22', 2)
plot('S22', 3)
plot('S22', 4)
plot('S22', 5)
plot('S22', 6)
plot('S22', 0)

#Figure S23
plot('S23', 1)
plot('S23', 2)
plot('S23', 3)
plot('S23', 4)
plot('S23', 5)
plot('S23', 6)
plot('S23', 0)

#Figure S24
plot('S24', 1)
plot('S24', 2)
plot('S24', 3)
plot('S24', 4)
plot('S24', 5)
plot('S24', 6)
plot('S24', 0)

#Figure S25
plot('S25', 1)
plot('S25', 2)
plot('S25', 3)
plot('S25', 4)
plot('S25', 5)
plot('S25', 6)
plot('S25', 0)

#Figure S26
plot('S26', 1)
plot('S26', 2)
plot('S26', 3)
plot('S26', 4)
plot('S26', 5)
plot('S26', 6)
plot('S26', 0)

#Figure S27
plot('S27', 1)
plot('S27', 2)
plot('S27', 3)
plot('S27', 4)
plot('S27', 5)
plot('S27', 6)
plot('S27', 0)

#Figure S28
plot('S28', 1)
plot('S28', 2)
plot('S28', 3)
plot('S28', 4)
plot('S28', 5)
plot('S28', 6)
plot('S28', 0)

#Figure S29
plot('S29', 1)
plot('S29', 2)
plot('S29', 3)
plot('S29', 4)
plot('S29', 5)
plot('S29', 6)
plot('S29', 0)

#Figure S30
plot('S30', 1)
plot('S30', 2)
plot('S30', 3)
plot('S30', 4)
plot('S30', 5)
plot('S30', 6)
plot('S30', 0)

#Figure S31
plot('S31', 1)
plot('S31', 2)
plot('S31', 3)
plot('S31', 4)
plot('S31', 5)
plot('S31', 6)
plot('S31', 0)

#Figure S32
plot('S32', 1)
plot('S32', 2)
plot('S32', 3)
plot('S32', 4)
plot('S32', 5)
plot('S32', 6)
plot('S32', 0)

#Figure S33
plot('S33', 1)
plot('S33', 2)
plot('S33', 3)
plot('S33', 4)
plot('S33', 5)
plot('S33', 0)

#Figure S34
plot('S34', 1)
plot('S34', 2)
plot('S34', 3)
plot('S34', 4)
plot('S34', 5)
plot('S34', 6)
plot('S34', 0)

#Figure S35
plot('S35', 1)
plot('S35', 2)
plot('S35', 3)
plot('S35', 4)
plot('S35', 5)
plot('S35', 0)

#Figure S36
plot('S36', 1)
plot('S36', 2)
plot('S36', 3)
plot('S36', 4)
plot('S36', 5)
plot('S36', 6)
plot('S36', 0)

#Figure S37
plot('S37', 1)
plot('S37', 2)
plot('S37', 3)
plot('S37', 4)
plot('S37', 5)
plot('S37', 6)
plot('S37', 0)

#Figure S38
plot('S38', 1)
plot('S38', 2)
plot('S38', 3)
plot('S38', 4)
plot('S38', 5)
plot('S38', 6)
plot('S38', 0)

#Figure S39
plot('S39', 1)
plot('S39', 2)
plot('S39', 3)
plot('S39', 4)
plot('S39', 5)
plot('S39', 6)
plot('S39', 0)

#Figure S40
plot('S40', 1)
plot('S40', 2)
plot('S40', 3)
plot('S40', 4)
plot('S40', 5)
plot('S40', 6)
plot('S40', 0)

#Figure S41
plot('S41', 1)
plot('S41', 2)
plot('S41', 3)
plot('S41', 4)
plot('S41', 5)
plot('S41', 6)
plot('S41', 0)

#Figure S42
plot('S42', 1)
plot('S42', 2)
plot('S42', 3)
plot('S42', 4)
plot('S42', 5)
plot('S42', 6)
plot('S42', 0)

#Figure S43
plot('S43', 1)
plot('S43', 2)
plot('S43', 3)
plot('S43', 4)
plot('S43', 5)
plot('S43', 6)
plot('S43', 0)

#Figure S44
plot('S44', 1)
plot('S44', 2)
plot('S44', 3)
plot('S44', 4)
plot('S44', 5)
plot('S44', 6)
plot('S44', 0)

#Figure S45
plot('S45', 1)
plot('S45', 2)
plot('S45', 3)
plot('S45', 4)
plot('S45', 5)
plot('S45', 6)
plot('S45', 0)