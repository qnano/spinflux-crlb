# -*- coding: utf-8 -*-
"""
Simulates all SpinFlux Cramér-Rao lower bound data underlying figures.
    - simulate(plotnumber) simulates the SpinFlux Cramér-Rao lower bound (CRLB) 
    data underlying the figure specified in plotnumber.
    
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
"""

from spinflux_simulation_utils import simulate

simulate('SMLM')
simulate('3-S4')
simulate('3-heatmap')
simulate('4de-S6')
simulate('4bc-heatmap')
simulate('4ij-S10')
simulate('4gh-heatmap')
simulate('5-S16')
simulate('5-heatmap')
simulate('S2')
simulate('S5')
simulate('S7')
simulate('S8')
simulate('S9')
simulate('S11')
simulate('S12')
simulate('S13')
simulate('S14')
simulate('S15')
simulate('S17')
simulate('S18')
simulate('S19')
simulate('S20')
simulate('S21')
simulate('S22')
simulate('S23')
simulate('S24')
simulate('S25')
simulate('S26')
simulate('S27')
simulate('S28')
simulate('S29')
simulate('S30')
simulate('S31')
simulate('S32')
simulate('S33')
simulate('S34')
simulate('S35')
simulate('S36')
simulate('S37')
simulate('S38')
simulate('S39')
simulate('S40')
simulate('S41')
simulate('S42')
simulate('S43')
simulate('S44')
simulate('S45')