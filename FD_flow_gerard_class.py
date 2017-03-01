import numpy as np
import scipy as sp
import matploblib.pyplot as plt

class glaciermodel:
    ## --------------------------------
    #define parameters and constants
    #----------------------------------
    rho=917              # density of ice in kg/m^3
    g=9.8                  # gravity in m/s^2
    n=3                  # glen's flow law constant
    
    s_per_year = 365.35*24*3600        # number of seconds in a year (s yr^-1)
    A_T = 6.8e-24*s_per_year       # softness parameter for -5 degree ice (yr^-1 Pa^-3)
    fd = 2*A_T/(n+2)                           # constants lumped together, and put into standard units.  (Pa^-3 yr^-1)
    fs = 5.7e-20*s_per_year                # sliding parameter fromOerlemans; (Pa^-3 m^2 yr-1)
    
    def __init__(self):
        pass
        
    def 