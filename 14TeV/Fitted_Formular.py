#!bin/env python

def getWKFactor(x):
    if x>7000:
        x=7000
    return 1.12231+0.000396307*x-3.21092e-07*(x**2)+7.03988e-11*(x**3)-4.93892e-15*(x**4) 
