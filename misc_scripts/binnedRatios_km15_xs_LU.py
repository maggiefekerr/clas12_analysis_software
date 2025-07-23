#binnedRatios_km15_xs_LU.py
import sys
sys.stdout.reconfigure(line_buffering=True)
import os
import subprocess
import numpy as np
import pandas as pd
import concurrent.futures

import gepard
from gepard.fits import th_KM15

# this is specifically generating data that can be called to be used in the code where I am looking at the binned ratios
# could definitely be modified to be used elsewhere but alas

def xs(directory, name, Q2, xB, t, beamE, numBins):
    phiDivs = []
    for i in range(numBins):
        phiDivs.append((i*2*np.pi)/numBins)
    phiDivs.append(2*np.pi)

    phiList = []
    xsList = []
    
    t_km15 = -abs(t) # ensuring we are using actual t value (negative) for model

    for i in range(numBins):
        phi_rad = 0.
        if (i != numBins - 1):
            phi_rad = (0.5*(phiDivs[i] + phiDivs[i+1]))
            phiList.append(phi_rad)
        else:
            phi_rad = (0.5*(phiDivs[i] + 2*np.pi))
            phiList.append(phi_rad)
        phi_trento = np.pi - phi_rad

        pt = gepard.DataPoint( # defining a DataPoint
            xB                      = xB,
            t                       = t_km15,
            Q2                      = Q2,
            phi                     = phi_trento,
            observable              = 'XS',
            frame                   = 'Trento',
            process                 = 'ep2epgamma',
            exptype                 = 'fixed target',
            in1energy               = beamE,
            in1charge               = -1,
            in1polarization         = +1,
            in2polarizationvector   = "U"
        )

        pt.prepare()
        xsList.append(float(th_KM15.predict(pt)))

   file = "./ratios/data/" + directory + "/xs/-t_" + str(t) + "/" + name + ".txt"

    with open(file, "w") as f:
        for i in range(numBins):
            if i != (numBins - 1):
                f.write(str(phiList[i]) + " " + str(xsList[i]) + '\n')
            else:
                f.write(str(phiList[i]) + " " + str(xsList[i]))
    f.close()

if __name__ == "__main__":
    method = sys.argv[1]
    a = sys.argv[2]
    b = sys.argv[3]
    c = sys.argv[4]
    d = sys.argv[5]
    e = sys.argv[6]
    f = sys.argv[7]
    g = sys.argv[8]

    if method == "xs":
        xs(str(a), str(b), float(c), float(d), float(e), float(f), int(g))
    else:
        print(":(")