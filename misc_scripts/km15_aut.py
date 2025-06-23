import sys
sys.stdout.reconfigure(line_buffering=True)
import os
import subprocess
import numpy as np
import pandas as pd
import concurrent.futures

import gepard
from gepard.fits import th_KM15

# tranverse target spin asymmetry (TTSA) predictions using km15 model
def km15_aut(file, tpos, xB, Q2, numBins, beamE=10.604):
    
    t_km15 = -abs(tpos) # converting to negative for model
    inputPhiList = [20., 60., 100., 140., 180., 220., 260., 300., 340.]
    phiList = []
    autList = []

    for i in range(numBins):
        phi_rad = inputPhiList[i]*np.pi/180.
        phi_trento = np.pi - phi_rad

        pt = gepard.DataPoint( # defining a DataPoint
            xB                      = xB,
            t                       = t_km15,
            Q2                      = Q2,
            phi                     = phi_trento,
            observable              = 'AUT',
            frame                   = 'Trento',
            process                 = 'ep2epgamma',
            exptype                 = 'fixed target',
            in1energy               = beamE, # defaults to 10.604 if not given
            in1charge               = -1,
            in1polarization         = 0, # I think this should be zero because it is unpolarized, but maybe it's not? because isn't the beam to Hall B always polarized
            in2polarization         = -1,
            in2polarizationvector   = 'T',
            varphi                  = np.pi/2 - phi_trento # I **think** this is referring to the variable in its non-Trento config but I honestly don't know and I'm not sure why it can't compute it automatically
        )
        pt.prepare()
        phiList.append(phi_rad*180./np.pi)
        autList.append(float(th_KM15.predict(pt)))

    #file = "./output_-t_" + str(tpos) + "_xB_" + str(xB) + "_Q2_" + str(Q2) + ".txt"
    #print(file)
    with open(file, "w") as f:
        for i in range(len(phiList)):
            if i != (len(phiList) - 1):
                f.write(str(phiList[i]) + " " + str(autList[i]) + '\n')
                print(str(phiList[i]) + " " + str(autList[i]))
            else:
                f.write(str(phiList[i]) + " " + str(autList[i]))
                print(str(phiList[i]) + " " + str(autList[i]))

    f.close()

#km15_aut(0.198, 0.315, 4.433, 9, 10.604)
            
if __name__ == "__main__":
    method = sys.argv[1]
    a = sys.argv[2]
    b = sys.argv[3]
    c = sys.argv[4]
    d = sys.argv[5]
    e = sys.argv[6]
    f = sys.argv[7]

    if method == "km15_aut":
        km15_aut(str(a), float(b), float(c), float(d), int(e), float(f))
    else:
        print(":(")