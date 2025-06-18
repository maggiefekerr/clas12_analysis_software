# practice using gepard from Timothy's code

# useful to see what imports I need to get Python to run on my computer
import sys
sys.stdout.reconfigure(line_buffering=True)
import os
import subprocess
import numpy as np
import pandas as pd
import concurrent.futures

# gepard imports
import gepard as g
from gepard.fits import th_KM15

# inputs: Bjorken x, Q2, t_pos (positive value of t, which is actually negative), phi (in degrees), beam Energy, observable
def km15_model(xB, Q2, t_pos, phi_deg, beam_E=10.604, obs='XS'):

    t_km15 = -abs(t_pos) # converting to negative for model
    phi_rad = np.radians(phi_deg) # converting to radians
    phi_trento = np.pi - phi_rad # converting to the phi they use in model

    pt = g.DataPoint( # defining a DataPoint
        xB        = xB,
        t         = t_km15,
        Q2        = Q2,
        phi       = phi_trento,
        observable= obs, # defaults to cross section if not given
        frame     = 'trento',
        process   = 'ep2epgamma',
        exptype   = 'fixed target',
        in1energy = beam_E, # defaults to 10.604 if not given
        in1charge = -1,
        in1polarization = +1
    )
    print(th_KM15.predict(pt))
    pt.prepare()
    return th_KM15.predict(pt)

def km15_model_inb(xB, Q2, t_pos, phi, beam_E=10.604, obs='XS'):

    t_km15 = -abs(t_pos) # converting to negative for model
    phi_rad = np.radians(phi_deg) # converting to radians
    phi_trento = np.pi - phi_rad # converting to the phi they use in model

    pt = g.DataPoint( # defining a DataPoint
        xB        = xB,
        t         = t_km15,
        Q2        = Q2,
        phi       = phi_trento,
        observable= obs, # defaults to cross section if not given
        frame     = 'trento',
        process   = 'ep2epgamma',
        exptype   = 'fixed target',
        in1energy = beam_E, # defaults to 10.604 if not given
        in1charge = -1,
        in1polarization = -1
    )
    pt.prepare()
    return th_KM15.predict(pt)

def km15_model_outb(xB, Q2, t_pos, phi, beam_E=10.604, obs='XS'):

    t_km15 = -abs(t_pos) # converting to negative for model
    phi_rad = np.radians(phi_deg) # converting to radians
    phi_trento = np.pi - phi_rad # converting to the phi they use in model

    pt = g.DataPoint( # defining a DataPoint
        xB        = xB,
        t         = t_km15,
        Q2        = Q2,
        phi       = phi_trento,
        observable= 'ALU', # defaults to cross section if not given
        frame     = 'trento',
        process   = 'ep2epgamma',
        exptype   = 'fixed target',
        in1energy = beam_E, # defaults to 10.604 if not given
        in1charge = -1,
        in1polarization = +1
    )
    pt.prepare()
    print("hello")
    return th_KM15.predict(pt)

if __name__ == "__main__":
    method = sys.argv[1]
    a = sys.argv[2]
    b = sys.argv[3]
    c = sys.argv[4]
    d = sys.argv[5]
    e = sys.argv[6]
    f = sys.argv[7]

    if method == "km15_model":
        #km15_model(float(a), float(b), float(c))
        km15_model(float(a), float(b), float(c), float(d), float(e), str(f)) # at this stage we have to specify the observable because I don't know how to handle it otherwise actually I do I'm just too lazy

    elif method == "km15_model_inb":
        km15_model_inb(float(a), float(b), float(c), float(d), str(e))

    elif method == "km15_model_outb":
        km15_model_outb(float(a), float(b), float(c), float(d), str(e))

    else:
        print(":(")