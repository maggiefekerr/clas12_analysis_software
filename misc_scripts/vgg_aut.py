import sys
sys.stdout.reconfigure(line_buffering=True)
import os
import subprocess
import concurrent.futures
import numpy as np

# pointing at print-mode dvcsgen installation
os.environ['PATH'] = '/u/home/mkerr/dvcsgens/dvcsgen_print:' + os.environ.get('PATH','')
os.environ['CLASDVCS_PDF'] = '/u/home/mkerr/dvcsgens/dvcsgen_print'

def vgg_aut(file, xB, Q2, tpos, numBins, beamE=10.604):
    # calls dvcsgen in print mode and writes a file containing the phi values, pos & neg cross sections
    inputPhiList = [20., 60., 100., 140., 180., 220., 260., 300., 340.]
    phiList = []
    xsPosList = []
    xsNegList = []

    for i in range(numBins):
        phi_rad = inputPhiList[i]*np.pi/180.
        cmd = [
            'dvcsgen',
            '--beam', f'{beamE:.3f}',
            '--x',    str(xB), str(xB),
            '--q2',   str(Q2), str(Q2),
            '--t',    str(tpos), str(tpos),
            '--phi',  f'{phi_rad:.6f}',
            '--bh',   '3',
            '--gpd',  '101',
            '--ycol', '0.0001'
        ]
        cmd.append('--tpol')
        cmd.append('--globalfit')
        #print(cmd)
        proc = subprocess.run(cmd,
                            stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE,
                            text=True)
        if proc.returncode != 0:
            raise RuntimeError(f"dvcsgen failed:\n{proc.stderr}")
        lines = proc.stdout.splitlines()
        numeric = [ln for ln in lines if ln.strip()]
    
        print(numeric)
        phiList.append(inputPhiList[i])
        xsPosList.append((float)(numeric[-3]))
        xsNegList.append((float)(numeric[-2]))    

    with open(file, "w") as f:
        for i in range(numBins):
            if i != (numBins - 1):
                f.write(str(phiList[i]) + " " + str(xsPosList[i]) + " " + str(xsNegList[i]) + '\n')
            else:
                f.write(str(phiList[i]) + " " + str(xsPosList[i]) + " " + str(xsNegList[i]))
    f.close()

#vgg_aut("test.txt", 0.126, 1.759, 0.670, 9, 10.604)

if __name__ == "__main__":
    method = sys.argv[1]
    a = sys.argv[2]
    b = sys.argv[3]
    c = sys.argv[4]
    d = sys.argv[5]
    e = sys.argv[6]
    f = sys.argv[7]

    if method == "vgg_aut":
        vgg_aut(str(a), float(b), float(c), float(d), int(e), float(f))
    else:
        print(":(")