import os
import subprocess

# pointing at print-mode dvcsgen installation
os.environ['PATH'] = '/u/home/mkerr/dvcsgens/dvcsgen_print:' + os.environ.get('PATH','')
os.environ['CLASDVCS_PDF'] = '/u/home/mkerr/dvcsgens/dvcsgen_print'

# returns polarized cross section for +'ve bea, polarization
def vgg_model(xB, Q2, tpos, numBins, numDiv, beamE=10.604, bh=3, gpd=101, globalfit=True):
    # calls dvcsgen in print mode and writes a file containing the phi values, pos & neg cross sections
    phiList = []
    xsPosList = []
    xsNegList = []

    for i in range(numBins):
        for j in range(numDiv):
            phi_deg = (2.0*i +1.0)*0.5*(360.0/numBins) - 0.5*360.0/numBins + j*360/(numDiv*numBins) + 360/(2*numDiv*numBins)
            #(i*360/numBins) + (360.0/numBins)*(j/numDiv)
            cmd = [
                'dvcsgen',
                '--beam', f'{beamE:.3f}',
                '--x',    str(xB), str(xB),
                '--q2',   str(Q2), str(Q2),
                '--t',    str(tpos), str(tpos),
                '--phi',  f'{phi_deg:.6f}',
                '--bh',   str(bh),
                '--gpd',  str(gpd),
                '--ycol', '0.0001'
            ]
            if globalfit:
                cmd.append('--globalfit')
            proc = subprocess.run(cmd,
                                stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE,
                                text=True)
            if proc.returncode != 0:
                raise RuntimeError(f"dvcsgen failed:\n{proc.stderr}")
            lines = proc.stdout.splitlines()
            numeric = [ln for ln in lines if ln.strip()]
                
            phiList.append(phi_deg)
            xsPosList.append((float)(numeric[-2]))
            xsNegList.append((float)(numeric[-3]))    

    with open("vgg_xs_phi-xsPos-xsNeg.txt", "w") as f:
        for i in range(len(phiList)):
            if i != (len(phiList) - 1):
                f.write(str(phiList[i]) + " " + str(xsPosList[i]) + " " + str(xsNegList[i]) + '\n')
            else:
                print(str(phiList[i]) + " " + str(xsPosList[i]) + " " + str(xsNegList[i]))
                f.write(str(phiList[i]) + " " + str(xsPosList[i]) + " " + str(xsNegList[i]))
    f.close()

#vgg_model(0.126, 1.759, 0.670, 9, 10, 10.604)