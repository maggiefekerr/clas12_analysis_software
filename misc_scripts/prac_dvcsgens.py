import os
import subprocess

# pointing at print-mode dvcsgen installation
os.environ['PATH'] = '/u/home/mkerr/dvcsgens/dvcsgen_print:' + os.environ.get('PATH','')
os.environ['CLASDVCS_PDF'] = '/u/home/mkerr/dvcsgens/dvcsgen_print'

# returns polarized cross section for +'ve bea, polarization
def vgg_model_xs_pos(beam, x, Q2, t, phi_rad,
                                    bh=3, gpd=101, globalfit=True):
    cmd = [
        'dvcsgen',
        '--beam', f'{beam:.3f}',
        '--x',    str(x), str(x),
        '--q2',   str(Q2), str(Q2),
        '--t',    str(t), str(t),
        '--phi',  f'{phi_rad:.6f}',
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

    #sigma_minus = float(numeric[-3])
    sigma_plus  = float(numeric[-2])
    return sigma_plus#, sigma_minus


# returns polarized cross section for -'ve beam polarization
def vgg_model_xs_neg(x, Q2, t, phi_rad, beam=10.604, bh=3, gpd=101, globalfit=True):
    # calls dvcsgen in print mode and returns negatively polarized cross section

    cmd = [
        'dvcsgen',
        '--beam', f'{beam:.3f}',
        '--x',    str(x), str(x),
        '--q2',   str(Q2), str(Q2),
        '--t',    str(t), str(t),
        '--phi',  f'{phi_rad:.6f}',
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

    return float(numeric[-3])

a= vgg_model_xs_pos(
    beam    = 10.6,
    x       = 0.126,
    Q2      = 1.759,
    t       = 0.670,
    phi_rad = 90.0
)
print("sigma(+) =", a)
#print("sigma(−) =", b)