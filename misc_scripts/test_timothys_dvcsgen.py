import os
import subprocess

os.environ['PATH'] = '/u/home/mkerr/dvcsgens/dvcsgen_print:' + os.environ.get('PATH','')
os.environ['CLASDVCS_PDF'] = '/u/home/mkerr/dvcsgens/dvcsgen_print'

def dvcsgen_polarized_cross_section(beam, x, Q2, t, phi_rad,
                                    bh=3, gpd=101, globalfit=True):
    """
    Calls dvcsgen in print‐mode and returns (sigma_plus, sigma_minus).
    """
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

    sigma_minus = float(numeric[-3])
    sigma_plus  = float(numeric[-2])
    return sigma_plus, sigma_minus
#enddef

a, b = dvcsgen_polarized_cross_section(
    beam    = 10.6,
    x       = 0.126,
    Q2      = 1.759,
    t       = 0.670,
    phi_rad = 90.0
)
print("sigma(+) =", a)
print("sigma(−) =", b)