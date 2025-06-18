import os
import subprocess

# pointing at print-mode dvcsgen installation
os.environ['PATH'] = 'u/home/mkerr/dvcsgens/dvcsgen_print:' + os.environ.get('PATH','')
os.environ['CLASDVCS_PDF'] = '/u/home/mkerr/dvcsgens/dvcsgen_print'

# returns polarized cross section for +'ve bea, polarization
def vgg_model_xs_pos(xB, Q2, t_pos, phi_rad, beam_E=10.604, bh=3, gpd=101, globalfit=True):
    # calls dvcsgen in print mode and returns positively polarized cross section
    cmd = ['dvcsgen',
            '--beam', f'{beam:.6f}',
            '--x', str(x), str(x),
            '--q2', str(Q2), str(Q2),
            '--t', str(t_pos), str(t_pos),
            '--phi', f'{phi_rad:.6f}',
            '--bh', str(bh),
            '--gpd', str(gpd),
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

    return float(numeric[-2])


# returns polarized cross section for -'ve beam polarization
def vgg_model_xs_neg(xB, Q2, t_pos, phi_rad, beam_E=10.604, bh=3, gpd=101, globalfit=True):
    # calls dvcsgen in print mode and returns negatively polarized cross section

    cmd = ['dvcsgen',
            '--beam', f'{beam:.3f}',
            '--x', str(x), str(x),
            '--q2', str(Q2), str(Q2),
            '--t', str(t_pos), str(t_pos),
            '--phi', f'{phi_rad:.6f}',
            '--bh', str(bh),
            '--gpd', str(gpd),
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