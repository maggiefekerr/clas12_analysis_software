import os
import json
import math
import numpy as np
from collections import defaultdict

def load_contamination(period, contamination_dir="contamination"):
    file_path = os.path.join(contamination_dir, f"contamination_{period}.json")
    with open(file_path) as f:
        # Handle both tuple formats: "(0, 0, 0, 0)" and "0,0,0,0"
        return {tuple(map(int, k.strip("()").replace(" ", "").split(","))): v 
                for k, v in json.load(f).items()}

def load_bsa(period, channel, bsa_dir="bsa_results"):
    file_path = os.path.join(bsa_dir, f"raw_bsa_{channel}_{period}.json")
    try:
        with open(file_path) as f:
            return {tuple(map(int, k.strip("()").replace(" ", "").split(","))): v 
                    for k, v in json.load(f).items() 
                    if v.get("valid", True)}
    except:
        return {}

def propagate_errors(D, D_err, c, c_err, E, E_err):
    denominator = 1 - c
    # Check if the denominator is zero or nearly zero
    if abs(denominator) < 1e-12:  # tolerance to catch floating-point issues
        return 0, 1  # fail-safe: return 0 ± 1
    A = (D - c * E) / denominator
    dA_dD = 1 / denominator
    dA_dc = (D - E) / (denominator ** 2)
    dA_dE = -c / denominator
    variance = (dA_dD * D_err) ** 2 + (dA_dc * c_err) ** 2 + (dA_dE * E_err) ** 2
    A_err = math.sqrt(variance)
    return A, A_err

def process_period(period, contamination_dir, bsa_dir):
    # For DVCS periods, load corresponding eppi0 period data
    eppi0_period = period.replace("DVCS", "eppi0") if "DVCS" in period else period
    
    contamination = load_contamination(period, contamination_dir)
    dvcs_bsa = load_bsa(period, "dvcs", bsa_dir)
    eppi0_bsa = load_bsa(eppi0_period, "eppi0", bsa_dir)  # Changed here
    
    print(f"\nProcessing {period}:")
    print(f"Contamination bins: {len(contamination)}")
    print(f"DVCS BSA bins: {len(dvcs_bsa)}")
    print(f"EPPI0 BSA bins ({eppi0_period}): {len(eppi0_bsa)}")
    
    common_bins = set(dvcs_bsa.keys()) & set(eppi0_bsa.keys()) & set(contamination.keys())
    print(f"Common bins: {len(common_bins)}")

     # Temporary debug: compare first 5 bins
    print("\nSample DVCS BSA bins:", list(dvcs_bsa.keys())[:5])
    print("Sample EPPI0 BSA bins:", list(eppi0_bsa.keys())[:5])
    print("Sample Contamination bins:", list(contamination.keys())[:5])
    
    results = {}
    for bin_key in set(dvcs_bsa.keys()) & set(eppi0_bsa.keys()) & set(contamination.keys()):
        try:
            D = dvcs_bsa[bin_key]["bsa"]
            D_err = dvcs_bsa[bin_key]["bsa_err"]
            c = contamination[bin_key]["c_i"]
            c_err = contamination[bin_key]["c_i_err"]
            E = eppi0_bsa[bin_key]["bsa"]
            E_err = eppi0_bsa[bin_key]["bsa_err"]
            
            if abs(c) < 1e-6:
                A, A_err = D, D_err
            else:
                A, A_err = propagate_errors(D, D_err, c, c_err, E, E_err)
            
            results[bin_key] = {
                "bsa": round(A, 5),
                "bsa_err": round(A_err, 5),
                "valid": True,
                "components": {
                    "D": round(D, 5),
                    "D_err": round(D_err, 5),
                    "c": round(c, 5),
                    "c_err": round(c_err, 5),
                    "E": round(E, 5),
                    "E_err": round(E_err, 5)
                }
            }
        except KeyError:
            continue
    return results

def combine_periods(periods, final_dir):
    all_data = defaultdict(list)
    for period in periods:
        try:
            with open(os.path.join(final_dir, f"adjusted_bsa_{period}.json")) as f:
                data = {tuple(map(int, k.strip("()").split(", "))): v 
                       for k, v in json.load(f).items()}
            for bin_key, values in data.items():
                if values.get("valid", True):
                    all_data[bin_key].append(values)
        except:
            continue

    combined = {}
    for bin_key, measurements in all_data.items():
        if len(measurements) < 1: continue
        
        weights = [1/(m["bsa_err"]**2 + 1e-9) for m in measurements]
        total_weight = sum(weights)
        
        if total_weight > 0:
            avg_bsa = sum(w*m["bsa"] for w, m in zip(weights, measurements))/total_weight
            avg_err = 1/math.sqrt(total_weight)
        else:
            bsa_vals = [m["bsa"] for m in measurements]
            avg_bsa = np.mean(bsa_vals)
            avg_err = np.std(bsa_vals) if len(bsa_vals) > 1 else 0
        
        combined[bin_key] = {
            "bsa": round(avg_bsa, 5),
            "bsa_err": round(avg_err, 5),
            "n_periods": len(measurements),
            "valid": True
        }
    return combined

def determine_final_bsa(contamination_dir="contamination", bsa_dir="bsa_results", final_dir="final_results"):
    os.makedirs(final_dir, exist_ok=True)
    periods = ["DVCS_Fa18_inb", "DVCS_Fa18_out", "DVCS_Sp19_inb", "DVCS_Sp18_inb", "DVCS_Sp18_out"]
    
    for period in periods:
        try:
            result = process_period(period, contamination_dir, bsa_dir)
            print(f"{period} processed bins: {len(result)}")
            with open(os.path.join(final_dir, f"adjusted_bsa_{period}.json"), "w") as f:
                json.dump({str(k): v for k, v in result.items()}, f, indent=2)
        except Exception as e:
            print(f"Failed processing {period}: {str(e)}")
    
    # Combine periods
    combined = combine_periods(periods, final_dir)
    with open(os.path.join(final_dir, "combined_bsa.json"), "w") as f:
        json.dump({str(k): v for k, v in combined.items()}, f, indent=2)
    return combined