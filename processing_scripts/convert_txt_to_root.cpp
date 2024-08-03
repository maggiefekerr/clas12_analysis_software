// Include necessary ROOT headers
#include <TFile.h>
#include <TTree.h>
#include <iostream>
#include <fstream>
#include <cmath>
#include <sstream>
#include <sys/types.h> // Include sys/types.h for stat
#include <sys/stat.h>  // Include sys/stat.h for stat structure

using namespace std;

// Function to find the root directory of the "clas12_analysis_software" repository
std::string findPackageRoot() {
    char cwd[4096];  // Use a fixed size buffer if PATH_MAX is not defined
    if (getcwd(cwd, sizeof(cwd)) != NULL) {
        std::string currentDir(cwd);
        std::string searchDir = "clas12_analysis_software";

        // Check if the current directory or its parent directories contain "clas12_analysis_software"
        size_t found = currentDir.find(searchDir);
        if (found != std::string::npos) {
            return currentDir.substr(0, found + searchDir.length()) + "/";
        } else {
            std::cerr << "Warning: Unable to locate the 'clas12_analysis_software' directory.\n"
                      << "Using relative path for CSV file.\n";
            return "./";  // Default to current directory
        }
    } else {
        std::cerr << "Error: Unable to get current working directory.\n";
        return "./";
    }
}

// Function to calculate the Cartesian components from spherical coordinates
void sphericalToCartesian(double p, double phi, double theta, double& px, double& py, double& pz) {
    px = p * sin(theta) * cos(phi);
    py = p * sin(theta) * sin(phi);
    pz = p * cos(theta);
}

// Function to calculate the magnitude of a vector from its Cartesian components
double vectorMagnitude(double px, double py, double pz) {
    return sqrt(px * px + py * py + pz * pz);
}

// Function to calculate the polar angle from the Cartesian components
double vectorPolarAngle(double pz, double vectorMag) {
    return acos(pz / vectorMag);
}

struct RunInfo {
      int runnum;
      float total_charge;
      float positive_charge;
      float negative_charge;
      float target_polarization;
      float target_polarization_uncertainty;
};

// Declare a vector to store the run information
std::vector<RunInfo> run_info_list;

void load_run_info_from_csv(const std::string& filename) {
   // Open the input file with the given filename
  std::ifstream file(filename);
  if (!file) {
    std::cerr << "Failed to open file: " << filename << std::endl;
    return;
  }

  // Declare a string to store each line read from the file
  std::string line;

  // Loop through each line in the file until there are no more lines left to read
  while (std::getline(file, line)) {
    // If the line is empty or starts with a '#' (comment), skip to the next line
    if (line.empty() || line[0] == '#') { continue; }

    // Use a stringstream to split the line by commas
    std::stringstream ss(line);

    // Declare a struct to store the run information
    RunInfo run_info;

    // Declare a string to store each piece of information read from the stringstream
    std::string info;

    // Read the run number from the stringstream and convert it to an integer
    std::getline(ss, info, ',');
    run_info.runnum = std::stoi(info);

    // Read the total charge from the stringstream and convert it to a float
    std::getline(ss, info, ',');
    run_info.total_charge = std::stof(info);

    // Read the positive charge from the stringstream and convert it to a float
    std::getline(ss, info, ',');
    run_info.positive_charge = std::stof(info);

    // Read the negative charge from the stringstream and convert it to a float
    std::getline(ss, info, ',');
    run_info.negative_charge = std::stof(info);

    // Read the target polarization from the stringstream and convert it to a float
    std::getline(ss, info, ',');
    run_info.target_polarization = std::stof(info);

    // Read the target polarization from the stringstream and convert it to a float
    std::getline(ss, info, ',');
    run_info.target_polarization_uncertainty = std::stof(info);

    // Add the struct to the run_info_list vector
    run_info_list.push_back(run_info);
  }
}

// function to get tmin 
double gettmin(double x) {
    double mp = 0.938272; // proton mass in GeV
    return -pow((mp*x),2)/(1-x);
}

// function get t
double gett(double p, double theta) {
    double mp = 0.938272; // proton mass in GeV
    double E = mp; // target proton energy (written as E to help checking calculation but at rest)
    return 2*mp*(E - p) - 2*sqrt(mp*mp + E*E)*sqrt(mp*mp + p*p) +
          2*sqrt(mp*mp + E*E)*sqrt(mp*mp + p*p)*cos(theta);
}

// function to get the polarization value
double getPol(int runnum) {
  double pol = 0.86; 
    if (runnum == 11 ) { pol = 0.86; } // runnum == 11 indicates Monte Carlo in CLAS12
    else if (runnum >= 5032 && runnum < 5333) { pol = 0.8592; } 
    else if (runnum >= 5333 && runnum <= 5666) { pol = 0.8922; }
    else if (runnum >= 6616 && runnum <= 6783) { pol = 0.8453; }
    else if (runnum >= 6142 && runnum <= 6149) { pol = 0.81132; }
    else if (runnum >= 6150 && runnum <= 6188) { pol = 0.82137; }
    else if (runnum >= 6189 && runnum <= 6260) { pol = 0.83598; }
    else if (runnum >= 6261 && runnum <= 6339) { pol = 0.80770; }
    else if (runnum >= 6340 && runnum <= 6342) { pol = 0.85536; }
    else if (runnum >= 6344 && runnum <= 6399) { pol = 0.87038; }
    else if (runnum >= 6420 && runnum <= 6476) { pol = 0.88214; }
    else if (runnum >= 6479 && runnum <= 6532) { pol = 0.86580; }
    else if (runnum >= 6533 && runnum <= 6603) { pol = 0.87887; }
    else if (runnum >= 11013 && runnum <= 11309) { pol = 0.84983; }
    else if (runnum >= 11323 && runnum <= 11334) { pol = 0.87135; }
    else if (runnum >= 11335 && runnum <= 11387) { pol = 0.85048; }
    else if (runnum >= 11389 && runnum <= 11571) { pol = 0.84262; }
    else if (runnum >= 16042 && runnum <= 17811) { pol = 0.827; } // RGC +/- 0.01440
  return pol;
}

int main(int argc, char *argv[]) {
    // Check for correct number of command line arguments
    if (argc != 5) {
        cout << "Usage: " << argv[0];
        cout << " <input_text_file> <output_root_file> <hadron_count> <is_mc>" << endl;
        cout << " <hadron_count> = 0 for inclusive" << endl;
        cout << " <hadron_count> = 1 for single hadron" << endl;
        cout << " <hadron_count> = 2 for dihadron" << endl;
        cout << " <hadron_count> = 3 for trihadron" << endl;
        cout << " <hadron_count> = 4 for dvcs (epgamma)" << endl;
        cout << " <hadron_count> = 5 for calibration" << endl;
        cout << "Should rename <hadron_count> at some point." << endl;
        return 1;
    }
    
    // Open input text file for reading
    ifstream infile(argv[1]);
    if (!infile.is_open()) {
        cout << "Error opening input text file." << endl;
        return 1;
    }

    // Create a new ROOT file and TTree for output
    TFile *outfile = new TFile(argv[2], "RECREATE");
    TTree *tree = new TTree("PhysicsEvents", "Physics Events Tree");

    // Determine the hadron count from the command line argument
    int hadron_count = atoi(argv[3]);
    int is_mc = atoi(argv[4]);

    // Declare common variables
    int runnum, evnum, helicity;
    double beam_pol, target_pol, e_p, e_theta, e_phi, vz_e, Q2, W, Mx, Mx2, x, y;
    double t, tmin;
    double z, xF, pT, zeta, eta, phi, DepA, DepB, DepC, DepV, DepW;
    double p_p, p_theta, p_phi, vz_p;
    // Additional variables for one or two hadrons
    double p1_p, p1_theta, p1_phi, vz_p1, p2_p, p2_theta, p2_phi, vz_p2;
    double z1, z2, Mh, xF1, xF2, pT1, pT2, pTpT, zeta1, zeta2;
    double t1, t1min, t2, t2min, Mx1;
    double eta1, eta2, Delta_eta, eta1_gN, eta2_gN;
    double phi1, phi2, Delta_phi, phih, phiR, theta;
    double Emiss2, theta_gamma_gamma, pTmiss, Mxgammasquared;
    // Additional variables for three hadrons
    double p3_p, p3_theta, p3_phi, vz_p3;
    double z3, z12, z13, z23, Mh12, Mh13, Mh23, xF3, xF12, xF13, xF23;
    double t3, t3min, Mx3, Mx12, Mx13, Mx23;
    double pT3, pT12, pT13, pT23, zeta3, zeta12, zeta23, zeta13;
    double eta3, eta12, eta23, eta13;
    double phi3, phi12, phi13, phi23, Delta_phi12, Delta_phi13, Delta_phi23;

    // Additional variables for mc single hadron
    double mc_e_p, mc_e_theta, mc_e_phi, mc_vz_e, mc_p_p, mc_p_theta, mc_p_phi, mc_vz_p;
    double mc_Q2, mc_W, mc_Mx, mc_Mx2, mc_x, mc_y, mc_t, mc_tmin;
    double mc_z, mc_xF, mc_pT, mc_zeta, mc_eta, mc_phi;
    double mc_DepA, mc_DepB, mc_DepC, mc_DepV, mc_DepW;
    int matching_e_pid, matching_p1_pid, mc_p1_parent;
    // Additional variables for mc dihadron 
    double mc_p1_p, mc_p1_theta, mc_p1_phi, mc_vz_p1, mc_p2_p, mc_p2_theta, mc_p2_phi, mc_vz_p2;
    double mc_z1, mc_z2, mc_Mh, mc_xF1, mc_xF2, mc_pT1, mc_pT2, mc_pTpT, mc_zeta1, mc_zeta2;
    double mc_t1, mc_t1min, mc_t2, mc_t2min, mc_Mx1;
    double mc_eta1, mc_eta2, mc_Delta_eta, mc_eta1_gN, mc_eta2_gN;
    double mc_phi1, mc_phi2, mc_Delta_phi, mc_phih, mc_phiR, mc_theta;
    int matching_p2_pid, mc_p2_parent;

    // Additional variables for calibration scripts
    int config_run, config_event, event_helicity;
    int particle_pid, particle_status;
    double particle_px, particle_py, particle_pz, p; // theta and phi already defined
    double particle_vx, particle_vy, particle_vz, particle_beta, particle_chi2pid;
    int cal_sector;
    double cal_energy_1, cal_x_1, cal_y_1, cal_z_1, cal_lu_1, cal_lv_1, cal_lw_1;
    double cal_energy_4, cal_x_4, cal_y_4, cal_z_4, cal_lu_4, cal_lv_4, cal_lw_4;
    double cal_energy_7, cal_x_7, cal_y_7, cal_z_7, cal_lu_7, cal_lv_7, cal_lw_7;
    double cc_nphe_15, cc_nphe_16;
    int track_sector_5, track_ndf_5, track_sector_6, track_ndf_6;
    double track_chi2_5, track_chi2_6;
    double traj_x_6, traj_y_6, traj_z_6, traj_edge_6;
    double traj_x_18, traj_y_18, traj_z_18, traj_edge_18;
    double traj_x_36, traj_y_36, traj_z_36, traj_edge_36;
    double traj_x_1, traj_y_1, traj_z_1, traj_edge_1;
    double traj_x_3, traj_y_3, traj_z_3, traj_edge_3;
    double traj_x_5, traj_y_5, traj_z_5, traj_edge_5;
    double traj_x_7, traj_y_7, traj_z_7, traj_edge_7;
    double traj_x_12, traj_y_12, traj_z_12, traj_edge_12;
    double ft_energy, ft_x, ft_y, ft_z, ft_radius;

    // Case for zero hadrons (inclusive)
    if (hadron_count == 0 && is_mc == 0) {
        // Link TTree branches to variables for zero hadrons
        tree->Branch("runnum", &runnum, "runnum/I");
        tree->Branch("evnum", &evnum, "evnum/I");
        tree->Branch("helicity", &helicity, "helicity/I");
        tree->Branch("beam_pol", &beam_pol, "beam_pol/D");
        tree->Branch("target_pol", &target_pol, "target_pol/D");
        tree->Branch("e_p", &e_p, "e_p/D");
        tree->Branch("e_theta", &e_theta, "e_theta/D");
        tree->Branch("e_phi", &e_phi, "e_phi/D");
        tree->Branch("vz_e", &vz_e, "vz_e/D");
        tree->Branch("Q2", &Q2, "Q2/D");
        tree->Branch("W", &W, "W/D");
        tree->Branch("Mx", &Mx, "Mx/D");
        tree->Branch("Mx2", &Mx2, "Mx2/D");
        tree->Branch("x", &x, "x/D");
        tree->Branch("y", &y, "y/D");
        tree->Branch("t", &t, "t/D");
        tree->Branch("tmin", &tmin, "tmin/D");
        tree->Branch("DepA", &DepA, "DepA/D");
        tree->Branch("DepB", &DepB, "DepB/D");
        tree->Branch("DepC", &DepC, "DepC/D");
        tree->Branch("DepV", &DepV, "DepV/D");
        tree->Branch("DepW", &DepW, "DepW/D");
    }
    // Case for zero hadrons (inclusive) with mc
    if (hadron_count == 0 && is_mc == 1) {
        // Link TTree branches to variables for zero hadrons
        tree->Branch("e_p", &e_p, "e_p/D");
        tree->Branch("e_theta", &e_theta, "e_theta/D");
        tree->Branch("e_phi", &e_phi, "e_phi/D");
        tree->Branch("vz_e", &vz_e, "vz_e/D");
        tree->Branch("Q2", &Q2, "Q2/D");
        tree->Branch("W", &W, "W/D");
        tree->Branch("Mx", &Mx, "Mx/D");
        tree->Branch("Mx2", &Mx2, "Mx2/D");
        tree->Branch("x", &x, "x/D");
        tree->Branch("y", &y, "y/D");
        tree->Branch("t", &t, "t/D");
        tree->Branch("tmin", &tmin, "tmin/D");
        tree->Branch("DepA", &DepA, "DepA/D");
        tree->Branch("DepB", &DepB, "DepB/D");
        tree->Branch("DepC", &DepC, "DepC/D");
        tree->Branch("DepV", &DepV, "DepV/D");
        tree->Branch("DepW", &DepW, "DepW/D");
        //
        tree->Branch("mc_e_p", &mc_e_p, "mc_e_p/D");
        tree->Branch("mc_e_theta", &mc_e_theta, "mc_e_theta/D");
        tree->Branch("mc_e_phi", &mc_e_phi, "mc_e_phi/D");
        tree->Branch("mc_vz_e", &mc_vz_e, "mc_vz_e/D");
        tree->Branch("mc_Q2", &mc_Q2, "mc_Q2/D");
        tree->Branch("mc_W", &mc_W, "mc_W/D");
        tree->Branch("mc_Mx", &mc_Mx, "mc_Mx/D");
        tree->Branch("mc_Mx2", &mc_Mx2, "mc_Mx2/D");
        tree->Branch("mc_x", &mc_x, "mc_x/D");
        tree->Branch("mc_y", &mc_y, "mc_y/D");
        tree->Branch("mc_t", &mc_t, "mc_t/D");
        tree->Branch("mc_tmin", &mc_tmin, "mc_tmin/D");
        tree->Branch("mc_DepA", &mc_DepA, "mc_DepA/D");
        tree->Branch("mc_DepB", &mc_DepB, "mc_DepB/D");
        tree->Branch("mc_DepC", &mc_DepC, "mc_DepC/D");
        tree->Branch("mc_DepV", &mc_DepV, "mc_DepV/D");
        tree->Branch("mc_DepW", &mc_DepW, "mc_DepW/D");
    }
    // Case for one hadron
    else if (hadron_count == 1 && is_mc == 0) {
        // Link TTree branches to variables for one hadron
        tree->Branch("runnum", &runnum, "runnum/I");
        tree->Branch("evnum", &evnum, "evnum/I");
        tree->Branch("helicity", &helicity, "helicity/I");
        tree->Branch("beam_pol", &beam_pol, "beam_pol/D");
        tree->Branch("target_pol", &target_pol, "target_pol/D");
        tree->Branch("e_p", &e_p, "e_p/D");
        tree->Branch("e_theta", &e_theta, "e_theta/D");
        tree->Branch("e_phi", &e_phi, "e_phi/D");
        tree->Branch("vz_e", &vz_e, "vz_e/D");
        tree->Branch("p_p", &p_p, "p_p/D");
        tree->Branch("p_theta", &p_theta, "p_theta/D");
        tree->Branch("p_phi", &p_phi, "p_phi/D");
        tree->Branch("vz_p", &vz_p, "vz_p/D");
        tree->Branch("Q2", &Q2, "Q2/D");
        tree->Branch("W", &W, "W/D");
        tree->Branch("Mx", &Mx, "Mx/D");
        tree->Branch("Mx2", &Mx2, "Mx2/D");
        tree->Branch("x", &x, "x/D");
        tree->Branch("y", &y, "y/D");
        tree->Branch("t", &t, "t/D");
        tree->Branch("tmin", &tmin, "tmin/D");
        tree->Branch("z", &z, "z/D");
        tree->Branch("xF", &xF, "xF/D");
        tree->Branch("pT", &pT, "pT/D");
        tree->Branch("zeta", &zeta, "zeta/D");
        tree->Branch("eta", &eta, "eta/D");
        tree->Branch("phi", &phi, "phi/D");
        tree->Branch("DepA", &DepA, "DepA/D");
        tree->Branch("DepB", &DepB, "DepB/D");
        tree->Branch("DepC", &DepC, "DepC/D");
        tree->Branch("DepV", &DepV, "DepV/D");
        tree->Branch("DepW", &DepW, "DepW/D");
    }
    // Case for one hadron and is monte carlo
    else if (hadron_count == 1 && is_mc == 1) {
        // Link TTree branches to variables for one hadron
        tree->Branch("e_p", &e_p, "e_p/D");
        tree->Branch("e_theta", &e_theta, "e_theta/D");
        tree->Branch("e_phi", &e_phi, "e_phi/D");
        tree->Branch("vz_e", &vz_e, "vz_e/D");
        tree->Branch("p_p", &p_p, "p_p/D");
        tree->Branch("p_theta", &p_theta, "p_theta/D");
        tree->Branch("p_phi", &p_phi, "p_phi/D");
        tree->Branch("vz_p", &vz_p, "vz_p/D");
        tree->Branch("Q2", &Q2, "Q2/D");
        tree->Branch("W", &W, "W/D");
        tree->Branch("Mx", &Mx, "Mx/D");
        tree->Branch("Mx2", &Mx2, "Mx2/D");
        tree->Branch("x", &x, "x/D");
        tree->Branch("y", &y, "y/D");
        tree->Branch("t", &t, "t/D");
        tree->Branch("tmin", &tmin, "tmin/D");
        tree->Branch("z", &z, "z/D");
        tree->Branch("xF", &xF, "xF/D");
        tree->Branch("pT", &pT, "pT/D");
        tree->Branch("zeta", &zeta, "zeta/D");
        tree->Branch("eta", &eta, "eta/D");
        tree->Branch("phi", &phi, "phi/D");
        tree->Branch("DepA", &DepA, "DepA/D");
        tree->Branch("DepB", &DepB, "DepB/D");
        tree->Branch("DepC", &DepC, "DepC/D");
        tree->Branch("DepV", &DepV, "DepV/D");
        tree->Branch("DepW", &DepW, "DepW/D");
        //
        tree->Branch("mc_e_p", &mc_e_p, "mc_e_p/D");
        tree->Branch("mc_e_theta", &mc_e_theta, "mc_e_theta/D");
        tree->Branch("mc_e_phi", &mc_e_phi, "mc_e_phi/D");
        tree->Branch("mc_vz_e", &mc_vz_e, "mc_vz_e/D");
        tree->Branch("mc_p_p", &mc_p_p, "mc_p_p/D");
        tree->Branch("mc_p_theta", &mc_p_theta, "mc_p_theta/D");
        tree->Branch("mc_p_phi", &mc_p_phi, "mc_p_phi/D");
        tree->Branch("mc_vz_p", &mc_vz_p, "mc_vz_p/D");
        tree->Branch("mc_Q2", &mc_Q2, "mc_Q2/D");
        tree->Branch("mc_W", &mc_W, "mc_W/D");
        tree->Branch("mc_Mx", &mc_Mx, "mc_Mx/D");
        tree->Branch("mc_Mx2", &mc_Mx2, "mc_Mx2/D");
        tree->Branch("mc_x", &mc_x, "mc_x/D");
        tree->Branch("mc_y", &mc_y, "mc_y/D");
        tree->Branch("mc_t", &mc_t, "mc_t/D");
        tree->Branch("mc_tmin", &mc_tmin, "mc_tmin/D");
        tree->Branch("mc_z", &mc_z, "mc_z/D");
        tree->Branch("mc_xF", &mc_xF, "mc_xF/D");
        tree->Branch("mc_pT", &mc_pT, "mc_pT/D");
        tree->Branch("mc_zeta", &mc_zeta, "mc_zeta/D");
        tree->Branch("mc_eta", &mc_eta, "mc_eta/D");
        tree->Branch("mc_phi", &mc_phi, "mc_phi/D");
        tree->Branch("mc_DepA", &mc_DepA, "mc_DepA/D");
        tree->Branch("mc_DepB", &mc_DepB, "mc_DepB/D");
        tree->Branch("mc_DepC", &mc_DepC, "mc_DepC/D");
        tree->Branch("mc_DepV", &mc_DepV, "mc_DepV/D");
        tree->Branch("mc_DepW", &mc_DepW, "mc_DepW/D");
        //
        tree->Branch("matching_e_pid", &matching_e_pid, "matching_e_pid/I");
        tree->Branch("matching_p1_pid", &matching_p1_pid, "matching_p1_pid/I");
        tree->Branch("mc_p1_parent", &mc_p1_parent, "mc_p1_parent/I");
    }
    // Case for two hadrons (dihadrons)
    else if (hadron_count == 2 && is_mc == 0) {
        // Link TTree branches to variables for two hadrons
        tree->Branch("runnum", &runnum, "runnum/I");
        tree->Branch("evnum", &evnum, "evnum/I");
        tree->Branch("helicity", &helicity, "helicity/I");
        tree->Branch("beam_pol", &beam_pol, "beam_pol/D");
        tree->Branch("target_pol", &target_pol, "target_pol/D");
        tree->Branch("e_p", &e_p, "e_p/D");
        tree->Branch("e_theta", &e_theta, "e_theta/D");
        tree->Branch("e_phi", &e_phi, "e_phi/D");
        tree->Branch("vz_e", &vz_e, "vz_e/D");
        tree->Branch("p1_p", &p1_p, "p1_p/D");
        tree->Branch("p1_theta", &p1_theta, "p1_theta/D");
        tree->Branch("p1_phi", &p1_phi, "p1_phi/D");
        tree->Branch("vz_p1", &vz_p1, "vz_p1/D");
        tree->Branch("p2_p", &p2_p, "p2_p/D");
        tree->Branch("p2_theta", &p2_theta, "p2_theta/D");
        tree->Branch("p2_phi", &p2_phi, "p2_phi/D");
        tree->Branch("vz_p2", &vz_p2, "vz_p2/D");
        tree->Branch("Q2", &Q2, "Q2/D");
        tree->Branch("W", &W, "W/D");
        tree->Branch("Mx", &Mx, "Mx/D");
        tree->Branch("Mx1", &Mx1, "Mx1/D");
        tree->Branch("Mx2", &Mx2, "Mx2/D");
        tree->Branch("x", &x, "x/D");
        tree->Branch("y", &y, "y/D");
        tree->Branch("t", &t, "t/D");
        tree->Branch("t1", &t1, "t1/D");
        tree->Branch("t2", &t2, "t2/D");
        tree->Branch("tmin", &tmin, "tmin/D");
        tree->Branch("z", &z, "z/D");
        tree->Branch("z1", &z1, "z1/D");
        tree->Branch("z2", &z2, "z2/D");
        tree->Branch("Mh", &Mh, "Mh/D");
        tree->Branch("xF", &xF, "xF/D");
        tree->Branch("xF1", &xF1, "xF1/D");
        tree->Branch("xF2", &xF2, "xF2/D");
        tree->Branch("pT", &pT, "pT/D");
        tree->Branch("pT1", &pT1, "pT1/D");
        tree->Branch("pT2", &pT2, "pT2/D");
        tree->Branch("pTpT", &pTpT, "pTpT/D");
        tree->Branch("zeta", &zeta, "zeta/D");
        tree->Branch("zeta1", &zeta1, "zeta1/D");
        tree->Branch("zeta2", &zeta2, "zeta2/D");
        tree->Branch("eta", &eta, "eta/D");
        tree->Branch("eta1", &eta1, "eta1/D");
        tree->Branch("eta2", &eta2, "eta2/D");
        tree->Branch("Delta_eta", &Delta_eta, "Delta_eta/D");
        tree->Branch("eta1_gN", &eta1_gN, "eta1_gN/D");
        tree->Branch("eta2_gN", &eta2_gN, "eta2_gN/D");
        tree->Branch("phi1", &phi1, "phi1/D");
        tree->Branch("phi2", &phi2, "phi2/D");
        tree->Branch("Delta_phi", &Delta_phi, "Delta_phi/D");
        tree->Branch("phi", &phi, "phi/D");
        tree->Branch("phiR", &phiR, "phiR/D");
        tree->Branch("theta", &theta, "theta/D");
        tree->Branch("DepA", &DepA, "DepA/D");
        tree->Branch("DepB", &DepB, "DepB/D");
        tree->Branch("DepC", &DepC, "DepC/D");
        tree->Branch("DepV", &DepV, "DepV/D");
        tree->Branch("DepW", &DepW, "DepW/D");
    }
    // Case for two hadrons (dihadrons) and is mc
    else if (hadron_count == 2 && is_mc == 1) {
        // Link TTree branches to variables for two hadrons
        tree->Branch("e_p", &e_p, "e_p/D");
        tree->Branch("e_theta", &e_theta, "e_theta/D");
        tree->Branch("e_phi", &e_phi, "e_phi/D");
        tree->Branch("vz_e", &vz_e, "vz_e/D");
        tree->Branch("p1_p", &p1_p, "p1_p/D");
        tree->Branch("p1_theta", &p1_theta, "p1_theta/D");
        tree->Branch("p1_phi", &p1_phi, "p1_phi/D");
        tree->Branch("vz_p1", &vz_p1, "vz_p1/D");
        tree->Branch("p2_p", &p2_p, "p2_p/D");
        tree->Branch("p2_theta", &p2_theta, "p2_theta/D");
        tree->Branch("p2_phi", &p2_phi, "p2_phi/D");
        tree->Branch("vz_p2", &vz_p2, "vz_p2/D");
        tree->Branch("Q2", &Q2, "Q2/D");
        tree->Branch("W", &W, "W/D");
        tree->Branch("Mx", &Mx, "Mx/D");
        tree->Branch("Mx1", &Mx1, "Mx1/D");
        tree->Branch("Mx2", &Mx2, "Mx2/D");
        tree->Branch("x", &x, "x/D");
        tree->Branch("y", &y, "y/D");
        tree->Branch("t", &t, "t/D");
        tree->Branch("t1", &t1, "t1/D");
        tree->Branch("t2", &t2, "t2/D");
        tree->Branch("tmin", &tmin, "tmin/D");
        tree->Branch("z", &z, "z/D");
        tree->Branch("z1", &z1, "z1/D");
        tree->Branch("z2", &z2, "z2/D");
        tree->Branch("Mh", &Mh, "Mh/D");
        tree->Branch("xF", &xF, "xF/D");
        tree->Branch("xF1", &xF1, "xF1/D");
        tree->Branch("xF2", &xF2, "xF2/D");
        tree->Branch("pT", &pT, "pT/D");
        tree->Branch("pT1", &pT1, "pT1/D");
        tree->Branch("pT2", &pT2, "pT2/D");
        tree->Branch("pTpT", &pTpT, "pTpT/D");
        tree->Branch("zeta", &zeta, "zeta/D");
        tree->Branch("zeta1", &zeta1, "zeta1/D");
        tree->Branch("zeta2", &zeta2, "zeta2/D");
        tree->Branch("eta", &eta, "eta/D");
        tree->Branch("eta1", &eta1, "eta1/D");
        tree->Branch("eta2", &eta2, "eta2/D");
        tree->Branch("Delta_eta", &Delta_eta, "Delta_eta/D");
        tree->Branch("eta1_gN", &eta1_gN, "eta1_gN/D");
        tree->Branch("eta2_gN", &eta2_gN, "eta2_gN/D");
        tree->Branch("phi1", &phi1, "phi1/D");
        tree->Branch("phi2", &phi2, "phi2/D");
        tree->Branch("Delta_phi", &Delta_phi, "Delta_phi/D");
        tree->Branch("phi", &phi, "phi/D");
        tree->Branch("phiR", &phiR, "phiR/D");
        tree->Branch("theta", &theta, "theta/D");
        tree->Branch("DepA", &DepA, "DepA/D");
        tree->Branch("DepB", &DepB, "DepB/D");
        tree->Branch("DepC", &DepC, "DepC/D");
        tree->Branch("DepV", &DepV, "DepV/D");
        tree->Branch("DepW", &DepW, "DepW/D");
        //
        tree->Branch("mc_e_p", &mc_e_p, "mc_e_p/D");
        tree->Branch("mc_e_theta", &mc_e_theta, "mc_e_theta/D");
        tree->Branch("mc_e_phi", &mc_e_phi, "mc_e_phi/D");
        tree->Branch("mc_vz_e", &mc_vz_e, "mc_vz_e/D");
        tree->Branch("mc_p1_p", &mc_p1_p, "mc_p1_p/D");
        tree->Branch("mc_p1_theta", &mc_p1_theta, "mc_p1_theta/D");
        tree->Branch("mc_p1_phi", &mc_p1_phi, "mc_p1_phi/D");
        tree->Branch("mc_vz_p1", &mc_vz_p1, "mc_vz_p1/D");
        tree->Branch("mc_p2_p", &mc_p2_p, "mc_p2_p/D");
        tree->Branch("mc_p2_theta", &mc_p2_theta, "mc_p2_theta/D");
        tree->Branch("mc_p2_phi", &mc_p2_phi, "mc_p2_phi/D");
        tree->Branch("mc_vz_p2", &mc_vz_p2, "mc_vz_p2/D");
        tree->Branch("mc_Q2", &mc_Q2, "mc_Q2/D");
        tree->Branch("mc_W", &mc_W, "mc_W/D");
        tree->Branch("mc_Mx", &mc_Mx, "mc_Mx/D");
        tree->Branch("mc_Mx1", &mc_Mx1, "mc_Mx1/D");
        tree->Branch("mc_Mx2", &mc_Mx2, "mc_Mx2/D");
        tree->Branch("mc_x", &mc_x, "mc_x/D");
        tree->Branch("mc_y", &mc_y, "mc_y/D");
        tree->Branch("mc_t", &mc_t, "mc_t/D");
        tree->Branch("mc_t1", &mc_t1, "mc_t1/D");
        tree->Branch("mc_t2", &mc_t2, "mc_t2/D");
        tree->Branch("mc_tmin", &mc_tmin, "mc_tmin/D");
        tree->Branch("mc_z", &mc_z, "mc_z/D");
        tree->Branch("mc_z1", &mc_z1, "mc_z1/D");
        tree->Branch("mc_z2", &mc_z2, "mc_z2/D");
        tree->Branch("mc_Mh", &mc_Mh, "mc_Mh/D");
        tree->Branch("mc_xF", &mc_xF, "mc_xF/D");
        tree->Branch("mc_xF1", &mc_xF1, "mc_xF1/D");
        tree->Branch("mc_xF2", &mc_xF2, "mc_xF2/D");
        tree->Branch("mc_pT", &mc_pT, "mc_pT/D");
        tree->Branch("mc_pT1", &mc_pT1, "mc_pT1/D");
        tree->Branch("mc_pT2", &mc_pT2, "mc_pT2/D");
        tree->Branch("mc_pTpT", &mc_pTpT, "mc_pTpT/D");
        tree->Branch("mc_zeta", &mc_zeta, "mc_zeta/D");
        tree->Branch("mc_zeta1", &mc_zeta1, "mc_zeta1/D");
        tree->Branch("mc_zeta2", &mc_zeta2, "mc_zeta2/D");
        tree->Branch("mc_eta", &mc_eta, "mc_eta/D");
        tree->Branch("mc_eta1", &mc_eta1, "mc_eta1/D");
        tree->Branch("mc_eta2", &mc_eta2, "mc_eta2/D");
        tree->Branch("mc_Delta_eta", &mc_Delta_eta, "mc_Delta_eta/D");
        tree->Branch("mc_eta1_gN", &mc_eta1_gN, "mc_eta1_gN/D");
        tree->Branch("mc_eta2_gN", &mc_eta2_gN, "mc_eta2_gN/D");
        tree->Branch("mc_phi1", &mc_phi1, "mc_phi1/D");
        tree->Branch("mc_phi2", &mc_phi2, "mc_phi2/D");
        tree->Branch("mc_Delta_phi", &mc_Delta_phi, "mc_Delta_phi/D");
        tree->Branch("mc_phi", &mc_phi, "mc_phi/D");
        tree->Branch("mc_phiR", &mc_phiR, "mc_phiR/D");
        tree->Branch("mc_theta", &mc_theta, "mc_theta/D");
        tree->Branch("mc_DepA", &mc_DepA, "mc_DepA/D");
        tree->Branch("mc_DepB", &mc_DepB, "mc_DepB/D");
        tree->Branch("mc_DepC", &mc_DepC, "mc_DepC/D");
        tree->Branch("mc_DepV", &mc_DepV, "mc_DepV/D");
        tree->Branch("mc_DepW", &mc_DepW, "mc_DepW/D");
        //
        tree->Branch("matching_e_pid", &matching_e_pid, "matching_e_pid/I");
        tree->Branch("matching_p1_pid", &matching_p1_pid, "matching_p1_pid/I");
        tree->Branch("matching_p2_pid", &matching_p2_pid, "matching_p2_pid/I");
        tree->Branch("mc_p1_parent", &mc_p1_parent, "mc_p1_parent/I");
        tree->Branch("mc_p2_parent", &mc_p2_parent, "mc_p2_parent/I");
    }
    // Case for two hadrons (trihadrons)
    else if (hadron_count == 3 && is_mc == 0) {

        // Link TTree branches to variables for three hadrons
        tree->Branch("runnum", &runnum, "runnum/I");
        tree->Branch("evnum", &evnum, "evnum/I");
        tree->Branch("helicity", &helicity, "helicity/I");
        tree->Branch("beam_pol", &beam_pol, "beam_pol/D");
        tree->Branch("target_pol", &target_pol, "target_pol/D");
        tree->Branch("e_p", &e_p, "e_p/D");
        tree->Branch("e_theta", &e_theta, "e_theta/D");
        tree->Branch("e_phi", &e_phi, "e_phi/D");
        tree->Branch("vz_e", &vz_e, "vz_e/D");
        tree->Branch("p1_p", &p1_p, "p1_p/D");
        tree->Branch("p1_theta", &p1_theta, "p1_theta/D");
        tree->Branch("p1_phi", &p1_phi, "p1_phi/D");
        tree->Branch("vz_p1", &vz_p1, "vz_p1/D");
        tree->Branch("p2_p", &p2_p, "p2_p/D");
        tree->Branch("p2_theta", &p2_theta, "p2_theta/D");
        tree->Branch("p2_phi", &p2_phi, "p2_phi/D");
        tree->Branch("vz_p2", &vz_p2, "vz_p2/D");
        tree->Branch("p3_p", &p3_p, "p3_p/D");
        tree->Branch("p3_theta", &p3_theta, "p3_theta/D");
        tree->Branch("p3_phi", &p3_phi, "p3_phi/D");
        tree->Branch("vz_p3", &vz_p3, "vz_p3/D");
        tree->Branch("Q2", &Q2, "Q2/D");
        tree->Branch("W", &W, "W/D");
        tree->Branch("Mx", &Mx, "Mx/D");
        tree->Branch("Mx1", &Mx1, "Mx1/D");
        tree->Branch("Mx2", &Mx2, "Mx2/D");
        tree->Branch("Mx3", &Mx3, "Mx3/D");
        tree->Branch("Mx12", &Mx12, "Mx12/D");
        tree->Branch("Mx13", &Mx13, "Mx13/D");
        tree->Branch("Mx23", &Mx23, "Mx23/D");
        tree->Branch("x", &x, "x/D");
        tree->Branch("y", &y, "y/D");
        tree->Branch("t1", &t1, "t1/D");
        tree->Branch("t2", &t2, "t2/D");
        tree->Branch("t3", &t3, "t3/D");
        tree->Branch("tmin", &tmin, "tmin/D");
        tree->Branch("z", &z, "z/D");
        tree->Branch("z1", &z1, "z1/D");
        tree->Branch("z2", &z2, "z2/D");
        tree->Branch("z3", &z3, "z3/D");
        tree->Branch("z12", &z12, "z12/D");
        tree->Branch("z13", &z13, "z13/D");
        tree->Branch("z23", &z23, "z23/D");
        tree->Branch("zeta", &zeta, "zeta/D");
        tree->Branch("zeta1", &zeta1, "zeta1/D");
        tree->Branch("zeta2", &zeta2, "zeta2/D");
        tree->Branch("zeta3", &zeta3, "zeta3/D");
        tree->Branch("zeta12", &zeta12, "zeta12/D");
        tree->Branch("zeta13", &zeta13, "zeta13/D");
        tree->Branch("zeta23", &zeta23, "zeta23/D");
        tree->Branch("pT", &pT, "pT/D");
        tree->Branch("pT1", &pT1, "pT1/D");
        tree->Branch("pT2", &pT2, "pT2/D");
        tree->Branch("pT3", &pT3, "pT3/D");
        tree->Branch("pT12", &pT12, "pT12/D");
        tree->Branch("pT13", &pT13, "pT13/D");
        tree->Branch("pT23", &pT23, "pT23/D");
        tree->Branch("Mh", &Mh, "Mh/D");
        tree->Branch("Mh12", &Mh12, "Mh12/D");
        tree->Branch("Mh13", &Mh13, "Mh13/D");
        tree->Branch("Mh23", &Mh23, "Mh23/D");
        tree->Branch("xF", &xF, "xF/D");
        tree->Branch("xF1", &xF1, "xF1/D");
        tree->Branch("xF2", &xF2, "xF2/D");
        tree->Branch("xF3", &xF3, "xF3/D");
        tree->Branch("xF12", &xF12, "xF12/D");
        tree->Branch("xF13", &xF13, "xF13/D");
        tree->Branch("xF23", &xF23, "xF23/D");
        tree->Branch("eta", &eta, "eta/D");
        tree->Branch("eta1", &eta1, "eta1/D");
        tree->Branch("eta2", &eta2, "eta2/D");
        tree->Branch("eta3", &eta3, "eta3/D");
        tree->Branch("eta12", &eta12, "eta12/D");
        tree->Branch("eta13", &eta13, "eta13/D");
        tree->Branch("eta23", &eta23, "eta23/D");
        tree->Branch("phi1", &phi1, "phi1/D");
        tree->Branch("phi2", &phi2, "phi2/D");
        tree->Branch("phi3", &phi3, "phi3/D");
        tree->Branch("phi12", &phi12, "phi12/D");
        tree->Branch("phi13", &phi13, "phi13/D");
        tree->Branch("phi23", &phi23, "phi23/D");
        tree->Branch("phih", &phih, "phih/D");
        tree->Branch("phiR", &phiR, "phiR/D");
        tree->Branch("theta", &theta, "theta/D");
        tree->Branch("Delta_phi12", &Delta_phi12, "Delta_phi12/D");
        tree->Branch("Delta_phi13", &Delta_phi13, "Delta_phi13/D");
        tree->Branch("Delta_phi23", &Delta_phi23, "Delta_phi23/D");
        tree->Branch("DepA", &DepA, "DepA/D");
        tree->Branch("DepB", &DepB, "DepB/D");
        tree->Branch("DepC", &DepC, "DepC/D");
        tree->Branch("DepV", &DepV, "DepV/D");
        tree->Branch("DepW", &DepW, "DepW/D");
    }
    // Case for dvcs (not actually two hadrons)
    else if (hadron_count == 4 && is_mc == 0) { // dvcs, not actually 4 hadrons
        // Link TTree branches to variables for dvcs
        tree->Branch("runnum", &runnum, "runnum/I");
        tree->Branch("evnum", &evnum, "evnum/I");
        tree->Branch("helicity", &helicity, "helicity/I");
        tree->Branch("beam_pol", &beam_pol, "beam_pol/D");
        tree->Branch("target_pol", &target_pol, "target_pol/D");
        tree->Branch("e_p", &e_p, "e_p/D");
        tree->Branch("e_theta", &e_theta, "e_theta/D");
        tree->Branch("e_phi", &e_phi, "e_phi/D");
        tree->Branch("vz_e", &vz_e, "vz_e/D");
        tree->Branch("p1_p", &p1_p, "p1_p/D");
        tree->Branch("p1_theta", &p1_theta, "p1_theta/D");
        tree->Branch("p1_phi", &p1_phi, "p1_phi/D");
        tree->Branch("vz_p1", &vz_p1, "vz_p1/D");
        tree->Branch("p2_p", &p2_p, "p2_p/D");
        tree->Branch("p2_theta", &p2_theta, "p2_theta/D");
        tree->Branch("p2_phi", &p2_phi, "p2_phi/D");
        tree->Branch("vz_p2", &vz_p2, "vz_p2/D");
        tree->Branch("Q2", &Q2, "Q2/D");
        tree->Branch("W", &W, "W/D");
        tree->Branch("Mx", &Mx, "Mx/D");
        tree->Branch("Mx1", &Mx1, "Mx1/D");
        tree->Branch("Mx2", &Mx2, "Mx2/D");
        tree->Branch("x", &x, "x/D");
        tree->Branch("y", &y, "y/D");
        tree->Branch("t", &t, "t/D");
        tree->Branch("t1", &t1, "t1/D");
        tree->Branch("t2", &t2, "t2/D");
        tree->Branch("tmin", &tmin, "tmin/D");
        tree->Branch("z", &z, "z/D");
        tree->Branch("z1", &z1, "z1/D");
        tree->Branch("z2", &z2, "z2/D");
        tree->Branch("Mh", &Mh, "Mh/D");
        tree->Branch("xF", &xF, "xF/D");
        tree->Branch("xF1", &xF1, "xF1/D");
        tree->Branch("xF2", &xF2, "xF2/D");
        tree->Branch("pT", &pT, "pT/D");
        tree->Branch("pT1", &pT1, "pT1/D");
        tree->Branch("pT2", &pT2, "pT2/D");
        tree->Branch("pTpT", &pTpT, "pTpT/D");
        tree->Branch("zeta", &zeta, "zeta/D");
        tree->Branch("zeta1", &zeta1, "zeta1/D");
        tree->Branch("zeta2", &zeta2, "zeta2/D");
        tree->Branch("eta", &eta, "eta/D");
        tree->Branch("eta1", &eta1, "eta1/D");
        tree->Branch("eta2", &eta2, "eta2/D");
        tree->Branch("Delta_eta", &Delta_eta, "Delta_eta/D");
        tree->Branch("eta1_gN", &eta1_gN, "eta1_gN/D");
        tree->Branch("eta2_gN", &eta2_gN, "eta2_gN/D");
        tree->Branch("phi1", &phi1, "phi1/D");
        tree->Branch("phi2", &phi2, "phi2/D");
        tree->Branch("Delta_phi", &Delta_phi, "Delta_phi/D");
        tree->Branch("phi", &phi, "phi/D");
        tree->Branch("phiR", &phiR, "phiR/D");
        tree->Branch("theta", &theta, "theta/D");
        tree->Branch("DepA", &DepA, "DepA/D");
        tree->Branch("DepB", &DepB, "DepB/D");
        tree->Branch("DepC", &DepC, "DepC/D");
        tree->Branch("DepV", &DepV, "DepV/D");
        tree->Branch("DepW", &DepW, "DepW/D");
        tree->Branch("Emiss2", &Emiss2, "Emiss2/D");
        tree->Branch("theta_gamma_gamma", &theta_gamma_gamma, "theta_gamma_gamma/D");
        tree->Branch("pTmiss", &pTmiss, "pTmiss/D");
        tree->Branch("Mxgammasquared", &Mxgammasquared, "Mxgammasquared/D");
    }
    // Case for calibration script (not actually five hadrons)
    if (hadron_count == 5 && is_mc == 0) {
        // Link TTree branches to variables for four hadrons
        tree->Branch("config_run", &config_run, "config_run/I");
        tree->Branch("config_event", &config_event, "config_event/I");
        tree->Branch("event_helicity", &event_helicity, "event_helicity/I");
        tree->Branch("particle_pid", &particle_pid, "particle_pid/I");
        tree->Branch("particle_status", &particle_status, "particle_status/I");
        tree->Branch("particle_px", &particle_px, "particle_px/D");
        tree->Branch("particle_py", &particle_py, "particle_py/D");
        tree->Branch("particle_pz", &particle_pz, "particle_pz/D");
        tree->Branch("p", &p, "p/D");
        tree->Branch("theta", &theta, "theta/D");
        tree->Branch("phi", &phi, "phi/D");
        tree->Branch("particle_vx", &particle_vx, "particle_vx/D");
        tree->Branch("particle_vy", &particle_vy, "particle_vy/D");
        tree->Branch("particle_vz", &particle_vz, "particle_vz/D");
        tree->Branch("particle_beta", &particle_beta, "particle_beta/D");
        tree->Branch("particle_chi2pid", &particle_chi2pid, "particle_chi2pid/D");

        tree->Branch("cal_sector", &cal_sector, "cal_sector/I");
        tree->Branch("cal_energy_1", &cal_energy_1, "cal_energy_1/D");
        tree->Branch("cal_x_1", &cal_x_1, "cal_x_1/D");
        tree->Branch("cal_y_1", &cal_y_1, "cal_y_1/D");
        tree->Branch("cal_z_1", &cal_z_1, "cal_z_1/D");
        tree->Branch("cal_lu_1", &cal_lu_1, "cal_lu_1/D");
        tree->Branch("cal_lv_1", &cal_lv_1, "cal_lv_1/D");
        tree->Branch("cal_lw_1", &cal_lw_1, "cal_lw_1/D");
        tree->Branch("cal_energy_4", &cal_energy_4, "cal_energy_4/D");
        tree->Branch("cal_x_4", &cal_x_4, "cal_x_4/D");
        tree->Branch("cal_y_4", &cal_y_4, "cal_y_4/D");
        tree->Branch("cal_z_4", &cal_z_4, "cal_z_4/D");
        tree->Branch("cal_lu_4", &cal_lu_4, "cal_lu_4/D");
        tree->Branch("cal_lv_4", &cal_lv_4, "cal_lv_4/D");
        tree->Branch("cal_lw_4", &cal_lw_4, "cal_lw_4/D");
        tree->Branch("cal_energy_7", &cal_energy_7, "cal_energy_7/D");
        tree->Branch("cal_x_7", &cal_x_7, "cal_x_7/D");
        tree->Branch("cal_y_7", &cal_y_7, "cal_y_7/D");
        tree->Branch("cal_z_7", &cal_z_7, "cal_z_7/D");
        tree->Branch("cal_lu_7", &cal_lu_7, "cal_lu_7/D");
        tree->Branch("cal_lv_7", &cal_lv_7, "cal_lv_7/D");
        tree->Branch("cal_lw_7", &cal_lw_7, "cal_lw_7/D");

        tree->Branch("cc_nphe_15", &cc_nphe_15, "cc_nphe_15/D");
        tree->Branch("cc_nphe_16", &cc_nphe_16, "cc_nphe_16/D");

        tree->Branch("track_sector_5", &track_sector_5, "track_sector_5/I");
        tree->Branch("track_chi2_5", &track_chi2_5, "track_chi2_5/D");
        tree->Branch("track_ndf_5", &track_ndf_5, "track_ndf_5/I");
        tree->Branch("track_sector_6", &track_sector_6, "track_sector_6/I");
        tree->Branch("track_chi2_6", &track_chi2_6, "track_chi2_6/D");
        tree->Branch("track_ndf_6", &track_ndf_6, "track_ndf_6/I");

        tree->Branch("traj_x_6", &traj_x_6, "traj_x_6/D");
        tree->Branch("traj_y_6", &traj_y_6, "traj_y_6/D");
        tree->Branch("traj_z_6", &traj_z_6, "traj_z_6/D");
        tree->Branch("traj_edge_6", &traj_edge_6, "traj_edge_6/D");
        tree->Branch("traj_x_18", &traj_x_18, "traj_x_18/D");
        tree->Branch("traj_y_18", &traj_y_18, "traj_y_18/D");
        tree->Branch("traj_z_18", &traj_z_18, "traj_z_18/D");
        tree->Branch("traj_edge_18", &traj_edge_18, "traj_edge_18/D");
        tree->Branch("traj_x_36", &traj_x_36, "traj_x_36/D");
        tree->Branch("traj_y_36", &traj_y_36, "traj_y_36/D");
        tree->Branch("traj_z_36", &traj_z_36, "traj_z_36/D");
        tree->Branch("traj_edge_36", &traj_edge_36, "traj_edge_36/D");

        tree->Branch("traj_x_1", &traj_x_1, "traj_x_1/D");
        tree->Branch("traj_y_1", &traj_y_1, "traj_y_1/D");
        tree->Branch("traj_z_1", &traj_z_1, "traj_z_1/D");
        tree->Branch("traj_edge_1", &traj_edge_1, "traj_edge_1/D");
        tree->Branch("traj_x_3", &traj_x_3, "traj_x_3/D");
        tree->Branch("traj_y_3", &traj_y_3, "traj_y_3/D");
        tree->Branch("traj_z_3", &traj_z_3, "traj_z_3/D");
        tree->Branch("traj_edge_3", &traj_edge_3, "traj_edge_3/D");
        tree->Branch("traj_x_5", &traj_x_5, "traj_x_5/D");
        tree->Branch("traj_y_5", &traj_y_5, "traj_y_5/D");
        tree->Branch("traj_z_5", &traj_z_5, "traj_z_5/D");
        tree->Branch("traj_edge_5", &traj_edge_5, "traj_edge_5/D");
        tree->Branch("traj_x_7", &traj_x_7, "traj_x_7/D");
        tree->Branch("traj_y_7", &traj_y_7, "traj_y_7/D");
        tree->Branch("traj_z_7", &traj_z_7, "traj_z_7/D");
        tree->Branch("traj_edge_7", &traj_edge_7, "traj_edge_7/D");
        tree->Branch("traj_x_12", &traj_x_12, "traj_x_12/D");
        tree->Branch("traj_y_12", &traj_y_12, "traj_y_12/D");
        tree->Branch("traj_z_12", &traj_z_12, "traj_z_12/D");
        tree->Branch("traj_edge_12", &traj_edge_12, "traj_edge_12/D");

        tree->Branch("ft_energy", &ft_energy, "ft_energy/D");
        tree->Branch("ft_x", &ft_x, "ft_x/D");
        tree->Branch("ft_y", &ft_y, "ft_y/D");
        tree->Branch("ft_z", &ft_z, "ft_z/D");
        tree->Branch("ft_radius", &ft_radius, "ft_radius/D");
    }

    // Find the root directory of the repository
    std::string package_location = findPackageRoot();
    // Define the CSV path relative to the package root
    std::string csv_location = "analysis_scripts/asymmetry_extraction/imports/clas12_run_info.csv";
    // Load run info from the CSV file
    load_run_info_from_csv(package_location + csv_location);
    

    // Loop to read each line from the text file and fill the TTree based on hadron_count
    if (hadron_count == 0 && is_mc == 0) {
        while (infile >> runnum >> evnum >> helicity >> e_p >> e_theta >> e_phi >> vz_e >> 
            Q2 >> W >> Mx >> Mx2 >> x >> y >> DepA >> DepB >> DepC >> DepV >> DepW) {

            beam_pol = getPol(runnum);
            if (runnum < 16000) { target_pol = 0; }
            else { 
                for (const auto& run_info : run_info_list) {
                    if (run_info.runnum == runnum) {
                        target_pol = run_info.target_polarization;
                        break;
                    }
                }
            }

            t = gett(e_p, e_theta); // for inclusive we calculate t with electron kinematics
            tmin = gettmin(x);  

            tree->Fill(); // Fill the tree with the read data
        }
    } 
    if (hadron_count == 0 && is_mc == 1) {
        while (infile >> e_p >> mc_e_p >> e_theta >> mc_e_theta >>
            e_phi >> mc_e_phi >> vz_e >> mc_vz_e >> Q2 >> mc_Q2 >> W >> mc_W >> Mx >> mc_Mx >>
            Mx2 >> mc_Mx2 >> x >> mc_x >> y >> mc_y >> DepA >> mc_DepA >> DepB >> mc_DepB >> 
            DepC >> mc_DepC >> DepV >> mc_DepV >> DepW >> mc_DepW >> matching_e_pid) {

            beam_pol = getPol(runnum);
            if (runnum < 16000) { target_pol = 0; }
            else { 
                for (const auto& run_info : run_info_list) {
                    if (run_info.runnum == runnum) {
                        target_pol = run_info.target_polarization;
                        break;
                    }
                }
            }

            t = gett(e_p, e_theta); // for inclusive we calculate t with electron kinematics
            tmin = gettmin(x);  
            mc_t = gett(mc_e_p, mc_e_theta);
            mc_tmin = gettmin(mc_x);  

            tree->Fill(); // Fill the tree with the read data
        }
    } 
    if (hadron_count == 1 && is_mc == 0) {
        // double helicity_double;
        while (infile >> runnum >> evnum >> helicity >> e_p >> e_theta >> e_phi >> vz_e >> 
            p_p >> p_theta >> p_phi >> vz_p >> Q2 >> W >> Mx >> Mx2 >> x >> y >> z >> xF >> 
            pT >> zeta >> eta >> phi >> DepA >> DepB >> DepC >> DepV >> DepW) {

            beam_pol = getPol(runnum);
            if (runnum < 16000) { target_pol = 0; }
            else { 
                for (const auto& run_info : run_info_list) {
                    if (run_info.runnum == runnum) {
                        target_pol = run_info.target_polarization;
                        break;
                    }
                }
            }

            t = gett(p_p, p_theta); // for SIDIS we calculate t with proton kinematics
            tmin = gettmin(x); 


            tree->Fill(); // Fill the tree with the read data
        }
    }
    if (hadron_count == 1 && is_mc == 1) {
        while (infile >> e_p >> mc_e_p >> e_theta >> mc_e_theta >> e_phi >> mc_e_phi >> vz_e >> 
            mc_vz_e >> p_p >> mc_p_p >> p_theta >> mc_p_theta >> p_phi >> mc_p_phi >> vz_p >>
            mc_vz_p >> Q2 >> mc_Q2 >> W >> mc_W >> Mx >> mc_Mx >> Mx2 >> mc_Mx2 >> x >> mc_x >> 
            y >> mc_y >> z >> mc_z >> xF >> mc_xF >> pT >> mc_pT >> zeta >> mc_zeta >> eta >> 
            mc_eta >> phi >> mc_phi >> DepA >> mc_DepA >> DepB >> mc_DepB >> DepC >> mc_DepC >> 
            DepV >> mc_DepV >> DepW >> mc_DepW >> matching_e_pid >> matching_p1_pid >>
            mc_p1_parent) {

            runnum = 11;
            beam_pol = 0;
            target_pol = 0;

            t = gett(p_p, p_theta); // for SIDIS we calculate t with proton kinematics
            mc_t = gett(mc_p_p, mc_p_theta);
            tmin = gettmin(x);
            mc_tmin = gettmin(mc_x); 

            tree->Fill(); // Fill the tree with the read data
        }
    } 
    if (hadron_count == 2 && is_mc == 0) {
        while (infile >> runnum >> evnum >> helicity >> e_p >> e_theta >> e_phi >> vz_e >> 
            p1_p >> p1_theta >> p1_phi >> vz_p1 >> p2_p >> p2_theta >> p2_phi >> vz_p2 >> 
            Q2 >> W >> Mx >> Mx1 >> Mx2 >> x >> y >> z >> z1 >> z2 >> Mh >> xF >> xF1 >> xF2 >> 
            pT >> pT1 >> pT2 >> pTpT >> zeta >> zeta1 >> zeta2 >> eta >> eta1 >> eta2 >> Delta_eta>> 
            eta1_gN >> eta2_gN >> phi1 >> phi2 >> Delta_phi >> phi >> phiR >> theta >> 
            DepA >> DepB >> DepC >> DepV >> DepW) {

            beam_pol = getPol(runnum);
            if (runnum < 16000) { target_pol = 0; }
            else { 
                for (const auto& run_info : run_info_list) {
                    if (run_info.runnum == runnum) {
                        target_pol = run_info.target_polarization;
                        break;
                    }
                }
            }

            // Convert spherical coordinates to Cartesian coordinates for each hadron
            double p1_x, p1_y, p1_z;
            sphericalToCartesian(p1_p, p1_phi, p1_theta, p1_x, p1_y, p1_z);
            
            double p2_x, p2_y, p2_z;
            sphericalToCartesian(p2_p, p2_phi, p2_theta, p2_x, p2_y, p2_z);
            
            // Calculate the total momentum components of the parent hadron
            double p_parent_x = p1_x + p2_x;
            double p_parent_y = p1_y + p2_y;
            double p_parent_z = p1_z + p2_z;
            
            // Calculate the magnitude of the parent hadron's momentum
            p_p = vectorMagnitude(p_parent_x, p_parent_y, p_parent_z);
            
            // Calculate the polar angle of the parent hadron's momentum
            p_theta = vectorPolarAngle(p_parent_z, p_p);

            t = gett(p_p, p_theta);
            t1 = gett(p1_p, p1_theta);
            t2 = gett(p2_p, p2_theta); 
            tmin = gettmin(x); 

            tree->Fill(); // Fill the tree with the read data
        }
    }
    
    if (hadron_count == 2 && is_mc == 1) {

        while (infile >> e_p >> mc_e_p >> e_theta >> mc_e_theta >> e_phi >> mc_e_phi >> vz_e >> 
            mc_vz_e >> p1_p >> mc_p1_p >> p1_theta >> mc_p1_theta >> p1_phi >> mc_p1_phi >> 
            vz_p1 >> mc_vz_p1 >> p2_p >> mc_p2_p >> p2_theta >> mc_p2_theta >> p2_phi >> 
            mc_p2_phi >> vz_p2 >> mc_vz_p2 >> Q2 >> mc_Q2 >> W >> mc_W >> Mx >> mc_Mx >> Mx1 >> 
            mc_Mx1 >> Mx2 >> mc_Mx2 >> x >> mc_x >> y >> mc_y >> z >> mc_z >> z1 >> mc_z1 >> 
            z2 >> mc_z2 >> Mh >> mc_Mh >> xF >> mc_xF >> xF1 >> mc_xF1 >> xF2 >> mc_xF2 >> 
            pT >> mc_pT >> pT1 >> mc_pT1 >> pT2 >> mc_pT2 >> pTpT >> mc_pTpT >> zeta >> mc_zeta >>
            zeta1 >> mc_zeta1 >> zeta2 >> mc_zeta2 >> eta >> mc_eta >> eta1 >> mc_eta1 >> eta2 >> 
            mc_eta2 >> Delta_eta >> mc_Delta_eta >> eta1_gN >> mc_eta1_gN >> eta2_gN >> 
            mc_eta2_gN >> phi1 >> mc_phi1 >> phi2 >> mc_phi2 >> Delta_phi >> mc_Delta_phi >> 
            phi >> mc_phi >> phiR >> mc_phiR >> theta >> mc_theta >> DepA >> mc_DepA >> DepB >> 
            mc_DepB >> DepC >> mc_DepC >> DepV >> mc_DepV >> DepW >> mc_DepW >> matching_e_pid >>
            matching_p1_pid >> matching_p2_pid >> mc_p1_parent >> mc_p2_parent) {

            // Convert spherical coordinates to Cartesian coordinates for each hadron
            double p1_x, p1_y, p1_z;
            sphericalToCartesian(p1_p, p1_phi, p1_theta, p1_x, p1_y, p1_z);
            double mc_p1_x, mc_p1_y, mc_p1_z;
            sphericalToCartesian(mc_p1_p, mc_p1_phi, mc_p1_theta, mc_p1_x, mc_p1_y, mc_p1_z);
            
            double p2_x, p2_y, p2_z;
            sphericalToCartesian(p2_p, p2_phi, p2_theta, p2_x, p2_y, p2_z);
            double mc_p2_x, mc_p2_y, mc_p2_z;
            sphericalToCartesian(mc_p2_p, mc_p2_phi, mc_p2_theta, mc_p2_x, mc_p2_y, mc_p2_z);
            
            // Calculate the total momentum components of the parent hadron
            double p_parent_x = p1_x + p2_x;
            double p_parent_y = p1_y + p2_y;
            double p_parent_z = p1_z + p2_z;
            double mc_p_parent_x = mc_p1_x + mc_p2_x;
            double mc_p_parent_y = mc_p1_y + mc_p2_y;
            double mc_p_parent_z = mc_p1_z + mc_p2_z;
            
            // Calculate the magnitude of the parent hadron's momentum
            p_p = vectorMagnitude(p_parent_x, p_parent_y, p_parent_z);
            mc_p_p = vectorMagnitude(mc_p_parent_x, mc_p_parent_y, mc_p_parent_z);
            
            // Calculate the polar angle of the parent hadron's momentum
            p_theta = vectorPolarAngle(p_parent_z, p_p);
            mc_p_theta = vectorPolarAngle(mc_p_parent_z, mc_p_p);

            t = gett(p_p, p_theta);
            t1 = gett(p1_p, p1_theta);
            t2 = gett(p2_p, p2_theta); 
            tmin = gettmin(x); 
            mc_t = gett(mc_p_p, mc_p_theta);
            mc_t1 = gett(mc_p1_p, mc_p1_theta);
            mc_t2 = gett(mc_p2_p, mc_p2_theta); 
            mc_tmin = gettmin(mc_x); 

            tree->Fill(); // Fill the tree with the read data
        }
    }
    if (hadron_count == 3 && is_mc == 0) {
        while (infile >> runnum >> evnum >> helicity >> e_p >> e_theta >> e_phi >> vz_e >> 
            p1_p >> p1_theta >> p1_phi >> vz_p1 >> p2_p >> p2_theta >> p2_phi >> vz_p2 >>
            p3_p >> p3_theta >> p3_phi >> vz_p3 >> 
            Q2 >> W >> Mx >> Mx1 >> Mx2 >> Mx3 >> Mx12 >> Mx13 >> Mx23 >> 
            x >> y >> z >> z1 >> z2 >> z3 >> z12 >> z13 >> z23 >> 
            zeta >> zeta1 >> zeta2 >> zeta3 >> zeta12 >> zeta13 >> zeta23 >>
            pT >> pT1 >> pT2 >> pT3 >> pT12 >> pT13 >> pT23 >> 
            Mh >> Mh12 >> Mh13 >> Mh23 >> 
            xF >> xF1 >> xF2 >> xF3 >> xF12 >> xF13 >> xF23 >> 
            eta >> eta1 >> eta2 >> eta3 >> eta12 >> eta13 >> eta23 >> 
            phi1 >> phi2 >> phi3 >> phi12 >> phi13 >> phi23 >> phih >> phiR >> theta >>
            Delta_phi12 >> Delta_phi13 >> Delta_phi23 >> 
            DepA >> DepB >> DepC >> DepV >> DepW) {

            beam_pol = getPol(runnum);
            if (runnum < 16000) { target_pol = 0; }
            else { 
                for (const auto& run_info : run_info_list) {
                    if (run_info.runnum == runnum) {
                        target_pol = run_info.target_polarization;
                        break;
                    }
                }
            }

            t1 = gett(p1_p, p1_theta);
            t2 = gett(p2_p, p2_theta); 
            t3 = gett(p3_p, p3_theta); 
            tmin = gettmin(x); 

            tree->Fill(); // Fill the tree with the read data
        }
    }

    if (hadron_count == 4 && is_mc == 0) {
        while (infile >> runnum >> evnum >> helicity >> e_p >> e_theta >> e_phi >> vz_e >> 
            p1_p >> p1_theta >> p1_phi >> vz_p1 >> p2_p >> p2_theta >> p2_phi >> vz_p2 >> 
            Q2 >> W >> Mx >> Mx1 >> Mx2 >> x >> y >> z >> z1 >> z2 >> Mh >> xF >> xF1 >> xF2 >> 
            pT >> pT1 >> pT2 >> pTpT >> zeta >> zeta1 >> zeta2 >> eta >> eta1 >> eta2 >> Delta_eta>> 
            eta1_gN >> eta2_gN >> phi1 >> phi2 >> Delta_phi >> phi >> phiR >> theta >> 
            DepA >> DepB >> DepC >> DepV >> DepW >> Emiss2 >> theta_gamma_gamma >> pTmiss >>
            Mxgammasquared) {

            beam_pol = getPol(runnum);
            if (runnum < 16000) { target_pol = 0; }
            else { 
                for (const auto& run_info : run_info_list) {
                    if (run_info.runnum == runnum) {
                        target_pol = run_info.target_polarization;
                        break;
                    }
                }
            }

            // Convert spherical coordinates to Cartesian coordinates for each hadron
            double p1_x, p1_y, p1_z;
            sphericalToCartesian(p1_p, p1_phi, p1_theta, p1_x, p1_y, p1_z);
            
            double p2_x, p2_y, p2_z;
            sphericalToCartesian(p2_p, p2_phi, p2_theta, p2_x, p2_y, p2_z);
            
            // Calculate the total momentum components of the parent hadron
            double p_parent_x = p1_x + p2_x;
            double p_parent_y = p1_y + p2_y;
            double p_parent_z = p1_z + p2_z;
            
            // Calculate the magnitude of the parent hadron's momentum
            p_p = vectorMagnitude(p_parent_x, p_parent_y, p_parent_z);
            
            // Calculate the polar angle of the parent hadron's momentum
            p_theta = vectorPolarAngle(p_parent_z, p_p);

            t = gett(p_p, p_theta);
            t1 = gett(p1_p, p1_theta);
            t2 = gett(p2_p, p2_theta); 
            tmin = gettmin(x); 

            tree->Fill(); // Fill the tree with the read data
        }
    }

    // calibration script
    if (hadron_count == 5 && is_mc == 0) {
        while (infile >> config_run >> config_event >> event_helicity >> particle_pid >>
                particle_px >> particle_py >> particle_pz >> p >> theta >> phi >> 
                particle_vx >> particle_vy >> particle_vz >> particle_beta >> 
                particle_chi2pid >> particle_status >> cal_sector >> cal_energy_1 >> 
                cal_x_1 >> cal_y_1 >> cal_z_1 >> cal_lu_1 >> cal_lv_1 >> cal_lw_1 >> 
                cal_energy_4 >> cal_x_4 >> cal_y_4 >> cal_z_4 >> cal_lu_4 >> cal_lv_4 >> cal_lw_4 >> 
                cal_energy_7 >> cal_x_7 >> cal_y_7 >> cal_z_7 >> cal_lu_7 >> cal_lv_7 >> cal_lw_7 >> 
                cc_nphe_15 >> cc_nphe_16 >> track_sector_5 >> track_chi2_5 >> track_ndf_5 >> 
                track_sector_6 >> track_chi2_6 >> track_ndf_6 >> traj_x_6 >> traj_y_6 >> traj_z_6 >> traj_edge_6 >> 
                traj_x_18 >> traj_y_18 >> traj_z_18 >> traj_edge_18 >> traj_x_36 >> traj_y_36 >> 
                traj_z_36 >> traj_edge_36 >> traj_x_1 >> traj_y_1 >> traj_z_1 >> traj_edge_1 >> 
                traj_x_3 >> traj_y_3 >> traj_z_3 >> traj_edge_3 >> traj_x_5 >> traj_y_5 >> 
                traj_z_5 >> traj_edge_5 >> traj_x_7 >> traj_y_7 >> traj_z_7 >> traj_edge_7 >> 
                traj_x_12 >> traj_y_12 >> traj_z_12 >> traj_edge_12 >> ft_energy >> ft_x >> 
                ft_y >> ft_z >> ft_radius) {

            tree->Fill(); // Fill the tree with the read data
        }
    }

    // Write the TTree to the ROOT file and close it
    tree->Write();
    outfile->Close();

    // Close the input text file
    infile.close();

    cout << "Output ROOT file: " << argv[2] << endl << endl;

    return 0;
}
