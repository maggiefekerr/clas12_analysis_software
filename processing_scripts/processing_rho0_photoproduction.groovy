/*
 * author Maggie F. E. Kerr
 * using processing_rho0_photoproduction.groovy 
 * and derivative files of
 * Timothy B. Hayward as template
 * 
 * rho0 photoproduction
 */

// import CLAS12 physics classes
import org.jlab.io.hipo.*;
import org.jlab.io.base.DataEvent;
import org.jlab.clas.physics.*;
import org.jlab.clas12.physics.*;

// import from hayward_coatjava_extensions
import extended_kinematic_fitters.*; 
import analyzers.*;

// filetype for gathering files in directory
import groovy.io.FileType;

// dilks CLAS QA analysis
import clasqa.QADB

public static void main(String[] args) {
    // Start time
	long startTime = System.currentTimeMillis();

    // Check if an argument is provided
	if (!args) {
	    // Print an error message and exit the program if the input directory is not specified
	    println("ERROR: Please enter a hipo file directory as the first argument");
	    System.exit(0);
	}

    // If the input directory is provided, iterate through each file recursively
	def hipo_list = []
	(args[0] as File).eachFileRecurse(FileType.FILES) 
		{ if (it.name.endsWith('.hipo')) hipo_list << it }
    
    String nucl_str = args.length >= 2 ? ((args[1].equals("2212") || args[1].equals("2112")) ? args[1] : "2212") : "2212";
    if (args.length < 2) println("WARNING: Specify either proton or neutron PDG PID for TCS type! Set to proton (2212).")
    if ((args[1] != "2212") && (args[1] !="2112")) println("WARNING: Specify either proton or neutron PDG PID for TCS type! Set to proton (2212).")
    println("Set target PID for rho0 production type = $nucl_str")
    int nucl_int = nucl_str.toInteger()

    // Set the PDG PID for rho0 decay product based on the provided 3rd argument or default to 211 (pi+)
	String p1_Str = args.length < 3 ? "211" : args[2];
	if (args.length < 3) println("WARNING: Specify a PDG PID for p1! Set to pi+ (211).");
	println("Set p1 PID = $p1_Str");
	int p1_int = p1_Str.toInteger(); // Convert p1_Str to integer

    String output_file = args.length < 4 ? "rho0_dummy_out.txt" : args[3]
    if (args.length < 4) println('WARNING: Specify an output file name. Set to "rho0_dummy_out.txt".')
    File file = new File(output_file)
    file.delete()
    BufferedWriter writer = new BufferedWriter(new FileWriter(file))

    int n_files = args.length < 5 || Integer.parseInt(args[4]) == 0 || Integer.parseInt(args[4]) > hipo_list.size()
        ? hipo_list.size() : Integer.parseInt(args[4])
    if (args.length < 5 || Integer.parseInt(args[4]) == 0 || Integer.parseInt(args[4]) > hipo_list.size()) {
        println("WARNING: Number of files not specified, set to 0, or number too large.")
        println("Setting # of files to be equal to number of files in the directory.")
        println("There are $hipo_list.size files.")
    }

    double beam_energy = args.length < 6 ? 10.6 : Double.parseDouble(args[5])
    if (args.length < 6) {
        println("No beam energy provided, defaulting to 10.6 GeV.")
        println("All MC will use 10.604 GeV. You must manually enter a beam energy to change this.")
    }

    Integer userProvidedRun = null
    if (args.length < 7 || Integer.parseInt(args[6]) == 0) {
        println("Run number not provided, will pull from hipo files.")
        println("Think carefully about this if you are processing MC.")
    } else {
        userProvidedRun = Integer.parseInt(args[6])
    }

    // Allow for QADB override (usually meaning you're processing MC)
    Integer userProvidedOverride = 0;
    if (args.length < 8) {
        println("No indication of QADB provided. Will use QADB.");
    } else {
        userProvidedOverride = Integer.parseInt(args[7]);
    }

    // ~~~~~~~~~~~~~~~~ prepare physics analysis ~~~~~~~~~~~~~~~~ //
    // declare physics event variables
    int helicity;
	int num_pos, num_neg, num_neutrals; 
    int elec_detector, plus_detector, min_detector, nucl_detector;
    int nucl_pid;

	double elec_pid_chi2, plus_pid_chi2, min_pid_chi2, nucl_pid_chi2;

    double elec_px, elec_py, elec_pz, elec_p, elec_e, elec_theta, elec_phi;
	double min_px, min_py, min_pz, min_p, min_e, min_theta, min_phi;
    double plus_px, plus_py, plus_pz, plus_p, plus_e, plus_theta, plus_phi;
    double nucl_px, nucl_py, nucl_pz, nucl_p, nucl_e, nucl_theta, nucl_phi;

    double elec_vx, elec_vy, elec_vz;
    double min_vx, min_vy, min_vz;
    double plus_vx, plus_vy, plus_vz;
	double nucl_vx, nucl_vy, nucl_vz;
    
	double elec_e_pcal, elec_e_ecin, elec_e_ecout;
	double plus_e_pcal, plus_e_ecin, plus_e_ecout;
	double elec_m2_pcal_u, elec_m2_pcal_v, elec_m2_pcal_w;
	double elec_m2_ecin_u, elec_m2_ecin_v, elec_m2_ecin_w;
	double elec_m2_ecout_u, elec_m2_ecout_v, elec_m2_ecout_w;
	double posi_m2_pcal_u, posi_m2_pcal_v, posi_m2_pcal_w;
	double posi_m2_ecin_u, posi_m2_ecin_v, posi_m2_ecin_w;
	double posi_m2_ecout_u, posi_m2_ecout_v, posi_m2_ecout_w;

    // load kinematic fitter/PID
	GenericKinematicFitter fitter = new rho0_fitter(10.6041);

    // set filter for final states
	EventFilter filter = new EventFilter("11:"+p1_Str+":-"+p1_Str+":"+nucl_str+":X+:X-:Xn");

    QADB qa = new QADB("latest");
	qa.checkForDefect('TotalOutlier')    
	qa.checkForDefect('TerminalOutlier')
	qa.checkForDefect('MarginalOutlier')
	qa.checkForDefect('SectorLoss')
	qa.checkForDefect('LowLiveTime')
	qa.checkForDefect('Misc')
	qa.checkForDefect('ChargeHigh')
	qa.checkForDefect('ChargeNegative')
	qa.checkForDefect('ChargeUnknown')
	qa.checkForDefect('PossiblyNoBeam')
	qa.checkForDefect('BSAWrong')
	qa.checkForDefect('BSAUnknown')
	qa.checkForDefect('TSAWrong')
	qa.checkForDefect('TSAUnknown')
	qa.checkForDefect('DSAWrong')
	qa.checkForDefect('DSAUnknown')
    [ // list of runs with `Misc` that should be allowed, generally empty target etc for dilution factor calculations
  		6736, 6737, 6738,
  		6739, 6740, 6741, 6742, 6743, 6744, 6746, 6747,
  		6748, 6749, 6750, 6751, 6753, 6754, 6755, 6756,
  		6757, 											 // RGA runs FADC failure sector 6
  		16194, 16089, 16185, 16308, 16184, 16307, 16309, // RGC Su22 He/ET
  		16872, 16975, 									 // RGC Fa22 He/ET
  		17763, 17764, 17765, 17766, 17767, 17768,		 // RGC Sp23 He/ET
  		17179, 17180, 17181, 17182, 17183, 17188, 17189, // RICH off/partially down
  		17252,
		19660, 19662, 19663, 19664, 19665, 19666, 19823, 19824, 19825, // lumi scans
		19835 // RICH PMTs down
	].each{ run -> qa.allowMiscBit(run) }

    
}