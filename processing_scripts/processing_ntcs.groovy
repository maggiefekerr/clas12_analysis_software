/*
 * author Maggie F. E. Kerr
 * using processing_dvcs.groovy 
 * and derivative files of
 * Timothy B. Hayward as template
 * 
 * nTCS
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

    // Set the output file name based on the provided 2nd argument or use the default name
	String output_file = args.length < 2 ? "ntcs_dummy_out.txt" : args[1];
	if (args.length < 2) 
	    println("WARNING: Specify an output file name. Set to \"ntcs_dummy_out.txt\".");
	File file = new File(output_file);
	file.delete();
	BufferedWriter writer = new BufferedWriter(new FileWriter(file));

    // Set the number of files to process based on the provided 3rd argument or list size
	// If the argument is "0", default to the full list size
	int n_files = args.length < 3 || Integer.parseInt(args[2]) == 0 || Integer.parseInt(args[2]) > hipo_list.size()
	    ? hipo_list.size() : Integer.parseInt(args[2]);
	if (args.length < 3 || Integer.parseInt(args[2]) == 0 || Integer.parseInt(args[2]) > hipo_list.size()) {
	    // Print warnings and information if the number of files is not specified, set to 0, or too large
	    println("WARNING: Number of files not specified, set to 0, or number too large.")
	    println("Setting # of files to be equal to number of files in the directory.");
	    println("There are $hipo_list.size files.");
	}

    // Set the beam energy based on the provided 4th argument or default to 10.6
	double beam_energy = args.length < 4 ? 10.6 : Double.parseDouble(args[3]);
	if (args.length < 4) {
	    println("No beam energy provided, defaulting to 10.6 GeV.");
	}

    // Set the user-provided run number if available
	Integer userProvidedRun = null
	if (args.length < 5) {
	    println("Run number not provided, will pull from hipo files.")
	    println("Think carefully about this if you are processing MC.")
	} else {
		userProvidedRun = Integer.parseInt(args[4]);
	}

    // Allow for QADB override (usually meaning you're processing MC)
	Integer userProvidedOverride = 0;
	if (args.length < 6) {
		println("No indication of QADB provided. Will use QADB.");
	} else {
		userProvidedOverride = Integer.parseInt(args[5]);
	}

    // ~~~~~~~~~~~~~~~~ prepare physics analysis ~~~~~~~~~~~~~~~~ //

    // declare physics event variables

    int helicity;
	int num_pos, num_neg, num_neutrals; 
    int elec_detector, posi_detector, neut_detector;
	double elec_chi2, posi_chi2, neut_chi2;
	double elec_px, elec_py, elec_pz, elec_p, elec_e, elec_theta, elec_phi;
    double posi_px, posi_py, posi_pz, posi_p, posi_e, posi_theta, posi_phi;
    double neut_px, neut_py, neut_pz, neut_p, neut_e, neut_theta, neut_phi;
    double elec_vx, elec_vy, elec_vz;
    double posi_vx, posi_vy, posi_vz;
	double neut_vx, neut_vy, neut_vz; // will need to readjust these

    // load kinematic fitter/PID
	GenericKinematicFitter fitter = new ntcs_fitter(10.6041);

    // set filter for final states
	EventFilter filter = new EventFilter("11:-11:2112");

    // setup QA database
	QADB qa = new QADB("latest");
	qa.checkForDefect('TotalOutlier')    
	qa.checkForDefect('TerminalOutlier')
	qa.checkForDefect('MarginalOutlier')
	qa.checkForDefect('SectorLoss')
	// qa.checkForDefect('LowLiveTime')
	qa.checkForDefect('Misc')
	qa.checkForDefect('ChargeHigh')
	qa.checkForDefect('ChargeNegative')
	qa.checkForDefect('ChargeUnknown')
	qa.checkForDefect('PossiblyNoBeam')
	[ // list of runs with `Misc` that should be allowed, generally empty target etc for dilution factor calculations
  		6736, 6737, 6738,
  		6739, 6740, 6741, 6742, 6743, 6744, 6746, 6747,
  		6748, 6749, 6750, 6751, 6753, 6754, 6755, 6756,
  		6757, 											 // RGA runs FADC failure sector 6
  		16194, 16089, 16185, 16308, 16184, 16307, 16309, // RGC Su22 He/ET
  		16872, 16975, 									 // RGC Fa22 He/ET
  		17763, 17764, 17765, 17766, 17767, 17768,		 // RGC Sp23 He/ET
  		17179, 17180, 17181, 17182, 17183, 17188, 17189, // RICH off/partially down
  		17252
	].each{ run -> qa.allowMiscBit(run) }

    // create a StringBuilder for accumulating lines
	StringBuilder batchLines = new StringBuilder();

    int num_events = 0;
	int max_lines = 1000;
	int lineCount = 0;

    for (current_file in 0..<n_files) {
        // limit to a certain number of files defined by n_files
		println("\n Opening file "+Integer.toString(current_file+1)
			+" out of "+n_files+".\n"); 
		
		HipoDataSource reader = new HipoDataSource();
		reader.open(hipo_list[current_file]); // open next hipo file
		HipoDataEvent event = reader.getNextEvent();

        while (reader.hasEvent()) {
            ++num_events;
		    if (num_events % 500000 == 0) { // not necessary, just updates output
		        print("processed: " + num_events + " events. ");
		    }

            // get run and event numbers
		    event = reader.getNextEvent();
		    // collect info for QA
		    int runnum = userProvidedRun ?: event.getBank("RUN::config").getInt('run', 0);
		    int evnum = event.getBank("RUN::config").getInt('event', 0);

            PhysicsEvent research_Event = fitter.getPhysicsEvent(event);

            // do not use the qa if it is MC (runnum = 11)
			boolean process_event = filter.isValid(research_Event) && 
												  (runnum == 11 ||  // MC
												   userProvidedOverride == 1 ||
												   qa.pass(runnum, evnum)); // just using QADB for now, if ever decide to do RG-A will need to revisit this
			if (process_event) {
				// get # of particles 
		        int elec_num = research_Event.countByPid(11);
		        int posi_num = research_Event.countByPid(-11);
				int neut_num = research_Event.countByPid(2112);

				// supply runnum and boolean for radiative simulation or not
				BeamEnergy Eb = new BeamEnergy(research_Event, runnum, false);
				// Use the input beam energy if runnum == 11, otherwise use Eb.Eb()
				double energy = (runnum == 11) ? beam_energy : Eb.Eb();
				nTCSParticles variables = new nTCSParticles(event, research_Event, energy); 
				// this is the class for defining all relevant kinematic variables
				if (variables.channel_test(variables)) {
					helicity = variables.get_helicity(); // helicity of event
	                elec_detector = variables.get_elec_detector();
	                posi_detector = variables.get_posi_detector();
					neut_detector = variables.get_neut_detector();
	                num_pos = variables.get_num_pos();
	                num_neg = variables.get_num_neg();
	                num_neutrals = variables.get_num_neutrals();

					// pid elec_chi2
					elec_chi2 = variables.get_elec_chi2pid();
					posi_chi2 = variables.get_posi_chi2pid();
					System.out.println(posi_chi2);
					neut_chi2 = variables.get_neut_chi2pid();

					// lab kinematics
					elec_px    = variables.get_elec_px();
					elec_py    = variables.get_elec_py(); 
					elec_pz    = variables.get_elec_pz(); 
					elec_p     = variables.get_elec_p(); 
					elec_e     = variables.get_elec_e(); 
					elec_theta = variables.get_elec_theta();
					elec_phi   = variables.get_elec_phi();
					posi_px    = variables.get_posi_px();
					posi_py    = variables.get_posi_py(); 
					posi_pz    = variables.get_posi_pz(); 
					posi_p     = variables.get_posi_p(); 
					posi_e     = variables.get_posi_e(); 
					posi_theta = variables.get_posi_theta();
					posi_phi   = variables.get_posi_phi();
					neut_px    = variables.get_neut_px();
					neut_py    = variables.get_neut_py(); 
					neut_pz    = variables.get_neut_pz(); 
					neut_p     = variables.get_neut_p(); 
					neut_e     = variables.get_neut_e(); 
					neut_theta = variables.get_neut_theta();
					neut_phi   = variables.get_neut_phi();

					// vertices
					elec_vx = variables.get_elec_vx();
					elec_vy = variables.get_elec_vy();
					elec_vz = variables.get_elec_vz();
					posi_vx = variables.get_posi_vx();
					posi_vy = variables.get_posi_vy();
					posi_vz = variables.get_posi_vz();
					neut_vx = variables.get_neut_vx();
					neut_vy = variables.get_neut_vy();
					neut_vz = variables.get_neut_vz();

					// Use a StringBuilder to append all data in a single call
					StringBuilder line = new StringBuilder();
					line.append(num_pos).append(" ")
						.append(num_neg).append(" ")
						.append(num_neutrals).append(" ")
						.append(runnum).append(" ")
	                	.append(evnum).append(" ")
	                	.append(helicity).append(" ")
	                	.append(elec_detector).append(" ")
	                	.append(posi_detector).append(" ")
						.append(neut_detector).append(" ")
						.append(elec_chi2).append(" ")
						.append(posi_chi2).append(" ")
						.append(neut_chi2).append(" ")
	                	.append(elec_p).append(" ")
	                	.append(elec_theta).append(" ")
	                	.append(elec_phi).append(" ")
	                	.append(elec_vz).append(" ")
	                	.append(posi_p).append(" ")
	                	.append(posi_theta).append(" ")
	                	.append(posi_phi).append(" ")
	                	.append(posi_vz).append(" ")
	                	.append(neut_p).append(" ")
	                	.append(neut_theta).append(" ")
	                	.append(neut_phi).append(" ")
	                	.append(neut_vz).append("\n");

					// print(line)

					// Append the line to the batchLines StringBuilder
	                batchLines.append(line.toString());
	                lineCount++; // Increment the line count

	                // If the line count reaches 1000, write to the file and reset
	                if (lineCount >= max_lines) {
	                    file.append(batchLines.toString());
	                    batchLines.setLength(0);
	                    lineCount = 0;
	                }
				}
			}
			reader.close();
        }

		// Write any remaining lines in the batchLines StringBuilder to the file
		if (batchLines.length() > 0) {
		    file.append(batchLines.toString());
		    batchLines.setLength(0);
		}

		/*println("1: fiducial_status, 2: num_pos, 3: num_neg, 4: num_neutrals, " +
	    "5: runnum, 6: evnum, 7: helicity, 8: detector1, 9: detector2, 10: e_p, 11: e_theta, 12: e_phi, 13: vz_e, " +
	    "14: p1_p, 15: p1_theta, 16: p1_phi, 17: vz_p1, 18: p2_p, 19: p2_theta, 20: p2_phi, 21: vz_p2, " +
	    "22: open_angle_ep, 23: open_angle_ep1, 24: open_angle_ep2, 25: open_angle_p1p2, " +
	    "26: Q2, 27: W, 28: Mx2, 29: Mx2_1, 30: Mx2_2, 31: x, 32: t, 33: t1, 34: t2, 35: tmin, 36: y, 37: z, " +
	    "38: z1, 39: z2, 40: Mh, 41: xF, 42: xF1, 43: xF2, 44: pT, 45: pT1, 46: pT2, 47: pTpT, " +
	    "48: xi, 49: xi1, 50: xi2, 51: eta, 52: eta1, 53: eta2, 54: Delta_eta, 55: eta1_gN, 56: eta2_gN, " +
	    "57: phi1, 58: phi2, 59: Delta_phi, 60: phih, 61: phiR, 62: theta, " +
	    "63: DepA, 64: DepB, 65: DepC, 66: DepV, 67: DepW, 68: Emiss2, 69: theta_gamma_gamma, " +
	    "70: pTmiss");*/ // not introducing this line until we have all the variables we want to use

		println("Analyzing ntcs.");
		println("output text file is: $file");
	}

	writer.close();

	// End time
	long endTime = System.currentTimeMillis()
	// Calculate the elapsed time
	long elapsedTime = endTime - startTime
	// Print the elapsed time in milliseconds
	println("Elapsed time: ${elapsedTime} ms");
}