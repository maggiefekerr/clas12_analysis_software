/*
 * author Maggie F. E. Kerr
 * using processing_dvcs.groovy 
 * and derivative files of
 * Timothy B. Hayward as template
 * 
 * TCS
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
    println("Set PID for TCS type = $nucl_str")
    int nucl_int = nucl_str.toInteger()

    String output_file = args.length < 3 ? "tcs_dummy_out.txt" : args[2]
    if (args.length < 3) println('WARNING: Specify an output file name. Set to "tcs_dummy_out.txt".')
    File file = new File(output_file)
    file.delete()
    BufferedWriter writer = new BufferedWriter(new FileWriter(file))

    int n_files = args.length < 4 || Integer.parseInt(args[3]) == 0 || Integer.parseInt(args[3]) > hipo_list.size()
        ? hipo_list.size() : Integer.parseInt(args[3])
    if (args.length < 4 || Integer.parseInt(args[3]) == 0 || Integer.parseInt(args[3]) > hipo_list.size()) {
        println("WARNING: Number of files not specified, set to 0, or number too large.")
        println("Setting # of files to be equal to number of files in the directory.")
        println("There are $hipo_list.size files.")
    }

    double beam_energy = args.length < 5 ? 10.6 : Double.parseDouble(args[4])
    if (args.length < 5) {
        println("No beam energy provided, defaulting to 10.6 GeV.")
        println("All MC will use 10.604 GeV. You must manually enter a beam energy to change this.")
    }

    Integer userProvidedRun = null
    if (args.length < 6) {
        println("Run number not provided, will pull from hipo files.")
        println("Think carefully about this if you are processing MC.")
    } else {
        userProvidedRun = Integer.parseInt(args[5])
    }

    // Allow for QADB override (usually meaning you're processing MC)
    Integer userProvidedOverride = 0;
    if (args.length < 7) {
        println("No indication of QADB provided. Will use QADB.");
    } else {
        userProvidedOverride = Integer.parseInt(args[6]);
    }

    // ~~~~~~~~~~~~~~~~ prepare physics analysis ~~~~~~~~~~~~~~~~ //
    // declare physics event variables
    int helicity;
	int num_pos, num_neg, num_neutrals; 
    int elec_detector, posi_detector, nucl_detector;
    int nucl_pid;
	double elec_chi2, posi_chi2, nucl_chi2;
	double elec_px, elec_py, elec_pz, elec_p, elec_e, elec_theta, elec_phi;
    double posi_px, posi_py, posi_pz, posi_p, posi_e, posi_theta, posi_phi;
    double nucl_px, nucl_py, nucl_pz, nucl_p, nucl_e, nucl_theta, nucl_phi;
    double elec_vx, elec_vy, elec_vz;
    double posi_vx, posi_vy, posi_vz;
	double nucl_vx, nucl_vy, nucl_vz;

    // load kinematic fitter/PID
	GenericKinematicFitter fitter = new tcs_fitter(10.6041);

    // set filter for final states
	EventFilter filter = new EventFilter("11:-11:"+nucl_str+":X+:X-:Xn");

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
				int nucl_num = research_Event.countByPid(nucl_int);

                // supply runnum and boolean for radiative simulation or not
				BeamEnergy Eb = new BeamEnergy(research_Event, runnum, false);
				// Use the input beam energy if runnum == 11, otherwise use Eb.Eb()
				double energy = (runnum == 11) ? beam_energy : Eb.Eb();
                TCSParticles variables = new TCSParticles(event, research_Event, energy, nucl_int, nucl_str);
                // this is the class for defining all relevant kinematic variables
                if (variables.channel_test(variables)) {
					helicity = variables.get_helicity(); // helicity of event
	                elec_detector = variables.get_elec_detector();
	                posi_detector = variables.get_posi_detector();
					nucl_detector = variables.get_nucl_detector();
	                num_pos = variables.get_num_pos();
	                num_neg = variables.get_num_neg();
	                num_neutrals = variables.get_num_neutrals();
                    nucl_pid = nucl_int;

                    // pid chi2
					elec_chi2 = variables.get_elec_chi2pid();
					posi_chi2 = variables.get_posi_chi2pid();
					nucl_chi2 = variables.get_nucl_chi2pid();

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
					nucl_px    = variables.get_nucl_px();
					nucl_py    = variables.get_nucl_py(); 
					nucl_pz    = variables.get_nucl_pz(); 
					nucl_p     = variables.get_nucl_p(); 
					nucl_e     = variables.get_nucl_e(); 
					nucl_theta = variables.get_nucl_theta();
					nucl_phi   = variables.get_nucl_phi();

                    // vertices
					elec_vx = variables.get_elec_vx();
					elec_vy = variables.get_elec_vy();
					elec_vz = variables.get_elec_vz();
					posi_vx = variables.get_posi_vx();
					posi_vy = variables.get_posi_vy();
					posi_vz = variables.get_posi_vz();
					nucl_vx = variables.get_nucl_vx();
					nucl_vy = variables.get_nucl_vy();
					nucl_vz = variables.get_nucl_vz();

                    // Use a StringBuilder to append all data in a single call
					StringBuilder line = new StringBuilder();
					line.append(num_pos).append(" ")
						.append(num_neg).append(" ")
						.append(num_neutrals).append(" ")
						.append(runnum).append(" ")
	                	.append(evnum).append(" ")
	                	.append(helicity).append(" ")
                        .append(nucl_pid).append(" ")
	                	.append(elec_detector).append(" ")
	                	.append(posi_detector).append(" ")
						.append(nucl_detector).append(" ")
						.append(elec_chi2).append(" ")
						.append(posi_chi2).append(" ")
						.append(nucl_chi2).append(" ")
	                	.append(elec_p).append(" ")
	                	.append(elec_theta).append(" ")
	                	.append(elec_phi).append(" ")
	                	.append(elec_vz).append(" ")
	                	.append(posi_p).append(" ")
	                	.append(posi_theta).append(" ")
	                	.append(posi_phi).append(" ")
	                	.append(posi_vz).append(" ")
	                	.append(nucl_p).append(" ")
	                	.append(nucl_theta).append(" ")
	                	.append(nucl_phi).append(" ")
	                	.append(nucl_vz).append("\n");
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

        println("Analyzing tcs.");
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