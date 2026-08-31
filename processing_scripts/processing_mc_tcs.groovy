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
    if (args.length < 6 || Integer.parseInt(args[5]) == 0) {
        println("Run number not provided, will pull from hipo files.")
        println("Think carefully about this if you are processing MC.")
    } else {
        userProvidedRun = Integer.parseInt(args[5])
    }

    println("Script will analyze generated banks and any reconstructed banks and save both.")

    // reconstructed variables
    int reconstructed; // variable to declare whether the generated event has reconstructed particles
    int nucl_pid;
    int num_pos, num_neg, num_neutrals; 
    int elec_detector, posi_detector, nucl_detector;
	double elec_chi2, posi_chi2, nucl_chi2;
	double elec_px, elec_py, elec_pz, elec_p, elec_e, elec_theta, elec_phi;
    double posi_px, posi_py, posi_pz, posi_p, posi_e, posi_theta, posi_phi;
    double nucl_px, nucl_py, nucl_pz, nucl_p, nucl_e, nucl_theta, nucl_phi;
    double elec_vx, elec_vy, elec_vz;
    double posi_vx, posi_vy, posi_vz;
	double nucl_vx, nucl_vy, nucl_vz;
	double elec_e_pcal, elec_e_ecin, elec_e_ecout;
	double posi_e_pcal, posi_e_ecin, posi_e_ecout;
	double elec_m2_pcal_u, elec_m2_pcal_v, elec_m2_pcal_w;
	double elec_m2_ecin_u, elec_m2_ecin_v, elec_m2_ecin_w;
	double elec_m2_ecout_u, elec_m2_ecout_v, elec_m2_ecout_w;
	double posi_m2_pcal_u, posi_m2_pcal_v, posi_m2_pcal_w;
	double posi_m2_ecin_u, posi_m2_ecin_v, posi_m2_ecin_w;
	double posi_m2_ecout_u, posi_m2_ecout_v, posi_m2_ecout_w;

    // "generated" variables (mc)
    double weight;
    double gen_elec_px, gen_elec_py, gen_elec_pz, gen_elec_p, gen_elec_e, gen_elec_theta, gen_elec_phi;
    double gen_posi_px, gen_posi_py, gen_posi_pz, gen_posi_p, gen_posi_e, gen_posi_theta, gen_posi_phi;
    double gen_nucl_px, gen_nucl_py, gen_nucl_pz, gen_nucl_p, gen_nucl_e, gen_nucl_theta, gen_nucl_phi;
    double gen_elec_vx, gen_elec_vy, gen_elec_vz;
    double gen_posi_vx, gen_posi_vy, gen_posi_vz;
	double gen_nucl_vx, gen_nucl_vy, gen_nucl_vz;

    // load kinematic fitter/PID
	GenericKinematicFitter rec_fitter = new tcs_fitter(10.6041);
    GenericKinematicFitter gen_fitter = new monte_carlo_fitter(10.6041);

    // set filter for final states
	EventFilter filter = new EventFilter("11:-11:"+nucl_str+":X+:X-:Xn");

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
            // instantiating variables to use -999 as a flag
            // reconstructed variables
            reconstructed = -999; // variable to declare whether the generated event has reconstructed particles
            nucl_pid = nucl_int;
            num_pos = -999;
            num_neg = -999;
            num_neutrals = -999; 
            elec_detector = -999; 
            posi_detector = -999;
            nucl_detector = -999;
            elec_chi2 = -999; 
            posi_chi2 = -999;
            nucl_chi2 = -999;

            elec_px = -999; 
            elec_py = -999; 
            elec_pz = -999; 
            elec_p = -999; 
            elec_e = -999; 
            elec_theta = -999; 
            elec_phi = -999;

            posi_px = -999; 
            posi_py = -999; 
            posi_pz = -999; 
            posi_p = -999; 
            posi_e = -999; 
            posi_theta = -999; 
            posi_phi = -999;
            
            nucl_px = -999; 
            nucl_py = -999; 
            nucl_pz = -999; 
            nucl_p = -999; 
            nucl_e = -999; 
            nucl_theta = -999; 
            nucl_phi = -999;

            elec_vx = -999; 
            elec_vy = -999; 
            elec_vz = -999;
            posi_vx = -999; 
            posi_vy = -999; 
            posi_vz = -999;
            nucl_vx = -999; 
            nucl_vy = -999; 
            nucl_vz = -999;

            elec_e_pcal = -999; 
            elec_e_ecin = -999; 
            elec_e_ecout = -999;
            posi_e_pcal = -999; 
            posi_e_ecin = -999; 
            posi_e_ecout = -999;

            elec_m2_pcal_u = -999;
            elec_m2_pcal_v = -999; 
            elec_m2_pcal_w = -999;
            elec_m2_ecin_u = -999;
            elec_m2_ecin_v = -999;
            elec_m2_ecin_w = -999;
            elec_m2_ecout_u = -999;
            elec_m2_ecout_v = -999;
            elec_m2_ecout_w = -999;

            posi_m2_pcal_u = -999;
            posi_m2_pcal_v = -999;
            posi_m2_pcal_w = -999;
            posi_m2_ecin_u = -999;
            posi_m2_ecin_v = -999;
            posi_m2_ecin_w = -999;
            posi_m2_ecout_u = -999;
            posi_m2_ecout_v = -999;
            posi_m2_ecout_w = -999;

            // generated variables
            weight = -999;

            gen_elec_px = -999;
            gen_elec_py = -999;
            gen_elec_pz = -999;
            gen_elec_p = -999;
            gen_elec_e = -999;
            gen_elec_theta = -999;
            gen_elec_phi = -999;

            gen_posi_px = -999;
            gen_posi_py = -999;
            gen_posi_pz = -999;
            gen_posi_p = -999;
            gen_posi_e = -999;
            gen_posi_theta = -999;
            gen_posi_phi = -999;

            gen_nucl_px = -999;
            gen_nucl_py = -999;
            gen_nucl_pz = -999;
            gen_nucl_p = -999;
            gen_nucl_e = -999;
            gen_nucl_theta = -999;
            gen_nucl_phi = -999;

            gen_elec_vx = -999;
            gen_elec_vy = -999;
            gen_elec_vz = -999;
            gen_posi_vx = -999;
            gen_posi_vy = -999;
            gen_posi_vz = -999;
            gen_nucl_vx = -999;
            gen_nucl_vy = -999;
            gen_nucl_vz = -999;

            ++num_events;
		    if (num_events % 500000 == 0) { // not necessary, just updates output
		        print("processed: " + num_events + " events. ");
		    }

            // get run and event numbers
		    event = reader.getNextEvent();
		    // collect info for QA
		    int runnum = userProvidedRun ?: event.getBank("RUN::config").getInt('run', 0);
		    int evnum = event.getBank("RUN::config").getInt('event', 0);

            PhysicsEvent rec_Event = rec_fitter.getPhysicsEvent(event);
		    PhysicsEvent gen_Event = gen_fitter.getPhysicsEvent(event);

            //if (runnum != 11) {
			//    throw new IllegalArgumentException("error: runnum != 11; this script is intended for use on MC (runnum ==11).")
			//}
            boolean process_event = filter.isValid(gen_Event)
			reconstructed = filter.isValid(rec_Event) ? 1 : 0;
			generated_cut = false; // not using right now

            if (true) {
                // get # of particles 
		        int elec_num = gen_Event.countByPid(11);
		        int posi_num = gen_Event.countByPid(-11);
				int nucl_num = gen_Event.countByPid(nucl_int);

                // supply runnum and boolean for radiative simulation or not
				BeamEnergy Eb = new BeamEnergy(gen_Event, runnum, false);
				// Use the input beam energy if runnum == 11, otherwise use Eb.Eb()
				//double energy = (runnum == 11) ? beam_energy : Eb.Eb();
                double energy = beam_energy
                TCSParticles variables = new TCSParticles(event, gen_Event, energy, nucl_int, nucl_str);
                generated_cut = variables.channel_test(variables);
                if (generated_cut) {
                    weight = variables.get_weight();

                    gen_elec_px = variables.get_elec_px();
                    gen_elec_py = variables.get_elec_py();
                    gen_elec_pz = variables.get_elec_pz();
                    gen_elec_p = variables.get_elec_p();
                    gen_elec_e = variables.get_elec_e();
                    gen_elec_theta = variables.get_elec_theta();
                    gen_elec_phi = variables.get_elec_phi();

                    gen_posi_px = variables.get_posi_px();
                    gen_posi_py = variables.get_posi_py();
                    gen_posi_pz = variables.get_posi_pz();
                    gen_posi_p = variables.get_posi_p();
                    gen_posi_e = variables.get_posi_e();
                    gen_posi_theta = variables.get_posi_theta();
                    gen_posi_phi = variables.get_posi_phi();

                    gen_nucl_px = variables.get_nucl_px();
                    gen_nucl_py = variables.get_nucl_py();
                    gen_nucl_pz = variables.get_nucl_pz();
                    gen_nucl_p = variables.get_nucl_p();
                    gen_nucl_e = variables.get_nucl_e();
                    gen_nucl_theta = variables.get_nucl_theta();
                    gen_nucl_phi = variables.get_nucl_phi();

                    gen_elec_vx = variables.get_elec_vx();
                    gen_elec_vy = variables.get_elec_vy();
                    gen_elec_vz = variables.get_elec_vz();
                    gen_posi_vx = variables.get_posi_vx();
                    gen_posi_vy = variables.get_posi_vy();
                    gen_posi_vz = variables.get_posi_vz();
                    gen_nucl_vx = variables.get_nucl_vx();
                    gen_nucl_vy = variables.get_nucl_vy();
                    gen_nucl_vz = variables.get_nucl_vz();
                }
            }

            if (process_event && reconstructed) {
                // get # of particles 
		        int elec_num = rec_Event.countByPid(11);
		        int posi_num = rec_Event.countByPid(-11);
				int nucl_num = rec_Event.countByPid(nucl_int);

                // supply runnum and boolean for radiative simulation or not
				BeamEnergy Eb = new BeamEnergy(rec_Event, runnum, false);
				// Use the input beam energy if runnum == 11, otherwise use Eb.Eb()
				double energy = (runnum == 11) ? beam_energy : Eb.Eb();
                TCSParticles variables = new TCSParticles(event, rec_Event, energy, nucl_int, nucl_str);
                if (variables.channel_test(variables)) {
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

                    // PCAL/ECAL readout
					elec_e_pcal = variables.get_elec_e_pcal();
					elec_e_ecin = variables.get_elec_e_ecin();
					elec_e_ecout = variables.get_elec_e_ecout();
					posi_e_pcal = variables.get_posi_e_pcal();
					posi_e_ecin = variables.get_posi_e_ecin();
					posi_e_ecout = variables.get_posi_e_ecout();
					elec_m2_pcal_u = variables.get_elec_m2_pcal_u();
					elec_m2_pcal_v = variables.get_elec_m2_pcal_v();
					elec_m2_pcal_w = variables.get_elec_m2_pcal_w();
					elec_m2_ecin_u = variables.get_elec_m2_ecin_u();
					elec_m2_ecin_v = variables.get_elec_m2_ecin_v();
					elec_m2_ecin_w = variables.get_elec_m2_ecin_w();
					elec_m2_ecout_u = variables.get_elec_m2_ecout_u();
					elec_m2_ecout_v = variables.get_elec_m2_ecout_v();
					elec_m2_ecout_w = variables.get_elec_m2_ecout_w();
					posi_m2_pcal_u = variables.get_posi_m2_pcal_u();
					posi_m2_pcal_v = variables.get_posi_m2_pcal_v();
					posi_m2_pcal_w = variables.get_posi_m2_pcal_w();
					posi_m2_ecin_u = variables.get_posi_m2_ecin_u();
					posi_m2_ecin_v = variables.get_posi_m2_ecin_v();
					posi_m2_ecin_w = variables.get_posi_m2_ecin_w();
					posi_m2_ecout_u = variables.get_posi_m2_ecout_u();
					posi_m2_ecout_v = variables.get_posi_m2_ecout_v();
					posi_m2_ecout_w = variables.get_posi_m2_ecout_w();
                }
            }

            if (generated_cut) {
                // Use a StringBuilder to append all data in a single call
	            StringBuilder line = new StringBuilder();
                // first the generated variables
                line.append(reconstructed).append(" ")
                    .append(weight).append(" ")
                    .append(gen_elec_p).append(" ")
                    .append(gen_elec_theta).append(" ")
                    .append(gen_elec_phi).append(" ")
                    .append(gen_elec_vz).append(" ")
                    .append(gen_posi_p).append(" ")
                    .append(gen_posi_theta).append(" ")
                    .append(gen_posi_phi).append(" ")
                    .append(gen_posi_vz).append(" ")
                    .append(gen_nucl_p).append(" ")
                    .append(gen_nucl_theta).append(" ")
                    .append(gen_nucl_phi).append(" ")
                    .append(gen_nucl_vz).append(" ")
                    .append(num_pos).append(" ")
                    .append(num_neg).append(" ")
                    .append(num_neutrals).append(" ")
                    .append(runnum).append(" ")
                    .append(evnum).append(" ")
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
                    .append(nucl_vz).append(" ")
                    .append(elec_e_pcal).append(" ")
                    .append(elec_e_ecin).append(" ")
                    .append(elec_e_ecout).append(" ")
                    .append(posi_e_pcal).append(" ")
                    .append(posi_e_ecin).append(" ")
                    .append(posi_e_ecout).append(" ")
                    .append(elec_m2_pcal_u).append(" ")
                    .append(elec_m2_pcal_v).append(" ")
                    .append(elec_m2_pcal_w).append(" ")
                    .append(elec_m2_ecin_u).append(" ")
                    .append(elec_m2_ecin_v).append(" ")
                    .append(elec_m2_ecin_w).append(" ")
                    .append(elec_m2_ecout_u).append(" ")
                    .append(elec_m2_ecout_v).append(" ")
                    .append(elec_m2_ecout_w).append(" ")
                    .append(posi_m2_pcal_u).append(" ")
                    .append(posi_m2_pcal_v).append(" ")
                    .append(posi_m2_pcal_w).append(" ")
                    .append(posi_m2_ecin_u).append(" ")
                    .append(posi_m2_ecin_v).append(" ")
                    .append(posi_m2_ecin_w).append(" ")
                    .append(posi_m2_ecout_u).append(" ")
                    .append(posi_m2_ecout_v).append(" ")
                    .append(posi_m2_ecout_w).append("\n");
                
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

        // Write any remaining lines in the batchLines StringBuilder to the file
		if (batchLines.length() > 0) {
		    file.append(batchLines.toString());
		    batchLines.setLength(0);
		}
        println("Analyzing tcs mc.");
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