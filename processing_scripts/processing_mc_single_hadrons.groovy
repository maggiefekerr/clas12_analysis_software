/*
 * author Timothy B. Hayward
 * 
 * SIDIS hadron mc (generated and reconstructed variables saved)
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

public static double phi_calculation (double x, double y) {
	// tracks are given with Cartesian values and so must be converted to cylindrical
	double phi = Math.toDegrees(Math.atan2(x,y));
	phi = phi - 90;
	if (phi < 0) {
		phi = 360 + phi;
	}
	phi = 360 - phi;
	return phi;	
}

public static double theta_calculation (double x, double y, double z) {
	// convert cartesian coordinates to polar angle
	double r = Math.pow(Math.pow(x,2)+Math.pow(y,2)+Math.pow(z,2),0.5);
	return (double) (180/Math.PI)*Math.acos(z/r);
}

public static void main(String[] args) {

	double scale = 3.0;

	// Start time
	long startTime = System.currentTimeMillis();

	// ~~~~~~~~~~~~~~~~ set up input paramaeters ~~~~~~~~~~~~~~~~ //

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

	// Set the PDG PID for p1 based on the provided 2nd argument or default to 211 (pi+)
	String p1_Str = args.length < 2 ? "211" : args[1];
	if (args.length < 2) println("WARNING: Specify a PDG PID for p1! Set to pi+ (211).");
	println("Set p1 PID = $p1_Str");
	int p1_int = p1_Str.toInteger(); // Convert p1_Str to integer

	// Set the output file name based on the provided 3rd argument or use the default name
	String output_file = args.length < 3 ? "hadron_dummy_out.txt" : args[2];
	if (args.length < 3) 
	    println("WARNING: Specify an output file name. Set to \"hadron_dummy_out.txt\".");
	File file = new File(output_file);
	file.delete();
	BufferedWriter writer = new BufferedWriter(new FileWriter(file));

	// Set the number of files to process based on the provided 4th argument
	// use the size of the hipo_list if no argument provided
	int n_files = args.length < 4 || Integer.parseInt(args[3]) > hipo_list.size()
	    ? hipo_list.size() : Integer.parseInt(args[3]);
	if (args.length < 4 || Integer.parseInt(args[3]) > hipo_list.size()) {
	    // Print warnings and information if the number of files is not specified or too large
	    println("WARNING: Number of files not specified or number too large.")
	    println("Setting # of files to be equal to number of files in the directory.");
	    println("There are $hipo_list.size files.");
	}

	int hadron_pair_counts = 0;
	GenericKinematicFitter research_fitter = new analysis_fitter(10.6041);
	GenericKinematicFitter mc_fitter = new monte_carlo_fitter(10.6041);
	EventFilter filter = new EventFilter("11:"+p1_Str+":"+":X+:X-:Xn");

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

		while(reader.hasEvent()){
			++num_events; 
			if (num_events%1000000 == 0) { 
				print("processed: "+num_events+" events.     ");
			}

			// get run and event numbers
			event = reader.getNextEvent();
		    int runnum = event.getBank("RUN::config").getInt('run',0);
		    int evnum = event.getBank("RUN::config").getInt('event',0);

		    PhysicsEvent research_Event = research_fitter.getPhysicsEvent(event);
		    PhysicsEvent mc_Event = mc_fitter.getPhysicsEvent(event);

			if (filter.isValid(research_Event)) {
			
				HipoDataBank recBank = (HipoDataBank) event.getBank("REC::Event");
				HipoDataBank lundBank = (HipoDataBank) event.getBank("MC::Lund");
				HipoDataBank mcBank = (HipoDataBank) event.getBank("MC::Particle");

				int num_p1 = research_Event.countByPid(p1_Str.toInteger()); 

				for (int current_p1 = 0; current_p1 < num_p1; current_p1++) {

					Particle exp_e = research_Event.getParticleByPid(11,0);
					Particle exp_p1 = research_Event.getParticleByPid(p1_Str.toInteger(),current_p1);

					BeamEnergy Eb = new BeamEnergy(runnum, false);
					BeamEnergy mc_Eb = new BeamEnergy(runnum, false);
					Hadron variables = new Hadron(event, research_Event, 
						p1_Str.toInteger(), current_p1, Eb.Eb());
					Hadron mc_variables = new Hadron(event, mc_Event, 
						p1_Str.toInteger(), current_p1, mc_Eb.Eb());

					if (variables.channel_test(variables)) {

						// lab kinematics data
						double e_p = variables.e_p();
						double e_theta = variables.e_theta();
						double e_phi = variables.e_phi();
						double vz_e = variables.vz_e();
						double p_p = variables.p_p();
						double p_theta = variables.p_theta();
						double p_phi = variables.p_phi();
						double vz_p = variables.vz_p();

						// lab kinematics MC
						double mc_e_p = mc_variables.e_p();
						double mc_e_theta = mc_variables.e_theta();
						double mc_e_phi = mc_variables.e_phi();
						double mc_vz_e = mc_variables.vz_e();
						double mc_p_p = mc_variables.p_p();
						double mc_p_theta = mc_variables.p_theta();
						double mc_p_phi = mc_variables.p_phi();
						double mc_vz_p = mc_variables.vz_p();

						// DIS variables data
						double Q2 = variables.Q2();
						double W = variables.W();
						double y = variables.y();
						double Mx = variables.Mx();
						double Mx2 = variables.Mx2();

						// DIS variables MC
						double mc_Q2 = mc_variables.Q2();
						double mc_W = mc_variables.W(); 
						double mc_y = mc_variables.y();
						double mc_Mx = mc_variables.Mx();
						double mc_Mx2 = mc_variables.Mx2();

						// SIDIS variables data
						double x = variables.x();
						double z = variables.z();
						double xF = variables.xF();
						double pT = variables.pT();
						double eta = variables.eta();
						double zeta = variables.zeta();

						// SIDIS variables MC
						double mc_x = mc_variables.x();
						double mc_z = mc_variables.z();
						double mc_xF = mc_variables.xF();
						double mc_pT = mc_variables.pT();
						double mc_eta = mc_variables.eta();
						double mc_zeta = mc_variables.zeta();

						// depolarization factors data
		                double Depolarization_A = variables.Depolarization_A();
		                double Depolarization_B = variables.Depolarization_B();
		                double Depolarization_C = variables.Depolarization_C();
		                double Depolarization_V = variables.Depolarization_V();
				    	double Depolarization_W = variables.Depolarization_W();

				    	// depolarization factors MC
		                double mc_Depolarization_A = mc_variables.Depolarization_A();
		                double mc_Depolarization_B = mc_variables.Depolarization_B();
		                double mc_Depolarization_C = mc_variables.Depolarization_C();
		                double mc_Depolarization_V = mc_variables.Depolarization_V();
				    	double mc_Depolarization_W = mc_variables.Depolarization_W();

						// angles data
						double trento_phi = variables.phi();

						// angles MC
						double mc_trento_phi = mc_variables.phi();

						boolean matching_e = false;
						boolean matching_p1 = false;

						int matching_p1_pid = 0;
						int mc_p1_parent_index = 0;
						for (int current_part = 0; current_part < lundBank.rows(); current_part++) {
							int pid = lundBank.getInt("pid", current_part);
							if (matching_p1) { continue; }
							double mc_px = lundBank.getFloat("px", current_part);
							double mc_py = lundBank.getFloat("py", current_part);
							double mc_pz = lundBank.getFloat("pz", current_part);

							double mc_phi = phi_calculation(mc_px, mc_py);
							double mc_theta_lab = theta_calculation(mc_px, mc_py, mc_pz);

							double exp_phi = phi_calculation(exp_p1.px(), exp_p1.py());
							double exp_theta = theta_calculation(exp_p1.px(), exp_p1.py(), 
								exp_p1.pz());

							matching_p1 = Math.abs(exp_phi - mc_phi) < scale*3.0 && 
								Math.abs(exp_theta - mc_theta_lab) < scale*1.0;
							if (matching_p1) {
								matching_p1_pid = pid;
								mc_p1_parent_index = lundBank.getInt("parent", current_part)-1;
							}
						}
						int mc_p1_parent = lundBank.getInt("pid", mc_p1_parent_index);


						matching_e = false;
						matching_p1 = false;

						int matching_e_pid = 0;
						int mc_e_parent_index = 0;
						for (int current_part = 0; current_part < mcBank.rows(); current_part++) {
							int pid = mcBank.getInt("pid", current_part);
							if (matching_e) { continue; }
							double mc_px = mcBank.getFloat("px", current_part);
							double mc_py = mcBank.getFloat("py", current_part);
							double mc_pz = mcBank.getFloat("pz", current_part);

							double mc_phi = phi_calculation(mc_px, mc_py);
							double mc_theta_lab = theta_calculation(mc_px, mc_py, mc_pz);

							double exp_phi = phi_calculation(exp_e.px(), exp_e.py());
							double exp_theta = theta_calculation(exp_e.px(), exp_e.py(), exp_e.pz());

							matching_e = Math.abs(exp_phi - mc_phi) < scale*3.0 && 
								Math.abs(exp_theta - mc_theta_lab) < scale*1.0;
							if (matching_e) {
								matching_e_pid = pid;
							}
						}

						matching_p1_pid = 0;
						mc_p1_parent_index = 0;
						for (int current_part = 0; current_part < mcBank.rows(); current_part++) {
							int pid = mcBank.getInt("pid", current_part);
							if (matching_p1) { continue; }
							double mc_px = mcBank.getFloat("px", current_part);
							double mc_py = mcBank.getFloat("py", current_part);
							double mc_pz = mcBank.getFloat("pz", current_part);

							double mc_phi = phi_calculation(mc_px, mc_py);
							double mc_theta_lab = theta_calculation(mc_px, mc_py, mc_pz);

							double exp_phi = phi_calculation(exp_p1.px(), exp_p1.py());
							double exp_theta = theta_calculation(exp_p1.px(), exp_p1.py(), 
								exp_p1.pz());

							matching_p1 = Math.abs(exp_phi - mc_phi) < scale*3.0 && 
								Math.abs(exp_theta - mc_theta_lab) < scale*1.0;
							if (matching_p1) {
								matching_p1_pid = pid;
							}
						}

						// Use a StringBuilder to append all data in a single call
		                StringBuilder line = new StringBuilder();
						line.append(e_p).append(" ")
							.append(mc_e_p).append(" ")
							.append(e_theta).append(" ")
							.append(mc_e_theta).append(" ")
							.append(e_phi).append(" ")
							.append(mc_e_phi).append(" ")
							.append(vz_e).append(" ")
							.append(mc_vz_e).append(" ")
							.append(p_p).append(" ")
							.append(mc_p_p).append(" ")
							.append(p_theta).append(" ")
							.append(mc_p_theta).append(" ")
							.append(p_phi).append(" ")
							.append(mc_p_phi).append(" ")
							.append(vz_p).append(" ")
							.append(mc_vz_p).append(" ")
							.append(Q2).append(" ")
							.append(mc_Q2).append(" ")
							.append(W).append(" ")
							.append(mc_W).append(" ")
							.append(Mx).append(" ")
							.append(mc_Mx).append(" ")
							.append(Mx2).append(" ")
							.append(mc_Mx2).append(" ")
							.append(x).append(" ")
							.append(mc_x).append(" ")
							.append(y).append(" ")
							.append(mc_y).append(" ")
							.append(z).append(" ")
							.append(mc_z).append(" ")
							.append(xF).append(" ")
							.append(mc_xF).append(" ")
							.append(pT).append(" ")
							.append(mc_pT).append(" ")
							.append(zeta).append(" ")
							.append(mc_zeta).append(" ")
							.append(eta).append(" ")
							.append(mc_eta).append(" ")
							.append(trento_phi).append(" ")
							.append(mc_trento_phi).append(" ")
							.append(Depolarization_A).append(" ")
							.append(mc_Depolarization_A).append(" ")
							.append(Depolarization_B).append(" ")
							.append(mc_Depolarization_B).append(" ")
							.append(Depolarization_C).append(" ")
							.append(mc_Depolarization_C).append(" ")
							.append(Depolarization_V).append(" ")
							.append(mc_Depolarization_V).append(" ")
							.append(Depolarization_W).append(" ")
							.append(mc_Depolarization_W).append(" ")
							.append(matching_e_pid).append(" ")
							.append(matching_p1_pid).append(" ")
							.append(mc_p1_parent).append("\n");

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
			}
		reader.close();
		}

		// Write any remaining lines in the batchLines StringBuilder to the file
		if (batchLines.length() > 0) {
		    file.append(batchLines.toString());
		    batchLines.setLength(0);
		}

		println(); println(); // NUMBERS ARE WRONG, NEED PHIS AND VERTICES
		print("1: e_p, 3: e_theta, 5: p1_p, 7: p1_theta, ");
		print("9: Q2, 11: W, 13: Mx, 15: x, 17: y, 19: z, 21: xF, 23: pT, 25: zeta, ");
		print("27: eta, 29: phi, ");
		print("31: Dep_A, 33: Dep_B, 35: Dep_C, 37: Dep_V, 39: Dep_W, ")
		print("41: matching e pid, 42: matching p1 pid, 43: p1 parent id\n");

		println(); println();
		println("Set p1 PID = "+p1_Str+"\n");
		println("output file is: "+file);
		println();
	}

}