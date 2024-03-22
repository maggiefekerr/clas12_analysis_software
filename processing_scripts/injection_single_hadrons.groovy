/*
 * author Timothy B. Hayward
 * 
 * SIDIS hadron 
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

def helicity_assignment(double Q2, double x, double PT, double z, double zeta, double xF, 
	double phi, double A, double B, double C, double V, double W) {
	double Pb = 0.83534; // injected beam polarization
	double Pt = 0.76200; // injected target polarization
	double Df = 0.1158; // injected dilution factor

	// injected asymmetry values, can depend on parameters or not

	// // TEST 1, null injection
	// double ALUsinphi = 0; 
	// double AULsinphi = 0;
	// double AULsin2phi = 0;
	// double ALL = 0;
	// double ALLcosphi = 0;

	// // TEST 2, single constant BSA
	// double ALUsinphi = -0.05; 
	// double AULsinphi = 0;
	// double AULsin2phi = 0;
	// double ALL = 0;
	// double ALLcosphi = 0;

	// // TEST 3, mixed constant injection
	// double ALUsinphi = -0.05; 
	// double AULsinphi = -0.10;
	// double AULsin2phi = 0.05;
	// double ALL = 0.40;
	// double ALLcosphi = 0.10;

	// // TEST 4, data-like xB injection
	// double ALUsinphi = -0.021; 
	// double AULsinphi = -0.025;
	// double AULsin2phi = -0.014;
	// double ALL = -0.046+2.017*x-2.986*x*x+2.761*x*x*x;
	// double ALLcosphi = 0.0191;

	// // TEST 5, data-like PT injection
	// double ALUsinphi = 0.002-0.056*PT; 
	// double AULsinphi = 0.009-0.088*PT;
	// double AULsin2phi = -0.013;
	// double ALL = 0.296;
	// double ALLcosphi = 0.100-0.020*PT;

	// // TEST 6, data-like xF injection
	// double ALUsinphi = -0.017+0.052*xF+0.103*xF*xF; 
	// double AULsinphi = -0.017+0.102*xF+0.204*xF*xF;
	// double AULsin2phi = -0.020-0.119*xF-0.248*xF*xF;
	// double ALL = 0.308+0.018*xF-0.250*xF*xF;
	// double ALLcosphi = 0.019+0.090*xF+0.109*xF*xF;

	// TEST "7", UU studies
	double AUUcosphi = -0.1;
	double AUUcos2phi = 0.05;
	double ALUsinphi = -0.01;
	double AULsinphi = 0;
	double AULsin2phi = 0;
	double ALL = 0;
	double ALLcosphi = 0;

	// TEST "8", UU studies
	double AUUcosphi = -0.15;
	double AUUcos2phi = -0.05;
	double ALUsinphi = -0.01;
	double AULsinphi = 0;
	double AULsin2phi = 0;
	double ALL = 0;
	double ALLcosphi = 0;

	// TEST "9", UU studies
	double AUUcosphi = -0.3*PT;
	double AUUcos2phi = 0.2*PT;
	double ALUsinphi = 0.002-0.056*PT;
	double AULsinphi = 0;
	double AULsin2phi = 0;
	double ALL = 0;
	double ALLcosphi = 0;


	int hb, ht;
	boolean weight_check = true;

	while(weight_check) {
		hb = new Random().nextBoolean() ? 1 : -1; // beam helicity
		ht = new Random().nextBoolean() ? 1 : -1; // target helicity
		double weight = 1 + 
			(V/A)*AUUcosphi*Math.cos(phi) +
			(B/A)*AUUcos2phi*Math.cos(2*phi) +
			hb*Pb*(W/A)*ALUsinphi*Math.sin(phi) + 
			ht*Pt*Df*(V/A)*AULsinphi*Math.sin(phi) +
			ht*Pt*Df*(B/A)*AULsin2phi*Math.sin(2*phi) + 
			hb*Pb*ht*Pt*Df*(C/A)*ALL + 
			hb*Pb*ht*Pt*Df*(W/A)*ALLcosphi*Math.cos(phi); 
		def randomValue = new Random().nextDouble() * 2;
		if (weight > randomValue) { weight_check = false; }
	}

	return [hb, ht]; // beam and target helicities
}

public static void main(String[] args) {

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

	// ~~~~~~~~~~~~~~~~ prepare physics analysis ~~~~~~~~~~~~~~~~ //

	// declare physics event variables
	int helicity;
	double e_p, e_theta, e_phi, p_phi, p_p, p_theta;
	double Q2, W, y, Mx, Mx2, x, z, xF, pT, eta, zeta, phi, vz_e, vz_p;
	double Depolarization_A, Depolarization_B, Depolarization_C;
	double Depolarization_V, Depolarization_W;

	// load my kinematic fitter/PID
	GenericKinematicFitter fitter = new analysis_fitter(10.6041); 
	GenericKinematicFitter mc_fitter = new monte_carlo_fitter(10.6041);
	// GenericKinematicFitter fitter = new event_builder_fitter(10.6041); 
	// GenericKinematicFitter fitter = new proton_energy_loss_corrections_fitter(10.6041); 
	
	// set filter for final states
	EventFilter filter = new EventFilter("11:"+p1_Str+":X+:X-:Xn"); 
	
	// setup QA database
	QADB qa = new QADB();

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
		    if (num_events % 100000 == 0) { // not necessary, just updates output
		        print("processed: " + num_events + " events. ");
		    }

		    // get run and event numbers
		    event = reader.getNextEvent();
		    // collect info for QA
		    int runnum = event.getBank("RUN::config").getInt('run', 0);
		    int evnum = event.getBank("RUN::config").getInt('event', 0);

		    PhysicsEvent research_Event = fitter.getPhysicsEvent(event);
		    PhysicsEvent mc_Event = mc_fitter.getPhysicsEvent(event);

		    // do not use the qa if it is MC (runnum = 11) 
		    // do not use the qa if the run is from RGC (until QA is produced!)
		    boolean process_event = filter.isValid(research_Event);

		    if (process_event) {

		    	HipoDataBank recBank = (HipoDataBank) event.getBank("REC::Event");
				HipoDataBank lundBank = (HipoDataBank) event.getBank("MC::Lund");
				HipoDataBank mcBank = (HipoDataBank) event.getBank("MC::Particle");

		        // get # of particles w/ pid1
		        int num_p1 = research_Event.countByPid(p1_Str.toInteger());
		        int mc_num_p1 = mc_Event.countByPid(p1_Str.toInteger()); 

		        Hadron mc_variables = new Hadron(event, mc_Event, p1_int, 0, 10.6);
		        mc_Q2 = mc_variables.Q2();
		        mc_x = mc_variables.x();
		        mc_pT = mc_variables.pT();
		        mc_z = mc_variables.z();
		        mc_zeta = mc_variables.zeta();
		        mc_xF = mc_variables.xF();
		        mc_phi = mc_variables.phi();
		        mc_Depolarization_A = mc_variables.Depolarization_A();
                mc_Depolarization_B = mc_variables.Depolarization_B();
                mc_Depolarization_C = mc_variables.Depolarization_C();
                mc_Depolarization_V = mc_variables.Depolarization_V();
		    	mc_Depolarization_W = mc_variables.Depolarization_W();
		        int[] helicities = helicity_assignment(mc_Q2, mc_x, mc_pT, mc_z, mc_zeta, mc_xF,
		        	mc_phi, mc_Depolarization_A, mc_Depolarization_B, mc_Depolarization_C,
		        	mc_Depolarization_V, mc_Depolarization_W);
		        int hb = helicities[0]; double ht = helicities[1];

		        // cycle over all hadrons
		        for (int current_p1 = 0; current_p1 < 1; current_p1++) { 

		            Hadron variables = new Hadron(event, research_Event,
		                    p1_int, current_p1, 10.6);
		            // this is my class for defining all relevant kinematic variables

		            if (variables.channel_test(variables)) {

		                // lab kinematics
		                e_p = variables.e_p(); // lab frame momentum
		                e_theta = variables.e_theta(); // lab polar angle
		                e_phi = variables.e_phi(); // lab azimuthal angle
		                p_phi = variables.p_phi(); // lab azimuthal angle
		                p_p = variables.p_p(); // lab momentum
		                p_theta = variables.p_theta(); // lab polar angle

		                // DIS variables
		                Q2 = variables.Q2(); // exchanged virtual photon energy
		                W = variables.W(); // hadronic mass
		                x = variables.x(); // Bjorken-x
		                y = variables.y(); // E_scat/E_beam
		                Mx = variables.Mx(); // missing mass
		                Mx2 = variables.Mx2(); // missing mass square

		                // SIDIS variables
		                z = variables.z(); // fractional hadron energy wrt virtual photon
		                xF = variables.xF(); // Feynman-x
		                pT = variables.pT(); // transverse momentum of hadron
		                eta = variables.eta(); // rapidity
		                zeta = variables.zeta(); // fractional longitudinal momentum of hadron

		                // angles
		                phi = variables.phi(); // trento phi of the hadron

		                // vertices
		                vz_e = variables.vz_e();
		                vz_p = variables.vz_p();

		                // depolarization factors
		                Depolarization_A = variables.Depolarization_A();
		                Depolarization_B = variables.Depolarization_B();
		                Depolarization_C = variables.Depolarization_C();
		                Depolarization_V = variables.Depolarization_V();
				    	Depolarization_W = variables.Depolarization_W();

		                // Use a StringBuilder to append all data in a single call
		                StringBuilder line = new StringBuilder();
		                line.append(runnum).append(" ")
		                	.append(evnum).append(" ")
		                	.append(hb).append(" ")
		                	.append(e_p).append(" ")
		                	.append(e_theta).append(" ")
		                	.append(e_phi).append(" ")
		                	.append(vz_e).append(" ")
		                	.append(p_p).append(" ")
		                	.append(p_theta).append(" ")
		                	.append(p_phi).append(" ")
		                	.append(vz_p).append(" ")
		                	.append(Q2).append(" ")
		                	.append(W).append(" ")
		                	.append(Mx).append(" ")
		                	.append(Mx2).append(" ")
		                	.append(x).append(" ")
		                	.append(y).append(" ")
		                	.append(z).append(" ")
		                	.append(xF).append(" ")
		                	.append(pT).append(" ")
		                	.append(zeta).append(" ")
		                	.append(ht).append(" ")
		                	.append(phi).append(" ")
		                	.append(Depolarization_A).append(" ")
		                    .append(Depolarization_B).append(" ")
		                    .append(Depolarization_C).append(" ")
		                    .append(Depolarization_V).append(" ")
		                    .append(Depolarization_W).append("\n");

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

		println("\n1:runnum, 2:evnum, 3:helicity, 4:e_p, 5:e_theta, 6:e_phi, 7:vz_e,"+
		"8:p_p, 9:p_theta, 10:p_phi, 11:vz_p, 12:Q2, 13:W, 14:Mx, 15: Mx2, 16:x, 17:y, 18:z,"+
		"19:xF, 20:pT, 21:zeta, 22:eta, 23:phi (trento), "+
		"24:DepA, 25:DepB, 26:DepC, 27:DepV, 28:DepW\n");

		println("Set p1 PID = $p1_Str");
		println("output file is: $file");
	}

	writer.close();

	// End time
	long endTime = System.currentTimeMillis()
	// Calculate the elapsed time
	long elapsedTime = endTime - startTime
	// Print the elapsed time in milliseconds
	println("Elapsed time: ${elapsedTime} ms");

}