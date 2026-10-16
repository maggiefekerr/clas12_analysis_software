/*
 * author Timothy B. Hayward & Maggie F. E. Kerr
 *
 * SIDIS hadron
 */

// import CLAS12 physics classes
import org.jlab.io.hipo.*
import org.jlab.io.base.DataEvent
import org.jlab.clas.physics.*
import org.jlab.clas12.physics.*

// import from hayward_coatjava_extensions
import extended_kinematic_fitters.*
import analyzers.*

// filetype for gathering files in directory
import groovy.io.FileType

// dilks CLAS QA analysis
import clasqa.QADB // access QADB

// -------------------- small helpers --------------------
public static double phi_calculation(double x, double y) {
    double phi = Math.toDegrees(Math.atan2(x, y))
    phi = phi - 90
    if (phi < 0) phi = 360 + phi
    phi = 360 - phi
    return phi
}

public static double theta_calculation(double x, double y, double z) {
    double r = Math.pow(Math.pow(x, 2) + Math.pow(y, 2) + Math.pow(z, 2), 0.5)
    return (double) (180 / Math.PI) * Math.acos(z / r)
}

// -------------------- main --------------------
public static void main(String[] args) {

    long startTime = System.currentTimeMillis()

    // ~~~~~~~~~~~~~~~~ set up input parameters ~~~~~~~~~~~~~~~~ //
    if (!args) {
        println("ERROR: Please enter a hipo file directory as the first argument")
        System.exit(0)
    }
    def hipo_list = []
    (args[0] as File).eachFileRecurse(FileType.FILES) { if (it.name.endsWith('.hipo')) hipo_list << it }

    // Set the PDG PID for the target particle based on the provided 2nd argument or default to 2212 (proton)
	String nucl_str = args.length >= 2 ? ((args[1].equals("2212") || args[1].equals("2112")) ? args[1] : "2212") : "2212";
    if (args.length < 2) println("WARNING: Specify either proton or neutron PDG PID for target type! Set to proton (2212).")
    if ((args[1] != "2212") && (args[1] !="2112")) println("WARNING: Specify either proton or neutron PDG PID for target type! Set to proton (2212).")
    println("Set PID for target type = $nucl_str")
    int nucl_int = nucl_str.toInteger()

    String p1_Str = args.length < 3 ? "211" : args[2]
    if (args.length < 3) println("WARNING: Specify a PDG PID for p1! Set to pi+ (211).")
    println("Set p1 PID = $p1_Str")
    int p1_int = p1_Str.toInteger()

    String output_file = args.length < 4 ? "hadron_dummy_out.txt" : args[3]
    if (args.length < 4) println('WARNING: Specify an output file name. Set to "hadron_dummy_out.txt".')
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

    // ~~~~~~~~~~~~~~~~ prepare physics analysis ~~~~~~~~~~~~~~~~ //
    int helicity, detector, sector
    int num_pos, num_neg, num_neutrals
    double e_p, e_theta, e_phi, p_phi, p_p, p_theta, open_angle
    double Q2, W, y, Mx2, x, t, tmin, z, xF, pT, eta, xi, phi, vz_e, vz_p
    double Depolarization_A, Depolarization_B, Depolarization_C, Depolarization_V, Depolarization_W
    int qadb_status
    int e_fd_cut, e_vertex_cut, e_sf_cut, e_diag_cut, e_pcal_fid_cut, e_dc_fid_cut
    int p_chi2pid_cut, p_vertex_cut, p_dc_fid_cut

    // Variables for fiducial, electron identification studies
    // Calorimeter variables
    double pcal_lu, pcal_lv, pcal_lw; // lu, lv, lw from pcal hit
    double ecin_lu, ecin_lv, ecin_lw; // lu, lv, lw from ecin hit
    double ecout_lu, ecout_lv, ecout_lw; // lu, lv, lw from ecout hit
    double pcal_e, ecin_e, ecout_e; // energy from pcal, ecin, ecout hit
    double pcal_x, pcal_y, pcal_z; // x, y, z from pcal hit
    double ecin_x, ecin_y, ecin_z; // x, y, z from ecin hit
    double ecout_x, ecout_y, ecout_z; // x, y, z from ecout hit
    // DC
    double e_dc_track_chi2, p_dc_track_chi2; // chi2 of track for different particles in DC
    int e_dc_track_ndf, p_dc_track_ndf; // ndf of track for different particles in DC
    double e_dc_edge_r1, e_dc_edge_r2, e_dc_edge_r3; // electron edges for DC regions
    double p_dc_edge_r1, p_dc_edge_r2, p_dc_edge_r3; // p edges for DC regions
    double e_dc_r1_x, e_dc_r1_y, e_dc_r1_z, e_dc_r2_x, e_dc_r2_y, e_dc_r2_z, e_dc_r3_x, e_dc_r3_y, e_dc_r3_z; // electron x, y, z for DC regions
    double p_dc_r1_x, p_dc_r1_y, p_dc_r1_z, p_dc_r2_x, p_dc_r2_y, p_dc_r2_z, p_dc_r3_x, p_dc_r3_y, p_dc_r3_z; // p x, y, z for DC regions
    // CVT
    double p_cvt_track_chi2; // chi2 of track for particle in CVT
    int p_cvt_track_ndf; // ndf of track for particle in CVT
    double p_cvt_edge_l1, p_cvt_edge_l3, p_cvt_edge_l5, p_cvt_edge_l7, p_cvt_edge_l12; // p edges for CVT layers
    double p_cvt_l1_x, p_cvt_l1_y, p_cvt_l1_z, p_cvt_l3_x, p_cvt_l3_y, p_cvt_l3_z, p_cvt_l5_x, p_cvt_l5_y, p_cvt_l5_z; // p x, y, z for CVT layers
    double p_cvt_l7_x, p_cvt_l7_y, p_cvt_l7_z, p_cvt_l12_x, p_cvt_l12_y, p_cvt_l12_z; // p x, y, z for CVT layers

    GenericKinematicFitter fitter = new analysis_fitter(10.6041);
    // GenericKinematicFitter fitter = new monte_carlo_fitter(10.6041);
    // GenericKinematicFitter fitter = new event_builder_fitter(10.6041);
    EventFilter filter = new EventFilter("11:" + p1_Str + ":X+:X-:Xn")

    QADB qa = new QADB("latest")
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
    [
        5046, 5047, 5051, 5128, 5129, 5130, 5158, 5159,
        5160, 5163, 5165, 5166, 5167, 5168, 5169, 5180,
        5181, 5182, 5183, 5400, 5448, 5495, 5496, 5505,
        5567, 5610, 5617, 5621, 5623, 6736, 6737, 6738,
        6739, 6740, 6741, 6742, 6743, 6744, 6746, 6747,
        6748, 6749, 6750, 6751, 6753, 6754, 6755, 6756,
        6757, 6378,
        16194, 16089, 16185, 16308, 16184, 16307, 16309,
        16872, 16975,
        17763, 17764, 17765, 17766, 17767, 17768,
        17179, 17180, 17181, 17182, 17183, 17188, 17189,
        17252,
        6378 											 // RG-B Helicity Flip but OK
    ].each { run -> qa.allowMiscBit(run) }

    StringBuilder batchLines = new StringBuilder()

    int num_events = 0
    int max_lines = 1000
    int lineCount = 0

    for (current_file in 0..<n_files) {
        println("\n Opening file " + Integer.toString(current_file + 1) + " out of " + n_files + ".\n")

        HipoDataSource reader = new HipoDataSource()
        reader.open(hipo_list[current_file])

        HipoDataEvent event = reader.getNextEvent()

        while (reader.hasEvent()) {
            ++num_events
            if (num_events % 500000 == 0) print("processed: " + num_events + " events. ")

            event = reader.getNextEvent()
            int runnum = event.getBank("RUN::config").getInt('run', 0)
            if ((userProvidedRun!=null)&&(userProvidedRun!=runnum)) {
                continue
            }
            // if (runnum > 16600 && runnum < 16700) break // Hall C bleedthrough
            int evnum = event.getBank("RUN::config").getInt('event', 0)
            // print(runnum+" "+evnum+"\n")
            PhysicsEvent research_Event = fitter.getPhysicsEvent(event)

            boolean process_event = filter.isValid(research_Event)
            
            // processing all events but recording the qadb status
            if (qa.pass(runnum,evnum)) { qadb_status = 1 }
            else                       { qadb_status = -1}

            // --- Toggle here ---
            // false -> baseline (no inverse-ISR)
            // true  -> apply inverse-ISR (subtract R from q inside TwoParticles)
            BeamEnergy Eb = new BeamEnergy(research_Event, runnum, /*isRadiative=*/false)

            double Ebeam = (runnum == 11) ? beam_energy : Eb.Eb()

            // If radiative sampling requested, prepare one-shot correction for TwoParticles
            boolean applyInverseISR = Eb.isRadiativeApplied() && Eb.getEgammaGeV() > 0.0
            double Egamma = 0.0, isrTheta = 0.0, isrPhi = 0.0

            if (process_event) {
                int num_p1 = research_Event.countByPid(p1_Str.toInteger())

                for (int current_p1 = 0; current_p1 < num_p1; current_p1++) {

                	if (applyInverseISR) {
		                Egamma  = Eb.getEgammaGeV()
		                isrTheta = analyzers.ISRThetaKernel.sampleThetaRad(Egamma) // radians
		                isrPhi   = 2.0 * Math.PI * Math.random()
		                TwoParticles.setNextInverseISRPhoton(Egamma, isrTheta, isrPhi)
		            } // endif
                    TwoParticles variables = new TwoParticles(event, research_Event, p1_int, current_p1, Ebeam, nucl_int, nucl_str)

                    if (TwoParticles.channel_test(variables)) {
                        fiducial_status = variables.get_fiducial_status()
                        helicity = variables.get_helicity()
                        detector = variables.get_detector()
                        sector = variables.get_sector()
                        num_pos = variables.get_num_pos()
                        num_neg = variables.get_num_neg()
                        num_neutrals = variables.get_num_neutrals()

                        // lab kinematics
                        e_p = variables.e_p()
                        e_theta = variables.e_theta()
                        e_phi = variables.e_phi()
                        p_phi = variables.p_phi()
                        p_p = variables.p_p()
                        p_theta = variables.p_theta()
                        open_angle = variables.open_angle()

                        // DIS/SIDIS (already corrected if toggle on)
                        Q2 = variables.Q2()
                        W = variables.W()
                        x = variables.x()
                        t = variables.t()
                        tmin = variables.tmin()
                        y = variables.y()
                        Mx2 = variables.Mx2()
                        if (Mx2 > 2) { continue; } // added to reduce size of files
                        z = variables.z()
                        xF = variables.xF()
                        pT = variables.pT()
                        eta = variables.eta()
                        xi = variables.xi()
                        phi = variables.phi()
                        vz_e = variables.vz_e()
                        vz_p = variables.vz_p()
                        Depolarization_A = variables.Depolarization_A()
                        Depolarization_B = variables.Depolarization_B()
                        Depolarization_C = variables.Depolarization_C()
                        Depolarization_V = variables.Depolarization_V()
                        Depolarization_W = variables.Depolarization_W()

                        // status of different cuts
                        e_fd_cut = variables.e_fd_cut()
                        e_vertex_cut = variables.e_vertex_cut()
                        e_sf_cut = variables.e_sf_cut()
                        e_diag_cut = variables.e_diag_cut()
                        e_pcal_fid_cut = variables.e_pcal_fid_cut()
                        e_dc_fid_cut = variables.e_dc_fid_cut()
                        p_chi2pid_cut = variables.p_chi2pid_cut()
                        p_vertex_cut = variables.p_vertex_cut()
                        p_dc_fid_cut = variables.p_dc_fid_cut()

                        // variables for fiducial and pid cuts
                        // Calorimeter
                        pcal_lu = variables.pcal_lu();
                        pcal_lv = variables.pcal_lv();
                        pcal_lw = variables.pcal_lw();
                        ecin_lu = variables.ecin_lu();
                        ecin_lv = variables.ecin_lv();
                        ecin_lw = variables.ecin_lw();
                        ecout_lu = variables.ecout_lu();
                        ecout_lv = variables.ecout_lv();
                        ecout_lw = variables.ecout_lw();
                        pcal_e = variables.pcal_e();
                        ecin_e = variables.ecin_e();
                        ecout_e = variables.ecout_e();
                        pcal_x = variables.pcal_x();
                        pcal_y = variables.pcal_y();
                        pcal_z = variables.pcal_z();
                        ecin_x = variables.ecin_x();
                        ecin_y = variables.ecin_y();
                        ecin_z = variables.ecin_z();
                        ecout_x = variables.ecout_x();
                        ecout_y = variables.ecout_y();
                        ecout_z = variables.ecout_z();
                        // DC
                        e_dc_track_chi2 = variables.e_dc_track_chi2();
                        p_dc_track_chi2 = variables.p_dc_track_chi2();
                        e_dc_track_ndf = variables.e_dc_track_ndf();
                        p_dc_track_ndf = variables.p_dc_track_ndf();
                        e_dc_edge_r1 = variables.e_dc_edge_r1();
                        e_dc_edge_r2 = variables.e_dc_edge_r2();
                        e_dc_edge_r3 = variables.e_dc_edge_r3();
                        p_dc_edge_r1 = variables.p_dc_edge_r1();
                        p_dc_edge_r2 = variables.p_dc_edge_r2();
                        p_dc_edge_r3 = variables.p_dc_edge_r3();
                        e_dc_r1_x = variables.e_dc_r1_x();
                        e_dc_r1_y = variables.e_dc_r1_y();
                        e_dc_r1_z = variables.e_dc_r1_z();
                        e_dc_r2_x = variables.e_dc_r2_x();
                        e_dc_r2_y = variables.e_dc_r2_y();
                        e_dc_r2_z = variables.e_dc_r2_z();
                        e_dc_r3_x = variables.e_dc_r3_x();
                        e_dc_r3_y = variables.e_dc_r3_y();
                        e_dc_r3_z = variables.e_dc_r3_z();
                        p_dc_r1_x = variables.p_dc_r1_x();
                        p_dc_r1_y = variables.p_dc_r1_y();
                        p_dc_r1_z = variables.p_dc_r1_z();
                        p_dc_r2_x = variables.p_dc_r2_x();
                        p_dc_r2_y = variables.p_dc_r2_y();
                        p_dc_r2_z = variables.p_dc_r2_z();
                        p_dc_r3_x = variables.p_dc_r3_x();
                        p_dc_r3_y = variables.p_dc_r3_y();
                        p_dc_r3_z = variables.p_dc_r3_z();
                        // CVT
                        p_cvt_track_chi2 = variables.p_cvt_track_chi2();
                        p_cvt_track_ndf = variables.p_cvt_track_ndf();
                        p_cvt_edge_l1 = variables.p_cvt_edge_l1();
                        p_cvt_edge_l3 = variables.p_cvt_edge_l3();
                        p_cvt_edge_l5 = variables.p_cvt_edge_l5();
                        p_cvt_edge_l7 = variables.p_cvt_edge_l7();
                        p_cvt_edge_l12 = variables.p_cvt_edge_l12();
                        p_cvt_l1_x = variables.p_cvt_l1_x();
                        p_cvt_l1_y = variables.p_cvt_l1_y();
                        p_cvt_l1_z = variables.p_cvt_l1_z();
                        p_cvt_l3_x = variables.p_cvt_l3_x();
                        p_cvt_l3_y = variables.p_cvt_l3_y();
                        p_cvt_l3_z = variables.p_cvt_l3_z();
                        p_cvt_l5_x = variables.p_cvt_l5_x();
                        p_cvt_l5_y = variables.p_cvt_l5_y();
                        p_cvt_l5_z = variables.p_cvt_l5_z();
                        p_cvt_l7_x = variables.p_cvt_l7_x();
                        p_cvt_l7_y = variables.p_cvt_l7_y();
                        p_cvt_l7_z = variables.p_cvt_l7_z();
                        p_cvt_l12_x = variables.p_cvt_l12_x();
                        p_cvt_l12_y = variables.p_cvt_l12_y();
                        p_cvt_l12_z = variables.p_cvt_l12_z();

                        // Output the sampled R (or zeros if baseline)
                        double isrTheta_deg = Math.toDegrees(isrTheta)
                        double isrPhi_deg   = Math.toDegrees(isrPhi)

                        StringBuilder line = new StringBuilder()
                        line.append(fiducial_status).append(" ")
                            .append(num_pos).append(" ")
                            .append(num_neg).append(" ")
                            .append(num_neutrals).append(" ")
                            .append(runnum).append(" ")
                            .append(evnum).append(" ")
                            .append(helicity).append(" ")
                            .append(detector).append(" ")
                            .append(e_p).append(" ")
                            .append(e_theta).append(" ")
                            .append(e_phi).append(" ")
                            .append(vz_e).append(" ")
                            .append(p_p).append(" ")
                            .append(p_theta).append(" ")
                            .append(p_phi).append(" ")
                            .append(vz_p).append(" ")
                            .append(open_angle).append(" ")
                            .append(applyInverseISR ? Egamma : 0.0).append(" ")
                            .append(applyInverseISR ? isrTheta_deg : 0.0).append(" ")
                            .append(applyInverseISR ? isrPhi_deg   : 0.0).append(" ")
                            .append(Q2).append(" ")
                            .append(W).append(" ")
                            .append(Mx2).append(" ")
                            .append(x).append(" ")
                            .append(t).append(" ")
                            .append(tmin).append(" ")
                            .append(y).append(" ")
                            .append(z).append(" ")
                            .append(xF).append(" ")
                            .append(pT).append(" ")
                            .append(xi).append(" ")
                            .append(eta).append(" ")
                            .append(phi).append(" ")
                            .append(Depolarization_A).append(" ")
                            .append(Depolarization_B).append(" ")
                            .append(Depolarization_C).append(" ")
                            .append(Depolarization_V).append(" ")
                            .append(Depolarization_W).append(" ")
                            .append(qadb_status).append(" ")
                            .append(e_fd_cut).append(" ")
                            .append(e_vertex_cut).append(" ")
                            .append(e_sf_cut).append(" ")
                            .append(e_diag_cut).append(" ")
                            .append(e_pcal_fid_cut).append(" ")
                            .append(e_dc_fid_cut).append(" ")
                            .append(p_chi2pid_cut).append(" ")
                            .append(p_vertex_cut).append(" ")
                            .append(p_dc_fid_cut).append(" ")
                            .append(e_chi2).append(" ")
                            .append(p_chi2).append(" ")
                            .append(pcal_lu).append(" ")
                            .append(pcal_lv).append(" ")
                            .append(pcal_lw).append(" ")
                            .append(ecin_lu).append(" ")
                            .append(ecin_lv).append(" ")
                            .append(ecin_lw).append(" ")
                            .append(ecout_lu).append(" ")
                            .append(ecout_lv).append(" ")
                            .append(ecout_lw).append(" ")
                            .append(pcal_e).append(" ")
                            .append(ecin_e).append(" ")
                            .append(ecout_e).append(" ")
                            .append(pcal_x).append(" ")
                            .append(pcal_y).append(" ")
                            .append(pcal_z).append(" ")
                            .append(ecin_x).append(" ")
                            .append(ecin_y).append(" ")
                            .append(ecin_z).append(" ")
                            .append(ecout_x).append(" ")
                            .append(ecout_y).append(" ")
                            .append(ecout_z).append(" ")
                            .append(e_dc_track_chi2).append(" ")
                            .append(p_dc_track_chi2).append(" ")
                            .append(e_dc_track_ndf).append(" ")
                            .append(p_dc_track_ndf).append(" ")
                            .append(e_dc_edge_r1).append(" ")
                            .append(e_dc_edge_r2).append(" ")
                            .append(e_dc_edge_r3).append(" ")
                            .append(p_dc_edge_r1).append(" ")
                            .append(p_dc_edge_r2).append(" ")
                            .append(p_dc_edge_r3).append(" ")
                            .append(e_dc_r1_x).append(" ")
                            .append(e_dc_r1_y).append(" ")
                            .append(e_dc_r1_z).append(" ")
                            .append(e_dc_r2_x).append(" ")
                            .append(e_dc_r2_y).append(" ")
                            .append(e_dc_r2_z).append(" ")
                            .append(e_dc_r3_x).append(" ")
                            .append(e_dc_r3_y).append(" ")
                            .append(e_dc_r3_z).append(" ")
                            .append(p_dc_r1_x).append(" ")
                            .append(p_dc_r1_y).append(" ")
                            .append(p_dc_r1_z).append(" ")
                            .append(p_dc_r2_x).append(" ")
                            .append(p_dc_r2_y).append(" ")
                            .append(p_dc_r2_z).append(" ")
                            .append(p_dc_r3_x).append(" ")
                            .append(p_dc_r3_y).append(" ")
                            .append(p_dc_r3_z).append(" ")
                            .append(p_cvt_track_chi2).append(" ")
                            .append(p_cvt_track_ndf).append(" ")
                            .append(p_cvt_edge_l1).append(" ")
                            .append(p_cvt_edge_l3).append(" ")
                            .append(p_cvt_edge_l5).append(" ")
                            .append(p_cvt_edge_l7).append(" ")
                            .append(p_cvt_edge_l12).append(" ")
                            .append(p_cvt_l1_x).append(" ")
                            .append(p_cvt_l1_y).append(" ")
                            .append(p_cvt_l1_z).append(" ")
                            .append(p_cvt_l3_x).append(" ")
                            .append(p_cvt_l3_y).append(" ")
                            .append(p_cvt_l3_z).append(" ")
                            .append(p_cvt_l5_x).append(" ")
                            .append(p_cvt_l5_y).append(" ")
                            .append(p_cvt_l5_z).append(" ")
                            .append(p_cvt_l7_x).append(" ")
                            .append(p_cvt_l7_y).append(" ")
                            .append(p_cvt_l7_z).append(" ")
                            .append(p_cvt_l12_x).append(" ")
                            .append(p_cvt_l12_y).append(" ")
                            .append(p_cvt_l12_z).append("\n");
                        batchLines.append(line.toString())
                        lineCount++

                        if (lineCount >= 1000) {
                            file.append(batchLines.toString())
                            batchLines.setLength(0)
                            lineCount = 0
                        }
                    } // channel test
                } // loop hadrons
            } // process_event
        } // while
        reader.close()

        if (batchLines.length() > 0) {
            file.append(batchLines.toString())
            batchLines.setLength(0)
        }

        println(
          "1:  fiducial_status,  " +
          "2:  num_pos,          " +
          "3:  num_neg,          " +
          "4:  num_neutrals,     " +
          "5:  runnum,           " +
          "6:  evnum,            " +
          "7:  helicity,         " +
          "8:  detector,         " +
          "9:  e_p,              " +
          "10: e_theta,          " +
          "11: e_phi,            " +
          "12: vz_e,             " +
          "13: p_p,              " +
          "14: p_theta,          " +
          "15: p_phi,            " +
          "16: vz_p,             " +
          "17: open_angle,       " +
          "18: Egamma (GeV),     " +
          "19: isrTheta (deg),   " +
          "20: isrPhi (deg),     " +
          "21: Q2,               " +
          "22: W,                " +
          "23: Mx2,              " +
          "24: x,                " +
          "25: t,                " +
          "26: tmin,             " +
          "27: y,                " +
          "28: z,                " +
          "29: xF,               " +
          "30: pT,               " +
          "31: xi,               " +
          "32: eta,              " +
          "33: phi (trento),     " +
          "34: DepA,             " +
          "35: DepB,             " +
          "36: DepC,             " +
          "37: DepV,             " +
          "38: DepW              " +
          "39: sector            " +
          "40: qadb_status       " +
          "41: e_fd_cut          " +
          "42: e_vertex_cut      " +
          "43: e_sf_cut          " +
          "44: e_diag_cut        " +
          "45: e_pcal_fid_cut    " +
          "46: e_dc_fid_cut      " +
          "47: p_chi2pid_cut     " +
          "48: p_vertex_cut      " +
          "49: p_dc_fid_cut,     + values for fiducial cuts"
        )

        println("Set p1 PID = $p1_Str")
        println("output text file is: $file")
    } // files

    writer.close()
    long endTime = System.currentTimeMillis()
    println("Elapsed time: ${endTime - startTime} ms")
}