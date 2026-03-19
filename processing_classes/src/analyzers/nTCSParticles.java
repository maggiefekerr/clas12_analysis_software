 package analyzers;

/**
 *
 * @author Maggie F. E. Kerr
 */
import extended_kinematic_fitters.fiducial_cuts;
import extended_kinematic_fitters.generic_tests;
import extended_kinematic_fitters.momentum_corrections;
import org.jlab.clas.physics.Particle;
import org.jlab.clas.physics.PhysicsEvent;
import org.jlab.io.base.DataEvent;
import org.jlab.io.hipo.HipoDataBank;
import org.jlab.clas.physics.*;

public class nTCSParticles {

    protected byte helicity;
    protected int runnum;

    protected int elec_detector = -1;
    protected int posi_detector = -1;
    protected int neut_detector = -1;

    protected int num_electrons, num_piplus, num_piminus, num_kplus, num_kminus, num_protons, num_neutrons, num_particles;
    protected int num_pos, num_neg, num_neutrals;
    protected int num_positrons, num_antiprotons, num_antineutrons;

    protected double elec_chi2, posi_chi2, neut_chi2; // pid chi2 values

    protected double elec_px, elec_py, elec_pz, elec_p, elec_e, elec_theta, elec_phi; // electron kinematics
    protected double posi_px, posi_py, posi_pz, posi_p, posi_e, posi_theta, posi_phi; // positron kinematics
    protected double neut_px, neut_py, neut_pz, neut_p, neut_e, neut_theta, neut_phi; // neutron kinematics

    protected double elec_vx, elec_vy, elec_vz; // electron vertex
    protected double posi_vx, posi_vy, posi_vz; // positron vertex
    protected double neut_vx, neut_vy, neut_vz; // neutron vertex

    // Can and will add more variables later of course but want to keep it simple for now that I am just looking at the
    // event selection :)

    public static boolean channel_test(nTCSParticles variables) {
        if (variables.helicity == 0 && variables.runnum != 11) {
            return false;
        }
        // DVCS version includes some exclusivity conditions, just keeping it simple for now 
        return true;
    }

    public static int getIndex(HipoDataBank rec_Bank, int input_pid, int input_index) {
        int index = -1;
        for (int particle_Index = 0; particle_Index < rec_Bank.rows(); particle_Index++) {
            int pid = rec_Bank.getInt("pid", particle_Index);
            if (pid == input_pid) {
                index++;
            }
            if (index == input_index) {
                return particle_Index;
            }
        }
        return -1;
    }

    public nTCSParticles(DataEvent event, PhysicsEvent recEvent, double Eb) {
        // Alterring to remove PID numbers as input as just for nTCS

        kinematic_variables kinematic_variables = new kinematic_variables();

        // load banks
        HipoDataBank eventBank = (HipoDataBank) event.getBank("REC::Event");
        HipoDataBank configBank = (HipoDataBank) event.getBank("RUN::config");
        HipoDataBank rec_Bank = (HipoDataBank) event.getBank("REC::Particle");
        HipoDataBank cal_Bank = (HipoDataBank) event.getBank("REC::Calorimeter");
        HipoDataBank traj_Bank = (HipoDataBank) event.getBank("REC::Traj");

        helicity = eventBank.getByte("helicity", 0);
        runnum = configBank.getInt("run", 0); // used for beam energy and polarization

        num_electrons = recEvent.countByPid(11); // returns number of electrons
        num_positrons = recEvent.countByPid(-11); // returns number of positrons
        num_piplus = recEvent.countByPid(211);
        num_piminus = recEvent.countByPid(-211);
        num_kplus = recEvent.countByPid(321);
        num_kminus = recEvent.countByPid(-321);
        num_protons = recEvent.countByPid(2212);
        num_antiprotons = recEvent.countByPid(-2212);
        num_neutrons = recEvent.countByPid(2112);
        num_antineutrons = recEvent.countByPid(-2112);
        num_particles = num_electrons + num_piplus + num_piminus + num_kplus + num_kminus + num_protons;
        num_pos = num_positrons + num_piplus + num_kplus + num_protons;
        num_neg = num_electrons + num_piminus + num_kminus + num_antiprotons;
        num_neutrals = recEvent.countByPid(22) + num_neutrons + num_antineutrons;

        generic_tests generic_tests = new generic_tests();
        fiducial_cuts fiducial_cuts = new fiducial_cuts();

        int elec_rec_index = getIndex(rec_Bank, 11, 0);
        int posi_rec_index = getIndex(rec_Bank, -11, 0);
        int neut_rec_index = getIndex(rec_Bank, 2112, 0);

        elec_chi2 = rec_Bank.getFloat("chi2pid", elec_rec_index);
        posi_chi2 = rec_Bank.getFloat("chi2pid", posi_rec_index);
        neut_chi2 = rec_Bank.getFloat("chi2pid", neut_rec_index);

        // Fiducial cuts & fiducial status, will need to add this in later

        // electron detector
        if (generic_tests.forward_tagger_cut(elec_rec_index, rec_Bank)) {
            elec_detector = 0; // Forward Tagger
        } else if (generic_tests.forward_detector_cut(elec_rec_index, rec_Bank)) {
            elec_detector = 1; // Forward Detector
        } else if (generic_tests.central_detector_cut(elec_rec_index, rec_Bank)) {
            elec_detector = 2; // Central Detector
        }
        
        // positron detector
        if (generic_tests.forward_tagger_cut(posi_rec_index, rec_Bank)) {
            posi_detector = 0; // Forward Tagger
        } else if (generic_tests.forward_detector_cut(posi_rec_index, rec_Bank)) {
            posi_detector = 1; // Forward Detector
        } else if (generic_tests.central_detector_cut(posi_rec_index, rec_Bank)) {
            posi_detector = 2; // Central Detector
        }

        // neutron detector
        if (generic_tests.forward_tagger_cut(neut_rec_index, rec_Bank)) {
            neut_detector = 0; // Forward Tagger
        } else if (generic_tests.forward_detector_cut(neut_rec_index, rec_Bank)) {
            neut_detector = 1; // Forward Detector
        } else if (generic_tests.central_detector_cut(neut_rec_index, rec_Bank)) {
            neut_detector = 2; // Central Detector
        }

        // Set up Lorentz vectors
        // target
        LorentzVector target_lv = new LorentzVector();
        target_lv.setPxPyPzM(0,0,0,kinematic_variables.particle_mass(2112));
        // beam electron (not sure how relevant or necessary here but for completeness)
        LorentzVector beam_lv = new LorentzVector();
        beam_lv.setPxPyPzM(0, 0, Math.pow(Eb * Eb - kinematic_variables.particle_mass(11) * kinematic_variables.particle_mass(11), 0.5),
                           kinematic_variables.particle_mass(11));
        // electron
        String electron_string = "[11,0]"; // using found index value for all in case e+ listed first
        Particle scattered_electron = recEvent.getParticle(electron_string);
        LorentzVector elec_lv = new LorentzVector();
        elec_lv.setPxPyPzM(scattered_electron.px(), scattered_electron.py(),
                           scattered_electron.pz(), kinematic_variables.particle_mass(11));
        // positron
        String positron_string = "[-11,0]"; // using found index value for all in case e+ listed first
        Particle scattered_positron = recEvent.getParticle(positron_string);
        LorentzVector posi_lv = new LorentzVector();
        //System.out.println(Double.toString(scattered_positron.pz()));
        posi_lv.setPxPyPzM(scattered_positron.px(), scattered_positron.py(),
                           scattered_positron.pz(), kinematic_variables.particle_mass(-11));
        // neutron
        String neutron_string = "[2112,0]"; // using found index value for all in case e+ listed first
        Particle scattered_neutron = recEvent.getParticle(neutron_string);
        LorentzVector neut_lv = new LorentzVector();
        neut_lv.setPxPyPzM(scattered_neutron.px(), scattered_neutron.py(),
                           scattered_neutron.pz(), kinematic_variables.particle_mass(2112));

        // positions of electron, positron, neutron
        elec_vx = scattered_electron.vx();
        posi_vx = scattered_positron.vx();
        neut_vx = scattered_neutron.vx();
        elec_vy = scattered_electron.vy();
        posi_vy = scattered_positron.vy();
        neut_vy = scattered_neutron.vy();
        elec_vz = scattered_electron.vz();
        posi_vz = scattered_positron.vz();
        neut_vz = scattered_neutron.vz();

        // Initialize momentum corrections at some point

        // kinematics of electron
        elec_px    = elec_lv.px();
        elec_py    = elec_lv.py();
        elec_px    = elec_lv.pz();
        elec_p     = elec_lv.p();
        elec_e     = elec_lv.e();
        elec_theta = scattered_electron.theta();
        elec_phi   = scattered_electron.phi();
        if (elec_phi < 0) {
            elec_phi = 2 * Math.PI + elec_phi;
        }
        // kinematics of positron
        posi_px    = posi_lv.px();
        posi_py    = posi_lv.py();
        posi_pz    = posi_lv.pz();
        posi_p     = posi_lv.p();
        posi_e     = posi_lv.e();
        posi_theta = scattered_positron.theta();
        posi_phi   = scattered_positron.phi();
        if (posi_phi < 0) {
            posi_phi = 2 * Math.PI + posi_phi;
        }
        // kinematics of neutron
        neut_px    = neut_lv.px();
        neut_py    = neut_lv.py();
        neut_pz    = neut_lv.pz();
        neut_p     = neut_lv.p();
        neut_e     = neut_lv.e();
        neut_theta = scattered_neutron.theta();
        neut_phi   = scattered_neutron.phi();
        if (neut_phi < 0) {
            neut_phi = 2 * Math.PI + neut_phi;
        }
    }

    public int get_helicity() { // -1, 0, or 1. 0 equals unassigned by EventBuilder
        if (runnum >= 4326 && runnum <= 5666) {
            return -1 * helicity;
        }  else if (runnum >= 6616 && runnum <= 6783) {
            return -1 * helicity;
        } else if (runnum >= 6120 && runnum <= 6604) {
            return -1 * helicity;
        } else if (runnum >= 11093 && runnum <= 11283) {
            return helicity;
        } else if (runnum >= 11284 && runnum < 11300) {
            return -1 * helicity;
        } else if (runnum >= 11323 && runnum <= 11571) {
            return helicity;
        }
//        System.out.println("runnum not found, assigning helicity flip");
        return helicity;
    } // returns helicity

    public int get_runnum() {
        return runnum;
    } // returns run number

    public int get_elec_detector() {
        return elec_detector;
    } // returns integer value representing the detector of electron track

    public int get_posi_detector() {
        return posi_detector;
    } // returns integer value representing the detector of positron track

    public int get_neut_detector() {
        return neut_detector;
    } // returns integer value representing the detector of neutron track

    public int get_num_pos() {
        return num_pos;
    } // returns number of positively charged particles
    
    public int get_num_neg() {
        return num_neg;
    } // returns number of negatively charged particles
    
    public int get_num_neutrals() {
        return num_neutrals;
    } // returns number of neutrally charged particles

    public int get_num_electrons() {
        return num_electrons;
    } // returns number of electrons

    public int get_num_piplus() {
        return num_piplus;
    } // returns number of piplus

    public int get_num_piminus() {
        return num_piminus;
    } // returns number of piminus

    public int get_num_kplus() {
        return num_kplus;
    }// returns number of kplus

    public int get_num_kminus() {
        return num_kminus;
    } // returns number of kminus

    public int get_num_protons() {
        return num_protons;
    } // returns number of protons

    public int get_num_neutrons() {
        return num_neutrons;
    } // returns number of neutrons

    public int get_num_positrons() {
        return num_positrons;
    } // returns number of positrons

    public double get_elec_chi2pid() {
        return elec_chi2;
    } // returns chi2 value of electron pid

    public double get_posi_chi2pid() {
        return posi_chi2;
    } // returns chi2 value of positron pid

    public double get_neut_chi2pid() {
        return neut_chi2;
    } // returns chi2 value of neutron pid

    public double get_elec_px(){
        return elec_px;
    } // returns electron px

    public double get_elec_py() {
        return elec_py;
    } // returns electron py

    public double get_elec_pz(){
        return elec_pz;
    } // returns electron pz

    public double get_elec_p() {
        return elec_p;
    } // returns electron p

    public double get_elec_e() {
        return elec_e;
    } // returns electron e

    public double get_elec_theta() {
        return elec_theta;
    } // returns electron theta

    public double get_elec_phi() {
        return elec_phi;
    } // returns electron phi

    public double get_posi_px(){
        return posi_px;
    } // returns positron px

    public double get_posi_py() {
        return posi_py;
    } // returns positron py

    public double get_posi_pz(){
        return posi_pz;
    } // returns positron pz

    public double get_posi_p() {
        return posi_p;
    } // returns positron p

    public double get_posi_e() {
        return posi_e;
    } // returns positron e

    public double get_posi_theta() {
        return posi_theta;
    } // returns positron theta

    public double get_posi_phi() {
        return posi_phi;
    } // returns positron phi

    public double get_neut_px(){
        return neut_px;
    } // returns neutron px

    public double get_neut_py() {
        return neut_py;
    } // returns neutron py

    public double get_neut_pz(){
        return neut_pz;
    } // returns neutron pz

    public double get_neut_p() {
        return neut_p;
    } // returns neutron p

    public double get_neut_e() {
        return neut_e;
    } // returns neutron e

    public double get_neut_theta() {
        return neut_theta;
    } // returns neutron theta

    public double get_neut_phi() {
        return neut_phi;
    } // returns neutron phi

    public double get_elec_vx() {
        return elec_vx;
    } // returns electron vx

    public double get_elec_vy() {
        return elec_vy;
    } // returns electron vy

    public double get_elec_vz() {
        return elec_vz;
    } // returns electron vz

    public double get_posi_vx() {
        return posi_vx;
    } // returns positron vx

    public double get_posi_vy() {
        return posi_vy;
    } // returns positron vy

    public double get_posi_vz() {
        return posi_vz;
    } // returns positron vz

    public double get_neut_vx() {
        return neut_vx;
    } // returns neutron vx

    public double get_neut_vy() {
        return neut_vy;
    } // returns neutron vy

    public double get_neut_vz() {
        return neut_vz;
    } // returns neutron vz

}