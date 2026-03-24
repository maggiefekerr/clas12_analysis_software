/**
 *
 * @author Maggie F. E. Kerr
 */
package extended_kinematic_fitters;

import org.jlab.clas.physics.GenericKinematicFitter;
import org.jlab.clas.physics.Particle;
import org.jlab.clas.physics.PhysicsEvent;
import org.jlab.io.base.DataEvent;
import org.jlab.io.hipo.HipoDataBank;

import org.jlab.clas.physics.*;

public class xtcs_fitter extends GenericKinematicFitter {

    protected final Double mybeam;

    public xtcs_fitter(double beam) {
        super(beam);
        mybeam = beam;
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////
    @Override
    public PhysicsEvent getPhysicsEvent(DataEvent event) {
        
        generic_tests generic_tests = new generic_tests();
        if (generic_tests.banks_test(event)) {
            PhysicsEvent physEvent = new PhysicsEvent();
            // load the hipo banks
            // assumption is we are using trains which would require all of these banks to exist
            HipoDataBank rec_Bank = (HipoDataBank) event.getBank("REC::Particle");
            HipoDataBank cal_Bank = (HipoDataBank) event.getBank("REC::Calorimeter");
            HipoDataBank cc_Bank = (HipoDataBank) event.getBank("REC::Cherenkov");
            HipoDataBank traj_Bank = (HipoDataBank) event.getBank("REC::Traj");
            HipoDataBank track_Bank = (HipoDataBank) event.getBank("REC::Track");
            HipoDataBank run_Bank = (HipoDataBank) event.getBank("RUN::config");
            HipoDataBank ft_Bank = null;
            if (event.hasBank("REC::ForwardTagger")) {
                ft_Bank = (HipoDataBank) event.getBank("REC::ForwardTagger");
            }
            double vz_e = -999;

            LorentzVector lv_e = new LorentzVector();
            if (rec_Bank.getInt("pid", 0) == 11) {
                // trigger particle was an electron
                // highest momentum electron listed first (used for DIS calculations)
                float px = rec_Bank.getFloat("px", 0);
                float py = rec_Bank.getFloat("py", 0);
                float pz = rec_Bank.getFloat("pz", 0);
                double p = Math.sqrt(px * px + py * py + pz * pz);
                lv_e.setPxPyPzM(px, py, pz, 0.0005109989461);
                vz_e = rec_Bank.getFloat("vz", 0);

            } else {
                return physEvent;
            } // trigger particle was not an electron

            for (int particle_Index = 0; particle_Index < rec_Bank.rows(); particle_Index++) {
                int pid = rec_Bank.getInt("pid", particle_Index);
                float px = rec_Bank.getFloat("px", particle_Index);
                float py = rec_Bank.getFloat("py", particle_Index);
                float pz = rec_Bank.getFloat("pz", particle_Index);
                float vx = rec_Bank.getFloat("vx", particle_Index);
                float vy = rec_Bank.getFloat("vy", particle_Index);
                float vz = rec_Bank.getFloat("vz", particle_Index);
                float chi2pid = rec_Bank.getFloat("chi2pid", particle_Index);
                double p = Math.sqrt(px * px + py * py + pz * pz);

                int sector = generic_tests.sector(particle_Index, track_Bank); // 0 FT/CD, 1-6 FD

                int runnum = run_Bank.getInt("run", 0);

                boolean inbending = false;
                boolean outbending = false;
                if (run_Bank.getFloat("torus", 0) == 1) {
                    outbending = true;
                } else {
                    inbending = true;
                }

                if (pid == 11) {
                    float[] momentum = {px, py, pz};
                    px = momentum[0];
                    py = momentum[1];
                    pz = momentum[2];
                    Particle elec = new Particle(pid, px, py, pz, vx, vy, vz_e);
                    physEvent.addParticle(elec);
                }

                if (pid == -11) {
                    float[] momentum = {px, py, pz};
                    px = momentum[0];
                    py = momentum[1];
                    pz = momentum[2];
                    Particle posi = new Particle(pid, px, py, pz, vx, vy, vz);
                    physEvent.addParticle(posi);
                }

                if (pid == 2212 && (p > 0.0001)) {
                    float[] momentum = {px, py, pz};
                    px = momentum[0];
                    py = momentum[1];
                    pz = momentum[2];
                    Particle prot = new Particle(pid, px, py, pz, vx, vy, vz);
                    physEvent.addParticle(prot);
                }

                if (pid == 2112 && (p > 0.0001)) {
                    float[] momentum = {px, py, pz};
                    px = momentum[0];
                    py = momentum[1];
                    pz = momentum[2];
                    Particle neut = new Particle(pid, px, py, pz, vx, vy, vz);
                    physEvent.addParticle(neut);
                }
            }
            return physEvent;
        }
        return new PhysicsEvent(this.mybeam);
    }
}