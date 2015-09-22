package org.simbrain.special_projects.plasticity_proj;

import java.io.PrintWriter;
import java.lang.reflect.Field;
import java.util.concurrent.ThreadLocalRandom;

import org.simbrain.network.core.SpikingNeuronUpdateRule;
import org.simbrain.network.core.Synapse;
import org.simbrain.network.synapse_update_rules.STDPRule;
import org.simbrain.util.SimbrainConstants.Polarity;

public class SymmetricSTDPRule extends STDPRule{
    /** Time constant for LTP. */
    protected double tau_plus = TAU_PLUS_DEFAULT;

    /** Time constant for LTD. */
    protected double tau_minus = TAU_MINUS_DEFAULT;

    /**
     * Learning rate for LTP case. Controls magnitude of LTP changes.
     */
    protected double W_plus = W_PLUS_DEFAULT;

    /**
     * Learning rate for LTP case. Controls magnitude of LTD changes.
     */
    protected double W_minus = W_MINUS_DEFAULT;

    /** General learning rate. */
    protected double learningRate = LEARNING_RATE_DEFAULT;

    @Override
    public void init(Synapse synapse) {
    }

    @Override
    public String getDescription() {
        return "STDP";
    }

    @Override
    public SymmetricSTDPRule deepCopy() {
	ThreadLocalRandom tlr = ThreadLocalRandom.current();
        SymmetricSTDPRule duplicateSynapse = new  SymmetricSTDPRule();
        duplicateSynapse.setTau_minus(this.getTau_minus() + (tlr.nextGaussian() * (getTau_minus()/10)));
        duplicateSynapse.setTau_plus(this.getTau_plus() + (tlr.nextGaussian() * (getTau_plus()/10))) ;
        duplicateSynapse.setW_minus(this.getW_minus());
        duplicateSynapse.setW_plus(this.getW_plus());
        duplicateSynapse.setLearningRate(this.getLearningRate());
        duplicateSynapse.setHebbian(hebbian);
        duplicateSynapse.setSymmetric(isSymmetric());
        duplicateSynapse.setSymTau(symTau);
        duplicateSynapse.setWtSoftening(wtSoftening);
        return duplicateSynapse;
    }
    
    private boolean hebbian = true;

    private double delta_w = 0;
    
    private boolean symmetric = true;

    private double symTau = 20;

    public double getSymTau() {
	return symTau;
    }

    public void setSymTau(final double symTau) {
	this.symTau = symTau;
    }
    
    @Override
    public void update(Synapse synapse) {
        final double str = synapse.getStrength();
        if (synapse.spkArrived() || synapse.getTarget().isSpike()) {
            double delta_t = (synapse.getLastSpkArrival()
                    - ((SpikingNeuronUpdateRule) synapse
                            .getTarget().getUpdateRule()).getLastSpikeTime())
                            * (hebbian ? 1 : -1); // Reverse time window for
            double wm;
            double tm;
            double eVal;
            if (symmetric) {
                double ptfr = ((IPIFRule) synapse.getTarget()
                        .getUpdateRule()).getEstFR();
                double psfr = ((IPIFRule) synapse.getSource()
                        .getUpdateRule()).getEstFR();
                if (ptfr > 50) {
                    ptfr = 50;
                }
                if (psfr > 50) {
                    psfr = 50;
                }
                double pfr = ptfr + psfr;
                eVal = Math.exp((pfr - 100) / symTau);

                if (pfr < 100) {
                    //wm = (((-W_plus - W_minus) / 100) * pfr) + W_minus;
                    //tm = (((tau_plus - tau_minus) / 100) * pfr) + tau_minus;
                    wm = ((-W_plus - W_minus) * eVal)
                            + W_minus;
                    tm = ((tau_plus - tau_minus) * eVal) + tau_minus;
                } else {
                    wm = -W_plus;
                    tm = tau_plus;
                }
            } else {
                wm = W_minus;
                tm = tau_minus;
                eVal = 0;
            }
            if (delta_t < 0) {
                delta_w = learningRate * (W_plus * Math.exp(delta_t / tau_plus) - (eVal * 0.5))
                        * Math.exp(-Math.abs(str) / wtSoftening);
            } else if (delta_t > 0) {
                delta_w = (-wm * Math.exp(-delta_t / tm) - (eVal * 0.5))
                        * learningRate; 
            }
        }
        // Subtracts deltaW if inhibitory adds otherwise
        int sign = (int) Math.signum(str);
        if (sign != 0) {
            synapse.setStrength(str + (sign * delta_w));
        } else {
            if (synapse.getSource().getPolarity() == Polarity.EXCITATORY) {
                synapse.setStrength(str + delta_w);
            } else {
                synapse.setStrength(str - delta_w);
            }
        }
    }

    private float wtSoftening = 60;
    
    /**
     * @return the tau_plus
     */
    public double getTau_plus() {
        return tau_plus;
    }

    /**
     * @param tauPlus
     *            the tau_plus to set
     */
    public void setTau_plus(double tauPlus) {
        tau_plus = tauPlus;
    }

    /**
     * @return the tau_minus
     */
    public double getTau_minus() {
        return tau_minus;
    }

    /**
     * @param tauMinus
     *            the tau_minus to set
     */
    public void setTau_minus(double tauMinus) {
        tau_minus = tauMinus;
    }

    /**
     * @return the w_plus
     */
    public double getW_plus() {
        return W_plus;
    }

    /**
     * @param wPlus
     *            the w_plus to set
     */
    public void setW_plus(double wPlus) {
        W_plus = wPlus;
    }

    /**
     * @return the w_minus
     */
    public double getW_minus() {
        return W_minus;
    }

    /**
     * @param wMinus
     *            the w_minus to set
     */
    public void setW_minus(double wMinus) {
        W_minus = wMinus;
    }

    /**
     * @return the learningRate
     */
    public double getLearningRate() {
        return learningRate;
    }

    /**
     * @param learningRate
     *            the learningRate to set
     */
    public void setLearningRate(double learningRate) {
        this.learningRate = learningRate;
    }

    public boolean isHebbian() {
        return hebbian;
    }

    public void setHebbian(boolean hebbian) {
        this.hebbian = hebbian;
    }

    public boolean isSymmetric() {
        return symmetric;
    }

    public void setSymmetric(boolean symmetric) {
        this.symmetric = symmetric;
    }

    public float getWtSoftening() {
        return wtSoftening;
    }

    public void setWtSoftening(float wtSoftening) {
        this.wtSoftening = wtSoftening;
    }
    
    public void reportAllValsToFile(PrintWriter pw)
            throws IllegalAccessException {
        pw.println("\t\t**[" + this.toString() + "]**");
        for (Field f : SymmetricSTDPRule.class.getDeclaredFields()) {
            pw.println("\t" + f.getName() + ": " + f.get(this));
        }
    }
    
    public String toString() {
        return (symmetric ? "Symmetric" : "Non-Symmetric") + "_"
                + (hebbian ? "Hebbian" : "Anti-Hebbian") + "_" + "STDP";
    }
}
