package org.simbrain.special_projects.plasticity_proj;

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
        SymmetricSTDPRule duplicateSynapse = new  SymmetricSTDPRule();
        duplicateSynapse.setTau_minus(this.getTau_minus());
        duplicateSynapse.setTau_plus(this.getTau_plus());
        duplicateSynapse.setW_minus(this.getW_minus());
        duplicateSynapse.setW_plus(this.getW_plus());
        duplicateSynapse.setLearningRate(this.getLearningRate());
        duplicateSynapse.setHebbian(hebbian);
        return duplicateSynapse;
    }
    
    private boolean hebbian = true;

    private double delta_w = 0;
    
    @Override
    public void update(Synapse synapse) {
        final double str = synapse.getStrength();
        if (synapse.getSource().isSpike() || synapse.getTarget().isSpike()) {
            double delta_t = ((((SpikingNeuronUpdateRule) synapse
                    .getSource().getUpdateRule()).getLastSpikeTime())
                    - ((SpikingNeuronUpdateRule) synapse
                            .getTarget().getUpdateRule()).getLastSpikeTime())
                            * (hebbian ? 1 : -1); // Reverse time window for
            delta_t += ThreadLocalRandom.current().nextGaussian()
                    * (delta_t / 10);
            double ptfr = ((IPIFRule) synapse.getTarget()
                    .getUpdateRule()).getEstFR();
            double psfr = ((IPIFRule) synapse.getSource()
                    .getUpdateRule()).getEstFR();
            double pfr = ptfr + psfr;
            double wm;
            double tm;
            if (pfr < 100) {
                wm = (((-W_plus - W_minus) / 100) * pfr) + W_minus;
                tm = (((tau_plus - tau_minus) / 100) * pfr) + tau_minus;
            } else {
                wm = -W_plus;
                tm = tau_plus;
            }

            if (delta_t < 0) {
                delta_w = W_plus * Math.exp(delta_t / tau_plus)
                        * learningRate * Math.exp(-Math.abs(str) / 200);
            } else if (delta_t > 0) {
                delta_w = -wm * Math.exp(-delta_t / tm)
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
}
