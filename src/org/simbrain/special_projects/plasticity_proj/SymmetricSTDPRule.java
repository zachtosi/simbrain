package org.simbrain.special_projects.plasticity_proj;

import java.util.concurrent.ThreadLocalRandom;

import org.simbrain.network.core.SpikingNeuronUpdateRule;
import org.simbrain.network.core.Synapse;
import org.simbrain.network.synapse_update_rules.STDPRule;

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
        duplicateSynapse.setDampen(dampen);
        duplicateSynapse.setDampeningFactor(dampeningFactor);
        duplicateSynapse.setNoisy(noisy);
        return duplicateSynapse;
    }
    
//    public enum WindowType {
//        HEBBIAN {
//            @Override
//            public double wtChange(double originalWtChange,
//                    boolean excite) {
//                
//                return 0;
//            }
//        },
//        ANTI_HEBBIAN {
//            @Override
//            public double wtChange(double originalWtChange,
//                    boolean excite) {
//                // TODO Auto-generated method stub
//                return 0;
//            }
//        },
//        SYMMETRIC{
//            @Override
//            public double wtChange(double originalWtChange,
//                    boolean excite) {
//                return Math.abs(originalWtChange);
//            }
//        };
//        public abstract double wtChange(double originalWtChange,
//                boolean excite);
    // }

    private boolean hebbian = true;

    private boolean dampen = false;

    private double dampeningFactor = 0.1;

    private boolean noisy = true;

    private double noiseStd = 0.07;

    public void setNoisy(boolean noisy) {
        this.noisy = noisy;
    }

    public boolean isNoisy() {
        return noisy;
    }

    private double delta_w = 0;
    
    public static double coef = 1.0 / (2*Math.PI*8000);
    
    @Override
    public void update(Synapse synapse) {
        final double str = synapse.getStrength();
        final double delta_t = ((((SpikingNeuronUpdateRule) synapse
                .getSource().getUpdateRule()).getLastSpikeTime())
                - ((SpikingNeuronUpdateRule) synapse
                        .getTarget().getUpdateRule()).getLastSpikeTime())
                        * (hebbian ? 1 : -1); // Reverse time window for
        double pfr = ((IPIFRule) synapse.getTarget()
                .getUpdateRule()).getEstFR();
        double wm;
        double tm;
        if (pfr < 55) {
            wm = (((-W_plus - W_minus) / 55) * pfr) + W_minus;
            tm = (((tau_plus - tau_minus) / 55) * pfr) + tau_minus;
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
        if (noisy) {
            delta_w *= (1 + (ThreadLocalRandom.current()
                    .nextGaussian() * noiseStd));
        }
        if (dampen && synapse.getSource().isSpike()) {
            delta_w -= dampeningFactor;
        }
        // Subtracts deltaW if inhibitory adds otherwise
        int sign = (int) Math.signum(str);
        if (sign != 0) {
            synapse.setStrength(str + (sign * delta_w));
        } else {
            synapse.setStrength(str + Math.abs(delta_w));
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
//        if (this.hebbian != hebbian) {
//            double holder = W_plus;
//            W_plus = W_minus;
//            W_minus = holder;
//        }
        this.hebbian = hebbian;
    }

    public boolean isDampen() {
        return dampen;
    }

    public void setDampen(boolean dampen) {
        this.dampen = dampen;
    }

    public double getDampeningFactor() {
        return dampeningFactor;
    }

    public void setDampeningFactor(double dampeningFactor) {
        this.dampeningFactor = dampeningFactor;
    }
}
