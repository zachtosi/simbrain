package org.simbrain.special_projects.plasticity_proj;

import java.io.PrintWriter;
import java.lang.reflect.Field;

import org.simbrain.network.core.SpikingNeuronUpdateRule;
import org.simbrain.network.core.Synapse;
import org.simbrain.network.synapse_update_rules.STDPRule;
import org.simbrain.util.SimbrainConstants.Polarity;

public class MexicanHatSTDPRule extends STDPRule {
	public static final double PI_4th_RT = Math.pow(Math.PI, 0.25);
	
	private double delta_w;
	
	public double sigma = 30;
	
	private double a = 25;
	
	public MexicanHatSTDPRule deepCopy() {
		MexicanHatSTDPRule mhcpy = new MexicanHatSTDPRule();
		mhcpy.learningRate = this.learningRate;
		mhcpy.sigma = sigma;
		mhcpy.a = a;
		mhcpy.W_plus = W_plus;
		mhcpy.W_minus = W_minus;
		return mhcpy;
	}
	
	 public void update(Synapse synapse) {
	        final double str = synapse.getStrength();
	        if (synapse.spkArrived() || synapse.getTarget().isSpike()) {
	            double delta_t = (((SpikingNeuronUpdateRule) synapse
	                    .getSource().getUpdateRule()).getLastSpikeTime()
	                    - ((SpikingNeuronUpdateRule) synapse
	                            .getTarget().getUpdateRule()).getLastSpikeTime());
	            delta_w = this.learningRate * mexican_hat(delta_t, sigma, a);
	            if(delta_w > 0) {
	            	delta_w *= W_plus;//Math.exp(-Math.abs(str) / 15.0) * W_plus;
	            } else {
	            	delta_w *= W_minus;
	            }
	        
	        }
	        double was = Math.abs(str);
	        // Subtracts deltaW if inhibitory adds otherwise
	        int sign = (int) Math.signum(str);
	        if (sign != 0) {
	            synapse.changeStrength((sign * delta_w));
	        } else {
	            if (synapse.getSource().getPolarity() == Polarity.EXCITATORY) {
	                synapse.changeStrength(delta_w);
	            } else {
	                synapse.changeStrength(-delta_w);
	            }
	        }
//	        double is = Math.abs(synapse.getStrength());
//	        if (was < 1 && is >= 1) {
//	            synapse.setDelay(synapse.getDelay() - 2);
//	        }
//	        if (was > 1 && is <= 1) {
//	            synapse.setDelay(synapse.getDelay() + 2);
//	        }
//	        if (was < 5 && is >= 5) {
//	            synapse.setDelay(synapse.getDelay() - 2);
//	        }
//	        if (was > 5 && is <= 5) {
//	            synapse.setDelay(synapse.getDelay() + 2);
//	        }
//	        if (was < 10 && is >= 10) {
//	            synapse.setDelay(synapse.getDelay() - 2);
//	        }
//	        if (was > 10 && is <= 10) {
//	            synapse.setDelay(synapse.getDelay() + 2);
//	        }
//	        if (was < 15 && is >= 15) {
//	            synapse.setDelay(synapse.getDelay() - 2);
//	        }
//	        if (was > 15 && is <= 15) {
//	            synapse.setDelay(synapse.getDelay() + 2);
//	        }
	    }
	
	//TODO: Generalize this to SimbrainMath
	public static double mexican_hat(double x, double sigma, double a) {
		return a * (2.0 / (Math.sqrt(3*sigma) * PI_4th_RT)) * (1 - ((x*x)/(sigma*sigma)))
				* Math.exp(-(x*x) / (2 * sigma * sigma));
	}

    public void reportAllValsToFile(PrintWriter pw)
            throws IllegalAccessException {
        pw.println("\t\t**[" + this.toString() + "]**");
        for (Field f : MexicanHatSTDPRule.class.getDeclaredFields()) {
            pw.println("\t" + f.getName() + ": " + f.get(this));
        }
    }
    
    public String toString() {
        return "Mexican Hat STDP";
    }
	
}
