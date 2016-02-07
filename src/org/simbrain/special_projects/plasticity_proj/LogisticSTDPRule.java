package org.simbrain.special_projects.plasticity_proj;

import java.io.PrintWriter;
import java.lang.reflect.Field;

import org.simbrain.network.core.SpikingNeuronUpdateRule;
import org.simbrain.network.core.Synapse;
import org.simbrain.network.synapse_update_rules.STDPRule;
import org.simbrain.util.SimbrainConstants.Polarity;

public class LogisticSTDPRule extends STDPRule {
	public static final double PI_4th_RT = Math.pow(Math.PI, 0.25);
	
	private double delta_w;
	
	public LogisticSTDPRule deepCopy() {
		LogisticSTDPRule mhcpy = new LogisticSTDPRule();
		mhcpy.learningRate = this.learningRate;
		mhcpy.W_plus = W_plus;
		mhcpy.W_minus = W_minus;
		return mhcpy;
	}
	
	 public void update(Synapse synapse) {
	        final double str = synapse.getStrength();
	        if (synapse.spkArrived() || synapse.getTarget().isSpike()) {
	            double delta_t = (((SpikingNeuronUpdateRule) synapse
                        .getTarget().getUpdateRule()).getLastSpikeTime()
	                    - ((SpikingNeuronUpdateRule) synapse
	                            .getTarget().getUpdateRule()).getLastSpikeTime());
	            
	            delta_w = learningRate * double2ndDerLogit(delta_t);
//	            if(delta_w > 0) {
//	            	delta_w *= Math.exp(-Math.abs(str) / 20.0) * W_plus;
//	            } else {
//	            	delta_w *= W_minus;
//	            }
	        
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
	
	//TODO: Generalize this to SimbrainMath
	public static double double2ndDerLogit(double x) {
		double dnomFull = 1.0/(1.0 + (0.25*x*x));
		double dnomNum = Math.exp(x/4) + 1;
		dnomNum = dnomNum*dnomNum*dnomNum;
		double numer = -4 * Math.exp(x/4) * (Math.exp(x/4) - 1);
		double ret = (numer/dnomNum)/dnomFull;
		if (Double.isNaN(ret) || Double.isInfinite(ret)) {
			ret = 0;
			//System.out.println(numer + " " + dnomNum + " " + dnomFull);
		}
		return ret;
	}

    public void reportAllValsToFile(PrintWriter pw)
            throws IllegalAccessException {
        pw.println("\t\t**[" + this.toString() + "]**");
        for (Field f : LogisticSTDPRule.class.getDeclaredFields()) {
            pw.println("\t" + f.getName() + ": " + f.get(this));
        }
    }
    
    public String toString() {
        return "Logistic STDP";
    }
}
