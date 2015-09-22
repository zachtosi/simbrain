package org.simbrain.special_projects.plasticity_proj;

import java.io.PrintWriter;
import java.lang.reflect.Field;

import org.simbrain.network.core.Synapse;
import org.simbrain.network.core.SynapseUpdateRule;
import org.simbrain.util.SimbrainConstants.Polarity;

public class ModifiedBCMRule extends SynapseUpdateRule {

    private double learningRate = 5E-3;

    @Override
    public void init(Synapse synapse) {
        // TODO Auto-generated method stub

    }
    
    
    int dummy = 0;
    @Override
    public void update(Synapse synapse) {
        if (synapse.getParentNetwork().getTime() < 1000) {
            return;
        }
//        double x = ((IPIFRule) synapse.getSource().getUpdateRule())
//                .getApproxFR();
//        double y = ((IPIFRule) synapse.getTarget().getUpdateRule())
//                .getApproxFR();
//        double theta = ((IPIFRule) synapse.getTarget().getUpdateRule())
//                .getMeanApproxFR();
//        double delta_w = learningRate * x * y * (y - theta);
        float x = ((IPIFRule) synapse.getSource().getUpdateRule())
                .getMeanApproxFR();
        float y = ((IPIFRule) synapse.getTarget().getUpdateRule())
                .getMeanApproxFR();
        float theta = ((IPIFRule) synapse.getTarget().getUpdateRule())
                .getMafrLong();
        double delta_w = 0;
        if (y > theta * (0.1f)) {
            delta_w = x * (y - theta);
        } else {
            delta_w = -x * y * (0.9) / 0.1;
        }
        delta_w *= learningRate;
        double str = synapse.getStrength();
        if (delta_w > 0) {
            delta_w *= 2;
            if (str > 0){
                delta_w *= 1.5;
            }
        } else {
            delta_w *= 0.75;
        }
        // if (str < 0 && delta_w < 0) {
        // delta_w *= 2;
        // }
//         if (delta_w > 0) {
//             delta_w *= Math.exp(-Math.abs(str) / 40);
//         }
        //
        
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

    @Override
    public ModifiedBCMRule deepCopy() {
        ModifiedBCMRule cpy = new ModifiedBCMRule();
        cpy.setLearningRate(getLearningRate());
        return cpy;
    }

    @Override
    public String getDescription() {
        return "Modified BCM Rule";
    }

    public double getLearningRate() {
        return learningRate;
    }

    public void setLearningRate(double learningRate) {
        this.learningRate = learningRate;
    }
    
    public void reportAllValsToFile(PrintWriter pw)
            throws IllegalAccessException {
        pw.println("\t\t**[" + this.toString() + "]**");
        for (Field f : ModifiedBCMRule.class.getDeclaredFields()) {
            pw.println("\t" + f.getName() + ": " + f.get(this));
        }
    }
    
    public String toString() {
        return "ModifiedBCMRule";
    }

}
