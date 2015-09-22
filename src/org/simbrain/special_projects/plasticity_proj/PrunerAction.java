package org.simbrain.special_projects.plasticity_proj;

import java.io.PrintWriter;
import java.lang.reflect.Field;
import java.util.concurrent.ThreadLocalRandom;

import org.simbrain.network.core.NetworkUpdateAction;
import org.simbrain.network.core.Synapse;
import org.simbrain.network.groups.SynapseGroup;

public class PrunerAction implements NetworkUpdateAction {

    private final SynapseGroup synGrp;
    
    private double pruneInterval = 10000; // 10s
    
    private double pruneThresholdEx = 0.05;
    
    private double pruneThresholdIn = 0.05;
    
    private double minVal = 1E-5;
    
    public PrunerAction(SynapseGroup synGrp) {
        this.synGrp = synGrp;
    }
    
    @Override
    public String getDescription() {
        return synGrp.getLabel() + " synapse pruner";
    }

    @Override
    public String getLongDescription() {
        return "Prunes synapses marked for deletion from a synapse group";
    }

    @Override
    public void invoke() {
        final double time = synGrp.getParentNetwork().getTime();
        final double timeStep = synGrp.getParentNetwork().getTimeStep();
        final int srcSize = synGrp.getSourceNeuronGroup().size();
        final int tarSize = synGrp.getTargetNeuronGroup().size();
        ThreadLocalRandom tlr = ThreadLocalRandom.current();
        if (((int) (time / timeStep)) % ((int)(pruneInterval / timeStep)) == 0) {
            System.out.println("\t" + synGrp.size());
            double max = Double.MIN_VALUE;
            for (Synapse s : synGrp.getExcitatorySynapses()) {
                double str = s.getStrength();
                if (str > max) {
                    max = str;
                }
            }
            for (Synapse s : synGrp.getExcitatorySynapses()) {
                if (s.getStrength() < minVal) {
                    synGrp.removeSynapse(s);
                    continue;
                }
                if (Math.abs(s.getStrength()) < max * pruneThresholdEx) {
                    double removeP = (((double) s.getSource().getFanOut().size()
                            / tarSize) + ((double) s.getTarget().getFanIn()
                                    .size() / srcSize)) / 2.0;
                    if (tlr.nextDouble() < removeP) {
                        synGrp.removeSynapse(s);
                    }
                }
            }
            max = Double.MAX_VALUE;
            for (Synapse s : synGrp.getInhibitorySynapses()) {
                double str = s.getStrength();
                if (str < max) {
                    max = str;
                }
            }
            for (Synapse s : synGrp.getInhibitorySynapses()) {
                if (s.getStrength() > -minVal) {
                    synGrp.removeSynapse(s);
                    continue;
                }
                if (s.getStrength() > max * pruneThresholdIn) {
                    double removeP = (((double) s.getSource().getFanOut().size()
                            / tarSize) + ((double) s.getTarget().getFanIn()
                                    .size() / srcSize)) / 2.0;
                    if (tlr.nextDouble() < removeP) {
                        synGrp.removeSynapse(s);
                    }
                }
            }
            System.out.print("\t" + synGrp.size() + " ");
            System.out.println((double) synGrp.size() / (synGrp
                    .getSourceNeuronGroup().size() * (synGrp
                    .getTargetNeuronGroup().size()-1)));
        }
        
    }
    
    public void reportAllValsToFile(PrintWriter pw)
            throws IllegalAccessException {
        pw.println("\t\t**[" + this.toString() + "]**");
        for (Field f : PrunerAction.class.getDeclaredFields()) {
            pw.println("\t" + f.getName() + ": " + f.get(this));
        }
    }
    
    public String toString() {
        return "Pruner Action";
    }
}
