package org.simbrain.special_projects.plasticity_proj;

import java.io.PrintWriter;
import java.lang.reflect.Field;
import java.util.concurrent.ThreadLocalRandom;

import org.simbrain.network.core.NetworkUpdateAction;
import org.simbrain.network.core.Neuron;
import org.simbrain.network.core.Synapse;
import org.simbrain.network.groups.NeuronGroup;
import org.simbrain.network.groups.SynapseGroup;
import org.simbrain.util.SimbrainConstants.Polarity;

public class PrunerAction implements NetworkUpdateAction {

    private final SynapseGroup synGrp;
    
    private double pruneInterval = 5000; // 5s
    
    private double pruneThresholdEx = 0.02;
    
    private double pruneThresholdIn = 0.02;
    
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
        NeuronGroup ng = synGrp.getSourceNeuronGroup();
        int ex = 0;
        int in = 0;
        for (Neuron n : ng.getNeuronList()) {
        	if (n.getPolarity() == Polarity.EXCITATORY) {
        		ex++;
        	}
        	if (n.getPolarity() == Polarity.INHIBITORY) {
        		in++;
        	}
        }
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
                	int exin = 0;
                	for (int k = 0; k < s.getTarget().getFanIn().size(); k++) {
                		if (s.getTarget().getFanIn().get(k).getStrength() > 0)
                				//&& s.getTarget().getFanIn()
                				//.get(k).getAuxVal() == 0) 
                			exin++;
                	}
                	double exInRat = (double)exin/ex;
                	double outRat = (double)(s.getSource().getFanOut().size()+5)/tarSize;
                	double removeP = exInRat * outRat +.01; //* (((IPIFRule)s.getTarget().getUpdateRule()).getThreshold()/15) + .01;
                    //double removeP = (0.5 * ((((double) s.getSource().getFanOut().size()
                    //        / tarSize) + ((double) s.getTarget().getFanIn()
                    //                .size() / srcSize)) / 2.0)) + 0.01;
                	//double removeP = (0.5 * ((double) s.getSource().getFanOut().size() / tarSize)) + 0.01;
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
                	int inin = 0;
                	for (int k = 0; k < s.getTarget().getFanIn().size(); k++) {
                		if (s.getTarget().getFanIn().get(k).getStrength() < 0)
                				//&& s.getTarget().getFanIn()
                				//.get(k).getAuxVal() == 0) 
                			inin++;
                	}
                	double inInRat = (double)inin/in;
                	double outRat = (double)(s.getSource().getFanOut().size()+5)/tarSize;
                	double removeP = inInRat * outRat +.01; //* (15/((IPIFRule)s.getTarget().getUpdateRule()).getThreshold()) + .01;
//                	double removeP = (0.5 * ((((double) s.getSource().getFanOut().size()
//                            / tarSize) + ((double) s.getTarget().getFanIn()
//                                    .size() / srcSize)) / 2.0)) + 0.01;
                    //double removeP = (0.5 * ((double) s.getSource().getFanOut().size() / tarSize)) + 0.01;
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
