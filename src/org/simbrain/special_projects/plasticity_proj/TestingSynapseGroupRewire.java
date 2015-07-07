package org.simbrain.special_projects.plasticity_proj;

import java.util.Arrays;

import org.simbrain.network.connections.Sparse;
import org.simbrain.network.core.Network;
import org.simbrain.network.groups.NeuronGroup;
import org.simbrain.network.groups.SynapseGroup;
import org.simbrain.util.SimbrainConstants.Polarity;
import org.simbrain.util.math.ProbDistribution;
import org.simbrain.util.randomizer.PolarizedRandomizer;

public class TestingSynapseGroupRewire {

    
    public static void main(String[] args) {
        
        Network net = new Network(); // Create the network
        // Create a neuron group with 1000 neurons
        NeuronGroup ng = new NeuronGroup(net, 100);
        // Create a sparse connector which specifies a synapse density of 0.1
        // i.e. of all possible synaptic connections create only 1/10 of them
        Sparse sp = new Sparse(0.1, false, false);
        // Set up the randomizers used to determine synaptic strength
        PolarizedRandomizer prEx = new PolarizedRandomizer(Polarity.EXCITATORY,
                ProbDistribution.LOGNORMAL);
        prEx.setParam1(1.5); // Location
        prEx.setParam2(0.75); // Scale
        PolarizedRandomizer prIn = new PolarizedRandomizer(Polarity.INHIBITORY,
                ProbDistribution.LOGNORMAL);
        prIn.setParam1(3.0);
        prIn.setParam2(2.5);
        // Create a group of synapses connecting ng to itself, using the
        // sparse connector defined above where 50% of the synapses are 
        // excitatory, using the randomizers defined above to determine
        // their strengths.
        SynapseGroup synGrp = SynapseGroup.createSynapseGroup(ng, ng, sp, 0.5,
                prEx, prIn);
        
        // Store all the parameters before rewiring that should remain the same
        // afterwards
        double trueSparsity = (double) synGrp.size() / (ng.size()
                * (ng.size() - 1));
        double [] exWts = synGrp.getExcitatoryStrengths();
        double [] inWts = synGrp.getInhibitoryStrengths();
        int [] inDegs = new int[ng.size()];
        int [] outDegs = new int[ng.size()];
        for (int i = 0, n = ng.size(); i < n; i++) {
            inDegs[i] = ng.getNeuronList().get(i).getInDegree();
            outDegs[i] = ng.getNeuronList().get(i).getOutDegree();
        }
        
        System.out.println("Begin Rewiring...:");
        long start = System.nanoTime();
        // Perform rewiring...
        MainSimulation.rewireSynapseGroup(synGrp, 100);
        long end = System.nanoTime();
        System.out.println("COMPLETE");
        System.out.println("Time: " + (end - start) / 1E9 + "s");
        
        // Store all the parameters after rewiring
        double trueSparsityRw = (double) synGrp.size() / (ng.size()
                * (ng.size() - 1));
        double [] exWtsRw = synGrp.getExcitatoryStrengths();
        double [] inWtsRw = synGrp.getInhibitoryStrengths();
        int [] inDegsRw = new int[ng.size()];
        int [] outDegsRw = new int[ng.size()];
        for (int i = 0, n = ng.size(); i < n; i++) {
            inDegsRw[i] = ng.getNeuronList().get(i).getInDegree();
            outDegsRw[i] = ng.getNeuronList().get(i).getOutDegree();
        }
        
        // Sort all arrays making them directly comparable.
        Arrays.sort(exWts);
        Arrays.sort(inWts);
        Arrays.sort(inDegs);
        Arrays.sort(outDegs);
        Arrays.sort(exWtsRw);
        Arrays.sort(inWtsRw);
        Arrays.sort(inDegsRw);
        Arrays.sort(outDegsRw);


        System.out.println("Testing that no synapses were lost...");
        if (trueSparsity == trueSparsityRw) {
            System.out.println("PASSED!");
        } else {
            System.out.print("FAILED! ");
            System.out.println("Original: " + trueSparsity + " " + "Rw: "
            + trueSparsityRw);
        }
        
        boolean exWtsPassed = true;
        boolean inWtsPassed = true;
        boolean inDegsPassed = true;
        boolean outDegsPassed = true;
        
        for (int i = 0, n = ng.size(); i < n; i++) {
            if (exWtsPassed) {
                exWtsPassed = exWts[i] == exWtsRw[i];
            } 
            if (inWtsPassed) {
                inWtsPassed = inWts[i] == inWtsRw[i];
            }
            if (inDegsPassed) {
                inDegsPassed = inDegs[i] == inDegsRw[i];
            }
            if (outDegsPassed) {
                outDegsPassed = outDegs[i] == outDegsRw[i];
            }
        }
        System.out.print("Excitatory strengths: ");
        if (exWtsPassed) {
            System.out.println("PASSED!");
        } else {
            System.out.println("FAILED!");
            System.out.println(Arrays.toString(exWts));
            System.out.println(Arrays.toString(exWtsRw));
        }
        System.out.print("Inhibitory strengths: ");
        if (inWtsPassed) {
            System.out.println("PASSED!");
        } else {
            System.out.println("FAILED!");
            System.out.println(Arrays.toString(inWts));
            System.out.println(Arrays.toString(inWtsRw));
        }
        System.out.print("In Degree: ");
        if (inDegsPassed) {
            System.out.println("PASSED!");
        } else {
            System.out.println("FAILED!");
            System.out.println(Arrays.toString(inDegs));
            System.out.println(Arrays.toString(inDegsRw));
        }
        System.out.print("Out Degree: ");
        if (outDegsPassed) {
            System.out.println("PASSED!");
        } else {
            System.out.println("FAILED!");
            System.out.println(Arrays.toString(outDegs));
            System.out.println(Arrays.toString(outDegsRw));
        }
    }
    
}
