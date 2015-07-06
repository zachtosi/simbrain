package org.simbrain.special_projects.plasticity_proj;

import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.time.LocalDateTime;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Random;
import java.util.Scanner;

import javax.swing.JButton;
import javax.swing.JFrame;
import javax.swing.JPanel;

import org.simbrain.network.connections.AllToAll;
import org.simbrain.network.connections.ConnectNeurons;
import org.simbrain.network.connections.Sparse;
import org.simbrain.network.core.Network;
import org.simbrain.network.core.Neuron;
import org.simbrain.network.core.SpikingNeuronUpdateRule;
import org.simbrain.network.core.Synapse;
import org.simbrain.network.groups.NeuronGroup;
import org.simbrain.network.groups.SynapseGroup;
import org.simbrain.network.neuron_update_rules.SpikingThresholdRule;
import org.simbrain.network.synapse_update_rules.STDPRule;
import org.simbrain.network.synapse_update_rules.StaticSynapseRule;
import org.simbrain.network.synapse_update_rules.spikeresponders.ConvolvedJumpAndDecay;
import org.simbrain.network.synapse_update_rules.spikeresponders.UDF;
import org.simbrain.network.update_actions.ConcurrentBufferedUpdate;
import org.simbrain.util.SimbrainConstants.Polarity;
import org.simbrain.util.math.ProbDistribution;
import org.simbrain.util.randomizer.PolarizedRandomizer;
import org.simbrain.util.randomizer.Randomizer;

import com.jmatio.io.MatFileWriter;
import com.jmatio.types.MLDouble;

public class MainSimulation {
    /*
     * _______________________________________________________________________ |
     * General Sim parameters |
     * _______________________________________________________________________|
     */
    private static int NUM_NEURONS = 540;
    private static int GRID_SPACE = 100; // influences delay times
    private static int NUM_INPUTS = 4;
    private static String IN_FILE_NAME = "./ZSim_Inputs/Generalization_Test_Inputs"
            + "_4ch_20Hz_200ms_05ms_bins1.csv";
    private static String TH_FILE_NAME = "LAIP_"
            + LocalDateTime.now().toString() + "_Thresholds.mat";
    private static String PF_FILE_NAME = "LAIP_"
            + LocalDateTime.now().toString() + "_PrefFRs.mat";
    private static String FR_FILE_NAME = "LAIP_"
            + LocalDateTime.now().toString() + "_FiringRates.mat";
    private static double[][] TEST_DATA;
    private static double[][] TH_VALS;
    private static double[][] PF_VALS;
    private static int REC_INTERVAL = 1; // seconds
    private static int TEST_RUNS = 40;

    /*
     * _______________________________________________________________________ |
     * Network parameters |
     * _______________________________________________________________________|
     */
    private static Network NETWORK;
    private static NeuronGroup LAIP_RES;
    private static NeuronGroup INPUT_LAYER;
    private final static double TIME_STEP = 0.5; // ms
    static {
        // REC_INTERVAL secs for timeStep/iteration iterations
        // Translates seconds into number of network update iterations
        // based on the network time step.
        REC_INTERVAL = (int) (REC_INTERVAL / (TIME_STEP / 1000));
    }
    private static final double INHIB_RATIO = 0.2; // 20%

    /*
     * _______________________________________________________________________ |
     * Synapse and synaptic plasticity parameters |
     * _______________________________________________________________________|
     */
    private static SynapseGroup LAIP_RES_SYNS;
    private static final boolean II_STDP_HEB = true;
    private static final boolean IE_STDP_HEB = true;
    private static final boolean USE_ADAPTIVE_SYMMETRY = true;
    private static final double I_LR = 0.00005;
    private static final double E_LR = 0.00005;
    private static final double II_W_PLUS = 4.8;
    private static final double II_W_MINUS = 1;
    private static final double IE_W_PLUS = 5.4;
    private static final double IE_W_MINUS = 1;
    private static final double E_W_PLUS = 4.8;
    private static final double E_W_MINUS = 1;
    private static final double II_TAU_PLUS = 20;
    private static final double II_TAU_MINUS = 50;
    private static final double IE_TAU_PLUS = 20;
    private static final double IE_TAU_MINUS = 50;
    private static final double E_TAU_PLUS = 20;
    private static final double E_TAU_MINUS = 100;

    /*
     * Values that all synapses take on at the very begining. This should be low
     * so as to emulate a total lack of connectivity. Synaptic growth or
     * depression through STDP combined with SynapsePruning causes some of these
     * synapses to "grow" by taking on values > 0 and others to wither.
     */
    private static final double startingWtExFl = 0.01;
    private static final double startingWtExCl = 0.011;
    private static final double startingWtInFl = 0.01;
    private static final double startingWtInCl = 0.011;

    /*
     * _______________________________________________________________________ |
     * Input parameters |
     * _______________________________________________________________________|
     */
    private static final ProbDistribution IN_SYN_DIST = ProbDistribution.LOGNORMAL;
    private static final double DIST_PARAM_1 = 1.5; // "Location" if dist is
                                                    // lognormal
    private static final double DIST_PARAM_2 = 0.75; // "Scale" if dist is
                                                     // lognormal
    // 1 meaning all possible connections exist from the input neurons to the
    // output neurons, and 0 meaning that no synapses will be constructed
    // between the input and the output.
    private static double IN_DENSITY = 1; // 0.05 when using > 4 inputs

    /*
     * _______________________________________________________________________ |
     * Neuron and neuronal plasticity parameters |
     * _______________________________________________________________________|
     */
    private static final double NOISE_MEAN = 0; // nA
    private static final double NOISE_STD = 0.2; // nA
    private static final double INHIB_TAU = 20; // ms
    private static final double EXCITE_TAU = 30; // ms
    private static final double INHIB_REF = 2.0; // ms
    private static final double EXCITE_REF = 3.0; // ms
    private static final double I_BG = 13.5; // nA

    /*
     * Homeostatic and Intrinsic plasticity (IP) parameters
     */
    private static final boolean USING_IP = true; // Turns on or off all
                                                  // neuronal plasticity

    private static final int DEFAULT_INIT_PFR = 1; // Hz or Spks/s
    // private static final double FR_ESTIMATE_TAU = 1000; // ms
    private static final double INITIAL_HOMEOSTATIC_CONST = 1E7;
    private static final double FINAL_HOMEOSTATIC_CONST = 1E5;
    private static final double HOMEOSTATIC_COOLING_RATE = 1E-4;
    private static final double INITIAL_IP_CONST = .1;
    private static final double FINAL_IP_CONST = 1E-7;
    private static final double INTRINSIC_COOLING_RATE = 1E-4;
    private static final double ALPHA_IP = 1;
    private static final double BETA_IP = 60;
    private static final double LOW_FR_BOUNDARY = 5; // Hz or Spks/s

    /**
     * ________________________________________________________________________
     * MAIN
     * 
     * @param args
     *            ________________________________________________________________________
     */
    public static void main(String[] args) {
        buildNetwork();
        TH_VALS = new double[1 + (TEST_RUNS * TEST_DATA.length / REC_INTERVAL)]
                [NUM_NEURONS];
        PF_VALS = new double[1 + (TEST_RUNS * TEST_DATA.length / REC_INTERVAL)]
                [NUM_NEURONS];
        JPanel dum = new JPanel();
        JButton dumpButton = new JButton("Dump Stats");
        dumpButton.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                try {
                    printNetworkStats();
                } catch (IOException e1) {
                    e1.printStackTrace();
                    System.out.println("PROGRAM IS STILL RUNNING!");
                }
            }
        });
        dum.add(dumpButton);
        JFrame frame = new JFrame();
        frame.setContentPane(dum);
        frame.pack();
        frame.setVisible(true);
        // initializeSynapseWatching();
        printParams();
        System.out.println(((IPIFRule) LAIP_RES.getNeuronList()
                .get(0).getUpdateRule()).getResetPotential());
        Integer[][] delays = new Integer[LAIP_RES.size()][LAIP_RES.size()];

        delays = (Integer[][]) LAIP_RES_SYNS.getMatrixOfValues(
                (Synapse s) -> s.getDelay(), delays);
        double[][] ddlys = new double[LAIP_RES.size()][LAIP_RES.size()];
        for (int i = 0; i < delays.length; i++) {
            for (int j = 0; j < delays.length; j++) {
                if (i != j)
                    ddlys[i][j] = delays[i][j];
            }
        }

        MLDouble dlys = new MLDouble("Delays", ddlys);
        try {
            new MatFileWriter("SynDelays"
                    + ".mat", Collections.singletonList(dlys));
        } catch (IOException ie) {
            ie.printStackTrace();
        }
        runAnnealingProcedure();
        // printSynChanges(2);
        System.exit(0);
    }

    /**
     * Run the network forward in time using LAIP (or not)
     */
    public static void runAnnealingProcedure() {
        System.out.println(TEST_RUNS * TEST_DATA.length / REC_INTERVAL);
        int c = 0;
        double[] frs = new double[NUM_NEURONS];
        for (Neuron neuron : LAIP_RES.getNeuronList()) {
            ((IPIFRule) neuron.getUpdateRule())
                    .setUsingIP(USING_IP);
        }
        try {
            for (int i = 0, n = TEST_RUNS * TEST_DATA.length; i < n; i++) {
                c = i;

                // Every "REC_INTERVAL" usually 1s, record some data
                if (i % REC_INTERVAL == 0) {
                    List<Neuron> neurons = LAIP_RES.getNeuronList();
                    for (int j = 0; j < NUM_NEURONS; j++) {
                        TH_VALS[i / REC_INTERVAL][j] = ((IPIFRule)
                                neurons.get(j).getUpdateRule()).getThreshold();
                        PF_VALS[i / REC_INTERVAL][j] = ((IPIFRule)
                                neurons.get(j).getUpdateRule()).getPrefFR();
                    }
                }
                
                if (i % TEST_DATA.length == 0) {
                    // Create a new file for spike times
                    LAIP_RES.stopRecording();
                    LAIP_RES.startRecording(new File(
                            ((double) i / REC_INTERVAL)
                                    + "-"
                                    + ((double) (i + TEST_DATA.length)
                                            / REC_INTERVAL) + "s.csv"));
                }

                // Turn off STDP & IP for last test run
                if (i == (TEST_RUNS - 1) * TEST_DATA.length) {
                    LAIP_RES_SYNS.setLearningRule(new StaticSynapseRule(),
                            Polarity.BOTH);
                    for (Neuron neuron : LAIP_RES.getNeuronList()) {
                        ((IPIFRule) neuron.getUpdateRule())
                                .setUsingIP(false);
                    }
                }

                // Update the network
                synchronized (NETWORK) {
                    NETWORK.update();
                }
                
                // Print out which simulated second we're on
                if (i % REC_INTERVAL == 0) {
                    System.out.println(i / REC_INTERVAL);
                }

            }
        } finally {
            LAIP_RES.stopRecording();
            int k = 0;
            for (Neuron n : LAIP_RES.getNeuronList()) {
                frs[k++] = ((SpikingNeuronUpdateRule) n.getUpdateRule())
                        .getFrequency(n);
            }
            System.out.println(c + " " + TEST_RUNS * TEST_DATA.length);
            LAIP_RES.stopRecording();
            MLDouble spks = new MLDouble("FiringRates", frs, 1);
            MLDouble th = new MLDouble("thresholds",
                    TH_VALS);
            MLDouble pf = new MLDouble("PFs",
                    PF_VALS);
            MLDouble mat = new MLDouble("wts",
                    LAIP_RES_SYNS.getWeightMatrix());
            try {
                new MatFileWriter(FR_FILE_NAME,
                        Collections.singletonList(spks));
                new MatFileWriter(TH_FILE_NAME,
                        Collections.singletonList(th));
                new MatFileWriter(PF_FILE_NAME,
                        Collections.singletonList(pf));
                new MatFileWriter("LAIPWtMat.mat",
                        Collections.singletonList(mat));
            } catch (IOException e) {
                e.printStackTrace();
            }
        }
    }

    /**
     * Prints various network stats: -weight matrix -firing rates -preferred
     * firing rates -thresholds to .mat files
     * 
     * @throws IOException
     */
    public static void printNetworkStats() throws IOException {
        synchronized (NETWORK) {
            double time = LAIP_RES.getParentNetwork().getTime();
            MLDouble mat = new MLDouble("wtMat",
                    LAIP_RES_SYNS.getWeightMatrix());
            new MatFileWriter("LAIPWtMatCheckIn"
                    + time + ".mat", Collections.singletonList(mat));
            int k = 0;
            double[] frs = new double[LAIP_RES.size()];
            double[] pfrs = new double[LAIP_RES.size()];
            double[] ths = new double[LAIP_RES.size()];
            System.out.println(((IPIFRule) LAIP_RES
                    .getNeuronList().get(0).getUpdateRule()).getLearningRate());
            System.out.println(((IPIFRule) LAIP_RES
                    .getNeuronList().get(0).getUpdateRule()).getIpConst());
            for (Neuron neuron : LAIP_RES.getNeuronList()) {
                frs[k] = ((SpikingNeuronUpdateRule) neuron
                        .getUpdateRule())
                        .getFrequency(neuron);
                ((SpikingNeuronUpdateRule) neuron.getUpdateRule())
                        .resetCounter(neuron);
                pfrs[k] = ((IPIFRule) neuron
                        .getUpdateRule())
                        .getPrefFR();
                ths[k] = ((IPIFRule) neuron
                        .getUpdateRule())
                        .getThreshold();
                ((SpikingNeuronUpdateRule) neuron.getUpdateRule())
                        .resetCounter(neuron);
                k++;
            }
            MLDouble spks = new MLDouble("FiringRates", frs, 1);
            new MatFileWriter(time + FR_FILE_NAME,
                    Collections.singletonList(spks));
            MLDouble pref = new MLDouble("PrefFRs", pfrs, 1);
            new MatFileWriter(time + "PrefFRs.mat",
                    Collections.singletonList(pref));
            MLDouble thresh = new MLDouble("Thresholds", ths, 1);
            new MatFileWriter(time + "Thresholds.mat",
                    Collections.singletonList(thresh));
        }
    }

    // public static void printSynChanges(int no) {
    // ArrayList<Synapse> taggedSyns = new ArrayList<Synapse>();
    // for (Synapse s : LAIP_RES_SYNS.getAllSynapses()) {
    // if (s.tag)
    // taggedSyns.add(s);
    // }
    // double [][] strChanges = new double[taggedSyns.size()]
    // [taggedSyns.get(0).strChanges.length];
    // System.out.println(strChanges.length);
    // for (int i = 0, n = strChanges.length; i < n; i++) {
    // for (int j = 0, m = strChanges[i].length; j < m; j++) {
    // strChanges[i][j] = taggedSyns.get(i).strChanges[j];
    // }
    // }
    // MLDouble synCh = new MLDouble("synCh", strChanges);
    // try {
    // new MatFileWriter("synChanges" + no + ".mat", Collections
    // .singletonList(synCh));
    // } catch (IOException ie) {
    // ie.printStackTrace();
    // }
    // }

    public static void printParams() {
        try (PrintWriter pw = new PrintWriter(new FileWriter("Parameters.txt"));
                Scanner sc = new Scanner(System.in);) {
            System.out.println("Special notes?");
            String notes = sc.nextLine();
            pw.println("Network Parameters");
            pw.println("RES Neurons: " + NUM_NEURONS);
            pw.println("Input Neurons: " + NUM_INPUTS);
            pw.println("Time-step: " + TIME_STEP);
            pw.println("Inhib %: " + INHIB_RATIO);
            pw.println("IP On: " + USING_IP);
            pw.println();
            pw.println("Synapse Parameters");
            pw.println("II_STDP Hebbian: " + II_STDP_HEB);
            pw.println("IE_STDP Hebbian: " + IE_STDP_HEB);
            pw.println("E-STDP Learning Rate: " + E_LR);
            pw.println("I-STDP Learning Rate: " + I_LR);
            pw.println("II_W_PLUS: " + II_W_PLUS);
            pw.println("II_W_MINUS: " + II_W_MINUS);
            pw.println("IE_W_PLUS: " + IE_W_PLUS);
            pw.println("IE_W_MINUS: " + IE_W_MINUS);
            pw.println("E_W_PLUS: " + E_W_PLUS);
            pw.println("E_W_MINUS: " + E_W_MINUS);
            pw.println("Start wt val:" + startingWtExFl);
            pw.println();
            pw.println("Neuron Parameters");
            IPIFRule example = (IPIFRule) LAIP_RES
                    .getNeuronList().get(0).getUpdateRule();
            pw.println("Alpha: " + example.getAlpha());
            pw.println("StartIP Const.: " + example.getIpConst());
            pw.println("Start Learning Rate: " + example.getLearningRate());
            pw.println("Beta: " + example.getBeta());
            pw.println("Small FR Threshold: " + example.getLowFRBoundary());
            pw.println("Tau_a: " + example.getTauA());
            pw.println();
            pw.println(notes);
        } catch (IOException ie) {
            ie.printStackTrace();
        }
    }

    /**
     * Constructs the network
     */
    private static void buildNetwork() {
        NETWORK = new Network();
        NETWORK.setTimeStep(TIME_STEP);

        List<Neuron> neurons = new ArrayList<Neuron>(NUM_NEURONS);
        Randomizer rand = new Randomizer();
        rand.setPdf(ProbDistribution.NORMAL);
        rand.setParam1(NOISE_MEAN);
        rand.setParam2(NOISE_STD);
        Random randi = new Random();
        for (int i = 0; i < NUM_NEURONS; i++) {
            Neuron neuron = new Neuron(NETWORK);
            randomizeNeuronLocation(randi, neuron);
            neurons.add(neuron);

            // Set up neurons
            IPIFRule lif = new IPIFRule();
            // lif.setTauA(FR_ESTIMATE_TAU);
            lif.setIpConst(INITIAL_HOMEOSTATIC_CONST);
            lif.setIpConstFinal(FINAL_HOMEOSTATIC_CONST);
            lif.setLearningRate(INITIAL_IP_CONST);
            lif.setLearningRateFinal(FINAL_IP_CONST);
            lif.setIntrinsicCooling(INTRINSIC_COOLING_RATE);
            lif.setHomeostaticCooling(HOMEOSTATIC_COOLING_RATE);
            lif.setAlpha(ALPHA_IP);
            lif.setBeta(BETA_IP);
            lif.setLowFRBoundary(LOW_FR_BOUNDARY);
            lif.setPrefFR(DEFAULT_INIT_PFR);

            if (Math.random() < INHIB_RATIO) {
                neuron.setPolarity(Polarity.INHIBITORY);
                lif.setRefractoryPeriod(INHIB_REF);
                lif.setTimeConstant(INHIB_TAU);
            } else {
                neuron.setPolarity(Polarity.EXCITATORY);
                lif.setRefractoryPeriod(EXCITE_REF);
                lif.setTimeConstant(EXCITE_TAU);
            }
            lif.setBackgroundCurrent(I_BG);
            lif.setAddNoise(true);
            lif.setNoiseGenerator(rand);
            neuron.setUpdateRule(lif);
            neuron.randomize();
        }
        LAIP_RES = new NeuronGroup(NETWORK, neurons);
        LAIP_RES.setLabel("LAIP_Reservoir");

        NETWORK.addGroup(LAIP_RES);

        PolarizedRandomizer exRand = new PolarizedRandomizer(
                Polarity.EXCITATORY,
                ProbDistribution.UNIFORM);
        PolarizedRandomizer inRand = new PolarizedRandomizer(
                Polarity.INHIBITORY,
                ProbDistribution.UNIFORM);
        exRand.setParam1(startingWtExFl);
        exRand.setParam2(startingWtExCl);
        inRand.setParam1(startingWtInFl);
        inRand.setParam2(startingWtInCl);

        AllToAll con = new AllToAll(false);
        LAIP_RES_SYNS = SynapseGroup.createSynapseGroup(LAIP_RES, LAIP_RES,
                con, 0.5, exRand, inRand);
        LAIP_RES_SYNS.setSpikeResponder(new UDF(), Polarity.BOTH);
        LAIP_RES_SYNS.setLabel("Recurrent Synapses");

        for (Synapse s : LAIP_RES_SYNS.getAllSynapses()) {
            double dist = Network
                    .getEuclideanDist(s.getSource(), s.getTarget());
            double delay = Math.round(dist / (0.6 * NUM_NEURONS));

            s.setDelay((int) delay);
            if (Math.random() < 0.01) {
                System.out.println(delay * 0.5);
            }
        }

        STDPRule dummy = new STDPRule();
        LAIP_RES_SYNS.setLearningRule(dummy, Polarity.INHIBITORY);

        STDPRule iiSTDP = USE_ADAPTIVE_SYMMETRY ? new SymmetricSTDPRule()
                : new STDPRule();
        iiSTDP.setLearningRate(I_LR);
        iiSTDP.setTau_minus(II_TAU_MINUS);
        iiSTDP.setTau_plus(II_TAU_PLUS);
        iiSTDP.setW_plus(II_W_PLUS);
        iiSTDP.setW_minus(II_W_MINUS);
        iiSTDP.setHebbian(II_STDP_HEB);

        STDPRule ieSTDP = USE_ADAPTIVE_SYMMETRY ? new SymmetricSTDPRule()
                : new STDPRule();
        ieSTDP.setLearningRate(I_LR);
        ieSTDP.setTau_minus(IE_TAU_MINUS);
        ieSTDP.setTau_plus(IE_TAU_PLUS);
        ieSTDP.setW_plus(IE_W_PLUS);
        ieSTDP.setW_minus(IE_W_MINUS);
        ieSTDP.setHebbian(IE_STDP_HEB);
        for (Neuron n : neurons) {
            if (n.getPolarity() == Polarity.INHIBITORY) {
                for (Synapse s : n.getFanOut().values()) {
                    if (s.getTarget().getPolarity() == Polarity.EXCITATORY) {
                        s.setLearningRule(ieSTDP.deepCopy());
                    } else {
                        s.setLearningRule(iiSTDP.deepCopy());
                    }
                }
            }
        }

        STDPRule eeSTDP = USE_ADAPTIVE_SYMMETRY ? new SymmetricSTDPRule()
                : new STDPRule();
        eeSTDP.setLearningRate(E_LR);
        eeSTDP.setTau_minus(E_TAU_MINUS);
        eeSTDP.setTau_plus(E_TAU_PLUS);
        eeSTDP.setW_plus(E_W_PLUS);
        eeSTDP.setW_minus(E_W_MINUS);
        LAIP_RES_SYNS.setLearningRule(eeSTDP, Polarity.EXCITATORY);
        NETWORK.addGroup(LAIP_RES_SYNS);

        // Clears all the serial buffered updates
        NETWORK.getUpdateManager().clear();

        INPUT_LAYER = new NeuronGroup(NETWORK, NUM_INPUTS);
        INPUT_LAYER.setNeuronType(new SpikingThresholdRule());
        INPUT_LAYER.setTestData(getData(IN_FILE_NAME));
        INPUT_LAYER.setInputMode(true);
        NETWORK.addGroup(INPUT_LAYER);
        ConnectNeurons sCon;
        if (IN_DENSITY < 1) {
            sCon = new Sparse(IN_DENSITY, true, false);
        } else {
            sCon = new AllToAll(false);
        }
        PolarizedRandomizer inputRand = new PolarizedRandomizer(Polarity
                .EXCITATORY, IN_SYN_DIST);
        inputRand.setParam1(DIST_PARAM_1);
        inputRand.setParam2(DIST_PARAM_2);
        SynapseGroup in_sg = SynapseGroup.createSynapseGroup(INPUT_LAYER,
                LAIP_RES, sCon, 1.0, inputRand, inRand);
        in_sg.setSpikeResponder(new ConvolvedJumpAndDecay(), Polarity.BOTH);
        NETWORK.addGroup(in_sg);
        NETWORK.getUpdateManager().removeGroupAction(in_sg);

        LAIP_RES_SYNS.setUpperBound(100000, Polarity.EXCITATORY);
        LAIP_RES_SYNS.setLowerBound(0, Polarity.EXCITATORY);
        LAIP_RES_SYNS.setLowerBound(-100000, Polarity.INHIBITORY);
        LAIP_RES_SYNS.setUpperBound(0, Polarity.INHIBITORY);
        NETWORK.fireSynapsesUpdated();

        // Add the prune synapses action
        PrunerAction pruner = new PrunerAction(LAIP_RES_SYNS);
        NETWORK.getUpdateManager().addAction(pruner);
        NETWORK.getUpdateManager().addAction(ConcurrentBufferedUpdate
                .createConcurrentBufferedUpdate(NETWORK));
        NETWORK.updateTimeType();

    }

    // private static void initializeSynapseWatching() {
    // for (Synapse s : LAIP_RES_SYNS.getAllSynapses()) {
    // if (Math.random() < 0.0025) {
    // s.tag = true;
    // s.strChanges = new float[TEST_RUNS
    // * TEST_DATA.length / s.chkInterval];
    // }
    // }
    //
    //
    // }

    private static void randomizeNeuronLocation(Random randi, Neuron neuron) {
        neuron.setX(randi.nextInt((int) (Math.sqrt(NUM_NEURONS)
                * GRID_SPACE)));
        neuron.setY(randi.nextInt((int) (Math.sqrt(NUM_NEURONS)
                * GRID_SPACE)));
        neuron.setZ(randi.nextInt((int) (Math.sqrt(NUM_NEURONS)
                * GRID_SPACE)));
    }

    private static double[][] getData(String filename) {
        try (FileReader fr = new FileReader(filename);
                Scanner lineScan = new Scanner(fr);) {
            ArrayList<ArrayList<Double>> dynData =
                    new ArrayList<ArrayList<Double>>();
            while (lineScan.hasNextLine()) {
                Scanner reader = new Scanner(lineScan.nextLine());
                reader.useDelimiter(", ");
                ArrayList<Double> line = new ArrayList<Double>();
                while (reader.hasNextDouble()) {
                    line.add(reader.nextDouble());
                }
                dynData.add(line);
                reader.close();
            }
            double[][] data = new double[dynData.size()][];
            for (int i = 0, n = data.length; i < n; i++) {
                data[i] = new double[dynData.get(i).size()];
                for (int j = 0, m = data[i].length; j < m; j++) {
                    data[i][j] = dynData.get(i).get(j);
                }
            }
            TEST_DATA = data;
            return data;
        } catch (IOException ie) {
            ie.printStackTrace();
            System.err.println("File read in failed.");
            return null;
        }
    }

}
