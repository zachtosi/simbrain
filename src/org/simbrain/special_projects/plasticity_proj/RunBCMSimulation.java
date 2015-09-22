package org.simbrain.special_projects.plasticity_proj;

import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.lang.reflect.Field;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Date;
import java.util.HashSet;
import java.util.List;
import java.util.Random;
import java.util.Scanner;
import java.util.Set;
import java.util.concurrent.ThreadLocalRandom;

import javax.swing.JButton;
import javax.swing.JFrame;
import javax.swing.JPanel;

import org.simbrain.network.connections.AllToAll;
import org.simbrain.network.connections.ConnectNeurons;
import org.simbrain.network.connections.ConnectionUtilities.SynapseParameterGetter;
import org.simbrain.network.connections.Sparse;
import org.simbrain.network.core.Network;
import org.simbrain.network.core.NetworkUpdateAction;
import org.simbrain.network.core.Neuron;
import org.simbrain.network.core.SpikingNeuronUpdateRule;
import org.simbrain.network.core.Synapse;
import org.simbrain.network.groups.NeuronGroup;
import org.simbrain.network.groups.SynapseGroup;
import org.simbrain.network.neuron_update_rules.SpikingThresholdRule;
import org.simbrain.network.synapse_update_rules.StaticSynapseRule;
import org.simbrain.network.synapse_update_rules.spikeresponders.ConvolvedJumpAndDecay;
import org.simbrain.network.update_actions.ConcurrentBufferedUpdate;
import org.simbrain.util.SimbrainConstants.Polarity;
import org.simbrain.util.math.ProbDistribution;
import org.simbrain.util.randomizer.PolarizedRandomizer;
import org.simbrain.util.randomizer.Randomizer;

import com.jmatio.io.MatFileWriter;
import com.jmatio.types.MLArray;
import com.jmatio.types.MLDouble;

public class RunBCMSimulation {
    /*
     * _______________________________________________________________________|
     * General Sim parameters |
     * _______________________________________________________________________|
     */
    private static int NUM_NEURONS = 1000;
    private static int GRID_SPACE = 200; // influences delay times
    private static int NUM_INPUTS = 4;
    private static String IN_FILE_NAME = "./ZSim_Inputs/Generalization_Test"
            + "_Inputs_" + NUM_INPUTS + "ch_20Hz_200ms_05ms_bins0.csv";
    private static File OUTPUT_DIR;
//    static {
//        Date dNow = new Date( );
//        SimpleDateFormat ft = 
//        new SimpleDateFormat ("E yyyy.MM.dd 'at' hh:mm:ss a zzz");
//        String dirName = "results: " + ft.format(dNow);
//        File mainOut = new File("./Outputs");
//        if (!mainOut.exists()) {
//            if (!mainOut.mkdir()) {
//                System.err.println("FATAL ERROR: FAILED TO CREATE MAIN "
//                        + "OUTPUT DIRECTORY.");
//                System.exit(1);
//            }
//        }
//        OUTPUT_DIR = new File(mainOut, dirName);
//        if (!OUTPUT_DIR.mkdir()) {
//            System.err.println("FATAL ERROR: FAILED TO CREATE INSTANCE "
//                    + "OUTPUT DIRECTORY.");
//            System.exit(1);
//        }
//    }
    private static String TH_FILE_NAME = "LAIP_Thresholds.mat";
    private static String PF_FILE_NAME = "LAIP_PrefFRs.mat";
    private static String FR_FILE_NAME = "LAIP_FiringRates.mat";
    private static double[][] TEST_DATA;
    private static double[][] TH_VALS;
    private static double[][] PF_VALS;
    private static int REC_INTERVAL = 1; // seconds
    private static int TEST_RUNS = 30;
    private static ConcurrentBufferedUpdate NET_UPDATER;
    private static PrunerAction PRUNER;

    /*
     * _______________________________________________________________________|
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
     * _______________________________________________________________________|
     * Synapse and synaptic plasticity parameters |
     * _______________________________________________________________________|
     */
    private static SynapseGroup LAIP_RES_SYNS;
    private static SynapseGroup IN_SG;
//    private static final boolean II_STDP_HEB = true;
//    private static final boolean IE_STDP_HEB = true;
    private static final double I_LR = 5E-2; //2E-6
    private static final double E_LR = 5E-2; //2E-6
//    private static final double II_W_PLUS = 5;
//    private static final double II_W_MINUS = 1;
//    private static final double IE_W_PLUS = 6;
//    private static final double IE_W_MINUS = 1.5;
//    private static final double E_W_PLUS = 5;
//    private static final double E_W_MINUS = 1;
//    private static final double II_TAU_PLUS = 20;
//    private static final double II_TAU_MINUS = 100;
//    private static final double IE_TAU_PLUS = 20;
//    private static final double IE_TAU_MINUS = 40;
//    private static final double E_TAU_PLUS = 20;
//    private static final double E_TAU_MINUS = 100;
    private static final ConvolvedJumpAndDecay E_SPK_RESP
        = new ConvolvedJumpAndDecay();
    private static final ConvolvedJumpAndDecay I_SPK_RESP
        = new ConvolvedJumpAndDecay();
    static {
        E_SPK_RESP.setTimeConstant(3);
        I_SPK_RESP.setTimeConstant(6);
    }
    private static final ModifiedBCMRule E_SYN_PLASTICITY =
            new ModifiedBCMRule();
    private static final ModifiedBCMRule II_SYN_PLASTICITY =
            new ModifiedBCMRule();
    private static final ModifiedBCMRule IE_SYN_PLASTICITY =
            new ModifiedBCMRule();
    
//    private static final SymmetricSTDPRule E_SYN_PLASTICITY =
//            new SymmetricSTDPRule();
//    private static final SymmetricSTDPRule II_SYN_PLASTICITY =
//            new SymmetricSTDPRule();
//    private static final SymmetricSTDPRule IE_SYN_PLASTICITY =
//            new SymmetricSTDPRule(); 
    
    static {
        E_SYN_PLASTICITY.setLearningRate(E_LR);
        II_SYN_PLASTICITY.setLearningRate(I_LR);
        IE_SYN_PLASTICITY.setLearningRate(I_LR);
    }
    
//    static {
//        E_SYN_PLASTICITY.setLearningRate(E_LR);
//        E_SYN_PLASTICITY.setHebbian(true);
//        E_SYN_PLASTICITY.setSymmetric(true);
//        E_SYN_PLASTICITY.setSymTau(15);
//        E_SYN_PLASTICITY.setWtSoftening(50);
//        E_SYN_PLASTICITY.setTau_plus(E_TAU_PLUS);
//        E_SYN_PLASTICITY.setTau_minus(E_TAU_MINUS);
//        E_SYN_PLASTICITY.setW_minus(E_W_MINUS);
//        E_SYN_PLASTICITY.setW_plus(E_W_PLUS);
//
//        IE_SYN_PLASTICITY.setLearningRate(I_LR);
//        IE_SYN_PLASTICITY.setHebbian(true);
//        IE_SYN_PLASTICITY.setSymmetric(true);
//        IE_SYN_PLASTICITY.setSymTau(5);
//        IE_SYN_PLASTICITY.setWtSoftening(50);
//        IE_SYN_PLASTICITY.setTau_plus(IE_TAU_PLUS);
//        IE_SYN_PLASTICITY.setTau_minus(IE_TAU_MINUS);
//        IE_SYN_PLASTICITY.setW_minus(IE_W_MINUS);
//        IE_SYN_PLASTICITY.setW_plus(IE_W_PLUS);
//        
//        II_SYN_PLASTICITY.setLearningRate(I_LR);
//        II_SYN_PLASTICITY.setHebbian(true);
//        II_SYN_PLASTICITY.setSymmetric(true);
//        II_SYN_PLASTICITY.setSymTau(5);
//        II_SYN_PLASTICITY.setWtSoftening(50);
//        II_SYN_PLASTICITY.setTau_plus(II_TAU_PLUS);
//        II_SYN_PLASTICITY.setTau_minus(II_TAU_MINUS);
//        II_SYN_PLASTICITY.setW_minus(II_W_MINUS);
//        II_SYN_PLASTICITY.setW_plus(II_W_PLUS);
//    }

    /*
     * Values that all synapses take on at the very begining. This should be low
     * so as to emulate a total lack of connectivity. Synaptic growth or
     * depression through STDP combined with SynapsePruning causes some of these
     * synapses to "grow" by taking on values > 0 and others to wither.
     */
    private static final double startingWtExFl = .1; //0.01
    private static final double startingWtExCl = .11;
    private static final double startingWtInFl = .1; //0.01
    private static final double startingWtInCl = .11;

    /*
     * _______________________________________________________________________ |
     * Input parameters |
     * _______________________________________________________________________|
     */
    private static final ProbDistribution IN_SYN_DIST = ProbDistribution.LOGNORMAL;
    private static final double DIST_PARAM_1 = 3; // "Location" if dist is
                                                  // lognormal
    private static final double DIST_PARAM_2 = 4; // "Scale" if dist is
                                                  // lognormal
    // 1 meaning all possible connections exist from the input neurons to the
    // output neurons, and 0 meaning that no synapses will be constructed
    // between the input and the output.
    private static double IN_DENSITY = 1; // make this closer to 1 for less
    // inputs and closer to 0 for less inputs/more res neurons

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
    private static final double I_BG = 13; // nA
    private static final double RESET_POTENTIAL = 13; // mV, use 10 for larger
                                                        // sims

    /*
     * Homeostatic and Intrinsic plasticity (IP) parameters
     */
    private static final boolean USING_IP = true; // Turns on or off all
                                                  // neuronal plasticity

    private static final double DEFAULT_INIT_PFR = 0.1; // Hz or Spks/s
    // private static final double FR_ESTIMATE_TAU = 1000; // ms
    private static final double INITIAL_HOMEOSTATIC_CONST = 1E7;
    private static final double FINAL_HOMEOSTATIC_CONST = 1E5;
    private static final double HOMEOSTATIC_COOLING_RATE = 1E-4;
    private static final double INITIAL_IP_CONST = 1;// 1;
    private static final double FINAL_IP_CONST = 1E-7;
    private static final double INTRINSIC_COOLING_RATE = 1E-4;
    private static final double ALPHA_IP = 1;// 5;
    private static final double BETA_IP = 400;
    private static final double LOW_FR_BOUNDARY = 0.5; // Hz or Spks/s

    /**
     * ________________________________________________________________________
     * MAIN
     * 
     * @param args
     *            ________________________________________________________________________
     */
    public static void main(String[] args) {
        System.out.println(Runtime.getRuntime().availableProcessors());
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
                new SynapseParameterGetter<Integer>() {
                    @Override
                    public Integer getParameterFromSynapse(Synapse synapse) {
                        return synapse.getDelay();
                    }
                }, delays);
        double[][] ddlys = new double[LAIP_RES.size()][LAIP_RES.size()];
        for (int i = 0; i < delays.length; i++) {
            for (int j = 0; j < delays.length; j++) {
                if (i != j)
                    ddlys[i][j] = delays[i][j];
            }
        }

        MLArray dlys = new MLDouble("Delays", ddlys);
        MLArray inMat = new MLDouble("InMat", IN_SG.getWeightMatrix());
        try {
            MatFileWriter imw = new MatFileWriter();
            imw.write("InMat.mat", Collections.singletonList(inMat));
            MatFileWriter mfw = new MatFileWriter();
            mfw.write("SynDelays"
                    + ".mat", Collections.singletonList(dlys));
        } catch (IOException ie) {
            ie.printStackTrace();
        }
        for (NetworkUpdateAction nua : NETWORK.getUpdateManager()
                .getActionList()) {
            System.out.println(nua.getDescription());
        }
        runAnnealingProcedure();
        // printSynChanges(2);
        System.exit(0);
    }

    /**
     * Run the network forward in time using LAIP (or not)
     */
    public static void runAnnealingProcedure() {
        long start = System.nanoTime();
        System.out.println(TEST_RUNS * TEST_DATA.length / REC_INTERVAL);
        int c = 0;
        boolean prune = true;
        double[] frs = new double[NUM_NEURONS];
        for (Neuron neuron : LAIP_RES.getNeuronList()) {
            ((IPIFRule) neuron.getUpdateRule())
                    .setUsingIP(USING_IP);
        }
        LAIP_RES.stopRecording();
        try {
            for (int i = 0, n = TEST_RUNS * TEST_DATA.length; i < n; i++) {
                c = i;

                if (i % (TEST_DATA.length * 2) == 0
                        && i < (30 * TEST_DATA.length)) {
                    try {
                        printNetworkStats();
                    } catch (IOException ie) {
                        ie.printStackTrace();
                    }
                }

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

                if (i % TEST_DATA.length == 0 && i > 20 * TEST_DATA.length) {
                    LAIP_RES_SYNS.setLearningRule(new StaticSynapseRule(),
                            Polarity.BOTH);
                    if (prune) {
                        prune = false;
                        NETWORK.getUpdateManager().removeAction(PRUNER);
                    }
                    // Create a new file for spike times
                    LAIP_RES.stopRecording();
                    LAIP_RES.startRecording(new File(
                            ((double) i / REC_INTERVAL)
                                    + "-"
                                    + ((double) (i + TEST_DATA.length)
                                    / REC_INTERVAL) + "s.csv"));
                }

                if (prune && LAIP_RES_SYNS.size() <
                        (NUM_NEURONS * (NUM_NEURONS - 1) * 0.05)) {
                    prune = false;
                    NETWORK.getUpdateManager().removeAction(PRUNER);
                    for (Neuron neu : LAIP_RES.getNeuronList()) {
                        ((IPIFRule) neu.getUpdateRule()).setFullNorm(true);
                    }
                }

                // Turn off STDP & IP for last test runs
                if (i == 24 * TEST_DATA.length && !prune) {
                    printNetworkStats();
                    System.out.println("REMOVING INPUT. ADJUSTING HOMEOSTASIS");
                    NETWORK.getUpdateManager().clear();
                    INPUT_LAYER.setInputMode(false);
                    NETWORK.removeGroup(INPUT_LAYER);
                    NETWORK.getUpdateManager()
                            .addAction(ConcurrentBufferedUpdate
                                    .createConcurrentBufferedUpdate(NETWORK));
                    for (Neuron neu : LAIP_RES.getNeuronList()) {
                        ((IPIFRule) neu.getUpdateRule()).setIpConst(1E5);
                        ((IPIFRule) neu.getUpdateRule()).setLearningRate(1E-8);
                        ((IPIFRule) neu.getUpdateRule()).setAnnealing(false);
                        ((IPIFRule) neu.getUpdateRule()).setResistance(1);
                        ((IPIFRule) neu.getUpdateRule())
                                .setBackgroundCurrent(13.5);
                        // if (neu.getPolarity() == Polarity.INHIBITORY) {
                        // ((IPIFRule) neu.getUpdateRule()).setUsingIP(false);
                        // }
                    }
                }
                // // Re-wire the synapse group for the last run.
                // if (i == (TEST_RUNS - 2) * TEST_DATA.length) {
                // // rewireSynapseGroup(LAIP_RES_SYNS, 100);
                // printNetworkStats();
                // }

                // Update the network
                synchronized (NETWORK) {
                    NETWORK.update();
                }
                
                for (int j = 0; j < NUM_NEURONS; j++) {
                    ((IPIFRule) LAIP_RES.getNeuronListUnsafe()
                            .get(j).getUpdateRule()).pushBuffers();
                }

                // Print out which simulated second we're on
                if (i % REC_INTERVAL == 0) {
                    System.out.println(i / REC_INTERVAL);
                }

            }
        } catch (IOException ie) {
            // TODO: something better here.
            ie.printStackTrace();
        } finally {
            long end = System.nanoTime();
            System.out.println((end - start) / 1E9 + " secs.");
            LAIP_RES.stopRecording();
            int k = 0;
            for (Neuron n : LAIP_RES.getNeuronList()) {
                frs[k++] = ((SpikingNeuronUpdateRule) n.getUpdateRule())
                        .getFrequency(n);
            }
            System.out.println(c + " " + TEST_RUNS * TEST_DATA.length);
            LAIP_RES.stopRecording();
            MLArray spks = new MLDouble("FiringRates", frs, 1);
            MLArray th = new MLDouble("thresholds",
                    TH_VALS);
            MLArray pf = new MLDouble("PFs",
                    PF_VALS);
            MLArray mat = new MLDouble("wts",
                    LAIP_RES_SYNS.getWeightMatrix());
            try {
                File f = new File(OUTPUT_DIR, FR_FILE_NAME);
                f.createNewFile();
                new MatFileWriter().write(f,Collections.singletonList(spks));
                File f1 = new File(OUTPUT_DIR, TH_FILE_NAME);
                f1.createNewFile();
                new MatFileWriter().write(f1, Collections.singletonList(th));
                File f2 = new File(OUTPUT_DIR, PF_FILE_NAME);
                f2.createNewFile();
                new MatFileWriter().write(f2, Collections.singletonList(pf));
                File f3 = new File(OUTPUT_DIR, "LAIPTwMat.mat");
                f3.createNewFile();
                new MatFileWriter().write(f3, Collections.singletonList(mat));
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
            MLArray mat = new MLDouble("wtMat",
                    LAIP_RES_SYNS.getWeightMatrix());
            File f1 = new File(OUTPUT_DIR, time + "LAIPWtMatCheckIn.mat");
            f1.createNewFile();
            new MatFileWriter(f1, Collections.singletonList(mat));
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
            MLArray spks = new MLDouble("FiringRates", frs, 1);
            File f2 = new File(OUTPUT_DIR, time + FR_FILE_NAME);
            f2.createNewFile();
            new MatFileWriter(f2, Collections.singletonList(spks));
            MLArray pref = new MLDouble("PrefFRs", pfrs, 1);
            File f3 = new File(OUTPUT_DIR, time + "PrefFRs.mat");
            f3.createNewFile();
            new MatFileWriter(f3, Collections.singletonList(pref));
            MLArray thresh = new MLDouble("Thresholds", ths, 1);
            File f4 = new File(OUTPUT_DIR, time + "Thresholds.mat");
            f4.createNewFile();
            new MatFileWriter(f4, Collections.singletonList(thresh));
        }
    }
    
    public static void printParams() {
        File paramF = new File(OUTPUT_DIR, "Parameters.txt");
       try {
           paramF.createNewFile();
       } catch (Exception e) {
           e.printStackTrace();
           System.exit(1);
       }
        try (   FileWriter fw = new FileWriter(paramF);
                PrintWriter pw = new PrintWriter(fw);) {
            pw.println("\t" + RunBCMSimulation.class.getName());
            pw.println("\n"
                    + "========================================================"
                    + "============================================\n"
                    + "\t\tGLOBAL PARAMS (may overlap w/ Neuron/Synapse Params)\n"
                    + "======================================================="
                    + "=============================================\n");
            try {
                for (Field f : RunBCMSimulation.class.getDeclaredFields()) {
                    pw.println("\t" + f.getName() + ": "
                            + f.get(new RunBCMSimulation()));
                }
                pw.println();
                pw.println("\n"
                        + "===================================================="
                        + "================================================\n"
                        + "\t\tNEURON PARAMS\n"
                        + "===================================================="
                        + "================================================\n");
                IPIFRule Iproto = null;
                IPIFRule Eproto = null;
                for (Neuron n : LAIP_RES.getNeuronList()) {
                    if (n.getPolarity() == Polarity.EXCITATORY) {
                        Eproto = (IPIFRule) n.getUpdateRule();
                    }
                    if (n.getPolarity() == Polarity.INHIBITORY) {
                        Iproto = (IPIFRule) n.getUpdateRule();
                    }
                    if (Iproto != null && Eproto != null) {
                        break;
                    }
                }
                pw.println("\n[ EXCITATORY ]\n");
                Eproto.reportAllValsToFile(pw);
                pw.println("\n[ INHIBITORY ]\n");
                Iproto.reportAllValsToFile(pw);
                pw.println();
                
                pw.println("\n"
                        + "===================================================="
                        + "================================================\n"
                        + "\t\tSYNAPSE PARAMS\n"
                        + "===================================================="
                        + "================================================\n");
                pw.println("\n[ EXCITATORY ]\n");
                pw.println(E_SPK_RESP.getClass().getName());
                E_SYN_PLASTICITY.reportAllValsToFile(pw);
                pw.println();
                pw.println("\n[ INHIBITORY-EXCITATORY ]\n");
                pw.println(I_SPK_RESP.getClass().getName());
                pw.println();
                IE_SYN_PLASTICITY.reportAllValsToFile(pw);
                pw.println("\n[ INHIBITORY-INHIBITORY ]\n");
                pw.println(I_SPK_RESP.getClass().getName());
                pw.println();
                II_SYN_PLASTICITY.reportAllValsToFile(pw);
                pw.println();
                
                pw.println("\n"
                        + "===================================================="
                        + "================================================\n"
                        + "\t\tPRUNER PARAMS\n"
                        + "===================================================="
                        + "================================================\n");
                PRUNER.reportAllValsToFile(pw);
                pw.println();
                
            } catch (IllegalAccessException iae) {
                iae.printStackTrace();
            }
            // pw.println(notes);
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
            lif.setResetPotential(RESET_POTENTIAL);
            lif.setUseFRApprox(true);

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
        LAIP_RES_SYNS.setSpikeResponder(E_SPK_RESP, Polarity.EXCITATORY);
        LAIP_RES_SYNS.setSpikeResponder(I_SPK_RESP, Polarity.INHIBITORY);
//        LAIP_RES_SYNS.setSpikeResponder(new UDF(), Polarity.BOTH);
        LAIP_RES_SYNS.setLabel("Recurrent Synapses");

        for (Synapse s : LAIP_RES_SYNS.getAllSynapses()) {
            // if (s.getTarget().getPolarity() == Polarity.INHIBITORY
            // && s.getSource().getPolarity() == Polarity.INHIBITORY) {
            // s.setDelay(0); //Gap Junctions
            // continue;
            // }
            double dist = Network
                    .getEuclideanDist(s.getSource(), s.getTarget());
            double delay = Math.round(dist / (450)) + 1;

            s.setDelay((int) delay);
            if (Math.random() < 0.01) {
                System.out.println(delay * 0.5);
            }
        }
        LAIP_RES_SYNS.setLearningRule(IE_SYN_PLASTICITY, Polarity.INHIBITORY);
        for (Neuron n : neurons) {
            if (n.getPolarity() == Polarity.INHIBITORY) {
                for (Synapse s : n.getFanOut().values()) {
                    if (s.getTarget().getPolarity() == Polarity.EXCITATORY) {
                        s.setLearningRule(IE_SYN_PLASTICITY.deepCopy());
                    } else {
                        s.setLearningRule(II_SYN_PLASTICITY.deepCopy());
                    }
                }
            }
        }
        LAIP_RES_SYNS.setLearningRule(E_SYN_PLASTICITY, Polarity.EXCITATORY);
        NETWORK.addGroup(LAIP_RES_SYNS);

        INPUT_LAYER = new NeuronGroup(NETWORK, NUM_INPUTS);
        INPUT_LAYER.setNeuronType(new SpikingThresholdRule());
        INPUT_LAYER.setTestData(getData(IN_FILE_NAME));
        INPUT_LAYER.setInputMode(true);
        INPUT_LAYER.setLabel("INPUT LAYER");
        NETWORK.addGroup(INPUT_LAYER);
        ConnectNeurons sCon;
        if (IN_DENSITY < 1.0) {
            sCon = new Sparse(IN_DENSITY, true, false);
        } else {
            sCon = new AllToAll(false);
        }
        PolarizedRandomizer inputRand = new PolarizedRandomizer(Polarity
                .EXCITATORY, IN_SYN_DIST);
        inputRand.setParam1(DIST_PARAM_1);
        inputRand.setParam2(DIST_PARAM_2);
        IN_SG = SynapseGroup.createSynapseGroup(INPUT_LAYER,
                LAIP_RES, sCon, 1.0, inputRand, inRand);
        IN_SG.setFrozen(true, Polarity.BOTH);
        ConvolvedJumpAndDecay cjd = new ConvolvedJumpAndDecay();
        IN_SG.setSpikeResponder(cjd, Polarity.BOTH);
        NETWORK.addGroup(IN_SG);
        for (Synapse s : IN_SG.getAllSynapses()) {
            s.setAuxVal(1);
        }

        LAIP_RES_SYNS.setUpperBound(100000, Polarity.EXCITATORY);
        LAIP_RES_SYNS.setLowerBound(0, Polarity.EXCITATORY);
        LAIP_RES_SYNS.setLowerBound(-100000, Polarity.INHIBITORY);
        LAIP_RES_SYNS.setUpperBound(0, Polarity.INHIBITORY);
        NETWORK.fireSynapsesUpdated();

        // Clears all the serial buffered updates
        NETWORK.getUpdateManager().clear();
        // Add the prune synapses action
        PRUNER = new PrunerAction(LAIP_RES_SYNS);
        NETWORK.getUpdateManager().addAction(PRUNER);
        NET_UPDATER = ConcurrentBufferedUpdate
                .createConcurrentBufferedUpdate(NETWORK);
        NETWORK.getUpdateManager().addAction(NET_UPDATER);

        NETWORK.updateTimeType();

        for (NeuronGroup ng : NET_UPDATER.getInputGroups()) {
            System.out.println(ng.getLabel());
        }

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
                * GRID_SPACE / 2)));
        neuron.setY(randi.nextInt((int) (Math.sqrt(NUM_NEURONS)
                * GRID_SPACE / 2)));
        neuron.setZ(randi.nextInt((int) (Math.sqrt(NUM_NEURONS)
                * GRID_SPACE * 2)));
    }

    public static void rewireSynapseGroup(SynapseGroup synGrp, int iter) {

        List<Neuron> targetNeurons = new ArrayList<Neuron>(synGrp
                .getTargetNeurons());
        List<Neuron> sourceNeurons = new ArrayList<Neuron>(synGrp
                .getSourceNeurons());
        Collections.shuffle(sourceNeurons);
        ThreadLocalRandom tlr = ThreadLocalRandom.current();
        for (int k = 0; k < iter; k++) {
            Collections.shuffle(sourceNeurons);
            for (int i = 0; i < sourceNeurons.size() - 1; i++) {
                Neuron n1 = sourceNeurons.get(i);
                Neuron n2 = sourceNeurons.get(i + 1);

                Set<Neuron> n1Targets = new HashSet<Neuron>(n1.getFanOut()
                        .keySet());
                n1Targets.retainAll(targetNeurons);
                n1Targets.remove(n2);
                Set<Neuron> n2Targets = new HashSet<Neuron>(n2.getFanOut()
                        .keySet());
                n2Targets.retainAll(targetNeurons);
                n2Targets.remove(n1);

                n2Targets.removeAll(n1.getFanOut().keySet());
                n1Targets.removeAll(n2.getFanOut().keySet());

                if (n1Targets.isEmpty() || n2Targets.isEmpty()) {
                    continue;
                }

                Neuron n3 = (Neuron) n1Targets.toArray()
                        [tlr.nextInt(n1Targets.size())];
                Neuron n4 = (Neuron) n2Targets.toArray()
                        [tlr.nextInt(n2Targets.size())];

                if (n3 == n4) {
                    continue;
                }
                Synapse toRemove2 = n2.getFanOut().get(n4);
                Synapse toRemove1 = n1.getFanOut().get(n3);

                synGrp.getParentNetwork().removeSynapse(toRemove1);
                synGrp.getParentNetwork().removeSynapse(toRemove2);
                Synapse s1 = Synapse.copySynapseWithNewSrcTarg(toRemove1, n1,
                        n4);
                Synapse s2 = Synapse.copySynapseWithNewSrcTarg(toRemove2, n2,
                        n3);
                synGrp.addSynapseUnsafe(s1);
                synGrp.addSynapseUnsafe(s2);

            }
        }
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
