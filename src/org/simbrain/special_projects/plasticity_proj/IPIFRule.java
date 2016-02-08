package org.simbrain.special_projects.plasticity_proj;

import java.io.PrintWriter;
import java.lang.reflect.Field;
import java.util.concurrent.ThreadLocalRandom;

import org.simbrain.network.core.Neuron;
import org.simbrain.network.core.SpikingNeuronUpdateRule;
import org.simbrain.network.neuron_update_rules.interfaces.NoisyUpdateRule;
import org.simbrain.util.SimbrainConstants.Polarity;
import org.simbrain.util.math.ProbDistribution;
import org.simbrain.util.randomizer.Randomizer;

public class IPIFRule extends SpikingNeuronUpdateRule implements
        NoisyUpdateRule {

    /** Resistance (M ohms). */
    private double resistance = 1;

    /** Time constant (ms) */
    private double timeConstant = 30;

    /** Threshold (mV) */
    private double threshold = 15;

    /** Reset potential (mV) */
    private double resetPotential = 10;

    /** Resting potential (mV) Default: 0.0 */
    private double restingPotential;

    /** Background Current (nA) . */
    private double backgroundCurrent = 13.5;

    /** Noise dialog. */
    private Randomizer noiseGenerator = new Randomizer();

    /** Add noise to neuron. */
    private boolean addNoise;

    /**
     * @deprecated here for backwards compatibility.
     */
    @Deprecated
    private boolean hasSpiked;

    private boolean usingIP = true;

    private double ipConst = 1E7;

    /** Hz */
    private double prefFR = 1;

    /** Hz */
    private double lowFRBoundary = 5;

    private double beta = 60;

    private double tauMinus = 100;

    private double tauPlus = 100;

    private double learningRate = .1;

    private double alpha = 1;

    private double refractoryPeriod = 1.0;

    private double ipConstFinal = 1E5;

    private double learningRateFinal = 1E-7;

    private double intrinsicCooling = 1E-4;

    private double homeostaticCooling = 1E-4;

    private double minPrefFR = 1E-4; // Hz or Spk/s

    private double estFR = 0;

    private double adapt = 0;

    private double tauA = 10000;

 //   private native float gaussConv(float x);

//    static {
//        File nativeFile = new File("libipifrule.jnilib");
//        if (!nativeFile.exists())
//            System.exit(1);
//        System.load(nativeFile.getAbsolutePath());
//    }

    private float approxFR = 0;

    private float meanApproxFR;

    private int count = 0;

    /**
     * {@inheritDoc}
     */
    public IPIFRule deepCopy() {
        IPIFRule ifn = new IPIFRule();
        ifn.setRestingPotential(getRestingPotential());
        ifn.setResetPotential(getResetPotential());
        ifn.setThreshold(getThreshold());
        ifn.setBackgroundCurrent(getBackgroundCurrent());
        ifn.setTimeConstant(getTimeConstant());
        ifn.setResistance(getResistance());
        ifn.setIncrement(getIncrement());
        ifn.setAddNoise(getAddNoise());
        ifn.setBeta(getBeta());
        ifn.setAlpha(getAlpha());
        ifn.setUsingIP(isUsingIP());
        ifn.setLowFRBoundary(getLowFRBoundary());
        ifn.setTauMinus(getTauMinus());
        ifn.setTauPlus(getTauPlus());
        ifn.setLearningRate(getLearningRate());
        ifn.setRefractoryPeriod(getRefractoryPeriod());
        ifn.noiseGenerator = new Randomizer(noiseGenerator);
        ifn.setPrefFR(getPrefFR());
        ifn.setIpConst(getIpConst());
        ifn.setTauA(tauA);
        ifn.setIpConstFinal(ipConstFinal);
        ifn.setLearningRateFinal(learningRateFinal);
        ifn.setMinPrefFR(minPrefFR);
        ifn.setIntrinsicCooling(intrinsicCooling);
        ifn.setHomeostaticCooling(homeostaticCooling);
        ifn.setUseFRApprox(useFRApprox);
        ifn.ss = ss;
        ifn.s = s;
        ifn.m = m;
        ifn.l = l;
        ifn.setFullNorm(fullNormEx);
        ifn.satA = satA;
        ifn.satB = satB;
        ifn.adaptiveTauA = adaptiveTauA;
        return ifn;
        
    }

    double adaptPar = .1;
    
    boolean init = false;
    
    /**
     * {@inheritDoc}
     */
    public void update(Neuron neuron) {
        // if (ThreadLocalRandom.current().nextDouble() < 0.1) {
        // return;
        // }
        
    	if (!init) {
    		if (neuron.getPolarity() == Polarity.EXCITATORY) {
    			adaptPar = .1;
    			satC = 2;
    		} else {
    			alpha = 10;
    			beta = 40;
    			adaptPar = .05;
    			satC = 0;
    			lowFRBoundary = 4;
    		}
    		init = true;
    	}
    	
        double iSyn = inputType.getInput(neuron);
        
        if (addNoise) {
            iSyn += ThreadLocalRandom.current().nextGaussian() * 0.05;
        }
        double timeStep = neuron.getNetwork().getTimeStep();

        double memPotential = neuron.getActivation();

        // Update the neuron's preferred firing rate and firing threshold
        // based on homeostatic and intrinsic plasticity
        performPlasticityChanges(neuron, timeStep);

        // Perform the update of the membrane potential.
        /*
         * Formula:
         * 
         * dV/dt = ( -(Vm - Vr) + Rm * (Isyn + Ibg) ) / tau
         * 
         * Vm > theta ? Vm <- Vreset ; spike
         * 
         * Vm: membrane potential Vr: resting potential* Rm: membrane resistance
         * Isyn: synaptic input current Ibg: background input current tau: time
         * constant Vreset: reset potential theta: threshold
         */
        double dVm =
                timeStep
                        * (-(memPotential - restingPotential) + resistance
                                * (iSyn + backgroundCurrent - (adaptPar * adapt)))
                        / timeConstant;

        memPotential += dVm;
        
        //if (neuron.getNetwork().getTime() < 50000) {
        //    if (ThreadLocalRandom.current().nextInt((int)(1000/neuron
        //            .getNetwork().getTimeStep())) == 1) {
        //        memPotential = threshold + 1;
        //    }
        //}

        // Determine if there was a spike
        if (memPotential >= threshold
                && neuron.getNetwork().getTime()
                > getLastSpikeTime() + refractoryPeriod)
        {
            // Send a spike
            neuron.setSpkBuffer(true);
            setHasSpiked(true, neuron);
            memPotential = resetPotential;
            // Update the firing rate estimator
            adapt += 1;
        } else {
            // Do not send a spike.
            neuron.setSpkBuffer(false);
            setHasSpiked(false, neuron);
        }
//        if (useFRApprox) {
//            double iSyn_theta = threshold - restingPotential
//                    - backgroundCurrent + (0.1 * adapt);
//            afrBuff = gaussConv((float) (iSyn - iSyn_theta));
//            ss_buff += 0.5f * (approxFR- ss);
//            s_buff += 0.5f * (ss - s);
//            m_buff += 0.1f * (s - m);
//            mafrBuff = (0.9f * s) + (0.1f * m);
//            
//
////          
//            if (meanApproxFR > l) {
//                l_buff += 0.6 * (meanApproxFR - l);
//            } else {
//                l_buff += 0.05 * (meanApproxFR - l);
//            }
//
////            if (loc_index == 1 && neuron.getNetwork().getTime() > 1000) {
////                SimbrainConstants.Debug_stream.print(approxFR + " ");
////                SimbrainConstants.Debug_stream.print(ss + " ");
////                SimbrainConstants.Debug_stream.print(s + " ");
////                SimbrainConstants.Debug_stream.print(m + " ");
////                SimbrainConstants.Debug_stream.print(meanApproxFR + " ");
////                SimbrainConstants.Debug_stream.print(l + " ");
////                SimbrainConstants.Debug_stream.print(getEstFR() + " ");
////                SimbrainConstants.Debug_stream.print(prefFR+ " ");
////                SimbrainConstants.Debug_stream.println(this.getSpikeCount());
////            }
//
////            mafrBuff = ((meanApproxFR * count) + (approxFR * approxFR))
////                    / ++count;
//        }
        neuron.setBuffer(memPotential);
    }

    float ss = 1;
    float s = 1;
    float m = 1;
    private float l = 1;
    float ss_buff = 1;
    float s_buff = 1;
    float m_buff = 1;
    private float l_buff = 1;
    
    
    public float getMafrLong() {
        return l;
    }
    
    private boolean useFRApprox = false;

    static int index = 0;
    
    int loc_index;
    
    public IPIFRule() {
        super();
        loc_index = index++;
    }
    
    public void setUseFRApprox(boolean useFRApprox) {
        this.useFRApprox = useFRApprox;
    }

    private float afrBuff = 0;

    private float mafrBuff = 0;

    public void pushBuffers() {
        this.approxFR = afrBuff;
        this.meanApproxFR = mafrBuff;
        ss = ss_buff;
        s = s_buff;
        m = m_buff;
        l = l_buff;
    }

    private boolean annealing = true;

    private double dPi = 0;

    private boolean fullNormEx = false;
    private boolean fullNormIn = false;
    
    private float satA = 5; // 10
    private float satB = 5; // 20
    private float satC = 2;
    
    private boolean adaptiveTauA = false;
    

    /**
     * Performs the neuronal plasticity changes that affect the preferred firing
     * rate and threshold.
     * 
     * @param neuron
     * @param timeStep
     */
    private void performPlasticityChanges(final Neuron neuron,
            final double timeStep) {

        // The neuron estimates its firing rate by a decaying amount of
        // arbitrary "resource", which decays faster if the neuron preferrs to
        // fire more frequently.
        //if (adaptiveTauA) {
          //  tauA = 10000 / (prefFR + 1);
        //}
        double pfr = prefFR > 50 ? 50 : prefFR;
    	if (neuron.getPolarity() == Polarity.INHIBITORY) {
    		tauA = 20000 / (pfr + 2);
    	}
        adapt += -timeStep * adapt / tauA;
        // final firing rate estimate
        estFR += timeStep * ((adapt / tauA) - estFR);
        // convert from kHz to Hz
        double estFRScale = 1000 * estFR;



        if (usingIP && neuron.getNetwork().getTime() > 1000) {
            // Synaptic normalization procedure. Ensures that the sum of
            // incoming synaptic connections cannot exceed a certain value.
            // Leaves synapses undisturbed if their sum is below this value
            // and normalizes them to this value otherwise.
            if (isAnnealing()) {
                if (usingIP && neuron.getNetwork().getTime() > 10000) {
                    double inExSum = 0;
                    double inInSum = 0;
                    // double outSum = 0;
                    int ex = 0;
                    int in = 0;
                    for (int i = 0, n = neuron.getFanIn().size(); i < n; i++) {
                    	if (neuron.getFanIn().get(i).getAuxVal() == 0) {
                    		if( neuron.getFanIn().get(i).getStrength() > 0) {
                    			inExSum += neuron.getFanIn().get(i).getStrength();
                    			ex++;
                    		}
                    		if( neuron.getFanIn().get(i).getStrength() < 0) {
                    			inInSum -= neuron.getFanIn().get(i).getStrength();
                    			in++;
                    		}
                    	}
                    }

                    // double saturationVal = (prefFR * 150) + 200;
//                    double saturationValEx = ((13.5/threshold)*((prefFR+10) / (estFRScale+10))  * (prefFR * satA))
//                            + satB;
                    //double saturationValEx = 75*Math.log((prefFR+4.25)/4.25)*((prefFR+10)/(estFRScale+10))+20;
                    //double saturationValEx = (5 * prefFR * (15/threshold) * (estFRScale+20)/(prefFR+20)) + 50; 

                    double saturationValEx = (threshold/15) * ((pfr*(satA-satC))+estFRScale*satC) + satB;
                    double sValEx = 20 * ex;
                    sValEx = (saturationValEx > sValEx) ? sValEx : saturationValEx;
                    if ((inExSum > sValEx || fullNormEx)) {
                    	if (!fullNormEx) {
                    		System.out.println(neuron.getPolarity().title() + " EXC\tTripped!");
                    	}
                    	fullNormEx = true;
                    	for (int i = 0, n = neuron.getFanIn().size(); i < n; i++) {
                    		if (neuron.getFanIn().get(i).getAuxVal() == 0) {
                    			if (neuron.getFanIn().get(i).getStrength() > 0) {
                    				double str = neuron.getFanIn().get(i).getStrength();
                    				double newstr = sValEx * str / inExSum;
                    				neuron.getFanIn().get(i).changeStrength(newstr - str);
                    			}
                    		}
                        }
                    }
                    
//                    double saturationValIn = ((threshold/13.5)*((estFRScale+10) / (prefFR+10))  * (prefFR * satA))
//                            + satB;
                    //double saturationValIn = (.5*Math.exp(.04*prefFR) + .75) * saturationValEx;
                    double saturationValIn = (15/threshold) * ((pfr*(satA-satC))+estFRScale*satC) + satB;
                    double sValIn = 20 * in;
                    sValIn = (saturationValIn > sValIn) ? sValIn : saturationValIn;
                    if ((inInSum > sValIn || fullNormIn)) {
                    	if (!fullNormIn) {
                    		System.out.println(neuron.getPolarity().title() + " INH\tTripped!");
                    	}
                    	fullNormIn = true;
                    	for (int i = 0, n = neuron.getFanIn().size(); i < n; i++) {
                    		if (neuron.getFanIn().get(i).getAuxVal() == 0) {
                    			if (neuron.getFanIn().get(i).getStrength() < 0) {
                    				double str = neuron.getFanIn().get(i).getStrength();
                    				double newstr = sValIn * str / inInSum;
                    				neuron.getFanIn().get(i).changeStrength(newstr - str);
                    			}
                    		}
                        }
                    }
                }
            }

            // if (prefFR > 1) {
            // // double outSat = (-3900 * Math.exp(-0.45 * prefFR))
            // // + (3500 * Math.exp(-0.041 * prefFR));
            // double outSat = saturationVal;
            // if (outSum > outSat) {
            // for (Synapse s : neuron.getFanOut().values()) {
            // s.forceSetStrength(outSat * s.getStrength()
            // / outSum);
            // }
            // }
            // }

            // Annealing procedure. ipConst (homeostatic plasticity) and
            // learningRate (intrinsic plasticity) change over time such that
            // initially neurons can easily change their preferred firing rate
            // but cannot easily change their firing threshold to maintain a
            // given preferred firing rate, but eventually the opposite is true
            // making it difficult for a neuron's preferred firing rate to
            // change but very easy for it to change its threshold so as to
            // maintain a given firing rate.
            if (isAnnealing()) {
            	ipConst -= timeStep * homeostaticCooling
            			* (ipConst - ipConstFinal);
            	learningRate -= timeStep * intrinsicCooling * (learningRate
            			- learningRateFinal);
            }

            // Determine the difference between the estimate and the preference
            double deltaFR = prefFR - estFRScale;
            // Use whichever is larger (estimated or preferred firing rate)
            // as an attenuating factor.
            //double fac = estFRScale > prefFR ? estFRScale + 1E-5
           // 		: prefFR + 1E-5;

            double thVal = 1;
            // Adjust the preferred firing rate so as to bring it closer to
            // the preferred firing rate
            double noise = ProbDistribution.NORMAL.nextRand(0, 0.05);
            if (estFRScale > prefFR) {
            	dPi = ((learningRate * Math.exp(-prefFR
            			/ (beta * lowFRBoundary)))) // + (0.9 * dPi))
            			* (1 + noise);
            	prefFR += dPi * timeStep * invLogitDer((prefFR/estFRScale) - 1);
            	//thVal = 2;
            } else {
            	if (prefFR <= lowFRBoundary) {
            		dPi = ((-learningRate * (prefFR / lowFRBoundary)))
            				// + (0.9 * dPi))
            				* (1 + noise);

            		//thVal = 0.5;
            	} else {
            		double wTerm = 1 + (Math.log(1
            				+ (alpha * ((prefFR / lowFRBoundary) - 1)))
            				/ alpha);
            		dPi = ((-learningRate * wTerm))// + (0.9 * dPi))
            				* (1 + noise);
            		//thVal = 1;
            	}
            	prefFR += dPi * timeStep * invLogitDer((prefFR/estFRScale) - 1);
            }
            if(prefFR < .001){
            	prefFR = .001;
            }
            // Alter threshold to maintain firing rate
            // homeostasis at preferred firing rate
            threshold += thVal * timeStep * threshold
            		* (Math.exp((-deltaFR) / (prefFR * ipConst)) - 1);
        }
    }
    
    private double invLogitDer(double a) {
    	//return 1;
    	if (Math.abs(a) > 2) return 1;
    	return -4*(Math.exp(a*4)/(Math.pow(Math.exp(a*4) + 1, 2) ) ) + 1;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public double getRandomValue() {
        // Equal chance of spiking or not spiking, taking on any value between
        // the resting potential and the threshold if not.
        return 2 * (threshold - restingPotential) * Math.random()
                + restingPotential;
    }

    /**
     * @return Returns the lowerValue.
     */
    public double getRestingPotential() {
        return restingPotential;
    }

    /**
     * @param restingPotential
     *            The restingPotential to set.
     */
    public void setRestingPotential(final double restingPotential) {
        this.restingPotential = restingPotential;
    }

    /**
     * @return Returns the upperValue.
     */
    public double getResistance() {
        return resistance;
    }

    /**
     * @param resistance
     *            The resistance to set.
     */
    public void setResistance(final double resistance) {
        this.resistance = resistance;
    }

    /**
     * @return Returns the lowerValue.
     */
    public boolean getAddNoise() {
        return addNoise;
    }

    /**
     * @param addNoise
     *            The addNoise to set.
     */
    public void setAddNoise(final boolean addNoise) {
        this.addNoise = addNoise;
    }

    /**
     * @param noise
     *            The noise to set.
     */
    public void setAddNoise(final Randomizer noise) {
        this.noiseGenerator = noise;
    }

    /**
     * @return Returns the noiseGenerator.
     */
    public Randomizer getNoiseGenerator() {
        return noiseGenerator;
    }

    /**
     * @param noiseGenerator
     *            The noiseGenerator to set.
     */
    public void setNoiseGenerator(final Randomizer noiseGenerator) {
        this.noiseGenerator = noiseGenerator;
    }

    /**
     * @return Returns the resetPotential.
     */
    public double getResetPotential() {
        return resetPotential;
    }

    /**
     * @param resetPotential
     *            The resetPotential to set.
     */
    public void setResetPotential(final double resetPotential) {
        this.resetPotential = resetPotential;
    }

    /**
     * @return Returns the background current
     */
    public double getBackgroundCurrent() {
        return backgroundCurrent;
    }

    /**
     * @param backgroundCurrent
     *            The background current to set
     */
    public void setBackgroundCurrent(double backgroundCurrent) {
        this.backgroundCurrent = backgroundCurrent;
    }

    /**
     * @return Returns the threshold.
     */
    public double getThreshold() {
        return threshold;
    }

    /**
     * @param threshold
     *            The threshold to set.
     */
    public void setThreshold(final double threshold) {
        this.threshold = threshold;
    }

    /**
     * @return Returns the timeConstant.
     */
    public double getTimeConstant() {
        return timeConstant;
    }

    /**
     * @param timeConstant
     *            The timeConstant to set.
     */
    public void setTimeConstant(final double timeConstant) {
        this.timeConstant = timeConstant;
    }

    @Override
    public String getDescription() {
        return "Integrate and Fire";
    }

    @Override
    public double getGraphicalLowerBound() {
        return resetPotential - 30;
    }

    @Override
    public double getGraphicalUpperBound() {
        return restingPotential + 30;
    }

    public boolean isUsingIP() {
        return usingIP;
    }

    public void setUsingIP(boolean usingIP) {
        this.usingIP = usingIP;
    }

    public double getIpConst() {
        return ipConst;
    }

    public void setIpConst(double ipConst) {
        this.ipConst = ipConst;
    }

    public double getPrefFR() {
        return prefFR;
    }

    public void setPrefFR(double prefFR) {
        this.prefFR = prefFR;
    }

    public double getLowFRBoundary() {
        return lowFRBoundary;
    }

    public void setLowFRBoundary(double lowFRBoundary) {
        this.lowFRBoundary = lowFRBoundary;
    }

    public double getTauMinus() {
        return tauMinus;
    }

    public void setTauMinus(double tauMinus) {
        this.tauMinus = tauMinus;
    }

    public double getTauPlus() {
        return tauPlus;
    }

    public void setTauPlus(double tauPlus) {
        this.tauPlus = tauPlus;
    }

    public double getLearningRate() {
        return learningRate;
    }

    public void setLearningRate(double learningRate) {
        this.learningRate = learningRate;
    }

    public double getRefractoryPeriod() {
        return refractoryPeriod;
    }

    public void setRefractoryPeriod(double refractoryPeriod) {
        this.refractoryPeriod = refractoryPeriod;
    }

    public double getBeta() {
        return beta;
    }

    public void setBeta(double beta) {
        this.beta = beta;
    }

    public double getAlpha() {
        return alpha;
    }

    public void setAlpha(double alpha) {
        this.alpha = alpha;
    }

    public double getEstFR() {
        return 1000.0 * estFR;
    }

    public void setEstFR(double estFR) {
        this.estFR = estFR;
    }

    public double getAdapt() {
        return adapt;
    }

    public void setAdapt(double adapt) {
        this.adapt = adapt;
    }

    public double getTauA() {
        return tauA;
    }

    public void setTauA(double tauA) {
        this.tauA = tauA;
    }

    public double getIpConstFinal() {
        return ipConstFinal;
    }

    public void setIpConstFinal(double ipConstFinal) {
        this.ipConstFinal = ipConstFinal;
    }

    public double getLearningRateFinal() {
        return learningRateFinal;
    }

    public void setLearningRateFinal(double learningRateFinal) {
        this.learningRateFinal = learningRateFinal;
    }

    public double getIntrinsicCooling() {
        return intrinsicCooling;
    }

    public void setIntrinsicCooling(double intrinsicCooling) {
        this.intrinsicCooling = intrinsicCooling;
    }

    public double getHomeostaticCooling() {
        return homeostaticCooling;
    }

    public void setHomeostaticCooling(double homeostaticCooling) {
        this.homeostaticCooling = homeostaticCooling;
    }

    public double getMinPrefFR() {
        return minPrefFR;
    }

    public void setMinPrefFR(double minPrefFR) {
        this.minPrefFR = minPrefFR;
    }

    public double getMaxFR() {
        return 1000 / refractoryPeriod;
    }

    public boolean isAnnealing() {
        return annealing;
    }

    public void setAnnealing(boolean annealing) {
        this.annealing = annealing;
    }

    public boolean isFullNorm() {
        return fullNormEx;
    }

    public void setFullNorm(boolean fullNorm) {
        this.fullNormEx = fullNorm;
    }

    public float getApproxFR() {
        return approxFR;
    }

    public void setApproxFR(float approxFR) {
        this.approxFR = approxFR;
    }

    public float getMeanApproxFR() {
        return meanApproxFR;
    }

    public void setMeanApproxFR(float meanApproxFR) {
        this.meanApproxFR = meanApproxFR;
    }
    
    public void reportAllValsToFile(PrintWriter pw)
            throws IllegalAccessException {
        pw.println("\t\t**[" + this.toString() + "]**");
        for (Field f : IPIFRule.class.getDeclaredFields()) {
            pw.println("\t" + f.getName() + ": " + f.get(this));
        }
    }

    public String toString() {
        return "IPIFRule";
    }
    
}
