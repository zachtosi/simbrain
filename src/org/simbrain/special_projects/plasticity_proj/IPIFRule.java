package org.simbrain.special_projects.plasticity_proj;

import java.util.concurrent.ThreadLocalRandom;

import org.simbrain.network.core.Neuron;
import org.simbrain.network.core.SpikingNeuronUpdateRule;
import org.simbrain.network.neuron_update_rules.interfaces.NoisyUpdateRule;
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

    private double tauA = 1000;

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
        return ifn;
    }

    /**
     * {@inheritDoc}
     */
    public void update(Neuron neuron) {

        double iSyn = inputType.getInput(neuron);

        if (addNoise) {
            iSyn += ThreadLocalRandom.current().nextGaussian() * 0.5;
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
                                * (iSyn + backgroundCurrent))
                        / timeConstant;

        memPotential += dVm;

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
        neuron.setBuffer(memPotential);
    }

    /**
     * Performs the neuronal plasticity changes that affect the preferred firing
     * rate and threshold.
     * 
     * @param neuron
     * @param timeStep
     */
    private void performPlasticityChanges(Neuron neuron, double timeStep) {

        // The neuron estimates its firing rate by a decaying amount of
        // arbitrary "resource", which decays faster if the neuron preferrs to
        // fire more frequently.
        tauA = 10000 / (prefFR + 1);
        adapt += -timeStep * adapt / tauA;
        // final firing rate estimate
        estFR += timeStep * ((adapt / tauA) - estFR);
        // convert from kHz to Hz
        double estFRScale = 1000 * estFR;

        // Determine the difference between the estimate and the preference
        double deltaFR = prefFR - estFRScale;
        // Use whichever is larger (estimated or preferred firing rate)
        // as an attenuating factor.
        double fac = estFRScale > prefFR ? estFRScale + 1E-5
                : prefFR + 1E-5;

        if (usingIP && neuron.getNetwork().getTime() > 1000) {

            // Annealing procedure. ipConst (homeostatic plasticity) and
            // learningRate (intrinsic plasticity) change over time such that
            // initially neurons can easily change their preferred firing rate
            // but cannot easily change their firing threshold to maintain a
            // given preferred firing rate, but eventually the opposite is true
            // making it difficult for a neuron's preferred firing rate to
            // change but very easy for it to change its threshold so as to
            // maintain a given firing rate.
            ipConst -= timeStep * homeostaticCooling * (ipConst - ipConstFinal);
            learningRate -= timeStep * intrinsicCooling * (learningRate
                    - learningRateFinal);

            // Synaptic normalization procedure. Ensures that the sum of
            // incoming synaptic connections cannot exceed a certain value.
            // Leaves synapses undisturbed if their sum is below this value
            // and normalizes them to this value otherwise.
            double sum = 0;
            for (int i = 0, n = neuron.getFanIn().size(); i < n; i++) {
                sum += Math.abs(neuron.getFanIn().get(i).getStrength());
            }
            double saturationVal = 400 * Math.exp(0.05 * prefFR) + 1000;
            if (sum > saturationVal) {
                for (int i = 0, n = neuron.getFanIn().size(); i < n; i++) {
                    double str = neuron.getFanIn().get(i).getStrength();
                    neuron.getFanIn().get(i).forceSetStrength(saturationVal
                            * str / sum);
                }
            }

            // Alter threshold to maintain firing rate
            // homeostasis at preferred firing rate
            threshold += timeStep * threshold
                    * (Math.exp((-deltaFR) / (fac * ipConst)) - 1);

            // Adjust the preferred firing rate so as to bring it closer to
            // the preferred firing rate
            double noise = ProbDistribution.NORMAL.nextRand(0, 0.1);
            if (estFRScale > prefFR) {
                double dPi = 0;
                dPi = 2 * learningRate * Math.exp(-prefFR
                        / (beta * lowFRBoundary)) * (1 + noise);
                prefFR += dPi * timeStep;
            } else {
                double dPi = 0.0;
                if (prefFR <= lowFRBoundary) {
                    dPi = learningRate * (prefFR / lowFRBoundary) * (1 + noise);
                } else {
                    double wTerm = 1 + (Math.log(1
                            + (alpha * ((prefFR / lowFRBoundary) - 1)))
                            / alpha);
                    dPi = learningRate * wTerm * (1 + noise);
                }
                prefFR -= dPi * timeStep;
            }

            // Neurons cannot fall below a minimum preferred firing rate
            if (prefFR < minPrefFR) {
                prefFR = minPrefFR;
            }
        }
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
        return estFR;
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

}
