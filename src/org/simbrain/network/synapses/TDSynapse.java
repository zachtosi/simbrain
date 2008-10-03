package org.simbrain.network.synapses;

import org.simbrain.network.interfaces.Neuron;
import org.simbrain.network.interfaces.Synapse;

public class TDSynapse extends Synapse {
    /** Learning rate. */
    public static final double DEFAULT_LEARNING_RATE = 1;
    /** reward discount factor */ 
    public static final double DEFAULT_GAMMA = 1;
    /** inputs from the last time step */
    public static final double DEFAULT_LAST_INPUT = 0;
    /** reward expectation from the last time step */
    public static final double DEFAULT_LAST_REWARD_EXPECTATION = 0;
    /** signal synapse to get the reward value */
    public static final SignalSynapse DEFAULT_REWARD_SYNAPSE = null;    
    
    /** Learning rate. */
    private double learningRate = DEFAULT_LEARNING_RATE;
    /** reward discount factor */ 
    private double gamma = DEFAULT_GAMMA;
    /** inputs from the last time step */
    private double lastInput = DEFAULT_LAST_INPUT;
    /** reward expectation from the last time step */
    private double lastRewardExpectation = DEFAULT_LAST_REWARD_EXPECTATION;
    /** signal synapse to get the reward value */
    private SignalSynapse rewardSynapse = DEFAULT_REWARD_SYNAPSE;    

    /**
     * Creates a weight of some value connecting two neurons.
     *
     * @param src source neuron
     * @param tar target neuron
     * @param val initial weight value
     * @param theId Id of the synapse
     */
    public TDSynapse(final Neuron src, final Neuron tar, final double val, final String theId) {
    	super(src, tar);
//        setSource(src);
//        setTarget(tar);
        strength = val;
        id = theId;
    }

    /**
     * This constructor is used when creating a neuron of one type from another neuron of another type Only values
     * common to different types of neuron are copied.
     * @param s Synapse to make of the type
     */
    public TDSynapse(final Synapse s) {
        super(s);
    }

    /**
     * @return Name of synapse type.
     */
    public static String getName() {
        return "TDSynapse";
    }

    /**
     * @return duplicate TDSynapse (used, e.g., in copy/paste).
     */
    public Synapse duplicate() {
        TDSynapse h = new TDSynapse(this.getSource(), this.getTarget());
        h.setLearningRate(getLearningRate());

        return super.duplicate(h);
    }

    /**
     * Creates a weight connecting source and target neurons.
     *
     * @param source source neuron
     * @param target target neuron
     */
    public TDSynapse(final Neuron source, final Neuron target) {
    	super(source, target);
//        setSource(source);
//        setTarget(target);
    }

    /**
     * Updates the strength of the synapse.
     */
    public void update() {
	
	if(rewardSynapse == null)
	    rewardSynapse = findSignalSynapse();	
	
        if(rewardSynapse != null){
            strength += learningRate * this.lastInput * 
        	(getTarget().getActivation() * gamma + rewardSynapse.getSource().getActivation() - this.lastRewardExpectation);
            strength = clip(strength);
        }
        this.lastInput = getSource().getActivation();
        this.lastRewardExpectation = getTarget().getActivation();
    }
    
    /**
     * Returns the first signal synapse discovered, null if there are none.
     *
     * @return the first signal synapse discovered, null if there are none.
     */
    private SignalSynapse findSignalSynapse() {
        SignalSynapse ret = null;
        for(Synapse synapse : this.getTarget().getFanIn()){
            if (synapse instanceof SignalSynapse) {
                return (SignalSynapse) synapse;
            }
        }
        return ret;
    }    
    
    /**
     * Clear the state of input and output neuron for
     * the last time step
     */
    public void ResetPreviousNeuronState(){
	this.lastInput = 0;
        this.lastRewardExpectation = 0;
    }    
    
    /**
     * @return Returns the momentum.
     */
    public double getLearningRate() {
        return learningRate;
    }

    /**
     * @param momentum The momentum to set.
     */
    public void setLearningRate(final double momentum) {
        this.learningRate = momentum;
    }

    /**
     * @return Returns gamma (reward discount factor)
     */
    public double getGamma() {
        return gamma;
    }

    /**
     * @param gamma The reward discount factor to set
     */
    public void setGamma(double gamma) {
        this.gamma = gamma;
    }
}