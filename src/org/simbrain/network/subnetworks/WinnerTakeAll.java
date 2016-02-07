/*
 * Part of Simbrain--a java-based neural network kit
 * Copyright (C) 2005,2007 The Authors.  See http://www.simbrain.net/credits
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
 */
package org.simbrain.network.subnetworks;

import java.util.Random;

import org.simbrain.network.core.Network;
import org.simbrain.network.core.Neuron;
import org.simbrain.network.groups.NeuronGroup;
import org.simbrain.network.neuron_update_rules.LinearRule;

/**
 * <b>WinnerTakeAll</b>.The neuron with the highest weighted input in a
 * winner-take-all network takes on an upper value, all other neurons take on
 * the lower value. In case of a tie the node which wins is arbitrary (the first
 * in an internally maintained list).
 */
public class WinnerTakeAll extends NeuronGroup {

    /** Default initial number of units. */
    private static final int DEFAULT_NUM_UNITS = 5;

    /** Number of neurons. */
    private int numUnits = DEFAULT_NUM_UNITS;

    /** Winning value. */
    private double winValue = 1;

    /** Losing value. */
    private double loseValue = 0;

    /** If true, sometimes set the winner randomly. */
    private boolean useRandom;

    /** Probability of setting the winner randomly, when useRandom is true. */
    private double randomProb = .1;

    /**
     * Copy constructor.
     *
     * @param newRoot new root net
     * @param oldNet old network
     */
    public WinnerTakeAll(Network newRoot, WinnerTakeAll oldNet) {
        super(newRoot, oldNet);
        setLoseValue(oldNet.getLoseValue());
        setWinValue(oldNet.getWinValue());
        setUseRandom(oldNet.isUseRandom());
        setRandomProb(oldNet.getRandomProb());
        setLabel("WTA Group (copy)");
    }

    /**
     * Creates a new winner take all network.
     *
     * @param root the network containing this subnetwork
     * @param numNeurons Number of neurons in new network
     */
    public WinnerTakeAll(final Network root, final int numNeurons) {
        super(root);
        for (int i = 0; i < numNeurons; i++) {
            // TODO: Prevent invalid states like this?
            this.addNeuron(new Neuron(root, new LinearRule()));
        }
        setLabel("Winner take all network");
    }

    @Override
    public WinnerTakeAll deepCopy(Network newNetwork) {
    	return new WinnerTakeAll(newNetwork, this);
    }
    
    @Override
    public String getTypeDescription() {
        return "Winner Take All Group";
    }

    @Override
    public void update() {

        // Determine the winning neuron
        int winnerIndex;
        if (useRandom) {
            if (Math.random() < randomProb) {
                winnerIndex = getRandomWinnerIndex();
            } else {
                winnerIndex = getWinningIndex();
            }
        } else {
            winnerIndex = getWinningIndex();
        }

        // Set neuron values
        for (int i = 0; i < getNeuronList().size(); i++) {
            if (i == winnerIndex) {
                getNeuronList().get(i).setActivation(winValue);
            } else {
                getNeuronList().get(i).setActivation(loseValue);
            }
        }
    }

    /**
     *
     * Returns index of random winning neuron.
     *
     * @return index of random winner
     */
    private int getRandomWinnerIndex() {
        return new Random().nextInt(getNeuronList().size());
    }

    /**
     * Returns the index of the input node with the greatest net input.
     *
     * @return winning node's index
     */
    private int getWinningIndex() {
        int winnerIndex = 0;
        double max = Double.NEGATIVE_INFINITY;
        double lastVal = getNeuronList().get(0).getWeightedInputs();
        boolean tie = true;
        for (int i = 0; i < getNeuronList().size(); i++) {
            Neuron n = getNeuronList().get(i);
            double val = n.getWeightedInputs();
            if (val != lastVal) {
                tie = false;
            }
            lastVal = val;
            if (val > max) {
                winnerIndex = i;
                max = n.getWeightedInputs();
            }
        }
        // Break ties randomly
        // (TODO: Add a field so use can decide if they want this)
        if (tie) {
            winnerIndex = getRandomWinnerIndex();
        }
        return winnerIndex;
    }

    /**
     * @return Returns the loseValue.
     */
    public double getLoseValue() {
        return loseValue;
    }

    /**
     * @param loseValue The loseValue to set.
     */
    public void setLoseValue(final double loseValue) {
        this.loseValue = loseValue;
    }

    /**
     * @return Returns the winValue.
     */
    public double getWinValue() {
        return winValue;
    }

    /**
     * @param winValue The winValue to set.
     */
    public void setWinValue(final double winValue) {
        this.winValue = winValue;
    }

    /**
     * @return Number of neurons.
     */
    public int getNumUnits() {
        return numUnits;
    }

    /**
     * @return the useRandom
     */
    public boolean isUseRandom() {
        return useRandom;
    }

    /**
     * @param useRandom the useRandom to set
     */
    public void setUseRandom(boolean useRandom) {
        this.useRandom = useRandom;
    }

    /**
     * @return the randomProb
     */
    public double getRandomProb() {
        return randomProb;
    }

    /**
     * @param randomProb the randomProb to set
     */
    public void setRandomProb(double randomProb) {
        this.randomProb = randomProb;
    }
}
