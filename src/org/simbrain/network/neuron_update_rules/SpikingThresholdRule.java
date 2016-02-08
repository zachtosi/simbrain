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
package org.simbrain.network.neuron_update_rules;

import java.util.Random;

import org.simbrain.network.core.Neuron;
import org.simbrain.network.core.SpikingNeuronUpdateRule;
import org.simbrain.network.neuron_update_rules.interfaces.NoisyUpdateRule;
import org.simbrain.util.randomizer.Randomizer;

/**
 * A simple spiking neuron that fires when weighted inputs exceed a threshold.
 * TODO: Has no documentation.
 */
public class SpikingThresholdRule extends SpikingNeuronUpdateRule implements
    NoisyUpdateRule {

    /** The noise generating randomizer. */
    private Randomizer noiseGenerator = new Randomizer();

    /** Whether or not to add noise to the inputs .*/
    private boolean addNoise;

    public SpikingThresholdRule() {
        setThreshold(.5);
    }
    
    @Override
    public SpikingThresholdRule deepCopy() {
        SpikingThresholdRule neuron = new SpikingThresholdRule();
        neuron.setThreshold(getThreshold());
        return neuron;
    }

    @Override
    public void update(Neuron neuron) {
        final double input = inputType.getInput(neuron)
                + (addNoise ? noiseGenerator.getRandom() : 0);
        if (input >= getThreshold()) {
            neuron.setSpkBuffer(true);
            setHasSpiked(true, neuron);
            neuron.setBuffer(1);
        } else {
            neuron.setSpkBuffer(false);
            setHasSpiked(false, neuron);
            neuron.setBuffer(0); // Make this a separate variable?
        }

    }

    /**
     * {@inheritDoc}
     */
    @Override
    public double getRandomValue() {
        Random rand = new Random();
        return rand.nextBoolean() ? 1 : 0;
    }

    @Override
    public String getDescription() {
        return "Spiking Threshold";
    }

    @Override
    public Randomizer getNoiseGenerator() {
        return noiseGenerator;
    }

    @Override
    public void setNoiseGenerator(Randomizer rand) {
        this.noiseGenerator = rand;
    }

    @Override
    public boolean getAddNoise() {
        return addNoise;
    }

    @Override
    public void setAddNoise(boolean noise) {
        this.addNoise = noise;
    }

}
