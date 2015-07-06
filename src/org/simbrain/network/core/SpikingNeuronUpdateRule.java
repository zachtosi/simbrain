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
package org.simbrain.network.core;

import org.simbrain.network.core.Network.TimeType;

/**
 * <b>SpikingNeuron</b> is the superclass for spiking neuron types (e.g.
 * integrate and fire) with functions common to spiking neurons. For example a
 * boolean hasSpiked field is used in the gui to indicate that this neuron has
 * spiked.
 *
 * @author Jeff Yoshimi
 * @author Zach Tosi
 *
 */
public abstract class SpikingNeuronUpdateRule extends NeuronUpdateRule {

    {
        inputType = InputType.SYNAPTIC;
    }

    /** Time of last spike. */
    private double lastSpikeTime;
    
    private int spikeCounter = 0;
    
    private double startTime = 0;

    @Override
    public void clear(Neuron neuron) {
        super.clear(neuron);
        setLastSpikeTime(0);
    }

    /**
     * {@inheritDoc}
     */
    public TimeType getTimeType() {
        return TimeType.CONTINUOUS;
    }

    /**
     * {@inheritDoc}
     */
    public abstract void update(Neuron neuron);

    /**
     * @param hasSpiked
     *            the hasSpiked to set
     * @param neuron
     *            the neuron which has (or has not) spiked.
     */
    public void setHasSpiked(final boolean hasSpiked, final Neuron neuron) {
        if (hasSpiked) {
            lastSpikeTime = neuron.getNetwork().getTime();
            spikeCounter++;
        }
    }

    /**
     * @return the lastSpikeTime
     */
    public double getLastSpikeTime() {
        return lastSpikeTime;
    }

    public void resetCounter(final Neuron neuron) {
        spikeCounter = 0;
        startTime = neuron.getNetwork().getTime();
    }

    public int getSpikeCount() {
        return spikeCounter;
    }

    public double getFrequency(final Neuron neuron) {
        return 1000 * (spikeCounter / (neuron.getNetwork().getTime()
                - startTime));
    }

    /**
     * @param lastSpikeTime
     *            the lastSpikeTime to set
     */
    public void setLastSpikeTime(double lastSpikeTime) {
        this.lastSpikeTime = lastSpikeTime;
    }

    /**
     * A helper method which identifies this and all subclasses as variations of
     * spiking neurons. While instanceof is often bad practice this is a faster
     * way of determining if a neuron is spiking without using instanceof.
     * While normally this would still be bad practice, this is often used by
     * GUI components which are separate from the logical code.
     *
     * @return TRUE: Any subclass of SpikingNeuronUpdate rule, must by
     *         definition be a spiking neuron.
     */
    @Override
    public final boolean isSpikingNeuron() {
        return true;
    }

}
