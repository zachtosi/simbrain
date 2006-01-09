
package org.simbrain.network.actions;

import org.simbrain.network.NetworkPanel;
import org.simnet.coupling.InteractionMode;

import org.simbrain.resource.ResourceManager;

/**
 * World to network interaction mode action.
 */
public final class WorldToNetworkInteractionModeAction
    extends InteractionModeAction {

    /**
     * Create a new world to network interaction mode action with the specified
     * network panel.
     *
     * @param networkPanel network panel, must not be null
     */
    public WorldToNetworkInteractionModeAction(final NetworkPanel networkPanel) {
        super("World to network", networkPanel, InteractionMode.WORLD_TO_NETWORK);

        putValue(SMALL_ICON, ResourceManager.getImageIcon("WorldToNet.gif"));
        putValue(SHORT_DESCRIPTION, "World is sending stimuli to the network");
    }
}