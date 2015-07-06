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

import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;

import org.simbrain.network.groups.Group;
import org.simbrain.network.listeners.GroupAdapter;
import org.simbrain.network.listeners.NetworkEvent;
import org.simbrain.network.update_actions.BufferedUpdate;
import org.simbrain.network.update_actions.ConcurrentBufferedUpdate;
import org.simbrain.network.update_actions.CustomUpdate;
import org.simbrain.network.update_actions.PriorityUpdate;
import org.simbrain.network.update_actions.UpdateGroup;

/**
 * Manage network updates. Maintains a list of actions that are updated in the
 * order in which they appear in the list when the network is iterated once (in
 * the GUI, when the step button is clicked).
 *
 * @author Jeff Yoshimi
 */
public class NetworkUpdateManager {

    /**
     * The list of update actions, in a specific order. One run through these
     * actions constitutes a single "update" in the network.
     */
    private final List<NetworkUpdateAction> actionList =
            new ArrayList<NetworkUpdateAction>();

    /**
     * List of listeners on this update manager.
     */
    private List<UpdateManagerListener> listeners =
            new ArrayList<UpdateManagerListener>();

    /** Reference to parent network. */
    private final Network network;

    /**
     * Construct a new update manager.
     *
     * @param network
     */
    public NetworkUpdateManager(Network network) {
        this.network = network;
        // Default update method
        addAction(new BufferedUpdate(network));
        addListeners();

    }

    /**
     * Perform any initialization required after opening a network from xml.
     * UpdateManager will have been created from a default no argument
     * constructor ands its fields populated using xstream.
     */
    public void postUnmarshallingInit() {
        listeners = new ArrayList<UpdateManagerListener>();
        addListeners();
        Iterator<NetworkUpdateAction> actions = actionList.iterator();
        // TODO: Hack-y solution. Revisit this.
        while (actions.hasNext()) {
            NetworkUpdateAction nua = actions.next();
            if (nua instanceof ConcurrentBufferedUpdate) {
                actions.remove();
                actionList.add(ConcurrentBufferedUpdate
                        .createConcurrentBufferedUpdate(network));
                break;
            }
        }

        for (NetworkUpdateAction action : getActionList()) {
            if (action instanceof CustomUpdate) {
                ((CustomUpdate) action).init();
            }
        }
    }

    /**
     * Update manager listen for relevant changes in network. In particular
     * group update actions are added or removed as groups are added or removed.
     */
    private void addListeners() {
        network.addGroupListener(new GroupAdapter() {

            public void groupAdded(NetworkEvent<Group> e) {
                if (e.getObject().isTopLevelGroup()) {
                    addAction(new UpdateGroup(e.getObject()));
                }
            }

            public void groupRemoved(NetworkEvent<Group> e) {
                // Find corresponding group update action and remove it
                removeGroupAction(e.getObject());
            }

        });
    }

    /**
     * Returns a list of network update actions that can be added.
     *
     * @return available action list
     */
    public List<NetworkUpdateAction> getAvailableActionList() {
        final List<NetworkUpdateAction> availableActionList =
                new ArrayList<NetworkUpdateAction>();

        // By default these guys are always available
        availableActionList.add(new BufferedUpdate(network));
        availableActionList.add(new PriorityUpdate(network));
        availableActionList.add(ConcurrentBufferedUpdate
                .createConcurrentBufferedUpdate(network));

        // Add update actions for all groups available
        for (Group group : network.getGroupList()) {
            if (group.isTopLevelGroup()) {
                availableActionList.add(new UpdateGroup(group));
            }
        }

        return availableActionList;
    }

    /**
     * Remove action (if one exists) associated with the provided group.
     *
     * @param group
     *            the group being removed
     */
    public void removeGroupAction(Group group) {
        NetworkUpdateAction toDelete = null;
        for (NetworkUpdateAction action : actionList) {
            if (action instanceof UpdateGroup) {
                if (((UpdateGroup) action).getGroup() == group) {
                    toDelete = action;
                }
            }
        }
        if (toDelete != null) {
            removeAction(toDelete);
        }
    }

    /**
     * Listen for updates to the update manager.
     *
     * @param listener
     *            the listener to add
     */
    public void addListener(UpdateManagerListener listener) {
        listeners.add(listener);
    }

    /**
     * Remove listener.
     *
     * @param listener
     *            the listener to remove
     */
    public void removeListener(UpdateManagerListener listener) {
        listeners.remove(listener);
    }

    /**
     * @return the actionList
     */
    public List<NetworkUpdateAction> getActionList() {
        return actionList;
    }

    /**
     * Swap elements at the specified location.
     *
     * @param index1
     *            index of first element
     * @param index2
     *            index of second element
     */
    public void swapElements(final int index1, final int index2) {
        Collections.swap(actionList, index1, index2);
        for (UpdateManagerListener listener : listeners) {
            listener.actionOrderChanged();
        }
    }

    /**
     * Add an action to the list.
     *
     * @param action
     *            the action to add.
     */
    public void addAction(NetworkUpdateAction action) {
        actionList.add(action);
        for (UpdateManagerListener listener : listeners) {
            listener.actionAdded(action);
        }
    }

    /**
     * Completely remove an action.
     *
     * @param action
     *            the action to completely remove
     */
    public void removeAction(NetworkUpdateAction action) {
        actionList.remove(action);
        for (UpdateManagerListener listener : listeners) {
            listener.actionRemoved(action);
        }
    }

    /**
     * Listen from changes to update manager.
     */
    public interface UpdateManagerListener {

        /**
         * An action was added.
         * 
         * @param action
         *            the action to add
         */
        void actionAdded(NetworkUpdateAction action);

        /**
         * An action was removed.
         * 
         * @param action
         *            the action to remove
         */
        void actionRemoved(NetworkUpdateAction action);

        /** The action order was changed. */
        void actionOrderChanged();
    }

    /**
     * Remove all actions completely.
     */
    public void clear() {
        for (NetworkUpdateAction action : actionList) {
            for (UpdateManagerListener l : listeners) {
                l.actionRemoved(action);
            }
        }
        actionList.clear();
    }

}
