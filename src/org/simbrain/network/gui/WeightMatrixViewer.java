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
package org.simbrain.network.gui;

import java.awt.BorderLayout;
import java.awt.FlowLayout;
import java.util.ArrayList;
import java.util.List;

import javax.swing.JOptionPane;
import javax.swing.JPanel;

import org.simbrain.network.core.Neuron;
import org.simbrain.network.core.Synapse;
import org.simbrain.network.util.SimnetUtils;
import org.simbrain.util.table.NumericTable;
import org.simbrain.util.table.SimbrainJTable;
import org.simbrain.util.table.SimbrainJTableScrollPanel;

/**
 * Widget to display the synaptic connections between two layers of neurons as a
 * matrix, in a jtable.
 *
 * TODO: Better display of non-existent connections (perhaps by disabling those
 * cells for now). What would be super cool is gray for no synapse. and then as
 * added or deleted ungray it.
 *
 * @author jyoshimi
 */
public class WeightMatrixViewer extends SimbrainJTableScrollPanel {

    /** JTable contained in scroller. */
    private SimbrainJTable table;

    /**
     * Embed the scrollpanel in a widget with a toolbar.
     *
     * @param scroller the scroller to embed
     * @return the formatted jpanel.
     */
    public static JPanel getWeightMatrixPanel(WeightMatrixViewer scroller) {
        JPanel panel = new JPanel(new BorderLayout());
        panel.add("Center", scroller);
        JPanel toolbar = new JPanel(new FlowLayout(FlowLayout.LEFT));
        toolbar.add(scroller.getTable().getToolbarRandomize());
        toolbar.add(scroller.getTable().getToolbarCSV(false, false));
        panel.add("North", toolbar);
        return panel;
    }

    /**
     * Construct a weight matrix viewer using a specified list of source and
     * target neurons.
     *
     * @param sourceList the source neurons
     * @param targetList the target neurons
     * @param panel the parent network panel
     */
    public WeightMatrixViewer(List<Neuron> sourceList, List<Neuron> targetList,
            NetworkPanel panel) {
        init(sourceList, targetList, panel);
    }

    /**
     * Create a panel for viewing the matrices connecting a set of source and
     * target neuron lists.
     *
     * @param panel the panel from which to draw the matrix.
     */
    public WeightMatrixViewer(NetworkPanel panel) {

        // Get source and target lists
        ArrayList<Neuron> sourceList = panel.getSourceModelNeurons();
        ArrayList<Neuron> targetList = panel.getSelectedModelNeurons();
        init(sourceList, targetList, panel);
    }

    /**
     * Initialize the weight matrix viewer.
     *
     * @param sourceList the source neurons
     * @param targetList the target neurons
     * @param panel the network panel
     */
    private void init(List<Neuron> sourceList, List<Neuron> targetList,
            NetworkPanel panel) {
        // By default the lists are sorted horizontally.
        // TODO: Allow for vertical sorting, or for some appropriate sorting
        // when displaying an adjacency matrix

        // By default the lists are sorted horizontally.
        // TODO: Allow for vertical sorting, or for some appropriate sorting
        // when displaying an adjacency matrix
        // Collections.sort(sourceList, OrientationComparator.X_ORDER);
        // Collections.sort(targetList, OrientationComparator.X_ORDER);

        // Populate data in simbrain table
        Synapse[][] weights = SimnetUtils.getWeightMatrix(sourceList,
                targetList);
        //displayWarningIfEmptyCells(weights);
        WeightMatrix weightMatrix = new WeightMatrix(weights);
        table = SimbrainJTable.createTable(weightMatrix);
        table.disableTableModificationMenus();

        // Create names for row headings
        List<String> rowHeaders = new ArrayList<String>();
        for (Neuron neuron : sourceList) {
            rowHeaders.add(new String(neuron.getId()));
        }

        // Create names for column headings
        List<String> colHeaders = new ArrayList<String>();
        for (Neuron neuron : targetList) {
            colHeaders.add(new String(neuron.getId()));
        }
        table.setColumnHeadings(colHeaders);
        table.setRowHeadings(rowHeaders);
        table.getData().fireTableStructureChanged();

        // Set the table
        this.setTable(table);

    }

    /**
     * Display a warning message if there are empty weights.
     *
     * @param weights weight matrix to check
     */
    private void displayWarningIfEmptyCells(Synapse[][] weights) {
        String warningMessage = "Only fully connected source-target pairs \n"
                + "are supported.  Some zeros in the matrix \n"
                + "correspond to non-existent weights and \n"
                + "cannot be modified in the viewer.";
        for (int i = 0; i < weights.length; i++) {
            for (int j = 0; j < weights[0].length; j++) {
                if (weights[i][j] == null) {
                    JOptionPane.showMessageDialog(null, warningMessage,
                            "Weight Matrix Error", JOptionPane.WARNING_MESSAGE);
                    return;
                }
            }
        }
    }

    /**
     * Matrix of synapses to be viewed in a SimbrainJTable.
     *
     * A matrix representation of synapses is passed in, and as the table data
     * are changed the synapses are directly updated.
     *
     * Note that the "mutable" features of numerictable are all passed over.
     *
     */
    private class WeightMatrix extends NumericTable {

        /** Underlying data. */
        private Synapse[][] weights;

        /**
         * @param weights the weights to set
         */
        public WeightMatrix(Synapse[][] weights) {
            super(weights.length, weights[0].length);
            this.weights = weights;
        }

        // TOOD: Explain below. Possibly move to  a special immutable numeric class
        @Override
        public void setValue(final int row, final int col, final Double value,
                final boolean fireEvent) {
            if (col == 0) {
                return;
            }
            if (weights[row][col-1] != null) {
                weights[row][col-1].forceSetStrength(value);
            }
            if (fireEvent) {
                fireTableDataChanged();
            }
        }

        @Override
        public void setValue(int row, int col, Double value) {
            setValue(row, col, value, true);
        }


        @Override
        public Double getValueAt(int row, int col) {
            if (col == 0) {
                return null;
            }
            if (weights[row][col - 1] != null) {
                return weights[row][col - 1].getStrength();
            } else {
                return null;
            }
        }


        @Override
        public void setLogicalValue(int row, int column, Double value,
                boolean fireEvent) {
            setValue(row, column+1, value, fireEvent);
        }


        @Override
        public Double getLogicalValueAt(int row, int col) {
            return getValueAt(row, col+1);
        }



    }

    /**
     * @return the table
     */
    public SimbrainJTable getTable() {
        return table;
    }

}
