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
package org.simbrain.util;

import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.Insets;

import javax.swing.JComponent;
import javax.swing.JLabel;
import javax.swing.JPanel;

/**
 * <b>LabelledItemPanel</b> provides a panel for laying out labeled elements
 * neatly with all the labels and elements aligned down the screen.
 *
 * @author David Fraser
 * @author Michael Harris
 * @author Zach Tosi
 *
 */
public class LabelledItemPanel extends JPanel {

    /** The row to add the next labeled item to. */
    private int myNextItemRow = 0;

    /**
     * This method is the default constructor.
     */
    public LabelledItemPanel() {
        init();
    }

    /**
     * Initializes the panel and layout manager.
     */
    private void init() {
        setLayout(new GridBagLayout());

        // Create a blank label to use as a vertical fill so that the
        // label/item pairs are aligned to the top of the panel and are not
        // grouped in the centre if the parent component is taller than
        // the preferred size of the panel.
        GridBagConstraints constraints = new GridBagConstraints();
        constraints.gridx = 0;
        constraints.gridy = 99;
        constraints.insets = new Insets(10, 0, 0, 0);
        constraints.weighty = 1.0;
        constraints.fill = GridBagConstraints.VERTICAL;

        JLabel verticalFillLabel = new JLabel();

        add(verticalFillLabel, constraints);
    }

    /**
     * Add a labeled item to the panel. The item is added to the row below the
     * last item added.
     *
     * @param labelText The label text for the item.
     * @param item The item to be added.
     */
    public void addItem(final String labelText, final JComponent item) {
        // Create the label and its constraints
        JLabel label = new JLabel(labelText);

        GridBagConstraints labelConstraints = new GridBagConstraints();

        labelConstraints.gridx = 0;
        labelConstraints.gridy = myNextItemRow;
        labelConstraints.insets = new Insets(10, 10, 0, 0);
        labelConstraints.anchor = GridBagConstraints.NORTHEAST;
        labelConstraints.fill = GridBagConstraints.NONE;

        add(label, labelConstraints);

        // Add the component with its constraints
        GridBagConstraints itemConstraints = new GridBagConstraints();

        itemConstraints.gridx = 1;
        itemConstraints.gridy = myNextItemRow;
        itemConstraints.insets = new Insets(10, 10, 0, 10);
        itemConstraints.weightx = 1.0;
        itemConstraints.anchor = GridBagConstraints.WEST;
        itemConstraints.fill = GridBagConstraints.HORIZONTAL;

        add(item, itemConstraints);

        myNextItemRow++;
    }

    /**
     * Provides support for multi-column labeled item panels. Adds a label and
     * item to the panel on the current row, at the specified column.
     *
     * @param label The label text for the item.
     * @param item The item to be added.
     * @param col desired grid bag layer column
     */
    public void addItem(final JLabel label, final JComponent item, int col) {

        GridBagConstraints labelConstraints = new GridBagConstraints();

        labelConstraints.gridx = col;
        labelConstraints.gridy = myNextItemRow;
        labelConstraints.insets = new Insets(10, 10, 0, 0);
        labelConstraints.anchor = GridBagConstraints.NORTHEAST;
        labelConstraints.fill = GridBagConstraints.NONE;

        add(label, labelConstraints);

        // Add the component with its constraints
        GridBagConstraints itemConstraints = new GridBagConstraints();

        itemConstraints.gridx = col + 1;
        itemConstraints.gridy = myNextItemRow;
        itemConstraints.insets = new Insets(10, 10, 0, 10);
        itemConstraints.weightx = 1.0;
        itemConstraints.anchor = GridBagConstraints.WEST;
        itemConstraints.fill = GridBagConstraints.HORIZONTAL;

        add(item, itemConstraints);

    }

    /**
     * Adds a labeled item to the panel on the current myNextItemRow, at the
     * specified column.
     *
     * @param name The label text for the item.
     * @param item The item to be added.
     * @param col desired grid bag layor column
     */
    public void addItem(final String name, final JComponent item, int col) {

        JLabel label = new JLabel(name);

        GridBagConstraints labelConstraints = new GridBagConstraints();

        labelConstraints.gridx = col;
        labelConstraints.gridy = myNextItemRow;
        labelConstraints.insets = new Insets(10, 10, 0, 0);
        labelConstraints.anchor = GridBagConstraints.NORTHEAST;
        labelConstraints.fill = GridBagConstraints.NONE;

        add(label, labelConstraints);

        // Add the component with its constraints
        GridBagConstraints itemConstraints = new GridBagConstraints();

        itemConstraints.gridx = col + 1;
        itemConstraints.gridy = myNextItemRow;
        itemConstraints.insets = new Insets(10, 10, 0, 10);
        itemConstraints.weightx = 1.0;
        itemConstraints.anchor = GridBagConstraints.WEST;
        itemConstraints.fill = GridBagConstraints.HORIZONTAL;

        add(item, itemConstraints);

    }

    /**
     * Modification of addItem which takes a label, rather than text, as an
     * argument.
     *
     * @param label Label to be added
     * @param item SimbrainComponent to be added
     */
    public void addItemLabel(final JLabel label, final JComponent item) {
        GridBagConstraints labelConstraints = new GridBagConstraints();

        labelConstraints.gridx = 0;
        labelConstraints.gridy = myNextItemRow;
        labelConstraints.insets = new Insets(10, 10, 0, 0);
        labelConstraints.anchor = GridBagConstraints.NORTHEAST;
        labelConstraints.fill = GridBagConstraints.NONE;

        add(label, labelConstraints);

        // Add the component with its constraints
        GridBagConstraints itemConstraints = new GridBagConstraints();

        itemConstraints.gridx = 1;
        itemConstraints.gridy = myNextItemRow;
        itemConstraints.insets = new Insets(10, 10, 0, 10);
        itemConstraints.weightx = 1.0;
        itemConstraints.anchor = GridBagConstraints.WEST;
        itemConstraints.fill = GridBagConstraints.HORIZONTAL;

        add(item, itemConstraints);

        myNextItemRow++;
    }

    /**
     * A function which adds an item without a label. This function works with
     * the already established constraints and coordinates of the
     * LabelledItemPanel class to keep the code clean.
     *
     * @param item the desired item
     * @param col the column in which it is to be deposited
     */
    public void addItem(final JComponent item, int col) {
        GridBagConstraints itemConstraints = new GridBagConstraints();

        itemConstraints.gridx = col;
        itemConstraints.gridy = myNextItemRow;
        itemConstraints.insets = new Insets(10, 10, 0, 10);
        itemConstraints.weightx = 1.0;
        itemConstraints.anchor = GridBagConstraints.EAST;
        itemConstraints.fill = GridBagConstraints.HORIZONTAL;

        add(item, itemConstraints);
    }

    /**
     * Returns the current row.
     *
     * @return the next row where an item will be placed.
     */
    public int getMyNextItemRow() {
        return myNextItemRow;
    }

    /**
     * Sets the current row where the next item will be placed.
     *
     * @param myNextItemRow the desired row
     */
    public void setMyNextItemRow(int myNextItemRow) {
        this.myNextItemRow = myNextItemRow;
    }
}
