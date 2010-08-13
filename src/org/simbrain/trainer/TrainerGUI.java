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
package org.simbrain.trainer;

import java.awt.BorderLayout;
import java.awt.Dimension;
import java.awt.FlowLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.File;
import java.util.Iterator;

import javax.swing.BorderFactory;
import javax.swing.JButton;
import javax.swing.JComboBox;
import javax.swing.JComponent;
import javax.swing.JLabel;
import javax.swing.JMenu;
import javax.swing.JMenuBar;
import javax.swing.JMenuItem;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JSplitPane;
import javax.swing.JTable;
import javax.swing.JTextField;
import javax.swing.table.DefaultTableModel;

import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;
import org.simbrain.network.NetworkComponent;
import org.simbrain.network.builders.LayeredNetworkBuilder;
import org.simbrain.network.groups.NeuronGroup;
import org.simbrain.network.interfaces.Group;
import org.simbrain.network.interfaces.Neuron;
import org.simbrain.network.interfaces.RootNetwork;
import org.simbrain.network.listeners.GroupListener;
import org.simbrain.network.listeners.NetworkEvent;
import org.simbrain.resource.ResourceManager;
import org.simbrain.util.LabelledItemPanel;
import org.simbrain.util.SFileChooser;
import org.simbrain.util.Utils;
import org.simbrain.util.table.SimbrainDataTable;
import org.simbrain.util.table.SimbrainJTable;
import org.simbrain.workspace.Workspace;
import org.simbrain.workspace.WorkspaceComponent;
import org.simbrain.workspace.WorkspaceListener;
import org.simbrain.workspace.gui.GenericFrame;
import org.simbrain.workspace.gui.GenericJFrame;

/**
 * GUI for supervised learning in Simbrain, using back-propagation, LMS, and
 * (eventually) other algorithms. A front end for the trainer class.
 *
 * @author ericneilj
 * @author jeff yoshimi
 * @see org.simbrain.trainer.Trainer
 */
public class TrainerGUI extends JPanel {

    /** Parent frame. */
    GenericFrame parentFrame;

    /** Choices of training algorithms. */
	private String[] trainingAlgorithms = {"Backprop  ", "Least Mean Squares"};

	/** Network selection combo box. */
	private JComboBox cbNetworkChooser = new JComboBox();

	/** Input layer combo box. */
	private JComboBox cbInputLayer = new JComboBox();

	/** Table displaying input data. */
	private SimbrainJTable inputDataTable;

	/** Output layer combo box. */
    private JComboBox cbOutputLayer = new JComboBox();

    /** Left scroll pane. */
    JScrollPane leftScroll;

    /** Right scroll pane. */
    JScrollPane rightScroll;

    /** Table displaying training data. */
    private SimbrainJTable trainingDataTable;

	/** Reference to trainer object. */
	private Trainer trainer;

	/** Reference to workspace object. */
	private Workspace workspace;

	/** Current network. */
	private RootNetwork currentNetwork;

	/** Data for the error graph. */
    private XYSeries graphData;

    /** Text field for setting number of iterations to run. */
    private JTextField tfIterations;

    /** Error label. */
    private JLabel rmsError = new JLabel("Error: --- ");

    /** Update completed boolean value. */
    private boolean updateCompleted = true;

	/**
	 * Default constructor.
	 */
	public TrainerGUI(Workspace workspace, GenericFrame frame) {

	    // Initial setup
		this.workspace = workspace;
		workspace.addListener(workspaceListener);
		this.parentFrame = frame;

		// Initialize combo box action listeners
		cbNetworkChooser.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent arg0) {
                updateCurrentNetwork();
                updateLayerBoxes();
            }
		});
        cbInputLayer.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent arg0) {
                updateInputTable();
            }
        });
        cbOutputLayer.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent arg0) {
                updateOutputTable();
            }
        });

		//Main Panel
		JPanel mainPanel = new JPanel();
		mainPanel.setLayout(new BorderLayout());
		//mainPanel.setBorder(BorderFactory.createTitledBorder("Backprop properties"));

		//Top Panel
		JPanel topPanel = new JPanel();
		topPanel.setLayout(new BorderLayout());
		//topPanel.setPreferredSize(new Dimension(800, 200));


		// Top items
		JPanel  topItems = new JPanel(new FlowLayout(FlowLayout.LEFT));
	    topItems.setBorder(BorderFactory.createEtchedBorder());
	    LabelledItemPanel netSelect = new LabelledItemPanel();
        netSelect.setLayout(new FlowLayout(FlowLayout.LEFT));
        //netSelect.setPreferredSize(new Dimension(140, 220));
        JLabel netSelectLabel = new JLabel("Select Root Network");
        topItems.add(netSelectLabel);
        topItems.add(cbNetworkChooser);

        topItems.add(new JLabel("Training Algorithm"));
        JComboBox cbTrainingAlgorithm = new JComboBox(trainingAlgorithms);
        topItems.add(cbTrainingAlgorithm);
        JButton properties = new JButton(TrainerGuiActions.getPropertiesDialogAction(this));
        topItems.add(properties);
        topPanel.add("North", topItems);
    	topPanel.add("Center", getGraphPanel());

		// Split Pane (Contains two data tables)
		JSplitPane splitPane = new JSplitPane(JSplitPane.HORIZONTAL_SPLIT);
		splitPane.setBorder(null);
		splitPane.setResizeWeight(.5); //keeps divider centered on resize
		JPanel leftPanel = new JPanel();
		leftPanel.setBorder(BorderFactory.createTitledBorder("Input data"));
		JPanel rightPanel = new JPanel();
		rightPanel.setBorder(BorderFactory.createTitledBorder("Training data"));
		leftPanel.setLayout(new BorderLayout());
		rightPanel.setLayout(new BorderLayout());
		splitPane.setLeftComponent(leftPanel);
		splitPane.setRightComponent(rightPanel);

	    // Input Data
        inputDataTable = new SimbrainJTable(new SimbrainDataTable(5, 4));
        leftScroll = new JScrollPane(inputDataTable);
        leftScroll.setPreferredSize(new Dimension(100,100));
        JPanel leftMenuPanel = new JPanel();
		leftMenuPanel.setLayout(new FlowLayout(FlowLayout.LEADING));
		JLabel inputLabel = new JLabel("Input Layer:");
		leftMenuPanel.add(inputLabel);
		leftMenuPanel.add(cbInputLayer);
        leftMenuPanel.add(inputDataTable.getToolbarCSV());
        leftPanel.add("North", leftMenuPanel);
        leftPanel.add("Center", leftScroll);

        // Training Data
        trainingDataTable = new SimbrainJTable(new SimbrainDataTable(5,4));
        rightScroll = new JScrollPane(trainingDataTable);
        rightScroll.setPreferredSize(new Dimension(100,100));
		JPanel rightMenuPanel = new JPanel();
		rightMenuPanel.setLayout(new FlowLayout(FlowLayout.LEADING));
		JLabel outputLabel = new JLabel("Output Layer:");
		rightMenuPanel.add(outputLabel);
		rightMenuPanel.add(cbOutputLayer);
		rightMenuPanel.add(trainingDataTable.getToolbarCSV());
		//rightMenuPanel.add(trainingDataTable.getToolbarRandomize());
		rightPanel.add("North", rightMenuPanel);
		rightPanel.add("Center", rightScroll);

		//Bottom Button Panel
		JPanel bottomPanel = new JPanel();
		bottomPanel.add(new JButton("Cancel"));
		bottomPanel.add(new JButton("Ok"));

		// Put it all together
		mainPanel.add("North", topPanel);
		mainPanel.add("Center", splitPane);
		add(mainPanel);

		// Initialize selection box
        resetNetworkSelectionBox();

        // Initialize menus
        createMenus();
	}

	/**
	 * Create the graph panel.
	 *
	 * @return the graph panel
	 */
	private JPanel getGraphPanel() {

	    //Graph
        JPanel graphPanel = new JPanel();
        graphPanel.setBorder(BorderFactory.createTitledBorder("Error / Trainer Controls"));
        graphPanel.setLayout(new BorderLayout());
        XYSeriesCollection series = new XYSeriesCollection();
        graphData = new XYSeries(1);
        series.addSeries(graphData);
        JFreeChart chart = ChartFactory.createXYLineChart(
                null, // Title
                "Iterations", // x-axis Label
                "Error", // y-axis Label
                series, // Dataset
                PlotOrientation.VERTICAL, // Plot Orientation
                false, // Show Legend
                true, // Use tooltips
                false // Configure chart to generate URLs?
            );
        chart.setBackgroundPaint(null);
        //chart.getXYPlot().getRangeAxis().setUpperBound(1); TODO: Make autorange a dialog option
        ChartPanel centerPanel = new ChartPanel(chart);
        centerPanel.setPreferredSize(new Dimension(centerPanel
                .getPreferredSize().width, 200));

        // Make button panel
        JPanel buttonPanel = new JPanel();

        // Init
        JButton initButton = new JButton("Init");
        buttonPanel.add(initButton);
        initButton.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent arg0) {
                currentNetwork.randomizeBiases(-1, 1);
                trainer.init();
                trainer.setInputData(inputDataTable.getData().asArray());
                trainer.setTrainingData(trainingDataTable.getData().asArray());
            }
        });

        // Run
        buttonPanel.add(new JButton(TrainerGuiActions.getRunAction(this)));


        // Batch
        JButton batchButton = new JButton("Batch");
        buttonPanel.add(batchButton);
        batchButton.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent arg0) {
                if (trainer != null) {
                    trainer.train(Integer.parseInt(tfIterations.getText()));
                    rmsError.setText("Error:" +  Utils.round(trainer.getCurrentError(), 6));
                    // TODO (above): repeated code
                }
            }
        });

        // Iterations
        tfIterations = new JTextField("300");
        buttonPanel.add(new JLabel("Iterations"));
        buttonPanel.add(tfIterations);

        // Error
        buttonPanel.add(rmsError);

        // Clear Button (Not used)
        JButton clearButton = new JButton("Clear");
        clearButton.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent arg0) {
                graphData.clear();
            }
        });

        // Randomize
        buttonPanel.add(new JButton(TrainerGuiActions.getRandomizeNetworkAction(this)));


        // Finish up panel
        graphPanel.add("Center", centerPanel);
        graphPanel.add("South", buttonPanel);
        return graphPanel;
	}

	/**
	 * Create menus.
	 */
	private void createMenus() {
        JMenuBar menuBar = new JMenuBar();

        // File Menu
        JMenu fileMenu = new JMenu("File");
        fileMenu.add(new JMenuItem("Open")); // TODO
        fileMenu.add(new JMenuItem("Save")); // TODO
        fileMenu.add(new JMenuItem("Save as...")); // TODO
        menuBar.add(fileMenu);

        // File Menu
        JMenu buildMenu = new JMenu("Build");
        JMenuItem threeLayerItem = new JMenuItem(TrainerGuiActions
                .getBuildThreeLayerAction(this));
        buildMenu.add(threeLayerItem);
        JMenuItem multiLayerItem = new JMenuItem(TrainerGuiActions
                .getBuildMultiLayerAction(this));
        buildMenu.add(multiLayerItem);
        menuBar.add(buildMenu);

        parentFrame.setJMenuBar(menuBar);

	}

	/**
	 * Load data into a JTable.
	 *
	 * @param table table to laod data in to
	 */
	private void loadData(JTable table) {

        SFileChooser chooser = new SFileChooser(".", "Comma Separated Values",
                "csv");
        File theFile = chooser.showOpenDialog();
        if (theFile == null) {
            return;
        }

        // Load data in to trainer
        if (trainer != null) {
            if (table == inputDataTable) {
                trainer.setInputData(theFile);
            } else if (table == trainingDataTable) {
                trainer.setTrainingData(theFile);
            }
        }

        double[][] doubleVals = Utils.getDoubleMatrix(theFile);
        DefaultTableModel model = (DefaultTableModel) table.getModel();
        model.setNumRows(doubleVals.length);
        model.setColumnCount(doubleVals[0].length + 1);
        for (int i = 0; i < doubleVals.length; i++) {
            for (int j = 0; j < doubleVals[i].length + 1; j++) {
                if (j == 0) {
                    model.setValueAt(i + 1, i, j); // Row number column
                } else {
                    model.setValueAt(doubleVals[i][j - 1], i, j);
                }
            }
        }
    }

    /**
     * User has changed the current network in the network selection combo box.
     * Make appropriate changes.
     */
	private void updateCurrentNetwork() {

	    Object object = cbNetworkChooser.getSelectedItem();
	    if (object instanceof NetworkComponent) {
	        currentNetwork = ((NetworkComponent)object).getRootNetwork();
	    } else {
	        currentNetwork = null;
	        return;
	    }
        updateInputTable();
        updateOutputTable();

        // Initialize trainer
        if (trainer == null) {
            trainer = new BackpropTrainer(currentNetwork);
            trainer.listeners.add(new TrainerListener() {

                public void errorUpdated(double error) {
                    graphData.add(trainer.getIteration() , error);
                }

            });
        } else {
            trainer.setNetwork(currentNetwork);
        }

	    //TODO: Remove old listener
	    //previousNetwork.removeGroupListener(previousListener)
	    currentNetwork.addGroupListener(new GroupListener() {

            public void groupAdded(NetworkEvent<Group> e) {
                updateLayerBoxes();
             }

            public void groupChanged(NetworkEvent<Group> networkEvent) {
                updateLayerBoxes();
             }

            public void groupRemoved(NetworkEvent<Group> e) {
                updateLayerBoxes();
            }
	    });

	}

	/**
	 *  Reset the network selection combo box.
	 */
	private void resetNetworkSelectionBox() {

	        cbNetworkChooser.removeAllItems();
	        for (WorkspaceComponent component : workspace
                    .getComponentList(NetworkComponent.class)) {
	            cbNetworkChooser.addItem(component);
	        }
	        // TODO: This does not seem to work.  Test: Set box to network 3, add a network, it resets to 1
	        if (currentNetwork != null) {
	            cbNetworkChooser.setSelectedItem(currentNetwork);
	        } else {
	            if (cbNetworkChooser.getItemCount() >= 1)  {
	                cbNetworkChooser.setSelectedIndex(1);
	            }
	        }

	        updateInputTable();
	        updateOutputTable();
	}


    /**
     * Update input data table.
     */
    private void updateInputTable() {
        NeuronGroup group = (NeuronGroup) cbInputLayer.getSelectedItem();
        if (group != null) {
            if (trainer != null) {
                trainer.setInputLayer(group);
            }
            updateTable(group, inputDataTable);
        }
    }

    /**
     * Update training data table.
     */
    private void updateOutputTable() {
        NeuronGroup group = (NeuronGroup) cbOutputLayer.getSelectedItem();
        if (group != null) {
            if (trainer != null) {
                trainer.setOutputLayer(group);
            }
            updateTable(group, trainingDataTable);
        }
    }

    /**
     * Update the input layer and output layer combo boxes (when groups are
     * added, removed, or changed in the current network).
     */
	private void updateLayerBoxes() {

	    if (currentNetwork != null) {
	        cbInputLayer.removeAllItems();
	        for (Group group : currentNetwork.getGroupList()) {
	            if (group instanceof NeuronGroup) {
	                cbInputLayer.addItem(group);
	            }
	        }
	        cbOutputLayer.removeAllItems();
	        for (Group group : currentNetwork.getGroupList()) {
                if (group instanceof NeuronGroup) {
                    cbOutputLayer.addItem(group);
                    //TODO: Abstract this in an appropriate way
                    if (group.getLabel().equalsIgnoreCase("Output Layer")) {
                        cbOutputLayer.setSelectedItem(group);
                    }
                }
	        }
	        updateInputTable();
	        updateOutputTable();
	    }
	}

	/**
	 * Update table using the supplied neuron group.
	 */
	private void updateTable(Group group, SimbrainJTable table) {

	    int groupSize = group.getNeuronList().size();
	    int tableSize = table.getColumnCount();

	    // Make sure table and group have same number of columns
	    if (groupSize != tableSize) {
            table.getData().modifyRowsColumns(table.getData().getRowCount(),
                    groupSize, 0);
	    }

	    // Rename column headings
	    // Note the for loop starts at column 1 (column 0 is the "#" value)

	    Iterator<Neuron> neuronIterator = group.getNeuronList().iterator();
        for (int i = 1; i < tableSize; i++) { 
            if (neuronIterator.hasNext()) {
                table.getColumnModel().getColumn(i).setHeaderValue(neuronIterator.next().getDescription());
            } else {
                // Table columns not assigned to a neuron
                // TODO: Current not working properly (when reducing columns of table)
                // table.getColumnModel().getColumn(i).setHeaderValue("--");
            }
	    }
        table.getTableHeader().resizeAndRepaint();

        int rows = table.getData().getRowCount();
        int cols = table.getData().getColumnCount();
        //((JComponent)table).getParent().getParent().setPreferredSize(new Dimension(200 + cols * 100, rows *  25));
        if (((JComponent)table).getParent() != null) {
            ((JScrollPane)((JComponent)table).getParent().getParent()).revalidate();            
        }

	}

    /**
     * Listen to the workspace. When components are added update the network
     * selection combo box.
     */
    private final WorkspaceListener workspaceListener = new WorkspaceListener() {

        /**
         * Clear the Simbrain desktop.
         */
        public void workspaceCleared() {
            resetNetworkSelectionBox();
        }

        @SuppressWarnings("unchecked")
        public void componentAdded(final WorkspaceComponent workspaceComponent) {
            resetNetworkSelectionBox();
        }

        @SuppressWarnings("unchecked")
        public void componentRemoved(final WorkspaceComponent workspaceComponent) {
            if (workspaceComponent instanceof NetworkComponent) {
                if (((NetworkComponent)workspaceComponent).getRootNetwork() == currentNetwork) {
                    currentNetwork = null;
                }
            }
            resetNetworkSelectionBox();
        }
    };

    /**
     * Iterate the trainer one time and update graphics.
     */
    public void iterate() {
        trainer.train(1);
        rmsError.setText("Error:" +  Utils.round(trainer.getCurrentError(), 6));
        graphData.add(trainer.getIteration() ,trainer.getCurrentError());
    }

    /**
     * @return boolean updated completed.
     */
    public boolean isUpdateCompleted() {
        return updateCompleted;
    }

    /**
     * Sets updated completed value.
     *
     * @param updateCompleted Updated completed value to be set
     */
    public void setUpdateCompleted(final boolean updateCompleted) {
        this.updateCompleted = updateCompleted;
    }

	/**
	 * Test GUI.
	 *
	 * @param args
	 */
	public static void main(String[] args) {
	    Workspace workspace = new Workspace();

	    // Make network 1
        RootNetwork network = new RootNetwork();
        LayeredNetworkBuilder builder = new LayeredNetworkBuilder();
        int[] nodesPerLayer = new int[]{2,4,4,1};
        builder.setNodesPerLayer(nodesPerLayer);
        builder.buildNetwork(network);
        NetworkComponent networkComponent = new NetworkComponent("Net 1", network);
        workspace.addWorkspaceComponent(networkComponent);

        // Make network 2
        RootNetwork network2 = new RootNetwork();
        LayeredNetworkBuilder builder2 = new LayeredNetworkBuilder();
        int[] nodesPerLayer2 = new int[]{6, 4, 8};
        builder2.setNodesPerLayer(nodesPerLayer2);
        builder2.buildNetwork(network2);
        NetworkComponent networkComponent2 = new NetworkComponent("Net 2", network2);
        workspace.addWorkspaceComponent(networkComponent2);

        GenericJFrame topFrame = new GenericJFrame();
		TrainerGUI trainer = new TrainerGUI(workspace, topFrame);
		topFrame.setContentPane(trainer);
        topFrame.pack();
        topFrame.setVisible(true);
	}

    /**
     * @return the currentNetwork
     */
    public RootNetwork getCurrentNetwork() {
        return currentNetwork;
    }

    /**
     * @return the inputDataTable
     */
    public SimbrainJTable getInputDataTable() {
        return inputDataTable;
    }

    /**
     * @return the trainingDataTable
     */
    public SimbrainJTable getTrainingDataTable() {
        return trainingDataTable;
    }

    /**
     * @return the trainer
     */
    public Trainer getTrainer() {
        return trainer;
    }

    /**
     * @return the workspace
     */
    public Workspace getWorkspace() {
        return workspace;
    }

}
