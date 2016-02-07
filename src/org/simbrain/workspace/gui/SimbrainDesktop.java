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
package org.simbrain.workspace.gui;

import java.awt.BorderLayout;
import java.awt.Dimension;
import java.awt.Point;
import java.awt.Rectangle;
import java.awt.Toolkit;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.ComponentAdapter;
import java.awt.event.ComponentEvent;
import java.awt.event.ComponentListener;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import java.awt.event.WindowListener;
import java.beans.PropertyVetoException;
import java.io.File;
import java.lang.reflect.Constructor;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.Map;
import java.util.Properties;
import java.util.Stack;
import java.util.Vector;

import javax.swing.AbstractAction;
import javax.swing.Action;
import javax.swing.BorderFactory;
import javax.swing.JButton;
import javax.swing.JDesktopPane;
import javax.swing.JFrame;
import javax.swing.JInternalFrame;
import javax.swing.JLabel;
import javax.swing.JMenu;
import javax.swing.JMenuBar;
import javax.swing.JMenuItem;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JPopupMenu;
import javax.swing.JSplitPane;
import javax.swing.JTabbedPane;
import javax.swing.JToolBar;
import javax.swing.SwingUtilities;
import javax.swing.WindowConstants;
import javax.swing.event.InternalFrameAdapter;
import javax.swing.event.InternalFrameEvent;
import javax.swing.event.InternalFrameListener;
import javax.swing.event.MenuEvent;
import javax.swing.event.MenuListener;

import org.apache.log4j.Logger;
import org.simbrain.console.ConsoleComponent;
import org.simbrain.console.ConsoleDesktopComponent;
import org.simbrain.docviewer.DocViewerComponent;
import org.simbrain.docviewer.DocViewerDesktopComponent;
import org.simbrain.network.NetworkComponent;
import org.simbrain.network.desktop.NetworkDesktopComponent;
import org.simbrain.plot.barchart.BarChartComponent;
import org.simbrain.plot.barchart.BarChartGui;
import org.simbrain.plot.histogram.HistogramComponent;
import org.simbrain.plot.histogram.HistogramGui;
import org.simbrain.plot.piechart.PieChartComponent;
import org.simbrain.plot.piechart.PieChartGui;
import org.simbrain.plot.projection.ProjectionComponent;
import org.simbrain.plot.projection.ProjectionGui;
import org.simbrain.plot.rasterchart.RasterPlotComponent;
import org.simbrain.plot.rasterchart.RasterPlotGui;
import org.simbrain.plot.scatterplot.ScatterPlotComponent;
import org.simbrain.plot.scatterplot.ScatterPlotGui;
import org.simbrain.plot.timeseries.TimeSeriesPlotComponent;
import org.simbrain.plot.timeseries.TimeSeriesPlotGui;
import org.simbrain.resource.ResourceManager;
import org.simbrain.util.SFileChooser;
import org.simbrain.util.StandardDialog;
import org.simbrain.util.Utils;
import org.simbrain.util.genericframe.GenericFrame;
import org.simbrain.util.genericframe.GenericJFrame;
import org.simbrain.util.genericframe.GenericJInternalFrame;
import org.simbrain.util.widgets.ShowHelpAction;
import org.simbrain.util.widgets.ToggleButton;
import org.simbrain.workspace.Coupling;
import org.simbrain.workspace.Workspace;
import org.simbrain.workspace.WorkspaceComponent;
import org.simbrain.workspace.WorkspaceListener;
import org.simbrain.workspace.WorkspaceSerializer;
import org.simbrain.workspace.updater.InterceptingEventQueue;
import org.simbrain.workspace.updater.WorkspaceUpdaterListener;
import org.simbrain.world.dataworld.DataWorldComponent;
import org.simbrain.world.dataworld.DataWorldDesktopComponent;
import org.simbrain.world.game.GameComponent;
import org.simbrain.world.game.GameDesktopComponent;
import org.simbrain.world.odorworld.OdorWorldComponent;
import org.simbrain.world.odorworld.OdorWorldDesktopComponent;
import org.simbrain.world.textworld.DisplayComponent;
import org.simbrain.world.textworld.DisplayComponentDesktopGui;
import org.simbrain.world.textworld.ReaderComponent;
import org.simbrain.world.textworld.ReaderComponentDesktopGui;
import org.simbrain.world.visionworld.VisionWorldComponent;
import org.simbrain.world.visionworld.VisionWorldDesktopComponent;

import bsh.Interpreter;
import bsh.util.JConsole;

/**
 * Creates a Swing-based environment for working with a workspace.
 *
 * Also provides wrappers for GUI elements called from a terminal.
 *
 * @author Matt Watson
 * @author Jeff Yoshimi
 */
public class SimbrainDesktop {

    /** The x offset for popup menus. */
    private static final int MENU_X_OFFSET = 5;

    /** The y offset for popup menus. */
    private static final int MENU_Y_OFFSET = 53;

    /** The default serial version ID. */
    private static final long serialVersionUID = 1L;

    /** Log4j logger. */
    private static final Logger LOGGER = Logger.getLogger(Workspace.class);

    /** Initial indent of entire workspace. */
    private static final int WORKSPACE_INSET = 80;

    /** After placing one simbrain window how far away to put the next one. */
    private static final int DEFAULT_WINDOW_OFFSET = 30;

    /**
     * Reference to the last internal frames that were focused, so that they can
     * get the focus when the next one is closed.
     */
    private static final Stack<GuiComponent<?>> lastFocusedStack = new Stack<GuiComponent<?>>();

    /** TODO: Create Javadoc comment. */
    private static final Map<Workspace, SimbrainDesktop> INSTANCES = new HashMap<Workspace, SimbrainDesktop>();

    // TODO: Review. Part of a hack solution in NeuronGroupNode
    // and SynapseGroup dialog, useful anyway?
    public static Map<Workspace, SimbrainDesktop> getInstances() {
		return INSTANCES;
	}

	/** Desktop pane. */
    private JDesktopPane desktop;

    /** Cached context menu. */
    private JPopupMenu contextMenu;

    /** Workspace tool bar. */
    private JToolBar wsToolBar = new JToolBar();

    /** Whether the bottom dock is visible. */
    private boolean dockVisible = true;

    /** the frame that will hold the workspace. */
    private JFrame frame;

    /** The bottom dock. */
    private JTabbedPane bottomDock;

    /** Pane splitter for bottom dock. */
    private JSplitPane horizontalSplitter;

    /** The workspace this desktop wraps. */
    private final Workspace workspace;

    /** Boundary of workspace. */
    private Rectangle workspaceBounds;

    /** Workspace action manager. */
    private WorkspaceActionManager actionManager;
    
    /** Menu for re-centering windows. */
    private static JMenu recenterMenu;

    /** Interpreter for terminal. */
    Interpreter interpreter;

    /** Time indicator. */
    private JLabel timeLabel = new JLabel();

    /** "Throbber" to indicate a simulation is running. */
    private JLabel runningLabel = new JLabel();

    /** Associates workspace components with their corresponding gui components. */
    private Map<WorkspaceComponent, GuiComponent<?>> guiComponents = new LinkedHashMap<WorkspaceComponent, GuiComponent<?>>();

    /** Listener on the workspace. */
    private final WorkspaceListener workspaceListener = new WorkspaceListener() {

        /**
         * Clear the Simbrain desktop.
         */
        public void workspaceCleared() {
            desktop.removeAll();
            desktop.repaint();
            frame.setTitle("Simbrain");
            updateTimeLabel();
        }

        /**
         * Add a new <c>SimbrainComponent</c>.
         */
        @SuppressWarnings("unchecked")
        public void componentAdded(final WorkspaceComponent workspaceComponent) {
            addDesktopComponent(workspaceComponent);
        }

        @SuppressWarnings("unchecked")
        public void componentRemoved(final WorkspaceComponent workspaceComponent) {
            GuiComponent<?> component = guiComponents.get(workspaceComponent);
            guiComponents.remove(component);
            component.getParentFrame().dispose();
            if (!lastFocusedStack.isEmpty()) {
                lastFocusedStack.remove(component);
            }
            moveLastFocusedComponentToFront();
        }

        /**
         * {@inheritDoc}
         */
        public void newWorkspaceOpened() {
            frame.setTitle(workspace.getCurrentFile().getName());
            updateTimeLabel();
        }

    };

    /** Listens for workspace updater events. */
    private final WorkspaceUpdaterListener updaterListener = new WorkspaceUpdaterListener() {

        /**
         * {@inheritDoc}
         */
        public void changeNumThreads() {
        }

        /**
         * {@inheritDoc}
         */
        public void changedUpdateController() {
        }

        /**
         * {@inheritDoc}
         */
        public void updatedCouplings(int update) {
        }

        /**
         * {@inheritDoc}
         */
        public void updatingStarted() {
            StandardDialog.setSimulationRunning(true);
        }

        /**
         * {@inheritDoc}
         */
        public void updatingFinished() {
            StandardDialog.setSimulationRunning(false);
        }

        /**
         * {@inheritDoc}
         */
        public void workspaceUpdated() {
            updateTimeLabel();
        }
    };

    // TODO this should be addressed at a higher level
    public static SimbrainDesktop getDesktop(final Workspace workspace) {
        return INSTANCES.get(workspace);
    }

    /**
     * Default constructor.
     *
     * @param workspace The workspace for this desktop.
     */
    public SimbrainDesktop(final Workspace workspace) {

        INSTANCES.put(workspace, this);
        this.workspace = workspace;
        frame = new JFrame("Simbrain");
        actionManager = new WorkspaceActionManager(this);
        createAndAttachMenus();
        wsToolBar = createToolBar();
        createContextMenu();
        workspace.addListener(workspaceListener);
        workspace.getUpdater().addUpdaterListener(updaterListener);
        SimbrainDesktop.registerComponents();
        Dimension screenSize = Toolkit.getDefaultToolkit().getScreenSize();
        workspaceBounds = new Rectangle(WORKSPACE_INSET, WORKSPACE_INSET,
                screenSize.width - (WORKSPACE_INSET * 2), screenSize.height
                        - (WORKSPACE_INSET * 2));

        // Set the bottom dock to visible or not based on the properties file.
        Properties properties = Utils.getSimbrainProperties();
        if (properties.containsKey("showBottomDock")) {
            dockVisible = Boolean.parseBoolean(properties
                    .getProperty("showBottomDock"));
        }

        // Set up Desktop
        desktop = new JDesktopPane();
        desktop.addMouseListener(mouseListener);
        desktop.addKeyListener(new WorkspaceKeyAdapter(workspace));
        desktop.setPreferredSize(new Dimension(screenSize.width
                - (WORKSPACE_INSET * 2), screenSize.height
                - (WORKSPACE_INSET * 3)));

        // Create the Tabbed Pane for bottom of the desktop
        bottomDock = new JTabbedPane();
        bottomDock.addTab("Components", null, new ComponentPanel(this),
                "Show workspace components");

        // List of current couplings for populating couplings panel.
        Vector<Coupling<?>> couplings = new Vector<Coupling<?>>(workspace
                .getCouplingManager().getCouplings());
        bottomDock.addTab("Couplings", null, new CouplingListPanel(this,
                couplings), "Show current couplings");
        bottomDock.addTab("Terminal", null, this.getTerminalPanel(),
                "Simbrain terminal");
        bottomDock.addTab("Updater", null,
                new ThreadViewerPanel(this.getWorkspace()),
                "Simbrain thread viewer");
        // Set up the main panel
        horizontalSplitter = new JSplitPane(JSplitPane.VERTICAL_SPLIT);
        horizontalSplitter.setDividerLocation(getDividerLocation());
        horizontalSplitter.setTopComponent(desktop);
        horizontalSplitter.setBottomComponent(bottomDock);
        JPanel mainPanel = new JPanel(new BorderLayout());
        mainPanel.add(wsToolBar, "North");
        mainPanel.add(horizontalSplitter, "Center");
        if (dockVisible == false) {
            horizontalSplitter.getBottomComponent().setVisible(false);
        }

        // Set up Frame
        frame.setBounds(workspaceBounds);
        frame.setContentPane(mainPanel);
        frame.pack();
        frame.addWindowListener(windowListener);
        frame.addKeyListener(new WorkspaceKeyAdapter(workspace));

        // Start terminal
        new Thread(interpreter).start();

        // Make dragging a little faster but perhaps uglier.
        // desktop.setDragMode(JDesktopPane.OUTLINE_DRAG_MODE);
    }

    /**
     * Create mappings from guiComponents to their GUI wrappers.
     */
    private static void registerComponents() {
        // TODO use a configuration file
        registerComponent(BarChartComponent.class, BarChartGui.class);
        registerComponent(DocViewerComponent.class,
                DocViewerDesktopComponent.class);
        registerComponent(ConsoleComponent.class, ConsoleDesktopComponent.class);
        registerComponent(DisplayComponent.class,
                DisplayComponentDesktopGui.class);
        registerComponent(DataWorldComponent.class,
                DataWorldDesktopComponent.class);
        // registerComponent(MidiWorldComponent.class,
        // MidiWorldDesktopComponent.class);
        registerComponent(HistogramComponent.class, HistogramGui.class);
        registerComponent(NetworkComponent.class, NetworkDesktopComponent.class);
        registerComponent(OdorWorldComponent.class,
                OdorWorldDesktopComponent.class);
        registerComponent(PieChartComponent.class, PieChartGui.class);
        registerComponent(ProjectionComponent.class, ProjectionGui.class);
        registerComponent(ReaderComponent.class,
                ReaderComponentDesktopGui.class);
        registerComponent(ScatterPlotComponent.class, ScatterPlotGui.class);
        registerComponent(TimeSeriesPlotComponent.class,
                TimeSeriesPlotGui.class);
        registerComponent(RasterPlotComponent.class, RasterPlotGui.class);
        registerComponent(VisionWorldComponent.class,
                VisionWorldDesktopComponent.class);
        registerComponent(GameComponent.class, GameDesktopComponent.class);
    }

    /** Listener for swing component changes. */
    private final ComponentListener componentListener = new ComponentAdapter() {
        /**
         * Responds to component moved events.
         *
         * @param event SimbrainComponent event
         */
        public void componentMoved(final ComponentEvent event) {

            // Prevent window from being moved outside of visible area
            int x = (int) event.getComponent().getBounds().getX();
            int y = (int) event.getComponent().getBounds().getY();
            int width = (int) event.getComponent().getBounds().getWidth();
            int height = (int) event.getComponent().getBounds().getHeight();
            if (x < desktop.getVisibleRect().getX()) {
                event.getComponent().setBounds(0, y, width, height);
            }
            if (y < desktop.getVisibleRect().getY()) {
                event.getComponent().setBounds(x, 0, width, height);
            }

            // Workspace has changed
            workspace.setWorkspaceChanged(true);
        }

        /**
         * Responds to component resized events.
         *
         * @param arg0 SimbrainComponent event
         */
        public void componentResized(final ComponentEvent arg0) {
            // System.out.println("Component resized");
            workspace.setWorkspaceChanged(true);
        }
    };

    /**
     * Takes the last gui component opened and moves it to the front of the
     * simbrain desktop, place it in focus.
     */
    private void moveLastFocusedComponentToFront() {
        if (!lastFocusedStack.isEmpty()) {
            GuiComponent<?> lastFocused = lastFocusedStack.peek();
            if (lastFocused != null) {
                try {
                    ((JInternalFrame) lastFocused.getParentFrame())
                            .setSelected(true);
                } catch (Exception e) {
                    e.printStackTrace();
                }
            }
        }
    }

    /**
     * @return Terminal panel.
     */
    private JConsole getTerminalPanel() {
        JConsole console = new JConsole();
        interpreter = ConsoleDesktopComponent.getSimbrainInterpreter(console,
                this.getWorkspace());
        try {
            interpreter.set("desktop", this);
        } catch (Exception e) {
            e.printStackTrace();
        }
        console.setPreferredSize(new Dimension(400, 300));
        return console;
    }

    /**
     * Print text to terminal.
     *
     * @param toPrint text to print
     */
    public void printToTerminal(final String toPrint) {
        interpreter.println(toPrint);
    }

    /**
     * Returns the workspace.
     *
     * @return the workspace.
     */
    public Workspace getWorkspace() {
        return workspace;
    }

    /**
     * Returns the main frame for the desktop.
     *
     * @return the main frame for the desktop.
     */
    public JFrame getFrame() {
        return frame;
    }

    /**
     * Creates the workspace tool bar.
     *
     * @return JToolBar tool bar created
     */
    private JToolBar createToolBar() {
        JToolBar bar = new JToolBar();

        bar.add(actionManager.getOpenWorkspaceAction());
        bar.add(actionManager.getSaveWorkspaceAction());
        bar.addSeparator();
        bar.add(actionManager.getGlobalUpdateAction());
        bar.add(new ToggleButton(actionManager.getGlobalControlActions()));

        bar.addSeparator();
        bar.add(actionManager.getOpenCouplingManagerAction());

        bar.addSeparator();
        bar.add(actionManager.getNewNetworkAction());

        /* World menu button. */
        JButton button = new JButton();
        button.setIcon(ResourceManager.getImageIcon("World.png"));
        final JPopupMenu worldMenu = new JPopupMenu();
        for (Action action : actionManager.getNewWorldActions()) {
            worldMenu.add(action);
        }
        button.addActionListener(new ActionListener() {
            public void actionPerformed(final ActionEvent e) {
                JButton button = (JButton) e.getSource();
                worldMenu.show(button, 0, button.getHeight());
            }
        });
        button.setComponentPopupMenu(worldMenu);
        bar.add(button);

        /* Chart menu button. */
        button = new JButton();
        button.setIcon(ResourceManager.getImageIcon("BarChart.png"));
        final JPopupMenu gaugeMenu = new JPopupMenu();
        for (Action action : actionManager.getPlotActions()) {
            gaugeMenu.add(action);
        }
        button.addActionListener(new ActionListener() {
            public void actionPerformed(final ActionEvent e) {
                JButton button = (JButton) e.getSource();
                gaugeMenu.show(button, 0, button.getHeight());
            }
        });
        button.setComponentPopupMenu(gaugeMenu);
        bar.add(button);
        bar.add(actionManager.getNewConsoleAction());

        // Initialize time label
        timeLabel.setBorder(BorderFactory.createEmptyBorder(0, 10, 0, 10));
        timeLabel.addMouseListener(new MouseAdapter() {
            // Reset time if user double clicks on label.
            public void mousePressed(final MouseEvent event) {
                if (event.getClickCount() == 2) {
                    workspace.getUpdater().resetTime();
                    updateTimeLabel();
                }
            }
        });
        runningLabel.setIcon(ResourceManager.getImageIcon("Throbber.gif"));
        runningLabel.setVisible(false);
        updateTimeLabel();
        bar.add(timeLabel);
        bar.add(runningLabel);

        return bar;
    }

    /**
     * Create and attach workspace menus.
     */
    private void createAndAttachMenus() {
        JMenuBar menuBar = new JMenuBar();
        menuBar.add(createFileMenu());
        menuBar.add(createViewMenu());
        menuBar.add(createInsertMenu());
        menuBar.add(createScriptMenu());
        menuBar.add(createCoupleMenu());
        menuBar.add(createHelpMenu());
        frame.setJMenuBar(menuBar);
    }

    /**
     * Create script menu.
     *
     * @return script JMenu
     */
    private JMenu createScriptMenu() {
        JMenu scriptMenu = new JMenu("Scripts");
        // scriptMenu.add(actionManager.getRunScriptAction());
        scriptMenu.add(actionManager.getShowScriptEditorAction());
        scriptMenu.addSeparator();
        scriptMenu.addMenuListener(menuListener);
        for (Action action : actionManager.getScriptActions(this)) {
            scriptMenu.add(action);
        }
        return scriptMenu;
    }

    /**
     * Create the workspace file menu.
     *
     * @return file menu
     */
    private JMenu createFileMenu() {
        JMenu fileMenu = new JMenu("File");
        fileMenu.addMenuListener(menuListener);
        for (Action action : actionManager.getOpenSaveWorkspaceActions()) {
            fileMenu.add(action);
        }
        fileMenu.addSeparator();
        fileMenu.add(actionManager.getClearWorkspaceAction());
        fileMenu.addSeparator();
        fileMenu.add(actionManager.getOpenNetworkAction());

        JMenu worldSubMenu = new JMenu("Open World");
        for (Action action : actionManager.getOpenWorldActions()) {
            worldSubMenu.add(action);
        }
        fileMenu.add(worldSubMenu);
        fileMenu.addSeparator();
        fileMenu.add(actionManager.getShowUpdaterDialog());
        fileMenu.add(actionManager.getShowPropertyDialogAction());
        fileMenu.addSeparator();
        fileMenu.add(actionManager.getQuitWorkspaceAction());
        return fileMenu;
    }

    /**
     * Create the workspace view menu.
     *
     * @return view menu
     */
    private JMenu createViewMenu() {
        JMenu viewMenu = new JMenu("View");
        viewMenu.add(actionManager.getPropertyTabAction());
        viewMenu.addSeparator();
        initializeRecenterWindowMenu();
        viewMenu.add(recenterMenu);
        return viewMenu;
    }
    
    /**
     * Create a menu for re-centering windows in case they get "lost" off-screen.
     */
    private void initializeRecenterWindowMenu() {
        recenterMenu = new JMenu("Recenter window");
        final Workspace workspace = this.getWorkspace();
        for (WorkspaceComponent component : workspace.getComponentList()) {
            JMenuItem item = createRecenterWindowMenuItem(component.getName());
            item.setName(component.getName());
            recenterMenu.add(item);
        }
    }  

    /** Create new Menu item to be added to re-center-window menu */
    private JMenuItem createRecenterWindowMenuItem(String title) {
        return new JMenuItem(new AbstractAction(title) {
            @Override
            public void actionPerformed(ActionEvent e) {
                //find which window we need to bring to front
                for (WorkspaceComponent component : workspace.getComponentList()) {
                    if (component.getName().equals(getValue(AbstractAction.NAME))) {
                        GuiComponent c = guiComponents.get(component);
                        // Re-centers windows to top left at 10,10
                        c.getParentFrame().setLocation(10, 10);
                        c.getParentFrame().toFront();
                    }
                }
            }
        });
    }

    /**
     * Create the workspace insert menu.
     *
     * @return insert menu
     */
    private JMenu createInsertMenu() {
        JMenu insertMenu = new JMenu("Insert");
        insertMenu.add(actionManager.getNewNetworkAction());
        // insertMenu.add(new OpenEditorAction(this)); //TODO: Move this action
        // manager
        JMenu newGaugeSubMenu = new JMenu("New Plot");
        for (Action action : actionManager.getPlotActions()) {
            newGaugeSubMenu.add(action);
        }
        insertMenu.add(newGaugeSubMenu);
        JMenu newWorldSubMenu = new JMenu("New World");
        for (Action action : actionManager.getNewWorldActions()) {
            newWorldSubMenu.add(action);
        }
        insertMenu.add(newWorldSubMenu);
        insertMenu.addSeparator();
        insertMenu.add(actionManager.getNewDocViewerAction());
        insertMenu.add(actionManager.getNewConsoleAction());
        return insertMenu;
    }

    /**
     * Create the workspace couplings menu.
     *
     * @return couplings menu
     */
    private JMenu createCoupleMenu() {
        JMenu coupleMenu = new JMenu("Couplings");
        coupleMenu.add(actionManager.getOpenCouplingManagerAction());
        coupleMenu.add(actionManager.getOpenCouplingListAction());
        return coupleMenu;
    }

    /**
     * Create the workspace help menu.
     *
     * @return help menu
     */
    private JMenu createHelpMenu() {
        JMenu helpMenu = new JMenu("Help");
        helpMenu.add(new ShowHelpAction("Main Help", "SimbrainDocs.html"));
        helpMenu.addSeparator();
        helpMenu.add(new ShowHelpAction("Quick start", "Pages/QuickStart.html"));
        helpMenu.add(new ShowHelpAction("Keyboard Shortcuts",
                "KeyboardShortcuts.html"));
        helpMenu.add(new ShowHelpAction("Credits", "SimbrainCredits.html"));
        return helpMenu;
    }

    /**
     * Create a new context menu for this network panel.
     */
    private void createContextMenu() {
        contextMenu = new JPopupMenu();
        contextMenu.add(actionManager.getNewNetworkAction());
        JMenu newGaugeSubMenu = new JMenu("New Plot");
        for (Action action : actionManager.getPlotActions()) {
            newGaugeSubMenu.add(action);
        }
        contextMenu.add(newGaugeSubMenu);
        JMenu newWorldSubMenu = new JMenu("New World");
        for (Action action : actionManager.getNewWorldActions()) {
            newWorldSubMenu.add(action);
        }
        contextMenu.add(newWorldSubMenu);
        contextMenu.addSeparator();
        contextMenu.add(actionManager.getNewDocViewerAction());
        contextMenu.add(actionManager.getNewConsoleAction());

    }

    /**
     * This nasty declaration creates a map of the workspace guiComponents to
     * their associated wrapper class.
     */
    private static final Map<Class<? extends WorkspaceComponent>,
            Class<? extends GuiComponent<?>>> wrappers =
            new HashMap<Class<? extends WorkspaceComponent>,
            Class<? extends GuiComponent<?>>>();

    /**
     * Registers a gui wrapper class association with a component class.
     *
     * @param component
     *            The component class.
     * @param gui
     *            The gui class.
     */
    private static void registerComponent(
            final Class<? extends WorkspaceComponent> component,
            final Class<? extends GuiComponent<?>> gui) {
        wrappers.put(component, gui);
    }

    /**
     * Returns the desktop component corresponding to a workspace component.
     *
     * @param component component to check with
     * @return component guicomponent
     */
    public GuiComponent<?> getDesktopComponent(
            final WorkspaceComponent component) {
        return guiComponents.get(component);
    }

    /**
     * Returns the desktop component corresponding to a named workspace
     * component.
     *
     * @param componentName name of desktop component to return
     * @return component desktop component, or null if none found
     */
    public GuiComponent<?> getDesktopComponent(final String componentName) {
        WorkspaceComponent wc = workspace.getComponent(componentName);
        if (wc != null) {
            return guiComponents.get(wc);
        } else {
            return null;
        }
    }

    /**
     * Utility class for adding internal frames, which are not wrappers for
     * WorkspaceComponents. Wraps GUI Component in a JInternalFrame for Desktop.
     */
    private static class DesktopInternalFrame extends GenericJInternalFrame {

        /** Reference to workspace component. */
        private WorkspaceComponent workspaceComponent;

        /** Gui Component. */
        private GuiComponent guiComponent;

        /**
         * Construct an internal frame.
         *
         * @param workspaceComponent workspace component.
         */
        public DesktopInternalFrame(final WorkspaceComponent workspaceComponent) {
            init();
            this.workspaceComponent = workspaceComponent;
        }

        /**
         * Initialize the frame.
         */
        private void init() {
            setResizable(true);
            setMaximizable(true);
            setIconifiable(true);
            setClosable(true);
            setDefaultCloseOperation(WindowConstants.DO_NOTHING_ON_CLOSE);
            addInternalFrameListener(new WindowFrameListener());
        }

        /**
         * Set the Gui Component.
         *
         * @param guiComponent the component to set.
         */
        public void setGuiComponent(final GuiComponent guiComponent) {
            this.guiComponent = guiComponent;
        }

        /**
         * Manage cleanup when a component is closed.
         */
        private class WindowFrameListener extends InternalFrameAdapter {
            @Override
            public void internalFrameActivated(final InternalFrameEvent e) {
                // TODO: Does not work properly. Should be used so that
                // the last focused stack tracks changes in focus and not just
                // open / close events.
                // lastFocusedStack.remove(guiComponent);
                // lastFocusedStack.push(guiComponent);
            }

            @Override
            public void internalFrameOpened(InternalFrameEvent e) {
                super.internalFrameOpened(e);

                String title = workspaceComponent.getName();
                if (recenterMenu != null) {
                    JMenuItem item =  guiComponent
                            .getDesktop().createRecenterWindowMenuItem(title);
                    item.setName(title);
                    recenterMenu.add(item);
                }
            }

            @Override
            public void internalFrameClosing(final InternalFrameEvent e) {
                if (workspaceComponent.hasChangedSinceLastSave()) {
                    boolean hasCancelled = guiComponent.showHasChangedDialog();
                    if (hasCancelled) {
                        return;
                    }
                }
                guiComponent.close();
            }

            @Override
            public void internalFrameClosed(InternalFrameEvent e) {
                super.internalFrameClosed(e);
                if (recenterMenu != null) {
                    for (int i = 0; i < recenterMenu.getItemCount(); i++) {
                        if (recenterMenu.getItem(i).getName()
                                .equals(e.getInternalFrame().getTitle())) {
                            recenterMenu.remove(i);
                            break;
                        }
                    }
                }
            }

        }
    } // End DesktopInternalFrame

    /**
     * Add internal frame.
     *
     * @param internalFrame the frame to add.
     */
    public void addInternalFrame(final JInternalFrame internalFrame) {
        internalFrame.addInternalFrameListener(new InternalFrameListener() {

            @Override
            public void internalFrameActivated(InternalFrameEvent arg0) {
                // TODO Auto-generated method stub

            }

            @Override
            public void internalFrameClosed(InternalFrameEvent arg0) {
            }

            @Override
            public void internalFrameClosing(InternalFrameEvent arg0) {
                moveLastFocusedComponentToFront();
            }

            @Override
            public void internalFrameDeactivated(InternalFrameEvent arg0) {
                // TODO Auto-generated method stub

            }

            @Override
            public void internalFrameDeiconified(InternalFrameEvent arg0) {
                // TODO Auto-generated method stub

            }

            @Override
            public void internalFrameIconified(InternalFrameEvent arg0) {
                // TODO Auto-generated method stub

            }

            @Override
            public void internalFrameOpened(InternalFrameEvent arg0) {
                // TODO Auto-generated method stub

            }

        });
        desktop.add(internalFrame);
    }

    /**
     * Registers instance of guiComponents.
     *
     * @param workspaceComponent Workspace component
     * @param guiComponent GUI component
     */
    public void registerComponentInstance(
            final WorkspaceComponent workspaceComponent,
            final GuiComponent guiComponent) {
        guiComponent.setDesktop(this);
        guiComponents.put(workspaceComponent, guiComponent);
    }

    /**
     * Add a new SimbrainComponent.
     *
     * @param workspaceComponent Workspace Component
     */
    @SuppressWarnings("unchecked")
    public void addDesktopComponent(final WorkspaceComponent workspaceComponent) {
        LOGGER.trace("Adding workspace component: " + workspaceComponent);

        final DesktopInternalFrame componentFrame = new DesktopInternalFrame(
                workspaceComponent);
        GuiComponent guiComponent = createDesktopComponent(componentFrame,
                workspaceComponent);
        componentFrame.setGuiComponent(guiComponent);

        // Either add the window at a default location, or relative to the last
        // added window. Note that this is overridden when individual
        // components are opened
        if (guiComponents.size() == 0) {
            componentFrame.setBounds(DEFAULT_WINDOW_OFFSET,
                    DEFAULT_WINDOW_OFFSET, (int) guiComponent
                            .getPreferredSize().getWidth(), (int) guiComponent
                            .getPreferredSize().getHeight());
        } else {
            int highestComponentNumber = workspace.getComponentList().size();
            componentFrame
                    .setBounds(
                            (int) ((highestComponentNumber * DEFAULT_WINDOW_OFFSET) % (desktop
                                    .getWidth() - guiComponent
                                    .getPreferredSize().getWidth())),
                            (int) ((highestComponentNumber * DEFAULT_WINDOW_OFFSET) % (desktop
                                    .getHeight() - guiComponent
                                    .getPreferredSize().getHeight())),
                            (int) guiComponent.getPreferredSize().getWidth(),
                            (int) guiComponent.getPreferredSize().getHeight());
        }

        // Other initialization
        componentFrame.addComponentListener(componentListener);
        componentFrame.setContentPane(guiComponent);
        registerComponentInstance(workspaceComponent, guiComponent);
        componentFrame.setVisible(true);
        componentFrame.setTitle(workspaceComponent.getName());
        desktop.add(componentFrame);
        guiComponent.postAddInit();
        lastFocusedStack.push(guiComponent);
        // System.out.println(lastOpened.getName());

        // Forces last component of the desktop to the front
        try {
            ((JInternalFrame) componentFrame).setSelected(true);
        } catch (PropertyVetoException e) {
            e.printStackTrace();
        }

    }

    /**
     * Creates an instance of the proper wrapper class around the provided
     * instance.
     *
     * @param component The component to wrap.
     * @param parentFrame The frame of this component
     * @return A new desktop component wrapping the provided component.
     */
    @SuppressWarnings("unchecked")
    static GuiComponent<?> createDesktopComponent(
            final GenericFrame parentFrame, final WorkspaceComponent component) {
        Class<? extends WorkspaceComponent> componentClass = component
                .getClass();
        Class<? extends GuiComponent<?>> guiClass = wrappers
                .get(componentClass);

        if (guiClass == null) {
            throw new IllegalArgumentException(
                    "no desktop component registered for "
                            + component.getClass());
        }

        try {
            GenericFrame genericFrame = parentFrame != null ? parentFrame
                    : new DesktopInternalFrame(component);

            Constructor<? extends GuiComponent<?>> constructor = guiClass
                    .getConstructor(GenericFrame.class, componentClass);
            return constructor.newInstance(genericFrame, component);
        } catch (RuntimeException e) {
            throw e;
        } catch (Exception e) {
            throw new RuntimeException(e);
        }
    }

    /**
     * Shows the dialog for opening a workspace file.
     */
    public void openWorkspace() {
        SFileChooser simulationChooser = new SFileChooser(
                workspace.getCurrentDirectory(), "Zip Archive", "zip");
        File simFile = simulationChooser.showOpenDialog();
        if (simFile != null) {
            workspace.openWorkspace(simFile);
            workspace.setCurrentDirectory(simulationChooser
                    .getCurrentLocation());
            workspace.setCurrentFile(simFile);
        }
    }

    /**
     * Show Gui View of a workspace component. Used from terminal.
     *
     * @param component component to view
     */
    public static void showJFrame(final WorkspaceComponent component) {

        SimbrainDesktop.registerComponents();
        GenericJFrame theFrame = new GenericJFrame();
        GuiComponent<?> desktopComponent = createDesktopComponent(theFrame,
                component);
        theFrame.setResizable(true);
        theFrame.setVisible(true);
        theFrame.setBounds(100, 100, 200, 200);
        theFrame.setContentPane(desktopComponent);
        desktopComponent.postAddInit();
    }

    /**
     * Show a save-as dialog.
     */
    public void saveAs() {

        // Create the file chooser
        SFileChooser chooser = new SFileChooser(
                workspace.getCurrentDirectory(), "Zip Archive", "zip");

        // Set the file
        File theFile;
        if (workspace.getCurrentFile() != null) {
            theFile = chooser.showSaveDialog(workspace.getCurrentFile());
        } else {
            // Default workspace
            theFile = chooser.showSaveDialog("workspace");
        }

        // Save the file by setting the current file
        if (theFile != null) {
            workspace.setCurrentFile(theFile);
            workspace.setCurrentDirectory(chooser.getCurrentLocation());
            save(theFile);
        }
    }

    /**
     * If changes exist, show a change dialog, otherwise just save the current
     * file.
     */
    public void save() {

        // Ignore the save command if there are no changes
        if (workspace.changesExist()) {
            if (workspace.getCurrentFile() != null) {
                save(workspace.getCurrentFile());
            } else {
                saveAs(); // Show save-as if there is no current file.
            }
        }
    }

    /**
     * Save a specified file.
     *
     * @param file file to save.
     */
    private void save(File file) {
        if (file != null) {
            frame.setTitle(file.getName());
            WorkspaceSerializer.save(file, workspace);
        }
    }

    /**
     * Clear desktop of all components. Show a save-as dialog if there have been
     * changes.
     */
    public void clearDesktop() {

        // If there have been changes, show a save-as dialog
        if (workspace.changesExist()) {
            int s = showHasChangedDialog();
            if (s == JOptionPane.OK_OPTION) {
                save();
                workspace.clearWorkspace();
            } else if (s == JOptionPane.NO_OPTION) {
                workspace.clearWorkspace();
            } else if (s == JOptionPane.CANCEL_OPTION) {
                return;
            }
        } else {
            // If there have been no changes, just clear away!
            workspace.clearWorkspace();
        }
    }

    /**
     * Create the GUI and show it. For thread safety, this method should be
     * invoked from the event-dispatching thread.
     */
    private void createAndShowGUI() {
        /*
         * Make sure we have nice window decorations.
         * JFrame.setDefaultLookAndFeelDecorated(true); Create and set up the
         * window.
         */
        frame.setDefaultCloseOperation(JFrame.DO_NOTHING_ON_CLOSE);

        /** Open a default workspace */
        // openWorkspace(workspace.getCurrentFile());

        /* Display the window. */
        frame.setVisible(true);

    }

    /**
     * Simbrain main method. Creates a single instance of the Simulation class
     *
     * @param args currently not used
     */
    public static void main(final String[] args) {

        final Workspace workspace = new Workspace();

        InterceptingEventQueue eventQueue = new InterceptingEventQueue(
                workspace);

        workspace.setTaskSynchronizationManager(eventQueue);

        Toolkit.getDefaultToolkit().getSystemEventQueue().push(eventQueue);

        SwingUtilities.invokeLater(new Runnable() {
            public void run() {
                new SimbrainDesktop(workspace).createAndShowGUI();
            }
        });
    }

    /**
     * Checks to see if anything has changed and then offers to save if true.
     *
     * @return the JOptionPane pane result
     */
    private int showHasChangedDialog() {
        Object[] options = { "Save", "Don't Save", "Cancel" };
        return JOptionPane.showOptionDialog(frame,
                "The workspace has changed since last save,"
                        + "\nWould you like to save these changes?",
                "Workspace Has Changed", JOptionPane.YES_NO_OPTION,
                JOptionPane.WARNING_MESSAGE, null, options, options[0]);
    }

    /**
     * Quit application.
     *
     * @param forceQuit should quit be forced.
     */
    public void quit(final boolean forceQuit) {

        if (workspace.changesExist() && (!forceQuit)
                && (workspace.getComponentList().size() > 0)) {
            int s = showHasChangedDialog();
            if (s == JOptionPane.OK_OPTION) {
                save();
                quit(true);
            } else if (s == JOptionPane.NO_OPTION) {
                quit(true);
            } else if (s == JOptionPane.CANCEL_OPTION) {
                return;
            }
        } else {
            workspace.removeAllComponents();
            System.exit(0);
        }
    }

    /** Listener for mouse presses. */
    private final MouseListener mouseListener = new MouseAdapter() {
        /**
         * Responds to mouse events.
         *
         * @param mouseEvent Mouse Event
         */
        public void mousePressed(final MouseEvent mouseEvent) {
            Point lastClickedPoint = mouseEvent.getPoint();
            // System.out.println("desktop-->" + lastClickedPoint); //TODO: Make
            // this visible somehow
            if (mouseEvent.isControlDown()
                    || (mouseEvent.getButton() == MouseEvent.BUTTON3)) {
                contextMenu.show(frame, (int) lastClickedPoint.getX()
                        + MENU_X_OFFSET, (int) lastClickedPoint.getY()
                        + MENU_Y_OFFSET);
            }
        }
    };

    /** listener for window closing events. */
    private final WindowListener windowListener = new WindowAdapter() {
        /**
         * Responds to window closing events.
         *
         * @param arg0 Window event
         */
        public void windowClosing(final WindowEvent arg0) {
            quit(false);
        }
    };

    /** listens to menu events for setting save enabled. */
    private final MenuListener menuListener = new MenuListener() {
        /**
         * Responds to menu selected events.
         *
         * @param arg0 Menu event
         */
        public void menuSelected(final MenuEvent arg0) {
            if (workspace.changesExist()) {
                actionManager.getSaveWorkspaceAction().setEnabled(true);
            } else {
                actionManager.getSaveWorkspaceAction().setEnabled(false);
            }
        }

        /**
         * Responds to menu deslected events.
         *
         * @param arg0 Menu event
         */
        public void menuDeselected(final MenuEvent arg0) {
            /* no implementation */
        }

        /**
         * Responds to menu canceled events.
         *
         * @param arg0 Menu event
         */
        public void menuCanceled(final MenuEvent arg0) {
            /* no implementation */
        }
    };

    /**
     * Provisional Code for toggling tab dock's visibility.
     */
    public void toggleDock() {
        if (dockVisible) {
            dockVisible = false;
            horizontalSplitter.getBottomComponent().setVisible(false);
        } else {
            dockVisible = true;
            horizontalSplitter.getBottomComponent().setVisible(true);
            horizontalSplitter.setDividerLocation(getDividerLocation());
        }

    }

    /**
     * Update time label.
     */
    public void updateTimeLabel() {
        timeLabel.setText("Time:" + workspace.getTime());
        if (workspace.getUpdater().isRunning()) {
            runningLabel.setVisible(true);
            // SimbrainDesktop.this.desktop.setBackground(Color.red);
        } else {
            runningLabel.setVisible(false);
            // SimbrainDesktop.this.desktop.setBackground(Color.blue);
        }
    }

    /**
     * Helper method for determining where the bottom tab should be placed.
     *
     * @return the location
     */
    private int getDividerLocation() {
        return (int) (3 * (workspaceBounds.getHeight() / 4));
    }

    /**
     * Returns the width of the visible portion of the desktop.
     *
     * @return visible width.
     */
    public double getWidth() {
        return desktop.getVisibleRect().getWidth();
    }

    /**
     * Returns the height of the visible portion of the desktop.
     *
     * @return the visible height
     */
    public double getHeight() {
        return desktop.getVisibleRect().getHeight();
    }

    /**
     * Returns the internal desktop object. Sometimes useful in scripts.
     * @return
     */
    public JDesktopPane getDesktop() {
        return desktop;
    }

}
