/**
 * NERDSS GUI - Main Application
 * 
 * This file serves as the entry point for the NERDSS GUI application.
 * It initializes the application and manages the main window.
 * 
 * The application allows users to:
 * - Create and configure molecules with binding sites
 * - Define reactions between molecules
 * - Set simulation parameters
 * - Visualize molecular structures and reactions in 3D
 * - Export configuration files for simulations
 */

import { MoleculeManager } from './MoleculeManager.js';
import { ReactionManager } from './ReactionManager.js';
import { ParameterManager } from './ParameterManager.js';
import { BoundaryManager } from './BoundaryManager.js';
import { FileManager } from './FileManager.js';
import { MoleculeTab } from './ui/MoleculeTab.js';
import { ReactionTab } from './ui/ReactionTab.js';
import { ParameterTab } from './ui/ParameterTab.js';
import { BoundaryTab } from './ui/BoundaryTab.js';
import { Camera } from './graphics/Camera.js';
import { Renderer } from './graphics/Renderer.js';
import { Molecule } from './models/Molecule.js';
import { Reaction } from './models/Reaction.js';

/**
 * Main application class that coordinates all components of the NERDSS GUI
 */
class NERDSSApp {
    /**
     * Initialize the application
     */
    constructor() {
        // Application state managers
        this.moleculeManager = new MoleculeManager();
        this.reactionManager = new ReactionManager();
        this.parameterManager = new ParameterManager();
        this.boundaryManager = new BoundaryManager();
        this.fileManager = new FileManager(
            this.moleculeManager, 
            this.reactionManager, 
            this.parameterManager,
            this.boundaryManager
        );
        
        // Create UI tabs
        this.moleculeTab = new MoleculeTab(this.moleculeManager);
        this.reactionTab = new ReactionTab(this.reactionManager, this.moleculeManager);
        this.parameterTab = new ParameterTab(this.parameterManager);
        this.boundaryTab = new BoundaryTab(this.boundaryManager);
        
        // Initialize 3D visualization components
        this.renderer = new Renderer();
        this.moleculeCamera = new Camera();
        this.reactionCamera = new Camera();
        
        this.init();
    }
    
    /**
     * Initialize the application UI and event handlers
     */
    init() {
        this.initUI();
        this.initEventHandlers();
        this.setupTabbedPane();
        console.log('NERDSS GUI initialized');
    }
    
    /**
     * Set up the main UI components
     */
    initUI() {
        // Create main window and menu bar
        this.createMainWindow();
        this.createMenuBar();
        
        // Add tabs to the tabbed pane
        this.tabbedPane = document.getElementById('tabbed-pane');
        this.tabbedPane.appendChild(this.moleculeTab.getTabElement());
        this.tabbedPane.appendChild(this.reactionTab.getTabElement());
        this.tabbedPane.appendChild(this.parameterTab.getTabElement());
        this.tabbedPane.appendChild(this.boundaryTab.getTabElement());
        
        // Initialize the 3D visualization canvas
        this.initMoleculeCanvas();
        this.initReactionCanvas();
    }
    
    /**
     * Create the main application window
     */
    createMainWindow() {
        document.title = 'NERDSS GUI';
        
        // Create main container
        const container = document.createElement('div');
        container.id = 'nerdss-container';
        container.className = 'main-container';
        document.body.appendChild(container);
        
        // Create tabbed pane
        const tabbedPane = document.createElement('div');
        tabbedPane.id = 'tabbed-pane';
        tabbedPane.className = 'tabbed-pane';
        container.appendChild(tabbedPane);
        
        // Create status bar
        const statusBar = document.createElement('div');
        statusBar.id = 'status-bar';
        statusBar.className = 'status-bar';
        statusBar.textContent = 'Ready';
        container.appendChild(statusBar);
    }
    
    /**
     * Create the application menu bar
     */
    createMenuBar() {
        const menuBar = document.createElement('div');
        menuBar.className = 'menu-bar';
        
        // Create File menu
        const fileMenu = document.createElement('div');
        fileMenu.className = 'menu';
        
        const fileMenuButton = document.createElement('button');
        fileMenuButton.textContent = 'File';
        fileMenu.appendChild(fileMenuButton);
        
        const fileMenuDropdown = document.createElement('div');
        fileMenuDropdown.className = 'menu-dropdown';
        
        // Add export option
        const exportOption = document.createElement('div');
        exportOption.className = 'menu-item';
        exportOption.textContent = 'Export';
        exportOption.addEventListener('click', () => this.fileManager.exportFiles());
        fileMenuDropdown.appendChild(exportOption);
        
        // Add exit option
        const exitOption = document.createElement('div');
        exitOption.className = 'menu-item';
        exitOption.textContent = 'Exit';
        exitOption.addEventListener('click', () => this.exitApplication());
        fileMenuDropdown.appendChild(exitOption);
        
        fileMenu.appendChild(fileMenuDropdown);
        menuBar.appendChild(fileMenu);
        
        // Add menu bar to the document
        document.body.insertBefore(menuBar, document.body.firstChild);
    }
    
    /**
     * Initialize the molecule visualization canvas
     */
    initMoleculeCanvas() {
        const canvas = document.createElement('canvas');
        canvas.id = 'molecule-canvas';
        canvas.width = 400;
        canvas.height = 300;
        this.moleculeTab.addCanvas(canvas);
        
        this.renderer.initMoleculeRenderer(canvas, this.moleculeCamera);
        this.moleculeCamera.installTrackball(canvas);
    }
    
    /**
     * Initialize the reaction visualization canvas
     */
    initReactionCanvas() {
        const canvas = document.createElement('canvas');
        canvas.id = 'reaction-canvas';
        canvas.width = 400;
        canvas.height = 300;
        this.reactionTab.addCanvas(canvas);
        
        this.renderer.initReactionRenderer(canvas, this.reactionCamera);
        this.reactionCamera.installTrackball(canvas);
    }
    
    /**
     * Set up the tabbed pane for switching between views
     */
    setupTabbedPane() {
        // Implement tabbed pane functionality
        const tabs = document.querySelectorAll('.tab-button');
        const tabContents = document.querySelectorAll('.tab-content');
        
        tabs.forEach(tab => {
            tab.addEventListener('click', () => {
                // Remove active class from all tabs
                tabs.forEach(t => t.classList.remove('active'));
                tabContents.forEach(content => content.classList.remove('active'));
                
                // Add active class to current tab
                tab.classList.add('active');
                const tabId = tab.getAttribute('data-tab');
                document.getElementById(tabId).classList.add('active');
                
                // Update 3D visualization if needed
                if (tabId === 'molecule-tab') {
                    this.renderer.refreshMoleculeView();
                } else if (tabId === 'reaction-tab') {
                    this.renderer.refreshReactionView();
                }
            });
        });
        
        // Activate the first tab by default
        if (tabs.length > 0) {
            tabs[0].click();
        }
    }
    
    /**
     * Initialize global event handlers
     */
    initEventHandlers() {
        // Handle window resize
        window.addEventListener('resize', () => {
            this.renderer.resizeCanvases();
        });
        
        // Setup communication between tabs
        this.moleculeTab.on('molecule-added', (molecule) => {
            this.reactionTab.updateMoleculeList();
        });
        
        this.moleculeTab.on('molecule-updated', (molecule) => {
            this.reactionTab.updateMoleculeList();
        });
        
        this.moleculeTab.on('molecule-deleted', (moleculeId) => {
            this.reactionTab.updateMoleculeList();
            this.reactionManager.removeMoleculeFromReactions(moleculeId);
        });
    }
    
    /**
     * Exit the application after confirmation
     */
    exitApplication() {
        const confirmed = confirm('Are you sure you want to exit?');
        if (confirmed) {
            window.close();
            // For web applications, might redirect to a different page
            // window.location.href = 'exit.html';
        }
    }
    
    /**
     * Application entry point
     */
    static main() {
        // Create and initialize the application
        const app = new NERDSSApp();
        
        // Register global error handler
        window.addEventListener('error', (event) => {
            console.error('Application error:', event.error);
            alert(`An error occurred: ${event.error.message}`);
        });
    }
}

// Start the application when DOM is fully loaded
document.addEventListener('DOMContentLoaded', () => {
    NERDSSApp.main();
});

// Export the NERDSSApp class for testing or external access
export { NERDSSApp };