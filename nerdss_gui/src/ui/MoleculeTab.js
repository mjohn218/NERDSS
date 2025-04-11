/**
 * MoleculeTab.js
 * 
 * This file defines the MoleculeTab class, which provides the UI for
 * creating, editing, and managing molecules in the NERDSS GUI.
 */

import { Molecule, BindingSite } from '../models/Molecule.js';
import { Vector3 } from '../math/Vector3.js';

/**
 * MoleculeTab class for the molecule editing UI
 */
class MoleculeTab {
    /**
     * Create a new MoleculeTab
     * @param {MoleculeManager} moleculeManager - Reference to the molecule manager
     */
    constructor(moleculeManager) {
        this.moleculeManager = moleculeManager;
        this.tabElement = null;
        this.moleculeCanvas = null;
        this.bindingSiteTable = null;
        this.moleculeTable = null;
        this.currentMolecule = null;
        this.currentBindingSiteIndex = -1;
        
        // Create the tab
        this.createTab();
        
        // Set up event listeners
        this.setupEventListeners();
    }
    
    /**
     * Get the tab element
     * @returns {HTMLElement} The tab element
     */
    getTabElement() {
        return this.tabElement;
    }
    
    /**
     * Add a canvas to the tab for rendering
     * @param {HTMLCanvasElement} canvas - Canvas element
     */
    addCanvas(canvas) {
        this.moleculeCanvas = canvas;
        const canvasContainer = this.tabElement.querySelector('#molecule-canvas-container');
        if (canvasContainer) {
            canvasContainer.innerHTML = '';
            canvasContainer.appendChild(canvas);
        }
    }
    
    /**
     * Create the tab element and its contents
     * @private
     */
    createTab() {
        // Create the tab element
        this.tabElement = document.createElement('div');
        this.tabElement.id = 'molecule-tab';
        this.tabElement.className = 'tab-content';
        
        // Create tab header with controls
        const header = document.createElement('div');
        header.className = 'tab-header';
        header.innerHTML = `
            <div class="tab-title">Molecules</div>
            <div class="tab-controls">
                <button id="add-molecule-btn" class="btn">Add Molecule</button>
                <button id="update-molecule-btn" class="btn">Update Molecule</button>
                <button id="delete-molecule-btn" class="btn">Delete Molecule</button>
            </div>
        `;
        this.tabElement.appendChild(header);
        
        // Create the main content area with a two-column layout
        const content = document.createElement('div');
        content.className = 'tab-content-area';
        
        // Left column: molecule editor
        const leftColumn = document.createElement('div');
        leftColumn.className = 'left-column';
        leftColumn.innerHTML = `
            <div class="panel molecule-editor-panel">
                <div class="panel-header">Molecule Properties</div>
                <div class="panel-content">
                    <div class="form-group">
                        <label for="molecule-name">Molecule Name</label>
                        <input id="molecule-name" type="text" value="MolecX">
                    </div>
                    <div class="form-group">
                        <label for="molecule-count">Molecule Count</label>
                        <input id="molecule-count" type="number" value="100" min="1">
                    </div>
                    <div class="form-group">
                        <label for="molecule-radius">Molecule Radius (nm)</label>
                        <input id="molecule-radius" type="number" value="10.0" min="0">
                        <button id="calc-diff-btn" class="btn small">Calculate Diffusion</button>
                    </div>
                    <div class="form-group">
                        <label for="trans-diff">Translational Diffusion (μm²/s)</label>
                        <input id="trans-diff" type="number" value="0.0" min="0">
                    </div>
                    <div class="form-group">
                        <label for="rot-diff">Rotational Diffusion (1/μs)</label>
                        <input id="rot-diff" type="number" value="0.0" min="0">
                    </div>
                    <div class="form-group">
                        <div class="checkbox-group">
                            <input id="is-membrane" type="checkbox">
                            <label for="is-membrane">Anchored to Membrane</label>
                        </div>
                    </div>
                    <div class="form-group">
                        <div class="checkbox-group">
                            <input id="is-lipid" type="checkbox">
                            <label for="is-lipid">Is Lipid</label>
                        </div>
                    </div>
                    <div class="form-group">
                        <div class="checkbox-group">
                            <input id="is-implicit-lipid" type="checkbox">
                            <label for="is-implicit-lipid">Is Implicit Lipid</label>
                        </div>
                    </div>
                </div>
            </div>
            
            <div class="panel binding-site-panel">
                <div class="panel-header">Binding Site Properties</div>
                <div class="panel-content">
                    <div class="form-group">
                        <label for="site-name">Site Name</label>
                        <input id="site-name" type="text" value="Site1">
                    </div>
                    <div class="form-group">
                        <label>Coordinates Relative to C.o.M. (nm)</label>
                        <div class="coord-inputs">
                            <div class="coord-input">
                                <label for="site-x">x:</label>
                                <input id="site-x" type="number" value="3.0" step="0.1">
                            </div>
                            <div class="coord-input">
                                <label for="site-y">y:</label>
                                <input id="site-y" type="number" value="4.0" step="0.1">
                            </div>
                            <div class="coord-input">
                                <label for="site-z">z:</label>
                                <input id="site-z" type="number" value="5.0" step="0.1">
                            </div>
                        </div>
                    </div>
                    <div class="form-group binding-site-btns">
                        <button id="add-site-btn" class="btn small">Add</button>
                        <button id="update-site-btn" class="btn small">Update</button>
                        <button id="delete-site-btn" class="btn small">Delete</button>
                    </div>
                </div>
            </div>
            
            <div id="molecule-error" class="error-message"></div>
        `;
        
        // Right column: tables and visualization
        const rightColumn = document.createElement('div');
        rightColumn.className = 'right-column';
        rightColumn.innerHTML = `
            <div class="panel binding-sites-panel">
                <div class="panel-header">Binding Sites</div>
                <div class="panel-content">
                    <div class="table-container">
                        <table id="binding-site-table" class="data-table">
                            <thead>
                                <tr>
                                    <th>Name</th>
                                    <th>x (nm)</th>
                                    <th>y (nm)</th>
                                    <th>z (nm)</th>
                                    <th>#States</th>
                                    <th>States</th>
                                </tr>
                            </thead>
                            <tbody>
                                <tr>
                                    <td>Center of Mass</td>
                                    <td>0.0</td>
                                    <td>0.0</td>
                                    <td>0.0</td>
                                    <td>0</td>
                                    <td>-</td>
                                </tr>
                            </tbody>
                        </table>
                    </div>
                </div>
            </div>
            
            <div class="panel molecule-viz-panel">
                <div class="panel-header">3D Visualization</div>
                <div class="panel-content">
                    <div id="molecule-canvas-container" class="canvas-container"></div>
                </div>
            </div>
            
            <div class="panel molecules-list-panel">
                <div class="panel-header">Molecules</div>
                <div class="panel-content">
                    <div class="table-container">
                        <table id="molecule-table" class="data-table">
                            <thead>
                                <tr>
                                    <th>Name</th>
                                    <th>Count</th>
                                    <th>On Membrane</th>
                                    <th>Trans. Diff. (μm²/s)</th>
                                    <th>Rot. Diff. (1/μs)</th>
                                    <th>#Sites</th>
                                </tr>
                            </thead>
                            <tbody>
                                <!-- Molecules will be populated here -->
                            </tbody>
                        </table>
                    </div>
                </div>
            </div>
        `;
        
        // Append columns to content
        content.appendChild(leftColumn);
        content.appendChild(rightColumn);
        
        // Append content to tab
        this.tabElement.appendChild(content);
        
        // Store references to tables
        this.bindingSiteTable = this.tabElement.querySelector('#binding-site-table');
        this.moleculeTable = this.tabElement.querySelector('#molecule-table');
    }
    
    /**
     * Set up event listeners for the tab
     * @private
     */
    setupEventListeners() {
        // Add molecule button
        const addMoleculeBtn = this.tabElement.querySelector('#add-molecule-btn');
        if (addMoleculeBtn) {
            addMoleculeBtn.addEventListener('click', () => this.addMolecule());
        }
        
        // Update molecule button
        const updateMoleculeBtn = this.tabElement.querySelector('#update-molecule-btn');
        if (updateMoleculeBtn) {
            updateMoleculeBtn.addEventListener('click', () => this.updateMolecule());
        }
        
        // Delete molecule button
        const deleteMoleculeBtn = this.tabElement.querySelector('#delete-molecule-btn');
        if (deleteMoleculeBtn) {
            deleteMoleculeBtn.addEventListener('click', () => this.deleteMolecule());
        }
        
        // Add binding site button
        const addSiteBtn = this.tabElement.querySelector('#add-site-btn');
        if (addSiteBtn) {
            addSiteBtn.addEventListener('click', () => this.addBindingSite());
        }
        
        // Update binding site button
        const updateSiteBtn = this.tabElement.querySelector('#update-site-btn');
        if (updateSiteBtn) {
            updateSiteBtn.addEventListener('click', () => this.updateBindingSite());
        }
        
        // Delete binding site button
        const deleteSiteBtn = this.tabElement.querySelector('#delete-site-btn');
        if (deleteSiteBtn) {
            deleteSiteBtn.addEventListener('click', () => this.deleteBindingSite());
        }
        
        // Calculate diffusion button
        const calcDiffBtn = this.tabElement.querySelector('#calc-diff-btn');
        if (calcDiffBtn) {
            calcDiffBtn.addEventListener('click', () => this.calculateDiffusion());
        }
        
        // Molecule table row click
        if (this.moleculeTable) {
            this.moleculeTable.addEventListener('click', event => {
                const row = event.target.closest('tr');
                if (row && row.dataset.moleculeId) {
                    this.selectMolecule(row.dataset.moleculeId);
                }
            });
        }
        
        // Binding site table row click
        if (this.bindingSiteTable) {
            this.bindingSiteTable.addEventListener('click', event => {
                const row = event.target.closest('tr');
                if (row && row.dataset.siteIndex) {
                    this.selectBindingSite(parseInt(row.dataset.siteIndex));
                }
            });
        }
        
        // Molecule manager events
        if (this.moleculeManager) {
            this.moleculeManager.on('molecules-changed', () => this.refreshMoleculeTable());
            this.moleculeManager.on('selection-changed', selection => this.handleSelectionChanged(selection));
            this.moleculeManager.on('error', errorMsg => this.showError(errorMsg));
        }
    }
    
    /**
     * Handle selection change event from the molecule manager
     * @param {Object} selection - Selected molecule and binding site
     * @private
     */
    handleSelectionChanged(selection) {
        // Update UI to reflect the new selection
        this.currentMolecule = selection.molecule;
        this.currentBindingSiteIndex = selection.bindingSiteIndex;
        
        // Update form fields
        this.populateMoleculeForm(selection.molecule);
        
        // Highlight selected rows in tables
        this.highlightSelectedMolecule(selection.moleculeId);
        this.highlightSelectedBindingSite(selection.bindingSiteIndex);
        
        // Update binding site table
        this.refreshBindingSiteTable();
        
        // Update binding site form if a site is selected
        if (selection.bindingSite) {
            this.populateBindingSiteForm(selection.bindingSite);
        } else {
            this.clearBindingSiteForm();
        }
    }
    
    /**
     * Populate the molecule form with data from a molecule
     * @param {Molecule} molecule - Molecule to load
     * @private
     */
    populateMoleculeForm(molecule) {
        if (!molecule) {
            this.clearMoleculeForm();
            return;
        }
        
        const nameInput = this.tabElement.querySelector('#molecule-name');
        const countInput = this.tabElement.querySelector('#molecule-count');
        const transDiffInput = this.tabElement.querySelector('#trans-diff');
        const rotDiffInput = this.tabElement.querySelector('#rot-diff');
        const membraneCheckbox = this.tabElement.querySelector('#is-membrane');
        const lipidCheckbox = this.tabElement.querySelector('#is-lipid');
        const implicitLipidCheckbox = this.tabElement.querySelector('#is-implicit-lipid');
        
        if (nameInput) nameInput.value = molecule.name;
        if (countInput) countInput.value = molecule.count;
        if (transDiffInput) transDiffInput.value = molecule.diffusionTranslational;
        if (rotDiffInput) rotDiffInput.value = molecule.diffusionRotational;
        if (membraneCheckbox) membraneCheckbox.checked = molecule.membraneAnchored;
        if (lipidCheckbox) lipidCheckbox.checked = molecule.isLipid;
        if (implicitLipidCheckbox) implicitLipidCheckbox.checked = molecule.isImplicitLipid;
    }
    
    /**
     * Clear the molecule form
     * @private
     */
    clearMoleculeForm() {
        const nameInput = this.tabElement.querySelector('#molecule-name');
        const countInput = this.tabElement.querySelector('#molecule-count');
        const radiusInput = this.tabElement.querySelector('#molecule-radius');
        const transDiffInput = this.tabElement.querySelector('#trans-diff');
        const rotDiffInput = this.tabElement.querySelector('#rot-diff');
        const membraneCheckbox = this.tabElement.querySelector('#is-membrane');
        const lipidCheckbox = this.tabElement.querySelector('#is-lipid');
        const implicitLipidCheckbox = this.tabElement.querySelector('#is-implicit-lipid');
        
        if (nameInput) nameInput.value = 'MolecX';
        if (countInput) countInput.value = '100';
        if (radiusInput) radiusInput.value = '10.0';
        if (transDiffInput) transDiffInput.value = '0.0';
        if (rotDiffInput) rotDiffInput.value = '0.0';
        if (membraneCheckbox) membraneCheckbox.checked = false;
        if (lipidCheckbox) lipidCheckbox.checked = false;
        if (implicitLipidCheckbox) implicitLipidCheckbox.checked = false;
    }
    
    /**
     * Populate the binding site form with data from a binding site
     * @param {BindingSite} bindingSite - Binding site to load
     * @private
     */
    populateBindingSiteForm(bindingSite) {
        if (!bindingSite) {
            this.clearBindingSiteForm();
            return;
        }
        
        const nameInput = this.tabElement.querySelector('#site-name');
        const xInput = this.tabElement.querySelector('#site-x');
        const yInput = this.tabElement.querySelector('#site-y');
        const zInput = this.tabElement.querySelector('#site-z');
        
        if (nameInput) nameInput.value = bindingSite.name;
        if (xInput) xInput.value = bindingSite.position.x;
        if (yInput) yInput.value = bindingSite.position.y;
        if (zInput) zInput.value = bindingSite.position.z;
    }
    
    /**
     * Clear the binding site form
     * @private
     */
    clearBindingSiteForm() {
        const nameInput = this.tabElement.querySelector('#site-name');
        const xInput = this.tabElement.querySelector('#site-x');
        const yInput = this.tabElement.querySelector('#site-y');
        const zInput = this.tabElement.querySelector('#site-z');
        
        if (nameInput) nameInput.value = 'Site1';
        if (xInput) xInput.value = '3.0';
        if (yInput) yInput.value = '4.0';
        if (zInput) zInput.value = '5.0';
    }
    
    /**
     * Add a new molecule
     * @private
     */
    addMolecule() {
        // Clear any previous errors
        this.clearError();
        
        // Get form values
        const name = this.tabElement.querySelector('#molecule-name').value;
        const count = parseInt(this.tabElement.querySelector('#molecule-count').value);
        const transDiff = parseFloat(this.tabElement.querySelector('#trans-diff').value);
        const rotDiff = parseFloat(this.tabElement.querySelector('#rot-diff').value);
        const isOnMembrane = this.tabElement.querySelector('#is-membrane').checked;
        const isLipid = this.tabElement.querySelector('#is-lipid').checked;
        const isImplicitLipid = this.tabElement.querySelector('#is-implicit-lipid').checked;
        
        // Validate inputs
        if (!name || name.trim() === '') {
            this.showError('Molecule name cannot be empty');
            return;
        }
        
        if (isNaN(count) || count <= 0) {
            this.showError('Molecule count must be a positive integer');
            return;
        }
        
        if (isNaN(transDiff) || transDiff < 0) {
            this.showError('Translational diffusion must be a non-negative number');
            return;
        }
        
        if (isNaN(rotDiff) || rotDiff < 0) {
            this.showError('Rotational diffusion must be a non-negative number');
            return;
        }
        
        // Create a new molecule
        const molecule = new Molecule(
            name,
            count,
            transDiff,
            rotDiff,
            isOnMembrane,
            isLipid,
            isImplicitLipid
        );
        
        // Add binding sites from the binding site table (skip Center of Mass)
        const siteRows = Array.from(this.bindingSiteTable.querySelectorAll('tbody tr')).slice(1);
        for (const row of siteRows) {
            const cells = row.cells;
            const siteName = cells[0].textContent;
            const x = parseFloat(cells[1].textContent);
            const y = parseFloat(cells[2].textContent);
            const z = parseFloat(cells[3].textContent);
            const stateCount = parseInt(cells[4].textContent);
            
            // Parse states
            let stateNames = [];
            if (cells[5].textContent !== '-') {
                stateNames = cells[5].textContent.split(',').map(s => s.trim());
            }
            
            molecule.addBindingSite(new BindingSite(
                siteName,
                new Vector3(x, y, z),
                stateCount,
                stateNames
            ));
        }
        
        // Add the molecule
        if (this.moleculeManager.addMolecule(molecule)) {
            // Clear forms
            this.clearMoleculeForm();
            this.clearBindingSiteForm();
            
            // Reset binding site table
            this.resetBindingSiteTable();
            
            // Select the new molecule
            const newMoleculeId = molecule.id;
            this.moleculeManager.setSelection(newMoleculeId);
        }
    }
    
    /**
     * Update the selected molecule
     * @private
     */
    updateMolecule() {
        // Clear any previous errors
        this.clearError();
        
        // Make sure a molecule is selected
        if (!this.currentMolecule) {
            this.showError('No molecule selected');
            return;
        }
        
        // Get form values
        const name = this.tabElement.querySelector('#molecule-name').value;
        const count = parseInt(this.tabElement.querySelector('#molecule-count').value);
        const transDiff = parseFloat(this.tabElement.querySelector('#trans-diff').value);
        const rotDiff = parseFloat(this.tabElement.querySelector('#rot-diff').value);
        const isOnMembrane = this.tabElement.querySelector('#is-membrane').checked;
        const isLipid = this.tabElement.querySelector('#is-lipid').checked;
        const isImplicitLipid = this.tabElement.querySelector('#is-implicit-lipid').checked;
        
        // Validate inputs
        if (!name || name.trim() === '') {
            this.showError('Molecule name cannot be empty');
            return;
        }
        
        if (isNaN(count) || count <= 0) {
            this.showError('Molecule count must be a positive integer');
            return;
        }
        
        if (isNaN(transDiff) || transDiff < 0) {
            this.showError('Translational diffusion must be a non-negative number');
            return;
        }
        
        if (isNaN(rotDiff) || rotDiff < 0) {
            this.showError('Rotational diffusion must be a non-negative number');
            return;
        }
        
        // Collect binding sites
        const bindingSites = [];
        const siteRows = Array.from(this.bindingSiteTable.querySelectorAll('tbody tr')).slice(1); // Skip COM
        for (const row of siteRows) {
            const cells = row.cells;
            const siteName = cells[0].textContent;
            const x = parseFloat(cells[1].textContent);
            const y = parseFloat(cells[2].textContent);
            const z = parseFloat(cells[3].textContent);
            const stateCount = parseInt(cells[4].textContent);
            
            // Parse states
            let stateNames = [];
            if (cells[5].textContent !== '-') {
                stateNames = cells[5].textContent.split(',').map(s => s.trim());
            }
            
            bindingSites.push({
                name: siteName,
                x: x,
                y: y,
                z: z,
                stateCount: stateCount,
                stateNames: stateNames
            });
        }
        
        // Update the molecule
        if (this.moleculeManager.updateMolecule(this.currentMolecule.id, {
            name: name,
            count: count,
            diffusionTranslational: transDiff,
            diffusionRotational: rotDiff,
            membraneAnchored: isOnMembrane,
            isLipid: isLipid,
            isImplicitLipid: isImplicitLipid,
            bindingSites: bindingSites
        })) {
            // Keep selection on the updated molecule
            this.moleculeManager.setSelection(this.currentMolecule.id);
        }
    }
    
    /**
     * Delete the selected molecule
     * @private
     */
    deleteMolecule() {
        if (!this.currentMolecule) {
            this.showError('No molecule selected');
            return;
        }
        
        // Ask for confirmation
        if (confirm(`Are you sure you want to delete the molecule '${this.currentMolecule.name}'?`)) {
            const success = this.moleculeManager.removeMolecule(this.currentMolecule.id);
            if (success) {
                this.clearMoleculeForm();
                this.clearBindingSiteForm();
                this.resetBindingSiteTable();
                this.currentMolecule = null;
                this.currentBindingSiteIndex = -1;
            }
        }
    }
    
    /**
     * Add a binding site to the current molecule
     * @private
     */
    addBindingSite() {
        // Make sure a molecule is selected
        if (!this.currentMolecule) {
            this.showError('No molecule selected');
            return;
        }
        
        // Get form values
        const name = this.tabElement.querySelector('#site-name').value;
        const x = parseFloat(this.tabElement.querySelector('#site-x').value);
        const y = parseFloat(this.tabElement.querySelector('#site-y').value);
        const z = parseFloat(this.tabElement.querySelector('#site-z').value);
        
        // Validate inputs
        if (!name || name.trim() === '') {
            this.showError('Binding site name cannot be empty');
            return;
        }
        
        if (isNaN(x) || isNaN(y) || isNaN(z)) {
            this.showError('Coordinates must be numbers');
            return;
        }
        
        // Check for interface state information
        const stateInfo = this.promptForStateInfo();
        if (stateInfo === null) {
            // User cancelled
            return;
        }
        
        const { stateCount, stateNames } = stateInfo;
        
        // Add row to the binding site table
        const tbody = this.bindingSiteTable.querySelector('tbody');
        const row = tbody.insertRow();
        row.dataset.siteIndex = this.currentMolecule.bindingSites.length; // This will be the new index
        
        // Add cells
        row.insertCell().textContent = name;
        row.insertCell().textContent = x.toFixed(1);
        row.insertCell().textContent = y.toFixed(1);
        row.insertCell().textContent = z.toFixed(1);
        row.insertCell().textContent = stateCount;
        row.insertCell().textContent = stateCount > 0 ? stateNames.join(', ') : '-';
        
        // Clear binding site form
        this.clearBindingSiteForm();
    }
    
    /**
     * Update the selected binding site
     * @private
     */
    updateBindingSite() {
        // Make sure a molecule and binding site are selected
        if (!this.currentMolecule) {
            this.showError('No molecule selected');
            return;
        }
        
        if (this.currentBindingSiteIndex < 0) {
            this.showError('No binding site selected');
            return;
        }
        
        // Get form values
        const name = this.tabElement.querySelector('#site-name').value;
        const x = parseFloat(this.tabElement.querySelector('#site-x').value);
        const y = parseFloat(this.tabElement.querySelector('#site-y').value);
        const z = parseFloat(this.tabElement.querySelector('#site-z').value);
        
        // Validate inputs
        if (!name || name.trim() === '') {
            this.showError('Binding site name cannot be empty');
            return;
        }
        
        if (isNaN(x) || isNaN(y) || isNaN(z)) {
            this.showError('Coordinates must be numbers');
            return;
        }
        
        // Get current site for state info
        const currentSite = this.currentMolecule.bindingSites[this.currentBindingSiteIndex];
        
        // Check for interface state information
        const defaultStateCount = currentSite ? currentSite.stateCount : 0;
        const defaultStateNames = currentSite ? currentSite.stateNames : [];
        
        const stateInfo = this.promptForStateInfo(defaultStateCount, defaultStateNames);
        if (stateInfo === null) {
            // User cancelled
            return;
        }
        
        const { stateCount, stateNames } = stateInfo;
        
        // Update row in the binding site table
        const rows = this.bindingSiteTable.querySelectorAll('tbody tr');
        if (rows.length > this.currentBindingSiteIndex + 1) { // +1 for COM row
            const row = rows[this.currentBindingSiteIndex + 1];
            
            // Update cells
            row.cells[0].textContent = name;
            row.cells[1].textContent = x.toFixed(1);
            row.cells[2].textContent = y.toFixed(1);
            row.cells[3].textContent = z.toFixed(1);
            row.cells[4].textContent = stateCount;
            row.cells[5].textContent = stateCount > 0 ? stateNames.join(', ') : '-';
        }
        
        // Clear binding site form
        this.clearBindingSiteForm();
    }
    
    /**
     * Delete the selected binding site
     * @private
     */
    deleteBindingSite() {
        // Make sure a molecule and binding site are selected
        if (!this.currentMolecule) {
            this.showError('No molecule selected');
            return;
        }
        
        if (this.currentBindingSiteIndex < 0) {
            this.showError('No binding site selected');
            return;
        }
        
        // Ask for confirmation
        const site = this.currentMolecule.bindingSites[this.currentBindingSiteIndex];
        if (confirm(`Are you sure you want to delete the binding site '${site.name}'?`)) {
            // Remove row from the binding site table
            const rows = this.bindingSiteTable.querySelectorAll('tbody tr');
            if (rows.length > this.currentBindingSiteIndex + 1) { // +1 for COM row
                const row = rows[this.currentBindingSiteIndex + 1];
                row.remove();
            }
            
            // Update data-site-index attributes for remaining rows
            const remainingRows = this.bindingSiteTable.querySelectorAll('tbody tr');
            for (let i = 1; i < remainingRows.length; i++) { // Start at 1 to skip COM
                remainingRows[i].dataset.siteIndex = i - 1;
            }
            
            // Clear binding site form
            this.clearBindingSiteForm();
            
            // Reset selection
            this.currentBindingSiteIndex = -1;
        }
    }
    
    /**
     * Calculate diffusion coefficients based on molecule radius
     * @private
     */
    calculateDiffusion() {
        const radiusInput = this.tabElement.querySelector('#molecule-radius');
        if (!radiusInput) return;
        
        const radius = parseFloat(radiusInput.value);
        if (isNaN(radius) || radius <= 0) {
            this.showError('Radius must be a positive number');
            return;
        }
        
        // Calculate diffusion coefficients
        const diffusion = this.moleculeManager.calculateDiffusionCoefficients(radius);
        
        // Update form fields
        const transDiffInput = this.tabElement.querySelector('#trans-diff');
        const rotDiffInput = this.tabElement.querySelector('#rot-diff');
        
        if (transDiffInput) transDiffInput.value = diffusion.translational.toFixed(4);
        if (rotDiffInput) rotDiffInput.value = diffusion.rotational.toFixed(4);
    }
    
    /**
     * Select a molecule
     * @param {string} moleculeId - ID of the molecule to select
     * @private
     */
    selectMolecule(moleculeId) {
        this.moleculeManager.setSelection(moleculeId);
    }
    
    /**
     * Select a binding site
     * @param {number} siteIndex - Index of the binding site to select
     * @private
     */
    selectBindingSite(siteIndex) {
        if (!this.currentMolecule) return;
        
        this.moleculeManager.setSelection(this.currentMolecule.id, siteIndex);
    }
    
    /**
     * Refresh the molecule table
     * @private
     */
    refreshMoleculeTable() {
        const tbody = this.moleculeTable.querySelector('tbody');
        tbody.innerHTML = '';
        
        const molecules = this.moleculeManager.getAllMolecules();
        
        for (const molecule of molecules) {
            const row = tbody.insertRow();
            row.dataset.moleculeId = molecule.id;
            
            // Add cells
            row.insertCell().textContent = molecule.name;
            row.insertCell().textContent = molecule.count;
            row.insertCell().textContent = molecule.membraneAnchored ? 'Yes' : 'No';
            row.insertCell().textContent = molecule.diffusionTranslational.toFixed(4);
            row.insertCell().textContent = molecule.diffusionRotational.toFixed(4);
            row.insertCell().textContent = molecule.bindingSites.length;
            
            // Highlight if this is the selected molecule
            if (this.currentMolecule && molecule.id === this.currentMolecule.id) {
                row.classList.add('selected');
            }
        }
    }
    
    /**
     * Refresh the binding site table
     * @private
     */
    refreshBindingSiteTable() {
        // Don't refresh if no molecule is selected
        if (!this.currentMolecule) {
            this.resetBindingSiteTable();
            return;
        }
        
        const tbody = this.bindingSiteTable.querySelector('tbody');
        tbody.innerHTML = '';
        
        // Add center of mass row
        const comRow = tbody.insertRow();
        comRow.dataset.siteIndex = -1;
        comRow.insertCell().textContent = 'Center of Mass';
        comRow.insertCell().textContent = '0.0';
        comRow.insertCell().textContent = '0.0';
        comRow.insertCell().textContent = '0.0';
        comRow.insertCell().textContent = '0';
        comRow.insertCell().textContent = '-';
        
        // Add binding sites
        for (let i = 0; i < this.currentMolecule.bindingSites.length; i++) {
            const site = this.currentMolecule.bindingSites[i];
            
            const row = tbody.insertRow();
            row.dataset.siteIndex = i;
            
            // Add cells
            row.insertCell().textContent = site.name;
            row.insertCell().textContent = site.position.x.toFixed(1);
            row.insertCell().textContent = site.position.y.toFixed(1);
            row.insertCell().textContent = site.position.z.toFixed(1);
            row.insertCell().textContent = site.stateCount;
            row.insertCell().textContent = site.stateCount > 0 ? site.stateNames.join(', ') : '-';
            
            // Highlight if this is the selected binding site
            if (i === this.currentBindingSiteIndex) {
                row.classList.add('selected');
            }
        }
    }
    
    /**
     * Reset the binding site table to show only center of mass
     * @private
     */
    resetBindingSiteTable() {
        const tbody = this.bindingSiteTable.querySelector('tbody');
        tbody.innerHTML = '';
        
        // Add center of mass row
        const comRow = tbody.insertRow();
        comRow.dataset.siteIndex = -1;
        comRow.insertCell().textContent = 'Center of Mass';
        comRow.insertCell().textContent = '0.0';
        comRow.insertCell().textContent = '0.0';
        comRow.insertCell().textContent = '0.0';
        comRow.insertCell().textContent = '0';
        comRow.insertCell().textContent = '-';
    }
    
    /**
     * Highlight the selected molecule in the table
     * @param {string} moleculeId - ID of the selected molecule
     * @private
     */
    highlightSelectedMolecule(moleculeId) {
        const rows = this.moleculeTable.querySelectorAll('tbody tr');
        
        for (const row of rows) {
            if (row.dataset.moleculeId === moleculeId) {
                row.classList.add('selected');
            } else {
                row.classList.remove('selected');
            }
        }
    }
    
    /**
     * Highlight the selected binding site in the table
     * @param {number} siteIndex - Index of the selected binding site
     * @private
     */
    highlightSelectedBindingSite(siteIndex) {
        const rows = this.bindingSiteTable.querySelectorAll('tbody tr');
        
        for (const row of rows) {
            if (parseInt(row.dataset.siteIndex) === siteIndex) {
                row.classList.add('selected');
            } else {
                row.classList.remove('selected');
            }
        }
    }
    
    /**
     * Prompt for binding site state information
     * @param {number} [defaultStateCount=0] - Default number of states
     * @param {string[]} [defaultStateNames=[]] - Default state names
     * @returns {Object|null} Object with stateCount and stateNames, or null if cancelled
     * @private
     */
    promptForStateInfo(defaultStateCount = 0, defaultStateNames = []) {
        // Create dialog for state information
        const dialogHtml = `
            <div class="dialog-content">
                <div class="form-group">
                    <label for="state-count">Number of interface states</label>
                    <input id="state-count" type="number" min="0" max="10" value="${defaultStateCount}">
                </div>
                <div class="form-group">
                    <label for="state-names">List of interface states (comma-separated)</label>
                    <input id="state-names" type="text" value="${defaultStateNames.join(', ')}">
                </div>
            </div>
        `;
        
        // Show dialog
        const result = prompt("Can this site exist in additional states?\nEnter number of states and names (comma-separated):", defaultStateNames.join(', '));
        if (result === null) {
            return null;
        }
        
        // Parse result
        let stateNames = [];
        let stateCount = 0;
        
        if (result.trim() !== '') {
            stateNames = result.split(',').map(s => s.trim());
            stateCount = stateNames.length;
        }
        
        return { stateCount, stateNames };
    }
    
    /**
     * Show an error message
     * @param {string} message - Error message to display
     * @private
     */
    showError(message) {
        const errorElement = this.tabElement.querySelector('#molecule-error');
        if (errorElement) {
            errorElement.textContent = message;
            errorElement.style.display = 'block';
        }
    }
    
    /**
     * Clear error messages
     * @private
     */
    clearError() {
        const errorElement = this.tabElement.querySelector('#molecule-error');
        if (errorElement) {
            errorElement.textContent = '';
            errorElement.style.display = 'none';
        }
    }
}

export { MoleculeTab };