/**
 * MoleculeManager.js
 * 
 * This file defines the MoleculeManager class, which manages the collection
 * of molecules in the NERDSS application. It provides methods for adding,
 * updating, removing, and retrieving molecules.
 */

import { Molecule, BindingSite } from './models/Molecule.js';
import { EventEmitter } from './utils/EventEmitter.js';
import { Vector3 } from './math/Vector3.js';

/**
 * Manages the collection of molecules in the application
 * @extends EventEmitter
 */
class MoleculeManager extends EventEmitter {
    /**
     * Create a new MoleculeManager
     */
    constructor() {
        super();
        this.molecules = new Map();
        this.selectedMoleculeId = null;
        this.selectedBindingSiteIndex = -1;
    }
    
    /**
     * Get all molecules as an array
     * @returns {Molecule[]} Array of all molecules
     */
    getAllMolecules() {
        return Array.from(this.molecules.values());
    }
    
    /**
     * Get a molecule by ID
     * @param {string} id - Molecule ID
     * @returns {Molecule|undefined} The molecule or undefined if not found
     */
    getMolecule(id) {
        return this.molecules.get(id);
    }
    
    /**
     * Get molecule by name
     * @param {string} name - Molecule name
     * @returns {Molecule|undefined} The first molecule with the given name or undefined
     */
    getMoleculeByName(name) {
        for (const molecule of this.molecules.values()) {
            if (molecule.name === name) {
                return molecule;
            }
        }
        return undefined;
    }
    
    /**
     * Add a new molecule to the collection
     * @param {Molecule} molecule - Molecule to add
     * @returns {boolean} True if added successfully, false if duplicate name or validation error
     */
    addMolecule(molecule) {
        // Check if molecule with same name already exists
        if (this.getMoleculeByName(molecule.name)) {
            this.emit('error', `Molecule with name '${molecule.name}' already exists`);
            return false;
        }
        
        // Validate molecule
        if (!this.validateMolecule(molecule)) {
            return false;
        }
        
        // Add molecule
        this.molecules.set(molecule.id, molecule);
        
        // Emit events
        this.emit('molecule-added', molecule);
        this.emit('molecules-changed');
        
        return true;
    }
    
    /**
     * Create and add a new molecule with specified properties
     * @param {Object} properties - Molecule properties
     * @returns {Molecule|null} The added molecule or null if failed
     */
    createMolecule(properties) {
        const molecule = new Molecule(
            properties.name,
            properties.count,
            properties.diffusionTranslational,
            properties.diffusionRotational,
            properties.membraneAnchored,
            properties.isLipid,
            properties.isImplicitLipid
        );
        
        // Add binding sites if provided
        if (properties.bindingSites) {
            for (const site of properties.bindingSites) {
                molecule.addBindingSite(new BindingSite(
                    site.name,
                    new Vector3(site.x, site.y, site.z),
                    site.stateCount,
                    site.stateNames
                ));
            }
        }
        
        if (this.addMolecule(molecule)) {
            return molecule;
        }
        
        return null;
    }
    
    /**
     * Update an existing molecule
     * @param {string} id - ID of the molecule to update
     * @param {Object} updates - Properties to update
     * @returns {boolean} True if updated successfully
     */
    updateMolecule(id, updates) {
        const molecule = this.getMolecule(id);
        if (!molecule) {
            this.emit('error', `Molecule with ID '${id}' not found`);
            return false;
        }
        
        // Check for name change and verify it doesn't conflict
        if (updates.name && updates.name !== molecule.name) {
            const existingMolecule = this.getMoleculeByName(updates.name);
            if (existingMolecule && existingMolecule.id !== id) {
                this.emit('error', `Molecule with name '${updates.name}' already exists`);
                return false;
            }
            molecule.name = updates.name;
        }
        
        // Update other properties if provided
        if (updates.count !== undefined) molecule.count = updates.count;
        if (updates.diffusionTranslational !== undefined) molecule.diffusionTranslational = updates.diffusionTranslational;
        if (updates.diffusionRotational !== undefined) molecule.diffusionRotational = updates.diffusionRotational;
        if (updates.membraneAnchored !== undefined) molecule.membraneAnchored = updates.membraneAnchored;
        if (updates.isLipid !== undefined) molecule.isLipid = updates.isLipid;
        if (updates.isImplicitLipid !== undefined) molecule.isImplicitLipid = updates.isImplicitLipid;
        
        // Update binding sites if provided
        if (updates.bindingSites) {
            // Replace all binding sites
            molecule.bindingSites = [];
            for (const site of updates.bindingSites) {
                molecule.addBindingSite(new BindingSite(
                    site.name,
                    new Vector3(site.x, site.y, site.z),
                    site.stateCount,
                    site.stateNames
                ));
            }
        }
        
        // Validate updated molecule
        if (!this.validateMolecule(molecule)) {
            return false;
        }
        
        // Emit events
        this.emit('molecule-updated', molecule);
        this.emit('molecules-changed');
        
        return true;
    }
    
    /**
     * Remove a molecule from the collection
     * @param {string} id - ID of the molecule to remove
     * @returns {boolean} True if removed successfully
     */
    removeMolecule(id) {
        if (!this.molecules.has(id)) {
            this.emit('error', `Molecule with ID '${id}' not found`);
            return false;
        }
        
        const molecule = this.molecules.get(id);
        this.molecules.delete(id);
        
        // Clear selection if the removed molecule was selected
        if (this.selectedMoleculeId === id) {
            this.selectedMoleculeId = null;
            this.selectedBindingSiteIndex = -1;
        }
        
        // Emit events
        this.emit('molecule-deleted', id);
        this.emit('molecules-changed');
        
        return true;
    }
    
    /**
     * Add a binding site to a molecule
     * @param {string} moleculeId - ID of the molecule
     * @param {Object} siteProperties - Binding site properties
     * @returns {boolean} True if added successfully
     */
    addBindingSite(moleculeId, siteProperties) {
        const molecule = this.getMolecule(moleculeId);
        if (!molecule) {
            this.emit('error', `Molecule with ID '${moleculeId}' not found`);
            return false;
        }
        
        try {
            molecule.addBindingSiteByProperties(
                siteProperties.name,
                siteProperties.x,
                siteProperties.y,
                siteProperties.z,
                siteProperties.stateCount,
                siteProperties.stateNames
            );
            
            // Emit events
            this.emit('molecule-updated', molecule);
            this.emit('molecules-changed');
            
            return true;
        } catch (error) {
            this.emit('error', error.message);
            return false;
        }
    }
    
    /**
     * Update a binding site in a molecule
     * @param {string} moleculeId - ID of the molecule
     * @param {string} siteName - Name of the binding site to update
     * @param {Object} updates - Properties to update
     * @returns {boolean} True if updated successfully
     */
    updateBindingSite(moleculeId, siteName, updates) {
        const molecule = this.getMolecule(moleculeId);
        if (!molecule) {
            this.emit('error', `Molecule with ID '${moleculeId}' not found`);
            return false;
        }
        
        try {
            const success = molecule.updateBindingSite(siteName, updates);
            if (!success) {
                this.emit('error', `Binding site '${siteName}' not found in molecule '${molecule.name}'`);
                return false;
            }
            
            // Emit events
            this.emit('molecule-updated', molecule);
            this.emit('molecules-changed');
            
            return true;
        } catch (error) {
            this.emit('error', error.message);
            return false;
        }
    }
    
    /**
     * Remove a binding site from a molecule
     * @param {string} moleculeId - ID of the molecule
     * @param {string} siteName - Name of the binding site to remove
     * @returns {boolean} True if removed successfully
     */
    removeBindingSite(moleculeId, siteName) {
        const molecule = this.getMolecule(moleculeId);
        if (!molecule) {
            this.emit('error', `Molecule with ID '${moleculeId}' not found`);
            return false;
        }
        
        const success = molecule.removeBindingSite(siteName);
        if (!success) {
            this.emit('error', `Binding site '${siteName}' not found in molecule '${molecule.name}'`);
            return false;
        }
        
        // Emit events
        this.emit('molecule-updated', molecule);
        this.emit('molecules-changed');
        
        return true;
    }
    
    /**
     * Set the currently selected molecule and binding site
     * @param {string} moleculeId - ID of the molecule to select
     * @param {number} [bindingSiteIndex=-1] - Index of the binding site to select (-1 for none)
     */
    setSelection(moleculeId, bindingSiteIndex = -1) {
        const oldMoleculeId = this.selectedMoleculeId;
        const oldBindingSiteIndex = this.selectedBindingSiteIndex;
        
        this.selectedMoleculeId = moleculeId;
        this.selectedBindingSiteIndex = bindingSiteIndex;
        
        // Only emit event if selection actually changed
        if (oldMoleculeId !== moleculeId || oldBindingSiteIndex !== bindingSiteIndex) {
            this.emit('selection-changed', {
                moleculeId,
                bindingSiteIndex,
                molecule: moleculeId ? this.getMolecule(moleculeId) : null,
                bindingSite: moleculeId && bindingSiteIndex >= 0 ? 
                    this.getMolecule(moleculeId)?.getBindingSiteByIndex(bindingSiteIndex) : null
            });
        }
    }
    
    /**
     * Get the currently selected molecule and binding site
     * @returns {Object} Object with moleculeId, bindingSiteIndex, molecule, and bindingSite
     */
    getSelection() {
        return {
            moleculeId: this.selectedMoleculeId,
            bindingSiteIndex: this.selectedBindingSiteIndex,
            molecule: this.selectedMoleculeId ? this.getMolecule(this.selectedMoleculeId) : null,
            bindingSite: this.selectedMoleculeId && this.selectedBindingSiteIndex >= 0 ? 
                this.getMolecule(this.selectedMoleculeId)?.getBindingSiteByIndex(this.selectedBindingSiteIndex) : null
        };
    }
    
    /**
     * Calculate diffusion coefficients for a molecule based on its radius
     * @param {number} radius - Molecule radius in nm
     * @returns {Object} Object with translational and rotational diffusion coefficients
     */
    calculateDiffusionCoefficients(radius) {
        return Molecule.prototype.calculateDiffusionCoefficients.call(null, radius);
    }
    
    /**
     * Validate a molecule to ensure it has required properties
     * @param {Molecule} molecule - Molecule to validate
     * @returns {boolean} True if valid
     * @private
     */
    validateMolecule(molecule) {
        // Name is required
        if (!molecule.name || molecule.name.trim() === '') {
            this.emit('error', 'Molecule name is required');
            return false;
        }
        
        // Count must be a positive integer
        if (!Number.isInteger(molecule.count) || molecule.count <= 0) {
            this.emit('error', 'Molecule count must be a positive integer');
            return false;
        }
        
        // Diffusion coefficients must be non-negative numbers
        if (isNaN(molecule.diffusionTranslational) || molecule.diffusionTranslational < 0) {
            this.emit('error', 'Translational diffusion coefficient must be a non-negative number');
            return false;
        }
        
        if (isNaN(molecule.diffusionRotational) || molecule.diffusionRotational < 0) {
            this.emit('error', 'Rotational diffusion coefficient must be a non-negative number');
            return false;
        }
        
        // Validate binding site names are unique
        const siteNames = new Set();
        for (const site of molecule.bindingSites) {
            if (siteNames.has(site.name)) {
                this.emit('error', `Duplicate binding site name '${site.name}'`);
                return false;
            }
            siteNames.add(site.name);
        }
        
        return true;
    }
    
    /**
     * Import molecules from MolecX.mol file content
     * @param {string} content - Content of the MolecX.mol file
     * @returns {Molecule|null} The imported molecule or null if import failed
     */
    importFromMolFile(content) {
        try {
            const molecule = Molecule.fromMolFile(content);
            
            if (this.validateMolecule(molecule)) {
                // Check for duplicate name
                if (this.getMoleculeByName(molecule.name)) {
                    // Generate a unique name by appending a number
                    let counter = 1;
                    let uniqueName = `${molecule.name}_${counter}`;
                    
                    while (this.getMoleculeByName(uniqueName)) {
                        counter++;
                        uniqueName = `${molecule.name}_${counter}`;
                    }
                    
                    molecule.name = uniqueName;
                }
                
                // Add the molecule
                this.molecules.set(molecule.id, molecule);
                
                // Emit events
                this.emit('molecule-added', molecule);
                this.emit('molecules-changed');
                
                return molecule;
            }
        } catch (error) {
            this.emit('error', `Error importing molecule: ${error.message}`);
        }
        
        return null;
    }
    
    /**
     * Export a molecule to MolecX.mol file format
     * @param {string} id - ID of the molecule to export
     * @returns {string|null} The generated file content or null if export failed
     */
    exportToMolFile(id) {
        const molecule = this.getMolecule(id);
        if (!molecule) {
            this.emit('error', `Molecule with ID '${id}' not found`);
            return null;
        }
        
        return molecule.toMolFile();
    }
    
    /**
     * Clear all molecules
     */
    clear() {
        this.molecules.clear();
        this.selectedMoleculeId = null;
        this.selectedBindingSiteIndex = -1;
        
        // Emit event
        this.emit('molecules-changed');
    }
    
    /**
     * Serialize the molecules to JSON
     * @returns {string} JSON string of all molecules
     */
    toJSON() {
        const moleculesArray = Array.from(this.molecules.values()).map(m => m.toObject());
        return JSON.stringify(moleculesArray);
    }
    
    /**
     * Load molecules from JSON
     * @param {string} json - JSON string of molecules
     * @returns {boolean} True if loaded successfully
     */
    fromJSON(json) {
        try {
            const moleculesArray = JSON.parse(json);
            
            // Clear existing molecules
            this.clear();
            
            // Add molecules from JSON
            for (const obj of moleculesArray) {
                const molecule = Molecule.fromObject(obj);
                if (this.validateMolecule(molecule)) {
                    this.molecules.set(molecule.id, molecule);
                }
            }
            
            // Emit event
            this.emit('molecules-changed');
            
            return true;
        } catch (error) {
            this.emit('error', `Error loading molecules from JSON: ${error.message}`);
            return false;
        }
    }
}

export { MoleculeManager };