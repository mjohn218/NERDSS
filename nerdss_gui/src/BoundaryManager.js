/**
 * BoundaryManager.js
 * 
 * This file defines the BoundaryManager class, which manages the simulation
 * boundaries and boundary conditions in the NERDSS application. This includes
 * box dimensions, boundary types, and related settings.
 */

import { EventEmitter } from './utils/EventEmitter.js';

/**
 * Boundary type enum
 * @readonly
 * @enum {string}
 */
const BoundaryType = {
    PERIODIC: 'periodic',
    REFLECTIVE: 'reflect'
};

/**
 * Manages the simulation boundaries and boundary conditions
 * @extends EventEmitter
 */
class BoundaryManager extends EventEmitter {
    /**
     * Create a new BoundaryManager with default values
     */
    constructor() {
        super();
        
        // Initialize with default boundary settings
        this.boundaries = {
            // Box dimensions (nm)
            boxSizeX: 20.0,
            boxSizeY: 20.0,
            boxSizeZ: 20.0,
            
            // Boundary types
            boundaryTypeX: BoundaryType.PERIODIC,
            boundaryTypeY: BoundaryType.PERIODIC,
            boundaryTypeZ: BoundaryType.PERIODIC,
            
            // Lipid settings
            isLipid: false,
            isImplicitLipid: false
        };
    }
    
    /**
     * Get all boundary settings
     * @returns {Object} All boundary settings
     */
    getAllBoundaries() {
        return { ...this.boundaries };
    }
    
    /**
     * Get a specific boundary setting
     * @param {string} name - Name of the setting
     * @returns {*} Setting value or undefined if setting doesn't exist
     */
    getBoundary(name) {
        return this.boundaries[name];
    }
    
    /**
     * Set a boundary setting value
     * @param {string} name - Name of the setting
     * @param {*} value - New value for the setting
     * @returns {boolean} True if the setting was set successfully
     */
    setBoundary(name, value) {
        // Validate setting
        if (!this.validateBoundary(name, value)) {
            return false;
        }
        
        // Update setting
        this.boundaries[name] = value;
        
        // Emit events
        this.emit('boundary-changed', { name, value });
        this.emit('boundaries-changed');
        
        return true;
    }
    
    /**
     * Update multiple boundary settings at once
     * @param {Object} updates - Object with setting names as keys and new values as values
     * @returns {boolean} True if all settings were updated successfully
     */
    updateBoundaries(updates) {
        let allUpdatesSuccessful = true;
        const changedBoundaries = {};
        
        // Validate all settings first
        for (const [name, value] of Object.entries(updates)) {
            if (!this.validateBoundary(name, value)) {
                allUpdatesSuccessful = false;
            }
        }
        
        // If any validation failed, don't update any settings
        if (!allUpdatesSuccessful) {
            return false;
        }
        
        // Update all settings
        for (const [name, value] of Object.entries(updates)) {
            this.boundaries[name] = value;
            changedBoundaries[name] = value;
        }
        
        // Emit events
        for (const [name, value] of Object.entries(changedBoundaries)) {
            this.emit('boundary-changed', { name, value });
        }
        
        this.emit('boundaries-changed');
        
        return true;
    }
    
    /**
     * Validate a boundary setting value
     * @param {string} name - Name of the setting
     * @param {*} value - Value to validate
     * @returns {boolean} True if the value is valid for the setting
     * @private
     */
    validateBoundary(name, value) {
        switch (name) {
            case 'boxSizeX':
            case 'boxSizeY':
            case 'boxSizeZ':
                if (typeof value !== 'number' || isNaN(value) || value <= 0) {
                    this.emit('error', `Box size ${name.slice(-1)} must be a positive number`);
                    return false;
                }
                break;
                
            case 'boundaryTypeX':
            case 'boundaryTypeY':
            case 'boundaryTypeZ':
                if (value !== BoundaryType.PERIODIC && value !== BoundaryType.REFLECTIVE) {
                    this.emit('error', `Boundary type ${name.slice(-1)} must be either 'periodic' or 'reflect'`);
                    return false;
                }
                break;
                
            case 'isLipid':
            case 'isImplicitLipid':
                if (typeof value !== 'boolean') {
                    this.emit('error', `${name} must be a boolean`);
                    return false;
                }
                break;
                
            default:
                this.emit('error', `Unknown boundary setting: ${name}`);
                return false;
        }
        
        return true;
    }
    
    /**
     * Reset all boundary settings to default values
     */
    resetToDefaults() {
        this.boundaries = {
            boxSizeX: 20.0,
            boxSizeY: 20.0,
            boxSizeZ: 20.0,
            boundaryTypeX: BoundaryType.PERIODIC,
            boundaryTypeY: BoundaryType.PERIODIC,
            boundaryTypeZ: BoundaryType.PERIODIC,
            isLipid: false,
            isImplicitLipid: false
        };
        
        // Emit event
        this.emit('boundaries-changed');
    }
    
    /**
     * Set box dimensions
     * @param {number} x - X dimension (nm)
     * @param {number} y - Y dimension (nm)
     * @param {number} z - Z dimension (nm)
     * @returns {boolean} True if dimensions were set successfully
     */
    setBoxDimensions(x, y, z) {
        return this.updateBoundaries({
            boxSizeX: x,
            boxSizeY: y,
            boxSizeZ: z
        });
    }
    
    /**
     * Set boundary types
     * @param {string} x - X boundary type ('periodic' or 'reflect')
     * @param {string} y - Y boundary type ('periodic' or 'reflect')
     * @param {string} z - Z boundary type ('periodic' or 'reflect')
     * @returns {boolean} True if boundary types were set successfully
     */
    setBoundaryTypes(x, y, z) {
        return this.updateBoundaries({
            boundaryTypeX: x,
            boundaryTypeY: y,
            boundaryTypeZ: z
        });
    }
    
    /**
     * Get box dimensions
     * @returns {Object} Object with x, y, z properties
     */
    getBoxDimensions() {
        return {
            x: this.boundaries.boxSizeX,
            y: this.boundaries.boxSizeY,
            z: this.boundaries.boxSizeZ
        };
    }
    
    /**
     * Get boundary types
     * @returns {Object} Object with x, y, z properties
     */
    getBoundaryTypes() {
        return {
            x: this.boundaries.boundaryTypeX,
            y: this.boundaries.boundaryTypeY,
            z: this.boundaries.boundaryTypeZ
        };
    }
    
    /**
     * Serialize the boundary settings to JSON
     * @returns {string} JSON string of all boundary settings
     */
    toJSON() {
        return JSON.stringify(this.boundaries);
    }
    
    /**
     * Load boundary settings from JSON
     * @param {string} json - JSON string of boundary settings
     * @returns {boolean} True if loaded successfully
     */
    fromJSON(json) {
        try {
            const boundaries = JSON.parse(json);
            
            // Validate all settings
            for (const [name, value] of Object.entries(boundaries)) {
                if (!this.validateBoundary(name, value)) {
                    return false;
                }
            }
            
            // Update settings
            this.boundaries = { ...boundaries };
            
            // Emit event
            this.emit('boundaries-changed');
            
            return true;
        } catch (error) {
            this.emit('error', `Error loading boundaries from JSON: ${error.message}`);
            return false;
        }
    }
    
    /**
     * Generate model.inp file boundaries section content
     * @returns {string} Model.inp file boundaries section content
     */
    generateModelInputContent() {
        let content = "start boundaries\n";
        
        content += `\tWaterBox = [${this.boundaries.boxSizeX}, ${this.boundaries.boxSizeY}, ${this.boundaries.boxSizeZ}]\n`;
        content += `\tisLipid = ${this.boundaries.isLipid ? 'true' : 'false'}\n`;
        content += `\tisImplicitLipid = ${this.boundaries.isImplicitLipid ? 'true' : 'false'}\n`;
        content += `\txBCtype = ${this.boundaries.boundaryTypeX}\n`;
        content += `\tyBCtype = ${this.boundaries.boundaryTypeY}\n`;
        content += `\tzBCtype = ${this.boundaries.boundaryTypeZ}\n`;
        
        content += "end boundaries\n";
        
        return content;
    }
}

// Export the BoundaryManager class and BoundaryType enum
export { BoundaryManager, BoundaryType };