/**
 * ParameterManager.js
 * 
 * This file defines the ParameterManager class, which manages the simulation
 * parameters in the NERDSS application. This includes time step, number of
 * iterations, output frequencies, and other simulation control parameters.
 */

import { EventEmitter } from './utils/EventEmitter.js';

/**
 * Manages the simulation parameters
 * @extends EventEmitter
 */
class ParameterManager extends EventEmitter {
    /**
     * Create a new ParameterManager with default values
     */
    constructor() {
        super();
        
        // Initialize with default parameters
        this.parameters = {
            // Time parameters
            timeStep: 10.0,        // Simulation time step (μs)
            numIterations: 1000,   // Number of time steps to run
            
            // Output frequencies
            statOutputFreq: 100,   // Frequency to output statistics (time steps)
            configOutputFreq: 100, // Frequency to output configuration (time steps)
            restartOutputFreq: 100, // Frequency to output restart files (time steps)
            
            // Restart parameters
            isRestart: false,      // Whether this is a restart simulation
            
            // Additional parameters
            overlapSepLimit: 0.6,  // Overlap separation limit
        };
    }
    
    /**
     * Get all parameters
     * @returns {Object} All simulation parameters
     */
    getAllParameters() {
        return { ...this.parameters };
    }
    
    /**
     * Get a specific parameter
     * @param {string} name - Name of the parameter
     * @returns {*} Parameter value or undefined if parameter doesn't exist
     */
    getParameter(name) {
        return this.parameters[name];
    }
    
    /**
     * Set a parameter value
     * @param {string} name - Name of the parameter
     * @param {*} value - New value for the parameter
     * @returns {boolean} True if the parameter was set successfully
     */
    setParameter(name, value) {
        // Validate parameter
        if (!this.validateParameter(name, value)) {
            return false;
        }
        
        // Update parameter
        this.parameters[name] = value;
        
        // Emit events
        this.emit('parameter-changed', { name, value });
        this.emit('parameters-changed');
        
        return true;
    }
    
    /**
     * Update multiple parameters at once
     * @param {Object} updates - Object with parameter names as keys and new values as values
     * @returns {boolean} True if all parameters were updated successfully
     */
    updateParameters(updates) {
        let allUpdatesSuccessful = true;
        const changedParameters = {};
        
        // Validate all parameters first
        for (const [name, value] of Object.entries(updates)) {
            if (!this.validateParameter(name, value)) {
                allUpdatesSuccessful = false;
            }
        }
        
        // If any validation failed, don't update any parameters
        if (!allUpdatesSuccessful) {
            return false;
        }
        
        // Update all parameters
        for (const [name, value] of Object.entries(updates)) {
            this.parameters[name] = value;
            changedParameters[name] = value;
        }
        
        // Emit events
        for (const [name, value] of Object.entries(changedParameters)) {
            this.emit('parameter-changed', { name, value });
        }
        
        this.emit('parameters-changed');
        
        return true;
    }
    
    /**
     * Validate a parameter value
     * @param {string} name - Name of the parameter
     * @param {*} value - Value to validate
     * @returns {boolean} True if the value is valid for the parameter
     * @private
     */
    validateParameter(name, value) {
        switch (name) {
            case 'timeStep':
                if (typeof value !== 'number' || isNaN(value) || value <= 0) {
                    this.emit('error', 'Time step must be a positive number');
                    return false;
                }
                break;
                
            case 'numIterations':
                if (!Number.isInteger(value) || value <= 0) {
                    this.emit('error', 'Number of iterations must be a positive integer');
                    return false;
                }
                break;
                
            case 'statOutputFreq':
            case 'configOutputFreq':
            case 'restartOutputFreq':
                if (!Number.isInteger(value) || value < 0) {
                    this.emit('error', 'Output frequency must be a non-negative integer');
                    return false;
                }
                break;
                
            case 'isRestart':
                if (typeof value !== 'boolean') {
                    this.emit('error', 'Restart flag must be a boolean');
                    return false;
                }
                break;
                
            case 'overlapSepLimit':
                if (typeof value !== 'number' || isNaN(value) || value <= 0) {
                    this.emit('error', 'Overlap separation limit must be a positive number');
                    return false;
                }
                break;
                
            default:
                this.emit('error', `Unknown parameter: ${name}`);
                return false;
        }
        
        return true;
    }
    
    /**
     * Reset all parameters to default values
     */
    resetToDefaults() {
        this.parameters = {
            timeStep: 10.0,
            numIterations: 1000,
            statOutputFreq: 100,
            configOutputFreq: 100,
            restartOutputFreq: 100,
            isRestart: false,
            overlapSepLimit: 0.6,
        };
        
        // Emit event
        this.emit('parameters-changed');
    }
    
    /**
     * Serialize the parameters to JSON
     * @returns {string} JSON string of all parameters
     */
    toJSON() {
        return JSON.stringify(this.parameters);
    }
    
    /**
     * Load parameters from JSON
     * @param {string} json - JSON string of parameters
     * @returns {boolean} True if loaded successfully
     */
    fromJSON(json) {
        try {
            const parameters = JSON.parse(json);
            
            // Validate all parameters
            for (const [name, value] of Object.entries(parameters)) {
                if (!this.validateParameter(name, value)) {
                    return false;
                }
            }
            
            // Update parameters
            this.parameters = { ...parameters };
            
            // Emit event
            this.emit('parameters-changed');
            
            return true;
        } catch (error) {
            this.emit('error', `Error loading parameters from JSON: ${error.message}`);
            return false;
        }
    }
    
    /**
     * Generate model.inp file parameters section content
     * @returns {string} Model.inp file parameters section content
     */
    generateModelInputContent() {
        let content = "start parameters\n";
        
        content += `\tnItr = ${this.parameters.numIterations}\n`;
        content += `\ttimeStep = ${this.parameters.timeStep}\n\n`;
        
        content += `\ttimeWrite = ${this.parameters.configOutputFreq}\n`;
        content += `\ttrajWrite = ${this.parameters.statOutputFreq}\n`;
        content += `\trestartWrite = ${this.parameters.restartOutputFreq}\n`;
        content += `\tfromRestart = ${this.parameters.isRestart ? 'true' : 'false'}\n`;
        content += `\toverlapSepLimit = ${this.parameters.overlapSepLimit}\n`;
        
        content += "end parameters\n";
        
        return content;
    }
}

export { ParameterManager };