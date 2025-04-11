/**
 * ReactionManager.js
 * 
 * This file defines the ReactionManager class, which manages the collection
 * of reactions in the NERDSS application. It provides methods for adding,
 * updating, removing, and retrieving reactions.
 */

import { Reaction, ReactionParticipant, Complex } from './models/Reaction.js';
import { EventEmitter } from './utils/EventEmitter.js';
import { Vector3 } from './math/Vector3.js';

/**
 * Manages the collection of reactions in the application
 * @extends EventEmitter
 */
class ReactionManager extends EventEmitter {
    /**
     * Create a new ReactionManager
     * @param {MoleculeManager} moleculeManager - Reference to the molecule manager
     */
    constructor(moleculeManager) {
        super();
        this.reactions = new Map();
        this.moleculeManager = moleculeManager;
        this.selectedReactionId = null;
        
        // Listen for molecule changes
        if (moleculeManager) {
            moleculeManager.on('molecule-deleted', this.handleMoleculeDeleted.bind(this));
            moleculeManager.on('molecule-updated', this.handleMoleculeUpdated.bind(this));
        }
    }
    
    /**
     * Set the molecule manager reference
     * @param {MoleculeManager} moleculeManager - Reference to the molecule manager
     */
    setMoleculeManager(moleculeManager) {
        this.moleculeManager = moleculeManager;
        
        // Add event listeners
        moleculeManager.on('molecule-deleted', this.handleMoleculeDeleted.bind(this));
        moleculeManager.on('molecule-updated', this.handleMoleculeUpdated.bind(this));
    }
    
    /**
     * Get all reactions as an array
     * @returns {Reaction[]} Array of all reactions
     */
    getAllReactions() {
        return Array.from(this.reactions.values());
    }
    
    /**
     * Get a reaction by ID
     * @param {string} id - Reaction ID
     * @returns {Reaction|undefined} The reaction or undefined if not found
     */
    getReaction(id) {
        return this.reactions.get(id);
    }
    
    /**
     * Add a new reaction to the collection
     * @param {Reaction} reaction - Reaction to add
     * @returns {boolean} True if added successfully
     */
    addReaction(reaction) {
        // Validate reaction
        if (!this.validateReaction(reaction)) {
            return false;
        }
        
        // Add reaction
        this.reactions.set(reaction.id, reaction);
        
        // Emit events
        this.emit('reaction-added', reaction);
        this.emit('reactions-changed');
        
        return true;
    }
    
    /**
     * Create and add a new reaction with specified properties
     * @param {Object} properties - Reaction properties
     * @returns {Reaction|null} The added reaction or null if failed
     */
    createReaction(properties) {
        // Create reactants
        const reactants = properties.reactants.map(r => 
            new ReactionParticipant(r.moleculeId, r.bindingSiteName, r.state)
        );
        
        // Create products (either participants or a complex)
        let products;
        if (properties.productType === 'complex') {
            products = new Complex(
                properties.products.participants.map(p => 
                    new ReactionParticipant(p.moleculeId, p.bindingSiteName, p.state)
                ),
                properties.products.bonds
            );
        } else {
            products = properties.products.map(p => 
                new ReactionParticipant(p.moleculeId, p.bindingSiteName, p.state)
            );
        }
        
        // Create the reaction
        const reaction = new Reaction({
            reactants,
            products,
            onRate: properties.onRate,
            onRateUnit: properties.onRateUnit,
            offRate: properties.offRate,
            offRateUnit: properties.offRateUnit,
            reversible: properties.reversible,
            membraneBound: properties.membraneBound,
            geometry: {
                sigma: properties.geometry.sigma,
                angles: properties.geometry.angles,
                norm1: new Vector3(
                    properties.geometry.norm1.x,
                    properties.geometry.norm1.y,
                    properties.geometry.norm1.z
                ),
                norm2: new Vector3(
                    properties.geometry.norm2.x,
                    properties.geometry.norm2.y,
                    properties.geometry.norm2.z
                )
            }
        });
        
        if (this.addReaction(reaction)) {
            return reaction;
        }
        
        return null;
    }
    
    /**
     * Update an existing reaction
     * @param {string} id - ID of the reaction to update
     * @param {Object} updates - Properties to update
     * @returns {boolean} True if updated successfully
     */
    updateReaction(id, updates) {
        const reaction = this.getReaction(id);
        if (!reaction) {
            this.emit('error', `Reaction with ID '${id}' not found`);
            return false;
        }
        
        // Update reactants if provided
        if (updates.reactants) {
            reaction.reactants = updates.reactants.map(r => 
                new ReactionParticipant(r.moleculeId, r.bindingSiteName, r.state)
            );
        }
        
        // Update products if provided
        if (updates.products) {
            if (updates.productType === 'complex') {
                reaction.products = new Complex(
                    updates.products.participants.map(p => 
                        new ReactionParticipant(p.moleculeId, p.bindingSiteName, p.state)
                    ),
                    updates.products.bonds
                );
            } else {
                reaction.products = updates.products.map(p => 
                    new ReactionParticipant(p.moleculeId, p.bindingSiteName, p.state)
                );
            }
        }
        
        // Update rates
        if (updates.onRate !== undefined) reaction.onRate = updates.onRate;
        if (updates.onRateUnit !== undefined) reaction.onRateUnit = updates.onRateUnit;
        if (updates.offRate !== undefined) reaction.offRate = updates.offRate;
        if (updates.offRateUnit !== undefined) reaction.offRateUnit = updates.offRateUnit;
        
        // Update reaction properties
        if (updates.reversible !== undefined) reaction.reversible = updates.reversible;
        if (updates.membraneBound !== undefined) reaction.membraneBound = updates.membraneBound;
        
        // Update geometry
        if (updates.geometry) {
            if (updates.geometry.sigma !== undefined) reaction.geometry.sigma = updates.geometry.sigma;
            
            if (updates.geometry.angles) {
                reaction.geometry.angles = [...updates.geometry.angles];
            }
            
            if (updates.geometry.norm1) {
                reaction.geometry.norm1 = new Vector3(
                    updates.geometry.norm1.x,
                    updates.geometry.norm1.y,
                    updates.geometry.norm1.z
                );
            }
            
            if (updates.geometry.norm2) {
                reaction.geometry.norm2 = new Vector3(
                    updates.geometry.norm2.x,
                    updates.geometry.norm2.y,
                    updates.geometry.norm2.z
                );
            }
        }
        
        // Validate updated reaction
        if (!this.validateReaction(reaction)) {
            return false;
        }
        
        // Emit events
        this.emit('reaction-updated', reaction);
        this.emit('reactions-changed');
        
        return true;
    }
    
    /**
     * Remove a reaction from the collection
     * @param {string} id - ID of the reaction to remove
     * @returns {boolean} True if removed successfully
     */
    removeReaction(id) {
        if (!this.reactions.has(id)) {
            this.emit('error', `Reaction with ID '${id}' not found`);
            return false;
        }
        
        const reaction = this.reactions.get(id);
        this.reactions.delete(id);
        
        // Clear selection if the removed reaction was selected
        if (this.selectedReactionId === id) {
            this.selectedReactionId = null;
        }
        
        // Emit events
        this.emit('reaction-deleted', id);
        this.emit('reactions-changed');
        
        return true;
    }
    
    /**
     * Set the currently selected reaction
     * @param {string} reactionId - ID of the reaction to select
     */
    setSelection(reactionId) {
        const oldReactionId = this.selectedReactionId;
        this.selectedReactionId = reactionId;
        
        // Only emit event if selection actually changed
        if (oldReactionId !== reactionId) {
            this.emit('selection-changed', {
                reactionId,
                reaction: reactionId ? this.getReaction(reactionId) : null
            });
        }
    }
    
    /**
     * Get the currently selected reaction
     * @returns {Object} Object with reactionId and reaction
     */
    getSelection() {
        return {
            reactionId: this.selectedReactionId,
            reaction: this.selectedReactionId ? this.getReaction(this.selectedReactionId) : null
        };
    }
    
    /**
     * Remove all reactions involving a molecule
     * @param {string} moleculeId - ID of the molecule
     */
    removeMoleculeFromReactions(moleculeId) {
        const reactionsToRemove = [];
        
        // Find reactions involving the molecule
        for (const [id, reaction] of this.reactions.entries()) {
            let containsMolecule = false;
            
            // Check reactants
            for (const reactant of reaction.reactants) {
                if (reactant.moleculeId === moleculeId) {
                    containsMolecule = true;
                    break;
                }
            }
            
            // Check products
            if (!containsMolecule && Array.isArray(reaction.products)) {
                for (const product of reaction.products) {
                    if (product.moleculeId === moleculeId) {
                        containsMolecule = true;
                        break;
                    }
                }
            } else if (!containsMolecule && reaction.products instanceof Complex) {
                for (const participant of reaction.products.participants) {
                    if (participant.moleculeId === moleculeId) {
                        containsMolecule = true;
                        break;
                    }
                }
            }
            
            if (containsMolecule) {
                reactionsToRemove.push(id);
            }
        }
        
        // Remove affected reactions
        for (const id of reactionsToRemove) {
            this.removeReaction(id);
        }
    }
    
    /**
     * Handle a molecule deletion event from the MoleculeManager
     * @param {string} moleculeId - ID of the deleted molecule
     * @private
     */
    handleMoleculeDeleted(moleculeId) {
        this.removeMoleculeFromReactions(moleculeId);
    }
    
    /**
     * Handle a molecule update event from the MoleculeManager
     * @param {Molecule} molecule - The updated molecule
     * @private
     */
    handleMoleculeUpdated(molecule) {
        // Currently no special handling needed
        // If we needed to update references to the molecule, we would do it here
    }
    
    /**
     * Validate a reaction to ensure it has required properties and references valid molecules
     * @param {Reaction} reaction - Reaction to validate
     * @returns {boolean} True if valid
     * @private
     */
    validateReaction(reaction) {
        // Validate reactants
        if (!reaction.reactants || reaction.reactants.length === 0) {
            // Allow zero reactants for creation reactions
        } else {
            for (const reactant of reaction.reactants) {
                // Check if molecule exists
                if (this.moleculeManager && !this.moleculeManager.getMolecule(reactant.moleculeId)) {
                    this.emit('error', `Reactant references non-existent molecule ID: ${reactant.moleculeId}`);
                    return false;
                }
            }
        }
        
        // Validate products
        if (!reaction.products) {
            // Allow undefined products for destruction reactions
        } else if (Array.isArray(reaction.products)) {
            for (const product of reaction.products) {
                // Check if molecule exists
                if (this.moleculeManager && !this.moleculeManager.getMolecule(product.moleculeId)) {
                    this.emit('error', `Product references non-existent molecule ID: ${product.moleculeId}`);
                    return false;
                }
            }
        } else if (reaction.products instanceof Complex) {
            for (const participant of reaction.products.participants) {
                // Check if molecule exists
                if (this.moleculeManager && !this.moleculeManager.getMolecule(participant.moleculeId)) {
                    this.emit('error', `Product complex references non-existent molecule ID: ${participant.moleculeId}`);
                    return false;
                }
            }
        }
        
        // Validate kinetic parameters
        if (isNaN(reaction.onRate) || reaction.onRate < 0) {
            this.emit('error', 'On-rate must be a non-negative number');
            return false;
        }
        
        if (reaction.reversible && (isNaN(reaction.offRate) || reaction.offRate < 0)) {
            this.emit('error', 'Off-rate must be a non-negative number for reversible reactions');
            return false;
        }
        
        // Validate geometry
        if (isNaN(reaction.geometry.sigma) || reaction.geometry.sigma <= 0) {
            this.emit('error', 'Sigma must be a positive number');
            return false;
        }
        
        // Validate reaction type consistency
        const reactantCount = reaction.reactants.length;
        const productCount = Array.isArray(reaction.products) ? reaction.products.length : 
                             (reaction.products instanceof Complex ? 1 : 0);
                             
        if (reactantCount === 2 && productCount === 1) {
            // Binding reaction - need to validate angles
            for (let i = 0; i < 5; i++) {
                if (isNaN(reaction.geometry.angles[i])) {
                    this.emit('error', `Invalid angle at index ${i}`);
                    return false;
                }
            }
        }
        
        return true;
    }
    
    /**
     * Clear all reactions
     */
    clear() {
        this.reactions.clear();
        this.selectedReactionId = null;
        
        // Emit event
        this.emit('reactions-changed');
    }
    
    /**
     * Serialize the reactions to JSON
     * @returns {string} JSON string of all reactions
     */
    toJSON() {
        const reactionsArray = Array.from(this.reactions.values()).map(r => r.toObject());
        return JSON.stringify(reactionsArray);
    }
    
    /**
     * Load reactions from JSON
     * @param {string} json - JSON string of reactions
     * @returns {boolean} True if loaded successfully
     */
    fromJSON(json) {
        try {
            const reactionsArray = JSON.parse(json);
            
            // Clear existing reactions
            this.clear();
            
            // Add reactions from JSON
            for (const obj of reactionsArray) {
                const reaction = Reaction.fromObject(obj);
                if (this.validateReaction(reaction)) {
                    this.reactions.set(reaction.id, reaction);
                }
            }
            
            // Emit event
            this.emit('reactions-changed');
            
            return true;
        } catch (error) {
            this.emit('error', `Error loading reactions from JSON: ${error.message}`);
            return false;
        }
    }
    
    /**
     * Generate model.inp file content for all reactions
     * @returns {string} Model.inp file reaction section content
     */
    generateModelInputContent() {
        if (!this.moleculeManager) {
            this.emit('error', 'Cannot generate model input without MoleculeManager reference');
            return '';
        }
        
        // Create a map of molecule IDs to molecule objects
        const moleculeMap = {};
        for (const molecule of this.moleculeManager.getAllMolecules()) {
            moleculeMap[molecule.id] = molecule;
        }
        
        let content = "start reactions\n";
        
        for (const reaction of this.reactions.values()) {
            content += reaction.toModelInputFormat(moleculeMap);
        }
        
        content += "end reactions";
        
        return content;
    }
}

export { ReactionManager };