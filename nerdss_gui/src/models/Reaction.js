/**
 * Reaction.js
 * 
 * This file defines the Reaction class, which represents a chemical reaction
 * between molecules in the NERDSS simulation.
 * 
 * A reaction has:
 * - Reactants and products (molecules and their binding sites)
 * - Kinetic parameters (on-rate, off-rate)
 * - Geometric parameters (sigma, angles)
 * - Properties related to reaction dimension (2D/3D)
 */

import { Vector3 } from '../math/Vector3.js';
import { Quaternion } from '../math/Quaternion.js';

/**
 * Represents a reaction participant (reactant or product)
 */
class ReactionParticipant {
    /**
     * Create a reaction participant
     * @param {string} moleculeId - ID of the molecule
     * @param {string} bindingSiteName - Name of the binding site
     * @param {string} [state] - Optional state of the binding site
     */
    constructor(moleculeId, bindingSiteName, state = null) {
        this.moleculeId = moleculeId;
        this.bindingSiteName = bindingSiteName;
        this.state = state;
    }
    
    /**
     * Create a deep copy of this participant
     * @returns {ReactionParticipant} New reaction participant with same values
     */
    clone() {
        return new ReactionParticipant(
            this.moleculeId,
            this.bindingSiteName,
            this.state
        );
    }
    
    /**
     * Convert to a simple object for serialization
     * @returns {Object} Plain object representation
     */
    toObject() {
        return {
            moleculeId: this.moleculeId,
            bindingSiteName: this.bindingSiteName,
            state: this.state
        };
    }
    
    /**
     * Create from a plain object
     * @param {Object} obj - Object with participant properties
     * @returns {ReactionParticipant} New reaction participant
     */
    static fromObject(obj) {
        return new ReactionParticipant(
            obj.moleculeId,
            obj.bindingSiteName,
            obj.state
        );
    }
}

/**
 * Represents a complex (bound molecules) in a reaction
 */
class Complex {
    /**
     * Create a reaction complex
     * @param {ReactionParticipant[]} participants - Participants in this complex
     * @param {Object[]} bonds - Bonds between participants (indices to participants array)
     */
    constructor(participants = [], bonds = []) {
        this.participants = participants;
        this.bonds = bonds;
    }
    
    /**
     * Add a participant to this complex
     * @param {ReactionParticipant} participant - Participant to add
     * @returns {Complex} This complex for chaining
     */
    addParticipant(participant) {
        this.participants.push(participant);
        return this;
    }
    
    /**
     * Add a bond between participants
     * @param {number} index1 - Index of first participant
     * @param {number} index2 - Index of second participant
     * @returns {Complex} This complex for chaining
     */
    addBond(index1, index2) {
        this.bonds.push({ from: index1, to: index2 });
        return this;
    }
    
    /**
     * Create a deep copy of this complex
     * @returns {Complex} New complex with same values
     */
    clone() {
        return new Complex(
            this.participants.map(p => p.clone()),
            [...this.bonds]
        );
    }
    
    /**
     * Convert to a simple object for serialization
     * @returns {Object} Plain object representation
     */
    toObject() {
        return {
            participants: this.participants.map(p => p.toObject()),
            bonds: this.bonds
        };
    }
    
    /**
     * Create from a plain object
     * @param {Object} obj - Object with complex properties
     * @returns {Complex} New complex
     */
    static fromObject(obj) {
        return new Complex(
            obj.participants.map(p => ReactionParticipant.fromObject(p)),
            obj.bonds
        );
    }
}

/**
 * Represents a reaction in the simulation
 */
class Reaction {
    /**
     * Create a new reaction
     * @param {Object} options - Configuration options
     * @param {ReactionParticipant[]} options.reactants - Reactants (for unimolecular/bimolecular)
     * @param {(ReactionParticipant[]|Complex[])} options.products - Products
     * @param {number} options.onRate - Forward reaction rate
     * @param {string} options.onRateUnit - Unit for the forward rate
     * @param {number} options.offRate - Backward reaction rate (for reversible reactions)
     * @param {string} options.offRateUnit - Unit for the backward rate
     * @param {boolean} options.reversible - Whether the reaction is reversible
     * @param {boolean} options.membraneBound - Whether the reaction occurs on a membrane (2D)
     * @param {Object} options.geometry - Geometric parameters
     * @param {number} options.geometry.sigma - Reaction radius (nm)
     * @param {number[]} options.geometry.angles - Array of angles [theta1, theta2, phi1, phi2, omega]
     * @param {Vector3} options.geometry.norm1 - Normal vector for first reactant
     * @param {Vector3} options.geometry.norm2 - Normal vector for second reactant
     */
    constructor(options) {
        this.id = Date.now().toString(36) + Math.random().toString(36).substring(2);
        
        // Reaction participants
        this.reactants = options.reactants || [];
        this.products = options.products || [];
        
        // Kinetic parameters
        this.onRate = options.onRate || 0;
        this.onRateUnit = options.onRateUnit || "M⁻¹s⁻¹";
        this.offRate = options.offRate || 0;
        this.offRateUnit = options.offRateUnit || "s⁻¹";
        
        // Reaction properties
        this.reversible = options.reversible || false;
        this.membraneBound = options.membraneBound || false;
        
        // Geometric parameters
        this.geometry = {
            sigma: options.geometry?.sigma || 1.0,
            angles: options.geometry?.angles || [0, 0, 0, 0, 0],
            norm1: options.geometry?.norm1 || new Vector3(0, 0, 1),
            norm2: options.geometry?.norm2 || new Vector3(0, 0, 1)
        };
    }
    
    /**
     * Get the number of reactants
     * @returns {number} Number of reactants
     */
    getReactantCount() {
        return this.reactants.length;
    }
    
    /**
     * Get the number of products
     * @returns {number} Number of products
     */
    getProductCount() {
        return this.products.length;
    }
    
    /**
     * Get the reaction type based on reactant and product counts
     * @returns {string} Reaction type (binding, dissociation, catalysis, etc.)
     */
    getReactionType() {
        const reactantCount = this.getReactantCount();
        const productCount = this.getProductCount();
        
        if (reactantCount === 0 && productCount > 0) {
            return "creation";
        } else if (reactantCount > 0 && productCount === 0) {
            return "destruction";
        } else if (reactantCount === 1 && productCount === 1) {
            return "conversion";
        } else if (reactantCount === 2 && productCount === 1) {
            return "binding";
        } else if (reactantCount === 1 && productCount === 2) {
            return "dissociation";
        } else if (reactantCount === 2 && productCount === 2) {
            return "catalysis";
        } else {
            return "custom";
        }
    }
    
    /**
     * Get the correct rate units based on reaction type and dimensionality
     * @returns {Object} Object with onRateUnit and offRateUnit
     */
    getCorrectRateUnits() {
        const reactantCount = this.getReactantCount();
        const is2D = this.membraneBound;
        
        let onRateUnit, offRateUnit;
        
        // On-rate units depend on reaction order and dimensionality
        if (reactantCount === 0) {
            onRateUnit = is2D ? "nm⁻²μs⁻¹" : "μM·s⁻¹";
        } else if (reactantCount === 1) {
            onRateUnit = "s⁻¹";
        } else if (reactantCount === 2) {
            onRateUnit = is2D ? "nm²·μs⁻¹" : "μM⁻¹·s⁻¹";
        }
        
        // Off-rate units (for reversible reactions)
        if (this.reversible) {
            const productCount = this.getProductCount();
            
            if (productCount === 0) {
                offRateUnit = is2D ? "nm⁻²μs⁻¹" : "μM·s⁻¹";
            } else if (productCount === 1) {
                offRateUnit = "s⁻¹";
            } else if (productCount === 2) {
                offRateUnit = is2D ? "nm²·μs⁻¹" : "μM⁻¹·s⁻¹";
            }
        } else {
            offRateUnit = "s⁻¹";
        }
        
        return { onRateUnit, offRateUnit };
    }
    
    /**
     * Update rate units based on reaction properties
     * This should be called after changing reactant/product counts or dimensionality
     */
    updateRateUnits() {
        const { onRateUnit, offRateUnit } = this.getCorrectRateUnits();
        this.onRateUnit = onRateUnit;
        this.offRateUnit = offRateUnit;
    }
    
    /**
     * Set the reaction angles
     * @param {number} theta1 - Theta angle for first reactant (degrees)
     * @param {number} theta2 - Theta angle for second reactant (degrees)
     * @param {number} phi1 - Phi angle for first reactant (degrees)
     * @param {number} phi2 - Phi angle for second reactant (degrees)
     * @param {number} omega - Omega angle (degrees)
     * @returns {Reaction} This reaction for chaining
     */
    setAngles(theta1, theta2, phi1, phi2, omega) {
        this.geometry.angles = [
            theta1, 
            theta2, 
            phi1, 
            phi2, 
            omega
        ];
        return this;
    }
    
    /**
     * Convert angles from degrees to radians
     * @returns {number[]} Array of angles in radians
     */
    getAnglesInRadians() {
        return this.geometry.angles.map(angle => angle * Math.PI / 180);
    }
    
    /**
     * Update normal vectors from angles
     * This calculates appropriate normal vectors based on current angles
     */
    updateNormals() {
        const [theta1, theta2, phi1, phi2, omega] = this.getAnglesInRadians();
        
        // Create rotation quaternions
        const qTheta1 = new Quaternion().setFromAxisAngle(new Vector3(0, 1, 0), theta1);
        const qPhi1 = new Quaternion().setFromAxisAngle(new Vector3(1, 0, 0), phi1);
        const qTheta2 = new Quaternion().setFromAxisAngle(new Vector3(0, 1, 0), theta2);
        const qPhi2 = new Quaternion().setFromAxisAngle(new Vector3(1, 0, 0), phi2);
        
        // Apply rotations to initial normal vector
        const initialNormal = new Vector3(0, 0, 1);
        
        // Calculate norm1
        this.geometry.norm1 = initialNormal.clone();
        qTheta1.multiplyVector3(this.geometry.norm1);
        qPhi1.multiplyVector3(this.geometry.norm1);
        
        // Calculate norm2
        this.geometry.norm2 = initialNormal.clone();
        qTheta2.multiplyVector3(this.geometry.norm2);
        qPhi2.multiplyVector3(this.geometry.norm2);
        
        // Apply omega rotation to norm2
        const qOmega = new Quaternion().setFromAxisAngle(new Vector3(0, 0, 1), omega);
        qOmega.multiplyVector3(this.geometry.norm2);
    }
    
    /**
     * Format the reaction as a BNGL string
     * @param {Object} moleculeMap - Map of molecule IDs to molecule objects
     * @returns {string} BNGL-formatted reaction string
     */
    toBNGLString(moleculeMap) {
        // Format reactants
        let reactantsStr = this.reactants.map(reactant => {
            const molecule = moleculeMap[reactant.moleculeId];
            if (!molecule) return "UnknownMolecule";
            
            return molecule.toBNGL(
                reactant.bindingSiteName,
                reactant.state
            );
        }).join(" + ");
        
        // Format products for simple cases
        let productsStr = "";
        if (Array.isArray(this.products) && this.products.length > 0 && 
            this.products[0] instanceof ReactionParticipant) {
            // Simple products (dissociated molecules)
            productsStr = this.products.map(product => {
                const molecule = moleculeMap[product.moleculeId];
                if (!molecule) return "UnknownMolecule";
                
                return molecule.toBNGL(
                    product.bindingSiteName,
                    product.state
                );
            }).join(" + ");
        } else if (this.getReactionType() === "binding") {
            // For binding reactions, create a complex string
            const mol1 = moleculeMap[this.reactants[0].moleculeId];
            const mol2 = moleculeMap[this.reactants[1].moleculeId];
            
            if (mol1 && mol2) {
                productsStr = `${mol1.toBNGL(this.reactants[0].bindingSiteName, null, "1")}.${mol2.toBNGL(this.reactants[1].bindingSiteName, null, "1")}`;
            } else {
                productsStr = "UnknownComplex";
            }
        }
        
        // Create the full reaction string
        let reactionStr = `${reactantsStr} ${this.reversible ? "<->" : "->"} ${productsStr}`;
        
        return reactionStr;
    }
    
    /**
     * Check if this reaction is identical to another
     * @param {Reaction} other - The other reaction to compare with
     * @returns {boolean} True if the reactions are identical
     */
    isIdenticalTo(other) {
        // Check basic properties
        if (this.reversible !== other.reversible ||
            this.membraneBound !== other.membraneBound ||
            this.reactants.length !== other.reactants.length ||
            this.products.length !== other.products.length ||
            this.onRate !== other.onRate ||
            this.offRate !== other.offRate ||
            this.geometry.sigma !== other.geometry.sigma) {
            return false;
        }
        
        // Check angles
        for (let i = 0; i < 5; i++) {
            if (this.geometry.angles[i] !== other.geometry.angles[i]) {
                return false;
            }
        }
        
        // Check reactants and products
        // This is a simplification - full implementation would need deeper comparison
        // of reactants/products and their properties
        
        return true;
    }
    
    /**
     * Create a deep copy of this reaction
     * @returns {Reaction} New reaction with same values
     */
    clone() {
        const clone = new Reaction({
            reactants: this.reactants.map(r => r.clone()),
            products: Array.isArray(this.products) ? 
                this.products.map(p => p.clone()) : 
                this.products.clone(),
            onRate: this.onRate,
            onRateUnit: this.onRateUnit,
            offRate: this.offRate,
            offRateUnit: this.offRateUnit,
            reversible: this.reversible,
            membraneBound: this.membraneBound,
            geometry: {
                sigma: this.geometry.sigma,
                angles: [...this.geometry.angles],
                norm1: this.geometry.norm1.clone(),
                norm2: this.geometry.norm2.clone()
            }
        });
        
        clone.id = this.id;
        return clone;
    }
    
    /**
     * Convert to a simple object for serialization
     * @returns {Object} Plain object representation
     */
    toObject() {
        return {
            id: this.id,
            reactants: this.reactants.map(r => r.toObject()),
            products: Array.isArray(this.products) ? 
                this.products.map(p => p.toObject()) : 
                this.products.toObject(),
            onRate: this.onRate,
            onRateUnit: this.onRateUnit,
            offRate: this.offRate,
            offRateUnit: this.offRateUnit,
            reversible: this.reversible,
            membraneBound: this.membraneBound,
            geometry: {
                sigma: this.geometry.sigma,
                angles: [...this.geometry.angles],
                norm1: {
                    x: this.geometry.norm1.x,
                    y: this.geometry.norm1.y,
                    z: this.geometry.norm1.z
                },
                norm2: {
                    x: this.geometry.norm2.x,
                    y: this.geometry.norm2.y,
                    z: this.geometry.norm2.z
                }
            }
        };
    }
    
    /**
     * Create from a plain object
     * @param {Object} obj - Object with reaction properties
     * @returns {Reaction} New reaction
     */
    static fromObject(obj) {
        const reaction = new Reaction({
            reactants: obj.reactants.map(r => ReactionParticipant.fromObject(r)),
            products: Array.isArray(obj.products) ? 
                obj.products.map(p => ReactionParticipant.fromObject(p)) : 
                Complex.fromObject(obj.products),
            onRate: obj.onRate,
            onRateUnit: obj.onRateUnit,
            offRate: obj.offRate,
            offRateUnit: obj.offRateUnit,
            reversible: obj.reversible,
            membraneBound: obj.membraneBound,
            geometry: {
                sigma: obj.geometry.sigma,
                angles: obj.geometry.angles,
                norm1: new Vector3(
                    obj.geometry.norm1.x,
                    obj.geometry.norm1.y,
                    obj.geometry.norm1.z
                ),
                norm2: new Vector3(
                    obj.geometry.norm2.x,
                    obj.geometry.norm2.y,
                    obj.geometry.norm2.z
                )
            }
        });
        
        if (obj.id) {
            reaction.id = obj.id;
        }
        
        return reaction;
    }
    
    /**
     * Export reaction to NERDSS model.inp file format
     * @param {Object} moleculeMap - Map of molecule IDs to molecule objects
     * @returns {string} Model.inp formatted reaction section
     */
    toModelInputFormat(moleculeMap) {
        const reactionStr = this.toBNGLString(moleculeMap);
        let content = `\t${reactionStr}\n`;
        content += `\tonRate = ${this.onRate}\n`;
        content += `\toffRate = ${this.offRate}\n`;
        content += `\tsigma = ${this.geometry.sigma}\n`;
        
        const norm1Str = `[${this.geometry.norm1.x.toFixed(7)}, ${this.geometry.norm1.y.toFixed(7)}, ${this.geometry.norm1.z.toFixed(7)}]`;
        const norm2Str = `[${this.geometry.norm2.x.toFixed(7)}, ${this.geometry.norm2.y.toFixed(7)}, ${this.geometry.norm2.z.toFixed(7)}]`;
        
        content += `\tnorm1 = ${norm1Str}\n`;
        content += `\tnorm2 = ${norm2Str}\n`;
        content += "\tbindRadSameCom = 1.5 #scales sigma to define distance\n";
        content += "\tloopCoopFactor = 0.001\n";
        content += `\tlength3Dto2D = ${2 * this.geometry.sigma} #default 2*sigma\n`;
        
        // Format the angles for the file
        const [theta1, theta2, phi1, phi2, omega] = this.geometry.angles;
        
        // Handle -999 (undefined) values for phi angles
        const phi1Str = phi1 === -999 ? "-" : theta1.toString();
        const phi2Str = phi2 === -999 ? "-" : theta2.toString();
        
        content += `\tassocAngles = [${theta1}, ${theta2}, ${phi1Str}, ${phi2Str}, ${omega}]\n\n`;
        
        return content;
    }
}

// Export classes
export { Reaction, ReactionParticipant, Complex };