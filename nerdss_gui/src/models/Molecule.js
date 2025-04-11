/**
 * Molecule.js
 * 
 * This file defines the Molecule class, which represents a molecular structure
 * with binding sites (interfaces) in the NERDSS simulation.
 * 
 * A molecule has:
 * - Basic properties (name, count, diffusion coefficients)
 * - Binding sites (interfaces) with positions and states
 * - Membrane anchoring properties
 */

import { Vector3 } from '../math/Vector3.js';

/**
 * Represents a binding site (interface) on a molecule
 */
class BindingSite {
    /**
     * Create a new binding site
     * @param {string} name - Name of the binding site
     * @param {Vector3} position - Position relative to molecule center
     * @param {number} stateCount - Number of states this binding site can have
     * @param {string[]} stateNames - Names of the possible states
     */
    constructor(name, position, stateCount = 0, stateNames = []) {
        this.name = name;
        this.position = position;
        this.stateCount = stateCount;
        this.stateNames = stateNames.length > 0 ? stateNames : ['-'];
        
        // Calculate normal vector for the binding site (for phi calculations)
        this.normal = this.calculateNormal();
    }
    
    /**
     * Calculate the normal vector for this binding site
     * Used for reaction orientation calculations
     * @returns {Vector3} Normal vector
     */
    calculateNormal() {
        // Create a reference Z vector
        const zAxis = new Vector3(0, 0, 1);
        
        // Copy the position vector to avoid modifying it
        const positionVector = this.position.clone();
        
        // Get normal by cross product with z-axis
        const normal = positionVector.cross(zAxis);
        
        // If the normal has zero length (site on z-axis), use x-axis
        if (normal.lengthSquared() < 0.0001) {
            return new Vector3(1, 0, 0);
        }
        
        return normal.normalize();
    }
    
    /**
     * Create a string representation of the binding site states
     * @returns {string} Comma-separated list of state names
     */
    getStateNamesString() {
        return this.stateNames.join(', ');
    }
    
    /**
     * Clone this binding site
     * @returns {BindingSite} A new binding site with the same properties
     */
    clone() {
        return new BindingSite(
            this.name,
            this.position.clone(),
            this.stateCount,
            [...this.stateNames]
        );
    }
    
    /**
     * Convert the binding site to a plain object for serialization
     * @returns {Object} Plain object representation
     */
    toObject() {
        return {
            name: this.name,
            position: {
                x: this.position.x,
                y: this.position.y,
                z: this.position.z
            },
            stateCount: this.stateCount,
            stateNames: [...this.stateNames]
        };
    }
    
    /**
     * Create a binding site from a plain object
     * @param {Object} obj - Object with binding site properties
     * @returns {BindingSite} New binding site
     */
    static fromObject(obj) {
        return new BindingSite(
            obj.name,
            new Vector3(obj.position.x, obj.position.y, obj.position.z),
            obj.stateCount,
            obj.stateNames
        );
    }
}

/**
 * Represents a molecule in the simulation
 */
class Molecule {
    /**
     * Create a new molecule
     * @param {string} name - Molecule name
     * @param {number} count - Number of copies in simulation
     * @param {number} diffusionTranslational - Translational diffusion coefficient
     * @param {number} diffusionRotational - Rotational diffusion coefficient
     * @param {boolean} membraneAnchored - Whether the molecule is anchored to a membrane
     * @param {boolean} isLipid - Whether the molecule is a lipid
     * @param {boolean} isImplicitLipid - Whether the molecule is an implicit lipid
     */
    constructor(name, count, diffusionTranslational, diffusionRotational, 
                membraneAnchored = false, isLipid = false, isImplicitLipid = false) {
        this.id = Date.now().toString(36) + Math.random().toString(36).substring(2);
        this.name = name;
        this.count = count;
        this.diffusionTranslational = diffusionTranslational;
        this.diffusionRotational = diffusionRotational;
        this.membraneAnchored = membraneAnchored;
        this.isLipid = isLipid;
        this.isImplicitLipid = isImplicitLipid;
        this.bindingSites = [];
        
        // Center of mass is always at the origin (0,0,0)
        this.centerOfMass = new Vector3(0, 0, 0);
    }
    
    /**
     * Add a binding site to the molecule
     * @param {BindingSite} bindingSite - The binding site to add
     * @returns {Molecule} This molecule for chaining
     */
    addBindingSite(bindingSite) {
        // Check for duplicate name
        if (this.bindingSites.some(site => site.name === bindingSite.name)) {
            throw new Error(`Binding site with name '${bindingSite.name}' already exists`);
        }
        
        this.bindingSites.push(bindingSite);
        return this;
    }
    
    /**
     * Add a binding site to the molecule by providing individual properties
     * @param {string} name - Binding site name
     * @param {number} x - X coordinate
     * @param {number} y - Y coordinate
     * @param {number} z - Z coordinate
     * @param {number} stateCount - Number of states
     * @param {string[]} stateNames - Names of the states
     * @returns {Molecule} This molecule for chaining
     */
    addBindingSiteByProperties(name, x, y, z, stateCount = 0, stateNames = []) {
        const position = new Vector3(x, y, z);
        const bindingSite = new BindingSite(name, position, stateCount, stateNames);
        return this.addBindingSite(bindingSite);
    }
    
    /**
     * Update an existing binding site
     * @param {string} name - Name of the binding site to update
     * @param {Object} updates - Properties to update
     * @returns {boolean} True if the binding site was found and updated
     */
    updateBindingSite(name, updates) {
        const index = this.bindingSites.findIndex(site => site.name === name);
        if (index === -1) {
            return false;
        }
        
        const site = this.bindingSites[index];
        
        if (updates.name !== undefined && updates.name !== site.name) {
            // Check for duplicate name
            if (this.bindingSites.some(s => s.name === updates.name)) {
                throw new Error(`Binding site with name '${updates.name}' already exists`);
            }
            site.name = updates.name;
        }
        
        if (updates.position !== undefined) {
            site.position.copy(updates.position);
            site.normal = site.calculateNormal(); // Recalculate normal
        } else if (updates.x !== undefined || updates.y !== undefined || updates.z !== undefined) {
            if (updates.x !== undefined) site.position.x = updates.x;
            if (updates.y !== undefined) site.position.y = updates.y;
            if (updates.z !== undefined) site.position.z = updates.z;
            site.normal = site.calculateNormal(); // Recalculate normal
        }
        
        if (updates.stateCount !== undefined) {
            site.stateCount = updates.stateCount;
        }
        
        if (updates.stateNames !== undefined) {
            site.stateNames = updates.stateNames;
        }
        
        return true;
    }
    
    /**
     * Remove a binding site from the molecule
     * @param {string} name - Name of the binding site to remove
     * @returns {boolean} True if the binding site was found and removed
     */
    removeBindingSite(name) {
        const initialLength = this.bindingSites.length;
        this.bindingSites = this.bindingSites.filter(site => site.name !== name);
        return this.bindingSites.length < initialLength;
    }
    
    /**
     * Get a binding site by name
     * @param {string} name - Name of the binding site
     * @returns {BindingSite|null} The binding site or null if not found
     */
    getBindingSite(name) {
        return this.bindingSites.find(site => site.name === name) || null;
    }
    
    /**
     * Get a binding site by index
     * @param {number} index - Index of the binding site
     * @returns {BindingSite|null} The binding site or null if index is out of bounds
     */
    getBindingSiteByIndex(index) {
        return (index >= 0 && index < this.bindingSites.length) 
            ? this.bindingSites[index] 
            : null;
    }
    
    /**
     * Calculate diffusion coefficients based on molecule radius
     * @param {number} radius - Molecule radius in nm
     * @returns {Object} Object with translational and rotational diffusion coefficients
     */
    calculateDiffusionCoefficients(radius) {
        if (radius <= 0) {
            return { 
                translational: 0, 
                rotational: 0 
            };
        }
        
        // Constants for diffusion calculation
        const kB = 1.3806e-23;  // Boltzmann constant
        const T = 298;          // Temperature (Kelvin)
        const eta = 0.00089;    // Viscosity (Pa·s)
        
        // Calculate translational diffusion (um²/s)
        const translational = kB * T / (6 * Math.PI * eta * radius * 1e-9) * 1e12;
        
        // Calculate rotational diffusion (1/us)
        const rotational = 1e-6 * kB * T / (8 * Math.PI * eta * Math.pow(radius, 3) * 1e-27);
        
        return {
            translational,
            rotational
        };
    }
    
    /**
     * Get the total count of elements (1 + #interfaces + #states)
     * Used for traversing the molecule in tree views
     * @returns {number} Total number of elements
     */
    getTotalElementCount() {
        let total = 1; // The molecule itself
        
        for (const site of this.bindingSites) {
            total += 1; // The binding site
            total += site.stateCount; // States of the binding site
        }
        
        return total;
    }
    
    /**
     * Clone this molecule
     * @returns {Molecule} A new molecule with the same properties
     */
    clone() {
        const clone = new Molecule(
            this.name,
            this.count,
            this.diffusionTranslational,
            this.diffusionRotational,
            this.membraneAnchored,
            this.isLipid,
            this.isImplicitLipid
        );
        
        // Clone all binding sites
        for (const site of this.bindingSites) {
            clone.addBindingSite(site.clone());
        }
        
        return clone;
    }
    
    /**
     * Check if this molecule is identical to another
     * @param {Molecule} other - The other molecule to compare with
     * @returns {boolean} True if the molecules are identical
     */
    isIdenticalTo(other) {
        // Check basic properties
        if (this.name !== other.name || 
            this.count !== other.count ||
            this.membraneAnchored !== other.membraneAnchored ||
            this.diffusionTranslational !== other.diffusionTranslational ||
            this.diffusionRotational !== other.diffusionRotational ||
            this.isLipid !== other.isLipid ||
            this.isImplicitLipid !== other.isImplicitLipid ||
            this.bindingSites.length !== other.bindingSites.length) {
            return false;
        }
        
        // Check all binding sites
        for (let i = 0; i < this.bindingSites.length; i++) {
            const thisSite = this.bindingSites[i];
            const otherSite = other.bindingSites[i];
            
            if (thisSite.name !== otherSite.name ||
                thisSite.stateCount !== otherSite.stateCount ||
                !thisSite.position.equals(otherSite.position) ||
                thisSite.stateNames.length !== otherSite.stateNames.length) {
                return false;
            }
            
            // Check state names
            for (let j = 0; j < thisSite.stateNames.length; j++) {
                if (thisSite.stateNames[j] !== otherSite.stateNames[j]) {
                    return false;
                }
            }
        }
        
        return true;
    }
    
    /**
     * Convert the molecule to a string for BNGL format
     * Used for reaction definitions
     * @param {string} [bindingSite] - Optional binding site to include
     * @param {string} [state] - Optional state to include 
     * @param {string} [bond] - Optional bond to include
     * @returns {string} BNGL representation of the molecule
     */
    toBNGL(bindingSite, state, bond) {
        let result = this.name + '(';
        
        if (bindingSite) {
            result += bindingSite;
            
            if (state) {
                result += '~' + state;
            }
            
            if (bond) {
                result += '!' + bond;
            }
        }
        
        result += ')';
        return result;
    }
    
    /**
     * Convert the molecule to a plain object for serialization
     * @returns {Object} Plain object representation
     */
    toObject() {
        return {
            id: this.id,
            name: this.name,
            count: this.count,
            diffusionTranslational: this.diffusionTranslational,
            diffusionRotational: this.diffusionRotational,
            membraneAnchored: this.membraneAnchored,
            isLipid: this.isLipid,
            isImplicitLipid: this.isImplicitLipid,
            bindingSites: this.bindingSites.map(site => site.toObject())
        };
    }
    
    /**
     * Create a molecule from a plain object
     * @param {Object} obj - Object with molecule properties
     * @returns {Molecule} New molecule
     */
    static fromObject(obj) {
        const molecule = new Molecule(
            obj.name,
            obj.count,
            obj.diffusionTranslational,
            obj.diffusionRotational,
            obj.membraneAnchored,
            obj.isLipid,
            obj.isImplicitLipid
        );
        
        if (obj.id) {
            molecule.id = obj.id;
        }
        
        if (obj.bindingSites) {
            for (const siteObj of obj.bindingSites) {
                molecule.addBindingSite(BindingSite.fromObject(siteObj));
            }
        }
        
        return molecule;
    }
    
    /**
     * Create a molecule from MolecX.mol file format
     * @param {string} molFileContent - Content of a MolecX.mol file
     * @returns {Molecule} New molecule
     */
    static fromMolFile(molFileContent) {
        const lines = molFileContent.split(/\r?\n/);
        let name = "";
        let count = 0;
        let isLipid = false;
        let diffTranslational = 0;
        let diffRotational = 0;
        const bindingSites = [];
        
        // Parse header section
        for (let i = 0; i < lines.length; i++) {
            const line = lines[i].trim();
            
            if (line.startsWith("Name")) {
                name = line.split("=")[1].trim();
            } else if (line.startsWith("copies")) {
                count = parseInt(line.split("=")[1].trim(), 10);
            } else if (line.startsWith("isLipid")) {
                isLipid = line.split("=")[1].trim().toLowerCase() === "true";
            } else if (line.startsWith("D")) {
                // Parse translational diffusion
                const diffValues = line.split("=")[1].trim()
                    .replace(/[\[\]]/g, "")
                    .split(",")
                    .map(v => parseFloat(v.trim()));
                
                diffTranslational = diffValues[0]; // Use x component
            } else if (line.startsWith("Dr")) {
                // Parse rotational diffusion
                const diffValues = line.split("=")[1].trim()
                    .replace(/[\[\]]/g, "")
                    .split(",")
                    .map(v => parseFloat(v.trim()));
                
                diffRotational = diffValues[0]; // Use x component
            } else if (line.startsWith("#Coords")) {
                // Start of coordinates section
                for (let j = i + 1; j < lines.length; j++) {
                    const coordLine = lines[j].trim();
                    
                    if (coordLine === "" || coordLine.startsWith("#")) {
                        break;
                    }
                    
                    const parts = coordLine.split(/\s+/);
                    if (parts.length >= 4) {
                        const siteName = parts[0];
                        const x = parseFloat(parts[1]);
                        const y = parseFloat(parts[2]);
                        const z = parseFloat(parts[3]);
                        
                        // Skip center of mass
                        if (siteName !== "COM" && siteName !== "Cm") {
                            bindingSites.push({
                                name: siteName,
                                position: new Vector3(x, y, z),
                                stateCount: 0,
                                stateNames: ["-"]
                            });
                        }
                    }
                }
                
                break; // Stop parsing after coordinates
            }
        }
        
        // Create the molecule
        const molecule = new Molecule(
            name,
            count,
            diffTranslational,
            diffRotational,
            false, // membraneAnchored
            isLipid,
            false // isImplicitLipid
        );
        
        // Add binding sites
        for (const site of bindingSites) {
            molecule.addBindingSite(new BindingSite(
                site.name,
                site.position,
                site.stateCount,
                site.stateNames
            ));
        }
        
        return molecule;
    }
    
    /**
     * Export molecule to MolecX.mol file format
     * @returns {string} Content for a MolecX.mol file
     */
    toMolFile() {
        // Format numbers with 4 decimal places
        const format = (num) => num.toFixed(4);
        
        let content = "##\n";
        content += `# ${this.name} molecule information file\n`;
        content += "##\n\n";
        content += `Name\t= ${this.name}\n`;
        content += `copies\t= ${this.count}\n`;
        content += `isLipid = ${this.isLipid ? 'true' : 'false'}\n`;
        content += "checkOverlap = true\n";
        content += "mass = 1 #default value\n\n";
        
        content += "#translational diffusion\n";
        content += `D = [${format(this.diffusionTranslational)},${format(this.diffusionTranslational)},${format(0.0)}]\n\n`;
        
        content += "#rotational diffusion\n";
        content += `Dr = [${format(this.diffusionRotational)},${format(this.diffusionRotational)},${format(0.0)}]\n\n`;
        
        content += "#Coords\tx\ty\tz\n";
        content += "COM\t0.0000\t0.0000\t0.0000\n";
        
        for (const site of this.bindingSites) {
            content += `${site.name}\t${format(site.position.x)}\t${format(site.position.y)}\t${format(site.position.z)}\n`;
        }
        
        return content;
    }
}

// Export the classes
export { Molecule, BindingSite };