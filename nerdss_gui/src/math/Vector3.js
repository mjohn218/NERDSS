/**
 * Vector3.js
 * 
 * A 3D vector class for mathematical operations in 3D space.
 * Used for molecular coordinates, normals, and transformations.
 */

/**
 * Vector3 class representing a 3D vector
 */
class Vector3 {
    /**
     * Create a new 3D vector
     * @param {number} x - X component
     * @param {number} y - Y component
     * @param {number} z - Z component
     */
    constructor(x = 0, y = 0, z = 0) {
        this.x = x;
        this.y = y;
        this.z = z;
    }
    
    /**
     * Set the vector components
     * @param {number} x - X component
     * @param {number} y - Y component
     * @param {number} z - Z component
     * @returns {Vector3} This vector for chaining
     */
    set(x, y, z) {
        this.x = x;
        this.y = y;
        this.z = z;
        return this;
    }
    
    /**
     * Copy values from another vector
     * @param {Vector3} v - Vector to copy from
     * @returns {Vector3} This vector for chaining
     */
    copy(v) {
        this.x = v.x;
        this.y = v.y;
        this.z = v.z;
        return this;
    }
    
    /**
     * Create a clone of this vector
     * @returns {Vector3} A new vector with the same values
     */
    clone() {
        return new Vector3(this.x, this.y, this.z);
    }
    
    /**
     * Add another vector to this one
     * @param {Vector3} v - Vector to add
     * @returns {Vector3} This vector for chaining
     */
    add(v) {
        this.x += v.x;
        this.y += v.y;
        this.z += v.z;
        return this;
    }
    
    /**
     * Subtract another vector from this one
     * @param {Vector3} v - Vector to subtract
     * @returns {Vector3} This vector for chaining
     */
    subtract(v) {
        this.x -= v.x;
        this.y -= v.y;
        this.z -= v.z;
        return this;
    }
    
    /**
     * Multiply this vector by a scalar
     * @param {number} scalar - Value to multiply by
     * @returns {Vector3} This vector for chaining
     */
    scale(scalar) {
        this.x *= scalar;
        this.y *= scalar;
        this.z *= scalar;
        return this;
    }
    
    /**
     * Calculate the dot product with another vector
     * @param {Vector3} v - Other vector
     * @returns {number} Dot product
     */
    dot(v) {
        return this.x * v.x + this.y * v.y + this.z * v.z;
    }
    
    /**
     * Calculate the cross product with another vector
     * @param {Vector3} v - Other vector
     * @returns {Vector3} This vector for chaining
     */
    cross(v) {
        const x = this.y * v.z - this.z * v.y;
        const y = this.z * v.x - this.x * v.z;
        const z = this.x * v.y - this.y * v.x;
        
        this.x = x;
        this.y = y;
        this.z = z;
        
        return this;
    }
    
    /**
     * Calculate the length (magnitude) of the vector
     * @returns {number} Length
     */
    length() {
        return Math.sqrt(this.lengthSquared());
    }
    
    /**
     * Calculate the squared length (magnitude) of the vector
     * More efficient than length() when only comparing lengths
     * @returns {number} Squared length
     */
    lengthSquared() {
        return this.x * this.x + this.y * this.y + this.z * this.z;
    }
    
    /**
     * Normalize the vector (set length to 1)
     * @returns {Vector3} This vector for chaining
     */
    normalize() {
        const length = this.length();
        
        if (length > 0) {
            this.scale(1 / length);
        }
        
        return this;
    }
    
    /**
     * Calculate the angle (in radians) between this vector and another
     * @param {Vector3} v - Other vector
     * @returns {number} Angle in radians
     */
    angle(v) {
        const denominator = Math.sqrt(this.lengthSquared() * v.lengthSquared());
        
        if (denominator === 0) {
            return 0;
        }
        
        // Ensure the dot product is in range [-1, 1]
        const dot = Math.max(-1, Math.min(1, this.dot(v) / denominator));
        return Math.acos(dot);
    }
    
    /**
     * Project this vector onto a plane defined by a normal vector
     * @param {Vector3} normal - Normal vector of the plane
     * @returns {Vector3} This vector for chaining
     */
    projectOnPlane(normal) {
        const normalClone = normal.clone().normalize();
        const dot = this.dot(normalClone);
        
        // Calculate projection along normal and subtract from this vector
        this.x -= normalClone.x * dot;
        this.y -= normalClone.y * dot;
        this.z -= normalClone.z * dot;
        
        return this;
    }
    
    /**
     * Rotate the vector around an axis by an angle
     * @param {Vector3} axis - Axis of rotation (will be normalized)
     * @param {number} angle - Angle in radians
     * @returns {Vector3} This vector for chaining
     */
    rotateAroundAxis(axis, angle) {
        // Implementation of Rodrigues' rotation formula
        const normalizedAxis = axis.clone().normalize();
        const cosTheta = Math.cos(angle);
        const sinTheta = Math.sin(angle);
        
        // v_rot = v * cos(θ) + (axis × v) * sin(θ) + axis * (axis · v) * (1 - cos(θ))
        const dot = this.dot(normalizedAxis);
        const cross = normalizedAxis.clone().cross(this);
        
        this.x = this.x * cosTheta + cross.x * sinTheta + normalizedAxis.x * dot * (1 - cosTheta);
        this.y = this.y * cosTheta + cross.y * sinTheta + normalizedAxis.y * dot * (1 - cosTheta);
        this.z = this.z * cosTheta + cross.z * sinTheta + normalizedAxis.z * dot * (1 - cosTheta);
        
        return this;
    }
    
    /**
     * Check if this vector is equal to another
     * @param {Vector3} v - Other vector
     * @param {number} [epsilon=1e-6] - Tolerance for floating-point comparison
     * @returns {boolean} True if vectors are equal within epsilon
     */
    equals(v, epsilon = 1e-6) {
        return (
            Math.abs(this.x - v.x) < epsilon &&
            Math.abs(this.y - v.y) < epsilon &&
            Math.abs(this.z - v.z) < epsilon
        );
    }
    
    /**
     * Convert vector to array [x, y, z]
     * @returns {number[]} Array with x, y, z components
     */
    toArray() {
        return [this.x, this.y, this.z];
    }
    
    /**
     * Get a string representation of the vector
     * @returns {string} String representation
     */
    toString() {
        return `Vector3(${this.x}, ${this.y}, ${this.z})`;
    }
    
    /**
     * Create a Vector3 from an array
     * @param {number[]} array - Array with x, y, z components
     * @returns {Vector3} New Vector3
     */
    static fromArray(array) {
        return new Vector3(array[0], array[1], array[2]);
    }
}

export { Vector3 };