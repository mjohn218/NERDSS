/**
 * Quaternion.js
 * 
 * A quaternion class for representing 3D rotations.
 * Quaternions provide a more stable way to handle rotations compared to Euler angles,
 * avoiding issues like gimbal lock.
 */

import { Vector3 } from './Vector3.js';

/**
 * Quaternion class representing a rotation in 3D space
 */
class Quaternion {
    /**
     * Create a new quaternion
     * @param {number} x - X component
     * @param {number} y - Y component
     * @param {number} z - Z component
     * @param {number} w - W component (real part)
     */
    constructor(x = 0, y = 0, z = 0, w = 1) {
        this.x = x;
        this.y = y;
        this.z = z;
        this.w = w;
    }
    
    /**
     * Set quaternion components directly
     * @param {number} x - X component
     * @param {number} y - Y component
     * @param {number} z - Z component
     * @param {number} w - W component
     * @returns {Quaternion} This quaternion for chaining
     */
    set(x, y, z, w) {
        this.x = x;
        this.y = y;
        this.z = z;
        this.w = w;
        return this;
    }
    
    /**
     * Copy values from another quaternion
     * @param {Quaternion} q - Quaternion to copy from
     * @returns {Quaternion} This quaternion for chaining
     */
    copy(q) {
        this.x = q.x;
        this.y = q.y;
        this.z = q.z;
        this.w = q.w;
        return this;
    }
    
    /**
     * Create a clone of this quaternion
     * @returns {Quaternion} A new quaternion with the same values
     */
    clone() {
        return new Quaternion(this.x, this.y, this.z, this.w);
    }
    
    /**
     * Set quaternion from axis-angle representation
     * @param {Vector3} axis - Axis of rotation (will be normalized)
     * @param {number} angle - Angle of rotation in radians
     * @returns {Quaternion} This quaternion for chaining
     */
    setFromAxisAngle(axis, angle) {
        // Normalize the axis
        const normalizedAxis = axis.clone().normalize();
        
        // Calculate the quaternion components
        const halfAngle = angle / 2;
        const sinHalfAngle = Math.sin(halfAngle);
        
        this.x = normalizedAxis.x * sinHalfAngle;
        this.y = normalizedAxis.y * sinHalfAngle;
        this.z = normalizedAxis.z * sinHalfAngle;
        this.w = Math.cos(halfAngle);
        
        return this;
    }
    
    /**
     * Calculate the magnitude (length) of the quaternion
     * @returns {number} Quaternion magnitude
     */
    length() {
        return Math.sqrt(this.lengthSquared());
    }
    
    /**
     * Calculate the squared magnitude of the quaternion
     * @returns {number} Squared magnitude
     */
    lengthSquared() {
        return this.x * this.x + this.y * this.y + this.z * this.z + this.w * this.w;
    }
    
    /**
     * Normalize the quaternion (set length to 1)
     * @returns {Quaternion} This quaternion for chaining
     */
    normalize() {
        let length = this.length();
        
        if (length === 0) {
            this.x = 0;
            this.y = 0;
            this.z = 0;
            this.w = 1;
        } else {
            length = 1 / length;
            this.x *= length;
            this.y *= length;
            this.z *= length;
            this.w *= length;
        }
        
        return this;
    }
    
    /**
     * Calculate the conjugate of this quaternion
     * @returns {Quaternion} This quaternion for chaining
     */
    conjugate() {
        this.x *= -1;
        this.y *= -1;
        this.z *= -1;
        return this;
    }
    
    /**
     * Calculate the inverse of this quaternion
     * @returns {Quaternion} This quaternion for chaining
     */
    invert() {
        // For unit quaternions, the inverse is the conjugate
        return this.conjugate().normalize();
    }
    
    /**
     * Multiply this quaternion by another quaternion (composition of rotations)
     * @param {Quaternion} q - Quaternion to multiply with
     * @returns {Quaternion} This quaternion for chaining
     */
    multiply(q) {
        const x = this.x, y = this.y, z = this.z, w = this.w;
        
        this.x = w * q.x + x * q.w + y * q.z - z * q.y;
        this.y = w * q.y + y * q.w + z * q.x - x * q.z;
        this.z = w * q.z + z * q.w + x * q.y - y * q.x;
        this.w = w * q.w - x * q.x - y * q.y - z * q.z;
        
        return this;
    }
    
    /**
     * Premultiply this quaternion by another quaternion
     * @param {Quaternion} q - Quaternion to premultiply with
     * @returns {Quaternion} This quaternion for chaining
     */
    premultiply(q) {
        return this.clone().copy(q).multiply(this);
    }
    
    /**
     * Rotate a vector using this quaternion
     * @param {Vector3} v - Vector to rotate
     * @returns {Vector3} Rotated vector (modified in place)
     */
    multiplyVector3(v) {
        // Implementation of quaternion-vector multiplication: q*v*q^(-1)
        // This is optimized to avoid creating temporary quaternions
        
        const x = v.x, y = v.y, z = v.z;
        const qx = this.x, qy = this.y, qz = this.z, qw = this.w;
        
        // Calculate quaternion * vector
        const ix = qw * x + qy * z - qz * y;
        const iy = qw * y + qz * x - qx * z;
        const iz = qw * z + qx * y - qy * x;
        const iw = -qx * x - qy * y - qz * z;
        
        // Calculate result * quaternion^-1
        v.x = ix * qw + iw * -qx + iy * -qz - iz * -qy;
        v.y = iy * qw + iw * -qy + iz * -qx - ix * -qz;
        v.z = iz * qw + iw * -qz + ix * -qy - iy * -qx;
        
        return v;
    }
    
    /**
     * Calculate the angle (in radians) between this quaternion and another
     * @param {Quaternion} q - Other quaternion
     * @returns {number} Angle in radians
     */
    angleTo(q) {
        const dot = this.dot(q);
        return Math.acos(2 * dot * dot - 1);
    }
    
    /**
     * Calculate the dot product with another quaternion
     * @param {Quaternion} q - Other quaternion
     * @returns {number} Dot product
     */
    dot(q) {
        return this.x * q.x + this.y * q.y + this.z * q.z + this.w * q.w;
    }
    
    /**
     * Spherical linear interpolation between this quaternion and another
     * @param {Quaternion} q - Target quaternion
     * @param {number} t - Interpolation parameter (0-1)
     * @returns {Quaternion} This quaternion for chaining
     */
    slerp(q, t) {
        if (t === 0) return this;
        if (t === 1) return this.copy(q);
        
        const x = this.x, y = this.y, z = this.z, w = this.w;
        
        // Calculate cosine of angle between quaternions
        let cosHalfTheta = x * q.x + y * q.y + z * q.z + w * q.w;
        
        // Choose the shorter path
        if (cosHalfTheta < 0) {
            cosHalfTheta = -cosHalfTheta;
            q = new Quaternion(-q.x, -q.y, -q.z, -q.w);
        }
        
        // If quaternions are very close, use linear interpolation
        if (cosHalfTheta >= 0.999) {
            this.x = x + t * (q.x - x);
            this.y = y + t * (q.y - y);
            this.z = z + t * (q.z - z);
            this.w = w + t * (q.w - w);
            return this.normalize();
        }
        
        // Calculate interpolation parameters
        const halfTheta = Math.acos(cosHalfTheta);
        const sinHalfTheta = Math.sqrt(1 - cosHalfTheta * cosHalfTheta);
        
        if (Math.abs(sinHalfTheta) < 0.001) {
            this.x = 0.5 * (x + q.x);
            this.y = 0.5 * (y + q.y);
            this.z = 0.5 * (z + q.z);
            this.w = 0.5 * (w + q.w);
            return this.normalize();
        }
        
        // Calculate interpolation ratios
        const ratioA = Math.sin((1 - t) * halfTheta) / sinHalfTheta;
        const ratioB = Math.sin(t * halfTheta) / sinHalfTheta;
        
        // Calculate interpolated quaternion
        this.x = x * ratioA + q.x * ratioB;
        this.y = y * ratioA + q.y * ratioB;
        this.z = z * ratioA + q.z * ratioB;
        this.w = w * ratioA + q.w * ratioB;
        
        return this;
    }
    
    /**
     * Convert to a rotation matrix (3x3)
     * @returns {Array<Array<number>>} 3x3 rotation matrix
     */
    toRotationMatrix() {
        const x = this.x, y = this.y, z = this.z, w = this.w;
        const x2 = x + x, y2 = y + y, z2 = z + z;
        const xx = x * x2, xy = x * y2, xz = x * z2;
        const yy = y * y2, yz = y * z2, zz = z * z2;
        const wx = w * x2, wy = w * y2, wz = w * z2;
        
        return [
            [1 - (yy + zz), xy - wz, xz + wy],
            [xy + wz, 1 - (xx + zz), yz - wx],
            [xz - wy, yz + wx, 1 - (xx + yy)]
        ];
    }
    
    /**
     * Extract Euler angles from quaternion (in radians)
     * @returns {Object} Object with x, y, z properties (radians)
     */
    toEulerAngles() {
        // Convert to Euler angles (roll, pitch, yaw)
        const x = this.x, y = this.y, z = this.z, w = this.w;
        
        // Roll (x-axis rotation)
        const sinr_cosp = 2 * (w * x + y * z);
        const cosr_cosp = 1 - 2 * (x * x + y * y);
        const roll = Math.atan2(sinr_cosp, cosr_cosp);
        
        // Pitch (y-axis rotation)
        const sinp = 2 * (w * y - z * x);
        let pitch;
        if (Math.abs(sinp) >= 1) {
            // Use 90 degrees if out of range
            pitch = Math.sign(sinp) * Math.PI / 2;
        } else {
            pitch = Math.asin(sinp);
        }
        
        // Yaw (z-axis rotation)
        const siny_cosp = 2 * (w * z + x * y);
        const cosy_cosp = 1 - 2 * (y * y + z * z);
        const yaw = Math.atan2(siny_cosp, cosy_cosp);
        
        return { x: roll, y: pitch, z: yaw };
    }
    
    /**
     * Set quaternion from Euler angles (in radians)
     * @param {number} roll - Rotation around x-axis (radians)
     * @param {number} pitch - Rotation around y-axis (radians)
     * @param {number} yaw - Rotation around z-axis (radians)
     * @returns {Quaternion} This quaternion for chaining
     */
    setFromEulerAngles(roll, pitch, yaw) {
        // Calculate quaternion components using Euler angles
        const cr = Math.cos(roll * 0.5);
        const sr = Math.sin(roll * 0.5);
        const cp = Math.cos(pitch * 0.5);
        const sp = Math.sin(pitch * 0.5);
        const cy = Math.cos(yaw * 0.5);
        const sy = Math.sin(yaw * 0.5);
        
        this.w = cr * cp * cy + sr * sp * sy;
        this.x = sr * cp * cy - cr * sp * sy;
        this.y = cr * sp * cy + sr * cp * sy;
        this.z = cr * cp * sy - sr * sp * cy;
        
        return this;
    }
    
    /**
     * Get a string representation of the quaternion
     * @returns {string} String representation
     */
    toString() {
        return `Quaternion(${this.x}, ${this.y}, ${this.z}, ${this.w})`;
    }
}

export { Quaternion };