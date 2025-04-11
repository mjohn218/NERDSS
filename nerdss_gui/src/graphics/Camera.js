/**
 * Camera.js
 * 
 * This file defines the Camera class for 3D visualization in the NERDSS GUI.
 * It handles camera positioning, view transformations, and user interactions
 * like rotation, zooming, and panning.
 */

import { Vector3 } from '../math/Vector3.js';
import { Quaternion } from '../math/Quaternion.js';

/**
 * Camera class for 3D visualization
 */
class Camera {
    /**
     * Create a new camera
     */
    constructor() {
        // Camera position and orientation
        this.eye = new Vector3(0, 0, 30);     // Eye position
        this.target = new Vector3(0, 0, 0);   // Look-at target
        this.up = new Vector3(0, 1, 0);       // Up vector
        
        // View volume parameters
        this.xMin = -10;
        this.xMax = 10;
        this.yMin = -10;
        this.yMax = 10;
        this.zMin = -10;
        this.zMax = 10;
        
        this.orthographic = false;      // Orthographic or perspective projection
        this.preserveAspect = true;     // Preserve aspect ratio
        
        // Variables for trackball interaction
        this.trackball = {
            dragging: false,
            prevPos: null,
            element: null
        };
        
        // Bind methods that are used as event handlers
        this.onMouseDown = this.onMouseDown.bind(this);
        this.onMouseMove = this.onMouseMove.bind(this);
        this.onMouseUp = this.onMouseUp.bind(this);
        this.onWheel = this.onWheel.bind(this);
    }
    
    /**
     * Set the view limits
     * @param {number} xMin - Minimum X value
     * @param {number} xMax - Maximum X value
     * @param {number} yMin - Minimum Y value
     * @param {number} yMax - Maximum Y value
     * @param {number} zMin - Minimum Z value
     * @param {number} zMax - Maximum Z value
     */
    setLimits(xMin, xMax, yMin, yMax, zMin, zMax) {
        this.xMin = xMin;
        this.xMax = xMax;
        this.yMin = yMin;
        this.yMax = yMax;
        this.zMin = zMin;
        this.zMax = zMax;
    }
    
    /**
     * Set the camera scale (convenience method)
     * @param {number} scale - Scale factor
     */
    setScale(scale) {
        const absScale = Math.abs(scale);
        this.setLimits(-absScale, absScale, -absScale, absScale, -2 * absScale, 2 * absScale);
    }
    
    /**
     * Set the camera to look at a specific point
     * @param {number} eyeX - Eye X position
     * @param {number} eyeY - Eye Y position
     * @param {number} eyeZ - Eye Z position
     * @param {number} targetX - Target X position
     * @param {number} targetY - Target Y position
     * @param {number} targetZ - Target Z position
     * @param {number} upX - Up vector X component
     * @param {number} upY - Up vector Y component
     * @param {number} upZ - Up vector Z component
     */
    lookAt(eyeX, eyeY, eyeZ, targetX, targetY, targetZ, upX, upY, upZ) {
        this.eye.set(eyeX, eyeY, eyeZ);
        this.target.set(targetX, targetY, targetZ);
        this.up.set(upX, upY, upZ);
    }
    
    /**
     * Get the view parameters
     * @returns {Object} Object with eye, target, and up properties
     */
    getViewParameters() {
        return {
            eye: this.eye.clone(),
            target: this.target.clone(),
            up: this.up.clone()
        };
    }
    
    /**
     * Get the view limits
     * @returns {number[]} Array of [xMin, xMax, yMin, yMax, zMin, zMax]
     */
    getLimits() {
        return [this.xMin, this.xMax, this.yMin, this.yMax, this.zMin, this.zMax];
    }
    
    /**
     * Apply the camera to a WebGL context
     * @param {WebGLRenderingContext} gl - WebGL rendering context
     */
    apply(gl) {
        // Get viewport dimensions
        const viewport = gl.getViewport ? gl.getViewport() : [0, 0, gl.canvas.width, gl.canvas.height];
        const viewWidth = viewport[2];
        const viewHeight = viewport[3];
        
        // Calculate actual view volume dimensions (respecting aspect ratio if required)
        let xMinActual = this.xMin;
        let xMaxActual = this.xMax;
        let yMinActual = this.yMin;
        let yMaxActual = this.yMax;
        
        if (this.preserveAspect) {
            const viewAspect = viewHeight / viewWidth;
            const windowWidth = xMaxActual - xMinActual;
            const windowHeight = yMaxActual - yMinActual;
            const windowAspect = windowHeight / windowWidth;
            
            if (windowAspect > viewAspect) {
                // Expand width
                const extra = (windowAspect / viewAspect - 1.0) * windowWidth / 2.0;
                xMinActual -= extra;
                xMaxActual += extra;
            } else if (viewAspect > windowAspect) {
                // Expand height
                const extra = (viewAspect / windowAspect - 1.0) * windowHeight / 2.0;
                yMinActual -= extra;
                yMaxActual += extra;
            }
        }
        
        // Calculate view distance (distance from eye to target)
        const viewVector = new Vector3().copy(this.target).subtract(this.eye);
        const viewDistance = viewVector.length();
        
        // Set up projection matrix
        gl.matrixMode(gl.PROJECTION);
        gl.loadIdentity();
        
        if (this.orthographic) {
            // Orthographic projection
            gl.ortho(
                xMinActual, xMaxActual,
                yMinActual, yMaxActual,
                viewDistance - this.zMax,
                viewDistance - this.zMin
            );
        } else {
            // Perspective projection
            let near = viewDistance - this.zMax;
            if (near < 0.1) near = 0.1;
            
            const centerX = (xMinActual + xMaxActual) / 2;
            const centerY = (yMinActual + yMaxActual) / 2;
            const newWidth = (near / viewDistance) * (xMaxActual - xMinActual);
            const newHeight = (near / viewDistance) * (yMaxActual - yMinActual);
            const x1 = centerX - newWidth / 2;
            const x2 = centerX + newWidth / 2;
            const y1 = centerY - newHeight / 2;
            const y2 = centerY + newHeight / 2;
            
            gl.frustum(
                x1, x2, y1, y2,
                near,
                viewDistance - this.zMin
            );
        }
        
        // Set up modelview matrix
        gl.matrixMode(gl.MODELVIEW);
        gl.loadIdentity();
        
        // Apply lookAt transformation
        this.gluLookAt(gl, 
            this.eye.x, this.eye.y, this.eye.z,
            this.target.x, this.target.y, this.target.z,
            this.up.x, this.up.y, this.up.z
        );
    }
    
    /**
     * Implementation of gluLookAt for WebGL
     * @param {WebGLRenderingContext} gl - WebGL rendering context
     * @param {number} eyeX - Eye X position
     * @param {number} eyeY - Eye Y position
     * @param {number} eyeZ - Eye Z position
     * @param {number} centerX - Target X position
     * @param {number} centerY - Target Y position
     * @param {number} centerZ - Target Z position
     * @param {number} upX - Up vector X component
     * @param {number} upY - Up vector Y component
     * @param {number} upZ - Up vector Z component
     */
    gluLookAt(gl, eyeX, eyeY, eyeZ, centerX, centerY, centerZ, upX, upY, upZ) {
        // Calculate look vector (normalized)
        const forward = new Vector3(
            centerX - eyeX,
            centerY - eyeY,
            centerZ - eyeZ
        ).normalize();
        
        // Calculate side vector (normalized)
        const side = new Vector3().copy(forward).cross(new Vector3(upX, upY, upZ)).normalize();
        
        // Recalculate up vector
        const up = new Vector3().copy(side).cross(forward);
        
        // Create rotation matrix
        const rotMat = [
            side.x, up.x, -forward.x, 0,
            side.y, up.y, -forward.y, 0,
            side.z, up.z, -forward.z, 0,
            0, 0, 0, 1
        ];
        
        // Apply rotation
        gl.multMatrixf(rotMat);
        
        // Apply translation
        gl.translatef(-eyeX, -eyeY, -eyeZ);
    }
    
    /**
     * Install a trackball for user interaction on a canvas element
     * @param {HTMLCanvasElement} canvas - Canvas element to attach trackball to
     */
    installTrackball(canvas) {
        // Remove any existing trackball
        if (this.trackball.element && this.trackball.element !== canvas) {
            this.trackball.element.removeEventListener('mousedown', this.onMouseDown);
            this.trackball.element.removeEventListener('wheel', this.onWheel);
        }
        
        // Set up new trackball
        this.trackball.element = canvas;
        
        if (canvas) {
            canvas.addEventListener('mousedown', this.onMouseDown);
            canvas.addEventListener('wheel', this.onWheel);
        }
    }
    
    /**
     * Mouse down event handler for trackball
     * @param {MouseEvent} event - Mouse event
     */
    onMouseDown(event) {
        if (this.trackball.dragging) return;
        
        this.trackball.dragging = true;
        this.trackball.prevPos = this.mousePointToRay(event.clientX, event.clientY);
        
        document.addEventListener('mousemove', this.onMouseMove);
        document.addEventListener('mouseup', this.onMouseUp);
    }
    
    /**
     * Mouse move event handler for trackball
     * @param {MouseEvent} event - Mouse event
     */
    onMouseMove(event) {
        if (!this.trackball.dragging) return;
        
        const currentPos = this.mousePointToRay(event.clientX, event.clientY);
        this.applyTransvection(this.trackball.prevPos, currentPos);
        this.trackball.prevPos = currentPos;
        
        // Request redraw
        if (this.trackball.element) {
            if (typeof this.trackball.element.requestRedraw === 'function') {
                this.trackball.element.requestRedraw();
            } else {
                // Trigger a redraw through canvas update event
                const event = new CustomEvent('camera-update');
                this.trackball.element.dispatchEvent(event);
            }
        }
    }
    
    /**
     * Mouse up event handler for trackball
     * @param {MouseEvent} event - Mouse event
     */
    onMouseUp(event) {
        if (!this.trackball.dragging) return;
        
        this.trackball.dragging = false;
        document.removeEventListener('mousemove', this.onMouseMove);
        document.removeEventListener('mouseup', this.onMouseUp);
    }
    
    /**
     * Mouse wheel event handler for zooming
     * @param {WheelEvent} event - Wheel event
     */
    onWheel(event) {
        event.preventDefault();
        
        // Adjust scale based on wheel direction
        const delta = event.deltaY > 0 ? 1.1 : 0.9;
        const limits = this.getLimits();
        
        this.setScale(limits[1] * delta);
        
        // Request redraw
        if (this.trackball.element) {
            if (typeof this.trackball.element.requestRedraw === 'function') {
                this.trackball.element.requestRedraw();
            } else {
                // Trigger a redraw through canvas update event
                const event = new CustomEvent('camera-update');
                this.trackball.element.dispatchEvent(event);
            }
        }
    }
    
    /**
     * Convert mouse point to 3D ray
     * @param {number} clientX - Mouse X coordinate
     * @param {number} clientY - Mouse Y coordinate
     * @returns {Vector3} Ray direction vector
     * @private
     */
    mousePointToRay(clientX, clientY) {
        const canvas = this.trackball.element;
        if (!canvas) return new Vector3(0, 0, 1);
        
        const rect = canvas.getBoundingClientRect();
        const x = clientX - rect.left;
        const y = clientY - rect.top;
        
        const centerX = canvas.width / 2;
        const centerY = canvas.height / 2;
        const scale = 0.8 * Math.min(centerX, centerY);
        
        const dx = x - centerX;
        const dy = centerY - y;  // Flip Y since screen coords are top-down
        
        // Calculate ray that passes through mouse point
        let dz;
        const norm = Math.sqrt(dx*dx + dy*dy);
        
        if (norm >= scale) {
            // Point outside the trackball sphere
            dz = 0;
        } else {
            // Point on the trackball sphere - use Pythagoras to get Z
            dz = Math.sqrt(scale*scale - dx*dx - dy*dy);
        }
        
        // Normalize the ray
        const length = Math.sqrt(dx*dx + dy*dy + dz*dz);
        return new Vector3(dx/length, dy/length, dz/length);
    }
    
    /**
     * Apply a transvection (reflection) transformation
     * Used for trackball rotation
     * @param {Vector3} v1 - First ray
     * @param {Vector3} v2 - Second ray
     * @private
     */
    applyTransvection(v1, v2) {
        // Get current camera directions
        const forwardDir = new Vector3().copy(this.target).subtract(this.eye).normalize();
        const rightDir = new Vector3().copy(this.up).cross(forwardDir).normalize();
        const upDir = new Vector3().copy(forwardDir).cross(rightDir).normalize();
        
        // Transform the rays to camera space
        const v1Camera = this.transformToViewCoords(v1, rightDir, upDir, forwardDir);
        const v2Camera = this.transformToViewCoords(v2, rightDir, upDir, forwardDir);
        
        // Calculate reflection axis
        const axis = new Vector3().copy(v1Camera).add(v2Camera).normalize();
        
        // Calculate rotation amount
        const angle = v1Camera.angle(v2Camera);
        
        // Apply rotation to current camera vectors
        const rotation = new Quaternion().setFromAxisAngle(axis, 2 * angle);
        
        // Rotate the view direction
        rotation.multiplyVector3(forwardDir);
        
        // Rotate the up vector
        rotation.multiplyVector3(upDir);
        
        // Recalculate eye position based on the rotated view direction
        const viewDistance = this.eye.clone().subtract(this.target).length();
        this.eye = this.target.clone().add(forwardDir.clone().scale(-viewDistance));
        
        // Update up vector
        this.up = upDir;
    }
    
    /**
     * Transform a vector to view coordinates
     * @param {Vector3} v - Vector to transform
     * @param {Vector3} rightDir - Right direction vector in world space
     * @param {Vector3} upDir - Up direction vector in world space
     * @param {Vector3} forwardDir - Forward direction vector in world space
     * @returns {Vector3} Transformed vector
     * @private
     */
    transformToViewCoords(v, rightDir, upDir, forwardDir) {
        const result = new Vector3();
        
        result.x = v.x * rightDir.x + v.y * upDir.x + v.z * forwardDir.x;
        result.y = v.x * rightDir.y + v.y * upDir.y + v.z * forwardDir.y;
        result.z = v.x * rightDir.z + v.y * upDir.z + v.z * forwardDir.z;
        
        return result;
    }
    
    /**
     * Reset the camera to default position
     */
    reset() {
        this.eye = new Vector3(0, 0, 30);
        this.target = new Vector3(0, 0, 0);
        this.up = new Vector3(0, 1, 0);
        this.setScale(10);
    }
}

export { Camera };