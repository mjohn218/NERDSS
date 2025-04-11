/**
 * Renderer.js
 * 
 * This file defines the Renderer class for 3D visualization of molecules
 * and reactions in the NERDSS GUI. It handles WebGL rendering of molecular
 * structures, binding sites, and reaction geometries.
 */

import { Vector3 } from '../math/Vector3.js';
import { Quaternion } from '../math/Quaternion.js';
import { ShapeGenerator } from './ShapeGenerator.js';

/**
 * Renderer class for 3D visualization
 */
class Renderer {
    /**
     * Create a new renderer
     */
    constructor() {
        // WebGL contexts
        this.moleculeContext = null;
        this.reactionContext = null;
        
        // Cameras
        this.moleculeCamera = null;
        this.reactionCamera = null;
        
        // Canvas elements
        this.moleculeCanvas = null;
        this.reactionCanvas = null;
        
        // Data references
        this.moleculeManager = null;
        this.reactionManager = null;
        
        // Display lists and buffers
        this.sphereDisplayList = null;
        this.cylinderDisplayList = null;
        this.sphereVertexBuffer = null;
        this.sphereNormalBuffer = null;
        this.cylinderVertexBuffer = null;
        this.cylinderNormalBuffer = null;
        
        // WebGL buffer IDs
        this.vertexVboId = null;
        this.normalVboId = null;
        
        // Shapes for rendering
        this.shapeGenerator = null;
        
        // Animation state
        this.animationFrameId = null;
        
        // Bind methods used as callbacks
        this.renderMoleculeView = this.renderMoleculeView.bind(this);
        this.renderReactionView = this.renderReactionView.bind(this);
    }
    
    /**
     * Set data managers
     * @param {MoleculeManager} moleculeManager - Molecule manager
     * @param {ReactionManager} reactionManager - Reaction manager
     */
    setManagers(moleculeManager, reactionManager) {
        this.moleculeManager = moleculeManager;
        this.reactionManager = reactionManager;
        
        // Set up event listeners for data changes
        if (moleculeManager) {
            moleculeManager.on('molecules-changed', this.refreshMoleculeView.bind(this));
            moleculeManager.on('selection-changed', this.refreshMoleculeView.bind(this));
        }
        
        if (reactionManager) {
            reactionManager.on('reactions-changed', this.refreshReactionView.bind(this));
            reactionManager.on('selection-changed', this.refreshReactionView.bind(this));
        }
    }
    
    /**
     * Initialize the molecule renderer
     * @param {HTMLCanvasElement} canvas - Canvas element for rendering
     * @param {Camera} camera - Camera for the view
     */
    initMoleculeRenderer(canvas, camera) {
        this.moleculeCanvas = canvas;
        this.moleculeCamera = camera;
        
        // Initialize WebGL context
        try {
            this.moleculeContext = canvas.getContext('webgl') || canvas.getContext('experimental-webgl');
            if (!this.moleculeContext) {
                throw new Error('WebGL not supported');
            }
        } catch (error) {
            console.error('Error initializing WebGL:', error);
            this.showWebGLError(canvas);
            return;
        }
        
        // Initialize WebGL state
        this.initWebGL(this.moleculeContext);
        
        // Set up shape generator
        this.shapeGenerator = new ShapeGenerator(this.moleculeContext);
        
        // Create display lists and buffers for shapes
        this.createShapeDisplayLists(this.moleculeContext);
        
        // Set up animation loop
        canvas.addEventListener('camera-update', this.refreshMoleculeView.bind(this));
        this.startMoleculeRenderLoop();
    }
    
    /**
     * Initialize the reaction renderer
     * @param {HTMLCanvasElement} canvas - Canvas element for rendering
     * @param {Camera} camera - Camera for the view
     */
    initReactionRenderer(canvas, camera) {
        this.reactionCanvas = canvas;
        this.reactionCamera = camera;
        
        // Initialize WebGL context
        try {
            this.reactionContext = canvas.getContext('webgl') || canvas.getContext('experimental-webgl');
            if (!this.reactionContext) {
                throw new Error('WebGL not supported');
            }
        } catch (error) {
            console.error('Error initializing WebGL:', error);
            this.showWebGLError(canvas);
            return;
        }
        
        // Initialize WebGL state
        this.initWebGL(this.reactionContext);
        
        // Use the same shape generator for both views
        if (!this.shapeGenerator) {
            this.shapeGenerator = new ShapeGenerator(this.reactionContext);
        }
        
        // Create display lists and buffers for shapes
        this.createShapeDisplayLists(this.reactionContext);
        
        // Set up animation loop
        canvas.addEventListener('camera-update', this.refreshReactionView.bind(this));
        this.startReactionRenderLoop();
    }
    
    /**
     * Initialize WebGL state
     * @param {WebGLRenderingContext} gl - WebGL context
     * @private
     */
    initWebGL(gl) {
        // Enable depth testing
        gl.enable(gl.DEPTH_TEST);
        
        // Enable lighting
        gl.enable(gl.LIGHTING);
        gl.enable(gl.LIGHT0);
        
        // Set up material properties
        gl.enable(gl.COLOR_MATERIAL);
        gl.colorMaterial(gl.FRONT_AND_BACK, gl.AMBIENT_AND_DIFFUSE);
        
        // Set up light position
        const lightPos = [10.0, 10.0, 10.0, 0.0]; // Directional light
        gl.lightfv(gl.LIGHT0, gl.POSITION, lightPos);
        
        // Set up ambient and diffuse light colors
        const ambientLight = [0.3, 0.3, 0.3, 1.0];
        const diffuseLight = [0.7, 0.7, 0.7, 1.0];
        gl.lightfv(gl.LIGHT0, gl.AMBIENT, ambientLight);
        gl.lightfv(gl.LIGHT0, gl.DIFFUSE, diffuseLight);
        
        // Set up default material properties
        const defaultAmbient = [0.2, 0.2, 0.2, 1.0];
        const defaultDiffuse = [0.8, 0.8, 0.8, 1.0];
        const defaultSpecular = [0.0, 0.0, 0.0, 1.0];
        gl.materialfv(gl.FRONT_AND_BACK, gl.AMBIENT, defaultAmbient);
        gl.materialfv(gl.FRONT_AND_BACK, gl.DIFFUSE, defaultDiffuse);
        gl.materialfv(gl.FRONT_AND_BACK, gl.SPECULAR, defaultSpecular);
        gl.materialf(gl.FRONT_AND_BACK, gl.SHININESS, 0.0);
        
        // Set up smooth shading
        gl.shadeModel(gl.SMOOTH);
        
        // Set clear color to white
        gl.clearColor(1.0, 1.0, 1.0, 1.0);
    }
    
    /**
     * Create display lists and buffers for shapes
     * @param {WebGLRenderingContext} gl - WebGL context
     * @private
     */
    createShapeDisplayLists(gl) {
        // Create display list for sphere
        this.sphereDisplayList = gl.genLists(1);
        gl.newList(this.sphereDisplayList, gl.COMPILE);
        this.shapeGenerator.drawSphere(0.4, 32, 16);
        gl.endList();
        
        // Create display list for cylinder
        this.cylinderDisplayList = gl.genLists(1);
        gl.newList(this.cylinderDisplayList, gl.COMPILE);
        this.shapeGenerator.drawCylinder(0.2, 1.0, 16, 1);
        gl.endList();
        
        // Create vertex and normal buffers for sphere
        const sphereData = this.shapeGenerator.generateSphereData(0.4, 32, 16);
        this.sphereVertexBuffer = gl.createBuffer();
        gl.bindBuffer(gl.ARRAY_BUFFER, this.sphereVertexBuffer);
        gl.bufferData(gl.ARRAY_BUFFER, new Float32Array(sphereData.vertices), gl.STATIC_DRAW);
        
        this.sphereNormalBuffer = gl.createBuffer();
        gl.bindBuffer(gl.ARRAY_BUFFER, this.sphereNormalBuffer);
        gl.bufferData(gl.ARRAY_BUFFER, new Float32Array(sphereData.normals), gl.STATIC_DRAW);
        
        // Create vertex and normal buffers for cylinder
        const cylinderData = this.shapeGenerator.generateCylinderData(0.2, 1.0, 16, 1);
        this.cylinderVertexBuffer = gl.createBuffer();
        gl.bindBuffer(gl.ARRAY_BUFFER, this.cylinderVertexBuffer);
        gl.bufferData(gl.ARRAY_BUFFER, new Float32Array(cylinderData.vertices), gl.STATIC_DRAW);
        
        this.cylinderNormalBuffer = gl.createBuffer();
        gl.bindBuffer(gl.ARRAY_BUFFER, this.cylinderNormalBuffer);
        gl.bufferData(gl.ARRAY_BUFFER, new Float32Array(cylinderData.normals), gl.STATIC_DRAW);
        
        // Reset buffer binding
        gl.bindBuffer(gl.ARRAY_BUFFER, null);
    }
    
    /**
     * Start the render loop for the molecule view
     * @private
     */
    startMoleculeRenderLoop() {
        const renderLoop = () => {
            this.renderMoleculeView();
            this.moleculeAnimationFrameId = requestAnimationFrame(renderLoop);
        };
        
        renderLoop();
    }
    
    /**
     * Start the render loop for the reaction view
     * @private
     */
    startReactionRenderLoop() {
        const renderLoop = () => {
            this.renderReactionView();
            this.reactionAnimationFrameId = requestAnimationFrame(renderLoop);
        };
        
        renderLoop();
    }
    
    /**
     * Stop all render loops
     */
    stopRenderLoops() {
        if (this.moleculeAnimationFrameId) {
            cancelAnimationFrame(this.moleculeAnimationFrameId);
            this.moleculeAnimationFrameId = null;
        }
        
        if (this.reactionAnimationFrameId) {
            cancelAnimationFrame(this.reactionAnimationFrameId);
            this.reactionAnimationFrameId = null;
        }
    }
    
    /**
     * Render the molecule view
     */
    renderMoleculeView() {
        if (!this.moleculeContext || !this.moleculeCamera) return;
        
        const gl = this.moleculeContext;
        
        // Clear the canvas
        gl.clear(gl.COLOR_BUFFER_BIT | gl.DEPTH_BUFFER_BIT);
        
        // Apply camera
        this.moleculeCamera.apply(gl);
        
        // Rotate for better view
        gl.rotatef(-45, 0.0, 1.0, 0.0);
        
        // Draw coordinate axes
        this.drawAxes(gl);
        
        // Render molecules if available
        if (this.moleculeManager) {
            this.renderMolecules(gl);
        }
        
        // Ensure all commands are sent to the GPU
        gl.flush();
    }
    
    /**
     * Render the reaction view
     */
    renderReactionView() {
        if (!this.reactionContext || !this.reactionCamera) return;
        
        const gl = this.reactionContext;
        
        // Clear the canvas
        gl.clear(gl.COLOR_BUFFER_BIT | gl.DEPTH_BUFFER_BIT);
        
        // Apply camera
        this.reactionCamera.apply(gl);
        
        // Render selected reaction if available
        if (this.reactionManager && this.moleculeManager) {
            this.renderReaction(gl);
        }
        
        // Ensure all commands are sent to the GPU
        gl.flush();
    }
    
    /**
     * Render all molecules in the molecule view
     * @param {WebGLRenderingContext} gl - WebGL context
     * @private
     */
    renderMolecules(gl) {
        // Get the selected molecule
        const selection = this.moleculeManager.getSelection();
        const selectedMoleculeId = selection.moleculeId;
        const selectedBindingSiteIndex = selection.bindingSiteIndex;
        
        // Get all molecules
        const molecules = this.moleculeManager.getAllMolecules();
        if (!molecules || molecules.length === 0) return;
        
        // Only render the selected molecule if one is selected
        const moleculesToRender = selectedMoleculeId ? 
            [this.moleculeManager.getMolecule(selectedMoleculeId)] : 
            molecules;
        
        for (const molecule of moleculesToRender) {
            // Draw center of mass
            gl.pushMatrix();
            this.setColor(gl, 1.0, 1.0, 1.0); // White
            gl.callList(this.sphereDisplayList);
            gl.popMatrix();
            
            // Draw binding sites
            for (let i = 0; i < molecule.bindingSites.length; i++) {
                const site = molecule.bindingSites[i];
                const isSelected = (selectedMoleculeId === molecule.id && selectedBindingSiteIndex === i);
                
                gl.pushMatrix();
                
                // Translate to binding site position
                gl.translatef(site.position.x, site.position.y, site.position.z);
                
                // Draw binding site
                if (isSelected) {
                    this.setColor(gl, 1.0, 0.0, 0.0); // Red for selected
                } else {
                    this.setColor(gl, 1.0, 1.0, 0.0); // Yellow for normal
                }
                gl.callList(this.sphereDisplayList);
                
                gl.popMatrix();
                
                // Draw connecting cylinder from center to binding site
                gl.pushMatrix();
                
                // Calculate rotation to align cylinder with binding site direction
                const direction = site.position.clone().normalize();
                const length = site.position.length();
                
                // Default cylinder is along z-axis, need to rotate it to align with direction
                const zAxis = new Vector3(0, 0, 1);
                const rotationAxis = new Vector3().copy(zAxis).cross(direction);
                const angle = Math.acos(zAxis.dot(direction)) * 180 / Math.PI;
                
                // Apply rotation and scale to align with binding site
                if (rotationAxis.length() > 0.001) {
                    gl.rotatef(angle, rotationAxis.x, rotationAxis.y, rotationAxis.z);
                }
                
                // Draw the cylinder
                this.setColor(gl, 1.0, 1.0, 0.0); // Yellow for connections
                gl.scalef(1, 1, length);
                gl.callList(this.cylinderDisplayList);
                
                gl.popMatrix();
            }
        }
    }
    
    /**
     * Render the selected reaction in the reaction view
     * @param {WebGLRenderingContext} gl - WebGL context
     * @private
     */
    renderReaction(gl) {
        // Get the selected reaction
        const selection = this.reactionManager.getSelection();
        if (!selection.reactionId) return;
        
        const reaction = selection.reaction;
        if (!reaction) return;
        
        // Get reactants (for bimolecular reactions)
        if (reaction.reactants.length !== 2) {
            // Currently only supporting visualization of bimolecular reactions
            this.drawTextMessage(gl, "Only bimolecular reactions can be visualized");
            return;
        }
        
        // Find the molecules
        const molecule1Id = reaction.reactants[0].moleculeId;
        const molecule2Id = reaction.reactants[1].moleculeId;
        
        const molecule1 = this.moleculeManager.getMolecule(molecule1Id);
        const molecule2 = this.moleculeManager.getMolecule(molecule2Id);
        
        if (!molecule1 || !molecule2) {
            this.drawTextMessage(gl, "Reactant molecules not found");
            return;
        }
        
        // Find the binding sites
        const site1Name = reaction.reactants[0].bindingSiteName;
        const site2Name = reaction.reactants[1].bindingSiteName;
        
        const bindingSite1 = molecule1.getBindingSite(site1Name);
        const bindingSite2 = molecule2.getBindingSite(site2Name);
        
        if (!bindingSite1 || !bindingSite2) {
            this.drawTextMessage(gl, "Binding sites not found");
            return;
        }
        
        // Get reaction geometry
        const sigma = reaction.geometry.sigma;
        const angles = reaction.geometry.angles;
        const norm1 = reaction.geometry.norm1;
        const norm2 = reaction.geometry.norm2;
        
        // Draw the first molecule
        gl.pushMatrix();
        
        // Draw center of mass of first molecule
        this.setColor(gl, 1.0, 1.0, 1.0); // White
        gl.callList(this.sphereDisplayList);
        
        // Draw binding site
        gl.pushMatrix();
        gl.translatef(bindingSite1.position.x, bindingSite1.position.y, bindingSite1.position.z);
        this.setColor(gl, 0.0, 0.0, 0.0); // Black
        gl.callList(this.sphereDisplayList);
        gl.popMatrix();
        
        // Draw connection from center to binding site
        this.drawCylinderBetweenPoints(
            gl, 
            new Vector3(0, 0, 0), 
            bindingSite1.position, 
            0.25, 
            [1.0, 1.0, 0.0] // Yellow
        );
        
        gl.popMatrix();
        
        // Draw the second molecule
        gl.pushMatrix();
        
        // Position the second molecule at distance sigma along x-axis
        gl.translatef(sigma, 0, 0);
        
        // Draw center of mass of second molecule
        this.setColor(gl, 1.0, 1.0, 1.0); // White
        gl.callList(this.sphereDisplayList);
        
        // Draw binding site
        gl.pushMatrix();
        gl.translatef(bindingSite2.position.x, bindingSite2.position.y, bindingSite2.position.z);
        this.setColor(gl, 0.0, 0.0, 0.0); // Black
        gl.callList(this.sphereDisplayList);
        gl.popMatrix();
        
        // Draw connection from center to binding site
        this.drawCylinderBetweenPoints(
            gl, 
            new Vector3(0, 0, 0), 
            bindingSite2.position, 
            0.25, 
            [1.0, 1.0, 0.0] // Yellow
        );
        
        gl.popMatrix();
        
        // Draw dashed line representing the sigma distance
        this.drawDashedLine(
            gl,
            new Vector3(0, 0, 0),
            new Vector3(sigma, 0, 0),
            0.0, 0.0, 0.0 // Black
        );
        
        // Display angles
        this.drawAnglesText(gl, angles);
    }
    
    /**
     * Draw a cylinder between two points
     * @param {WebGLRenderingContext} gl - WebGL context
     * @param {Vector3} start - Start point
     * @param {Vector3} end - End point
     * @param {number} radius - Cylinder radius
     * @param {number[]} color - RGB color array [r, g, b]
     * @private
     */
    drawCylinderBetweenPoints(gl, start, end, radius, color) {
        gl.pushMatrix();
        
        // Calculate direction and length
        const direction = end.clone().subtract(start);
        const length = direction.length();
        direction.normalize();
        
        // Default cylinder is along z-axis, need to rotate it to align with direction
        const zAxis = new Vector3(0, 0, 1);
        const rotationAxis = new Vector3().copy(zAxis).cross(direction);
        const angle = Math.acos(zAxis.dot(direction)) * 180 / Math.PI;
        
        // Translate to start position
        gl.translatef(start.x, start.y, start.z);
        
        // Apply rotation and scale to align with binding site
        if (rotationAxis.length() > 0.001) {
            gl.rotatef(angle, rotationAxis.x, rotationAxis.y, rotationAxis.z);
        }
        
        // Set color
        this.setColor(gl, color[0], color[1], color[2]);
        
        // Scale cylinder to the required length and radius
        gl.scalef(radius, radius, length);
        
        // Draw the cylinder
        gl.callList(this.cylinderDisplayList);
        
        gl.popMatrix();
    }
    
    /**
     * Draw a dashed line between two points
     * @param {WebGLRenderingContext} gl - WebGL context
     * @param {Vector3} start - Start point
     * @param {Vector3} end - End point
     * @param {number} r - Red component
     * @param {number} g - Green component
     * @param {number} b - Blue component
     * @private
     */
    drawDashedLine(gl, start, end, r, g, b) {
        gl.pushMatrix();
        
        gl.enable(gl.LINE_STIPPLE);
        const factor = 1;
        const pattern = 0x5555; // Dashed pattern
        gl.lineStipple(factor, pattern);
        gl.lineWidth(4.0);
        
        this.setColor(gl, r, g, b);
        
        gl.begin(gl.LINES);
        gl.vertex3f(start.x, start.y, start.z);
        gl.vertex3f(end.x, end.y, end.z);
        gl.end();
        
        gl.disable(gl.LINE_STIPPLE);
        
        gl.popMatrix();
    }
    
    /**
     * Draw text displaying angles
     * @param {WebGLRenderingContext} gl - WebGL context
     * @param {number[]} angles - Array of angles [theta1, theta2, phi1, phi2, omega]
     * @private
     */
    drawAnglesText(gl, angles) {
        // Disable lighting for text rendering
        gl.disable(gl.LIGHTING);
        
        // Set up orthographic projection for text
        gl.matrixMode(gl.PROJECTION);
        gl.pushMatrix();
        gl.loadIdentity();
        gl.ortho(0, this.reactionCanvas.width, 0, this.reactionCanvas.height, -1, 1);
        
        gl.matrixMode(gl.MODELVIEW);
        gl.pushMatrix();
        gl.loadIdentity();
        
        // Set text color
        gl.color3f(0.0, 0.0, 0.0); // Black
        
        // Position for text
        const x = 10;
        let y = this.reactionCanvas.height - 20;
        
        // Draw angle values
        this.drawText(gl, `Theta 1: ${angles[0].toFixed(1)}°`, x, y);
        y -= 20;
        this.drawText(gl, `Theta 2: ${angles[1].toFixed(1)}°`, x, y);
        y -= 20;
        this.drawText(gl, `Phi 1: ${angles[2] === -999 ? '-' : angles[2].toFixed(1) + '°'}`, x, y);
        y -= 20;
        this.drawText(gl, `Phi 2: ${angles[3] === -999 ? '-' : angles[3].toFixed(1) + '°'}`, x, y);
        y -= 20;
        this.drawText(gl, `Omega: ${angles[4].toFixed(1)}°`, x, y);
        
        // Restore matrices
        gl.matrixMode(gl.PROJECTION);
        gl.popMatrix();
        
        gl.matrixMode(gl.MODELVIEW);
        gl.popMatrix();
        
        // Re-enable lighting
        gl.enable(gl.LIGHTING);
    }
    
    /**
     * Draw text in the canvas
     * @param {WebGLRenderingContext} gl - WebGL context
     * @param {string} text - Text to draw
     * @param {number} x - X position
     * @param {number} y - Y position
     * @private
     */
    drawText(gl, text, x, y) {
        // Since WebGL doesn't have direct text rendering,
        // we'll use a 2D canvas overlay for this in a real implementation
        // This is a placeholder
        console.log(`Draw text: ${text} at (${x}, ${y})`);
    }
    
    /**
     * Draw a message in the center of the canvas
     * @param {WebGLRenderingContext} gl - WebGL context
     * @param {string} message - Message to display
     * @private
     */
    drawTextMessage(gl, message) {
        // Disable lighting for text rendering
        gl.disable(gl.LIGHTING);
        
        // Set up orthographic projection for text
        gl.matrixMode(gl.PROJECTION);
        gl.pushMatrix();
        gl.loadIdentity();
        gl.ortho(0, gl.canvas.width, 0, gl.canvas.height, -1, 1);
        
        gl.matrixMode(gl.MODELVIEW);
        gl.pushMatrix();
        gl.loadIdentity();
        
        // Set text color
        gl.color3f(0.0, 0.0, 0.0); // Black
        
        // Position in center of canvas
        const x = gl.canvas.width / 2;
        const y = gl.canvas.height / 2;
        
        // Draw message
        this.drawText(gl, message, x, y);
        
        // Restore matrices
        gl.matrixMode(gl.PROJECTION);
        gl.popMatrix();
        
        gl.matrixMode(gl.MODELVIEW);
        gl.popMatrix();
        
        // Re-enable lighting
        gl.enable(gl.LIGHTING);
    }
    
    /**
     * Draw coordinate axes
     * @param {WebGLRenderingContext} gl - WebGL context
     * @private
     */
    drawAxes(gl) {
        // Draw X axis (red)
        gl.pushMatrix();
        gl.rotatef(90, 0, 1, 0);
        this.setColor(gl, 1.0, 0.0, 0.0);
        this.drawAxis(gl);
        gl.popMatrix();
        
        // Draw Y axis (green)
        gl.pushMatrix();
        gl.rotatef(-90, 1, 0, 0);
        this.setColor(gl, 0.0, 1.0, 0.0);
        this.drawAxis(gl);
        gl.popMatrix();
        
        // Draw Z axis (blue)
        gl.pushMatrix();
        this.setColor(gl, 0.0, 0.0, 1.0);
        this.drawAxis(gl);
        gl.popMatrix();
    }
    
    /**
     * Draw a single coordinate axis with arrow
     * @param {WebGLRenderingContext} gl - WebGL context
     * @private
     */
    drawAxis(gl) {
        // Draw cylinder for axis body
        gl.pushMatrix();
        gl.scalef(0.2, 0.2, 4); // Length 4, radius 0.2
        gl.callList(this.cylinderDisplayList);
        gl.popMatrix();
        
        // Draw cone for arrow head
        gl.pushMatrix();
        gl.translatef(0, 0, 4); // Position at end of cylinder
        gl.scalef(0.45, 0.45, 0.55); // Scale cone
        this.shapeGenerator.drawCone(0.5, 1.0, 12, 1);
        gl.popMatrix();
    }
    
    /**
     * Set the current color for rendering
     * @param {WebGLRenderingContext} gl - WebGL context
     * @param {number} r - Red component (0-1)
     * @param {number} g - Green component (0-1)
     * @param {number} b - Blue component (0-1)
     * @private
     */
    setColor(gl, r, g, b) {
        gl.color3f(r, g, b);
    }
    
    /**
     * Refresh the molecule view
     */
    refreshMoleculeView() {
        if (this.moleculeCanvas) {
            this.renderMoleculeView();
        }
    }
    
    /**
     * Refresh the reaction view
     */
    refreshReactionView() {
        if (this.reactionCanvas) {
            this.renderReactionView();
        }
    }
    
    /**
     * Resize the canvases when the window size changes
     */
    resizeCanvases() {
        if (this.moleculeCanvas) {
            const container = this.moleculeCanvas.parentElement;
            if (container) {
                this.moleculeCanvas.width = container.clientWidth;
                this.moleculeCanvas.height = container.clientHeight;
                this.moleculeContext.viewport(0, 0, this.moleculeCanvas.width, this.moleculeCanvas.height);
                this.refreshMoleculeView();
            }
        }
        
        if (this.reactionCanvas) {
            const container = this.reactionCanvas.parentElement;
            if (container) {
                this.reactionCanvas.width = container.clientWidth;
                this.reactionCanvas.height = container.clientHeight;
                this.reactionContext.viewport(0, 0, this.reactionCanvas.width, this.reactionCanvas.height);
                this.refreshReactionView();
            }
        }
    }
    
    /**
     * Display an error message when WebGL is not supported
     * @param {HTMLCanvasElement} canvas - Canvas element
     * @private
     */
    showWebGLError(canvas) {
        // Create a 2D context for error message
        const ctx = canvas.getContext('2d');
        if (!ctx) return;
        
        // Clear canvas
        ctx.fillStyle = 'white';
        ctx.fillRect(0, 0, canvas.width, canvas.height);
        
        // Draw error message
        ctx.fillStyle = 'red';
        ctx.font = '16px Arial';
        ctx.textAlign = 'center';
        ctx.textBaseline = 'middle';
        ctx.fillText('WebGL not supported by your browser.', canvas.width / 2, canvas.height / 2 - 20);
        ctx.fillText('Please try a different browser or update your graphics drivers.', canvas.width / 2, canvas.height / 2 + 20);
    }
    
    /**
     * Clean up resources when the renderer is no longer needed
     */
    dispose() {
        // Stop render loops
        this.stopRenderLoops();
        
        // Clean up event listeners
        if (this.moleculeCanvas) {
            this.moleculeCanvas.removeEventListener('camera-update', this.refreshMoleculeView);
        }
        
        if (this.reactionCanvas) {
            this.reactionCanvas.removeEventListener('camera-update', this.refreshReactionView);
        }
        
        // Delete WebGL resources
        if (this.moleculeContext) {
            this.deleteWebGLResources(this.moleculeContext);
        }
        
        if (this.reactionContext) {
            this.deleteWebGLResources(this.reactionContext);
        }
    }
    
    /**
     * Delete WebGL resources
     * @param {WebGLRenderingContext} gl - WebGL context
     * @private
     */
    deleteWebGLResources(gl) {
        // Delete buffers
        if (this.sphereVertexBuffer) gl.deleteBuffer(this.sphereVertexBuffer);
        if (this.sphereNormalBuffer) gl.deleteBuffer(this.sphereNormalBuffer);
        if (this.cylinderVertexBuffer) gl.deleteBuffer(this.cylinderVertexBuffer);
        if (this.cylinderNormalBuffer) gl.deleteBuffer(this.cylinderNormalBuffer);
        
        // Delete display lists
        if (this.sphereDisplayList) gl.deleteList(this.sphereDisplayList);
        if (this.cylinderDisplayList) gl.deleteList(this.cylinderDisplayList);
    }
}

export { Renderer };