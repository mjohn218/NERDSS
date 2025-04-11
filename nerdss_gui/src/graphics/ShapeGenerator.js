/**
 * ShapeGenerator.js
 * 
 * This file defines the ShapeGenerator class, which provides methods for
 * generating 3D shapes for the NERDSS GUI. It includes primitives like
 * spheres, cylinders, and cones for visualizing molecular structures.
 */

/**
 * Class for generating 3D shapes for WebGL rendering
 */
class ShapeGenerator {
    /**
     * Create a new shape generator
     * @param {WebGLRenderingContext} gl - WebGL rendering context
     */
    constructor(gl) {
        this.gl = gl;
    }
    
    /**
     * Draw a UV sphere (sphere with texture coordinates)
     * @param {number} radius - Sphere radius
     * @param {number} slices - Number of slices (longitude divisions)
     * @param {number} stacks - Number of stacks (latitude divisions)
     * @param {boolean} [generateTexCoords=false] - Whether to generate texture coordinates
     */
    drawSphere(radius, slices, stacks, generateTexCoords = false) {
        if (radius <= 0) {
            throw new Error('Radius must be positive');
        }
        if (slices < 3) {
            throw new Error('Number of slices must be at least 3');
        }
        if (stacks < 2) {
            throw new Error('Number of stacks must be at least 2');
        }
        
        const gl = this.gl;
        
        for (let j = 0; j < stacks; j++) {
            const latitude1 = (Math.PI / stacks) * j - Math.PI / 2;
            const latitude2 = (Math.PI / stacks) * (j + 1) - Math.PI / 2;
            const sinLat1 = Math.sin(latitude1);
            const cosLat1 = Math.cos(latitude1);
            const sinLat2 = Math.sin(latitude2);
            const cosLat2 = Math.cos(latitude2);
            
            gl.begin(gl.QUAD_STRIP);
            
            for (let i = 0; i <= slices; i++) {
                const longitude = (2 * Math.PI / slices) * i;
                const sinLong = Math.sin(longitude);
                const cosLong = Math.cos(longitude);
                
                const x1 = cosLong * cosLat1;
                const y1 = sinLong * cosLat1;
                const z1 = sinLat1;
                
                const x2 = cosLong * cosLat2;
                const y2 = sinLong * cosLat2;
                const z2 = sinLat2;
                
                // First vertex (bottom of quad)
                gl.normal3f(x1, y1, z1);
                if (generateTexCoords) {
                    gl.texCoord2f(i / slices, j / stacks);
                }
                gl.vertex3f(radius * x1, radius * y1, radius * z1);
                
                // Second vertex (top of quad)
                gl.normal3f(x2, y2, z2);
                if (generateTexCoords) {
                    gl.texCoord2f(i / slices, (j + 1) / stacks);
                }
                gl.vertex3f(radius * x2, radius * y2, radius * z2);
            }
            
            gl.end();
        }
    }
    
    /**
     * Generate vertex and normal data for a sphere
     * @param {number} radius - Sphere radius
     * @param {number} slices - Number of slices
     * @param {number} stacks - Number of stacks
     * @returns {Object} Object with vertices and normals arrays
     */
    generateSphereData(radius, slices, stacks) {
        const vertices = [];
        const normals = [];
        
        for (let j = 0; j < stacks; j++) {
            const latitude1 = (Math.PI / stacks) * j - Math.PI / 2;
            const latitude2 = (Math.PI / stacks) * (j + 1) - Math.PI / 2;
            const sinLat1 = Math.sin(latitude1);
            const cosLat1 = Math.cos(latitude1);
            const sinLat2 = Math.sin(latitude2);
            const cosLat2 = Math.cos(latitude2);
            
            for (let i = 0; i <= slices; i++) {
                const longitude = (2 * Math.PI / slices) * i;
                const sinLong = Math.sin(longitude);
                const cosLong = Math.cos(longitude);
                
                // First vertex (bottom of quad)
                const x1 = cosLong * cosLat1;
                const y1 = sinLong * cosLat1;
                const z1 = sinLat1;
                
                // Second vertex (top of quad)
                const x2 = cosLong * cosLat2;
                const y2 = sinLong * cosLat2;
                const z2 = sinLat2;
                
                // Add first vertex
                normals.push(x1, y1, z1);
                vertices.push(radius * x1, radius * y1, radius * z1);
                
                // Add second vertex
                normals.push(x2, y2, z2);
                vertices.push(radius * x2, radius * y2, radius * z2);
            }
        }
        
        return { vertices, normals };
    }
    
    /**
     * Draw a cylinder
     * @param {number} radius - Cylinder radius
     * @param {number} height - Cylinder height
     * @param {number} slices - Number of slices (divisions around circumference)
     * @param {number} [stacks=1] - Number of stacks (divisions along height)
     * @param {boolean} [caps=true] - Whether to include end caps
     * @param {boolean} [generateTexCoords=false] - Whether to generate texture coordinates
     */
    drawCylinder(radius, height, slices, stacks = 1, caps = true, generateTexCoords = false) {
        if (radius <= 0) {
            throw new Error('Radius must be positive');
        }
        if (height <= 0) {
            throw new Error('Height must be positive');
        }
        if (slices < 3) {
            throw new Error('Number of slices must be at least 3');
        }
        if (stacks < 1) {
            throw new Error('Number of stacks must be at least 1');
        }
        
        const gl = this.gl;
        
        // Draw the cylinder body
        for (let j = 0; j < stacks; j++) {
            const z1 = height * j / stacks;
            const z2 = height * (j + 1) / stacks;
            
            gl.begin(gl.QUAD_STRIP);
            
            for (let i = 0; i <= slices; i++) {
                const angle = 2 * Math.PI * i / slices;
                const sin = Math.sin(angle);
                const cos = Math.cos(angle);
                
                // Normal points outward from cylinder axis
                gl.normal3f(cos, sin, 0);
                
                // First vertex (bottom of quad)
                if (generateTexCoords) {
                    gl.texCoord2f(i / slices, j / stacks);
                }
                gl.vertex3f(radius * cos, radius * sin, z1);
                
                // Second vertex (top of quad)
                if (generateTexCoords) {
                    gl.texCoord2f(i / slices, (j + 1) / stacks);
                }
                gl.vertex3f(radius * cos, radius * sin, z2);
            }
            
            gl.end();
        }
        
        // Draw caps if requested
        if (caps) {
            // Draw bottom cap (z = 0)
            gl.begin(gl.TRIANGLE_FAN);
            gl.normal3f(0, 0, -1);
            gl.vertex3f(0, 0, 0); // Center of cap
            
            for (let i = 0; i <= slices; i++) {
                const angle = 2 * Math.PI * i / slices;
                const sin = Math.sin(angle);
                const cos = Math.cos(angle);
                
                if (generateTexCoords) {
                    gl.texCoord2f(0.5 + 0.5 * cos, 0.5 + 0.5 * sin);
                }
                gl.vertex3f(radius * cos, radius * sin, 0);
            }
            
            gl.end();
            
            // Draw top cap (z = height)
            gl.begin(gl.TRIANGLE_FAN);
            gl.normal3f(0, 0, 1);
            gl.vertex3f(0, 0, height); // Center of cap
            
            for (let i = slices; i >= 0; i--) {
                const angle = 2 * Math.PI * i / slices;
                const sin = Math.sin(angle);
                const cos = Math.cos(angle);
                
                if (generateTexCoords) {
                    gl.texCoord2f(0.5 + 0.5 * cos, 0.5 + 0.5 * sin);
                }
                gl.vertex3f(radius * cos, radius * sin, height);
            }
            
            gl.end();
        }
    }
    
    /**
     * Generate vertex and normal data for a cylinder
     * @param {number} radius - Cylinder radius
     * @param {number} height - Cylinder height
     * @param {number} slices - Number of slices
     * @param {number} [stacks=1] - Number of stacks
     * @param {boolean} [caps=true] - Whether to include end caps
     * @returns {Object} Object with vertices and normals arrays
     */
    generateCylinderData(radius, height, slices, stacks = 1, caps = true) {
        const vertices = [];
        const normals = [];
        
        // Generate the cylinder body
        for (let j = 0; j < stacks; j++) {
            const z1 = height * j / stacks;
            const z2 = height * (j + 1) / stacks;
            
            for (let i = 0; i <= slices; i++) {
                const angle = 2 * Math.PI * i / slices;
                const sin = Math.sin(angle);
                const cos = Math.cos(angle);
                
                // Normal for both vertices
                const nx = cos;
                const ny = sin;
                const nz = 0;
                
                // First vertex (bottom of quad)
                normals.push(nx, ny, nz);
                vertices.push(radius * cos, radius * sin, z1);
                
                // Second vertex (top of quad)
                normals.push(nx, ny, nz);
                vertices.push(radius * cos, radius * sin, z2);
            }
        }
        
        // Generate caps if requested
        if (caps) {
            // Bottom cap (z = 0)
            // Center vertex
            normals.push(0, 0, -1);
            vertices.push(0, 0, 0);
            
            // Rim vertices
            for (let i = 0; i <= slices; i++) {
                const angle = 2 * Math.PI * i / slices;
                const sin = Math.sin(angle);
                const cos = Math.cos(angle);
                
                normals.push(0, 0, -1);
                vertices.push(radius * cos, radius * sin, 0);
            }
            
            // Top cap (z = height)
            // Center vertex
            normals.push(0, 0, 1);
            vertices.push(0, 0, height);
            
            // Rim vertices
            for (let i = 0; i <= slices; i++) {
                const angle = 2 * Math.PI * (slices - i) / slices;
                const sin = Math.sin(angle);
                const cos = Math.cos(angle);
                
                normals.push(0, 0, 1);
                vertices.push(radius * cos, radius * sin, height);
            }
        }
        
        return { vertices, normals };
    }
    
    /**
     * Draw a cone
     * @param {number} radius - Base radius
     * @param {number} height - Cone height
     * @param {number} slices - Number of slices (divisions around circumference)
     * @param {number} [stacks=1] - Number of stacks (divisions along height)
     * @param {boolean} [base=true] - Whether to include the base cap
     * @param {boolean} [generateTexCoords=false] - Whether to generate texture coordinates
     */
    drawCone(radius, height, slices, stacks = 1, base = true, generateTexCoords = false) {
        if (radius <= 0) {
            throw new Error('Radius must be positive');
        }
        if (height <= 0) {
            throw new Error('Height must be positive');
        }
        if (slices < 3) {
            throw new Error('Number of slices must be at least 3');
        }
        if (stacks < 1) {
            throw new Error('Number of stacks must be at least 1');
        }
        
        const gl = this.gl;
        
        // Draw the cone body
        for (let j = 0; j < stacks; j++) {
            const z1 = height * j / stacks;
            const z2 = height * (j + 1) / stacks;
            const r1 = radius * (1 - z1 / height);
            const r2 = radius * (1 - z2 / height);
            
            gl.begin(gl.QUAD_STRIP);
            
            for (let i = 0; i <= slices; i++) {
                const angle = 2 * Math.PI * i / slices;
                const sin = Math.sin(angle);
                const cos = Math.cos(angle);
                
                // Calculate normal for the cone surface
                // This is more complex than for a cylinder because the radius changes
                const nx1 = cos * height;
                const ny1 = sin * height;
                const nz1 = radius;
                const len1 = Math.sqrt(nx1 * nx1 + ny1 * ny1 + nz1 * nz1);
                
                // First vertex (bottom of quad)
                gl.normal3f(nx1 / len1, ny1 / len1, nz1 / len1);
                if (generateTexCoords) {
                    gl.texCoord2f(i / slices, j / stacks);
                }
                gl.vertex3f(r1 * cos, r1 * sin, z1);
                
                // Second vertex (top of quad)
                if (generateTexCoords) {
                    gl.texCoord2f(i / slices, (j + 1) / stacks);
                }
                gl.vertex3f(r2 * cos, r2 * sin, z2);
            }
            
            gl.end();
        }
        
        // Draw base if requested
        if (base) {
            gl.begin(gl.TRIANGLE_FAN);
            gl.normal3f(0, 0, -1);
            gl.vertex3f(0, 0, 0); // Center of base
            
            for (let i = 0; i <= slices; i++) {
                const angle = 2 * Math.PI * i / slices;
                const sin = Math.sin(angle);
                const cos = Math.cos(angle);
                
                if (generateTexCoords) {
                    gl.texCoord2f(0.5 + 0.5 * cos, 0.5 + 0.5 * sin);
                }
                gl.vertex3f(radius * cos, radius * sin, 0);
            }
            
            gl.end();
        }
    }
    
    /**
     * Draw a square in the XY plane
     * @param {number} size - Side length
     * @param {boolean} [generateTexCoords=false] - Whether to generate texture coordinates
     */
    drawSquare(size, generateTexCoords = false) {
        const gl = this.gl;
        const halfSize = size / 2;
        
        gl.begin(gl.POLYGON);
        gl.normal3f(0, 0, 1);
        
        if (generateTexCoords) gl.texCoord2f(0, 0);
        gl.vertex2f(-halfSize, -halfSize);
        
        if (generateTexCoords) gl.texCoord2f(1, 0);
        gl.vertex2f(halfSize, -halfSize);
        
        if (generateTexCoords) gl.texCoord2f(1, 1);
        gl.vertex2f(halfSize, halfSize);
        
        if (generateTexCoords) gl.texCoord2f(0, 1);
        gl.vertex2f(-halfSize, halfSize);
        
        gl.end();
    }
    
    /**
     * Draw a circle in the XY plane
     * @param {number} radius - Circle radius
     * @param {number} slices - Number of slices (divisions around circumference)
     * @param {boolean} [generateTexCoords=false] - Whether to generate texture coordinates
     */
    drawCircle(radius, slices, generateTexCoords = false) {
        if (radius <= 0) {
            throw new Error('Radius must be positive');
        }
        if (slices < 3) {
            throw new Error('Number of slices must be at least 3');
        }
        
        const gl = this.gl;
        
        gl.begin(gl.TRIANGLE_FAN);
        gl.normal3f(0, 0, 1);
        gl.vertex3f(0, 0, 0); // Center of circle
        
        for (let i = 0; i <= slices; i++) {
            const angle = 2 * Math.PI * i / slices;
            const sin = Math.sin(angle);
            const cos = Math.cos(angle);
            
            if (generateTexCoords) {
                gl.texCoord2f(0.5 + 0.5 * cos, 0.5 + 0.5 * sin);
            }
            gl.vertex3f(radius * cos, radius * sin, 0);
        }
        
        gl.end();
    }
    
    /**
     * Draw a torus (donut shape)
     * @param {number} outerRadius - Outer radius of the torus
     * @param {number} innerRadius - Inner radius of the torus (tube radius)
     * @param {number} sides - Number of sides around the tube
     * @param {number} rings - Number of rings around the torus
     * @param {boolean} [generateTexCoords=false] - Whether to generate texture coordinates
     */
    drawTorus(outerRadius, innerRadius, sides, rings, generateTexCoords = false) {
        if (outerRadius <= innerRadius) {
            throw new Error('Outer radius must be greater than inner radius');
        }
        if (innerRadius <= 0) {
            throw new Error('Inner radius must be positive');
        }
        if (sides < 3) {
            throw new Error('Number of sides must be at least 3');
        }
        if (rings < 3) {
            throw new Error('Number of rings must be at least 3');
        }
        
        const gl = this.gl;
        const centerRadius = (outerRadius + innerRadius) / 2;
        const tubeRadius = outerRadius - centerRadius;
        
        for (let i = 0; i < rings; i++) {
            const ringAngle1 = 2 * Math.PI * i / rings;
            const ringAngle2 = 2 * Math.PI * (i + 1) / rings;
            const ringCos1 = Math.cos(ringAngle1);
            const ringSin1 = Math.sin(ringAngle1);
            const ringCos2 = Math.cos(ringAngle2);
            const ringSin2 = Math.sin(ringAngle2);
            
            gl.begin(gl.QUAD_STRIP);
            
            for (let j = 0; j <= sides; j++) {
                const sideAngle = 2 * Math.PI * j / sides;
                const sideCos = Math.cos(sideAngle);
                const sideSin = Math.sin(sideAngle);
                
                // First vertex
                const nx1 = ringCos1 * sideCos;
                const ny1 = ringSin1 * sideCos;
                const nz1 = sideSin;
                
                const x1 = ringCos1 * (centerRadius + tubeRadius * sideCos);
                const y1 = ringSin1 * (centerRadius + tubeRadius * sideCos);
                const z1 = tubeRadius * sideSin;
                
                gl.normal3f(nx1, ny1, nz1);
                if (generateTexCoords) {
                    gl.texCoord2f(i / rings, j / sides);
                }
                gl.vertex3f(x1, y1, z1);
                
                // Second vertex
                const nx2 = ringCos2 * sideCos;
                const ny2 = ringSin2 * sideCos;
                const nz2 = sideSin;
                
                const x2 = ringCos2 * (centerRadius + tubeRadius * sideCos);
                const y2 = ringSin2 * (centerRadius + tubeRadius * sideCos);
                const z2 = tubeRadius * sideSin;
                
                gl.normal3f(nx2, ny2, nz2);
                if (generateTexCoords) {
                    gl.texCoord2f((i + 1) / rings, j / sides);
                }
                gl.vertex3f(x2, y2, z2);
            }
            
            gl.end();
        }
    }
    
    /**
     * Draw a cube
     * @param {number} size - Side length
     * @param {boolean} [generateTexCoords=false] - Whether to generate texture coordinates
     */
    drawCube(size, generateTexCoords = false) {
        const gl = this.gl;
        const halfSize = size / 2;
        
        // Front face (z = halfSize)
        gl.begin(gl.QUADS);
        gl.normal3f(0, 0, 1);
        if (generateTexCoords) gl.texCoord2f(0, 0);
        gl.vertex3f(-halfSize, -halfSize, halfSize);
        if (generateTexCoords) gl.texCoord2f(1, 0);
        gl.vertex3f(halfSize, -halfSize, halfSize);
        if (generateTexCoords) gl.texCoord2f(1, 1);
        gl.vertex3f(halfSize, halfSize, halfSize);
        if (generateTexCoords) gl.texCoord2f(0, 1);
        gl.vertex3f(-halfSize, halfSize, halfSize);
        gl.end();
        
        // Back face (z = -halfSize)
        gl.begin(gl.QUADS);
        gl.normal3f(0, 0, -1);
        if (generateTexCoords) gl.texCoord2f(0, 0);
        gl.vertex3f(-halfSize, -halfSize, -halfSize);
        if (generateTexCoords) gl.texCoord2f(0, 1);
        gl.vertex3f(-halfSize, halfSize, -halfSize);
        if (generateTexCoords) gl.texCoord2f(1, 1);
        gl.vertex3f(halfSize, halfSize, -halfSize);
        if (generateTexCoords) gl.texCoord2f(1, 0);
        gl.vertex3f(halfSize, -halfSize, -halfSize);
        gl.end();
        
        // Left face (x = -halfSize)
        gl.begin(gl.QUADS);
        gl.normal3f(-1, 0, 0);
        if (generateTexCoords) gl.texCoord2f(0, 0);
        gl.vertex3f(-halfSize, -halfSize, -halfSize);
        if (generateTexCoords) gl.texCoord2f(1, 0);
        gl.vertex3f(-halfSize, -halfSize, halfSize);
        if (generateTexCoords) gl.texCoord2f(1, 1);
        gl.vertex3f(-halfSize, halfSize, halfSize);
        if (generateTexCoords) gl.texCoord2f(0, 1);
        gl.vertex3f(-halfSize, halfSize, -halfSize);
        gl.end();
        
        // Right face (x = halfSize)
        gl.begin(gl.QUADS);
        gl.normal3f(1, 0, 0);
        if (generateTexCoords) gl.texCoord2f(0, 0);
        gl.vertex3f(halfSize, -halfSize, -halfSize);
        if (generateTexCoords) gl.texCoord2f(0, 1);
        gl.vertex3f(halfSize, halfSize, -halfSize);
        if (generateTexCoords) gl.texCoord2f(1, 1);
        gl.vertex3f(halfSize, halfSize, halfSize);
        if (generateTexCoords) gl.texCoord2f(1, 0);
        gl.vertex3f(halfSize, -halfSize, halfSize);
        gl.end();
        
        // Top face (y = halfSize)
        gl.begin(gl.QUADS);
        gl.normal3f(0, 1, 0);
        if (generateTexCoords) gl.texCoord2f(0, 0);
        gl.vertex3f(-halfSize, halfSize, -halfSize);
        if (generateTexCoords) gl.texCoord2f(0, 1);
        gl.vertex3f(-halfSize, halfSize, halfSize);
        if (generateTexCoords) gl.texCoord2f(1, 1);
        gl.vertex3f(halfSize, halfSize, halfSize);
        if (generateTexCoords) gl.texCoord2f(1, 0);
        gl.vertex3f(halfSize, halfSize, -halfSize);
        gl.end();
        
        // Bottom face (y = -halfSize)
        gl.begin(gl.QUADS);
        gl.normal3f(0, -1, 0);
        if (generateTexCoords) gl.texCoord2f(0, 0);
        gl.vertex3f(-halfSize, -halfSize, -halfSize);
        if (generateTexCoords) gl.texCoord2f(1, 0);
        gl.vertex3f(halfSize, -halfSize, -halfSize);
        if (generateTexCoords) gl.texCoord2f(1, 1);
        gl.vertex3f(halfSize, -halfSize, halfSize);
        if (generateTexCoords) gl.texCoord2f(0, 1);
        gl.vertex3f(-halfSize, -halfSize, halfSize);
        gl.end();
    }
}

export { ShapeGenerator };