/**
 * FileManager.js
 * 
 * This file defines the FileManager class, which handles file operations
 * in the NERDSS application. It provides methods for importing and exporting
 * molecule files, parameter files, and generating model input files.
 */

import { EventEmitter } from './utils/EventEmitter.js';

/**
 * Manages file operations in the application
 * @extends EventEmitter
 */
class FileManager extends EventEmitter {
    /**
     * Create a new FileManager
     * @param {MoleculeManager} moleculeManager - Reference to the molecule manager
     * @param {ReactionManager} reactionManager - Reference to the reaction manager
     * @param {ParameterManager} parameterManager - Reference to the parameter manager
     * @param {BoundaryManager} boundaryManager - Reference to the boundary manager
     */
    constructor(moleculeManager, reactionManager, parameterManager, boundaryManager) {
        super();
        this.moleculeManager = moleculeManager;
        this.reactionManager = reactionManager;
        this.parameterManager = parameterManager;
        this.boundaryManager = boundaryManager;
    }
    
    /**
     * Set manager references
     * @param {MoleculeManager} moleculeManager - Reference to the molecule manager
     * @param {ReactionManager} reactionManager - Reference to the reaction manager
     * @param {ParameterManager} parameterManager - Reference to the parameter manager
     * @param {BoundaryManager} boundaryManager - Reference to the boundary manager
     */
    setManagers(moleculeManager, reactionManager, parameterManager, boundaryManager) {
        this.moleculeManager = moleculeManager;
        this.reactionManager = reactionManager;
        this.parameterManager = parameterManager;
        this.boundaryManager = boundaryManager;
    }
    
    /**
     * Import a molecule from a MolecX.mol file
     * @param {File|Blob} file - File or Blob object containing the MolecX.mol file
     * @returns {Promise<Object>} Promise resolving to the imported molecule or rejecting with an error
     */
    async importMolecule(file) {
        return new Promise((resolve, reject) => {
            const reader = new FileReader();
            
            reader.onload = (event) => {
                try {
                    const content = event.target.result;
                    const molecule = this.moleculeManager.importFromMolFile(content);
                    
                    if (molecule) {
                        resolve(molecule);
                    } else {
                        reject(new Error('Failed to import molecule'));
                    }
                } catch (error) {
                    reject(error);
                }
            };
            
            reader.onerror = () => {
                reject(new Error('Error reading file'));
            };
            
            reader.readAsText(file);
        });
    }
    
    /**
     * Export a molecule to a MolecX.mol file
     * @param {string} moleculeId - ID of the molecule to export
     * @returns {Blob} Blob containing the MolecX.mol file content
     */
    exportMolecule(moleculeId) {
        const content = this.moleculeManager.exportToMolFile(moleculeId);
        
        if (!content) {
            this.emit('error', 'Failed to export molecule');
            return null;
        }
        
        return new Blob([content], { type: 'text/plain' });
    }
    
    /**
     * Generate and download a model.inp file
     * @returns {Blob} Blob containing the model.inp file content
     */
    generateModelInput() {
        if (!this.parameterManager || !this.boundaryManager || 
            !this.moleculeManager || !this.reactionManager) {
            this.emit('error', 'Missing manager reference(s)');
            return null;
        }
        
        let content = "# Input file\n\n";
        
        // Add parameters section
        content += this.parameterManager.generateModelInputContent();
        content += "\n";
        
        // Add boundaries section
        content += this.boundaryManager.generateModelInputContent();
        content += "\n";
        
        // Add molecules section
        content += "start molecules\n";
        for (const molecule of this.moleculeManager.getAllMolecules()) {
            content += `\t${molecule.name}\n`;
        }
        content += "end molecules\n\n";
        
        // Add reactions section
        content += this.reactionManager.generateModelInputContent();
        
        return new Blob([content], { type: 'text/plain' });
    }
    
    /**
     * Export files to the user's file system
     * This will create a directory with all required files
     */
    exportFiles() {
        // Check if managers are available
        if (!this.moleculeManager || !this.reactionManager || 
            !this.parameterManager || !this.boundaryManager) {
            this.emit('error', 'Cannot export files without manager references');
            return;
        }
        
        // Check if file system API is available (only in secure contexts)
        if (window.showDirectoryPicker) {
            this.exportFilesUsingFileSystemAPI();
        } else {
            // Fall back to downloadable zip file or individual downloads
            this.exportFilesAsDownloads();
        }
    }
    
    /**
     * Export files using the File System Access API
     * This method is only available in supported browsers
     * @private
     */
    async exportFilesUsingFileSystemAPI() {
        try {
            // Request directory access
            const dirHandle = await window.showDirectoryPicker({
                mode: 'readwrite',
                startIn: 'downloads',
                id: 'nerdss-export'
            });
            
            // Generate model.inp file
            const modelInpContent = this.generateModelInput();
            const modelInpFile = await dirHandle.getFileHandle('model.inp', { create: true });
            const modelInpWriter = await modelInpFile.createWritable();
            await modelInpWriter.write(modelInpContent);
            await modelInpWriter.close();
            
            // Generate molecule files
            for (const molecule of this.moleculeManager.getAllMolecules()) {
                const molContent = this.moleculeManager.exportToMolFile(molecule.id);
                const molFile = await dirHandle.getFileHandle(`${molecule.name}.mol`, { create: true });
                const molWriter = await molFile.createWritable();
                await molWriter.write(new Blob([molContent], { type: 'text/plain' }));
                await molWriter.close();
            }
            
            this.emit('files-exported');
            
        } catch (error) {
            console.error('Error exporting files:', error);
            this.emit('error', `Error exporting files: ${error.message}`);
            
            // Fall back to downloads if permission was denied or API failed
            this.exportFilesAsDownloads();
        }
    }
    
    /**
     * Export files as downloads
     * This is a fallback for browsers that don't support the File System Access API
     * @private
     */
    exportFilesAsDownloads() {
        // Generate and download model.inp
        this.downloadFile(this.generateModelInput(), 'model.inp');
        
        // Download each molecule file
        for (const molecule of this.moleculeManager.getAllMolecules()) {
            const molContent = this.moleculeManager.exportToMolFile(molecule.id);
            this.downloadFile(
                new Blob([molContent], { type: 'text/plain' }), 
                `${molecule.name}.mol`
            );
        }
        
        this.emit('files-exported');
    }
    
    /**
     * Trigger a download of a file
     * @param {Blob} blob - Blob containing the file content
     * @param {string} filename - Name of the file to download
     * @private
     */
    downloadFile(blob, filename) {
        const url = URL.createObjectURL(blob);
        const a = document.createElement('a');
        a.href = url;
        a.download = filename;
        a.style.display = 'none';
        document.body.appendChild(a);
        a.click();
        
        // Clean up
        setTimeout(() => {
            document.body.removeChild(a);
            URL.revokeObjectURL(url);
        }, 100);
    }
    
    /**
     * Import a project from a JSON file
     * @param {File|Blob} file - File or Blob object containing the project JSON
     * @returns {Promise<boolean>} Promise resolving to true if import was successful
     */
    async importProject(file) {
        return new Promise((resolve, reject) => {
            const reader = new FileReader();
            
            reader.onload = (event) => {
                try {
                    const json = event.target.result;
                    const project = JSON.parse(json);
                    
                    // Import molecules
                    if (project.molecules) {
                        this.moleculeManager.fromJSON(JSON.stringify(project.molecules));
                    }
                    
                    // Import reactions
                    if (project.reactions) {
                        this.reactionManager.fromJSON(JSON.stringify(project.reactions));
                    }
                    
                    // Import parameters
                    if (project.parameters) {
                        this.parameterManager.fromJSON(JSON.stringify(project.parameters));
                    }
                    
                    // Import boundaries
                    if (project.boundaries) {
                        this.boundaryManager.fromJSON(JSON.stringify(project.boundaries));
                    }
                    
                    this.emit('project-imported');
                    resolve(true);
                    
                } catch (error) {
                    console.error('Error importing project:', error);
                    this.emit('error', `Error importing project: ${error.message}`);
                    reject(error);
                }
            };
            
            reader.onerror = () => {
                reject(new Error('Error reading file'));
            };
            
            reader.readAsText(file);
        });
    }
    
    /**
     * Export the current project as a JSON file
     * @returns {Blob} Blob containing the project JSON
     */
    exportProject() {
        const project = {
            molecules: JSON.parse(this.moleculeManager.toJSON()),
            reactions: JSON.parse(this.reactionManager.toJSON()),
            parameters: JSON.parse(this.parameterManager.toJSON()),
            boundaries: JSON.parse(this.boundaryManager.toJSON())
        };
        
        return new Blob([JSON.stringify(project, null, 2)], { type: 'application/json' });
    }
    
    /**
     * Save the current project to a file
     */
    saveProject() {
        const projectBlob = this.exportProject();
        this.downloadFile(projectBlob, 'nerdss-project.json');
    }
    
    /**
     * Load a project from a file
     * This opens a file picker to let the user select a project file
     */
    loadProject() {
        // Create a file input element
        const fileInput = document.createElement('input');
        fileInput.type = 'file';
        fileInput.accept = '.json';
        fileInput.style.display = 'none';
        document.body.appendChild(fileInput);
        
        fileInput.onchange = async (event) => {
            const file = event.target.files[0];
            if (file) {
                try {
                    await this.importProject(file);
                } catch (error) {
                    console.error('Error loading project:', error);
                    this.emit('error', `Error loading project: ${error.message}`);
                }
            }
            
            // Clean up
            document.body.removeChild(fileInput);
        };
        
        fileInput.click();
    }
}

export { FileManager };