/**
 * EventEmitter.js
 * 
 * A simple event emitter implementation for handling custom events
 * in the application. This follows the observer pattern to allow
 * different components to subscribe to and react to events.
 */

/**
 * Event emitter class for implementing the observer pattern
 */
class EventEmitter {
    /**
     * Create a new event emitter
     */
    constructor() {
        // Map of event types to arrays of listener functions
        this.events = new Map();
    }
    
    /**
     * Register an event listener
     * @param {string} eventName - Name of the event to listen for
     * @param {Function} listener - Function to call when event is triggered
     * @returns {EventEmitter} This emitter for chaining
     */
    on(eventName, listener) {
        if (!this.events.has(eventName)) {
            this.events.set(eventName, []);
        }
        
        this.events.get(eventName).push(listener);
        return this;
    }
    
    /**
     * Register a one-time event listener
     * @param {string} eventName - Name of the event to listen for
     * @param {Function} listener - Function to call when event is triggered
     * @returns {EventEmitter} This emitter for chaining
     */
    once(eventName, listener) {
        const onceWrapper = (...args) => {
            this.off(eventName, onceWrapper);
            listener.apply(this, args);
        };
        
        // Store the original listener to allow for removal
        onceWrapper.originalListener = listener;
        
        return this.on(eventName, onceWrapper);
    }
    
    /**
     * Remove an event listener
     * @param {string} eventName - Name of the event
     * @param {Function} listener - Listener function to remove
     * @returns {EventEmitter} This emitter for chaining
     */
    off(eventName, listener) {
        if (!this.events.has(eventName)) {
            return this;
        }
        
        const listeners = this.events.get(eventName);
        const index = listeners.findIndex(l => 
            l === listener || 
            (l.originalListener && l.originalListener === listener)
        );
        
        if (index !== -1) {
            listeners.splice(index, 1);
            
            // Clean up empty event arrays
            if (listeners.length === 0) {
                this.events.delete(eventName);
            }
        }
        
        return this;
    }
    
    /**
     * Remove all listeners for an event
     * @param {string} [eventName] - Name of the event (if omitted, removes all listeners for all events)
     * @returns {EventEmitter} This emitter for chaining
     */
    removeAllListeners(eventName) {
        if (eventName) {
            this.events.delete(eventName);
        } else {
            this.events.clear();
        }
        
        return this;
    }
    
    /**
     * Emit an event
     * @param {string} eventName - Name of the event to emit
     * @param {...any} args - Arguments to pass to listeners
     * @returns {boolean} True if event had listeners
     */
    emit(eventName, ...args) {
        if (!this.events.has(eventName)) {
            return false;
        }
        
        const listeners = this.events.get(eventName).slice(); // Create a copy to prevent issues if listeners modify the array
        
        for (const listener of listeners) {
            try {
                listener.apply(this, args);
            } catch (error) {
                console.error(`Error in event listener for '${eventName}':`, error);
                
                // Emit an error event
                if (eventName !== 'error') {
                    this.emit('error', error);
                }
            }
        }
        
        return listeners.length > 0;
    }
    
    /**
     * Get the number of listeners for an event
     * @param {string} eventName - Name of the event
     * @returns {number} Number of listeners
     */
    listenerCount(eventName) {
        return this.events.has(eventName) ? this.events.get(eventName).length : 0;
    }
    
    /**
     * Get all listeners for an event
     * @param {string} eventName - Name of the event
     * @returns {Function[]} Array of listener functions
     */
    listeners(eventName) {
        return this.events.has(eventName) ? this.events.get(eventName).slice() : [];
    }
    
    /**
     * Get all registered event names
     * @returns {string[]} Array of event names
     */
    eventNames() {
        return Array.from(this.events.keys());
    }
}

export { EventEmitter };