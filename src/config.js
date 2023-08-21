// =============================================================================
// Settings
// =============================================================================
export const LOGGING = 1;
export const MAX_FILE_SIZE_TO_CACHE = 50 * 1024 * 1024; // 50MB
export const DIR_TUTORIAL = "/root/tutorial";

// =============================================================================
// Constants
// =============================================================================

// Logging
export const LOGGING_NONE = 0;
export const LOGGING_INFO = 1;
export const LOGGING_DEBUG = 2;

// Tutorials
export const DIR_TUTORIAL_SHORT = DIR_TUTORIAL.replace("/root", "~");

// v86 Emulator
export const BUS_SERIAL_INPUT = "serial0-input";
export const BUS_SERIAL_OUTPUT = "serial0-output-char";
// Use a different serial port for communication between JavaScript and v86 (UART1)
export const BUS_SERIAL_APP_INPUT = "serial1-input";
export const BUS_SERIAL_APP_OUTPUT = "serial1-output-char";
