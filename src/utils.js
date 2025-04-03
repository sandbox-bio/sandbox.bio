import { get } from "svelte/store";
import localforage from "localforage";
import { env } from "$env/dynamic/public";
import { user } from "$stores/user";


// =============================================================================
// Local state management
// =============================================================================

export const STATE_FS = "fs";
export const STATE_PLAYGROUND = "sandbox";
export const STATE_STUDIO = "studio";
export const STATE_QUIZ = "quiz";
export const STATE_IDE = "ide";

export class LocalState {
    static async getKey(state, description = "") {
        const prefix = get(user)?.email || "guest";
        return `${prefix}:${state}:${description}`;
    }

    // -------------------------------------------------------------------------
    // General state (sandbox, studio)
    // -------------------------------------------------------------------------
    static async get(state) {
        const key = await this.getKey(state);
        return await localforage.getItem(key);
    }

    static async set(state, value) {
        const key = await this.getKey(state);
        return await localforage.setItem(key, value);
    }

    // -------------------------------------------------------------------------
    // IDE
    // -------------------------------------------------------------------------
    static async getIDE(fn) {
        const key = await this.getKey(STATE_IDE, fn);
        return await localforage.getItem(key);
    }

    static async setIDE(fn, value) {
        const key = await this.getKey(STATE_IDE, fn);
        return await localforage.setItem(key, value);
    }

    // -------------------------------------------------------------------------
    // File system
    // -------------------------------------------------------------------------
    static async getFS(tutorial) {
        if (!tutorial) return [];

        const key = await this.getKey(STATE_FS, tutorial);
        return (await localforage.getItem(key)) || [];
    }

    static async setFS(tutorial, value) {
        if (!tutorial) throw "Stopped saving FS state because moved away from terminal.";

        const key = await this.getKey(STATE_FS, tutorial);
        return await localforage.setItem(key, value);
    }

    // -------------------------------------------------------------------------
    // Quiz
    // -------------------------------------------------------------------------
    static getQuizKey(tutorial, step, exerciseId) {
        return `${tutorial}-${step}-${exerciseId || "default"}`;
    }

    static async getQuiz(tutorial, step, exerciseId) {
        const key = await this.getKey(STATE_QUIZ, this.getQuizKey(tutorial, step, exerciseId));
        return (await localforage.getItem(key)) || null;
    }

    static async setQuiz(tutorial, step, exerciseId, value) {
        const key = await this.getKey(STATE_QUIZ, this.getQuizKey(tutorial, step, exerciseId));
        return await localforage.setItem(key, value);
    }
}

// =============================================================================
// Utility functions
// =============================================================================

// Get table name based on environment. Only use this for public.* tables
export function t(tableName) {
    if (env.PUBLIC_ENVIRONMENT !== "prd") return `${tableName}_stg`;
    return tableName;
}


export function strToChars(str) {
    const chars = str.split("");
    return chars.map((d) => d.charCodeAt(0));
}
