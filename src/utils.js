import localforage from "localforage";

const KEY_FILESYSTEM = "fs";

export class LocalState {
	static async getFS(tutorial) {
		// FIXME: add user info
		const key = `guest:${KEY_FILESYSTEM}:${tutorial}`;
		return (await localforage.getItem(key)) || [];
	}

	static async setFS(tutorial, value) {
		// FIXME: add user info
		const key = `guest:${KEY_FILESYSTEM}:${tutorial}`;
		return localforage.setItem(key, value);
	}
}
