const config = {
	webServer: {
		command: "npm run preview",
		port: 4173
	},
	workers: 2,
	timeout: 60000,
	preserveOutput: "never",
	testDir: "tests",
	testMatch: /(.+\.)?(test|spec)\.[jt]s/
};

export default config;
