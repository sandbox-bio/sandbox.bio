const { defineConfig } = require('cypress')

module.exports = defineConfig({
  fixturesFolder: false,
  defaultCommandTimeout: 40000,
  video: false,
  screenshotOnRunFailure: false,
  port: 11111,
  fileServerFolder: "public/",
  e2e: {
    setupNodeEvents(on, config) {},
    specPattern: 'tests/**/*.cy.{js,jsx,ts,tsx}',
    supportFile: false,
  },
})
