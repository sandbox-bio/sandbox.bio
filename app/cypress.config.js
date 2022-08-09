const { defineConfig } = require('cypress')

module.exports = defineConfig({
  fixturesFolder: false,
  defaultCommandTimeout: 100000,
  video: false,
  screenshotOnRunFailure: false,
  port: 11111,
  e2e: {
    setupNodeEvents(on, config) {},
    specPattern: 'tests/**/*.cy.{js,jsx,ts,tsx}',
    supportFile: false,
  },
})
