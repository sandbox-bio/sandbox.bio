const { defineConfig } = require("cypress")
const webpackPreprocessor = require("@cypress/webpack-preprocessor")

module.exports = defineConfig({
  fixturesFolder: false,
  defaultCommandTimeout: 40000,
  video: false,
  screenshotOnRunFailure: false,
  port: 11111,
  fileServerFolder: "public/",
  e2e: {
    setupNodeEvents(on, config) {
      // Custom Webpack settings so we can import man/*.txt files without error when running tests
      on("file:preprocessor", webpackPreprocessor({
        webpackOptions: {
          module: {
            rules: [{ test: /\.txt$/, use: "raw-loader" }],
          },
        }
      }))
    },
    specPattern: "tests/**/*.cy.{js,jsx,ts,tsx}",
    supportFile: false,
  },
})
