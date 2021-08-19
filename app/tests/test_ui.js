// Test that the exercise functionality works as expected

describe("Test UI", () => {
	before(() => {
		cy.visit("http://localhost:5000/tutorials?id=bedtools-intro&step=14");
		cy.get("body").get(".spinner-border").should("not.exist");
	});

	it("Exercise criteria checker", () => {
		// Initially, there should be no checkmarks
		cy.get("#tutorial-wrapper").find(".bi-check-circle-fill").should("have.length", 0);

		// Solve 1/2 of the criteria
		cy.get(".xterm-screen").type("echo 1 > notexons.bed{enter}");
		cy.get("#tutorial-wrapper").find(".bi-check-circle-fill").should("have.length", 1);

		// Undo
		cy.get(".xterm-screen").type("rm notexons.bed{enter}");
		cy.get("#tutorial-wrapper").find(".bi-check-circle-fill").should("have.length", 0);

		// Solve 2/2 of the criteria
		cy.get(".xterm-screen").type("bedtools complement -i exons.bed -g genome.txt > notexons.bed{enter}");
		cy.get("#tutorial-wrapper").find(".bi-check-circle-fill").should("have.length", 2);

		// Undo
		cy.get(".xterm-screen").type("echo 2 > notexons.bed{enter}");
		cy.get("#tutorial-wrapper").find(".bi-check-circle-fill").should("have.length", 1);
		cy.get(".xterm-screen").type("rm notexons.bed{enter}");
		cy.get("#tutorial-wrapper").find(".bi-check-circle-fill").should("have.length", 0);
	});
});
