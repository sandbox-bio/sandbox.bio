// TODO: 
// * npm run test needs npm run dev?
// * how to retrieve output of terminal once execute a command?
// * example test: "samtools view | head -n5 | wc -l > somefile"

// const baseUrl = "http://localhost:5000/";
// describe("Input validation", () => {
// 	it("Empty constructor", () => {
// 		cy.intercept("GET", "/v2/bedtools/2.29.2/bedtools.wasm").as("yes")
// 		cy.visit(baseUrl);
// 		cy.wait("@yes", { timeout: 10e3 }).then(d => {
// 			console.log("DONE")
// 			cy.get(".xterm-helper-textarea").type("pwd\n")	
// 			expect(5).to.equal(5);
// 		})
// 	});
// });
