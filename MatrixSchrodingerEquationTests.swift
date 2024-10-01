
import XCTest

// Test class for the 1D Schrödinger Equation Solver using matrix method
final class MatrixSchrodingerEquationTests: XCTestCase {

    // Test case to check if the potential is loaded correctly
    func testLoadPotential() {
        let filePath = "test_potential.txt" // Assume this file is created for testing purposes
        let potential = loadPotential(from: filePath)
        
        XCTAssertFalse(potential.isEmpty, "Potential should not be empty.")
    }

    // Test case to verify the particle in a box wavefunction at certain points
    func testParticleInBoxWavefunction() {
        let wavefunctionAtZero = particleInBoxWavefunction(n: 1, x: 0.0)
        let wavefunctionAtBoxLength = particleInBoxWavefunction(n: 1, x: boxLength)
        
        XCTAssertEqual(wavefunctionAtZero, 0.0, "The wavefunction should be 0 at x = 0.")
        XCTAssertEqual(wavefunctionAtBoxLength, 0.0, "The wavefunction should be 0 at x = L.")
    }
    
    // Test case to verify that the Hamiltonian matrix is constructed correctly
    func testConstructHamiltonian() {
        let potential = [1.0, 2.0, 1.5, 1.0, 0.5] // Simple test potential
        let dx = boxLength / Double(potential.count)
        let H = constructHamiltonian(potential: potential, dx: dx)
        
        XCTAssertEqual(H.count, numBasis, "The Hamiltonian matrix should have numBasis rows.")
        XCTAssertEqual(H[0].count, numBasis, "The Hamiltonian matrix should have numBasis columns.")
    }

    // Test case to verify that the matrix solver returns valid eigenvalues
    func testSolveHamiltonian() {
        let potential = [1.0, 2.0, 1.5, 1.0, 0.5] // Simple test potential
        let dx = boxLength / Double(potential.count)
        let H = constructHamiltonian(potential: potential, dx: dx)
        let (eigenvalues, eigenvectors) = solveHamiltonian(H: H)
        
        XCTAssertEqual(eigenvalues.count, numBasis, "There should be numBasis eigenvalues.")
        XCTAssertEqual(eigenvectors.count, numBasis, "There should be numBasis eigenvectors.")
    }
    
    // Test case to verify the eigenvalue order (eigenvalues should be in ascending order)
    func testEigenvalueOrder() {
        let potential = [1.0, 2.0, 1.5, 1.0, 0.5] // Simple test potential
        let dx = boxLength / Double(potential.count)
        let H = constructHamiltonian(potential: potential, dx: dx)
        let (eigenvalues, _) = solveHamiltonian(H: H)
        
        for i in 0..<(eigenvalues.count - 1) {
            XCTAssertLessThanOrEqual(eigenvalues[i], eigenvalues[i + 1], "Eigenvalues should be in ascending order.")
        }
    }
}

// Run the tests
MatrixSchrodingerEquationTests.defaultTestSuite.run()
