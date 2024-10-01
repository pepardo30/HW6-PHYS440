
import Foundation

// Main function to run the matrix Schrödinger equation simulation
func main() {
    print("Running Matrix Schrödinger Equation Simulation...")
    
    // Example usage: Load a potential, solve the Schrödinger equation, and print the result
    let filePath = "potential.txt" // Path to the potential file
    let potential = loadPotential(from: filePath)
    
    if potential.isEmpty {
        print("No potential data loaded.")
        return
    }
    
    // Construct the Hamiltonian matrix and solve for eigenvalues
    let dx = boxLength / Double(potential.count)
    let H = constructHamiltonian(potential: potential, dx: dx)
    let (eigenvalues, eigenvectors) = solveHamiltonian(H: H)
    
    // Print the first few eigenvalues (energy levels)
    print("Eigenvalues (Energy Levels):")
    for i in 0..<min(10, eigenvalues.count) {
        print("E[\(i)] = \(eigenvalues[i])")
    }
    
    // Print the first eigenvector (ground state wavefunction)
    print("
Ground state wavefunction:")
    for i in 0..<min(10, eigenvectors[0].count) {
        print("ψ[\(i)] = \(eigenvectors[0][i])")
    }
}

main()
