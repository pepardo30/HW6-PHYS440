
import Foundation
import Accelerate // For linear algebra operations

// Constants
let hbar = 1.0 // Reduced Planck constant (can be set to 1 in natural units)
let mass = 1.0 // Particle mass (can be set to 1 in natural units)
let boxLength = 1.0 // Length of the 1D box
let numBasis = 100 // Number of basis states (particle in a box states)
let gridPoints = 1000 // Number of grid points for the potential

// Function to read the potential from a file
func loadPotential(from filePath: String) -> [Double] {
    do {
        let data = try String(contentsOfFile: filePath, encoding: .utf8)
        return data.split(separator: "
").compactMap { Double($0) }
    } catch {
        print("Error loading file: \(error)")
        return []
    }
}

// Particle in a box wavefunction
func particleInBoxWavefunction(n: Int, x: Double) -> Double {
    return sqrt(2 / boxLength) * sin(Double(n) * Double.pi * x / boxLength)
}

// Construct the Hamiltonian matrix
func constructHamiltonian(potential: [Double], dx: Double) -> [[Double]] {
    var H = Array(repeating: Array(repeating: 0.0, count: numBasis), count: numBasis)
    
    // Kinetic energy term (diagonal)
    for n in 1...numBasis {
        H[n - 1][n - 1] = (Double(n * n) * Double.pi * Double.pi * hbar * hbar) / (2 * mass * boxLength * boxLength)
    }
    
    // Potential energy term (off-diagonal elements)
    for i in 0..<numBasis {
        for j in 0..<numBasis {
            var integral = 0.0
            for k in 0..<potential.count {
                let x = Double(k) * dx
                let Vx = potential[k]
                integral += particleInBoxWavefunction(n: i + 1, x: x) * particleInBoxWavefunction(n: j + 1, x: x) * Vx * dx
            }
            H[i][j] += integral
        }
    }
    
    return H
}

// Solve the Hamiltonian matrix using LAPACK for eigenvalues and eigenvectors
func solveHamiltonian(H: [[Double]]) -> (eigenvalues: [Double], eigenvectors: [[Double]]) {
    var flatH = H.flatMap { $0 } // Flatten the matrix into a 1D array for LAPACK
    
    var N = __CLPK_integer(numBasis)
    var lwork = __CLPK_integer(-1)
    var work = [Double](repeating: 0.0, count: 1)
    var info: __CLPK_integer = 0
    var w = [Double](repeating: 0.0, count: numBasis)
    
    // Query for optimal workspace
    dgeev_("V", "V", &N, &flatH, &N, &w, nil, &N, nil, &N, &work, &lwork, &info)
    
    lwork = __CLPK_integer(work[0])
    work = [Double](repeating: 0.0, count: Int(lwork))
    
    var vl = [Double](repeating: 0.0, count: numBasis * numBasis)
    var vr = [Double](repeating: 0.0, count: numBasis * numBasis)
    
    // Solve for eigenvalues and eigenvectors
    dgeev_("V", "V", &N, &flatH, &N, &w, &vl, &N, &vr, &N, &work, &lwork, &info)
    
    // Extract eigenvectors
    var eigenvectors = [[Double]]()
    for i in 0..<numBasis {
        var eigenvector = [Double]()
        for j in 0..<numBasis {
            eigenvector.append(vr[i * numBasis + j])
        }
        eigenvectors.append(eigenvector)
    }
    
    return (w, eigenvectors)
}

// Main function to load the potential, construct the Hamiltonian, and solve the matrix equation
func main() {
    let filePath = "potential.txt" // Specify the path to the potential file
    let potential = loadPotential(from: filePath)
    
    if potential.isEmpty {
        print("No potential data loaded.")
        return
    }
    
    let dx = boxLength / Double(potential.count)
    
    // Construct the Hamiltonian matrix
    let H = constructHamiltonian(potential: potential, dx: dx)
    
    // Solve for eigenvalues and eigenvectors
    let (eigenvalues, eigenvectors) = solveHamiltonian(H: H)
    
    // Print the first few eigenvalues (energy levels)
    print("Eigenvalues (Energy Levels):")
    for i in 0..<10 {
        print("E[\(i)] = \(eigenvalues[i])")
    }
    
    // Print the first eigenvector (ground state wavefunction)
    print("
Ground state wavefunction:")
    for i in 0..<eigenvectors[0].count {
        print("ψ[\(i)] = \(eigenvectors[0][i])")
    }
}

main()
