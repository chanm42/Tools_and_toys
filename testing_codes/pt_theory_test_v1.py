import numpy as np

# Small example of first-order non-degenerate matrix perturbation theory
# for a 5x5 symmetric matrix.

def first_order_perturbation(A0, E, idx):
    # Compute first-order eigenvalue correction for eigenpair idx
    evals, evecs = np.linalg.eigh(A0)
    v = evecs[:, idx]
    ans = v.T @ E @ v
    return ans


def main():
    np.random.seed(42)

    # Unperturbed base matrix A0 with simple, well-separated eigenvalues.
    A0 = np.diag([1.0, 2.0, 3.0, 4.0, 5.0])

    # Symmetric perturbation matrix E with small entries.
    E = np.array([
        [0.00, 0.02, 0.01, 0.00, 0.00],
        [0.02, 0.00, 0.03, 0.01, 0.00],
        [0.01, 0.03, 0.00, 0.02, 0.01],
        [0.00, 0.01, 0.02, 0.00, 0.03],
        [0.00, 0.00, 0.01, 0.03, 0.00],
    ])

    epsilon = 0.1
    A = A0 + epsilon * E

    evals0, evecs0 = np.linalg.eigh(A0)
    evals, evecs = np.linalg.eigh(A)

    print("A0 (base matrix):")
    print(A0)
    print("\nE (perturbation matrix):")
    print(E)
    print(f"\nPerturbation strength epsilon = {epsilon}\n")

    print("Unperturbed eigenvalues:")
    print(evals0)
    print("\nPerturbed eigenvalues:")
    print(evals)

    print("\nFirst-order eigenvalue corrections:")
    corrections = [first_order_perturbation(A0, E, i) for i in range(5)]
    print(np.array(corrections))

    print("\nPerturbation theory predicts new eigenvalues:")
    predicted = evals0 + epsilon * np.array(corrections)
    print(predicted)

    print("\nDifference between exact perturbed eigenvalues and first-order prediction:")
    print(evals - predicted)

main()
