"""General Huckel solver

Usage:    Edit the mol assignment for the desired adjaceny matrix that
            describes the molecule. For the given adjacency matrix the program
                calculates the Huckel energies and their degeneracies.
                    Results are displayed in a table.
"""

import numpy as np


def get_evals(matrix):
    """
    Get eigenvalues of Huckel matrix from the adjacency matrix

    """

    hMatrix = np.multiply(-1, matrix)
    evals, evecs = np.linalg.eig(hMatrix)
    evals.sort()
    return evals


def degen(evals):
    """
    Returns number of degeneracies for each eigenvalue

    """
    eps = 10e-5  # Threshold for energy differences to be considered degenerate
    deg = 1
    result = []
    for i in range(evals.size - 1):
        if (evals[i+1] - eps) < evals[i] < (evals[i+1] + eps):
            deg += 1
        else:
            result.append((evals[i], deg))
            deg = 1
    result.append((evals[-1], deg))
    return result

# Runs the program for now


mol = np.matrix(
        '0 1 0 0 0 0 \
        ;1 0 1 0 0 0 \
        ;0 1 0 1 0 0 \
        ;0 0 1 0 1 0 \
        ;0 0 0 1 0 1 \
        ;0 0 0 0 1 0 ')  # Initialzie test allyl molecule
evals = get_evals(mol)
result = degen(evals)

print(f'{"Energy":10}', f'{"Degeneracy":10}')
for a, b in result:
    print(f"{a:<10.3f}", b)












