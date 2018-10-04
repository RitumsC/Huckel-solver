"""General Huckel solver

Usage:    Edit the mol assignment for the desired adjaceny matrix that
            describes the molecule. For the given adjacency matrix the program
                calculates the Huckel energies and their degeneracies.
                    Results are displayed in a table.
"""

import numpy as np
import sys


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
    eps = 10e-7  # Threshold for energy differences to be considered degenerate
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


def adj_linear(n):
    """
    Creates adjacency matrix for linear polyene consiting of n atoms

    """
    matrix = np.eye(n, n, 1) + np.eye(n, n, -1)

    return matrix


def adj_cyclic(n):
    """
    Creates adjacency matrix for cyclic polyene consiting of n atoms

    """
    matrix = adj_linear(n) + np.eye(n, n, n-1) + np.eye(n, n, -n+1)

    return matrix


def return_result(matrix):
    """
    Prints the result given the adjacency matrix
    """
    evals = get_evals(matrix)
    result = degen(evals)
    print(f'{"Energy":10}', f'{"Degeneracy":10}')
    for a, b in result:
        print(f"A{-a:+.3f}*B", f'{b:^10}')


def main():
    arguments = sys.argv[1:]
    if arguments[0] == 'linear':
        try:
            print("The Huckel energies "
                  f"for linear polyene with {arguments[1]} carbons.\n")
            return_result(adj_linear(int(arguments[1])))
        except IndexError:
            print('You must supply the numer of atoms in the chain')
            print('e.g. : "python3 huckel.py linear 5"')
    elif arguments[0] == 'cyclic':
        try:
            print("The Huckel energies "
                  f"for cyclic polyene with {arguments[1]} carbons.\n")
            return_result(adj_cyclic(int(arguments[1])))
        except IndexError:
            print('You must supply the numer of atoms in the cycle')
            print('e.g. : "python3 huckel.py cyclic 5"')
    else:
        try:
            print("The Huckel energies "
                  f"for  {arguments[0]}.\n")
            return_result(np.loadtxt(arguments[0]+'.txt'))
        except OSError:
            print('No adjacency matrix supplied,\
                  check the argument filename')


if __name__ == '__main__':
    main()
