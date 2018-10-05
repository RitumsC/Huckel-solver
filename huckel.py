"""General Huckel solver

Run in command line:
    $python3 huckel.py mode n

Arguments:
    mode
        Decides which task to do, posibilities:
           linear
           chain
           name of text file that contains the adjacency matrix of atoms
    n
        Number of atoms to be considered, only used for linear/chain mode.

For Platonic solids and Buckinsterfullerene the adjaceny matrix are prepared
and kept in the same git repository.

Example command lines:
     $python3 huckel.py linear 5
     $python3 huckel.py cylic 4
     $python3 huckel.py cube
"""

import numpy as np
import sys


def get_evals(adjMatrix):
    """Get eigenvalues of Huckel matrix from the adjacency matrix

    hMatrix = adjMatrix * (-1)
    """
    evals, evecs = np.linalg.eig(adjMatrix)
    evals.sort()

    return (-1)*evals


def degen(evals):
    """Returns list of pairs listing degeneracies for each eigenvalue

        Must take in sorted evals
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
    """Creates adjacency matrix for linear polyene consiting of n atoms

    """
    matrix = np.eye(n, n, 1) + np.eye(n, n, -1)

    return matrix


def adj_cyclic(n):
    """Creates adjacency matrix for cyclic polyene consiting of n atoms

    """
    matrix = adj_linear(n) + np.eye(n, n, n-1) + np.eye(n, n, -n+1)

    return matrix


def return_result(adjMatrix):
    """Prints the result given the adjacency matrix

    """
    evals = get_evals(adjMatrix)
    degEvals = degen(evals)
    print(f'{"Energy":10}', f'{"Degeneracy"}')
    for a, b in degEvals:
        print(f"A{-a.real:+.3f}*B", f'{b:^10}')


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
            print("No adjacency matrix supplied, "
                  "check typos or create new adjacency matrix")


if __name__ == '__main__':
    main()
