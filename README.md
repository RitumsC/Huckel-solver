# Huckel-solver
A general Huckel solver as a part of programming practical option for part II chemistry tripos.

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
