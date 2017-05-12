import numpy as np

from fdispersion import fgetc6, fgetc6_periodic

def getc6(coords, natoms, atnum):

    fcoords = np.array(coords, order="F").T
    fatnum = np.array(atnum, dtype=np.int32)

    return fgetc6(fcoords, natoms, fatnum)

def getc6_pbc(coords, natoms, atnum, latvecs):

    fcoords = np.array(coords, order="F").T
    flatvecs = np.array(latvecs, order="F").T
    fatnum = np.array(atnum, dtype=np.int32)

    return fgetc6_periodic(fcoords, natoms, fatnum, flatvecs)


if __name__ == "__main__":

    from dna import coords as dnacoords
    from dna import natoms as dnanatoms
    from dna import atnum as dnaatnum

    c6 = getc6(dnacoords, dnanatoms, dnaatnum)

    print "C6"
    print c6

    from dna import coords as dnacoords
    from dna import natoms as dnanatoms
    from dna import atnum as dnaatnum
    from dna import latvecs as dnalatvecs

    c6p = getc6_pbc(dnacoords, dnanatoms, dnaatnum, dnalatvecs)

    print "C6 periodic"
    print c6[np.diag_indices_from(c6)]
    print c6p[np.diag_indices_from(c6p)]
