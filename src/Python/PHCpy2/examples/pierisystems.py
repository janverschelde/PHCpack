"""
Script to generate polynomial systems in the homotopy to solve
hypersurface Pieri problems in the numerical Schubert calculus.
The script needs sympy for the symbolic representation of the variable
entries in a matrix.  To compute determinants, numpy.linalg.det is used.
"""

import sympy as sp
from numpy import ndarray, zeros
from numpy.linalg import det, qr
from random import uniform
from math import sin, cos, pi

def random_matrix(rows, cols):
    """
    Returns a numpy array of shape rows and cols
    to store the data for a random complex matrix
    of as many rows as the value of rows and
    as many columns as the value of cols.
    The entries are random complex numbers
    of modulus one.  For numerical conditioning,
    the columns of the random matrix are orthonormalized.
    """
    numbers = ndarray(shape=(rows, cols), dtype=complex)
    for i in range(rows):
        for j in range(cols):
            angle = uniform(0, 2*pi)
            numbers[i, j] = complex(cos(angle), sin(angle))
    result, upper = qr(numbers)
    return result

def double2linear(nbrows, nbcols, row, col):
    """
    Returns one single index to locate the variable
    at position (row, col) in a variable matrix with
    as many rows as the value of nbrows and
    as many columns as the value of nbcols.
    The counting of indices starts at one
    and the matrix has a banded shape, with the ones
    in the identity matrix starting on top.
    The first variables start at the last column.
    """
    nvc = nbrows - nbcols    # number of variables per column
    prv = (nbcols - col)*nvc # number of preceding variables
    idx = row - col          # row index w.r.t. identity
    return prv + idx

def variables(rows, cols, pivots, linear=False):
    """
    Returns a symbolic matrix with as many rows
    as the value of rows and as many columns as
    the value of cols.  The pivots determine the
    position of the last nonzero row.
    The length of the list pivots must equal the
    value of the parameter cols.
    If linear, then the names of the variables
    will following a linear indexing scheme as
    defined by the function double2linear.
    """
    result = sp.Matrix(rows, cols, lambda i, j: 0)
    for i in range(cols):
        result[i, i] = 1
    for k in range(cols):
        for i in range(k+1, pivots[k]+1):
            if linear:
                idx = double2linear(rows, cols, i+1, k+1)
                name = 'x' + str(idx)
            else:
                name = 'x' + '_' + str(i+1) + '_' + str(k+1)
            sp.var(name)
            result[i, k] = eval(name)
    return result

def enumerate_all_minors(mat, rows, cols, crnt, accu):
    """
    Prints the row indices that define all minors
    of the matrix mat of dimension equal to cols.
    The matrix must have as many rows as the value of rows
    and as many columns as the value of cols.
    The variable crnt is the current column and accu
    accumulates the current selection of rows.
    """
    if(crnt >= cols):
        print accu
    else:
        start = -1
        if(len(accu) > 0):
            start = accu[-1]
        for i in range(start+1, rows):
            if(mat[i, crnt] != 0):
                enumerate_all_minors(mat, rows, cols, crnt+1, accu+[i])

def complement(cff, rows, cols, idx):
    """
    Returns the determinant of the complement of the selected
    rows as defined by the list idx of the matrix cff.
    """
    dim = rows - cols
    mat = ndarray(shape=(dim, dim), dtype=complex)
    cnt = -1
    for row in range(rows):
        if(not (row in idx)):
            cnt = cnt + 1
            for col in range(dim):
                mat[cnt, col] = cff[row, col]
    return det(mat)

def make_term(mat, rows, cols, idx, cffs):
    """
    Returns the string representation of the term defined by
    the selection of rows in the list idx.
    If the coefficient is zero, then the empty string is returned.
    """
    coeff = complement(cffs, rows, cols, idx)
    if(coeff == 0.0):
        result = ''
    else:
        cffre = '%.16e' % coeff.real
        cffim = '%+.16e' % coeff.imag
        result = '(' + cffre + cffim + '*i)'
        for k in range(cols):
            if(mat[idx[k], k] != 1):
                if(result == ''):
                    result = str(mat[idx[k],k])
                else:
                    result = result + '*' + str(mat[idx[k],k])
    return result

def enumerate_all_terms(mat, rows, cols, crnt, accu, cffs, equ):
    """
    Computes all terms in the Pieri condition imposed on the matrix mat
    with coefficients in the matrix cffs.
    The matrix mat has as many rows as the value of rows
    and as many columns as the value of cols.
    The matrix cffs has as many rows as the value of rows
    and as many columns as the value of rows - cols.
    The function assumes the matrix to have a banded shape.
    The variable crnt is the current column and accu
    accumulates the current selection of rows.
    The terms are accumulated as strings in the list equ.
    """
    if(crnt >= cols):
        term = make_term(mat, rows, cols, accu, cffs)
        equ.append(term)
    else:
        start = -1
        if(len(accu) > 0):
            start = accu[-1]
        for i in range(start+1, rows):
            if(mat[i, crnt] != 0):
                enumerate_all_terms(mat, rows, cols, crnt+1, accu+[i], \
                                    cffs, equ)

def all_minors(mat, rows, cols):
    """
    Prints the row indices that define all minors
    of the matrix mat of dimension equal to cols.
    The matrix must have as many rows as the value of rows
    and as many columns as the value of cols.
    """
    enumerate_all_minors(mat, rows, cols, 0, [])

def all_terms(mat, rows, cols, cffs):
    """
    Returns all terms in the Pieri condition imposed on the variable
    matrix mat by the coefficients in the matrix cffs.
    The matrix must have as many rows as the value of rows
    and as many columns as the value of cols.
    The numerical matrix cffs has as many rows as the value of rows
    and as many columns as the value of rows - cols.
    """
    result = []
    enumerate_all_terms(mat, rows, cols, 0, [], cffs, result) 
    return result

def sum_terms(equ):
    """
    Given in equ a list of string representing terms,
    returns the sum of the terms as a string.
    """
    result = ''
    if(len(equ) > 0):
        if(equ[0] != ''):
            result = equ[0]
        for i in range(1, len(equ)):
            if(equ[i] != ''):
                result = result + ' + ' + equ[i]
    return result

def number_of_variables(pivots):
    """
    Returns the number of variables in the variable matrix
    defined by the pivots.
    """
    result = 0
    for k in range(len(pivots)):
        result = result + pivots[k] - k;
    return result

def ask_inputs():
    """
    Prompts the user for the number of rows, 
    the number of columns, and pivots.
    These inputs are returned as a tuple.
    """
    rawnbr = raw_input("Give the number of rows : ")
    rawnbc = raw_input("Give the number of columns : ")
    rows = int(rawnbr)
    cols = int(rawnbc)
    pvts = [k for k in range(cols)]
    # pvts[cols-1] = rows-1 # rows - 1 is largest linear system
    pvts[cols-1] = rows-2   # rows - 2 allows for special plane
    print '-> the start pivots :', pvts
    print '-> number of variables :', number_of_variables(pvts)
    answer = raw_input('Do you want to give other pivots ? (y/n) ')
    if(answer == 'y'):
        rawpvts = raw_input('Give a list of increasing numbers : ')
        pvts = eval(rawpvts)
        print '-> your pivots :', pvts
        print '-> number of variables :', number_of_variables(pvts)
    return rows, cols, pvts

def show_variable_matrix():
    """
    Prompts the user for the number of rows and columns,
    the pivot information of the matrix of variables.
    Then shows the variable matrix and all its minors,
    along with the equation for a random coefficient matrix.
    """
    rows, cols, pvts = ask_inputs()
    answer = raw_input('Double indexing of variables ? (y/n) ')
    linear = (answer != 'y')
    varmat = variables(rows, cols, pvts, linear)
    print '-> the matrix of variables :'
    print varmat
    answer = raw_input('Enumerate all minors ? (y/n) ')
    if(answer == 'y'):
        all_minors(varmat, rows, cols)
    answer = raw_input('Enumerate all terms ? (y/n) ')
    dim = rows - cols
    cff = random_matrix(rows, dim)
    equ = all_terms(varmat, rows, cols, cff)
    print '-> all terms in the Pieri condition :'
    print sum_terms(equ)

def rightmost_nonfull_pivot(dim, pivots):
    """
    Let dim be the dimension of the ambient space.
    The list of pivots is full if equal to 
    [dim - L, dim - L + 1, .., dim - 1], L = len(pivots).
    If the list of pivots is full, then -1 is returned,
    otherwise return the index of the rightmost entry
    in pivots that can still be increased.
    """
    L = len(pivots) 
    for k in range(L-1,-1,-1):
        if(pivots[k] < dim - L + k):
            return k
    return -1

def update_pivots(dim, pivots):
    """
    Returns the updated pivots after incrementing
    the rightmost nonfull pivot in pivots.
    If the pivots are full, then pivots are returned.
    """
    piv = rightmost_nonfull_pivot(dim, pivots)
    if(piv < 0):
        result = pivots
    else:
        result = [x for x in pivots]
        result[piv] = result[piv] + 1
    return result

def extend_solutions(sols, dim, pivots):
    """
    Adds one extra zero coordinate to all solutions in sols,
    corresponding to the rightmost nonfull pivot.
    """
    from phcpy.solutions import strsol2dict, variables, make_solution
    piv = rightmost_nonfull_pivot(dim, pivots)
    print '-> in extend_solutions, piv =', piv
    if(piv < 0):
        result = sols
    else:
        rown = pivots[piv] + 1  # row number of new variable 
        name = 'x_' + str(rown+1) + '_' + str(piv+1)
        print '-> the name of the new variable :', name
        result = []
        for sol in sols:
            dicsol = strsol2dict(sol)
            solvar = variables(dicsol)
            solval = [dicsol[var] for var in solvar] 
            solvar.append(name)
            solval.append(complex(0.0))
            result.append(make_solution(solvar, solval))
    return result

def special_plane(dim, pivots):
    """
    Returns a matrix with as many rows as the value of dim
    and as many columns as the value of dim - len(pivots).
    The plane is a special plane so the Pieri condition is
    satisfied when the variable with the rightmost nonfull
    pivot is set to zero.
    """
    cols = dim - len(pivots)
    result = zeros(shape=(dim, cols), dtype=complex)
    ind = 0
    for k in range(dim):
        if(not (k in pivots)):
            (result[k, ind], ind) = (1.0, ind+1)
    return result

def pieri_system(rows, cols, pivots, planes, linear=False):
    """
    Generates a Pieri system for the numbers of rows in rows,
    the number of columns in cols, and the given pivots.
    The list planes must contain at least as many matrices
    as the number of variables in the pivots.
    The number of rows of the matrices in planes is the value
    of rows and the number of columns of the matrices in planes
    must equal to rows - cols.
    If linear, then variables in the matrix are named with 
    a strict linear index instead of a double index.
    """
    mat = variables(rows, cols, pivots, linear)
    nvr = number_of_variables(pivots)
    dim = rows - cols
    result = []
    for k in range(nvr):
        cff = planes[k]
        equ = all_terms(mat, rows, cols, cff)
        print '-> number of terms in equation %d : %d' % (k+1, len(equ))
        pol = sum_terms(equ)
        result.append(pol + ';')
    return result

def random_planes(rows, cols, pivots):
    """
    Generates as many random planes as the number of variables
    described by the list of pivots.
    The dimension of the planes is the value of cols
    in a space of dimension equal to the value of rows.
    """
    nvr = number_of_variables(pivots)
    dim = rows - cols
    planes = []
    for k in range(nvr):
        cff = random_matrix(rows, dim)
        planes.append(cff)
    return planes

def random_pieri_system(rows, cols, pivots, linear=False):
    """
    Generates a Pieri system for the numbers of rows in rows,
    the number of columns in cols, and the given pivots.
    The conditions are formulated for random planes.
    If linear, then variables in the matrix are named with 
    a strict linear index instead of a double index.
    """
    planes = random_planes(rows, cols, pivots)
    return pieri_system(rows, cols, pivots, planes, linear)

def start_pieri_system(rows, cols, pivots, planes, linear=False):
    """
    Generates a Pieri system where the last equation
    is for a special plane.  For all other equations,
    the planes in the list planes will be used.
    If linear, then variables in the matrix are named with 
    a strict linear index instead of a double index.
    """
    newpiv = update_pivots(rows, pivots)
    mat = variables(rows, cols, newpiv, linear)
    nvr = number_of_variables(pivots)
    dim = rows - cols
    result = []
    for k in range(nvr):
        cff = planes[k]
        equ = all_terms(mat, rows, cols, cff)
        print '-> number of terms in equation %d : %d' % (k+1, len(equ))
        pol = sum_terms(equ)
        result.append(pol + ';')
    lastcff = special_plane(rows, newpiv)
    print '-> the special plane :'
    print lastcff
    lastequ = all_terms(mat, rows, cols, lastcff)
    lastpol = sum_terms(lastequ)
    result.append(lastpol + ';')
    return result

def random_start_pieri_system(rows, cols, pivots, linear=False):
    """
    Generates a Pieri system where the last equation
    is for a special plane.  For all other equations,
    random matrices will be generated.
    If linear, then variables in the matrix are named with 
    a strict linear index instead of a double index.
    """
    planes = random_planes(rows, cols, pivots)
    return start_pieri_system(rows, cols, pivots, planes, linear) 

def solve_pieri_system(eqs):
    """
    Applies the black box solver of phcpy to solve the equations
    defined by the strings in the list eqs.  Note that, because
    of its particular structure, this will only work well if the 
    system is linear.
    """
    from phcpy.solver import solve
    sols = solve(eqs)
    print '-> the solutions :'
    for sol in sols:
        print sol
    return sols

def write_to_name_file(name, eqs, sols, dimspace, dimplane):
    """
    Opens a file with the given name for writing and the writes the
    system defined by the string representations of polynomials to file.
    If the list sols is not empty, then the solutions will be
    writen as well to file.  The parameters dimspace and dimplane
    respectively represent the dimensions of the space and 
    the solution planes.
    """
    file = open(name, 'w')
    file.write(str(len(eqs)) + '\n')
    for pol in eqs:
        file.write(pol + '\n')
    file.write('\n')
    file.write('TITLE : Pieri conditions on %d-planes in %d-space\n' \
               % (dimplane, dimspace))
    file.write('\n')
    if(len(sols) > 0):
        cnt = 0
        file.write('THE SOLUTIONS :\n')
        file.write(str(len(sols)) + ' ' + str(len(eqs)) + '\n')
        file.write('=====================================================\n')
        for sol in sols:
            cnt = cnt + 1
            file.write('solution ' + str(cnt) + ' : \n')
            file.write(sol + '\n')

def write_to_file(eqs, sols, dimspace, dimplane):
    """
    Prompts the user for a file name and then writes the system
    defined by the string representations of polynomials to file.
    If the list sols is not empty, then the solutions will be
    writen as well to file.  The parameters dimspace and dimplane
    respectively represent the dimensions of the space and 
    the solution planes.
    """
    name = raw_input('Give a name of a file to write to : ')
    write_to_name_file(name, eqs, sols, dimspace, dimplane)

def show_pieri_system():
    """
    Prompts the user for the number of rows and columns,
    the pivot information of the matrix of variables.
    Then shows the Pieri system for as many random matrices
    as there are variables in the matrix.
    """
    rows, cols, pvts = ask_inputs()
    answer = raw_input('Double indexing of variables ? (y/n) ')
    linear = (answer != 'y')
    answer = raw_input('Make start Pieri system ? (y/n) ')
    if(answer == 'y'):
        eqs = random_start_pieri_system(rows, cols, pvts, linear)
    else:
        eqs = random_pieri_system(rows, cols, pvts, linear)
    answer = raw_input('Do you want to see the polynomials ? (y/n) ')
    if(answer == 'y'):
        print '-> the polynomials in the Pieri system :'
        for pol in eqs:
            print pol
    answer = raw_input('Apply solve of phcpy ? (y/n) ')
    if(answer == 'y'):
        sols = solve_pieri_system(eqs)
    else:
        sols = []
    answer = raw_input('Write the system to file ? (y/n) ')
    if(answer == 'y'):
        write_to_file(eqs, sols, rows, cols)

def track_sequence_of_pieri_systems():
    """
    Makes a sequence of Pieri systems
    and runs the path trackers on one path.
    """
    from phcpy.trackers import track
    nbrows, nbcols, pvts = ask_inputs()
    answer = raw_input('Double indexing of variables ? (y/n) ')
    linear = (answer != 'y')
    dim = nbrows - nbcols
    planes = random_planes(nbrows, nbcols, pvts)
    start = start_pieri_system(nbrows, nbcols, pvts, planes, linear)
    startsols = solve_pieri_system(start)
    newpvts = update_pivots(nbrows, pvts)
    planes.append(random_matrix(nbrows, dim))
    target = pieri_system(nbrows, nbcols, newpvts, planes, linear)
    sols = track(target, start, startsols)
    print '-> the solutions after track :'
    for sol in sols:
        print sol
    while True:
        answer = raw_input('Continue to next level ? (y/n) ')
        if(answer != 'y'):
            break
        else:
            pvts = [piv for piv in newpvts]
            newpvts = update_pivots(nbrows, pvts)
            if(pvts == newpvts):
                print '-> no extension of pivots possible: no next level'
                break
            start = start_pieri_system(nbrows, nbcols, pvts, planes, linear)
            print '-> the new start system :'
            for pol in start:
                print pol
            startsols = extend_solutions(sols, nbrows, pvts)
            print '-> the extended solutions :'
            for sol in startsols:
                print sol
            planes.append(random_matrix(nbrows, dim))
            target = pieri_system(nbrows, nbcols, newpvts, planes, linear)
            print '-> the new target system :'
            for pol in target:
                print pol
            sols = track(target, start, startsols)
            print '-> the solutions after track :'
            for sol in sols:
                print sol
    answer = raw_input('Write system and solutions to file ? (y/n) ')
    if(answer == 'y'):
        write_to_file(target, sols, nbrows, nbcols)

def write_sequence_of_pieri_systems():
    """
    Makes a sequence of Pieri systems
    and write each pair of start and target system to file.
    """
    nbrows, nbcols, pvts = ask_inputs()
    answer = raw_input('Double indexing of variables ? (y/n) ')
    linear = (answer != 'y')
    answer = raw_input('Continue through all levels ? (y/n) ')
    contin = (answer == 'y')
    dim = nbrows - nbcols
    planes = random_planes(nbrows, nbcols, pvts)
    start = start_pieri_system(nbrows, nbcols, pvts, planes, linear)
    name = 'pieri' + str(nbrows) + str(nbcols) + 'start0'
    write_to_name_file(name, start, [], nbrows, nbcols)
    newpvts = update_pivots(nbrows, pvts)
    planes.append(random_matrix(nbrows, dim))
    target = pieri_system(nbrows, nbcols, newpvts, planes, linear)
    name = 'pieri' + str(nbrows) + str(nbcols) + 'target0'
    write_to_name_file(name, target, [], nbrows, nbcols)
    cnt = 0
    while True:
        if(not contin):
            answer = raw_input('Continue to next level ? (y/n) ')
            if(answer != 'y'):
                break
        else:
            cnt = cnt + 1
            pvts = [piv for piv in newpvts]
            newpvts = update_pivots(nbrows, pvts)
            if(pvts == newpvts):
                print '-> no extension of pivots possible: no next level'
                break
            start = start_pieri_system(nbrows, nbcols, pvts, planes, linear)
            name = 'pieri' + str(nbrows) + str(nbcols) + 'start' + str(cnt)
            write_to_name_file(name, start, [], nbrows, nbcols)
            planes.append(random_matrix(nbrows, dim))
            target = pieri_system(nbrows, nbcols, newpvts, planes, linear)
            name = 'pieri' + str(nbrows) + str(nbcols) + 'target' + str(cnt)
            write_to_name_file(name, target, [], nbrows, nbcols)

def main():
    """
    Collects the test programs.
    """
    print 'MENU to test setup of Pieri systems :'
    print '  1. show a matrix of variables'
    print '  2. make a Pieri system for pivots'
    print '  3. run path trackers between Pieri systems'
    print '  4. generate sequence of Pieri systems to file'
    answer = raw_input('Type 1, 2, 3, or 4 to select a test : ')
    if(answer == '1'):
        show_variable_matrix()
    elif(answer == '2'):
        show_pieri_system()
    elif(answer == '3'):
        track_sequence_of_pieri_systems()
    elif(answer == '4'):
        write_sequence_of_pieri_systems()
    else:
        print 'Wrong answer, please try again ...'

if __name__ == '__main__':
    main()
