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

def variables(rows, cols, pivots):
    """
    Returns a symbolic matrix with as many rows
    as the value of rows and as many columns as
    the value of cols.  The pivots determine the
    position of the last nonzero row.
    The length of the list pivots must equal the
    value of the parameter cols.
    """
    result = sp.Matrix(rows, cols, lambda i, j: 0)
    for i in range(cols):
        result[i, i] = 1
    for k in range(cols):
        for i in range(k+1, pivots[k]+1):
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
    pvts[cols-1] = rows-1
    print 'the start pivots :', pvts
    print 'number of variables :', number_of_variables(pvts)
    answer = raw_input('Do you want to give other pivots ? (y/n) ')
    if(answer == 'y'):
        rawpvts = raw_input('Give a list of increasing numbers : ')
        pvts = eval(rawpvts)
        print 'your pivots :', pvts
        print 'number of variables :', number_of_variables(pvts)
    return rows, cols, pvts

def show_variable_matrix():
    """
    Prompts the user for the number of rows and columns,
    the pivot information of the matrix of variables.
    Then shows the variable matrix and all its minors,
    along with the equation for a random coefficient matrix.
    """
    rows, cols, pvts = ask_inputs()
    varmat = variables(rows, cols, pvts)
    print 'the matrix of variables :'
    print varmat
    answer = raw_input('Enumerate all minors ? (y/n) ')
    if(answer == 'y'):
        all_minors(varmat, rows, cols)
    answer = raw_input('Enumerate all terms ? (y/n) ')
    dim = rows - cols
    cff = random_matrix(rows, dim)
    equ = all_terms(varmat, rows, cols, cff)
    print 'all terms in the Pieri condition :'
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

def pieri_system(rows, cols, pivots):
    """
    Generates a Pieri system for the numbers of rows in rows,
    the number of columns in cols, and the given pivots.
    """
    mat = variables(rows, cols, pivots)
    nvr = number_of_variables(pivots)
    dim = rows - cols
    result = []
    for k in range(nvr):
        cff = random_matrix(rows, dim)
        equ = all_terms(mat, rows, cols, cff)
        print 'number of terms in equation %d : %d' % (k+1, len(equ))
        pol = sum_terms(equ)
        result.append(pol + ';')
    return result

def start_pieri_system(rows, cols, pivots):
    """
    Generates a Pieri system where the last equation
    is for a special plane.
    """
    newpiv = update_pivots(rows, pivots)
    mat = variables(rows, cols, newpiv)
    nvr = number_of_variables(pivots)
    dim = rows - cols
    result = []
    for k in range(nvr):
        cff = random_matrix(rows, dim)
        equ = all_terms(mat, rows, cols, cff)
        print 'number of terms in equation %d : %d' % (k+1, len(equ))
        pol = sum_terms(equ)
        result.append(pol + ';')
    lastcff = special_plane(rows, newpiv)
    print 'the special plane :'
    print lastcff
    lastequ = all_terms(mat, rows, cols, lastcff)
    lastpol = sum_terms(lastequ)
    result.append(lastpol + ';')
    return result

def solve_pieri_system(eqs):
    """
    Applies the black box solver of phcpy to solve the equations
    defined by the strings in the list eqs.  Note that, because
    of its particular structure, this will only work well if the 
    system is linear.
    """
    from phcpy.solver import solve
    sols = solve(eqs)
    print 'the solutions :'
    for sol in sols:
        print sol
    return sols

def show_pieri_system():
    """
    Prompts the user for the number of rows and columns,
    the pivot information of the matrix of variables.
    Then shows the Pieri system for as many random matrices
    as there are variables in the matrix.
    """
    rows, cols, pvts = ask_inputs()
    answer = raw_input('start Pieri system ? (y/n) ')
    if(answer == 'y'):
        eqs = start_pieri_system(rows, cols, pvts)
    else:
        eqs = pieri_system(rows, cols, pvts)
    answer = raw_input('see the polynomials ? (y/n) ')
    if(answer == 'y'):
        print 'the polynomials in the Pieri system :'
        for pol in eqs:
            print pol
    answer = raw_input('apply solve of phcpy ? (y/n) ')
    if(answer == 'y'):
        sols = solve_pieri_system(eqs)

def main():
    """
    Collects the test programs.
    """
    # show_variable_matrix()
    show_pieri_system()

if __name__ == '__main__':
    main()
