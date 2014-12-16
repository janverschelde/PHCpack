"""
Script to generate polynomial systems in the homotopy to solve
hypersurface Pieri problems in the numerical Schubert calculus.
The script needs sympy for the symbolic representation of the variable
entries in a matrix.  To compute determinants, numpy.linalg.det is used.
"""

import sympy as sp
from numpy import ndarray
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
    """
    coeff = complement(cffs, rows, cols, idx)
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
        result = equ[0]
        for i in range(1, len(equ)):
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
        result.append(pol)
    return result

def show_pieri_system():
    """
    Prompts the user for the number of rows and columns,
    the pivot information of the matrix of variables.
    Then shows the Pieri system for as many random matrices
    as there are variables in the matrix.
    """
    rows, cols, pvts = ask_inputs()
    eqs = pieri_system(rows, cols, pvts)
    answer = raw_input('see the polynomials ? (y/n) ')
    if(answer == 'y'):
        print 'the polynomials in the Pieri system :'
        for pol in eqs:
            print pol

def main():
    """
    Collects the test programs.
    """
    # show_variable_matrix()
    show_pieri_system()

if __name__ == '__main__':
    main()
