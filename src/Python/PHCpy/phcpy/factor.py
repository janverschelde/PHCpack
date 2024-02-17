"""
Given a witness set representation of a pure dimensional solution set,
the functions in this module separate the generic points in the witness set
according to the irreducible components of the solution set.
"""
from ctypes import c_int32, c_double, pointer
from random import uniform
from cmath import exp, pi
from phcpy.version import get_phcfun, int4a2nbr
from phcpy.solutions import make_solution, endmultiplicity
from phcpy.solutions import get_double_solutions
from phcpy.solutions import set_double_solutions
from phcpy.solutions import clear_double_solutions
from phcpy.solutions import get_double_double_solutions
from phcpy.solutions import set_double_double_solutions
from phcpy.solutions import clear_double_double_solutions
from phcpy.solutions import get_quad_double_solutions
from phcpy.solutions import set_quad_double_solutions
from phcpy.solutions import clear_quad_double_solutions
from phcpy.solutions import filter_zero_coordinates
from phcpy.families import cyclic
from phcpy.solver import solve
from phcpy.sets import double_embed, double_double_embed
from phcpy.sets import quad_double_embed
from phcpy.sets import set_double_witness_set
from phcpy.sets import set_double_double_witness_set
from phcpy.sets import set_quad_double_witness_set
from phcpy.sets import set_double_laurent_witness_set
from phcpy.sets import set_double_double_laurent_witness_set
from phcpy.sets import set_quad_double_laurent_witness_set

def set_double_verbose(verbose=True, vrblvl=0):
    """
    Sets the state of the monodromy algorithm in double precision
    to verbose if vrblvl > 0, otherwise it is set to remain mute.
    """
    if vrblvl > 0:
        print('in set_double_verbose, verbose :', verbose)
    phc = get_phcfun(vrblvl-1)
    aaa = pointer(c_int32(0))
    bbb = pointer(c_int32(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl-1)
    if vrblvl > 0:
        print('-> set_double_monodromy_verbose calls phc', end='')
    if verbose:
        retval = phc(630, aaa, bbb, ccc, vrb)
    else:
        retval = phc(39, aaa, bbb, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    return retval

def set_double_double_verbose(verbose=True, vrblvl=0):
    """
    Sets the state of the monodromy algorithm in double double precision
    to verbose if vrblvl > 0, otherwise it is set to remain mute.
    """
    if vrblvl > 0:
        print('in set_double_double_verbose, verbose :', verbose)
    phc = get_phcfun(vrblvl-1)
    aaa = pointer(c_int32(0))
    bbb = pointer(c_int32(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl-1)
    if vrblvl > 0:
        print('-> set_double_double_verbose calls phc', end='')
    if verbose:
        retval = phc(660, aaa, bbb, ccc, vrb)
    else:
        retval = phc(658, aaa, bbb, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    return retval

def set_quad_double_verbose(verbose=True, vrblvl=0):
    """
    Sets the state of the monodromy algorithm in quad double precision
    to verbose if vrblvl > 0, otherwise it is set to remain mute.
    """
    if vrblvl > 0:
        print('in set_quad_double_verbose, verbose :', verbose)
    phc = get_phcfun(vrblvl-1)
    aaa = pointer(c_int32(0))
    bbb = pointer(c_int32(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl-1)
    if vrblvl > 0:
        print('-> set_quad_double_verbose calls phc', end='')
    if verbose:
        retval = phc(690, aaa, bbb, ccc, vrb)
    else:
        retval = phc(688, aaa, bbb, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    return retval

def double_assign_labels(nvr, vrblvl=0):
    """
    Assigns a unique label to the multiplicity field
    for each solution set in double precision.
    """
    if vrblvl > 0:
        print('in double_assign_labels, nvr :', nvr)
    sols = get_double_solutions(vrblvl-1)
    newsols = []
    for (idx, sol) in enumerate(sols):
        data = sol.split('m : 1')
        idxm = f"m : {idx+1}"
        newsol = idxm.join(data)
        newsols.append(newsol)
    clear_double_solutions(vrblvl-1)
    set_double_solutions(nvr, newsols, vrblvl-1)
    return 0

def double_double_assign_labels(nvr, vrblvl=0):
    """
    Assigns a unique label to the multiplicity field
    for each solution set in double double precision.
    """
    if vrblvl > 0:
        print('in double_double_assign_labels, nvr :', nvr)
    sols = get_double_double_solutions(vrblvl-1)
    newsols = []
    for (idx, sol) in enumerate(sols):
        data = sol.split('m : 1')
        idxm = f"m : {idx+1}"
        newsol = idxm.join(data)
        newsols.append(newsol)
    clear_double_double_solutions(vrblvl-1)
    set_double_double_solutions(nvr, newsols, vrblvl-1)
    return 0

def quad_double_assign_labels(nvr, vrblvl=0):
    """
    Assigns a unique label to the multiplicity field
    for each solution set in quad double precision.
    """
    if vrblvl > 0:
        print('in quad_double_assign_labels, nvr :', nvr)
    sols = get_quad_double_solutions(vrblvl-1)
    newsols = []
    for (idx, sol) in enumerate(sols):
        data = sol.split('m : 1')
        idxm = f"m : {idx+1}"
        newsol = idxm.join(data)
        newsols.append(newsol)
    clear_quad_double_solutions(vrblvl-1)
    set_quad_double_solutions(nvr, newsols, vrblvl-1)
    return 0

def initialize_double_sampler(dim, vrblvl=0):
    """
    Initializes the sampling machine with a witness set,
    set in double precision.
    """
    if vrblvl > 0:
        print('in initialize_double_sampler, dim :', dim)
    phc = get_phcfun(vrblvl-1)
    adim = pointer(c_int32(dim))
    bbb = pointer(c_int32(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl-1)
    if vrblvl > 0:
        print('-> initialize_double_sampler calls phc', end='')
    retval = phc(42, adim, bbb, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    return retval

def initialize_double_double_sampler(dim, vrblvl=0):
    """
    Initializes the sampling machine with a witness set,
    set in double double precision.
    """
    if vrblvl > 0:
        print('in initialize_double_double_sampler, dim :', dim)
    phc = get_phcfun(vrblvl-1)
    adim = pointer(c_int32(dim))
    bbb = pointer(c_int32(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl-1)
    if vrblvl > 0:
        print('-> initialize_double_double_sampler calls phc', end='')
    retval = phc(632, adim, bbb, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    return retval

def initialize_quad_double_sampler(dim, vrblvl=0):
    """
    Initializes the sampling machine with a witness set,
    set in quad double precision.
    """
    if vrblvl > 0:
        print('in initialize_quad_double_sampler, dim :', dim)
    phc = get_phcfun(vrblvl-1)
    adim = pointer(c_int32(dim))
    bbb = pointer(c_int32(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl-1)
    if vrblvl > 0:
        print('-> initialize_quad_double_sampler calls phc', end='')
    retval = phc(662, adim, bbb, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    return retval

def initialize_double_laurent_sampler(dim, vrblvl=0):
    """
    Initializes the sampling machine with a witness set,
    for a laurent system set in double precision.
    """
    if vrblvl > 0:
        print('in initialize_double_laurent_sampler, dim :', dim)
    phc = get_phcfun(vrblvl-1)
    adim = pointer(c_int32(dim))
    bbb = pointer(c_int32(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl-1)
    if vrblvl > 0:
        print('-> initialize_double_laurent_sampler calls phc', end='')
    retval = phc(804, adim, bbb, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    return retval

def initialize_double_double_laurent_sampler(dim, vrblvl=0):
    """
    Initializes the sampling machine with a witness set,
    for a laurent system set in double double precision.
    """
    if vrblvl > 0:
        print('in initialize_double_double_laurent_sampler, dim :', dim)
    phc = get_phcfun(vrblvl-1)
    adim = pointer(c_int32(dim))
    bbb = pointer(c_int32(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl-1)
    if vrblvl > 0:
        print('-> initialize_double_double_laurent_sampler calls phc', end='')
    retval = phc(805, adim, bbb, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    return retval

def initialize_quad_double_laurent_sampler(dim, vrblvl=0):
    """
    Initializes the sampling machine with a witness set,
    for a laurent system set in quad double precision.
    """
    if vrblvl > 0:
        print('in initialize_quad_double_laurent_sampler, dim :', dim)
    phc = get_phcfun(vrblvl-1)
    adim = pointer(c_int32(dim))
    bbb = pointer(c_int32(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl-1)
    if vrblvl > 0:
        print('-> initialize_quad_double_laurent_sampler calls phc', end='')
    retval = phc(806, adim, bbb, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    return retval

def initialize_double_monodromy(nbloops, deg, dim, vrblvl=0):
    """
    Initializes to run nbloops monodromy loops to factor
    a positive dimension solution set of dimension dim and of degree deg
    in double precision.
    """
    if vrblvl > 0:
        print('in initialize_double_monodromy', end='')
        print(', nbloops :', nbloops, ', deg :', deg, ', dim :', dim)
    phc = get_phcfun(vrblvl-1)
    anbloops = pointer(c_int32(nbloops))
    bpars = int4a2nbr([deg, dim], vrblvl=vrblvl-1)
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl-1)
    if vrblvl > 0:
        print('-> initialize_double_monodromy calls phc', end='')
    retval = phc(50, anbloops, bpars, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    return retval

def initialize_double_double_monodromy(nbloops, deg, dim, vrblvl=0):
    """
    Initializes to run nbloops monodromy loops to factor
    a positive dimension solution set of dimension dim and of degree deg
    in double double precision.
    """
    if vrblvl > 0:
        print('in initialize_double_double_monodromy', end='')
        print(', nbloops :', nbloops, ', deg :', deg, ', dim :', dim)
    phc = get_phcfun(vrblvl-1)
    anbloops = pointer(c_int32(nbloops))
    bpars = int4a2nbr([deg, dim], vrblvl=vrblvl-1)
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl-1)
    if vrblvl > 0:
        print('-> initialize_double_double_monodromy calls phc', end='')
    retval = phc(640, anbloops, bpars, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    return retval

def initialize_quad_double_monodromy(nbloops, deg, dim, vrblvl=0):
    """
    Initializes to run nbloops monodromy loops to factor
    a positive dimension solution set of dimension dim and of degree deg
    in double double precision.
    """
    if vrblvl > 0:
        print('in initialize_quad_double_monodromy', end='')
        print(', nbloops :', nbloops, ', deg :', deg, ', dim :', dim)
    phc = get_phcfun(vrblvl-1)
    anbloops = pointer(c_int32(nbloops))
    bpars = int4a2nbr([deg, dim], vrblvl=vrblvl-1)
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl-1)
    if vrblvl > 0:
        print('-> initialize_double_double_monodromy calls phc', end='')
    retval = phc(670, anbloops, bpars, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    return retval

def preset_double_solutions(vrblvl=0):
    """
    Initializes the start solutions for monodromy loops to the
    solutions set in double precision.
    """
    if vrblvl > 0:
        print('in preset_double_solutions ...')
    phc = get_phcfun(vrblvl-1)
    aaa = pointer(c_int32(0))
    bbb = pointer(c_int32(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl-1)
    if vrblvl > 0:
        print('-> preset_double_solutions calls phc', end='')
    retval = phc(51, aaa, bbb, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    return retval

def preset_double_double_solutions(vrblvl=0):
    """
    Initializes the start solutions for monodromy loops to the
    solutions set in double double precision.
    """
    if vrblvl > 0:
        print('in preset_double_double_solutions ...')
    phc = get_phcfun(vrblvl-1)
    aaa = pointer(c_int32(0))
    bbb = pointer(c_int32(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl-1)
    if vrblvl > 0:
        print('-> preset_double_double_solutions calls phc', end='')
    retval = phc(641, aaa, bbb, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    return retval

def preset_quad_double_solutions(vrblvl=0):
    """
    Initializes the start solutions for monodromy loops to the
    solutions set in quad double precision.
    """
    if vrblvl > 0:
        print('in preset_quad_double_solutions ...')
    phc = get_phcfun(vrblvl-1)
    aaa = pointer(c_int32(0))
    bbb = pointer(c_int32(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl-1)
    if vrblvl > 0:
        print('-> preset_quad_double_solutions calls phc', end='')
    retval = phc(671, aaa, bbb, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    return retval

def reset_double_solutions(vrblvl=0):
    """
    Resets the start solutions for monodromy loops in double precision
    to the solutions used to initialize the start solutions with.
    """
    if vrblvl > 0:
        print('in reset_double_solutions ...')
    phc = get_phcfun(vrblvl-1)
    aaa = pointer(c_int32(0))
    bbb = pointer(c_int32(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl-1)
    if vrblvl > 0:
        print('-> reset_double_solutions calls phc', end='')
    retval = phc(48, aaa, bbb, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    return retval

def reset_double_double_solutions(vrblvl=0):
    """
    Resets the start solutions for monodromy loops in double double precision
    to the solutions used to initialize the start solutions with.
    """
    if vrblvl > 0:
        print('in reset_double_double_solutions ...')
    phc = get_phcfun(vrblvl-1)
    aaa = pointer(c_int32(0))
    bbb = pointer(c_int32(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl-1)
    if vrblvl > 0:
        print('-> reset_double_double_solutions calls phc', end='')
    retval = phc(638, aaa, bbb, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    return retval

def reset_quad_double_solutions(vrblvl=0):
    """
    Resets the start solutions for monodromy loops in quad double precision
    to the solutions used to initialize the start solutions with.
    """
    if vrblvl > 0:
        print('in reset_quad_double_solutions ...')
    phc = get_phcfun(vrblvl-1)
    aaa = pointer(c_int32(0))
    bbb = pointer(c_int32(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl-1)
    if vrblvl > 0:
        print('-> reset_quad_double_solutions calls phc', end='')
    retval = phc(668, aaa, bbb, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    return retval

def set_double_slice(equ, idx, cff, vrblvl=0):
    """
    Sets the coefficients of slicing equation with index equ
    at position idx to the value in double precision:
    cff[0] + cff[1]*complex(0, 1).
    """
    if vrblvl > 0:
        print('in set_double_slice, equ :', equ, end='')
        print(', in idx :', idx, ', cff :', cff)
    phc = get_phcfun(vrblvl-1)
    aequ = pointer(c_int32(equ))
    bidx = pointer(c_int32(idx))
    cffs = (c_double * 2)()
    cffs[0] = c_double(cff[0])
    cffs[1] = c_double(cff[1])
    ccff = pointer(cffs)
    vrb = c_int32(vrblvl-1)
    if vrblvl > 0:
        print('-> set_double_slice calls phc', end='')
    retval = phc(43, aequ, bidx, ccff, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    return retval

def set_double_double_slice(equ, idx, cff, vrblvl=0):
    """
    Sets the coefficients of slicing equation with index equ
    at position idx to the value in double double precision:
    (cff[0], cff[1]) + (cff[2], cff[3])*complex(0, 1).
    """
    if vrblvl > 0:
        print('in set_double_double_slice, equ :', equ, end='')
        print(', in idx :', idx, ', cff :', cff)
    phc = get_phcfun(vrblvl-1)
    aequ = pointer(c_int32(equ))
    bidx = pointer(c_int32(idx))
    cffs = (c_double * 4)()
    cffs[0] = c_double(cff[0])
    cffs[1] = c_double(cff[1])
    cffs[2] = c_double(cff[2])
    cffs[3] = c_double(cff[3])
    ccff = pointer(cffs)
    vrb = c_int32(vrblvl-1)
    if vrblvl > 0:
        print('-> set_double_double_slice calls phc', end='')
    retval = phc(633, aequ, bidx, ccff, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    return retval

def set_quad_double_slice(equ, idx, cff, vrblvl=0):
    """
    Sets the coefficients of slicing equation with index equ
    at position idx to the value in quad double precision:
    (cff[0], cff[1], cff[2], cff[3])
    + (cff[4], cff[5], cff[6], cff[7])*complex(0, 1).
    """
    if vrblvl > 0:
        print('in set_double_double_slice, equ :', equ, end='')
        print(', in idx :', idx, ', cff :', cff)
    phc = get_phcfun(vrblvl-1)
    aequ = pointer(c_int32(equ))
    bidx = pointer(c_int32(idx))
    cffs = (c_double * 8)()
    cffs[0] = c_double(cff[0])
    cffs[1] = c_double(cff[1])
    cffs[2] = c_double(cff[2])
    cffs[3] = c_double(cff[3])
    cffs[4] = c_double(cff[4])
    cffs[5] = c_double(cff[5])
    cffs[6] = c_double(cff[6])
    cffs[7] = c_double(cff[7])
    ccff = pointer(cffs)
    vrb = c_int32(vrblvl-1)
    if vrblvl > 0:
        print('-> set_quad_double_slice calls phc', end='')
    retval = phc(663, aequ, bidx, ccff, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    return retval

def set_double_trace(first, vrblvl=0):
    """
    Sets the constant coefficient of the first slice,
    for use in the linear trace test in double precision.
    The integer first indicates if it is the first time or not.
    """
    if vrblvl > 0:
        print('in set_double_trace, first :', first)
    if first == 1:
        cff = [-1.0, 0.0]
    else:
        cff = [+1.0, 0.0]
    return set_double_slice(0, 0, cff, vrblvl)

def set_double_double_trace(first, vrblvl=0):
    """
    Sets the constant coefficient of the first slice,
    for use in the linear trace test in double double precision.
    The integer first indicates if it is the first time or not.
    """
    if vrblvl > 0:
        print('in set_double_double_trace, first :', first)
    if first == 1:
        cff = [-1.0, 0.0, 0.0, 0.0]
    else:
        cff = [+1.0, 0.0, 0.0, 0.0]
    return set_double_double_slice(0, 0, cff, vrblvl)

def set_quad_double_trace(first, vrblvl=0):
    """
    Sets the constant coefficient of the first slice,
    for use in the linear trace test in quad double precision.
    The integer first indicates if it is the first time or not.
    """
    if vrblvl > 0:
        print('in set_quad_double_trace, first :', first)
    if first == 1:
        cff = [-1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    else:
        cff = [+1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    return set_quad_double_slice(0, 0, cff, vrblvl)

def set_double_gammas(dim, vrblvl=0):
    """
    Sets the gamma constants in double precision
    for the sampler in the monodromy loops.
    Generates as many random complex constants as dim.
    """
    if vrblvl > 0:
        print('in set_double_gammas, dim :', dim)
    phc = get_phcfun(vrblvl-1)
    aidx = pointer(c_int32(0))
    bbb = pointer(c_int32(0))
    gamma = (c_double * 2)()
    cgamma = pointer(gamma)
    vrb = c_int32(vrblvl-1)
    retval = 0
    for i in range(dim):
        angle = uniform(0, 2*pi)
        rndgamma = exp(angle*complex(0, 1))
        gamma[0] = c_double(rndgamma.real)
        gamma[1] = c_double(rndgamma.imag)
        if vrblvl > 0:
            print('random gamma :', rndgamma)
        if vrblvl > 0:
            print('-> set_double_gammas calls phc', end='')
        aidx[0] = c_int32(i)
        retval = phc(44, aidx, bbb, cgamma, vrb)
        if vrblvl > 0:
            print(', return value :', retval)
    return retval

def set_double_double_gammas(dim, vrblvl=0):
    """
    Sets the gamma constants in double_double precision
    for the sampler in the monodromy loops.
    Generates as many random complex constants as dim.
    """
    if vrblvl > 0:
        print('in set_double_double_gammas, dim :', dim)
    phc = get_phcfun(vrblvl-1)
    aidx = pointer(c_int32(0))
    bbb = pointer(c_int32(0))
    gamma = (c_double * 4)()
    cgamma = pointer(gamma)
    vrb = c_int32(vrblvl-1)
    retval = 0
    for i in range(dim):
        angle = uniform(0, 2*pi)
        rndgamma = exp(angle*complex(0, 1))
        gamma[0] = c_double(rndgamma.real)
        gamma[1] = c_double(0.0)
        gamma[2] = c_double(rndgamma.imag)
        gamma[3] = c_double(0.0)
        if vrblvl > 0:
            print('random gamma :', rndgamma)
        if vrblvl > 0:
            print('-> set_double_double_gammas calls phc', end='')
        aidx[0] = c_int32(i)
        retval = phc(634, aidx, bbb, cgamma, vrb)
        if vrblvl > 0:
            print(', return value :', retval)
    return retval

def set_quad_double_gammas(dim, vrblvl=0):
    """
    Sets the gamma constants in quad_double precision
    for the sampler in the monodromy loops.
    Generates as many random complex constants as dim.
    """
    if vrblvl > 0:
        print('in set_quad_double_gammas, dim :', dim)
    phc = get_phcfun(vrblvl-1)
    aidx = pointer(c_int32(0))
    bbb = pointer(c_int32(0))
    gamma = (c_double * 8)()
    cgamma = pointer(gamma)
    vrb = c_int32(vrblvl-1)
    retval = 0
    for i in range(dim):
        angle = uniform(0, 2*pi)
        rndgamma = exp(angle*complex(0, 1))
        gamma[0] = c_double(rndgamma.real)
        gamma[1] = c_double(0.0)
        gamma[2] = c_double(0.0)
        gamma[3] = c_double(0.0)
        gamma[4] = c_double(rndgamma.imag)
        gamma[5] = c_double(0.0)
        gamma[6] = c_double(0.0)
        gamma[7] = c_double(0.0)
        if vrblvl > 0:
            print('random gamma :', rndgamma)
        if vrblvl > 0:
            print('-> set_quad_double_gammas calls phc', end='')
        aidx[0] = c_int32(i)
        retval = phc(664, aidx, bbb, cgamma, vrb)
        if vrblvl > 0:
            print(', return value :', retval)
    return retval

def swap_double_slices(vrblvl=0):
    """
    Swaps the slices with new ones to turn back in one monodromy loop
    in double precision.
    """
    if vrblvl > 0:
        print('in swap_double_slices ...')
    phc = get_phcfun(vrblvl-1)
    aaa = pointer(c_int32(0))
    bbb = pointer(c_int32(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl-1)
    if vrblvl > 0:
        print('-> swap_double_slices calls phc', end='')
    retval = phc(46, aaa, bbb, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    return retval

def swap_double_double_slices(vrblvl=0):
    """
    Swaps the slices with new ones to turn back in one monodromy loop
    in double double precision.
    """
    if vrblvl > 0:
        print('in swap_double_double_slices ...')
    phc = get_phcfun(vrblvl-1)
    aaa = pointer(c_int32(0))
    bbb = pointer(c_int32(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl-1)
    if vrblvl > 0:
        print('-> swap_double_double_slices calls phc', end='')
    retval = phc(636, aaa, bbb, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    return retval

def swap_quad_double_slices(vrblvl=0):
    """
    Swaps the slices with new ones to turn back in one monodromy loop
    in quad double precision.
    """
    if vrblvl > 0:
        print('in swap_quad_double_slices ...')
    phc = get_phcfun(vrblvl-1)
    aaa = pointer(c_int32(0))
    bbb = pointer(c_int32(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl-1)
    if vrblvl > 0:
        print('-> swap_quad_double_slices calls phc', end='')
    retval = phc(666, aaa, bbb, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    return retval

def new_double_slices(dim, nvr, vrblvl=0):
    """
    Generates dim random hyperplanes in double precision
    where nvr is the number of variables.
    """
    if vrblvl > 0:
        print('in new_double_slices, dim :', dim, end='')
        print(', nvr :', nvr)
    fail = 0
    for row in range(dim):
        for col in range(nvr+1):
            angle = uniform(0, 2*pi)
            rancff = exp(angle*complex(0, 1))
            cff = [rancff.real, rancff.imag]
            fail = fail + set_double_slice(row, col, cff, vrblvl)
    return fail

def new_double_double_slices(dim, nvr, vrblvl=0):
    """
    Generates dim random hyperplanes in double double precision
    where nvr is the number of variables.
    """
    if vrblvl > 0:
        print('in new_double_double_slices, dim :', dim, end='')
        print(', nvr :', nvr)
    fail = 0
    for row in range(dim):
        for col in range(nvr+1):
            angle = uniform(0, 2*pi)
            rancff = exp(angle*complex(0, 1))
            cff = [rancff.real, 0.0, rancff.imag, 0.0]
            fail = fail + set_double_double_slice(row, col, cff, vrblvl)
    return fail

def new_quad_double_slices(dim, nvr, vrblvl=0):
    """
    Generates dim random hyperplanes in quad double precision
    where nvr is the number of variables.
    """
    if vrblvl > 0:
        print('in new_quad_double_slices, dim :', dim, end='')
        print(', nvr :', nvr)
    fail = 0
    for row in range(dim):
        for col in range(nvr+1):
            angle = uniform(0, 2*pi)
            rancff = exp(angle*complex(0, 1))
            cff = [rancff.real, 0.0, 0.0, 0.0, rancff.imag, 0.0, 0.0, 0.0]
            fail = fail + set_quad_double_slice(row, col, cff, vrblvl)
    return fail

def double_witness_sample(vrblvl=0):
    """
    Computes a new witness set for a new set of slices,
    in double precision.
    """
    if vrblvl > 0:
        print('in double_witness_sample ...')
    phc = get_phcfun(vrblvl-1)
    aaa = pointer(c_int32(0))
    bbb = pointer(c_int32(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl-1)
    if vrblvl > 0:
        print('-> double_witness_sample calls phc', end='')
    retval = phc(45, aaa, bbb, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    return retval

def double_double_witness_sample(vrblvl=0):
    """
    Computes a new witness set for a new set of slices,
    in double double precision.
    """
    if vrblvl > 0:
        print('in double_double_witness_sample ...')
    phc = get_phcfun(vrblvl-1)
    aaa = pointer(c_int32(0))
    bbb = pointer(c_int32(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl-1)
    if vrblvl > 0:
        print('-> double_double_witness_sample calls phc', end='')
    retval = phc(635, aaa, bbb, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    return retval

def quad_double_witness_sample(vrblvl=0):
    """
    Computes a new witness set for a new set of slices,
    in quad double precision.
    """
    if vrblvl > 0:
        print('in quad_double_witness_sample ...')
    phc = get_phcfun(vrblvl-1)
    aaa = pointer(c_int32(0))
    bbb = pointer(c_int32(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl-1)
    if vrblvl > 0:
        print('-> quad_double_witness_sample calls phc', end='')
    retval = phc(665, aaa, bbb, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    return retval

def copy_double_witness_set(vrblvl=0):
    """
    Copies the witness set from the sampler to the system set
    in double precision.
    """
    if vrblvl > 0:
        print('in copy_double_witness_set ...')
    phc = get_phcfun(vrblvl-1)
    aaa = pointer(c_int32(0))
    bbb = pointer(c_int32(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl-1)
    if vrblvl > 0:
        print('-> copy_double_witness_set calls phc', end='')
    retval = phc(47, aaa, bbb, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    return retval

def copy_double_double_witness_set(vrblvl=0):
    """
    Copies the witness set from the sampler to the system set
    in double double precision.
    """
    if vrblvl > 0:
        print('in copy_double_double_witness_set ...')
    phc = get_phcfun(vrblvl-1)
    aaa = pointer(c_int32(0))
    bbb = pointer(c_int32(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl-1)
    if vrblvl > 0:
        print('-> copy_double_double_witness_set calls phc', end='')
    retval = phc(637, aaa, bbb, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    return retval

def copy_quad_double_witness_set(vrblvl=0):
    """
    Copies the witness set from the sampler to the system set
    in quad double precision.
    """
    if vrblvl > 0:
        print('in copy_quad_double_witness_set ...')
    phc = get_phcfun(vrblvl-1)
    aaa = pointer(c_int32(0))
    bbb = pointer(c_int32(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl-1)
    if vrblvl > 0:
        print('-> copy_quad_double_witness_set calls phc', end='')
    retval = phc(667, aaa, bbb, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    return retval

def copy_double_laurent_witness_set(vrblvl=0):
    """
    Copies the witness set from the sampler to the laurent system set
    in double precision.
    """
    if vrblvl > 0:
        print('in copy_double_witness_set ...')
    phc = get_phcfun(vrblvl-1)
    aaa = pointer(c_int32(0))
    bbb = pointer(c_int32(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl-1)
    if vrblvl > 0:
        print('-> copy_double_laurent_witness_set calls phc', end='')
    retval = phc(807, aaa, bbb, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    return retval

def copy_double_double_laurent_witness_set(vrblvl=0):
    """
    Copies the witness set from the sampler to the laurent system set
    in double double precision.
    """
    if vrblvl > 0:
        print('in copy_double_double_witness_set ...')
    phc = get_phcfun(vrblvl-1)
    aaa = pointer(c_int32(0))
    bbb = pointer(c_int32(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl-1)
    if vrblvl > 0:
        print('-> copy_double_double_laurent_witness_set calls phc', end='')
    retval = phc(808, aaa, bbb, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    return retval

def copy_quad_double_laurent_witness_set(vrblvl=0):
    """
    Copies the witness set from the sampler to the laurent system set
    in quad double precision.
    """
    if vrblvl > 0:
        print('in copy_quad_double_witness_set ...')
    phc = get_phcfun(vrblvl-1)
    aaa = pointer(c_int32(0))
    bbb = pointer(c_int32(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl-1)
    if vrblvl > 0:
        print('-> copy_quad_double_laurent_witness_set calls phc', end='')
    retval = phc(809, aaa, bbb, ccc, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    return retval

def double_witness_track(islaurent=False, vrblvl=0):
    """
    Tracks as many paths as set in double precision,
    as many as the size of the witness set.
    If the system is laurent, then islaurent must be True.
    """
    if vrblvl > 0:
        print('in double_witness_track, islaurent :', islaurent)
    double_witness_sample(vrblvl)
    swap_double_slices(vrblvl)
    if islaurent:
        copy_double_laurent_witness_set(vrblvl)
    else:
        copy_double_witness_set(vrblvl)

def double_double_witness_track(islaurent=False, vrblvl=0):
    """
    Tracks as many paths as set in double double precision,
    as many as the size of the witness set.
    If the system is laurent, then islaurent must be True.
    """
    if vrblvl > 0:
        print('in double_double_witness_track, islaurent :', islaurent)
    double_double_witness_sample(vrblvl)
    swap_double_double_slices(vrblvl)
    if islaurent:
        copy_double_double_laurent_witness_set(vrblvl)
    else:
        copy_double_double_witness_set(vrblvl)

def quad_double_witness_track(islaurent=False, vrblvl=0):
    """
    Tracks as many paths as set in quad double precision,
    as many as the size of the witness set.
    If the system is laurent, then islaurent must be True.
    """
    if vrblvl > 0:
        print('in quad_double_witness_track, islaurent :', islaurent)
    quad_double_witness_sample(vrblvl)
    swap_quad_double_slices(vrblvl)
    if islaurent:
        copy_quad_double_laurent_witness_set(vrblvl)
    else:
        copy_quad_double_witness_set(vrblvl)

def double_trace_grid_diagnostics(vrblvl=0):
    """
    Returns the maximal error on the samples in the trace grid
    and the minimal distance between the samples in the trace grid,
    computed in double precision.
    """
    if vrblvl > 0:
        print('in double_trace_grid_diagnostics ...')
    phc = get_phcfun(vrblvl-1)
    aaa = pointer(c_int32(0))
    bbb = pointer(c_int32(0))
    pars = (c_double * 2)()
    cpars = pointer(pars)
    vrb = c_int32(vrblvl-1)
    if vrblvl > 0:
        print('-> double_trace_grid_diagnostics calls phc', end='')
    retval = phc(56, aaa, bbb, cpars, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    vals = cpars[:2]
    err = vals[0][0]
    dis = vals[0][1]
    if vrblvl > 0:
        print('err :', err, ', dis :', dis)
    return (err, dis)

def double_double_trace_grid_diagnostics(vrblvl=0):
    """
    Returns the maximal error on the samples in the trace grid
    and the minimal distance between the samples in the trace grid,
    computed in double double precision.
    """
    if vrblvl > 0:
        print('in double_double_trace_grid_diagnostics ...')
    phc = get_phcfun(vrblvl-1)
    aaa = pointer(c_int32(0))
    bbb = pointer(c_int32(0))
    pars = (c_double * 2)()
    cpars = pointer(pars)
    vrb = c_int32(vrblvl-1)
    if vrblvl > 0:
        print('-> double_double_trace_grid_diagnostics calls phc', end='')
    retval = phc(646, aaa, bbb, cpars, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    vals = cpars[:2]
    err = vals[0][0]
    dis = vals[0][1]
    if vrblvl > 0:
        print('err :', err, ', dis :', dis)
    return (err, dis)

def quad_double_trace_grid_diagnostics(vrblvl=0):
    """
    Returns the maximal error on the samples in the trace grid
    and the minimal distance between the samples in the trace grid,
    computed in quad double precision.
    """
    if vrblvl > 0:
        print('in quad_double_trace_grid_diagnostics ...')
    phc = get_phcfun(vrblvl-1)
    aaa = pointer(c_int32(0))
    bbb = pointer(c_int32(0))
    pars = (c_double * 2)()
    cpars = pointer(pars)
    vrb = c_int32(vrblvl-1)
    if vrblvl > 0:
        print('-> quad_double_trace_grid_diagnostics calls phc', end='')
    retval = phc(676, aaa, bbb, cpars, vrb)
    if vrblvl > 0:
        print(', return value :', retval)
    vals = cpars[:2]
    err = vals[0][0]
    dis = vals[0][1]
    if vrblvl > 0:
        print('err :', err, ', dis :', dis)
    return (err, dis)

def double_trace_sum_difference(labels, vrblvl=0):
    """
    Returns the difference between the actual sum at the samples
    defined by the labels to the generic points of a factor
    and the trace sum, in double precision.
    """
    if vrblvl > 0:
        print('in double_trace_sum_difference, labels :', labels)
    phc = get_phcfun(vrblvl-1)
    adim = pointer(c_int32(len(labels)))
    fac = (c_int32 * len(labels))()
    for (idx, label) in enumerate(labels):
        fac[idx] = c_int32(label)
    bfac = pointer(fac)
    cdif = pointer(c_double(0.0))
    vrb = c_int32(vrblvl-1)
    if vrblvl > 0:
        print('-> double_trace_sum_difference calls phc', end='')
    retval = phc(57, adim, bfac, cdif, vrb)
    result = cdif[0]
    if vrblvl > 0:
        print(', return value :', retval)
        print('trace sum difference :', result)
    return result

def double_double_trace_sum_difference(labels, vrblvl=0):
    """
    Returns the difference between the actual sum at the samples
    defined by the labels to the generic points of a factor
    and the trace sum, in double double precision.
    """
    if vrblvl > 0:
        print('in double_double_trace_sum_difference, labels :', labels)
    phc = get_phcfun(vrblvl-1)
    adim = pointer(c_int32(len(labels)))
    fac = (c_int32 * len(labels))()
    for (idx, label) in enumerate(labels):
        fac[idx] = c_int32(label)
    bfac = pointer(fac)
    cdif = pointer(c_double(0.0))
    vrb = c_int32(vrblvl-1)
    if vrblvl > 0:
        print('-> double_double_trace_sum_difference calls phc', end='')
    retval = phc(647, adim, bfac, cdif, vrb)
    result = cdif[0]
    if vrblvl > 0:
        print(', return value :', retval)
        print('trace sum difference :', result)
    return result

def quad_double_trace_sum_difference(labels, vrblvl=0):
    """
    Returns the difference between the actual sum at the samples
    defined by the labels to the generic points of a factor
    and the trace sum, in quad double precision.
    """
    if vrblvl > 0:
        print('in quad_double_trace_sum_difference, labels :', labels)
    phc = get_phcfun(vrblvl-1)
    adim = pointer(c_int32(len(labels)))
    fac = (c_int32 * len(labels))()
    for (idx, label) in enumerate(labels):
        fac[idx] = c_int32(label)
    bfac = pointer(fac)
    cdif = pointer(c_double(0.0))
    vrb = c_int32(vrblvl-1)
    if vrblvl > 0:
        print('-> quad_double_trace_sum_difference calls phc', end='')
    retval = phc(677, adim, bfac, cdif, vrb)
    result = cdif[0]
    if vrblvl > 0:
        print(', return value :', retval)
        print('trace sum difference :', result)
    return result

def double_loop_permutation(deg, vrblvl=0):
    """
    Returns the permutation using the solution most recently computed,
    for a set of degree deg, after a loop in double precision.
    """
    if vrblvl > 0:
        print('in double_loop_permutation, deg :', deg)
    phc = get_phcfun(vrblvl-1)
    adeg = pointer(c_int32(deg))
    perm = (c_int32 * deg)()
    bperm = pointer(perm)
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl-1)
    if vrblvl > 0:
        print('-> double_loop permutation calls phc', end='')
    retval = phc(52, adeg, bperm, ccc, vrb)
    vals = bperm[0:deg]
    result = []
    for idx in range(deg):
        result.append(int(vals[0][idx]))
    if vrblvl > 0:
        print(', return value :', retval)
        print('permutation :', result)
    return result

def double_double_loop_permutation(deg, vrblvl=0):
    """
    Returns the permutation using the solution most recently computed,
    for a set of degree deg, after a loop in double double precision.
    """
    if vrblvl > 0:
        print('in double_double_loop_permutation, deg :', deg)
    phc = get_phcfun(vrblvl-1)
    adeg = pointer(c_int32(deg))
    perm = (c_int32 * deg)()
    bperm = pointer(perm)
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl-1)
    if vrblvl > 0:
        print('-> double_double_loop permutation calls phc', end='')
    retval = phc(642, adeg, bperm, ccc, vrb)
    vals = bperm[0:deg]
    result = []
    for idx in range(deg):
        result.append(int(vals[0][idx]))
    if vrblvl > 0:
        print(', return value :', retval)
        print('permutation :', result)
    return result

def quad_double_loop_permutation(deg, vrblvl=0):
    """
    Returns the permutation using the solution most recently computed,
    for a set of degree deg, after a loop in quad double precision.
    """
    if vrblvl > 0:
        print('in quad_double_loop_permutation, deg :', deg)
    phc = get_phcfun(vrblvl-1)
    adeg = pointer(c_int32(deg))
    perm = (c_int32 * deg)()
    bperm = pointer(perm)
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl-1)
    if vrblvl > 0:
        print('-> quad_double_loop permutation calls phc', end='')
    retval = phc(672, adeg, bperm, ccc, vrb)
    vals = bperm[0:deg]
    result = []
    for idx in range(deg):
        result.append(int(vals[0][idx]))
    if vrblvl > 0:
        print(', return value :', retval)
        print('permutation :', result)
    return result

def double_factor_count(vrblvl=0):
    """
    Returns the number of factors computed in double precision.
    """
    if vrblvl > 0:
        print('in double_factor_count ...')
    phc = get_phcfun(vrblvl-1)
    anbr = pointer(c_int32(0))
    bbb = pointer(c_int32(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl-1)
    if vrblvl > 0:
        print('-> double_factor_count calls phc', end='')
    retval = phc(68, anbr, bbb, ccc, vrb)
    result = anbr[0]
    if vrblvl > 0:
        print(', return value :', retval)
        print('number of factors :', result)
    return result

def double_double_factor_count(vrblvl=0):
    """
    Returns the number of factors computed in double double precision.
    """
    if vrblvl > 0:
        print('in double_double_factor_count ...')
    phc = get_phcfun(vrblvl-1)
    anbr = pointer(c_int32(0))
    bbb = pointer(c_int32(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl-1)
    if vrblvl > 0:
        print('-> double_double_factor_count calls phc', end='')
    retval = phc(656, anbr, bbb, ccc, vrb)
    result = anbr[0]
    if vrblvl > 0:
        print(', return value :', retval)
        print('number of factors :', result)
    return result

def quad_double_factor_count(vrblvl=0):
    """
    Returns the number of factors computed in quad double precision.
    """
    if vrblvl > 0:
        print('in quad_double_factor_count ...')
    phc = get_phcfun(vrblvl-1)
    anbr = pointer(c_int32(0))
    bbb = pointer(c_int32(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl-1)
    if vrblvl > 0:
        print('-> quad_double_factor_count calls phc', end='')
    retval = phc(686, anbr, bbb, ccc, vrb)
    result = anbr[0]
    if vrblvl > 0:
        print(', return value :', retval)
        print('number of factors :', result)
    return result

def update_double_decomposition(deg, perm, vrblvl=0):
    """
    Updates the decomposition with a permutation of deg elements
    computed in double precision.
    Returns the tuple with the previous and the new number of factors.
    """
    if vrblvl > 0:
        print('in update_double_decomposition, deg :', deg)
        print('permutation :', perm)
    phc = get_phcfun(vrblvl-1)
    nbrs = (c_int32 * 2)()
    nbrs[0] = len(perm)
    anbr = pointer(nbrs)
    iperm = (c_int32 * len(perm))()
    for (idx, nbr) in enumerate(perm):
        iperm[idx] = c_int32(nbr)
    bperm = pointer(iperm)
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl-1)
    if vrblvl > 0:
        print('-> update_double_decomposition_count calls phc', end='')
    retval = phc(53, anbr, bperm, ccc, vrb)
    vals = anbr[:2]
    result = (vals[0][0], vals[0][1])
    if vrblvl > 0:
        print(', return value :', retval)
        print('number of factors :', result)
    return result

def update_double_double_decomposition(deg, perm, vrblvl=0):
    """
    Updates the decomposition with a permutation of deg elements
    computed in double double precision.
    Returns the tuple with the previous and the new number of factors.
    """
    if vrblvl > 0:
        print('in update_double_double_decomposition, deg :', deg)
        print('permutation :', perm)
    phc = get_phcfun(vrblvl-1)
    nbrs = (c_int32 * 2)()
    nbrs[0] = len(perm)
    anbr = pointer(nbrs)
    iperm = (c_int32 * len(perm))()
    for (idx, nbr) in enumerate(perm):
        iperm[idx] = c_int32(nbr)
    bperm = pointer(iperm)
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl-1)
    if vrblvl > 0:
        print('-> update_double_double_decomposition_count calls phc', end='')
    retval = phc(643, anbr, bperm, ccc, vrb)
    vals = anbr[:2]
    result = (vals[0][0], vals[0][1])
    if vrblvl > 0:
        print(', return value :', retval)
        print('number of factors :', result)
    return result

def update_quad_double_decomposition(deg, perm, vrblvl=0):
    """
    Updates the decomposition with a permutation of deg elements
    computed in quad double precision.
    Returns the tuple with the previous and the new number of factors.
    """
    if vrblvl > 0:
        print('in update_quad_double_decomposition, deg :', deg)
        print('permutation :', perm)
    phc = get_phcfun(vrblvl-1)
    nbrs = (c_int32 * 2)()
    nbrs[0] = len(perm)
    anbr = pointer(nbrs)
    iperm = (c_int32 * len(perm))()
    for (idx, nbr) in enumerate(perm):
        iperm[idx] = c_int32(nbr)
    bperm = pointer(iperm)
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl-1)
    if vrblvl > 0:
        print('-> update_quad_double_decomposition_count calls phc', end='')
    retval = phc(673, anbr, bperm, ccc, vrb)
    vals = anbr[:2]
    result = (vals[0][0], vals[0][1])
    if vrblvl > 0:
        print(', return value :', retval)
        print('number of factors :', result)
    return result

def double_trace_test(vrblvl=0):
    """
    Runs the trace test on the decompostion in double precision,
    returns True if certified, otherwise returns False.
    """
    if vrblvl > 0:
        print('in double_trace_test ...')
    phc = get_phcfun(vrblvl-1)
    adone = pointer(c_int32(0))
    bbb = pointer(c_int32(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl-1)
    if vrblvl > 0:
        print('-> double_trace_test calls phc', end='')
    retval = phc(55, adone, bbb, ccc, vrb)
    result = adone[0] == 1
    if vrblvl > 0:
        print(', return value :', retval)
        print('certified :', result)
    return result

def double_double_trace_test(vrblvl=0):
    """
    Runs the trace test on the decompostion in double double precision,
    returns True if certified, otherwise returns False.
    """
    if vrblvl > 0:
        print('in double_double_trace_test ...')
    phc = get_phcfun(vrblvl-1)
    adone = pointer(c_int32(0))
    bbb = pointer(c_int32(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl-1)
    if vrblvl > 0:
        print('-> double_double_trace_test calls phc', end='')
    retval = phc(645, adone, bbb, ccc, vrb)
    result = adone[0] == 1
    if vrblvl > 0:
        print(', return value :', retval)
        print('certified :', result)
    return result

def quad_double_trace_test(vrblvl=0):
    """
    Runs the trace test on the decompostion in quad double precision,
    returns True if certified, otherwise returns False.
    """
    if vrblvl > 0:
        print('in quad_double_trace_test ...')
    phc = get_phcfun(vrblvl-1)
    adone = pointer(c_int32(0))
    bbb = pointer(c_int32(0))
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl-1)
    if vrblvl > 0:
        print('-> quad_double_trace_test calls phc', end='')
    retval = phc(675, adone, bbb, ccc, vrb)
    result = adone[0] == 1
    if vrblvl > 0:
        print(', return value :', retval)
        print('certified :', result)
    return result

def double_witness_points(idx, deg, vrblvl=0):
    """
    Given an index idx of an irreducible component,
    computed in double precision,
    returns the labels of the witness points that span the component.
    The input deg is the upper bound on the degree of a factor.
    The degree of the factor is the length of the returned list.
    """
    if vrblvl > 0:
        print('in double_witness_points, idx :', idx, end='')
        print(', deg :', deg)
    phc = get_phcfun(vrblvl-1)
    aidx = pointer(c_int32(idx))
    witpts = (c_int32 * deg)()
    bwit = pointer(witpts)
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl-1)
    if vrblvl > 0:
        print('-> double_witness_points calls phc', end='')
    retval = phc(69, aidx, bwit, ccc, vrb)
    wdeg = aidx[0]
    if vrblvl > 0:
        print(', the return value :', retval)
        print('the degree of factor', idx, ':', wdeg)
    result = []
    vals = bwit[:wdeg]
    for index in range(wdeg):
        result.append(vals[0][index])
    if vrblvl > 0:
        print('labels of witness points :', result)
    return result

def double_double_witness_points(idx, deg, vrblvl=0):
    """
    Given an index idx of an irreducible component,
    computed in double double precision,
    returns the labels of the witness points that span the component.
    The input deg is the upper bound on the degree of a factor.
    The degree of the factor is the length of the returned list.
    """
    if vrblvl > 0:
        print('in double_double_witness_points, idx :', idx, end='')
        print(', deg :', deg)
    phc = get_phcfun(vrblvl-1)
    aidx = pointer(c_int32(idx))
    witpts = (c_int32 * deg)()
    bwit = pointer(witpts)
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl-1)
    if vrblvl > 0:
        print('-> double_double_witness_points calls phc', end='')
    retval = phc(657, aidx, bwit, ccc, vrb)
    wdeg = aidx[0]
    if vrblvl > 0:
        print(', the return value :', retval)
        print('the degree of factor', idx, ':', wdeg)
    result = []
    vals = bwit[:wdeg]
    for index in range(wdeg):
        result.append(vals[0][index])
    if vrblvl > 0:
        print('labels of witness points :', result)
    return result

def quad_double_witness_points(idx, deg, vrblvl=0):
    """
    Given an index idx of an irreducible component,
    computed in quad double precision,
    returns the labels of the witness points that span the component.
    The input deg is the upper bound on the degree of a factor.
    The degree of the factor is the length of the returned list.
    """
    if vrblvl > 0:
        print('in quad_double_witness_points, idx :', idx, end='')
        print(', deg :', deg)
    phc = get_phcfun(vrblvl-1)
    aidx = pointer(c_int32(idx))
    witpts = (c_int32 * deg)()
    bwit = pointer(witpts)
    ccc = pointer(c_double(0.0))
    vrb = c_int32(vrblvl-1)
    if vrblvl > 0:
        print('-> quad_double_witness_points calls phc', end='')
    retval = phc(687, aidx, bwit, ccc, vrb)
    wdeg = aidx[0]
    if vrblvl > 0:
        print(', the return value :', retval)
        print('the degree of factor', idx, ':', wdeg)
    result = []
    vals = bwit[:wdeg]
    for index in range(wdeg):
        result.append(vals[0][index])
    if vrblvl > 0:
        print('labels of witness points :', result)
    return result

def double_decomposition(deg, vrblvl=0):
    """
    Returns the decomposition as a list of labels of witness points
    on the components, computed in double precision.
    """
    if vrblvl > 0:
        print('in double_decomposition, deg :', deg)
    nbr = double_factor_count(vrblvl)
    result = []
    for idx in range(1, nbr+1):
        witpts = double_witness_points(idx, deg, vrblvl)
        trace = double_trace_sum_difference(witpts, vrblvl)
        result.append((witpts, trace))
    return result

def double_double_decomposition(deg, vrblvl=0):
    """
    Returns the decomposition as a list of labels of witness points
    on the components, computed in double double precision.
    """
    if vrblvl > 0:
        print('in double_double_decomposition, deg :', deg)
    nbr = double_double_factor_count(vrblvl)
    result = []
    for idx in range(1, nbr+1):
        witpts = double_double_witness_points(idx, deg, vrblvl)
        trace = double_double_trace_sum_difference(witpts, vrblvl)
        result.append((witpts, trace))
    return result

def quad_double_decomposition(deg, vrblvl=0):
    """
    Returns the decomposition as a list of labels of witness points
    on the components, computed in quad double precision.
    """
    if vrblvl > 0:
        print('in quad_double_decomposition, deg :', deg)
    nbr = quad_double_factor_count(vrblvl)
    result = []
    for idx in range(1, nbr+1):
        witpts = quad_double_witness_points(idx, deg, vrblvl)
        trace = quad_double_trace_sum_difference(witpts, vrblvl)
        result.append((witpts, trace))
    return result

def write_factorization(deco):
    """
    Writes the decomposition in deco.
    """
    for (idx, factor) in enumerate(deco):
        print('  factor', idx+1, ':', factor)

def double_monodromy_breakup(embsys, esols, dim, \
    islaurent=False, verbose=False, nbloops=20, vrblvl=0):
    r"""
    Applies the monodromy breakup algorithm in double precision
    to factor the *dim*-dimensional algebraic set represented by the
    embedded system *embsys* and its solutions *esols*.
    If the embedded polynomial system is a laurent system,
    then *islaurent* must be True.
    If *verbose* is False, then no output is written.
    The value of *nbloops* equals the maximum number of loops.
    """
    if vrblvl > 0:
        print('in double_monodromy_breakup', end='')
        print(', dim :', dim, ', nbloops :', nbloops)
        print(', islaurent :', islaurent, ', verbose :', verbose)
        print('the embedded polynomials :')
        for pol in embsys:
            print(pol)
        print('the generic points :')
        for (idx, sol) in enumerate(esols):
            print('Solution', idx+1, ':')
            print(sol)
    if verbose:
        print('... running monodromy loops in double precision ...')
    deg = len(esols)
    nvr = len(embsys)
    set_double_verbose(verbose, vrblvl)
    if islaurent:
        set_double_laurent_witness_set(nvr, dim, embsys, esols, vrblvl-1)
        double_assign_labels(nvr, vrblvl)
        initialize_double_laurent_sampler(dim, vrblvl)
    else:
        set_double_witness_set(nvr, dim, embsys, esols, vrblvl-1)
        double_assign_labels(nvr, vrblvl)
        initialize_double_sampler(dim, vrblvl)
    initialize_double_monodromy(nbloops, deg, dim, vrblvl)
    preset_double_solutions(vrblvl)
    if verbose:
        print('... initializing the grid for the linear trace ...')
    for i in range(1, 3):
        set_double_trace(i, vrblvl)
        set_double_gammas(nvr, vrblvl)
        double_witness_track(islaurent, vrblvl)
        preset_double_solutions(vrblvl)
        reset_double_solutions(vrblvl)
        swap_double_slices(vrblvl)
    (err, dis) = double_trace_grid_diagnostics(vrblvl)
    if verbose:
        print('The diagnostics of the trace grid :')
        print('  largest error on the samples :', err)
        print('  smallest distance between the samples :', dis)
    for i in range(1, nbloops+1):
        if verbose:
            print(f'... starting loop {i} ...')
        new_double_slices(dim, nvr, vrblvl)
        set_double_gammas(nvr)
        double_witness_track(islaurent, vrblvl)
        clear_double_solutions(vrblvl-1)
        set_double_gammas(nvr)
        double_witness_track(islaurent, vrblvl)
        preset_double_solutions(vrblvl)
        perm = double_loop_permutation(deg, vrblvl)
        if verbose:
            print('new permutation :', perm)
        nb0 = double_factor_count(vrblvl)
        nf0, nf1 = update_double_decomposition(deg, perm, vrblvl)
        if vrblvl > 0:
            print('nf0 :', nf0, ', nf1 :', nf1)
        nb1 = double_factor_count(vrblvl)
        if verbose:
            print(f'number of factors : {nb0} -> {nb1}')
            deco = double_decomposition(deg, vrblvl)
            print('the decomposition :')
            write_factorization(deco)
        done = double_trace_test(vrblvl)
        if done:
            break
        reset_double_solutions(vrblvl)
    return double_decomposition(deg, vrblvl)

def double_double_monodromy_breakup(embsys, esols, dim, \
    islaurent=False, verbose=False, nbloops=20, vrblvl=0):
    r"""
    Applies the monodromy breakup algorithm in double double precision
    to factor the *dim*-dimensional algebraic set represented by the
    embedded system *embsys* and its solutions *esols*.
    If the embedded polynomial system is a laurent system,
    then *islaurent* must be True.
    If *verbose* is False, then no output is written.
    The value of *nbloops* equals the maximum number of loops.
    """
    if vrblvl > 0:
        print('in double_double_monodromy_breakup', end='')
        print(', dim :', dim, ', nbloops :', nbloops)
        print(', islaurent :', islaurent, ', verbose :', verbose)
        print('the embedded polynomials :')
        for pol in embsys:
            print(pol)
        print('the generic points :')
        for (idx, sol) in enumerate(esols):
            print('Solution', idx+1, ':')
            print(sol)
    if verbose:
        print('... running monodromy loops in double double precision ...')
    deg = len(esols)
    nvr = len(embsys)
    set_double_double_verbose(verbose, vrblvl)
    if islaurent:
        set_double_double_laurent_witness_set(nvr, dim, embsys, esols, \
            vrblvl-1)
        double_double_assign_labels(nvr, vrblvl)
        initialize_double_double_laurent_sampler(dim, vrblvl)
    else:
        set_double_double_witness_set(nvr, dim, embsys, esols, vrblvl-1)
        double_double_assign_labels(nvr, vrblvl)
        initialize_double_double_sampler(dim, vrblvl)
    initialize_double_double_monodromy(nbloops, deg, dim, vrblvl)
    preset_double_double_solutions(vrblvl)
    if verbose:
        print('... initializing the grid for the linear trace ...')
    for i in range(1, 3):
        set_double_double_trace(i, vrblvl)
        set_double_double_gammas(nvr, vrblvl)
        double_double_witness_track(islaurent, vrblvl)
        preset_double_double_solutions(vrblvl)
        reset_double_double_solutions(vrblvl)
        swap_double_double_slices(vrblvl)
    (err, dis) = double_double_trace_grid_diagnostics(vrblvl)
    if verbose:
        print('The diagnostics of the trace grid :')
        print('  largest error on the samples :', err)
        print('  smallest distance between the samples :', dis)
    for i in range(1, nbloops+1):
        if verbose:
            print(f'... starting loop {i} ...')
        new_double_double_slices(dim, nvr, vrblvl)
        set_double_double_gammas(nvr)
        double_double_witness_track(islaurent, vrblvl)
        clear_double_double_solutions(vrblvl-1)
        set_double_double_gammas(nvr)
        double_double_witness_track(islaurent, vrblvl)
        preset_double_double_solutions(vrblvl)
        perm = double_double_loop_permutation(deg, vrblvl)
        if verbose:
            print('new permutation :', perm)
        nb0 = double_double_factor_count(vrblvl)
        nf0, nf1 = update_double_double_decomposition(deg, perm, vrblvl)
        if vrblvl > 0:
            print('nf0 :', nf0, ', nf1 :', nf1)
        nb1 = double_double_factor_count(vrblvl)
        if verbose:
            print(f'number of factors : {nb0} -> {nb1}')
            deco = double_double_decomposition(deg, vrblvl)
            print('the decomposition :')
            write_factorization(deco)
        done = double_double_trace_test(vrblvl)
        if done:
            break
        reset_double_double_solutions(vrblvl)
    return double_double_decomposition(deg, vrblvl)

def quad_double_monodromy_breakup(embsys, esols, dim, \
    islaurent=False, verbose=False, nbloops=20, vrblvl=0):
    r"""
    Applies the monodromy breakup algorithm in quad double precision
    to factor the *dim*-dimensional algebraic set represented by the
    embedded system *embsys* and its solutions *esols*.
    If the embedded polynomial system is a laurent system,
    then *islaurent* must be True.
    If *verbose* is False, then no output is written.
    The value of *nbloops* equals the maximum number of loops.
    """
    if vrblvl > 0:
        print('in quad_double_monodromy_breakup', end='')
        print(', dim :', dim, ', nbloops :', nbloops)
        print(', islaurent :', islaurent, ', verbose :', verbose)
        print('the embedded polynomials :')
        for pol in embsys:
            print(pol)
        print('the generic points :')
        for (idx, sol) in enumerate(esols):
            print('Solution', idx+1, ':')
            print(sol)
    if verbose:
        print('... running monodromy loops in quad double precision ...')
    deg = len(esols)
    nvr = len(embsys)
    set_quad_double_verbose(verbose, vrblvl)
    if islaurent:
        set_quad_double_laurent_witness_set(nvr, dim, embsys, esols, \
            vrblvl-1)
        quad_double_assign_labels(nvr, vrblvl)
        initialize_quad_double_laurent_sampler(dim, vrblvl)
    else:
        set_quad_double_witness_set(nvr, dim, embsys, esols, vrblvl-1)
        quad_double_assign_labels(nvr, vrblvl)
        initialize_quad_double_sampler(dim, vrblvl)
    initialize_quad_double_monodromy(nbloops, deg, dim, vrblvl)
    preset_quad_double_solutions(vrblvl)
    if verbose:
        print('... initializing the grid for the linear trace ...')
    for i in range(1, 3):
        set_quad_double_trace(i, vrblvl)
        set_quad_double_gammas(nvr, vrblvl)
        quad_double_witness_track(islaurent, vrblvl)
        preset_quad_double_solutions(vrblvl)
        reset_quad_double_solutions(vrblvl)
        swap_quad_double_slices(vrblvl)
    (err, dis) = quad_double_trace_grid_diagnostics(vrblvl)
    if verbose:
        print('The diagnostics of the trace grid :')
        print('  largest error on the samples :', err)
        print('  smallest distance between the samples :', dis)
    for i in range(1, nbloops+1):
        if verbose:
            print(f'... starting loop {i} ...')
        new_quad_double_slices(dim, nvr, vrblvl)
        set_quad_double_gammas(nvr)
        quad_double_witness_track(islaurent, vrblvl)
        clear_quad_double_solutions(vrblvl-1)
        set_quad_double_gammas(nvr)
        quad_double_witness_track(islaurent, vrblvl)
        preset_quad_double_solutions(vrblvl)
        perm = quad_double_loop_permutation(deg, vrblvl)
        if verbose:
            print('new permutation :', perm)
        nb0 = quad_double_factor_count(vrblvl)
        nf0, nf1 = update_quad_double_decomposition(deg, perm, vrblvl)
        if vrblvl > 0:
            print('nf0 :', nf0, ', nf1 :', nf1)
        nb1 = quad_double_factor_count(vrblvl)
        if verbose:
            print(f'number of factors : {nb0} -> {nb1}')
            deco = quad_double_decomposition(deg, vrblvl)
            print('the decomposition :')
            write_factorization(deco)
        done = quad_double_trace_test(vrblvl)
        if done:
            break
        reset_quad_double_solutions(vrblvl)
    return quad_double_decomposition(deg, vrblvl)

def test_double_assign_labels(vrblvl=0):
    """
    Tests the assigning of labels for solutions
    set in double precision.
    """
    if vrblvl > 0:
        print('in test_double_assign_labels ...')
    names = ['x', 'y']
    sols = []
    for idx in range(3):
        sol = make_solution(names, [idx, 0.0])
        sols.append(sol)
    if vrblvl > 0:
        print('the test solutions :')
        for (idx, sol) in enumerate(sols):
            print('Solution', idx+1, ':')
            print(sol)
    set_double_solutions(2, sols, vrblvl-1)
    double_assign_labels(2, vrblvl)
    newsols = get_double_solutions(vrblvl-1)
    if vrblvl > 0:
        print('after assigning labels :')
    fail = 0
    for (idx, sol) in enumerate(newsols):
        if vrblvl > 0:
            print('Solution', idx+1, ':')
            print(sol)
        (_, mval) = endmultiplicity(sol, vrblvl-1)
        if vrblvl > 0:
            print('mval :', mval, mval == idx+1)
        fail = fail + int(mval != idx+1)
    return fail

def test_double_double_assign_labels(vrblvl=0):
    """
    Tests the assigning of labels for solutions
    set in double double precision.
    """
    if vrblvl > 0:
        print('in test_double_double_assign_labels ...')
    names = ['x', 'y']
    sols = []
    for idx in range(3):
        sol = make_solution(names, [idx, 0.0])
        sols.append(sol)
    if vrblvl > 0:
        print('the test solutions :')
        for (idx, sol) in enumerate(sols):
            print('Solution', idx+1, ':')
            print(sol)
    set_double_double_solutions(2, sols, vrblvl-1)
    double_double_assign_labels(2, vrblvl)
    newsols = get_double_double_solutions(vrblvl-1)
    if vrblvl > 0:
        print('after assigning labels :')
    fail = 0
    for (idx, sol) in enumerate(newsols):
        if vrblvl > 0:
            print('Solution', idx+1, ':')
            print(sol)
        (_, mval) = endmultiplicity(sol, vrblvl-1)
        if vrblvl > 0:
            print('mval :', mval, mval == idx+1)
        fail = fail + int(mval != idx+1)
    return fail

def test_quad_double_assign_labels(vrblvl=0):
    """
    Tests the assigning of labels for solutions
    set in quad double precision.
    """
    if vrblvl > 0:
        print('in test_quad_double_assign_labels ...')
    names = ['x', 'y']
    sols = []
    for idx in range(3):
        sol = make_solution(names, [idx, 0.0])
        sols.append(sol)
    if vrblvl > 0:
        print('the test solutions :')
        for (idx, sol) in enumerate(sols):
            print('Solution', idx+1, ':')
            print(sol)
    set_quad_double_solutions(2, sols, vrblvl-1)
    quad_double_assign_labels(2, vrblvl)
    newsols = get_quad_double_solutions(vrblvl-1)
    if vrblvl > 0:
        print('after assigning labels :')
    fail = 0
    for (idx, sol) in enumerate(newsols):
        if vrblvl > 0:
            print('Solution', idx+1, ':')
            print(sol)
        (_, mval) = endmultiplicity(sol, vrblvl-1)
        if vrblvl > 0:
            print('mval :', mval, mval == idx+1)
        fail = fail + int(mval != idx+1)
    return fail

def test_double_monodromy(vrblvl=0):
    """
    Tests the monodromy breakup in double precision.
    """
    if vrblvl > 0:
        print('in test_double_monodromy ...')
    cyc4 = cyclic(4)
    cyc4e1 = double_embed(4, 1, cyc4, vrblvl=vrblvl-1)
    clear_double_solutions(vrblvl-1)
    c4sols = solve(cyc4e1, vrblvl=vrblvl-1)
    esols = filter_zero_coordinates(c4sols, varname='zz1', tol=1.0e-8, \
        oper='select', vrblvl=vrblvl-1)
    print('the embedded cyclic 4-roots system :')
    for pol in cyc4e1:
        print(pol)
    print('the generic points :')
    for (idx, sol) in enumerate(esols):
        print('Solution', idx+1, ':')
        print(sol)
    fail = int(len(esols) != 4)
    deco = double_monodromy_breakup(cyc4e1, esols, dim=1, \
        verbose=True, nbloops=15, vrblvl=vrblvl)
    if vrblvl > 0:
        print('the decomposition :')
        for (idx, factor) in enumerate(deco):
            print('  factor', idx+1, ':', factor)
    fail = fail + int(len(deco) != 2)
    return fail

def test_double_double_monodromy(vrblvl=0):
    """
    Tests the monodromy breakup in double double precision.
    """
    if vrblvl > 0:
        print('in test_double_double_monodromy ...')
    cyc4 = cyclic(4)
    cyc4e1 = double_double_embed(4, 1, cyc4, vrblvl=vrblvl-1)
    clear_double_double_solutions(vrblvl-1)
    c4sols = solve(cyc4e1, precision='dd', vrblvl=vrblvl-1)
    esols = filter_zero_coordinates(c4sols, varname='zz1', tol=1.0e-8, \
        oper='select', vrblvl=vrblvl-1)
    print('the embedded cyclic 4-roots system :')
    for pol in cyc4e1:
        print(pol)
    print('the generic points :')
    for (idx, sol) in enumerate(esols):
        print('Solution', idx+1, ':')
        print(sol)
    fail = int(len(esols) != 4)
    deco = double_double_monodromy_breakup(cyc4e1, esols, dim=1, \
        verbose=True, nbloops=15, vrblvl=vrblvl)
    if vrblvl > 0:
        print('the decomposition :')
        for (idx, factor) in enumerate(deco):
            print('  factor', idx+1, ':', factor)
    fail = fail + int(len(deco) != 2)
    return fail

def test_quad_double_monodromy(vrblvl=0):
    """
    Tests the monodromy breakup in quad double precision.
    """
    if vrblvl > 0:
        print('in test_quad_double_monodromy ...')
    cyc4 = cyclic(4)
    cyc4e1 = quad_double_embed(4, 1, cyc4, vrblvl=vrblvl-1)
    clear_quad_double_solutions(vrblvl-1)
    c4sols = solve(cyc4e1, precision='qd', vrblvl=vrblvl-1)
    esols = filter_zero_coordinates(c4sols, varname='zz1', tol=1.0e-8, \
        oper='select', vrblvl=vrblvl-1)
    print('the embedded cyclic 4-roots system :')
    for pol in cyc4e1:
        print(pol)
    print('the generic points :')
    for (idx, sol) in enumerate(esols):
        print('Solution', idx+1, ':')
        print(sol)
    fail = int(len(esols) != 4)
    deco = quad_double_monodromy_breakup(cyc4e1, esols, dim=1, \
        verbose=True, nbloops=15, vrblvl=vrblvl)
    if vrblvl > 0:
        print('the decomposition :')
        write_factorization(deco)
    fail = fail + int(len(deco) != 2)
    return fail

def main():
    """
    Runs some tests.
    """
    lvl = 1
    fail = test_double_assign_labels(lvl)
    fail = fail + test_double_double_assign_labels(lvl)
    fail = fail + test_quad_double_assign_labels(lvl)
    fail = fail + test_double_monodromy(lvl)
    fail = fail + test_double_double_monodromy(lvl)
    fail = fail + test_quad_double_monodromy(lvl)
    if fail == 0:
        print('=> All tests passed.')
    else:
        print('Number of failed tests :', fail)

if __name__=='__main__':
    main()
