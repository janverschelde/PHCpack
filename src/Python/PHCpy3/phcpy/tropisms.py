"""
The module tropisms exports function to manage numerically computed
tropisms in double, double double, or quad double precision.
"""

def standard_initialize(nbt, dim, wnd, dir, err):
    r"""
    Initializes the direction vectors computed in double precision,
    along with estimates for their winding numbers and errors.
    On entry are the following five parameters:

    *nbt*: the number of direction vectors;

    *dim*: the number of coordinates in each vector;

    *wnd*: a list of integer values for the winding numbers, as many as *nbt*;

    *dir*: a list of lists of doubles with the coordinates of the directions,
    each inner list has *dim* doubles and *nbt* vectors are given;

    *err*: a list of *nbt* doubles.
    """
    from phcpy.phcpy2c3 import py2c_numbtrop_standard_initialize as store
    flat = []
    for vec in dir:
        flat = flat + vec
    data = wnd + flat + err
    store(nbt, dim, str(data))

def dobldobl_initialize(nbt, dim, wnd, dir, err):
    r"""
    Initializes the direction vectors computed in double double precision,
    along with estimates for their winding numbers and errors.
    On entry are the following five parameters:

    *nbt*: the number of direction vectors;

    *dim*: the number of coordinates in each vector;

    *wnd*: a list of integer values for the winding numbers, as many as *nbt*;

    *dir*: a list of lists of double doubles with the coordinates of the 
    directions, each inner list has *dim* double doubles and *nbt* vectors;

    *err*: a list of *nbt* double doubles.
    """
    from phcpy.phcpy2c3 import py2c_numbtrop_dobldobl_initialize as store
    flat = []
    for vec in dir:
        flat = flat + vec
    data = wnd + flat + err
    store(nbt, dim, str(data))

def quaddobl_initialize(nbt, dim, wnd, dir, err):
    r"""
    Initializes the direction vectors computed in quad double precision,
    along with estimates for their winding numbers and errors.
    On entry are the following five parameters:

    *nbt*: the number of direction vectors;

    *dim*: the number of coordinates in each vector;

    *wnd*: a list of integer values for the winding numbers, as many as *nbt*;

    *dir*: a list of lists of quad doubles with the coordinates of the 
    directions, each inner list has *dim* quad doubles and *nbt* vectors;

    *err*: a list of *nbt* double doubles.
    """
    from phcpy.phcpy2c3 import py2c_numbtrop_quaddobl_initialize as store
    flat = []
    for vec in dir:
        flat = flat + vec
    data = wnd + flat + err
    store(nbt, dim, str(data))

def standard_size():
    """
    Returns the number of tropisms stored in standard double precision.
    """
    from phcpy.phcpy2c3 import py2c_numbtrop_standard_size as get_size
    return get_size()

def dobldobl_size():
    """
    Returns the number of tropisms stored in double double precision.
    """
    from phcpy.phcpy2c3 import py2c_numbtrop_dobldobl_size as get_size
    return get_size()

def quaddobl_size():
    """
    Returns the number of tropisms stored in quad double precision.
    """
    from phcpy.phcpy2c3 import py2c_numbtrop_quaddobl_size as get_size
    return get_size()

def standard_dimension():
    """
    Returns the dimension of the tropisms stored in standard double precision.
    """
    from phcpy.phcpy2c3 import py2c_numbtrop_standard_dimension as get_dim
    return get_dim()

def dobldobl_dimension():
    """
    Returns the dimension of the tropisms stored in double double precision.
    """
    from phcpy.phcpy2c3 import py2c_numbtrop_dobldobl_dimension as get_dim
    return get_dim()

def quaddobl_dimension():
    """
    Returns the dimension of the tropisms stored in quad double precision.
    """
    from phcpy.phcpy2c3 import py2c_numbtrop_quaddobl_dimension as get_dim
    return get_dim()

def standard_retrieve(nbt, dim):
    r"""
    Given on input the number of tropisms in *nbt* and the dimension in *dim*,
    returns a tuple of three lists: the winding numbers, coordinates of
    the direction vectors, and the errrors; in standard double precision.
    """
    from ast import literal_eval
    from phcpy.phcpy2c3 import py2c_numbtrop_standard_retrieve as load
    (fail, strdata) = load(nbt, dim)
    data = literal_eval(strdata)
    wnd = [int(data[k]) for k in range(nbt)]
    dirs = []
    for i in range(nbt):
        dirs.append([data[nbt+i*dim+j] for j in range(dim)])
    err = [data[nbt*(dim+1)+k] for k in range(nbt)]
    return (wnd, dirs, err)

def dobldobl_retrieve(nbt, dim):
    r"""
    Given on input the number of tropisms in *nbt* and the dimension in *dim*,
    returns a tuple of three lists: the winding numbers, coordinates of
    the direction vectors, and the errrors; in double double precision.
    """
    from ast import literal_eval
    from phcpy.phcpy2c3 import py2c_numbtrop_dobldobl_retrieve as load
    (fail, strdata) = load(nbt, dim)
    data = literal_eval(strdata)
    wnd = [int(data[k]) for k in range(nbt)]
    dirs = []
    for i in range(nbt):
        dirs.append([data[nbt+i*2*dim+j] for j in range(2*dim)])
    err = [data[nbt*(2*dim+1)+k] for k in range(2*nbt)]
    return (wnd, dirs, err)

def quaddobl_retrieve(nbt, dim):
    r"""
    Given on input the number of tropisms in *nbt* and the dimension in *dim*,
    returns a tuple of three lists: the winding numbers, coordinates of
    the direction vectors, and the errrors; in quad double precision.
    """
    from ast import literal_eval
    from phcpy.phcpy2c3 import py2c_numbtrop_quaddobl_retrieve as load
    (fail, strdata) = load(nbt, dim)
    data = literal_eval(strdata)
    wnd = [int(data[k]) for k in range(nbt)]
    dirs = []
    for i in range(nbt):
        dirs.append([data[nbt+i*4*dim+j] for j in range(4*dim)])
    err = [data[nbt*(4*dim+1)+k] for k in range(4*nbt)]
    return (wnd, dirs, err)

def retrieve_standard_tropism(dim, idx):
    r"""
    Returns the winding number, coordinates of the direction, and its error,
    stored in double precision, of dimension *dim*, and index *idx*.
    The index must be in the range 1..standard_size().
    Observe that the index counter starts at one and not at zero.
    """
    from ast import literal_eval
    from phcpy.phcpy2c3 \
    import py2c_numbtrop_standard_retrieve_tropism as load
    (fail, wnd, strdata) = load(dim, idx)
    data = literal_eval(strdata)
    dir = [data[k] for k in range(dim)]
    err = data[dim]
    return (wnd, dir, err)

def retrieve_dobldobl_tropism(dim, idx):
    r"""
    Returns the winding number, coordinates of the direction, and its error,
    stored in double double precision, of dimension *dim*, and index *idx*.
    The index must be in the range 1..dobldobl_size().
    Observe that the index counter starts at one and not at zero.
    """
    from ast import literal_eval
    from phcpy.phcpy2c3 \
    import py2c_numbtrop_dobldobl_retrieve_tropism as load
    (fail, wnd, strdata) = load(dim, idx)
    data = literal_eval(strdata)
    dir = [data[k] for k in range(2*dim)]
    err = [data[2*dim], data[2*dim+1]]
    return (wnd, dir, err)

def retrieve_quaddobl_tropism(dim, idx):
    r"""
    Returns the winding number, coordinates of the direction, and its error,
    stored in quad double precision, of dimension *dim*, and index *idx*.
    The index must be in the range 1..quaddobl_size().
    Observe that the index counter starts at one and not at zero.
    """
    from ast import literal_eval
    from phcpy.phcpy2c3 \
    import py2c_numbtrop_quaddobl_retrieve_tropism as load
    (fail, wnd, strdata) = load(dim, idx)
    data = literal_eval(strdata)
    dir = [data[k] for k in range(4*dim)]
    err = [data[4*dim], data[4*dim+1], data[4*dim+2], data[4*dim+3]]
    return (wnd, dir, err)

def store_standard_tropism(dim, idx, wnd, dir, err):
    r"""
    Stores the tropism, given in standard double precision,
    with dim doubles as coordinates in the list *dir*, the error in *err*,
    and the winding number *wnd*, at position *idx*.
    The index *idx* must be in the range 1..standard_size().
    """
    from phcpy.phcpy2c3 import py2c_numbtrop_store_standard_tropism as store
    data = [x for x in dir]
    data.append(err)
    strdata = str(data)
    store(dim, idx, wnd, strdata)

def store_dobldobl_tropism(dim, idx, wnd, dir, err):
    r"""
    Stores the tropism, given in double double precision,
    with dim doubles as coordinates in the list *dir*, the error in *err*,
    and the winding number *wnd*, at position *idx*.
    The index *idx* must be in the range 1..dobldobl_size().
    """
    from phcpy.phcpy2c3 import py2c_numbtrop_store_dobldobl_tropism as store
    data = [x for x in dir]
    data.append(err[0])
    data.append(err[1])
    strdata = str(data)
    store(dim, idx, wnd, strdata)

def store_quaddobl_tropism(dim, idx, wnd, dir, err):
    r"""
    Stores the tropism, given in quad double precision,
    with dim doubles as coordinates in the list *dir*, the error in *err*,
    and the winding number *wnd*, at position *idx*.
    The index *idx* must be in the range 1..quaddobl_size().
    """
    from phcpy.phcpy2c3 import py2c_numbtrop_store_quaddobl_tropism as store
    data = [x for x in dir]
    data.append(err[0])
    data.append(err[1])
    data.append(err[2])
    data.append(err[3])
    strdata = str(data)
    store(dim, idx, wnd, strdata)

def standard_clear():
    """
    Clears the tropisms stored in standard double precision.
    """
    from phcpy.phcpy2c3 import py2c_numbtrop_standard_clear
    py2c_numbtrop_standard_clear()

def dobldobl_clear():
    """
    Clears the tropisms stored in double double precision.
    """
    from phcpy.phcpy2c3 import py2c_numbtrop_dobldobl_clear
    py2c_numbtrop_dobldobl_clear()

def quaddobl_clear():
    """
    Clears the tropisms stored in quad double precision.
    """
    from phcpy.phcpy2c3 import py2c_numbtrop_quaddobl_clear
    py2c_numbtrop_quaddobl_clear()

def test_standard_store_load():
    """
    Tests the storing and loading of numerically computed tropisms,
    in standard double precision.
    """
    print('testing storing and loading numerical tropisms...')
    nbr = int(input('Give the number of tropisms : '))
    dim = int(input('Give the dimension : '))
    from random import randint
    wnd = [randint(10, 99) for _ in range(nbr)]
    from random import uniform as u
    dirs = []
    for _ in range(nbr):
        dirs.append([u(-1, +1) for _ in range(dim)])
    errs = [u(0,1)/1000.0 for _ in range(nbr)]
    print('random winding numbers :', wnd)
    for k in range(len(dirs)):
        print('direction', k+1, ':', dirs[k])
    print('random errors :', errs)
    standard_initialize(nbr, dim, wnd, dirs, errs)
    (retwnd, retdirs, reterrs) = standard_retrieve(nbr, dim)
    retsize = standard_size()
    retdim = standard_dimension()
    print('retrieved number of tropisms :', retsize)
    print('retrieved dimension of the tropisms :', retdim)
    print('retrieved winding numbers :', retwnd)
    print('retrieved directions :')
    for k in range(len(retdirs)):
        print('direction', k+1, ':', retdirs[k])
    print('retrieved errors :', reterrs)
    while True:
        idx = int(input('Give an index (0 to exit): '))
        if idx < 1:
            break
        (idxwnd, idxdir, idxerr) = retrieve_standard_tropism(dim, idx)
        print('-> winding number for tropism', idx, ':', idxwnd)
        print('-> tropism', idx, 'has coordinates :', idxdir)
        print('-> the error :', idxerr)
        ranwnd = randint(10, 99)
        randir = [u(-1, 1) for _ in range(dim)]
        ranerr = u(0, 1)/1000.0
        print('store tropism %d with winding number %d' % (idx,ranwnd))
        print('with new coordinates :', randir)
        print('and error', ranerr);
        store_standard_tropism(dim, idx, ranwnd, randir, ranerr)
    standard_clear()

def test_dobldobl_store_load():
    """
    Tests the storing and loading of numerically computed tropisms,
    in double double precision.
    """
    print('testing storing and loading numerical tropisms...')
    nbr = int(input('Give the number of tropisms : '))
    dim = int(input('Give the dimension : '))
    from random import randint
    wnd = [randint(10, 99) for _ in range(nbr)]
    from random import uniform as u
    dirs = []
    for _ in range(nbr):
        dirs.append([u(-1, +1) for _ in range(2*dim)])
    errs = [u(0,1)/1000.0 for _ in range(2*nbr)]
    print('random winding numbers :', wnd)
    for k in range(len(dirs)):
        print('direction', k+1, ':', dirs[k])
    print('random errors :', errs)
    dobldobl_initialize(nbr, dim, wnd, dirs, errs)
    (retwnd, retdirs, reterrs) = dobldobl_retrieve(nbr, dim)
    retsize = dobldobl_size()
    retdim = dobldobl_dimension()
    print('retrieved number of tropisms :', retsize)
    print('retrieved dimension of the tropisms :', retdim)
    print('retrieved winding numbers :', retwnd)
    print('retrieved directions :')
    for k in range(len(retdirs)):
        print('direction', k+1, ':', retdirs[k])
    print('retrieved errors :', reterrs)
    while True:
        idx = int(input('Give an index (0 to exit): '))
        if idx < 1:
            break
        (idxwnd, idxdir, idxerr) = retrieve_dobldobl_tropism(dim, idx)
        print('-> winding number for tropism', idx, ':', idxwnd)
        print('-> tropism', idx, 'has coordinates :', idxdir)
        print('-> the error :', idxerr)
        ranwnd = randint(10, 99)
        randir = [u(-1, 1) for _ in range(2*dim)]
        ranerr = [u(0, 1)/1000.0, u(0, 1)/1000.0]
        print('store tropism %d with winding number %d' % (idx,ranwnd))
        print('with new coordinates :', randir)
        print('and error', ranerr);
        store_dobldobl_tropism(dim, idx, ranwnd, randir, ranerr)
    dobldobl_clear()

def test_quaddobl_store_load():
    """
    Tests the storing and loading of numerically computed tropisms,
    in quad double precision.
    """
    print('testing storing and loading numerical tropisms...')
    nbr = int(input('Give the number of tropisms : '))
    dim = int(input('Give the dimension : '))
    from random import randint
    wnd = [randint(10, 99) for _ in range(nbr)]
    from random import uniform as u
    dirs = []
    for _ in range(nbr):
        dirs.append([u(-1, +1) for _ in range(4*dim)])
    errs = [u(0,1)/1000.0 for _ in range(4*nbr)]
    print('random winding numbers :', wnd)
    for k in range(len(dirs)):
        print('direction', k+1, ':', dirs[k])
    print('random errors :', errs)
    quaddobl_initialize(nbr, dim, wnd, dirs, errs)
    (retwnd, retdirs, reterrs) = quaddobl_retrieve(nbr, dim)
    retsize = quaddobl_size()
    retdim = quaddobl_dimension()
    print('retrieved number of tropisms :', retsize)
    print('retrieved dimension of the tropisms :', retdim)
    print('retrieved winding numbers :', retwnd)
    print('retrieved directions :')
    for k in range(len(retdirs)):
        print('direction', k+1, ':', retdirs[k])
    print('retrieved errors :', reterrs)
    while True:
        idx = int(input('Give an index (0 to exit): '))
        if idx < 1:
            break
        (idxwnd, idxdir, idxerr) = retrieve_quaddobl_tropism(dim, idx)
        print('-> winding number for tropism', idx, ':', idxwnd)
        print('-> tropism', idx, 'has coordinates :', idxdir)
        print('-> the error :', idxerr)
        ranwnd = randint(10, 99)
        randir = [u(-1, 1) for _ in range(4*dim)]
        ranerr = [u(0, 1)/1000.0, u(0, 1)/1000.0]
        print('store tropism %d with winding number %d' % (idx,ranwnd))
        print('with new coordinates :', randir)
        print('and error', ranerr);
        store_quaddobl_tropism(dim, idx, ranwnd, randir, ranerr)
    quaddobl_clear()

def management_test():
    """
    Tests the management of numerically computed tropisms.
    """
    print('testing management of numerical tropisms...')
    prc = input('give the precision (d, dd, or qd) : ')
    if prc == 'd':
        test_standard_store_load()
    elif prc == 'dd':
        test_dobldobl_store_load()
    else:
        test_quaddobl_store_load()

def test(precision='d'):
    """
    Tests the numerical computation of a tropism.
    """
    pols = ['x + x*y + y^3;', 'x;']
    from phcpy.solver import solve
    start = ['x + (-5.89219623474258E-02 - 9.98262591883082E-01*i)*y^2' \
         + ' + ( 5.15275165825429E-01 + 8.57024797472965E-01*i)*x*y;', \
        'x+(-9.56176913648087E-01 - 2.92789531586119E-01*i);']
    startsols = solve(start, verbose=False)
    print('computed', len(startsols), 'start solutions')
    from phcpy.tuning import order_endgame_extrapolator_set as set
    set(4)
    if precision == 'd':
        from phcpy.trackers import standard_double_track as track
    elif precision == 'dd':
        from phcpy.trackers import double_double_track as track
    elif precision == 'qd':
        from phcpy.trackers import quad_double_track as track
    else:
        print('wrong level of precision')
        return
    sols = track(pols, start, startsols)
    print('the solutions at the end :')
    for sol in sols:
        print(sol)
    if precision == 'd':
        print('size of the tropisms container :', standard_size())
        (wnd, dirs, errs) = standard_retrieve(len(sols), len(pols))
    elif precision == 'dd':
        print('size of the tropisms container :', dobldobl_size())
        (wnd, dirs, errs) = dobldobl_retrieve(len(sols), len(pols))
    else:
        print('size of the tropisms container :', quaddobl_size())
        (wnd, dirs, errs) = quaddobl_retrieve(len(sols), len(pols))
    print('the winding numbers :', wnd)
    print('the directions :')
    for dir in dirs:
        print(dir)
    print('the errors :', errs)
    # check which tuning parameters were used
    # from phcpy.tuning import tune_track_parameters as tune
    # tune(interactive=True)

if __name__ == "__main__":
    test('d')
