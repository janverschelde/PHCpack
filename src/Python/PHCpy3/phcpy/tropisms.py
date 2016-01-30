"""
The module tropisms exports function to manage numerically computed
tropisms in double, double double, or quad double precision.
"""

def standard_initialize(nbt, dim, wnd, dir, err):
    """
    Initializes the direction vectors computed in double precision,
    along with estimates for their winding numbers and errors.
    On entry are the following five parameters:
    nbt : the number of direction vectors;
    dim : the number of coordinates in each vector;
    wnd : a list of integer values for the winding numbers, as many as nbt;
    dir : a list of lists of doubles with the coordinates of the directions,
    each inner list has dim doubles and nbt vectors are given;
    err : a list of nbt doubles.
    """
    from phcpy.phcpy2c3 import py2c_numbtrop_standard_initialize as store
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

def standard_dimension():
    """
    Returns the dimension of the tropisms stored in standard double precision.
    """
    from phcpy.phcpy2c3 import py2c_numbtrop_standard_dimension as get_dim
    return get_dim()

def standard_retrieve(nbt, dim):
    """
    Given on input the number of tropisms in nbt and the dimension in dim,
    returns a tuple of three lists: the winding numbers, coordinates of
    the direction vectors, and the errrors.
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

def standard_retrieve_tropism(dim, idx):
    """
    Returns the winding number, coordinates of the direction, and its error,
    stored in double precision, of dimensin dim, and index idx.
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

def test_standard_store_load():
    """
    Tests the storing and loading of numerically computed tropisms,
    in standard double precision.
    """
    print('testing storing and loading numerical tropisms...')
    nbr = int(input('Give the number of tropisms : '))
    dim = int(input('Give the dimension : '))
    from random import randint
    wnd = [randint(1, 10) for _ in range(nbr)]
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
       (idxwnd, idxdir, idxerr) = standard_retrieve_tropism(dim, idx)
       print('winding number for tropism', idx, ':', idxwnd)
       print('tropism', idx, 'has coordinates :', idxdir)
       print('the error :', idx)

def test():
    """
    Tests the management of numerically computed tropisms.
    """
    print('testing management of numerical tropisms...')
    test_standard_store_load()

if __name__ == "__main__":
    test()
