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
    print('retrieved winding numbers :', retwnd)
    print('retrieved directions :')
    for k in range(len(retdirs)):
        print('direction', k+1, ':', retdirs[k])
    print('retrieved errors :', reterrs)

def test():
    """
    Tests the management of numerically computed tropisms.
    """
    print('testing management of numerical tropisms...')
    test_standard_store_load()

if __name__ == "__main__":
    test()
