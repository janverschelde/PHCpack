with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with C_Integer_Arrays;                  use C_Integer_Arrays;
with C_Double_Arrays;                   use C_Double_Arrays;

function use_sweep ( job : integer32;
                     a : C_intarrs.Pointer;
                     b : C_intarrs.Pointer;
                     c : C_dblarrs.Pointer;
                     vrblvl : integer32 := 0 ) return integer32;

-- DESCRIPTION :
--   Provides a gateway to the sweep homotopy and the procedures that
--   track solution paths defined by the sweep homotopy.
--   A sweep homotopy makes a convex combination between sequences of
--   values for the parameters.

-- ON ENTRY :
--   job    =   0 : defines the indices of the variables that serve
--                  as parameters in the sweep homotopy numerically, given 
--                  in a[0] the number of equations,
--                  in a[1] the total number of variables,
--                  in a[2] the number m of parameters,
--                  in b a list of m integers with indices
--                  in the range 1..n defining the parameters;
--          =   1 : defines the parameters in the sweep homotopy symbolically,
--                  given in a[0] the number of equations,
--                  in a[1] the total number of variables,
--                  in a[2] the number m of parameters, and
--                  in a[3] the number of characters stored in b,
--                  the symbols in b are separated by one space,
--                  and the symbol table must be initialized properly
--                  so it contains all the symbols listed in b;
--          =   2 : returns in a the number of equations;
--          =   3 : returns in a the number of variables;
--          =   4 : returns in a the number of parameters;
--          =   5 : returns the indices of the parameters in a;
--          =   6 : returns in b the string representations of the
--                  parameters in the sweep homotopy, given on input in a
--                  an upper bound on the length of the string to store
--                  the sequence of symbols (separated by one space);
--          =   7 : clears the definition of the parameters in the sweep.
--          =   8 : sets the start or target values of the parameters,
--                  in standard double, double double, or quad double
--                  precision, on input are
--                  in a[0] : the precision level, 0, 1, or 2 for standard
--                  double, double double, or quad double respectively,
--                  in a[1] : 0 or 1, for start or target respectively,
--                  in b[0] : the number of coefficients to store the values
--                  real and imaginary parts of the start or target,
--                  irrespective the precision (d, dd, or qd), and
--                  in c : the values of the start or target parameters;
--          =   9 : gets the start or target values of the parameters,
--                  in standard double, double double, or quad double
--                  precision, on input are
--                  in a[0] : the precision level, 0, 1, or 2 for standard
--                  double, double double, or quad double respectively,
--                  in a[1] : 0 or 1, for start or target respectively,
--                  in b[0] : the number of coefficients to store the values
--                  real and imaginary parts of the start or target,
--                  irrespective the precision (d, dd, or qd), and
--                  on return in c are the values of the start or target;
--          =  10 : runs a complex convex-parameter sweep between the
--                  defined values of the parameter for start and target,
--                  with on input 
--                  in a[0] : the precision level, 0, 1, or 2 for standard
--                  double, double double, or quad double respectively,
--                  in a[1] : 0, 1, or 2 for the choice of the gamma,
--                  0: randomly generated gamma, should be the default,
--                  1: no gamma constant, gamma = 1 (which could fail),
--                  2: user given gamma, with real and complex part in
--                  c[0] and c[1] respectively;
--                  for this to work, the proper systems and solutions
--                  containers must be initialized, on return, the new
--                  solutions are stored in the solutions container;
--          =  11 : runs a real natural parameter sweep between the
--                  defined values of the parameter for start and target,
--                  with on input 
--                  in a : the precision level, 0, 1, or 2 for standard
--                  double, double double, or quad double respectively,
--                  for this to work, the proper systems and solutions
--                  containers must be initialized, on return, the new
--                  solutions are stored in the solutions container.

-- ON RETURN :
--   0 if all went well.
