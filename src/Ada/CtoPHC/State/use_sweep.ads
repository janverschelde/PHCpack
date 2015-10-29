with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with C_Integer_Arrays;                  use C_Integer_Arrays;
with C_Double_Arrays;                   use C_Double_Arrays;

function use_sweep ( job : integer32;
                     a : C_intarrs.Pointer;
                     b : C_intarrs.Pointer;
                     c : C_dblarrs.Pointer ) return integer32;

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
