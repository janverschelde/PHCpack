with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with C_Integer_Arrays;                  use C_Integer_Arrays;
with C_Double_Arrays;                   use C_Double_Arrays;

function use_avvcon ( job : integer32;
                       a : C_intarrs.Pointer;
                       b : C_intarrs.Pointer;
                       c : C_dblarrs.Pointer;
                       vrblvl : integer32 := 0 ) return integer32;

-- DESCRIPTION :
--   Provides a gateway to the container for arrays of vectors of vectors.

-- ON ENTRY :
--   job    =   0 : initialize the arrays of vectors of real vectors
--                  with the dimension in a[0] and the number of vectors
--                  in b[k] for the k-th array;
--   job    =   1 : initialize the arrays of vectors of complex vectors
--                  with the dimension in a[0] and the number of vectors
--                  in b[k] for the k-th array;
--   job    =   2 : sets the j-th vector of the i-th array to the a[0]
--                  real numbers in c given on input,
--                  where i is given by b[0] and j is given by b[1];
--   job    =   3 : sets the j-th vector of the i-th array to the a[0]
--                  complex numbers in c given on input,
--                  where i is given by b[0] and j is given by b[1];
--   job    =   4 : returns in a[0] the number of double vecvecs arrays;
--   job    =   5 : returns in a[0] the number of complex vecvecs arrays;
--   job    =   6 : given in a[0] the index of a double vecvecs array,
--                  returns in b[0] the size of that array;
--   job    =   7 : given in a[0] the index of a complex vecvecs array,
--                  returns in b[0] the size of that array;
--   job    =   8 : given in a[0] the index of a double vecvecs array,
--                  and in a[1] the index of a vector,
--                  returns in b[0] the size of that vector;
--   job    =   9 : given in a[0] the index of a complex vecvecs array,
--                  and in a[1] the index of a vector,
--                  returns in b[0] the size of that vector.
--   job    =  10 : gets the j-th vector of the i-th array to the a[0]
--                  real numbers in c returned as output,
--                  where i is given by b[0] and j is given by b[1];
--   job    =  11 : gets the j-th vector of the i-th array to the a[0]
--                  complex numbers in c returned as output,
--                  where i is given by b[0] and j is given by b[1];
--   job    =  12 : clears the arrays of vectors of real vectors;
--   job    =  13 : clears the arrays of vectors of complex vectors.
