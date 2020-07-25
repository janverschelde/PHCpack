with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with C_Integer_Arrays;                  use C_Integer_Arrays;
with C_Double_Arrays;                   use C_Double_Arrays;

function use_mapcon ( job : integer32;
                      a : C_intarrs.Pointer;
		      b : C_intarrs.Pointer;
                      c : C_dblarrs.Pointer;
                      vrblvl : integer32 := 0 ) return integer32;

-- DESCRIPTION :
--   Provides a gateway to the container for monomial maps.

-- ON ENTRY :
--   job    =   0 : solves binomial system in Laurent system container,
--                  the option puretopdim is given in a[0];
--          =   1 : writes the maps stored in the container to screen.
--          =   2 : clears the maps stored in the container;
--          =   3 : returns in a the top dimension of the maps;
--          =   4 : returns in b the number of maps of dimension given in a;
--          =   5 : given in a[0] the dimension and in a[1] the index of
--                  a map, returns in b the degree of the map.
--          =   6 : given in a[0] the dimension of the map,
--                        in a[1] the index of of the map
--                        in a[2] the number n of variables,
--                  returns in c the coefficients of the map,
--                  as 2*n doubles, representing real and imaginary parts;
--          =   7 : given in a[0] the dimension of the map,
--                        in a[1] the index of of the map
--                        in a[2] the number n of variables,
--                  returns in b the exponents of the map, as 
--                  a flattened vector of dim vectors of range 1..n;
--          =   8 : given in a[0] the dimension of the map,
--                        in a[1] the index of of the map
--                        in a[2] the number n of variables,
--                  returns in b the exponents of the map, as 
--                  a flattened vector of dim vectors of range 1..n, and
--                  returns in c the coefficients of the map,
--                  as 2*n doubles, representing real and imaginary parts.

-- ON RETURN :
--   0 if the operation was successful, otherwise something went wrong,
--   e.g.: indices to monomial out of range, or job not in the proper range.
