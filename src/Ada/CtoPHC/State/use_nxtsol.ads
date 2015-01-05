with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with C_Integer_Arrays;                  use C_Integer_Arrays;
with C_Double_Arrays;                   use C_Double_Arrays;

function use_nxtsol ( job : integer32;
                      a : C_intarrs.Pointer;
                      b : C_intarrs.Pointer;
                      c : C_dblarrs.Pointer ) return integer32;

-- DESCRIPTION :
--   Provides a gateway to the path trackers in PHCpack with generators.

-- ON ENTRY :
--   job    =   0 : initialize homotopy in standard double precision
--                  with the target and start systems stored in containers;
--          =   1 : initialize homotopy in double double precision
--                  with the target and start systems stored in containers;
--          =   2 : initialize homotopy in quad double precision
--                  with the target and start systems stored in containers;
--          =   3 : takes solution at position a[0] in the standard solution
--                  container and initializes the standard path tracker;
--          =   4 : takes solution at position a[0] in the double double
--                  container and initializes the double double path tracker;
--          =   5 : takes solution at position a[0] in the quad double
--                  container and initializes the quad double path tracker;
--          =   6 : applies one predictor-corrector step in standard double
--                  precision and replaces the solution in the standard
--                  solutions container at position equal to the value of a[0];
--          =   7 : applies one predictor-corrector step in double double
--                  precision and replaces the solution in the double double
--                  solutions container at position equal to the value of a[0];
--          =   8 : applies one predictor-corrector step in quad double
--                  precision and replaces the solution in the quad double
--                  solutions container at position equal to the value of a[0];
--         =    9 : clears standard path tracker;
--         =   10 : clears double double path tracker;
--         =   11 : clears quad double path tracker;
--         =   12 : initialize homotopy for arbitrary multiprecision,
--                  with the target and start systems stored in container,
--                  using the number of decimal places in the working
--                  precision as given by the value in a[0];
--         =   13 : takes solution at position a[0] in the multiprecision 
--                  solution container and initializes the path tracker.
--         =   14 : applies one predictor-corrector step in multiple
--                  precision and replaces the solution in the multiprecision
--                  solutions container at position equal to the value of a[0];
--         =   15 : clears the multiprecision path tracker.
--         =   16 : initialize homotopy for variable precision tracking,
--                  with the target and start systems stored as strings,
--                  the function expects three values in a:
--                  a[0] = whether a fixed gamma is used or not (1 or 0),
--                  a[1] = the total number of characters in the string b,
--                  a[2] = the start of the second (start) system in b;
--                  b holds the string representations of two systems,
--                  respectively the target and start system in the homotopy;
--         =   17 : takes solution given as string representation in b,
--                  in a[0] = the number of characters in the string b,
--                  in a[1] = the number of variables in the solution
--                  represented by b, to initialize the variable precision
--                  path tracker with;
--         =   18 : applies one predictor-corrector step in variable
--                  precision and expects four values in a on entry:
--                  a[0] = the wanted number of accurate decimal places;
--                  a[1] = the maximum precision to be used;
--                  a[2] = the maximum number of corrector steps;
--                  a[3] = whether intermediate output is wanted or not,
--                  where 1 is true and 0 is false; and
--                  returns in a[0] the number of characters in the string
--                  of the solution and returns in b the solution string;
--         =   19 : clears the variable precision path tracker.
--
-- ON RETURN :
--   0 if the operation was successful, otherwise something went wrong,
--   e.g.: job not in the right range.
