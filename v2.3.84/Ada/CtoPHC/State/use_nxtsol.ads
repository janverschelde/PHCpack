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
--                  solutions container at position equal to the value of a[0].
--
-- ON RETURN :
--   0 if the operation was successful, otherwise something went wrong,
--   e.g.: job not in the right range.
