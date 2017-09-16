with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with C_Integer_Arrays;                  use C_Integer_Arrays;
with C_Double_Arrays;                   use C_Double_Arrays;

function use_reduction ( job : integer32;
                         a : C_intarrs.Pointer;
                         b : C_intarrs.Pointer;
                         c : C_dblarrs.Pointer ) return integer32;

-- DESCRIPTION :
--   Gateway to the reduction methods on the coefficient matrix.

-- ON ENTRY :
--   job =     1 : reduces the system in the standard container,
--                 applying row reduction on the coefficient matrix:
--                 on input in a is either 0 or 1:
--                 0 : triangulation of the coefficient matrix,
--                 1 : makes the coefficient matrix diagonal;
--                 the system in the container is replaced by the
--                 reduced system.

-- ON RETURN :
--   0 if the operation was successful, otherwise something went wrong,
--   or job not in the right range.
