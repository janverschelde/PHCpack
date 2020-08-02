with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with C_Integer_Arrays;                  use C_Integer_Arrays;
with C_Double_Arrays;                   use C_Double_Arrays;

function use_reduction ( job : integer32;
                         a : C_intarrs.Pointer;
                         b : C_intarrs.Pointer;
                         c : C_dblarrs.Pointer;
                         vrblvl : integer32 := 0 ) return integer32;

-- DESCRIPTION :
--   Gateway to the reduction methods on the coefficient matrix.

-- ON ENTRY :
--   job =     1 : reduces the system in the standard container,
--                 applying row reduction in standard double precision,
--                 on the coefficient matrix:
--                 on input in a is either 0 or 1:
--                 0 : triangulation of the coefficient matrix,
--                 1 : makes the coefficient matrix diagonal;
--                 the system in the standard container is replaced
--                 by the reduced system;
--       =     2 : reduces the system in the dobldobl container,
--                 applying row reduction in double double precision,
--                 on the coefficient matrix:
--                 on input in a is either 0 or 1:
--                 0 : triangulation of the coefficient matrix,
--                 1 : makes the coefficient matrix diagonal;
--                 the system in the dobldobl container is replaced
--                 by the reduced system;
--       =     3 : reduces the system in the quaddobl container,
--                 applying row reduction in quad double precision,
--                 on the coefficient matrix:
--                 on input in a is either 0 or 1:
--                 0 : triangulation of the coefficient matrix,
--                 1 : makes the coefficient matrix diagonal;
--                 the system in the quaddobl container is replaced
--                 by the reduced system.
--       =     4 : applies nonlinear reduction to the systems in the
--                 standard double container, on input in a are three limits:
--                 (1) the maximum number of equal degree replacements,
--                 (2) the maximum number of computed S-polynomials, and
--                 (3) the maximum number of computed R-polynomials;
--                 on return in b are three counters:
--                 (1) the number of equal degree replacements,
--                 (2) the number of computed S-polynomials, and
--                 (3) the number of computed R-polynomials;
--                 the system in the standard container is replaced
--                 by the reduced system.

-- ON RETURN :
--   0 if the operation was successful, otherwise something went wrong,
--   or job not in the right range.
