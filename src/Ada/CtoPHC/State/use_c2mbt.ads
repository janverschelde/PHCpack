with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with C_Integer_Arrays;                  use C_Integer_Arrays;
with C_Double_Arrays;                   use C_Double_Arrays;

function use_c2mbt ( job : integer32;
                     a : C_intarrs.Pointer;
		     b : C_intarrs.Pointer;
                     c : C_dblarrs.Pointer ) return integer32;

-- DESCRIPTION :
--   Provides a gateway from C to the operations in PHCpack
--   to run a homotopy membership test to decide whether a given point
--   belongs to a positive dimensional algebraic set,
--   represented by a witness set.

-- ON ENTRY :
--   job    =  0 : runs the membership test in standard double precision,
--                 for a witness set stored in the standard containers,
--                 on input in a is the verbose flag:
--                 0 for no output, or 1 for diagnostic intermediate output,
--                 on input in b[0] is the dimension n of the test point,
--                 which equals the number of complex coordinates,
--                 and in b[1] is the dimension of the witness set,
--                 on input in c are first the two tolerances:
--                 first on the residual for the evaluation,
--                 second on the tolerance on the membership,
--                 followed then by the coordinates of the test point,
--                 if a test point has dimension n, then c on input
--                 is expected to have 2 + 2*n doubles;
--                 on return in a is the result of the evaluation test:
--                 0 if failure, or 1 if success,
--                 on return in b is the result of the membership test:
--                 0 if failure, or 1 if success;
--   job    =  1 : runs the membership test in double double precision,
--                 for a witness set stored in the dobldobl containers,
--                 on input in a is the verbose flag:
--                 0 for no output, or 1 for diagnostic intermediate output,
--                 on input in b is the dimension n of the test point,
--                 and in b[1] is the dimension of the witness set,
--                 which equals the number of complex coordinates,
--                 on input in c are first the two tolerances:
--                 first on the residual for the evaluation,
--                 second on the tolerance on the membership,
--                 followed then by the coordinates of the test point,
--                 if a test point has dimension n, then c on input
--                 is expected to have 2 + 4*n doubles;
--                 on return in a is the result of the evaluation test:
--                 0 if failure, or 1 if success,
--                 on return in b is the result of the membership test:
--                 0 if failure, or 1 if success;
--   job    =  2 : runs the membership test in quad double precision,
--                 for a witness set stored in the quaddobl containers,
--                 on input in a is the verbose flag:
--                 0 for no output, or 1 for diagnostic intermediate output,
--                 on input in b is the dimension n of the test point,
--                 and in b[1] is the dimension of the witness set,
--                 which equals the number of complex coordinates,
--                 on input in c are first the two tolerances:
--                 first on the residual for the evaluation,
--                 second on the tolerance on the membership,
--                 followed then by the coordinates of the test point,
--                 if a test point has dimension n, then c on input
--                 is expected to have 2 + 8*n doubles;
--                 on return in a is the result of the evaluation test:
--                 0 if failure, or 1 if success,
--                 on return in b is the result of the membership test:
--                 0 if failure, or 1 if success.

-- ON RETURN :
--   0 if the operation was successful, otherwise something went wrong,
--   e.g.: job not in the right range.
