with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with C_Integer_Arrays;                  use C_Integer_Arrays;
with C_Double_Arrays;                   use C_Double_Arrays;

function use_multip ( job : integer32;
                      a : C_intarrs.Pointer;
                      b : C_intarrs.Pointer;
                      c : C_dblarrs.Pointer;
                      vrblvl : integer32 := 0 ) return integer32;

-- DESCRIPTION :
--   Exports the computation of the multiplicity structure
--   for an isolated solution of a polynomial system.
--   Both the system and the solution are expected in the
--   systems and solutions container prior to calling use_multip.
--   Three different levels of precision are supported:
--   standard double, double double, and quad double.

-- ON ENTRY :
--   job    =   0 : computes in standard double precision,
--                  on entry in a[0] is the maximum deflation order
--                  in b[0] is 0 or 1, as a flag to indicate whether
--                  to be silent or verbose (to screen) respectively
--                  (note that b must hold space for order+1 integers),
--                  in c[0] is the tolerance for the numerical rank,
--                  on return in a is the computed multiplicity
--                  and in b are order+1 integers with the values
--                  of the Hilbert function;
--   job    =   1 : computes in double double precision,
--                  on entry in a[0] is the maximum deflation order
--                  in b[0] is 0 or 1, as a flag to indicate whether
--                  to be silent or verbose (to screen) respectively
--                  (note that b must hold space for order+1 integers),
--                  in c[0] is the tolerance for the numerical rank,
--                  on return in a is the computed multiplicity
--                  and in b are order+1 integers with the values
--                  of the Hilbert function;
--   job    =   2 : computes in quad double precision,
--                  on entry in a[0] is the maximum deflation order
--                  in b[0] is 0 or 1, as a flag to indicate whether
--                  to be silent or verbose (to screen) respectively
--                  (note that b must hold space for order+1 integers),
--                  in c[0] is the tolerance for the numerical rank,
--                  on return in a is the computed multiplicity
--                  and in b are order+1 integers with the values
--                  of the Hilbert function.

-- ON RETURN :
--   0 if all went well.
