with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with C_Integer_Arrays;                  use C_Integer_Arrays;
with C_Double_Arrays;                   use C_Double_Arrays;

function use_numbtrop ( job : integer32;
                        a : C_intarrs.Pointer;
                        b : C_intarrs.Pointer;
                        c : C_dblarrs.Pointer;
                        vrblvl : integer32 := 0 ) return integer32;

-- DESCRIPTION :
--   Gateway to the container to manage numerically computed tropisms.

-- ON ENTRY :
--   job   =   1 : initialize the tropisms in standard double precision,
--                 on input in a[0] are the number of tropisms, and
--                 in a[1] the length of the tropisms, in b are then
--                 the winding numbers, as many as the value of a[0],
--                 in c are all the floating-point coefficients,
--                 as many as a[0]*(a[1]+1), for the coordinates of the
--                 tropisms and the errors on the numerical tropisms; 
--             2 : initialize the tropisms in double double precision,
--                 on input in a[0] are the number of tropisms, and
--                 in a[1] the length of the tropisms, in b are then
--                 the winding numbers, as many as the value of a[0],
--                 in c are all the floating-point coefficients,
--                 as many as 2*a[0]*(a[1]+1), for the coordinates of the
--                 tropisms and the errors on the numerical tropisms; 
--             3 : initialize the tropisms in quad double precision,
--                 on input in a[0] are the number of tropisms, and
--                 in a[1] the length of the tropisms, in b are then
--                 the winding numbers, as many as the value of a[0],
--                 in c are all the floating-point coefficients,
--                 as many as 4*a[0](a[1]+1), for the coordinates of the
--                 tropisms and the errors on the numerical tropisms; 
--             4 : store a tropism in standard double precision,
--                 on input in a[0] is the length of the tropism, and
--                 in a[1] the index of the tropism, 
--                 in b is the winding number,
--                 in c are all the floating-point coefficients,
--                 as many as a[0] + 1, with the coordinates of the 
--                 tropisms and the error;
--             5 : store a tropism in double double precision,
--                 on input in a[0] is the length of the tropism, and
--                 in a[1] the index of the tropism, 
--                 in b is the winding number,
--                 in c are all the floating-point coefficients,
--                 as many as 2*a[0] + 2, with the coordinates of the 
--                 tropisms and the error;
--             6 : store a tropism in quad double precision,
--                 on input in a[0] is the length of the tropism, and
--                 in a[1] the index of the tropism, 
--                 in b is the winding number,
--                 in c are all the floating-point coefficients,
--                 as many as 4*a[0] + 4, with the coordinates of the 
--                 tropisms and the error;
--             7 : retrieve all tropisms in standard double precision,
--                 on return in a[0] are the number of tropisms, and
--                 in a[1] the length of each tropism,
--                 in b are the winding numbers for each tropism, and
--                 in c the a[0]*(a[1]+1) coefficients;
--             8 : retrieve all tropisms in double double precision,
--                 on return in a[0] are the number of tropisms, and
--                 in a[1] the length of each tropism,
--                 in b are the winding numbers for each tropism, and
--                 in c the 2*a[0]*(a[1]+1) coefficients;
--             9 : retrieve all tropisms in quad double precision,
--                 on return in a[0] are the number of tropisms, and
--                 in a[1] the length of each tropism,
--                 in b are the winding numbers for each tropism, and
--                 in c the 4*a[0]*(a[1]+1) coefficients;
--            10 : the number of tropisms in standard double precision
--                 is returned in a;
--            11 : the number of tropisms in double double precision
--                 is returned in a;
--            12 : the number of tropisms in quad double precision
--                 is returned in a;
--            13 : retrieve one tropism in standard double precision,
--                 on input in a[0] is the length of the tropism, and
--                 in a[1] is the index of the tropism,
--                 on return in b is the winding number and in c are
--                 the coefficiennts, as many as a[0] + 1, of the
--                 coordinates of the tropism and the error,
--            14 : retrieve one tropism in double double precision,
--                 on input in a[0] is the length of the tropism, and
--                 in a[1] is the index of the tropism,
--                 on return in b is the winding number and in c are
--                 the coefficiennts, as many as 2*a[0] + 2, of the
--                 coordinates of the tropism and the error,
--            15 : retrieve one tropism in quad double precision,
--                 on input in a[0] is the length of the tropism, and
--                 in a[1] is the index of the tropism,
--                 on return in b is the winding number and in c are
--                 the coefficiennts, as many as 4*a[0] + 4, of the
--                 coordinates of the tropism and the error,
--            16 : clear all tropisms in standard double precision;
--            17 : clear all tropisms in double double precision;
--            18 : clear all tropisms in quad double precision;
--            19 : the length of the tropisms in standard double precision
--                 is returned in a;
--            20 : the length of the tropisms in double double precision
--                 is returned in a;
--            21 : the length of the tropisms in quad double precision
--                 is returned in a.

-- ON RETURN :
--   0 if the operation was successful, otherwise something went wrong,
--   or job not in the right range.
