with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with C_Integer_Arrays;                  use C_Integer_Arrays;
with C_Double_Arrays;                   use C_Double_Arrays;

function use_witsols ( job : integer32;
                       a : C_intarrs.Pointer;
                       b : C_intarrs.Pointer;
                       c : C_dblarrs.Pointer;
                       vrblvl : integer32 := 0 ) return integer32;

-- DESCRIPTION :
--   Provides a gateway from C to the cascade homotopies, the filter,
--   and factor operations for a numerical irreducible decomposition.
--   The verbose level is given by the parameter vrblvl.

-- ON ENTRY :
--   job =  0 : numerical irreducible decomposition with standard doubles,
--              on the system in the standard polynomial systems container,
--              the input parameters are
--              a[0] : the top dimension of the solution set,
--              a[1] : the number of tasks, 0 for no multitasking,
--              a[2] : 0 or 1 if the witness supersets need to be filtered,
--              a[3] : 0 or 1 if the witness sets need to be factored,
--              b[0] : the verbose flag, 0 for false, 1 for true.
--              If verbose, then intermediate results are written to screen.
--   job =  1 : numerical irreducible decomposition with standard doubles,
--              on the system in the standard Laurent systems container,
--              the input parameters are
--              a[0] : the top dimension of the solution set,
--              a[1] : the number of tasks, 0 for no multitasking,
--              a[2] : 0 or 1 if the witness supersets need to be filtered,
--              a[3] : 0 or 1 if the witness sets need to be factored,
--              b[0] : the verbose flag, 0 for false, 1 for true.
--              If verbose, then intermediate results are written to screen.
--   job =  2 : numerical irreducible decomposition with double doubles,
--              on the system in the dobldobl polynomial systems container,
--              the input parameters are
--              a[0] : the top dimension of the solution set,
--              a[1] : the number of tasks, 0 for no multitasking,
--              a[2] : 0 or 1 if the witness supersets need to be filtered,
--              a[3] : 0 or 1 if the witness sets need to be factored,
--              b[0] : the verbose flag, 0 for false, 1 for true.
--              If verbose, then intermediate results are written to screen.
--   job =  3 : numerical irreducible decomposition with double doubles,
--              on the system in the dobldobl Laurent systems container,
--              the input parameters are
--              a[0] : the top dimension of the solution set,
--              a[1] : the number of tasks, 0 for no multitasking,
--              a[2] : 0 or 1 if the witness supersets need to be filtered,
--              a[3] : 0 or 1 if the witness sets need to be factored,
--              b[0] : the verbose flag, 0 for false, 1 for true.
--              If verbose, then intermediate results are written to screen.
--   job =  4 : numerical irreducible decomposition with quad doubles,
--              on the system in the quaddobl polynomial systems container,
--              the input parameters are
--              a[0] : the top dimension of the solution set,
--              a[1] : the number of tasks, 0 for no multitasking,
--              a[2] : 0 or 1 if the witness supersets need to be filtered,
--              a[3] : 0 or 1 if the witness sets need to be factored,
--              b[0] : the verbose flag, 0 for false, 1 for true.
--              If verbose, then intermediate results are written to screen.
--   job =  5 : numerical irreducible decomposition with quad doubles,
--              on the system in the quaddobl Laurent systems container,
--              the input parameters are
--              a[0] : the top dimension of the solution set,
--              a[1] : the number of tasks, 0 for no multitasking,
--              a[2] : 0 or 1 if the witness supersets need to be filtered,
--              a[3] : 0 or 1 if the witness sets need to be factored,
--              b[0] : the verbose flag, 0 for false, 1 for true.
--              If verbose, then intermediate results are written to screen.
--   job =  6 : copies the witness set from the stored witness solutions
--              into the standard polynomial system and solutions containers,
--              on input in a[0] is the dimension of the witness set.
--   job =  7 : copies the witness set from the stored witness solutions
--              into the standard Laurent system and solutions containers,
--              on input in a[0] is the dimension of the witness set.
--   job =  8 : copies the witness set from the stored witness solutions
--              into the dobldobl polynomial system and solutions containers,
--              on input in a[0] is the dimension of the witness set.
--   job =  9 : copies the witness set from the stored witness solutions
--              into the dobldobl Laurent system and solutions containers,
--              on input in a[0] is the dimension of the witness set.
--   job = 10 : copies the witness set from the stored witness solutions
--              into the quaddobl polynomial system and solutions containers,
--              on input in a[0] is the dimension of the witness set.
--   job = 11 : copies the witness set from the stored witness solutions
--              into the quaddobl Laurent system and solutions containers,
--              on input in a[0] is the dimension of the witness set.
--   job = 12 : clears the witness solutions for standard double precision.
--   job = 13 : clears the witness solutions for double double precision.
--   job = 14 : clears the witness solutions for quad double precision.
--   job = 15 : given in b[0], the number of characters in the string,
--              return the string representation of the factors in b,
--              returns in a[0] the number of characters in b.

-- ON RETURN :
--   0 if the operation was successful, otherwise something went wrong,
--   e.g.: job not in the right range.
