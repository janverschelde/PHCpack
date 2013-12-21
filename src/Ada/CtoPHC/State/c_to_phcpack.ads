with Standard_Integer_Numbers;          use Standard_Integer_Numbers;

function C_to_PHCpack ( job : integer32 ) return integer32;

-- DESCRIPTION :
--   Interactive gateway to the core operations in PHCpack.
--   Calling this routine with n = 0 displays the menu.

-- ON ENTRY :
--   job      =  0 : displays the menu of available operations;
--            =  1 : read target polynomial system;
--            =  2 : write target polynomial system;
--            =  3 : read start polynomial system;
--            =  4 : write start polynomial system;
--            =  5 : write start solutions;
--            =  6 : solve by homotopy continuation;
--            =  7 : write the target solutions;
--            =  8 : clear the data in PHCpack_Operations;
--            =  9 : define the output file;
--            = 11 : read double double target polynomial system;
--            = 12 : write double double target polynomial system;
--            = 13 : read double double start polynomial system;
--            = 14 : write double double start polynomial system;
--            = 15 : write double double start solutions;
--            = 16 : solve by double double homotopy continuation;
--            = 17 : write the double double target solutions;
--            = 18 : clear the double double data;
--            = 21 : read quad double target polynomial system;
--            = 22 : write quad double target polynomial system;
--            = 23 : read quad double start polynomial system;
--            = 24 : write quad double start polynomial system;
--            = 25 : write quad double start solutions;
--            = 26 : solve by quad double homotopy continuation;
--            = 27 : write the quad double target solutions;
--            = 28 : clear the quad double data.

-- ON RETURN :
--   0 if the operation was successful, otherwise something went wrong,
--   e.g.: job not in the range 0..9, 11..18, or 21..28.
