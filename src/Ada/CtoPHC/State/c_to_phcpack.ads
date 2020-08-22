with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;

function C_to_PHCpack
           ( job : integer32;
             number_of_tasks : natural32;
             vrblvl : integer32 := 0 ) return integer32;

-- DESCRIPTION :
--   Interactive gateway to the basic path tracking operations in PHCpack.
--   Calling this routine with n = 0 displays the menu.

-- ON ENTRY :
--   job      =  0 : displays the menu of available operations;
--            =  1 : read target polynomial system;
--            =  2 : write target polynomial system;
--            =  3 : read start polynomial system;
--            =  4 : write start polynomial system;
--            =  5 : write start solutions;
--            =  6 : solve by standard homotopy continuation,
--                   using as many tasks as the value of number_of_tasks,
--                   if that value is zero, then no multitasking is used;
--            =  7 : write the target solutions;
--            =  8 : clear the data in PHCpack_Operations;
--            =  9 : define the output file;
--            = 11 : read double double target polynomial system;
--            = 12 : write double double target polynomial system;
--            = 13 : read double double start polynomial system;
--            = 14 : write double double start polynomial system;
--            = 15 : write double double start solutions;
--            = 16 : solve by double double homotopy continuation,
--                   using as many tasks as the value of number_of_tasks,
--                   if that value is zero, then no multitasking is used;
--            = 17 : write the double double target solutions;
--            = 18 : clear the double double data;
--            = 21 : read quad double target polynomial system;
--            = 22 : write quad double target polynomial system;
--            = 23 : read quad double start polynomial system;
--            = 24 : write quad double start polynomial system;
--            = 25 : write quad double start solutions;
--            = 26 : solve by quad double homotopy continuation,
--                   using as many tasks as the value of number_of_tasks,
--                   if that value is zero, then no multitasking is used;
--            = 27 : write the quad double target solutions;
--            = 28 : clear the quad double data;
--            = 29 : read standard double Laurent start system;
--            = 30 : write standard double Laurent start system;
--            = 31 : read standard double Laurent target system;
--            = 32 : write standard double Laurent target system;
--            = 33 : read double double Laurent start system;
--            = 34 : write double double Laurent start system;
--            = 35 : read double double Laurent target system;
--            = 36 : write double double Laurent target system;
--            = 37 : read quad double Laurent start system;
--            = 38 : write quad double Laurent start system;
--            = 39 : read quad double Laurent target system;
--            = 40 : write quad double Laurent target system;
--            = 41 : clear data for Laurent homotopy with standard doubles;
--            = 42 : clear data for Laurent homotopy with double doubles;
--            = 43 : clear data for Laurent homotopy with quad doubles;
--            = 44 : solve stored Laurent target system using the stored
--                   Laurent system in standard double precision,
--                   using as many tasks as the value of number_of_tasks,
--                   if that value is zero, then no multitasking is used;
--            = 45 : solve stored Laurent target system using the stored
--                   Laurent system in double double precision,
--                   using as many tasks as the value of number_of_tasks,
--                   if that value is zero, then no multitasking is used;
--            = 46 : solve stored Laurent target system using the stored
--                   Laurent system in quad double precision,
--                   using as many tasks as the value of number_of_tasks,
--                   if that value is zero, then no multitasking is used.

-- ON RETURN :
--   0 if the operation was successful, otherwise something went wrong,
--   e.g.: job not in the range 0..9, 11..18, or 21..43.
