with text_io;                             use text_io;
with Standard_Floating_Numbers;           use Standard_Floating_Numbers;
with Standard_Complex_Solutions;          use Standard_Complex_Solutions;

procedure Filter_and_Split_Solutions
              ( file : in file_type; sols : in Solution_List;
                n,k : in natural; tol : in double_float;
                sols0,sols1 : out Solution_List );

-- DESCRIPTION :
--   Removes the spurious solutions from the solution list sols and
--   breaks it up into two lists, with (n+k)-th variable as zero component
--   and nonzero component.

-- ON ENTRY :
--   file       to write the results on;
--   sols       computed solutions at the end of the paths;
--   n          original dimension;
--   k          number of random slices and variables added;
--   tol        everything that has absolute value less than tol is zero.

-- ON RETURN :
--   sols0      solutions with n-th component less than tol;
--   sols1      solutions with n-th component larger than tol.
