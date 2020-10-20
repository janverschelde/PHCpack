with Test_mtNewton_Convolutions;

procedure ts_mtnewton is

-- DESCRIPTION :
--   Calls the test procedures on Newton's method for power series
--   with the reverse mode of algorithmic differentation
--   and linearization to solve the matrix series equations,
--   in double, double double, and quad double arithmetic,
--   with multitasking for shared memory parallel computers.

begin
  Test_mtNewton_Convolutions.Main;
end ts_mtnewton;
