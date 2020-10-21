with Test_mtHessian_Circuits;

procedure ts_mthesscrc is

-- DESCRIPTION :
--   Calls the main test on the Hessian criterion,
--   for systems of complex circuits,
--   in double, double double, and quad double arithmetic,
--   with multitasking for shared memory parallel computers.

begin
  Test_mtHessian_Circuits.Main;
end ts_mthesscrc;
