with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Number_of_Cores;
with Test_Vectored_Hexa_Doubles;

procedure ts_mtvhda is

-- DESCRIPTION :
--   Code to run a noninteractive wall clock time test
--   on vectorized hexa double arithmetic.

  nt : constant integer32 := Number_of_Cores;

begin
  Test_Vectored_Hexa_Doubles.Wall_Time_Parallel_Test(nt);
end ts_mtvhda;
