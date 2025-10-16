with Ada.Text_IO;
with Test_Weighted_Assignment;

procedure ts_hunwas is

begin
  Ada.Text_IO.put_line("Testing the Hungarian algorithm ...");
  Test_Weighted_Assignment.main;
end ts_hunwas;
