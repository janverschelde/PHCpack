with text_io;                           use text_io;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with Number_of_Cores;

procedure ts_corecount is

-- DESCRIPTION :
--   Writes the number of cores available on the system.

  procedure Main is

    cnt : constant integer32 := Number_of_Cores;

  begin
    put("Number of cores : "); put(cnt,1); new_line;
  end Main;

begin
  Main;
end ts_corecount;
