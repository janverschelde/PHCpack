with Ada.Text_IO;
with Ada.Command_Line;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Multitasked_Geometric_Products;

procedure ts_mtgeoprod is

-- DESCRIPTION :
--   Runs tests on the multitasked geometric inner products,
--   expecting up to two optional arguments at the command line:
--   (1) the number of threads,
--   (2) the precision: 2 for double double, 4 for quad double.

  argcnt : constant integer32 := integer32(Ada.Command_Line.Argument_Count);

  procedure Double_Run is

  -- DESCRIPTION :
  --   Runs a test in double precision.

  begin
    if argcnt = 0 then
      Multitasked_Geometric_Products.Double_Sequential_Run;
    else
      declare
        arg1 : constant String := Ada.Command_Line.Argument(1);
        nprc : constant integer32 := integer32'value(arg1);
      begin
        Ada.Text_IO.Put_Line("Running with" & nprc'Image & " threads ...");
        Multitasked_Geometric_Products.Double_Parallel_Run(nprc);
      end;
    end if;
  end Double_Run;

  procedure Double_Double_Run is

  -- DESCRIPTION :
  --   Runs a test in double double precision.
  
  begin
    if argcnt = 0 then
      Multitasked_Geometric_Products.Double_Double_Sequential_Run;
    else
      declare
        arg1 : constant String := Ada.Command_Line.Argument(1);
        nprc : constant integer32 := integer32'value(arg1);
      begin
        Ada.Text_IO.Put_Line("Running with" & nprc'Image & " threads ...");
        Multitasked_Geometric_Products.Double_Double_Parallel_Run(nprc);
      end;
    end if;
  end Double_Double_Run;

  procedure Quad_Double_Run is

  -- DESCRIPTION :
  --   Runs a test in quad double precision.
  
  begin
    if argcnt = 0 then
      Multitasked_Geometric_Products.Quad_Double_Sequential_Run;
    else
      declare
        arg1 : constant String := Ada.Command_Line.Argument(1);
        nprc : constant integer32 := integer32'value(arg1);
      begin
        Ada.Text_IO.Put_Line("Running with" & nprc'Image & " threads ...");
        Multitasked_Geometric_Products.Quad_Double_Parallel_Run(nprc);
      end;
    end if;
  end Quad_Double_Run;

begin
  if argcnt <= 1 then
    Ada.Text_IO.Put_Line("Running in double precision ...");
    Double_Run;
  else
    declare
      arg2 : constant String := Ada.Command_Line.Argument(2);
      prec : constant integer32 := integer32'value(arg2);
    begin
      if prec = 2 then
        Ada.Text_IO.Put_Line("Running in double double precision ...");
        Double_Double_Run;
      elsif prec = 4 then
        Ada.Text_IO.Put_Line("Running in quad double precision ...");
        Quad_Double_Run;
      end if;
    end;
  end if;
end ts_mtgeoprod;
