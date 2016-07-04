with text_io;                            use text_io;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Complex_Numbers;
with Standard_Random_Numbers;
with Standard_Complex_Poly_Systems;
with Standard_Complex_Poly_Systems_io;   use Standard_Complex_Poly_Systems_io;
with Standard_Homotopy;
with Standard_Series_Poly_Systems;
with Series_and_Polynomials;
with Series_and_Polynomials_io;          use Series_and_Polynomials_io;
with Series_and_Homotopies;

procedure ts_serhom is

-- DESCRIPTION :
--   Test on the operations in the package series_and_homotopies.

  procedure Standard_Test_Creation ( nq : in integer32 ) is

  -- DESCRIPTION :
  --   Test on converting Standard_Homotopy.Homotopy_System of
  --   range 1..nq to a series system.
  --   The output of put_line(h) is the same as put(s,1),
  --   except for that the last variable in h has become
  --   the first variable in s.

    h : Standard_Complex_Poly_Systems.Poly_Sys(1..nq)
      := Standard_Homotopy.Homotopy_System;
    s : Standard_Series_Poly_Systems.Poly_Sys(1..nq);

  begin
    new_line;
    put_line("The homotopy system :"); put_line(h);
    s := Series_and_Homotopies.Create(h,nq+1,true);
    put_line("The series system :"); put(s,1);
  end Standard_Test_Creation;

  procedure Main is

  -- DESCRIPTION :
  --   Test on converting a homotopy system to a series system.

    target,start : Standard_Complex_Poly_Systems.Link_to_Poly_Sys;
    k : constant natural32 := 2;
    gamma : constant Standard_Complex_Numbers.Complex_Number
          := Standard_Random_Numbers.Random1;

  begin
    new_line;
    put_line("Testing the creation of a homotopies as a series system ...");
    new_line;
    put_line("Reading the target system ..."); get(target);
    new_line;
    put_line("Reading the start system ..."); get(start);
    Standard_Homotopy.Create(target.all,start.all,k,gamma);
    Standard_Test_Creation(target'last);
  end Main;

begin
  Main;
end ts_serhom;
