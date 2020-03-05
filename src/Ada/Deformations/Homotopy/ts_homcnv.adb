with text_io;                           use text_io;
with Communications_with_User;          use Communications_with_User;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Complex_Numbers;
with DoblDobl_Complex_Numbers;
with QuadDobl_Complex_Numbers;
with Standard_Complex_Vectors;
with Standard_Complex_Vectors_io;       use Standard_Complex_Vectors_io;
with DoblDobl_Complex_Vectors;
with DoblDobl_Complex_Vectors_io;       use DoblDobl_Complex_Vectors_io;
with QuadDobl_Complex_Vectors;
with QuadDobl_Complex_Vectors_io;       use QuadDobl_Complex_Vectors_io;
with Standard_Integer_VecVecs_io;
with Standard_Complex_Poly_Systems;
with DoblDobl_Complex_Poly_Systems;
with QuadDobl_Complex_Poly_Systems;
with Standard_Complex_Solutions;
with DoblDobl_Complex_Solutions;
with QuadDobl_Complex_Solutions;
with Standard_Homotopy;
with DoblDobl_Homotopy;
with QuadDobl_Homotopy;
with Artificial_Parameter_Homotopy_io;
with Standard_CSeries_Poly_Systems;
with Standard_CSeries_Poly_Systems_io;  use Standard_CSeries_Poly_Systems_io;
with DoblDobl_CSeries_Poly_Systems;
with DoblDobl_CSeries_Poly_Systems_io;  use DoblDobl_CSeries_Poly_Systems_io;
with QuadDobl_CSeries_Poly_Systems;
with QuadDobl_CSeries_Poly_Systems_io;  use QuadDobl_CSeries_Poly_Systems_io;
with Series_and_Homotopies;
with Complex_Series_and_Polynomials;
with Standard_Speelpenning_Convolutions;
with DoblDobl_Speelpenning_Convolutions;
with QuadDobl_Speelpenning_Convolutions;
with System_Convolution_Circuits;

procedure ts_homcnv is

-- DESCRIPTION :
--   Test on making convolution circuits for a homotopy,
--   in double, double double, or quad double precision.

  procedure Standard_Test
               ( nq : in integer32;
                 sols : in Standard_Complex_Solutions.Solution_List ) is

  -- DESCRIPTION :
  --   Converts the stored homotopy into convolution circuits
  --   and then verifies that the solutions are proper start solutions.
  --   The number of equations is given in nq.

    ls : constant Standard_Complex_Solutions.Link_to_Solution
       := Standard_Complex_Solutions.Head_Of(sols);
    nv : constant integer32 := ls.n;
    h : constant Standard_Complex_Poly_Systems.Poly_Sys(1..nq)
      := Standard_Homotopy.Homotopy_System;
    s : Standard_CSeries_Poly_Systems.Poly_Sys(1..nq)
      := Series_and_Homotopies.Create(h,nv+1);
    c : Standard_Speelpenning_Convolutions.Circuits(1..nq);
    deg : integer32 := 0;
    zero : constant Standard_Complex_Numbers.Complex_Number
         := Standard_Complex_Numbers.Create(0.0);
    eva : Standard_Complex_Vectors.Vector(1..nq);

  begin
    put_line("The homotopy as a series system : "); put(s);
    put("Give the degree of the series : "); get(deg);
    Complex_Series_and_Polynomials.Set_Degree(s,deg);
    c := System_Convolution_Circuits.Make_Convolution_Circuits(s);
    put_line("The exponents in the circuits :");
    for k in c'range loop
      Standard_Integer_VecVecs_io.put(c(k).xps);
    end loop;
    eva := Standard_Speelpenning_Convolutions.Eval(c,ls.v,zero);
    put_line("Value at the circuit at the first start solution :");
    put_line(eva);
  end Standard_Test;

  procedure DoblDobl_Test
               ( nq : in integer32;
                 sols : in DoblDobl_Complex_Solutions.Solution_List ) is

  -- DESCRIPTION :
  --   Converts the stored homotopy into convolution circuits
  --   and then verifies that a solution is a proper start solution.
  --   The number of equations is given in nq.

    ls : constant DoblDobl_Complex_Solutions.Link_to_Solution
       := DoblDobl_Complex_Solutions.Head_Of(sols);
    nv : constant integer32 := ls.n;
    h : constant DoblDobl_Complex_Poly_Systems.Poly_Sys(1..nq)
      := DoblDobl_Homotopy.Homotopy_System;
    s : DoblDobl_CSeries_Poly_Systems.Poly_Sys(1..nq)
      := Series_and_Homotopies.Create(h,nv+1);
    c : DoblDobl_Speelpenning_Convolutions.Circuits(1..nq);
    deg : integer32 := 0;
    zero : constant DoblDobl_Complex_Numbers.Complex_Number
         := DoblDobl_Complex_Numbers.Create(integer(0));
    eva : DoblDobl_Complex_Vectors.Vector(1..nq);

  begin
    put_line("The homotopy as a series system : "); put(s);
    put("Give the degree of the series : "); get(deg);
    Complex_Series_and_Polynomials.Set_Degree(s,deg);
    c := System_Convolution_Circuits.Make_Convolution_Circuits(s);
    put_line("The exponents in the circuits :");
    for k in c'range loop
      Standard_Integer_VecVecs_io.put(c(k).xps);
    end loop;
    eva := DoblDobl_Speelpenning_Convolutions.Eval(c,ls.v,zero);
    put_line("Value at the circuit at the first start solution :");
    put_line(eva);
  end DoblDobl_Test;

  procedure QuadDobl_Test
               ( nq : in integer32;
                 sols : in QuadDobl_Complex_Solutions.Solution_List ) is

  -- DESCRIPTION :
  --   Converts the stored homotopy into convolution circuits
  --   and then verifies that a solution is a proper start solution.
  --   The number of equations is given in nq.

    ls : constant QuadDobl_Complex_Solutions.Link_to_Solution
       := QuadDobl_Complex_Solutions.Head_Of(sols);
    nv : constant integer32 := ls.n;
    h : constant QuadDobl_Complex_Poly_Systems.Poly_Sys(1..nq)
      := QuadDobl_Homotopy.Homotopy_System;
    s : QuadDobl_CSeries_Poly_Systems.Poly_Sys(1..nq)
      := Series_and_Homotopies.Create(h,nv+1);
    c : QuadDobl_Speelpenning_Convolutions.Circuits(1..nq);
    deg : integer32 := 0;
    zero : constant QuadDobl_Complex_Numbers.Complex_Number
         := QuadDobl_Complex_Numbers.Create(integer(0));
    eva : QuadDobl_Complex_Vectors.Vector(1..nq);

  begin
    put_line("The homotopy as a series system : "); put(s);
    put("Give the degree of the series : "); get(deg);
    Complex_Series_and_Polynomials.Set_Degree(s,deg);
    c := System_Convolution_Circuits.Make_Convolution_Circuits(s);
    put_line("The exponents in the circuits :");
    for k in c'range loop
      Standard_Integer_VecVecs_io.put(c(k).xps);
    end loop;
    eva := QuadDobl_Speelpenning_Convolutions.Eval(c,ls.v,zero);
    put_line("Value at the circuit at the first start solution :");
    put_line(eva);
  end QuadDobl_Test;

  procedure Standard_Main is

  -- DESCRIPTION :
  --   Prompts the user for an artificial-parameter homotopy for testing
  --   the making of convolution circuits, in double precision.

    target,start : Standard_Complex_Poly_Systems.Link_to_Poly_Sys;
    sols : Standard_Complex_Solutions.Solution_List;
    tpow : constant natural32 := 2;
    gamma : constant Standard_Complex_Numbers.Complex_Number
          := Standard_Complex_Numbers.Create(1.0);

  begin
    Artificial_Parameter_Homotopy_io.get(target,start,sols);
    Standard_Homotopy.Create(target.all,start.all,tpow,gamma);
    Standard_Test(target'last,sols);
  end Standard_Main;

  procedure DoblDobl_Main is

  -- DESCRIPTION :
  --   Prompts the user for an artificial-parameter homotopy for testing
  --   the making of convolution circuits, in double double precision.

    target,start : DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
    sols : DoblDobl_Complex_Solutions.Solution_List;
    tpow : constant natural32 := 2;
    gamma : constant DoblDobl_Complex_Numbers.Complex_Number
          := DoblDobl_Complex_Numbers.Create(integer32(1));

  begin
    Artificial_Parameter_Homotopy_io.get(target,start,sols);
    DoblDobl_Homotopy.Create(target.all,start.all,tpow,gamma);
    DoblDobl_Test(target'last,sols);
  end DoblDobl_Main;

  procedure QuadDobl_Main is

  -- DESCRIPTION :
  --   Prompts the user for an artificial-parameter homotopy for testing
  --   the making of convolution circuits, in quad double precision.

    target,start : QuadDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
    sols : QuadDobl_Complex_Solutions.Solution_List;
    tpow : constant natural32 := 2;
    gamma : constant QuadDobl_Complex_Numbers.Complex_Number
          := QuadDobl_Complex_Numbers.Create(integer32(1));

  begin
    Artificial_Parameter_Homotopy_io.get(target,start,sols);
    QuadDobl_Homotopy.Create(target.all,start.all,tpow,gamma);
    QuadDobl_Test(target'last,sols);
  end QuadDobl_Main;

  procedure Main is

  -- DESCRIPTION :
  --   Prompts the user for the working precision and then launches
  --   the test corresponding to that precision.

    precision : character;

  begin
    new_line;
    put_line("Testing convolution circuits for homotopies ...");
    new_line;
    put_line("MENU for the working precision :");
    put_line("  0. double precision");
    put_line("  1. double double precision");
    put_line("  2. quad double precision");
    put("Type 0, 1, or 2 to select the working precision : ");
    Ask_Alternative(precision,"012");
    case precision is
      when '0' => Standard_Main;
      when '1' => DoblDobl_Main;
      when '2' => QuadDobl_Main;
      when others => null;
    end case;
  end Main;

begin
  Main;
end ts_homcnv;
