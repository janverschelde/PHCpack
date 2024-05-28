with text_io;                            use text_io;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;
with Standard_Complex_Numbers;
with Standard_Random_Numbers;
with Standard_Natural_Vectors;
with Standard_Floating_Vectors;
with Standard_Floating_Vectors_io;       use Standard_Floating_Vectors_io;
with Symbol_Table;
with Standard_Complex_Polynomials;       use Standard_Complex_Polynomials;
with Standard_Complex_Poly_Systems;      use Standard_Complex_Poly_Systems;
with Standard_Complex_Poly_Systems_io;   use Standard_Complex_Poly_Systems_io;
with Double_Taylor_Developments;         use Double_Taylor_Developments;

procedure ts_fliphom is

-- DESCRIPTION :
--   Prototypes a polyhedral homotopy on extending an initial form system
--   with one monomial, as geometrically in a bistellar flip.

  function Random_Polynomial
             ( cff : Standard_Floating_Vectors.Vector ) return Poly is

  -- DESCRIPTION :
  --   Returns the polynomial supported on (0,0), (1,0), (0,1), (1,1),
  --   with random complex coefficients, using the coefficients cff
  --   of the Taylor developments of the coefficient of (1,1).

    res : Poly;
    mon : Term;

  begin
    mon.dg := new Standard_Natural_Vectors.Vector(1..3);
    mon.dg(1) := 0;
    mon.dg(2) := 0;
    mon.dg(3) := 0;
    mon.cf := Standard_Random_Numbers.Random1;
    res := Create(mon);
    mon.dg(1) := 0;
    mon.dg(2) := 1;
    mon.dg(3) := 0;
    mon.cf := Standard_Random_Numbers.Random1;
    Add(res, mon);
    mon.dg(1) := 0;
    mon.dg(2) := 0;
    mon.dg(3) := 1;
    mon.cf := Standard_Random_Numbers.Random1;
    Add(res, mon);
    for k in cff'range loop
      mon.dg(1) := natural32(k);
      mon.dg(2) := 1;
      mon.dg(3) := 1;
      mon.cf := Standard_Complex_Numbers.Create(cff(k));
      Add(res, mon);
    end loop;
    Clear(mon);
    return res;
  end Random_Polynomial;

  function Random_Square_System 
             ( cff : Standard_Floating_Vectors.Vector ) return Poly_Sys is

  -- DESCRIPTION :
  --   Returns the system supported on (0,0), (1,0), (0,1), (1,1),
  --   with random complex coefficients, using the coefficients of
  --   the Taylor developments for the coefficient of (1,1).

    res : Poly_Sys(1..2);

  begin
    res(1) := Random_Polynomial(cff);
    res(2) := Random_Polynomial(cff);
    return res;
  end Random_Square_System;

  procedure Double_Test
              ( deg : in integer32; alpha, point : in double_float ) is

  -- DESCRIPTION :
  --   Generates a random square system and builds a homotopy.

    cff : constant Standard_Floating_Vectors.Vector(0..deg)
        := Double_Taylor_Coefficients(deg,alpha,point);
    sys : constant Poly_Sys := Random_Square_System(cff);

  begin
    Symbol_Table.Init(3);
    Symbol_Table.Add_String("t");
    Symbol_Table.Add_String("x");
    Symbol_Table.Add_String("y");
    put_line("The coefficients of the Taylor series :");
    put_line(cff);
    put_line("The homotopy :");
    put_line(sys);
  end Double_Test;

  procedure Main is

  -- DESCRIPTION :
  --   Prompts for the main parameters: truncation degree
  --   and exponent of the power of the continuation parameter.

    deg : integer32 := 0;
    alpha,point : double_float := 0.0;

  begin
    new_line;
    put("Give the truncation degree : "); get(deg);
    put("Give the positive real power : "); get(alpha);
    put("Give the positive real point : "); get(point);
    new_line;
    put("-> the truncation degree : "); put(deg,1); new_line;
    put("-> power of the monomial :"); put(alpha); new_line;
    put("-> point of development  :"); put(point); new_line;
    Double_Test(deg,alpha,point);
  end Main;

begin
  Main;
end ts_fliphom;
