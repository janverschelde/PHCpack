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
with Standard_Integer_Vectors;
with Standard_Complex_Vectors;
with Standard_Complex_Vectors_io;        use Standard_Complex_Vectors_io;
with Standard_Complex_Matrices;
with Standard_Complex_Matrices_io;       use Standard_Complex_Matrices_io;
with Standard_Complex_Linear_Solvers;    use Standard_Complex_Linear_Solvers;
with Symbol_Table;
with Standard_Complex_Polynomials;       use Standard_Complex_Polynomials;
with Standard_Complex_Poly_Systems;      use Standard_Complex_Poly_Systems;
with Standard_Complex_Poly_Systems_io;   use Standard_Complex_Poly_Systems_io;
with Standard_Complex_Poly_SysFun;
with Standard_Complex_Solutions;
with Standard_Homotopy;
with Standard_CSeries_Poly_Systems;
with Standard_CSeries_Poly_Systems_io;   use Standard_CSeries_Poly_Systems_io;
-- with Standard_Fabry_on_Homotopy;
with Double_Taylor_Developments;         use Double_Taylor_Developments;
with Double_Taylor_Homotopies;           use Double_Taylor_Homotopies;
with Double_Taylor_Homotopies_io;        use Double_Taylor_Homotopies_io;
with Taylor_Homotopy_Series;

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

  procedure Solve_Linear_System
              ( p : in Poly_Sys; s : out Standard_Complex_Vectors.Vector ) is 

  -- DESCRIPTION :
  --   Extracts the linear part of the system in p
  --   and return the solution of the linear system in s.

    mat : Standard_Complex_Matrices.Matrix(1..2,1..2);
    rhs : Standard_Complex_Vectors.Vector(1..2);
    deg : Standard_Complex_Polynomials.Degrees
        := new Standard_Natural_Vectors.Vector(1..3);
    piv : Standard_Integer_Vectors.Vector(1..2);
    info : integer32;
    wrk : Standard_Complex_Vectors.Vector(1..3);
    eva : Standard_Complex_Vectors.Vector(1..2);

    use Standard_Complex_Numbers;

  begin
    deg(1) := 0;
    deg(2) := 0;
    deg(3) := 0;
    rhs(1) := -Coeff(p(1),deg);
    rhs(2) := -Coeff(p(2),deg);
    put_line("The right hand side vector :");
    put_line(rhs);
    deg(2) := 1;
    mat(1,1) := Coeff(p(1),deg);
    mat(2,1) := Coeff(p(2),deg);
    deg(2) := 0;
    deg(3) := 1;
    mat(1,2) := Coeff(p(1),deg);
    mat(2,2) := Coeff(p(2),deg);
    put_line("The coefficient matrix :"); put(mat);
    lufac(mat,2,piv,info);
    lusolve(mat,2,piv,rhs);
    Clear(deg);
    wrk(1) := Standard_Complex_Numbers.Create(0.0);
    wrk(2) := rhs(1);
    wrk(3) := rhs(2);
    put_line("The solution :");
    put_line(rhs);
    eva := Standard_Complex_Poly_SysFun.Eval(p,wrk);
    put_line("The value at the system :");
    put_line(eva);
    s := rhs;
  end Solve_Linear_System;

  function Make_Taylor_Monomial_Vector
             ( p : Poly; deg : integer32; alpha,point : double_float )
             return Taylor_Monomial_Vector is

   -- len : constant integer32 := integer32(Number_of_Terms(p));
    len : constant integer32 := 4;
    res : Taylor_Monomial_Vector(1..len);
    idx : integer32 := 0;

    procedure Visit_Term ( t : in Term; continue : out boolean ) is

      tm : Taylor_Monomial(2,deg);
      dg : Standard_Integer_Vectors.Vector(1..2);

    begin
      if t.dg(1) = 0 then
        idx := idx + 1; 
        dg(1) := integer32(t.dg(2));
        dg(2) := integer32(t.dg(3));
        if dg(1) = 1 and dg(2) = 1
         then tm := Make(deg,alpha,point,t.cf,dg);
         else tm := Make(deg,0.0,point,t.cf,dg);
        end if;
        res(idx) := new Taylor_Monomial'(tm);
      end if;
      continue := true;
    end Visit_Term;
    procedure Visit_Terms is new Visiting_Iterator(Visit_Term);

  begin
    Visit_Terms(p);
    return res;
  end Make_Taylor_Monomial_Vector;

  function Make_Taylor_Homotopy
             ( p : Poly_Sys; deg : integer32; alpha,point : double_float )
             return Taylor_Homotopy is

  -- DESCRIPTION :
  --   Makes a Taylor homotopy, using t^alpha developed at the point,
  --   as power series truncated at degree deg,
  --   for the monomial with support (1, 1).

    res : Taylor_Homotopy(p'range);

  begin
    for i in p'range loop
      declare
        tmv : constant Taylor_Monomial_Vector
            := Make_Taylor_Monomial_Vector(p(i),deg,alpha,point);
      begin
        res(i) := new Taylor_Monomial_Vector'(tmv);
      end;
    end loop;
    return res;
  end Make_Taylor_Homotopy;

  procedure Double_Test
              ( deg : in integer32; alpha, point : in double_float ) is

  -- DESCRIPTION :
  --   Generates a random square system and builds a homotopy.

    cff : constant Standard_Floating_Vectors.Vector(0..deg)
        := Double_Taylor_Coefficients(deg,alpha,point);
    sys : constant Poly_Sys := Random_Square_System(cff);
    sol : Standard_Complex_Solutions.Solution(2);
    sols : Standard_Complex_Solutions.Solution_List;
    thm : Taylor_Homotopy(sys'range)
        := Make_Taylor_Homotopy(sys,deg,alpha,point);
    hom : Standard_CSeries_Poly_Systems.Poly_Sys(sys'range)
        := Taylor_Homotopy_Series.Make(thm);

  begin
    Symbol_Table.Init(3);
    Symbol_Table.Add_String("t");
    Symbol_Table.Add_String("x");
    Symbol_Table.Add_String("y");
    put_line("The coefficients of the Taylor series :");
    put_line(cff);
    put_line("The homotopy :");
    put_line(sys);
    Standard_Homotopy.Create(sys,1);
    Solve_Linear_System(sys,sol.v);
    sol.m := 1;
    sol.t := Standard_Complex_Numbers.Create(point);
    Standard_Complex_Solutions.Add(sols,sol);
   -- Standard_Fabry_on_Homotopy.Run(0,2,1,deg,sols);
    put_line("The Taylor homotopy :"); put(thm);
    put_line("The Taylor homotopy as series system :"); put(hom);
    Clear(thm);
    Standard_CSeries_Poly_Systems.Clear(hom);
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
