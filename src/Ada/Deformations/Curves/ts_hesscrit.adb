with text_io;                           use text_io;
with Communications_with_User;          use Communications_with_User;
with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;       use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;         use Standard_Floating_Numbers;
with Standard_Floating_Numbers_io;      use Standard_Floating_Numbers_io;
with Double_Double_Numbers;             use Double_Double_Numbers;
with Double_Double_Numbers_io;          use Double_Double_Numbers_io;
with Quad_Double_Numbers;               use Quad_Double_Numbers;
with Quad_Double_Numbers_io;            use Quad_Double_Numbers_io;
with Standard_Complex_Vectors;
with DoblDobl_Complex_Vectors;
with QuadDobl_Complex_Vectors;
with Standard_Complex_Polynomials;
with Standard_Complex_Poly_Systems;
with Standard_Complex_Poly_Systems_io;  use Standard_Complex_Poly_Systems_io;
with DoblDobl_Complex_Polynomials;
with DoblDobl_Complex_Poly_Systems;
with DoblDobl_Complex_Poly_Systems_io;  use DoblDobl_Complex_Poly_Systems_io;
with QuadDobl_Complex_Polynomials;
with QuadDobl_Complex_Poly_Systems;
with QuadDobl_Complex_Poly_Systems_io;  use QuadDobl_Complex_Poly_Systems_io;
with Standard_Complex_Jaco_Matrices;
with Standard_Complex_Hessians;
with DoblDobl_Complex_Jaco_Matrices;
with DoblDobl_Complex_Hessians;
with QuadDobl_Complex_Jaco_Matrices;
with QuadDobl_Complex_Hessians;
with Standard_Complex_Solutions;
with DoblDobl_Complex_Solutions;
with QuadDobl_Complex_Solutions;
with Symbol_Table;
with Standard_Homotopy;
with DoblDobl_Homotopy;
with QuadDobl_Homotopy;
with Homotopy_Series_Readers;
with Singular_Values_of_Hessians;

procedure ts_hesscrit is

-- DESCRIPTION :
--   Interactive development testing of the Hessian criterion
--   for solution curves defined by polynomial homotopies.

  procedure Standard_Evaluate
              ( jm : in Standard_Complex_Jaco_Matrices.Jaco_Mat;
                hess : in Standard_Complex_Hessians.Array_of_Hessians;
                sol : in Standard_Complex_Solutions.Solution ) is

  -- DESCRIPTION :
  --   Evaluates the Hessians at the solution vector in sol,
  --   with the value for sol.t added.

    xt : Standard_Complex_Vectors.Vector(1..sol.n+1);
    eta : double_float;

    use Singular_Values_of_Hessians;

  begin
    xt(sol.v'range) := sol.v;
    xt(xt'last) := sol.t;
    eta := Standard_Distance(jm,hess,xt);
    put("eta : "); put(eta,2); new_line;
  end Standard_Evaluate;

  procedure DoblDobl_Evaluate
              ( jm : in DoblDobl_Complex_Jaco_Matrices.Jaco_Mat;
                hess : in DoblDobl_Complex_Hessians.Array_of_Hessians;
                sol : in DoblDobl_Complex_Solutions.Solution ) is

  -- DESCRIPTION :
  --   Evaluates the Hessians at the solution vector in sol,
  --   with the value for sol.t added.

    xt : DoblDobl_Complex_Vectors.Vector(1..sol.n+1);
    eta : double_double;

    use Singular_Values_of_Hessians;

  begin
    xt(sol.v'range) := sol.v;
    xt(xt'last) := sol.t;
    eta := DoblDobl_Distance(jm,hess,xt);
    put("eta : "); put(eta,2); new_line;
  end DoblDobl_Evaluate;

  procedure QuadDobl_Evaluate
              ( jm : in QuadDobl_Complex_Jaco_Matrices.Jaco_Mat;
                hess : in QuadDobl_Complex_Hessians.Array_of_Hessians;
                sol : in QuadDobl_Complex_Solutions.Solution ) is

  -- DESCRIPTION :
  --   Evaluates the Hessians at the solution vector in sol,
  --   with the value for sol.t added.

    xt : QuadDobl_Complex_Vectors.Vector(1..sol.n+1);
    eta : quad_double;

    use Singular_Values_of_Hessians;

  begin
    xt(sol.v'range) := sol.v;
    xt(xt'last) := sol.t;
    eta := QuadDobl_Distance(jm,hess,xt);
    put("eta : "); put(eta,2); new_line;
  end QuadDobl_Evaluate;

  procedure Standard_Evaluate
              ( jm : in Standard_Complex_Jaco_Matrices.Jaco_Mat;
                hess : in Standard_Complex_Hessians.Array_of_Hessians;
                sols : in Standard_Complex_Solutions.Solution_List ) is

  -- DESCRIPTION :
  --   Evaluates the Hessians at the solutions in sols.

    tmp : Standard_Complex_Solutions.Solution_List := sols;
    ls : Standard_Complex_Solutions.Link_to_Solution;

  begin
    while not Standard_Complex_Solutions.Is_Null(tmp) loop
      ls := Standard_Complex_Solutions.Head_Of(tmp);
      Standard_Evaluate(jm,hess,ls.all);
      tmp := Standard_Complex_Solutions.Tail_Of(tmp);
    end loop;
  end Standard_Evaluate;

  procedure DoblDobl_Evaluate
              ( jm : in DoblDobl_Complex_Jaco_Matrices.Jaco_Mat;
                hess : in DoblDobl_Complex_Hessians.Array_of_Hessians;
                sols : in DoblDobl_Complex_Solutions.Solution_List ) is

  -- DESCRIPTION :
  --   Evaluates the Hessians at the solutions in sols.

    tmp : DoblDobl_Complex_Solutions.Solution_List := sols;
    ls : DoblDobl_Complex_Solutions.Link_to_Solution;

  begin
    while not DoblDobl_Complex_Solutions.Is_Null(tmp) loop
      ls := DoblDobl_Complex_Solutions.Head_Of(tmp);
      DoblDobl_Evaluate(jm,hess,ls.all);
      tmp := DoblDobl_Complex_Solutions.Tail_Of(tmp);
    end loop;
  end DoblDobl_Evaluate;

  procedure QuadDobl_Evaluate
              ( jm : in QuadDobl_Complex_Jaco_Matrices.Jaco_Mat;
                hess : in QuadDobl_Complex_Hessians.Array_of_Hessians;
                sols : in QuadDobl_Complex_Solutions.Solution_List ) is

  -- DESCRIPTION :
  --   Evaluates the Hessians at the solutions in sols.

    tmp : QuadDobl_Complex_Solutions.Solution_List := sols;
    ls : QuadDobl_Complex_Solutions.Link_to_Solution;

  begin
    while not QuadDobl_Complex_Solutions.Is_Null(tmp) loop
      ls := QuadDobl_Complex_Solutions.Head_Of(tmp);
      QuadDobl_Evaluate(jm,hess,ls.all);
      tmp := QuadDobl_Complex_Solutions.Tail_Of(tmp);
    end loop;
  end QuadDobl_Evaluate;

  procedure Standard_Test
              ( nbeq : in integer32;
                sols : in Standard_Complex_Solutions.Solution_List ) is

  -- DESCRIPTION :
  --   For the polynomial system defined in Standard_Homotopy,
  --   computes the array of Hessians.

    p : constant Standard_Complex_Poly_Systems.Poly_Sys(1..nbeq)
      := Standard_Homotopy.Homotopy_System;
    n : constant integer32
      := integer32(Standard_Complex_Polynomials.Number_of_Unknowns(p(p'first)));
    n1 : constant integer32 := n-1; -- minus homotopy parameter
    s : Symbol_Table.Symbol;
    h : constant Standard_Complex_Hessians.Array_of_Hessians(p'range)
      := Standard_Complex_Hessians.Create(p,n);
    jm : Standard_Complex_Jaco_Matrices.Jaco_Mat(1..nbeq,1..n1);

  begin
    for i in 1..nbeq loop
      for j in 1..n1 loop
        jm(i,j) := Standard_Complex_Polynomials.Diff(p(i),j);
      end loop;
    end loop;
    s := (s'range => ' ');
    s(s'first) := 't';
    Symbol_Table.Enlarge(1);
    Symbol_Table.Add(s);
    put_line("The homotopy system : ");
    put(p);
    Standard_Evaluate(jm,h,sols);
  end Standard_Test;

  procedure DoblDobl_Test
              ( nbeq : in integer32;
                sols : in DoblDobl_Complex_Solutions.Solution_List ) is

  -- DESCRIPTION :
  --   For the polynomial system defined in DoblDobl_Homotopy,
  --   computes the array of Hessians.

    p : constant DoblDobl_Complex_Poly_Systems.Poly_Sys(1..nbeq)
      := DoblDobl_Homotopy.Homotopy_System;
    n : constant integer32
      := integer32(DoblDobl_Complex_Polynomials.Number_of_Unknowns(p(p'first)));
    n1 : constant integer32 := n-1; -- minus homotopy parameter
    s : Symbol_Table.Symbol;
    h : constant DoblDobl_Complex_Hessians.Array_of_Hessians(p'range)
      := DoblDobl_Complex_Hessians.Create(p,n);
    jm : DoblDobl_Complex_Jaco_Matrices.Jaco_Mat(1..nbeq,1..n1);

  begin
    for i in 1..nbeq loop
      for j in 1..n1 loop
        jm(i,j) := DoblDobl_Complex_Polynomials.Diff(p(i),j);
      end loop;
    end loop;
    s := (s'range => ' ');
    s(s'first) := 't';
    Symbol_Table.Enlarge(1);
    Symbol_Table.Add(s);
    put_line("The homotopy system : ");
    put(p);
    DoblDobl_Evaluate(jm,h,sols);
  end DoblDobl_Test;

  procedure QuadDobl_Test
              ( nbeq : in integer32;
                sols : in QuadDobl_Complex_Solutions.Solution_List ) is

  -- DESCRIPTION :
  --   For the polynomial system defined in QuadDobl_Homotopy,
  --   computes the array of Hessians.

    p : constant QuadDobl_Complex_Poly_Systems.Poly_Sys(1..nbeq)
      := QuadDobl_Homotopy.Homotopy_System;
    n : constant integer32
      := integer32(QuadDobl_Complex_Polynomials.Number_of_Unknowns(p(p'first)));
    n1 : constant integer32 := n-1; -- minus homotopy parameter
    s : Symbol_Table.Symbol;
    h : constant QuadDobl_Complex_Hessians.Array_of_Hessians(p'range)
      := QuadDobl_Complex_Hessians.Create(p,n);
    jm : QuadDobl_Complex_Jaco_Matrices.Jaco_Mat(1..nbeq,1..n1);

  begin
    for i in 1..nbeq loop
      for j in 1..n1 loop
        jm(i,j) := QuadDobl_Complex_Polynomials.Diff(p(i),j);
      end loop;
    end loop;
    s := (s'range => ' ');
    s(s'first) := 't';
    Symbol_Table.Enlarge(1);
    Symbol_Table.Add(s);
    put_line("The homotopy system : ");
    put(p);
    QuadDobl_Evaluate(jm,h,sols);
  end QuadDobl_Test;

  procedure Standard_Main is

  -- DESCRIPTION :
  --   Prompts the user for a homotopy,
  --   for running a test in double precision.

    nbeq : integer32 := 0;
    sols : Standard_Complex_Solutions.Solution_List;
    len : natural32;

  begin
    new_line;
    put_line("Reading an artificial parameter homotopy ...");
    Homotopy_Series_Readers.Standard_Reader(nbeq,sols);
    new_line;
    put("Number of equations : "); put(nbeq,1); new_line;
    len := Standard_Complex_Solutions.Length_Of(sols);
    put("Number of solutions : "); put(len,1); new_line;
    Standard_Test(nbeq,sols);
  end Standard_Main;

  procedure DoblDobl_Main is

  -- DESCRIPTION :
  --   Prompts the user for a homotopy,
  --   for running a test in double double precision.

    nbeq : integer32 := 0;
    sols : DoblDobl_Complex_Solutions.Solution_List;
    len : natural32;

  begin
    new_line;
    put_line("Reading an artificial parameter homotopy ...");
    Homotopy_Series_Readers.DoblDobl_Reader(nbeq,sols);
    new_line;
    put("Number of equations : "); put(nbeq,1); new_line;
    len := DoblDobl_Complex_Solutions.Length_Of(sols);
    put("Number of solutions : "); put(len,1); new_line;
    DoblDobl_Test(nbeq,sols);
  end DoblDobl_Main;

  procedure QuadDobl_Main is

  -- DESCRIPTION :
  --   Prompts the user for a homotopy,
  --   for running a test in quad double precision.

    nbeq : integer32 := 0;
    sols : QuadDobl_Complex_Solutions.Solution_List;
    len : natural32;

  begin
    new_line;
    put_line("Reading an artificial parameter homotopy ...");
    Homotopy_Series_Readers.QuadDobl_Reader(nbeq,sols);
    new_line;
    put("Number of equations : "); put(nbeq,1); new_line;
    len := QuadDobl_Complex_Solutions.Length_Of(sols);
    put("Number of solutions : "); put(len,1); new_line;
    QuadDobl_Test(nbeq,sols);
  end QuadDobl_Main;

  procedure Main is

  -- DESCRIPTION :
  --   Prompts the user for the precision and then launches
  --   the corresponding test.

    ans : character;

  begin
    new_line;
    put_line("MENU for the working precision :");
    put_line("  0. standard double precision;");
    put_line("  1. double double precision;");
    put_line("  2. quad double precision.");
    put("Type 0, 1, or 2 to select the precision : ");
    Ask_Alternative(ans,"012");
    case ans is
      when '0' => Standard_Main;
      when '1' => DoblDobl_Main;
      when '2' => QuadDobl_Main;
      when others => null;
    end case;
  end Main;

begin
  Main;
end ts_hesscrit;
