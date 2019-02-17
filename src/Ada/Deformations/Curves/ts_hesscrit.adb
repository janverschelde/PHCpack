with text_io;                           use text_io;
with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;       use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;         use Standard_Floating_Numbers;
with Standard_Floating_Numbers_io;      use Standard_Floating_Numbers_io;
with Standard_Complex_Vectors;
with Standard_Complex_Polynomials;
with Standard_Complex_Poly_Systems;
with Standard_Complex_Poly_Systems_io;  use Standard_Complex_Poly_Systems_io;
with Standard_Complex_Jaco_Matrices;
with Standard_Complex_Hessians;
with Standard_Complex_Solutions;
with Symbol_Table;
with Standard_Homotopy;
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

  procedure Main is

  -- DESCRIPTION :
  --   Prompts the user for a homotopy.

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
  end Main;

begin
  Main;
end ts_hesscrit;
