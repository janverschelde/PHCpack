with text_io;                           use text_io;
with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;       use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with Standard_Complex_Vectors;
with Standard_Complex_Vectors_io;       use Standard_Complex_Vectors_io;
with Standard_Complex_Matrices;
with Standard_Complex_Singular_Values;
with Standard_Complex_Polynomials;
with Standard_Complex_Poly_Systems;
with Standard_Complex_Poly_Systems_io;  use Standard_Complex_Poly_Systems_io;
with Standard_Complex_Hessians;
with Standard_Complex_Solutions;
with Symbol_Table;
with Standard_Homotopy;
with Homotopy_Series_Readers;

procedure ts_hesscrit is

-- DESCRIPTION :
--   Interactive development testing of the Hessian criterion
--   for solution curves defined by polynomial homotopies.

  procedure Singular_Values
              ( A : in out Standard_Complex_Matrices.Matrix ) is

  -- DESCRIPTION :
  --   Computes the singular values of A.

    n : constant integer32 := A'last(1);
    p : constant integer32 := A'last(2);
    mm : constant integer32 := Standard_Complex_Singular_Values.Min0(n+1,p);
    s : Standard_Complex_Vectors.Vector(1..mm);
    e : Standard_Complex_Vectors.Vector(1..p);
    u : Standard_Complex_Matrices.Matrix(1..n,1..n);
    v : Standard_Complex_Matrices.Matrix(1..p,1..p);
    job : constant integer32 := 11;
    info : integer32;

  begin
    Standard_Complex_Singular_Values.SVD(A,n,p,s,e,u,v,job,info);
    put_line("The singular values : "); put_line(s);
  end Singular_Values;

  procedure Standard_Evaluate
              ( hess : in Standard_Complex_Hessians.Array_of_Hessians;
                sol : in Standard_Complex_Solutions.Solution ) is

  -- DESCRIPTION :
  --   Evaluates the Hessians at the solution vector in sol,
  --   with the value for sol.t added.

    xt : Standard_Complex_Vectors.Vector(1..sol.n+1);

  begin
    xt(sol.v'range) := sol.v;
    xt(xt'last) := sol.t;
    for i in hess'range loop
      declare
        hs : constant Standard_Complex_Hessians.Link_to_hessian := hess(i);
        mh : Standard_Complex_Matrices.Matrix(hs'range(1),hs'range(2));
      begin
        mh := Standard_Complex_Hessians.Eval(hs,xt);
        Singular_Values(mh);
      end;
    end loop;
  end Standard_Evaluate;

  procedure Standard_Evaluate
              ( hess : in Standard_Complex_Hessians.Array_of_Hessians;
                sols : in Standard_Complex_Solutions.Solution_List ) is

  -- DESCRIPTION :
  --   Evaluates the Hessians at the solutions in sols.

    tmp : Standard_Complex_Solutions.Solution_List := sols;
    ls : Standard_Complex_Solutions.Link_to_Solution;

  begin
    while not Standard_Complex_Solutions.Is_Null(tmp) loop
      ls := Standard_Complex_Solutions.Head_Of(tmp);
      Standard_Evaluate(hess,ls.all);
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
    s : Symbol_Table.Symbol;
    h : Standard_Complex_Hessians.Array_of_Hessians(p'range)
      := Standard_Complex_Hessians.Create(p,n);

  begin
    s := (s'range => ' ');
    s(s'first) := 't';
    Symbol_Table.Enlarge(1);
    Symbol_Table.Add(s);
    put_line("The homotopy system : ");
    put(p);
    Standard_Evaluate(h,sols);
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
