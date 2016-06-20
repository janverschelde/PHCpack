with text_io;                           use text_io;
with Communications_with_User;          use Communications_with_User;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with Standard_Integer_Vectors;
with Symbol_Table;
with Standard_Complex_Poly_Systems;
with Standard_Complex_Poly_Systems_io;
with Standard_Dense_Series_Vectors;
with Standard_Dense_Series_Matrices;
with Standard_Linear_Series_Solvers;    use Standard_Linear_Series_Solvers;
with Standard_Series_Polynomials;
with Series_and_Polynomials;
with Series_and_Polynomials_io;
with Standard_Series_Poly_Systems;
with Standard_Series_Poly_SysFun;
with Standard_Series_Jaco_Matrices;

procedure ts_sersys is

-- DESCRIPTION :
--   Tests the methods on systems of series polynomials.

  procedure Read_Series_Vector
              ( v : out Standard_Dense_Series_Vectors.Link_to_Vector;
                dim,idx : in integer32 ) is

  -- DESCRIPTION :
  --   Prompts the user for a vector of series.
  --   The dim is the number of variables in the system where the
  --   series will be evaluated.
  --   The idx is the index of the variable used as series variable.

    use Standard_Series_Polynomials;

    vx : integer32 := 1;

  begin
    if idx = 0 then
      Symbol_Table.Enlarge(1);
      vx := dim+1;
    end if;
    new_line;
    put("Reading a vector of "); put(dim,1); put_line(" series ...");
    Series_and_Polynomials_io.get(v,vx);
    put_line("The vector of series :");
    Series_and_Polynomials_io.put(v.all);
  end Read_Series_Vector;

  procedure Test_Evaluation
              ( p : in Standard_Series_Poly_Systems.Poly_Sys;
                idx : in integer32 ) is

  -- DESCRIPTION :
  --   Prompts the user for a vector of series
  --   and evaluates the system p.
  --   The idx is the index of the variable used as series variable.

    use Standard_Series_Polynomials;

    n : constant integer32 := integer32(Number_of_Unknowns(p(p'first)));
    x : Standard_Dense_Series_Vectors.Link_to_Vector;
    y : Standard_Dense_Series_Vectors.Vector(p'range);
    order : integer32 := 0;

  begin
    Read_Series_Vector(x,n,idx);
    new_line;
    put("Give the order of series : "); get(order);
    Series_and_Polynomials.Set_Order(x.all,order);
    put_line("Evaluating the series ...");
    y := Standard_Series_Poly_SysFun.Eval(p,x.all);
    put_line("The value of the system at the given series :");
    Series_and_Polynomials_io.put(y);
  end Test_Evaluation;

  procedure Write ( A : Standard_Dense_Series_Matrices.Matrix ) is

  -- DESCRIPTION :
  --   Writes the series in the matrix A to screen.

  begin
    for i in A'range(1) loop
      for j in A'range(2) loop
        put("A["); put(i,1); put(","); put(j,1); put_line("] :");
        Series_and_Polynomials_io.put(A(i,j)); new_line;
      end loop;
    end loop;
  end Write;

  procedure Test_Newton_Step
              ( p : in Standard_Series_Poly_Systems.Poly_Sys;
                idx : in integer32 ) is

  -- DESCRIPTION :
  --   Prompts the user for a vector of series
  --   and evaluates the system p.
  --   The idx is the index of the variable used as series variable.

    use Standard_Series_Polynomials;

    n : constant integer32 := integer32(Number_of_Unknowns(p(p'first)));
    x : Standard_Dense_Series_Vectors.Link_to_Vector;
    dx : Standard_Dense_Series_Vectors.Vector(1..n);
    px : Standard_Dense_Series_Vectors.Vector(p'range);
    jp : Standard_Series_Jaco_Matrices.Jaco_Mat(p'range,1..n)
       := Standard_Series_Jaco_Matrices.Create(p);
    jm : Standard_Dense_Series_Matrices.Matrix(p'range,1..n);
    order : integer32 := 0;
    ipvt : Standard_Integer_Vectors.Vector(1..n);
    info : integer32;

  begin
    Read_Series_Vector(x,n,idx);
    new_line;
    put("Give the order of series : "); get(order);
    Series_and_Polynomials.Set_Order(x.all,order);
    put_line("Evaluating the series ...");
    px := Standard_Series_Poly_SysFun.Eval(p,x.all);
    put_line("The value of the system at the given series :");
    Series_and_Polynomials_io.put(px);
    Standard_Dense_Series_Vectors.Min(px);
    Series_and_Polynomials.Set_Order(px,order);
    jm := Standard_Series_Jaco_Matrices.Eval(jp,x.all);
    Series_and_Polynomials.Set_Order(jm,order);
    put_line("The Jacobian matrix : ");
    Write(jm);
    LUfac(jm,n,ipvt,info);
    if info /= 0 then
      put("info = "); put(info,1); new_line;
    else
      dx := px;
      LUsolve(jm,n,ipvt,dx);
      put_line("The update of the Newton step :");
      Series_and_Polynomials_io.put(dx);
      Standard_Dense_Series_Vectors.Add(x.all,dx);
      put_line("After adding the update to the current series :");
      Series_and_Polynomials_io.put(x.all);
      px := Standard_Series_Poly_SysFun.Eval(p,x.all);
      put_line("After evaluation in the original system :");
      Series_and_Polynomials_io.put(px);
    end if;
  end Test_Newton_Step;

  procedure Main is

  -- DESCRIPTION :
  --   Prompts the user for a system of series polynomials,
  --   displays the menu of test operations, asks for a choice,
  --   and then launches the corresponding test.

    ls : Standard_Series_Poly_Systems.Link_to_Poly_Sys;
    ix : integer32 := 0;
    ans : character;

  begin
    new_line;
    put("Give the index of the series variable : "); get(ix);
    new_line;
    Series_and_Polynomials_io.get(ls,ix);
    new_line;
    put_line("The polynomial system : ");
    Series_and_Polynomials_io.put(ls.all,ix);
    new_line;
    put_line("MENU to test systems of series polynomials :");
    put_line("  0. test evaluation in a given series;");
    put_line("  1. run one Newton step.");
    put("Type 0 or 1 to select a test : ");
    Ask_Alternative(ans,"01");
    case ans is
      when '0' => Test_Evaluation(ls.all,ix);
      when '1' => Test_Newton_Step(ls.all,ix);
      when others => null;
    end case;
  end Main;

begin
  Main;
end ts_sersys;
