with text_io;                           use text_io;
with Communications_with_User;          use Communications_with_User;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;         use Standard_Floating_Numbers;
with Standard_Integer_Vectors;
with Symbol_Table;
with Standard_Complex_Series;
with Standard_Complex_Series_Vectors;
with Standard_Complex_Series_Matrices;
with Standard_Series_Linear_Solvers;    use Standard_Series_Linear_Solvers;
with Standard_CSeries_Polynomials;
with Standard_CSeries_Poly_Systems;
with Standard_CSeries_Poly_SysFun;
with Standard_CSeries_Jaco_Matrices;
with DoblDobl_Complex_Series;
with DoblDobl_Complex_Series_Vectors;
with DoblDobl_Complex_Series_Matrices;
with DoblDobl_Series_Linear_Solvers;    use DoblDobl_Series_Linear_Solvers;
with DoblDobl_CSeries_Polynomials;
with DoblDobl_CSeries_Poly_Systems;
with DoblDobl_CSeries_Poly_SysFun;
with DoblDobl_CSeries_Jaco_Matrices;
with QuadDobl_Complex_Series;
with QuadDobl_Complex_Series_Vectors;
with QuadDobl_Complex_Series_Matrices;
with QuadDobl_Series_Linear_Solvers;    use QuadDobl_Series_Linear_Solvers;
with QuadDobl_CSeries_Polynomials;
with QuadDobl_CSeries_Poly_Systems;
with QuadDobl_CSeries_Poly_SysFun;
with QuadDobl_CSeries_Jaco_Matrices;
with Complex_Series_and_Polynomials;
with Complex_Series_and_Polynomials_io;

procedure ts_csersys is

-- DESCRIPTION :
--   Tests the methods on systems of series polynomials.

  procedure Read_Series_Vector
              ( v : out Standard_Complex_Series_Vectors.Link_to_Vector;
                dim,idx : in integer32 ) is

  -- DESCRIPTION :
  --   Prompts the user for a vector of series, in standard double precision.
  --   The dim is the number of variables in the system where the
  --   series will be evaluated.
  --   The idx is the index of the variable used as series variable.

    use Standard_CSeries_Polynomials;

    vx : integer32 := 1;

  begin
    if idx = 0 then
      Symbol_Table.Enlarge(1);
      vx := dim+1;
    end if;
    new_line;
    put("Reading a vector of "); put(dim,1); put_line(" series ...");
    Complex_Series_and_Polynomials_io.get(v,vx);
    put_line("The vector of series :");
    Complex_Series_and_Polynomials_io.put(v.all);
  end Read_Series_Vector;

  procedure Read_Series_Vector
              ( v : out DoblDobl_Complex_Series_Vectors.Link_to_Vector;
                dim,idx : in integer32 ) is

  -- DESCRIPTION :
  --   Prompts the user for a vector of series, in double double precision.
  --   The dim is the number of variables in the system where the
  --   series will be evaluated.
  --   The idx is the index of the variable used as series variable.

    use DoblDobl_CSeries_Polynomials;

    vx : integer32 := 1;

  begin
    if idx = 0 then
      Symbol_Table.Enlarge(1);
      vx := dim+1;
    end if;
    new_line;
    put("Reading a vector of "); put(dim,1); put_line(" series ...");
    Complex_Series_and_Polynomials_io.get(v,vx);
    put_line("The vector of series :");
    Complex_Series_and_Polynomials_io.put(v.all);
  end Read_Series_Vector;

  procedure Read_Series_Vector
              ( v : out QuadDobl_Complex_Series_Vectors.Link_to_Vector;
                dim,idx : in integer32 ) is

  -- DESCRIPTION :
  --   Prompts the user for a vector of series, in quad double precision.
  --   The dim is the number of variables in the system where the
  --   series will be evaluated.
  --   The idx is the index of the variable used as series variable.

    use QuadDobl_CSeries_Polynomials;

    vx : integer32 := 1;

  begin
    if idx = 0 then
      Symbol_Table.Enlarge(1);
      vx := dim+1;
    end if;
    new_line;
    put("Reading a vector of "); put(dim,1); put_line(" series ...");
    Complex_Series_and_Polynomials_io.get(v,vx);
    put_line("The vector of series :");
    Complex_Series_and_Polynomials_io.put(v.all);
  end Read_Series_Vector;

  procedure Test_Evaluation
              ( p : in Standard_CSeries_Poly_Systems.Poly_Sys;
                idx : in integer32 ) is

  -- DESCRIPTION :
  --   Prompts the user for a vector of series
  --   and evaluates the system p, in standard double precision.
  --   The idx is the index of the variable used as series variable.

    use Standard_CSeries_Polynomials;

    n : constant integer32 := integer32(Number_of_Unknowns(p(p'first)));
    x : Standard_Complex_Series_Vectors.Link_to_Vector;
    y : Standard_Complex_Series_Vectors.Vector(p'range);
    degree : integer32 := 0;

  begin
    Read_Series_Vector(x,n,idx);
    new_line;
    put("Give the degree of series : "); get(degree);
    Complex_Series_and_Polynomials.Set_Degree(x.all,degree);
    put_line("Evaluating the series ...");
    y := Standard_CSeries_Poly_SysFun.Eval(p,x.all);
    put_line("The value of the system at the given series :");
    Complex_Series_and_Polynomials_io.put(y);
  end Test_Evaluation;

  procedure Test_Evaluation
              ( p : in DoblDobl_CSeries_Poly_Systems.Poly_Sys;
                idx : in integer32 ) is

  -- DESCRIPTION :
  --   Prompts the user for a vector of series
  --   and evaluates the system p, in double double precision.
  --   The idx is the index of the variable used as series variable.

    use DoblDobl_CSeries_Polynomials;

    n : constant integer32 := integer32(Number_of_Unknowns(p(p'first)));
    x : DoblDobl_Complex_Series_Vectors.Link_to_Vector;
    y : DoblDobl_Complex_Series_Vectors.Vector(p'range);
    degree : integer32 := 0;

  begin
    Read_Series_Vector(x,n,idx);
    new_line;
    put("Give the degree of series : "); get(degree);
    Complex_Series_and_Polynomials.Set_Degree(x.all,degree);
    put_line("Evaluating the series ...");
    y := DoblDobl_CSeries_Poly_SysFun.Eval(p,x.all);
    put_line("The value of the system at the given series :");
    Complex_Series_and_Polynomials_io.put(y);
  end Test_Evaluation;

  procedure Test_Evaluation
              ( p : in QuadDobl_CSeries_Poly_Systems.Poly_Sys;
                idx : in integer32 ) is

  -- DESCRIPTION :
  --   Prompts the user for a vector of series
  --   and evaluates the system p, in double double precision.
  --   The idx is the index of the variable used as series variable.

    use QuadDobl_CSeries_Polynomials;

    n : constant integer32 := integer32(Number_of_Unknowns(p(p'first)));
    x : QuadDobl_Complex_Series_Vectors.Link_to_Vector;
    y : QuadDobl_Complex_Series_Vectors.Vector(p'range);
    degree : integer32 := 0;

  begin
    Read_Series_Vector(x,n,idx);
    new_line;
    put("Give the degree of series : "); get(degree);
    Complex_Series_and_Polynomials.Set_Degree(x.all,degree);
    put_line("Evaluating the series ...");
    y := QuadDobl_CSeries_Poly_SysFun.Eval(p,x.all);
    put_line("The value of the system at the given series :");
    Complex_Series_and_Polynomials_io.put(y);
  end Test_Evaluation;

  procedure Write ( A : Standard_Complex_Series_Matrices.Matrix ) is

  -- DESCRIPTION :
  --   Writes the series in the matrix A to screen.

    use Standard_Complex_Series;

  begin
    for i in A'range(1) loop
      for j in A'range(2) loop
        put("A["); put(i,1); put(","); put(j,1); put_line("] :");
        if A(i,j) /= null then
          Complex_Series_and_Polynomials_io.put(A(i,j).all);
          new_line;
        end if;
      end loop;
    end loop;
  end Write;

  procedure Write ( A : DoblDobl_Complex_Series_Matrices.Matrix ) is

  -- DESCRIPTION :
  --   Writes the series in the matrix A to screen.

    use DoblDobl_Complex_Series;

  begin
    for i in A'range(1) loop
      for j in A'range(2) loop
        put("A["); put(i,1); put(","); put(j,1); put_line("] :");
        if A(i,j) /= null then
          Complex_Series_and_Polynomials_io.put(A(i,j).all);
          new_line;
        end if;
      end loop;
    end loop;
  end Write;

  procedure Write ( A : QuadDobl_Complex_Series_Matrices.Matrix ) is

  -- DESCRIPTION :
  --   Writes the series in the matrix A to screen.

    use QuadDobl_Complex_Series;

  begin
    for i in A'range(1) loop
      for j in A'range(2) loop
        put("A["); put(i,1); put(","); put(j,1); put_line("] :");
        if A(i,j) /= null then
          Complex_Series_and_Polynomials_io.put(A(i,j).all);
          new_line;
        end if;
      end loop;
    end loop;
  end Write;

  procedure Test_Newton_Step
              ( p : in Standard_CSeries_Poly_Systems.Poly_Sys;
                idx : in integer32 ) is

  -- DESCRIPTION :
  --   Prompts the user for a vector of series and applies one
  --   Newton step on p, in standard double precision.
  --   The idx is the index of the variable used as series variable.

    use Standard_CSeries_Polynomials;

    n : constant integer32 := integer32(Number_of_Unknowns(p(p'first)));
    x : Standard_Complex_Series_Vectors.Link_to_Vector;
    dx : Standard_Complex_Series_Vectors.Vector(1..n);
    px : Standard_Complex_Series_Vectors.Vector(p'range);
    jp : constant Standard_CSeries_Jaco_Matrices.Jaco_Mat(p'range,1..n)
       := Standard_CSeries_Jaco_Matrices.Create(p);
    jm : Standard_Complex_Series_Matrices.Matrix(p'range,1..n);
    degree : integer32 := 0;
    ipvt : Standard_Integer_Vectors.Vector(1..n);
    info : integer32;
    tol : constant double_float := 1.0e-12;

  begin
    Read_Series_Vector(x,n,idx);
    new_line;
    put("Give the degree of series : "); get(degree);
    Complex_Series_and_Polynomials.Set_Degree(x.all,degree);
    put_line("Evaluating the series ...");
    px := Standard_CSeries_Poly_SysFun.Eval(p,x.all);
    put_line("The value of the system at the given series :");
    Complex_Series_and_Polynomials_io.put(px);
    Standard_Complex_Series_Vectors.Min(px);
    Complex_Series_and_Polynomials.Set_Degree(px,degree);
    jm := Standard_CSeries_Jaco_Matrices.Eval(jp,x.all);
    Complex_Series_and_Polynomials.Set_Degree(jm,degree);
    put_line("The Jacobian matrix : ");
    Write(jm);
    LUfac(jm,n,ipvt,info);
    if info /= 0 then
      put("info = "); put(info,1); new_line;
    else
      dx := px;
      LUsolve(jm,n,ipvt,dx);
      put_line("The update of the Newton step :");
      Complex_Series_and_Polynomials.Filter(dx,tol);
      Complex_Series_and_Polynomials_io.put(dx);
      Standard_Complex_Series_Vectors.Add(x.all,dx);
      put_line("After adding the update to the current series :");
      Complex_Series_and_Polynomials.Filter(x.all,tol);
      Complex_Series_and_Polynomials_io.put(x.all);
      px := Standard_CSeries_Poly_SysFun.Eval(p,x.all);
      put_line("After evaluation in the original system :");
      Complex_Series_and_Polynomials.Filter(px,tol);
      Complex_Series_and_Polynomials_io.put(px);
    end if;
  end Test_Newton_Step;

  procedure Test_Newton_Step
              ( p : in DoblDobl_CSeries_Poly_Systems.Poly_Sys;
                idx : in integer32 ) is

  -- DESCRIPTION :
  --   Prompts the user for a vector of series and applies one
  --   Newton step on p, in double double precision.
  --   The idx is the index of the variable used as series variable.

    use DoblDobl_CSeries_Polynomials;

    n : constant integer32 := integer32(Number_of_Unknowns(p(p'first)));
    x : DoblDobl_Complex_Series_Vectors.Link_to_Vector;
    dx : DoblDobl_Complex_Series_Vectors.Vector(1..n);
    px : DoblDobl_Complex_Series_Vectors.Vector(p'range);
    jp : constant DoblDobl_CSeries_Jaco_Matrices.Jaco_Mat(p'range,1..n)
       := DoblDobl_CSeries_Jaco_Matrices.Create(p);
    jm : DoblDobl_Complex_Series_Matrices.Matrix(p'range,1..n);
    degree : integer32 := 0;
    ipvt : Standard_Integer_Vectors.Vector(1..n);
    info : integer32;
    tol : constant double_float := 1.0e-20;

  begin
    Read_Series_Vector(x,n,idx);
    new_line;
    put("Give the degree of series : "); get(degree);
    Complex_Series_and_Polynomials.Set_Degree(x.all,degree);
    put_line("Evaluating the series ...");
    px := DoblDobl_CSeries_Poly_SysFun.Eval(p,x.all);
    put_line("The value of the system at the given series :");
    Complex_Series_and_Polynomials_io.put(px);
    DoblDobl_Complex_Series_Vectors.Min(px);
    Complex_Series_and_Polynomials.Set_Degree(px,degree);
    jm := DoblDobl_CSeries_Jaco_Matrices.Eval(jp,x.all);
    Complex_Series_and_Polynomials.Set_Degree(jm,degree);
    put_line("The Jacobian matrix : ");
    Write(jm);
    LUfac(jm,n,ipvt,info);
    if info /= 0 then
      put("info = "); put(info,1); new_line;
    else
      dx := px;
      LUsolve(jm,n,ipvt,dx);
      put_line("The update of the Newton step :");
      Complex_Series_and_Polynomials.Filter(dx,tol);
      Complex_Series_and_Polynomials_io.put(dx);
      DoblDobl_Complex_Series_Vectors.Add(x.all,dx);
      put_line("After adding the update to the current series :");
      Complex_Series_and_Polynomials.Filter(x.all,tol);
      Complex_Series_and_Polynomials_io.put(x.all);
      px := DoblDobl_CSeries_Poly_SysFun.Eval(p,x.all);
      put_line("After evaluation in the original system :");
      Complex_Series_and_Polynomials.Filter(px,tol);
      Complex_Series_and_Polynomials_io.put(px);
    end if;
  end Test_Newton_Step;

  procedure Test_Newton_Step
              ( p : in QuadDobl_CSeries_Poly_Systems.Poly_Sys;
                idx : in integer32 ) is

  -- DESCRIPTION :
  --   Prompts the user for a vector of series and applies one
  --   Newton step on p, in quad double precision.
  --   The idx is the index of the variable used as series variable.

    use QuadDobl_CSeries_Polynomials;

    n : constant integer32 := integer32(Number_of_Unknowns(p(p'first)));
    x : QuadDobl_Complex_Series_Vectors.Link_to_Vector;
    dx : QuadDobl_Complex_Series_Vectors.Vector(1..n);
    px : QuadDobl_Complex_Series_Vectors.Vector(p'range);
    jp : constant QuadDobl_CSeries_Jaco_Matrices.Jaco_Mat(p'range,1..n)
       := QuadDobl_CSeries_Jaco_Matrices.Create(p);
    jm : QuadDobl_Complex_Series_Matrices.Matrix(p'range,1..n);
    degree : integer32 := 0;
    ipvt : Standard_Integer_Vectors.Vector(1..n);
    info : integer32;
    tol : constant double_float := 1.0e-20;

  begin
    Read_Series_Vector(x,n,idx);
    new_line;
    put("Give the degree of series : "); get(degree);
    Complex_Series_and_Polynomials.Set_Degree(x.all,degree);
    put_line("Evaluating the series ...");
    px := QuadDobl_CSeries_Poly_SysFun.Eval(p,x.all);
    put_line("The value of the system at the given series :");
    Complex_Series_and_Polynomials_io.put(px);
    QuadDobl_Complex_Series_Vectors.Min(px);
    Complex_Series_and_Polynomials.Set_Degree(px,degree);
    jm := QuadDobl_CSeries_Jaco_Matrices.Eval(jp,x.all);
    Complex_Series_and_Polynomials.Set_Degree(jm,degree);
    put_line("The Jacobian matrix : ");
    Write(jm);
    LUfac(jm,n,ipvt,info);
    if info /= 0 then
      put("info = "); put(info,1); new_line;
    else
      dx := px;
      LUsolve(jm,n,ipvt,dx);
      put_line("The update of the Newton step :");
      Complex_Series_and_Polynomials.Filter(dx,tol);
      Complex_Series_and_Polynomials_io.put(dx);
      QuadDobl_Complex_Series_Vectors.Add(x.all,dx);
      put_line("After adding the update to the current series :");
      Complex_Series_and_Polynomials.Filter(x.all,tol);
      Complex_Series_and_Polynomials_io.put(x.all);
      px := QuadDobl_CSeries_Poly_SysFun.Eval(p,x.all);
      put_line("After evaluation in the original system :");
      Complex_Series_and_Polynomials.Filter(px,tol);
      Complex_Series_and_Polynomials_io.put(px);
    end if;
  end Test_Newton_Step;

  procedure Standard_Test is

  -- DESCRIPTION :
  --   Prompts the user for a system of series polynomials,
  --   in standard double precision,
  --   displays the menu of test operations, asks for a choice,
  --   and then launches the corresponding test.

    ls : Standard_CSeries_Poly_Systems.Link_to_Poly_Sys;
    ix : integer32 := 0;
    ans : character;

  begin
    new_line;
    put("Give the index of the series variable : "); get(ix);
    new_line;
    Complex_Series_and_Polynomials_io.get(ls,ix);
    new_line;
    put_line("The polynomial system : ");
    Complex_Series_and_Polynomials_io.put(ls.all,ix);
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
  end Standard_Test;

  procedure DoblDobl_Test is

  -- DESCRIPTION :
  --   Prompts the user for a system of series polynomials,
  --   in double double precision,
  --   displays the menu of test operations, asks for a choice,
  --   and then launches the corresponding test.

    ls : DoblDobl_CSeries_Poly_Systems.Link_to_Poly_Sys;
    ix : integer32 := 0;
    ans : character;

  begin
    new_line;
    put("Give the index of the series variable : "); get(ix);
    new_line;
    Complex_Series_and_Polynomials_io.get(ls,ix);
    new_line;
    put_line("The polynomial system : ");
    Complex_Series_and_Polynomials_io.put(ls.all,ix);
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
  end DoblDobl_Test;

  procedure QuadDobl_Test is

  -- DESCRIPTION :
  --   Prompts the user for a system of series polynomials,
  --   in quad double precision,
  --   displays the menu of test operations, asks for a choice,
  --   and then launches the corresponding test.

    ls : QuadDobl_CSeries_Poly_Systems.Link_to_Poly_Sys;
    ix : integer32 := 0;
    ans : character;

  begin
    new_line;
    put("Give the index of the series variable : "); get(ix);
    new_line;
    Complex_Series_and_Polynomials_io.get(ls,ix);
    new_line;
    put_line("The polynomial system : ");
    Complex_Series_and_Polynomials_io.put(ls.all,ix);
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
  end QuadDobl_Test;

  procedure Main is

  -- DESCRIPTION :
  --   Displays the menu for the working precision and
  --   calls the test corresponding to the choice made.

    ans : character;

  begin
    new_line;
    put_line("MENU for the working precision :");
    put_line("  0. standard double precision;");
    put_line("  1. double double precision; or");
    put_line("  2. quad double precision.");
    put("Type 0, 1, or 2 to select the working precision : ");
    Ask_Alternative(ans,"012");
    case ans is
      when '0' => Standard_Test;
      when '1' => DoblDobl_Test;
      when '2' => QuadDobl_Test;
      when others => null;
    end case;
  end Main;

begin
  Main;
end ts_csersys;
