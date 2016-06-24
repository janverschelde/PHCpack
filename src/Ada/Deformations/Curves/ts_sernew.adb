with text_io;                           use text_io;
with Communications_with_User;          use Communications_with_User;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;         use Standard_Floating_Numbers;
with Symbol_Table;
with Standard_Dense_Series_Vectors;
with Standard_Series_Polynomials;
with Series_and_Polynomials;
with Series_and_Polynomials_io;
with Standard_Series_Poly_Systems;
with Standard_Series_Poly_SysFun;
with Standard_Newton_Series;            use Standard_Newton_Series;

procedure ts_sernew is

-- DESCRIPTION :
--   Test on the application of Newton's method to computed series solutions.

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

  procedure Test_LU_Newton_Step
              ( p : in Standard_Series_Poly_Systems.Poly_Sys;
                order : in integer32;
                x : in out Standard_Dense_Series_Vectors.Vector ) is

  -- DESCRIPTION :
  --  Does one step with Newton's method on the system p,
  --  calculating with series x of the given order.

    info : integer32;
    tol : constant double_float := 1.0E-12;
    eva : Standard_Dense_Series_Vectors.Vector(p'range);
    ans : character;

  begin
    new_line;
    put("Do you want intermediate output ? (y/n) ");
    Ask_Yes_or_No(ans);
    Series_and_Polynomials.Set_Order(x,order);
    if ans = 'y'
     then LU_Newton_Step(standard_output,p,order,x,info);
     else LU_Newton_Step(p,order,x,info);
    end if;
    if info /= 0 then
      put("info = "); put(info,1); new_line;
    else
      Series_and_Polynomials.Filter(x,tol);
      put_line("The updated power series solution :");
      Series_and_Polynomials_io.put(x);
      eva := Standard_Series_Poly_SysFun.Eval(p,x);
      Series_and_Polynomials.Filter(eva,tol);
      put_line("The evaluated solution :");
      Series_and_Polynomials_io.put(eva);
    end if;
  end Test_LU_Newton_Step;

  procedure Test_LU_Newton_Steps
              ( p : in Standard_Series_Poly_Systems.Poly_Sys;
                order : in out integer32; nbrit : in integer32;
                x : in out Standard_Dense_Series_Vectors.Vector ) is

  -- DESCRIPTION :
  --  Does as many steps with Newton's method on the system p,
  --  as the value of nbrit, calculating with series x of the given order.

    info : integer32;
    tol : constant double_float := 1.0E-12;
    eva : Standard_Dense_Series_Vectors.Vector(p'range);
    ans : character;

  begin
    new_line;
    put("Do you want intermediate output ? (y/n) ");
    Ask_Yes_or_No(ans);
    Series_and_Polynomials.Set_Order(x,order);
    if ans = 'y'
     then LU_Newton_Steps(standard_output,p,order,nbrit,x,info);
     else LU_Newton_Steps(p,order,nbrit,x,info);
    end if;
    if info /= 0 then
      put("info = "); put(info,1); new_line;
    else
      Series_and_Polynomials.Filter(x,tol);
      put_line("The updated power series solution :");
      Series_and_Polynomials_io.put(x);
      eva := Standard_Series_Poly_SysFun.Eval(p,x);
      Series_and_Polynomials.Filter(eva,tol);
      put_line("The evaluated solution :");
      Series_and_Polynomials_io.put(eva);
    end if;
  end Test_LU_Newton_Steps;

  procedure Test_QR_Newton_Step
              ( p : in Standard_Series_Poly_Systems.Poly_Sys;
                order : in integer32;
                x : in out Standard_Dense_Series_Vectors.Vector ) is

  -- DESCRIPTION :
  --  Does one step with Newton's method on the system p,
  --  calculating with series x of the given order.

    info : integer32;
    tol : constant double_float := 1.0E-12;
    eva : Standard_Dense_Series_Vectors.Vector(p'range);
    ans : character;

  begin
    new_line;
    put("Do you want intermediate output ? (y/n) ");
    Ask_Yes_or_No(ans);
    Series_and_Polynomials.Set_Order(x,order);
    if ans = 'y'
     then QR_Newton_Step(standard_output,p,order,x,info);
     else QR_Newton_Step(p,order,x,info);
    end if;
    if info /= 0 then
      put("info = "); put(info,1); new_line;
    else
      Series_and_Polynomials.Filter(x,tol);
      put_line("The updated power series solution :");
      Series_and_Polynomials_io.put(x);
      eva := Standard_Series_Poly_SysFun.Eval(p,x);
      Series_and_Polynomials.Filter(eva,tol);
      put_line("The evaluated solution :");
      Series_and_Polynomials_io.put(eva);
    end if;
  end Test_QR_Newton_Step;

  procedure Test_QR_Newton_Steps
              ( p : in Standard_Series_Poly_Systems.Poly_Sys;
                order : in out integer32; nbrit : in integer32;
                x : in out Standard_Dense_Series_Vectors.Vector ) is

  -- DESCRIPTION :
  --  Does as many steps with Newton's method on the system p,
  --  as the value of nbrit, calculating with series x of the given order.

    info : integer32;
    tol : constant double_float := 1.0E-12;
    eva : Standard_Dense_Series_Vectors.Vector(p'range);
    ans : character;

  begin
    new_line;
    put("Do you want intermediate output ? (y/n) ");
    Ask_Yes_or_No(ans);
    Series_and_Polynomials.Set_Order(x,order);
    if ans = 'y'
     then QR_Newton_Steps(standard_output,p,order,nbrit,x,info);
     else QR_Newton_Steps(p,order,nbrit,x,info);
    end if;
    if info /= 0 then
      put("info = "); put(info,1); new_line;
    else
      Series_and_Polynomials.Filter(x,tol);
      put_line("The updated power series solution :");
      Series_and_Polynomials_io.put(x);
      eva := Standard_Series_Poly_SysFun.Eval(p,x);
      Series_and_Polynomials.Filter(eva,tol);
      put_line("The evaluated solution :");
      Series_and_Polynomials_io.put(eva);
    end if;
  end Test_QR_Newton_Steps;

  procedure Test_LU_Newton is

  -- DESCRIPTION :
  --   Prompts the user for a system of series coefficients
  --   and for an initial series approximation.

    ls : Standard_Series_Poly_Systems.Link_to_Poly_Sys;
    sol : Standard_Dense_Series_Vectors.Link_to_Vector;
    idx,order,nbr : integer32 := 0;
    dim : integer32;

    use Standard_Series_Polynomials;

  begin
    new_line;
    put("Give the index of the series variable : "); get(idx);
    new_line;
    Series_and_Polynomials_io.get(ls,idx);
    new_line;
    put_line("The polynomial system : ");
    dim := integer32(Number_of_Unknowns(ls(ls'first)));
    Series_and_Polynomials_io.put(ls.all,idx);
    Read_Series_Vector(sol,dim,idx);
    new_line;
    put("Give the start order of the computations : "); get(order);
    new_line;
    put("Give the number of Newton steps : "); get(nbr);
    if nbr = 1
     then Test_LU_Newton_Step(ls.all,order,sol.all);
     else Test_LU_Newton_Steps(ls.all,order,nbr,sol.all);
    end if;
  end Test_LU_Newton;

  procedure Test_QR_Newton is

  -- DESCRIPTION :
  --   Prompts the user for a system of series coefficients
  --   and for an initial series approximation.
  --   Then the test on the QR Newton method are done.

    ls : Standard_Series_Poly_Systems.Link_to_Poly_Sys;
    sol : Standard_Dense_Series_Vectors.Link_to_Vector;
    idx,order,nbr : integer32 := 0;
    dim : integer32;

    use Standard_Series_Polynomials;

  begin
    new_line;
    put("Give the index of the series variable : "); get(idx);
    new_line;
    Series_and_Polynomials_io.get(ls,idx);
    new_line;
    dim := integer32(Number_of_Unknowns(ls(ls'first)));
    put("The number of variables : "); put(dim,1); new_line;
    put_line("The polynomial system : ");
    Series_and_Polynomials_io.put(ls.all,idx);
    Read_Series_Vector(sol,dim,idx);
    new_line;
    put("Give the start order of the computations : "); get(order);
    new_line;
    put("Give the number of Newton steps : "); get(nbr);
    if nbr = 1
     then Test_QR_Newton_Step(ls.all,order,sol.all);
     else Test_QR_Newton_Steps(ls.all,order,nbr,sol.all);
    end if;
  end Test_QR_Newton;

  procedure Main is

  -- DESCRIPTION :
  --   Displays the test menu and then calls the selected test.

    ans : character;

  begin
    new_line;
    put_line("MENU to test Newton's method on truncated power series :");
    put_line("  1. test LU factorization for square systems;");
    put_line("  2. test QR decomposition for overdetermined systems.");
    put("Type 1 or 2 to select the test : ");
    Ask_Alternative(ans,"12");
    case ans is
      when '1' => Test_LU_Newton;
      when '2' => Test_QR_Newton;
      when others => null;
    end case;
  end Main;

begin
  Main;
end ts_sernew;
