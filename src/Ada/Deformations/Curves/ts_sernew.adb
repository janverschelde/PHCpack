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
with Standard_Newton_Series;
with Standard_Newton_Matrix_Series;
with DoblDobl_Dense_Series_Vectors;
with DoblDobl_Series_Polynomials;
with DoblDobl_Series_Poly_Systems;
with DoblDobl_Series_Poly_SysFun;
with DoblDobl_Newton_Series;
with DoblDobl_Newton_Matrix_Series;
with QuadDobl_Dense_Series_Vectors;
with QuadDobl_Series_Polynomials;
with QuadDobl_Series_Poly_Systems;
with QuadDobl_Series_Poly_SysFun;
with QuadDobl_Newton_Series;
with QuadDobl_Newton_Matrix_Series;

procedure ts_sernew is

-- DESCRIPTION :
--   Test on the application of Newton's method to computed series solutions.

  procedure Read_Series_Vector
              ( v : out Standard_Dense_Series_Vectors.Link_to_Vector;
                dim,idx : in integer32 ) is

  -- DESCRIPTION :
  --   Prompts the user for a vector of series.
  --   The dim is the number of variables in the system where the
  --   series will be evaluated, in standard double precision.
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

  procedure Read_Series_Vector
              ( v : out DoblDobl_Dense_Series_Vectors.Link_to_Vector;
                dim,idx : in integer32 ) is

  -- DESCRIPTION :
  --   Prompts the user for a vector of series.
  --   The dim is the number of variables in the system where the
  --   series will be evaluated, in double double precision.
  --   The idx is the index of the variable used as series variable.

    use DoblDobl_Series_Polynomials;

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

  procedure Read_Series_Vector
              ( v : out QuadDobl_Dense_Series_Vectors.Link_to_Vector;
                dim,idx : in integer32 ) is

  -- DESCRIPTION :
  --   Prompts the user for a vector of series.
  --   The dim is the number of variables in the system where the
  --   series will be evaluated, in quad double precision.
  --   The idx is the index of the variable used as series variable.

    use QuadDobl_Series_Polynomials;

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
                degree : in integer32;
                x : in Standard_Dense_Series_Vectors.Vector ) is

  -- DESCRIPTION :
  --  Does one step with Newton's method on the system p,
  --  calculating with series x of the given degree,
  --  in standard double precision.

    info : integer32;
    tol : constant double_float := 1.0E-12;
    eva : Standard_Dense_Series_Vectors.Vector(p'range);
    ans : character;
    otp : boolean;
    z : Standard_Dense_Series_Vectors.Vector(x'range);

  begin
    new_line;
    put("Do you want intermediate output ? (y/n) ");
    Ask_Yes_or_No(ans); otp := (ans = 'y');
    new_line;
    put("Run on series matrices ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans = 'y' then
      z := x;
      Series_and_Polynomials.Set_Degree(z,degree);
      if otp then
        Standard_Newton_Series.LU_Newton_Step(standard_output,p,degree,z,info);
      else
        Standard_Newton_Series.LU_Newton_Step(p,degree,z,info);
      end if;
      if info /= 0 then
        put("info = "); put(info,1); new_line;
      else
        Series_and_Polynomials.Filter(z,tol);
        put_line("The updated power series solution :");
        Series_and_Polynomials_io.put(z);
        eva := Standard_Series_Poly_SysFun.Eval(p,z);
        Series_and_Polynomials.Filter(eva,tol);
        put_line("The evaluated solution :");
        Series_and_Polynomials_io.put(eva);
      end if;
    end if;
    new_line;
    put("Run on matrix series ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans = 'y' then
      z := x;
      Series_and_Polynomials.Set_Degree(z,degree);
      if otp then
        Standard_Newton_Matrix_Series.LU_Newton_Step
          (standard_output,p,degree,z,info);
      else
        Standard_Newton_Matrix_Series.LU_Newton_Step(p,degree,z,info);
      end if;
      if info /= 0 then
        put("info = "); put(info,1); new_line;
      else
        Series_and_Polynomials.Filter(z,tol);
        put_line("The updated power series solution :");
        Series_and_Polynomials_io.put(z);
        eva := Standard_Series_Poly_SysFun.Eval(p,z);
        Series_and_Polynomials.Filter(eva,tol);
        put_line("The evaluated solution :");
        Series_and_Polynomials_io.put(eva);
      end if;
    end if;
  end Test_LU_Newton_Step;

  procedure Test_LU_Newton_Step
              ( p : in DoblDobl_Series_Poly_Systems.Poly_Sys;
                degree : in integer32;
                x : in out DoblDobl_Dense_Series_Vectors.Vector ) is

  -- DESCRIPTION :
  --  Does one step with Newton's method on the system p,
  --  calculating with series x of the given degree,
  --  in double double precision.

    info : integer32;
    tol : constant double_float := 1.0E-20;
    eva : DoblDobl_Dense_Series_Vectors.Vector(p'range);
    ans : character;

  begin
    new_line;
    put("Do you want intermediate output ? (y/n) ");
    Ask_Yes_or_No(ans);
    Series_and_Polynomials.Set_Degree(x,degree);
    if ans = 'y' then
      DoblDobl_Newton_Series.LU_Newton_Step(standard_output,p,degree,x,info);
    else
      DoblDobl_Newton_Series.LU_Newton_Step(p,degree,x,info);
    end if;
    if info /= 0 then
      put("info = "); put(info,1); new_line;
    else
      Series_and_Polynomials.Filter(x,tol);
      put_line("The updated power series solution :");
      Series_and_Polynomials_io.put(x);
      eva := DoblDobl_Series_Poly_SysFun.Eval(p,x);
      Series_and_Polynomials.Filter(eva,tol);
      put_line("The evaluated solution :");
      Series_and_Polynomials_io.put(eva);
    end if;
  end Test_LU_Newton_Step;

  procedure Test_LU_Newton_Step
              ( p : in QuadDobl_Series_Poly_Systems.Poly_Sys;
                degree : in integer32;
                x : in out QuadDobl_Dense_Series_Vectors.Vector ) is

  -- DESCRIPTION :
  --  Does one step with Newton's method on the system p,
  --  calculating with series x of the given degree,
  --  in quad double precision.

    info : integer32;
    tol : constant double_float := 1.0E-32;
    eva : QuadDobl_Dense_Series_Vectors.Vector(p'range);
    ans : character;

  begin
    new_line;
    put("Do you want intermediate output ? (y/n) ");
    Ask_Yes_or_No(ans);
    Series_and_Polynomials.Set_Degree(x,degree);
    if ans = 'y' then
      QuadDobl_Newton_Series.LU_Newton_Step(standard_output,p,degree,x,info);
    else
      QuadDobl_Newton_Series.LU_Newton_Step(p,degree,x,info);
    end if;
    if info /= 0 then
      put("info = "); put(info,1); new_line;
    else
      Series_and_Polynomials.Filter(x,tol);
      put_line("The updated power series solution :");
      Series_and_Polynomials_io.put(x);
      eva := QuadDobl_Series_Poly_SysFun.Eval(p,x);
      Series_and_Polynomials.Filter(eva,tol);
      put_line("The evaluated solution :");
      Series_and_Polynomials_io.put(eva);
    end if;
  end Test_LU_Newton_Step;

  procedure LU_Newton_on_Series_Matrices
              ( p : in Standard_Series_Poly_Systems.Poly_Sys;
                otp : in boolean;
                degree : in integer32; nbrit : in integer32;
                x : in Standard_Dense_Series_Vectors.Vector ) is

  -- DESCRIPTION :
  --   Runs Newton's method on matrices with series as entries,
  --   on the system p, starting at x, of the given degree,
  --   with as many iterations as the value of nbrit,
  --   in standard double precision.
  --   If otp, then extra output is written to screen.

    info : integer32;
    tol : constant double_float := 1.0E-12;
    eva : Standard_Dense_Series_Vectors.Vector(p'range);
    deg : integer32 := degree;
    z : Standard_Dense_Series_Vectors.Vector(x'range) := x;

  begin
    Series_and_Polynomials.Set_Degree(z,deg);
    if otp then
      Standard_Newton_Series.LU_Newton_Steps
        (standard_output,p,deg,nbrit,z,info);
    else
      Standard_Newton_Series.LU_Newton_Steps(p,deg,nbrit,z,info);
    end if;
    if info /= 0 then
      put("info = "); put(info,1); new_line;
    else
      Series_and_Polynomials.Filter(z,tol);
      put_line("The updated power series solution :");
      Series_and_Polynomials_io.put(z);
      eva := Standard_Series_Poly_SysFun.Eval(p,z);
      Series_and_Polynomials.Filter(eva,tol);
      put_line("The evaluated solution :");
      Series_and_Polynomials_io.put(eva);
    end if;
  end LU_Newton_on_Series_Matrices;

  procedure LU_Newton_on_Matrix_Series
              ( p : in Standard_Series_Poly_Systems.Poly_Sys;
                otp : in boolean;
                degree : in integer32; nbrit : in integer32;
                x : in Standard_Dense_Series_Vectors.Vector ) is

  -- DESCRIPTION :
  --   Runs Newton's method on matrix series,
  --   on the system p, starting at x, of the given degree,
  --   with as many iterations as the value of nbrit,
  --   in standard double precision.
  --   If otp, then extra output is written to screen.

    info : integer32;
    tol : constant double_float := 1.0E-12;
    eva : Standard_Dense_Series_Vectors.Vector(p'range);
    deg : integer32 := degree;
    z : Standard_Dense_Series_Vectors.Vector(x'range) := x;

  begin
    Series_and_Polynomials.Set_Degree(z,deg);
    if otp then
      Standard_Newton_Matrix_Series.LU_Newton_Steps
        (standard_output,p,deg,nbrit,z,info);
    else
      Standard_Newton_Matrix_Series.LU_Newton_Steps(p,deg,nbrit,z,info);
    end if;
    if info /= 0 then
      put("info = "); put(info,1); new_line;
    else
      Series_and_Polynomials.Filter(z,tol);
      put_line("The updated power series solution :");
      Series_and_Polynomials_io.put(z);
      eva := Standard_Series_Poly_SysFun.Eval(p,z);
      Series_and_Polynomials.Filter(eva,tol);
      put_line("The evaluated solution :");
      Series_and_Polynomials_io.put(eva);
    end if;
  end LU_Newton_on_Matrix_Series;

  procedure Test_LU_Newton_Steps
              ( p : in Standard_Series_Poly_Systems.Poly_Sys;
                degree : in integer32; nbrit : in integer32;
                x : in Standard_Dense_Series_Vectors.Vector ) is

  -- DESCRIPTION :
  --  Does as many steps with Newton's method on the system p,
  --  as the value of nbrit, calculating with series x of the given degree,
  --  in standard double precision.

    ans : character;
    otp : boolean;

  begin
    new_line;
    put("Do you want intermediate output ? (y/n) ");
    Ask_Yes_or_No(ans); otp := (ans = 'y');
    new_line;
    put("Run on series matrices ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans = 'y'
     then LU_Newton_on_Series_Matrices(p,otp,degree,nbrit,x);
    end if;
    new_line;
    put("Run on matrix series ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans = 'y'
     then LU_Newton_on_Matrix_Series(p,otp,degree,nbrit,x);
    end if;
  end Test_LU_Newton_Steps;

  procedure Test_LU_Newton_Steps
              ( p : in DoblDobl_Series_Poly_Systems.Poly_Sys;
                degree : in out integer32; nbrit : in integer32;
                x : in out DoblDobl_Dense_Series_Vectors.Vector ) is

  -- DESCRIPTION :
  --  Does as many steps with Newton's method on the system p,
  --  as the value of nbrit, calculating with series x of the given degree,
  --  in double double precision.

    info : integer32;
    tol : constant double_float := 1.0E-20;
    eva : DoblDobl_Dense_Series_Vectors.Vector(p'range);
    ans : character;

  begin
    new_line;
    put("Do you want intermediate output ? (y/n) ");
    Ask_Yes_or_No(ans);
    Series_and_Polynomials.Set_Degree(x,degree);
    if ans = 'y' then
      DoblDobl_Newton_Series.LU_Newton_Steps
        (standard_output,p,degree,nbrit,x,info);
    else
      DoblDobl_Newton_Series.LU_Newton_Steps(p,degree,nbrit,x,info);
    end if;
    if info /= 0 then
      put("info = "); put(info,1); new_line;
    else
      Series_and_Polynomials.Filter(x,tol);
      put_line("The updated power series solution :");
      Series_and_Polynomials_io.put(x);
      eva := DoblDobl_Series_Poly_SysFun.Eval(p,x);
      Series_and_Polynomials.Filter(eva,tol);
      put_line("The evaluated solution :");
      Series_and_Polynomials_io.put(eva);
    end if;
  end Test_LU_Newton_Steps;

  procedure Test_LU_Newton_Steps
              ( p : in QuadDobl_Series_Poly_Systems.Poly_Sys;
                degree : in out integer32; nbrit : in integer32;
                x : in out QuadDobl_Dense_Series_Vectors.Vector ) is

  -- DESCRIPTION :
  --  Does as many steps with Newton's method on the system p,
  --  as the value of nbrit, calculating with series x of the given degree,
  --  in quad double precision.

    info : integer32;
    tol : constant double_float := 1.0E-20;
    eva : QuadDobl_Dense_Series_Vectors.Vector(p'range);
    ans : character;

  begin
    new_line;
    put("Do you want intermediate output ? (y/n) ");
    Ask_Yes_or_No(ans);
    Series_and_Polynomials.Set_Degree(x,degree);
    if ans = 'y' then
      QuadDobl_Newton_Series.LU_Newton_Steps
        (standard_output,p,degree,nbrit,x,info);
    else
      QuadDobl_Newton_Series.LU_Newton_Steps(p,degree,nbrit,x,info);
    end if;
    if info /= 0 then
      put("info = "); put(info,1); new_line;
    else
      Series_and_Polynomials.Filter(x,tol);
      put_line("The updated power series solution :");
      Series_and_Polynomials_io.put(x);
      eva := QuadDobl_Series_Poly_SysFun.Eval(p,x);
      Series_and_Polynomials.Filter(eva,tol);
      put_line("The evaluated solution :");
      Series_and_Polynomials_io.put(eva);
    end if;
  end Test_LU_Newton_Steps;

  procedure Test_QR_Newton_Step
              ( p : in Standard_Series_Poly_Systems.Poly_Sys;
                degree : in integer32;
                x : in Standard_Dense_Series_Vectors.Vector ) is

  -- DESCRIPTION :
  --  Does one step with Newton's method on the system p,
  --  calculating with series x of the given degree,
  --  in standard double precision.

    info : integer32;
    tol : constant double_float := 1.0E-12;
    eva : Standard_Dense_Series_Vectors.Vector(p'range);
    ans : character;
    otp : boolean;
    z : Standard_Dense_Series_Vectors.Vector(x'range);

  begin
    new_line;
    put("Do you want intermediate output ? (y/n) ");
    Ask_Yes_or_No(ans); otp := (ans = 'y');
    new_line;
    put("Run on series matrices ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans = 'y' then
      z := x;
      Series_and_Polynomials.Set_Degree(z,degree);
      if otp then
        Standard_Newton_Series.QR_Newton_Step
          (standard_output,p,degree,z,info);
      else
        Standard_Newton_Series.QR_Newton_Step(p,degree,z,info);
      end if;
      if info /= 0 then
        put("info = "); put(info,1); new_line;
      else
        Series_and_Polynomials.Filter(z,tol);
        put_line("The updated power series solution :");
        Series_and_Polynomials_io.put(z);
        eva := Standard_Series_Poly_SysFun.Eval(p,z);
        Series_and_Polynomials.Filter(eva,tol);
        put_line("The evaluated solution :");
        Series_and_Polynomials_io.put(eva);
      end if;
    end if;
    new_line;
    put("Run on matrix series ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans = 'y' then
      z := x;
      Series_and_Polynomials.Set_Degree(z,degree);
      if otp then
        Standard_Newton_Matrix_Series.QR_Newton_Step
          (standard_output,p,degree,z,info);
      else
        Standard_Newton_Matrix_Series.QR_Newton_Step(p,degree,z,info);
      end if;
      if info /= 0 then
        put("info = "); put(info,1); new_line;
      else
        Series_and_Polynomials.Filter(z,tol);
        put_line("The updated power series solution :");
        Series_and_Polynomials_io.put(z);
        eva := Standard_Series_Poly_SysFun.Eval(p,z);
        Series_and_Polynomials.Filter(eva,tol);
        put_line("The evaluated solution :");
        Series_and_Polynomials_io.put(eva);
      end if;
    end if;
  end Test_QR_Newton_Step;

  procedure Test_QR_Newton_Step
              ( p : in DoblDobl_Series_Poly_Systems.Poly_Sys;
                degree : in integer32;
                x : in out DoblDobl_Dense_Series_Vectors.Vector ) is

  -- DESCRIPTION :
  --  Does one step with Newton's method on the system p,
  --  calculating with series x of the given degree,
  --  in double double precision.

    info : integer32;
    tol : constant double_float := 1.0E-12;
    eva : DoblDobl_Dense_Series_Vectors.Vector(p'range);
    ans : character;

  begin
    new_line;
    put("Do you want intermediate output ? (y/n) ");
    Ask_Yes_or_No(ans);
    Series_and_Polynomials.Set_Degree(x,degree);
    if ans = 'y' then
      DoblDobl_Newton_Series.QR_Newton_Step(standard_output,p,degree,x,info);
    else
      DoblDobl_Newton_Series.QR_Newton_Step(p,degree,x,info);
    end if;
    if info /= 0 then
      put("info = "); put(info,1); new_line;
    else
      Series_and_Polynomials.Filter(x,tol);
      put_line("The updated power series solution :");
      Series_and_Polynomials_io.put(x);
      eva := DoblDobl_Series_Poly_SysFun.Eval(p,x);
      Series_and_Polynomials.Filter(eva,tol);
      put_line("The evaluated solution :");
      Series_and_Polynomials_io.put(eva);
    end if;
  end Test_QR_Newton_Step;

  procedure Test_QR_Newton_Step
              ( p : in QuadDobl_Series_Poly_Systems.Poly_Sys;
                degree : in integer32;
                x : in out QuadDobl_Dense_Series_Vectors.Vector ) is

  -- DESCRIPTION :
  --  Does one step with Newton's method on the system p,
  --  calculating with series x of the given degree,
  --  in quad double precision.

    info : integer32;
    tol : constant double_float := 1.0E-12;
    eva : QuadDobl_Dense_Series_Vectors.Vector(p'range);
    ans : character;

  begin
    new_line;
    put("Do you want intermediate output ? (y/n) ");
    Ask_Yes_or_No(ans);
    Series_and_Polynomials.Set_Degree(x,degree);
    if ans = 'y' then
      QuadDobl_Newton_Series.QR_Newton_Step(standard_output,p,degree,x,info);
    else
      QuadDobl_Newton_Series.QR_Newton_Step(p,degree,x,info);
    end if;
    if info /= 0 then
      put("info = "); put(info,1); new_line;
    else
      Series_and_Polynomials.Filter(x,tol);
      put_line("The updated power series solution :");
      Series_and_Polynomials_io.put(x);
      eva := QuadDobl_Series_Poly_SysFun.Eval(p,x);
      Series_and_Polynomials.Filter(eva,tol);
      put_line("The evaluated solution :");
      Series_and_Polynomials_io.put(eva);
    end if;
  end Test_QR_Newton_Step;

  procedure Test_QR_Newton_Steps
              ( p : in Standard_Series_Poly_Systems.Poly_Sys;
                degree : in out integer32; nbrit : in integer32;
                x : in out Standard_Dense_Series_Vectors.Vector ) is

  -- DESCRIPTION :
  --  Does as many steps with Newton's method on the system p,
  --  as the value of nbrit, calculating with series x of the given degree,
  --  in standard double precision.

    info : integer32;
    tol : constant double_float := 1.0E-12;
    eva : Standard_Dense_Series_Vectors.Vector(p'range);
    ans : character;

  begin
    new_line;
    put("Do you want intermediate output ? (y/n) ");
    Ask_Yes_or_No(ans);
    Series_and_Polynomials.Set_Degree(x,degree);
    if ans = 'y' then
      Standard_Newton_Series.QR_Newton_Steps
        (standard_output,p,degree,nbrit,x,info);
    else
      Standard_Newton_Series.QR_Newton_Steps(p,degree,nbrit,x,info);
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

  procedure Test_QR_Newton_Steps
              ( p : in DoblDobl_Series_Poly_Systems.Poly_Sys;
                degree : in out integer32; nbrit : in integer32;
                x : in out DoblDobl_Dense_Series_Vectors.Vector ) is

  -- DESCRIPTION :
  --  Does as many steps with Newton's method on the system p,
  --  as the value of nbrit, calculating with series x of the given degree,
  --  in double double precision.

    info : integer32;
    tol : constant double_float := 1.0E-20;
    eva : DoblDobl_Dense_Series_Vectors.Vector(p'range);
    ans : character;

  begin
    new_line;
    put("Do you want intermediate output ? (y/n) ");
    Ask_Yes_or_No(ans);
    Series_and_Polynomials.Set_Degree(x,degree);
    if ans = 'y' then
      DoblDobl_Newton_Series.QR_Newton_Steps
        (standard_output,p,degree,nbrit,x,info);
    else
      DoblDobl_Newton_Series.QR_Newton_Steps(p,degree,nbrit,x,info);
    end if;
    if info /= 0 then
      put("info = "); put(info,1); new_line;
    else
      Series_and_Polynomials.Filter(x,tol);
      put_line("The updated power series solution :");
      Series_and_Polynomials_io.put(x);
      eva := DoblDobl_Series_Poly_SysFun.Eval(p,x);
      Series_and_Polynomials.Filter(eva,tol);
      put_line("The evaluated solution :");
      Series_and_Polynomials_io.put(eva);
    end if;
  end Test_QR_Newton_Steps;

  procedure Test_QR_Newton_Steps
              ( p : in QuadDobl_Series_Poly_Systems.Poly_Sys;
                degree : in out integer32; nbrit : in integer32;
                x : in out QuadDobl_Dense_Series_Vectors.Vector ) is

  -- DESCRIPTION :
  --  Does as many steps with Newton's method on the system p,
  --  as the value of nbrit, calculating with series x of the given degree,
  --  in quad double precision.

    info : integer32;
    tol : constant double_float := 1.0E-20;
    eva : QuadDobl_Dense_Series_Vectors.Vector(p'range);
    ans : character;

  begin
    new_line;
    put("Do you want intermediate output ? (y/n) ");
    Ask_Yes_or_No(ans);
    Series_and_Polynomials.Set_Degree(x,degree);
    if ans = 'y' then
      QuadDobl_Newton_Series.QR_Newton_Steps
        (standard_output,p,degree,nbrit,x,info);
    else
      QuadDobl_Newton_Series.QR_Newton_Steps(p,degree,nbrit,x,info);
    end if;
    if info /= 0 then
      put("info = "); put(info,1); new_line;
    else
      Series_and_Polynomials.Filter(x,tol);
      put_line("The updated power series solution :");
      Series_and_Polynomials_io.put(x);
      eva := QuadDobl_Series_Poly_SysFun.Eval(p,x);
      Series_and_Polynomials.Filter(eva,tol);
      put_line("The evaluated solution :");
      Series_and_Polynomials_io.put(eva);
    end if;
  end Test_QR_Newton_Steps;

  procedure Standard_Test_LU_Newton is

  -- DESCRIPTION :
  --   Prompts the user for a system of series coefficients
  --   and for an initial series approximation.
  --   Computations are done in standard double precision.

    ls : Standard_Series_Poly_Systems.Link_to_Poly_Sys;
    sol : Standard_Dense_Series_Vectors.Link_to_Vector;
    idx,degree,nbr : integer32 := 0;
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
    put("Give the start degree of the computations : "); get(degree);
    new_line;
    put("Give the number of Newton steps : "); get(nbr);
    if nbr = 1
     then Test_LU_Newton_Step(ls.all,degree,sol.all);
     else Test_LU_Newton_Steps(ls.all,degree,nbr,sol.all);
    end if;
  end Standard_Test_LU_Newton;

  procedure DoblDobl_Test_LU_Newton is

  -- DESCRIPTION :
  --   Prompts the user for a system of series coefficients
  --   and for an initial series approximation.
  --   Computations are done in double double precision.

    ls : DoblDobl_Series_Poly_Systems.Link_to_Poly_Sys;
    sol : DoblDobl_Dense_Series_Vectors.Link_to_Vector;
    idx,degree,nbr : integer32 := 0;
    dim : integer32;

    use DoblDobl_Series_Polynomials;

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
    put("Give the start degree of the computations : "); get(degree);
    new_line;
    put("Give the number of Newton steps : "); get(nbr);
    if nbr = 1
     then Test_LU_Newton_Step(ls.all,degree,sol.all);
     else Test_LU_Newton_Steps(ls.all,degree,nbr,sol.all);
    end if;
  end DoblDobl_Test_LU_Newton;

  procedure QuadDobl_Test_LU_Newton is

  -- DESCRIPTION :
  --   Prompts the user for a system of series coefficients
  --   and for an initial series approximation.
  --   Computations are done in double double precision.

    ls : QuadDobl_Series_Poly_Systems.Link_to_Poly_Sys;
    sol : QuadDobl_Dense_Series_Vectors.Link_to_Vector;
    idx,degree,nbr : integer32 := 0;
    dim : integer32;

    use QuadDobl_Series_Polynomials;

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
    put("Give the start degree of the computations : "); get(degree);
    new_line;
    put("Give the number of Newton steps : "); get(nbr);
    if nbr = 1
     then Test_LU_Newton_Step(ls.all,degree,sol.all);
     else Test_LU_Newton_Steps(ls.all,degree,nbr,sol.all);
    end if;
  end QuadDobl_Test_LU_Newton;

  procedure Standard_Test_QR_Newton is

  -- DESCRIPTION :
  --   Prompts the user for a system of series coefficients
  --   and for an initial series approximation.
  --   Then the test on the QR Newton method are done,
  --   in standard double precision.

    ls : Standard_Series_Poly_Systems.Link_to_Poly_Sys;
    sol : Standard_Dense_Series_Vectors.Link_to_Vector;
    idx,degree,nbr : integer32 := 0;
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
    put("Give the start degree of the computations : "); get(degree);
    new_line;
    put("Give the number of Newton steps : "); get(nbr);
    if nbr = 1
     then Test_QR_Newton_Step(ls.all,degree,sol.all);
     else Test_QR_Newton_Steps(ls.all,degree,nbr,sol.all);
    end if;
  end Standard_Test_QR_Newton;

  procedure DoblDobl_Test_QR_Newton is

  -- DESCRIPTION :
  --   Prompts the user for a system of series coefficients
  --   and for an initial series approximation.
  --   Then the test on the QR Newton method are done,
  --   in double double precision.

    ls : DoblDobl_Series_Poly_Systems.Link_to_Poly_Sys;
    sol : DoblDobl_Dense_Series_Vectors.Link_to_Vector;
    idx,degree,nbr : integer32 := 0;
    dim : integer32;

    use DoblDobl_Series_Polynomials;

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
    put("Give the start degree of the computations : "); get(degree);
    new_line;
    put("Give the number of Newton steps : "); get(nbr);
    if nbr = 1
     then Test_QR_Newton_Step(ls.all,degree,sol.all);
     else Test_QR_Newton_Steps(ls.all,degree,nbr,sol.all);
    end if;
  end DoblDobl_Test_QR_Newton;

  procedure QuadDobl_Test_QR_Newton is

  -- DESCRIPTION :
  --   Prompts the user for a system of series coefficients
  --   and for an initial series approximation.
  --   Then the test on the QR Newton method are done,
  --   in quad double precision.

    ls : QuadDobl_Series_Poly_Systems.Link_to_Poly_Sys;
    sol : QuadDobl_Dense_Series_Vectors.Link_to_Vector;
    idx,degree,nbr : integer32 := 0;
    dim : integer32;

    use QuadDobl_Series_Polynomials;

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
    put("Give the start degree of the computations : "); get(degree);
    new_line;
    put("Give the number of Newton steps : "); get(nbr);
    if nbr = 1
     then Test_QR_Newton_Step(ls.all,degree,sol.all);
     else Test_QR_Newton_Steps(ls.all,degree,nbr,sol.all);
    end if;
  end QuadDobl_Test_QR_Newton;

  procedure Main is

  -- DESCRIPTION :
  --   Displays the test menu and then calls the selected test.

    ans,prc : character;

  begin
    new_line;
    put_line("MENU to test Newton's method on truncated power series :");
    put_line("  1. test LU factorization for square systems;");
    put_line("  2. test QR decomposition for overdetermined systems.");
    put("Type 1 or 2 to select the test : ");
    Ask_Alternative(ans,"12");
    new_line;
    put_line("MENU to set the working precision :");
    put_line("  0. standard double precision;");
    put_line("  1. double double precision; or");
    put_line("  2. quad double precision.");
    put("Type 0, 1, or 2 to select the working precision : ");
    Ask_Alternative(prc,"012");
    case ans is
      when '1' =>
        case prc is 
          when '0' => Standard_Test_LU_Newton;
          when '1' => DoblDobl_Test_LU_Newton;
          when '2' => QuadDobl_Test_LU_Newton;
          when others => null;
        end case;
      when '2' =>
        case prc is
          when '0' => Standard_Test_QR_Newton;
          when '1' => DoblDobl_Test_QR_Newton;
          when '2' => QuadDobl_Test_QR_Newton;
          when others => null;
        end case;
      when others => null;
    end case;
  end Main;

begin
  Main;
end ts_sernew;
