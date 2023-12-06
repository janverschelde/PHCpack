with text_io;                           use text_io;
with Communications_with_User;          use Communications_with_User;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;         use Standard_Floating_Numbers;
with Hexa_Double_Numbers;               use Hexa_Double_Numbers;
with Hexa_Double_Numbers_io;            use Hexa_Double_Numbers_io;
with HexaDobl_Complex_Numbers;
with HexaDobl_Complex_Numbers_io;       use HexaDobl_Complex_Numbers_io;
with Symbol_Table;
with Complex_Series_and_Polynomials;
with Complex_Series_and_Polynomials_io;
with HexaDobl_Complex_Series_Vectors;
with HexaDobl_CSeries_Polynomials;
with HexaDobl_CSeries_Poly_Systems;
with HexaDobl_CSeries_Poly_SysFun;
with HexaDobl_Newton_Matrix_Series;

package body Test_HD_Newton_Matrix_Series is

  procedure Read_Series_Vector
              ( v : out HexaDobl_Complex_Series_Vectors.Link_to_Vector;
                dim,idx : in integer32 ) is

  -- DESCRIPTION :
  --   Prompts the user for a vector of series.
  --   The dim is the number of variables in the system where the
  --   series will be evaluated, in double double precision.
  --   The idx is the index of the variable used as series variable.

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

  procedure Test_LU_Newton_Step
              ( p : in HexaDobl_CSeries_Poly_Systems.Poly_Sys;
                degree : in integer32;
                x : in HexaDobl_Complex_Series_Vectors.Vector ) is

  -- DESCRIPTION :
  --  Does one step with Newton's method on the system p,
  --  calculating with series x of the given degree,
  --  in double double precision.

    info : integer32;
    tol : constant double_float := 1.0E-20;
    eva : HexaDobl_Complex_Series_Vectors.Vector(p'range);
    ans : character;
    otp : boolean;
    z : HexaDobl_Complex_Series_Vectors.Vector(x'range);

  begin
    new_line;
    put("Do you want intermediate output ? (y/n) ");
    Ask_Yes_or_No(ans); otp := (ans = 'y');
    new_line;
    z := x;
    Complex_Series_and_Polynomials.Set_Degree(z,degree);
    if otp then
      HexaDobl_Newton_Matrix_Series.LU_Newton_Step
        (standard_output,p,degree,z,info);
    else
      HexaDobl_Newton_Matrix_Series.LU_Newton_Step(p,degree,z,info);
    end if;
    if info /= 0 then
      put("info = "); put(info,1); new_line;
    else
      Complex_Series_and_Polynomials.Filter(z,tol);
      put_line("The updated power series solution :");
      Complex_Series_and_Polynomials_io.put(z);
      eva := HexaDobl_CSeries_Poly_SysFun.Eval(p,z);
      Complex_Series_and_Polynomials.Filter(eva,tol);
      put_line("The evaluated solution :");
      Complex_Series_and_Polynomials_io.put(eva);
    end if;
  end Test_LU_Newton_Step;

  procedure LU_Newton_on_Matrix_Series
              ( p : in HexaDobl_CSeries_Poly_Systems.Poly_Sys;
                otp : in boolean;
                degree,maxdeg : in integer32; nbrit : in integer32;
                x : in HexaDobl_Complex_Series_Vectors.Vector ) is

  -- DESCRIPTION :
  --   Runs Newton's method on matrix series,
  --   on the system p, starting at x, of the given degree,
  --   with as many iterations as the value of nbrit,
  --   in double double precision.
  --   If otp, then extra output is written to screen.

    info : integer32;
    tol : constant double_float := 1.0E-12;
    eva : HexaDobl_Complex_Series_Vectors.Vector(p'range);
    deg : integer32 := degree;
    z : HexaDobl_Complex_Series_Vectors.Vector(x'range) := x;

  begin
    Complex_Series_and_Polynomials.Set_Degree(z,deg);
    if otp then
      HexaDobl_Newton_Matrix_Series.LU_Newton_Steps
        (standard_output,p,deg,maxdeg,nbrit,z,info);
    else
      HexaDobl_Newton_Matrix_Series.LU_Newton_Steps(p,deg,maxdeg,nbrit,z,info);
    end if;
    if info /= 0 then
      put("info = "); put(info,1); new_line;
    else
      Complex_Series_and_Polynomials.Filter(z,tol);
      put_line("The updated power series solution :");
      Complex_Series_and_Polynomials_io.put(z);
      eva := HexaDobl_CSeries_Poly_SysFun.Eval(p,z);
      Complex_Series_and_Polynomials.Filter(eva,tol);
      put_line("The evaluated solution :");
      Complex_Series_and_Polynomials_io.put(eva);
    end if;
  end LU_Newton_on_Matrix_Series;

  procedure Test_LU_Newton_Steps
              ( p : in HexaDobl_CSeries_Poly_Systems.Poly_Sys;
                degree,maxdeg : in integer32; nbrit : in integer32;
                x : in HexaDobl_Complex_Series_Vectors.Vector ) is

  -- DESCRIPTION :
  --  Does as many steps with Newton's method on the system p,
  --  as the value of nbrit, calculating with series x of the given degree,
  --  in double double precision.

    ans : character;
    otp : boolean;

  begin
    new_line;
    put("Do you want intermediate output ? (y/n) ");
    Ask_Yes_or_No(ans); otp := (ans = 'y');
    new_line;
    LU_Newton_on_Matrix_Series(p,otp,degree,maxdeg,nbrit,x);
  end Test_LU_Newton_Steps;

  procedure Test_QR_Newton_Step
              ( p : in HexaDobl_CSeries_Poly_Systems.Poly_Sys;
                degree : in integer32;
                x : in HexaDobl_Complex_Series_Vectors.Vector ) is

  -- DESCRIPTION :
  --  Does one step with Newton's method on the system p,
  --  calculating with series x of the given degree,
  --  in double double precision.

    info : integer32;
    tol : constant double_float := 1.0E-12;
    eva : HexaDobl_Complex_Series_Vectors.Vector(p'range);
    ans : character;
    otp : boolean;
    z : HexaDobl_Complex_Series_Vectors.Vector(x'range);

  begin
    new_line;
    put("Do you want intermediate output ? (y/n) ");
    Ask_Yes_or_No(ans); otp := (ans = 'y');
    new_line;
    z := x;
    Complex_Series_and_Polynomials.Set_Degree(z,degree);
    if otp then
      HexaDobl_Newton_Matrix_Series.QR_Newton_Step
        (standard_output,p,degree,z,info);
    else
      HexaDobl_Newton_Matrix_Series.QR_Newton_Step(p,degree,z,info);
    end if;
    if info /= 0 then
      put("info = "); put(info,1); new_line;
    else
      Complex_Series_and_Polynomials.Filter(z,tol);
      put_line("The updated power series solution :");
      Complex_Series_and_Polynomials_io.put(z);
      eva := HexaDobl_CSeries_Poly_SysFun.Eval(p,z);
      Complex_Series_and_Polynomials.Filter(eva,tol);
      put_line("The evaluated solution :");
      Complex_Series_and_Polynomials_io.put(eva);
    end if;
  end Test_QR_Newton_Step;

  procedure QR_Newton_on_Matrix_Series
              ( p : in HexaDobl_CSeries_Poly_Systems.Poly_Sys;
                otp : in boolean;
                degree,maxdeg : in integer32; nbrit : in integer32;
                x : in HexaDobl_Complex_Series_Vectors.Vector ) is

  -- DESCRIPTION :
  --  Runs Newton's method on matrix series.
  --  Does as many steps with Newton's method on the system p,
  --  as the value of nbrit, calculating with series x of the given degree,
  --  in double double precision.
  --  If otp, then extra output is written to screen.

    info : integer32;
    tol : constant double_float := 1.0E-20;
    eva : HexaDobl_Complex_Series_Vectors.Vector(p'range);
    deg : integer32 := degree;
    z : HexaDobl_Complex_Series_Vectors.Vector(x'range) := x;

  begin
    Complex_Series_and_Polynomials.Set_Degree(z,degree);
    if otp then
      HexaDobl_Newton_Matrix_Series.QR_Newton_Steps
        (standard_output,p,deg,maxdeg,nbrit,z,info);
    else
      HexaDobl_Newton_Matrix_Series.QR_Newton_Steps
        (p,deg,maxdeg,nbrit,z,info);
    end if;
    if info /= 0 then
      put("info = "); put(info,1); new_line;
    else
      Complex_Series_and_Polynomials.Filter(z,tol);
      put_line("The updated power series solution :");
      Complex_Series_and_Polynomials_io.put(z);
      eva := HexaDobl_CSeries_Poly_SysFun.Eval(p,z);
      Complex_Series_and_Polynomials.Filter(eva,tol);
      put_line("The evaluated solution :");
      Complex_Series_and_Polynomials_io.put(eva);
    end if;
  end QR_Newton_on_Matrix_Series;

  procedure Test_QR_Newton_Steps
              ( p : in HexaDobl_CSeries_Poly_Systems.Poly_Sys;
                degree,maxdeg : in integer32; nbrit : in integer32;
                x : in HexaDobl_Complex_Series_Vectors.Vector ) is

  -- DESCRIPTION :
  --  Does as many steps with Newton's method on the system p,
  --  as the value of nbrit, calculating with series x of the given degree,
  --  in double double precision.

    ans : character;
    otp : boolean;

  begin
    new_line;
    put("Do you want intermediate output ? (y/n) ");
    Ask_Yes_or_No(ans); otp := (ans = 'y');
    new_line;
    QR_Newton_on_Matrix_Series(p,otp,degree,maxdeg,nbrit,x);
  end Test_QR_Newton_Steps;

  procedure Test_SVD_Newton_Step
              ( p : in HexaDobl_CSeries_Poly_Systems.Poly_Sys;
                degree : in integer32; 
                x : in HexaDobl_Complex_Series_Vectors.Vector ) is

  -- DESCRIPTION :
  --  Does one step with Newton's method on the system p,
  --  calculating with series x of the given degree,
  --  in double double precision.

    ans : character;
    otp : boolean;
    info : integer32;
    tol : constant double_float := 1.0E-12;
    eva : HexaDobl_Complex_Series_Vectors.Vector(p'range);
    deg : constant integer32 := degree;
    z : HexaDobl_Complex_Series_Vectors.Vector(x'range) := x;
    rcond : hexa_double;
    one : constant hexa_double := create(1.0);

  begin
    new_line;
    put("Do you want intermediate output ? (y/n) ");
    Ask_Yes_or_No(ans); otp := (ans = 'y');
    new_line;
    Complex_Series_and_Polynomials.Set_Degree(z,degree);
    if otp then
      HexaDobl_Newton_Matrix_Series.SVD_Newton_Step
        (standard_output,p,deg,z,info,rcond);
    else
      HexaDobl_Newton_Matrix_Series.SVD_Newton_Step(p,deg,z,info,rcond);
    end if;
    put("rcond = "); put(rcond,3); new_line;
    if one + rcond /= one then
      Complex_Series_and_Polynomials.Filter(z,tol);
      put_line("The updated power series solution :");
      Complex_Series_and_Polynomials_io.put(z);
      eva := HexaDobl_CSeries_Poly_SysFun.Eval(p,z);
      Complex_Series_and_Polynomials.Filter(eva,tol);
      put_line("The evaluated solution :");
      Complex_Series_and_Polynomials_io.put(eva);
    end if;
  end Test_SVD_Newton_Step;

  procedure Test_SVD_Newton_Steps
              ( p : in HexaDobl_CSeries_Poly_Systems.Poly_Sys;
                degree,maxdeg : in integer32; nbrit : in integer32;
                x : in HexaDobl_Complex_Series_Vectors.Vector ) is

  -- DESCRIPTION :
  --  Does as many steps with Newton's method on the system p,
  --  as the value of nbrit, calculating with series x of the given degree,
  --  in double double precision.

    ans : character;
    otp : boolean;
    info : integer32;
    tol : constant double_float := 1.0E-12;
    eva : HexaDobl_Complex_Series_Vectors.Vector(p'range);
    deg : integer32 := degree;
    z : HexaDobl_Complex_Series_Vectors.Vector(x'range) := x;
    rcond : hexa_double;
    one : constant hexa_double := create(1.0);

  begin
    new_line;
    put("Do you want intermediate output ? (y/n) ");
    Ask_Yes_or_No(ans); otp := (ans = 'y');
    new_line;
    Complex_Series_and_Polynomials.Set_Degree(z,degree);
    if otp then
      HexaDobl_Newton_Matrix_Series.SVD_Newton_Steps
        (standard_output,p,deg,maxdeg,nbrit,z,info,rcond);
    else
      HexaDobl_Newton_Matrix_Series.SVD_Newton_Steps
        (p,deg,maxdeg,nbrit,z,info,rcond);
    end if;
    put("rcond = "); put(rcond,3); new_line;
    if one + rcond /= one then
      Complex_Series_and_Polynomials.Filter(z,tol);
      put_line("The updated power series solution :");
      Complex_Series_and_Polynomials_io.put(z);
      eva := HexaDobl_CSeries_Poly_SysFun.Eval(p,z);
      Complex_Series_and_Polynomials.Filter(eva,tol);
      put_line("The evaluated solution :");
      Complex_Series_and_Polynomials_io.put(eva);
    end if;
  end Test_SVD_Newton_Steps;

  procedure Test_Echelon_Newton_Step
              ( p : in HexaDobl_CSeries_Poly_Systems.Poly_Sys;
                degree : in integer32;
                x : in HexaDobl_Complex_Series_Vectors.Vector ) is

  -- DESCRIPTION :
  --  Does one step with the Echelon Newton's method on the system p,
  --  calculating with series x of the given degree,
  --  in double double precision.

    det : HexaDobl_Complex_Numbers.Complex_Number;
    tol : constant double_float := 1.0E-12;
    eva : HexaDobl_Complex_Series_Vectors.Vector(p'range);
    ans : character;
    otp : boolean;
    z : HexaDobl_Complex_Series_Vectors.Vector(x'range);

  begin
    new_line;
    put("Do you want intermediate output ? (y/n) ");
    Ask_Yes_or_No(ans); otp := (ans = 'y');
    new_line;
    z := x;
    Complex_Series_and_Polynomials.Set_Degree(z,degree);
    if otp then
      HexaDobl_Newton_Matrix_Series.Echelon_Newton_Step
        (standard_output,p,degree,z,det);
    else
      HexaDobl_Newton_Matrix_Series.Echelon_Newton_Step(p,degree,z,det);
    end if;
    put("det : "); put(det); new_line;
    Complex_Series_and_Polynomials.Filter(z,tol);
    put_line("The updated power series solution :");
    Complex_Series_and_Polynomials_io.put(z);
    eva := HexaDobl_CSeries_Poly_SysFun.Eval(p,z);
    Complex_Series_and_Polynomials.Filter(eva,tol);
    put_line("The evaluated solution :");
    Complex_Series_and_Polynomials_io.put(eva);
  end Test_Echelon_Newton_Step;

  procedure Echelon_Newton_on_Matrix_Series
              ( p : in HexaDobl_CSeries_Poly_Systems.Poly_Sys;
                otp : in boolean;
                degree,maxdeg : in integer32; nbrit : in integer32;
                x : in HexaDobl_Complex_Series_Vectors.Vector ) is

  -- DESCRIPTION :
  --   Runs the Echelon Newton's method on matrix series,
  --   on the system p, starting at x, of the given degree,
  --   with as many iterations as the value of nbrit,
  --   in double double precision.
  --   If otp, then extra output is written to screen.

    det : HexaDobl_Complex_Numbers.Complex_Number;
    tol : constant double_float := 1.0E-12;
    eva : HexaDobl_Complex_Series_Vectors.Vector(p'range);
    deg : integer32 := degree;
    z : HexaDobl_Complex_Series_Vectors.Vector(x'range) := x;

  begin
    Complex_Series_and_Polynomials.Set_Degree(z,deg);
    if otp then
      HexaDobl_Newton_Matrix_Series.Echelon_Newton_Steps
        (standard_output,p,deg,maxdeg,nbrit,z,det);
    else
      HexaDobl_Newton_Matrix_Series.Echelon_Newton_Steps
        (p,deg,maxdeg,nbrit,z,det);
    end if;
    put("det : "); put(det); new_line;
    Complex_Series_and_Polynomials.Filter(z,tol);
    put_line("The updated power series solution :");
    Complex_Series_and_Polynomials_io.put(z);
    eva := HexaDobl_CSeries_Poly_SysFun.Eval(p,z);
    Complex_Series_and_Polynomials.Filter(eva,tol);
    put_line("The evaluated solution :");
    Complex_Series_and_Polynomials_io.put(eva);
  end Echelon_Newton_on_Matrix_Series;

  procedure Test_Echelon_Newton_Steps
              ( p : in HexaDobl_CSeries_Poly_Systems.Poly_Sys;
                degree,maxdeg : in integer32; nbrit : in integer32;
                x : in HexaDobl_Complex_Series_Vectors.Vector ) is

  -- DESCRIPTION :
  --  Does as many steps with Newton's method on the system p,
  --  as the value of nbrit, calculating with series x of the given degree,
  --  in double double precision.

    ans : character;
    otp : boolean;

  begin
    new_line;
    put("Do you want intermediate output ? (y/n) ");
    Ask_Yes_or_No(ans); otp := (ans = 'y');
    Echelon_Newton_on_Matrix_Series(p,otp,degree,maxdeg,nbrit,x);
  end Test_Echelon_Newton_Steps;

  procedure Test_LU_Newton is

    ls : HexaDobl_CSeries_Poly_Systems.Link_to_Poly_Sys;
    sol : HexaDobl_Complex_Series_Vectors.Link_to_Vector;
    idx,degree,maxdeg,nbr,dim : integer32 := 0;

    use HexaDobl_CSeries_Polynomials;

  begin
    new_line;
    put("Give the index of the series variable : "); get(idx);
    new_line;
    Complex_Series_and_Polynomials_io.get(ls,idx);
    new_line;
    put_line("The polynomial system : ");
    dim := integer32(Number_of_Unknowns(ls(ls'first)));
    Complex_Series_and_Polynomials_io.put(ls.all,idx);
    Read_Series_Vector(sol,dim,idx);
    new_line;
    put("Give the start degree of the computations : "); get(degree);
    new_line;
    put("Give the number of Newton steps : "); get(nbr);
    if nbr = 1 then
      Test_LU_Newton_Step(ls.all,degree,sol.all);
    else
      put("Give the maximal degree of the series : "); get(maxdeg);
      Test_LU_Newton_Steps(ls.all,degree,maxdeg,nbr,sol.all);
    end if;
  end Test_LU_Newton;

  procedure Test_QR_Newton is

    ls : HexaDobl_CSeries_Poly_Systems.Link_to_Poly_Sys;
    sol : HexaDobl_Complex_Series_Vectors.Link_to_Vector;
    idx,degree,maxdeg,nbr,dim : integer32 := 0;

    use HexaDobl_CSeries_Polynomials;

  begin
    new_line;
    put("Give the index of the series variable : "); get(idx);
    new_line;
    Complex_Series_and_Polynomials_io.get(ls,idx);
    new_line;
    dim := integer32(Number_of_Unknowns(ls(ls'first)));
    put("The number of variables : "); put(dim,1); new_line;
    put_line("The polynomial system : ");
    Complex_Series_and_Polynomials_io.put(ls.all,idx);
    Read_Series_Vector(sol,dim,idx);
    new_line;
    put("Give the start degree of the computations : "); get(degree);
    new_line;
    put("Give the number of Newton steps : "); get(nbr);
    if nbr = 1 then
      Test_QR_Newton_Step(ls.all,degree,sol.all);
    else
      put("Give the maximal degree of the series : "); get(maxdeg);
      Test_QR_Newton_Steps(ls.all,degree,maxdeg,nbr,sol.all);
    end if;
  end Test_QR_Newton;

  procedure Test_SVD_Newton is

    ls : HexaDobl_CSeries_Poly_Systems.Link_to_Poly_Sys;
    sol : HexaDobl_Complex_Series_Vectors.Link_to_Vector;
    idx,degree,maxdeg,nbr,dim : integer32 := 0;

    use HexaDobl_CSeries_Polynomials;

  begin
    new_line;
    put("Give the index of the series variable : "); get(idx);
    new_line;
    Complex_Series_and_Polynomials_io.get(ls,idx);
    new_line;
    dim := integer32(Number_of_Unknowns(ls(ls'first)));
    put("The number of variables : "); put(dim,1); new_line;
    put_line("The polynomial system : ");
    Complex_Series_and_Polynomials_io.put(ls.all,idx);
    Read_Series_Vector(sol,dim,idx);
    new_line;
    put("Give the start degree of the computations : "); get(degree);
    new_line;
    put("Give the number of Newton steps : "); get(nbr);
    if nbr = 1 then
      Test_SVD_Newton_Step(ls.all,degree,sol.all);
    else
      put("Give the maximal degree of the series : "); get(maxdeg);
      Test_SVD_Newton_Steps(ls.all,degree,maxdeg,nbr,sol.all);
    end if;
  end Test_SVD_Newton;

  procedure Test_Echelon_Newton is

    ls : HexaDobl_CSeries_Poly_Systems.Link_to_Poly_Sys;
    sol : HexaDobl_Complex_Series_Vectors.Link_to_Vector;
    idx,degree,maxdeg,nbr,dim : integer32 := 0;

    use HexaDobl_CSeries_Polynomials;

  begin
    new_line;
    put("Give the index of the series variable : "); get(idx);
    new_line;
    Complex_Series_and_Polynomials_io.get(ls,idx);
    new_line;
    dim := integer32(Number_of_Unknowns(ls(ls'first)));
    put("The number of variables : "); put(dim,1); new_line;
    put_line("The polynomial system : ");
    Complex_Series_and_Polynomials_io.put(ls.all,idx);
    Read_Series_Vector(sol,dim,idx);
    new_line;
    put("Give the start degree of the computations : "); get(degree);
    new_line;
    put("Give the number of Newton steps : "); get(nbr);
    if nbr = 1 then
      Test_Echelon_Newton_Step(ls.all,degree,sol.all);
    else
      put("Give the maximal degree of the series : "); get(maxdeg);
      Test_Echelon_Newton_Steps(ls.all,degree,maxdeg,nbr,sol.all);
    end if;
  end Test_Echelon_Newton;

  procedure Main is

    ans : character;

  begin
    new_line;
    put_line("MENU to test Newton's method on truncated power series :");
    put_line("  1. test LU factorization for square systems;");
    put_line("  2. test QR decomposition for overdetermined systems;");
    put_line("  3. test SVD decomposition om matrix series;");
    put_line("  4. test lower triangular echelon form for any system.");
    put("Type 1, 2, 3, or 4 to select the test : ");
    Ask_Alternative(ans,"1234");
    case ans is
      when '1' => Test_LU_Newton;
      when '2' => Test_QR_Newton;
      when '3' => Test_SVD_Newton;
      when '4' => Test_Echelon_Newton;
      when others => null;
    end case;
  end Main;

end Test_HD_Newton_Matrix_Series;
