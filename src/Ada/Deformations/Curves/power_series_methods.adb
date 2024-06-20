with Communications_with_User;           use Communications_with_User;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;
with Double_Double_Numbers_io;           use Double_Double_Numbers_io;
with Triple_Double_Numbers_io;           use Triple_Double_Numbers_io;
with Quad_Double_Numbers_io;             use Quad_Double_Numbers_io;
with Penta_Double_Numbers_io;            use Penta_Double_Numbers_io;
with Octo_Double_Numbers_io;             use Octo_Double_Numbers_io;
with Deca_Double_Numbers_io;             use Deca_Double_Numbers_io;
with Standard_Complex_Numbers_io;        use Standard_Complex_Numbers_io;
with DoblDobl_Complex_Numbers_io;        use DoblDobl_Complex_Numbers_io;
with TripDobl_Complex_Numbers_io;        use TripDobl_Complex_Numbers_io;
with QuadDobl_Complex_Numbers_io;        use QuadDobl_Complex_Numbers_io;
with Standard_Complex_Vectors;
with Standard_Complex_Solutions_io;
with DoblDobl_Complex_Vectors;
with DoblDobl_Complex_Solutions_io;
with TripDobl_Complex_Vectors;
with TripDobl_Complex_Solutions_io;
with QuadDobl_Complex_Vectors;
with QuadDobl_Complex_Solutions_io;
with PentDobl_Complex_Vectors;
with PentDobl_Complex_Solutions_io;
with OctoDobl_Complex_Vectors;
with OctoDobl_Complex_Solutions_io;
with DecaDobl_Complex_Vectors;
with DecaDobl_Complex_Solutions_io;
with Complex_Series_and_Polynomials;
with Complex_Series_and_Polynomials_io;
with Standard_Newton_Matrix_Series;
with DoblDobl_Newton_Matrix_Series;
with TripDobl_Newton_Matrix_Series;
with QuadDobl_Newton_Matrix_Series;
with PentDobl_Newton_Matrix_Series;
with OctoDobl_Newton_Matrix_Series;
with DecaDobl_Newton_Matrix_Series;

package body Power_Series_Methods is

  procedure Run_LU_Newton
             ( maxdeg,nbrit : in integer32;
               p : in Standard_CSeries_Poly_Systems.Poly_Sys;
               s : in out Standard_Complex_Series_Vectors.Vector;
               info : out integer32; verbose : in boolean := false;
               vrblvl : in integer32 := 0 ) is
  begin
    if vrblvl > 0
     then put_line("-> in power_series_methods.Run_LU_Newton 1 ...");
    end if;
    Run_LU_Newton(standard_output,maxdeg,nbrit,p,s,info,verbose,vrblvl-1);
  end Run_LU_Newton;

  procedure Run_LU_Newton
             ( file : in file_type; maxdeg,nbrit : in integer32;
               p : in Standard_CSeries_Poly_Systems.Poly_Sys;
               s : in out Standard_Complex_Series_Vectors.Vector;
               info : out integer32; verbose : in boolean := false;
               vrblvl : in integer32 := 0 ) is

    use Standard_Newton_Matrix_Series; -- use Standard_Newton_Series;

    order : integer32 := 1;
    tol : constant double_float := 1.0E-12;
    eva : Standard_Complex_Series_Vectors.Vector(p'range);

  begin
    if vrblvl > 0
     then put_line("-> in power_series_methods.Run_LU_Newton 2 ...");
    end if;
    if not verbose then
      LU_Newton_Steps(p,order,maxdeg,nbrit,s,info,vrblvl-1);
    else
      LU_Newton_Steps(file,p,order,maxdeg,nbrit,s,info,vrblvl-1);
      if info /= 0 then
        put(file,"info = "); put(file,info,1); new_line(file);
      else
        Complex_Series_and_Polynomials.Filter(s,tol);
        put_line(file,"The updated power series solution :");
        Complex_Series_and_Polynomials_io.put(file,s);
        eva := Standard_CSeries_Poly_SysFun.Eval(p,s);
        Complex_Series_and_Polynomials.Filter(eva,tol);
        put_line(file,"The evaluated solution :");
        Complex_Series_and_Polynomials_io.put(file,eva);
        Standard_Complex_Series_Vectors.Clear(eva);
      end if;
    end if;
  end Run_LU_Newton;

  procedure Run_LU_Newton
             ( maxdeg,nbrit : in integer32;
               p : in DoblDobl_CSeries_Poly_Systems.Poly_Sys;
               s : in out DoblDobl_Complex_Series_Vectors.Vector;
               info : out integer32; verbose : in boolean := false;
               vrblvl : in integer32 := 0 ) is
  begin
    if vrblvl > 0
     then put_line("-> in power_series_methods.Run_LU_Newton 3 ...");
    end if;
    Run_LU_Newton(standard_output,maxdeg,nbrit,p,s,info,verbose,vrblvl-1);
  end Run_LU_Newton;

  procedure Run_LU_Newton
             ( file : in file_type; maxdeg,nbrit : in integer32;
               p : in DoblDobl_CSeries_Poly_Systems.Poly_Sys;
               s : in out DoblDobl_Complex_Series_Vectors.Vector;
               info : out integer32; verbose : in boolean := false;
               vrblvl : in integer32 := 0 ) is

    use DoblDobl_Newton_Matrix_Series; -- use DoblDobl_Newton_Series;

    order : integer32 := 1;
    tol : constant double_float := 1.0E-24;
    eva : DoblDobl_Complex_Series_Vectors.Vector(p'range);

  begin
    if vrblvl > 0
     then put_line("-> in power_series_methods.Run_LU_Newton 4 ...");
    end if;
    if not verbose then
      LU_Newton_Steps(p,order,maxdeg,nbrit,s,info,vrblvl-1);
    else
      LU_Newton_Steps(file,p,order,maxdeg,nbrit,s,info,vrblvl-1);
      if info /= 0 then
        put(file,"info = "); put(file,info,1); new_line(file);
      else
        Complex_Series_and_Polynomials.Filter(s,tol);
        put_line(file,"The updated power series solution :");
        Complex_Series_and_Polynomials_io.put(file,s);
        eva := DoblDobl_CSeries_Poly_SysFun.Eval(p,s);
        Complex_Series_and_Polynomials.Filter(eva,tol);
        put_line(file,"The evaluated solution :");
        Complex_Series_and_Polynomials_io.put(file,eva);
        DoblDobl_Complex_Series_Vectors.Clear(eva);
      end if;
    end if;
  end Run_LU_Newton;

  procedure Run_LU_Newton
             ( maxdeg,nbrit : in integer32;
               p : in TripDobl_CSeries_Poly_Systems.Poly_Sys;
               s : in out TripDobl_Complex_Series_Vectors.Vector;
               info : out integer32; verbose : in boolean := false;
               vrblvl : in integer32 := 0 ) is
  begin
    if vrblvl > 0
     then put_line("-> in power_series_methods.Run_LU_Newton 5 ...");
    end if;
    Run_LU_Newton(standard_output,maxdeg,nbrit,p,s,info,verbose,vrblvl-1);
  end Run_LU_Newton;

  procedure Run_LU_Newton
             ( file : in file_type; maxdeg,nbrit : in integer32;
               p : in TripDobl_CSeries_Poly_Systems.Poly_Sys;
               s : in out TripDobl_Complex_Series_Vectors.Vector;
               info : out integer32; verbose : in boolean := false;
               vrblvl : in integer32 := 0 ) is

    use TripDobl_Newton_Matrix_Series; -- use TripDobl_Newton_Series;

    order : integer32 := 1;
    tol : constant double_float := 1.0E-24;
    eva : TripDobl_Complex_Series_Vectors.Vector(p'range);

  begin
    if vrblvl > 0
     then put_line("-> in power_series_methods.Run_LU_Newton 6 ...");
    end if;
    if not verbose then
      LU_Newton_Steps(p,order,maxdeg,nbrit,s,info,vrblvl-1);
    else
      LU_Newton_Steps(file,p,order,maxdeg,nbrit,s,info,vrblvl-1);
      if info /= 0 then
        put(file,"info = "); put(file,info,1); new_line(file);
      else
        Complex_Series_and_Polynomials.Filter(s,tol);
        put_line(file,"The updated power series solution :");
        Complex_Series_and_Polynomials_io.put(file,s);
        eva := TripDobl_CSeries_Poly_SysFun.Eval(p,s);
        Complex_Series_and_Polynomials.Filter(eva,tol);
        put_line(file,"The evaluated solution :");
        Complex_Series_and_Polynomials_io.put(file,eva);
        TripDobl_Complex_Series_Vectors.Clear(eva);
      end if;
    end if;
  end Run_LU_Newton;

  procedure Run_LU_Newton
             ( maxdeg,nbrit : in integer32;
               p : in QuadDobl_CSeries_Poly_Systems.Poly_Sys;
               s : in out QuadDobl_Complex_Series_Vectors.Vector;
               info : out integer32; verbose : in boolean := false;
               vrblvl : in integer32 := 0 ) is
  begin
    if vrblvl > 0
     then put_line("-> in power_series_methods.Run_LU_Newton 7 ...");
    end if;
    Run_LU_Newton(standard_output,maxdeg,nbrit,p,s,info,verbose,vrblvl-1);
  end Run_LU_Newton;

  procedure Run_LU_Newton
             ( file : in file_type; maxdeg,nbrit : in integer32;
               p : in QuadDobl_CSeries_Poly_Systems.Poly_Sys;
               s : in out QuadDobl_Complex_Series_Vectors.Vector;
               info : out integer32; verbose : in boolean := false;
               vrblvl : in integer32 := 0 ) is

    use QuadDobl_Newton_Matrix_Series; -- use QuadDobl_Newton_Series;

    order : integer32 := 1;
    tol : constant double_float := 1.0E-32;
    eva : QuadDobl_Complex_Series_Vectors.Vector(p'range);

  begin
    if vrblvl > 0
     then put_line("-> in power_series_methods.Run_LU_Newton 8 ...");
    end if;
    if not verbose then
      LU_Newton_Steps(p,order,maxdeg,nbrit,s,info,vrblvl-1);
    else
      LU_Newton_Steps(file,p,order,maxdeg,nbrit,s,info,vrblvl-1);
      if info /= 0 then
        put(file,"info = "); put(file,info,1); new_line(file);
      else
        Complex_Series_and_Polynomials.Filter(s,tol);
        put_line(file,"The updated power series solution :");
        Complex_Series_and_Polynomials_io.put(file,s);
        eva := QuadDobl_CSeries_Poly_SysFun.Eval(p,s);
        Complex_Series_and_Polynomials.Filter(eva,tol);
        put_line(file,"The evaluated solution :");
        Complex_Series_and_Polynomials_io.put(file,eva);
        QuadDobl_Complex_Series_Vectors.Clear(eva);
      end if;
    end if;
  end Run_LU_Newton;

  procedure Run_LU_Newton
             ( maxdeg,nbrit : in integer32;
               p : in Standard_CSeries_Poly_Systems.Poly_Sys;
               s : in out Standard_Complex_Series_Vectors.Vector;
               rcond : out double_float; verbose : in boolean := false;
               vrblvl : in integer32 := 0 ) is
  begin
    if vrblvl > 0
     then put_line("-> in power_series_methods.Run_LU_Newton 9 ...");
    end if;
    Run_LU_Newton(standard_output,maxdeg,nbrit,p,s,rcond,verbose,vrblvl-1);
  end Run_LU_Newton;

  procedure Run_LU_Newton
             ( file : in file_type; maxdeg,nbrit : in integer32;
               p : in Standard_CSeries_Poly_Systems.Poly_Sys;
               s : in out Standard_Complex_Series_Vectors.Vector;
               rcond : out double_float; verbose : in boolean := false;
               vrblvl : in integer32 := 0 ) is

    use Standard_Newton_Matrix_Series; -- use Standard_Newton_Series;

    order : integer32 := 1;
    tol : constant double_float := 1.0E-12;
    eva : Standard_Complex_Series_Vectors.Vector(p'range);

  begin
    if vrblvl > 0
     then put_line("-> in power_series_methods.Run_LU_Newton 10 ...");
    end if;
    if not verbose then
      LU_Newton_Steps(p,order,maxdeg,nbrit,s,rcond,vrblvl-1);
    else
      LU_Newton_Steps(file,p,order,maxdeg,nbrit,s,rcond,vrblvl-1);
      put(file,"rcond : "); put(file,rcond,3); new_line(file);
      if 1.0 + rcond /= 1.0 then
        Complex_Series_and_Polynomials.Filter(s,tol);
        put_line(file,"The updated power series solution :");
        Complex_Series_and_Polynomials_io.put(file,s);
        eva := Standard_CSeries_Poly_SysFun.Eval(p,s);
        Complex_Series_and_Polynomials.Filter(eva,tol);
        put_line(file,"The evaluated solution :");
        Complex_Series_and_Polynomials_io.put(file,eva);
        Standard_Complex_Series_Vectors.Clear(eva);
      end if;
    end if;
  end Run_LU_Newton;

  procedure Run_LU_Newton
             ( maxdeg,nbrit : in integer32;
               p : in DoblDobl_CSeries_Poly_Systems.Poly_Sys;
               s : in out DoblDobl_Complex_Series_Vectors.Vector;
               rcond : out double_double; verbose : in boolean := false;
               vrblvl : in integer32 := 0 ) is
  begin
    if vrblvl > 0
     then put_line("-> in power_series_methods.Run_LU_Newton 11 ...");
    end if;
    Run_LU_Newton(standard_output,maxdeg,nbrit,p,s,rcond,verbose,vrblvl-1);
  end Run_LU_Newton;

  procedure Run_LU_Newton
             ( file : in file_type; maxdeg,nbrit : in integer32;
               p : in DoblDobl_CSeries_Poly_Systems.Poly_Sys;
               s : in out DoblDobl_Complex_Series_Vectors.Vector;
               rcond : out double_double; verbose : in boolean := false;
               vrblvl : in integer32 := 0 ) is

    use DoblDobl_Newton_Matrix_Series; -- use DoblDobl_Newton_Series;

    order : integer32 := 1;
    tol : constant double_float := 1.0E-24;
    eva : DoblDobl_Complex_Series_Vectors.Vector(p'range);
    one : constant double_double := create(1.0);

  begin
    if vrblvl > 0
     then put_line("-> in power_series_methods.Run_LU_Newton 12 ...");
    end if;
    if not verbose then
      LU_Newton_Steps(p,order,maxdeg,nbrit,s,rcond,vrblvl-1);
    else
      LU_Newton_Steps(file,p,order,maxdeg,nbrit,s,rcond,vrblvl-1);
      put(file,"rcond : "); put(file,rcond,3); new_line(file);
      if one + rcond /= one then
        Complex_Series_and_Polynomials.Filter(s,tol);
        put_line(file,"The updated power series solution :");
        Complex_Series_and_Polynomials_io.put(file,s);
        eva := DoblDobl_CSeries_Poly_SysFun.Eval(p,s);
        Complex_Series_and_Polynomials.Filter(eva,tol);
        put_line(file,"The evaluated solution :");
        Complex_Series_and_Polynomials_io.put(file,eva);
        DoblDobl_Complex_Series_Vectors.Clear(eva);
      end if;
    end if;
  end Run_LU_Newton;

  procedure Run_LU_Newton
             ( maxdeg,nbrit : in integer32;
               p : in TripDobl_CSeries_Poly_Systems.Poly_Sys;
               s : in out TripDobl_Complex_Series_Vectors.Vector;
               rcond : out triple_double; verbose : in boolean := false;
               vrblvl : in integer32 := 0 ) is
  begin
    if vrblvl > 0
     then put_line("-> in power_series_methods.Run_LU_Newton 13 ...");
    end if;
    Run_LU_Newton(standard_output,maxdeg,nbrit,p,s,rcond,verbose,vrblvl-1);
  end Run_LU_Newton;

  procedure Run_LU_Newton
             ( file : in file_type; maxdeg,nbrit : in integer32;
               p : in TripDobl_CSeries_Poly_Systems.Poly_Sys;
               s : in out TripDobl_Complex_Series_Vectors.Vector;
               rcond : out triple_double; verbose : in boolean := false;
               vrblvl : in integer32 := 0 ) is

    use TripDobl_Newton_Matrix_Series; -- use TripDobl_Newton_Series;

    order : integer32 := 1;
    tol : constant double_float := 1.0E-24;
    eva : TripDobl_Complex_Series_Vectors.Vector(p'range);
    one : constant triple_double := create(1.0);

  begin
    if vrblvl > 0
     then put_line("-> in power_series_methods.Run_LU_Newton 14 ...");
    end if;
    if not verbose then
      LU_Newton_Steps(p,order,maxdeg,nbrit,s,rcond,vrblvl-1);
    else
      LU_Newton_Steps(file,p,order,maxdeg,nbrit,s,rcond,vrblvl-1);
      put(file,"rcond : "); put(file,rcond,3); new_line(file);
      if one + rcond /= one then
        Complex_Series_and_Polynomials.Filter(s,tol);
        put_line(file,"The updated power series solution :");
        Complex_Series_and_Polynomials_io.put(file,s);
        eva := TripDobl_CSeries_Poly_SysFun.Eval(p,s);
        Complex_Series_and_Polynomials.Filter(eva,tol);
        put_line(file,"The evaluated solution :");
        Complex_Series_and_Polynomials_io.put(file,eva);
        TripDobl_Complex_Series_Vectors.Clear(eva);
      end if;
    end if;
  end Run_LU_Newton;

  procedure Run_LU_Newton
             ( maxdeg,nbrit : in integer32;
               p : in QuadDobl_CSeries_Poly_Systems.Poly_Sys;
               s : in out QuadDobl_Complex_Series_Vectors.Vector;
               rcond : out quad_double; verbose : in boolean := false;
               vrblvl : in integer32 := 0 ) is
  begin
    if vrblvl > 0
     then put_line("-> in power_series_methods.Run_LU_Newton 15 ...");
    end if;
    Run_LU_Newton(standard_output,maxdeg,nbrit,p,s,rcond,verbose,vrblvl-1);
  end Run_LU_Newton;

  procedure Run_LU_Newton
             ( file : in file_type; maxdeg,nbrit : in integer32;
               p : in QuadDobl_CSeries_Poly_Systems.Poly_Sys;
               s : in out QuadDobl_Complex_Series_Vectors.Vector;
               rcond : out quad_double; verbose : in boolean := false;
               vrblvl : in integer32 := 0 ) is

    use QuadDobl_Newton_Matrix_Series; -- use QuadDobl_Newton_Series;

    order : integer32 := 1;
    tol : constant double_float := 1.0E-32;
    eva : QuadDobl_Complex_Series_Vectors.Vector(p'range);
    one : constant quad_double := create(1.0);

  begin
    if vrblvl > 0
     then put_line("-> in power_series_methods.Run_LU_Newton 16 ...");
    end if;
    if not verbose then
      LU_Newton_Steps(p,order,maxdeg,nbrit,s,rcond,vrblvl-1);
    else
      LU_Newton_Steps(file,p,order,maxdeg,nbrit,s,rcond,vrblvl-1);
      put(file,"rcond : "); put(file,rcond,3); new_line(file);
      if one + rcond /= one then
        Complex_Series_and_Polynomials.Filter(s,tol);
        put_line(file,"The updated power series solution :");
        Complex_Series_and_Polynomials_io.put(file,s);
        eva := QuadDobl_CSeries_Poly_SysFun.Eval(p,s);
        Complex_Series_and_Polynomials.Filter(eva,tol);
        put_line(file,"The evaluated solution :");
        Complex_Series_and_Polynomials_io.put(file,eva);
        QuadDobl_Complex_Series_Vectors.Clear(eva);
      end if;
    end if;
  end Run_LU_Newton;

-- LU NEWTON ON COEFFICIENT-PARAMETER HOMOTOPIES :

  procedure Run_LU_Newton
             ( maxdeg,nbrit : in integer32;
               f : in Standard_CSeries_Poly_SysFun.Eval_Coeff_Poly_Sys;
               c : in Standard_Complex_Series_VecVecs.VecVec;
               ejm : in Standard_CSeries_Jaco_Matrices.Eval_Coeff_Jaco_Mat;
               mlt : in Standard_CSeries_Jaco_Matrices.Mult_Factors;
               s : in out Standard_Complex_Series_Vectors.Vector;
               info : out integer32; verbose : in boolean := false;
               vrblvl : in integer32 := 0 ) is
  begin
    if vrblvl > 0
     then put_line("-> in power_series_methods.Run_LU_Newton 17 ...");
    end if;
    Run_LU_Newton
      (standard_output,maxdeg,nbrit,f,c,ejm,mlt,s,info,verbose,vrblvl-1);
  end Run_LU_Newton;

  procedure Run_LU_Newton
             ( file : in file_type; maxdeg,nbrit : in integer32;
               f : in Standard_CSeries_Poly_SysFun.Eval_Coeff_Poly_Sys;
               c : in Standard_Complex_Series_VecVecs.VecVec;
               ejm : in Standard_CSeries_Jaco_Matrices.Eval_Coeff_Jaco_Mat;
               mlt : in Standard_CSeries_Jaco_Matrices.Mult_Factors;
               s : in out Standard_Complex_Series_Vectors.Vector;
               info : out integer32; verbose : in boolean := false;
               vrblvl : in integer32 := 0 ) is

    use Standard_Newton_Matrix_Series; -- use Standard_Newton_Series;

    order : integer32 := 1;
    tol : constant double_float := 1.0E-12;
    eva : Standard_Complex_Series_Vectors.Vector(f'range);

  begin
    if vrblvl > 0
     then put_line("-> in power_series_methods.Run_LU_Newton 18 ...");
    end if;
    if not verbose then
      LU_Newton_Steps(f,c,ejm,mlt,order,maxdeg,nbrit,s,info,vrblvl-1);
    else
      LU_Newton_Steps(file,f,c,ejm,mlt,order,maxdeg,nbrit,s,info,vrblvl-1);
      if info /= 0 then
        put(file,"info = "); put(file,info,1); new_line(file);
      else
        Complex_Series_and_Polynomials.Filter(s,tol);
        put_line(file,"The updated power series solution :");
        Complex_Series_and_Polynomials_io.put(file,s);
        eva := Standard_CSeries_Poly_SysFun.Eval(f,c,s);
        Complex_Series_and_Polynomials.Filter(eva,tol);
        put_line(file,"The evaluated solution :");
        Complex_Series_and_Polynomials_io.put(file,eva);
        Standard_Complex_Series_Vectors.Clear(eva);
      end if;
    end if;
  end Run_LU_Newton;

  procedure Run_LU_Newton
             ( maxdeg,nbrit : in integer32;
               f : in DoblDobl_CSeries_Poly_SysFun.Eval_Coeff_Poly_Sys;
               c : in DoblDobl_Complex_Series_VecVecs.VecVec;
               ejm : in DoblDobl_CSeries_Jaco_Matrices.Eval_Coeff_Jaco_Mat;
               mlt : in DoblDobl_CSeries_Jaco_Matrices.Mult_Factors;
               s : in out DoblDobl_Complex_Series_Vectors.Vector;
               info : out integer32; verbose : in boolean := false;
               vrblvl : in integer32 := 0 ) is
  begin
    if vrblvl > 0
     then put_line("-> in power_series_methods.Run_LU_Newton 19 ...");
    end if;
    Run_LU_Newton
      (standard_output,maxdeg,nbrit,f,c,ejm,mlt,s,info,verbose,vrblvl-1);
  end Run_LU_Newton;

  procedure Run_LU_Newton
             ( file : in file_type; maxdeg,nbrit : in integer32;
               f : in DoblDobl_CSeries_Poly_SysFun.Eval_Coeff_Poly_Sys;
               c : in DoblDobl_Complex_Series_VecVecs.VecVec;
               ejm : in DoblDobl_CSeries_Jaco_Matrices.Eval_Coeff_Jaco_Mat;
               mlt : in DoblDobl_CSeries_Jaco_Matrices.Mult_Factors;
               s : in out DoblDobl_Complex_Series_Vectors.Vector;
               info : out integer32; verbose : in boolean := false;
               vrblvl : in integer32 := 0 ) is

    use DoblDobl_Newton_Matrix_Series; -- use DoblDobl_Newton_Series;

    order : integer32 := 1;
    tol : constant double_float := 1.0E-12;
    eva : DoblDobl_Complex_Series_Vectors.Vector(f'range);

  begin
    if vrblvl > 0
     then put_line("-> in power_series_methods.Run_LU_Newton 20 ...");
    end if;
    if not verbose then
      LU_Newton_Steps(f,c,ejm,mlt,order,maxdeg,nbrit,s,info,vrblvl-1);
    else
      LU_Newton_Steps(file,f,c,ejm,mlt,order,maxdeg,nbrit,s,info,vrblvl-1);
      if info /= 0 then
        put(file,"info = "); put(file,info,1); new_line(file);
      else
        Complex_Series_and_Polynomials.Filter(s,tol);
        put_line(file,"The updated power series solution :");
        Complex_Series_and_Polynomials_io.put(file,s);
        eva := DoblDobl_CSeries_Poly_SysFun.Eval(f,c,s);
        Complex_Series_and_Polynomials.Filter(eva,tol);
        put_line(file,"The evaluated solution :");
        Complex_Series_and_Polynomials_io.put(file,eva);
        DoblDobl_Complex_Series_Vectors.Clear(eva);
      end if;
    end if;
  end Run_LU_Newton;

  procedure Run_LU_Newton
             ( maxdeg,nbrit : in integer32;
               f : in TripDobl_CSeries_Poly_SysFun.Eval_Coeff_Poly_Sys;
               c : in TripDobl_Complex_Series_VecVecs.VecVec;
               ejm : in TripDobl_CSeries_Jaco_Matrices.Eval_Coeff_Jaco_Mat;
               mlt : in TripDobl_CSeries_Jaco_Matrices.Mult_Factors;
               s : in out TripDobl_Complex_Series_Vectors.Vector;
               info : out integer32; verbose : in boolean := false;
               vrblvl : in integer32 := 0 ) is
  begin
    if vrblvl > 0
     then put_line("-> in power_series_methods.Run_LU_Newton 21 ...");
    end if;
    Run_LU_Newton
      (standard_output,maxdeg,nbrit,f,c,ejm,mlt,s,info,verbose,vrblvl-1);
  end Run_LU_Newton;

  procedure Run_LU_Newton
             ( file : in file_type; maxdeg,nbrit : in integer32;
               f : in TripDobl_CSeries_Poly_SysFun.Eval_Coeff_Poly_Sys;
               c : in TripDobl_Complex_Series_VecVecs.VecVec;
               ejm : in TripDobl_CSeries_Jaco_Matrices.Eval_Coeff_Jaco_Mat;
               mlt : in TripDobl_CSeries_Jaco_Matrices.Mult_Factors;
               s : in out TripDobl_Complex_Series_Vectors.Vector;
               info : out integer32; verbose : in boolean := false;
               vrblvl : in integer32 := 0 ) is

    use TripDobl_Newton_Matrix_Series;

    order : integer32 := 1;
    tol : constant double_float := 1.0E-12;
    eva : TripDobl_Complex_Series_Vectors.Vector(f'range);

  begin
    if vrblvl > 0
     then put_line("-> in power_series_methods.Run_LU_Newton 22 ...");
    end if;
    if not verbose then
      LU_Newton_Steps(f,c,ejm,mlt,order,maxdeg,nbrit,s,info,vrblvl-1);
    else
      LU_Newton_Steps(file,f,c,ejm,mlt,order,maxdeg,nbrit,s,info,vrblvl-1);
      if info /= 0 then
        put(file,"info = "); put(file,info,1); new_line(file);
      else
        Complex_Series_and_Polynomials.Filter(s,tol);
        put_line(file,"The updated power series solution :");
        Complex_Series_and_Polynomials_io.put(file,s);
        eva := TripDobl_CSeries_Poly_SysFun.Eval(f,c,s);
        Complex_Series_and_Polynomials.Filter(eva,tol);
        put_line(file,"The evaluated solution :");
        Complex_Series_and_Polynomials_io.put(file,eva);
        TripDobl_Complex_Series_Vectors.Clear(eva);
      end if;
    end if;
  end Run_LU_Newton;

  procedure Run_LU_Newton
             ( maxdeg,nbrit : in integer32;
               f : in QuadDobl_CSeries_Poly_SysFun.Eval_Coeff_Poly_Sys;
               c : in QuadDobl_Complex_Series_VecVecs.VecVec;
               ejm : in QuadDobl_CSeries_Jaco_Matrices.Eval_Coeff_Jaco_Mat;
               mlt : in QuadDobl_CSeries_Jaco_Matrices.Mult_Factors;
               s : in out QuadDobl_Complex_Series_Vectors.Vector;
               info : out integer32; verbose : in boolean := false;
               vrblvl : in integer32 := 0 ) is
  begin
    if vrblvl > 0
     then put_line("-> in power_series_methods.Run_LU_Newton 23 ...");
    end if;
    Run_LU_Newton
      (standard_output,maxdeg,nbrit,f,c,ejm,mlt,s,info,verbose,vrblvl-1);
  end Run_LU_Newton;

  procedure Run_LU_Newton
             ( file : in file_type; maxdeg,nbrit : in integer32;
               f : in QuadDobl_CSeries_Poly_SysFun.Eval_Coeff_Poly_Sys;
               c : in QuadDobl_Complex_Series_VecVecs.VecVec;
               ejm : in QuadDobl_CSeries_Jaco_Matrices.Eval_Coeff_Jaco_Mat;
               mlt : in QuadDobl_CSeries_Jaco_Matrices.Mult_Factors;
               s : in out QuadDobl_Complex_Series_Vectors.Vector;
               info : out integer32; verbose : in boolean := false;
               vrblvl : in integer32 := 0 ) is

    use QuadDobl_Newton_Matrix_Series; -- use QuadDobl_Newton_Series;

    order : integer32 := 1;
    tol : constant double_float := 1.0E-12;
    eva : QuadDobl_Complex_Series_Vectors.Vector(f'range);

  begin
    if vrblvl > 0
     then put_line("-> in power_series_methods.Run_LU_Newton 24 ...");
    end if;
    if not verbose then
      LU_Newton_Steps(f,c,ejm,mlt,order,maxdeg,nbrit,s,info,vrblvl-1);
    else
      LU_Newton_Steps(file,f,c,ejm,mlt,order,maxdeg,nbrit,s,info,vrblvl-1);
      if info /= 0 then
        put(file,"info = "); put(file,info,1); new_line(file);
      else
        Complex_Series_and_Polynomials.Filter(s,tol);
        put_line(file,"The updated power series solution :");
        Complex_Series_and_Polynomials_io.put(file,s);
        eva := QuadDobl_CSeries_Poly_SysFun.Eval(f,c,s);
        Complex_Series_and_Polynomials.Filter(eva,tol);
        put_line(file,"The evaluated solution :");
        Complex_Series_and_Polynomials_io.put(file,eva);
        QuadDobl_Complex_Series_Vectors.Clear(eva);
      end if;
    end if;
  end Run_LU_Newton;

  procedure Run_QR_Newton
             ( maxdeg,nbrit : in integer32;
               p : in Standard_CSeries_Poly_Systems.Poly_Sys;
               s : in out Standard_Complex_Series_Vectors.Vector;
               verbose : in boolean := false;
               vrblvl : in integer32 := 0 ) is
  begin
    if vrblvl > 0
     then put_line("-> in power_series_methods.Run_QR_Newton 1 ...");
    end if;
    Run_QR_Newton(standard_output,maxdeg,nbrit,p,s,verbose,vrblvl-1);
  end Run_QR_Newton;

  procedure Run_QR_Newton
             ( file : in file_type; maxdeg,nbrit : in integer32;
               p : in Standard_CSeries_Poly_Systems.Poly_Sys;
               s : in out Standard_Complex_Series_Vectors.Vector;
               verbose : in boolean := false;
               vrblvl : in integer32 := 0 ) is

    use Standard_Newton_Matrix_Series; -- use Standard_Newton_Series;

    order : integer32 := 1;
    info : integer32;
    tol : constant double_float := 1.0E-12;
    eva : Standard_Complex_Series_Vectors.Vector(p'range);

  begin
    if vrblvl > 0
     then put_line("-> in power_series_methods.Run_QR_Newton 2 ...");
    end if;
    if not verbose then
      QR_Newton_Steps(p,order,maxdeg,nbrit,s,info,vrblvl-1);
    else
      QR_Newton_Steps(file,p,order,maxdeg,nbrit,s,info,vrblvl-1);
      if info /= 0 then
        put(file,"info = "); put(file,info,1); new_line(file);
      else
        Complex_Series_and_Polynomials.Filter(s,tol);
        put_line(file,"The updated power series solution :");
        Complex_Series_and_Polynomials_io.put(file,s);
        eva := Standard_CSeries_Poly_SysFun.Eval(p,s);
        Complex_Series_and_Polynomials.Filter(eva,tol);
        put_line(file,"The evaluated solution :");
        Complex_Series_and_Polynomials_io.put(file,eva);
        Standard_Complex_Series_Vectors.Clear(eva);
      end if;
    end if;
  end Run_QR_Newton;

  procedure Run_QR_Newton
             ( maxdeg,nbrit : in integer32;
               p : in DoblDobl_CSeries_Poly_Systems.Poly_Sys;
               s : in out DoblDobl_Complex_Series_Vectors.Vector;
               verbose : in boolean := false;
               vrblvl : in integer32 := 0 ) is
  begin
    if vrblvl > 0
     then put_line("-> in power_series_methods.Run_QR_Newton 3 ...");
    end if;
    Run_QR_Newton(standard_output,maxdeg,nbrit,p,s,verbose,vrblvl-1);
  end Run_QR_Newton;

  procedure Run_QR_Newton
             ( file : in file_type; maxdeg,nbrit : in integer32;
               p : in DoblDobl_CSeries_Poly_Systems.Poly_Sys;
               s : in out DoblDobl_Complex_Series_Vectors.Vector;
               verbose : in boolean := false;
               vrblvl : in integer32 := 0 ) is

    use DoblDobl_Newton_Matrix_Series;

    order : integer32 := 1;
    info : integer32;
    tol : constant double_float := 1.0E-24;
    eva : DoblDobl_Complex_Series_Vectors.Vector(p'range);

  begin
    if vrblvl > 0
     then put_line("-> in power_series_methods.Run_QR_Newton 4 ...");
    end if;
    if not verbose then
      QR_Newton_Steps(p,order,maxdeg,nbrit,s,info,vrblvl-1);
    else
      QR_Newton_Steps(file,p,order,maxdeg,nbrit,s,info,vrblvl-1);
      if info /= 0 then
        put(file,"info = "); put(file,info,1); new_line(file);
      else
        Complex_Series_and_Polynomials.Filter(s,tol);
        put_line(file,"The updated power series solution :");
        Complex_Series_and_Polynomials_io.put(file,s);
        eva := DoblDobl_CSeries_Poly_SysFun.Eval(p,s);
        Complex_Series_and_Polynomials.Filter(eva,tol);
        put_line(file,"The evaluated solution :");
        Complex_Series_and_Polynomials_io.put(file,eva);
        DoblDobl_Complex_Series_Vectors.Clear(eva);
      end if;
    end if;
  end Run_QR_Newton;

  procedure Run_QR_Newton
             ( maxdeg,nbrit : in integer32;
               p : in TripDobl_CSeries_Poly_Systems.Poly_Sys;
               s : in out TripDobl_Complex_Series_Vectors.Vector;
               verbose : in boolean := false;
               vrblvl : in integer32 := 0 ) is
  begin
    if vrblvl > 0
     then put_line("-> in power_series_methods.Run_QR_Newton 5 ...");
    end if;
    Run_QR_Newton(standard_output,maxdeg,nbrit,p,s,verbose,vrblvl-1);
  end Run_QR_Newton;

  procedure Run_QR_Newton
             ( file : in file_type; maxdeg,nbrit : in integer32;
               p : in TripDobl_CSeries_Poly_Systems.Poly_Sys;
               s : in out TripDobl_Complex_Series_Vectors.Vector;
               verbose : in boolean := false;
               vrblvl : in integer32 := 0 ) is

    use TripDobl_Newton_Matrix_Series;

    order : integer32 := 1;
    info : integer32;
    tol : constant double_float := 1.0E-24;
    eva : TripDobl_Complex_Series_Vectors.Vector(p'range);

  begin
    if vrblvl > 0
     then put_line("-> in power_series_methods.Run_QR_Newton 6 ...");
    end if;
    if not verbose then
      QR_Newton_Steps(p,order,maxdeg,nbrit,s,info,vrblvl-1);
    else
      QR_Newton_Steps(file,p,order,maxdeg,nbrit,s,info,vrblvl-1);
      if info /= 0 then
        put(file,"info = "); put(file,info,1); new_line(file);
      else
        Complex_Series_and_Polynomials.Filter(s,tol);
        put_line(file,"The updated power series solution :");
        Complex_Series_and_Polynomials_io.put(file,s);
        eva := TripDobl_CSeries_Poly_SysFun.Eval(p,s);
        Complex_Series_and_Polynomials.Filter(eva,tol);
        put_line(file,"The evaluated solution :");
        Complex_Series_and_Polynomials_io.put(file,eva);
        TripDobl_Complex_Series_Vectors.Clear(eva);
      end if;
    end if;
  end Run_QR_Newton;

  procedure Run_QR_Newton
             ( maxdeg,nbrit : in integer32;
               p : in QuadDobl_CSeries_Poly_Systems.Poly_Sys;
               s : in out QuadDobl_Complex_Series_Vectors.Vector;
               verbose : in boolean := false;
               vrblvl : in integer32 := 0 ) is
  begin
    if vrblvl > 0
     then put_line("-> in power_series_methods.Run_QR_Newton 7 ...");
    end if;
    Run_QR_Newton(standard_output,maxdeg,nbrit,p,s,verbose,vrblvl-1);
  end Run_QR_Newton;

  procedure Run_QR_Newton
             ( file : in file_type; maxdeg,nbrit : in integer32;
               p : in QuadDobl_CSeries_Poly_Systems.Poly_Sys;
               s : in out QuadDobl_Complex_Series_Vectors.Vector;
               verbose : in boolean := false;
               vrblvl : in integer32 := 0 ) is

    use QuadDobl_Newton_Matrix_Series; -- use QuadDobl_Newton_Series;

    order : integer32 := 1;
    info : integer32;
    tol : constant double_float := 1.0E-32;
    eva : QuadDobl_Complex_Series_Vectors.Vector(p'range);

  begin
    if vrblvl > 0
     then put_line("-> in power_series_methods.Run_QR_Newton 8 ...");
    end if;
    if not verbose then
      QR_Newton_Steps(p,order,maxdeg,nbrit,s,info,vrblvl-1);
    else
      QR_Newton_Steps(file,p,order,maxdeg,nbrit,s,info,vrblvl-1);
      if info /= 0 then
        put(file,"info = "); put(file,info,1); new_line(file);
      else
        Complex_Series_and_Polynomials.Filter(s,tol);
        put_line(file,"The updated power series solution :");
        Complex_Series_and_Polynomials_io.put(file,s);
        eva := QuadDobl_CSeries_Poly_SysFun.Eval(p,s);
        Complex_Series_and_Polynomials.Filter(eva,tol);
        put_line(file,"The evaluated solution :");
        Complex_Series_and_Polynomials_io.put(file,eva);
        QuadDobl_Complex_Series_Vectors.Clear(eva);
      end if;
    end if;
  end Run_QR_Newton;

-- QR NEWTON ON COEFFICIENT-PARAMETER HOMOTOPIES :

  procedure Run_QR_Newton
              ( maxdeg,nbrit : in integer32;
                f : in Standard_CSeries_Poly_SysFun.Eval_Coeff_Poly_Sys;
                c : in Standard_Complex_Series_VecVecs.VecVec;
                ejm : in Standard_CSeries_Jaco_Matrices.Eval_Coeff_Jaco_Mat;
                mlt : in Standard_CSeries_Jaco_Matrices.Mult_Factors;
                s : in out Standard_Complex_Series_Vectors.Vector;
                verbose : in boolean := false;
                vrblvl : in integer32 := 0 ) is
  begin
    if vrblvl > 0
     then put_line("-> in power_series_methods.Run_QR_Newton 9 ...");
    end if;
    Run_QR_Newton(standard_output,maxdeg,nbrit,f,c,ejm,mlt,s,verbose,vrblvl-1);
  end Run_QR_Newton;

  procedure Run_QR_Newton
              ( file : in file_type; maxdeg,nbrit : in integer32;
                f : in Standard_CSeries_Poly_SysFun.Eval_Coeff_Poly_Sys;
                c : in Standard_Complex_Series_VecVecs.VecVec;
                ejm : in Standard_CSeries_Jaco_Matrices.Eval_Coeff_Jaco_Mat;
                mlt : in Standard_CSeries_Jaco_Matrices.Mult_Factors;
                s : in out Standard_Complex_Series_Vectors.Vector;
                verbose : in boolean := false;
                vrblvl : in integer32 := 0 ) is

    use Standard_Newton_Matrix_Series;

    order : integer32 := 1;
    info : integer32;
    tol : constant double_float := 1.0E-12;
    eva : Standard_Complex_Series_Vectors.Vector(f'range);

  begin
    if vrblvl > 0
     then put_line("-> in power_series_methods.Run_QR_Newton 10 ...");
    end if;
    if not verbose then
      QR_Newton_Steps(f,c,ejm,mlt,order,maxdeg,nbrit,s,info,vrblvl-1);
    else
      QR_Newton_Steps(file,f,c,ejm,mlt,order,maxdeg,nbrit,s,info,vrblvl-1);
      if info /= 0 then
        put(file,"info = "); put(file,info,1); new_line(file);
      else
        Complex_Series_and_Polynomials.Filter(s,tol);
        put_line(file,"The updated power series solution :");
        Complex_Series_and_Polynomials_io.put(file,s);
        eva := Standard_CSeries_Poly_SysFun.Eval(f,c,s);
        Complex_Series_and_Polynomials.Filter(eva,tol);
        put_line(file,"The evaluated solution :");
        Complex_Series_and_Polynomials_io.put(file,eva);
        Standard_Complex_Series_Vectors.Clear(eva);
      end if;
    end if;
  end Run_QR_Newton;

  procedure Run_QR_Newton
              ( maxdeg,nbrit : in integer32;
                f : in DoblDobl_CSeries_Poly_SysFun.Eval_Coeff_Poly_Sys;
                c : in DoblDobl_Complex_Series_VecVecs.VecVec;
                ejm : in DoblDobl_CSeries_Jaco_Matrices.Eval_Coeff_Jaco_Mat;
                mlt : in DoblDobl_CSeries_Jaco_Matrices.Mult_Factors;
                s : in out DoblDobl_Complex_Series_Vectors.Vector;
                verbose : in boolean := false;
                vrblvl : in integer32 := 0 ) is
  begin
    if vrblvl > 0
     then put_line("-> in power_series_methods.Run_QR_Newton 11 ...");
    end if;
    Run_QR_Newton(standard_output,maxdeg,nbrit,f,c,ejm,mlt,s,verbose,vrblvl-1);
  end Run_QR_Newton;

  procedure Run_QR_Newton
              ( file : in file_type; maxdeg,nbrit : in integer32;
                f : in DoblDobl_CSeries_Poly_SysFun.Eval_Coeff_Poly_Sys;
                c : in DoblDobl_Complex_Series_VecVecs.VecVec;
                ejm : in DoblDobl_CSeries_Jaco_Matrices.Eval_Coeff_Jaco_Mat;
                mlt : in DoblDobl_CSeries_Jaco_Matrices.Mult_Factors;
                s : in out DoblDobl_Complex_Series_Vectors.Vector;
                verbose : in boolean := false;
                vrblvl : in integer32 := 0 ) is

    use DoblDobl_Newton_Matrix_Series;

    order : integer32 := 1;
    info : integer32;
    tol : constant double_float := 1.0E-24;
    eva : DoblDobl_Complex_Series_Vectors.Vector(f'range);

  begin
    if vrblvl > 0
     then put_line("-> in power_series_methods.Run_QR_Newton 12 ...");
    end if;
    if not verbose then
      QR_Newton_Steps(f,c,ejm,mlt,order,maxdeg,nbrit,s,info,vrblvl-1);
    else
      QR_Newton_Steps(file,f,c,ejm,mlt,order,maxdeg,nbrit,s,info,vrblvl-1);
      if info /= 0 then
        put(file,"info = "); put(file,info,1); new_line(file);
      else
        Complex_Series_and_Polynomials.Filter(s,tol);
        put_line(file,"The updated power series solution :");
        Complex_Series_and_Polynomials_io.put(file,s);
        eva := DoblDobl_CSeries_Poly_SysFun.Eval(f,c,s);
        Complex_Series_and_Polynomials.Filter(eva,tol);
        put_line(file,"The evaluated solution :");
        Complex_Series_and_Polynomials_io.put(file,eva);
        DoblDobl_Complex_Series_Vectors.Clear(eva);
      end if;
    end if;
  end Run_QR_Newton;

  procedure Run_QR_Newton
              ( maxdeg,nbrit : in integer32;
                f : in TripDobl_CSeries_Poly_SysFun.Eval_Coeff_Poly_Sys;
                c : in TripDobl_Complex_Series_VecVecs.VecVec;
                ejm : in TripDobl_CSeries_Jaco_Matrices.Eval_Coeff_Jaco_Mat;
                mlt : in TripDobl_CSeries_Jaco_Matrices.Mult_Factors;
                s : in out TripDobl_Complex_Series_Vectors.Vector;
                verbose : in boolean := false;
                vrblvl : in integer32 := 0 ) is
  begin
    if vrblvl > 0
     then put_line("-> in power_series_methods.Run_QR_Newton 13 ...");
    end if;
    Run_QR_Newton(standard_output,maxdeg,nbrit,f,c,ejm,mlt,s,verbose,vrblvl-1);
  end Run_QR_Newton;

  procedure Run_QR_Newton
              ( file : in file_type; maxdeg,nbrit : in integer32;
                f : in TripDobl_CSeries_Poly_SysFun.Eval_Coeff_Poly_Sys;
                c : in TripDobl_Complex_Series_VecVecs.VecVec;
                ejm : in TripDobl_CSeries_Jaco_Matrices.Eval_Coeff_Jaco_Mat;
                mlt : in TripDobl_CSeries_Jaco_Matrices.Mult_Factors;
                s : in out TripDobl_Complex_Series_Vectors.Vector;
                verbose : in boolean := false;
                vrblvl : in integer32 := 0 ) is

    use TripDobl_Newton_Matrix_Series;

    order : integer32 := 1;
    info : integer32;
    tol : constant double_float := 1.0E-24;
    eva : TripDobl_Complex_Series_Vectors.Vector(f'range);

  begin
    if vrblvl > 0
     then put_line("-> in power_series_methods.Run_QR_Newton 14 ...");
    end if;
    if not verbose then
      QR_Newton_Steps(f,c,ejm,mlt,order,maxdeg,nbrit,s,info,vrblvl-1);
    else
      QR_Newton_Steps(file,f,c,ejm,mlt,order,maxdeg,nbrit,s,info,vrblvl-1);
      if info /= 0 then
        put(file,"info = "); put(file,info,1); new_line(file);
      else
        Complex_Series_and_Polynomials.Filter(s,tol);
        put_line(file,"The updated power series solution :");
        Complex_Series_and_Polynomials_io.put(file,s);
        eva := TripDobl_CSeries_Poly_SysFun.Eval(f,c,s);
        Complex_Series_and_Polynomials.Filter(eva,tol);
        put_line(file,"The evaluated solution :");
        Complex_Series_and_Polynomials_io.put(file,eva);
        TripDobl_Complex_Series_Vectors.Clear(eva);
      end if;
    end if;
  end Run_QR_Newton;

  procedure Run_QR_Newton
              ( maxdeg,nbrit : in integer32;
                f : in QuadDobl_CSeries_Poly_SysFun.Eval_Coeff_Poly_Sys;
                c : in QuadDobl_Complex_Series_VecVecs.VecVec;
                ejm : in QuadDobl_CSeries_Jaco_Matrices.Eval_Coeff_Jaco_Mat;
                mlt : in QuadDobl_CSeries_Jaco_Matrices.Mult_Factors;
                s : in out QuadDobl_Complex_Series_Vectors.Vector;
                verbose : in boolean := false;
                vrblvl : in integer32 := 0 ) is
  begin
    if vrblvl > 0
     then put_line("-> in power_series_methods.Run_QR_Newton 15 ...");
    end if;
    Run_QR_Newton(standard_output,maxdeg,nbrit,f,c,ejm,mlt,s,verbose,vrblvl-1);
  end Run_QR_Newton;

  procedure Run_QR_Newton
              ( file : in file_type; maxdeg,nbrit : in integer32;
                f : in QuadDobl_CSeries_Poly_SysFun.Eval_Coeff_Poly_Sys;
                c : in QuadDobl_Complex_Series_VecVecs.VecVec;
                ejm : in QuadDobl_CSeries_Jaco_Matrices.Eval_Coeff_Jaco_Mat;
                mlt : in QuadDobl_CSeries_Jaco_Matrices.Mult_Factors;
                s : in out QuadDobl_Complex_Series_Vectors.Vector;
                verbose : in boolean := false;
                vrblvl : in integer32 := 0 ) is

    use QuadDobl_Newton_Matrix_Series;

    order : integer32 := 1;
    info : integer32;
    tol : constant double_float := 1.0E-48;
    eva : QuadDobl_Complex_Series_Vectors.Vector(f'range);

  begin
    if vrblvl > 0
     then put_line("-> in power_series_methods.Run_QR_Newton 16 ...");
    end if;
    if not verbose then
      QR_Newton_Steps(f,c,ejm,mlt,order,maxdeg,nbrit,s,info,vrblvl-1);
    else
      QR_Newton_Steps(file,f,c,ejm,mlt,order,maxdeg,nbrit,s,info,vrblvl-1);
      if info /= 0 then
        put(file,"info = "); put(file,info,1); new_line(file);
      else
        Complex_Series_and_Polynomials.Filter(s,tol);
        put_line(file,"The updated power series solution :");
        Complex_Series_and_Polynomials_io.put(file,s);
        eva := QuadDobl_CSeries_Poly_SysFun.Eval(f,c,s);
        Complex_Series_and_Polynomials.Filter(eva,tol);
        put_line(file,"The evaluated solution :");
        Complex_Series_and_Polynomials_io.put(file,eva);
        QuadDobl_Complex_Series_Vectors.Clear(eva);
      end if;
    end if;
  end Run_QR_Newton;

-- SVD NEWTON :

  procedure Run_SVD_Newton
              ( maxdeg,nbrit : in integer32;
                p : in Standard_CSeries_Poly_Systems.Poly_Sys;
                s : in out Standard_Complex_Series_Vectors.Vector;
                rcond : out double_float; verbose : in boolean := false;
                vrblvl : in integer32 := 0 ) is
  begin
    if vrblvl > 0
     then put_line("-> in power_series_methods.Run_SVD_Newton 1 ...");
    end if;
    Run_SVD_Newton(standard_output,maxdeg,nbrit,p,s,rcond,verbose,vrblvl-1);
  end Run_SVD_Newton;

  procedure Run_SVD_Newton
              ( maxdeg,nbrit : in integer32;
                p : in Standard_CSeries_Poly_Systems.Poly_Sys;
                s : in out Standard_Complex_Series_Vectors.Vector;
                rcond : out double_float;
                evp : out Standard_Complex_Series_Vectors.Vector;
                verbose : in boolean := false;
                vrblvl : in integer32 := 0 ) is
  begin
    if vrblvl > 0
     then put_line("-> in power_series_methods.Run_SVD_Newton 2 ...");
    end if;
    Run_SVD_Newton(standard_output,maxdeg,nbrit,p,s,rcond,verbose,vrblvl-1);
    if 1.0 + rcond /= 1.0
     then evp := Standard_CSeries_Poly_SysFun.Eval(p,s);
    end if;
  end Run_SVD_Newton;

  procedure Run_SVD_Newton
              ( file : in file_type; maxdeg,nbrit : in integer32;
                p : in Standard_CSeries_Poly_Systems.Poly_Sys;
                s : in out Standard_Complex_Series_Vectors.Vector;
                rcond : out double_float; verbose : in boolean := false;
                vrblvl : in integer32 := 0 ) is

    use Standard_Newton_Matrix_Series; -- use Standard_Newton_Series;

    order : integer32 := 1;
    info : integer32;
    tol : constant double_float := 1.0E-12;
    eva : Standard_Complex_Series_Vectors.Vector(p'range);

  begin
    if vrblvl > 0
     then put_line("-> in power_series_methods.Run_SVD_Newton 3 ...");
    end if;
    if not verbose then
      SVD_Newton_Steps(p,order,maxdeg,nbrit,s,info,rcond,vrblvl-1);
    else
      SVD_Newton_Steps(file,p,order,maxdeg,nbrit,s,info,rcond,vrblvl-1);
      put(file,"rcond = "); put(file,rcond,3); new_line(file);
      if 1.0 + rcond /= 1.0 then
        Complex_Series_and_Polynomials.Filter(s,tol);
        put_line(file,"The updated power series solution :");
        Complex_Series_and_Polynomials_io.put(file,s);
        eva := Standard_CSeries_Poly_SysFun.Eval(p,s);
        Complex_Series_and_Polynomials.Filter(eva,tol);
        put_line(file,"The evaluated solution :");
        Complex_Series_and_Polynomials_io.put(file,eva);
        Standard_Complex_Series_Vectors.Clear(eva);
      end if;
    end if;
  end Run_SVD_Newton;

  procedure Run_SVD_Newton
              ( file : in file_type; maxdeg,nbrit : in integer32;
                p : in Standard_CSeries_Poly_Systems.Poly_Sys;
                s : in out Standard_Complex_Series_Vectors.Vector;
                rcond : out double_float;
                evp : out Standard_Complex_Series_Vectors.Vector;
                verbose : in boolean := false;
                vrblvl : in integer32 := 0 ) is

    use Standard_Newton_Matrix_Series; -- use Standard_Newton_Series;

    order : integer32 := 1;
    info : integer32;
    tol : constant double_float := 1.0E-12;

  begin
    if vrblvl > 0
     then put_line("-> in power_series_methods.Run_SVD_Newton 4 ...");
    end if;
    if not verbose then
      Run_SVD_Newton(maxdeg,nbrit,p,s,rcond,evp,verbose,vrblvl-1);
    else
      SVD_Newton_Steps(file,p,order,maxdeg,nbrit,s,info,rcond,vrblvl-1);
      put(file,"rcond = "); put(file,rcond,3); new_line(file);
      if 1.0 + rcond /= 1.0 then
        Complex_Series_and_Polynomials.Filter(s,tol);
        put_line(file,"The updated power series solution :");
        Complex_Series_and_Polynomials_io.put(file,s);
        evp := Standard_CSeries_Poly_SysFun.Eval(p,s);
        Complex_Series_and_Polynomials.Filter(evp,tol);
        put_line(file,"The evaluated solution :");
        Complex_Series_and_Polynomials_io.put(file,evp);
      end if;
    end if;
  end Run_SVD_Newton;

  procedure Run_SVD_Newton
             ( maxdeg,nbrit : in integer32;
               p : in DoblDobl_CSeries_Poly_Systems.Poly_Sys;
               s : in out DoblDobl_Complex_Series_Vectors.Vector;
               rcond : out double_double; verbose : in boolean := false;
               vrblvl : in integer32 := 0 ) is
  begin
    if vrblvl > 0
     then put_line("-> in power_series_methods.Run_SVD_Newton 5 ...");
    end if;
    Run_SVD_Newton(standard_output,maxdeg,nbrit,p,s,rcond,verbose,vrblvl-1);
  end Run_SVD_Newton;

  procedure Run_SVD_Newton
             ( file : in file_type; maxdeg,nbrit : in integer32;
               p : in DoblDobl_CSeries_Poly_Systems.Poly_Sys;
               s : in out DoblDobl_Complex_Series_Vectors.Vector;
               rcond : out double_double; verbose : in boolean := false;
               vrblvl : in integer32 := 0 ) is

    use DoblDobl_Newton_Matrix_Series; -- use DoblDobl_Newton_Series;

    order : integer32 := 1;
    info : integer32;
    tol : constant double_float := 1.0E-24;
    eva : DoblDobl_Complex_Series_Vectors.Vector(p'range);
    one : constant double_double := create(1.0);

  begin
    if vrblvl > 0
     then put_line("-> in power_series_methods.Run_SVD_Newton 6 ...");
    end if;
    if not verbose then
      SVD_Newton_Steps(p,order,maxdeg,nbrit,s,info,rcond,vrblvl-1);
    else
      SVD_Newton_Steps(file,p,order,maxdeg,nbrit,s,info,rcond,vrblvl-1);
      put(file,"rcond = "); put(file,rcond,3); new_line(file);
      if one + rcond /= one then 
        Complex_Series_and_Polynomials.Filter(s,tol);
        put_line(file,"The updated power series solution :");
        Complex_Series_and_Polynomials_io.put(file,s);
        eva := DoblDobl_CSeries_Poly_SysFun.Eval(p,s);
        Complex_Series_and_Polynomials.Filter(eva,tol);
        put_line(file,"The evaluated solution :");
        Complex_Series_and_Polynomials_io.put(file,eva);
        DoblDobl_Complex_Series_Vectors.Clear(eva);
      end if;
    end if;
  end Run_SVD_Newton;

  procedure Run_SVD_Newton
             ( maxdeg,nbrit : in integer32;
               p : in TripDobl_CSeries_Poly_Systems.Poly_Sys;
               s : in out TripDobl_Complex_Series_Vectors.Vector;
               rcond : out triple_double; verbose : in boolean := false;
               vrblvl : in integer32 := 0 ) is
  begin
    if vrblvl > 0
     then put_line("-> in power_series_methods.Run_SVD_Newton 7 ...");
    end if;
    Run_SVD_Newton(standard_output,maxdeg,nbrit,p,s,rcond,verbose,vrblvl-1);
  end Run_SVD_Newton;

  procedure Run_SVD_Newton
             ( file : in file_type; maxdeg,nbrit : in integer32;
               p : in TripDobl_CSeries_Poly_Systems.Poly_Sys;
               s : in out TripDobl_Complex_Series_Vectors.Vector;
               rcond : out triple_double; verbose : in boolean := false;
               vrblvl : in integer32 := 0 ) is

    use TripDobl_Newton_Matrix_Series; -- use TripDobl_Newton_Series;

    order : integer32 := 1;
    info : integer32;
    tol : constant double_float := 1.0E-24;
    eva : TripDobl_Complex_Series_Vectors.Vector(p'range);
    one : constant triple_double := create(1.0);

  begin
    if vrblvl > 0
     then put_line("-> in power_series_methods.Run_SVD_Newton 8 ...");
    end if;
    if not verbose then
      SVD_Newton_Steps(p,order,maxdeg,nbrit,s,info,rcond,vrblvl-1);
    else
      SVD_Newton_Steps(file,p,order,maxdeg,nbrit,s,info,rcond,vrblvl-1);
      put(file,"rcond = "); put(file,rcond,3); new_line(file);
      if one + rcond /= one then 
        Complex_Series_and_Polynomials.Filter(s,tol);
        put_line(file,"The updated power series solution :");
        Complex_Series_and_Polynomials_io.put(file,s);
        eva := TripDobl_CSeries_Poly_SysFun.Eval(p,s);
        Complex_Series_and_Polynomials.Filter(eva,tol);
        put_line(file,"The evaluated solution :");
        Complex_Series_and_Polynomials_io.put(file,eva);
        TripDobl_Complex_Series_Vectors.Clear(eva);
      end if;
    end if;
  end Run_SVD_Newton;

  procedure Run_SVD_Newton
             ( maxdeg,nbrit : in integer32;
               p : in QuadDobl_CSeries_Poly_Systems.Poly_Sys;
               s : in out QuadDobl_Complex_Series_Vectors.Vector;
               rcond : out quad_double; verbose : in boolean := false;
               vrblvl : in integer32 := 0 ) is
  begin
    if vrblvl > 0
     then put_line("-> in power_series_methods.Run_SVD_Newton 9 ...");
    end if;
    Run_SVD_Newton(standard_output,maxdeg,nbrit,p,s,rcond,verbose,vrblvl-1);
  end Run_SVD_Newton;

  procedure Run_SVD_Newton
             ( file : in file_type; maxdeg,nbrit : in integer32;
               p : in QuadDobl_CSeries_Poly_Systems.Poly_Sys;
               s : in out QuadDobl_Complex_Series_Vectors.Vector;
               rcond : out quad_double; verbose : in boolean := false;
               vrblvl : in integer32 := 0 ) is

    use QuadDobl_Newton_Matrix_Series; -- use QuadDobl_Newton_Series;

    order : integer32 := 1;
    info : integer32;
    tol : constant double_float := 1.0E-32;
    eva : QuadDobl_Complex_Series_Vectors.Vector(p'range);
    one : constant quad_double := create(1.0);

  begin
    if vrblvl > 0
     then put_line("-> in power_series_methods.Run_SVD_Newton 10 ...");
    end if;
    if not verbose then
      SVD_Newton_Steps(p,order,maxdeg,nbrit,s,info,rcond,vrblvl-1);
    else
      SVD_Newton_Steps(file,p,order,maxdeg,nbrit,s,info,rcond,vrblvl-1);
      put(file,"rcond = "); put(file,rcond,3); new_line(file);
      if one + rcond /= one then
        Complex_Series_and_Polynomials.Filter(s,tol);
        put_line(file,"The updated power series solution :");
        Complex_Series_and_Polynomials_io.put(file,s);
        eva := QuadDobl_CSeries_Poly_SysFun.Eval(p,s);
        Complex_Series_and_Polynomials.Filter(eva,tol);
        put_line(file,"The evaluated solution :");
        Complex_Series_and_Polynomials_io.put(file,eva);
        QuadDobl_Complex_Series_Vectors.Clear(eva);
      end if;
    end if;
  end Run_SVD_Newton;

  procedure Run_SVD_Newton
             ( maxdeg,nbrit : in integer32;
               p : in PentDobl_CSeries_Poly_Systems.Poly_Sys;
               s : in out PentDobl_Complex_Series_Vectors.Vector;
               rcond : out penta_double; verbose : in boolean := false;
               vrblvl : in integer32 := 0 ) is
  begin
    if vrblvl > 0
     then put_line("-> in power_series_methods.Run_SVD_Newton 11 ...");
    end if;
    Run_SVD_Newton(standard_output,maxdeg,nbrit,p,s,rcond,verbose,vrblvl-1);
  end Run_SVD_Newton;

  procedure Run_SVD_Newton
             ( file : in file_type; maxdeg,nbrit : in integer32;
               p : in PentDobl_CSeries_Poly_Systems.Poly_Sys;
               s : in out PentDobl_Complex_Series_Vectors.Vector;
               rcond : out penta_double; verbose : in boolean := false;
               vrblvl : in integer32 := 0 ) is

    use PentDobl_Newton_Matrix_Series; -- use PentDobl_Newton_Series;

    order : integer32 := 1;
    info : integer32;
    tol : constant double_float := 1.0E-40;
    eva : PentDobl_Complex_Series_Vectors.Vector(p'range);
    one : constant penta_double := create(1.0);

  begin
    if vrblvl > 0
     then put_line("-> in power_series_methods.Run_SVD_Newton 12 ...");
    end if;
    if not verbose then
      SVD_Newton_Steps(p,order,maxdeg,nbrit,s,info,rcond,vrblvl-1);
    else
      SVD_Newton_Steps(file,p,order,maxdeg,nbrit,s,info,rcond,vrblvl-1);
      put(file,"rcond = "); put(file,rcond,3); new_line(file);
      if one + rcond /= one then
        Complex_Series_and_Polynomials.Filter(s,tol);
        put_line(file,"The updated power series solution :");
        Complex_Series_and_Polynomials_io.put(file,s);
        eva := PentDobl_CSeries_Poly_SysFun.Eval(p,s);
        Complex_Series_and_Polynomials.Filter(eva,tol);
        put_line(file,"The evaluated solution :");
        Complex_Series_and_Polynomials_io.put(file,eva);
        PentDobl_Complex_Series_Vectors.Clear(eva);
      end if;
    end if;
  end Run_SVD_Newton;

  procedure Run_SVD_Newton
             ( maxdeg,nbrit : in integer32;
               p : in OctoDobl_CSeries_Poly_Systems.Poly_Sys;
               s : in out OctoDobl_Complex_Series_Vectors.Vector;
               rcond : out octo_double; verbose : in boolean := false;
               vrblvl : in integer32 := 0 ) is
  begin
    if vrblvl > 0
     then put_line("-> in power_series_methods.Run_SVD_Newton 13 ...");
    end if;
    Run_SVD_Newton(standard_output,maxdeg,nbrit,p,s,rcond,verbose,vrblvl-1);
  end Run_SVD_Newton;

  procedure Run_SVD_Newton
             ( file : in file_type; maxdeg,nbrit : in integer32;
               p : in OctoDobl_CSeries_Poly_Systems.Poly_Sys;
               s : in out OctoDobl_Complex_Series_Vectors.Vector;
               rcond : out octo_double; verbose : in boolean := false;
               vrblvl : in integer32 := 0 ) is

    use OctoDobl_Newton_Matrix_Series; -- use OctoDobl_Newton_Series;

    order : integer32 := 1;
    info : integer32;
    tol : constant double_float := 1.0E-40;
    eva : OctoDobl_Complex_Series_Vectors.Vector(p'range);
    one : constant octo_double := create(1.0);

  begin
    if vrblvl > 0
     then put_line("-> in power_series_methods.Run_SVD_Newton 14 ...");
    end if;
    if not verbose then
      SVD_Newton_Steps(p,order,maxdeg,nbrit,s,info,rcond,vrblvl-1);
    else
      SVD_Newton_Steps(file,p,order,maxdeg,nbrit,s,info,rcond,vrblvl-1);
      put(file,"rcond = "); put(file,rcond,3); new_line(file);
      if one + rcond /= one then
        Complex_Series_and_Polynomials.Filter(s,tol);
        put_line(file,"The updated power series solution :");
        Complex_Series_and_Polynomials_io.put(file,s);
        eva := OctoDobl_CSeries_Poly_SysFun.Eval(p,s);
        Complex_Series_and_Polynomials.Filter(eva,tol);
        put_line(file,"The evaluated solution :");
        Complex_Series_and_Polynomials_io.put(file,eva);
        OctoDobl_Complex_Series_Vectors.Clear(eva);
      end if;
    end if;
  end Run_SVD_Newton;

  procedure Run_SVD_Newton
             ( maxdeg,nbrit : in integer32;
               p : in DecaDobl_CSeries_Poly_Systems.Poly_Sys;
               s : in out DecaDobl_Complex_Series_Vectors.Vector;
               rcond : out deca_double; verbose : in boolean := false;
               vrblvl : in integer32 := 0 ) is
  begin
    if vrblvl > 0
     then put_line("-> in power_series_methods.Run_SVD_Newton 15 ...");
    end if;
    Run_SVD_Newton(standard_output,maxdeg,nbrit,p,s,rcond,verbose,vrblvl-1);
  end Run_SVD_Newton;

  procedure Run_SVD_Newton
             ( file : in file_type; maxdeg,nbrit : in integer32;
               p : in DecaDobl_CSeries_Poly_Systems.Poly_Sys;
               s : in out DecaDobl_Complex_Series_Vectors.Vector;
               rcond : out deca_double; verbose : in boolean := false;
               vrblvl : in integer32 := 0 ) is

    use DecaDobl_Newton_Matrix_Series; -- use DecaDobl_Newton_Series;

    order : integer32 := 1;
    info : integer32;
    tol : constant double_float := 1.0E-40;
    eva : DecaDobl_Complex_Series_Vectors.Vector(p'range);
    one : constant deca_double := create(1.0);

  begin
    if vrblvl > 0
     then put_line("-> in power_series_methods.Run_SVD_Newton 16 ...");
    end if;
    if not verbose then
      SVD_Newton_Steps(p,order,maxdeg,nbrit,s,info,rcond,vrblvl-1);
    else
      SVD_Newton_Steps(file,p,order,maxdeg,nbrit,s,info,rcond,vrblvl-1);
      put(file,"rcond = "); put(file,rcond,3); new_line(file);
      if one + rcond /= one then
        Complex_Series_and_Polynomials.Filter(s,tol);
        put_line(file,"The updated power series solution :");
        Complex_Series_and_Polynomials_io.put(file,s);
        eva := DecaDobl_CSeries_Poly_SysFun.Eval(p,s);
        Complex_Series_and_Polynomials.Filter(eva,tol);
        put_line(file,"The evaluated solution :");
        Complex_Series_and_Polynomials_io.put(file,eva);
        DecaDobl_Complex_Series_Vectors.Clear(eva);
      end if;
    end if;
  end Run_SVD_Newton;

  procedure Run_Echelon_Newton
              ( maxdeg,nbrit : in integer32;
                p : in Standard_CSeries_Poly_Systems.Poly_Sys;
                s : in out Standard_Complex_Series_Vectors.Vector;
                det : out Standard_Complex_Numbers.Complex_Number;
                verbose : in boolean := false;
                vrblvl : in integer32 := 0 ) is
  begin
    if vrblvl > 0
     then put_line("-> in power_series_methods.Run_Echelon_Newton 1 ...");
    end if;
    Run_Echelon_Newton(standard_output,maxdeg,nbrit,p,s,det,verbose,vrblvl-1);
  end Run_Echelon_Newton;

  procedure Run_Echelon_Newton
              ( file : in file_type; maxdeg,nbrit : in integer32;
                p : in Standard_CSeries_Poly_Systems.Poly_Sys;
                s : in out Standard_Complex_Series_Vectors.Vector;
                det : out Standard_Complex_Numbers.Complex_Number;
                verbose : in boolean := false;
                vrblvl : in integer32 := 0 ) is

    use Standard_Newton_Matrix_Series;

    order : integer32 := 1;
    tol : constant double_float := 1.0E-12;
    eva : Standard_Complex_Series_Vectors.Vector(p'range);

  begin
    if vrblvl > 0
     then put_line("-> in power_series_methods.Run_Echelon_Newton 2 ...");
    end if;
    if not verbose then
      Echelon_Newton_Steps(p,order,maxdeg,nbrit,s,det,vrblvl-1);
    else
      Echelon_Newton_Steps(file,p,order,maxdeg,nbrit,s,det,vrblvl-1);
      put(file,"det : "); put(file,det); new_line(file);
      Complex_Series_and_Polynomials.Filter(s,tol);
      put_line(file,"The updated power series solution :");
      Complex_Series_and_Polynomials_io.put(file,s);
      eva := Standard_CSeries_Poly_SysFun.Eval(p,s);
      Complex_Series_and_Polynomials.Filter(eva,tol);
      put_line(file,"The evaluated solution :");
      Complex_Series_and_Polynomials_io.put(file,eva);
      Standard_Complex_Series_Vectors.Clear(eva);
    end if;
  exception
    when others =>
      put_line("exception in Power_Series_Methods.Run_Echelon_Newton");
      raise;
  end Run_Echelon_Newton;

  procedure Run_Echelon_Newton
              ( maxdeg,nbrit : in integer32;
                p : in DoblDobl_CSeries_Poly_Systems.Poly_Sys;
                s : in out DoblDobl_Complex_Series_Vectors.Vector;
                det : out DoblDobl_Complex_Numbers.Complex_Number;
                verbose : in boolean := false;
                vrblvl : in integer32 := 0 ) is
  begin
    if vrblvl > 0
     then put_line("-> in power_series_methods.Run_Echelon_Newton 3 ...");
    end if;
    Run_Echelon_Newton(standard_output,maxdeg,nbrit,p,s,det,verbose,vrblvl-1);
  end Run_Echelon_Newton;

  procedure Run_Echelon_Newton
              ( file : in file_type; maxdeg,nbrit : in integer32;
                p : in DoblDobl_CSeries_Poly_Systems.Poly_Sys;
                s : in out DoblDobl_Complex_Series_Vectors.Vector;
                det : out DoblDobl_Complex_Numbers.Complex_Number;
                verbose : in boolean := false;
                vrblvl : in integer32 := 0 ) is

    use DoblDobl_Newton_Matrix_Series;

    order : integer32 := 1;
    tol : constant double_float := 1.0E-12;
    eva : DoblDobl_Complex_Series_Vectors.Vector(p'range);

  begin
    if vrblvl > 0
     then put_line("-> in power_series_methods.Run_Echelon_Newton 4 ...");
    end if;
    if not verbose then
      Echelon_Newton_Steps(p,order,maxdeg,nbrit,s,det,vrblvl-1);
    else
      Echelon_Newton_Steps(file,p,order,maxdeg,nbrit,s,det,vrblvl-1);
      put(file,"det : "); put(file,det); new_line(file);
      Complex_Series_and_Polynomials.Filter(s,tol);
      put_line(file,"The updated power series solution :");
      Complex_Series_and_Polynomials_io.put(file,s);
      eva := DoblDobl_CSeries_Poly_SysFun.Eval(p,s);
      Complex_Series_and_Polynomials.Filter(eva,tol);
      put_line(file,"The evaluated solution :");
      Complex_Series_and_Polynomials_io.put(file,eva);
      DoblDobl_Complex_Series_Vectors.Clear(eva);
    end if;
  end Run_Echelon_Newton;

  procedure Run_Echelon_Newton
              ( maxdeg,nbrit : in integer32;
                p : in TripDobl_CSeries_Poly_Systems.Poly_Sys;
                s : in out TripDobl_Complex_Series_Vectors.Vector;
                det : out TripDobl_Complex_Numbers.Complex_Number;
                verbose : in boolean := false;
                vrblvl : in integer32 := 0 ) is
  begin
    if vrblvl > 0
     then put_line("-> in power_series_methods.Run_Echelon_Newton 5 ...");
    end if;
    Run_Echelon_Newton(standard_output,maxdeg,nbrit,p,s,det,verbose,vrblvl-1);
  end Run_Echelon_Newton;

  procedure Run_Echelon_Newton
              ( file : in file_type; maxdeg,nbrit : in integer32;
                p : in TripDobl_CSeries_Poly_Systems.Poly_Sys;
                s : in out TripDobl_Complex_Series_Vectors.Vector;
                det : out TripDobl_Complex_Numbers.Complex_Number;
                verbose : in boolean := false;
                vrblvl : in integer32 := 0 ) is

    use TripDobl_Newton_Matrix_Series;

    order : integer32 := 1;
    tol : constant double_float := 1.0E-12;
    eva : TripDobl_Complex_Series_Vectors.Vector(p'range);

  begin
    if vrblvl > 0
     then put_line("-> in power_series_methods.Run_Echelon_Newton 6 ...");
    end if;
    if not verbose then
      Echelon_Newton_Steps(p,order,maxdeg,nbrit,s,det,vrblvl-1);
    else
      Echelon_Newton_Steps(file,p,order,maxdeg,nbrit,s,det,vrblvl-1);
      put(file,"det : "); put(file,det); new_line(file);
      Complex_Series_and_Polynomials.Filter(s,tol);
      put_line(file,"The updated power series solution :");
      Complex_Series_and_Polynomials_io.put(file,s);
      eva := TripDobl_CSeries_Poly_SysFun.Eval(p,s);
      Complex_Series_and_Polynomials.Filter(eva,tol);
      put_line(file,"The evaluated solution :");
      Complex_Series_and_Polynomials_io.put(file,eva);
      TripDobl_Complex_Series_Vectors.Clear(eva);
    end if;
  end Run_Echelon_Newton;

  procedure Run_Echelon_Newton
              ( maxdeg,nbrit : in integer32;
                p : in QuadDobl_CSeries_Poly_Systems.Poly_Sys;
                s : in out QuadDobl_Complex_Series_Vectors.Vector;
                det : out QuadDobl_Complex_Numbers.Complex_Number;
                verbose : in boolean := false;
                vrblvl : in integer32 := 0 ) is
  begin
    if vrblvl > 0
     then put_line("-> in power_series_methods.Run_Echelon_Newton 7 ...");
    end if;
    Run_Echelon_Newton(standard_output,maxdeg,nbrit,p,s,det,verbose,vrblvl-1);
  end Run_Echelon_Newton;

  procedure Run_Echelon_Newton
              ( file : in file_type; maxdeg,nbrit : in integer32;
                p : in QuadDobl_CSeries_Poly_Systems.Poly_Sys;
                s : in out QuadDobl_Complex_Series_Vectors.Vector;
                det : out QuadDobl_Complex_Numbers.Complex_Number;
                verbose : in boolean := false;
                vrblvl : in integer32 := 0 ) is

    use QuadDobl_Newton_Matrix_Series;

    order : integer32 := 1;
    tol : constant double_float := 1.0E-12;
    eva : QuadDobl_Complex_Series_Vectors.Vector(p'range);

  begin
    if vrblvl > 0
     then put_line("-> in power_series_methods.Run_Echelon_Newton 8 ...");
    end if;
    if not verbose then
      Echelon_Newton_Steps(p,order,maxdeg,nbrit,s,det,vrblvl-1);
    else
      Echelon_Newton_Steps(file,p,order,maxdeg,nbrit,s,det,vrblvl-1);
      put(file,"det : "); put(file,det); new_line(file);
      Complex_Series_and_Polynomials.Filter(s,tol);
      put_line(file,"The updated power series solution :");
      Complex_Series_and_Polynomials_io.put(file,s);
      eva := QuadDobl_CSeries_Poly_SysFun.Eval(p,s);
      Complex_Series_and_Polynomials.Filter(eva,tol);
      put_line(file,"The evaluated solution :");
      Complex_Series_and_Polynomials_io.put(file,eva);
      QuadDobl_Complex_Series_Vectors.Clear(eva);
    end if;
  end Run_Echelon_Newton;

  procedure Run_LU_Newton
              ( maxdeg,nbrit : in integer32;
                p : in Standard_CSeries_Poly_Systems.Poly_Sys;
                v : in Standard_Complex_Series_VecVecs.VecVec;
                verbose : in boolean := false;
                pause : in boolean := false;
                vrblvl : in integer32 := 0 ) is

    ans : character;
    info : integer32;

  begin
    if vrblvl > 0
     then put_line("-> in power_series_methods.Run_LU_Newton 25 ...");
    end if;
    for i in v'range loop
      if verbose then
        put("Running on solution "); put(i,1); put_line(" ...");
        declare
          lvi : constant Standard_Complex_Series_Vectors.Link_to_Vector := v(i);
          sol : Standard_Complex_Vectors.Vector(lvi'range);
        begin
          for k in sol'range loop
            sol(k) := lvi(k).cff(0);
          end loop;
          Standard_Complex_Solutions_io.put_vector(sol);
        end;
      end if;
      Run_LU_Newton(maxdeg,nbrit,p,v(i).all,info,verbose,vrblvl-1);
      if verbose then
        if pause then
          put("Continue ? (y/n) "); Ask_Yes_or_No(ans);
          exit when (ans /= 'y');
        end if;
      end if;
    end loop;
  end Run_LU_Newton;

  procedure Run_LU_Newton
              ( maxdeg,nbrit : in integer32;
                p : in DoblDobl_CSeries_Poly_Systems.Poly_Sys;
                v : in DoblDobl_Complex_Series_VecVecs.VecVec;
                verbose : in boolean := false;
                pause : in boolean := false;
                vrblvl : in integer32 := 0 ) is

    ans : character;
    info : integer32;

  begin
    if vrblvl > 0
     then put_line("-> in power_series_methods.Run_LU_Newton 26 ...");
    end if;
    for i in v'range loop
      if verbose then
        put("Running on solution "); put(i,1); put_line(" ...");
        declare
          lvi : constant DoblDobl_Complex_Series_Vectors.Link_to_Vector := v(i);
          sol : DoblDobl_Complex_Vectors.Vector(lvi'range);
        begin
          for k in sol'range loop
            sol(k) := lvi(k).cff(0);
          end loop;
          DoblDobl_Complex_Solutions_io.put_vector(sol);
        end;
      end if;
      Run_LU_Newton(maxdeg,nbrit,p,v(i).all,info,verbose,vrblvl-1);
      if verbose then
        if pause then
          put("Continue ? (y/n) "); Ask_Yes_or_No(ans);
          exit when (ans /= 'y');
        end if;
      end if;
    end loop;
  end Run_LU_Newton;

  procedure Run_LU_Newton
              ( maxdeg,nbrit : in integer32;
                p : in TripDobl_CSeries_Poly_Systems.Poly_Sys;
                v : in TripDobl_Complex_Series_VecVecs.VecVec;
                verbose : in boolean := false;
                pause : in boolean := false;
                vrblvl : in integer32 := 0 ) is

    ans : character;
    info : integer32;

  begin
    if vrblvl > 0
     then put_line("-> in power_series_methods.Run_LU_Newton 27 ...");
    end if;
    for i in v'range loop
      if verbose then
        put("Running on solution "); put(i,1); put_line(" ...");
        declare
          lvi : constant TripDobl_Complex_Series_Vectors.Link_to_Vector := v(i);
          sol : TripDobl_Complex_Vectors.Vector(lvi'range);
        begin
          for k in sol'range loop
            sol(k) := lvi(k).cff(0);
          end loop;
          TripDobl_Complex_Solutions_io.put_vector(sol);
        end;
      end if;
      Run_LU_Newton(maxdeg,nbrit,p,v(i).all,info,verbose,vrblvl-1);
      if verbose then
        if pause then
          put("Continue ? (y/n) "); Ask_Yes_or_No(ans);
          exit when (ans /= 'y');
        end if;
      end if;
    end loop;
  end Run_LU_Newton;

  procedure Run_LU_Newton
              ( maxdeg,nbrit : in integer32;
                p : in QuadDobl_CSeries_Poly_Systems.Poly_Sys;
                v : in QuadDobl_Complex_Series_VecVecs.VecVec;
                verbose : in boolean := false;
                pause : in boolean := false;
                vrblvl : in integer32 := 0 ) is

    ans : character;
    info : integer32;

  begin
    if vrblvl > 0
     then put_line("-> in power_series_methods.Run_LU_Newton 28 ...");
    end if;
    for i in v'range loop
      if verbose then
        put("Running on solution "); put(i,1); put_line(" ...");
        declare
          lvi : constant QuadDobl_Complex_Series_Vectors.Link_to_Vector := v(i);
          sol : QuadDobl_Complex_Vectors.Vector(lvi'range);
        begin
          for k in sol'range loop
            sol(k) := lvi(k).cff(0);
          end loop;
          QuadDobl_Complex_Solutions_io.put_vector(sol);
        end;
      end if;
      Run_LU_Newton(maxdeg,nbrit,p,v(i).all,info,verbose,vrblvl-1);
      if verbose then
        if pause then
          put("Continue ? (y/n) "); Ask_Yes_or_No(ans);
          exit when (ans /= 'y');
        end if;
      end if;
    end loop;
  end Run_LU_Newton;

  procedure Run_LU_Newton
              ( file : in file_type; maxdeg,nbrit : in integer32;
                p : in Standard_CSeries_Poly_Systems.Poly_Sys;
                v : in Standard_Complex_Series_VecVecs.VecVec;
                verbose : in boolean := false;
                vrblvl : in integer32 := 0 ) is

    info : integer32;

  begin
    if vrblvl > 0
     then put_line("-> in power_series_methods.Run_LU_Newton 29 ...");
    end if;
    for i in v'range loop
      if verbose then
        put(file,"Running on solution ");
        put(file,i,1); put_line(file," ...");
        declare
          lvi : constant Standard_Complex_Series_Vectors.Link_to_Vector := v(i);
          sol : Standard_Complex_Vectors.Vector(lvi'range);
        begin
          for k in sol'range loop
            sol(k) := lvi(k).cff(0);
          end loop;
          Standard_Complex_Solutions_io.put_vector(file,sol);
        end;
      end if;
      Run_LU_Newton(file,maxdeg,nbrit,p,v(i).all,info,verbose,vrblvl-1);
    end loop;
  end Run_LU_Newton;

  procedure Run_LU_Newton
              ( file : in file_type; maxdeg,nbrit : in integer32;
                p : in DoblDobl_CSeries_Poly_Systems.Poly_Sys;
                v : in DoblDobl_Complex_Series_VecVecs.VecVec;
                verbose : in boolean := false;
                vrblvl : in integer32 := 0 ) is

    info : integer32;

  begin
    if vrblvl > 0
     then put_line("-> in power_series_methods.Run_LU_Newton 30 ...");
    end if;
    for i in v'range loop
      if verbose then
        put(file,"Running on solution ");
        put(file,i,1); put_line(file," ...");
        declare
          lvi : constant DoblDobl_Complex_Series_Vectors.Link_to_Vector := v(i);
          sol : DoblDobl_Complex_Vectors.Vector(lvi'range);
        begin
          for k in sol'range loop
            sol(k) := lvi(k).cff(0);
          end loop;
          DoblDobl_Complex_Solutions_io.put_vector(file,sol);
        end;
      end if;
      Run_LU_Newton(file,maxdeg,nbrit,p,v(i).all,info,verbose,vrblvl-1);
    end loop;
  end Run_LU_Newton;

  procedure Run_LU_Newton
              ( file : in file_type; maxdeg,nbrit : in integer32;
                p : in TripDobl_CSeries_Poly_Systems.Poly_Sys;
                v : in TripDobl_Complex_Series_VecVecs.VecVec;
                verbose : in boolean := false;
                vrblvl : in integer32 := 0 ) is

    info : integer32;

  begin
    if vrblvl > 0
     then put_line("-> in power_series_methods.Run_LU_Newton 31 ...");
    end if;
    for i in v'range loop
      if verbose then
        put(file,"Running on solution ");
        put(file,i,1); put_line(file," ...");
        declare
          lvi : constant TripDobl_Complex_Series_Vectors.Link_to_Vector := v(i);
          sol : TripDobl_Complex_Vectors.Vector(lvi'range);
        begin
          for k in sol'range loop
            sol(k) := lvi(k).cff(0);
          end loop;
          TripDobl_Complex_Solutions_io.put_vector(file,sol);
        end;
      end if;
      Run_LU_Newton(file,maxdeg,nbrit,p,v(i).all,info,verbose,vrblvl-1);
    end loop;
  end Run_LU_Newton;

  procedure Run_LU_Newton
              ( file : in file_type; maxdeg,nbrit : in integer32;
                p : in QuadDobl_CSeries_Poly_Systems.Poly_Sys;
                v : in QuadDobl_Complex_Series_VecVecs.VecVec;
                verbose : in boolean := false;
                vrblvl : in integer32 := 0 ) is

    info : integer32;

  begin
    if vrblvl > 0
     then put_line("-> in power_series_methods.Run_LU_Newton 32 ...");
    end if;
    for i in v'range loop
      if verbose then
        put(file,"Running on solution ");
        put(file,i,1); put_line(file," ...");
        declare
          lvi : constant QuadDobl_Complex_Series_Vectors.Link_to_Vector := v(i);
          sol : QuadDobl_Complex_Vectors.Vector(lvi'range);
        begin
          for k in sol'range loop
            sol(k) := lvi(k).cff(0);
          end loop;
          QuadDobl_Complex_Solutions_io.put_vector(file,sol);
        end;
      end if;
      Run_LU_Newton(file,maxdeg,nbrit,p,v(i).all,info,verbose,vrblvl-1);
    end loop;
  end Run_LU_Newton;

  procedure Run_QR_Newton
              ( maxdeg,nbrit : in integer32;
                p : in Standard_CSeries_Poly_Systems.Poly_Sys;
                v : in Standard_Complex_Series_VecVecs.VecVec;
                verbose : in boolean := false;
                pause : in boolean := false;
                vrblvl : in integer32 := 0 ) is

    ans : character;

  begin
    if vrblvl > 0
     then put_line("-> in power_series_methods.Run_QR_Newton 17 ...");
    end if;
    for i in v'range loop
      if verbose then
        put("Running on solution "); put(i,1); put_line(" ...");
        declare
          lvi : constant Standard_Complex_Series_Vectors.Link_to_Vector := v(i);
          sol : Standard_Complex_Vectors.Vector(lvi'range);
        begin
          for k in sol'range loop
            sol(k) := lvi(k).cff(0);
          end loop;
          Standard_Complex_Solutions_io.put_vector(sol);
        end;
      end if;
      Run_QR_Newton(maxdeg,nbrit,p,v(i).all,verbose,vrblvl-1);
      if verbose then
        if pause then
          put("Continue ? (y/n) "); Ask_Yes_or_No(ans);
          exit when (ans /= 'y');
        end if;
      end if;
    end loop;
  end Run_QR_Newton;

  procedure Run_QR_Newton
              ( maxdeg,nbrit : in integer32;
                p : in DoblDobl_CSeries_Poly_Systems.Poly_Sys;
                v : in DoblDobl_Complex_Series_VecVecs.VecVec;
                verbose : in boolean := false;
                pause : in boolean := false;
                vrblvl : in integer32 := 0 ) is

    ans : character;

  begin
    if vrblvl > 0
     then put_line("-> in power_series_methods.Run_QR_Newton 18 ...");
    end if;
    for i in v'range loop
      if verbose then
        put("Running on solution "); put(i,1); put_line(" ...");
        declare
          lvi : constant DoblDobl_Complex_Series_Vectors.Link_to_Vector := v(i);
          sol : DoblDobl_Complex_Vectors.Vector(lvi'range);
        begin
          for k in sol'range loop
            sol(k) := lvi(k).cff(0);
          end loop;
          DoblDobl_Complex_Solutions_io.put_vector(sol);
        end;
      end if;
      Run_QR_Newton(maxdeg,nbrit,p,v(i).all,verbose,vrblvl-1);
      if verbose then
        if pause then
          put("Continue ? (y/n) "); Ask_Yes_or_No(ans);
          exit when (ans /= 'y');
        end if;
      end if;
    end loop;
  end Run_QR_Newton;

  procedure Run_QR_Newton
              ( maxdeg,nbrit : in integer32;
                p : in TripDobl_CSeries_Poly_Systems.Poly_Sys;
                v : in TripDobl_Complex_Series_VecVecs.VecVec;
                verbose : in boolean := false;
                pause : in boolean := false;
                vrblvl : in integer32 := 0 ) is

    ans : character;

  begin
    if vrblvl > 0
     then put_line("-> in power_series_methods.Run_QR_Newton 19 ...");
    end if;
    for i in v'range loop
      if verbose then
        put("Running on solution "); put(i,1); put_line(" ...");
        declare
          lvi : constant TripDobl_Complex_Series_Vectors.Link_to_Vector := v(i);
          sol : TripDobl_Complex_Vectors.Vector(lvi'range);
        begin
          for k in sol'range loop
            sol(k) := lvi(k).cff(0);
          end loop;
          TripDobl_Complex_Solutions_io.put_vector(sol);
        end;
      end if;
      Run_QR_Newton(maxdeg,nbrit,p,v(i).all,verbose,vrblvl-1);
      if verbose then
        if pause then
          put("Continue ? (y/n) "); Ask_Yes_or_No(ans);
          exit when (ans /= 'y');
        end if;
      end if;
    end loop;
  end Run_QR_Newton;

  procedure Run_QR_Newton
              ( maxdeg,nbrit : in integer32;
                p : in QuadDobl_CSeries_Poly_Systems.Poly_Sys;
                v : in QuadDobl_Complex_Series_VecVecs.VecVec;
                verbose : in boolean := false;
                pause : in boolean := false;
                vrblvl : in integer32 := 0 ) is

    ans : character;

  begin
    if vrblvl > 0
     then put_line("-> in power_series_methods.Run_QR_Newton 20 ...");
    end if;
    for i in v'range loop
      if verbose then
        put("Running on solution "); put(i,1); put_line(" ...");
        declare
          lvi : constant QuadDobl_Complex_Series_Vectors.Link_to_Vector := v(i);
          sol : QuadDobl_Complex_Vectors.Vector(lvi'range);
        begin
          for k in sol'range loop
            sol(k) := lvi(k).cff(0);
          end loop;
          QuadDobl_Complex_Solutions_io.put_vector(sol);
        end;
      end if;
      Run_QR_Newton(maxdeg,nbrit,p,v(i).all,verbose,vrblvl-1);
      if verbose then
        if pause then
          put("Continue ? (y/n) "); Ask_Yes_or_No(ans);
          exit when (ans /= 'y');
        end if;
      end if;
    end loop;
  end Run_QR_Newton;

  procedure Run_QR_Newton
              ( file : in file_type; maxdeg,nbrit : in integer32;
                p : in Standard_CSeries_Poly_Systems.Poly_Sys;
                v : in Standard_Complex_Series_VecVecs.VecVec;
                verbose : in boolean := false;
                vrblvl : in integer32 := 0 ) is
  begin
    if vrblvl > 0
     then put_line("-> in power_series_methods.Run_QR_Newton 21 ...");
    end if;
    for i in v'range loop
      if verbose then
        put(file,"Running on solution ");
        put(file,i,1); put_line(file," ...");
        declare
          lvi : constant Standard_Complex_Series_Vectors.Link_to_Vector := v(i);
          sol : Standard_Complex_Vectors.Vector(lvi'range);
        begin
          for k in sol'range loop
            sol(k) := lvi(k).cff(0);
          end loop;
          Standard_Complex_Solutions_io.put_vector(file,sol);
        end;
      end if;
      Run_QR_Newton(file,maxdeg,nbrit,p,v(i).all,verbose,vrblvl-1);
    end loop;
  end Run_QR_Newton;

  procedure Run_QR_Newton
              ( file : in file_type; maxdeg,nbrit : in integer32;
                p : in DoblDobl_CSeries_Poly_Systems.Poly_Sys;
                v : in DoblDobl_Complex_Series_VecVecs.VecVec;
                verbose : in boolean := false;
                vrblvl : in integer32 := 0 ) is
  begin
    if vrblvl > 0
     then put_line("-> in power_series_methods.Run_QR_Newton 22 ...");
    end if;
    for i in v'range loop
      if verbose then
        put(file,"Running on solution ");
        put(file,i,1); put_line(file," ...");
        declare
          lvi : constant DoblDobl_Complex_Series_Vectors.Link_to_Vector := v(i);
          sol : DoblDobl_Complex_Vectors.Vector(lvi'range);
        begin
          for k in sol'range loop
            sol(k) := lvi(k).cff(0);
          end loop;
          DoblDobl_Complex_Solutions_io.put_vector(file,sol);
        end;
      end if;
      Run_QR_Newton(file,maxdeg,nbrit,p,v(i).all,verbose,vrblvl-1);
    end loop;
  end Run_QR_Newton;

  procedure Run_QR_Newton
              ( file : in file_type; maxdeg,nbrit : in integer32;
                p : in TripDobl_CSeries_Poly_Systems.Poly_Sys;
                v : in TripDobl_Complex_Series_VecVecs.VecVec;
                verbose : in boolean := false;
                vrblvl : in integer32 := 0 ) is
  begin
    if vrblvl > 0
     then put_line("-> in power_series_methods.Run_QR_Newton 23 ...");
    end if;
    for i in v'range loop
      if verbose then
        put(file,"Running on solution ");
        put(file,i,1); put_line(file," ...");
        declare
          lvi : constant TripDobl_Complex_Series_Vectors.Link_to_Vector := v(i);
          sol : TripDobl_Complex_Vectors.Vector(lvi'range);
        begin
          for k in sol'range loop
            sol(k) := lvi(k).cff(0);
          end loop;
          TripDobl_Complex_Solutions_io.put_vector(file,sol);
        end;
      end if;
      Run_QR_Newton(file,maxdeg,nbrit,p,v(i).all,verbose,vrblvl-1);
    end loop;
  end Run_QR_Newton;

  procedure Run_QR_Newton
              ( file : in file_type; maxdeg,nbrit : in integer32;
                p : in QuadDobl_CSeries_Poly_Systems.Poly_Sys;
                v : in QuadDobl_Complex_Series_VecVecs.VecVec;
                verbose : in boolean := false;
                vrblvl : in integer32 := 0 ) is
  begin
    if vrblvl > 0
     then put_line("-> in power_series_methods.Run_QR_Newton 24 ...");
    end if;
    for i in v'range loop
      if verbose then
        put(file,"Running on solution ");
        put(file,i,1); put_line(file," ...");
        declare
          lvi : constant QuadDobl_Complex_Series_Vectors.Link_to_Vector := v(i);
          sol : QuadDobl_Complex_Vectors.Vector(lvi'range);
        begin
          for k in sol'range loop
            sol(k) := lvi(k).cff(0);
          end loop;
          QuadDobl_Complex_Solutions_io.put_vector(file,sol);
        end;
      end if;
      Run_QR_Newton(file,maxdeg,nbrit,p,v(i).all,verbose,vrblvl-1);
    end loop;
  end Run_QR_Newton;

  procedure Run_SVD_Newton
              ( maxdeg,nbrit : in integer32;
                p : in Standard_CSeries_Poly_Systems.Poly_Sys;
                v : in Standard_Complex_Series_VecVecs.VecVec;
                verbose : in boolean := false;
                pause : in boolean := false;
                vrblvl : in integer32 := 0 ) is

    ans : character;
    rcond : double_float;

  begin
    if vrblvl > 0
     then put_line("-> in power_series_methods.Run_SVD_Newton 17 ...");
    end if;
    for i in v'range loop
      if verbose then
        put("Running on solution "); put(i,1); put_line(" ...");
        declare
          lvi : constant Standard_Complex_Series_Vectors.Link_to_Vector := v(i);
          sol : Standard_Complex_Vectors.Vector(lvi'range);
        begin
          for k in sol'range loop
            sol(k) := lvi(k).cff(0);
          end loop;
          Standard_Complex_Solutions_io.put_vector(sol);
        end;
      end if;
      Run_SVD_Newton(maxdeg,nbrit,p,v(i).all,rcond,verbose,vrblvl-1);
      if verbose then
        if pause then
          put("Continue ? (y/n) "); Ask_Yes_or_No(ans);
          exit when (ans /= 'y');
        end if;
      end if;
    end loop;
  end Run_SVD_Newton;

  procedure Run_SVD_Newton
              ( maxdeg,nbrit : in integer32;
                p : in DoblDobl_CSeries_Poly_Systems.Poly_Sys;
                v : in DoblDobl_Complex_Series_VecVecs.VecVec;
                verbose : in boolean := false;
                pause : in boolean := false;
                vrblvl : in integer32 := 0 ) is

    ans : character;
    rcond : double_double;

  begin
    if vrblvl > 0
     then put_line("-> in power_series_methods.Run_SVD_Newton 18 ...");
    end if;
    for i in v'range loop
      if verbose then
        put("Running on solution "); put(i,1); put_line(" ...");
        declare
          lvi : constant DoblDobl_Complex_Series_Vectors.Link_to_Vector := v(i);
          sol : DoblDobl_Complex_Vectors.Vector(lvi'range);
        begin
          for k in sol'range loop
            sol(k) := lvi(k).cff(0);
          end loop;
          DoblDobl_Complex_Solutions_io.put_vector(sol);
        end;
      end if;
      Run_SVD_Newton(maxdeg,nbrit,p,v(i).all,rcond,verbose,vrblvl-1);
      if verbose then
        if pause then
          put("Continue ? (y/n) "); Ask_Yes_or_No(ans);
          exit when (ans /= 'y');
        end if;
      end if;
    end loop;
  end Run_SVD_Newton;

  procedure Run_SVD_Newton
              ( maxdeg,nbrit : in integer32;
                p : in TripDobl_CSeries_Poly_Systems.Poly_Sys;
                v : in TripDobl_Complex_Series_VecVecs.VecVec;
                verbose : in boolean := false;
                pause : in boolean := false;
                vrblvl : in integer32 := 0 ) is

    ans : character;
    rcond : triple_double;

  begin
    if vrblvl > 0
     then put_line("-> in power_series_methods.Run_SVD_Newton 19 ...");
    end if;
    for i in v'range loop
      if verbose then
        put("Running on solution "); put(i,1); put_line(" ...");
        declare
          lvi : constant TripDobl_Complex_Series_Vectors.Link_to_Vector := v(i);
          sol : TripDobl_Complex_Vectors.Vector(lvi'range);
        begin
          for k in sol'range loop
            sol(k) := lvi(k).cff(0);
          end loop;
          TripDobl_Complex_Solutions_io.put_vector(sol);
        end;
      end if;
      Run_SVD_Newton(maxdeg,nbrit,p,v(i).all,rcond,verbose,vrblvl-1);
      if verbose then
        if pause then
          put("Continue ? (y/n) "); Ask_Yes_or_No(ans);
          exit when (ans /= 'y');
        end if;
      end if;
    end loop;
  end Run_SVD_Newton;

  procedure Run_SVD_Newton
              ( maxdeg,nbrit : in integer32;
                p : in QuadDobl_CSeries_Poly_Systems.Poly_Sys;
                v : in QuadDobl_Complex_Series_VecVecs.VecVec;
                verbose : in boolean := false;
                pause : in boolean := false;
                vrblvl : in integer32 := 0 ) is

    ans : character;
    rcond : quad_double;

  begin
    if vrblvl > 0
     then put_line("-> in power_series_methods.Run_SVD_Newton 20 ...");
    end if;
    for i in v'range loop
      if verbose then
        put("Running on solution "); put(i,1); put_line(" ...");
        declare
          lvi : constant QuadDobl_Complex_Series_Vectors.Link_to_Vector := v(i);
          sol : QuadDobl_Complex_Vectors.Vector(lvi'range);
        begin
          for k in sol'range loop
            sol(k) := lvi(k).cff(0);
          end loop;
          QuadDobl_Complex_Solutions_io.put_vector(sol);
        end;
      end if;
      Run_SVD_Newton(maxdeg,nbrit,p,v(i).all,rcond,verbose,vrblvl-1);
      if verbose then
        if pause then
          put("Continue ? (y/n) "); Ask_Yes_or_No(ans);
          exit when (ans /= 'y');
        end if;
      end if;
    end loop;
  end Run_SVD_Newton;

  procedure Run_SVD_Newton
              ( maxdeg,nbrit : in integer32;
                p : in PentDobl_CSeries_Poly_Systems.Poly_Sys;
                v : in PentDobl_Complex_Series_VecVecs.VecVec;
                verbose : in boolean := false;
                pause : in boolean := false;
                vrblvl : in integer32 := 0 ) is

    ans : character;
    rcond : penta_double;

  begin
    if vrblvl > 0
     then put_line("-> in power_series_methods.Run_SVD_Newton 21 ...");
    end if;
    for i in v'range loop
      if verbose then
        put("Running on solution "); put(i,1); put_line(" ...");
        declare
          lvi : constant PentDobl_Complex_Series_Vectors.Link_to_Vector := v(i);
          sol : PentDobl_Complex_Vectors.Vector(lvi'range);
        begin
          for k in sol'range loop
            sol(k) := lvi(k).cff(0);
          end loop;
          PentDobl_Complex_Solutions_io.put_vector(sol);
        end;
      end if;
      Run_SVD_Newton(maxdeg,nbrit,p,v(i).all,rcond,verbose,vrblvl-1);
      if verbose then
        if pause then
          put("Continue ? (y/n) "); Ask_Yes_or_No(ans);
          exit when (ans /= 'y');
        end if;
      end if;
    end loop;
  end Run_SVD_Newton;

  procedure Run_SVD_Newton
              ( maxdeg,nbrit : in integer32;
                p : in OctoDobl_CSeries_Poly_Systems.Poly_Sys;
                v : in OctoDobl_Complex_Series_VecVecs.VecVec;
                verbose : in boolean := false;
                pause : in boolean := false;
                vrblvl : in integer32 := 0 ) is

    ans : character;
    rcond : octo_double;

  begin
    if vrblvl > 0
     then put_line("-> in power_series_methods.Run_SVD_Newton 22 ...");
    end if;
    for i in v'range loop
      if verbose then
        put("Running on solution "); put(i,1); put_line(" ...");
        declare
          lvi : constant OctoDobl_Complex_Series_Vectors.Link_to_Vector := v(i);
          sol : OctoDobl_Complex_Vectors.Vector(lvi'range);
        begin
          for k in sol'range loop
            sol(k) := lvi(k).cff(0);
          end loop;
          OctoDobl_Complex_Solutions_io.put_vector(sol);
        end;
      end if;
      Run_SVD_Newton(maxdeg,nbrit,p,v(i).all,rcond,verbose,vrblvl-1);
      if verbose then
        if pause then
          put("Continue ? (y/n) "); Ask_Yes_or_No(ans);
          exit when (ans /= 'y');
        end if;
      end if;
    end loop;
  end Run_SVD_Newton;

  procedure Run_SVD_Newton
              ( maxdeg,nbrit : in integer32;
                p : in DecaDobl_CSeries_Poly_Systems.Poly_Sys;
                v : in DecaDobl_Complex_Series_VecVecs.VecVec;
                verbose : in boolean := false;
                pause : in boolean := false;
                vrblvl : in integer32 := 0 ) is

    ans : character;
    rcond : deca_double;

  begin
    if vrblvl > 0
     then put_line("-> in power_series_methods.Run_SVD_Newton 23 ...");
    end if;
    for i in v'range loop
      if verbose then
        put("Running on solution "); put(i,1); put_line(" ...");
        declare
          lvi : constant DecaDobl_Complex_Series_Vectors.Link_to_Vector := v(i);
          sol : DecaDobl_Complex_Vectors.Vector(lvi'range);
        begin
          for k in sol'range loop
            sol(k) := lvi(k).cff(0);
          end loop;
          DecaDobl_Complex_Solutions_io.put_vector(sol);
        end;
      end if;
      Run_SVD_Newton(maxdeg,nbrit,p,v(i).all,rcond,verbose,vrblvl-1);
      if verbose then
        if pause then
          put("Continue ? (y/n) "); Ask_Yes_or_No(ans);
          exit when (ans /= 'y');
        end if;
      end if;
    end loop;
  end Run_SVD_Newton;

  procedure Run_SVD_Newton
              ( file : in file_type; maxdeg,nbrit : in integer32;
                p : in Standard_CSeries_Poly_Systems.Poly_Sys;
                v : in Standard_Complex_Series_VecVecs.VecVec;
                verbose : in boolean := false;
                vrblvl : in integer32 := 0 ) is

    rcond : double_float;

  begin
    if vrblvl > 0
     then put_line("-> in power_series_methods.Run_SVD_Newton 24 ...");
    end if;
    for i in v'range loop
      if verbose then
        put(file,"Running on solution ");
        put(file,i,1); put_line(file," ...");
        declare
          lvi : constant Standard_Complex_Series_Vectors.Link_to_Vector := v(i);
          sol : Standard_Complex_Vectors.Vector(lvi'range);
        begin
          for k in sol'range loop
            sol(k) := lvi(k).cff(0);
          end loop;
          Standard_Complex_Solutions_io.put_vector(file,sol);
        end;
      end if;
      Run_SVD_Newton(file,maxdeg,nbrit,p,v(i).all,rcond,verbose,vrblvl-1);
    end loop;
  end Run_SVD_Newton;

  procedure Run_SVD_Newton
              ( file : in file_type; maxdeg,nbrit : in integer32;
                p : in DoblDobl_CSeries_Poly_Systems.Poly_Sys;
                v : in DoblDobl_Complex_Series_VecVecs.VecVec;
                verbose : in boolean := false;
                vrblvl : in integer32 := 0 ) is

    rcond : double_double;

  begin
    if vrblvl > 0
     then put_line("-> in power_series_methods.Run_SVD_Newton 25 ...");
    end if;
    for i in v'range loop
      if verbose then
        put(file,"Running on solution ");
        put(file,i,1); put_line(file," ...");
        declare
          lvi : constant DoblDobl_Complex_Series_Vectors.Link_to_Vector := v(i);
          sol : DoblDobl_Complex_Vectors.Vector(lvi'range);
        begin
          for k in sol'range loop
            sol(k) := lvi(k).cff(0);
          end loop;
          DoblDobl_Complex_Solutions_io.put_vector(file,sol);
        end;
      end if;
      Run_SVD_Newton(file,maxdeg,nbrit,p,v(i).all,rcond,verbose,vrblvl-1);
    end loop;
  end Run_SVD_Newton;

  procedure Run_SVD_Newton
              ( file : in file_type; maxdeg,nbrit : in integer32;
                p : in TripDobl_CSeries_Poly_Systems.Poly_Sys;
                v : in TripDobl_Complex_Series_VecVecs.VecVec;
                verbose : in boolean := false;
                vrblvl : in integer32 := 0 ) is

    rcond : triple_double;

  begin
    if vrblvl > 0
     then put_line("-> in power_series_methods.Run_SVD_Newton 26 ...");
    end if;
    for i in v'range loop
      if verbose then
        put(file,"Running on solution ");
        put(file,i,1); put_line(file," ...");
        declare
          lvi : constant TripDobl_Complex_Series_Vectors.Link_to_Vector := v(i);
          sol : TripDobl_Complex_Vectors.Vector(lvi'range);
        begin
          for k in sol'range loop
            sol(k) := lvi(k).cff(0);
          end loop;
          TripDobl_Complex_Solutions_io.put_vector(file,sol);
        end;
      end if;
      Run_SVD_Newton(file,maxdeg,nbrit,p,v(i).all,rcond,verbose,vrblvl-1);
    end loop;
  end Run_SVD_Newton;

  procedure Run_SVD_Newton
              ( file : in file_type; maxdeg,nbrit : in integer32;
                p : in QuadDobl_CSeries_Poly_Systems.Poly_Sys;
                v : in QuadDobl_Complex_Series_VecVecs.VecVec;
                verbose : in boolean := false;
                vrblvl : in integer32 := 0 ) is

    rcond : quad_double;

  begin
    if vrblvl > 0
     then put_line("-> in power_series_methods.Run_SVD_Newton 27 ...");
    end if;
    for i in v'range loop
      if verbose then
        put(file,"Running on solution ");
        put(file,i,1); put_line(file," ...");
        declare
          lvi : constant QuadDobl_Complex_Series_Vectors.Link_to_Vector := v(i);
          sol : QuadDobl_Complex_Vectors.Vector(lvi'range);
        begin
          for k in sol'range loop
            sol(k) := lvi(k).cff(0);
          end loop;
          QuadDobl_Complex_Solutions_io.put_vector(file,sol);
        end;
      end if;
      Run_SVD_Newton(file,maxdeg,nbrit,p,v(i).all,rcond,verbose,vrblvl-1);
    end loop;
  end Run_SVD_Newton;

  procedure Run_SVD_Newton
              ( file : in file_type; maxdeg,nbrit : in integer32;
                p : in PentDobl_CSeries_Poly_Systems.Poly_Sys;
                v : in PentDobl_Complex_Series_VecVecs.VecVec;
                verbose : in boolean := false;
                vrblvl : in integer32 := 0 ) is

    rcond : penta_double;

  begin
    if vrblvl > 0
     then put_line("-> in power_series_methods.Run_SVD_Newton 28 ...");
    end if;
    for i in v'range loop
      if verbose then
        put(file,"Running on solution ");
        put(file,i,1); put_line(file," ...");
        declare
          lvi : constant PentDobl_Complex_Series_Vectors.Link_to_Vector := v(i);
          sol : PentDobl_Complex_Vectors.Vector(lvi'range);
        begin
          for k in sol'range loop
            sol(k) := lvi(k).cff(0);
          end loop;
          PentDobl_Complex_Solutions_io.put_vector(file,sol);
        end;
      end if;
      Run_SVD_Newton(file,maxdeg,nbrit,p,v(i).all,rcond,verbose,vrblvl-1);
    end loop;
  end Run_SVD_Newton;

  procedure Run_SVD_Newton
              ( file : in file_type; maxdeg,nbrit : in integer32;
                p : in OctoDobl_CSeries_Poly_Systems.Poly_Sys;
                v : in OctoDobl_Complex_Series_VecVecs.VecVec;
                verbose : in boolean := false;
                vrblvl : in integer32 := 0 ) is

    rcond : octo_double;

  begin
    if vrblvl > 0
     then put_line("-> in power_series_methods.Run_SVD_Newton 29 ...");
    end if;
    for i in v'range loop
      if verbose then
        put(file,"Running on solution ");
        put(file,i,1); put_line(file," ...");
        declare
          lvi : constant OctoDobl_Complex_Series_Vectors.Link_to_Vector := v(i);
          sol : OctoDobl_Complex_Vectors.Vector(lvi'range);
        begin
          for k in sol'range loop
            sol(k) := lvi(k).cff(0);
          end loop;
          OctoDobl_Complex_Solutions_io.put_vector(file,sol);
        end;
      end if;
      Run_SVD_Newton(file,maxdeg,nbrit,p,v(i).all,rcond,verbose,vrblvl-1);
    end loop;
  end Run_SVD_Newton;

  procedure Run_SVD_Newton
              ( file : in file_type; maxdeg,nbrit : in integer32;
                p : in DecaDobl_CSeries_Poly_Systems.Poly_Sys;
                v : in DecaDobl_Complex_Series_VecVecs.VecVec;
                verbose : in boolean := false;
                vrblvl : in integer32 := 0 ) is

    rcond : deca_double;

  begin
    if vrblvl > 0
     then put_line("-> in power_series_methods.Run_SVD_Newton 30 ...");
    end if;
    for i in v'range loop
      if verbose then
        put(file,"Running on solution ");
        put(file,i,1); put_line(file," ...");
        declare
          lvi : constant DecaDobl_Complex_Series_Vectors.Link_to_Vector := v(i);
          sol : DecaDobl_Complex_Vectors.Vector(lvi'range);
        begin
          for k in sol'range loop
            sol(k) := lvi(k).cff(0);
          end loop;
          DecaDobl_Complex_Solutions_io.put_vector(file,sol);
        end;
      end if;
      Run_SVD_Newton(file,maxdeg,nbrit,p,v(i).all,rcond,verbose,vrblvl-1);
    end loop;
  end Run_SVD_Newton;

  procedure Run_Echelon_Newton
              ( maxdeg,nbrit : in integer32;
                p : in Standard_CSeries_Poly_Systems.Poly_Sys;
                v : in Standard_Complex_Series_VecVecs.VecVec;
                verbose : in boolean := false;
                pause : in boolean := false;
                vrblvl : in integer32 := 0 ) is

    ans : character;
    det : Standard_Complex_Numbers.Complex_Number;

  begin
    if vrblvl > 0
     then put_line("-> in power_series_methods.Run_Echelon_Newton 9 ...");
    end if;
    for i in v'range loop
      if verbose then
        put("Running on solution "); put(i,1); put_line(" ...");
        declare
          lvi : constant Standard_Complex_Series_Vectors.Link_to_Vector := v(i);
          sol : Standard_Complex_Vectors.Vector(lvi'range);
        begin
          for k in sol'range loop
            sol(k) := lvi(k).cff(0);
          end loop;
          Standard_Complex_Solutions_io.put_vector(sol);
        end;
      end if;
      Run_Echelon_Newton(maxdeg,nbrit,p,v(i).all,det,verbose,vrblvl-1);
      if verbose then
        if pause then
          put("Continue ? (y/n) "); Ask_Yes_or_No(ans);
          exit when (ans /= 'y');
        end if;
      end if;
    end loop;
  end Run_Echelon_Newton;

  procedure Run_Echelon_Newton
              ( maxdeg,nbrit : in integer32;
                p : in DoblDobl_CSeries_Poly_Systems.Poly_Sys;
                v : in DoblDobl_Complex_Series_VecVecs.VecVec;
                verbose : in boolean := false;
                pause : in boolean := false;
                vrblvl : in integer32 := 0 ) is

    ans : character;
    det : DoblDobl_Complex_Numbers.Complex_Number;

  begin
    if vrblvl > 0
     then put_line("-> in power_series_methods.Run_Echelon_Newton 10 ...");
    end if;
    for i in v'range loop
      if verbose then
        put("Running on solution "); put(i,1); put_line(" ...");
        declare
          lvi : constant DoblDobl_Complex_Series_Vectors.Link_to_Vector := v(i);
          sol : DoblDobl_Complex_Vectors.Vector(lvi'range);
        begin
          for k in sol'range loop
            sol(k) := lvi(k).cff(0);
          end loop;
          DoblDobl_Complex_Solutions_io.put_vector(sol);
        end;
      end if;
      Run_Echelon_Newton(maxdeg,nbrit,p,v(i).all,det,verbose,vrblvl-1);
      if verbose then
        if pause then
          put("Continue ? (y/n) "); Ask_Yes_or_No(ans);
          exit when (ans /= 'y');
        end if;
      end if;
    end loop;
  end Run_Echelon_Newton;

  procedure Run_Echelon_Newton
              ( maxdeg,nbrit : in integer32;
                p : in TripDobl_CSeries_Poly_Systems.Poly_Sys;
                v : in TripDobl_Complex_Series_VecVecs.VecVec;
                verbose : in boolean := false;
                pause : in boolean := false;
                vrblvl : in integer32 := 0 ) is

    ans : character;
    det : TripDobl_Complex_Numbers.Complex_Number;

  begin
    if vrblvl > 0
     then put_line("-> in power_series_methods.Run_Echelon_Newton 11 ...");
    end if;
    for i in v'range loop
      if verbose then
        put("Running on solution "); put(i,1); put_line(" ...");
        declare
          lvi : constant TripDobl_Complex_Series_Vectors.Link_to_Vector := v(i);
          sol : TripDobl_Complex_Vectors.Vector(lvi'range);
        begin
          for k in sol'range loop
            sol(k) := lvi(k).cff(0);
          end loop;
          TripDobl_Complex_Solutions_io.put_vector(sol);
        end;
      end if;
      Run_Echelon_Newton(maxdeg,nbrit,p,v(i).all,det,verbose,vrblvl-1);
      if verbose then
        if pause then
          put("Continue ? (y/n) "); Ask_Yes_or_No(ans);
          exit when (ans /= 'y');
        end if;
      end if;
    end loop;
  end Run_Echelon_Newton;

  procedure Run_Echelon_Newton
              ( maxdeg,nbrit : in integer32;
                p : in QuadDobl_CSeries_Poly_Systems.Poly_Sys;
                v : in QuadDobl_Complex_Series_VecVecs.VecVec;
                verbose : in boolean := false;
                pause : in boolean := false;
                vrblvl : in integer32 := 0 ) is

    ans : character;
    det : QuadDobl_Complex_Numbers.Complex_Number;

  begin
    if vrblvl > 0
     then put_line("-> in power_series_methods.Run_Echelon_Newton 12 ...");
    end if;
    for i in v'range loop
      if verbose then
        put("Running on solution "); put(i,1); put_line(" ...");
        declare
          lvi : constant QuadDobl_Complex_Series_Vectors.Link_to_Vector := v(i);
          sol : QuadDobl_Complex_Vectors.Vector(lvi'range);
        begin
          for k in sol'range loop
            sol(k) := lvi(k).cff(0);
          end loop;
          QuadDobl_Complex_Solutions_io.put_vector(sol);
        end;
      end if;
      Run_Echelon_Newton(maxdeg,nbrit,p,v(i).all,det,verbose,vrblvl-1);
      if verbose then
        if pause then
          put("Continue ? (y/n) "); Ask_Yes_or_No(ans);
          exit when (ans /= 'y');
        end if;
      end if;
    end loop;
  end Run_Echelon_Newton;

  procedure Run_Echelon_Newton
              ( file : in file_type; maxdeg,nbrit : in integer32;
                p : in Standard_CSeries_Poly_Systems.Poly_Sys;
                v : in Standard_Complex_Series_VecVecs.VecVec;
                verbose : in boolean := false;
                vrblvl : in integer32 := 0 ) is

    det : Standard_Complex_Numbers.Complex_Number;

  begin
    if vrblvl > 0
     then put_line("-> in power_series_methods.Run_Echelon_Newton 13 ...");
    end if;
    for i in v'range loop
      if verbose then
        put(file,"Running on solution ");
        put(file,i,1); put_line(file," ...");
        declare
          lvi : constant Standard_Complex_Series_Vectors.Link_to_Vector := v(i);
          sol : Standard_Complex_Vectors.Vector(lvi'range);
        begin
          for k in sol'range loop
            sol(k) := lvi(k).cff(0);
          end loop;
          Standard_Complex_Solutions_io.put_vector(file,sol);
        end;
      end if;
      Run_Echelon_Newton(file,maxdeg,nbrit,p,v(i).all,det,verbose,vrblvl-1);
    end loop;
  end Run_Echelon_Newton;

  procedure Run_Echelon_Newton
              ( file : in file_type; maxdeg,nbrit : in integer32;
                p : in DoblDobl_CSeries_Poly_Systems.Poly_Sys;
                v : in DoblDobl_Complex_Series_VecVecs.VecVec;
                verbose : in boolean := false;
                vrblvl : in integer32 := 0 ) is

    det : DoblDobl_Complex_Numbers.Complex_Number;

  begin
    if vrblvl > 0
     then put_line("-> in power_series_methods.Run_Echelon_Newton 14 ...");
    end if;
    for i in v'range loop
      if verbose then
        put(file,"Running on solution ");
        put(file,i,1); put_line(file," ...");
        declare
          lvi : constant DoblDobl_Complex_Series_Vectors.Link_to_Vector := v(i);
          sol : DoblDobl_Complex_Vectors.Vector(lvi'range);
        begin
          for k in sol'range loop
            sol(k) := lvi(k).cff(0);
          end loop;
          DoblDobl_Complex_Solutions_io.put_vector(file,sol);
        end;
      end if;
      Run_Echelon_Newton(file,maxdeg,nbrit,p,v(i).all,det,verbose,vrblvl-1);
    end loop;
  end Run_Echelon_Newton;

  procedure Run_Echelon_Newton
              ( file : in file_type; maxdeg,nbrit : in integer32;
                p : in TripDobl_CSeries_Poly_Systems.Poly_Sys;
                v : in TripDobl_Complex_Series_VecVecs.VecVec;
                verbose : in boolean := false;
                vrblvl : in integer32 := 0 ) is

    det : TripDobl_Complex_Numbers.Complex_Number;

  begin
    if vrblvl > 0
     then put_line("-> in power_series_methods.Run_Echelon_Newton 15 ...");
    end if;
    for i in v'range loop
      if verbose then
        put(file,"Running on solution ");
        put(file,i,1); put_line(file," ...");
        declare
          lvi : constant TripDobl_Complex_Series_Vectors.Link_to_Vector := v(i);
          sol : TripDobl_Complex_Vectors.Vector(lvi'range);
        begin
          for k in sol'range loop
            sol(k) := lvi(k).cff(0);
          end loop;
          TripDobl_Complex_Solutions_io.put_vector(file,sol);
        end;
      end if;
      Run_Echelon_Newton(file,maxdeg,nbrit,p,v(i).all,det,verbose,vrblvl-1);
    end loop;
  end Run_Echelon_Newton;

  procedure Run_Echelon_Newton
              ( file : in file_type; maxdeg,nbrit : in integer32;
                p : in QuadDobl_CSeries_Poly_Systems.Poly_Sys;
                v : in QuadDobl_Complex_Series_VecVecs.VecVec;
                verbose : in boolean := false;
                vrblvl : in integer32 := 0 ) is

    det : QuadDobl_Complex_Numbers.Complex_Number;

  begin
    if vrblvl > 0
     then put_line("-> in power_series_methods.Run_Echelon_Newton 16 ...");
    end if;
    for i in v'range loop
      if verbose then
        put(file,"Running on solution ");
        put(file,i,1); put_line(file," ...");
        declare
          lvi : constant QuadDobl_Complex_Series_Vectors.Link_to_Vector := v(i);
          sol : QuadDobl_Complex_Vectors.Vector(lvi'range);
        begin
          for k in sol'range loop
            sol(k) := lvi(k).cff(0);
          end loop;
          QuadDobl_Complex_Solutions_io.put_vector(file,sol);
        end;
      end if;
      Run_Echelon_Newton(file,maxdeg,nbrit,p,v(i).all,det,verbose,vrblvl-1);
    end loop;
  end Run_Echelon_Newton;

end Power_Series_Methods;
