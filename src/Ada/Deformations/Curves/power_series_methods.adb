with Communications_with_User;           use Communications_with_User;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;
with Double_Double_Numbers_io;           use Double_Double_Numbers_io;
with Quad_Double_Numbers_io;             use Quad_Double_Numbers_io;
with Standard_Complex_Numbers_io;        use Standard_Complex_Numbers_io;
with DoblDobl_Complex_Numbers_io;        use DoblDobl_Complex_Numbers_io;
with QuadDobl_Complex_Numbers_io;        use QuadDobl_Complex_Numbers_io;
with Standard_Complex_Vectors;
with Standard_Complex_Solutions_io;
with DoblDobl_Complex_Vectors;
with DoblDobl_Complex_Solutions_io;
with QuadDobl_Complex_Vectors;
with QuadDobl_Complex_Solutions_io;
with Complex_Series_and_Polynomials;
with Complex_Series_and_Polynomials_io;
with Standard_Newton_Matrix_Series; -- with Standard_Newton_Series;
with DoblDobl_Newton_Matrix_Series; -- with DoblDobl_Newton_Series;
with QuadDobl_Newton_Matrix_Series; -- with QuadDobl_Newton_Series;

package body Power_Series_Methods is

  procedure Run_LU_Newton
             ( maxdeg,nbrit : in integer32;
               p : in Standard_CSeries_Poly_Systems.Poly_Sys;
               s : in out Standard_Complex_Series_Vectors.Vector;
               info : out integer32; verbose : in boolean := false ) is
  begin
    Run_LU_Newton(standard_output,maxdeg,nbrit,p,s,info,verbose);
  end Run_LU_Newton;

  procedure Run_LU_Newton
             ( file : in file_type; maxdeg,nbrit : in integer32;
               p : in Standard_CSeries_Poly_Systems.Poly_Sys;
               s : in out Standard_Complex_Series_Vectors.Vector;
               info : out integer32; verbose : in boolean := false ) is

    use Standard_Newton_Matrix_Series; -- use Standard_Newton_Series;

    order : integer32 := 1;
    tol : constant double_float := 1.0E-12;
    eva : Standard_Complex_Series_Vectors.Vector(p'range);

  begin
    if not verbose then
      LU_Newton_Steps(p,order,maxdeg,nbrit,s,info);
    else
      LU_Newton_Steps(file,p,order,maxdeg,nbrit,s,info);
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
               info : out integer32; verbose : in boolean := false ) is
  begin
    Run_LU_Newton(standard_output,maxdeg,nbrit,p,s,info,verbose);
  end Run_LU_Newton;

  procedure Run_LU_Newton
             ( file : in file_type; maxdeg,nbrit : in integer32;
               p : in DoblDobl_CSeries_Poly_Systems.Poly_Sys;
               s : in out DoblDobl_Complex_Series_Vectors.Vector;
               info : out integer32; verbose : in boolean := false ) is

    use DoblDobl_Newton_Matrix_Series; -- use DoblDobl_Newton_Series;

    order : integer32 := 1;
    tol : constant double_float := 1.0E-24;
    eva : DoblDobl_Complex_Series_Vectors.Vector(p'range);

  begin
    if not verbose then
      LU_Newton_Steps(p,order,maxdeg,nbrit,s,info);
    else
      LU_Newton_Steps(file,p,order,maxdeg,nbrit,s,info);
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
               p : in QuadDobl_CSeries_Poly_Systems.Poly_Sys;
               s : in out QuadDobl_Complex_Series_Vectors.Vector;
               info : out integer32; verbose : in boolean := false ) is
  begin
    Run_LU_Newton(standard_output,maxdeg,nbrit,p,s,info,verbose);
  end Run_LU_Newton;

  procedure Run_LU_Newton
             ( file : in file_type; maxdeg,nbrit : in integer32;
               p : in QuadDobl_CSeries_Poly_Systems.Poly_Sys;
               s : in out QuadDobl_Complex_Series_Vectors.Vector;
               info : out integer32; verbose : in boolean := false ) is

    use QuadDobl_Newton_Matrix_Series; -- use QuadDobl_Newton_Series;

    order : integer32 := 1;
    tol : constant double_float := 1.0E-32;
    eva : QuadDobl_Complex_Series_Vectors.Vector(p'range);

  begin
    if not verbose then
      LU_Newton_Steps(p,order,maxdeg,nbrit,s,info);
    else
      LU_Newton_Steps(file,p,order,maxdeg,nbrit,s,info);
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
               rcond : out double_float; verbose : in boolean := false ) is
  begin
    Run_LU_Newton(standard_output,maxdeg,nbrit,p,s,rcond,verbose);
  end Run_LU_Newton;

  procedure Run_LU_Newton
             ( file : in file_type; maxdeg,nbrit : in integer32;
               p : in Standard_CSeries_Poly_Systems.Poly_Sys;
               s : in out Standard_Complex_Series_Vectors.Vector;
               rcond : out double_float; verbose : in boolean := false ) is

    use Standard_Newton_Matrix_Series; -- use Standard_Newton_Series;

    order : integer32 := 1;
    tol : constant double_float := 1.0E-12;
    eva : Standard_Complex_Series_Vectors.Vector(p'range);

  begin
    if not verbose then
      LU_Newton_Steps(p,order,maxdeg,nbrit,s,rcond);
    else
      LU_Newton_Steps(file,p,order,maxdeg,nbrit,s,rcond);
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
               rcond : out double_double; verbose : in boolean := false ) is
  begin
    Run_LU_Newton(standard_output,maxdeg,nbrit,p,s,rcond,verbose);
  end Run_LU_Newton;

  procedure Run_LU_Newton
             ( file : in file_type; maxdeg,nbrit : in integer32;
               p : in DoblDobl_CSeries_Poly_Systems.Poly_Sys;
               s : in out DoblDobl_Complex_Series_Vectors.Vector;
               rcond : out double_double; verbose : in boolean := false ) is

    use DoblDobl_Newton_Matrix_Series; -- use DoblDobl_Newton_Series;

    order : integer32 := 1;
    tol : constant double_float := 1.0E-24;
    eva : DoblDobl_Complex_Series_Vectors.Vector(p'range);
    one : constant double_double := create(1.0);

  begin
    if not verbose then
      LU_Newton_Steps(p,order,maxdeg,nbrit,s,rcond);
    else
      LU_Newton_Steps(file,p,order,maxdeg,nbrit,s,rcond);
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
               p : in QuadDobl_CSeries_Poly_Systems.Poly_Sys;
               s : in out QuadDobl_Complex_Series_Vectors.Vector;
               rcond : out quad_double; verbose : in boolean := false ) is
  begin
    Run_LU_Newton(standard_output,maxdeg,nbrit,p,s,rcond,verbose);
  end Run_LU_Newton;

  procedure Run_LU_Newton
             ( file : in file_type; maxdeg,nbrit : in integer32;
               p : in QuadDobl_CSeries_Poly_Systems.Poly_Sys;
               s : in out QuadDobl_Complex_Series_Vectors.Vector;
               rcond : out quad_double; verbose : in boolean := false ) is

    use QuadDobl_Newton_Matrix_Series; -- use QuadDobl_Newton_Series;

    order : integer32 := 1;
    tol : constant double_float := 1.0E-32;
    eva : QuadDobl_Complex_Series_Vectors.Vector(p'range);
    one : constant quad_double := create(1.0);

  begin
    if not verbose then
      LU_Newton_Steps(p,order,maxdeg,nbrit,s,rcond);
    else
      LU_Newton_Steps(file,p,order,maxdeg,nbrit,s,rcond);
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
               info : out integer32; verbose : in boolean := false ) is
  begin
    Run_LU_Newton(standard_output,maxdeg,nbrit,f,c,ejm,mlt,s,info,verbose);
  end Run_LU_Newton;

  procedure Run_LU_Newton
             ( file : in file_type; maxdeg,nbrit : in integer32;
               f : in Standard_CSeries_Poly_SysFun.Eval_Coeff_Poly_Sys;
               c : in Standard_Complex_Series_VecVecs.VecVec;
               ejm : in Standard_CSeries_Jaco_Matrices.Eval_Coeff_Jaco_Mat;
               mlt : in Standard_CSeries_Jaco_Matrices.Mult_Factors;
               s : in out Standard_Complex_Series_Vectors.Vector;
               info : out integer32; verbose : in boolean := false ) is

    use Standard_Newton_Matrix_Series; -- use Standard_Newton_Series;

    order : integer32 := 1;
    tol : constant double_float := 1.0E-12;
    eva : Standard_Complex_Series_Vectors.Vector(f'range);

  begin
    if not verbose then
      LU_Newton_Steps(f,c,ejm,mlt,order,maxdeg,nbrit,s,info);
    else
      LU_Newton_Steps(file,f,c,ejm,mlt,order,maxdeg,nbrit,s,info);
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
               info : out integer32; verbose : in boolean := false ) is
  begin
    Run_LU_Newton(standard_output,maxdeg,nbrit,f,c,ejm,mlt,s,info,verbose);
  end Run_LU_Newton;

  procedure Run_LU_Newton
             ( file : in file_type; maxdeg,nbrit : in integer32;
               f : in DoblDobl_CSeries_Poly_SysFun.Eval_Coeff_Poly_Sys;
               c : in DoblDobl_Complex_Series_VecVecs.VecVec;
               ejm : in DoblDobl_CSeries_Jaco_Matrices.Eval_Coeff_Jaco_Mat;
               mlt : in DoblDobl_CSeries_Jaco_Matrices.Mult_Factors;
               s : in out DoblDobl_Complex_Series_Vectors.Vector;
               info : out integer32; verbose : in boolean := false ) is

    use DoblDobl_Newton_Matrix_Series; -- use DoblDobl_Newton_Series;

    order : integer32 := 1;
    tol : constant double_float := 1.0E-12;
    eva : DoblDobl_Complex_Series_Vectors.Vector(f'range);

  begin
    if not verbose then
      LU_Newton_Steps(f,c,ejm,mlt,order,maxdeg,nbrit,s,info);
    else
      LU_Newton_Steps(file,f,c,ejm,mlt,order,maxdeg,nbrit,s,info);
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
               f : in QuadDobl_CSeries_Poly_SysFun.Eval_Coeff_Poly_Sys;
               c : in QuadDobl_Complex_Series_VecVecs.VecVec;
               ejm : in QuadDobl_CSeries_Jaco_Matrices.Eval_Coeff_Jaco_Mat;
               mlt : in QuadDobl_CSeries_Jaco_Matrices.Mult_Factors;
               s : in out QuadDobl_Complex_Series_Vectors.Vector;
               info : out integer32; verbose : in boolean := false ) is
  begin
    Run_LU_Newton(standard_output,maxdeg,nbrit,f,c,ejm,mlt,s,info,verbose);
  end Run_LU_Newton;

  procedure Run_LU_Newton
             ( file : in file_type; maxdeg,nbrit : in integer32;
               f : in QuadDobl_CSeries_Poly_SysFun.Eval_Coeff_Poly_Sys;
               c : in QuadDobl_Complex_Series_VecVecs.VecVec;
               ejm : in QuadDobl_CSeries_Jaco_Matrices.Eval_Coeff_Jaco_Mat;
               mlt : in QuadDobl_CSeries_Jaco_Matrices.Mult_Factors;
               s : in out QuadDobl_Complex_Series_Vectors.Vector;
               info : out integer32; verbose : in boolean := false ) is

    use QuadDobl_Newton_Matrix_Series; -- use QuadDobl_Newton_Series;

    order : integer32 := 1;
    tol : constant double_float := 1.0E-12;
    eva : QuadDobl_Complex_Series_Vectors.Vector(f'range);

  begin
    if not verbose then
      LU_Newton_Steps(f,c,ejm,mlt,order,maxdeg,nbrit,s,info);
    else
      LU_Newton_Steps(file,f,c,ejm,mlt,order,maxdeg,nbrit,s,info);
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
               verbose : in boolean := false ) is
  begin
    Run_QR_Newton(standard_output,maxdeg,nbrit,p,s,verbose);
  end Run_QR_Newton;

  procedure Run_QR_Newton
             ( file : in file_type; maxdeg,nbrit : in integer32;
               p : in Standard_CSeries_Poly_Systems.Poly_Sys;
               s : in out Standard_Complex_Series_Vectors.Vector;
               verbose : in boolean := false ) is

    use Standard_Newton_Matrix_Series; -- use Standard_Newton_Series;

    order : integer32 := 1;
    info : integer32;
    tol : constant double_float := 1.0E-12;
    eva : Standard_Complex_Series_Vectors.Vector(p'range);

  begin
    if not verbose then
      QR_Newton_Steps(p,order,maxdeg,nbrit,s,info);
    else
      QR_Newton_Steps(file,p,order,maxdeg,nbrit,s,info);
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
               verbose : in boolean := false ) is
  begin
    Run_QR_Newton(standard_output,maxdeg,nbrit,p,s,verbose);
  end Run_QR_Newton;

  procedure Run_QR_Newton
             ( file : in file_type; maxdeg,nbrit : in integer32;
               p : in DoblDobl_CSeries_Poly_Systems.Poly_Sys;
               s : in out DoblDobl_Complex_Series_Vectors.Vector;
               verbose : in boolean := false ) is

    use DoblDobl_Newton_Matrix_Series; -- use DoblDobl_Newton_Series;

    order : integer32 := 1;
    info : integer32;
    tol : constant double_float := 1.0E-24;
    eva : DoblDobl_Complex_Series_Vectors.Vector(p'range);

  begin
    if not verbose then
      QR_Newton_Steps(p,order,maxdeg,nbrit,s,info);
    else
      QR_Newton_Steps(file,p,order,maxdeg,nbrit,s,info);
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
               p : in QuadDobl_CSeries_Poly_Systems.Poly_Sys;
               s : in out QuadDobl_Complex_Series_Vectors.Vector;
               verbose : in boolean := false ) is
  begin
    Run_QR_Newton(standard_output,maxdeg,nbrit,p,s,verbose);
  end Run_QR_Newton;

  procedure Run_QR_Newton
             ( file : in file_type; maxdeg,nbrit : in integer32;
               p : in QuadDobl_CSeries_Poly_Systems.Poly_Sys;
               s : in out QuadDobl_Complex_Series_Vectors.Vector;
               verbose : in boolean := false ) is

    use QuadDobl_Newton_Matrix_Series; -- use QuadDobl_Newton_Series;

    order : integer32 := 1;
    info : integer32;
    tol : constant double_float := 1.0E-32;
    eva : QuadDobl_Complex_Series_Vectors.Vector(p'range);

  begin
    if not verbose then
      QR_Newton_Steps(p,order,maxdeg,nbrit,s,info);
    else
      QR_Newton_Steps(file,p,order,maxdeg,nbrit,s,info);
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
                verbose : in boolean := false ) is
  begin
    Run_QR_Newton(standard_output,maxdeg,nbrit,f,c,ejm,mlt,s,verbose);
  end Run_QR_Newton;

  procedure Run_QR_Newton
              ( file : in file_type; maxdeg,nbrit : in integer32;
                f : in Standard_CSeries_Poly_SysFun.Eval_Coeff_Poly_Sys;
                c : in Standard_Complex_Series_VecVecs.VecVec;
                ejm : in Standard_CSeries_Jaco_Matrices.Eval_Coeff_Jaco_Mat;
                mlt : in Standard_CSeries_Jaco_Matrices.Mult_Factors;
                s : in out Standard_Complex_Series_Vectors.Vector;
                verbose : in boolean := false ) is

    use Standard_Newton_Matrix_Series; -- use Standard_Newton_Series;

    order : integer32 := 1;
    info : integer32;
    tol : constant double_float := 1.0E-12;
    eva : Standard_Complex_Series_Vectors.Vector(f'range);

  begin
    if not verbose then
      QR_Newton_Steps(f,c,ejm,mlt,order,maxdeg,nbrit,s,info);
    else
      QR_Newton_Steps(file,f,c,ejm,mlt,order,maxdeg,nbrit,s,info);
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
                verbose : in boolean := false ) is
  begin
    Run_QR_Newton(standard_output,maxdeg,nbrit,f,c,ejm,mlt,s,verbose);
  end Run_QR_Newton;

  procedure Run_QR_Newton
              ( file : in file_type; maxdeg,nbrit : in integer32;
                f : in DoblDobl_CSeries_Poly_SysFun.Eval_Coeff_Poly_Sys;
                c : in DoblDobl_Complex_Series_VecVecs.VecVec;
                ejm : in DoblDobl_CSeries_Jaco_Matrices.Eval_Coeff_Jaco_Mat;
                mlt : in DoblDobl_CSeries_Jaco_Matrices.Mult_Factors;
                s : in out DoblDobl_Complex_Series_Vectors.Vector;
                verbose : in boolean := false ) is

    use DoblDobl_Newton_Matrix_Series; -- use DoblDobl_Newton_Series;

    order : integer32 := 1;
    info : integer32;
    tol : constant double_float := 1.0E-24;
    eva : DoblDobl_Complex_Series_Vectors.Vector(f'range);

  begin
    if not verbose then
      QR_Newton_Steps(f,c,ejm,mlt,order,maxdeg,nbrit,s,info);
    else
      QR_Newton_Steps(file,f,c,ejm,mlt,order,maxdeg,nbrit,s,info);
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
                f : in QuadDobl_CSeries_Poly_SysFun.Eval_Coeff_Poly_Sys;
                c : in QuadDobl_Complex_Series_VecVecs.VecVec;
                ejm : in QuadDobl_CSeries_Jaco_Matrices.Eval_Coeff_Jaco_Mat;
                mlt : in QuadDobl_CSeries_Jaco_Matrices.Mult_Factors;
                s : in out QuadDobl_Complex_Series_Vectors.Vector;
                verbose : in boolean := false ) is
  begin
    Run_QR_Newton(standard_output,maxdeg,nbrit,f,c,ejm,mlt,s,verbose);
  end Run_QR_Newton;

  procedure Run_QR_Newton
              ( file : in file_type; maxdeg,nbrit : in integer32;
                f : in QuadDobl_CSeries_Poly_SysFun.Eval_Coeff_Poly_Sys;
                c : in QuadDobl_Complex_Series_VecVecs.VecVec;
                ejm : in QuadDobl_CSeries_Jaco_Matrices.Eval_Coeff_Jaco_Mat;
                mlt : in QuadDobl_CSeries_Jaco_Matrices.Mult_Factors;
                s : in out QuadDobl_Complex_Series_Vectors.Vector;
                verbose : in boolean := false ) is

    use QuadDobl_Newton_Matrix_Series; -- use QuadDobl_Newton_Series;

    order : integer32 := 1;
    info : integer32;
    tol : constant double_float := 1.0E-48;
    eva : QuadDobl_Complex_Series_Vectors.Vector(f'range);

  begin
    if not verbose then
      QR_Newton_Steps(f,c,ejm,mlt,order,maxdeg,nbrit,s,info);
    else
      QR_Newton_Steps(file,f,c,ejm,mlt,order,maxdeg,nbrit,s,info);
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
               rcond : out double_float; verbose : in boolean := false ) is
  begin
    Run_SVD_Newton(standard_output,maxdeg,nbrit,p,s,rcond,verbose);
  end Run_SVD_Newton;

  procedure Run_SVD_Newton
             ( file : in file_type; maxdeg,nbrit : in integer32;
               p : in Standard_CSeries_Poly_Systems.Poly_Sys;
               s : in out Standard_Complex_Series_Vectors.Vector;
               rcond : out double_float; verbose : in boolean := false ) is

    use Standard_Newton_Matrix_Series; -- use Standard_Newton_Series;

    order : integer32 := 1;
    info : integer32;
    tol : constant double_float := 1.0E-12;
    eva : Standard_Complex_Series_Vectors.Vector(p'range);

  begin
    if not verbose then
      SVD_Newton_Steps(p,order,maxdeg,nbrit,s,info,rcond);
    else
      SVD_Newton_Steps(file,p,order,maxdeg,nbrit,s,info,rcond);
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
             ( maxdeg,nbrit : in integer32;
               p : in DoblDobl_CSeries_Poly_Systems.Poly_Sys;
               s : in out DoblDobl_Complex_Series_Vectors.Vector;
               rcond : out double_double; verbose : in boolean := false ) is
  begin
    Run_SVD_Newton(standard_output,maxdeg,nbrit,p,s,rcond,verbose);
  end Run_SVD_Newton;

  procedure Run_SVD_Newton
             ( file : in file_type; maxdeg,nbrit : in integer32;
               p : in DoblDobl_CSeries_Poly_Systems.Poly_Sys;
               s : in out DoblDobl_Complex_Series_Vectors.Vector;
               rcond : out double_double; verbose : in boolean := false ) is

    use DoblDobl_Newton_Matrix_Series; -- use DoblDobl_Newton_Series;

    order : integer32 := 1;
    info : integer32;
    tol : constant double_float := 1.0E-24;
    eva : DoblDobl_Complex_Series_Vectors.Vector(p'range);
    one : constant double_double := create(1.0);

  begin
    if not verbose then
      SVD_Newton_Steps(p,order,maxdeg,nbrit,s,info,rcond);
    else
      SVD_Newton_Steps(file,p,order,maxdeg,nbrit,s,info,rcond);
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
               p : in QuadDobl_CSeries_Poly_Systems.Poly_Sys;
               s : in out QuadDobl_Complex_Series_Vectors.Vector;
               rcond : out quad_double; verbose : in boolean := false ) is
  begin
    Run_SVD_Newton(standard_output,maxdeg,nbrit,p,s,rcond,verbose);
  end Run_SVD_Newton;

  procedure Run_SVD_Newton
             ( file : in file_type; maxdeg,nbrit : in integer32;
               p : in QuadDobl_CSeries_Poly_Systems.Poly_Sys;
               s : in out QuadDobl_Complex_Series_Vectors.Vector;
               rcond : out quad_double; verbose : in boolean := false ) is

    use QuadDobl_Newton_Matrix_Series; -- use QuadDobl_Newton_Series;

    order : integer32 := 1;
    info : integer32;
    tol : constant double_float := 1.0E-32;
    eva : QuadDobl_Complex_Series_Vectors.Vector(p'range);
    one : constant quad_double := create(1.0);

  begin
    if not verbose then
      SVD_Newton_Steps(p,order,maxdeg,nbrit,s,info,rcond);
    else
      SVD_Newton_Steps(file,p,order,maxdeg,nbrit,s,info,rcond);
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

  procedure Run_Echelon_Newton
              ( maxdeg,nbrit : in integer32;
                p : in Standard_CSeries_Poly_Systems.Poly_Sys;
                s : in out Standard_Complex_Series_Vectors.Vector;
                det : out Standard_Complex_Numbers.Complex_Number;
                verbose : in boolean := false ) is
  begin
    Run_Echelon_Newton(standard_output,maxdeg,nbrit,p,s,det,verbose);
  end Run_Echelon_Newton;

  procedure Run_Echelon_Newton
              ( file : in file_type; maxdeg,nbrit : in integer32;
                p : in Standard_CSeries_Poly_Systems.Poly_Sys;
                s : in out Standard_Complex_Series_Vectors.Vector;
                det : out Standard_Complex_Numbers.Complex_Number;
                verbose : in boolean := false ) is

    use Standard_Newton_Matrix_Series;

    order : integer32 := 1;
    tol : constant double_float := 1.0E-12;
    eva : Standard_Complex_Series_Vectors.Vector(p'range);

  begin
    if not verbose then
      Echelon_Newton_Steps(p,order,maxdeg,nbrit,s,det);
    else
      Echelon_Newton_Steps(file,p,order,maxdeg,nbrit,s,det);
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
                verbose : in boolean := false ) is
  begin
    Run_Echelon_Newton(standard_output,maxdeg,nbrit,p,s,det,verbose);
  end Run_Echelon_Newton;

  procedure Run_Echelon_Newton
              ( file : in file_type; maxdeg,nbrit : in integer32;
                p : in DoblDobl_CSeries_Poly_Systems.Poly_Sys;
                s : in out DoblDobl_Complex_Series_Vectors.Vector;
                det : out DoblDobl_Complex_Numbers.Complex_Number;
                verbose : in boolean := false ) is

    use DoblDobl_Newton_Matrix_Series;

    order : integer32 := 1;
    tol : constant double_float := 1.0E-12;
    eva : DoblDobl_Complex_Series_Vectors.Vector(p'range);

  begin
    if not verbose then
      Echelon_Newton_Steps(p,order,maxdeg,nbrit,s,det);
    else
      Echelon_Newton_Steps(file,p,order,maxdeg,nbrit,s,det);
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
                p : in QuadDobl_CSeries_Poly_Systems.Poly_Sys;
                s : in out QuadDobl_Complex_Series_Vectors.Vector;
                det : out QuadDobl_Complex_Numbers.Complex_Number;
                verbose : in boolean := false ) is
  begin
    Run_Echelon_Newton(standard_output,maxdeg,nbrit,p,s,det,verbose);
  end Run_Echelon_Newton;

  procedure Run_Echelon_Newton
              ( file : in file_type; maxdeg,nbrit : in integer32;
                p : in QuadDobl_CSeries_Poly_Systems.Poly_Sys;
                s : in out QuadDobl_Complex_Series_Vectors.Vector;
                det : out QuadDobl_Complex_Numbers.Complex_Number;
                verbose : in boolean := false ) is

    use QuadDobl_Newton_Matrix_Series;

    order : integer32 := 1;
    tol : constant double_float := 1.0E-12;
    eva : QuadDobl_Complex_Series_Vectors.Vector(p'range);

  begin
    if not verbose then
      Echelon_Newton_Steps(p,order,maxdeg,nbrit,s,det);
    else
      Echelon_Newton_Steps(file,p,order,maxdeg,nbrit,s,det);
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
                pause : in boolean := false ) is

    ans : character;
    info : integer32;

  begin
    for i in v'range loop
      if verbose then
        put("Running on solution "); put(i,1); put_line(" ...");
        declare
          lvi : Standard_Complex_Series_Vectors.Link_to_Vector := v(i);
          sol : Standard_Complex_Vectors.Vector(lvi'range);
        begin
          for k in sol'range loop
            sol(k) := lvi(k).cff(0);
          end loop;
          Standard_Complex_Solutions_io.put_vector(sol);
        end;
      end if;
      Run_LU_Newton(maxdeg,nbrit,p,v(i).all,info,verbose);
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
                pause : in boolean := false ) is

    ans : character;
    info : integer32;

  begin
    for i in v'range loop
      if verbose then
        put("Running on solution "); put(i,1); put_line(" ...");
        declare
          lvi : DoblDobl_Complex_Series_Vectors.Link_to_Vector := v(i);
          sol : DoblDobl_Complex_Vectors.Vector(lvi'range);
        begin
          for k in sol'range loop
            sol(k) := lvi(k).cff(0);
          end loop;
          DoblDobl_Complex_Solutions_io.put_vector(sol);
        end;
      end if;
      Run_LU_Newton(maxdeg,nbrit,p,v(i).all,info,verbose);
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
                pause : in boolean := false ) is

    ans : character;
    info : integer32;

  begin
    for i in v'range loop
      if verbose then
        put("Running on solution "); put(i,1); put_line(" ...");
        declare
          lvi : QuadDobl_Complex_Series_Vectors.Link_to_Vector := v(i);
          sol : QuadDobl_Complex_Vectors.Vector(lvi'range);
        begin
          for k in sol'range loop
            sol(k) := lvi(k).cff(0);
          end loop;
          QuadDobl_Complex_Solutions_io.put_vector(sol);
        end;
      end if;
      Run_LU_Newton(maxdeg,nbrit,p,v(i).all,info,verbose);
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
                verbose : in boolean := false ) is

    info : integer32;

  begin
    for i in v'range loop
      if verbose then
        put(file,"Running on solution ");
        put(file,i,1); put_line(file," ...");
        declare
          lvi : Standard_Complex_Series_Vectors.Link_to_Vector := v(i);
          sol : Standard_Complex_Vectors.Vector(lvi'range);
        begin
          for k in sol'range loop
            sol(k) := lvi(k).cff(0);
          end loop;
          Standard_Complex_Solutions_io.put_vector(file,sol);
        end;
      end if;
      Run_LU_Newton(file,maxdeg,nbrit,p,v(i).all,info,verbose);
    end loop;
  end Run_LU_Newton;

  procedure Run_LU_Newton
              ( file : in file_type; maxdeg,nbrit : in integer32;
                p : in DoblDobl_CSeries_Poly_Systems.Poly_Sys;
                v : in DoblDobl_Complex_Series_VecVecs.VecVec;
                verbose : in boolean := false ) is

    info : integer32;

  begin
    for i in v'range loop
      if verbose then
        put(file,"Running on solution ");
        put(file,i,1); put_line(file," ...");
        declare
          lvi : DoblDobl_Complex_Series_Vectors.Link_to_Vector := v(i);
          sol : DoblDobl_Complex_Vectors.Vector(lvi'range);
        begin
          for k in sol'range loop
            sol(k) := lvi(k).cff(0);
          end loop;
          DoblDobl_Complex_Solutions_io.put_vector(file,sol);
        end;
      end if;
      Run_LU_Newton(file,maxdeg,nbrit,p,v(i).all,info,verbose);
    end loop;
  end Run_LU_Newton;

  procedure Run_LU_Newton
              ( file : in file_type; maxdeg,nbrit : in integer32;
                p : in QuadDobl_CSeries_Poly_Systems.Poly_Sys;
                v : in QuadDobl_Complex_Series_VecVecs.VecVec;
                verbose : in boolean := false ) is

    info : integer32;

  begin
    for i in v'range loop
      if verbose then
        put(file,"Running on solution ");
        put(file,i,1); put_line(file," ...");
        declare
          lvi : QuadDobl_Complex_Series_Vectors.Link_to_Vector := v(i);
          sol : QuadDobl_Complex_Vectors.Vector(lvi'range);
        begin
          for k in sol'range loop
            sol(k) := lvi(k).cff(0);
          end loop;
          QuadDobl_Complex_Solutions_io.put_vector(file,sol);
        end;
      end if;
      Run_LU_Newton(file,maxdeg,nbrit,p,v(i).all,info,verbose);
    end loop;
  end Run_LU_Newton;

  procedure Run_QR_Newton
              ( maxdeg,nbrit : in integer32;
                p : in Standard_CSeries_Poly_Systems.Poly_Sys;
                v : in Standard_Complex_Series_VecVecs.VecVec;
                verbose : in boolean := false;
                pause : in boolean := false ) is

    ans : character;

  begin
    for i in v'range loop
      if verbose then
        put("Running on solution "); put(i,1); put_line(" ...");
        declare
          lvi : Standard_Complex_Series_Vectors.Link_to_Vector := v(i);
          sol : Standard_Complex_Vectors.Vector(lvi'range);
        begin
          for k in sol'range loop
            sol(k) := lvi(k).cff(0);
          end loop;
          Standard_Complex_Solutions_io.put_vector(sol);
        end;
      end if;
      Run_QR_Newton(maxdeg,nbrit,p,v(i).all,verbose);
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
                pause : in boolean := false ) is

    ans : character;

  begin
    for i in v'range loop
      if verbose then
        put("Running on solution "); put(i,1); put_line(" ...");
        declare
          lvi : DoblDobl_Complex_Series_Vectors.Link_to_Vector := v(i);
          sol : DoblDobl_Complex_Vectors.Vector(lvi'range);
        begin
          for k in sol'range loop
            sol(k) := lvi(k).cff(0);
          end loop;
          DoblDobl_Complex_Solutions_io.put_vector(sol);
        end;
      end if;
      Run_QR_Newton(maxdeg,nbrit,p,v(i).all,verbose);
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
                pause : in boolean := false ) is

    ans : character;

  begin
    for i in v'range loop
      if verbose then
        put("Running on solution "); put(i,1); put_line(" ...");
        declare
          lvi : QuadDobl_Complex_Series_Vectors.Link_to_Vector := v(i);
          sol : QuadDobl_Complex_Vectors.Vector(lvi'range);
        begin
          for k in sol'range loop
            sol(k) := lvi(k).cff(0);
          end loop;
          QuadDobl_Complex_Solutions_io.put_vector(sol);
        end;
      end if;
      Run_QR_Newton(maxdeg,nbrit,p,v(i).all,verbose);
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
                verbose : in boolean := false ) is
  begin
    for i in v'range loop
      if verbose then
        put(file,"Running on solution ");
        put(file,i,1); put_line(file," ...");
        declare
          lvi : Standard_Complex_Series_Vectors.Link_to_Vector := v(i);
          sol : Standard_Complex_Vectors.Vector(lvi'range);
        begin
          for k in sol'range loop
            sol(k) := lvi(k).cff(0);
          end loop;
          Standard_Complex_Solutions_io.put_vector(file,sol);
        end;
      end if;
      Run_QR_Newton(file,maxdeg,nbrit,p,v(i).all,verbose);
    end loop;
  end Run_QR_Newton;

  procedure Run_QR_Newton
              ( file : in file_type; maxdeg,nbrit : in integer32;
                p : in DoblDobl_CSeries_Poly_Systems.Poly_Sys;
                v : in DoblDobl_Complex_Series_VecVecs.VecVec;
                verbose : in boolean := false ) is
  begin
    for i in v'range loop
      if verbose then
        put(file,"Running on solution ");
        put(file,i,1); put_line(file," ...");
        declare
          lvi : DoblDobl_Complex_Series_Vectors.Link_to_Vector := v(i);
          sol : DoblDobl_Complex_Vectors.Vector(lvi'range);
        begin
          for k in sol'range loop
            sol(k) := lvi(k).cff(0);
          end loop;
          DoblDobl_Complex_Solutions_io.put_vector(file,sol);
        end;
      end if;
      Run_QR_Newton(file,maxdeg,nbrit,p,v(i).all,verbose);
    end loop;
  end Run_QR_Newton;

  procedure Run_QR_Newton
              ( file : in file_type; maxdeg,nbrit : in integer32;
                p : in QuadDobl_CSeries_Poly_Systems.Poly_Sys;
                v : in QuadDobl_Complex_Series_VecVecs.VecVec;
                verbose : in boolean := false ) is
  begin
    for i in v'range loop
      if verbose then
        put(file,"Running on solution ");
        put(file,i,1); put_line(file," ...");
        declare
          lvi : QuadDobl_Complex_Series_Vectors.Link_to_Vector := v(i);
          sol : QuadDobl_Complex_Vectors.Vector(lvi'range);
        begin
          for k in sol'range loop
            sol(k) := lvi(k).cff(0);
          end loop;
          QuadDobl_Complex_Solutions_io.put_vector(file,sol);
        end;
      end if;
      Run_QR_Newton(file,maxdeg,nbrit,p,v(i).all,verbose);
    end loop;
  end Run_QR_Newton;

  procedure Run_SVD_Newton
              ( maxdeg,nbrit : in integer32;
                p : in Standard_CSeries_Poly_Systems.Poly_Sys;
                v : in Standard_Complex_Series_VecVecs.VecVec;
                verbose : in boolean := false;
                pause : in boolean := false ) is

    ans : character;
    rcond : double_float;

  begin
    for i in v'range loop
      if verbose then
        put("Running on solution "); put(i,1); put_line(" ...");
        declare
          lvi : Standard_Complex_Series_Vectors.Link_to_Vector := v(i);
          sol : Standard_Complex_Vectors.Vector(lvi'range);
        begin
          for k in sol'range loop
            sol(k) := lvi(k).cff(0);
          end loop;
          Standard_Complex_Solutions_io.put_vector(sol);
        end;
      end if;
      Run_SVD_Newton(maxdeg,nbrit,p,v(i).all,rcond,verbose);
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
                pause : in boolean := false ) is

    ans : character;
    rcond : double_double;

  begin
    for i in v'range loop
      if verbose then
        put("Running on solution "); put(i,1); put_line(" ...");
        declare
          lvi : DoblDobl_Complex_Series_Vectors.Link_to_Vector := v(i);
          sol : DoblDobl_Complex_Vectors.Vector(lvi'range);
        begin
          for k in sol'range loop
            sol(k) := lvi(k).cff(0);
          end loop;
          DoblDobl_Complex_Solutions_io.put_vector(sol);
        end;
      end if;
      Run_SVD_Newton(maxdeg,nbrit,p,v(i).all,rcond,verbose);
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
                pause : in boolean := false ) is

    ans : character;
    rcond : quad_double;

  begin
    for i in v'range loop
      if verbose then
        put("Running on solution "); put(i,1); put_line(" ...");
        declare
          lvi : QuadDobl_Complex_Series_Vectors.Link_to_Vector := v(i);
          sol : QuadDobl_Complex_Vectors.Vector(lvi'range);
        begin
          for k in sol'range loop
            sol(k) := lvi(k).cff(0);
          end loop;
          QuadDobl_Complex_Solutions_io.put_vector(sol);
        end;
      end if;
      Run_SVD_Newton(maxdeg,nbrit,p,v(i).all,rcond,verbose);
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
                verbose : in boolean := false ) is

    rcond : double_float;

  begin
    for i in v'range loop
      if verbose then
        put(file,"Running on solution ");
        put(file,i,1); put_line(file," ...");
        declare
          lvi : Standard_Complex_Series_Vectors.Link_to_Vector := v(i);
          sol : Standard_Complex_Vectors.Vector(lvi'range);
        begin
          for k in sol'range loop
            sol(k) := lvi(k).cff(0);
          end loop;
          Standard_Complex_Solutions_io.put_vector(file,sol);
        end;
      end if;
      Run_SVD_Newton(file,maxdeg,nbrit,p,v(i).all,rcond,verbose);
    end loop;
  end Run_SVD_Newton;

  procedure Run_SVD_Newton
              ( file : in file_type; maxdeg,nbrit : in integer32;
                p : in DoblDobl_CSeries_Poly_Systems.Poly_Sys;
                v : in DoblDobl_Complex_Series_VecVecs.VecVec;
                verbose : in boolean := false ) is

    rcond : double_double;

  begin
    for i in v'range loop
      if verbose then
        put(file,"Running on solution ");
        put(file,i,1); put_line(file," ...");
        declare
          lvi : DoblDobl_Complex_Series_Vectors.Link_to_Vector := v(i);
          sol : DoblDobl_Complex_Vectors.Vector(lvi'range);
        begin
          for k in sol'range loop
            sol(k) := lvi(k).cff(0);
          end loop;
          DoblDobl_Complex_Solutions_io.put_vector(file,sol);
        end;
      end if;
      Run_SVD_Newton(file,maxdeg,nbrit,p,v(i).all,rcond,verbose);
    end loop;
  end Run_SVD_Newton;

  procedure Run_SVD_Newton
              ( file : in file_type; maxdeg,nbrit : in integer32;
                p : in QuadDobl_CSeries_Poly_Systems.Poly_Sys;
                v : in QuadDobl_Complex_Series_VecVecs.VecVec;
                verbose : in boolean := false ) is

    rcond : quad_double;

  begin
    for i in v'range loop
      if verbose then
        put(file,"Running on solution ");
        put(file,i,1); put_line(file," ...");
        declare
          lvi : QuadDobl_Complex_Series_Vectors.Link_to_Vector := v(i);
          sol : QuadDobl_Complex_Vectors.Vector(lvi'range);
        begin
          for k in sol'range loop
            sol(k) := lvi(k).cff(0);
          end loop;
          QuadDobl_Complex_Solutions_io.put_vector(file,sol);
        end;
      end if;
      Run_SVD_Newton(file,maxdeg,nbrit,p,v(i).all,rcond,verbose);
    end loop;
  end Run_SVD_Newton;

  procedure Run_Echelon_Newton
              ( maxdeg,nbrit : in integer32;
                p : in Standard_CSeries_Poly_Systems.Poly_Sys;
                v : in Standard_Complex_Series_VecVecs.VecVec;
                verbose : in boolean := false;
                pause : in boolean := false ) is

    ans : character;
    det : Standard_Complex_Numbers.Complex_Number;

  begin
    for i in v'range loop
      if verbose then
        put("Running on solution "); put(i,1); put_line(" ...");
        declare
          lvi : Standard_Complex_Series_Vectors.Link_to_Vector := v(i);
          sol : Standard_Complex_Vectors.Vector(lvi'range);
        begin
          for k in sol'range loop
            sol(k) := lvi(k).cff(0);
          end loop;
          Standard_Complex_Solutions_io.put_vector(sol);
        end;
      end if;
      Run_Echelon_Newton(maxdeg,nbrit,p,v(i).all,det,verbose);
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
                pause : in boolean := false ) is

    ans : character;
    det : DoblDobl_Complex_Numbers.Complex_Number;

  begin
    for i in v'range loop
      if verbose then
        put("Running on solution "); put(i,1); put_line(" ...");
        declare
          lvi : DoblDobl_Complex_Series_Vectors.Link_to_Vector := v(i);
          sol : DoblDobl_Complex_Vectors.Vector(lvi'range);
        begin
          for k in sol'range loop
            sol(k) := lvi(k).cff(0);
          end loop;
          DoblDobl_Complex_Solutions_io.put_vector(sol);
        end;
      end if;
      Run_Echelon_Newton(maxdeg,nbrit,p,v(i).all,det,verbose);
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
                pause : in boolean := false ) is

    ans : character;
    det : QuadDobl_Complex_Numbers.Complex_Number;

  begin
    for i in v'range loop
      if verbose then
        put("Running on solution "); put(i,1); put_line(" ...");
        declare
          lvi : QuadDobl_Complex_Series_Vectors.Link_to_Vector := v(i);
          sol : QuadDobl_Complex_Vectors.Vector(lvi'range);
        begin
          for k in sol'range loop
            sol(k) := lvi(k).cff(0);
          end loop;
          QuadDobl_Complex_Solutions_io.put_vector(sol);
        end;
      end if;
      Run_Echelon_Newton(maxdeg,nbrit,p,v(i).all,det,verbose);
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
                verbose : in boolean := false ) is

    det : Standard_Complex_Numbers.Complex_Number;

  begin
    for i in v'range loop
      if verbose then
        put(file,"Running on solution ");
        put(file,i,1); put_line(file," ...");
        declare
          lvi : Standard_Complex_Series_Vectors.Link_to_Vector := v(i);
          sol : Standard_Complex_Vectors.Vector(lvi'range);
        begin
          for k in sol'range loop
            sol(k) := lvi(k).cff(0);
          end loop;
          Standard_Complex_Solutions_io.put_vector(file,sol);
        end;
      end if;
      Run_Echelon_Newton(file,maxdeg,nbrit,p,v(i).all,det,verbose);
    end loop;
  end Run_Echelon_Newton;

  procedure Run_Echelon_Newton
              ( file : in file_type; maxdeg,nbrit : in integer32;
                p : in DoblDobl_CSeries_Poly_Systems.Poly_Sys;
                v : in DoblDobl_Complex_Series_VecVecs.VecVec;
                verbose : in boolean := false ) is

    det : DoblDobl_Complex_Numbers.Complex_Number;

  begin
    for i in v'range loop
      if verbose then
        put(file,"Running on solution ");
        put(file,i,1); put_line(file," ...");
        declare
          lvi : DoblDobl_Complex_Series_Vectors.Link_to_Vector := v(i);
          sol : DoblDobl_Complex_Vectors.Vector(lvi'range);
        begin
          for k in sol'range loop
            sol(k) := lvi(k).cff(0);
          end loop;
          DoblDobl_Complex_Solutions_io.put_vector(file,sol);
        end;
      end if;
      Run_Echelon_Newton(file,maxdeg,nbrit,p,v(i).all,det,verbose);
    end loop;
  end Run_Echelon_Newton;

  procedure Run_Echelon_Newton
              ( file : in file_type; maxdeg,nbrit : in integer32;
                p : in QuadDobl_CSeries_Poly_Systems.Poly_Sys;
                v : in QuadDobl_Complex_Series_VecVecs.VecVec;
                verbose : in boolean := false ) is

    det : QuadDobl_Complex_Numbers.Complex_Number;

  begin
    for i in v'range loop
      if verbose then
        put(file,"Running on solution ");
        put(file,i,1); put_line(file," ...");
        declare
          lvi : QuadDobl_Complex_Series_Vectors.Link_to_Vector := v(i);
          sol : QuadDobl_Complex_Vectors.Vector(lvi'range);
        begin
          for k in sol'range loop
            sol(k) := lvi(k).cff(0);
          end loop;
          QuadDobl_Complex_Solutions_io.put_vector(file,sol);
        end;
      end if;
      Run_Echelon_Newton(file,maxdeg,nbrit,p,v(i).all,det,verbose);
    end loop;
  end Run_Echelon_Newton;

end Power_Series_Methods;
