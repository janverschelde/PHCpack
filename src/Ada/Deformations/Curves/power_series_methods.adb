with Communications_with_User;           use Communications_with_User;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;
with Double_Double_Numbers_io;           use Double_Double_Numbers_io;
with Quad_Double_Numbers_io;             use Quad_Double_Numbers_io;
with Standard_Complex_Vectors;
with Standard_Complex_Solutions_io;
with DoblDobl_Complex_Vectors;
with DoblDobl_Complex_Solutions_io;
with QuadDobl_Complex_Vectors;
with QuadDobl_Complex_Solutions_io;
with Standard_Series_Poly_SysFun;
with DoblDobl_Series_Poly_SysFun;
with QuadDobl_Series_Poly_SysFun;
with Series_and_Polynomials;
with Series_and_Polynomials_io;
with Standard_Newton_Matrix_Series; -- with Standard_Newton_Series;
with DoblDobl_Newton_Matrix_Series; -- with DoblDobl_Newton_Series;
with QuadDobl_Newton_Matrix_Series; -- with QuadDobl_Newton_Series;

package body Power_Series_Methods is

  procedure Run_LU_Newton
             ( nbrit : in integer32;
               p : in Standard_Series_Poly_Systems.Poly_Sys;
               s : in out Standard_Dense_Series_Vectors.Vector;
               info : out integer32; verbose : in boolean := false ) is
  begin
    Run_LU_Newton(standard_output,nbrit,p,s,info,verbose);
  end Run_LU_Newton;

  procedure Run_LU_Newton
             ( file : in file_type; nbrit : in integer32;
               p : in Standard_Series_Poly_Systems.Poly_Sys;
               s : in out Standard_Dense_Series_Vectors.Vector;
               info : out integer32; verbose : in boolean := false ) is

    use Standard_Newton_Matrix_Series; -- use Standard_Newton_Series;

    order : integer32 := 1;
    tol : constant double_float := 1.0E-12;
    eva : Standard_Dense_Series_Vectors.Vector(p'range);

  begin
    if not verbose then
      LU_Newton_Steps(p,order,nbrit,s,info);
    else
      LU_Newton_Steps(file,p,order,nbrit,s,info);
      if info /= 0 then
        put(file,"info = "); put(file,info,1); new_line(file);
      else
        Series_and_Polynomials.Filter(s,tol);
        put_line(file,"The updated power series solution :");
        Series_and_Polynomials_io.put(file,s);
        eva := Standard_Series_Poly_SysFun.Eval(p,s);
        Series_and_Polynomials.Filter(eva,tol);
        put_line(file,"The evaluated solution :");
        Series_and_Polynomials_io.put(file,eva);
      end if;
    end if;
  end Run_LU_Newton;

  procedure Run_LU_Newton
             ( nbrit : in integer32;
               p : in DoblDobl_Series_Poly_Systems.Poly_Sys;
               s : in out DoblDobl_Dense_Series_Vectors.Vector;
               info : out integer32; verbose : in boolean := false ) is
  begin
    Run_LU_Newton(standard_output,nbrit,p,s,info,verbose);
  end Run_LU_Newton;

  procedure Run_LU_Newton
             ( file : in file_type; nbrit : in integer32;
               p : in DoblDobl_Series_Poly_Systems.Poly_Sys;
               s : in out DoblDobl_Dense_Series_Vectors.Vector;
               info : out integer32; verbose : in boolean := false ) is

    use DoblDobl_Newton_Matrix_Series; -- use DoblDobl_Newton_Series;

    order : integer32 := 1;
    tol : constant double_float := 1.0E-24;
    eva : DoblDobl_Dense_Series_Vectors.Vector(p'range);

  begin
    if not verbose then
      LU_Newton_Steps(p,order,nbrit,s,info);
    else
      LU_Newton_Steps(file,p,order,nbrit,s,info);
      if info /= 0 then
        put(file,"info = "); put(file,info,1); new_line(file);
      else
        Series_and_Polynomials.Filter(s,tol);
        put_line(file,"The updated power series solution :");
        Series_and_Polynomials_io.put(file,s);
        eva := DoblDobl_Series_Poly_SysFun.Eval(p,s);
        Series_and_Polynomials.Filter(eva,tol);
        put_line(file,"The evaluated solution :");
        Series_and_Polynomials_io.put(file,eva);
      end if;
    end if;
  end Run_LU_Newton;

  procedure Run_LU_Newton
             ( nbrit : in integer32;
               p : in QuadDobl_Series_Poly_Systems.Poly_Sys;
               s : in out QuadDobl_Dense_Series_Vectors.Vector;
               info : out integer32; verbose : in boolean := false ) is
  begin
    Run_LU_Newton(standard_output,nbrit,p,s,info,verbose);
  end Run_LU_Newton;

  procedure Run_LU_Newton
             ( file : in file_type; nbrit : in integer32;
               p : in QuadDobl_Series_Poly_Systems.Poly_Sys;
               s : in out QuadDobl_Dense_Series_Vectors.Vector;
               info : out integer32; verbose : in boolean := false ) is

    use QuadDobl_Newton_Matrix_Series; -- use QuadDobl_Newton_Series;

    order : integer32 := 1;
    tol : constant double_float := 1.0E-32;
    eva : QuadDobl_Dense_Series_Vectors.Vector(p'range);

  begin
    if not verbose then
      LU_Newton_Steps(p,order,nbrit,s,info);
    else
      LU_Newton_Steps(file,p,order,nbrit,s,info);
      if info /= 0 then
        put(file,"info = "); put(file,info,1); new_line(file);
      else
        Series_and_Polynomials.Filter(s,tol);
        put_line(file,"The updated power series solution :");
        Series_and_Polynomials_io.put(file,s);
        eva := QuadDobl_Series_Poly_SysFun.Eval(p,s);
        Series_and_Polynomials.Filter(eva,tol);
        put_line(file,"The evaluated solution :");
        Series_and_Polynomials_io.put(file,eva);
      end if;
    end if;
  end Run_LU_Newton;

  procedure Run_LU_Newton
             ( nbrit : in integer32;
               p : in Standard_Series_Poly_Systems.Poly_Sys;
               s : in out Standard_Dense_Series_Vectors.Vector;
               rcond : out double_float; verbose : in boolean := false ) is
  begin
    Run_LU_Newton(standard_output,nbrit,p,s,rcond,verbose);
  end Run_LU_Newton;

  procedure Run_LU_Newton
             ( file : in file_type; nbrit : in integer32;
               p : in Standard_Series_Poly_Systems.Poly_Sys;
               s : in out Standard_Dense_Series_Vectors.Vector;
               rcond : out double_float; verbose : in boolean := false ) is

    use Standard_Newton_Matrix_Series; -- use Standard_Newton_Series;

    order : integer32 := 1;
    tol : constant double_float := 1.0E-12;
    eva : Standard_Dense_Series_Vectors.Vector(p'range);

  begin
    if not verbose then
      LU_Newton_Steps(p,order,nbrit,s,rcond);
    else
      LU_Newton_Steps(file,p,order,nbrit,s,rcond);
      put(file,"rcond : "); put(file,rcond,3); new_line(file);
      if 1.0 + rcond /= 1.0 then
        Series_and_Polynomials.Filter(s,tol);
        put_line(file,"The updated power series solution :");
        Series_and_Polynomials_io.put(file,s);
        eva := Standard_Series_Poly_SysFun.Eval(p,s);
        Series_and_Polynomials.Filter(eva,tol);
        put_line(file,"The evaluated solution :");
        Series_and_Polynomials_io.put(file,eva);
      end if;
    end if;
  end Run_LU_Newton;

  procedure Run_LU_Newton
             ( nbrit : in integer32;
               p : in DoblDobl_Series_Poly_Systems.Poly_Sys;
               s : in out DoblDobl_Dense_Series_Vectors.Vector;
               rcond : out double_double; verbose : in boolean := false ) is
  begin
    Run_LU_Newton(standard_output,nbrit,p,s,rcond,verbose);
  end Run_LU_Newton;

  procedure Run_LU_Newton
             ( file : in file_type; nbrit : in integer32;
               p : in DoblDobl_Series_Poly_Systems.Poly_Sys;
               s : in out DoblDobl_Dense_Series_Vectors.Vector;
               rcond : out double_double; verbose : in boolean := false ) is

    use DoblDobl_Newton_Matrix_Series; -- use DoblDobl_Newton_Series;

    order : integer32 := 1;
    tol : constant double_float := 1.0E-24;
    eva : DoblDobl_Dense_Series_Vectors.Vector(p'range);
    one : constant double_double := create(1.0);

  begin
    if not verbose then
      LU_Newton_Steps(p,order,nbrit,s,rcond);
    else
      LU_Newton_Steps(file,p,order,nbrit,s,rcond);
      put(file,"rcond : "); put(file,rcond,3); new_line(file);
      if one + rcond /= one then
        Series_and_Polynomials.Filter(s,tol);
        put_line(file,"The updated power series solution :");
        Series_and_Polynomials_io.put(file,s);
        eva := DoblDobl_Series_Poly_SysFun.Eval(p,s);
        Series_and_Polynomials.Filter(eva,tol);
        put_line(file,"The evaluated solution :");
        Series_and_Polynomials_io.put(file,eva);
      end if;
    end if;
  end Run_LU_Newton;

  procedure Run_LU_Newton
             ( nbrit : in integer32;
               p : in QuadDobl_Series_Poly_Systems.Poly_Sys;
               s : in out QuadDobl_Dense_Series_Vectors.Vector;
               rcond : out quad_double; verbose : in boolean := false ) is
  begin
    Run_LU_Newton(standard_output,nbrit,p,s,rcond,verbose);
  end Run_LU_Newton;

  procedure Run_LU_Newton
             ( file : in file_type; nbrit : in integer32;
               p : in QuadDobl_Series_Poly_Systems.Poly_Sys;
               s : in out QuadDobl_Dense_Series_Vectors.Vector;
               rcond : out quad_double; verbose : in boolean := false ) is

    use QuadDobl_Newton_Matrix_Series; -- use QuadDobl_Newton_Series;

    order : integer32 := 1;
    tol : constant double_float := 1.0E-32;
    eva : QuadDobl_Dense_Series_Vectors.Vector(p'range);
    one : constant quad_double := create(1.0);

  begin
    if not verbose then
      LU_Newton_Steps(p,order,nbrit,s,rcond);
    else
      LU_Newton_Steps(file,p,order,nbrit,s,rcond);
      put(file,"rcond : "); put(file,rcond,3); new_line(file);
      if one + rcond /= one then
        Series_and_Polynomials.Filter(s,tol);
        put_line(file,"The updated power series solution :");
        Series_and_Polynomials_io.put(file,s);
        eva := QuadDobl_Series_Poly_SysFun.Eval(p,s);
        Series_and_Polynomials.Filter(eva,tol);
        put_line(file,"The evaluated solution :");
        Series_and_Polynomials_io.put(file,eva);
      end if;
    end if;
  end Run_LU_Newton;

  procedure Run_QR_Newton
             ( nbrit : in integer32;
               p : in Standard_Series_Poly_Systems.Poly_Sys;
               s : in out Standard_Dense_Series_Vectors.Vector;
               verbose : in boolean := false ) is
  begin
    Run_QR_Newton(standard_output,nbrit,p,s,verbose);
  end Run_QR_Newton;

  procedure Run_QR_Newton
             ( file : in file_type; nbrit : in integer32;
               p : in Standard_Series_Poly_Systems.Poly_Sys;
               s : in out Standard_Dense_Series_Vectors.Vector;
               verbose : in boolean := false ) is

    use Standard_Newton_Matrix_Series; -- use Standard_Newton_Series;

    order : integer32 := 1;
    info : integer32;
    tol : constant double_float := 1.0E-12;
    eva : Standard_Dense_Series_Vectors.Vector(p'range);

  begin
    if not verbose then
      QR_Newton_Steps(p,order,nbrit,s,info);
    else
      QR_Newton_Steps(file,p,order,nbrit,s,info);
      if info /= 0 then
        put(file,"info = "); put(file,info,1); new_line(file);
      else
        Series_and_Polynomials.Filter(s,tol);
        put_line(file,"The updated power series solution :");
        Series_and_Polynomials_io.put(file,s);
        eva := Standard_Series_Poly_SysFun.Eval(p,s);
        Series_and_Polynomials.Filter(eva,tol);
        put_line(file,"The evaluated solution :");
        Series_and_Polynomials_io.put(file,eva);
      end if;
    end if;
  end Run_QR_Newton;

  procedure Run_QR_Newton
             ( nbrit : in integer32;
               p : in DoblDobl_Series_Poly_Systems.Poly_Sys;
               s : in out DoblDobl_Dense_Series_Vectors.Vector;
               verbose : in boolean := false ) is
  begin
    Run_QR_Newton(standard_output,nbrit,p,s,verbose);
  end Run_QR_Newton;

  procedure Run_QR_Newton
             ( file : in file_type; nbrit : in integer32;
               p : in DoblDobl_Series_Poly_Systems.Poly_Sys;
               s : in out DoblDobl_Dense_Series_Vectors.Vector;
               verbose : in boolean := false ) is

    use DoblDobl_Newton_Matrix_Series; -- use DoblDobl_Newton_Series;

    order : integer32 := 1;
    info : integer32;
    tol : constant double_float := 1.0E-24;
    eva : DoblDobl_Dense_Series_Vectors.Vector(p'range);

  begin
    if not verbose then
      QR_Newton_Steps(p,order,nbrit,s,info);
    else
      QR_Newton_Steps(file,p,order,nbrit,s,info);
      if info /= 0 then
        put(file,"info = "); put(file,info,1); new_line(file);
      else
        Series_and_Polynomials.Filter(s,tol);
        put_line(file,"The updated power series solution :");
        Series_and_Polynomials_io.put(file,s);
        eva := DoblDobl_Series_Poly_SysFun.Eval(p,s);
        Series_and_Polynomials.Filter(eva,tol);
        put_line(file,"The evaluated solution :");
        Series_and_Polynomials_io.put(file,eva);
      end if;
    end if;
  end Run_QR_Newton;

  procedure Run_QR_Newton
             ( nbrit : in integer32;
               p : in QuadDobl_Series_Poly_Systems.Poly_Sys;
               s : in out QuadDobl_Dense_Series_Vectors.Vector;
               verbose : in boolean := false ) is
  begin
    Run_QR_Newton(standard_output,nbrit,p,s,verbose);
  end Run_QR_Newton;

  procedure Run_QR_Newton
             ( file : in file_type; nbrit : in integer32;
               p : in QuadDobl_Series_Poly_Systems.Poly_Sys;
               s : in out QuadDobl_Dense_Series_Vectors.Vector;
               verbose : in boolean := false ) is

    use QuadDobl_Newton_Matrix_Series; -- use QuadDobl_Newton_Series;

    order : integer32 := 1;
    info : integer32;
    tol : constant double_float := 1.0E-32;
    eva : QuadDobl_Dense_Series_Vectors.Vector(p'range);

  begin
    if not verbose then
      QR_Newton_Steps(p,order,nbrit,s,info);
    else
      QR_Newton_Steps(file,p,order,nbrit,s,info);
      if info /= 0 then
        put(file,"info = "); put(file,info,1); new_line(file);
      else
        Series_and_Polynomials.Filter(s,tol);
        put_line(file,"The updated power series solution :");
        Series_and_Polynomials_io.put(file,s);
        eva := QuadDobl_Series_Poly_SysFun.Eval(p,s);
        Series_and_Polynomials.Filter(eva,tol);
        put_line(file,"The evaluated solution :");
        Series_and_Polynomials_io.put(file,eva);
      end if;
    end if;
  end Run_QR_Newton;

  procedure Run_SVD_Newton
             ( nbrit : in integer32;
               p : in Standard_Series_Poly_Systems.Poly_Sys;
               s : in out Standard_Dense_Series_Vectors.Vector;
               rcond : out double_float; verbose : in boolean := false ) is
  begin
    Run_SVD_Newton(standard_output,nbrit,p,s,rcond,verbose);
  end Run_SVD_Newton;

  procedure Run_SVD_Newton
             ( file : in file_type; nbrit : in integer32;
               p : in Standard_Series_Poly_Systems.Poly_Sys;
               s : in out Standard_Dense_Series_Vectors.Vector;
               rcond : out double_float; verbose : in boolean := false ) is

    use Standard_Newton_Matrix_Series; -- use Standard_Newton_Series;

    order : integer32 := 1;
    info : integer32;
    tol : constant double_float := 1.0E-12;
    eva : Standard_Dense_Series_Vectors.Vector(p'range);

  begin
    if not verbose then
      SVD_Newton_Steps(p,order,nbrit,s,info,rcond);
    else
      SVD_Newton_Steps(file,p,order,nbrit,s,info,rcond);
      put(file,"rcond = "); put(file,rcond,3); new_line(file);
      if 1.0 + rcond /= 1.0 then
        Series_and_Polynomials.Filter(s,tol);
        put_line(file,"The updated power series solution :");
        Series_and_Polynomials_io.put(file,s);
        eva := Standard_Series_Poly_SysFun.Eval(p,s);
        Series_and_Polynomials.Filter(eva,tol);
        put_line(file,"The evaluated solution :");
        Series_and_Polynomials_io.put(file,eva);
      end if;
    end if;
  end Run_SVD_Newton;

  procedure Run_SVD_Newton
             ( nbrit : in integer32;
               p : in DoblDobl_Series_Poly_Systems.Poly_Sys;
               s : in out DoblDobl_Dense_Series_Vectors.Vector;
               rcond : out double_double; verbose : in boolean := false ) is
  begin
    Run_SVD_Newton(standard_output,nbrit,p,s,rcond,verbose);
  end Run_SVD_Newton;

  procedure Run_SVD_Newton
             ( file : in file_type; nbrit : in integer32;
               p : in DoblDobl_Series_Poly_Systems.Poly_Sys;
               s : in out DoblDobl_Dense_Series_Vectors.Vector;
               rcond : out double_double; verbose : in boolean := false ) is

    use DoblDobl_Newton_Matrix_Series; -- use DoblDobl_Newton_Series;

    order : integer32 := 1;
    info : integer32;
    tol : constant double_float := 1.0E-24;
    eva : DoblDobl_Dense_Series_Vectors.Vector(p'range);
    one : constant double_double := create(1.0);

  begin
    if not verbose then
      SVD_Newton_Steps(p,order,nbrit,s,info,rcond);
    else
      SVD_Newton_Steps(file,p,order,nbrit,s,info,rcond);
      put(file,"rcond = "); put(file,rcond,3); new_line(file);
      if one + rcond /= one then 
        Series_and_Polynomials.Filter(s,tol);
        put_line(file,"The updated power series solution :");
        Series_and_Polynomials_io.put(file,s);
        eva := DoblDobl_Series_Poly_SysFun.Eval(p,s);
        Series_and_Polynomials.Filter(eva,tol);
        put_line(file,"The evaluated solution :");
        Series_and_Polynomials_io.put(file,eva);
      end if;
    end if;
  end Run_SVD_Newton;

  procedure Run_SVD_Newton
             ( nbrit : in integer32;
               p : in QuadDobl_Series_Poly_Systems.Poly_Sys;
               s : in out QuadDobl_Dense_Series_Vectors.Vector;
               rcond : out quad_double; verbose : in boolean := false ) is
  begin
    Run_SVD_Newton(standard_output,nbrit,p,s,rcond,verbose);
  end Run_SVD_Newton;

  procedure Run_SVD_Newton
             ( file : in file_type; nbrit : in integer32;
               p : in QuadDobl_Series_Poly_Systems.Poly_Sys;
               s : in out QuadDobl_Dense_Series_Vectors.Vector;
               rcond : out quad_double; verbose : in boolean := false ) is

    use QuadDobl_Newton_Matrix_Series; -- use QuadDobl_Newton_Series;

    order : integer32 := 1;
    info : integer32;
    tol : constant double_float := 1.0E-32;
    eva : QuadDobl_Dense_Series_Vectors.Vector(p'range);
    one : constant quad_double := create(1.0);

  begin
    if not verbose then
      SVD_Newton_Steps(p,order,nbrit,s,info,rcond);
    else
      SVD_Newton_Steps(file,p,order,nbrit,s,info,rcond);
      put(file,"rcond = "); put(file,rcond,3); new_line(file);
      if one + rcond /= one then
        Series_and_Polynomials.Filter(s,tol);
        put_line(file,"The updated power series solution :");
        Series_and_Polynomials_io.put(file,s);
        eva := QuadDobl_Series_Poly_SysFun.Eval(p,s);
        Series_and_Polynomials.Filter(eva,tol);
        put_line(file,"The evaluated solution :");
        Series_and_Polynomials_io.put(file,eva);
      end if;
    end if;
  end Run_SVD_Newton;

  procedure Run_LU_Newton
              ( nbrit : in integer32;
                p : in Standard_Series_Poly_Systems.Poly_Sys;
                v : in Standard_Dense_Series_VecVecs.VecVec;
                verbose : in boolean := false;
                pause : in boolean := false ) is

    ans : character;
    info : integer32;

  begin
    for i in v'range loop
      if verbose then
        put("Running on solution "); put(i,1); put_line(" ...");
        declare
          lvi : Standard_Dense_Series_Vectors.Link_to_Vector := v(i);
          sol : Standard_Complex_Vectors.Vector(lvi'range);
        begin
          for k in sol'range loop
            sol(k) := lvi(k).cff(0);
          end loop;
          Standard_Complex_Solutions_io.put_vector(sol);
        end;
      end if;
      Run_LU_Newton(nbrit,p,v(i).all,info,verbose);
      if verbose then
        if pause then
          put("Continue ? (y/n) "); Ask_Yes_or_No(ans);
          exit when (ans /= 'y');
        end if;
      end if;
    end loop;
  end Run_LU_Newton;

  procedure Run_LU_Newton
              ( nbrit : in integer32;
                p : in DoblDobl_Series_Poly_Systems.Poly_Sys;
                v : in DoblDobl_Dense_Series_VecVecs.VecVec;
                verbose : in boolean := false;
                pause : in boolean := false ) is

    ans : character;
    info : integer32;

  begin
    for i in v'range loop
      if verbose then
        put("Running on solution "); put(i,1); put_line(" ...");
        declare
          lvi : DoblDobl_Dense_Series_Vectors.Link_to_Vector := v(i);
          sol : DoblDobl_Complex_Vectors.Vector(lvi'range);
        begin
          for k in sol'range loop
            sol(k) := lvi(k).cff(0);
          end loop;
          DoblDobl_Complex_Solutions_io.put_vector(sol);
        end;
      end if;
      Run_LU_Newton(nbrit,p,v(i).all,info,verbose);
      if verbose then
        if pause then
          put("Continue ? (y/n) "); Ask_Yes_or_No(ans);
          exit when (ans /= 'y');
        end if;
      end if;
    end loop;
  end Run_LU_Newton;

  procedure Run_LU_Newton
              ( nbrit : in integer32;
                p : in QuadDobl_Series_Poly_Systems.Poly_Sys;
                v : in QuadDobl_Dense_Series_VecVecs.VecVec;
                verbose : in boolean := false;
                pause : in boolean := false ) is

    ans : character;
    info : integer32;

  begin
    for i in v'range loop
      if verbose then
        put("Running on solution "); put(i,1); put_line(" ...");
        declare
          lvi : QuadDobl_Dense_Series_Vectors.Link_to_Vector := v(i);
          sol : QuadDobl_Complex_Vectors.Vector(lvi'range);
        begin
          for k in sol'range loop
            sol(k) := lvi(k).cff(0);
          end loop;
          QuadDobl_Complex_Solutions_io.put_vector(sol);
        end;
      end if;
      Run_LU_Newton(nbrit,p,v(i).all,info,verbose);
      if verbose then
        if pause then
          put("Continue ? (y/n) "); Ask_Yes_or_No(ans);
          exit when (ans /= 'y');
        end if;
      end if;
    end loop;
  end Run_LU_Newton;

  procedure Run_LU_Newton
              ( file : in file_type; nbrit : in integer32;
                p : in Standard_Series_Poly_Systems.Poly_Sys;
                v : in Standard_Dense_Series_VecVecs.VecVec;
                verbose : in boolean := false ) is

    info : integer32;

  begin
    for i in v'range loop
      if verbose then
        put(file,"Running on solution ");
        put(file,i,1); put_line(file," ...");
        declare
          lvi : Standard_Dense_Series_Vectors.Link_to_Vector := v(i);
          sol : Standard_Complex_Vectors.Vector(lvi'range);
        begin
          for k in sol'range loop
            sol(k) := lvi(k).cff(0);
          end loop;
          Standard_Complex_Solutions_io.put_vector(file,sol);
        end;
      end if;
      Run_LU_Newton(file,nbrit,p,v(i).all,info,verbose);
    end loop;
  end Run_LU_Newton;

  procedure Run_LU_Newton
              ( file : in file_type; nbrit : in integer32;
                p : in DoblDobl_Series_Poly_Systems.Poly_Sys;
                v : in DoblDobl_Dense_Series_VecVecs.VecVec;
                verbose : in boolean := false ) is

    info : integer32;

  begin
    for i in v'range loop
      if verbose then
        put(file,"Running on solution ");
        put(file,i,1); put_line(file," ...");
        declare
          lvi : DoblDobl_Dense_Series_Vectors.Link_to_Vector := v(i);
          sol : DoblDobl_Complex_Vectors.Vector(lvi'range);
        begin
          for k in sol'range loop
            sol(k) := lvi(k).cff(0);
          end loop;
          DoblDobl_Complex_Solutions_io.put_vector(file,sol);
        end;
      end if;
      Run_LU_Newton(file,nbrit,p,v(i).all,info,verbose);
    end loop;
  end Run_LU_Newton;

  procedure Run_LU_Newton
              ( file : in file_type; nbrit : in integer32;
                p : in QuadDobl_Series_Poly_Systems.Poly_Sys;
                v : in QuadDobl_Dense_Series_VecVecs.VecVec;
                verbose : in boolean := false ) is

    info : integer32;

  begin
    for i in v'range loop
      if verbose then
        put(file,"Running on solution ");
        put(file,i,1); put_line(file," ...");
        declare
          lvi : QuadDobl_Dense_Series_Vectors.Link_to_Vector := v(i);
          sol : QuadDobl_Complex_Vectors.Vector(lvi'range);
        begin
          for k in sol'range loop
            sol(k) := lvi(k).cff(0);
          end loop;
          QuadDobl_Complex_Solutions_io.put_vector(file,sol);
        end;
      end if;
      Run_LU_Newton(file,nbrit,p,v(i).all,info,verbose);
    end loop;
  end Run_LU_Newton;

  procedure Run_QR_Newton
              ( nbrit : in integer32;
                p : in Standard_Series_Poly_Systems.Poly_Sys;
                v : in Standard_Dense_Series_VecVecs.VecVec;
                verbose : in boolean := false;
                pause : in boolean := false ) is

    ans : character;

  begin
    for i in v'range loop
      if verbose then
        put("Running on solution "); put(i,1); put_line(" ...");
        declare
          lvi : Standard_Dense_Series_Vectors.Link_to_Vector := v(i);
          sol : Standard_Complex_Vectors.Vector(lvi'range);
        begin
          for k in sol'range loop
            sol(k) := lvi(k).cff(0);
          end loop;
          Standard_Complex_Solutions_io.put_vector(sol);
        end;
      end if;
      Run_QR_Newton(nbrit,p,v(i).all,verbose);
      if verbose then
        if pause then
          put("Continue ? (y/n) "); Ask_Yes_or_No(ans);
          exit when (ans /= 'y');
        end if;
      end if;
    end loop;
  end Run_QR_Newton;

  procedure Run_QR_Newton
              ( nbrit : in integer32;
                p : in DoblDobl_Series_Poly_Systems.Poly_Sys;
                v : in DoblDobl_Dense_Series_VecVecs.VecVec;
                verbose : in boolean := false;
                pause : in boolean := false ) is

    ans : character;

  begin
    for i in v'range loop
      if verbose then
        put("Running on solution "); put(i,1); put_line(" ...");
        declare
          lvi : DoblDobl_Dense_Series_Vectors.Link_to_Vector := v(i);
          sol : DoblDobl_Complex_Vectors.Vector(lvi'range);
        begin
          for k in sol'range loop
            sol(k) := lvi(k).cff(0);
          end loop;
          DoblDobl_Complex_Solutions_io.put_vector(sol);
        end;
      end if;
      Run_QR_Newton(nbrit,p,v(i).all,verbose);
      if verbose then
        if pause then
          put("Continue ? (y/n) "); Ask_Yes_or_No(ans);
          exit when (ans /= 'y');
        end if;
      end if;
    end loop;
  end Run_QR_Newton;

  procedure Run_QR_Newton
              ( nbrit : in integer32;
                p : in QuadDobl_Series_Poly_Systems.Poly_Sys;
                v : in QuadDobl_Dense_Series_VecVecs.VecVec;
                verbose : in boolean := false;
                pause : in boolean := false ) is

    ans : character;

  begin
    for i in v'range loop
      if verbose then
        put("Running on solution "); put(i,1); put_line(" ...");
        declare
          lvi : QuadDobl_Dense_Series_Vectors.Link_to_Vector := v(i);
          sol : QuadDobl_Complex_Vectors.Vector(lvi'range);
        begin
          for k in sol'range loop
            sol(k) := lvi(k).cff(0);
          end loop;
          QuadDobl_Complex_Solutions_io.put_vector(sol);
        end;
      end if;
      Run_QR_Newton(nbrit,p,v(i).all,verbose);
      if verbose then
        if pause then
          put("Continue ? (y/n) "); Ask_Yes_or_No(ans);
          exit when (ans /= 'y');
        end if;
      end if;
    end loop;
  end Run_QR_Newton;

  procedure Run_QR_Newton
              ( file : in file_type; nbrit : in integer32;
                p : in Standard_Series_Poly_Systems.Poly_Sys;
                v : in Standard_Dense_Series_VecVecs.VecVec;
                verbose : in boolean := false ) is
  begin
    for i in v'range loop
      if verbose then
        put(file,"Running on solution ");
        put(file,i,1); put_line(file," ...");
        declare
          lvi : Standard_Dense_Series_Vectors.Link_to_Vector := v(i);
          sol : Standard_Complex_Vectors.Vector(lvi'range);
        begin
          for k in sol'range loop
            sol(k) := lvi(k).cff(0);
          end loop;
          Standard_Complex_Solutions_io.put_vector(file,sol);
        end;
      end if;
      Run_QR_Newton(file,nbrit,p,v(i).all,verbose);
    end loop;
  end Run_QR_Newton;

  procedure Run_QR_Newton
              ( file : in file_type; nbrit : in integer32;
                p : in DoblDobl_Series_Poly_Systems.Poly_Sys;
                v : in DoblDobl_Dense_Series_VecVecs.VecVec;
                verbose : in boolean := false ) is
  begin
    for i in v'range loop
      if verbose then
        put(file,"Running on solution ");
        put(file,i,1); put_line(file," ...");
        declare
          lvi : DoblDobl_Dense_Series_Vectors.Link_to_Vector := v(i);
          sol : DoblDobl_Complex_Vectors.Vector(lvi'range);
        begin
          for k in sol'range loop
            sol(k) := lvi(k).cff(0);
          end loop;
          DoblDobl_Complex_Solutions_io.put_vector(file,sol);
        end;
      end if;
      Run_QR_Newton(file,nbrit,p,v(i).all,verbose);
    end loop;
  end Run_QR_Newton;

  procedure Run_QR_Newton
              ( file : in file_type; nbrit : in integer32;
                p : in QuadDobl_Series_Poly_Systems.Poly_Sys;
                v : in QuadDobl_Dense_Series_VecVecs.VecVec;
                verbose : in boolean := false ) is
  begin
    for i in v'range loop
      if verbose then
        put(file,"Running on solution ");
        put(file,i,1); put_line(file," ...");
        declare
          lvi : QuadDobl_Dense_Series_Vectors.Link_to_Vector := v(i);
          sol : QuadDobl_Complex_Vectors.Vector(lvi'range);
        begin
          for k in sol'range loop
            sol(k) := lvi(k).cff(0);
          end loop;
          QuadDobl_Complex_Solutions_io.put_vector(file,sol);
        end;
      end if;
      Run_QR_Newton(file,nbrit,p,v(i).all,verbose);
    end loop;
  end Run_QR_Newton;

  procedure Run_SVD_Newton
              ( nbrit : in integer32;
                p : in Standard_Series_Poly_Systems.Poly_Sys;
                v : in Standard_Dense_Series_VecVecs.VecVec;
                verbose : in boolean := false;
                pause : in boolean := false ) is

    ans : character;
    rcond : double_float;

  begin
    for i in v'range loop
      if verbose then
        put("Running on solution "); put(i,1); put_line(" ...");
        declare
          lvi : Standard_Dense_Series_Vectors.Link_to_Vector := v(i);
          sol : Standard_Complex_Vectors.Vector(lvi'range);
        begin
          for k in sol'range loop
            sol(k) := lvi(k).cff(0);
          end loop;
          Standard_Complex_Solutions_io.put_vector(sol);
        end;
      end if;
      Run_SVD_Newton(nbrit,p,v(i).all,rcond,verbose);
      if verbose then
        if pause then
          put("Continue ? (y/n) "); Ask_Yes_or_No(ans);
          exit when (ans /= 'y');
        end if;
      end if;
    end loop;
  end Run_SVD_Newton;

  procedure Run_SVD_Newton
              ( nbrit : in integer32;
                p : in DoblDobl_Series_Poly_Systems.Poly_Sys;
                v : in DoblDobl_Dense_Series_VecVecs.VecVec;
                verbose : in boolean := false;
                pause : in boolean := false ) is

    ans : character;
    rcond : double_double;

  begin
    for i in v'range loop
      if verbose then
        put("Running on solution "); put(i,1); put_line(" ...");
        declare
          lvi : DoblDobl_Dense_Series_Vectors.Link_to_Vector := v(i);
          sol : DoblDobl_Complex_Vectors.Vector(lvi'range);
        begin
          for k in sol'range loop
            sol(k) := lvi(k).cff(0);
          end loop;
          DoblDobl_Complex_Solutions_io.put_vector(sol);
        end;
      end if;
      Run_SVD_Newton(nbrit,p,v(i).all,rcond,verbose);
      if verbose then
        if pause then
          put("Continue ? (y/n) "); Ask_Yes_or_No(ans);
          exit when (ans /= 'y');
        end if;
      end if;
    end loop;
  end Run_SVD_Newton;

  procedure Run_SVD_Newton
              ( nbrit : in integer32;
                p : in QuadDobl_Series_Poly_Systems.Poly_Sys;
                v : in QuadDobl_Dense_Series_VecVecs.VecVec;
                verbose : in boolean := false;
                pause : in boolean := false ) is

    ans : character;
    rcond : quad_double;

  begin
    for i in v'range loop
      if verbose then
        put("Running on solution "); put(i,1); put_line(" ...");
        declare
          lvi : QuadDobl_Dense_Series_Vectors.Link_to_Vector := v(i);
          sol : QuadDobl_Complex_Vectors.Vector(lvi'range);
        begin
          for k in sol'range loop
            sol(k) := lvi(k).cff(0);
          end loop;
          QuadDobl_Complex_Solutions_io.put_vector(sol);
        end;
      end if;
      Run_SVD_Newton(nbrit,p,v(i).all,rcond,verbose);
      if verbose then
        if pause then
          put("Continue ? (y/n) "); Ask_Yes_or_No(ans);
          exit when (ans /= 'y');
        end if;
      end if;
    end loop;
  end Run_SVD_Newton;

  procedure Run_SVD_Newton
              ( file : in file_type; nbrit : in integer32;
                p : in Standard_Series_Poly_Systems.Poly_Sys;
                v : in Standard_Dense_Series_VecVecs.VecVec;
                verbose : in boolean := false ) is

    rcond : double_float;

  begin
    for i in v'range loop
      if verbose then
        put(file,"Running on solution ");
        put(file,i,1); put_line(file," ...");
        declare
          lvi : Standard_Dense_Series_Vectors.Link_to_Vector := v(i);
          sol : Standard_Complex_Vectors.Vector(lvi'range);
        begin
          for k in sol'range loop
            sol(k) := lvi(k).cff(0);
          end loop;
          Standard_Complex_Solutions_io.put_vector(file,sol);
        end;
      end if;
      Run_SVD_Newton(file,nbrit,p,v(i).all,rcond,verbose);
    end loop;
  end Run_SVD_Newton;

  procedure Run_SVD_Newton
              ( file : in file_type; nbrit : in integer32;
                p : in DoblDobl_Series_Poly_Systems.Poly_Sys;
                v : in DoblDobl_Dense_Series_VecVecs.VecVec;
                verbose : in boolean := false ) is

    rcond : double_double;

  begin
    for i in v'range loop
      if verbose then
        put(file,"Running on solution ");
        put(file,i,1); put_line(file," ...");
        declare
          lvi : DoblDobl_Dense_Series_Vectors.Link_to_Vector := v(i);
          sol : DoblDobl_Complex_Vectors.Vector(lvi'range);
        begin
          for k in sol'range loop
            sol(k) := lvi(k).cff(0);
          end loop;
          DoblDobl_Complex_Solutions_io.put_vector(file,sol);
        end;
      end if;
      Run_SVD_Newton(file,nbrit,p,v(i).all,rcond,verbose);
    end loop;
  end Run_SVD_Newton;

  procedure Run_SVD_Newton
              ( file : in file_type; nbrit : in integer32;
                p : in QuadDobl_Series_Poly_Systems.Poly_Sys;
                v : in QuadDobl_Dense_Series_VecVecs.VecVec;
                verbose : in boolean := false ) is

    rcond : quad_double;

  begin
    for i in v'range loop
      if verbose then
        put(file,"Running on solution ");
        put(file,i,1); put_line(file," ...");
        declare
          lvi : QuadDobl_Dense_Series_Vectors.Link_to_Vector := v(i);
          sol : QuadDobl_Complex_Vectors.Vector(lvi'range);
        begin
          for k in sol'range loop
            sol(k) := lvi(k).cff(0);
          end loop;
          QuadDobl_Complex_Solutions_io.put_vector(file,sol);
        end;
      end if;
      Run_SVD_Newton(file,nbrit,p,v(i).all,rcond,verbose);
    end loop;
  end Run_SVD_Newton;

end Power_Series_Methods;
