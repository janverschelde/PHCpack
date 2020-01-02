with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;         use Standard_Floating_Numbers;
with Quad_Double_Numbers_io;            use Quad_Double_Numbers_io;
with QuadDobl_Complex_Numbers_io;       use QuadDobl_Complex_Numbers_io;
with QuadDobl_Complex_Series_io;
with QuadDobl_CSeries_Vector_Norms;
with QuadDobl_Complex_Series_Matrices;
with QuadDobl_Complex_Vector_Series;
with QuadDobl_Complex_Matrix_Series;
with QuadDobl_Series_Matrix_Solvers;    use QuadDobl_Series_Matrix_Solvers;
with QuadDobl_CSeries_Polynomials;
with Complex_Series_and_Polynomials;
with Standard_Newton_Matrix_Series; -- for the Double_Degree_with_Threshold

package body QuadDobl_Newton_Matrix_Series is

-- ONE NEWTON STEP WITH LU WITHOUT CONDITION NUMBER ESTIMATE :

  procedure LU_Newton_Step
              ( p : in QuadDobl_CSeries_Poly_Systems.Poly_Sys;
                jp : in QuadDobl_CSeries_Jaco_Matrices.Jaco_Mat;
                degree : in integer32;
                x : in out QuadDobl_Complex_Series_Vectors.Vector;
                info : out integer32; vrblvl : in integer32 := 0 ) is

    dx : QuadDobl_Complex_Series_Vectors.Vector(x'range);
    px : QuadDobl_Complex_Series_Vectors.Vector(p'range);
    jm : QuadDobl_Complex_Series_Matrices.Matrix(p'range,x'range);
    xp,xd : QuadDobl_Complex_Vector_Series.Vector(degree);
    mj : QuadDobl_Complex_Matrix_Series.Matrix(degree);

  begin
    if vrblvl > 0 then
      put_line("-> in quaddobl_newton_matrix_series.LU_Newton_Step 1 ...");
    end if;
    px := QuadDobl_CSeries_Poly_SysFun.Eval(p,x);
    QuadDobl_Complex_Series_Vectors.Min(px);
    Complex_Series_and_Polynomials.Set_degree(px,degree);
    jm := QuadDobl_CSeries_Jaco_Matrices.Eval(jp,x);
    Complex_Series_and_Polynomials.Set_degree(jm,degree);
    mj := QuadDobl_Complex_Matrix_Series.Create(jm);
    xp := QuadDobl_Complex_Vector_Series.Create(px);
    Solve_by_lufac(mj,xp,info,xd);
    if info = 0 then
      dx := QuadDobl_Complex_Vector_Series.Create(xd);
      QuadDobl_Complex_Series_Vectors.Add(x,dx);
      QuadDobl_Complex_Series_Vectors.Clear(dx);
    end if;
    QuadDobl_Complex_Series_Vectors.Clear(px);
    QuadDobl_Complex_Series_Matrices.Clear(jm);
    QuadDobl_Complex_Matrix_Series.Clear(mj);
    QuadDobl_Complex_Vector_Series.Clear(xp);
    QuadDobl_Complex_Vector_Series.Clear(xd);
  end LU_Newton_Step;

  procedure LU_Newton_Step
              ( p : in QuadDobl_CSeries_Poly_Systems.Poly_Sys;
                degree : in integer32;
                x : in out QuadDobl_Complex_Series_Vectors.Vector;
                info : out integer32; vrblvl : in integer32 := 0 ) is

    jp : QuadDobl_CSeries_Jaco_Matrices.Jaco_Mat(p'range,x'range)
       := QuadDobl_CSeries_Jaco_Matrices.Create(p);

  begin
    if vrblvl > 0 then
      put_line("-> in quaddobl_newton_matrix_series.LU_Newton_Step 2 ...");
    end if;
    LU_Newton_Step(p,jp,degree,x,info,vrblvl-1);
    QuadDobl_CSeries_Jaco_Matrices.Clear(jp);
  end LU_Newton_Step;

  procedure LU_Newton_Step
              ( file : in file_type;
                p : in QuadDobl_CSeries_Poly_Systems.Poly_Sys;
                jp : in QuadDobl_CSeries_Jaco_Matrices.Jaco_Mat;
                degree : in integer32;
                x : in out QuadDobl_Complex_Series_Vectors.Vector;
                info : out integer32; vrblvl : in integer32 := 0 ) is

    dx : QuadDobl_Complex_Series_Vectors.Vector(x'range);
    px : QuadDobl_Complex_Series_Vectors.Vector(p'range);
    jm : QuadDobl_Complex_Series_Matrices.Matrix(p'range,x'range);
    xp,xd : QuadDobl_Complex_Vector_Series.Vector(degree);
    mj : QuadDobl_Complex_Matrix_Series.Matrix(degree);
    nrm : quad_double;

  begin
    if vrblvl > 0 then
      put_line("-> in quaddobl_newton_matrix_series.LU_Newton_Step 3 ...");
    end if;
    px := QuadDobl_CSeries_Poly_SysFun.Eval(p,x);
    put_line(file,"The evaluated series :");
    for i in px'range loop
      QuadDobl_Complex_Series_io.put(file,px(i)); new_line(file);
    end loop;
    nrm := QuadDobl_CSeries_Vector_Norms.Max_Norm(px);
    put(file,"The max norm of the evaluation : ");
    put(file,nrm,3); new_line(file);
    QuadDobl_Complex_Series_Vectors.Min(px);
    Complex_Series_and_Polynomials.Set_degree(px,degree);
    jm := QuadDobl_CSeries_Jaco_Matrices.Eval(jp,x);
    Complex_Series_and_Polynomials.Set_degree(jm,degree);
    mj := QuadDobl_Complex_Matrix_Series.Create(jm);
    xp := QuadDobl_Complex_Vector_Series.Create(px);
    Solve_by_lufac(mj,xp,info,xd);
    if info /= 0 then
      put(file,"LUfac info : "); put(file,info,1); new_line(file);
    else
      dx := QuadDobl_Complex_Vector_Series.Create(xd);
      put_line(file,"The update to the series :");
      for i in dx'range loop
        QuadDobl_Complex_Series_io.put(file,dx(i)); new_line(file);
      end loop;
      nrm := QuadDobl_CSeries_Vector_Norms.Max_Norm(dx);
      put(file,"The max norm of the update : ");
      put(file,nrm,3); new_line(file);
      QuadDobl_Complex_Series_Vectors.Add(x,dx);
      QuadDobl_Complex_Series_Vectors.Clear(dx);
    end if;
    QuadDobl_Complex_Series_Vectors.Clear(px);
    QuadDobl_Complex_Series_Matrices.Clear(jm);
    QuadDobl_Complex_Matrix_Series.Clear(mj);
    QuadDobl_Complex_Vector_Series.Clear(xp);
    QuadDobl_Complex_Vector_Series.Clear(xd);
  end LU_Newton_Step;

  procedure LU_Newton_Step
              ( file : in file_type;
                p : in QuadDobl_CSeries_Poly_Systems.Poly_Sys;
                degree : in integer32;
                x : in out QuadDobl_Complex_Series_Vectors.Vector;
                info : out integer32; vrblvl : in integer32 := 0 ) is

    jp : QuadDobl_CSeries_Jaco_Matrices.Jaco_Mat(p'range,x'range)
       := QuadDobl_CSeries_Jaco_Matrices.Create(p);

  begin
    if vrblvl > 0 then
      put_line("-> in quaddobl_newton_matrix_series.LU_Newton_Step 4 ...");
    end if;
    LU_Newton_Step(file,p,jp,degree,x,info,vrblvl-1);
    QuadDobl_CSeries_Jaco_Matrices.Clear(jp);
  end LU_Newton_Step;

-- ONE NEWTON STEP WITH LU ON COEFFICIENT-PARAMETER HOMOTOPIES :

  procedure LU_Newton_Step
              ( f : in QuadDobl_CSeries_Poly_SysFun.Eval_Coeff_Poly_Sys;
                c : in QuadDobl_Complex_Series_VecVecs.VecVec;
                ejm : in QuadDobl_CSeries_Jaco_Matrices.Eval_Coeff_Jaco_Mat;
                mlt : in QuadDobl_CSeries_Jaco_Matrices.Mult_Factors;
                degree : in integer32;
                x : in out QuadDobl_Complex_Series_Vectors.Vector;
                info : out integer32; vrblvl : in integer32 := 0 ) is

    dx : QuadDobl_Complex_Series_Vectors.Vector(x'range);
    px : QuadDobl_Complex_Series_Vectors.Vector(f'range);
    jm : QuadDobl_Complex_Series_Matrices.Matrix(f'range,x'range);
    xp,xd : QuadDobl_Complex_Vector_Series.Vector(degree);
    mj : QuadDobl_Complex_Matrix_Series.Matrix(degree);

  begin
    if vrblvl > 0 then
      put_line("-> in quaddobl_newton_matrix_series.LU_Newton_Step 5 ...");
    end if;
    px := QuadDobl_CSeries_Poly_SysFun.Eval(f,c,x);
    QuadDobl_Complex_Series_Vectors.Min(px);
    Complex_Series_and_Polynomials.Set_degree(px,degree);
    jm := QuadDobl_CSeries_Jaco_Matrices.Eval(ejm,mlt,c,x);
    Complex_Series_and_Polynomials.Set_degree(jm,degree);
    mj := QuadDobl_Complex_Matrix_Series.Create(jm);
    xp := QuadDobl_Complex_Vector_Series.Create(px);
    Solve_by_lufac(mj,xp,info,xd);
    if info = 0 then
      dx := QuadDobl_Complex_Vector_Series.Create(xd);
      QuadDobl_Complex_Series_Vectors.Add(x,dx);
      QuadDobl_Complex_Series_Vectors.Clear(dx);
    end if;
    QuadDobl_Complex_Series_Vectors.Clear(px);
    QuadDobl_Complex_Series_Matrices.Clear(jm);
    QuadDobl_Complex_Matrix_Series.Clear(mj);
    QuadDobl_Complex_Vector_Series.Clear(xp);
    QuadDobl_Complex_Vector_Series.Clear(xd);
  end LU_Newton_Step;

  procedure LU_Newton_Step
              ( file : in file_type;
                f : in QuadDobl_CSeries_Poly_SysFun.Eval_Coeff_Poly_Sys;
                c : in QuadDobl_Complex_Series_VecVecs.VecVec;
                ejm : in QuadDobl_CSeries_Jaco_Matrices.Eval_Coeff_Jaco_Mat;
                mlt : in QuadDobl_CSeries_Jaco_Matrices.Mult_Factors;
                degree : in integer32;
                x : in out QuadDobl_Complex_Series_Vectors.Vector;
                info : out integer32; vrblvl : in integer32 := 0 ) is

    dx : QuadDobl_Complex_Series_Vectors.Vector(x'range);
    px : QuadDobl_Complex_Series_Vectors.Vector(f'range);
    jm : QuadDobl_Complex_Series_Matrices.Matrix(f'range,x'range);
    xp,xd : QuadDobl_Complex_Vector_Series.Vector(degree);
    mj : QuadDobl_Complex_Matrix_Series.Matrix(degree);
    nrm : quad_double;

  begin
    if vrblvl > 0 then
      put_line("-> in quaddobl_newton_matrix_series.LU_Newton_Step 6 ...");
    end if;
    px := QuadDobl_CSeries_Poly_SysFun.Eval(f,c,x);
    put_line(file,"The evaluated series :");
    for i in px'range loop
      QuadDobl_Complex_Series_io.put(file,px(i)); new_line(file);
    end loop;
    nrm := QuadDobl_CSeries_Vector_Norms.Max_Norm(px);
    put(file,"The max norm of the evaluation : ");
    put(file,nrm,3); new_line(file);
    QuadDobl_Complex_Series_Vectors.Min(px);
    Complex_Series_and_Polynomials.Set_degree(px,degree);
    jm := QuadDobl_CSeries_Jaco_Matrices.Eval(ejm,mlt,c,x);
    Complex_Series_and_Polynomials.Set_degree(jm,degree);
    mj := QuadDobl_Complex_Matrix_Series.Create(jm);
    xp := QuadDobl_Complex_Vector_Series.Create(px);
    Solve_by_lufac(mj,xp,info,xd);
    if info /= 0 then
      put(file,"LUfac info : "); put(file,info,1); new_line(file);
    else
      dx := QuadDobl_Complex_Vector_Series.Create(xd);
      put_line(file,"The update to the series :");
      for i in dx'range loop
        QuadDobl_Complex_Series_io.put(file,dx(i)); new_line(file);
      end loop;
      nrm := QuadDobl_CSeries_Vector_Norms.Max_Norm(px);
      put(file,"The max norm of the evaluation : ");
      put(file,nrm,3); new_line(file);
      QuadDobl_Complex_Series_Vectors.Add(x,dx);
      QuadDobl_Complex_Series_Vectors.Clear(dx);
    end if;
    QuadDobl_Complex_Series_Vectors.Clear(px);
    QuadDobl_Complex_Series_Matrices.Clear(jm);
    QuadDobl_Complex_Matrix_Series.Clear(mj);
    QuadDobl_Complex_Vector_Series.Clear(xp);
    QuadDobl_Complex_Vector_Series.Clear(xd);
  end LU_Newton_Step;

-- ONE NEWTON STEP WITH LU WITH CONDITION NUMBER ESTIMATE :

  procedure LU_Newton_Step
              ( p : in QuadDobl_CSeries_Poly_Systems.Poly_Sys;
                jp : in QuadDobl_CSeries_Jaco_Matrices.Jaco_Mat;
                degree : in integer32;
                x : in out QuadDobl_Complex_Series_Vectors.Vector;
                rcond : out quad_double; vrblvl : in integer32 := 0 ) is

    dx : QuadDobl_Complex_Series_Vectors.Vector(x'range);
    px : QuadDobl_Complex_Series_Vectors.Vector(p'range);
    jm : QuadDobl_Complex_Series_Matrices.Matrix(p'range,x'range);
    xp,xd : QuadDobl_Complex_Vector_Series.Vector(degree);
    mj : QuadDobl_Complex_Matrix_Series.Matrix(degree);
    one : constant quad_double := create(1.0);

  begin
    if vrblvl > 0 then
      put_line("-> in quaddobl_newton_matrix_series.LU_Newton_Step 7 ...");
    end if;
    px := QuadDobl_CSeries_Poly_SysFun.Eval(p,x);
    QuadDobl_Complex_Series_Vectors.Min(px);
    Complex_Series_and_Polynomials.Set_degree(px,degree);
    jm := QuadDobl_CSeries_Jaco_Matrices.Eval(jp,x);
    Complex_Series_and_Polynomials.Set_degree(jm,degree);
    mj := QuadDobl_Complex_Matrix_Series.Create(jm);
    xp := QuadDobl_Complex_Vector_Series.Create(px);
    Solve_by_lufco(mj,xp,rcond,xd);
    if one + rcond /= one then
      dx := QuadDobl_Complex_Vector_Series.Create(xd);
      QuadDobl_Complex_Series_Vectors.Add(x,dx);
      QuadDobl_Complex_Series_Vectors.Clear(dx);
    end if;
    QuadDobl_Complex_Series_Vectors.Clear(px);
    QuadDobl_Complex_Series_Matrices.Clear(jm);
    QuadDobl_Complex_Matrix_Series.Clear(mj);
    QuadDobl_Complex_Vector_Series.Clear(xp);
    QuadDobl_Complex_Vector_Series.Clear(xd);
  end LU_Newton_Step;

  procedure LU_Newton_Step
              ( p : in QuadDobl_CSeries_Poly_Systems.Poly_Sys;
                degree : in integer32;
                x : in out QuadDobl_Complex_Series_Vectors.Vector;
                rcond : out quad_double; vrblvl : in integer32 := 0 ) is

    jp : QuadDobl_CSeries_Jaco_Matrices.Jaco_Mat(p'range,x'range)
       := QuadDobl_CSeries_Jaco_Matrices.Create(p);

  begin
    if vrblvl > 0 then
      put_line("-> in quaddobl_newton_matrix_series.LU_Newton_Step 8 ...");
    end if;
    LU_Newton_Step(p,jp,degree,x,rcond,vrblvl-1);
    QuadDobl_CSeries_Jaco_Matrices.Clear(jp);
  end LU_Newton_Step;

  procedure LU_Newton_Step
              ( file : in file_type;
                p : in QuadDobl_CSeries_Poly_Systems.Poly_Sys;
                jp : in QuadDobl_CSeries_Jaco_Matrices.Jaco_Mat;
                degree : in integer32;
                x : in out QuadDobl_Complex_Series_Vectors.Vector;
                rcond : out quad_double; vrblvl : in integer32 := 0 ) is

    dx : QuadDobl_Complex_Series_Vectors.Vector(x'range);
    px : QuadDobl_Complex_Series_Vectors.Vector(p'range);
    jm : QuadDobl_Complex_Series_Matrices.Matrix(p'range,x'range);
    xp,xd : QuadDobl_Complex_Vector_Series.Vector(degree);
    mj : QuadDobl_Complex_Matrix_Series.Matrix(degree);
    one : constant quad_double := create(1.0);
    nrm : quad_double;

  begin
    if vrblvl > 0 then
      put_line("-> in quaddobl_newton_matrix_series.LU_Newton_Step 9 ...");
    end if;
    px := QuadDobl_CSeries_Poly_SysFun.Eval(p,x);
    put_line(file,"The evaluated series :");
    for i in px'range loop
      QuadDobl_Complex_Series_io.put(file,px(i)); new_line(file);
    end loop;
    nrm := QuadDobl_CSeries_Vector_Norms.Max_Norm(px);
    put(file,"The max norm of the evaluation : ");
    put(file,nrm,3); new_line(file);
    QuadDobl_Complex_Series_Vectors.Min(px);
    Complex_Series_and_Polynomials.Set_degree(px,degree);
    jm := QuadDobl_CSeries_Jaco_Matrices.Eval(jp,x);
    Complex_Series_and_Polynomials.Set_degree(jm,degree);
    mj := QuadDobl_Complex_Matrix_Series.Create(jm);
    xp := QuadDobl_Complex_Vector_Series.Create(px);
    Solve_by_lufco(mj,xp,rcond,xd);
    put(file,"LUfco rcond : "); put(file,rcond,3); new_line(file);
    if one + rcond /= one then
      dx := QuadDobl_Complex_Vector_Series.Create(xd);
      put_line(file,"The update to the series :");
      for i in dx'range loop
        QuadDobl_Complex_Series_io.put(file,dx(i)); new_line(file);
      end loop;
      nrm := QuadDobl_CSeries_Vector_Norms.Max_Norm(dx);
      put(file,"The max norm of the update : ");
      put(file,nrm,3); new_line(file);
      QuadDobl_Complex_Series_Vectors.Add(x,dx);
      QuadDobl_Complex_Series_Vectors.Clear(dx);
    end if;
    QuadDobl_Complex_Series_Vectors.Clear(px);
    QuadDobl_Complex_Series_Matrices.Clear(jm);
    QuadDobl_Complex_Matrix_Series.Clear(mj);
    QuadDobl_Complex_Vector_Series.Clear(xp);
    QuadDobl_Complex_Vector_Series.Clear(xd);
  end LU_Newton_Step;

  procedure LU_Newton_Step
              ( file : in file_type;
                p : in QuadDobl_CSeries_Poly_Systems.Poly_Sys;
                degree : in integer32;
                x : in out QuadDobl_Complex_Series_Vectors.Vector;
                rcond : out quad_double; vrblvl : in integer32 := 0 ) is

    jp : QuadDobl_CSeries_Jaco_Matrices.Jaco_Mat(p'range,x'range)
       := QuadDobl_CSeries_Jaco_Matrices.Create(p);

  begin
    if vrblvl > 0 then
      put_line("-> in quaddobl_newton_matrix_series.LU_Newton_Step 10 ...");
    end if;
    LU_Newton_Step(file,p,jp,degree,x,rcond,vrblvl-1);
    QuadDobl_CSeries_Jaco_Matrices.Clear(jp);
  end LU_Newton_Step;

-- ONE NEWTON STEP WITH QR :

  procedure QR_Newton_Step
              ( p : in QuadDobl_CSeries_Poly_Systems.Poly_Sys;
                jp : in QuadDobl_CSeries_Jaco_Matrices.Jaco_Mat;
                degree : in integer32;
                x : in out QuadDobl_Complex_Series_Vectors.Vector;
                info : out integer32; vrblvl : in integer32 := 0 ) is

    wp : QuadDobl_CSeries_Poly_Systems.Poly_Sys(p'range);
    dx : QuadDobl_Complex_Series_Vectors.Vector(x'range);
    px : QuadDobl_Complex_Series_Vectors.Vector(p'range);
    jm : QuadDobl_Complex_Series_Matrices.Matrix(p'range,x'range);
    mj : QuadDobl_Complex_Matrix_Series.Matrix(degree);
    xp,xd : QuadDobl_Complex_Vector_Series.Vector(degree);
    nrm : quad_double;
    tol : constant double_float := 1.0E-48;
    wjp : QuadDobl_CSeries_Jaco_Matrices.Jaco_Mat(jp'range(1),jp'range(2));

  begin
    if vrblvl > 0 then
      put_line("-> in quaddobl_newton_matrix_series.QR_Newton_Step 1 ...");
    end if;
    QuadDobl_CSeries_Poly_Systems.Copy(p,wp);
    Complex_Series_and_Polynomials.Set_degree(wp,degree);
    px := QuadDobl_CSeries_Poly_SysFun.Eval(wp,x);
    nrm := QuadDobl_CSeries_Vector_Norms.Max_Norm(px);
    if nrm > tol then
      for i in jp'range(1) loop
        for j in jp'range(2) loop
          QuadDobl_CSeries_Polynomials.Copy(jp(i,j),wjp(i,j));
        end loop;
      end loop;
      Complex_Series_and_Polynomials.Set_degree(wjp,degree);
      QuadDobl_Complex_Series_Vectors.Min(px);
      Complex_Series_and_Polynomials.Set_degree(px,degree);
      jm := QuadDobl_CSeries_Jaco_Matrices.Eval(wjp,x);
      Complex_Series_and_Polynomials.Set_degree(jm,degree);
      mj := QuadDobl_Complex_Matrix_Series.Create(jm);
      xp := QuadDobl_Complex_Vector_Series.Create(px);
      Solve_by_QRLS(mj,xp,info,xd);
      if info = 0 then
        dx := QuadDobl_Complex_Vector_Series.Create(xd);
        QuadDobl_Complex_Series_Vectors.Add(x,dx);
        QuadDobl_Complex_Series_Vectors.Clear(dx);
      end if;
      QuadDobl_CSeries_Jaco_Matrices.Clear(wjp);
    end if;
    QuadDobl_Complex_Series_Vectors.Clear(px);
    QuadDobl_Complex_Series_Matrices.Clear(jm);
    QuadDobl_CSeries_Poly_Systems.Clear(wp);
    QuadDobl_Complex_Matrix_Series.Clear(mj);
    QuadDobl_Complex_Vector_Series.Clear(xp);
    QuadDobl_Complex_Vector_Series.Clear(xd);
  end QR_Newton_Step;

  procedure QR_Newton_Step
              ( p : in QuadDobl_CSeries_Poly_Systems.Poly_Sys;
                degree : in integer32;
                x : in out QuadDobl_Complex_Series_Vectors.Vector;
                info : out integer32; vrblvl : in integer32 := 0 ) is

    jp : QuadDobl_CSeries_Jaco_Matrices.Jaco_Mat(p'range,x'range)
       := QuadDobl_CSeries_Jaco_Matrices.Create(p);

  begin
    if vrblvl > 0 then
      put_line("-> in quaddobl_newton_matrix_series.QR_Newton_Step 2 ...");
    end if;
    QR_Newton_Step(p,jp,degree,x,info,vrblvl-1);
    QuadDobl_CSeries_Jaco_Matrices.Clear(jp);
  end QR_Newton_Step;

  procedure QR_Newton_Step
              ( file : in file_type;
                p : in QuadDobl_CSeries_Poly_Systems.Poly_Sys;
                jp : in QuadDobl_CSeries_Jaco_Matrices.Jaco_Mat;
                degree : in integer32;
                x : in out QuadDobl_Complex_Series_Vectors.Vector;
                info : out integer32; vrblvl : in integer32 := 0 ) is

    wp : QuadDobl_CSeries_Poly_Systems.Poly_Sys(p'range);
    dx : QuadDobl_Complex_Series_Vectors.Vector(x'range);
    px : QuadDobl_Complex_Series_Vectors.Vector(p'range);
    jm : QuadDobl_Complex_Series_Matrices.Matrix(p'range,x'range);
    mj : QuadDobl_Complex_Matrix_Series.Matrix(degree);
    xp,xd : QuadDobl_Complex_Vector_Series.Vector(degree);
    nrm : quad_double;
    tol : constant double_float := 1.0E-48;
    wjp : QuadDobl_CSeries_Jaco_Matrices.Jaco_Mat(jm'range(1),jm'range(2));

  begin
    if vrblvl > 0 then
      put_line("-> in quaddobl_newton_matrix_series.QR_Newton_Step 3 ...");
    end if;
    QuadDobl_CSeries_Poly_Systems.Copy(p,wp);
    Complex_Series_and_Polynomials.Set_degree(wp,degree);
    px := QuadDobl_CSeries_Poly_SysFun.Eval(wp,x);
    put_line(file,"The evaluated series :");
    for i in px'range loop
      QuadDobl_Complex_Series_io.put(file,px(i)); new_line(file);
    end loop;
    nrm := QuadDobl_CSeries_Vector_Norms.Max_Norm(px);
    put(file,"The max norm of the evaluation : ");
    put(file,nrm,3); new_line(file);
    if nrm > tol then
      for i in jp'range(1) loop
        for j in jp'range(2) loop
          QuadDobl_CSeries_Polynomials.Copy(jp(i,j),wjp(i,j));
        end loop;
      end loop;
      Complex_Series_and_Polynomials.Set_degree(wjp,degree);
      QuadDobl_Complex_Series_Vectors.Min(px);
      Complex_Series_and_Polynomials.Set_degree(px,degree);
      jm := QuadDobl_CSeries_Jaco_Matrices.Eval(jp,x);
      Complex_Series_and_Polynomials.Set_degree(jm,degree);
      mj := QuadDobl_Complex_Matrix_Series.Create(jm);
      xp := QuadDobl_Complex_Vector_Series.Create(px);
      Solve_by_QRLS(mj,xp,info,xd);
      if info /= 0 then
        put(file,"QRLS info : "); put(file,info,1); new_line(file);
      else
        dx := QuadDobl_Complex_Vector_Series.Create(xd);
        put_line(file,"The update to the series :");
        for i in dx'range loop
          QuadDobl_Complex_Series_io.put(file,dx(i)); new_line(file);
        end loop;
        nrm := QuadDobl_CSeries_Vector_Norms.Max_Norm(dx);
        put(file,"The max norm of the update : ");
        put(file,nrm,3); new_line(file);
        QuadDobl_Complex_Series_Vectors.Add(x,dx);
        QuadDobl_Complex_Series_Vectors.Clear(dx);
      end if;
      QuadDobl_CSeries_Jaco_Matrices.Clear(wjp);
    end if;
    QuadDobl_Complex_Series_Vectors.Clear(px);
    QuadDobl_Complex_Series_Matrices.Clear(jm);
    QuadDobl_CSeries_Poly_Systems.Clear(wp);
    QuadDobl_Complex_Matrix_Series.Clear(mj);
    QuadDobl_Complex_Vector_Series.Clear(xp);
    QuadDobl_Complex_Vector_Series.Clear(xd);
  end QR_Newton_Step;

  procedure QR_Newton_Step
              ( file : in file_type;
                p : in QuadDobl_CSeries_Poly_Systems.Poly_Sys;
                degree : in integer32;
                x : in out QuadDobl_Complex_Series_Vectors.Vector;
                info : out integer32; vrblvl : in integer32 := 0 ) is

    jp : QuadDobl_CSeries_Jaco_Matrices.Jaco_Mat(p'range,x'range)
       := QuadDobl_CSeries_Jaco_Matrices.Create(p);

  begin
    if vrblvl > 0 then
      put_line("-> in quaddobl_newton_matrix_series.QR_Newton_Step 4 ...");
    end if;
    QR_Newton_Step(file,p,jp,degree,x,info,vrblvl-1);
    QuadDobl_CSeries_Jaco_Matrices.Clear(jp);
  end QR_Newton_Step;

-- ONE NEWTON STEP WITH QR ON COEFFICIENT-PARAMETER HOMOTOPIES :

  procedure QR_Newton_Step
              ( f : in QuadDobl_CSeries_Poly_SysFun.Eval_Coeff_Poly_Sys;
                c : in QuadDobl_Complex_Series_VecVecs.VecVec;
                ejm : in QuadDobl_CSeries_Jaco_Matrices.Eval_Coeff_Jaco_Mat;
                mlt : in QuadDobl_CSeries_Jaco_Matrices.Mult_Factors;
                degree : in integer32;
                x : in out QuadDobl_Complex_Series_Vectors.Vector;
                info : out integer32; vrblvl : in integer32 := 0 ) is

    dx : QuadDobl_Complex_Series_Vectors.Vector(x'range);
    px : QuadDobl_Complex_Series_Vectors.Vector(f'range);
    jm : QuadDobl_Complex_Series_Matrices.Matrix(f'range,x'range);
    mj : QuadDobl_Complex_Matrix_Series.Matrix(degree);
    xp,xd : QuadDobl_Complex_Vector_Series.Vector(degree);
    nrm : quad_double;
    tol : constant double_float := 1.0E-48;

  begin
    if vrblvl > 0 then
      put_line("-> in quaddobl_newton_matrix_series.QR_Newton_Step 5 ...");
    end if;
    px := QuadDobl_CSeries_Poly_SysFun.Eval(f,c,x);
    nrm := QuadDobl_CSeries_Vector_Norms.Max_Norm(px);
    if nrm > tol then
      QuadDobl_Complex_Series_Vectors.Min(px);
      Complex_Series_and_Polynomials.Set_degree(px,degree);
      jm := QuadDobl_CSeries_Jaco_Matrices.Eval(ejm,mlt,c,x);
      Complex_Series_and_Polynomials.Set_degree(jm,degree);
      mj := QuadDobl_Complex_Matrix_Series.Create(jm);
      xp := QuadDobl_Complex_Vector_Series.Create(px);
      Solve_by_QRLS(mj,xp,info,xd);
      if info = 0 then
        dx := QuadDobl_Complex_Vector_Series.Create(xd);
        QuadDobl_Complex_Series_Vectors.Add(x,dx);
        QuadDobl_Complex_Series_Vectors.Clear(dx);
      end if;
    end if;
    QuadDobl_Complex_Series_Vectors.Clear(px);
    QuadDobl_Complex_Series_Matrices.Clear(jm);
    QuadDobl_Complex_Matrix_Series.Clear(mj);
    QuadDobl_Complex_Vector_Series.Clear(xp);
    QuadDobl_Complex_Vector_Series.Clear(xd);
  end QR_Newton_Step;

  procedure QR_Newton_Step
              ( file : in file_type;
                f : in QuadDobl_CSeries_Poly_SysFun.Eval_Coeff_Poly_Sys;
                c : in QuadDobl_Complex_Series_VecVecs.VecVec;
                ejm : in QuadDobl_CSeries_Jaco_Matrices.Eval_Coeff_Jaco_Mat;
                mlt : in QuadDobl_CSeries_Jaco_Matrices.Mult_Factors;
                degree : in integer32;
                x : in out QuadDobl_Complex_Series_Vectors.Vector;
                info : out integer32; vrblvl : in integer32 := 0 ) is

    dx : QuadDobl_Complex_Series_Vectors.Vector(x'range);
    px : QuadDobl_Complex_Series_Vectors.Vector(f'range);
    jm : QuadDobl_Complex_Series_Matrices.Matrix(f'range,x'range);
    mj : QuadDobl_Complex_Matrix_Series.Matrix(degree);
    xp,xd : QuadDobl_Complex_Vector_Series.Vector(degree);
    nrm : quad_double;
    tol : constant double_float := 1.0E-48;

  begin
    if vrblvl > 0 then
      put_line("-> in quaddobl_newton_matrix_series.QR_Newton_Step 6 ...");
    end if;
    px := QuadDobl_CSeries_Poly_SysFun.Eval(f,c,x);
    put_line(file,"The evaluated series :");
    for i in px'range loop
      QuadDobl_Complex_Series_io.put(file,px(i)); new_line(file);
    end loop;
    nrm := QuadDobl_CSeries_Vector_Norms.Max_Norm(px);
    put(file,"The max norm of the evaluation : ");
    put(file,nrm,3); new_line(file);
    if nrm > tol then
      QuadDobl_Complex_Series_Vectors.Min(px);
      Complex_Series_and_Polynomials.Set_degree(px,degree);
      jm := QuadDobl_CSeries_Jaco_Matrices.Eval(ejm,mlt,c,x);
      Complex_Series_and_Polynomials.Set_degree(jm,degree);
      mj := QuadDobl_Complex_Matrix_Series.Create(jm);
      xp := QuadDobl_Complex_Vector_Series.Create(px);
      Solve_by_QRLS(mj,xp,info,xd);
      if info /= 0 then
        put(file,"QRLS info : "); put(file,info,1); new_line(file);
      else
        dx := QuadDobl_Complex_Vector_Series.Create(xd);
        put_line(file,"The update to the series :");
        for i in dx'range loop
          QuadDobl_Complex_Series_io.put(file,dx(i)); new_line(file);
        end loop;
        nrm := QuadDobl_CSeries_Vector_Norms.Max_Norm(dx);
        put(file,"The max norm of the update : ");
        put(file,nrm,3); new_line(file);
        QuadDobl_Complex_Series_Vectors.Add(x,dx);
        QuadDobl_Complex_Series_Vectors.Clear(dx);
      end if;
    end if;
    QuadDobl_Complex_Series_Vectors.Clear(px);
    QuadDobl_Complex_Series_Matrices.Clear(jm);
    QuadDobl_Complex_Matrix_Series.Clear(mj);
    QuadDobl_Complex_Vector_Series.Clear(xp);
    QuadDobl_Complex_Vector_Series.Clear(xd);
  end QR_Newton_Step;

-- ONE NEWTON STEP WITH SVD :

  procedure SVD_Newton_Step
              ( p : in QuadDobl_CSeries_Poly_Systems.Poly_Sys;
                jp : in QuadDobl_CSeries_Jaco_Matrices.Jaco_Mat;
                degree : in integer32;
                x : in out QuadDobl_Complex_Series_Vectors.Vector;
                info : out integer32; rcond : out quad_double;
                vrblvl : in integer32 := 0 ) is

    wp : QuadDobl_CSeries_Poly_Systems.Poly_Sys(p'range);
    dx : QuadDobl_Complex_Series_Vectors.Vector(x'range);
    px : QuadDobl_Complex_Series_Vectors.Vector(p'range);
    jm : QuadDobl_Complex_Series_Matrices.Matrix(p'range,x'range);
    mj : QuadDobl_Complex_Matrix_Series.Matrix(degree);
    xp,xd : QuadDobl_Complex_Vector_Series.Vector(degree);
    nrm : quad_double;
    tol : constant double_float := 1.0E-26;
    wjp : QuadDobl_CSeries_Jaco_Matrices.Jaco_Mat(jp'range(1),jp'range(2));
    one : constant quad_double := create(1.0);

  begin
    if vrblvl > 0 then
      put_line("-> in quaddobl_newton_matrix_series.SVD_Newton_Step 1 ...");
    end if;
    QuadDobl_CSeries_Poly_Systems.Copy(p,wp);
    Complex_Series_and_Polynomials.Set_degree(wp,degree);
    px := QuadDobl_CSeries_Poly_SysFun.Eval(wp,x);
    nrm := QuadDobl_CSeries_Vector_Norms.Max_Norm(px);
    if nrm > tol then
      for i in jp'range(1) loop
        for j in jp'range(2) loop
          QuadDobl_CSeries_Polynomials.Copy(jp(i,j),wjp(i,j));
        end loop;
      end loop;
      Complex_Series_and_Polynomials.Set_degree(wjp,degree);
      QuadDobl_Complex_Series_Vectors.Min(px);
      Complex_Series_and_Polynomials.Set_degree(px,degree);
      jm := QuadDobl_CSeries_Jaco_Matrices.Eval(wjp,x);
      Complex_Series_and_Polynomials.Set_degree(jm,degree);
      mj := QuadDobl_Complex_Matrix_Series.Create(jm);
      xp := QuadDobl_Complex_Vector_Series.Create(px);
      Solve_by_SVD(mj,xp,info,rcond,xd);
      if one + rcond /= one then
        dx := QuadDobl_Complex_Vector_Series.Create(xd);
        QuadDobl_Complex_Series_Vectors.Add(x,dx);
        QuadDobl_Complex_Series_Vectors.Clear(dx);
      end if;
      QuadDobl_CSeries_Jaco_Matrices.Clear(wjp);
    end if;
    QuadDobl_Complex_Series_Vectors.Clear(px);
    QuadDobl_Complex_Series_Matrices.Clear(jm);
    QuadDobl_CSeries_Poly_Systems.Clear(wp);
    QuadDobl_Complex_Matrix_Series.Clear(mj);
    QuadDobl_Complex_Vector_Series.Clear(xp);
    QuadDobl_Complex_Vector_Series.Clear(xd);
  end SVD_Newton_Step;

  procedure SVD_Newton_Step
              ( p : in QuadDobl_CSeries_Poly_Systems.Poly_Sys;
                degree : in integer32;
                x : in out QuadDobl_Complex_Series_Vectors.Vector;
                info : out integer32; rcond : out quad_double;
                vrblvl : in integer32 := 0 ) is

    jp : QuadDobl_CSeries_Jaco_Matrices.Jaco_Mat(p'range,x'range)
       := QuadDobl_CSeries_Jaco_Matrices.Create(p);

  begin
    if vrblvl > 0 then
      put_line("-> in quaddobl_newton_matrix_series.SVD_Newton_Step 2 ...");
    end if;
    SVD_Newton_Step(p,jp,degree,x,info,rcond,vrblvl-1);
    QuadDobl_CSeries_Jaco_Matrices.Clear(jp);
  end SVD_Newton_Step;

  procedure SVD_Newton_Step
              ( file : in file_type;
                p : in QuadDobl_CSeries_Poly_Systems.Poly_Sys;
                jp : in QuadDobl_CSeries_Jaco_Matrices.Jaco_Mat;
                degree : in integer32;
                x : in out QuadDobl_Complex_Series_Vectors.Vector;
                info : out integer32; rcond : out quad_double;
                vrblvl : in integer32 := 0 ) is

    wp : QuadDobl_CSeries_Poly_Systems.Poly_Sys(p'range);
    dx : QuadDobl_Complex_Series_Vectors.Vector(x'range);
    px : QuadDobl_Complex_Series_Vectors.Vector(p'range);
    jm : QuadDobl_Complex_Series_Matrices.Matrix(p'range,x'range);
    mj : QuadDobl_Complex_Matrix_Series.Matrix(degree);
    xp,xd : QuadDobl_Complex_Vector_Series.Vector(degree);
    nrm : quad_double;
    tol : constant double_float := 1.0E-26;
    wjp : QuadDobl_CSeries_Jaco_Matrices.Jaco_Mat(jm'range(1),jm'range(2));
    one : constant quad_double := create(1.0);

  begin
    if vrblvl > 0 then
      put_line("-> in quaddobl_newton_matrix_series.SVD_Newton_Step 3 ...");
    end if;
    QuadDobl_CSeries_Poly_Systems.Copy(p,wp);
    Complex_Series_and_Polynomials.Set_degree(wp,degree);
    px := QuadDobl_CSeries_Poly_SysFun.Eval(wp,x);
    put_line(file,"The evaluated series :");
    for i in px'range loop
      QuadDobl_Complex_Series_io.put(file,px(i)); new_line(file);
    end loop;
    nrm := QuadDobl_CSeries_Vector_Norms.Max_Norm(px);
    put(file,"The max norm of the evaluation : ");
    put(file,nrm,3); new_line(file);
    if nrm > tol then
      for i in jp'range(1) loop
        for j in jp'range(2) loop
          QuadDobl_CSeries_Polynomials.Copy(jp(i,j),wjp(i,j));
        end loop;
      end loop;
      Complex_Series_and_Polynomials.Set_degree(wjp,degree);
      QuadDobl_Complex_Series_Vectors.Min(px);
      Complex_Series_and_Polynomials.Set_degree(px,degree);
      jm := QuadDobl_CSeries_Jaco_Matrices.Eval(jp,x);
      Complex_Series_and_Polynomials.Set_degree(jm,degree);
      mj := QuadDobl_Complex_Matrix_Series.Create(jm);
      xp := QuadDobl_Complex_Vector_Series.Create(px);
      Solve_by_SVD(mj,xp,info,rcond,xd);
      put(file,"SVD rcond : "); put(file,rcond,3); new_line(file);
      if one + rcond /= one then
        dx := QuadDobl_Complex_Vector_Series.Create(xd);
        put_line(file,"The update to the series :");
        for i in dx'range loop
          QuadDobl_Complex_Series_io.put(file,dx(i)); new_line(file);
        end loop;
        nrm := QuadDobl_CSeries_Vector_Norms.Max_Norm(dx);
        put(file,"The max norm of the update : ");
        put(file,nrm,3); new_line(file);
        QuadDobl_Complex_Series_Vectors.Add(x,dx);
        QuadDobl_Complex_Series_Vectors.Clear(dx);
      end if;
      QuadDobl_CSeries_Jaco_Matrices.Clear(wjp);
    end if;
    QuadDobl_Complex_Series_Vectors.Clear(px);
    QuadDobl_Complex_Series_Matrices.Clear(jm);
    QuadDobl_CSeries_Poly_Systems.Clear(wp);
    QuadDobl_Complex_Matrix_Series.Clear(mj);
    QuadDobl_Complex_Vector_Series.Clear(xp);
    QuadDobl_Complex_Vector_Series.Clear(xd);
  end SVD_Newton_Step;

  procedure SVD_Newton_Step
              ( file : in file_type;
                p : in QuadDobl_CSeries_Poly_Systems.Poly_Sys;
                degree : in integer32;
                x : in out QuadDobl_Complex_Series_Vectors.Vector;
                info : out integer32; rcond : out quad_double;
                vrblvl : in integer32 := 0 ) is

    jp : QuadDobl_CSeries_Jaco_Matrices.Jaco_Mat(p'range,x'range)
       := QuadDobl_CSeries_Jaco_Matrices.Create(p);

  begin
    if vrblvl > 0 then
      put_line("-> in quaddobl_newton_matrix_series.SVD_Newton_Step 4 ...");
    end if;
    SVD_Newton_Step(file,p,jp,degree,x,info,rcond,vrblvl-1);
    QuadDobl_CSeries_Jaco_Matrices.Clear(jp);
  end SVD_Newton_Step;

-- ONE NEWTON STEP WITH ECHELON FORM :

  procedure Echelon_Newton_Step
              ( p : in QuadDobl_CSeries_Poly_Systems.Poly_Sys;
                jp : in QuadDobl_CSeries_Jaco_Matrices.Jaco_Mat;
                degree : in integer32;
                x : in out QuadDobl_Complex_Series_Vectors.Vector;
                det : out Complex_Number; vrblvl : in integer32 := 0 ) is

    dx : QuadDobl_Complex_Series_Vectors.Vector(x'range);
    px : QuadDobl_Complex_Series_Vectors.Vector(p'range);
    jm : QuadDobl_Complex_Series_Matrices.Matrix(p'range,x'range);
    xp,xd,xdn : QuadDobl_Complex_Vector_Series.Vector(degree);
    mj : QuadDobl_Complex_Matrix_Series.Matrix(degree);

  begin
    if vrblvl > 0 then
      put_line("-> in quaddobl_newton_matrix_series.Echelon_Newton_Step 1 ...");
    end if;
    px := QuadDobl_CSeries_Poly_SysFun.Eval(p,x);
    QuadDobl_Complex_Series_Vectors.Min(px);
    Complex_Series_and_Polynomials.Set_degree(px,degree);
    jm := QuadDobl_CSeries_Jaco_Matrices.Eval(jp,x);
    Complex_Series_and_Polynomials.Set_degree(jm,degree);
    mj := QuadDobl_Complex_Matrix_Series.Create(jm);
    xp := QuadDobl_Complex_Vector_Series.Create(px);
    Echelon_Solve(mj,xp,det,xd,xdn);
    dx := QuadDobl_Complex_Vector_Series.Create(xd);
    QuadDobl_Complex_Series_Vectors.Add(x,dx);
    QuadDobl_Complex_Matrix_Series.Clear(mj);
    QuadDobl_Complex_Vector_Series.Clear(xp);
    QuadDobl_Complex_Vector_Series.Clear(xd);
  end Echelon_Newton_Step;

  procedure Echelon_Newton_Step
              ( p : in QuadDobl_CSeries_Poly_Systems.Poly_Sys;
                degree : in integer32;
                x : in out QuadDobl_Complex_Series_Vectors.Vector;
                det : out Complex_Number; vrblvl : in integer32 := 0 ) is

    jp : QuadDobl_CSeries_Jaco_Matrices.Jaco_Mat(p'range,x'range)
       := QuadDobl_CSeries_Jaco_Matrices.Create(p);

  begin
    if vrblvl > 0 then
      put_line("-> in quaddobl_newton_matrix_series.Echelon_Newton_Step 2 ...");
    end if;
    Echelon_Newton_Step(p,jp,degree,x,det,vrblvl-1);
    QuadDobl_CSeries_Jaco_Matrices.Clear(jp);
  end Echelon_Newton_Step;

  procedure Echelon_Newton_Step
              ( file : in file_type;
                p : in QuadDobl_CSeries_Poly_Systems.Poly_Sys;
                jp : in QuadDobl_CSeries_Jaco_Matrices.Jaco_Mat;
                degree : in integer32;
                x : in out QuadDobl_Complex_Series_Vectors.Vector;
                det : out Complex_Number; vrblvl : in integer32 := 0 ) is

    dx : QuadDobl_Complex_Series_Vectors.Vector(x'range);
    px : QuadDobl_Complex_Series_Vectors.Vector(p'range);
    jm : QuadDobl_Complex_Series_Matrices.Matrix(p'range,x'range);
    xp,xd,xdn : QuadDobl_Complex_Vector_Series.Vector(degree);
    mj : QuadDobl_Complex_Matrix_Series.Matrix(degree);
    nrm : quad_double;

  begin
    if vrblvl > 0 then
      put_line("-> in quaddobl_newton_matrix_series.Echelon_Newton_Step 3 ...");
    end if;
    px := QuadDobl_CSeries_Poly_SysFun.Eval(p,x);
    put_line(file,"The evaluated series :");
    for i in px'range loop
      QuadDobl_Complex_Series_io.put(file,px(i)); new_line(file);
    end loop;
    nrm := QuadDobl_CSeries_Vector_Norms.Max_Norm(px);
    put(file,"The max norm of the evaluation : ");
    put(file,nrm,3); new_line(file);
    QuadDobl_Complex_Series_Vectors.Min(px);
    Complex_Series_and_Polynomials.Set_degree(px,degree);
    jm := QuadDobl_CSeries_Jaco_Matrices.Eval(jp,x);
    Complex_Series_and_Polynomials.Set_degree(jm,degree);
    mj := QuadDobl_Complex_Matrix_Series.Create(jm);
    xp := QuadDobl_Complex_Vector_Series.Create(px);
    Echelon_Solve(mj,xp,det,xd,xdn);
    put(file,"n.deg : "); put(file,xdn.deg,1); 
    put(file,"  det : "); put(file,det); new_line(file);
    dx := QuadDobl_Complex_Vector_Series.Create(xd);
    QuadDobl_Complex_Series_Vectors.Add(x,dx);
    put_line(file,"The update to the series :");
    for i in dx'range loop
      QuadDobl_Complex_Series_io.put(file,dx(i)); new_line(file);
    end loop;
    nrm := QuadDobl_CSeries_Vector_Norms.Max_Norm(dx);
    put(file,"The max norm of the update : ");
    put(file,nrm,3); new_line(file);
    QuadDobl_Complex_Matrix_Series.Clear(mj);
    QuadDobl_Complex_Vector_Series.Clear(xp);
    QuadDobl_Complex_Vector_Series.Clear(xd);
  end Echelon_Newton_Step;

  procedure Echelon_Newton_Step
              ( file : in file_type;
                p : in QuadDobl_CSeries_Poly_Systems.Poly_Sys;
                degree : in integer32;
                x : in out QuadDobl_Complex_Series_Vectors.Vector;
                det : out Complex_Number; vrblvl : in integer32 := 0 ) is

    jp : QuadDobl_CSeries_Jaco_Matrices.Jaco_Mat(p'range,x'range)
       := QuadDobl_CSeries_Jaco_Matrices.Create(p);

  begin
    if vrblvl > 0 then
      put_line("-> in quaddobl_newton_matrix_series.Echelon_Newton_Step 4 ...");
    end if;
    Echelon_Newton_Step(file,p,jp,degree,x,det,vrblvl-1);
    QuadDobl_CSeries_Jaco_Matrices.Clear(jp);
  end Echelon_Newton_Step;

-- MANY NEWTON STEPS WITH LU WITHOUT CONDITION NUMBER ESTIMATE :

  procedure LU_Newton_Steps
              ( p : in QuadDobl_CSeries_Poly_Systems.Poly_Sys;
                jp : in QuadDobl_CSeries_Jaco_Matrices.Jaco_Mat;
                degree : in out integer32; maxdeg,nbrit : in integer32;
                x : in out QuadDobl_Complex_Series_Vectors.Vector;
                info : out integer32; vrblvl : in integer32 := 0 ) is
  begin
    if vrblvl > 0 then
      put_line("-> in quaddobl_newton_matrix_series.LU_Newton_Steps 1 ...");
    end if;
    for i in 1..nbrit loop
      LU_Newton_Step(p,jp,degree,x,info,vrblvl-1);
      exit when (info /= 0); -- stop if Jacobian matrix is singular
      exit when (i = nbrit); -- do not double degree after last step
      Standard_Newton_Matrix_Series.Double_Degree_with_Threshold(degree,maxdeg);
    end loop;
  end LU_Newton_Steps;

  procedure LU_Newton_Steps
              ( p : in QuadDobl_CSeries_Poly_Systems.Poly_Sys;
                degree : in out integer32; maxdeg,nbrit : in integer32;
                x : in out QuadDobl_Complex_Series_Vectors.Vector;
                info : out integer32; vrblvl : in integer32 := 0 ) is

    jp : QuadDobl_CSeries_Jaco_Matrices.Jaco_Mat(p'range,x'range)
       := QuadDobl_CSeries_Jaco_Matrices.Create(p);

  begin
    if vrblvl > 0 then
      put_line("-> in quaddobl_newton_matrix_series.LU_Newton_Steps 2 ...");
    end if;
    LU_Newton_Steps(p,jp,degree,maxdeg,nbrit,x,info,vrblvl-1);
    QuadDobl_CSeries_Jaco_Matrices.Clear(jp);
  end LU_Newton_Steps;

  procedure LU_Newton_Steps
              ( file : in file_type;
                p : in QuadDobl_CSeries_Poly_Systems.Poly_Sys;
                jp : in QuadDobl_CSeries_Jaco_Matrices.Jaco_Mat;
                degree : in out integer32; maxdeg,nbrit : in integer32;
                x : in out QuadDobl_Complex_Series_Vectors.Vector;
                info : out integer32; vrblvl : in integer32 := 0 ) is
  begin
    if vrblvl > 0 then
      put_line("-> in quaddobl_newton_matrix_series.LU_Newton_Steps 3 ...");
    end if;
    for i in 1..nbrit loop
      put(file,"LU Newton step "); put(file,i,1); put_line(file," :");
      LU_Newton_Step(file,p,jp,degree,x,info,vrblvl-1);
      exit when (info /= 0); -- stop if Jacobian matrix is singular
      exit when (i = nbrit); -- do not double degree after last step
      Standard_Newton_Matrix_Series.Double_Degree_with_Threshold(degree,maxdeg);
    end loop;
  end LU_Newton_Steps;

  procedure LU_Newton_Steps
              ( file : in file_type;
                p : in QuadDobl_CSeries_Poly_Systems.Poly_Sys;
                degree : in out integer32; maxdeg,nbrit : in integer32;
                x : in out QuadDobl_Complex_Series_Vectors.Vector;
                info : out integer32; vrblvl : in integer32 := 0 ) is

    jp : QuadDobl_CSeries_Jaco_Matrices.Jaco_Mat(p'range,x'range)
       := QuadDobl_CSeries_Jaco_Matrices.Create(p);

  begin
    if vrblvl > 0 then
      put_line("-> in quaddobl_newton_matrix_series.LU_Newton_Steps 4 ...");
    end if;
    LU_Newton_Steps(file,p,jp,degree,maxdeg,nbrit,x,info,vrblvl-1);
    QuadDobl_CSeries_Jaco_Matrices.Clear(jp);
  end LU_Newton_Steps;

-- MANY NEWTON STEPS WITH LU ON COEFFICIENT-PARAMETER HOMOTOPIES :

  procedure LU_Newton_Steps
              ( f : in QuadDobl_CSeries_Poly_SysFun.Eval_Coeff_Poly_Sys;
                c : in QuadDobl_Complex_Series_VecVecs.VecVec;
                ejm : in QuadDobl_CSeries_Jaco_Matrices.Eval_Coeff_Jaco_Mat;
                mlt : in QuadDobl_CSeries_Jaco_Matrices.Mult_Factors;
                degree : in out integer32; maxdeg,nbrit : in integer32;
                x : in out QuadDobl_Complex_Series_Vectors.Vector;
                info : out integer32; vrblvl : in integer32 := 0 ) is
  begin
    if vrblvl > 0 then
      put_line("-> in quaddobl_newton_matrix_series.LU_Newton_Steps 5 ...");
    end if;
    for i in 1..nbrit loop
      LU_Newton_Step(f,c,ejm,mlt,degree,x,info,vrblvl-1);
      exit when (info /= 0); -- stop if Jacobian matrix is singular
      exit when (i = nbrit); -- do not double degree after last step
      Standard_Newton_Matrix_Series.Double_Degree_with_Threshold(degree,maxdeg);
    end loop;
  end LU_Newton_Steps;

  procedure LU_Newton_Steps
              ( file : in file_type;
                f : in QuadDobl_CSeries_Poly_SysFun.Eval_Coeff_Poly_Sys;
                c : in QuadDobl_Complex_Series_VecVecs.VecVec;
                ejm : in QuadDobl_CSeries_Jaco_Matrices.Eval_Coeff_Jaco_Mat;
                mlt : in QuadDobl_CSeries_Jaco_Matrices.Mult_Factors;
                degree : in out integer32; maxdeg,nbrit : in integer32;
                x : in out QuadDobl_Complex_Series_Vectors.Vector;
                info : out integer32; vrblvl : in integer32 := 0 ) is
  begin
    if vrblvl > 0 then
      put_line("-> in quaddobl_newton_matrix_series.LU_Newton_Steps 6 ...");
    end if;
    for i in 1..nbrit loop
      put(file,"LU Newton step "); put(file,i,1); put_line(file," :");
      LU_Newton_Step(file,f,c,ejm,mlt,degree,x,info,vrblvl-1);
      exit when (info /= 0); -- stop if Jacobian matrix is singular
      exit when (i = nbrit); -- do not double degree after last step
      Standard_Newton_Matrix_Series.Double_Degree_with_Threshold(degree,maxdeg);
    end loop;
  end LU_Newton_Steps;

-- MANY NEWTON STEPS WITH LU WITH CONDITION NUMBER ESTIMATE :

  procedure LU_Newton_Steps
              ( p : in QuadDobl_CSeries_Poly_Systems.Poly_Sys;
                jp : in QuadDobl_CSeries_Jaco_Matrices.Jaco_Mat;
                degree : in out integer32; maxdeg,nbrit : in integer32;
                x : in out QuadDobl_Complex_Series_Vectors.Vector;
                rcond : out quad_double; vrblvl : in integer32 := 0 ) is

    one : constant quad_double := create(1.0);

  begin
    if vrblvl > 0 then
      put_line("-> in quaddobl_newton_matrix_series.LU_Newton_Steps 7 ...");
    end if;
    for i in 1..nbrit loop
      LU_Newton_Step(p,jp,degree,x,rcond,vrblvl-1);
      exit when (one + rcond = one); -- stop if Jacobian matrix is singular
      exit when (i = nbrit); -- do not double degree after last step
      Standard_Newton_Matrix_Series.Double_Degree_with_Threshold(degree,maxdeg);
    end loop;
  end LU_Newton_Steps;

  procedure LU_Newton_Steps
              ( p : in QuadDobl_CSeries_Poly_Systems.Poly_Sys;
                degree : in out integer32; maxdeg,nbrit : in integer32;
                x : in out QuadDobl_Complex_Series_Vectors.Vector;
                rcond : out quad_double; vrblvl : in integer32 := 0 ) is

    jp : QuadDobl_CSeries_Jaco_Matrices.Jaco_Mat(p'range,x'range)
       := QuadDobl_CSeries_Jaco_Matrices.Create(p);

  begin
    if vrblvl > 0 then
      put_line("-> in quaddobl_newton_matrix_series.LU_Newton_Steps 8 ...");
    end if;
    LU_Newton_Steps(p,jp,degree,maxdeg,nbrit,x,rcond,vrblvl-1);
    QuadDobl_CSeries_Jaco_Matrices.Clear(jp);
  end LU_Newton_Steps;

  procedure LU_Newton_Steps
              ( file : in file_type;
                p : in QuadDobl_CSeries_Poly_Systems.Poly_Sys;
                jp : in QuadDobl_CSeries_Jaco_Matrices.Jaco_Mat;
                degree : in out integer32; maxdeg,nbrit : in integer32;
                x : in out QuadDobl_Complex_Series_Vectors.Vector;
                rcond : out quad_double; vrblvl : in integer32 := 0 ) is

    one : constant quad_double := create(1.0);

  begin
    if vrblvl > 0 then
      put_line("-> in quaddobl_newton_matrix_series.LU_Newton_Steps 9 ...");
    end if;
    for i in 1..nbrit loop
      put(file,"LU Newton step "); put(file,i,1); put_line(file," :");
      LU_Newton_Step(file,p,jp,degree,x,rcond,vrblvl-1);
      exit when (one + rcond = one); -- stop if Jacobian matrix is singular
      exit when (i = nbrit); -- do not double degree after last step
      Standard_Newton_Matrix_Series.Double_Degree_with_Threshold(degree,maxdeg);
    end loop;
  end LU_Newton_Steps;

  procedure LU_Newton_Steps
              ( file : in file_type;
                p : in QuadDobl_CSeries_Poly_Systems.Poly_Sys;
                degree : in out integer32; maxdeg,nbrit : in integer32;
                x : in out QuadDobl_Complex_Series_Vectors.Vector;
                rcond : out quad_double; vrblvl : in integer32 := 0 ) is

    jp : QuadDobl_CSeries_Jaco_Matrices.Jaco_Mat(p'range,x'range)
       := QuadDobl_CSeries_Jaco_Matrices.Create(p);

  begin
    if vrblvl > 0 then
      put_line("-> in quaddobl_newton_matrix_series.LU_Newton_Steps 10 ...");
    end if;
    LU_Newton_Steps(file,p,jp,degree,maxdeg,nbrit,x,rcond,vrblvl-1);
    QuadDobl_CSeries_Jaco_Matrices.Clear(jp);
  end LU_Newton_Steps;

-- MANY NEWTON STEPS WITH QR :

  procedure QR_Newton_Steps
              ( p : in QuadDobl_CSeries_Poly_Systems.Poly_Sys;
                jp : in QuadDobl_CSeries_Jaco_Matrices.Jaco_Mat;
                degree : in out integer32; maxdeg,nbrit : in integer32;
                x : in out QuadDobl_Complex_Series_Vectors.Vector;
                info : out integer32; vrblvl : in integer32 := 0 ) is
  begin
    if vrblvl > 0 then
      put_line("-> in quaddobl_newton_matrix_series.QR_Newton_Steps 1 ...");
    end if;
    for i in 1..nbrit loop
      QR_Newton_Step(p,jp,degree,x,info,vrblvl-1);
      exit when (info /= 0); -- stop if Jacobian matrix is singular
      exit when (i = nbrit); -- do not double degree after last step
      Standard_Newton_Matrix_Series.Double_Degree_with_Threshold(degree,maxdeg);
    end loop;
  end QR_Newton_Steps;

  procedure QR_Newton_Steps
              ( p : in QuadDobl_CSeries_Poly_Systems.Poly_Sys;
                degree : in out integer32; maxdeg,nbrit : in integer32;
                x : in out QuadDobl_Complex_Series_Vectors.Vector;
                info : out integer32; vrblvl : in integer32 := 0 ) is

    jp : QuadDobl_CSeries_Jaco_Matrices.Jaco_Mat(p'range,x'range)
       := QuadDobl_CSeries_Jaco_Matrices.Create(p);

  begin
    if vrblvl > 0 then
      put_line("-> in quaddobl_newton_matrix_series.QR_Newton_Steps 2 ...");
    end if;
    QR_Newton_Steps(p,jp,degree,maxdeg,nbrit,x,info,vrblvl-1);
    QuadDobl_CSeries_Jaco_Matrices.Clear(jp);
  end QR_Newton_Steps;

  procedure QR_Newton_Steps
              ( file : in file_type;
                p : in QuadDobl_CSeries_Poly_Systems.Poly_Sys;
                jp : in QuadDobl_CSeries_Jaco_Matrices.Jaco_Mat;
                degree : in out integer32; maxdeg,nbrit : in integer32;
                x : in out QuadDobl_Complex_Series_Vectors.Vector;
                info : out integer32; vrblvl : in integer32 := 0 ) is
  begin
    if vrblvl > 0 then
      put_line("-> in quaddobl_newton_matrix_series.QR_Newton_Steps 3 ...");
    end if;
    for i in 1..nbrit loop
      put(file,"QR Newton step "); put(file,i,1); put_line(file," :");
      QR_Newton_Step(file,p,jp,degree,x,info,vrblvl-1);
      exit when (info /= 0); -- stop if Jacobian matrix is singular
      exit when (i = nbrit); -- do not double degree after last step
      Standard_Newton_Matrix_Series.Double_Degree_with_Threshold(degree,maxdeg);
    end loop;
  end QR_Newton_Steps;

  procedure QR_Newton_Steps
              ( file : in file_type;
                p : in QuadDobl_CSeries_Poly_Systems.Poly_Sys;
                degree : in out integer32; maxdeg,nbrit : in integer32;
                x : in out QuadDobl_Complex_Series_Vectors.Vector;
                info : out integer32; vrblvl : in integer32 := 0 ) is

    jp : QuadDobl_CSeries_Jaco_Matrices.Jaco_Mat(p'range,x'range)
       := QuadDobl_CSeries_Jaco_Matrices.Create(p);

  begin
    if vrblvl > 0 then
      put_line("-> in quaddobl_newton_matrix_series.QR_Newton_Steps 4 ...");
    end if;
    QR_Newton_Steps(file,p,jp,degree,maxdeg,nbrit,x,info,vrblvl-1);
    QuadDobl_CSeries_Jaco_Matrices.Clear(jp);
  end QR_Newton_Steps;

-- MANY QR NEWTON STEPS ON COEFFICIENT-PARAMETER HOMOTOPIES :

  procedure QR_Newton_Steps
              ( f : in QuadDobl_CSeries_Poly_SysFun.Eval_Coeff_Poly_Sys;
                c : in QuadDobl_Complex_Series_VecVecs.VecVec;
                ejm : in QuadDobl_CSeries_Jaco_Matrices.Eval_Coeff_Jaco_Mat;
                mlt : in QuadDobl_CSeries_Jaco_Matrices.Mult_Factors;
                degree : in out integer32; maxdeg,nbrit : in integer32;
                x : in out QuadDobl_Complex_Series_Vectors.Vector;
                info : out integer32; vrblvl : in integer32 := 0 ) is
  begin
    if vrblvl > 0 then
      put_line("-> in quaddobl_newton_matrix_series.QR_Newton_Steps 5 ...");
    end if;
    for i in 1..nbrit loop
      QR_Newton_Step(f,c,ejm,mlt,degree,x,info,vrblvl-1);
      exit when (info /= 0); -- stop if Jacobian matrix is singular
      exit when (i = nbrit); -- do not double degree after last step
      Standard_Newton_Matrix_Series.Double_Degree_with_Threshold(degree,maxdeg);
    end loop;
  end QR_Newton_Steps;

  procedure QR_Newton_Steps
              ( file : in file_type;
                f : in QuadDobl_CSeries_Poly_SysFun.Eval_Coeff_Poly_Sys;
                c : in QuadDobl_Complex_Series_VecVecs.VecVec;
                ejm : in QuadDobl_CSeries_Jaco_Matrices.Eval_Coeff_Jaco_Mat;
                mlt : in QuadDobl_CSeries_Jaco_Matrices.Mult_Factors;
                degree : in out integer32; maxdeg,nbrit : in integer32;
                x : in out QuadDobl_Complex_Series_Vectors.Vector;
                info : out integer32; vrblvl : in integer32 := 0 ) is
  begin
    if vrblvl > 0 then
      put_line("-> in quaddobl_newton_matrix_series.QR_Newton_Steps 6 ...");
    end if;
    for i in 1..nbrit loop
      put(file,"QR Newton step "); put(file,i,1); put_line(file," :");
      QR_Newton_Step(file,f,c,ejm,mlt,degree,x,info,vrblvl-1);
      exit when (info /= 0); -- stop if Jacobian matrix is singular
      exit when (i = nbrit); -- do not double degree after last step
      Standard_Newton_Matrix_Series.Double_Degree_with_Threshold(degree,maxdeg);
    end loop;
  end QR_Newton_Steps;

-- MANY NEWTON STEPS WITH SVD :

  procedure SVD_Newton_Steps
              ( p : in QuadDobl_CSeries_Poly_Systems.Poly_Sys;
                jp : in QuadDobl_CSeries_Jaco_Matrices.Jaco_Mat;
                degree : in out integer32; maxdeg,nbrit : in integer32;
                x : in out QuadDobl_Complex_Series_Vectors.Vector;
                info : out integer32; rcond : out quad_double;
                vrblvl : in integer32 := 0 ) is

    one : constant quad_double := create(1.0);

  begin
    if vrblvl > 0 then
      put_line("-> in quaddobl_newton_matrix_series.SVD_Newton_Steps 1 ...");
    end if;
    for i in 1..nbrit loop
      SVD_Newton_Step(p,jp,degree,x,info,rcond,vrblvl-1);
      exit when (one + rcond = one); -- stop if Jacobian matrix is singular
      exit when (i = nbrit); -- do not double degree after last step
      Standard_Newton_Matrix_Series.Double_Degree_with_Threshold(degree,maxdeg);
    end loop;
  end SVD_Newton_Steps;

  procedure SVD_Newton_Steps
              ( p : in QuadDobl_CSeries_Poly_Systems.Poly_Sys;
                degree : in out integer32; maxdeg,nbrit : in integer32;
                x : in out QuadDobl_Complex_Series_Vectors.Vector;
                info : out integer32; rcond : out quad_double;
                vrblvl : in integer32 := 0 ) is

    jp : QuadDobl_CSeries_Jaco_Matrices.Jaco_Mat(p'range,x'range)
       := QuadDobl_CSeries_Jaco_Matrices.Create(p);

  begin
    if vrblvl > 0 then
      put_line("-> in quaddobl_newton_matrix_series.SVD_Newton_Steps 2 ...");
    end if;
    SVD_Newton_Steps(p,jp,degree,maxdeg,nbrit,x,info,rcond,vrblvl-1);
    QuadDobl_CSeries_Jaco_Matrices.Clear(jp);
  end SVD_Newton_Steps;

  procedure SVD_Newton_Steps
              ( file : in file_type;
                p : in QuadDobl_CSeries_Poly_Systems.Poly_Sys;
                jp : in QuadDobl_CSeries_Jaco_Matrices.Jaco_Mat;
                degree : in out integer32; maxdeg,nbrit : in integer32;
                x : in out QuadDobl_Complex_Series_Vectors.Vector;
                info : out integer32; rcond : out quad_double;
                vrblvl : in integer32 := 0 ) is

    one : constant quad_double := create(1.0);

  begin
    if vrblvl > 0 then
      put_line("-> in quaddobl_newton_matrix_series.SVD_Newton_Steps 3 ...");
    end if;
    for i in 1..nbrit loop
      put(file,"SVD Newton step "); put(file,i,1); put_line(file," :");
      SVD_Newton_Step(file,p,jp,degree,x,info,rcond,vrblvl-1);
      exit when (one + rcond = one); -- stop if Jacobian matrix is singular
      exit when (i = nbrit); -- do not double degree after last step
      Standard_Newton_Matrix_Series.Double_Degree_with_Threshold(degree,maxdeg);
    end loop;
  end SVD_Newton_Steps;

  procedure SVD_Newton_Steps
              ( file : in file_type;
                p : in QuadDobl_CSeries_Poly_Systems.Poly_Sys;
                degree : in out integer32; maxdeg,nbrit : in integer32;
                x : in out QuadDobl_Complex_Series_Vectors.Vector;
                info : out integer32; rcond : out quad_double;
                vrblvl : in integer32 := 0 ) is

    jp : QuadDobl_CSeries_Jaco_Matrices.Jaco_Mat(p'range,x'range)
       := QuadDobl_CSeries_Jaco_Matrices.Create(p);

  begin
    if vrblvl > 0 then
      put_line("-> in quaddobl_newton_matrix_series.SVD_Newton_Steps 4 ...");
    end if;
    SVD_Newton_Steps(file,p,jp,degree,maxdeg,nbrit,x,info,rcond,vrblvl-1);
    QuadDobl_CSeries_Jaco_Matrices.Clear(jp);
  end SVD_Newton_Steps;

-- MANY NEWTON STEPS WITH ECHELON FORM :

  procedure Echelon_Newton_Steps
              ( p : in QuadDobl_CSeries_Poly_Systems.Poly_Sys;
                jp : in QuadDobl_CSeries_Jaco_Matrices.Jaco_Mat;
                degree : in out integer32; maxdeg,nbrit : in integer32;
                x : in out QuadDobl_Complex_Series_Vectors.Vector;
                det : out Complex_Number; vrblvl : in integer32 := 0 ) is
  begin
    if vrblvl > 0 then
      put_line
        ("-> in quaddobl_newton_matrix_series.Echelon_Newton_Steps 1 ...");
    end if;
    for i in 1..nbrit loop
      Echelon_Newton_Step(p,jp,degree,x,det,vrblvl-1);
      exit when (i = nbrit); -- do not double degree after last step
      Standard_Newton_Matrix_Series.Double_Degree_with_Threshold(degree,maxdeg);
    end loop;
  end Echelon_Newton_Steps;

  procedure Echelon_Newton_Steps
              ( p : in QuadDobl_CSeries_Poly_Systems.Poly_Sys;
                degree : in out integer32; maxdeg,nbrit : in integer32;
                x : in out QuadDobl_Complex_Series_Vectors.Vector;
                det : out Complex_Number; vrblvl : in integer32 := 0 ) is

    jp : QuadDobl_CSeries_Jaco_Matrices.Jaco_Mat(p'range,x'range)
       := QuadDobl_CSeries_Jaco_Matrices.Create(p);

  begin
    if vrblvl > 0 then
      put_line
        ("-> in quaddobl_newton_matrix_series.Echelon_Newton_Steps 2 ...");
    end if;
    Echelon_Newton_Steps(p,jp,degree,maxdeg,nbrit,x,det,vrblvl-1);
    QuadDobl_CSeries_Jaco_Matrices.Clear(jp);
  end Echelon_Newton_Steps;

  procedure Echelon_Newton_Steps
              ( file : in file_type;
                p : in QuadDobl_CSeries_Poly_Systems.Poly_Sys;
                jp : in QuadDobl_CSeries_Jaco_Matrices.Jaco_Mat;
                degree : in out integer32; maxdeg,nbrit : in integer32;
                x : in out QuadDobl_Complex_Series_Vectors.Vector;
                det : out Complex_Number; vrblvl : in integer32 := 0 ) is
  begin
    if vrblvl > 0 then
      put_line
        ("-> in quaddobl_newton_matrix_series.Echelon_Newton_Steps 3 ...");
    end if;
    for i in 1..nbrit loop
      put(file,"Echelon Newton step "); put(file,i,1); put_line(file," :");
      Echelon_Newton_Step(file,p,jp,degree,x,det,vrblvl-1);
      exit when (i = nbrit); -- do not double degree after last step
      Standard_Newton_Matrix_Series.Double_Degree_with_Threshold(degree,maxdeg);
    end loop;
  end Echelon_Newton_Steps;

  procedure Echelon_Newton_Steps
              ( file : in file_type;
                p : in QuadDobl_CSeries_Poly_Systems.Poly_Sys;
                degree : in out integer32; maxdeg,nbrit : in integer32;
                x : in out QuadDobl_Complex_Series_Vectors.Vector;
                det : out Complex_Number; vrblvl : in integer32 := 0 ) is

    jp : QuadDobl_CSeries_Jaco_Matrices.Jaco_Mat(p'range,x'range)
       := QuadDobl_CSeries_Jaco_Matrices.Create(p);

  begin
    if vrblvl > 0 then
      put_line
        ("-> in quaddobl_newton_matrix_series.Echelon_Newton_Steps 4 ...");
    end if;
    Echelon_Newton_Steps(file,p,jp,degree,maxdeg,nbrit,x,det,vrblvl-1);
    QuadDobl_CSeries_Jaco_Matrices.Clear(jp);
  end Echelon_Newton_Steps;

end QuadDobl_Newton_Matrix_Series;
