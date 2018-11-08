with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers_io;      use Standard_Floating_Numbers_io;
with Standard_Complex_Numbers_io;       use Standard_Complex_Numbers_io;
with Standard_Complex_Series;
with Standard_Complex_Series_io;
with Standard_CSeries_Vector_Norms;
with Standard_Complex_Series_Matrices;
with Standard_Complex_Vector_Series;
with Standard_Complex_Matrix_Series;
with Standard_Series_Matrix_Solvers;    use Standard_Series_Matrix_Solvers;
with Standard_CSeries_Polynomials;
with Complex_Series_and_Polynomials;
with Standard_CSeries_Poly_SysFun;

package body Standard_Newton_Matrix_Series is

-- ONE NEWTON STEP WITH LU WITHOUT CONDITION NUMBER ESTIMATE :

  procedure LU_Newton_Step
              ( p : in Standard_CSeries_Poly_Systems.Poly_Sys;
                jp : in Standard_CSeries_Jaco_Matrices.Jaco_Mat;
                degree : in integer32;
                x : in out Standard_Complex_Series_Vectors.Vector;
                info : out integer32 ) is

    dx : Standard_Complex_Series_Vectors.Vector(x'range);
    px : Standard_Complex_Series_Vectors.Vector(p'range);
    jm : Standard_Complex_Series_Matrices.Matrix(p'range,x'range);
    xp,xd : Standard_Complex_Vector_Series.Vector(degree);
    mj : Standard_Complex_Matrix_Series.Matrix(degree);

  begin
    px := Standard_CSeries_Poly_SysFun.Eval(p,x);
    Standard_Complex_Series_Vectors.Min(px);
    Complex_Series_and_Polynomials.Set_degree(px,degree);
    jm := Standard_CSeries_Jaco_Matrices.Eval(jp,x);
    Complex_Series_and_Polynomials.Set_degree(jm,degree);
    mj := Standard_Complex_Matrix_Series.Create(jm);
    xp := Standard_Complex_Vector_Series.Create(px);
    Solve_by_lufac(mj,xp,info,xd);
    if info = 0 then
      dx := Standard_Complex_Vector_Series.Create(xd);
      Standard_Complex_Series_Vectors.Add(x,dx);
    end if;
    Standard_Complex_Matrix_Series.Clear(mj);
    Standard_Complex_Vector_Series.Clear(xp);
    Standard_Complex_Vector_Series.Clear(xd);
  end LU_Newton_Step;

  procedure LU_Newton_Step
              ( p : in Standard_CSeries_Poly_Systems.Poly_Sys;
                degree : in integer32;
                x : in out Standard_Complex_Series_Vectors.Vector;
                info : out integer32 ) is

    jp : Standard_CSeries_Jaco_Matrices.Jaco_Mat(p'range,x'range)
       := Standard_CSeries_Jaco_Matrices.Create(p);

  begin
    LU_Newton_Step(p,jp,degree,x,info);
    Standard_CSeries_Jaco_Matrices.Clear(jp);
  end LU_Newton_Step;

  procedure LU_Newton_Step
              ( file : in file_type;
                p : in Standard_CSeries_Poly_Systems.Poly_Sys;
                jp : in Standard_CSeries_Jaco_Matrices.Jaco_Mat;
                degree : in integer32;
                x : in out Standard_Complex_Series_Vectors.Vector;
                info : out integer32 ) is

    dx : Standard_Complex_Series_Vectors.Vector(x'range);
    px : Standard_Complex_Series_Vectors.Vector(p'range);
    jm : Standard_Complex_Series_Matrices.Matrix(p'range,x'range);
    xp,xd : Standard_Complex_Vector_Series.Vector(degree);
    mj : Standard_Complex_Matrix_Series.Matrix(degree);

  begin
    px := Standard_CSeries_Poly_SysFun.Eval(p,x);
    put_line(file,"The evaluated series :");
    for i in px'range loop
      Standard_Complex_Series_io.put(file,px(i)); new_line(file);
    end loop;
    Standard_Complex_Series_Vectors.Min(px);
    Complex_Series_and_Polynomials.Set_degree(px,degree);
    jm := Standard_CSeries_Jaco_Matrices.Eval(jp,x);
    Complex_Series_and_Polynomials.Set_degree(jm,degree);
    mj := Standard_Complex_Matrix_Series.Create(jm);
    xp := Standard_Complex_Vector_Series.Create(px);
    Solve_by_lufac(mj,xp,info,xd);
    if info /= 0 then
      put(file,"LUfac info : "); put(file,info,1); new_line(file);
    else
      dx := Standard_Complex_Vector_Series.Create(xd);
      put_line(file,"The update to the series :");
      for i in dx'range loop
        Standard_Complex_Series_io.put(file,dx(i)); new_line(file);
      end loop;
      Standard_Complex_Series_Vectors.Add(x,dx);
    end if;
    Standard_Complex_Matrix_Series.Clear(mj);
    Standard_Complex_Vector_Series.Clear(xp);
    Standard_Complex_Vector_Series.Clear(xd);
  end LU_Newton_Step;

  procedure LU_Newton_Step
              ( file : in file_type;
                p : in Standard_CSeries_Poly_Systems.Poly_Sys;
                degree : in integer32;
                x : in out Standard_Complex_Series_Vectors.Vector;
                info : out integer32 ) is

    jp : Standard_CSeries_Jaco_Matrices.Jaco_Mat(p'range,x'range)
       := Standard_CSeries_Jaco_Matrices.Create(p);

  begin
    LU_Newton_Step(file,p,jp,degree,x,info);
    Standard_CSeries_Jaco_Matrices.Clear(jp);
  end LU_Newton_Step;

-- ONE NEWTON STEP WITH LU WITH CONDITION NUMBER ESTIMATE :

  procedure LU_Newton_Step
              ( p : in Standard_CSeries_Poly_Systems.Poly_Sys;
                jp : in Standard_CSeries_Jaco_Matrices.Jaco_Mat;
                degree : in integer32;
                x : in out Standard_Complex_Series_Vectors.Vector;
                rcond : out double_float ) is

    dx : Standard_Complex_Series_Vectors.Vector(x'range);
    px : Standard_Complex_Series_Vectors.Vector(p'range);
    jm : Standard_Complex_Series_Matrices.Matrix(p'range,x'range);
    xp,xd : Standard_Complex_Vector_Series.Vector(degree);
    mj : Standard_Complex_Matrix_Series.Matrix(degree);

  begin
    px := Standard_CSeries_Poly_SysFun.Eval(p,x);
    Standard_Complex_Series_Vectors.Min(px);
    Complex_Series_and_Polynomials.Set_degree(px,degree);
    jm := Standard_CSeries_Jaco_Matrices.Eval(jp,x);
    Complex_Series_and_Polynomials.Set_degree(jm,degree);
    mj := Standard_Complex_Matrix_Series.Create(jm);
    xp := Standard_Complex_Vector_Series.Create(px);
    Solve_by_lufco(mj,xp,rcond,xd);
    if 1.0 + rcond /= 1.0 then
      dx := Standard_Complex_Vector_Series.Create(xd);
      Standard_Complex_Series_Vectors.Add(x,dx);
    end if;
    Standard_Complex_Matrix_Series.Clear(mj);
    Standard_Complex_Vector_Series.Clear(xp);
    Standard_Complex_Vector_Series.Clear(xd);
  end LU_Newton_Step;

  procedure LU_Newton_Step
              ( p : in Standard_CSeries_Poly_Systems.Poly_Sys;
                degree : in integer32;
                x : in out Standard_Complex_Series_Vectors.Vector;
                rcond : out double_float ) is

    jp : Standard_CSeries_Jaco_Matrices.Jaco_Mat(p'range,x'range)
       := Standard_CSeries_Jaco_Matrices.Create(p);

  begin
    LU_Newton_Step(p,jp,degree,x,rcond);
    Standard_CSeries_Jaco_Matrices.Clear(jp);
  end LU_Newton_Step;

  procedure LU_Newton_Step
              ( file : in file_type;
                p : in Standard_CSeries_Poly_Systems.Poly_Sys;
                jp : in Standard_CSeries_Jaco_Matrices.Jaco_Mat;
                degree : in integer32;
                x : in out Standard_Complex_Series_Vectors.Vector;
                rcond : out double_float ) is

    dx : Standard_Complex_Series_Vectors.Vector(x'range);
    px : Standard_Complex_Series_Vectors.Vector(p'range);
    jm : Standard_Complex_Series_Matrices.Matrix(p'range,x'range);
    xp,xd : Standard_Complex_Vector_Series.Vector(degree);
    mj : Standard_Complex_Matrix_Series.Matrix(degree);

  begin
    px := Standard_CSeries_Poly_SysFun.Eval(p,x);
    put_line(file,"The evaluated series :");
    for i in px'range loop
      Standard_Complex_Series_io.put(file,px(i)); new_line(file);
    end loop;
    Standard_Complex_Series_Vectors.Min(px);
    Complex_Series_and_Polynomials.Set_degree(px,degree);
    jm := Standard_CSeries_Jaco_Matrices.Eval(jp,x);
    Complex_Series_and_Polynomials.Set_degree(jm,degree);
    mj := Standard_Complex_Matrix_Series.Create(jm);
    xp := Standard_Complex_Vector_Series.Create(px);
    Solve_by_lufco(mj,xp,rcond,xd);
    put(file,"LUfco rcond : "); put(file,rcond,3); new_line(file);
    if 1.0 + rcond /= 1.0 then
      dx := Standard_Complex_Vector_Series.Create(xd);
      put_line(file,"The update to the series :");
      for i in dx'range loop
        Standard_Complex_Series_io.put(file,dx(i)); new_line(file);
      end loop;
      Standard_Complex_Series_Vectors.Add(x,dx);
    end if;
    Standard_Complex_Matrix_Series.Clear(mj);
    Standard_Complex_Vector_Series.Clear(xp);
    Standard_Complex_Vector_Series.Clear(xd);
  end LU_Newton_Step;

  procedure LU_Newton_Step
              ( file : in file_type;
                p : in Standard_CSeries_Poly_Systems.Poly_Sys;
                degree : in integer32;
                x : in out Standard_Complex_Series_Vectors.Vector;
                rcond : out double_float ) is

    jp : Standard_CSeries_Jaco_Matrices.Jaco_Mat(p'range,x'range)
       := Standard_CSeries_Jaco_Matrices.Create(p);

  begin
    LU_Newton_Step(file,p,jp,degree,x,rcond);
    Standard_CSeries_Jaco_Matrices.Clear(jp);
  end LU_Newton_Step;

-- ONE NEWTON STEP WITH QR :

  procedure QR_Newton_Step
              ( p : in Standard_CSeries_Poly_Systems.Poly_Sys;
                jp : in Standard_CSeries_Jaco_Matrices.Jaco_Mat;
                degree : in integer32;
                x : in out Standard_Complex_Series_Vectors.Vector;
                info : out integer32 ) is

    wp : Standard_CSeries_Poly_Systems.Poly_Sys(p'range);
    dx : Standard_Complex_Series_Vectors.Vector(x'range);
    px : Standard_Complex_Series_Vectors.Vector(p'range);
    jm : Standard_Complex_Series_Matrices.Matrix(p'range,x'range);
    mj : Standard_Complex_Matrix_Series.Matrix(degree);
    xp,xd : Standard_Complex_Vector_Series.Vector(degree);
    nrm : double_float;
    tol : constant double_float := 1.0E-13;
    wjp : Standard_CSeries_Jaco_Matrices.Jaco_Mat(jp'range(1),jp'range(2));

  begin
    Standard_CSeries_Poly_Systems.Copy(p,wp);
    Complex_Series_and_Polynomials.Set_degree(wp,degree);
    px := Standard_CSeries_Poly_SysFun.Eval(wp,x);
    nrm := Standard_CSeries_Vector_Norms.Max_Norm(px);
    if nrm > tol then
      for i in jp'range(1) loop
        for j in jp'range(2) loop
          Standard_CSeries_Polynomials.Copy(jp(i,j),wjp(i,j));
        end loop;
      end loop;
      Complex_Series_and_Polynomials.Set_degree(wjp,degree);
      Standard_Complex_Series_Vectors.Min(px);
      Complex_Series_and_Polynomials.Set_degree(px,degree);
      jm := Standard_CSeries_Jaco_Matrices.Eval(wjp,x);
      Complex_Series_and_Polynomials.Set_degree(jm,degree);
      mj := Standard_Complex_Matrix_Series.Create(jm);
      xp := Standard_Complex_Vector_Series.Create(px);
      Solve_by_QRLS(mj,xp,info,xd);
      if info = 0 then
        dx := Standard_Complex_Vector_Series.Create(xd);
        Standard_Complex_Series_Vectors.Add(x,dx);
      end if;
      Standard_CSeries_Jaco_Matrices.Clear(wjp);
    end if;
    Standard_CSeries_Poly_Systems.Clear(wp);
    Standard_Complex_Matrix_Series.Clear(mj);
    Standard_Complex_Vector_Series.Clear(xp);
    Standard_Complex_Vector_Series.Clear(xd);
  end QR_Newton_Step;

  procedure QR_Newton_Step
              ( p : in Standard_CSeries_Poly_Systems.Poly_Sys;
                degree : in integer32;
                x : in out Standard_Complex_Series_Vectors.Vector;
                info : out integer32 ) is

    jp : Standard_CSeries_Jaco_Matrices.Jaco_Mat(p'range,x'range)
       := Standard_CSeries_Jaco_Matrices.Create(p);

  begin
    QR_Newton_Step(p,jp,degree,x,info);
    Standard_CSeries_Jaco_Matrices.Clear(jp);
  end QR_Newton_Step;

  procedure QR_Newton_Step
              ( file : in file_type;
                p : in Standard_CSeries_Poly_Systems.Poly_Sys;
                jp : in Standard_CSeries_Jaco_Matrices.Jaco_Mat;
                degree : in integer32;
                x : in out Standard_Complex_Series_Vectors.Vector;
                info : out integer32 ) is

    wp : Standard_CSeries_Poly_Systems.Poly_Sys(p'range);
    dx : Standard_Complex_Series_Vectors.Vector(x'range);
    px : Standard_Complex_Series_Vectors.Vector(p'range);
    jm : Standard_Complex_Series_Matrices.Matrix(p'range,x'range);
    mj : Standard_Complex_Matrix_Series.Matrix(degree);
    xp,xd : Standard_Complex_Vector_Series.Vector(degree);
    nrm : double_float;
    tol : constant double_float := 1.0E-13;
    wjp : Standard_CSeries_Jaco_Matrices.Jaco_Mat(jm'range(1),jm'range(2));

  begin
    Standard_CSeries_Poly_Systems.Copy(p,wp);
    Complex_Series_and_Polynomials.Set_degree(wp,degree);
    px := Standard_CSeries_Poly_SysFun.Eval(wp,x);
    put_line(file,"The evaluated series :");
    for i in px'range loop
      Standard_Complex_Series_io.put(file,px(i)); new_line(file);
    end loop;
    nrm := Standard_CSeries_Vector_Norms.Max_Norm(px);
    put(file,"The max norm of the evaluation : ");
    put(file,nrm,3); new_line(file);
    if nrm > tol then
      for i in jp'range(1) loop
        for j in jp'range(2) loop
          Standard_CSeries_Polynomials.Copy(jp(i,j),wjp(i,j));
        end loop;
      end loop;
      Complex_Series_and_Polynomials.Set_degree(wjp,degree);
      Standard_Complex_Series_Vectors.Min(px);
      Complex_Series_and_Polynomials.Set_degree(px,degree);
      jm := Standard_CSeries_Jaco_Matrices.Eval(jp,x);
      Complex_Series_and_Polynomials.Set_degree(jm,degree);
      mj := Standard_Complex_Matrix_Series.Create(jm);
      xp := Standard_Complex_Vector_Series.Create(px);
      Solve_by_QRLS(mj,xp,info,xd);
      if info /= 0 then
        put(file,"QRLS info : "); put(file,info,1); new_line(file);
      else
        dx := Standard_Complex_Vector_Series.Create(xd);
        put_line(file,"The update to the series :");
        for i in dx'range loop
          Standard_Complex_Series_io.put(file,dx(i)); new_line(file);
        end loop;
        Standard_Complex_Series_Vectors.Add(x,dx);
      end if;
      Standard_CSeries_Jaco_Matrices.Clear(wjp);
    end if;
    Standard_CSeries_Poly_Systems.Clear(wp);
    Standard_Complex_Matrix_Series.Clear(mj);
    Standard_Complex_Vector_Series.Clear(xp);
    Standard_Complex_Vector_Series.Clear(xd);
  end QR_Newton_Step;

  procedure QR_Newton_Step
              ( file : in file_type;
                p : in Standard_CSeries_Poly_Systems.Poly_Sys;
                degree : in integer32;
                x : in out Standard_Complex_Series_Vectors.Vector;
                info : out integer32 ) is

    jp : Standard_CSeries_Jaco_Matrices.Jaco_Mat(p'range,x'range)
       := Standard_CSeries_Jaco_Matrices.Create(p);

  begin
    QR_Newton_Step(file,p,jp,degree,x,info);
    Standard_CSeries_Jaco_Matrices.Clear(jp);
  end QR_Newton_Step;

-- ONE NEWTON STEP WITH SVD :

  procedure SVD_Newton_Step
              ( p : in Standard_CSeries_Poly_Systems.Poly_Sys;
                jp : in Standard_CSeries_Jaco_Matrices.Jaco_Mat;
                degree : in integer32;
                x : in out Standard_Complex_Series_Vectors.Vector;
                info : out integer32; rcond : out double_float ) is

    wp : Standard_CSeries_Poly_Systems.Poly_Sys(p'range);
    dx : Standard_Complex_Series_Vectors.Vector(x'range);
    px : Standard_Complex_Series_Vectors.Vector(p'range);
    jm : Standard_Complex_Series_Matrices.Matrix(p'range,x'range);
    mj : Standard_Complex_Matrix_Series.Matrix(degree);
    xp,xd : Standard_Complex_Vector_Series.Vector(degree);
    nrm : double_float;
    tol : constant double_float := 1.0E-13;
    wjp : Standard_CSeries_Jaco_Matrices.Jaco_Mat(jp'range(1),jp'range(2));

  begin
    Standard_CSeries_Poly_Systems.Copy(p,wp);
    Complex_Series_and_Polynomials.Set_degree(wp,degree);
    px := Standard_CSeries_Poly_SysFun.Eval(wp,x);
    nrm := Standard_CSeries_Vector_Norms.Max_Norm(px);
    if nrm > tol then
      for i in jp'range(1) loop
        for j in jp'range(2) loop
          Standard_CSeries_Polynomials.Copy(jp(i,j),wjp(i,j));
        end loop;
      end loop;
      Complex_Series_and_Polynomials.Set_degree(wjp,degree);
      Standard_Complex_Series_Vectors.Min(px);
      Complex_Series_and_Polynomials.Set_degree(px,degree);
      jm := Standard_CSeries_Jaco_Matrices.Eval(wjp,x);
      Complex_Series_and_Polynomials.Set_degree(jm,degree);
      mj := Standard_Complex_Matrix_Series.Create(jm);
      xp := Standard_Complex_Vector_Series.Create(px);
      Solve_by_SVD(mj,xp,info,rcond,xd);
      if 1.0 + rcond /= 1.0 then
        dx := Standard_Complex_Vector_Series.Create(xd);
        Standard_Complex_Series_Vectors.Add(x,dx);
      end if;
      Standard_CSeries_Jaco_Matrices.Clear(wjp);
    end if;
    Standard_CSeries_Poly_Systems.Clear(wp);
    Standard_Complex_Matrix_Series.Clear(mj);
    Standard_Complex_Vector_Series.Clear(xp);
    Standard_Complex_Vector_Series.Clear(xd);
  end SVD_Newton_Step;

  procedure SVD_Newton_Step
              ( p : in Standard_CSeries_Poly_Systems.Poly_Sys;
                degree : in integer32;
                x : in out Standard_Complex_Series_Vectors.Vector;
                info : out integer32; rcond : out double_float ) is

    jp : Standard_CSeries_Jaco_Matrices.Jaco_Mat(p'range,x'range)
       := Standard_CSeries_Jaco_Matrices.Create(p);

  begin
    SVD_Newton_Step(p,jp,degree,x,info,rcond);
    Standard_CSeries_Jaco_Matrices.Clear(jp);
  end SVD_Newton_Step;

  procedure SVD_Newton_Step
              ( file : in file_type;
                p : in Standard_CSeries_Poly_Systems.Poly_Sys;
                jp : in Standard_CSeries_Jaco_Matrices.Jaco_Mat;
                degree : in integer32;
                x : in out Standard_Complex_Series_Vectors.Vector;
                info : out integer32; rcond : out double_float ) is

    wp : Standard_CSeries_Poly_Systems.Poly_Sys(p'range);
    dx : Standard_Complex_Series_Vectors.Vector(x'range);
    px : Standard_Complex_Series_Vectors.Vector(p'range);
    jm : Standard_Complex_Series_Matrices.Matrix(p'range,x'range);
    mj : Standard_Complex_Matrix_Series.Matrix(degree);
    xp,xd : Standard_Complex_Vector_Series.Vector(degree);
    nrm : double_float;
    tol : constant double_float := 1.0E-13;
    wjp : Standard_CSeries_Jaco_Matrices.Jaco_Mat(jm'range(1),jm'range(2));

  begin
    Standard_CSeries_Poly_Systems.Copy(p,wp);
    Complex_Series_and_Polynomials.Set_degree(wp,degree);
    px := Standard_CSeries_Poly_SysFun.Eval(wp,x);
    put_line(file,"The evaluated series :");
    for i in px'range loop
      Standard_Complex_Series_io.put(file,px(i)); new_line(file);
    end loop;
    nrm := Standard_CSeries_Vector_Norms.Max_Norm(px);
    put(file,"The max norm of the evaluation : ");
    put(file,nrm,3); new_line(file);
    if nrm > tol then
      for i in jp'range(1) loop
        for j in jp'range(2) loop
          Standard_CSeries_Polynomials.Copy(jp(i,j),wjp(i,j));
        end loop;
      end loop;
      Complex_Series_and_Polynomials.Set_degree(wjp,degree);
      Standard_Complex_Series_Vectors.Min(px);
      Complex_Series_and_Polynomials.Set_degree(px,degree);
      jm := Standard_CSeries_Jaco_Matrices.Eval(jp,x);
      Complex_Series_and_Polynomials.Set_degree(jm,degree);
      mj := Standard_Complex_Matrix_Series.Create(jm);
      xp := Standard_Complex_Vector_Series.Create(px);
      Solve_by_SVD(mj,xp,info,rcond,xd);
      put(file,"SVD rcond : "); put(file,rcond,3); new_line(file);
      if 1.0 + rcond /= 1.0 then
        dx := Standard_Complex_Vector_Series.Create(xd);
        put_line(file,"The update to the series :");
        for i in dx'range loop
          Standard_Complex_Series_io.put(file,dx(i)); new_line(file);
        end loop;
        Standard_Complex_Series_Vectors.Add(x,dx);
      end if;
      Standard_CSeries_Jaco_Matrices.Clear(wjp);
    end if;
    Standard_CSeries_Poly_Systems.Clear(wp);
    Standard_Complex_Matrix_Series.Clear(mj);
    Standard_Complex_Vector_Series.Clear(xp);
    Standard_Complex_Vector_Series.Clear(xd);
  end SVD_Newton_Step;

  procedure SVD_Newton_Step
              ( file : in file_type;
                p : in Standard_CSeries_Poly_Systems.Poly_Sys;
                degree : in integer32;
                x : in out Standard_Complex_Series_Vectors.Vector;
                info : out integer32; rcond : out double_float ) is

    jp : Standard_CSeries_Jaco_Matrices.Jaco_Mat(p'range,x'range)
       := Standard_CSeries_Jaco_Matrices.Create(p);

  begin
    SVD_Newton_Step(file,p,jp,degree,x,info,rcond);
    Standard_CSeries_Jaco_Matrices.Clear(jp);
  end SVD_Newton_Step;

-- ONE NEWTON STEP WITH ECHELON FORM :

  procedure Echelon_Newton_Step
              ( p : in Standard_CSeries_Poly_Systems.Poly_Sys;
                jp : in Standard_CSeries_Jaco_Matrices.Jaco_Mat;
                degree : in integer32;
                x : in out Standard_Complex_Series_Vectors.Vector;
                det : out Complex_Number ) is

    dx : Standard_Complex_Series_Vectors.Vector(x'range);
    px : Standard_Complex_Series_Vectors.Vector(p'range);
    jm : Standard_Complex_Series_Matrices.Matrix(p'range,x'range);
    xp,xd,xdn : Standard_Complex_Vector_Series.Vector(degree);
    mj : Standard_Complex_Matrix_Series.Matrix(degree);

  begin
    px := Standard_CSeries_Poly_SysFun.Eval(p,x);
    Standard_Complex_Series_Vectors.Min(px);
    Complex_Series_and_Polynomials.Set_degree(px,degree);
    jm := Standard_CSeries_Jaco_Matrices.Eval(jp,x);
    Complex_Series_and_Polynomials.Set_degree(jm,degree);
    mj := Standard_Complex_Matrix_Series.Create(jm);
    xp := Standard_Complex_Vector_Series.Create(px);
    Echelon_Solve(mj,xp,det,xd,xdn);
    dx := Standard_Complex_Vector_Series.Create(xd);
    Standard_Complex_Series_Vectors.Add(x,dx);
    Standard_Complex_Matrix_Series.Clear(mj);
    Standard_Complex_Vector_Series.Clear(xp);
    Standard_Complex_Vector_Series.Clear(xd);
  end Echelon_Newton_Step;

  procedure Echelon_Newton_Step
              ( p : in Standard_CSeries_Poly_Systems.Poly_Sys;
                degree : in integer32;
                x : in out Standard_Complex_Series_Vectors.Vector;
                det : out Complex_Number ) is

    jp : Standard_CSeries_Jaco_Matrices.Jaco_Mat(p'range,x'range)
       := Standard_CSeries_Jaco_Matrices.Create(p);

  begin
    Echelon_Newton_Step(p,jp,degree,x,det);
    Standard_CSeries_Jaco_Matrices.Clear(jp);
  end Echelon_Newton_Step;

  procedure Echelon_Newton_Step
              ( file : in file_type;
                p : in Standard_CSeries_Poly_Systems.Poly_Sys;
                jp : in Standard_CSeries_Jaco_Matrices.Jaco_Mat;
                degree : in integer32;
                x : in out Standard_Complex_Series_Vectors.Vector;
                det : out Complex_Number ) is

    dx : Standard_Complex_Series_Vectors.Vector(x'range);
    px : Standard_Complex_Series_Vectors.Vector(p'range);
    jm : Standard_Complex_Series_Matrices.Matrix(p'range,x'range);
    xp,xd,xdn : Standard_Complex_Vector_Series.Vector(degree);
    mj : Standard_Complex_Matrix_Series.Matrix(degree);

  begin
    px := Standard_CSeries_Poly_SysFun.Eval(p,x);
    put_line(file,"The evaluated series :");
    for i in px'range loop
      Standard_Complex_Series_io.put(file,px(i)); new_line(file);
    end loop;
    Standard_Complex_Series_Vectors.Min(px);
    Complex_Series_and_Polynomials.Set_degree(px,degree);
    put_line(file,"Evaluating the Jacobian matrix ...");
    jm := Standard_CSeries_Jaco_Matrices.Eval(jp,x);
    put_line(file,"Done evaluating the Jacobian matrix.");
    Complex_Series_and_Polynomials.Set_degree(jm,degree);
    put_line(file,"Creating matrices and vectors ...");
    mj := Standard_Complex_Matrix_Series.Create(jm);
    xp := Standard_Complex_Vector_Series.Create(px);
    put_line(file,"Calling Echelon_Solve ...");
    Echelon_Solve(mj,xp,det,xd,xdn);
    put(file,"n.deg : "); put(file,xdn.deg,1); 
    put(file,"  det : "); put(file,det); new_line(file);
    dx := Standard_Complex_Vector_Series.Create(xd);
    Standard_Complex_Series_Vectors.Add(x,dx);
    put_line(file,"The update to the series :");
    for i in dx'range loop
      Standard_Complex_Series_io.put(file,dx(i)); new_line(file);
    end loop;
    Standard_Complex_Matrix_Series.Clear(mj);
    Standard_Complex_Vector_Series.Clear(xp);
    Standard_Complex_Vector_Series.Clear(xd);
  exception
    when others =>
      put("exception in ");
      put_line("Standard_Newton_Matrix_Series.Echelon_Newton_Step");
      raise;
  end Echelon_Newton_Step;

  procedure Echelon_Newton_Step
              ( file : in file_type;
                p : in Standard_CSeries_Poly_Systems.Poly_Sys;
                degree : in integer32;
                x : in out Standard_Complex_Series_Vectors.Vector;
                det : out Complex_Number ) is

    jp : Standard_CSeries_Jaco_Matrices.Jaco_Mat(p'range,x'range)
       := Standard_CSeries_Jaco_Matrices.Create(p);

  begin
    Echelon_Newton_Step(file,p,jp,degree,x,det);
    Standard_CSeries_Jaco_Matrices.Clear(jp);
  end Echelon_Newton_Step;

-- MANY NEWTON STEPS WITH LU WITHOUT CONDITION NUMBER ESTIMATE :

  procedure LU_Newton_Steps
              ( p : in Standard_CSeries_Poly_Systems.Poly_Sys;
                jp : in Standard_CSeries_Jaco_Matrices.Jaco_Mat;
                degree : in out integer32; nbrit : in integer32;
                x : in out Standard_Complex_Series_Vectors.Vector;
                info : out integer32 ) is

  begin
    for i in 1..nbrit loop
      LU_Newton_Step(p,jp,degree,x,info);
      exit when (info /= 0); -- stop if Jacobian matrix is singular
      exit when (i = nbrit); -- do not double degree after last step
      degree := 2*degree;
    end loop;
  end LU_Newton_Steps;

  procedure LU_Newton_Steps
              ( p : in Standard_CSeries_Poly_Systems.Poly_Sys;
                degree : in out integer32; nbrit : in integer32;
                x : in out Standard_Complex_Series_Vectors.Vector;
                info : out integer32 ) is

    jp : Standard_CSeries_Jaco_Matrices.Jaco_Mat(p'range,x'range)
       := Standard_CSeries_Jaco_Matrices.Create(p);

  begin
    LU_Newton_Steps(p,jp,degree,nbrit,x,info);
    Standard_CSeries_Jaco_Matrices.Clear(jp);
  end LU_Newton_Steps;

  procedure LU_Newton_Steps
              ( file : in file_type;
                p : in Standard_CSeries_Poly_Systems.Poly_Sys;
                jp : in Standard_CSeries_Jaco_Matrices.Jaco_Mat;
                degree : in out integer32; nbrit : in integer32;
                x : in out Standard_Complex_Series_Vectors.Vector;
                info : out integer32 ) is
  begin
    for i in 1..nbrit loop
      put(file,"LU Newton step "); put(file,i,1); put_line(file," :");
      LU_Newton_Step(file,p,jp,degree,x,info);
      exit when (info /= 0); -- stop if Jacobian matrix is singular
      exit when (i = nbrit); -- do not double degree after last step
      degree := 2*degree;
    end loop;
  end LU_Newton_Steps;

  procedure LU_Newton_Steps
              ( file : in file_type;
                p : in Standard_CSeries_Poly_Systems.Poly_Sys;
                degree : in out integer32; nbrit : in integer32;
                x : in out Standard_Complex_Series_Vectors.Vector;
                info : out integer32 ) is

    jp : Standard_CSeries_Jaco_Matrices.Jaco_Mat(p'range,x'range)
       := Standard_CSeries_Jaco_Matrices.Create(p);

  begin
    LU_Newton_Steps(file,p,jp,degree,nbrit,x,info);
    Standard_CSeries_Jaco_Matrices.Clear(jp);
  end LU_Newton_Steps;

-- MANY NEWTON STEPS WITH LU WITHOUT CONDITION NUMBER ESTIMATE :

  procedure LU_Newton_Steps
              ( p : in Standard_CSeries_Poly_Systems.Poly_Sys;
                jp : in Standard_CSeries_Jaco_Matrices.Jaco_Mat;
                degree : in out integer32; nbrit : in integer32;
                x : in out Standard_Complex_Series_Vectors.Vector;
                rcond : out double_float ) is
  begin
    for i in 1..nbrit loop
      LU_Newton_Step(p,jp,degree,x,rcond);
      exit when (1.0 + rcond = 1.0); -- stop if Jacobian matrix is singular
      exit when (i = nbrit); -- do not double degree after last step
      degree := 2*degree;
    end loop;
  end LU_Newton_Steps;

  procedure LU_Newton_Steps
              ( p : in Standard_CSeries_Poly_Systems.Poly_Sys;
                degree : in out integer32; nbrit : in integer32;
                x : in out Standard_Complex_Series_Vectors.Vector;
                rcond : out double_float ) is

    jp : Standard_CSeries_Jaco_Matrices.Jaco_Mat(p'range,x'range)
       := Standard_CSeries_Jaco_Matrices.Create(p);

  begin
    LU_Newton_Steps(p,jp,degree,nbrit,x,rcond);
    Standard_CSeries_Jaco_Matrices.Clear(jp);
  end LU_Newton_Steps;

  procedure LU_Newton_Steps
              ( file : in file_type;
                p : in Standard_CSeries_Poly_Systems.Poly_Sys;
                jp : in Standard_CSeries_Jaco_Matrices.Jaco_Mat;
                degree : in out integer32; nbrit : in integer32;
                x : in out Standard_Complex_Series_Vectors.Vector;
                rcond : out double_float ) is
  begin
    for i in 1..nbrit loop
      put(file,"LU Newton step "); put(file,i,1); put_line(file," :");
      LU_Newton_Step(file,p,jp,degree,x,rcond);
      exit when (1.0 + rcond = 1.0); -- stop if Jacobian matrix is singular
      exit when (i = nbrit); -- do not double degree after last step
      degree := 2*degree;
    end loop;
  end LU_Newton_Steps;

  procedure LU_Newton_Steps
              ( file : in file_type;
                p : in Standard_CSeries_Poly_Systems.Poly_Sys;
                degree : in out integer32; nbrit : in integer32;
                x : in out Standard_Complex_Series_Vectors.Vector;
                rcond : out double_float ) is

    jp : Standard_CSeries_Jaco_Matrices.Jaco_Mat(p'range,x'range)
       := Standard_CSeries_Jaco_Matrices.Create(p);

  begin
    LU_Newton_Steps(file,p,jp,degree,nbrit,x,rcond);
    Standard_CSeries_Jaco_Matrices.Clear(jp);
  end LU_Newton_Steps;

-- MANY NEWTON STEPS WITH QR :

  procedure QR_Newton_Steps
              ( p : in Standard_CSeries_Poly_Systems.Poly_Sys;
                jp : in Standard_CSeries_Jaco_Matrices.Jaco_Mat;
                degree : in out integer32; nbrit : in integer32;
                x : in out Standard_Complex_Series_Vectors.Vector;
                info : out integer32 ) is
  begin
    for i in 1..nbrit loop
      QR_Newton_Step(p,jp,degree,x,info);
      exit when (info /= 0); -- stop if Jacobian matrix is singular
      exit when (i = nbrit); -- do not double degree after last step
      degree := 2*degree;
    end loop;
  end QR_Newton_Steps;

  procedure QR_Newton_Steps
              ( p : in Standard_CSeries_Poly_Systems.Poly_Sys;
                degree : in out integer32; nbrit : in integer32;
                x : in out Standard_Complex_Series_Vectors.Vector;
                info : out integer32 ) is

    jp : Standard_CSeries_Jaco_Matrices.Jaco_Mat(p'range,x'range)
       := Standard_CSeries_Jaco_Matrices.Create(p);

  begin
    QR_Newton_Steps(p,jp,degree,nbrit,x,info);
    Standard_CSeries_Jaco_Matrices.Clear(jp);
  end QR_Newton_Steps;

  procedure QR_Newton_Steps
              ( file : in file_type;
                p : in Standard_CSeries_Poly_Systems.Poly_Sys;
                jp : in Standard_CSeries_Jaco_Matrices.Jaco_Mat;
                degree : in out integer32; nbrit : in integer32;
                x : in out Standard_Complex_Series_Vectors.Vector;
                info : out integer32 ) is
  begin
    for i in 1..nbrit loop
      put(file,"QR Newton step "); put(file,i,1); put_line(file," :");
      QR_Newton_Step(file,p,jp,degree,x,info);
      exit when (info /= 0); -- stop if Jacobian matrix is singular
      exit when (i = nbrit); -- do not double degree after last step
      degree := 2*degree;
    end loop;
  end QR_Newton_Steps;

  procedure QR_Newton_Steps
              ( file : in file_type;
                p : in Standard_CSeries_Poly_Systems.Poly_Sys;
                degree : in out integer32; nbrit : in integer32;
                x : in out Standard_Complex_Series_Vectors.Vector;
                info : out integer32 ) is

    jp : Standard_CSeries_Jaco_Matrices.Jaco_Mat(p'range,x'range)
       := Standard_CSeries_Jaco_Matrices.Create(p);

  begin
    QR_Newton_Steps(file,p,jp,degree,nbrit,x,info);
    Standard_CSeries_Jaco_Matrices.Clear(jp);
  end QR_Newton_Steps;

-- MANY NEWTON STEPS WITH SVD :

  procedure SVD_Newton_Steps
              ( p : in Standard_CSeries_Poly_Systems.Poly_Sys;
                jp : in Standard_CSeries_Jaco_Matrices.Jaco_Mat;
                degree : in out integer32; nbrit : in integer32;
                x : in out Standard_Complex_Series_Vectors.Vector;
                info : out integer32; rcond : out double_float ) is
  begin
    for i in 1..nbrit loop
      SVD_Newton_Step(p,jp,degree,x,info,rcond);
      exit when (1.0 + rcond = 1.0); -- stop if Jacobian matrix is singular
      exit when (i = nbrit); -- do not double degree after last step
      degree := 2*degree;
    end loop;
  end SVD_Newton_Steps;

  procedure SVD_Newton_Steps
              ( p : in Standard_CSeries_Poly_Systems.Poly_Sys;
                degree : in out integer32; nbrit : in integer32;
                x : in out Standard_Complex_Series_Vectors.Vector;
                info : out integer32; rcond : out double_float ) is

    jp : Standard_CSeries_Jaco_Matrices.Jaco_Mat(p'range,x'range)
       := Standard_CSeries_Jaco_Matrices.Create(p);

  begin
    SVD_Newton_Steps(p,jp,degree,nbrit,x,info,rcond);
    Standard_CSeries_Jaco_Matrices.Clear(jp);
  end SVD_Newton_Steps;

  procedure SVD_Newton_Steps
              ( file : in file_type;
                p : in Standard_CSeries_Poly_Systems.Poly_Sys;
                jp : in Standard_CSeries_Jaco_Matrices.Jaco_Mat;
                degree : in out integer32; nbrit : in integer32;
                x : in out Standard_Complex_Series_Vectors.Vector;
                info : out integer32; rcond : out double_float ) is
  begin
    for i in 1..nbrit loop
      put(file,"SVD Newton step "); put(file,i,1); put_line(file," :");
      SVD_Newton_Step(file,p,jp,degree,x,info,rcond);
      exit when (1.0 + rcond = 1.0); -- stop if Jacobian matrix is singular
      exit when (i = nbrit); -- do not double degree after last step
      degree := 2*degree;
    end loop;
  end SVD_Newton_Steps;

  procedure SVD_Newton_Steps
              ( file : in file_type;
                p : in Standard_CSeries_Poly_Systems.Poly_Sys;
                degree : in out integer32; nbrit : in integer32;
                x : in out Standard_Complex_Series_Vectors.Vector;
                info : out integer32; rcond : out double_float ) is

    jp : Standard_CSeries_Jaco_Matrices.Jaco_Mat(p'range,x'range)
       := Standard_CSeries_Jaco_Matrices.Create(p);

  begin
    SVD_Newton_Steps(file,p,jp,degree,nbrit,x,info,rcond);
    Standard_CSeries_Jaco_Matrices.Clear(jp);
  end SVD_Newton_Steps;

-- MANY NEWTON STEPS WITH ECHELON FORM :

  procedure Echelon_Newton_Steps
              ( p : in Standard_CSeries_Poly_Systems.Poly_Sys;
                jp : in Standard_CSeries_Jaco_Matrices.Jaco_Mat;
                degree : in out integer32; nbrit : in integer32;
                x : in out Standard_Complex_Series_Vectors.Vector;
                det : out Complex_Number ) is
  begin
    for i in 1..nbrit loop
      Echelon_Newton_Step(p,jp,degree,x,det);
      exit when (i = nbrit); -- do not double degree after last step
      degree := 2*degree;
    end loop;
  end Echelon_Newton_Steps;

  procedure Echelon_Newton_Steps
              ( p : in Standard_CSeries_Poly_Systems.Poly_Sys;
                degree : in out integer32; nbrit : in integer32;
                x : in out Standard_Complex_Series_Vectors.Vector;
                det : out Complex_Number ) is

    jp : Standard_CSeries_Jaco_Matrices.Jaco_Mat(p'range,x'range)
       := Standard_CSeries_Jaco_Matrices.Create(p);

  begin
    Echelon_Newton_Steps(p,jp,degree,nbrit,x,det);
    Standard_CSeries_Jaco_Matrices.Clear(jp);
  end Echelon_Newton_Steps;

  procedure Echelon_Newton_Steps
              ( file : in file_type;
                p : in Standard_CSeries_Poly_Systems.Poly_Sys;
                jp : in Standard_CSeries_Jaco_Matrices.Jaco_Mat;
                degree : in out integer32; nbrit : in integer32;
                x : in out Standard_Complex_Series_Vectors.Vector;
                det : out Complex_Number ) is
  begin
    for i in 1..nbrit loop
      put(file,"Echelon Newton step "); put(file,i,1); put_line(file," :");
      Echelon_Newton_Step(file,p,jp,degree,x,det);
      exit when (i = nbrit); -- do not double degree after last step
      degree := 2*degree;
    end loop;
  exception
    when others =>
      put("exception ");
      put_line("in Standard_Newton_Matrix_Series.Echelon_Newton_Steps");
      raise;
  end Echelon_Newton_Steps;

  procedure Echelon_Newton_Steps
              ( file : in file_type;
                p : in Standard_CSeries_Poly_Systems.Poly_Sys;
                degree : in out integer32; nbrit : in integer32;
                x : in out Standard_Complex_Series_Vectors.Vector;
                det : out Complex_Number ) is

    jp : Standard_CSeries_Jaco_Matrices.Jaco_Mat(p'range,x'range)
       := Standard_CSeries_Jaco_Matrices.Create(p);

  begin
    Echelon_Newton_Steps(file,p,jp,degree,nbrit,x,det);
    Standard_CSeries_Jaco_Matrices.Clear(jp);
  end Echelon_Newton_Steps;

end Standard_Newton_Matrix_Series;
