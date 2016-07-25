with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;         use Standard_Floating_Numbers;
with Quad_Double_Numbers;               use Quad_Double_Numbers;
with Quad_Double_Numbers_io;            use Quad_Double_Numbers_io;
with QuadDobl_Dense_Series;
with QuadDobl_Dense_Series_io;
with QuadDobl_Series_Vector_Norms;
with QuadDobl_Dense_Series_Matrices;
with QuadDobl_Dense_Vector_Series;
with QuadDobl_Dense_Matrix_Series;
with QuadDobl_Matrix_Series_Solvers;    use QuadDobl_Matrix_Series_Solvers;
with QuadDobl_Series_Polynomials;
with Series_and_Polynomials;
with QuadDobl_Series_Poly_SysFun;

package body QuadDobl_Newton_Matrix_Series is

-- ONE NEWTON STEP :

  procedure LU_Newton_Step
              ( p : in QuadDobl_Series_Poly_Systems.Poly_Sys;
                jp : in QuadDobl_Series_Jaco_Matrices.Jaco_Mat;
                degree : in integer32;
                x : in out QuadDobl_Dense_Series_Vectors.Vector;
                info : out integer32 ) is

    dx : QuadDobl_Dense_Series_Vectors.Vector(x'range);
    px : QuadDobl_Dense_Series_Vectors.Vector(p'range);
    jm : QuadDobl_Dense_Series_Matrices.Matrix(p'range,x'range);
    xp,xd : QuadDobl_Dense_Vector_Series.Vector;
    mj : QuadDobl_Dense_Matrix_Series.Matrix;

  begin
    px := QuadDobl_Series_Poly_SysFun.Eval(p,x);
    QuadDobl_Dense_Series_Vectors.Min(px);
    Series_and_Polynomials.Set_degree(px,degree);
    jm := QuadDobl_Series_Jaco_Matrices.Eval(jp,x);
    Series_and_Polynomials.Set_degree(jm,degree);
    mj := QuadDobl_Dense_Matrix_Series.Create(jm);
    xp := QuadDobl_Dense_Vector_Series.Create(px);
    Solve_by_lufac(mj,xp,info,xd);
    if info = 0 then
      dx := QuadDobl_Dense_Vector_Series.Create(xd);
      QuadDobl_Dense_Series_Vectors.Add(x,dx);
    end if;
    QuadDobl_Dense_Matrix_Series.Clear(mj);
    QuadDobl_Dense_Vector_Series.Clear(xp);
    QuadDobl_Dense_Vector_Series.Clear(xd);
  end LU_Newton_Step;

  procedure LU_Newton_Step
              ( p : in QuadDobl_Series_Poly_Systems.Poly_Sys;
                degree : in integer32;
                x : in out QuadDobl_Dense_Series_Vectors.Vector;
                info : out integer32 ) is

    jp : QuadDobl_Series_Jaco_Matrices.Jaco_Mat(p'range,x'range)
       := QuadDobl_Series_Jaco_Matrices.Create(p);

  begin
    LU_Newton_Step(p,jp,degree,x,info);
    QuadDobl_Series_Jaco_Matrices.Clear(jp);
  end LU_Newton_Step;

  procedure LU_Newton_Step
              ( file : in file_type;
                p : in QuadDobl_Series_Poly_Systems.Poly_Sys;
                jp : in QuadDobl_Series_Jaco_Matrices.Jaco_Mat;
                degree : in integer32;
                x : in out QuadDobl_Dense_Series_Vectors.Vector;
                info : out integer32 ) is

    dx : QuadDobl_Dense_Series_Vectors.Vector(x'range);
    px : QuadDobl_Dense_Series_Vectors.Vector(p'range);
    jm : QuadDobl_Dense_Series_Matrices.Matrix(p'range,x'range);
    xp,xd : QuadDobl_Dense_Vector_Series.Vector;
    mj : QuadDobl_Dense_Matrix_Series.Matrix;

  begin
    px := QuadDobl_Series_Poly_SysFun.Eval(p,x);
    put_line(file,"The evaluated series :");
    for i in px'range loop
      QuadDobl_Dense_Series_io.put(file,px(i)); new_line(file);
    end loop;
    QuadDobl_Dense_Series_Vectors.Min(px);
    Series_and_Polynomials.Set_degree(px,degree);
    jm := QuadDobl_Series_Jaco_Matrices.Eval(jp,x);
    Series_and_Polynomials.Set_degree(jm,degree);
    mj := QuadDobl_Dense_Matrix_Series.Create(jm);
    xp := QuadDobl_Dense_Vector_Series.Create(px);
    Solve_by_lufac(mj,xp,info,xd);
    if info /= 0 then
      put(file,"LUfac info : "); put(file,info,1); new_line(file);
    else
      dx := QuadDobl_Dense_Vector_Series.Create(xd);
      put_line(file,"The update to the series :");
      for i in dx'range loop
        QuadDobl_Dense_Series_io.put(file,dx(i)); new_line(file);
      end loop;
      QuadDobl_Dense_Series_Vectors.Add(x,dx);
    end if;
    QuadDobl_Dense_Matrix_Series.Clear(mj);
    QuadDobl_Dense_Vector_Series.Clear(xp);
    QuadDobl_Dense_Vector_Series.Clear(xd);
  end LU_Newton_Step;

  procedure LU_Newton_Step
              ( file : in file_type;
                p : in QuadDobl_Series_Poly_Systems.Poly_Sys;
                degree : in integer32;
                x : in out QuadDobl_Dense_Series_Vectors.Vector;
                info : out integer32 ) is

    jp : QuadDobl_Series_Jaco_Matrices.Jaco_Mat(p'range,x'range)
       := QuadDobl_Series_Jaco_Matrices.Create(p);

  begin
    LU_Newton_Step(file,p,jp,degree,x,info);
    QuadDobl_Series_Jaco_Matrices.Clear(jp);
  end LU_Newton_Step;

  procedure QR_Newton_Step
              ( p : in QuadDobl_Series_Poly_Systems.Poly_Sys;
                jp : in QuadDobl_Series_Jaco_Matrices.Jaco_Mat;
                degree : in integer32;
                x : in out QuadDobl_Dense_Series_Vectors.Vector;
                info : out integer32 ) is

    wp : QuadDobl_Series_Poly_Systems.Poly_Sys(p'range);
    dx : QuadDobl_Dense_Series_Vectors.Vector(x'range);
    px : QuadDobl_Dense_Series_Vectors.Vector(p'range);
    jm : QuadDobl_Dense_Series_Matrices.Matrix(p'range,x'range);
    mj : QuadDobl_Dense_Matrix_Series.Matrix;
    xp,xd : QuadDobl_Dense_Vector_Series.Vector;
    nrm : quad_double;
    tol : constant double_float := 1.0E-48;
    wjp : QuadDobl_Series_Jaco_Matrices.Jaco_Mat(jp'range(1),jp'range(2));

  begin
    QuadDobl_Series_Poly_Systems.Copy(p,wp);
    Series_and_Polynomials.Set_degree(wp,degree);
    px := QuadDobl_Series_Poly_SysFun.Eval(wp,x);
    nrm := QuadDobl_Series_Vector_Norms.Max_Norm(px);
    if nrm > tol then
      for i in jp'range(1) loop
        for j in jp'range(2) loop
          QuadDobl_Series_Polynomials.Copy(jp(i,j),wjp(i,j));
        end loop;
      end loop;
      Series_and_Polynomials.Set_degree(wjp,degree);
      QuadDobl_Dense_Series_Vectors.Min(px);
      Series_and_Polynomials.Set_degree(px,degree);
      jm := QuadDobl_Series_Jaco_Matrices.Eval(wjp,x);
      Series_and_Polynomials.Set_degree(jm,degree);
      mj := QuadDobl_Dense_Matrix_Series.Create(jm);
      xp := QuadDobl_Dense_Vector_Series.Create(px);
      Solve_by_QRLS(mj,xp,info,xd);
      if info = 0 then
        dx := QuadDobl_Dense_Vector_Series.Create(xd);
        QuadDobl_Dense_Series_Vectors.Add(x,dx);
      end if;
      QuadDobl_Series_Jaco_Matrices.Clear(wjp);
    end if;
    QuadDobl_Series_Poly_Systems.Clear(wp);
    QuadDobl_Dense_Matrix_Series.Clear(mj);
    QuadDobl_Dense_Vector_Series.Clear(xp);
    QuadDobl_Dense_Vector_Series.Clear(xd);
  end QR_Newton_Step;

  procedure QR_Newton_Step
              ( p : in QuadDobl_Series_Poly_Systems.Poly_Sys;
                degree : in integer32;
                x : in out QuadDobl_Dense_Series_Vectors.Vector;
                info : out integer32 ) is

    jp : QuadDobl_Series_Jaco_Matrices.Jaco_Mat(p'range,x'range)
       := QuadDobl_Series_Jaco_Matrices.Create(p);

  begin
    QR_Newton_Step(p,jp,degree,x,info);
    QuadDobl_Series_Jaco_Matrices.Clear(jp);
  end QR_Newton_Step;

  procedure QR_Newton_Step
              ( file : in file_type;
                p : in QuadDobl_Series_Poly_Systems.Poly_Sys;
                jp : in QuadDobl_Series_Jaco_Matrices.Jaco_Mat;
                degree : in integer32;
                x : in out QuadDobl_Dense_Series_Vectors.Vector;
                info : out integer32 ) is

    wp : QuadDobl_Series_Poly_Systems.Poly_Sys(p'range);
    dx : QuadDobl_Dense_Series_Vectors.Vector(x'range);
    px : QuadDobl_Dense_Series_Vectors.Vector(p'range);
    jm : QuadDobl_Dense_Series_Matrices.Matrix(p'range,x'range);
    mj : QuadDobl_Dense_Matrix_Series.Matrix;
    xp,xd : QuadDobl_Dense_Vector_Series.Vector;
    nrm : quad_double;
    tol : constant double_float := 1.0E-48;
    wjp : QuadDobl_Series_Jaco_Matrices.Jaco_Mat(jm'range(1),jm'range(2));

  begin
    QuadDobl_Series_Poly_Systems.Copy(p,wp);
    Series_and_Polynomials.Set_degree(wp,degree);
    px := QuadDobl_Series_Poly_SysFun.Eval(wp,x);
    put_line(file,"The evaluated series :");
    for i in px'range loop
      QuadDobl_Dense_Series_io.put(file,px(i)); new_line(file);
    end loop;
    nrm := QuadDobl_Series_Vector_Norms.Max_Norm(px);
    put(file,"The max norm of the evaluation : ");
    put(file,nrm,3); new_line(file);
    if nrm > tol then
      for i in jp'range(1) loop
        for j in jp'range(2) loop
          QuadDobl_Series_Polynomials.Copy(jp(i,j),wjp(i,j));
        end loop;
      end loop;
      Series_and_Polynomials.Set_degree(wjp,degree);
      QuadDobl_Dense_Series_Vectors.Min(px);
      Series_and_Polynomials.Set_degree(px,degree);
      jm := QuadDobl_Series_Jaco_Matrices.Eval(jp,x);
      Series_and_Polynomials.Set_degree(jm,degree);
      mj := QuadDobl_Dense_Matrix_Series.Create(jm);
      xp := QuadDobl_Dense_Vector_Series.Create(px);
      Solve_by_QRLS(mj,xp,info,xd);
      if info /= 0 then
        put(file,"QRLS info : "); put(file,info,1); new_line(file);
      else
        dx := QuadDobl_Dense_Vector_Series.Create(xd);
        put_line(file,"The update to the series :");
        for i in dx'range loop
          QuadDobl_Dense_Series_io.put(file,dx(i)); new_line(file);
        end loop;
        QuadDobl_Dense_Series_Vectors.Add(x,dx);
      end if;
      QuadDobl_Series_Jaco_Matrices.Clear(wjp);
    end if;
    QuadDobl_Series_Poly_Systems.Clear(wp);
    QuadDobl_Dense_Matrix_Series.Clear(mj);
    QuadDobl_Dense_Vector_Series.Clear(xp);
    QuadDobl_Dense_Vector_Series.Clear(xd);
  end QR_Newton_Step;

  procedure QR_Newton_Step
              ( file : in file_type;
                p : in QuadDobl_Series_Poly_Systems.Poly_Sys;
                degree : in integer32;
                x : in out QuadDobl_Dense_Series_Vectors.Vector;
                info : out integer32 ) is

    jp : QuadDobl_Series_Jaco_Matrices.Jaco_Mat(p'range,x'range)
       := QuadDobl_Series_Jaco_Matrices.Create(p);

  begin
    QR_Newton_Step(file,p,jp,degree,x,info);
    QuadDobl_Series_Jaco_Matrices.Clear(jp);
  end QR_Newton_Step;

-- MANY NEWTON STEPS :

  procedure LU_Newton_Steps
              ( p : in QuadDobl_Series_Poly_Systems.Poly_Sys;
                jp : in QuadDobl_Series_Jaco_Matrices.Jaco_Mat;
                degree : in out integer32; nbrit : in integer32;
                x : in out QuadDobl_Dense_Series_Vectors.Vector;
                info : out integer32 ) is

  begin
    for i in 1..nbrit loop
      LU_Newton_Step(p,jp,degree,x,info);
      exit when (info /= 0); -- stop if Jacobian matrix is singular
      exit when (i = nbrit); -- do not double degree after last step
      degree := 2*degree;
      if degree > QuadDobl_Dense_Series.max_deg
       then degree := QuadDobl_Dense_Series.max_deg;
      end if;
    end loop;
  end LU_Newton_Steps;

  procedure LU_Newton_Steps
              ( p : in QuadDobl_Series_Poly_Systems.Poly_Sys;
                degree : in out integer32; nbrit : in integer32;
                x : in out QuadDobl_Dense_Series_Vectors.Vector;
                info : out integer32 ) is

    jp : QuadDobl_Series_Jaco_Matrices.Jaco_Mat(p'range,x'range)
       := QuadDobl_Series_Jaco_Matrices.Create(p);

  begin
    LU_Newton_Steps(p,jp,degree,nbrit,x,info);
    QuadDobl_Series_Jaco_Matrices.Clear(jp);
  end LU_Newton_Steps;

  procedure LU_Newton_Steps
              ( file : in file_type;
                p : in QuadDobl_Series_Poly_Systems.Poly_Sys;
                jp : in QuadDobl_Series_Jaco_Matrices.Jaco_Mat;
                degree : in out integer32; nbrit : in integer32;
                x : in out QuadDobl_Dense_Series_Vectors.Vector;
                info : out integer32 ) is

  begin
    for i in 1..nbrit loop
      put(file,"LU Newton step "); put(file,i,1); put_line(file," :");
      LU_Newton_Step(file,p,jp,degree,x,info);
      exit when (info /= 0); -- stop if Jacobian matrix is singular
      exit when (i = nbrit); -- do not double degree after last step
      degree := 2*degree;
      if degree > QuadDobl_Dense_Series.max_deg
       then degree := QuadDobl_Dense_Series.max_deg;
      end if;
    end loop;
  end LU_Newton_Steps;

  procedure LU_Newton_Steps
              ( file : in file_type;
                p : in QuadDobl_Series_Poly_Systems.Poly_Sys;
                degree : in out integer32; nbrit : in integer32;
                x : in out QuadDobl_Dense_Series_Vectors.Vector;
                info : out integer32 ) is

    jp : QuadDobl_Series_Jaco_Matrices.Jaco_Mat(p'range,x'range)
       := QuadDobl_Series_Jaco_Matrices.Create(p);

  begin
    LU_Newton_Steps(file,p,jp,degree,nbrit,x,info);
    QuadDobl_Series_Jaco_Matrices.Clear(jp);
  end LU_Newton_Steps;

  procedure QR_Newton_Steps
              ( p : in QuadDobl_Series_Poly_Systems.Poly_Sys;
                jp : in QuadDobl_Series_Jaco_Matrices.Jaco_Mat;
                degree : in out integer32; nbrit : in integer32;
                x : in out QuadDobl_Dense_Series_Vectors.Vector;
                info : out integer32 ) is

  begin
    for i in 1..nbrit loop
      QR_Newton_Step(p,jp,degree,x,info);
      exit when (info /= 0); -- stop if Jacobian matrix is singular
      exit when (i = nbrit); -- do not double degree after last step
      degree := 2*degree;
      if degree > QuadDobl_Dense_Series.max_deg
       then degree := QuadDobl_Dense_Series.max_deg;
      end if;
    end loop;
  end QR_Newton_Steps;

  procedure QR_Newton_Steps
              ( p : in QuadDobl_Series_Poly_Systems.Poly_Sys;
                degree : in out integer32; nbrit : in integer32;
                x : in out QuadDobl_Dense_Series_Vectors.Vector;
                info : out integer32 ) is

    jp : QuadDobl_Series_Jaco_Matrices.Jaco_Mat(p'range,x'range)
       := QuadDobl_Series_Jaco_Matrices.Create(p);

  begin
    QR_Newton_Steps(p,jp,degree,nbrit,x,info);
    QuadDobl_Series_Jaco_Matrices.Clear(jp);
  end QR_Newton_Steps;

  procedure QR_Newton_Steps
              ( file : in file_type;
                p : in QuadDobl_Series_Poly_Systems.Poly_Sys;
                jp : in QuadDobl_Series_Jaco_Matrices.Jaco_Mat;
                degree : in out integer32; nbrit : in integer32;
                x : in out QuadDobl_Dense_Series_Vectors.Vector;
                info : out integer32 ) is

  begin
    for i in 1..nbrit loop
      put(file,"QR Newton step "); put(file,i,1); put_line(file," :");
      QR_Newton_Step(file,p,jp,degree,x,info);
      exit when (info /= 0); -- stop if Jacobian matrix is singular
      exit when (i = nbrit); -- do not double degree after last step
      degree := 2*degree;
      if degree > QuadDobl_Dense_Series.max_deg
       then degree := QuadDobl_Dense_Series.max_deg;
      end if;
    end loop;
  end QR_Newton_Steps;

  procedure QR_Newton_Steps
              ( file : in file_type;
                p : in QuadDobl_Series_Poly_Systems.Poly_Sys;
                degree : in out integer32; nbrit : in integer32;
                x : in out QuadDobl_Dense_Series_Vectors.Vector;
                info : out integer32 ) is

    jp : QuadDobl_Series_Jaco_Matrices.Jaco_Mat(p'range,x'range)
       := QuadDobl_Series_Jaco_Matrices.Create(p);

  begin
    QR_Newton_Steps(file,p,jp,degree,nbrit,x,info);
    QuadDobl_Series_Jaco_Matrices.Clear(jp);
  end QR_Newton_Steps;

end QuadDobl_Newton_Matrix_Series;
