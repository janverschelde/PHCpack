with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with Standard_Integer_Vectors;
with QuadDobl_Dense_Series;
with QuadDobl_Dense_Series_io;
with QuadDobl_Dense_Series_Matrices;
with QuadDobl_Linear_Series_Solvers;    use QuadDobl_Linear_Series_Solvers;
with QuadDobl_Least_Squares_Series;     use QuadDobl_Least_Squares_Series;
with Series_and_Polynomials;
with QuadDobl_Series_Poly_SysFun;

package body QuadDobl_Newton_Series is

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
    n : constant integer32 := p'last;
    ipvt : Standard_Integer_Vectors.Vector(1..n);

  begin
    px := QuadDobl_Series_Poly_SysFun.Eval(p,x);
    QuadDobl_Dense_Series_Vectors.Min(px);
    Series_and_Polynomials.Set_degree(px,degree);
    jm := QuadDobl_Series_Jaco_Matrices.Eval(jp,x);
    Series_and_Polynomials.Set_degree(jm,degree);
    LUfac(jm,n,ipvt,info);
    if info = 0 then
      dx := px;
      LUsolve(jm,n,ipvt,dx);
      QuadDobl_Dense_Series_Vectors.Add(x,dx);
    end if;
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
    n : constant integer32 := p'last;
    ipvt : Standard_Integer_Vectors.Vector(1..n);

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
    LUfac(jm,n,ipvt,info);
    if info /= 0 then
      put(file,"LUfac info : "); put(file,info,1); new_line(file);
    else
      dx := px;
      LUsolve(jm,n,ipvt,dx);
      put_line(file,"The update to the series :");
      for i in dx'range loop
        QuadDobl_Dense_Series_io.put(file,dx(i)); new_line(file);
      end loop;
      QuadDobl_Dense_Series_Vectors.Add(x,dx);
    end if;
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

    dx : QuadDobl_Dense_Series_Vectors.Vector(x'range);
    px : QuadDobl_Dense_Series_Vectors.Vector(p'range);
    jm : QuadDobl_Dense_Series_Matrices.Matrix(p'range,x'range);
    ipvt : Standard_Integer_Vectors.Vector(x'range);
    qraux : QuadDobl_Dense_Series_Vectors.Vector(x'range);
    n : constant integer32 := jm'last(1);
    m : constant integer32 := jm'last(2);
    rsd,dum,dum2,dum3 : QuadDobl_Dense_Series_Vectors.Vector(1..n);

  begin
    px := QuadDobl_Series_Poly_SysFun.Eval(p,x);
    QuadDobl_Dense_Series_Vectors.Min(px);
    Series_and_Polynomials.Set_degree(px,degree);
    jm := QuadDobl_Series_Jaco_Matrices.Eval(jp,x);
    Series_and_Polynomials.Set_degree(jm,degree);
    QRD(jm,qraux,ipvt,false);
    QRLS(jm,n,m,qraux,px,dum2,dum3,dx,rsd,dum,110,info);
    if info = 0
     then QuadDobl_Dense_Series_Vectors.Add(x,dx);
    end if;
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

    dx : QuadDobl_Dense_Series_Vectors.Vector(x'range);
    px : QuadDobl_Dense_Series_Vectors.Vector(p'range);
    jm : QuadDobl_Dense_Series_Matrices.Matrix(p'range,x'range);
    ipvt : Standard_Integer_Vectors.Vector(x'range);
    qraux : QuadDobl_Dense_Series_Vectors.Vector(x'range);
    n : constant integer32 := jm'last(1);
    m : constant integer32 := jm'last(2);
    rsd,dum,dum2,dum3 : QuadDobl_Dense_Series_Vectors.Vector(1..n);

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
    QRD(jm,qraux,ipvt,false);
    QRLS(jm,n,m,qraux,px,dum2,dum3,dx,rsd,dum,110,info);
    if info /= 0 then
      put(file,"QRLS info : "); put(file,info,1); new_line(file);
    else
      put_line(file,"The update to the series :");
      for i in dx'range loop
        QuadDobl_Dense_Series_io.put(file,dx(i)); new_line(file);
      end loop;
      QuadDobl_Dense_Series_Vectors.Add(x,dx);
    end if;
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

end QuadDobl_Newton_Series;
