with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;         use Standard_Floating_Numbers;
with Standard_Floating_Numbers_io;      use Standard_Floating_Numbers_io;
with Standard_Integer_Vectors;
with Standard_Dense_Series;
with Standard_Dense_Series_io;
with Standard_Series_Vector_Norms;
with Standard_Dense_Series_Matrices;
with Standard_Linear_Series_Solvers;    use Standard_Linear_Series_Solvers;
with Standard_Least_Squares_Series;     use Standard_Least_Squares_Series;
with Standard_Series_Polynomials;
with Series_and_Polynomials;
with Standard_Series_Poly_SysFun;

package body Standard_Newton_Series is

-- ONE NEWTON STEP :

  procedure LU_Newton_Step
              ( p : in Standard_Series_Poly_Systems.Poly_Sys;
                jp : in Standard_Series_Jaco_Matrices.Jaco_Mat;
                order : in integer32;
                x : in out Standard_Dense_Series_Vectors.Vector;
                info : out integer32 ) is

    dx : Standard_Dense_Series_Vectors.Vector(x'range);
    px : Standard_Dense_Series_Vectors.Vector(p'range);
    jm : Standard_Dense_Series_Matrices.Matrix(p'range,x'range);
    n : constant integer32 := p'last;
    ipvt : Standard_Integer_Vectors.Vector(1..n);

  begin
    px := Standard_Series_Poly_SysFun.Eval(p,x);
    Standard_Dense_Series_Vectors.Min(px);
    Series_and_Polynomials.Set_Order(px,order);
    jm := Standard_Series_Jaco_Matrices.Eval(jp,x);
    Series_and_Polynomials.Set_Order(jm,order);
    LUfac(jm,n,ipvt,info);
    if info = 0 then
      dx := px;
      LUsolve(jm,n,ipvt,dx);
      Standard_Dense_Series_Vectors.Add(x,dx);
    end if;
  end LU_Newton_Step;

  procedure LU_Newton_Step
              ( p : in Standard_Series_Poly_Systems.Poly_Sys;
                order : in integer32;
                x : in out Standard_Dense_Series_Vectors.Vector;
                info : out integer32 ) is

    jp : Standard_Series_Jaco_Matrices.Jaco_Mat(p'range,x'range)
       := Standard_Series_Jaco_Matrices.Create(p);

  begin
    LU_Newton_Step(p,jp,order,x,info);
    Standard_Series_Jaco_Matrices.Clear(jp);
  end LU_Newton_Step;

  procedure LU_Newton_Step
              ( file : in file_type;
                p : in Standard_Series_Poly_Systems.Poly_Sys;
                jp : in Standard_Series_Jaco_Matrices.Jaco_Mat;
                order : in integer32;
                x : in out Standard_Dense_Series_Vectors.Vector;
                info : out integer32 ) is

    dx : Standard_Dense_Series_Vectors.Vector(x'range);
    px : Standard_Dense_Series_Vectors.Vector(p'range);
    jm : Standard_Dense_Series_Matrices.Matrix(p'range,x'range);
    n : constant integer32 := p'last;
    ipvt : Standard_Integer_Vectors.Vector(1..n);

  begin
    px := Standard_Series_Poly_SysFun.Eval(p,x);
    put_line(file,"The evaluated series :");
    for i in px'range loop
      Standard_Dense_Series_io.put(file,px(i)); new_line(file);
    end loop;
    Standard_Dense_Series_Vectors.Min(px);
    Series_and_Polynomials.Set_Order(px,order);
    jm := Standard_Series_Jaco_Matrices.Eval(jp,x);
    Series_and_Polynomials.Set_Order(jm,order);
    LUfac(jm,n,ipvt,info);
    if info /= 0 then
      put(file,"LUfac info : "); put(file,info,1); new_line(file);
    else
      dx := px;
      LUsolve(jm,n,ipvt,dx);
      put_line(file,"The update to the series :");
      for i in dx'range loop
        Standard_Dense_Series_io.put(file,dx(i)); new_line(file);
      end loop;
      Standard_Dense_Series_Vectors.Add(x,dx);
    end if;
  end LU_Newton_Step;

  procedure LU_Newton_Step
              ( file : in file_type;
                p : in Standard_Series_Poly_Systems.Poly_Sys;
                order : in integer32;
                x : in out Standard_Dense_Series_Vectors.Vector;
                info : out integer32 ) is

    jp : Standard_Series_Jaco_Matrices.Jaco_Mat(p'range,x'range)
       := Standard_Series_Jaco_Matrices.Create(p);

  begin
    LU_Newton_Step(file,p,jp,order,x,info);
    Standard_Series_Jaco_Matrices.Clear(jp);
  end LU_Newton_Step;

  procedure QR_Newton_Step
              ( p : in Standard_Series_Poly_Systems.Poly_Sys;
                jp : in Standard_Series_Jaco_Matrices.Jaco_Mat;
                order : in integer32;
                x : in out Standard_Dense_Series_Vectors.Vector;
                info : out integer32 ) is

    wp : Standard_Series_Poly_Systems.Poly_Sys(p'range);
    dx : Standard_Dense_Series_Vectors.Vector(x'range);
    px : Standard_Dense_Series_Vectors.Vector(p'range);
    jm : Standard_Dense_Series_Matrices.Matrix(p'range,x'range);
    ipvt : Standard_Integer_Vectors.Vector(x'range);
    qraux : Standard_Dense_Series_Vectors.Vector(x'range);
    n : constant integer32 := jm'last(1);
    m : constant integer32 := jm'last(2);
    rsd,dum,dum2,dum3 : Standard_Dense_Series_Vectors.Vector(1..n);
    nrm : double_float;
    tol : constant double_float := 1.0E-13;
    wjp : Standard_Series_Jaco_Matrices.Jaco_Mat(jp'range(1),jp'range(2));

  begin
    Standard_Series_Poly_Systems.Copy(p,wp);
    Series_and_Polynomials.Set_Order(wp,order);
    px := Standard_Series_Poly_SysFun.Eval(wp,x);
    nrm := Standard_Series_Vector_Norms.Max_Norm(px);
    if nrm > tol then
      for i in jp'range(1) loop
        for j in jp'range(2) loop
          Standard_Series_Polynomials.Copy(jp(i,j),wjp(i,j));
        end loop;
      end loop;
      Series_and_Polynomials.Set_Order(wjp,order);
      Standard_Dense_Series_Vectors.Min(px);
      Series_and_Polynomials.Set_Order(px,order);
      jm := Standard_Series_Jaco_Matrices.Eval(wjp,x);
      Series_and_Polynomials.Set_Order(jm,order);
      QRD(jm,qraux,ipvt,false);
      QRLS(jm,n,m,qraux,px,dum2,dum3,dx,rsd,dum,110,info);
      if info = 0
       then Standard_Dense_Series_Vectors.Add(x,dx);
      end if;
      Standard_Series_Jaco_Matrices.Clear(wjp);
    end if;
    Standard_Series_Poly_Systems.Clear(wp);
  end QR_Newton_Step;

  procedure QR_Newton_Step
              ( p : in Standard_Series_Poly_Systems.Poly_Sys;
                order : in integer32;
                x : in out Standard_Dense_Series_Vectors.Vector;
                info : out integer32 ) is

    jp : Standard_Series_Jaco_Matrices.Jaco_Mat(p'range,x'range)
       := Standard_Series_Jaco_Matrices.Create(p);

  begin
    QR_Newton_Step(p,jp,order,x,info);
    Standard_Series_Jaco_Matrices.Clear(jp);
  end QR_Newton_Step;

  procedure QR_Newton_Step
              ( file : in file_type;
                p : in Standard_Series_Poly_Systems.Poly_Sys;
                jp : in Standard_Series_Jaco_Matrices.Jaco_Mat;
                order : in integer32;
                x : in out Standard_Dense_Series_Vectors.Vector;
                info : out integer32 ) is

    wp : Standard_Series_Poly_Systems.Poly_Sys(p'range);
    dx : Standard_Dense_Series_Vectors.Vector(x'range);
    px : Standard_Dense_Series_Vectors.Vector(p'range);
    jm : Standard_Dense_Series_Matrices.Matrix(p'range,x'range);
    ipvt : Standard_Integer_Vectors.Vector(x'range);
    qraux : Standard_Dense_Series_Vectors.Vector(x'range);
    n : constant integer32 := jm'last(1);
    m : constant integer32 := jm'last(2);
    rsd,dum,dum2,dum3 : Standard_Dense_Series_Vectors.Vector(1..n);
    nrm : double_float;
    tol : constant double_float := 1.0E-13;
    wjp : Standard_Series_Jaco_Matrices.Jaco_Mat(jm'range(1),jm'range(2));

  begin
    Standard_Series_Poly_Systems.Copy(p,wp);
    Series_and_Polynomials.Set_Order(wp,order);
    px := Standard_Series_Poly_SysFun.Eval(wp,x);
    put_line(file,"The evaluated series :");
    for i in px'range loop
      Standard_Dense_Series_io.put(file,px(i)); new_line(file);
    end loop;
    nrm := Standard_Series_Vector_Norms.Max_Norm(px);
    put(file,"The max norm of the evaluation : ");
    put(file,nrm,3); new_line(file);
    if nrm > tol then
      for i in jp'range(1) loop
        for j in jp'range(2) loop
          Standard_Series_Polynomials.Copy(jp(i,j),wjp(i,j));
        end loop;
      end loop;
      Series_and_Polynomials.Set_Order(wjp,order);
      Standard_Dense_Series_Vectors.Min(px);
      Series_and_Polynomials.Set_Order(px,order);
      jm := Standard_Series_Jaco_Matrices.Eval(jp,x);
      Series_and_Polynomials.Set_Order(jm,order);
      QRD(jm,qraux,ipvt,false);
      QRLS(jm,n,m,qraux,px,dum2,dum3,dx,rsd,dum,110,info);
      if info /= 0 then
        put(file,"QRLS info : "); put(file,info,1); new_line(file);
      else
        put_line(file,"The update to the series :");
        for i in dx'range loop
          Standard_Dense_Series_io.put(file,dx(i)); new_line(file);
        end loop;
        Standard_Dense_Series_Vectors.Add(x,dx);
      end if;
      Standard_Series_Jaco_Matrices.Clear(wjp);
    end if;
    Standard_Series_Poly_Systems.Clear(wp);
  end QR_Newton_Step;

  procedure QR_Newton_Step
              ( file : in file_type;
                p : in Standard_Series_Poly_Systems.Poly_Sys;
                order : in integer32;
                x : in out Standard_Dense_Series_Vectors.Vector;
                info : out integer32 ) is

    jp : Standard_Series_Jaco_Matrices.Jaco_Mat(p'range,x'range)
       := Standard_Series_Jaco_Matrices.Create(p);

  begin
    QR_Newton_Step(file,p,jp,order,x,info);
    Standard_Series_Jaco_Matrices.Clear(jp);
  end QR_Newton_Step;

-- MANY NEWTON STEPS :

  procedure LU_Newton_Steps
              ( p : in Standard_Series_Poly_Systems.Poly_Sys;
                jp : in Standard_Series_Jaco_Matrices.Jaco_Mat;
                order : in out integer32; nbrit : in integer32;
                x : in out Standard_Dense_Series_Vectors.Vector;
                info : out integer32 ) is

  begin
    for i in 1..nbrit loop
      LU_Newton_Step(p,jp,order,x,info);
      exit when (info /= 0); -- stop if Jacobian matrix is singular
      exit when (i = nbrit); -- do not double order after last step
      order := 2*order;
      if order > Standard_Dense_Series.max_order
       then order := Standard_Dense_Series.max_order;
      end if;
    end loop;
  end LU_Newton_Steps;

  procedure LU_Newton_Steps
              ( p : in Standard_Series_Poly_Systems.Poly_Sys;
                order : in out integer32; nbrit : in integer32;
                x : in out Standard_Dense_Series_Vectors.Vector;
                info : out integer32 ) is

    jp : Standard_Series_Jaco_Matrices.Jaco_Mat(p'range,x'range)
       := Standard_Series_Jaco_Matrices.Create(p);

  begin
    LU_Newton_Steps(p,jp,order,nbrit,x,info);
    Standard_Series_Jaco_Matrices.Clear(jp);
  end LU_Newton_Steps;

  procedure LU_Newton_Steps
              ( file : in file_type;
                p : in Standard_Series_Poly_Systems.Poly_Sys;
                jp : in Standard_Series_Jaco_Matrices.Jaco_Mat;
                order : in out integer32; nbrit : in integer32;
                x : in out Standard_Dense_Series_Vectors.Vector;
                info : out integer32 ) is

  begin
    for i in 1..nbrit loop
      put(file,"LU Newton step "); put(file,i,1); put_line(file," :");
      LU_Newton_Step(file,p,jp,order,x,info);
      exit when (info /= 0); -- stop if Jacobian matrix is singular
      exit when (i = nbrit); -- do not double order after last step
      order := 2*order;
      if order > Standard_Dense_Series.max_order
       then order := Standard_Dense_Series.max_order;
      end if;
    end loop;
  end LU_Newton_Steps;

  procedure LU_Newton_Steps
              ( file : in file_type;
                p : in Standard_Series_Poly_Systems.Poly_Sys;
                order : in out integer32; nbrit : in integer32;
                x : in out Standard_Dense_Series_Vectors.Vector;
                info : out integer32 ) is

    jp : Standard_Series_Jaco_Matrices.Jaco_Mat(p'range,x'range)
       := Standard_Series_Jaco_Matrices.Create(p);

  begin
    LU_Newton_Steps(file,p,jp,order,nbrit,x,info);
    Standard_Series_Jaco_Matrices.Clear(jp);
  end LU_Newton_Steps;

  procedure QR_Newton_Steps
              ( p : in Standard_Series_Poly_Systems.Poly_Sys;
                jp : in Standard_Series_Jaco_Matrices.Jaco_Mat;
                order : in out integer32; nbrit : in integer32;
                x : in out Standard_Dense_Series_Vectors.Vector;
                info : out integer32 ) is

  begin
    for i in 1..nbrit loop
      QR_Newton_Step(p,jp,order,x,info);
      exit when (info /= 0); -- stop if Jacobian matrix is singular
      exit when (i = nbrit); -- do not double order after last step
      order := 2*order;
      if order > Standard_Dense_Series.max_order
       then order := Standard_Dense_Series.max_order;
      end if;
    end loop;
  end QR_Newton_Steps;

  procedure QR_Newton_Steps
              ( p : in Standard_Series_Poly_Systems.Poly_Sys;
                order : in out integer32; nbrit : in integer32;
                x : in out Standard_Dense_Series_Vectors.Vector;
                info : out integer32 ) is

    jp : Standard_Series_Jaco_Matrices.Jaco_Mat(p'range,x'range)
       := Standard_Series_Jaco_Matrices.Create(p);

  begin
    QR_Newton_Steps(p,jp,order,nbrit,x,info);
    Standard_Series_Jaco_Matrices.Clear(jp);
  end QR_Newton_Steps;

  procedure QR_Newton_Steps
              ( file : in file_type;
                p : in Standard_Series_Poly_Systems.Poly_Sys;
                jp : in Standard_Series_Jaco_Matrices.Jaco_Mat;
                order : in out integer32; nbrit : in integer32;
                x : in out Standard_Dense_Series_Vectors.Vector;
                info : out integer32 ) is

  begin
    for i in 1..nbrit loop
      put(file,"QR Newton step "); put(file,i,1); put_line(file," :");
      QR_Newton_Step(file,p,jp,order,x,info);
      exit when (info /= 0); -- stop if Jacobian matrix is singular
      exit when (i = nbrit); -- do not double order after last step
      order := 2*order;
      if order > Standard_Dense_Series.max_order
       then order := Standard_Dense_Series.max_order;
      end if;
    end loop;
  end QR_Newton_Steps;

  procedure QR_Newton_Steps
              ( file : in file_type;
                p : in Standard_Series_Poly_Systems.Poly_Sys;
                order : in out integer32; nbrit : in integer32;
                x : in out Standard_Dense_Series_Vectors.Vector;
                info : out integer32 ) is

    jp : Standard_Series_Jaco_Matrices.Jaco_Mat(p'range,x'range)
       := Standard_Series_Jaco_Matrices.Create(p);

  begin
    QR_Newton_Steps(file,p,jp,order,nbrit,x,info);
    Standard_Series_Jaco_Matrices.Clear(jp);
  end QR_Newton_Steps;

end Standard_Newton_Series;
