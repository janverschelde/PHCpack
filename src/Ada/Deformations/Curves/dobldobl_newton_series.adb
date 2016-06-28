with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with Standard_Integer_Vectors;
with DoblDobl_Dense_Series;
with DoblDobl_Dense_Series_io;
with DoblDobl_Dense_Series_Matrices;
with DoblDobl_Linear_Series_Solvers;    use DoblDobl_Linear_Series_Solvers;
with DoblDobl_Least_Squares_Series;     use DoblDobl_Least_Squares_Series;
with Series_and_Polynomials;
with DoblDobl_Series_Poly_SysFun;

package body DoblDobl_Newton_Series is

-- ONE NEWTON STEP :

  procedure LU_Newton_Step
              ( p : in DoblDobl_Series_Poly_Systems.Poly_Sys;
                jp : in DoblDobl_Series_Jaco_Matrices.Jaco_Mat;
                order : in integer32;
                x : in out DoblDobl_Dense_Series_Vectors.Vector;
                info : out integer32 ) is

    dx : DoblDobl_Dense_Series_Vectors.Vector(x'range);
    px : DoblDobl_Dense_Series_Vectors.Vector(p'range);
    jm : DoblDobl_Dense_Series_Matrices.Matrix(p'range,x'range);
    n : constant integer32 := p'last;
    ipvt : Standard_Integer_Vectors.Vector(1..n);

  begin
    px := DoblDobl_Series_Poly_SysFun.Eval(p,x);
    DoblDobl_Dense_Series_Vectors.Min(px);
    Series_and_Polynomials.Set_Order(px,order);
    jm := DoblDobl_Series_Jaco_Matrices.Eval(jp,x);
    Series_and_Polynomials.Set_Order(jm,order);
    LUfac(jm,n,ipvt,info);
    if info = 0 then
      dx := px;
      LUsolve(jm,n,ipvt,dx);
      DoblDobl_Dense_Series_Vectors.Add(x,dx);
    end if;
  end LU_Newton_Step;

  procedure LU_Newton_Step
              ( p : in DoblDobl_Series_Poly_Systems.Poly_Sys;
                order : in integer32;
                x : in out DoblDobl_Dense_Series_Vectors.Vector;
                info : out integer32 ) is

    jp : DoblDobl_Series_Jaco_Matrices.Jaco_Mat(p'range,x'range)
       := DoblDobl_Series_Jaco_Matrices.Create(p);

  begin
    LU_Newton_Step(p,jp,order,x,info);
    DoblDobl_Series_Jaco_Matrices.Clear(jp);
  end LU_Newton_Step;

  procedure LU_Newton_Step
              ( file : in file_type;
                p : in DoblDobl_Series_Poly_Systems.Poly_Sys;
                jp : in DoblDobl_Series_Jaco_Matrices.Jaco_Mat;
                order : in integer32;
                x : in out DoblDobl_Dense_Series_Vectors.Vector;
                info : out integer32 ) is

    dx : DoblDobl_Dense_Series_Vectors.Vector(x'range);
    px : DoblDobl_Dense_Series_Vectors.Vector(p'range);
    jm : DoblDobl_Dense_Series_Matrices.Matrix(p'range,x'range);
    n : constant integer32 := p'last;
    ipvt : Standard_Integer_Vectors.Vector(1..n);

  begin
    px := DoblDobl_Series_Poly_SysFun.Eval(p,x);
    put_line(file,"The evaluated series :");
    for i in px'range loop
      DoblDobl_Dense_Series_io.put(file,px(i)); new_line(file);
    end loop;
    DoblDobl_Dense_Series_Vectors.Min(px);
    Series_and_Polynomials.Set_Order(px,order);
    jm := DoblDobl_Series_Jaco_Matrices.Eval(jp,x);
    Series_and_Polynomials.Set_Order(jm,order);
    LUfac(jm,n,ipvt,info);
    if info /= 0 then
      put(file,"LUfac info : "); put(file,info,1); new_line(file);
    else
      dx := px;
      LUsolve(jm,n,ipvt,dx);
      put_line(file,"The update to the series :");
      for i in dx'range loop
        DoblDobl_Dense_Series_io.put(file,dx(i)); new_line(file);
      end loop;
      DoblDobl_Dense_Series_Vectors.Add(x,dx);
    end if;
  end LU_Newton_Step;

  procedure LU_Newton_Step
              ( file : in file_type;
                p : in DoblDobl_Series_Poly_Systems.Poly_Sys;
                order : in integer32;
                x : in out DoblDobl_Dense_Series_Vectors.Vector;
                info : out integer32 ) is

    jp : DoblDobl_Series_Jaco_Matrices.Jaco_Mat(p'range,x'range)
       := DoblDobl_Series_Jaco_Matrices.Create(p);

  begin
    LU_Newton_Step(file,p,jp,order,x,info);
    DoblDobl_Series_Jaco_Matrices.Clear(jp);
  end LU_Newton_Step;

  procedure QR_Newton_Step
              ( p : in DoblDobl_Series_Poly_Systems.Poly_Sys;
                jp : in DoblDobl_Series_Jaco_Matrices.Jaco_Mat;
                order : in integer32;
                x : in out DoblDobl_Dense_Series_Vectors.Vector;
                info : out integer32 ) is

    dx : DoblDobl_Dense_Series_Vectors.Vector(x'range);
    px : DoblDobl_Dense_Series_Vectors.Vector(p'range);
    jm : DoblDobl_Dense_Series_Matrices.Matrix(p'range,x'range);
    ipvt : Standard_Integer_Vectors.Vector(x'range);
    qraux : DoblDobl_Dense_Series_Vectors.Vector(x'range);
    n : constant integer32 := jm'last(1);
    m : constant integer32 := jm'last(2);
    rsd,dum,dum2,dum3 : DoblDobl_Dense_Series_Vectors.Vector(1..n);

  begin
    px := DoblDobl_Series_Poly_SysFun.Eval(p,x);
    DoblDobl_Dense_Series_Vectors.Min(px);
    Series_and_Polynomials.Set_Order(px,order);
    jm := DoblDobl_Series_Jaco_Matrices.Eval(jp,x);
    Series_and_Polynomials.Set_Order(jm,order);
    QRD(jm,qraux,ipvt,false);
    QRLS(jm,n,m,qraux,px,dum2,dum3,dx,rsd,dum,110,info);
    if info = 0
     then DoblDobl_Dense_Series_Vectors.Add(x,dx);
    end if;
  end QR_Newton_Step;

  procedure QR_Newton_Step
              ( p : in DoblDobl_Series_Poly_Systems.Poly_Sys;
                order : in integer32;
                x : in out DoblDobl_Dense_Series_Vectors.Vector;
                info : out integer32 ) is

    jp : DoblDobl_Series_Jaco_Matrices.Jaco_Mat(p'range,x'range)
       := DoblDobl_Series_Jaco_Matrices.Create(p);

  begin
    QR_Newton_Step(p,jp,order,x,info);
    DoblDobl_Series_Jaco_Matrices.Clear(jp);
  end QR_Newton_Step;

  procedure QR_Newton_Step
              ( file : in file_type;
                p : in DoblDobl_Series_Poly_Systems.Poly_Sys;
                jp : in DoblDobl_Series_Jaco_Matrices.Jaco_Mat;
                order : in integer32;
                x : in out DoblDobl_Dense_Series_Vectors.Vector;
                info : out integer32 ) is

    dx : DoblDobl_Dense_Series_Vectors.Vector(x'range);
    px : DoblDobl_Dense_Series_Vectors.Vector(p'range);
    jm : DoblDobl_Dense_Series_Matrices.Matrix(p'range,x'range);
    ipvt : Standard_Integer_Vectors.Vector(x'range);
    qraux : DoblDobl_Dense_Series_Vectors.Vector(x'range);
    n : constant integer32 := jm'last(1);
    m : constant integer32 := jm'last(2);
    rsd,dum,dum2,dum3 : DoblDobl_Dense_Series_Vectors.Vector(1..n);

  begin
    px := DoblDobl_Series_Poly_SysFun.Eval(p,x);
    DoblDobl_Dense_Series_Vectors.Min(px);
    Series_and_Polynomials.Set_Order(px,order);
    jm := DoblDobl_Series_Jaco_Matrices.Eval(jp,x);
    Series_and_Polynomials.Set_Order(jm,order);
    QRD(jm,qraux,ipvt,false);
    QRLS(jm,n,m,qraux,px,dum2,dum3,dx,rsd,dum,110,info);
    if info /= 0 then
      put(file,"QRLS info : "); put(file,info,1); new_line(file);
    else
      put_line(file,"The update to the series :");
      for i in dx'range loop
        DoblDobl_Dense_Series_io.put(file,dx(i)); new_line(file);
      end loop;
      DoblDobl_Dense_Series_Vectors.Add(x,dx);
    end if;
  end QR_Newton_Step;

  procedure QR_Newton_Step
              ( file : in file_type;
                p : in DoblDobl_Series_Poly_Systems.Poly_Sys;
                order : in integer32;
                x : in out DoblDobl_Dense_Series_Vectors.Vector;
                info : out integer32 ) is

    jp : DoblDobl_Series_Jaco_Matrices.Jaco_Mat(p'range,x'range)
       := DoblDobl_Series_Jaco_Matrices.Create(p);

  begin
    QR_Newton_Step(file,p,jp,order,x,info);
    DoblDobl_Series_Jaco_Matrices.Clear(jp);
  end QR_Newton_Step;

-- MANY NEWTON STEPS :

  procedure LU_Newton_Steps
              ( p : in DoblDobl_Series_Poly_Systems.Poly_Sys;
                jp : in DoblDobl_Series_Jaco_Matrices.Jaco_Mat;
                order : in out integer32; nbrit : in integer32;
                x : in out DoblDobl_Dense_Series_Vectors.Vector;
                info : out integer32 ) is

  begin
    for i in 1..nbrit loop
      LU_Newton_Step(p,jp,order,x,info);
      exit when (info /= 0); -- stop if Jacobian matrix is singular
      exit when (i = nbrit); -- do not double order after last step
      order := 2*order;
      if order > DoblDobl_Dense_Series.max_order
       then order := DoblDobl_Dense_Series.max_order;
      end if;
    end loop;
  end LU_Newton_Steps;

  procedure LU_Newton_Steps
              ( p : in DoblDobl_Series_Poly_Systems.Poly_Sys;
                order : in out integer32; nbrit : in integer32;
                x : in out DoblDobl_Dense_Series_Vectors.Vector;
                info : out integer32 ) is

    jp : DoblDobl_Series_Jaco_Matrices.Jaco_Mat(p'range,x'range)
       := DoblDobl_Series_Jaco_Matrices.Create(p);

  begin
    LU_Newton_Steps(p,jp,order,nbrit,x,info);
    DoblDobl_Series_Jaco_Matrices.Clear(jp);
  end LU_Newton_Steps;

  procedure LU_Newton_Steps
              ( file : in file_type;
                p : in DoblDobl_Series_Poly_Systems.Poly_Sys;
                jp : in DoblDobl_Series_Jaco_Matrices.Jaco_Mat;
                order : in out integer32; nbrit : in integer32;
                x : in out DoblDobl_Dense_Series_Vectors.Vector;
                info : out integer32 ) is

  begin
    for i in 1..nbrit loop
      put(file,"LU Newton step "); put(file,i,1); put_line(file," :");
      LU_Newton_Step(file,p,jp,order,x,info);
      exit when (info /= 0); -- stop if Jacobian matrix is singular
      exit when (i = nbrit); -- do not double order after last step
      order := 2*order;
      if order > DoblDobl_Dense_Series.max_order
       then order := DoblDobl_Dense_Series.max_order;
      end if;
    end loop;
  end LU_Newton_Steps;

  procedure LU_Newton_Steps
              ( file : in file_type;
                p : in DoblDobl_Series_Poly_Systems.Poly_Sys;
                order : in out integer32; nbrit : in integer32;
                x : in out DoblDobl_Dense_Series_Vectors.Vector;
                info : out integer32 ) is

    jp : DoblDobl_Series_Jaco_Matrices.Jaco_Mat(p'range,x'range)
       := DoblDobl_Series_Jaco_Matrices.Create(p);

  begin
    LU_Newton_Steps(file,p,jp,order,nbrit,x,info);
    DoblDobl_Series_Jaco_Matrices.Clear(jp);
  end LU_Newton_Steps;

  procedure QR_Newton_Steps
              ( p : in DoblDobl_Series_Poly_Systems.Poly_Sys;
                jp : in DoblDobl_Series_Jaco_Matrices.Jaco_Mat;
                order : in out integer32; nbrit : in integer32;
                x : in out DoblDobl_Dense_Series_Vectors.Vector;
                info : out integer32 ) is

  begin
    for i in 1..nbrit loop
      QR_Newton_Step(p,jp,order,x,info);
      exit when (info /= 0); -- stop if Jacobian matrix is singular
      exit when (i = nbrit); -- do not double order after last step
      order := 2*order;
      if order > DoblDobl_Dense_Series.max_order
       then order := DoblDobl_Dense_Series.max_order;
      end if;
    end loop;
  end QR_Newton_Steps;

  procedure QR_Newton_Steps
              ( p : in DoblDobl_Series_Poly_Systems.Poly_Sys;
                order : in out integer32; nbrit : in integer32;
                x : in out DoblDobl_Dense_Series_Vectors.Vector;
                info : out integer32 ) is

    jp : DoblDobl_Series_Jaco_Matrices.Jaco_Mat(p'range,x'range)
       := DoblDobl_Series_Jaco_Matrices.Create(p);

  begin
    QR_Newton_Steps(p,jp,order,nbrit,x,info);
    DoblDobl_Series_Jaco_Matrices.Clear(jp);
  end QR_Newton_Steps;

  procedure QR_Newton_Steps
              ( file : in file_type;
                p : in DoblDobl_Series_Poly_Systems.Poly_Sys;
                jp : in DoblDobl_Series_Jaco_Matrices.Jaco_Mat;
                order : in out integer32; nbrit : in integer32;
                x : in out DoblDobl_Dense_Series_Vectors.Vector;
                info : out integer32 ) is

  begin
    for i in 1..nbrit loop
      put(file,"QR Newton step "); put(file,i,1); put_line(file," :");
      QR_Newton_Step(file,p,jp,order,x,info);
      exit when (info /= 0); -- stop if Jacobian matrix is singular
      exit when (i = nbrit); -- do not double order after last step
      order := 2*order;
      if order > DoblDobl_Dense_Series.max_order
       then order := DoblDobl_Dense_Series.max_order;
      end if;
    end loop;
  end QR_Newton_Steps;

  procedure QR_Newton_Steps
              ( file : in file_type;
                p : in DoblDobl_Series_Poly_Systems.Poly_Sys;
                order : in out integer32; nbrit : in integer32;
                x : in out DoblDobl_Dense_Series_Vectors.Vector;
                info : out integer32 ) is

    jp : DoblDobl_Series_Jaco_Matrices.Jaco_Mat(p'range,x'range)
       := DoblDobl_Series_Jaco_Matrices.Create(p);

  begin
    QR_Newton_Steps(file,p,jp,order,nbrit,x,info);
    DoblDobl_Series_Jaco_Matrices.Clear(jp);
  end QR_Newton_Steps;

end DoblDobl_Newton_Series;
