with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;         use Standard_Floating_Numbers;
with Standard_Floating_Numbers_io;      use Standard_Floating_Numbers_io;
with Standard_Integer_Vectors;
with Standard_Complex_Series_io;
with Standard_CSeries_Vector_Norms;
with Standard_Complex_Series_Matrices;
with Standard_Series_Linear_Solvers;    use Standard_Series_Linear_Solvers;
with Standard_Series_Least_Squares;     use Standard_Series_Least_Squares;
with Standard_CSeries_Polynomials;
with Complex_Series_and_Polynomials;
with Standard_CSeries_Poly_SysFun;

package body Standard_Newton_Series is

-- ONE NEWTON STEP :

  procedure LU_Newton_Step
              ( p : in Standard_CSeries_Poly_Systems.Poly_Sys;
                jp : in Standard_CSeries_Jaco_Matrices.Jaco_Mat;
                degree : in integer32;
                x : in out Standard_Complex_Series_Vectors.Vector;
                info : out integer32 ) is

    dx : Standard_Complex_Series_Vectors.Vector(x'range);
    px : Standard_Complex_Series_Vectors.Vector(p'range);
    jm : Standard_Complex_Series_Matrices.Matrix(p'range,x'range);
    n : constant integer32 := p'last;
    ipvt : Standard_Integer_Vectors.Vector(1..n);

  begin
    px := Standard_CSeries_Poly_SysFun.Eval(p,x);
    Standard_Complex_Series_Vectors.Min(px);
    Complex_Series_and_Polynomials.Set_Degree(px,degree);
    jm := Standard_CSeries_Jaco_Matrices.Eval(jp,x);
    Complex_Series_and_Polynomials.Set_Degree(jm,degree);
    LUfac(jm,n,ipvt,info);
    if info = 0 then
      dx := px;
      LUsolve(jm,n,ipvt,dx);
      Standard_Complex_Series_Vectors.Add(x,dx);
    end if;
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
    n : constant integer32 := p'last;
    ipvt : Standard_Integer_Vectors.Vector(1..n);

  begin
    px := Standard_CSeries_Poly_SysFun.Eval(p,x);
    put_line(file,"The evaluated series :");
    for i in px'range loop
      Standard_Complex_Series_io.put(file,px(i)); new_line(file);
    end loop;
    Standard_Complex_Series_Vectors.Min(px);
    Complex_Series_and_Polynomials.Set_Degree(px,degree);
    jm := Standard_CSeries_Jaco_Matrices.Eval(jp,x);
    Complex_Series_and_Polynomials.Set_Degree(jm,degree);
    LUfac(jm,n,ipvt,info);
    if info /= 0 then
      put(file,"LUfac info : "); put(file,info,1); new_line(file);
    else
      dx := px;
      LUsolve(jm,n,ipvt,dx);
      put_line(file,"The update to the series :");
      for i in dx'range loop
        Standard_Complex_Series_io.put(file,dx(i)); new_line(file);
      end loop;
      Standard_Complex_Series_Vectors.Add(x,dx);
    end if;
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
    ipvt : Standard_Integer_Vectors.Vector(x'range);
    qraux : Standard_Complex_Series_Vectors.Vector(x'range);
    n : constant integer32 := jm'last(1);
    m : constant integer32 := jm'last(2);
    rsd,dum,dum2,dum3 : Standard_Complex_Series_Vectors.Vector(1..n);
    nrm : double_float;
    tol : constant double_float := 1.0E-13;
    wjp : Standard_CSeries_Jaco_Matrices.Jaco_Mat(jp'range(1),jp'range(2));

  begin
    Standard_CSeries_Poly_Systems.Copy(p,wp);
    Complex_Series_and_Polynomials.Set_Degree(wp,degree);
    px := Standard_CSeries_Poly_SysFun.Eval(wp,x);
    nrm := Standard_CSeries_Vector_Norms.Max_Norm(px);
    if nrm > tol then
      for i in jp'range(1) loop
        for j in jp'range(2) loop
          Standard_CSeries_Polynomials.Copy(jp(i,j),wjp(i,j));
        end loop;
      end loop;
      Complex_Series_and_Polynomials.Set_Degree(wjp,degree);
      Standard_Complex_Series_Vectors.Min(px);
      Complex_Series_and_Polynomials.Set_Degree(px,degree);
      jm := Standard_CSeries_Jaco_Matrices.Eval(wjp,x);
      Complex_Series_and_Polynomials.Set_Degree(jm,degree);
      QRD(jm,qraux,ipvt,false);
      QRLS(jm,n,m,qraux,px,dum2,dum3,dx,rsd,dum,110,info);
      if info = 0
       then Standard_Complex_Series_Vectors.Add(x,dx);
      end if;
      Standard_CSeries_Jaco_Matrices.Clear(wjp);
    end if;
    Standard_CSeries_Poly_Systems.Clear(wp);
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
    ipvt : Standard_Integer_Vectors.Vector(x'range);
    qraux : Standard_Complex_Series_Vectors.Vector(x'range);
    n : constant integer32 := jm'last(1);
    m : constant integer32 := jm'last(2);
    rsd,dum,dum2,dum3 : Standard_Complex_Series_Vectors.Vector(1..n);
    nrm : double_float;
    tol : constant double_float := 1.0E-13;
    wjp : Standard_CSeries_Jaco_Matrices.Jaco_Mat(jm'range(1),jm'range(2));

  begin
    Standard_CSeries_Poly_Systems.Copy(p,wp);
    Complex_Series_and_Polynomials.Set_Degree(wp,degree);
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
      Complex_Series_and_Polynomials.Set_Degree(wjp,degree);
      Standard_Complex_Series_Vectors.Min(px);
      Complex_Series_and_Polynomials.Set_Degree(px,degree);
      jm := Standard_CSeries_Jaco_Matrices.Eval(jp,x);
      Complex_Series_and_Polynomials.Set_Degree(jm,degree);
      QRD(jm,qraux,ipvt,false);
      QRLS(jm,n,m,qraux,px,dum2,dum3,dx,rsd,dum,110,info);
      if info /= 0 then
        put(file,"QRLS info : "); put(file,info,1); new_line(file);
      else
        put_line(file,"The update to the series :");
        for i in dx'range loop
          Standard_Complex_Series_io.put(file,dx(i)); new_line(file);
        end loop;
        Standard_Complex_Series_Vectors.Add(x,dx);
      end if;
      Standard_CSeries_Jaco_Matrices.Clear(wjp);
    end if;
    Standard_CSeries_Poly_Systems.Clear(wp);
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

-- MANY NEWTON STEPS :

  procedure LU_Newton_Steps
              ( p : in Standard_CSeries_Poly_Systems.Poly_Sys;
                jp : in Standard_CSeries_Jaco_Matrices.Jaco_Mat;
                degree : in out integer32; maxdeg,nbrit : in integer32;
                x : in out Standard_Complex_Series_Vectors.Vector;
                info : out integer32 ) is

  begin
    for i in 1..nbrit loop
      LU_Newton_Step(p,jp,degree,x,info);
      exit when (info /= 0); -- stop if Jacobian matrix is singular
      exit when (i = nbrit); -- do not double degree after last step
      if degree < maxdeg then
        degree := 2*degree;
        if degree > maxdeg
         then degree := maxdeg;
        end if;
      end if;
    end loop;
  end LU_Newton_Steps;

  procedure LU_Newton_Steps
              ( p : in Standard_CSeries_Poly_Systems.Poly_Sys;
                degree : in out integer32; maxdeg,nbrit : in integer32;
                x : in out Standard_Complex_Series_Vectors.Vector;
                info : out integer32 ) is

    jp : Standard_CSeries_Jaco_Matrices.Jaco_Mat(p'range,x'range)
       := Standard_CSeries_Jaco_Matrices.Create(p);

  begin
    LU_Newton_Steps(p,jp,degree,maxdeg,nbrit,x,info);
    Standard_CSeries_Jaco_Matrices.Clear(jp);
  end LU_Newton_Steps;

  procedure LU_Newton_Steps
              ( file : in file_type;
                p : in Standard_CSeries_Poly_Systems.Poly_Sys;
                jp : in Standard_CSeries_Jaco_Matrices.Jaco_Mat;
                degree : in out integer32; maxdeg,nbrit : in integer32;
                x : in out Standard_Complex_Series_Vectors.Vector;
                info : out integer32 ) is

  begin
    for i in 1..nbrit loop
      put(file,"LU Newton step "); put(file,i,1); put_line(file," :");
      LU_Newton_Step(file,p,jp,degree,x,info);
      exit when (info /= 0); -- stop if Jacobian matrix is singular
      exit when (i = nbrit); -- do not double degree after last step
      if degree < maxdeg then
        degree := 2*degree;
        if degree > maxdeg
         then degree := maxdeg;
        end if;
      end if;
    end loop;
  end LU_Newton_Steps;

  procedure LU_Newton_Steps
              ( file : in file_type;
                p : in Standard_CSeries_Poly_Systems.Poly_Sys;
                degree : in out integer32; maxdeg,nbrit : in integer32;
                x : in out Standard_Complex_Series_Vectors.Vector;
                info : out integer32 ) is

    jp : Standard_CSeries_Jaco_Matrices.Jaco_Mat(p'range,x'range)
       := Standard_CSeries_Jaco_Matrices.Create(p);

  begin
    LU_Newton_Steps(file,p,jp,degree,maxdeg,nbrit,x,info);
    Standard_CSeries_Jaco_Matrices.Clear(jp);
  end LU_Newton_Steps;

  procedure QR_Newton_Steps
              ( p : in Standard_CSeries_Poly_Systems.Poly_Sys;
                jp : in Standard_CSeries_Jaco_Matrices.Jaco_Mat;
                degree : in out integer32; maxdeg,nbrit : in integer32;
                x : in out Standard_Complex_Series_Vectors.Vector;
                info : out integer32 ) is

  begin
    for i in 1..nbrit loop
      QR_Newton_Step(p,jp,degree,x,info);
      exit when (info /= 0); -- stop if Jacobian matrix is singular
      exit when (i = nbrit); -- do not double degree after last step
      if degree < maxdeg then
        degree := 2*degree;
        if degree > maxdeg
         then degree := maxdeg;
        end if;
      end if;
    end loop;
  end QR_Newton_Steps;

  procedure QR_Newton_Steps
              ( p : in Standard_CSeries_Poly_Systems.Poly_Sys;
                degree : in out integer32; maxdeg,nbrit : in integer32;
                x : in out Standard_Complex_Series_Vectors.Vector;
                info : out integer32 ) is

    jp : Standard_CSeries_Jaco_Matrices.Jaco_Mat(p'range,x'range)
       := Standard_CSeries_Jaco_Matrices.Create(p);

  begin
    QR_Newton_Steps(p,jp,degree,maxdeg,nbrit,x,info);
    Standard_CSeries_Jaco_Matrices.Clear(jp);
  end QR_Newton_Steps;

  procedure QR_Newton_Steps
              ( file : in file_type;
                p : in Standard_CSeries_Poly_Systems.Poly_Sys;
                jp : in Standard_CSeries_Jaco_Matrices.Jaco_Mat;
                degree : in out integer32; maxdeg,nbrit : in integer32;
                x : in out Standard_Complex_Series_Vectors.Vector;
                info : out integer32 ) is

  begin
    for i in 1..nbrit loop
      put(file,"QR Newton step "); put(file,i,1); put_line(file," :");
      QR_Newton_Step(file,p,jp,degree,x,info);
      exit when (info /= 0); -- stop if Jacobian matrix is singular
      exit when (i = nbrit); -- do not double degree after last step
      if degree < maxdeg then
        degree := 2*degree;
        if degree > maxdeg
         then degree := maxdeg;
        end if;
      end if;
    end loop;
  end QR_Newton_Steps;

  procedure QR_Newton_Steps
              ( file : in file_type;
                p : in Standard_CSeries_Poly_Systems.Poly_Sys;
                degree : in out integer32; maxdeg,nbrit : in integer32;
                x : in out Standard_Complex_Series_Vectors.Vector;
                info : out integer32 ) is

    jp : Standard_CSeries_Jaco_Matrices.Jaco_Mat(p'range,x'range)
       := Standard_CSeries_Jaco_Matrices.Create(p);

  begin
    QR_Newton_Steps(file,p,jp,degree,maxdeg,nbrit,x,info);
    Standard_CSeries_Jaco_Matrices.Clear(jp);
  end QR_Newton_Steps;

end Standard_Newton_Series;
