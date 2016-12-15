with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;         use Standard_Floating_Numbers;
with Double_Double_Numbers_io;          use Double_Double_Numbers_io;
with DoblDobl_Complex_Numbers_io;       use DoblDobl_Complex_Numbers_io;
with DoblDobl_Dense_Series;
with DoblDobl_Dense_Series_io;
with DoblDobl_Series_Vector_Norms;
with DoblDobl_Dense_Series_Matrices;
with DoblDobl_Dense_Vector_Series;
with DoblDobl_Dense_Matrix_Series;
with DoblDobl_Matrix_Series_Solvers;    use DoblDobl_Matrix_Series_Solvers;
with DoblDobl_Series_Polynomials;
with Series_and_Polynomials;
with DoblDobl_Series_Poly_SysFun;

package body DoblDobl_Newton_Matrix_Series is

-- ONE NEWTON STEP WITH LU WITHOUT CONDITION NUMBER ESTIMATE :

  procedure LU_Newton_Step
              ( p : in DoblDobl_Series_Poly_Systems.Poly_Sys;
                jp : in DoblDobl_Series_Jaco_Matrices.Jaco_Mat;
                degree : in integer32;
                x : in out DoblDobl_Dense_Series_Vectors.Vector;
                info : out integer32 ) is

    dx : DoblDobl_Dense_Series_Vectors.Vector(x'range);
    px : DoblDobl_Dense_Series_Vectors.Vector(p'range);
    jm : DoblDobl_Dense_Series_Matrices.Matrix(p'range,x'range);
    xp,xd : DoblDobl_Dense_Vector_Series.Vector;
    mj : DoblDobl_Dense_Matrix_Series.Matrix;

  begin
    px := DoblDobl_Series_Poly_SysFun.Eval(p,x);
    DoblDobl_Dense_Series_Vectors.Min(px);
    Series_and_Polynomials.Set_degree(px,degree);
    jm := DoblDobl_Series_Jaco_Matrices.Eval(jp,x);
    Series_and_Polynomials.Set_degree(jm,degree);
    mj := DoblDobl_Dense_Matrix_Series.Create(jm);
    xp := DoblDobl_Dense_Vector_Series.Create(px);
    Solve_by_lufac(mj,xp,info,xd);
    if info = 0 then
      dx := DoblDobl_Dense_Vector_Series.Create(xd);
      DoblDobl_Dense_Series_Vectors.Add(x,dx);
    end if;
    DoblDobl_Dense_Matrix_Series.Clear(mj);
    DoblDobl_Dense_Vector_Series.Clear(xp);
    DoblDobl_Dense_Vector_Series.Clear(xd);
  end LU_Newton_Step;

  procedure LU_Newton_Step
              ( p : in DoblDobl_Series_Poly_Systems.Poly_Sys;
                degree : in integer32;
                x : in out DoblDobl_Dense_Series_Vectors.Vector;
                info : out integer32 ) is

    jp : DoblDobl_Series_Jaco_Matrices.Jaco_Mat(p'range,x'range)
       := DoblDobl_Series_Jaco_Matrices.Create(p);

  begin
    LU_Newton_Step(p,jp,degree,x,info);
    DoblDobl_Series_Jaco_Matrices.Clear(jp);
  end LU_Newton_Step;

  procedure LU_Newton_Step
              ( file : in file_type;
                p : in DoblDobl_Series_Poly_Systems.Poly_Sys;
                jp : in DoblDobl_Series_Jaco_Matrices.Jaco_Mat;
                degree : in integer32;
                x : in out DoblDobl_Dense_Series_Vectors.Vector;
                info : out integer32 ) is

    dx : DoblDobl_Dense_Series_Vectors.Vector(x'range);
    px : DoblDobl_Dense_Series_Vectors.Vector(p'range);
    jm : DoblDobl_Dense_Series_Matrices.Matrix(p'range,x'range);
    xp,xd : DoblDobl_Dense_Vector_Series.Vector;
    mj : DoblDobl_Dense_Matrix_Series.Matrix;

  begin
    px := DoblDobl_Series_Poly_SysFun.Eval(p,x);
    put_line(file,"The evaluated series :");
    for i in px'range loop
      DoblDobl_Dense_Series_io.put(file,px(i)); new_line(file);
    end loop;
    DoblDobl_Dense_Series_Vectors.Min(px);
    Series_and_Polynomials.Set_degree(px,degree);
    jm := DoblDobl_Series_Jaco_Matrices.Eval(jp,x);
    Series_and_Polynomials.Set_degree(jm,degree);
    mj := DoblDobl_Dense_Matrix_Series.Create(jm);
    xp := DoblDobl_Dense_Vector_Series.Create(px);
    Solve_by_lufac(mj,xp,info,xd);
    if info /= 0 then
      put(file,"LUfac info : "); put(file,info,1); new_line(file);
    else
      dx := DoblDobl_Dense_Vector_Series.Create(xd);
      put_line(file,"The update to the series :");
      for i in dx'range loop
        DoblDobl_Dense_Series_io.put(file,dx(i)); new_line(file);
      end loop;
      DoblDobl_Dense_Series_Vectors.Add(x,dx);
    end if;
    DoblDobl_Dense_Matrix_Series.Clear(mj);
    DoblDobl_Dense_Vector_Series.Clear(xp);
    DoblDobl_Dense_Vector_Series.Clear(xd);
  end LU_Newton_Step;

  procedure LU_Newton_Step
              ( file : in file_type;
                p : in DoblDobl_Series_Poly_Systems.Poly_Sys;
                degree : in integer32;
                x : in out DoblDobl_Dense_Series_Vectors.Vector;
                info : out integer32 ) is

    jp : DoblDobl_Series_Jaco_Matrices.Jaco_Mat(p'range,x'range)
       := DoblDobl_Series_Jaco_Matrices.Create(p);

  begin
    LU_Newton_Step(file,p,jp,degree,x,info);
    DoblDobl_Series_Jaco_Matrices.Clear(jp);
  end LU_Newton_Step;

-- ONE NEWTON STEP WITH LU WITH CONDITION NUMBER ESTIMATE :

  procedure LU_Newton_Step
              ( p : in DoblDobl_Series_Poly_Systems.Poly_Sys;
                jp : in DoblDobl_Series_Jaco_Matrices.Jaco_Mat;
                degree : in integer32;
                x : in out DoblDobl_Dense_Series_Vectors.Vector;
                rcond : out double_double ) is

    dx : DoblDobl_Dense_Series_Vectors.Vector(x'range);
    px : DoblDobl_Dense_Series_Vectors.Vector(p'range);
    jm : DoblDobl_Dense_Series_Matrices.Matrix(p'range,x'range);
    xp,xd : DoblDobl_Dense_Vector_Series.Vector;
    mj : DoblDobl_Dense_Matrix_Series.Matrix;
    one : constant double_double := create(1.0);

  begin
    px := DoblDobl_Series_Poly_SysFun.Eval(p,x);
    DoblDobl_Dense_Series_Vectors.Min(px);
    Series_and_Polynomials.Set_degree(px,degree);
    jm := DoblDobl_Series_Jaco_Matrices.Eval(jp,x);
    Series_and_Polynomials.Set_degree(jm,degree);
    mj := DoblDobl_Dense_Matrix_Series.Create(jm);
    xp := DoblDobl_Dense_Vector_Series.Create(px);
    Solve_by_lufco(mj,xp,rcond,xd);
    if one + rcond /= one then
      dx := DoblDobl_Dense_Vector_Series.Create(xd);
      DoblDobl_Dense_Series_Vectors.Add(x,dx);
    end if;
    DoblDobl_Dense_Matrix_Series.Clear(mj);
    DoblDobl_Dense_Vector_Series.Clear(xp);
    DoblDobl_Dense_Vector_Series.Clear(xd);
  end LU_Newton_Step;

  procedure LU_Newton_Step
              ( p : in DoblDobl_Series_Poly_Systems.Poly_Sys;
                degree : in integer32;
                x : in out DoblDobl_Dense_Series_Vectors.Vector;
                rcond : out double_double ) is

    jp : DoblDobl_Series_Jaco_Matrices.Jaco_Mat(p'range,x'range)
       := DoblDobl_Series_Jaco_Matrices.Create(p);

  begin
    LU_Newton_Step(p,jp,degree,x,rcond);
    DoblDobl_Series_Jaco_Matrices.Clear(jp);
  end LU_Newton_Step;

  procedure LU_Newton_Step
              ( file : in file_type;
                p : in DoblDobl_Series_Poly_Systems.Poly_Sys;
                jp : in DoblDobl_Series_Jaco_Matrices.Jaco_Mat;
                degree : in integer32;
                x : in out DoblDobl_Dense_Series_Vectors.Vector;
                rcond : out double_double ) is

    dx : DoblDobl_Dense_Series_Vectors.Vector(x'range);
    px : DoblDobl_Dense_Series_Vectors.Vector(p'range);
    jm : DoblDobl_Dense_Series_Matrices.Matrix(p'range,x'range);
    xp,xd : DoblDobl_Dense_Vector_Series.Vector;
    mj : DoblDobl_Dense_Matrix_Series.Matrix;
    one : constant double_double := create(1.0);

  begin
    px := DoblDobl_Series_Poly_SysFun.Eval(p,x);
    put_line(file,"The evaluated series :");
    for i in px'range loop
      DoblDobl_Dense_Series_io.put(file,px(i)); new_line(file);
    end loop;
    DoblDobl_Dense_Series_Vectors.Min(px);
    Series_and_Polynomials.Set_degree(px,degree);
    jm := DoblDobl_Series_Jaco_Matrices.Eval(jp,x);
    Series_and_Polynomials.Set_degree(jm,degree);
    mj := DoblDobl_Dense_Matrix_Series.Create(jm);
    xp := DoblDobl_Dense_Vector_Series.Create(px);
    Solve_by_lufco(mj,xp,rcond,xd);
    put(file,"LUfco rcond : "); put(file,rcond,3); new_line(file);
    if one + rcond /= one then
      dx := DoblDobl_Dense_Vector_Series.Create(xd);
      put_line(file,"The update to the series :");
      for i in dx'range loop
        DoblDobl_Dense_Series_io.put(file,dx(i)); new_line(file);
      end loop;
      DoblDobl_Dense_Series_Vectors.Add(x,dx);
    end if;
    DoblDobl_Dense_Matrix_Series.Clear(mj);
    DoblDobl_Dense_Vector_Series.Clear(xp);
    DoblDobl_Dense_Vector_Series.Clear(xd);
  end LU_Newton_Step;

  procedure LU_Newton_Step
              ( file : in file_type;
                p : in DoblDobl_Series_Poly_Systems.Poly_Sys;
                degree : in integer32;
                x : in out DoblDobl_Dense_Series_Vectors.Vector;
                rcond : out double_double ) is

    jp : DoblDobl_Series_Jaco_Matrices.Jaco_Mat(p'range,x'range)
       := DoblDobl_Series_Jaco_Matrices.Create(p);

  begin
    LU_Newton_Step(file,p,jp,degree,x,rcond);
    DoblDobl_Series_Jaco_Matrices.Clear(jp);
  end LU_Newton_Step;

-- ONE NEWTON STEP WITH QR :

  procedure QR_Newton_Step
              ( p : in DoblDobl_Series_Poly_Systems.Poly_Sys;
                jp : in DoblDobl_Series_Jaco_Matrices.Jaco_Mat;
                degree : in integer32;
                x : in out DoblDobl_Dense_Series_Vectors.Vector;
                info : out integer32 ) is

    wp : DoblDobl_Series_Poly_Systems.Poly_Sys(p'range);
    dx : DoblDobl_Dense_Series_Vectors.Vector(x'range);
    px : DoblDobl_Dense_Series_Vectors.Vector(p'range);
    jm : DoblDobl_Dense_Series_Matrices.Matrix(p'range,x'range);
    mj : DoblDobl_Dense_Matrix_Series.Matrix;
    xp,xd : DoblDobl_Dense_Vector_Series.Vector;
    nrm : double_double;
    tol : constant double_float := 1.0E-26;
    wjp : DoblDobl_Series_Jaco_Matrices.Jaco_Mat(jp'range(1),jp'range(2));

  begin
    DoblDobl_Series_Poly_Systems.Copy(p,wp);
    Series_and_Polynomials.Set_degree(wp,degree);
    px := DoblDobl_Series_Poly_SysFun.Eval(wp,x);
    nrm := DoblDobl_Series_Vector_Norms.Max_Norm(px);
    if nrm > tol then
      for i in jp'range(1) loop
        for j in jp'range(2) loop
          DoblDobl_Series_Polynomials.Copy(jp(i,j),wjp(i,j));
        end loop;
      end loop;
      Series_and_Polynomials.Set_degree(wjp,degree);
      DoblDobl_Dense_Series_Vectors.Min(px);
      Series_and_Polynomials.Set_degree(px,degree);
      jm := DoblDobl_Series_Jaco_Matrices.Eval(wjp,x);
      Series_and_Polynomials.Set_degree(jm,degree);
      mj := DoblDobl_Dense_Matrix_Series.Create(jm);
      xp := DoblDobl_Dense_Vector_Series.Create(px);
      Solve_by_QRLS(mj,xp,info,xd);
      if info = 0 then
        dx := DoblDobl_Dense_Vector_Series.Create(xd);
        DoblDobl_Dense_Series_Vectors.Add(x,dx);
      end if;
      DoblDobl_Series_Jaco_Matrices.Clear(wjp);
    end if;
    DoblDobl_Series_Poly_Systems.Clear(wp);
    DoblDobl_Dense_Matrix_Series.Clear(mj);
    DoblDobl_Dense_Vector_Series.Clear(xp);
    DoblDobl_Dense_Vector_Series.Clear(xd);
  end QR_Newton_Step;

  procedure QR_Newton_Step
              ( p : in DoblDobl_Series_Poly_Systems.Poly_Sys;
                degree : in integer32;
                x : in out DoblDobl_Dense_Series_Vectors.Vector;
                info : out integer32 ) is

    jp : DoblDobl_Series_Jaco_Matrices.Jaco_Mat(p'range,x'range)
       := DoblDobl_Series_Jaco_Matrices.Create(p);

  begin
    QR_Newton_Step(p,jp,degree,x,info);
    DoblDobl_Series_Jaco_Matrices.Clear(jp);
  end QR_Newton_Step;

  procedure QR_Newton_Step
              ( file : in file_type;
                p : in DoblDobl_Series_Poly_Systems.Poly_Sys;
                jp : in DoblDobl_Series_Jaco_Matrices.Jaco_Mat;
                degree : in integer32;
                x : in out DoblDobl_Dense_Series_Vectors.Vector;
                info : out integer32 ) is

    wp : DoblDobl_Series_Poly_Systems.Poly_Sys(p'range);
    dx : DoblDobl_Dense_Series_Vectors.Vector(x'range);
    px : DoblDobl_Dense_Series_Vectors.Vector(p'range);
    jm : DoblDobl_Dense_Series_Matrices.Matrix(p'range,x'range);
    mj : DoblDobl_Dense_Matrix_Series.Matrix;
    xp,xd : DoblDobl_Dense_Vector_Series.Vector;
    nrm : double_double;
    tol : constant double_float := 1.0E-26;
    wjp : DoblDobl_Series_Jaco_Matrices.Jaco_Mat(jm'range(1),jm'range(2));

  begin
    DoblDobl_Series_Poly_Systems.Copy(p,wp);
    Series_and_Polynomials.Set_degree(wp,degree);
    px := DoblDobl_Series_Poly_SysFun.Eval(wp,x);
    put_line(file,"The evaluated series :");
    for i in px'range loop
      DoblDobl_Dense_Series_io.put(file,px(i)); new_line(file);
    end loop;
    nrm := DoblDobl_Series_Vector_Norms.Max_Norm(px);
    put(file,"The max norm of the evaluation : ");
    put(file,nrm,3); new_line(file);
    if nrm > tol then
      for i in jp'range(1) loop
        for j in jp'range(2) loop
          DoblDobl_Series_Polynomials.Copy(jp(i,j),wjp(i,j));
        end loop;
      end loop;
      Series_and_Polynomials.Set_degree(wjp,degree);
      DoblDobl_Dense_Series_Vectors.Min(px);
      Series_and_Polynomials.Set_degree(px,degree);
      jm := DoblDobl_Series_Jaco_Matrices.Eval(jp,x);
      Series_and_Polynomials.Set_degree(jm,degree);
      mj := DoblDobl_Dense_Matrix_Series.Create(jm);
      xp := DoblDobl_Dense_Vector_Series.Create(px);
      Solve_by_QRLS(mj,xp,info,xd);
      if info /= 0 then
        put(file,"QRLS info : "); put(file,info,1); new_line(file);
      else
        dx := DoblDobl_Dense_Vector_Series.Create(xd);
        put_line(file,"The update to the series :");
        for i in dx'range loop
          DoblDobl_Dense_Series_io.put(file,dx(i)); new_line(file);
        end loop;
        DoblDobl_Dense_Series_Vectors.Add(x,dx);
      end if;
      DoblDobl_Series_Jaco_Matrices.Clear(wjp);
    end if;
    DoblDobl_Series_Poly_Systems.Clear(wp);
    DoblDobl_Dense_Matrix_Series.Clear(mj);
    DoblDobl_Dense_Vector_Series.Clear(xp);
    DoblDobl_Dense_Vector_Series.Clear(xd);
  end QR_Newton_Step;

  procedure QR_Newton_Step
              ( file : in file_type;
                p : in DoblDobl_Series_Poly_Systems.Poly_Sys;
                degree : in integer32;
                x : in out DoblDobl_Dense_Series_Vectors.Vector;
                info : out integer32 ) is

    jp : DoblDobl_Series_Jaco_Matrices.Jaco_Mat(p'range,x'range)
       := DoblDobl_Series_Jaco_Matrices.Create(p);

  begin
    QR_Newton_Step(file,p,jp,degree,x,info);
    DoblDobl_Series_Jaco_Matrices.Clear(jp);
  end QR_Newton_Step;

-- ONE NEWTON STEP WITH SVD :

  procedure SVD_Newton_Step
              ( p : in DoblDobl_Series_Poly_Systems.Poly_Sys;
                jp : in DoblDobl_Series_Jaco_Matrices.Jaco_Mat;
                degree : in integer32;
                x : in out DoblDobl_Dense_Series_Vectors.Vector;
                info : out integer32; rcond : out double_double ) is

    wp : DoblDobl_Series_Poly_Systems.Poly_Sys(p'range);
    dx : DoblDobl_Dense_Series_Vectors.Vector(x'range);
    px : DoblDobl_Dense_Series_Vectors.Vector(p'range);
    jm : DoblDobl_Dense_Series_Matrices.Matrix(p'range,x'range);
    mj : DoblDobl_Dense_Matrix_Series.Matrix;
    xp,xd : DoblDobl_Dense_Vector_Series.Vector;
    nrm : double_double;
    tol : constant double_float := 1.0E-26;
    wjp : DoblDobl_Series_Jaco_Matrices.Jaco_Mat(jp'range(1),jp'range(2));
    one : double_double := create(1.0);

  begin
    DoblDobl_Series_Poly_Systems.Copy(p,wp);
    Series_and_Polynomials.Set_degree(wp,degree);
    px := DoblDobl_Series_Poly_SysFun.Eval(wp,x);
    nrm := DoblDobl_Series_Vector_Norms.Max_Norm(px);
    if nrm > tol then
      for i in jp'range(1) loop
        for j in jp'range(2) loop
          DoblDobl_Series_Polynomials.Copy(jp(i,j),wjp(i,j));
        end loop;
      end loop;
      Series_and_Polynomials.Set_degree(wjp,degree);
      DoblDobl_Dense_Series_Vectors.Min(px);
      Series_and_Polynomials.Set_degree(px,degree);
      jm := DoblDobl_Series_Jaco_Matrices.Eval(wjp,x);
      Series_and_Polynomials.Set_degree(jm,degree);
      mj := DoblDobl_Dense_Matrix_Series.Create(jm);
      xp := DoblDobl_Dense_Vector_Series.Create(px);
      Solve_by_SVD(mj,xp,info,rcond,xd);
      if one + rcond /= one then
        dx := DoblDobl_Dense_Vector_Series.Create(xd);
        DoblDobl_Dense_Series_Vectors.Add(x,dx);
      end if;
      DoblDobl_Series_Jaco_Matrices.Clear(wjp);
    end if;
    DoblDobl_Series_Poly_Systems.Clear(wp);
    DoblDobl_Dense_Matrix_Series.Clear(mj);
    DoblDobl_Dense_Vector_Series.Clear(xp);
    DoblDobl_Dense_Vector_Series.Clear(xd);
  end SVD_Newton_Step;

  procedure SVD_Newton_Step
              ( p : in DoblDobl_Series_Poly_Systems.Poly_Sys;
                degree : in integer32;
                x : in out DoblDobl_Dense_Series_Vectors.Vector;
                info : out integer32; rcond : out double_double ) is

    jp : DoblDobl_Series_Jaco_Matrices.Jaco_Mat(p'range,x'range)
       := DoblDobl_Series_Jaco_Matrices.Create(p);

  begin
    SVD_Newton_Step(p,jp,degree,x,info,rcond);
    DoblDobl_Series_Jaco_Matrices.Clear(jp);
  end SVD_Newton_Step;

  procedure SVD_Newton_Step
              ( file : in file_type;
                p : in DoblDobl_Series_Poly_Systems.Poly_Sys;
                jp : in DoblDobl_Series_Jaco_Matrices.Jaco_Mat;
                degree : in integer32;
                x : in out DoblDobl_Dense_Series_Vectors.Vector;
                info : out integer32; rcond : out double_double ) is

    wp : DoblDobl_Series_Poly_Systems.Poly_Sys(p'range);
    dx : DoblDobl_Dense_Series_Vectors.Vector(x'range);
    px : DoblDobl_Dense_Series_Vectors.Vector(p'range);
    jm : DoblDobl_Dense_Series_Matrices.Matrix(p'range,x'range);
    mj : DoblDobl_Dense_Matrix_Series.Matrix;
    xp,xd : DoblDobl_Dense_Vector_Series.Vector;
    nrm : double_double;
    tol : constant double_float := 1.0E-26;
    wjp : DoblDobl_Series_Jaco_Matrices.Jaco_Mat(jm'range(1),jm'range(2));
    one : constant double_double := create(1.0);

  begin
    DoblDobl_Series_Poly_Systems.Copy(p,wp);
    Series_and_Polynomials.Set_degree(wp,degree);
    px := DoblDobl_Series_Poly_SysFun.Eval(wp,x);
    put_line(file,"The evaluated series :");
    for i in px'range loop
      DoblDobl_Dense_Series_io.put(file,px(i)); new_line(file);
    end loop;
    nrm := DoblDobl_Series_Vector_Norms.Max_Norm(px);
    put(file,"The max norm of the evaluation : ");
    put(file,nrm,3); new_line(file);
    if nrm > tol then
      for i in jp'range(1) loop
        for j in jp'range(2) loop
          DoblDobl_Series_Polynomials.Copy(jp(i,j),wjp(i,j));
        end loop;
      end loop;
      Series_and_Polynomials.Set_degree(wjp,degree);
      DoblDobl_Dense_Series_Vectors.Min(px);
      Series_and_Polynomials.Set_degree(px,degree);
      jm := DoblDobl_Series_Jaco_Matrices.Eval(jp,x);
      Series_and_Polynomials.Set_degree(jm,degree);
      mj := DoblDobl_Dense_Matrix_Series.Create(jm);
      xp := DoblDobl_Dense_Vector_Series.Create(px);
      Solve_by_SVD(mj,xp,info,rcond,xd);
      put(file,"SVD rcond : "); put(file,rcond,3); new_line(file);
      if one + rcond /= one then
        dx := DoblDobl_Dense_Vector_Series.Create(xd);
        put_line(file,"The update to the series :");
        for i in dx'range loop
          DoblDobl_Dense_Series_io.put(file,dx(i)); new_line(file);
        end loop;
        DoblDobl_Dense_Series_Vectors.Add(x,dx);
      end if;
      DoblDobl_Series_Jaco_Matrices.Clear(wjp);
    end if;
    DoblDobl_Series_Poly_Systems.Clear(wp);
    DoblDobl_Dense_Matrix_Series.Clear(mj);
    DoblDobl_Dense_Vector_Series.Clear(xp);
    DoblDobl_Dense_Vector_Series.Clear(xd);
  end SVD_Newton_Step;

  procedure SVD_Newton_Step
              ( file : in file_type;
                p : in DoblDobl_Series_Poly_Systems.Poly_Sys;
                degree : in integer32;
                x : in out DoblDobl_Dense_Series_Vectors.Vector;
                info : out integer32; rcond : out double_double ) is

    jp : DoblDobl_Series_Jaco_Matrices.Jaco_Mat(p'range,x'range)
       := DoblDobl_Series_Jaco_Matrices.Create(p);

  begin
    SVD_Newton_Step(file,p,jp,degree,x,info,rcond);
    DoblDobl_Series_Jaco_Matrices.Clear(jp);
  end SVD_Newton_Step;

-- ONE NEWTON STEP WITH ECHELON FORM :

  procedure Echelon_Newton_Step
              ( p : in DoblDobl_Series_Poly_Systems.Poly_Sys;
                jp : in DoblDobl_Series_Jaco_Matrices.Jaco_Mat;
                degree : in integer32;
                x : in out DoblDobl_Dense_Series_Vectors.Vector;
                det : out Complex_Number ) is

    dx : DoblDobl_Dense_Series_Vectors.Vector(x'range);
    px : DoblDobl_Dense_Series_Vectors.Vector(p'range);
    jm : DoblDobl_Dense_Series_Matrices.Matrix(p'range,x'range);
    xp,xd,xdn : DoblDobl_Dense_Vector_Series.Vector;
    mj : DoblDobl_Dense_Matrix_Series.Matrix;

  begin
    px := DoblDobl_Series_Poly_SysFun.Eval(p,x);
    DoblDobl_Dense_Series_Vectors.Min(px);
    Series_and_Polynomials.Set_degree(px,degree);
    jm := DoblDobl_Series_Jaco_Matrices.Eval(jp,x);
    Series_and_Polynomials.Set_degree(jm,degree);
    mj := DoblDobl_Dense_Matrix_Series.Create(jm);
    xp := DoblDobl_Dense_Vector_Series.Create(px);
    Echelon_Solve(mj,xp,det,xd,xdn);
    dx := DoblDobl_Dense_Vector_Series.Create(xd);
    DoblDobl_Dense_Series_Vectors.Add(x,dx);
    DoblDobl_Dense_Matrix_Series.Clear(mj);
    DoblDobl_Dense_Vector_Series.Clear(xp);
    DoblDobl_Dense_Vector_Series.Clear(xd);
  end Echelon_Newton_Step;

  procedure Echelon_Newton_Step
              ( p : in DoblDobl_Series_Poly_Systems.Poly_Sys;
                degree : in integer32;
                x : in out DoblDobl_Dense_Series_Vectors.Vector;
                det : out Complex_Number ) is

    jp : DoblDobl_Series_Jaco_Matrices.Jaco_Mat(p'range,x'range)
       := DoblDobl_Series_Jaco_Matrices.Create(p);

  begin
    Echelon_Newton_Step(p,jp,degree,x,det);
    DoblDobl_Series_Jaco_Matrices.Clear(jp);
  end Echelon_Newton_Step;

  procedure Echelon_Newton_Step
              ( file : in file_type;
                p : in DoblDobl_Series_Poly_Systems.Poly_Sys;
                jp : in DoblDobl_Series_Jaco_Matrices.Jaco_Mat;
                degree : in integer32;
                x : in out DoblDobl_Dense_Series_Vectors.Vector;
                det : out Complex_Number ) is

    dx : DoblDobl_Dense_Series_Vectors.Vector(x'range);
    px : DoblDobl_Dense_Series_Vectors.Vector(p'range);
    jm : DoblDobl_Dense_Series_Matrices.Matrix(p'range,x'range);
    xp,xd,xdn : DoblDobl_Dense_Vector_Series.Vector;
    mj : DoblDobl_Dense_Matrix_Series.Matrix;

  begin
    px := DoblDobl_Series_Poly_SysFun.Eval(p,x);
    put_line(file,"The evaluated series :");
    for i in px'range loop
      DoblDobl_Dense_Series_io.put(file,px(i)); new_line(file);
    end loop;
    DoblDobl_Dense_Series_Vectors.Min(px);
    Series_and_Polynomials.Set_degree(px,degree);
    jm := DoblDobl_Series_Jaco_Matrices.Eval(jp,x);
    Series_and_Polynomials.Set_degree(jm,degree);
    mj := DoblDobl_Dense_Matrix_Series.Create(jm);
    xp := DoblDobl_Dense_Vector_Series.Create(px);
    Echelon_Solve(mj,xp,det,xd,xdn);
    put(file,"n.deg : "); put(file,xdn.deg,1); 
    put(file,"  det : "); put(file,det); new_line(file);
    dx := DoblDobl_Dense_Vector_Series.Create(xd);
    DoblDobl_Dense_Series_Vectors.Add(x,dx);
    put_line(file,"The update to the series :");
    for i in dx'range loop
      DoblDobl_Dense_Series_io.put(file,dx(i)); new_line(file);
    end loop;
    DoblDobl_Dense_Matrix_Series.Clear(mj);
    DoblDobl_Dense_Vector_Series.Clear(xp);
    DoblDobl_Dense_Vector_Series.Clear(xd);
  end Echelon_Newton_Step;

  procedure Echelon_Newton_Step
              ( file : in file_type;
                p : in DoblDobl_Series_Poly_Systems.Poly_Sys;
                degree : in integer32;
                x : in out DoblDobl_Dense_Series_Vectors.Vector;
                det : out Complex_Number ) is

    jp : DoblDobl_Series_Jaco_Matrices.Jaco_Mat(p'range,x'range)
       := DoblDobl_Series_Jaco_Matrices.Create(p);

  begin
    Echelon_Newton_Step(file,p,jp,degree,x,det);
    DoblDobl_Series_Jaco_Matrices.Clear(jp);
  end Echelon_Newton_Step;

-- MANY NEWTON STEPS WITH LU WITHOUT CONDITION NUMBER ESTIMATE :

  procedure LU_Newton_Steps
              ( p : in DoblDobl_Series_Poly_Systems.Poly_Sys;
                jp : in DoblDobl_Series_Jaco_Matrices.Jaco_Mat;
                degree : in out integer32; nbrit : in integer32;
                x : in out DoblDobl_Dense_Series_Vectors.Vector;
                info : out integer32 ) is

  begin
    for i in 1..nbrit loop
      LU_Newton_Step(p,jp,degree,x,info);
      exit when (info /= 0); -- stop if Jacobian matrix is singular
      exit when (i = nbrit); -- do not double degree after last step
      degree := 2*degree;
      if degree > DoblDobl_Dense_Series.max_deg
       then degree := DoblDobl_Dense_Series.max_deg;
      end if;
    end loop;
  end LU_Newton_Steps;

  procedure LU_Newton_Steps
              ( p : in DoblDobl_Series_Poly_Systems.Poly_Sys;
                degree : in out integer32; nbrit : in integer32;
                x : in out DoblDobl_Dense_Series_Vectors.Vector;
                info : out integer32 ) is

    jp : DoblDobl_Series_Jaco_Matrices.Jaco_Mat(p'range,x'range)
       := DoblDobl_Series_Jaco_Matrices.Create(p);

  begin
    LU_Newton_Steps(p,jp,degree,nbrit,x,info);
    DoblDobl_Series_Jaco_Matrices.Clear(jp);
  end LU_Newton_Steps;

  procedure LU_Newton_Steps
              ( file : in file_type;
                p : in DoblDobl_Series_Poly_Systems.Poly_Sys;
                jp : in DoblDobl_Series_Jaco_Matrices.Jaco_Mat;
                degree : in out integer32; nbrit : in integer32;
                x : in out DoblDobl_Dense_Series_Vectors.Vector;
                info : out integer32 ) is

  begin
    for i in 1..nbrit loop
      put(file,"LU Newton step "); put(file,i,1); put_line(file," :");
      LU_Newton_Step(file,p,jp,degree,x,info);
      exit when (info /= 0); -- stop if Jacobian matrix is singular
      exit when (i = nbrit); -- do not double degree after last step
      degree := 2*degree;
      if degree > DoblDobl_Dense_Series.max_deg
       then degree := DoblDobl_Dense_Series.max_deg;
      end if;
    end loop;
  end LU_Newton_Steps;

  procedure LU_Newton_Steps
              ( file : in file_type;
                p : in DoblDobl_Series_Poly_Systems.Poly_Sys;
                degree : in out integer32; nbrit : in integer32;
                x : in out DoblDobl_Dense_Series_Vectors.Vector;
                info : out integer32 ) is

    jp : DoblDobl_Series_Jaco_Matrices.Jaco_Mat(p'range,x'range)
       := DoblDobl_Series_Jaco_Matrices.Create(p);

  begin
    LU_Newton_Steps(file,p,jp,degree,nbrit,x,info);
    DoblDobl_Series_Jaco_Matrices.Clear(jp);
  end LU_Newton_Steps;

-- MANY NEWTON STEPS WITH LU WITH CONDITION NUMBER ESTIMATE :

  procedure LU_Newton_Steps
              ( p : in DoblDobl_Series_Poly_Systems.Poly_Sys;
                jp : in DoblDobl_Series_Jaco_Matrices.Jaco_Mat;
                degree : in out integer32; nbrit : in integer32;
                x : in out DoblDobl_Dense_Series_Vectors.Vector;
                rcond : out double_double ) is

    one : constant double_double := create(1.0);

  begin
    for i in 1..nbrit loop
      LU_Newton_Step(p,jp,degree,x,rcond);
      exit when (one + rcond = one); -- stop if Jacobian matrix is singular
      exit when (i = nbrit); -- do not double degree after last step
      degree := 2*degree;
      if degree > DoblDobl_Dense_Series.max_deg
       then degree := DoblDobl_Dense_Series.max_deg;
      end if;
    end loop;
  end LU_Newton_Steps;

  procedure LU_Newton_Steps
              ( p : in DoblDobl_Series_Poly_Systems.Poly_Sys;
                degree : in out integer32; nbrit : in integer32;
                x : in out DoblDobl_Dense_Series_Vectors.Vector;
                rcond : out double_double ) is

    jp : DoblDobl_Series_Jaco_Matrices.Jaco_Mat(p'range,x'range)
       := DoblDobl_Series_Jaco_Matrices.Create(p);

  begin
    LU_Newton_Steps(p,jp,degree,nbrit,x,rcond);
    DoblDobl_Series_Jaco_Matrices.Clear(jp);
  end LU_Newton_Steps;

  procedure LU_Newton_Steps
              ( file : in file_type;
                p : in DoblDobl_Series_Poly_Systems.Poly_Sys;
                jp : in DoblDobl_Series_Jaco_Matrices.Jaco_Mat;
                degree : in out integer32; nbrit : in integer32;
                x : in out DoblDobl_Dense_Series_Vectors.Vector;
                rcond : out double_double ) is

    one : constant double_double := create(1.0);

  begin
    for i in 1..nbrit loop
      put(file,"LU Newton step "); put(file,i,1); put_line(file," :");
      LU_Newton_Step(file,p,jp,degree,x,rcond);
      exit when (one + rcond = one); -- stop if Jacobian matrix is singular
      exit when (i = nbrit); -- do not double degree after last step
      degree := 2*degree;
      if degree > DoblDobl_Dense_Series.max_deg
       then degree := DoblDobl_Dense_Series.max_deg;
      end if;
    end loop;
  end LU_Newton_Steps;

  procedure LU_Newton_Steps
              ( file : in file_type;
                p : in DoblDobl_Series_Poly_Systems.Poly_Sys;
                degree : in out integer32; nbrit : in integer32;
                x : in out DoblDobl_Dense_Series_Vectors.Vector;
                rcond : out double_double ) is

    jp : DoblDobl_Series_Jaco_Matrices.Jaco_Mat(p'range,x'range)
       := DoblDobl_Series_Jaco_Matrices.Create(p);

  begin
    LU_Newton_Steps(file,p,jp,degree,nbrit,x,rcond);
    DoblDobl_Series_Jaco_Matrices.Clear(jp);
  end LU_Newton_Steps;

-- MANY NEWTON STEPS WITH QR :

  procedure QR_Newton_Steps
              ( p : in DoblDobl_Series_Poly_Systems.Poly_Sys;
                jp : in DoblDobl_Series_Jaco_Matrices.Jaco_Mat;
                degree : in out integer32; nbrit : in integer32;
                x : in out DoblDobl_Dense_Series_Vectors.Vector;
                info : out integer32 ) is

  begin
    for i in 1..nbrit loop
      QR_Newton_Step(p,jp,degree,x,info);
      exit when (info /= 0); -- stop if Jacobian matrix is singular
      exit when (i = nbrit); -- do not double degree after last step
      degree := 2*degree;
      if degree > DoblDobl_Dense_Series.max_deg
       then degree := DoblDobl_Dense_Series.max_deg;
      end if;
    end loop;
  end QR_Newton_Steps;

  procedure QR_Newton_Steps
              ( p : in DoblDobl_Series_Poly_Systems.Poly_Sys;
                degree : in out integer32; nbrit : in integer32;
                x : in out DoblDobl_Dense_Series_Vectors.Vector;
                info : out integer32 ) is

    jp : DoblDobl_Series_Jaco_Matrices.Jaco_Mat(p'range,x'range)
       := DoblDobl_Series_Jaco_Matrices.Create(p);

  begin
    QR_Newton_Steps(p,jp,degree,nbrit,x,info);
    DoblDobl_Series_Jaco_Matrices.Clear(jp);
  end QR_Newton_Steps;

  procedure QR_Newton_Steps
              ( file : in file_type;
                p : in DoblDobl_Series_Poly_Systems.Poly_Sys;
                jp : in DoblDobl_Series_Jaco_Matrices.Jaco_Mat;
                degree : in out integer32; nbrit : in integer32;
                x : in out DoblDobl_Dense_Series_Vectors.Vector;
                info : out integer32 ) is

  begin
    for i in 1..nbrit loop
      put(file,"QR Newton step "); put(file,i,1); put_line(file," :");
      QR_Newton_Step(file,p,jp,degree,x,info);
      exit when (info /= 0); -- stop if Jacobian matrix is singular
      exit when (i = nbrit); -- do not double degree after last step
      degree := 2*degree;
      if degree > DoblDobl_Dense_Series.max_deg
       then degree := DoblDobl_Dense_Series.max_deg;
      end if;
    end loop;
  end QR_Newton_Steps;

  procedure QR_Newton_Steps
              ( file : in file_type;
                p : in DoblDobl_Series_Poly_Systems.Poly_Sys;
                degree : in out integer32; nbrit : in integer32;
                x : in out DoblDobl_Dense_Series_Vectors.Vector;
                info : out integer32 ) is

    jp : DoblDobl_Series_Jaco_Matrices.Jaco_Mat(p'range,x'range)
       := DoblDobl_Series_Jaco_Matrices.Create(p);

  begin
    QR_Newton_Steps(file,p,jp,degree,nbrit,x,info);
    DoblDobl_Series_Jaco_Matrices.Clear(jp);
  end QR_Newton_Steps;

-- MANY NEWTON STEPS WITH SVD :

  procedure SVD_Newton_Steps
              ( p : in DoblDobl_Series_Poly_Systems.Poly_Sys;
                jp : in DoblDobl_Series_Jaco_Matrices.Jaco_Mat;
                degree : in out integer32; nbrit : in integer32;
                x : in out DoblDobl_Dense_Series_Vectors.Vector;
                info : out integer32; rcond : out double_double ) is

    one : constant double_double := create(1.0);

  begin
    for i in 1..nbrit loop
      SVD_Newton_Step(p,jp,degree,x,info,rcond);
      exit when (one + rcond = one); -- stop if Jacobian matrix is singular
      exit when (i = nbrit); -- do not double degree after last step
      degree := 2*degree;
      if degree > DoblDobl_Dense_Series.max_deg
       then degree := DoblDobl_Dense_Series.max_deg;
      end if;
    end loop;
  end SVD_Newton_Steps;

  procedure SVD_Newton_Steps
              ( p : in DoblDobl_Series_Poly_Systems.Poly_Sys;
                degree : in out integer32; nbrit : in integer32;
                x : in out DoblDobl_Dense_Series_Vectors.Vector;
                info : out integer32; rcond : out double_double ) is

    jp : DoblDobl_Series_Jaco_Matrices.Jaco_Mat(p'range,x'range)
       := DoblDobl_Series_Jaco_Matrices.Create(p);

  begin
    SVD_Newton_Steps(p,jp,degree,nbrit,x,info,rcond);
    DoblDobl_Series_Jaco_Matrices.Clear(jp);
  end SVD_Newton_Steps;

  procedure SVD_Newton_Steps
              ( file : in file_type;
                p : in DoblDobl_Series_Poly_Systems.Poly_Sys;
                jp : in DoblDobl_Series_Jaco_Matrices.Jaco_Mat;
                degree : in out integer32; nbrit : in integer32;
                x : in out DoblDobl_Dense_Series_Vectors.Vector;
                info : out integer32; rcond : out double_double ) is

    one : constant double_double := create(1.0);

  begin
    for i in 1..nbrit loop
      put(file,"SVD Newton step "); put(file,i,1); put_line(file," :");
      SVD_Newton_Step(file,p,jp,degree,x,info,rcond);
      exit when (one + rcond = one); -- stop if Jacobian matrix is singular
      exit when (i = nbrit); -- do not double degree after last step
      degree := 2*degree;
      if degree > DoblDobl_Dense_Series.max_deg
       then degree := DoblDobl_Dense_Series.max_deg;
      end if;
    end loop;
  end SVD_Newton_Steps;

  procedure SVD_Newton_Steps
              ( file : in file_type;
                p : in DoblDobl_Series_Poly_Systems.Poly_Sys;
                degree : in out integer32; nbrit : in integer32;
                x : in out DoblDobl_Dense_Series_Vectors.Vector;
                info : out integer32; rcond : out double_double ) is

    jp : DoblDobl_Series_Jaco_Matrices.Jaco_Mat(p'range,x'range)
       := DoblDobl_Series_Jaco_Matrices.Create(p);

  begin
    SVD_Newton_Steps(file,p,jp,degree,nbrit,x,info,rcond);
    DoblDobl_Series_Jaco_Matrices.Clear(jp);
  end SVD_Newton_Steps;

-- MANY NEWTON STEPS WITH ECHELON FORM :

  procedure Echelon_Newton_Steps
              ( p : in DoblDobl_Series_Poly_Systems.Poly_Sys;
                jp : in DoblDobl_Series_Jaco_Matrices.Jaco_Mat;
                degree : in out integer32; nbrit : in integer32;
                x : in out DoblDobl_Dense_Series_Vectors.Vector;
                det : out Complex_Number ) is
  begin
    for i in 1..nbrit loop
      Echelon_Newton_Step(p,jp,degree,x,det);
      exit when (i = nbrit); -- do not double degree after last step
      degree := 2*degree;
      if degree > DoblDobl_Dense_Series.max_deg
       then degree := DoblDobl_Dense_Series.max_deg;
      end if;
    end loop;
  end Echelon_Newton_Steps;

  procedure Echelon_Newton_Steps
              ( p : in DoblDobl_Series_Poly_Systems.Poly_Sys;
                degree : in out integer32; nbrit : in integer32;
                x : in out DoblDobl_Dense_Series_Vectors.Vector;
                det : out Complex_Number ) is

    jp : DoblDobl_Series_Jaco_Matrices.Jaco_Mat(p'range,x'range)
       := DoblDobl_Series_Jaco_Matrices.Create(p);

  begin
    Echelon_Newton_Steps(p,jp,degree,nbrit,x,det);
    DoblDobl_Series_Jaco_Matrices.Clear(jp);
  end Echelon_Newton_Steps;

  procedure Echelon_Newton_Steps
              ( file : in file_type;
                p : in DoblDobl_Series_Poly_Systems.Poly_Sys;
                jp : in DoblDobl_Series_Jaco_Matrices.Jaco_Mat;
                degree : in out integer32; nbrit : in integer32;
                x : in out DoblDobl_Dense_Series_Vectors.Vector;
                det : out Complex_Number ) is
  begin
    for i in 1..nbrit loop
      put(file,"LU Newton step "); put(file,i,1); put_line(file," :");
      Echelon_Newton_Step(file,p,jp,degree,x,det);
      exit when (i = nbrit); -- do not double degree after last step
      degree := 2*degree;
      if degree > DoblDobl_Dense_Series.max_deg
       then degree := DoblDobl_Dense_Series.max_deg;
      end if;
    end loop;
  end Echelon_Newton_Steps;

  procedure Echelon_Newton_Steps
              ( file : in file_type;
                p : in DoblDobl_Series_Poly_Systems.Poly_Sys;
                degree : in out integer32; nbrit : in integer32;
                x : in out DoblDobl_Dense_Series_Vectors.Vector;
                det : out Complex_Number ) is

    jp : DoblDobl_Series_Jaco_Matrices.Jaco_Mat(p'range,x'range)
       := DoblDobl_Series_Jaco_Matrices.Create(p);

  begin
    Echelon_Newton_Steps(file,p,jp,degree,nbrit,x,det);
    DoblDobl_Series_Jaco_Matrices.Clear(jp);
  end Echelon_Newton_Steps;

end DoblDobl_Newton_Matrix_Series;
