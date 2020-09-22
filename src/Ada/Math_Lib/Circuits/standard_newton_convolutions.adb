with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;
with Standard_Complex_Numbers;
with Standard_Complex_VecVecs_io;        use Standard_Complex_VecVecs_io;
with Standard_Series_Matrix_Solvers;

package body Standard_Newton_Convolutions is

  function Series_Coefficients
             ( v : Standard_Complex_Vectors.Vector;
               d : integer32 )
             return Standard_Complex_VecVecs.VecVec is

    res : Standard_Complex_VecVecs.VecVec(v'range);

  begin
    for k in res'range loop
      declare
        cff : Standard_Complex_Vectors.Vector(0..d);
      begin
        cff(0) := v(k);
        cff(1..d) := (1..d => Standard_Complex_Numbers.Create(0.0));
        res(k) := new Standard_Complex_Vectors.Vector'(cff);
      end;
    end loop;
    return res;
  end Series_Coefficients;

  procedure Minus ( v : in Standard_Complex_VecVecs.VecVec ) is

    cf : Standard_Complex_Vectors.Link_to_Vector;

  begin
    for i in v'range loop
      cf := v(i);
      for j in cf'range loop
        Standard_Complex_Numbers.Min(cf(j));
      end loop;
    end loop;
  end Minus;

  procedure Minus ( deg : in integer32;
                    v : in Standard_Complex_VecVecs.VecVec ) is

    cf : Standard_Complex_Vectors.Link_to_Vector;

  begin
    for i in v'first..deg loop
      cf := v(i);
      for j in cf'range loop
        Standard_Complex_Numbers.Min(cf(j));
      end loop;
    end loop;
  end Minus;

  procedure Power_Divide
	      ( x : in Standard_Complex_VecVecs.VecVec;
	        f : in double_float ) is

    fac : double_float := f;
    lnk : Standard_Complex_Vectors.Link_to_Vector;

  begin
    for k in 1..x'last loop
      lnk := x(k);
      for i in lnk'range loop
        Standard_Complex_Numbers.Div(lnk(i),fac);
      end loop;
      fac := f*fac;
    end loop;
  end Power_Divide;

  procedure Update ( x,y : in Standard_Complex_VecVecs.VecVec ) is

    xcf,ycf : Standard_Complex_Vectors.Link_to_Vector;

  begin
    for i in x'range loop
      xcf := x(i);
      ycf := y(i);
      for j in xcf'range loop
        Standard_Complex_Numbers.Add(xcf(j),ycf(j));
      end loop;
    end loop;
  end Update;

  procedure Update ( deg : in integer32;
                     x,y : in Standard_Complex_VecVecs.VecVec ) is

    xcf,ycf : Standard_Complex_Vectors.Link_to_Vector;

  begin
    for i in x'range loop
      xcf := x(i);
      ycf := y(i);
      for j in xcf'first..deg loop
        Standard_Complex_Numbers.Add(xcf(j),ycf(j));
      end loop;
    end loop;
  end Update;

  procedure Update ( idx,deg : in integer32;
                     x,y : in Standard_Complex_VecVecs.VecVec ) is

    xcf,ycf : Standard_Complex_Vectors.Link_to_Vector;

  begin
    for i in x'range loop
      xcf := x(i);
      ycf := y(i);
      for j in idx..deg loop
        Standard_Complex_Numbers.Add(xcf(j),ycf(j));
      end loop;
    end loop;
  end Update;

  function Max ( v : Standard_Complex_Vectors.Link_to_Vector )
               return double_float is

    res : double_float := Standard_Complex_Numbers.AbsVal(v(v'first));
    val : double_float;

  begin
    for k in v'first+1..v'last loop
      val := Standard_Complex_Numbers.AbsVal(v(k));
      if val > res
       then res := val;
      end if;
    end loop;
    return res;
  end Max;

  function Max ( deg : integer32;
                 v : Standard_Complex_Vectors.Link_to_Vector )
               return double_float is

    res : double_float := Standard_Complex_Numbers.AbsVal(v(v'first));
    val : double_float;

  begin
    for k in v'first+1..deg loop
      val := Standard_Complex_Numbers.AbsVal(v(k));
      if val > res
       then res := val;
      end if;
    end loop;
    return res;
  end Max;

  function Max ( v : Standard_Complex_VecVecs.VecVec ) return double_float is

    res : double_float := Max(v(v'first));
    val : double_float;

  begin
    for k in v'first+1..v'last loop
      val := Max(v(k));
      if val > res
       then res := val;
      end if;
    end loop;
    return res;
  end Max;

  function Max ( deg : integer32;
                 v : Standard_Complex_VecVecs.VecVec ) return double_float is

    res : double_float := Max(deg,v(v'first));
    val : double_float;

  begin
    for k in v'first+1..v'last loop
      val := Max(deg,v(k));
      if val > res
       then res := val;
      end if;
    end loop;
    return res;
  end Max;

  procedure MaxIdx ( v : in Standard_Complex_VecVecs.VecVec;
                     tol : in double_float;
                     maxval : out double_float; idx : out integer32 ) is

    val : double_float;

  begin
    maxval := Max(v(v'first));
    if maxval > tol then
      idx := v'first-1;
    else
      for k in v'first+1..v'last loop
        val := Max(v(k));
        if val < tol then
          maxval := val;
        else
          idx := k-1; return;
        end if;
      end loop;
    end if;
    idx := v'last;
  end MaxIdx;

  procedure MaxIdx ( deg : in integer32;
                     v : in Standard_Complex_VecVecs.VecVec;
                     tol : in double_float;
                     maxval : out double_float; idx : out integer32 ) is

    val : double_float;

  begin
    maxval := Max(v(v'first));
    if maxval > tol then
      idx := v'first-1;
    else
      for k in v'first+1..deg loop
        val := Max(v(k));
        if val < tol then
          maxval := val;
        else
          idx := k-1; return;
        end if;
      end loop;
    end if;
    idx := v'last;
  end MaxIdx;

-- ONE NEWTON STEP WITH LU WITHOUT CONDITION NUMBER ESTIMATE :

  procedure LU_Newton_Step
              ( s : in Standard_Speelpenning_Convolutions.Link_to_System;
                scf : in Standard_Complex_VecVecs.VecVec;
                absdx : out double_float; info : out integer32;
                ipvt : out Standard_Integer_Vectors.Vector;
                wrk : in Standard_Complex_Vectors.Link_to_Vector;
                scaledx : in boolean := true;
                vrblvl : in integer32 := 0 ) is
  begin
    if vrblvl > 0 then
      put_line("-> in standard_newton_convolutions.LU_Newton_Step 1 ...");
    end if;
    Standard_Speelpenning_Convolutions.Compute(s.pwt,s.mxe,scf);
    Standard_Speelpenning_Convolutions.EvalDiff(s,scf);
    Minus(s.vy);
    Standard_Series_Matrix_Solvers.Solve_by_lufac(s.vm,s.vy,ipvt,info,wrk);
    if scaledx
     then Power_Divide(s.vy,1.0);
    end if;
    Standard_Speelpenning_Convolutions.Delinearize(s.vy,s.yv);
    absdx := max(s.yv);
    Update(scf,s.yv);
  end LU_Newton_Step;

  procedure LU_Newton_Step
              ( file : in file_type;
                s : in Standard_Speelpenning_Convolutions.Link_to_System;
                scf : in Standard_Complex_VecVecs.VecVec;
                absdx : out double_float; info : out integer32;
                ipvt : out Standard_Integer_Vectors.Vector;
                wrk : in Standard_Complex_Vectors.Link_to_Vector;
                scaledx : in boolean := true;
                vrblvl : in integer32 := 0 ) is
  begin
    if vrblvl > 0 then
      put_line("-> in standard_newton_convolutions.LU_Newton_Step 2 ...");
    end if;
    put_line(file,"scf :"); put_line(file,scf);
    Standard_Speelpenning_Convolutions.Compute(s.pwt,s.mxe,scf);
    Standard_Speelpenning_Convolutions.EvalDiff(s,scf);
    put_line(file,"vy :"); put_line(file,s.vy);
    Minus(s.vy);
    Standard_Series_Matrix_Solvers.Solve_by_lufac(s.vm,s.vy,ipvt,info,wrk);
    put_line(file,"dx :"); put_line(file,s.vy);
    if scaledx then
      Power_Divide(s.vy,1.0);
      put_line(file,"scaled dx :"); put_line(file,s.vy);
    end if;
    Standard_Speelpenning_Convolutions.Delinearize(s.vy,s.yv);
    absdx := Max(s.yv);
    put(file,"max |dx| :"); put(file,absdx,3); new_line(file);
    Update(scf,s.yv);
  end LU_Newton_Step;

-- ONE NEWTON STEP WITH LU WITH CONDITION NUMBER ESTIMATE :

  procedure LU_Newton_Step
              ( s : in Standard_Speelpenning_Convolutions.Link_to_System;
                scf : in Standard_Complex_VecVecs.VecVec;
                absdx,rcond : out double_float;
                ipvt : out Standard_Integer_Vectors.Vector;
                wrk : in Standard_Complex_Vectors.Link_to_Vector;
                scaledx : in boolean := true;
                vrblvl : in integer32 := 0 ) is
  begin
    if vrblvl > 0 then
      put_line("-> in standard_newton_convolutions.LU_Newton_Step 3 ...");
    end if;
    Standard_Speelpenning_Convolutions.Compute(s.pwt,s.mxe,scf);
    Standard_Speelpenning_Convolutions.EvalDiff(s,scf);
    Minus(s.vy);
    Standard_Series_Matrix_Solvers.Solve_by_lufco(s.vm,s.vy,ipvt,rcond,wrk);
    if scaledx
     then Power_Divide(s.vy,1.0);
    end if;
    Standard_Speelpenning_Convolutions.Delinearize(s.vy,s.yv);
    absdx := Max(s.yv);
    Update(scf,s.yv);
  end LU_Newton_Step;

  procedure LU_Newton_Step
              ( file : in file_type;
                s : in Standard_Speelpenning_Convolutions.Link_to_System;
                scf : in Standard_Complex_VecVecs.VecVec;
                absdx,rcond : out double_float;
                ipvt : out Standard_Integer_Vectors.Vector;
                wrk : in Standard_Complex_Vectors.Link_to_Vector;
                scaledx : in boolean := true;
                vrblvl : in integer32 := 0 ) is
  begin
    if vrblvl > 0 then
      put_line("-> in standard_newton_convolutions.LU_Newton_Step 4 ...");
    end if;
    put_line(file,"scf :"); put_line(file,scf);
    Standard_Speelpenning_Convolutions.Compute(s.pwt,s.mxe,scf);
    Standard_Speelpenning_Convolutions.EvalDiff(s,scf);
    put_line(file,"vy :"); put_line(file,s.vy);
    Minus(s.vy);
    Standard_Series_Matrix_Solvers.Solve_by_lufco(s.vm,s.vy,ipvt,rcond,wrk);
    put_line(file,"dx :"); put_line(file,s.vy);
    if scaledx then
      Power_Divide(s.vy,1.0);
      put_line(file,"scaled dx :"); put_line(file,s.vy);
    end if;
    Standard_Speelpenning_Convolutions.Delinearize(s.vy,s.yv);
    absdx := Max(s.yv);
    put(file,"max |dx| :"); put(file,absdx,3); new_line(file);
    Update(scf,s.yv);
  end LU_Newton_Step;

-- ONE NEWTON STEP WITH QR :

  procedure QR_Newton_Step
              ( s : in Standard_Speelpenning_Convolutions.Link_to_System;
                scf,dx,xd : in Standard_Complex_VecVecs.VecVec;
                absdx : out double_float;
                qraux : out Standard_Complex_Vectors.Vector;
                w1,w2,w3,w4,w5 : in out Standard_Complex_Vectors.Vector;
                info : out integer32;
                ipvt : out Standard_Integer_Vectors.Vector;
                wrk : in Standard_Complex_Vectors.Link_to_Vector;
                scaledx : in boolean := true;
                vrblvl : in integer32 := 0 ) is
  begin
    if vrblvl > 0 then
      put_line("-> in newton_convolutions.QR_Newton_Step 1 ...");
    end if;
    Standard_Speelpenning_Convolutions.Compute(s.pwt,s.mxe,scf);
    Standard_Speelpenning_Convolutions.EvalDiff(s,scf);
    Minus(s.vy);
    Standard_Series_Matrix_Solvers.Solve_by_QRLS
      (s.vm,s.vy,xd,qraux,w1,w2,w3,w4,w5,ipvt,info,wrk);
    if scaledx
     then Power_Divide(xd,1.0);
    end if;
    Standard_Speelpenning_Convolutions.Delinearize(xd,dx);
    absdx := Max(dx);
    Update(scf,dx);
  end QR_Newton_Step;

  procedure QR_Newton_Step
              ( file : in file_type;
                s : in Standard_Speelpenning_Convolutions.Link_to_System;
                scf,dx,xd : in Standard_Complex_VecVecs.VecVec;
                absdx : out double_float;
                qraux : out Standard_Complex_Vectors.Vector;
                w1,w2,w3,w4,w5 : in out Standard_Complex_Vectors.Vector;
                info : out integer32;
                ipvt : out Standard_Integer_Vectors.Vector;
                wrk : in Standard_Complex_Vectors.Link_to_Vector;
                scaledx : in boolean := true;
                vrblvl : in integer32 := 0 ) is
  begin
    if vrblvl > 0 then
      put_line("-> in newton_convolutions.QR_Newton_Step 2 ...");
    end if;
    put_line(file,"scf :"); put_line(file,scf);
    Standard_Speelpenning_Convolutions.Compute(s.pwt,s.mxe,scf);
    Standard_Speelpenning_Convolutions.EvalDiff(s,scf);
    put_line(file,"vy :"); put_line(file,s.vy);
    Minus(s.vy);
    Standard_Series_Matrix_Solvers.Solve_by_QRLS
      (s.vm,s.vy,xd,qraux,w1,w2,w3,w4,w5,ipvt,info,wrk);
    put_line(file,"dx :"); put_line(file,xd);
    if scaledx then
      Power_Divide(xd,1.0);
      put(file,"scaled dx :"); put_line(file,xd);
    end if;
    Standard_Speelpenning_Convolutions.Delinearize(xd,dx);
    absdx := Max(dx);
    put(file,"max |dx| :"); put(file,absdx,3); new_line(file);
    Update(scf,dx);
  end QR_Newton_Step;

-- ONE NEWTON STEP WITH SVD :

  procedure SVD_Newton_Step
              ( s : in Standard_Speelpenning_Convolutions.Link_to_System;
                scf,dx,xd : in Standard_Complex_VecVecs.VecVec;
                absdx : out double_float;
                svl : out Standard_Complex_Vectors.Vector;
                U,V : out Standard_Complex_Matrices.Matrix;
                info : out integer32; rcond : out double_float;
                ewrk : in Standard_Complex_Vectors.Link_to_Vector;
                wrkv : in Standard_Complex_Vectors.Link_to_Vector;
                scaledx : in boolean := true;
                vrblvl : in integer32 := 0 ) is
  begin
    if vrblvl > 0 then
      put_line("-> in newton_convolutions.SVD_Newton_Step 1 ...");
    end if;
    Standard_Speelpenning_Convolutions.Compute(s.pwt,s.mxe,scf);
    Standard_Speelpenning_Convolutions.EvalDiff(s,scf);
    Minus(s.vy);
    Standard_Series_Matrix_Solvers.Solve_by_SVD
      (s.vm,s.vy,xd,svl,U,V,info,rcond,ewrk,wrkv);
    if scaledx
     then Power_Divide(xd,1.0);
    end if;
    Standard_Speelpenning_Convolutions.Delinearize(xd,dx);
    absdx := Max(dx);
    Update(scf,dx);
  end SVD_Newton_Step;

  procedure SVD_Newton_Step
              ( file : in file_type;
                s : in Standard_Speelpenning_Convolutions.Link_to_System;
                scf,dx,xd : in Standard_Complex_VecVecs.VecVec;
                absdx : out double_float;
                svl : out Standard_Complex_Vectors.Vector;
                U,V : out Standard_Complex_Matrices.Matrix;
                info : out integer32; rcond : out double_float;
                ewrk : in Standard_Complex_Vectors.Link_to_Vector;
                wrkv : in Standard_Complex_Vectors.Link_to_Vector;
                scaledx : in boolean := true;
                vrblvl : in integer32 := 0 ) is
  begin
    if vrblvl > 0 then
      put_line("-> in newton_convolutions.SVD_Newton_Step 2 ...");
    end if;
    put_line(file,"scf :"); put_line(file,scf);
    Standard_Speelpenning_Convolutions.Compute(s.pwt,s.mxe,scf);
    Standard_Speelpenning_Convolutions.EvalDiff(s,scf);
    put_line(file,"vy :"); put_line(file,s.vy);
    Minus(s.vy);
    Standard_Series_Matrix_Solvers.Solve_by_SVD
      (s.vm,s.vy,xd,svl,U,V,info,rcond,ewrk,wrkv);
    put_line(file,"dx :"); put_line(file,xd);
    if scaledx then
      Power_Divide(xd,1.0);
      put(file,"scaled dx :"); put_line(file,xd);
    end if;
    Standard_Speelpenning_Convolutions.Delinearize(xd,dx);
    absdx := Max(dx);
    put(file,"max |dx| : "); put(file,absdx,3); new_line(file);
    Update(scf,dx);
  end SVD_Newton_Step;

end Standard_Newton_Convolutions;
