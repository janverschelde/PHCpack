with Quad_Double_Numbers_io;             use Quad_Double_Numbers_io;
with QuadDobl_Complex_Numbers;
with QuadDobl_Complex_VecVecs_io;        use QUadDobl_Complex_VecVecs_io;
with QuadDobl_Series_Matrix_Solvers;

package body QuadDobl_Newton_Convolutions is

  function Series_Coefficients
             ( v : QuadDobl_Complex_Vectors.Vector;
               d : integer32 )
             return QuadDobl_Complex_VecVecs.VecVec is

    res : QuadDobl_Complex_VecVecs.VecVec(v'range);
    zero : constant quad_double := create(0.0);

  begin
    for k in res'range loop
      declare
        cff : QuadDobl_Complex_Vectors.Vector(0..d);
      begin
        cff(0) := v(k);
        cff(1..d) := (1..d => QuadDobl_Complex_Numbers.Create(zero));
        res(k) := new QuadDobl_Complex_Vectors.Vector'(cff);
      end;
    end loop;
    return res;
  end Series_Coefficients;

  procedure Minus ( v : in QuadDobl_Complex_VecVecs.VecVec ) is

    cf : QuadDobl_Complex_Vectors.Link_to_Vector;

  begin
    for i in v'range loop
      cf := v(i);
      for j in cf'range loop
        QuadDobl_Complex_Numbers.Min(cf(j));
      end loop;
    end loop;
  end Minus;

  procedure Minus ( deg : in integer32;
                    v : in QuadDobl_Complex_VecVecs.VecVec ) is

    cf : QuadDobl_Complex_Vectors.Link_to_Vector;

  begin
    for i in v'first..deg loop
      cf := v(i);
      for j in cf'range loop
        QuadDobl_Complex_Numbers.Min(cf(j));
      end loop;
    end loop;
  end Minus;

  procedure Power_Divide
	      ( x : in QuadDobl_Complex_VecVecs.VecVec;
	        f : in quad_double ) is

    fac : quad_double := f;
    lnk : QuadDobl_Complex_Vectors.Link_to_Vector;

  begin
    for k in 2..x'last loop
      lnk := x(k);
      for i in lnk'range loop
        QuadDobl_Complex_Numbers.Div(lnk(i),fac);
      end loop;
      fac := f*fac;
    end loop;
  end Power_Divide;

  procedure Update ( x,y : in QuadDobl_Complex_VecVecs.VecVec ) is

    xcf,ycf : QuadDobl_Complex_Vectors.Link_to_Vector;

  begin
    for i in x'range loop
      xcf := x(i);
      ycf := y(i);
      for j in xcf'range loop
        QuadDobl_Complex_Numbers.Add(xcf(j),ycf(j));
      end loop;
    end loop;
  end Update;

  procedure Update ( deg : in integer32;
                     x,y : in QuadDobl_Complex_VecVecs.VecVec ) is

    xcf,ycf : QuadDobl_Complex_Vectors.Link_to_Vector;

  begin
    for i in x'range loop
      xcf := x(i);
      ycf := y(i);
      for j in xcf'first..deg loop
        QuadDobl_Complex_Numbers.Add(xcf(j),ycf(j));
      end loop;
    end loop;
  end Update;

  function Max ( v : QuadDobl_Complex_Vectors.Link_to_Vector )
               return quad_double is

    res : quad_double := QuadDobl_Complex_Numbers.AbsVal(v(v'first));
    val : quad_double;

  begin
    for k in v'first+1..v'last loop
      val := QuadDobl_Complex_Numbers.AbsVal(v(k));
      if val > res
       then res := val;
      end if;
    end loop;
    return res;
  end Max;

  function Max ( deg : integer32;
                 v : QuadDobl_Complex_Vectors.Link_to_Vector )
               return quad_double is

    res : quad_double := QuadDobl_Complex_Numbers.AbsVal(v(v'first));
    val : quad_double;

  begin
    for k in v'first+1..deg loop
      val := QuadDobl_Complex_Numbers.AbsVal(v(k));
      if val > res
       then res := val;
      end if;
    end loop;
    return res;
  end Max;

  function Max ( v : QuadDobl_Complex_VecVecs.VecVec ) return quad_double is

    res : quad_double := Max(v(v'first));
    val : quad_double;

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
                 v : QuadDobl_Complex_VecVecs.VecVec ) return quad_double is

    res : quad_double := Max(deg,v(v'first));
    val : quad_double;

  begin
    for k in v'first+1..v'last loop
      val := Max(deg,v(k));
      if val > res
       then res := val;
      end if;
    end loop;
    return res;
  end Max;

  procedure MaxIdx ( v : in QuadDobl_Complex_VecVecs.VecVec;
                     tol : in double_float;
                     maxval : out quad_double; idx : out integer32 ) is

    val : quad_double;

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
                     v : in QuadDobl_Complex_VecVecs.VecVec;
                     tol : in double_float;
                     maxval : out quad_double; idx : out integer32 ) is

    val : quad_double;

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
              ( s : in QuadDobl_Speelpenning_Convolutions.Link_to_System;
                scf : in QuadDobl_Complex_VecVecs.VecVec;
                absdx : out quad_double; info : out integer32;
                ipvt : out Standard_Integer_Vectors.Vector;
                wrk : in QuadDobl_Complex_Vectors.Link_to_Vector;
                scaledx : in boolean := true;
                vrblvl : in integer32 := 0 ) is

    fac : constant quad_double := create(1.0);

  begin
    if vrblvl > 0 then
      put_line("-> in quaddobl_newton_convolutions.LU_Newton_Step 1 ...");
    end if;
    QuadDobl_Speelpenning_Convolutions.Compute(s.pwt,s.mxe,scf);
    QuadDobl_Speelpenning_Convolutions.EvalDiff(s,scf);
    Minus(s.vy);
    QuadDobl_Series_Matrix_Solvers.Solve_by_lufac(s.vm,s.vy,ipvt,info,wrk);
    if scaledx
     then Power_Divide(s.vy,fac);
    end if;
    QuadDobl_Speelpenning_Convolutions.Delinearize(s.vy,s.yv);
    absdx := Max(s.yv);
    Update(scf,s.yv);
  end LU_Newton_Step;

  procedure LU_Newton_Step
              ( file : in file_type;
                s : in QuadDobl_Speelpenning_Convolutions.Link_to_System;
                scf : in QuadDobl_Complex_VecVecs.VecVec;
                absdx : out quad_double; info : out integer32;
                ipvt : out Standard_Integer_Vectors.Vector;
                wrk : in QuadDobl_Complex_Vectors.Link_to_Vector;
                scaledx : in boolean := true;
                vrblvl : in integer32 := 0 ) is

    fac : constant quad_double := create(1.0);

  begin
    if vrblvl > 0 then
      put_line("-> in quaddobl_newton_convolutions.LU_Newton_Step 2 ...");
    end if;
    put_line(file,"scf :"); put_line(file,scf);
    QuadDobl_Speelpenning_Convolutions.Compute(s.pwt,s.mxe,scf);
    QuadDobl_Speelpenning_Convolutions.EvalDiff(s,scf);
    put_line(file,"vy :"); put_line(file,s.vy);
    Minus(s.vy);
    QuadDobl_Series_Matrix_Solvers.Solve_by_lufac(s.vm,s.vy,ipvt,info,wrk);
    put_line(file,"dx :"); put_line(file,s.vy);
    if scaledx then
      Power_Divide(s.vy,fac);
      put_line(file,"scaled dx :"); put_line(file,s.vy);
    end if;
    QuadDobl_Speelpenning_Convolutions.Delinearize(s.vy,s.yv);
    absdx := Max(s.yv);
    put(file,"max |dx| : "); put(file,absdx,3); new_line(file);
    Update(scf,s.yv);
  end LU_Newton_Step;

-- ONE NEWTON STEP WITH LU WITH CONDITION NUMBER ESTIMATE :

  procedure LU_Newton_Step
              ( s : in QuadDobl_Speelpenning_Convolutions.Link_to_System;
                scf : in QuadDobl_Complex_VecVecs.VecVec;
                absdx,rcond : out quad_double;
                ipvt : out Standard_Integer_Vectors.Vector;
                wrk : in QuadDobl_Complex_Vectors.Link_to_Vector;
                scaledx : in boolean := true;
                vrblvl : in integer32 := 0 ) is

    fac : constant quad_double := create(1.0);

  begin
    if vrblvl > 0 then
      put_line("-> in quaddobl_newton_convolutions.LU_Newton_Step 3 ...");
    end if;
    QuadDobl_Speelpenning_Convolutions.Compute(s.pwt,s.mxe,scf);
    QuadDobl_Speelpenning_Convolutions.EvalDiff(s,scf);
    Minus(s.vy);
    QuadDobl_Series_Matrix_Solvers.Solve_by_lufco(s.vm,s.vy,ipvt,rcond,wrk);
    if scaledx
     then Power_Divide(s.vy,fac);
    end if;
    QuadDobl_Speelpenning_Convolutions.Delinearize(s.vy,s.yv);
    absdx := Max(s.yv);
    Update(scf,s.yv);
  end LU_Newton_Step;

  procedure LU_Newton_Step
              ( file : in file_type;
                s : in QuadDobl_Speelpenning_Convolutions.Link_to_System;
                scf : in QuadDobl_Complex_VecVecs.VecVec;
                absdx,rcond : out quad_double;
                ipvt : out Standard_Integer_Vectors.Vector;
                wrk : in QuadDobl_Complex_Vectors.Link_to_Vector;
                scaledx : in boolean := true;
                vrblvl : in integer32 := 0 ) is

    fac : constant quad_double := create(1.0);

  begin
    if vrblvl > 0 then
      put_line("-> in quaddobl_newton_convolutions.LU_Newton_Step 4 ...");
    end if;
    put_line(file,"scf :"); put_line(file,scf);
    QuadDobl_Speelpenning_Convolutions.Compute(s.pwt,s.mxe,scf);
    QuadDobl_Speelpenning_Convolutions.EvalDiff(s,scf);
    put_line(file,"vy :"); put_line(file,s.vy);
    Minus(s.vy);
    QuadDobl_Series_Matrix_Solvers.Solve_by_lufco(s.vm,s.vy,ipvt,rcond,wrk);
    put_line(file,"dx :"); put_line(file,s.vy);
    if scaledx then
      Power_Divide(s.vy,fac);
      put_line(file,"scaled dx :"); put_line(file,s.vy);
    end if;
    QuadDobl_Speelpenning_Convolutions.Delinearize(s.vy,s.yv);
    absdx := Max(s.yv);
    put(file,"max |dx| :"); put(file,absdx,3); new_line(file);
    Update(scf,s.yv);
  end LU_Newton_Step;

-- ONE NEWTON STEP WITH QR :

  procedure QR_Newton_Step
              ( s : in QuadDobl_Speelpenning_Convolutions.Link_to_System;
                scf,dx,xd : in QuadDobl_Complex_VecVecs.VecVec;
                absdx : out quad_double;
                qraux : out QuadDobl_Complex_Vectors.Vector;
                w1,w2,w3,w4,w5 : in out QuadDobl_Complex_Vectors.Vector;
                info : out integer32;
                ipvt : out Standard_Integer_Vectors.Vector;
                wrk : in QuadDobl_Complex_Vectors.Link_to_Vector;
                scaledx : in boolean := true;
                vrblvl : in integer32 := 0 ) is

    fac : constant quad_double := create(1.0);

  begin
    if vrblvl > 0 then
      put_line("-> in quaddobl_newton_convolutions.QR_Newton_Step 1 ...");
    end if;
    QuadDobl_Speelpenning_Convolutions.Compute(s.pwt,s.mxe,scf);
    QuadDobl_Speelpenning_Convolutions.EvalDiff(s,scf);
    Minus(s.vy);
    QuadDobl_Series_Matrix_Solvers.Solve_by_QRLS
      (s.vm,s.vy,xd,qraux,w1,w2,w3,w4,w5,ipvt,info,wrk);
    if scaledx
     then Power_Divide(xd,fac);
    end if;
    QuadDobl_Speelpenning_Convolutions.Delinearize(xd,dx);
    absdx := Max(dx);
    Update(scf,dx);
  end QR_Newton_Step;

  procedure QR_Newton_Step
              ( file : in file_type;
                s : in QuadDobl_Speelpenning_Convolutions.Link_to_System;
                scf,dx,xd : in QuadDobl_Complex_VecVecs.VecVec;
                absdx : out quad_double;
                qraux : out QuadDobl_Complex_Vectors.Vector;
                w1,w2,w3,w4,w5 : in out QuadDobl_Complex_Vectors.Vector;
                info : out integer32;
                ipvt : out Standard_Integer_Vectors.Vector;
                wrk : in QuadDobl_Complex_Vectors.Link_to_Vector;
                scaledx : in boolean := true;
                vrblvl : in integer32 := 0 ) is

    fac : constant quad_double := create(1.0);

  begin
    if vrblvl > 0 then
      put_line("-> in quaddobl_newton_convolutions.QR_Newton_Step 2 ...");
    end if;
    put_line(file,"scf :"); put_line(file,scf);
    QuadDobl_Speelpenning_Convolutions.Compute(s.pwt,s.mxe,scf);
    QuadDobl_Speelpenning_Convolutions.EvalDiff(s,scf);
    put_line(file,"vy :"); put_line(file,s.vy);
    Minus(s.vy);
    QuadDobl_Series_Matrix_Solvers.Solve_by_QRLS
      (s.vm,s.vy,xd,qraux,w1,w2,w3,w4,w5,ipvt,info,wrk);
    put_line(file,"dx :"); put_line(file,xd);
    if scaledx then
      Power_Divide(xd,fac);
      put(file,"scaled dx :"); put_line(file,xd);
    end if;
    QuadDobl_Speelpenning_Convolutions.Delinearize(xd,dx);
    absdx := Max(dx);
    put(file,"max |dx| : "); put(file,absdx,3); new_line(file);
    Update(scf,dx);
  end QR_Newton_Step;

-- ONE NEWTON STEP WITH SVD :

  procedure SVD_Newton_Step
              ( s : in QuadDobl_Speelpenning_Convolutions.Link_to_System;
                scf,dx,xd : in QuadDobl_Complex_VecVecs.VecVec;
                absdx : out quad_double;
                svl : out QuadDobl_Complex_Vectors.Vector;
                U,V : out QuadDobl_Complex_Matrices.Matrix;
                info : out integer32; rcond : out quad_double;
                ewrk : in QuadDobl_Complex_Vectors.Link_to_Vector;
                wrkv : in QuadDobl_Complex_Vectors.Link_to_Vector;
                scaledx : in boolean := true;
                vrblvl : in integer32 := 0 ) is

    fac : constant quad_double := create(1.0);

  begin
    if vrblvl > 0 then
      put_line("-> in quaddobl_newton_convolutions.SVD_Newton_Step 3 ...");
    end if;
    QuadDobl_Speelpenning_Convolutions.Compute(s.pwt,s.mxe,scf);
    QuadDobl_Speelpenning_Convolutions.EvalDiff(s,scf);
    Minus(s.vy);
    QuadDobl_Series_Matrix_Solvers.Solve_by_SVD
      (s.vm,s.vy,xd,svl,U,V,info,rcond,ewrk,wrkv);
    if scaledx
     then Power_Divide(xd,fac);
    end if;
    QuadDobl_Speelpenning_Convolutions.Delinearize(xd,dx);
    absdx := Max(dx);
    Update(scf,dx);
  end SVD_Newton_Step;

  procedure SVD_Newton_Step
              ( file : in file_type;
                s : in QuadDobl_Speelpenning_Convolutions.Link_to_System;
                scf,dx,xd : in QuadDobl_Complex_VecVecs.VecVec;
                absdx : out quad_double;
                svl : out QuadDobl_Complex_Vectors.Vector;
                U,V : out QuadDobl_Complex_Matrices.Matrix;
                info : out integer32; rcond : out quad_double;
                ewrk : in QuadDobl_Complex_Vectors.Link_to_Vector;
                wrkv : in QuadDobl_Complex_Vectors.Link_to_Vector;
                scaledx : in boolean := true;
                vrblvl : in integer32 := 0 ) is

    fac : constant quad_double := create(1.0);

  begin
    if vrblvl > 0 then
      put_line("-> in quaddobl_newton_convolutions.SVD_Newton_Step 4 ...");
    end if;
    put_line(file,"scf :"); put_line(file,scf);
    QuadDobl_Speelpenning_Convolutions.Compute(s.pwt,s.mxe,scf);
    QuadDobl_Speelpenning_Convolutions.EvalDiff(s,scf);
    put_line(file,"vy :"); put_line(file,s.vy);
    Minus(s.vy);
    QuadDobl_Series_Matrix_Solvers.Solve_by_SVD
      (s.vm,s.vy,xd,svl,U,V,info,rcond,ewrk,wrkv);
    put_line(file,"dx :"); put_line(file,xd);
    if scaledx then
      Power_Divide(xd,fac);
      put(file,"scaled dx :"); put_line(file,xd);
    end if;
    QuadDobl_Speelpenning_Convolutions.Delinearize(xd,dx);
    absdx := Max(dx);
    put(file,"max |dx| : "); put(file,absdx,3); new_line(file);
    Update(scf,dx);
  end SVD_Newton_Step;

end QuadDobl_Newton_Convolutions;
