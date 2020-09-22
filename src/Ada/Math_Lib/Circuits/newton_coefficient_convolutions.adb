with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;
with Double_Double_Numbers_io;           use Double_Double_Numbers_io;
with Standard_Complex_Vectors_io;        use Standard_Complex_Vectors_io;
with Standard_Complex_VecVecs_io;        use Standard_Complex_VecVecs_io;
with DoblDobl_Complex_VecVecs_io;        use DoblDobl_Complex_VecVecs_io;
with Standard_Vector_Splitters;
with Standard_Matrix_Splitters;
with DoblDobl_Vector_Splitters;
with Standard_Inlined_Linear_Solvers;
with Standard_Inlined_Linearization;
with Standard_Series_Matrix_Solvers;
with DoblDobl_Series_Matrix_Solvers;
with Standard_Speelpenning_Convolutions;
with DoblDobl_Speelpenning_Convolutions;
with Standard_Newton_Convolutions;
with DoblDobl_Newton_Convolutions;

package body Newton_Coefficient_Convolutions is

  procedure Tolerance_Index
              ( idx,deg : in integer32;
                v : in Standard_Complex_VecVecs.VecVec;
                tol : in double_float;
                idxtoldx : out integer32; absdx : out double_float ) is

    val : double_float;
    done : boolean := false;

  begin
    absdx := Standard_Newton_Convolutions.Max(v(idx));
    if absdx > tol
     then idxtoldx := idx-1; done := true;
     else idxtoldx := deg;
    end if;
    for k in idx+1..deg loop
      val := Standard_Newton_Convolutions.Max(v(k));
      if not done then
        if val > tol
         then idxtoldx := k-1; done := true;
        end if;
      end if;
      if val > absdx
       then absdx := val;
      end if;
    end loop;
  end Tolerance_Index;

-- ONE INLINED NEWTON STEP WITHOUT CONDITION NUMBER ESTIMATE :

  procedure Inlined_LU_Newton_Step
              ( s : in Standard_Coefficient_Convolutions.Link_to_System;
                scf : in Standard_Complex_VecVecs.VecVec;
                rx,ix : in Standard_Floating_VecVecs.Link_to_VecVec;
                absdx : out double_float; info : out integer32;
                ipvt : out Standard_Integer_Vectors.Vector;
                rc,ic : in Standard_Floating_VecVecs.Link_to_VecVec;
                rv,iv : in Standard_Floating_VecVecVecs.Link_to_VecVecVec;
                rb,ib : in Standard_Floating_VecVecs.Link_to_VecVec;
                ry,iy : in Standard_Floating_Vectors.Link_to_Vector;
                scaledx : in boolean := true;
                vrblvl : in integer32 := 0 ) is
  begin
    if vrblvl > 0 then
      put("-> in newton_coefficient_convolutions.");
      put_line("Inlined_LU_Newton_Step 1 ...");
    end if;
    Standard_Vector_Splitters.Complex_Parts(scf,rx,ix);
    Standard_Coefficient_Convolutions.Compute(s.rpwt,s.ipwt,s.mxe,rx,ix);
    Standard_Coefficient_Convolutions.EvalDiff(s,rx.all,ix.all);
    Standard_Newton_Convolutions.Minus(s.vy);
    Standard_Matrix_Splitters.Complex_Parts(s.vm(0).all,rc,ic);
    Standard_Vector_Splitters.Complex_Parts(s.vy,rb,ib);
    Standard_Matrix_Splitters.Split_Rows(s.vm,rv,iv);
    Standard_Inlined_Linearization.Inlined_Solve_by_lufac
      (s.dim,rc,ic,rv,iv,rb,ib,ipvt,info,ry,iy);
    Standard_Vector_Splitters.Complex_Merge(rb,ib,s.vy);
    if scaledx
     then Standard_Newton_Convolutions.Power_Divide(s.vy,1.0);
    end if;
    Standard_Coefficient_Convolutions.Delinearize(s.vy,s.yv);
    absdx := Standard_Newton_Convolutions.Max(s.yv);
    Standard_Newton_Convolutions.Update(scf,s.yv);
  end Inlined_LU_Newton_Step;

  procedure Inlined_LU_Newton_Step
              ( deg : in integer32;
                s : in Standard_Coefficient_Convolutions.Link_to_System;
                scf : in Standard_Complex_VecVecs.VecVec;
                rx,ix : in Standard_Floating_VecVecs.Link_to_VecVec;
                absdx : out double_float; info : out integer32;
                ipvt : out Standard_Integer_Vectors.Vector;
                rc,ic : in Standard_Floating_VecVecs.Link_to_VecVec;
                rv,iv : in Standard_Floating_VecVecVecs.Link_to_VecVecVec;
                rb,ib : in Standard_Floating_VecVecs.Link_to_VecVec;
                ry,iy : in Standard_Floating_Vectors.Link_to_Vector;
                scaledx : in boolean := true;
                vrblvl : in integer32 := 0 ) is
  begin
    if vrblvl > 0 then
      put("-> in newton_coefficient_convolutions.");
      put_line("Inlined_LU_Newton_Step 2 ...");
    end if;
    Standard_Vector_Splitters.Complex_Parts(deg,scf,rx,ix);
    Standard_Coefficient_Convolutions.Compute(deg,s.rpwt,s.ipwt,s.mxe,rx,ix);
    Standard_Coefficient_Convolutions.EvalDiff(deg,s,rx.all,ix.all);
    Standard_Newton_Convolutions.Minus(deg,s.vy);
    Standard_Matrix_Splitters.Complex_Parts(s.vm(0).all,rc,ic);
    for k in 0..deg loop
      Standard_Vector_Splitters.Complex_Parts(s.vy(k),rb(k),ib(k));
    end loop;
    Standard_Matrix_Splitters.Split_Rows(s.vm,rv,iv);
    Standard_Inlined_Linearization.Inlined_Solve_by_lufac
      (deg,s.dim,rc,ic,rv,iv,rb,ib,ipvt,info,ry,iy);
    for k in 0..deg loop
      Standard_Vector_Splitters.Complex_Merge(rb(k),ib(k),s.vy(k));
    end loop;
    if scaledx
     then Standard_Newton_Convolutions.Power_Divide(s.vy,1.0);
    end if;
    Standard_Coefficient_Convolutions.Delinearize(deg,s.vy,s.yv);
    absdx := Standard_Newton_Convolutions.Max(deg,s.yv);
    Standard_Newton_Convolutions.Update(deg,scf,s.yv);
  end Inlined_LU_Newton_Step;

  procedure Inlined_LU_Newton_Step
              ( idx,deg : in integer32;
                s : in Standard_Coefficient_Convolutions.Link_to_System;
                scf : in Standard_Complex_VecVecs.VecVec;
                rx,ix : in Standard_Floating_VecVecs.Link_to_VecVec;
                toldx : in double_float; idxtoldx : out integer32;
                absdx : out double_float; info : out integer32;
                ipvt : in out Standard_Integer_Vectors.Vector;
                rc,ic : in Standard_Floating_VecVecs.Link_to_VecVec;
                rv,iv : in Standard_Floating_VecVecVecs.Link_to_VecVecVec;
                rb,ib : in Standard_Floating_VecVecs.Link_to_VecVec;
                ry,iy : in Standard_Floating_Vectors.Link_to_Vector;
                scaledx : in boolean := true;
                vrblvl : in integer32 := 0 ) is

    rlnk,ilnk : Standard_Floating_Vectors.Link_to_Vector;

  begin
    if vrblvl > 0 then
      put("-> in newton_coefficient_convolutions.");
      put_line("Inlined_LU_Newton_Step 3 ...");
    end if;
    Standard_Vector_Splitters.Complex_Parts(deg,scf,rx,ix);
    Standard_Coefficient_Convolutions.Compute(deg,s.rpwt,s.ipwt,s.mxe,rx,ix);
    Standard_Coefficient_Convolutions.EvalDiff(deg,s,rx.all,ix.all);
    Standard_Newton_Convolutions.Minus(deg,s.vy);
    for k in idx..deg loop
      Standard_Vector_Splitters.Complex_Parts(s.vy(k),rb(k),ib(k));
    end loop;
    for k in rv'first..deg loop
      Standard_Matrix_Splitters.Split_Rows(s.vm(k),rv(k),iv(k));
    end loop;
    if idx = 0 then
      Standard_Matrix_Splitters.Complex_Parts(s.vm(0).all,rc,ic);
      Standard_Inlined_Linearization.Inlined_Solve_by_lufac
        (deg,s.dim,rc,ic,rv,iv,rb,ib,ipvt,info,ry,iy);
    else
      info := 0; -- assuming the previous info was zero
      for k in idx..deg loop
        rlnk := rb(k); ilnk := ib(k);
        for i in idx..(k-1) loop -- all rhs vectors prior to idx are zero
          Standard_Inlined_Linearization.Row_Matrix_Multiply
            (rv(k-i),iv(k-i),rb(i),ib(i),ry,iy);
          for j in rlnk'range loop         -- subtract A(k-1)*x(k) from b(k)
            rlnk(j) := rlnk(j) - ry(j);
            ilnk(j) := ilnk(j) - iy(j);
          end loop;
        end loop;
        Standard_Inlined_Linear_Solvers.lusolve(rc,ic,s.dim,ipvt,rb(k),ib(k));
      end loop;
    end if;
    for k in idx..deg loop
      Standard_Vector_Splitters.Complex_Merge(rb(k),ib(k),s.vy(k));
    end loop;
    if scaledx
     then Standard_Newton_Convolutions.Power_Divide(s.vy,1.0);
    end if;
    Standard_Coefficient_Convolutions.Delinearize(deg,s.vy,s.yv);
    Tolerance_Index(idx,deg,s.vy,toldx,idxtoldx,absdx);
    Standard_Newton_Convolutions.Update(idx,deg,scf,s.yv);
  end Inlined_LU_Newton_Step;

  procedure Inlined_LU_Newton_Step
              ( file : in file_type;
                s : in Standard_Coefficient_Convolutions.Link_to_System;
                scf : in Standard_Complex_VecVecs.VecVec;
                rx,ix : in Standard_Floating_VecVecs.Link_to_VecVec;
                absdx : out double_float; info : out integer32;
                ipvt : out Standard_Integer_Vectors.Vector;
                rc,ic : in Standard_Floating_VecVecs.Link_to_VecVec;
                rv,iv : in Standard_Floating_VecVecVecs.Link_to_VecVecVec;
                rb,ib : in Standard_Floating_VecVecs.Link_to_VecVec;
                ry,iy : in Standard_Floating_Vectors.Link_to_Vector;
                scaledx : in boolean := true;
                vrblvl : in integer32 := 0 ) is
  begin
    if vrblvl > 0 then
      put("-> in newton_coefficient_convolutions.");
      put_line("Inlined_LU_Newton_Step 4 ...");
    end if;
    put_line(file,"scf :"); put_line(file,scf);
    Standard_Vector_Splitters.Complex_Parts(scf,rx,ix);
    Standard_Coefficient_Convolutions.Compute(s.rpwt,s.ipwt,s.mxe,rx,ix);
    Standard_Coefficient_Convolutions.EvalDiff(s,rx.all,ix.all);
    put_line(file,"vy :"); put_line(file,s.vy);
    Standard_Newton_Convolutions.Minus(s.vy);
    Standard_Matrix_Splitters.Complex_Parts(s.vm(0).all,rc,ic);
    Standard_Vector_Splitters.Complex_Parts(s.vy,rb,ib);
    Standard_Matrix_Splitters.Split_Rows(s.vm,rv,iv);
    Standard_Inlined_Linearization.Inlined_Solve_by_lufac
       (s.dim,rc,ic,rv,iv,rb,ib,ipvt,info,ry,iy);
    put(file,"info : "); put(file,info,1); new_line(file);
    Standard_Vector_Splitters.Complex_Merge(rb,ib,s.vy);
    put_line(file,"dx :"); put_line(file,s.vy);
    if scaledx then
      Standard_Newton_Convolutions.Power_Divide(s.vy,1.0);
      put_line(file,"scaled dx :"); put_line(file,s.vy);
    end if;
    Standard_Coefficient_Convolutions.Delinearize(s.vy,s.yv);
    absdx := Standard_Newton_Convolutions.Max(s.yv);
    put(file,"max |dx| :"); put(file,absdx,3); new_line(file);
    Standard_Newton_Convolutions.Update(scf,s.yv);
  end Inlined_LU_Newton_Step;

  procedure Inlined_LU_Newton_Step
              ( file : in file_type; deg : in integer32;
                s : in Standard_Coefficient_Convolutions.Link_to_System;
                scf : in Standard_Complex_VecVecs.VecVec;
                rx,ix : in Standard_Floating_VecVecs.Link_to_VecVec;
                absdx : out double_float; info : out integer32;
                ipvt : out Standard_Integer_Vectors.Vector;
                rc,ic : in Standard_Floating_VecVecs.Link_to_VecVec;
                rv,iv : in Standard_Floating_VecVecVecs.Link_to_VecVecVec;
                rb,ib : in Standard_Floating_VecVecs.Link_to_VecVec;
                ry,iy : in Standard_Floating_Vectors.Link_to_Vector;
                scaledx : in boolean := true;
                vrblvl : in integer32 := 0 ) is
  begin
    if vrblvl > 0 then
      put("-> in newton_coefficient_convolutions.");
      put_line("Inlined_LU_Newton_Step 5 ...");
    end if;
    put_line(file,"scf :"); put_line(file,scf);
    Standard_Vector_Splitters.Complex_Parts(deg,scf,rx,ix);
    Standard_Coefficient_Convolutions.Compute(deg,s.rpwt,s.ipwt,s.mxe,rx,ix);
    Standard_Coefficient_Convolutions.EvalDiff(deg,s,rx.all,ix.all);
    put_line(file,"vy :");
    for k in 0..deg loop
      put_line(file,s.vy(k)); new_line(file);
    end loop;
    Standard_Newton_Convolutions.Minus(deg,s.vy);
    Standard_Matrix_Splitters.Complex_Parts(s.vm(0).all,rc,ic);
    for k in 0..deg loop -- s.vy, rb, and ib are linearized
      Standard_Vector_Splitters.Complex_Parts(s.vy(k),rb(k),ib(k));
    end loop;
    for k in rv'first..deg loop
      Standard_Matrix_Splitters.Split_Rows(s.vm(k),rv(k),iv(k));
    end loop;
    Standard_Inlined_Linearization.Inlined_Solve_by_lufac
      (deg,s.dim,rc,ic,rv,iv,rb,ib,ipvt,info,ry,iy);
    put(file,"info : "); put(file,info,1); new_line(file);
    for k in 0..deg loop
      Standard_Vector_Splitters.Complex_Merge(rb(k),ib(k),s.vy(k));
    end loop;
    put_line(file,"dx :");
    for k in 0..deg loop
      put_line(file,s.vy(k)); new_line(file);
    end loop;
    if scaledx then
      Standard_Newton_Convolutions.Power_Divide(s.vy,1.0);
      put_line(file,"scaled dx :"); put_line(file,s.vy);
    end if;
    Standard_Coefficient_Convolutions.Delinearize(deg,s.vy,s.yv);
    absdx := Standard_Newton_Convolutions.Max(deg,s.yv);
    put(file,"max |dx| :"); put(file,absdx,3); new_line(file);
    Standard_Newton_Convolutions.Update(deg,scf,s.yv);
  end Inlined_LU_Newton_Step;

  procedure Inlined_LU_Newton_Step
              ( file : in file_type; idx,deg : in integer32;
                s : in Standard_Coefficient_Convolutions.Link_to_System;
                scf : in Standard_Complex_VecVecs.VecVec;
                rx,ix : in Standard_Floating_VecVecs.Link_to_VecVec;
                toldx : in double_float; idxtoldx : out integer32;
                absdx : out double_float; info : out integer32;
                ipvt : in out Standard_Integer_Vectors.Vector;
                rc,ic : in Standard_Floating_VecVecs.Link_to_VecVec;
                rv,iv : in Standard_Floating_VecVecVecs.Link_to_VecVecVec;
                rb,ib : in Standard_Floating_VecVecs.Link_to_VecVec;
                ry,iy : in Standard_Floating_Vectors.Link_to_Vector;
                scaledx : in boolean := true;
                vrblvl : in integer32 := 0 ) is

    rlnk,ilnk : Standard_Floating_Vectors.Link_to_Vector;

  begin
    if vrblvl > 0 then
      put("-> in newton_coefficient_convolutions.");
      put_line("Inlined_LU_Newton_Step 6 ...");
    end if;
    put_line(file,"scf :"); put_line(file,scf);
    Standard_Vector_Splitters.Complex_Parts(deg,scf,rx,ix);
    Standard_Coefficient_Convolutions.Compute(deg,s.rpwt,s.ipwt,s.mxe,rx,ix);
    Standard_Coefficient_Convolutions.EvalDiff(deg,s,rx.all,ix.all);
    put_line(file,"vy :");
    for k in 0..deg loop
      put_line(file,s.vy(k)); new_line(file);
    end loop;
    Standard_Newton_Convolutions.Minus(deg,s.vy);
    for k in idx..deg loop -- s.vy, rb, and ib are linearized
      Standard_Vector_Splitters.Complex_Parts(s.vy(k),rb(k),ib(k));
    end loop;
    for k in rv'first..deg loop
      Standard_Matrix_Splitters.Split_Rows(s.vm(k),rv(k),iv(k));
    end loop;
    if idx = 0 then
      Standard_Matrix_Splitters.Complex_Parts(s.vm(0).all,rc,ic);
      Standard_Inlined_Linearization.Inlined_Solve_by_lufac
        (deg,s.dim,rc,ic,rv,iv,rb,ib,ipvt,info,ry,iy);
      put(file,"info : "); put(file,info,1); new_line(file);
    else
      info := 0; -- assuming previous info was zero ...
      for k in idx..deg loop
        rlnk := rb(k); ilnk := ib(k);
        for i in idx..(k-1) loop -- all rhs vectors prior to idx are zero
          Standard_Inlined_Linearization.Row_Matrix_Multiply
            (rv(k-i),iv(k-i),rb(i),ib(i),ry,iy);
          for j in rlnk'range loop         -- subtract A(k-i)*x(i) from b(k)
            rlnk(j) := rlnk(j) - ry(j);
            ilnk(j) := ilnk(j) - iy(j);
          end loop;
        end loop;
        Standard_Inlined_Linear_Solvers.lusolve(rc,ic,s.dim,ipvt,rb(k),ib(k));
      end loop;
    end if;
    for k in idx..deg loop
      Standard_Vector_Splitters.Complex_Merge(rb(k),ib(k),s.vy(k));
    end loop;
    put_line(file,"dx :");
    for k in 0..deg loop
      put_line(file,s.vy(k)); new_line(file);
    end loop;
    if scaledx then
      Standard_Newton_Convolutions.Power_Divide(s.vy,1.0);
      put_line(file,"scaled dx :"); put_line(file,s.vy);
    end if;
    Standard_Coefficient_Convolutions.Delinearize(deg,s.vy,s.yv);
    Tolerance_Index(idx,deg,s.vy,toldx,idxtoldx,absdx);
    put(file,"max |dx| :"); put(file,absdx,3);
    put(file,"  tolerance index : "); put(file,idxtoldx,1); new_line(file);
    Standard_Newton_Convolutions.Update(idx,deg,scf,s.yv);
  end Inlined_LU_Newton_Step;

-- ONE INLINED NEWTON STEP WITH CONDITION NUMBER ESTIMATE :

  procedure Inlined_LU_Newton_Step
              ( s : in Standard_Coefficient_Convolutions.Link_to_System;
                scf : in Standard_Complex_VecVecs.VecVec;
                rx,ix : in Standard_Floating_VecVecs.Link_to_VecVec;
                absdx,rcond : out double_float;
                ipvt : out Standard_Integer_Vectors.Vector;
                rc,ic : in Standard_Floating_VecVecs.Link_to_VecVec;
                rv,iv : in Standard_Floating_VecVecVecs.Link_to_VecVecVec;
                rb,ib : in Standard_Floating_VecVecs.Link_to_VecVec;
                ry,iy : in Standard_Floating_Vectors.Link_to_Vector;
                scaledx : in boolean := true;
                vrblvl : in integer32 := 0 ) is
  begin
    if vrblvl > 0 then
      put("-> in newton_coefficient_convolutions.");
      put_line("Inlined_LU_Newton_Step 7 ...");
    end if;
    Standard_Vector_Splitters.Complex_Parts(scf,rx,ix);
    Standard_Coefficient_Convolutions.Compute(s.rpwt,s.ipwt,s.mxe,rx,ix);
    Standard_Coefficient_Convolutions.EvalDiff(s,rx.all,ix.all);
    Standard_Newton_Convolutions.Minus(s.vy);
    Standard_Matrix_Splitters.Complex_Parts(s.vm(0).all,rc,ic);
    Standard_Vector_Splitters.Complex_Parts(s.vy,rb,ib);
    Standard_Matrix_Splitters.Split_Rows(s.vm,rv,iv);
    Standard_Inlined_Linearization.Inlined_Solve_by_lufco
       (s.dim,rc,ic,rv,iv,rb,ib,ipvt,rcond,ry,iy);
    Standard_Vector_Splitters.Complex_Merge(rb,ib,s.vy);
    if scaledx
     then Standard_Newton_Convolutions.Power_Divide(s.vy,1.0);
    end if;
    Standard_Coefficient_Convolutions.Delinearize(s.vy,s.yv);
    absdx := Standard_Newton_Convolutions.Max(s.yv);
    Standard_Newton_Convolutions.Update(scf,s.yv);
  end Inlined_LU_Newton_Step;

  procedure Inlined_LU_Newton_Step
              ( deg : in integer32;
                s : in Standard_Coefficient_Convolutions.Link_to_System;
                scf : in Standard_Complex_VecVecs.VecVec;
                rx,ix : in Standard_Floating_VecVecs.Link_to_VecVec;
                absdx,rcond : out double_float;
                ipvt : out Standard_Integer_Vectors.Vector;
                rc,ic : in Standard_Floating_VecVecs.Link_to_VecVec;
                rv,iv : in Standard_Floating_VecVecVecs.Link_to_VecVecVec;
                rb,ib : in Standard_Floating_VecVecs.Link_to_VecVec;
                ry,iy : in Standard_Floating_Vectors.Link_to_Vector;
                scaledx : in boolean := true;
                vrblvl : in integer32 := 0 ) is
  begin
    if vrblvl > 0 then
      put("-> in newton_coefficient_convolutions.");
      put_line("Inlined_LU_Newton_Step 8 ...");
    end if;
    Standard_Vector_Splitters.Complex_Parts(deg,scf,rx,ix);
    Standard_Coefficient_Convolutions.Compute(deg,s.rpwt,s.ipwt,s.mxe,rx,ix);
    Standard_Coefficient_Convolutions.EvalDiff(deg,s,rx.all,ix.all);
    Standard_Newton_Convolutions.Minus(deg,s.vy);
    Standard_Matrix_Splitters.Complex_Parts(s.vm(0).all,rc,ic);
    for k in 0..deg loop
      Standard_Vector_Splitters.Complex_Parts(s.vy(k),rb(k),ib(k));
    end loop;
    Standard_Matrix_Splitters.Split_Rows(s.vm,rv,iv);
    Standard_Inlined_Linearization.Inlined_Solve_by_lufco
       (deg,s.dim,rc,ic,rv,iv,rb,ib,ipvt,rcond,ry,iy);
    for k in 0..deg loop
      Standard_Vector_Splitters.Complex_Merge(rb(k),ib(k),s.vy(k));
    end loop;
    if scaledx
     then Standard_Newton_Convolutions.Power_Divide(s.vy,1.0);
    end if;
    Standard_Coefficient_Convolutions.Delinearize(deg,s.vy,s.yv);
    absdx := Standard_Newton_Convolutions.Max(deg,s.yv);
    Standard_Newton_Convolutions.Update(deg,scf,s.yv);
  end Inlined_LU_Newton_Step;

  procedure Inlined_LU_Newton_Step
              ( idx,deg : in integer32;
                s : in Standard_Coefficient_Convolutions.Link_to_System;
                scf : in Standard_Complex_VecVecs.VecVec;
                rx,ix : in Standard_Floating_VecVecs.Link_to_VecVec;
                toldx : in double_float; idxtoldx : out integer32;
                absdx,rcond : out double_float;
                ipvt : in out Standard_Integer_Vectors.Vector;
                rc,ic : in Standard_Floating_VecVecs.Link_to_VecVec;
                rv,iv : in Standard_Floating_VecVecVecs.Link_to_VecVecVec;
                rb,ib : in Standard_Floating_VecVecs.Link_to_VecVec;
                ry,iy : in Standard_Floating_Vectors.Link_to_Vector;
                scaledx : in boolean := true;
                vrblvl : in integer32 := 0 ) is

    rlnk,ilnk : Standard_Floating_Vectors.Link_to_Vector;

  begin
    if vrblvl > 0 then
      put("-> in newton_coefficient_convolutions.");
      put_line("Inlined_LU_Newton_Step 9 ...");
    end if;
    Standard_Vector_Splitters.Complex_Parts(deg,scf,rx,ix);
    Standard_Coefficient_Convolutions.Compute(deg,s.rpwt,s.ipwt,s.mxe,rx,ix);
    Standard_Coefficient_Convolutions.EvalDiff(deg,s,rx.all,ix.all);
    Standard_Newton_Convolutions.Minus(deg,s.vy);
    for k in idx..deg loop
      Standard_Vector_Splitters.Complex_Parts(s.vy(k),rb(k),ib(k));
    end loop;
    for k in rv'first..deg loop
      Standard_Matrix_Splitters.Split_Rows(s.vm(k),rv(k),iv(k));
    end loop;
    if idx = 0 then
      Standard_Matrix_Splitters.Complex_Parts(s.vm(0).all,rc,ic);
      Standard_Inlined_Linearization.Inlined_Solve_by_lufco
        (deg,s.dim,rc,ic,rv,iv,rb,ib,ipvt,rcond,ry,iy);
    else
      rcond := 1.0; -- assuming perfect condition
      for k in idx..deg loop
        rlnk := rb(k); ilnk := ib(k);
        for i in idx..(k-1) loop -- all rhs vectors prior to idx are zero
          Standard_Inlined_Linearization.Row_Matrix_Multiply
            (rv(k-i),iv(k-i),rb(i),ib(i),ry,iy);
          for j in rlnk'range loop         -- subtract A(k-1)*x(k) from b(k)
            rlnk(j) := rlnk(j) - ry(j);
            ilnk(j) := ilnk(j) - iy(j);
          end loop;
        end loop;
        Standard_Inlined_Linear_Solvers.lusolve(rc,ic,s.dim,ipvt,rb(k),ib(k));
      end loop;
    end if;
    for k in idx..deg loop
      Standard_Vector_Splitters.Complex_Merge(rb(k),ib(k),s.vy(k));
    end loop;
    if scaledx
     then Standard_Newton_Convolutions.Power_Divide(s.vy,1.0);
    end if;
    Standard_Coefficient_Convolutions.Delinearize(deg,s.vy,s.yv);
    Tolerance_Index(idx,deg,s.vy,toldx,idxtoldx,absdx);
    Standard_Newton_Convolutions.Update(idx,deg,scf,s.yv);
  end Inlined_LU_Newton_Step;

  procedure Inlined_LU_Newton_Step
              ( file : in file_type;
                s : in Standard_Coefficient_Convolutions.Link_to_System;
                scf : in Standard_Complex_VecVecs.VecVec;
                rx,ix : in Standard_Floating_VecVecs.Link_to_VecVec;
                absdx,rcond : out double_float;
                ipvt : out Standard_Integer_Vectors.Vector;
                rc,ic : in Standard_Floating_VecVecs.Link_to_VecVec;
                rv,iv : in Standard_Floating_VecVecVecs.Link_to_VecVecVec;
                rb,ib : in Standard_Floating_VecVecs.Link_to_VecVec;
                ry,iy : in Standard_Floating_Vectors.Link_to_Vector;
                scaledx : in boolean := true;
                vrblvl : in integer32 := 0 ) is
  begin
    if vrblvl > 0 then
      put("-> in newton_coefficient_convolutions.");
      put_line("Inlined_LU_Newton_Step 10 ...");
    end if;
    put_line(file,"scf :"); put_line(file,scf);
    Standard_Vector_Splitters.Complex_Parts(scf,rx,ix);
    Standard_Coefficient_Convolutions.Compute(s.rpwt,s.ipwt,s.mxe,rx,ix);
    Standard_Coefficient_Convolutions.EvalDiff(s,rx.all,ix.all);
    put_line(file,"vy :"); put_line(file,s.vy);
    Standard_Newton_Convolutions.Minus(s.vy);
    Standard_Matrix_Splitters.Complex_Parts(s.vm(0).all,rc,ic);
    Standard_Vector_Splitters.Complex_Parts(s.vy,rb,ib);
    Standard_Matrix_Splitters.Split_Rows(s.vm,rv,iv);
    Standard_Inlined_Linearization.Inlined_Solve_by_lufco
       (s.dim,rc,ic,rv,iv,rb,ib,ipvt,rcond,ry,iy);
    put(file,"rcond :"); put(file,rcond,3); new_line(file);
    Standard_Vector_Splitters.Complex_Merge(rb,ib,s.vy);
    put_line(file,"dx :"); put_line(file,s.vy);
    if scaledx then
      Standard_Newton_Convolutions.Power_Divide(s.vy,1.0);
      put_line(file,"scaled dx :"); put_line(file,s.vy);
    end if;
    Standard_Coefficient_Convolutions.Delinearize(s.vy,s.yv);
    absdx := Standard_Newton_Convolutions.Max(s.yv);
    put(file,"max |dx| :"); put(file,absdx,3); new_line(file);
    Standard_Newton_Convolutions.Update(scf,s.yv);
  end Inlined_LU_Newton_Step;

  procedure Inlined_LU_Newton_Step
              ( file : in file_type; deg : in integer32;
                s : in Standard_Coefficient_Convolutions.Link_to_System;
                scf : in Standard_Complex_VecVecs.VecVec;
                rx,ix : in Standard_Floating_VecVecs.Link_to_VecVec;
                absdx,rcond : out double_float;
                ipvt : out Standard_Integer_Vectors.Vector;
                rc,ic : in Standard_Floating_VecVecs.Link_to_VecVec;
                rv,iv : in Standard_Floating_VecVecVecs.Link_to_VecVecVec;
                rb,ib : in Standard_Floating_VecVecs.Link_to_VecVec;
                ry,iy : in Standard_Floating_Vectors.Link_to_Vector;
                scaledx : in boolean := true;
                vrblvl : in integer32 := 0 ) is
  begin
    if vrblvl > 0 then
      put("-> in newton_coefficient_convolutions.");
      put_line("Inlined_LU_Newton_Step 11 ...");
    end if;
    put_line(file,"scf :"); put_line(file,scf);
    Standard_Vector_Splitters.Complex_Parts(deg,scf,rx,ix);
    Standard_Coefficient_Convolutions.Compute(deg,s.rpwt,s.ipwt,s.mxe,rx,ix);
    Standard_Coefficient_Convolutions.EvalDiff(deg,s,rx.all,ix.all);
    put_line(file,"vy :");
    for k in 0..deg loop
      put_line(file,s.vy(k)); new_line(file);
    end loop;
    Standard_Newton_Convolutions.Minus(deg,s.vy);
    Standard_Matrix_Splitters.Complex_Parts(s.vm(0).all,rc,ic);
    for k in 0..deg loop -- s.vy, rb, and ib are linearized
      Standard_Vector_Splitters.Complex_Parts(s.vy(k),rb(k),ib(k));
    end loop;
    for k in rv'first..deg loop
      Standard_Matrix_Splitters.Split_Rows(s.vm(k),rv(k),iv(k));
    end loop;
    Standard_Inlined_Linearization.Inlined_Solve_by_lufco
       (deg,s.dim,rc,ic,rv,iv,rb,ib,ipvt,rcond,ry,iy);
    put(file,"rcond :"); put(file,rcond,3); new_line(file);
    for k in 0..deg loop
      Standard_Vector_Splitters.Complex_Merge(rb(k),ib(k),s.vy(k));
    end loop;
    put_line(file,"dx :");
    for k in 0..deg loop
      put_line(file,s.vy(k)); new_line(file);
    end loop;
    if scaledx then
      Standard_Newton_Convolutions.Power_Divide(s.vy,1.0);
      put_line(file,"scaled dx :"); put_line(file,s.vy);
    end if;
    Standard_Coefficient_Convolutions.Delinearize(deg,s.vy,s.yv);
    absdx := Standard_Newton_Convolutions.Max(deg,s.yv);
    put(file,"max |dx| :"); put(file,absdx,3); new_line(file);
    Standard_Newton_Convolutions.Update(deg,scf,s.yv);
  end Inlined_LU_Newton_Step;

  procedure Inlined_LU_Newton_Step
              ( file : in file_type; idx,deg : in integer32;
                s : in Standard_Coefficient_Convolutions.Link_to_System;
                scf : in Standard_Complex_VecVecs.VecVec;
                rx,ix : in Standard_Floating_VecVecs.Link_to_VecVec;
                toldx : in double_float; idxtoldx : out integer32;
                absdx,rcond : out double_float;
                ipvt : in out Standard_Integer_Vectors.Vector;
                rc,ic : in Standard_Floating_VecVecs.Link_to_VecVec;
                rv,iv : in Standard_Floating_VecVecVecs.Link_to_VecVecVec;
                rb,ib : in Standard_Floating_VecVecs.Link_to_VecVec;
                ry,iy : in Standard_Floating_Vectors.Link_to_Vector;
                scaledx : in boolean := true;
                vrblvl : in integer32 := 0 ) is

    rlnk,ilnk : Standard_Floating_Vectors.Link_to_Vector;

  begin
    if vrblvl > 0 then
      put("-> in newton_coefficient_convolutions.");
      put_line("Inlined_LU_Newton_Step 12 ...");
    end if;
    put_line(file,"scf :"); put_line(file,scf);
    Standard_Vector_Splitters.Complex_Parts(deg,scf,rx,ix);
    Standard_Coefficient_Convolutions.Compute(deg,s.rpwt,s.ipwt,s.mxe,rx,ix);
    Standard_Coefficient_Convolutions.EvalDiff(deg,s,rx.all,ix.all);
    put_line(file,"vy :");
    for k in 0..deg loop
      put_line(file,s.vy(k)); new_line(file);
    end loop;
    Standard_Newton_Convolutions.Minus(deg,s.vy);
    for k in idx..deg loop -- s.vy, rb, and ib are linearized
      Standard_Vector_Splitters.Complex_Parts(s.vy(k),rb(k),ib(k));
    end loop;
    for k in rv'first..deg loop
      Standard_Matrix_Splitters.Split_Rows(s.vm(k),rv(k),iv(k));
    end loop;
    if idx = 0 then
      Standard_Matrix_Splitters.Complex_Parts(s.vm(0).all,rc,ic);
      Standard_Inlined_Linearization.Inlined_Solve_by_lufco
        (deg,s.dim,rc,ic,rv,iv,rb,ib,ipvt,rcond,ry,iy);
      put(file,"rcond :"); put(file,rcond,3); new_line(file);
    else
      rcond := 1.0; -- assuming perfect condition ...
      for k in idx..deg loop
        rlnk := rb(k); ilnk := ib(k);
        for i in idx..(k-1) loop -- all rhs vectors prior to idx are zero
          Standard_Inlined_Linearization.Row_Matrix_Multiply
            (rv(k-i),iv(k-i),rb(i),ib(i),ry,iy);
          for j in rlnk'range loop         -- subtract A(k-1)*x(k) from b(k)
            rlnk(j) := rlnk(j) - ry(j);
            ilnk(j) := ilnk(j) - iy(j);
          end loop;
        end loop;
        Standard_Inlined_Linear_Solvers.lusolve(rc,ic,s.dim,ipvt,rb(k),ib(k));
      end loop;
    end if;
    for k in idx..deg loop
      Standard_Vector_Splitters.Complex_Merge(rb(k),ib(k),s.vy(k));
    end loop;
    put_line(file,"dx :");
    for k in 0..deg loop
      put_line(file,s.vy(k)); new_line(file);
    end loop;
    if scaledx then
      Standard_Newton_Convolutions.Power_Divide(s.vy,1.0);
      put_line(file,"scaled dx :"); put_line(file,s.vy);
    end if;
    Standard_Coefficient_Convolutions.Delinearize(deg,s.vy,s.yv);
    absdx := Standard_Newton_Convolutions.Max(deg,s.yv);
    Tolerance_Index(idx,deg,s.vy,toldx,idxtoldx,absdx);
    put(file,"max |dx| :"); put(file,absdx,3);
    put(file,"  tolerance index : "); put(file,idxtoldx,1); new_line(file);
    Standard_Newton_Convolutions.Update(idx,deg,scf,s.yv);
  end Inlined_LU_Newton_Step;

-- ONE NEWTON STEP WITH LU WITHOUT CONDITION NUMBER ESTIMATE :

  procedure LU_Newton_Step
              ( s : in Standard_Coefficient_Convolutions.Link_to_System;
                scf : in Standard_Complex_VecVecs.VecVec;
                rx,ix : in Standard_Floating_VecVecs.Link_to_VecVec;
                absdx : out double_float; info : out integer32;
                ipvt : out Standard_Integer_Vectors.Vector;
                wrk : in Standard_Complex_Vectors.Link_to_Vector;
                scaledx : in boolean := true;
                vrblvl : in integer32 := 0 ) is
  begin
    if vrblvl > 0 then
      put_line("-> in newton_coefficient_convolutions.LU_Newton_Step 1 ...");
    end if;
    Standard_Vector_Splitters.Complex_Parts(scf,rx,ix);
    Standard_Coefficient_Convolutions.Compute(s.rpwt,s.ipwt,s.mxe,rx,ix);
    Standard_Coefficient_Convolutions.EvalDiff(s,rx.all,ix.all);
    Standard_Newton_Convolutions.Minus(s.vy);
    Standard_Series_Matrix_Solvers.Solve_by_lufac(s.vm,s.vy,ipvt,info,wrk);
    if scaledx
     then Standard_Newton_Convolutions.Power_Divide(s.vy,1.0);
    end if;
    Standard_Coefficient_Convolutions.Delinearize(s.vy,s.yv);
    absdx := Standard_Newton_Convolutions.Max(s.yv);
    Standard_Newton_Convolutions.Update(scf,s.yv);
  end LU_Newton_Step;

  procedure LU_Newton_Step
              ( deg : in integer32;
                s : in Standard_Coefficient_Convolutions.Link_to_System;
                scf : in Standard_Complex_VecVecs.VecVec;
                rx,ix : in Standard_Floating_VecVecs.Link_to_VecVec;
                absdx : out double_float; info : out integer32;
                ipvt : out Standard_Integer_Vectors.Vector;
                wrk : in Standard_Complex_Vectors.Link_to_Vector;
                scaledx : in boolean := true;
                vrblvl : in integer32 := 0 ) is
  begin
    if vrblvl > 0 then
      put_line("-> in newton_coefficient_convolutions.LU_Newton_Step 2 ...");
    end if;
    Standard_Vector_Splitters.Complex_Parts(deg,scf,rx,ix);
    Standard_Coefficient_Convolutions.Compute(deg,s.rpwt,s.ipwt,s.mxe,rx,ix);
    Standard_Coefficient_Convolutions.EvalDiff(deg,s,rx.all,ix.all);
    Standard_Newton_Convolutions.Minus(deg,s.vy);
    Standard_Series_Matrix_Solvers.Solve_by_lufac(deg,s.vm,s.vy,ipvt,info,wrk);
    if scaledx
     then Standard_Newton_Convolutions.Power_Divide(s.vy,1.0);
    end if;
    Standard_Coefficient_Convolutions.Delinearize(deg,s.vy,s.yv);
    absdx := Standard_Newton_Convolutions.Max(deg,s.yv);
    Standard_Newton_Convolutions.Update(deg,scf,s.yv);
  end LU_Newton_Step;

  procedure LU_Newton_Step
              ( file : in file_type;
                s : in Standard_Coefficient_Convolutions.Link_to_System;
                scf : in Standard_Complex_VecVecs.VecVec;
                rx,ix : in Standard_Floating_VecVecs.Link_to_VecVec;
                absdx : out double_float; info : out integer32;
                ipvt : out Standard_Integer_Vectors.Vector;
                wrk : in Standard_Complex_Vectors.Link_to_Vector;
                scaledx : in boolean := true;
                vrblvl : in integer32 := 0 ) is
  begin
    if vrblvl > 0 then
      put_line("-> in newton_coefficient_convolutions.LU_Newton_Step 3 ...");
    end if;
    put_line(file,"scf :"); put_line(file,scf);
    Standard_Vector_Splitters.Complex_Parts(scf,rx,ix);
    Standard_Coefficient_Convolutions.Compute(s.rpwt,s.ipwt,s.mxe,rx,ix);
    Standard_Coefficient_Convolutions.EvalDiff(s,rx.all,ix.all);
    put_line(file,"vy :"); put_line(file,s.vy);
    Standard_Newton_Convolutions.Minus(s.vy);
    Standard_Series_Matrix_Solvers.Solve_by_lufac(s.vm,s.vy,ipvt,info,wrk);
    put_line(file,"dx :"); put_line(file,s.vy);
    if scaledx then
      Standard_Newton_Convolutions.Power_Divide(s.vy,1.0);
      put_line(file,"scaled dx :"); put_line(file,s.vy);
    end if;
    Standard_Speelpenning_Convolutions.Delinearize(s.vy,s.yv);
    absdx := Standard_Newton_Convolutions.Max(s.yv);
    put(file,"max |dx| :"); put(file,absdx,3); new_line(file);
    Standard_Newton_Convolutions.Update(scf,s.yv);
  end LU_Newton_Step;

  procedure LU_Newton_Step
              ( file : in file_type; deg : in integer32;
                s : in Standard_Coefficient_Convolutions.Link_to_System;
                scf : in Standard_Complex_VecVecs.VecVec;
                rx,ix : in Standard_Floating_VecVecs.Link_to_VecVec;
                absdx : out double_float; info : out integer32;
                ipvt : out Standard_Integer_Vectors.Vector;
                wrk : in Standard_Complex_Vectors.Link_to_Vector;
                scaledx : in boolean := true;
                vrblvl : in integer32 := 0 ) is
  begin
    if vrblvl > 0 then
      put_line("-> in newton_coefficient_convolutions.LU_Newton_Step 4 ...");
    end if;
    put_line(file,"scf :"); put_line(file,scf);
    Standard_Vector_Splitters.Complex_Parts(deg,scf,rx,ix);
    Standard_Coefficient_Convolutions.Compute(deg,s.rpwt,s.ipwt,s.mxe,rx,ix);
    Standard_Coefficient_Convolutions.EvalDiff(deg,s,rx.all,ix.all);
    put_line(file,"vy :");
    for k in 0..deg loop
      put_line(file,s.vy(k)); new_line(file);
    end loop;
    Standard_Newton_Convolutions.Minus(deg,s.vy);
    Standard_Series_Matrix_Solvers.Solve_by_lufac(deg,s.vm,s.vy,ipvt,info,wrk);
    put_line(file,"dx :");
    for k in 0..deg loop
      put_line(file,s.vy(k)); new_line(file);
    end loop;
    if scaledx then
      Standard_Newton_Convolutions.Power_Divide(s.vy,1.0);
      put_line(file,"scaled dx :"); put_line(file,s.vy);
    end if;
    Standard_Speelpenning_Convolutions.Delinearize(deg,s.vy,s.yv);
    absdx := Standard_Newton_Convolutions.Max(deg,s.yv);
    put(file,"max |dx| :"); put(file,absdx,3); new_line(file);
    Standard_Newton_Convolutions.Update(deg,scf,s.yv);
  end LU_Newton_Step;

  procedure LU_Newton_Step
              ( s : in DoblDobl_Coefficient_Convolutions.Link_to_System;
                scf : in DoblDobl_Complex_VecVecs.VecVec;
                rhx,ihx : in Standard_Floating_VecVecs.Link_to_VecVec;
                rlx,ilx : in Standard_Floating_VecVecs.Link_to_VecVec;
                absdx : out double_double; info : out integer32;
                ipvt : out Standard_Integer_Vectors.Vector;
                wrk : in DoblDobl_Complex_Vectors.Link_to_Vector;
                scaledx : in boolean := true;
                vrblvl : in integer32 := 0 ) is

    fac : constant double_double := create(1.0);

  begin
    if vrblvl > 0 then
      put_line("-> in newton_coefficient_convolutions.LU_Newton_Step 5 ...");
    end if;
    DoblDobl_Vector_Splitters.Complex_Parts(scf,rhx,ihx,rlx,ilx);
    DoblDobl_Coefficient_Convolutions.Compute
      (s.rhpwt,s.ihpwt,s.rlpwt,s.ilpwt,s.mxe,rhx,ihx,rlx,ilx);
    DoblDobl_Coefficient_Convolutions.EvalDiff
      (s,rhx.all,ihx.all,rlx.all,ilx.all);
    DoblDobl_Newton_Convolutions.Minus(s.vy);
    DoblDobl_Series_Matrix_Solvers.Solve_by_lufac(s.vm,s.vy,ipvt,info,wrk);
    if scaledx
     then DoblDobl_Newton_Convolutions.Power_Divide(s.vy,fac);
    end if;
    DoblDobl_Coefficient_Convolutions.Delinearize(s.vy,s.yv);
    absdx := DoblDobl_Newton_Convolutions.Max(s.yv);
    DoblDobl_Newton_Convolutions.Update(scf,s.yv);
  end LU_Newton_Step;

  procedure LU_Newton_Step
              ( file : in file_type;
                s : in DoblDobl_Coefficient_Convolutions.Link_to_System;
                scf : in DoblDobl_Complex_VecVecs.VecVec;
                rhx,ihx : in Standard_Floating_VecVecs.Link_to_VecVec;
                rlx,ilx : in Standard_Floating_VecVecs.Link_to_VecVec;
                absdx : out double_double; info : out integer32;
                ipvt : out Standard_Integer_Vectors.Vector;
                wrk : in DoblDobl_Complex_Vectors.Link_to_Vector;
                scaledx : in boolean := true;
                vrblvl : in integer32 := 0 ) is

    fac : constant double_double := create(1.0);

  begin
    if vrblvl > 0 then
      put_line("-> in newton_coefficient_convolutions.LU_Newton_Step 6 ...");
    end if;
    put_line(file,"scf :"); put_line(file,scf);
    DoblDobl_Vector_Splitters.Complex_Parts(scf,rhx,ihx,rlx,ilx);
    DoblDobl_Coefficient_Convolutions.Compute
      (s.rhpwt,s.ihpwt,s.rlpwt,s.ilpwt,s.mxe,rhx,ihx,rlx,ilx);
    DoblDobl_Coefficient_Convolutions.EvalDiff
      (s,rhx.all,ihx.all,rlx.all,ilx.all);
    put_line(file,"vy :"); put_line(file,s.vy);
    DoblDobl_Newton_Convolutions.Minus(s.vy);
    DoblDobl_Series_Matrix_Solvers.Solve_by_lufac(s.vm,s.vy,ipvt,info,wrk);
    put_line(file,"dx :"); put_line(file,s.vy);
    if scaledx then
      DoblDobl_Newton_Convolutions.Power_Divide(s.vy,fac);
      put_line(file,"scaled dx :"); put_line(file,s.vy);
    end if;
    DoblDobl_Coefficient_Convolutions.Delinearize(s.vy,s.yv);
    absdx := DoblDobl_Newton_Convolutions.Max(s.yv);
    put(file,"max |dx| : "); put(file,absdx,3); new_line(file);
    DoblDobl_Newton_Convolutions.Update(scf,s.yv);
  end LU_Newton_Step;

-- ONE NEWTON STEP WITH LU WITH CONDITION NUMBER ESTIMATE :

  procedure LU_Newton_Step
              ( s : in Standard_Coefficient_Convolutions.Link_to_System;
                scf : in Standard_Complex_VecVecs.VecVec;
                rx,ix : in Standard_Floating_VecVecs.Link_to_VecVec;
                absdx,rcond : out double_float;
                ipvt : out Standard_Integer_Vectors.Vector;
                wrk : in Standard_Complex_Vectors.Link_to_Vector;
                scaledx : in boolean := true;
                vrblvl : in integer32 := 0 ) is
  begin
    if vrblvl > 0 then
      put_line("-> in newton_coefficient_convolutions.LU_Newton_Step 7 ...");
    end if;
    Standard_Vector_Splitters.Complex_Parts(scf,rx,ix);
    Standard_Coefficient_Convolutions.Compute(s.rpwt,s.ipwt,s.mxe,rx,ix);
    Standard_Coefficient_Convolutions.EvalDiff(s,rx.all,ix.all);
    Standard_Newton_Convolutions.Minus(s.vy);
    Standard_Series_Matrix_Solvers.Solve_by_lufco(s.vm,s.vy,ipvt,rcond,wrk);
    if scaledx
     then Standard_Newton_Convolutions.Power_Divide(s.vy,1.0);
    end if;
    Standard_Coefficient_Convolutions.Delinearize(s.vy,s.yv);
    absdx := Standard_Newton_Convolutions.Max(s.yv);
    Standard_Newton_Convolutions.Update(scf,s.yv);
  end LU_Newton_Step;

  procedure LU_Newton_Step
              ( deg : in integer32;
                s : in Standard_Coefficient_Convolutions.Link_to_System;
                scf : in Standard_Complex_VecVecs.VecVec;
                rx,ix : in Standard_Floating_VecVecs.Link_to_VecVec;
                absdx,rcond : out double_float;
                ipvt : out Standard_Integer_Vectors.Vector;
                wrk : in Standard_Complex_Vectors.Link_to_Vector;
                scaledx : in boolean := true;
                vrblvl : in integer32 := 0 ) is
  begin
    if vrblvl > 0 then
      put_line("-> in newton_coefficient_convolutions.LU_Newton_Step 8 ...");
    end if;
    Standard_Vector_Splitters.Complex_Parts(deg,scf,rx,ix);
    Standard_Coefficient_Convolutions.Compute(deg,s.rpwt,s.ipwt,s.mxe,rx,ix);
    Standard_Coefficient_Convolutions.EvalDiff(deg,s,rx.all,ix.all);
    Standard_Newton_Convolutions.Minus(deg,s.vy);
    Standard_Series_Matrix_Solvers.Solve_by_lufco
      (deg,s.vm,s.vy,ipvt,rcond,wrk);
    if scaledx
     then Standard_Newton_Convolutions.Power_Divide(s.vy,1.0);
    end if;
    Standard_Coefficient_Convolutions.Delinearize(deg,s.vy,s.yv);
    absdx := Standard_Newton_Convolutions.Max(deg,s.yv);
    Standard_Newton_Convolutions.Update(deg,scf,s.yv);
  end LU_Newton_Step;

  procedure LU_Newton_Step
              ( file : in file_type;
                s : in Standard_Coefficient_Convolutions.Link_to_System;
                scf : in Standard_Complex_VecVecs.VecVec;
                rx,ix : in Standard_Floating_VecVecs.Link_to_VecVec;
                absdx,rcond : out double_float;
                ipvt : out Standard_Integer_Vectors.Vector;
                wrk : in Standard_Complex_Vectors.Link_to_Vector;
                scaledx : in boolean := true;
                vrblvl : in integer32 := 0 ) is
  begin
    if vrblvl > 0 then
      put_line("-> in newton_coefficient_convolutions.LU_Newton_Step 9 ...");
    end if;
    put_line(file,"scf :"); put_line(file,scf);
    Standard_Vector_Splitters.Complex_Parts(scf,rx,ix);
    Standard_Coefficient_Convolutions.Compute(s.rpwt,s.ipwt,s.mxe,rx,ix);
    Standard_Coefficient_Convolutions.EvalDiff(s,rx.all,ix.all);
    put_line(file,"vy :"); put_line(file,s.vy);
    Standard_Newton_Convolutions.Minus(s.vy);
    Standard_Series_Matrix_Solvers.Solve_by_lufco(s.vm,s.vy,ipvt,rcond,wrk);
    put_line(file,"dx :"); put_line(file,s.vy);
    if scaledx then
      Standard_Newton_Convolutions.Power_Divide(s.vy,1.0);
      put_line(file,"scaled dx :"); put_line(file,s.vy);
    end if;
    Standard_Coefficient_Convolutions.Delinearize(s.vy,s.yv);
    absdx := Standard_Newton_Convolutions.Max(s.yv);
    put(file,"max |dx| :"); put(file,absdx,3); new_line(file);
    Standard_Newton_Convolutions.Update(scf,s.yv);
  end LU_Newton_Step;

  procedure LU_Newton_Step
              ( file : in file_type; deg : in integer32;
                s : in Standard_Coefficient_Convolutions.Link_to_System;
                scf : in Standard_Complex_VecVecs.VecVec;
                rx,ix : in Standard_Floating_VecVecs.Link_to_VecVec;
                absdx,rcond : out double_float;
                ipvt : out Standard_Integer_Vectors.Vector;
                wrk : in Standard_Complex_Vectors.Link_to_Vector;
                scaledx : in boolean := true;
                vrblvl : in integer32 := 0 ) is
  begin
    if vrblvl > 0 then
      put_line("-> in newton_coefficient_convolutions.LU_Newton_Step 10 ...");
    end if;
    put_line(file,"scf :"); put_line(file,scf);
    Standard_Vector_Splitters.Complex_Parts(deg,scf,rx,ix);
    Standard_Coefficient_Convolutions.Compute(deg,s.rpwt,s.ipwt,s.mxe,rx,ix);
    Standard_Coefficient_Convolutions.EvalDiff(deg,s,rx.all,ix.all);
    put_line(file,"vy :");
    for k in 0..deg loop
      put_line(file,s.vy(k)); new_line(file);
    end loop;
    Standard_Newton_Convolutions.Minus(deg,s.vy);
    Standard_Series_Matrix_Solvers.Solve_by_lufco
      (deg,s.vm,s.vy,ipvt,rcond,wrk);
    put_line(file,"dx :");
    for k in 0..deg loop
      put_line(file,s.vy(k)); new_line(file);
    end loop;
    if scaledx then
      Standard_Newton_Convolutions.Power_Divide(s.vy,1.0);
      put_line(file,"scaled dx :"); put_line(file,s.vy);
    end if;
    Standard_Coefficient_Convolutions.Delinearize(deg,s.vy,s.yv);
    absdx := Standard_Newton_Convolutions.Max(deg,s.yv);
    put(file,"max |dx| :"); put(file,absdx,3); new_line(file);
    Standard_Newton_Convolutions.Update(deg,scf,s.yv);
  end LU_Newton_Step;

  procedure LU_Newton_Step
              ( s : in DoblDobl_Coefficient_Convolutions.Link_to_System;
                scf : in DoblDobl_Complex_VecVecs.VecVec;
                rhx,ihx : in Standard_Floating_VecVecs.Link_to_VecVec;
                rlx,ilx : in Standard_Floating_VecVecs.Link_to_VecVec;
                absdx,rcond  : out double_double;
                ipvt : out Standard_Integer_Vectors.Vector;
                wrk : in DoblDobl_Complex_Vectors.Link_to_Vector;
                scaledx : in boolean := true;
                vrblvl : in integer32 := 0 ) is

    fac : constant double_double := create(1.0);

  begin
    if vrblvl > 0 then
      put_line("-> in newton_coefficient_convolutions.LU_Newton_Step 7 ...");
    end if;
    DoblDobl_Vector_Splitters.Complex_Parts(scf,rhx,ihx,rlx,ilx);
    DoblDobl_Coefficient_Convolutions.Compute
      (s.rhpwt,s.ihpwt,s.rlpwt,s.ilpwt,s.mxe,rhx,ihx,rlx,ilx);
    DoblDobl_Coefficient_Convolutions.EvalDiff
      (s,rhx.all,ihx.all,rlx.all,ilx.all);
    DoblDobl_Newton_Convolutions.Minus(s.vy);
    DoblDobl_Series_Matrix_Solvers.Solve_by_lufco(s.vm,s.vy,ipvt,rcond,wrk);
    if scaledx
     then DoblDobl_Newton_Convolutions.Power_Divide(s.vy,fac);
    end if;
    DoblDobl_Speelpenning_Convolutions.Delinearize(s.vy,s.yv);
    absdx := DoblDobl_Newton_Convolutions.Max(s.vy);
    DoblDobl_Newton_Convolutions.Update(scf,s.yv);
  end LU_Newton_Step;

  procedure LU_Newton_Step
              ( file : in file_type;
                s : in DoblDobl_Coefficient_Convolutions.Link_to_System;
                scf : in DoblDobl_Complex_VecVecs.VecVec;
                rhx,ihx : in Standard_Floating_VecVecs.Link_to_VecVec;
                rlx,ilx : in Standard_Floating_VecVecs.Link_to_VecVec;
                absdx,rcond  : out double_double;
                ipvt : out Standard_Integer_Vectors.Vector;
                wrk : in DoblDobl_Complex_Vectors.Link_to_Vector;
                scaledx : in boolean := true;
                vrblvl : in integer32 := 0 ) is

    fac : constant double_double := create(1.0);

  begin
    if vrblvl > 0 then
      put_line("-> in newton_coefficient_convolutions.LU_Newton_Step 8 ...");
    end if;
    put_line(file,"scf :"); put_line(file,scf);
    DoblDobl_Vector_Splitters.Complex_Parts(scf,rhx,ihx,rlx,ilx);
    DoblDobl_Coefficient_Convolutions.Compute
      (s.rhpwt,s.ihpwt,s.rlpwt,s.ilpwt,s.mxe,rhx,ihx,rlx,ilx);
    DoblDobl_Coefficient_Convolutions.EvalDiff
      (s,rhx.all,ihx.all,rlx.all,ilx.all);
    put_line(file,"vy :"); put_line(file,s.vy);
    DoblDobl_Newton_Convolutions.Minus(s.vy);
    DoblDobl_Series_Matrix_Solvers.Solve_by_lufco(s.vm,s.vy,ipvt,rcond,wrk);
    put_line(file,"dx :"); put_line(file,s.vy);
    if scaledx then
      DoblDobl_Newton_Convolutions.Power_Divide(s.vy,fac);
      put_line(file,"scaled dx :"); put_line(file,s.vy);
    end if;
    DoblDobl_Speelpenning_Convolutions.Delinearize(s.vy,s.yv);
    absdx := DoblDobl_Newton_Convolutions.Max(s.vy);
    put(file,"max |dx| :"); put(file,absdx,3); new_line(file);
    DoblDobl_Newton_Convolutions.Update(scf,s.yv);
  end LU_Newton_Step;

-- ONE NEWTON STEP WITH QR :

  procedure QR_Newton_Step
              ( s : in Standard_Coefficient_Convolutions.Link_to_System;
                scf,dx,xd : in Standard_Complex_VecVecs.VecVec;
                rx,ix : in Standard_Floating_VecVecs.Link_to_VecVec;
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
      put_line("-> in newton_coefficient_convolutions.QR_Newton_Step 1 ...");
    end if;
    Standard_Vector_Splitters.Complex_Parts(scf,rx,ix);
    Standard_Coefficient_Convolutions.Compute(s.rpwt,s.ipwt,s.mxe,rx,ix);
    Standard_Coefficient_Convolutions.EvalDiff(s,rx.all,ix.all);
    Standard_Newton_Convolutions.Minus(s.vy);
    Standard_Series_Matrix_Solvers.Solve_by_QRLS
      (s.vm,s.vy,xd,qraux,w1,w2,w3,w4,w5,ipvt,info,wrk);
    if scaledx
     then Standard_Newton_Convolutions.Power_Divide(xd,1.0);
    end if;
    Standard_Coefficient_Convolutions.Delinearize(xd,dx);
    absdx := Standard_Newton_Convolutions.Max(dx);
    Standard_Newton_Convolutions.Update(scf,dx);
  end QR_Newton_Step;

  procedure QR_Newton_Step
              ( deg : in integer32;
                s : in Standard_Coefficient_Convolutions.Link_to_System;
                scf,dx,xd : in Standard_Complex_VecVecs.VecVec;
                rx,ix : in Standard_Floating_VecVecs.Link_to_VecVec;
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
      put_line("-> in newton_coefficient_convolutions.QR_Newton_Step 2 ...");
    end if;
    Standard_Vector_Splitters.Complex_Parts(deg,scf,rx,ix);
    Standard_Coefficient_Convolutions.Compute(deg,s.rpwt,s.ipwt,s.mxe,rx,ix);
    Standard_Coefficient_Convolutions.EvalDiff(deg,s,rx.all,ix.all);
    Standard_Newton_Convolutions.Minus(deg,s.vy);
    Standard_Series_Matrix_Solvers.Solve_by_QRLS
      (deg,s.vm,s.vy,xd,qraux,w1,w2,w3,w4,w5,ipvt,info,wrk);
    if scaledx
     then Standard_Newton_Convolutions.Power_Divide(xd,1.0);
    end if;
    Standard_Coefficient_Convolutions.Delinearize(deg,xd,dx);
    absdx := Standard_Newton_Convolutions.Max(deg,dx);
    Standard_Newton_Convolutions.Update(deg,scf,dx);
  end QR_Newton_Step;

  procedure QR_Newton_Step
              ( file : in file_type;
                s : in Standard_Coefficient_Convolutions.Link_to_System;
                scf,dx,xd : in Standard_Complex_VecVecs.VecVec;
                rx,ix : in Standard_Floating_VecVecs.Link_to_VecVec;
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
      put_line("-> in newton_coefficient_convolutions.QR_Newton_Step 3 ...");
    end if;
    put_line(file,"scf :"); put_line(file,scf);
    Standard_Vector_Splitters.Complex_Parts(scf,rx,ix);
    Standard_Coefficient_Convolutions.Compute(s.rpwt,s.ipwt,s.mxe,rx,ix);
    Standard_Coefficient_Convolutions.EvalDiff(s,rx.all,ix.all);
    put_line(file,"vy :"); put_line(file,s.vy);
    Standard_Newton_Convolutions.Minus(s.vy);
    Standard_Series_Matrix_Solvers.Solve_by_QRLS
      (s.vm,s.vy,xd,qraux,w1,w2,w3,w4,w5,ipvt,info,wrk);
    put_line(file,"dx :"); put_line(file,xd);
    if scaledx then
      Standard_Newton_Convolutions.Power_Divide(xd,1.0);
      put(file,"scaled dx :"); put_line(file,xd);
    end if;
    Standard_Coefficient_Convolutions.Delinearize(xd,dx);
    absdx := Standard_Newton_Convolutions.Max(dx);
    put(file,"max |dx| :"); put(file,absdx,3); new_line(file);
    Standard_Newton_Convolutions.Update(scf,dx);
  end QR_Newton_Step;

  procedure QR_Newton_Step
              ( file : in file_type; deg : in integer32;
                s : in Standard_Coefficient_Convolutions.Link_to_System;
                scf,dx,xd : in Standard_Complex_VecVecs.VecVec;
                rx,ix : in Standard_Floating_VecVecs.Link_to_VecVec;
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
      put_line("-> in newton_coefficient_convolutions.QR_Newton_Step 4 ...");
    end if;
    put_line(file,"scf :"); put_line(file,scf);
    Standard_Vector_Splitters.Complex_Parts(deg,scf,rx,ix);
    Standard_Coefficient_Convolutions.Compute(deg,s.rpwt,s.ipwt,s.mxe,rx,ix);
    Standard_Coefficient_Convolutions.EvalDiff(deg,s,rx.all,ix.all);
    put_line(file,"vy :");
    for k in 0..deg loop
      put_line(file,s.vy(k)); new_line(file);
    end loop;
    Standard_Newton_Convolutions.Minus(deg,s.vy);
    Standard_Series_Matrix_Solvers.Solve_by_QRLS
      (deg,s.vm,s.vy,xd,qraux,w1,w2,w3,w4,w5,ipvt,info,wrk);
    put_line(file,"dx :");
    for k in 0..deg loop
      put_line(file,xd(k)); new_line(file);
    end loop;
    if scaledx then
      Standard_Newton_Convolutions.Power_Divide(xd,1.0);
      put(file,"scaled dx :"); put_line(file,xd);
    end if;
    Standard_Coefficient_Convolutions.Delinearize(deg,xd,dx);
    absdx := Standard_Newton_Convolutions.Max(deg,dx);
    put(file,"max |dx| :"); put(file,absdx,3); new_line(file);
    Standard_Newton_Convolutions.Update(deg,scf,dx);
  end QR_Newton_Step;

  procedure QR_Newton_Step
              ( s : in DoblDobl_Coefficient_Convolutions.Link_to_System;
                scf,dx,xd : in DoblDobl_Complex_VecVecs.VecVec;
                rhx,ihx : in Standard_Floating_VecVecs.Link_to_VecVec;
                rlx,ilx : in Standard_Floating_VecVecs.Link_to_VecVec;
                absdx : out double_double;
                qraux : out DoblDobl_Complex_Vectors.Vector;
                w1,w2,w3,w4,w5 : in out DoblDobl_Complex_Vectors.Vector;
                info : out integer32;
                ipvt : out Standard_Integer_Vectors.Vector;
                wrk : in DoblDobl_Complex_Vectors.Link_to_Vector;
                scaledx : in boolean := true;
                vrblvl : in integer32 := 0 ) is

    fac : constant double_double := create(1.0);

  begin
    if vrblvl > 0 then
      put_line("-> in newton_coefficient_convolutions.QR_Newton_Step 5 ...");
    end if;
    DoblDobl_Vector_Splitters.Complex_Parts(scf,rhx,ihx,rlx,ilx);
    DoblDobl_Coefficient_Convolutions.Compute
      (s.rhpwt,s.ihpwt,s.rlpwt,s.ilpwt,s.mxe,rhx,ihx,rlx,ilx);
    DoblDobl_Coefficient_Convolutions.EvalDiff
      (s,rhx.all,ihx.all,rlx.all,ilx.all);
    DoblDobl_Newton_Convolutions.Minus(s.vy);
    DoblDobl_Series_Matrix_Solvers.Solve_by_QRLS
      (s.vm,s.vy,xd,qraux,w1,w2,w3,w4,w5,ipvt,info,wrk);
    if scaledx
     then DoblDobl_Newton_Convolutions.Power_Divide(xd,fac);
    end if;
    DoblDobl_Speelpenning_Convolutions.Delinearize(xd,dx);
    absdx := DoblDobl_Newton_Convolutions.Max(dx);
    DoblDobl_Newton_Convolutions.Update(scf,dx);
  end QR_Newton_Step;

  procedure QR_Newton_Step
              ( file : in file_type;
                s : in DoblDobl_Coefficient_Convolutions.Link_to_System;
                scf,dx,xd : in DoblDobl_Complex_VecVecs.VecVec;
                rhx,ihx : in Standard_Floating_VecVecs.Link_to_VecVec;
                rlx,ilx : in Standard_Floating_VecVecs.Link_to_VecVec;
                absdx : out double_double;
                qraux : out DoblDobl_Complex_Vectors.Vector;
                w1,w2,w3,w4,w5 : in out DoblDobl_Complex_Vectors.Vector;
                info : out integer32;
                ipvt : out Standard_Integer_Vectors.Vector;
                wrk : in DoblDobl_Complex_Vectors.Link_to_Vector;
                scaledx : in boolean := true;
                vrblvl : in integer32 := 0 ) is

    fac : constant double_double := create(1.0);

  begin
    if vrblvl > 0 then
      put_line("-> in newton_coefficient_convolutions.QR_Newton_Step 6 ...");
    end if;
    put_line(file,"scf :"); put_line(file,scf);
    DoblDobl_Vector_Splitters.Complex_Parts(scf,rhx,ihx,rlx,ilx);
    DoblDobl_Coefficient_Convolutions.Compute
      (s.rhpwt,s.ihpwt,s.rlpwt,s.ilpwt,s.mxe,rhx,ihx,rlx,ilx);
    DoblDobl_Coefficient_Convolutions.EvalDiff
      (s,rhx.all,ihx.all,rlx.all,ilx.all);
    put_line(file,"vy :"); put_line(file,s.vy);
    DoblDobl_Newton_Convolutions.Minus(s.vy);
    DoblDobl_Series_Matrix_Solvers.Solve_by_QRLS
      (s.vm,s.vy,xd,qraux,w1,w2,w3,w4,w5,ipvt,info,wrk);
    if scaledx then
      DoblDobl_Newton_Convolutions.Power_Divide(xd,fac);
      put(file,"scaled dx :"); put_line(file,xd);
    end if;
    DoblDobl_Speelpenning_Convolutions.Delinearize(xd,dx);
    absdx := DoblDobl_Newton_Convolutions.Max(dx);
    put(file,"max |dx| :"); put(file,absdx,3); new_line(file);
    DoblDobl_Newton_Convolutions.Update(scf,dx);
  end QR_Newton_Step;

-- ONE NEWTON STEP WITH SVD :

  procedure SVD_Newton_Step
              ( s : in Standard_Coefficient_Convolutions.Link_to_System;
                scf,dx,xd : in Standard_Complex_VecVecs.VecVec;
                rx,ix : in Standard_Floating_VecVecs.Link_to_VecVec;
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
      put_line("-> in newton_coefficient_convolutions.SVD_Newton_Step 1 ...");
    end if;
    Standard_Vector_Splitters.Complex_Parts(scf,rx,ix);
    Standard_Coefficient_Convolutions.Compute(s.rpwt,s.ipwt,s.mxe,rx,ix);
    Standard_Coefficient_Convolutions.EvalDiff(s,rx.all,ix.all);
    Standard_Newton_Convolutions.Minus(s.vy);
    Standard_Series_Matrix_Solvers.Solve_by_SVD
      (s.vm,s.vy,xd,svl,U,V,info,rcond,ewrk,wrkv);
    if scaledx
     then Standard_Newton_Convolutions.Power_Divide(xd,1.0);
    end if;
    Standard_Coefficient_Convolutions.Delinearize(xd,dx);
    absdx := Standard_Newton_Convolutions.Max(dx);
    Standard_Newton_Convolutions.Update(scf,dx);
  end SVD_Newton_Step;

  procedure SVD_Newton_Step
              ( deg : in integer32;
                s : in Standard_Coefficient_Convolutions.Link_to_System;
                scf,dx,xd : in Standard_Complex_VecVecs.VecVec;
                rx,ix : in Standard_Floating_VecVecs.Link_to_VecVec;
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
      put_line("-> in newton_coefficient_convolutions.SVD_Newton_Step 2 ...");
    end if;
    Standard_Vector_Splitters.Complex_Parts(deg,scf,rx,ix);
    Standard_Coefficient_Convolutions.Compute(deg,s.rpwt,s.ipwt,s.mxe,rx,ix);
    Standard_Coefficient_Convolutions.EvalDiff(deg,s,rx.all,ix.all);
    Standard_Newton_Convolutions.Minus(deg,s.vy);
    Standard_Series_Matrix_Solvers.Solve_by_SVD
      (deg,s.vm,s.vy,xd,svl,U,V,info,rcond,ewrk,wrkv);
    if scaledx
     then Standard_Newton_Convolutions.Power_Divide(xd,1.0);
    end if;
    Standard_Coefficient_Convolutions.Delinearize(deg,xd,dx);
    absdx := Standard_Newton_Convolutions.Max(deg,dx);
    Standard_Newton_Convolutions.Update(deg,scf,dx);
  end SVD_Newton_Step;

  procedure SVD_Newton_Step
              ( file : in file_type;
                s : in Standard_Coefficient_Convolutions.Link_to_System;
                scf,dx,xd : in Standard_Complex_VecVecs.VecVec;
                rx,ix : in Standard_Floating_VecVecs.Link_to_VecVec;
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
      put_line("-> in newton_coefficient_convolutions.SVD_Newton_Step 3 ...");
    end if;
    put_line(file,"scf :"); put_line(file,scf);
    Standard_Vector_Splitters.Complex_Parts(scf,rx,ix);
    Standard_Coefficient_Convolutions.Compute(s.rpwt,s.ipwt,s.mxe,rx,ix);
    Standard_Coefficient_Convolutions.EvalDiff(s,rx.all,ix.all);
    put_line(file,"vy :"); put_line(file,s.vy);
    Standard_Newton_Convolutions.Minus(s.vy);
    Standard_Series_Matrix_Solvers.Solve_by_SVD
      (s.vm,s.vy,xd,svl,U,V,info,rcond,ewrk,wrkv);
    put_line(file,"dx :"); put_line(file,xd);
    if scaledx then
      Standard_Newton_Convolutions.Power_Divide(xd,1.0);
      put(file,"scaled dx :"); put_line(file,xd);
    end if;
    Standard_Coefficient_Convolutions.Delinearize(xd,dx);
    absdx := Standard_Newton_Convolutions.Max(dx);
    put(file,"max |dx| : "); put(file,absdx,3); new_line(file);
    Standard_Newton_Convolutions.Update(scf,dx);
  end SVD_Newton_Step;

  procedure SVD_Newton_Step
              ( file : in file_type; deg : in integer32;
                s : in Standard_Coefficient_Convolutions.Link_to_System;
                scf,dx,xd : in Standard_Complex_VecVecs.VecVec;
                rx,ix : in Standard_Floating_VecVecs.Link_to_VecVec;
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
      put_line("-> in newton_coefficient_convolutions.SVD_Newton_Step 4 ...");
    end if;
    put_line(file,"scf :"); put_line(file,scf);
    Standard_Vector_Splitters.Complex_Parts(deg,scf,rx,ix);
    Standard_Coefficient_Convolutions.Compute(deg,s.rpwt,s.ipwt,s.mxe,rx,ix);
    Standard_Coefficient_Convolutions.EvalDiff(deg,s,rx.all,ix.all);
    put_line(file,"vy :");
    for k in 0..deg loop
      put_line(file,s.vy(k)); new_line(file);
    end loop;
    Standard_Newton_Convolutions.Minus(deg,s.vy);
    Standard_Series_Matrix_Solvers.Solve_by_SVD
      (deg,s.vm,s.vy,xd,svl,U,V,info,rcond,ewrk,wrkv);
    put_line(file,"dx :");
    for k in 0..deg loop
      put_line(file,xd(k)); new_line(file);
    end loop;
    if scaledx then
      Standard_Newton_Convolutions.Power_Divide(xd,1.0);
      put(file,"scaled dx :"); put_line(file,xd);
    end if;
    Standard_Coefficient_Convolutions.Delinearize(deg,xd,dx);
    absdx := Standard_Newton_Convolutions.Max(deg,dx);
    put(file,"max |dx| : "); put(file,absdx,3); new_line(file);
    Standard_Newton_Convolutions.Update(deg,scf,dx);
  end SVD_Newton_Step;

  procedure SVD_Newton_Step
              ( s : in DoblDobl_Coefficient_Convolutions.Link_to_System;
                scf,dx,xd : in DoblDobl_Complex_VecVecs.VecVec;
                rhx,ihx : in Standard_Floating_VecVecs.Link_to_VecVec;
                rlx,ilx : in Standard_Floating_VecVecs.Link_to_VecVec;
                absdx : out double_double;
                svl : out DoblDobl_Complex_Vectors.Vector;
                U,V : out DoblDobl_Complex_Matrices.Matrix;
                info : out integer32; rcond : out double_double;
                ewrk : in DoblDobl_Complex_Vectors.Link_to_Vector;
                wrkv : in DoblDobl_Complex_Vectors.Link_to_Vector;
                scaledx : in boolean := true;
                vrblvl : in integer32 := 0 ) is

    fac : constant double_double := create(1.0);

  begin
    if vrblvl > 0 then
      put_line("-> in newton_coefficient_convolutions.SVD_Newton_Step 5 ...");
    end if;
    DoblDobl_Vector_Splitters.Complex_Parts(scf,rhx,ihx,rlx,ilx);
    DoblDobl_Coefficient_Convolutions.Compute
      (s.rhpwt,s.ihpwt,s.rlpwt,s.ilpwt,s.mxe,rhx,ihx,rlx,ilx);
    DoblDobl_Coefficient_Convolutions.EvalDiff
      (s,rhx.all,ihx.all,rlx.all,ilx.all);
    DoblDobl_Newton_Convolutions.Minus(s.vy);
    DoblDobl_Series_Matrix_Solvers.Solve_by_SVD
      (s.vm,s.vy,xd,svl,U,V,info,rcond,ewrk,wrkv);
    if scaledx
     then DoblDobl_Newton_Convolutions.Power_Divide(xd,fac);
    end if;
    DoblDobl_Speelpenning_Convolutions.Delinearize(xd,dx);
    absdx := DoblDobl_Newton_Convolutions.Max(dx);
    DoblDobl_Newton_Convolutions.Update(scf,dx);
  end SVD_Newton_Step;

  procedure SVD_Newton_Step
              ( file : in file_type;
                s : in DoblDobl_Coefficient_Convolutions.Link_to_System;
                scf,dx,xd : in DoblDobl_Complex_VecVecs.VecVec;
                rhx,ihx : in Standard_Floating_VecVecs.Link_to_VecVec;
                rlx,ilx : in Standard_Floating_VecVecs.Link_to_VecVec;
                absdx : out double_double;
                svl : out DoblDobl_Complex_Vectors.Vector;
                U,V : out DoblDobl_Complex_Matrices.Matrix;
                info : out integer32; rcond : out double_double;
                ewrk : in DoblDobl_Complex_Vectors.Link_to_Vector;
                wrkv : in DoblDobl_Complex_Vectors.Link_to_Vector;
                scaledx : in boolean := true;
                vrblvl : in integer32 := 0 ) is

    fac : constant double_double := create(1.0);

  begin
    if vrblvl > 0 then
      put_line("-> in newton_coefficient_convolutions.SVD_Newton_Step 6 ...");
    end if;
    put_line(file,"scf :"); put_line(file,scf);
    DoblDobl_Vector_Splitters.Complex_Parts(scf,rhx,ihx,rlx,ilx);
    DoblDobl_Coefficient_Convolutions.Compute
      (s.rhpwt,s.ihpwt,s.rlpwt,s.ilpwt,s.mxe,rhx,ihx,rlx,ilx);
    DoblDobl_Coefficient_Convolutions.EvalDiff
      (s,rhx.all,ihx.all,rlx.all,ilx.all);
    put_line(file,"vy :"); put_line(file,s.vy);
    DoblDobl_Newton_Convolutions.Minus(s.vy);
    DoblDobl_Series_Matrix_Solvers.Solve_by_SVD
      (s.vm,s.vy,xd,svl,U,V,info,rcond,ewrk,wrkv);
    if scaledx then
      DoblDobl_Newton_Convolutions.Power_Divide(xd,fac);
      put(file,"scaled dx :"); put_line(file,xd);
    end if;
    DoblDobl_Speelpenning_Convolutions.Delinearize(xd,dx);
    absdx := DoblDobl_Newton_Convolutions.Max(dx);
    put(file,"max |dx| : "); put(file,absdx,3); new_line(file);
    DoblDobl_Newton_Convolutions.Update(scf,dx);
  end SVD_Newton_Step;

end Newton_Coefficient_Convolutions;
