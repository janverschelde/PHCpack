with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;
with Double_Double_Numbers_io;           use Double_Double_Numbers_io;
with Triple_Double_Numbers_io;           use Triple_Double_Numbers_io;
with Quad_Double_Numbers_io;             use Quad_Double_Numbers_io;
with Penta_Double_Numbers_io;            use Penta_Double_Numbers_io;
with Octo_Double_Numbers_io;             use Octo_Double_Numbers_io;
with Deca_Double_Numbers_io;             use Deca_Double_Numbers_io;
with Standard_Newton_Convolutions;
with DoblDobl_Newton_Convolutions;
with TripDobl_Newton_Convolutions;
with QuadDobl_Newton_Convolutions;
with PentDobl_Newton_Convolutions;
with OctoDobl_Newton_Convolutions;
with DecaDobl_Newton_Convolutions;
with Multitasked_AlgoDiff_Convolutions;  use Multitasked_AlgoDiff_Convolutions;
with Multitasked_Series_Linearization;   use Multitasked_Series_Linearization;

package body Multitasked_Newton_Convolutions is

-- ONE NEWTON STEP WITH LU WITHOUT CONDITION NUMBER ESTIMATE :

  procedure Multitasked_LU_Newton_Step
              ( nbt : in integer32;
                s : in Standard_Speelpenning_Convolutions.Link_to_System;
                x : in Standard_Complex_VecVecs.VecVec;
                absdx : out double_float; info : out integer32;
                ipvt : out Standard_Integer_Vectors.Vector;
                wrk : in Standard_Complex_VecVecs.VecVec;
                output : in boolean := false ) is
  begin
    Standard_Multitasked_EvalDiff(nbt,s.crc,x,s.mxe,s.pwt,s.vy,s.vm,output);
    Standard_Newton_Convolutions.Minus(s.vy);
    Multitasked_Solve_by_lufac(nbt,s.vm,s.vy,ipvt,info,wrk,output);
    Standard_Speelpenning_Convolutions.Delinearize(s.vy,s.yv);
    absdx := Standard_Newton_Convolutions.Max(s.yv);
    Standard_Newton_Convolutions.Update(x,s.yv);
  end Multitasked_LU_Newton_Step;

  procedure Multitasked_LU_Newton_Step
              ( nbt : in integer32;
                s : in DoblDobl_Speelpenning_Convolutions.Link_to_System;
                x : in DoblDobl_Complex_VecVecs.VecVec;
                absdx : out double_double; info : out integer32;
                ipvt : out Standard_Integer_Vectors.Vector;
                wrk : in DoblDobl_Complex_VecVecs.VecVec;
                output : in boolean := false ) is
  begin
    DoblDobl_Multitasked_EvalDiff(nbt,s.crc,x,s.mxe,s.pwt,s.vy,s.vm,output);
    DoblDobl_Newton_Convolutions.Minus(s.vy);
    Multitasked_Solve_by_lufac(nbt,s.vm,s.vy,ipvt,info,wrk,output);
    DoblDobl_Speelpenning_Convolutions.Delinearize(s.vy,s.yv);
    absdx := DoblDobl_Newton_Convolutions.Max(s.yv);
    DoblDobl_Newton_Convolutions.Update(x,s.yv);
  end Multitasked_LU_Newton_Step;

  procedure Multitasked_LU_Newton_Step
              ( nbt : in integer32;
                s : in TripDobl_Speelpenning_Convolutions.Link_to_System;
                x : in TripDobl_Complex_VecVecs.VecVec;
                absdx : out triple_double; info : out integer32;
                ipvt : out Standard_Integer_Vectors.Vector;
                wrk : in TripDobl_Complex_VecVecs.VecVec;
                output : in boolean := false ) is
  begin
    TripDobl_Multitasked_EvalDiff(nbt,s.crc,x,s.mxe,s.pwt,s.vy,s.vm,output);
    TripDobl_Newton_Convolutions.Minus(s.vy);
    Multitasked_Solve_by_lufac(nbt,s.vm,s.vy,ipvt,info,wrk,output);
    TripDobl_Speelpenning_Convolutions.Delinearize(s.vy,s.yv);
    absdx := TripDobl_Newton_Convolutions.Max(s.yv);
    TripDobl_Newton_Convolutions.Update(x,s.yv);
  end Multitasked_LU_Newton_Step;

  procedure Multitasked_LU_Newton_Step
              ( nbt : in integer32;
                s : in QuadDobl_Speelpenning_Convolutions.Link_to_System;
                x : in QuadDobl_Complex_VecVecs.VecVec;
                absdx : out quad_double; info : out integer32;
                ipvt : out Standard_Integer_Vectors.Vector;
                wrk : in QuadDobl_Complex_VecVecs.VecVec;
                output : in boolean := false ) is
  begin
    QuadDobl_Multitasked_EvalDiff(nbt,s.crc,x,s.mxe,s.pwt,s.vy,s.vm,output);
    QuadDobl_Newton_Convolutions.Minus(s.vy);
    Multitasked_Solve_by_lufac(nbt,s.vm,s.vy,ipvt,info,wrk,output);
    QuadDobl_Speelpenning_Convolutions.Delinearize(s.vy,s.yv);
    absdx := QuadDobl_Newton_Convolutions.Max(s.yv);
    QuadDobl_Newton_Convolutions.Update(x,s.yv);
  end Multitasked_LU_Newton_Step;

  procedure Multitasked_LU_Newton_Step
              ( nbt : in integer32;
                s : in PentDobl_Speelpenning_Convolutions.Link_to_System;
                x : in PentDobl_Complex_VecVecs.VecVec;
                absdx : out penta_double; info : out integer32;
                ipvt : out Standard_Integer_Vectors.Vector;
                wrk : in PentDobl_Complex_VecVecs.VecVec;
                output : in boolean := false ) is
  begin
    PentDobl_Multitasked_EvalDiff(nbt,s.crc,x,s.mxe,s.pwt,s.vy,s.vm,output);
    PentDobl_Newton_Convolutions.Minus(s.vy);
    Multitasked_Solve_by_lufac(nbt,s.vm,s.vy,ipvt,info,wrk,output);
    PentDobl_Speelpenning_Convolutions.Delinearize(s.vy,s.yv);
    absdx := PentDobl_Newton_Convolutions.Max(s.yv);
    PentDobl_Newton_Convolutions.Update(x,s.yv);
  end Multitasked_LU_Newton_Step;

  procedure Multitasked_LU_Newton_Step
              ( nbt : in integer32;
                s : in OctoDobl_Speelpenning_Convolutions.Link_to_System;
                x : in OctoDobl_Complex_VecVecs.VecVec;
                absdx : out octo_double; info : out integer32;
                ipvt : out Standard_Integer_Vectors.Vector;
                wrk : in OctoDobl_Complex_VecVecs.VecVec;
                output : in boolean := false ) is
  begin
    OctoDobl_Multitasked_EvalDiff(nbt,s.crc,x,s.mxe,s.pwt,s.vy,s.vm,output);
    OctoDobl_Newton_Convolutions.Minus(s.vy);
    Multitasked_Solve_by_lufac(nbt,s.vm,s.vy,ipvt,info,wrk,output);
    OctoDobl_Speelpenning_Convolutions.Delinearize(s.vy,s.yv);
    absdx := OctoDobl_Newton_Convolutions.Max(s.yv);
    OctoDobl_Newton_Convolutions.Update(x,s.yv);
  end Multitasked_LU_Newton_Step;

  procedure Multitasked_LU_Newton_Step
              ( nbt : in integer32;
                s : in DecaDobl_Speelpenning_Convolutions.Link_to_System;
                x : in DecaDobl_Complex_VecVecs.VecVec;
                absdx : out deca_double; info : out integer32;
                ipvt : out Standard_Integer_Vectors.Vector;
                wrk : in DecaDobl_Complex_VecVecs.VecVec;
                output : in boolean := false ) is
  begin
    DecaDobl_Multitasked_EvalDiff(nbt,s.crc,x,s.mxe,s.pwt,s.vy,s.vm,output);
    DecaDobl_Newton_Convolutions.Minus(s.vy);
    Multitasked_Solve_by_lufac(nbt,s.vm,s.vy,ipvt,info,wrk,output);
    DecaDobl_Speelpenning_Convolutions.Delinearize(s.vy,s.yv);
    absdx := DecaDobl_Newton_Convolutions.Max(s.yv);
    DecaDobl_Newton_Convolutions.Update(x,s.yv);
  end Multitasked_LU_Newton_Step;

-- ONE NEWTON STEP WITH LU WITH CONDITION NUMBER ESTIMATE :

  procedure Multitasked_LU_Newton_Step
              ( nbt : in integer32;
                s : in Standard_Speelpenning_Convolutions.Link_to_System;
                x : in Standard_Complex_VecVecs.VecVec;
                absdx : out double_float; rcond : out double_float;
                ipvt : out Standard_Integer_Vectors.Vector;
                wrk : in Standard_Complex_VecVecs.VecVec;
                output : in boolean := false ) is
  begin
    Standard_Multitasked_EvalDiff(nbt,s.crc,x,s.mxe,s.pwt,s.vy,s.vm,output);
    Standard_Newton_Convolutions.Minus(s.vy);
    Multitasked_Solve_by_lufco(nbt,s.vm,s.vy,ipvt,rcond,wrk,output);
    Standard_Speelpenning_Convolutions.Delinearize(s.vy,s.yv);
    absdx := Standard_Newton_Convolutions.Max(s.yv);
    Standard_Newton_Convolutions.Update(x,s.yv);
  end Multitasked_LU_Newton_Step;

  procedure Multitasked_LU_Newton_Step
              ( nbt : in integer32;
                s : in DoblDobl_Speelpenning_Convolutions.Link_to_System;
                x : in DoblDobl_Complex_VecVecs.VecVec;
                absdx : out double_double; rcond : out double_double;
                ipvt : out Standard_Integer_Vectors.Vector;
                wrk : in DoblDobl_Complex_VecVecs.VecVec;
                output : in boolean := false ) is
  begin
    DoblDobl_Multitasked_EvalDiff(nbt,s.crc,x,s.mxe,s.pwt,s.vy,s.vm,output);
    DoblDobl_Newton_Convolutions.Minus(s.vy);
    Multitasked_Solve_by_lufco(nbt,s.vm,s.vy,ipvt,rcond,wrk,output);
    DoblDobl_Speelpenning_Convolutions.Delinearize(s.vy,s.yv);
    absdx := DoblDobl_Newton_Convolutions.Max(s.yv);
    DoblDobl_Newton_Convolutions.Update(x,s.yv);
  end Multitasked_LU_Newton_Step;

  procedure Multitasked_LU_Newton_Step
              ( nbt : in integer32;
                s : in TripDobl_Speelpenning_Convolutions.Link_to_System;
                x : in TripDobl_Complex_VecVecs.VecVec;
                absdx : out triple_double; rcond : out triple_double;
                ipvt : out Standard_Integer_Vectors.Vector;
                wrk : in TripDobl_Complex_VecVecs.VecVec;
                output : in boolean := false ) is
  begin
    TripDobl_Multitasked_EvalDiff(nbt,s.crc,x,s.mxe,s.pwt,s.vy,s.vm,output);
    TripDobl_Newton_Convolutions.Minus(s.vy);
    Multitasked_Solve_by_lufco(nbt,s.vm,s.vy,ipvt,rcond,wrk,output);
    TripDobl_Speelpenning_Convolutions.Delinearize(s.vy,s.yv);
    absdx := TripDobl_Newton_Convolutions.Max(s.yv);
    TripDobl_Newton_Convolutions.Update(x,s.yv);
  end Multitasked_LU_Newton_Step;

  procedure Multitasked_LU_Newton_Step
              ( nbt : in integer32;
                s : in QuadDobl_Speelpenning_Convolutions.Link_to_System;
                x : in QuadDobl_Complex_VecVecs.VecVec;
                absdx : out quad_double; rcond : out quad_double;
                ipvt : out Standard_Integer_Vectors.Vector;
                wrk : in QuadDobl_Complex_VecVecs.VecVec;
                output : in boolean := false ) is
  begin
    QuadDobl_Multitasked_EvalDiff(nbt,s.crc,x,s.mxe,s.pwt,s.vy,s.vm,output);
    QuadDobl_Newton_Convolutions.Minus(s.vy);
    Multitasked_Solve_by_lufco(nbt,s.vm,s.vy,ipvt,rcond,wrk,output);
    QuadDobl_Speelpenning_Convolutions.Delinearize(s.vy,s.yv);
    absdx := QuadDobl_Newton_Convolutions.Max(s.yv);
    QuadDobl_Newton_Convolutions.Update(x,s.yv);
  end Multitasked_LU_Newton_Step;

  procedure Multitasked_LU_Newton_Step
              ( nbt : in integer32;
                s : in PentDobl_Speelpenning_Convolutions.Link_to_System;
                x : in PentDobl_Complex_VecVecs.VecVec;
                absdx : out penta_double; rcond : out penta_double;
                ipvt : out Standard_Integer_Vectors.Vector;
                wrk : in PentDobl_Complex_VecVecs.VecVec;
                output : in boolean := false ) is
  begin
    PentDobl_Multitasked_EvalDiff(nbt,s.crc,x,s.mxe,s.pwt,s.vy,s.vm,output);
    PentDobl_Newton_Convolutions.Minus(s.vy);
    Multitasked_Solve_by_lufco(nbt,s.vm,s.vy,ipvt,rcond,wrk,output);
    PentDobl_Speelpenning_Convolutions.Delinearize(s.vy,s.yv);
    absdx := PentDobl_Newton_Convolutions.Max(s.yv);
    PentDobl_Newton_Convolutions.Update(x,s.yv);
  end Multitasked_LU_Newton_Step;

  procedure Multitasked_LU_Newton_Step
              ( nbt : in integer32;
                s : in OctoDobl_Speelpenning_Convolutions.Link_to_System;
                x : in OctoDobl_Complex_VecVecs.VecVec;
                absdx : out octo_double; rcond : out octo_double;
                ipvt : out Standard_Integer_Vectors.Vector;
                wrk : in OctoDobl_Complex_VecVecs.VecVec;
                output : in boolean := false ) is
  begin
    OctoDobl_Multitasked_EvalDiff(nbt,s.crc,x,s.mxe,s.pwt,s.vy,s.vm,output);
    OctoDobl_Newton_Convolutions.Minus(s.vy);
    Multitasked_Solve_by_lufco(nbt,s.vm,s.vy,ipvt,rcond,wrk,output);
    OctoDobl_Speelpenning_Convolutions.Delinearize(s.vy,s.yv);
    absdx := OctoDobl_Newton_Convolutions.Max(s.yv);
    OctoDobl_Newton_Convolutions.Update(x,s.yv);
  end Multitasked_LU_Newton_Step;

  procedure Multitasked_LU_Newton_Step
              ( nbt : in integer32;
                s : in DecaDobl_Speelpenning_Convolutions.Link_to_System;
                x : in DecaDobl_Complex_VecVecs.VecVec;
                absdx : out deca_double; rcond : out deca_double;
                ipvt : out Standard_Integer_Vectors.Vector;
                wrk : in DecaDobl_Complex_VecVecs.VecVec;
                output : in boolean := false ) is
  begin
    DecaDobl_Multitasked_EvalDiff(nbt,s.crc,x,s.mxe,s.pwt,s.vy,s.vm,output);
    DecaDobl_Newton_Convolutions.Minus(s.vy);
    Multitasked_Solve_by_lufco(nbt,s.vm,s.vy,ipvt,rcond,wrk,output);
    DecaDobl_Speelpenning_Convolutions.Delinearize(s.vy,s.yv);
    absdx := DecaDobl_Newton_Convolutions.Max(s.yv);
    DecaDobl_Newton_Convolutions.Update(x,s.yv);
  end Multitasked_LU_Newton_Step;

-- SEVERAL NEWTON STEPS WITH LU WITHOUT CONDITION NUMBER ESTIMATE :

  procedure Multitasked_LU_Newton_Steps
              ( nbt : in integer32;
                s : in Standard_Speelpenning_Convolutions.Link_to_System;
                x : in Standard_Complex_VecVecs.VecVec;
                maxit : in integer32; nbrit : out integer32;
                tol : in double_float; absdx : out double_float; 
                fail : out boolean; info : out integer32;
                ipvt : out Standard_Integer_Vectors.Vector;
                wrk : in Standard_Complex_VecVecs.VecVec;
                output : in boolean := false ) is
  begin
    fail := true; nbrit := maxit;
    for k in 1..maxit loop
      Multitasked_LU_Newton_Step(nbt,s,x,absdx,info,ipvt,wrk,output);
      if absdx <= tol
       then fail := false; nbrit := k; exit;
      end if;
    end loop;
  end Multitasked_LU_Newton_Steps;

  procedure Multitasked_LU_Newton_Steps
              ( file : in file_type; nbt : in integer32;
                s : in Standard_Speelpenning_Convolutions.Link_to_System;
                x : in Standard_Complex_VecVecs.VecVec;
                maxit : in integer32; nbrit : out integer32;
                tol : in double_float; absdx : out double_float; 
                fail : out boolean; info : out integer32;
                ipvt : out Standard_Integer_Vectors.Vector;
                wrk : in Standard_Complex_VecVecs.VecVec;
                output : in boolean := false ) is
  begin
    fail := true; nbrit := maxit;
    for k in 1..maxit loop
      put(file,"Step "); put(file,k,1); put_line(file," :");
      Multitasked_LU_Newton_Step(nbt,s,x,absdx,info,ipvt,wrk,output);
      put(file,"  info : "); put(file,info,1);
      put(file,"  absdx :"); put(file,absdx,3); new_line(file);
      if absdx <= tol
       then fail := false; nbrit := k; exit;
      end if;
    end loop;
  end Multitasked_LU_Newton_Steps;

  procedure Multitasked_LU_Newton_Steps
              ( nbt : in integer32;
                s : in DoblDobl_Speelpenning_Convolutions.Link_to_System;
                x : in DoblDobl_Complex_VecVecs.VecVec;
                maxit : in integer32; nbrit : out integer32;
                tol : in double_double; absdx : out double_double; 
                fail : out boolean; info : out integer32;
                ipvt : out Standard_Integer_Vectors.Vector;
                wrk : in DoblDobl_Complex_VecVecs.VecVec;
                output : in boolean := false ) is
  begin
    fail := true; nbrit := maxit;
    for k in 1..maxit loop
      Multitasked_LU_Newton_Step(nbt,s,x,absdx,info,ipvt,wrk,output);
      if absdx <= tol
       then fail := false; nbrit := k; exit;
      end if;
    end loop;
  end Multitasked_LU_Newton_Steps;

  procedure Multitasked_LU_Newton_Steps
              ( file : in file_type; nbt : in integer32;
                s : in DoblDobl_Speelpenning_Convolutions.Link_to_System;
                x : in DoblDobl_Complex_VecVecs.VecVec;
                maxit : in integer32; nbrit : out integer32;
		tol : in double_double; absdx : out double_double; 
		fail : out boolean; info : out integer32;
                ipvt : out Standard_Integer_Vectors.Vector;
                wrk : in DoblDobl_Complex_VecVecs.VecVec;
                output : in boolean := false ) is
  begin
    fail := true; nbrit := maxit;
    for k in 1..maxit loop
      put(file,"Step "); put(file,k,1); put_line(file," :");
      Multitasked_LU_Newton_Step(nbt,s,x,absdx,info,ipvt,wrk,output);
      put(file,"  info : "); put(file,info,1);
      put(file,"  absdx : "); put(file,absdx,3); new_line(file);
      if absdx <= tol
       then fail := false; nbrit := k; exit;
      end if;
    end loop;
  end Multitasked_LU_Newton_Steps;

  procedure Multitasked_LU_Newton_Steps
              ( nbt : in integer32;
                s : in TripDobl_Speelpenning_Convolutions.Link_to_System;
                x : in TripDobl_Complex_VecVecs.VecVec;
                maxit : in integer32; nbrit : out integer32;
                tol : in triple_double; absdx : out triple_double; 
                fail : out boolean; info : out integer32;
                ipvt : out Standard_Integer_Vectors.Vector;
                wrk : in TripDobl_Complex_VecVecs.VecVec;
                output : in boolean := false ) is
  begin
    fail := true; nbrit := maxit;
    for k in 1..maxit loop
      Multitasked_LU_Newton_Step(nbt,s,x,absdx,info,ipvt,wrk,output);
      if absdx <= tol
       then fail := false; nbrit := k; exit;
      end if;
    end loop;
  end Multitasked_LU_Newton_Steps;

  procedure Multitasked_LU_Newton_Steps
              ( file : in file_type; nbt : in integer32;
                s : in TripDobl_Speelpenning_Convolutions.Link_to_System;
                x : in TripDobl_Complex_VecVecs.VecVec;
                maxit : in integer32; nbrit : out integer32;
		tol : in triple_double; absdx : out triple_double; 
		fail : out boolean; info : out integer32;
                ipvt : out Standard_Integer_Vectors.Vector;
                wrk : in TripDobl_Complex_VecVecs.VecVec;
                output : in boolean := false ) is
  begin
    fail := true; nbrit := maxit;
    for k in 1..maxit loop
      put(file,"Step "); put(file,k,1); put_line(file," :");
      Multitasked_LU_Newton_Step(nbt,s,x,absdx,info,ipvt,wrk,output);
      put(file,"  info : "); put(file,info,1);
      put(file,"  absdx : "); put(file,absdx,3); new_line(file);
      if absdx <= tol
       then fail := false; nbrit := k; exit;
      end if;
    end loop;
  end Multitasked_LU_Newton_Steps;

  procedure Multitasked_LU_Newton_Steps
              ( nbt : in integer32;
                s : in QuadDobl_Speelpenning_Convolutions.Link_to_System;
                x : in QuadDobl_Complex_VecVecs.VecVec;
                maxit : in integer32; nbrit : out integer32;
		tol : in quad_double; absdx : out quad_double; 
		fail : out boolean; info : out integer32;
                ipvt : out Standard_Integer_Vectors.Vector;
                wrk : in QuadDobl_Complex_VecVecs.VecVec;
                output : in boolean := false ) is
  begin
    fail := true; nbrit := maxit;
    for k in 1..maxit loop
      Multitasked_LU_Newton_Step(nbt,s,x,absdx,info,ipvt,wrk,output);
      if absdx <= tol
       then fail := false; nbrit := k; exit;
      end if;
    end loop;
  end Multitasked_LU_Newton_Steps;

  procedure Multitasked_LU_Newton_Steps
              ( file : in file_type; nbt : in integer32;
                s : in QuadDobl_Speelpenning_Convolutions.Link_to_System;
                x : in QuadDobl_Complex_VecVecs.VecVec;
                maxit : in integer32; nbrit : out integer32;
		tol : in quad_double; absdx : out quad_double; 
		fail : out boolean; info : out integer32;
                ipvt : out Standard_Integer_Vectors.Vector;
                wrk : in QuadDobl_Complex_VecVecs.VecVec;
                output : in boolean := false ) is
  begin
    fail := true; nbrit := maxit;
    for k in 1..maxit loop
      put(file,"Step "); put(file,k,1); put_line(file," :");
      Multitasked_LU_Newton_Step(nbt,s,x,absdx,info,ipvt,wrk,output);
      put(file,"  info : "); put(file,info,1);
      put(file,"  absdx : "); put(file,absdx,3); new_line(file);
      if absdx <= tol
       then fail := false; nbrit := k; exit;
      end if;
    end loop;
  end Multitasked_LU_Newton_Steps;

  procedure Multitasked_LU_Newton_Steps
              ( nbt : in integer32;
                s : in PentDobl_Speelpenning_Convolutions.Link_to_System;
                x : in PentDobl_Complex_VecVecs.VecVec;
                maxit : in integer32; nbrit : out integer32;
		tol : in penta_double; absdx : out penta_double; 
		fail : out boolean; info : out integer32;
                ipvt : out Standard_Integer_Vectors.Vector;
                wrk : in PentDobl_Complex_VecVecs.VecVec;
                output : in boolean := false ) is
  begin
    fail := true; nbrit := maxit;
    for k in 1..maxit loop
      Multitasked_LU_Newton_Step(nbt,s,x,absdx,info,ipvt,wrk,output);
      if absdx <= tol
       then fail := false; nbrit := k; exit;
      end if;
    end loop;
  end Multitasked_LU_Newton_Steps;

  procedure Multitasked_LU_Newton_Steps
              ( file : in file_type; nbt : in integer32;
                s : in PentDobl_Speelpenning_Convolutions.Link_to_System;
                x : in PentDobl_Complex_VecVecs.VecVec;
                maxit : in integer32; nbrit : out integer32;
		tol : in penta_double; absdx : out penta_double; 
		fail : out boolean; info : out integer32;
                ipvt : out Standard_Integer_Vectors.Vector;
                wrk : in PentDobl_Complex_VecVecs.VecVec;
                output : in boolean := false ) is
  begin
    fail := true; nbrit := maxit;
    for k in 1..maxit loop
      put(file,"Step "); put(file,k,1); put_line(file," :");
      Multitasked_LU_Newton_Step(nbt,s,x,absdx,info,ipvt,wrk,output);
      put(file,"  info : "); put(file,info,1);
      put(file,"  absdx : "); put(file,absdx,3); new_line(file);
      if absdx <= tol
       then fail := false; nbrit := k; exit;
      end if;
    end loop;
  end Multitasked_LU_Newton_Steps;

  procedure Multitasked_LU_Newton_Steps
              ( nbt : in integer32;
                s : in OctoDobl_Speelpenning_Convolutions.Link_to_System;
                x : in OctoDobl_Complex_VecVecs.VecVec;
                maxit : in integer32; nbrit : out integer32;
		tol : in octo_double; absdx : out octo_double; 
		fail : out boolean; info : out integer32;
                ipvt : out Standard_Integer_Vectors.Vector;
                wrk : in OctoDobl_Complex_VecVecs.VecVec;
                output : in boolean := false ) is
  begin
    fail := true; nbrit := maxit;
    for k in 1..maxit loop
      Multitasked_LU_Newton_Step(nbt,s,x,absdx,info,ipvt,wrk,output);
      if absdx <= tol
       then fail := false; nbrit := k; exit;
      end if;
    end loop;
  end Multitasked_LU_Newton_Steps;

  procedure Multitasked_LU_Newton_Steps
              ( file : in file_type; nbt : in integer32;
                s : in OctoDobl_Speelpenning_Convolutions.Link_to_System;
                x : in OctoDobl_Complex_VecVecs.VecVec;
                maxit : in integer32; nbrit : out integer32;
		tol : in octo_double; absdx : out octo_double; 
		fail : out boolean; info : out integer32;
                ipvt : out Standard_Integer_Vectors.Vector;
                wrk : in OctoDobl_Complex_VecVecs.VecVec;
                output : in boolean := false ) is
  begin
    fail := true; nbrit := maxit;
    for k in 1..maxit loop
      put(file,"Step "); put(file,k,1); put_line(file," :");
      Multitasked_LU_Newton_Step(nbt,s,x,absdx,info,ipvt,wrk,output);
      put(file,"  info : "); put(file,info,1);
      put(file,"  absdx : "); put(file,absdx,3); new_line(file);
      if absdx <= tol
       then fail := false; nbrit := k; exit;
      end if;
    end loop;
  end Multitasked_LU_Newton_Steps;

  procedure Multitasked_LU_Newton_Steps
              ( nbt : in integer32;
                s : in DecaDobl_Speelpenning_Convolutions.Link_to_System;
                x : in DecaDobl_Complex_VecVecs.VecVec;
                maxit : in integer32; nbrit : out integer32;
		tol : in deca_double; absdx : out deca_double; 
		fail : out boolean; info : out integer32;
                ipvt : out Standard_Integer_Vectors.Vector;
                wrk : in DecaDobl_Complex_VecVecs.VecVec;
                output : in boolean := false ) is
  begin
    fail := true; nbrit := maxit;
    for k in 1..maxit loop
      Multitasked_LU_Newton_Step(nbt,s,x,absdx,info,ipvt,wrk,output);
      if absdx <= tol
       then fail := false; nbrit := k; exit;
      end if;
    end loop;
  end Multitasked_LU_Newton_Steps;

  procedure Multitasked_LU_Newton_Steps
              ( file : in file_type; nbt : in integer32;
                s : in DecaDobl_Speelpenning_Convolutions.Link_to_System;
                x : in DecaDobl_Complex_VecVecs.VecVec;
                maxit : in integer32; nbrit : out integer32;
		tol : in deca_double; absdx : out deca_double; 
		fail : out boolean; info : out integer32;
                ipvt : out Standard_Integer_Vectors.Vector;
                wrk : in DecaDobl_Complex_VecVecs.VecVec;
                output : in boolean := false ) is
  begin
    fail := true; nbrit := maxit;
    for k in 1..maxit loop
      put(file,"Step "); put(file,k,1); put_line(file," :");
      Multitasked_LU_Newton_Step(nbt,s,x,absdx,info,ipvt,wrk,output);
      put(file,"  info : "); put(file,info,1);
      put(file,"  absdx : "); put(file,absdx,3); new_line(file);
      if absdx <= tol
       then fail := false; nbrit := k; exit;
      end if;
    end loop;
  end Multitasked_LU_Newton_Steps;

-- SEVERAL NEWTON STEPS WITH LU WITH CONDITION NUMBER ESTIMATE :

  procedure Multitasked_LU_Newton_Steps
              ( nbt : in integer32;
                s : in Standard_Speelpenning_Convolutions.Link_to_System;
                x : in Standard_Complex_VecVecs.VecVec;
                maxit : in integer32; nbrit : out integer32;
                tol : in double_float; absdx : out double_float; 
                fail : out boolean; rcond : out double_float;
                ipvt : out Standard_Integer_Vectors.Vector;
                wrk : in Standard_Complex_VecVecs.VecVec;
                output : in boolean := false ) is
  begin
    fail := true; nbrit := maxit;
    for k in 1..maxit loop
      Multitasked_LU_Newton_Step(nbt,s,x,absdx,rcond,ipvt,wrk,output);
      if absdx <= tol
       then fail := false; nbrit := k; exit;
      end if;
    end loop;
  end Multitasked_LU_Newton_Steps;

  procedure Multitasked_LU_Newton_Steps
              ( file : in file_type; nbt : in integer32;
                s : in Standard_Speelpenning_Convolutions.Link_to_System;
                x : in Standard_Complex_VecVecs.VecVec;
                maxit : in integer32; nbrit : out integer32;
                tol : in double_float; absdx : out double_float; 
                fail : out boolean; rcond : out double_float;
                ipvt : out Standard_Integer_Vectors.Vector;
                wrk : in Standard_Complex_VecVecs.VecVec;
                output : in boolean := false ) is
  begin
    fail := true; nbrit := maxit;
    for k in 1..maxit loop
      put(file,"Step "); put(file,k,1); put_line(file," :");
      Multitasked_LU_Newton_Step(nbt,s,x,absdx,rcond,ipvt,wrk,output);
      put(file,"  rcond :"); put(file,rcond,3);
      put(file,"  absdx :"); put(file,absdx,3); new_line(file);
      if absdx <= tol
       then fail := false; nbrit := k; exit;
      end if;
    end loop;
  end Multitasked_LU_Newton_Steps;

  procedure Multitasked_LU_Newton_Steps
              ( nbt : in integer32;
                s : in DoblDobl_Speelpenning_Convolutions.Link_to_System;
                x : in DoblDobl_Complex_VecVecs.VecVec;
                maxit : in integer32; nbrit : out integer32;
                tol : in double_double; absdx : out double_double; 
                fail : out boolean; rcond : out double_double;
                ipvt : out Standard_Integer_Vectors.Vector;
                wrk : in DoblDobl_Complex_VecVecs.VecVec;
                output : in boolean := false ) is
  begin
    fail := true; nbrit := maxit;
    for k in 1..maxit loop
      Multitasked_LU_Newton_Step(nbt,s,x,absdx,rcond,ipvt,wrk,output);
      if absdx <= tol
       then fail := false; nbrit := k; exit;
      end if;
    end loop;
  end Multitasked_LU_Newton_Steps;

  procedure Multitasked_LU_Newton_Steps
              ( file : in file_type; nbt : in integer32;
                s : in DoblDobl_Speelpenning_Convolutions.Link_to_System;
                x : in DoblDobl_Complex_VecVecs.VecVec;
                maxit : in integer32; nbrit : out integer32;
		tol : in double_double; absdx : out double_double; 
		fail : out boolean; rcond : out double_double;
                ipvt : out Standard_Integer_Vectors.Vector;
                wrk : in DoblDobl_Complex_VecVecs.VecVec;
                output : in boolean := false ) is
  begin
    fail := true; nbrit := maxit;
    for k in 1..maxit loop
      put(file,"Step "); put(file,k,1); put_line(file," :");
      Multitasked_LU_Newton_Step(nbt,s,x,absdx,rcond,ipvt,wrk,output);
      put(file,"  rcond : "); put(file,rcond,3);
      put(file,"  absdx : "); put(file,absdx,3); new_line(file);
      if absdx <= tol
       then fail := false; nbrit := k; exit;
      end if;
    end loop;
  end Multitasked_LU_Newton_Steps;

  procedure Multitasked_LU_Newton_Steps
              ( nbt : in integer32;
                s : in TripDobl_Speelpenning_Convolutions.Link_to_System;
                x : in TripDobl_Complex_VecVecs.VecVec;
                maxit : in integer32; nbrit : out integer32;
                tol : in triple_double; absdx : out triple_double; 
                fail : out boolean; rcond : out triple_double;
                ipvt : out Standard_Integer_Vectors.Vector;
                wrk : in TripDobl_Complex_VecVecs.VecVec;
                output : in boolean := false ) is
  begin
    fail := true; nbrit := maxit;
    for k in 1..maxit loop
      Multitasked_LU_Newton_Step(nbt,s,x,absdx,rcond,ipvt,wrk,output);
      if absdx <= tol
       then fail := false; nbrit := k; exit;
      end if;
    end loop;
  end Multitasked_LU_Newton_Steps;

  procedure Multitasked_LU_Newton_Steps
              ( file : in file_type; nbt : in integer32;
                s : in TripDobl_Speelpenning_Convolutions.Link_to_System;
                x : in TripDobl_Complex_VecVecs.VecVec;
                maxit : in integer32; nbrit : out integer32;
		tol : in triple_double; absdx : out triple_double; 
		fail : out boolean; rcond : out triple_double;
                ipvt : out Standard_Integer_Vectors.Vector;
                wrk : in TripDobl_Complex_VecVecs.VecVec;
                output : in boolean := false ) is
  begin
    fail := true; nbrit := maxit;
    for k in 1..maxit loop
      put(file,"Step "); put(file,k,1); put_line(file," :");
      Multitasked_LU_Newton_Step(nbt,s,x,absdx,rcond,ipvt,wrk,output);
      put(file,"  rcond : "); put(file,rcond,3);
      put(file,"  absdx : "); put(file,absdx,3); new_line(file);
      if absdx <= tol
       then fail := false; nbrit := k; exit;
      end if;
    end loop;
  end Multitasked_LU_Newton_Steps;

  procedure Multitasked_LU_Newton_Steps
              ( nbt : in integer32;
                s : in QuadDobl_Speelpenning_Convolutions.Link_to_System;
                x : in QuadDobl_Complex_VecVecs.VecVec;
                maxit : in integer32; nbrit : out integer32;
		tol : in quad_double; absdx : out quad_double; 
		fail : out boolean; rcond : out quad_double;
                ipvt : out Standard_Integer_Vectors.Vector;
                wrk : in QuadDobl_Complex_VecVecs.VecVec;
                output : in boolean := false ) is
  begin
    fail := true; nbrit := maxit;
    for k in 1..maxit loop
      Multitasked_LU_Newton_Step(nbt,s,x,absdx,rcond,ipvt,wrk,output);
      if absdx <= tol
       then fail := false; nbrit := k; exit;
      end if;
    end loop;
  end Multitasked_LU_Newton_Steps;

  procedure Multitasked_LU_Newton_Steps
              ( file : in file_type; nbt : in integer32;
                s : in QuadDobl_Speelpenning_Convolutions.Link_to_System;
                x : in QuadDobl_Complex_VecVecs.VecVec;
                maxit : in integer32; nbrit : out integer32;
		tol : in quad_double; absdx : out quad_double; 
		fail : out boolean; rcond : out quad_double;
                ipvt : out Standard_Integer_Vectors.Vector;
                wrk : in QuadDobl_Complex_VecVecs.VecVec;
                output : in boolean := false ) is
  begin
    fail := true; nbrit := maxit;
    for k in 1..maxit loop
      put(file,"Step "); put(file,k,1); put_line(file," :");
      Multitasked_LU_Newton_Step(nbt,s,x,absdx,rcond,ipvt,wrk,output);
      put(file,"  rcond : "); put(file,rcond,3);
      put(file,"  absdx : "); put(file,absdx,3); new_line(file);
      if absdx <= tol
       then fail := false; nbrit := k; exit;
      end if;
    end loop;
  end Multitasked_LU_Newton_Steps;

  procedure Multitasked_LU_Newton_Steps
              ( nbt : in integer32;
                s : in PentDobl_Speelpenning_Convolutions.Link_to_System;
                x : in PentDobl_Complex_VecVecs.VecVec;
                maxit : in integer32; nbrit : out integer32;
		tol : in penta_double; absdx : out penta_double; 
		fail : out boolean; rcond : out penta_double;
                ipvt : out Standard_Integer_Vectors.Vector;
                wrk : in PentDobl_Complex_VecVecs.VecVec;
                output : in boolean := false ) is
  begin
    fail := true; nbrit := maxit;
    for k in 1..maxit loop
      Multitasked_LU_Newton_Step(nbt,s,x,absdx,rcond,ipvt,wrk,output);
      if absdx <= tol
       then fail := false; nbrit := k; exit;
      end if;
    end loop;
  end Multitasked_LU_Newton_Steps;

  procedure Multitasked_LU_Newton_Steps
              ( file : in file_type; nbt : in integer32;
                s : in PentDobl_Speelpenning_Convolutions.Link_to_System;
                x : in PentDobl_Complex_VecVecs.VecVec;
                maxit : in integer32; nbrit : out integer32;
		tol : in penta_double; absdx : out penta_double; 
		fail : out boolean; rcond : out penta_double;
                ipvt : out Standard_Integer_Vectors.Vector;
                wrk : in PentDobl_Complex_VecVecs.VecVec;
                output : in boolean := false ) is
  begin
    fail := true; nbrit := maxit;
    for k in 1..maxit loop
      put(file,"Step "); put(file,k,1); put_line(file," :");
      Multitasked_LU_Newton_Step(nbt,s,x,absdx,rcond,ipvt,wrk,output);
      put(file,"  rcond : "); put(file,rcond,3);
      put(file,"  absdx : "); put(file,absdx,3); new_line(file);
      if absdx <= tol
       then fail := false; nbrit := k; exit;
      end if;
    end loop;
  end Multitasked_LU_Newton_Steps;

  procedure Multitasked_LU_Newton_Steps
              ( nbt : in integer32;
                s : in OctoDobl_Speelpenning_Convolutions.Link_to_System;
                x : in OctoDobl_Complex_VecVecs.VecVec;
                maxit : in integer32; nbrit : out integer32;
		tol : in octo_double; absdx : out octo_double; 
		fail : out boolean; rcond : out octo_double;
                ipvt : out Standard_Integer_Vectors.Vector;
                wrk : in OctoDobl_Complex_VecVecs.VecVec;
                output : in boolean := false ) is
  begin
    fail := true; nbrit := maxit;
    for k in 1..maxit loop
      Multitasked_LU_Newton_Step(nbt,s,x,absdx,rcond,ipvt,wrk,output);
      if absdx <= tol
       then fail := false; nbrit := k; exit;
      end if;
    end loop;
  end Multitasked_LU_Newton_Steps;

  procedure Multitasked_LU_Newton_Steps
              ( file : in file_type; nbt : in integer32;
                s : in OctoDobl_Speelpenning_Convolutions.Link_to_System;
                x : in OctoDobl_Complex_VecVecs.VecVec;
                maxit : in integer32; nbrit : out integer32;
		tol : in octo_double; absdx : out octo_double; 
		fail : out boolean; rcond : out octo_double;
                ipvt : out Standard_Integer_Vectors.Vector;
                wrk : in OctoDobl_Complex_VecVecs.VecVec;
                output : in boolean := false ) is
  begin
    fail := true; nbrit := maxit;
    for k in 1..maxit loop
      put(file,"Step "); put(file,k,1); put_line(file," :");
      Multitasked_LU_Newton_Step(nbt,s,x,absdx,rcond,ipvt,wrk,output);
      put(file,"  rcond : "); put(file,rcond,3);
      put(file,"  absdx : "); put(file,absdx,3); new_line(file);
      if absdx <= tol
       then fail := false; nbrit := k; exit;
      end if;
    end loop;
  end Multitasked_LU_Newton_Steps;

  procedure Multitasked_LU_Newton_Steps
              ( nbt : in integer32;
                s : in DecaDobl_Speelpenning_Convolutions.Link_to_System;
                x : in DecaDobl_Complex_VecVecs.VecVec;
                maxit : in integer32; nbrit : out integer32;
		tol : in deca_double; absdx : out deca_double; 
		fail : out boolean; rcond : out deca_double;
                ipvt : out Standard_Integer_Vectors.Vector;
                wrk : in DecaDobl_Complex_VecVecs.VecVec;
                output : in boolean := false ) is
  begin
    fail := true; nbrit := maxit;
    for k in 1..maxit loop
      Multitasked_LU_Newton_Step(nbt,s,x,absdx,rcond,ipvt,wrk,output);
      if absdx <= tol
       then fail := false; nbrit := k; exit;
      end if;
    end loop;
  end Multitasked_LU_Newton_Steps;

  procedure Multitasked_LU_Newton_Steps
              ( file : in file_type; nbt : in integer32;
                s : in DecaDobl_Speelpenning_Convolutions.Link_to_System;
                x : in DecaDobl_Complex_VecVecs.VecVec;
                maxit : in integer32; nbrit : out integer32;
		tol : in deca_double; absdx : out deca_double; 
		fail : out boolean; rcond : out deca_double;
                ipvt : out Standard_Integer_Vectors.Vector;
                wrk : in DecaDobl_Complex_VecVecs.VecVec;
                output : in boolean := false ) is
  begin
    fail := true; nbrit := maxit;
    for k in 1..maxit loop
      put(file,"Step "); put(file,k,1); put_line(file," :");
      Multitasked_LU_Newton_Step(nbt,s,x,absdx,rcond,ipvt,wrk,output);
      put(file,"  rcond : "); put(file,rcond,3);
      put(file,"  absdx : "); put(file,absdx,3); new_line(file);
      if absdx <= tol
       then fail := false; nbrit := k; exit;
      end if;
    end loop;
  end Multitasked_LU_Newton_Steps;

end Multitasked_Newton_Convolutions; 
