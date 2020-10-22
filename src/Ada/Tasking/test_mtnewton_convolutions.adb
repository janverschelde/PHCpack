with duration_io;
with Ada.Calendar;
with Communications_with_User;           use Communications_with_User;
with Time_Stamps;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;
with Double_Double_Numbers;              use Double_Double_Numbers;
with Double_Double_Numbers_io;           use Double_Double_Numbers_io;
with Triple_Double_Numbers;              use Triple_Double_Numbers;
with Triple_Double_Numbers_io;           use Triple_Double_Numbers_io;
with Quad_Double_Numbers;                use Quad_Double_Numbers;
with Quad_Double_Numbers_io;             use Quad_Double_Numbers_io;
with Penta_Double_Numbers;               use Penta_Double_Numbers;
with Penta_Double_Numbers_io;            use Penta_Double_Numbers_io;
with Octo_Double_Numbers;                use Octo_Double_Numbers;
with Octo_Double_Numbers_io;             use Octo_Double_Numbers_io;
with Deca_Double_Numbers;                use Deca_Double_Numbers;
with Deca_Double_Numbers_io;             use Deca_Double_Numbers_io;
with Standard_Complex_Numbers;
with DoblDobl_Complex_Numbers;
with TripDobl_Complex_Numbers;
with QuadDobl_Complex_Numbers;
with PentDobl_Complex_Numbers;
with OctoDobl_Complex_Numbers;
with DecaDobl_Complex_Numbers;
with Standard_Integer_Vectors_io;        use Standard_Integer_Vectors_io;
with QuadDobl_Complex_Vectors_cv;
with Standard_Complex_Solutions;
with Standard_System_and_Solutions_io;
with DoblDobl_Complex_Solutions;
with DoblDobl_System_and_Solutions_io;
with TripDobl_Complex_Solutions;
with TripDobl_System_and_Solutions_io;
with QuadDobl_System_and_Solutions_io;
with PentDobl_Complex_Solutions;
with PentDobl_System_and_Solutions_io;
with OctoDobl_Complex_Solutions;
with OctoDobl_System_and_Solutions_io;
with DecaDobl_Complex_Solutions;
with DecaDobl_System_and_Solutions_io;
with System_Convolution_Circuits;        use System_Convolution_Circuits;
with Homotopy_Convolution_Circuits;      use Homotopy_Convolution_Circuits;
with Random_Convolution_Circuits;        use Random_Convolution_Circuits;
with Standard_Newton_Convolutions;
with DoblDobl_Newton_Convolutions;
with TripDobl_Newton_Convolutions;
with QuadDobl_Newton_Convolutions;
with PentDobl_Newton_Convolutions;
with OctoDobl_Newton_Convolutions;
with DecaDobl_Newton_Convolutions;
with Convergence_Radius_Estimates;
with Multitasked_Series_Linearization;
with Multitasked_Newton_Convolutions;    use Multitasked_Newton_Convolutions;

package body Test_mtNewton_Convolutions is

  procedure Apply_Fabry
              ( c : in Standard_Complex_VecVecs.VecVec;
                verbose : in boolean := true ) is

    z : Standard_Complex_Numbers.Complex_Number;
    r,e : double_float;
    fail : boolean;

  begin
    Convergence_Radius_Estimates.Fabry(c,z,r,e,fail,0,verbose);
  end Apply_Fabry;

  procedure Apply_Fabry
              ( c : in DoblDobl_Complex_VecVecs.VecVec;
                verbose : in boolean := true ) is

    z : DoblDobl_Complex_Numbers.Complex_Number;
    r,e : double_double;
    fail : boolean;

  begin
    Convergence_Radius_Estimates.Fabry(c,z,r,e,fail,0,verbose);
  end Apply_Fabry;

  procedure Apply_Fabry
              ( c : in TripDobl_Complex_VecVecs.VecVec;
                verbose : in boolean := true ) is

    z : TripDobl_Complex_Numbers.Complex_Number;
    r,e : triple_double;
    fail : boolean;

  begin
    Convergence_Radius_Estimates.Fabry(c,z,r,e,fail,0,verbose);
  end Apply_Fabry;

  procedure Apply_Fabry
              ( c : in QuadDobl_Complex_VecVecs.VecVec;
                verbose : in boolean := true ) is

    z : QuadDobl_Complex_Numbers.Complex_Number;
    r,e : quad_double;
    fail : boolean;

  begin
    Convergence_Radius_Estimates.Fabry(c,z,r,e,fail,0,verbose);
  end Apply_Fabry;

  procedure Apply_Fabry
              ( c : in PentDobl_Complex_VecVecs.VecVec;
                verbose : in boolean := true ) is

    z : PentDobl_Complex_Numbers.Complex_Number;
    r,e : penta_double;
    fail : boolean;

  begin
    Convergence_Radius_Estimates.Fabry(c,z,r,e,fail,0,verbose);
  end Apply_Fabry;

  procedure Apply_Fabry
              ( c : in OctoDobl_Complex_VecVecs.VecVec;
                verbose : in boolean := true ) is

    z : OctoDobl_Complex_Numbers.Complex_Number;
    r,e : octo_double;
    fail : boolean;

  begin
    Convergence_Radius_Estimates.Fabry(c,z,r,e,fail,0,verbose);
  end Apply_Fabry;

  procedure Apply_Fabry
              ( c : in DecaDobl_Complex_VecVecs.VecVec;
                verbose : in boolean := true ) is

    z : DecaDobl_Complex_Numbers.Complex_Number;
    r,e : deca_double;
    fail : boolean;

  begin
    Convergence_Radius_Estimates.Fabry(c,z,r,e,fail,0,verbose);
  end Apply_Fabry;

  procedure Standard_Run
              ( nbt,dim,maxit : in integer32;
                s : in Standard_Speelpenning_Convolutions.Link_to_System;
                scf : in Standard_Complex_VecVecs.VecVec;
                serelp,mltelp,speedup,efficiency : in out Duration;
                output,estco : in boolean; verbose : in boolean := true ) is

    wks : Standard_Complex_VecVecs.VecVec(1..nbt)
        := Multitasked_Series_Linearization.Allocate_Work_Space(nbt,dim);
    ipvt : Standard_Integer_Vectors.Vector(1..dim);
    info,nbrit : integer32 := 0;
    fail : boolean;
    tol : constant double_float := 1.0E-12;
    rcond,absdx : double_float;
    seristart,seristop,multstart,multstop : Ada.Calendar.Time;

    use Ada.Calendar; -- for the difference operation on Duration

  begin
    if verbose then
      new_line;
      put("Running with "); put(nbt,1); put_line(" tasks ...");
    end if;
    if nbt = 1
     then seristart := Ada.Calendar.Clock;
     else multstart := Ada.Calendar.Clock;
    end if;
    if verbose then
      if estco then
        Multitasked_LU_Newton_Steps
         (standard_output,nbt,s,scf,maxit,nbrit,tol,absdx,fail,
          rcond,ipvt,wks,output);
      else
        Multitasked_LU_Newton_Steps
         (standard_output,nbt,s,scf,maxit,nbrit,tol,absdx,fail,
          info,ipvt,wks,output);
      end if;
    else
      if estco then
        Multitasked_LU_Newton_Steps
         (nbt,s,scf,maxit,nbrit,tol,absdx,fail,rcond,ipvt,wks);
      else
        Multitasked_LU_Newton_Steps
         (nbt,s,scf,maxit,nbrit,tol,absdx,fail,info,ipvt,wks);
      end if;
    end if;
    if nbt = 1 then
      seristop := Ada.Calendar.Clock;
      serelp := seristop - seristart;
    else
      multstop := Ada.Calendar.Clock;
      mltelp := multstop - multstart;
    end if;
    if verbose then
      put("#steps : "); put(nbrit,1); put("  absdx :"); put(absdx,3);
      if fail
       then put("  failed to reach tolerance");
       else put("  succeeded to reach tolerance");
      end if;
      put(tol,3); new_line;
      Apply_Fabry(scf,verbose);
    end if;
    if nbt = 1 then
      if verbose then
        Time_Stamps.Write_Elapsed_Time(standard_output,seristart,seristop);
      end if;
    else
      if verbose then
        Time_Stamps.Write_Elapsed_Time(standard_output,multstart,multstop);
      end if;
      if serelp + 1.0 /= 1.0 then
        speedup := serelp/mltelp;
        efficiency := speedup/duration(nbt);
        efficiency := duration(100)*efficiency;
        if verbose then
          put("The speedup : "); duration_io.put(speedup,1,3);
          put("  efficiency : "); duration_io.put(efficiency,2,2); new_line;
        end if;
      end if;
    end if;
    Standard_Complex_VecVecs.Clear(wks);
  end Standard_Run;

  procedure DoblDobl_Run
              ( nbt,dim,maxit : in integer32;
                s : in DoblDobl_Speelpenning_Convolutions.Link_to_System;
                scf : in DoblDobl_Complex_VecVecs.VecVec;
                serelp,mltelp,speedup,efficiency : in out Duration;
                output,estco : in boolean; verbose : in boolean := true ) is

    wks : DoblDobl_Complex_VecVecs.VecVec(1..nbt)
        := Multitasked_Series_Linearization.Allocate_Work_Space(nbt,dim);
    ipvt : Standard_Integer_Vectors.Vector(1..dim);
    info,nbrit : integer32 := 0;
    fail : boolean;
    tol : constant double_double := create(1.0E-24);
    rcond,absdx : double_double;
    seristart,seristop,multstart,multstop : Ada.Calendar.Time;

    use Ada.Calendar; -- for the difference operation on Duration

  begin
    if verbose then
      new_line;
      put("Running with "); put(nbt,1); put_line(" tasks ...");
    end if;
    if nbt = 1
     then seristart := Ada.Calendar.Clock;
     else multstart := Ada.Calendar.Clock;
    end if;
    if verbose then
      if estco then
        Multitasked_LU_Newton_Steps
          (standard_output,nbt,s,scf,maxit,nbrit,tol,absdx,fail,
           rcond,ipvt,wks,output);
      else
        Multitasked_LU_Newton_Steps
          (standard_output,nbt,s,scf,maxit,nbrit,tol,absdx,fail,
           info,ipvt,wks,output);
      end if;
    else
      if estco then
        Multitasked_LU_Newton_Steps
          (nbt,s,scf,maxit,nbrit,tol,absdx,fail,rcond,ipvt,wks);
      else
        Multitasked_LU_Newton_Steps
          (nbt,s,scf,maxit,nbrit,tol,absdx,fail,info,ipvt,wks);
      end if;
    end if;
    if nbt = 1 then
      seristop := Ada.Calendar.Clock;
      serelp := seristop - seristart;
    else
      multstop := Ada.Calendar.Clock;
      mltelp := multstop - multstart;
    end if;
    if verbose then
      put("#steps : "); put(nbrit,1); put("  absdx :"); put(absdx,3);
      if fail
       then put("  failed to reach tolerance ");
       else put("  succeeded to reach tolerance ");
      end if;
      put(tol,3); new_line;
      Apply_Fabry(scf,verbose);
    end if;
    if nbt = 1 then
      if verbose then
        Time_Stamps.Write_Elapsed_Time(standard_output,seristart,seristop);
      end if;
    else
      if verbose then
        Time_Stamps.Write_Elapsed_Time(standard_output,multstart,multstop);
      end if;
      if serelp + 1.0 /= 1.0 then
        speedup := serelp/mltelp;
        efficiency := speedup/duration(nbt);
        efficiency := duration(100)*efficiency;
        if verbose then
          put("The speedup : "); duration_io.put(speedup,1,3);
          put("  efficiency : "); duration_io.put(efficiency,2,2); new_line;
        end if;
      end if;
    end if;
    DoblDobl_Complex_VecVecs.Clear(wks);
  end DoblDobl_Run;

  procedure TripDobl_Run
              ( nbt,dim,maxit : in integer32;
                s : in TripDobl_Speelpenning_Convolutions.Link_to_System;
                scf : in TripDobl_Complex_VecVecs.VecVec;
                serelp,mltelp,speedup,efficiency : in out Duration;
                output,estco : in boolean; verbose : in boolean := true ) is

    wks : TripDobl_Complex_VecVecs.VecVec(1..nbt)
        := Multitasked_Series_Linearization.Allocate_Work_Space(nbt,dim);
    ipvt : Standard_Integer_Vectors.Vector(1..dim);
    info,nbrit : integer32 := 0;
    fail : boolean;
    tol : constant triple_double := create(1.0E-36);
    rcond,absdx : triple_double;
    seristart,seristop,multstart,multstop : Ada.Calendar.Time;

    use Ada.Calendar; -- for the difference operation on Duration

  begin
    if verbose then
      new_line;
      put("Running with "); put(nbt,1); put_line(" tasks ...");
    end if;
    if nbt = 1
     then seristart := Ada.Calendar.Clock;
     else multstart := Ada.Calendar.Clock;
    end if;
    if verbose then
      if estco then
        Multitasked_LU_Newton_Steps
          (standard_output,nbt,s,scf,maxit,nbrit,tol,absdx,fail,
           rcond,ipvt,wks,output);
      else
        Multitasked_LU_Newton_Steps
          (standard_output,nbt,s,scf,maxit,nbrit,tol,absdx,fail,
           info,ipvt,wks,output);
      end if;
    else
      if estco then
        Multitasked_LU_Newton_Steps
          (nbt,s,scf,maxit,nbrit,tol,absdx,fail,rcond,ipvt,wks);
      else
        Multitasked_LU_Newton_Steps
          (nbt,s,scf,maxit,nbrit,tol,absdx,fail,info,ipvt,wks);
      end if;
    end if;
    if nbt = 1 then
      seristop := Ada.Calendar.Clock;
      serelp := seristop - seristart;
    else
      multstop := Ada.Calendar.Clock;
      mltelp := multstop - multstart;
    end if;
    if verbose then
      put("#steps : "); put(nbrit,1); put("  absdx :"); put(absdx,3);
      if fail
       then put("  failed to reach tolerance ");
       else put("  succeeded to reach tolerance ");
      end if;
      put(tol,3); new_line;
      Apply_Fabry(scf,verbose);
    end if;
    if nbt = 1 then
      if verbose then
        Time_Stamps.Write_Elapsed_Time(standard_output,seristart,seristop);
      end if;
    else
      if verbose then
        Time_Stamps.Write_Elapsed_Time(standard_output,multstart,multstop);
      end if;
      if serelp + 1.0 /= 1.0 then
        speedup := serelp/mltelp;
        efficiency := speedup/duration(nbt);
        efficiency := duration(100)*efficiency;
        if verbose then
          put("The speedup : "); duration_io.put(speedup,1,3);
          put("  efficiency : "); duration_io.put(efficiency,2,2); new_line;
        end if;
      end if;
    end if;
    TripDobl_Complex_VecVecs.Clear(wks);
  end TripDobl_Run;

  procedure QuadDobl_Run
              ( nbt,dim,maxit : in integer32;
                s : in QuadDobl_Speelpenning_Convolutions.Link_to_System;
                scf : in QuadDobl_Complex_VecVecs.VecVec;
                serelp,mltelp,speedup,efficiency : in out Duration;
                output,estco : in boolean; verbose : in boolean := true ) is

    wks : QuadDobl_Complex_VecVecs.VecVec(1..nbt)
        := Multitasked_Series_Linearization.Allocate_Work_Space(nbt,dim);
    ipvt : Standard_Integer_Vectors.Vector(1..dim);
    info,nbrit : integer32 := 0;
    fail : boolean;
    tol : constant quad_double := create(1.0E-48);
    rcond,absdx : quad_double;
    seristart,seristop,multstart,multstop : Ada.Calendar.Time;

    use Ada.Calendar; -- for the difference operation on Duration

  begin
    if verbose then
      new_line;
      put("Running with "); put(nbt,1); put_line(" tasks ...");
    end if;
    if nbt = 1
     then seristart := Ada.Calendar.Clock;
     else multstart := Ada.Calendar.Clock;
    end if;
    if verbose then
      if estco then
        Multitasked_LU_Newton_Steps
          (standard_output,nbt,s,scf,maxit,nbrit,tol,absdx,fail,
           rcond,ipvt,wks,output);
      else
        Multitasked_LU_Newton_Steps
          (standard_output,nbt,s,scf,maxit,nbrit,tol,absdx,fail,
           info,ipvt,wks,output);
      end if;
    else
      if estco then
        Multitasked_LU_Newton_Steps
          (nbt,s,scf,maxit,nbrit,tol,absdx,fail,rcond,ipvt,wks);
      else
        Multitasked_LU_Newton_Steps
          (nbt,s,scf,maxit,nbrit,tol,absdx,fail,info,ipvt,wks);
      end if;
    end if;
    if nbt = 1 then
      seristop := Ada.Calendar.Clock;
      serelp := seristop - seristart;
    else
      multstop := Ada.Calendar.Clock;
      mltelp := multstop - multstart;
    end if;
    if verbose then
      put("#steps : "); put(nbrit,1); put("  absdx :"); put(absdx,3);
      if fail
       then put("  failed to reach tolerance ");
       else put("  succeeded to reach tolerance ");
      end if;
      put(tol,3); new_line;
      Apply_Fabry(scf,verbose);
    end if;
    if nbt = 1 then
      if verbose then
        Time_Stamps.Write_Elapsed_Time(standard_output,seristart,seristop);
      end if;
    else
      if verbose then
        Time_Stamps.Write_Elapsed_Time(standard_output,multstart,multstop);
      end if;
      if serelp + 1.0 /= 1.0 then
        speedup := serelp/mltelp;
        efficiency := speedup/duration(nbt);
        efficiency := duration(100)*efficiency;
        if verbose then
          put("The speedup : "); duration_io.put(speedup,1,3);
          put("  efficiency : "); duration_io.put(efficiency,2,2); new_line;
        end if;
      end if;
    end if;
    QuadDobl_Complex_VecVecs.Clear(wks);
  end QuadDobl_Run;

  procedure PentDobl_Run
              ( nbt,dim,maxit : in integer32;
                s : in PentDobl_Speelpenning_Convolutions.Link_to_System;
                scf : in PentDobl_Complex_VecVecs.VecVec;
                serelp,mltelp,speedup,efficiency : in out Duration;
                output,estco : in boolean; verbose : in boolean := true ) is

    wks : PentDobl_Complex_VecVecs.VecVec(1..nbt)
        := Multitasked_Series_Linearization.Allocate_Work_Space(nbt,dim);
    ipvt : Standard_Integer_Vectors.Vector(1..dim);
    info,nbrit : integer32 := 0;
    fail : boolean;
    tol : constant penta_double := create(1.0E-60);
    rcond,absdx : penta_double;
    seristart,seristop,multstart,multstop : Ada.Calendar.Time;

    use Ada.Calendar; -- for the difference operation on Duration

  begin
    if verbose then
      new_line;
      put("Running with "); put(nbt,1); put_line(" tasks ...");
    end if;
    if nbt = 1
     then seristart := Ada.Calendar.Clock;
     else multstart := Ada.Calendar.Clock;
    end if;
    if verbose then
      if estco then
        Multitasked_LU_Newton_Steps
          (standard_output,nbt,s,scf,maxit,nbrit,tol,absdx,fail,
           rcond,ipvt,wks,output);
      else
        Multitasked_LU_Newton_Steps
          (standard_output,nbt,s,scf,maxit,nbrit,tol,absdx,fail,
           info,ipvt,wks,output);
      end if;
    else
      if estco then
        Multitasked_LU_Newton_Steps
          (nbt,s,scf,maxit,nbrit,tol,absdx,fail,rcond,ipvt,wks);
      else
        Multitasked_LU_Newton_Steps
          (nbt,s,scf,maxit,nbrit,tol,absdx,fail,info,ipvt,wks);
      end if;
    end if;
    if nbt = 1 then
      seristop := Ada.Calendar.Clock;
      serelp := seristop - seristart;
    else
      multstop := Ada.Calendar.Clock;
      mltelp := multstop - multstart;
    end if;
    if verbose then
      put("#steps : "); put(nbrit,1); put("  absdx :"); put(absdx,3);
      if fail
       then put("  failed to reach tolerance ");
       else put("  succeeded to reach tolerance ");
      end if;
      put(tol,3); new_line;
      Apply_Fabry(scf,verbose);
    end if;
    if nbt = 1 then
      if verbose then
        Time_Stamps.Write_Elapsed_Time(standard_output,seristart,seristop);
      end if;
    else
      if verbose then
        Time_Stamps.Write_Elapsed_Time(standard_output,multstart,multstop);
      end if;
      if serelp + 1.0 /= 1.0 then
        speedup := serelp/mltelp;
        efficiency := speedup/duration(nbt);
        efficiency := duration(100)*efficiency;
        if verbose then
          put("The speedup : "); duration_io.put(speedup,1,3);
          put("  efficiency : "); duration_io.put(efficiency,2,2); new_line;
        end if;
      end if;
    end if;
    PentDobl_Complex_VecVecs.Clear(wks);
  end PentDobl_Run;

  procedure OctoDobl_Run
              ( nbt,dim,maxit : in integer32;
                s : in OctoDobl_Speelpenning_Convolutions.Link_to_System;
                scf : in OctoDobl_Complex_VecVecs.VecVec;
                serelp,mltelp,speedup,efficiency : in out Duration;
                output,estco : in boolean; verbose : in boolean := true ) is

    wks : OctoDobl_Complex_VecVecs.VecVec(1..nbt)
        := Multitasked_Series_Linearization.Allocate_Work_Space(nbt,dim);
    ipvt : Standard_Integer_Vectors.Vector(1..dim);
    info,nbrit : integer32 := 0;
    fail : boolean;
    tol : constant octo_double := create(1.0E-96);
    rcond,absdx : octo_double;
    seristart,seristop,multstart,multstop : Ada.Calendar.Time;

    use Ada.Calendar; -- for the difference operation on Duration

  begin
    if verbose then
      new_line;
      put("Running with "); put(nbt,1); put_line(" tasks ...");
    end if;
    if nbt = 1
     then seristart := Ada.Calendar.Clock;
     else multstart := Ada.Calendar.Clock;
    end if;
    if verbose then
      if estco then
        Multitasked_LU_Newton_Steps
          (standard_output,nbt,s,scf,maxit,nbrit,tol,absdx,fail,
           rcond,ipvt,wks,output);
      else
        Multitasked_LU_Newton_Steps
          (standard_output,nbt,s,scf,maxit,nbrit,tol,absdx,fail,
           info,ipvt,wks,output);
      end if;
    else
      if estco then
        Multitasked_LU_Newton_Steps
          (nbt,s,scf,maxit,nbrit,tol,absdx,fail,rcond,ipvt,wks);
      else
        Multitasked_LU_Newton_Steps
          (nbt,s,scf,maxit,nbrit,tol,absdx,fail,info,ipvt,wks);
      end if;
    end if;
    if nbt = 1 then
      seristop := Ada.Calendar.Clock;
      serelp := seristop - seristart;
    else
      multstop := Ada.Calendar.Clock;
      mltelp := multstop - multstart;
    end if;
    if verbose then
      put("#steps : "); put(nbrit,1); put("  absdx :"); put(absdx,3);
      if fail
       then put("  failed to reach tolerance ");
       else put("  succeeded to reach tolerance ");
      end if;
      put(tol,3); new_line;
      Apply_Fabry(scf,verbose);
    end if;
    if nbt = 1 then
      if verbose then
        Time_Stamps.Write_Elapsed_Time(standard_output,seristart,seristop);
      end if;
    else
      if verbose then
        Time_Stamps.Write_Elapsed_Time(standard_output,multstart,multstop);
      end if;
      if serelp + 1.0 /= 1.0 then
        speedup := serelp/mltelp;
        efficiency := speedup/duration(nbt);
        efficiency := duration(100)*efficiency;
        if verbose then
          put("The speedup : "); duration_io.put(speedup,1,3);
          put("  efficiency : "); duration_io.put(efficiency,2,2); new_line;
        end if;
      end if;
    end if;
    OctoDobl_Complex_VecVecs.Clear(wks);
  end OctoDobl_Run;

  procedure DecaDobl_Run
              ( nbt,dim,maxit : in integer32;
                s : in DecaDobl_Speelpenning_Convolutions.Link_to_System;
                scf : in DecaDobl_Complex_VecVecs.VecVec;
                serelp,mltelp,speedup,efficiency : in out Duration;
                output,estco : in boolean; verbose : in boolean := true ) is

    wks : DecaDobl_Complex_VecVecs.VecVec(1..nbt)
        := Multitasked_Series_Linearization.Allocate_Work_Space(nbt,dim);
    ipvt : Standard_Integer_Vectors.Vector(1..dim);
    info,nbrit : integer32 := 0;
    fail : boolean;
    tol : constant deca_double := create(1.0E-108);
    rcond,absdx : deca_double;
    seristart,seristop,multstart,multstop : Ada.Calendar.Time;

    use Ada.Calendar; -- for the difference operation on Duration

  begin
    if verbose then
      new_line;
      put("Running with "); put(nbt,1); put_line(" tasks ...");
    end if;
    if nbt = 1
     then seristart := Ada.Calendar.Clock;
     else multstart := Ada.Calendar.Clock;
    end if;
    if verbose then
      if estco then
        Multitasked_LU_Newton_Steps
          (standard_output,nbt,s,scf,maxit,nbrit,tol,absdx,fail,
           rcond,ipvt,wks,output);
      else
        Multitasked_LU_Newton_Steps
          (standard_output,nbt,s,scf,maxit,nbrit,tol,absdx,fail,
           info,ipvt,wks,output);
      end if;
    else
      if estco then
        Multitasked_LU_Newton_Steps
          (nbt,s,scf,maxit,nbrit,tol,absdx,fail,rcond,ipvt,wks);
      else
        Multitasked_LU_Newton_Steps
          (nbt,s,scf,maxit,nbrit,tol,absdx,fail,info,ipvt,wks);
      end if;
    end if;
    if nbt = 1 then
      seristop := Ada.Calendar.Clock;
      serelp := seristop - seristart;
    else
      multstop := Ada.Calendar.Clock;
      mltelp := multstop - multstart;
    end if;
    if verbose then
      put("#steps : "); put(nbrit,1); put("  absdx :"); put(absdx,3);
      if fail
       then put("  failed to reach tolerance ");
       else put("  succeeded to reach tolerance ");
      end if;
      put(tol,3); new_line;
      Apply_Fabry(scf,verbose);
    end if;
    if nbt = 1 then
      if verbose then
        Time_Stamps.Write_Elapsed_Time(standard_output,seristart,seristop);
      end if;
    else
      if verbose then
        Time_Stamps.Write_Elapsed_Time(standard_output,multstart,multstop);
      end if;
      if serelp + 1.0 /= 1.0 then
        speedup := serelp/mltelp;
        efficiency := speedup/duration(nbt);
        efficiency := duration(100)*efficiency;
        if verbose then
          put("The speedup : "); duration_io.put(speedup,1,3);
          put("  efficiency : "); duration_io.put(efficiency,2,2); new_line;
        end if;
      end if;
    end if;
    DecaDobl_Complex_VecVecs.Clear(wks);
  end DecaDobl_Run;

  procedure Standard_Run_Loop
              ( p : in Standard_Complex_Poly_Systems.Link_to_Poly_Sys;
	        sol : in Standard_Complex_Vectors.Vector;
	        deg : in integer32 ) is

    use Standard_Speelpenning_Convolutions;

    c : constant Circuits(p'range)
      := Make_Convolution_Circuits(p.all,natural32(deg));
    dim : constant integer32 := sol'last;
    s : constant Link_to_System := Create(c,dim,deg);
    scf : Standard_Complex_VecVecs.VecVec(1..dim);
    maxit,nbt : integer32 := 0;
    ans : character;
    otp,est : boolean;
    seri_elapsed,mult_elapsed,speedup,efficiency : Duration := 0.0;

  begin
    Add_Parameter_to_Constant(s);
    new_line;
    put("Give the maximum number of iterations : "); get(maxit);
    loop
      put("Give the number of tasks (0 to exit) : "); get(nbt);
      exit when (nbt = 0);
      put("Output during multitasking ? (y/n) "); Ask_Yes_or_No(ans);
      otp := (ans = 'y');
      put("Estimate condition number ? (y/n) "); Ask_Yes_or_No(ans);
      est := (ans = 'y');
      scf := Standard_Newton_Convolutions.Series_Coefficients(sol,deg);
      Standard_Run(nbt,dim,maxit,s,scf,seri_elapsed,mult_elapsed,
                   speedup,efficiency,otp,est);
      Standard_Complex_VecVecs.Clear(scf);
    end loop;
  end Standard_Run_Loop;

  procedure DoblDobl_Run_Loop
              ( p : in DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
	        sol : in DoblDobl_Complex_Vectors.Vector;
	        deg : in integer32 ) is

    use DoblDobl_Speelpenning_Convolutions;

    c : constant Circuits(p'range)
      := Make_Convolution_Circuits(p.all,natural32(deg));
    dim : constant integer32 := sol'last;
    s : constant Link_to_System := Create(c,dim,deg);
    scf : DoblDobl_Complex_VecVecs.VecVec(1..dim);
    maxit,nbt : integer32 := 0;
    ans : character;
    otp,est : boolean;
    seri_elapsed,mult_elapsed,speedup,efficiency : Duration := 0.0;

  begin
    Add_Parameter_to_Constant(s);
    new_line;
    put("Give the maximum number of iterations : "); get(maxit);
    loop
      put("Give the number of tasks (0 to exit) : "); get(nbt);
      exit when (nbt = 0);
      put("Output during multitasking ? (y/n) "); Ask_Yes_or_No(ans);
      otp := (ans = 'y');
      put("Estimate condition number ? (y/n) "); Ask_Yes_or_No(ans);
      est := (ans = 'y');
      scf := DoblDobl_Newton_Convolutions.Series_Coefficients(sol,deg);
      DoblDobl_Run(nbt,dim,maxit,s,scf,seri_elapsed,mult_elapsed,
                   speedup,efficiency,otp,est);
      DoblDobl_Complex_VecVecs.Clear(scf);
    end loop;
  end DoblDobl_Run_Loop;

  procedure TripDobl_Run_Loop
              ( p : in TripDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
	        sol : in TripDobl_Complex_Vectors.Vector;
	        deg : in integer32 ) is

    use TripDobl_Speelpenning_Convolutions;

    c : constant Circuits(p'range)
      := Make_Convolution_Circuits(p.all,natural32(deg));
    dim : constant integer32 := sol'last;
    s : constant Link_to_System := Create(c,dim,deg);
    scf : TripDobl_Complex_VecVecs.VecVec(1..dim);
    maxit,nbt : integer32 := 0;
    ans : character;
    otp,est : boolean;
    seri_elapsed,mult_elapsed,speedup,efficiency : Duration := 0.0;

  begin
    Add_Parameter_to_Constant(s);
    new_line;
    put("Give the maximum number of iterations : "); get(maxit);
    loop
      put("Give the number of tasks (0 to exit) : "); get(nbt);
      exit when (nbt = 0);
      put("Output during multitasking ? (y/n) "); Ask_Yes_or_No(ans);
      otp := (ans = 'y');
      put("Estimate condition number ? (y/n) "); Ask_Yes_or_No(ans);
      est := (ans = 'y');
      scf := TripDobl_Newton_Convolutions.Series_Coefficients(sol,deg);
      TripDobl_Run(nbt,dim,maxit,s,scf,seri_elapsed,mult_elapsed,
                   speedup,efficiency,otp,est);
      TripDobl_Complex_VecVecs.Clear(scf);
    end loop;
  end TripDobl_Run_Loop;

  procedure QuadDobl_Run_Loop
              ( p : in QuadDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
	        sol : in QuadDobl_Complex_Vectors.Vector;
	        deg : in integer32 ) is

    use QuadDobl_Speelpenning_Convolutions;

    c : constant Circuits(p'range)
      := Make_Convolution_Circuits(p.all,natural32(deg));
    dim : constant integer32 := sol'last;
    s : constant Link_to_System := Create(c,dim,deg);
    scf : QuadDobl_Complex_VecVecs.VecVec(1..dim);
    maxit,nbt : integer32 := 0;
    ans : character;
    otp,est : boolean;
    seri_elapsed,mult_elapsed,speedup,efficiency : Duration := 0.0;

  begin
    Add_Parameter_to_Constant(s);
    new_line;
    put("Give the maximum number of iterations : "); get(maxit);
    loop
      put("Give the number of tasks (0 to exit) : "); get(nbt);
      exit when (nbt = 0);
      put("Output during multitasking ? (y/n) "); Ask_Yes_or_No(ans);
      otp := (ans = 'y');
      put("Estimate condition number ? (y/n) "); Ask_Yes_or_No(ans);
      est := (ans = 'y');
      scf := QuadDobl_Newton_Convolutions.Series_Coefficients(sol,deg);
      QuadDobl_Run(nbt,dim,maxit,s,scf,seri_elapsed,mult_elapsed,
                   speedup,efficiency,otp,est);
      QuadDobl_Complex_VecVecs.Clear(scf);
    end loop;
  end QuadDobl_Run_Loop;

  procedure PentDobl_Run_Loop
              ( p : in PentDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
	        sol : in PentDobl_Complex_Vectors.Vector;
	        deg : in integer32 ) is

    use PentDobl_Speelpenning_Convolutions;

    c : constant Circuits(p'range)
      := Make_Convolution_Circuits(p.all,natural32(deg));
    dim : constant integer32 := sol'last;
    s : constant Link_to_System := Create(c,dim,deg);
    scf : PentDobl_Complex_VecVecs.VecVec(1..dim);
    maxit,nbt : integer32 := 0;
    ans : character;
    otp,est : boolean;
    seri_elapsed,mult_elapsed,speedup,efficiency : Duration := 0.0;

  begin
    Add_Parameter_to_Constant(s);
    new_line;
    put("Give the maximum number of iterations : "); get(maxit);
    loop
      put("Give the number of tasks (0 to exit) : "); get(nbt);
      exit when (nbt = 0);
      put("Output during multitasking ? (y/n) "); Ask_Yes_or_No(ans);
      otp := (ans = 'y');
      put("Estimate condition number ? (y/n) "); Ask_Yes_or_No(ans);
      est := (ans = 'y');
      scf := PentDobl_Newton_Convolutions.Series_Coefficients(sol,deg);
      PentDobl_Run(nbt,dim,maxit,s,scf,seri_elapsed,mult_elapsed,
                   speedup,efficiency,otp,est);
      PentDobl_Complex_VecVecs.Clear(scf);
    end loop;
  end PentDobl_Run_Loop;

  procedure OctoDobl_Run_Loop
              ( p : in OctoDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
	        sol : in OctoDobl_Complex_Vectors.Vector;
	        deg : in integer32 ) is

    use OctoDobl_Speelpenning_Convolutions;

    c : constant Circuits(p'range)
      := Make_Convolution_Circuits(p.all,natural32(deg));
    dim : constant integer32 := sol'last;
    s : constant Link_to_System := Create(c,dim,deg);
    scf : OctoDobl_Complex_VecVecs.VecVec(1..dim);
    maxit,nbt : integer32 := 0;
    ans : character;
    otp,est : boolean;
    seri_elapsed,mult_elapsed,speedup,efficiency : Duration := 0.0;

  begin
    Add_Parameter_to_Constant(s);
    new_line;
    put("Give the maximum number of iterations : "); get(maxit);
    loop
      put("Give the number of tasks (0 to exit) : "); get(nbt);
      exit when (nbt = 0);
      put("Output during multitasking ? (y/n) "); Ask_Yes_or_No(ans);
      otp := (ans = 'y');
      put("Estimate condition number ? (y/n) "); Ask_Yes_or_No(ans);
      est := (ans = 'y');
      scf := OctoDobl_Newton_Convolutions.Series_Coefficients(sol,deg);
      OctoDobl_Run(nbt,dim,maxit,s,scf,seri_elapsed,mult_elapsed,
                   speedup,efficiency,otp,est);
      OctoDobl_Complex_VecVecs.Clear(scf);
    end loop;
  end OctoDobl_Run_Loop;

  procedure DecaDobl_Run_Loop
              ( p : in DecaDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
	        sol : in DecaDobl_Complex_Vectors.Vector;
	        deg : in integer32 ) is

    use DecaDobl_Speelpenning_Convolutions;

    c : constant Circuits(p'range)
      := Make_Convolution_Circuits(p.all,natural32(deg));
    dim : constant integer32 := sol'last;
    s : constant Link_to_System := Create(c,dim,deg);
    scf : DecaDobl_Complex_VecVecs.VecVec(1..dim);
    maxit,nbt : integer32 := 0;
    ans : character;
    otp,est : boolean;
    seri_elapsed,mult_elapsed,speedup,efficiency : Duration := 0.0;

  begin
    Add_Parameter_to_Constant(s);
    new_line;
    put("Give the maximum number of iterations : "); get(maxit);
    loop
      put("Give the number of tasks (0 to exit) : "); get(nbt);
      exit when (nbt = 0);
      put("Output during multitasking ? (y/n) "); Ask_Yes_or_No(ans);
      otp := (ans = 'y');
      put("Estimate condition number ? (y/n) "); Ask_Yes_or_No(ans);
      est := (ans = 'y');
      scf := DecaDobl_Newton_Convolutions.Series_Coefficients(sol,deg);
      DecaDobl_Run(nbt,dim,maxit,s,scf,seri_elapsed,mult_elapsed,
                   speedup,efficiency,otp,est);
      DecaDobl_Complex_VecVecs.Clear(scf);
    end loop;
  end DecaDobl_Run_Loop;

  procedure Standard_Test is

    lp : Standard_Complex_Poly_Systems.Link_to_Poly_Sys;
    sols : Standard_Complex_Solutions.Solution_List;
    nbr,dim : natural32;
    ls : Standard_Complex_Solutions.Link_to_Solution;
    deg : integer32 := 0;

  begin
    new_line;
    put_line("Reading a polynomial system with solutions ...");
    Standard_System_and_Solutions_io.get(lp,sols);
    nbr := Standard_Complex_Solutions.Length_Of(sols);
    ls := Standard_Complex_Solutions.Head_Of(sols);
    dim := natural32(ls.n);
    new_line;
    put("Read "); put(nbr,1); put(" solutions in dimension ");
    put(dim,1); put_line(".");
    new_line;
    put("Give the degree of the series : "); get(deg);
    Standard_Run_Loop(lp,ls.v,deg);
  end Standard_Test;

  procedure DoblDobl_Test is

    lp : DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
    sols : DoblDobl_Complex_Solutions.Solution_List;
    nbr,dim : natural32;
    ls : DoblDobl_Complex_Solutions.Link_to_Solution;
    deg : integer32 := 0;

  begin
    new_line;
    put_line("Reading a polynomial system with solutions ...");
    DoblDobl_System_and_Solutions_io.get(lp,sols);
    nbr := DoblDobl_Complex_Solutions.Length_Of(sols);
    ls := DoblDobl_Complex_Solutions.Head_Of(sols);
    dim := natural32(ls.n);
    new_line;
    put("Read "); put(nbr,1); put(" solutions in dimension ");
    put(dim,1); put_line(".");
    new_line;
    put("Give the degree of the series : "); get(deg);
    DoblDobl_Run_Loop(lp,ls.v,deg);
  end DoblDobl_Test;

  procedure TripDobl_Test is

    lp : TripDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
    sols : TripDobl_Complex_Solutions.Solution_List;
    nbr,dim : natural32;
    ls : TripDobl_Complex_Solutions.Link_to_Solution;
    deg : integer32 := 0;

  begin
    new_line;
    put_line("Reading a polynomial system with solutions ...");
    TripDobl_System_and_Solutions_io.get(lp,sols);
    nbr := TripDobl_Complex_Solutions.Length_Of(sols);
    ls := TripDobl_Complex_Solutions.Head_Of(sols);
    dim := natural32(ls.n);
    new_line;
    put("Read "); put(nbr,1); put(" solutions in dimension ");
    put(dim,1); put_line(".");
    new_line;
    put("Give the degree of the series : "); get(deg);
    TripDobl_Run_Loop(lp,ls.v,deg);
  end TripDobl_Test;

  procedure QuadDobl_Test is

    lp : QuadDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
    sols : QuadDobl_Complex_Solutions.Solution_List;
    nbr,dim : natural32;
    ls : QuadDobl_Complex_Solutions.Link_to_Solution;
    deg : integer32 := 0;

  begin
    new_line;
    put_line("Reading a polynomial system with solutions ...");
    QuadDobl_System_and_Solutions_io.get(lp,sols);
    nbr := QuadDobl_Complex_Solutions.Length_Of(sols);
    ls := QuadDobl_Complex_Solutions.Head_Of(sols);
    dim := natural32(ls.n);
    new_line;
    put("Read "); put(nbr,1); put(" solutions in dimension ");
    put(dim,1); put_line(".");
    new_line;
    put("Give the degree of the series : "); get(deg);
    QuadDobl_Run_Loop(lp,ls.v,deg);
  end QuadDobl_Test;

  procedure PentDobl_Test is

    lp : PentDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
    sols : PentDobl_Complex_Solutions.Solution_List;
    nbr,dim : natural32;
    ls : PentDobl_Complex_Solutions.Link_to_Solution;
    deg : integer32 := 0;

  begin
    new_line;
    put_line("Reading a polynomial system with solutions ...");
    PentDobl_System_and_Solutions_io.get(lp,sols);
    nbr := PentDobl_Complex_Solutions.Length_Of(sols);
    ls := PentDobl_Complex_Solutions.Head_Of(sols);
    dim := natural32(ls.n);
    new_line;
    put("Read "); put(nbr,1); put(" solutions in dimension ");
    put(dim,1); put_line(".");
    new_line;
    put("Give the degree of the series : "); get(deg);
    PentDobl_Run_Loop(lp,ls.v,deg);
  end PentDobl_Test;

  procedure OctoDobl_Test is

    lp : OctoDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
    sols : OctoDobl_Complex_Solutions.Solution_List;
    nbr,dim : natural32;
    ls : OctoDobl_Complex_Solutions.Link_to_Solution;
    deg : integer32 := 0;

  begin
    new_line;
    put_line("Reading a polynomial system with solutions ...");
    OctoDobl_System_and_Solutions_io.get(lp,sols);
    nbr := OctoDobl_Complex_Solutions.Length_Of(sols);
    ls := OctoDobl_Complex_Solutions.Head_Of(sols);
    dim := natural32(ls.n);
    new_line;
    put("Read "); put(nbr,1); put(" solutions in dimension ");
    put(dim,1); put_line(".");
    new_line;
    put("Give the degree of the series : "); get(deg);
    OctoDobl_Run_Loop(lp,ls.v,deg);
  end OctoDobl_Test;

  procedure DecaDobl_Test is

    lp : DecaDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
    sols : DecaDobl_Complex_Solutions.Solution_List;
    nbr,dim : natural32;
    ls : DecaDobl_Complex_Solutions.Link_to_Solution;
    deg : integer32 := 0;

  begin
    new_line;
    put_line("Reading a polynomial system with solutions ...");
    DecaDobl_System_and_Solutions_io.get(lp,sols);
    nbr := DecaDobl_Complex_Solutions.Length_Of(sols);
    ls := DecaDobl_Complex_Solutions.Head_Of(sols);
    dim := natural32(ls.n);
    new_line;
    put("Read "); put(nbr,1); put(" solutions in dimension ");
    put(dim,1); put_line(".");
    new_line;
    put("Give the degree of the series : "); get(deg);
    DecaDobl_Run_Loop(lp,ls.v,deg);
  end DecaDobl_Test;

  procedure Standard_Random_Test
              ( dim,deg,nbr,pwr : in integer32 ) is

    use Standard_Complex_VecVecs;
    use Standard_Speelpenning_Convolutions;

    s : Link_to_System;
    x : Link_to_VecVec;
    scf : VecVec(1..dim);
    maxit,nbt : integer32 := 0;
    ans : character;
    otp,est : boolean;
    seri_elapsed,mult_elapsed,speedup,efficiency : Duration := 0.0;

  begin
    Standard_Random_Newton_Homotopy(dim,deg,nbr,pwr,s,x);
    new_line;
    put("Give the maximum number of iterations : "); get(maxit);
    loop
      put("Give the number of tasks (0 to exit) : "); get(nbt);
      exit when (nbt = 0);
      put("Output during multitasking ? (y/n) "); Ask_Yes_or_No(ans);
      otp := (ans = 'y');
      put("Estimate condition number ? (y/n) "); Ask_Yes_or_No(ans);
      est := (ans = 'y');
      Standard_Complex_VecVecs.Copy(x.all,scf);
      Standard_Run(nbt,dim,maxit,s,scf,seri_elapsed,mult_elapsed,
                   speedup,efficiency,otp,est);
      Standard_Complex_VecVecs.Clear(scf);
    end loop;
  end Standard_Random_Test;

  procedure DoblDobl_Random_Test
              ( dim,deg,nbr,pwr : in integer32 ) is

    use DoblDobl_Complex_VecVecs;
    use DoblDobl_Speelpenning_Convolutions;

    s : Link_to_System;
    x : Link_to_VecVec;
    scf : VecVec(1..dim);
    maxit,nbt : integer32 := 0;
    ans : character;
    otp,est : boolean;
    seri_elapsed,mult_elapsed,speedup,efficiency : Duration := 0.0;

  begin
    DoblDobl_Random_Newton_Homotopy(dim,deg,nbr,pwr,s,x);
    new_line;
    put("Give the maximum number of iterations : "); get(maxit);
    loop
      put("Give the number of tasks (0 to exit) : "); get(nbt);
      exit when (nbt = 0);
      put("Output during multitasking ? (y/n) "); Ask_Yes_or_No(ans);
      otp := (ans = 'y');
      put("Estimate condition number ? (y/n) "); Ask_Yes_or_No(ans);
      est := (ans = 'y');
      DoblDobl_Complex_VecVecs.Copy(x.all,scf);
      DoblDobl_Run(nbt,dim,maxit,s,scf,seri_elapsed,mult_elapsed,
                   speedup,efficiency,otp,est);
      DoblDobl_Complex_VecVecs.Clear(scf);
    end loop;
  end DoblDobl_Random_Test;

  procedure TripDobl_Random_Test
              ( dim,deg,nbr,pwr : in integer32 ) is

    use TripDobl_Complex_VecVecs;
    use TripDobl_Speelpenning_Convolutions;

    s : Link_to_System;
    x : Link_to_VecVec;
    scf : VecVec(1..dim);
    maxit,nbt : integer32 := 0;
    ans : character;
    otp,est : boolean;
    seri_elapsed,mult_elapsed,speedup,efficiency : Duration := 0.0;

  begin
    TripDobl_Random_Newton_Homotopy(dim,deg,nbr,pwr,s,x);
    new_line;
    put("Give the maximum number of iterations : "); get(maxit);
    loop
      put("Give the number of tasks (0 to exit) : "); get(nbt);
      exit when (nbt = 0);
      put("Output during multitasking ? (y/n) "); Ask_Yes_or_No(ans);
      otp := (ans = 'y');
      put("Estimate condition number ? (y/n) "); Ask_Yes_or_No(ans);
      est := (ans = 'y');
      TripDobl_Complex_VecVecs.Copy(x.all,scf);
      TripDobl_Run(nbt,dim,maxit,s,scf,seri_elapsed,mult_elapsed,
                   speedup,efficiency,otp,est);
      TripDobl_Complex_VecVecs.Clear(scf);
    end loop;
  end TripDobl_Random_Test;

  procedure QuadDobl_Random_Test
              ( dim,deg,nbr,pwr : in integer32 ) is

    use QuadDobl_Complex_VecVecs;
    use QuadDobl_Speelpenning_Convolutions;

    s : Link_to_System;
    x : Link_to_VecVec;
    scf : VecVec(1..dim);
    maxit,nbt : integer32 := 0;
    ans : character;
    otp,est : boolean;
    seri_elapsed,mult_elapsed,speedup,efficiency : Duration := 0.0;

  begin
    QuadDobl_Random_Newton_Homotopy(dim,deg,nbr,pwr,s,x);
    new_line;
    put("Give the maximum number of iterations : "); get(maxit);
    loop
      put("Give the number of tasks (0 to exit) : "); get(nbt);
      exit when (nbt = 0);
      put("Output during multitasking ? (y/n) "); Ask_Yes_or_No(ans);
      otp := (ans = 'y');
      put("Estimate condition number ? (y/n) "); Ask_Yes_or_No(ans);
      est := (ans = 'y');
      QuadDobl_Complex_VecVecs.Copy(x.all,scf);
      QuadDobl_Run(nbt,dim,maxit,s,scf,seri_elapsed,mult_elapsed,
                   speedup,efficiency,otp,est);
      QuadDobl_Complex_VecVecs.Clear(scf);
    end loop;
  end QuadDobl_Random_Test;

  procedure PentDobl_Random_Test
              ( dim,deg,nbr,pwr : in integer32 ) is

    use PentDobl_Complex_VecVecs;
    use PentDobl_Speelpenning_Convolutions;

    s : Link_to_System;
    x : Link_to_VecVec;
    scf : VecVec(1..dim);
    maxit,nbt : integer32 := 0;
    ans : character;
    otp,est : boolean;
    seri_elapsed,mult_elapsed,speedup,efficiency : Duration := 0.0;

  begin
    PentDobl_Random_Newton_Homotopy(dim,deg,nbr,pwr,s,x);
    new_line;
    put("Give the maximum number of iterations : "); get(maxit);
    loop
      put("Give the number of tasks (0 to exit) : "); get(nbt);
      exit when (nbt = 0);
      put("Output during multitasking ? (y/n) "); Ask_Yes_or_No(ans);
      otp := (ans = 'y');
      put("Estimate condition number ? (y/n) "); Ask_Yes_or_No(ans);
      est := (ans = 'y');
      PentDobl_Complex_VecVecs.Copy(x.all,scf);
      PentDobl_Run(nbt,dim,maxit,s,scf,seri_elapsed,mult_elapsed,
                   speedup,efficiency,otp,est);
      PentDobl_Complex_VecVecs.Clear(scf);
    end loop;
  end PentDobl_Random_Test;

  procedure OctoDobl_Random_Test
              ( dim,deg,nbr,pwr : in integer32 ) is

    use OctoDobl_Complex_VecVecs;
    use OctoDobl_Speelpenning_Convolutions;

    s : Link_to_System;
    x : Link_to_VecVec;
    scf : VecVec(1..dim);
    maxit,nbt : integer32 := 0;
    ans : character;
    otp,est : boolean;
    seri_elapsed,mult_elapsed,speedup,efficiency : Duration := 0.0;

  begin
    OctoDobl_Random_Newton_Homotopy(dim,deg,nbr,pwr,s,x);
    new_line;
    put("Give the maximum number of iterations : "); get(maxit);
    loop
      put("Give the number of tasks (0 to exit) : "); get(nbt);
      exit when (nbt = 0);
      put("Output during multitasking ? (y/n) "); Ask_Yes_or_No(ans);
      otp := (ans = 'y');
      put("Estimate condition number ? (y/n) "); Ask_Yes_or_No(ans);
      est := (ans = 'y');
      OctoDobl_Complex_VecVecs.Copy(x.all,scf);
      OctoDobl_Run(nbt,dim,maxit,s,scf,seri_elapsed,mult_elapsed,
                   speedup,efficiency,otp,est);
      OctoDobl_Complex_VecVecs.Clear(scf);
    end loop;
  end OctoDobl_Random_Test;

  procedure DecaDobl_Random_Test
              ( dim,deg,nbr,pwr : in integer32 ) is

    use DecaDobl_Complex_VecVecs;
    use DecaDobl_Speelpenning_Convolutions;

    s : Link_to_System;
    x : Link_to_VecVec;
    scf : VecVec(1..dim);
    maxit,nbt : integer32 := 0;
    ans : character;
    otp,est : boolean;
    seri_elapsed,mult_elapsed,speedup,efficiency : Duration := 0.0;

  begin
    DecaDobl_Random_Newton_Homotopy(dim,deg,nbr,pwr,s,x);
    new_line;
    put("Give the maximum number of iterations : "); get(maxit);
    loop
      put("Give the number of tasks (0 to exit) : "); get(nbt);
      exit when (nbt = 0);
      put("Output during multitasking ? (y/n) "); Ask_Yes_or_No(ans);
      otp := (ans = 'y');
      put("Estimate condition number ? (y/n) "); Ask_Yes_or_No(ans);
      est := (ans = 'y');
      DecaDobl_Complex_VecVecs.Copy(x.all,scf);
      DecaDobl_Run(nbt,dim,maxit,s,scf,seri_elapsed,mult_elapsed,
                   speedup,efficiency,otp,est);
      DecaDobl_Complex_VecVecs.Clear(scf);
    end loop;
  end DecaDobl_Random_Test;

  procedure Standard_Benchmark
              ( file : in file_type; nbruns,inc,maxit : in integer32;
                nbtseq : in Standard_Integer_Vectors.Link_to_Vector;
                s : in Standard_Speelpenning_Convolutions.Link_to_System;
                x : in Standard_Complex_VecVecs.Link_to_VecVec;
                verbose : in boolean := false ) is

    scf : Standard_Complex_VecVecs.VecVec(1..s.dim);
    nbt : integer32 := 2;
    seri_elapsed,mult_elapsed,speedup,efficiency : Duration := 0.0;

  begin
    put_line(file,"double precision");
    Standard_Complex_VecVecs.Copy(x.all,scf);
    Standard_Run(1,s.dim,maxit,s,scf,seri_elapsed,mult_elapsed,
                 speedup,efficiency,false,false,verbose);
    Standard_Complex_VecVecs.Clear(scf);
    put(file,"  1 : ");
    duration_io.put(file,seri_elapsed,1,3); new_line(file); flush(file);
    if nbruns /= 0 then
      for k in 1..nbruns loop
        Standard_Complex_VecVecs.Copy(x.all,scf);
        Standard_Run(nbt,s.dim,maxit,s,scf,seri_elapsed,mult_elapsed,
                     speedup,efficiency,false,false,verbose);
        Standard_Complex_VecVecs.Clear(scf);
        put(file,nbt,3);
        put(file," : "); duration_io.put(file,mult_elapsed,1,3);
        put(file," : "); duration_io.put(file,speedup,1,3);
        put(file," : "); duration_io.put(file,efficiency,2,2);
        new_line(file); flush(file);
        nbt := nbt + inc;
      end loop;
   else
      for k in nbtseq'range loop
        nbt := nbtseq(k);
        Standard_Complex_VecVecs.Copy(x.all,scf);
        Standard_Run(nbt,s.dim,maxit,s,scf,seri_elapsed,mult_elapsed,
                     speedup,efficiency,false,false,verbose);
        Standard_Complex_VecVecs.Clear(scf);
        put(file,nbt,3);
        put(file," : "); duration_io.put(file,mult_elapsed,1,3);
        put(file," : "); duration_io.put(file,speedup,1,3);
        put(file," : "); duration_io.put(file,efficiency,2,2);
        new_line(file); flush(file);
     end loop;
    end if;
  end Standard_Benchmark;

  procedure DoblDobl_Benchmark
              ( file : in file_type; nbruns,inc,maxit : in integer32;
                nbtseq : in Standard_Integer_Vectors.Link_to_Vector;
                s : in DoblDobl_Speelpenning_Convolutions.Link_to_System;
                x : in DoblDobl_Complex_VecVecs.Link_to_VecVec;
                verbose : in boolean := false ) is

    scf : DoblDobl_Complex_VecVecs.VecVec(1..s.dim);
    nbt : integer32 := 2;
    seri_elapsed,mult_elapsed,speedup,efficiency : Duration := 0.0;

  begin
    put_line(file,"double double precision");
    DoblDobl_Complex_VecVecs.Copy(x.all,scf);
    DoblDobl_Run(1,s.dim,maxit,s,scf,seri_elapsed,mult_elapsed,
                 speedup,efficiency,false,false,verbose);
    DoblDobl_Complex_VecVecs.Clear(scf);
    put(file,"  1 : ");
    duration_io.put(file,seri_elapsed,1,3); new_line(file); flush(file);
    if nbruns /= 0 then
      for k in 1..nbruns loop
        DoblDobl_Complex_VecVecs.Copy(x.all,scf);
        DoblDobl_Run(nbt,s.dim,maxit,s,scf,seri_elapsed,mult_elapsed,
                     speedup,efficiency,false,false,verbose);
        DoblDobl_Complex_VecVecs.Clear(scf);
        put(file,nbt,3);
        put(file," : "); duration_io.put(file,mult_elapsed,1,3);
        put(file," : "); duration_io.put(file,speedup,1,3);
        put(file," : "); duration_io.put(file,efficiency,2,2);
        new_line(file); flush(file);
        nbt := nbt + inc;
      end loop;
    else
      for k in nbtseq'range loop
        nbt := nbtseq(k);
        DoblDobl_Complex_VecVecs.Copy(x.all,scf);
        DoblDobl_Run(nbt,s.dim,maxit,s,scf,seri_elapsed,mult_elapsed,
                     speedup,efficiency,false,false,verbose);
        DoblDobl_Complex_VecVecs.Clear(scf);
        put(file,nbt,3);
        put(file," : "); duration_io.put(file,mult_elapsed,1,3);
        put(file," : "); duration_io.put(file,speedup,1,3);
        put(file," : "); duration_io.put(file,efficiency,2,2);
        new_line(file); flush(file);
      end loop;
    end if;
  end DoblDobl_Benchmark;

  procedure TripDobl_Benchmark
              ( file : in file_type; nbruns,inc,maxit : in integer32;
                nbtseq : in Standard_Integer_Vectors.Link_to_Vector;
                s : in TripDobl_Speelpenning_Convolutions.Link_to_System;
                x : in TripDobl_Complex_VecVecs.Link_to_VecVec;
                verbose : in boolean := false ) is

    scf : TripDobl_Complex_VecVecs.VecVec(1..s.dim);
    nbt : integer32 := 2;
    seri_elapsed,mult_elapsed,speedup,efficiency : Duration := 0.0;

  begin
    put_line(file,"triple double precision");
    TripDobl_Complex_VecVecs.Copy(x.all,scf);
    TripDobl_Run(1,s.dim,maxit,s,scf,seri_elapsed,mult_elapsed,
                 speedup,efficiency,false,false,verbose);
    TripDobl_Complex_VecVecs.Clear(scf);
    put(file,"  1 : ");
    duration_io.put(file,seri_elapsed,1,3); new_line(file); flush(file);
    if nbruns /= 0 then
      for k in 1..nbruns loop
        TripDobl_Complex_VecVecs.Copy(x.all,scf);
        TripDobl_Run(nbt,s.dim,maxit,s,scf,seri_elapsed,mult_elapsed,
                     speedup,efficiency,false,false,verbose);
        TripDobl_Complex_VecVecs.Clear(scf);
        put(file,nbt,3);
        put(file," : "); duration_io.put(file,mult_elapsed,1,3);
        put(file," : "); duration_io.put(file,speedup,1,3);
        put(file," : "); duration_io.put(file,efficiency,2,2);
        new_line(file); flush(file);
        nbt := nbt + inc;
      end loop;
    else
      for k in nbtseq'range loop
        nbt := nbtseq(k);
        TripDobl_Complex_VecVecs.Copy(x.all,scf);
        TripDobl_Run(nbt,s.dim,maxit,s,scf,seri_elapsed,mult_elapsed,
                     speedup,efficiency,false,false,verbose);
        TripDobl_Complex_VecVecs.Clear(scf);
        put(file,nbt,3);
        put(file," : "); duration_io.put(file,mult_elapsed,1,3);
        put(file," : "); duration_io.put(file,speedup,1,3);
        put(file," : "); duration_io.put(file,efficiency,2,2);
        new_line(file); flush(file);
      end loop;
    end if;
  end TripDobl_Benchmark;

  procedure QuadDobl_Benchmark
              ( file : in file_type; nbruns,inc,maxit : in integer32;
                nbtseq : in Standard_Integer_Vectors.Link_to_Vector;
                s : in QuadDobl_Speelpenning_Convolutions.Link_to_System;
                x : in QuadDobl_Complex_VecVecs.Link_to_VecVec;
                verbose : in boolean := false ) is

    scf : QuadDobl_Complex_VecVecs.VecVec(1..s.dim);
    nbt : integer32 := 2;
    seri_elapsed,mult_elapsed,speedup,efficiency : Duration := 0.0;

  begin
    put_line(file,"quad double precision");
    QuadDobl_Complex_VecVecs.Copy(x.all,scf);
    QuadDobl_Run(1,s.dim,maxit,s,scf,seri_elapsed,mult_elapsed,
                 speedup,efficiency,false,false,verbose);
    QuadDobl_Complex_VecVecs.Clear(scf);
    put(file,"  1 : ");
    duration_io.put(file,seri_elapsed,1,3); new_line(file); flush(file);
    if nbruns /= 0 then
      for k in 1..nbruns loop
        QuadDobl_Complex_VecVecs.Copy(x.all,scf);
        QuadDobl_Run(nbt,s.dim,maxit,s,scf,seri_elapsed,mult_elapsed,
                     speedup,efficiency,false,false,verbose);
        QuadDobl_Complex_VecVecs.Clear(scf);
        put(file,nbt,3);
        put(file," : "); duration_io.put(file,mult_elapsed,1,3);
        put(file," : "); duration_io.put(file,speedup,1,3);
        put(file," : "); duration_io.put(file,efficiency,2,2);
        new_line(file); flush(file);
        nbt := nbt + inc;
      end loop;
    else
      for k in nbtseq'range loop
        nbt := nbtseq(k);
        QuadDobl_Complex_VecVecs.Copy(x.all,scf);
        QuadDobl_Run(nbt,s.dim,maxit,s,scf,seri_elapsed,mult_elapsed,
                     speedup,efficiency,false,false,verbose);
        QuadDobl_Complex_VecVecs.Clear(scf);
        put(file,nbt,3);
        put(file," : "); duration_io.put(file,mult_elapsed,1,3);
        put(file," : "); duration_io.put(file,speedup,1,3);
        put(file," : "); duration_io.put(file,efficiency,2,2);
        new_line(file); flush(file);
      end loop;
    end if;
  end QuadDobl_Benchmark;

  function Prompt_for_Sequence
             ( max : in integer32 )
             return Standard_Integer_Vectors.Link_to_Vector is

  -- DESCRIPTION :
  --   Prompts the user for a sequence of numbers
  --   and the sequence is returned in a vector.
  --   The length of the sequence may not exceed max.

    res : Standard_Integer_Vectors.Link_to_Vector;
    seq : Standard_Integer_Vectors.Vector(1..max);
    cnt,nbr : integer32 := 0;
    ans : character;

  begin
    loop
      put_line("Reading a sequence of numbers ...");
      for k in 1..max loop
        put("-> give a number (0 to exit) : "); get(nbr);
        exit when (nbr = 0);
        cnt := cnt + 1;
        seq(cnt) := nbr;
      end loop;
      put("The sequence : "); put(seq(1..cnt)); new_line;
      put("Okay to exit ? (y/n) "); Ask_Yes_or_No(ans);
      exit when (ans = 'y');
      put_line("Entering the sequence again ...");
      cnt := 0;
    end loop;
    res := new Standard_Integer_Vectors.Vector'(seq(1..cnt));
    return res;
  end Prompt_For_Sequence;

  procedure Benchmark ( dim,deg,nbr,pwr : in integer32 ) is

    qds : QuadDobl_Speelpenning_Convolutions.Link_to_System;
    tds : TripDobl_Speelpenning_Convolutions.Link_to_System;
    dds : DoblDobl_Speelpenning_Convolutions.Link_to_System;
    d_s : Standard_Speelpenning_Convolutions.Link_to_System;
    qdx : QuadDobl_Complex_VecVecs.Link_to_VecVec;
    tdx : TripDobl_Complex_VecVecs.Link_to_VecVec;
    ddx : DoblDobl_Complex_VecVecs.Link_to_VecVec;
    d_x : Standard_Complex_VecVecs.Link_to_VecVec;
    file : file_type;
    nbruns,inc,maxit : integer32 := 0;
    nbtseq : Standard_Integer_Vectors.Link_to_Vector;
 
  begin
    new_line;
    put("Give the number of multitasked runs (0 for sequence) : ");
    get(nbruns);
    if nbruns /= 0
     then put("Give the increment on the tasks : "); get(inc);
     else nbtseq := Prompt_for_Sequence(44);
    end if;
    new_line;
    put("Give the maximum number of iterations : "); get(maxit);
    skip_line;
    new_line;
    put_line("Reading the name of the output file ...");
    Read_Name_and_Create_File(file);
    new_line;
    put_line("See the output file for results ...");
    new_line;
    QuadDobl_Random_Newton_Homotopy(dim,deg,nbr,pwr,qds,qdx);
    tds := System_Convolution_Circuits.to_triple_double(qds);
    tdx := QuadDobl_Complex_Vectors_cv.to_triple_double(qdx);
    dds := System_Convolution_Circuits.to_double_double(qds);
    ddx := QuadDobl_Complex_Vectors_cv.to_double_double(qdx);
    d_s := System_Convolution_Circuits.to_double(qds);
    d_x := QuadDobl_Complex_Vectors_cv.to_double(qdx);
    put(file,"dimension : "); put(file,dim,1);
    put(file,"  degree : "); put(file,deg,1);
    put(file,"  largest power : "); put(file,pwr,1); new_line(file);
    put(file,"maximum number of iterations : ");
    put(file,maxit,1); new_line(file);
    Standard_Benchmark(file,nbruns,inc,maxit,nbtseq,d_s,d_x);
    DoblDobl_Benchmark(file,nbruns,inc,maxit,nbtseq,dds,ddx);
    TripDobl_Benchmark(file,nbruns,inc,maxit,nbtseq,tds,tdx);
    QuadDobl_Benchmark(file,nbruns,inc,maxit,nbtseq,qds,qdx);
  end Benchmark;

  procedure Benchmark
              ( p : in QuadDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
                sols : in QuadDobl_Complex_Solutions.Solution_List;
                dim,deg : in integer32 ) is

    file : file_type;
    nbruns,inc,maxit : integer32 := 0;

    use QuadDobl_Speelpenning_Convolutions;

    c : constant QuadDobl_Speelpenning_Convolutions.Circuits(p'range)
      := Make_Convolution_Circuits(p.all,natural32(deg));
    qds : constant QuadDobl_Speelpenning_Convolutions.Link_to_System
        := Create(c,dim,deg);
    dds : DoblDobl_Speelpenning_Convolutions.Link_to_System;
    tds : TripDobl_Speelpenning_Convolutions.Link_to_System;
    d_s : Standard_Speelpenning_Convolutions.Link_to_System;
    scf : QuadDobl_Complex_VecVecs.VecVec(1..dim);
    ls : constant QuadDobl_Complex_Solutions.Link_to_Solution
       := QuadDobl_Complex_Solutions.Head_Of(sols);
    qdx : QuadDobl_Complex_VecVecs.Link_to_VecVec;
    tdx : TripDobl_Complex_VecVecs.Link_to_VecVec;
    ddx : DoblDobl_Complex_VecVecs.Link_to_VecVec;
    d_x : Standard_Complex_VecVecs.Link_to_VecVec;
    nbtseq : Standard_Integer_Vectors.Link_to_Vector;

  begin
    Add_Parameter_to_Constant(qds); -- make Newton homotopy
    dds := System_Convolution_Circuits.to_double_double(qds);
    tds := System_Convolution_Circuits.to_triple_double(qds);
    d_s := System_Convolution_Circuits.to_double(qds);
    scf := QuadDobl_Newton_Convolutions.Series_Coefficients(ls.v,deg);
    qdx := new QuadDobl_Complex_VecVecs.VecVec'(scf);
    ddx := QuadDobl_Complex_Vectors_cv.to_double_double(qdx);
    tdx := QuadDobl_Complex_Vectors_cv.to_triple_double(qdx);
    d_x := QuadDobl_Complex_Vectors_cv.to_double(qdx);
    new_line;
    put("Give the number of multitasked runs (0 for sequence) : ");
    get(nbruns);
    if nbruns /= 0
     then put("Give the increment on the tasks : "); get(inc);
     else nbtseq := Prompt_for_Sequence(44);
    end if;
    new_line;
    put("Give the maximum number of iterations : "); get(maxit);
    skip_line;
    new_line;
    put_line("Reading the name of the output file ...");
    Read_Name_and_Create_File(file);
    new_line;
    put_line("See the output file for results ...");
    new_line;
    put(file,"dimension : "); put(file,dim,1);
    put(file,"  degree : "); put(file,deg,1); new_line(file);
    put(file,"maximum number of iterations : ");
    put(file,maxit,1); new_line(file);
    Standard_Benchmark(file,nbruns,inc,maxit,nbtseq,d_s,d_x);
    DoblDobl_Benchmark(file,nbruns,inc,maxit,nbtseq,dds,ddx);
    TripDobl_Benchmark(file,nbruns,inc,maxit,nbtseq,tds,tdx);
    QuadDobl_Benchmark(file,nbruns,inc,maxit,nbtseq,qds,qdx);
  end Benchmark;

  procedure Prompt_for_Dimensions
              ( dim,deg,nbr,pwr : in out integer32 ) is
  begin
    new_line;
    put("Give the dimension : "); get(dim);
    put("Give the degree of the power series : "); get(deg);
    put("Give the number of terms in each circuit : "); get(nbr);
    put("Give the largest power of the variables : "); get(pwr);
  end Prompt_for_Dimensions;           

  procedure Main is

    prc,ans : character;
    dim,deg,nbr,pwr : integer32 := 0;
    lp : QuadDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
    sols : QuadDobl_Complex_Solutions.Solution_List;
    ls : QuadDobl_Complex_Solutions.Link_to_Solution;

  begin
    new_line;
    put_line("Testing Newton's method on power series ...");
    new_line;
    put("Benchmark for all precisions ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans = 'y' then
      new_line;
      put("Generate a random problem ? (y/n) "); Ask_Yes_or_No(ans);
      if ans = 'y' then
        Prompt_for_Dimensions(dim,deg,nbr,pwr);
        Benchmark(dim,deg,nbr,pwr);
      else
        new_line;
        put_line("Reading a polynomial system with solutions ...");
        QuadDobl_System_and_Solutions_io.get(lp,sols);
        nbr := integer32(QuadDobl_Complex_Solutions.Length_Of(sols));
        ls := QuadDobl_Complex_Solutions.Head_Of(sols);
        dim := ls.n;
        new_line;
        put("Read "); put(nbr,1); put(" solutions in dimension ");
        put(dim,1); put_line(".");
        new_line;
        put("Give the degree of the power series : "); get(deg);
        Benchmark(lp,sols,dim,deg);
      end if;
    else
      new_line;
      put_line("MENU for the working precision :");
      put_line("  1. double precision");
      put_line("  2. double double precision");
      put_line("  3. triple double precision");
      put_line("  4. quad double precision");
      put_line("  5. penta double precision");
      put_line("  6. octo double precision");
      put_line("  7. deca double precision");
      put("Type 1, 2, 3, 4, 5, 6, or 7 to select the precision : ");
      Ask_Alternative(prc,"1234567");
      new_line;
      put("Generate a random problem ? (y/n) "); Ask_Yes_or_No(ans);
      if ans ='y' then
        Prompt_for_Dimensions(dim,deg,nbr,pwr);
        case prc is
          when '1' => Standard_Random_Test(dim,deg,nbr,pwr);
          when '2' => DoblDobl_Random_Test(dim,deg,nbr,pwr);
          when '3' => TripDobl_Random_Test(dim,deg,nbr,pwr);
          when '4' => QuadDobl_Random_Test(dim,deg,nbr,pwr);
          when '5' => PentDobl_Random_Test(dim,deg,nbr,pwr);
          when '6' => OctoDobl_Random_Test(dim,deg,nbr,pwr);
          when '7' => DecaDobl_Random_Test(dim,deg,nbr,pwr);
          when others => null;
        end case;
      else
        case prc is
          when '1' => Standard_Test;
          when '2' => DoblDobl_Test;
          when '3' => TripDobl_Test;
          when '4' => QuadDobl_Test;
          when '5' => PentDobl_Test;
          when '6' => OctoDobl_Test;
          when '7' => DecaDobl_Test;
          when others => null;
        end case;
      end if;
    end if;
  end Main;

end Test_mtNewton_Convolutions;
