with duration_io;
with Ada.Calendar;
with Communications_with_User;           use Communications_with_User;
with Time_Stamps;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
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
with Standard_Integer_Vectors_io;        use Standard_Integer_Vectors_io;
with Standard_Floating_Vectors;
with Standard_Floating_VecVecs;
with Standard_Floating_VecVecVecs;
with Standard_Complex_VecMats;
with DoblDobl_Complex_VecMats;
with TripDobl_Complex_VecMats;
with QuadDobl_Complex_VecMats;
with PentDobl_Complex_VecMats;
with OctoDobl_Complex_VecMats;
with DecaDobl_Complex_VecMats;
with DecaDobl_Complex_Vectors_cv;
with Standard_Complex_Series_Vectors;
with Standard_Random_Series_Vectors;
with DoblDobl_Complex_Series_Vectors;
with DoblDobl_Random_Series_Vectors;
with TripDobl_Complex_Series_Vectors;
with TripDobl_Random_Series_Vectors;
with QuadDobl_Complex_Series_Vectors;
with QuadDobl_Random_Series_Vectors;
with PentDobl_Complex_Series_Vectors;
with PentDobl_Random_Series_Vectors;
with OctoDobl_Complex_Series_Vectors;
with OctoDobl_Random_Series_Vectors;
with DecaDobl_Complex_Series_Vectors;
with DecaDobl_Random_Series_Vectors;
with Series_Coefficient_Vectors;         use Series_Coefficient_Vectors;
with Standard_Complex_Polynomials;
with Standard_Complex_Poly_Systems;
with Standard_Complex_Poly_Systems_io;   use Standard_Complex_Poly_Systems_io;
with DoblDobl_Complex_Polynomials;
with DoblDobl_Complex_Poly_Systems;
with DoblDobl_Complex_Poly_Systems_io;   use DoblDobl_Complex_Poly_Systems_io;
with TripDobl_Complex_Polynomials;
with TripDobl_Complex_Poly_Systems;
with TripDobl_Complex_Poly_Systems_io;   use TripDobl_Complex_Poly_Systems_io;
with QuadDobl_Complex_Polynomials;
with QuadDobl_Complex_Poly_Systems;
with QuadDobl_Complex_Poly_Systems_io;   use QuadDobl_Complex_Poly_Systems_io;
with PentDobl_Complex_Polynomials;
with PentDobl_Complex_Poly_Systems;
with PentDobl_Complex_Poly_Systems_io;   use PentDobl_Complex_Poly_Systems_io;
with OctoDobl_Complex_Polynomials;
with OctoDobl_Complex_Poly_Systems;
with OctoDobl_Complex_Poly_Systems_io;   use OctoDobl_Complex_Poly_Systems_io;
with DecaDobl_Complex_Polynomials;
with DecaDobl_Complex_Poly_Systems;
with DecaDobl_Complex_Poly_Systems_io;   use DecaDobl_Complex_Poly_Systems_io;
with Standard_Vector_Splitters;
with DoblDobl_Vector_Splitters;
with QuadDobl_Vector_Splitters;
with Standard_Convolution_Splitters;
with DoblDobl_Convolution_Splitters;
with QuadDobl_Convolution_Splitters;
with System_Convolution_Circuits;        use System_Convolution_Circuits;
with Random_Convolution_Circuits;        use Random_Convolution_Circuits;
with Evaluation_Differentiation_Errors;  use Evaluation_Differentiation_Errors;
with Multitasked_AlgoDiff_Convolutions;  use Multitasked_AlgoDiff_Convolutions;

package body Test_mtAlgoDiff_Convolutions is

  procedure Standard_Test
              ( c : in Standard_Speelpenning_Convolutions.Circuits;
                dim,deg : in integer32; nbt : in out integer32;
                output : in boolean := true ) is

    use Standard_Speelpenning_Convolutions;

    x : constant Standard_Complex_Series_Vectors.Vector(1..dim)
      := Standard_Random_Series_Vectors.Random_Series_Vector(1,dim,deg);
    xcff : constant Standard_Complex_VecVecs.VecVec(1..dim)
         := Standard_Series_Coefficients(x);
    mxe : constant Standard_Integer_Vectors.Vector(1..dim)
        := Exponent_Maxima(c,dim);
    pwt : constant Link_to_VecVecVec := Allocate(mxe,deg);
    yd : constant Standard_Complex_VecVecs.VecVec(1..dim+1)
       := Allocate_Coefficients(dim+1,deg);
    vy1 : constant Standard_Complex_VecVecs.VecVec(0..deg)
        := Linearized_Allocation(dim,deg);
    vy2 : constant Standard_Complex_VecVecs.VecVec(0..deg)
        := Linearized_Allocation(dim,deg);
    vm1 : constant Standard_Complex_VecMats.VecMat(0..deg)
        := Allocate_Coefficients(dim,dim,deg);
    vm2 : constant Standard_Complex_VecMats.VecMat(0..deg)
        := Allocate_Coefficients(dim,dim,deg);
    err : double_float;
    multstart,multstop,seristart,seristop : Ada.Calendar.Time;
    seri_elapsed,mult_elapsed,speedup,efficiency : Duration;
    ans : character;

    use Ada.Calendar; -- for the difference operator on Duration

  begin
    put_line("Running on one task ...");
    seristart := Ada.Calendar.Clock;
    Compute(pwt,mxe,xcff);
    EvalDiff(c,xcff,pwt,yd,vy1,vm1);
    seristop := Ada.Calendar.Clock;
    seri_elapsed := seristop - seristart;
    put_line("-> Elapsed time without multitasking : ");
    Time_Stamps.Write_Elapsed_Time(standard_output,seristart,seristop);
    loop
      new_line;
      put("Running with "); put(nbt,1); put_line(" tasks ...");
      multstart := Ada.Calendar.Clock;
      Standard_Multitasked_EvalDiff(nbt,c,xcff,mxe,pwt,vy2,vm2,output);
      multstop := Ada.Calendar.Clock;
      mult_elapsed := multstop - multstart;
      err := Difference(vy1,vy2);
      put("  the error of evaluation : "); put(err,3); new_line;
      err := Difference(vm1,vm2);
      put("  the error of differentiation : "); put(err,3); new_line;
      put_line("-> Elapsed time on multitasking : ");
      Time_Stamps.Write_Elapsed_Time(standard_output,multstart,multstop);
      if seri_elapsed + 1.0 /= 1.0 then
        speedup := seri_elapsed/mult_elapsed;
        efficiency := speedup/duration(nbt);
        efficiency := duration(100)*efficiency;
        put("The speedup : "); duration_io.put(speedup,1,3);
        put("  efficiency : "); duration_io.put(efficiency,2,2); new_line;
      end if;
      put("Continue with different number of tasks ? (y/n) ");
      Ask_Yes_or_No(ans);
      exit when (ans /= 'y');
      put("Give the number of tasks : "); get(nbt);
    end loop;
  end Standard_Test;

  procedure Standard_Coefficient_Test
              ( c : in Standard_Coefficient_Convolutions.Circuits;
                dim,deg : in integer32; nbt : in out integer32;
                static : in boolean; output : in boolean := true ) is

    use Standard_Floating_VecVecVecs;
    use Standard_Coefficient_Convolutions;

    s : constant Standard_Coefficient_Convolutions.Link_to_System
      := Standard_Coefficient_Convolutions.Create(c,dim,deg);
    x : constant Standard_Complex_Series_Vectors.Vector(1..dim)
      := Standard_Random_Series_Vectors.Random_Series_Vector(1,dim,deg);
    xcff : constant Standard_Complex_VecVecs.VecVec(1..dim)
         := Standard_Series_Coefficients(x);
    rx : constant Standard_Floating_VecVecs.Link_to_VecVec
       := Standard_Vector_Splitters.Allocate_Floating_Coefficients(dim,deg);
    ix : constant Standard_Floating_VecVecs.Link_to_VecVec
       := Standard_Vector_Splitters.Allocate_Floating_Coefficients(dim,deg);
    rpwt : constant Link_to_VecVecVec := Allocate(s.mxe,deg);
    ipwt : constant Link_to_VecVecVec := Allocate(s.mxe,deg);
    vy2 : constant Standard_Complex_VecVecs.VecVec(0..deg)
        := Linearized_Allocation(dim,deg);
    vm2 : constant Standard_Complex_VecMats.VecMat(0..deg)
        := Allocate_Coefficients(dim,dim,deg);
    multstart,multstop,seristart,seristop : Ada.Calendar.Time;
    seri_elapsed,mult_elapsed,speedup,efficiency : Duration;
    err : double_float;
    ans : character;

    use Ada.Calendar; -- for the difference operator on Duration

  begin
    Standard_Vector_Splitters.Complex_Parts(xcff,rx,ix);
    put_line("Running on one task ...");
    seristart := Ada.Calendar.Clock;
    Compute(s.rpwt,s.ipwt,s.mxe,rx,ix);
    EvalDiff(c,rx.all,ix.all,s.rpwt,s.ipwt,s.ryd,s.iyd,s.vy,s.vm);
    seristop := Ada.Calendar.Clock;
    seri_elapsed := seristop - seristart;
    put_line("-> Elapsed time without multitasking : ");
    Time_Stamps.Write_Elapsed_Time(standard_output,seristart,seristop);
    loop
      new_line;
      put("Running with "); put(nbt,1); put_line(" tasks ...");
      multstart := Ada.Calendar.Clock;
      Standard_Multitasked_EvalDiff
        (nbt,s.crc,rx,ix,s.mxe,rpwt,ipwt,vy2,vm2,static,output);
      multstop := Ada.Calendar.Clock;
      mult_elapsed := multstop - multstart;
      err := Difference(s.vy,vy2);
      put("  the error of evaluation : "); put(err,3); new_line;
      err := Difference(s.vm,vm2);
      put("  the error of differentiation : "); put(err,3); new_line;
      put_line("-> Elapsed time on multitasking : ");
      Time_Stamps.Write_Elapsed_Time(standard_output,multstart,multstop);
      if seri_elapsed + 1.0 /= 1.0 then
        speedup := seri_elapsed/mult_elapsed;
        efficiency := speedup/duration(nbt);
        efficiency := duration(100)*efficiency;
        put("The speedup : "); duration_io.put(speedup,1,3);
        put("  efficiency : "); duration_io.put(efficiency,2,2); new_line;
      end if;
      put("Continue with different number of tasks ? (y/n) ");
      Ask_Yes_or_No(ans);
      exit when (ans /= 'y');
      put("Give the number of tasks : "); get(nbt);
    end loop;
  end Standard_Coefficient_Test;

  procedure DoblDobl_Test
              ( c : in DoblDobl_Speelpenning_Convolutions.Circuits;
                dim,deg : in integer32; nbt : in out integer32;
                output : in boolean := true ) is

    use DoblDobl_Speelpenning_Convolutions;

    x : constant DoblDobl_Complex_Series_Vectors.Vector(1..dim)
      := DoblDobl_Random_Series_Vectors.Random_Series_Vector(1,dim,deg);
    xcff : constant DoblDobl_Complex_VecVecs.VecVec(1..dim)
         := DoblDobl_Series_Coefficients(x);
    mxe : constant Standard_Integer_Vectors.Vector(1..dim)
        := Exponent_Maxima(c,dim);
    pwt : constant Link_to_VecVecVec := Allocate(mxe,deg);
    yd : constant DoblDobl_Complex_VecVecs.VecVec(1..dim+1)
       := Allocate_Coefficients(dim+1,deg);
    vy1 : constant DoblDobl_Complex_VecVecs.VecVec(0..deg)
        := Linearized_Allocation(dim,deg);
    vy2 : constant DoblDobl_Complex_VecVecs.VecVec(0..deg)
        := Linearized_Allocation(dim,deg);
    vm1 : constant DoblDobl_Complex_VecMats.VecMat(0..deg)
        := Allocate_Coefficients(dim,dim,deg);
    vm2 : constant DoblDobl_Complex_VecMats.VecMat(0..deg)
        := Allocate_Coefficients(dim,dim,deg);
    err : double_double;
    multstart,multstop,seristart,seristop : Ada.Calendar.Time;
    seri_elapsed,mult_elapsed,speedup,efficiency : Duration;
    ans : character;

    use Ada.Calendar; -- for the difference operator on Duration

  begin
    put_line("Running on one task ...");
    seristart := Ada.Calendar.Clock;
    Compute(pwt,mxe,xcff);
    EvalDiff(c,xcff,pwt,yd,vy1,vm1);
    seristop := Ada.Calendar.Clock;
    seri_elapsed := seristop - seristart;
    new_line;
    put_line("-> Elapsed time without multitasking : ");
    Time_Stamps.Write_Elapsed_Time(standard_output,seristart,seristop);
    loop
      new_line;
      put("Running with "); put(nbt,1); put_line(" tasks ...");
      multstart := Ada.Calendar.Clock;
      DoblDobl_Multitasked_EvalDiff(nbt,c,xcff,mxe,pwt,vy2,vm2,output);
      multstop := Ada.Calendar.Clock;
      mult_elapsed := multstop - multstart;
      err := Difference(vy1,vy2);
      put("  the error of evaluation : "); put(err,3); new_line;
      err := Difference(vm1,vm2);
      put("  the error of differentiation : "); put(err,3); new_line;
      put_line("-> Elapsed time on multitasking : ");
      Time_Stamps.Write_Elapsed_Time(standard_output,multstart,multstop);
      if seri_elapsed + 1.0 /= 1.0 then
        speedup := seri_elapsed/mult_elapsed;
        efficiency := speedup/duration(nbt);
        efficiency := duration(100)*efficiency;
        put("The speedup : "); duration_io.put(speedup,1,3);
        put("  efficiency : "); duration_io.put(efficiency,2,2); new_line;
      end if;
      put("Continue with different number of tasks ? (y/n) ");
      Ask_Yes_or_No(ans);
      exit when (ans /= 'y');
      put("Give the number of tasks : "); get(nbt);
    end loop;
  end DoblDobl_Test;

  procedure DoblDobl_Coefficient_Test
              ( c : in DoblDobl_Coefficient_Convolutions.Circuits;
                dim,deg : in integer32; nbt : in out integer32;
                static : in boolean; output : in boolean := true ) is

    use DoblDobl_Coefficient_Convolutions;

    s : constant DoblDobl_Coefficient_Convolutions.Link_to_System
      := DoblDobl_Coefficient_Convolutions.Create(c,dim,deg);
    x : constant DoblDobl_Complex_Series_Vectors.Vector(1..dim)
      := DoblDobl_Random_Series_Vectors.Random_Series_Vector(1,dim,deg);
    xcff : constant DoblDobl_Complex_VecVecs.VecVec(1..dim)
         := DoblDobl_Series_Coefficients(x);
    rhx : constant Standard_Floating_VecVecs.Link_to_VecVec
        := Standard_Vector_Splitters.Allocate_Floating_Coefficients(dim,deg);
    ihx : constant Standard_Floating_VecVecs.Link_to_VecVec
        := Standard_Vector_Splitters.Allocate_Floating_Coefficients(dim,deg);
    rlx : constant Standard_Floating_VecVecs.Link_to_VecVec
        := Standard_Vector_Splitters.Allocate_Floating_Coefficients(dim,deg);
    ilx : constant Standard_Floating_VecVecs.Link_to_VecVec
        := Standard_Vector_Splitters.Allocate_Floating_Coefficients(dim,deg);
    rhpwt : constant Standard_Floating_VecVecVecs.Link_to_VecVecVec
          := Standard_Coefficient_Convolutions.Allocate(s.mxe,deg);
    rlpwt : constant Standard_Floating_VecVecVecs.Link_to_VecVecVec
          := Standard_Coefficient_Convolutions.Allocate(s.mxe,deg);
    ihpwt : constant Standard_Floating_VecVecVecs.Link_to_VecVecVec
          := Standard_Coefficient_Convolutions.Allocate(s.mxe,deg);
    ilpwt : constant Standard_Floating_VecVecVecs.Link_to_VecVecVec
          := Standard_Coefficient_Convolutions.Allocate(s.mxe,deg);
    vy2 : constant DoblDobl_Complex_VecVecs.VecVec(0..deg)
        := Linearized_Allocation(dim,deg);
    vm2 : constant DoblDobl_Complex_VecMats.VecMat(0..deg)
        := Allocate_Coefficients(dim,dim,deg);
    multstart,multstop,seristart,seristop : Ada.Calendar.Time;
    seri_elapsed,mult_elapsed,speedup,efficiency : Duration;
    err : double_double;
    ans : character;

    use Ada.Calendar; -- for the difference operator on Duration

  begin
    DoblDobl_Vector_Splitters.Complex_Parts(xcff,rhx,ihx,rlx,ilx);
    put_line("Running on one task ...");
    seristart := Ada.Calendar.Clock;
    Compute(s.rhpwt,s.ihpwt,s.rlpwt,s.ilpwt,s.mxe,rhx,ihx,rlx,ilx);
    EvalDiff(c,rhx.all,ihx.all,rlx.all,ilx.all,s.rhpwt,s.ihpwt,
             s.rlpwt,s.ilpwt,s.rhyd,s.ihyd,s.rlyd,s.ilyd,s.vy,s.vm);
    seristop := Ada.Calendar.Clock;
    seri_elapsed := seristop - seristart;
    put_line("-> Elapsed time without multitasking : ");
    Time_Stamps.Write_Elapsed_Time(standard_output,seristart,seristop);
    loop
      new_line;
      put("Running with "); put(nbt,1); put_line(" tasks ...");
      multstart := Ada.Calendar.Clock;
      DoblDobl_Multitasked_EvalDiff
        (nbt,s.crc,rhx,ihx,rlx,ilx,s.mxe,rhpwt,ihpwt,rlpwt,ilpwt,
         vy2,vm2,static,output);
      multstop := Ada.Calendar.Clock;
      mult_elapsed := multstop - multstart;
      err := Difference(s.vy,vy2);
      put("  the error of evaluation : "); put(err,3); new_line;
      err := Difference(s.vm,vm2);
      put("  the error of differentiation : "); put(err,3); new_line;
      put_line("-> Elapsed time on multitasking : ");
      Time_Stamps.Write_Elapsed_Time(standard_output,multstart,multstop);
      if seri_elapsed + 1.0 /= 1.0 then
        speedup := seri_elapsed/mult_elapsed;
        efficiency := speedup/duration(nbt);
        efficiency := duration(100)*efficiency;
        put("The speedup : "); duration_io.put(speedup,1,3);
        put("  efficiency : "); duration_io.put(efficiency,2,2); new_line;
      end if;
      put("Continue with different number of tasks ? (y/n) ");
      Ask_Yes_or_No(ans);
      exit when (ans /= 'y');
      put("Give the number of tasks : "); get(nbt);
    end loop;
  end DoblDobl_Coefficient_Test;

  procedure TripDobl_Test
              ( c : in TripDobl_Speelpenning_Convolutions.Circuits;
                dim,deg : in integer32; nbt : in out integer32;
                output : in boolean := true ) is

    use TripDobl_Speelpenning_Convolutions;

    x : constant TripDobl_Complex_Series_Vectors.Vector(1..dim)
      := TripDobl_Random_Series_Vectors.Random_Series_Vector(1,dim,deg);
    xcff : constant TripDobl_Complex_VecVecs.VecVec(1..dim)
         := TripDobl_Series_Coefficients(x);
    mxe : constant Standard_Integer_Vectors.Vector(1..dim)
        := Exponent_Maxima(c,dim);
    pwt : constant Link_to_VecVecVec := Allocate(mxe,deg);
    yd : constant TripDobl_Complex_VecVecs.VecVec(1..dim+1)
       := Allocate_Coefficients(dim+1,deg);
    vy1 : constant TripDobl_Complex_VecVecs.VecVec(0..deg)
        := Linearized_Allocation(dim,deg);
    vy2 : constant TripDobl_Complex_VecVecs.VecVec(0..deg)
        := Linearized_Allocation(dim,deg);
    vm1 : constant TripDobl_Complex_VecMats.VecMat(0..deg)
        := Allocate_Coefficients(dim,dim,deg);
    vm2 : constant TripDobl_Complex_VecMats.VecMat(0..deg)
        := Allocate_Coefficients(dim,dim,deg);
    err : triple_double;
    multstart,multstop,seristart,seristop : Ada.Calendar.Time;
    seri_elapsed,mult_elapsed,speedup,efficiency : Duration;
    ans : character;

    use Ada.Calendar; -- for the difference operator on Duration

  begin
    put_line("Running on one task ...");
    seristart := Ada.Calendar.Clock;
    Compute(pwt,mxe,xcff);
    EvalDiff(c,xcff,pwt,yd,vy1,vm1);
    seristop := Ada.Calendar.Clock;
    seri_elapsed := seristop - seristart;
    put_line("-> Elapsed time without multitasking : ");
    Time_Stamps.Write_Elapsed_Time(standard_output,seristart,seristop);
    loop
      new_line;
      put("Running with "); put(nbt,1); put_line(" tasks ...");
      multstart := Ada.Calendar.Clock;
      TripDobl_Multitasked_EvalDiff(nbt,c,xcff,mxe,pwt,vy2,vm2,output);
      multstop := Ada.Calendar.Clock;
      mult_elapsed := multstop - multstart;
      err := Difference(vy1,vy2);
      put("  the error of evaluation : "); put(err,3); new_line;
      err := Difference(vm1,vm2);
      put("  the error of differentiation : "); put(err,3); new_line;
      new_line;
      put_line("-> Elapsed time on multitasking : ");
      Time_Stamps.Write_Elapsed_Time(standard_output,multstart,multstop);
      if seri_elapsed + 1.0 /= 1.0 then
        speedup := seri_elapsed/mult_elapsed;
        efficiency := speedup/duration(nbt);
        efficiency := duration(100)*efficiency;
        put("The speedup : "); duration_io.put(speedup,1,3);
        put("  efficiency : "); duration_io.put(efficiency,2,2); new_line;
      end if;
      put("Continue with different number of tasks ? (y/n) ");
      Ask_Yes_or_No(ans);
      exit when (ans /= 'y');
      put("Give the number of tasks : "); get(nbt);
    end loop;
  end TripDobl_Test;

  procedure QuadDobl_Test
              ( c : in QuadDobl_Speelpenning_Convolutions.Circuits;
                dim,deg : in integer32; nbt : in out integer32;
                output : in boolean := true ) is

    use QuadDobl_Speelpenning_Convolutions;

    x : constant QuadDobl_Complex_Series_Vectors.Vector(1..dim)
      := QuadDobl_Random_Series_Vectors.Random_Series_Vector(1,dim,deg);
    xcff : constant QuadDobl_Complex_VecVecs.VecVec(1..dim)
         := QuadDobl_Series_Coefficients(x);
    mxe : constant Standard_Integer_Vectors.Vector(1..dim)
        := Exponent_Maxima(c,dim);
    pwt : constant Link_to_VecVecVec := Allocate(mxe,deg);
    yd : constant QuadDobl_Complex_VecVecs.VecVec(1..dim+1)
       := Allocate_Coefficients(dim+1,deg);
    vy1 : constant QuadDobl_Complex_VecVecs.VecVec(0..deg)
        := Linearized_Allocation(dim,deg);
    vy2 : constant QuadDobl_Complex_VecVecs.VecVec(0..deg)
        := Linearized_Allocation(dim,deg);
    vm1 : constant QuadDobl_Complex_VecMats.VecMat(0..deg)
        := Allocate_Coefficients(dim,dim,deg);
    vm2 : constant QuadDobl_Complex_VecMats.VecMat(0..deg)
        := Allocate_Coefficients(dim,dim,deg);
    err : quad_double;
    multstart,multstop,seristart,seristop : Ada.Calendar.Time;
    seri_elapsed,mult_elapsed,speedup,efficiency : Duration;
    ans : character;

    use Ada.Calendar; -- for the difference operator on Duration

  begin
    put_line("Running on one task ...");
    seristart := Ada.Calendar.Clock;
    Compute(pwt,mxe,xcff);
    EvalDiff(c,xcff,pwt,yd,vy1,vm1);
    seristop := Ada.Calendar.Clock;
    seri_elapsed := seristop - seristart;
    put_line("-> Elapsed time without multitasking : ");
    Time_Stamps.Write_Elapsed_Time(standard_output,seristart,seristop);
    loop
      new_line;
      put("Running with "); put(nbt,1); put_line(" tasks ...");
      multstart := Ada.Calendar.Clock;
      QuadDobl_Multitasked_EvalDiff(nbt,c,xcff,mxe,pwt,vy2,vm2,output);
      multstop := Ada.Calendar.Clock;
      mult_elapsed := multstop - multstart;
      err := Difference(vy1,vy2);
      put("  the error of evaluation : "); put(err,3); new_line;
      err := Difference(vm1,vm2);
      put("  the error of differentiation : "); put(err,3); new_line;
      new_line;
      put_line("-> Elapsed time on multitasking : ");
      Time_Stamps.Write_Elapsed_Time(standard_output,multstart,multstop);
      if seri_elapsed + 1.0 /= 1.0 then
        speedup := seri_elapsed/mult_elapsed;
        efficiency := speedup/duration(nbt);
        efficiency := duration(100)*efficiency;
        put("The speedup : "); duration_io.put(speedup,1,3);
        put("  efficiency : "); duration_io.put(efficiency,2,2); new_line;
      end if;
      put("Continue with different number of tasks ? (y/n) ");
      Ask_Yes_or_No(ans);
      exit when (ans /= 'y');
      put("Give the number of tasks : "); get(nbt);
    end loop;
  end QuadDobl_Test;

  procedure QuadDobl_Coefficient_Test
              ( c : in QuadDobl_Coefficient_Convolutions.Circuits;
                dim,deg : in integer32; nbt : in out integer32;
                static : in boolean; output : in boolean := true ) is

    use QuadDobl_Coefficient_Convolutions;

    s : constant QuadDobl_Coefficient_Convolutions.Link_to_System
      := QuadDobl_Coefficient_Convolutions.Create(c,dim,deg);
    x : constant QuadDobl_Complex_Series_Vectors.Vector(1..dim)
      := QuadDobl_Random_Series_Vectors.Random_Series_Vector(1,dim,deg);
    xcff : constant QuadDobl_Complex_VecVecs.VecVec(1..dim)
         := QuadDobl_Series_Coefficients(x);
    degdim : constant integer32 := 4*(deg+1)-1;
    xr : constant Standard_Floating_VecVecs.Link_to_VecVec
       := Standard_Vector_Splitters.Allocate_Floating_Coefficients(dim,degdim);
    xi : constant Standard_Floating_VecVecs.Link_to_VecVec
       := Standard_Vector_Splitters.Allocate_Floating_Coefficients(dim,degdim);
    rpwt : constant Standard_Floating_VecVecVecs.Link_to_VecVecVec
         := Standard_Coefficient_Convolutions.Allocate(s.mxe,degdim);
    ipwt : constant Standard_Floating_VecVecVecs.Link_to_VecVecVec
         := Standard_Coefficient_Convolutions.Allocate(s.mxe,degdim);
    vy2 : constant QuadDobl_Complex_VecVecs.VecVec(0..deg)
        := Linearized_Allocation(dim,deg);
    vm2 : constant QuadDobl_Complex_VecMats.VecMat(0..deg)
        := Allocate_Coefficients(dim,dim,deg);
    multstart,multstop,seristart,seristop : Ada.Calendar.Time;
    seri_elapsed,mult_elapsed,speedup,efficiency : Duration;
    err : quad_double;
    ans : character;
    ur,vr,wr : constant Standard_Floating_Vectors.Vector(0..3)
             := (0..3 => 0.0);
    u : constant Standard_Floating_Vectors.Link_to_Vector
      := new Standard_Floating_Vectors.Vector'(ur);
    v : constant Standard_Floating_Vectors.Link_to_Vector
      := new Standard_Floating_Vectors.Vector'(vr);
    w : constant Standard_Floating_Vectors.Link_to_Vector
      := new Standard_Floating_Vectors.Vector'(wr);

    use Ada.Calendar; -- for the difference operator on Duration

  begin
    QuadDobl_Vector_Splitters.Complex_Parts(xcff,xr,xi);
    put("mxe : "); put(s.mxe); new_line;
    put_line("Running on one task ...");
    seristart := Ada.Calendar.Clock;
    Compute(s.rpwt,s.ipwt,s.mxe,xr,xi,u,v,w);
    EvalDiff(c,xr.all,xi.all,s.rpwt,s.ipwt,s.ryd,s.iyd,s.vy,s.vm,u,v,w);
    seristop := Ada.Calendar.Clock;
    seri_elapsed := seristop - seristart;
    put_line("-> Elapsed time without multitasking : ");
    Time_Stamps.Write_Elapsed_Time(standard_output,seristart,seristop);
    loop
      new_line;
      put("Running with "); put(nbt,1); put_line(" tasks ...");
      multstart := Ada.Calendar.Clock;
      QuadDobl_Multitasked_EvalDiff
        (nbt,s.crc,xr,xi,s.mxe,rpwt,ipwt,vy2,vm2,static,output);
      multstop := Ada.Calendar.Clock;
      mult_elapsed := multstop - multstart;
      err := Difference(s.vy,vy2);
      put("  the error of evaluation : "); put(err,3); new_line;
      err := Difference(s.vm,vm2);
      put("  the error of differentiation : "); put(err,3); new_line;
      put_line("-> Elapsed time on multitasking : ");
      Time_Stamps.Write_Elapsed_Time(standard_output,multstart,multstop);
      if seri_elapsed + 1.0 /= 1.0 then
        speedup := seri_elapsed/mult_elapsed;
        efficiency := speedup/duration(nbt);
        efficiency := duration(100)*efficiency;
        put("The speedup : "); duration_io.put(speedup,1,3);
        put("  efficiency : "); duration_io.put(efficiency,2,2); new_line;
      end if;
      put("Continue with different number of tasks ? (y/n) ");
      Ask_Yes_or_No(ans);
      exit when (ans /= 'y');
      put("Give the number of tasks : "); get(nbt);
    end loop;
  end QuadDobl_Coefficient_Test;

  procedure PentDobl_Test
              ( c : in PentDobl_Speelpenning_Convolutions.Circuits;
                dim,deg : in integer32; nbt : in out integer32;
                output : in boolean := true ) is

    use PentDobl_Speelpenning_Convolutions;

    x : constant PentDobl_Complex_Series_Vectors.Vector(1..dim)
      := PentDobl_Random_Series_Vectors.Random_Series_Vector(1,dim,deg);
    xcff : constant PentDobl_Complex_VecVecs.VecVec(1..dim)
         := PentDobl_Series_Coefficients(x);
    mxe : constant Standard_Integer_Vectors.Vector(1..dim)
        := Exponent_Maxima(c,dim);
    pwt : constant Link_to_VecVecVec := Allocate(mxe,deg);
    yd : constant PentDobl_Complex_VecVecs.VecVec(1..dim+1)
       := Allocate_Coefficients(dim+1,deg);
    vy1 : constant PentDobl_Complex_VecVecs.VecVec(0..deg)
        := Linearized_Allocation(dim,deg);
    vy2 : constant PentDobl_Complex_VecVecs.VecVec(0..deg)
        := Linearized_Allocation(dim,deg);
    vm1 : constant PentDobl_Complex_VecMats.VecMat(0..deg)
        := Allocate_Coefficients(dim,dim,deg);
    vm2 : constant PentDobl_Complex_VecMats.VecMat(0..deg)
        := Allocate_Coefficients(dim,dim,deg);
    err : penta_double;
    multstart,multstop,seristart,seristop : Ada.Calendar.Time;
    seri_elapsed,mult_elapsed,speedup,efficiency : Duration;
    ans : character;

    use Ada.Calendar; -- for the difference operator on Duration

  begin
    put_line("Running on one task ...");
    seristart := Ada.Calendar.Clock;
    Compute(pwt,mxe,xcff);
    EvalDiff(c,xcff,pwt,yd,vy1,vm1);
    seristop := Ada.Calendar.Clock;
    seri_elapsed := seristop - seristart;
    put_line("-> Elapsed time without multitasking : ");
    Time_Stamps.Write_Elapsed_Time(standard_output,seristart,seristop);
    loop
      new_line;
      put("Running with "); put(nbt,1); put_line(" tasks ...");
      multstart := Ada.Calendar.Clock;
      PentDobl_Multitasked_EvalDiff(nbt,c,xcff,mxe,pwt,vy2,vm2,output);
      multstop := Ada.Calendar.Clock;
      mult_elapsed := multstop - multstart;
      err := Difference(vy1,vy2);
      put("  the error of evaluation : "); put(err,3); new_line;
      err := Difference(vm1,vm2);
      put("  the error of differentiation : "); put(err,3); new_line;
      new_line;
      put_line("-> Elapsed time on multitasking : ");
      Time_Stamps.Write_Elapsed_Time(standard_output,multstart,multstop);
      if seri_elapsed + 1.0 /= 1.0 then
        speedup := seri_elapsed/mult_elapsed;
        efficiency := speedup/duration(nbt);
        efficiency := duration(100)*efficiency;
        put("The speedup : "); duration_io.put(speedup,1,3);
        put("  efficiency : "); duration_io.put(efficiency,2,2); new_line;
      end if;
      put("Continue with different number of tasks ? (y/n) ");
      Ask_Yes_or_No(ans);
      exit when (ans /= 'y');
      put("Give the number of tasks : "); get(nbt);
    end loop;
  end PentDobl_Test;

  procedure OctoDobl_Test
              ( c : in OctoDobl_Speelpenning_Convolutions.Circuits;
                dim,deg : in integer32; nbt : in out integer32;
                output : in boolean := true ) is

    use OctoDobl_Speelpenning_Convolutions;

    x : constant OctoDobl_Complex_Series_Vectors.Vector(1..dim)
      := OctoDobl_Random_Series_Vectors.Random_Series_Vector(1,dim,deg);
    xcff : constant OctoDobl_Complex_VecVecs.VecVec(1..dim)
         := OctoDobl_Series_Coefficients(x);
    mxe : constant Standard_Integer_Vectors.Vector(1..dim)
        := Exponent_Maxima(c,dim);
    pwt : constant Link_to_VecVecVec := Allocate(mxe,deg);
    yd : constant OctoDobl_Complex_VecVecs.VecVec(1..dim+1)
       := Allocate_Coefficients(dim+1,deg);
    vy1 : constant OctoDobl_Complex_VecVecs.VecVec(0..deg)
        := Linearized_Allocation(dim,deg);
    vy2 : constant OctoDobl_Complex_VecVecs.VecVec(0..deg)
        := Linearized_Allocation(dim,deg);
    vm1 : constant OctoDobl_Complex_VecMats.VecMat(0..deg)
        := Allocate_Coefficients(dim,dim,deg);
    vm2 : constant OctoDobl_Complex_VecMats.VecMat(0..deg)
        := Allocate_Coefficients(dim,dim,deg);
    err : octo_double;
    multstart,multstop,seristart,seristop : Ada.Calendar.Time;
    seri_elapsed,mult_elapsed,speedup,efficiency : Duration;
    ans : character;

    use Ada.Calendar; -- for the difference operator on Duration

  begin
    put_line("Running on one task ...");
    seristart := Ada.Calendar.Clock;
    Compute(pwt,mxe,xcff);
    EvalDiff(c,xcff,pwt,yd,vy1,vm1);
    seristop := Ada.Calendar.Clock;
    seri_elapsed := seristop - seristart;
    put_line("-> Elapsed time without multitasking : ");
    Time_Stamps.Write_Elapsed_Time(standard_output,seristart,seristop);
    loop
      new_line;
      put("Running with "); put(nbt,1); put_line(" tasks ...");
      multstart := Ada.Calendar.Clock;
      OctoDobl_Multitasked_EvalDiff(nbt,c,xcff,mxe,pwt,vy2,vm2,output);
      multstop := Ada.Calendar.Clock;
      mult_elapsed := multstop - multstart;
      err := Difference(vy1,vy2);
      put("  the error of evaluation : "); put(err,3); new_line;
      err := Difference(vm1,vm2);
      put("  the error of differentiation : "); put(err,3); new_line;
      new_line;
      put_line("-> Elapsed time on multitasking : ");
      Time_Stamps.Write_Elapsed_Time(standard_output,multstart,multstop);
      if seri_elapsed + 1.0 /= 1.0 then
        speedup := seri_elapsed/mult_elapsed;
        efficiency := speedup/duration(nbt);
        efficiency := duration(100)*efficiency;
        put("The speedup : "); duration_io.put(speedup,1,3);
        put("  efficiency : "); duration_io.put(efficiency,2,2); new_line;
      end if;
      put("Continue with different number of tasks ? (y/n) ");
      Ask_Yes_or_No(ans);
      exit when (ans /= 'y');
      put("Give the number of tasks : "); get(nbt);
    end loop;
  end OctoDobl_Test;

  procedure DecaDobl_Test
              ( c : in DecaDobl_Speelpenning_Convolutions.Circuits;
                dim,deg : in integer32; nbt : in out integer32;
                output : in boolean := true ) is

    use DecaDobl_Speelpenning_Convolutions;

    x : constant DecaDobl_Complex_Series_Vectors.Vector(1..dim)
      := DecaDobl_Random_Series_Vectors.Random_Series_Vector(1,dim,deg);
    xcff : constant DecaDobl_Complex_VecVecs.VecVec(1..dim)
         := DecaDobl_Series_Coefficients(x);
    mxe : constant Standard_Integer_Vectors.Vector(1..dim)
        := Exponent_Maxima(c,dim);
    pwt : constant Link_to_VecVecVec := Allocate(mxe,deg);
    yd : constant DecaDobl_Complex_VecVecs.VecVec(1..dim+1)
       := Allocate_Coefficients(dim+1,deg);
    vy1 : constant DecaDobl_Complex_VecVecs.VecVec(0..deg)
        := Linearized_Allocation(dim,deg);
    vy2 : constant DecaDobl_Complex_VecVecs.VecVec(0..deg)
        := Linearized_Allocation(dim,deg);
    vm1 : constant DecaDobl_Complex_VecMats.VecMat(0..deg)
        := Allocate_Coefficients(dim,dim,deg);
    vm2 : constant DecaDobl_Complex_VecMats.VecMat(0..deg)
        := Allocate_Coefficients(dim,dim,deg);
    err : deca_double;
    multstart,multstop,seristart,seristop : Ada.Calendar.Time;
    seri_elapsed,mult_elapsed,speedup,efficiency : Duration;
    ans : character;

    use Ada.Calendar; -- for the difference operator on Duration

  begin
    put_line("Running on one task ...");
    seristart := Ada.Calendar.Clock;
    Compute(pwt,mxe,xcff);
    EvalDiff(c,xcff,pwt,yd,vy1,vm1);
    seristop := Ada.Calendar.Clock;
    seri_elapsed := seristop - seristart;
    put_line("-> Elapsed time without multitasking : ");
    Time_Stamps.Write_Elapsed_Time(standard_output,seristart,seristop);
    loop
      new_line;
      put("Running with "); put(nbt,1); put_line(" tasks ...");
      multstart := Ada.Calendar.Clock;
      DecaDobl_Multitasked_EvalDiff(nbt,c,xcff,mxe,pwt,vy2,vm2,output);
      multstop := Ada.Calendar.Clock;
      mult_elapsed := multstop - multstart;
      err := Difference(vy1,vy2);
      put("  the error of evaluation : "); put(err,3); new_line;
      err := Difference(vm1,vm2);
      put("  the error of differentiation : "); put(err,3); new_line;
      new_line;
      put_line("-> Elapsed time on multitasking : ");
      Time_Stamps.Write_Elapsed_Time(standard_output,multstart,multstop);
      if seri_elapsed + 1.0 /= 1.0 then
        speedup := seri_elapsed/mult_elapsed;
        efficiency := speedup/duration(nbt);
        efficiency := duration(100)*efficiency;
        put("The speedup : "); duration_io.put(speedup,1,3);
        put("  efficiency : "); duration_io.put(efficiency,2,2); new_line;
      end if;
      put("Continue with different number of tasks ? (y/n) ");
      Ask_Yes_or_No(ans);
      exit when (ans /= 'y');
      put("Give the number of tasks : "); get(nbt);
    end loop;
  end DecaDobl_Test;

  procedure Standard_Random_Test
              ( dim,deg,nbr,pwr : in integer32; nbt : in out integer32;
                output : in boolean := true ) is

    c : constant Standard_Speelpenning_Convolutions.Circuits
      := Standard_Random_Convolution_Circuits(dim,deg,nbr,pwr);
    cc : constant Standard_Coefficient_Convolutions.Circuits
       := Standard_Convolution_Splitters.Split(c);
    ans : character;
    static : boolean;

  begin
    put("Run on circuits with splitted coefficient vectors ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans = 'y' then
      put("Static load balancing ? (y/n) "); Ask_Yes_or_No(ans);
      static := (ans = 'y');
      Standard_Coefficient_Test(cc,dim,deg,nbt,static,output);
    else
      Standard_Test(c,dim,deg,nbt,output);
    end if;
  end Standard_Random_Test;

  procedure DoblDobl_Random_Test
              ( dim,deg,nbr,pwr : in integer32; nbt : in out integer32;
                output : in boolean := true ) is

    c : constant DoblDobl_Speelpenning_Convolutions.Circuits
      := DoblDobl_Random_Convolution_Circuits(dim,deg,nbr,pwr);
    cc : constant DoblDobl_Coefficient_Convolutions.Circuits
       := DoblDobl_Convolution_Splitters.Split(c);
    ans : character;
    static : boolean;

  begin
    put("Run on circuits with splitted coefficient vectors ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans = 'y' then
      put("Static load balancing ? (y/n) "); Ask_Yes_or_No(ans);
      static := (ans = 'y');
      DoblDobl_Coefficient_Test(cc,dim,deg,nbt,static,output);
    else
      DoblDobl_Test(c,dim,deg,nbt,output);
    end if;
  end DoblDobl_Random_Test;

  procedure TripDobl_Random_Test
              ( dim,deg,nbr,pwr : in integer32; nbt : in out integer32;
                output : in boolean := true ) is

    c : constant TripDobl_Speelpenning_Convolutions.Circuits
      := TripDobl_Random_Convolution_Circuits(dim,deg,nbr,pwr);

  begin
    TripDobl_Test(c,dim,deg,nbt,output);
  end TripDobl_Random_Test;

  procedure QuadDobl_Random_Test
              ( dim,deg,nbr,pwr : in integer32; nbt : in out integer32;
                output : in boolean := true ) is

    c : constant QuadDobl_Speelpenning_Convolutions.Circuits
      := QuadDobl_Random_Convolution_Circuits(dim,deg,nbr,pwr);
    cc : constant QuadDobl_Coefficient_Convolutions.Circuits
       := QuadDobl_Convolution_Splitters.Split(c);
    ans : character;
    static : boolean;

  begin
    put("Run on circuits with splitted coefficient vectors ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans = 'y' then
      put("Static load balancing ? (y/n) "); Ask_Yes_or_No(ans);
      static := (ans = 'y');
      QuadDobl_Coefficient_Test(cc,dim,deg,nbt,static,output);
    else
      QuadDobl_Test(c,dim,deg,nbt,output);
    end if;
  end QuadDobl_Random_Test;

  procedure PentDobl_Random_Test
              ( dim,deg,nbr,pwr : in integer32; nbt : in out integer32;
                output : in boolean := true ) is

    c : constant PentDobl_Speelpenning_Convolutions.Circuits
      := PentDobl_Random_Convolution_Circuits(dim,deg,nbr,pwr);

  begin
    PentDobl_Test(c,dim,deg,nbt,output);
  end PentDobl_Random_Test;

  procedure OctoDobl_Random_Test
              ( dim,deg,nbr,pwr : in integer32; nbt : in out integer32;
                output : in boolean := true ) is

    c : constant OctoDobl_Speelpenning_Convolutions.Circuits
      := OctoDobl_Random_Convolution_Circuits(dim,deg,nbr,pwr);

  begin
    OctoDobl_Test(c,dim,deg,nbt,output);
  end OctoDobl_Random_Test;

  procedure DecaDobl_Random_Test
              ( dim,deg,nbr,pwr : in integer32; nbt : in out integer32;
                output : in boolean := true ) is

    c : constant DecaDobl_Speelpenning_Convolutions.Circuits
      := DecaDobl_Random_Convolution_Circuits(dim,deg,nbr,pwr);

  begin
    DecaDobl_Test(c,dim,deg,nbt,output);
  end DecaDobl_Random_Test;

  procedure Standard_Benchmark
              ( file : in file_type; deg,nbruns,inc : in integer32;
                c : in Standard_Speelpenning_Convolutions.Circuits;
                x : in Standard_Complex_VecVecs.VecVec;
                mxe : in Standard_Integer_Vectors.Vector;
                verbose : in boolean := true ) is

    use Standard_Speelpenning_Convolutions;

    dim : constant integer32 := x'last;
    pwt : constant Link_to_VecVecVec := Allocate(mxe,deg);
    yd : constant Standard_Complex_VecVecs.VecVec(1..dim+1)
       := Allocate_Coefficients(dim+1,deg);
    vy1 : constant Standard_Complex_VecVecs.VecVec(0..deg)
        := Linearized_Allocation(dim,deg);
    vy2 : constant Standard_Complex_VecVecs.VecVec(0..deg)
        := Linearized_Allocation(dim,deg);
    vm1 : constant Standard_Complex_VecMats.VecMat(0..deg)
        := Allocate_Coefficients(dim,dim,deg);
    vm2 : constant Standard_Complex_VecMats.VecMat(0..deg)
        := Allocate_Coefficients(dim,dim,deg);
    multstart,multstop,seristart,seristop : Ada.Calendar.Time;
    seri_elapsed,mult_elapsed,speedup,efficiency : Duration;
    nbt : integer32 := 2;

    use Ada.Calendar;

  begin
    put_line(file,"double precision"); flush(file);
    if verbose
     then put_line("Running on one task ...");
    end if;
    seristart := Ada.Calendar.Clock;
    Compute(pwt,mxe,x);
    EvalDiff(c,x,pwt,yd,vy1,vm1);
    seristop := Ada.Calendar.Clock;
    seri_elapsed := seristop - seristart;
    put(file,"  1 : ");
    duration_io.put(file,seri_elapsed,1,3); new_line(file); flush(file);
    if verbose then
      put_line("-> Elapsed time without multitasking : ");
      Time_Stamps.Write_Elapsed_Time(standard_output,seristart,seristop);
    end if;   
    for k in 1..nbruns loop
      if verbose then
        new_line;
        put("Running with "); put(nbt,1); put_line(" tasks ...");
      end if;
      multstart := Ada.Calendar.Clock;
      Standard_Multitasked_EvalDiff(nbt,c,x,mxe,pwt,vy2,vm2,false);
      multstop := Ada.Calendar.Clock;
      mult_elapsed := multstop - multstart;
      if verbose then
        new_line;
        put_line("-> Elapsed time on multitasking : ");
        Time_Stamps.Write_Elapsed_Time(standard_output,multstart,multstop);
      end if;
      if seri_elapsed + 1.0 /= 1.0 then
        speedup := seri_elapsed/mult_elapsed;
        efficiency := speedup/duration(nbt);
        efficiency := duration(100)*efficiency;
        if verbose then
          put("The speedup : "); duration_io.put(speedup,1,3);
          put("  efficiency : "); duration_io.put(efficiency,2,2); new_line;
        end if;
      end if;
      put(file,nbt,3);
      put(file," : "); duration_io.put(file,mult_elapsed,1,3);
      put(file," : "); duration_io.put(file,speedup,1,3);
      put(file," : "); duration_io.put(file,efficiency,1,3);
      new_line(file); flush(file);
      nbt := nbt + inc;
    end loop;
  end Standard_Benchmark;

  procedure Standard_Coefficient_Benchmark
              ( file : in file_type; deg,nbruns,inc : in integer32;
                c : in Standard_Coefficient_Convolutions.Circuits;
                x : in Standard_Complex_VecVecs.VecVec;
                mxe : in Standard_Integer_Vectors.Vector;
                verbose : in boolean := true ) is

    use Standard_Floating_VecVecVecs;
    use Standard_Vector_Splitters;
    use Standard_Coefficient_Convolutions;

    dim : constant integer32 := x'last;
    rx : constant Standard_Floating_VecVecs.Link_to_VecVec
       := Allocate_Floating_Coefficients(dim,deg);
    ix : constant Standard_Floating_VecVecs.Link_to_VecVec
       := Allocate_Floating_Coefficients(dim,deg);
    rpwt : constant Link_to_VecVecVec := Allocate(mxe,deg);
    ipwt : constant Link_to_VecVecVec := Allocate(mxe,deg);
    ryd : constant Standard_Floating_VecVecs.Link_to_VecVec
        := Allocate_Floating_Coefficients(dim+1,deg);
    iyd : constant Standard_Floating_VecVecs.Link_to_VecVec
        := Allocate_Floating_Coefficients(dim+1,deg);
    vy1 : constant Standard_Complex_VecVecs.VecVec(0..deg)
        := Linearized_Allocation(dim,deg);
    vy2 : constant Standard_Complex_VecVecs.VecVec(0..deg)
        := Linearized_Allocation(dim,deg);
    vm1 : constant Standard_Complex_VecMats.VecMat(0..deg)
        := Allocate_Coefficients(dim,dim,deg);
    vm2 : constant Standard_Complex_VecMats.VecMat(0..deg)
        := Allocate_Coefficients(dim,dim,deg);
    multstart,multstop,seristart,seristop : Ada.Calendar.Time;
    seri_elapsed,mult_elapsed,speedup,efficiency : Duration;
    nbt : integer32 := 2;

    use Ada.Calendar;

  begin
    Standard_Vector_Splitters.Complex_Parts(x,rx,ix);
    put_line(file,"double precision"); flush(file);
    if verbose
     then put_line("Running on one task ...");
    end if;
    seristart := Ada.Calendar.Clock;
    Compute(rpwt,ipwt,mxe,rx,ix);
    EvalDiff(c,rx.all,ix.all,rpwt,ipwt,ryd.all,iyd.all,vy1,vm1);
    seristop := Ada.Calendar.Clock;
    seri_elapsed := seristop - seristart;
    put(file,"  1 : ");
    duration_io.put(file,seri_elapsed,1,3); new_line(file); flush(file);
    if verbose then
      put_line("-> Elapsed time without multitasking : ");
      Time_Stamps.Write_Elapsed_Time(standard_output,seristart,seristop);
    end if;   
    for k in 1..nbruns loop
      if verbose then
        new_line;
        put("Running with "); put(nbt,1); put_line(" tasks ...");
      end if;
      multstart := Ada.Calendar.Clock;
      Standard_Multitasked_EvalDiff
        (nbt,c,rx,ix,mxe,rpwt,ipwt,vy2,vm2,false,false);
      multstop := Ada.Calendar.Clock;
      mult_elapsed := multstop - multstart;
      if verbose then
        new_line;
        put_line("-> Elapsed time on multitasking : ");
        Time_Stamps.Write_Elapsed_Time(standard_output,multstart,multstop);
      end if;
      if seri_elapsed + 1.0 /= 1.0 then
        speedup := seri_elapsed/mult_elapsed;
        efficiency := speedup/duration(nbt);
        efficiency := duration(100)*efficiency;
        if verbose then
          put("The speedup : "); duration_io.put(speedup,1,3);
          put("   efficiency : "); duration_io.put(efficiency,2,2); new_line;
        end if;
      end if;
      put(file,nbt,3);
      put(file," : "); duration_io.put(file,mult_elapsed,1,3);
      put(file," : "); duration_io.put(file,speedup,1,3);
      put(file," : "); duration_io.put(file,efficiency,2,2);
      new_line(file); flush(file);
      nbt := nbt + inc;
    end loop;
  end Standard_Coefficient_Benchmark;

  procedure DoblDobl_Benchmark
              ( file : in file_type; deg,nbruns,inc : in integer32;
                c : in DoblDobl_Speelpenning_Convolutions.Circuits;
                x : in DoblDobl_Complex_VecVecs.VecVec;
                mxe : in Standard_Integer_Vectors.Vector;
                verbose : in boolean := true ) is

    use DoblDobl_Speelpenning_Convolutions;

    dim : constant integer32 := x'last;
    pwt : constant Link_to_VecVecVec := Allocate(mxe,deg);
    yd : constant DoblDobl_Complex_VecVecs.VecVec(1..dim+1)
       := Allocate_Coefficients(dim+1,deg);
    vy1 : constant DoblDobl_Complex_VecVecs.VecVec(0..deg)
        := Linearized_Allocation(dim,deg);
    vy2 : constant DoblDobl_Complex_VecVecs.VecVec(0..deg)
        := Linearized_Allocation(dim,deg);
    vm1 : constant DoblDobl_Complex_VecMats.VecMat(0..deg)
        := Allocate_Coefficients(dim,dim,deg);
    vm2 : constant DoblDobl_Complex_VecMats.VecMat(0..deg)
        := Allocate_Coefficients(dim,dim,deg);
    multstart,multstop,seristart,seristop : Ada.Calendar.Time;
    seri_elapsed,mult_elapsed,speedup,efficiency : Duration;
    nbt : integer32 := 2;

    use Ada.Calendar;

  begin
    put_line(file,"double double precision"); flush(file);
    if verbose
     then put_line("Running on one task ...");
    end if;
    seristart := Ada.Calendar.Clock;
    Compute(pwt,mxe,x);
    EvalDiff(c,x,pwt,yd,vy1,vm1);
    seristop := Ada.Calendar.Clock;
    seri_elapsed := seristop - seristart;
    put(file,"  1 : ");
    duration_io.put(file,seri_elapsed,1,3); new_line(file); flush(file);
    if verbose then
      put_line("-> Elapsed time without multitasking : ");
      Time_Stamps.Write_Elapsed_Time(standard_output,seristart,seristop);
    end if;   
    for k in 1..nbruns loop
      if verbose then
        new_line;
        put("Running with "); put(nbt,1); put_line(" tasks ...");
      end if;
      multstart := Ada.Calendar.Clock;
      DoblDobl_Multitasked_EvalDiff(nbt,c,x,mxe,pwt,vy2,vm2,false);
      multstop := Ada.Calendar.Clock;
      mult_elapsed := multstop - multstart;
      if verbose then
        new_line;
        put_line("-> Elapsed time on multitasking : ");
        Time_Stamps.Write_Elapsed_Time(standard_output,multstart,multstop);
      end if;
      if seri_elapsed + 1.0 /= 1.0 then
        speedup := seri_elapsed/mult_elapsed;
        efficiency := speedup/duration(nbt);
        efficiency := duration(100)*efficiency;
        if verbose then
          put("The speedup : "); duration_io.put(speedup,1,3);
          put("  efficiency : "); duration_io.put(speedup,2,2); new_line;
        end if;
      end if;
      put(file,nbt,3);
      put(file," : "); duration_io.put(file,mult_elapsed,1,3);
      put(file," : "); duration_io.put(file,speedup,1,3);
      put(file," : "); duration_io.put(file,efficiency,2,2);
      new_line(file); flush(file);
      nbt := nbt + inc;
    end loop;
  end DoblDobl_Benchmark;

  procedure DoblDobl_Coefficient_Benchmark
              ( file : in file_type; deg,nbruns,inc : in integer32;
                c : in DoblDobl_Coefficient_Convolutions.Circuits;
                x : in DoblDobl_Complex_VecVecs.VecVec;
                mxe : in Standard_Integer_Vectors.Vector;
                verbose : in boolean := true ) is

    use Standard_Vector_Splitters;
    use DoblDobl_Coefficient_Convolutions;

    dim : constant integer32 := x'last;
    rhx : constant Standard_Floating_VecVecs.Link_to_VecVec
        := Standard_Vector_Splitters.Allocate_Floating_Coefficients(dim,deg);
    ihx : constant Standard_Floating_VecVecs.Link_to_VecVec
        := Standard_Vector_Splitters.Allocate_Floating_Coefficients(dim,deg);
    rlx : constant Standard_Floating_VecVecs.Link_to_VecVec
        := Standard_Vector_Splitters.Allocate_Floating_Coefficients(dim,deg);
    ilx : constant Standard_Floating_VecVecs.Link_to_VecVec
        := Standard_Vector_Splitters.Allocate_Floating_Coefficients(dim,deg);
    rhpwt : constant Standard_Floating_VecVecVecs.Link_to_VecVecVec
          := Standard_Coefficient_Convolutions.Allocate(mxe,deg);
    rlpwt : constant Standard_Floating_VecVecVecs.Link_to_VecVecVec
          := Standard_Coefficient_Convolutions.Allocate(mxe,deg);
    ihpwt : constant Standard_Floating_VecVecVecs.Link_to_VecVecVec
          := Standard_Coefficient_Convolutions.Allocate(mxe,deg);
    ilpwt : constant Standard_Floating_VecVecVecs.Link_to_VecVecVec
          := Standard_Coefficient_Convolutions.Allocate(mxe,deg);
   -- yd : constant DoblDobl_Complex_VecVecs.VecVec(1..dim+1)
   --    := Allocate_Coefficients(dim+1,deg);
    rhyd : constant Standard_Floating_VecVecs.Link_to_VecVec
         := Allocate_Floating_Coefficients(dim+1,deg);
    ihyd : constant Standard_Floating_VecVecs.Link_to_VecVec
         := Allocate_Floating_Coefficients(dim+1,deg);
    rlyd : constant Standard_Floating_VecVecs.Link_to_VecVec
         := Allocate_Floating_Coefficients(dim+1,deg);
    ilyd : constant Standard_Floating_VecVecs.Link_to_VecVec
         := Allocate_Floating_Coefficients(dim+1,deg);
    vy1 : constant DoblDobl_Complex_VecVecs.VecVec(0..deg)
        := Linearized_Allocation(dim,deg);
    vy2 : constant DoblDobl_Complex_VecVecs.VecVec(0..deg)
        := Linearized_Allocation(dim,deg);
    vm1 : constant DoblDobl_Complex_VecMats.VecMat(0..deg)
        := Allocate_Coefficients(dim,dim,deg);
    vm2 : constant DoblDobl_Complex_VecMats.VecMat(0..deg)
        := Allocate_Coefficients(dim,dim,deg);
    multstart,multstop,seristart,seristop : Ada.Calendar.Time;
    seri_elapsed,mult_elapsed,speedup,efficiency : Duration;
    nbt : integer32 := 2;

    use Ada.Calendar;

  begin
    DoblDobl_Vector_Splitters.Complex_Parts(x,rhx,ihx,rlx,ilx);
    put_line(file,"double double precision"); flush(file);
    if verbose
     then put_line("Running on one task ...");
    end if;
    seristart := Ada.Calendar.Clock;
    Compute(rhpwt,ihpwt,rlpwt,ilpwt,mxe,rhx,ihx,rlx,ilx);
    EvalDiff(c,rhx.all,ihx.all,rlx.all,ilx.all,rhpwt,ihpwt,rlpwt,ilpwt,
             rhyd.all,ihyd.all,rlyd.all,ilyd.all,vy1,vm1);
    seristop := Ada.Calendar.Clock;
    seri_elapsed := seristop - seristart;
    put(file,"  1 : ");
    duration_io.put(file,seri_elapsed,1,3); new_line(file); flush(file);
    if verbose then
      put_line("-> Elapsed time without multitasking : ");
      Time_Stamps.Write_Elapsed_Time(standard_output,seristart,seristop);
    end if;   
    for k in 1..nbruns loop
      if verbose then
        new_line;
        put("Running with "); put(nbt,1); put_line(" tasks ...");
      end if;
      multstart := Ada.Calendar.Clock;
      DoblDobl_Multitasked_EvalDiff
        (nbt,c,rhx,ihx,rlx,ilx,mxe,rhpwt,ihpwt,rlpwt,ilpwt,vy2,vm2,
         false,false);
      multstop := Ada.Calendar.Clock;
      mult_elapsed := multstop - multstart;
      if verbose then
        new_line;
        put_line("-> Elapsed time on multitasking : ");
        Time_Stamps.Write_Elapsed_Time(standard_output,multstart,multstop);
      end if;
      if seri_elapsed + 1.0 /= 1.0 then
        speedup := seri_elapsed/mult_elapsed;
        efficiency := speedup/duration(nbt);
        efficiency := duration(100)*efficiency;
        if verbose then
          put("The speedup : "); duration_io.put(speedup,1,3);
          put("  efficiency : "); duration_io.put(speedup,2,2); new_line;
        end if;
      end if;
      put(file,nbt,3);
      put(file," : "); duration_io.put(file,mult_elapsed,1,3);
      put(file," : "); duration_io.put(file,speedup,1,3);
      put(file," : "); duration_io.put(file,efficiency,2,2);
      new_line(file); flush(file);
      nbt := nbt + inc;
    end loop;
  end DoblDobl_Coefficient_Benchmark;

  procedure TripDobl_Benchmark
              ( file : in file_type; deg,nbruns,inc : in integer32;
                c : in TripDobl_Speelpenning_Convolutions.Circuits;
                x : in TripDobl_Complex_VecVecs.VecVec;
                mxe : in Standard_Integer_Vectors.Vector;
                verbose : in boolean := true ) is

    use TripDobl_Speelpenning_Convolutions;

    dim : constant integer32 := x'last;
    pwt : constant Link_to_VecVecVec := Allocate(mxe,deg);
    yd : constant TripDobl_Complex_VecVecs.VecVec(1..dim+1)
       := Allocate_Coefficients(dim+1,deg);
    vy1 : constant TripDobl_Complex_VecVecs.VecVec(0..deg)
        := Linearized_Allocation(dim,deg);
    vy2 : constant TripDobl_Complex_VecVecs.VecVec(0..deg)
        := Linearized_Allocation(dim,deg);
    vm1 : constant TripDobl_Complex_VecMats.VecMat(0..deg)
        := Allocate_Coefficients(dim,dim,deg);
    vm2 : constant TripDobl_Complex_VecMats.VecMat(0..deg)
        := Allocate_Coefficients(dim,dim,deg);
    multstart,multstop,seristart,seristop : Ada.Calendar.Time;
    seri_elapsed,mult_elapsed,speedup,efficiency : Duration;
    nbt : integer32 := 2;

    use Ada.Calendar;

  begin
    put_line(file,"triple double precision"); flush(file);
    if verbose
     then put_line("Running on one task ...");
    end if;
    seristart := Ada.Calendar.Clock;
    Compute(pwt,mxe,x);
    EvalDiff(c,x,pwt,yd,vy1,vm1);
    seristop := Ada.Calendar.Clock;
    seri_elapsed := seristop - seristart;
    put(file,"  1 : ");
    duration_io.put(file,seri_elapsed,1,3); new_line(file); flush(file);
    if verbose then
      put_line("-> Elapsed time without multitasking : ");
      Time_Stamps.Write_Elapsed_Time(standard_output,seristart,seristop);
    end if;   
    for k in 1..nbruns loop
      if verbose then
        new_line;
        put("Running with "); put(nbt,1); put_line(" tasks ...");
      end if;
      multstart := Ada.Calendar.Clock;
      TripDobl_Multitasked_EvalDiff(nbt,c,x,mxe,pwt,vy2,vm2,false);
      multstop := Ada.Calendar.Clock;
      mult_elapsed := multstop - multstart;
      if verbose then
        new_line;
        put_line("-> Elapsed time on multitasking : ");
        Time_Stamps.Write_Elapsed_Time(standard_output,multstart,multstop);
      end if;
      if seri_elapsed + 1.0 /= 1.0 then
        speedup := seri_elapsed/mult_elapsed;
        efficiency := speedup/duration(nbt);
        efficiency := duration(100)*efficiency;
        if verbose then
          put("The speedup : "); duration_io.put(speedup,1,3);
          put("  efficiency : "); duration_io.put(speedup,2,2); new_line;
        end if;
      end if;
      put(file,nbt,3);
      put(file," : "); duration_io.put(file,mult_elapsed,1,3);
      put(file," : "); duration_io.put(file,speedup,1,3);
      put(file," : "); duration_io.put(file,efficiency,2,2);
      new_line(file); flush(file);
      nbt := nbt + inc;
    end loop;
  end TripDobl_Benchmark;

  procedure QuadDobl_Benchmark
              ( file : in file_type; deg,nbruns,inc : in integer32;
                c : in QuadDobl_Speelpenning_Convolutions.Circuits;
                x : in QuadDobl_Complex_VecVecs.VecVec;
                mxe : in Standard_Integer_Vectors.Vector;
                verbose : in boolean := true ) is

    use QuadDobl_Speelpenning_Convolutions;

    dim : constant integer32 := x'last;
    pwt : constant Link_to_VecVecVec := Allocate(mxe,deg);
    yd : constant QuadDobl_Complex_VecVecs.VecVec(1..dim+1)
       := Allocate_Coefficients(dim+1,deg);
    vy1 : constant QuadDobl_Complex_VecVecs.VecVec(0..deg)
        := Linearized_Allocation(dim,deg);
    vy2 : constant QuadDobl_Complex_VecVecs.VecVec(0..deg)
        := Linearized_Allocation(dim,deg);
    vm1 : constant QuadDobl_Complex_VecMats.VecMat(0..deg)
        := Allocate_Coefficients(dim,dim,deg);
    vm2 : constant QuadDobl_Complex_VecMats.VecMat(0..deg)
        := Allocate_Coefficients(dim,dim,deg);
    multstart,multstop,seristart,seristop : Ada.Calendar.Time;
    seri_elapsed,mult_elapsed,speedup,efficiency : Duration;
    nbt : integer32 := 2;

    use Ada.Calendar;

  begin
    put_line(file,"quad double precision"); flush(file);
    if verbose
     then put_line("Running on one task ...");
    end if;
    seristart := Ada.Calendar.Clock;
    Compute(pwt,mxe,x);
    EvalDiff(c,x,pwt,yd,vy1,vm1);
    seristop := Ada.Calendar.Clock;
    seri_elapsed := seristop - seristart;
    put(file,"  1 : ");
    duration_io.put(file,seri_elapsed,1,3); new_line(file); flush(file);
    if verbose then
      put_line("-> Elapsed time without multitasking : ");
      Time_Stamps.Write_Elapsed_Time(standard_output,seristart,seristop);
    end if;   
    for k in 1..nbruns loop
      if verbose then
        new_line;
        put("Running with "); put(nbt,1); put_line(" tasks ...");
      end if;
      multstart := Ada.Calendar.Clock;
      QuadDobl_Multitasked_EvalDiff(nbt,c,x,mxe,pwt,vy2,vm2,false);
      multstop := Ada.Calendar.Clock;
      mult_elapsed := multstop - multstart;
      if verbose then
        new_line;
        put_line("-> Elapsed time on multitasking : ");
        Time_Stamps.Write_Elapsed_Time(standard_output,multstart,multstop);
      end if;
      if seri_elapsed + 1.0 /= 1.0 then
        speedup := seri_elapsed/mult_elapsed;
        efficiency := speedup/duration(nbt);
        efficiency := duration(100)*efficiency;
        if verbose then
          put("The speedup : "); duration_io.put(speedup,1,3);
          put("  efficiency : "); duration_io.put(speedup,2,2); new_line;
        end if;
      end if;
      put(file,nbt,3);
      put(file," : "); duration_io.put(file,mult_elapsed,1,3);
      put(file," : "); duration_io.put(file,speedup,1,3);
      put(file," : "); duration_io.put(file,efficiency,2,2);
      new_line(file); flush(file);
      nbt := nbt + inc;
    end loop;
  end QuadDobl_Benchmark;

  procedure QuadDobl_Coefficient_Benchmark
              ( file : in file_type; deg,nbruns,inc : in integer32;
                c : in QuadDobl_Coefficient_Convolutions.Circuits;
                x : in QuadDobl_Complex_VecVecs.VecVec;
                mxe : in Standard_Integer_Vectors.Vector;
                verbose : in boolean := true ) is

    use Standard_Vector_Splitters;
    use QuadDobl_Coefficient_Convolutions;

    dim : constant integer32 := x'last;
   -- pwt : constant Link_to_VecVecVec := Allocate(mxe,deg);
    degdim : constant integer32 := 4*(deg+1)-1;
    xr : constant Standard_Floating_VecVecs.Link_to_VecVec
       := Standard_Vector_Splitters.Allocate_Floating_Coefficients(dim,degdim);
    xi : constant Standard_Floating_VecVecs.Link_to_VecVec
       := Standard_Vector_Splitters.Allocate_Floating_Coefficients(dim,degdim);
    rpwt : constant Standard_Floating_VecVecVecs.Link_to_VecVecVec
         := Standard_Coefficient_Convolutions.Allocate(mxe,degdim);
    ipwt : constant Standard_Floating_VecVecVecs.Link_to_VecVecVec
         := Standard_Coefficient_Convolutions.Allocate(mxe,degdim);
    ryd : constant Standard_Floating_VecVecs.Link_to_VecVec
        := Allocate_Floating_Coefficients(dim+1,degdim);
    iyd : constant Standard_Floating_VecVecs.Link_to_VecVec
        := Allocate_Floating_Coefficients(dim+1,degdim);
   -- yd : constant QuadDobl_Complex_VecVecs.VecVec(1..dim+1)
   --    := Allocate_Coefficients(dim+1,deg);
    vy1 : constant QuadDobl_Complex_VecVecs.VecVec(0..deg)
        := Linearized_Allocation(dim,deg);
    vy2 : constant QuadDobl_Complex_VecVecs.VecVec(0..deg)
        := Linearized_Allocation(dim,deg);
    vm1 : constant QuadDobl_Complex_VecMats.VecMat(0..deg)
        := Allocate_Coefficients(dim,dim,deg);
    vm2 : constant QuadDobl_Complex_VecMats.VecMat(0..deg)
        := Allocate_Coefficients(dim,dim,deg);
    multstart,multstop,seristart,seristop : Ada.Calendar.Time;
    seri_elapsed,mult_elapsed,speedup,efficiency : Duration;
    nbt : integer32 := 2;
    ur,vr,wr : constant Standard_Floating_Vectors.Vector(0..3)
             := (0..3 => 0.0);
    u : constant Standard_Floating_Vectors.Link_to_Vector
      := new Standard_Floating_Vectors.Vector'(ur);
    v : constant Standard_Floating_Vectors.Link_to_Vector
      := new Standard_Floating_Vectors.Vector'(vr);
    w : constant Standard_Floating_Vectors.Link_to_Vector
      := new Standard_Floating_Vectors.Vector'(wr);

    use Ada.Calendar;

  begin
    QuadDobl_Vector_Splitters.Complex_Parts(x,xr,xi);
    put_line(file,"quad double precision"); flush(file);
    if verbose
     then put_line("Running on one task ...");
    end if;
    seristart := Ada.Calendar.Clock;
    Compute(rpwt,ipwt,mxe,xr,xi,u,v,w);
    EvalDiff(c,xr.all,xi.all,rpwt,ipwt,ryd.all,iyd.all,vy1,vm1,u,v,w);
    seristop := Ada.Calendar.Clock;
    seri_elapsed := seristop - seristart;
    put(file,"  1 : ");
    duration_io.put(file,seri_elapsed,1,3); new_line(file); flush(file);
    if verbose then
      put_line("-> Elapsed time without multitasking : ");
      Time_Stamps.Write_Elapsed_Time(standard_output,seristart,seristop);
    end if;   
    for k in 1..nbruns loop
      if verbose then
        new_line;
        put("Running with "); put(nbt,1); put_line(" tasks ...");
      end if;
      multstart := Ada.Calendar.Clock;
      QuadDobl_Multitasked_EvalDiff
        (nbt,c,xr,xi,mxe,rpwt,ipwt,vy2,vm2,false,false);
      multstop := Ada.Calendar.Clock;
      mult_elapsed := multstop - multstart;
      if verbose then
        new_line;
        put_line("-> Elapsed time on multitasking : ");
        Time_Stamps.Write_Elapsed_Time(standard_output,multstart,multstop);
      end if;
      if seri_elapsed + 1.0 /= 1.0 then
        speedup := seri_elapsed/mult_elapsed;
        efficiency := speedup/duration(nbt);
        efficiency := duration(100)*efficiency;
        if verbose then
          put("The speedup : "); duration_io.put(speedup,1,3);
          put("  efficiency : "); duration_io.put(speedup,2,2); new_line;
        end if;
      end if;
      put(file,nbt,3);
      put(file," : "); duration_io.put(file,mult_elapsed,1,3);
      put(file," : "); duration_io.put(file,speedup,1,3);
      put(file," : "); duration_io.put(file,efficiency,2,2);
      new_line(file); flush(file);
      nbt := nbt + inc;
    end loop;
  end QuadDobl_Coefficient_Benchmark;

  procedure PentDobl_Benchmark
              ( file : in file_type; deg,nbruns,inc : in integer32;
                c : in PentDobl_Speelpenning_Convolutions.Circuits;
                x : in PentDobl_Complex_VecVecs.VecVec;
                mxe : in Standard_Integer_Vectors.Vector;
                verbose : in boolean := true ) is

    use PentDobl_Speelpenning_Convolutions;

    dim : constant integer32 := x'last;
    pwt : constant Link_to_VecVecVec := Allocate(mxe,deg);
    yd : constant PentDobl_Complex_VecVecs.VecVec(1..dim+1)
       := Allocate_Coefficients(dim+1,deg);
    vy1 : constant PentDobl_Complex_VecVecs.VecVec(0..deg)
        := Linearized_Allocation(dim,deg);
    vy2 : constant PentDobl_Complex_VecVecs.VecVec(0..deg)
        := Linearized_Allocation(dim,deg);
    vm1 : constant PentDobl_Complex_VecMats.VecMat(0..deg)
        := Allocate_Coefficients(dim,dim,deg);
    vm2 : constant PentDobl_Complex_VecMats.VecMat(0..deg)
        := Allocate_Coefficients(dim,dim,deg);
    multstart,multstop,seristart,seristop : Ada.Calendar.Time;
    seri_elapsed,mult_elapsed,speedup,efficiency : Duration;
    nbt : integer32 := 2;

    use Ada.Calendar;

  begin
    put_line(file,"penta double precision"); flush(file);
    if verbose
     then put_line("Running on one task ...");
    end if;
    seristart := Ada.Calendar.Clock;
    Compute(pwt,mxe,x);
    EvalDiff(c,x,pwt,yd,vy1,vm1);
    seristop := Ada.Calendar.Clock;
    seri_elapsed := seristop - seristart;
    put(file,"  1 : ");
    duration_io.put(file,seri_elapsed,1,3); new_line(file); flush(file);
    if verbose then
      put_line("-> Elapsed time without multitasking : ");
      Time_Stamps.Write_Elapsed_Time(standard_output,seristart,seristop);
    end if;   
    for k in 1..nbruns loop
      if verbose then
        new_line;
        put("Running with "); put(nbt,1); put_line(" tasks ...");
      end if;
      multstart := Ada.Calendar.Clock;
      PentDobl_Multitasked_EvalDiff(nbt,c,x,mxe,pwt,vy2,vm2,false);
      multstop := Ada.Calendar.Clock;
      mult_elapsed := multstop - multstart;
      if verbose then
        new_line;
        put_line("-> Elapsed time on multitasking : ");
        Time_Stamps.Write_Elapsed_Time(standard_output,multstart,multstop);
      end if;
      if seri_elapsed + 1.0 /= 1.0 then
        speedup := seri_elapsed/mult_elapsed;
        efficiency := speedup/duration(nbt);
        efficiency := duration(100)*efficiency;
        if verbose then
          put("The speedup : "); duration_io.put(speedup,1,3);
          put("  efficiency : "); duration_io.put(speedup,2,2); new_line;
        end if;
      end if;
      put(file,nbt,3);
      put(file," : "); duration_io.put(file,mult_elapsed,1,3);
      put(file," : "); duration_io.put(file,speedup,1,3);
      put(file," : "); duration_io.put(file,efficiency,2,2);
      new_line(file); flush(file);
      nbt := nbt + inc;
    end loop;
  end PentDobl_Benchmark;

  procedure OctoDobl_Benchmark
              ( file : in file_type; deg,nbruns,inc : in integer32;
                c : in OctoDobl_Speelpenning_Convolutions.Circuits;
                x : in OctoDobl_Complex_VecVecs.VecVec;
                mxe : in Standard_Integer_Vectors.Vector;
                verbose : in boolean := true ) is

    use OctoDobl_Speelpenning_Convolutions;

    dim : constant integer32 := x'last;
    pwt : constant Link_to_VecVecVec := Allocate(mxe,deg);
    yd : constant OctoDobl_Complex_VecVecs.VecVec(1..dim+1)
       := Allocate_Coefficients(dim+1,deg);
    vy1 : constant OctoDobl_Complex_VecVecs.VecVec(0..deg)
        := Linearized_Allocation(dim,deg);
    vy2 : constant OctoDobl_Complex_VecVecs.VecVec(0..deg)
        := Linearized_Allocation(dim,deg);
    vm1 : constant OctoDobl_Complex_VecMats.VecMat(0..deg)
        := Allocate_Coefficients(dim,dim,deg);
    vm2 : constant OctoDobl_Complex_VecMats.VecMat(0..deg)
        := Allocate_Coefficients(dim,dim,deg);
    multstart,multstop,seristart,seristop : Ada.Calendar.Time;
    seri_elapsed,mult_elapsed,speedup,efficiency : Duration;
    nbt : integer32 := 2;

    use Ada.Calendar;

  begin
    put_line(file,"octo double precision"); flush(file);
    if verbose
     then put_line("Running on one task ...");
    end if;
    seristart := Ada.Calendar.Clock;
    Compute(pwt,mxe,x);
    EvalDiff(c,x,pwt,yd,vy1,vm1);
    seristop := Ada.Calendar.Clock;
    seri_elapsed := seristop - seristart;
    put(file,"  1 : ");
    duration_io.put(file,seri_elapsed,1,3); new_line(file); flush(file);
    if verbose then
      put_line("-> Elapsed time without multitasking : ");
      Time_Stamps.Write_Elapsed_Time(standard_output,seristart,seristop);
    end if;   
    for k in 1..nbruns loop
      if verbose then
        new_line;
        put("Running with "); put(nbt,1); put_line(" tasks ...");
      end if;
      multstart := Ada.Calendar.Clock;
      OctoDobl_Multitasked_EvalDiff(nbt,c,x,mxe,pwt,vy2,vm2,false);
      multstop := Ada.Calendar.Clock;
      mult_elapsed := multstop - multstart;
      if verbose then
        new_line;
        put_line("-> Elapsed time on multitasking : ");
        Time_Stamps.Write_Elapsed_Time(standard_output,multstart,multstop);
      end if;
      if seri_elapsed + 1.0 /= 1.0 then
        speedup := seri_elapsed/mult_elapsed;
        efficiency := speedup/duration(nbt);
        efficiency := duration(100)*efficiency;
        if verbose then
          put("The speedup : "); duration_io.put(speedup,1,3);
          put("  efficiency : "); duration_io.put(speedup,2,2); new_line;
        end if;
      end if;
      put(file,nbt,3);
      put(file," : "); duration_io.put(file,mult_elapsed,1,3);
      put(file," : "); duration_io.put(file,speedup,1,3);
      put(file," : "); duration_io.put(file,efficiency,2,2);
      new_line(file); flush(file);
      nbt := nbt + inc;
    end loop;
  end OctoDobl_Benchmark;

  procedure DecaDobl_Benchmark
              ( file : in file_type; deg,nbruns,inc : in integer32;
                c : in DecaDobl_Speelpenning_Convolutions.Circuits;
                x : in DecaDobl_Complex_VecVecs.VecVec;
                mxe : in Standard_Integer_Vectors.Vector;
                verbose : in boolean := true ) is

    use DecaDobl_Speelpenning_Convolutions;

    dim : constant integer32 := x'last;
    pwt : constant Link_to_VecVecVec := Allocate(mxe,deg);
    yd : constant DecaDobl_Complex_VecVecs.VecVec(1..dim+1)
       := Allocate_Coefficients(dim+1,deg);
    vy1 : constant DecaDobl_Complex_VecVecs.VecVec(0..deg)
        := Linearized_Allocation(dim,deg);
    vy2 : constant DecaDobl_Complex_VecVecs.VecVec(0..deg)
        := Linearized_Allocation(dim,deg);
    vm1 : constant DecaDobl_Complex_VecMats.VecMat(0..deg)
        := Allocate_Coefficients(dim,dim,deg);
    vm2 : constant DecaDobl_Complex_VecMats.VecMat(0..deg)
        := Allocate_Coefficients(dim,dim,deg);
    multstart,multstop,seristart,seristop : Ada.Calendar.Time;
    seri_elapsed,mult_elapsed,speedup,efficiency : Duration;
    nbt : integer32 := 2;

    use Ada.Calendar;

  begin
    put_line(file,"deca double precision"); flush(file);
    if verbose
     then put_line("Running on one task ...");
    end if;
    seristart := Ada.Calendar.Clock;
    Compute(pwt,mxe,x);
    EvalDiff(c,x,pwt,yd,vy1,vm1);
    seristop := Ada.Calendar.Clock;
    seri_elapsed := seristop - seristart;
    put(file,"  1 : ");
    duration_io.put(file,seri_elapsed,1,3); new_line(file); flush(file);
    if verbose then
      put_line("-> Elapsed time without multitasking : ");
      Time_Stamps.Write_Elapsed_Time(standard_output,seristart,seristop);
    end if;   
    for k in 1..nbruns loop
      if verbose then
        new_line;
        put("Running with "); put(nbt,1); put_line(" tasks ...");
      end if;
      multstart := Ada.Calendar.Clock;
      DecaDobl_Multitasked_EvalDiff(nbt,c,x,mxe,pwt,vy2,vm2,false);
      multstop := Ada.Calendar.Clock;
      mult_elapsed := multstop - multstart;
      if verbose then
        new_line;
        put_line("-> Elapsed time on multitasking : ");
        Time_Stamps.Write_Elapsed_Time(standard_output,multstart,multstop);
      end if;
      if seri_elapsed + 1.0 /= 1.0 then
        speedup := seri_elapsed/mult_elapsed;
        efficiency := speedup/duration(nbt);
        efficiency := duration(100)*efficiency;
        if verbose then
          put("The speedup : "); duration_io.put(speedup,1,3);
          put("  efficiency : "); duration_io.put(speedup,2,2); new_line;
        end if;
      end if;
      put(file,nbt,3);
      put(file," : "); duration_io.put(file,mult_elapsed,1,3);
      put(file," : "); duration_io.put(file,speedup,1,3);
      put(file," : "); duration_io.put(file,efficiency,2,2);
      new_line(file); flush(file);
      nbt := nbt + inc;
    end loop;
  end DecaDobl_Benchmark;

  procedure Benchmark ( dim,deg,nbr,pwr : in integer32 ) is

    da_c : constant DecaDobl_Speelpenning_Convolutions.Circuits
         := DecaDobl_Random_Convolution_Circuits(dim,deg,nbr,pwr);
    od_c : constant OctoDobl_Speelpenning_Convolutions.Circuits
         := to_octo_double(da_c);
    pd_c : constant PentDobl_Speelpenning_Convolutions.Circuits
         := to_penta_double(da_c);
    qd_c : constant QuadDobl_Speelpenning_Convolutions.Circuits
         := to_quad_double(da_c);
    qd_cc : constant QuadDobl_Coefficient_Convolutions.Circuits
          := QuadDobl_Convolution_Splitters.Split(qd_c);
    td_c : constant TripDobl_Speelpenning_Convolutions.Circuits
         := to_triple_double(qd_c);
    dd_c : constant DoblDobl_Speelpenning_Convolutions.Circuits
         := to_double_double(qd_c);
    dd_cc : constant DoblDobl_Coefficient_Convolutions.Circuits
          := DoblDobl_Convolution_Splitters.Split(dd_c);
    d_c : constant Standard_Speelpenning_Convolutions.Circuits
        := to_double(qd_c);
    d_cc : constant Standard_Coefficient_Convolutions.Circuits
         := Standard_Convolution_Splitters.Split(d_c);
    da_x : constant DecaDobl_Complex_Series_Vectors.Vector(1..dim)
         := DecaDobl_Random_Series_Vectors.Random_Series_Vector(1,dim,deg);
    da_xcff : constant DecaDobl_Complex_VecVecs.VecVec(1..dim)
            := DecaDobl_Series_Coefficients(da_x);
    od_xcff : constant OctoDobl_Complex_VecVecs.VecVec(1..dim)
            := DecaDobl_Complex_Vectors_cv.to_octo_double(da_xcff);
    pd_xcff : constant PentDobl_Complex_VecVecs.VecVec(1..dim)
            := DecaDobl_Complex_Vectors_cv.to_penta_double(da_xcff);
    qd_xcff : constant QuadDobl_Complex_VecVecs.VecVec(1..dim)
            := DecaDobl_Complex_Vectors_cv.to_quad_double(da_xcff);
    td_xcff : constant TripDobl_Complex_VecVecs.VecVec(1..dim)
            := DecaDobl_Complex_Vectors_cv.to_triple_double(da_xcff);
    dd_xcff : constant DoblDobl_Complex_VecVecs.VecVec(1..dim)
            := DecaDobl_Complex_Vectors_cv.to_double_double(da_xcff);
    d_xcff : constant Standard_Complex_VecVecs.VecVec(1..dim)
           := DecaDobl_Complex_Vectors_cv.to_double(da_xcff);
    mxe : constant Standard_Integer_Vectors.Vector(1..dim)
        := QuadDobl_Speelpenning_Convolutions.Exponent_Maxima(qd_c,dim);
    nbruns,inc : integer32 := 0;
    file : file_type;
    rtp : character;

  begin
    new_line;
    put_line("MENU for benchmarked runs :");
    put_line("  1. complex convolutions only");
    put_line("  2. coefficient convolutions only");
    put_line("  3. both complex and coefficient convolutions");
    put("Type 1, 2, or 3 to select the type of runs : ");
    Ask_Alternative(rtp,"123");
    new_line;
    put("Give the number of multitasked runs : "); get(nbruns);
    put("Give the increment on the tasks : "); get(inc);
    skip_line;
    new_line;
    put_line("Reading the name of the output file ...");
    Read_Name_and_Create_File(file);
    new_line;
    put_line("See the output file for results ...");
    new_line;
    case rtp is
      when '1' =>
        Standard_Benchmark(file,deg,nbruns,inc,d_c,d_xcff,mxe,false);
        DoblDobl_Benchmark(file,deg,nbruns,inc,dd_c,dd_xcff,mxe,false);
        TripDobl_Benchmark(file,deg,nbruns,inc,td_c,td_xcff,mxe,false);
        QuadDobl_Benchmark(file,deg,nbruns,inc,qd_c,qd_xcff,mxe,false);
        PentDobl_Benchmark(file,deg,nbruns,inc,pd_c,pd_xcff,mxe,false);
        OctoDobl_Benchmark(file,deg,nbruns,inc,od_c,od_xcff,mxe,false);
        DecaDobl_Benchmark(file,deg,nbruns,inc,da_c,da_xcff,mxe,false);
      when '2' =>
        put_line(file,"runs on coefficient convolution circuits ...");
        Standard_Coefficient_Benchmark
          (file,deg,nbruns,inc,d_cc,d_xcff,mxe,false);
        DoblDobl_Coefficient_Benchmark
          (file,deg,nbruns,inc,dd_cc,dd_xcff,mxe,false);
        QuadDobl_Coefficient_Benchmark
          (file,deg,nbruns,inc,qd_cc,qd_xcff,mxe,false);
      when '3' =>
        Standard_Benchmark(file,deg,nbruns,inc,d_c,d_xcff,mxe,false);
        put_line(file,"runs on coefficient convolution circuits ...");
        Standard_Coefficient_Benchmark
          (file,deg,nbruns,inc,d_cc,d_xcff,mxe,false);
        DoblDobl_Benchmark(file,deg,nbruns,inc,dd_c,dd_xcff,mxe,false);
        put_line(file,"runs on coefficient convolution circuits ...");
        DoblDobl_Coefficient_Benchmark
          (file,deg,nbruns,inc,dd_cc,dd_xcff,mxe,false);
        QuadDobl_Benchmark(file,deg,nbruns,inc,qd_c,qd_xcff,mxe,false);
        put_line(file,"runs on coefficient convolution circuits ...");
        QuadDobl_Coefficient_Benchmark
          (file,deg,nbruns,inc,qd_cc,qd_xcff,mxe,false);
        PentDobl_Benchmark(file,deg,nbruns,inc,pd_c,pd_xcff,mxe,false);
        OctoDobl_Benchmark(file,deg,nbruns,inc,od_c,od_xcff,mxe,false);
        DecaDobl_Benchmark(file,deg,nbruns,inc,da_c,da_xcff,mxe,false);
      when others => null;
    end case;
  end Benchmark;

  function Prompt_for_Precision return character is

    res : character;

  begin
    new_line;
    put_line("MENU for the working precision :");
    put_line("  1. double precision");
    put_line("  2. double double precision");
    put_line("  3. triple double precision");
    put_line("  4. quad double precision");
    put_line("  5. penta double precision");
    put_line("  6. octo double precision");
    put("Type 1, 2, 3, 4, 5, 6, or 7 to select the precision : ");
    Ask_Alternative(res,"1234567");
    return res;
  end Prompt_for_Precision;

  procedure Standard_Test_Problem is

    lp : Standard_Complex_Poly_Systems.Link_to_Poly_Sys;
    deg,nbt : integer32 := 0;
    answer : character;
    cffcrc,static,output : boolean;
    dim : natural32;

  begin
    new_line;
    put_line("Reading a polynomial system ...");
    get(lp);
    dim := Standard_Complex_Polynomials.Number_of_Unknowns(lp(lp'first));
    new_line;
    put("Give the degree of the series : "); get(deg);
    put("Give the number of tasks : "); get(nbt);
    put("Output during multitasking ? (y/n) "); Ask_Yes_or_No(answer);
    output := (answer = 'y');
    put("Run on circuits with splitted coefficient vectors ? (y/n) ");
    Ask_Yes_or_No(answer); cffcrc := (answer = 'y');
    if answer = 'y' then
      put("Static load balancing ? (y/n) "); Ask_Yes_or_No(answer);
      static := (answer = 'y');
    end if;
    declare
      c : constant Standard_Speelpenning_Convolutions.Circuits(lp'range)
        := Make_Convolution_Circuits(lp.all,natural32(deg));
      cc : constant Standard_Coefficient_Convolutions.Circuits
         := Standard_Convolution_Splitters.Split(c);
    begin
      if cffcrc then
        Standard_Coefficient_Test(cc,integer32(dim),deg,nbt,static,output);
      else
        Standard_Test(c,integer32(dim),deg,nbt,output);
      end if;
    end;
  end Standard_Test_Problem;

  procedure DoblDobl_Test_Problem is

    use DoblDobl_Speelpenning_Convolutions;

    lp : DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
    deg,nbt : integer32 := 0;
    dim : natural32;
    answer : character;
    output,cffcrc,static : boolean;

  begin
    new_line;
    put_line("Reading a polynomial system ...");
    get(lp);
    dim := DoblDobl_Complex_Polynomials.Number_of_Unknowns(lp(lp'first));
    new_line;
    put("Give the degree of the series : "); get(deg);
    put("Give the number of tasks : "); get(nbt);
    put("Output during multitasking ? (y/n) "); Ask_Yes_or_No(answer);
    output := (answer = 'y');
    put("Run on circuits with splitted coefficient vectors ? (y/n) ");
    Ask_Yes_or_No(answer); cffcrc := (answer = 'y');
    if answer = 'y' then
      put("Static load balancing ? (y/n) "); Ask_Yes_or_No(answer);
      static := (answer = 'y');
    end if;
    declare
      c : constant DoblDobl_Speelpenning_Convolutions.Circuits(lp'range)
        := Make_Convolution_Circuits(lp.all,natural32(deg));
      cc : constant DoblDobl_Coefficient_Convolutions.Circuits
         := DoblDobl_Convolution_Splitters.Split(c);
    begin
      if cffcrc then
        DoblDobl_Coefficient_Test(cc,integer32(dim),deg,nbt,static,output);
      else
        DoblDobl_Test(c,integer32(dim),deg,nbt,output);
      end if;
    end;
  end DoblDobl_Test_Problem;

  procedure TripDobl_Test_Problem is

    use TripDobl_Speelpenning_Convolutions;

    lp : TripDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
    deg,nbt : integer32 := 0;
    dim : natural32;
    answer : character;
    output : boolean;

  begin
    new_line;
    put_line("Reading a polynomial system ...");
    get(lp);
    dim := TripDobl_Complex_Polynomials.Number_of_Unknowns(lp(lp'first));
    new_line;
    put("Give the degree of the series : "); get(deg);
    put("Give the number of tasks : "); get(nbt);
    put("Output during multitasking ? (y/n) "); Ask_Yes_or_No(answer);
    output := (answer = 'y');
    declare
      c : constant TripDobl_Speelpenning_Convolutions.Circuits(lp'range)
        := Make_Convolution_Circuits(lp.all,natural32(deg));
    begin
      TripDobl_Test(c,integer32(dim),deg,nbt,output);
    end;
  end TripDobl_Test_Problem;

  procedure QuadDobl_Test_Problem is

    lp : QuadDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
    deg,nbt : integer32 := 0;
    dim : natural32;
    answer : character;
    output,cffcrc,static : boolean;

  begin
    new_line;
    put_line("Reading a polynomial system ...");
    get(lp);
    dim := QuadDobl_Complex_Polynomials.Number_of_Unknowns(lp(lp'first));
    new_line;
    put("Give the degree of the series : "); get(deg);
    put("Give the number of tasks : "); get(nbt);
    put("Output during multitasking ? (y/n) "); Ask_Yes_or_No(answer);
    output := (answer = 'y');
    put("Run on circuits with splitted coefficient vectors ? (y/n) ");
    Ask_Yes_or_No(answer); cffcrc := (answer = 'y');
    if answer = 'y' then
      put("Static load balancing ? (y/n) "); Ask_Yes_or_No(answer);
      static := (answer = 'y');
    end if;
    declare
      c : constant QuadDobl_Speelpenning_Convolutions.Circuits(lp'range)
        := Make_Convolution_Circuits(lp.all,natural32(deg));
      cc : constant QuadDobl_Coefficient_Convolutions.Circuits
         := QuadDobl_Convolution_Splitters.Split(c);
    begin
      if cffcrc then
        QuadDobl_Coefficient_Test(cc,integer32(dim),deg,nbt,static,output);
      else
        QuadDobl_Test(c,integer32(dim),deg,nbt,output);
      end if;
    end;
  end QuadDobl_Test_Problem;

  procedure PentDobl_Test_Problem is

    lp : PentDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
    deg,nbt : integer32 := 0;
    dim : natural32;
    answer : character;
    output : boolean;

  begin
    new_line;
    put_line("Reading a polynomial system ...");
    get(lp);
    dim := PentDobl_Complex_Polynomials.Number_of_Unknowns(lp(lp'first));
    new_line;
    put("Give the degree of the series : "); get(deg);
    put("Give the number of tasks : "); get(nbt);
    put("Output during multitasking ? (y/n) "); Ask_Yes_or_No(answer);
    output := (answer = 'y');
    declare
      c : constant PentDobl_Speelpenning_Convolutions.Circuits(lp'range)
        := Make_Convolution_Circuits(lp.all,natural32(deg));
    begin
      PentDobl_Test(c,integer32(dim),deg,nbt,output);
    end;
  end PentDobl_Test_Problem;

  procedure OctoDobl_Test_Problem is

    lp : OctoDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
    deg,nbt : integer32 := 0;
    dim : natural32;
    answer : character;
    output : boolean;

  begin
    new_line;
    put_line("Reading a polynomial system ...");
    get(lp);
    dim := OctoDobl_Complex_Polynomials.Number_of_Unknowns(lp(lp'first));
    new_line;
    put("Give the degree of the series : "); get(deg);
    put("Give the number of tasks : "); get(nbt);
    put("Output during multitasking ? (y/n) "); Ask_Yes_or_No(answer);
    output := (answer = 'y');
    declare
      c : constant OctoDobl_Speelpenning_Convolutions.Circuits(lp'range)
        := Make_Convolution_Circuits(lp.all,natural32(deg));
    begin
      OctoDobl_Test(c,integer32(dim),deg,nbt,output);
    end;
  end OctoDobl_Test_Problem;

  procedure DecaDobl_Test_Problem is

    lp : DecaDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
    deg,nbt : integer32 := 0;
    dim : natural32;
    answer : character;
    output : boolean;

  begin
    new_line;
    put_line("Reading a polynomial system ...");
    get(lp);
    dim := DecaDobl_Complex_Polynomials.Number_of_Unknowns(lp(lp'first));
    new_line;
    put("Give the degree of the series : "); get(deg);
    put("Give the number of tasks : "); get(nbt);
    put("Output during multitasking ? (y/n) "); Ask_Yes_or_No(answer);
    output := (answer = 'y');
    declare
      c : constant DecaDobl_Speelpenning_Convolutions.Circuits(lp'range)
        := Make_Convolution_Circuits(lp.all,natural32(deg));
    begin
      DecaDobl_Test(c,integer32(dim),deg,nbt,output);
    end;
  end DecaDobl_Test_Problem;

  procedure Main is

    dim,deg,nbr,pwr,nbt : integer32 := 0;
    answer,precision : character;

  begin
    new_line;
    put_line("Testing the multitasked evaluation and differentiation ...");
    new_line;
    put("Test on random problems ? (y/n) ");
    Ask_Yes_or_No(answer);
    if answer /= 'y' then
      precision := Prompt_for_Precision;
      case precision is
        when '1' => Standard_Test_Problem;
        when '2' => DoblDobl_Test_Problem;
        when '3' => TripDobl_Test_Problem;
        when '4' => QuadDobl_Test_Problem;
        when '5' => PentDobl_Test_Problem;
        when '6' => OctoDobl_Test_Problem;
        when '7' => DecaDobl_Test_Problem;
        when others => null;
      end case;
    else
      new_line;
      put("Give the dimension : "); get(dim);
      put("Give the degree of the series : "); get(deg);
      put("Give the number of monomials : "); get(nbr);
      put("Give the highest power of each variable : "); get(pwr);
      put("Benchmarking for all precisions ? (y/n) "); Ask_Yes_or_No(answer);
      if answer  = 'y' then
        Benchmark(dim,deg,nbr,pwr);
      else
        put("Give the number of tasks : "); get(nbt);
        put("Output during multitasking ? (y/n) "); Ask_Yes_or_No(answer);
        precision := Prompt_for_Precision;
        new_line;
        case precision is
          when '1' => Standard_Random_Test(dim,deg,nbr,pwr,nbt,answer = 'y');
          when '2' => DoblDobl_Random_Test(dim,deg,nbr,pwr,nbt,answer = 'y');
          when '3' => TripDobl_Random_Test(dim,deg,nbr,pwr,nbt,answer = 'y');
          when '4' => QuadDobl_Random_Test(dim,deg,nbr,pwr,nbt,answer = 'y');
          when '5' => PentDobl_Random_Test(dim,deg,nbr,pwr,nbt,answer = 'y');
          when '6' => OctoDobl_Random_Test(dim,deg,nbr,pwr,nbt,answer = 'y');
          when '7' => DecaDobl_Random_Test(dim,deg,nbr,pwr,nbt,answer = 'y');
          when others => null;
        end case;
      end if;
    end if;
  end Main;

end Test_mtAlgoDiff_Convolutions;
