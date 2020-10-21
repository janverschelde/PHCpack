with duration_io;
with Ada.Calendar;
with Communications_with_User;           use Communications_with_User;
with Time_Stamps;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Double_Double_Numbers;              use Double_Double_Numbers;
with Quad_Double_Numbers;                use Quad_Double_Numbers;
with Standard_Complex_Numbers;
with DoblDobl_Complex_Numbers;
with QuadDobl_Complex_Numbers;
with Standard_Complex_Vectors;
with Standard_Complex_Vectors_io;        use Standard_Complex_Vectors_io;
with DoblDobl_Complex_Vectors;
with DoblDobl_Complex_Vectors_io;        use DoblDobl_Complex_Vectors_io;
with QuadDobl_Complex_Vectors;
with QuadDobl_Complex_Vectors_io;        use QuadDobl_Complex_Vectors_io;
with System_Convolution_Circuits;        use System_Convolution_Circuits;
with Random_Convolution_Circuits;        use Random_Convolution_Circuits;
with Shift_Convolution_Circuits;         use Shift_Convolution_Circuits;
with Multitasked_Shift_Circuits;         use Multitasked_Shift_Circuits;

package body Test_mtShift_Convolutions is

  procedure Standard_Test
              ( dim,deg,nbr,pwr : in integer32; nbt : in out integer32;
                output : in boolean := true ) is

    use Standard_Complex_Numbers;
    use Standard_Speelpenning_Convolutions;

    c : constant Circuits
      := Standard_Random_Convolution_Circuits(dim,deg,nbr,pwr);
    cwork : Circuits(c'range);
    multstart,multstop,seristart,seristop : Ada.Calendar.Time;
    seri_elapsed,mult_elapsed,speedup,efficiency : Duration;
    ans : character;
    dt : constant double_float := 0.1;
    zero : constant Complex_Number := Create(0.0);
    wrk : constant Standard_Complex_Vectors.Link_to_Vector
        := new Standard_Complex_Vectors.Vector'(0..deg => zero);
    cff1,cff2 : Standard_Complex_Vectors.Vector(0..deg);

    use Ada.Calendar; -- for the difference operator on Duration

  begin
    put_line("Running on one task ...");
    Copy(c,cwork);
    seristart := Ada.Calendar.Clock;
    Shift(cwork,wrk,dt);
    seristop := Ada.Calendar.Clock;
    seri_elapsed := seristop - seristart;
    put_line("-> Elapsed time without multitasking : ");
    Time_Stamps.Write_Elapsed_Time(standard_output,seristart,seristop);
    if output then
      Standard_Complex_Vectors.Copy(cwork(1).cff(1).all,cff1);
      put_line("The first shifted coefficients :"); put_line(cff1);
    end if;
    loop
      new_line;
      put("Running with "); put(nbt,1); put_line(" tasks ...");
      Copy(c,cwork);
      multstart := Ada.Calendar.Clock;
      Standard_Multitasking(nbt,deg,cwork,dt,output);
      multstop := Ada.Calendar.Clock;
      mult_elapsed := multstop - multstart;
      Time_Stamps.Write_Elapsed_Time(standard_output,multstart,multstop);
      if seri_elapsed + 1.0 /= 1.0 then
        speedup := seri_elapsed/mult_elapsed;
        put("The speedup : "); duration_io.put(speedup,1,3);
        efficiency := speedup/duration(nbt);
        efficiency := duration(100)*efficiency;
        put("  efficiency : "); duration_io.put(efficiency,2,2); new_line;
      end if;
      Standard_Complex_Vectors.Copy(cwork(1).cff(1).all,cff2);
      if output then
        Standard_Complex_Vectors.Copy(cwork(1).cff(1).all,cff2);
        put_line("The first shifted coefficients :"); put_line(cff2);
      end if;
      put("Continue with different number of tasks ? (y/n) ");
      Ask_Yes_or_No(ans);
      exit when (ans /= 'y');
      put("Give the number of tasks : "); get(nbt);
    end loop;
  end Standard_Test;

  procedure DoblDobl_Test
              ( dim,deg,nbr,pwr : in integer32; nbt : in out integer32;
                output : in boolean := true ) is

    use DoblDobl_Complex_Numbers;
    use DoblDobl_Speelpenning_Convolutions;

    c : constant Circuits
      := DoblDobl_Random_Convolution_Circuits(dim,deg,nbr,pwr);
    cwork : Circuits(c'range);
    multstart,multstop,seristart,seristop : Ada.Calendar.Time;
    seri_elapsed,mult_elapsed,speedup,efficiency : Duration;
    ans : character;
    dt : constant double_double := create(0.1);
    zero : constant Complex_Number := Create(integer(0));
    wrk : constant DoblDobl_Complex_Vectors.Link_to_Vector
        := new DoblDobl_Complex_Vectors.Vector'(0..deg => zero);
    cff1,cff2 : DoblDobl_Complex_Vectors.Vector(0..deg);

    use Ada.Calendar; -- for the difference operator on Duration

  begin
    put_line("Running on one task ...");
    Copy(c,cwork);
    seristart := Ada.Calendar.Clock;
    Shift(cwork,wrk,dt);
    seristop := Ada.Calendar.Clock;
    seri_elapsed := seristop - seristart;
    put_line("-> Elapsed time without multitasking : ");
    Time_Stamps.Write_Elapsed_Time(standard_output,seristart,seristop);
    if output then
      DoblDobl_Complex_Vectors.Copy(cwork(1).cff(1).all,cff1);
      put_line("The first shifted coefficients :"); put_line(cff1);
    end if;
    loop
      new_line;
      put("Running with "); put(nbt,1); put_line(" tasks ...");
      Copy(c,cwork);
      multstart := Ada.Calendar.Clock;
      DoblDobl_Multitasking(nbt,deg,cwork,dt,output);
      multstop := Ada.Calendar.Clock;
      mult_elapsed := multstop - multstart;
      Time_Stamps.Write_Elapsed_Time(standard_output,multstart,multstop);
      if seri_elapsed + 1.0 /= 1.0 then
        speedup := seri_elapsed/mult_elapsed;
        put("The speedup : "); duration_io.put(speedup,1,3);
        efficiency := speedup/duration(nbt);
        efficiency := duration(100)*efficiency;
        put("  efficiency : "); duration_io.put(efficiency,2,2); new_line;
      end if;
      DoblDobl_Complex_Vectors.Copy(cwork(1).cff(1).all,cff2);
      if output then
        DoblDobl_Complex_Vectors.Copy(cwork(1).cff(1).all,cff2);
        put_line("The first shifted coefficients :"); put_line(cff2);
      end if;
      put("Continue with different number of tasks ? (y/n) ");
      Ask_Yes_or_No(ans);
      exit when (ans /= 'y');
      put("Give the number of tasks : "); get(nbt);
    end loop;
  end DoblDobl_Test;

  procedure QuadDobl_Test
              ( dim,deg,nbr,pwr : in integer32; nbt : in out integer32;
                output : in boolean := true ) is

    use QuadDobl_Complex_Numbers;
    use QuadDobl_Speelpenning_Convolutions;

    c : constant Circuits
      := QuadDobl_Random_Convolution_Circuits(dim,deg,nbr,pwr);
    cwork : Circuits(c'range);
    multstart,multstop,seristart,seristop : Ada.Calendar.Time;
    seri_elapsed,mult_elapsed,speedup,efficiency : Duration;
    ans : character;
    dt : constant quad_double := create(0.1);
    zero : constant Complex_Number := Create(integer(0));
    wrk : constant QuadDobl_Complex_Vectors.Link_to_Vector
        := new QuadDobl_Complex_Vectors.Vector'(0..deg => zero);
    cff1,cff2 : QuadDobl_Complex_Vectors.Vector(0..deg);

    use Ada.Calendar; -- for the difference operator on Duration

  begin
    put_line("Running on one task ...");
    Copy(c,cwork);
    seristart := Ada.Calendar.Clock;
    Shift(cwork,wrk,dt);
    seristop := Ada.Calendar.Clock;
    seri_elapsed := seristop - seristart;
    put_line("-> Elapsed time without multitasking : ");
    Time_Stamps.Write_Elapsed_Time(standard_output,seristart,seristop);
    if output then
      QuadDobl_Complex_Vectors.Copy(cwork(1).cff(1).all,cff1);
      put_line("The first shifted coefficients :"); put_line(cff1);
    end if;
    loop
      new_line;
      put("Running with "); put(nbt,1); put_line(" tasks ...");
      Copy(c,cwork);
      multstart := Ada.Calendar.Clock;
      QuadDobl_Multitasking(nbt,deg,cwork,dt,output);
      multstop := Ada.Calendar.Clock;
      mult_elapsed := multstop - multstart;
      Time_Stamps.Write_Elapsed_Time(standard_output,multstart,multstop);
      if seri_elapsed + 1.0 /= 1.0 then
        speedup := seri_elapsed/mult_elapsed;
        put("The speedup : "); duration_io.put(speedup,1,3); new_line;
        efficiency := speedup/duration(nbt);
        efficiency := duration(100)*efficiency;
        put("  efficiency : "); duration_io.put(efficiency,2,2); new_line;
      end if;
      QuadDobl_Complex_Vectors.Copy(cwork(1).cff(1).all,cff2);
      if output then
        QuadDobl_Complex_Vectors.Copy(cwork(1).cff(1).all,cff2);
        put_line("The first shifted coefficients :"); put_line(cff2);
      end if;
      put("Continue with different number of tasks ? (y/n) ");
      Ask_Yes_or_No(ans);
      exit when (ans /= 'y');
      put("Give the number of tasks : "); get(nbt);
    end loop;
  end QuadDobl_Test;

  procedure Standard_Benchmark
              ( file : in file_type; deg,nbruns,inc : in integer32;
                c : in Standard_Speelpenning_Convolutions.Circuits;
                verbose : in boolean := true ) is

    use Standard_Complex_Numbers;
    use Standard_Speelpenning_Convolutions;

    cwork : Circuits(c'range);
    multstart,multstop,seristart,seristop : Ada.Calendar.Time;
    seri_elapsed,mult_elapsed,speedup,efficiency : Duration;
    nbt : integer32 := 2;
    t : constant double_float := 0.1;
    zero : constant Complex_Number := Create(0.0);
    wrk : constant Standard_Complex_Vectors.Link_to_Vector
        := new Standard_Complex_Vectors.Vector'(0..deg => zero);

    use Ada.Calendar;

  begin
    put_line(file,"double precision");
    if verbose
     then put_line("Running on one task ...");
    end if;
    Copy(c,cwork);
    seristart := Ada.Calendar.Clock;
    Shift(cwork,wrk,t);
    seristop := Ada.Calendar.Clock;
    seri_elapsed := seristop - seristart;
    put(file,"  1 : ");
    duration_io.put(file,seri_elapsed,1,3); new_line(file);
    if verbose then
      put_line("-> Elapsed time without multitasking : ");
      Time_Stamps.Write_Elapsed_Time(standard_output,seristart,seristop);
    end if;   
    for k in 1..nbruns loop
      if verbose then
        new_line;
        put("Running with "); put(nbt,1); put_line(" tasks ...");
      end if;
      Copy(c,cwork);
      multstart := Ada.Calendar.Clock;
      Standard_Multitasking(nbt,deg,cwork,t,false);
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
  end Standard_Benchmark;

  procedure DoblDobl_Benchmark
              ( file : in file_type; deg,nbruns,inc : in integer32;
                c : in DoblDobl_Speelpenning_Convolutions.Circuits;
                verbose : in boolean := true ) is

    use DoblDobl_Complex_Numbers;
    use DoblDobl_Speelpenning_Convolutions;

    cwork : Circuits(c'range);
    multstart,multstop,seristart,seristop : Ada.Calendar.Time;
    seri_elapsed,mult_elapsed,speedup,efficiency : Duration;
    nbt : integer32 := 2;
    t : constant double_double := create(0.1);
    zero : constant Complex_Number := Create(integer(0));
    wrk : constant DoblDobl_Complex_Vectors.Link_to_Vector
        := new DoblDobl_Complex_Vectors.Vector'(0..deg => zero);

    use Ada.Calendar;

  begin
    put_line(file,"double double precision");
    if verbose
     then put_line("Running on one task ...");
    end if;
    Copy(c,cwork);
    seristart := Ada.Calendar.Clock;
    Shift(cwork,wrk,t);
    seristop := Ada.Calendar.Clock;
    seri_elapsed := seristop - seristart;
    put(file,"  1 : ");
    duration_io.put(file,seri_elapsed,1,3); new_line(file);
    if verbose then
      put_line("-> Elapsed time without multitasking : ");
      Time_Stamps.Write_Elapsed_Time(standard_output,seristart,seristop);
    end if;   
    for k in 1..nbruns loop
      if verbose then
        new_line;
        put("Running with "); put(nbt,1); put_line(" tasks ...");
      end if;
      Copy(c,cwork);
      multstart := Ada.Calendar.Clock;
      DoblDobl_Multitasking(nbt,deg,cwork,t,false);
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
      put(file," : "); duration_io.put(file,efficiency,2,2);
      new_line(file); flush(file);
      nbt := nbt + inc;
    end loop;
  end DoblDobl_Benchmark;

  procedure QuadDobl_Benchmark
              ( file : in file_type; deg,nbruns,inc : in integer32;
                c : in QuadDobl_Speelpenning_Convolutions.Circuits;
                verbose : in boolean := true ) is

    use QuadDobl_Complex_Numbers;
    use QuadDobl_Speelpenning_Convolutions;

    cwork : Circuits(c'range);
    multstart,multstop,seristart,seristop : Ada.Calendar.Time;
    seri_elapsed,mult_elapsed,speedup,efficiency : Duration;
    nbt : integer32 := 2;
    t : constant quad_double := create(0.1);
    zero : constant Complex_Number := Create(integer(0));
    wrk : constant QuadDobl_Complex_Vectors.Link_to_Vector
        := new QuadDobl_Complex_Vectors.Vector'(0..deg => zero);

    use Ada.Calendar;

  begin
    put_line(file,"quad double precision");
    if verbose
     then put_line("Running on one task ...");
    end if;
    Copy(c,cwork);
    seristart := Ada.Calendar.Clock;
    Shift(cwork,wrk,t);
    seristop := Ada.Calendar.Clock;
    seri_elapsed := seristop - seristart;
    put(file,"  1 : ");
    duration_io.put(file,seri_elapsed,1,3); new_line(file);
    if verbose then
      put_line("-> Elapsed time without multitasking : ");
      Time_Stamps.Write_Elapsed_Time(standard_output,seristart,seristop);
    end if;   
    for k in 1..nbruns loop
      if verbose then
        new_line;
        put("Running with "); put(nbt,1); put_line(" tasks ...");
      end if;
      Copy(c,cwork);
      multstart := Ada.Calendar.Clock;
      QuadDobl_Multitasking(nbt,deg,cwork,t,false);
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
      put(file," : "); duration_io.put(file,efficiency,2,2);
      new_line(file); flush(file);
      nbt := nbt + inc;
    end loop;
  end QuadDobl_Benchmark;

  procedure Benchmark ( dim,deg,nbr,pwr : in integer32 ) is

    qd_c : constant QuadDobl_Speelpenning_Convolutions.Circuits
         := QuadDobl_Random_Convolution_Circuits(dim,deg,nbr,pwr);
    dd_c : constant DoblDobl_Speelpenning_Convolutions.Circuits
         := to_double_double(qd_c);
    d_c : constant Standard_Speelpenning_Convolutions.Circuits
        := to_double(qd_c);
    nbruns,inc : integer32 := 0;
    file : file_type;

  begin
    put("Give the number of multitasked runs : "); get(nbruns);
    put("Give the increment on the tasks : "); get(inc);
    skip_line;
    new_line;
    put_line("Reading the name of the output file ...");
    Read_Name_and_Create_File(file);
    new_line;
    put_line("See the output file for results ...");
    new_line;
    put(file,"dimension : "); put(file,dim,1);
    put(file,"  degree : "); put(file,deg,1); new_line(file);
    Standard_Benchmark(file,deg,nbruns,inc,d_c,false);
    DoblDobl_Benchmark(file,deg,nbruns,inc,dd_c,false);
    QuadDobl_Benchmark(file,deg,nbruns,inc,qd_c,false);
  end Benchmark;

  procedure Main is

    dim,deg,nbr,pwr,nbt : integer32 := 0;
    answer,precision : character;

  begin
    new_line;
    put_line("Testing multitasked shifting of coefficients ...");
    new_line;
    put("Give the dimension : "); get(dim);
    put("Give the degree of the series : "); get(deg);
    put("Give the number of monomials : "); get(nbr);
    put("Give the highest power of each variable : "); get(pwr);
    new_line;
    put("Benchmark for all precisions ? (y/n) "); Ask_Yes_or_No(answer);
    if answer = 'y' then
      Benchmark(dim,deg,nbr,pwr);
    else
      put("Give the number of tasks : "); get(nbt);
      put("Output during multitasking ? (y/n) "); Ask_Yes_or_No(answer);
      new_line;
      put_line("MENU for the working precision :");
      put_line("  0. double precision");
      put_line("  1. double double precision");
      put_line("  2. quad double precision");
      put("Type 0, 1, or 2 to select the working precision : ");
      Ask_Alternative(precision,"012");
      case precision is
        when '0' => Standard_Test(dim,deg,nbr,pwr,nbt,answer='y');
        when '1' => DoblDobl_Test(dim,deg,nbr,pwr,nbt,answer='y');
        when '2' => QuadDobl_Test(dim,deg,nbr,pwr,nbt,answer='y');
        when others => null;
      end case;
    end if;
  end Main;

end Test_mtShift_Convolutions;
