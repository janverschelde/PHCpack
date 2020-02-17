with text_io;                            use text_io;
with duration_io;
with Ada.Calendar;
with Communications_with_User;           use Communications_with_User;
with Time_Stamps;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
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
with Standard_Speelpenning_Convolutions;
with DoblDobl_Speelpenning_Convolutions;
with QuadDobl_Speelpenning_Convolutions;
with Random_Convolution_Circuits;        use Random_Convolution_Circuits;
with Shift_Convolution_Circuits;         use Shift_Convolution_Circuits;
with Multitasked_Shift_Circuits;         use Multitasked_Shift_Circuits;

procedure ts_mtshiftcnv is

-- DESCRIPTION :
--   Development of multitasked code to shift coefficients of circuits.

  procedure Standard_Test
              ( dim,deg,nbr,pwr : in integer32; nbt : in out integer32;
                output : in boolean := true ) is

  -- DESCRIPTION :
  --   Generates a random convolution circuit and applies multitasking
  --   to shift the coefficients in double precision.
  --   Compares with the one task run if output.

  -- ON ENTRY :
  --   dim      dimension of the exponent vectors;
  --   deg      degree of the power series;
  --   nbr      number of products;
  --   pwr      largest power of the variables;
  --   nbt      the number of tasks;
  --   output   flag for output during multitasking.

    use Standard_Complex_Numbers,Standard_Speelpenning_Convolutions;

    c : constant Circuits
      := Standard_Random_Convolution_Circuits(dim,deg,nbr,pwr);
    cwork : Circuits(c'range);
    multstart,multstop,seristart,seristop : Ada.Calendar.Time;
    seri_elapsed,mult_elapsed,speedup : Duration;
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
        put("The speedup : ");
        duration_io.put(speedup,1,3); new_line;
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

  -- DESCRIPTION :
  --   Generates a random convolution circuit and applies multitasking
  --   to shift the coefficients in double double precision.
  --   Compares with the one task run if output.

  -- ON ENTRY :
  --   dim      dimension of the exponent vectors;
  --   deg      degree of the power series;
  --   nbr      number of products;
  --   pwr      largest power of the variables;
  --   nbt      the number of tasks;
  --   output   flag for output during multitasking.

    use DoblDobl_Complex_Numbers,DoblDobl_Speelpenning_Convolutions;

    c : constant Circuits
      := DoblDobl_Random_Convolution_Circuits(dim,deg,nbr,pwr);
    cwork : Circuits(c'range);
    multstart,multstop,seristart,seristop : Ada.Calendar.Time;
    seri_elapsed,mult_elapsed,speedup : Duration;
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
        put("The speedup : ");
        duration_io.put(speedup,1,3); new_line;
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

  -- DESCRIPTION :
  --   Generates a random convolution circuit and applies multitasking
  --   to shift the coefficients in double double precision.
  --   Compares with the one task run if output.

  -- ON ENTRY :
  --   dim      dimension of the exponent vectors;
  --   deg      degree of the power series;
  --   nbr      number of products;
  --   pwr      largest power of the variables;
  --   nbt      the number of tasks;
  --   output   flag for output during multitasking.

    use QuadDobl_Complex_Numbers,QuadDobl_Speelpenning_Convolutions;

    c : constant Circuits
      := QuadDobl_Random_Convolution_Circuits(dim,deg,nbr,pwr);
    cwork : Circuits(c'range);
    multstart,multstop,seristart,seristop : Ada.Calendar.Time;
    seri_elapsed,mult_elapsed,speedup : Duration;
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
        put("The speedup : ");
        duration_io.put(speedup,1,3); new_line;
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

  procedure Main is

  -- DESCRIPTION :
  --   Prompts the user for the parameters of the test
  --   and then launches the test.

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
  end Main;

begin
  Main;
end ts_mtshiftcnv;
