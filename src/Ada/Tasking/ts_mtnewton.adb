with text_io;                            use text_io;
with duration_io;
with Ada.Calendar;
with Communications_with_User;           use Communications_with_User;
with Time_Stamps;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;
with Double_Double_Numbers;              use Double_Double_Numbers;
with Double_Double_Numbers_io;           use Double_Double_Numbers_io;
with Quad_Double_Numbers;                use Quad_Double_Numbers;
with Quad_Double_Numbers_io;             use Quad_Double_Numbers_io;
with Standard_Integer_Vectors;
with Standard_Complex_Vectors;
with Standard_Complex_VecVecs;
with DoblDobl_Complex_Vectors;
with DoblDobl_Complex_VecVecs;
with QuadDobl_Complex_Vectors;
with QuadDobl_Complex_VecVecs;
with Standard_Complex_Poly_Systems;
with DoblDobl_Complex_Poly_Systems;
with QuadDobl_Complex_Poly_Systems;
with Standard_Complex_Solutions;
with Standard_System_and_Solutions_io;
with DoblDobl_Complex_Solutions;
with DoblDobl_System_and_Solutions_io;
with QuadDobl_Complex_Solutions;
with QuadDobl_System_and_Solutions_io;
with Standard_Speelpenning_Convolutions;
with DoblDobl_Speelpenning_Convolutions;
with QuadDobl_Speelpenning_Convolutions;
with System_Convolution_Circuits;        use System_Convolution_Circuits;
with Homotopy_Convolution_Circuits;      use Homotopy_Convolution_Circuits;
with Random_Convolution_Circuits;        use Random_Convolution_Circuits;
with Newton_Convolutions;
with Multitasked_Series_Linearization;
with Multitasked_Newton_Convolutions;    use Multitasked_Newton_Convolutions;

procedure ts_mtnewton is

-- DESCRIPTION :
--   Tests the development of Newton's method on power series
--   with the reverse mode of algorithmic differentation
--   and linearization to solve the matrix series equations,
--   in double, double double, and quad double arithmetic,
--   with multitasking for shared memory parallel computers.

  procedure Standard_Run
              ( nbt,dim,maxit : in integer32;
                s : in Standard_Speelpenning_Convolutions.Link_to_System;
                scf : in Standard_Complex_VecVecs.VecVec;
                serelp,mltelp : in out Duration;
                output : in boolean ) is

  -- DESCRIPTION :
  --   Runs Newton's method with nbt tasks in double precision.

  -- ON ENTRY :
  --   nbt      the number of tasks;
  --   dim      number of equations and variables;
  --   maxit    maximum number of iterations;
  --   s        system of convolution circuits;
  --   scf      its leading coefficients are solution vector;
  --   serelp   the previous elapsed wall clock time of a serial run;
  --   mltelp   the previous elapsed wall clock time of a multitasked run;
  --   output   if true, then the run is verbose, else silent;

  -- ON RETURN :
  --   serelp   updated elapsed wall clock time of a serial run,
  --            if nbt = 1 and the user did not want multitasking;
  --   mltelp   updated elapsed wall clock time of a multitasked run,
  --            if nbt > 1.

    wks : Standard_Complex_VecVecs.VecVec(1..nbt)
        := Multitasked_Series_Linearization.Allocate_Work_Space(nbt,dim);
    ipvt : Standard_Integer_Vectors.Vector(1..dim);
    info,nbrit : integer32 := 0;
    fail : boolean;
    tol : constant double_float := 1.0E-12;
    absdx : double_float;
    seristart,seristop,multstart,multstop : Ada.Calendar.Time;
    speedup : Duration := 0.0;

    use Ada.Calendar; -- for the difference operation on Duration

  begin
    new_line;
    put("Running with "); put(nbt,1); put_line(" tasks ...");
    if nbt = 1
     then seristart := Ada.Calendar.Clock;
     else multstart := Ada.Calendar.Clock;
    end if;
    Multitasked_LU_Newton_Steps
      (standard_output,nbt,s,scf,maxit,nbrit,tol,absdx,fail,
       info,ipvt,wks,output);
    if nbt = 1 then
      seristop := Ada.Calendar.Clock;
      serelp := seristop - seristart;
    else
      multstop := Ada.Calendar.Clock;
      mltelp := multstop - multstart;
    end if;
    put("#steps : "); put(nbrit,1); put("  absdx :"); put(absdx,3);
    if fail
     then put("  failed to reach tolerance");
     else put("  succeeded to reach tolerance");
    end if;
    put(tol,3); new_line;
    if nbt = 1 then
      Time_Stamps.Write_Elapsed_Time(standard_output,seristart,seristop);
    else
      Time_Stamps.Write_Elapsed_Time(standard_output,multstart,multstop);
      if serelp + 1.0 /= 1.0 then
        speedup := serelp/mltelp;
        put("The speedup : ");
        duration_io.put(speedup,1,3); new_line;
      end if;
    end if;
    Standard_Complex_VecVecs.Clear(wks);
  end Standard_Run;

  procedure DoblDobl_Run
              ( nbt,dim,maxit : in integer32;
                s : in DoblDobl_Speelpenning_Convolutions.Link_to_System;
                scf : in DoblDobl_Complex_VecVecs.VecVec;
                serelp,mltelp : in out Duration;
                output : in boolean ) is

  -- DESCRIPTION :
  --   Runs Newton's method with nbt tasks in double precision.

  -- ON ENTRY :
  --   nbt      the number of tasks;
  --   dim      number of equations and variables;
  --   maxit    maximum number of iterations;
  --   s        system of convolution circuits;
  --   scf      its leading coefficients are the solution vector;
  --   serelp   the previous elapsed wall clock time of a serial run;
  --   mltelp   the previous elapsed wall clock time of a multitasked run;
  --   output   if true, then the run is verbose, else silent;

  -- ON RETURN :
  --   serelp   updated elapsed wall clock time of a serial run,
  --            if nbt = 1 and the user did not want multitasking;
  --   mltelp   updated elapsed wall clock time of a multitasked run,
  --            if nbt > 1.

    wks : DoblDobl_Complex_VecVecs.VecVec(1..nbt)
        := Multitasked_Series_Linearization.Allocate_Work_Space(nbt,dim);
    ipvt : Standard_Integer_Vectors.Vector(1..dim);
    info,nbrit : integer32 := 0;
    fail : boolean;
    tol : constant double_double := create(1.0E-24);
    absdx : double_double;
    seristart,seristop,multstart,multstop : Ada.Calendar.Time;
    speedup : Duration := 0.0;

    use Ada.Calendar; -- for the difference operation on Duration

  begin
    new_line;
    put("Running with "); put(nbt,1); put_line(" tasks ...");
    if nbt = 1
     then seristart := Ada.Calendar.Clock;
     else multstart := Ada.Calendar.Clock;
    end if;
    Multitasked_LU_Newton_Steps
      (standard_output,nbt,s,scf,maxit,nbrit,tol,absdx,fail,
       info,ipvt,wks,output);
    if nbt = 1 then
      seristop := Ada.Calendar.Clock;
      serelp := seristop - seristart;
    else
      multstop := Ada.Calendar.Clock;
      mltelp := multstop - multstart;
    end if;
    put("#steps : "); put(nbrit,1); put("  absdx :"); put(absdx,3);
    if fail
     then put("  failed to reach tolerance ");
     else put("  succeeded to reach tolerance ");
    end if;
    put(tol,3); new_line;
    if nbt = 1 then
      Time_Stamps.Write_Elapsed_Time(standard_output,seristart,seristop);
    else
      Time_Stamps.Write_Elapsed_Time(standard_output,multstart,multstop);
      if serelp + 1.0 /= 1.0 then
        speedup := serelp/mltelp;
        put("The speedup : ");
        duration_io.put(speedup,1,3); new_line;
      end if;
    end if;
    DoblDobl_Complex_VecVecs.Clear(wks);
  end DoblDobl_Run;

  procedure QuadDobl_Run
              ( nbt,dim,maxit : in integer32;
                s : in QuadDobl_Speelpenning_Convolutions.Link_to_System;
                scf : in QuadDobl_Complex_VecVecs.VecVec;
                serelp,mltelp : in out Duration;
                output : in boolean ) is

  -- DESCRIPTION :
  --   Runs Newton's method with nbt tasks in double precision.

  -- ON ENTRY :
  --   nbt      the number of tasks;
  --   dim      number of equations and variables;
  --   deg      degree of the truncated series;
  --   maxit    maximum number of iterations;
  --   s        system of convolution circuits;
  --   scf      its leading coefficients are the solution vector;
  --   serelp   the previous elapsed wall clock time of a serial run;
  --   mltelp   the previous elapsed wall clock time of a multitasked run;
  --   output   if true, then the run is verbose, else silent;

  -- ON RETURN :
  --   serelp   updated elapsed wall clock time of a serial run,
  --            if nbt = 1 and the user did not want multitasking;
  --   mltelp   updated elapsed wall clock time of a multitasked run,
  --            if nbt > 1.

    wks : QuadDobl_Complex_VecVecs.VecVec(1..nbt)
        := Multitasked_Series_Linearization.Allocate_Work_Space(nbt,dim);
    ipvt : Standard_Integer_Vectors.Vector(1..dim);
    info,nbrit : integer32 := 0;
    fail : boolean;
    tol : constant quad_double := create(1.0E-48);
    absdx : quad_double;
    seristart,seristop,multstart,multstop : Ada.Calendar.Time;
    speedup : Duration := 0.0;

    use Ada.Calendar; -- for the difference operation on Duration

  begin
    new_line;
    put("Running with "); put(nbt,1); put_line(" tasks ...");
    if nbt = 1
     then seristart := Ada.Calendar.Clock;
     else multstart := Ada.Calendar.Clock;
    end if;
    Multitasked_LU_Newton_Steps
      (standard_output,nbt,s,scf,maxit,nbrit,tol,absdx,fail,
       info,ipvt,wks,output);
    if nbt = 1 then
      seristop := Ada.Calendar.Clock;
      serelp := seristop - seristart;
    else
      multstop := Ada.Calendar.Clock;
      mltelp := multstop - multstart;
    end if;
    put("#steps : "); put(nbrit,1); put("  absdx :"); put(absdx,3);
    if fail
     then put("  failed to reach tolerance ");
     else put("  succeeded to reach tolerance ");
    end if;
    put(tol,3); new_line;
    if nbt = 1 then
      Time_Stamps.Write_Elapsed_Time(standard_output,seristart,seristop);
    else
      Time_Stamps.Write_Elapsed_Time(standard_output,multstart,multstop);
      if serelp + 1.0 /= 1.0 then
        speedup := serelp/mltelp;
        put("The speedup : ");
        duration_io.put(speedup,1,3); new_line;
      end if;
    end if;
    QuadDobl_Complex_VecVecs.Clear(wks);
  end QuadDobl_Run;

  procedure Standard_Run_Loop
              ( p : in Standard_Complex_Poly_Systems.Link_to_Poly_Sys;
	        sol : in Standard_Complex_Vectors.Vector;
	        deg : in integer32 ) is

  -- DESCRIPTION :
  --   Runs Newton's method in double precision on a solution sol
  --   of the system p, with power series of degree deg.

    use Standard_Speelpenning_Convolutions;

    c : constant Circuits(p'range)
      := Make_Convolution_Circuits(p.all,natural32(deg));
    s : constant Link_to_System := Create(c,p'last,deg);
    dim : constant integer32 := sol'last;
    scf : Standard_Complex_VecVecs.VecVec(1..dim);
    maxit,nbt : integer32 := 0;
    ans : character;
    otp : boolean;
    seri_elapsed,mult_elapsed : Duration := 0.0;

  begin
    Add_Parameter_to_Constant(s);
    new_line;
    put("Give the maximum number of iterations : "); get(maxit);
    loop
      put("Give the number of tasks (0 to exit) : "); get(nbt);
      exit when (nbt = 0);
      put("Output during multitasking ? (y/n) "); Ask_Yes_or_No(ans);
      otp := (ans = 'y');
      scf := Newton_Convolutions.Series_Coefficients(sol,deg);
      Standard_Run(nbt,dim,maxit,s,scf,seri_elapsed,mult_elapsed,otp);
      Standard_Complex_VecVecs.Clear(scf);
    end loop;
  end Standard_Run_Loop;

  procedure DoblDobl_Run_Loop
              ( p : in DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
	        sol : in DoblDobl_Complex_Vectors.Vector;
	        deg : in integer32 ) is

  -- DESCRIPTION :
  --   Runs Newton's method in double double precision on a solution sol
  --   of the system p, with power series of degree deg.

    use DoblDobl_Speelpenning_Convolutions;

    c : constant Circuits(p'range)
      := Make_Convolution_Circuits(p.all,natural32(deg));
    s : constant Link_to_System := Create(c,p'last,deg);
    dim : constant integer32 := sol'last;
    scf : DoblDobl_Complex_VecVecs.VecVec(1..dim);
    maxit,nbt : integer32 := 0;
    ans : character;
    otp : boolean;
    seri_elapsed,mult_elapsed : Duration := 0.0;

  begin
    Add_Parameter_to_Constant(s);
    new_line;
    put("Give the maximum number of iterations : "); get(maxit);
    loop
      put("Give the number of tasks (0 to exit) : "); get(nbt);
      exit when (nbt = 0);
      put("Output during multitasking ? (y/n) "); Ask_Yes_or_No(ans);
      otp := (ans = 'y');
      scf := Newton_Convolutions.Series_Coefficients(sol,deg);
      DoblDobl_Run(nbt,dim,maxit,s,scf,seri_elapsed,mult_elapsed,otp);
      DoblDobl_Complex_VecVecs.Clear(scf);
    end loop;
  end DoblDobl_Run_Loop;

  procedure QuadDobl_Run_Loop
              ( p : in QuadDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
	        sol : in QuadDobl_Complex_Vectors.Vector;
	        deg : in integer32 ) is

  -- DESCRIPTION :
  --   Runs Newton's method in quad double precision on a solution sol
  --   of the system p, with power series of degree deg.

    use QuadDobl_Speelpenning_Convolutions;

    c : constant Circuits(p'range)
      := Make_Convolution_Circuits(p.all,natural32(deg));
    s : constant Link_to_System := Create(c,p'last,deg);
    dim : constant integer32 := sol'last;
    scf : QuadDobl_Complex_VecVecs.VecVec(1..dim);
    maxit,nbt : integer32 := 0;
    ans : character;
    otp : boolean;
    seri_elapsed,mult_elapsed : Duration := 0.0;

  begin
    Add_Parameter_to_Constant(s);
    new_line;
    put("Give the maximum number of iterations : "); get(maxit);
    loop
      put("Give the number of tasks (0 to exit) : "); get(nbt);
      exit when (nbt = 0);
      put("Output during multitasking ? (y/n) "); Ask_Yes_or_No(ans);
      otp := (ans = 'y');
      scf := Newton_Convolutions.Series_Coefficients(sol,deg);
      QuadDobl_Run(nbt,dim,maxit,s,scf,seri_elapsed,mult_elapsed,otp);
      QuadDobl_Complex_VecVecs.Clear(scf);
    end loop;
  end QuadDobl_Run_Loop;

  procedure Standard_Test is

  -- DESCRIPTION :
  --   Prompts the user for a polynomial system with solutions
  --   and launches the test in double precision.

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

  -- DESCRIPTION :
  --   Prompts the user for a polynomial system with solutions
  --   and launches the test in double double precision.

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

  procedure QuadDobl_Test is

  -- DESCRIPTION :
  --   Prompts the user for a polynomial system with solutions
  --   and launches the test in quad double precision.

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

  procedure Standard_Random_Test
              ( dim,deg,nbr,pwr : in integer32 ) is

  -- DESCRIPTION :
  --   Generates a random Newton homotopy in double precision.

  -- ON ENTRY :
  --   dim      dimension of the exponent vectors;
  --   deg      degree of the power series;
  --   nbr      number of products;
  --   pwr      largest power of the variables.

    use Standard_Complex_VecVecs;
    use Standard_Speelpenning_Convolutions;

    s : Link_to_System;
    x : Link_to_VecVec;
    scf : VecVec(1..dim);
    maxit,nbt : integer32 := 0;
    ans : character;
    otp : boolean;
    seri_elapsed,mult_elapsed : Duration := 0.0;

  begin
    Standard_Random_Newton_Homotopy(dim,deg,nbr,pwr,s,x);
    new_line;
    put("Give the maximum number of iterations : "); get(maxit);
    loop
      put("Give the number of tasks (0 to exit) : "); get(nbt);
      exit when (nbt = 0);
      put("Output during multitasking ? (y/n) "); Ask_Yes_or_No(ans);
      otp := (ans = 'y');
      Standard_Complex_VecVecs.Copy(x.all,scf);
      Standard_Run(nbt,dim,maxit,s,scf,seri_elapsed,mult_elapsed,otp);
      Standard_Complex_VecVecs.Clear(scf);
    end loop;
  end Standard_Random_Test;

  procedure DoblDobl_Random_Test
              ( dim,deg,nbr,pwr : in integer32 ) is

  -- DESCRIPTION :
  --   Generates a random Newton homotopy in double double precision.

  -- ON ENTRY :
  --   dim      dimension of the exponent vectors;
  --   deg      degree of the power series;
  --   nbr      number of products;
  --   pwr      largest power of the variables.

    use DoblDobl_Complex_VecVecs;
    use DoblDobl_Speelpenning_Convolutions;

    s : Link_to_System;
    x : Link_to_VecVec;
    scf : VecVec(1..dim);
    maxit,nbt : integer32 := 0;
    ans : character;
    otp : boolean;
    seri_elapsed,mult_elapsed : Duration := 0.0;

  begin
    DoblDobl_Random_Newton_Homotopy(dim,deg,nbr,pwr,s,x);
    new_line;
    put("Give the maximum number of iterations : "); get(maxit);
    loop
      put("Give the number of tasks (0 to exit) : "); get(nbt);
      exit when (nbt = 0);
      put("Output during multitasking ? (y/n) "); Ask_Yes_or_No(ans);
      otp := (ans = 'y');
      DoblDobl_Complex_VecVecs.Copy(x.all,scf);
      DoblDobl_Run(nbt,dim,maxit,s,scf,seri_elapsed,mult_elapsed,otp);
      DoblDobl_Complex_VecVecs.Clear(scf);
    end loop;
  end DoblDobl_Random_Test;

  procedure QuadDobl_Random_Test
              ( dim,deg,nbr,pwr : in integer32 ) is

  -- DESCRIPTION :
  --   Generates a random Newton homotopy in quad double precision.

  -- ON ENTRY :
  --   dim      dimension of the exponent vectors;
  --   deg      degree of the power series;
  --   nbr      number of products;
  --   pwr      largest power of the variables.

    use QuadDobl_Complex_VecVecs;
    use QuadDobl_Speelpenning_Convolutions;

    s : Link_to_System;
    x : Link_to_VecVec;
    scf : VecVec(1..dim);
    maxit,nbt : integer32 := 0;
    ans : character;
    otp : boolean;
    seri_elapsed,mult_elapsed : Duration := 0.0;

  begin
    QuadDobl_Random_Newton_Homotopy(dim,deg,nbr,pwr,s,x);
    new_line;
    put("Give the maximum number of iterations : "); get(maxit);
    loop
      put("Give the number of tasks (0 to exit) : "); get(nbt);
      exit when (nbt = 0);
      put("Output during multitasking ? (y/n) "); Ask_Yes_or_No(ans);
      otp := (ans = 'y');
      QuadDobl_Complex_VecVecs.Copy(x.all,scf);
      QuadDobl_Run(nbt,dim,maxit,s,scf,seri_elapsed,mult_elapsed,otp);
      QuadDobl_Complex_VecVecs.Clear(scf);
    end loop;
  end QuadDobl_Random_Test;

  procedure Main is

  -- DESCRIPTION :
  --   Launches the test after prompting for the precision.

    prc,ans : character;
    dim,deg,nbr,pwr : integer32 := 0;

  begin
    new_line;
    put_line("Testing Newton's method on power series ...");
    new_line;
    put_line("MENU for the working precision :");
    put_line("  0. standard double precision");
    put_line("  1. double double precision");
    put_line("  2. quad double precision");
    put("Type 0, 1, or 2 to select the precision : ");
    Ask_Alternative(prc,"012");
    new_line;
    put("Generate a random problem ? (y/n) "); Ask_Yes_or_No(ans);
    if ans ='y' then
      new_line;
      put("Give the dimension : "); get(dim);
      put("Give the degree of the power series : "); get(deg);
      put("Give the number of terms in each circuit : "); get(nbr);
      put("Give the largest power of the variables : "); get(pwr);
      case prc is
        when '0' => Standard_Random_Test(dim,deg,nbr,pwr);
        when '1' => DoblDobl_Random_Test(dim,deg,nbr,pwr);
        when '2' => QuadDobl_Random_Test(dim,deg,nbr,pwr);
        when others => null;
      end case;
    else
      case prc is
        when '0' => Standard_Test;
        when '1' => DoblDobl_Test;
        when '2' => QuadDobl_Test;
        when others => null;
      end case;
    end if;
  end Main;

begin
  Main;
end ts_mtnewton;
