with text_io,integer_io;                 use text_io,integer_io;
with Communications_with_User;           use Communications_with_User;
with Timing_Package;                     use Timing_Package;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;
with Standard_Natural_Vectors;
with Standard_Natural_Vectors_io;        use Standard_Natural_Vectors_io;
with Standard_Complex_Vectors;
with Standard_Complex_Matrices;
with Standard_Complex_Linear_Solvers;
with Standard_Complex_Norms_Equals;
with Standard_Random_Vectors;
with Standard_Random_Matrices;
with Double_Double_Numbers;              use Double_Double_Numbers;
with Double_Double_Numbers_io;           use Double_Double_Numbers_io;
with DoblDobl_Complex_Vectors;
with DoblDobl_Complex_Matrices;
with DoblDobl_Complex_Linear_Solvers;
with DoblDobl_Complex_Vector_Norms;
with DoblDobl_Random_Vectors;
with DoblDobl_Random_Matrices;
with Quad_Double_Numbers;                use Quad_Double_Numbers;
with Quad_Double_Numbers_io;             use Quad_Double_Numbers_io;
with QuadDobl_Complex_Vectors;
with QuadDobl_Complex_Matrices;
with QuadDobl_Complex_Linear_Solvers;
with QuadDobl_Complex_Vector_Norms;
with QuadDobl_Random_Vectors;
with QuadDobl_Random_Matrices;
with Multitasking_Linear_Solvers;

procedure ts_mtlinsol is

-- DESCRIPTION :
--   Interactive development of a linear system solver with multitasking.

  procedure Standard_Test ( n,t : in positive; output : in boolean ) is

  -- DESCRIPTION :
  --   Generates a random n-by-n linear system and tests the
  --   linear solver using t tasks.

    a : constant Standard_Complex_Matrices.Matrix(1..n,1..n)
      := Standard_Random_Matrices.Random_Matrix(n,n);
    c : Standard_Complex_Matrices.Matrix(1..n,1..n) := a;
    b : constant Standard_Complex_Vectors.Vector(1..n)
      := Standard_Random_Vectors.Random_Vector(1,n);
    x : Standard_Complex_Vectors.Vector(1..n) := b;
    r : Standard_Complex_Vectors.Vector(1..n);
    residue : double_float;
    ipvt : Standard_Natural_Vectors.Vector(1..n);
    info : integer;

    use Standard_Complex_Vectors,Standard_Complex_Matrices;

  begin
   -- Standard_Complex_Linear_Solvers.lufac(c,n,ipvt,info);
    if output
     then Multitasking_Linear_Solvers.reporting_lufac(t,c,n,ipvt,info);
     else Multitasking_Linear_Solvers.silent_lufac(t,c,n,ipvt,info);
    end if;
    put("info after lufac : "); put(info,1); new_line;
    put("ipvt = "); put(ipvt); new_line;
   -- Standard_Complex_Linear_Solvers.lusolve(c,n,ipvt,x);
    if output
     then Multitasking_Linear_Solvers.reporting_lusolve(t,c,n,ipvt,x);
     else Multitasking_Linear_Solvers.silent_lusolve(t,c,n,ipvt,x);
    end if;
    r := b - a*x;
    residue := Standard_Complex_Norms_Equals.Max_Norm(r);
    put("residue after lusolve : "); put(residue,3); new_line;
  end Standard_Test;

  procedure Standard_Performance_lusolve
              ( n,t : in positive; m : in natural ) is

  -- DESCRIPTION :
  --   Generates a random n-by-n linear system and runs lusolve
  --   with t tasks m times after a lufac on a random system.

    a : constant Standard_Complex_Matrices.Matrix(1..n,1..n)
      := Standard_Random_Matrices.Random_Matrix(n,n);
    c : Standard_Complex_Matrices.Matrix(1..n,1..n) := a;
    b : constant Standard_Complex_Vectors.Vector(1..n)
      := Standard_Random_Vectors.Random_Vector(1,n);
    x : Standard_Complex_Vectors.Vector(1..n) := b;
    ipvt : Standard_Natural_Vectors.Vector(1..n);
    info : integer;
    timer : timing_widget;

  begin
    Standard_Complex_Linear_Solvers.lufac(c,n,ipvt,info);
    put("info after lufac : "); put(info,1); new_line;
    tstart(timer);
    for i in 1..m loop
      Multitasking_Linear_Solvers.silent_lusolve(t,c,n,ipvt,x);
    end loop;
    tstop(timer);
    new_line;
    print_times(standard_output,timer,"standard performance lusolve");
  end Standard_Performance_lusolve;

  procedure Standard_Performance_lufac
              ( n,t : in positive; m : in natural ) is

  -- DESCRIPTION :
  --   Generates a random n-by-n linear system and tests the
  --   linear solver using t tasks.

    a : constant Standard_Complex_Matrices.Matrix(1..n,1..n)
      := Standard_Random_Matrices.Random_Matrix(n,n);
    c : Standard_Complex_Matrices.Matrix(1..n,1..n);
    ipvt : Standard_Natural_Vectors.Vector(1..n);
    info : integer;
    timer : timing_widget;

  begin
    tstart(timer);
    for i in 1..m loop
      c := a;
      Multitasking_Linear_Solvers.silent_lufac(t,c,n,ipvt,info);
    end loop;
    tstop(timer);
    new_line;
    print_times(standard_output,timer,"standard performance lufac");
  end Standard_Performance_lufac;

  procedure DoblDobl_Test ( n,t : in positive; output : in boolean ) is

  -- DESCRIPTION :
  --   Generates a random n-by-n linear system and tests the
  --   linear solver using t tasks.

    a : constant DoblDobl_Complex_Matrices.Matrix(1..n,1..n)
      := DoblDobl_Random_Matrices.Random_Matrix(n,n);
    c : DoblDobl_Complex_Matrices.Matrix(1..n,1..n) := a;
    b : constant DoblDobl_Complex_Vectors.Vector(1..n)
      := DoblDobl_Random_Vectors.Random_Vector(1,n);
    x : DoblDobl_Complex_Vectors.Vector(1..n) := b;
    r : DoblDobl_Complex_Vectors.Vector(1..n);
    residue : double_double;
    ipvt : Standard_Natural_Vectors.Vector(1..n);
    info : integer;

    use DoblDobl_Complex_Vectors,DoblDobl_Complex_Matrices;

  begin
   -- DoblDobl_Complex_Linear_Solvers.lufac(c,n,ipvt,info);
    if output
     then Multitasking_Linear_Solvers.reporting_lufac(t,c,n,ipvt,info);
     else Multitasking_Linear_Solvers.silent_lufac(t,c,n,ipvt,info);
    end if;
    put("info after lufac : "); put(info,1); new_line;
    put("ipvt = "); put(ipvt); new_line;
   -- DoblDobl_Complex_Linear_Solvers.lusolve(c,n,ipvt,x);
    if output
     then Multitasking_Linear_Solvers.reporting_lusolve(t,c,n,ipvt,x);
     else Multitasking_Linear_Solvers.silent_lusolve(t,c,n,ipvt,x);
    end if;
    r := b - a*x;
    residue := DoblDobl_Complex_Vector_Norms.Max_Norm(r);
    put("residue after lusolve : "); put(residue,3); new_line;
  end DoblDobl_Test;

  procedure DoblDobl_Performance_lusolve
              ( n,t : in positive; m : in natural ) is

  -- DESCRIPTION :
  --   Generates a random n-by-n linear system and runs lusolve
  --   with t tasks m times after a lufac on a random system.

    a : constant DoblDobl_Complex_Matrices.Matrix(1..n,1..n)
      := DoblDobl_Random_Matrices.Random_Matrix(n,n);
    c : DoblDobl_Complex_Matrices.Matrix(1..n,1..n) := a;
    b : constant DoblDobl_Complex_Vectors.Vector(1..n)
      := DoblDobl_Random_Vectors.Random_Vector(1,n);
    x : DoblDobl_Complex_Vectors.Vector(1..n) := b;
    ipvt : Standard_Natural_Vectors.Vector(1..n);
    info : integer;
    timer : timing_widget;

  begin
    DoblDobl_Complex_Linear_Solvers.lufac(c,n,ipvt,info);
    put("info after lufac : "); put(info,1); new_line;
    tstart(timer);
    for i in 1..m loop
      Multitasking_Linear_Solvers.silent_lusolve(t,c,n,ipvt,x);
    end loop;
    tstop(timer);
    new_line;
    print_times(standard_output,timer,"dobldobl performance lusolve");
  end DoblDobl_Performance_lusolve;

  procedure DoblDobl_Performance_lufac
              ( n,t : in positive; m : in natural ) is

  -- DESCRIPTION :
  --   Generates a random n-by-n linear system and tests the
  --   linear solver using t tasks.

    a : constant DoblDobl_Complex_Matrices.Matrix(1..n,1..n)
      := DoblDobl_Random_Matrices.Random_Matrix(n,n);
    c : DoblDobl_Complex_Matrices.Matrix(1..n,1..n);
    ipvt : Standard_Natural_Vectors.Vector(1..n);
    info : integer;
    timer : timing_widget;

  begin
    tstart(timer);
    for i in 1..m loop
      c := a;
      Multitasking_Linear_Solvers.silent_lufac(t,c,n,ipvt,info);
    end loop;
    tstop(timer);
    new_line;
    print_times(standard_output,timer,"dobldobl performance lufac");
  end DoblDobl_Performance_lufac;

  procedure QuadDobl_Test ( n,t : in positive; output : in boolean ) is

  -- DESCRIPTION :
  --   Generates a random n-by-n linear system and tests the
  --   linear solver using t tasks.

    a : constant QuadDobl_Complex_Matrices.Matrix(1..n,1..n)
      := QuadDobl_Random_Matrices.Random_Matrix(n,n);
    c : QuadDobl_Complex_Matrices.Matrix(1..n,1..n) := a;
    b : constant QuadDobl_Complex_Vectors.Vector(1..n)
      := QuadDobl_Random_Vectors.Random_Vector(1,n);
    x : QuadDobl_Complex_Vectors.Vector(1..n) := b;
    r : QuadDobl_Complex_Vectors.Vector(1..n);
    residue : quad_double;
    ipvt : Standard_Natural_Vectors.Vector(1..n);
    info : integer;

    use QuadDobl_Complex_Vectors,QuadDobl_Complex_Matrices;

  begin
   -- QuadDobl_Complex_Linear_Solvers.lufac(c,n,ipvt,info);
    if output
     then Multitasking_Linear_Solvers.reporting_lufac(t,c,n,ipvt,info);
     else Multitasking_Linear_Solvers.silent_lufac(t,c,n,ipvt,info);
    end if;
    put("info after lufac : "); put(info,1); new_line;
    put("ipvt = "); put(ipvt); new_line;
   -- QuadDobl_Complex_Linear_Solvers.lusolve(c,n,ipvt,x);
    if output
     then Multitasking_Linear_Solvers.reporting_lusolve(t,c,n,ipvt,x);
     else Multitasking_Linear_Solvers.silent_lusolve(t,c,n,ipvt,x);
    end if;
    r := b - a*x;
    residue := QuadDobl_Complex_Vector_Norms.Max_Norm(r);
    put("residue after lusolve : "); put(residue,3); new_line;
  end QuadDobl_Test;

  procedure QuadDobl_Performance_lusolve
              ( n,t : in positive; m : in natural ) is

  -- DESCRIPTION :
  --   Generates a random n-by-n linear system and runs lusolve
  --   with t tasks m times after a lufac on a random system.

    a : constant QuadDobl_Complex_Matrices.Matrix(1..n,1..n)
      := QuadDobl_Random_Matrices.Random_Matrix(n,n);
    c : QuadDobl_Complex_Matrices.Matrix(1..n,1..n) := a;
    b : constant QuadDobl_Complex_Vectors.Vector(1..n)
      := QuadDobl_Random_Vectors.Random_Vector(1,n);
    x : QuadDobl_Complex_Vectors.Vector(1..n) := b;
    ipvt : Standard_Natural_Vectors.Vector(1..n);
    info : integer;
    timer : timing_widget;

  begin
    QuadDobl_Complex_Linear_Solvers.lufac(c,n,ipvt,info);
    put("info after lufac : "); put(info,1); new_line;
    tstart(timer);
    for i in 1..m loop
      Multitasking_Linear_Solvers.silent_lusolve(t,c,n,ipvt,x);
    end loop;
    tstop(timer);
    new_line;
    print_times(standard_output,timer,"dobldobl performance lusolve");
  end QuadDobl_Performance_lusolve;

  procedure QuadDobl_Performance_lufac
              ( n,t : in positive; m : in natural ) is

  -- DESCRIPTION :
  --   Generates a random n-by-n linear system and tests the
  --   linear solver using t tasks.

    a : constant QuadDobl_Complex_Matrices.Matrix(1..n,1..n)
      := QuadDobl_Random_Matrices.Random_Matrix(n,n);
    c : QuadDobl_Complex_Matrices.Matrix(1..n,1..n);
    ipvt : Standard_Natural_Vectors.Vector(1..n);
    info : integer;
    timer : timing_widget;

  begin
    tstart(timer);
    for i in 1..m loop
      c := a;
      Multitasking_Linear_Solvers.silent_lufac(t,c,n,ipvt,info);
    end loop;
    tstop(timer);
    new_line;
    print_times(standard_output,timer,"quaddobl performance lufac");
  end QuadDobl_Performance_lufac;

  procedure Main is

    n,t : positive;
    ans,choice : character;
    output : boolean;
    m : natural;
 
  begin
    new_line;
    put_line("Testing linear system solvers with multitasking.");
    new_line;
    put_line("MENU to test multitasking linear system solving : ");
    put_line("  1. use standard complex arithmetic;");
    put_line("  2. use double double complex arithmetic.");
    put_line("  3. use quad double complex arithmetic;");
    put_line("  4. test performance of standard complex lusolve;");
    put_line("  5. test performance of double double complex lusolve;");
    put_line("  6. test performance of quad double complex lusolve;");
    put_line("  7. test performance of standard complex lufac;");
    put_line("  8. test performance of double double complex lufac;");
    put_line("  9. test performance of quad double complex lufac.");
    put("Type 1, 2, 3, 4, 5, 6, 7, 8, or 9 to make your choice : ");
    Ask_Alternative(choice,"123456789");
    new_line;
    put("Give the dimension : "); get(n);
    put("Give the number of tasks : "); get(t);
    if choice = '1' or choice = '2' or choice = '3' then
      put("Output during computations ? (y/n) ");
      Ask_Yes_or_No(ans);
      output := (ans = 'y');
    else
      put("Give frequency for performance : "); get(m);
    end if;
    case choice is
      when '1' => Standard_Test(n,t,output);
      when '2' => DoblDobl_Test(n,t,output);
      when '3' => QuadDobl_Test(n,t,output);
      when '4' => Standard_Performance_lusolve(n,t,m);
      when '5' => DoblDobl_Performance_lusolve(n,t,m);
      when '6' => QuadDobl_Performance_lusolve(n,t,m);
      when '7' => Standard_Performance_lufac(n,t,m);
      when '8' => DoblDobl_Performance_lufac(n,t,m);
      when '9' => QuadDobl_Performance_lufac(n,t,m);
      when others => null;
    end case;
  end Main;

begin
  Main;
end ts_mtlinsol;
