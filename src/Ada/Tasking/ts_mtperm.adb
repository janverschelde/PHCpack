with Ada.Calendar;
with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User; 
with Timing_Package;                     use Timing_Package;
with Time_Stamps;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Random_Numbers;
with Standard_Integer_Vectors;
with Standard_Integer_Vectors_io;        use Standard_Integer_Vectors_io;
with Standard_Integer_VecVecs;
with Standard_Integer_Matrices;
with Standard_Integer_Matrices_io;       use Standard_Integer_Matrices_io;
with Boolean_Matrices;
with Boolean_Matrices_io;                use Boolean_Matrices_io;
with Standard_Random_Matrices;
with Integer_Permanents;                 use Integer_Permanents;
with Multitasking_Integer_Permanents;

procedure ts_mtperm is

-- DESCRIPTION :
--   Interactive test on the permanent computation of an integer matrix.

  procedure Integer_Random_Test ( dim : in natural32 ) is

  -- DESCRIPTION :
  --   Generates a random integer matrix of the given dimension dim
  --   and computes its permanent.

    idm : constant integer32 := integer32(dim);
    mat : constant Standard_Integer_Matrices.Matrix(1..idm,1..idm)
        := Standard_Random_Matrices.Random_Matrix(dim,dim,0,9);
    per : integer64;
    ans : character;
    timer : Timing_Widget;

  begin
    put("A random "); put(dim,1); put("-by-"); put(dim,1);
    put_line(" matrix : "); put(mat);
    put("Do you want output of the factors ? (y/n) ");
    Ask_Yes_or_No(ans);
    tstart(timer);
    if ans = 'y' 
     then per := Permanent(standard_output,mat);
     else per := Permanent(mat);
    end if;
    tstop(timer);
    put("The permanent : "); put(per,1); new_line;
    new_line;
    print_times(standard_output,timer,"Permanent Computation");
  end Integer_Random_Test;

  procedure Boolean_Random_Test ( dim : in natural32 ) is

  -- DESCRIPTION :
  --   Generates a random Boolean matrix of the given dimension dim
  --   and computes its permanent.

    idm : constant integer32 := integer32(dim);
    mat : Boolean_Matrices.Matrix(1..idm,1..idm)
        := Standard_Random_Matrices.Random_Matrix(dim);
    per : integer64;
    ans : character;
    timer : Timing_Widget;

  begin
    put("A random "); put(dim,1); put("-by-"); put(dim,1);
    put_line(" matrix : "); put(mat);
    put("Do you want output of the factors ? (y/n) ");
    Ask_Yes_or_No(ans);
    tstart(timer);
    if ans = 'y' 
     then per := Permanent(standard_output,mat);
     else per := Permanent(mat);
    end if;
    tstop(timer);
    put("The permanent : "); put(per,1); new_line;
    new_line;
    print_times(standard_output,timer,"Permanent Computation");
  end Boolean_Random_Test;

  procedure Random_Test is

  -- DESCRIPTION :
  --   Prompts the user for the dimension
  --   and asks whether the matrix is Boolean or not.

    dim : natural32 := 0;
    ans : character;

  begin
    new_line;
    put("Give the dimension : "); get(dim);
    put("Generate a Boolean matrix ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans = 'y'
     then Boolean_Random_Test(dim);
     else Integer_Random_Test(dim);
    end if;
  end Random_Test;

  procedure Start_Integer_Random_Test ( dim : in natural32 ) is

  -- DESCRIPTION :
  --   Generates a random integer matrix of the given dimension dim
  --   and computes its permanent, starting with the column indices
  --   for a number of rows first.

    idm : constant integer32 := integer32(dim);
    mat : constant Standard_Integer_Matrices.Matrix(1..idm,1..idm)
        := Standard_Random_Matrices.Random_Matrix(dim,dim,0,9);
    nbr : integer32 := 0;

  begin
    put("A random "); put(dim,1); put("-by-"); put(dim,1);
    put_line(" matrix : "); put(mat);
    put("Give the number of rows : "); get(nbr);
    declare
      nbstc : constant integer32
            := Number_of_Start_Columns(integer32(dim),nbr);
      stcol : Standard_Integer_VecVecs.VecVec(1..nbstc);
      cols : Standard_Integer_Vectors.Vector(mat'range(1))
           := (mat'range(1) => 0);
      cnts : Standard_Integer_Vectors.Vector(mat'range(1))
           := (mat'range(1) => 1);
      idx : integer32 := 0;
      per : integer64;
      ans : character;
    begin
      Start_Columns(mat'first(1),nbr,mat,cols,cnts,idx,stcol);
      put_line("The start column indices : ");
      for k in 1..idx loop
        put(stcol(k)); new_line;
      end loop;
      put("Do you want output of the factors ? (y/n) ");
      Ask_Yes_or_No(ans);
      if ans = 'y'
       then per := Start_Permanent(standard_output,nbr+1,stcol(1..idx),mat);
       else per := Start_Permanent(nbr+1,stcol(1..idx),mat);
      end if;
      put("The permanent : "); put(per,1); new_line;
      per := Permanent(mat);
      put("The permanent : "); put(per,1); new_line;
    end;
  end Start_Integer_Random_Test;

  procedure Start_Boolean_Random_Test ( dim : in natural32 ) is

  -- DESCRIPTION :
  --   Generates a random Boolean matrix of the given dimension dim
  --   and computes its permanent, starting with the column indices
  --   for a number of rows first.

    idm : constant integer32 := integer32(dim);
    mat : Boolean_Matrices.Matrix(1..idm,1..idm)
        := Standard_Random_Matrices.Random_Matrix(dim);
    nbr : integer32 := 0;

  begin
    put("A random "); put(dim,1); put("-by-"); put(dim,1);
    put_line(" matrix : "); put(mat);
    put("Give the number of rows : "); get(nbr);
    declare
      nbstc : constant integer32
            := Number_of_Start_Columns(integer32(dim),nbr);
      stcol : Standard_Integer_VecVecs.VecVec(1..nbstc);
      cols : Standard_Integer_Vectors.Vector(mat'range(1))
           := (mat'range(1) => 0);
      cnts : Standard_Integer_Vectors.Vector(mat'range(1))
           := (mat'range(1) => 1);
      idx : integer32 := 0;
      per : integer64;
      ans : character;
    begin
      Start_Columns(mat'first(1),nbr,mat,cols,cnts,idx,stcol);
      put_line("The start column indices : ");
      for k in 1..idx loop
        put(stcol(k)); new_line;
      end loop;
      put("Do you want output of the factors ? (y/n) ");
      Ask_Yes_or_No(ans);
      if ans = 'y'
       then per := Start_Permanent(standard_output,nbr+1,stcol(1..idx),mat);
       else per := Start_Permanent(nbr+1,stcol(1..idx),mat);
      end if;
      put("The permanent : "); put(per,1); new_line;
      per := Permanent(mat);
      put("The permanent : "); put(per,1); new_line;
    end;
  end Start_Boolean_Random_Test;

  procedure Start_Random_Test is

  -- DESCRIPTION :
  --   Prompts the user for the dimension
  --   and asks whether the matrix is Boolean or not.

    dim : natural32 := 0;
    ans : character;

  begin
    new_line;
    put("Give the dimension : "); get(dim);
    put("Generate a Boolean matrix ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans = 'y'
     then Start_Boolean_Random_Test(dim);
     else Start_Integer_Random_Test(dim);
    end if;
  end Start_Random_Test;

  procedure Multitasked_Boolean_Test ( dim : in natural32 ) is

  -- DESCRIPTION :
  --   Generates a random integer matrix of the given dimension dim
  --   and computes its permanent.
  --   Elapsed wall clock times in seconds are written.

    idm : constant integer32 := integer32(dim);
    mat : constant Boolean_Matrices.Matrix(1..idm,1..idm)
        := Standard_Random_Matrices.Random_Matrix(dim);
    per : integer64;
    nbrows,ntasks : integer32 := 0;
    ans : character;
    verb : boolean;
    multstart,multstop,seristart,seristop : Ada.Calendar.Time;

  begin
    put("A random "); put(dim,1); put("-by-"); put(dim,1);
    put_line(" matrix : "); put(mat);
    put("Give the number of first rows : "); get(nbrows);
    put("Give the number of tasks : "); get(ntasks);
    put("Output during computations ? "); Ask_Yes_or_No(ans);
    verb := (ans = 'y');
    multstart := Ada.Calendar.Clock;
    per := Multitasking_Integer_Permanents.Permanent(mat,nbrows,ntasks,verb);
    multstop := Ada.Calendar.Clock;
    put("The permanent : "); put(per,1); new_line;
    seristart := Ada.Calendar.Clock;
    per := Permanent(mat);
    seristop := Ada.Calendar.Clock;
    put("The permanent : "); put(per,1); new_line;
    new_line;
    put_line("Elapsed time on multitasking : ");
    Time_Stamps.Write_Elapsed_Time(standard_output,multstart,multstop);
    new_line;
    put_line("Elapsed time without multitasking : ");
    Time_Stamps.Write_Elapsed_Time(standard_output,seristart,seristop);
  end Multitasked_Boolean_Test;

  procedure Multitasked_Integer_Test ( dim : in natural32 ) is

  -- DESCRIPTION :
  --   Generates a random integer matrix of the given dimension dim
  --   and computes its permanent.

    idm : constant integer32 := integer32(dim);
    mat : constant Standard_Integer_Matrices.Matrix(1..idm,1..idm)
        := Standard_Random_Matrices.Random_Matrix(dim,dim,0,9);
    per : integer64;
    nbrows,ntasks : integer32 := 0;
    ans : character;
    verb : boolean;
    multstart,multstop,seristart,seristop : Ada.Calendar.Time;

  begin
    put("A random "); put(dim,1); put("-by-"); put(dim,1);
    put_line(" matrix : "); put(mat);
    put("Give the number of first rows : "); get(nbrows);
    put("Give the number of tasks : "); get(ntasks);
    put("Output during computations ? "); Ask_Yes_or_No(ans);
    verb := (ans = 'y');
    multstart := Ada.Calendar.Clock;
    per := Multitasking_Integer_Permanents.Permanent(mat,nbrows,ntasks,verb);
    multstop := Ada.Calendar.Clock;
    put("The permanent : "); put(per,1); new_line;
    seristart := Ada.Calendar.Clock;
    per := Permanent(mat);
    seristop := Ada.Calendar.Clock;
    put("The permanent : "); put(per,1); new_line;
    new_line;
    put_line("Elapsed time on multitasking : ");
    Time_Stamps.Write_Elapsed_Time(standard_output,multstart,multstop);
    new_line;
    put_line("Elapsed time without multitasking : ");
    Time_Stamps.Write_Elapsed_Time(standard_output,seristart,seristop);
  end Multitasked_Integer_Test;

  procedure Multitasked_Test is

  -- DESCRIPTION :
  --   Prompts the user for the dimension
  --   and asks whether the matrix is Boolean or not.

    dim : natural32 := 0;
    ans : character;

  begin
    new_line;
    put("Give the dimension : "); get(dim);
    put("Generate a Boolean matrix ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans = 'y'
     then Multitasked_Boolean_Test(dim);
     else Multitasked_Integer_Test(dim);
    end if;
  end Multitasked_Test;

  procedure Doubling_Multitasked_Boolean_Test ( dim : in natural32 ) is

  -- DESCRIPTION :
  --   Generates a random integer matrix of the given dimension dim
  --   and computes its permanent.
  --   Elapsed wall clock times in seconds are written.

    idm : constant integer32 := integer32(dim);
    mat : constant Boolean_Matrices.Matrix(1..idm,1..idm)
        := Standard_Random_Matrices.Random_Matrix(dim);
    per : integer64;
    nbrows,ntasks,maxnbtasks : integer32 := 0;
    ans : character;
    verb : boolean;
    timestart,timestop : Ada.Calendar.Time;

  begin
    put("A random "); put(dim,1); put("-by-"); put(dim,1);
    put_line(" matrix : "); put(mat);
    put("Give the number of first rows : "); get(nbrows);
    put("Give the largest number of tasks : "); get(maxnbtasks);
    put("Output during computations ? "); Ask_Yes_or_No(ans);
    verb := (ans = 'y');
    timestart := Ada.Calendar.Clock;
    per := Permanent(mat);
    timestop := Ada.Calendar.Clock;
    put("The permanent : "); put(per,1); new_line;
    put_line("Elapsed time without multitasking : ");
    Time_Stamps.Write_Elapsed_Time(standard_output,timestart,timestop);
    ntasks := 2;
    while ntasks <= maxnbtasks loop
      timestart := Ada.Calendar.Clock;
      per := Multitasking_Integer_Permanents.Permanent(mat,nbrows,ntasks,verb);
      timestop := Ada.Calendar.Clock;
      put("The permanent : "); put(per,1); new_line;
      put("With "); put(ntasks,1); put(" tasks.  ");
      Time_Stamps.Write_Elapsed_Time(standard_output,timestart,timestop);
      ntasks := 2*ntasks;
    end loop;
  end Doubling_Multitasked_Boolean_Test;

  procedure Doubling_Multitasked_Integer_Test ( dim : in natural32 ) is

  -- DESCRIPTION :
  --   Generates a random integer matrix of the given dimension dim
  --   and computes its permanent.

    idm : constant integer32 := integer32(dim);
    mat : constant Standard_Integer_Matrices.Matrix(1..idm,1..idm)
        := Standard_Random_Matrices.Random_Matrix(dim,dim,0,9);
    per : integer64;
    nbrows,ntasks,maxnbtasks : integer32 := 0;
    ans : character;
    verb : boolean;
    timestart,timestop : Ada.Calendar.Time;

  begin
    put("A random "); put(dim,1); put("-by-"); put(dim,1);
    put_line(" matrix : "); put(mat);
    put("Give the number of first rows : "); get(nbrows);
    put("Give the maximum number of tasks : "); get(maxnbtasks);
    put("Output during computations ? "); Ask_Yes_or_No(ans);
    verb := (ans = 'y');
    timestart := Ada.Calendar.Clock;
    per := Permanent(mat);
    timestop := Ada.Calendar.Clock;
    put("The permanent : "); put(per,1); new_line;
    put_line("Elapsed time without multitasking : ");
    Time_Stamps.Write_Elapsed_Time(standard_output,timestart,timestop);
    ntasks := 2;
    while ntasks <= maxnbtasks loop
      timestart := Ada.Calendar.Clock;
      per := Multitasking_Integer_Permanents.Permanent(mat,nbrows,ntasks,verb);
      timestop := Ada.Calendar.Clock;
      put("The permanent : "); put(per,1); new_line;
      put("With "); put(ntasks,1); put(" tasks.  ");
      Time_Stamps.Write_Elapsed_Time(standard_output,timestart,timestop);
      ntasks := 2*ntasks;
    end loop;
  end Doubling_Multitasked_Integer_Test;

  procedure Doubling_Multitasked_Test is

  -- DESCRIPTION :
  --   Prompts the user for the dimension
  --   and asks whether the matrix is Boolean or not.
  --   Starting with the serial run,
  --   the number of tasks doubles in each stage.

    dim : natural32 := 0;
    ans : character;

  begin
    new_line;
    put("Give the dimension : "); get(dim);
    put("Generate a Boolean matrix ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans = 'y'
     then Doubling_Multitasked_Boolean_Test(dim);
     else Doubling_Multitasked_Integer_Test(dim);
    end if;
  end Doubling_Multitasked_Test;

  procedure Main is

    ans : character;

  begin
    new_line;
    put_line("MENU to test a permament computation : ");
    put_line("  1. test on a randomly generated matrix;");
    put_line("  2. compute start columns for the first rows;");
    put_line("  3. test multitasked permanent computation;");
    put_line("  4. double the number of tasks in each stage.");
    put("Type 1, 2, 3, or 4 to choose a test : ");
    Ask_Alternative(ans,"1234");
    case ans is
      when '1' => Random_Test;
      when '2' => Start_Random_Test;
      when '3' => Multitasked_Test;
      when '4' => Doubling_Multitasked_Test;
      when others => null;
    end case;
  end Main;

begin
  Main;
end ts_mtperm;
