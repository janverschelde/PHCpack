with text_io;                          use text_io;
with Communications_with_User;         use Communications_with_User;
with Standard_Integer_Numbers;         use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;      use Standard_Integer_Numbers_io;
with Standard_Random_Numbers;
with Multitasking;

procedure ts_tasking is

-- DESCRIPTION :
--   Prompts the user for the number of threads
--   and launches as many tasks as the number given.

  procedure Busy_Sum ( n,m : in integer32 ) is

  -- DESCRIPTION :
  --   Just to keep the tasks very busy, summing n*m numbers.

    sum : integer32;

  begin
    for i in 1..n loop
      sum := 0;
      for j in 1..m loop
        sum := sum + 1;
      end loop;
    end loop;
  end Busy_Sum;

  procedure Start_Silent_Workers ( n : in integer32 ) is

    procedure a_job ( i,n : in integer32 ) is

      w : constant integer32 := 100000;

    begin
      Busy_Sum(w/n,w);
    end a_job;
    procedure do_jobs is new Multitasking.Silent_Workers(a_job);

  begin
    do_jobs(n);
  end Start_Silent_Workers;

  procedure Start_Reporting_Workers ( n : in integer32 ) is

    procedure a_job ( i,n : in integer32 ) is

      w : constant integer32 := 100000;

    begin
      Busy_Sum(w/n,w);
    end a_job;
    procedure do_jobs is new Multitasking.Reporting_Workers(a_job);

  begin
    do_jobs(n);
  end Start_Reporting_Workers;

  procedure Start_Reporting_Looping_Workers ( n : in integer32 ) is

    secret : constant integer32 := 1;

    procedure guess ( i,n : in integer32; continue : out boolean ) is

      r : integer32;
      use Multitasking;

    begin
      r := Standard_Random_Numbers.Random(0,3);
      if r = secret then
        put_line("*** Worker " & to_string(i) & " guessed " 
                               & to_string(r) & " and quits ***");
        continue := false;
      else
        put_line("Worker " & to_string(i) & " guessed " 
                           & to_string(r) & " and continues.");
        continue := true;
      end if;
    end guess;
    procedure do_jobs is new Multitasking.Reporting_Looping_Workers(guess);

  begin
    put("the secret is "); put(secret,1); new_line;
    do_jobs(n);
  end Start_Reporting_Looping_Workers;

  procedure Start_Silent_Looping_Workers ( n : in integer32 ) is

    secret : constant integer32 := 1;

    procedure guess ( i,n : in integer32; continue : out boolean ) is

      r : integer32;
      use Multitasking;

    begin
      r := Standard_Random_Numbers.Random(0,3);
      if r = secret then
        put_line("*** Worker " & to_string(i) & " guessed " 
                               & to_string(r) & " and quits ***");
        continue := false;
      else
        put_line("Worker " & to_string(i) & " guessed " 
                           & to_string(r) & " and continues.");
        continue := true;
      end if;
    end guess;
    procedure do_jobs is new Multitasking.Silent_Looping_Workers(guess);

  begin
    put("the secret is "); put(secret,1); new_line;
    do_jobs(n);
  end Start_Silent_Looping_Workers;


  procedure Main is

    n : integer32 := 0;
    ans,choice : character;

  begin
    new_line;
    put_line("MENU for testing tasking : ");
    put_line("  1. test a busy sum ...");
    put_line("  2. test looping workers");
    put("Type 1 or 2 to make your choice : ");
    Ask_Alternative(choice,"12");
    new_line;
    put("Give the number of workers : "); get(n);
    put("Do you want intermediate output ? ");
    Ask_Yes_or_No(ans);
    if choice = '1' then
      if ans = 'y' 
       then Start_Reporting_Workers(n);
       else Start_Silent_Workers(n);
      end if;
    else
      if ans = 'y'
       then Start_Reporting_Looping_Workers(n);
       else Start_Silent_Looping_Workers(n);
      end if;
    end if;
  end Main;

begin
  Main;
end ts_tasking;
