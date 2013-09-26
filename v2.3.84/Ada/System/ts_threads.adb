with text_io,integer_io;                use text_io,integer_io;
with System;
with GNAT.Threads;                      use GNAT.Threads;
with Communications_with_User;          use Communications_with_User;
with Multithreading;

procedure ts_threads is

-- DESCRIPTION :
--   Test on how to use threads with GNAT.

  function busy_divisor ( id : Void_Ptr; parm : Void_Ptr ) return Void_Ptr is

    n : constant natural := 100000007;
    r : natural;
    p : constant integer := integer(parm.all);

  begin
    put("hi there from thread "); put(p,1); new_line;
    for k in 1..n loop
      for k in 1..n loop
        r := n mod k;
        if r = 0 then
          put("Thread "); put(p,1);
          put(" found divisor "); put(k,1); new_line;
        end if;
      end loop;
    end loop;
    return null;
  end busy_divisor;

  procedure Busy is

    a1,a2 : System.Address;
    p1 : constant Void_Ptr := new integer'(1);
    p2 : constant Void_Ptr := new integer'(2);

  begin
    put_line("creating two threads ...");
    a1 := Create_Thread(busy_divisor'address,p1,10000,0);
    a2 := Create_Thread(busy_divisor'address,p2,10000,0);
    put_line("threads will remain busy ...");
  end Busy;

  procedure Enumerate is

    a1,a2 : System.Address;
    p1 : constant Void_Ptr := new integer'(1);
    p2 : constant Void_Ptr := new integer'(2);
    done1 : boolean := false;
    done2 : boolean := false;
    n : natural;

    function divisor ( id : Void_Ptr; parm : Void_Ptr ) return Void_Ptr is

      k,r : natural;
      p : constant integer := integer(parm.all);

    begin
      put("hi there from thread "); put(p,1); new_line;
      delay 5.0; -- wait for other thread to start as well
      for j in 1..10000000 loop
        if p = 1 
         then k := 2;
         else k := 3;
        end if;
        while k <= n loop
          r := n mod k;
          if r = 0 then
            put("Thread "); put(p,1);
            put(" found divisor "); put(k,1); new_line;
          end if;
          k := k + 2;
        end loop;
      end loop;
      if p = 1
       then done1 := true;
       else done2 := true;
      end if;
      return null;
    end divisor;

  begin
    new_line;
    put("Give a natural number : "); get(n);
    put_line("starting two threads ...");
    a1 := Create_Thread(divisor'address,p1,10000,0);
    a2 := Create_Thread(divisor'address,p2,10000,0);
    put_line("threads with remain busy ...");
    delay 10000.0;
   -- while not done1 and not done2 loop
   --   null;
   -- end loop;
   -- put_line("both threads are done");
  end Enumerate;

  procedure Test_Multithreading is

    done1 : boolean := false;
    done2 : boolean := false;
    d : integer;
    n,t : natural;

    procedure divisor ( id : in integer ) is

      k,r : natural;

    begin
      put("Thread "); put(id,1);
      put(" looks at number : "); put(n,1); new_line;
     -- delay 5.0; -- allow other thread to grab cpu...
      for j in 1..10000000 loop
        if id = 1 
         then k := 2;
         else k := 3;
        end if;
        while k <= n loop
          r := n mod k;
          if r = 0 then
            put("Thread "); put(id,1);
            put(" found divisor "); put(k,1); new_line;
          end if;
          k := k + 2;
        end loop;
       -- delay 10.0;
        put("Thread "); put(id,1); put_line(" completed one cycle");
       -- delay 10.0; -- take a break
      end loop;
      if id = 1
       then done1 := true;
       else done2 := true;
      end if;
    end divisor;
    procedure Start is new Multithreading.Start_Threads(divisor);

  begin
    new_line;
    put("Enter the delay : "); get(d);
    put("Give a natural number : "); get(n);
    put("Give number of threads : "); get(t);
    put("starting "); put(t,1); put(" threads ...");
    Start(t);
    delay 10.0; -- wait till threads have said hello
    put_line("will enter a busy waiting loop ...");
    delay duration(d);
    if done1
     then put_line("thread 1 is done");
    end if;
    if done2
     then put_line("thread 2 is done");
    end if;
  end Test_Multithreading;

  procedure Main is

    ans : character;

  begin
    new_line;
    put_line("Multithreading with GNAT...");
    new_line;
    put_line("MENU to test using threads : ");
    put_line("  0. test the busyness of two threads;");
    put_line("  1. have two threads find even and odd divisors;"); 
    put_line("  2. do the same with multithreading package.");
    put("Type 0, 1, or 2 to make your choice : ");
    Ask_Alternative(ans,"012");
    if ans = '0' then
      Busy;
    elsif ans = '1' then
      Enumerate;
    else
      Test_Multithreading;
    end if;
  end Main;

begin
  Main;
end ts_threads;
