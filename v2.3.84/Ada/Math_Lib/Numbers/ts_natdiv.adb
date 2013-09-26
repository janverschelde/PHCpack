with text_io;                            use text_io;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with System;
with GNAT.Threads;                       use GNAT.Threads;
with Communications_with_User;           use Communications_with_User;
with Timing_Package;                     use Timing_Package;
with Multprec_Natural_Numbers;           use Multprec_Natural_Numbers;
with Multprec_Natural_Numbers_io;        use Multprec_Natural_Numbers_io;
with Multithreading;

procedure ts_natdiv is

-- DESCRIPTION :
--   Test on computing all divisors of a multiprecision natural number.
--   In addition to testing the soundness of the multiprecision arithmetic
--   its main purpose is to provide a nontrivial hard computation.

  procedure Test_Divisor 
              ( n,d : in Natural_Number; test : in boolean;
                isdiv : out boolean; bug : out boolean ) is

  -- DESCRIPTION :
  --   If d divides n, then quotient q and d are written to screen.
  --   If test, then the test d*q = n is done.

    q,dq,r : Natural_Number;

  begin
    bug := false;
    Div(n,d,q,r);
    if Equal(r,0) then
      isdiv := true;
      put(d); put("*"); put(q);
      if not test then
        new_line;
      else
        dq := d*q;
        put(" = "); put(dq);
        if Equal(dq,n) then
          put(" = "); put(n); put_line("  Okay");
        else
          put(" /= "); put(n); put_line("  BUG !???");
          bug := true;
        end if;
        Clear(dq);
      end if;
    else
      isdiv := false;
    end if;
    Clear(q); Clear(r);
  end Test_Divisor;

  procedure Enumerate_Divisors
               ( n : in Natural_Number; test : in boolean ) is

  -- DESCRIPTION :
  --   Writes all divisors and quotients of n to screen, just by plainly
  --   enumerating all candidate divisors.
  --   If test, then q*d = n is tested explicitly.

    timer : Timing_Widget;
    d : Natural_Number;
    isdiv,bug : boolean;

  begin
    d := Create(natural32(2));
    put_line("Checking divisors ...");
    tstart(timer);
    while d < n loop
      Test_Divisor(n,d,test,isdiv,bug);
      exit when bug;
      Add(d,1);
    end loop;
    tstop(timer);
    new_line;
    print_times(standard_output,timer,"checking divisors");
  end Enumerate_Divisors;

  type Array_of_Natural_Numbers is 
    array ( natural32 range <> ) of Natural_Number;

  procedure Divisors ( n,a,b : in Natural_Number; cnt : out natural32;
                       dn : out Array_of_Natural_Numbers;
                       cputime : out duration ) is

  -- DESCRIPTION :
  --   Enumerates all candidate divisors of n, starting from a
  --   and ending at b-1 (so b is not included).

  -- REQUIRED : dn'first equals 1 and dn'last is large enough...

  -- ON ENTRY :
  --   n         a natural number;
  --   a         start of the range of candidate divisors;
  --   b         strict upper bound for the divisor.

  -- ON RETURN :
  --   cnt       number of divisors of n in range a..b-1;  
  --   dn        dn(1..cnt) contains divisors of n;
  --   cputime   cpu time spent on the computation of the divisors.

    timer : Timing_Widget;
    ind : natural32 := 0;
    d,q,r : Natural_Number;

  begin
    tstart(timer);
    Copy(a,d);
    while d < b loop
      Div(n,d,q,r);
      if Equal(r,0) then
        ind := ind + 1;
        Copy(d,dn(ind));
      end if;
      Add(d,1);
      Clear(q); Clear(r);
    end loop;
    cnt := ind;
    Clear(d);
    tstop(timer);
    cputime := Elapsed_User_Time(timer);
  end Divisors;

  procedure Store_all_Divisors ( n : in Natural_Number ) is

    dp : constant natural32 := Decimal_Places(n);
    dn : Array_of_Natural_Numbers(1..10*dp);
    a : Natural_Number := Create(natural32(2));
    cnt : natural32;
    usertime : duration;

  begin
    put("Number of decimal places : "); put(dp,1); new_line;
    Divisors(n,a,n,cnt,dn,usertime);
    put("Found "); put(cnt,1); put(" divisors :");
    for i in 1..cnt loop
      put(" "); put(dn(i));
    end loop;
    new_line;
    put("user cpu time : "); print_time(standard_output,usertime); new_line;
  end Store_all_Divisors;

  procedure Two_Tasked_Divider ( n : in Natural_Number ) is

  -- DESCRIPTION :
  --   Uses Ada tasking ...

    dp : constant natural32 := Decimal_Places(n);
    n2 : Natural_Number := n/2;
    two : Natural_Number := Create(natural32(2));

    task First_Half;
    task Second_Half;
   
    task body First_Half is
 
      dn : Array_of_Natural_Numbers(1..5*dp);
      cnt : natural32;
      t : duration;
     
    begin
      put("Starting first half ...");
      Divisors(n,two,n2,cnt,dn,t);
      put("... done with first half."); new_line;
      put("user cpu time : "); print_time(standard_output,t); new_line;
    end First_Half;
   
    task body Second_Half is
 
      dn : Array_of_Natural_Numbers(1..5*dp);
      cnt : natural32;
      t : duration;
     
    begin
      put("Starting second half ...");
      Divisors(n,n2,n,cnt,dn,t);
      put("... done with second half."); new_line;
      put("user cpu time : "); print_time(standard_output,t); new_line;
    end Second_Half;

  begin
    null;
  end Two_Tasked_Divider;

  procedure Tasked_Divider ( n : in Natural_Number ) is

    timer : Timing_Widget;
    
  begin
    tstart(timer);
    Two_Tasked_Divider(n);
    tstop(timer);
    print_times(standard_output,timer,"divided with two tasks");
  end Tasked_Divider;

  procedure Threaded_Divider ( n : in Natural_Number ) is

    dp : constant natural32 := Decimal_Places(n);
    n2 : Natural_Number := n/2;
    two : Natural_Number := Create(natural32(2));
    dn0,dn1 : Array_of_Natural_Numbers(1..5*dp);
    cnt0,cnt1 : natural32;
    t0,t1 : duration;

    procedure Thread_Code ( id : in integer ) is
    begin
      if id = 1
       then Divisors(n,two,n2,cnt0,dn0,t0);
       else Divisors(n,n2,n,cnt1,dn1,t1);
      end if;
    end Thread_Code;
    procedure Divider is new Multithreading.Start_Threads(Thread_Code);

  begin
    Divider(2);
  end Threaded_Divider;

  procedure Two_Threaded_Enumerator ( n : in Natural_Number ) is

    timer : Timing_Widget;
    done1 : boolean := false;
    done2 : boolean := false;

    function divide ( id : Void_Ptr; p : Void_Ptr ) return Void_Ptr is

      i : constant integer32 := integer32(p.all);
      d : Natural_Number;
      isdiv,bug : boolean;

    begin
      put("Hi from thread "); put(i,1); new_line;
     -- delay 5.0; -- wait for other thread to start
      if i = 1
       then d := Create(natural32(2));
       else d := Create(natural32(3));
      end if;
      while d < n loop
         Test_Divisor(n,d,false,isdiv,bug);
         if isdiv
          then put("-> found by thread "); put(i,1); new_line;
         end if;
         Add(d,2); 
      end loop;
      if i = 1
       then done1 := true;
       else done2 := true;
      end if;
      return null;
    end divide;

  begin
    tstart(timer);
    declare
      a1,a2 : System.address;
      p1 : constant Void_Ptr := new integer'(1);
      p2 : constant Void_Ptr := new integer'(2);
    begin
      a1 := Create_Thread(divide'address,p1,10000,0);
      a2 := Create_Thread(divide'address,p2,10000,0);
    end;
   -- while not done1 and not done2 loop
   --   delay duration(10);
   --   put_line("in busy loop...");
   -- end loop;
    tstop(timer);
    new_line;
    print_times(standard_output,timer,"two threaded enumerator");
  end Two_Threaded_Enumerator;

  procedure Factor_in_Primes ( n : in Natural_Number ) is

    timer : Timing_Widget;
    m,d,q,r : Natural_Number;

  begin
    tstart(timer);
    Copy(n,m);
    d := Create(natural32(2));
    put(n); put(" = ");
    while d < m loop
      Div(m,d,q,r);
      if Equal(r,0) then
        put(d); put(" ");
        Clear(m); m := q;
        Clear(d); d := Create(natural32(2));
      else
        Clear(q);
        Add(d,1);
      end if;
      Clear(r);
    end loop;
    put(m); new_line;
    tstop(timer);
    new_line;
    print_times(standard_output,timer,"factor in primes");
  end Factor_in_Primes;

  procedure Main is

    n : Natural_Number;
    ans : character;

  begin
    new_line;
    put_line("Testing division of natural numbers...");
    new_line;
    put("Give a natural number : "); get(n);
    put("Your number is "); put(n); new_line;
    new_line;
    put_line("MENU for testing division on multiprecision natural numbers :");
    put_line("  1. plainly enumerate all candidate divisors and check;");
    put_line("  2. store all divisors into an array;");
    put_line("  3. use two threads to enumerate all candidate divisors;");
    put_line("  4. compute factorization into prime numbers.");
    put("Type 1, 2, 3, or 4 to make your choice : ");
    Ask_Alternative(ans,"1234");
    new_line;
    if ans = '1' then
      Enumerate_Divisors(n,true);
    elsif ans = '2' then
      Store_all_Divisors(n);
    elsif ans = '3' then
      Tasked_Divider(n);
     -- Threaded_Divider(n);
     -- Two_Threaded_Enumerator(n);
    else
      Factor_in_Primes(n);
    end if;
  end Main;

begin
  Main;
end ts_natdiv;
