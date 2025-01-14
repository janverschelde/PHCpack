with text_io;                            use text_io;
with Timing_Package;                     use Timing_Package;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;

procedure ts_leibniz is

-- This example is based on section 3.2.2 on loop unrolling in 
-- Scientific Programming and Computer Architecture
-- by Divakar Viswanath, Springer-Verlag, 2017.

-- The branching in the straighforward code below prevents
-- a pipelined execution of the floating-point operations.

  function leibniz1 ( n : integer32 ) return double_float is

  -- DESCRIPTION :
  --   Straighforward implementation of the Leibniz series to
  --   approximate pi/4 as 1 - 1/3 + 1/5 - 1/7 + 1/9 - 1/11 + ...

    s : double_float := 1.0;

  begin
    for i in 1..n loop
      if i mod 2 = 1
       then s := s - 1.0/double_float(2*i + 1);
       else s := s + 1.0/double_float(2*i + 1);
      end if;
    end loop;
    return s;
  end leibniz1;

  function leibniz2 ( n : integer32 ) return double_float is

  -- DESCRIPTION :
  --   Reformulated implementation of the Leibniz series to
  --   approximate pi/4 as 1 + 1/5 + 1/9 + ... - 1/3 - 1/7 - 1/11 + ...

    s : double_float := 1.0;
    i : integer32;

  begin
    i := 2;
    while i <= n loop
       s := s + 1.0/double_float(2*i + 1);
       i := i + 2;
    end loop;
    i := 1;
    while i <= n loop
       s := s - 1.0/double_float(2*i + 1);
       i := i + 2;
    end loop;
    return s;
  end leibniz2;

  procedure Main is

  -- DESCRIPTION :
  --   Runs the tests.

    nbr : constant integer32 := 10**8;
    pi1,pi2 : double_float;
    timer1,timer2 : Timing_Widget;

  begin
    tstart(timer1);
    pi1 := leibniz1(nbr);
    tstop(timer1);
    put("pi ~"); put(4.0*pi1); new_line;
    tstart(timer2);
    pi2 := leibniz2(nbr);
    tstop(timer2);
    put("pi ~"); put(4.0*pi2); new_line;
    new_line;
    print_times(standard_output,timer1,"straightforward Leibniz rule");
    new_line;
    print_times(standard_output,timer2,"reformulated Leibniz rule");
  end Main;

begin
  Main;
end ts_leibniz;
