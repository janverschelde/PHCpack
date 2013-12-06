with text_io,integer_io;            use text_io,integer_io;
with Timing_Package;                use Timing_Package;

procedure ts_timer is

-- DESCRIPTION :
--   Test of timing package on computing sums of natural numbers.
--   Do "time ts_timer" when running this program to see whether
--   the UNIX timer gives the same results as the printed times.

  n,m,acc : natural;
  timer,totaltimer : timing_widget;

begin
  new_line;
  put_line("Test of timer on computation of sums and products.");
  new_line;
  tstart(totaltimer);
  put("Give number of sums     : "); get(n);
  put("Give number of products : "); get(m);
  tstart(timer);
  acc := 0;
  for i in 1..n loop
    acc := acc + 1;
  end loop;
  tstop(timer);
  new_line;
  put("1 + 1 + .. + 1 : "); put(acc); new_line;
  print_times(timer,"calculating the sum");
  new_line;
  tstart(timer);
  acc := 1;
  for i in 1..m loop
    acc:= acc*acc;
  end loop;
  tstop(timer);
  put("1 * 1 * .. * 1 : "); put(acc); new_line;
  print_times(timer,"calculating the product");
  new_line;
  tstop(totaltimer);
  print_times(totaltimer,"the whole program");
  new_line;
end ts_timer;
