with text_io,integer_io;               use text_io,integer_io;
with Standard_Natural_Numbers;         use Standard_Natural_Numbers;
with Standard_Integer_Numbers;         use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;      use Standard_Integer_Numbers_io;
with Multitasking,Semaphore;

procedure ts_jobs is

-- DESCRIPTION :
--   Prompts the user for the number of threads, the number of jobs,
--   and launches as many tasks as the number given.

  cnt : natural32 := 0;

  procedure Start_Workers ( n,m : in integer32 ) is

  -- DESCRIPTION :
  --   Starts n tasks to do m jobs.

    procedure do_next_job ( i,n : in integer32 ) is

    -- DESCRIPTION :
    --   Increases the job counter and does the next job if
    --   the counter is still less than or equal to m.
    --   After each critical section, a message is printed.

      s : Semaphore.Lock;
      myjob : natural32;

    begin
      loop
        Semaphore.Request(s);
        cnt := cnt + 1; myjob := cnt;
        Semaphore.Release(s);
        if integer32(myjob) <= m then
           put_line("Worker " & Multitasking.to_string(natural32(i))
            & " will do job " & Multitasking.to_string(myjob) & ".");
           delay 1.0;
        else
          put_line("No more jobs for worker "
            & Multitasking.to_string(natural32(i)) & ".");
        end if;
        exit when (integer32(myjob) > m);
      end loop;
    end do_next_job;
    procedure do_jobs is new Multitasking.Reporting_Workers(do_next_job);

  begin
    do_jobs(n);
  end Start_Workers;

  procedure Main is

    n,m : integer32 := 0;

  begin
    put("Give the number of workers : "); get(n);
    put("Give the number of jobs : "); get(m);
    Start_Workers(n,m);
  end Main;

begin
  Main;
end ts_jobs;
