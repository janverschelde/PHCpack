with text_io;
with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;       use Standard_Natural_Numbers_io;
with Multitasking,Semaphore;

procedure ts_mutex is

-- DESCRIPTION :
--   Exploration of critical sections using semaphores.
--   Two tasks take turns incrementing a counter.
--   It does not work when the tasks do input/output in
--   their critical sections!

  use Multitasking; -- for the routine to_string

  procedure Main ( max : in natural32 ) is

    s : Semaphore.Lock;
    turn : natural32 := 1;
    cnt : natural32 := 0;
    task t1;
    task t2;
   
    task body t1 is
    begin
      text_io.put_line("Task One says hello, count = " & to_string(cnt));
      loop
        if turn = 1 then
          text_io.put_line("Task One is requesting, count = "
            & to_string(cnt));
          Semaphore.Request(s);
         -- text_io.put_line("Task One requested, value = "
         --    & to_string(Semaphore.Value(s)));
          cnt := cnt + 1; turn := 2;
         -- text_io.put_line("Task One increased, count = " 
         --    & to_string(cnt) & ", turn = " & to_string(turn));
         -- text_io.put_line("Task One about to release, value = "
         --    & to_string(Semaphore.Value(s)));
          Semaphore.Release(s);
          text_io.put_line("Task One released, count = "
             & to_string(cnt));
        end if;
        exit when (cnt > max);
      end loop;
    end t1;

    task body t2 is
    begin
      text_io.put_line("Task Two says hello, count = " & to_string(cnt));
      loop
        if turn = 2 then
          text_io.put_line("Task Two is requesting, count = "
            & to_string(cnt));
          Semaphore.Request(s);
         -- text_io.put_line("Task Two requested, value = "
         --    & to_string(Semaphore.Value(s)));
          cnt := cnt + 1; turn := 1;
         -- text_io.put_line("Task Two increased, count = " 
         --   & to_string(cnt) & ", turn = " & to_string(turn));
         -- text_io.put_line("Task Two about to release, value = "
         --    & to_string(Semaphore.Value(s)));
          Semaphore.Release(s);
          text_io.put_line("Task Two released, count = "
             & to_string(cnt));
        end if;
        exit when (cnt > max);
      end loop;
    end t2;
   
  begin
    text_io.put_line("Launched two tasks flipping turns.");
  end Main;

begin
  declare
    n : natural32 := 0;
  begin
    text_io.put("Give number of turns : ");
    get(n);
    Main(n);
  end;
end ts_mutex;
