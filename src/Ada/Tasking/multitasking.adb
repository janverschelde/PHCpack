with text_io;                            use text_io;

package body Multitasking is

  function all_true ( n : integer32; b : boolean_array ) return boolean is
  begin
    for i in 1..n loop
      if not b(i)
       then return false;
      end if;
    end loop;
    return true;
  end all_true;

  function all_false ( n : integer32; b : boolean_array ) return boolean is
  begin
    for i in 1..n loop
      if b(i)
       then return false;
      end if;
    end loop;
    return true;
  end all_false;

  function to_string ( n : natural32 ) return string is
  begin
    case n is
      when 0 => return "0";
      when 1 => return "1";
      when 2 => return "2";
      when 3 => return "3";
      when 4 => return "4";
      when 5 => return "5";
      when 6 => return "6";
      when 7 => return "7";
      when 8 => return "8";
      when 9 => return "9";
      when others => return to_string(n/10) & to_string(n mod 10);
    end case;
  end to_string;

  function to_string ( n : integer32 ) return string is
  begin
    return to_string(natural32(n));
  end to_string;

  procedure show ( n : in integer32; b : in boolean_array; s : in string ) is
  begin
    for i in 1..n loop
      if b(i)
       then put_line("-> " & s & " for task " & to_string(natural32(i)) 
                           & " is true");
       else put_line("-> " & s & " for task " & to_string(natural32(i))
                           & " is false");
      end if;
    end loop;
  end show;

-- WORKERS executing one job :

  procedure Silent_Workers ( n : in integer32 ) is

    task type Worker ( id,n : integer32 );

    task body Worker is

    -- i : constant string := to_string(natural32(id));

    begin
      Job(id,n);
    end Worker; 

    procedure Launch_Workers ( i,n : in integer32 ) is

    -- DESCRIPTION :
    --   A recursive procedure to start n tasks, starting at i.

    -- REQUIRED : i > 0.

      w : Worker(i,n);

    begin
      if i < n
       then Launch_Workers(i+1,n);
      end if;
    end Launch_Workers;

  begin
    Launch_Workers(1,n);
  end Silent_Workers;

  procedure Reporting_Workers ( n : in integer32 ) is

    task type Worker ( id,n : integer32 );

    task body Worker is

      i : constant string := to_string(natural32(id));

    begin
      put_line("Worker " & i & " will get busy ...");
      Job(id,n);
      put_line("... worker " & i & " returns from business.");
    end Worker; 

    procedure Launch_Workers ( i,n : in integer32 ) is

    -- DESCRIPTION :
    --   A recursive procedure to start n tasks, starting at i.

    -- REQUIRED : i > 0.

      w : Worker(i,n);
      id : constant string := to_string(natural32(i));
      ns : constant string := to_string(natural32(n));

    begin
      if i >= n then
        put_line("All " & ns & " workers have been launched.");
      else
        put_line("Launched worker " & id & " ...");
        Launch_Workers(i+1,n);
      end if;
    end Launch_Workers;

  begin
    Launch_Workers(1,n);
  end Reporting_Workers;

-- WORKERS executing jobs in a while loop :

  procedure Silent_Looping_Workers ( n : in integer32 ) is

    busy : boolean_array(1..n) := (1..n => true);
    done : boolean_array(1..n) := (1..n => false);
    jobcnt : natural32 := 0;

    task type Worker ( id,n : integer32 );

    task body Worker is

      -- i : constant string := to_string(natural32(id));
      c : boolean := true;
      cnt : natural32 := 0;
       
    begin
      while c or ((id = 1) and not all_true(n,done)) loop
        if (id = 1 and not done(id)) or (id > 1) then
          Job(id,n,c);
        end if;
        if not c then          -- id = 1 regulates and must wait to quit
          done(id) := true;
          if (id = 1 and all_true(n,done)) or (id > 1)
           then busy(id) := false; exit;
          end if;
        end if;
        busy(id) := false;
        cnt := cnt + 1;
        if id = 1 then
          loop
            exit when all_false(n,busy);
          end loop;
          for j in 1..n loop  -- tasks that quit do not get busy
            if not done(j)
             then busy(j) := true;
            end if;
          end loop;
          jobcnt := jobcnt + 1; 
        else
          loop
            exit when (jobcnt = cnt);
          end loop;
        end if;
      end loop;
    end Worker; 

    procedure Launch_Workers ( i,n : in integer32 ) is

    -- DESCRIPTION :
    --   A recursive procedure to start n tasks, starting at i.

    -- REQUIRED : i > 0.

      w : Worker(i,n);

    begin
      if i < n
       then Launch_Workers(i+1,n);
      end if;
    end Launch_Workers;

  begin
    Launch_Workers(1,n);
  end Silent_Looping_Workers;

  procedure Reporting_Looping_Workers ( n : in integer32 ) is

     busy : boolean_array(1..n) := (1..n => true);
     done : boolean_array(1..n) := (1..n => false);
     jobcnt : natural32 := 0;

     task type Worker ( id,n : integer32 );

     task body Worker is

       i : constant string := to_string(natural32(id));
       c : boolean := true;
       cnt : natural32 := 0;

     begin
       put_line("Worker " & i & " will get busy ...");
       while c or ((id = 1) and not all_true(n,done)) loop
         if ((id = 1) and not done(id)) or (id > 1) then
           Job(id,n,c);
           cnt := cnt + 1;
           put_line("worker " & i & " has done iteration " & to_string(cnt));
         end if;
         if not c then
           show(n,done,"done");
           done(id) := true;
           if (id = 1 and all_true(n,done)) or (id > 1) then
             put_line("worker " & i & " quits");
             busy(id) := false;
             exit;
           end if;
         end if;
         busy(id) := false;
         if id = 1 then
           put_line("worker 1 waits for all tasks to finish ...");
           show(n,busy,"busy");
           loop
             exit when all_false(n,busy);
           end loop;
           for j in 1..n loop   -- tasks that quit do not become busy
             if not done(j) 
              then busy(j) := true;
             end if;
           end loop;
           jobcnt := jobcnt + 1;
           put_line("worker 1 updates job counter to " & to_string(jobcnt));
         else
           put_line("worker " & i & " waits for update of job counter ...");
           loop
             exit when (jobcnt = cnt);
           end loop;
           put_line("worker " & i & " finished waiting for job counter");
         end if;
       end loop;
       put_line("... worker " & i & " returns from business.");
     end Worker; 

     procedure Launch_Workers ( i,n : in integer32 ) is

     -- DESCRIPTION :
     --   A recursive procedure to start n tasks, starting at i.

     -- REQUIRED : i > 0.

       w : Worker(i,n);
       id : constant string := to_string(natural32(i));
       ns : constant string := to_string(natural32(n));

     begin
       if i >= n then
         put_line("All " & ns & " workers have been launched.");
       else
         put_line("Launched worker " & id & " ...");
         Launch_Workers(i+1,n);
       end if;
     end Launch_Workers;

  begin
    Launch_Workers(1,n);
  end Reporting_Looping_Workers;

end Multitasking;
