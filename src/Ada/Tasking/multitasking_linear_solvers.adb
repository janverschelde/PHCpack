with text_io;                            use text_io;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Double_Double_Numbers;              use Double_Double_Numbers;
with Quad_Double_Numbers;                use Quad_Double_Numbers;
with Standard_Complex_Numbers;
with DoblDobl_Complex_Numbers;
with QuadDobl_Complex_Numbers;
with Multitasking,Semaphore;             use Multitasking;

package body multitasking_linear_solvers is

-- AUXILIARY for lufac :

  function cabs ( c : Standard_Complex_Numbers.Complex_Number )
                return double_float is

    use Standard_Complex_Numbers;

  begin
    return (ABS(REAL_PART(c)) + ABS(IMAG_PART(c)));
  end cabs;

  function cabs ( c : DoblDobl_Complex_Numbers.Complex_Number )
                return double_double is

    use DoblDobl_Complex_Numbers;

  begin
    return (ABS(REAL_PART(c)) + ABS(IMAG_PART(c)));
  end cabs;

  function cabs ( c : QuadDobl_Complex_Numbers.Complex_Number )
                return quad_double is

    use QuadDobl_Complex_Numbers;

  begin
    return (ABS(REAL_PART(c)) + ABS(IMAG_PART(c)));
  end cabs;

-- TARGET PROCEDURES for STANDARD complex :

  procedure silent_lufac
              ( t : in integer32;
                a : in out Standard_Complex_Matrices.Matrix;
                n : in integer32;
                ipvt : out Standard_Integer_Vectors.Vector;
                info : out integer32 ) is

    nm1 : constant integer32 := n-1;
    go_pivot : integer32 := 0;
    i_row,j_col : integer32;
    mult_done : boolean := false;
    elim_done : boolean := false;
    use Standard_Complex_Numbers;
    temp : Complex_Number;
    sem : Semaphore.Lock;
    k : integer32 := 1;
    ll : integer32;

    procedure do_job ( id,nt : in integer32 ) is

      smax : double_float;
      kp1,my_row,my_col : integer32;

    begin
      if nm1 >= 1 then
        while k <= nm1 loop
          kp1 := k + 1;
          while go_pivot < k loop
            Semaphore.Request(sem);               -- find the pivot index ll
            if go_pivot < k then
              ll := k; smax := cabs(a(k,k));
              for i in kp1..n loop
                if cabs(a(i,k)) > smax then
                  ll := i;
                  smax := cabs(a(i,k));
                end if;
              end loop;
              ipvt(k) := ll;
              if smax = 0.0 then    -- this column is already triangularized
                info := k;
              else
                if ll /= k then                  -- interchange if necessary
                  temp := a(ll,k);
                  a(ll,k) := a(k,k);
                  a(k,k) := temp;
                end if;
                temp := -Create(1.0);                 -- compute multipliers
                temp := temp/a(k,k);
              end if;
              j_col := kp1;  -- for j in kp1..n loop in the serial algorithm
              i_row := kp1;  -- counter in the j-loop
              mult_done := false;
              elim_done := false;
              go_pivot := k;
            end if;
            Semaphore.Release(sem);
          end loop;
          while not mult_done loop    -- compute multipliers
            Semaphore.Request(sem);
            if not mult_done then
              if i_row = n then         -- last update in critical section
                a(i_row,k) := temp*a(i_row,k);
                i_row := k;             -- for elimination stage
                mult_done := true;
              else
                my_row := i_row;
                i_row := i_row + 1;     -- update job counter
              end if;
            end if;
            Semaphore.Release(sem);
            exit when mult_done;
            a(my_row,k) := temp*a(my_row,k);
          end loop;
          while not elim_done loop   -- row elimination with column indexing
            Semaphore.Request(sem);
            if not elim_done then
              if i_row = k then        -- swap job in critical section
                temp := a(ll,j_col);   -- set the multiplier to temp
                if ll /= k then        -- swap if necessary
                  a(ll,j_col) := a(k,j_col);
                  a(k,j_col) := temp;
                end if;
                my_row := i_row;       -- swapper will skip a step
              elsif i_row = n then     -- the last update with the temp
                a(i_row,j_col) := a(i_row,j_col) + temp*a(i_row,k);
                if j_col = n then
                  k := k + 1;
                  elim_done := true;
                else
                  my_row := k;     -- make sure to skip a step
                end if;
              else
                my_row := i_row;       -- grab next update job
                my_col := j_col;
              end if;
              if not elim_done then
                i_row := i_row + 1;    -- update current job counter
                if i_row > n then
                  j_col := j_col + 1;
                  i_row := k;
                end if;
              end if;
            end if;
            Semaphore.Release(sem);
            exit when elim_done;
            if (my_row /= k)
             then a(my_row,my_col) := a(my_row,my_col) + temp*a(my_row,k);
            end if;
          end loop;
        end loop;
      end if;
      if id = 1 then
        ipvt(n) := n;
        if AbsVal(a(n,n)) = 0.0
         then info := n;
        end if;
      end if;
    end do_job;
    procedure do_jobs is new Silent_Workers(do_job);

  begin
    info := 0;
    if n > t
     then do_jobs(t);
     else do_jobs(n);
    end if;
  end silent_lufac;

  function Different
             ( a,b : in Standard_Complex_Numbers.Complex_Number )
             return boolean is

  -- DESCRIPTION :
  --   Returns true if a and b different by more than 1.0E-8.

    use Standard_Complex_Numbers;
    d : constant Complex_Number := a - b;
    tol : constant double_float := 1.0E-8;

  begin
    if AbsVal(d) > tol
     then return true;
     else return false;
    end if;
  end Different;

  procedure reporting_lufac
              ( t : in integer32;
                lu : in Standard_Complex_Matrices.Matrix;
                a : in out Standard_Complex_Matrices.Matrix;
                n : in integer32;
                ipvt : out Standard_Integer_Vectors.Vector;
                info : out integer32 ) is

    nm1 : constant integer32 := n-1;
    go_pivot : integer32 := 0;
    mult_done : boolean := false;
    elim_done : boolean := false;
    i_row,j_col : integer32;
    use Standard_Complex_Numbers;
    temp : Complex_Number;
    sem : Semaphore.Lock;
    k : integer32 := 1;
    ll : integer32;

    procedure do_job ( id,nt : integer32 ) is

      smax : double_float;
      kp1,my_row,my_col : integer32;

    begin
      if nm1 >= 1 then
        while k <= nm1 loop
          kp1 := k + 1;
          while go_pivot < k loop
            Semaphore.Request(sem);               -- find the pivot index ll
            if go_pivot < k then
              put_line("*** task " & to_string(id)
                                   & " computes pivot at column "
                                   & to_string(k));
              ll := k; smax := cabs(a(k,k));
              for i in kp1..n loop
                if cabs(a(i,k)) > smax then
                  ll := i;
                  smax := cabs(a(i,k));
                end if;
              end loop;
              ipvt(k) := ll;
              if smax = 0.0 then    -- this column is already triangularized
                info := k;
              else
                if ll /= k then                  -- interchange if necessary
                  temp := a(ll,k);
                  a(ll,k) := a(k,k);
                  a(k,k) := temp;
                end if;
                temp := -Create(1.0)/a(k,k);          -- compute multipliers
              end if;
              i_row := kp1;  -- double indexing for the job counter
              j_col := kp1;  -- for j in kp1..n and i in kp1..n
              mult_done := false;
              elim_done := false;
              go_pivot := k;
              put_line("task " & to_string(id) & " sets go_pivot to "
                               & to_string(go_pivot));
            end if;
            Semaphore.Release(sem);
          end loop;
          while not mult_done loop
            Semaphore.Request(sem);
            if not mult_done then
              if i_row = n then
                put_line("** task " & to_string(id) & " last multiplies row "
                                    & to_string(i_row));
                a(i_row,k) := temp*a(i_row,k);
                i_row := k;          -- for elimination stage
                mult_done := true;
              else
                my_row := i_row;
                i_row := i_row + 1;
              end if;
            end if;
            Semaphore.Release(sem);
            exit when mult_done;
            put_line("task " & to_string(id) & " multiplies "
                             & to_string(my_row));
            a(my_row,k) := temp*a(my_row,k);
          end loop;
          while not elim_done loop
            Semaphore.Request(sem);
            if not elim_done then
              put_line("task " & to_string(id)
                               & " requests at i_row = " & to_string(i_row)
                               & " and at j_col = " & to_string(j_col));
              if i_row = k then
                put_line("task " & to_string(id)
                                 & " sets multiplier and swaps");
                temp := a(ll,j_col);    -- set multiplier to temp
                if ll /= k then         -- swap if necessary
                  a(ll,j_col) := a(k,j_col);
                  a(k,j_col) := temp;
                end if;
                my_row := i_row;       -- swapper will skip next job
              elsif i_row = n then
                put_line("** task " & to_string(id) & " does last update");
                a(i_row,j_col) := a(i_row,j_col) + temp*a(i_row,k);
                if j_col = n then
                  put_line("*** update done at column " & to_string(k));
                  k := k+1;
                  elim_done := true;   -- no next job
                else
                  my_row := k;        -- make sure to skip a job
                end if;
              else
                put_line("task " & to_string(id) & " grabs next job");
                my_row := i_row;    -- next update out of critical section
                my_col := j_col;
              end if;
              if not elim_done then
                i_row := i_row + 1;    -- update current job counter
                if i_row > n then
                  j_col := j_col + 1;
                  i_row := k;
                end if;
              end if;
            end if;
            Semaphore.Release(sem);
            exit when elim_done;
            if (my_row /= k) then
              put_line("task " & to_string(id)
                               & " updates row " & to_string(my_row)
                               & " updates column " & to_string(my_col));
              a(my_row,my_col) := a(my_row,my_col) + temp*a(my_row,k);
            end if;
          end loop;
          if k > 1 then
            for i in 1..n loop
              if Different(lu(i,k-1),a(i,k-1)) then
                put_line("BUG at row " & to_string(i)
                         & " and column " & to_string(k-1)); -- return;
              end if;
            end loop;
          end if;
        end loop;
      end if;
      if id = 1 then
        ipvt(n) := n;
        if AbsVal(a(n,n)) = 0.0
         then info := n;
        end if;
      end if;
    end do_job;
    procedure do_jobs is new Reporting_Workers(do_job);

  begin
    info := 0;
    if n > t
     then do_jobs(t);
     else do_jobs(n);
    end if;
  end reporting_lufac;

  procedure silent_lusolve
              ( t : in integer32;
                a : in Standard_Complex_Matrices.Matrix;
                n : in integer32;
                ipvt : in Standard_Integer_Vectors.Vector;
                b : in out Standard_Complex_Vectors.Vector ) is

    nm1 : constant integer32 := n-1;
    go_forward : integer32 := 0;
    go_back : integer32 := n+1;
    busy : boolean_array(1..t) := (1..t => false);
    done : boolean_array(1..t) := (1..t => false);
    nt : integer32;
    use Standard_Complex_Numbers;
    temp : Complex_Number;

    procedure do_job ( id,nt : in integer32 ) is

      ll,kb : integer32;
 
    begin
      if nm1 >= 1 then                              -- solve L*y = b
        for k in 1..nm1 loop
          if id = 1 then
            loop
              exit when all_false(nt,busy); 
            end loop;
            ll := ipvt(k);
            temp := b(ll);
            if ll /= k then
              b(ll) := b(k);
              b(k) := temp;
            end if;
            busy := (1..nt => true);
            go_forward := k;
          else
            loop
              exit when (go_forward = k);
            end loop;
          end if;
          for i in (k+1)..n loop
            if i mod nt = id - 1 then
              b(i) := b(i) + temp*a(i,k);
              --if i = ipvt(k)
              -- then b(k) := b(k) + b(ipvt(k))*a(k,k);
              -- else b(i) := b(i) + b(ipvt(k))*a(i,k);
              --end if;
            end if;
          end loop;
          busy(id) := false;
        end loop;
      end if;
      done(id) := true;
      loop
        exit when all_true(nt,done);
      end loop;
      for k in 1..n loop                            -- solve U*x = y
        kb := n+1-k;
        if id = 1 then
          loop
            exit when all_false(nt,busy);
          end loop;
          b(kb) := b(kb) / a(kb,kb);
          busy(1..nt) := (1..nt => true);
          go_back := kb;
        else
          loop
            exit when (go_back = kb);
          end loop;
        end if;
        for j in 1..(kb-1) loop
          if j mod nt = id - 1 then
            b(j) := b(j) - b(kb)*a(j,kb);
          end if;
        end loop;
        busy(id) := false;
      end loop;
    end do_job;
    procedure do_jobs is new Silent_Workers(do_job);

  begin
    if t > n
     then nt := n; do_jobs(n);
     else nt := t; do_jobs(t);
    end if;
  end silent_lusolve;

  procedure reporting_lusolve
              ( t : in integer32;
                a : in Standard_Complex_Matrices.Matrix;
                n : in integer32;
                ipvt : in Standard_Integer_Vectors.Vector;
                b : in out Standard_Complex_Vectors.Vector ) is

    nm1 : constant integer32 := n-1;
    go_forward : integer32 := 0;
    go_back : integer32 := n+1;
    busy : boolean_array(1..t) := (1..t => false);
    done : boolean_array(1..t) := (1..t => false);
    nt : integer32;
    use Standard_Complex_Numbers;
    temp : Complex_Number;

    procedure do_job ( id,nt : in integer32 ) is

      ll,kb : integer32;
 
    begin
      if nm1 >= 1 then                                  -- solve L*y = b
        for k in 1..nm1 loop
          if id = 1 then
            show(nt,busy,"busy");
            put_line("task 1 waits till none are busy ...");
            loop
              exit when all_false(nt,busy); 
            end loop;
            put_line("task 1 pivots at k = " & to_string(k));
            ll := ipvt(k);
            temp := b(ll);
            if ll /= k then
              b(ll) := b(k);
              b(k) := temp;
            end if;
            busy(1..nt) := (1..nt => true);
            go_forward := k;
            put_line("**** task 1 sets go to " & to_string(go_forward));
          else
            put_line("task " & to_string(id) & " waits to update"
                  & " at go = " & to_string(go_forward));
            loop
              exit when (go_forward = k);
            end loop;
            put_line("task " & to_string(id) & " finished wait for go");
          end if;
          for i in (k+1)..n loop
            if i mod nt = id - 1 then
              put_line("task " & to_string(id) & " updates " & to_string(i));
              b(i) := b(i) + temp*a(i,k);
            end if;
          end loop;
          busy(id) := false;
          put_line("task " & to_string(id) & " done at k = " & to_string(k));
        end loop;
      end if;
      done(id) := true;
      put_line("task " & to_string(id) & " is done with L*y = b");
      loop
        exit when all_true(nt,done);
      end loop;
      for k in 1..n loop                            -- solve U*x = y
        kb := n+1-k;
        if id = 1 then
          show(nt,busy,"busy");
          put_line("task 1 waits till none are busy ...");
          loop
            exit when all_false(nt,busy);
          end loop;
          put_line("task 1 sets temp at kb = " & to_string(kb));
          b(kb) := b(kb) / a(kb,kb);
          put_line("**** task 1 sets go to " & to_string(go_back));
          busy := (1..nt => true);
          go_back := kb;
        else
          put_line("task " & to_string(id) & " waits to update"
                & " at go = " & to_string(go_back));
          loop
            exit when (go_back = kb);
          end loop;
          put_line("task " & to_string(id) & " finished wait for go");
        end if;
        for j in 1..(kb-1) loop
          if j mod nt = id - 1 then
            put_line("task " & to_string(id) & " updates " & to_string(j));
            b(j) := b(j) - b(kb)*a(j,kb);
          end if;
        end loop;
        busy(id) := false;
        put_line("task " & to_string(id) & " done at kb = " & to_string(kb));
      end loop;
    end do_job;
    procedure do_jobs is new Reporting_Workers(do_job);

  begin
    if t > n
     then nt := n; do_jobs(n);
     else nt := t; do_jobs(t);
    end if;
  end reporting_lusolve;

-- TARGET PROCEDURES for DOUBLE DOUBLE complex :

  procedure silent_lufac
              ( t : in integer32;
                a : in out DoblDobl_Complex_Matrices.Matrix;
                n : in integer32;
                ipvt : out Standard_Integer_Vectors.Vector;
                info : out integer32 ) is

    nm1 : constant integer32 := n-1;
    go_pivot : integer32 := 0;
    ll,i_row,j_col : integer32;
    mult_done : boolean := false;
    elim_done : boolean := false;
    zero : constant double_double := Create(0.0);
    use DoblDobl_Complex_Numbers;
    temp : Complex_Number;
    sem : Semaphore.Lock;
    k : integer32 := 1;

    procedure do_job ( id,nt : in integer32 ) is

      kp1,my_row,my_col : integer32;
      smax : double_double;

    begin
      if nm1 >= 1 then
        while k <= nm1 loop
          kp1 := k + 1;   
          while go_pivot < k loop                -- find the pivot index ll
            Semaphore.Request(sem);
            if go_pivot < k then
              ll := k; smax := cabs(a(k,k));
              for i in kp1..n loop
                if cabs(a(i,k)) > smax then
                  ll := i;
                  smax := cabs(a(i,k));
                end if;
              end loop;
              ipvt(k) := ll;
              if smax = zero then  -- this column is already triangularized
                info := k;
              else
                if ll /= k then                 -- interchange if necessary
                  temp := a(ll,k);
                  a(ll,k) := a(k,k);
                  a(k,k) := temp;
                end if;
                temp := -Create(double_double_numbers.create(1.0));
                temp := temp/a(k,k);
              end if;
              i_row := kp1; -- initialization of the double job index 
              j_col := kp1; -- for j in kp1..n in the serial algorithm
              mult_done := false;
              elim_done := false;
              go_pivot := k;
            end if;
            Semaphore.Release(sem);
          end loop;
          while not mult_done loop      -- compute multipliers
            Semaphore.Request(sem);
            if not mult_done then
              if i_row = n then         -- last update in critical section
                a(i_row,k) := temp*a(i_row,k);
                i_row := k;             -- for elimination stage
                mult_done := true;
              else
                my_row := i_row;
                i_row := i_row + 1;     -- update job counter
              end if;
            end if;
            Semaphore.Release(sem);
            exit when mult_done;
            a(my_row,k) := temp*a(my_row,k);
          end loop;
          while not elim_done loop  -- row elimination with column indexing
            Semaphore.Request(sem);
            if not elim_done then
              if i_row = k then        -- swap job in critical section
                temp := a(ll,j_col);   -- set the multiplier to temp
                if ll /= k then        -- swap if necessary
                  a(ll,j_col) := a(k,j_col);
                  a(k,j_col) := temp;
                end if;
                my_row := i_row;       -- swapper will skip a step
              elsif i_row = n then     -- the last update with the temp
                a(i_row,j_col) := a(i_row,j_col) + temp*a(i_row,k);
                if j_col = n then
                  k := k + 1;
                  elim_done := true;
                else
                  my_row := k;     -- make sure to skip a step
                end if;
              else
                my_row := i_row;       -- grab next update job
                my_col := j_col;
              end if;
              if not elim_done then
                i_row := i_row + 1;    -- update current job counter
                if i_row > n then
                  j_col := j_col + 1;
                  i_row := k;
                end if;
              end if;
            end if;
            Semaphore.Release(sem);
            exit when elim_done;
            if (my_row /= k)
             then a(my_row,my_col) := a(my_row,my_col) + temp*a(my_row,k);
            end if;
          end loop;
        end loop;
      end if;
      if id = 1 then
        ipvt(n) := n;
        if AbsVal(a(n,n)) = zero
         then info := n;
        end if;
      end if;
    end do_job;
    procedure do_jobs is new Silent_Workers(do_job);

  begin
    info := 0;
    if n > t
     then do_jobs(t);
     else do_jobs(n);
    end if;
  end silent_lufac;

  procedure reporting_lufac
              ( t : in integer32;
                a : in out DoblDobl_Complex_Matrices.Matrix;
                n : in integer32;
                ipvt : out Standard_Integer_Vectors.Vector;
                info : out integer32 ) is

    nm1 : constant integer32 := n-1;
    go : integer32 := 0;
    go_elim : integer32 := 0;
    busy : boolean_array(1..t) := (1..t => false);
    nt : integer32;
    zero : constant double_double := Create(0.0);
    use DoblDobl_Complex_Numbers;
    temp : Complex_Number;

    procedure do_job ( id,nt : in integer32 ) is

      kp1,ll : integer32;
      smax : double_double;

    begin
      if nm1 >= 1 then
        for k in 1..nm1 loop
          kp1 := k + 1;                            -- find the pivot index ll
          if id = 1 then
            show(nt,busy,"busy");
            put_line("task 1 waits till none are busy to find pivot...");
            loop
              exit when all_false(nt,busy);
            end loop;
            put_line("task 1 finds pivot at k = " & to_string(k));
            ll := k; smax := cabs(a(k,k));
            for i in kp1..n loop
              if cabs(a(i,k)) > smax then
                ll := i;
                smax := cabs(a(i,k));
              end if;
            end loop;
            ipvt(k) := ll;
            if smax = zero then      -- this column is already triangularized
              info := k;
            else
              if ll /= k then                     -- interchange if necessary
                temp := a(ll,k);
                a(ll,k) := a(k,k);
                a(k,k) := temp;
              end if;
             -- temp := -Create(1.0)/a(k,k);           -- compute multipliers
              temp := -Create(double_double_numbers.create(1.0));
              temp := temp/a(k,k);
            end if;
            busy(1..nt) := (1..nt => true);
            go := k;
            put_line("**** task 1 sets go to " & to_string(go));
          else
            put_line("task " & to_string(id) & " waits to update"
                  & " at go = " & to_string(go));
            loop
              exit when (go = k);
            end loop;
            put_line("task " & to_string(id) & " finished wait for go");
          end if;
          for i in kp1..n loop
            if i mod nt = id - 1 then
              put_line("task " & to_string(id) & " updates " & to_string(i));
              a(i,k) := temp*a(i,k);
            end if;
          end loop;
          busy(id) := false;
          for j in kp1..n loop       -- row elimination with column indexing
            if id = 1 then
              show(nt,busy,"busy");
              put_line("task 1 waits till none are busy to swap...");
              loop
                exit when all_false(nt,busy);
              end loop;
              put_line("task 1 swaps at k = " & to_string(k));
              temp := a(ll,j);
              if ll /= k then
                a(ll,j) := a(k,j);
                a(k,j) := temp;
              end if;
              busy(1..nt) := (1..nt => true);
              go_elim := k*(n-1) + j;
              put_line("**** task 1 sets go_elim to " & to_string(go_elim));
            else
              put_line("task " & to_string(id) & " waits to update"
                    & " at go_elim = " & to_string(go_elim));
              loop
                exit when (go_elim = (k*(n-1)+j));
              end loop;
              put_line("task " & to_string(id) & " finished wait for go");
            end if;
            for i in kp1..n loop
              if i mod nt = id - 1 then
                put_line("task " & to_string(id) & " updates " & to_string(i));
                a(i,j) := a(i,j) + temp*a(i,k);
              end if;
            end loop;
            busy(id) := false;
          end loop;
        end loop;
      end if;
      if id = 1 then
        ipvt(n) := n;
        if AbsVal(a(n,n)) = zero
         then info := n;
        end if;
      end if;
    end do_job;
    procedure do_jobs is new Reporting_Workers(do_job);

  begin
    info := 0;
    if n > t
     then nt := t; do_jobs(t);
     else nt := n; do_jobs(n);
    end if;
  end reporting_lufac;

  procedure silent_lusolve
              ( t : in integer32;
                a : in DoblDobl_Complex_Matrices.Matrix;
                n : in integer32;
                ipvt : in Standard_Integer_Vectors.Vector;
                b : in out DoblDobl_Complex_Vectors.Vector ) is

    nm1 : constant integer32 := n-1;
    go_forward : integer32 := 0;
    go_back : integer32 := n+1;
    busy : boolean_array(1..t) := (1..t => false);
    done : boolean_array(1..t) := (1..t => false);
    nt : integer32;
    use DoblDobl_Complex_Numbers;
    temp : Complex_Number;

    procedure do_job ( id,nt : in integer32 ) is

      ll,kb : integer32;
 
    begin
      if nm1 >= 1 then                              -- solve L*y = b
        for k in 1..nm1 loop
          if id = 1 then
            loop
              exit when all_false(nt,busy); 
            end loop;
            ll := ipvt(k);
            temp := b(ll);
            if ll /= k then
              b(ll) := b(k);
              b(k) := temp;
            end if;
            busy := (1..nt => true);
            go_forward := k;
          else
            loop
              exit when (go_forward = k);
            end loop;
          end if;
          for i in (k+1)..n loop
            if i mod nt = id - 1 then
              b(i) := b(i) + temp*a(i,k);
            end if;
          end loop;
          busy(id) := false;
        end loop;
      end if;
      done(id) := true;
      loop
        exit when all_true(nt,done);
      end loop;
      for k in 1..n loop                            -- solve U*x = y
        kb := n+1-k;
        if id = 1 then
          loop
            exit when all_false(nt,busy);
          end loop;
          b(kb) := b(kb) / a(kb,kb);
          busy := (1..nt => true);
          go_back := kb;
        else
          loop
            exit when (go_back = kb);
          end loop;
        end if;
        for j in 1..(kb-1) loop
          if j mod nt = id - 1 then
            b(j) := b(j) - b(kb)*a(j,kb);
          end if;
        end loop;
        busy(id) := false;
      end loop;
    end do_job;
    procedure do_jobs is new Silent_Workers(do_job);

  begin
    if t > n
     then nt := n; do_jobs(n);
     else nt := t; do_jobs(t);
    end if;
  end silent_lusolve;

  procedure reporting_lusolve
              ( t : in integer32;
                a : in DoblDobl_Complex_Matrices.Matrix;
                n : in integer32;
                ipvt : in Standard_Integer_Vectors.Vector;
                b : in out DoblDobl_Complex_Vectors.Vector ) is

    nm1 : constant integer32 := n-1;
    go_forward : integer32 := 0;
    go_back : integer32 := n+1;
    busy : boolean_array(1..t) := (1..t => false);
    done : boolean_array(1..t) := (1..t => false);
    nt : integer32;
    use DoblDobl_Complex_Numbers;
    temp : Complex_Number;

    procedure do_job ( id,nt : in integer32 ) is

      ll,kb : integer32;
 
    begin
      if nm1 >= 1 then                                  -- solve L*y = b
        for k in 1..nm1 loop
          if id = 1 then
            show(nt,busy,"busy");
            put_line("task 1 waits till none are busy ...");
            loop
              exit when all_false(nt,busy); 
            end loop;
            put_line("task 1 pivots at k = " & to_string(k));
            ll := ipvt(k);
            temp := b(ll);
            if ll /= k then
              b(ll) := b(k);
              b(k) := temp;
            end if;
            busy := (1..nt => true);
            go_forward := k;
            put_line("**** task 1 sets go to " & to_string(go_forward));
          else
            put_line("task " & to_string(id) & " waits to update"
                  & " at go = " & to_string(go_forward));
            loop
              exit when (go_forward = k);
            end loop;
            put_line("task " & to_string(id) & " finished wait for go");
          end if;
          for i in (k+1)..n loop
            if i mod nt = id - 1 then
              put_line("task " & to_string(id) & " updates " & to_string(i));
              b(i) := b(i) + temp*a(i,k);
            end if;
          end loop;
          busy(id) := false;
          put_line("task " & to_string(id) & " done at k = " & to_string(k));
        end loop;
      end if;
      done(id) := true;
      put_line("task " & to_string(id) & " is done with L*y = b");
      loop
        exit when all_true(nt,done);
      end loop;
      for k in 1..n loop                            -- solve U*x = y
        kb := n+1-k;
        if id = 1 then
          show(nt,busy,"busy");
          put_line("task 1 waits till none are busy ...");
          loop
            exit when all_false(nt,busy);
          end loop;
          put_line("task 1 sets temp at kb = " & to_string(kb));
          b(kb) := b(kb) / a(kb,kb);
          put_line("**** task 1 sets go to " & to_string(go_back));
          busy := (1..nt => true);
          go_back := kb;
        else
          put_line("task " & to_string(id) & " waits to update"
                & " at go = " & to_string(go_back));
          loop
            exit when (go_back = kb);
          end loop;
          put_line("task " & to_string(id) & " finished wait for go");
        end if;
        for j in 1..(kb-1) loop
          if j mod nt = id - 1 then
            put_line("task " & to_string(id) & " updates " & to_string(j));
            b(j) := b(j) - b(kb)*a(j,kb);
          end if;
        end loop;
        busy(id) := false;
        put_line("task " & to_string(id) & " done at kb = " & to_string(kb));
      end loop;
    end do_job;
    procedure do_jobs is new Reporting_Workers(do_job);

  begin
    if t > n
     then nt := n; do_jobs(n);
     else nt := t; do_jobs(t);
    end if;
  end reporting_lusolve;

-- TARGET PROCEDURES for QUAD DOUBLE complex :

  procedure silent_lufac
              ( t : in integer32;
                a : in out QuadDobl_Complex_Matrices.Matrix;
                n : in integer32;
                ipvt : out Standard_Integer_Vectors.Vector;
                info : out integer32 ) is

    nm1 : constant integer32 := n-1;
    go_pivot : integer32 := 0;
    ll,i_row,j_col : integer32;
    mult_done : boolean := false;
    elim_done : boolean := false;
    zero : constant quad_double := Create(0.0);
    use QuadDobl_Complex_Numbers;
    temp : Complex_Number;
    sem : Semaphore.Lock;
    k : integer32 := 1;

    procedure do_job ( id,nt : in integer32 ) is

      kp1,my_row,my_col : integer32;
      smax : quad_double;

    begin
      if nm1 >= 1 then
        while k <= nm1 loop
          kp1 := k + 1; 
          while go_pivot < k loop                 -- find the pivot index ll
            Semaphore.Request(sem);
            if go_pivot < k then
              ll := k; smax := cabs(a(k,k));
              for i in kp1..n loop
                if cabs(a(i,k)) > smax then
                  ll := i;
                  smax := cabs(a(i,k));
                end if;
              end loop;
              ipvt(k) := ll;
              if smax = zero then  -- this column is already triangularized
                info := k;
              else
                if ll /= k then                 -- interchange if necessary
                  temp := a(ll,k);
                  a(ll,k) := a(k,k);
                  a(k,k) := temp;
                end if;
                temp := -Create(quad_double_numbers.create(1.0));
                temp := temp/a(k,k);
              end if;
              i_row := kp1;      -- initialization of the double job index 
              j_col := kp1;      -- for j in kp1..n in the serial algorithm
              mult_done := false;
              elim_done := false;
              go_pivot := k;
            end if;
            Semaphore.Release(sem);
          end loop;
          while not mult_done loop      -- compute multipliers
            Semaphore.Request(sem);
            if not mult_done then
              if i_row = n then         -- last update in critical section
                a(i_row,k) := temp*a(i_row,k);
                i_row := k;             -- for elimination stage
                mult_done := true;
              else
                my_row := i_row;
                i_row := i_row + 1;     -- update job counter
              end if;
            end if;
            Semaphore.Release(sem);
            exit when mult_done;
            a(my_row,k) := temp*a(my_row,k);
          end loop;
          while not elim_done loop  -- row elimination with column indexing
            Semaphore.Request(sem);
            if not elim_done then
              if i_row = k then        -- swap job in critical section
                temp := a(ll,j_col);   -- set the multiplier to temp
                if ll /= k then        -- swap if necessary
                  a(ll,j_col) := a(k,j_col);
                  a(k,j_col) := temp;
                end if;
                my_row := i_row;       -- swapper will skip a step
              elsif i_row = n then     -- the last update with the temp
                a(i_row,j_col) := a(i_row,j_col) + temp*a(i_row,k);
                if j_col = n then
                  k := k + 1;
                  elim_done := true;
                else
                  my_row := k;     -- make sure to skip a step
                end if;
              else
                my_row := i_row;       -- grab next update job
                my_col := j_col;
              end if;
              if not elim_done then
                i_row := i_row + 1;    -- update current job counter
                if i_row > n then
                  j_col := j_col + 1;
                  i_row := k;
                end if;
              end if;
            end if;
            Semaphore.Release(sem);
            exit when elim_done;
            if (my_row /= k)
             then a(my_row,my_col) := a(my_row,my_col) + temp*a(my_row,k);
            end if;
          end loop;
        end loop;
      end if;
      if id = 1 then
        ipvt(n) := n;
        if AbsVal(a(n,n)) = zero
         then info := n;
        end if;
      end if;
    end do_job;
    procedure do_jobs is new Silent_Workers(do_job);

  begin
    info := 0;
    if n > t
     then do_jobs(t);
     else do_jobs(n);
    end if;
  end silent_lufac;

  procedure reporting_lufac
              ( t : in integer32;
                a : in out QuadDobl_Complex_Matrices.Matrix;
                n : in integer32;
                ipvt : out Standard_Integer_Vectors.Vector;
                info : out integer32 ) is

    nm1 : constant integer32 := n-1;
    go : integer32 := 0;
    go_elim : integer32 := 0;
    busy : boolean_array(1..t) := (1..t => false);
    nt : integer32;
    zero : constant quad_double := Create(0.0);
    use QuadDobl_Complex_Numbers;
    temp : Complex_Number;

    procedure do_job ( id,nt : in integer32 ) is

      kp1,ll : integer32;
      smax : quad_double;

    begin
      if nm1 >= 1 then
        for k in 1..nm1 loop
          kp1 := k + 1;                            -- find the pivot index ll
          if id = 1 then
            show(nt,busy,"busy");
            put_line("task 1 waits till none are busy to find pivot...");
            loop
              exit when all_false(nt,busy);
            end loop;
            put_line("task 1 finds pivot at k = " & to_string(k));
            ll := k; smax := cabs(a(k,k));
            for i in kp1..n loop
              if cabs(a(i,k)) > smax then
                ll := i;
                smax := cabs(a(i,k));
              end if;
            end loop;
            ipvt(k) := ll;
            if smax = zero then      -- this column is already triangularized
              info := k;
            else
              if ll /= k then                     -- interchange if necessary
                temp := a(ll,k);
                a(ll,k) := a(k,k);
                a(k,k) := temp;
              end if;
             -- temp := -Create(1.0)/a(k,k);           -- compute multipliers
              temp := -Create(quad_double_numbers.create(1.0));
              temp := temp/a(k,k);
            end if;
            busy(1..nt) := (1..nt => true);
            go := k;
            put_line("**** task 1 sets go to " & to_string(go));
          else
            put_line("task " & to_string(id) & " waits to update"
                  & " at go = " & to_string(go));
            loop
              exit when (go = k);
            end loop;
            put_line("task " & to_string(id) & " finished wait for go");
          end if;
          for i in kp1..n loop
            if i mod nt = id - 1 then
              put_line("task " & to_string(id) & " updates " & to_string(i));
              a(i,k) := temp*a(i,k);
            end if;
          end loop;
          busy(id) := false;
          for j in kp1..n loop       -- row elimination with column indexing
            if id = 1 then
              show(nt,busy,"busy");
              put_line("task 1 waits till none are busy to swap...");
              loop
                exit when all_false(nt,busy);
              end loop;
              put_line("task 1 swaps at k = " & to_string(k));
              temp := a(ll,j);
              if ll /= k then
                a(ll,j) := a(k,j);
                a(k,j) := temp;
              end if;
              busy(1..nt) := (1..nt => true);
              go_elim := k*(n-1) + j;
              put_line("**** task 1 sets go_elim to " & to_string(go_elim));
            else
              put_line("task " & to_string(id) & " waits to update"
                    & " at go_elim = " & to_string(go_elim));
              loop
                exit when (go_elim = (k*(n-1)+j));
              end loop;
              put_line("task " & to_string(id) & " finished wait for go");
            end if;
            for i in kp1..n loop
              if i mod nt = id - 1 then
                put_line("task " & to_string(id) & " updates " & to_string(i));
                a(i,j) := a(i,j) + temp*a(i,k);
              end if;
            end loop;
            busy(id) := false;
          end loop;
        end loop;
      end if;
      if id = 1 then
        ipvt(n) := n;
        if AbsVal(a(n,n)) = zero
         then info := n;
        end if;
      end if;
    end do_job;
    procedure do_jobs is new Reporting_Workers(do_job);

  begin
    info := 0;
    if n > t
     then nt := t; do_jobs(t);
     else nt := n; do_jobs(n);
    end if;
  end reporting_lufac;

  procedure silent_lusolve
              ( t : in integer32;
                a : in QuadDobl_Complex_Matrices.Matrix;
                n : in integer32;
                ipvt : in Standard_Integer_Vectors.Vector;
                b : in out QuadDobl_Complex_Vectors.Vector ) is

    nm1 : constant integer32 := n-1;
    go_forward : integer32 := 0;
    go_back : integer32 := n+1;
    busy : boolean_array(1..t) := (1..t => false);
    done : boolean_array(1..t) := (1..t => false);
    nt : integer32;
    use QuadDobl_Complex_Numbers;
    temp : Complex_Number;

    procedure do_job ( id,nt : in integer32 ) is

      ll,kb : integer32;
 
    begin
      if nm1 >= 1 then                              -- solve L*y = b
        for k in 1..nm1 loop
          if id = 1 then
            loop
              exit when all_false(nt,busy); 
            end loop;
            ll := ipvt(k);
            temp := b(ll);
            if ll /= k then
              b(ll) := b(k);
              b(k) := temp;
            end if;
            busy := (1..nt => true);
            go_forward := k;
          else
            loop
              exit when (go_forward = k);
            end loop;
          end if;
          for i in (k+1)..n loop
            if i mod nt = id - 1 then
              b(i) := b(i) + temp*a(i,k);
            end if;
          end loop;
          busy(id) := false;
        end loop;
      end if;
      done(id) := true;
      loop
        exit when all_true(nt,done);
      end loop;
      for k in 1..n loop                            -- solve U*x = y
        kb := n+1-k;
        if id = 1 then
          loop
            exit when all_false(nt,busy);
          end loop;
          b(kb) := b(kb) / a(kb,kb);
          busy := (1..nt => true);
          go_back := kb;
        else
          loop
            exit when (go_back = kb);
          end loop;
        end if;
        for j in 1..(kb-1) loop
          if j mod nt = id - 1 then
            b(j) := b(j) - b(kb)*a(j,kb);
          end if;
        end loop;
        busy(id) := false;
      end loop;
    end do_job;
    procedure do_jobs is new Silent_Workers(do_job);

  begin
    if t > n
     then nt := n; do_jobs(n);
     else nt := t; do_jobs(t);
    end if;
  end silent_lusolve;

  procedure reporting_lusolve
              ( t : in integer32;
                a : in QuadDobl_Complex_Matrices.Matrix;
                n : in integer32;
                ipvt : in Standard_Integer_Vectors.Vector;
                b : in out QuadDobl_Complex_Vectors.Vector ) is

    nm1 : constant integer32 := n-1;
    go_forward : integer32 := 0;
    go_back : integer32 := n+1;
    busy : boolean_array(1..t) := (1..t => false);
    done : boolean_array(1..t) := (1..t => false);
    nt : integer32;
    use QuadDobl_Complex_Numbers;
    temp : Complex_Number;

    procedure do_job ( id,nt : in integer32 ) is

      ll,kb : integer32;
 
    begin
      if nm1 >= 1 then                                  -- solve L*y = b
        for k in 1..nm1 loop
          if id = 1 then
            show(nt,busy,"busy");
            put_line("task 1 waits till none are busy ...");
            loop
              exit when all_false(nt,busy); 
            end loop;
            put_line("task 1 pivots at k = " & to_string(k));
            ll := ipvt(k);
            temp := b(ll);
            if ll /= k then
              b(ll) := b(k);
              b(k) := temp;
            end if;
            busy := (1..nt => true);
            go_forward := k;
            put_line("**** task 1 sets go to " & to_string(go_forward));
          else
            put_line("task " & to_string(id) & " waits to update"
                  & " at go = " & to_string(go_forward));
            loop
              exit when (go_forward = k);
            end loop;
            put_line("task " & to_string(id) & " finished wait for go");
          end if;
          for i in (k+1)..n loop
            if i mod nt = id - 1 then
              put_line("task " & to_string(id) & " updates " & to_string(i));
              b(i) := b(i) + temp*a(i,k);
            end if;
          end loop;
          busy(id) := false;
          put_line("task " & to_string(id) & " done at k = " & to_string(k));
        end loop;
      end if;
      done(id) := true;
      put_line("task " & to_string(id) & " is done with L*y = b");
      loop
        exit when all_true(nt,done);
      end loop;
      for k in 1..n loop                            -- solve U*x = y
        kb := n+1-k;
        if id = 1 then
          show(nt,busy,"busy");
          put_line("task 1 waits till none are busy ...");
          loop
            exit when all_false(nt,busy);
          end loop;
          put_line("task 1 sets temp at kb = " & to_string(kb));
          b(kb) := b(kb) / a(kb,kb);
          put_line("**** task 1 sets go to " & to_string(go_back));
          busy := (1..nt => true);
          go_back := kb;
        else
          put_line("task " & to_string(id) & " waits to update"
                & " at go = " & to_string(go_back));
          loop
            exit when (go_back = kb);
          end loop;
          put_line("task " & to_string(id) & " finished wait for go");
        end if;
        for j in 1..(kb-1) loop
          if j mod nt = id - 1 then
            put_line("task " & to_string(id) & " updates " & to_string(j));
            b(j) := b(j) - b(kb)*a(j,kb);
          end if;
        end loop;
        busy(id) := false;
        put_line("task " & to_string(id) & " done at kb = " & to_string(kb));
      end loop;
    end do_job;
    procedure do_jobs is new Reporting_Workers(do_job);

  begin
    if t > n
     then nt := n; do_jobs(n);
     else nt := t; do_jobs(t);
    end if;
  end reporting_lusolve;

end multitasking_linear_solvers;
