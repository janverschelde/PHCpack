with text_io,integer_io;                 use text_io,integer_io;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Double_Double_Numbers;              use Double_Double_Numbers;
with Quad_Double_Numbers;                use Quad_Double_Numbers;
with Standard_Complex_Numbers;
with DoblDobl_Complex_Numbers;
with QuadDobl_Complex_Numbers;
with Multitasking;                       use Multitasking;

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
              ( t : in positive;
                a : in out Standard_Complex_Matrices.Matrix;
                n : in positive;
                ipvt : out Standard_Natural_Vectors.Vector;
                info : out natural ) is

    nm1 : constant integer := n-1;
    go : natural := 0;
    go_elim : natural := 0;
    busy : boolean_array(1..t) := (1..t => false);
    nt : natural;
    use Standard_Complex_Numbers;
    temp : Complex_Number;

    procedure do_job ( id,nt : in natural ) is

      kp1,ll : integer;
      smax : double_float;

    begin
      if nm1 >= 1 then
        for k in 1..nm1 loop
          kp1 := k + 1;                            -- find the pivot index ll
          if id = 1 then
            loop
              exit when all_false(nt,busy);
            end loop;
            ll := k; smax := cabs(a(k,k));
            for i in kp1..n loop
              if cabs(a(i,k)) > smax then
                ll := i;
                smax := cabs(a(i,k));
              end if;
            end loop;
            ipvt(k) := ll;
            if smax = 0.0 then       -- this column is already triangularized
              info := k;
            else
              if ll /= k then                     -- interchange if necessary
                temp := a(ll,k);
                a(ll,k) := a(k,k);
                a(k,k) := temp;
              end if;
             -- temp := -Create(1.0)/a(k,k);           -- compute multipliers
              temp := -Create(1);
              temp := temp/a(k,k);
            end if;
            busy(1..nt) := (1..nt => true);
            go := k;
          else
            loop
              exit when (go = k);
            end loop;
          end if;
          for i in kp1..n loop
            if i mod nt = id - 1 then
              a(i,k) := temp*a(i,k);
            end if;
          end loop;
          busy(id) := false;
          for j in kp1..n loop       -- row elimination with column indexing
            if id = 1 then
              loop
                exit when all_false(nt,busy);
              end loop;
              temp := a(ll,j);
              if ll /= k then
                a(ll,j) := a(k,j);
                a(k,j) := temp;
              end if;
              busy(1..nt) := (1..nt => true);
              go_elim := k*(n-1) + j;
            else
              loop
                exit when (go_elim = (k*(n-1)+j));
              end loop;
            end if;
            for i in kp1..n loop
              if i mod nt = id - 1 then
                a(i,j) := a(i,j) + temp*a(i,k);
              end if;
            end loop;
            busy(id) := false;
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
     then nt := t; do_jobs(t);
     else nt := n; do_jobs(n);
    end if;
  end silent_lufac;

  procedure reporting_lufac
              ( t : in positive;
                a : in out Standard_Complex_Matrices.Matrix;
                n : in positive;
                ipvt : out Standard_Natural_Vectors.Vector;
                info : out natural ) is

    nm1 : constant integer := n-1;
    go : natural := 0;
    go_elim : natural := 0;
    busy : boolean_array(1..t) := (1..t => false);
    nt : natural;
    use Standard_Complex_Numbers;
    temp : Complex_Number;

    procedure do_job ( id,nt : in natural ) is

      kp1,ll : integer;
      smax : double_float;

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
            if smax = 0.0 then       -- this column is already triangularized
              info := k;
            else
              if ll /= k then                     -- interchange if necessary
                temp := a(ll,k);
                a(ll,k) := a(k,k);
                a(k,k) := temp;
              end if;
              temp := -Create(1.0)/a(k,k);             -- compute multipliers
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
        if AbsVal(a(n,n)) = 0.0
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
              ( t : in positive;
                a : in Standard_Complex_Matrices.Matrix;
                n : in positive;
                ipvt : in Standard_Natural_Vectors.Vector;
                b : in out Standard_Complex_Vectors.Vector ) is

    nm1 : constant natural := n-1;
    go_forward : natural := 0;
    go_back : natural := n+1;
    busy : boolean_array(1..t) := (1..t => false);
    done : boolean_array(1..t) := (1..t => false);
    nt : natural;
    use Standard_Complex_Numbers;
    temp : Complex_Number;

    procedure do_job ( id,nt : in natural ) is

      ll,kb : integer;
 
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
              ( t : in positive;
                a : in Standard_Complex_Matrices.Matrix;
                n : in positive;
                ipvt : in Standard_Natural_Vectors.Vector;
                b : in out Standard_Complex_Vectors.Vector ) is

    nm1 : constant natural := n-1;
    go_forward : natural := 0;
    go_back : natural := n+1;
    busy : boolean_array(1..t) := (1..t => false);
    done : boolean_array(1..t) := (1..t => false);
    nt : natural;
    use Standard_Complex_Numbers;
    temp : Complex_Number;

    procedure do_job ( id,nt : in natural ) is

      ll,kb : integer;
 
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
              ( t : in positive;
                a : in out DoblDobl_Complex_Matrices.Matrix;
                n : in positive;
                ipvt : out Standard_Natural_Vectors.Vector;
                info : out natural ) is

    nm1 : constant integer := n-1;
    go : natural := 0;
    go_elim : natural := 0;
    busy : boolean_array(1..t) := (1..t => false);
    nt : natural;
    zero : constant double_double := Create(0);
    use DoblDobl_Complex_Numbers;
    temp : Complex_Number;

    procedure do_job ( id,nt : in natural ) is

      kp1,ll : integer;
      smax : double_double;

    begin
      if nm1 >= 1 then
        for k in 1..nm1 loop
          kp1 := k + 1;                            -- find the pivot index ll
          if id = 1 then
            loop
              exit when all_false(nt,busy);
            end loop;
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
             -- temp := -Create(1.0)/a(k,k);          -- compute multipliers
              temp := -Create(1);
              temp := temp/a(k,k);
            end if;
            busy(1..nt) := (1..nt => true);
            go := k;
          else
            loop
              exit when (go = k);
            end loop;
          end if;
          for i in kp1..n loop
            if i mod nt = id - 1 then
              a(i,k) := temp*a(i,k);
            end if;
          end loop;
          busy(id) := false;
          for j in kp1..n loop       -- row elimination with column indexing
            if id = 1 then
              loop
                exit when all_false(nt,busy);
              end loop;
              temp := a(ll,j);
              if ll /= k then
                a(ll,j) := a(k,j);
                a(k,j) := temp;
              end if;
              busy(1..nt) := (1..nt => true);
              go_elim := k*(n-1) + j;
            else
              loop
                exit when (go_elim = (k*(n-1)+j));
              end loop;
            end if;
            for i in kp1..n loop
              if i mod nt = id - 1 then
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
    procedure do_jobs is new Silent_Workers(do_job);

  begin
    info := 0;
    if n > t
     then nt := t; do_jobs(t);
     else nt := n; do_jobs(n);
    end if;
  end silent_lufac;

  procedure reporting_lufac
              ( t : in positive;
                a : in out DoblDobl_Complex_Matrices.Matrix;
                n : in positive;
                ipvt : out Standard_Natural_Vectors.Vector;
                info : out natural ) is

    nm1 : constant integer := n-1;
    go : natural := 0;
    go_elim : natural := 0;
    busy : boolean_array(1..t) := (1..t => false);
    nt : natural;
    zero : constant double_double := Create(0);
    use DoblDobl_Complex_Numbers;
    temp : Complex_Number;

    procedure do_job ( id,nt : in natural ) is

      kp1,ll : integer;
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
              temp := -Create(1);
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
              ( t : in positive;
                a : in DoblDobl_Complex_Matrices.Matrix;
                n : in positive;
                ipvt : in Standard_Natural_Vectors.Vector;
                b : in out DoblDobl_Complex_Vectors.Vector ) is

    nm1 : constant natural := n-1;
    go_forward : natural := 0;
    go_back : natural := n+1;
    busy : boolean_array(1..t) := (1..t => false);
    done : boolean_array(1..t) := (1..t => false);
    nt : natural;
    use DoblDobl_Complex_Numbers;
    temp : Complex_Number;

    procedure do_job ( id,nt : in natural ) is

      ll,kb : integer;
 
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
              ( t : in positive;
                a : in DoblDobl_Complex_Matrices.Matrix;
                n : in positive;
                ipvt : in Standard_Natural_Vectors.Vector;
                b : in out DoblDobl_Complex_Vectors.Vector ) is

    nm1 : constant natural := n-1;
    go_forward : natural := 0;
    go_back : natural := n+1;
    busy : boolean_array(1..t) := (1..t => false);
    done : boolean_array(1..t) := (1..t => false);
    nt : natural;
    use DoblDobl_Complex_Numbers;
    temp : Complex_Number;

    procedure do_job ( id,nt : in natural ) is

      ll,kb : integer;
 
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
              ( t : in positive;
                a : in out QuadDobl_Complex_Matrices.Matrix;
                n : in positive;
                ipvt : out Standard_Natural_Vectors.Vector;
                info : out natural ) is

    nm1 : constant integer := n-1;
    go : natural := 0;
    go_elim : natural := 0;
    busy : boolean_array(1..t) := (1..t => false);
    nt : natural;
    zero : constant quad_double := Create(0);
    use QuadDobl_Complex_Numbers;
    temp : Complex_Number;

    procedure do_job ( id,nt : in natural ) is

      kp1,ll : integer;
      smax : quad_double;

    begin
      if nm1 >= 1 then
        for k in 1..nm1 loop
          kp1 := k + 1;                            -- find the pivot index ll
          if id = 1 then
            loop
              exit when all_false(nt,busy);
            end loop;
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
             -- temp := -Create(1.0)/a(k,k);          -- compute multipliers
              temp := -Create(1);
              temp := temp/a(k,k);
            end if;
            busy(1..nt) := (1..nt => true);
            go := k;
          else
            loop
              exit when (go = k);
            end loop;
          end if;
          for i in kp1..n loop
            if i mod nt = id - 1 then
              a(i,k) := temp*a(i,k);
            end if;
          end loop;
          busy(id) := false;
          for j in kp1..n loop       -- row elimination with column indexing
            if id = 1 then
              loop
                exit when all_false(nt,busy);
              end loop;
              temp := a(ll,j);
              if ll /= k then
                a(ll,j) := a(k,j);
                a(k,j) := temp;
              end if;
              busy(1..nt) := (1..nt => true);
              go_elim := k*(n-1) + j;
            else
              loop
                exit when (go_elim = (k*(n-1)+j));
              end loop;
            end if;
            for i in kp1..n loop
              if i mod nt = id - 1 then
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
    procedure do_jobs is new Silent_Workers(do_job);

  begin
    info := 0;
    if n > t
     then nt := t; do_jobs(t);
     else nt := n; do_jobs(n);
    end if;
  end silent_lufac;

  procedure reporting_lufac
              ( t : in positive;
                a : in out QuadDobl_Complex_Matrices.Matrix;
                n : in positive;
                ipvt : out Standard_Natural_Vectors.Vector;
                info : out natural ) is

    nm1 : constant integer := n-1;
    go : natural := 0;
    go_elim : natural := 0;
    busy : boolean_array(1..t) := (1..t => false);
    nt : natural;
    zero : constant quad_double := Create(0);
    use QuadDobl_Complex_Numbers;
    temp : Complex_Number;

    procedure do_job ( id,nt : in natural ) is

      kp1,ll : integer;
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
              temp := -Create(1);
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
              ( t : in positive;
                a : in QuadDobl_Complex_Matrices.Matrix;
                n : in positive;
                ipvt : in Standard_Natural_Vectors.Vector;
                b : in out QuadDobl_Complex_Vectors.Vector ) is

    nm1 : constant natural := n-1;
    go_forward : natural := 0;
    go_back : natural := n+1;
    busy : boolean_array(1..t) := (1..t => false);
    done : boolean_array(1..t) := (1..t => false);
    nt : natural;
    use QuadDobl_Complex_Numbers;
    temp : Complex_Number;

    procedure do_job ( id,nt : in natural ) is

      ll,kb : integer;
 
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
              ( t : in positive;
                a : in QuadDobl_Complex_Matrices.Matrix;
                n : in positive;
                ipvt : in Standard_Natural_Vectors.Vector;
                b : in out QuadDobl_Complex_Vectors.Vector ) is

    nm1 : constant natural := n-1;
    go_forward : natural := 0;
    go_back : natural := n+1;
    busy : boolean_array(1..t) := (1..t => false);
    done : boolean_array(1..t) := (1..t => false);
    nt : natural;
    use QuadDobl_Complex_Numbers;
    temp : Complex_Number;

    procedure do_job ( id,nt : in natural ) is

      ll,kb : integer;
 
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
