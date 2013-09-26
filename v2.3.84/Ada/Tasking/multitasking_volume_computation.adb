with text_io;                            use text_io;
with Timing_Package;                     use Timing_Package;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Multitasking,Semaphore;
with Standard_Floating_Vectors;
with Standard_Integer_Linear_Solvers;
with Lists_of_Floating_Vectors;          use Lists_of_Floating_Vectors;

package body Multitasking_Volume_Computation is

  procedure Extract_Support_Vectors
              ( mic : in Mixed_Cell;
                A : out Standard_Integer_Matrices.Matrix ) is

    tmp : List;
    lead,next : Standard_Floating_Vectors.Link_to_Vector;
    row : integer32 := A'first(1) - 1;

  begin
    for i in mic.pts'range loop
      lead := Head_Of(mic.pts(i));
      tmp := Tail_Of(mic.pts(i)); 
      while not Is_Null(tmp) loop
        next := Head_Of(tmp);
        row := row + 1;
        for j in A'range(2) loop
          A(row,j) := integer32(next(j)) - integer32(lead(j));
        end loop;
        tmp := Tail_Of(tmp);
      end loop;
    end loop;
  end Extract_Support_Vectors;

  procedure Reporting_Static_Multithreaded_Mixed_Volume
              ( nt,dim : in integer32;
                mcc : in out Mixed_Subdivision;
                vol : out Standard_Integer_Vectors.Vector;
                mixvol : out natural32 ) is

    timer : Timing_Widget;

    procedure Job ( i,n : integer32 ) is

      tmp : Mixed_Subdivision := mcc;
      mic : Mixed_Cell;
      mv : integer32;
      A : Standard_Integer_Matrices.Matrix(1..dim,1..dim);

    begin
      put_line("hello from task " & Multitasking.to_string(natural32(i)));
      for k in 1..Length_Of(mcc) loop
        if integer32(k) mod n = i-1 then
          put_line("task " & Multitasking.to_string(natural32(i))
                           & " computes cell "
                           & Multitasking.to_string(k));
          mic := Head_Of(tmp);
          Extract_Support_Vectors(mic,A);
          mv := Standard_Integer_Linear_Solvers.Det(A);
          if mv > 0
           then vol(i) := vol(i) + mv;
           else vol(i) := vol(i) - mv;
          end if;
        end if;
        tmp := Tail_Of(tmp);
      end loop;
    end Job;
    procedure do_jobs is new Multitasking.Reporting_Workers(Job);

  begin
    vol := (vol'range => 0);
    tstart(timer);
    do_jobs(nt);
    mixvol := natural32(Standard_Integer_Vectors.sum(vol));
    tstop(timer);
    new_line;
    put("the mixed volume : "); put(mixvol,1); new_line;
    new_line;
    print_times(standard_output,timer,"multithreaded mixed volume");
  end Reporting_Static_Multithreaded_Mixed_Volume;

  procedure Reporting_Dynamic_Multithreaded_Mixed_Volume
              ( nt,dim : in integer32;
                mcc : in out Mixed_Subdivision;
                vol : out Standard_Integer_Vectors.Vector;
                mixvol : out natural32 ) is

    timer : Timing_Widget;
    ptr : Mixed_Subdivision;
    cnt : natural32 := 0;
    s : Semaphore.Lock;

    procedure Next_Cell ( i,n : integer32 ) is

      myjob : natural32;
      myptr : Mixed_Subdivision;
      mic : Mixed_Cell;
      mv : integer32;
      A : Standard_Integer_Matrices.Matrix(1..dim,1..dim);

    begin
      put_line("hello from task " & Multitasking.to_string(natural32(i)));
      loop
        Semaphore.Request(s);
        if cnt = 0 then
          cnt := 1;
          ptr := mcc;
        else
          cnt := cnt + 1;
          if not Is_Null(ptr)
           then ptr := Tail_Of(ptr);
          end if;
        end if;
        myptr := ptr;
        myjob := cnt;
        Semaphore.Release(s);
        exit when Is_Null(myptr);
        put_line("task " & Multitasking.to_string(natural32(i))
                         & " computes cell "
                         & Multitasking.to_string(myjob));
        mic := Head_Of(myptr);
        Extract_Support_Vectors(mic,A);
        mv := Standard_Integer_Linear_Solvers.Det(A);
        if mv > 0
         then vol(i) := vol(i) + mv;
         else vol(i) := vol(i) - mv;
        end if;
      end loop;
    end Next_Cell;
    procedure do_jobs is new Multitasking.Silent_Workers(Next_Cell);

  begin
    vol := (vol'range => 0);
    tstart(timer);
    do_jobs(nt);
    mixvol := natural32(Standard_Integer_Vectors.sum(vol));
    tstop(timer);
    new_line;
    put("the mixed volume : "); put(mixvol,1); new_line;
    new_line;
    print_times(standard_output,timer,"multithreaded mixed volume");
  end Reporting_Dynamic_Multithreaded_Mixed_Volume;

  procedure Reporting_Static_Multithreaded_Mixed_Volume
              ( nt,dim : in integer32;
                mcc : in out Mixed_Subdivision;
                mixvol : out natural32 ) is

    vol : Standard_Integer_Vectors.Vector(1..nt);

  begin
    Reporting_Static_Multithreaded_Mixed_Volume(nt,dim,mcc,vol,mixvol);
  end Reporting_Static_Multithreaded_Mixed_Volume;

  procedure Reporting_Dynamic_Multithreaded_Mixed_Volume
              ( nt,dim : in integer32;
                mcc : in out Mixed_Subdivision;
                mixvol : out natural32 ) is

    vol : Standard_Integer_Vectors.Vector(1..nt);

  begin
    Reporting_Dynamic_Multithreaded_Mixed_Volume(nt,dim,mcc,vol,mixvol);
  end Reporting_Dynamic_Multithreaded_Mixed_Volume;

  procedure Silent_Static_Multithreaded_Mixed_Volume
              ( nt,dim : in integer32;
                mcc : in out Mixed_Subdivision;
                vol : out Standard_Integer_Vectors.Vector;
                mixvol : out natural32 ) is

    procedure Job ( i,n : integer32 ) is

      tmp : Mixed_Subdivision := mcc;
      mic : Mixed_Cell;
      mv : integer32;
      A : Standard_Integer_Matrices.Matrix(1..dim,1..dim);

    begin
      for k in 1..Length_Of(mcc) loop
        if integer32(k) mod n = i-1 then
          mic := Head_Of(tmp);
          Extract_Support_Vectors(mic,A);
          mv := Standard_Integer_Linear_Solvers.Det(A);
          if mv > 0
           then vol(i) := vol(i) + mv;
           else vol(i) := vol(i) - mv;
          end if;
        end if;
        tmp := Tail_Of(tmp);
      end loop;
    end Job;
    procedure do_jobs is new Multitasking.Silent_Workers(Job);

  begin
    vol := (vol'range => 0);
    do_jobs(nt);
    mixvol := natural32(Standard_Integer_Vectors.sum(vol));
  end Silent_Static_Multithreaded_Mixed_Volume;

  procedure Silent_Dynamic_Multithreaded_Mixed_Volume
              ( nt,dim : in integer32;
                mcc : in out Mixed_Subdivision;
                vol : out Standard_Integer_Vectors.Vector;
                mixvol : out natural32 ) is

    ptr : Mixed_Subdivision;
    first : boolean := true;
    s : Semaphore.Lock;

    procedure Next_Cell ( i,n : integer32 ) is

      myptr : Mixed_Subdivision;
      mic : Mixed_Cell;
      mv : integer32;
      A : Standard_Integer_Matrices.Matrix(1..dim,1..dim);

    begin
      loop
        Semaphore.Request(s);
        if first then
          first := false;
          ptr := mcc;
        else
          if not Is_Null(ptr)
           then ptr := Tail_Of(ptr);
          end if;
        end if;
        myptr := ptr;
        Semaphore.Release(s);
        exit when Is_Null(myptr);
        mic := Head_Of(myptr);
        Extract_Support_Vectors(mic,A);
        mv := Standard_Integer_Linear_Solvers.Det(A);
        if mv > 0
         then vol(i) := vol(i) + mv;
         else vol(i) := vol(i) - mv;
        end if;
      end loop;
    end Next_Cell;
    procedure do_jobs is new Multitasking.Silent_Workers(Next_Cell);

  begin
    vol := (vol'range => 0);
    do_jobs(nt);
    mixvol := natural32(Standard_Integer_Vectors.sum(vol));
  end Silent_Dynamic_Multithreaded_Mixed_Volume;

  procedure Silent_Static_Multithreaded_Mixed_Volume
              ( nt,dim : in integer32;
                mcc : in out Mixed_Subdivision;
                mixvol : out natural32 ) is

    vol : Standard_Integer_Vectors.Vector(1..nt);

  begin
    Silent_Static_Multithreaded_Mixed_Volume(nt,dim,mcc,vol,mixvol);
  end Silent_Static_Multithreaded_Mixed_Volume;

  procedure Silent_Dynamic_Multithreaded_Mixed_Volume
              ( nt,dim : in integer32;
                mcc : in out Mixed_Subdivision;
                mixvol : out natural32 ) is

    vol : Standard_Integer_Vectors.Vector(1..nt);

  begin
    Silent_Dynamic_Multithreaded_Mixed_Volume(nt,dim,mcc,vol,mixvol);
  end Silent_Dynamic_Multithreaded_Mixed_Volume;

end Multitasking_Volume_Computation;
