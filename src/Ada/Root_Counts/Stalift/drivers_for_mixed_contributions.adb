with Timing_Package;                     use Timing_Package;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Lists_of_Integer_Vectors;           use Lists_of_Integer_Vectors;
with Contributions_to_Mixed_Volume;      use Contributions_to_Mixed_Volume;

package body Drivers_for_Mixed_Contributions is

  procedure Count_Zero_Contributions
              ( file : in file_type; L,nl : in Array_of_Lists;
                nred : out natural32 ) is

  -- DESCRIPTION :
  --   Counts the number of zero contributions per component and reports
  --   the cardinalities on file.

    res,inc : natural32 := 0;

  begin
    new_line(file);
    put_line(file,"#Eliminated points per component : ");
    for i in l'range loop
      inc := Length_Of(l(i)) - Length_Of(nl(i));
      put(file,"  "); put(file,inc,1);
      res := res + inc;
    end loop;
    new_line(file);
    nred := res;
  end Count_Zero_Contributions;

-- TARGET ROUTINES :

  procedure Once_Simple_Sweep
              ( file : in file_type; L : in out Array_of_Lists;
                nred : out natural32 ) is

    nl : Array_of_Lists(L'range);
    timer : Timing_Widget;

  begin
    tstart(timer);
    nl := Simple_Sweep(L);
    tstop(timer);
    Count_Zero_Contributions(file,L,nl,nred);
    new_line(file);
    print_times(file,timer,"one simple sweep");
    new_line(file);
    Copy(nl,L); Deep_Clear(nl);
  end Once_Simple_Sweep;

  procedure Once_Exhaustive_Sweep
              ( file : in file_type; L : in out Array_of_Lists;
                nred : out natural32 ) is

    nl : Array_of_Lists(l'range);
    timer : Timing_Widget;

  begin
    tstart(timer);
    nl := Exhaustive_Sweep(L);
    tstop(timer);
    Count_Zero_Contributions(file,L,nl,nred);
    new_line(file);
    print_times(file,timer,"exhaustive sweep");
    new_line(file);
    Copy(nl,L); Deep_Clear(nl);
  end Once_Exhaustive_Sweep;

  procedure Full_Simple_Sweep
              ( file : in file_type; L : in out Array_of_Lists;
                nred : out natural32 ) is

    totnred,wrknred : natural32;
    timer : Timing_Widget;

  begin
    totnred := 0;
    tstart(timer);
    loop
      Once_Simple_Sweep(file,L,wrknred);
      exit when (wrknred = 0);
      totnred := totnred + wrknred;
    end loop;
    tstop(timer);
    new_line(file);
    print_times(file,timer,"full simple sweep");
    new_line(file);
    nred := totnred;
  end Full_Simple_Sweep;

  procedure Full_Exhaustive_Sweep
              ( file : in file_type; L : in out Array_of_Lists;
                nred : out natural32 ) is

    totnred,wrknred : natural32;
    timer : Timing_Widget;

  begin
    totnred := 0;
    tstart(timer);
    loop
      Once_Exhaustive_Sweep(file,L,wrknred);
      exit when (wrknred = 0);
      totnred := totnred + wrknred;
    end loop;
    tstop(timer);
    new_line(file);
    print_times(file,timer,"full exhaustive sweep");
    new_line(file);
    nred := totnred;
  end Full_Exhaustive_Sweep;

end Drivers_for_Mixed_Contributions;
