with System;
with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with Timing_Package;                     use Timing_Package;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Integer_Vectors;
with Lists_of_Floating_Vectors;          use Lists_of_Floating_Vectors;
with Floating_Mixed_Subdivisions;        use Floating_Mixed_Subdivisions;
with Floating_Mixed_Subdivisions_io;     use Floating_Mixed_Subdivisions_io;
with Mixed_Volume_Computation;
with Multitasking_Volume_Computation;    use Multitasking_Volume_Computation;

procedure ts_mtvol is

-- DESCRIPTION :
--   Development of multithreaded computation of the mixed volume.

  procedure Mixed_Volume
              ( n : in integer32; mix : in Standard_Integer_Vectors.Vector;
                mcc : in out Mixed_Subdivision ) is

  -- DESCRIPTION :
  --   Given a mixed cell configuration in mcc for a tuple of point sets
  --   with mixture in mix in n-space, the mixed volume is computed
  --   with one single thread.

    vol : natural32;
    timer : Timing_Widget;

  begin
    tstart(timer);
    Mixed_Volume_Computation.Mixed_Volume(n,mix,mcc,vol);
    tstop(timer);
    new_line;
    put("the mixed volume : "); put(vol,1); new_line;
    new_line;
    print_times(standard_output,timer,"single thread mixed volume");
  end Mixed_Volume;

  procedure Main is

    file : file_type;
    nt,n : integer32 := 0;
    m,vol : natural32;
    mix : Standard_Integer_Vectors.Link_to_Vector;
    mcc : Mixed_Subdivision;
    ans : character;
    static : boolean;
    timer : Timing_Widget;

  begin
    new_line;
    put_line("Reading a file for a mixed cell configuration ...");
    Read_Name_and_Open_File(file);
    get(file,natural32(n),m,mix,mcc);
    close(file);
    new_line;
    put("Read "); put(Length_Of(mcc),1); put_line(" cells.");
    Mixed_Volume(n,mix.all,mcc);
    new_line;
    put("Give the number of threads : "); get(nt);
    put("-> static or dynamic load balancing ? (s/d) ");
    Ask_Alternative(ans,"sd");
    static := (ans = 's');
    put("Output during computations ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans = 'y' then
      if static
       then Reporting_Static_Multithreaded_Mixed_Volume(nt,n,mcc,vol);
       else Reporting_Dynamic_Multithreaded_Mixed_Volume(nt,n,mcc,vol);
      end if;
    else
      tstart(timer);
      if static 
       then Silent_Static_Multithreaded_Mixed_Volume(nt,n,mcc,vol);
       else Silent_Dynamic_Multithreaded_Mixed_Volume(nt,n,mcc,vol);
      end if;
      tstop(timer);
      new_line;
      put("The mixed volume : "); put(vol,1); new_line;
      new_line;
      print_times(standard_output,timer,"multithreaded mixed volume");
    end if;
  end Main;

begin
  Main;
end ts_mtvol;
