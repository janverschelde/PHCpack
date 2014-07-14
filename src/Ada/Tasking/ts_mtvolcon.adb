with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with Timing_Package;                     use Timing_Package;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;
with Double_Double_Numbers_io;           use Double_Double_Numbers_io;
with Quad_Double_Numbers_io;             use Quad_Double_Numbers_io;
with Standard_Integer_Vectors;
with Standard_Floating_Vectors;
with Double_Double_Vectors;
with Quad_Double_Vectors;
with Lists_of_Floating_Vectors;          use Lists_of_Floating_Vectors;
with Standard_Complex_Laur_Systems;
with Standard_Complex_Laur_Systems_io;   use Standard_Complex_Laur_Systems_io;
with DoblDobl_Complex_Laur_Systems;
with DoblDobl_Complex_Laur_Systems_io;   use DoblDobl_Complex_Laur_Systems_io;
with QuadDobl_Complex_Laur_Systems;
with QuadDobl_Complex_Laur_Systems_io;   use QuadDobl_Complex_Laur_Systems_io;
with Standard_Complex_Solutions;
with Standard_Complex_Solutions_io;      use Standard_Complex_Solutions_io;
with DoblDobl_Complex_Solutions;
with DoblDobl_Complex_Solutions_io;      use DoblDobl_Complex_Solutions_io;
with QuadDobl_Complex_Solutions;
with QuadDobl_Complex_Solutions_io;      use QuadDobl_Complex_Solutions_io;
with Floating_Mixed_Subdivisions;        use Floating_Mixed_Subdivisions;
with Floating_Mixed_Subdivisions_io;     use Floating_Mixed_Subdivisions_io;
with Continuation_Parameters;
with Polyhedral_Start_Systems;           use Polyhedral_Start_Systems;
with Multitasking_Polyhedral_Starters;   use Multitasking_Polyhedral_Starters;
with Multitasking_Polyhedral_Trackers;   use Multitasking_Polyhedral_Trackers;

procedure ts_mtvolcon is

-- DESCRIPTION :
--   Development of multithreaded polyhedral continuation.

  procedure Standard_Solve_Start_Systems
              ( nt : in integer32;
                q : in Standard_Complex_Laur_Systems.Laur_Sys;
                m : in integer32; mix : in Standard_Integer_Vectors.Vector;
                mcc : in out Mixed_Subdivision ) is

  -- DESCRIPTION :
  --   Solves all initial form systems defined by the mixed-cell
  --   configuration in mcc to start the polyhedral homotopies in q,
  --   in standard double precision using nt tasks.

    use Standard_Complex_Solutions;

    mv : natural32;
    sols : Array_of_Solution_Lists(1..nt);
    res : Standard_Floating_Vectors.Vector(1..nt);
    ans : character;
    timer : Timing_Widget;

  begin
    put("Intermediate output during computation wanted ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans = 'y' then    
      Reporting_Multithreaded_Solve_Start_Systems(nt,q,m,mix,mcc,mv,sols,res);
    else
      tstart(timer);
      Silent_Multithreaded_Solve_Start_Systems(nt,q,m,mix,mcc,mv,sols,res);
      tstop(timer);
      new_line;
      print_times(standard_output,timer,"multithreaded start solutions");
    end if;
    new_line;
    put("residuals : ");
    for i in res'range loop
      put(res(i),3);
    end loop;
    new_line;
  end Standard_Solve_Start_Systems;

  procedure Standard_Compute_Start_Solutions
              ( q : in Standard_Complex_Laur_Systems.Laur_Sys;
                m : in integer32;
                mix : in Standard_Integer_Vectors.Vector;
                mcc : in out Mixed_Subdivision ) is

  -- DESCRIPTION :
  --   First solves all initial form systems defined by the mixed-cell
  --   configuration in mcc by one single thread, then prompts the user
  --   for the number of tasks for multithreaded solving,
  --   in standard double precision.

    timer : Timing_Widget;
    nt : integer32 := 0;

  begin
    if m = q'last then
      tstart(timer);
      Fully_Mixed_Start_Systems(q,mcc);
      tstop(timer);
      new_line;
      print_times(standard_output,timer,"single thread start solutions");
    else
      tstart(timer);
      Semi_Mixed_Start_Systems(q,natural32(m),mix,mcc);
      tstop(timer);
      new_line;
      print_times(standard_output,timer,"single thread start solutions");
    end if;
    new_line;
    put("Give the number of threads : "); get(nt);
    Standard_Solve_Start_Systems(nt,q,m,mix,mcc);
  end Standard_Compute_Start_Solutions;

  procedure DoblDobl_Solve_Start_Systems
              ( nt : in integer32;
                q : in DoblDobl_Complex_Laur_Systems.Laur_Sys;
                m : in integer32; mix : in Standard_Integer_Vectors.Vector;
                mcc : in out Mixed_Subdivision ) is

  -- DESCRIPTION :
  --   Solves all initial form systems defined by the mixed-cell
  --   configuration in mcc to start the polyhedral homotopies in q,
  --   in double double precision using nt tasks.

    use DoblDobl_Complex_Solutions;

    mv : natural32;
    sols : Array_of_Solution_Lists(1..nt);
    res : Double_Double_Vectors.Vector(1..nt);
    ans : character;
    timer : Timing_Widget;

  begin
    put("Intermediate output during computation wanted ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans = 'y' then    
      Reporting_Multithreaded_Solve_Start_Systems(nt,q,m,mix,mcc,mv,sols,res);
    else
      tstart(timer);
      Silent_Multithreaded_Solve_Start_Systems(nt,q,m,mix,mcc,mv,sols,res);
      tstop(timer);
      new_line;
      print_times(standard_output,timer,"multithreaded start solutions");
    end if;
    new_line;
    put("residuals : ");
    for i in res'range loop
      put(" "); put(res(i),3);
    end loop;
    new_line;
  end DoblDobl_Solve_Start_Systems;

  procedure QuadDobl_Solve_Start_Systems
              ( nt : in integer32;
                q : in QuadDobl_Complex_Laur_Systems.Laur_Sys;
                m : in integer32; mix : in Standard_Integer_Vectors.Vector;
                mcc : in out Mixed_Subdivision ) is

  -- DESCRIPTION :
  --   Solves all initial form systems defined by the mixed-cell
  --   configuration in mcc to start the polyhedral homotopies in q,
  --   in double double precision using nt tasks.

    use QuadDobl_Complex_Solutions;

    mv : natural32;
    sols : Array_of_Solution_Lists(1..nt);
    res : Quad_Double_Vectors.Vector(1..nt);
    ans : character;
    timer : Timing_Widget;

  begin
    put("Intermediate output during computation wanted ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans = 'y' then    
      Reporting_Multithreaded_Solve_Start_Systems(nt,q,m,mix,mcc,mv,sols,res);
    else
      tstart(timer);
      Silent_Multithreaded_Solve_Start_Systems(nt,q,m,mix,mcc,mv,sols,res);
      tstop(timer);
      new_line;
      print_times(standard_output,timer,"multithreaded start solutions");
    end if;
    new_line;
    put("residuals : ");
    for i in res'range loop
      put(" "); put(res(i),3);
    end loop;
    new_line;
  end QuadDobl_Solve_Start_Systems;

  procedure DoblDobl_Compute_Start_Solutions
              ( q : in DoblDobl_Complex_Laur_Systems.Laur_Sys;
                m : in integer32;
                mix : in Standard_Integer_Vectors.Vector;
                mcc : in out Mixed_Subdivision ) is

  -- DESCRIPTION :
  --   First solves all initial form systems defined by the mixed-cell
  --   configuration in mcc by one single thread, then prompts the user
  --   for the number of tasks for multithreaded solving,
  --   in standard double precision.

    timer : Timing_Widget;
    nt : integer32 := 0;

  begin
    if m = q'last then
      tstart(timer);
      Fully_Mixed_Start_Systems(q,mcc);
      tstop(timer);
      new_line;
      print_times(standard_output,timer,"single thread start solutions");
    else
      tstart(timer);
      Semi_Mixed_Start_Systems(q,natural32(m),mix,mcc);
      tstop(timer);
      new_line;
      print_times(standard_output,timer,"single thread start solutions");
    end if;
    new_line;
    put("Give the number of threads : "); get(nt);
    DoblDobl_Solve_Start_Systems(nt,q,m,mix,mcc);
  end DoblDobl_Compute_Start_Solutions;

  procedure QuadDobl_Compute_Start_Solutions
              ( q : in QuadDobl_Complex_Laur_Systems.Laur_Sys;
                m : in integer32;
                mix : in Standard_Integer_Vectors.Vector;
                mcc : in out Mixed_Subdivision ) is

  -- DESCRIPTION :
  --   First solves all initial form systems defined by the mixed-cell
  --   configuration in mcc by one single thread, then prompts the user
  --   for the number of tasks for multithreaded solving,
  --   in standard double precision.

    timer : Timing_Widget;
    nt : integer32 := 0;

  begin
    if m = q'last then
      tstart(timer);
      Fully_Mixed_Start_Systems(q,mcc);
      tstop(timer);
      new_line;
      print_times(standard_output,timer,"single thread start solutions");
    else
      tstart(timer);
      Semi_Mixed_Start_Systems(q,natural32(m),mix,mcc);
      tstop(timer);
      new_line;
      print_times(standard_output,timer,"single thread start solutions");
    end if;
    new_line;
    put("Give the number of threads : "); get(nt);
    QuadDobl_Solve_Start_Systems(nt,q,m,mix,mcc);
  end QuadDobl_Compute_Start_Solutions;

  procedure Standard_Track_Paths
              ( q : in Standard_Complex_Laur_Systems.Laur_Sys;
                n,m : in integer32;
                mix : in Standard_Integer_Vectors.Vector;
                mcc : in out Mixed_Subdivision ) is

  -- DESCRIPTION :
  --   Prompts the suesr for an output file, the number of tasks,
  --   and whether intermediate output is needed before launching
  --   the multitasked polyhedral trackers in standard double precision.

    use Standard_Complex_Solutions;

    nt : integer32 := 0;
    ans : character;
    sols : Solution_List;
    file : file_type;
 
  begin
    put_line("Reading the name of the output file ...");
    Read_Name_and_Create_File(file);
    put_line(file,q);
    new_line;
    put("Give the number of threads : "); get(nt);
    put("Do you want output to screen ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans = 'n'
     then Silent_Multitasking_Path_Tracker(q,nt,n,m,mix,mcc,sols);
     else Reporting_Multitasking_Path_Tracker(file,q,nt,n,m,mix,mcc,sols);
    end if;
    new_line(file);
    put_line(file,"THE SOLUTIONS :");
    put(file,Length_Of(sols),natural32(Head_Of(sols).n),sols);
  end Standard_Track_Paths;

  procedure DoblDobl_Track_Paths
              ( q : in DoblDobl_Complex_Laur_Systems.Laur_Sys;
                n,m : in integer32;
                mix : in Standard_Integer_Vectors.Vector;
                mcc : in out Mixed_Subdivision ) is

  -- DESCRIPTION :
  --   Prompts the suesr for an output file, the number of tasks,
  --   and whether intermediate output is needed before launching
  --   the multitasked polyhedral trackers in double double precision.

    use DoblDobl_Complex_Solutions;

    nt : integer32 := 0;
    ans : character;
    sols : Solution_List;
    file : file_type;
 
  begin
    put_line("Reading the name of the output file ...");
    Read_Name_and_Create_File(file);
    put_line(file,q);
    new_line;
    put("Give the number of threads : "); get(nt);
    put("Do you want output to screen ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans = 'n'
     then Silent_Multitasking_Path_Tracker(q,nt,n,m,mix,mcc,sols);
     else Reporting_Multitasking_Path_Tracker(file,q,nt,n,m,mix,mcc,sols);
    end if;
    new_line(file);
    put_line(file,"THE SOLUTIONS :");
    put(file,Length_Of(sols),natural32(Head_Of(sols).n),sols);
  end DoblDobl_Track_Paths;

  procedure QuadDobl_Track_Paths
              ( q : in QuadDobl_Complex_Laur_Systems.Laur_Sys;
                n,m : in integer32;
                mix : in Standard_Integer_Vectors.Vector;
                mcc : in out Mixed_Subdivision ) is

  -- DESCRIPTION :
  --   Prompts the suesr for an output file, the number of tasks,
  --   and whether intermediate output is needed before launching
  --   the multitasked polyhedral trackers in double double precision.

    use QuadDobl_Complex_Solutions;

    nt : integer32 := 0;
    ans : character;
    sols : Solution_List;
    file : file_type;
 
  begin
    put_line("Reading the name of the output file ...");
    Read_Name_and_Create_File(file);
    put_line(file,q);
    new_line;
    put("Give the number of threads : "); get(nt);
    put("Do you want output to screen ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans = 'n'
     then Silent_Multitasking_Path_Tracker(q,nt,n,m,mix,mcc,sols);
     else Reporting_Multitasking_Path_Tracker(file,q,nt,n,m,mix,mcc,sols);
    end if;
    new_line(file);
    put_line(file,"THE SOLUTIONS :");
    put(file,Length_Of(sols),natural32(Head_Of(sols).n),sols);
  end QuadDobl_Track_Paths;

  procedure Standard_Test is

  -- DESCRIPTION :
  --   Prompts the user for a random coefficient polynomial system,
  --   a mixed-cell configuration, and then tests the multitasked
  --   polyhedral homotopies in standard double precision.

    file : file_type;
    q : Standard_Complex_Laur_Systems.Link_to_Laur_Sys;
    n,m : integer32;
    mix : Standard_Integer_Vectors.Link_to_Vector;
    mcc : Mixed_Subdivision;
    ans : character;

  begin
    new_line;
    put_line("Reading a random coefficient system ...");
    Read_Name_and_Open_File(file);
    get(file,q);
    close(file);
    new_line;
    put_line("Reading a file for a mixed cell configuration ...");
    Read_Name_and_Open_File(file);
    get(file,natural32(n),natural32(m),mix,mcc);
    close(file);
    new_line;
    put("Read "); put(Length_Of(mcc),1); put_line(" cells.");
    new_line;
    put_line("MENU to test polyhedral homotopies : ");
    put_line("  1. solve only the start systems;");
    put_line("  2. solve the start systems and track all paths.");
    put("Type 1 or 2 to make your choice : ");
    Ask_Alternative(ans,"12");
    new_line;
    if ans = '1' then
      Standard_Compute_Start_Solutions(q.all,m,mix.all,mcc);
    else
      Continuation_Parameters.Tune(0);
      Standard_Track_Paths(q.all,n,m,mix.all,mcc);
    end if;
  end Standard_Test;

  procedure DoblDobl_Test is

  -- DESCRIPTION :
  --   Prompts the user for a random coefficient polynomial system,
  --   a mixed-cell configuration, and then tests the multitasked
  --   polyhedral homotopies in double double precision.

    file : file_type;
    q : DoblDobl_Complex_Laur_Systems.Link_to_Laur_Sys;
    n,m : integer32;
    mix : Standard_Integer_Vectors.Link_to_Vector;
    mcc : Mixed_Subdivision;
    ans : character;

  begin
    new_line;
    put_line("Reading a random coefficient system ...");
    Read_Name_and_Open_File(file);
    get(file,q);
    close(file);
    new_line;
    put_line("Reading a file for a mixed cell configuration ...");
    Read_Name_and_Open_File(file);
    get(file,natural32(n),natural32(m),mix,mcc);
    close(file);
    new_line;
    put("Read "); put(Length_Of(mcc),1); put_line(" cells.");
    new_line;
    put_line("MENU to test polyhedral homotopies : ");
    put_line("  1. solve only the start systems;");
    put_line("  2. solve the start systems and track all paths.");
    put("Type 1 or 2 to make your choice : ");
    Ask_Alternative(ans,"12");
    new_line;
    if ans = '1' then
      DoblDobl_Compute_Start_Solutions(q.all,m,mix.all,mcc);
    else
      Continuation_Parameters.Tune(0,32);
      DoblDobl_Track_Paths(q.all,n,m,mix.all,mcc);
    end if;
  end DoblDobl_Test;

  procedure QuadDobl_Test is

  -- DESCRIPTION :
  --   Prompts the user for a random coefficient polynomial system,
  --   a mixed-cell configuration, and then tests the multitasked
  --   polyhedral homotopies in quad double precision.

    file : file_type;
    q : QuadDobl_Complex_Laur_Systems.Link_to_Laur_Sys;
    n,m : integer32;
    mix : Standard_Integer_Vectors.Link_to_Vector;
    mcc : Mixed_Subdivision;
    ans : character;

  begin
    new_line;
    put_line("Reading a random coefficient system ...");
    Read_Name_and_Open_File(file);
    get(file,q);
    close(file);
    new_line;
    put_line("Reading a file for a mixed cell configuration ...");
    Read_Name_and_Open_File(file);
    get(file,natural32(n),natural32(m),mix,mcc);
    close(file);
    new_line;
    put("Read "); put(Length_Of(mcc),1); put_line(" cells.");
    new_line;
    put_line("MENU to test polyhedral homotopies : ");
    put_line("  1. solve only the start systems;");
    put_line("  2. solve the start systems and track all paths.");
    put("Type 1 or 2 to make your choice : ");
    Ask_Alternative(ans,"12");
    new_line;
    if ans = '1' then
      QuadDobl_Compute_Start_Solutions(q.all,m,mix.all,mcc);
    else
      Continuation_Parameters.Tune(0,64);
      QuadDobl_Track_Paths(q.all,n,m,mix.all,mcc);
    end if;
  end QuadDobl_Test;

  procedure Main is

    ans : character;

  begin
    new_line;
    put_line("MENU to test multitasked polyhedral homotopies :");
    put_line("  1. in standard double precision;");
    put_line("  2. in double double precision;");
    put_line("  3. in quad double precision.");
    put("Type 1, 2, or 3 to select the precision : ");
    Ask_Alternative(ans,"123");
    case ans is
      when '1' => Standard_Test;
      when '2' => DoblDobl_Test;
      when '3' => QuadDobl_Test;
      when others => null;
    end case;
  end Main;

begin
  Main;
end ts_mtvolcon;
