with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with Timing_Package;                     use Timing_Package;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;
with Standard_Complex_Numbers;           use Standard_Complex_Numbers;
with Standard_Complex_Numbers_io;        use Standard_Complex_Numbers_io;
with Standard_Integer_Vectors;
with Standard_Integer_Vectors_io;        use Standard_Integer_Vectors_io;
with Standard_Integer_VecVecs;
with Standard_Floating_Vectors;
with Standard_Complex_Vectors;
with Standard_Complex_VecVecs;
with Standard_Integer_Matrices;
with Standard_Integer_Linear_Solvers;
with Lists_of_Floating_Vectors;          use Lists_of_Floating_Vectors;
with Arrays_of_Floating_Vector_Lists;
with Standard_Complex_Laur_Systems;      use Standard_Complex_Laur_Systems;
with Standard_Complex_Laur_Systems_io;   use Standard_Complex_Laur_Systems_io;
with Standard_Tableau_Formats;
with Standard_Complex_Solutions_io;      use Standard_Complex_Solutions_io;
with Floating_Mixed_Subdivisions;        use Floating_Mixed_Subdivisions;
with Floating_Mixed_Subdivisions_io;     use Floating_Mixed_Subdivisions_io;
with Standard_Complex_Solutions;         use Standard_Complex_Solutions;
with Standard_Radial_Solvers;
with Standard_Binomial_Systems;
with Standard_Binomial_Solvers;
with Continuation_Parameters;
with Multitasking;
with Multitasking_Volume_Computation;    use Multitasking_Volume_Computation;
with Polyhedral_Start_Systems;           use Polyhedral_Start_Systems;
with Multitasking_Polyhedral_Trackers;   use Multitasking_Polyhedral_Trackers;

procedure ts_mtvolcon is

-- DESCRIPTION :
--   Development of multithreaded polyhedral continuation.

  procedure Solve_Start_Systems
              ( nt : in integer32; q : in Laur_Sys;
                m : in integer32; mix : in Standard_Integer_Vectors.Vector;
                mcc : in out Mixed_Subdivision ) is

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
  end Solve_Start_Systems;

  procedure Compute_Start_Solutions
              ( q : in Laur_Sys; m : in integer32;
                mix : in Standard_Integer_Vectors.Vector;
                mcc : in out Mixed_Subdivision ) is

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
    Solve_Start_Systems(nt,q,m,mix,mcc);
  end Compute_Start_Solutions;

  procedure Track_Paths
              ( q : in Laur_Sys; n,m : in integer32;
                mix : in Standard_Integer_Vectors.Vector;
                mcc : in out Mixed_Subdivision ) is

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
  end Track_Paths;

  procedure Main is

    file : file_type;
    q : Link_to_Laur_Sys;
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
      Compute_Start_Solutions(q.all,m,mix.all,mcc);
    else
      Continuation_Parameters.Tune(2);
      Track_Paths(q.all,n,m,mix.all,mcc);
    end if;
  end Main;

begin
  Main;
end ts_mtvolcon;
