with text_io;                           use text_io;
with Communications_with_User;          use Communications_with_User;
with Timing_Package;                    use Timing_Package;
with Standard_Natural_Numbers_io;       use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with Standard_Complex_Numbers;          use Standard_Complex_Numbers;
with Standard_Random_Numbers;           use Standard_Random_Numbers;
with Standard_Complex_Poly_Systems;     use Standard_Complex_Poly_Systems;
with Standard_Complex_Poly_Systems_io;  use Standard_Complex_Poly_Systems_io;
with Standard_Homotopy;
with Standard_Complex_Solutions;        use Standard_Complex_Solutions;
with Standard_System_and_Solutions_io;  use Standard_System_and_Solutions_io;
with Main_Poly_Continuation;            use Main_Poly_Continuation;
with Multitasking_Continuation;         use Multitasking_Continuation;

procedure ts_mtcont is

-- DESCRIPTION :
--   Test on multitasking polynomial continuation.

  procedure Main is

    p,q : Link_to_Poly_Sys;
    s : Solution_List;
    n : integer32 := 0;
    file : file_type;
    timer : Timing_Widget;
    gamma : constant Complex_Number := Random1;
    ans : character;

  begin
    new_line;
    put_line("Test on multitasking path tracking ...");
    new_line;
    put_line("Reading a target system ..."); get(p);
    new_line;
    put_line("Reading a start system and its solutions ...");
    get(q,s);
    new_line;
    put("Read "); put(Length_Of(s),1); put_line(" solutions.");
    new_line;
    put("Give the number of threads : "); get(n);
    new_line;
    skip_line;
    put_line("Reading the name of the output file...");
    Read_Name_and_Create_File(file);
    Standard_Homotopy.Create(p.all,q.all,2,gamma);
    new_line;
    Driver_for_Continuation_Parameters(file);
    new_line;
    put("Do you want to monitor the progress on screen ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans = 'y' then
      tstart(timer);
      Reporting_Multitasking_Path_Tracker(s,n);
      tstop(timer);
    else
      tstart(timer);
      Silent_Multitasking_Path_Tracker(s,n);
      tstop(timer);
    end if;
    put(file,p.all,s);
    new_line(file);
    print_times(file,timer,"multitasking continuation");
  end Main;

begin
  Main;
end ts_mtcont;
