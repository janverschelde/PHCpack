with text_io;                           use text_io;
with Communications_with_User;          use Communications_with_User;
with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with Standard_Complex_Poly_Systems;
with Standard_Complex_Poly_Systems_io;  use Standard_Complex_Poly_Systems_io;
with Drivers_to_Track_Standard_Paths;
with Drivers_to_Track_DoblDobl_Paths;
with Drivers_to_Track_QuadDobl_Paths;
with Standard_Continuation_Data;
with Standard_Continuation_Data_io;
with DoblDobl_Continuation_Data;
with DoblDobl_Continuation_Data_io;
with QuadDobl_Continuation_Data;
with QuadDobl_Continuation_Data_io;

procedure ts_track is

-- DESCRIPTION :
--   Calls the main interactive driver to track solution paths,
--   for a cheater homotopy where all solutions of the start system
--   are read one after the other from file.

  procedure Call_Standard_Trackers
              ( target,start,output : in file_type ) is

    use Standard_Complex_Poly_Systems;
    use Drivers_to_Track_Standard_Paths;
    use Standard_Continuation_Data;

    procedure Write ( s : in Solu_Info ) is
    begin
      Standard_Continuation_Data_io.Write_Solution(standard_output,s);
    end Write;

    p : Link_to_Poly_Sys;
    ans : character;

  begin
    Read_Target_System(target,p);
    put(output,p'last,1); new_line(output);
    put(output,p.all);
    new_line;
    put("Do you want to see intermediate solutions ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans = 'y'
     then Read_Systems_and_Track(p.all,start,output,'3',Write'access);
     else Read_Systems_and_Track(p.all,start,output,'3');
    end if;
  end Call_Standard_Trackers;

  procedure Call_DoblDobl_Trackers
              ( target,start,output : in file_type ) is

    use Drivers_to_Track_DoblDobl_Paths;
    use DoblDobl_Continuation_Data;

    procedure Write ( s : in Solu_Info ) is
    begin
      DoblDobl_Continuation_Data_io.Write_Solution(standard_output,s);
    end Write;

    ans : character;

  begin
    new_line;
    put("Do you want to see intermediate solutions ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans = 'y'
     then Read_Systems_and_Track(target,start,output,Write'access);
     else Read_Systems_and_Track(target,start,output);
    end if;
  end Call_DoblDobl_Trackers;

  procedure Call_QuadDobl_Trackers
              ( target,start,output : in file_type ) is

    use Drivers_to_Track_QuadDobl_Paths;
    use QuadDobl_Continuation_Data;

    procedure Write ( s : in Solu_Info ) is
    begin
      QuadDobl_Continuation_Data_io.Write_Solution(standard_output,s);
    end Write;

    ans : character;

  begin
    new_line;
    put("Do you want to see intermediate solutions ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans = 'y'
     then Read_Systems_and_Track(target,start,output,Write'access);
     else Read_Systems_and_Track(target,start,output);
    end if;
  end Call_QuadDobl_Trackers;

  procedure Main is

    target_file,start_file,output_file : file_type;
    ans : character;

  begin
    new_line;
    put_line("Path tracking with incremental read/write of solutions.");
    new_line;
    put_line("Reading the name of the file for the target system.");
    Read_Name_and_Open_File(target_file);
    new_line;
    put_line("Reading the name of the file for the start system.");
    Read_Name_and_Open_File(start_file);
    new_line;
    put_line("Reading the name of the output file.");
    Read_Name_and_Create_File(output_file);
    new_line;
    put_line("MENU for arithmetic during path tracking : ");
    put_line("  0. standard double floating point arithmetic;");
    put_line("  1. double double floating point arithmetic; or");
    put_line("  2. quad double floating point arithmetic.");
    put("Type 0, 1, or 2 to select the arithmetic : ");
    Ask_Alternative(ans,"012");
    case ans is
      when '0' => Call_Standard_Trackers(target_file,start_file,output_file);
      when '1' => Call_DoblDobl_Trackers(target_file,start_file,output_file);
      when '2' => Call_QuadDobl_Trackers(target_file,start_file,output_file);
      when others => null;
    end case;
  end Main;

begin
  Main;
end ts_track;
