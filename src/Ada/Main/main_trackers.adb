with String_Splitters;                  use String_Splitters;
with Communications_with_User;          use Communications_with_User;
with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with Standard_Complex_Poly_Systems;
with Standard_Complex_Poly_Systems_io;  use Standard_Complex_Poly_Systems_io;
with DoblDobl_Complex_Poly_Systems;
with DoblDobl_Complex_Poly_Systems_io;  use DoblDobl_Complex_Poly_Systems_io;
with QuadDobl_Complex_Poly_Systems;
with QuadDobl_Complex_Poly_Systems_io;  use QuadDobl_Complex_Poly_Systems_io;
with Drivers_to_Track_Standard_Paths;
with Drivers_to_Track_DoblDobl_Paths;
with Drivers_to_Track_QuadDobl_Paths;
with Jumpstart_Polyhedral_Homotopies;   use Jumpstart_Polyhedral_Homotopies;
with Jumpstart_Diagonal_Homotopies;     use Jumpstart_Diagonal_Homotopies;

package body Main_Trackers is

  function Ask_for_Start_Type return character is

    res : character;

  begin
    new_line;
    put_line("MENU for type of start system or homotopy : ");
    put_line("  1. start system is based on total degree;");
    put_line("  2. a linear-product start system will be given;");
    put_line("  3. start system and start solutions are provided;");
    put_line("  4. polyhedral continuation on a generic system;");
    put_line("  5. diagonal homotopy to intersect algebraic sets;");
    put_line("  6. descend one level down in a cascade of homotopies;");
    put_line("  7. remove last slack variable in a witness set.");
    put("Type 1, 2, 3, 4, 5, 6, or 7 to select type : ");
    Ask_Alternative(res,"1234567");
    return res;
  end Ask_for_Start_Type;

  procedure Standard_Track
              ( target,start,output : in file_type;
                kind : in character ) is

    use Drivers_to_Track_Standard_Paths;
    tgtsys : Standard_Complex_Poly_Systems.Link_to_Poly_Sys;

  begin 
    Read_Target_System(target,tgtsys);
    put(output,tgtsys'last,1);
    new_line(output);
    put(output,tgtsys.all);
    Read_Systems_and_Track(tgtsys.all,start,output,kind);
  end Standard_Track;

  procedure DoblDobl_Track
              ( target,start,output : in file_type; kind : in character ) is

    use Drivers_to_Track_DoblDobl_Paths;
    tgtsys : DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys;

  begin 
    Read_Target_System(target,tgtsys);
    put(output,tgtsys'last,1);
    new_line(output);
    put(output,tgtsys.all);
    Read_Systems_and_Track(tgtsys.all,start,output,kind);
  end DoblDobl_Track;

  procedure QuadDobl_Track
              ( target,start,output : in file_type; kind : in character  ) is

    use Drivers_to_Track_QuadDobl_Paths;
    tgtsys : QuadDobl_Complex_Poly_Systems.Link_to_Poly_Sys;

  begin 
    Read_Target_System(target,tgtsys);
    put(output,tgtsys'last,1);
    new_line(output);
    put(output,tgtsys.all);
    Read_Systems_and_Track(tgtsys.all,start,output,kind);
  end QuadDobl_Track;

  procedure Track ( target,start,output : in file_type; kd : in character ) is

    ans : character;

  begin
    new_line;
    put_line("MENU for arithmetic during path tracking : ");
    put_line("  0. standard double floating point arithmetic;");
    put_line("  1. double double floating point arithmetic; or");
    put_line("  2. quad double floating point arithmetic.");
    put("Type 0, 1, or 2 to select the arithmetic : ");
    Ask_Alternative(ans,"012");
    case ans is
      when '0' => Standard_Track(target,start,output,kd);
      when '1' => DoblDobl_Track(target,start,output,kd);
      when '2' => QuadDobl_Track(target,start,output,kd);
      when others => null;
    end case;
  end Track;

  procedure Main ( targetfilename,startfilename,outfilename : in string;
                   verbose : in integer32 := 0 ) is

    target_file,start_file,output_file : file_type;
    target_name,start_name,output_name : Link_to_String;
    target : Standard_Complex_Poly_Systems.Link_to_Poly_Sys;
    start_type : character;

    use Drivers_to_Track_Standard_Paths;

  begin
    if verbose > 0 then
      put("At verbose level "); put(verbose,1);
      put_line(", in main_trackers.Main ...");
    end if;
    start_type := Ask_for_Start_Type;
    if targetfilename = "" then
      new_line;
      put_line("Reading the name of the file for the target system.");
      Read_Name_and_Open_File(target_file);
    else
      Open_Input_File(target_file,targetfilename,target_name);
    end if;
   -- postpone reading of target system in case of cheater homotopy
    if start_type /= '1' and start_type /= '2' and start_type /= '3'
     then Read_Target_System(target_file,target);
    end if;
    if outfilename = "" then
      new_line;
      put_line("Reading the name of the output file.");
      Read_Name_and_Create_File(output_file);
    else
      Create_Output_File(output_file,outfilename,output_name);
    end if;
    if (start_type = '5') then
      Jumpstart_Diagonal_Homotopy(target_file,output_file,target.all);
    elsif (start_type = '6') then
      Jumpstart_Cascade_Homotopy(target_file,output_file,target.all);
    elsif (start_type = '7') then
      Remove_Last_Slack_Variable(target_file,output_file,target.all);
    else
      if start_type /= '1' and start_type /= '2' and start_type /= '3' then
        put(output_file,target'last,1);
        new_line(output_file);
        put(output_file,target.all);
      end if;
      if start_type = '4' then
        if startfilename = "" then
          new_line;
          put_line("Reading file name for regular mixed-cell configuration.");
          Read_Name_and_Open_File(start_file);
        else 
          Open_Input_File(start_file,startfilename,start_name);
        end if;
        Read_Cells_and_Track(target.all,start_file,output_file);
      else
        if startfilename = "" then
          new_line;
          put_line("Reading the name of the file for the start system.");
          Read_Name_and_Open_File(start_file);
        else
          Open_Input_File(start_file,startfilename,start_name);
        end if;
        Track(target_file,start_file,output_file,start_type);
      end if;
    end if;
    close(target_file);
  end Main;

end Main_Trackers;
