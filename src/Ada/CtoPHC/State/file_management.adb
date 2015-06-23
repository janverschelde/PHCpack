with Communications_with_User;          use Communications_with_User;
with Standard_Natural_Numbers_io;       use Standard_Natural_Numbers_io;

package body File_Management is

-- INTERNAL DATA :

  input_file,output_file,wfile1,wfile2 : file_type;

-- OPERATIONS :

  procedure Silent_Open_Input_File ( filename : in string ) is
  begin
    Open(input_file,in_file,filename);
  end Silent_Open_Input_File;

  procedure Silent_Open_Input_File
              ( k : in natural32; filename : in string ) is
  begin
    if k = 1 then
      Open(wfile1,in_file,filename);
    elsif k = 2 then
      Open(wfile2,in_file,filename);
    end if;
  end Silent_Open_Input_File;

  procedure Silent_Open_Input_File is
  begin
    Read_Name_and_Open_File(input_file);
  end Silent_Open_Input_File;

  procedure Silent_Open_Input_File ( k : in natural32 ) is
  begin
    if k = 1 then
      Read_Name_and_Open_File(wfile1);
    elsif k = 2 then
      Read_Name_and_Open_File(wfile2);
    end if;
  end Silent_Open_Input_File;

  procedure Open_Input_File is
  begin
    put_line("Reading the name of the input file for solutions.");
    Read_Name_and_Open_File(input_file);
  end Open_Input_File;

  procedure Open_Input_File ( k : in natural32 ) is
  begin
    put("Reading the name of the input file for witness set ");
    put(k,1); put_line(".");
    if k = 1 then
      Read_Name_and_Open_File(wfile1);
    elsif k = 2 then
      Read_Name_and_Open_File(wfile2);
    end if;
  end Open_Input_File;

  procedure Create_Output_File is
  begin
    put_line("Reading the name of the output file for solutions.");
    Read_Name_and_Create_File(output_file);
  end Create_Output_File;

  function Solution_Input_File return file_type is
  begin
    return input_file;
  end Solution_Input_File;

  function Solution_Input_File ( k : natural32 ) return file_type is
  begin
    if k = 1 then
      return wfile1;
    elsif k = 2 then
      return wfile2;
    else
      return input_file;
    end if;
  end Solution_Input_File;

  function Solution_Output_File return file_type is
  begin
    return output_file;
  end Solution_Output_File;

  procedure Reset_Input_File ( k : in natural32 ) is
  begin
    if k = 1 then
      Reset(wfile1);
    elsif k = 2 then
      Reset(wfile2);
    end if;
  end Reset_Input_File;

  procedure Close_Input_File is
  begin
    close(input_file);
  end Close_Input_File;

  procedure Close_Input_File ( k : in natural32 ) is
  begin
    if k = 0 then
      close(input_file);
    elsif k = 1 then
      close(wfile1);
    elsif k = 2 then 
      close(wfile2);
    end if;
  end Close_Input_File;

  procedure Close_Output_File is
  begin
    close(output_file);
  end Close_Output_File;

end File_Management;
