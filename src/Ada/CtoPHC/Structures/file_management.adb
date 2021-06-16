with Communications_with_User;          use Communications_with_User;
with Standard_Natural_Numbers_io;       use Standard_Natural_Numbers_io;

package body File_Management is

-- INTERNAL DATA :

  link_to_infile,link_to_outfile : link_to_file_type;
  link_to_wfile1,link_to_wfile2 : link_to_file_type;

-- OPERATIONS :

  procedure Silent_Open_Input_File ( filename : in string ) is
  begin
    link_to_infile := new file_type;
    Open(link_to_infile.all,in_file,filename);
  end Silent_Open_Input_File;

  procedure Silent_Open_Input_File
              ( k : in natural32; filename : in string ) is
  begin
    if k = 1 then
      link_to_wfile1 := new file_type;
      Open(link_to_wfile1.all,in_file,filename);
    elsif k = 2 then
      link_to_wfile2 := new file_type;
      Open(link_to_wfile2.all,in_file,filename);
    end if;
  end Silent_Open_Input_File;

  procedure Silent_Open_Input_File is
  begin
    link_to_infile := new file_type;
    Read_Name_and_Open_File(link_to_infile.all);
  end Silent_Open_Input_File;

  procedure Silent_Open_Input_File ( k : in natural32 ) is
  begin
    if k = 1 then
      link_to_wfile1 := new file_type;
      Read_Name_and_Open_File(link_to_wfile1.all);
    elsif k = 2 then
      link_to_wfile2 := new file_type;
      Read_Name_and_Open_File(link_to_wfile2.all);
    end if;
  end Silent_Open_Input_File;

  procedure Open_Input_File is
  begin
    link_to_infile := new file_type;
    put_line("Reading the name of the input file...");
    Read_Name_and_Open_File(link_to_infile.all);
  end Open_Input_File;

  procedure Open_Input_File ( k : in natural32 ) is
  begin
    put("Reading the name of the input file for witness set ");
    put(k,1); put_line(".");
    if k = 1 then
      link_to_wfile1 := new file_type;
      Read_Name_and_Open_File(link_to_wfile1.all);
    elsif k = 2 then
      link_to_wfile2 := new file_type;
      Read_Name_and_Open_File(link_to_wfile2.all);
    end if;
  end Open_Input_File;

  procedure Create_Output_File is
  begin
    link_to_outfile := new file_type;
    put_line("Reading the name of the output file...");
    Read_Name_and_Create_File(link_to_outfile.all);
  end Create_Output_File;

  function Link_to_Input return link_to_file_type is
  begin
    return link_to_infile;
  end Link_to_Input;

  function Link_to_Input ( k : natural32 ) return link_to_file_type is
  begin
    if k = 1 then
      return link_to_wfile1;
    elsif k = 2 then
      return link_to_wfile2;
    else
      return link_to_infile;
    end if;
  end Link_to_Input;

  function Link_to_Output return link_to_file_type is
  begin
    return link_to_outfile;
  end Link_to_Output;

  procedure Reset_Input_File ( k : in natural32 ) is
  begin
    if k = 1 then
      Reset(link_to_wfile1.all);
    elsif k = 2 then
      Reset(link_to_wfile2.all);
    end if;
  end Reset_Input_File;

  procedure Close_Input_File is
  begin
    close(link_to_infile.all);
  end Close_Input_File;

  procedure Close_Input_File ( k : in natural32 ) is
  begin
    if k = 0 then
      close(link_to_infile.all);
    elsif k = 1 then
      close(link_to_wfile1.all);
    elsif k = 2 then 
      close(link_to_wfile2.all);
    end if;
  end Close_Input_File;

  procedure Close_Output_File is
  begin
    close(link_to_outfile.all);
  end Close_Output_File;

end File_Management;
