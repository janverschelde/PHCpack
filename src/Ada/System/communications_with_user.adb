package body Communications_with_User is

-- AUXILIARY :

  function Is_In ( s : string; ch : character ) return boolean is

  -- DESCRIPTION :
  --   Returns true if the character occurs in the string, false otherwise.

  begin
    for i in s'range loop
      if s(i) = ch
       then return true;
      end if;
    end loop;
    return false;
  end Is_In;

-- TARGET ROUTINES :

  procedure Ask ( ans : out character ) is

    ch : character;

  begin
    loop
      get(ch); skip_line;
      exit when Valid_Alternative(ch);
      put("Invalid alternative.  Please try again : ");
    end loop;
    ans := ch;
  end Ask;

  procedure Ask_Yes_or_No ( ans : out character ) is

    function Yes_or_No ( alt : character ) return boolean is
    begin
      if alt = 'y' or else alt = 'n'
       then return true;
       else return false;
      end if;
    end Yes_or_No;
    procedure Yes_or_No_Ask is new Ask (Yes_or_No);

  begin
    Yes_or_No_Ask(ans);
  end Ask_Yes_or_No;

  procedure Ask_Alternative ( ans : out character; alternatives : in string ) is

    function Is_Valid ( alt : character ) return boolean is
    begin
      return Is_In(alternatives,alt);
    end Is_Valid;
    procedure Ask_Alt is new Ask ( Is_Valid );

  begin
    Ask_Alt(ans);
  end Ask_Alternative;

  procedure Ask_Alternative
                ( ans : in out string; alternatives : string;
                  prefix : in character ) is

    ok : boolean := false;
    tmp : string(1..10);
    ind,cnt : natural;

  begin
    loop
      get_line(tmp,cnt);
      ans := "  ";
      ind := 1;
      while (ans(ans'first) = ' ') and (ind <= cnt) loop
        ans(ans'first) := tmp(ind);
        ind := ind+1;
      end loop;
      if ans(ans'first) = prefix then
        while (ans(ans'first+1) = ' ') and (ind <= cnt) loop
          ans(ans'first+1) := tmp(ind);
          ind := ind+1;
        end loop;
        if Is_In(alternatives,ans(ans'first+1))
         then ok := true;
         else put("Invalid alternative.  Please try again : ");
        end if;
      else
        if Is_In(alternatives,ans(ans'first))
         then ok := true;
         else put("Invalid alternative.  Please try again : ");
        end if;
      end if;
      exit when ok;
    end loop;
  end Ask_Alternative;

  procedure Report_Location_Error ( filename : in string ) is
  begin
    put("The file "); put(filename);
    put_line(" could not be found, please try again ...");
  end Report_Location_Error;

  procedure Read_Name_and_Open_File ( file : in out file_type ) is

    name : constant string := Read_String;

  begin
    Open(file,in_file,name);
  exception
    when NAME_ERROR => 
       Report_Location_Error(name);
       Read_Name_and_Open_File(file);
    when USE_ERROR =>
       put_line("File is not readable, please try again ...");
       Read_Name_and_Open_File(file);
  end Read_Name_and_Open_File;

  procedure Read_Name_and_Open_File 
              ( file : in out file_type; name : out Link_to_String ) is

    filename : constant string := Read_String;

  begin
    Open(file,in_file,filename);
    name := new string'(filename);
  exception
    when NAME_ERROR => 
       Report_Location_Error(filename);
       Read_Name_and_Open_File(file,name);
    when USE_ERROR =>
       put_line("File is not readable, please try again ...");
       Read_Name_and_Open_File(file,name);
  end Read_Name_and_Open_File;

  procedure Read_Name_and_Create_File ( file : in out file_type ) is

    filename : constant string := Read_String;
    ans : character;
    temp : file_type;

    procedure Retry is
    begin
      Create(file,out_file,filename);
    exception
      when USE_ERROR =>
        put_line("Could not create file, file already in use.");
        put_line("Please, try again...");
        Read_Name_and_Create_File(file);
      when NAME_ERROR =>
        put_line("Could not create file, perhaps wrong directory ?");
        put_line("Please, try again...");
        Read_Name_and_Create_File(file);
    end Retry;

  begin
    Open(temp,in_file,filename);
    Close(temp);
    put("There exists already a file named "); put_line(filename);
    put("Do you want to destroy this file ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans = 'y'
     then create(file,out_file,filename);
     else Read_Name_and_Create_File(file);
    end if;
  exception
    when others => Retry;
  end Read_Name_and_Create_File;

  procedure Read_Name_and_Create_File
              ( file : in out file_type; name : out Link_to_String ) is

    filename : constant string := Read_String;
    ans : character;
    temp : file_type;

    procedure Retry is
    begin
      Create(file,out_file,filename);
      name := new string'(filename);
    exception
      when USE_ERROR =>
        put_line("Could not create file, file already in use.");
        put_line("Please, try again...");
        Read_Name_and_Create_File(file,name);
      when NAME_ERROR =>
        put_line("Could not create file, perhaps wrong directory ?");
        put_line("Please, try again...");
        Read_Name_and_Create_File(file,name);
    end Retry;

  begin
    Open(temp,in_file,filename);
    Close(temp);
    put("There exists already a file named "); put_line(filename);
    put("Do you want to destroy this file ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans = 'y'
     then create(file,out_file,filename);
          name := new string'(filename);
     else Read_Name_and_Create_File(file,name);
    end if;
  exception
    when others => Retry;
  end Read_Name_and_Create_File;

  procedure Read_Name_and_Append_File ( file : in out file_type ) is

    name : constant string := Read_String;

  begin
    Open(file,append_file,name);
  exception
    when NAME_ERROR =>
       Report_Location_Error(name);
       Read_Name_and_Open_File(file);
    when USE_ERROR =>
       put_line("File is not readable, please try again ...");
       Read_Name_and_Open_File(file);
  end Read_Name_and_Append_File;

  procedure Open_Input_File
               ( file : in out file_type; filename : in string ) is
  begin
    Open(file,in_file,filename);
  exception
    when NAME_ERROR =>
       Report_Location_Error(filename);
       Read_Name_and_Open_File(file);
    when USE_ERROR =>
       put("The file "); put(filename);
       put_line(" is not readable, please try again ...");
       Read_Name_and_Open_File(file);
  end Open_Input_File;

  procedure Open_Input_File
               ( file : in out file_type; filename : in string;
                 name : out Link_to_String ) is
  begin
    Open(file,in_file,filename);
    name := new string'(filename);
  exception
    when NAME_ERROR =>
       Report_Location_Error(filename);
       Read_Name_and_Open_File(file,name);
    when USE_ERROR =>
       put("The file "); put(filename);
       put_line(" is not readable, please try again ...");
       Read_Name_and_Open_File(file,name);
  end Open_Input_File;

  procedure Create_Output_File
                 ( file : in out file_type; filename : in string ) is

    ans : character;
    temp : file_type;

    procedure Retry is
    begin
      Create(file,out_file,filename);
    exception
      when USE_ERROR =>
        put("Could not create file "); put(filename);
        put_line(", file already in use.");
        put_line("Please, try again ...");
        Read_Name_and_Create_File(file);
      when NAME_ERROR =>
        put("Could not create file "); put(filename);
        put_line(", perhaps wrong directory ?");
        put_line("Please, try again ...");
        Read_Name_and_Create_File(file);
    end Retry;

  begin
    if filename = "" then
      new_line;
      put_line("Reading the name of the output file.");
      Read_Name_and_Create_File(file);
    else
      Open(temp,in_file,filename); Close(temp);
      new_line;
      put("There exists already a file named "); put_line(filename);
      put("Do you want to destroy this file ? (y/n) ");
      Ask_Yes_or_No(ans);
      if ans = 'y'
       then create(file,out_file,filename);
       else Read_Name_and_Create_File(file);
      end if;
    end if;
  exception
    when others => Retry;
  end Create_Output_File;

  procedure Create_Output_File
                 ( file : in out file_type; filename : in string;
                   name : out Link_to_String ) is

    ans : character;
    temp : file_type;

    procedure Retry is
    begin
      Create(file,out_file,filename);
      name := new string'(filename);
    exception
      when USE_ERROR =>
        put("Could not create file "); put(filename);
        put_line(", file already in use.");
        put_line("Please, try again ...");
        Read_Name_and_Create_File(file,name);
      when NAME_ERROR =>
        put("Could not create file "); put(filename);
        put_line(", perhaps wrong directory ?");
        put_line("Please, try again ...");
        Read_Name_and_Create_File(file,name);
    end Retry;

  begin
    if filename = "" then
      new_line;
      put_line("Reading the name of the output file.");
      Read_Name_and_Create_File(file,name);
    else
      Open(temp,in_file,filename); Close(temp);
      new_line;
      put("There exists already a file named ");
      put_line(filename);
      put("Do you want to destroy this file ? (y/n) ");
      Ask_Yes_or_No(ans);
      if ans = 'y' then
        create(file,out_file,filename);
        name := new string'(filename);
      else
        Read_Name_and_Create_File(file,name);
      end if;
    end if;
  exception
    when others => Retry;
  end Create_Output_File;

  procedure Open_Append_File
               ( file : in out file_type; filename : in string ) is
  begin
    Open(file,append_file,filename);
  exception
    when NAME_ERROR =>
       Report_Location_Error(filename);
       Read_Name_and_Append_File(file);
    when USE_ERROR =>
       put("The file "); put(filename);
       put_line(" is not readable, please try again ...");
       Read_Name_and_Append_File(file);
  end Open_Append_File;

  function Prompt_for_Precision return character is

    res : character;

  begin
    new_line;
    put_line("MENU for the working precision :");
    put_line("  0. hardware double precision;");
    put_line("  1. double double precision; or");
    put_line("  2. quad double precision.");
    put("Type 0, 1, or 2 to select the precision : ");
    Ask_Alternative(res,"012");
    return res;
  end Prompt_for_Precision;

  procedure End_of_Input_Message is
  begin
    new_line;
    put_line("No more input expected.  See output file for results.");
    new_line;
  end End_of_Input_Message;

end Communications_with_User;
