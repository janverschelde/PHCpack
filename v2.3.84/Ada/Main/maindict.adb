with text_io;                            use text_io;
with File_Scanning;                      use File_Scanning;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Complex_Solutions;         use Standard_Complex_Solutions;
with Standard_Complex_Solutions_io;      use Standard_Complex_Solutions_io;
with Standard_Dictionary_Solutions_io;

procedure maindict ( infilename,outfilename : in string ) is

-- READING THE INPUT :

  procedure Scan_Solutions
               ( filename : in string; sols : in out Solution_List;
                 solsonfile : out boolean; python_format : out boolean ) is

  -- DESCRIPTION :
  --   Checks whether the given file name corresponds to a file with
  --   the solutions in a correct format.  If so, then solsonfile is
  --   true on return and sols contain the solutions.
  --   When the solutions are in Python format, then python_format is true.

    file : file_type;
    found : boolean;

  begin
    solsonfile := false;
    if filename /= "" then
      declare
      begin
        Open(file,in_file,filename);
      exception
        when others =>
          put("Could not open the file with name ");
          put(filename); put_line("."); return;
      end;
      Scan_and_Skip(file,"THE SOLUTIONS",found);
      if found then
        try_get(file,sols);
        python_format := false;
      else
        Reset(file);
        Standard_Dictionary_Solutions_io.get(file,sols);
        python_format := true;
      end if;
      Close(file);
      solsonfile := (Length_Of(sols) > 0);
    end if;
  exception
    when others =>
      put("Something wrong with solutions on "); put(filename);
      put_line("."); return;
  end Scan_Solutions;

  procedure Write ( file : in file_type; python_format : in boolean;
                    sols : in Solution_List ) is

  -- DESCRIPTION :
  --   Writes the solutions to the file.

  -- REQUIRED : file is opened in the right mode.

  begin
    if python_format
     then put(file,Length_Of(sols),natural32(Head_Of(sols).n),sols);
     else Standard_Dictionary_Solutions_io.put(file,sols);
    end if;
  end Write;

  procedure Write_Solutions ( python_format : in boolean;
                              sols : in Solution_List ) is

  -- DESCRIPTION :
  --   If the file with name outfilename does not exist,
  --   then it is created and used to write the solutions on;
  --   otherwise, the solutions are appended to the file with
  --   name outfilename.

    temp,file : file_type;

  begin
    Open(temp,in_file,outfilename); Close(temp);
    Open(file,append_file,outfilename);
    new_line(file);
    put_line(file,"THE SOLUTIONS :"); 
    Write(file,python_format,sols); Close(file);
  exception
    when others => Create(file,out_file,outfilename);
                   Write(file,python_format,sols); Close(file);
  end Write_Solutions;

  procedure Main is

    sols : Solution_List;
    solsonfile,python_format : boolean;

  begin
    Scan_Solutions(infilename,sols,solsonfile,python_format);
    if outfilename = "" then
      if python_format
       then put(Standard_Output,Length_Of(sols),
                natural32(Head_Of(sols).n),sols);
       else Standard_Dictionary_Solutions_io.put(sols);
      end if;
    elsif not Is_Null(sols) then
      Write_Solutions(python_format,sols);
    end if;
  end Main;

begin
  Main;
end maindict;
