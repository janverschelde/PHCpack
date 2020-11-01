with File_Scanning;                      use File_Scanning;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Complex_Solutions_io;      use Standard_Complex_Solutions_io;
with Standard_Dictionary_Solutions_io;

package body Main_Dictionary_Solutions is

  procedure Scan_Solutions
               ( filename : in string; sols : in out Solution_List;
                 solsonfile : out boolean; python_format : out boolean ) is

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
  begin
    if python_format
     then put(file,Length_Of(sols),natural32(Head_Of(sols).n),sols);
     else Standard_Dictionary_Solutions_io.put(file,sols);
    end if;
  end Write;

  procedure Write_Solutions
              ( outfilename : in string; python_format : in boolean;
                sols : in Solution_List ) is

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

  procedure Main ( infilename,outfilename : in string;
                   verbose : in integer32 := 0 ) is

    sols : Solution_List;
    solsonfile,python_format : boolean;

  begin
    if verbose > 0
     then put_line("-> in main_dictionary_solutions.Main ...");
    end if;
    Scan_Solutions(infilename,sols,solsonfile,python_format);
    if outfilename = "" then
      if python_format
       then put(Standard_Output,Length_Of(sols),
                natural32(Head_Of(sols).n),sols);
       else Standard_Dictionary_Solutions_io.put(sols);
      end if;
    elsif not Is_Null(sols) then
      Write_Solutions(outfilename,python_format,sols);
    end if;
  exception
    when others
      => new_line;
         put_line("Use as phc -x input output");
         new_line;
  end Main;

end Main_Dictionary_Solutions;
