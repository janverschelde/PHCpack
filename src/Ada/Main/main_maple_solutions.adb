with File_Scanning;                      use File_Scanning;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
-- with Maple, the default format uses multiprecision numbers
--with Standard_Complex_Solutions;         use Standard_Complex_Solutions;
--with Standard_Complex_Solutions_io;      use Standard_Complex_Solutions_io;
--with Standard_Maple_Solutions_io;
with Multprec_Complex_Solutions_io;      use Multprec_Complex_Solutions_io;
with Multprec_Maple_Solutions_io;

package body Main_Maple_Solutions is

  procedure Scan_Solutions
              ( filename : in string; sols : in out Solution_List;
                solsonfile : out boolean; maple_format : out boolean ) is

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
        get(file,sols);
        maple_format := false;
      else
        Reset(file);
        Multprec_Maple_Solutions_io.get(file,sols);
       -- Standard_Maple_Solutions_io.get(file,sols);
        maple_format := true;
      end if;
      Close(file);
      solsonfile := (Length_Of(sols) > 0);
    end if;
  exception
    when others =>
      put("Something wrong with solutions on the file '"); put(filename);
      put_line("'.");
      if maple_format then
        put_line("The solution were expected to be in Maple format!?");
      else  
        put_line("The solution were expected to be in PHCpack format!?");
      end if;
      put("Solution in PHCpack format must start with ");
      put_line("'THE SOLUTIONS :'");
      put_line("on the first line of the inpt file.");
      return;
  end Scan_Solutions;

  procedure Write ( file : in file_type; maple_format : in boolean;
                    sols : in Solution_List ) is
  begin
    if maple_format
     then put(file,Length_Of(sols),natural32(Head_Of(sols).n),sols);
     else Multprec_Maple_Solutions_io.put(file,sols);
    -- else Standard_Maple_Solutions_io.put(file,sols);
    end if;
  end Write;

  procedure Write_Solutions
              ( outfilename : in string; maple_format : in boolean;
                sols : in Solution_List ) is

    temp,file : file_type;

  begin
    Open(temp,in_file,outfilename); Close(temp);
    Open(file,append_file,outfilename);
    new_line(file);
    put_line(file,"THE SOLUTIONS :"); 
    Write(file,maple_format,sols); Close(file);
  exception
    when others => Create(file,out_file,outfilename);
                   Write(file,maple_format,sols); Close(file);
  end Write_Solutions;

  procedure Main ( infilename,outfilename : in string;
                   verbose : in integer32 := 0 ) is

    sols : Solution_List;
    solsonfile,maple_format : boolean;

  begin
    if verbose > 0 then
      put("At verbose level "); put(verbose,1);
      put_line(", in main_maple_solutions.Main ...");
    end if;
    Scan_Solutions(infilename,sols,solsonfile,maple_format);
    if outfilename = "" then
      if maple_format
       then put(Standard_Output,Length_Of(sols),
                natural32(Head_Of(sols).n),sols);
       else Multprec_Maple_Solutions_io.put(sols);
      -- else Standard_Maple_Solutions_io.put(sols);
      end if;
    elsif not Is_Null(sols) then
      Write_Solutions(outfilename,maple_format,sols);
    end if;
  exception
    when others
      => new_line;
         put_line("Use as phc -z input output");
         new_line;
  end Main;

end Main_Maple_Solutions;
