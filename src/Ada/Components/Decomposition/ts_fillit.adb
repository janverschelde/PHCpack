with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with File_Scanning;                      use File_Scanning;
with String_Splitters;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Characters_and_Numbers;             use Characters_and_Numbers;
with Standard_Complex_Poly_Systems;      use Standard_Complex_Poly_Systems;
with Standard_Complex_Poly_Systems_io;   use Standard_Complex_Poly_Systems_io;
with Standard_Complex_Solutions;         use Standard_Complex_Solutions;
with Standard_Complex_Solutions_io;      use Standard_Complex_Solutions_io;
with Witness_Sets_io;                    use Witness_Sets_io;
-- with Homotopy_Membership_Tests;          use Homotopy_Membership_Tests;

procedure ts_fillit is

-- DESCRIPTION :
--   This procedure is the interactive development of the filtration
--   of the junk from the list of candidate witness points.
--   The name of the procedure is a contraction of filter and split,
--   as the output is spread over as many files as the dimension + 1.

  type File_Array is array ( integer32 range <> ) of file_type;
  type Array_of_Solution_Lists is
    array ( integer32 range <> ) of Solution_List;

  procedure Output_File_Creation
              ( k : in integer32; fa : in out File_Array ) is

  -- DESCRIPTION :
  --   Creates interactively a sequence of k+1 files, adding numbers
  --   to the file name the user gives.

  -- REQUIRED :
  --   The range of the file array is 0..k.

  begin
    new_line;
    put_line("Reading the name of the sequence of output files.");
    declare
      name : constant string := String_Splitters.Read_String;
    begin
      for i in 0..k loop
        declare
          sn : constant string := Convert(i);
          fn : constant string := name & sn;
        begin 
          put("  creating the file ");
          put(fn); put_line(" ...");
          Create_Output_File(fa(i),fn);
        end;
      end loop;
    end;
  end Output_file_Creation;

  procedure Write_Solutions ( file : in file_type;
                              sols : in Solution_List ) is
  begin
    if not Is_Null(sols) then
      new_line(file);
      put_line(file,"THE SOLUTIONS :");
      new_line(file);
      put(file,Length_Of(sols),natural32(Head_Of(sols).n),sols);
    end if;
  end Write_Solutions;

  procedure Remove_Junk ( logfile : in file_type; k : in integer32;
                          files : in out File_Array;
                          embsys : in Array_of_Poly_Sys;
                          embsols : in Array_of_Solution_Lists ) is

   -- restol : constant double_float := 1.0E-10;
   -- homtol : constant double_float := 1.0E-6;
 
  begin
    put_line(files(k),embsys(k).all);
    Write_Solutions(files(k),embsols(k));
    close(files(k));
    for i in reverse 0..k-1 loop
      put_line(files(i),embsys(i).all);
      Write_Solutions(files(i),embsols(i));
      close(files(i));
    end loop;
  end Remove_Junk;

  procedure Main is

  -- DESCRIPTION :
  --   This procedures reads the input file, creates the output files,
  --   and then launches the filtering procedure.

    k : integer32 := 0;
    infile,logfile : file_type;
    lp : Link_to_Poly_Sys;
    sols : Solution_List;

  begin
    new_line;
    put_line("Filter junk points from the candidate witness point lists.");
    new_line;
    put_line("Reading the name of the output log file for diagnostics.");
    Read_Name_and_Create_File(logfile);
    new_line;
    put_line("Reading the name of the file for the embedding sequence.");
    Read_Name_and_Open_File(infile);
    Standard_Read_Embedding(infile,lp,sols,natural32(k));
    declare
      embsys : Array_of_Poly_Sys(0..k);
      embsols : Array_of_Solution_Lists(0..k);
      found : boolean;
      files : File_Array(0..k); 
    begin
      embsys(k) := lp;
      embsols(k) := sols;
      for i in reverse 0..k-1 loop
        get(infile,embsys(i));
        Scan_and_Skip(infile,"SOLUTIONS",found);
        if found
         then get(infile,embsols(i));
        end if;
      end loop;
      Output_File_Creation(k,files);
      Remove_Junk(logfile,k,files,embsys,embsols);
    end;
  end Main;

begin
  Main;
end ts_fillit;
