with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Complex_Poly_Systems;
with Standard_Complex_Laur_Systems;
with Standard_Complex_Solutions;
with Standard_Complex_Solutions_io;      use Standard_Complex_Solutions_io;
with Standard_Solution_Filters;          use Standard_Solution_Filters;
with DoblDobl_Complex_Poly_Systems;
with DoblDobl_Complex_Laur_Systems;
with DoblDobl_Complex_Solutions;
with DoblDobl_Complex_Solutions_io;      use DoblDobl_Complex_Solutions_io;
with DoblDobl_Solution_Filters;          use DoblDobl_Solution_Filters;
with QuadDobl_Complex_Poly_Systems;
with QuadDobl_Complex_Laur_Systems;
with QuadDobl_Complex_Solutions;
with QuadDobl_Complex_Solutions_io;      use QuadDobl_Complex_Solutions_io;
with QuadDobl_Solution_Filters;          use QuadDobl_Solution_Filters;
with Standard_Tracked_Solutions_io;
with DoblDobl_Tracked_Solutions_io;
with QuadDobl_Tracked_Solutions_io;

procedure ts_trasols_io is

-- DESCRIPTION :
--   Development of processing of the output file of a path tracker.

  procedure Write_to_File
              ( sols : in Standard_Complex_Solutions.Solution_List ) is

  -- DESCRIPTION :
  --   Prompts the user for a file name to write the solutions on.

    ans : character;
    file : file_type;
 
  begin
    put("Selected ");
    put(Standard_Complex_Solutions.Length_Of(sols),1);
    put_line(" failed solutions.");
    if not Standard_Complex_Solutions.Is_Null(sols) then
      new_line;
      put("Write start solutions corresponding to failed paths ");
      put("to file ? (y/n) ");
      Ask_Yes_or_No(ans);
      if ans = 'y' then
        new_line;
        put_line("Reading the name of an output file ...");
        Read_Name_and_Create_File(file);
        put(file,Standard_Complex_Solutions.Length_Of(sols),
            natural32(Standard_Complex_Solutions.Head_Of(sols).n),sols);
        close(file);
      end if;
    end if;
  end Write_to_File;

  procedure Write_to_File
              ( sols : in DoblDobl_Complex_Solutions.Solution_List ) is

  -- DESCRIPTION :
  --   Prompts the user for a file name to write the solutions on.

    ans : character;
    file : file_type;
 
  begin
    put("Selected ");
    put(DoblDobl_Complex_Solutions.Length_Of(sols),1);
    put_line(" failed solutions.");
    if not DoblDobl_Complex_Solutions.Is_Null(sols) then
      new_line;
      put("Write start solutions corresponding to failed paths ");
      put("to file ? (y/n) ");
      Ask_Yes_or_No(ans);
      if ans = 'y' then
        new_line;
        put_line("Reading the name of an output file ...");
        Read_Name_and_Create_File(file);
        put(file,DoblDobl_Complex_Solutions.Length_Of(sols),
            natural32(DoblDobl_Complex_Solutions.Head_Of(sols).n),sols);
        close(file);
      end if;
    end if;
  end Write_to_File;

  procedure Write_to_File
              ( sols : in QuadDobl_Complex_Solutions.Solution_List ) is

  -- DESCRIPTION :
  --   Prompts the user for a file name to write the solutions on.

    ans : character;
    file : file_type;
 
  begin
    put("Selected ");
    put(QuadDobl_Complex_Solutions.Length_Of(sols),1);
    put_line(" failed solutions.");
    if not QuadDobl_Complex_Solutions.Is_Null(sols) then
      new_line;
      put("Write start solutions corresponding to failed paths ");
      put("to file ? (y/n) ");
      Ask_Yes_or_No(ans);
      if ans = 'y' then
        new_line;
        put_line("Reading the name of an output file ...");
        Read_Name_and_Create_File(file);
        put(file,QuadDobl_Complex_Solutions.Length_Of(sols),
            natural32(QuadDobl_Complex_Solutions.Head_Of(sols).n),sols);
        close(file);
      end if;
    end if;
  end Write_to_File;

  procedure Standard_Poly_Read
              ( file : in file_type; verbose : in boolean := false ) is

  -- DESCRIPTION :
  --   Reads the target system, the start system, the start solution,
  --   and the solutions of the target system (in this order) from
  --   the file.  If verbose, then extra output is written to screen.
  --   The systems are ordinary polynomial systems and coefficients
  --   are read in standard double precision.

    lp,lq : Standard_Complex_Poly_Systems.Link_to_Poly_Sys;
    qsols,psols,failed : Standard_Complex_Solutions.Solution_List;
    tol : constant double_float := 1.0E-8;

  begin
    Standard_Tracked_Solutions_io.get(file,lp,lq,psols,qsols,verbose);
    failed := Select_Failed_Solutions(psols,qsols,tol,verbose);
    Write_to_File(failed);
  end Standard_Poly_Read;

  procedure Standard_Laur_Read
              ( file : in file_type; verbose : in boolean := false ) is

  -- DESCRIPTION :
  --   Reads the target system, the start system, the start solution,
  --   and the solutions of the target system (in this order) from
  --   the file.  If verbose, then extra output is written to screen.
  --   The systems are Laurent polynomial systems and coefficients
  --   are read in standard double precision.

    lp,lq : Standard_Complex_Laur_Systems.Link_to_Laur_Sys;
    qsols,psols,failed : Standard_Complex_Solutions.Solution_List;
    tol : constant double_float := 1.0E-8;

  begin
    Standard_Tracked_Solutions_io.get(file,lp,lq,psols,qsols,verbose);
    failed := Select_Failed_Solutions(psols,qsols,tol,verbose);
    Write_to_File(failed);
  end Standard_Laur_Read;

  procedure DoblDobl_Poly_Read
              ( file : in file_type; verbose : in boolean := false ) is

  -- DESCRIPTION :
  --   Reads the target system, the start system, the start solution,
  --   and the solutions of the target system (in this order) from
  --   the file.  If verbose, then extra output is written to screen.
  --   The systems are ordinary polynomial systems and coefficients
  --   are read in double double precision.

    lp,lq : DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
    qsols,psols,failed : DoblDobl_Complex_Solutions.Solution_List;
    tol : constant double_float := 1.0E-8;

  begin
    DoblDobl_Tracked_Solutions_io.get(file,lp,lq,psols,qsols,verbose);
    failed := Select_Failed_Solutions(psols,qsols,tol,verbose);
    Write_to_File(failed);
  end DoblDobl_Poly_Read;

  procedure DoblDobl_Laur_Read
              ( file : in file_type; verbose : in boolean := false ) is

  -- DESCRIPTION :
  --   Reads the target system, the start system, the start solution,
  --   and the solutions of the target system (in this order) from
  --   the file.  If verbose, then extra output is written to screen.
  --   The systems are Laurent polynomial systems and coefficients
  --   are read in double double precision.

    lp,lq : DoblDobl_Complex_Laur_Systems.Link_to_Laur_Sys;
    qsols,psols,failed : DoblDobl_Complex_Solutions.Solution_List;
    tol : constant double_float := 1.0E-8;

  begin
    DoblDobl_Tracked_Solutions_io.get(file,lp,lq,psols,qsols,verbose);
    failed := Select_Failed_Solutions(psols,qsols,tol,verbose);
    Write_to_File(failed);
  end DoblDobl_Laur_Read;

  procedure QuadDobl_Poly_Read
              ( file : in file_type; verbose : in boolean := false ) is

  -- DESCRIPTION :
  --   Reads the target system, the start system, the start solution,
  --   and the solutions of the target system (in this order) from
  --   the file.  If verbose, then extra output is written to screen.
  --   The systems are ordinary polynomial systems and coefficients
  --   are read in quad double precision.

    lp,lq : QuadDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
    qsols,psols,failed : QuadDobl_Complex_Solutions.Solution_List;
    tol : constant double_float := 1.0E-8;

  begin
    QuadDobl_Tracked_Solutions_io.get(file,lp,lq,psols,qsols,verbose);
    failed := Select_Failed_Solutions(psols,qsols,tol,verbose);
    Write_to_File(failed);
  end QuadDobl_Poly_Read;

  procedure QuadDobl_Laur_Read
              ( file : in file_type; verbose : in boolean := false ) is

  -- DESCRIPTION :
  --   Reads the target system, the start system, the start solution,
  --   and the solutions of the target system (in this order) from
  --   the file.  If verbose, then extra output is written to screen.
  --   The systems are Laurent polynomial systems and coefficients
  --   are read in quad double precision.

    lp,lq : QuadDobl_Complex_Laur_Systems.Link_to_Laur_Sys;
    qsols,psols,failed : QuadDobl_Complex_Solutions.Solution_List;
    tol : constant double_float := 1.0E-8;

  begin
    QuadDobl_Tracked_Solutions_io.get(file,lp,lq,psols,qsols,verbose);
    failed := Select_Failed_Solutions(psols,qsols,tol,verbose);
    Write_to_File(failed);
  end QuadDobl_Laur_Read;

  procedure Main is

  -- DESCRIPTION :
  --   Prompts the user for a file name, the output of a path tracker.

    infile : file_type;
    prc,ans : character;

  begin
    new_line;
    put_line("MENU for the working precision :");
    put_line("  0. standard double precision");
    put_line("  1. double double precision");
    put_line("  2. quad double precision");
    put("Type 0, 1, or 2 to select the precision :");
    Ask_Alternative(prc,"012");
    new_line;
    put_line("Reading the name of a file with output of a tracker ...");
    Read_Name_and_Open_File(infile);
    new_line;
    put("Laurent polynomial systems ? (y/n) ");
    Ask_Yes_or_No(ans);
    new_line;
    put_line("Reading the contents of the file ...");
    case prc is
      when '0' => 
        if ans = 'y'
         then Standard_Laur_Read(infile,true);
         else Standard_Poly_Read(infile,true);
        end if;
      when '1' => 
        if ans = 'y'
         then DoblDobl_Laur_Read(infile,true);
         else DoblDobl_Poly_Read(infile,true);
        end if;
      when '2' => 
        if ans = 'y'
         then QuadDobl_Laur_Read(infile,true);
         else QuadDobl_Poly_Read(infile,true);
        end if;
      when others => null;
    end case;
  end Main;

begin
  Main;
end ts_trasols_io;
