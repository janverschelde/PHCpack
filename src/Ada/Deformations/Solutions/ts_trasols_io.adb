with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Complex_Poly_Systems;
with Standard_Complex_Laur_Systems;
with Standard_Complex_Solutions;
with Standard_Solution_Filters;          use Standard_Solution_Filters;
with DoblDobl_Complex_Poly_Systems;
with DoblDobl_Complex_Laur_Systems;
with DoblDobl_Complex_Solutions;
with DoblDobl_Solution_Filters;          use DoblDobl_Solution_Filters;
with QuadDobl_Complex_Poly_Systems;
with QuadDobl_Complex_Laur_Systems;
with QuadDobl_Complex_Solutions;
with QuadDobl_Solution_Filters;          use QuadDobl_Solution_Filters;
with Standard_Tracked_Solutions_io;
with DoblDobl_Tracked_Solutions_io;
with QuadDobl_Tracked_Solutions_io;
with Drivers_for_Failed_Paths;           use Drivers_for_Failed_Paths;

procedure ts_trasols_io is

-- DESCRIPTION :
--   Development of processing of the output file of a path tracker.

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
    Write_to_File(lq.all,failed);
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
    Write_to_File(lq.all,failed);
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
    Write_to_File(lq.all,failed);
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
    Write_to_File(lq.all,failed);
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
    Write_to_File(lq.all,failed);
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
    Write_to_File(lq.all,failed);
  end QuadDobl_Laur_Read;

  procedure Main_Test ( precision : in character ) is

  -- DESCRIPTION :
  --   For the precision equals to '0', '1', or '2',
  --   calls the test in double, double double, or quad double precision.

    infile : file_type;
    ans : character;

  begin
    new_line;
    put_line("Reading the name of a file with output of a tracker ...");
    Read_Name_and_Open_File(infile);
    new_line;
    put("Laurent polynomial systems ? (y/n) ");
    Ask_Yes_or_No(ans);
    new_line;
    put_line("Reading the contents of the file ...");
    case precision is
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
  end Main_Test;

  procedure Driver_Test ( precision : in character ) is

  -- DESCRIPTION :
  --   For the precision equals to '0', '1', or '2', tests the
  --   driver in double, double double, or quad double precision.

  begin
    case precision is
      when '0' => Standard_Scan_Failed_Paths("","");
      when '1' => DoblDobl_Scan_Failed_Paths("","");
      when '2' => QuadDobl_Scan_Failed_Paths("","");
      when others => null;
    end case;
  end Driver_Test;

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
    put("Type 0, 1, or 2 to select the precision : ");
    Ask_Alternative(prc,"012");
    new_line;
    put("Test the main driver ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans = 'y'
     then Driver_Test(prc);
     else Main_Test(prc);
    end if;
  end Main;

begin
  Main;
end ts_trasols_io;
