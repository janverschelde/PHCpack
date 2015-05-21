with text_io;                            use text_io;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Communications_with_User;           use Communications_with_User;
with Standard_Complex_Poly_Systems;
with Standard_Complex_Poly_Systems_io;   use Standard_Complex_Poly_Systems_io;
with Standard_Complex_Solutions;
with DoblDobl_Complex_Poly_Systems;
with DoblDobl_Complex_Poly_Systems_io;   use DoblDobl_Complex_Poly_Systems_io;
with DoblDobl_Complex_Solutions;
with QuadDobl_Complex_Poly_Systems;
with QuadDobl_Complex_Poly_Systems_io;   use QuadDobl_Complex_Poly_Systems_io;
with QuadDobl_Complex_Solutions;
with Drivers_for_Static_Lifting;         use Drivers_for_Static_Lifting;
with Drivers_for_MixedVol_Algorithm;     use Drivers_for_MixedVol_Algorithm;

procedure ts_drivstal is

-- DESCRIPTION :
--   This procedure calls the driver to static lifting.

  procedure Standard_Test is

  -- DESCRIPTION :
  --   Prompts the user for a polynomial system and calls the drivers
  --   to execute the polyhedral homotopies in standard double precision.

    use Standard_Complex_Poly_Systems,Standard_Complex_Solutions;

    lp,lq : Link_to_Poly_Sys;
    qsols,qsols0 : Solution_List;
    nt,mv,smv,tmv : natural32 := 0;
    file : file_type;
    ans : character;

  begin
    get(lp);
    new_line;
    put_line("Reading the name of the output file.");
    Read_Name_and_Create_File(file);
    lq := new Poly_Sys(lp'range);
    new_line;
    put("Use MixedVol algorithm ? (y/n) ");
    Ask_Yes_or_No(ans);
    new_line;
    put("Give the number of tasks (0 for no multitasking) : ");
    get(nt);
    if ans = 'y' then
      new_line;
      put("Use multiprecision arithmetic for Hermite normal form ? (y/n) ");
      Ask_Yes_or_No(ans);
      Driver_for_MixedVol_Algorithm
        (file,integer32(nt),lp.all,true,false,lq.all,qsols,qsols0,
         mv,smv,tmv,(ans = 'y'));
    else
      Driver_for_Mixed_Volume_Computation
        (file,integer32(nt),lp.all,true,lq.all,qsols,qsols0,mv,smv,tmv);
    end if;
  end Standard_Test;

  procedure DoblDobl_Test is

  -- DESCRIPTION :
  --   Prompts the user for a polynomial system and calls the drivers
  --   to execute the polyhedral homotopies in standard double precision.

    use DoblDobl_Complex_Poly_Systems,DoblDobl_Complex_Solutions;

    lp,lq : Link_to_Poly_Sys;
    qsols,qsols0 : Solution_List;
    nt,mv,smv,tmv : natural32 := 0;
    file : file_type;
    ans : character;

  begin
    get(lp);
    new_line;
    put_line("Reading the name of the output file.");
    Read_Name_and_Create_File(file);
    lq := new Poly_Sys(lp'range);
   -- new_line;
   -- put("Use MixedVol algorithm ? (y/n) ");
   -- Ask_Yes_or_No(ans);
   -- if ans = 'y' then
    new_line;
    put("Use multiprecision arithmetic for Hermite normal form ? (y/n) ");
    Ask_Yes_or_No(ans);
    new_line;
    put("Give the number of tasks (0 for no multitasking) : ");
    get(nt);
    Driver_for_MixedVol_Algorithm
      (file,integer32(nt),lp.all,true,false,lq.all,qsols,qsols0,
       mv,smv,tmv,(ans = 'y'));
   -- else
   --   Driver_for_Mixed_Volume_Computation
   --     (file,integer32(nt),lp.all,true,lq.all,qsols,qsols0,mv,smv,tmv);
   -- end if;
  end DoblDobl_Test;

  procedure QuadDobl_Test is

  -- DESCRIPTION :
  --   Prompts the user for a polynomial system and calls the drivers
  --   to execute the polyhedral homotopies in standard double precision.

    use QuadDobl_Complex_Poly_Systems,QuadDobl_Complex_Solutions;

    lp,lq : Link_to_Poly_Sys;
    qsols,qsols0 : Solution_List;
    nt,mv,smv,tmv : natural32 := 0;
    file : file_type;
    ans : character;

  begin
    get(lp);
    new_line;
    put_line("Reading the name of the output file.");
    Read_Name_and_Create_File(file);
    lq := new Poly_Sys(lp'range);
   -- new_line;
   -- put("Use MixedVol algorithm ? (y/n) ");
   -- Ask_Yes_or_No(ans);
   -- if ans = 'y' then
    new_line;
    put("Use multiprecision arithmetic for Hermite normal form ? (y/n) ");
    Ask_Yes_or_No(ans);
    new_line;
    put("Give the number of tasks (0 for no multitasking) : ");
    get(nt);
    Driver_for_MixedVol_Algorithm
      (file,integer32(nt),lp.all,true,false,lq.all,qsols,qsols0,
       mv,smv,tmv,(ans = 'y'));
   -- else
   --   Driver_for_Mixed_Volume_Computation
   --     (file,integer32(nt),lp.all,true,lq.all,qsols,qsols0,mv,smv,tmv);
   -- end if;
  end QuadDobl_Test;

  procedure Main is

  -- DESCRIPTION :
  --   The user is first prompted for the level of precision,
  --   before the system is entered.

    ans : character;

  begin
    new_line;
    put_line("MENU to call drivers for polyhedral homotopies :");
    put_line("  1. path tracking in standard double precision;");
    put_line("  2. path tracking in double double precision;");
    put_line("  3. path tracking in quad double precision.");
    put("Type 1, 2, or 3 to select the precision : ");
    Ask_Alternative(ans,"123");
    new_line;
    case ans is
      when '1' => Standard_Test;
      when '2' => DoblDobl_Test;
      when '3' => QuadDobl_Test;
      when others => null;
    end case;
  end Main;

begin
  Main;
end ts_drivstal;
