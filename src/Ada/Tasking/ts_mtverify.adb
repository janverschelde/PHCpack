with text_io;                           use text_io;
with Communications_with_User;          use Communications_with_User;
with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;       use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;         use Standard_Floating_Numbers;
with Standard_Complex_Poly_Systems;
with Standard_Complex_Laur_Systems;
with Standard_Complex_Laur_Systems_io;  use Standard_Complex_Laur_Systems_io;
with Standard_Laur_Poly_Convertors;
with DoblDobl_Complex_Poly_Systems;
with DoblDobl_Complex_Laur_Systems;
with DoblDobl_Complex_Laur_Systems_io;  use DoblDobl_Complex_Laur_Systems_io;
with DoblDobl_Laur_Poly_Convertors;
with QuadDobl_Complex_Poly_Systems;
with QuadDobl_Complex_Laur_Systems;
with QuadDobl_Complex_Laur_Systems_io;  use QuadDobl_Complex_Laur_Systems_io;
with QuadDobl_Laur_Poly_Convertors;
with Standard_Complex_Solutions;
with Standard_System_and_Solutions_io;
with DoblDobl_Complex_Solutions;
with DoblDobl_System_and_Solutions_io;
with QuadDobl_Complex_Solutions;
with QuadDobl_System_and_Solutions_io;
with Root_Refining_Parameters;          use Root_Refining_Parameters;
with Multitasking_Root_Refiners;        use Multitasking_Root_Refiners;

procedure ts_mtverify is

-- DESCRIPTION :
--   Verification of a list of solutions with multitasking.

  procedure Standard_Verify is

  -- DESCRIPTION :
  --   Applies the verification in standard double precision arithmetic.

    use Standard_Complex_Poly_Systems;
    use Standard_Complex_Laur_Systems;
    use Standard_Complex_Solutions;

    lq : Link_to_Laur_Sys;
    sols : Solution_List;
    nt : integer32 := 0;
    file : file_type;
    ans : character;
    epsxa,epsfa,tolsing : double_float;
    numbit,maxit : natural32 := 0;
    deflate,wout,reporting : boolean;

  begin
    Standard_System_and_Solutions_io.get(lq,sols);
    new_line;
    put("Read "); put(Length_Of(sols),1); put_line(" solutions.");
    new_line;
    put_line("Reading the name of the output file ...");
    Read_Name_and_Create_File(file);
    put(file,lq'last,1); new_line(file); put(file,lq.all);
    new_line;
    put("Give the number of tasks : "); get(nt);
    Standard_Default_Root_Refining_Parameters
      (epsxa,epsfa,tolsing,maxit,deflate,wout);
    Standard_Menu_Root_Refining_Parameters
      (file,epsxa,epsfa,tolsing,maxit,deflate,wout);
    new_line;
    put("Do you want output to screen during refinement ? (y/n) ");
    Ask_Yes_or_No(ans);
    reporting := (ans = 'y');
    new_line;
    if Standard_Laur_Poly_Convertors.Is_Genuine_Laurent(lq.all) then
      put_line("The system is a genuine Laurent polynomial system.");
      if reporting then
        Reporting_Multitasking_Root_Refiner
          (file,nt,lq.all,sols,epsxa,epsfa,tolsing,numbit,maxit,deflate);
      else
        Silent_Multitasking_Root_Refiner
          (file,nt,lq.all,sols,epsxa,epsfa,tolsing,numbit,maxit,deflate);
      end if;
    else
      put_line("The system is a positive Laurent polynomial system.");
      declare
        use Standard_Laur_Poly_Convertors;
        p : Poly_Sys(lq'range) := Positive_Laurent_Polynomial_System(lq.all);
      begin
        if reporting then
          Reporting_Multitasking_Root_Refiner
            (file,nt,p,sols,epsxa,epsfa,tolsing,numbit,maxit,deflate);
        else
          Silent_Multitasking_Root_Refiner
            (file,nt,p,sols,epsxa,epsfa,tolsing,numbit,maxit,deflate);
        end if;
        Clear(p);
      end;
    end if;
  end Standard_Verify;

  procedure DoblDobl_Verify is

  -- DESCRIPTION :
  --   Applies the verification in double double precision arithmetic.

    use DoblDobl_Complex_Poly_Systems;
    use DoblDobl_Complex_Laur_Systems;
    use DoblDobl_Complex_Solutions;

    lq : Link_to_Laur_Sys;
    sols : Solution_List;
    nt : integer32 := 0;
    file : file_type;
    ans : character;
    epsxa,epsfa,tolsing : double_float;
    numbit,maxit : natural32 := 0;
    deflate,wout,reporting : boolean;

  begin
    DoblDobl_System_and_Solutions_io.get(lq,sols);
    new_line;
    put("Read "); put(Length_Of(sols),1); put_line(" solutions.");
    new_line;
    put_line("Reading the name of the output file ...");
    Read_Name_and_Create_File(file);
    put(file,lq'last,1); new_line(file); put(file,lq.all);
    new_line;
    put("Give the number of tasks : "); get(nt);
    Standard_Default_Root_Refining_Parameters
      (epsxa,epsfa,tolsing,maxit,deflate,wout);
    Standard_Menu_Root_Refining_Parameters
      (file,epsxa,epsfa,tolsing,maxit,deflate,wout);
    new_line;
    put("Do you want output to screen during refinement ? (y/n) ");
    Ask_Yes_or_No(ans);
    reporting := (ans = 'y');
    new_line;
    if DoblDobl_Laur_Poly_Convertors.Is_Genuine_Laurent(lq.all) then
      put_line("The system is a genuine Laurent polynomial system.");
      if reporting then
        Reporting_Multitasking_Root_Refiner
          (file,nt,lq.all,sols,epsxa,epsfa,tolsing,numbit,maxit,deflate);
      else
        Silent_Multitasking_Root_Refiner
          (file,nt,lq.all,sols,epsxa,epsfa,tolsing,numbit,maxit,deflate);
      end if;
    else
      put_line("The system is a positive Laurent polynomial system.");
      declare
        use DoblDobl_Laur_Poly_Convertors;
        p : Poly_Sys(lq'range) := Positive_Laurent_Polynomial_System(lq.all);
      begin
        if reporting then
          Reporting_Multitasking_Root_Refiner
            (file,nt,p,sols,epsxa,epsfa,tolsing,numbit,maxit,deflate);
        else
          Silent_Multitasking_Root_Refiner
            (file,nt,p,sols,epsxa,epsfa,tolsing,numbit,maxit,deflate);
        end if;
        Clear(p);
      end;
    end if;
  end DoblDobl_Verify;

  procedure QuadDobl_Verify is

  -- DESCRIPTION :
  --   Applies the verification in double double precision arithmetic.

    use QuadDobl_Complex_Poly_Systems;
    use QuadDobl_Complex_Laur_Systems;
    use QuadDobl_Complex_Solutions;

    lq : Link_to_Laur_Sys;
    sols : Solution_List;
    nt : integer32 := 0;
    file : file_type;
    ans : character;
    epsxa,epsfa,tolsing : double_float;
    numbit,maxit : natural32 := 0;
    deflate,wout,reporting : boolean;

  begin
    QuadDobl_System_and_Solutions_io.get(lq,sols);
    new_line;
    put("Read "); put(Length_Of(sols),1); put_line(" solutions.");
    new_line;
    put_line("Reading the name of the output file ...");
    Read_Name_and_Create_File(file);
    put(file,lq'last,1); new_line(file); put(file,lq.all);
    new_line;
    put("Give the number of tasks : "); get(nt);
    Standard_Default_Root_Refining_Parameters
      (epsxa,epsfa,tolsing,maxit,deflate,wout);
    Standard_Menu_Root_Refining_Parameters
      (file,epsxa,epsfa,tolsing,maxit,deflate,wout);
    new_line;
    put("Do you want output to screen during refinement ? (y/n) ");
    Ask_Yes_or_No(ans);
    reporting := (ans = 'y');
    new_line;
    if QuadDobl_Laur_Poly_Convertors.Is_Genuine_Laurent(lq.all) then
      put_line("The system is a genuine Laurent polynomial system.");
      if reporting then
        Reporting_Multitasking_Root_Refiner
          (file,nt,lq.all,sols,epsxa,epsfa,tolsing,numbit,maxit,deflate);
      else
        Silent_Multitasking_Root_Refiner
          (file,nt,lq.all,sols,epsxa,epsfa,tolsing,numbit,maxit,deflate);
      end if;
    else
      put_line("The system is a positive Laurent polynomial system.");
      declare
        use QuadDobl_Laur_Poly_Convertors;
        p : Poly_Sys(lq'range) := Positive_Laurent_Polynomial_System(lq.all);
      begin
        if reporting then
          Reporting_Multitasking_Root_Refiner
            (file,nt,p,sols,epsxa,epsfa,tolsing,numbit,maxit,deflate);
        else
          Silent_Multitasking_Root_Refiner
            (file,nt,p,sols,epsxa,epsfa,tolsing,numbit,maxit,deflate);
        end if;
        Clear(p);
      end;
    end if;
  end QuadDobl_Verify;

  procedure Main is

  -- DESCRIPTION :
  --   Prompts the user for the level of precision
  --   and then calls the proper verification procedure.

    ans : character;

  begin
    new_line;
    put_line("MENU to set the precision in the multitasking root refiner.");
    put_line("  0. standard double precision; or");
    put_line("  1. double double precision; or");
    put_line("  2. quad double precision.");
    put("Type 0, 1, or 2 to select the precision : ");
    Ask_Alternative(ans,"012");
    case ans is
      when '0' => Standard_Verify;
      when '1' => DoblDobl_Verify;
      when '2' => QuadDobl_Verify;
      when others => null;
    end case;
  end Main;

begin
  Main;
end ts_mtverify;
