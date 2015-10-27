with text_io;                           use text_io;
with Communications_with_User;          use Communications_with_User;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Complex_Poly_Systems;
with DoblDobl_Complex_Poly_Systems;
with QuadDobl_Complex_Poly_Systems;
with Standard_Complex_to_Real_Poly;
with DoblDobl_Complex_to_Real_Poly;
with QuadDobl_Complex_to_Real_Poly;
with Standard_Complex_Solutions;
with DoblDobl_Complex_Solutions;
with QuadDobl_Complex_Solutions;
with Standard_Parameter_Systems;
with DoblDobl_Parameter_Systems;
with QuadDobl_Parameter_Systems;
with Parameter_Homotopy_Continuation;   use Parameter_Homotopy_Continuation;

procedure ts_parcon is

-- DESCRIPTION :
--   Stand-alone calling procedure to coefficient parameter continuation.

  procedure Standard_Main is

  -- DESCRIPTION :
  --   Test on parameter continuation in standard double precision.

    use Standard_Complex_Poly_Systems;
    use Standard_Complex_Solutions;
    use Standard_Parameter_Systems;

    file : file_type;
    lp : Link_to_Poly_Sys;
    sols : Solution_List;
    nb_equ,nb_unk,nb_par : integer32;
    ans : character;
    isreal : boolean;

  begin
    new_line;
    ans := Parameter_Homotopy_Continuation.Show_Menu;
    new_line;
    if ans = '1' then
      put_line("Running coefficient-parameter homotopy continuation...");
      Read_Parameter_Homotopy(file,lp,sols,nb_equ,nb_unk,nb_par);
      Coefficient_Parameter_Homotopy_Continuation
        (file,lp.all,sols,nb_equ,nb_unk,nb_par);
    else
      put_line("Running real sweep to target or first singularity...");
      Read_Parameter_Homotopy(file,lp,sols,nb_equ,nb_unk,nb_par);
      isreal := Standard_Complex_to_Real_Poly.Is_Real(lp.all);
      Sweep(file,isreal,lp.all,sols,nb_equ,nb_unk,nb_par);
    end if;
  end Standard_Main;

  procedure DoblDobl_Main is

  -- DESCRIPTION :
  --   Test on parameter continuation in double double precision.

    use DoblDobl_Complex_Poly_Systems;
    use DoblDobl_Complex_Solutions;
    use DoblDobl_Parameter_Systems;

    file : file_type;
    lp : Link_to_Poly_Sys;
    sols : Solution_List;
    nb_equ,nb_unk,nb_par : integer32;
    ans : character;
    isreal : boolean;

  begin
    new_line;
    ans := Parameter_Homotopy_Continuation.Show_Menu;
    new_line;
    if ans = '1' then
      put_line("Running coefficient-parameter homotopy continuation...");
      Read_Parameter_Homotopy(file,lp,sols,nb_equ,nb_unk,nb_par);
      Coefficient_Parameter_Homotopy_Continuation
        (file,lp.all,sols,nb_equ,nb_unk,nb_par);
    else
      put_line("Running real sweep to target or first singularity...");
      Read_Parameter_Homotopy(file,lp,sols,nb_equ,nb_unk,nb_par);
      isreal := DoblDobl_Complex_to_Real_Poly.Is_Real(lp.all);
      Sweep(file,isreal,lp.all,sols,nb_equ,nb_unk,nb_par);
    end if;
  end DoblDobl_Main;

  procedure QuadDobl_Main is

  -- DESCRIPTION :
  --   Test on parameter continuation in double double precision.

    use QuadDobl_Complex_Poly_Systems;
    use QuadDobl_Complex_Solutions;
    use QuadDobl_Parameter_Systems;

    file : file_type;
    lp : Link_to_Poly_Sys;
    sols : Solution_List;
    nb_equ,nb_unk,nb_par : integer32;
    ans : character;
    isreal : boolean;

  begin
    new_line;
    ans := Parameter_Homotopy_Continuation.Show_Menu;
    new_line;
    if ans = '1' then
      put_line("Running coefficient-parameter homotopy continuation...");
      Read_Parameter_Homotopy(file,lp,sols,nb_equ,nb_unk,nb_par);
      Coefficient_Parameter_Homotopy_Continuation
        (file,lp.all,sols,nb_equ,nb_unk,nb_par);
    else
      put_line("Running real sweep to target or first singularity...");
      Read_Parameter_Homotopy(file,lp,sols,nb_equ,nb_unk,nb_par);
      isreal := quadDobl_Complex_to_Real_Poly.Is_Real(lp.all);
      Sweep(file,isreal,lp.all,sols,nb_equ,nb_unk,nb_par);
    end if;
  end QuadDobl_Main;

  procedure Main is

  -- DESCRIPTION :
  --   Prompts the user for the level of precision
  --   and then calls the corresponding test.

    ans : character;

  begin
    new_line;
    put_line("MENU to select the working precision :");
    put_line("  0. standard double precision;");
    put_line("  1. double double precision; or");
    put_line("  2. quad double precision.");
    put("Type 0, 1, or 2 to select the precision : ");
    Ask_Alternative(ans,"012");
    case ans is
      when '0' => Standard_Main;
      when '1' => DoblDobl_Main;
      when '2' => QuadDobl_Main;
      when others => null;
    end case;
  end Main;

begin
  Main;
end ts_parcon; 
