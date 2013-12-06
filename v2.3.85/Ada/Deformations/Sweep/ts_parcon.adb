with text_io;                           use text_io;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Complex_Poly_Systems;     use Standard_Complex_Poly_Systems;
with Standard_Complex_to_Real_Poly;
with Standard_Complex_Solutions;        use Standard_Complex_Solutions;
with Systems_with_Parameters;           use Systems_with_Parameters;
with Parameter_Homotopy_Continuation;   use Parameter_Homotopy_Continuation;

procedure ts_parcon is

-- DESCRIPTION :
--   Stand-alone calling procedure to coefficient parameter continuation.

  procedure Main is

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
  end Main;

begin
  Main;
end ts_parcon; 
