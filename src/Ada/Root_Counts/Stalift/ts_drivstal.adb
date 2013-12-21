with text_io;                            use text_io;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Communications_with_User;           use Communications_with_User;
with Standard_Complex_Poly_Systems;      use Standard_Complex_Poly_Systems;
with Standard_Complex_Poly_Systems_io;   use Standard_Complex_Poly_Systems_io;
with Standard_Complex_Solutions;         use Standard_Complex_Solutions;
with Drivers_for_Static_Lifting;         use Drivers_for_Static_Lifting;
with Drivers_for_MixedVol_Algorithm;     use Drivers_for_MixedVol_Algorithm;

procedure ts_drivstal is

-- DESCRIPTION :
--   This procedure calls the driver to static lifting.

  procedure Main is

    lp,lq : Link_to_Poly_Sys;
    qsols,qsols0 : Solution_List;
    nt,mv,smv,tmv : natural32 := 0;
    file : file_type;
    ans : character;

  begin
    new_line;
    put_line("Test on driver for static lifting.");
    new_line;
    get(lp);
    new_line;
    put_line("Reading the name of the output file.");
    Read_Name_and_Create_File(file);
    lq := new Poly_Sys(lp'range);
    new_line;
    put("Use MixedVol algorithm ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans = 'y' then
      new_line;
      put("Use multiprecision arithmetic for Hermite normal form ? (y/n) ");
      Ask_Yes_or_No(ans);
      Driver_for_MixedVol_Algorithm
        (file,0,lp.all,true,lq.all,qsols,qsols0,mv,smv,tmv,(ans = 'y'));
    else
      new_line;
      put("Give the number of tasks (0 for no multitasking) : ");
      get(nt);
      Driver_for_Mixed_Volume_Computation
        (file,integer32(nt),lp.all,true,lq.all,qsols,qsols0,mv,smv,tmv);
    end if;
  end Main;

begin
  Main;
end ts_drivstal;
