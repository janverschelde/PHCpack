with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with DoblDobl_Complex_Poly_Systems;
with DoblDobl_Complex_Laur_Systems;
with DoblDobl_Laur_Poly_Convertors;
with DoblDobl_Complex_Solutions;
with DoblDobl_Complex_Solutions_io;      use DoblDobl_Complex_Solutions_io;
with DoblDobl_Root_Refiners;             use DoblDobl_Root_Refiners;
with Drivers_to_DD_QD_Root_Refiners;     use Drivers_to_DD_QD_Root_Refiners;

procedure ts_ddnewt is

-- DESCRIPTION :
--  Test on Newton's method in double double arithmetic.

  procedure Main is

    ans : character;
    p : DoblDobl_Complex_Laur_Systems.Link_to_Laur_Sys;
    s : DoblDobl_Complex_Solutions.Solution_List;

  begin
    new_line;
    put_line("MENU for Newton's method in double double arithmetic.");
    put_line("  1. start from system with standard complex numbers;");
    put_line("  2. give system and solution in multiprecision numbers.");
    put("Type 1 or 2 to make a choice : ");
    Ask_Alternative(ans,"12");
    if ans = '1'
     then Standard_to_DoblDobl_Complex(p,s);
     else Multprec_to_DoblDobl_Complex(p,s);
    end if;
    put("read "); 
    put(DoblDobl_Complex_Solutions.Length_Of(s),1);
    put_line(" solutions");
    if DoblDobl_Laur_Poly_Convertors.Is_Genuine_Laurent(p.all) then
      DoblDobl_Root_Refiner(p.all,s);
    else
      declare
        q : DoblDobl_Complex_Poly_Systems.Poly_Sys(p'range)
          := DoblDobl_Laur_Poly_Convertors.Laurent_to_Polynomial_System(p.all);
      begin
        DoblDobl_Root_Refiner(q,s);
      end;
    end if;
    put_line("The refined solutions :");
    put(standard_output,DoblDobl_Complex_Solutions.Length_Of(s),
        natural32(DoblDobl_Complex_Solutions.Head_Of(s).n),s);
  end Main;

begin
  Main;
end ts_ddnewt;
