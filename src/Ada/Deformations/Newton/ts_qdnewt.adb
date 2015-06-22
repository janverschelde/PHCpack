with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with QuadDobl_Complex_Poly_Systems;
with QuadDobl_Complex_Laur_Systems;
with QuadDobl_Laur_Poly_Convertors;
with QuadDobl_Complex_Solutions;
with QuadDobl_Complex_Solutions_io;      use QuadDobl_Complex_Solutions_io;
with QuadDobl_Root_Refiners;             use QuadDobl_Root_Refiners;
with Drivers_to_DD_QD_Root_Refiners;     use Drivers_to_DD_QD_Root_Refiners;

procedure ts_qdnewt is

-- DESCRIPTION :
--   Development of Newton's method in quad double arithmetic.

  procedure Main is

    ans : character;
    p : QuadDobl_Complex_Laur_Systems.Link_to_Laur_Sys;
    s : QuadDobl_Complex_Solutions.Solution_List;

  begin
    new_line;
    put_line("MENU for Newton's method in quad double arithmetic.");
    put_line("  1. start from system with standard complex numbers;");
    put_line("  2. give system and solution in multiprecision numbers.");
    put("Type 1 or 2 to make a choice : ");
    Ask_Alternative(ans,"12");
    if ans = '1'
     then Standard_to_QuadDobl_Complex(p,s);
     else Multprec_to_QuadDobl_Complex(p,s);
    end if;
    put("read "); 
    put(QuadDobl_Complex_Solutions.Length_Of(s),1);
    put_line(" solutions");
    if QuadDobl_Laur_Poly_Convertors.Is_Genuine_Laurent(p.all) then
      QuadDobl_Root_Refiner(p.all,s);
    else
      declare
        q : QuadDobl_Complex_Poly_Systems.Poly_Sys(p'range)
          := QuadDobl_Laur_Poly_Convertors.Laurent_to_Polynomial_System(p.all);
      begin
        QuadDobl_Root_Refiner(q,s);
      end;
    end if;
    put_line("The solutions after refinement :");
    Write(standard_output,s);
  end Main;

begin
  Main;
end ts_qdnewt;
