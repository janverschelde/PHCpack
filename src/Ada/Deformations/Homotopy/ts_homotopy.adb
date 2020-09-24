with text_io;                           use text_io;
with Communications_with_User;          use Communications_with_User;
with Test_Standard_Laur_Homotopy;
with Test_Standard_Poly_Homotopy;
with Test_DoblDobl_Poly_Homotopy;
with Test_TripDobl_Poly_Homotopy;
with Test_QuadDobl_Poly_Homotopy;
with Test_PentDobl_Poly_Homotopy;
with Test_OctoDobl_Poly_Homotopy;
with Test_DecaDobl_Poly_Homotopy;

procedure ts_homotopy is

-- DESCRIPTION :
--   Tests polynomial homotopies.

  procedure Main is

  -- DESCRIPTION :
  --   Prompts for the choice between Laurent and polynomial homotopies.
  --   If not Laurent, then prompts for the precision.

    ans : character;

  begin
    new_line;
    put_line("Testing polynomial homotopies ...");
    new_line;
    put("Laurent polynomial systems ? (y/n) "); Ask_Yes_or_No(ans);
    if ans = 'y' then
      Test_Standard_Laur_Homotopy.Main;
    else
      new_line;
      put_line("MENU for the working precision :");
      put_line("  1. double precision");
      put_line("  2. double double precision");
      put_line("  3. triple double precision");
      put_line("  4. quad double precision");
      put_line("  5. penta double precision");
      put_line("  6. octo double precision");
      put_line("  7. deca double precision");
      put("Type 1, 2, 3, 4, 5, 6, or 7 to select the precision : ");
      Ask_Alternative(ans,"1234567");
      case ans is
        when '1' => Test_Standard_Poly_Homotopy.Main;
        when '2' => Test_DoblDobl_Poly_Homotopy.Main;
        when '3' => Test_TripDobl_Poly_Homotopy.Main;
        when '4' => Test_QuadDobl_Poly_Homotopy.Main;
        when '5' => Test_PentDobl_Poly_Homotopy.Main;
        when '6' => Test_OctoDobl_Poly_Homotopy.Main;
        when '7' => Test_DecaDobl_Poly_Homotopy.Main;
        when others => null;
      end case;
    end if;
  end Main;

begin
  Main;
end ts_homotopy;
