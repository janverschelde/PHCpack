with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with Test_Multprec_DoblDobl_Casts;
with Test_Multprec_TripDobl_Casts;
with Test_Multprec_QuadDobl_Casts;
with Test_Multprec_PentDobl_Casts;
with Test_Multprec_OctoDobl_Casts;
with Test_Multprec_DecaDobl_Casts;
with Test_Multprec_HexaDobl_Casts;

procedure ts_mpcasts is

-- DESCRIPTION :
--   Calls the main interactive test on type casts between
--   multiprecision and multiple double precision numbers.

  ans : character;

begin
  new_line;
  put_line("MENU for multiprecision to multiple double type casts :");
  put_line("  2. double double precision");
  put_line("  3. triple double precision");
  put_line("  4. quad double precision");
  put_line("  5. penta double precision");
  put_line("  6. octo double precision");
  put_line("  7. deca double precision");
  put_line("  8. deca double precision");
  put("Type 2, 3, 4, 5, 6, 7, or 8 to select a precistion : ");
  Ask_Alternative(ans,"23458A");
  case ans is
    when '2' => Test_Multprec_DoblDobl_Casts.Main;
    when '3' => Test_Multprec_TripDobl_Casts.Main;
    when '4' => Test_Multprec_QuadDobl_Casts.Main;
    when '5' => Test_Multprec_PentDobl_Casts.Main;
    when '6' => Test_Multprec_OctoDobl_Casts.Main;
    when '7' => Test_Multprec_DecaDobl_Casts.Main;
    when '8' => Test_Multprec_HexaDobl_Casts.Main;
    when others => null;
  end case;
end ts_mpcasts;
