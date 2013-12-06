with text_io;                             use text_io;
with Standard_Natural_Numbers;            use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;         use Standard_Natural_Numbers_io;
with Multprec_Floating_Numbers;           use Multprec_Floating_Numbers;
with Multprec_Floating_Numbers_io;        use Multprec_Floating_Numbers_io;
with Multprec_Floating_Constants;

procedure ts_mpcst is

-- DESCRIPTION :
--   Development of multiprecision constants.

  procedure Main is

    d : natural32 := 0;
    f,g : Floating_Number;

  begin
    put("Give number of decimal places : "); get(d);
    f := Multprec_Floating_Constants.Pi(d);
    g := Multprec_Floating_Constants.TwoPi(d);
    put("Pi with "); put(d,1); put(" decimal places : "); put(f); new_line;
    put("2*Pi with "); put(d,1); put(" decimal places : "); put(g);
  end Main;

begin
  Main;
end ts_mpcst;
