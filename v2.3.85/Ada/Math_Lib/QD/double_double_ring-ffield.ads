with Double_Double_Ring;
with Double_Double_Numbers;               use Double_Double_Numbers;
with Abstract_Ring.Field;

package Double_Double_Ring.FField is
  new Double_Double_Ring.Field("<",">",AbsVal,"/",Div);

-- DESCRIPTION :
--   Defines the extension of the ring of double double numbers
--   to a field.
