with Penta_Double_Ring;
with Penta_Double_Numbers;                use Penta_Double_Numbers;
with Abstract_Ring.Field;

package Penta_Double_Ring.FField is
  new Penta_Double_Ring.Field("<",">",AbsVal,"/",Div);

-- DESCRIPTION :
--   Extends the ring of penta double numbers to a field.
