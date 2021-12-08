with Hexa_Double_Ring;
with Hexa_Double_Numbers;                 use Hexa_Double_Numbers;
with Abstract_Ring.Field;

package Hexa_Double_Ring.FField is
  new Hexa_Double_Ring.Field("<",">",AbsVal,"/",Div);

-- DESCRIPTION :
--   Extends the ring of hexa double numbers to a field.
