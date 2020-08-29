with Deca_Double_Ring;
with Deca_Double_Numbers;                 use Deca_Double_Numbers;
with Abstract_Ring.Field;

package Deca_Double_Ring.FField is
  new Deca_Double_Ring.Field("<",">",AbsVal,"/",Div);

-- DESCRIPTION :
--   Extends the ring of deca double numbers to a field.
