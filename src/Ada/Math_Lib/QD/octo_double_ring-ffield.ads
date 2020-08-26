with Octo_Double_Ring;
with Octo_Double_Numbers;                 use Octo_Double_Numbers;
with Abstract_Ring.Field;

package Octo_Double_Ring.FField is
  new Octo_Double_Ring.Field("<",">",AbsVal,"/",Div);

-- DESCRIPTION :
--   Extends the ring of octo double numbers to a field.
