with Triple_Double_Ring;
with Triple_Double_Numbers;               use Triple_Double_Numbers;
with Abstract_Ring.Field;

package Triple_Double_Ring.FField is
  new Triple_Double_Ring.Field("<",">",AbsVal,"/",Div);

-- DESCRIPTION :
--   Extends the ring of triple double numbers to a field.
