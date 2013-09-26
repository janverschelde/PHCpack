with Quad_Double_Ring;
with Quad_Double_Numbers;                 use Quad_Double_Numbers;
with Abstract_Ring.Field;

package Quad_Double_Ring.FField is
  new Quad_Double_Ring.Field("<",">",AbsVal,"/",Div);

-- DESCRIPTION :
--   Defines the extension of the ring of quad double numbers
--   to a field.
