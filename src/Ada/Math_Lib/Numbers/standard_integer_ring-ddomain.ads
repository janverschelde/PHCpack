with Standard_Integer_Ring;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Abstract_Ring.Domain;

package Standard_Integer_Ring.DDomain is
  new Standard_Integer_Ring.Domain("<",">","/",Rmd,Rmd,Div,Div,Div);

-- DESCRIPTION :
--   Defines the extension of the ring of standard integer numbers
--   to an Euclidean domain.
