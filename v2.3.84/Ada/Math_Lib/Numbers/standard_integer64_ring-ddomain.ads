with Standard_Integer64_Ring;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Abstract_Ring.Domain;

package Standard_Integer64_Ring.DDomain is
  new Standard_Integer64_Ring.Domain("<",">","/",Rmd,Rmd,Div,Div,Div);

-- DESCRIPTION :
--   Defines the extension of the ring of standard integer numbers
--   of type long_long_integer to an Euclidean domain.
