with Multprec_Integer64_Ring;
with Multprec_Integer64_Numbers;         use Multprec_Integer64_Numbers;
with Abstract_Ring.Domain;

package Multprec_Integer64_Ring.DDomain is
  new Multprec_Integer64_Ring.Domain("<",">","/",Rmd,Rmd,Div,Div,Div);

-- DESCRIPTION :
--   Defines the extension of the ring of multi-precision integer numbers
--   to an Euclidean domain.
