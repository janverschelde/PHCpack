with Standard_Complex_Poly_Systems;      use Standard_Complex_Poly_Systems;

package Reduction_of_Nonsquare_Systems is

-- DESCRIPTION :
--   This package provides two different approaches for reducing
--   a overconstrained system to a square system.

  function Random_Square ( p : in Poly_Sys ) return Poly_Sys;

  -- DESCRIPTION :
  --   The system q (q'range = 1..N) on return will be square, like
  --     q(i) = p(i) + l(i,N+1)*p(N+1) + .. + l(i,n)*p(n), for i=1,2,..,N.
  --   This type of reduction is recommended when the solution of any
  --   subsystem of p (p'range = 1,..,n) has a connected solution component.

  function Reduced_Square ( p : in Poly_Sys ) return Poly_Sys;

  -- DESCRIPTION :
  --   The system q (q'range = 1,..,N) on return will be square, like
  --     q(i) = Rpoly(..(Rpoly(p(i),p(N+1)), p(n)), for i=1,2,..,N,
  --   i.e.: the additional polynomial will be used for reducing the system.
  --   This type of reduction is recommended when the solution set of the
  --   first N equation of p consists of isolated points that do not all
  --   satisfy the remaining equations.

end Reduction_of_Nonsquare_Systems;
