with Standard_Floating_Vectors;
with Standard_Complex_Vectors;
with Standard_Integer_Matrices;
with Multprec_Integer_Matrices;
with Standard_Complex_Solutions;       use Standard_Complex_Solutions;

package Standard_Radial_Solvers is

-- DESCRIPTION :
--   For numerical stability, it is best to separate the radii from
--   the arguments in the righthand side of a binomial system when
--   solving it.  The routines in this package help to compute the 
--   magnitude of the solution of a binomial system.

  function Radii ( c : Standard_Complex_Vectors.Vector )
                 return Standard_Floating_Vectors.Vector;

  -- DESCRIPTION :
  --   Returns the vector of radii for all complex numbers in c.

  function Scale ( c : Standard_Complex_Vectors.Vector;
                   r : Standard_Floating_Vectors.Vector )
                 return Standard_Complex_Vectors.Vector;

  -- DESCRIPTION :
  --   Return c/r, divisions executed componentwise.

  function Log10 ( r : Standard_Floating_Vectors.Vector ) 
                 return Standard_Floating_Vectors.Vector;

  -- DESCRIPTION :
  --   Returns the logarithm of all the real numbers in r.

  function Exp10 ( r : Standard_Floating_Vectors.Vector ) 
                 return Standard_Floating_Vectors.Vector;

  -- DESCRIPTION :
  --   Returns the vector 10^r.

  function Radial_Upper_Solve
              ( U : Standard_Integer_Matrices.Matrix;
                logr : Standard_Floating_Vectors.Vector )
              return Standard_Floating_Vectors.Vector;
  function Radial_Upper_Solve
              ( U : Multprec_Integer_Matrices.Matrix;
                logr : Standard_Floating_Vectors.Vector )
              return Standard_Floating_Vectors.Vector;

  -- DESCRIPTION :
  --   Solves U*log(x) = logr, for U an upper triangular exponent matrix.

  function Multiply ( A : in Standard_Integer_Matrices.Matrix;
                      x : in Standard_Floating_Vectors.Vector )
                    return Standard_Floating_Vectors.Vector;
  function Multiply ( A : in Multprec_Integer_Matrices.Matrix;
                      x : in Standard_Floating_Vectors.Vector )
                    return Standard_Floating_Vectors.Vector;

  -- DESCRIPTION : returns A*x, taking the transpose of A.

  function Eval ( A : in Standard_Integer_Matrices.Matrix;
                  x : in Standard_Floating_Vectors.Vector )
                return Standard_Floating_Vectors.Vector;
  function Eval ( A : in Multprec_Integer_Matrices.Matrix;
                  x : in Standard_Floating_Vectors.Vector )
                return Standard_Floating_Vectors.Vector;

  -- DESCRIPTION : returns x^A.

  procedure Multiply ( s : in out Standard_Complex_Vectors.Vector;
                       r : in Standard_Floating_Vectors.Vector );
  procedure Multiply ( s : in out Solution;
                       r : in Standard_Floating_Vectors.Vector );
  procedure Multiply ( s : in out Solution_List;
                       r : in Standard_Floating_Vectors.Vector );

  -- DESCRIPTION :
  --   Multiplies every solution vector s componentwise by r.

end Standard_Radial_Solvers;
