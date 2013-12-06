with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Natural_Vectors;
with Sets_of_Unknowns;                   use Sets_of_Unknowns;
with Partitions_of_Sets_of_Unknowns;     use Partitions_of_Sets_of_Unknowns;
with Standard_Complex_Poly_Systems;      use Standard_Complex_Poly_Systems;

package Degree_Structure is

-- DESCRIPTION :
--   This package analyzes the structure of a given 
--   polynomial system, with respect to its degrees. 

-- CREATORS :

  procedure Create ( p : in Poly_Sys );

  -- DESCRIPTION : 
  --   For each equation of the polynomial system p, a heuristic is applied
  --   for constructing a good tuple of partitions.
  --   This degree structure will be contained by the package itself.

  procedure Put ( p : in Poly_Sys;
                  i,m : in natural32; z : in Partition );

  -- DESCRIPTION : 
  --   A partition for the set of unknowns of the ith polynomial 
  --   of the system p is added.

-- SELECTORS :

  function Empty return boolean;

  -- DESCRIPTION :
  --   is true initially and after the Clear procedure.

  function Get ( i : natural32 ) return natural32;

  -- DESCRIPTION :
  --   returns the number of sets for the partition of
  --   the i-th equation of a polynomial system p for
  --   which the Create operation must be performed.

  procedure Get ( i : in natural32; z : in out Partition;
                  d : out Standard_Natural_Vectors.Vector );

  -- DESCRIPTION :
  --   The partition for the i-the equation of the polynomial 
  --   system p is returned together with its degrees;
  --   for j in d'range :  d(j) = Degree(p(i),z(j)).

  -- REQUIRED :
  --   The Create operation or Put operations have been performed on a 
  --   polynomial system p having at least i equations. 

-- COMPUTING THE GENERALIZED BEZOUT NUMBER : 

  function Admissible ( z : Partition; n : natural32 ) return boolean;

  -- DESCRIPTION : 
  --   return true if the following holds: 
  --    `Any union of k sets of z contains at least k unknowns',
  --   for k in 2..n. 

  -- NOTE :
  --   z is not necessary a partition, it is rather a set of n sets 
  --   of unknowns.

  function Admissible ( z : Partition; n : natural32; s : Set )
                      return boolean;

  -- DESCRIPTION : 
  --   given admissible(z,n) = true, this function returns 
  --   admissible(union(z,s),n+1).

  function Admissible ( z : Partition; k,n : natural32; s : Set )
                      return boolean;

  -- DESCRIPTION : 
  --   returns true if admissible(union(z,s),n+1), with respect to k.

  function Generalized_Bezout_Number return natural32;

  -- DESCRIPTION : 
  --   This function returns a Bezout number based on the
  --   degree structure of the system.

  function Generalized_Bezout_Number ( p : in Poly_Sys ) return natural32;

  -- DESCRIPTION :
  --   Calling this function is equivalent to :
  --      Degree_Structure.Create(p);
  --      return Degree_Structure.Generalized_Bezout_Number; 

-- DESTRUCTOR :

  procedure Clear;

  -- DESCRIPTION :
  --   All allocated memory space is freed. 

end Degree_Structure;
