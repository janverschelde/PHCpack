with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Integer_Vectors;
with Standard_Integer_Matrices;
with Standard_Complex_Laur_Systems;     use Standard_Complex_Laur_Systems;
with Standard_Complex_Solutions;        use Standard_Complex_Solutions;
with Standard_Monomial_Maps;            use Standard_Monomial_Maps;

package Standard_Monomial_Map_Solvers is

-- DESCRIPTION :
--   A standard monomial map solver is a solver of a binomial system
--   that returns an array or a list of monomial maps.

  function Free_Monomial_Map
              ( n,d : integer32;
                free : Standard_Integer_Vectors.Vector )
              return Monomial_Map;

  -- DESCRIPTION :
  --   Returns a d-dimensional monomial map of d free variables.

  -- REQUIRED : d equals the number of ones in free.

  -- ON ENTRY :
  --   n        the number of variables;
  --   d        the number of free variables;
  --   free     if the k-th variable is free, then f(k) = 1.

  function Monomial_Map_Solution
              ( n,d : integer32;
                M : Standard_Integer_Matrices.Matrix; c : Solution )
              return Monomial_Map;
  function Monomial_Map_Solution
              ( n,d : integer32;
                s : Standard_Integer_Vectors.Vector;
                M : Standard_Integer_Matrices.Matrix; c : Solution )
              return Monomial_Map;
  function Monomial_Map_Solution
              ( n,d,e : integer32;
                s,f : Standard_Integer_Vectors.Vector;
                M : Standard_Integer_Matrices.Matrix; c : Solution )
              return Monomial_Map;

  -- DESCRIPTION :
  --   Returns the map in n variables and d+e parameters defined by
  --   the tropisms in the first d rows of the matrix M and the solutions
  --   to the reduced binomial system in c, taking into account the
  --   variables set to zero, as defined in s, and the free variables
  --   as defined in f.  The number of free variables equals e.
  --   If not provided, then e is consider zero.

  function Monomial_Map_Solutions 
              ( n,d : integer32;
                M : Standard_Integer_Matrices.Matrix; c : Solution_List )
              return Monomial_Map_Array;
  function Monomial_Map_Solutions 
              ( n,d : integer32;
                s : in Standard_Integer_Vectors.Vector;
                M : Standard_Integer_Matrices.Matrix; c : Solution_List )
              return Monomial_Map_Array;
  function Monomial_Map_Solutions 
              ( n,d,e : integer32;
                s,f : in Standard_Integer_Vectors.Vector;
                M : Standard_Integer_Matrices.Matrix; c : Solution_List )
              return Monomial_Map_Array;

  -- DESCRIPTION :
  --   Returns an array of monomial maps to represent the solutions
  --   with coefficients in c and tropisms in the first d rows of M.

  -- REQUIRED : Length_Of(c) > 0.

  function Toric_Solve ( p : Laur_Sys ) return Link_to_Monomial_Map_Array;

  -- DESCRIPTION :
  --   Solves the binomial system defined by p
  --   and returns the solutions as monomial maps.
  --   and stores the solutions as monomial maps.

  function Affine_Solve
              ( p : Laur_Sys;
                s : Standard_Integer_Vectors.Vector )
              return Link_to_Monomial_Map_Array;

  -- DESCRIPTION :
  --   Solves the binomial system defined by p
  --   and returns the solution maps, taking into account the
  --   selected variables set to zero, as defined by s.

  function Affine_Solve
              ( p : Laur_Sys; nbf : integer32;
                s,f : Standard_Integer_Vectors.Vector )
              return Link_to_Monomial_Map_Array;

  -- DESCRIPTION :
  --   Solves the binomial system defined by p
  --   and returns the solution maps, taking into account the
  --   selected variables set to zero, as defined by s,
  --   and the free variables as defined by f.
  --   The number of free variables equals nbf.
  
end Standard_Monomial_Map_Solvers;
