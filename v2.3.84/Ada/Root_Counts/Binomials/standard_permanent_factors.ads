with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Vectors;
with Standard_Integer_Matrices;
with Standard_Complex_Laur_Systems;      use Standard_Complex_Laur_Systems;
with Standard_Monomial_Maps;             use Standard_Monomial_Maps;

package Standard_Permanent_Factors is

-- DESCRIPTION :
--   The enumeration of all affine solution sets is connected to the
--   enumeration of all factors that contributed to a permanent.
--   This package offers drivers to enumerate all selections of subsets
--   of variables to be set to zero and tools to process such selections.
--   The routines come in two flavors:
--   (1) without output data structures, just printing to screen;
--   (2) with monomial maps on output.

  procedure Solve_Affine_Subsystem 
               ( output : in boolean; p : in Laur_Sys;
                 eqcnt,s0cnt : in integer32;
                 eq,s0,free : in Standard_Integer_Vectors.Vector );

  -- DESCRIPTION :
  --   Solves the subsystem of p as defined by eq and shows the solutions.

  -- REQUIRED : eqcnt > 0.

  -- ON ENTRY :
  --   output    true if additional diagnostics during the solving
  --             procedure must be shown, false otherwise;
  --   p         a binomial system;
  --   eqcnt     number of equations to select from p;
  --   s0cnt     number of variables to be set to zero;
  --   eq        defines the equations of p that will be selected:
  --             eq(k) = 1 if the k-th equation of p will be selected,
  --             eq(k) = 0 if the k-th equation of p does not matter;
  --   s0        defines the variables that will be set to zero:
  --             s0(k) = 1 if the k-th variable will be set to zero,
  --             s0(k) = 0 otherwise;
  --   free      defines the free variables who no longer occur in the
  --             system after the selected variables as defined by s0
  --             are set to zero: free(k) = 1 indicates that the k-th
  --             variable is free, free(k) = 0 otherwise.

  function Monomial_Map_Solution
               ( output : in boolean; p : in Laur_Sys;
                 eqcnt,s0cnt : in integer32;
                 eq,s0,free : in Standard_Integer_Vectors.Vector )
             return Link_to_Monomial_Map_Array;

  -- DESCRIPTION :
  --   Selects the subsystem of p and return the solution map.

  -- REQUIRED : eqcnt > 0.

  -- ON ENTRY :
  --   output    true if the solutions should be shown on screen;
  --   p         a binomial system;
  --   eqcnt     number of equations to select from p;
  --   s0cnt     number of variables to be set to zero;
  --   eq        defines the equations of p that will be selected:
  --             eq(k) = 1 if the k-th equation of p will be selected,
  --             eq(k) = 0 if the k-th equation of p does not matter;
  --   s0        defines the variables that will be set to zero:
  --             s0(k) = 1 if the k-th variable will be set to zero,
  --             s0(k) = 0 otherwise;
  --   free      defines the free variables who no longer occur in the
  --             system after the selected variables as defined by s0
  --             are set to zero: free(k) = 1 indicates that the k-th
  --             variable is free, free(k) = 0 otherwise.

  -- ON RETURN :
  --   An array of maps (possibly empty if overconstrained) as solutions.

  procedure Show_Selection
               ( p : in Laur_Sys; A : in Standard_Integer_Matrices.Matrix;
                 s0,s1 : in Standard_Integer_Vectors.Vector;
                 cnt,s0cnt,eqcnt : in integer32; fail : out boolean );

  -- DESCRIPTION :
  --   Shows the selection of variables made and calls the solver.

  -- ON ENTRY :
  --   p         a binomial system;
  --   A         the incidence matrix of the binomial system p;
  --   s0        variables set to zero;
  --   s1        variables to be nonzero;
  --   cnt       index of the selection;
  --   s0cnt     number of variables in s0;
  --   eqcnt     number of remaining equations.

  -- ON RETURN :
  --   fail      true if the remaining equations after setting variables
  --             to zero as indicated by s0 are no longer binomials,
  --             false otherwise.

  procedure Append_Solution_Maps
               ( output : in boolean;
                 p : in Laur_Sys; A : in Standard_Integer_Matrices.Matrix;
                 s0,s1 : in Standard_Integer_Vectors.Vector;
                 cnt,s0cnt,eqcnt : in integer32; fail : out boolean;
                 first,last : in out Monomial_Map_List );

  -- DESCRIPTION :
  --   Concatenates the contributing factors to the permanent selection
  --   as monomial maps to the already accumulated affine solution sets.

  -- ON ENTRY :
  --   output    true if the solutions and other diagnostics are to
  --             be shown on screen, otherwise no output;
  --   p         a binomial system;
  --   A         the incidence matrix of the binomial system p;
  --   s0        variables set to zero;
  --   s1        variables to be nonzero;
  --   cnt       index of the selection;
  --   s0cnt     number of variables in s0;
  --   eqcnt     number of remaining equations;
  --   first     first element of the list of already computed maps;
  --   last      last computed solution map.

  -- ON RETURN :
  --   fail      true if the remaining equations after setting variables
  --             to zero as indicated by s0 are no longer binomials,
  --             false otherwise;
  --   first     updated list of solution maps;
  --   last      last computed solution map.

  procedure Selection_Iterator
               ( p : in Laur_Sys; A : in Standard_Integer_Matrices.Matrix;
                 max : in integer32; puredim : in boolean );
  procedure Selection_Iterator
               ( p : in Laur_Sys; A : in Standard_Integer_Matrices.Matrix;
                 max : in integer32; puredim,output : in boolean;
                 sols : out Monomial_Map_List );

  -- DESCRIPTION :
  --   Calls an iterator to enumerate all selections of variables
  --   to be set to zero for all candidate affine solution sets,
  --   for a binomial system with incidence matrix in A.

  procedure Recursive_Enumerator
               ( p : in Laur_Sys; A : in Standard_Integer_Matrices.Matrix;
                 max : in integer32; puredim : in boolean );
  procedure Recursive_Enumerator
               ( p : in Laur_Sys; A : in Standard_Integer_Matrices.Matrix;
                 max : in integer32; puredim,output : in boolean;
                 sols : out Monomial_Map_List );

  -- DESCRIPTION :
  --   Recursively enumerates all selections of variables
  --   to be set to zero for all candidate affine solution sets,
  --   for a binomial system with incidence matrix in A.

  procedure Silent_Affine_Solutions_with_Recursion
               ( p : in Laur_Sys; sols : out Monomial_Map_List;
                 fail : out boolean );
  procedure Interactive_Affine_Solutions_with_Recursion
               ( p : in Laur_Sys; sols : out Monomial_Map_List;
                 fail : out boolean );

  -- DESCRIPTION :
  --   Searches for affine solutions of a binomial system
  --   with a recursive enumeration procedure.
  --   If the system is not binomial, then fail is true on return,
  --   otherwise fail is false on return.
  --   This procedure is a driver to the recursive enumerator.
  --   The silent version does not prompt for settings and stays mute,
  --   while the interactive version communicates with the user.
 
  procedure Silent_Affine_Solutions_with_Iterator
               ( p : in Laur_Sys; puretopdim : in boolean;
                 sols : out Monomial_Map_List; fail : out boolean );
  procedure Interactive_Affine_Solutions_with_Iterator
               ( p : in Laur_Sys; sols : out Monomial_Map_List;
                 fail : out boolean );

  -- DESCRIPTION :
  --   Searches for affine solutions of a binomial system,
  --   with an interator.  If the system is not binomial,
  --   then fail is true on return, otherwise fail is false.
  --   This procedure is a driver to the selection iterator.
  --   The silent version does not prompt for settings and stays mute,
  --   while the interactive version communicates with the user.

end Standard_Permanent_Factors;
