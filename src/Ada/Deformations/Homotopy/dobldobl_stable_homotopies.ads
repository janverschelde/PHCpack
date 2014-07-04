with Standard_Integer_Numbers;            use Standard_Integer_Numbers;
with Standard_Integer_Vectors;
with DoblDobl_Complex_Vectors;
with DoblDobl_Complex_Polynomials;
with DoblDobl_Complex_Laurentials;
with DoblDobl_Complex_Laur_Systems;       use DoblDobl_Complex_Laur_Systems;
with DoblDobl_Complex_Poly_Systems;       use DoblDobl_Complex_Poly_Systems;
with DoblDobl_Complex_Solutions;          use DoblDobl_Complex_Solutions;

package DoblDobl_Stable_Homotopies is

-- DESCRIPTION :
--   This package provides operations to create homotopies that
--   allow the stable computation of solutions with zero components,
--   as predicted by stable mixed volumes.

-- TYPE OF ZERO :
--   is described by an integer vector z and number of zeroes nbz.
--   For every entry i of z, the signs means the following:
--     z(i) < 0 : the i-th component will go to infinity (spurious),
--     z(i) = 0 : the i-th component will go to zero,
--     z(i) > 0 : the i-th component will most likely be nonzero.
--   In case all entries of z are nonnegative, the nbz counts the
--   number of zero components in the solution, otherwise nbz = -1
--   for a spurious solution.

  function Number_of_Zeroes
             ( z : Standard_Integer_Vectors.Vector ) return integer32;

  -- DESCRIPTION :
  --   For z on input the zero type of a mixed cell,
  --   returns -1 if the zero type is for a spurious cell,
  --            0 if all components will be nonzero, otherwise
  --   returns the number of zero components in z.

  procedure Zero_Type
              ( v : in DoblDobl_Complex_Vectors.Vector; nz : out integer32;
                z : out Standard_Integer_Vectors.Vector );

  -- DESCRIPTION :
  --   Returns in z the type of zero components of the vector v.
  --   The number of zero components equals nz on return.

  function Remove_Zeroes
             ( s : in Solution; nz : in integer32;
               z : in Standard_Integer_Vectors.Vector ) return Solution;

  -- DESCRIPTION :
  --   Removes nz zeroes from the solution s,
  --   according to the zero type in z.

  -- REQUIRED : nz < s.n, otherwise the solution is entirely zero.

  function Origin ( n,m : integer32 ) return Solution;

  -- DESCRIPTION :
  --   Returns the origin as a solution of dimension n and multiplicity m.
  --   The full case of this zero type involves no continuation.

  function Insert_Zeroes
              ( v : DoblDobl_Complex_Vectors.Vector;
                z : Standard_Integer_Vectors.Vector )
              return DoblDobl_Complex_Vectors.Vector;

  -- DESCRIPTION :
  --   Returns a solution vector of z'range with the zeroes at those
  --   corresponding entries in the zero type z.  The nonzero values
  --   are taken from the vector v.

  -- REQUIRED : v'last >= number of nonzero values in z.

  function Insert_Zeroes
              ( s : Solution; z : Standard_Integer_Vectors.Vector )
              return Solution;
  function Insert_Zeroes
              ( s : Solution_List; z : Standard_Integer_Vectors.Vector )
              return Solution_List;

  -- DESCRIPTION :
  --   The solution(s) on return has(have) in its(their) vector(s) the 
  --   appropriate zeroes inserted, according to the zero type in z.

  function Is_Same ( s1,s2 : Link_to_Solution ) return boolean;

  -- DESCRIPTION :
  --   Returns true if both solutions have the same zero type
  --   and are the same otherwise.

  procedure Merge_and_Concat
              ( first,last : in out Solution_List; sols : in Solution_List );

  -- DESCRIPTION :
  --   Merges the solutions sols into the list headed by first
  --   and whose last element is pointed to by last.
  --   For solutions with the same zero type, the multiplicity 
  --   of the corresponding solution in first is increased
  --   by the multiplicity of the corresponding solution in sols.

  function Vanish_by_Zeroes
             ( t : DoblDobl_Complex_Polynomials.Term;
               z : Standard_Integer_Vectors.Vector;
               nbz : integer32 ) return boolean;
  function Vanish_by_Zeroes
             ( t : DoblDobl_Complex_Laurentials.Term;
               z : Standard_Integer_Vectors.Vector;
               nbz : integer32 ) return boolean;

  -- DESCRIPTION :
  --   If nbz <= 0, then false is returned, otherwise, for nbz > 0,
  --   true is returned if the term t vanishes along the zero type z.
  --   A term vanishes if it has a nonzero exponent for a corresponding
  --   zero entry in z.

  function Substitute_Zeroes
             ( t : DoblDobl_Complex_Polynomials.Term;
               z : Standard_Integer_Vectors.Vector;
               nbz : integer32 ) return DoblDobl_Complex_Polynomials.Term;
  function Substitute_Zeroes
             ( t : DoblDobl_Complex_Laurentials.Term;
               z : Standard_Integer_Vectors.Vector;
               nbz : integer32 ) return DoblDobl_Complex_Laurentials.Term;

  -- DESCRIPTION :
  --   If nbz <= 0, t is returned, otherwise if t vanishes along z,
  --   the coefficient of the term on return will equal zero,
  --   if t does not vanish along z, it will have nbz less variables.

  function Substitute_Zeroes
             ( p : DoblDobl_Complex_Polynomials.Poly;
               z : Standard_Integer_Vectors.Vector;
               nbz : integer32 ) return DoblDobl_Complex_Polynomials.Poly;
  function Substitute_Zeroes
             ( p : DoblDobl_Complex_Laurentials.Poly;
               z : Standard_Integer_Vectors.Vector;
               nbz : integer32 ) return DoblDobl_Complex_Laurentials.Poly;

  -- DESCRIPTION :
  --   Applies Substitute_Zeroes to all terms in p.

  function Substitute_Zeroes
             ( p : Poly_Sys; z : Standard_Integer_Vectors.Vector;
               nbz : integer32 ) return Poly_Sys;
  function Substitute_Zeroes
             ( p : Laur_Sys; z : Standard_Integer_Vectors.Vector;
               nbz : integer32 ) return Laur_Sys;

  -- DESCRIPTION :
  --   Substitutes the zeroes in the polynomials in p.

  function Filter ( p : Poly_Sys ) return Poly_Sys;
  function Filter ( p : Laur_Sys ) return Laur_Sys;

  -- DESCRIPTION :
  --   Removes from p all polynomials with degree <= 0.

end DoblDobl_Stable_Homotopies;
