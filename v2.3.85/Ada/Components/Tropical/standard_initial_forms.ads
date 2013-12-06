with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Vectors;           use Standard_Integer_Vectors;
with Standard_Integer_Matrices;          use Standard_Integer_Matrices;
with Standard_Complex_Laurentials;
with Standard_Complex_Polynomials;
with Standard_Complex_Poly_Systems;      use Standard_Complex_Poly_Systems;
with Standard_Complex_Laur_Systems;      use Standard_Complex_Laur_Systems;

package Standard_Initial_Forms is

-- DESCRIPTION :
--   The initial form of a (Laurent) polynomial with respect to a direction v
--   is supported on the face of the Newton polytope determined by v.

  function Degree ( t : Standard_Complex_Polynomials.Term; v : Vector )
                  return integer32;
  function Degree ( t : Standard_Complex_Laurentials.Term; v : Vector )
                  return integer32;

  -- DESCRIPTION :
  --   Returns the inner product of the exponents in t with the vector v.

  function Degree ( p : Standard_Complex_Polynomials.Poly; v : Vector )
                  return integer32;
  function Degree ( p : Standard_Complex_Laurentials.Poly; v : Vector )
                  return integer32;

  -- DESCRIPTION :
  --   Returns the minimal inner product of the exponent vectors in p
  --   with the vector v.

  function Form ( p : Standard_Complex_Polynomials.Poly;
                  v : Vector; d : integer32 )
                return Standard_Complex_Polynomials.Poly;
  function Form ( p : Standard_Complex_Laurentials.Poly;
                  v : Vector; d : integer32 )
                return Standard_Complex_Laurentials.Poly;

  -- DESCRIPTION :
  --   Returns all terms in p whose degree in v equals d.
  --   If d = degree(p,v), then the polynomial on return
  --   is the initial form of p.

  function Initial ( p : Standard_Complex_Polynomials.Poly; v : Vector )
                   return Standard_Complex_Polynomials.Poly;
  function Initial ( p : Standard_Complex_Laurentials.Poly; v : Vector )
                   return Standard_Complex_Laurentials.Poly;

  -- DESCRIPTION :
  --   Returns the initial form of p in the direction v.

  function Initial ( p : Poly_Sys; v : Vector ) return Poly_Sys;
  function Initial ( s : Laur_Sys; v : Vector ) return Laur_Sys;

  -- DESCRIPTION :
  --   Returns the initial form system of p in the direction of v.

  function Transform ( t : Standard_Complex_Laurentials.Term; m : Matrix )
                     return Standard_Complex_Laurentials.Term;
  function Transform ( p : Standard_Complex_Laurentials.Poly; m : Matrix )
                     return Standard_Complex_Laurentials.Poly;
  function Transform ( s : Laur_Sys; m : Matrix ) return Laur_Sys;

  -- DESCRIPTION :
  --   Returns the transformed term, polynomial, or system,
  --   after multiplying the exponent vectors by the unimodular m.

  function Eliminate ( t : Standard_Complex_Laurentials.Term; k : integer32 )
                     return Standard_Complex_Laurentials.Term;
  function Eliminate ( p : Standard_Complex_Laurentials.Poly; k : integer32 )
                     return Standard_Complex_Laurentials.Poly;
  function Eliminate ( s : Laur_Sys; k : integer32 ) return Laur_Sys;

  -- DESCRIPTION :
  --   Removes the k-th variable from every term, where k is the pivot
  --   in the direction v that defines the unimodular transformation m.

end Standard_Initial_Forms;
