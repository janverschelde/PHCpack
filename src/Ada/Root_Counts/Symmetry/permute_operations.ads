with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Natural_Vectors;
with Standard_Integer_Vectors;
with Standard_Floating_Vectors;
with Standard_Complex_Vectors;
with DoblDobl_Complex_Vectors;
with QuadDobl_Complex_Vectors;
with Standard_Complex_Polynomials;
with Standard_Complex_Laurentials;
with Standard_Complex_Poly_Systems;
with Standard_Complex_Laur_Systems;
with DoblDobl_Complex_Laurentials;
with DoblDobl_Complex_Polynomials;
with DoblDobl_Complex_Poly_Systems;
with DoblDobl_Complex_Laur_Systems;
with QuadDobl_Complex_Laurentials;
with QuadDobl_Complex_Polynomials;
with QuadDobl_Complex_Poly_Systems;
with QuadDobl_Complex_Laur_Systems;
with Permutations;                       use Permutations;

package Permute_Operations is

-- DESCRIPTION :
--   This package provides permute operations on vectors,
--   on polynomials and on systems of polynomials.

  function "*" ( p : Permutation; v : Standard_Natural_Vectors.Vector )
	       return Standard_Natural_Vectors.Vector;

  function "*" ( p : Permutation; v : Standard_Integer_Vectors.Vector )
	       return Standard_Integer_Vectors.Vector;

  function "*" ( p : Permutation; v : Standard_Floating_Vectors.Vector )
	       return Standard_Floating_Vectors.Vector;

  function "*" ( p : Permutation; v : Standard_Complex_Vectors.Vector )
	       return Standard_Complex_Vectors.Vector;

  function "*" ( p : Permutation; v : DoblDobl_Complex_Vectors.Vector )
	       return DoblDobl_Complex_Vectors.Vector;

  function "*" ( p : Permutation; v : QuadDobl_Complex_Vectors.Vector )
	       return QuadDobl_Complex_Vectors.Vector;

  -- DESCRIPTION :
  --   Returns the result of the permutation of p on the vector v.

  -- REQUIRED :
  --   p'range = v'range

  function Permutable ( v1,v2 : Standard_Natural_Vectors.Vector )
                      return boolean;
  function Permutable ( v1,v2 : Standard_Integer_Vectors.Vector )
                      return boolean;
  function Permutable ( v1,v2 : Standard_Floating_Vectors.Vector )
                      return boolean;
  function Permutable ( v1,v2 : Standard_Complex_Vectors.Vector )
                      return boolean;
  function Permutable ( v1,v2 : Standard_Floating_Vectors.Vector;
                        tol : double_float ) return boolean;
  function Permutable ( v1,v2 : Standard_Complex_Vectors.Vector;
                        tol : double_float ) return boolean;

  -- DESCRIPTION :
  --   Returns true if there exists a permutation between the two vectors.
  --   If provided, tol is the tolerance for comparing two numeric values.

  function Sign_Permutable ( v1,v2 : Standard_Natural_Vectors.Vector )
                           return boolean;
  function Sign_Permutable ( v1,v2 : Standard_Integer_Vectors.Vector )
                           return boolean;
  function Sign_Permutable ( v1,v2 : Standard_Floating_Vectors.Vector )
                           return boolean;
  function Sign_Permutable ( v1,v2 : Standard_Complex_Vectors.Vector )
                           return boolean;
  function Sign_Permutable ( v1,v2 : Standard_Floating_Vectors.Vector;
                             tol : double_float ) return boolean;
  function Sign_Permutable ( v1,v2 : Standard_Complex_Vectors.Vector; 
                             tol : double_float ) return boolean;

  -- DESCRIPTION :
  --   Also permutations where the sign of one of the components can
  --   be changed, are checked.

  function "*" ( p : Permutation; t : Standard_Complex_Polynomials.Term )
	       return Standard_Complex_Polynomials.Term;
  function "*" ( p : Permutation; t : DoblDobl_Complex_Polynomials.Term )
	       return DoblDobl_Complex_Polynomials.Term;
  function "*" ( p : Permutation; t : QuadDobl_Complex_Polynomials.Term )
	       return QuadDobl_Complex_Polynomials.Term;

  -- DESCRIPTION :
  --   Permutes the unknowns in the term t according to the permutation p,
  --   in standard double, double double or quad double precision.

  function "*" ( p : Permutation; s : Standard_Complex_Polynomials.Poly )
	       return Standard_Complex_Polynomials.Poly;
  function "*" ( p : Permutation; s : DoblDobl_Complex_Polynomials.Poly )
	       return DoblDobl_Complex_Polynomials.Poly;
  function "*" ( p : Permutation; s : QuadDobl_Complex_Polynomials.Poly )
	       return QuadDobl_Complex_Polynomials.Poly;

  -- DESCRIPTION :
  --   Permutes the unknowns in s according to the permutation p,
  --   in standard double, double double or quad double precision.

  function "*" ( p : Permutation; t : Standard_Complex_Laurentials.Term )
	       return Standard_Complex_Laurentials.Term;
  function "*" ( p : Permutation; t : DoblDobl_Complex_Laurentials.Term )
	       return DoblDobl_Complex_Laurentials.Term;
  function "*" ( p : Permutation; t : QuadDobl_Complex_Laurentials.Term )
	       return QuadDobl_Complex_Laurentials.Term;
  function "*" ( p : Permutation; s : Standard_Complex_Laurentials.Poly )
               return Standard_Complex_Laurentials.Poly;
  function "*" ( p : Permutation; s : DoblDobl_Complex_Laurentials.Poly )
               return DoblDobl_Complex_Laurentials.Poly;
  function "*" ( p : Permutation; s : QuadDobl_Complex_Laurentials.Poly )
               return QuadDobl_Complex_Laurentials.Poly;

  -- DESCRIPTION :
  --   Permutes the unknowns in the term t or the polynonomial s,
  --   according to the permutation p.

  function "*" ( s : Standard_Complex_Poly_Systems.Poly_Sys;
                 p : Permutation )
               return Standard_Complex_Poly_Systems.Poly_Sys;
  function "*" ( s : DoblDobl_Complex_Poly_Systems.Poly_Sys;
                 p : Permutation )
               return DoblDobl_Complex_Poly_Systems.Poly_Sys;
  function "*" ( s : QuadDobl_Complex_Poly_Systems.Poly_Sys;
                 p : Permutation )
               return QuadDobl_Complex_Poly_Systems.Poly_Sys;

  -- DESCRIPTION :
  --   s*p permutes the unknowns in the individual polynomials,
  --   in standard double, double double, or quad double precision.

  function "*" ( p : Permutation;
                 s : Standard_Complex_Poly_Systems.Poly_Sys )
               return Standard_Complex_Poly_Systems.Poly_Sys;
  function "*" ( p : Permutation;
                 s : DoblDobl_Complex_Poly_Systems.Poly_Sys )
               return DoblDobl_Complex_Poly_Systems.Poly_Sys;
  function "*" ( p : Permutation;
                 s : QuadDobl_Complex_Poly_Systems.Poly_Sys )
               return QuadDobl_Complex_Poly_Systems.Poly_Sys;

  -- DESCRIPTION :
  --   p*s permutes the equations in the system,
  --   in standard double, double double, or quad double precision.
  --   Watch out for sharing by this second type of operation!

  function "*" ( s : Standard_Complex_Laur_Systems.Laur_Sys;
                 p : Permutation )
               return Standard_Complex_Laur_Systems.Laur_Sys;
  function "*" ( p : Permutation;
                 s : Standard_Complex_Laur_Systems.Laur_Sys )
               return Standard_Complex_Laur_Systems.Laur_Sys;
  function "*" ( s : DoblDobl_Complex_Laur_Systems.Laur_Sys;
                 p : Permutation )
               return DoblDobl_Complex_Laur_Systems.Laur_Sys;
  function "*" ( p : Permutation;
                 s : DoblDobl_Complex_Laur_Systems.Laur_Sys )
               return DoblDobl_Complex_Laur_Systems.Laur_Sys;
  function "*" ( s : QuadDobl_Complex_Laur_Systems.Laur_Sys;
                 p : Permutation )
               return QuadDobl_Complex_Laur_Systems.Laur_Sys;
  function "*" ( p : Permutation;
                 s : QuadDobl_Complex_Laur_Systems.Laur_Sys )
               return QuadDobl_Complex_Laur_Systems.Laur_Sys;

  -- DESCRIPTION :
  --   s*p permutes the unknowns in the individual polynomials.
  --   p*s permutes the equations in the system.
  --   Watch out for sharing by this second type of operation!

end Permute_Operations;
