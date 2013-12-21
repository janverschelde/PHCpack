with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Vectors;           use Standard_Integer_Vectors;
with Standard_Complex_Polynomials;       use Standard_Complex_Polynomials;
with Standard_Integer32_Triangulations;  use Standard_Integer32_Triangulations;
with Integer_Mixed_Subdivisions;         use Integer_Mixed_Subdivisions;

package Minkowski_Polynomials is

-- DESCRIPTION :
--   This package allows the computation of the Minkowski-polynomial
--   of a tuple of polytopes (P1,P2,..,Pr).  This polynomial is the
--   expansion of the volume of a positive linear combination of the
--   polytopes in the tuple: vol_n(l1*P1 + l2*P2 + .. + lr*Pr), which
--   is a homogeneous polynomial of degree n in the coefficients l1,l2,..,lr,
--   according to Minkowski's theorem.

  function Minkowski_Polynomial ( n,r : natural32 ) return Poly;

  -- DESCRIPTION :
  --   Returns the structure of the Minkowski-polynomial, given the
  --   dimension and the number of different polytopes in the tuple.

  -- ON ENTRY :
  --   n           dimension of the polytopes before lifting;
  --   r           number of different polytopes in the tuple.

  -- ON RETURN :
  --   The structure of the Minkowski-polynomial, with all coefficients
  --   equal to one.

  procedure Minkowski_Polynomial
              ( p : in out Poly; t : in Triangulation; n : in natural32;
                mix : in Vector; mixsub : out Mixed_Subdivision );

  -- DESCRIPTION :
  --   Computes the coefficients of the Minkowski-polynomial, given its
  --   structure and based on the triangulation of the Cayley-polytope.
  --   On return, one also obtains the mixed subdivision, corresponding
  --   the given type of mixture.

  generic
    with procedure Process ( submix : in Vector; sub : in Mixed_Subdivision;
                             vol : out natural32 );
  procedure Minkowski_Polynomial_Subdivisions
                 ( p : in out Poly; t : in Triangulation; n : in natural32;
                   mix : in Vector; mixsub : out Mixed_Subdivision );

  -- DESCRIPTION :
  --   Computes the coefficients of the Minkowski-polynomial, given its
  --   structure and the triangulation of the Cayley polytope.
  --   The generic procedure returns the subdivision for a given type of
  --   mixture and asks to compute its volume.
  --   On return, one also obtains the mixed subdivision, corresponding
  --   the given type of mixture.

end Minkowski_Polynomials;
