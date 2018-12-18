with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Monomial_Vectors;
with DoblDobl_Monomial_Vectors;
with QuadDobl_Monomial_Vectors;

package Random_Monomial_Vectors is

-- DESCRIPTION :
--   A random monomial vector is a vector of random monomials.
--   Each random monomial is defined by a random complex coefficient
--   and a vector of random natural numbers.

  function Standard_Random_Monomial_Vector
             ( size,dim : integer32; expmax : natural32;
               verbose : boolean := false )  
             return Standard_Monomial_Vectors.Monomial_Vector;

  -- DESCRIPTION :
  --   On return is a vector of range 1..size.  Each entry contains a random
  --   monomial of dimension dim with exponents in the range 0..expmax.
  --   The coefficients are double complex numbers.
  --   If verbose, then the coefficient and exponents are shown.

  -- REQUIRED : expmax > 0.

  function DoblDobl_Random_Monomial_Vector
             ( size,dim : integer32; expmax : natural32;
               verbose : boolean := false )  
             return DoblDobl_Monomial_Vectors.Monomial_Vector;

  -- DESCRIPTION :
  --   On return is a vector of range 1..size.  Each entry contains a random
  --   monomial of dimension dim with exponents in the range 0..expmax.
  --   The coefficients are double double complex numbers.
  --   If verbose, then the coefficient and exponents are shown.

  -- REQUIRED : expmax > 0.

  function QuadDobl_Random_Monomial_Vector
             ( size,dim : integer32; expmax : natural32;
               verbose : boolean := false )  
             return QuadDobl_Monomial_Vectors.Monomial_Vector;

  -- DESCRIPTION :
  --   On return is a vector of range 1..size.  Each entry contains a random
  --   monomial of dimension dim with exponents in the range 0..expmax.
  --   The coefficients are quad double complex numbers.
  --   If verbose, then the coefficient and exponents are shown.

  -- REQUIRED : expmax > 0.

  function Standard_Random_Polynomial
             ( size,dim : integer32; expmax : natural32;
               verbose : boolean := false )  
             return Standard_Monomial_Vectors.Polynomial;
  function Standard_Random_Polynomial
             ( size,dim : integer32; expmax : natural32;
               verbose : boolean := false )  
             return Standard_Monomial_Vectors.Link_to_Polynomial;

  -- DESCRIPTION :
  --   On return is a polynomial with size+1 monomials,
  --   the +1 corresponds to the random constant coefficient.
  --   The monomials are generated with Standard_Random_Monomial_Vector.

  -- REQUIRED : expmax > 0.

  function DoblDobl_Random_Polynomial
             ( size,dim : integer32; expmax : natural32;
               verbose : boolean := false )  
             return DoblDobl_Monomial_Vectors.Polynomial;
  function DoblDobl_Random_Polynomial
             ( size,dim : integer32; expmax : natural32;
               verbose : boolean := false )  
             return DoblDobl_Monomial_Vectors.Link_to_Polynomial;

  -- DESCRIPTION :
  --   On return is a polynomial with size+1 monomials,
  --   the +1 corresponds to the random constant coefficient.
  --   The monomials are generated with DoblDobl_Random_Monomial_Vector.

  -- REQUIRED : expmax > 0.

  function QuadDobl_Random_Polynomial
             ( size,dim : integer32; expmax : natural32;
               verbose : boolean := false )  
             return QuadDobl_Monomial_Vectors.Polynomial;
  function QuadDobl_Random_Polynomial
             ( size,dim : integer32; expmax : natural32;
               verbose : boolean := false )  
             return QuadDobl_Monomial_Vectors.Link_to_Polynomial;

  -- DESCRIPTION :
  --   On return is a polynomial with size+1 monomials,
  --   the +1 corresponds to the random constant coefficient.
  --   The monomials are generated with QuadDobl_Random_Monomial_Vector.

  -- REQUIRED : expmax > 0.
 
end Random_Monomial_Vectors;
