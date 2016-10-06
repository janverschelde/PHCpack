with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Natural_Vectors;
with Standard_Integer_Vectors;
with Standard_Floating_Vectors;
with Standard_Complex_Vectors;
with Standard_Complex_VecVecs;
with Lists_of_Integer_Vectors;
with Lists_of_Floating_Vectors;
with Arrays_of_Integer_Vector_Lists;
with Arrays_of_Floating_Vector_Lists;
with Standard_Complex_Polynomials;
with Standard_Complex_Laurentials;
with Standard_Complex_Poly_Systems;
with Standard_Complex_Laur_Systems;
with DoblDobl_Complex_Polynomials;
with DoblDobl_Complex_Laurentials;
with DoblDobl_Complex_Poly_Systems;
with DoblDobl_Complex_Laur_Systems;
with QuadDobl_Complex_Polynomials;
with QuadDobl_Complex_Laurentials;
with QuadDobl_Complex_Poly_Systems;
with QuadDobl_Complex_Laur_Systems;

package Supports_of_Polynomial_Systems is

-- DESCRIPTION :
--   The support of the polynomial collects the exponent vectors of those
--   monomials which occur with a nonzero coefficient of the polynomial.
--   The supports of a polynomial system are a tuple of supports of each
--   polynomial in the system.  The multitude of routines is due to the
--   versatility of types of systems: ordinary polynomials or Laurent
--   polynomials; and the integer and floating-points lists used to
--   represent the supports.

  function Create ( p : Standard_Complex_Polynomials.Poly )
                  return Lists_of_Integer_Vectors.List;
  function Create ( p : Standard_Complex_Laurentials.Poly )
                  return Lists_of_Integer_Vectors.List;
  function Create ( p : DoblDobl_Complex_Polynomials.Poly )
                  return Lists_of_Integer_Vectors.List;
  function Create ( p : DoblDobl_Complex_Laurentials.Poly )
                  return Lists_of_Integer_Vectors.List;
  function Create ( p : QuadDobl_Complex_Polynomials.Poly )
                  return Lists_of_Integer_Vectors.List;
  function Create ( p : QuadDobl_Complex_Laurentials.Poly )
                  return Lists_of_Integer_Vectors.List;

  -- DESCRIPTION : Returns the support of p.

  function Random_Complex_Polynomial
             ( s : Lists_of_Integer_Vectors.List )
             return Standard_Complex_Polynomials.Poly;

  -- DESCRIPTION :
  --   Returns a polynomial with random complex coefficients on
  --   the unit circle and with the given support in s.

  -- REQUIRED:
  --   All exponents in s are either zero or positive.

  function Random_Complex_Laurential
             ( s : Lists_of_Integer_Vectors.List )
             return Standard_Complex_Laurentials.Poly;

  -- DESCRIPTION :
  --   Returns a Laurent polynomial with random complex coefficients on
  --   the unit circle and with the given support in s.

  function Random_Complex_Laurent_System
             ( dim : integer32; s : Lists_of_Integer_Vectors.List )
             return Standard_Complex_Laur_Systems.Laur_Sys;

  -- DESCRIPTION :
  --   Returns a Laurent polynomial system of range 1..dim
  --   with random complex coefficients and support equal to s.

  function Random_Complex_Laurent_System
             ( dim : integer32; mix : Standard_Integer_Vectors.Vector;
               sup : Arrays_of_Integer_Vector_Lists.Array_of_Lists )
             return Standard_Complex_Laur_Systems.Laur_Sys;

  -- DESCRIPTION :
  --   Returns a Laurent polynomial system of range 1..dim
  --   with random complex coefficients and support equal to sup,
  --   and type of mixture as prescribed in the vector mix.

  -- REQUIRED : sup'range = mix'range and sum(mix) = dim.

  function Is_Equal ( v1 : Standard_Floating_Vectors.Vector;
                      v2 : Standard_Natural_Vectors.Vector ) return boolean;
  function Is_Equal ( v1 : Standard_Floating_Vectors.Vector;
                      v2 : Standard_Integer_Vectors.Vector ) return boolean;

  -- DESCRIPTION :
  --   The vectors v1 and v2 are equal if and only if all elements of
  --   the vector v2 converted into double_float equal v1.
  --   If v1 is longer than v2, then these elements will be ignored.

  function Is_In ( s : Lists_of_Floating_Vectors.List;
                   v : Standard_Natural_Vectors.Vector ) return boolean;
  function Is_In ( s : Lists_of_Floating_Vectors.List;
                   v : Standard_Integer_Vectors.Vector ) return boolean;

  -- DESCRIPTION :
  --   Returns true if the vector v belongs to the list s.

  function Select_Terms ( p : Standard_Complex_Polynomials.Poly;
                          s : Lists_of_Integer_Vectors.List )
                        return Standard_Complex_Polynomials.Poly;
  function Select_Terms ( p : Standard_Complex_Laurentials.Poly;
                          s : Lists_of_Integer_Vectors.List )
                        return Standard_Complex_Laurentials.Poly;
  function Select_Terms ( p : Standard_Complex_Polynomials.Poly;
                          s : Lists_of_Floating_Vectors.List )
                        return Standard_Complex_Polynomials.Poly;
  function Select_Terms ( p : Standard_Complex_Laurentials.Poly;
                          s : Lists_of_Floating_Vectors.List )
                        return Standard_Complex_Laurentials.Poly;

  function Select_Terms ( p : DoblDobl_Complex_Polynomials.Poly;
                          s : Lists_of_Integer_Vectors.List )
                        return DoblDobl_Complex_Polynomials.Poly;
  function Select_Terms ( p : DoblDobl_Complex_Polynomials.Poly;
                          s : Lists_of_Floating_Vectors.List )
                        return DoblDobl_Complex_Polynomials.Poly;
  function Select_Terms ( p : DoblDobl_Complex_Laurentials.Poly;
                          s : Lists_of_Integer_Vectors.List )
                        return DoblDobl_Complex_Laurentials.Poly;
  function Select_Terms ( p : DoblDobl_Complex_Laurentials.Poly;
                          s : Lists_of_Floating_Vectors.List )
                        return DoblDobl_Complex_Laurentials.Poly;

  function Select_Terms ( p : QuadDobl_Complex_Polynomials.Poly;
                          s : Lists_of_Integer_Vectors.List )
                        return QuadDobl_Complex_Polynomials.Poly;
  function Select_Terms ( p : QuadDobl_Complex_Polynomials.Poly;
                          s : Lists_of_Floating_Vectors.List )
                        return QuadDobl_Complex_Polynomials.Poly;
  function Select_Terms ( p : QuadDobl_Complex_Laurentials.Poly;
                          s : Lists_of_Integer_Vectors.List )
                        return QuadDobl_Complex_Laurentials.Poly;
  function Select_Terms ( p : QuadDobl_Complex_Laurentials.Poly;
                          s : Lists_of_Floating_Vectors.List )
                        return QuadDobl_Complex_Laurentials.Poly;

  -- DESCRIPTION :
  --   Returns those terms in p whose exponent vector occurs in the list s.

  function Create ( p : Standard_Complex_Poly_Systems.Poly_Sys )
                  return Arrays_of_Integer_Vector_Lists.Array_of_Lists;
  function Create ( p : Standard_Complex_Laur_Systems.Laur_Sys )
                  return Arrays_of_Integer_Vector_Lists.Array_of_Lists;
  function Create ( p : DoblDobl_Complex_Poly_Systems.Poly_Sys )
                  return Arrays_of_Integer_Vector_Lists.Array_of_Lists;
  function Create ( p : DoblDobl_Complex_Laur_Systems.Laur_Sys )
                  return Arrays_of_Integer_Vector_Lists.Array_of_Lists;
  function Create ( p : QuadDobl_Complex_Poly_Systems.Poly_Sys )
                  return Arrays_of_Integer_Vector_Lists.Array_of_Lists;
  function Create ( p : QuadDobl_Complex_Laur_Systems.Laur_Sys )
                  return Arrays_of_Integer_Vector_Lists.Array_of_Lists;

  -- DESCRIPTION :
  --   Returns the supports of the polynomial system.

  function Select_Terms ( p : Standard_Complex_Poly_Systems.Poly_Sys;
                          s : Arrays_of_Integer_Vector_Lists.Array_of_Lists )
                        return Standard_Complex_Poly_Systems.Poly_Sys;
  function Select_Terms ( p : Standard_Complex_Poly_Systems.Poly_Sys;
                          m : Standard_Integer_Vectors.Vector;
                          s : Arrays_of_Integer_Vector_Lists.Array_of_Lists )
                        return Standard_Complex_Poly_Systems.Poly_Sys;
  function Select_Terms ( p : Standard_Complex_Laur_Systems.Laur_Sys;
                          s : Arrays_of_Integer_Vector_Lists.Array_of_Lists )
                        return Standard_Complex_Laur_Systems.Laur_Sys;
  function Select_Terms ( p : Standard_Complex_Laur_Systems.Laur_Sys;
                          m : Standard_Integer_Vectors.Vector;
                          s : Arrays_of_Integer_Vector_Lists.Array_of_Lists )
                        return Standard_Complex_Laur_Systems.Laur_Sys;

  -- DESCRIPTION :
  --   Returns those terms in each polynomial p(i) whose exponent vector
  --   occurs in the list s(i), for i in p'range = s'range.

  function Select_Terms ( p : Standard_Complex_Poly_Systems.Poly_Sys;
                          s : Arrays_of_Floating_Vector_Lists.Array_of_Lists )
                        return Standard_Complex_Poly_Systems.Poly_Sys;
  function Select_Terms ( p : Standard_Complex_Poly_Systems.Poly_Sys;
                          m : Standard_Integer_Vectors.Vector; 
                          s : Arrays_of_Floating_Vector_Lists.Array_of_Lists )
                        return Standard_Complex_Poly_Systems.Poly_Sys;
  function Select_Terms ( p : Standard_Complex_Laur_Systems.Laur_Sys;
                          s : Arrays_of_Floating_Vector_Lists.Array_of_Lists )
                        return Standard_Complex_Laur_Systems.Laur_Sys;
  function Select_Terms ( p : Standard_Complex_Laur_Systems.Laur_Sys;
                          m : Standard_Integer_Vectors.Vector;
                          s : Arrays_of_Floating_Vector_Lists.Array_of_Lists )
                        return Standard_Complex_Laur_Systems.Laur_Sys;

  function Select_Terms ( p : DoblDobl_Complex_Poly_Systems.Poly_Sys;
                          s : Arrays_of_Integer_Vector_Lists.Array_of_Lists )
                        return DoblDobl_Complex_Poly_Systems.Poly_Sys;
  function Select_Terms ( p : DoblDobl_Complex_Poly_Systems.Poly_Sys;
                          m : Standard_Integer_Vectors.Vector;
                          s : Arrays_of_Integer_Vector_Lists.Array_of_Lists )
                        return DoblDobl_Complex_Poly_Systems.Poly_Sys;
  function Select_Terms ( p : DoblDobl_Complex_Poly_Systems.Poly_Sys;
                          s : Arrays_of_Floating_Vector_Lists.Array_of_Lists )
                        return DoblDobl_Complex_Poly_Systems.Poly_Sys;
  function Select_Terms ( p : DoblDobl_Complex_Poly_Systems.Poly_Sys;
                          m : Standard_Integer_Vectors.Vector; 
                          s : Arrays_of_Floating_Vector_Lists.Array_of_Lists )
                        return DoblDobl_Complex_Poly_Systems.Poly_Sys;
  function Select_Terms ( p : DoblDobl_Complex_Laur_Systems.Laur_Sys;
                          s : Arrays_of_Floating_Vector_Lists.Array_of_Lists )
                        return DoblDobl_Complex_Laur_Systems.Laur_Sys;
  function Select_Terms ( p : DoblDobl_Complex_Laur_Systems.Laur_Sys;
                          s : Arrays_of_Integer_Vector_Lists.Array_of_Lists )
                        return DoblDobl_Complex_Laur_Systems.Laur_Sys;
  function Select_Terms ( p : DoblDobl_Complex_Laur_Systems.Laur_Sys;
                          m : Standard_Integer_Vectors.Vector;
                          s : Arrays_of_Floating_Vector_Lists.Array_of_Lists )
                        return DoblDobl_Complex_Laur_Systems.Laur_Sys;

  function Select_Terms ( p : QuadDobl_Complex_Poly_Systems.Poly_Sys;
                          s : Arrays_of_Integer_Vector_Lists.Array_of_Lists )
                        return QuadDobl_Complex_Poly_Systems.Poly_Sys;
  function Select_Terms ( p : QuadDobl_Complex_Poly_Systems.Poly_Sys;
                          m : Standard_Integer_Vectors.Vector;
                          s : Arrays_of_Integer_Vector_Lists.Array_of_Lists )
                        return QuadDobl_Complex_Poly_Systems.Poly_Sys;
  function Select_Terms ( p : QuadDobl_Complex_Poly_Systems.Poly_Sys;
                          s : Arrays_of_Floating_Vector_Lists.Array_of_Lists )
                        return QuadDobl_Complex_Poly_Systems.Poly_Sys;
  function Select_Terms ( p : QuadDobl_Complex_Poly_Systems.Poly_Sys;
                          m : Standard_Integer_Vectors.Vector; 
                          s : Arrays_of_Floating_Vector_Lists.Array_of_Lists )
                        return QuadDobl_Complex_Poly_Systems.Poly_Sys;
  function Select_Terms ( p : QuadDobl_Complex_Laur_Systems.Laur_Sys;
                          s : Arrays_of_Floating_Vector_Lists.Array_of_Lists )
                        return QuadDobl_Complex_Laur_Systems.Laur_Sys;
  function Select_Terms ( p : QuadDobl_Complex_Laur_Systems.Laur_Sys;
                          s : Arrays_of_Integer_Vector_Lists.Array_of_Lists )
                        return QuadDobl_Complex_Laur_Systems.Laur_Sys;
  function Select_Terms ( p : QuadDobl_Complex_Laur_Systems.Laur_Sys;
                          m : Standard_Integer_Vectors.Vector;
                          s : Arrays_of_Floating_Vector_Lists.Array_of_Lists )
                        return QuadDobl_Complex_Laur_Systems.Laur_Sys;

  -- DESCRIPTION :
  --   Returns those terms of the polynomials in p whose exponent vectors
  --   occur in the lists s.  If the type of mixture m is omitted,
  --   then s'range must equal p'range, otherwise the type of mixture
  --   is taken into account when selecting terms from p.

  function Select_Lifted
               ( p : Standard_Complex_Poly_Systems.Poly_Sys;
                 mix : Standard_Integer_Vectors.Vector;
                 lifsup : Arrays_of_Floating_Vector_Lists.Array_of_Lists )
               return Standard_Complex_Poly_Systems.Poly_Sys;
  function Select_Lifted
               ( p : Standard_Complex_Laur_Systems.Laur_Sys;
                 mix : Standard_Integer_Vectors.Vector;
                 lifsup : Arrays_of_Floating_Vector_Lists.Array_of_Lists )
               return Standard_Complex_Laur_Systems.Laur_Sys;

  -- DESCRIPTION :
  --   Returns the system whose terms have exponents that occur in the
  --   lifted supports.

-- PROCEDURES without extra memory allocation :

  procedure Select_Coefficients
              ( p : in Standard_Complex_Polynomials.Poly;
                s : in Lists_of_Integer_Vectors.List;
                dim : in natural32;
                deg : in Standard_Complex_Polynomials.Degrees;
                cff : out Standard_Complex_Vectors.Vector );
  procedure Select_Coefficients
              ( p : in Standard_Complex_Laurentials.Poly;
                s : in Lists_of_Integer_Vectors.List;
                dim : in natural32;
                deg : in Standard_Complex_Laurentials.Degrees;
                cff : out Standard_Complex_Vectors.Vector );
  procedure Select_Coefficients
              ( p : in Standard_Complex_Polynomials.Poly;
                s : in Lists_of_Floating_Vectors.List;
                dim : in natural32;
                deg : in Standard_Complex_Polynomials.Degrees;
                cff : out Standard_Complex_Vectors.Vector );
  procedure Select_Coefficients
              ( p : in Standard_Complex_Laurentials.Poly;
                s : in Lists_of_Floating_Vectors.List;
                dim : in natural32;
                deg : in Standard_Complex_Laurentials.Degrees;
                cff : out Standard_Complex_Vectors.Vector );

  -- DESCRIPTION :
  --   Returns the coefficients of p with exponent vectors in s.

  -- REQUIRED : deg'range = 1..dim and cff'range = 1..Length_Of(s).

  -- ON ENTRY :
  --   p        a polynomial with complex coefficients in dim variables;
  --   s        a list of exponent vectors of range at least 1..dim,
  --            note that lifting values will be ignored;
  --   dim      the number of variables in the polynomial p;
  --   deg      work space to select the coefficients of p,
  --            points to a natural vector of range 1..dim.
 
  -- ON RETURN :
  --   cff      cff(k) is the coefficient in p that has as exponents
  --            the k-th element in s.

  procedure Select_Coefficients
              ( p : in Standard_Complex_Poly_Systems.Poly_Sys;
                s : in Arrays_of_Integer_Vector_Lists.Array_of_Lists;
                dim : in natural32;
                deg : in Standard_Complex_Polynomials.Degrees;
                cff : in Standard_Complex_VecVecs.VecVec );
  procedure Select_Coefficients
              ( p : in Standard_Complex_Laur_Systems.Laur_Sys;
                s : in Arrays_of_Integer_Vector_Lists.Array_of_Lists;
                dim : in natural32;
                deg : in Standard_Complex_Laurentials.Degrees;
                cff : in Standard_Complex_VecVecs.VecVec );
  procedure Select_Coefficients
              ( p : in Standard_Complex_Poly_Systems.Poly_Sys;
                s : in Arrays_of_Floating_Vector_Lists.Array_of_Lists;
                dim : in natural32;
                deg : in Standard_Complex_Polynomials.Degrees;
                cff : in Standard_Complex_VecVecs.VecVec );
  procedure Select_Coefficients
              ( p : in Standard_Complex_Laur_Systems.Laur_Sys;
                s : in Arrays_of_Floating_Vector_Lists.Array_of_Lists;
                dim : in natural32;
                deg : in Standard_Complex_Laurentials.Degrees;
                cff : in Standard_Complex_VecVecs.VecVec );

  -- DESCRIPTION :
  --   Returns the coefficients of p with exponents vectors in s,
  --   for the fully mixed case, i.e.: s'range = p'range.

  -- REQUIRED : deg'range = 1..dim and cff'range(i) = 1..Length_Of(s(i)),
  --   for i in s'range.

  -- ON ENTRY :
  --   p        a system with complex coefficients in dim variables;
  --   s        a list of exponent vectors of range at least 1..dim,
  --            note that lifting values will be ignored;
  --   dim      the number of variables in the polynomial p;
  --   deg      work space to select the coefficients of p,
  --            points to a natural vector of range 1..dim.
 
  -- ON RETURN :
  --   cff      cff(i)(k) is the coefficient in p that has as exponents
  --            the k-th element in s(i), for i in s'range.

  procedure Select_Coefficients
              ( p : in Standard_Complex_Poly_Systems.Poly_Sys;
                m : in Standard_Integer_Vectors.Vector;
                s : in Arrays_of_Integer_Vector_Lists.Array_of_Lists;
                dim : in natural32;
                deg : in Standard_Complex_Polynomials.Degrees;
                cff : in Standard_Complex_VecVecs.VecVec );
  procedure Select_Coefficients
              ( p : in Standard_Complex_Laur_Systems.Laur_Sys;
                m : in Standard_Integer_Vectors.Vector;
                s : in Arrays_of_Integer_Vector_Lists.Array_of_Lists;
                dim : in natural32;
                deg : in Standard_Complex_Laurentials.Degrees;
                cff : in Standard_Complex_VecVecs.VecVec );
  procedure Select_Coefficients
              ( p : in Standard_Complex_Poly_Systems.Poly_Sys;
                m : in Standard_Integer_Vectors.Vector;
                s : in Arrays_of_Floating_Vector_Lists.Array_of_Lists;
                dim : in natural32;
                deg : in Standard_Complex_Polynomials.Degrees;
                cff : in Standard_Complex_VecVecs.VecVec );
  procedure Select_Coefficients
              ( p : in Standard_Complex_Laur_Systems.Laur_Sys;
                m : in Standard_Integer_Vectors.Vector;
                s : in Arrays_of_Floating_Vector_Lists.Array_of_Lists;
                dim : in natural32;
                deg : in Standard_Complex_Laurentials.Degrees;
                cff : in Standard_Complex_VecVecs.VecVec );

  -- DESCRIPTION :
  --   Returns the coefficients of p with exponents vectors in s,
  --   for the semi-mixed case, i.e.: s'range = m'range.

  -- REQUIRED : deg'range = 1..dim and cff'range(i) = 1..Length_Of(s(j)),
  --   for j corresponding to the correct i along the mixture in m.

  -- ON ENTRY :
  --   p        a system with complex coefficients in dim variables,
  --            polynomials with the same support appear consecutively;
  --   m        m(i) indicates how many times the i-th support occurs
  --   s        a list of exponent vectors of range at least 1..dim,
  --            note that lifting values will be ignored;
  --   dim      the number of variables in the polynomial p;
  --   deg      work space to select the coefficients of p,
  --            points to a natural vector of range 1..dim.
 
  -- ON RETURN :
  --   cff      cff(i)(k) is the coefficient in p(i) that has as exponents
  --            the k-th element in s(j), for j corresponding the mixture.

end Supports_of_Polynomial_Systems;
