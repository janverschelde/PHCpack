with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with DoblDobl_Complex_Numbers;          use DoblDobl_Complex_Numbers;
with Standard_Natural_Vectors;
with Standard_Natural_VecVecs;
with Standard_Integer_Vectors;
with Standard_Integer_VecVecs;
with DoblDobl_Complex_Vectors;
with DoblDobl_Complex_VecVecs;
with DoblDobl_Complex_Matrices;
with Lists_of_Integer_Vectors;          use Lists_of_Integer_Vectors;
with DoblDobl_Complex_Polynomials;
with DoblDobl_Complex_Laurentials;
with DoblDobl_Complex_Poly_Systems;     use DoblDobl_Complex_Poly_Systems;
with DoblDobl_Complex_Laur_Systems;     use DoblDobl_Complex_Laur_Systems;

package DoblDobl_Polynomial_Flatteners is

-- DESCRIPTION :
--   A flattened representation of a polynomial system consists
--   of a coefficient matrix and a monomial vector.
--   The monomial vector is defined by a vector of exponents.
--   For a dense polynomial system, we use a matrix for the coefficients.
--   Columns of the coefficient matrix are indexed by the exponents
--   and the rows contain the coefficients of the polynomials in the system.
--   For a sparse polynomial, we use a vector of coefficient vectors,
--   with corresponding indices to the exponent vectors.

-- CONSTRUCTORS for SUPPORTS :

  function Number_of_Terms ( p : Poly_Sys ) return natural32;
  function Number_of_Terms ( p : Laur_Sys ) return natural32;

  -- DESCRIPTION :
  --   Returns the sum of the number of terms of the polynomials in p,
  --   as an upper bound on the number of columns in the coefficient matrix.

  procedure Update_Supports ( s,s_last : in out List;
                              p : in DoblDobl_Complex_Polynomials.Poly );
  procedure Update_Supports ( s,s_last : in out List;
                              p : in DoblDobl_Complex_Laurentials.Poly );

  -- DESCRIPTION :
  --   Updates the supports in s with the exponent vectors in p.

  -- ON ENTRY :
  --   s        current list of supports, all vectors are distinct;
  --   s_last   points to the last element in the list s;
  --   p        a polynomial in several variables.

  -- ON RETURN :
  --   s        updated list of supports, all vectors are distinct
  --            and contain the exponent vectors of p;
  --   s_last   points to the last element in the list s.

  function Distinct_Supports ( p : Poly_Sys ) return List;
  function Distinct_Supports ( p : Laur_Sys ) return List;

  -- DESCRIPTION :
  --   Returns a list of integer vectors with the support of p,
  --   in which every exponent vector is listed only once.

-- CONSTRUCTORS for the DENSE case :

  procedure Update_Coefficient_Matrix
              ( A : in out DoblDobl_Complex_Matrices.Matrix;
                i : in integer32;
                s : in Standard_Integer_VecVecs.VecVec;
                p : in DoblDobl_Complex_Polynomials.Poly );
  procedure Update_Coefficient_Matrix
              ( A : in out DoblDobl_Complex_Matrices.Matrix;
                i : in integer32;
                s : in Standard_Integer_VecVecs.VecVec;
                p : in DoblDobl_Complex_Laurentials.Poly );

  -- DESCRIPTION :
  --   Updates the i-th row of the coefficient matrix A
  --   with the coefficients of the polynomial p.
  --   In s are the lexicographically sorted exponent vectors.

  function Coefficient_Matrix 
             ( p : Poly_Sys; s : Standard_Integer_VecVecs.VecVec ) 
             return DoblDobl_Complex_Matrices.Matrix;
  function Coefficient_Matrix 
             ( p : Laur_Sys; s : Standard_Integer_VecVecs.VecVec ) 
             return DoblDobl_Complex_Matrices.Matrix;

  -- DESCRIPTION :
  --   Computes the coefficient matrix representation of the 
  --   polynomial system in p, given the lexicographically sorted
  --   list of all distinct exponent vectors in s.

  procedure Flatten ( p : in Poly_Sys;
                      s : out Standard_Integer_VecVecs.Link_to_VecVec;
                      c : out DoblDobl_Complex_Matrices.Link_to_Matrix );
  procedure Flatten ( p : in Laur_Sys;
                      s : out Standard_Integer_VecVecs.Link_to_VecVec;
                      c : out DoblDobl_Complex_Matrices.Link_to_Matrix );

  -- DESCRIPTION :
  --   Returns the flattened representation of the polynomial system p
  --   in the vector of exponents s and coefficient matrix c.

-- CONSTRUCTORS for the SPARSE case :

  procedure Coefficients_of_Support
              ( p : in DoblDobl_Complex_Polynomials.Poly;
                s : in Standard_Integer_VecVecs.VecVec;
                c : out DoblDobl_Complex_Vectors.Link_to_Vector;
                k : out Standard_Natural_Vectors.Link_to_Vector );
  procedure Coefficients_of_Support
              ( p : in DoblDobl_Complex_Laurentials.Poly;
                s : in Standard_Integer_VecVecs.VecVec;
                c : out DoblDobl_Complex_Vectors.Link_to_Vector;
                k : out Standard_Natural_Vectors.Link_to_Vector );

  -- DESCRIPTION :
  --   Returns in c the coefficients of p, ordered along exponents in s.
  --   Indices mapping coefficients to exponents are in k.

  -- ON ENTRY :
  --   p        a (Laurent) polynomial in several variables;
  --   s        lexicographically ordered vector of exponents in p.

  -- ON RETURN :
  --   c        nonzero coefficients of p occurring with exponents in s;
  --   k        p is the sum of c(i)*x^s(k(i)), for i in c'range.

  procedure Coefficients_of_Supports
              ( p : in Poly_Sys;
                s : in Standard_Integer_VecVecs.VecVec;
                c : out DoblDobl_Complex_VecVecs.VecVec;
                k : out Standard_Natural_VecVecs.VecVec );
  procedure Coefficients_of_Supports
              ( p : in Laur_Sys;
                s : in Standard_Integer_VecVecs.VecVec;
                c : out DoblDobl_Complex_VecVecs.VecVec;
                k : out Standard_Natural_VecVecs.VecVec );

  -- DESCRIPTION :
  --   Calls Coefficients_of_Support to all polynomial in p,
  --   defining the coefficients in c and indices in k.

  -- REQUIRED : c'range = k'range = p'range.

-- EVALUATORS :

  function Eval ( v : Standard_Integer_Vectors.Vector;
                  x : DoblDobl_Complex_Vectors.Vector )
                return Complex_Number;

  -- DESCRIPTION :
  --   Returns the value of the monomial defined by the exponent
  --   vector v at the vector x.

  -- REQUIRED : x'range = v'range.

  function Compressed_Eval
             ( v : Standard_Integer_Vectors.Vector;
               x : DoblDobl_Complex_Vectors.Vector )
             return Complex_Number;

  -- DESCRIPTION :
  --   Given a compressed vector of monomial indices and powers,
  --   returns the value of the monomial at x.

  function Eval ( v : Standard_Integer_VecVecs.VecVec;
                  x : DoblDobl_Complex_Vectors.Vector )
                return DoblDobl_Complex_Vectors.Vector;

  -- DESCRIPTION :
  --   Evaluates the monomial vector defined by the exponents in v
  --   at the given vector x.

  -- REQUIRED : x'range must match the range of the vectors in v.

  function Factored_Eval
             ( v : Standard_Integer_VecVecs.VecVec;
               x : DoblDobl_Complex_Vectors.Vector )
             return DoblDobl_Complex_Vectors.Vector;

  -- DESCRIPTION :
  --   Evaluates the monomial vectors in factored form in v
  --   at the given x.

  -- REQUIRED : x'range fits in the range of the vectors in v.

  function Factored_Compressed_Eval
             ( v : Standard_Integer_VecVecs.VecVec;
               x : DoblDobl_Complex_Vectors.Vector )
             return DoblDobl_Complex_Vectors.Vector;

  -- DESCRIPTION :
  --   Evaluates the monomial vectors in factored compressed form in v
  --   at the given x.

  function Eval ( A : in DoblDobl_Complex_Matrices.Matrix;
                  v : in Standard_Integer_VecVecs.VecVec;
                  x : in DoblDobl_Complex_Vectors.Vector ) 
                return DoblDobl_Complex_Vectors.Vector;

  -- DESCRIPTION :
  --   Returns the value of the flattened polynomial system representation
  --   defined by the coefficient matrix A and the vector v of exponents 
  --   at the vector x.

  function Eval ( c : DoblDobl_Complex_Vectors.Vector;
                  v : Standard_Integer_VecVecs.VecVec;
                  k : Standard_Natural_Vectors.Vector;
                  x : in DoblDobl_Complex_Vectors.Vector )
                return Complex_Number;

  -- DESCRIPTION :
  --   Evaluates the polynomial with flattened coefficient representation
  --   in (c,v,k), as c(i)*x^v(k(i)), for i in c'range, at x.

  function Eval ( c,vx : DoblDobl_Complex_Vectors.Vector;
                  k : Standard_Natural_Vectors.Vector )
                return Complex_Number;

  -- DESCRIPTION :
  --   Given the values of the monomials at x in v and coefficients in c
  --   with indices mapping coefficients to monomials in k, on return is
  --   the value of the polynomial c(i)*x^vx(k(i)), for i in c'range.

  function Eval ( c : DoblDobl_Complex_VecVecs.VecVec;
                  v : Standard_Integer_VecVecs.VecVec;
                  k : Standard_Natural_VecVecs.VecVec;
                  x : DoblDobl_Complex_Vectors.Vector )
                return DoblDobl_Complex_Vectors.Vector;

  -- DESCRIPTION :
  --   Evaluates the flattened sparse representation (c,v,k) at x,
  --   evaluating first the monomials in v at x and then calling Eval
  --   on every coefficient vector c(i) with corresponding indices in k(i).

  function Factored_Eval
             ( c : DoblDobl_Complex_VecVecs.VecVec;
               v : Standard_Integer_VecVecs.VecVec;
               k : Standard_Natural_VecVecs.VecVec;
               x : DoblDobl_Complex_Vectors.Vector )
             return DoblDobl_Complex_Vectors.Vector;

  -- DESCRIPTION :
  --   Assumes the factored form of the monomials in v to evaluate
  --   the flattened coefficient representation in (c,v,k) at x.

  function Factored_Compressed_Eval
             ( c : DoblDobl_Complex_VecVecs.VecVec;
               v : Standard_Integer_VecVecs.VecVec;
               k : Standard_Natural_VecVecs.VecVec;
               x : DoblDobl_Complex_Vectors.Vector )
             return DoblDobl_Complex_Vectors.Vector;

  -- DESCRIPTION :
  --   Assumes the compressed factored form of the monomials in v to evaluate
  --   the flattened coefficient representation in (c,v,k) at x.

-- TEST PROCEDURES :

  procedure Spy ( A : in DoblDobl_Complex_Matrices.Matrix; 
                  v : Standard_Integer_VecVecs.VecVec );

  -- DESCRIPTION :
  --   Creates a "spy plot" of a complex matrix A.
  --   The columns of A are indexed by exponent vectors in v
  --   and the rows are coefficient vectors of polynomials.
  --   A nonzero coefficient is represented by a *
  --   otherwise 0 is written to screen.
  --   Because there are typically many more monomials than polynomials
  --   in the system, the columns are swapped with rows in the display.

  procedure Test_Eval ( p : in Poly_Sys;
                        A : in DoblDobl_Complex_Matrices.Matrix;
                        v : in Standard_Integer_VecVecs.VecVec );
  procedure Test_Eval ( p : in Laur_Sys;
                        A : in DoblDobl_Complex_Matrices.Matrix;
                        v : in Standard_Integer_VecVecs.VecVec );

  -- DESCRIPTION :
  --   Computes the evaluation of p at a randomly generated point
  --   via the flattened matrix representation of p.
  --   Compares the computed value with the value of p at the point.

  procedure Test_Eval ( p : in Poly_Sys;
                        c : in DoblDobl_Complex_VecVecs.VecVec;
                        v,fv,cfv : in Standard_Integer_VecVecs.VecVec;
                        k : in Standard_Natural_VecVecs.VecVec );
  procedure Test_Eval ( p : in Laur_Sys;
                        c : in DoblDobl_Complex_VecVecs.VecVec;
                        v : in Standard_Integer_VecVecs.VecVec;
                        k : in Standard_Natural_VecVecs.VecVec );

  -- DESCRIPTION :
  --   Computes the evaluation of p at a randomly generated point
  --   via the flattened sparse representation of p in (c,v,k),
  --   where fv is the factored form of v and cfv is the compressed form.
  --   Compares the computed value with the value of p at the point.

end DoblDobl_Polynomial_Flatteners;
