with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with DoblDobl_Complex_Matrices;          use DoblDobl_Complex_Matrices;
with DoblDobl_Complex_Polynomials;       use DoblDobl_Complex_Polynomials;
with DoblDobl_Complex_Poly_Systems;      use DoblDobl_Complex_Poly_Systems;

package DoblDobl_Linear_Reduction is

-- DESCRIPTION :
--   The application of row reduction on the coefficient matrix of
--   a polynomial system may lead to a system with a lower total degree.
--   The operations in this package are utilities to build the
--   coefficient matrix and to reconstruct the system after row reduction.
--   The row reduction is executed in double double precision.

  mach_eps : constant double_float := 1.0E-30;

  type Degrees_Array is array ( integer32 range <> ) of Degrees;
  type Terms_Array   is array ( integer32 range <> ) of Term;
  type Boolean_Array is array ( integer32 range <> ) of boolean;

  procedure Pop_First_Term ( p : in out Poly; t : in out Term );

  -- DESCRIPTION :
  --   The term on return is the leading term of p.
  --   This term is removed from p.

  procedure Leading_Terms ( p : in out Poly_Sys; ta : in out Terms_Array );

  -- DESCRIPTION :
  --   Puts the leading terms of the polynomials in p in the array ta.
  --   The leading terms are removed afterwards.

  procedure Find_Max ( ta : in Terms_Array; index : in out Boolean_Array;
                       stop : in out boolean );

  procedure Update ( p : in out Poly_Sys; n : in natural32;
                     ta : in out Terms_Array; da : in out Degrees_Array;
                     nda,cnt : in out natural32; mat : in out Matrix;
                     stop : in out boolean );

  -- DESCRIPTION :
  --   Updates the coefficient matrix.

  procedure Coefficient_Matrix
                ( p : in Poly_Sys; mat : in out Matrix;
                  da : in out Degrees_Array; nda : in out natural32;
                  diagonal : in out boolean );

  -- DESCRIPTION :
  --   Constructs the coefficient matrix of the polynomial system.
  --   Stops when the system is diagonal.

  -- REQUIRED :
  --   mat'range(1) = p'range, mat'range(2) = 1..Sum_Number_of_Terms(p),
  --   da'range = mat'range(2).

  -- ON ENTRY :
  --   p          a polynomial system.

  -- ON RETURN :
  --   mat        coefficient matrix, up to column nda filled up;
  --   da         da(1..nda) collects the different terms in the system;
  --   nda        number of different terms;
  --   diagonal   true if the leading terms are all different.

  procedure Coefficient_Matrix
                ( p : in Poly_Sys; mat : in out Matrix;
                  da : in out Degrees_Array; nda : in out natural32 );

  -- DESCRIPTION :
  --   Constructs the coefficient matrix of the polynomial system.

  -- REQUIRED :
  --   mat'range(1) = p'range, mat'range(2) = 1..Sum_Number_of_Terms(p),
  --   da'range = mat'range(2).

  -- ON ENTRY :
  --   p          a polynomial system.

  -- ON RETURN :
  --   mat        coefficient matrix, up to column nda filled up;
  --   da         da(1..nda) collects the different terms in the system;
  --   nda        number of different terms.

  procedure Make_Polynomial_System
                ( p : in out Poly_Sys; mat : in Matrix;
                  da : in Degrees_Array; nda : in natural32;
                  inconsistent,infinite : out boolean );

  -- DESCRIPTION :
  --   Reconstructs the polynomial p from the coefficient matrix in mat
  --   and the degrees in the array da.

  function Sum_Number_of_Terms ( p : Poly_Sys ) return natural32;

  -- DESCRIPTION :
  --   Returns the sum of the number of terms of every polynomial in p.

end DoblDobl_Linear_Reduction;
