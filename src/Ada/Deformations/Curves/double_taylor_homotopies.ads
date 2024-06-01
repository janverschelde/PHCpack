with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Complex_Numbers;
with Standard_Integer_Vectors;
with Standard_Floating_Vectors;

package Double_Taylor_Homotopies is

-- DESCRIPTION :
--   A Taylor monomial homotopy, or Taylor homotopy for short,
--   is a homotopy where the coefficients of the monomials are 
--   Taylor developments of functions of the continuation parameter t.
--   One example of such function is t^e, where e is a real exponent.

  type Taylor_Monomial ( dim,deg : integer32 ) is record
    cst : Standard_Complex_Numbers.Complex_Number;
    pwr : double_float; -- power of t
    cff : Standard_Floating_Vectors.Vector(0..deg);
    exp : Standard_Integer_Vectors.Vector(1..dim);
  end record;
  type Link_to_Taylor_Monomial is access Taylor_Monomial;

  type Taylor_Monomial_Vector is
    array ( integer32 range <> ) of Link_to_Taylor_Monomial;
  type Link_to_Taylor_Monomial_Vector is access Taylor_Monomial_Vector;

  type Taylor_Homotopy is
    array ( integer32 range <> ) of Link_to_Taylor_Monomial_Vector;
  type Link_to_Taylor_Homotopy is access Taylor_Homotopy;

  function Make ( deg : integer32; alpha,point : double_float;
                  cst : Standard_Complex_Numbers.Complex_Number;
                  exp : Standard_Integer_Vectors.Vector )
                return Taylor_Monomial;
  function Make ( deg : integer32; alpha,point : double_float;
                  cst : Standard_Complex_Numbers.Complex_Number;
                  exp : Standard_Integer_Vectors.Vector )
                return Link_to_Taylor_Monomial;

  -- DESCRIPTION :
  --   Given the truncation degree deg, the power alpha,
  --   the point at which to develop the Taylor series,
  --   a complex constant cst, and the exponents exp,
  --   returns a Taylor monomial, wrapping the computation
  --   of the Taylor series development.

  procedure Clear ( tm : in out Link_to_Taylor_Monomial );
  procedure Clear ( tmv : in out Taylor_Monomial_Vector );
  procedure Clear ( tmv : in out Link_to_Taylor_Monomial_Vector );
  procedure Clear ( th : in out Taylor_Homotopy );
  procedure Clear ( th : in out Link_to_Taylor_Homotopy );

  -- DESCRIPTION :
  --   Deallocates the data occupied by tm.

end Double_Taylor_Homotopies;
