with Standard_Integer_Numbers;            use Standard_Integer_Numbers;
with QuadDobl_Complex_Numbers;
with Standard_Integer_Vectors;
with QuadDobl_Complex_Vectors;
with QuadDobl_Complex_Matrices;
with QuadDobl_Complex_Polynomials;
with QuadDobl_Complex_Poly_Systems;
with QuadDobl_Complex_Circuits;
with QuadDobl_Speelpenning_Convolutions;

package QuadDobl_Circuit_Makers is

-- DESCRIPTION :
--   A circuit maker is one of the following kinds:
--   (1) make a random circuit for testing purposes,
--   (2) convert between circuit and polynomial types,
--   (3) make a circuit from a convolution circuit.
--   The package also offers useful procedures in testing.

  function Random_Indices
             ( dim : integer32 )
             return Standard_Integer_Vectors.Vector;

  -- DESCRIPTION :
  --   Returns a vector of indices between 1 and dim
  --   of variables participating in a product.
  --   The size of the vector on return is at least two.

  function Random_Indices
             ( dim,size : integer32 )
             return Standard_Integer_Vectors.Vector;

  -- DESCRIPTION :
  --   Returns a vector of indices of size,
  --   with numbers in 1..dim, without duplicates.
  --   The vector on return is sorted.

  -- REQUIRED : dim >= size.

  function Random_Complex_Circuit
             ( nbr,dim : integer32 )
             return QuadDobl_Complex_Circuits.Circuit;

  -- DESCRIPTION :
  --   Returns a random complex circuit with products of dimension dim
  --   and as many nonconstant coefficients as the number nbr.
  --   Every exponent vector in the circuit has at least size two.
  --   All exponents are either zero or one.

  function Random_Complex_Circuit
             ( nbr,dim,pwr : integer32 )
             return QuadDobl_Complex_Circuits.Circuit;
  function Random_Complex_Circuit
             ( nbr,dim,pwr : integer32 )
             return QuadDobl_Complex_Circuits.Link_to_Circuit;

  -- DESCRIPTION :
  --   Returns a random complex circuit with products of dimension dim,
  --   as many nonconstant coefficients as the number nbr,
  --   and pwr as the value for the higest power.

  function Random_Complex_Circuits
             ( neq,nbr,dim,pwr : integer32 )
             return QuadDobl_Complex_Circuits.Circuits;

  -- DESCRIPTION :
  --   Returns a sequence of neq circuits, with nbr terms,
  --   of dimension dim and with highest power equal to pwr.

  function Random_Complex_System
             ( neq,nbr,dim,pwr : integer32 )
             return QuadDobl_Complex_Circuits.System;
  function Random_Complex_System
             ( neq,nbr,dim,pwr : integer32 )
             return QuadDobl_Complex_Circuits.Link_to_System;

  -- DESCRIPTION :
  --   Returns a system of neq circuits, with nbr terms,
  --   of dimension dim and with highest power equal to pwr.

  function Make_Polynomial
             ( c : QuadDobl_Complex_Circuits.Circuit;
               index : boolean := false )
             return QuadDobl_Complex_Polynomials.Poly;

  -- DESCRIPTION :
  --   Returns the polynomial equivalent to the circuit c.
  --   If index, then c.xps should be considered as the index
  --   vector of the participating variables in each term,
  --   instead of the exponent vector.

  function Gradient ( p : QuadDobl_Complex_Polynomials.Poly;
                      x : QuadDobl_Complex_Vectors.Vector )
                    return QuadDobl_Complex_Vectors.Vector;

  -- DESCRIPTION :
  --   Straighforward computation of the gradient of p at x,
  --   for testing purposes.

  function Hessian ( p : QuadDobl_Complex_Polynomials.Poly;
                     x : QuadDobl_Complex_Vectors.Vector )
                   return QuadDobl_Complex_Matrices.Matrix;

  -- DESCRIPTION :
  --   Returns the Hessian matrix of p at x, for testing purposes.

  function Constant_Coefficient
             ( p : QuadDobl_Complex_Polynomials.Poly )
             return QuadDobl_Complex_Numbers.Complex_Number;

  -- DESCRIPTION :
  --   Returns the constant coefficient of the polynomial p.

  function Is_NonZero ( c : QuadDobl_Complex_Numbers.Complex_Number )
                   return integer32;

  -- DESCRIPTION :
  --   Returns 1 if the coefficient is nonzero,
  --   returns 0 otherwise.

  function Make_Complex_Circuit
             ( p : QuadDobl_Complex_Polynomials.Poly;
               verbose : in boolean := true )
             return QuadDobl_Complex_Circuits.Circuit;

  -- DESCRIPTION :
  --   Returns the circuit representation of the polynomial p.
  --   If verbose, then extra information is written to screen.
    
  function Make_Complex_System
             ( p : in QuadDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
               verbose : in boolean := true )
             return QuadDobl_Complex_Circuits.Link_to_System;

  -- DESCRIPTION :
  --   Returns the system of circuits defined by p.
  --   If verbose, then the tableau format of the system is written.

  procedure Write_Matrix ( A : in QuadDobl_Complex_Matrices.Matrix );

  -- DESCRIPTION :
  --   Writes the matrix A component-wise, with explicit indexing,
  --   for visual inspection of small matrices while testing.

-- FROM CONVOLUTION CIRCUITS TO COMPLEX CIRCUITS :

  function Make_Complex_Circuit
             ( c : QuadDobl_Speelpenning_Convolutions.Circuit )
             return QuadDobl_Complex_Circuits.Circuit;
  function Make_Complex_Circuit
             ( c : QuadDobl_Speelpenning_Convolutions.Link_to_Circuit )
             return QuadDobl_Complex_Circuits.Link_to_Circuit;

  -- DESCRIPTION :
  --   The circuit on return has the leading coefficient of every 
  --   power series coefficient of the given c and stores copies
  --   of the exponents, indices, and factors, as stored in c.
  --   All data structures are allocated.

  function Make_Complex_System
             ( s : QuadDobl_Speelpenning_Convolutions.System )
             return QuadDobl_Complex_Circuits.System;
  function Make_Complex_System
             ( s : QuadDobl_Speelpenning_Convolutions.Link_to_System )
             return QuadDobl_Complex_Circuits.Link_to_System;

  -- DESCRIPTION :
  --   Takes the leading coefficient for every power series coefficient
  --   in s.crc on input as the coefficient for the circuits on output.

end QuadDobl_Circuit_Makers;
