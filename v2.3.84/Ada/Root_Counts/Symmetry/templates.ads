with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Natural_Vectors;           use Standard_Natural_Vectors;
with Lists_of_Integer_Vectors;           use Lists_of_Integer_Vectors;

package Templates is

-- DESCRIPTION :
--   This package provides the basic data abstraction to be used for the
--   construction of a symmetric linear-product start system.

-- DATA STRUCTURES :

  -- A template for an n-dimensional system is made of a table
  -- of natural vectors h(*).
  -- If h(i) = 0, then this coefficient h(i) will be zero,
  -- elsif h(i) = j /= 0,
  --     then this coefficient h(i) will be the j-th random number.

-- CONSTRUCTORS :

  procedure Create ( n : in natural32 );

  -- DESCRIPTION :
  --   Allocates memory space to contain a template for
  --   an n-dimensional polynomial system.

  procedure Add_Hyperplane ( i : in natural32; h : in Vector );

  -- DESCRIPTION :
  --   The hyperplane h is added to the i-the equation of the
  --   random product system.

  -- REQUIRED :                               __n_
  --   i <= n                                 \
  --   h : Vector(0..n) representing  h(0) +   >    h(j) x
  --                                          /___        j
  --                                           j=1

  procedure Change_Hyperplane ( i,j : in natural32; h : in Vector );

  -- DESCRIPTION :
  --   The (i,j)-th hyperplane will be changed into h.

-- SELECTORS :

  function Number_of_Hyperplanes ( i : natural32 ) return natural32;

  -- DESCRIPTION :
  --   returns the number of added hyperplanes for the i-th equation

  procedure Get_Hyperplane ( i,j : in natural32; h : out Vector );

  -- DESCRIPTION :
  --   returns the j-th hyperplane h for the i-th equation

  procedure Polynomial_System ( n,nbfree : in natural32 );

  -- DESCRIPTION :
  --   Based on the template, an n-dimensional random product
  --   polynomial system will be generated.
  --   The parameter nbfree indicates the number of free coefficients.
  --   After calling this routine, the package Random_Product_System
  --   will contain the data for a polynomial system.

  function Verify ( n : natural32; lp : List ) return natural32;

  -- DESCRIPTION :
  --   Computes the number of finite nonsingular solutions 
  --   of the final symmetric polynomial system.
  --   The list of positions lp indicates where the acceptable
  --   classes in the structure can be found.
  --   The structure is degenerate if this number does not
  --   correspond with the generalized Bezout number.

-- DESTRUCTOR :

  procedure Clear;

  -- DESCRIPTION :
  --   This procedure frees all memory space used by the template.

end Templates;
