with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Complex_Numbers;          use Standard_Complex_Numbers;
with Standard_Integer_Vectors;
with Standard_Integer_VecVecs;
with Standard_Floating_Vectors;
with Standard_Floating_VecVecs;
with Standard_Complex_Vectors;
with Standard_Complex_VecVecs;
with Standard_Complex_Laurentials;
with Standard_Complex_Laur_Systems;

package Double_Puiseux_Structures is

-- DESCRIPTION :
--   In extracting the data vectors from the containers, 
--   the coefficient vectors are reindexed to start from zero,
--   for the constant coefficients of the series coefficients.
--   Other functionality provided is the normalization of
--   binomial and product homotopies into the (hcf, hct, hdg) format
--   used by the Newton steps procedures.

  procedure Indexing_Series 
              ( powers : out Standard_Floating_VecVecs.Link_to_Array_of_VecVecs;
                coeffs : out Standard_Complex_VecVecs.Link_to_Array_of_VecVecs;
                vrblvl : in integer32 := 0);

  -- DESCRIPTION :
  --   Copies the arrays of vectors of vectors stored in the containers
  --   into powers and coefficients, adjusting the start of the indices
  --   of the coefficient vectors and shifting the powers.

  procedure Show_Data ( vrblvl : in integer32 := 0);

  -- DESCRIPTION :
  --   Retrieves the Laurent polynomial system and the corresponding
  --   coefficients of the series coefficients and the arrays of vectors.
  --   Writes the series system to screen as a test.

  function Is_Zero ( v : Standard_Integer_Vectors.Link_to_Vector )
                   return boolean;

  -- DESCRIPTION :
  --   Returns true if all elements in v are zero.

  function Is_Variable ( deg : Standard_Integer_Vectors.Link_to_Vector;
                         var : integer32 ) return boolean;

  -- DESCRIPTION :
  --   Returns true if the degrees in deg correspond to the variable
  --   with index given in var.

  function Is_Diagonal_Binomial_Homotopy
              ( p : Standard_Complex_Laur_Systems.Link_to_Laur_Sys;
                vrblvl : integer32 := 0 ) return boolean;

  -- DESCRIPTION :
  --   A binomial contains in its k-th equation the monomials x(k)
  --   and a constant.  This function runs through the monomials of p
  --   and checks if this condition is satified.

  function Extract_Leading_Coefficients
             ( cffs : Standard_Complex_VecVecs.Link_to_Array_of_VecVecs;
               hdg : Standard_Integer_VecVecs.Array_of_VecVecs;
               vrblvl : integer32 := 0 )
             return Standard_Complex_VecVecs.VecVec;

  -- DESCRIPTION :
  --   Given properly indexed coefficients of a Laurent homotopy,
  --   extracts the coefficients at index 1.
  --   Although this is not correct when the monomial is a variable,
  --   it gets fixed when normalizing the binomial homotopy.
  --   If the corresponding degrees in hdg are zero,
  --   that is: we have a constant, then the power should be zero.

  function Extract_Leading_Powers
             ( pwrs : Standard_Floating_VecVecs.Link_to_Array_of_VecVecs;
               hdg : Standard_Integer_VecVecs.Array_of_VecVecs;
               vrblvl : integer32 := 0 )
             return Standard_Floating_VecVecs.VecVec;

  -- DESCRIPTION :
  --   Given properly indexed powers of a Laurent homotopy,
  --   extracts the powers at index 1.
  --   Although this is not correct when the monomial is a variable,
  --   it gets fixed when normalizing the binomial homotopy.
  --   If the corresponding degrees in hdg are zero,
  --   that is: we have a constant, then the power should be zero.

  procedure Normalize_Binomial_Homotopy 
              ( hdg : in out Standard_Integer_VecVecs.Array_of_VecVecs;
                hcf : in out Standard_Complex_VecVecs.VecVec;
                hct : in out Standard_Floating_VecVecs.VecVec;
                vrblvl : in integer32 := 0 );

  -- DESCRIPTION :
  --   In a binomial homotopy, the first term is a variable,
  --   with index equal to the index of the polynomial,
  --   and the second term is the constant term, the solution series.

  function Truncation_Index
              ( pwrs : Standard_Floating_VecVecs.Link_to_VecVec )
              return integer32;

  -- DESCRIPTION :
  --   Given the powers of the series coefficients,
  --   returns the index of the last power,
  --   which corresponds to the truncation index,
  --   or the equivalent of Taylor series degree.

  function Is_Zero ( c : Complex_Number ) return boolean;

  -- DESCRIPTION :
  --   Returns true if both real and imaginary part of c are zero.

  procedure Product_Monomials
              ( p : in Standard_Complex_Laurentials.Poly;
                cffs : in Standard_Complex_VecVecs.Link_to_VecVec;
                pwrs : in Standard_Floating_VecVecs.Link_to_VecVec;
                hdg : out Standard_Integer_VecVecs.Link_to_VecVec;
                hcf : out Standard_Complex_Vectors.Link_to_Vector;
                hct : out Standard_Floating_Vectors.Link_to_Vector;
                vrblvl : in integer32 := 0 );

  -- DESCRIPTION :
  --   Given in p are the monomials of a product homotopy,
  --   where the coefficients of product are series,
  --   defines in hdg, hcf, hct the degrees, coefficients and powers
  --   of the expanded monomial version.

  procedure Normalize_Product_Homotopy
              ( p : in Standard_Complex_Laur_Systems.Link_to_Laur_Sys;
                cffs : in Standard_Complex_VecVecs.Link_to_Array_of_VecVecs;
                pwrs : in Standard_Floating_VecVecs.Link_to_Array_of_VecVecs;
                hdg : out Standard_Integer_VecVecs.Array_of_VecVecs;
                hcf : out Standard_Complex_VecVecs.VecVec;
                hct : out Standard_Floating_VecVecs.VecVec;
                vrblvl : in integer32 := 0 );

  -- DESCRIPTION :
  --   Extracts the degrees, coefficients and powers of t
  --   of a product Laurent homotopy with monomials given in p,
  --   coefficients of the series in cffs and powers in pwrs.
  --   The products of the polynomials are expanded,
  --   but the coefficients are still power series.

end Double_Puiseux_Structures;
