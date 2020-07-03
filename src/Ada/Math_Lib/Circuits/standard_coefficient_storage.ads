with Standard_Floating_VecVecs;
with Standard_Floating_VecVecVecs;
with Standard_Coefficient_Convolutions;  use Standard_Coefficient_Convolutions;

package Standard_Coefficient_Storage is

-- DESCRIPTION :
--   The coefficients of the power series of circuits are shifted
--   during the tracking on one path.  This package allocates space,
--   stores, and restores the coefficients of the circuits.

  procedure Allocate_and_Store
              ( c : in Link_to_Circuit;
                rcf : out Standard_Floating_VecVecs.Link_to_VecVec;
                icf : out Standard_Floating_VecVecs.Link_to_VecVec );

  -- DESCRIPTION :
  --   Allocates space for the real and imaginary parts of the coefficients
  --   of the circuit c and stores the coefficients in vectors of vectors.

  -- ON ENTRY :
  --   c        a convolution circuit with well defined coefficients.

  -- ON RETURN :
  --   rcf      real parts of the coefficients of the series in c,
  --            with the real part of the constant series at index 0;
  --   icf      imaginary parts of the coefficients of the series in c,
  --            with the imaginary part of the constant series at index 0.

  procedure Allocate_and_Store
              ( c : in Circuits;
                rcf : out Standard_Floating_VecVecVecs.Link_to_VecVecVec;
                icf : out Standard_Floating_VecVecVecs.Link_to_VecVecVec );

  -- DESCRIPTION :
  --   Allocates space for the real and imaginary parts of the coefficients
  --   of the circuits in c and stores all coefficients of the circuits
  --   in c into vectors of vectors of vectors.

  -- ON ENTRY :
  --   c        convolution circuits with well defined coefficients.

  -- ON RETURN :
  --   rcf      vector of vectors of vectors, where rcf(k) contains
  --            the real parts of the coefficients of the series in c(k),
  --            with the real part of the constant series at index 0,
  --            for all k in c'range;
  --   icf      vector of vectors of vectors, where icf(k) contains
  --            the imaginary parts of the coefficients of the series in c(k),
  --            with the imaginary part of the constant series at index 0,
  --            for all k in c'range.

  procedure Restore
              ( rcf : in Standard_Floating_VecVecs.Link_to_VecVec;
                icf : in Standard_Floating_VecVecs.Link_to_VecVec;
                c : in Link_to_Circuit );

  -- DESCRIPTION :
  --   Given allocated space for all coefficients of the circuit c,
  --   restores the values to their original values.

  -- REQUIRED : rcf'range = icf'range = 0..c.nbr,
  --   and for all k in rcf'range, rcf(k) and icf(k) is allocated.

  -- ON ENTRY :
  --   rcf      real parts of the coefficients of series for c,
  --            with the real part of the constant series at index 0;
  --   icf      imaginary parts of the coefficients of series for c,
  --            with the imaginary part of the constant series at index 0.

  -- ON RETURN :
  --   c        a convolution circuit with coefficients restored
  --            with the values in rcf and icf.

  procedure Restore
              ( rcf : in Standard_Floating_VecVecVecs.Link_to_VecVecVec;
                icf : in Standard_Floating_VecVecVecs.Link_to_VecVecVec;
                c : in Circuits );

  -- DESCRIPTION :
  --   Restores all coefficients in c to the values
  --   stored in the rcf and icf.

end Standard_Coefficient_Storage;
