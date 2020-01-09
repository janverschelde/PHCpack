with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Vectors;
with Standard_Complex_Vectors;
with Standard_Complex_VecVecs;
with DoblDobl_Complex_Vectors;
with DoblDobl_Complex_VecVecs;
with QuadDobl_Complex_Vectors;
with QuadDobl_Complex_VecVecs;
with Standard_Speelpenning_Convolutions;
with DoblDobl_Speelpenning_Convolutions;
with QuadDobl_Speelpenning_Convolutions;

package Newton_Convolutions is

-- DESCRIPTION :
--   Newton's method on power series computed with convolution circuits

  function Series_Coefficients
             ( v : Standard_Complex_Vectors.Vector;
               d : integer32 )
             return Standard_Complex_VecVecs.VecVec;
  function Series_Coefficients
             ( v : DoblDobl_Complex_Vectors.Vector;
               d : integer32 )
             return DoblDobl_Complex_VecVecs.VecVec;
  function Series_Coefficients
             ( v : QuadDobl_Complex_Vectors.Vector;
               d : integer32 )
             return QuadDobl_Complex_VecVecs.VecVec;

  -- DESCRIPTION :
  --   Returns the coefficients of a vector of series,
  --   with leading coefficients the components of a solution v.
  --   The series are of degree d.

  procedure Minus ( v : in Standard_Complex_VecVecs.VecVec );
  procedure Minus ( v : in DoblDobl_Complex_VecVecs.VecVec );
  procedure Minus ( v : in QuadDobl_Complex_VecVecs.VecVec );

  -- DESCRIPTION :
  --   Flips the sign of all coefficients in v.

  procedure Update ( x,y : in Standard_Complex_VecVecs.VecVec );
  procedure Update ( x,y : in DoblDobl_Complex_VecVecs.VecVec );
  procedure Update ( x,y : in QuadDobl_Complex_VecVecs.VecVec );

  -- DESCRIPTION :
  --   Adds to all coefficients of x the corresponding coefficient of y.

  -- REQUIRED : x'range = y'range, and for all k in x'range
  --   x(k)'range = y(k)'range.

  procedure Standard_Newton_Step
              ( s : in Standard_Speelpenning_Convolutions.Link_to_System;
                scf : in Standard_Complex_VecVecs.VecVec;
                info : out integer32;
                ipvt : in out Standard_Integer_Vectors.Vector;
                wrk : in Standard_Complex_Vectors.Link_to_Vector );
  procedure DoblDobl_Newton_Step
              ( s : in DoblDobl_Speelpenning_Convolutions.Link_to_System;
                scf : in DoblDobl_Complex_VecVecs.VecVec;
                info : out integer32;
                ipvt : in out Standard_Integer_Vectors.Vector;
                wrk : in DoblDobl_Complex_Vectors.Link_to_Vector );
  procedure QuadDobl_Newton_Step
              ( s : in QuadDobl_Speelpenning_Convolutions.Link_to_System;
                scf : in QuadDobl_Complex_VecVecs.VecVec;
                info : out integer32;
                ipvt : in out Standard_Integer_Vectors.Vector;
                wrk : in QuadDobl_Complex_Vectors.Link_to_Vector );

  -- DESCRIPTION :
  --   Applies one Newton step on the convolution circuits c,
  --   departing from the series coefficients in s,
  --   in double, double double, or quad double precision.

  -- ON ENTRY :
  --   s        system of convolution circuits;
  --   scf      vector of coefficients of power series;
  --   ipvt     vector for the pivoting information in the LU factorization;
  --   wrk      work space for the matrix series solver.

  -- ON RETURN :
  --   info     info from the LU factorization;
  --   ipvt     pivoting of the LU factorization on the lead matrix.

end Newton_Convolutions;
