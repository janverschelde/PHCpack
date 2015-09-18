with text_io;                           use text_io;
with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Natural_Vectors;          
with Standard_Natural_VecVecs;
with Standard_Complex_Vectors;
with DoblDobl_Complex_Vectors;
with QuadDobl_Complex_Vectors;
with Standard_Complex_Polynomials;
with Standard_Complex_Poly_Systems;
with DoblDobl_Complex_Polynomials;
with DoblDobl_Complex_Poly_Systems;
with QuadDobl_Complex_Polynomials;
with QuadDobl_Complex_Poly_Systems;

package Drivers_to_Factor_Polynomials is

-- DESCRIPTION :
--   This package provides the main drivers to factor polynomials
--   in several variables with complex coefficients.

  procedure Factor ( monodromy : in boolean; n : in natural32;
                     p : in Standard_Complex_Polynomials.Poly;
                     b,v : out Standard_Complex_Vectors.Vector;
                     wp : out Standard_Complex_Vectors.Link_to_Vector;
                     mw,mf : out Standard_Natural_Vectors.Link_to_Vector;
                     deco : out Standard_Natural_VecVecs.Link_to_VecVec;
                     fail : out boolean;
                     factors
                       : out Standard_Complex_Poly_Systems.Link_to_Poly_Sys );
  procedure Factor ( monodromy : in boolean; n : in natural32;
                     p : in DoblDobl_Complex_Polynomials.Poly;
                     b,v : out DoblDobl_Complex_Vectors.Vector;
                     wp : out DoblDobl_Complex_Vectors.Link_to_Vector;
                     mw,mf : out Standard_Natural_Vectors.Link_to_Vector;
                     deco : out Standard_Natural_VecVecs.Link_to_VecVec;
                     fail : out boolean;
                     factors
                       : out DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys );
  procedure Factor ( monodromy : in boolean; n : in natural32;
                     p : in QuadDobl_Complex_Polynomials.Poly;
                     b,v : out QuadDobl_Complex_Vectors.Vector;
                     wp : out QuadDobl_Complex_Vectors.Link_to_Vector;
                     mw,mf : out Standard_Natural_Vectors.Link_to_Vector;
                     deco : out Standard_Natural_VecVecs.Link_to_VecVec;
                     fail : out boolean;
                     factors
                       : out QuadDobl_Complex_Poly_Systems.Link_to_Poly_Sys );

  -- DESCRIPTION :
  --   Decomposes a polynomial into irreducible factor without output.

  -- ON ENTRY :
  --   monodromy     flag to indicate if monodromy needs to be used;
  --   n             number of variables in the polynomial;
  --   p             polynomial in n variables with complex coefficients.

  -- ON RETURN :
  --   b,v           represents the line b + t*v;
  --   wp            witness points on the line b + t*v;
  --   mw            mw(i) is multiplicity of i-th witness point wp(i);
  --   mf            mf(i) is multiplicity of i-th factor deco(i);
  --   deco          groups witness points according to irreducible factors;
  --   fail          true if failed to compute witness points,
  --                 or if certificate with traces failed;
  --   factors       irreducible factors of the polynomial.
  
  procedure Driver_to_Factor 
               ( monodromy : in boolean; n : in natural32;
                 p : in Standard_Complex_Polynomials.Poly;
                 fail : out boolean;
                 factors : out Standard_Complex_Poly_Systems.Link_to_Poly_Sys );
  procedure Driver_to_Factor 
               ( monodromy : in boolean; n : in natural32;
                 p : in DoblDobl_Complex_Polynomials.Poly;
                 fail : out boolean;
                 factors : out DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys );
  procedure Driver_to_Factor 
               ( monodromy : in boolean; n : in natural32;
                 p : in QuadDobl_Complex_Polynomials.Poly;
                 fail : out boolean;
                 factors : out QuadDobl_Complex_Poly_Systems.Link_to_Poly_Sys );
  procedure Driver_to_Factor 
               ( file : in file_type; output,monodromy : in boolean;
	         n : in natural32;
                 p : in Standard_Complex_Polynomials.Poly;
		 fail : out boolean;
                 factors : out Standard_Complex_Poly_Systems.Link_to_Poly_Sys );
  procedure Driver_to_Factor 
               ( file : in file_type; output,monodromy : in boolean;
	         n : in natural32;
                 p : in DoblDobl_Complex_Polynomials.Poly;
		 fail : out boolean;
                 factors : out DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys );
  procedure Driver_to_Factor 
               ( file : in file_type; output,monodromy : in boolean;
	         n : in natural32;
                 p : in QuadDobl_Complex_Polynomials.Poly;
		 fail : out boolean;
                 factors : out QuadDobl_Complex_Poly_Systems.Link_to_Poly_Sys );

  -- DESCRIPTION :
  --   Decomposes a polynomial into irreducible factors, writing
  --   intermediate output and diagnostics on screen or to a file.

  -- ON ENTRY :
  --   file      to write intermediate output and diagnostics;
  --   output    flag to require intermediate output during continuation;
  --   monodromy flag to indicate whether monodromy needs to be used;
  --   n         number of variables of the polynomial;
  --   p         polynomials to factor.

  -- ON RETURN :
  --   fail      true if failed to compute witness points;
  --   factors   factors of the polynomial.

  procedure Write_Factors
               ( filename : in string;
                 factors : in Standard_Complex_Poly_Systems.Poly_Sys );
  procedure Write_Factors
               ( filename : in string;
                 factors : in DoblDobl_Complex_Poly_Systems.Poly_Sys );
  procedure Write_Factors
               ( filename : in string;
                 factors : in QuadDobl_Complex_Poly_Systems.Poly_Sys );

  -- DESCRIPTION :
  --   Writes the factors to files which start with filename.

  procedure Driver_to_Factor_Polynomial;
  procedure Driver_to_Factor_Polynomial ( filename : in string );

  -- DESCRIPTION :
  --   Interactive driver to factor a polynomial in several variables,
  --   with complex coefficients.
  --   The second routine writes the factors to separate files,
  --   all starting with the same filename, but with different suffix.

end Drivers_to_Factor_Polynomials;
