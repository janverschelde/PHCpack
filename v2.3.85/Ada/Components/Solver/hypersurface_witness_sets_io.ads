with text_io;                           use text_io;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Floating_Numbers;         use Standard_Floating_Numbers;
with Standard_Complex_Vectors;          use Standard_Complex_Vectors;
with Standard_Complex_Polynomials;      use Standard_Complex_Polynomials;
with Standard_Complex_Poly_Systems;     use Standard_Complex_Poly_Systems;
with Standard_Complex_Solutions;        use Standard_Complex_Solutions;

package Hypersurface_Witness_Sets_io is

-- DESCRIPTION :
--   This package collects the function to compute and write a witness set
--   for one polynomial in several variables.

  function Embedded_System 
               ( n : in integer32; p : in Poly; b,v,t : in Vector )
               return Poly_Sys;

  -- DESCRIPTION :
  --   Returns the embedded system in the extrinsic format,
  --   converting the intrinsic line x = b + t*v.

  procedure Write_Witness_Set
               ( n : in integer32; p : in Poly; b,v,t : in Vector );
  procedure Write_Witness_Set
               ( file : in file_type;
                 n : in integer32; p : in Poly; b,v,t : in Vector );

  -- DESCRIPTION :
  --   Writes the witness set for p given intrinsically as b + t*v
  --   in extrinsic format to file.
  --   Prompts the user to provide a file name if file is not provided.

  procedure Call_Root_Finder ( p : in Poly );

  -- DESCRIPTION :
  --   Interactive driver to compute a witness set for p
  --   and to write the witness set to file.
  --   The user is prompted for information.

  procedure Call_Root_Finder
               ( file : in file_type; p : in Poly; output : in boolean;
                 eps : in double_float; fail : out boolean );

  -- DESCRIPTION :
  --   Computes a witness set for p and writes it to file.

  -- ON ENTRY :
  --   file      to write the results on;
  --   p         a polynomial in several variables;
  --   output    set to true if intermediate output is desired;
  --   eps       accuracy requirement for univariate root finder.

  -- ON RETURN :
  --   fail      true if accuracy is not met, false otherwise.

  procedure Silent_Root_Finder
               ( p : in Poly; eps : in double_float; fail : out boolean;
                 e : out Link_to_Poly_Sys; esols : out Solution_List );

  -- DESCRIPTION :
  --   Computes a witness set for p and writes it to file.

  -- ON ENTRY :
  --   p         a polynomial in several variables;
  --   eps       accuracy requirement for univariate root finder.

  -- ON RETURN :
  --   fail      true if accuracy is not met, false otherwise;
  --   e         embedded polynomial system;
  --   esols     witness points.

end Hypersurface_Witness_Sets_io;
