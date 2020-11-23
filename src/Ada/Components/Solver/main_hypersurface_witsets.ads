with text_io;                           use text_io;
with String_Splitters;                  use String_Splitters;
with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Complex_Polynomials;
with Standard_Complex_Laurentials;
with DoblDobl_Complex_Polynomials;
with DoblDobl_Complex_Laurentials;
with QuadDobl_Complex_Polynomials;
with QuadDobl_Complex_Laurentials;

package Main_Hypersurface_Witsets is

-- DESCRIPTION :
--   Wraps the witness set representation for hypersurfaces.

  procedure Make_Witness_Set
              ( file : in file_type; name : in string;
                n : in natural32;
                p : in Standard_Complex_Polynomials.Poly );

  -- DESCRIPTION :
  --   Makes a witness set for the hypersurface defined by p,
  --   where p is an ordinary polynomial in several variables,
  --   with coefficients in double precision.

  -- ON ENTRY :
  --   file     file which must be opened for output;
  --   name     name of the input file, used as prefix
  --            for the witness set;
  --   n        number of variables in the polynomial p;
  --   p        an ordinary polynomial in n veriables.

  procedure Make_Witness_Set
              ( file : in file_type; name : in string;
                n : in natural32;
                p : in DoblDobl_Complex_Polynomials.Poly );

  -- DESCRIPTION :
  --   Generates a random offset and random direction vector
  --   and then computes all roots of a univariate problem
  --   to make a witness set for the hypersurface defined by p,
  --   in double double precision.

  -- ON ENTRY :
  --   file     file which must be opened for output;
  --   name     name of the input file, used as prefix
  --            for the witness set;
  --   n        number of variables in the polynomial p;
  --   p        an ordinary polynomial p in n veriables.

  procedure Make_Witness_Set
              ( file : in file_type; name : in string;
                n : in natural32;
                p : in QuadDobl_Complex_Polynomials.Poly );

  -- DESCRIPTION :
  --   Generates a random offset and random direction vector
  --   and then computes all roots of a univariate problem
  --   to make a witness set for the hypersurface defined by p,
  --   in quad double precision.

  -- ON ENTRY :
  --   file     file which must be opened for output;
  --   name     name of the input file, used as prefix
  --            for the witness set;
  --   n        number of variables in the polynomial p;
  --   p        an ordinary polynomial p in n veriables.

  procedure Make_Witness_Set
              ( file : in file_type; name : in string;
                n : in natural32;
                p : in Standard_Complex_Laurentials.Poly );

  -- DESCRIPTION :
  --   Makes a witness set for the hypersurface defined by p,
  --   where p is a Laurent polynomial,
  --   with coefficients in double precision.

  -- ON ENTRY :
  --   file     file which must be opened for output;
  --   name     name of the input file, used as prefix
  --            for the witness set;
  --   n        number of variables in the polynomial p;
  --   p        a Laurent polynomial in n veriables.

  procedure Make_Witness_Set
              ( file : in file_type; name : in string;
                n : in natural32;
                p : in DoblDobl_Complex_Laurentials.Poly );

  -- DESCRIPTION :
  --   Makes a witness set for the hypersurface defined by p,
  --   where p is a Laurent polynomial,
  --   with coefficients in double double precision.

  -- ON ENTRY :
  --   file     file which must be opened for output;
  --   name     name of the input file, used as prefix
  --            for the witness set;
  --   n        number of variables in the polynomial p;
  --   p        a Laurent polynomial in n veriables.

  procedure Make_Witness_Set
              ( file : in file_type; name : in string;
                n : in natural32;
                p : in QuadDobl_Complex_Laurentials.Poly );

  -- DESCRIPTION :
  --   Makes a witness set for the hypersurface defined by p,
  --   where p is a Laurent polynomial,
  --   with coefficients in quad double precision.

  -- ON ENTRY :
  --   file     file which must be opened for output;
  --   name     name of the input file, used as prefix
  --            for the witness set;
  --   n        number of variables in the polynomial p;
  --   p        a Laurent polynomial in n veriables.

  procedure Read_Input_Polynomial
              ( polyfile : in string; name : out Link_to_String;
                n : out natural32;
                p : out Standard_Complex_Laurentials.Poly );

  -- DESCRIPTION :
  --   Read the polynomial from polyfile, or 
  --   if polyfile is the empty string, then prompts the user
  --   for a file name to read the polynomial from.

  -- ON ENTRY :
  --   polyfile is the name of the file for the polynomial.

  -- ON RETURN :
  --   name     name of the input file;
  --   n        number of variables in the polynomial;
  --   p        the polynomial read.

  procedure Read_Input_Polynomial
              ( polyfile : in string; name : out Link_to_String;
                n : out natural32;
                p : out DoblDobl_Complex_Laurentials.Poly );

  -- DESCRIPTION :
  --   Read the polynomial from polyfile, or 
  --   if polyfile is the empty string, then prompts the user
  --   for a file name to read the polynomial from.

  -- ON ENTRY :
  --   polyfile is the name of the file for the polynomial.

  -- ON RETURN :
  --   name     name of the input file;
  --   n        number of variables in the polynomial;
  --   p        the polynomial read.

  procedure Read_Input_Polynomial
              ( polyfile : in string; name : out Link_to_String;
                n : out natural32;
                p : out QuadDobl_Complex_Laurentials.Poly );

  -- DESCRIPTION :
  --   Read the polynomial from polyfile, or 
  --   if polyfile is the empty string, then prompts the user
  --   for a file name to read the polynomial from.

  -- ON ENTRY :
  --   polyfile is the name of the file for the polynomial.

  -- ON RETURN :
  --   name     name of the input file;
  --   n        number of variables in the polynomial;
  --   p        the polynomial read.

  procedure Standard_Main
              ( polyfile,logfile : in string; vrblvl : in integer32 := 0 );

  -- DESCRIPTION :
  --   Constructs a witness set of a hypersurface in double precision,
  --   called by phc -l

  -- ON ENTRY :
  --   polyfile     file with one (Laurent) polynomial in several variables;
  --   logfile      file name to write diagnostics on;
  --   vrblvl       the verbose level.

  -- ON RETURN :
  --   polyfile_wk, a file with a witness set for the hypersurface,
  --   with k the dimension of the hypersurface.

  procedure DoblDobl_Main
              ( polyfile,logfile : in string; vrblvl : in integer32 := 0 );

  -- DESCRIPTION :
  --   Constructs a witness set of a hypersurface in double double precision,
  --   called by phc -l2.

  -- ON ENTRY :
  --   polyfile     file with one (Laurent) polynomial in several variables;
  --   logfile      file name to write diagnostics on;
  --   vrblvl       the verbose level.

  -- ON RETURN :
  --   polyfile_wk, a file with a witness set for the hypersurface,
  --   with k the dimension of the hypersurface.

  procedure QuadDobl_Main
              ( polyfile,logfile : in string; vrblvl : in integer32 := 0 );

  -- DESCRIPTION :
  --   Constructs a witness set of a hypersurface in quad double precision,
  --   called by phc -l4.

  -- ON ENTRY :
  --   polyfile     file with one (Laurent) polynomial in several variables;
  --   logfile      file name to write diagnostics on;
  --   vrblvl       the verbose level.

  -- ON RETURN :
  --   polyfile_wk, a file with a witness set for the hypersurface,
  --   with k the dimension of the hypersurface.

end Main_Hypersurface_Witsets;
