with text_io;                          use text_io;
with Standard_Natural_Numbers;         use Standard_Natural_Numbers;
with Standard_Natural_Vectors;
with Standard_Natural_VecVecs;
with Standard_Complex_Poly_Systems;
with Standard_Complex_Laur_Systems;
with DoblDobl_Complex_Poly_Systems;
with DoblDobl_Complex_Laur_Systems;
with QuadDobl_Complex_Poly_Systems;
with QuadDobl_Complex_Laur_Systems;
with Standard_Complex_Solutions;
with DoblDobl_Complex_Solutions;
with QuadDobl_Complex_Solutions;
with Sample_Point_Lists;               use Sample_Point_Lists;
with DoblDobl_Sample_Lists;            use DoblDobl_Sample_Lists;
with QuadDobl_Sample_Lists;            use QuadDobl_Sample_Lists;

package Drivers_to_Breakup_Solutions is

-- DESCRIPTION :
--   The drivers in this package prompt the user for 
--   a witness set representation for a pure dimensional solution set,
--   ask to select the method (monodromy or enumeration) and then
--   decompose the solution sets into irreducible factors.

  function Create ( file : file_type;
                    p : Standard_Complex_Poly_Systems.Poly_Sys;
                    sols : Standard_Complex_Solutions.Solution_List;
                    dim : natural32 )
                  return Array_of_Standard_Sample_Lists;
  function Create ( file : file_type;
                    p : DoblDobl_Complex_Poly_Systems.Poly_Sys;
                    sols : DoblDobl_Complex_Solutions.Solution_List;
                    dim : natural32 )
                  return Array_of_DoblDobl_Sample_Lists;
  function Create ( file : file_type;
                    p : QuadDobl_Complex_Poly_Systems.Poly_Sys;
                    sols : QuadDobl_Complex_Solutions.Solution_List;
                    dim : natural32 )
                  return Array_of_QuadDobl_Sample_Lists;

  -- DESCRIPTION :
  --   Returns a grid of sample points needed for linear traces
  --   in the combinatorial enumeration factorization method,
  --   for a witness set defined by an ordinary polynomial system,
  --   in standard double, double double, and quad double precision.

  -- ON ENTRY :
  --   file     file which must be opened for output;
  --   ep       embedded ordinary polynomial system of the witness set;
  --   sols     generic points in the witness set;
  --   dim      dimension of the solution set.

  function Create ( file : file_type;
                    p : Standard_Complex_Laur_Systems.Laur_Sys;
                    sols : Standard_Complex_Solutions.Solution_List;
                    dim : natural32 )
                  return Array_of_Standard_Sample_Lists;
  function Create ( file : file_type;
                    p : DoblDobl_Complex_Laur_Systems.Laur_Sys;
                    sols : DoblDobl_Complex_Solutions.Solution_List;
                    dim : natural32 )
                  return Array_of_DoblDobl_Sample_Lists;
  function Create ( file : file_type;
                    p : QuadDobl_Complex_Laur_Systems.Laur_Sys;
                    sols : QuadDobl_Complex_Solutions.Solution_List;
                    dim : natural32 )
                  return Array_of_QuadDobl_Sample_Lists;

  -- DESCRIPTION :
  --   Returns a grid of sample points needed for linear traces,
  --   in the combinatorial enumeration factorization method,
  --   for a witness set defined by a Laurent polynomial system,
  --   in standard double, double double, and quad double precision.

  -- ON ENTRY :
  --   file     file which must be opened for output;
  --   ep       embedded Laurent polynomial system of the witness set;
  --   sols     generic points in the witness set;
  --   dim      dimension of the solution set.

  procedure Standard_Enumerate_Decomposition
              ( file : in file_type;
                ep : in Standard_Complex_Poly_Systems.Poly_Sys;
                sols : in Standard_Complex_Solutions.Solution_List;
                dim : in natural32 );
  procedure DoblDobl_Enumerate_Decomposition
              ( file : in file_type;
                ep : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                sols : in DoblDobl_Complex_Solutions.Solution_List;
                dim : in natural32 );
  procedure QuadDobl_Enumerate_Decomposition
              ( file : in file_type;
                ep : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                sols : in QuadDobl_Complex_Solutions.Solution_List;
                dim : in natural32 );

  -- DESCRIPTION :
  --   Factors a component using linear traces combinatorially,
  --   in standard double, double double, or quad double precision,
  --   for a witness set defined by an ordinary polynomial system.
  --   All output is written to the single output file.

  -- ON ENTRY :
  --   file     file which must be opened for output;
  --   ep       embedded polynomial system of the witness set;
  --   sols     generic points in the witness set;
  --   dim      dimension of the solution set.

  procedure Standard_Enumerate_Decomposition
              ( file : in file_type;
                ep : in Standard_Complex_Laur_Systems.Laur_Sys;
                sols : in Standard_Complex_Solutions.Solution_List;
                dim : in natural32 );
  procedure DoblDobl_Enumerate_Decomposition
              ( file : in file_type;
                ep : in DoblDobl_Complex_Laur_Systems.Laur_Sys;
                sols : in DoblDobl_Complex_Solutions.Solution_List;
                dim : in natural32 );
  procedure QuadDobl_Enumerate_Decomposition
              ( file : in file_type;
                ep : in QuadDobl_Complex_Laur_Systems.Laur_Sys;
                sols : in QuadDobl_Complex_Solutions.Solution_List;
                dim : in natural32 );

  -- DESCRIPTION :
  --   Factors a component using linear traces combinatorially,
  --   in standard double, double double, or quad double precision,
  --   for a witness set defined by a Laurent polynomial system.
  --   All output is written to the single output file.

  -- ON ENTRY :
  --   file     file which must be opened for output;
  --   ep       embedded polynomial system of the witness set;
  --   sols     generic points in the witness set;
  --   dim      dimension of the solution set.

  function Append_fk ( name : string; k : natural32 ) return string;

  -- DESCRIPTION :
  --   Returns the name_fk, where the k is the string representation
  --   of the given natural number k.

  function Select_Witness_Set_for_Factor
             ( sols : in Standard_Complex_Solutions.Solution_List;
               f : in Standard_Natural_Vectors.Vector )
             return Standard_Complex_Solutions.Solution_List;
  function Select_Witness_Set_for_Factor
             ( sols : in DoblDobl_Complex_Solutions.Solution_List;
               f : in Standard_Natural_Vectors.Vector )
             return DoblDobl_Complex_Solutions.Solution_List;
  function Select_Witness_Set_for_Factor
             ( sols : in QuadDobl_Complex_Solutions.Solution_List;
               f : in Standard_Natural_Vectors.Vector )
             return QuadDobl_Complex_Solutions.Solution_List;

  -- DESCRIPTION :
  --   Returns those solutions whose index occurs in f.

  procedure Write_Witness_Sets_for_Factors
              ( name : in string; 
                p : in Standard_Complex_Poly_Systems.Poly_Sys;
                sols : in Standard_Complex_Solutions.Solution_List;
                dim : in natural32; f : in Standard_Natural_VecVecs.VecVec );
  procedure Write_Witness_Sets_for_Factors
              ( name : in string; 
                p : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                sols : in DoblDobl_Complex_Solutions.Solution_List;
                dim : in natural32; f : in Standard_Natural_VecVecs.VecVec );
  procedure Write_Witness_Sets_for_Factors
              ( name : in string; 
                p : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                sols : in QuadDobl_Complex_Solutions.Solution_List;
                dim : in natural32; f : in Standard_Natural_VecVecs.VecVec );

  -- DESCRIPTION  :
  --   Writes a witness set for every irreducible factor of the
  --   embedded polynomial system p with corresponding solutions in sols.
  --   The factorization is given in f.

  procedure Write_Witness_Sets_for_Factors
              ( name : in string; 
                p : in Standard_Complex_Laur_Systems.Laur_Sys;
                sols : in Standard_Complex_Solutions.Solution_List;
                dim : in natural32; f : in Standard_Natural_VecVecs.VecVec );
  procedure Write_Witness_Sets_for_Factors
              ( name : in string; 
                p : in DoblDobl_Complex_Laur_Systems.Laur_Sys;
                sols : in DoblDobl_Complex_Solutions.Solution_List;
                dim : in natural32; f : in Standard_Natural_VecVecs.VecVec );
  procedure Write_Witness_Sets_for_Factors
              ( name : in string; 
                p : in QuadDobl_Complex_Laur_Systems.Laur_Sys;
                sols : in QuadDobl_Complex_Solutions.Solution_List;
                dim : in natural32; f : in Standard_Natural_VecVecs.VecVec );

  -- DESCRIPTION  :
  --   Writes a witness set for every irreducible factor of the
  --   embedded Laurent system p with corresponding solutions in sols.
  --   The factorization is given in f.

  procedure Standard_Enumerate_Decomposition
              ( file : in file_type; name : in string;
                ep : in Standard_Complex_Poly_Systems.Poly_Sys;
                sols : in Standard_Complex_Solutions.Solution_List;
                dim : in natural32 );
  procedure DoblDobl_Enumerate_Decomposition
              ( file : in file_type; name : in string;
                ep : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                sols : in DoblDobl_Complex_Solutions.Solution_List;
                dim : in natural32 );
  procedure QuadDobl_Enumerate_Decomposition
              ( file : in file_type; name : in string;
                ep : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                sols : in QuadDobl_Complex_Solutions.Solution_List;
                dim : in natural32 );

  -- DESCRIPTION :
  --   Factors a component using linear traces combinatorially,
  --   in standard double, double double, or quad double precision,
  --   for a witness set defined by an ordinary polynomial system.
  --   Each irreducible factor is written as a witness set to a 
  --   separate file, with the same prefix name.

  -- ON ENTRY :
  --   file     file which must be opened for output;
  --   name     prefix for the output file names for the factors;
  --   ep       embedded polynomial system of the witness set;
  --   sols     generic points in the witness set;
  --   dim      dimension of the solution set.

  procedure Standard_Enumerate_Decomposition
              ( file : in file_type; name : in string;
                ep : in Standard_Complex_Laur_Systems.Laur_Sys;
                sols : in Standard_Complex_Solutions.Solution_List;
                dim : in natural32 );
  procedure DoblDobl_Enumerate_Decomposition
              ( file : in file_type; name : in string;
                ep : in DoblDobl_Complex_Laur_Systems.Laur_Sys;
                sols : in DoblDobl_Complex_Solutions.Solution_List;
                dim : in natural32 );
  procedure QuadDobl_Enumerate_Decomposition
              ( file : in file_type; name : in string;
                ep : in QuadDobl_Complex_Laur_Systems.Laur_Sys;
                sols : in QuadDobl_Complex_Solutions.Solution_List;
                dim : in natural32 );

  -- DESCRIPTION :
  --   Factors a component using linear traces combinatorially,
  --   in standard double, double double, or quad double precision,
  --   for a witness set defined by a Laurent system.
  --   Each irreducible factor is written as a witness set to a 
  --   separate file, with the same prefix name.

  -- ON ENTRY :
  --   file     file which must be opened for output;
  --   name     prefix for the output file names for the factors;
  --   ep       embedded polynomial system of the witness set;
  --   sols     generic points in the witness set;
  --   dim      dimension of the solution set.

  procedure Standard_Monodromy_Decomposition
              ( file : in file_type; name : in string;
                ep : in Standard_Complex_Poly_Systems.Poly_Sys;
                sols : in Standard_Complex_Solutions.Solution_List;
                dim : in natural32 );
  procedure DoblDobl_Monodromy_Decomposition
              ( file : in file_type; name : in string;
                ep : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                sols : in DoblDobl_Complex_Solutions.Solution_List;
                dim : in natural32 );
  procedure QuadDobl_Monodromy_Decomposition
              ( file : in file_type; name : in string;
                ep : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                sols : in QuadDobl_Complex_Solutions.Solution_List;
                dim : in natural32 );

  -- DESCRIPTION :
  --   Factors a pure dimensional solution set using monodromy loops
  --   in standard double, double double, or quad double precision.

  -- ON ENTRY :
  --   file     file which must be opened for output;
  --   name     prefix for the output file names for the factors;
  --   ep       embedded polynomial system of the witness set;
  --   sols     generic points in the witness set;
  --   dim      dimension of the solution set.

  function Prompt_for_Method return character;

  -- DESCRIPTION :
  --   Displays the menu with the choice between monodromy
  --   and combinatorial enumeration.
  --   Returns '1' for monodromy and '2' for enumeration.

  procedure Standard_Breakup
              ( file : in file_type; name : in string;
                p : in Standard_Complex_Poly_Systems.Poly_Sys;
                sols : in Standard_Complex_Solutions.Solution_List;
                dim : in natural32 );
  procedure DoblDobl_Breakup
              ( file : in file_type; name : in string;
                p : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                sols : in DoblDobl_Complex_Solutions.Solution_List;
                dim : in natural32 );
  procedure QuadDobl_Breakup
              ( file : in file_type; name : in string;
                p : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                sols : in QuadDobl_Complex_Solutions.Solution_List;
                dim : in natural32 );

  -- DESCRIPTION :
  --   Prompts the user for method (combinatorial or monodromy)
  --   and then factors the given witness set,
  --   defined by an ordinary polynomial system,
  --   in standard double, double double, or quad double precision.

  -- ON ENTRY :
  --   file     file which must be opened for output;
  --   name     prefix for the output file names for the factors;
  --   p        embedded polynomial system of the witness set;
  --   sols     generic points in the witness set;
  --   dim      dimension of the solution set.

  procedure Standard_Breakup
              ( file : in file_type; name : in string;
                p : in Standard_Complex_Laur_Systems.Laur_Sys;
                sols : in Standard_Complex_Solutions.Solution_List;
                dim : in natural32 );
  procedure DoblDobl_Breakup
              ( file : in file_type; name : in string;
                p : in DoblDobl_Complex_Laur_Systems.Laur_Sys;
                sols : in DoblDobl_Complex_Solutions.Solution_List;
                dim : in natural32 );
  procedure QuadDobl_Breakup
              ( file : in file_type; name : in string;
                p : in QuadDobl_Complex_Laur_Systems.Laur_Sys;
                sols : in QuadDobl_Complex_Solutions.Solution_List;
                dim : in natural32 );

  -- DESCRIPTION :
  --   Prompts the user for method (combinatorial or monodromy)
  --   and then factors the given witness set,
  --   defined by a Laurent polynomial system,
  --   in standard double, double double, or quad double precision.

  -- ON ENTRY :
  --   file     file which must be opened for output;
  --   name     prefix for the output file names for the factors;
  --   p        embedded polynomial system of the witness set;
  --   sols     generic points in the witness set;
  --   dim      dimension of the solution set.

  procedure Standard_Breakup ( infilename,outfilename : in string );
  procedure DoblDobl_Breakup ( infilename,outfilename : in string );
  procedure QuadDobl_Breakup ( infilename,outfilename : in string );

  -- DESCRIPTION :
  --   Reads a witness set from the file with name infilename if not empty,
  --   otherwise prompts the user for a filtered witness set a computes
  --   its numerical irreducible decomposition,
  --   in standard double, double double, or quad double precision.
  --   Reads the embedded polynomial system, display the factorization menu
  --   and either does monodromy loops or enumerates factors.  

end Drivers_to_Breakup_Solutions;
