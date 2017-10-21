with text_io;                           use text_io;
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

package Monodromy_Homotopies_io is

-- DESCRIPTION :
--   This packages provide output routines to write the results of
--   monodromy homotopies, a numerical irreducible decomposition.

  procedure Write_Factor
              ( file : in file_type;
                eqs : in Standard_Complex_Poly_Systems.Poly_Sys;
                pts : in Standard_Complex_Solutions.Solution_List;
                fac : in Standard_Natural_Vectors.Link_to_Vector );
  procedure Write_Factor
              ( file : in file_type;
                eqs : in Standard_Complex_Laur_Systems.Laur_Sys;
                pts : in Standard_Complex_Solutions.Solution_List;
                fac : in Standard_Natural_Vectors.Link_to_Vector );
  procedure Write_Factor
              ( file : in file_type;
                eqs : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                pts : in DoblDobl_Complex_Solutions.Solution_List;
                fac : in Standard_Natural_Vectors.Link_to_Vector );
  procedure Write_Factor
              ( file : in file_type;
                eqs : in DoblDobl_Complex_Laur_Systems.Laur_Sys;
                pts : in DoblDobl_Complex_Solutions.Solution_List;
                fac : in Standard_Natural_Vectors.Link_to_Vector );
  procedure Write_Factor
              ( file : in file_type;
                eqs : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                pts : in QuadDobl_Complex_Solutions.Solution_List;
                fac : in Standard_Natural_Vectors.Link_to_Vector );
  procedure Write_Factor
              ( file : in file_type;
                eqs : in QuadDobl_Complex_Laur_Systems.Laur_Sys;
                pts : in QuadDobl_Complex_Solutions.Solution_List;
                fac : in Standard_Natural_Vectors.Link_to_Vector );

  -- DESCRIPTION :
  --   Writes the witness set representation of the factor fac
  --   to file, or to standard_output.

  -- ON ENTRY :
  --   file     must be either standard_output or a file opened for output;
  --   eqs      the equations of the witness set;
  --   pts      generic points in the witness set;
  --   fac      fac(k) is the k-th generic point in pts which belongs
  --            to the factor represented by the witness set.

  procedure Write_Factors
              ( file : in file_type;
                eqs : in Standard_Complex_Poly_Systems.Poly_Sys;
                pts : in Standard_Complex_Solutions.Solution_List;
                fac : in Standard_Natural_VecVecs.Link_to_VecVec );
  procedure Write_Factors
              ( file : in file_type;
                eqs : in Standard_Complex_Laur_Systems.Laur_Sys;
                pts : in Standard_Complex_Solutions.Solution_List;
                fac : in Standard_Natural_VecVecs.Link_to_VecVec );
  procedure Write_Factors
              ( file : in file_type;
                eqs : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                pts : in DoblDobl_Complex_Solutions.Solution_List;
                fac : in Standard_Natural_VecVecs.Link_to_VecVec );
  procedure Write_Factors
              ( file : in file_type;
                eqs : in DoblDobl_Complex_Laur_Systems.Laur_Sys;
                pts : in DoblDobl_Complex_Solutions.Solution_List;
                fac : in Standard_Natural_VecVecs.Link_to_VecVec );
  procedure Write_Factors
              ( file : in file_type;
                eqs : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                pts : in QuadDobl_Complex_Solutions.Solution_List;
                fac : in Standard_Natural_VecVecs.Link_to_VecVec );
  procedure Write_Factors
              ( file : in file_type;
                eqs : in QuadDobl_Complex_Laur_Systems.Laur_Sys;
                pts : in QuadDobl_Complex_Solutions.Solution_List;
                fac : in Standard_Natural_VecVecs.Link_to_VecVec );

  -- DESCRIPTION :
  --   Writes the witness set representation of the factors in fac
  --   to file, or to standard_output.

  -- ON ENTRY :
  --   file     must be either standard_output or a file opened for output;
  --   eqs      the equations of the witness set;
  --   pts      generic points in the witness set;
  --   fac      irreducible factors of the witness set.

  procedure Write_Components
              ( file : in file_type;
                eqs : in Standard_Complex_Poly_Systems.Array_of_Poly_Sys;
                pts : in Standard_Complex_Solutions.Array_of_Solution_Lists );
  procedure Write_Components
              ( file : in file_type;
                eqs : in Standard_Complex_Laur_Systems.Array_of_Laur_Sys;
                pts : in Standard_Complex_Solutions.Array_of_Solution_Lists );
  procedure Write_Components
              ( file : in file_type;
                eqs : in DoblDobl_Complex_Poly_Systems.Array_of_Poly_Sys;
                pts : in DoblDobl_Complex_Solutions.Array_of_Solution_Lists );
  procedure Write_Components
              ( file : in file_type;
                eqs : in DoblDobl_Complex_Laur_Systems.Array_of_Laur_Sys;
                pts : in DoblDobl_Complex_Solutions.Array_of_Solution_Lists );
  procedure Write_Components
              ( file : in file_type;
                eqs : in QuadDobl_Complex_Poly_Systems.Array_of_Poly_Sys;
                pts : in QuadDobl_Complex_Solutions.Array_of_Solution_Lists );
  procedure Write_Components
              ( file : in file_type;
                eqs : in QuadDobl_Complex_Laur_Systems.Array_of_Laur_Sys;
                pts : in QuadDobl_Complex_Solutions.Array_of_Solution_Lists );

  -- DESCRIPTION :
  --   Writes the pure dimensional components to file,
  --   or to standard_output.

  -- ON ENTRY :
  --   file     must be either standard_output or a file opened for output;
  --   eqs      the equations of the witness sets;
  --   pts      generic points in the witness sets.

  procedure Write_Decomposition
              ( file : in file_type;
                eqs : in Standard_Complex_Poly_Systems.Array_of_Poly_Sys;
                pts : in Standard_Complex_Solutions.Array_of_Solution_Lists;
                fac : in Standard_Natural_VecVecs.Array_of_VecVecs );
  procedure Write_Decomposition
              ( file : in file_type;
                eqs : in Standard_Complex_Laur_Systems.Array_of_Laur_Sys;
                pts : in Standard_Complex_Solutions.Array_of_Solution_Lists;
                fac : in Standard_Natural_VecVecs.Array_of_VecVecs );
  procedure Write_Decomposition
              ( file : in file_type;
                eqs : in DoblDobl_Complex_Poly_Systems.Array_of_Poly_Sys;
                pts : in DoblDobl_Complex_Solutions.Array_of_Solution_Lists;
                fac : in Standard_Natural_VecVecs.Array_of_VecVecs );
  procedure Write_Decomposition
              ( file : in file_type;
                eqs : in DoblDobl_Complex_Laur_Systems.Array_of_Laur_Sys;
                pts : in DoblDobl_Complex_Solutions.Array_of_Solution_Lists;
                fac : in Standard_Natural_VecVecs.Array_of_VecVecs );
  procedure Write_Decomposition
              ( file : in file_type;
                eqs : in QuadDobl_Complex_Poly_Systems.Array_of_Poly_Sys;
                pts : in QuadDobl_Complex_Solutions.Array_of_Solution_Lists;
                fac : in Standard_Natural_VecVecs.Array_of_VecVecs );
  procedure Write_Decomposition
              ( file : in file_type;
                eqs : in QuadDobl_Complex_Laur_Systems.Array_of_Laur_Sys;
                pts : in QuadDobl_Complex_Solutions.Array_of_Solution_Lists;
                fac : in Standard_Natural_VecVecs.Array_of_VecVecs );

  -- DESCRIPTION :
  --   Writes the numerical irreducible decomposition to file,
  --   or to standard_output.

  -- ON ENTRY :
  --   file     must be either standard_output or a file opened for output;
  --   eqs      the equations of the witness sets;
  --   pts      generic points in the witness sets;
  --   fac      numerical irreducible decomposition.

end Monodromy_Homotopies_io;
