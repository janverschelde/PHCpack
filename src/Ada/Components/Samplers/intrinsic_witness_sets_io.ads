with text_io;                            use text_io;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Complex_Matrices;
with DoblDobl_Complex_Matrices;
with QuadDobl_Complex_Matrices;
with Standard_Complex_VecMats;
with Standard_Complex_Poly_Systems;
with Standard_Complex_Laur_Systems;
with DoblDobl_Complex_Poly_Systems;
with DoblDobl_Complex_Laur_Systems;
with QuadDobl_Complex_Poly_Systems;
with QuadDobl_Complex_Laur_Systems;
with Standard_Complex_Solutions;
with DoblDobl_Complex_Solutions;
with QuadDobl_Complex_Solutions;
with Standard_Continuation_Data;
with DoblDobl_Continuation_Data;
with QuadDobl_Continuation_Data;

package Intrinsic_Witness_Sets_io is

-- DESCRIPTION :
--   This package provides procedures to write a witness set given
--   in intrinsic coordinates to file.  A witness stone is used as
--   stepping stone in the equation-by-equation solver.

-- INPUT/OUTPUT of WITNESS STONES :

  procedure Read_Witness_Stone
              ( d,k : out natural32;
                s : out Standard_Complex_Solutions.Solution_List;
                p : out Standard_Complex_Matrices.Link_to_Matrix );
  procedure Read_Witness_Stone
              ( f : out Standard_Complex_Poly_Systems.Link_to_Poly_Sys;
                d,k : out natural32;
                s : out Standard_Complex_Solutions.Solution_List;
                p : out Standard_Complex_Matrices.Link_to_Matrix );

  -- DESCRIPTION :
  --   Prompts the user for a file name, and then reads a witness set.

  -- ON RETURN :
  --   f        embedding of the defining polynomial system,
  --            if solving is done without straight-line program;
  --   d        dimension of the solution set;
  --   k        number of equations in the original polynomial system;
  --   s        intrinsic coordinates of witness points;
  --   p        orthonormal basis of the plane defining the witness set.
  
  procedure Write_Witness_Stone
              ( file : in file_type; filename : in string; nv,d : in natural32;
                s : in Standard_Complex_Solutions.Solution_List;
                p : in Standard_Complex_Matrices.Matrix );
  procedure Write_Witness_Stone
              ( file : in file_type; filename : in string; nv,d : in natural32;
                f : in Standard_Complex_Poly_Systems.Poly_Sys;
                s : in Standard_Complex_Solutions.Solution_List;
                p : in Standard_Complex_Matrices.Matrix );

  -- DESCRIPTION :
  --   Writes a witness stone, just like "Write_Witness_Set" above.
  --   A witness stone is an abbreviation for stepping stone of a
  --   witness set, obtained after solving some equations in a system.
  --   When f is omitted, the hyperplanes defined by p are written.

-- OUTPUT OF ONE WITNESS SET :

  procedure Write_Witness_Set
              ( file : in file_type; nv,d : in natural32;
                f : in Standard_Complex_Poly_Systems.Poly_Sys;
                s : in Standard_Complex_Solutions.Solution_List;
                p : in Standard_Complex_Matrices.Matrix );
  procedure Write_Witness_Set
              ( file : in file_type; nv,d : in natural32;
                f : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                s : in DoblDobl_Complex_Solutions.Solution_List;
                p : in DoblDobl_Complex_Matrices.Matrix );
  procedure Write_Witness_Set
              ( file : in file_type; nv,d : in natural32;
                f : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                s : in QuadDobl_Complex_Solutions.Solution_List;
                p : in QuadDobl_Complex_Matrices.Matrix );

  -- DESCRIPTION :
  --   Writes the witness set defined by f, s, and p to file,
  --   in standard double, double double, or quad double precision.

  procedure Write_Recentered_Witness_Set
              ( file : in file_type; nv,d : in natural32;
                f : in Standard_Complex_Poly_Systems.Poly_Sys;
                s : in Standard_Complex_Solutions.Solution_List;
                p : in Standard_Complex_Matrices.Matrix );
  procedure Write_Recentered_Witness_Set
              ( file : in file_type; nv,d : in natural32;
                f : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                s : in DoblDobl_Complex_Solutions.Solution_List;
                p : in DoblDobl_Complex_Matrices.Matrix );
  procedure Write_Recentered_Witness_Set
              ( file : in file_type; nv,d : in natural32;
                f : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                s : in QuadDobl_Complex_Solutions.Solution_List;
                p : in QuadDobl_Complex_Matrices.Matrix );

  -- DESCRIPTION :
  --   Writes the witness set defined by f, s, and p to file.
  --   The solutions in s are already in extrinsic coordinates,
  --   as computed in a recentered version of intrinsic tracking.

  procedure Write_Witness_Set_to_File
              ( filename : in string; n,k : in natural32;
                f : in Standard_Complex_Poly_Systems.Poly_Sys;
                p : in Standard_Complex_Matrices.Matrix;
                sols : in Standard_Continuation_Data.Solu_Info_Array );
  procedure Write_Witness_Set_to_File
              ( filename : in string; n,k : in natural32;
                f : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                p : in DoblDobl_Complex_Matrices.Matrix;
                sols : in DoblDobl_Continuation_Data.Solu_Info_Array );
  procedure Write_Witness_Set_to_File
              ( filename : in string; n,k : in natural32;
                f : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                p : in QuadDobl_Complex_Matrices.Matrix;
                sols : in QuadDobl_Continuation_Data.Solu_Info_Array );

  -- DESCRIPTION :
  --   Prompts the user for a file name and then writes the
  --   witness set to file, calling Write_Recentered_Witness_Set.

  procedure Write_Witness_Set
              ( file : in file_type; filename : in string;
                nv,d : in natural32;
                s : in Standard_Complex_Solutions.Solution_List;
                p : in Standard_Complex_Matrices.Matrix );
  procedure Write_Witness_Set
              ( file : in file_type; filename : in string;
                nv,d : in natural32;
                s : in DoblDobl_Complex_Solutions.Solution_List;
                p : in DoblDobl_Complex_Matrices.Matrix );
  procedure Write_Witness_Set
              ( file : in file_type; filename : in string;
                nv,d : in natural32;
                s : in QuadDobl_Complex_Solutions.Solution_List;
                p : in QuadDobl_Complex_Matrices.Matrix );

  -- DESCRIPTION :
  --   Writes the solutions in a witness set to file with name in filename
  --   concatenated with _wd, where d is the dimension of the set,
  --   in standard double, double double, or quad double precision.

  procedure Write_Witness_Set
              ( file : in file_type; filename : in string;
                nv,d : in natural32;
                f : in Standard_Complex_Poly_Systems.Poly_Sys;
                s : in Standard_Complex_Solutions.Solution_List;
                p : in Standard_Complex_Matrices.Matrix );
  procedure Write_Witness_Set
              ( file : in file_type; filename : in string;
                nv,d : in natural32;
                f : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                s : in DoblDobl_Complex_Solutions.Solution_List;
                p : in DoblDobl_Complex_Matrices.Matrix );
  procedure Write_Witness_Set
              ( file : in file_type; filename : in string;
                nv,d : in natural32;
                f : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                s : in QuadDobl_Complex_Solutions.Solution_List;
                p : in QuadDobl_Complex_Matrices.Matrix );

  -- DESCRIPTION :
  --   Writes a witness set to file with name given in filename
  --   concatenated with _wd, where d is the dimension of the set,
  --   in standard double, double double, or quad double precision,
  --   for a set defined by an ordinary polynomial system.

  procedure Write_Witness_Set
              ( file : in file_type; filename : in string;
                nv,d : in natural32;
                f : in Standard_Complex_Laur_Systems.Laur_Sys;
                s : in Standard_Complex_Solutions.Solution_List;
                p : in Standard_Complex_Matrices.Matrix );
  procedure Write_Witness_Set
              ( file : in file_type; filename : in string;
                nv,d : in natural32;
                f : in DoblDobl_Complex_Laur_Systems.Laur_Sys;
                s : in DoblDobl_Complex_Solutions.Solution_List;
                p : in DoblDobl_Complex_Matrices.Matrix );
  procedure Write_Witness_Set
              ( file : in file_type; filename : in string;
                nv,d : in natural32;
                f : in QuadDobl_Complex_Laur_Systems.Laur_Sys;
                s : in QuadDobl_Complex_Solutions.Solution_List;
                p : in QuadDobl_Complex_Matrices.Matrix );

  -- DESCRIPTION :
  --   Writes a witness set to file with name given in filename
  --   concatenated with _wd, where d is the dimension of the set,
  --   in standard double, double double, or quad double precision,
  --   for a set defined by a Laurent polynomial system.

-- OUTPUT OF AN ARRAY OF WITNESS SETS :

  procedure Write_Witness_Sets
              ( file : in file_type; filename : in string; nv : in natural32;
                witset : in Standard_Complex_Solutions.Array_of_Solution_Lists;
                planes : in Standard_Complex_VecMats.VecMat );
  procedure Write_Witness_Sets
              ( file : in file_type; filename : in string; nv : in natural32;
                f : in Standard_Complex_Poly_Systems.Poly_Sys;
                witset : in Standard_Complex_Solutions.Array_of_Solution_Lists;
                planes : in Standard_Complex_VecMats.VecMat );

  -- DESCRIPTION :
  --   For every nonzero witness set, a new file with extension _wd
  --   is created with on it the witness set of dimension d.

  -- ON ENTRY :
  --   file     is the main output file, already opened for output;
  --   filename is the name of the main output file;
  --   nv       number of variables;
  --   f        polynomial system defining the solution sets;
  --   witset   array of witness sets, in intrinsic format, 
  --            witset(i) contains witness points on component of
  --            dimension nv - i;
  --   planes   generators for planes defining witness sets,
  --            for every witset(i), there corresponds a planes(i).

end Intrinsic_Witness_Sets_io;
