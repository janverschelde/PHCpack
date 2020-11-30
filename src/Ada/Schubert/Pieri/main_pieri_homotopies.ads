with text_io;                            use text_io;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Natural_Matrices;
with Standard_Complex_VecMats;           use Standard_Complex_VecMats;
with Standard_Complex_Poly_Systems;      use Standard_Complex_Poly_Systems;
with Standard_Complex_Poly_Matrices;
with Standard_Complex_Solutions;         use Standard_Complex_Solutions;

package Main_Pieri_Homotopies is

-- DESCRIPTION :
--   The two main procedures setup Pieri homotopies 
--   to intersect m-planes in n-space. 

  function Create_Hypersurface_System
             ( n : natural32; locmap : Standard_Natural_Matrices.Matrix;
               xpm : Standard_Complex_Poly_Matrices.Matrix; planes : VecMat )
             return Poly_Sys;

  -- DESCRIPTION :
  --   Returns the polynomial system that collects the hypersurface
  --   intersection conditions for meeting the given m-planes.
  --   The system is localized according to the given localization map.

  function Solution_Planes ( locmap : Standard_Natural_Matrices.Matrix;
                             vm : VecMat ) return Solution_List;

  -- DESCRIPTION :
  --   Returns the representation of the vector of planes as a solution list.

  function Square ( n : natural32; p : Poly_Sys ) return Poly_Sys;

  -- DESCRIPTION :
  --   Returns a n-by-n system, by adding random linear combinations of 
  --   the polynomials p(i), i > n, to the first n polynomials.

  procedure Path_Following_with_Cheater
               ( file : in file_type; start,target : in Poly_Sys;
                 sols : in out Solution_List; report : in boolean );

  -- DESCRIPTION :
  --   Calls the standard continuation routines to solve a specific
  --   target system, starting at the solutions of the start system.
  --   This is the usual linear cheater between start and target system,
  --   although it will take nonsquare inputs, it should only be used
  --   for hypersurface intersection conditions.

  -- REQUIRED : not Is_Null(sols).

  procedure Main ( n,d : in natural32 );
  procedure Main ( file : in file_type; n,d : in natural32 );

  -- DESCRIPTION :
  --   Finds all d-planes that meet linear subspaces in n-space.

  -- ON ENTRY :
  --   file     output file (optional) must be opened for writing;
  --   n        dimension of the ambient space;
  --   d        dimension of the output planes.

end Main_Pieri_Homotopies;
