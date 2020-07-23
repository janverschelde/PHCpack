with text_io;                            use text_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Complex_Poly_Systems;
with Standard_Complex_Laur_Systems;
with DoblDobl_Complex_Poly_Systems;
with DoblDobl_Complex_Laur_Systems;
with QuadDobl_Complex_Poly_Systems;
with QuadDobl_Complex_Laur_Systems;
with Standard_Complex_Solutions;
with DoblDobl_Complex_Solutions;
with QuadDobl_Complex_Solutions;

package Black_Box_Simplex_Solvers is

-- DESCRIPTION :
--   This package provides a blackbox interface to the simplex system
--   solvers in PHCpack, to be called by phc -b.
 
  procedure Black_Box_Simplex_Solver
              ( p : in Standard_Complex_Poly_Systems.Poly_Sys;
                sols : out Standard_Complex_Solutions.Solution_List;
                fail : out boolean; verbose : integer32 := 0 );
  procedure Black_Box_Simplex_Solver
              ( p : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                sols : out DoblDobl_Complex_Solutions.Solution_List;
                fail : out boolean; verbose : integer32 := 0 );
  procedure Black_Box_Simplex_Solver
              ( p : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                sols : out QuadDobl_Complex_Solutions.Solution_List;
                fail : out boolean; verbose : integer32 := 0 );

  procedure Black_Box_Simplex_Solver
              ( p : in Standard_Complex_Laur_Systems.Laur_Sys;
                sols : out Standard_Complex_Solutions.Solution_List;
                fail : out boolean; verbose : integer32 := 0 );
  procedure Black_Box_Simplex_Solver
              ( p : in DoblDobl_Complex_Laur_Systems.Laur_Sys;
                sols : out DoblDobl_Complex_Solutions.Solution_List;
                fail : out boolean; verbose : integer32 := 0 );
  procedure Black_Box_Simplex_Solver
              ( p : in QuadDobl_Complex_Laur_Systems.Laur_Sys;
                sols : out QuadDobl_Complex_Solutions.Solution_List;
                fail : out boolean; verbose : integer32 := 0 );

  procedure Black_Box_Simplex_Solver
              ( file : in file_type;
                p : in Standard_Complex_Poly_Systems.Poly_Sys;
                sols : out Standard_Complex_Solutions.Solution_List;
                fail : out boolean; verbose : integer32 := 0 );
  procedure Black_Box_Simplex_Solver
              ( file : in file_type;
                p : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                sols : out DoblDobl_Complex_Solutions.Solution_List;
                fail : out boolean; verbose : integer32 := 0 );
  procedure Black_Box_Simplex_Solver
              ( file : in file_type;
                p : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                sols : out QuadDobl_Complex_Solutions.Solution_List;
                fail : out boolean; verbose : integer32 := 0 );

  procedure Black_Box_Simplex_Solver
              ( file : in file_type;
                p : in Standard_Complex_Laur_Systems.Laur_Sys;
                sols : out Standard_Complex_Solutions.Solution_List;
                fail : out boolean; verbose : integer32 := 0 );
  procedure Black_Box_Simplex_Solver
              ( file : in file_type;
                p : in DoblDobl_Complex_Laur_Systems.Laur_Sys;
                sols : out DoblDobl_Complex_Solutions.Solution_List;
                fail : out boolean; verbose : integer32 := 0 );
  procedure Black_Box_Simplex_Solver
              ( file : in file_type;
                p : in QuadDobl_Complex_Laur_Systems.Laur_Sys;
                sols : out QuadDobl_Complex_Solutions.Solution_List;
                fail : out boolean; verbose : integer32 := 0 );

  -- DESCRIPTION :
  --   Returns the solutions of p if p is simplex,
  --   otherwise fail is true on return.

  -- ON ENTRY :
  --   file     for results of the reporting root refiner;
  --   p        a (Laurent) polynomial system, presumed to be a simplex.

  -- ON RETURN :
  --   sols     solutions of p if not fail;
  --   fail     true if p is not a simplex system, otherwise false.

end Black_Box_Simplex_Solvers;
