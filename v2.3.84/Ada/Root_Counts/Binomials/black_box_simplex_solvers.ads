with text_io;                            use text_io;
with Standard_Complex_Poly_Systems;      use Standard_Complex_Poly_Systems;
with Standard_Complex_Laur_Systems;      use Standard_Complex_Laur_Systems;
with Standard_Complex_Solutions;         use Standard_Complex_Solutions;

package Black_Box_Simplex_Solvers is

-- DESCRIPTION :
--   This package provides a blackbox interface to the simplex system
--   solvers in PHCpack, to be called by phc -b.
 
  procedure Black_Box_Simplex_Solver
              ( p : in Poly_Sys; sols : out Solution_List;
                fail : out boolean );
  procedure Black_Box_Simplex_Solver
              ( p : in Laur_Sys; sols : out Solution_List;
                fail : out boolean );
  procedure Black_Box_Simplex_Solver
              ( file : in file_type;
                p : in Poly_Sys; sols : out Solution_List;
                fail : out boolean );
  procedure Black_Box_Simplex_Solver
              ( file : in file_type;
                p : in Laur_Sys; sols : out Solution_List;
                fail : out boolean );

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
