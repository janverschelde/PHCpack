with text_io;                            use text_io;
with Standard_Complex_Poly_Systems;
with Standard_Complex_Laur_Systems;
with Standard_Monomial_Maps;             use Standard_Monomial_Maps;

package Black_Box_Binomial_Solvers is

-- DESCRIPTION :
--   This package provides a blackbox interface to the binomial system
--   solvers in PHCpack, to be called by phc -b.
 
  procedure Black_Box_Binomial_Solver
              ( p : in Standard_Complex_Poly_Systems.Poly_Sys;
                pure : in boolean;
                sols : out Link_to_Array_of_Monomial_Map_Lists;
                fail : out boolean );
  procedure Black_Box_Binomial_Solver
              ( p : in Standard_Complex_Laur_Systems.Laur_Sys;
                pure : in boolean;
                sols : out Link_to_Array_of_Monomial_Map_Lists;
                fail : out boolean );
  procedure Black_Box_Binomial_Solver
              ( file : in file_type;
                p : in Standard_Complex_Poly_Systems.Poly_Sys;
                sols : out Link_to_Array_of_Monomial_Map_Lists;
                fail : out boolean );
  procedure Black_Box_Binomial_Solver
              ( file : in file_type;
                p : in Standard_Complex_Laur_Systems.Laur_Sys;
                sols : out Link_to_Array_of_Monomial_Map_Lists;
                fail : out boolean );

  -- DESCRIPTION :
  --   Returns the solutions of p if p is simplex,
  --   otherwise fail is true on return.

  -- ON ENTRY :
  --   file     for intermediate diagnostic results to monitor progress;
  --   p        a (Laurent) polynomial system, presumed to be binomial;
  --   pure     true if only the pure top dimensional sets are wanted;
  --            if omitted, then the user is prompted for this value.

  -- ON RETURN :
  --   sols     solution maps of p if not fail;
  --   fail     true if p is not a simplex system, otherwise false.

end Black_Box_Binomial_Solvers;
