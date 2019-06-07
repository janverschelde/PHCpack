with text_io;                           use text_io;
with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Floating_Numbers;         use Standard_Floating_Numbers;
with Standard_Integer_Vectors;
with Standard_Floating_Vectors;
with Standard_Complex_Vectors;
with DoblDobl_Complex_Vectors;
with QuadDobl_Complex_Vectors;
with Standard_Complex_Laurentials;
with Standard_Complex_Laur_Systems;
with DoblDobl_Complex_Laurentials;
with DoblDobl_Complex_Laur_Systems;
with QuadDobl_Complex_Laurentials;
with QuadDobl_Complex_Laur_Systems;
with Lists_of_Floating_Vectors;
with Arrays_of_Floating_Vector_Lists;
with Floating_Mixed_Subdivisions;       use Floating_Mixed_Subdivisions;
with Standard_Complex_Solutions;
with DoblDobl_Complex_Solutions;
with QuadDobl_Complex_Solutions;

package Stable_Polyhedral_Continuation is

-- DESCRIPTION :
--   This package collects dedicated operations to deal with the extra stable
--   mixed cells arising in the computation of the stable mixed volume.
--   The main functionality of this package is to allow the computation of
--   solutions with zero components of random coefficient systems, by
--   polyhedral continuation for the nonzero components of lower dimensional
--   polynomial systems derived from substituting the zeroes.

-- PREPARATORY FUNCTIONS :
--   The Eliminate_Zeroes removes those components that match zero types,
--   whereas Substitute_Zeroes take into account that an entire monomial
--   vanishes as soon as there is one variable present of zero type.

  function Eliminate_Zeroes
             ( v : Standard_Floating_Vectors.Vector;
               z : Standard_Integer_Vectors.Vector; nbz : integer32 )
             return Standard_Floating_Vectors.Vector;

  -- DESCRIPTION :
  --   Removes components of v that match zero types.

  function Vanish_by_Zeroes
             ( x : Standard_Floating_Vectors.Vector;
               z : Standard_Integer_Vectors.Vector; nbz : integer32 )
             return boolean;

  -- DESCRIPTION :
  --   If nbz <= 0, then false is returned, otherwise for nbz > 0,
  --   true is returned if the point x vanishes along the zero type z.
  --   A point x vanishes if it contains a nonzero entry corresponding
  --   to a zero entry in z.

  function Substitute_Zeroes
             ( pts : Lists_of_Floating_Vectors.List; nbz : integer32;
               z : Standard_Integer_Vectors.Vector )
             return Lists_of_Floating_Vectors.List;

  -- DESCRIPTION :
  --   If nbz > 0, then the points in pts vanish or have those entries
  --   that correspond to the zero type in z removed.

  function Substitute_Zeroes
             ( pts : Arrays_of_Floating_Vector_Lists.Array_of_Lists;
               nbz : integer32; z : Standard_Integer_Vectors.Vector )
             return Arrays_of_Floating_Vector_Lists.Array_of_Lists;

  -- DESCRIPTION :
  --   Applies the Substitute_Zeroes to each list in pts.

  function Substitute_Zeroes
             ( mic : Mixed_Cell; nbz : integer32;
               z : Standard_Integer_Vectors.Vector ) return Mixed_Cell;

  -- DESCRIPTION :
  --   If nbz > 0, then the cell on return has in its inner normal
  --   the corresponding zero components removed.  Either the points
  --   that span the cell vanish or have also the corresponding zero
  --   components removed.

  function Filter ( pts : Arrays_of_Floating_Vector_Lists.Array_of_Lists )
                  return Arrays_of_Floating_Vector_Lists.Array_of_Lists;

  -- DESCRIPTION :
  --   Removes all lists from pts with one or no elements.

  function Filter ( pts : Arrays_of_Floating_Vector_Lists.Array_of_Lists;
                    mix : Standard_Integer_Vectors.Link_to_Vector )
                  return Standard_Integer_Vectors.Link_to_Vector;

  -- DESCRIPTION :
  --   Returns a new type of mixture, according to Filter(pts).

  function Filter ( mic : Mixed_Cell ) return Mixed_Cell;

  -- DESCRIPTION :
  --   Removes all lists from mic with one or no elements.

-- DEFINING POLYHEDRAL HOMOTOPIES :
--   A polyhedral homotopy is defined by a lifting and mixed cells.

  procedure Assign_Multiplicity
               ( s : in out Standard_Complex_Solutions.Solution_List;
                 vol : in natural32 );
  procedure Assign_Multiplicity
               ( s : in out DoblDobl_Complex_Solutions.Solution_List;
                 vol : in natural32 );
  procedure Assign_Multiplicity
               ( s : in out QuadDobl_Complex_Solutions.Solution_List;
                 vol : in natural32 );

  -- DESCRIPTION :
  --   Assigns the multiplicity of the solutions, using the volume of
  --   the corresponding mixed cell.  If vol = Length_Of(s), then all
  --   solutions have multiplicity one, otherwise the multiplicity of
  --   each solution equals vol/Length(s).

  procedure Silent_Polyhedral_Continuation
              ( q : in Standard_Complex_Laur_Systems.Laur_Sys;
                mix : in Standard_Integer_Vectors.Link_to_Vector;
                lif : in Arrays_of_Floating_Vector_Lists.Array_of_Lists;
                mic : in Mixed_Cell; nbz : in integer32; vol : in natural32;
                ztp : in Standard_Integer_Vectors.Vector;
                sols : out Standard_Complex_Solutions.Solution_List;
                verbose : in integer32 := 0 );
  procedure Silent_Polyhedral_Continuation
              ( q : in DoblDobl_Complex_Laur_Systems.Laur_Sys;
                mix : in Standard_Integer_Vectors.Link_to_Vector;
                lif : in Arrays_of_Floating_Vector_Lists.Array_of_Lists;
                mic : in Mixed_Cell; nbz : in integer32; vol : in natural32;
                ztp : in Standard_Integer_Vectors.Vector;
                sols : out DoblDobl_Complex_Solutions.Solution_List;
                verbose : in integer32 := 0 );
  procedure Silent_Polyhedral_Continuation
              ( q : in QuadDobl_Complex_Laur_Systems.Laur_Sys;
                mix : in Standard_Integer_Vectors.Link_to_Vector;
                lif : in Arrays_of_Floating_Vector_Lists.Array_of_Lists;
                mic : in Mixed_Cell; nbz : in integer32; vol : in natural32;
                ztp : in Standard_Integer_Vectors.Vector;
                sols : out QuadDobl_Complex_Solutions.Solution_List;
                verbose : in integer32 := 0 );
  procedure Reporting_Polyhedral_Continuation
              ( file : in file_type;
                q : in Standard_Complex_Laur_Systems.Laur_Sys;
                mix : in Standard_Integer_Vectors.Link_to_Vector;
                lif : in Arrays_of_Floating_Vector_Lists.Array_of_Lists;
                mic : in Mixed_Cell; nbz : in integer32; vol : in natural32;
                ztp : in Standard_Integer_Vectors.Vector;
                sols : out Standard_Complex_Solutions.Solution_List;
                verbose : in integer32 := 0 );
  procedure Reporting_Polyhedral_Continuation
              ( file : in file_type;
                q : in DoblDobl_Complex_Laur_Systems.Laur_Sys;
                mix : in Standard_Integer_Vectors.Link_to_Vector;
                lif : in Arrays_of_Floating_Vector_Lists.Array_of_Lists;
                mic : in Mixed_Cell; nbz : in integer32; vol : in natural32;
                ztp : in Standard_Integer_Vectors.Vector;
                sols : out DoblDobl_Complex_Solutions.Solution_List;
                verbose : in integer32 := 0 );
  procedure Reporting_Polyhedral_Continuation
              ( file : in file_type;
                q : in QuadDobl_Complex_Laur_Systems.Laur_Sys;
                mix : in Standard_Integer_Vectors.Link_to_Vector;
                lif : in Arrays_of_Floating_Vector_Lists.Array_of_Lists;
                mic : in Mixed_Cell; nbz : in integer32; vol : in natural32;
                ztp : in Standard_Integer_Vectors.Vector;
                sols : out QuadDobl_Complex_Solutions.Solution_List;
                verbose : in integer32 := 0 );

  -- DESCRIPTION :
  --   Performs polyhedral continuation to compute the solutions with
  --   zero components for one extra stable mixed cell.

  -- REQUIRED : nbz < q'last, i.e.: not all components are zero.

  -- ON ENTRY :
  --   file     for intermediate output and diagnostics;
  --   q        a random coefficient system;
  --   mix      type of mixture of the supports;
  --   lif      lifted supports;
  --   mic      one extra stable mixed cell;
  --   nbz      number of zeroes in the solutions;
  --   vol      volume of the cell;
  --   verbose  the verbose level.

  -- ON RETURN :
  --   sols     solutions of q with zero components, 
  --            corresponding to the mixed cell.

  procedure Silent_Polyhedral_Continuation
              ( q : in Standard_Complex_Laur_Systems.Laur_Sys;
                b : in double_float; 
                mix : in Standard_Integer_Vectors.Link_to_Vector;
                lif : in Arrays_of_Floating_Vector_Lists.Array_of_Lists;
                mic : in out Mixed_Cell; -- k : in integer32;
                sols : out Standard_Complex_Solutions.Solution_List;
                verbose : in integer32 := 0 );
  procedure Silent_Polyhedral_Continuation
              ( q : in DoblDobl_Complex_Laur_Systems.Laur_Sys;
                b : in double_float; 
                mix : in Standard_Integer_Vectors.Link_to_Vector;
                lif : in Arrays_of_Floating_Vector_Lists.Array_of_Lists;
                mic : in out Mixed_Cell; -- k : in integer32;
                sols : out DoblDobl_Complex_Solutions.Solution_List;
                verbose : in integer32 := 0 );
  procedure Silent_Polyhedral_Continuation
              ( q : in QuadDobl_Complex_Laur_Systems.Laur_Sys;
                b : in double_float; 
                mix : in Standard_Integer_Vectors.Link_to_Vector;
                lif : in Arrays_of_Floating_Vector_Lists.Array_of_Lists;
                mic : in out Mixed_Cell; -- k : in integer32;
                sols : out QuadDobl_Complex_Solutions.Solution_List;
                verbose : in integer32 := 0 );
  procedure Reporting_Polyhedral_Continuation
              ( file : in file_type;
                q : in Standard_Complex_Laur_Systems.Laur_Sys;
                b : in double_float; 
                mix : in Standard_Integer_Vectors.Link_to_Vector;
                lif : in Arrays_of_Floating_Vector_Lists.Array_of_Lists;
                mic : in out Mixed_Cell; k : in integer32;
                sols : in out Standard_Complex_Solutions.Solution_List;
                verbose : in integer32 := 0 );
  procedure Reporting_Polyhedral_Continuation
              ( file : in file_type;
                q : in DoblDobl_Complex_Laur_Systems.Laur_Sys;
                b : in double_float; 
                mix : in Standard_Integer_Vectors.Link_to_Vector;
                lif : in Arrays_of_Floating_Vector_Lists.Array_of_Lists;
                mic : in out Mixed_Cell; k : in integer32;
                sols : in out DoblDobl_Complex_Solutions.Solution_List;
                verbose : in integer32 := 0 );
  procedure Reporting_Polyhedral_Continuation
              ( file : in file_type;
                q : in QuadDobl_Complex_Laur_Systems.Laur_Sys;
                b : in double_float; 
                mix : in Standard_Integer_Vectors.Link_to_Vector;
                lif : in Arrays_of_Floating_Vector_Lists.Array_of_Lists;
                mic : in out Mixed_Cell; k : in integer32;
                sols : in out QuadDobl_Complex_Solutions.Solution_List;
                verbose : in integer32 := 0 );

  -- DESCRIPTION :
  --   Performs polyhedral continuation to compute the solutions with
  --   zero components for one extra stable mixed cell.

  -- ON ENTRY :
  --   file     for intermediate output and diagnostics;
  --   q        a random coefficient system;
  --   b        lifting bound used for the artificial origins;
  --   mix      type of mixture of the supports;
  --   lif      lifted supports;
  --   mic      one extra stable mixed cell;
  --   k        counter of cell;
  --   verbose  the verbose level.
 
  -- ON RETURN :
  --   mic      might be refined for mixed volume computation;
  --   sols     solutions corresponding to the mixed cell.

  procedure Silent_Polyhedral_Continuation
              ( q : in Standard_Complex_Laur_Systems.Laur_Sys;
                b : in double_float; 
                mix : in Standard_Integer_Vectors.Link_to_Vector;
                lif : in Arrays_of_Floating_Vector_Lists.Array_of_Lists;
                mcc : in Mixed_Subdivision;
                sols : in out Standard_Complex_Solutions.Solution_List;
                verbose : in integer32 := 0 );
  procedure Silent_Polyhedral_Continuation
              ( q : in DoblDobl_Complex_Laur_Systems.Laur_Sys;
                b : in double_float; 
                mix : in Standard_Integer_Vectors.Link_to_Vector;
                lif : in Arrays_of_Floating_Vector_Lists.Array_of_Lists;
                mcc : in Mixed_Subdivision;
                sols : in out DoblDobl_Complex_Solutions.Solution_List;
                verbose : in integer32 := 0 );
  procedure Silent_Polyhedral_Continuation
              ( q : in QuadDobl_Complex_Laur_Systems.Laur_Sys;
                b : in double_float; 
                mix : in Standard_Integer_Vectors.Link_to_Vector;
                lif : in Arrays_of_Floating_Vector_Lists.Array_of_Lists;
                mcc : in Mixed_Subdivision;
                sols : in out QuadDobl_Complex_Solutions.Solution_List;
                verbose : in integer32 := 0 );
  procedure Reporting_Polyhedral_Continuation
              ( file : in file_type; 
                q : in Standard_Complex_Laur_Systems.Laur_Sys;
                b : in double_float; 
                mix : in Standard_Integer_Vectors.Link_to_Vector;
                lif : in Arrays_of_Floating_Vector_Lists.Array_of_Lists;
                mcc : in Mixed_Subdivision;
                sols : in out Standard_Complex_Solutions.Solution_List;
                verbose : in integer32 := 0 );
  procedure Reporting_Polyhedral_Continuation
              ( file : in file_type; 
                q : in DoblDobl_Complex_Laur_Systems.Laur_Sys;
                b : in double_float; 
                mix : in Standard_Integer_Vectors.Link_to_Vector;
                lif : in Arrays_of_Floating_Vector_Lists.Array_of_Lists;
                mcc : in Mixed_Subdivision;
                sols : in out DoblDobl_Complex_Solutions.Solution_List;
                verbose : in integer32 := 0 );
  procedure Reporting_Polyhedral_Continuation
              ( file : in file_type; 
                q : in QuadDobl_Complex_Laur_Systems.Laur_Sys;
                b : in double_float; 
                mix : in Standard_Integer_Vectors.Link_to_Vector;
                lif : in Arrays_of_Floating_Vector_Lists.Array_of_Lists;
                mcc : in Mixed_Subdivision;
                sols : in out QuadDobl_Complex_Solutions.Solution_List;
                verbose : in integer32 := 0 );

  -- DESCRIPTION :
  --   Performs polyhedral continuation to compute the solutions with
  --   zero components for a list of extra stable mixed cells.

  -- ON ENTRY :
  --   file     for intermediate output and diagnostics;
  --   q        a random coefficient system;
  --   b        lifting bound used for the artificial origins;
  --   mix      type of mixture of the supports;
  --   lif      lifted supports;
  --   mcc      extra stable mixed cells;
  --   verbose  the verbose level.

  -- ON RETURN :
  --   sols     solutions with zero components of q.
 
end Stable_Polyhedral_Continuation;
