with text_io;                            use text_io;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Lists_of_Integer_Vectors;           use Lists_of_Integer_Vectors;
with Standard_Complex_Poly_Systems;
with Multprec_Complex_Poly_Systems;
with Standard_Complex_Solutions;
with Irreducible_Component_Lists;        use Irreducible_Component_Lists;

package Irreducible_Decompositions is

-- DESCRIPTION :
--   This package offers a data abstraction to represent irreducible
--   decompositions of solution sets of polynomial systems.

  type Standard_Irreducible_Decomposition is private;
  type Multprec_Irreducible_Decomposition is private;

-- CREATORS :

  function Create ( k : integer32 ) return Standard_Irreducible_Decomposition;

  -- DESCRIPTION :
  --   Returns a structure ready to hold a sequence of embedded systems,
  --   of range 0..k where k is the top dimension.

  function Create ( p : Standard_Complex_Poly_Systems.Array_of_Poly_Sys )
                  return Standard_Irreducible_Decomposition;
  function Create ( p : Standard_Complex_Poly_Systems.Array_of_Poly_Sys )
                  return Multprec_Irreducible_Decomposition;

  -- DESCRIPTION :
  --   Returns a decomposition initialized with a sequence of embedded
  --   systems, of range 0..k, where k is the top dimension.

  function Create ( dc : Standard_Irreducible_Decomposition )
                  return Multprec_Irreducible_Decomposition;

  -- DESCRIPTION :
  --   Transfer the embedding and generic points over to the multi-precision
  --   version of the irreducible_decomposition.  There is sharing.

  procedure Add_Original
               ( dc : in out Multprec_Irreducible_Decomposition;
                 p : in Multprec_Complex_Poly_Systems.Link_to_Poly_Sys );

  -- DESCRIPTION :
  --   Adds the original polynomial system with multi-precision coefficients
  --   to the internal data structure.

  procedure Add_Embedding
               ( dc : in out Standard_Irreducible_Decomposition;
                 i : in integer32;
                 ep : in Standard_Complex_Poly_Systems.Link_to_Poly_Sys );

  -- DESCRIPTION :
  --   The decomposition is updated with an embedded system at level i,
  --   which means that i slices and i slack variables have been added
  --   to the original square polynomial system.
  --   Only the pointer is copied, so be aware of sharing.

  -- REQUIRED : dc has been created and i <= Top_Dimension(dc).

  procedure Add_Generic_Points
               ( dc : in out Standard_Irreducible_Decomposition;
                 i : in integer32;
                 gp : in Standard_Complex_Solutions.Solution_List );
  procedure Add_Generic_Points
               ( dc : in out Multprec_Irreducible_Decomposition;
                 i : in integer32;
                 gp : in Standard_Complex_Solutions.Solution_List );

  -- DESCRIPTION :
  --   The decomposition is updated with a list of generic points on
  --   the i-dimensional solution set.
  --   Also this procedure does not make a deep copy of the solution list,
  --   so be aware of sharing.

  -- REQUIRED : dc has been created and i <= Top_Dimension(dc).

  procedure Breakup ( file : in file_type; full_output : in boolean;
                      dc : in out Standard_Irreducible_Decomposition;
                      i : in integer32; method : in natural32;
                      stoptol,membtol : in double_float;
                      fp,fp_last : in out List );

  procedure Breakup ( file : in file_type; full_output : in boolean;
                      dc : in out Multprec_Irreducible_Decomposition;
                      i : in integer32; method,size : in natural32;
                      stoptol,membtol : in double_float;
                      fp,fp_last : in out List );

  -- DESCRIPTION :
  --   Breaks up the i-th dimensional solution set of dc.

  -- ON ENTRY :
  --   file      to write intermedate results and diagnostics to;
  --   full_output indicates whether all intermediate diagnositics of
  --             the sampling or whether only a summary will be written;
  --   dc        contains embeddings and solution sets;
  --   i         dimension of solution set of dc to be broken up;
  --   method    interpolating method used in breaking up
  --              = 0 : massive interpolate, comparison only at end,
  --              = 1 : incremental interpolate with linear projections,
  --              = 2 : interpolate with exploitation of span of component,
  --              = 3 : interpolate with central projections;
  --   size      size of multi-precision numbers;
  --   stoptol   tolerance to decide to stop interpolating;
  --   membtol   tolerance to determine membership.

  -- ON RETURN :
  --   dc        contains interpolation filters;
  --   fp        each component contains vector with its degree at entry i;
  --   fp_last   pointer to last element of the list fp.

  procedure Monodromy_Breakup
               ( file : in file_type;
                 dc : in out Standard_Irreducible_Decomposition;
                 i : in integer32; threshold : in natural32;
                 tol : in double_float;
                 fp,fp_last : in out List );

  -- DESCRIPTION :
  --   Applies the actions of the monodromy group to breakup the
  --   i-dimensional solution set of dc into irreducible components.

  -- ON ENTRY :
  --   file      for intermediate output and diagnostics;
  --   dc        contains embeddings and solution sets;
  --   i         dimension of current solution set to be decomposed;
  --   threshold is stabilizing threshold for the monodromy actions breakup;
  --   tol       tolerance to decide for equality of solution vectors.

  -- ON RETURN :
  --   dc        contains the grouping of points according to components;
  --   fp        each component contains vector with its degree at entry i;
  --   fp_last   pointer to last element of the list fp.

  procedure Filter ( file : in file_type;
                     dc : in out Standard_Irreducible_Decomposition;
                     i : in integer32; membtol : in double_float;
                     fp,fp_last : in out List; junkcnt : out natural32 );

  procedure Filter ( file : in file_type;
                     dc : in out Multprec_Irreducible_Decomposition;
                     i : in integer32; membtol : in double_float;
                     fp,fp_last : in out List; junkcnt : out natural32 );

  -- DESCRIPTION :
  --   Filters out the junk from the points at the i dimensional component.
  --   If a point belongs to a component of dimension higher than i,
  --   then Filter classifies the point as junk and removes it from dc.
  --   The list fp counts the junk points classified on each component.
  --   The output parameter junkcnt totals the number of junk points.

  procedure Homotopy_Filter
               ( file : in file_type;
                 dc : in out Standard_Irreducible_Decomposition;
                 i : in integer32; membtol : in double_float;
                 fp,fp_last : in out List; junkcnt : out natural32 );

  -- DESCRIPTION :
  --   Uses the homotopy membership test to perform the filtering.

-- SELECTORS :

  function Top_Dimension ( dc : Standard_Irreducible_Decomposition ) 
                         return integer32;
  function Top_Dimension ( dc : Multprec_Irreducible_Decomposition ) 
                         return integer32;

  -- DESCRIPTION :
  --   Returns the top dimension of the decomposition.

  function Embedding ( dc : Standard_Irreducible_Decomposition;
                       i : integer32 )
                     return Standard_Complex_Poly_Systems.Link_to_Poly_Sys;

  -- DESCRIPTION :
  --   Returns the i-th embedded system in the decomposition.

  function Original ( dc : Multprec_Irreducible_Decomposition )
                    return Multprec_Complex_Poly_Systems.Link_to_Poly_Sys;

  -- DESCRIPTION :
  --   Returns the original polynomial system with multi-precision
  --   coefficients.

  function Generic_Points ( dc : Standard_Irreducible_Decomposition;
                            i : integer32 )
                          return Standard_Complex_Solutions.Solution_List;
  function Generic_Points ( dc : Multprec_Irreducible_Decomposition;
                            i : integer32 )
                          return Standard_Complex_Solutions.Solution_List;

  -- DESCRIPTION :
  --   Returns the generic points on the i-dimensional components.

  function Components ( dc : Standard_Irreducible_Decomposition;
                        i : integer32 )
                      return Standard_Irreducible_Component_List;
  function Components ( dc : Multprec_Irreducible_Decomposition;
                        i : integer32 )
                      return Multprec_Irreducible_Component_List;

  -- DESCRIPTION :
  --   Returns the component list of the decomposition at level i.

-- DESTRUCTORS :

  procedure Clear ( dc : in out Standard_Irreducible_Decomposition );
  procedure Clear ( dc : in out Multprec_Irreducible_Decomposition );

  -- DESCRIPTION :
  --   Destroys the allocated memory.

private

  type Standard_Irreducible_Decomposition_Rep;
  type Standard_Irreducible_Decomposition is
    access Standard_Irreducible_Decomposition_Rep;

  type Multprec_Irreducible_Decomposition_Rep;
  type Multprec_Irreducible_Decomposition is
    access Multprec_Irreducible_Decomposition_Rep;

end Irreducible_Decompositions;
