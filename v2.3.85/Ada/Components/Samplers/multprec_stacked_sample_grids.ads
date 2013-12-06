with text_io;                            use text_io;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Multprec_Floating_Numbers;          use Multprec_Floating_Numbers;
with Multprec_Complex_Vectors;           use Multprec_Complex_Vectors;
with Multprec_Complex_VecVecs;           use Multprec_Complex_VecVecs;
with Sample_Points;                      use Sample_Points;
with Sample_Point_Lists;                 use Sample_Point_Lists;

package Multprec_Stacked_Sample_Grids is

-- DESCRIPTION :
--   A stacked sample grid has samples on parallel slices to apply
--   divided differences to interpolate a multivariate polynomial.
--   The grids managed by this package contain samples in arbitrary
--   precision floating-point numbers.

-- DATA STRUCTURES :

  type Link_to_Array_of_Multprec_Sample_Lists is
    access Array_of_Multprec_Sample_Lists;

  type Stacked_Sample_Grid;
  type Link_to_Stacked_Sample_Grid is access Stacked_Sample_Grid;
  type Array_of_Stacked_Sample_Grids is
    array ( integer32 range <> ) of Link_to_Stacked_Sample_Grid;

  type Stacked_Sample_Grid ( k,d : integer32 ) is record
    n : natural32;       -- n = #variables; k = #slices; d = degree 
    hyp : VecVec(1..k);  -- slices of the original sample list
    pts : Vector(0..d);  -- current interpolation points
    case k is
      when 1 => g : Link_to_Array_of_Multprec_Sample_Lists;
      when others => a : Array_of_Stacked_Sample_Grids(0..d);
                     spt : Multprec_Sample;
    end case;
  end record;

  -- CONVENTIONS :
  --   for k = 1 : g points to a rectangular sample grid;
  --   for k > 1 : a(i) contains samples for degree d of dimension k-1,
  --               spt contains sample for last point in grid.
  -- NOTICE :
  --   For the divided differences in the bootstrapping Newton,
  --   the range of "a" that is used goes from 1 to d.  With a
  --   little trick we save d samples at every level using "spt".

-- CREATORS :

  function Create ( file : file_type; sps : Standard_Sample_List;
                    size : natural32 ) return Stacked_Sample_Grid;

  -- DESCRIPTION :
  --   Returns a grid with samples to interpolate a polynomial of degree d
  --   in n variables to represent a component with generic points in sps.
  --   The parameter size controls the size of the multi-precision numbers.

  -- NOTE:
  --   d = Length_Of(sps), n = Number_of_Variables(Head_Of(sps)),
  --   and dimension = Number_of_Slices(Head_Of(sps)).

  function Create_Full ( file : file_type; sps : Standard_Sample_List;
                         size : natural32 ) return Stacked_Sample_Grid;

  -- DESCRIPTION :
  --   This procedure create the full grid, as required by the basic
  --   version of the trace form of the interpolating polynomial.
  --   Instead of creating the spt sample, it fills up grid.a(0).
  --   The same note as above on d, n, and dimension apply here.

-- DIAGNOSTICS :

  procedure Write_Errors 
              ( file : in file_type; grid : in Stacked_Sample_Grid );

  -- DESCRIPTION :
  --   Writes the "err" value for every sample in the grid on file.

  function Maximal_Error ( grid : Stacked_Sample_Grid ) return Floating_Number;

  -- DESCRIPTION :
  --   Returns the maximal "err" value for every sample in the grid.

  function Minimal_Distance ( grid : Stacked_Sample_Grid )
                            return Floating_Number;

  -- DESCRIPTION :
  --   Returns the minimal distance between the samples in the lists
  --   of the grid.

  procedure Write_Grid_Values
              ( file : in file_type; grid : in Stacked_Sample_Grid );

  -- DESCRIPTION :
  --   Writes the values of the samples in the grid at the directions
  --   perpendicular to the slicing hyperplanes.

  procedure Write_Full_Grid_Values
              ( file : in file_type; grid : in Stacked_Sample_Grid );

  -- DESCRIPTION :
  --   Writes the values of the samples in the grid at the directions
  --   perpendicular to the slicing hyperplanes.  This routine should
  --   be applied to the output of Create_Full, since the order of
  --   generating the slices is different from the order in Create.

-- SELECTOR :

  function Grid_Size ( n,d : natural32 ) return natural32;

  -- DESCRIPTION :
  --   Returns the number of samples in the grid to interpolate a
  --   polynomial in n variables of degree d.

  function Full_Grid_Size ( n,d : natural32 ) return natural32;

  -- DESCRIPTION :
  --   Returns d*(d+1)^(n-1) which is the size of the full grid, required
  --   by the basic version of the trace form of the interpolator.

-- DESTRUCTOR :

  procedure Clear ( L : in out Link_to_Array_of_Multprec_Sample_Lists );
  procedure Clear ( s : in out Link_to_Stacked_Sample_Grid );
  procedure Clear ( s : in out Array_of_Stacked_Sample_Grids );
  procedure Clear ( s : in out Stacked_Sample_Grid );

  -- DESCRIPTION :
  --   Deallocation of memory with "Deep_Clear".

end Multprec_Stacked_Sample_Grids;
