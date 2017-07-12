with text_io;                           use text_io;
with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Generic_Lists;
with DoblDobl_Sample_Lists;             use DoblDobl_Sample_Lists;

package DoblDobl_Sample_Grids is

-- DESCRIPTION :
--   This package provides an implementation for grids of sample points,
--   computed with double double precision.
--   A grid contains for every component a list of sample points and
--   is represented as a list of sample point lists.

  package Lists_of_DoblDobl_Sample_Lists is
    new Generic_Lists(DoblDobl_Sample_List);

-- DATA TYPES :

  type DoblDobl_Sample_Grid is new Lists_of_DoblDobl_Sample_Lists.List;

-- THE STATE IS POLYNOMIAL OR LAURENT :

  procedure Set_Polynomial_Type ( laurent : in boolean );

  -- DESCRIPTION :
  --   If laurent, then the witness set is assumed to be defined
  --   by a Laurent polynomial system, otherwise, the witness set
  --   is defined by an ordinary polynomial system.
  --   This state determines the type of Sampling_Machine in the
  --   packages Sample_Points and Sample_Point_Lists.

-- CREATORS AS TYPE CONVERTORS :

  function Create ( grid : DoblDobl_Sample_Grid )
                  return Array_of_DoblDobl_Sample_Lists;

  procedure Create ( samples : in Array_of_DoblDobl_Sample_Lists;
                     grid,grid_last : in out DoblDobl_Sample_Grid );

  -- DESCRIPTION :
  --   With these creators, grids can have dual representation:
  --   as lists during buildup, and as arrays when dimensions are fixed.

-- SAMPLERS :

  procedure Sample ( samples : in DoblDobl_Sample_List; nb : in natural32;
                     grid,grid_last : in out DoblDobl_Sample_Grid );
  procedure Sample ( file : in file_type; full_output : in boolean;
                     samples : in DoblDobl_Sample_List; nb : in natural32;
                     grid,grid_last : in out DoblDobl_Sample_Grid );

  -- DESCRIPTION :
  --   Creates a grid of sample points, generating from every point in
  --   the list nb new sample points.  If a file is given as argument,
  --   then diagnostics will be written on that file, completely if
  --   full_output is true, and only in summary if full_output is false.

-- DESTRUCTORS :

  procedure Shallow_Clear ( g : in out DoblDobl_Sample_Grid );

  -- DESCRIPTION :
  --   Only destroys the encapsulating list structure, but not the
  --   lists with the sample points.

  procedure Deep_Clear ( g : in out DoblDobl_Sample_Grid );

  -- DESCRIPTION :
  --   Deallocation of all memory resources occupied by the grid.

end DoblDobl_Sample_Grids;
