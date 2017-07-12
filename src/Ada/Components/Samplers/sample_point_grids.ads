with text_io;                           use text_io;
with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Generic_Lists;
with Sample_Point_Lists;                use Sample_Point_Lists;

package Sample_Point_Grids is

-- DESCRIPTION :
--   This package provides an implementation for grids of sample points.
--   A grid contains for every component a list of sample points and
--   is represented as a list of sample point lists.

  package Lists_of_Standard_Sample_Lists is
    new Generic_Lists(Standard_Sample_List);
  package Lists_of_Multprec_Sample_Lists is
    new Generic_Lists(Multprec_Sample_List);

-- DATA TYPES :

  type Standard_Sample_Grid is new Lists_of_Standard_Sample_Lists.List;
  type Multprec_Sample_Grid is new Lists_of_Multprec_Sample_Lists.List;

-- THE STATE IS POLYNOMIAL OR LAURENT :

  procedure Set_Polynomial_Type ( laurent : in boolean );

  -- DESCRIPTION :
  --   If laurent, then the witness set is assumed to be defined
  --   by a Laurent polynomial system, otherwise, the witness set
  --   is defined by an ordinary polynomial system.
  --   This state determines the type of Sampling_Machine in the
  --   packages Sample_Points and Sample_Point_Lists.

-- CREATORS AS TYPE CONVERTORS :

  function Create ( grid : Standard_Sample_Grid )
                  return Array_of_Standard_Sample_Lists;
  function Create ( grid : Multprec_Sample_Grid )
                  return Array_of_Multprec_Sample_Lists;

  procedure Create ( samples : in Array_of_Standard_Sample_Lists;
                     grid,grid_last : in out Standard_Sample_Grid );
  procedure Create ( samples : in Array_of_Multprec_Sample_Lists;
                     grid,grid_last : in out Multprec_Sample_Grid );

  -- DESCRIPTION :
  --   With these creators, grids can have dual representation:
  --   as lists during buildup, and as arrays when dimensions are fixed.

-- SAMPLERS :

  procedure Sample ( samples : in Standard_Sample_List; nb : in natural32;
                     grid,grid_last : in out Standard_Sample_Grid );
  procedure Sample ( file : in file_type; full_output : in boolean;
                     samples : in Standard_Sample_List; nb : in natural32;
                     grid,grid_last : in out Standard_Sample_Grid );
  procedure Sample ( samples : in Standard_Sample_List; nb : in natural32;
                     grid,grid_last : in out Multprec_Sample_Grid );
  procedure Sample ( file : in file_type; full_output : in boolean;
                     samples : in Standard_Sample_List; nb : in natural32;
                     grid,grid_last : in out Multprec_Sample_Grid );

  -- DESCRIPTION :
  --   Creates a grid of sample points, generating from every point in
  --   the list nb new sample points.  If a file is given as argument,
  --   then diagnostics will be written on that file, completely if
  --   full_output is true, and only in summary if full_output is false.

-- DESTRUCTORS :

  procedure Shallow_Clear ( g : in out Standard_Sample_Grid );
  procedure Shallow_Clear ( g : in out Multprec_Sample_Grid );

  -- DESCRIPTION :
  --   Only destroys the encapsulating list structure, but not the
  --   lists with the sample points.

  procedure Deep_Clear ( g : in out Standard_Sample_Grid );
  procedure Deep_Clear ( g : in out Multprec_Sample_Grid );

  -- DESCRIPTION :
  --   Deallocation of all memory resources occupied by the grid.

end Sample_Point_Grids;
