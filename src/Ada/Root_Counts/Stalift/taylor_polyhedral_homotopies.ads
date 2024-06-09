with text_io;                            use text_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Integer_Vectors;
with Standard_Floating_VecVecs;
with Standard_Complex_VecVecs;
with Lists_of_Floating_Vectors;          use Lists_of_Floating_Vectors;
with Arrays_of_Floating_Vector_Lists;
with Floating_Mixed_Subdivisions;        use Floating_Mixed_Subdivisions;
with Double_Taylor_Homotopies;           use Double_Taylor_Homotopies;

package Taylor_Polyhedral_Homotopies is

-- DESCRIPTION :
--   A Taylor polyhedral homotopy applies a Taylor development to
--   the continuation parameter raised to a real positive power,
--   given a mixed-cell configuration and a random coefficient system.

  function Offset ( mic : Mixed_Cell; idx : integer32 )
                  return double_float;

  -- DESCRIPTION :
  --   Returns the value of the offset in the homotopy as the product
  --   of the compenent of the mixed cell with index idx.

  procedure Make_Powers
              ( file : in file_type;
                cfq : in Standard_Complex_VecVecs.VecVec;
                mix : in Standard_Integer_Vectors.Link_to_Vector;
                lif : in Arrays_of_Floating_Vector_Lists.Array_of_Lists;
                mic : in Mixed_Cell; 
                pow : out Standard_Floating_VecVecs.VecVec );

  -- DESCRIPTION :
  --   Defines the powers of the continuation parameter t
  --   in the homotopy for the mixed cell and the lifting,
  --   using the coefficients in the system q.

  -- ON ENTRY :
  --   file     for intermediate output;
  --   cfq      coefficient vectors of the random coefficient system;
  --   mix      type of mixture of the cell;
  --   lif      lifting values for all supports;
  --   mic      mixed cell defines the initial form system.
  --   pow      space for cfq'last pointers to floating-point vectors.

  -- ON RETURN :
  --   pow      powers of t in the homotopy.

  function Smallest_Nonzero_Power
             ( pow : Standard_Floating_VecVecs.VecVec )
             return double_float;

  -- DESCRIPTION :
  --   Returns the smallest nonzero number in pow.

  procedure Scale_Powers
              ( pow : in out Standard_Floating_VecVecs.VecVec;
                val : in double_float );

  -- DESCRIPTION :
  --   Multiplies every nonzero value in pow by 1/val.

  -- REQUIRED : val /= 0.0.

  procedure Make_Homotopy
              ( file : in file_type;
                deg : in integer32;
                pnt : in double_float;
                cfq : in Standard_Complex_VecVecs.VecVec;
                mix : in Standard_Integer_Vectors.Link_to_Vector;
                lif : in Arrays_of_Floating_Vector_Lists.Array_of_Lists;
                pow : in Standard_Floating_VecVecs.VecVec;
                thm : out Taylor_Homotopy );

  -- DESCRIPTION :
  --   Given the powers of the continuation parameter t,
  --   defines the Taylor monomial homotopy.

  -- ON ENTRY :
  --   file     opened for output;
  --   deg      truncation degree of the series;
  --   pnt      expansion point for the Taylor series;
  --   cfq      coefficient vectors of the random coefficient system;
  --   mix      type of mixture of the cell;
  --   lif      lifting values for all supports;
  --   pow      powers of t in the homotopy.

  -- ON RETURN :
  --   thm      the Taylor homotopy for the powers.

  procedure Make_Homotopy
              ( file : in file_type;
                deg : in integer32;
                pnt : in double_float;
                cfq : in Standard_Complex_VecVecs.VecVec;
                mix : in Standard_Integer_Vectors.Link_to_Vector;
                lif : in Arrays_of_Floating_Vector_Lists.Array_of_Lists;
                mic : in Mixed_Cell; thm : out Taylor_Homotopy );

  -- DESCRIPTION :
  --   Defines the homotopy for the mixed cell and the lifting,
  --   using the coefficients in the system q.

  -- ON ENTRY :
  --   file     opened for output;
  --   deg      truncation degree of the series;
  --   pnt      expansion point for the Taylor series;
  --   cfq      coefficient vectors of the random coefficient system;
  --   mix      type of mixture of the cell;
  --   lif      lifting values for all supports;
  --   mic      a mixed cell.

  -- ON RETURN :
  --   thm      the Taylor homotopy for the powers.

end Taylor_Polyhedral_Homotopies;
