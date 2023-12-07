with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Complex_Vectors;
with Standard_Complex_Solutions;
with DoblDobl_Complex_Vectors;
with DoblDobl_Complex_Solutions;
with TripDobl_Complex_Vectors;
with TripDobl_Complex_Solutions;
with QuadDobl_Complex_Vectors;
with QuadDobl_Complex_Solutions;
with PentDobl_Complex_Vectors;
with PentDobl_Complex_Solutions;
with OctoDobl_Complex_Vectors;
with OctoDobl_Complex_Solutions;
with DecaDobl_Complex_Vectors;
with DecaDobl_Complex_Solutions;
with HexaDobl_Complex_Vectors;
with HexaDobl_Complex_Solutions;
with Standard_Complex_Series_Vectors;
with Standard_Complex_Series_VecVecs;
with DoblDobl_Complex_Series_Vectors;
with DoblDobl_Complex_Series_VecVecs;
with TripDobl_Complex_Series_Vectors;
with TripDobl_Complex_Series_VecVecs;
with QuadDobl_Complex_Series_Vectors;
with QuadDobl_Complex_Series_VecVecs;
with PentDobl_Complex_Series_Vectors;
with PentDobl_Complex_Series_VecVecs;
with OctoDobl_Complex_Series_Vectors;
with OctoDobl_Complex_Series_VecVecs;
with DecaDobl_Complex_Series_Vectors;
with DecaDobl_Complex_Series_VecVecs;
with HexaDobl_Complex_Series_Vectors;
with HexaDobl_Complex_Series_VecVecs;

package Series_and_Solutions is

-- DESCRIPTION :
--   A solution of a polynomial system defines the constant term
--   of a power series.

  function Create ( sol : Standard_Complex_Vectors.Vector;
                    idx : integer32 ) 
                  return Standard_Complex_Series_Vectors.Vector;
  function Create ( sol : DoblDobl_Complex_Vectors.Vector;
                    idx : integer32 ) 
                  return DoblDobl_Complex_Series_Vectors.Vector;
  function Create ( sol : TripDobl_Complex_Vectors.Vector;
                    idx : integer32 ) 
                  return TripDobl_Complex_Series_Vectors.Vector;
  function Create ( sol : QuadDobl_Complex_Vectors.Vector;
                    idx : integer32 ) 
                  return QuadDobl_Complex_Series_Vectors.Vector;
  function Create ( sol : PentDobl_Complex_Vectors.Vector;
                    idx : integer32 ) 
                  return PentDobl_Complex_Series_Vectors.Vector;
  function Create ( sol : OctoDobl_Complex_Vectors.Vector;
                    idx : integer32 ) 
                  return OctoDobl_Complex_Series_Vectors.Vector;
  function Create ( sol : DecaDobl_Complex_Vectors.Vector;
                    idx : integer32 ) 
                  return DecaDobl_Complex_Series_Vectors.Vector;
  function Create ( sol : HexaDobl_Complex_Vectors.Vector;
                    idx : integer32 ) 
                  return HexaDobl_Complex_Series_Vectors.Vector;

  -- DESCRIPTION :
  --   Takes the coordinates of the vector sol and returns a vector of 
  --   series which have the coordinates as their first, constant term.
  --   If idx = 0, then all coordinates of the solution are copied
  --   to the result on return, otherwise, the coordinate with idx
  --   is skipped in the vector of series on return.

  function Create ( sol : Standard_Complex_Solutions.Solution;
                    idx : integer32 ) 
                  return Standard_Complex_Series_Vectors.Vector;
  function Create ( sol : DoblDobl_Complex_Solutions.Solution;
                    idx : integer32 ) 
                  return DoblDobl_Complex_Series_Vectors.Vector;
  function Create ( sol : TripDobl_Complex_Solutions.Solution;
                    idx : integer32 ) 
                  return TripDobl_Complex_Series_Vectors.Vector;
  function Create ( sol : QuadDobl_Complex_Solutions.Solution;
                    idx : integer32 ) 
                  return QuadDobl_Complex_Series_Vectors.Vector;
  function Create ( sol : PentDobl_Complex_Solutions.Solution;
                    idx : integer32 ) 
                  return PentDobl_Complex_Series_Vectors.Vector;
  function Create ( sol : OctoDobl_Complex_Solutions.Solution;
                    idx : integer32 ) 
                  return OctoDobl_Complex_Series_Vectors.Vector;
  function Create ( sol : DecaDobl_Complex_Solutions.Solution;
                    idx : integer32 ) 
                  return DecaDobl_Complex_Series_Vectors.Vector;
  function Create ( sol : HexaDobl_Complex_Solutions.Solution;
                    idx : integer32 ) 
                  return HexaDobl_Complex_Series_Vectors.Vector;

  -- DESCRIPTION :
  --   Takes the coordinates of the solution and returns a vector of 
  --   series which have the coordinates as their first, constant term.
  --   If idx = 0, then all coordinates of the solution are copied
  --   to the result on return, otherwise, the coordinate with idx
  --   is skipped in the vector of series on return.

  function Create ( sols : Standard_Complex_Solutions.Solution_List;
                    idx : integer32 ) 
                  return Standard_Complex_Series_VecVecs.VecVec;
  function Create ( sols : DoblDobl_Complex_Solutions.Solution_List;
                    idx : integer32 ) 
                  return DoblDobl_Complex_Series_VecVecs.VecVec;
  function Create ( sols : TripDobl_Complex_Solutions.Solution_List;
                    idx : integer32 ) 
                  return TripDobl_Complex_Series_VecVecs.VecVec;
  function Create ( sols : QuadDobl_Complex_Solutions.Solution_List;
                    idx : integer32 ) 
                  return QuadDobl_Complex_Series_VecVecs.VecVec;
  function Create ( sols : PentDobl_Complex_Solutions.Solution_List;
                    idx : integer32 ) 
                  return PentDobl_Complex_Series_VecVecs.VecVec;
  function Create ( sols : OctoDobl_Complex_Solutions.Solution_List;
                    idx : integer32 ) 
                  return OctoDobl_Complex_Series_VecVecs.VecVec;
  function Create ( sols : DecaDobl_Complex_Solutions.Solution_List;
                    idx : integer32 ) 
                  return DecaDobl_Complex_Series_VecVecs.VecVec;
  function Create ( sols : HexaDobl_Complex_Solutions.Solution_List;
                    idx : integer32 ) 
                  return HexaDobl_Complex_Series_VecVecs.VecVec;

  -- DESCRIPTION :
  --   Returns a vector of vectors of series.  The series are defined
  --   by the coordinates of the solutions in the list.
  --   The range of the vecvec on return is 1..Length_Of(sols).
  --   The coordinate with the given index idx is skipped.

end Series_and_Solutions;
