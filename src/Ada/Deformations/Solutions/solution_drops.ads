with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Complex_Solutions;
with DoblDobl_Complex_Solutions;
with TripDobl_Complex_Solutions;
with QuadDobl_Complex_Solutions;
with PentDobl_Complex_Solutions;
with OctoDobl_Complex_Solutions;
with DecaDobl_Complex_Solutions;
with HexaDobl_Complex_Solutions;

package Solution_Drops is

-- DESCRIPTION :
--   This package exports operations to drop a coordinate of a solution.

  function Drop ( s : Standard_Complex_Solutions.Solution; k : natural32 )
                return Standard_Complex_Solutions.Solution;
  function Drop ( s : DoblDobl_Complex_Solutions.Solution; k : natural32 )
                return DoblDobl_Complex_Solutions.Solution;
  function Drop ( s : TripDobl_Complex_Solutions.Solution; k : natural32 )
                return TripDobl_Complex_Solutions.Solution;
  function Drop ( s : QuadDobl_Complex_Solutions.Solution; k : natural32 )
                return QuadDobl_Complex_Solutions.Solution;
  function Drop ( s : PentDobl_Complex_Solutions.Solution; k : natural32 )
                return PentDobl_Complex_Solutions.Solution;
  function Drop ( s : OctoDobl_Complex_Solutions.Solution; k : natural32 )
                return OctoDobl_Complex_Solutions.Solution;
  function Drop ( s : DecaDobl_Complex_Solutions.Solution; k : natural32 )
                return DecaDobl_Complex_Solutions.Solution;
  function Drop ( s : HexaDobl_Complex_Solutions.Solution; k : natural32 )
                return HexaDobl_Complex_Solutions.Solution;

  -- DESCRIPTION :
  --   Returns a solution with the k-th coordinate removed.
 
  -- REQUIRED : s.n > 1 and k is in the range 1..s.n.

  function Drop ( s : Standard_Complex_Solutions.Solution_List;
                  k : natural32 )
                return Standard_Complex_Solutions.Solution_List;
  function Drop ( s : DoblDobl_Complex_Solutions.Solution_List;
                  k : natural32 )
                return DoblDobl_Complex_Solutions.Solution_List;
  function Drop ( s : TripDobl_Complex_Solutions.Solution_List;
                  k : natural32 )
                return TripDobl_Complex_Solutions.Solution_List;
  function Drop ( s : QuadDobl_Complex_Solutions.Solution_List;
                  k : natural32 )
                return QuadDobl_Complex_Solutions.Solution_List;
  function Drop ( s : PentDobl_Complex_Solutions.Solution_List;
                  k : natural32 )
                return PentDobl_Complex_Solutions.Solution_List;
  function Drop ( s : OctoDobl_Complex_Solutions.Solution_List;
                  k : natural32 )
                return OctoDobl_Complex_Solutions.Solution_List;
  function Drop ( s : DecaDobl_Complex_Solutions.Solution_List;
                  k : natural32 )
                return DecaDobl_Complex_Solutions.Solution_List;
  function Drop ( s : HexaDobl_Complex_Solutions.Solution_List;
                  k : natural32 )
                return HexaDobl_Complex_Solutions.Solution_List;

  -- DESCRIPTION :
  --   Returns the solutions in s with the k-th coordinate removed.
 
  -- REQUIRED : s.n > 1 and k is in the range 1..s.n.

end Solution_Drops;
