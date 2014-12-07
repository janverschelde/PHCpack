with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Complex_Solutions;
with DoblDobl_Complex_Solutions;
with QuadDobl_Complex_Solutions;
with Multprec_Complex_Solutions;

package Varbprec_Complex_Solutions is

-- DESCRIPTION :
--   Offers conversions of solutions and lists of solutions
--   over various levels of precision.

  function Multprec_to_Standard_Solution
             ( s : Multprec_Complex_Solutions.Solution )
             return Standard_Complex_Solutions.Solution;
  function Multprec_to_Standard_Solutions
             ( s : Multprec_Complex_Solutions.Solution_List )
             return Standard_Complex_Solutions.Solution_List;

  -- DESCRIPTION :
  --   Converts the multiprecision solutions into standard precision.

  function Multprec_to_DoblDobl_Solution
             ( s : Multprec_Complex_Solutions.Solution )
             return DoblDobl_Complex_Solutions.Solution;
  function Multprec_to_DoblDobl_Solutions
             ( s : Multprec_Complex_Solutions.Solution_List )
             return DoblDobl_Complex_Solutions.Solution_List;

  -- DESCRIPTION :
  --   Converts the multiprecision solutions into double double precision.

  function Multprec_to_QuadDobl_Solution
             ( s : Multprec_Complex_Solutions.Solution )
             return QuadDobl_Complex_Solutions.Solution;
  function Multprec_to_QuadDobl_Solutions
             ( s : Multprec_Complex_Solutions.Solution_List )
             return QuadDobl_Complex_Solutions.Solution_List;

  -- DESCRIPTION :
  --   Converts the multiprecision solutions into quad double precision.

  function Standard_to_Multprec_Solution
             ( s : Standard_Complex_Solutions.Solution; size : natural32 )
             return Multprec_Complex_Solutions.Solution;
  function Standard_to_Multprec_Solutions
             ( s : Standard_Complex_Solutions.Solution_List; size : natural32 )
             return Multprec_Complex_Solutions.Solution_List;

  -- DESCRIPTION :
  --   Converts the solutions from standard to arbitrary multiprecision,
  --   using size as the size of the multiprecision numbers;

  function DoblDobl_to_Multprec_Solution
             ( s : DoblDobl_Complex_Solutions.Solution; size : natural32 )
             return Multprec_Complex_Solutions.Solution;
  function DoblDobl_to_Multprec_Solutions
             ( s : DoblDobl_Complex_Solutions.Solution_List; size : natural32 )
             return Multprec_Complex_Solutions.Solution_List;

  -- DESCRIPTION :
  --   Converts the solutions from double double to arbitrary multiprecision,
  --   using size as the size of the multiprecision numbers;

  function QuadDobl_to_Multprec_Solution
             ( s : QuadDobl_Complex_Solutions.Solution; size : natural32 )
             return Multprec_Complex_Solutions.Solution;
  function QuadDobl_to_Multprec_Solutions
             ( s : QuadDobl_Complex_Solutions.Solution_List; size : natural32 )
             return Multprec_Complex_Solutions.Solution_List;

  -- DESCRIPTION :
  --   Converts the solutions from quad double to arbitrary multiprecision,
  --   using size as the size of the multiprecision numbers;

end Varbprec_Complex_Solutions;
