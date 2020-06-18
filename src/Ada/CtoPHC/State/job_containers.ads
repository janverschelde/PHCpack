with Standard_Integer_Numbers;          use Standard_Integer_Numbers;

package Job_Containers is

-- DESCRIPTION :
--   Defines the functions to handle the jobs for the C gateway,
--   in particular those functions that involve the copying of data
--   from and into the containers for systems and solutions.
--   The functions return 0 if all went well, or else the job code.
--   If the verbose level vrbvlv is positive, then the name of the
--   function is written to screen, to track bugs faster.

-- COPYING POLYNOMIAL SYSTEMS :

  function Standard_Target_Poly_System_to_Container
             ( vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Copies the target system to the systems container,
  --   for polynomial systems in standard double precision.

  function DoblDobl_Target_Poly_System_to_Container
             ( vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Copies the target system to the systems container,
  --   for polynomial systems in double double precision.

  function QuadDobl_Target_Poly_System_to_Container
             ( vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Copies the target system to the systems container,
  --   for polynomial systems in quad double precision.

  function Multprec_Target_Poly_System_to_Container
             ( vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Copies the target system to the systems container,
  --   for polynomial systems in arbitrary multiprecision.

  function Standard_Container_Poly_System_to_Target
             ( vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Copies the system in the container to the target system,
  --   for polynomial systems in standard double precision.

  function DoblDobl_Container_Poly_System_to_Target
             ( vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Copies the system in the container to the target system,
  --   for polynomial systems in double double precision.

  function QuadDobl_Container_Poly_System_to_Target
             ( vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Copies the system in the container to the target system,
  --   for polynomial systems in quad double precision.

  function Multprec_Container_Poly_System_to_Target
             ( vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Copies the system in the container to the target system,
  --   for polynomial systems in arbitrary multiprecision.

  function Standard_Start_Poly_System_to_Container
             ( vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Copies the start system to the systems container,
  --   for polynomial systems in standard double precision.

  function DoblDobl_Start_Poly_System_to_Container
             ( vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Copies the start system to the systems container,
  --   for polynomial systems in double double precision.

  function QuadDobl_Start_Poly_System_to_Container
             ( vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Copies the start system to the systems container,
  --   for polynomial systems in quad double precision.

  function Multprec_Start_Poly_System_to_Container
             ( vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Copies the start system to the systems container,
  --   for polynomial systems in arbitrary multiprecision.

  function Standard_Container_Poly_System_to_Start
             ( vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Copies the system in the container to the start system,
  --   for polynomial systems in standard double precision.

  function DoblDobl_Container_Poly_System_to_Start
             ( vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Copies the system in the container to the start system,
  --   for polynomial systems in double double precision.

  function QuadDobl_Container_Poly_System_to_Start
             ( vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Copies the system in the container to the start system,
  --   for polynomial systems in quad double precision.

  function Multprec_Container_Poly_System_to_Start
             ( vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Copies the system in the container to the start system,
  --   for polynomial systems in arbitrary multiprecision.

-- COPYING LAURENT POLYNOMIAL SYSTEMS :

  function Standard_Target_Laur_System_to_Container
             ( vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Copies the target system to the systems container,
  --   for Laurent polynomial systems in standard double precision.

  function DoblDobl_Target_Laur_System_to_Container
             ( vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Copies the target system to the systems container,
  --   for Laurent polynomial systems in double double precision.

  function QuadDobl_Target_Laur_System_to_Container
             ( vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Copies the target system to the systems container,
  --   for Laurent polynomial systems in quad double precision.

  function Standard_Container_Laur_System_to_Target
             ( vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Copies the system in the container to the target system,
  --   for Laurent polynomial systems in standard double precision.

  function DoblDobl_Container_Laur_System_to_Target
             ( vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Copies the system in the container to the target system,
  --   for Laurent polynomial systems in double double precision.

  function QuadDobl_Container_Laur_System_to_Target
             ( vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Copies the system in the container to the target system,
  --   for Laurent polynomial systems in quad double precision.

  function Standard_Start_Laur_System_to_Container
             ( vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Copies the start system to the systems container,
  --   for Laurent polynomial systems in standard double precision.

  function DoblDobl_Start_Laur_System_to_Container
             ( vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Copies the start system to the systems container,
  --   for Laurent polynomial systems in double double precision.

  function QuadDobl_Start_Laur_System_to_Container
             ( vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Copies the start system to the systems container,
  --   for Laurent polynomial systems in quad double precision.

  function Standard_Container_Laur_System_to_Start
             ( vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Copies the system in the container to the start system,
  --   for Laurent polynomial systems in standard double precision.

  function DoblDobl_Container_Laur_System_to_Start
             ( vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Copies the system in the container to the start system,
  --   for Laurent polynomial systems in double double precision.

  function QuadDobl_Container_Laur_System_to_Start
             ( vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Copies the system in the container to the start system,
  --   for Laurent polynomial systems in quad double precision.

-- COPYING SOLUTIONS :

  function Standard_Target_Solutions_to_Container
             ( vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Copies the target solutions to the solutions container,
  --   for solutions in standard double precision.

  function DoblDobl_Target_Solutions_to_Container
             ( vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Copies the target solutions to the solutions container,
  --   for solutions in double double precision.

  function QuadDobl_Target_Solutions_to_Container
             ( vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Copies the target solutions to the solutions container,
  --   for solutions in quad double precision.

  function Multprec_Target_Solutions_to_Container
             ( vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Copies the target solutions to the solutions container,
  --   for solutions in arbitrary multiprecision.

  function Standard_Container_Solutions_to_Target
             ( vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Copies the solutions in the container to the target solutions,
  --   for solutions in standard double precision.

  function DoblDobl_Container_Solutions_to_Target
             ( vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Copies the solutions in the container to the target solutions,
  --   for solutions in double double precision.

  function QuadDobl_Container_Solutions_to_Target
             ( vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Copies the solutions in the container to the target solutions,
  --   for solutions in quad double precision.

  function Multprec_Container_Solutions_to_Target
             ( vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Copies the solutions in the container to the target solutions,
  --   for solutions in arbitrary multiprecision.

  function Standard_Start_Solutions_to_Container
             ( vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Copies the start solutions to the solutions container,
  --   for solutions in standard double precision.

  function DoblDobl_Start_Solutions_to_Container
             ( vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Copies the start solutions to the solutions container,
  --   for solutions in double double precision.

  function QuadDobl_Start_Solutions_to_Container
             ( vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Copies the start solutions to the solutions container,
  --   for solutions in quad double precision.

  function Multprec_Start_Solutions_to_Container
             ( vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Copies the start solutions to the solutions container,
  --   for solutions in arbitrary multiprecision.

  function Standard_Container_Solutions_to_Start
             ( vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Copies the solutions in the container to the start solutions,
  --   for solutions in standard double precision.

  function DoblDobl_Container_Solutions_to_Start
             ( vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Copies the solutions in the container to the start solutions,
  --   for solutions in double double precision.

  function QuadDobl_Container_Solutions_to_Start
             ( vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Copies the solutions in the container to the start solutions,
  --   for solutions in quad double precision.

  function Multprec_Container_Solutions_to_Start
             ( vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Copies the solutions in the container to the start solutions,
  --   for solutions in arbitrary multiprecision.

end Job_Containers;
