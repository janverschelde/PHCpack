with Standard_Integer_Numbers;          use Standard_Integer_Numbers;

package Job_Containers is

-- DESCRIPTION :
--   Defines the functions to handle the jobs for the C gateway,
--   in particular those functions that involve the copying of data
--   from and into the containers for systems and solutions.
--   The functions return 0 if all went well, or else the job code.
--   If the verbose level vrbvlv is positive, then the name of the
--   function is written to screen, to track bugs faster.

  function Standard_Copy_Target_System_to_Container
             ( vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Copies the target system to the systems container,
  --   for systems in standard double precision.

  function DoblDobl_Copy_Target_System_to_Container
             ( vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Copies the target system to the systems container,
  --   for systems in double double precision.

  function QuadDobl_Copy_Target_System_to_Container
             ( vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Copies the target system to the systems container,
  --   for systems in quad double precision.

  function Multprec_Copy_Target_System_to_Container
             ( vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Copies the target system to the systems container,
  --   for systems in arbitrary multiprecision.

  function Standard_Copy_Container_System_to_Target
             ( vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Copies the system in the container to the target system,
  --   for systems in standard double precision.

  function DoblDobl_Copy_Container_System_to_Target
             ( vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Copies the system in the container to the target system,
  --   for systems in double double precision.

  function QuadDobl_Copy_Container_System_to_Target
             ( vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Copies the system in the container to the target system,
  --   for systems in quad double precision.

  function Multprec_Copy_Container_System_to_Target
             ( vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Copies the system in the container to the target system,
  --   for systems in arbitrary multiprecision.

  function Standard_Copy_Start_System_to_Container
             ( vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Copies the start system to the systems container,
  --   for systems in standard double precision.

  function DoblDobl_Copy_Start_System_to_Container
             ( vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Copies the start system to the systems container,
  --   for systems in double double precision.

  function QuadDobl_Copy_Start_System_to_Container
             ( vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Copies the start system to the systems container,
  --   for systems in quad double precision.

  function Multprec_Copy_Start_System_to_Container
             ( vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Copies the start system to the systems container,
  --   for systems in arbitrary multiprecision.

  function Standard_Copy_Container_System_to_Start
             ( vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Copies the system in the container to the start system,
  --   for systems in standard double precision.

  function DoblDobl_Copy_Container_System_to_Start
             ( vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Copies the system in the container to the start system,
  --   for systems in double double precision.

  function QuadDobl_Copy_Container_System_to_Start
             ( vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Copies the system in the container to the start system,
  --   for systems in quad double precision.

  function Multprec_Copy_Container_System_to_Start
             ( vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Copies the system in the container to the start system,
  --   for systems in arbitrary multiprecision.

  function Standard_Copy_Target_Solutions_to_Container
             ( vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Copies the target solutions to the solutions container,
  --   for solutions in standard double precision.

  function DoblDobl_Copy_Target_Solutions_to_Container
             ( vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Copies the target solutions to the solutions container,
  --   for solutions in double double precision.

  function QuadDobl_Copy_Target_Solutions_to_Container
             ( vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Copies the target solutions to the solutions container,
  --   for solutions in quad double precision.

  function Multprec_Copy_Target_Solutions_to_Container
             ( vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Copies the target solutions to the solutions container,
  --   for solutions in arbitrary multiprecision.

  function Standard_Copy_Container_Solutions_to_Target
             ( vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Copies the solutions in the container to the target solutions,
  --   for solutions in standard double precision.

  function DoblDobl_Copy_Container_Solutions_to_Target
             ( vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Copies the solutions in the container to the target solutions,
  --   for solutions in double double precision.

  function QuadDobl_Copy_Container_Solutions_to_Target
             ( vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Copies the solutions in the container to the target solutions,
  --   for solutions in quad double precision.

  function Multprec_Copy_Container_Solutions_to_Target
             ( vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Copies the solutions in the container to the target solutions,
  --   for solutions in arbitrary multiprecision.

  function Standard_Copy_Start_Solutions_to_Container
             ( vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Copies the start solutions to the solutions container,
  --   for solutions in standard double precision.

  function DoblDobl_Copy_Start_Solutions_to_Container
             ( vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Copies the start solutions to the solutions container,
  --   for solutions in double double precision.

  function QuadDobl_Copy_Start_Solutions_to_Container
             ( vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Copies the start solutions to the solutions container,
  --   for solutions in quad double precision.

  function Multprec_Copy_Start_Solutions_to_Container
             ( vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Copies the start solutions to the solutions container,
  --   for solutions in arbitrary multiprecision.

  function Standard_Copy_Container_Solutions_to_Start
             ( vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Copies the solutions in the container to the start solutions,
  --   for solutions in standard double precision.

  function DoblDobl_Copy_Container_Solutions_to_Start
             ( vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Copies the solutions in the container to the start solutions,
  --   for solutions in double double precision.

  function QuadDobl_Copy_Container_Solutions_to_Start
             ( vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Copies the solutions in the container to the start solutions,
  --   for solutions in quad double precision.

  function Multprec_Copy_Container_Solutions_to_Start
             ( vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Copies the solutions in the container to the start solutions,
  --   for solutions in arbitrary multiprecision.

end Job_Containers;
