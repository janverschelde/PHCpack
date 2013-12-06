with DoblDobl_Complex_Laur_Systems;
with QuadDobl_Complex_Laur_Systems;
with DoblDobl_Complex_Solutions;
with QuadDobl_Complex_Solutions;

package Drivers_to_dd_qd_Root_Refiners is

-- DESCRIPTION :
--   Gives access to calling the root refiners in double double
--   or quad double complex arithmetic.

  procedure Standard_to_DoblDobl_Complex
              ( p : out DoblDobl_Complex_Laur_Systems.Link_to_Laur_Sys;
                s : out DoblDobl_Complex_Solutions.Solution_List );

  -- DESCRIPTION :
  --   Prompts the user for a system and solution list 
  --   with standard complex coefficients.
  --   On return is the given system converted to double double numbers.

  procedure Multprec_to_DoblDobl_Complex
              ( p : out DoblDobl_Complex_Laur_Systems.Link_to_Laur_Sys;
                s : out DoblDobl_Complex_Solutions.Solution_List );

  -- DESCRIPTION :
  --   Prompts the user for a system and solution list 
  --   with standard complex coefficients.
  --   On return is the given system converted to double double numbers.

  procedure Standard_to_QuadDobl_Complex
              ( p : out QuadDobl_Complex_Laur_Systems.Link_to_Laur_Sys;
                s : out QuadDobl_Complex_Solutions.Solution_List );

  -- DESCRIPTION :
  --   Prompts the user for a system and solution list 
  --   with standard complex coefficients.
  --   On return is the given system converted to double double numbers.

  procedure Multprec_to_QuadDobl_Complex
              ( p : out QuadDobl_Complex_Laur_Systems.Link_to_Laur_Sys;
                s : out QuadDobl_Complex_Solutions.Solution_List );

  -- DESCRIPTION :
  --   Prompts the user for a system and solution list 
  --   with standard complex coefficients.
  --   On return is the given system converted to quad double numbers.

  procedure DD_QD_Root_Refinement;

  -- DESCRIPTION :
  --   Called by mainvali.

end Drivers_to_dd_qd_Root_Refiners;
