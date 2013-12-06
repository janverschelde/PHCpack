with text_io;                            use text_io;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Irreducible_Component_Lists;        use Irreducible_Component_Lists;

package Drivers_to_Component_Creators is

-- DESCRIPTION :
--   This package contains facilities to tune the parameters and to
--   write the final summary of the breakup of an equidimensional
--   solution set into irreducible components.

  procedure Display_Filter_Settings
               ( file : in file_type; full_output : in boolean;
                 stoptol,membtol : in double_float );

  -- DESCRIPTION :
  --   Displays the current settings of the filter.

  -- ON ENTRY :
  --   file      if Standard_Output, then display will be on screen,
  --             otherwise, the settings are written to the file;
  --   full_output is flag to indicate whether all intermediate diagnostics
  --             of the sampling need to be written, or just a summary;
  --   stoptol   tolerance to stop interpolating in the filter construction;
  --   membtol   tolerance to decide membership of any point to component.

  procedure Standard_Tuner
               ( file : in file_type; full_output : in out boolean;
                 stoptol,membtol : in out double_float );
  procedure Multprec_Tuner
               ( file : in file_type; full_output : in out boolean;
                 size : in natural32; stoptol,membtol : in out double_float );

  -- DESCRIPTION :
  --   Given default values for tolerances of stop and membership test,
  --   this routine allows interactive tuning.

  procedure Write_Summary
               ( file : in file_type; k : in natural32;
                 deco : in Standard_Irreducible_Component_List );
  procedure Write_Summary
               ( file : in file_type; k : in natural32;
                 deco : in Multprec_Irreducible_Component_List );

  -- DESCRIPTION :
  --   Writes a summary of the list of irreducible components on file.

end Drivers_to_Component_Creators;
