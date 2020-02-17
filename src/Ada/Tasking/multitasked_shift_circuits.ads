with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Double_Double_Numbers;              use Double_Double_Numbers;
with Quad_Double_Numbers;                use Quad_Double_Numbers;
with Standard_Speelpenning_Convolutions;
with DoblDobl_Speelpenning_Convolutions;
with QuadDobl_Speelpenning_Convolutions;

package Multitasked_Shift_Circuits is

-- DESCRIPTION :
--   Provides multitasked code to shift coefficients of convolution circuits.

  procedure Standard_Multitasking
              ( nbtasks,deg : in integer32;
                c : in out Standard_Speelpenning_Convolutions.Circuits;
                t : in double_float;
                output : in boolean := false );
  procedure DoblDobl_Multitasking
              ( nbtasks,deg : in integer32;
                c : in out DoblDobl_Speelpenning_Convolutions.Circuits;
                t : in double_double;
                output : in boolean := false );
  procedure QuadDobl_Multitasking
              ( nbtasks,deg : in integer32;
                c : in out QuadDobl_Speelpenning_Convolutions.Circuits;
                t : in quad_double;
                output : in boolean := false );

  -- DESCRIPTION :
  --   Shifts the coefficients in the circuits with multitasking,
  --   in double, double double, or quad double precision.

  -- ON ENTRY :
  --   nbtasks  the number of tasks;
  --   deg      degree of the power series;
  --   c        a sequence of convolution circuits;
  --   t        value to use in the shift; 
  --   output   true if the tasks are verbose,
  --            false if no output during the multitasking.

  -- ON RETURN :
  --   c        circuits with shifted coefficients.

end Multitasked_Shift_Circuits;
