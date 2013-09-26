with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Complex_Vectors;           use Standard_Complex_Vectors;

package Polyhedral_Coefficient_Predictors is

-- DESCRIPTION :
--   This package offers predictors for path tracking with an exponential
--   function of the continuation parameter in the polyhedral homotopy.

  procedure Predictor ( s,ds : in out double_float );

  -- DESCRIPTION :
  --   Determines the next value the continuation parameter s := s + ds,
  --   keeping s >= 0.

  -- ON ENTRY :
  --   s        current value for the continuation parameter, s < 0;
  --   ds       value for the step size, increment for s.

  -- ON RETURN :
  --   s        new value for the continuation parameter;
  --   ds       in case s has become zero, ds is likely to be equal to
  --            the previous value of s, taken in absolute value.

  procedure Step_Control
              ( fail : in boolean; s,ds : in out double_float;
                max_ds : in double_float; cnt : in out natural32;
                x,backup : in out Standard_Complex_Vectors.Vector );

  procedure Step_Control
              ( fail : in boolean; s,ds : in out double_float;
                max_ds : in double_float; cnt : in out natural32;
                x,backup,px : in out Standard_Complex_Vectors.Vector );

  -- DESCRIPTION :
  --   Depending on the success of the corrector, this routine determines
  --   the next value for the increment ds and sets the backup value to x,
  --   or in case of fail restores x to the backup value.

  -- ON ENTRY :
  --   fail     outcome of the corrector, if fail then we need to backup
  --            and reduce the step size, otherwise we update the backup
  --            and may increase the step size;
  --   s        current value of the continuation parameter;
  --   ds       current value for the step size;
  --   max_ds   maximum value for the step size;
  --   cnt      countes the number of successive successes of the corrector;
  --   x        current value of the solution vector;
  --   backup   backup value for the solution;
  --   px       previous value for the solution.
 
  -- ON RETURN :
  --   s        if fail, then s is restored to its previous value;
  --   ds       if fail, then ds is decreased, otherwise if the value
  --            for cnt is right, then ds is increased;
  --   cnt      updated counter for consecutive successes;
  --   x        if fail, then x is restored to backup;
  --   backup   if not fail, then backup equals x;
  --   px       if not fail, then px equals the backup value on entry.

  procedure Secant_Predictor 
              ( x : in out Standard_Complex_Vectors.Vector;
                px : in Standard_Complex_Vectors.Vector;
                h : in double_float );

  -- DESCRIPTION :
  --   Predicts the next value for x along the path by adding to x
  --   h times the difference with its previous value.

  -- ON ENTRY :
  --   x        current solution vector;
  --   px       previous solution vector;
  --   h        step size in t.

  -- ON RETURN :
  --   x        predicted solution vector.

end Polyhedral_Coefficient_Predictors;
