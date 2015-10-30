with text_io;                           use text_io;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Complex_Numbers;
with DoblDobl_Complex_Numbers;
with QuadDobl_Complex_Numbers;
with Standard_Integer_Vectors;
with Standard_Complex_Vectors;
with DoblDobl_Complex_Vectors;
with QuadDobl_Complex_Vectors;
with Standard_Complex_Poly_Systems;
with DoblDobl_Complex_Poly_Systems;
with DoblDobl_Complex_Jaco_Matrices;
with QuadDobl_Complex_Poly_Systems;
with Standard_Complex_Solutions;
with DoblDobl_Complex_Solutions;
with QuadDobl_Complex_Solutions;

package Complex_Convex_Continuation is

-- DESCRIPTION :
--   A special case of Coefficient Parameter Homotopy Continuation 
--   is when the parameter space is convex and we can move from one
--   set of parameter values to another via convex-linear complex arcs,
--   made generic by the gamma trick, so singularities can be avoided.
--   The definition of start and target values for the parameters
--   leads then naturally to interpolation between the two.
--   This interpolation may be seen as the default function to instantiate
--   the general parameter continuation routine with.
--   For efficiency -- to reduce the number of variables -- interpolating
--   between start and target is not expressed by explicit equations,
--   but via direct substitution in the increment-and-fix continuation.

  function Interpolate ( a,b : Standard_Complex_Vectors.Vector;
                         t : Standard_Complex_Numbers.Complex_Number )
                       return Standard_Complex_Vectors.Vector;
  function Interpolate ( a,b : DoblDobl_Complex_Vectors.Vector;
                         t : DoblDobl_Complex_Numbers.Complex_Number )
                       return DoblDobl_Complex_Vectors.Vector;
  function Interpolate ( a,b : QuadDobl_Complex_Vectors.Vector;
                         t : QuadDobl_Complex_Numbers.Complex_Number )
                       return QuadDobl_Complex_Vectors.Vector;

  -- DESCRIPTION :
  --   Returns (1-t)*a + t*b, for t between 0 and 1.  Is the 1st
  --   default to use in the instantiation of Parameter_Continuation.

  function Circulate ( a,b : Standard_Complex_Vectors.Vector;
                       gamma,t : Standard_Complex_Numbers.Complex_Number )
                     return Standard_Complex_Vectors.Vector;
  function Circulate ( a,b : DoblDobl_Complex_Vectors.Vector;
                       gamma,t : DoblDobl_Complex_Numbers.Complex_Number )
                     return DoblDobl_Complex_Vectors.Vector;
  function Circulate ( a,b : QuadDobl_Complex_Vectors.Vector;
                       gamma,t : QuadDobl_Complex_Numbers.Complex_Number )
                     return QuadDobl_Complex_Vectors.Vector;

  -- DESCRIPTION :
  --   Returns (1-s)*a + s*b, for s = t + gamma*t*(1-t).  This move
  --   through the parameter space ensures complex arithmetic,
  --   in standard double, double double, or quad double arithmetic.

  function Differentiate ( a,b : Standard_Complex_Vectors.Vector )
                         return Standard_Complex_Vectors.Vector;
  function Differentiate ( a,b : DoblDobl_Complex_Vectors.Vector )
                         return DoblDobl_Complex_Vectors.Vector;
  function Differentiate ( a,b : QuadDobl_Complex_Vectors.Vector )
                         return QuadDobl_Complex_Vectors.Vector;

  -- DESCRIPTION :
  --   Just returns b-a, in standard double, double double, or quad
  --   double precision, as this is the derivative of the Interpolate
  --   function of the parameters.  Is the 2nd default to use in the
  --   instantiation of Parameter_Continuation.

  generic
    with function Evaluate_Parameters 
                    ( t : Standard_Complex_Numbers.Complex_Number )
                    return Standard_Complex_Vectors.Vector;
    -- returns value of parameters at t
    with function Differentiate_Parameters
                    ( t : Standard_Complex_Numbers.Complex_Number )
                    return Standard_Complex_Vectors.Vector;
    -- returns derivatives of parameters with respect to t
  procedure Standard_Reporting_Parameter_Continuation
              ( file : in file_type;
                n : in integer32;
                p : in Standard_Complex_Poly_Systems.Poly_Sys;
                pars : in Standard_Integer_Vectors.Vector;
                vars : in Standard_Integer_Vectors.Vector;
		sols : in out Standard_Complex_Solutions.Solution_List;
                output : in boolean );
  generic
    with function Evaluate_Parameters 
                    ( t : DoblDobl_Complex_Numbers.Complex_Number )
                    return DoblDobl_Complex_Vectors.Vector;
    -- returns value of parameters at t
    with function Differentiate_Parameters
                    ( t : DoblDobl_Complex_Numbers.Complex_Number )
                    return DoblDobl_Complex_Vectors.Vector;
    -- returns derivatives of parameters with respect to t
  procedure DoblDobl_Reporting_Parameter_Continuation
              ( file : in file_type;
                n : in integer32;
                p : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                pars : in Standard_Integer_Vectors.Vector;
                vars : in Standard_Integer_Vectors.Vector;
		sols : in out DoblDobl_Complex_Solutions.Solution_List;
                output : in boolean );
  generic
    with function Evaluate_Parameters 
                    ( t : QuadDobl_Complex_Numbers.Complex_Number )
                    return QuadDobl_Complex_Vectors.Vector;
    -- returns value of parameters at t
    with function Differentiate_Parameters
                    ( t : QuadDobl_Complex_Numbers.Complex_Number )
                    return QuadDobl_Complex_Vectors.Vector;
    -- returns derivatives of parameters with respect to t
  procedure QuadDobl_Reporting_Parameter_Continuation
              ( file : in file_type;
                n : in integer32;
                p : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                pars : in Standard_Integer_Vectors.Vector;
                vars : in Standard_Integer_Vectors.Vector;
		sols : in out QuadDobl_Complex_Solutions.Solution_List;
                output : in boolean );

  -- DESCRIPTION :
  --   Executes the parameter continuation for the homotopy defined by
  --   the system p and for parameters defined by the generics functions,
  --   in standard double, double double, or quad double precision.
  --   The reporting versions allow for the writing of diagnostics
  --   to file, depending on the setting of the output code level.

  -- ON ENTRY :
  --   file     for intermediate output and diagnostics;
  --   n        number of unknowns and parameters in p.
  --   p        polynomial system with parameters;
  --   pars     indices to the parameter variables in p;
  --   vars     indices to the unknown variables in p;
  --   sols     start solutions for start values of parameters;
  --   output   true if intermediate output during continuation.

  -- ON RETURN :
  --   sols     solutions for target values of parameters.

  generic
    with function Evaluate_Parameters 
                    ( t : Standard_Complex_Numbers.Complex_Number )
                    return Standard_Complex_Vectors.Vector;
    -- returns value of parameters at t
    with function Differentiate_Parameters
                    ( t : Standard_Complex_Numbers.Complex_Number )
                    return Standard_Complex_Vectors.Vector;
    -- returns derivatives of parameters with respect to t
  procedure Standard_Silent_Parameter_Continuation
              ( n : in integer32;
                p : in Standard_Complex_Poly_Systems.Poly_Sys;
                pars : in Standard_Integer_Vectors.Vector;
                vars : in Standard_Integer_Vectors.Vector;
		sols : in out Standard_Complex_Solutions.Solution_List );
  generic
    with function Evaluate_Parameters 
                    ( t : DoblDobl_Complex_Numbers.Complex_Number )
                    return DoblDobl_Complex_Vectors.Vector;
    -- returns value of parameters at t
    with function Differentiate_Parameters
                    ( t : DoblDobl_Complex_Numbers.Complex_Number )
                    return DoblDobl_Complex_Vectors.Vector;
    -- returns derivatives of parameters with respect to t
  procedure DoblDobl_Silent_Parameter_Continuation
              ( n : in integer32;
                p : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                pars : in Standard_Integer_Vectors.Vector;
                vars : in Standard_Integer_Vectors.Vector;
		sols : in out DoblDobl_Complex_Solutions.Solution_List );
  generic
    with function Evaluate_Parameters 
                    ( t : QuadDobl_Complex_Numbers.Complex_Number )
                    return QuadDobl_Complex_Vectors.Vector;
    -- returns value of parameters at t
    with function Differentiate_Parameters
                    ( t : QuadDobl_Complex_Numbers.Complex_Number )
                    return QuadDobl_Complex_Vectors.Vector;
    -- returns derivatives of parameters with respect to t
  procedure QuadDobl_Silent_Parameter_Continuation
              ( n : in integer32;
                p : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                pars : in Standard_Integer_Vectors.Vector;
                vars : in Standard_Integer_Vectors.Vector;
		sols : in out QuadDobl_Complex_Solutions.Solution_List );

  -- DESCRIPTION :
  --   Executes the parameter continuation for the homotopy defined by
  --   the system p and for parameters defined by the generics functions,
  --   in standard double, double double, or quad double precision.
  --   The silent versions do not write any output to file or screen.

  -- ON ENTRY :
  --   n        number of unknowns and parameters in p.
  --   p        polynomial system with parameters;
  --   pars     indices to the parameter variables in p;
  --   vars     indices to the unknown variables in p;
  --   sols     start solutions for start values of parameters;
  --   output   true if intermediate output during continuation.

  -- ON RETURN :
  --   sols     solutions for target values of parameters.

end Complex_Convex_Continuation;
