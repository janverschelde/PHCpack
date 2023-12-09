with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Double_Double_Numbers;              use Double_Double_Numbers;
with Triple_Double_Numbers;              use Triple_Double_Numbers;
with Quad_Double_Numbers;                use Quad_Double_Numbers;
with Standard_Complex_Numbers;
with DoblDobl_Complex_Numbers;
with TripDobl_Complex_Numbers;
with QuadDobl_Complex_Numbers;
with PentDobl_Complex_Numbers;
with OctoDobl_Complex_Numbers;
with DecaDobl_Complex_Numbers;
with HexaDobl_Complex_Numbers;
with Standard_Complex_Polynomials;
with Standard_Complex_Poly_Systems;
with DoblDobl_Complex_Polynomials;
with DoblDobl_Complex_Poly_Systems;
with TripDobl_Complex_Polynomials;
with TripDobl_Complex_Poly_Systems;
with QuadDobl_Complex_Polynomials;
with QuadDobl_Complex_Poly_Systems;
with PentDobl_Complex_Polynomials;
with PentDobl_Complex_Poly_Systems;
with OctoDobl_Complex_Polynomials;
with OctoDobl_Complex_Poly_Systems;
with DecaDobl_Complex_Polynomials;
with DecaDobl_Complex_Poly_Systems;
with HexaDobl_Complex_Polynomials;
with HexaDobl_Complex_Poly_Systems;
with Standard_CSeries_Polynomials;
with Standard_CSeries_Poly_Systems;
with DoblDobl_CSeries_Polynomials;
with DoblDobl_CSeries_Poly_Systems;
with TripDobl_CSeries_Polynomials;
with TripDobl_CSeries_Poly_Systems;
with QuadDobl_CSeries_Polynomials;
with QuadDobl_CSeries_Poly_Systems;
with PentDobl_CSeries_Polynomials;
with PentDobl_CSeries_Poly_Systems;
with OctoDobl_CSeries_Polynomials;
with OctoDobl_CSeries_Poly_Systems;
with DecaDobl_CSeries_Polynomials;
with DecaDobl_CSeries_Poly_Systems;
with HexaDobl_CSeries_Polynomials;
with HexaDobl_CSeries_Poly_Systems;

package Series_and_Homotopies is

-- DESCRIPTION :
--   A homotopy in one parameter is naturally encoded as a polynomial
--   system with truncated power series as coefficients,
--   in double, double double, triple double, quad double, penta double,
--   octo double, deca double, and hexa double precision.

  function Create ( h : in Standard_Complex_Poly_Systems.Poly_Sys;
                    idx : in integer32; verbose : boolean := false )
                  return Standard_CSeries_Poly_Systems.Poly_Sys;
  function Create ( h : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                    idx : in integer32; verbose : boolean := false )
                  return DoblDobl_CSeries_Poly_Systems.Poly_Sys;
  function Create ( h : in TripDobl_Complex_Poly_Systems.Poly_Sys;
                    idx : in integer32; verbose : boolean := false )
                  return TripDobl_CSeries_Poly_Systems.Poly_Sys;
  function Create ( h : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                    idx : in integer32; verbose : boolean := false )
                  return QuadDobl_CSeries_Poly_Systems.Poly_Sys;
  function Create ( h : in PentDobl_Complex_Poly_Systems.Poly_Sys;
                    idx : in integer32; verbose : boolean := false )
                  return PentDobl_CSeries_Poly_Systems.Poly_Sys;
  function Create ( h : in OctoDobl_Complex_Poly_Systems.Poly_Sys;
                    idx : in integer32; verbose : boolean := false )
                  return OctoDobl_CSeries_Poly_Systems.Poly_Sys;
  function Create ( h : in DecaDobl_Complex_Poly_Systems.Poly_Sys;
                    idx : in integer32; verbose : boolean := false )
                  return DecaDobl_CSeries_Poly_Systems.Poly_Sys;
  function Create ( h : in HexaDobl_Complex_Poly_Systems.Poly_Sys;
                    idx : in integer32; verbose : boolean := false )
                  return HexaDobl_CSeries_Poly_Systems.Poly_Sys;

  -- DESCRIPTION :
  --   Given in h the output of Standard_Homotopy.Homotopy_System
  --   and in idx the index to the homotopy continuation parameter,
  --   on return is polynomial system with coefficients power series
  --   in the homotopy continuation parameter.
  --   Additional information is written to screen, if verbose.

  function Eval ( p : Standard_CSeries_Polynomials.Poly;
                  t : double_float )
                return Standard_Complex_Polynomials.Poly;
  function Eval ( p : Standard_CSeries_Polynomials.Poly;
                  t : Standard_Complex_Numbers.Complex_Number )
                return Standard_Complex_Polynomials.Poly;
  function Eval ( p : DoblDobl_CSeries_Polynomials.Poly;
                  t : double_double )
                return DoblDobl_Complex_Polynomials.Poly;
  function Eval ( p : DoblDobl_CSeries_Polynomials.Poly;
                  t : DoblDobl_Complex_Numbers.Complex_Number )
                return DoblDobl_Complex_Polynomials.Poly;
  function Eval ( p : TripDobl_CSeries_Polynomials.Poly;
                  t : triple_double )
                return TripDobl_Complex_Polynomials.Poly;
  function Eval ( p : TripDobl_CSeries_Polynomials.Poly;
                  t : TripDobl_Complex_Numbers.Complex_Number )
                return TripDobl_Complex_Polynomials.Poly;
  function Eval ( p : QuadDobl_CSeries_Polynomials.Poly;
                  t : quad_double )
                return QuadDobl_Complex_Polynomials.Poly;
  function Eval ( p : QuadDobl_CSeries_Polynomials.Poly;
                  t : QuadDobl_Complex_Numbers.Complex_Number )
                return QuadDobl_Complex_Polynomials.Poly;

  -- DESCRIPTION :
  --   Every series coefficient in p is evaluated at t.
  --   On return is a polynomial with complex coefficients
  --   in double, double double, triple double, or quad double precision.

  function Eval ( h : Standard_CSeries_Poly_Systems.Poly_Sys;
                  t : double_float )
                return Standard_Complex_Poly_Systems.Poly_Sys;
  function Eval ( h : Standard_CSeries_Poly_Systems.Poly_Sys;
                  t : Standard_Complex_Numbers.Complex_Number )
                return Standard_Complex_Poly_Systems.Poly_Sys;
  function Eval ( h : DoblDobl_CSeries_Poly_Systems.Poly_Sys;
                  t : double_double )
                return DoblDobl_Complex_Poly_Systems.Poly_Sys;
  function Eval ( h : DoblDobl_CSeries_Poly_Systems.Poly_Sys;
                  t : DoblDobl_Complex_Numbers.Complex_Number )
                return DoblDobl_Complex_Poly_Systems.Poly_Sys;
  function Eval ( h : TripDobl_CSeries_Poly_Systems.Poly_Sys;
                  t : triple_double )
                return TripDobl_Complex_Poly_Systems.Poly_Sys;
  function Eval ( h : TripDobl_CSeries_Poly_Systems.Poly_Sys;
                  t : TripDobl_Complex_Numbers.Complex_Number )
                return TripDobl_Complex_Poly_Systems.Poly_Sys;
  function Eval ( h : QuadDobl_CSeries_Poly_Systems.Poly_Sys;
                  t : quad_double )
                return QuadDobl_Complex_Poly_Systems.Poly_Sys;
  function Eval ( h : QuadDobl_CSeries_Poly_Systems.Poly_Sys;
                  t : QuadDobl_Complex_Numbers.Complex_Number )
                return QuadDobl_Complex_Poly_Systems.Poly_Sys;

  -- DESCRIPTION :
  --   Every series coefficient in h is evaluated at t.
  --   On return is the evaluated polynomial system,
  --   with complex coefficients in double, double double,
  --   triple double, or quad double precision.

  function Shift ( p : Standard_CSeries_Polynomials.Poly;
                   c : double_float )
                 return Standard_CSeries_Polynomials.Poly;
  function Shift ( p : Standard_CSeries_Polynomials.Poly;
                   c : Standard_Complex_Numbers.Complex_Number )
                 return Standard_CSeries_Polynomials.Poly;
  function Shift ( p : DoblDobl_CSeries_Polynomials.Poly;
                   c : double_double )
                 return DoblDobl_CSeries_Polynomials.Poly;
  function Shift ( p : DoblDobl_CSeries_Polynomials.Poly;
                   c : DoblDobl_Complex_Numbers.Complex_Number )
                 return DoblDobl_CSeries_Polynomials.Poly;
  function Shift ( p : QuadDobl_CSeries_Polynomials.Poly;
                   c : quad_double )
                 return QuadDobl_CSeries_Polynomials.Poly;
  function Shift ( p : QuadDobl_CSeries_Polynomials.Poly;
                   c : QuadDobl_Complex_Numbers.Complex_Number )
                 return QuadDobl_CSeries_Polynomials.Poly;

  -- DESCRIPTION :
  --   The series parameter t in the series coefficient of p is
  --   replaced by t-c, so that: Eval(p,c) = Eval(Shift(p,c),0).

  procedure Shift ( p : in out Standard_CSeries_Polynomials.Poly;
                    c : in double_float );
  procedure Shift ( p : in out Standard_CSeries_Polynomials.Poly;
                    c : in Standard_Complex_Numbers.Complex_Number );
  procedure Shift ( p : in out DoblDobl_CSeries_Polynomials.Poly;
                    c : in double_double );
  procedure Shift ( p : in out DoblDobl_CSeries_Polynomials.Poly;
                    c : in DoblDobl_Complex_Numbers.Complex_Number );
  procedure Shift ( p : in out QuadDobl_CSeries_Polynomials.Poly;
                    c : in quad_double );
  procedure Shift ( p : in out QuadDobl_CSeries_Polynomials.Poly;
                    c : in QuadDobl_Complex_Numbers.Complex_Number );

  -- DESCRIPTION :
  --   On return, p = Shift(p,c).

  function Shift ( p : Standard_CSeries_Poly_Systems.Poly_Sys;
                   c : double_float )
                 return Standard_CSeries_Poly_Systems.Poly_Sys;
  function Shift ( p : Standard_CSeries_Poly_Systems.Poly_Sys;
                   c : Standard_Complex_Numbers.Complex_Number )
                 return Standard_CSeries_Poly_Systems.Poly_Sys;
  function Shift ( p : DoblDobl_CSeries_Poly_Systems.Poly_Sys;
                   c : double_double )
                 return DoblDobl_CSeries_Poly_Systems.Poly_Sys;
  function Shift ( p : DoblDobl_CSeries_Poly_Systems.Poly_Sys;
                   c : DoblDobl_Complex_Numbers.Complex_Number )
                 return DoblDobl_CSeries_Poly_Systems.Poly_Sys;
  function Shift ( p : QuadDobl_CSeries_Poly_Systems.Poly_Sys;
                   c : quad_double )
                 return QuadDobl_CSeries_Poly_Systems.Poly_Sys;
  function Shift ( p : QuadDobl_CSeries_Poly_Systems.Poly_Sys;
                   c : QuadDobl_Complex_Numbers.Complex_Number )
                 return QuadDobl_CSeries_Poly_Systems.Poly_Sys;

  -- DESCRIPTION :
  --   The series parameter t in the series coefficient of p is
  --   replaced by t-c, so that: Eval(p,c) = Eval(Shift(p,c),0).

  procedure Shift ( p : in out Standard_CSeries_Poly_Systems.Poly_Sys;
                    c : in double_float );
  procedure Shift ( p : in out Standard_CSeries_Poly_Systems.Poly_Sys;
                    c : in Standard_Complex_Numbers.Complex_Number );
  procedure Shift ( p : in out DoblDobl_CSeries_Poly_Systems.Poly_Sys;
                    c : in double_double );
  procedure Shift ( p : in out DoblDobl_CSeries_Poly_Systems.Poly_Sys;
                    c : in DoblDobl_Complex_Numbers.Complex_Number );
  procedure Shift ( p : in out QuadDobl_CSeries_Poly_Systems.Poly_Sys;
                    c : in quad_double );
  procedure Shift ( p : in out QuadDobl_CSeries_Poly_Systems.Poly_Sys;
                    c : in QuadDobl_Complex_Numbers.Complex_Number );

  -- DESCRIPTION :
  --   On return, p = Shift(p,c).

end Series_and_Homotopies;
