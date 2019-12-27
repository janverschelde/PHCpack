with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Complex_Polynomials;
with Standard_Complex_Poly_Systems;
with DoblDobl_Complex_Polynomials;
with DoblDobl_Complex_Poly_Systems;
with QuadDobl_Complex_Polynomials;
with QuadDobl_Complex_Poly_Systems;
with Standard_CSeries_Polynomials;
with Standard_CSeries_Poly_Systems;
with DoblDobl_CSeries_Polynomials;
with DoblDobl_CSeries_Poly_Systems;
with QuadDobl_CSeries_Polynomials;
with QuadDobl_CSeries_Poly_Systems;
with Standard_Speelpenning_Convolutions;
with DoblDobl_Speelpenning_Convolutions;
with QuadDobl_Speelpenning_Convolutions;

package System_Convolution_Circuits is

-- DESCRIPTION :
--   The functions in this package return convolution circuits,
--   with polynomials and systems of polynomials on input.

  function Make_Convolution_Circuit
             ( p : Standard_Complex_Polynomials.Poly;
               d : natural32 )
             return Standard_Speelpenning_Convolutions.Convolution_Circuit;
  function Make_Convolution_Circuit
             ( p : DoblDobl_Complex_Polynomials.Poly;
               d : natural32 )
             return DoblDobl_Speelpenning_Convolutions.Convolution_Circuit;
  function Make_Convolution_Circuit
             ( p : QuadDobl_Complex_Polynomials.Poly;
               d : natural32 )
             return QuadDobl_Speelpenning_Convolutions.Convolution_Circuit;

  -- DESCRIPTION :
  --   Returns the convolution circuit representation of the polynomial p,
  --   in double, double double, or quad double precision.
  --   The degree of the series is given by the natural number d.

  function Make_Convolution_Circuit
             ( p : Standard_CSeries_Polynomials.Poly )
             return Standard_Speelpenning_Convolutions.Convolution_Circuit;
  function Make_Convolution_Circuit
             ( p : DoblDobl_CSeries_Polynomials.Poly )
             return DoblDobl_Speelpenning_Convolutions.Convolution_Circuit;
  function Make_Convolution_Circuit
             ( p : QuadDobl_CSeries_Polynomials.Poly )
             return QuadDobl_Speelpenning_Convolutions.Convolution_Circuit;

  -- DESCRIPTION :
  --   Returns the convolution circuit representation of the polynomial p,
  --   in double, double double, or quad double precision.

  function Make_Convolution_Circuits
             ( p : Standard_Complex_Poly_Systems.Poly_Sys;
               d : natural32 )
             return Standard_Speelpenning_Convolutions.Convolution_Circuits;
  function Make_Convolution_Circuits
             ( p : DoblDobl_Complex_Poly_Systems.Poly_Sys;
               d : natural32 )
             return DoblDobl_Speelpenning_Convolutions.Convolution_Circuits;
  function Make_Convolution_Circuits
             ( p : QuadDobl_Complex_Poly_Systems.Poly_Sys;
               d : natural32 )
             return QuadDobl_Speelpenning_Convolutions.Convolution_Circuits;

  -- DESCRIPTION :
  --   Returns the convolution circuits representation of the polynomials in p,
  --   in double, double double, or quad double precision.
  --   The degree of the series is given by the natural number d.

  function Make_Convolution_Circuits
             ( p : Standard_CSeries_Poly_Systems.Poly_Sys )
             return Standard_Speelpenning_Convolutions.Convolution_Circuits;
  function Make_Convolution_Circuits
             ( p : DoblDobl_CSeries_Poly_Systems.Poly_Sys )
             return DoblDobl_Speelpenning_Convolutions.Convolution_Circuits;
  function Make_Convolution_Circuits
             ( p : QuadDobl_CSeries_Poly_Systems.Poly_Sys )
             return QuadDobl_Speelpenning_Convolutions.Convolution_Circuits;

  -- DESCRIPTION :
  --   Returns the convolution circuits representation of the polynomials in p,
  --   in double, double double, or quad double precision.

end System_Convolution_Circuits;
