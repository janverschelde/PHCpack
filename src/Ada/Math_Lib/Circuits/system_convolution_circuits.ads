with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Complex_Polynomials;
with Standard_Complex_Poly_Systems;
with DoblDobl_Complex_Polynomials;
with DoblDobl_Complex_Poly_Systems;
with TripDobl_Complex_Polynomials;
with TripDobl_Complex_Poly_Systems;
with QuadDobl_Complex_Polynomials;
with QuadDobl_Complex_Poly_Systems;
with Standard_CSeries_Polynomials;
with Standard_CSeries_Poly_Systems;
with DoblDobl_CSeries_Polynomials;
with DoblDobl_CSeries_Poly_Systems;
with TripDobl_CSeries_Polynomials;
with TripDobl_CSeries_Poly_Systems;
with QuadDobl_CSeries_Polynomials;
with QuadDobl_CSeries_Poly_Systems;
with Standard_Speelpenning_Convolutions;
with DoblDobl_Speelpenning_Convolutions;
with TripDobl_Speelpenning_Convolutions;
with QuadDobl_Speelpenning_Convolutions;

package System_Convolution_Circuits is

-- DESCRIPTION :
--   The functions in this package return convolution circuits,
--   with polynomials and systems of polynomials on input.

  function Make_Convolution_Circuit
             ( p : Standard_Complex_Polynomials.Poly;
               d : natural32 )
             return Standard_Speelpenning_Convolutions.Circuit;
  function Make_Convolution_Circuit
             ( p : DoblDobl_Complex_Polynomials.Poly;
               d : natural32 )
             return DoblDobl_Speelpenning_Convolutions.Circuit;
  function Make_Convolution_Circuit
             ( p : TripDobl_Complex_Polynomials.Poly;
               d : natural32 )
             return TripDobl_Speelpenning_Convolutions.Circuit;
  function Make_Convolution_Circuit
             ( p : QuadDobl_Complex_Polynomials.Poly;
               d : natural32 )
             return QuadDobl_Speelpenning_Convolutions.Circuit;

  -- DESCRIPTION :
  --   Returns the convolution circuit representation of the polynomial p,
  --   in double, double double, triple double, or quad double precision.
  --   The degree of the series is given by the natural number d.

  function Make_Convolution_Circuit
             ( p : Standard_CSeries_Polynomials.Poly )
             return Standard_Speelpenning_Convolutions.Circuit;
  function Make_Convolution_Circuit
             ( p : DoblDobl_CSeries_Polynomials.Poly )
             return DoblDobl_Speelpenning_Convolutions.Circuit;
  function Make_Convolution_Circuit
             ( p : TripDobl_CSeries_Polynomials.Poly )
             return TripDobl_Speelpenning_Convolutions.Circuit;
  function Make_Convolution_Circuit
             ( p : QuadDobl_CSeries_Polynomials.Poly )
             return QuadDobl_Speelpenning_Convolutions.Circuit;

  -- DESCRIPTION :
  --   Returns the convolution circuit representation of the polynomial p,
  --   in double, double double, triple double, or quad double precision.

  function Make_Convolution_Circuits
             ( p : Standard_Complex_Poly_Systems.Poly_Sys;
               d : natural32 )
             return Standard_Speelpenning_Convolutions.Circuits;
  function Make_Convolution_Circuits
             ( p : DoblDobl_Complex_Poly_Systems.Poly_Sys;
               d : natural32 )
             return DoblDobl_Speelpenning_Convolutions.Circuits;
  function Make_Convolution_Circuits
             ( p : TripDobl_Complex_Poly_Systems.Poly_Sys;
               d : natural32 )
             return TripDobl_Speelpenning_Convolutions.Circuits;
  function Make_Convolution_Circuits
             ( p : QuadDobl_Complex_Poly_Systems.Poly_Sys;
               d : natural32 )
             return QuadDobl_Speelpenning_Convolutions.Circuits;

  -- DESCRIPTION :
  --   Returns the convolution circuits representation of the polynomials in p,
  --   in double, double double, triple double, or quad double precision.
  --   The degree of the series is given by the natural number d.

  function Make_Convolution_Circuits
             ( p : Standard_CSeries_Poly_Systems.Poly_Sys )
             return Standard_Speelpenning_Convolutions.Circuits;
  function Make_Convolution_Circuits
             ( p : DoblDobl_CSeries_Poly_Systems.Poly_Sys )
             return DoblDobl_Speelpenning_Convolutions.Circuits;
  function Make_Convolution_Circuits
             ( p : TripDobl_CSeries_Poly_Systems.Poly_Sys )
             return TripDobl_Speelpenning_Convolutions.Circuits;
  function Make_Convolution_Circuits
             ( p : QuadDobl_CSeries_Poly_Systems.Poly_Sys )
             return QuadDobl_Speelpenning_Convolutions.Circuits;

  -- DESCRIPTION :
  --   Returns the convolution circuits representation of the polynomials in p,
  --   in double, double double, triple double, or quad double precision.

  function Make_Convolution_System
             ( p : Standard_CSeries_Poly_Systems.Poly_Sys;
               d : natural32 )
             return Standard_Speelpenning_Convolutions.System;
  function Make_Convolution_System
             ( p : Standard_CSeries_Poly_Systems.Poly_Sys;
               d : natural32 )
             return Standard_Speelpenning_Convolutions.Link_to_System;
  function Make_Convolution_System
             ( p : DoblDobl_CSeries_Poly_Systems.Poly_Sys;
               d : natural32 )
             return DoblDobl_Speelpenning_Convolutions.System;
  function Make_Convolution_System
             ( p : DoblDobl_CSeries_Poly_Systems.Poly_Sys;
               d : natural32 )
             return DoblDobl_Speelpenning_Convolutions.Link_to_System;
  function Make_Convolution_System
             ( p : TripDobl_CSeries_Poly_Systems.Poly_Sys;
               d : natural32 )
             return TripDobl_Speelpenning_Convolutions.System;
  function Make_Convolution_System
             ( p : TripDobl_CSeries_Poly_Systems.Poly_Sys;
               d : natural32 )
             return TripDobl_Speelpenning_Convolutions.Link_to_System;
  function Make_Convolution_System
             ( p : QuadDobl_CSeries_Poly_Systems.Poly_Sys;
               d : natural32 )
             return QuadDobl_Speelpenning_Convolutions.System;
  function Make_Convolution_System
             ( p : QuadDobl_CSeries_Poly_Systems.Poly_Sys;
               d : natural32 )
             return QuadDobl_Speelpenning_Convolutions.Link_to_System;

  -- DESCRIPTION :
  --   Wraps the making of the convolution circuits of the polynomials in p,
  --   allocating all work space, for degree d of the power series,
  --   in double, double double, triple_double, or quad double precision.

  function to_double
	     ( c : DoblDobl_Speelpenning_Convolutions.Circuit )
	     return Standard_Speelpenning_Convolutions.Circuit;
  function to_double
	     ( c : QuadDobl_Speelpenning_Convolutions.Circuit )
	     return Standard_Speelpenning_Convolutions.Circuit;
  function to_double
	     ( c : DoblDobl_Speelpenning_Convolutions.Link_to_Circuit )
             return Standard_Speelpenning_Convolutions.Link_to_Circuit;
  function to_double
	     ( c : QuadDobl_Speelpenning_Convolutions.Link_to_Circuit )
             return Standard_Speelpenning_Convolutions.Link_to_Circuit;
  function to_double
	     ( c : DoblDobl_Speelpenning_Convolutions.Circuits )
	     return Standard_Speelpenning_Convolutions.Circuits;
  function to_double
	     ( c : QuadDobl_Speelpenning_Convolutions.Circuits )
	     return Standard_Speelpenning_Convolutions.Circuits;

  -- DESCRIPTION :
  --   Turns the coefficients in c to double precision
  --   and copies all data of c, and allocates the work space.

  function to_double_double
	     ( c : QuadDobl_Speelpenning_Convolutions.Circuit )
	     return DoblDobl_Speelpenning_Convolutions.Circuit;
  function to_double_double
	     ( c : QuadDobl_Speelpenning_Convolutions.Link_to_Circuit )
             return DoblDobl_Speelpenning_Convolutions.Link_to_Circuit;
  function to_double_double
	     ( c : QuadDobl_Speelpenning_Convolutions.Circuits )
	     return DoblDobl_Speelpenning_Convolutions.Circuits;

  -- DESCRIPTION :
  --   Turns the coefficients in c to double double precision
  --   and copies all data of c, and allocates the work space.

  function to_double
             ( s : QuadDobl_Speelpenning_Convolutions.System )
             return Standard_Speelpenning_Convolutions.System;
  function to_double
             ( s : QuadDobl_Speelpenning_Convolutions.Link_to_System )
             return Standard_Speelpenning_Convolutions.Link_to_System;

  -- DESCRIPTION :
  --   Turns the coefficients in s to double precision,
  --   copies all data of s, and allocates the work space.

  function to_double_double
             ( s : QuadDobl_Speelpenning_Convolutions.System )
             return DoblDobl_Speelpenning_Convolutions.System;
  function to_double_double
             ( s : QuadDobl_Speelpenning_Convolutions.Link_to_System )
             return DoblDobl_Speelpenning_Convolutions.Link_to_System;

  -- DESCRIPTION :
  --   Turns the coefficients in s to double double precision,
  --   copies all data of s, and allocates the work space.

end System_Convolution_Circuits;
