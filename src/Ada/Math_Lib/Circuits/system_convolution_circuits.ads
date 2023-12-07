with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
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
with Standard_Speelpenning_Convolutions;
with DoblDobl_Speelpenning_Convolutions;
with TripDobl_Speelpenning_Convolutions;
with QuadDobl_Speelpenning_Convolutions;
with PentDobl_Speelpenning_Convolutions;
with OctoDobl_Speelpenning_Convolutions;
with DecaDobl_Speelpenning_Convolutions;
with HexaDobl_Speelpenning_Convolutions;

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
  function Make_Convolution_Circuit
             ( p : PentDobl_Complex_Polynomials.Poly;
               d : natural32 )
             return PentDobl_Speelpenning_Convolutions.Circuit;
  function Make_Convolution_Circuit
             ( p : OctoDobl_Complex_Polynomials.Poly;
               d : natural32 )
             return OctoDobl_Speelpenning_Convolutions.Circuit;
  function Make_Convolution_Circuit
             ( p : DecaDobl_Complex_Polynomials.Poly;
               d : natural32 )
             return DecaDobl_Speelpenning_Convolutions.Circuit;
  function Make_Convolution_Circuit
             ( p : HexaDobl_Complex_Polynomials.Poly;
               d : natural32 )
             return HexaDobl_Speelpenning_Convolutions.Circuit;

  -- DESCRIPTION :
  --   Returns the convolution circuit representation of the polynomial p,
  --   in double, double double, triple double, quad double, penta double,
  --   octo double, deca double, or hexa double precision.
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
  function Make_Convolution_Circuit
             ( p : PentDobl_CSeries_Polynomials.Poly )
             return PentDobl_Speelpenning_Convolutions.Circuit;
  function Make_Convolution_Circuit
             ( p : OctoDobl_CSeries_Polynomials.Poly )
             return OctoDobl_Speelpenning_Convolutions.Circuit;
  function Make_Convolution_Circuit
             ( p : DecaDobl_CSeries_Polynomials.Poly )
             return DecaDobl_Speelpenning_Convolutions.Circuit;
  function Make_Convolution_Circuit
             ( p : HexaDobl_CSeries_Polynomials.Poly )
             return HexaDobl_Speelpenning_Convolutions.Circuit;

  -- DESCRIPTION :
  --   Returns the convolution circuit representation of the polynomial p,
  --   in double, double double, triple double, quad double, penta double,
  --   octo double, deca double, or hexa double precision.

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
  function Make_Convolution_Circuits
             ( p : PentDobl_Complex_Poly_Systems.Poly_Sys;
               d : natural32 )
             return PentDobl_Speelpenning_Convolutions.Circuits;
  function Make_Convolution_Circuits
             ( p : OctoDobl_Complex_Poly_Systems.Poly_Sys;
               d : natural32 )
             return OctoDobl_Speelpenning_Convolutions.Circuits;
  function Make_Convolution_Circuits
             ( p : DecaDobl_Complex_Poly_Systems.Poly_Sys;
               d : natural32 )
             return DecaDobl_Speelpenning_Convolutions.Circuits;
  function Make_Convolution_Circuits
             ( p : HexaDobl_Complex_Poly_Systems.Poly_Sys;
               d : natural32 )
             return HexaDobl_Speelpenning_Convolutions.Circuits;

  -- DESCRIPTION :
  --   Returns the convolution circuits representation of the polynomials
  --   in p, in double, double double, triple double, quad double,
  --   penta double, octo double, deca double, or hexa double precision.
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
  function Make_Convolution_Circuits
             ( p : PentDobl_CSeries_Poly_Systems.Poly_Sys )
             return PentDobl_Speelpenning_Convolutions.Circuits;
  function Make_Convolution_Circuits
             ( p : OctoDobl_CSeries_Poly_Systems.Poly_Sys )
             return OctoDobl_Speelpenning_Convolutions.Circuits;
  function Make_Convolution_Circuits
             ( p : DecaDobl_CSeries_Poly_Systems.Poly_Sys )
             return DecaDobl_Speelpenning_Convolutions.Circuits;
  function Make_Convolution_Circuits
             ( p : HexaDobl_CSeries_Poly_Systems.Poly_Sys )
             return HexaDobl_Speelpenning_Convolutions.Circuits;

  -- DESCRIPTION :
  --   Returns the convolution circuits representation of the polynomials
  --   in p, in double, double double, triple double, quad double,
  --   penta double, octo double, deca double, or hexa double precision.

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
  function Make_Convolution_System
             ( p : PentDobl_CSeries_Poly_Systems.Poly_Sys;
               d : natural32 )
             return PentDobl_Speelpenning_Convolutions.System;
  function Make_Convolution_System
             ( p : PentDobl_CSeries_Poly_Systems.Poly_Sys;
               d : natural32 )
             return PentDobl_Speelpenning_Convolutions.Link_to_System;
  function Make_Convolution_System
             ( p : OctoDobl_CSeries_Poly_Systems.Poly_Sys;
               d : natural32 )
             return OctoDobl_Speelpenning_Convolutions.System;
  function Make_Convolution_System
             ( p : OctoDobl_CSeries_Poly_Systems.Poly_Sys;
               d : natural32 )
             return OctoDobl_Speelpenning_Convolutions.Link_to_System;
  function Make_Convolution_System
             ( p : DecaDobl_CSeries_Poly_Systems.Poly_Sys;
               d : natural32 )
             return DecaDobl_Speelpenning_Convolutions.System;
  function Make_Convolution_System
             ( p : DecaDobl_CSeries_Poly_Systems.Poly_Sys;
               d : natural32 )
             return DecaDobl_Speelpenning_Convolutions.Link_to_System;
  function Make_Convolution_System
             ( p : HexaDobl_CSeries_Poly_Systems.Poly_Sys;
               d : natural32 )
             return HexaDobl_Speelpenning_Convolutions.System;
  function Make_Convolution_System
             ( p : HexaDobl_CSeries_Poly_Systems.Poly_Sys;
               d : natural32 )
             return HexaDobl_Speelpenning_Convolutions.Link_to_System;

  -- DESCRIPTION :
  --   Wraps the making of the convolution circuits of the polynomials in p,
  --   allocating all work space, for degree d of the power series,
  --   in double, double double, triple_double, quad double,
  --   penta double, octo double, deca double, or hexa double precision.

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

  function to_triple_double
	     ( c : QuadDobl_Speelpenning_Convolutions.Circuit )
	     return TripDobl_Speelpenning_Convolutions.Circuit;
  function to_triple_double
	     ( c : QuadDobl_Speelpenning_Convolutions.Link_to_Circuit )
             return TripDobl_Speelpenning_Convolutions.Link_to_Circuit;
  function to_triple_double
	     ( c : QuadDobl_Speelpenning_Convolutions.Circuits )
	     return TripDobl_Speelpenning_Convolutions.Circuits;

  -- DESCRIPTION :
  --   Turns the coefficients in c to triple double precision
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

  function to_triple_double
             ( s : QuadDobl_Speelpenning_Convolutions.System )
             return TripDobl_Speelpenning_Convolutions.System;
  function to_triple_double
             ( s : QuadDobl_Speelpenning_Convolutions.Link_to_System )
             return TripDobl_Speelpenning_Convolutions.Link_to_System;

  -- DESCRIPTION :
  --   Turns the coefficients in s to triple double precision,
  --   copies all data of s, and allocates the work space.

  function to_double
	     ( c : DecaDobl_Speelpenning_Convolutions.Circuit )
	     return Standard_Speelpenning_Convolutions.Circuit;
  function to_double_double
	     ( c : DecaDobl_Speelpenning_Convolutions.Circuit )
	     return DoblDobl_Speelpenning_Convolutions.Circuit;
  function to_triple_double
	     ( c : DecaDobl_Speelpenning_Convolutions.Circuit )
	     return TripDobl_Speelpenning_Convolutions.Circuit;
  function to_quad_double
	     ( c : DecaDobl_Speelpenning_Convolutions.Circuit )
	     return QuadDobl_Speelpenning_Convolutions.Circuit;
  function to_penta_double
	     ( c : DecaDobl_Speelpenning_Convolutions.Circuit )
	     return PentDobl_Speelpenning_Convolutions.Circuit;
  function to_octo_double
	     ( c : DecaDobl_Speelpenning_Convolutions.Circuit )
	     return OctoDobl_Speelpenning_Convolutions.Circuit;

  -- DESCRIPTION :
  --   The circuit on return has the same coefficients as c,
  --   but of a lower precision.

  function to_double
	     ( c : DecaDobl_Speelpenning_Convolutions.Link_to_Circuit )
             return Standard_Speelpenning_Convolutions.Link_to_Circuit;
  function to_double_double
	     ( c : DecaDobl_Speelpenning_Convolutions.Link_to_Circuit )
             return DoblDobl_Speelpenning_Convolutions.Link_to_Circuit;
  function to_triple_double
	     ( c : DecaDobl_Speelpenning_Convolutions.Link_to_Circuit )
             return TripDobl_Speelpenning_Convolutions.Link_to_Circuit;
  function to_quad_double
	     ( c : DecaDobl_Speelpenning_Convolutions.Link_to_Circuit )
             return QuadDobl_Speelpenning_Convolutions.Link_to_Circuit;
  function to_penta_double
	     ( c : DecaDobl_Speelpenning_Convolutions.Link_to_Circuit )
             return PentDobl_Speelpenning_Convolutions.Link_to_Circuit;
  function to_octo_double
	     ( c : DecaDobl_Speelpenning_Convolutions.Link_to_Circuit )
             return OctoDobl_Speelpenning_Convolutions.Link_to_Circuit;

  -- DESCRIPTION :
  --   The circuit on return has the same coefficients as c,
  --   but of a lower precision.

  function to_double
	     ( c : DecaDobl_Speelpenning_Convolutions.Circuits )
             return Standard_Speelpenning_Convolutions.Circuits;
  function to_double_double
	     ( c : DecaDobl_Speelpenning_Convolutions.Circuits )
             return DoblDobl_Speelpenning_Convolutions.Circuits;
  function to_triple_double
	     ( c : DecaDobl_Speelpenning_Convolutions.Circuits )
             return TripDobl_Speelpenning_Convolutions.Circuits;
  function to_quad_double
	     ( c : DecaDobl_Speelpenning_Convolutions.Circuits )
             return QuadDobl_Speelpenning_Convolutions.Circuits;
  function to_penta_double
	     ( c : DecaDobl_Speelpenning_Convolutions.Circuits )
             return PentDobl_Speelpenning_Convolutions.Circuits;
  function to_octo_double
	     ( c : DecaDobl_Speelpenning_Convolutions.Circuits )
             return OctoDobl_Speelpenning_Convolutions.Circuits;

  -- DESCRIPTION :
  --   The circuits on return have the same coefficients as c,
  --   but of a lower precision.

  function to_double
             ( s : DecaDobl_Speelpenning_Convolutions.System )
             return Standard_Speelpenning_Convolutions.System;
  function to_double_double
             ( s : DecaDobl_Speelpenning_Convolutions.System )
             return DoblDobl_Speelpenning_Convolutions.System;
  function to_triple_double
             ( s : DecaDobl_Speelpenning_Convolutions.System )
             return TripDobl_Speelpenning_Convolutions.System;
  function to_quad_double
             ( s : DecaDobl_Speelpenning_Convolutions.System )
             return QuadDobl_Speelpenning_Convolutions.System;
  function to_penta_double
             ( s : DecaDobl_Speelpenning_Convolutions.System )
             return PentDobl_Speelpenning_Convolutions.System ;
  function to_octo_double
             ( s : DecaDobl_Speelpenning_Convolutions.System )
             return OctoDobl_Speelpenning_Convolutions.System;

  -- DESCRIPTION :
  --   The system on return has the same coefficients as s,
  --   but of a lower precision.

  function to_double
             ( s : DecaDobl_Speelpenning_Convolutions.Link_to_System )
             return Standard_Speelpenning_Convolutions.Link_to_System;
  function to_double_double
             ( s : DecaDobl_Speelpenning_Convolutions.Link_to_System )
             return DoblDobl_Speelpenning_Convolutions.Link_to_System;
  function to_triple_double
             ( s : DecaDobl_Speelpenning_Convolutions.Link_to_System )
             return TripDobl_Speelpenning_Convolutions.Link_to_System;
  function to_quad_double
             ( s : DecaDobl_Speelpenning_Convolutions.Link_to_System )
             return QuadDobl_Speelpenning_Convolutions.Link_to_System;
  function to_penta_double
             ( s : DecaDobl_Speelpenning_Convolutions.Link_to_System )
             return PentDobl_Speelpenning_Convolutions.Link_to_System ;
  function to_octo_double
             ( s : DecaDobl_Speelpenning_Convolutions.Link_to_System )
             return OctoDobl_Speelpenning_Convolutions.Link_to_System;

  -- DESCRIPTION :
  --   The system on return has the same coefficients as s,
  --   but of a lower precision.

end System_Convolution_Circuits;
