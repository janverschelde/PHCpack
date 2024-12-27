with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Double_Double_Numbers;              use Double_Double_Numbers;
with Quad_Double_Numbers;                use Quad_Double_Numbers;
with Octo_Double_Numbers;                use Octo_Double_Numbers;
with Hexa_Double_Numbers;                use Hexa_Double_Numbers;

package Sign_Balancers is

-- DESCRIPTION :
--   A multiple double number is an unevaluated sum of doubles.
--   The sign balancers in this package ensure that the lower doubles
--   have the same sign as the leading double, as this is convenient
--   for the accurate computation of the product.

  function Different_Sign ( x,y : double_float ) return boolean;

  -- DESCRIPTION :
  --   Returns true if x and y have different signs, false otherwise.
  --   This function should be used to avoid underflow of x*y
  --   which causes the test x*y < 0.0 to fail.

  function Is_Sign_Balanced ( x : double_double ) return boolean;
  function Is_Sign_Balanced ( x : quad_double ) return boolean;
  function Is_Sign_Balanced ( x : octo_double ) return boolean;
  function Is_Sign_Balanced ( x : hexa_double ) return boolean;

  -- DESCRIPTION :
  --   Returns true of all parts of x have the same sign.

  procedure Sign_Balance ( hi,lo : in out double_float;
                           verbose : in boolean := true );

  -- DESCRIPTION :
  --   Given hi*lo < 0.0, balances the sign by redistributing
  --   the bits from hi to lo.
  --   If verbose, prints results of intermediate computations.

  procedure Sign_Balance ( x : in out double_double;
                           verbose : in boolean := true );

  -- DESCRIPTION :
  --   If the high and the low part of x have different signs,
  --   then the bits of x are redistributed so the double double
  --   representation of x is sign balanced.
  --   If verbose, prints results of intermediate computations.

  procedure Sign_Balance ( hihi,lohi,hilo,lolo : in out double_float;
                           verbose : in boolean := true );

  -- DESCRIPTION :
  --   Balances the doubles of a quad double so all parts have the same sign.

  procedure Sign_Balance
              ( hihihi,lohihi,hilohi,lolohi : in out double_float;
                hihilo,lohilo,hilolo,lololo : in out double_float;
                verbose : in boolean := true );

  -- DESCRIPTION :
  --   Balances the doubles of an octo double so all parts have the same sign.

  procedure Sign_Balance
              ( hihihihi,lohihihi,hilohihi,lolohihi : in out double_float;
                hihilohi,lohilohi,hilolohi,lololohi : in out double_float;
                hihihilo,lohihilo,hilohilo,lolohilo : in out double_float;
                hihilolo,lohilolo,hilololo,lolololo : in out double_float;
                verbose : in boolean := true );

  -- DESCRIPTION :
  --   Balances the doubles of a hexa double so all parts have the same sign.

  procedure Sign_Balance ( x : in out quad_double;
                           verbose : in boolean := true );
  procedure Sign_Balance ( x : in out octo_double;
                           verbose : in boolean := true );
  procedure Sign_Balance ( x : in out hexa_double;
                           verbose : in boolean := true );

  -- DESCRIPTION :
  --   Redistributes the bits so all parts of have the same sign.

end Sign_Balancers;
