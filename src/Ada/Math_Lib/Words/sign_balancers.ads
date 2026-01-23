with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
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

-- EQUILIZE signs via one last bit deduction/augmentation

  function One_Last_Bit ( nbr : double_float ) return double_float;

  -- DESCRIPTION :
  --   Given a double float number nbr, computes the one last bit of
  --   its fraction, using the exponent of nbr in the returned double.
  --   Subtracting the returned number from nbr results in a double
  --   where the fraction is one less than the fraction of nbr.

  procedure Equalize_Signs ( hi,lo : in out double_float;
                             vrblvl : in integer32 := 0 );

  -- DESCRIPTION :
  --   Given hi*lo < 0.0, equalizes the signs of hi and lo
  --   by subtracting the one last bit of hi and adding it to lo.
  --   If vrblvl > 0, prints results of intermediate computations.

  procedure Equalize_Signs ( hihi,lohi,hilo,lolo : in out double_float;
                             vrblvl : in integer32 := 0 );

  -- DESCRIPTION :
  --   Equalizes the signs of the parts of a quad double.

  procedure Equalize_Signs ( x : in out double_double;
                             vrblvl : in integer32 := 0 );

  -- DESCRIPTION :
  --   If the high and the low part of x have different signs,
  --   then the one last bit of the high part x is reduced and
  --   augmented to the low part of x so the double double
  --   representation of x is sign balanced.
  --   If vrblvl > 0, prints results of intermediate computations.

  procedure Equalize_Signs ( x : in out quad_double;
                             vrblvl : in integer32 := 0 );

  -- DESCRIPTION :
  --   Equalizes the signs of the quad double x.

-- BALANCE via redistribution of bits

  procedure Sign_Balance ( hi,lo : in out double_float;
                           vrblvl : in integer32 := 0 );

  -- DESCRIPTION :
  --   Given hi*lo < 0.0, balances the sign by redistributing
  --   the bits from hi to lo.
  --   If vrblvl > 0, prints results of intermediate computations.

  procedure Sign_Balance ( x : in out double_double;
                           vrblvl : in integer32 := 0 );

  -- DESCRIPTION :
  --   If the high and the low part of x have different signs,
  --   then the bits of x are redistributed so the double double
  --   representation of x is sign balanced.
  --   If vrblvl > 0, prints results of intermediate computations.

  procedure Sign_Balance ( hihi,lohi,hilo,lolo : in out double_float;
                           vrblvl : in integer32 := 0 );

  -- DESCRIPTION :
  --   Balances the doubles of a quad double so all parts have the same sign.

  procedure Sign_Balance
              ( hihihi,lohihi,hilohi,lolohi : in out double_float;
                hihilo,lohilo,hilolo,lololo : in out double_float;
                vrblvl : in integer32 := 0 );

  -- DESCRIPTION :
  --   Balances the doubles of an octo double so all parts have the same sign.

  procedure Sign_Balance
              ( hihihihi,lohihihi,hilohihi,lolohihi : in out double_float;
                hihilohi,lohilohi,hilolohi,lololohi : in out double_float;
                hihihilo,lohihilo,hilohilo,lolohilo : in out double_float;
                hihilolo,lohilolo,hilololo,lolololo : in out double_float;
                vrblvl : in integer32 := 0 );

  -- DESCRIPTION :
  --   Balances the doubles of a hexa double so all parts have the same sign.

  procedure Sign_Balance ( x : in out quad_double;
                           vrblvl : in integer32 := 0 );
  procedure Sign_Balance ( x : in out octo_double;
                           vrblvl : in integer32 := 0 );
  procedure Sign_Balance ( x : in out hexa_double;
                           vrblvl : in integer32 := 0 );

  -- DESCRIPTION :
  --   Redistributes the bits so all parts of have the same sign.

end Sign_Balancers;
