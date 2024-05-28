with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Floating_Vectors;

package Double_Taylor_Developments is

-- DESCRIPTION :
--   Computes Taylor series developments of monomials
--   with positive real exponents.

  function Double_Taylor_Coefficients 
             ( deg : integer32; alpha, point : double_float )
             return Standard_Floating_Vectors.Vector;

  -- DESCRIPTION :
  --   Returns the vector of coefficients of the Taylor series
  --   in the shifted basis.

  -- ON ENTRY :
  --   deg     truncation degree of the Taylor series;
  --   alpha   positive real power of the monomial;
  --   point   base point of the development.

  function Double_Taylor_Value
             ( cff : Standard_Floating_Vectors.Vector;
               arg, point : double_float ) return double_float;

  -- DESCRIPTION :
  --   Returns the value of the Taylor expansion about the point,
  --   given the coefficients of the series in the shifted basis.

  -- ON ENTRY :
  --   deg     truncation degree of the Taylor series;
  --   arg     argument where to evaluate the series;
  --   point   base point of the development.

  function Double_Taylor_Value
             ( cff : Standard_Floating_Vectors.Vector;
               arg : double_float ) return double_float;

  -- DESCRIPTION :
  --   Returns the value of the Taylor expansion about the point,
  --   given the coefficients of the series in the monomial basis.

  -- ON ENTRY :
  --   deg     truncation degree of the Taylor series;
  --   arg     argument where to evaluate the series.

  function Double_Taylor_Expansion
             ( cff : Standard_Floating_Vectors.Vector;
               point : double_float )
             return Standard_Floating_Vectors.Vector;

  -- DESCRIPTION :
  --   Returns the Taylor expansion about the point,
  --   expressed in the standard monomial basis.

  -- ON ENTRY :
  --   deg     truncation degree of the Taylor series;
  --   point   base point of the development.

end Double_Taylor_Developments;
