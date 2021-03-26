with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Integer_Vectors;
with Standard_Complex_VecVecs;
with Standard_Complex_Laur_Systems;     use Standard_Complex_Laur_Systems;
with Standard_Complex_Laur_JacoMats;    use Standard_Complex_Laur_JacoMats;
with Double_Lseries_Polynomials;        use Double_Lseries_Polynomials;

package Test_Double_Lseries_Polynomials is

-- DESCRIPTION :
--   Tests the evaluation and differentiation of a polynomial
--   at a sequence of Laurent series, in double precision.

  procedure Test ( dim,nbr,deg,pwr,low,upp : in integer32 );

  -- DESCRIPTION :
  --   Generates random data and runs tests.

  -- ON ENTRY :
  --   dim      the dimension is the number of variables;
  --   nbr      number of monomials in the polynomial;
  --   deg      degree of the series;
  --   pwr      largest power for every variable;
  --   low      lower bound on leading exponents of the series;
  --   upp      upper bound on leading exponents of the series.

  procedure Set_Parameters
              ( p : in Link_to_Laur_Sys;
                neq,dim,tdx,deg : out integer32 );

  -- DESCRIPTION :
  --   Given a Laurent system, sets the characteristic parameters.

  -- ON RETURN :
  --   neq      number of equations;
  --   dim      number of variables;
  --   tdx      index of the symbol t for the series;
  --   deg      truncation degree of the series.

  procedure Test_Input;

  -- DESCRIPTION :
  --   Prompts for a Laurent polynomial system and then constructs the data
  --   to evaluate the system at a vector of Laurent series.

  procedure Chop ( tdx : in integer32;
                   xlead : in Standard_Integer_Vectors.Vector;
                   xcffs : in Standard_Complex_VecVecs.Link_to_VecVec;
                   ylead : out Standard_Integer_Vectors.Vector;
                   ycffs : in Standard_Complex_VecVecs.Link_to_VecVec );

  -- DESCRIPTION :
  --   Copies all data from (xlead, xcffs) to (ylead, ycffs),
  --   except the data at the nonzero index tdx.
 
  -- REQUIRED :
  --   The range of ylead and ycffs should be one less than xlead and ycffs.

  procedure Test_Jacobian_Evaluation
              ( dim,nvr,tdx,deg : integer32; jp : in Jaco_Mat;
                tva : in Table_Vector_Array );

  -- DESCRIPTION :
  --   Tests the evaluation of the Jacobian matrix
  --   at a vector of a random series.

  -- ON ENTRY :
  --   dim      original number of variables, including t;
  --   nvr      number of variables, not including t;
  --   deg      degree of the power series;
  --   jp       symbolic form of the Jacobian matrix;
  --   tva      table representation of the Jacobian.    

  procedure Test_Jacobian;

  -- DESCRIPTION :
  --    Prompts for a Laurent polynomial system
  --    and evaluates its Jacobian matrix.

  procedure Main;

  -- DESCRIPTION :
  --   Prompts for the parameters of the tests and then runs tests.

end Test_Double_Lseries_Polynomials;
