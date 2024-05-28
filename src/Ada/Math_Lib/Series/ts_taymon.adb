with text_io;                            use text_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;
with Standard_Mathematical_Functions;    use Standard_Mathematical_Functions;
with Standard_Floating_Vectors;
with Standard_Floating_Vectors_io;       use Standard_Floating_Vectors_io;
with Double_Taylor_Developments;         use Double_Taylor_Developments;

procedure ts_taymon is

-- DESCRIPTION :
--   Tests on Taylor series developments of monomials
--   with positive real exponents.

  procedure Double_Test
              ( deg : in integer32; alpha, point : in double_float ) is

  -- DESCRIPTION :
  --   Runs a test on the double Taylor coefficients.
  --
  -- ON ENTRY :
  --   deg      truncation degree of the Taylor series;
  --   alpha    positive real power of the monomial;
  --   point    base point of the development.

    cff : constant Standard_Floating_Vectors.Vector(0..deg)
        := Double_Taylor_Coefficients(deg,alpha,point);
    tcf : constant Standard_Floating_Vectors.Vector(0..deg)
        := Double_Taylor_Expansion(cff,point);
    apt : double_float := point - 1.0e-2; -- Standard_Random_Numbers.Random;
    mval,tval1,tval2,err : double_float;

  begin
    new_line;
    put_line("The coefficients in the shifted basis :");
    put_line(cff);
    put_line("The coefficients in the monomial basis :");
    put_line(tcf);
    if apt < 0.0
     then apt := -apt;
    end if;
    put("The argument to evaluate the series :"); put(apt); new_line;
    tval1 := Double_Taylor_Value(cff,apt,point);
    put("The value in the shifted basis  :"); put(tval1); new_line;
    tval2 := Double_Taylor_Value(tcf,apt);
    put("The value in the monomial basis :"); put(tval2); new_line;
    mval := apt**alpha;
    put("The value of the monomial       :"); put(mval); new_line;
    err := abs(mval - tval1);
    put("-> the error on the values :"); put(err,3); new_line;
  end Double_Test;

  procedure Main is

  -- DESCRIPTION :
  --   Prompts for the input parameters.

    deg : integer32 := 0;
    alpha,point : double_float := 0.0;

  begin
    new_line;
    put("Give the truncation degree : "); get(deg);
    put("Give the positive real power : "); get(alpha);
    put("Give the positive real point : "); get(point);
    new_line;
    put("-> the truncation degree : "); put(deg,1); new_line;
    put("-> power of the monomial :"); put(alpha); new_line;
    put("-> point of development  :"); put(point); new_line;
    Double_Test(deg,alpha,point);
  end Main;

begin
  Main;
end ts_taymon;
