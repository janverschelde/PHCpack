with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;
with Standard_Complex_Numbers_io;        use Standard_Complex_Numbers_io;
with Standard_Random_Numbers;            use Standard_Random_Numbers;
with Standard_Mathematical_Functions;    use Standard_Mathematical_Functions;
with Standard_Complex_Vectors_io;        use Standard_Complex_Vectors_io;
with Standard_Complex_VecVecs;
with Standard_Complex_Polynomials;
with Standard_Complex_Polynomials_io;    use Standard_Complex_Polynomials_io;
with Standard_Polynomial_Interpolators;  use Standard_Polynomial_Interpolators;

procedure ts_adpolin is

-- DESCRIPTION :
--   Adaptive precision multivariate polynomial interpolation.

  function Interpolation_Condition_Number
              ( p : Standard_Complex_Polynomials.Poly;
                v : Standard_Complex_VecVecs.VecVec ) return double_float is

  -- DESCRIPTION :
  --   Returns an estimate for the condition number of the interpolation
  --   problem, using the samples in v.

    res : double_float;
    ip : Standard_Complex_Polynomials.Poly;

  begin
    Interpolate(p,v,ip,res);
    put_line("The interpolating polynomial : "); put_line(ip);
    Standard_Complex_Polynomials.Clear(ip);
    return res;
  end Interpolation_Condition_Number;

  function Interpolation_Condition_Number
              ( p : Standard_Complex_Polynomials.Poly;
                m : natural32; scalfac : double_float ) return double_float is

  -- DESCRIPTION :
  --   Returns an estimate for the condition number of the interpolation
  --   problem, using m samples of the polynomial.
  --   

    vv : Standard_Complex_VecVecs.VecVec := Sample(p,m,scalfac);
    rcond : constant double_float := Interpolation_Condition_Number(p,vv);

  begin
    Standard_Complex_VecVecs.Clear(vv);
    return rcond;
  end Interpolation_Condition_Number;

  procedure Determine_Working_Precision
              ( sp : in Standard_Complex_Polynomials.Poly;
                desacc : in natural32; scalfac : in double_float ) is

    nbterms,m,loss,wrkacc,size : natural32 := 0;
    rcond : double_float;

  begin
    new_line;
    nbterms := Standard_Complex_Polynomials.Number_of_Terms(sp);
    put("The number of terms : "); put(nbterms,1); new_line;
    m := nbterms-1;
    rcond := Interpolation_Condition_Number(sp,m,scalfac);
    put("Estimate for inverse of condition# : ");
    put(rcond,3); new_line;
    loss := natural32(LOG10(rcond));
    put("Loss of decimal places : "); put(loss,1); new_line;
    new_line;
    wrkacc := desacc - loss;
    size := wrkacc/8;
    if wrkacc mod 8 /= 0
     then size := size + 1;
    end if;
    put("With the desired accuracy (#decimal places) : ");
    put(desacc,1); new_line;
    put("Working precision (#decimal places) will be : ");
    put((size+1)*8,1); new_line;
  end Determine_Working_Precision;

  procedure Main is

    d,n,cff,desacc : natural32 := 0;
    sp : Standard_Complex_Polynomials.Poly;
    scalfac : double_float;
    ans : character;

  begin
    new_line;
    put_line("Multivariate Polynomial Interpolation with adaptive precision");
    loop
      new_line;                                            -- input section
      put("Give the degree : "); get(d);
      put("Give the number of variables : "); get(n);
      put_line("Choose one of the following : ");
      put_line("  1. all coefficients equal to one;");
      put_line("  2. random real coefficients;");
      put_line("  3. random complex coefficients.");
      put("Type 1,2, or 3 to make your choice : "); get(cff);
      put("Give real number, magnitude for samples : "); get(scalfac);
      put("Give desired accuracy (#decimal places) : "); get(desacc);  
      sp := Create(d,n,cff);                               -- end of input
      Determine_Working_Precision(sp,desacc,scalfac);
      new_line;
      put("Do you want more tests ? (y/n) "); Ask_Yes_or_No(ans);
      exit when (ans /= 'y');
    end loop;
  end Main;

begin
  Main;
end ts_adpolin;
