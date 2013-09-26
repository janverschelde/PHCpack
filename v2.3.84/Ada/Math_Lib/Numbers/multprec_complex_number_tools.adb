with Multprec_Floating_Numbers;          use Multprec_Floating_Numbers;

package body Multprec_Complex_Number_Tools is

  function Round ( c : Multprec_Complex_Numbers.Complex_Number )
                 return Standard_Complex_Numbers.Complex_Number is

    res : Standard_Complex_Numbers.Complex_Number;
    cre : Floating_Number := Multprec_Complex_Numbers.REAL_PART(c);
    cim : Floating_Number := Multprec_Complex_Numbers.IMAG_PART(c);

  begin
    res := Standard_Complex_Numbers.Create(Round(cre),Round(cim));
    Clear(cre); Clear(cim);
    return res;
  end Round;

  function Create ( f : double_float )
                  return Multprec_Complex_Numbers.Complex_Number is

    ff : constant Floating_Number := Create(f);

  begin
    return Multprec_Complex_Numbers.Create(ff);
  end Create;

  function Create ( re,im : double_float )
                  return Multprec_Complex_Numbers.Complex_Number is

    f_re : constant Floating_Number := Create(re);
    f_im : constant Floating_Number := Create(im);

  begin
    return Multprec_Complex_Numbers.Create(f_re,f_im);
  end Create;

  function Create ( c : Standard_Complex_Numbers.Complex_Number )
                  return Multprec_Complex_Numbers.Complex_Number is

    cre : constant double_float := Standard_Complex_Numbers.REAL_PART(c);
    cim : constant double_float := Standard_Complex_Numbers.IMAG_PART(c);
    rre : Floating_Number := Multprec_Floating_Numbers.Create(cre);
    rim : Floating_Number := Multprec_Floating_Numbers.Create(cim);
    res : constant Multprec_Complex_Numbers.Complex_Number
        := Multprec_Complex_Numbers.Create(rre,rim);

  begin
    Clear(rre); Clear(rim);
    return res;
  end Create;

  procedure Set_Size ( c : in out Multprec_Complex_Numbers.Complex_Number;
                       size : in natural32 ) is

    cre : Floating_Number := Multprec_Complex_Numbers.REAL_PART(c);
    cim : Floating_Number := Multprec_Complex_Numbers.IMAG_PART(c);

  begin
    Set_Size(cre,size);
    Set_Size(cim,size);
    Multprec_Complex_Numbers.Clear(c);
    c := Multprec_Complex_Numbers.Create(cre,cim);
    Clear(cre); Clear(cim);
  end Set_Size;

end Multprec_Complex_Number_Tools;
