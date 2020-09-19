with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Multprec_Floating_Numbers;          use Multprec_Floating_Numbers;
with Octo_Double_Numbers;                use Octo_Double_Numbers;
with Multprec_OctoDobl_Convertors;       use Multprec_OctoDobl_Convertors;

package body OctoDobl_Complex_Numbers_cv is

  function Standard_to_OctoDobl_Complex
             ( c : Standard_Complex_Numbers.Complex_Number )
             return OctoDobl_Complex_Numbers.Complex_Number is

    res : OctoDobl_Complex_Numbers.Complex_Number;
    cre : constant double_float := Standard_Complex_Numbers.REAL_PART(c);
    cim : constant double_float := Standard_Complex_Numbers.IMAG_PART(c);
    rre : constant octo_double := create(cre);
    rim : constant octo_double := create(cim);

  begin
    res := OctoDobl_Complex_Numbers.Create(rre,rim);
    return res;
  end Standard_to_OctoDobl_Complex;

  function Multprec_to_OctoDobl_Complex
             ( c : Multprec_Complex_Numbers.Complex_Number )
             return OctoDobl_Complex_Numbers.Complex_Number is

    res : OctoDobl_Complex_Numbers.Complex_Number;
    cre : Floating_Number := Multprec_Complex_Numbers.REAL_PART(c);
    cim : Floating_Number := Multprec_Complex_Numbers.IMAG_PART(c);
    rre : constant octo_double := to_octo_double(cre);
    rim : constant octo_double := to_octo_double(cim);

  begin
    Multprec_Floating_Numbers.Clear(cre);
    Multprec_Floating_Numbers.Clear(cim);
    res := OctoDobl_Complex_Numbers.Create(rre,rim);
    return res;
  end Multprec_to_OctoDobl_Complex;

  function OctoDobl_Complex_to_Standard
             ( c : OctoDobl_Complex_Numbers.Complex_Number )
             return Standard_Complex_Numbers.Complex_Number is

    res : Standard_Complex_Numbers.Complex_Number;
    cre : constant octo_double := OctoDobl_Complex_Numbers.REAL_PART(c);
    cim : constant octo_double := OctoDobl_Complex_Numbers.IMAG_PART(c);
    rre : constant double_float := to_double(cre);
    rim : constant double_float := to_double(cim);

  begin
    res := Standard_Complex_Numbers.Create(rre,rim);
    return res;
  end OctoDobl_Complex_to_Standard;

  function OctoDobl_Complex_to_Multprec
             ( c : OctoDobl_Complex_Numbers.Complex_Number )
             return Multprec_Complex_Numbers.Complex_Number is

    res : Multprec_Complex_Numbers.Complex_Number;
    cre : constant octo_double := OctoDobl_Complex_Numbers.REAL_PART(c);
    cim : constant octo_double := OctoDobl_Complex_Numbers.IMAG_PART(c);
    rre : Floating_Number := to_floating_number(cre);
    rim : Floating_Number := to_floating_number(cim);

  begin
    res := Multprec_Complex_Numbers.Create(rre,rim);
    Multprec_Floating_Numbers.Clear(rre);
    Multprec_Floating_Numbers.Clear(rim);
    return res;
  end OctoDobl_Complex_to_Multprec;

end OctoDobl_Complex_Numbers_cv;
