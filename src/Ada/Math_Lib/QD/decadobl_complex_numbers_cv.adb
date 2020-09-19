with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Multprec_Floating_Numbers;          use Multprec_Floating_Numbers;
with Deca_Double_Numbers;                use Deca_Double_Numbers;
with Multprec_DecaDobl_Convertors;       use Multprec_DecaDobl_Convertors;

package body DecaDobl_Complex_Numbers_cv is

  function Standard_to_DecaDobl_Complex
             ( c : Standard_Complex_Numbers.Complex_Number )
             return DecaDobl_Complex_Numbers.Complex_Number is

    res : DecaDobl_Complex_Numbers.Complex_Number;
    cre : constant double_float := Standard_Complex_Numbers.REAL_PART(c);
    cim : constant double_float := Standard_Complex_Numbers.IMAG_PART(c);
    rre : constant deca_double := create(cre);
    rim : constant deca_double := create(cim);

  begin
    res := DecaDobl_Complex_Numbers.Create(rre,rim);
    return res;
  end Standard_to_DecaDobl_Complex;

  function Multprec_to_DecaDobl_Complex
             ( c : Multprec_Complex_Numbers.Complex_Number )
             return DecaDobl_Complex_Numbers.Complex_Number is

    res : DecaDobl_Complex_Numbers.Complex_Number;
    cre : Floating_Number := Multprec_Complex_Numbers.REAL_PART(c);
    cim : Floating_Number := Multprec_Complex_Numbers.IMAG_PART(c);
    rre : constant deca_double := to_deca_double(cre);
    rim : constant deca_double := to_deca_double(cim);

  begin
    Multprec_Floating_Numbers.Clear(cre);
    Multprec_Floating_Numbers.Clear(cim);
    res := DecaDobl_Complex_Numbers.Create(rre,rim);
    return res;
  end Multprec_to_DecaDobl_Complex;

  function DecaDobl_Complex_to_Standard
             ( c : DecaDobl_Complex_Numbers.Complex_Number )
             return Standard_Complex_Numbers.Complex_Number is

    res : Standard_Complex_Numbers.Complex_Number;
    cre : constant deca_double := DecaDobl_Complex_Numbers.REAL_PART(c);
    cim : constant deca_double := DecaDobl_Complex_Numbers.IMAG_PART(c);
    rre : constant double_float := to_double(cre);
    rim : constant double_float := to_double(cim);

  begin
    res := Standard_Complex_Numbers.Create(rre,rim);
    return res;
  end DecaDobl_Complex_to_Standard;

  function DecaDobl_Complex_to_Multprec
             ( c : DecaDobl_Complex_Numbers.Complex_Number )
             return Multprec_Complex_Numbers.Complex_Number is

    res : Multprec_Complex_Numbers.Complex_Number;
    cre : constant deca_double := DecaDobl_Complex_Numbers.REAL_PART(c);
    cim : constant deca_double := DecaDobl_Complex_Numbers.IMAG_PART(c);
    rre : Floating_Number := to_floating_number(cre);
    rim : Floating_Number := to_floating_number(cim);

  begin
    res := Multprec_Complex_Numbers.Create(rre,rim);
    Multprec_Floating_Numbers.Clear(rre);
    Multprec_Floating_Numbers.Clear(rim);
    return res;
  end DecaDobl_Complex_to_Multprec;

end DecaDobl_Complex_Numbers_cv;
