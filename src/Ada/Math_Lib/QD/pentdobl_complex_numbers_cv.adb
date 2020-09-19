with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Multprec_Floating_Numbers;          use Multprec_Floating_Numbers;
with Penta_Double_Numbers;               use Penta_Double_Numbers;
with Multprec_PentDobl_Convertors;       use Multprec_PentDobl_Convertors;

package body PentDobl_Complex_Numbers_cv is

  function Standard_to_PentDobl_Complex
             ( c : Standard_Complex_Numbers.Complex_Number )
             return PentDobl_Complex_Numbers.Complex_Number is

    res : PentDobl_Complex_Numbers.Complex_Number;
    cre : constant double_float := Standard_Complex_Numbers.REAL_PART(c);
    cim : constant double_float := Standard_Complex_Numbers.IMAG_PART(c);
    rre : constant penta_double := create(cre);
    rim : constant penta_double := create(cim);

  begin
    res := PentDobl_Complex_Numbers.Create(rre,rim);
    return res;
  end Standard_to_PentDobl_Complex;

  function Multprec_to_PentDobl_Complex
             ( c : Multprec_Complex_Numbers.Complex_Number )
             return PentDobl_Complex_Numbers.Complex_Number is

    res : PentDobl_Complex_Numbers.Complex_Number;
    cre : Floating_Number := Multprec_Complex_Numbers.REAL_PART(c);
    cim : Floating_Number := Multprec_Complex_Numbers.IMAG_PART(c);
    rre : constant penta_double := to_penta_double(cre);
    rim : constant penta_double := to_penta_double(cim);

  begin
    Multprec_Floating_Numbers.Clear(cre);
    Multprec_Floating_Numbers.Clear(cim);
    res := PentDobl_Complex_Numbers.Create(rre,rim);
    return res;
  end Multprec_to_PentDobl_Complex;

  function PentDobl_Complex_to_Standard
             ( c : PentDobl_Complex_Numbers.Complex_Number )
             return Standard_Complex_Numbers.Complex_Number is

    res : Standard_Complex_Numbers.Complex_Number;
    cre : constant penta_double := PentDobl_Complex_Numbers.REAL_PART(c);
    cim : constant penta_double := PentDobl_Complex_Numbers.IMAG_PART(c);
    rre : constant double_float := to_double(cre);
    rim : constant double_float := to_double(cim);

  begin
    res := Standard_Complex_Numbers.Create(rre,rim);
    return res;
  end PentDobl_Complex_to_Standard;

  function PentDobl_Complex_to_Multprec
             ( c : PentDobl_Complex_Numbers.Complex_Number )
             return Multprec_Complex_Numbers.Complex_Number is

    res : Multprec_Complex_Numbers.Complex_Number;
    cre : constant penta_double := PentDobl_Complex_Numbers.REAL_PART(c);
    cim : constant penta_double := PentDobl_Complex_Numbers.IMAG_PART(c);
    rre : Floating_Number := to_floating_number(cre);
    rim : Floating_Number := to_floating_number(cim);

  begin
    res := Multprec_Complex_Numbers.Create(rre,rim);
    Multprec_Floating_Numbers.Clear(rre);
    Multprec_Floating_Numbers.Clear(rim);
    return res;
  end PentDobl_Complex_to_Multprec;

end PentDobl_Complex_Numbers_cv;
