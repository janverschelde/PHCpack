with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Multprec_Floating_Numbers;          use Multprec_Floating_Numbers;
with Double_Double_Numbers;              use Double_Double_Numbers;
with Multprec_DoblDobl_Convertors;       use Multprec_DoblDobl_Convertors;

package body DoblDobl_Complex_Numbers_cv is

  function Standard_to_DoblDobl_Complex
             ( c : Standard_Complex_Numbers.Complex_Number )
             return DoblDobl_Complex_Numbers.Complex_Number is

    res : DoblDobl_Complex_Numbers.Complex_Number;
    cre : constant double_float := Standard_Complex_Numbers.REAL_PART(c);
    cim : constant double_float := Standard_Complex_Numbers.IMAG_PART(c);
    rre : constant double_double := create(cre);
    rim : constant double_double := create(cim);

  begin
    res := DoblDobl_Complex_Numbers.Create(rre,rim);
    return res;
  end Standard_to_DoblDobl_Complex;

  function Multprec_to_DoblDobl_Complex
             ( c : Multprec_Complex_Numbers.Complex_Number )
             return DoblDobl_Complex_Numbers.Complex_Number is

    res : DoblDobl_Complex_Numbers.Complex_Number;
    cre : Floating_Number := Multprec_Complex_Numbers.REAL_PART(c);
    cim : Floating_Number := Multprec_Complex_Numbers.IMAG_PART(c);
    rre : constant double_double := to_double_double(cre);
    rim : constant double_double := to_double_double(cim);

  begin
    Multprec_Floating_Numbers.Clear(cre);
    Multprec_Floating_Numbers.Clear(cim);
    res := DoblDobl_Complex_Numbers.Create(rre,rim);
    return res;
  end Multprec_to_DoblDobl_Complex;

  function DoblDobl_Complex_to_Standard
             ( c : DoblDobl_Complex_Numbers.Complex_Number )
             return Standard_Complex_Numbers.Complex_Number is

    res : Standard_Complex_Numbers.Complex_Number;
    cre : constant double_double := DoblDobl_Complex_Numbers.REAL_PART(c);
    cim : constant double_double := DoblDobl_Complex_Numbers.IMAG_PART(c);
    rre : constant double_float := to_double(cre);
    rim : constant double_float := to_double(cim);

  begin
    res := Standard_Complex_Numbers.Create(rre,rim);
    return res;
  end DoblDobl_Complex_to_Standard;

  function DoblDobl_Complex_to_Multprec
             ( c : DoblDobl_Complex_Numbers.Complex_Number )
             return Multprec_Complex_Numbers.Complex_Number is

    res : Multprec_Complex_Numbers.Complex_Number;
    cre : constant double_double := DoblDobl_Complex_Numbers.REAL_PART(c);
    cim : constant double_double := DoblDobl_Complex_Numbers.IMAG_PART(c);
    rre : Floating_Number := to_floating_number(cre);
    rim : Floating_Number := to_floating_number(cim);

  begin
    res := Multprec_Complex_Numbers.Create(rre,rim);
    Multprec_Floating_Numbers.Clear(rre);
    Multprec_Floating_Numbers.Clear(rim);
    return res;
  end DoblDobl_Complex_to_Multprec;

end DoblDobl_Complex_Numbers_cv;
