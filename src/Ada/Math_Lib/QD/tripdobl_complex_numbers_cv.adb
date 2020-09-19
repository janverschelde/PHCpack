with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Multprec_Floating_Numbers;          use Multprec_Floating_Numbers;
with Triple_Double_Numbers;              use Triple_Double_Numbers;
with Multprec_TripDobl_Convertors;       use Multprec_TripDobl_Convertors;

package body TripDobl_Complex_Numbers_cv is

  function Standard_to_TripDobl_Complex
             ( c : Standard_Complex_Numbers.Complex_Number )
             return TripDobl_Complex_Numbers.Complex_Number is

    res : TripDobl_Complex_Numbers.Complex_Number;
    cre : constant double_float := Standard_Complex_Numbers.REAL_PART(c);
    cim : constant double_float := Standard_Complex_Numbers.IMAG_PART(c);
    rre : constant triple_double := create(cre);
    rim : constant triple_double := create(cim);

  begin
    res := TripDobl_Complex_Numbers.Create(rre,rim);
    return res;
  end Standard_to_TripDobl_Complex;

  function Multprec_to_TripDobl_Complex
             ( c : Multprec_Complex_Numbers.Complex_Number )
             return TripDobl_Complex_Numbers.Complex_Number is

    res : TripDobl_Complex_Numbers.Complex_Number;
    cre : Floating_Number := Multprec_Complex_Numbers.REAL_PART(c);
    cim : Floating_Number := Multprec_Complex_Numbers.IMAG_PART(c);
    rre : constant triple_double := to_triple_double(cre);
    rim : constant triple_double := to_triple_double(cim);

  begin
    Multprec_Floating_Numbers.Clear(cre);
    Multprec_Floating_Numbers.Clear(cim);
    res := TripDobl_Complex_Numbers.Create(rre,rim);
    return res;
  end Multprec_to_TripDobl_Complex;

  function TripDobl_Complex_to_Standard
             ( c : TripDobl_Complex_Numbers.Complex_Number )
             return Standard_Complex_Numbers.Complex_Number is

    res : Standard_Complex_Numbers.Complex_Number;
    cre : constant triple_double := TripDobl_Complex_Numbers.REAL_PART(c);
    cim : constant triple_double := TripDobl_Complex_Numbers.IMAG_PART(c);
    rre : constant double_float := to_double(cre);
    rim : constant double_float := to_double(cim);

  begin
    res := Standard_Complex_Numbers.Create(rre,rim);
    return res;
  end TripDobl_Complex_to_Standard;

  function TripDobl_Complex_to_Multprec
             ( c : TripDobl_Complex_Numbers.Complex_Number )
             return Multprec_Complex_Numbers.Complex_Number is

    res : Multprec_Complex_Numbers.Complex_Number;
    cre : constant triple_double := TripDobl_Complex_Numbers.REAL_PART(c);
    cim : constant triple_double := TripDobl_Complex_Numbers.IMAG_PART(c);
    rre : Floating_Number := to_floating_number(cre);
    rim : Floating_Number := to_floating_number(cim);

  begin
    res := Multprec_Complex_Numbers.Create(rre,rim);
    Multprec_Floating_Numbers.Clear(rre);
    Multprec_Floating_Numbers.Clear(rim);
    return res;
  end TripDobl_Complex_to_Multprec;

end TripDobl_Complex_Numbers_cv;
