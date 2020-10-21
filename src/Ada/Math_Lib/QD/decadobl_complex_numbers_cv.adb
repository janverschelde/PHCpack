with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Multprec_Floating_Numbers;          use Multprec_Floating_Numbers;
with Double_Double_Numbers;              use Double_Double_Numbers;
with Triple_Double_Numbers;              use Triple_Double_Numbers;
with Quad_Double_Numbers;                use Quad_Double_Numbers;
with Penta_Double_Numbers;               use Penta_Double_Numbers;
with Octo_Double_Numbers;                use Octo_Double_Numbers;
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

  function DecaDobl_Complex_to_DoblDobl
             ( c : DecaDobl_Complex_Numbers.Complex_Number )
             return DoblDobl_Complex_Numbers.Complex_Number is

    res : DoblDobl_Complex_Numbers.Complex_Number;
    cre : constant deca_double := DecaDobl_Complex_Numbers.REAL_PART(c);
    cim : constant deca_double := DecaDobl_Complex_Numbers.IMAG_PART(c);
    rre : constant double_double := to_double_double(cre);
    rim : constant double_double := to_double_double(cim);

  begin
    res := DoblDobl_Complex_Numbers.Create(rre,rim);
    return res;
  end DecaDobl_Complex_to_DoblDobl;

  function DecaDobl_Complex_to_TripDobl
             ( c : DecaDobl_Complex_Numbers.Complex_Number )
             return TripDobl_Complex_Numbers.Complex_Number is

    res : TripDobl_Complex_Numbers.Complex_Number;
    cre : constant deca_double := DecaDobl_Complex_Numbers.REAL_PART(c);
    cim : constant deca_double := DecaDobl_Complex_Numbers.IMAG_PART(c);
    rre : constant triple_double := to_triple_double(cre);
    rim : constant triple_double := to_triple_double(cim);

  begin
    res := TripDobl_Complex_Numbers.Create(rre,rim);
    return res;
  end DecaDobl_Complex_to_TripDobl;

  function DecaDobl_Complex_to_QuadDobl
             ( c : DecaDobl_Complex_Numbers.Complex_Number )
             return QuadDobl_Complex_Numbers.Complex_Number is

    res : QuadDobl_Complex_Numbers.Complex_Number;
    cre : constant deca_double := DecaDobl_Complex_Numbers.REAL_PART(c);
    cim : constant deca_double := DecaDobl_Complex_Numbers.IMAG_PART(c);
    rre : constant quad_double := to_quad_double(cre);
    rim : constant quad_double := to_quad_double(cim);

  begin
    res := QuadDobl_Complex_Numbers.Create(rre,rim);
    return res;
  end DecaDobl_Complex_to_QuadDobl;

  function DecaDobl_Complex_to_PentDobl
             ( c : DecaDobl_Complex_Numbers.Complex_Number )
             return PentDobl_Complex_Numbers.Complex_Number is

    res : PentDobl_Complex_Numbers.Complex_Number;
    cre : constant deca_double := DecaDobl_Complex_Numbers.REAL_PART(c);
    cim : constant deca_double := DecaDobl_Complex_Numbers.IMAG_PART(c);
    rre : constant penta_double := to_penta_double(cre);
    rim : constant penta_double := to_penta_double(cim);

  begin
    res := PentDobl_Complex_Numbers.Create(rre,rim);
    return res;
  end DecaDobl_Complex_to_PentDobl;

  function DecaDobl_Complex_to_OctoDobl
             ( c : DecaDobl_Complex_Numbers.Complex_Number )
             return OctoDobl_Complex_Numbers.Complex_Number is

    res : OctoDobl_Complex_Numbers.Complex_Number;
    cre : constant deca_double := DecaDobl_Complex_Numbers.REAL_PART(c);
    cim : constant deca_double := DecaDobl_Complex_Numbers.IMAG_PART(c);
    rre : constant octo_double := to_octo_double(cre);
    rim : constant octo_double := to_octo_double(cim);

  begin
    res := OctoDobl_Complex_Numbers.Create(rre,rim);
    return res;
  end DecaDobl_Complex_to_OctoDobl;

end DecaDobl_Complex_Numbers_cv;
