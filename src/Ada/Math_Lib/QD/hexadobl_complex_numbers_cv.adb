with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Multprec_Floating_Numbers;          use Multprec_Floating_Numbers;
with Double_Double_Numbers;              use Double_Double_Numbers;
with Triple_Double_Numbers;              use Triple_Double_Numbers;
with Quad_Double_Numbers;                use Quad_Double_Numbers;
with Penta_Double_Numbers;               use Penta_Double_Numbers;
with Octo_Double_Numbers;                use Octo_Double_Numbers;
with Deca_Double_Numbers;                use Deca_Double_Numbers;
with Hexa_Double_Numbers;                use Hexa_Double_Numbers;
with Multprec_HexaDobl_Convertors;       use Multprec_HexaDobl_Convertors;

package body HexaDobl_Complex_Numbers_cv is

  function Standard_to_HexaDobl_Complex
             ( c : Standard_Complex_Numbers.Complex_Number )
             return HexaDobl_Complex_Numbers.Complex_Number is

    res : HexaDobl_Complex_Numbers.Complex_Number;
    cre : constant double_float := Standard_Complex_Numbers.REAL_PART(c);
    cim : constant double_float := Standard_Complex_Numbers.IMAG_PART(c);
    rre : constant hexa_double := create(cre);
    rim : constant hexa_double := create(cim);

  begin
    res := HexaDobl_Complex_Numbers.Create(rre,rim);
    return res;
  end Standard_to_HexaDobl_Complex;

  function Multprec_to_HexaDobl_Complex
             ( c : Multprec_Complex_Numbers.Complex_Number )
             return HexaDobl_Complex_Numbers.Complex_Number is

    res : HexaDobl_Complex_Numbers.Complex_Number;
    cre : Floating_Number := Multprec_Complex_Numbers.REAL_PART(c);
    cim : Floating_Number := Multprec_Complex_Numbers.IMAG_PART(c);
    rre : constant Hexa_double := to_hexa_double(cre);
    rim : constant Hexa_double := to_hexa_double(cim);

  begin
    Multprec_Floating_Numbers.Clear(cre);
    Multprec_Floating_Numbers.Clear(cim);
    res := HexaDobl_Complex_Numbers.Create(rre,rim);
    return res;
  end Multprec_to_HexaDobl_Complex;

  function HexaDobl_Complex_to_Standard
             ( c : hexaDobl_Complex_Numbers.Complex_Number )
             return Standard_Complex_Numbers.Complex_Number is

    res : Standard_Complex_Numbers.Complex_Number;
    cre : constant hexa_double := HexaDobl_Complex_Numbers.REAL_PART(c);
    cim : constant hexa_double := HexaDobl_Complex_Numbers.IMAG_PART(c);
    rre : constant double_float := to_double(cre);
    rim : constant double_float := to_double(cim);

  begin
    res := Standard_Complex_Numbers.Create(rre,rim);
    return res;
  end HexaDobl_Complex_to_Standard;

  function HexaDobl_Complex_to_Multprec
             ( c : HexaDobl_Complex_Numbers.Complex_Number )
             return Multprec_Complex_Numbers.Complex_Number is

    res : Multprec_Complex_Numbers.Complex_Number;
    cre : constant hexa_double := HexaDobl_Complex_Numbers.REAL_PART(c);
    cim : constant hexa_double := HexaDobl_Complex_Numbers.IMAG_PART(c);
    rre : Floating_Number := to_floating_number(cre);
    rim : Floating_Number := to_floating_number(cim);

  begin
    res := Multprec_Complex_Numbers.Create(rre,rim);
    Multprec_Floating_Numbers.Clear(rre);
    Multprec_Floating_Numbers.Clear(rim);
    return res;
  end HexaDobl_Complex_to_Multprec;

  function HexaDobl_Complex_to_DoblDobl
             ( c : HexaDobl_Complex_Numbers.Complex_Number )
             return DoblDobl_Complex_Numbers.Complex_Number is

    res : DoblDobl_Complex_Numbers.Complex_Number;
    cre : constant hexa_double := HexaDobl_Complex_Numbers.REAL_PART(c);
    cim : constant hexa_double := HexaDobl_Complex_Numbers.IMAG_PART(c);
    rre : constant double_double := to_double_double(cre);
    rim : constant double_double := to_double_double(cim);

  begin
    res := DoblDobl_Complex_Numbers.Create(rre,rim);
    return res;
  end HexaDobl_Complex_to_DoblDobl;

  function HexaDobl_Complex_to_TripDobl
             ( c : HexaDobl_Complex_Numbers.Complex_Number )
             return TripDobl_Complex_Numbers.Complex_Number is

    res : TripDobl_Complex_Numbers.Complex_Number;
    cre : constant hexa_double := HexaDobl_Complex_Numbers.REAL_PART(c);
    cim : constant hexa_double := HexaDobl_Complex_Numbers.IMAG_PART(c);
    rre : constant triple_double := to_triple_double(cre);
    rim : constant triple_double := to_triple_double(cim);

  begin
    res := TripDobl_Complex_Numbers.Create(rre,rim);
    return res;
  end HexaDobl_Complex_to_TripDobl;

  function HexaDobl_Complex_to_QuadDobl
             ( c : HexaDobl_Complex_Numbers.Complex_Number )
             return QuadDobl_Complex_Numbers.Complex_Number is

    res : QuadDobl_Complex_Numbers.Complex_Number;
    cre : constant hexa_double := HexaDobl_Complex_Numbers.REAL_PART(c);
    cim : constant hexa_double := HexaDobl_Complex_Numbers.IMAG_PART(c);
    rre : constant quad_double := to_quad_double(cre);
    rim : constant quad_double := to_quad_double(cim);

  begin
    res := QuadDobl_Complex_Numbers.Create(rre,rim);
    return res;
  end HexaDobl_Complex_to_QuadDobl;

  function HexaDobl_Complex_to_PentDobl
             ( c : HexaDobl_Complex_Numbers.Complex_Number )
             return PentDobl_Complex_Numbers.Complex_Number is

    res : PentDobl_Complex_Numbers.Complex_Number;
    cre : constant hexa_double := HexaDobl_Complex_Numbers.REAL_PART(c);
    cim : constant hexa_double := HexaDobl_Complex_Numbers.IMAG_PART(c);
    rre : constant penta_double := to_penta_double(cre);
    rim : constant penta_double := to_penta_double(cim);

  begin
    res := PentDobl_Complex_Numbers.Create(rre,rim);
    return res;
  end HexaDobl_Complex_to_PentDobl;

  function HexaDobl_Complex_to_OctoDobl
             ( c : HexaDobl_Complex_Numbers.Complex_Number )
             return OctoDobl_Complex_Numbers.Complex_Number is

    res : OctoDobl_Complex_Numbers.Complex_Number;
    cre : constant hexa_double := HexaDobl_Complex_Numbers.REAL_PART(c);
    cim : constant hexa_double := HexaDobl_Complex_Numbers.IMAG_PART(c);
    rre : constant octo_double := to_octo_double(cre);
    rim : constant octo_double := to_octo_double(cim);

  begin
    res := OctoDobl_Complex_Numbers.Create(rre,rim);
    return res;
  end HexaDobl_Complex_to_OctoDobl;

  function HexaDobl_Complex_to_DecaDobl
             ( c : HexaDobl_Complex_Numbers.Complex_Number )
             return DecaDobl_Complex_Numbers.Complex_Number is

    res : DecaDobl_Complex_Numbers.Complex_Number;
    cre : constant hexa_double := HexaDobl_Complex_Numbers.REAL_PART(c);
    cim : constant hexa_double := HexaDobl_Complex_Numbers.IMAG_PART(c);
    rre : constant deca_double := to_deca_double(cre);
    rim : constant deca_double := to_deca_double(cim);

  begin
    res := DecaDobl_Complex_Numbers.Create(rre,rim);
    return res;
  end HexaDobl_Complex_to_DecaDobl;

end HexaDobl_Complex_Numbers_cv;
