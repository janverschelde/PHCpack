with text_io;                            use text_io;
with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;
with Double_Double_Numbers_io;           use Double_Double_Numbers_io;
with Quad_Double_Numbers_io;             use Quad_Double_Numbers_io;
with Octo_Double_Numbers_io;             use Octo_Double_Numbers_io;
with Hexa_Double_Numbers_io;             use Hexa_Double_Numbers_io;
with DoblDobl_Complex_Numbers;
with QuadDobl_Complex_Numbers;
with OctoDobl_Complex_Numbers;
with HexaDobl_Complex_Numbers;
with DoblDobl_Random_Numbers;
with QuadDobl_Random_Numbers;
with OctoDobl_Random_Numbers;
with HexaDobl_Random_Numbers;
with Sign_Balancers;                     use Sign_Balancers;

package body Test_Sign_Balancers is

  procedure Test_Sign_Balance ( nbr : in double_double ) is

    x : double_double := nbr;
    xb,err : double_double;

  begin
    put("x : "); put(x);
    if Is_Sign_Balanced(x) then
      put_line(" sign balanced");
    else
      put_line(" not sign balanced");
      put("x hi : "); put(hi_part(x)); new_line;
      put("x lo : "); put(lo_part(x)); new_line;
      xb := nbr;
      Sign_Balance(x);
      put("x hi : "); put(hi_part(x)); new_line;
      put("x lo : "); put(lo_part(x)); new_line;
      put("org x : "); put(xb); new_line;
      put("new x : "); put(x); 
      if Is_Sign_Balanced(x)
       then put_line(" sign balanced");
       else put_line(" NOT sign balanced, bug!");
      end if;
      err := abs(xb - x);
      put("error : "); put(err,2); new_line;
    end if;
  end Test_Sign_Balance;

  procedure Test_Sign_Balance ( nbr : in quad_double ) is

    x : quad_double := nbr;
    xb,err : quad_double;

  begin
    put("x : "); put(x);
    if Is_Sign_Balanced(x) then
      put_line(" sign balanced");
    else
      put_line(" not sign balanced");
      put("x hihi : "); put(hihi_part(x)); new_line;
      put("x lohi : "); put(lohi_part(x)); new_line;
      put("x hilo : "); put(hilo_part(x)); new_line;
      put("x lolo : "); put(lolo_part(x)); new_line;
      xb := nbr;
      Sign_Balance(x);
      put("x hihi : "); put(hihi_part(x)); new_line;
      put("x lohi : "); put(lohi_part(x)); new_line;
      put("x hilo : "); put(hilo_part(x)); new_line;
      put("x lolo : "); put(lolo_part(x)); new_line;
      put("org x : "); put(xb); new_line;
      put("new x : "); put(x); 
      if Is_Sign_Balanced(x)
       then put_line(" sign balanced");
       else put_line(" NOT sign balanced, bug!");
      end if;
      err := abs(xb - x);
      put("error : "); put(err,2); new_line;
    end if;
  end Test_Sign_Balance;

  procedure Test_Sign_Balance ( nbr : in octo_double ) is

    x : octo_double := nbr;
    xb,err : octo_double;

  begin
    put("x : "); put(x);
    if Is_Sign_Balanced(x) then
      put_line(" sign balanced");
    else
      put_line(" not sign balanced");
      put("x hihihi : "); put(hihihi_part(x)); new_line;
      put("x lohihi : "); put(lohihi_part(x)); new_line;
      put("x hilohi : "); put(hilohi_part(x)); new_line;
      put("x lolohi : "); put(lolohi_part(x)); new_line;
      put("x hihilo : "); put(hihilo_part(x)); new_line;
      put("x lohilo : "); put(lohilo_part(x)); new_line;
      put("x hilolo : "); put(hilolo_part(x)); new_line;
      put("x lololo : "); put(lololo_part(x)); new_line;
      xb := nbr;
      Sign_Balance(x);
      put("x hihihi : "); put(hihihi_part(x)); new_line;
      put("x lohihi : "); put(lohihi_part(x)); new_line;
      put("x hilohi : "); put(hilohi_part(x)); new_line;
      put("x lolohi : "); put(lolohi_part(x)); new_line;
      put("x hihilo : "); put(hihilo_part(x)); new_line;
      put("x lohilo : "); put(lohilo_part(x)); new_line;
      put("x hilolo : "); put(hilolo_part(x)); new_line;
      put("x lololo : "); put(lololo_part(x)); new_line;
      put("org x : "); put(xb); new_line;
      put("new x : "); put(x); 
      if Is_Sign_Balanced(x)
       then put_line(" sign balanced");
       else put_line(" NOT sign balanced, bug!");
      end if;
      err := abs(xb - x);
      put("error : "); put(err,2); new_line;
    end if;
  end Test_Sign_Balance;

  procedure Test_Sign_Balance ( nbr : in hexa_double ) is

    x : hexa_double := nbr;
    xb,err : hexa_double;

  begin
    put("x : "); put(x);
    if Is_Sign_Balanced(x) then
      put_line(" sign balanced");
    else
      put_line(" not sign balanced");
      put("x hihihihi : "); put(hihihihi_part(x)); new_line;
      put("x lohihihi : "); put(lohihihi_part(x)); new_line;
      put("x hilohihi : "); put(hilohihi_part(x)); new_line;
      put("x lolohihi : "); put(lolohihi_part(x)); new_line;
      put("x hihilohi : "); put(hihilohi_part(x)); new_line;
      put("x lohilohi : "); put(lohilohi_part(x)); new_line;
      put("x hilolohi : "); put(hilolohi_part(x)); new_line;
      put("x lololohi : "); put(lololohi_part(x)); new_line;
      put("x hihihilo : "); put(hihihilo_part(x)); new_line;
      put("x lohihilo : "); put(lohihilo_part(x)); new_line;
      put("x hilohilo : "); put(hilohilo_part(x)); new_line;
      put("x lolohilo : "); put(lolohilo_part(x)); new_line;
      put("x hihilolo : "); put(hihilolo_part(x)); new_line;
      put("x lohilolo : "); put(lohilolo_part(x)); new_line;
      put("x hilololo : "); put(hilololo_part(x)); new_line;
      put("x lolololo : "); put(lolololo_part(x)); new_line;
      xb := nbr;
      Sign_Balance(x);
      put("x hihihihi : "); put(hihihihi_part(x)); new_line;
      put("x lohihihi : "); put(lohihihi_part(x)); new_line;
      put("x hilohihi : "); put(hilohihi_part(x)); new_line;
      put("x lolohihi : "); put(lolohihi_part(x)); new_line;
      put("x hihilohi : "); put(hihilohi_part(x)); new_line;
      put("x lohilohi : "); put(lohilohi_part(x)); new_line;
      put("x hilolohi : "); put(hilolohi_part(x)); new_line;
      put("x lololohi : "); put(lololohi_part(x)); new_line;
      put("x hihihilo : "); put(hihihilo_part(x)); new_line;
      put("x lohihilo : "); put(lohihilo_part(x)); new_line;
      put("x hilohilo : "); put(hilohilo_part(x)); new_line;
      put("x lolohilo : "); put(lolohilo_part(x)); new_line;
      put("x hihilolo : "); put(hihilolo_part(x)); new_line;
      put("x lohilolo : "); put(lohilolo_part(x)); new_line;
      put("x hilololo : "); put(hilololo_part(x)); new_line;
      put("x lolololo : "); put(lolololo_part(x)); new_line;
      put("org x : "); put(xb); new_line;
      put("new x : "); put(x); 
      if Is_Sign_Balanced(x)
       then put_line(" sign balanced");
       else put_line(" NOT sign balanced, bug!");
      end if;
      err := abs(xb - x);
      put("error : "); put(err,2); new_line;
    end if;
  end Test_Sign_Balance;

  procedure Test_Sign_DD_Balance is

    rnd : constant DoblDobl_Complex_Numbers.Complex_Number
        := DoblDobl_Random_Numbers.Random1;
    x : constant double_double := DoblDobl_Complex_Numbers.REAL_PART(rnd);
    y : constant double_double := DoblDobl_Complex_Numbers.IMAG_PART(rnd);

  begin
    Test_Sign_Balance(x);
    Test_Sign_Balance(y);
  end Test_Sign_DD_Balance;

  procedure Test_Sign_QD_Balance is

    rnd : constant QuadDobl_Complex_Numbers.Complex_Number
        := QuadDobl_Random_Numbers.Random1;
    x : constant quad_double := QuadDobl_Complex_Numbers.REAL_PART(rnd);
    y : constant quad_double := QuadDobl_Complex_Numbers.IMAG_PART(rnd);

  begin
    Test_Sign_Balance(x);
    Test_Sign_Balance(y);
  end Test_Sign_QD_Balance;

  procedure Test_Sign_OD_Balance is

    rnd : constant OctoDobl_Complex_Numbers.Complex_Number
        := OctoDobl_Random_Numbers.Random1;
    x : constant octo_double := OctoDobl_Complex_Numbers.REAL_PART(rnd);
    y : constant octo_double := OctoDobl_Complex_Numbers.IMAG_PART(rnd);

  begin
    Test_Sign_Balance(x);
    Test_Sign_Balance(y);
  end Test_Sign_OD_Balance;

  procedure Test_Sign_HD_Balance is

    rnd : constant HexaDobl_Complex_Numbers.Complex_Number
        := HexaDobl_Random_Numbers.Random1;
    x : constant hexa_double := HexaDobl_Complex_Numbers.REAL_PART(rnd);
    y : constant hexa_double := HexaDobl_Complex_Numbers.IMAG_PART(rnd);

  begin
    Test_Sign_Balance(x);
    Test_Sign_Balance(y);
  end Test_Sign_HD_Balance;

  procedure Main is
  begin
    Test_Sign_DD_Balance;
    new_line;
    Test_Sign_QD_Balance;
    new_line;
    Test_Sign_OD_Balance;
    new_line;
    Test_Sign_HD_Balance;
  end Main;

end Test_Sign_Balancers;
