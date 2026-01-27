with text_io;                            use text_io;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;
with Bits_of_Doubles;
with Sign_Balancers;                     use Sign_Balancers;

package body Vectored_Octo_Doubles is

-- BASIC PROCEDURES :

  function Sign_Balance
             ( x : Octo_Double_Vectors.Vector; verbose : boolean := true )
             return Octo_Double_Vectors.Vector is

    res : Octo_Double_Vectors.Vector(x'range);

  begin
    for i in x'range loop
      res(i) := x(i);
      if not Is_Sign_Balanced(res(i)) then
        if verbose
         then Sign_Balance(res(i),vrblvl=>1);
         else Sign_Balance(res(i),vrblvl=>0);
        end if;
      end if;
    end loop;
    return res;
  end Sign_Balance;

  procedure Signed_Quarter
              ( x,y : in Octo_Double_Vectors.Vector;
                xs00,xs01,xs02,xs03 : out Standard_Floating_Vectors.Vector;
                xs04,xs05,xs06,xs07 : out Standard_Floating_Vectors.Vector;
                xs08,xs09,xs10,xs11 : out Standard_Floating_Vectors.Vector;
                xs12,xs13,xs14,xs15 : out Standard_Floating_Vectors.Vector;
                xs16,xs17,xs18,xs19 : out Standard_Floating_Vectors.Vector;
                xs20,xs21,xs22,xs23 : out Standard_Floating_Vectors.Vector;
                xs24,xs25,xs26,xs27 : out Standard_Floating_Vectors.Vector;
                xs28,xs29,xs30,xs31 : out Standard_Floating_Vectors.Vector;
                ys00,ys01,ys02,ys03 : out Standard_Floating_Vectors.Vector;
                ys04,ys05,ys06,ys07 : out Standard_Floating_Vectors.Vector;
                ys08,ys09,ys10,ys11 : out Standard_Floating_Vectors.Vector;
                ys12,ys13,ys14,ys15 : out Standard_Floating_Vectors.Vector;
                ys16,ys17,ys18,ys19 : out Standard_Floating_Vectors.Vector;
                ys20,ys21,ys22,ys23 : out Standard_Floating_Vectors.Vector;
                ys24,ys25,ys26,ys27 : out Standard_Floating_Vectors.Vector;
                ys28,ys29,ys30,ys31 : out Standard_Floating_Vectors.Vector;
                xd00,xd01,xd02,xd03 : out Standard_Floating_Vectors.Vector;
                xd04,xd05,xd06,xd07 : out Standard_Floating_Vectors.Vector;
                xd08,xd09,xd10,xd11 : out Standard_Floating_Vectors.Vector;
                xd12,xd13,xd14,xd15 : out Standard_Floating_Vectors.Vector;
                xd16,xd17,xd18,xd19 : out Standard_Floating_Vectors.Vector;
                xd20,xd21,xd22,xd23 : out Standard_Floating_Vectors.Vector;
                xd24,xd25,xd26,xd27 : out Standard_Floating_Vectors.Vector;
                xd28,xd29,xd30,xd31 : out Standard_Floating_Vectors.Vector;
                yd00,yd01,yd02,yd03 : out Standard_Floating_Vectors.Vector;
                yd04,yd05,yd06,yd07 : out Standard_Floating_Vectors.Vector;
                yd08,yd09,yd10,yd11 : out Standard_Floating_Vectors.Vector;
                yd12,yd13,yd14,yd15 : out Standard_Floating_Vectors.Vector;
                yd16,yd17,yd18,yd19 : out Standard_Floating_Vectors.Vector;
                yd20,yd21,yd22,yd23 : out Standard_Floating_Vectors.Vector;
                yd24,yd25,yd26,yd27 : out Standard_Floating_Vectors.Vector;
                yd28,yd29,yd30,yd31 : out Standard_Floating_Vectors.Vector;
                ns,nd : out integer32 ) is

    nbr : octo_double;
    flt : double_float;
    x00,x01,x02,x03,x04,x05,x06,x07 : double_float;
    x08,x09,x10,x11,x12,x13,x14,x15 : double_float;
    x16,x17,x18,x19,x20,x21,x22,x23 : double_float;
    x24,x25,x26,x27,x28,x29,x30,x31 : double_float;
    y00,y01,y02,y03,y04,y05,y06,y07 : double_float;
    y08,y09,y10,y11,y12,y13,y14,y15 : double_float;
    y16,y17,y18,y19,y20,y21,y22,y23 : double_float;
    y24,y25,y26,y27,y28,y29,y30,y31 : double_float;

  begin
    ns := 0; nd := 0;
    for i in x'range loop
      nbr := x(i);
      flt := hihihi_part(nbr); Bits_of_Doubles.Split(flt,x00,x01,x02,x03);
      flt := lohihi_part(nbr); Bits_of_Doubles.Split(flt,x04,x05,x06,x07);
      flt := hilohi_part(nbr); Bits_of_Doubles.Split(flt,x08,x09,x10,x11);
      flt := lolohi_part(nbr); Bits_of_Doubles.Split(flt,x12,x13,x14,x15);
      flt := hihilo_part(nbr); Bits_of_Doubles.Split(flt,x16,x17,x18,x19);
      flt := lohilo_part(nbr); Bits_of_Doubles.Split(flt,x20,x21,x22,x23);
      flt := hilolo_part(nbr); Bits_of_Doubles.Split(flt,x24,x25,x26,x27);
      flt := lololo_part(nbr); Bits_of_Doubles.Split(flt,x28,x29,x30,x31);
      nbr := y(i);
      flt := hihihi_part(nbr); Bits_of_Doubles.Split(flt,y00,y01,y02,y03);
      flt := lohihi_part(nbr); Bits_of_Doubles.Split(flt,y04,y05,y06,y07);
      flt := hilohi_part(nbr); Bits_of_Doubles.Split(flt,y08,y09,y10,y11);
      flt := lolohi_part(nbr); Bits_of_Doubles.Split(flt,y12,y13,y14,y15);
      flt := hihilo_part(nbr); Bits_of_Doubles.Split(flt,y16,y17,y18,y19);
      flt := lohilo_part(nbr); Bits_of_Doubles.Split(flt,y20,y21,y22,y23);
      flt := hilolo_part(nbr); Bits_of_Doubles.Split(flt,y24,y25,y26,y27);
      flt := lololo_part(nbr); Bits_of_Doubles.Split(flt,y28,y29,y30,y31);
      if not Sign_Balancers.Different_Sign(x00,y00) then
        ns := ns + 1;
        xs00(ns) := x00; xs01(ns) := x01;
        xs02(ns) := x02; xs03(ns) := x03;
        xs04(ns) := x04; xs05(ns) := x05;
        xs06(ns) := x06; xs07(ns) := x07;
        xs08(ns) := x08; xs09(ns) := x09;
        xs10(ns) := x10; xs11(ns) := x11;
        xs12(ns) := x12; xs13(ns) := x13;
        xs14(ns) := x14; xs15(ns) := x15;
        xs16(ns) := x16; xs17(ns) := x17;
        xs18(ns) := x18; xs19(ns) := x19;
        xs20(ns) := x20; xs21(ns) := x21;
        xs22(ns) := x22; xs23(ns) := x23;
        xs24(ns) := x24; xs25(ns) := x25;
        xs26(ns) := x26; xs27(ns) := x27;
        xs28(ns) := x28; xs29(ns) := x29;
        xs30(ns) := x30; xs31(ns) := x31;
        ys00(ns) := y00; ys01(ns) := y01;
        ys02(ns) := y02; ys03(ns) := y03;
        ys04(ns) := y04; ys05(ns) := y05;
        ys06(ns) := y06; ys07(ns) := y07;
        ys08(ns) := y08; ys09(ns) := y09;
        ys10(ns) := y10; ys11(ns) := y11;
        ys12(ns) := y12; ys13(ns) := y13;
        ys14(ns) := y14; ys15(ns) := y15;
        ys16(ns) := y16; ys17(ns) := y17;
        ys18(ns) := y18; ys19(ns) := y19;
        ys20(ns) := y20; ys21(ns) := y21;
        ys22(ns) := y22; ys23(ns) := y23;
        ys24(ns) := y24; ys25(ns) := y25;
        ys26(ns) := y26; ys27(ns) := y27;
        ys28(ns) := y28; ys29(ns) := y29;
        ys30(ns) := y30; ys31(ns) := y31;
      else -- Different_Sign(x00,y00)
        nd := nd + 1;
        xd00(nd) := x00; xd01(nd) := x01;
        xd02(nd) := x02; xd03(nd) := x03;
        xd04(nd) := x04; xd05(nd) := x05;
        xd06(nd) := x06; xd07(nd) := x07;
        xd08(nd) := x08; xd09(nd) := x09;
        xd10(nd) := x10; xd11(nd) := x11;
        xd12(nd) := x12; xd13(nd) := x13;
        xd14(nd) := x14; xd15(nd) := x15;
        xd16(nd) := x16; xd17(nd) := x17;
        xd18(nd) := x18; xd19(nd) := x19;
        xd20(nd) := x20; xd21(nd) := x21;
        xd22(nd) := x22; xd23(nd) := x23;
        xd24(nd) := x24; xd25(nd) := x25;
        xd26(nd) := x26; xd27(nd) := x27;
        xd28(nd) := x28; xd29(nd) := x29;
        xd30(nd) := x30; xd31(nd) := x31;
        yd00(nd) := y00; yd01(nd) := y01;
        yd02(nd) := y02; yd03(nd) := y03;
        yd04(nd) := y04; yd05(nd) := y05;
        yd06(nd) := y06; yd07(nd) := y07;
        yd08(nd) := y08; yd09(nd) := y09;
        yd10(nd) := y10; yd11(nd) := y11;
        yd12(nd) := y12; yd13(nd) := y13;
        yd14(nd) := y14; yd15(nd) := y15;
        yd16(nd) := y16; yd17(nd) := y17;
        yd18(nd) := y18; yd19(nd) := y19;
        yd20(nd) := y20; yd21(nd) := y21;
        yd22(nd) := y22; yd23(nd) := y23;
        yd24(nd) := y24; yd25(nd) := y25;
        yd26(nd) := y26; yd27(nd) := y27;
        yd28(nd) := y28; yd29(nd) := y29;
        yd30(nd) := y30; yd31(nd) := y31;
      end if;
    end loop;
  end Signed_Quarter;

  procedure Balanced_Quarter_Product
              ( dim : in integer32;
                x00,x01,x02,x03 : in Standard_Floating_Vectors.Vector;
                x04,x05,x06,x07 : in Standard_Floating_Vectors.Vector;
                x08,x09,x10,x11 : in Standard_Floating_Vectors.Vector;
                x12,x13,x14,x15 : in Standard_Floating_Vectors.Vector;
                x16,x17,x18,x19 : in Standard_Floating_Vectors.Vector;
                x20,x21,x22,x23 : in Standard_Floating_Vectors.Vector;
                x24,x25,x26,x27 : in Standard_Floating_Vectors.Vector;
                x28,x29,x30,x31 : in Standard_Floating_Vectors.Vector;
                y00,y01,y02,y03 : in Standard_Floating_Vectors.Vector;
                y04,y05,y06,y07 : in Standard_Floating_Vectors.Vector;
                y08,y09,y10,y11 : in Standard_Floating_Vectors.Vector;
                y12,y13,y14,y15 : in Standard_Floating_Vectors.Vector;
                y16,y17,y18,y19 : in Standard_Floating_Vectors.Vector;
                y20,y21,y22,y23 : in Standard_Floating_Vectors.Vector;
                y24,y25,y26,y27 : in Standard_Floating_Vectors.Vector;
                y28,y29,y30,y31 : in Standard_Floating_Vectors.Vector;
                s00,s01,s02,s03,s04,s05,s06,s07 : out double_float;
                s08,s09,s10,s11,s12,s13,s14,s15 : out double_float;
                s16,s17,s18,s19,s20,s21,s22,s23 : out double_float;
                s24,s25,s26,s27,s28,s29,s30,s31 : out double_float ) is
  begin
    s00 := 0.0; s01 := 0.0; s02 := 0.0; s03 := 0.0;
    s04 := 0.0; s05 := 0.0; s06 := 0.0; s07 := 0.0;
    s08 := 0.0; s09 := 0.0; s10 := 0.0; s11 := 0.0;
    s12 := 0.0; s13 := 0.0; s14 := 0.0; s15 := 0.0;
    s16 := 0.0; s17 := 0.0; s18 := 0.0; s19 := 0.0;
    s20 := 0.0; s21 := 0.0; s22 := 0.0; s23 := 0.0;
    s24 := 0.0; s25 := 0.0; s26 := 0.0; s27 := 0.0;
    s28 := 0.0; s29 := 0.0; s30 := 0.0; s31 := 0.0;
    for i in 1..dim loop
      s00 := s00 + x00(i)*y00(i);
      s01 := s01 + x00(i)*y01(i) + x01(i)*y00(i);
      s02 := s02 + x00(i)*y02(i) + x01(i)*y01(i) + x02(i)*y00(i);
      s03 := s03 + x00(i)*y03(i) + x01(i)*y02(i) + x02(i)*y01(i)
                 + x03(i)*y00(i);
      s04 := s04 + x00(i)*y04(i) + x01(i)*y03(i) + x02(i)*y02(i)
                 + x03(i)*y01(i) + x04(i)*y00(i);
      s05 := s05 + x00(i)*y05(i) + x01(i)*y04(i) + x02(i)*y03(i)
                 + x03(i)*y02(i) + x04(i)*y01(i) + x05(i)*y00(i);
      s06 := s06 + x00(i)*y06(i) + x01(i)*y05(i) + x02(i)*y04(i)
                 + x03(i)*y03(i) + x04(i)*y02(i) + x05(i)*y01(i)
                 + x06(i)*y00(i);
      s07 := s07 + x00(i)*y07(i) + x01(i)*y06(i) + x02(i)*y05(i)
                 + x03(i)*y04(i) + x04(i)*y03(i) + x05(i)*y02(i)
                 + x06(i)*y01(i) + x07(i)*y00(i);
      s08 := s08 + x00(i)*y08(i) + x01(i)*y07(i) + x02(i)*y06(i)
                 + x03(i)*y05(i) + x04(i)*y04(i) + x05(i)*y03(i)
                 + x06(i)*y02(i) + x07(i)*y01(i) + x08(i)*y00(i);
      s09 := s09 + x00(i)*y09(i) + x01(i)*y08(i) + x02(i)*y07(i)
                 + x03(i)*y06(i) + x04(i)*y05(i) + x05(i)*y04(i)
                 + x06(i)*y03(i) + x07(i)*y02(i) + x08(i)*y01(i)
                 + x09(i)*y00(i);
      s10 := s10 + x00(i)*y10(i) + x01(i)*y09(i) + x02(i)*y08(i)
                 + x03(i)*y07(i) + x04(i)*y06(i) + x05(i)*y05(i)
                 + x06(i)*y04(i) + x07(i)*y03(i) + x08(i)*y02(i)
                 + x09(i)*y01(i) + x10(i)*y00(i);
      s11 := s11 + x00(i)*y11(i) + x01(i)*y10(i) + x02(i)*y09(i)
                 + x03(i)*y08(i) + x04(i)*y07(i) + x05(i)*y06(i)
                 + x06(i)*y05(i) + x07(i)*y04(i) + x08(i)*y03(i)
                 + x09(i)*y02(i) + x10(i)*y01(i) + x11(i)*y00(i);
      s12 := s12 + x00(i)*y12(i) + x01(i)*y11(i) + x02(i)*y10(i)
                 + x03(i)*y09(i) + x04(i)*y08(i) + x05(i)*y07(i)
                 + x06(i)*y06(i) + x07(i)*y05(i) + x08(i)*y04(i)
                 + x09(i)*y03(i) + x10(i)*y02(i) + x11(i)*y01(i)
                 + x12(i)*y00(i);
      s13 := s13 + x00(i)*y13(i) + x01(i)*y12(i) + x02(i)*y11(i)
                 + x03(i)*y10(i) + x04(i)*y09(i) + x05(i)*y08(i)
                 + x06(i)*y07(i) + x07(i)*y06(i) + x08(i)*y05(i)
                 + x09(i)*y04(i) + x10(i)*y03(i) + x11(i)*y02(i)
                 + x12(i)*y01(i) + x13(i)*y00(i);
      s14 := s14 + x00(i)*y14(i) + x01(i)*y13(i) + x02(i)*y12(i)
                 + x03(i)*y11(i) + x04(i)*y10(i) + x05(i)*y09(i)
                 + x06(i)*y08(i) + x07(i)*y07(i) + x08(i)*y06(i)
                 + x09(i)*y05(i) + x10(i)*y04(i) + x11(i)*y03(i)
                 + x12(i)*y02(i) + x13(i)*y01(i) + x14(i)*y00(i);
      s15 := s15 + x00(i)*y15(i) + x01(i)*y14(i) + x02(i)*y13(i)
                 + x03(i)*y12(i) + x04(i)*y11(i) + x05(i)*y10(i)
                 + x06(i)*y09(i) + x07(i)*y08(i) + x08(i)*y07(i)
                 + x09(i)*y06(i) + x10(i)*y05(i) + x11(i)*y04(i)
                 + x12(i)*y03(i) + x13(i)*y02(i) + x14(i)*y01(i)
                 + x15(i)*y00(i);
      s16 := s16 + x00(i)*y16(i) + x01(i)*y15(i) + x02(i)*y14(i)
                 + x03(i)*y13(i) + x04(i)*y12(i) + x05(i)*y11(i)
                 + x06(i)*y10(i) + x07(i)*y09(i) + x08(i)*y08(i)
                 + x09(i)*y07(i) + x10(i)*y06(i) + x11(i)*y05(i)
                 + x12(i)*y04(i) + x13(i)*y03(i) + x14(i)*y02(i)
                 + x15(i)*y01(i) + x16(i)*y00(i);
      s17 := s17 + x00(i)*y17(i) + x01(i)*y16(i) + x02(i)*y15(i)
                 + x03(i)*y14(i) + x04(i)*y13(i) + x05(i)*y12(i)
                 + x06(i)*y11(i) + x07(i)*y10(i) + x08(i)*y09(i)
                 + x09(i)*y08(i) + x10(i)*y07(i) + x11(i)*y06(i)
                 + x12(i)*y05(i) + x13(i)*y04(i) + x14(i)*y03(i)
                 + x15(i)*y02(i) + x16(i)*y01(i) + x17(i)*y00(i);
      s18 := s18 + x00(i)*y18(i) + x01(i)*y17(i) + x02(i)*y16(i)
                 + x03(i)*y15(i) + x04(i)*y14(i) + x05(i)*y13(i)
                 + x06(i)*y12(i) + x07(i)*y11(i) + x08(i)*y10(i)
                 + x09(i)*y09(i) + x10(i)*y08(i) + x11(i)*y07(i)
                 + x12(i)*y06(i) + x13(i)*y05(i) + x14(i)*y04(i)
                 + x15(i)*y03(i) + x16(i)*y02(i) + x17(i)*y01(i)
                 + x18(i)*y00(i);
      s19 := s19 + x00(i)*y19(i) + x01(i)*y18(i) + x02(i)*y17(i)
                 + x03(i)*y16(i) + x04(i)*y15(i) + x05(i)*y14(i)
                 + x06(i)*y13(i) + x07(i)*y12(i) + x08(i)*y11(i)
                 + x09(i)*y10(i) + x10(i)*y09(i) + x11(i)*y08(i)
                 + x12(i)*y07(i) + x13(i)*y06(i) + x14(i)*y05(i)
                 + x15(i)*y04(i) + x16(i)*y03(i) + x17(i)*y02(i)
                 + x18(i)*y01(i) + x19(i)*y00(i);
      s20 := s20 + x00(i)*y20(i) + x01(i)*y19(i) + x02(i)*y18(i)
                 + x03(i)*y17(i) + x04(i)*y16(i) + x05(i)*y15(i)
                 + x06(i)*y14(i) + x07(i)*y13(i) + x08(i)*y12(i)
                 + x09(i)*y11(i) + x10(i)*y10(i) + x11(i)*y09(i)
                 + x12(i)*y08(i) + x13(i)*y07(i) + x14(i)*y06(i)
                 + x15(i)*y05(i) + x16(i)*y04(i) + x17(i)*y03(i)
                 + x18(i)*y02(i) + x19(i)*y01(i) + x20(i)*y00(i);
      s21 := s21 + x00(i)*y21(i) + x01(i)*y20(i) + x02(i)*y19(i)
                 + x03(i)*y18(i) + x04(i)*y17(i) + x05(i)*y16(i)
                 + x06(i)*y15(i) + x07(i)*y14(i) + x08(i)*y13(i)
                 + x09(i)*y12(i) + x10(i)*y11(i) + x11(i)*y10(i)
                 + x12(i)*y09(i) + x13(i)*y08(i) + x14(i)*y07(i)
                 + x15(i)*y06(i) + x16(i)*y05(i) + x17(i)*y04(i)
                 + x18(i)*y03(i) + x19(i)*y02(i) + x20(i)*y01(i)
                 + x21(i)*y00(i);
      s22 := s22 + x00(i)*y22(i) + x01(i)*y21(i) + x02(i)*y20(i)
                 + x03(i)*y19(i) + x04(i)*y18(i) + x05(i)*y17(i)
                 + x06(i)*y16(i) + x07(i)*y15(i) + x08(i)*y14(i)
                 + x09(i)*y13(i) + x10(i)*y12(i) + x11(i)*y11(i)
                 + x12(i)*y10(i) + x13(i)*y09(i) + x14(i)*y08(i)
                 + x15(i)*y07(i) + x16(i)*y06(i) + x17(i)*y05(i)
                 + x18(i)*y04(i) + x19(i)*y03(i) + x20(i)*y02(i)
                 + x21(i)*y01(i) + x22(i)*y00(i);
      s23 := s23 + x00(i)*y23(i) + x01(i)*y22(i) + x02(i)*y21(i)
                 + x03(i)*y20(i) + x04(i)*y19(i) + x05(i)*y18(i)
                 + x06(i)*y17(i) + x07(i)*y16(i) + x08(i)*y15(i)
                 + x09(i)*y14(i) + x10(i)*y13(i) + x11(i)*y12(i)
                 + x12(i)*y11(i) + x13(i)*y10(i) + x14(i)*y09(i)
                 + x15(i)*y08(i) + x16(i)*y07(i) + x17(i)*y06(i)
                 + x18(i)*y05(i) + x19(i)*y04(i) + x20(i)*y03(i)
                 + x21(i)*y02(i) + x22(i)*y01(i) + x23(i)*y00(i);
      s24 := s24 + x00(i)*y24(i) + x01(i)*y23(i) + x02(i)*y22(i)
                 + x03(i)*y21(i) + x04(i)*y20(i) + x05(i)*y19(i)
                 + x06(i)*y18(i) + x07(i)*y17(i) + x08(i)*y16(i)
                 + x09(i)*y15(i) + x10(i)*y14(i) + x11(i)*y13(i)
                 + x12(i)*y12(i) + x13(i)*y11(i) + x14(i)*y10(i)
                 + x15(i)*y09(i) + x16(i)*y08(i) + x17(i)*y07(i)
                 + x18(i)*y06(i) + x19(i)*y05(i) + x20(i)*y04(i)
                 + x21(i)*y03(i) + x22(i)*y02(i) + x23(i)*y01(i)
                 + x24(i)*y00(i);
      s25 := s25 + x00(i)*y25(i) + x01(i)*y24(i) + x02(i)*y23(i)
                 + x03(i)*y22(i) + x04(i)*y21(i) + x05(i)*y20(i)
                 + x06(i)*y19(i) + x07(i)*y18(i) + x08(i)*y17(i)
                 + x09(i)*y16(i) + x10(i)*y15(i) + x11(i)*y14(i)
                 + x12(i)*y13(i) + x13(i)*y12(i) + x14(i)*y11(i)
                 + x15(i)*y10(i) + x16(i)*y09(i) + x17(i)*y08(i)
                 + x18(i)*y07(i) + x19(i)*y06(i) + x20(i)*y05(i)
                 + x21(i)*y04(i) + x22(i)*y03(i) + x23(i)*y02(i)
                 + x24(i)*y01(i) + x25(i)*y00(i);
      s26 := s26 + x00(i)*y26(i) + x01(i)*y25(i) + x02(i)*y24(i)
                 + x03(i)*y23(i) + x04(i)*y22(i) + x05(i)*y21(i)
                 + x06(i)*y20(i) + x07(i)*y19(i) + x08(i)*y18(i)
                 + x09(i)*y17(i) + x10(i)*y16(i) + x11(i)*y15(i)
                 + x12(i)*y14(i) + x13(i)*y13(i) + x14(i)*y12(i)
                 + x15(i)*y11(i) + x16(i)*y10(i) + x17(i)*y09(i)
                 + x18(i)*y08(i) + x19(i)*y07(i) + x20(i)*y06(i)
                 + x21(i)*y05(i) + x22(i)*y04(i) + x23(i)*y03(i)
                 + x24(i)*y02(i) + x25(i)*y01(i) + x26(i)*y00(i);
      s27 := s27 + x00(i)*y27(i) + x01(i)*y26(i) + x02(i)*y25(i)
                 + x03(i)*y24(i) + x04(i)*y23(i) + x05(i)*y22(i)
                 + x06(i)*y21(i) + x07(i)*y20(i) + x08(i)*y19(i)
                 + x09(i)*y18(i) + x10(i)*y17(i) + x11(i)*y16(i)
                 + x12(i)*y15(i) + x13(i)*y14(i) + x14(i)*y13(i)
                 + x15(i)*y12(i) + x16(i)*y11(i) + x17(i)*y10(i)
                 + x18(i)*y09(i) + x19(i)*y08(i) + x20(i)*y07(i)
                 + x21(i)*y06(i) + x22(i)*y05(i) + x23(i)*y04(i)
                 + x24(i)*y03(i) + x25(i)*y02(i) + x26(i)*y01(i)
                 + x27(i)*y00(i);
      s28 := s28 + x00(i)*y28(i) + x01(i)*y27(i) + x02(i)*y26(i)
                 + x03(i)*y25(i) + x04(i)*y24(i) + x05(i)*y23(i)
                 + x06(i)*y22(i) + x07(i)*y21(i) + x08(i)*y20(i)
                 + x09(i)*y19(i) + x10(i)*y18(i) + x11(i)*y17(i)
                 + x12(i)*y16(i) + x13(i)*y15(i) + x14(i)*y14(i)
                 + x15(i)*y13(i) + x16(i)*y12(i) + x17(i)*y11(i)
                 + x18(i)*y10(i) + x19(i)*y09(i) + x20(i)*y08(i)
                 + x21(i)*y07(i) + x22(i)*y06(i) + x23(i)*y05(i)
                 + x24(i)*y04(i) + x25(i)*y03(i) + x26(i)*y02(i)
                 + x27(i)*y01(i) + x28(i)*y00(i);
      s29 := s29 + x00(i)*y29(i) + x01(i)*y28(i) + x02(i)*y27(i)
                 + x03(i)*y26(i) + x04(i)*y25(i) + x05(i)*y24(i)
                 + x06(i)*y23(i) + x07(i)*y22(i) + x08(i)*y21(i)
                 + x09(i)*y20(i) + x10(i)*y19(i) + x11(i)*y18(i)
                 + x12(i)*y17(i) + x13(i)*y16(i) + x14(i)*y15(i)
                 + x15(i)*y14(i) + x16(i)*y13(i) + x17(i)*y12(i)
                 + x18(i)*y11(i) + x19(i)*y10(i) + x20(i)*y09(i)
                 + x21(i)*y08(i) + x22(i)*y07(i) + x23(i)*y06(i)
                 + x24(i)*y05(i) + x25(i)*y04(i) + x26(i)*y03(i)
                 + x27(i)*y02(i) + x28(i)*y01(i) + x29(i)*y00(i);
      s30 := s30 + x00(i)*y30(i) + x01(i)*y29(i) + x02(i)*y28(i)
                 + x03(i)*y27(i) + x04(i)*y26(i) + x05(i)*y25(i)
                 + x06(i)*y24(i) + x07(i)*y23(i) + x08(i)*y22(i)
                 + x09(i)*y21(i) + x10(i)*y20(i) + x11(i)*y19(i)
                 + x12(i)*y18(i) + x13(i)*y17(i) + x14(i)*y16(i)
                 + x15(i)*y15(i) + x16(i)*y14(i) + x17(i)*y13(i)
                 + x18(i)*y12(i) + x19(i)*y11(i) + x20(i)*y10(i)
                 + x21(i)*y09(i) + x22(i)*y08(i) + x23(i)*y07(i)
                 + x24(i)*y06(i) + x25(i)*y05(i) + x26(i)*y04(i)
                 + x27(i)*y03(i) + x28(i)*y02(i) + x29(i)*y01(i)
                 + x30(i)*y00(i);
      s31 := s31 + x00(i)*y31(i) + x01(i)*y30(i) + x02(i)*y29(i)
                 + x03(i)*y28(i) + x04(i)*y27(i) + x05(i)*y26(i)
                 + x06(i)*y25(i) + x07(i)*y24(i) + x08(i)*y23(i)
                 + x09(i)*y22(i) + x10(i)*y21(i) + x11(i)*y20(i)
                 + x12(i)*y19(i) + x13(i)*y18(i) + x14(i)*y17(i)
                 + x15(i)*y16(i) + x16(i)*y15(i) + x17(i)*y14(i)
                 + x18(i)*y13(i) + x19(i)*y12(i) + x20(i)*y11(i)
                 + x21(i)*y10(i) + x22(i)*y09(i) + x23(i)*y08(i)
                 + x24(i)*y07(i) + x25(i)*y06(i) + x26(i)*y05(i)
                 + x27(i)*y04(i) + x28(i)*y03(i) + x29(i)*y02(i)
                 + x30(i)*y01(i) + x31(i)*y00(i);
    end loop;
  end Balanced_Quarter_Product;

  procedure Write_Subsums
              ( s00,s01,s02,s03,s04,s05,s06,s07 : in double_float ) is

    use Bits_of_Doubles;

  begin
    put("s00 : "); put(s00);
    put(", n00 : "); put(Last_Zero_Count(s00),1); new_line;
    put("s01 : "); put(s01);
    put(", n01 : "); put(Last_Zero_Count(s01),1); new_line;
    put("s02 : "); put(s02);
    put(", n02 : "); put(Last_Zero_Count(s02),1); new_line;
    put("s03 : "); put(s03);
    put(", n03 : "); put(Last_Zero_Count(s03),1); new_line;
    put("s04 : "); put(s04);
    put(", n04 : "); put(Last_Zero_Count(s04),1); new_line;
    put("s05 : "); put(s05);
    put(", n05 : "); put(Last_Zero_Count(s05),1); new_line;
    put("s06 : "); put(s06);
    put(", n06 : "); put(Last_Zero_Count(s06),1); new_line;
    put("s07 : "); put(s07);
    put(", n07 : "); put(Last_Zero_Count(s07),1); new_line;
  end Write_Subsums;

  procedure Write_Subsums
              ( s00,s01,s02,s03,s04,s05,s06,s07 : in double_float;
                s08,s09,s10,s11,s12,s13,s14,s15 : in double_float ) is

    use Bits_of_Doubles;

  begin
    put("s00 : "); put(s00);
    put(", n00 : "); put(Last_Zero_Count(s00),1); new_line;
    put("s01 : "); put(s01);
    put(", n01 : "); put(Last_Zero_Count(s01),1); new_line;
    put("s02 : "); put(s02);
    put(", n02 : "); put(Last_Zero_Count(s02),1); new_line;
    put("s03 : "); put(s03);
    put(", n03 : "); put(Last_Zero_Count(s03),1); new_line;
    put("s04 : "); put(s04);
    put(", n04 : "); put(Last_Zero_Count(s04),1); new_line;
    put("s05 : "); put(s05);
    put(", n05 : "); put(Last_Zero_Count(s05),1); new_line;
    put("s06 : "); put(s06);
    put(", n06 : "); put(Last_Zero_Count(s06),1); new_line;
    put("s07 : "); put(s07);
    put(", n07 : "); put(Last_Zero_Count(s07),1); new_line;
    put("s08 : "); put(s08);
    put(", n08 : "); put(Last_Zero_Count(s08),1); new_line;
    put("s09 : "); put(s09);
    put(", n09 : "); put(Last_Zero_Count(s09),1); new_line;
    put("s10 : "); put(s10);
    put(", n10 : "); put(Last_Zero_Count(s10),1); new_line;
    put("s11 : "); put(s11);
    put(", n11 : "); put(Last_Zero_Count(s11),1); new_line;
    put("s12 : "); put(s12);
    put(", n12 : "); put(Last_Zero_Count(s12),1); new_line;
    put("s13 : "); put(s13);
    put(", n13 : "); put(Last_Zero_Count(s13),1); new_line;
    put("s14 : "); put(s14);
    put(", n14 : "); put(Last_Zero_Count(s14),1); new_line;
    put("s15 : "); put(s15);
    put(", n15 : "); put(Last_Zero_Count(s15),1); new_line;
  end Write_Subsums;

  procedure Write_Subsums
              ( s00,s01,s02,s03,s04,s05,s06,s07 : in double_float;
                s08,s09,s10,s11,s12,s13,s14,s15 : in double_float;
                s16,s17,s18,s19,s20,s21,s22,s23 : in double_float ) is

    use Bits_of_Doubles;

  begin
    put("s00 : "); put(s00);
    put(", n00 : "); put(Last_Zero_Count(s00),1); new_line;
    put("s01 : "); put(s01);
    put(", n01 : "); put(Last_Zero_Count(s01),1); new_line;
    put("s02 : "); put(s02);
    put(", n02 : "); put(Last_Zero_Count(s02),1); new_line;
    put("s03 : "); put(s03);
    put(", n03 : "); put(Last_Zero_Count(s03),1); new_line;
    put("s04 : "); put(s04);
    put(", n04 : "); put(Last_Zero_Count(s04),1); new_line;
    put("s05 : "); put(s05);
    put(", n05 : "); put(Last_Zero_Count(s05),1); new_line;
    put("s06 : "); put(s06);
    put(", n06 : "); put(Last_Zero_Count(s06),1); new_line;
    put("s07 : "); put(s07);
    put(", n07 : "); put(Last_Zero_Count(s07),1); new_line;
    put("s08 : "); put(s08);
    put(", n08 : "); put(Last_Zero_Count(s08),1); new_line;
    put("s09 : "); put(s09);
    put(", n09 : "); put(Last_Zero_Count(s09),1); new_line;
    put("s10 : "); put(s10);
    put(", n10 : "); put(Last_Zero_Count(s10),1); new_line;
    put("s11 : "); put(s11);
    put(", n11 : "); put(Last_Zero_Count(s11),1); new_line;
    put("s12 : "); put(s12);
    put(", n12 : "); put(Last_Zero_Count(s12),1); new_line;
    put("s13 : "); put(s13);
    put(", n13 : "); put(Last_Zero_Count(s13),1); new_line;
    put("s14 : "); put(s14);
    put(", n14 : "); put(Last_Zero_Count(s14),1); new_line;
    put("s15 : "); put(s15);
    put(", n15 : "); put(Last_Zero_Count(s15),1); new_line;
    put("s16 : "); put(s16);
    put(", n16 : "); put(Last_Zero_Count(s16),1); new_line;
    put("s17 : "); put(s17);
    put(", n17 : "); put(Last_Zero_Count(s17),1); new_line;
    put("s18 : "); put(s18);
    put(", n18 : "); put(Last_Zero_Count(s18),1); new_line;
    put("s19 : "); put(s19);
    put(", n19 : "); put(Last_Zero_Count(s19),1); new_line;
    put("s20 : "); put(s20);
    put(", n20 : "); put(Last_Zero_Count(s20),1); new_line;
    put("s21 : "); put(s21);
    put(", n21 : "); put(Last_Zero_Count(s21),1); new_line;
    put("s22 : "); put(s22);
    put(", n22 : "); put(Last_Zero_Count(s22),1); new_line;
    put("s23 : "); put(s23);
    put(", n23 : "); put(Last_Zero_Count(s23),1); new_line;
  end Write_Subsums;

  procedure Write_Subsums
              ( s00,s01,s02,s03,s04,s05,s06,s07 : in double_float;
                s08,s09,s10,s11,s12,s13,s14,s15 : in double_float;
                s16,s17,s18,s19,s20,s21,s22,s23 : in double_float;
                s24,s25,s26,s27,s28,s29,s30,s31 : in double_float ) is

    use Bits_of_Doubles;

  begin
    put("s00 : "); put(s00);
    put(", n00 : "); put(Last_Zero_Count(s00),1); new_line;
    put("s01 : "); put(s01);
    put(", n01 : "); put(Last_Zero_Count(s01),1); new_line;
    put("s02 : "); put(s02);
    put(", n02 : "); put(Last_Zero_Count(s02),1); new_line;
    put("s03 : "); put(s03);
    put(", n03 : "); put(Last_Zero_Count(s03),1); new_line;
    put("s04 : "); put(s04);
    put(", n04 : "); put(Last_Zero_Count(s04),1); new_line;
    put("s05 : "); put(s05);
    put(", n05 : "); put(Last_Zero_Count(s05),1); new_line;
    put("s06 : "); put(s06);
    put(", n06 : "); put(Last_Zero_Count(s06),1); new_line;
    put("s07 : "); put(s07);
    put(", n07 : "); put(Last_Zero_Count(s07),1); new_line;
    put("s08 : "); put(s08);
    put(", n08 : "); put(Last_Zero_Count(s08),1); new_line;
    put("s09 : "); put(s09);
    put(", n09 : "); put(Last_Zero_Count(s09),1); new_line;
    put("s10 : "); put(s10);
    put(", n10 : "); put(Last_Zero_Count(s10),1); new_line;
    put("s11 : "); put(s11);
    put(", n11 : "); put(Last_Zero_Count(s11),1); new_line;
    put("s12 : "); put(s12);
    put(", n12 : "); put(Last_Zero_Count(s12),1); new_line;
    put("s13 : "); put(s13);
    put(", n13 : "); put(Last_Zero_Count(s13),1); new_line;
    put("s14 : "); put(s14);
    put(", n14 : "); put(Last_Zero_Count(s14),1); new_line;
    put("s15 : "); put(s15);
    put(", n15 : "); put(Last_Zero_Count(s15),1); new_line;
    put("s16 : "); put(s16);
    put(", n16 : "); put(Last_Zero_Count(s16),1); new_line;
    put("s17 : "); put(s17);
    put(", n17 : "); put(Last_Zero_Count(s17),1); new_line;
    put("s18 : "); put(s18);
    put(", n18 : "); put(Last_Zero_Count(s18),1); new_line;
    put("s19 : "); put(s19);
    put(", n19 : "); put(Last_Zero_Count(s19),1); new_line;
    put("s20 : "); put(s20);
    put(", n20 : "); put(Last_Zero_Count(s20),1); new_line;
    put("s21 : "); put(s21);
    put(", n21 : "); put(Last_Zero_Count(s21),1); new_line;
    put("s22 : "); put(s22);
    put(", n22 : "); put(Last_Zero_Count(s22),1); new_line;
    put("s23 : "); put(s23);
    put(", n23 : "); put(Last_Zero_Count(s23),1); new_line;
    put("s24 : "); put(s24);
    put(", n24 : "); put(Last_Zero_Count(s24),1); new_line;
    put("s25 : "); put(s25);
    put(", n25 : "); put(Last_Zero_Count(s25),1); new_line;
    put("s26 : "); put(s26);
    put(", n26 : "); put(Last_Zero_Count(s26),1); new_line;
    put("s27 : "); put(s27);
    put(", n27 : "); put(Last_Zero_Count(s27),1); new_line;
    put("s28 : "); put(s28);
    put(", n28 : "); put(Last_Zero_Count(s28),1); new_line;
    put("s29 : "); put(s29);
    put(", n29 : "); put(Last_Zero_Count(s29),1); new_line;
    put("s30 : "); put(s30);
    put(", n30 : "); put(Last_Zero_Count(s30),1); new_line;
    put("s31 : "); put(s31);
    put(", n31 : "); put(Last_Zero_Count(s31),1); new_line;
  end Write_Subsums;

  function to_octo_double
              ( s00,s01,s02,s03,s04,s05,s06,s07 : double_float;
                verbose : boolean := true ) return Octo_double is

    res : octo_double;

  begin
    if verbose then
      write_subsums(s00,s01,s02,s03,s04,s05,s06,s07);
    end if;
    res := create(s07);
    res := res + create(s06);
    res := res + create(s05);
    res := res + create(s04);
    res := res + create(s03);
    res := res + create(s02);
    res := res + create(s01);
    res := res + create(s00);
    return res;
  end to_octo_double;

  function to_octo_double
              ( s00,s01,s02,s03,s04,s05,s06,s07 : double_float;
                s08,s09,s10,s11,s12,s13,s14,s15 : double_float;
                verbose : boolean := true ) return Octo_double is

    res : octo_double;

  begin
    if verbose then
      write_subsums(s00,s01,s02,s03,s04,s05,s06,s07,
                    s08,s09,s10,s11,s12,s13,s14,s15);
    end if;
    res := create(s15);
    res := res + create(s14);
    res := res + create(s13);
    res := res + create(s12);
    res := res + create(s11);
    res := res + create(s10);
    res := res + create(s09);
    res := res + create(s08);
    res := res + create(s07);
    res := res + create(s06);
    res := res + create(s05);
    res := res + create(s04);
    res := res + create(s03);
    res := res + create(s02);
    res := res + create(s01);
    res := res + create(s00);
    return res;
  end to_octo_double;

  function to_octo_double
              ( s00,s01,s02,s03,s04,s05,s06,s07 : double_float;
                s08,s09,s10,s11,s12,s13,s14,s15 : double_float;
                s16,s17,s18,s19,s20,s21,s22,s23 : double_float;
                verbose : boolean := true ) return Octo_double is

    res : octo_double;

  begin
    if verbose then
      write_subsums(s00,s01,s02,s03,s04,s05,s06,s07,
                    s08,s09,s10,s11,s12,s13,s14,s15,
                    s16,s17,s18,s19,s20,s21,s22,s23);
    end if;
    res := create(s23);
    res := res + create(s22);
    res := res + create(s21);
    res := res + create(s20);
    res := res + create(s19);
    res := res + create(s18);
    res := res + create(s17);
    res := res + create(s16);
    res := res + create(s15);
    res := res + create(s14);
    res := res + create(s13);
    res := res + create(s12);
    res := res + create(s11);
    res := res + create(s10);
    res := res + create(s09);
    res := res + create(s08);
    res := res + create(s07);
    res := res + create(s06);
    res := res + create(s05);
    res := res + create(s04);
    res := res + create(s03);
    res := res + create(s02);
    res := res + create(s01);
    res := res + create(s00);
    return res;
  end to_octo_double;

  function to_octo_double
              ( s00,s01,s02,s03,s04,s05,s06,s07 : double_float;
                s08,s09,s10,s11,s12,s13,s14,s15 : double_float;
                s16,s17,s18,s19,s20,s21,s22,s23 : double_float;
                s24,s25,s26,s27,s28,s29,s30,s31 : double_float;
                verbose : boolean := true ) return Octo_double is

    res : octo_double;

  begin
    if verbose then
      write_subsums(s00,s01,s02,s03,s04,s05,s06,s07,
                    s08,s09,s10,s11,s12,s13,s14,s15,
                    s16,s17,s18,s19,s20,s21,s22,s23,
                    s24,s25,s26,s27,s28,s29,s30,s31);
    end if;
    res := create(s31);
    res := res + create(s30);
    res := res + create(s29);
    res := res + create(s28);
    res := res + create(s27);
    res := res + create(s26);
    res := res + create(s25);
    res := res + create(s24);
    res := res + create(s23);
    res := res + create(s22);
    res := res + create(s21);
    res := res + create(s20);
    res := res + create(s19);
    res := res + create(s18);
    res := res + create(s17);
    res := res + create(s16);
    res := res + create(s15);
    res := res + create(s14);
    res := res + create(s13);
    res := res + create(s12);
    res := res + create(s11);
    res := res + create(s10);
    res := res + create(s09);
    res := res + create(s08);
    res := res + create(s07);
    res := res + create(s06);
    res := res + create(s05);
    res := res + create(s04);
    res := res + create(s03);
    res := res + create(s02);
    res := res + create(s01);
    res := res + create(s00);
    return res;
  end to_octo_double;

-- SIGN AWARE WRAPPERS :

  function Product ( x,y : Octo_Double_Vectors.Vector;
                     verbose : boolean := true ) return octo_double is

    res : octo_double := create(0.0);
    xb : constant Octo_Double_Vectors.Vector(x'range)
       := Sign_Balance(x,verbose);
    yb : constant Octo_Double_Vectors.Vector(y'range)
       := Sign_Balance(y,verbose);
    xs00,xs01,xs02,xs03 : Standard_Floating_Vectors.Vector(x'range);
    xs04,xs05,xs06,xs07 : Standard_Floating_Vectors.Vector(x'range);
    xs08,xs09,xs10,xs11 : Standard_Floating_Vectors.Vector(x'range);
    xs12,xs13,xs14,xs15 : Standard_Floating_Vectors.Vector(x'range);
    xs16,xs17,xs18,xs19 : Standard_Floating_Vectors.Vector(x'range);
    xs20,xs21,xs22,xs23 : Standard_Floating_Vectors.Vector(x'range);
    xs24,xs25,xs26,xs27 : Standard_Floating_Vectors.Vector(x'range);
    xs28,xs29,xs30,xs31 : Standard_Floating_Vectors.Vector(x'range);
    ys00,ys01,ys02,ys03 : Standard_Floating_Vectors.Vector(y'range);
    ys04,ys05,ys06,ys07 : Standard_Floating_Vectors.Vector(y'range);
    ys08,ys09,ys10,ys11 : Standard_Floating_Vectors.Vector(y'range);
    ys12,ys13,ys14,ys15 : Standard_Floating_Vectors.Vector(y'range);
    ys16,ys17,ys18,ys19 : Standard_Floating_Vectors.Vector(y'range);
    ys20,ys21,ys22,ys23 : Standard_Floating_Vectors.Vector(y'range);
    ys24,ys25,ys26,ys27 : Standard_Floating_Vectors.Vector(y'range);
    ys28,ys29,ys30,ys31 : Standard_Floating_Vectors.Vector(y'range);
    xd00,xd01,xd02,xd03 : Standard_Floating_Vectors.Vector(x'range);
    xd04,xd05,xd06,xd07 : Standard_Floating_Vectors.Vector(x'range);
    xd08,xd09,xd10,xd11 : Standard_Floating_Vectors.Vector(x'range);
    xd12,xd13,xd14,xd15 : Standard_Floating_Vectors.Vector(x'range);
    xd16,xd17,xd18,xd19 : Standard_Floating_Vectors.Vector(x'range);
    xd20,xd21,xd22,xd23 : Standard_Floating_Vectors.Vector(x'range);
    xd24,xd25,xd26,xd27 : Standard_Floating_Vectors.Vector(x'range);
    xd28,xd29,xd30,xd31 : Standard_Floating_Vectors.Vector(x'range);
    yd00,yd01,yd02,yd03 : Standard_Floating_Vectors.Vector(y'range);
    yd04,yd05,yd06,yd07 : Standard_Floating_Vectors.Vector(y'range);
    yd08,yd09,yd10,yd11 : Standard_Floating_Vectors.Vector(y'range);
    yd12,yd13,yd14,yd15 : Standard_Floating_Vectors.Vector(y'range);
    yd16,yd17,yd18,yd19 : Standard_Floating_Vectors.Vector(y'range);
    yd20,yd21,yd22,yd23 : Standard_Floating_Vectors.Vector(y'range);
    yd24,yd25,yd26,yd27 : Standard_Floating_Vectors.Vector(y'range);
    yd28,yd29,yd30,yd31 : Standard_Floating_Vectors.Vector(y'range);
    ns,nd : integer32;
    s00,s01,s02,s03,s04,s05,s06,s07 : double_float;
    s08a,s09a,s10a,s11a,s12a,s13a,s14a,s15a : double_float;
    s08b,s09b,s10b,s11b,s12b,s13b,s14b,s15b : double_float;
    s16a,s17a,s18a,s19a,s20a,s21a,s22a,s23a : double_float;
    s16b,s17b,s18b,s19b,s20b,s21b,s22b,s23b : double_float;
    s16c,s17c,s18c,s19c,s20c,s21c,s22c,s23c : double_float;
    s24a,s25a,s26a,s27a,s28a,s29a,s30a,s31a : double_float;
    s24b,s25b,s26b,s27b,s28b,s29b,s30b,s31b : double_float;
    s24c,s25c,s26c,s27c,s28c,s29c,s30c,s31c : double_float;
    s24d,s25d,s26d,s27d,s28d,s29d,s30d,s31d : double_float;

  begin
    Signed_Quarter(xb,yb,xs00,xs01,xs02,xs03,xs04,xs05,xs06,xs07,
                         xs08,xs09,xs10,xs11,xs12,xs13,xs14,xs15,
                         xs16,xs17,xs18,xs19,xs20,xs21,xs22,xs23,
                         xs24,xs25,xs26,xs27,xs28,xs29,xs30,xs31,
                         ys00,ys01,ys02,ys03,ys04,ys05,ys06,ys07,
                         ys08,ys09,ys10,ys11,ys12,ys13,ys14,ys15,
                         ys16,ys17,ys18,ys19,ys20,ys21,ys22,ys23,
                         ys24,ys25,ys26,ys27,ys28,ys29,ys30,ys31,
                         xd00,xd01,xd02,xd03,xd04,xd05,xd06,xd07,
                         xd08,xd09,xd10,xd11,xd12,xd13,xd14,xd15,
                         xd16,xd17,xd18,xd19,xd20,xd21,xd22,xd23,
                         xd24,xd25,xd26,xd27,xd28,xd29,xd30,xd31,
                         yd00,yd01,yd02,yd03,yd04,yd05,yd06,yd07,
                         yd08,yd09,yd10,yd11,yd12,yd13,yd14,yd15,
                         yd16,yd17,yd18,yd19,yd20,yd21,yd22,yd23,
                         yd24,yd25,yd26,yd27,yd28,yd29,yd30,yd31,ns,nd);
    if verbose then
      put("#s : "); put(ns,1); 
      put(", #d : "); put(nd,1); new_line;
    end if;
    s00 := 0.0; s01 := 0.0; s02 := 0.0; s03 := 0.0;
    s04 := 0.0; s05 := 0.0; s06 := 0.0; s07 := 0.0;
    s08a := 0.0; s09a := 0.0; s10a := 0.0; s11a := 0.0;
    s12a := 0.0; s13a := 0.0; s14a := 0.0; s15a := 0.0;
    s08b := 0.0; s09b := 0.0; s10b := 0.0; s11b := 0.0;
    s12b := 0.0; s13b := 0.0; s14b := 0.0; s15b := 0.0;
    s16a := 0.0; s17a := 0.0; s18a := 0.0; s19a := 0.0;
    s20a := 0.0; s21a := 0.0; s22a := 0.0; s23a := 0.0;
    s16b := 0.0; s17b := 0.0; s18b := 0.0; s19b := 0.0;
    s20b := 0.0; s21b := 0.0; s22b := 0.0; s23b := 0.0;
    s16c := 0.0; s17c := 0.0; s18c := 0.0; s19c := 0.0;
    s20c := 0.0; s21c := 0.0; s22c := 0.0; s23c := 0.0;
    s24a := 0.0; s25a := 0.0; s26a := 0.0; s27a := 0.0;
    s28a := 0.0; s29a := 0.0; s30a := 0.0; s31a := 0.0;
    s24b := 0.0; s25b := 0.0; s26b := 0.0; s27b := 0.0;
    s28b := 0.0; s29b := 0.0; s30b := 0.0; s31b := 0.0;
    s24c := 0.0; s25c := 0.0; s26c := 0.0; s27c := 0.0;
    s28c := 0.0; s29c := 0.0; s30c := 0.0; s31c := 0.0;
    s24d := 0.0; s25d := 0.0; s26d := 0.0; s27d := 0.0;
    s28d := 0.0; s29d := 0.0; s30d := 0.0; s31d := 0.0;
    for i in 1..ns loop
      s00 := s00 + xs00(i)*ys00(i);
      s01 := s01 + xs00(i)*ys01(i) + xs01(i)*ys00(i);
      s02 := s02 + xs00(i)*ys02(i) + xs01(i)*ys01(i) + xs02(i)*ys00(i);
      s03 := s03 + xs00(i)*ys03(i) + xs01(i)*ys02(i) + xs02(i)*ys01(i)
                 + xs03(i)*ys00(i);
      s04 := s04 + xs00(i)*ys04(i) + xs01(i)*ys03(i) + xs02(i)*ys02(i)
                 + xs03(i)*ys01(i) + xs04(i)*ys00(i);
      s05 := s05 + xs00(i)*ys05(i) + xs01(i)*ys04(i) + xs02(i)*ys03(i)
                 + xs03(i)*ys02(i) + xs04(i)*ys01(i) + xs05(i)*ys00(i);
      s06 := s06 + xs00(i)*ys06(i) + xs01(i)*ys05(i) + xs02(i)*ys04(i)
                 + xs03(i)*ys03(i) + xs04(i)*ys02(i) + xs05(i)*ys01(i)
                 + xs06(i)*ys00(i);
      s07 := s07 + xs00(i)*ys07(i) + xs01(i)*ys06(i) + xs02(i)*ys05(i)
                 + xs03(i)*ys04(i) + xs04(i)*ys03(i) + xs05(i)*ys02(i)
                 + xs06(i)*ys01(i) + xs07(i)*ys00(i);
      s08a := s08a + xs00(i)*ys08(i) + xs01(i)*ys07(i) + xs02(i)*ys06(i)
                   + xs03(i)*ys05(i) + xs04(i)*ys04(i);
      s08b := s08b + xs05(i)*ys03(i) + xs06(i)*ys02(i) + xs07(i)*ys01(i)
                   + xs08(i)*ys00(i);
      s09a := s09a + xs00(i)*ys09(i) + xs01(i)*ys08(i) + xs02(i)*ys07(i)
                   + xs03(i)*ys06(i) + xs04(i)*ys05(i);
      s09b := s09b + xs05(i)*ys04(i) + xs06(i)*ys03(i) + xs07(i)*ys02(i)
                   + xs08(i)*ys01(i) + xs09(i)*ys00(i);
      s10a := s10a + xs00(i)*ys10(i) + xs01(i)*ys09(i) + xs02(i)*ys08(i)
                   + xs03(i)*ys07(i) + xs04(i)*ys06(i) + xs05(i)*ys05(i);
      s10b := s10b + xs06(i)*ys04(i) + xs07(i)*ys03(i) + xs08(i)*ys02(i)
                   + xs09(i)*ys01(i) + xs10(i)*ys00(i);
      s11a := s11a + xs00(i)*ys11(i) + xs01(i)*ys10(i) + xs02(i)*ys09(i)
                   + xs03(i)*ys08(i) + xs04(i)*ys07(i) + xs05(i)*ys06(i);
      s11b := s11b + xs06(i)*ys05(i) + xs07(i)*ys04(i) + xs08(i)*ys03(i)
                   + xs09(i)*ys02(i) + xs10(i)*ys01(i) + xs11(i)*ys00(i);
      s12a := s12a + xs00(i)*ys12(i) + xs01(i)*ys11(i) + xs02(i)*ys10(i)
                   + xs03(i)*ys09(i) + xs04(i)*ys08(i) + xs05(i)*ys07(i)
                   + xs06(i)*ys06(i);
      s12b := s12b + xs07(i)*ys05(i) + xs08(i)*ys04(i) + xs09(i)*ys03(i)
                   + xs10(i)*ys02(i) + xs11(i)*ys01(i) + xs12(i)*ys00(i);
      s13a := s13a + xs00(i)*ys13(i) + xs01(i)*ys12(i) + xs02(i)*ys11(i)
                   + xs03(i)*ys10(i) + xs04(i)*ys09(i) + xs05(i)*ys08(i)
                   + xs06(i)*ys07(i);
      s13b := s13b + xs07(i)*ys06(i) + xs08(i)*ys05(i) + xs09(i)*ys04(i)
                   + xs10(i)*ys03(i) + xs11(i)*ys02(i) + xs12(i)*ys01(i)
                   + xs13(i)*ys00(i);
      s14a := s14a + xs00(i)*ys14(i) + xs01(i)*ys13(i) + xs02(i)*ys12(i)
                   + xs03(i)*ys11(i) + xs04(i)*ys10(i) + xs05(i)*ys09(i)
                   + xs06(i)*ys08(i) + xs07(i)*ys07(i);
      s14b := s14b + xs08(i)*ys06(i) + xs09(i)*ys05(i) + xs10(i)*ys04(i)
                   + xs11(i)*ys03(i) + xs12(i)*ys02(i) + xs13(i)*ys01(i)
                   + xs14(i)*ys00(i);
      s15a := s15a + xs00(i)*ys15(i) + xs01(i)*ys14(i) + xs02(i)*ys13(i)
                   + xs03(i)*ys12(i) + xs04(i)*ys11(i) + xs05(i)*ys10(i)
                   + xs06(i)*ys09(i) + xs07(i)*ys08(i);
      s15b := s15b + xs08(i)*ys07(i) + xs09(i)*ys06(i) + xs10(i)*ys05(i)
                   + xs11(i)*ys04(i) + xs12(i)*ys03(i) + xs13(i)*ys02(i)
                   + xs14(i)*ys01(i) + xs15(i)*ys00(i);
      s16a := s16a + xs00(i)*ys16(i) + xs01(i)*ys15(i) + xs02(i)*ys14(i)
                   + xs03(i)*ys13(i) + xs04(i)*ys12(i) + xs05(i)*ys11(i);
      s16b := s16b + xs06(i)*ys10(i) + xs07(i)*ys09(i) + xs08(i)*ys08(i)
                   + xs09(i)*ys07(i) + xs10(i)*ys06(i) + xs11(i)*ys05(i);
      s16c := s16c + xs12(i)*ys04(i) + xs13(i)*ys03(i) + xs14(i)*ys02(i)
                   + xs15(i)*ys01(i) + xs16(i)*ys00(i);
      s17a := s17a + xs00(i)*ys17(i) + xs01(i)*ys16(i) + xs02(i)*ys15(i)
                   + xs03(i)*ys14(i) + xs04(i)*ys13(i) + xs05(i)*ys12(i);
      s17b := s17b + xs06(i)*ys11(i) + xs07(i)*ys10(i) + xs08(i)*ys09(i)
                   + xs09(i)*ys08(i) + xs10(i)*ys07(i) + xs11(i)*ys06(i);
      s17c := s17c + xs12(i)*ys05(i) + xs13(i)*ys04(i) + xs14(i)*ys03(i)
                   + xs15(i)*ys02(i) + xs16(i)*ys01(i) + xs17(i)*ys00(i);
      s18a := s18a + xs00(i)*ys18(i) + xs01(i)*ys17(i) + xs02(i)*ys16(i)
                   + xs03(i)*ys15(i) + xs04(i)*ys14(i) + xs05(i)*ys13(i)
                   + xs06(i)*ys12(i);
      s18b := s18b + xs07(i)*ys11(i) + xs08(i)*ys10(i) + xs09(i)*ys09(i)
                   + xs10(i)*ys08(i) + xs11(i)*ys07(i) + xs12(i)*ys06(i);
      s18c := s18c + xs13(i)*ys05(i) + xs14(i)*ys04(i) + xs15(i)*ys03(i)
                   + xs16(i)*ys02(i) + xs17(i)*ys01(i) + xs18(i)*ys00(i);
      s19a := s19a + xs00(i)*ys19(i) + xs01(i)*ys18(i) + xs02(i)*ys17(i)
                   + xs03(i)*ys16(i) + xs04(i)*ys15(i) + xs05(i)*ys14(i)
                   + xs06(i)*ys13(i);
      s19b := s19b + xs07(i)*ys12(i) + xs08(i)*ys11(i) + xs09(i)*ys10(i)
                   + xs10(i)*ys09(i) + xs11(i)*ys08(i) + xs12(i)*ys07(i)
                   + xs13(i)*ys06(i);
      s19c := s19c + xs14(i)*ys05(i) + xs15(i)*ys04(i) + xs16(i)*ys03(i)
                   + xs17(i)*ys02(i) + xs18(i)*ys01(i) + xs19(i)*ys00(i);
      s20a := s20a + xs00(i)*ys20(i) + xs01(i)*ys19(i) + xs02(i)*ys18(i)
                   + xs03(i)*ys17(i) + xs04(i)*ys16(i) + xs05(i)*ys15(i)
                   + xs06(i)*ys14(i);
      s20b := s20b + xs07(i)*ys13(i) + xs08(i)*ys12(i) + xs09(i)*ys11(i)
                   + xs10(i)*ys10(i) + xs11(i)*ys09(i) + xs12(i)*ys08(i)
                   + xs13(i)*ys07(i);
      s20c := s20c + xs14(i)*ys06(i) + xs15(i)*ys05(i) + xs16(i)*ys04(i)
                   + xs17(i)*ys03(i) + xs18(i)*ys02(i) + xs19(i)*ys01(i)
                   + xs20(i)*ys00(i);
      s21a := s21a + xs00(i)*ys21(i) + xs01(i)*ys20(i) + xs02(i)*ys19(i)
                   + xs03(i)*ys18(i) + xs04(i)*ys17(i) + xs05(i)*ys16(i)
                   + xs06(i)*ys15(i) + xs07(i)*ys14(i);
      s21b := s21b + xs08(i)*ys13(i) + xs09(i)*ys12(i) + xs10(i)*ys11(i)
                   + xs11(i)*ys10(i) + xs12(i)*ys09(i) + xs13(i)*ys08(i)
                   + xs14(i)*ys07(i);
      s21c := s21c + xs15(i)*ys06(i) + xs16(i)*ys05(i) + xs17(i)*ys04(i)
                   + xs18(i)*ys03(i) + xs19(i)*ys02(i) + xs20(i)*ys01(i)
                   + xs21(i)*ys00(i);
      s22a := s22a + xs00(i)*ys22(i) + xs01(i)*ys21(i) + xs02(i)*ys20(i)
                   + xs03(i)*ys19(i) + xs04(i)*ys18(i) + xs05(i)*ys17(i)
                   + xs06(i)*ys16(i) + xs07(i)*ys15(i);
      s22b := s22b + xs08(i)*ys14(i) + xs09(i)*ys13(i) + xs10(i)*ys12(i)
                   + xs11(i)*ys11(i) + xs12(i)*ys10(i) + xs13(i)*ys09(i)
                   + xs14(i)*ys08(i) + xs15(i)*ys07(i);
      s22c := s22c + xs16(i)*ys06(i) + xs17(i)*ys05(i) + xs18(i)*ys04(i)
                   + xs19(i)*ys03(i) + xs20(i)*ys02(i) + xs21(i)*ys01(i)
                   + xs22(i)*ys00(i);
      s23a := s23a + xs00(i)*ys23(i) + xs01(i)*ys22(i) + xs02(i)*ys21(i)
                   + xs03(i)*ys20(i) + xs04(i)*ys19(i) + xs05(i)*ys18(i)
                   + xs06(i)*ys17(i) + xs07(i)*ys16(i);
      s23b := s23b + xs08(i)*ys15(i) + xs09(i)*ys14(i) + xs10(i)*ys13(i)
                   + xs11(i)*ys12(i) + xs12(i)*ys11(i) + xs13(i)*ys10(i)
                   + xs14(i)*ys09(i) + xs15(i)*ys08(i);
      s23c := s23c + xs16(i)*ys07(i) + xs17(i)*ys06(i) + xs18(i)*ys05(i)
                   + xs19(i)*ys04(i) + xs20(i)*ys03(i) + xs21(i)*ys02(i)
                   + xs22(i)*ys01(i) + xs23(i)*ys00(i);
      s24a := s24a + xs00(i)*ys24(i) + xs01(i)*ys23(i) + xs02(i)*ys22(i)
                   + xs03(i)*ys21(i) + xs04(i)*ys20(i) + xs05(i)*ys19(i)
                   + xs06(i)*ys18(i);
      s24b := s24b + xs07(i)*ys17(i) + xs08(i)*ys16(i) + xs09(i)*ys15(i)
                   + xs10(i)*ys14(i) + xs11(i)*ys13(i) + xs12(i)*ys12(i);
      s24c := s24c + xs13(i)*ys11(i) + xs14(i)*ys10(i) + xs15(i)*ys09(i)
                   + xs16(i)*ys08(i) + xs17(i)*ys07(i) + xs18(i)*ys06(i);
      s24d := s24d + xs19(i)*ys05(i) + xs20(i)*ys04(i) + xs21(i)*ys03(i)
                   + xs22(i)*ys02(i) + xs23(i)*ys01(i) + xs24(i)*ys00(i);
      s25a := s25a + xs00(i)*ys25(i) + xs01(i)*ys24(i) + xs02(i)*ys23(i)
                   + xs03(i)*ys22(i) + xs04(i)*ys21(i) + xs05(i)*ys20(i)
                   + xs06(i)*ys19(i);
      s25b := s25b + xs07(i)*ys18(i) + xs08(i)*ys17(i) + xs09(i)*ys16(i)
                   + xs10(i)*ys15(i) + xs11(i)*ys14(i) + xs12(i)*ys13(i)
                   + xs13(i)*ys12(i);
      s25c := s25c + xs14(i)*ys11(i) + xs15(i)*ys10(i) + xs16(i)*ys09(i)
                   + xs17(i)*ys08(i) + xs18(i)*ys07(i) + xs19(i)*ys06(i);
      s25d := s25d + xs20(i)*ys05(i) + xs21(i)*ys04(i) + xs22(i)*ys03(i)
                   + xs23(i)*ys02(i) + xs24(i)*ys01(i) + xs25(i)*ys00(i);
      s26a := s26a + xs00(i)*ys26(i) + xs01(i)*ys25(i) + xs02(i)*ys24(i)
                   + xs03(i)*ys23(i) + xs04(i)*ys22(i) + xs05(i)*ys21(i)
                   + xs06(i)*ys20(i);
      s26b := s26b + xs07(i)*ys19(i) + xs08(i)*ys18(i) + xs09(i)*ys17(i)
                   + xs10(i)*ys16(i) + xs11(i)*ys15(i) + xs12(i)*ys14(i)
                   + xs13(i)*ys13(i);
      s26c := s26c + xs14(i)*ys12(i) + xs15(i)*ys11(i) + xs16(i)*ys10(i)
                   + xs17(i)*ys09(i) + xs18(i)*ys08(i) + xs19(i)*ys07(i)
                   + xs20(i)*ys06(i);
      s26d := s26d + xs21(i)*ys05(i) + xs22(i)*ys04(i) + xs23(i)*ys03(i)
                   + xs24(i)*ys02(i) + xs25(i)*ys01(i) + xs26(i)*ys00(i);
      s27a := s27a + xs00(i)*ys27(i) + xs01(i)*ys26(i) + xs02(i)*ys25(i)
                   + xs03(i)*ys24(i) + xs04(i)*ys23(i) + xs05(i)*ys22(i)
                   + xs06(i)*ys21(i);
      s27b := s27b + xs07(i)*ys20(i) + xs08(i)*ys19(i) + xs09(i)*ys18(i)
                   + xs10(i)*ys17(i) + xs11(i)*ys16(i) + xs12(i)*ys15(i)
                   + xs13(i)*ys14(i);
      s27c := s27c + xs14(i)*ys13(i) + xs15(i)*ys12(i) + xs16(i)*ys11(i)
                   + xs17(i)*ys10(i) + xs18(i)*ys09(i) + xs19(i)*ys08(i)
                   + xs20(i)*ys07(i);
      s27d := s27d + xs21(i)*ys06(i) + xs22(i)*ys05(i) + xs23(i)*ys04(i)
                   + xs24(i)*ys03(i) + xs25(i)*ys02(i) + xs26(i)*ys01(i)
                   + xs27(i)*ys00(i);
      s28a := s28a + xs00(i)*ys28(i) + xs01(i)*ys27(i) + xs02(i)*ys26(i)
                   + xs03(i)*ys25(i) + xs04(i)*ys24(i) + xs05(i)*ys23(i)
                   + xs06(i)*ys22(i) + xs07(i)*ys21(i);
      s28b := s28b + xs08(i)*ys20(i) + xs09(i)*ys19(i) + xs10(i)*ys18(i)
                   + xs11(i)*ys17(i) + xs12(i)*ys16(i) + xs13(i)*ys15(i)
                   + xs14(i)*ys14(i);
      s28c := s28c + xs15(i)*ys13(i) + xs16(i)*ys12(i) + xs17(i)*ys11(i)
                   + xs18(i)*ys10(i) + xs19(i)*ys09(i) + xs20(i)*ys08(i)
                   + xs21(i)*ys07(i);
      s28d := s28d + xs22(i)*ys06(i) + xs23(i)*ys05(i) + xs24(i)*ys04(i)
                   + xs25(i)*ys03(i) + xs26(i)*ys02(i) + xs27(i)*ys01(i)
                   + xs28(i)*ys00(i);
      s29a := s29a + xs00(i)*ys29(i) + xs01(i)*ys28(i) + xs02(i)*ys27(i)
                   + xs03(i)*ys26(i) + xs04(i)*ys25(i) + xs05(i)*ys24(i)
                   + xs06(i)*ys23(i) + xs07(i)*ys22(i);
      s29b := s29b + xs08(i)*ys21(i) + xs09(i)*ys20(i) + xs10(i)*ys19(i)
                   + xs11(i)*ys18(i) + xs12(i)*ys17(i) + xs13(i)*ys16(i)
                   + xs14(i)*ys15(i) + xs15(i)*ys14(i);
      s29c := s29c + xs16(i)*ys13(i) + xs17(i)*ys12(i) + xs18(i)*ys11(i)
                   + xs19(i)*ys10(i) + xs20(i)*ys09(i) + xs21(i)*ys08(i)
                   + xs22(i)*ys07(i);
      s29d := s29d + xs23(i)*ys06(i) + xs24(i)*ys05(i) + xs25(i)*ys04(i)
                   + xs26(i)*ys03(i) + xs27(i)*ys02(i) + xs28(i)*ys01(i)
                   + xs29(i)*ys00(i);
      s30a := s30a + xs00(i)*ys30(i) + xs01(i)*ys29(i) + xs02(i)*ys28(i)
                   + xs03(i)*ys27(i) + xs04(i)*ys26(i) + xs05(i)*ys25(i)
                   + xs06(i)*ys24(i) + xs07(i)*ys23(i);
      s30b := s30b + xs08(i)*ys22(i) + xs09(i)*ys21(i) + xs10(i)*ys20(i)
                   + xs11(i)*ys19(i) + xs12(i)*ys18(i) + xs13(i)*ys17(i)
                   + xs14(i)*ys16(i) + xs15(i)*ys15(i);
      s30c := s30c + xs16(i)*ys14(i) + xs17(i)*ys13(i) + xs18(i)*ys12(i)
                   + xs19(i)*ys11(i) + xs20(i)*ys10(i) + xs21(i)*ys09(i)
                   + xs22(i)*ys08(i) + xs23(i)*ys07(i);
      s30d := s30d + xs24(i)*ys06(i) + xs25(i)*ys05(i) + xs26(i)*ys04(i)
                   + xs27(i)*ys03(i) + xs28(i)*ys02(i) + xs29(i)*ys01(i)
                   + xs30(i)*ys00(i);
      s31a := s31a + xs00(i)*ys31(i) + xs01(i)*ys30(i) + xs02(i)*ys29(i)
                   + xs03(i)*ys28(i) + xs04(i)*ys27(i) + xs05(i)*ys26(i)
                   + xs06(i)*ys25(i) + xs07(i)*ys24(i);
      s31b := s31b + xs08(i)*ys23(i) + xs09(i)*ys22(i) + xs10(i)*ys21(i)
                   + xs11(i)*ys20(i) + xs12(i)*ys19(i) + xs13(i)*ys18(i)
                   + xs14(i)*ys17(i) + xs15(i)*ys16(i);
      s31c := s31c + xs16(i)*ys15(i) + xs17(i)*ys14(i) + xs18(i)*ys13(i)
                   + xs19(i)*ys12(i) + xs20(i)*ys11(i) + xs21(i)*ys10(i)
                   + xs22(i)*ys09(i) + xs23(i)*ys08(i);
      s31d := s31d + xs24(i)*ys07(i) + xs25(i)*ys06(i) + xs26(i)*ys05(i)
                   + xs27(i)*ys04(i) + xs28(i)*ys03(i) + xs29(i)*ys02(i)
                   + xs30(i)*ys01(i) + xs31(i)*ys00(i);
    end loop;
    if ns > 0 then
      res := to_octo_double(s00,s01,s02,s03,s04,s05,s06,s07,
                            s08a,s09a,s10a,s11a,s12a,s13a,s14a,s15a,
                            s16a,s17a,s18a,s19a,s20a,s21a,s22a,s23a,
                            s24a,s25a,s26a,s27a,s28a,s29a,s30a,s31a,
                            verbose=>true)
           + to_octo_double(s08b,s09b,s10b,s11b,s12b,s13b,s14b,s15b,
                            s16b,s17b,s18b,s19b,s20b,s21b,s22b,s23b,
                            s24b,s25b,s26b,s27b,s28b,s29b,s30b,s31b,
                            verbose=>true)
           + to_octo_double(s16c,s17c,s18c,s19c,s20c,s21c,s22c,s23c,
                            s24c,s25c,s26c,s27c,s28c,s29c,s30c,s31c,
                            verbose=>true)
           + to_octo_double(s24d,s25d,s26d,s27d,s28d,s29d,s30d,s31d,
                            verbose=>true);
    end if;
    s00 := 0.0; s01 := 0.0; s02 := 0.0; s03 := 0.0;
    s04 := 0.0; s05 := 0.0; s06 := 0.0; s07 := 0.0;
    s08a := 0.0; s09a := 0.0; s10a := 0.0; s11a := 0.0;
    s12a := 0.0; s13a := 0.0; s14a := 0.0; s15a := 0.0;
    s08b := 0.0; s09b := 0.0; s10b := 0.0; s11b := 0.0;
    s12b := 0.0; s13b := 0.0; s14b := 0.0; s15b := 0.0;
    s16a := 0.0; s17a := 0.0; s18a := 0.0; s19a := 0.0;
    s20a := 0.0; s21a := 0.0; s22a := 0.0; s23a := 0.0;
    s16b := 0.0; s17b := 0.0; s18b := 0.0; s19b := 0.0;
    s20b := 0.0; s21b := 0.0; s22b := 0.0; s23b := 0.0;
    s16c := 0.0; s17c := 0.0; s18c := 0.0; s19c := 0.0;
    s20c := 0.0; s21c := 0.0; s22c := 0.0; s23c := 0.0;
    s24a := 0.0; s25a := 0.0; s26a := 0.0; s27a := 0.0;
    s28a := 0.0; s29a := 0.0; s30a := 0.0; s31a := 0.0;
    s24b := 0.0; s25b := 0.0; s26b := 0.0; s27b := 0.0;
    s28b := 0.0; s29b := 0.0; s30b := 0.0; s31b := 0.0;
    s24c := 0.0; s25c := 0.0; s26c := 0.0; s27c := 0.0;
    s28c := 0.0; s29c := 0.0; s30c := 0.0; s31c := 0.0;
    s24d := 0.0; s25d := 0.0; s26d := 0.0; s27d := 0.0;
    s28d := 0.0; s29d := 0.0; s30d := 0.0; s31d := 0.0;
    for i in 1..nd loop
      s00 := s00 + xd00(i)*yd00(i);
      s01 := s01 + xd00(i)*yd01(i) + xd01(i)*yd00(i);
      s02 := s02 + xd00(i)*yd02(i) + xd01(i)*yd01(i) + xd02(i)*yd00(i);
      s03 := s03 + xd00(i)*yd03(i) + xd01(i)*yd02(i) + xd02(i)*yd01(i)
                 + xd03(i)*yd00(i);
      s04 := s04 + xd00(i)*yd04(i) + xd01(i)*yd03(i) + xd02(i)*yd02(i)
                 + xd03(i)*yd01(i) + xd04(i)*yd00(i);
      s05 := s05 + xd00(i)*yd05(i) + xd01(i)*yd04(i) + xd02(i)*yd03(i)
                 + xd03(i)*yd02(i) + xd04(i)*yd01(i) + xd05(i)*yd00(i);
      s06 := s06 + xd00(i)*yd06(i) + xd01(i)*yd05(i) + xd02(i)*yd04(i)
                 + xd03(i)*yd03(i) + xd04(i)*yd02(i) + xd05(i)*yd01(i)
                 + xd06(i)*yd00(i);
      s07 := s07 + xd00(i)*yd07(i) + xd01(i)*yd06(i) + xd02(i)*yd05(i)
                 + xd03(i)*yd04(i) + xd04(i)*yd03(i) + xd05(i)*yd02(i)
                 + xd06(i)*yd01(i) + xd07(i)*yd00(i);
      s08a := s08a + xd00(i)*yd08(i) + xd01(i)*yd07(i) + xd02(i)*yd06(i)
                   + xd03(i)*yd05(i) + xd04(i)*yd04(i);
      s08b := s08b + xd05(i)*yd03(i) + xd06(i)*yd02(i) + xd07(i)*yd01(i)
                   + xd08(i)*yd00(i);
      s09a := s09a + xd00(i)*yd09(i) + xd01(i)*yd08(i) + xd02(i)*yd07(i)
                   + xd03(i)*yd06(i) + xd04(i)*yd05(i);
      s09b := s09b + xd05(i)*yd04(i) + xd06(i)*yd03(i) + xd07(i)*yd02(i)
                   + xd08(i)*yd01(i) + xd09(i)*yd00(i);
      s10a := s10a + xd00(i)*yd10(i) + xd01(i)*yd09(i) + xd02(i)*yd08(i)
                   + xd03(i)*yd07(i) + xd04(i)*yd06(i) + xd05(i)*yd05(i);
      s10b := s10b + xd06(i)*yd04(i) + xd07(i)*yd03(i) + xd08(i)*yd02(i)
                   + xd09(i)*yd01(i) + xd10(i)*yd00(i);
      s11a := s11a + xd00(i)*yd11(i) + xd01(i)*yd10(i) + xd02(i)*yd09(i)
                   + xd03(i)*yd08(i) + xd04(i)*yd07(i) + xd05(i)*yd06(i);
      s11b := s11b + xd06(i)*yd05(i) + xd07(i)*yd04(i) + xd08(i)*yd03(i)
                   + xd09(i)*yd02(i) + xd10(i)*yd01(i) + xd11(i)*yd00(i);
      s12a := s12a + xd00(i)*yd12(i) + xd01(i)*yd11(i) + xd02(i)*yd10(i)
                   + xd03(i)*yd09(i) + xd04(i)*yd08(i) + xd05(i)*yd07(i)
                   + xd06(i)*yd06(i);
      s12b := s12b + xd07(i)*yd05(i) + xd08(i)*yd04(i) + xd09(i)*yd03(i)
                   + xd10(i)*yd02(i) + xd11(i)*yd01(i) + xd12(i)*yd00(i);
      s13a := s13a + xd00(i)*yd13(i) + xd01(i)*yd12(i) + xd02(i)*yd11(i)
                   + xd03(i)*yd10(i) + xd04(i)*yd09(i) + xd05(i)*yd08(i)
                   + xd06(i)*yd07(i);
      s13b := s13b + xd07(i)*yd06(i) + xd08(i)*yd05(i) + xd09(i)*yd04(i)
                   + xd10(i)*yd03(i) + xd11(i)*yd02(i) + xd12(i)*yd01(i)
                   + xd13(i)*yd00(i);
      s14a := s14a + xd00(i)*yd14(i) + xd01(i)*yd13(i) + xd02(i)*yd12(i)
                   + xd03(i)*yd11(i) + xd04(i)*yd10(i) + xd05(i)*yd09(i)
                   + xd06(i)*yd08(i) + xd07(i)*yd07(i);
      s14b := s14b + xd08(i)*yd06(i) + xd09(i)*yd05(i) + xd10(i)*yd04(i)
                   + xd11(i)*yd03(i) + xd12(i)*yd02(i) + xd13(i)*yd01(i)
                   + xd14(i)*yd00(i);
      s15a := s15a + xd00(i)*yd15(i) + xd01(i)*yd14(i) + xd02(i)*yd13(i)
                   + xd03(i)*yd12(i) + xd04(i)*yd11(i) + xd05(i)*yd10(i)
                   + xd06(i)*yd09(i) + xd07(i)*yd08(i);
      s15b := s15b + xd08(i)*yd07(i) + xd09(i)*yd06(i) + xd10(i)*yd05(i)
                   + xd11(i)*yd04(i) + xd12(i)*yd03(i) + xd13(i)*yd02(i)
                   + xd14(i)*yd01(i) + xd15(i)*yd00(i);
      s16a := s16a + xd00(i)*yd16(i) + xd01(i)*yd15(i) + xd02(i)*yd14(i)
                   + xd03(i)*yd13(i) + xd04(i)*yd12(i) + xd05(i)*yd11(i);
      s16b := s16b + xd06(i)*yd10(i) + xd07(i)*yd09(i) + xd08(i)*yd08(i)
                   + xd09(i)*yd07(i) + xd10(i)*yd06(i) + xd11(i)*yd05(i);
      s16c := s16c + xd12(i)*yd04(i) + xd13(i)*yd03(i) + xd14(i)*yd02(i)
                   + xd15(i)*yd01(i) + xd16(i)*yd00(i);
      s17a := s17a + xd00(i)*yd17(i) + xd01(i)*yd16(i) + xd02(i)*yd15(i)
                   + xd03(i)*yd14(i) + xd04(i)*yd13(i) + xd05(i)*yd12(i);
      s17b := s17b + xd06(i)*yd11(i) + xd07(i)*yd10(i) + xd08(i)*yd09(i)
                   + xd09(i)*yd08(i) + xd10(i)*yd07(i) + xd11(i)*yd06(i);
      s17c := s17c + xd12(i)*yd05(i) + xd13(i)*yd04(i) + xd14(i)*yd03(i)
                   + xd15(i)*yd02(i) + xd16(i)*yd01(i) + xd17(i)*yd00(i);
      s18a := s18a + xd00(i)*yd18(i) + xd01(i)*yd17(i) + xd02(i)*yd16(i)
                   + xd03(i)*yd15(i) + xd04(i)*yd14(i) + xd05(i)*yd13(i)
                   + xd06(i)*yd12(i);
      s18b := s18b + xd07(i)*yd11(i) + xd08(i)*yd10(i) + xd09(i)*yd09(i)
                   + xd10(i)*yd08(i) + xd11(i)*yd07(i) + xd12(i)*yd06(i);
      s18c := s18c + xd13(i)*yd05(i) + xd14(i)*yd04(i) + xd15(i)*yd03(i)
                   + xd16(i)*yd02(i) + xd17(i)*yd01(i) + xd18(i)*yd00(i);
      s19a := s19a + xd00(i)*yd19(i) + xd01(i)*yd18(i) + xd02(i)*yd17(i)
                   + xd03(i)*yd16(i) + xd04(i)*yd15(i) + xd05(i)*yd14(i)
                   + xd06(i)*yd13(i);
      s19b := s19b + xd07(i)*yd12(i) + xd08(i)*yd11(i) + xd09(i)*yd10(i)
                   + xd10(i)*yd09(i) + xd11(i)*yd08(i) + xd12(i)*yd07(i)
                   + xd13(i)*yd06(i);
      s19c := s19c + xd14(i)*yd05(i) + xd15(i)*yd04(i) + xd16(i)*yd03(i)
                   + xd17(i)*yd02(i) + xd18(i)*yd01(i) + xd19(i)*yd00(i);
      s20a := s20a + xd00(i)*yd20(i) + xd01(i)*yd19(i) + xd02(i)*yd18(i)
                   + xd03(i)*yd17(i) + xd04(i)*yd16(i) + xd05(i)*yd15(i)
                   + xd06(i)*yd14(i);
      s20b := s20b + xd07(i)*yd13(i) + xd08(i)*yd12(i) + xd09(i)*yd11(i)
                   + xd10(i)*yd10(i) + xd11(i)*yd09(i) + xd12(i)*yd08(i)
                   + xd13(i)*yd07(i);
      s20c := s20c + xd14(i)*yd06(i) + xd15(i)*yd05(i) + xd16(i)*yd04(i)
                   + xd17(i)*yd03(i) + xd18(i)*yd02(i) + xd19(i)*yd01(i)
                   + xd20(i)*yd00(i);
      s21a := s21a + xd00(i)*yd21(i) + xd01(i)*yd20(i) + xd02(i)*yd19(i)
                   + xd03(i)*yd18(i) + xd04(i)*yd17(i) + xd05(i)*yd16(i)
                   + xd06(i)*yd15(i) + xd07(i)*yd14(i);
      s21b := s21b + xd08(i)*yd13(i) + xd09(i)*yd12(i) + xd10(i)*yd11(i)
                   + xd11(i)*yd10(i) + xd12(i)*yd09(i) + xd13(i)*yd08(i)
                   + xd14(i)*yd07(i);
      s21c := s21c + xd15(i)*yd06(i) + xd16(i)*yd05(i) + xd17(i)*yd04(i)
                   + xd18(i)*yd03(i) + xd19(i)*yd02(i) + xd20(i)*yd01(i)
                   + xd21(i)*yd00(i);
      s22a := s22a + xd00(i)*yd22(i) + xd01(i)*yd21(i) + xd02(i)*yd20(i)
                   + xd03(i)*yd19(i) + xd04(i)*yd18(i) + xd05(i)*yd17(i)
                   + xd06(i)*yd16(i) + xd07(i)*yd15(i);
      s22b := s22b + xd08(i)*yd14(i) + xd09(i)*yd13(i) + xd10(i)*yd12(i)
                   + xd11(i)*yd11(i) + xd12(i)*yd10(i) + xd13(i)*yd09(i)
                   + xd14(i)*yd08(i) + xd15(i)*yd07(i);
      s22c := s22c + xd16(i)*yd06(i) + xd17(i)*yd05(i) + xd18(i)*yd04(i)
                   + xd19(i)*yd03(i) + xd20(i)*yd02(i) + xd21(i)*yd01(i)
                   + xd22(i)*yd00(i);
      s23a := s23a + xd00(i)*yd23(i) + xd01(i)*yd22(i) + xd02(i)*yd21(i)
                   + xd03(i)*yd20(i) + xd04(i)*yd19(i) + xd05(i)*yd18(i)
                   + xd06(i)*yd17(i) + xd07(i)*yd16(i);
      s23b := s23b + xd08(i)*yd15(i) + xd09(i)*yd14(i) + xd10(i)*yd13(i)
                   + xd11(i)*yd12(i) + xd12(i)*yd11(i) + xd13(i)*yd10(i)
                   + xd14(i)*yd09(i) + xd15(i)*yd08(i);
      s23c := s23c + xd16(i)*yd07(i) + xd17(i)*yd06(i) + xd18(i)*yd05(i)
                   + xd19(i)*yd04(i) + xd20(i)*yd03(i) + xd21(i)*yd02(i)
                   + xd22(i)*yd01(i) + xd23(i)*yd00(i);
      s24a := s24a + xd00(i)*yd24(i) + xd01(i)*yd23(i) + xd02(i)*yd22(i)
                   + xd03(i)*yd21(i) + xd04(i)*yd20(i) + xd05(i)*yd19(i)
                   + xd06(i)*yd18(i);
      s24b := s24b + xd07(i)*yd17(i) + xd08(i)*yd16(i) + xd09(i)*yd15(i)
                   + xd10(i)*yd14(i) + xd11(i)*yd13(i) + xd12(i)*yd12(i);
      s24c := s24c + xd13(i)*yd11(i) + xd14(i)*yd10(i) + xd15(i)*yd09(i)
                   + xd16(i)*yd08(i) + xd17(i)*yd07(i) + xd18(i)*yd06(i);
      s24d := s24d + xd19(i)*yd05(i) + xd20(i)*yd04(i) + xd21(i)*yd03(i)
                   + xd22(i)*yd02(i) + xd23(i)*yd01(i) + xd24(i)*yd00(i);
      s25a := s25a + xd00(i)*yd25(i) + xd01(i)*yd24(i) + xd02(i)*yd23(i)
                   + xd03(i)*yd22(i) + xd04(i)*yd21(i) + xd05(i)*yd20(i)
                   + xd06(i)*yd19(i);
      s25b := s25b + xd07(i)*yd18(i) + xd08(i)*yd17(i) + xd09(i)*yd16(i)
                   + xd10(i)*yd15(i) + xd11(i)*yd14(i) + xd12(i)*yd13(i)
                   + xd13(i)*yd12(i);
      s25c := s25c + xd14(i)*yd11(i) + xd15(i)*yd10(i) + xd16(i)*yd09(i)
                   + xd17(i)*yd08(i) + xd18(i)*yd07(i) + xd19(i)*yd06(i);
      s25d := s25d + xd20(i)*yd05(i) + xd21(i)*yd04(i) + xd22(i)*yd03(i)
                   + xd23(i)*yd02(i) + xd24(i)*yd01(i) + xd25(i)*yd00(i);
      s26a := s26a + xd00(i)*yd26(i) + xd01(i)*yd25(i) + xd02(i)*yd24(i)
                   + xd03(i)*yd23(i) + xd04(i)*yd22(i) + xd05(i)*yd21(i)
                   + xd06(i)*yd20(i);
      s26b := s26b + xd07(i)*yd19(i) + xd08(i)*yd18(i) + xd09(i)*yd17(i)
                   + xd10(i)*yd16(i) + xd11(i)*yd15(i) + xd12(i)*yd14(i)
                   + xd13(i)*yd13(i);
      s26c := s26c + xd14(i)*yd12(i) + xd15(i)*yd11(i) + xd16(i)*yd10(i)
                   + xd17(i)*yd09(i) + xd18(i)*yd08(i) + xd19(i)*yd07(i)
                   + xd20(i)*yd06(i);
      s26d := s26d + xd21(i)*yd05(i) + xd22(i)*yd04(i) + xd23(i)*yd03(i)
                   + xd24(i)*yd02(i) + xd25(i)*yd01(i) + xd26(i)*yd00(i);
      s27a := s27a + xd00(i)*yd27(i) + xd01(i)*yd26(i) + xd02(i)*yd25(i)
                   + xd03(i)*yd24(i) + xd04(i)*yd23(i) + xd05(i)*yd22(i)
                   + xd06(i)*yd21(i);
      s27b := s27b + xd07(i)*yd20(i) + xd08(i)*yd19(i) + xd09(i)*yd18(i)
                   + xd10(i)*yd17(i) + xd11(i)*yd16(i) + xd12(i)*yd15(i)
                   + xd13(i)*yd14(i);
      s27c := s27c + xd14(i)*yd13(i) + xd15(i)*yd12(i) + xd16(i)*yd11(i)
                   + xd17(i)*yd10(i) + xd18(i)*yd09(i) + xd19(i)*yd08(i)
                   + xd20(i)*yd07(i);
      s27d := s27d + xd21(i)*yd06(i) + xd22(i)*yd05(i) + xd23(i)*yd04(i)
                   + xd24(i)*yd03(i) + xd25(i)*yd02(i) + xd26(i)*yd01(i)
                   + xd27(i)*yd00(i);
      s28a := s28a + xd00(i)*yd28(i) + xd01(i)*yd27(i) + xd02(i)*yd26(i)
                   + xd03(i)*yd25(i) + xd04(i)*yd24(i) + xd05(i)*yd23(i)
                   + xd06(i)*yd22(i) + xd07(i)*yd21(i);
      s28b := s28b + xd08(i)*yd20(i) + xd09(i)*yd19(i) + xd10(i)*yd18(i)
                   + xd11(i)*yd17(i) + xd12(i)*yd16(i) + xd13(i)*yd15(i)
                   + xd14(i)*yd14(i);
      s28c := s28c + xd15(i)*yd13(i) + xd16(i)*yd12(i) + xd17(i)*yd11(i)
                   + xd18(i)*yd10(i) + xd19(i)*yd09(i) + xd20(i)*yd08(i)
                   + xd21(i)*yd07(i);
      s28d := s28d + xd22(i)*yd06(i) + xd23(i)*yd05(i) + xd24(i)*yd04(i)
                   + xd25(i)*yd03(i) + xd26(i)*yd02(i) + xd27(i)*yd01(i)
                   + xd28(i)*yd00(i);
      s29a := s29a + xd00(i)*yd29(i) + xd01(i)*yd28(i) + xd02(i)*yd27(i)
                   + xd03(i)*yd26(i) + xd04(i)*yd25(i) + xd05(i)*yd24(i)
                   + xd06(i)*yd23(i) + xd07(i)*yd22(i);
      s29b := s29b + xd08(i)*yd21(i) + xd09(i)*yd20(i) + xd10(i)*yd19(i)
                   + xd11(i)*yd18(i) + xd12(i)*yd17(i) + xd13(i)*yd16(i)
                   + xd14(i)*yd15(i) + xd15(i)*yd14(i);
      s29c := s29c + xd16(i)*yd13(i) + xd17(i)*yd12(i) + xd18(i)*yd11(i)
                   + xd19(i)*yd10(i) + xd20(i)*yd09(i) + xd21(i)*yd08(i)
                   + xd22(i)*yd07(i);
      s29d := s29d + xd23(i)*yd06(i) + xd24(i)*yd05(i) + xd25(i)*yd04(i)
                   + xd26(i)*yd03(i) + xd27(i)*yd02(i) + xd28(i)*yd01(i)
                   + xd29(i)*yd00(i);
      s30a := s30a + xd00(i)*yd30(i) + xd01(i)*yd29(i) + xd02(i)*yd28(i)
                   + xd03(i)*yd27(i) + xd04(i)*yd26(i) + xd05(i)*yd25(i)
                   + xd06(i)*yd24(i) + xd07(i)*yd23(i);
      s30b := s30b + xd08(i)*yd22(i) + xd09(i)*yd21(i) + xd10(i)*yd20(i)
                   + xd11(i)*yd19(i) + xd12(i)*yd18(i) + xd13(i)*yd17(i)
                   + xd14(i)*yd16(i) + xd15(i)*yd15(i);
      s30c := s30c + xd16(i)*yd14(i) + xd17(i)*yd13(i) + xd18(i)*yd12(i)
                   + xd19(i)*yd11(i) + xd20(i)*yd10(i) + xd21(i)*yd09(i)
                   + xd22(i)*yd08(i) + xd23(i)*yd07(i);
      s30d := s30d + xd24(i)*yd06(i) + xd25(i)*yd05(i) + xd26(i)*yd04(i)
                   + xd27(i)*yd03(i) + xd28(i)*yd02(i) + xd29(i)*yd01(i)
                   + xd30(i)*yd00(i);
      s31a := s31a + xd00(i)*yd31(i) + xd01(i)*yd30(i) + xd02(i)*yd29(i)
                   + xd03(i)*yd28(i) + xd04(i)*yd27(i) + xd05(i)*yd26(i)
                   + xd06(i)*yd25(i) + xd07(i)*yd24(i);
      s31b := s31b + xd08(i)*yd23(i) + xd09(i)*yd22(i) + xd10(i)*yd21(i)
                   + xd11(i)*yd20(i) + xd12(i)*yd19(i) + xd13(i)*yd18(i)
                   + xd14(i)*yd17(i) + xd15(i)*yd16(i);
      s31c := s31c + xd16(i)*yd15(i) + xd17(i)*yd14(i) + xd18(i)*yd13(i)
                   + xd19(i)*yd12(i) + xd20(i)*yd11(i) + xd21(i)*yd10(i)
                   + xd22(i)*yd09(i) + xd23(i)*yd08(i);
      s31d := s31d + xd24(i)*yd07(i) + xd25(i)*yd06(i) + xd26(i)*yd05(i)
                   + xd27(i)*yd04(i) + xd28(i)*yd03(i) + xd29(i)*yd02(i)
                   + xd30(i)*yd01(i) + xd31(i)*yd00(i);
    end loop;
    if nd > 0 then
      res := res + to_octo_double(s00,s01,s02,s03,s04,s05,s06,s07,
                                  s08a,s09a,s10a,s11a,s12a,s13a,s14a,s15a,
                                  s16a,s17a,s18a,s19a,s20a,s21a,s22a,s23a,
                                  s24a,s25a,s26a,s27a,s28a,s29a,s30a,s31a,
                                  verbose=>true)
                 + to_octo_double(s08b,s09b,s10b,s11b,s12b,s13b,s14b,s15b,
                                  s16b,s17b,s18b,s19b,s20b,s21b,s22b,s23b,
                                  s24b,s25b,s26b,s27b,s28b,s29b,s30b,s31b,
                                  verbose=>true)
                 + to_octo_double(s16c,s17c,s18c,s19c,s20c,s21c,s22c,s23c,
                                  s24c,s25c,s26c,s27c,s28c,s29c,s30c,s31c,
                                  verbose=>true)
                 + to_octo_double(s24d,s25d,s26d,s27d,s28d,s29d,s30d,s31d,
                                  verbose=>true);
    end if;
    return res;
  end Product;

  function Product ( x,y : OctoDobl_Complex_Vectors.Vector;
                     verbose : boolean := true ) return Complex_Number is

    res : Complex_Number;
    resre,resim : octo_double;
    xre,xim : Octo_Double_Vectors.Vector(x'range);
    yre,yim : Octo_Double_Vectors.Vector(y'range);

  begin
    for i in x'range loop
      xre(i) := OctoDobl_Complex_Numbers.REAL_PART(x(i));
      xim(i) := OctoDobl_Complex_Numbers.IMAG_PART(x(i));
      yre(i) := OctoDobl_Complex_Numbers.REAL_PART(y(i));
      yim(i) := OctoDobl_Complex_Numbers.IMAG_PART(y(i));
    end loop;
    resre := Product(xre,yre,verbose) - Product(xim,yim,verbose);
    resim := Product(xre,yim,verbose) + Product(xim,yre,verbose);
    res := create(resre,resim);
    return res;
  end Product;

end Vectored_Octo_Doubles;
