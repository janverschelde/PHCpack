with text_io;                            use text_io;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;
with Bits_of_Doubles;

package body Vectored_Hexa_Doubles is

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
                x32,x33,x34,x35 : in Standard_Floating_Vectors.Vector;
                x36,x37,x38,x39 : in Standard_Floating_Vectors.Vector;
                x40,x41,x42,x43 : in Standard_Floating_Vectors.Vector;
                x44,x45,x46,x47 : in Standard_Floating_Vectors.Vector;
                x48,x49,x50,x51 : in Standard_Floating_Vectors.Vector;
                x52,x53,x54,x55 : in Standard_Floating_Vectors.Vector;
                x56,x57,x58,x59 : in Standard_Floating_Vectors.Vector;
                x60,x61,x62,x63 : in Standard_Floating_Vectors.Vector;
                y00,y01,y02,y03 : in Standard_Floating_Vectors.Vector;
                y04,y05,y06,y07 : in Standard_Floating_Vectors.Vector;
                y08,y09,y10,y11 : in Standard_Floating_Vectors.Vector;
                y12,y13,y14,y15 : in Standard_Floating_Vectors.Vector;
                y16,y17,y18,y19 : in Standard_Floating_Vectors.Vector;
                y20,y21,y22,y23 : in Standard_Floating_Vectors.Vector;
                y24,y25,y26,y27 : in Standard_Floating_Vectors.Vector;
                y28,y29,y30,y31 : in Standard_Floating_Vectors.Vector;
                y32,y33,y34,y35 : in Standard_Floating_Vectors.Vector;
                y36,y37,y38,y39 : in Standard_Floating_Vectors.Vector;
                y40,y41,y42,y43 : in Standard_Floating_Vectors.Vector;
                y44,y45,y46,y47 : in Standard_Floating_Vectors.Vector;
                y48,y49,y50,y51 : in Standard_Floating_Vectors.Vector;
                y52,y53,y54,y55 : in Standard_Floating_Vectors.Vector;
                y56,y57,y58,y59 : in Standard_Floating_Vectors.Vector;
                y60,y61,y62,y63 : in Standard_Floating_Vectors.Vector;
                s00,s01,s02,s03,s04,s05,s06,s07 : out double_float;
                s08,s09,s10,s11,s12,s13,s14,s15 : out double_float;
                s16,s17,s18,s19,s20,s21,s22,s23 : out double_float;
                s24,s25,s26,s27,s28,s29,s30,s31 : out double_float;
                s32,s33,s34,s35,s36,s37,s38,s39 : out double_float;
                s40,s41,s42,s43,s44,s45,s46,s47 : out double_float;
                s48,s49,s50,s51,s52,s53,s54,s55 : out double_float;
                s56,s57,s58,s59,s60,s61,s62,s63 : out double_float ) is
  begin
    s00 := 0.0; s01 := 0.0; s02 := 0.0; s03 := 0.0;
    s04 := 0.0; s05 := 0.0; s06 := 0.0; s07 := 0.0;
    s08 := 0.0; s09 := 0.0; s10 := 0.0; s11 := 0.0;
    s12 := 0.0; s13 := 0.0; s14 := 0.0; s15 := 0.0;
    s16 := 0.0; s17 := 0.0; s18 := 0.0; s19 := 0.0;
    s20 := 0.0; s21 := 0.0; s22 := 0.0; s23 := 0.0;
    s24 := 0.0; s25 := 0.0; s26 := 0.0; s27 := 0.0;
    s28 := 0.0; s29 := 0.0; s30 := 0.0; s31 := 0.0;
    s32 := 0.0; s33 := 0.0; s34 := 0.0; s35 := 0.0;
    s36 := 0.0; s37 := 0.0; s38 := 0.0; s39 := 0.0;
    s40 := 0.0; s41 := 0.0; s42 := 0.0; s43 := 0.0;
    s44 := 0.0; s45 := 0.0; s46 := 0.0; s47 := 0.0;
    s48 := 0.0; s49 := 0.0; s50 := 0.0; s51 := 0.0;
    s52 := 0.0; s53 := 0.0; s54 := 0.0; s55 := 0.0;
    s56 := 0.0; s57 := 0.0; s58 := 0.0; s59 := 0.0;
    s60 := 0.0; s61 := 0.0; s62 := 0.0; s63 := 0.0;
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
      s32 := s32 + x00(i)*y32(i) + x01(i)*y31(i) + x02(i)*y30(i)
                 + x03(i)*y29(i) + x04(i)*y28(i) + x05(i)*y27(i)
                 + x06(i)*y26(i) + x07(i)*y25(i) + x08(i)*y24(i)
                 + x09(i)*y23(i) + x10(i)*y22(i) + x11(i)*y21(i)
                 + x12(i)*y20(i) + x13(i)*y19(i) + x14(i)*y18(i)
                 + x15(i)*y17(i) + x16(i)*y16(i) + x17(i)*y15(i)
                 + x18(i)*y14(i) + x19(i)*y13(i) + x20(i)*y12(i)
                 + x21(i)*y11(i) + x22(i)*y10(i) + x23(i)*y09(i)
                 + x24(i)*y08(i) + x25(i)*y07(i) + x26(i)*y06(i)
                 + x27(i)*y05(i) + x28(i)*y04(i) + x29(i)*y03(i)
                 + x30(i)*y02(i) + x31(i)*y01(i) + x32(i)*y00(i);
      s33 := s33 + x00(i)*y33(i) + x01(i)*y32(i) + x02(i)*y31(i)
                 + x03(i)*y30(i) + x04(i)*y29(i) + x05(i)*y28(i)
                 + x06(i)*y27(i) + x07(i)*y26(i) + x08(i)*y25(i)
                 + x09(i)*y24(i) + x10(i)*y23(i) + x11(i)*y22(i)
                 + x12(i)*y21(i) + x13(i)*y20(i) + x14(i)*y19(i)
                 + x15(i)*y18(i) + x16(i)*y17(i) + x17(i)*y16(i)
                 + x18(i)*y15(i) + x19(i)*y14(i) + x20(i)*y13(i)
                 + x21(i)*y12(i) + x22(i)*y11(i) + x23(i)*y10(i)
                 + x24(i)*y09(i) + x25(i)*y08(i) + x26(i)*y07(i)
                 + x27(i)*y06(i) + x28(i)*y05(i) + x29(i)*y04(i)
                 + x30(i)*y03(i) + x31(i)*y02(i) + x32(i)*y01(i)
                 + x33(i)*y00(i);
      s34 := s34 + x00(i)*y34(i) + x01(i)*y33(i) + x02(i)*y32(i)
                 + x03(i)*y31(i) + x04(i)*y30(i) + x05(i)*y29(i)
                 + x06(i)*y28(i) + x07(i)*y27(i) + x08(i)*y26(i)
                 + x09(i)*y25(i) + x10(i)*y24(i) + x11(i)*y23(i)
                 + x12(i)*y22(i) + x13(i)*y21(i) + x14(i)*y20(i)
                 + x15(i)*y19(i) + x16(i)*y18(i) + x17(i)*y17(i)
                 + x18(i)*y16(i) + x19(i)*y15(i) + x20(i)*y14(i)
                 + x21(i)*y13(i) + x22(i)*y12(i) + x23(i)*y11(i)
                 + x24(i)*y10(i) + x25(i)*y09(i) + x26(i)*y08(i)
                 + x27(i)*y07(i) + x28(i)*y06(i) + x29(i)*y05(i)
                 + x30(i)*y04(i) + x31(i)*y03(i) + x32(i)*y02(i)
                 + x33(i)*y01(i) + x34(i)*y00(i);
      s35 := s35 + x00(i)*y35(i) + x01(i)*y34(i) + x02(i)*y33(i)
                 + x03(i)*y32(i) + x04(i)*y31(i) + x05(i)*y30(i)
                 + x06(i)*y29(i) + x07(i)*y28(i) + x08(i)*y27(i)
                 + x09(i)*y26(i) + x10(i)*y25(i) + x11(i)*y24(i)
                 + x12(i)*y23(i) + x13(i)*y22(i) + x14(i)*y21(i)
                 + x15(i)*y20(i) + x16(i)*y19(i) + x17(i)*y18(i)
                 + x18(i)*y17(i) + x19(i)*y16(i) + x20(i)*y15(i)
                 + x21(i)*y14(i) + x22(i)*y13(i) + x23(i)*y12(i)
                 + x24(i)*y11(i) + x25(i)*y10(i) + x26(i)*y09(i)
                 + x27(i)*y08(i) + x28(i)*y07(i) + x29(i)*y06(i)
                 + x30(i)*y05(i) + x31(i)*y04(i) + x32(i)*y03(i)
                 + x33(i)*y02(i) + x34(i)*y01(i) + x35(i)*y00(i);
      s36 := s36 + x00(i)*y36(i) + x01(i)*y35(i) + x02(i)*y34(i)
                 + x03(i)*y33(i) + x04(i)*y32(i) + x05(i)*y31(i)
                 + x06(i)*y30(i) + x07(i)*y29(i) + x08(i)*y28(i)
                 + x09(i)*y27(i) + x10(i)*y26(i) + x11(i)*y25(i)
                 + x12(i)*y24(i) + x13(i)*y23(i) + x14(i)*y22(i)
                 + x15(i)*y21(i) + x16(i)*y20(i) + x17(i)*y19(i)
                 + x18(i)*y18(i) + x19(i)*y17(i) + x20(i)*y16(i)
                 + x21(i)*y15(i) + x22(i)*y14(i) + x23(i)*y13(i)
                 + x24(i)*y12(i) + x25(i)*y11(i) + x26(i)*y10(i)
                 + x27(i)*y09(i) + x28(i)*y08(i) + x29(i)*y07(i)
                 + x30(i)*y06(i) + x31(i)*y05(i) + x32(i)*y04(i)
                 + x33(i)*y03(i) + x34(i)*y02(i) + x35(i)*y01(i)
                 + x36(i)*y00(i);
      s37 := s37 + x00(i)*y37(i) + x01(i)*y36(i) + x02(i)*y35(i)
                 + x03(i)*y34(i) + x04(i)*y33(i) + x05(i)*y32(i)
                 + x06(i)*y31(i) + x07(i)*y30(i) + x08(i)*y29(i)
                 + x09(i)*y28(i) + x10(i)*y27(i) + x11(i)*y26(i)
                 + x12(i)*y25(i) + x13(i)*y24(i) + x14(i)*y23(i)
                 + x15(i)*y22(i) + x16(i)*y21(i) + x17(i)*y20(i)
                 + x18(i)*y19(i) + x19(i)*y18(i) + x20(i)*y17(i)
                 + x21(i)*y16(i) + x22(i)*y15(i) + x23(i)*y14(i)
                 + x24(i)*y13(i) + x25(i)*y12(i) + x26(i)*y11(i)
                 + x27(i)*y10(i) + x28(i)*y09(i) + x29(i)*y08(i)
                 + x30(i)*y07(i) + x31(i)*y06(i) + x32(i)*y05(i)
                 + x33(i)*y04(i) + x34(i)*y03(i) + x35(i)*y02(i)
                 + x36(i)*y01(i) + x37(i)*y00(i);
      s38 := s38 + x00(i)*y38(i) + x01(i)*y37(i) + x02(i)*y36(i)
                 + x03(i)*y35(i) + x04(i)*y34(i) + x05(i)*y33(i)
                 + x06(i)*y32(i) + x07(i)*y31(i) + x08(i)*y30(i)
                 + x09(i)*y29(i) + x10(i)*y28(i) + x11(i)*y27(i)
                 + x12(i)*y26(i) + x13(i)*y25(i) + x14(i)*y24(i)
                 + x15(i)*y23(i) + x16(i)*y22(i) + x17(i)*y21(i)
                 + x18(i)*y20(i) + x19(i)*y19(i) + x20(i)*y18(i)
                 + x21(i)*y17(i) + x22(i)*y16(i) + x23(i)*y15(i)
                 + x24(i)*y14(i) + x25(i)*y13(i) + x26(i)*y12(i)
                 + x27(i)*y11(i) + x28(i)*y10(i) + x29(i)*y09(i)
                 + x30(i)*y08(i) + x31(i)*y07(i) + x32(i)*y06(i)
                 + x33(i)*y05(i) + x34(i)*y04(i) + x35(i)*y03(i)
                 + x36(i)*y02(i) + x37(i)*y01(i) + x38(i)*y00(i);
      s39 := s39 + x00(i)*y39(i) + x01(i)*y38(i) + x02(i)*y37(i)
                 + x03(i)*y36(i) + x04(i)*y35(i) + x05(i)*y34(i)
                 + x06(i)*y33(i) + x07(i)*y32(i) + x08(i)*y31(i)
                 + x09(i)*y30(i) + x10(i)*y29(i) + x11(i)*y28(i)
                 + x12(i)*y27(i) + x13(i)*y26(i) + x14(i)*y25(i)
                 + x15(i)*y24(i) + x16(i)*y23(i) + x17(i)*y22(i)
                 + x18(i)*y21(i) + x19(i)*y20(i) + x20(i)*y19(i)
                 + x21(i)*y18(i) + x22(i)*y17(i) + x23(i)*y16(i)
                 + x24(i)*y15(i) + x25(i)*y14(i) + x26(i)*y13(i)
                 + x27(i)*y12(i) + x28(i)*y11(i) + x29(i)*y10(i)
                 + x30(i)*y09(i) + x31(i)*y08(i) + x32(i)*y07(i)
                 + x33(i)*y06(i) + x34(i)*y05(i) + x35(i)*y04(i)
                 + x36(i)*y03(i) + x37(i)*y02(i) + x38(i)*y01(i)
                 + x39(i)*y00(i);
      s40 := s40 + x00(i)*y40(i) + x01(i)*y39(i) + x02(i)*y38(i)
                 + x03(i)*y37(i) + x04(i)*y36(i) + x05(i)*y35(i)
                 + x06(i)*y34(i) + x07(i)*y33(i) + x08(i)*y32(i)
                 + x09(i)*y31(i) + x10(i)*y30(i) + x11(i)*y29(i)
                 + x12(i)*y28(i) + x13(i)*y27(i) + x14(i)*y26(i)
                 + x15(i)*y25(i) + x16(i)*y24(i) + x17(i)*y23(i)
                 + x18(i)*y22(i) + x19(i)*y21(i) + x20(i)*y20(i)
                 + x21(i)*y19(i) + x22(i)*y18(i) + x23(i)*y17(i)
                 + x24(i)*y16(i) + x25(i)*y15(i) + x26(i)*y14(i)
                 + x27(i)*y13(i) + x28(i)*y12(i) + x29(i)*y11(i)
                 + x30(i)*y10(i) + x31(i)*y09(i) + x32(i)*y08(i)
                 + x33(i)*y07(i) + x34(i)*y06(i) + x35(i)*y05(i)
                 + x36(i)*y04(i) + x37(i)*y03(i) + x38(i)*y02(i)
                 + x39(i)*y01(i) + x40(i)*y00(i);
      s41 := s41 + x00(i)*y41(i) + x01(i)*y40(i) + x02(i)*y39(i)
                 + x03(i)*y38(i) + x04(i)*y37(i) + x05(i)*y36(i)
                 + x06(i)*y35(i) + x07(i)*y34(i) + x08(i)*y33(i)
                 + x09(i)*y32(i) + x10(i)*y31(i) + x11(i)*y30(i)
                 + x12(i)*y29(i) + x13(i)*y28(i) + x14(i)*y27(i)
                 + x15(i)*y26(i) + x16(i)*y25(i) + x17(i)*y24(i)
                 + x18(i)*y23(i) + x19(i)*y22(i) + x20(i)*y21(i)
                 + x21(i)*y20(i) + x22(i)*y19(i) + x23(i)*y18(i)
                 + x24(i)*y17(i) + x25(i)*y16(i) + x26(i)*y15(i)
                 + x27(i)*y14(i) + x28(i)*y13(i) + x29(i)*y12(i)
                 + x30(i)*y11(i) + x31(i)*y10(i) + x32(i)*y09(i)
                 + x33(i)*y08(i) + x34(i)*y07(i) + x35(i)*y06(i)
                 + x36(i)*y05(i) + x37(i)*y04(i) + x38(i)*y03(i)
                 + x39(i)*y02(i) + x40(i)*y01(i) + x41(i)*y00(i);
      s42 := s42 + x00(i)*y42(i) + x01(i)*y41(i) + x02(i)*y40(i)
                 + x03(i)*y39(i) + x04(i)*y38(i) + x05(i)*y37(i)
                 + x06(i)*y36(i) + x07(i)*y35(i) + x08(i)*y34(i)
                 + x09(i)*y33(i) + x10(i)*y32(i) + x11(i)*y31(i)
                 + x12(i)*y30(i) + x13(i)*y29(i) + x14(i)*y28(i)
                 + x15(i)*y27(i) + x16(i)*y26(i) + x17(i)*y25(i)
                 + x18(i)*y24(i) + x19(i)*y23(i) + x20(i)*y22(i)
                 + x21(i)*y21(i) + x22(i)*y20(i) + x23(i)*y19(i)
                 + x24(i)*y18(i) + x25(i)*y17(i) + x26(i)*y16(i)
                 + x27(i)*y15(i) + x28(i)*y14(i) + x29(i)*y13(i)
                 + x30(i)*y12(i) + x31(i)*y11(i) + x32(i)*y10(i)
                 + x33(i)*y09(i) + x34(i)*y08(i) + x35(i)*y07(i)
                 + x36(i)*y06(i) + x37(i)*y05(i) + x38(i)*y04(i)
                 + x39(i)*y03(i) + x40(i)*y02(i) + x41(i)*y01(i)
                 + x42(i)*y00(i);
      s43 := s43 + x00(i)*y43(i) + x01(i)*y42(i) + x02(i)*y41(i)
                 + x03(i)*y40(i) + x04(i)*y39(i) + x05(i)*y38(i)
                 + x06(i)*y37(i) + x07(i)*y36(i) + x08(i)*y35(i)
                 + x09(i)*y34(i) + x10(i)*y33(i) + x11(i)*y32(i)
                 + x12(i)*y31(i) + x13(i)*y30(i) + x14(i)*y29(i)
                 + x15(i)*y28(i) + x16(i)*y27(i) + x17(i)*y26(i)
                 + x18(i)*y25(i) + x19(i)*y24(i) + x20(i)*y23(i)
                 + x21(i)*y22(i) + x22(i)*y21(i) + x23(i)*y20(i)
                 + x24(i)*y19(i) + x25(i)*y18(i) + x26(i)*y17(i)
                 + x27(i)*y16(i) + x28(i)*y15(i) + x29(i)*y14(i)
                 + x30(i)*y13(i) + x31(i)*y12(i) + x32(i)*y11(i)
                 + x33(i)*y10(i) + x34(i)*y09(i) + x35(i)*y08(i)
                 + x36(i)*y07(i) + x37(i)*y06(i) + x38(i)*y05(i)
                 + x39(i)*y04(i) + x40(i)*y03(i) + x41(i)*y02(i)
                 + x42(i)*y01(i) + x43(i)*y00(i);
      s44 := s44 + x00(i)*y44(i) + x01(i)*y43(i) + x02(i)*y42(i)
                 + x03(i)*y41(i) + x04(i)*y40(i) + x05(i)*y39(i)
                 + x06(i)*y38(i) + x07(i)*y37(i) + x08(i)*y36(i)
                 + x09(i)*y35(i) + x10(i)*y34(i) + x11(i)*y33(i)
                 + x12(i)*y32(i) + x13(i)*y31(i) + x14(i)*y30(i)
                 + x15(i)*y29(i) + x16(i)*y28(i) + x17(i)*y27(i)
                 + x18(i)*y26(i) + x19(i)*y25(i) + x20(i)*y24(i)
                 + x21(i)*y23(i) + x22(i)*y22(i) + x23(i)*y21(i)
                 + x24(i)*y20(i) + x25(i)*y19(i) + x26(i)*y18(i)
                 + x27(i)*y17(i) + x28(i)*y16(i) + x29(i)*y15(i)
                 + x30(i)*y14(i) + x31(i)*y13(i) + x32(i)*y12(i)
                 + x33(i)*y11(i) + x34(i)*y10(i) + x35(i)*y09(i)
                 + x36(i)*y08(i) + x37(i)*y07(i) + x38(i)*y06(i)
                 + x39(i)*y05(i) + x40(i)*y04(i) + x41(i)*y03(i)
                 + x42(i)*y02(i) + x43(i)*y01(i) + x44(i)*y00(i);
      s45 := s45 + x00(i)*y45(i) + x01(i)*y44(i) + x02(i)*y43(i)
                 + x03(i)*y42(i) + x04(i)*y41(i) + x05(i)*y40(i)
                 + x06(i)*y39(i) + x07(i)*y38(i) + x08(i)*y37(i)
                 + x09(i)*y36(i) + x10(i)*y35(i) + x11(i)*y34(i)
                 + x12(i)*y33(i) + x13(i)*y32(i) + x14(i)*y31(i)
                 + x15(i)*y30(i) + x16(i)*y29(i) + x17(i)*y28(i)
                 + x18(i)*y27(i) + x19(i)*y26(i) + x20(i)*y25(i)
                 + x21(i)*y24(i) + x22(i)*y23(i) + x23(i)*y22(i)
                 + x24(i)*y21(i) + x25(i)*y20(i) + x26(i)*y19(i)
                 + x27(i)*y18(i) + x28(i)*y17(i) + x29(i)*y16(i)
                 + x30(i)*y15(i) + x31(i)*y14(i) + x32(i)*y13(i)
                 + x33(i)*y12(i) + x34(i)*y11(i) + x35(i)*y10(i)
                 + x36(i)*y09(i) + x37(i)*y08(i) + x38(i)*y07(i)
                 + x39(i)*y06(i) + x40(i)*y05(i) + x41(i)*y04(i)
                 + x42(i)*y03(i) + x43(i)*y02(i) + x44(i)*y01(i)
                 + x45(i)*y00(i);
      s46 := s46 + x00(i)*y46(i) + x01(i)*y45(i) + x02(i)*y44(i)
                 + x03(i)*y43(i) + x04(i)*y42(i) + x05(i)*y41(i)
                 + x06(i)*y40(i) + x07(i)*y39(i) + x08(i)*y38(i)
                 + x09(i)*y37(i) + x10(i)*y36(i) + x11(i)*y35(i)
                 + x12(i)*y34(i) + x13(i)*y33(i) + x14(i)*y32(i)
                 + x15(i)*y31(i) + x16(i)*y30(i) + x17(i)*y29(i)
                 + x18(i)*y28(i) + x19(i)*y27(i) + x20(i)*y26(i)
                 + x21(i)*y25(i) + x22(i)*y24(i) + x23(i)*y23(i)
                 + x24(i)*y22(i) + x25(i)*y21(i) + x26(i)*y20(i)
                 + x27(i)*y19(i) + x28(i)*y18(i) + x29(i)*y17(i)
                 + x30(i)*y16(i) + x31(i)*y15(i) + x32(i)*y14(i)
                 + x33(i)*y13(i) + x34(i)*y12(i) + x35(i)*y11(i)
                 + x36(i)*y10(i) + x37(i)*y09(i) + x38(i)*y08(i)
                 + x39(i)*y07(i) + x40(i)*y06(i) + x41(i)*y05(i)
                 + x42(i)*y04(i) + x43(i)*y03(i) + x44(i)*y02(i)
                 + x45(i)*y01(i) + x46(i)*y00(i);
      s47 := s47 + x00(i)*y47(i) + x01(i)*y46(i) + x02(i)*y45(i)
                 + x03(i)*y44(i) + x04(i)*y43(i) + x05(i)*y42(i)
                 + x06(i)*y41(i) + x07(i)*y40(i) + x08(i)*y39(i)
                 + x09(i)*y38(i) + x10(i)*y37(i) + x11(i)*y36(i)
                 + x12(i)*y35(i) + x13(i)*y34(i) + x14(i)*y33(i)
                 + x15(i)*y32(i) + x16(i)*y31(i) + x17(i)*y30(i)
                 + x18(i)*y29(i) + x19(i)*y28(i) + x20(i)*y27(i)
                 + x21(i)*y26(i) + x22(i)*y25(i) + x23(i)*y24(i)
                 + x24(i)*y23(i) + x25(i)*y22(i) + x26(i)*y21(i)
                 + x27(i)*y20(i) + x28(i)*y19(i) + x29(i)*y18(i)
                 + x30(i)*y17(i) + x31(i)*y16(i) + x32(i)*y15(i)
                 + x33(i)*y14(i) + x34(i)*y13(i) + x35(i)*y12(i)
                 + x36(i)*y11(i) + x37(i)*y10(i) + x38(i)*y09(i)
                 + x39(i)*y08(i) + x40(i)*y07(i) + x41(i)*y06(i)
                 + x42(i)*y05(i) + x43(i)*y04(i) + x44(i)*y03(i)
                 + x45(i)*y02(i) + x46(i)*y01(i) + x47(i)*y00(i);
      s48 := s48 + x00(i)*y48(i) + x01(i)*y47(i) + x02(i)*y46(i)
                 + x03(i)*y45(i) + x04(i)*y44(i) + x05(i)*y43(i)
                 + x06(i)*y42(i) + x07(i)*y41(i) + x08(i)*y40(i)
                 + x09(i)*y39(i) + x10(i)*y38(i) + x11(i)*y37(i)
                 + x12(i)*y36(i) + x13(i)*y35(i) + x14(i)*y34(i)
                 + x15(i)*y33(i) + x16(i)*y32(i) + x17(i)*y31(i)
                 + x18(i)*y30(i) + x19(i)*y29(i) + x20(i)*y28(i)
                 + x21(i)*y27(i) + x22(i)*y26(i) + x23(i)*y25(i)
                 + x24(i)*y24(i) + x25(i)*y23(i) + x26(i)*y22(i)
                 + x27(i)*y21(i) + x28(i)*y20(i) + x29(i)*y19(i)
                 + x30(i)*y18(i) + x31(i)*y17(i) + x32(i)*y16(i)
                 + x33(i)*y15(i) + x34(i)*y14(i) + x35(i)*y13(i)
                 + x36(i)*y12(i) + x37(i)*y11(i) + x38(i)*y10(i)
                 + x39(i)*y09(i) + x40(i)*y08(i) + x41(i)*y07(i)
                 + x42(i)*y06(i) + x43(i)*y05(i) + x44(i)*y04(i)
                 + x45(i)*y03(i) + x46(i)*y02(i) + x47(i)*y01(i)
                 + x48(i)*y00(i);
      s49 := s49 + x00(i)*y49(i) + x01(i)*y48(i) + x02(i)*y47(i)
                 + x03(i)*y46(i) + x04(i)*y45(i) + x05(i)*y44(i)
                 + x06(i)*y43(i) + x07(i)*y42(i) + x08(i)*y41(i)
                 + x09(i)*y40(i) + x10(i)*y39(i) + x11(i)*y38(i)
                 + x12(i)*y37(i) + x13(i)*y36(i) + x14(i)*y35(i)
                 + x15(i)*y34(i) + x16(i)*y33(i) + x17(i)*y32(i)
                 + x18(i)*y31(i) + x19(i)*y30(i) + x20(i)*y29(i)
                 + x21(i)*y28(i) + x22(i)*y27(i) + x23(i)*y26(i)
                 + x24(i)*y25(i) + x25(i)*y24(i) + x26(i)*y23(i)
                 + x27(i)*y22(i) + x28(i)*y21(i) + x29(i)*y20(i)
                 + x30(i)*y19(i) + x31(i)*y18(i) + x32(i)*y17(i)
                 + x33(i)*y16(i) + x34(i)*y15(i) + x35(i)*y14(i)
                 + x36(i)*y13(i) + x37(i)*y12(i) + x38(i)*y11(i)
                 + x39(i)*y10(i) + x40(i)*y09(i) + x41(i)*y08(i)
                 + x42(i)*y07(i) + x43(i)*y06(i) + x44(i)*y05(i)
                 + x45(i)*y04(i) + x46(i)*y03(i) + x47(i)*y02(i)
                 + x48(i)*y01(i) + x49(i)*y00(i);
      s50 := s50 + x00(i)*y50(i) + x01(i)*y49(i) + x02(i)*y48(i)
                 + x03(i)*y47(i) + x04(i)*y46(i) + x05(i)*y45(i)
                 + x06(i)*y44(i) + x07(i)*y43(i) + x08(i)*y42(i)
                 + x09(i)*y41(i) + x10(i)*y40(i) + x11(i)*y39(i)
                 + x12(i)*y38(i) + x13(i)*y37(i) + x14(i)*y36(i)
                 + x15(i)*y35(i) + x16(i)*y34(i) + x17(i)*y33(i)
                 + x18(i)*y32(i) + x19(i)*y31(i) + x20(i)*y30(i)
                 + x21(i)*y29(i) + x22(i)*y28(i) + x23(i)*y27(i)
                 + x24(i)*y26(i) + x25(i)*y25(i) + x26(i)*y24(i)
                 + x27(i)*y23(i) + x28(i)*y22(i) + x29(i)*y21(i)
                 + x30(i)*y20(i) + x31(i)*y19(i) + x32(i)*y18(i)
                 + x33(i)*y17(i) + x34(i)*y16(i) + x35(i)*y15(i)
                 + x36(i)*y14(i) + x37(i)*y13(i) + x38(i)*y12(i)
                 + x39(i)*y11(i) + x40(i)*y10(i) + x41(i)*y09(i)
                 + x42(i)*y08(i) + x43(i)*y07(i) + x44(i)*y06(i)
                 + x45(i)*y05(i) + x46(i)*y04(i) + x47(i)*y03(i)
                 + x48(i)*y02(i) + x49(i)*y01(i) + x50(i)*y00(i);
      s51 := s51 + x00(i)*y51(i) + x01(i)*y50(i) + x02(i)*y49(i)
                 + x03(i)*y48(i) + x04(i)*y47(i) + x05(i)*y46(i)
                 + x06(i)*y45(i) + x07(i)*y44(i) + x08(i)*y43(i)
                 + x09(i)*y42(i) + x10(i)*y41(i) + x11(i)*y40(i)
                 + x12(i)*y39(i) + x13(i)*y38(i) + x14(i)*y37(i)
                 + x15(i)*y36(i) + x16(i)*y35(i) + x17(i)*y34(i)
                 + x18(i)*y33(i) + x19(i)*y32(i) + x20(i)*y31(i)
                 + x21(i)*y30(i) + x22(i)*y29(i) + x23(i)*y28(i)
                 + x24(i)*y27(i) + x25(i)*y26(i) + x26(i)*y25(i)
                 + x27(i)*y24(i) + x28(i)*y23(i) + x29(i)*y22(i)
                 + x30(i)*y21(i) + x31(i)*y20(i) + x32(i)*y19(i)
                 + x33(i)*y18(i) + x34(i)*y17(i) + x35(i)*y16(i)
                 + x36(i)*y15(i) + x37(i)*y14(i) + x38(i)*y13(i)
                 + x39(i)*y12(i) + x40(i)*y11(i) + x41(i)*y10(i)
                 + x42(i)*y09(i) + x43(i)*y08(i) + x44(i)*y07(i)
                 + x45(i)*y06(i) + x46(i)*y05(i) + x47(i)*y04(i)
                 + x48(i)*y03(i) + x49(i)*y02(i) + x50(i)*y01(i)
                 + x51(i)*y00(i);
      s52 := s52 + x00(i)*y52(i) + x01(i)*y51(i) + x02(i)*y50(i)
                 + x03(i)*y49(i) + x04(i)*y48(i) + x05(i)*y47(i)
                 + x06(i)*y46(i) + x07(i)*y45(i) + x08(i)*y44(i)
                 + x09(i)*y43(i) + x10(i)*y42(i) + x11(i)*y41(i)
                 + x12(i)*y40(i) + x13(i)*y39(i) + x14(i)*y38(i)
                 + x15(i)*y37(i) + x16(i)*y36(i) + x17(i)*y35(i)
                 + x18(i)*y34(i) + x19(i)*y33(i) + x20(i)*y32(i)
                 + x21(i)*y31(i) + x22(i)*y30(i) + x23(i)*y29(i)
                 + x24(i)*y28(i) + x25(i)*y27(i) + x26(i)*y26(i)
                 + x27(i)*y25(i) + x28(i)*y24(i) + x29(i)*y23(i)
                 + x30(i)*y22(i) + x31(i)*y21(i) + x32(i)*y20(i)
                 + x33(i)*y19(i) + x34(i)*y18(i) + x35(i)*y17(i)
                 + x36(i)*y16(i) + x37(i)*y15(i) + x38(i)*y14(i)
                 + x39(i)*y13(i) + x40(i)*y12(i) + x41(i)*y11(i)
                 + x42(i)*y10(i) + x43(i)*y09(i) + x44(i)*y08(i)
                 + x45(i)*y07(i) + x46(i)*y06(i) + x47(i)*y05(i)
                 + x48(i)*y04(i) + x49(i)*y03(i) + x50(i)*y02(i)
                 + x51(i)*y01(i) + x52(i)*y00(i);
      s53 := s53 + x00(i)*y53(i) + x01(i)*y52(i) + x02(i)*y51(i)
                 + x03(i)*y50(i) + x04(i)*y49(i) + x05(i)*y48(i)
                 + x06(i)*y47(i) + x07(i)*y46(i) + x08(i)*y45(i)
                 + x09(i)*y44(i) + x10(i)*y43(i) + x11(i)*y42(i)
                 + x12(i)*y41(i) + x13(i)*y40(i) + x14(i)*y39(i)
                 + x15(i)*y38(i) + x16(i)*y37(i) + x17(i)*y36(i)
                 + x18(i)*y35(i) + x19(i)*y34(i) + x20(i)*y33(i)
                 + x21(i)*y32(i) + x22(i)*y31(i) + x23(i)*y30(i)
                 + x24(i)*y29(i) + x25(i)*y28(i) + x26(i)*y27(i)
                 + x27(i)*y26(i) + x28(i)*y25(i) + x29(i)*y24(i)
                 + x30(i)*y23(i) + x31(i)*y22(i) + x32(i)*y21(i)
                 + x33(i)*y20(i) + x34(i)*y19(i) + x35(i)*y18(i)
                 + x36(i)*y17(i) + x37(i)*y16(i) + x38(i)*y15(i)
                 + x39(i)*y14(i) + x40(i)*y13(i) + x41(i)*y12(i)
                 + x42(i)*y11(i) + x43(i)*y10(i) + x44(i)*y09(i)
                 + x45(i)*y08(i) + x46(i)*y07(i) + x47(i)*y06(i)
                 + x48(i)*y05(i) + x49(i)*y04(i) + x50(i)*y03(i)
                 + x51(i)*y02(i) + x52(i)*y01(i) + x53(i)*y00(i);
      s54 := s54 + x00(i)*y54(i) + x01(i)*y53(i) + x02(i)*y52(i)
                 + x03(i)*y51(i) + x04(i)*y50(i) + x05(i)*y49(i)
                 + x06(i)*y48(i) + x07(i)*y47(i) + x08(i)*y46(i)
                 + x09(i)*y45(i) + x10(i)*y44(i) + x11(i)*y43(i)
                 + x12(i)*y42(i) + x13(i)*y41(i) + x14(i)*y40(i)
                 + x15(i)*y39(i) + x16(i)*y38(i) + x17(i)*y37(i)
                 + x18(i)*y36(i) + x19(i)*y35(i) + x20(i)*y34(i)
                 + x21(i)*y33(i) + x22(i)*y32(i) + x23(i)*y31(i)
                 + x24(i)*y30(i) + x25(i)*y29(i) + x26(i)*y28(i)
                 + x27(i)*y27(i) + x28(i)*y26(i) + x29(i)*y25(i)
                 + x30(i)*y24(i) + x31(i)*y23(i) + x32(i)*y22(i)
                 + x33(i)*y21(i) + x34(i)*y20(i) + x35(i)*y19(i)
                 + x36(i)*y18(i) + x37(i)*y17(i) + x38(i)*y16(i)
                 + x39(i)*y15(i) + x40(i)*y14(i) + x41(i)*y13(i)
                 + x42(i)*y12(i) + x43(i)*y11(i) + x44(i)*y10(i)
                 + x45(i)*y09(i) + x46(i)*y08(i) + x47(i)*y07(i)
                 + x48(i)*y06(i) + x49(i)*y05(i) + x50(i)*y04(i)
                 + x51(i)*y03(i) + x52(i)*y02(i) + x53(i)*y01(i)
                 + x54(i)*y00(i);
      s55 := s55 + x00(i)*y55(i) + x01(i)*y54(i) + x02(i)*y53(i)
                 + x03(i)*y52(i) + x04(i)*y51(i) + x05(i)*y50(i)
                 + x06(i)*y49(i) + x07(i)*y48(i) + x08(i)*y47(i)
                 + x09(i)*y46(i) + x10(i)*y45(i) + x11(i)*y44(i)
                 + x12(i)*y43(i) + x13(i)*y42(i) + x14(i)*y41(i)
                 + x15(i)*y40(i) + x16(i)*y39(i) + x17(i)*y38(i)
                 + x18(i)*y37(i) + x19(i)*y36(i) + x20(i)*y35(i)
                 + x21(i)*y34(i) + x22(i)*y33(i) + x23(i)*y32(i)
                 + x24(i)*y31(i) + x25(i)*y30(i) + x26(i)*y29(i)
                 + x27(i)*y28(i) + x28(i)*y27(i) + x29(i)*y26(i)
                 + x30(i)*y25(i) + x31(i)*y24(i) + x32(i)*y23(i)
                 + x33(i)*y22(i) + x34(i)*y21(i) + x35(i)*y20(i)
                 + x36(i)*y19(i) + x37(i)*y18(i) + x38(i)*y17(i)
                 + x39(i)*y16(i) + x40(i)*y15(i) + x41(i)*y14(i)
                 + x42(i)*y13(i) + x43(i)*y12(i) + x44(i)*y11(i)
                 + x45(i)*y10(i) + x46(i)*y09(i) + x47(i)*y08(i)
                 + x48(i)*y07(i) + x49(i)*y06(i) + x50(i)*y05(i)
                 + x51(i)*y04(i) + x52(i)*y03(i) + x53(i)*y02(i)
                 + x54(i)*y01(i) + x55(i)*y00(i);
      s56 := s56 + x00(i)*y56(i) + x01(i)*y55(i) + x02(i)*y54(i)
                 + x03(i)*y53(i) + x04(i)*y52(i) + x05(i)*y51(i)
                 + x06(i)*y50(i) + x07(i)*y49(i) + x08(i)*y48(i)
                 + x09(i)*y47(i) + x10(i)*y46(i) + x11(i)*y45(i)
                 + x12(i)*y44(i) + x13(i)*y43(i) + x14(i)*y42(i)
                 + x15(i)*y41(i) + x16(i)*y40(i) + x17(i)*y39(i)
                 + x18(i)*y38(i) + x19(i)*y37(i) + x20(i)*y36(i)
                 + x21(i)*y35(i) + x22(i)*y34(i) + x23(i)*y33(i)
                 + x24(i)*y32(i) + x25(i)*y31(i) + x26(i)*y30(i)
                 + x27(i)*y29(i) + x28(i)*y28(i) + x29(i)*y27(i)
                 + x30(i)*y26(i) + x31(i)*y25(i) + x32(i)*y24(i)
                 + x33(i)*y23(i) + x34(i)*y22(i) + x35(i)*y21(i)
                 + x36(i)*y20(i) + x37(i)*y19(i) + x38(i)*y18(i)
                 + x39(i)*y17(i) + x40(i)*y16(i) + x41(i)*y15(i)
                 + x42(i)*y14(i) + x43(i)*y13(i) + x44(i)*y12(i)
                 + x45(i)*y11(i) + x46(i)*y10(i) + x47(i)*y09(i)
                 + x48(i)*y08(i) + x49(i)*y07(i) + x50(i)*y06(i)
                 + x51(i)*y05(i) + x52(i)*y04(i) + x53(i)*y03(i)
                 + x54(i)*y02(i) + x55(i)*y01(i) + x56(i)*y00(i);
      s57 := s57 + x00(i)*y57(i) + x01(i)*y56(i) + x02(i)*y55(i)
                 + x03(i)*y54(i) + x04(i)*y53(i) + x05(i)*y52(i)
                 + x06(i)*y51(i) + x07(i)*y50(i) + x08(i)*y49(i)
                 + x09(i)*y48(i) + x10(i)*y47(i) + x11(i)*y46(i)
                 + x12(i)*y45(i) + x13(i)*y44(i) + x14(i)*y43(i)
                 + x15(i)*y42(i) + x16(i)*y41(i) + x17(i)*y40(i)
                 + x18(i)*y39(i) + x19(i)*y38(i) + x20(i)*y37(i)
                 + x21(i)*y36(i) + x22(i)*y35(i) + x23(i)*y34(i)
                 + x24(i)*y33(i) + x25(i)*y32(i) + x26(i)*y31(i)
                 + x27(i)*y30(i) + x28(i)*y29(i) + x29(i)*y28(i)
                 + x30(i)*y27(i) + x31(i)*y26(i) + x32(i)*y25(i)
                 + x33(i)*y24(i) + x34(i)*y23(i) + x35(i)*y22(i)
                 + x36(i)*y21(i) + x37(i)*y20(i) + x38(i)*y19(i)
                 + x39(i)*y18(i) + x40(i)*y17(i) + x41(i)*y16(i)
                 + x42(i)*y15(i) + x43(i)*y14(i) + x44(i)*y13(i)
                 + x45(i)*y12(i) + x46(i)*y11(i) + x47(i)*y10(i)
                 + x48(i)*y09(i) + x49(i)*y08(i) + x50(i)*y07(i)
                 + x51(i)*y06(i) + x52(i)*y05(i) + x53(i)*y04(i)
                 + x54(i)*y03(i) + x55(i)*y02(i) + x56(i)*y01(i)
                 + x57(i)*y00(i);
      s58 := s58 + x00(i)*y58(i) + x01(i)*y57(i) + x02(i)*y56(i)
                 + x03(i)*y55(i) + x04(i)*y54(i) + x05(i)*y53(i)
                 + x06(i)*y52(i) + x07(i)*y51(i) + x08(i)*y50(i)
                 + x09(i)*y49(i) + x10(i)*y48(i) + x11(i)*y47(i)
                 + x12(i)*y46(i) + x13(i)*y45(i) + x14(i)*y44(i)
                 + x15(i)*y43(i) + x16(i)*y42(i) + x17(i)*y41(i)
                 + x18(i)*y40(i) + x19(i)*y39(i) + x20(i)*y38(i)
                 + x21(i)*y37(i) + x22(i)*y36(i) + x23(i)*y35(i)
                 + x24(i)*y34(i) + x25(i)*y33(i) + x26(i)*y32(i)
                 + x27(i)*y31(i) + x28(i)*y30(i) + x29(i)*y29(i)
                 + x30(i)*y28(i) + x31(i)*y27(i) + x32(i)*y26(i)
                 + x33(i)*y25(i) + x34(i)*y24(i) + x35(i)*y23(i)
                 + x36(i)*y22(i) + x37(i)*y21(i) + x38(i)*y20(i)
                 + x39(i)*y19(i) + x40(i)*y18(i) + x41(i)*y17(i)
                 + x42(i)*y16(i) + x43(i)*y15(i) + x44(i)*y14(i)
                 + x45(i)*y13(i) + x46(i)*y12(i) + x47(i)*y11(i)
                 + x48(i)*y10(i) + x49(i)*y09(i) + x50(i)*y08(i)
                 + x51(i)*y07(i) + x52(i)*y06(i) + x53(i)*y05(i)
                 + x54(i)*y04(i) + x55(i)*y03(i) + x56(i)*y02(i)
                 + x57(i)*y01(i) + x58(i)*y00(i);
      s59 := s59 + x00(i)*y59(i) + x01(i)*y58(i) + x02(i)*y57(i)
                 + x03(i)*y56(i) + x04(i)*y55(i) + x05(i)*y54(i)
                 + x06(i)*y53(i) + x07(i)*y52(i) + x08(i)*y51(i)
                 + x09(i)*y50(i) + x10(i)*y49(i) + x11(i)*y48(i)
                 + x12(i)*y47(i) + x13(i)*y46(i) + x14(i)*y45(i)
                 + x15(i)*y44(i) + x16(i)*y43(i) + x17(i)*y42(i)
                 + x18(i)*y41(i) + x19(i)*y40(i) + x20(i)*y39(i)
                 + x21(i)*y38(i) + x22(i)*y37(i) + x23(i)*y36(i)
                 + x24(i)*y35(i) + x25(i)*y34(i) + x26(i)*y33(i)
                 + x27(i)*y32(i) + x28(i)*y31(i) + x29(i)*y30(i)
                 + x30(i)*y29(i) + x31(i)*y28(i) + x32(i)*y27(i)
                 + x33(i)*y26(i) + x34(i)*y25(i) + x35(i)*y24(i)
                 + x36(i)*y23(i) + x37(i)*y22(i) + x38(i)*y21(i)
                 + x39(i)*y20(i) + x40(i)*y19(i) + x41(i)*y18(i)
                 + x42(i)*y17(i) + x43(i)*y16(i) + x44(i)*y15(i)
                 + x45(i)*y14(i) + x46(i)*y13(i) + x47(i)*y12(i)
                 + x48(i)*y11(i) + x49(i)*y10(i) + x50(i)*y09(i)
                 + x51(i)*y08(i) + x52(i)*y07(i) + x53(i)*y06(i)
                 + x54(i)*y05(i) + x55(i)*y04(i) + x56(i)*y03(i)
                 + x57(i)*y02(i) + x58(i)*y01(i) + x59(i)*y00(i);
      s60 := s60 + x00(i)*y60(i) + x01(i)*y59(i) + x02(i)*y58(i)
                 + x03(i)*y57(i) + x04(i)*y56(i) + x05(i)*y55(i)
                 + x06(i)*y54(i) + x07(i)*y53(i) + x08(i)*y52(i)
                 + x09(i)*y51(i) + x10(i)*y50(i) + x11(i)*y49(i)
                 + x12(i)*y48(i) + x13(i)*y47(i) + x14(i)*y46(i)
                 + x15(i)*y45(i) + x16(i)*y44(i) + x17(i)*y43(i)
                 + x18(i)*y42(i) + x19(i)*y41(i) + x20(i)*y40(i)
                 + x21(i)*y39(i) + x22(i)*y38(i) + x23(i)*y37(i)
                 + x24(i)*y36(i) + x25(i)*y35(i) + x26(i)*y34(i)
                 + x27(i)*y33(i) + x28(i)*y32(i) + x29(i)*y31(i)
                 + x30(i)*y30(i) + x31(i)*y29(i) + x32(i)*y28(i)
                 + x33(i)*y27(i) + x34(i)*y26(i) + x35(i)*y25(i)
                 + x36(i)*y24(i) + x37(i)*y23(i) + x38(i)*y22(i)
                 + x39(i)*y21(i) + x40(i)*y20(i) + x41(i)*y19(i)
                 + x42(i)*y18(i) + x43(i)*y17(i) + x44(i)*y16(i)
                 + x45(i)*y15(i) + x46(i)*y14(i) + x47(i)*y13(i)
                 + x48(i)*y12(i) + x49(i)*y11(i) + x50(i)*y10(i)
                 + x51(i)*y09(i) + x52(i)*y08(i) + x53(i)*y07(i)
                 + x54(i)*y06(i) + x55(i)*y05(i) + x56(i)*y04(i)
                 + x57(i)*y03(i) + x58(i)*y02(i) + x59(i)*y01(i)
                 + x60(i)*y00(i);
      s61 := s61 + x00(i)*y61(i) + x01(i)*y60(i) + x02(i)*y59(i)
                 + x03(i)*y58(i) + x04(i)*y57(i) + x05(i)*y56(i)
                 + x06(i)*y55(i) + x07(i)*y54(i) + x08(i)*y53(i)
                 + x09(i)*y52(i) + x10(i)*y51(i) + x11(i)*y50(i)
                 + x12(i)*y49(i) + x13(i)*y48(i) + x14(i)*y47(i)
                 + x15(i)*y46(i) + x16(i)*y45(i) + x17(i)*y44(i)
                 + x18(i)*y43(i) + x19(i)*y42(i) + x20(i)*y41(i)
                 + x21(i)*y40(i) + x22(i)*y39(i) + x23(i)*y38(i)
                 + x24(i)*y37(i) + x25(i)*y36(i) + x26(i)*y35(i)
                 + x27(i)*y34(i) + x28(i)*y33(i) + x29(i)*y32(i)
                 + x30(i)*y31(i) + x31(i)*y30(i) + x32(i)*y29(i)
                 + x33(i)*y28(i) + x34(i)*y27(i) + x35(i)*y26(i)
                 + x36(i)*y25(i) + x37(i)*y24(i) + x38(i)*y23(i)
                 + x39(i)*y22(i) + x40(i)*y21(i) + x41(i)*y20(i)
                 + x42(i)*y19(i) + x43(i)*y18(i) + x44(i)*y17(i)
                 + x45(i)*y16(i) + x46(i)*y15(i) + x47(i)*y14(i)
                 + x48(i)*y13(i) + x49(i)*y12(i) + x50(i)*y11(i)
                 + x51(i)*y10(i) + x52(i)*y09(i) + x53(i)*y08(i)
                 + x54(i)*y07(i) + x55(i)*y06(i) + x56(i)*y05(i)
                 + x57(i)*y04(i) + x58(i)*y03(i) + x59(i)*y02(i)
                 + x60(i)*y01(i) + x61(i)*y00(i);
      s62 := s62 + x00(i)*y62(i) + x01(i)*y61(i) + x02(i)*y60(i)
                 + x03(i)*y59(i) + x04(i)*y58(i) + x05(i)*y57(i)
                 + x06(i)*y56(i) + x07(i)*y55(i) + x08(i)*y54(i)
                 + x09(i)*y53(i) + x10(i)*y52(i) + x11(i)*y51(i)
                 + x12(i)*y50(i) + x13(i)*y49(i) + x14(i)*y48(i)
                 + x15(i)*y47(i) + x16(i)*y46(i) + x17(i)*y45(i)
                 + x18(i)*y44(i) + x19(i)*y43(i) + x20(i)*y42(i)
                 + x21(i)*y41(i) + x22(i)*y40(i) + x23(i)*y39(i)
                 + x24(i)*y38(i) + x25(i)*y37(i) + x26(i)*y36(i)
                 + x27(i)*y35(i) + x28(i)*y34(i) + x29(i)*y33(i)
                 + x30(i)*y32(i) + x31(i)*y31(i) + x32(i)*y30(i)
                 + x33(i)*y29(i) + x34(i)*y28(i) + x35(i)*y27(i)
                 + x36(i)*y26(i) + x37(i)*y25(i) + x38(i)*y24(i)
                 + x39(i)*y23(i) + x40(i)*y22(i) + x41(i)*y21(i)
                 + x42(i)*y20(i) + x43(i)*y19(i) + x44(i)*y18(i)
                 + x45(i)*y17(i) + x46(i)*y16(i) + x47(i)*y15(i)
                 + x48(i)*y14(i) + x49(i)*y13(i) + x50(i)*y12(i)
                 + x51(i)*y11(i) + x52(i)*y10(i) + x53(i)*y09(i)
                 + x54(i)*y08(i) + x55(i)*y07(i) + x56(i)*y06(i)
                 + x57(i)*y05(i) + x58(i)*y04(i) + x59(i)*y03(i)
                 + x60(i)*y02(i) + x61(i)*y01(i) + x62(i)*y00(i);
      s63 := s63 + x00(i)*y63(i) + x01(i)*y62(i) + x02(i)*y61(i)
                 + x03(i)*y60(i) + x04(i)*y59(i) + x05(i)*y58(i)
                 + x06(i)*y57(i) + x07(i)*y56(i) + x08(i)*y55(i)
                 + x09(i)*y54(i) + x10(i)*y53(i) + x11(i)*y52(i)
                 + x12(i)*y51(i) + x13(i)*y50(i) + x14(i)*y49(i)
                 + x15(i)*y48(i) + x16(i)*y47(i) + x17(i)*y46(i)
                 + x18(i)*y45(i) + x19(i)*y44(i) + x20(i)*y43(i)
                 + x21(i)*y42(i) + x22(i)*y41(i) + x23(i)*y40(i)
                 + x24(i)*y39(i) + x25(i)*y38(i) + x26(i)*y37(i)
                 + x27(i)*y36(i) + x28(i)*y35(i) + x29(i)*y34(i)
                 + x30(i)*y33(i) + x31(i)*y32(i) + x32(i)*y31(i)
                 + x33(i)*y30(i) + x34(i)*y29(i) + x35(i)*y28(i)
                 + x36(i)*y27(i) + x37(i)*y26(i) + x38(i)*y25(i)
                 + x39(i)*y24(i) + x40(i)*y23(i) + x41(i)*y22(i)
                 + x42(i)*y21(i) + x43(i)*y20(i) + x44(i)*y19(i)
                 + x45(i)*y18(i) + x46(i)*y17(i) + x47(i)*y16(i)
                 + x48(i)*y15(i) + x49(i)*y14(i) + x50(i)*y13(i)
                 + x51(i)*y12(i) + x52(i)*y11(i) + x53(i)*y10(i)
                 + x54(i)*y09(i) + x55(i)*y08(i) + x56(i)*y07(i)
                 + x57(i)*y06(i) + x58(i)*y05(i) + x59(i)*y04(i)
                 + x60(i)*y03(i) + x61(i)*y02(i) + x62(i)*y01(i)
                 + x63(i)*y00(i);
    end loop;
  end Balanced_Quarter_Product;

  procedure Write_Subsums
              ( s00,s01,s02,s03,s04,s05,s06,s07 : in double_float;
                s08,s09,s10,s11,s12,s13,s14,s15 : in double_float;
                s16,s17,s18,s19,s20,s21,s22,s23 : in double_float;
                s24,s25,s26,s27,s28,s29,s30,s31 : in double_float;
                s32,s33,s34,s35,s36,s37,s38,s39 : in double_float;
                s40,s41,s42,s43,s44,s45,s46,s47 : in double_float;
                s48,s49,s50,s51,s52,s53,s54,s55 : in double_float;
                s56,s57,s58,s59,s60,s61,s62,s63 : in double_float ) is

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
    put("s32 : "); put(s32);
    put(", n32 : "); put(Last_Zero_Count(s32),1); new_line;
    put("s33 : "); put(s33);
    put(", n33 : "); put(Last_Zero_Count(s33),1); new_line;
    put("s34 : "); put(s34);
    put(", n34 : "); put(Last_Zero_Count(s34),1); new_line;
    put("s35 : "); put(s35);
    put(", n35 : "); put(Last_Zero_Count(s35),1); new_line;
    put("s36 : "); put(s36);
    put(", n36 : "); put(Last_Zero_Count(s36),1); new_line;
    put("s37 : "); put(s37);
    put(", n37 : "); put(Last_Zero_Count(s37),1); new_line;
    put("s38 : "); put(s38);
    put(", n38 : "); put(Last_Zero_Count(s38),1); new_line;
    put("s39 : "); put(s39);
    put(", n39 : "); put(Last_Zero_Count(s39),1); new_line;
    put("s40 : "); put(s40);
    put(", n40 : "); put(Last_Zero_Count(s40),1); new_line;
    put("s41 : "); put(s41);
    put(", n41 : "); put(Last_Zero_Count(s41),1); new_line;
    put("s42 : "); put(s42);
    put(", n42 : "); put(Last_Zero_Count(s42),1); new_line;
    put("s43 : "); put(s43);
    put(", n43 : "); put(Last_Zero_Count(s43),1); new_line;
    put("s44 : "); put(s44);
    put(", n44 : "); put(Last_Zero_Count(s44),1); new_line;
    put("s45 : "); put(s45);
    put(", n45 : "); put(Last_Zero_Count(s45),1); new_line;
    put("s46 : "); put(s46);
    put(", n46 : "); put(Last_Zero_Count(s46),1); new_line;
    put("s47 : "); put(s47);
    put(", n47 : "); put(Last_Zero_Count(s47),1); new_line;
    put("s48 : "); put(s48);
    put(", n48 : "); put(Last_Zero_Count(s48),1); new_line;
    put("s49 : "); put(s49);
    put(", n49 : "); put(Last_Zero_Count(s49),1); new_line;
    put("s50 : "); put(s50);
    put(", n50 : "); put(Last_Zero_Count(s50),1); new_line;
    put("s51 : "); put(s51);
    put(", n51 : "); put(Last_Zero_Count(s51),1); new_line;
    put("s52 : "); put(s52);
    put(", n52 : "); put(Last_Zero_Count(s52),1); new_line;
    put("s53 : "); put(s53);
    put(", n53 : "); put(Last_Zero_Count(s53),1); new_line;
    put("s54 : "); put(s54);
    put(", n54 : "); put(Last_Zero_Count(s54),1); new_line;
    put("s55 : "); put(s55);
    put(", n55 : "); put(Last_Zero_Count(s55),1); new_line;
    put("s56 : "); put(s56);
    put(", n56 : "); put(Last_Zero_Count(s56),1); new_line;
    put("s57 : "); put(s57);
    put(", n57 : "); put(Last_Zero_Count(s57),1); new_line;
    put("s58 : "); put(s58);
    put(", n58 : "); put(Last_Zero_Count(s58),1); new_line;
    put("s59 : "); put(s59);
    put(", n59 : "); put(Last_Zero_Count(s59),1); new_line;
    put("s60 : "); put(s60);
    put(", n60 : "); put(Last_Zero_Count(s60),1); new_line;
    put("s61 : "); put(s61);
    put(", n61 : "); put(Last_Zero_Count(s61),1); new_line;
    put("s62 : "); put(s62);
    put(", n62 : "); put(Last_Zero_Count(s62),1); new_line;
    put("s63 : "); put(s63);
    put(", n63 : "); put(Last_Zero_Count(s63),1); new_line;
  end Write_Subsums;

  function to_hexa_double
              ( s00,s01,s02,s03,s04,s05,s06,s07 : double_float;
                s08,s09,s10,s11,s12,s13,s14,s15 : double_float;
                s16,s17,s18,s19,s20,s21,s22,s23 : double_float;
                s24,s25,s26,s27,s28,s29,s30,s31 : double_float;
                s32,s33,s34,s35,s36,s37,s38,s39 : double_float;
                s40,s41,s42,s43,s44,s45,s46,s47 : double_float;
                s48,s49,s50,s51,s52,s53,s54,s55 : double_float;
                s56,s57,s58,s59,s60,s61,s62,s63 : double_float;
                verbose : boolean := true ) return hexa_double is

    res : hexa_double;

  begin
    if verbose then
      write_subsums(s00,s01,s02,s03,s04,s05,s06,s07,
                    s08,s09,s10,s11,s12,s13,s14,s15,
                    s16,s17,s18,s19,s20,s21,s22,s23,
                    s24,s25,s26,s27,s28,s29,s30,s31,
                    s32,s33,s34,s35,s36,s37,s38,s39,
                    s40,s41,s42,s43,s44,s45,s46,s47,
                    s48,s49,s50,s51,s52,s53,s54,s55,
                    s56,s57,s58,s59,s60,s61,s62,s63);
    end if;
    res := create(s63);
    res := res + create(s62);
    res := res + create(s61);
    res := res + create(s60);
    res := res + create(s59);
    res := res + create(s58);
    res := res + create(s57);
    res := res + create(s56);
    res := res + create(s55);
    res := res + create(s54);
    res := res + create(s53);
    res := res + create(s52);
    res := res + create(s51);
    res := res + create(s50);
    res := res + create(s49);
    res := res + create(s48);
    res := res + create(s47);
    res := res + create(s46);
    res := res + create(s45);
    res := res + create(s44);
    res := res + create(s43);
    res := res + create(s42);
    res := res + create(s41);
    res := res + create(s40);
    res := res + create(s39);
    res := res + create(s38);
    res := res + create(s37);
    res := res + create(s36);
    res := res + create(s35);
    res := res + create(s34);
    res := res + create(s33);
    res := res + create(s32);
    res := res + create(s31);
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
  end to_hexa_double;

end Vectored_Hexa_Doubles;
