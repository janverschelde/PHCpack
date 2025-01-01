with Standard_Random_Numbers;

package body Balanced_Quarter_Doubles is

  function Thirteen_Random_Bits return integer64 is

    res : integer64 := 1; -- first bit is one
    rnd : integer64;

  begin
    for i in 1..12 loop
      rnd := Standard_Random_Numbers.Random(0,1);
      res := 2*res + rnd;
    end loop;
    return res;
  end Thirteen_Random_Bits;

  function Random_Quarter ( e : in integer32 ) return double_float is

    frc : constant integer64 := Thirteen_Random_Bits;
    mrs : constant double_float := double_float(frc);
    res : constant double_float := double_float'compose(mrs,e);

  begin
    return res;
  end Random_Quarter;

  procedure Random ( x0,x1,x2,x3 : out double_float ) is
  begin
    x0 := Random_Quarter(0);
    x1 := Random_Quarter(-13);
    x2 := Random_Quarter(-26);
    x3 := Random_Quarter(-39);
  end Random;

  procedure Random ( x0,x1,x2,x3,x4,x5,x6,x7 : out double_float ) is
  begin
    Random(x0,x1,x2,x3);
    x4 := Random_Quarter(-52);
    x5 := Random_Quarter(-65);
    x6 := Random_Quarter(-78);
    x7 := Random_Quarter(-91);
  end Random;

  procedure Random ( x0,x1,x2,x3,x4,x5,x6,x7 : out double_float;
                     x8,x9,xA,xB,xC,xD,xE,xF : out double_float ) is
  begin
    Random(x0,x1,x2,x3,x4,x5,x6,x7);
    x8 := Random_Quarter(-104);
    x9 := Random_Quarter(-117);
    xA := Random_Quarter(-130);
    xB := Random_Quarter(-143);
    xC := Random_Quarter(-156);
    xD := Random_Quarter(-169);
    xE := Random_Quarter(-182);
    xF := Random_Quarter(-195);
  end Random;

  procedure Random ( x00,x01,x02,x03,x04,x05,x06,x07 : out double_float;
                     x08,x09,x10,x11,x12,x13,x14,x15 : out double_float;
                     x16,x17,x18,x19,x20,x21,x22,x23 : out double_float;
                     x24,x25,x26,x27,x28,x29,x30,x31 : out double_float ) is
  begin
    Random(x00,x01,x02,x03,x04,x05,x06,x07,x08,x09,x10,x11,x12,x13,x14,x15);
    x16 := Random_Quarter(-208);
    x17 := Random_Quarter(-221);
    x18 := Random_Quarter(-234);
    x19 := Random_Quarter(-247);
    x20 := Random_Quarter(-260);
    x21 := Random_Quarter(-273);
    x22 := Random_Quarter(-286);
    x23 := Random_Quarter(-299);
    x24 := Random_Quarter(-312);
    x25 := Random_Quarter(-325);
    x26 := Random_Quarter(-338);
    x27 := Random_Quarter(-351);
    x28 := Random_Quarter(-364);
    x29 := Random_Quarter(-377);
    x30 := Random_Quarter(-390);
    x31 := Random_Quarter(-403);
  end Random;

  procedure Random ( x00,x01,x02,x03,x04,x05,x06,x07 : out double_float;
                     x08,x09,x10,x11,x12,x13,x14,x15 : out double_float;
                     x16,x17,x18,x19,x20,x21,x22,x23 : out double_float;
                     x24,x25,x26,x27,x28,x29,x30,x31 : out double_float;
                     x32,x33,x34,x35,x36,x37,x38,x39 : out double_float;
                     x40,x41,x42,x43,x44,x45,x46,x47 : out double_float;
                     x48,x49,x50,x51,x52,x53,x54,x55 : out double_float;
                     x56,x57,x58,x59,x60,x61,x62,x63 : out double_float ) is
  begin
    Random(x00,x01,x02,x03,x04,x05,x06,x07,x08,x09,x10,x11,x12,x13,x14,x15,
           x16,x17,x18,x19,x20,x21,x22,x23,x24,x25,x26,x27,x28,x29,x30,x31);
    x32 := Random_Quarter(-416);
    x33 := Random_Quarter(-429);
    x34 := Random_Quarter(-442);
    x35 := Random_Quarter(-455);
    x36 := Random_Quarter(-468);
    x37 := Random_Quarter(-481);
    x38 := Random_Quarter(-494);
    x39 := Random_Quarter(-507);
    x40 := Random_Quarter(-520);
    x41 := Random_Quarter(-533);
    x42 := Random_Quarter(-546);
    x43 := Random_Quarter(-559);
    x44 := Random_Quarter(-572);
    x45 := Random_Quarter(-585);
    x46 := Random_Quarter(-598);
    x47 := Random_Quarter(-611);
    x48 := Random_Quarter(-624);
    x49 := Random_Quarter(-637);
    x50 := Random_Quarter(-650);
    x51 := Random_Quarter(-663);
    x52 := Random_Quarter(-676);
    x53 := Random_Quarter(-689);
    x54 := Random_Quarter(-702);
    x55 := Random_Quarter(-715);
    x56 := Random_Quarter(-728);
    x57 := Random_Quarter(-741);
    x58 := Random_Quarter(-754);
    x59 := Random_Quarter(-767);
    x60 := Random_Quarter(-780);
    x61 := Random_Quarter(-793);
    x62 := Random_Quarter(-806);
    x63 := Random_Quarter(-819);
  end Random;

  procedure Random ( dim : in integer32;
                     x0,x1,x2,x3 : out Standard_Floating_Vectors.Vector ) is
  begin
    for i in 1..dim loop
      Random(x0(i),x1(i),x2(i),x3(i));
    end loop;
  end Random;

  procedure Random ( dim : in integer32;
                     x0,x1,x2,x3 : out Standard_Floating_Vectors.Vector;
                     x4,x5,x6,x7 : out Standard_Floating_Vectors.Vector ) is
  begin
    for i in 1..dim loop
      Random(x0(i),x1(i),x2(i),x3(i),x4(i),x5(i),x6(i),x7(i));
    end loop;
  end Random;

  procedure Random ( dim : in integer32;
                     x0,x1,x2,x3 : out Standard_Floating_Vectors.Vector;
                     x4,x5,x6,x7 : out Standard_Floating_Vectors.Vector;
                     x8,x9,xA,xB : out Standard_Floating_Vectors.Vector;
                     xC,xD,xE,xF : out Standard_Floating_Vectors.Vector ) is
  begin
    for i in 1..dim loop
      Random(x0(i),x1(i),x2(i),x3(i),x4(i),x5(i),x6(i),x7(i),
             x8(i),x9(i),xA(i),xB(i),xC(i),xD(i),xE(i),xF(i));
    end loop;
  end Random;

  procedure Random ( dim : in integer32;
                     x00,x01,x02,x03 : out Standard_Floating_Vectors.Vector;
                     x04,x05,x06,x07 : out Standard_Floating_Vectors.Vector;
                     x08,x09,x10,x11 : out Standard_Floating_Vectors.Vector;
                     x12,x13,x14,x15 : out Standard_Floating_Vectors.Vector;
                     x16,x17,x18,x19 : out Standard_Floating_Vectors.Vector;
                     x20,x21,x22,x23 : out Standard_Floating_Vectors.Vector;
                     x24,x25,x26,x27 : out Standard_Floating_Vectors.Vector;
                     x28,x29,x30,x31 : out Standard_Floating_Vectors.Vector ) is
  begin
    for i in 1..dim loop
      Random(x00(i),x01(i),x02(i),x03(i),x04(i),x05(i),x06(i),x07(i),
             x08(i),x09(i),x10(i),x11(i),x12(i),x13(i),x14(i),x15(i),
             x16(i),x17(i),x18(i),x19(i),x20(i),x21(i),x22(i),x23(i),
             x24(i),x25(i),x26(i),x27(i),x28(i),x29(i),x30(i),x31(i));
    end loop;
  end Random;

  procedure Random ( dim : in integer32;
                     x00,x01,x02,x03 : out Standard_Floating_Vectors.Vector;
                     x04,x05,x06,x07 : out Standard_Floating_Vectors.Vector;
                     x08,x09,x10,x11 : out Standard_Floating_Vectors.Vector;
                     x12,x13,x14,x15 : out Standard_Floating_Vectors.Vector;
                     x16,x17,x18,x19 : out Standard_Floating_Vectors.Vector;
                     x20,x21,x22,x23 : out Standard_Floating_Vectors.Vector;
                     x24,x25,x26,x27 : out Standard_Floating_Vectors.Vector;
                     x28,x29,x30,x31 : out Standard_Floating_Vectors.Vector;
                     x32,x33,x34,x35 : out Standard_Floating_Vectors.Vector;
                     x36,x37,x38,x39 : out Standard_Floating_Vectors.Vector;
                     x40,x41,x42,x43 : out Standard_Floating_Vectors.Vector;
                     x44,x45,x46,x47 : out Standard_Floating_Vectors.Vector;
                     x48,x49,x50,x51 : out Standard_Floating_Vectors.Vector;
                     x52,x53,x54,x55 : out Standard_Floating_Vectors.Vector;
                     x56,x57,x58,x59 : out Standard_Floating_Vectors.Vector;
                     x60,x61,x62,x63 : out Standard_Floating_Vectors.Vector ) is
  begin
    for i in 1..dim loop
      Random(x00(i),x01(i),x02(i),x03(i),x04(i),x05(i),x06(i),x07(i),
             x08(i),x09(i),x10(i),x11(i),x12(i),x13(i),x14(i),x15(i),
             x16(i),x17(i),x18(i),x19(i),x20(i),x21(i),x22(i),x23(i),
             x24(i),x25(i),x26(i),x27(i),x28(i),x29(i),x30(i),x31(i),
             x32(i),x33(i),x34(i),x35(i),x36(i),x37(i),x38(i),x39(i),
             x40(i),x41(i),x42(i),x43(i),x44(i),x45(i),x46(i),x47(i),
             x48(i),x49(i),x50(i),x51(i),x52(i),x53(i),x54(i),x55(i),
             x56(i),x57(i),x58(i),x59(i),x60(i),x61(i),x62(i),x63(i));
    end loop;
  end Random;

  function Make_Double ( r0,r1,r2,r3 : double_float ) return double_float is
  begin
    return ((r3 + r2) + r1) + r0;
  end Make_Double;

  function Make_Double_Double
             ( r0,r1,r2,r3,r4,r5,r6,r7 : double_float )
             return double_double is

    res,reshi,reslo : double_double;

  begin
    reshi := ((Double_Double_Numbers.create(r3) + r2) + r1) + r0;
    reslo := ((Double_Double_Numbers.create(r7) + r6) + r5) + r4;
    res := reshi + reslo;
    return res;
  end Make_Double_Double;

  function Make_Quad_Double
             ( r0,r1,r2,r3,r4,r5,r6,r7 : double_float;
               r8,r9,rA,rB,rC,rD,rE,rF : double_float )
             return quad_double is

    res,reshihi,reslohi,reshilo,reslolo : quad_double;

  begin
    reshihi := ((Quad_Double_Numbers.create(r3) + r2) + r1) + r0;
    reslohi := ((Quad_Double_Numbers.create(r7) + r6) + r5) + r4;
    reshilo := ((Quad_Double_Numbers.create(rB) + rA) + r9) + r8;
    reslolo := ((Quad_Double_Numbers.create(rF) + rE) + rD) + rC;
    res := reshihi + reslohi + reshilo + reslolo;
    return res;
  end Make_Quad_Double;

  function Make_Octo_Double
             ( r00,r01,r02,r03,r04,r05,r06,r07 : double_float;
               r08,r09,r10,r11,r12,r13,r14,r15 : double_float;
               r16,r17,r18,r19,r20,r21,r22,r23 : double_float;
               r24,r25,r26,r27,r28,r29,r30,r31 : double_float )
             return octo_double is

    res : octo_double;
    reshihihi,reslohihi,reshilohi,reslolohi : octo_double;
    reshihilo,reslohilo,reshilolo,reslololo : octo_double;

  begin
    reshihihi := ((Octo_Double_Numbers.create(r03) + r02) + r01) + r00;
    reslohihi := ((Octo_Double_Numbers.create(r07) + r06) + r05) + r04;
    reshilohi := ((Octo_Double_Numbers.create(r11) + r10) + r09) + r08;
    reslolohi := ((Octo_Double_Numbers.create(r15) + r14) + r13) + r12;
    reshihilo := ((Octo_Double_Numbers.create(r19) + r18) + r17) + r16;
    reslohilo := ((Octo_Double_Numbers.create(r23) + r22) + r21) + r20;
    reshilolo := ((Octo_Double_Numbers.create(r27) + r26) + r25) + r24;
    reslololo := ((Octo_Double_Numbers.create(r31) + r30) + r29) + r28;
    res := reshihihi + reslohihi + reshilohi + reslolohi
         + reshihilo + reslohilo + reshilolo + reslololo;
    return res;
  end Make_Octo_Double;

  function Make_Hexa_Double
             ( r00,r01,r02,r03,r04,r05,r06,r07 : double_float;
               r08,r09,r10,r11,r12,r13,r14,r15 : double_float;
               r16,r17,r18,r19,r20,r21,r22,r23 : double_float;
               r24,r25,r26,r27,r28,r29,r30,r31 : double_float;
               r32,r33,r34,r35,r36,r37,r38,r39 : double_float;
               r40,r41,r42,r43,r44,r45,r46,r47 : double_float;
               r48,r49,r50,r51,r52,r53,r54,r55 : double_float;
               r56,r57,r58,r59,r60,r61,r62,r63 : double_float )
             return hexa_double is

    res : hexa_double;
    reshihihihi,reslohihihi,reshilohihi,reslolohihi : hexa_double;
    reshihilohi,reslohilohi,reshilolohi,reslololohi : hexa_double;
    reshihihilo,reslohihilo,reshilohilo,reslolohilo : hexa_double;
    reshihilolo,reslohilolo,reshilololo,reslolololo : hexa_double;

  begin
    reshihihihi := ((Hexa_Double_Numbers.create(r03) + r02) + r01) + r00;
    reslohihihi := ((Hexa_Double_Numbers.create(r07) + r06) + r05) + r04;
    reshilohihi := ((Hexa_Double_Numbers.create(r11) + r10) + r09) + r08;
    reslolohihi := ((Hexa_Double_Numbers.create(r15) + r14) + r13) + r12;
    reshihilohi := ((Hexa_Double_Numbers.create(r19) + r18) + r17) + r16;
    reslohilohi := ((Hexa_Double_Numbers.create(r23) + r22) + r21) + r20;
    reshilolohi := ((Hexa_Double_Numbers.create(r27) + r26) + r25) + r24;
    reslololohi := ((Hexa_Double_Numbers.create(r31) + r30) + r29) + r28;
    reshihihilo := ((Hexa_Double_Numbers.create(r35) + r34) + r33) + r32;
    reslohihilo := ((Hexa_Double_Numbers.create(r39) + r38) + r37) + r36;
    reshilohilo := ((Hexa_Double_Numbers.create(r43) + r42) + r41) + r40;
    reslolohilo := ((Hexa_Double_Numbers.create(r47) + r46) + r45) + r44;
    reshihilolo := ((Hexa_Double_Numbers.create(r51) + r50) + r49) + r48;
    reslohilolo := ((Hexa_Double_Numbers.create(r55) + r54) + r53) + r52;
    reshilololo := ((Hexa_Double_Numbers.create(r59) + r58) + r57) + r56;
    reslolololo := ((Hexa_Double_Numbers.create(r63) + r62) + r61) + r60;
    res := reshihihihi + reslohihihi + reshilohihi + reslolohihi
         + reshihilohi + reslohilohi + reshilolohi + reslololohi
         + reshihihilo + reslohihilo + reshilohilo + reslolohilo
         + reshihilolo + reslohilolo + reshilololo + reslolololo;
    return res;
  end Make_Hexa_Double;

  function Make_Doubles
             ( x0,x1,x2,x3 : Standard_Floating_Vectors.Vector )
             return Standard_Floating_Vectors.Vector is

    res : Standard_Floating_Vectors.Vector(x0'range);

  begin
    for i in x0'range loop
      res(i) := Make_Double(x0(i),x1(i),x2(i),x3(i));
    end loop;
    return res;
  end Make_Doubles;

  function Make_Double_Doubles
             ( x0,x1,x2,x3 : Standard_Floating_Vectors.Vector;
               x4,x5,x6,x7 : Standard_Floating_Vectors.Vector )
             return Double_Double_Vectors.Vector is

    res : Double_Double_Vectors.Vector(x0'range);

  begin
    for i in x0'range loop
      res(i) := Make_Double_Double
                  (x0(i),x1(i),x2(i),x3(i),x4(i),x5(i),x6(i),x7(i));
    end loop;
    return res;
  end Make_Double_Doubles;

  function Make_Quad_Doubles
             ( x0,x1,x2,x3 : Standard_Floating_Vectors.Vector;
               x4,x5,x6,x7 : Standard_Floating_Vectors.Vector;
               x8,x9,xA,xB : Standard_Floating_Vectors.Vector;
               xC,xD,xE,xF : Standard_Floating_Vectors.Vector )
             return Quad_Double_Vectors.Vector is

    res : Quad_Double_Vectors.Vector(x0'range);

  begin
    for i in x0'range loop
      res(i) := Make_Quad_Double
                  (x0(i),x1(i),x2(i),x3(i),x4(i),x5(i),x6(i),x7(i),
                   x8(i),x9(i),xA(i),xB(i),xC(i),xD(i),xE(i),xF(i));
    end loop;
    return res;
  end Make_Quad_Doubles;

  function Make_Octo_Doubles
             ( x00,x01,x02,x03 : Standard_Floating_Vectors.Vector;
               x04,x05,x06,x07 : Standard_Floating_Vectors.Vector;
               x08,x09,x10,x11 : Standard_Floating_Vectors.Vector;
               x12,x13,x14,x15 : Standard_Floating_Vectors.Vector;
               x16,x17,x18,x19 : Standard_Floating_Vectors.Vector;
               x20,x21,x22,x23 : Standard_Floating_Vectors.Vector;
               x24,x25,x26,x27 : Standard_Floating_Vectors.Vector;
               x28,x29,x30,x31 : Standard_Floating_Vectors.Vector )
             return Octo_Double_Vectors.Vector is

    res : Octo_Double_Vectors.Vector(x00'range);

  begin
    for i in x00'range loop
      res(i) := Make_Octo_Double
                  (x00(i),x01(i),x02(i),x03(i),x04(i),x05(i),x06(i),x07(i),
                   x08(i),x09(i),x10(i),x11(i),x12(i),x13(i),x14(i),x15(i),
                   x16(i),x17(i),x18(i),x19(i),x20(i),x21(i),x22(i),x23(i),
                   x24(i),x25(i),x26(i),x27(i),x28(i),x29(i),x30(i),x31(i));
    end loop;
    return res;
  end Make_Octo_Doubles;

  function Make_Hexa_Doubles
             ( x00,x01,x02,x03 : Standard_Floating_Vectors.Vector;
               x04,x05,x06,x07 : Standard_Floating_Vectors.Vector;
               x08,x09,x10,x11 : Standard_Floating_Vectors.Vector;
               x12,x13,x14,x15 : Standard_Floating_Vectors.Vector;
               x16,x17,x18,x19 : Standard_Floating_Vectors.Vector;
               x20,x21,x22,x23 : Standard_Floating_Vectors.Vector;
               x24,x25,x26,x27 : Standard_Floating_Vectors.Vector;
               x28,x29,x30,x31 : Standard_Floating_Vectors.Vector;
               x32,x33,x34,x35 : Standard_Floating_Vectors.Vector;
               x36,x37,x38,x39 : Standard_Floating_Vectors.Vector;
               x40,x41,x42,x43 : Standard_Floating_Vectors.Vector;
               x44,x45,x46,x47 : Standard_Floating_Vectors.Vector;
               x48,x49,x50,x51 : Standard_Floating_Vectors.Vector;
               x52,x53,x54,x55 : Standard_Floating_Vectors.Vector;
               x56,x57,x58,x59 : Standard_Floating_Vectors.Vector;
               x60,x61,x62,x63 : Standard_Floating_Vectors.Vector )
             return Hexa_Double_Vectors.Vector is

    res : Hexa_Double_Vectors.Vector(x00'range);

  begin
    for i in x00'range loop
      res(i) := Make_Hexa_Double
                  (x00(i),x01(i),x02(i),x03(i),x04(i),x05(i),x06(i),x07(i),
                   x08(i),x09(i),x10(i),x11(i),x12(i),x13(i),x14(i),x15(i),
                   x16(i),x17(i),x18(i),x19(i),x20(i),x21(i),x22(i),x23(i),
                   x24(i),x25(i),x26(i),x27(i),x28(i),x29(i),x30(i),x31(i),
                   x32(i),x33(i),x34(i),x35(i),x36(i),x37(i),x38(i),x39(i),
                   x40(i),x41(i),x42(i),x43(i),x44(i),x45(i),x46(i),x47(i),
                   x48(i),x49(i),x50(i),x51(i),x52(i),x53(i),x54(i),x55(i),
                   x56(i),x57(i),x58(i),x59(i),x60(i),x61(i),x62(i),x63(i));
    end loop;
    return res;
  end Make_Hexa_Doubles;

-- WRAPPERS :

  function Random return double_float is

    r0,r1,r2,r3 : double_float;

  begin
    Random(r0,r1,r2,r3);
    return Make_Double(r0,r1,r2,r3);
  end Random;

  function Random return double_double is

    r0,r1,r2,r3,r4,r5,r6,r7 : double_float;

  begin
    Random(r0,r1,r2,r3,r4,r5,r6,r7);
    return Make_Double_Double(r0,r1,r2,r3,r4,r5,r6,r7);
  end Random;

  function Random return quad_double is

    r0,r1,r2,r3,r4,r5,r6,r7,r8,r9,rA,rB,rC,rD,rE,rF : double_float;

  begin
    Random(r0,r1,r2,r3,r4,r5,r6,r7,r8,r9,rA,rB,rC,rD,rE,rF);
    return Make_Quad_Double(r0,r1,r2,r3,r4,r5,r6,r7,r8,r9,rA,rB,rC,rD,rE,rF);
  end Random;

  function Random return octo_double is

    r00,r01,r02,r03,r04,r05,r06,r07 : double_float;
    r08,r09,r10,r11,r12,r13,r14,r15 : double_float;
    r16,r17,r18,r19,r20,r21,r22,r23 : double_float;
    r24,r25,r26,r27,r28,r29,r30,r31 : double_float;

  begin
    Random(r00,r01,r02,r03,r04,r05,r06,r07,r08,r09,r10,r11,r12,r13,r14,r15,
           r16,r17,r18,r19,r20,r21,r22,r23,r24,r25,r26,r27,r28,r29,r30,r31);
    return Make_Octo_Double
             (r00,r01,r02,r03,r04,r05,r06,r07,r08,r09,r10,r11,r12,r13,r14,r15,
              r16,r17,r18,r19,r20,r21,r22,r23,r24,r25,r26,r27,r28,r29,r30,r31);
  end Random;

  function Random return hexa_double is

    r00,r01,r02,r03,r04,r05,r06,r07 : double_float;
    r08,r09,r10,r11,r12,r13,r14,r15 : double_float;
    r16,r17,r18,r19,r20,r21,r22,r23 : double_float;
    r24,r25,r26,r27,r28,r29,r30,r31 : double_float;
    r32,r33,r34,r35,r36,r37,r38,r39 : double_float;
    r40,r41,r42,r43,r44,r45,r46,r47 : double_float;
    r48,r49,r50,r51,r52,r53,r54,r55 : double_float;
    r56,r57,r58,r59,r60,r61,r62,r63 : double_float;

  begin
    Random(r00,r01,r02,r03,r04,r05,r06,r07,r08,r09,r10,r11,r12,r13,r14,r15,
           r16,r17,r18,r19,r20,r21,r22,r23,r24,r25,r26,r27,r28,r29,r30,r31,
           r32,r33,r34,r35,r36,r37,r38,r39,r40,r41,r42,r43,r44,r45,r46,r47,
           r48,r49,r50,r51,r52,r53,r54,r55,r56,r57,r58,r59,r60,r61,r62,r63);
    return Make_Hexa_Double
             (r00,r01,r02,r03,r04,r05,r06,r07,r08,r09,r10,r11,r12,r13,r14,r15,
              r16,r17,r18,r19,r20,r21,r22,r23,r24,r25,r26,r27,r28,r29,r30,r31,
              r32,r33,r34,r35,r36,r37,r38,r39,r40,r41,r42,r43,r44,r45,r46,r47,
              r48,r49,r50,r51,r52,r53,r54,r55,r56,r57,r58,r59,r60,r61,r62,r63);
  end Random;

  function Random 
             ( dim : integer32 ) return Standard_Floating_Vectors.Vector is

    res : Standard_Floating_Vectors.Vector(1..dim);

  begin
    for i in 1..dim loop
      res(i) := Random;
    end loop;
    return res;
  end Random;

  function Random ( dim : integer32 ) return Double_Double_Vectors.Vector is

    res : Double_Double_Vectors.Vector(1..dim);

  begin
    for i in 1..dim loop
      res(i) := Random;
    end loop;
    return res;
  end Random;

  function Random ( dim : integer32 ) return Quad_Double_Vectors.Vector is

    res : Quad_Double_Vectors.Vector(1..dim);

  begin
    for i in 1..dim loop
      res(i) := Random;
    end loop;
    return res;
  end Random;

  function Random ( dim : integer32 ) return Octo_Double_Vectors.Vector is

    res : Octo_Double_Vectors.Vector(1..dim);

  begin
    for i in 1..dim loop
      res(i) := Random;
    end loop;
    return res;
  end Random;

  function Random ( dim : integer32 ) return Hexa_Double_Vectors.Vector is

    res : Hexa_Double_Vectors.Vector(1..dim);

  begin
    for i in 1..dim loop
      res(i) := Random;
    end loop;
    return res;
  end Random;

end Balanced_Quarter_Doubles;
