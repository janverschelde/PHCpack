with text_io;                            use text_io;
with Timing_Package;                     use Timing_Package;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Random_Numbers;
with Standard_Floating_Vectors;
with Hexa_Double_Numbers;                use Hexa_Double_Numbers;
with Hexa_Double_Numbers_io;             use Hexa_Double_Numbers_io;
with Hexa_Double_Vectors;
with Balanced_Quarter_Doubles;
with Vectored_Hexa_Doubles;

package body Test_Vectored_Hexa_Doubles is

  procedure Test_Balanced_Product ( dim,freq : in integer32 ) is

    x00,x01,x02,x03 : Standard_Floating_Vectors.Vector(1..dim);
    x04,x05,x06,x07 : Standard_Floating_Vectors.Vector(1..dim);
    x08,x09,x10,x11 : Standard_Floating_Vectors.Vector(1..dim);
    x12,x13,x14,x15 : Standard_Floating_Vectors.Vector(1..dim);
    x16,x17,x18,x19 : Standard_Floating_Vectors.Vector(1..dim);
    x20,x21,x22,x23 : Standard_Floating_Vectors.Vector(1..dim);
    x24,x25,x26,x27 : Standard_Floating_Vectors.Vector(1..dim);
    x28,x29,x30,x31 : Standard_Floating_Vectors.Vector(1..dim);
    x32,x33,x34,x35 : Standard_Floating_Vectors.Vector(1..dim);
    x36,x37,x38,x39 : Standard_Floating_Vectors.Vector(1..dim);
    x40,x41,x42,x43 : Standard_Floating_Vectors.Vector(1..dim);
    x44,x45,x46,x47 : Standard_Floating_Vectors.Vector(1..dim);
    x48,x49,x50,x51 : Standard_Floating_Vectors.Vector(1..dim);
    x52,x53,x54,x55 : Standard_Floating_Vectors.Vector(1..dim);
    x56,x57,x58,x59 : Standard_Floating_Vectors.Vector(1..dim);
    x60,x61,x62,x63 : Standard_Floating_Vectors.Vector(1..dim);
    y00,y01,y02,y03 : Standard_Floating_Vectors.Vector(1..dim);
    y04,y05,y06,y07 : Standard_Floating_Vectors.Vector(1..dim);
    y08,y09,y10,y11 : Standard_Floating_Vectors.Vector(1..dim);
    y12,y13,y14,y15 : Standard_Floating_Vectors.Vector(1..dim);
    y16,y17,y18,y19 : Standard_Floating_Vectors.Vector(1..dim);
    y20,y21,y22,y23 : Standard_Floating_Vectors.Vector(1..dim);
    y24,y25,y26,y27 : Standard_Floating_Vectors.Vector(1..dim);
    y28,y29,y30,y31 : Standard_Floating_Vectors.Vector(1..dim);
    y32,y33,y34,y35 : Standard_Floating_Vectors.Vector(1..dim);
    y36,y37,y38,y39 : Standard_Floating_Vectors.Vector(1..dim);
    y40,y41,y42,y43 : Standard_Floating_Vectors.Vector(1..dim);
    y44,y45,y46,y47 : Standard_Floating_Vectors.Vector(1..dim);
    y48,y49,y50,y51 : Standard_Floating_Vectors.Vector(1..dim);
    y52,y53,y54,y55 : Standard_Floating_Vectors.Vector(1..dim);
    y56,y57,y58,y59 : Standard_Floating_Vectors.Vector(1..dim);
    y60,y61,y62,y63 : Standard_Floating_Vectors.Vector(1..dim);
    s00,s01,s02,s03,s04,s05,s06,s07 : double_float;
    s08,s09,s10,s11,s12,s13,s14,s15 : double_float;
    s16,s17,s18,s19,s20,s21,s22,s23 : double_float;
    s24,s25,s26,s27,s28,s29,s30,s31 : double_float;
    s32,s33,s34,s35,s36,s37,s38,s39 : double_float;
    s40,s41,s42,s43,s44,s45,s46,s47 : double_float;
    s48,s49,s50,s51,s52,s53,s54,s55 : double_float;
    s56,s57,s58,s59,s60,s61,s62,s63 : double_float;
    x,y : Hexa_Double_Vectors.Vector(1..dim);
    hdprd0,hdprd1,err : hexa_double;
    timer0,timer1 : Timing_Widget;

  begin
    put("Testing the balanced product on dimension ");
    put(dim,1); put(" and frequency "); put(freq,1); put_line(" ...");
    Balanced_Quarter_Doubles.Random
      (dim,x00,x01,x02,x03,x04,x05,x06,x07,x08,x09,x10,x11,x12,x13,x14,x15,
           x16,x17,x18,x19,x20,x21,x22,x23,x24,x25,x26,x27,x28,x29,x30,x31,
           x32,x33,x34,x35,x36,x37,x38,x39,x40,x41,x42,x43,x44,x45,x46,x47,
           x48,x49,x50,x51,x52,x53,x54,x55,x56,x57,x58,x59,x60,x61,x62,x63);
    Balanced_Quarter_Doubles.Random
      (dim,y00,y01,y02,y03,y04,y05,y06,y07,y08,y09,y10,y11,y12,y13,y14,y15,
           y16,y17,y18,y19,y20,y21,y22,y23,y24,y25,y26,y27,y28,y29,y30,y31,
           y32,y33,y34,y35,y36,y37,y38,y39,y40,y41,y42,y43,y44,y45,y46,y47,
           y48,y49,y50,y51,y52,y53,y54,y55,y56,y57,y58,y59,y60,y61,y62,y63);
    x := Balanced_Quarter_Doubles.Make_Hexa_Doubles
           (x00,x01,x02,x03,x04,x05,x06,x07,x08,x09,x10,x11,x12,x13,x14,x15,
            x16,x17,x18,x19,x20,x21,x22,x23,x24,x25,x26,x27,x28,x29,x30,x31,
            x32,x33,x34,x35,x36,x37,x38,x39,x40,x41,x42,x43,x44,x45,x46,x47,
            x48,x49,x50,x51,x52,x53,x54,x55,x56,x57,x58,x59,x60,x61,x62,x63);
    y := Balanced_Quarter_Doubles.Make_Hexa_Doubles
           (y00,y01,y02,y03,y04,y05,y06,y07,y08,y09,y10,y11,y12,y13,y14,y15,
            y16,y17,y18,y19,y20,y21,y22,y23,y24,y25,y26,y27,y28,y29,y30,y31,
            y32,y33,y34,y35,y36,y37,y38,y39,y40,y41,y42,y43,y44,y45,y46,y47,
            y48,y49,y50,y51,y52,y53,y54,y55,y56,y57,y58,y59,y60,y61,y62,y63);
    tstart(timer0);
    for i in 1..freq loop
      hdprd0 := create(integer32(0));
      for i in x'range loop
        hdprd0 := hdprd0 + x(i)*y(i);
      end loop;
    end loop;
    tstop(timer0);
    tstart(timer1);
    for i in 1..freq loop
      Vectored_Hexa_Doubles.Balanced_Quarter_Product
        (dim,x00,x01,x02,x03,x04,x05,x06,x07,x08,x09,x10,x11,x12,x13,x14,x15,
             x16,x17,x18,x19,x20,x21,x22,x23,x24,x25,x26,x27,x28,x29,x30,x31,
             x32,x33,x34,x35,x36,x37,x38,x39,x40,x41,x42,x43,x44,x45,x46,x47,
             x48,x49,x50,x51,x52,x53,x54,x55,x56,x57,x58,x59,x60,x61,x62,x63,
             y00,y01,y02,y03,y04,y05,y06,y07,y08,y09,y10,y11,y12,y13,y14,y15,
             y16,y17,y18,y19,y20,y21,y22,y23,y24,y25,y26,y27,y28,y29,y30,y31,
             y32,y33,y34,y35,y36,y37,y38,y39,y40,y41,y42,y43,y44,y45,y46,y47,
             y48,y49,y50,y51,y52,y53,y54,y55,y56,y57,y58,y59,y60,y61,y62,y63,
             s00,s01,s02,s03,s04,s05,s06,s07,s08,s09,s10,s11,s12,s13,s14,s15,
             s16,s17,s18,s19,s20,s21,s22,s23,s24,s25,s26,s27,s28,s29,s30,s31,
             s32,s33,s34,s35,s36,s37,s38,s39,s40,s41,s42,s43,s44,s45,s46,s47,
             s48,s49,s50,s51,s52,s53,s54,s55,s56,s57,s58,s59,s60,s61,s62,s63);
      if freq = 1 then
        hdprd1 := Vectored_Hexa_Doubles.to_hexa_double
          (s00,s01,s02,s03,s04,s05,s06,s07,s08,s09,s10,s11,s12,s13,s14,s15,
           s16,s17,s18,s19,s20,s21,s22,s23,s24,s25,s26,s27,s28,s29,s30,s31,
           s32,s33,s34,s35,s36,s37,s38,s39,s40,s41,s42,s43,s44,s45,s46,s47,
           s48,s49,s50,s51,s52,s53,s54,s55,s56,s57,s58,s59,s60,s61,s62,s63);
      else
        hdprd1 := Vectored_Hexa_Doubles.to_hexa_double
          (s00,s01,s02,s03,s04,s05,s06,s07,s08,s09,s10,s11,s12,s13,s14,s15,
           s16,s17,s18,s19,s20,s21,s22,s23,s24,s25,s26,s27,s28,s29,s30,s31,
           s32,s33,s34,s35,s36,s37,s38,s39,s40,s41,s42,s43,s44,s45,s46,s47,
           s48,s49,s50,s51,s52,s53,s54,s55,s56,s57,s58,s59,s60,s61,s62,s63,
           verbose=>false);
      end if;
    end loop;
    tstop(timer1);
    new_line;
    put("hd prd : "); put(hdprd0); new_line;
    put("hd sgn : "); put(hdprd1); new_line;
    err := hdprd0 - hdprd1;
    put(" error : "); put(err,2); new_line;
    new_line;
    print_times(standard_output,timer0,"hexa double inner product");
    new_line;
    print_times(standard_output,timer1,"vectored hexa double product");
  end Test_Balanced_Product;

  procedure Wall_Time_Test is

    dim : constant integer32 := 6144;
    freq : constant integer32 := 1024;
    x00,x01,x02,x03 : Standard_Floating_Vectors.Vector(1..dim);
    x04,x05,x06,x07 : Standard_Floating_Vectors.Vector(1..dim);
    x08,x09,x10,x11 : Standard_Floating_Vectors.Vector(1..dim);
    x12,x13,x14,x15 : Standard_Floating_Vectors.Vector(1..dim);
    x16,x17,x18,x19 : Standard_Floating_Vectors.Vector(1..dim);
    x20,x21,x22,x23 : Standard_Floating_Vectors.Vector(1..dim);
    x24,x25,x26,x27 : Standard_Floating_Vectors.Vector(1..dim);
    x28,x29,x30,x31 : Standard_Floating_Vectors.Vector(1..dim);
    x32,x33,x34,x35 : Standard_Floating_Vectors.Vector(1..dim);
    x36,x37,x38,x39 : Standard_Floating_Vectors.Vector(1..dim);
    x40,x41,x42,x43 : Standard_Floating_Vectors.Vector(1..dim);
    x44,x45,x46,x47 : Standard_Floating_Vectors.Vector(1..dim);
    x48,x49,x50,x51 : Standard_Floating_Vectors.Vector(1..dim);
    x52,x53,x54,x55 : Standard_Floating_Vectors.Vector(1..dim);
    x56,x57,x58,x59 : Standard_Floating_Vectors.Vector(1..dim);
    x60,x61,x62,x63 : Standard_Floating_Vectors.Vector(1..dim);
    y00,y01,y02,y03 : Standard_Floating_Vectors.Vector(1..dim);
    y04,y05,y06,y07 : Standard_Floating_Vectors.Vector(1..dim);
    y08,y09,y10,y11 : Standard_Floating_Vectors.Vector(1..dim);
    y12,y13,y14,y15 : Standard_Floating_Vectors.Vector(1..dim);
    y16,y17,y18,y19 : Standard_Floating_Vectors.Vector(1..dim);
    y20,y21,y22,y23 : Standard_Floating_Vectors.Vector(1..dim);
    y24,y25,y26,y27 : Standard_Floating_Vectors.Vector(1..dim);
    y28,y29,y30,y31 : Standard_Floating_Vectors.Vector(1..dim);
    y32,y33,y34,y35 : Standard_Floating_Vectors.Vector(1..dim);
    y36,y37,y38,y39 : Standard_Floating_Vectors.Vector(1..dim);
    y40,y41,y42,y43 : Standard_Floating_Vectors.Vector(1..dim);
    y44,y45,y46,y47 : Standard_Floating_Vectors.Vector(1..dim);
    y48,y49,y50,y51 : Standard_Floating_Vectors.Vector(1..dim);
    y52,y53,y54,y55 : Standard_Floating_Vectors.Vector(1..dim);
    y56,y57,y58,y59 : Standard_Floating_Vectors.Vector(1..dim);
    y60,y61,y62,y63 : Standard_Floating_Vectors.Vector(1..dim);
    s00,s01,s02,s03,s04,s05,s06,s07 : double_float;
    s08,s09,s10,s11,s12,s13,s14,s15 : double_float;
    s16,s17,s18,s19,s20,s21,s22,s23 : double_float;
    s24,s25,s26,s27,s28,s29,s30,s31 : double_float;
    s32,s33,s34,s35,s36,s37,s38,s39 : double_float;
    s40,s41,s42,s43,s44,s45,s46,s47 : double_float;
    s48,s49,s50,s51,s52,s53,s54,s55 : double_float;
    s56,s57,s58,s59,s60,s61,s62,s63 : double_float;
    x,y : Hexa_Double_Vectors.Vector(1..dim);
    hdprd1 : hexa_double;
    timer : Timing_Widget;

  begin
   -- Test_Balanced_Product(dim,freq);
    put("Running the balanced product on dimension ");
    put(dim,1); put(" and frequency "); put(freq,1); put_line(" ...");
    Balanced_Quarter_Doubles.Random
      (dim,x00,x01,x02,x03,x04,x05,x06,x07,x08,x09,x10,x11,x12,x13,x14,x15,
           x16,x17,x18,x19,x20,x21,x22,x23,x24,x25,x26,x27,x28,x29,x30,x31,
           x32,x33,x34,x35,x36,x37,x38,x39,x40,x41,x42,x43,x44,x45,x46,x47,
           x48,x49,x50,x51,x52,x53,x54,x55,x56,x57,x58,x59,x60,x61,x62,x63);
    Balanced_Quarter_Doubles.Random
      (dim,y00,y01,y02,y03,y04,y05,y06,y07,y08,y09,y10,y11,y12,y13,y14,y15,
           y16,y17,y18,y19,y20,y21,y22,y23,y24,y25,y26,y27,y28,y29,y30,y31,
           y32,y33,y34,y35,y36,y37,y38,y39,y40,y41,y42,y43,y44,y45,y46,y47,
           y48,y49,y50,y51,y52,y53,y54,y55,y56,y57,y58,y59,y60,y61,y62,y63);
    x := Balanced_Quarter_Doubles.Make_Hexa_Doubles
           (x00,x01,x02,x03,x04,x05,x06,x07,x08,x09,x10,x11,x12,x13,x14,x15,
            x16,x17,x18,x19,x20,x21,x22,x23,x24,x25,x26,x27,x28,x29,x30,x31,
            x32,x33,x34,x35,x36,x37,x38,x39,x40,x41,x42,x43,x44,x45,x46,x47,
            x48,x49,x50,x51,x52,x53,x54,x55,x56,x57,x58,x59,x60,x61,x62,x63);
    y := Balanced_Quarter_Doubles.Make_Hexa_Doubles
           (y00,y01,y02,y03,y04,y05,y06,y07,y08,y09,y10,y11,y12,y13,y14,y15,
            y16,y17,y18,y19,y20,y21,y22,y23,y24,y25,y26,y27,y28,y29,y30,y31,
            y32,y33,y34,y35,y36,y37,y38,y39,y40,y41,y42,y43,y44,y45,y46,y47,
            y48,y49,y50,y51,y52,y53,y54,y55,y56,y57,y58,y59,y60,y61,y62,y63);
    tstart(timer);
    for i in 1..freq loop
      Vectored_Hexa_Doubles.Balanced_Quarter_Product
        (dim,x00,x01,x02,x03,x04,x05,x06,x07,x08,x09,x10,x11,x12,x13,x14,x15,
             x16,x17,x18,x19,x20,x21,x22,x23,x24,x25,x26,x27,x28,x29,x30,x31,
             x32,x33,x34,x35,x36,x37,x38,x39,x40,x41,x42,x43,x44,x45,x46,x47,
             x48,x49,x50,x51,x52,x53,x54,x55,x56,x57,x58,x59,x60,x61,x62,x63,
             y00,y01,y02,y03,y04,y05,y06,y07,y08,y09,y10,y11,y12,y13,y14,y15,
             y16,y17,y18,y19,y20,y21,y22,y23,y24,y25,y26,y27,y28,y29,y30,y31,
             y32,y33,y34,y35,y36,y37,y38,y39,y40,y41,y42,y43,y44,y45,y46,y47,
             y48,y49,y50,y51,y52,y53,y54,y55,y56,y57,y58,y59,y60,y61,y62,y63,
             s00,s01,s02,s03,s04,s05,s06,s07,s08,s09,s10,s11,s12,s13,s14,s15,
             s16,s17,s18,s19,s20,s21,s22,s23,s24,s25,s26,s27,s28,s29,s30,s31,
             s32,s33,s34,s35,s36,s37,s38,s39,s40,s41,s42,s43,s44,s45,s46,s47,
             s48,s49,s50,s51,s52,s53,s54,s55,s56,s57,s58,s59,s60,s61,s62,s63);
      hdprd1 := Vectored_Hexa_Doubles.to_hexa_double
        (s00,s01,s02,s03,s04,s05,s06,s07,s08,s09,s10,s11,s12,s13,s14,s15,
         s16,s17,s18,s19,s20,s21,s22,s23,s24,s25,s26,s27,s28,s29,s30,s31,
         s32,s33,s34,s35,s36,s37,s38,s39,s40,s41,s42,s43,s44,s45,s46,s47,
         s48,s49,s50,s51,s52,s53,s54,s55,s56,s57,s58,s59,s60,s61,s62,s63,
         verbose=>false);
    end loop;
    tstop(timer);
    new_line;
    print_times(standard_output,timer,"vectored hexa double product");
  end Wall_Time_Test;

  procedure Main is

    seed : natural32 := 0;
    dim,freq : integer32 := 0;

  begin
    put("Give the seed (0 for none) : "); get(seed);
    if seed /= 0
     then Standard_Random_Numbers.Set_Seed(seed);
    end if;
    put("Give the dimension : "); get(dim);
    put("Give the frequency : "); get(freq);
    Test_Balanced_Product(dim,freq);
    put("Seed used : "); put(Standard_Random_Numbers.Get_Seed,1); new_line;
  end Main;

end Test_Vectored_Hexa_Doubles;
