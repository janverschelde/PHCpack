with text_io;                            use text_io;
with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;
with Double_Double_Basics;
with Bits_of_Doubles;

package body Vectored_Double_Doubles is

  procedure Split ( v : in Double_Double_Vectors.Vector;
                    v0,v1,v2,v3 : out Standard_Floating_Vectors.Vector ) is

    nbr : double_double;
    flt : double_float;

  begin
    for i in v'range loop
      nbr := v(i);
      flt := hi_part(nbr);
      Double_Double_Basics.Split(flt,v0(i),v1(i));
      flt := lo_part(nbr);
      Double_Double_Basics.Split(flt,v2(i),v3(i));
    end loop;
  end Split;

  procedure Split ( v : in DoblDobl_Complex_Vectors.Vector;
                    v0re,v1re : out Standard_Floating_Vectors.Vector;
                    v2re,v3re : out Standard_Floating_Vectors.Vector;
                    v0im,v1im : out Standard_Floating_Vectors.Vector;
                    v2im,v3im : out Standard_Floating_Vectors.Vector ) is

     nbr : double_double;
     flt : double_float;

  begin
    for i in v'range loop
      nbr := DoblDobl_Complex_Numbers.REAL_PART(v(i));
      flt := hi_part(nbr);
      Double_Double_Basics.Split(flt,v0re(i),v1re(i));
      flt := lo_part(nbr);
      Double_Double_Basics.Split(flt,v2re(i),v3re(i));
      nbr := DoblDobl_Complex_Numbers.IMAG_PART(v(i));
      flt := hi_part(nbr);
      Double_Double_Basics.Split(flt,v0im(i),v1im(i));
      flt := lo_part(nbr);
      Double_Double_Basics.Split(flt,v2im(i),v3im(i));
    end loop;
  end Split;

  procedure Quarter ( v : in Double_Double_Vectors.Vector;
                      v0,v1,v2,v3 : out Standard_Floating_Vectors.Vector;
                      v4,v5,v6,v7 : out Standard_Floating_Vectors.Vector ) is

    nbr : double_double;
    flt : double_float;

  begin
    for i in v'range loop
      nbr := v(i);
      flt := hi_part(nbr);
      Bits_of_Doubles.Split(flt,v0(i),v1(i),v2(i),v3(i));
      flt := lo_part(nbr);
      Bits_of_Doubles.Split(flt,v4(i),v5(i),v6(i),v7(i));
    end loop;
  end Quarter;

  procedure Quarter ( v : in DoblDobl_Complex_Vectors.Vector;
                      v0re,v1re : out Standard_Floating_Vectors.Vector;
                      v2re,v3re : out Standard_Floating_Vectors.Vector;
                      v4re,v5re : out Standard_Floating_Vectors.Vector;
                      v6re,v7re : out Standard_Floating_Vectors.Vector;
                      v0im,v1im : out Standard_Floating_Vectors.Vector;
                      v2im,v3im : out Standard_Floating_Vectors.Vector;
                      v4im,v5im : out Standard_Floating_Vectors.Vector;
                      v6im,v7im : out Standard_Floating_Vectors.Vector ) is

    nbr : double_double;
    flt : double_float;

  begin
    for i in v'range loop
      nbr := DoblDobl_Complex_Numbers.REAL_PART(v(i));
      flt := hi_part(nbr);
      Bits_of_Doubles.Split(flt,v0re(i),v1re(i),v2re(i),v3re(i));
      flt := lo_part(nbr);
      Bits_of_Doubles.Split(flt,v4re(i),v5re(i),v6re(i),v7re(i));
      nbr := DoblDobl_Complex_Numbers.IMAG_PART(v(i));
      flt := hi_part(nbr);
      Bits_of_Doubles.Split(flt,v0im(i),v1im(i),v2im(i),v3im(i));
      flt := lo_part(nbr);
      Bits_of_Doubles.Split(flt,v4im(i),v5im(i),v6im(i),v7im(i));
    end loop;
  end Quarter;

  procedure Sum ( v0,v1,v2,v3 : in Standard_Floating_Vectors.Vector;
                  s0,s1,s2,s3 : out double_float ) is
  begin
    s0 := 0.0;
    s1 := 0.0;
    s2 := 0.0;
    s3 := 0.0;
    for i in v0'range loop
      s0 := s0 + v0(i);
      s1 := s1 + v1(i);
      s2 := s2 + v2(i);
      s3 := s3 + v3(i);
    end loop;
  end Sum;

  procedure Sum ( v0re,v1re,v2re,v3re : in Standard_Floating_Vectors.Vector;
                  v0im,v1im,v2im,v3im : in Standard_Floating_Vectors.Vector;
                  s0re,s1re,s2re,s3re : out double_float;
                  s0im,s1im,s2im,s3im : out double_float ) is
  begin
    s0re := 0.0;
    s1re := 0.0;
    s2re := 0.0;
    s3re := 0.0;
    s0im := 0.0;
    s1im := 0.0;
    s2im := 0.0;
    s3im := 0.0;
    for i in v0re'range loop
      s0re := s0re + v0re(i);
      s1re := s1re + v1re(i);
      s2re := s2re + v2re(i);
      s3re := s3re + v3re(i);
      s0im := s0im + v0im(i);
      s1im := s1im + v1im(i);
      s2im := s2im + v2im(i);
      s3im := s3im + v3im(i);
    end loop;
  end Sum;

  procedure Product ( x0,x1,x2,x3 : in Standard_Floating_Vectors.Vector;
                      x4,x5,x6,x7 : in Standard_Floating_Vectors.Vector;
                      y0,y1,y2,y3 : in Standard_Floating_Vectors.Vector;
                      y4,y5,y6,y7 : in Standard_Floating_Vectors.Vector;
                      s0,s1,s2,s3,s4,s5,s6,s7 : out double_float ) is
  begin
    s0 := 0.0;
    s1 := 0.0;
    s2 := 0.0;
    s3 := 0.0;
    s4 := 0.0;
    s5 := 0.0;
    s6 := 0.0;
    s7 := 0.0;
    for i in x0'range loop
      s0 := s0 + x0(i)*y0(i);
      s1 := s1 + x0(i)*y1(i) + x1(i)*y0(i);
      s2 := s2 + x0(i)*y2(i) + x1(i)*y1(i) + x2(i)*y0(i);
      s3 := s3 + x0(i)*y3(i) + x1(i)*y2(i) + x2(i)*y1(i) + x3(i)*y0(i);
      s4 := s4 + x0(i)*y4(i) + x1(i)*y3(i) + x2(i)*y2(i) + x3(i)*y1(i)
               + x4(i)*y0(i);
      s5 := s5 + x0(i)*y5(i) + x1(i)*y4(i) + x2(i)*y3(i) + x3(i)*y2(i)
               + x4(i)*y1(i) + x5(i)*y0(i);
      s6 := s6 + x0(i)*y6(i) + x1(i)*y5(i) + x2(i)*y4(i) + x3(i)*y3(i)
               + x4(i)*y2(i) + x5(i)*y1(i) + x6(i)*y0(i);
      s7 := s7 + x0(i)*y7(i) + x1(i)*y6(i) + x2(i)*y5(i) + x3(i)*y4(i)
               + x4(i)*y3(i) + x5(i)*y2(i) + x6(i)*y1(i) + x7(i)*y0(i);
    end loop;
  end Product;

  procedure Product ( x0re,x1re : in Standard_Floating_Vectors.Vector;
                      x2re,x3re : in Standard_Floating_Vectors.Vector;
                      x4re,x5re : in Standard_Floating_Vectors.Vector;
                      x6re,x7re : in Standard_Floating_Vectors.Vector;
                      x0im,x1im : in Standard_Floating_Vectors.Vector;
                      x2im,x3im : in Standard_Floating_Vectors.Vector;
                      x4im,x5im : in Standard_Floating_Vectors.Vector;
                      x6im,x7im : in Standard_Floating_Vectors.Vector;
                      y0re,y1re : in Standard_Floating_Vectors.Vector;
                      y2re,y3re : in Standard_Floating_Vectors.Vector;
                      y4re,y5re : in Standard_Floating_Vectors.Vector;
                      y6re,y7re : in Standard_Floating_Vectors.Vector;
                      y0im,y1im : in Standard_Floating_Vectors.Vector;
                      y2im,y3im : in Standard_Floating_Vectors.Vector;
                      y4im,y5im : in Standard_Floating_Vectors.Vector;
                      y6im,y7im : in Standard_Floating_Vectors.Vector;
                      s0re,s1re,s2re,s3re : out double_float;
                      s4re,s5re,s6re,s7re : out double_float;
                      s0im,s1im,s2im,s3im : out double_float;
                      s4im,s5im,s6im,s7im : out double_float ) is
  begin
    s0re := 0.0; s1re := 0.0; s2re := 0.0; s3re := 0.0;
    s4re := 0.0; s5re := 0.0; s6re := 0.0; s7re := 0.0;
    s0im := 0.0; s1im := 0.0; s2im := 0.0; s3im := 0.0;
    s4im := 0.0; s5im := 0.0; s6im := 0.0; s7im := 0.0;
    for i in x0re'range loop
      s0re := s0re + x0re(i)*y0re(i) - x0im(i)*y0im(i);
      s1re := s1re + x0re(i)*y1re(i) - x0im(i)*y1im(i)
                   + x1re(i)*y0re(i) - x1im(i)*y0im(i);
      s2re := s2re + x0re(i)*y2re(i) - x0im(i)*y2im(i)
                   + x1re(i)*y1re(i) - x1im(i)*y1im(i)
                   + x2re(i)*y0re(i) - x2im(i)*y0im(i);
      s3re := s3re + x0re(i)*y3re(i) - x0im(i)*y3im(i)
                   + x1re(i)*y2re(i) - x1im(i)*y2im(i)
                   + x2re(i)*y1re(i) - x2im(i)*y1im(i)
                   + x3re(i)*y0re(i) - x3im(i)*y0im(i);
      s4re := s4re + x0re(i)*y4re(i) - x0im(i)*y4im(i)
                   + x1re(i)*y3re(i) - x1im(i)*y3im(i)
                   + x2re(i)*y2re(i) - x2im(i)*y2im(i)
                   + x3re(i)*y1re(i) - x3im(i)*y1im(i)
                   + x4re(i)*y0re(i) - x4im(i)*y0im(i);
      s5re := s5re + x0re(i)*y5re(i) - x0im(i)*y5im(i)
                   + x1re(i)*y4re(i) - x1im(i)*y4im(i)
                   + x2re(i)*y3re(i) - x2im(i)*y3im(i)
                   + x3re(i)*y2re(i) - x3im(i)*y2im(i)
                   + x4re(i)*y1re(i) - x4im(i)*y1im(i)
                   + x5re(i)*y0re(i) - x5im(i)*y0im(i);
      s6re := s6re + x0re(i)*y6re(i) - x0im(i)*y6im(i)
                   + x1re(i)*y5re(i) - x1im(i)*y5im(i)
                   + x2re(i)*y4re(i) - x2im(i)*y4im(i)
                   + x3re(i)*y3re(i) - x3im(i)*y3im(i)
                   + x4re(i)*y2re(i) - x4im(i)*y2im(i)
                   + x5re(i)*y1re(i) - x5im(i)*y1im(i)
                   + x6re(i)*y0re(i) - x6im(i)*y0im(i);
      s7re := s7re + x0re(i)*y7re(i) - x0im(i)*y7im(i)
                   + x1re(i)*y6re(i) - x1im(i)*y6im(i)
                   + x2re(i)*y5re(i) - x2im(i)*y5im(i)
                   + x3re(i)*y4re(i) - x3im(i)*y4im(i)
                   + x4re(i)*y3re(i) - x4im(i)*y3im(i)
                   + x5re(i)*y2re(i) - x5im(i)*y2im(i)
                   + x6re(i)*y1re(i) - x6im(i)*y1im(i)
                   + x7re(i)*y0re(i) - x7im(i)*y0im(i);
      s0im := s0im + x0re(i)*y0im(i) + x0im(i)*y0re(i);
      s1im := s1im + x0re(i)*y1im(i) + x1im(i)*y0re(i)
                   + x1re(i)*y0im(i) + x0im(i)*y1re(i);
      s2im := s2im + x0re(i)*y2im(i) + x0im(i)*y2re(i)
                   + x1re(i)*y1im(i) + x1im(i)*y1re(i)
                   + x2re(i)*y0im(i) + x2im(i)*y0re(i);
      s3im := s3im + x0re(i)*y3im(i) + x0im(i)*y3re(i)
                   + x1re(i)*y2im(i) + x1im(i)*y2re(i)
                   + x2re(i)*y1im(i) + x2im(i)*y1re(i)
                   + x3re(i)*y0im(i) + x3im(i)*y0re(i);
      s4im := s4im + x0re(i)*y4im(i) + x0im(i)*y4re(i)
                   + x1re(i)*y3im(i) + x1im(i)*y3re(i)
                   + x2re(i)*y2im(i) + x2im(i)*y2re(i)
                   + x3re(i)*y1im(i) + x3im(i)*y1re(i)
                   + x4re(i)*y0im(i) + x4im(i)*y0re(i);
      s5im := s5im + x0re(i)*y5im(i) + x0im(i)*y5re(i)
                   + x1re(i)*y4im(i) + x1im(i)*y4re(i)
                   + x2re(i)*y3im(i) + x2im(i)*y3re(i)
                   + x3re(i)*y2im(i) + x3im(i)*y2re(i)
                   + x4re(i)*y1im(i) + x4im(i)*y1re(i)
                   + x5re(i)*y0im(i) + x5im(i)*y0re(i);
      s6im := s6im + x0re(i)*y6im(i) + x0im(i)*y6re(i)
                   + x1re(i)*y5im(i) + x1im(i)*y5re(i)
                   + x2re(i)*y4im(i) + x2im(i)*y4re(i)
                   + x3re(i)*y3im(i) + x3im(i)*y3re(i)
                   + x4re(i)*y2im(i) + x4im(i)*y2re(i)
                   + x5re(i)*y1im(i) + x5im(i)*y1re(i)
                   + x6re(i)*y0im(i) + x6im(i)*y0re(i);
      s7im := s7im + x0re(i)*y7im(i) + x0im(i)*y7re(i)
                   + x1re(i)*y6im(i) + x1im(i)*y6re(i)
                   + x2re(i)*y5im(i) + x2im(i)*y5re(i)
                   + x3re(i)*y4im(i) + x3im(i)*y4re(i)
                   + x4re(i)*y3im(i) + x4im(i)*y3re(i)
                   + x5re(i)*y2im(i) + x5im(i)*y2re(i)
                   + x6re(i)*y1im(i) + x6im(i)*y1re(i)
                   + x7re(i)*y0im(i) + x7im(i)*y0re(i);
    end loop;
  end Product;

  function to_double_double
             ( s0,s1,s2,s3 : double_float;
               verbose : boolean := true ) return double_double is

    shi,slo,err : double_float;
    res : double_double;

  begin
    Double_Double_Basics.quick_two_sum(s0,s1,shi,err);
    res := create(shi,err);
    if verbose then
      put("shi : "); put(shi); new_line;
      put("err : "); put(err); new_line;
    end if;
    Double_Double_Basics.quick_two_sum(s2,s3,slo,err);
    res := res + create(slo,err);
    if verbose then
      put("slo : "); put(slo); new_line;
      put("err : "); put(err); new_line;
    end if;
    return res;
  end to_double_double;

  function to_Double_Double
             ( s0,s1,s2,s3,s4,s5,s6,s7 : double_float;
               verbose : boolean := true ) return double_double is

    res : double_double;
    z0,z1,z2,z3 : double_float;
    e1,e2,e3,e4 : double_float;

  begin
    Double_Double_Basics.two_sum(s0,s1,z0,e1);
    if verbose then
      put(" z0 : "); put(z0); new_line;
      put("err : "); put(e1); new_line;
    end if;
    Double_Double_Basics.two_sum(s2,s3,z1,e2);
    if verbose then
      put(" z1 : "); put(z1); new_line;
      put("err : "); put(e2); new_line;
    end if;
    Double_Double_Basics.two_sum(s4,s5,z2,e3);
    if verbose then
      put(" z2 : "); put(z2); new_line;
      put("err : "); put(e3); new_line;
    end if;
    Double_Double_Basics.two_sum(s6,s7,z3,e4);
    if verbose then
      put(" z3 : "); put(z3); new_line;
      put("err : "); put(e4); new_line;
    end if;
    res := Double_Double_Numbers.create(z3);
    res := res + z2;
    res := res + z1;
    res := res + z0;
    return res;
  end to_Double_Double;

  function to_Complex_Double_Double
             ( s0re,s1re,s2re,s3re,s0im,s1im,s2im,s3im : double_float;
               verbose : boolean := true ) return Complex_Number is

    res : Complex_Number;
    sre : constant double_double
        := to_Double_Double(s0re,s1re,s2re,s3re,verbose);
    sim : constant double_double
        := to_Double_Double(s0im,s1im,s2im,s3im,verbose);

  begin
    res := Create(sre,sim);
    return res;
  end to_Complex_Double_Double;

  function to_Complex_Double_Double
             ( s0re,s1re,s2re,s3re,s4re,s5re,s6re,s7re : double_float;
               s0im,s1im,s2im,s3im,s4im,s5im,s6im,s7im : double_float;
               verbose : boolean := true ) return Complex_Number is

    res : Complex_Number;
    sre : constant double_double
        := to_Double_Double(s0re,s1re,s2re,s3re,s4re,s5re,s6re,s7re,verbose);
    sim : constant double_double
        := to_Double_Double(s0im,s1im,s2im,s3im,s4im,s5im,s6im,s7im,verbose);

  begin
    res := Create(sre,sim);
    return res;
  end to_Complex_Double_Double;

end Vectored_Double_Doubles;
