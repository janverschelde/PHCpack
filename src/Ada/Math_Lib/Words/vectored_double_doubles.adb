with text_io;                            use text_io;
with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;
with Double_Double_Basics;

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

  function to_Complex_Double_Double
             ( s0re,s1re,s2re,s3re,s0im,s1im,s2im,s3im : double_float;
               verbose : boolean := true ) return Complex_Number is

    res : Complex_Number;
    sre : constant double_double
        := to_double_double(s0re,s1re,s2re,s3re,verbose);
    sim : constant double_double
        := to_double_double(s0im,s1im,s2im,s3im,verbose);

  begin
    res := Create(sre,sim);
    return res;
  end to_Complex_Double_Double;

end Vectored_Double_Doubles;
