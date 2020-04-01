with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Binomial_Coefficients;              use Binomial_Coefficients;

package body Shift_Convolution_Circuits is

  procedure Shift ( c,wrk : in out Standard_Complex_Vectors.Vector;
                    t : in double_float ) is

    use Standard_Complex_Numbers;

    bcf : double_float;
    sgn : integer32;

  begin
    for i in c'range loop
      wrk(i) := Create(0.0);
      if i mod 2 = 0
       then sgn := 1;
       else sgn := -1;
      end if;
      for j in 0..i loop
       -- bcf := double_float(sgn*binomial(i,j));
        bcf := double_float(sgn)*binomial(i,j);
        bcf := bcf*(t**(natural(i-j)));
        wrk(j) := wrk(j) + c(i)*bcf;
        sgn := -sgn;
      end loop;
    end loop;
    c := wrk;
  end Shift;

  procedure Shift ( c,wrk : in out Standard_Complex_Vectors.Vector;
                    t : in Standard_Complex_Numbers.Complex_Number ) is

    use Standard_Complex_Numbers;

    bcf : double_float;
    rcf : Complex_Number;
    sgn : integer32;

  begin
    for i in c'range loop
      wrk(i) := Create(0.0);
      if i mod 2 = 0
       then sgn := 1;
       else sgn := -1;
      end if;
      for j in 0..i loop
       -- bcf := double_float(sgn*binomial(i,j));
        bcf := double_float(sgn)*binomial(i,j);
        rcf := bcf*(t**(natural(i-j)));
        wrk(j) := wrk(j) + c(i)*rcf;
        sgn := -sgn;
      end loop;
    end loop;
    c := wrk;
  end Shift;

  procedure Shift ( c,wrk : in out DoblDobl_Complex_Vectors.Vector;
                    t : in double_double ) is

    use DoblDobl_Complex_Numbers;

    bcf : double_double;
    sgn : integer32;

  begin
    for i in c'range loop
      wrk(i) := Create(integer32(0));
      if i mod 2 = 0
       then sgn := 1;
       else sgn := -1;
      end if;
      for j in 0..i loop
       -- bcf := double_double_numbers.create(sgn*binomial(i,j));
        bcf := binomial(i,j);
        bcf := double_double_numbers.create(sgn)*bcf;
        bcf := bcf*(t**(natural(i-j)));
        wrk(j) := wrk(j) + c(i)*bcf;
        sgn := -sgn;
      end loop;
    end loop;
    c := wrk;
  end Shift;

  procedure Shift ( c,wrk : in out DoblDobl_Complex_Vectors.Vector;
                    t : in DoblDobl_Complex_Numbers.Complex_Number ) is

    use DoblDobl_Complex_Numbers;

    bcf : double_double;
    rcf : Complex_Number;
    sgn : integer32;

  begin
    for i in c'range loop
      wrk(i) := Create(integer32(0));
      if i mod 2 = 0
       then sgn := 1;
       else sgn := -1;
      end if;
      for j in 0..i loop
       -- bcf := double_double_numbers.create(sgn*binomial(i,j));
        bcf := binomial(i,j);
        bcf := double_double_numbers.create(sgn)*bcf;
        rcf := bcf*(t**(natural(i-j)));
        wrk(j) := wrk(j) + c(i)*rcf;
        sgn := -sgn;
      end loop;
    end loop;
    c := wrk;
  end Shift;

  procedure Shift ( c,wrk : in out QuadDobl_Complex_Vectors.Vector;
                    t : in quad_double ) is

    use QuadDobl_Complex_Numbers;

    bcf : quad_double;
    sgn : integer32;

  begin
    for i in c'range loop
      wrk(i) := Create(integer32(0));
      if i mod 2 = 0
       then sgn := 1;
       else sgn := -1;
      end if;
      for j in 0..i loop
       -- bcf := quad_double_numbers.create(sgn*binomial(i,j));
        bcf := binomial(i,j);
        bcf := quad_double_numbers.create(sgn)*bcf;
        bcf := bcf*(t**(natural(i-j)));
        wrk(j) := wrk(j) + c(i)*bcf;
        sgn := -sgn;
      end loop;
    end loop;
    c := wrk;
  end Shift;

  procedure Shift ( c,wrk : in out QuadDobl_Complex_Vectors.Vector;
                    t : in QuadDobl_Complex_Numbers.Complex_Number ) is

    use QuadDobl_Complex_Numbers;

    bcf : quad_double;
    rcf : Complex_Number;
    sgn : integer32;

  begin
    for i in c'range loop
      wrk(i) := Create(integer32(0));
      if i mod 2 = 0
       then sgn := 1;
       else sgn := -1;
      end if;
      for j in 0..i loop
       -- bcf := quad_double_numbers.create(sgn*binomial(i,j));
        bcf := binomial(i,j);
        bcf := quad_double_numbers.create(sgn)*bcf;
        rcf := bcf*(t**(natural(i-j)));
        wrk(j) := wrk(j) + c(i)*rcf;
        sgn := -sgn;
      end loop;
    end loop;
    c := wrk;
  end Shift;

  procedure Shift ( c,wrk : in Standard_Complex_Vectors.Link_to_Vector;
                    t : in double_float ) is

    use Standard_Complex_Vectors;

  begin
    if c /= null
     then Shift(c.all,wrk.all,t);
    end if;
  end Shift;

  procedure Shift ( c,wrk : in Standard_Complex_Vectors.Link_to_Vector;
                    t : in Standard_Complex_Numbers.Complex_Number ) is

    use Standard_Complex_Vectors;

  begin
    if c /= null
     then Shift(c.all,wrk.all,t);
    end if;
  end Shift;

  procedure Shift ( c,wrk : in DoblDobl_Complex_Vectors.Link_to_Vector;
                    t : in double_double ) is

    use DoblDobl_Complex_Vectors;

  begin
    if c /= null
     then Shift(c.all,wrk.all,t);
    end if;
  end Shift;

  procedure Shift ( c,wrk : in DoblDobl_Complex_Vectors.Link_to_Vector;
                    t : in DoblDobl_Complex_Numbers.Complex_Number ) is

    use DoblDobl_Complex_Vectors;

  begin
    if c /= null
     then Shift(c.all,wrk.all,t);
    end if;
  end Shift;

  procedure Shift ( c,wrk : in QuadDobl_Complex_Vectors.Link_to_Vector;
                    t : in quad_double ) is

    use QuadDobl_Complex_Vectors;

  begin
    if c /= null
     then Shift(c.all,wrk.all,t);
    end if;
  end Shift;

  procedure Shift ( c,wrk : in QuadDobl_Complex_Vectors.Link_to_Vector;
                    t : in QuadDobl_Complex_Numbers.Complex_Number ) is

    use QuadDobl_Complex_Vectors;

  begin
    if c /= null
     then Shift(c.all,wrk.all,t);
    end if;
  end Shift;

  procedure Shift ( c : in out Standard_Speelpenning_Convolutions.Circuit;
                    wrk : in Standard_Complex_Vectors.Link_to_Vector;
                    t : in double_float ) is

    use Standard_Complex_Vectors;
    lnk : Link_to_Vector;

  begin
    for k in c.cff'range loop
      lnk := c.cff(k);
      Shift(lnk,wrk,t);
    end loop;
    if c.cst /= null then
      if c.cst'last > 0
       then Shift(c.cst,wrk,t);
      end if;
    end if;
  end Shift;

  procedure Shift ( c : in out Standard_Speelpenning_Convolutions.Circuit;
                    wrk : in Standard_Complex_Vectors.Link_to_Vector;
                    t : in Standard_Complex_Numbers.Complex_Number ) is

    use Standard_Complex_Vectors;
    lnk : Link_to_Vector;

  begin
    for k in c.cff'range loop
      lnk := c.cff(k);
      Shift(lnk,wrk,t);
    end loop;
    if c.cst /= null then
      if c.cst'last > 0
       then Shift(c.cst,wrk,t);
      end if;
    end if;
  end Shift;

  procedure Shift ( c : in out DoblDobl_Speelpenning_Convolutions.Circuit;
                    wrk : in DoblDobl_Complex_Vectors.Link_to_Vector;
                    t : in double_double ) is

    use DoblDobl_Complex_Vectors;
    lnk : Link_to_Vector;

  begin
    for k in c.cff'range loop
      lnk := c.cff(k);
      Shift(lnk,wrk,t);
    end loop;
    if c.cst /= null then
      if c.cst'last > 0
       then Shift(c.cst,wrk,t);
      end if;
    end if;
  end Shift;

  procedure Shift ( c : in out DoblDobl_Speelpenning_Convolutions.Circuit;
                    wrk : in DoblDobl_Complex_Vectors.Link_to_Vector;
                    t : in DoblDobl_Complex_Numbers.Complex_Number ) is

    use DoblDobl_Complex_Vectors;
    lnk : Link_to_Vector;

  begin
    for k in c.cff'range loop
      lnk := c.cff(k);
      Shift(lnk,wrk,t);
    end loop;
    if c.cst /= null then
      if c.cst'last > 0
       then Shift(c.cst,wrk,t);
      end if;
    end if;
  end Shift;

  procedure Shift ( c : in out QuadDobl_Speelpenning_Convolutions.Circuit;
                    wrk : in QuadDobl_Complex_Vectors.Link_to_Vector;
                    t : in quad_double ) is

    use QuadDobl_Complex_Vectors;
    lnk : Link_to_Vector;

  begin
    for k in c.cff'range loop
      lnk := c.cff(k);
      Shift(lnk,wrk,t);
    end loop;
    if c.cst /= null then
      if c.cst'last > 0
       then Shift(c.cst,wrk,t);
      end if;
    end if;
  end Shift;

  procedure Shift ( c : in out QuadDobl_Speelpenning_Convolutions.Circuit;
                    wrk : in QuadDobl_Complex_Vectors.Link_to_Vector;
                    t : in QuadDobl_Complex_Numbers.Complex_Number ) is

    use QuadDobl_Complex_Vectors;
    lnk : Link_to_Vector;

  begin
    for k in c.cff'range loop
      lnk := c.cff(k);
      Shift(lnk,wrk,t);
    end loop;
    if c.cst /= null then
      if c.cst'last > 0
       then Shift(c.cst,wrk,t);
      end if;
    end if;
  end Shift;

  procedure Shift ( c : in Standard_Speelpenning_Convolutions.Link_to_Circuit;
                    wrk : in Standard_Complex_Vectors.Link_to_Vector;
                    t : in double_float ) is

    use Standard_Speelpenning_Convolutions;

  begin
    if c /= null
     then Shift(c.all,wrk,t);
    end if;
  end Shift;

  procedure Shift ( c : in Standard_Speelpenning_Convolutions.Link_to_Circuit;
                    wrk : in Standard_Complex_Vectors.Link_to_Vector;
                    t : in Standard_Complex_Numbers.Complex_Number ) is

    use Standard_Speelpenning_Convolutions;

  begin
    if c /= null
     then Shift(c.all,wrk,t);
    end if;
  end Shift;

  procedure Shift ( c : in DoblDobl_Speelpenning_Convolutions.Link_to_Circuit;
                    wrk : in DoblDobl_Complex_Vectors.Link_to_Vector;
                    t : in double_double ) is

    use DoblDobl_Speelpenning_Convolutions;

  begin
    if c /= null
     then Shift(c.all,wrk,t);
    end if;
  end Shift;

  procedure Shift ( c : in DoblDobl_Speelpenning_Convolutions.Link_to_Circuit;
                    wrk : in DoblDobl_Complex_Vectors.Link_to_Vector;
                    t : in DoblDobl_Complex_Numbers.Complex_Number ) is

    use DoblDobl_Speelpenning_Convolutions;

  begin
    if c /= null
     then Shift(c.all,wrk,t);
    end if;
  end Shift;

  procedure Shift ( c : in QuadDobl_Speelpenning_Convolutions.Link_to_Circuit;
                    wrk : in QuadDobl_Complex_Vectors.Link_to_Vector;
                    t : in quad_double ) is

    use QuadDobl_Speelpenning_Convolutions;

  begin
    if c /= null
     then Shift(c.all,wrk,t);
    end if;
  end Shift;

  procedure Shift ( c : in QuadDobl_Speelpenning_Convolutions.Link_to_Circuit;
                    wrk : in QuadDobl_Complex_Vectors.Link_to_Vector;
                    t : in QuadDobl_Complex_Numbers.Complex_Number ) is


    use QuadDobl_Speelpenning_Convolutions;

  begin
    if c /= null
     then Shift(c.all,wrk,t);
    end if;
  end Shift;

  procedure Shift ( c : in out Standard_Speelpenning_Convolutions.Circuits;
                    wrk : in Standard_Complex_Vectors.Link_to_Vector;
                    t : in double_float ) is
  begin
    for k in c'range loop
      Shift(c(k),wrk,t);
    end loop;
  end Shift;

  procedure Shift ( c : in Standard_Speelpenning_Convolutions.Link_to_Circuits;
                    wrk : in Standard_Complex_Vectors.Link_to_Vector;
                    t : in double_float ) is

    use Standard_Speelpenning_Convolutions;

  begin
    if c /= null
     then Shift(c.all,wrk,t);
    end if;
  end Shift;

  procedure Shift ( c : in out Standard_Speelpenning_Convolutions.Circuits;
                    wrk : in Standard_Complex_Vectors.Link_to_Vector;
                    t : in Standard_Complex_Numbers.Complex_Number ) is
  begin
    for k in c'range loop
      Shift(c(k),wrk,t);
    end loop;
  end Shift;

  procedure Shift ( c : in Standard_Speelpenning_Convolutions.Link_to_Circuits;
                    wrk : in Standard_Complex_Vectors.Link_to_Vector;
                    t : in Standard_Complex_Numbers.Complex_Number ) is

    use Standard_Speelpenning_Convolutions;

  begin
    if c /= null
     then Shift(c.all,wrk,t);
    end if;
  end Shift;

  procedure Shift ( c : in out DoblDobl_Speelpenning_Convolutions.Circuits;
                    wrk : in DoblDobl_Complex_Vectors.Link_to_Vector;
                    t : in double_double ) is
  begin
    for k in c'range loop
      Shift(c(k),wrk,t);
    end loop;
  end Shift;

  procedure Shift ( c : in DoblDobl_Speelpenning_Convolutions.Link_to_Circuits;
                    wrk : in DoblDobl_Complex_Vectors.Link_to_Vector;
                    t : in double_double ) is

    use DoblDobl_Speelpenning_Convolutions;

  begin
    if c /= null
     then Shift(c.all,wrk,t);
    end if;
  end Shift;

  procedure Shift ( c : in out DoblDobl_Speelpenning_Convolutions.Circuits;
                    wrk : in DoblDobl_Complex_Vectors.Link_to_Vector;
                    t : in DoblDobl_Complex_Numbers.Complex_Number ) is
  begin
    for k in c'range loop
      Shift(c(k),wrk,t);
    end loop;
  end Shift;

  procedure Shift ( c : in DoblDobl_Speelpenning_Convolutions.Link_to_Circuits;
                    wrk : in DoblDobl_Complex_Vectors.Link_to_Vector;
                    t : in DoblDobl_Complex_Numbers.Complex_Number ) is

    use DoblDobl_Speelpenning_Convolutions;

  begin
    if c /= null
     then Shift(c.all,wrk,t);
    end if;
  end Shift;

  procedure Shift ( c : in out QuadDobl_Speelpenning_Convolutions.Circuits;
                    wrk : in QuadDobl_Complex_Vectors.Link_to_Vector;
                    t : in quad_double ) is
  begin
    for k in c'range loop
      Shift(c(k),wrk,t);
    end loop;
  end Shift;

  procedure Shift ( c : in QuadDobl_Speelpenning_Convolutions.Link_to_Circuits;
                    wrk : in QuadDobl_Complex_Vectors.Link_to_Vector;
                    t : in quad_double ) is

    use QuadDobl_Speelpenning_Convolutions;

  begin
    if c /= null
     then Shift(c.all,wrk,t);
    end if;
  end Shift;

  procedure Shift ( c : in out QuadDobl_Speelpenning_Convolutions.Circuits;
                    wrk : in QuadDobl_Complex_Vectors.Link_to_Vector;
                    t : in QuadDobl_Complex_Numbers.Complex_Number ) is
  begin
    for k in c'range loop
      Shift(c(k),wrk,t);
    end loop;
  end Shift;

  procedure Shift ( c : in QuadDobl_Speelpenning_Convolutions.Link_to_Circuits;
                    wrk : in QuadDobl_Complex_Vectors.Link_to_Vector;
                    t : in QuadDobl_Complex_Numbers.Complex_Number ) is

    use QuadDobl_Speelpenning_Convolutions;

  begin
    if c /= null
     then Shift(c.all,wrk,t);
    end if;
  end Shift;

  procedure Shift ( s : in out Standard_Speelpenning_Convolutions.System;
                    wrk : in Standard_Complex_Vectors.Link_to_Vector;
                    t : in double_float ) is
  begin
    Shift(s.crc,wrk,t);
  end Shift;

  procedure Shift ( s : in Standard_Speelpenning_Convolutions.Link_to_System;
                    wrk : in Standard_Complex_Vectors.Link_to_Vector;
                    t : in double_float ) is

    use Standard_Speelpenning_Convolutions;

  begin
    if s /= null
     then Shift(s.crc,wrk,t);
    end if;
  end Shift;

  procedure Shift ( s : in out Standard_Speelpenning_Convolutions.System;
                    wrk : in Standard_Complex_Vectors.Link_to_Vector;
                    t : in Standard_Complex_Numbers.Complex_Number ) is
  begin
    Shift(s.crc,wrk,t);
  end Shift;

  procedure Shift ( s : in Standard_Speelpenning_Convolutions.Link_to_System;
                    wrk : in Standard_Complex_Vectors.Link_to_Vector;
                    t : in Standard_Complex_Numbers.Complex_Number ) is

    use Standard_Speelpenning_Convolutions;

  begin
    if s /= null
     then Shift(s.crc,wrk,t);
    end if;
  end Shift;

  procedure Shift ( s : in out DoblDobl_Speelpenning_Convolutions.System;
                    wrk : in DoblDobl_Complex_Vectors.Link_to_Vector;
                    t : in double_double ) is
  begin
    Shift(s.crc,wrk,t);
  end Shift;

  procedure Shift ( s : in DoblDobl_Speelpenning_Convolutions.Link_to_System;
                    wrk : in DoblDobl_Complex_Vectors.Link_to_Vector;
                    t : in double_double ) is

    use DoblDobl_Speelpenning_Convolutions;

  begin
    if s /= null
     then Shift(s.crc,wrk,t);
    end if;
  end Shift;

  procedure Shift ( s : in out DoblDobl_Speelpenning_Convolutions.System;
                    wrk : in DoblDobl_Complex_Vectors.Link_to_Vector;
                    t : in DoblDobl_Complex_Numbers.Complex_Number ) is
  begin
    Shift(s.crc,wrk,t);
  end Shift;

  procedure Shift ( s : in DoblDobl_Speelpenning_Convolutions.Link_to_System;
                    wrk : in DoblDobl_Complex_Vectors.Link_to_Vector;
                    t : in DoblDobl_Complex_Numbers.Complex_Number ) is

    use DoblDobl_Speelpenning_Convolutions;

  begin
    if s /= null
     then Shift(s.crc,wrk,t);
    end if;
  end Shift;

  procedure Shift ( s : in out QuadDobl_Speelpenning_Convolutions.System;
                    wrk : in QuadDobl_Complex_Vectors.Link_to_Vector;
                    t : in quad_double ) is
  begin
    Shift(s.crc,wrk,t);
  end Shift;

  procedure Shift ( s : in QuadDobl_Speelpenning_Convolutions.Link_to_System;
                    wrk : in QuadDobl_Complex_Vectors.Link_to_Vector;
                    t : in quad_double ) is

    use QuadDobl_Speelpenning_Convolutions;

  begin
    if s /= null
     then Shift(s.crc,wrk,t);
    end if;
  end Shift;

  procedure Shift ( s : in out QuadDobl_Speelpenning_Convolutions.System;
                    wrk : in QuadDobl_Complex_Vectors.Link_to_Vector;
                    t : in QuadDobl_Complex_Numbers.Complex_Number ) is
  begin
    Shift(s.crc,wrk,t);
  end Shift;

  procedure Shift ( s : in QuadDobl_Speelpenning_Convolutions.Link_to_System;
                    wrk : in QuadDobl_Complex_Vectors.Link_to_Vector;
                    t : in QuadDobl_Complex_Numbers.Complex_Number ) is

    use QuadDobl_Speelpenning_Convolutions;

  begin
    if s /= null
     then Shift(s.crc,wrk,t);
    end if;
  end Shift;

end Shift_Convolution_Circuits;
