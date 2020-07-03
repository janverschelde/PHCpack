with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Binomial_Coefficients;

package body Shift_Coefficient_Convolutions is

  procedure Powers_of_Shift
              ( pwt : in Standard_Floating_Vectors.Link_to_Vector;
                t : in double_float ) is
  begin
    pwt(0) := 1.0;
    pwt(1) := t;
    for k in 2..pwt'last loop
      pwt(k) := t*pwt(k-1);
    end loop;
  end Powers_of_Shift;

  procedure Powers_of_Shift
              ( rpwt,ipwt : in Standard_Floating_Vectors.Link_to_Vector;
                rpt,ipt : in double_float ) is

    pr,pi : double_float;

  begin
    rpwt(0) := 1.0; ipwt(0) := 0.0;
    rpwt(1) := rpt; ipwt(1) := ipt;
    for k in 2..rpwt'last loop  -- pwt(k) := t*pwt(k-1);
      pr := rpwt(k-1); pi := ipwt(k-1);
      rpwt(k) := pr*rpt - pi*ipt;
      ipwt(k) := pr*ipt + pi*rpt;
    end loop;
  end Powers_of_Shift;

  procedure Shift ( rcf,icf : in Standard_Floating_Vectors.Link_to_Vector;
                    rwk,iwk : in Standard_Floating_Vectors.Link_to_Vector;
                    pwt : in Standard_Floating_Vectors.Link_to_Vector ) is

    bcf : double_float;
    sgn : integer32;

  begin
    for i in rcf'range loop
      rwk(i) := 0.0;
      iwk(i) := 0.0;
      if i mod 2 = 0
       then sgn := 1;
       else sgn := -1;
      end if;
      for j in 0..i loop
        bcf := double_float(sgn)*Binomial_Coefficients.binomial(i,j);
        bcf := bcf*pwt(i-j);
        rwk(j) := rwk(j) + rcf(i)*bcf;
        iwk(j) := iwk(j) + icf(i)*bcf;
        sgn := -sgn;
      end loop;
    end loop;
    for i in rcf'range loop
      rcf(i) := rwk(i);
      icf(i) := iwk(i);
    end loop;
  end Shift;

  procedure Map ( rcf,icf : in Standard_Floating_Vectors.Link_to_Vector;
                  rsh,ish : in Standard_Floating_Vectors.Link_to_Vector;
                  pwt : in Standard_Floating_Vectors.Link_to_Vector ) is

    bcf : double_float;
    sgn : integer32;

  begin
    for i in rcf'range loop
      rsh(i) := 0.0;
      ish(i) := 0.0;
      if i mod 2 = 0
       then sgn := 1;
       else sgn := -1;
      end if;
      for j in 0..i loop
        bcf := double_float(sgn)*Binomial_Coefficients.binomial(i,j);
        bcf := bcf*pwt(i-j);
        rsh(j) := rsh(j) + rcf(i)*bcf;
        ish(j) := ish(j) + icf(i)*bcf;
        sgn := -sgn;
      end loop;
    end loop;
  end Map;

  procedure Shift ( rcf,icf : in Standard_Floating_Vectors.Link_to_Vector;
                    rwk,iwk : in Standard_Floating_Vectors.Link_to_Vector;
                    rpwt : in Standard_Floating_Vectors.Link_to_Vector;
                    ipwt : in Standard_Floating_Vectors.Link_to_Vector ) is

    bcf,pr,pi,zr,zi : double_float;
    sgn : integer32;

  begin
    for i in rcf'range loop
      rwk(i) := 0.0;
      iwk(i) := 0.0;
      if i mod 2 = 0
       then sgn := 1;
       else sgn := -1;
      end if;
      for j in 0..i loop
        bcf := double_float(sgn)*Binomial_Coefficients.binomial(i,j);
       -- bcf := bcf*(t**(natural(i-j)));
        pr := bcf*rpwt(i-j);
        pi := bcf*ipwt(i-j);
        zr := rcf(i)*pr - icf(i)*pi;
        zi := rcf(i)*pi + icf(i)*pr;
        rwk(j) := rwk(j) + zr;
        iwk(j) := iwk(j) + zi;
        sgn := -sgn;
      end loop;
    end loop;
    for i in rcf'range loop
      rcf(i) := rwk(i);
      icf(i) := iwk(i);
    end loop;
  end Shift;

  procedure Map ( rcf,icf : in Standard_Floating_Vectors.Link_to_Vector;
                  rsh,ish : in Standard_Floating_Vectors.Link_to_Vector;
                  rpwt : in Standard_Floating_Vectors.Link_to_Vector;
                  ipwt : in Standard_Floating_Vectors.Link_to_Vector ) is

    bcf,pr,pi,zr,zi : double_float;
    sgn : integer32;

  begin
    for i in rcf'range loop
      rsh(i) := 0.0;
      ish(i) := 0.0;
      if i mod 2 = 0
       then sgn := 1;
       else sgn := -1;
      end if;
      for j in 0..i loop
        bcf := double_float(sgn)*Binomial_Coefficients.binomial(i,j);
       -- bcf := bcf*(t**(natural(i-j)));
        pr := bcf*rpwt(i-j);
        pi := bcf*ipwt(i-j);
        zr := rcf(i)*pr - icf(i)*pi;
        zi := rcf(i)*pi + icf(i)*pr;
        rsh(j) := rsh(j) + zr;
        ish(j) := ish(j) + zi;
        sgn := -sgn;
      end loop;
    end loop;
  end Map;

  procedure Shift ( rcf,icf : in Standard_Floating_VecVecs.VecVec;
                    rwk,iwk : in Standard_Floating_Vectors.Link_to_Vector;
                    pwt : in Standard_Floating_Vectors.Link_to_Vector ) is
  begin
    for k in rcf'range loop
      Shift(rcf(k),icf(k),rwk,iwk,pwt);
    end loop;
  end Shift;

  procedure Map ( rcf,icf : in Standard_Floating_VecVecs.VecVec;
                  rsh,ish : in Standard_Floating_VecVecs.VecVec;
                  pwt : in Standard_Floating_Vectors.Link_to_Vector ) is
  begin
    for k in rcf'range loop
      Map(rcf(k),icf(k),rsh(k),ish(k),pwt);
    end loop;
  end Map;

  procedure Shift ( rcf,icf : in Standard_Floating_VecVecs.VecVec;
                    rwk,iwk : in Standard_Floating_Vectors.Link_to_Vector;
                    rpwt : in Standard_Floating_Vectors.Link_to_Vector;
                    ipwt : in Standard_Floating_Vectors.Link_to_Vector ) is
  begin
    for k in rcf'range loop
      Shift(rcf(k),icf(k),rwk,iwk,rpwt,ipwt);
    end loop;
  end Shift;

  procedure Map ( rcf,icf : in Standard_Floating_VecVecs.VecVec;
                  rsh,ish : in Standard_Floating_VecVecs.VecVec;
                  rpwt : in Standard_Floating_Vectors.Link_to_Vector;
                  ipwt : in Standard_Floating_Vectors.Link_to_Vector ) is
  begin
    for k in rcf'range loop
      Map(rcf(k),icf(k),rsh(k),ish(k),rpwt,ipwt);
    end loop;
  end Map;

-- SHIFTING CIRCUITS WITH A REAL/COMPLEX VALUE :

  procedure Shift ( c : in Standard_Coefficient_Convolutions.Circuit;
                    rwk,iwk : in Standard_Floating_Vectors.Link_to_Vector;
                    pwt : in Standard_Floating_Vectors.Link_to_Vector ) is

    use Standard_Floating_Vectors;

  begin
    Shift(c.rcf,c.icf,rwk,iwk,pwt);
    if c.rct /= null then
      Shift(c.rct,c.ict,rwk,iwk,pwt);
    end if;
  end Shift;

  procedure Shift ( c : in Standard_Coefficient_Convolutions.Circuit;
                    rwk,iwk : in Standard_Floating_Vectors.Link_to_Vector;
                    rpwt : in Standard_Floating_Vectors.Link_to_Vector;
                    ipwt : in Standard_Floating_Vectors.Link_to_Vector ) is

    use Standard_Floating_Vectors;

  begin
    Shift(c.rcf,c.icf,rwk,iwk,rpwt,ipwt);
    if c.rct /= null then
      Shift(c.rct,c.ict,rwk,iwk,rpwt,ipwt);
    end if;
  end Shift;

  procedure Shift ( c : in Standard_Coefficient_Convolutions.Link_to_Circuit;
                    rwk,iwk : in Standard_Floating_Vectors.Link_to_Vector;
                    pwt : in Standard_Floating_Vectors.Link_to_Vector ) is

    use Standard_Coefficient_Convolutions;

  begin
    if c /= null
     then Shift(c.all,rwk,iwk,pwt);
    end if;
  end Shift;

  procedure Shift ( c : in Standard_Coefficient_Convolutions.Link_to_Circuit;
                    rwk,iwk : in Standard_Floating_Vectors.Link_to_Vector;
                    rpwt : in Standard_Floating_Vectors.Link_to_Vector;
                    ipwt : in Standard_Floating_Vectors.Link_to_Vector ) is

    use Standard_Coefficient_Convolutions;

  begin
    if c /= null
     then Shift(c.all,rwk,iwk,rpwt,ipwt);
    end if;
  end Shift;

-- MAPPING COEFFICIENTS SHIFTED WITH REAL VALUE INTO A CIRCUIT :

  procedure Map ( rcf,icf : in Standard_Floating_VecVecs.VecVec;
                  c : in Standard_Coefficient_Convolutions.Circuit;
                  pwt : in Standard_Floating_Vectors.Link_to_Vector ) is

    use Standard_Floating_Vectors;

  begin
    if c.rct /= null then
      Map(rcf(0),icf(0),c.rct,c.ict,pwt);
    end if;
    for k in 1..c.nbr loop
      Map(rcf(k),icf(k),c.rcf(k),c.icf(k),pwt);
    end loop;
  end Map;

  procedure Map ( rcf,icf : in Standard_Floating_VecVecs.VecVec;
                  c : in Standard_Coefficient_Convolutions.Link_to_Circuit;
                  pwt : in Standard_Floating_Vectors.Link_to_Vector ) is

    use Standard_Floating_Vectors;
    use Standard_Coefficient_Convolutions;

  begin
    if c /= null then
      if c.rct /= null then
        Map(rcf(0),icf(0),c.rct,c.ict,pwt);
      end if;
      for k in 1..c.nbr loop
        Map(rcf(k),icf(k),c.rcf(k),c.icf(k),pwt);
      end loop;
    end if;
  end Map;

  procedure Map ( rcf,icf : in Standard_Floating_VecVecs.VecVec;
                  c : in Standard_Coefficient_Convolutions.Circuit;
                  rpwt : in Standard_Floating_Vectors.Link_to_Vector;
                  ipwt : in Standard_Floating_Vectors.Link_to_Vector ) is

    use Standard_Floating_Vectors;

  begin
    if c.rct /= null then
      Map(rcf(0),icf(0),c.rct,c.ict,rpwt,ipwt);
    end if;
    for k in 1..c.nbr loop
      Map(rcf(k),icf(k),c.rcf(k),c.icf(k),rpwt,ipwt);
    end loop;
  end Map;

  procedure Map ( rcf,icf : in Standard_Floating_VecVecs.VecVec;
                  c : in Standard_Coefficient_Convolutions.Link_to_Circuit;
                  rpwt : in Standard_Floating_Vectors.Link_to_Vector;
                  ipwt : in Standard_Floating_Vectors.Link_to_Vector ) is

    use Standard_Floating_Vectors;
    use Standard_Coefficient_Convolutions;

  begin
    if c /= null then
      if c.rct /= null then
        Map(rcf(0),icf(0),c.rct,c.ict,rpwt,ipwt);
      end if;
      for k in 1..c.nbr loop
        Map(rcf(k),icf(k),c.rcf(k),c.icf(k),rpwt,ipwt);
      end loop;
    end if;
  end Map;

-- SHIFTING MANY CIRCUITS :

  procedure Map ( rcf : in Standard_Floating_VecVecVecs.VecVecVec;
                  icf : in Standard_Floating_VecVecVecs.VecVecVec;
                  c : in Standard_Coefficient_Convolutions.Circuits;
                  pwt : in Standard_Floating_Vectors.Link_to_Vector ) is
  begin
    for k in c'range loop
      Map(rcf(k).all,icf(k).all,c(k),pwt);
    end loop;
  end Map;

  procedure Map ( rcf : in Standard_Floating_VecVecVecs.Link_to_VecVecVec;
                  icf : in Standard_Floating_VecVecVecs.Link_to_VecVecVec;
                  c : in Standard_Coefficient_Convolutions.Circuits;
                  pwt : in Standard_Floating_Vectors.Link_to_Vector ) is
  begin
    for k in c'range loop
      Map(rcf(k).all,icf(k).all,c(k),pwt);
    end loop;
  end Map;

  procedure Map ( rcf : in Standard_Floating_VecVecVecs.VecVecVec;
                  icf : in Standard_Floating_VecVecVecs.VecVecVec;
                  c : in Standard_Coefficient_Convolutions.Circuits;
                  rpwt : in Standard_Floating_Vectors.Link_to_Vector;
                  ipwt : in Standard_Floating_Vectors.Link_to_Vector ) is
  begin
    for k in c'range loop
      Map(rcf(k).all,icf(k).all,c(k),rpwt,ipwt);
    end loop;
  end Map;

  procedure Map ( rcf : in Standard_Floating_VecVecVecs.Link_to_VecVecVec;
                  icf : in Standard_Floating_VecVecVecs.Link_to_VecVecVec;
                  c : in Standard_Coefficient_Convolutions.Circuits;
                  rpwt : in Standard_Floating_Vectors.Link_to_Vector;
                  ipwt : in Standard_Floating_Vectors.Link_to_Vector ) is
  begin
    for k in c'range loop
      Map(rcf(k).all,icf(k).all,c(k),rpwt,ipwt);
    end loop;
  end Map;

  procedure Shift ( c : in Standard_Coefficient_Convolutions.Circuits;
                    rwk,iwk : in Standard_Floating_Vectors.Link_to_Vector;
                    pwt : in Standard_Floating_Vectors.Link_to_Vector ) is
  begin
    for k in c'range loop
      Shift(c(k),rwk,iwk,pwt);
    end loop;
  end Shift;

  procedure Shift ( c : in Standard_Coefficient_Convolutions.Circuits;
                    rwk,iwk : in Standard_Floating_Vectors.Link_to_Vector;
                    rpwt : in Standard_Floating_Vectors.Link_to_Vector;
                    ipwt : in Standard_Floating_Vectors.Link_to_Vector ) is
  begin
    for k in c'range loop
      Shift(c(k),rwk,iwk,rpwt,ipwt);
    end loop;
  end Shift;

-- SHIFTING SYSTEMS OF CIRCUITS :

  procedure Shift ( s : in Standard_Coefficient_Convolutions.System;
                    rwk,iwk : in Standard_Floating_Vectors.Link_to_Vector;
                    pwt : in Standard_Floating_Vectors.Link_to_Vector ) is
  begin
    Shift(s.crc,rwk,iwk,pwt);
  end Shift;

  procedure Shift ( s : in Standard_Coefficient_Convolutions.Link_to_System;
                    rwk,iwk : in Standard_Floating_Vectors.Link_to_Vector;
                    pwt : in Standard_Floating_Vectors.Link_to_Vector ) is

    use Standard_Coefficient_Convolutions;

  begin
    if s /= null
     then Shift(s.crc,rwk,iwk,pwt);
    end if;
  end Shift;

  procedure Shift ( s : in Standard_Coefficient_Convolutions.System;
                    rwk,iwk : in Standard_Floating_Vectors.Link_to_Vector;
                    rpwt : in Standard_Floating_Vectors.Link_to_Vector;
                    ipwt : in Standard_Floating_Vectors.Link_to_Vector ) is
  begin
    Shift(s.crc,rwk,iwk,rpwt,ipwt);
  end Shift;

  procedure Shift ( s : in Standard_Coefficient_Convolutions.Link_to_System;
                    rwk,iwk : in Standard_Floating_Vectors.Link_to_Vector;
                    rpwt : in Standard_Floating_Vectors.Link_to_Vector;
                    ipwt : in Standard_Floating_Vectors.Link_to_Vector ) is

    use Standard_Coefficient_Convolutions;

  begin
    if s /= null
     then Shift(s.crc,rwk,iwk,rpwt,ipwt);
    end if;
  end Shift;

  procedure Map ( rcf : in Standard_Floating_VecVecVecs.VecVecVec;
                  icf : in Standard_Floating_VecVecVecs.VecVecVec;
                  s : in Standard_Coefficient_Convolutions.System;
                  pwt : in Standard_Floating_Vectors.Link_to_Vector ) is
  begin
    Map(rcf,icf,s.crc,pwt);
  end Map;

  procedure Map ( rcf : in Standard_Floating_VecVecVecs.VecVecVec;
                  icf : in Standard_Floating_VecVecVecs.VecVecVec;
                  s : in Standard_Coefficient_Convolutions.Link_to_System;
                  pwt : in Standard_Floating_Vectors.Link_to_Vector ) is

    use Standard_Coefficient_Convolutions;

  begin
    if s /= null
     then Map(rcf,icf,s.crc,pwt);
    end if;
  end Map;

  procedure Map ( rcf : in Standard_Floating_VecVecVecs.Link_to_VecVecVec;
                  icf : in Standard_Floating_VecVecVecs.Link_to_VecVecVec;
                  s : in Standard_Coefficient_Convolutions.Link_to_System;
                  pwt : in Standard_Floating_Vectors.Link_to_Vector ) is

    use Standard_Coefficient_Convolutions;

  begin
    if s /= null
     then Map(rcf,icf,s.crc,pwt);
    end if;
  end Map;

  procedure Map ( rcf : in Standard_Floating_VecVecVecs.VecVecVec;
                  icf : in Standard_Floating_VecVecVecs.VecVecVec;
                  s : in Standard_Coefficient_Convolutions.System;
                  rpwt : in Standard_Floating_Vectors.Link_to_Vector;
                  ipwt : in Standard_Floating_Vectors.Link_to_Vector ) is
  begin
    Map(rcf,icf,s.crc,rpwt,ipwt);
  end Map;

  procedure Map ( rcf : in Standard_Floating_VecVecVecs.VecVecVec;
                  icf : in Standard_Floating_VecVecVecs.VecVecVec;
                  s : in Standard_Coefficient_Convolutions.Link_to_System;
                  rpwt : in Standard_Floating_Vectors.Link_to_Vector;
                  ipwt : in Standard_Floating_Vectors.Link_to_Vector ) is

    use Standard_Coefficient_Convolutions;

  begin
    if s /= null
     then Map(rcf,icf,s.crc,rpwt,ipwt);
    end if;
  end Map;

  procedure Map ( rcf : in Standard_Floating_VecVecVecs.Link_to_VecVecVec;
                  icf : in Standard_Floating_VecVecVecs.Link_to_VecVecVec;
                  s : in Standard_Coefficient_Convolutions.Link_to_System;
                  rpwt : in Standard_Floating_Vectors.Link_to_Vector;
                  ipwt : in Standard_Floating_Vectors.Link_to_Vector ) is

    use Standard_Coefficient_Convolutions;

  begin
    if s /= null
     then Map(rcf,icf,s.crc,rpwt,ipwt);
    end if;
  end Map;

end Shift_Coefficient_Convolutions;
