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

end Shift_Coefficient_Convolutions;
