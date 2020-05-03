with Standard_Integer_Numbers;            use Standard_Integer_Numbers;
with Standard_Floating_Numbers;           use Standard_Floating_Numbers;
with Standard_Complex_Numbers;

package body Standard_Coefficient_Convolutions is

  function Real_Part ( x : Standard_Complex_Vectors.Link_to_Vector )
                     return Standard_Floating_Vectors.Link_to_Vector is

    res : Standard_Floating_Vectors.Link_to_Vector;
    rpx : Standard_Floating_Vectors.Vector(x'range);

  begin
    for k in x'range loop
      rpx(k) := Standard_Complex_Numbers.REAL_PART(x(k));
    end loop;
    res := new Standard_Floating_Vectors.Vector'(rpx);
    return res;
  end Real_Part;

  function Imag_Part ( x : Standard_Complex_Vectors.Link_to_Vector )
                     return Standard_Floating_Vectors.Link_to_Vector is

    res : Standard_Floating_Vectors.Link_to_Vector;
    ipx : Standard_Floating_Vectors.Vector(x'range);

  begin
    for k in x'range loop
      ipx(k) := Standard_Complex_Numbers.IMAG_PART(x(k));
    end loop;
    res := new Standard_Floating_Vectors.Vector'(ipx);
    return res;
  end Imag_Part;

  function Make_Complex
             ( rpx,ipx : Standard_Floating_Vectors.Link_to_Vector )
             return Standard_Complex_Vectors.Link_to_Vector is

    res : Standard_Complex_Vectors.Link_to_Vector;
    cvx : Standard_Complex_Vectors.Vector(rpx'range);

  begin
    for k in cvx'range loop
      cvx(k) := Standard_Complex_Numbers.Create(rpx(k),ipx(k));
    end loop;
    res := new Standard_Complex_Vectors.Vector'(cvx);
    return res;
  end Make_Complex;

  procedure Multiply
              ( xr,xi,yr,yi : in Standard_Floating_Vectors.Link_to_Vector;
                zr,zi : in Standard_Floating_Vectors.Link_to_Vector ) is

    deg : constant integer32 := xr'last;
    rpa,ipa : double_float; -- accumulates real and imaginary parts
    xr0,xi0 : double_float; -- to hold values in xr and xi
    yr0,yi0 : double_float; -- to hold values in yr and yi
    idx : integer32;

  begin
   -- product(0) := first(0)*second(0);
    xr0 := xr(0); xi0 := xi(0); yr0 := yr(0); yi0 := yi(0);
    zr(0) := xr0*yr0 - xi0*yi0; zi(0) := xi0*yr0 + xr0*yi0;
    for k in 1..deg loop
     -- product(k) := first(0)*second(k);
      xr0 := xr(0); xi0 := xi(0); yr0 := yr(k); yi0 := yi(k);
      rpa := xr0*yr0 - xi0*yi0; ipa := xi0*yr0 + xr0*yi0;
      for i in 1..k loop
       -- product(k) := product(k) + first(i)*second(k-i);
        xr0 := xr(i); xi0 := xi(i);
        idx := k-i;
        yr0 := yr(idx); yi0 := yi(idx);
        rpa := rpa + xr0*yr0 - xi0*yi0; ipa := ipa + xi0*yr0 + xr0*yi0;
      end loop;
      zr(k) := rpa; zi(k) := ipa;
    end loop;
  end Multiply;

end Standard_Coefficient_Convolutions;
