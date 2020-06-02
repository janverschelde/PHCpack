with Standard_Integer_Numbers;            use Standard_Integer_Numbers;

package body Standard_Hessian_Updaters is

  procedure Speel1 ( H : in out Standard_Complex_Matrices.Matrix;
                     c : in Complex_Number;
                     xps : in Standard_Integer_Vectors.Vector;
                     idx : in Standard_Integer_Vectors.Vector;
                     fac : in Standard_Integer_Vectors.Vector;
                     x : in Standard_Complex_Vectors.Vector;
                     pwt : in Standard_Complex_VecVecs.VecVec ) is

    m1 : integer32;
    powfac : double_float; -- multiplier factor of two powers

  begin
    m1 := xps(fac(1)); -- the monomial is c*x**m1, m1 >= 2.
    powfac := double_float(m1*(m1-1));
    if m1 = 2 then
      H(idx(1),idx(1)) := H(idx(1),idx(1)) + c*powfac;
    elsif m1 = 3 then
      H(idx(1),idx(1)) := H(idx(1),idx(1)) + c*powfac*x(fac(1));
    else -- m > 3, if m = 4, then x**2 at pwt(fac(1))(1)
      H(idx(1),idx(1)) := H(idx(1),idx(1)) + c*powfac*(pwt(fac(1))(m1-3));
    end if;
  end Speel1;

  procedure Speel1 ( hrp : in Standard_Floating_VecVecs.VecVec;
                     hip : in Standard_Floating_VecVecs.VecVec;
                     rcff,icff : in double_float;
                     xps : in Standard_Integer_Vectors.Vector;
                     idx : in Standard_Integer_Vectors.Vector;
                     fac : in Standard_Integer_Vectors.Vector;
                     xr : in Standard_Floating_Vectors.Vector;
                     xi : in Standard_Floating_Vectors.Vector;
                     rpwt : in Standard_Floating_VecVecs.VecVec;
                     ipwt : in Standard_Floating_VecVecs.VecVec ) is

    powfac : double_float; -- multiplier factor of two powers
    hrprow,hiprow : Standard_Floating_Vectors.Link_to_Vector;
    idx1 : constant integer32 := idx(1);
    pr,pi,zr,zi : double_float;
    fac1 : constant integer32 := fac(1);
    m1 : constant integer32 := xps(fac1); -- we have c*x**m1, m1 >= 2

  begin
    powfac := double_float(m1*(m1-1));
    hrprow := hrp(idx1); hiprow := hip(idx1);
    if m1 = 2 then
     -- H(idx(1),idx(1)) := H(idx(1),idx(1)) + c*powfac;
      hrprow(idx1) := hrprow(idx1) + rcff*powfac;
      hiprow(idx1) := hiprow(idx1) + icff*powfac;
    elsif m1 = 3 then
     -- H(idx(1),idx(1)) := H(idx(1),idx(1)) + c*powfac*x(fac(1));
      pr := xr(fac1); pi := xi(fac1);
      zr := pr*rcff - pi*icff; zi := pr*icff + pi*rcff;
      hrprow(idx1) := hrprow(idx1) + powfac*zr;
      hiprow(idx1) := hiprow(idx1) + powfac*zi;
    else -- m > 3, if m = 4, then x**2 at pwt(fac(1))(1)
     -- H(idx(1),idx(1)) := H(idx(1),idx(1)) + c*powfac*(pwt(fac(1))(m1-3));
      pr := rpwt(fac1)(m1-3); pi := ipwt(fac1)(m1-3);
      zr := pr*rcff - pi*icff; zi := pr*icff + pi*rcff;
      hrprow(idx1) := hrprow(idx1) + powfac*zr;
      hiprow(idx1) := hiprow(idx1) + powfac*zi;
    end if;
  end Speel1;

  procedure Speel2 ( H : in out Standard_Complex_Matrices.Matrix;
                     c : in Complex_Number;
                     xps : in Standard_Integer_Vectors.Vector;
                     idx : in Standard_Integer_Vectors.Vector;
                     fac : in Standard_Integer_Vectors.Vector;
                     x : in Standard_Complex_Vectors.Vector;
                     pwt : in Standard_Complex_VecVecs.VecVec ) is

    m1,m2 : integer32;
    powfac : double_float; -- multiplier factor of two powers
    acc : Complex_Number;

  begin
    m1 := xps(fac(1)); powfac := double_float(m1*(m1-1));
    if fac'last = 1 then -- the other variable has no higher power
      if fac(1) = idx(1) -- do not forget to multiply with the other var
       then acc := c*powfac*x(idx(2));
       else acc := c*powfac*x(idx(1));
      end if;
    else -- the other variable appears with a higher power
      m2 := xps(fac(2));
      acc := c*powfac*pwt(fac(2))(m2-1);
    end if;
    if m1 = 2 then
      H(fac(1),fac(1)) := H(fac(1),fac(1)) + acc;
    elsif m1 = 3 then
      H(fac(1),fac(1)) := H(fac(1),fac(1)) + acc*x(fac(1));
    else
      H(fac(1),fac(1)) := H(fac(1),fac(1)) + acc*(pwt(fac(1))(m1-3));
    end if;
    if fac'last = 1 then
      powfac := double_float(m1);
      if m1 = 2 then
        H(idx(1),idx(2)) := H(idx(1),idx(2)) + c*powfac*x(fac(1));
      else
        H(idx(1),idx(2)) := H(idx(1),idx(2)) + c*powfac*(pwt(fac(1))(m1-2));
      end if;
    else
      powfac := double_float(m2*(m2-1));
      acc := c*powfac*(pwt(fac(1))(m1-1));
      if m2 = 2 then
        H(fac(2),fac(2)) := H(fac(2),fac(2)) + acc;
      elsif m2 = 3 then
        H(fac(2),fac(2)) := H(fac(2),fac(2)) + acc*x(fac(2));
      else
        H(fac(2),fac(2)) := H(fac(2),fac(2)) + acc*(pwt(fac(2))(m2-3));
      end if;
      powfac := double_float(m1*m2);
      acc := c*powfac;
      if m1 = 2
       then acc := acc*x(fac(1));
       else acc := acc*(pwt(fac(1))(m1-2));
      end if;
      if m2 = 2
       then acc := acc*x(fac(2));
       else acc := acc*(pwt(fac(2))(m2-2));
      end if;
      H(idx(1),idx(2)) := H(idx(1),idx(2)) + acc;
    end if;
  end Speel2;

  procedure Speel2 ( hrp : in Standard_Floating_VecVecs.VecVec;
                     hip : in Standard_Floating_VecVecs.VecVec;
                     rcff,icff : in double_float;
                     xps : in Standard_Integer_Vectors.Vector;
                     idx : in Standard_Integer_Vectors.Vector;
                     fac : in Standard_Integer_Vectors.Vector;
                     xr : in Standard_Floating_Vectors.Vector;
                     xi : in Standard_Floating_Vectors.Vector;
                     rpwt : in Standard_Floating_VecVecs.VecVec;
                     ipwt : in Standard_Floating_VecVecs.VecVec ) is

    m1,m2,fac2 : integer32;
    powfac : double_float; -- multiplier factor of two powers
    pr,pi,racc,iacc,zr,zi : double_float;
    idx1 : constant integer32 := idx(1);
    idx2 : constant integer32 := idx(2);
    fac1 : constant integer32 := fac(1);
    hrprow,hiprow : Standard_Floating_Vectors.Link_to_Vector;

  begin
    m1 := xps(fac1); powfac := double_float(m1*(m1-1));
    if fac'last = 1 then -- the other variable has no higher power
      if fac1 = idx1 then -- do not forget to multiply with the other var
       -- acc := c*powfac*x(idx(2));
        pr := xr(idx2); pi := xi(idx2);
        racc := rcff*pr - icff*pi; iacc := rcff*pi + icff*pr;
        racc := powfac*racc; iacc := powfac*iacc;
      else
       -- acc := c*powfac*x(idx(1));
        pr := xr(idx1); pi := xi(idx1);
        racc := rcff*pr - icff*pi; iacc := rcff*pi + icff*pr;
        racc := powfac*racc; iacc := powfac*iacc;
      end if;
    else -- the other variable appears with a higher power
      fac2 := fac(2);
      m2 := xps(fac2);
     -- acc := c*powfac*pwt(fac(2))(m2-1);
      pr := rpwt(fac2)(m2-1); pi := ipwt(fac2)(m2-1);
      racc := rcff*pr - icff*pi; iacc := rcff*pi + icff*pr;
      racc := powfac*racc; iacc := powfac*iacc;
    end if;
    hrprow := hrp(fac1); hiprow := hip(fac1);
    if m1 = 2 then
     -- H(fac(1),fac(1)) := H(fac(1),fac(1)) + acc;
      hrprow(fac1) := hrprow(fac1) + racc;
      hiprow(fac1) := hiprow(fac1) + iacc;
    elsif m1 = 3 then
     -- H(fac(1),fac(1)) := H(fac(1),fac(1)) + acc*x(fac(1));
      pr := xr(fac1); pi := xi(fac1);
      zr := racc*pr - iacc*pi; zi := racc*pi + iacc*pr;
      hrprow(fac1) := hrprow(fac1) + zr;
      hiprow(fac1) := hiprow(fac1) + zi;
    else
     -- H(fac(1),fac(1)) := H(fac(1),fac(1)) + acc*(pwt(fac(1))(m1-3));
      pr := rpwt(fac1)(m1-3); pi := ipwt(fac1)(m1-3);
      zr := racc*pr - iacc*pi; zi := racc*pi + iacc*pr;
      hrprow(fac1) := hrprow(fac1) + zr;
      hiprow(fac1) := hiprow(fac1) + zi;
    end if;
    if fac'last = 1 then
      powfac := double_float(m1);
      hrprow := hrp(idx1); hiprow := hip(idx1);
      if m1 = 2 then
       -- H(idx(1),idx(2)) := H(idx(1),idx(2)) + c*powfac*x(fac(1));
        pr := xr(fac1); pi := xi(fac1);
        zr := rcff*pr - icff*pi; zi := rcff*pi + icff*pr;
        hrprow(idx2) := hrprow(idx2) + powfac*zr;
        hiprow(idx2) := hiprow(idx2) + powfac*zi;
      else
       -- H(idx(1),idx(2)) := H(idx(1),idx(2)) + c*powfac*(pwt(fac(1))(m1-2));
        pr := rpwt(fac1)(m1-2); pi := ipwt(fac1)(m1-2);
        zr := rcff*pr - icff*pi; zi := rcff*pi + icff*pr;
        hrprow(idx2) := hrprow(idx2) + powfac*zr;
        hiprow(idx2) := hiprow(idx2) + powfac*zi;
      end if;
    else
      powfac := double_float(m2*(m2-1));
     -- acc := c*powfac*(pwt(fac(1))(m1-1));
      pr := rpwt(fac1)(m1-1); pi := ipwt(fac1)(m1-1);
      racc := rcff*pr - icff*pi; iacc := rcff*pi + icff*pr;
      racc := powfac*racc; iacc := powfac*iacc;
      hrprow := hrp(fac2); hiprow := hip(fac2);
      if m2 = 2 then
       -- H(fac(2),fac(2)) := H(fac(2),fac(2)) + acc;
        hrprow(fac2) := hrprow(fac2) + racc;
        hiprow(fac2) := hiprow(fac2) + iacc;
      elsif m2 = 3 then
       -- H(fac(2),fac(2)) := H(fac(2),fac(2)) + acc*x(fac(2));
        pr := xr(fac2); pi := xi(fac2);
        zr := racc*pr - iacc*pi; zi := racc*pi + iacc*pr;
        hrprow(fac2) := hrprow(fac2) + zr;
        hiprow(fac2) := hiprow(fac2) + zi;
      else
       -- H(fac(2),fac(2)) := H(fac(2),fac(2)) + acc*(pwt(fac(2))(m2-3));
        pr := rpwt(fac2)(m2-3); pi := ipwt(fac2)(m2-3);
        zr := racc*pr - iacc*pi; zi := racc*pi + iacc*pr;
        hrprow(fac2) := hrprow(fac2) + zr;
        hiprow(fac2) := hiprow(fac2) + zi;
      end if;
      powfac := double_float(m1*m2);
     -- acc := c*powfac;
      racc := powfac*rcff; iacc := powfac*icff;
      if m1 = 2 then
       -- acc := acc*x(fac(1));
        pr := xr(fac1); pi := xi(fac1);
        zr := racc*pr - iacc*pi; zi := racc*pi + iacc*pr;
        racc := zr; iacc := zi;
      else
       -- acc := acc*(pwt(fac(1))(m1-2));
        pr := rpwt(fac1)(m1-2); pi := ipwt(fac1)(m1-2);
        zr := racc*pr - iacc*pi; zi := racc*pi + iacc*pr;
        racc := zr; iacc := zi;
      end if;
      if m2 = 2 then
       -- acc := acc*x(fac(2));
        pr := xr(fac2); pi := xi(fac2);
        zr := racc*pr - iacc*pi; zi := racc*pi + iacc*pr;
        racc := zr; iacc := zi;
      else
       -- acc := acc*(pwt(fac(2))(m2-2));
        pr := rpwt(fac2)(m2-2); pi := ipwt(fac2)(m2-2);
        zr := racc*pr - iacc*pi; zi := racc*pi + iacc*pr;
        racc := zr; iacc := zi;
      end if;
     -- H(idx(1),idx(2)) := H(idx(1),idx(2)) + acc;
      hrprow := hrp(idx1); hiprow := hip(idx1);
      hrprow(idx2) := hrprow(idx2) + racc;
      hiprow(idx2) := hiprow(idx2) + iacc;
    end if;
  end Speel2;

  procedure Speel3 ( H : in out Standard_Complex_Matrices.Matrix;
                     c : in Complex_Number;
                     xps : in Standard_Integer_Vectors.Vector;
                     idx : in Standard_Integer_Vectors.Vector;
                     fac : in Standard_Integer_Vectors.Vector;
                     x : in Standard_Complex_Vectors.Vector;
                     pwt : in Standard_Complex_VecVecs.VecVec ) is

    m1,m2 : integer32;
    powfac : double_float; -- multiplier factor of two powers
    acc : Complex_Number;
    offdiagfac : Complex_Number; -- common off diagonal factor
    ondiagfac : Complex_Number;  -- common on diagonal factor

  begin
    offdiagfac := c; ondiagfac := c;  -- off/on diagonal factors
    for i in fac'range loop
      m1 := xps(fac(i));
      if m1 = 2 then
        offdiagfac := offdiagfac*x(fac(i));
      elsif m1 = 3 then
        offdiagfac := offdiagfac*(pwt(fac(i))(1));
        ondiagfac := ondiagfac*x(fac(i));
      else
        offdiagfac := offdiagfac*(pwt(fac(i))(m1-2));
        ondiagfac := ondiagfac*(pwt(fac(i))(m1-3));
      end if;
    end loop;
   -- compute the off diagonal elements of the Hessian
    powfac := double_float(xps(idx(1))*xps(idx(2)));
    H(idx(1),idx(2)) := H(idx(1),idx(2)) + offdiagfac*powfac*x(idx(3));
    powfac := double_float(xps(idx(1))*xps(idx(3)));
    H(idx(1),idx(3)) := H(idx(1),idx(3)) + offdiagfac*powfac*x(idx(2));
    powfac := double_float(xps(idx(2))*xps(idx(3)));
    H(idx(2),idx(3)) := H(idx(2),idx(3)) + offdiagfac*powfac*x(idx(1));
   -- ten cases for the on diagonal element of the Hessian
    if fac'last = 3 then -- all variables raised to higher power
      m1 := xps(fac(1)); powfac := double_float(m1*(m1-1));
      acc := ondiagfac*(pwt(fac(2))(1))*(pwt(fac(3))(1));
      H(fac(1),fac(1)) := H(fac(1),fac(1)) + acc*powfac;
      m1 := xps(fac(2)); powfac := double_float(m1*(m1-1));
      acc := ondiagfac*(pwt(fac(1))(1))*(pwt(fac(3))(1));
      H(fac(2),fac(2)) := H(fac(2),fac(2)) + acc*powfac;
      m1 := xps(fac(3)); powfac := double_float(m1*(m1-1));
      acc := ondiagfac*(pwt(fac(1))(1))*(pwt(fac(2))(1));
      H(fac(3),fac(3)) := H(fac(3),fac(3)) + acc*powfac;
    elsif fac'last = 1 then -- one variable raised to higher power
      m1 := xps(fac(1)); powfac := double_float(m1*(m1-1));
      if fac(1) = idx(1) then
        acc := ondiagfac*x(idx(2))*x(idx(3));
      elsif fac(1) = idx(2) then
        acc := ondiagfac*x(idx(1))*x(idx(3));
      else -- fac(1) = idx(3)
        acc := ondiagfac*x(idx(1))*x(idx(2));
      end if;
      H(fac(1),fac(1)) := H(fac(1),fac(1)) + acc*powfac;
    else -- fac'last = 2, two variables raised to higher power
      m1 := xps(fac(1)); powfac := double_float(m1*(m1-1));
      m2 := xps(fac(2));
      if fac(1) = idx(1) then
        if fac(2) = idx(2) then -- idx(3) has power 1
          acc := ondiagfac*(pwt(fac(2))(1))*x(idx(3));
          H(fac(1),fac(1)) := H(fac(1),fac(1)) + acc*powfac;
          powfac := double_float(m2*(m2-1));
          acc := ondiagfac*(pwt(fac(1))(1))*x(idx(3));
          H(fac(2),fac(2)) := H(fac(2),fac(2)) + acc*powfac;
        else -- idx(2) has power 1
          acc := ondiagfac*x(idx(2))*(pwt(fac(2))(1));
          H(fac(1),fac(1)) := H(fac(1),fac(1)) + acc*powfac;
          powfac := double_float(m2*(m2-1));
          acc := ondiagfac*x(idx(2))*(pwt(fac(1))(1));
          H(fac(2),fac(2)) := H(fac(2),fac(2)) + acc*powfac;
        end if;
      elsif fac(1) = idx(2) then
        if fac(2) = idx(1) then -- idx(3) has power 1
          acc := ondiagfac*(pwt(fac(2))(1))*x(idx(3));
          H(fac(1),fac(1)) := H(fac(1),fac(1)) + acc*powfac;
          powfac := double_float(m2*(m2-1));
          acc := ondiagfac*(pwt(fac(1))(1))*x(idx(3));
          H(fac(2),fac(2)) := H(fac(2),fac(2)) + acc*powfac;
        else -- idx(1) has power 1
          acc := ondiagfac*x(idx(1))*(pwt(fac(2))(1));
          H(fac(1),fac(1)) := H(fac(1),fac(1)) + acc*powfac;
          powfac := double_float(m2*(m2-1));
          acc := ondiagfac*x(idx(1))*(pwt(fac(1))(1));
          H(fac(2),fac(2)) := H(fac(2),fac(2)) + acc*powfac;
        end if;
      else --  fac(1) = idx(3)
        if fac(2) = idx(1) then -- idx(2) has power 1
          acc := ondiagfac*(pwt(fac(2))(1))*x(idx(2));
          H(fac(1),fac(1)) := H(fac(1),fac(1)) + acc*powfac;
          powfac := double_float(m2*(m2-1));
          acc := ondiagfac*(pwt(fac(1))(1))*x(idx(2));
          H(fac(2),fac(2)) := H(fac(2),fac(2)) + acc*powfac;
        else -- idx(1) has power 1
          acc := ondiagfac*x(idx(1))*(pwt(fac(2))(1));
          H(fac(1),fac(1)) := H(fac(1),fac(1)) + acc*powfac;
          powfac := double_float(m2*(m2-1));
          acc := ondiagfac*x(idx(1))*(pwt(fac(1))(1));
          H(fac(2),fac(2)) := H(fac(2),fac(2)) + acc*powfac;
        end if;
      end if;
    end if;
  end Speel3;

  procedure Speel3 ( hrp : in Standard_Floating_VecVecs.VecVec;
                     hip : in Standard_Floating_VecVecs.VecVec;
                     rcff,icff : in double_float;
                     xps : in Standard_Integer_Vectors.Vector;
                     idx : in Standard_Integer_Vectors.Vector;
                     fac : in Standard_Integer_Vectors.Vector;
                     xr : in Standard_Floating_Vectors.Vector;
                     xi : in Standard_Floating_Vectors.Vector;
                     rpwt : in Standard_Floating_VecVecs.VecVec;
                     ipwt : in Standard_Floating_VecVecs.VecVec ) is

    m1,m2,fidx1,fidx2,fidx3 : integer32;
    powfac : double_float; -- multiplier factor of two powers
    pr,pi,zr,zi,racc,iacc : double_float;
    roffdiagfac,ioffdiagfac : double_float; -- common off diagonal factor
    rondiagfac,iondiagfac : double_float;   -- common on diagonal factor
    idx1 : constant integer32 := idx(1);
    idx2 : constant integer32 := idx(2);
    idx3 : constant integer32 := idx(3);
    hrprow,hiprow : Standard_Floating_Vectors.Link_to_Vector;

  begin
   -- offdiagfac := c; ondiagfac := c;  -- off/on diagonal factors
    roffdiagfac := rcff; ioffdiagfac := icff;
    rondiagfac := rcff;  iondiagfac := icff;
    for i in fac'range loop
      fidx1 := fac(i); m1 := xps(fidx1);
      if m1 = 2 then
       -- offdiagfac := offdiagfac*x(fac(i));
        pr := xr(fidx1); pi := xi(fidx1);
        zr := roffdiagfac*pr - ioffdiagfac*pi;
        zi := roffdiagfac*pi + ioffdiagfac*pr;
        roffdiagfac := zr; ioffdiagfac := zi;
      elsif m1 = 3 then
       -- offdiagfac := offdiagfac*(pwt(fac(i))(1));
        pr := rpwt(fidx1)(1); pi := ipwt(fidx1)(1);
        zr := roffdiagfac*pr - ioffdiagfac*pi;
        zi := roffdiagfac*pi + ioffdiagfac*pr;
        roffdiagfac := zr; ioffdiagfac := zi;
       -- ondiagfac := ondiagfac*x(fac(i));
        pr := xr(fidx1); pi := xi(fidx1);
        zr := rondiagfac*pr - iondiagfac*pi;
        zi := rondiagfac*pi + iondiagfac*pr;
        rondiagfac := zr; iondiagfac := zi;
      else
       -- offdiagfac := offdiagfac*(pwt(fac(i))(m1-2));
        pr := rpwt(fidx1)(m1-2); pi := ipwt(fidx1)(m1-2);
        zr := roffdiagfac*pr - ioffdiagfac*pi;
        zi := roffdiagfac*pi + ioffdiagfac*pr;
        roffdiagfac := zr; ioffdiagfac := zi;
       -- ondiagfac := ondiagfac*(pwt(fac(i))(m1-3));
        pr := rpwt(fidx1)(m1-3); pi := ipwt(fidx1)(m1-3);
        zr := rondiagfac*pr - iondiagfac*pi;
        zi := rondiagfac*pi + iondiagfac*pr;
        rondiagfac := zr; iondiagfac := zi;
      end if;
    end loop;
   -- compute the off diagonal elements of the Hessian
    powfac := double_float(xps(idx1)*xps(idx2));
    hrprow := hrp(idx1); hiprow := hip(idx1); 
   -- H(idx(1),idx(2)) := H(idx(1),idx(2)) + offdiagfac*powfac*x(idx(3));
    pr := xr(idx3); pi := xi(idx3);
    zr := roffdiagfac*pr - ioffdiagfac*pi;
    zi := roffdiagfac*pi + ioffdiagfac*pr;
    hrprow(idx2) := hrprow(idx2) + powfac*zr;
    hiprow(idx2) := hiprow(idx2) + powfac*zi;
    powfac := double_float(xps(idx1)*xps(idx3));
   -- H(idx(1),idx(3)) := H(idx(1),idx(3)) + offdiagfac*powfac*x(idx(2));
    pr := xr(idx2); pi := xi(idx2);
    zr := roffdiagfac*pr - ioffdiagfac*pi;
    zi := roffdiagfac*pi + ioffdiagfac*pr;
    hrprow(idx3) := hrprow(idx3) + powfac*zr;
    hiprow(idx3) := hiprow(idx3) + powfac*zi;
    powfac := double_float(xps(idx2)*xps(idx3));
   -- H(idx(2),idx(3)) := H(idx(2),idx(3)) + offdiagfac*powfac*x(idx(1));
    hrprow := hrp(idx2); hiprow := hip(idx2);
    pr := xr(idx1); pi := xi(idx1);
    zr := roffdiagfac*pr - ioffdiagfac*pi;
    zi := roffdiagfac*pi + ioffdiagfac*pr;
    hrprow(idx3) := hrprow(idx3) + powfac*zr;
    hiprow(idx3) := hiprow(idx3) + powfac*zi;
   -- ten cases for the on diagonal elements of the Hessian
    if fac'last = 3 then -- all variables raised to higher power
      fidx1 := fac(1); fidx2 := fac(2); fidx3 := fac(3);
      m1 := xps(fidx1); powfac := double_float(m1*(m1-1));
     -- acc := ondiagfac*(pwt(fac(2))(1))*(pwt(fac(3))(1));
      pr := rpwt(fidx2)(1); pi := ipwt(fidx2)(1);
      zr := rondiagfac*pr - iondiagfac*pi;
      zi := rondiagfac*pi + iondiagfac*pr;
      pr := rpwt(fidx3)(1); pi := ipwt(fidx3)(1);
      racc := zr*pr - zi*pi; iacc := zr*pi + zi*pr;
     -- H(fac(1),fac(1)) := H(fac(1),fac(1)) + acc*powfac;
      hrprow := hrp(fidx1); hiprow := hip(fidx1);
      hrprow(fidx1) := hrprow(fidx1) + racc*powfac;
      hiprow(fidx1) := hiprow(fidx1) + iacc*powfac;
      m1 := xps(fidx2); powfac := double_float(m1*(m1-1));
     -- acc := ondiagfac*(pwt(fac(1))(1))*(pwt(fac(3))(1));
      pr := rpwt(fidx1)(1); pi := ipwt(fidx1)(1);
      zr := rondiagfac*pr - iondiagfac*pi;
      zi := rondiagfac*pi + iondiagfac*pr;
      pr := rpwt(fidx3)(1); pi := ipwt(fidx3)(1);
      racc := zr*pr - zi*pi; iacc := zr*pi + zi*pr;
     -- H(fac(2),fac(2)) := H(fac(2),fac(2)) + acc*powfac;
      hrprow := hrp(fidx2); hiprow := hip(fidx2);
      hrprow(fidx2) := hrprow(fidx2) + racc*powfac;
      hiprow(fidx2) := hiprow(fidx2) + iacc*powfac;
      m1 := xps(fidx3); powfac := double_float(m1*(m1-1));
     -- acc := ondiagfac*(pwt(fac(1))(1))*(pwt(fac(2))(1));
      pr := rpwt(fidx1)(1); pi := ipwt(fidx1)(1);
      zr := rondiagfac*pr - iondiagfac*pi;
      zi := rondiagfac*pi + iondiagfac*pr;
      pr := rpwt(fidx2)(1); pi := ipwt(fidx2)(1);
      racc := zr*pr - zi*pi; iacc := zr*pi + zi*pr;
     -- H(fac(3),fac(3)) := H(fac(3),fac(3)) + acc*powfac;
      hrprow := hrp(fidx3); hiprow := hip(fidx3);
      hrprow(fidx3) := hrprow(fidx3) + racc*powfac;
      hiprow(fidx3) := hiprow(fidx3) + iacc*powfac;
    elsif fac'last = 1 then -- one variable raised to higher power
      fidx1 := fac(1);
      m1 := xps(fidx1); powfac := double_float(m1*(m1-1));
      if fidx1 = idx1 then
       -- acc := ondiagfac*x(idx(2))*x(idx(3));
        pr := xr(idx2); pi := xi(idx2);
        zr := rondiagfac*pr - iondiagfac*pi;
        zi := rondiagfac*pi + iondiagfac*pr;
        pr := xr(idx3); pi := xi(idx3);
        racc := zr*pr - zi*pi; iacc := zr*pi + zi*pr;
      elsif fidx1 = idx2 then
       -- acc := ondiagfac*x(idx(1))*x(idx(3));
        pr := xr(idx1); pi := xi(idx1);
        zr := rondiagfac*pr - iondiagfac*pi;
        zi := rondiagfac*pi + iondiagfac*pr;
        pr := xr(idx3); pi := xi(idx3);
        racc := zr*pr - zi*pi; iacc := zr*pi + zi*pr;
      else -- fac(1) = idx(3)
       -- acc := ondiagfac*x(idx(1))*x(idx(2));
        pr := xr(idx1); pi := xi(idx1);
        zr := rondiagfac*pr - iondiagfac*pi;
        zi := rondiagfac*pi + iondiagfac*pr;
        pr := xr(idx2); pi := xi(idx2);
        racc := zr*pr - zi*pi; iacc := zr*pi + zi*pr;
      end if;
     -- H(fac(1),fac(1)) := H(fac(1),fac(1)) + acc*powfac;
      hrprow := hrp(fidx1); hiprow := hip(fidx1);
      hrprow(fidx1) := hrprow(fidx1) + racc*powfac;
      hiprow(fidx1) := hiprow(fidx1) + iacc*powfac;
    else -- fac'last = 2, two variables raised to higher power
      fidx1 := fac(1); fidx2 := fac(2);
      m1 := xps(fidx1); powfac := double_float(m1*(m1-1));
      m2 := xps(fidx2);
      if fidx1 = idx1 then
        if fidx2 = idx2 then -- idx(3) has power 1
         -- acc := ondiagfac*(pwt(fac(2))(1))*x(idx(3));
          pr := rpwt(fidx2)(1); pi := ipwt(fidx2)(1); 
          zr := rondiagfac*pr - iondiagfac*pi;
          zi := rondiagfac*pi + iondiagfac*pr;
          pr := xr(idx3); pi := xi(idx3); 
          racc := zr*pr - zi*pi; iacc := zr*pi + zi*pr;
         -- H(fac(1),fac(1)) := H(fac(1),fac(1)) + acc*powfac;
          hrprow := hrp(fidx1); hiprow := hip(fidx1);
          hrprow(fidx1) := hrprow(fidx1) + racc*powfac;
          hiprow(fidx1) := hiprow(fidx1) + iacc*powfac;
          powfac := double_float(m2*(m2-1));
         -- acc := ondiagfac*(pwt(fac(1))(1))*x(idx(3));
          pr := rpwt(fidx1)(1); pi := ipwt(fidx1)(1); 
          zr := rondiagfac*pr - iondiagfac*pi;
          zi := rondiagfac*pi + iondiagfac*pr;
          pr := xr(idx3); pi := xi(idx3); 
          racc := zr*pr - zi*pi; iacc := zr*pi + zi*pr;
         -- H(fac(2),fac(2)) := H(fac(2),fac(2)) + acc*powfac;
          hrprow := hrp(fidx2); hiprow := hip(fidx2);
          hrprow(fidx2) := hrprow(fidx2) + racc*powfac;
          hiprow(fidx2) := hiprow(fidx2) + iacc*powfac;
        else -- idx(2) has power 1
         -- acc := ondiagfac*x(idx(2))*(pwt(fac(2))(1));
          pr := rpwt(fidx2)(1); pi := ipwt(fidx2)(1); 
          zr := rondiagfac*pr - iondiagfac*pi;
          zi := rondiagfac*pi + iondiagfac*pr;
          pr := xr(idx2); pi := xi(idx2); 
          racc := zr*pr - zi*pi; iacc := zr*pi + zi*pr;
         -- H(fac(1),fac(1)) := H(fac(1),fac(1)) + acc*powfac;
          hrprow := hrp(fidx1); hiprow := hip(fidx1);
          hrprow(fidx1) := hrprow(fidx1) + racc*powfac;
          hiprow(fidx1) := hiprow(fidx1) + iacc*powfac;
          powfac := double_float(m2*(m2-1));
         -- acc := ondiagfac*x(idx(2))*(pwt(fac(1))(1));
          pr := rpwt(fidx1)(1); pi := ipwt(fidx1)(1); 
          zr := rondiagfac*pr - iondiagfac*pi;
          zi := rondiagfac*pi + iondiagfac*pr;
          pr := xr(idx2); pi := xi(idx2); 
          racc := zr*pr - zi*pi; iacc := zr*pi + zi*pr;
         -- H(fac(2),fac(2)) := H(fac(2),fac(2)) + acc*powfac;
          hrprow := hrp(fidx2); hiprow := hip(fidx2);
          hrprow(fidx2) := hrprow(fidx2) + racc*powfac;
          hiprow(fidx2) := hiprow(fidx2) + iacc*powfac;
        end if;
      elsif fidx1 = idx2 then
        if fidx2 = idx1 then -- idx(3) has power 1
         -- acc := ondiagfac*(pwt(fac(2))(1))*x(idx(3));
          pr := rpwt(fidx2)(1); pi := ipwt(fidx2)(1); 
          zr := rondiagfac*pr - iondiagfac*pi;
          zi := rondiagfac*pi + iondiagfac*pr;
          pr := xr(idx3); pi := xi(idx3); 
          racc := zr*pr - zi*pi; iacc := zr*pi + zi*pr;
         -- H(fac(1),fac(1)) := H(fac(1),fac(1)) + acc*powfac;
          hrprow := hrp(fidx1); hiprow := hip(fidx1);
          hrprow(fidx1) := hrprow(fidx1) + racc*powfac;
          hiprow(fidx1) := hiprow(fidx1) + iacc*powfac;
          powfac := double_float(m2*(m2-1));
         -- acc := ondiagfac*(pwt(fac(1))(1))*x(idx(3));
          pr := rpwt(fidx1)(1); pi := ipwt(fidx1)(1); 
          zr := rondiagfac*pr - iondiagfac*pi;
          zi := rondiagfac*pi + iondiagfac*pr;
          pr := xr(idx3); pi := xi(idx3); 
          racc := zr*pr - zi*pi; iacc := zr*pi + zi*pr;
         -- H(fac(2),fac(2)) := H(fac(2),fac(2)) + acc*powfac;
          hrprow := hrp(fidx2); hiprow := hip(fidx2);
          hrprow(fidx2) := hrprow(fidx2) + racc*powfac;
          hiprow(fidx2) := hiprow(fidx2) + iacc*powfac;
        else -- idx(1) has power 1
         -- acc := ondiagfac*x(idx(1))*(pwt(fac(2))(1));
          pr := rpwt(fidx2)(1); pi := ipwt(fidx2)(1); 
          zr := rondiagfac*pr - iondiagfac*pi;
          zi := rondiagfac*pi + iondiagfac*pr;
          pr := xr(idx1); pi := xi(idx1); 
          racc := zr*pr - zi*pi; iacc := zr*pi + zi*pr;
         -- H(fac(1),fac(1)) := H(fac(1),fac(1)) + acc*powfac;
          hrprow := hrp(fidx1); hiprow := hip(fidx1);
          hrprow(fidx1) := hrprow(fidx1) + racc*powfac;
          hiprow(fidx1) := hiprow(fidx1) + iacc*powfac;
          powfac := double_float(m2*(m2-1));
         -- acc := ondiagfac*x(idx(1))*(pwt(fac(1))(1));
          pr := rpwt(fidx1)(1); pi := ipwt(fidx1)(1); 
          zr := rondiagfac*pr - iondiagfac*pi;
          zi := rondiagfac*pi + iondiagfac*pr;
          pr := xr(idx1); pi := xi(idx1); 
          racc := zr*pr - zi*pi; iacc := zr*pi + zi*pr;
         -- H(fac(2),fac(2)) := H(fac(2),fac(2)) + acc*powfac;
          hrprow := hrp(fidx2); hiprow := hip(fidx2);
          hrprow(fidx2) := hrprow(fidx2) + racc*powfac;
          hiprow(fidx2) := hiprow(fidx2) + iacc*powfac;
        end if;
      else --  fac(1) = idx(3)
        if fidx2 = idx1 then -- idx(2) has power 1
         -- acc := ondiagfac*(pwt(fac(2))(1))*x(idx(2));
          pr := rpwt(fidx2)(1); pi := ipwt(fidx2)(1); 
          zr := rondiagfac*pr - iondiagfac*pi;
          zi := rondiagfac*pi + iondiagfac*pr;
          pr := xr(idx2); pi := xi(idx2); 
          racc := zr*pr - zi*pi; iacc := zr*pi + zi*pr;
         -- H(fac(1),fac(1)) := H(fac(1),fac(1)) + acc*powfac;
          hrprow := hrp(fidx1); hiprow := hip(fidx1);
          hrprow(fidx1) := hrprow(fidx1) + racc*powfac;
          hiprow(fidx1) := hiprow(fidx1) + iacc*powfac;
          powfac := double_float(m2*(m2-1));
         -- acc := ondiagfac*(pwt(fac(1))(1))*x(idx(2));
          pr := rpwt(fidx1)(1); pi := ipwt(fidx1)(1); 
          zr := rondiagfac*pr - iondiagfac*pi;
          zi := rondiagfac*pi + iondiagfac*pr;
          pr := xr(idx2); pi := xi(idx2); 
          racc := zr*pr - zi*pi; iacc := zr*pi + zi*pr;
         -- H(fac(2),fac(2)) := H(fac(2),fac(2)) + acc*powfac;
          hrprow := hrp(fidx2); hiprow := hip(fidx2);
          hrprow(fidx2) := hrprow(fidx2) + racc*powfac;
          hiprow(fidx2) := hiprow(fidx2) + iacc*powfac;
        else -- idx(1) has power 1
         -- acc := ondiagfac*x(idx(1))*(pwt(fac(2))(1));
          pr := rpwt(fidx2)(1); pi := ipwt(fidx2)(1); 
          zr := rondiagfac*pr - iondiagfac*pi;
          zi := rondiagfac*pi + iondiagfac*pr;
          pr := xr(idx1); pi := xi(idx1); 
          racc := zr*pr - zi*pi; iacc := zr*pi + zi*pr;
         -- H(fac(1),fac(1)) := H(fac(1),fac(1)) + acc*powfac;
          hrprow := hrp(fidx1); hiprow := hip(fidx1);
          hrprow(fidx1) := hrprow(fidx1) + racc*powfac;
          hiprow(fidx1) := hiprow(fidx1) + iacc*powfac;
          powfac := double_float(m2*(m2-1));
         -- acc := ondiagfac*x(idx(1))*(pwt(fac(1))(1));
          pr := rpwt(fidx1)(1); pi := ipwt(fidx1)(1); 
          zr := rondiagfac*pr - iondiagfac*pi;
          zi := rondiagfac*pi + iondiagfac*pr;
          pr := xr(idx1); pi := xi(idx1); 
          racc := zr*pr - zi*pi; iacc := zr*pi + zi*pr;
         -- H(fac(2),fac(2)) := H(fac(2),fac(2)) + acc*powfac;
          hrprow := hrp(fidx2); hiprow := hip(fidx2);
          hrprow(fidx2) := hrprow(fidx2) + racc*powfac;
          hiprow(fidx2) := hiprow(fidx2) + iacc*powfac;
        end if;
      end if;
    end if;
  end Speel3;

  procedure Speel4 ( H : in out Standard_Complex_Matrices.Matrix;
                     c : in Complex_Number;
                     xps : in Standard_Integer_Vectors.Vector;
                     idx : in Standard_Integer_Vectors.Vector;
                     fac : in Standard_Integer_Vectors.Vector;
                     x : in Standard_Complex_Vectors.Vector;
                     pwt : in Standard_Complex_VecVecs.VecVec ) is

    m1 : integer32;
    powfac : double_float; -- multiplier factor of two powers
    acc : Complex_Number;
    offdiagfac : Complex_Number; -- common off diagonal factor
    ondiagfac : Complex_Number;  -- common on diagonal factor
    fwd1 : Complex_Number; -- first forward product
    bck1 : Complex_Number; -- first backward prodcut

  begin
    offdiagfac := c; ondiagfac := c;  -- off/on diagonal factors
    for i in fac'range loop
      m1 := xps(fac(i));
      if m1 = 2 then
        offdiagfac := offdiagfac*x(fac(i));
      elsif m1 = 3 then
        offdiagfac := offdiagfac*(pwt(fac(i))(1));
        ondiagfac := ondiagfac*x(fac(i));
      else
        offdiagfac := offdiagfac*(pwt(fac(i))(m1-2));
        ondiagfac := ondiagfac*(pwt(fac(i))(m1-3));
      end if;
    end loop;
   -- the off diagonal elements use forward and backward products
    fwd1 := x(idx(1))*x(idx(2));
    bck1 := x(idx(4))*x(idx(3));
   -- the last element is a copy of fwd1, with a multiplier factor
    powfac := double_float(xps(idx(3))*xps(idx(4)));
    H(idx(3),idx(4)) := H(idx(3),idx(4)) + offdiagfac*powfac*fwd1;
   -- the first element is a copy of bck(1), with a multiplier factor
    powfac := double_float(xps(idx(1))*xps(idx(2)));
    H(idx(1),idx(2)) := H(idx(1),idx(2)) + offdiagfac*powfac*bck1;
   -- the other off diagonal elements
    acc := offdiagfac*x(idx(2));
    powfac := double_float(xps(idx(1))*xps(idx(3)));
    H(idx(1),idx(3)) := H(idx(1),idx(3)) + acc*powfac*x(idx(4));
    powfac := double_float(xps(idx(1))*xps(idx(4)));
    H(idx(1),idx(4)) := H(idx(1),idx(4)) + acc*powfac*x(idx(3));
    acc := offdiagfac*x(idx(1));
    powfac := double_float(xps(idx(2))*xps(idx(3)));
    H(idx(2),idx(3)) := H(idx(2),idx(3)) + acc*powfac*x(idx(4));
    powfac := double_float(xps(idx(2))*xps(idx(4)));
    H(idx(2),idx(4)) := H(idx(2),idx(4)) + acc*powfac*x(idx(3));
   -- compute the diagonal elements
    for k in fac'range loop
      m1 := xps(fac(k)); powfac := double_float(m1*(m1-1));
      acc := powfac*ondiagfac; -- acc is the cofactor
      for i in idx'range loop
        if idx(i) /= fac(k) then -- skip the current factor
          if xps(idx(i)) = 1 
           then acc := acc*x(idx(i));
           else acc := acc*(pwt(idx(i))(1));
          end if;
        end if;
      end loop;
      H(fac(k),fac(k)) := H(fac(k),fac(k)) + acc;
    end loop;
  end Speel4;

  procedure Speel4 ( hrp : in Standard_Floating_VecVecs.VecVec;
                     hip : in Standard_Floating_VecVecs.VecVec;
                     rcff,icff : in double_float;
                     xps : in Standard_Integer_Vectors.Vector;
                     idx : in Standard_Integer_Vectors.Vector;
                     fac : in Standard_Integer_Vectors.Vector;
                     xr : in Standard_Floating_Vectors.Vector;
                     xi : in Standard_Floating_Vectors.Vector;
                     rpwt : in Standard_Floating_VecVecs.VecVec;
                     ipwt : in Standard_Floating_VecVecs.VecVec ) is

    m1,fidx1 : integer32;
    powfac : double_float; -- multiplier factor of two powers
    pr,pi,zr,zi,racc,iacc : double_float;
    roffdiagfac,ioffdiagfac : double_float; -- common off diagonal factor
    rondiagfac,iondiagfac : double_float;   -- common on diagonal factor
    rfwd1,ifwd1,rbck1,ibck1 : double_float; -- forward & backward products
    idx1 : constant integer32 := idx(1);
    idx2 : constant integer32 := idx(2);
    idx3 : constant integer32 := idx(3);
    idx4 : constant integer32 := idx(4);
    hrprow,hiprow : Standard_Floating_Vectors.Link_to_Vector;

  begin
   -- offdiagfac := c; ondiagfac := c;  -- off/on diagonal factors
    roffdiagfac := rcff; ioffdiagfac := icff;
    rondiagfac := rcff;  iondiagfac := icff;
    for i in fac'range loop
      fidx1 := fac(i); m1 := xps(fidx1);
      if m1 = 2 then
       -- offdiagfac := offdiagfac*x(fac(i));
        pr := xr(fidx1); pi := xi(fidx1);
        zr := roffdiagfac*pr - ioffdiagfac*pi;
        zi := roffdiagfac*pi + ioffdiagfac*pr;
        roffdiagfac := zr; ioffdiagfac := zi;
      elsif m1 = 3 then
       -- offdiagfac := offdiagfac*(pwt(fac(i))(1));
        pr := rpwt(fidx1)(1); pi := ipwt(fidx1)(1);
        zr := roffdiagfac*pr - ioffdiagfac*pi;
        zi := roffdiagfac*pi + ioffdiagfac*pr;
        roffdiagfac := zr; ioffdiagfac := zi;
       -- ondiagfac := ondiagfac*x(fac(i));
        pr := xr(fidx1); pi := xi(fidx1);
        zr := rondiagfac*pr - iondiagfac*pi;
        zi := rondiagfac*pi + iondiagfac*pr;
        rondiagfac := zr; iondiagfac := zi;
      else
       -- offdiagfac := offdiagfac*(pwt(fac(i))(m1-2));
        pr := rpwt(fidx1)(m1-2); pi := ipwt(fidx1)(m1-2);
        zr := roffdiagfac*pr - ioffdiagfac*pi;
        zi := roffdiagfac*pi + ioffdiagfac*pr;
        roffdiagfac := zr; ioffdiagfac := zi;
       -- ondiagfac := ondiagfac*(pwt(fac(i))(m1-3));
        pr := rpwt(fidx1)(m1-3); pi := ipwt(fidx1)(m1-3);
        zr := rondiagfac*pr - iondiagfac*pi;
        zi := rondiagfac*pi + iondiagfac*pr;
        rondiagfac := zr; iondiagfac := zi;
      end if;
    end loop;
   -- the off diagonal elements use forward and backward products
   -- fwd1 := x(idx(1))*x(idx(2));
    pr := xr(idx1); pi := xi(idx1); zr := xr(idx2); zi := xi(idx2);
    rfwd1 := pr*zr - pi*zi; ifwd1 := pr*zi + pi*zr;
   -- bck1 := x(idx(4))*x(idx(3));
    pr := xr(idx4); pi := xi(idx4); zr := xr(idx3); zi := xi(idx3);
    rbck1 := pr*zr - pi*zi; ibck1 := pr*zi + pi*zr;
   -- the last element is a copy of fwd1, with a multiplier factor
    powfac := double_float(xps(idx3)*xps(idx4));
   -- H(idx(3),idx(4)) := H(idx(3),idx(4)) + offdiagfac*powfac*fwd1;
    hrprow := hrp(idx3); hiprow := hip(idx3);
    zr := roffdiagfac*rfwd1 - ioffdiagfac*ifwd1;
    zi := roffdiagfac*ifwd1 + ioffdiagfac*rfwd1;
    hrprow(idx4) := hrprow(idx4) + powfac*zr;
    hiprow(idx4) := hiprow(idx4) + powfac*zi;
   -- the first element is a copy of bck(1), with a multiplier factor
    powfac := double_float(xps(idx1)*xps(idx2));
   -- H(idx(1),idx(2)) := H(idx(1),idx(2)) + offdiagfac*powfac*bck1;
    hrprow := hrp(idx1); hiprow := hip(idx1);
    zr := roffdiagfac*rbck1 - ioffdiagfac*ibck1;
    zi := roffdiagfac*ibck1 + ioffdiagfac*rbck1;
    hrprow(idx2) := hrprow(idx2) + powfac*zr;
    hiprow(idx2) := hiprow(idx2) + powfac*zi;
   -- the other off diagonal elements
   -- acc := offdiagfac*x(idx(2));
    pr := xr(idx2); pi := xi(idx2);
    racc := roffdiagfac*pr - ioffdiagfac*pi;
    iacc := roffdiagfac*pi + ioffdiagfac*pr;
    powfac := double_float(xps(idx1)*xps(idx3));
   -- H(idx(1),idx(3)) := H(idx(1),idx(3)) + acc*powfac*x(idx(4));
    pr := xr(idx4); pi := xi(idx4);
    zr := racc*pr - iacc*pi; zi := racc*pi + iacc*pr;
    hrprow(idx3) := hrprow(idx3) + powfac*zr;
    hiprow(idx3) := hiprow(idx3) + powfac*zi;
    powfac := double_float(xps(idx1)*xps(idx4));
   -- H(idx(1),idx(4)) := H(idx(1),idx(4)) + acc*powfac*x(idx(3));
    pr := xr(idx3); pi := xi(idx3);
    zr := racc*pr - iacc*pi; zi := racc*pi + iacc*pr;
    hrprow(idx4) := hrprow(idx4) + powfac*zr;
    hiprow(idx4) := hiprow(idx4) + powfac*zi;
   -- acc := offdiagfac*x(idx(1));
    pr := xr(idx1); pi := xi(idx1);
    racc := roffdiagfac*pr - ioffdiagfac*pi;
    iacc := roffdiagfac*pi + ioffdiagfac*pr;
    powfac := double_float(xps(idx2)*xps(idx3));
   -- H(idx(2),idx(3)) := H(idx(2),idx(3)) + acc*powfac*x(idx(4));
    pr := xr(idx4); pi := xi(idx4);
    zr := racc*pr - iacc*pi; zi := racc*pi + iacc*pr;
    hrprow := hrp(idx2); hiprow := hip(idx2);
    hrprow(idx3) := hrprow(idx3) + powfac*zr;
    hiprow(idx3) := hiprow(idx3) + powfac*zi;
    powfac := double_float(xps(idx2)*xps(idx4));
   -- H(idx(2),idx(4)) := H(idx(2),idx(4)) + acc*powfac*x(idx(3));
    pr := xr(idx3); pi := xi(idx3);
    zr := racc*pr - iacc*pi; zi := racc*pi + iacc*pr;
    hrprow(idx4) := hrprow(idx4) + powfac*zr;
    hiprow(idx4) := hiprow(idx4) + powfac*zi;
   -- compute the diagonal elements
    for k in fac'range loop
      fidx1 := fac(k); m1 := xps(fidx1); powfac := double_float(m1*(m1-1));
     -- acc := powfac*ondiagfac; -- acc is the cofactor
      racc := powfac*rondiagfac; iacc := powfac*iondiagfac;
      for i in idx'range loop
        if idx(i) /= fidx1 then -- skip the current factor
          if xps(idx(i)) = 1  then
           -- acc := acc*x(idx(i));
            pr := xr(idx(i)); pi := xi(idx(i));
            zr := racc*pr - iacc*pi; zi := racc*pi + iacc*pr;
          else
           -- acc := acc*(pwt(idx(i))(1));
            pr := rpwt(idx(i))(1); pi := ipwt(idx(i))(1);
            zr := racc*pr - iacc*pi; zi := racc*pi + iacc*pr;
          end if;
          racc := zr; iacc := zi;
        end if;
      end loop;
     -- H(fac(k),fac(k)) := H(fac(k),fac(k)) + acc;
      hrprow := hrp(fidx1); hiprow := hip(fidx1);
      hrprow(fidx1) := hrprow(fidx1) + racc;
      hiprow(fidx1) := hiprow(fidx1) + iacc;
    end loop;
  end Speel4;

  procedure SpeelN ( H : in out Standard_Complex_Matrices.Matrix;
                     c : in Complex_Number;
                     xps : in Standard_Integer_Vectors.Vector;
                     idx : in Standard_Integer_Vectors.Vector;
                     fac : in Standard_Integer_Vectors.Vector;
                     x : in Standard_Complex_Vectors.Vector;
                     fwd : in Standard_Complex_Vectors.Link_to_Vector;
                     bck : in Standard_Complex_Vectors.Link_to_Vector;
                     pwt : in Standard_Complex_VecVecs.VecVec ) is

    sz : constant integer32 := idx'last;
    m1 : integer32;
    powfac : double_float; -- multiplier factor of two powers
    acc : Complex_Number;
    offdiagfac : Complex_Number; -- common off diagonal factor
    ondiagfac : Complex_Number;  -- common on diagonal factor

  begin
    offdiagfac := c; ondiagfac := c;  -- off/on diagonal factors
    for i in fac'range loop
      m1 := xps(fac(i));
      if m1 = 2 then
        offdiagfac := offdiagfac*x(fac(i));
      elsif m1 = 3 then
        offdiagfac := offdiagfac*(pwt(fac(i))(1));
        ondiagfac := ondiagfac*x(fac(i));
      else
        offdiagfac := offdiagfac*(pwt(fac(i))(m1-2));
        ondiagfac := ondiagfac*(pwt(fac(i))(m1-3));
      end if;
    end loop;
   -- the off diagonal elements use forward and backward products
   -- last element is copy of fwd(sz-3), multiplied with c
    powfac := double_float(xps(idx(sz-1))*xps(idx(sz)));
    H(idx(sz-1),idx(sz)) := H(idx(sz-1),idx(sz)) + offdiagfac*powfac*fwd(sz-3);
   -- first element is copy of bck(sz-3), multiplied with c
    powfac := double_float(xps(idx(1))*xps(idx(2)));
    H(idx(1),idx(2)) := H(idx(1),idx(2)) + offdiagfac*powfac*bck(sz-3);
   -- first row is special, starts with x(idx(2)) after diagonal
    acc := offdiagfac*x(idx(2));
    powfac := double_float(xps(idx(1))*xps(idx(3)));
    H(idx(1),idx(3)) := H(idx(1),idx(3)) + acc*powfac*bck(sz-4);
    for k in 4..sz-2 loop
      acc := acc*x(idx(k-1));
      powfac := double_float(xps(idx(1))*xps(idx(k)));
      H(idx(1),idx(k)) := H(idx(1),idx(k)) + acc*powfac*bck(sz-k-1);
    end loop;
    acc := acc*x(idx(sz-2));
    powfac := double_float(xps(idx(1))*xps(idx(sz-1)));
    H(idx(1),idx(sz-1)) := H(idx(1),idx(sz-1)) + acc*powfac*x(idx(sz));
    powfac := double_float(xps(idx(1))*xps(idx(sz)));
    H(idx(1),idx(sz)) := H(idx(1),idx(sz)) + acc*powfac*x(idx(sz-1));
   -- second row is special, starts with x(idx(1)) after diagonal
    acc := offdiagfac*x(idx(1));
    powfac := double_float(xps(idx(2))*xps(idx(3)));
    H(idx(2),idx(3)) := H(idx(2),idx(3)) + acc*powfac*bck(sz-4);
    for k in 4..sz-2 loop
      acc := acc*x(idx(k-1));
      powfac := double_float(xps(idx(2))*xps(idx(k)));
      H(idx(2),idx(k)) := H(idx(2),idx(k)) + acc*powfac*bck(sz-k-1);
    end loop;
    acc := acc*x(idx(sz-2));
    powfac := double_float(xps(idx(2))*xps(idx(sz-1)));
    H(idx(2),idx(sz-1)) := H(idx(2),idx(sz-1)) + acc*powfac*x(idx(sz));
    powfac := double_float(xps(idx(2))*xps(idx(sz)));
    H(idx(2),idx(sz)) := H(idx(2),idx(sz)) + acc*powfac*x(idx(sz-1));
   -- the row with index sz-2 has a general formula
    acc := offdiagfac*fwd(sz-4);
    powfac := double_float(xps(idx(sz-2))*xps(idx(sz-1)));
    H(idx(sz-2),idx(sz-1)) := H(idx(sz-2),idx(sz-1)) + acc*powfac*x(idx(sz));
    powfac := double_float(xps(idx(sz-2))*xps(idx(sz)));
    H(idx(sz-2),idx(sz)) := H(idx(sz-2),idx(sz)) + acc*powfac*x(idx(sz-1));
    for rw in 3..sz-3 loop  -- row rw starts with fwd(rw-2)
      acc := offdiagfac*fwd(rw-2);
      powfac := double_float(xps(idx(rw))*xps(idx(rw+1)));
      H(idx(rw),idx(rw+1)) := H(idx(rw),idx(rw+1)) + acc*powfac*bck(sz-rw-2);
      for k in rw+2..sz-2 loop
        acc := acc*x(idx(k-1));
        powfac := double_float(xps(idx(rw))*xps(idx(k)));
        H(idx(rw),idx(k)) := H(idx(rw),idx(k)) + acc*powfac*bck(sz-k-1);
      end loop;
      acc := acc*x(idx(sz-2));
      powfac := double_float(xps(idx(rw))*xps(idx(sz-1)));
      H(idx(rw),idx(sz-1)) := H(idx(rw),idx(sz-1)) + acc*powfac*x(idx(sz));
      powfac := double_float(xps(idx(rw))*xps(idx(sz)));
      H(idx(rw),idx(sz)) := H(idx(rw),idx(sz)) + acc*powfac*x(idx(sz-1));
    end loop;
   -- compute the diagonal elements
    for k in fac'range loop
      m1 := xps(fac(k)); powfac := double_float(m1*(m1-1));
      acc := powfac*ondiagfac; -- acc is the cofactor
      for i in idx'range loop
        if idx(i) /= fac(k) then -- skip the current factor
          if xps(idx(i)) = 1 
           then acc := acc*x(idx(i));
           else acc := acc*(pwt(idx(i))(1));
          end if;
        end if;
      end loop;
      H(fac(k),fac(k)) := H(fac(k),fac(k)) + acc;
    end loop;
   -- the above loop for the diagonal elements applies a loop
   -- for the cofactor, a similar triple loop with forward, backward,
   -- and cross porducts is possible for all fac'last cofactors
  end SpeelN;

  procedure SpeelN ( hrp : in Standard_Floating_VecVecs.VecVec;
                     hip : in Standard_Floating_VecVecs.VecVec;
                     rcff,icff : in double_float;
                     xps : in Standard_Integer_Vectors.Vector;
                     idx : in Standard_Integer_Vectors.Vector;
                     fac : in Standard_Integer_Vectors.Vector;
                     xr : in Standard_Floating_Vectors.Vector;
                     xi : in Standard_Floating_Vectors.Vector;
                     rfwd : in Standard_Floating_Vectors.Link_to_Vector;
                     ifwd : in Standard_Floating_Vectors.Link_to_Vector;
                     rbck : in Standard_Floating_Vectors.Link_to_Vector;
                     ibck : in Standard_Floating_Vectors.Link_to_Vector;
                     rpwt : in Standard_Floating_VecVecs.VecVec;
                     ipwt : in Standard_Floating_VecVecs.VecVec ) is

    sz : constant integer32 := idx'last;
    m1,idx1,fidx1 : integer32;
    powfac : double_float; -- multiplier factor of two powers
    pr,pi,zr,zi,racc,iacc : double_float;
    roffdiagfac,ioffdiagfac : double_float; -- common off diagonal factor
    rondiagfac,iondiagfac : double_float;   -- common on diagonal factor
    hrprow,hiprow : Standard_Floating_Vectors.Link_to_Vector;

  begin
   -- offdiagfac := c; ondiagfac := c;  -- off/on diagonal factors
    roffdiagfac := rcff; ioffdiagfac := icff;
    rondiagfac := rcff;  iondiagfac := icff;
    for i in fac'range loop
      fidx1 := fac(i); m1 := xps(fidx1);
      if m1 = 2 then
       -- offdiagfac := offdiagfac*x(fac(i));
        pr := xr(fidx1); pi := xi(fidx1);
        zr := roffdiagfac*pr - ioffdiagfac*pi;
        zi := roffdiagfac*pi + ioffdiagfac*pr;
        roffdiagfac := zr; ioffdiagfac := zi;
      elsif m1 = 3 then
       -- offdiagfac := offdiagfac*(pwt(fac(i))(1));
        pr := rpwt(fidx1)(1); pi := ipwt(fidx1)(1);
        zr := roffdiagfac*pr - ioffdiagfac*pi;
        zi := roffdiagfac*pi + ioffdiagfac*pr;
        roffdiagfac := zr; ioffdiagfac := zi;
       -- ondiagfac := ondiagfac*x(fac(i));
        pr := xr(fidx1); pi := xi(fidx1);
        zr := rondiagfac*pr - iondiagfac*pi;
        zi := rondiagfac*pi + iondiagfac*pr;
        rondiagfac := zr; iondiagfac := zi;
      else
       -- offdiagfac := offdiagfac*(pwt(fac(i))(m1-2));
        pr := rpwt(fidx1)(m1-2); pi := ipwt(fidx1)(m1-2);
        zr := roffdiagfac*pr - ioffdiagfac*pi;
        zi := roffdiagfac*pi + ioffdiagfac*pr;
        roffdiagfac := zr; ioffdiagfac := zi;
       -- ondiagfac := ondiagfac*(pwt(fac(i))(m1-3));
        pr := rpwt(fidx1)(m1-3); pi := ipwt(fidx1)(m1-3);
        zr := rondiagfac*pr - iondiagfac*pi;
        zi := rondiagfac*pi + iondiagfac*pr;
        rondiagfac := zr; iondiagfac := zi;
      end if;
    end loop;
   -- the off diagonal elements use forward and backward products
   -- last element is copy of fwd(sz-3), multiplied with c
    powfac := double_float(xps(idx(sz-1))*xps(idx(sz)));
   -- H(idx(sz-1),idx(sz))
   --   := H(idx(sz-1),idx(sz)) + offdiagfac*powfac*fwd(sz-3);
    pr := rfwd(sz-3); pi := ifwd(sz-3);
    zr := roffdiagfac*pr - ioffdiagfac*pi;
    zi := roffdiagfac*pi + ioffdiagfac*pr;
    idx1 := idx(sz-1); hrprow := hrp(idx1); hiprow := hip(idx1);
    idx1 := idx(sz);
    hrprow(idx1) := hrprow(idx1) + powfac*zr;
    hiprow(idx1) := hiprow(idx1) + powfac*zi;
   -- first element is copy of bck(sz-3), multiplied with c
    powfac := double_float(xps(idx(1))*xps(idx(2)));
   -- H(idx(1),idx(2)) := H(idx(1),idx(2)) + offdiagfac*powfac*bck(sz-3);
    pr := rbck(sz-3); pi := ibck(sz-3);
    zr := roffdiagfac*pr - ioffdiagfac*pi;
    zi := roffdiagfac*pi + ioffdiagfac*pr;
    idx1 := idx(1); hrprow := hrp(idx1); hiprow := hip(idx1);
    idx1 := idx(2);
    hrprow(idx1) := hrprow(idx1) + powfac*zr;
    hiprow(idx1) := hiprow(idx1) + powfac*zi;
   -- first row is special, starts with x(idx(2)) after diagonal
   -- acc := offdiagfac*x(idx(2));
    idx1 := idx(2); pr := xr(idx1); pi := xi(idx1);
    racc := roffdiagfac*pr - ioffdiagfac*pi;
    iacc := roffdiagfac*pi + ioffdiagfac*pr;
    powfac := double_float(xps(idx(1))*xps(idx(3)));
   -- H(idx(1),idx(3)) := H(idx(1),idx(3)) + acc*powfac*bck(sz-4);
    pr := rbck(sz-4); pi := ibck(sz-4);
    zr := racc*pr - iacc*pi; zi := racc*pi + iacc*pr;
   -- idx1 := idx(1); hrprow := hrp(idx1); hiprow := hip(idx1);
    idx1 := idx(3);
    hrprow(idx1) := hrprow(idx1) + powfac*zr;
    hiprow(idx1) := hiprow(idx1) + powfac*zi;
    for k in 4..sz-2 loop
     -- acc := acc*x(idx(k-1));
      idx1 := idx(k-1); pr := xr(idx1); pi := xi(idx1);
      zr := racc*pr - iacc*pi; zi := racc*pi + iacc*pr;
      racc := zr; iacc := zi;
      powfac := double_float(xps(idx(1))*xps(idx(k)));
     -- H(idx(1),idx(k)) := H(idx(1),idx(k)) + acc*powfac*bck(sz-k-1);
      pr := rbck(sz-k-1); pi := ibck(sz-k-1);
      zr := racc*pr - iacc*pi; zi := racc*pi + iacc*pr;
      idx1 := idx(k);
      hrprow(idx1) := hrprow(idx1) + powfac*zr;
      hiprow(idx1) := hiprow(idx1) + powfac*zi;
    end loop;
   -- acc := acc*x(idx(sz-2));
    idx1 := idx(sz-2); pr := xr(idx1); pi := xi(idx1);
    zr := racc*pr - iacc*pi; zi := racc*pi + iacc*pr;
    racc := zr; iacc := zi;
    powfac := double_float(xps(idx(1))*xps(idx(sz-1)));
   -- H(idx(1),idx(sz-1)) := H(idx(1),idx(sz-1)) + acc*powfac*x(idx(sz));
    idx1 := idx(sz); pr := xr(idx1); pi := xi(idx1);
    zr := racc*pr - iacc*pi; zi := racc*pi + iacc*pr;
    idx1 := idx(sz-1);
    hrprow(idx1) := hrprow(idx1) + powfac*zr;
    hiprow(idx1) := hiprow(idx1) + powfac*zi;
    powfac := double_float(xps(idx(1))*xps(idx(sz)));
   -- H(idx(1),idx(sz)) := H(idx(1),idx(sz)) + acc*powfac*x(idx(sz-1));
    idx1 := idx(sz-1); pr := xr(idx1); pi := xi(idx1);
    zr := racc*pr - iacc*pi; zi := racc*pi + iacc*pr;
    idx1 := idx(sz);
    hrprow(idx1) := hrprow(idx1) + powfac*zr;
    hiprow(idx1) := hiprow(idx1) + powfac*zi;
   -- second row is special, starts with x(idx(1)) after diagonal
   -- acc := offdiagfac*x(idx(1));
    idx1 := idx(1); pr := xr(idx1); pi := xi(idx1);
    racc := roffdiagfac*pr - ioffdiagfac*pi;
    iacc := roffdiagfac*pi + ioffdiagfac*pr;
    powfac := double_float(xps(idx(2))*xps(idx(3)));
   -- H(idx(2),idx(3)) := H(idx(2),idx(3)) + acc*powfac*bck(sz-4);
    pr := rbck(sz-4); pi := ibck(sz-4);
    zr := racc*pr - iacc*pi; zi := racc*pi + iacc*pr;
    idx1 := idx(2); hrprow := hrp(idx1); hiprow := hip(idx1);
    idx1 := idx(3);
    hrprow(idx1) := hrprow(idx1) + powfac*zr;
    hiprow(idx1) := hiprow(idx1) + powfac*zi;
    for k in 4..sz-2 loop
     -- acc := acc*x(idx(k-1));
      idx1 := idx(k-1); pr := xr(idx1); pi := xi(idx1);
      zr := racc*pr - iacc*pi; zi := racc*pi + iacc*pr;
      racc := zr; iacc := zi;
      powfac := double_float(xps(idx(2))*xps(idx(k)));
     -- H(idx(2),idx(k)) := H(idx(2),idx(k)) + acc*powfac*bck(sz-k-1);
      pr := rbck(sz-k-1); pi := ibck(sz-k-1);
      zr := racc*pr - iacc*pi; zi := racc*pi + iacc*pr;
      idx1 := idx(k);
      hrprow(idx1) := hrprow(idx1) + powfac*zr;
      hiprow(idx1) := hiprow(idx1) + powfac*zi;
    end loop;
   -- acc := acc*x(idx(sz-2));
    idx1 := idx(sz-2); pr := xr(idx1); pi := xi(idx1);
    zr := racc*pr - iacc*pi; zi := racc*pi + iacc*pr;
    racc := zr; iacc := zi;
    powfac := double_float(xps(idx(2))*xps(idx(sz-1)));
   -- H(idx(2),idx(sz-1)) := H(idx(2),idx(sz-1)) + acc*powfac*x(idx(sz));
    idx1 := idx(sz); pr := xr(idx1); pi := xi(idx1);
    zr := racc*pr - iacc*pi; zi := racc*pi + iacc*pr;
    idx1 := idx(sz-1);
    hrprow(idx1) := hrprow(idx1) + powfac*zr;
    hiprow(idx1) := hiprow(idx1) + powfac*zi;
    powfac := double_float(xps(idx(2))*xps(idx(sz)));
   -- H(idx(2),idx(sz)) := H(idx(2),idx(sz)) + acc*powfac*x(idx(sz-1));
    idx1 := idx(sz-1); pr := xr(idx1); pi := xi(idx1);
    zr := racc*pr - iacc*pi; zi := racc*pi + iacc*pr;
    idx1 := idx(sz);
    hrprow(idx1) := hrprow(idx1) + powfac*zr;
    hiprow(idx1) := hiprow(idx1) + powfac*zi;
   -- the row with index sz-2 has a general formula
   -- acc := offdiagfac*fwd(sz-4);
    pr := rfwd(sz-4); pi := ifwd(sz-4);
    racc := roffdiagfac*pr - ioffdiagfac*pi;
    iacc := roffdiagfac*pi + ioffdiagfac*pr;
    powfac := double_float(xps(idx(sz-2))*xps(idx(sz-1)));
   -- H(idx(sz-2),idx(sz-1)) := H(idx(sz-2),idx(sz-1)) + acc*powfac*x(idx(sz));
    idx1 := idx(sz); pr := xr(idx1); pi := xi(idx1);
    zr := racc*pr - iacc*pi; zi := racc*pi + iacc*pr;
    idx1 := idx(sz-2); hrprow := hrp(idx1); hiprow := hip(idx1);
    idx1 := idx(sz-1);
    hrprow(idx1) := hrprow(idx1) + powfac*zr;
    hiprow(idx1) := hiprow(idx1) + powfac*zi;
    powfac := double_float(xps(idx(sz-2))*xps(idx(sz)));
   -- H(idx(sz-2),idx(sz)) := H(idx(sz-2),idx(sz)) + acc*powfac*x(idx(sz-1));
    idx1 := idx(sz-1); pr := xr(idx1); pi := xi(idx1);
    zr := racc*pr - iacc*pi; zi := racc*pi + iacc*pr;
    idx1 := idx(sz);
    hrprow(idx1) := hrprow(idx1) + powfac*zr;
    hiprow(idx1) := hiprow(idx1) + powfac*zi;
    for rw in 3..sz-3 loop  -- row rw starts with fwd(rw-2)
     -- acc := offdiagfac*fwd(rw-2);
      pr := rfwd(rw-2); pi := ifwd(rw-2);
      racc := roffdiagfac*pr - ioffdiagfac*pi;
      iacc := roffdiagfac*pi + ioffdiagfac*pr;
      powfac := double_float(xps(idx(rw))*xps(idx(rw+1)));
     -- H(idx(rw),idx(rw+1)) := H(idx(rw),idx(rw+1)) + acc*powfac*bck(sz-rw-2);
      pr := rbck(sz-rw-2); pi := ibck(sz-rw-2);
      zr := racc*pr - iacc*pi; zi := racc*pi + iacc*pr;
      idx1 := idx(rw); hrprow := hrp(idx1); hiprow := hip(idx1);
      idx1 := idx(rw+1);
      hrprow(idx1) := hrprow(idx1) + powfac*zr;
      hiprow(idx1) := hiprow(idx1) + powfac*zi;
      for k in rw+2..sz-2 loop
       -- acc := acc*x(idx(k-1));
        idx1 := idx(k-1); pr := xr(idx1); pi := xi(idx1);
        zr := racc*pr - iacc*pi; zi := racc*pi + iacc*pr;
        racc := zr; iacc := zi;
        powfac := double_float(xps(idx(rw))*xps(idx(k)));
       -- H(idx(rw),idx(k)) := H(idx(rw),idx(k)) + acc*powfac*bck(sz-k-1);
        pr := rbck(sz-k-1); pi := ibck(sz-k-1);
        zr := racc*pr - iacc*pi; zi := racc*pi + iacc*pr;
        idx1 := idx(k);
        hrprow(idx1) := hrprow(idx1) + powfac*zr;
        hiprow(idx1) := hiprow(idx1) + powfac*zi;
      end loop;
     -- acc := acc*x(idx(sz-2));
      idx1 := idx(sz-2); pr := xr(idx1); pi := xi(idx1);
      zr := racc*pr - iacc*pi; zi := racc*pi + iacc*pr;
      racc := zr; iacc := zi;
      powfac := double_float(xps(idx(rw))*xps(idx(sz-1)));
     -- H(idx(rw),idx(sz-1)) := H(idx(rw),idx(sz-1)) + acc*powfac*x(idx(sz));
      idx1 := idx(sz); pr := xr(idx1); pi := xi(idx1);
      zr := racc*pr - iacc*pi; zi := racc*pi + iacc*pr;
      idx1 := idx(sz-1);
      hrprow(idx1) := hrprow(idx1) + powfac*zr;
      hiprow(idx1) := hiprow(idx1) + powfac*zi;
      powfac := double_float(xps(idx(rw))*xps(idx(sz)));
     -- H(idx(rw),idx(sz)) := H(idx(rw),idx(sz)) + acc*powfac*x(idx(sz-1));
      idx1 := idx(sz-1); pr := xr(idx1); pi := xi(idx1);
      zr := racc*pr - iacc*pi; zi := racc*pi + iacc*pr;
      idx1 := idx(sz);
      hrprow(idx1) := hrprow(idx1) + powfac*zr;
      hiprow(idx1) := hiprow(idx1) + powfac*zi;
    end loop;
   -- compute the diagonal elements
    for k in fac'range loop
      fidx1 := fac(k); m1 := xps(fidx1); powfac := double_float(m1*(m1-1));
     -- acc := powfac*ondiagfac; -- acc is the cofactor
      racc := powfac*rondiagfac; iacc := powfac*iondiagfac;
      for i in idx'range loop
        if idx(i) /= fidx1 then -- skip the current factor
          if xps(idx(i)) = 1  then
           -- acc := acc*x(idx(i));
            pr := xr(idx(i)); pi := xi(idx(i));
            zr := racc*pr - iacc*pi; zi := racc*pi + iacc*pr;
          else
           -- acc := acc*(pwt(idx(i))(1));
            pr := rpwt(idx(i))(1); pi := ipwt(idx(i))(1);
            zr := racc*pr - iacc*pi; zi := racc*pi + iacc*pr;
          end if;
          racc := zr; iacc := zi;
        end if;
      end loop;
     -- H(fac(k),fac(k)) := H(fac(k),fac(k)) + acc;
      hrprow := hrp(fidx1); hiprow := hip(fidx1);
      hrprow(fidx1) := hrprow(fidx1) + racc;
      hiprow(fidx1) := hiprow(fidx1) + iacc;
    end loop;
  end SpeelN;

end Standard_Hessian_Updaters;
