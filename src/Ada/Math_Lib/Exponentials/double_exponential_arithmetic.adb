with Standard_Complex_Numbers;           use Standard_Complex_Numbers;
-- with Standard_Random_Numbers;

package body Double_Exponential_Arithmetic is

  function Extension_Degree ( alpha,beta : double_float ) return integer32 is

    res : integer32 := 0;

  begin
    if beta = 0.0 then
      return 0;
    else
      res := integer32(alpha/beta);
      while double_float(res)*beta < alpha loop
        res := res + 1;
      end loop;
      return res;
    end if;
  end Extension_Degree;

  function Extension_Degree
	     ( alpha : double_float;
               beta : Standard_Floating_Vectors.Vector) return integer32 is

    res,deg : integer32 := 0;

  begin
    for i in beta'range loop
      deg := Extension_Degree(alpha,beta(i));
      if deg > res
       then res := deg;
      end if;
    end loop;
    return res;
  end Extension_Degree;

  procedure Normalize
              ( cff : in out Standard_Complex_Vectors.Vector;
                sxp : in out Standard_Floating_Vectors.Vector ) is

    cfftmp : Complex_Number;
    sxptmp : double_float;
    swapped : boolean;

  begin
    loop
      swapped := false;
      for i in sxp'first+1..sxp'last loop
        if sxp(i) < sxp(i-1) then
          sxptmp := sxp(i); sxp(i) := sxp(i-1); sxp(i-1) := sxptmp;
          cfftmp := cff(i); cff(i) := cff(i-1); cff(i-1) := cfftmp;
          swapped := true;
        end if;
      end loop;
      exit when not swapped;
    end loop;
  end Normalize;

  function Quadratic_Extend_Size ( deg : integer32 ) return integer32 is

   -- res : constant integer32 := 2*deg + deg*(deg-1)/2;
    res : constant integer32 := deg*(deg+3)/2;

  begin
    return res;
  end Quadratic_Extend_Size;

  procedure Quadratic_Extend
              ( deg,size : in integer32;
                cff : in Standard_Complex_Vectors.Vector;
                sxp : in Standard_Floating_Vectors.Vector;
                extcff : out Standard_Complex_Vectors.Vector;
                extsxp : out Standard_Floating_Vectors.Vector ) is

    idx : integer32;

  begin
    extcff(cff'range) := cff;
    extsxp(sxp'range) := sxp;
    extcff(deg+1..size) := (deg+1..size => create(1.0));
    idx := deg + 1;
    for i in 1..deg loop
      for j in i+1..deg loop
        extsxp(idx) := sxp(i) + sxp(j);
        idx := idx + 1;
      end loop;
    end loop;
    for i in 1..deg loop
      extsxp(idx) := 2.0*sxp(i);
      idx := idx + 1;
    end loop;
    Normalize(extcff,extsxp);
  end Quadratic_Extend;

  function Inverse ( cff : Standard_Complex_Vectors.Vector )
                    return Standard_Complex_Vectors.Vector is

    res : Standard_Complex_Vectors.Vector(cff'range);

  begin
    res(0) := 1.0/cff(0);
    for i in 1..res'last loop
      res(i) := -cff(1)*res(i-1);
      for j in 2..i loop
        res(i) := res(i) - cff(j)*res(i-j);
      end loop;
      res(i) := res(i)/cff(0);
    end loop;
    return res;
  end Inverse;

  function Convolute ( a,b : Standard_Complex_Vectors.Vector )
                      return Standard_Complex_Vectors.Vector is

    res : Standard_Complex_Vectors.Vector(0..a'last);

  begin
    for i in 0..res'last loop
      res(i) := a(0)*b(i);
      for j in 1..i loop
        res(i) := res(i) + a(j)*b(i-j);
      end loop;
    end loop;
    return res;
  end Convolute;

  procedure Zero_Padding
              ( deg : in integer32; idx : in out integer32;
                cff : in out Standard_Complex_Vectors.Vector;
                sxp : in out Standard_Floating_Vectors.Vector ) is

  -- DESCRIPTION :
  --   Adds extra zeroes at the end as long as idx <= deg.
  --   The corresponding exponents are raised by one each time.

  begin
    while idx <= deg loop
      cff(idx) := create(0.0);
      if idx = 0
       then sxp(idx) := 1.0;
       else sxp(idx) := sxp(idx-1) + 1.0;
      end if;
      idx := idx + 1;
    end loop;
  end Zero_Padding;

  procedure Add ( adeg,bdeg,cdeg : in integer32;
                  acf,bcf : in Standard_Complex_Vectors.Vector;
                  axp,bxp : in Standard_Floating_Vectors.Vector;
                  ccf : out Standard_Complex_Vectors.Vector;
                  cxp : out Standard_Floating_Vectors.Vector;
                  tol : in double_float := 1.0e-14 ) is

    aix : integer32 := acf'first;
    bix : integer32 := bcf'first;
    cix : integer32 := ccf'first;

  begin
    while cix <= cdeg loop
      if axp(aix) < bxp(bix) then
        ccf(cix) := acf(aix);
        cxp(cix) := axp(aix);
        aix := aix + 1;
        cix := cix + 1;
      elsif axp(aix) > bxp(bix) then
        ccf(cix) := bcf(bix);
        cxp(cix) := bxp(bix);
        bix := bix + 1;
        cix := cix + 1;
      else -- axp(aix) = bxp(bix) 
        ccf(cix) := acf(aix) + bcf(bix);
        cxp(cix) := axp(aix);
        aix := aix + 1;
        bix := bix + 1;
        if AbsVal(ccf(cix)) > tol
         then cix := cix + 1;
        end if;
      end if;
      exit when (aix > adeg) or (bix > bdeg);
    end loop;
    if cix <= cdeg then
      if aix <= adeg then
        while cix <= cdeg loop
          ccf(cix) := acf(aix);
          cxp(cix) := axp(aix);
          cix := cix + 1;
          aix := aix + 1;
          exit when (aix > adeg);
        end loop;
      elsif bix <= bdeg then
        while cix <= cdeg loop
          ccf(cix) := bcf(bix);
          cxp(cix) := bxp(bix);
          cix := cix + 1;
          bix := bix + 1;
          exit when (bix > bdeg);
        end loop;
      end if;
      Zero_Padding(cdeg,cix,ccf,cxp);
    end if;
  end Add;

  procedure Sub ( adeg,bdeg,cdeg : in integer32;
                  acf,bcf : in Standard_Complex_Vectors.Vector;
                  axp,bxp : in Standard_Floating_Vectors.Vector;
                  ccf : out Standard_Complex_Vectors.Vector;
                  cxp : out Standard_Floating_Vectors.Vector;
                  tol : in double_float := 1.0e-14 ) is

    aix : integer32 := acf'first;
    bix : integer32 := bcf'first;
    cix : integer32 := ccf'first;

  begin
    while cix <= cdeg loop
      if axp(aix) < bxp(bix) then
        ccf(cix) := acf(aix);
        cxp(cix) := axp(aix);
        aix := aix + 1;
        cix := cix + 1;
      elsif axp(aix) > bxp(bix) then
        ccf(cix) := -bcf(bix);
        cxp(cix) := bxp(bix);
        bix := bix + 1;
        cix := cix + 1;
      else -- axp(aix) = bxp(bix) 
        ccf(cix) := acf(aix) - bcf(bix);
        cxp(cix) := axp(aix);
        aix := aix + 1;
        bix := bix + 1;
        if AbsVal(ccf(cix)) > tol
         then cix := cix + 1;
        end if;
      end if;
      exit when (aix > adeg) or (bix > bdeg);
    end loop;
    if cix <= cdeg then
      if aix <= adeg then
        while cix <= cdeg loop
          ccf(cix) := acf(aix);
          cxp(cix) := axp(aix);
          cix := cix + 1;
          aix := aix + 1;
          exit when (aix > adeg);
        end loop;
      elsif bix <= bdeg then
        while cix <= cdeg loop
          ccf(cix) := -bcf(bix);
          cxp(cix) := bxp(bix);
          cix := cix + 1;
          bix := bix + 1;
          exit when (bix > bdeg);
        end loop;
      end if;
      Zero_Padding(cdeg,cix,ccf,cxp);
    end if;
  end Sub;

  procedure Mul ( adeg,bdeg,cdeg : in integer32;
                  acf,bcf : in Standard_Complex_Vectors.Vector;
                  axp,bxp : in Standard_Floating_Vectors.Vector;
                  ccf : out Standard_Complex_Vectors.Vector;
                  cxp : out Standard_Floating_Vectors.Vector;
                  prdcf,wrkcf : in out Standard_Complex_Vectors.Vector;
                  prdxp,wrkxp : in out Standard_Floating_Vectors.Vector;
                  tol : in double_float := 1.0e-14 ) is

    cix,pix,wix : integer32;

  begin
    ccf := (0..cdeg => Create(0.0));
    cxp := (0..cdeg => 0.0);
    prdcf := (prdcf'range => Create(0.0));
    prdxp := (prdxp'range => 0.0);
    wrkcf := (wrkcf'range => Create(0.0));
    wrkxp := (wrkxp'range => 0.0);
    for i in bcf'range loop   -- product of a(0) with b series
      ccf(i) := acf(0) * bcf(i);
      cxp(i) := axp(0) + bxp(i);
    end loop;
    for i in 1..adeg loop
      prdcf := (prdcf'range => create(0.0));
      prdxp := (prdxp'range => 0.0);
      for j in bcf'range loop -- product of a(i) with b series
        prdcf(j) := acf(i) * bcf(j);
        prdxp(j) := axp(i) + bxp(j);
      end loop;
      cix := ccf'first;
      pix := prdcf'first;
      wix := wrkcf'first;
      while wix <= bdeg + i*(bdeg+1) loop
        if cxp(cix) < prdxp(pix) then
          wrkcf(wix) := ccf(cix);
          wrkxp(wix) := cxp(cix);
          cix := cix + 1;
          wix := wix + 1;
        elsif cxp(cix) > prdxp(pix) then
          wrkcf(wix) := prdcf(pix);
          wrkxp(wix) := prdxp(pix);
          pix := pix + 1;
          wix := wix + 1;
        else -- cxp(cix) = prdxp(pix)
          wrkcf(wix) := ccf(cix) + prdcf(pix); -- accumulate products
          wrkxp(wix) := cxp(cix);              -- with same exponent
          cix := cix + 1;
          pix := pix + 1;
          if AbsVal(wrkcf(wix)) > tol
           then wix := wix + 1;
          end if;
        end if;
        exit when (pix > bdeg) or (cix > bdeg + (i-1)*(bdeg+1));
        exit when (cix > cdeg) or (wix > cdeg);
      end loop;
      if wix <= bdeg + i*(bdeg+1) then
        if pix <= bdeg then
          while wix <= bdeg + i*(bdeg+1) loop
            wrkcf(wix) := prdcf(pix);
            wrkxp(wix) := prdxp(pix);
            wix := wix + 1;
            pix := pix + 1;
            exit when (pix > bdeg);
          end loop;
        end if;
        if cix <= bdeg + i*(bdeg+1) then
          while wix <= bdeg + i*(bdeg+1) loop
            wrkcf(wix) := ccf(cix);
            wrkxp(wix) := cxp(cix);
            wix := wix + 1;
            cix := cix + 1;
            exit when (cix > bdeg + i*(bdeg+1)) or (cix > cdeg);
          end loop;
        end if;
      end if;
      cix := ccf'first;
      for j in 0..bdeg+i*(bdeg+1) loop
        exit when (cix > cdeg);
        ccf(cix) := wrkcf(j);
        cxp(cix) := wrkxp(j);
        if AbsVal(ccf(cix)) > tol
         then cix := cix + 1;
        end if;
      end loop;
    end loop;
    if cix <= cdeg
     then Zero_Padding(cdeg,cix,ccf,cxp);
    end if;
  end Mul;

  procedure Div ( adeg,bdeg,cdeg : in integer32;
                  acf,bcf : in Standard_Complex_Vectors.Vector;
                  axp,bxp : in Standard_Floating_Vectors.Vector;
                  ccf : out Standard_Complex_Vectors.Vector;
                  cxp : out Standard_Floating_Vectors.Vector;
                  invbcf,prdcf,wrkcf : in out Standard_Complex_Vectors.Vector;
                  prdxp,wrkxp : in out Standard_Floating_Vectors.Vector ) is
  begin
    invbcf := Inverse(bcf);
    Mul(adeg,bdeg,cdeg,acf,invbcf,axp,bxp,ccf,cxp,prdcf,wrkcf,prdxp,wrkxp);
  end Div;

end Double_Exponential_Arithmetic;
