with Standard_Complex_Numbers;           use Standard_Complex_Numbers;

package body Double_Exponential_Arithmetic is

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
      while cix <= cdeg loop
        ccf(cix) := create(0.0);
        cxp(cix) := cxp(cix-1) + 1.0;
        cix := cix + 1;
      end loop;
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
      while cix <= cdeg loop
        ccf(cix) := create(0.0);
        if cix = 0
         then cxp(cix) := 1.0;
         else cxp(cix) := cxp(cix-1) + 1.0;
        end if;
        cix := cix + 1;
      end loop;
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
    wrkcf := (wrkcf'range => Create(0.0));
    for i in bcf'range loop   -- product of a(0) with b series
      ccf(i) := acf(0)*bcf(i);
      cxp(i) := axp(0)+bxp(i);
    end loop;
    for i in 1..adeg loop
      for j in bcf'range loop -- product of a(i) with b series
        prdcf(j) := acf(i)*bcf(j);
        prdxp(j) := axp(i)+bxp(j);
      end loop;
      cix := ccf'first;
      pix := prdcf'first;
      wix := wrkcf'first;
      while wix <= (i+1)*bdeg loop
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
          wrkcf(wix) := ccf(cix) + prdcf(pix);
          wrkxp(wix) := cxp(cix);
          cix := cix + 1;
          pix := pix + 1;
          if AbsVal(wrkcf(wix)) > tol
           then wix := wix + 1;
          end if;
        end if;
        exit when (pix > bdeg) or (cix > i*bdeg);
        exit when (cix > cdeg) or (wix > cdeg);
      end loop;
      if wix <= (i+1)*bdeg then
        if pix <= bdeg then
          while wix <= (i+1)*bdeg loop
            wrkcf(wix) := prdcf(pix);
            wrkxp(wix) := prdxp(pix);
            wix := wix + 1;
            pix := pix + 1;
            exit when (pix > bdeg);
          end loop;
        elsif cix <= i*bdeg then
          while wix <= (i+1)*bdeg loop
            wrkcf(wix) := ccf(cix);
            wrkxp(wix) := cxp(cix);
            wix := wix + 1;
            cix := cix + 1;
            exit when (cix > i*bdeg) or (cix > cdeg);
          end loop;
        end if;
      end if;
      for j in 0..(i+1)*bdeg loop
        exit when (j > cdeg);
        ccf(j) := wrkcf(j);
        cxp(j) := wrkxp(j);
      end loop;
    end loop;
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
