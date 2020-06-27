with Double_Double_Numbers;              use Double_Double_Numbers;

package body DoblDobl_Vector_Splitters is

  procedure Split ( x : in Complex_Number;
                    rehi,imhi,relo,imlo : out double_float ) is

   nbr : double_double;

  begin
    nbr := REAL_PART(x); rehi := hi_part(nbr); relo := lo_part(nbr);
    nbr := IMAG_PART(x); imhi := hi_part(nbr); imlo := lo_part(nbr);
  end Split;

  procedure Split ( v : in DoblDobl_Complex_Vectors.Vector;
                    rehi : out Standard_Floating_Vectors.Vector;
                    imhi : out Standard_Floating_Vectors.Vector;
                    relo : out Standard_Floating_Vectors.Vector;
                    imlo : out Standard_Floating_Vectors.Vector ) is
  begin
    for k in v'range loop
      Split(v(k),rehi(k),imhi(k),relo(k),imlo(k));
    end loop;
  end Split;

  procedure Split_Complex
              ( x : in DoblDobl_Complex_Vectors.Link_to_Vector;
                rhpx,ihpx : out Standard_Floating_Vectors.Link_to_Vector;
                rlpx,ilpx : out Standard_Floating_Vectors.Link_to_Vector ) is

    rhx,ihx,rlx,ilx : Standard_Floating_Vectors.Vector(x'range);

  begin
    for k in x'range loop
      Split(x(k),rhx(k),ihx(k),rlx(k),ilx(k));
    end loop;
    rhpx := new Standard_Floating_Vectors.Vector'(rhx);
    ihpx := new Standard_Floating_Vectors.Vector'(ihx);
    rlpx := new Standard_Floating_Vectors.Vector'(rlx);
    ilpx := new Standard_Floating_Vectors.Vector'(ilx);
  end Split_Complex;

  procedure Split_Complex
              ( x : in DoblDobl_Complex_VecVecs.VecVec;
                rhpx,ihpx : out Standard_Floating_VecVecs.VecVec;
                rlpx,ilpx : out Standard_Floating_VecVecs.VecVec ) is

    use DoblDobl_Complex_Vectors;
 
  begin
    for k in x'range loop
      if x(k) /= null
       then Split_Complex(x(k),rhpx(k),ihpx(k),rlpx(k),ilpx(k));
      end if;
    end loop;
  end Split_Complex;

  procedure Split_Complex
              ( x : in DoblDobl_Complex_VecVecs.Link_to_VecVec;
                rhpx,ihpx : out Standard_Floating_VecVecs.Link_to_VecVec;
                rlpx,ilpx : out Standard_Floating_VecVecs.Link_to_VecVec ) is

    use DoblDobl_Complex_VecVecs,DoblDobl_Complex_Vectors;

  begin
    if x /= null then
      declare
        rhx,ihx,rlx,ilx : Standard_Floating_VecVecs.VecVec(x'range);
      begin
        for k in x'range loop
          if x(k) /= null
           then Split_Complex(x(k),rhx(k),ihx(k),rlx(k),ilx(k));
          end if;
        end loop;
        rhpx := new Standard_Floating_VecVecs.VecVec'(rhx);
        ihpx := new Standard_Floating_VecVecs.VecVec'(ihx);
        rlpx := new Standard_Floating_VecVecs.VecVec'(rlx);
        ilpx := new Standard_Floating_VecVecs.VecVec'(ilx);
      end;
    end if;
  end Split_Complex;

-- MEMORY ALLOCATORS :

  function Allocate_Complex_Coefficients
             ( deg : integer32 )
             return DoblDobl_Complex_Vectors.Link_to_Vector is

    zero : constant DoblDobl_Complex_Numbers.Complex_Number
         := DoblDobl_Complex_Numbers.Create(integer(0));
    cff : constant DoblDobl_Complex_Vectors.Vector(0..deg) := (0..deg => zero);
    res : constant DoblDobl_Complex_Vectors.Link_to_Vector
        := new DoblDobl_Complex_Vectors.Vector'(cff);

  begin
    return res;
  end Allocate_Complex_Coefficients;

  function Allocate_Complex_Coefficients
             ( dim,deg : integer32 )
             return DoblDobl_Complex_VecVecs.VecVec is

    res : DoblDobl_Complex_VecVecs.VecVec(1..dim);

  begin
    for k in res'range loop
      res(k) := Allocate_Complex_Coefficients(deg);
    end loop;
    return res;
  end Allocate_Complex_Coefficients;

  function Allocate_Complex_Coefficients
             ( dim,deg : integer32 )
             return DoblDobl_Complex_VecVecs.Link_to_VecVec is

    res : DoblDobl_Complex_VecVecs.Link_to_VecVec;
    cff : constant DoblDobl_Complex_VecVecs.VecVec(1..dim)
        := Allocate_Complex_Coefficients(dim,deg);

  begin
    res := new DoblDobl_Complex_VecVecs.VecVec'(cff);
    return res;
  end Allocate_Complex_Coefficients;

  function Allocate ( neq,dim : integer32; neqstart,dimstart : integer32 )
                    return DoblDobl_Complex_VecVecs.VecVec is

    res : DoblDobl_Complex_VecVecs.VecVec(neqstart..neq);

  begin
    for i in res'range loop
      declare
        v : constant DoblDobl_Complex_Vectors.Vector(dimstart..dim)
          := (dimstart..dim => DoblDobl_Complex_Numbers.Create(integer(0)));
      begin
        res(i) := new DoblDobl_Complex_Vectors.Vector'(v);
      end;
    end loop;
    return res;
  end Allocate;

  function Allocate ( neq,dim : integer32; neqstart,dimstart : integer32 )
                    return DoblDobl_Complex_VecVecs.Link_to_VecVec is

    res : DoblDobl_Complex_VecVecs.Link_to_VecVec;
    vv : DoblDobl_Complex_VecVecs.VecVec(neqstart..neq);

  begin
    for i in vv'range loop
      declare
        v : constant DoblDobl_Complex_Vectors.Vector(dimstart..dim)
          := (dimstart..dim => DoblDobl_Complex_Numbers.Create(integer(0)));
      begin
        vv(i) := new DoblDobl_Complex_Vectors.Vector'(v);
      end;
    end loop;
    res := new DoblDobl_Complex_VecVecs.VecVec'(vv);
    return res;
  end Allocate;

-- MERGE PROCEDURES :

  procedure Merge ( x : out Complex_Number;
                    rehi,imhi,relo,imlo : in double_float ) is

    realpart : constant double_double
             := Double_Double_Numbers.Create(rehi,relo);
    imagpart : constant double_double
             := Double_Double_Numbers.Create(imhi,imlo);

  begin
    x := DoblDobl_Complex_Numbers.Create(realpart,imagpart);
  end Merge;

  procedure Merge ( v : out DoblDobl_Complex_Vectors.Vector;
                    rehi : in Standard_Floating_Vectors.Vector;
                    imhi : in Standard_Floating_Vectors.Vector;
                    relo : in Standard_Floating_Vectors.Vector;
                    imlo : in Standard_Floating_Vectors.Vector ) is
  begin
    for k in v'range loop
      Merge(v(k),rehi(k),imhi(k),relo(k),imlo(k));
    end loop;
  end Merge;

-- PROCEDURES TO PART AND MERGE VECTORS OF VECTORS :

  procedure Complex_Parts
              ( x : in DoblDobl_Complex_Vectors.Link_to_Vector;
                rhpx,ihpx : in Standard_Floating_Vectors.Link_to_Vector;
                rlpx,ilpx : in Standard_Floating_Vectors.Link_to_Vector ) is

    dd : double_double;

  begin
    for k in x'range loop
      dd := DoblDobl_Complex_Numbers.REAL_PART(x(k));
      rhpx(k) := hi_part(dd); rlpx(k) := lo_part(dd);
      dd := DoblDobl_Complex_Numbers.IMAG_PART(x(k));
      ihpx(k) := hi_part(dd); ilpx(k) := lo_part(dd);
    end loop;
  end Complex_Parts;

  procedure Complex_Parts
              ( x : in DoblDobl_Complex_VecVecs.VecVec;
                rhpx,ihpx : in Standard_Floating_VecVecs.Link_to_VecVec;
                rlpx,ilpx : in Standard_Floating_VecVecs.Link_to_VecVec ) is
  begin
    for k in x'range loop
      Complex_Parts(x(k),rhpx(k),ihpx(k),rlpx(k),ilpx(k));
    end loop;
  end Complex_Parts;

  procedure Complex_Parts
              ( x : in DoblDobl_Complex_VecVecs.Link_to_VecVec;
                rhpx,ihpx : in Standard_Floating_VecVecs.Link_to_VecVec;
                rlpx,ilpx : in Standard_Floating_VecVecs.Link_to_VecVec ) is
  begin
    for k in x'range loop
      Complex_Parts(x(k),rhpx(k),ihpx(k),rlpx(k),ilpx(k));
    end loop;
  end Complex_Parts;

  procedure Complex_Parts
              ( deg : in integer32;
                x : in DoblDobl_Complex_Vectors.Link_to_Vector;
                rhpx,ihpx : in Standard_Floating_Vectors.Link_to_Vector;
                rlpx,ilpx : in Standard_Floating_Vectors.Link_to_Vector ) is

    dd : double_double;

  begin
    for k in x'first..deg loop
      dd := DoblDobl_Complex_Numbers.REAL_PART(x(k));
      rhpx(k) := hi_part(dd); rlpx(k) := lo_part(dd);
      dd := DoblDobl_Complex_Numbers.IMAG_PART(x(k));
      ihpx(k) := hi_part(dd); ilpx(k) := lo_part(dd);
    end loop;
  end Complex_Parts;

  procedure Complex_Parts
              ( deg : in integer32;
                x : in DoblDobl_Complex_VecVecs.VecVec;
                rhpx,ihpx : in Standard_Floating_VecVecs.Link_to_VecVec;
                rlpx,ilpx : in Standard_Floating_VecVecs.Link_to_VecVec ) is
  begin
    for k in x'range loop
      Complex_Parts(deg,x(k),rhpx(k),ihpx(k),rlpx(k),ilpx(k));
    end loop;
  end Complex_Parts;

  procedure Complex_Parts
              ( deg : in integer32;
                x : in DoblDobl_Complex_VecVecs.Link_to_VecVec;
                rhpx,ihpx : in Standard_Floating_VecVecs.Link_to_VecVec;
                rlpx,ilpx : in Standard_Floating_VecVecs.Link_to_VecVec ) is
  begin
    for k in x'range loop
      Complex_Parts(deg,x(k),rhpx(k),ihpx(k),rlpx(k),ilpx(k));
    end loop;
  end Complex_Parts;

  procedure Complex_Merge
              ( rhpx,ihpx : in Standard_Floating_Vectors.Link_to_Vector;
                rlpx,ilpx : in Standard_Floating_Vectors.Link_to_Vector;
                cvx : in DoblDobl_Complex_Vectors.Link_to_Vector ) is

    rdd,idd : double_double;

  begin
    for k in cvx'range loop
      rdd := Create(rhpx(k),rlpx(k));
      idd := Create(ihpx(k),ilpx(k));
      cvx(k) := DoblDobl_Complex_Numbers.Create(rdd,idd);
    end loop;
  end Complex_Merge;

  procedure Complex_Merge
              ( rhpx,ihpx : in Standard_Floating_VecVecs.Link_to_VecVec;
                rlpx,ilpx : in Standard_Floating_VecVecs.Link_to_VecVec;
                cvx : in DoblDobl_Complex_VecVecs.Link_to_VecVec ) is
  begin
    for k in cvx'range loop
      Complex_Merge(rhpx(k),ihpx(k),rlpx(k),ilpx(k),cvx(k));
    end loop;
  end Complex_Merge;

  procedure Complex_Merge
              ( rhpx,ihpx : in Standard_Floating_VecVecs.Link_to_VecVec;
                rlpx,ilpx : in Standard_Floating_VecVecs.Link_to_VecVec;
                cvx : in DoblDobl_Complex_VecVecs.VecVec ) is
  begin
    for k in cvx'range loop
      Complex_Merge(rhpx(k),ihpx(k),rlpx(k),ilpx(k),cvx(k));
    end loop;
  end Complex_Merge;

  procedure Complex_Merge
              ( deg : in integer32;
                rhpx,ihpx : in Standard_Floating_Vectors.Link_to_Vector;
                rlpx,ilpx : in Standard_Floating_Vectors.Link_to_Vector;
                cvx : in DoblDobl_Complex_Vectors.Link_to_Vector ) is

    rdd,idd : double_double;

  begin
    for k in cvx'first..deg loop
      rdd := Create(rhpx(k),rlpx(k));
      idd := Create(ihpx(k),ilpx(k));
      cvx(k) := DoblDobl_Complex_Numbers.Create(rdd,idd);
    end loop;
  end Complex_Merge;

  procedure Complex_Merge
              ( deg : in integer32;
                rhpx,ihpx : in Standard_Floating_VecVecs.Link_to_VecVec;
                rlpx,ilpx : in Standard_Floating_VecVecs.Link_to_VecVec;
                cvx : in DoblDobl_Complex_VecVecs.Link_to_VecVec ) is
  begin
    for k in cvx'range loop
      Complex_Merge(deg,rhpx(k),ihpx(k),rlpx(k),ilpx(k),cvx(k));
    end loop;
  end Complex_Merge;

  procedure Complex_Merge
              ( deg : in integer32;
                rhpx,ihpx : in Standard_Floating_VecVecs.Link_to_VecVec;
                rlpx,ilpx : in Standard_Floating_VecVecs.Link_to_VecVec;
                cvx : in DoblDobl_Complex_VecVecs.VecVec ) is
  begin
    for k in cvx'range loop
      Complex_Merge(deg,rhpx(k),ihpx(k),rlpx(k),ilpx(k),cvx(k));
    end loop;
  end Complex_Merge;

  procedure Add ( zrehi : in Standard_Floating_Vectors.Link_to_Vector;
                  zimhi : in Standard_Floating_Vectors.Link_to_Vector;
                  zrelo : in Standard_Floating_Vectors.Link_to_Vector;
                  zimlo : in Standard_Floating_Vectors.Link_to_Vector;
                  xrehi : in Standard_Floating_Vectors.Link_to_Vector;
                  ximhi : in Standard_Floating_Vectors.Link_to_Vector;
                  xrelo : in Standard_Floating_Vectors.Link_to_Vector;
                  ximlo : in Standard_Floating_Vectors.Link_to_Vector;
                  yrehi : in Standard_Floating_Vectors.Link_to_Vector;
                  yimhi : in Standard_Floating_Vectors.Link_to_Vector;
                  yrelo : in Standard_Floating_Vectors.Link_to_Vector;
                  yimlo : in Standard_Floating_Vectors.Link_to_Vector) is

    dim : constant integer32 := zrehi'last;
    a,b,bb,s,s1,s2,t1,t2 : double_float;

  begin
    for k in 1..dim loop
     -- first sum the real parts
     -- Double_Double_Basics.two_sum(xrehi(k),yrehi(k),s1,s2);
      a := xrehi(k); b := yrehi(k);
      s1 := a + b; bb := s1 - a; s2 := (a - (s1 - bb)) + (b - bb);
     -- Double_Double_Basics.two_sum(xrelo(k),yrelo(k),t1,t2);
      a := xrelo(k); b := yrelo(k);
      t1 := a + b; bb := t1 - a; t2 := (a - (t1 - bb)) + (b - bb);
      s2 := s2 + t1;
     -- Double_Double_Basics.quick_two_sum(s1,s2,s1,s2);
      a := s1; b := s2; s1 := a + b; s2 := b - (s1 - a);
      s2 := s2 + t2;
     -- Double_Double_Basics.quick_two_sum(s1,s2,zrehi(k),zrelo(k));
      s := s1 + s2; zrehi(k) := s; zrelo(k) := s2 - (s - s1);
     -- then sum the imaginary parts
     -- Double_Double_Basics.two_sum(ximhi(k),yimhi(k),s1,s2);
      a := ximhi(k); b := yimhi(k);
      s1 := a + b; bb := s1 - a; s2 := (a - (s1 - bb)) + (b - bb);
     -- Double_Double_Basics.two_sum(ximlo(k),yimlo(k),t1,t2);
      a := ximlo(k); b := yimlo(k);
      t1 := a + b; bb := t1 - a; t2 := (a - (t1 - bb)) + (b - bb);
      s2 := s2 + t1;
     -- Double_Double_Basics.quick_two_sum(s1,s2,s1,s2);
      a := s1; b := s2; s1 := a + b; s2 := b - (s1 - a);
      s2 := s2 + t2;
     -- Double_Double_Basics.quick_two_sum(s1,s2,zimhi(k),zimlo(k));
      s := s1 + s2; zimhi(k) := s; zimlo(k) := s2 - (s - s1);
    end loop;
  end Add;

  procedure Update ( zrehi : in Standard_Floating_Vectors.Link_to_Vector;
                     zimhi : in Standard_Floating_Vectors.Link_to_Vector;
                     zrelo : in Standard_Floating_Vectors.Link_to_Vector;
                     zimlo : in Standard_Floating_Vectors.Link_to_Vector;
                     xrehi : in Standard_Floating_Vectors.Link_to_Vector;
                     ximhi : in Standard_Floating_Vectors.Link_to_Vector;
                     xrelo : in Standard_Floating_Vectors.Link_to_Vector;
                     ximlo : in Standard_Floating_Vectors.Link_to_Vector ) is

    dim : constant integer32 := zrehi'last;
    a,b,bb,s,s1,s2,t1,t2 : double_float;

  begin
    for k in zrehi'first..dim loop
     -- first sum the real parts
     -- Double_Double_Basics.two_sum(xrehi(k),zrehi(k),s1,s2);
      a := xrehi(k); b := zrehi(k);
      s1 := a + b; bb := s1 - a; s2 := (a - (s1 - bb)) + (b - bb);
     -- Double_Double_Basics.two_sum(xrelo(k),zrelo(k),t1,t2);
      a := xrelo(k); b := zrelo(k);
      t1 := a + b; bb := t1 - a; t2 := (a - (t1 - bb)) + (b - bb);
      s2 := s2 + t1;
     -- Double_Double_Basics.quick_two_sum(s1,s2,s1,s2);
      a := s1; b := s2; s1 := a + b; s2 := b - (s1 - a);
      s2 := s2 + t2;
     -- Double_Double_Basics.quick_two_sum(s1,s2,zrehi(k),zrelo(k));
      s := s1 + s2; zrehi(k) := s; zrelo(k) := s2 - (s - s1);
     -- then sum the imaginary parts
     -- Double_Double_Basics.two_sum(ximhi(k),zimhi(k),s1,s2);
      a := ximhi(k); b := zimhi(k);
      s1 := a + b; bb := s1 - a; s2 := (a - (s1 - bb)) + (b - bb);
     -- Double_Double_Basics.two_sum(ximlo(k),zimlo(k),t1,t2);
      a := ximlo(k); b := zimlo(k);
      t1 := a + b; bb := t1 - a; t2 := (a - (t1 - bb)) + (b - bb);
      s2 := s2 + t1;
     -- Double_Double_Basics.quick_two_sum(s1,s2,s1,s2);
      a := s1; b := s2; s1 := a + b; s2 := b - (s1 - a);
      s2 := s2 + t2;
     -- Double_Double_Basics.quick_two_sum(s1,s2,zimhi(k),zimlo(k));
      s := s1 + s2; zimhi(k) := s; zimlo(k) := s2 - (s - s1);
    end loop;
  end Update;

  procedure Update_Product
                ( zrehi,zimhi,zrelo,zimlo : in out double_float;
                  xrehi,ximhi,xrelo,ximlo : in double_float;
                  yrehi,yimhi,yrelo,yimlo : in double_float ) is

    QD_SPLITTER : constant double_float := 134217729.0; -- 2^27 + 1
    QD_SPLIT_THRESH : constant double_float := 6.69692879491417e+299; -- 2^996

    a,b,bb,s1,s2,t1,t2,p1,p2,atemp,btemp,aa : double_float;
    a_hiprod,a_loprod,b_hiprod,b_loprod : double_float;
    c_hiprod,c_loprod,d_hiprod,d_loprod : double_float;
    a1_hi,a1_lo,b1_hi,b1_lo : double_float; -- 1st pair of splits
    a2_hi,a2_lo,b2_hi,b2_lo : double_float; -- 2nd pair of splits 

  begin
   -- prod := x*y, as complex product
   -- reprod := xre*yre - xim*yim, real part of product
   -- improd := xre*yim + xim*yre, imaginary part of product
   -- reprod has high part in rehiprod, low part in reloprod
   -- xre*yre is stored in a_hiprod and a_loprod
   -- Double_Double_Basics.two_prod(xrehi,yrehi,p1,p2);
   -- Double_Double_Basics.split(xrehi,a_hi,a_lo);
    if ( xrehi > QD_SPLIT_THRESH or xrehi < -QD_SPLIT_THRESH ) then
      aa := xrehi*3.7252902984619140625E-09;  -- 2^-28
      atemp := QD_SPLITTER * aa;
      a1_hi := atemp - (atemp - aa);
      a1_lo := xrehi - a1_hi;
      a1_hi := a1_hi*268435456.0;  -- 2^28
      a1_lo := 268435456.0;     -- 2^28
    else
      atemp := QD_SPLITTER * xrehi;
      a1_hi := atemp - (atemp - xrehi);
      a1_lo := xrehi - a1_hi;
    end if;
   -- Double_Double_Basics.split(yrehi,b1_hi,b1_lo);
    if ( yrehi > QD_SPLIT_THRESH or yrehi < -QD_SPLIT_THRESH ) then
      bb := yrehi*3.7252902984619140625E-09;  -- 2^-28
      btemp := QD_SPLITTER * bb;
      b1_hi := btemp - (btemp - bb);
      b1_lo := yrehi - b1_hi;
      b1_hi := b1_hi*268435456.0;  -- 2^28
      b1_lo := 268435456.0;     -- 2^28
    else
      btemp := QD_SPLITTER * yrehi;
      b1_hi := btemp - (btemp - yrehi);
      b1_lo := yrehi - b1_hi;
    end if;
   -- xim*yim is stored in b_hiprod and b_loprod
   -- Double_Double_Basics.two_prod(ximhi,yimhi,p1,p2);
   -- Double_Double_Basics.split(ximhi,a_hi,a_lo);
    if ( ximhi > QD_SPLIT_THRESH or ximhi < -QD_SPLIT_THRESH ) then
      aa := ximhi*3.7252902984619140625E-09;  -- 2^-28
      atemp := QD_SPLITTER * aa;
      a2_hi := atemp - (atemp - aa);
      a2_lo := ximhi - a2_hi;
      a2_hi := a2_hi*268435456.0;  -- 2^28
      a2_lo := 268435456.0;     -- 2^28
    else
      atemp := QD_SPLITTER * ximhi;
      a2_hi := atemp - (atemp - ximhi);
      a2_lo := ximhi - a2_hi;
    end if;
   -- Double_Double_Basics.split(yimhi,b_hi,b_lo);
    if ( yimhi > QD_SPLIT_THRESH or yimhi < -QD_SPLIT_THRESH ) then
      bb := yimhi*3.7252902984619140625E-09;  -- 2^-28
      btemp := QD_SPLITTER * bb;
      b2_hi := btemp - (btemp - bb);
      b2_lo := yimhi - b2_hi;
      b2_hi := b2_hi*268435456.0;  -- 2^28
      b2_lo := 268435456.0;     -- 2^28
    else
      btemp := QD_SPLITTER * yimhi;
      b2_hi := btemp - (btemp - yimhi);
      b2_lo := yimhi - b2_hi;
    end if;
    p1 := xrehi*yrehi;
    p2 := ((a1_hi*b1_hi - p1) + a1_hi*b1_lo + a1_lo*b1_hi) + a1_lo*b1_lo;
    p2 := p2 + (xrehi * yrelo + xrelo * yrehi);
   -- Double_Double_Basics.quick_two_sum(p1,p2,a_hiprod,a_loprod);
    a_hiprod := p1 + p2; a_loprod := p2 - (a_hiprod - p1);
    p1 := ximhi*yimhi;
    p2 := ((a2_hi*b2_hi - p1) + a2_hi*b2_lo + a2_lo*b2_hi) + a2_lo*b2_lo;
    p2 := p2 + (ximhi * yimlo + ximlo * yimhi);
   -- Double_Double_Basics.quick_two_sum(p1,p2,b_hiprod,b_loprod);
    b_hiprod := p1 + p2; b_loprod := p2 - (b_hiprod - p1);
   -- xre*yim is stored in c_hiprod and c_loprod
   -- Double_Double_Basics.two_prod(xrehi,yimhi,p1,p2);
    p1 := xrehi*yimhi;
   -- Double_Double_Basics.split(xrehi,a_hi,a_lo); -> a1_hi, a1_lo
   -- Double_Double_Basics.split(yimhi,b_hi,b_lo); -> b2_hi, b2_lo
    p2 := ((a1_hi*b2_hi - p1) + a1_hi*b2_lo + a1_lo*b2_hi) + a1_lo*b2_lo;
    p2 := p2 + (xrehi * yimlo + xrelo * yimhi);
   -- Double_Double_Basics.quick_two_sum(p1,p2,c_hiprod,c_loprod);
    c_hiprod := p1 + p2; c_loprod := p2 - (c_hiprod - p1);
   -- xim*yre is stored in hiprod and loprod
   -- Double_Double_Basics.two_prod(ximhi,yrehi,p1,p2);
    p1 := ximhi*yrehi;
   -- Double_Double_Basics.split(ximhi,a_hi,a_lo); -> a2_hi, a2_lo
   -- Double_Double_Basics.split(yrehi,b_hi,b_lo); -> b1_hi, b1_lo
    p2 := ((a2_hi*b1_hi - p1) + a2_hi*b1_lo + a2_lo*b1_hi) + a2_lo*b1_lo;
    p2 := p2 + (ximhi * yrelo + ximlo * yrehi);
   -- Double_Double_Basics.quick_two_sum(p1,p2,d_hiprod,d_loprod);
    d_hiprod := p1 + p2; d_loprod := p2 - (d_hiprod - p1);
   -- add xre*yre to zre, with high part in zrehi, low part in zrelo
   -- Double_Double_Basics.two_sum(a_hiprod,zrehi,s1,s2);
    s1 := a_hiprod + zrehi; bb := s1 - a_hiprod;
    s2 := (a_hiprod - (s1 - bb)) + (zrehi - bb);
   -- Double_Double_Basics.two_sum(a_loprod,zrelo,t1,t2);
    t1 := a_loprod + zrelo; bb := t1 - a_loprod;
    t2 := (a_loprod - (t1 - bb)) + (zrelo - bb);
    s2 := s2 + t1;
   -- Double_Double_Basics.quick_two_sum(s1,s2,s1,s2);
    a := s1; b := s2; s1 := a + b; s2 := b - (s1 - a);
    s2 := s2 + t2;
   -- Double_Double_Basics.quick_two_sum(s1,s2,zrehi,zrelo);
    zrehi := s1 + s2; zrelo := s2 - (zrehi - s1);
   -- subtract from zre xim(k)*yim(k), stored in b_hiprod, b_loprod
   -- Double_Double_Basics.two_diff(zrehi,b_hiprod,s1,s2);
    s1 := zrehi - b_hiprod; bb := s1 - zrehi;
    s2 := (zrehi - (s1 - bb)) - (b_hiprod + bb);
   -- Double_Double_Basics.two_diff(zrelo,b_loprod,t1,t2);
    t1 := zrelo - b_loprod; bb := t1 - zrelo;
    t2 := (zrelo - (t1 - bb)) - (b_loprod + bb);
    s2 := s2 + t1;
   -- Double_Double_Basics.quick_two_sum(s1,s2,s1,s2);
    a := s1; b := s2; s1 := a + b; s2 := b - (s1 - a);
    s2 := s2 + t2;
   -- Double_Double_Basics.quick_two_sum(s1,s2,zrehi,zrelo);
    zrehi := s1 + s2; zrelo := s2 - (zrehi - s1);
   -- add xre(k)*yim(k) to zim, with high part in zimhi, low part in zimlo
   -- Double_Double_Basics.two_sum(c_hiprod,zimhi,s1,s2);
    s1 := c_hiprod + zimhi; bb := s1 - c_hiprod;
    s2 := (c_hiprod - (s1 - bb)) + (zimhi - bb);
   -- Double_Double_Basics.two_sum(c_loprod,zimlo,t1,t2);
    t1 := c_loprod + zimlo; bb := t1 - c_loprod;
    t2 := (c_loprod - (t1 - bb)) + (zimlo - bb);
    s2 := s2 + t1;
   -- Double_Double_Basics.quick_two_sum(s1,s2,s1,s2);
    a := s1; b := s2; s1 := a + b; s2 := b - (s1 - a);
    s2 := s2 + t2;
   -- Double_Double_Basics.quick_two_sum(s1,s2,zimhi,zimlo);
    zimhi := s1 + s2; zimlo := s2 - (zimhi - s1);
   -- add to zim xim(k)*yre(k), stored in d_hiprod, d_loprod
   -- Double_Double_Basics.two_sum(d_hiprod,zimhi,s1,s2);
    s1 := d_hiprod + zimhi; bb := s1 - d_hiprod;
    s2 := (d_hiprod - (s1 - bb)) + (zimhi - bb);
   -- Double_Double_Basics.two_sum(d_loprod,zimlo,t1,t2);
    t1 := d_loprod + zimlo; bb := t1 - d_loprod;
    t2 := (d_loprod - (t1 - bb)) + (zimlo - bb);
    s2 := s2 + t1;
   -- Double_Double_Basics.quick_two_sum(s1,s2,s1,s2);
    a := s1; b := s2; s1 := a + b; s2 := b - (s1 - a);
    s2 := s2 + t2;
   -- Double_Double_Basics.quick_two_sum(s1,s2,zimhi,zimlo);
    zimhi := s1 + s2; zimlo := s2 - (zimhi - s1);
  end Update_Product;

  procedure Inner_Product
                ( zrehi,zimhi,zrelo,zimlo : out double_float;
                  xrehi : in Standard_Floating_Vectors.Link_to_Vector;
                  ximhi : in Standard_Floating_Vectors.Link_to_Vector;
                  xrelo : in Standard_Floating_Vectors.Link_to_Vector;
                  ximlo : in Standard_Floating_Vectors.Link_to_Vector;
                  yrehi : in Standard_Floating_Vectors.Link_to_Vector;
                  yimhi : in Standard_Floating_Vectors.Link_to_Vector;
                  yrelo : in Standard_Floating_Vectors.Link_to_Vector;
                  yimlo : in Standard_Floating_Vectors.Link_to_Vector) is

    dim : constant integer32 := xrehi'last;

  begin
    zrehi := 0.0; zimhi := 0.0; zrelo := 0.0; zimlo := 0.0;
    for k in 1..dim loop
      Update_Product(zrehi,zimhi,zrelo,zimlo,
                     xrehi(k),ximhi(k),xrelo(k),ximlo(k),
                     yrehi(k),yimhi(k),yrelo(k),yimlo(k));
    end loop;
  end Inner_Product;

  procedure Multiply
              ( xrehi : in Standard_Floating_Vectors.Link_to_Vector;
                ximhi : in Standard_Floating_Vectors.Link_to_Vector;
                xrelo : in Standard_Floating_Vectors.Link_to_Vector;
                ximlo : in Standard_Floating_Vectors.Link_to_Vector;
                yrehi : in Standard_Floating_Vectors.Link_to_Vector;
                yimhi : in Standard_Floating_Vectors.Link_to_Vector;
                yrelo : in Standard_Floating_Vectors.Link_to_Vector;
                yimlo : in Standard_Floating_Vectors.Link_to_Vector;
                zrehi : in Standard_Floating_Vectors.Link_to_Vector;
                zimhi : in Standard_Floating_Vectors.Link_to_Vector;
                zrelo : in Standard_Floating_Vectors.Link_to_Vector;
                zimlo : in Standard_Floating_Vectors.Link_to_Vector ) is

    deg : constant integer32 := xrehi'last;

  begin
   -- product(0) := first(0)*second(0);
    zrehi(0) := 0.0; zimhi(0) := 0.0; zrelo(0) := 0.0; zimlo(0) := 0.0;
    Update_Product(zrehi(0),zimhi(0),zrelo(0),zimlo(0),
                   xrehi(0),ximhi(0),xrelo(0),ximlo(0),
                   yrehi(0),yimhi(0),yrelo(0),yimlo(0));
    for k in 1..deg loop
     -- product(k) := first(0)*second(k);
      zrehi(k) := 0.0; zimhi(k) := 0.0; zrelo(k) := 0.0; zimlo(k) := 0.0;
      Update_Product(zrehi(k),zimhi(k),zrelo(k),zimlo(k),
                     xrehi(0),ximhi(0),xrelo(0),ximlo(0),
                     yrehi(k),yimhi(k),yrelo(k),yimlo(k));
      for i in 1..k loop
       -- product(k) := product(k) + first(i)*second(k-i);
        Update_Product(zrehi(k),zimhi(k),zrelo(k),zimlo(k),
                       xrehi(i),ximhi(i),xrelo(i),ximlo(i),
                       yrehi(k-i),yimhi(k-i),yrelo(k-i),yimlo(k-i));
      end loop;
    end loop;
  end Multiply;

end DoblDobl_Vector_Splitters;
