with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with Timing_Package;                     use Timing_Package;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Double_Double_Basics;
with Double_Double_Numbers;              use Double_Double_Numbers;
with DoblDobl_Complex_Numbers;           use DoblDobl_Complex_Numbers;
with DoblDobl_Complex_Numbers_io;        use DoblDobl_Complex_Numbers_io;
with Standard_Floating_Vectors;
with DoblDobl_Complex_Vectors;
with DoblDobl_Complex_Vectors_io;        use DoblDobl_Complex_Vectors_io;
with DoblDobl_Random_Vectors;

procedure ts_perfddvc is

-- DESCRIPTION :
--   Tests the performance on double double complex vector computations.

  procedure Split ( x : in Complex_Number;
                    rehi,imhi,relo,imlo : out double_float ) is

  -- DESCRIPTION :
  --   Splits the complex number in real and imaginary parts,
  --   in high and low parts.

  -- ON ENTRY :
  --   x        a double double complex number.

  -- ON RETURN :
  --   rehi     high part of the real part of x;
  --   imhi     high part of the imaginary part of x;
  --   relo     low part of the real part of x;
  --   imlo     low part of the imaginary part of x.

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
  
  -- DESCRIPTION :
  --   Splits the vector v into real and imaginary,
  --   into high and low parts.

  -- REQUIRED : v'range = rehi'range = imhi'range = relo'range = imlo'range.

  -- ON ENTRY :
  --   v        a vector of double double complex numbers.

  -- ON RETURN :
  --   rehi     high parts of the real parts of the numbers in v;
  --   imhi     high parts of the imaginary parts of the numbers in v;
  --   relo     low parts of the real parts of the numbers in v;
  --   imlo     low parts of the imaginary parts of the numbers in v.

  begin
    for k in v'range loop
      Split(v(k),rehi(k),imhi(k),relo(k),imlo(k));
    end loop;
  end Split;

  procedure Merge ( x : out Complex_Number;
                    rehi,imhi,relo,imlo : in double_float ) is

  -- DESCRIPTION :
  --   Merges the high and low parts of real and imaginary parts
  --   into one complex number.

  -- ON ENTRY :
  --   rehi     high part of the real part;
  --   imhi     high part of the imaginary part;
  --   relo     low part of the real part;
  --   imlo     low part of the imaginary part.

  -- ON RETURN :
  --   x        complex number with high part of real part in rehi,
  --            high part of imaginary part in imhi,
  --            low part of real part in imlo,
  --            low part of imaginary part in imlo.

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

  -- DESCRIPTION :
  --   Merges the high and low parts of real and imaginary parts
  --   into one complex vector v.

  -- REQUIRED : v'range = rehi'range = imhi'range = relo'range = imlo'range.

  -- ON ENTRY :
  --   rehi     high parts of the real parts for the numbers in v;
  --   imhi     high parts of the imaginary parts or the numbers in v;
  --   relo     low parts of the real parts for the numbers in v;
  --   imlo     low parts of the imaginary parts for the numbers in v.

  -- ON RETURN :
  --   v        a vector of double double complex numbers, with
  --            high part of real parts in rehi,
  --            high part of imaginary parts in imhi,
  --            low part of real parts in relo, and
  --            low part of imaginary parts in imlo,

  begin
    for k in v'range loop
      Merge(v(k),rehi(k),imhi(k),relo(k),imlo(k));
    end loop;
  end Merge;

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

  -- DESCRIPTION :
  --   Adds two double double complex vectors x and y to form
  --   the result z, in 4-vector representation.

  -- REQUIRED : all vectors have the same range.

  -- ON ENTRY :
  --   xrehi      high parts of the real parts for the numbers in x;
  --   ximhi      high parts of the imaginary parts or the numbers in x;
  --   xrelo      low parts of the real parts for the numbers in x;
  --   ximlo      low parts of the imaginary parts for the numbers in x;
  --   yrehi      high parts of the real parts for the numbers in y;
  --   yimhi      high parts of the imaginary parts or the numbers in y;
  --   yrelo      low parts of the real parts for the numbers in y;
  --   yimlo      low parts of the imaginary parts for the numbers in y.

  -- ON RETURN :
  --   zrehi      high parts of the real parts for the numbers in z;
  --   zimhi      high parts of the imaginary parts for the numbers in z;
  --   zrelo      low parts of the real parts for the numbers in z;
  --   zimlo      low parts of the imaginary parts for the numbers in z.

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

  function Inner_Product
             ( x,y : DoblDobl_Complex_Vectors.Vector )
             return Complex_Number is

  -- DESCRIPTION :
  --   Returns the inner product of the vectors x and y,
  --   without taking complex conjugates.

  -- REQUIRED : x'range = y'range.

    res : Complex_Number := Create(integer(0));

  begin
    for k in x'range loop
      res := res + x(k)*y(k);
    end loop;
    return res;   
  end Inner_Product;

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

  -- DESCRIPTION :
  --   Computes the inner product of two double double complex vectors
  --   x and y to form the result z, in 4-vector representation.

  -- REQUIRED : all vectors have the same range.

  -- ON ENTRY :
  --   xrehi      high parts of the real parts for the numbers in x;
  --   ximhi      high parts of the imaginary parts or the numbers in x;
  --   xrelo      low parts of the real parts for the numbers in x;
  --   ximlo      low parts of the imaginary parts for the numbers in x;
  --   yrehi      high parts of the real parts for the numbers in y;
  --   yimhi      high parts of the imaginary parts or the numbers in y;
  --   yrelo      low parts of the real parts for the numbers in y;
  --   yimlo      low parts of the imaginary parts for the numbers in y.

  -- ON RETURN :
  --   zrehi      high parts of the real parts of the inner product;
  --   zimhi      high parts of the imaginary parts of the inner product.
  --   zrelo      low parts of the real parts of the inner product;
  --   zimlo      low parts of the imaginary parts of the inner product.

    QD_SPLITTER : constant double_float := 134217729.0; -- 2^27 + 1
    QD_SPLIT_THRESH : constant double_float := 6.69692879491417e+299; -- 2^996

    dim : constant integer32 := xrehi'last;
    a,b,bb,s1,s2,t1,t2,p1,p2,atemp,btemp,aa : double_float;
    a_hiprod,a_loprod,b_hiprod,b_loprod : double_float;
    c_hiprod,c_loprod,d_hiprod,d_loprod : double_float;
    xrehik,xrelok,ximhik,ximlok : double_float;
    yrehik,yrelok,yimhik,yimlok : double_float;
    a1_hi,a1_lo,b1_hi,b1_lo : double_float; -- 1st pair of splits
    a2_hi,a2_lo,b2_hi,b2_lo : double_float; -- 2nd pair of splits 

  begin
    zrehi := 0.0; zimhi := 0.0; zrelo := 0.0; zimlo := 0.0;
    for k in 1..dim loop
     -- prod := x(k)*y(k), as complex product
     -- reprod := xre(k)*yre(k) - xim(k)*yim(k), real part of product
     -- improd := xre(k)*yim(k) + xim(k)*yre(k), imaginary part of product
     -- reprod has high part in rehiprod, low part in reloprod
      xrehik := xrehi(k); xrelok := xrelo(k);
      yrehik := yrehi(k); yrelok := yrelo(k);
      ximhik := ximhi(k); ximlok := ximlo(k);
      yimhik := yimhi(k); yimlok := yimlo(k);
     -- xre(k)*yre(k) is stored in a_hiprod and a_loprod
     -- Double_Double_Basics.two_prod(xrehik,yrehik,p1,p2);
     -- Double_Double_Basics.split(xrehik,a_hi,a_lo);
      if ( xrehik > QD_SPLIT_THRESH or xrehik < -QD_SPLIT_THRESH ) then
        aa := xrehik*3.7252902984619140625E-09;  -- 2^-28
        atemp := QD_SPLITTER * aa;
        a1_hi := atemp - (atemp - aa);
        a1_lo := xrehik - a1_hi;
        a1_hi := a1_hi*268435456.0;  -- 2^28
        a1_lo := 268435456.0;     -- 2^28
      else
        atemp := QD_SPLITTER * xrehik;
        a1_hi := atemp - (atemp - xrehik);
        a1_lo := xrehik - a1_hi;
      end if;
     -- Double_Double_Basics.split(yrehik,b1_hi,b1_lo);
      if ( yrehik > QD_SPLIT_THRESH or yrehik < -QD_SPLIT_THRESH ) then
        bb := yrehik*3.7252902984619140625E-09;  -- 2^-28
        btemp := QD_SPLITTER * bb;
        b1_hi := btemp - (btemp - bb);
        b1_lo := yrehik - b1_hi;
        b1_hi := b1_hi*268435456.0;  -- 2^28
        b1_lo := 268435456.0;     -- 2^28
      else
        btemp := QD_SPLITTER * yrehik;
        b1_hi := btemp - (btemp - yrehik);
        b1_lo := yrehik - b1_hi;
      end if;
     -- xim(k)*yim(k) is stored in b_hiprod and b_loprod
     -- Double_Double_Basics.two_prod(ximhik,yimhik,p1,p2);
     -- Double_Double_Basics.split(ximhik,a_hi,a_lo);
      if ( ximhik > QD_SPLIT_THRESH or ximhik < -QD_SPLIT_THRESH ) then
        aa := ximhik*3.7252902984619140625E-09;  -- 2^-28
        atemp := QD_SPLITTER * aa;
        a2_hi := atemp - (atemp - aa);
        a2_lo := ximhik - a2_hi;
        a2_hi := a2_hi*268435456.0;  -- 2^28
        a2_lo := 268435456.0;     -- 2^28
      else
        atemp := QD_SPLITTER * ximhik;
        a2_hi := atemp - (atemp - ximhik);
        a2_lo := ximhik - a2_hi;
      end if;
     -- Double_Double_Basics.split(yimhik,b_hi,b_lo);
      if ( yimhik > QD_SPLIT_THRESH or yimhik < -QD_SPLIT_THRESH ) then
        bb := yimhik*3.7252902984619140625E-09;  -- 2^-28
        btemp := QD_SPLITTER * bb;
        b2_hi := btemp - (btemp - bb);
        b2_lo := yimhik - b2_hi;
        b2_hi := b2_hi*268435456.0;  -- 2^28
        b2_lo := 268435456.0;     -- 2^28
      else
        btemp := QD_SPLITTER * yimhik;
        b2_hi := btemp - (btemp - yimhik);
        b2_lo := yimhik - b2_hi;
      end if;
      p1 := xrehik*yrehik;
      p2 := ((a1_hi*b1_hi - p1) + a1_hi*b1_lo + a1_lo*b1_hi) + a1_lo*b1_lo;
      p2 := p2 + (xrehik * yrelok + xrelok * yrehik);
     -- Double_Double_Basics.quick_two_sum(p1,p2,a_hiprod,a_loprod);
      a_hiprod := p1 + p2; a_loprod := p2 - (a_hiprod - p1);
      p1 := ximhik*yimhik;
      p2 := ((a2_hi*b2_hi - p1) + a2_hi*b2_lo + a2_lo*b2_hi) + a2_lo*b2_lo;
      p2 := p2 + (ximhik * yimlok + ximlok * yimhik);
     -- Double_Double_Basics.quick_two_sum(p1,p2,b_hiprod,b_loprod);
      b_hiprod := p1 + p2; b_loprod := p2 - (b_hiprod - p1);
     -- xre(k)*yim(k) is stored in c_hiprod and c_loprod
     -- Double_Double_Basics.two_prod(xrehik,yimhik,p1,p2);
      p1 := xrehik*yimhik;
     -- Double_Double_Basics.split(xrehik,a_hi,a_lo); -> a1_hi, a1_lo
     -- Double_Double_Basics.split(yimhik,b_hi,b_lo); -> b2_hi, b2_lo
      p2 := ((a1_hi*b2_hi - p1) + a1_hi*b2_lo + a1_lo*b2_hi) + a1_lo*b2_lo;
      p2 := p2 + (xrehik * yimlok + xrelok * yimhik);
     -- Double_Double_Basics.quick_two_sum(p1,p2,c_hiprod,c_loprod);
      c_hiprod := p1 + p2; c_loprod := p2 - (c_hiprod - p1);
     -- xim(k)*yre(k) is stored in hiprod and loprod
     -- Double_Double_Basics.two_prod(ximhik,yrehik,p1,p2);
      p1 := ximhik*yrehik;
     -- Double_Double_Basics.split(ximhik,a_hi,a_lo); -> a2_hi, a2_lo
     -- Double_Double_Basics.split(yrehik,b_hi,b_lo); -> b1_hi, b1_lo
      p2 := ((a2_hi*b1_hi - p1) + a2_hi*b1_lo + a2_lo*b1_hi) + a2_lo*b1_lo;
      p2 := p2 + (ximhik * yrelok + ximlok * yrehik);
     -- Double_Double_Basics.quick_two_sum(p1,p2,d_hiprod,d_loprod);
      d_hiprod := p1 + p2; d_loprod := p2 - (d_hiprod - p1);
     -- add xre(k)*yre(k) to zre, with high part in zrehi, low part in zrelo
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
    end loop;
  end Inner_Product;

  procedure Test_Add ( dim : in integer32 ) is

  -- DESCRIPTION :
  --   Generates two random vectors and tests their sum.

    x : constant DoblDobl_Complex_Vectors.Vector(1..dim)
      := DoblDobl_Random_Vectors.Random_Vector(1,dim);
    y : constant DoblDobl_Complex_Vectors.Vector(1..dim)
      := DoblDobl_Random_Vectors.Random_Vector(1,dim);
    z1,z2 : DoblDobl_Complex_Vectors.Vector(1..dim);
    xrehi : constant Standard_Floating_Vectors.Vector(1..dim)
          := (1..dim => 0.0);
    ximhi : constant Standard_Floating_Vectors.Vector(1..dim)
          := (1..dim => 0.0);
    xrelo : constant Standard_Floating_Vectors.Vector(1..dim)
          := (1..dim => 0.0);
    ximlo : constant Standard_Floating_Vectors.Vector(1..dim)
          := (1..dim => 0.0);
    yrehi : constant Standard_Floating_Vectors.Vector(1..dim)
          := (1..dim => 0.0);
    yimhi : constant Standard_Floating_Vectors.Vector(1..dim)
          := (1..dim => 0.0);
    yrelo : constant Standard_Floating_Vectors.Vector(1..dim)
          := (1..dim => 0.0);
    yimlo : constant Standard_Floating_Vectors.Vector(1..dim)
          := (1..dim => 0.0);
    zrehi : constant Standard_Floating_Vectors.Vector(1..dim)
          := (1..dim => 0.0);
    zimhi : constant Standard_Floating_Vectors.Vector(1..dim)
          := (1..dim => 0.0);
    zrelo : constant Standard_Floating_Vectors.Vector(1..dim)
          := (1..dim => 0.0);
    zimlo : constant Standard_Floating_Vectors.Vector(1..dim)
          := (1..dim => 0.0);
    xrh : constant Standard_Floating_Vectors.Link_to_Vector
        := new Standard_Floating_Vectors.Vector'(xrehi);
    xih : constant Standard_Floating_Vectors.Link_to_Vector
        := new Standard_Floating_Vectors.Vector'(ximhi);
    xrl : constant Standard_Floating_Vectors.Link_to_Vector
        := new Standard_Floating_Vectors.Vector'(xrelo);
    xil : constant Standard_Floating_Vectors.Link_to_Vector
        := new Standard_Floating_Vectors.Vector'(ximlo);
    yrh : constant Standard_Floating_Vectors.Link_to_Vector
        := new Standard_Floating_Vectors.Vector'(yrehi);
    yih : constant Standard_Floating_Vectors.Link_to_Vector
        := new Standard_Floating_Vectors.Vector'(yimhi);
    yrl : constant Standard_Floating_Vectors.Link_to_Vector
        := new Standard_Floating_Vectors.Vector'(yrelo);
    yil : constant Standard_Floating_Vectors.Link_to_Vector
        := new Standard_Floating_Vectors.Vector'(yimlo);
    zrh : constant Standard_Floating_Vectors.Link_to_Vector
        := new Standard_Floating_Vectors.Vector'(yrehi);
    zih : constant Standard_Floating_Vectors.Link_to_Vector
        := new Standard_Floating_Vectors.Vector'(yimhi);
    zrl : constant Standard_Floating_Vectors.Link_to_Vector
        := new Standard_Floating_Vectors.Vector'(yrelo);
    zil : constant Standard_Floating_Vectors.Link_to_Vector
        := new Standard_Floating_Vectors.Vector'(yimlo);

    use DoblDobl_Complex_Vectors;

  begin
    z1 := x + y;
    put_line("The sum of two random vectors :"); put_line(z1);
    Split(x,xrh.all,xih.all,xrl.all,xil.all);
    Split(y,yrh.all,yih.all,yrl.all,yil.all);
    Add(zrh,zih,zrl,zil,xrh,xih,xrl,xil,yrh,yih,yrl,yil);
    Merge(z2,zrh.all,zih.all,zrl.all,zil.all);
    put_line("The recomputed sum :"); put_line(z2);
  end Test_Add;

  procedure Test_Inner_Product ( dim : in integer32 ) is

  -- DESCRIPTION :
  --   Generates two random vectors and tests their inner product.

    x : constant DoblDobl_Complex_Vectors.Vector(1..dim)
      := DoblDobl_Random_Vectors.Random_Vector(1,dim);
    y : constant DoblDobl_Complex_Vectors.Vector(1..dim)
      := DoblDobl_Random_Vectors.Random_Vector(1,dim);
    z1,z2 : Complex_Number;
    xrehi,ximhi,xrelo,ximlo : Standard_Floating_Vectors.Vector(1..dim);
    yrehi,yimhi,yrelo,yimlo : Standard_Floating_Vectors.Vector(1..dim);
    xrh,xih,xrl,xil : Standard_Floating_Vectors.Link_to_Vector;
    yrh,yih,yrl,yil : Standard_Floating_Vectors.Link_to_Vector;
    zrehi,zimhi,zrelo,zimlo : double_float;

  begin
    z1 := Inner_Product(x,y);
    put_line("The inner product of two random vectors :");
    put(z1); new_line;
    Split(x,xrehi,ximhi,xrelo,ximlo);
    Split(y,yrehi,yimhi,yrelo,yimlo);
    xrh := new Standard_Floating_Vectors.Vector'(xrehi);
    xih := new Standard_Floating_Vectors.Vector'(ximhi);
    xrl := new Standard_Floating_Vectors.Vector'(xrelo);
    xil := new Standard_Floating_Vectors.Vector'(ximlo);
    yrh := new Standard_Floating_Vectors.Vector'(yrehi);
    yih := new Standard_Floating_Vectors.Vector'(yimhi);
    yrl := new Standard_Floating_Vectors.Vector'(yrelo);
    yil := new Standard_Floating_Vectors.Vector'(yimlo);
    Inner_Product(zrehi,zimhi,zrelo,zimlo,xrh,xih,xrl,xil,yrh,yih,yrl,yil);
    Merge(z2,zrehi,zimhi,zrelo,zimlo);
    put_line("The recomputed inner product :"); put(z2); new_line;
  end Test_Inner_Product;

  procedure Timing_Add ( dim,frq : in integer32 ) is

  -- DESCRIPTION :
  --   Generates random vectors and times their sum,
  --   for a frequency equal to frq.

    timer : Timing_Widget;
    x : constant DoblDobl_Complex_Vectors.Vector(1..dim)
      := DoblDobl_Random_Vectors.Random_Vector(1,dim);
    y : constant DoblDobl_Complex_Vectors.Vector(1..dim)
      := DoblDobl_Random_Vectors.Random_Vector(1,dim);
    z1,z2 : DoblDobl_Complex_Vectors.Vector(1..dim);
    xrehi : constant Standard_Floating_Vectors.Vector(1..dim)
          := (1..dim => 0.0);
    ximhi : constant Standard_Floating_Vectors.Vector(1..dim)
          := (1..dim => 0.0);
    xrelo : constant Standard_Floating_Vectors.Vector(1..dim)
          := (1..dim => 0.0);
    ximlo : constant Standard_Floating_Vectors.Vector(1..dim)
          := (1..dim => 0.0);
    yrehi : constant Standard_Floating_Vectors.Vector(1..dim)
          := (1..dim => 0.0);
    yimhi : constant Standard_Floating_Vectors.Vector(1..dim)
          := (1..dim => 0.0);
    yrelo : constant Standard_Floating_Vectors.Vector(1..dim)
          := (1..dim => 0.0);
    yimlo : constant Standard_Floating_Vectors.Vector(1..dim)
          := (1..dim => 0.0);
    zrehi : constant Standard_Floating_Vectors.Vector(1..dim)
          := (1..dim => 0.0);
    zimhi : constant Standard_Floating_Vectors.Vector(1..dim)
          := (1..dim => 0.0);
    zrelo : constant Standard_Floating_Vectors.Vector(1..dim)
          := (1..dim => 0.0);
    zimlo : constant Standard_Floating_Vectors.Vector(1..dim)
          := (1..dim => 0.0);
    xrh : constant Standard_Floating_Vectors.Link_to_Vector
        := new Standard_Floating_Vectors.Vector'(xrehi);
    xih : constant Standard_Floating_Vectors.Link_to_Vector
        := new Standard_Floating_Vectors.Vector'(ximhi);
    xrl : constant Standard_Floating_Vectors.Link_to_Vector
        := new Standard_Floating_Vectors.Vector'(xrelo);
    xil : constant Standard_Floating_Vectors.Link_to_Vector
        := new Standard_Floating_Vectors.Vector'(ximlo);
    yrh : constant Standard_Floating_Vectors.Link_to_Vector
        := new Standard_Floating_Vectors.Vector'(yrehi);
    yih : constant Standard_Floating_Vectors.Link_to_Vector
        := new Standard_Floating_Vectors.Vector'(yimhi);
    yrl : constant Standard_Floating_Vectors.Link_to_Vector
        := new Standard_Floating_Vectors.Vector'(yrelo);
    yil : constant Standard_Floating_Vectors.Link_to_Vector
        := new Standard_Floating_Vectors.Vector'(yimlo);
    zrh : constant Standard_Floating_Vectors.Link_to_Vector
        := new Standard_Floating_Vectors.Vector'(yrehi);
    zih : constant Standard_Floating_Vectors.Link_to_Vector
        := new Standard_Floating_Vectors.Vector'(yimhi);
    zrl : constant Standard_Floating_Vectors.Link_to_Vector
        := new Standard_Floating_Vectors.Vector'(yrelo);
    zil : constant Standard_Floating_Vectors.Link_to_Vector
        := new Standard_Floating_Vectors.Vector'(yimlo);

    use DoblDobl_Complex_Vectors;

  begin
    tstart(timer);
    for k in 1..frq loop
      z1 := x + y;
    end loop;
    tstop(timer);
    new_line;
    print_times(standard_output,timer,"complex addition");
    Split(x,xrh.all,xih.all,xrl.all,xil.all);
    Split(y,yrh.all,yih.all,yrl.all,yil.all);
    tstart(timer);
    for k in 1..frq loop
      Add(zrh,zih,zrl,zil,xrh,xih,xrl,xil,yrh,yih,yrl,yil);
    end loop;
    tstop(timer);
    new_line;
    print_times(standard_output,timer,"adding splitted vectors");
    tstart(timer);
    for k in 1..frq loop
      Split(x,xrh.all,xih.all,xrl.all,xil.all);
      Split(y,yrh.all,yih.all,yrl.all,yil.all);
      Add(zrh,zih,zrl,zil,xrh,xih,xrl,xil,yrh,yih,yrl,yil);
      Merge(z2,zrh.all,zih.all,zrl.all,zil.all);
    end loop;
    tstop(timer);
    new_line;
    print_times(standard_output,timer,"split, add, merge  vectors");
  end Timing_Add;

  procedure Timing_Inner_Product ( dim,frq : in integer32 ) is

  -- DESCRIPTION :
  --   Generates random vectors and times their sum,
  --   for a frequency equal to frq.

    timer : Timing_Widget;
    x : constant DoblDobl_Complex_Vectors.Vector(1..dim)
      := DoblDobl_Random_Vectors.Random_Vector(1,dim);
    y : constant DoblDobl_Complex_Vectors.Vector(1..dim)
      := DoblDobl_Random_Vectors.Random_Vector(1,dim);
    z1,z2 : Complex_Number;
    xrehi,ximhi,xrelo,ximlo : Standard_Floating_Vectors.Vector(1..dim);
    yrehi,yimhi,yrelo,yimlo : Standard_Floating_Vectors.Vector(1..dim);
    xrh,xih,xrl,xil : Standard_Floating_Vectors.Link_to_Vector;
    yrh,yih,yrl,yil : Standard_Floating_Vectors.Link_to_Vector;
    zrehi,zimhi,zrelo,zimlo : double_float;

    use DoblDobl_Complex_Vectors;

  begin
    tstart(timer);
    for k in 1..frq loop
      z1 := Inner_Product(x,y);
    end loop;
    tstop(timer);
    new_line;
    print_times(standard_output,timer,"complex inner product");
    Split(x,xrehi,ximhi,xrelo,ximlo);
    Split(y,yrehi,yimhi,yrelo,yimlo);
    xrh := new Standard_Floating_Vectors.Vector'(xrehi);
    xih := new Standard_Floating_Vectors.Vector'(ximhi);
    xrl := new Standard_Floating_Vectors.Vector'(xrelo);
    xil := new Standard_Floating_Vectors.Vector'(ximlo);
    yrh := new Standard_Floating_Vectors.Vector'(yrehi);
    yih := new Standard_Floating_Vectors.Vector'(yimhi);
    yrl := new Standard_Floating_Vectors.Vector'(yrelo);
    yil := new Standard_Floating_Vectors.Vector'(yimlo);
    tstart(timer);
    for k in 1..frq loop
      Inner_Product(zrehi,zimhi,zrelo,zimlo,xrh,xih,xrl,xil,yrh,yih,yrl,yil);
    end loop;
    tstop(timer);
    new_line;
    print_times(standard_output,timer,"inner product on splitted vectors");
    tstart(timer);
    for k in 1..frq loop
      Split(x,xrehi,ximhi,xrelo,ximlo);
      Split(y,yrehi,yimhi,yrelo,yimlo);
      Inner_Product(zrehi,zimhi,zrelo,zimlo,xrh,xih,xrl,xil,yrh,yih,yrl,yil);
      Merge(z2,zrehi,zimhi,zrelo,zimlo);
    end loop;
    tstop(timer);
    new_line;
    print_times(standard_output,timer,"split, *, merge  vectors");
  end Timing_Inner_Product;

  procedure Main is

  -- DESCRIPTION :
  --   Prompts the user for the dimensions of the vectors.

    dim,frq : integer32 := 0;
    tst,ans : character;

  begin
    new_line;
    put_line("MENU to test vector operations :");
    put_line("  1. test addition");
    put_line("  2. test inner product");
    put("Type 1 or 2 to select a test : "); Ask_Alternative(tst,"12");
    new_line;
    put("Give the dimension : "); get(dim);
    new_line;
    put("Interactive test ? (y/n) "); Ask_Yes_or_No(ans);
    case tst is
      when '1' =>
        if ans = 'y' then
          Test_Add(dim);
        else
          new_line;
          put("Give the frequency : "); get(frq);
          Timing_Add(dim,frq);
        end if;
      when '2' =>
        if ans = 'y' then
          Test_Inner_Product(dim);
        else
          new_line;
          put("Give the frequency : "); get(frq);
          Timing_Inner_Product(dim,frq);
        end if;
      when others => null;
    end case;
  end Main;

begin
  Main;
end ts_perfddvc;
