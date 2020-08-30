with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with Timing_Package;                     use Timing_Package;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with DoblDobl_Complex_Numbers_io;        use DoblDobl_Complex_Numbers_io;
with DoblDobl_Complex_Vectors_io;        use DoblDobl_Complex_Vectors_io;
with DoblDobl_Random_Vectors;
with DoblDobl_Vector_Splitters;          use DoblDobl_Vector_Splitters;

package body Test_DoblDobl_Performance is

  function Inner_Product
             ( x,y : DoblDobl_Complex_Vectors.Vector )
             return Complex_Number is

    res : Complex_Number := Create(integer(0));

  begin
    for k in x'range loop
      res := res + x(k)*y(k);
    end loop;
    return res;   
  end Inner_Product;

  procedure Inner_Product_Inlined
                ( zrehi,zimhi,zrelo,zimlo : out double_float;
                  xrehi : in Standard_Floating_Vectors.Link_to_Vector;
                  ximhi : in Standard_Floating_Vectors.Link_to_Vector;
                  xrelo : in Standard_Floating_Vectors.Link_to_Vector;
                  ximlo : in Standard_Floating_Vectors.Link_to_Vector;
                  yrehi : in Standard_Floating_Vectors.Link_to_Vector;
                  yimhi : in Standard_Floating_Vectors.Link_to_Vector;
                  yrelo : in Standard_Floating_Vectors.Link_to_Vector;
                  yimlo : in Standard_Floating_Vectors.Link_to_Vector) is

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
  end Inner_Product_Inlined;

  procedure Multiply
              ( first,second : in DoblDobl_Complex_Vectors.Link_to_Vector;
                product : in DoblDobl_Complex_Vectors.Link_to_Vector ) is

    deg : constant integer32 := first'last;

  begin
    product(0) := first(0)*second(0);
    for k in 1..deg loop
      product(k) := first(0)*second(k);
      for i in 1..k loop
        product(k) := product(k) + first(i)*second(k-i);
      end loop;
    end loop;
  end Multiply;

  procedure Test_Add ( dim : in integer32 ) is

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
        := new Standard_Floating_Vectors.Vector'(zrehi);
    zih : constant Standard_Floating_Vectors.Link_to_Vector
        := new Standard_Floating_Vectors.Vector'(zimhi);
    zrl : constant Standard_Floating_Vectors.Link_to_Vector
        := new Standard_Floating_Vectors.Vector'(zrelo);
    zil : constant Standard_Floating_Vectors.Link_to_Vector
        := new Standard_Floating_Vectors.Vector'(zimlo);

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

  procedure Test_Update ( dim : in integer32 ) is

    x : constant DoblDobl_Complex_Vectors.Vector(1..dim)
      := DoblDobl_Random_Vectors.Random_Vector(1,dim);
    z1 : DoblDobl_Complex_Vectors.Vector(1..dim)
       := DoblDobl_Random_Vectors.Random_Vector(1,dim);
    z2 : DoblDobl_Complex_Vectors.Vector(1..dim) := z1;
    xrehi : constant Standard_Floating_Vectors.Vector(1..dim)
          := (1..dim => 0.0);
    ximhi : constant Standard_Floating_Vectors.Vector(1..dim)
          := (1..dim => 0.0);
    xrelo : constant Standard_Floating_Vectors.Vector(1..dim)
          := (1..dim => 0.0);
    ximlo : constant Standard_Floating_Vectors.Vector(1..dim)
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
    zrh : constant Standard_Floating_Vectors.Link_to_Vector
        := new Standard_Floating_Vectors.Vector'(zrehi);
    zih : constant Standard_Floating_Vectors.Link_to_Vector
        := new Standard_Floating_Vectors.Vector'(zimhi);
    zrl : constant Standard_Floating_Vectors.Link_to_Vector
        := new Standard_Floating_Vectors.Vector'(zrelo);
    zil : constant Standard_Floating_Vectors.Link_to_Vector
        := new Standard_Floating_Vectors.Vector'(zimlo);

    use DoblDobl_Complex_Vectors;

  begin
    z1 := z1 + x;
    put_line("The update with two random vectors :"); put_line(z1);
    Split(z2,zrh.all,zih.all,zrl.all,zil.all);
    Split(x,xrh.all,xih.all,xrl.all,xil.all);
    Update(zrh,zih,zrl,zil,xrh,xih,xrl,xil);
    Merge(z2,zrh.all,zih.all,zrl.all,zil.all);
    put_line("The recomputed update :"); put_line(z2);
  end Test_Update;

  procedure Test_Inner_Product ( dim : in integer32 ) is

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

  procedure Test_Multiply ( dim : in integer32 ) is

    cx : constant DoblDobl_Complex_Vectors.Vector(0..dim)
       := DoblDobl_Random_Vectors.Random_Vector(0,dim);
    cy : constant DoblDobl_Complex_Vectors.Vector(0..dim)
       := DoblDobl_Random_Vectors.Random_Vector(0,dim);
    x : constant DoblDobl_Complex_Vectors.Link_to_Vector
      := new DoblDobl_Complex_Vectors.Vector'(cx);
    y : constant DoblDobl_Complex_Vectors.Link_to_Vector
      := new DoblDobl_Complex_Vectors.Vector'(cy);
    zero : constant Complex_Number := Create(integer(0));
    z1 : constant DoblDobl_Complex_Vectors.Link_to_Vector
       := new DoblDobl_Complex_Vectors.Vector'(0..dim => zero);
    z2 : DoblDobl_Complex_Vectors.Vector(0..dim);
    xrehi,ximhi,xrelo,ximlo : Standard_Floating_Vectors.Vector(0..dim);
    yrehi,yimhi,yrelo,yimlo : Standard_Floating_Vectors.Vector(0..dim);
    zrehi,zimhi,zrelo,zimlo : Standard_Floating_Vectors.Vector(0..dim);
    xrh,xih,xrl,xil : Standard_Floating_Vectors.Link_to_Vector;
    yrh,yih,yrl,yil : Standard_Floating_Vectors.Link_to_Vector;
    zrh,zih,zrl,zil : Standard_Floating_Vectors.Link_to_Vector;

  begin
    Multiply(x,y,z1);
    put_line("The convolution of two random vectors :"); put_line(z1);
    Split(cx,xrehi,ximhi,xrelo,ximlo);
    Split(cy,yrehi,yimhi,yrelo,yimlo);
    xrh := new Standard_Floating_Vectors.Vector'(xrehi);
    xih := new Standard_Floating_Vectors.Vector'(ximhi);
    xrl := new Standard_Floating_Vectors.Vector'(xrelo);
    xil := new Standard_Floating_Vectors.Vector'(ximlo);
    yrh := new Standard_Floating_Vectors.Vector'(yrehi);
    yih := new Standard_Floating_Vectors.Vector'(yimhi);
    yrl := new Standard_Floating_Vectors.Vector'(yrelo);
    yil := new Standard_Floating_Vectors.Vector'(yimlo);
    zrh := new Standard_Floating_Vectors.Vector'(zrehi);
    zih := new Standard_Floating_Vectors.Vector'(zimhi);
    zrl := new Standard_Floating_Vectors.Vector'(zrelo);
    zil := new Standard_Floating_Vectors.Vector'(zimlo);
    Multiply(xrh,xih,xrl,xil,yrh,yih,yrl,yil,zrh,zih,zrl,zil);
    Merge(z2,zrh.all,zih.all,zrl.all,zil.all);
    put_line("The recomputed convolution :"); put_line(z2);
  end Test_Multiply;

  procedure Timing_Add ( dim,frq : in integer32 ) is

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
        := new Standard_Floating_Vectors.Vector'(zrehi);
    zih : constant Standard_Floating_Vectors.Link_to_Vector
        := new Standard_Floating_Vectors.Vector'(zimhi);
    zrl : constant Standard_Floating_Vectors.Link_to_Vector
        := new Standard_Floating_Vectors.Vector'(zrelo);
    zil : constant Standard_Floating_Vectors.Link_to_Vector
        := new Standard_Floating_Vectors.Vector'(zimlo);

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

  procedure Timing_Update ( dim,frq : in integer32 ) is

    timer : Timing_Widget;
    z : constant DoblDobl_Complex_Vectors.Vector(1..dim)
      := DoblDobl_Random_Vectors.Random_Vector(1,dim);
    x : constant DoblDobl_Complex_Vectors.Vector(1..dim)
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
    zrh : constant Standard_Floating_Vectors.Link_to_Vector
        := new Standard_Floating_Vectors.Vector'(zrehi);
    zih : constant Standard_Floating_Vectors.Link_to_Vector
        := new Standard_Floating_Vectors.Vector'(zimhi);
    zrl : constant Standard_Floating_Vectors.Link_to_Vector
        := new Standard_Floating_Vectors.Vector'(zrelo);
    zil : constant Standard_Floating_Vectors.Link_to_Vector
        := new Standard_Floating_Vectors.Vector'(zimlo);

    use DoblDobl_Complex_Vectors;

  begin
    tstart(timer);
    for k in 1..frq loop
      z1 := z; z1 := z1 + x;
    end loop;
    tstop(timer);
    new_line;
    print_times(standard_output,timer,"complex update");
    Split(z,zrh.all,zih.all,zrl.all,zil.all);
    Split(x,xrh.all,xih.all,xrl.all,xil.all);
    tstart(timer);
    for k in 1..frq loop
      Update(zrh,zih,zrl,zil,xrh,xih,xrl,xil);
    end loop;
    tstop(timer);
    new_line;
    print_times(standard_output,timer,"updating splitted vectors");
    tstart(timer);
    for k in 1..frq loop
      Split(x,xrh.all,xih.all,xrl.all,xil.all);
      Split(z,zrh.all,zih.all,zrl.all,zil.all);
      Update(zrh,zih,zrl,zil,xrh,xih,xrl,xil);
      Merge(z2,zrh.all,zih.all,zrl.all,zil.all);
    end loop;
    tstop(timer);
    new_line;
    print_times(standard_output,timer,"split, add, merge  vectors");
  end Timing_Update;

  procedure Timing_Inner_Product ( dim,frq : in integer32 ) is

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
    print_times(standard_output,timer,"split, *, merge vectors");
  end Timing_Inner_Product;

  procedure Timing_Multiply ( dim,frq : in integer32 ) is

    timer : Timing_Widget;
    cx : constant DoblDobl_Complex_Vectors.Vector(0..dim)
       := DoblDobl_Random_Vectors.Random_Vector(0,dim);
    cy : constant DoblDobl_Complex_Vectors.Vector(0..dim)
       := DoblDobl_Random_Vectors.Random_Vector(0,dim);
    x : constant DoblDobl_Complex_Vectors.Link_to_Vector
      := new DoblDobl_Complex_Vectors.Vector'(cx);
    y : constant DoblDobl_Complex_Vectors.Link_to_Vector
      := new DoblDobl_Complex_Vectors.Vector'(cy);
    zero : constant Complex_Number := Create(integer(0));
    z1 : constant DoblDobl_Complex_Vectors.Link_to_Vector
       := new DoblDobl_Complex_Vectors.Vector'(0..dim => zero);
    z2 : DoblDobl_Complex_Vectors.Vector(0..dim);
    xrehi : constant Standard_Floating_Vectors.Vector(0..dim)
          := (0..dim => 0.0);
    ximhi : constant Standard_Floating_Vectors.Vector(0..dim)
          := (0..dim => 0.0);
    xrelo : constant Standard_Floating_Vectors.Vector(0..dim)
          := (0..dim => 0.0);
    ximlo : constant Standard_Floating_Vectors.Vector(0..dim)
          := (0..dim => 0.0);
    yrehi : constant Standard_Floating_Vectors.Vector(0..dim)
          := (0..dim => 0.0);
    yimhi : constant Standard_Floating_Vectors.Vector(0..dim)
          := (0..dim => 0.0);
    yrelo : constant Standard_Floating_Vectors.Vector(0..dim)
          := (0..dim => 0.0);
    yimlo : constant Standard_Floating_Vectors.Vector(0..dim)
          := (0..dim => 0.0);
    zrehi : constant Standard_Floating_Vectors.Vector(0..dim)
          := (0..dim => 0.0);
    zimhi : constant Standard_Floating_Vectors.Vector(0..dim)
          := (0..dim => 0.0);
    zrelo : constant Standard_Floating_Vectors.Vector(0..dim)
          := (0..dim => 0.0);
    zimlo : constant Standard_Floating_Vectors.Vector(0..dim)
          := (0..dim => 0.0);
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
        := new Standard_Floating_Vectors.Vector'(zrehi);
    zih : constant Standard_Floating_Vectors.Link_to_Vector
        := new Standard_Floating_Vectors.Vector'(zimhi);
    zrl : constant Standard_Floating_Vectors.Link_to_Vector
        := new Standard_Floating_Vectors.Vector'(zrelo);
    zil : constant Standard_Floating_Vectors.Link_to_Vector
        := new Standard_Floating_Vectors.Vector'(zimlo);

    use DoblDobl_Complex_Vectors;

  begin
    tstart(timer);
    for k in 1..frq loop
      Multiply(x,y,z1);
    end loop;
    tstop(timer);
    new_line;
    print_times(standard_output,timer,"complex convolution");
    Split(cx,xrh.all,xih.all,xrl.all,xil.all);
    Split(cy,yrh.all,yih.all,yrl.all,yil.all);
    tstart(timer);
    for k in 1..frq loop
      Multiply(xrh,xih,xrl,xil,yrh,yih,yrl,yil,zrh,zih,zrl,zil);
    end loop;
    tstop(timer);
    new_line;
    print_times(standard_output,timer,"convoluting splitted vectors");
    tstart(timer);
    for k in 1..frq loop
      Split(cx,xrh.all,xih.all,xrl.all,xil.all);
      Split(cy,yrh.all,yih.all,yrl.all,yil.all);
      Multiply(xrh,xih,xrl,xil,yrh,yih,yrl,yil,zrh,zih,zrl,zil);
      Merge(z2,zrh.all,zih.all,zrl.all,zil.all);
    end loop;
    tstop(timer);
    new_line;
    print_times(standard_output,timer,"split, *, merge vectors");
  end Timing_Multiply;

  procedure Main is

    dim,frq : integer32 := 0;
    tst,ans : character;

  begin
    new_line;
    put_line("MENU to test vector operations :");
    put_line("  1. test addition");
    put_line("  2. test update");
    put_line("  3. test inner product");
    put_line("  4. test convolution");
    put("Type 1, 2, 3, or 4 to select a test : ");
    Ask_Alternative(tst,"1234");
    new_line;
    put("Give the dimension : "); get(dim);
    new_line;
    put("Interactive test ? (y/n) "); Ask_Yes_or_No(ans);
    if ans /= 'y' then
      new_line;
      put("Give the frequency : "); get(frq);
    end if;
    if frq = 0 then
      case tst is
        when '1' => Test_Add(dim);
        when '2' => Test_Update(dim);
        when '3' => Test_Inner_Product(dim);
        when '4' => Test_Multiply(dim);
        when others => null;
      end case;
    else
      case tst is
        when '1' => Timing_Add(dim,frq);
        when '2' => Timing_Update(dim,frq);
        when '3' => Timing_Inner_Product(dim,frq);
        when '4' => Timing_Multiply(dim,frq);
        when others => null;
      end case;
    end if;
  end Main;

end Test_DoblDobl_Performance;
