with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with Timing_Package;                     use Timing_Package;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Double_Double_Basics;
with Quad_Double_Numbers;                use Quad_Double_Numbers;
with Quad_Double_Renormalizations;
with QuadDobl_Complex_Numbers;           use QuadDobl_Complex_Numbers;
with Standard_Floating_Vectors;
with QuadDobl_Complex_Vectors;
with QuadDobl_Complex_Vectors_io;        use QuadDobl_Complex_Vectors_io;
with QuadDobl_Random_Vectors;

procedure ts_perfqdvc is

-- DESCRIPTION :
--   Development of better performing computations with vectors of
--   complex numbers in quad double precision.

  procedure Split ( x : in Complex_Number;
                    rehihi,imhihi,relohi,imlohi : out double_float;
                    rehilo,imhilo,relolo,imlolo : out double_float ) is

  -- DESCRIPTION :
  --   Splits the complex number in quad double precision into 8 parts.

  -- ON ENTRY :
  --   x        a quad double complex number.

  -- ON RETURN :
  --   rehihi   highest double of the real part of x,
  --            or high_part(high_part(real_part(x)));
  --   imhihi   highest double of the imaginary part of x;
  --            or high_part(high_part(imag_part(x)));
  --   relohi   second highest double of the real part of x;
  --            or low_part(high_part(real_part(x)));
  --   imlohi   second highest double of the imaginary part of x;
  --            or low_part(high_part(imag_part(x)));
  --   rehilo   second lowest double of the real part of x;
  --            or high_part(low_part(real_part(x)));
  --   imhilo   second lowest double of the imaginary part of x;
  --            or high_part(low_part(imag_part(x)));
  --   relolo   lowest double of the real part of x;
  --            or low_part(low_part(imag_part(x)));
  --   imlolo   lowest double of the imaginary part of x;
  --            or low_part(low_part(real_part(x))).

    nbr : quad_double;

  begin
    nbr := REAL_PART(x);
    rehihi := hihi_part(nbr); relohi := lohi_part(nbr);
    rehilo := hilo_part(nbr); relolo := lolo_part(nbr);
    nbr := IMAG_PART(x);
    imhihi := hihi_part(nbr); imlohi := lohi_part(nbr);
    imhilo := hilo_part(nbr); imlolo := lolo_part(nbr);
  end Split;

  procedure Merge ( x : out Complex_Number;
                    rehihi,imhihi,relohi,imlohi : in double_float;
                    rehilo,imhilo,relolo,imlolo : in double_float ) is

  -- DESCRIPTION :
  --   Merges the 8 doubles into a quad double precision complex number.

  -- ON ENTRY :
  --   rehihi   highest double of the real part for x;
  --   imhihi   highest double of the imaginary part for x;
  --   relohi   second highest double of the real part for x;
  --   imlohi   second highest double of the imaginary part for x;
  --   rehilo   second lowest double of the real part for x;
  --   imhilo   second lowest double of the imaginary part for x;
  --   relolo   lowest double of the real part for x;
  --   imlolo   lowest double of the imaginary part for x.

  -- ON RETURN :
  --   x        a quad double complex number, with
  --            rehihi = high_part(high_part(real_part(x))),
  --            imhihi = high_part(high_part(imag_part(x))),
  --            relohi = low_part(high_part(real_part(x))),
  --            imlohi = low_part(high_part(imag_part(x))),
  --            rehilo = high_part(low_part(real_part(x))),
  --            imhilo = high_part(low_part(imag_part(x))),
  --            relolo = low_part(low_part(real_part(x))),
  --            imlolo = low_part(low_part(imag_part(x))).

    real_part : constant quad_double
              := Quad_Double_Numbers.create(rehihi,relohi,rehilo,relolo);
    imag_part : constant quad_double
              := Quad_Double_Numbers.create(imhihi,imlohi,imhilo,imlolo);

  begin
    x := QuadDobl_Complex_Numbers.Create(real_part,imag_part);
  end Merge;

  procedure Split ( x : in QuadDobl_Complex_Vectors.Vector;
                    xrhh,xihh : out Standard_Floating_Vectors.Vector;
                    xrlh,xilh : out Standard_Floating_Vectors.Vector;
                    xrhl,xihl : out Standard_Floating_Vectors.Vector;
                    xrll,xill : out Standard_Floating_Vectors.Vector ) is

  -- DESCRIPTION :
  --   Splits the vector x of quad double precision complex numbers
  --   into 8 vectors with real and imaginary parts,
  --   and for each part into four doubles.

  -- REQUIRED :
  --   x'range = xrhh'range = xihh'range = xrlh'range = xilh'range
  --           = xrhl'range = xihl'range = xrll'range = xill'range.

  -- ON ENTRY :
  --   x        a vector of quad double precision complex numbers.

  -- ON RETURN :
  --   xrhh     highest doubles of the real parts of x;
  --   xihh     highest doubles of the imaginary parts of x;
  --   xrlh     second highest doubles of the real parts of x;
  --   xilh     second highest doubles of the imaginary parts of x;
  --   xrhl     second lowest doubles of the real parts of x;
  --   xihl     second lowest doubles of the imaginary parts of x;
  --   xrll     lowest doubles of the real parts of x;
  --   xill     lowest doubles of the imaginary parts of x.

  begin
    for k in x'range loop
      Split(x(k),xrhh(k),xihh(k),xrlh(k),xilh(k),
                 xrhl(k),xihl(k),xrll(k),xill(k));
    end loop;
  end Split;

  procedure Merge ( x : out QuadDobl_Complex_Vectors.Vector;
                    xrhh,xihh : in Standard_Floating_Vectors.Vector;
                    xrlh,xilh : in Standard_Floating_Vectors.Vector;
                    xrhl,xihl : in Standard_Floating_Vectors.Vector;
                    xrll,xill : in Standard_Floating_Vectors.Vector ) is

  -- DESCRIPTION :
  --   Merges the 8 vectors of doubles into one vector of
  --   quad double precision complex numbers.

  -- REQUIRED :
  --   x'range = xrhh'range = xihh'range = xrlh'range = xilh'range
  --           = xrhl'range = xihl'range = xrll'range = xill'range.

  -- ON ENTRY :
  --   xrhh     highest doubles of the real parts for x;
  --   xihh     highest doubles of the imaginary parts for x;
  --   xrlh     second highest doubles of the real parts for x;
  --   xilh     second highest doubles of the imaginary parts for x;
  --   xrhl     second lowest doubles of the real parts for x;
  --   xihl     second lowest doubles of the imaginary parts for x;
  --   xrll     lowest doubles of the real parts for x;
  --   xill     lowest doubles of the imaginary parts for x.

  -- ON RETURN :
  --   x        a vector of quad double precision complex numbers, with
  --            xrhh = high_part(high_part(real_part(x))),
  --            xihh = high_part(high_part(imag_part(x))),
  --            xrlh = low_part(high_part(real_part(x))),
  --            xilh = low_part(high_part(imag_part(x))),
  --            xrhl = high_part(low_part(real_part(x))),
  --            xihl = high_part(low_part(imag_part(x))),
  --            xrll = low_part(low_part(real_part(x))),
  --            xill = low_part(low_part(imag_part(x))).

  begin
    for k in x'range loop
      Merge(x(k),xrhh(k),xihh(k),xrlh(k),xilh(k),
                 xrhl(k),xihl(k),xrll(k),xill(k));
    end loop;
  end Merge;

  procedure Add ( x,y : in Standard_Floating_Vectors.Vector;
                  z : out Standard_Floating_Vectors.Vector ) is

  -- DESCRIPTION :
  --   Adds two quad double numbers, x and y,
  --   and assigns the sum to the quad double number z.
  --   The vector representation for quad doubles is assumed.

  -- REQUIRED : x'range = y'range = z'range = 0..3.

  -- ON ENTRY :
  --   x        a quad double number, with
  --            x(0) the highest double of x,
  --            x(1) the second highest double of x,
  --            x(2) the second lowest double of x,
  --            x(3) the lowest double of x;
  --   y        a quad double number, with
  --            y(0) the highest double of y,
  --            y(1) the second highest double of y,
  --            y(2) the second lowest double of y,
  --            y(3) the lowest double of y.

  -- ON RETURN :
  --   z        the sum of x + y, with
  --            z(0) the highest double of x + y,
  --            z(1) the second highest double of x + y,
  --            z(2) the second lowest double of x + y,
  --            z(3) the lowest double of x + y.

    i,j,k : integer32;
    s,t : double_float;
    u,v : double_float; -- double-length accumulator

  begin
    z(0) := 0.0; z(1) := 0.0; z(2) := 0.0; z(3) := 0.0;
    i := 0; j := 0; k := 0;
    if abs(x(i)) > abs(y(j))
     then u := x(i); i := i+1;
     else u := y(j); j := j+1;
    end if;
    if abs(x(i)) > abs(y(j))
     then v := x(i); i := i+1;  
     else v := y(j); j := j+1;
    end if;
    Double_Double_Basics.quick_two_sum(u,v,u,v);
    while k < 4 loop
      if (i >= 4 and j >= 4) then
        z(k) := u;
        if k < 3
         then k := k+1; z(k) := v;
        end if;
        exit;
      end if;
      if i >= 4 then
        t := y(j); j := j+1;
      elsif j >= 4  then
        t := x(i); i := i+1;
      elsif abs(x(i)) > abs(y(j)) then
        t := x(i); i := i+1;
      else 
        t := y(j); j := j+1;
      end if;
      s := 0.0;
      Quad_Double_Renormalizations.quick_three_accum(u,v,s,t);
      if s /= 0.0
       then z(k) := s; k := k+1;
      end if;
    end loop;
    for k in i..3 loop                    -- add the rest
      z(3) := z(3) + x(k);
    end loop;
    for k in j..3 loop
      z(3) := z(3) + y(k);
    end loop;
    Quad_Double_Renormalizations.renorm4(z(0),z(1),z(2),z(3));
  end Add;

  procedure Add ( zrhh : in Standard_Floating_Vectors.Link_to_Vector;
                  zihh : in Standard_Floating_Vectors.Link_to_Vector;
                  zrlh : in Standard_Floating_Vectors.Link_to_Vector;
                  zilh : in Standard_Floating_Vectors.Link_to_Vector;
                  zrhl : in Standard_Floating_Vectors.Link_to_Vector;
                  zihl : in Standard_Floating_Vectors.Link_to_Vector;
                  zrll : in Standard_Floating_Vectors.Link_to_Vector;
                  zill : in Standard_Floating_Vectors.Link_to_Vector;
                  xrhh : in Standard_Floating_Vectors.Link_to_Vector;
                  xihh : in Standard_Floating_Vectors.Link_to_Vector;
                  xrlh : in Standard_Floating_Vectors.Link_to_Vector;
                  xilh : in Standard_Floating_Vectors.Link_to_Vector;
                  xrhl : in Standard_Floating_Vectors.Link_to_Vector;
                  xihl : in Standard_Floating_Vectors.Link_to_Vector;
                  xrll : in Standard_Floating_Vectors.Link_to_Vector;
                  xill : in Standard_Floating_Vectors.Link_to_Vector;
                  yrhh : in Standard_Floating_Vectors.Link_to_Vector;
                  yihh : in Standard_Floating_Vectors.Link_to_Vector;
                  yrlh : in Standard_Floating_Vectors.Link_to_Vector;
                  yilh : in Standard_Floating_Vectors.Link_to_Vector;
                  yrhl : in Standard_Floating_Vectors.Link_to_Vector;
                  yihl : in Standard_Floating_Vectors.Link_to_Vector;
                  yrll : in Standard_Floating_Vectors.Link_to_Vector;
                  yill : in Standard_Floating_Vectors.Link_to_Vector;
                  x,y,z : in out Standard_Floating_Vectors.Vector ) is

  -- DESCRIPTION :
  --   Adds two vectors of quad double complex numbers,
  --   in the 8-vector representation of doubles.

  -- REQUIRED : all vectors have the same range.

  -- ON ENTRY :
  --   xrhh     highest double of the real part of x;
  --   xihh     highest double of the imaginary part of x;
  --   xrlh     second highest double of the real part of x;
  --   xilh     second highest double of the imaginary part of x;
  --   xrhl     second lowest double of the real part of x;
  --   xihl     second lowest double of the imaginary part of x;
  --   xrll     lowest double of the real part of x;
  --   xill     lowest double of the imaginary part of x;
  --   yrhh     highest double of the real part of y;
  --   yihh     highest double of the imaginary part of y;
  --   yrlh     second highest double of the real part of y;
  --   yilh     second highest double of the imaginary part of y;
  --   yrhl     second lowest double of the real part of y;
  --   yihl     second lowest double of the imaginary part of y;
  --   yrll     lowest double of the real part of y;
  --   yill     lowest double of the imaginary part of y;
  --   x        a 4-vector work space of range 0..3;
  --   y        a 4-vector work space of range 0..3;
  --   z        a 4-vector work space of range 0..3.

  -- ON RETURN :
  --   zrhh     highest double of the real part of the sum;
  --   zihh     highest double of the imaginary part of the sum;
  --   zrlh     second highest double of the real part of the sum;
  --   zilh     second highest double of the imaginary part of the sum;
  --   zrhl     second lowest double of the real part of the sum;
  --   zihl     second lowest double of the imaginary part of the sum;
  --   zrll     lowest double of the real part of the sum;
  --   zill     lowest double of the imaginary part of the sum.

  begin
    for k in zrhh'range loop
      x(0) := xrhh(k); x(1) := xrlh(k); x(2) := xrhl(k); x(3) := xrll(k);
      y(0) := yrhh(k); y(1) := yrlh(k); y(2) := yrhl(k); y(3) := yrll(k);
      Add(x,y,z);
      zrhh(k) := z(0); zrlh(k) := z(1); zrhl(k) := z(2); zrll(k) := z(3); 
      x(0) := xihh(k); x(1) := xilh(k); x(2) := xihl(k); x(3) := xill(k);
      y(0) := yihh(k); y(1) := yilh(k); y(2) := yihl(k); y(3) := yill(k);
      Add(x,y,z);
      zihh(k) := z(0); zilh(k) := z(1); zihl(k) := z(2); zill(k) := z(3); 
    end loop;
  end Add;

  procedure Test_Add ( dim : in integer32 ) is

  -- DESCRIPTION :
  --   Generates two random vectors and tests their sum.

    x : constant QuadDobl_Complex_Vectors.Vector(1..dim)
      := QuadDobl_Random_Vectors.Random_Vector(1,dim);
    y : constant QuadDobl_Complex_Vectors.Vector(1..dim)
      := QuadDobl_Random_Vectors.Random_Vector(1,dim);
    z1,z2 : QuadDobl_Complex_Vectors.Vector(1..dim);
    xrhh,xihh,xrlh,xilh : Standard_Floating_Vectors.Vector(1..dim);
    xrhl,xihl,xrll,xill : Standard_Floating_Vectors.Vector(1..dim);
    yrhh,yihh,yrlh,yilh : Standard_Floating_Vectors.Vector(1..dim);
    yrhl,yihl,yrll,yill : Standard_Floating_Vectors.Vector(1..dim);
    zrhh,zihh,zrlh,zilh : constant Standard_Floating_Vectors.Vector(1..dim)
                        := (1..dim => 0.0);
    zrhl,zihl,zrll,zill : constant Standard_Floating_Vectors.Vector(1..dim)
                        := (1..dim => 0.0);
    urhh,uihh,urlh,uilh : Standard_Floating_Vectors.Link_to_Vector;
    urhl,uihl,urll,uill : Standard_Floating_Vectors.Link_to_Vector;
    vrhh,vihh,vrlh,vilh : Standard_Floating_Vectors.Link_to_Vector;
    vrhl,vihl,vrll,vill : Standard_Floating_Vectors.Link_to_Vector;
    wrhh,wihh,wrlh,wilh : Standard_Floating_Vectors.Link_to_Vector;
    wrhl,wihl,wrll,will : Standard_Floating_Vectors.Link_to_Vector;
    xwrk,ywrk,zwrk : Standard_Floating_Vectors.Vector(0..3);

    use QuadDobl_Complex_Vectors; -- for the + operator

  begin
    z1 := x + y;
    put_line("The sum of two random vectors :"); put_line(z1);
    Split(x,xrhh,xihh,xrlh,xilh,xrhl,xihl,xrll,xill);
    Split(y,yrhh,yihh,yrlh,yilh,yrhl,yihl,yrll,yill);
    urhh := new Standard_Floating_Vectors.Vector'(xrhh);
    uihh := new Standard_Floating_Vectors.Vector'(xihh);
    urlh := new Standard_Floating_Vectors.Vector'(xrlh);
    uilh := new Standard_Floating_Vectors.Vector'(xilh);
    urhl := new Standard_Floating_Vectors.Vector'(xrhl);
    uihl := new Standard_Floating_Vectors.Vector'(xihl);
    urll := new Standard_Floating_Vectors.Vector'(xrll);
    uill := new Standard_Floating_Vectors.Vector'(xill);
    vrhh := new Standard_Floating_Vectors.Vector'(yrhh);
    vihh := new Standard_Floating_Vectors.Vector'(yihh);
    vrlh := new Standard_Floating_Vectors.Vector'(yrlh);
    vilh := new Standard_Floating_Vectors.Vector'(yilh);
    vrhl := new Standard_Floating_Vectors.Vector'(yrhl);
    vihl := new Standard_Floating_Vectors.Vector'(yihl);
    vrll := new Standard_Floating_Vectors.Vector'(yrll);
    vill := new Standard_Floating_Vectors.Vector'(yill);
    wrhh := new Standard_Floating_Vectors.Vector'(zrhh);
    wihh := new Standard_Floating_Vectors.Vector'(zihh);
    wrlh := new Standard_Floating_Vectors.Vector'(zrlh);
    wilh := new Standard_Floating_Vectors.Vector'(zilh);
    wrhl := new Standard_Floating_Vectors.Vector'(zrhl);
    wihl := new Standard_Floating_Vectors.Vector'(zihl);
    wrll := new Standard_Floating_Vectors.Vector'(zrll);
    will := new Standard_Floating_Vectors.Vector'(zill);
    Add(wrhh,wihh,wrlh,wilh,wrhl,wihl,wrll,will,
        urhh,uihh,urlh,uilh,urhl,uihl,urll,uill,
        vrhh,vihh,vrlh,vilh,vrhl,vihl,vrll,vill,xwrk,ywrk,zwrk);
    Merge(z2,wrhh.all,wihh.all,wrlh.all,wilh.all,
             wrhl.all,wihl.all,wrll.all,will.all);
    put_line("The recomputed sum :"); put_line(z2);
  end Test_Add;

  procedure Timing_Add ( dim,frq : in integer32 ) is

  -- DESCRIPTION :
  --   Generates two random vectors and times the sum,
  --   for a frequency equal to the value of frq.

    timer : Timing_Widget;
    x : constant QuadDobl_Complex_Vectors.Vector(1..dim)
      := QuadDobl_Random_Vectors.Random_Vector(1,dim);
    y : constant QuadDobl_Complex_Vectors.Vector(1..dim)
      := QuadDobl_Random_Vectors.Random_Vector(1,dim);
    z1,z2 : QuadDobl_Complex_Vectors.Vector(1..dim);
    xrhh,xihh,xrlh,xilh : Standard_Floating_Vectors.Vector(1..dim);
    xrhl,xihl,xrll,xill : Standard_Floating_Vectors.Vector(1..dim);
    yrhh,yihh,yrlh,yilh : Standard_Floating_Vectors.Vector(1..dim);
    yrhl,yihl,yrll,yill : Standard_Floating_Vectors.Vector(1..dim);
    zrhh,zihh,zrlh,zilh : constant Standard_Floating_Vectors.Vector(1..dim)
                        := (1..dim => 0.0);
    zrhl,zihl,zrll,zill : constant Standard_Floating_Vectors.Vector(1..dim)
                        := (1..dim => 0.0);
    urhh,uihh,urlh,uilh : Standard_Floating_Vectors.Link_to_Vector;
    urhl,uihl,urll,uill : Standard_Floating_Vectors.Link_to_Vector;
    vrhh,vihh,vrlh,vilh : Standard_Floating_Vectors.Link_to_Vector;
    vrhl,vihl,vrll,vill : Standard_Floating_Vectors.Link_to_Vector;
    wrhh,wihh,wrlh,wilh : Standard_Floating_Vectors.Link_to_Vector;
    wrhl,wihl,wrll,will : Standard_Floating_Vectors.Link_to_Vector;
    xwrk,ywrk,zwrk : Standard_Floating_Vectors.Vector(0..3);

    use QuadDobl_Complex_Vectors; -- for the + operator

  begin
    tstart(timer);
    for k in 1..frq loop
      z1 := x + y;
    end loop;
    tstop(timer);
    new_line;
    print_times(standard_output,timer,"complex addition");
    Split(x,xrhh,xihh,xrlh,xilh,xrhl,xihl,xrll,xill);
    Split(y,yrhh,yihh,yrlh,yilh,yrhl,yihl,yrll,yill);
    urhh := new Standard_Floating_Vectors.Vector'(xrhh);
    uihh := new Standard_Floating_Vectors.Vector'(xihh);
    urlh := new Standard_Floating_Vectors.Vector'(xrlh);
    uilh := new Standard_Floating_Vectors.Vector'(xilh);
    urhl := new Standard_Floating_Vectors.Vector'(xrhl);
    uihl := new Standard_Floating_Vectors.Vector'(xihl);
    urll := new Standard_Floating_Vectors.Vector'(xrll);
    uill := new Standard_Floating_Vectors.Vector'(xill);
    vrhh := new Standard_Floating_Vectors.Vector'(yrhh);
    vihh := new Standard_Floating_Vectors.Vector'(yihh);
    vrlh := new Standard_Floating_Vectors.Vector'(yrlh);
    vilh := new Standard_Floating_Vectors.Vector'(yilh);
    vrhl := new Standard_Floating_Vectors.Vector'(yrhl);
    vihl := new Standard_Floating_Vectors.Vector'(yihl);
    vrll := new Standard_Floating_Vectors.Vector'(yrll);
    vill := new Standard_Floating_Vectors.Vector'(yill);
    wrhh := new Standard_Floating_Vectors.Vector'(zrhh);
    wihh := new Standard_Floating_Vectors.Vector'(zihh);
    wrlh := new Standard_Floating_Vectors.Vector'(zrlh);
    wilh := new Standard_Floating_Vectors.Vector'(zilh);
    wrhl := new Standard_Floating_Vectors.Vector'(zrhl);
    wihl := new Standard_Floating_Vectors.Vector'(zihl);
    wrll := new Standard_Floating_Vectors.Vector'(zrll);
    will := new Standard_Floating_Vectors.Vector'(zill);
    tstart(timer);
    for k in 1..frq loop
      Add(wrhh,wihh,wrlh,wilh,wrhl,wihl,wrll,will,
          urhh,uihh,urlh,uilh,urhl,uihl,urll,uill,
          vrhh,vihh,vrlh,vilh,vrhl,vihl,vrll,vill,xwrk,ywrk,zwrk);
    end loop;
    tstop(timer);
    new_line;
    print_times(standard_output,timer,"refitted addition");
    tstart(timer);
    for k in 1..frq loop
      Split(x,xrhh,xihh,xrlh,xilh,xrhl,xihl,xrll,xill);
      Split(y,yrhh,yihh,yrlh,yilh,yrhl,yihl,yrll,yill);
      Add(wrhh,wihh,wrlh,wilh,wrhl,wihl,wrll,will,
          urhh,uihh,urlh,uilh,urhl,uihl,urll,uill,
          vrhh,vihh,vrlh,vilh,vrhl,vihl,vrll,vill,xwrk,ywrk,zwrk);
      Merge(z2,wrhh.all,wihh.all,wrlh.all,wilh.all,
               wrhl.all,wihl.all,wrll.all,will.all);
    end loop;
    tstop(timer);
    new_line;
    print_times(standard_output,timer,"split, add, merge");
  end Timing_Add;

  procedure Main is

  -- DESCRIPTION :
  --   Prompts the user for the dimensions of the vectors.

    dim,frq : integer32 := 0;
    ans : character;

  begin
    new_line;
    put("Give the dimension : "); get(dim);
    new_line;
    put("Interactive test ? (y/n) "); Ask_Yes_or_No(ans);
    if ans = 'y' then
      Test_Add(dim);
    else
      new_line;
      put("Give the frequency : "); get(frq);
      Timing_Add(dim,frq);
    end if;
  end Main;

begin
  Main;
end ts_perfqdvc;
