with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with Timing_Package;                     use Timing_Package;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with QuadDobl_Complex_Numbers_io;        use QuadDobl_Complex_Numbers_io;
with Standard_Floating_Vectors;
with QuadDobl_Complex_Vectors_io;        use QuadDobl_Complex_Vectors_io;
with QuadDobl_Random_Vectors;
with QuadDobl_Vector_Splitters;          use QuadDobl_Vector_Splitters;

package body Test_QuadDobl_Performance is

  function Inner_Product
             ( x,y : QuadDobl_Complex_Vectors.Vector )
             return Complex_Number is

    res : Complex_Number := Create(integer(0));

  begin
    for k in x'range loop
      res := res + x(k)*y(k);
    end loop;
    return res;   
  end Inner_Product;

  procedure Multiply
              ( first,second : in QuadDobl_Complex_Vectors.Link_to_Vector;
                product : in QuadDobl_Complex_Vectors.Link_to_Vector ) is

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
    xwrk,ywrk,zwrk : constant Standard_Floating_Vectors.Vector(0..3)
                   := (0..3 => 0.0);
    u : constant Standard_Floating_Vectors.Link_to_Vector
      := new Standard_Floating_Vectors.Vector'(xwrk);
    v : constant Standard_Floating_Vectors.Link_to_Vector
      := new Standard_Floating_Vectors.Vector'(ywrk);
    w : constant Standard_Floating_Vectors.Link_to_Vector
      := new Standard_Floating_Vectors.Vector'(zwrk);

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
        vrhh,vihh,vrlh,vilh,vrhl,vihl,vrll,vill,u,v,w);
    Merge(z2,wrhh.all,wihh.all,wrlh.all,wilh.all,
             wrhl.all,wihl.all,wrll.all,will.all);
    put_line("The recomputed sum :"); put_line(z2);
  end Test_Add;

  procedure Test_Two_Add ( dim : in integer32 ) is

    x : constant QuadDobl_Complex_Vectors.Vector(1..dim)
      := QuadDobl_Random_Vectors.Random_Vector(1,dim);
    y : constant QuadDobl_Complex_Vectors.Vector(1..dim)
      := QuadDobl_Random_Vectors.Random_Vector(1,dim);
    z1,z2 : QuadDobl_Complex_Vectors.Vector(1..dim);
    xr,xi : Standard_Floating_Vectors.Vector(1..4*dim);
    yr,yi : Standard_Floating_Vectors.Vector(1..4*dim);
    zr,zi : constant Standard_Floating_Vectors.Vector(1..4*dim)
          := (1..4*dim => 0.0);
    ur,ui : Standard_Floating_Vectors.Link_to_Vector;
    vr,vi : Standard_Floating_Vectors.Link_to_Vector;
    wr,wi : Standard_Floating_Vectors.Link_to_Vector;

    use QuadDobl_Complex_Vectors; -- for the + operator

  begin
    z1 := x + y;
    put_line("The sum of two random vectors :"); put_line(z1);
    Two_Split(x,xr,xi); Two_Split(y,yr,yi);
    ur := new Standard_Floating_Vectors.Vector'(xr);
    ui := new Standard_Floating_Vectors.Vector'(xi);
    vr := new Standard_Floating_Vectors.Vector'(yr);
    vi := new Standard_Floating_Vectors.Vector'(yi);
    wr := new Standard_Floating_Vectors.Vector'(zr);
    wi := new Standard_Floating_Vectors.Vector'(zi);
    Two_Add(wr,wi,ur,ui,vr,vi);
    Two_Merge(z2,wr.all,wi.all);
    put_line("The recomputed sum :"); put_line(z2);
  end Test_Two_Add;

  procedure Test_Update ( dim : in integer32 ) is

    z : constant QuadDobl_Complex_Vectors.Vector(1..dim)
      := QuadDobl_Random_Vectors.Random_Vector(1,dim);
    x : constant QuadDobl_Complex_Vectors.Vector(1..dim)
      := QuadDobl_Random_Vectors.Random_Vector(1,dim);
    z1,z2 : QuadDobl_Complex_Vectors.Vector(1..dim);
    xr,xi : Standard_Floating_Vectors.Vector(1..4*dim);
    zr,zi : Standard_Floating_Vectors.Vector(1..4*dim);
    ur,ui : Standard_Floating_Vectors.Link_to_Vector;
    wr,wi : Standard_Floating_Vectors.Link_to_Vector;
    wrk : constant Standard_Floating_Vectors.Vector(0..3) := (0..3 => 0.0);
    lwrk : constant Standard_Floating_Vectors.Link_to_Vector
         := new Standard_Floating_Vectors.Vector'(wrk);

    use QuadDobl_Complex_Vectors; -- for the + operator

  begin
    z1 := z + x;
    put_line("The update on random vectors :"); put_line(z1);
    Two_Split(z,zr,zi); Two_Split(x,xr,xi);
    ur := new Standard_Floating_Vectors.Vector'(xr);
    ui := new Standard_Floating_Vectors.Vector'(xi);
    wr := new Standard_Floating_Vectors.Vector'(zr);
    wi := new Standard_Floating_Vectors.Vector'(zi);
    Update(wr,wi,ur,ui,lwrk);
    Two_Merge(z2,wr.all,wi.all);
    put_line("The recomputed update :"); put_line(z2);
  end Test_Update;

  procedure Test_Inner_Product ( dim : in integer32 ) is

    x : constant QuadDobl_Complex_Vectors.Vector(1..dim)
      := QuadDobl_Random_Vectors.Random_Vector(1,dim);
    y : constant QuadDobl_Complex_Vectors.Vector(1..dim)
      := QuadDobl_Random_Vectors.Random_Vector(1,dim);
    z1,z2 : Complex_Number;
    xr,xi : Standard_Floating_Vectors.Vector(1..4*dim);
    yr,yi : Standard_Floating_Vectors.Vector(1..4*dim);
    zrhh,zihh,zrlh,zilh,zrhl,zihl,zrll,zill : double_float; 
    ur,ui : Standard_Floating_Vectors.Link_to_Vector;
    vr,vi : Standard_Floating_Vectors.Link_to_Vector;
    xw,yw,zw : constant Standard_Floating_Vectors.Vector(0..3)
             := (0..3 => 0.0);
    u : constant Standard_Floating_Vectors.Link_to_Vector
      := new Standard_Floating_Vectors.Vector'(xw);
    v : constant Standard_Floating_Vectors.Link_to_Vector
      := new Standard_Floating_Vectors.Vector'(yw);
    w : constant Standard_Floating_Vectors.Link_to_Vector
      := new Standard_Floating_Vectors.Vector'(zw);

  begin
    z1 := Inner_Product(x,y);
    put_line("The inner product of two random vectors :");
    put(z1); new_line;
    Two_Split(x,xr,xi); Two_Split(y,yr,yi);
    ur := new Standard_Floating_Vectors.Vector'(xr);
    ui := new Standard_Floating_Vectors.Vector'(xi);
    vr := new Standard_Floating_Vectors.Vector'(yr);
    vi := new Standard_Floating_Vectors.Vector'(yi);
    Inner_Product(zrhh,zihh,zrlh,zilh,zrhl,zihl,zrll,zill,
                  ur,ui,vr,vi,u,v,w);
    Merge(z2,zrhh,zihh,zrlh,zilh,zrhl,zihl,zrll,zill);
    put_line("The recomputed inner product :");
    put(z2); new_line;
  end Test_Inner_Product;

  procedure Test_Multiply ( deg : in integer32 ) is

    cx : constant QuadDobl_Complex_Vectors.Vector(0..deg)
       := QuadDobl_Random_Vectors.Random_Vector(0,deg);
    cy : constant QuadDobl_Complex_Vectors.Vector(0..deg)
       := QuadDobl_Random_Vectors.Random_Vector(0,deg);
    x : constant QuadDobl_Complex_Vectors.Link_to_Vector
      := new QuadDobl_Complex_Vectors.Vector'(cx);
    y : constant QuadDobl_Complex_Vectors.Link_to_Vector
      := new QuadDobl_Complex_Vectors.Vector'(cy);
    zero : constant Complex_Number := Create(integer(0));
    z1 : constant QuadDobl_Complex_Vectors.Link_to_Vector
       := new QuadDobl_Complex_Vectors.Vector'(0..deg => zero);
    z2 : QuadDobl_Complex_Vectors.Vector(0..deg);
    dim : constant integer32 := 4*(deg+1)-1;
    xr,xi : Standard_Floating_Vectors.Vector(0..dim);
    yr,yi : Standard_Floating_Vectors.Vector(0..dim);
    zr,zi : constant Standard_Floating_Vectors.Vector(0..dim)
          := (0..dim => 0.0);
    ur,ui : Standard_Floating_Vectors.Link_to_Vector;
    vr,vi : Standard_Floating_Vectors.Link_to_Vector;
    wr,wi : Standard_Floating_Vectors.Link_to_Vector;
    xw,yw,zw : constant Standard_Floating_Vectors.Vector(0..3)
             := (0..3 => 0.0);
    u : constant Standard_Floating_Vectors.Link_to_Vector
      := new Standard_Floating_Vectors.Vector'(xw);
    v : constant Standard_Floating_Vectors.Link_to_Vector
      := new Standard_Floating_Vectors.Vector'(yw);
    w : constant Standard_Floating_Vectors.Link_to_Vector
      := new Standard_Floating_Vectors.Vector'(zw);

  begin
    Multiply(x,y,z1);
    put_line("The convolution of two random vectors :"); put_line(z1);
    Two_Split(cx,xr,xi);
    Two_Split(cy,yr,yi);
    ur := new Standard_Floating_Vectors.Vector'(xr);
    ui := new Standard_Floating_Vectors.Vector'(xi);
    vr := new Standard_Floating_Vectors.Vector'(yr);
    vi := new Standard_Floating_Vectors.Vector'(yi);
    wr := new Standard_Floating_Vectors.Vector'(zr);
    wi := new Standard_Floating_Vectors.Vector'(zi);
    Multiply(ur,ui,vr,vi,wr,wi,u,v,w);
    Two_Merge(z2,wr.all,wi.all);
    put_line("The recomputed convolution :"); put_line(z2);
  end Test_Multiply;

  procedure Timing_Add ( dim,frq : in integer32 ) is

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
    xwrk,ywrk,zwrk : constant Standard_Floating_Vectors.Vector(0..3)
                   := (0..3 => 0.0);
    u : constant Standard_Floating_Vectors.Link_to_Vector
      := new Standard_Floating_Vectors.Vector'(xwrk);
    v : constant Standard_Floating_Vectors.Link_to_Vector
      := new Standard_Floating_Vectors.Vector'(ywrk);
    w : constant Standard_Floating_Vectors.Link_to_Vector
      := new Standard_Floating_Vectors.Vector'(zwrk);

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
          vrhh,vihh,vrlh,vilh,vrhl,vihl,vrll,vill,u,v,w);
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
          vrhh,vihh,vrlh,vilh,vrhl,vihl,vrll,vill,u,v,w);
      Merge(z2,wrhh.all,wihh.all,wrlh.all,wilh.all,
               wrhl.all,wihl.all,wrll.all,will.all);
    end loop;
    tstop(timer);
    new_line;
    print_times(standard_output,timer,"split, add, merge");
  end Timing_Add;

  procedure Timing_Two_Add ( dim,frq : in integer32 ) is

    timer : Timing_Widget;
    x : constant QuadDobl_Complex_Vectors.Vector(1..dim)
      := QuadDobl_Random_Vectors.Random_Vector(1,dim);
    y : constant QuadDobl_Complex_Vectors.Vector(1..dim)
      := QuadDobl_Random_Vectors.Random_Vector(1,dim);
    z1,z2 : QuadDobl_Complex_Vectors.Vector(1..dim);
    xr,xi : Standard_Floating_Vectors.Vector(1..4*dim);
    yr,yi : Standard_Floating_Vectors.Vector(1..4*dim);
    zr,zi : constant Standard_Floating_Vectors.Vector(1..4*dim)
          := (1..4*dim => 0.0);
    ur,ui : Standard_Floating_Vectors.Link_to_Vector;
    vr,vi : Standard_Floating_Vectors.Link_to_Vector;
    wr,wi : Standard_Floating_Vectors.Link_to_Vector;

    use QuadDobl_Complex_Vectors; -- for the + operator

  begin
    tstart(timer);
    for k in 1..frq loop
      z1 := x + y;
    end loop;
    tstop(timer);
    new_line;
    print_times(standard_output,timer,"complex addition");
    Two_Split(x,xr,xi); Two_Split(y,yr,yi);
    ur := new Standard_Floating_Vectors.Vector'(xr);
    ui := new Standard_Floating_Vectors.Vector'(xi);
    vr := new Standard_Floating_Vectors.Vector'(yr);
    vi := new Standard_Floating_Vectors.Vector'(yi);
    wr := new Standard_Floating_Vectors.Vector'(zr);
    wi := new Standard_Floating_Vectors.Vector'(zi);
    tstart(timer);
    for k in 1..frq loop
      Two_Add(wr,wi,ur,ui,vr,vi);
    end loop;
    tstop(timer);
    new_line;
    print_times(standard_output,timer,"2-vector addition");
    tstart(timer);
    for k in 1..frq loop
      Two_Split(x,xr,xi); Two_Split(y,yr,yi);
      Two_Add(wr,wi,ur,ui,vr,vi);
      Two_Merge(z2,wr.all,wi.all);
    end loop;
    tstop(timer);
    new_line;
    print_times(standard_output,timer,"split, add, merge");
  end Timing_Two_Add;

  procedure Timing_Update ( dim,frq : in integer32 ) is

    timer : Timing_Widget;
    z : constant QuadDobl_Complex_Vectors.Vector(1..dim)
      := QuadDobl_Random_Vectors.Random_Vector(1,dim);
    x : constant QuadDobl_Complex_Vectors.Vector(1..dim)
      := QuadDobl_Random_Vectors.Random_Vector(1,dim);
    z1,z2 : QuadDobl_Complex_Vectors.Vector(1..dim);
    xr,xi : Standard_Floating_Vectors.Vector(1..4*dim);
    zr,zi : Standard_Floating_Vectors.Vector(1..4*dim);
    ur,ui : Standard_Floating_Vectors.Link_to_Vector;
    wr,wi : Standard_Floating_Vectors.Link_to_Vector;
    wrk : constant Standard_Floating_Vectors.Vector(0..3) := (0..3 => 0.0);
    lwrk : constant Standard_Floating_Vectors.Link_to_Vector
         := new Standard_Floating_Vectors.Vector'(wrk);

    use QuadDobl_Complex_Vectors; -- for the + operator

  begin
    tstart(timer);
    for k in 1..frq loop
      z1 := z + x;
    end loop;
    tstop(timer);
    new_line;
    print_times(standard_output,timer,"complex update");
    Two_Split(z,zr,zi); Two_Split(x,xr,xi);
    ur := new Standard_Floating_Vectors.Vector'(xr);
    ui := new Standard_Floating_Vectors.Vector'(xi);
    wr := new Standard_Floating_Vectors.Vector'(zr);
    wi := new Standard_Floating_Vectors.Vector'(zi);
    Two_Split(z,zr,zi); Two_Split(x,xr,xi);
    tstart(timer);
    for k in 1..frq loop
      Update(wr,wi,ur,ui,lwrk);
    end loop;
    tstop(timer);
    new_line;
    print_times(standard_output,timer,"2-vector update");
    tstart(timer);
    for k in 1..frq loop
      Two_Split(z,zr,zi); Two_Split(x,xr,xi);
      Update(wr,wi,ur,ui,lwrk);
      Two_Merge(z2,wr.all,wi.all);
    end loop;
    tstop(timer);
    new_line;
    print_times(standard_output,timer,"split, update, merge");
  end Timing_Update;

  procedure Timing_Inner_Product ( dim,frq : in integer32 ) is

    timer : Timing_Widget;
    x : constant QuadDobl_Complex_Vectors.Vector(1..dim)
      := QuadDobl_Random_Vectors.Random_Vector(1,dim);
    y : constant QuadDobl_Complex_Vectors.Vector(1..dim)
      := QuadDobl_Random_Vectors.Random_Vector(1,dim);
    z1,z2 : Complex_Number;
    xr,xi : Standard_Floating_Vectors.Vector(1..4*dim);
    yr,yi : Standard_Floating_Vectors.Vector(1..4*dim);
    zrhh,zihh,zrlh,zilh,zrhl,zihl,zrll,zill : double_float; 
    ur,ui : Standard_Floating_Vectors.Link_to_Vector;
    vr,vi : Standard_Floating_Vectors.Link_to_Vector;
    xw,yw,zw : constant Standard_Floating_Vectors.Vector(0..3) 
             := (0..3 => 0.0);
    u : constant Standard_Floating_Vectors.Link_to_Vector
      := new Standard_Floating_Vectors.Vector'(xw);
    v : constant Standard_Floating_Vectors.Link_to_Vector
      := new Standard_Floating_Vectors.Vector'(yw);
    w : constant Standard_Floating_Vectors.Link_to_Vector
      := new Standard_Floating_Vectors.Vector'(zw);

  begin
    tstart(timer);
    for k in 1..frq loop
      z1 := Inner_Product(x,y);
    end loop;
    tstop(timer);
    new_line;
    print_times(standard_output,timer,"complex inner product");
    Two_Split(x,xr,xi); Two_Split(y,yr,yi);
    ur := new Standard_Floating_Vectors.Vector'(xr);
    ui := new Standard_Floating_Vectors.Vector'(xi);
    vr := new Standard_Floating_Vectors.Vector'(yr);
    vi := new Standard_Floating_Vectors.Vector'(yi);
    tstart(timer);
    for k in 1..frq loop
      Inner_Product(zrhh,zihh,zrlh,zilh,zrhl,zihl,zrll,zill,
                    ur,ui,vr,vi,u,v,w);
    end loop;
    tstop(timer);
    new_line;
    print_times(standard_output,timer,"2-vector inner product");
    tstart(timer);
    for k in 1..frq loop
      Two_Split(x,xr,xi); Two_Split(y,yr,yi);
      Inner_Product(zrhh,zihh,zrlh,zilh,zrhl,zihl,zrll,zill,
                    ur,ui,vr,vi,u,v,w);
      Merge(z2,zrhh,zihh,zrlh,zilh,zrhl,zihl,zrll,zill);
    end loop;
    tstop(timer);
    new_line;
    print_times(standard_output,timer,"split, inner product, merge");
  end Timing_Inner_Product;

  procedure Timing_Multiply ( deg,frq : in integer32 ) is

    timer : Timing_Widget;
    cx : constant QuadDobl_Complex_Vectors.Vector(0..deg)
       := QuadDobl_Random_Vectors.Random_Vector(0,deg);
    cy : constant QuadDobl_Complex_Vectors.Vector(0..deg)
       := QuadDobl_Random_Vectors.Random_Vector(0,deg);
    x : constant QuadDobl_Complex_Vectors.Link_to_Vector
      := new QuadDobl_Complex_Vectors.Vector'(cx);
    y : constant QuadDobl_Complex_Vectors.Link_to_Vector
      := new QuadDobl_Complex_Vectors.Vector'(cy);
    zero : constant Complex_Number := Create(integer(0));
    z1 : constant QuadDobl_Complex_Vectors.Link_to_Vector
       := new QuadDobl_Complex_Vectors.Vector'(0..deg => zero);
    z2 : QuadDobl_Complex_Vectors.Vector(0..deg);
    dim : constant integer32 := 4*(deg+1)-1;
    xr,xi : Standard_Floating_Vectors.Vector(0..dim);
    yr,yi : Standard_Floating_Vectors.Vector(0..dim);
    zr,zi : constant Standard_Floating_Vectors.Vector(0..dim)
          := (0..dim => 0.0);
    ur,ui : Standard_Floating_Vectors.Link_to_Vector;
    vr,vi : Standard_Floating_Vectors.Link_to_Vector;
    wr,wi : Standard_Floating_Vectors.Link_to_Vector;
    xw,yw,zw : constant Standard_Floating_Vectors.Vector(0..3)
             := (0..3 => 0.0);
    u : constant Standard_Floating_Vectors.Link_to_Vector
      := new Standard_Floating_Vectors.Vector'(xw);
    v : constant Standard_Floating_Vectors.Link_to_Vector
      := new Standard_Floating_Vectors.Vector'(yw);
    w : constant Standard_Floating_Vectors.Link_to_Vector
      := new Standard_Floating_Vectors.Vector'(zw);

  begin
    tstart(timer);
    for k in 1..frq loop
      Multiply(x,y,z1);
    end loop;
    tstop(timer);
    new_line;
    print_times(standard_output,timer,"complex convolution");
    Two_Split(x.all,xr,xi); Two_Split(y.all,yr,yi);
    ur := new Standard_Floating_Vectors.Vector'(xr);
    ui := new Standard_Floating_Vectors.Vector'(xi);
    vr := new Standard_Floating_Vectors.Vector'(yr);
    vi := new Standard_Floating_Vectors.Vector'(yi);
    wr := new Standard_Floating_Vectors.Vector'(zr);
    wi := new Standard_Floating_Vectors.Vector'(zi);
    tstart(timer);
    for k in 1..frq loop
      Multiply(ur,ui,vr,vi,wr,wi,u,v,w);
    end loop;
    tstop(timer);
    new_line;
    print_times(standard_output,timer,"2-vector convolution");
    tstart(timer);
    for k in 1..frq loop
      Two_Split(x.all,xr,xi); Two_Split(y.all,yr,yi);
      Multiply(ur,ui,vr,vi,wr,wi,u,v,w);
      Two_Merge(z2,wr.all,wi.all);
    end loop;
    tstop(timer);
    new_line;
    print_times(standard_output,timer,"split, convolute, merge");
  end Timing_Multiply;

  procedure Main is

    dim,frq : integer32 := 0;
    ans,tst : character;

  begin
    new_line;
    put_line("MENU to test vector operations :");
    put_line("  1. test addition on 8-vector representation");
    put_line("  2. test addition on 2-vector representation");
    put_line("  3. test update");
    put_line("  4. test inner product");
    put_line("  5. test convolution");
    put("Type 1, 2, 3, 4, or 5 to select a test : ");
    Ask_Alternative(tst,"12345");
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
        when '2' => Test_Two_Add(dim);
        when '3' => Test_Update(dim);
        when '4' => Test_Inner_Product(dim);
        when '5' => Test_Multiply(dim);
        when others => null;
      end case;
    else
      case tst is
        when '1' => Timing_Add(dim,frq);
        when '2' => Timing_Two_Add(dim,frq);
        when '3' => Timing_Update(dim,frq);
        when '4' => Timing_Inner_Product(dim,frq);
        when '5' => Timing_Multiply(dim,frq);
        when others => null;
      end case;
    end if;
  end Main;

end Test_QuadDobl_Performance;
