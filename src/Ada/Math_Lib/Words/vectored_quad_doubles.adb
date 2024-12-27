with text_io;                            use text_io;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;
with Bits_of_Doubles;
with Sign_Balancers;                     use Sign_Balancers;

package body Vectored_Quad_Doubles is

-- BASIC PROCEDURES :

  function Sign_Balance
             ( x : Quad_Double_Vectors.Vector; verbose : boolean := true )
             return Quad_Double_Vectors.Vector is

    res : Quad_Double_Vectors.Vector(x'range);

  begin
    for i in x'range loop
      res(i) := x(i);
      if not Is_Sign_Balanced(res(i))
       then Sign_Balance(res(i),verbose);
      end if;
    end loop;
    return res;
  end Sign_Balance;

  procedure Signed_Quarter
              ( x,y : in Quad_Double_Vectors.Vector;
                xs0,xs1,xs2,xs3 : out Standard_Floating_Vectors.Vector;
                xs4,xs5,xs6,xs7 : out Standard_Floating_Vectors.Vector;
                xs8,xs9,xsA,xsB : out Standard_Floating_Vectors.Vector;
                xsC,xsD,xsE,xsF : out Standard_Floating_Vectors.Vector;
                ys0,ys1,ys2,ys3 : out Standard_Floating_Vectors.Vector;
                ys4,ys5,ys6,ys7 : out Standard_Floating_Vectors.Vector;
                ys8,ys9,ysA,ysB : out Standard_Floating_Vectors.Vector;
                ysC,ysD,ysE,ysF : out Standard_Floating_Vectors.Vector;
                xd0,xd1,xd2,xd3 : out Standard_Floating_Vectors.Vector;
                xd4,xd5,xd6,xd7 : out Standard_Floating_Vectors.Vector;
                xd8,xd9,xdA,xdB : out Standard_Floating_Vectors.Vector;
                xdC,xdD,xdE,xdF : out Standard_Floating_Vectors.Vector;
                yd0,yd1,yd2,yd3 : out Standard_Floating_Vectors.Vector;
                yd4,yd5,yd6,yd7 : out Standard_Floating_Vectors.Vector;
                yd8,yd9,ydA,ydB : out Standard_Floating_Vectors.Vector;
                ydC,ydD,ydE,ydF : out Standard_Floating_Vectors.Vector;
                ns,nd : out integer32 ) is

    nbr : quad_double;
    flt : double_float;
    x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,xA,xB,xC,xD,xE,xF : double_float;
    y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,yA,yB,yC,yD,yE,yF : double_float;

  begin
    ns := 0; nd := 0;
    for i in x'range loop
      nbr := x(i);
      flt := hihi_part(nbr); Bits_of_Doubles.Split(flt,x0,x1,x2,x3);
      flt := lohi_part(nbr); Bits_of_Doubles.Split(flt,x4,x5,x6,x7);
      flt := hilo_part(nbr); Bits_of_Doubles.Split(flt,x8,x9,xA,xB);
      flt := lolo_part(nbr); Bits_of_Doubles.Split(flt,xC,xD,xE,xF);
      nbr := y(i);
      flt := hihi_part(nbr); Bits_of_Doubles.Split(flt,y0,y1,y2,y3);
      flt := lohi_part(nbr); Bits_of_Doubles.Split(flt,y4,y5,y6,y7);
      flt := hilo_part(nbr); Bits_of_Doubles.Split(flt,y8,y9,yA,yB);
      flt := lolo_part(nbr); Bits_of_Doubles.Split(flt,yC,yD,yE,yF);
      if not Sign_Balancers.Different_Sign(x0,y0) then
        ns := ns + 1;
        xs0(ns) := x0; xs1(ns) := x1;
        xs2(ns) := x2; xs3(ns) := x3;
        xs4(ns) := x4; xs5(ns) := x5;
        xs6(ns) := x6; xs7(ns) := x7;
        xs8(ns) := x8; xs9(ns) := x9;
        xsA(ns) := xA; xsB(ns) := xB;
        xsC(ns) := xC; xsD(ns) := xD;
        xsE(ns) := xE; xsF(ns) := xF;
        ys0(ns) := y0; ys1(ns) := y1;
        ys2(ns) := y2; ys3(ns) := y3;
        ys4(ns) := y4; ys5(ns) := y5;
        ys6(ns) := y6; ys7(ns) := y7;
        ys8(ns) := y8; ys9(ns) := y9;
        ysA(ns) := yA; ysB(ns) := yB;
        ysC(ns) := yC; ysD(ns) := yD;
        ysE(ns) := yE; ysF(ns) := yF;
      else -- Different_Sign(x0,y0)
        nd := nd + 1;
        xd0(nd) := x0; xd1(nd) := x1;
        xd2(nd) := x2; xd3(nd) := x3;
        xd4(nd) := x4; xd5(nd) := x5;
        xd6(nd) := x6; xd7(nd) := x7;
        xd8(nd) := x8; xd9(nd) := x9;
        xdA(nd) := xA; xdB(nd) := xB;
        xdC(nd) := xC; xdD(nd) := xD;
        xdE(nd) := xE; xdF(nd) := xF;
        yd0(nd) := y0; yd1(nd) := y1;
        yd2(nd) := y2; yd3(nd) := y3;
        yd4(nd) := y4; yd5(nd) := y5;
        yd6(nd) := y6; yd7(nd) := y7;
        yd8(nd) := y8; yd9(nd) := y9;
        ydA(nd) := yA; ydB(nd) := yB;
        ydC(nd) := yC; ydD(nd) := yD;
        ydE(nd) := yE; ydF(nd) := yF;
      end if;
    end loop;
  end Signed_Quarter;

  procedure Write_Subsums
              ( s0,s1,s2,s3,s4,s5,s6,s7 : in double_float ) is

    use Bits_of_Doubles;

  begin
    put("s0 : "); put(s0);
    put(", n0 : "); put(Last_Zero_Count(s0),1); new_line;
    put("s1 : "); put(s1);
    put(", n1 : "); put(Last_Zero_Count(s1),1); new_line;
    put("s2 : "); put(s2);
    put(", n2 : "); put(Last_Zero_Count(s2),1); new_line;
    put("s3 : "); put(s3);
    put(", n3 : "); put(Last_Zero_Count(s3),1); new_line;
    put("s4 : "); put(s4);
    put(", n4 : "); put(Last_Zero_Count(s4),1); new_line;
    put("s5 : "); put(s5);
    put(", n5 : "); put(Last_Zero_Count(s5),1); new_line;
    put("s6 : "); put(s6);
    put(", n6 : "); put(Last_Zero_Count(s6),1); new_line;
    put("s7 : "); put(s7);
    put(", n7 : "); put(Last_Zero_Count(s7),1); new_line;
  end Write_Subsums;

  procedure Write_Subsums
              ( s0,s1,s2,s3,s4,s5,s6,s7 : in double_float;
                s8,s9,sA,sB,sC,sD,sE,sF : in double_float ) is

    use Bits_of_Doubles;

  begin
    put("s0 : "); put(s0);
    put(", n0 : "); put(Last_Zero_Count(s0),1); new_line;
    put("s1 : "); put(s1);
    put(", n1 : "); put(Last_Zero_Count(s1),1); new_line;
    put("s2 : "); put(s2);
    put(", n2 : "); put(Last_Zero_Count(s2),1); new_line;
    put("s3 : "); put(s3);
    put(", n3 : "); put(Last_Zero_Count(s3),1); new_line;
    put("s4 : "); put(s4);
    put(", n4 : "); put(Last_Zero_Count(s4),1); new_line;
    put("s5 : "); put(s5);
    put(", n5 : "); put(Last_Zero_Count(s5),1); new_line;
    put("s6 : "); put(s6);
    put(", n6 : "); put(Last_Zero_Count(s6),1); new_line;
    put("s7 : "); put(s7);
    put(", n7 : "); put(Last_Zero_Count(s7),1); new_line;
    put("s8 : "); put(s8);
    put(", n8 : "); put(Last_Zero_Count(s8),1); new_line;
    put("s9 : "); put(s9);
    put(", n9 : "); put(Last_Zero_Count(s9),1); new_line;
    put("sA : "); put(sA);
    put(", nA : "); put(Last_Zero_Count(sA),1); new_line;
    put("sB : "); put(sB);
    put(", nB : "); put(Last_Zero_Count(sB),1); new_line;
    put("sC : "); put(sC);
    put(", nC : "); put(Last_Zero_Count(sC),1); new_line;
    put("sD : "); put(sD);
    put(", nD : "); put(Last_Zero_Count(sD),1); new_line;
    put("sE : "); put(sE);
    put(", nE : "); put(Last_Zero_Count(sE),1); new_line;
    put("sF : "); put(sF);
    put(", nF : "); put(Last_Zero_Count(sF),1); new_line;
  end Write_Subsums;

  function to_quad_double
              ( s0,s1,s2,s3,s4,s5,s6,s7 : double_float;
                verbose : boolean := true ) return quad_double is

    res : quad_double;

  begin
    if verbose
     then write_subsums(s0,s1,s2,s3,s4,s5,s6,s7);
    end if;
    res := create(s7);
    res := res + create(s6);
    res := res + create(s5);
    res := res + create(s4);
    res := res + create(s3);
    res := res + create(s2);
    res := res + create(s1);
    res := res + create(s0);
    return res;
  end to_quad_double;

  function to_quad_double
              ( s0,s1,s2,s3,s4,s5,s6,s7 : double_float;
                s8,s9,sA,sB,sC,sD,sE,sF : double_float;
                verbose : boolean := true ) return quad_double is

    res : quad_double;

  begin
    if verbose
     then write_subsums(s0,s1,s2,s3,s4,s5,s6,s7,s8,s9,sA,sB,sC,sD,sE,sF);
    end if;
    res := create(sF);
    res := res + create(sE);
    res := res + create(sD);
    res := res + create(sC);
    res := res + create(sB);
    res := res + create(sA);
    res := res + create(s9);
    res := res + create(s8);
    res := res + create(s7);
    res := res + create(s6);
    res := res + create(s5);
    res := res + create(s4);
    res := res + create(s3);
    res := res + create(s2);
    res := res + create(s1);
    res := res + create(s0);
    return res;
  end to_quad_double;

-- SIGN AWARE WRAPPERS :

  function Product ( x,y : Quad_Double_Vectors.Vector;
                     verbose : boolean := true ) return quad_double is

    res : quad_double := create(0.0);
    xb : constant Quad_Double_Vectors.Vector(x'range)
       := Sign_Balance(x,verbose);
    yb : constant Quad_Double_Vectors.Vector(y'range)
       := Sign_Balance(y,verbose);
    xs0,xs1,xs2,xs3 : Standard_Floating_Vectors.Vector(x'range);
    xs4,xs5,xs6,xs7 : Standard_Floating_Vectors.Vector(x'range);
    xs8,xs9,xsA,xsB : Standard_Floating_Vectors.Vector(x'range);
    xsC,xsD,xsE,xsF : Standard_Floating_Vectors.Vector(x'range);
    ys0,ys1,ys2,ys3 : Standard_Floating_Vectors.Vector(y'range);
    ys4,ys5,ys6,ys7 : Standard_Floating_Vectors.Vector(y'range);
    ys8,ys9,ysA,ysB : Standard_Floating_Vectors.Vector(y'range);
    ysC,ysD,ysE,ysF : Standard_Floating_Vectors.Vector(y'range);
    xd0,xd1,xd2,xd3 : Standard_Floating_Vectors.Vector(x'range);
    xd4,xd5,xd6,xd7 : Standard_Floating_Vectors.Vector(x'range);
    xd8,xd9,xdA,xdB : Standard_Floating_Vectors.Vector(x'range);
    xdC,xdD,xdE,xdF : Standard_Floating_Vectors.Vector(x'range);
    yd0,yd1,yd2,yd3 : Standard_Floating_Vectors.Vector(y'range);
    yd4,yd5,yd6,yd7 : Standard_Floating_Vectors.Vector(y'range);
    yd8,yd9,ydA,ydB : Standard_Floating_Vectors.Vector(y'range);
    ydC,ydD,ydE,ydF : Standard_Floating_Vectors.Vector(y'range);
    ns,nd : integer32;
    s0,s1,s2,s3,s4,s5,s6,s7 : double_float;
    s8a,s9a,sAa,sBa,sCa,sDa,sEa,sFa : double_float;
    s8b,s9b,sAb,sBb,sCb,sDb,sEb,sFb : double_float;

  begin
    Signed_Quarter(xb,yb,xs0,xs1,xs2,xs3,xs4,xs5,xs6,xs7,
                         xs8,xs9,xsA,xsB,xsC,xsD,xsE,xsF,
                         ys0,ys1,ys2,ys3,ys4,ys5,ys6,ys7,
                         ys8,ys9,ysA,ysB,ysC,ysD,ysE,ysF,
                         xd0,xd1,xd2,xd3,xd4,xd5,xd6,xd7,
                         xd8,xd9,xdA,xdB,xdC,xdD,xdE,xdF,
                         yd0,yd1,yd2,yd3,yd4,yd5,yd6,yd7,
                         yd8,yd9,ydA,ydB,ydC,ydD,ydE,ydF,ns,nd);
    if verbose then
      put("#s : "); put(ns,1); 
      put(", #d : "); put(nd,1); new_line;
    end if;
    s0 := 0.0; s1 := 0.0; s2 := 0.0; s3 := 0.0;
    s4 := 0.0; s5 := 0.0; s6 := 0.0; s7 := 0.0;
    s8a := 0.0; s9a := 0.0; sAa := 0.0; sBa := 0.0;
    sCa := 0.0; sDa := 0.0; sEa := 0.0; sFa := 0.0;
    s8b := 0.0; s9b := 0.0; sAb := 0.0; sBb := 0.0;
    sCb := 0.0; sDb := 0.0; sEb := 0.0; sFb := 0.0;
    for i in 1..ns loop
      s0 := s0 + xs0(i)*ys0(i);
      s1 := s1 + xs0(i)*ys1(i) + xs1(i)*ys0(i);
      s2 := s2 + xs0(i)*ys2(i) + xs1(i)*ys1(i) + xs2(i)*ys0(i);
      s3 := s3 + xs0(i)*ys3(i) + xs1(i)*ys2(i) + xs2(i)*ys1(i)
               + xs3(i)*ys0(i);
      s4 := s4 + xs0(i)*ys4(i) + xs1(i)*ys3(i) + xs2(i)*ys2(i)
               + xs3(i)*ys1(i) + xs4(i)*ys0(i);
      s5 := s5 + xs0(i)*ys5(i) + xs1(i)*ys4(i) + xs2(i)*ys3(i)
               + xs3(i)*ys2(i) + xs4(i)*ys1(i) + xs5(i)*ys0(i);
      s6 := s6 + xs0(i)*ys6(i) + xs1(i)*ys5(i) + xs2(i)*ys4(i)
               + xs3(i)*ys3(i) + xs4(i)*ys2(i) + xs5(i)*ys1(i)
               + xs6(i)*ys0(i);
      s7 := s7 + xs0(i)*ys7(i) + xs1(i)*ys6(i) + xs2(i)*ys5(i)
               + xs3(i)*ys4(i) + xs4(i)*ys3(i) + xs5(i)*ys2(i)
               + xs6(i)*ys1(i) + xs7(i)*ys0(i);
      s8a := s8a + xs0(i)*ys8(i) + xs1(i)*ys7(i) + xs2(i)*ys6(i)
                 + xs3(i)*ys5(i) + xs4(i)*ys4(i);
      s8b := s8b + xs5(i)*ys3(i) + xs6(i)*ys2(i) + xs7(i)*ys1(i)
                 + xs8(i)*ys0(i);
      s9a := s9a + xs0(i)*ys9(i) + xs1(i)*ys8(i) + xs2(i)*ys7(i)
                 + xs3(i)*ys6(i) + xs4(i)*ys5(i);
      s9b := s9b + xs5(i)*ys4(i) + xs6(i)*ys3(i) + xs7(i)*ys2(i)
                 + xs8(i)*ys1(i) + xs9(i)*ys0(i);
      sAa := sAa + xs0(i)*ysA(i) + xs1(i)*ys9(i) + xs2(i)*ys8(i)
                 + xs3(i)*ys7(i) + xs4(i)*ys6(i) + xs5(i)*ys5(i);
      sAb := sAb + xs6(i)*ys4(i) + xs7(i)*ys3(i) + xs8(i)*ys2(i)
                 + xs9(i)*ys1(i) + xsA(i)*ys0(i);
      sBa := sBa + xs0(i)*ysB(i) + xs1(i)*ysA(i) + xs2(i)*ys9(i)
                 + xs3(i)*ys8(i) + xs4(i)*ys7(i) + xs5(i)*ys6(i);
      sBb := sBb + xs6(i)*ys5(i) + xs7(i)*ys4(i) + xs8(i)*ys3(i)
                 + xs9(i)*ys2(i) + xsA(i)*ys1(i) + xsB(i)*ys0(i);
      sCa := sCa + xs0(i)*ysC(i) + xs1(i)*ysB(i) + xs2(i)*ysA(i)
                 + xs3(i)*ys9(i) + xs4(i)*ys8(i) + xs5(i)*ys7(i)
                 + xs6(i)*ys6(i);
      sCb := sCb + xs7(i)*ys5(i) + xs8(i)*ys4(i) + xs9(i)*ys3(i)
                 + xsA(i)*ys2(i) + xsB(i)*ys1(i) + xsC(i)*ys0(i);
      sDa := sDa + xs0(i)*ysD(i) + xs1(i)*ysC(i) + xs2(i)*ysB(i)
                 + xs3(i)*ysA(i) + xs4(i)*ys9(i) + xs5(i)*ys8(i)
                 + xs6(i)*ys7(i);
      sDb := sDb + xs7(i)*ys6(i) + xs8(i)*ys5(i) + xs9(i)*ys4(i)
                 + xsA(i)*ys3(i) + xsB(i)*ys2(i) + xsC(i)*ys1(i)
                 + xsD(i)*ys0(i);
      sEa := sEa + xs0(i)*ysE(i) + xs1(i)*ysD(i) + xs2(i)*ysC(i)
                 + xs3(i)*ysB(i) + xs4(i)*ysA(i) + xs5(i)*ys9(i)
                 + xs6(i)*ys8(i) + xs7(i)*ys7(i);
      sEb := sEb + xs8(i)*ys6(i) + xs9(i)*ys5(i) + xsA(i)*ys4(i)
                 + xsB(i)*ys3(i) + xsC(i)*ys2(i) + xsD(i)*ys1(i)
                 + xsE(i)*ys0(i);
      sFa := sFa + xs0(i)*ysF(i) + xs1(i)*ysE(i) + xs2(i)*ysD(i)
                 + xs3(i)*ysC(i) + xs4(i)*ysB(i) + xs5(i)*ysA(i)
                 + xs6(i)*ys9(i) + xs7(i)*ys8(i);
      sFb := sFb + xs8(i)*ys7(i) + xs9(i)*ys6(i) + xsA(i)*ys5(i)
                 + xsB(i)*ys4(i) + xsC(i)*ys3(i) + xsD(i)*ys2(i)
                 + xsE(i)*ys1(i) + xsF(i)*ys0(i);
    end loop;
    res := to_quad_double(s0,s1,s2,s3,s4,s5,s6,s7,
                          s8a,s9a,sAa,sBa,sCa,sDa,sEa,sFa,verbose=>true)
         + to_quad_double(s8b,s9b,sAb,sBb,sCb,sDb,sEb,sFb,verbose=>true);
    s0 := 0.0; s1 := 0.0; s2 := 0.0; s3 := 0.0;
    s4 := 0.0; s5 := 0.0; s6 := 0.0; s7 := 0.0;
    s8a := 0.0; s9a := 0.0; sAa := 0.0; sBa := 0.0;
    sCa := 0.0; sDa := 0.0; sEa := 0.0; sFa := 0.0;
    s8b := 0.0; s9b := 0.0; sAb := 0.0; sBb := 0.0;
    sCb := 0.0; sDb := 0.0; sEb := 0.0; sFb := 0.0;
    for i in 1..nd loop
      s0 := s0 + xd0(i)*yd0(i);
      s1 := s1 + xd0(i)*yd1(i) + xd1(i)*yd0(i);
      s2 := s2 + xd0(i)*yd2(i) + xd1(i)*yd1(i) + xd2(i)*yd0(i);
      s3 := s3 + xd0(i)*yd3(i) + xd1(i)*yd2(i) + xd2(i)*yd1(i)
               + xd3(i)*yd0(i);
      s4 := s4 + xd0(i)*yd4(i) + xd1(i)*yd3(i) + xd2(i)*yd2(i)
               + xd3(i)*yd1(i) + xd4(i)*yd0(i);
      s5 := s5 + xd0(i)*yd5(i) + xd1(i)*yd4(i) + xd2(i)*yd3(i)
               + xd3(i)*yd2(i) + xd4(i)*yd1(i) + xd5(i)*yd0(i);
      s6 := s6 + xd0(i)*yd6(i) + xd1(i)*yd5(i) + xd2(i)*yd4(i)
               + xd3(i)*yd3(i) + xd4(i)*yd2(i) + xd5(i)*yd1(i)
               + xd6(i)*yd0(i);
      s7 := s7 + xd0(i)*yd7(i) + xd1(i)*yd6(i) + xd2(i)*yd5(i)
               + xd3(i)*yd4(i) + xd4(i)*yd3(i) + xd5(i)*yd2(i)
               + xd6(i)*yd1(i) + xd7(i)*yd0(i);
      s8a := s8a + xd0(i)*yd8(i) + xd1(i)*yd7(i) + xd2(i)*yd6(i)
                 + xd3(i)*yd5(i) + xd4(i)*yd4(i);
      s8b := s8b + xd5(i)*yd3(i) + xd6(i)*yd2(i) + xd7(i)*yd1(i)
                 + xd8(i)*yd0(i);
      s9a := s9a + xd0(i)*yd9(i) + xd1(i)*yd8(i) + xd2(i)*yd7(i)
                 + xd3(i)*yd6(i) + xd4(i)*yd5(i);
      s9b := s9b + xd5(i)*yd4(i) + xd6(i)*yd3(i) + xd7(i)*yd2(i)
                 + xd8(i)*yd1(i) + xd9(i)*yd0(i);
      sAa := sAa + xd0(i)*ydA(i) + xd1(i)*yd9(i) + xd2(i)*yd8(i)
                 + xd3(i)*yd7(i) + xd4(i)*yd6(i) + xd5(i)*yd5(i);
      sAb := sAb + xd6(i)*yd4(i) + xd7(i)*yd3(i) + xd8(i)*yd2(i)
                 + xd9(i)*yd1(i) + xdA(i)*yd0(i);
      sBa := sBa + xd0(i)*ydB(i) + xd1(i)*ydA(i) + xd2(i)*yd9(i)
                 + xd3(i)*yd8(i) + xd4(i)*yd7(i) + xd5(i)*yd6(i);
      sBb := sBb + xd6(i)*yd5(i) + xd7(i)*yd4(i) + xd8(i)*yd3(i)
                 + xd9(i)*yd2(i) + xdA(i)*yd1(i) + xdB(i)*yd0(i);
      sCa := sCa + xd0(i)*ydC(i) + xd1(i)*ydB(i) + xd2(i)*ydA(i)
                 + xd3(i)*yd9(i) + xd4(i)*yd8(i) + xd5(i)*yd7(i)
                 + xd6(i)*yd6(i);
      sCb := sCb + xd7(i)*yd5(i) + xd8(i)*yd4(i) + xd9(i)*yd3(i)
                 + xdA(i)*yd2(i) + xdB(i)*yd1(i) + xdC(i)*yd0(i);
      sDa := sDa + xd0(i)*ydD(i) + xd1(i)*ydC(i) + xd2(i)*ydB(i)
                 + xd3(i)*ydA(i) + xd4(i)*yd9(i) + xd5(i)*yd8(i)
                 + xd6(i)*yd7(i);
      sDb := sDb + xd7(i)*yd6(i) + xd8(i)*yd5(i) + xd9(i)*yd4(i)
                 + xdA(i)*yd3(i) + xdB(i)*yd2(i) + xdC(i)*yd1(i)
                 + xdD(i)*yd0(i);
      sEa := sEa + xd0(i)*ydE(i) + xd1(i)*ydD(i) + xd2(i)*ydC(i)
                 + xd3(i)*ydB(i) + xd4(i)*ydA(i) + xd5(i)*yd9(i)
                 + xd6(i)*yd8(i) + xd7(i)*yd7(i);
      sEb := sEb + xd8(i)*yd6(i) + xd9(i)*yd5(i) + xdA(i)*yd4(i)
                 + xdB(i)*yd3(i) + xdC(i)*yd2(i) + xdD(i)*yd1(i)
                 + xdE(i)*yd0(i);
      sFa := sFa + xd0(i)*ydF(i) + xd1(i)*ydE(i) + xd2(i)*ydD(i)
                 + xd3(i)*ydC(i) + xd4(i)*ydB(i) + xd5(i)*ydA(i)
                 + xd6(i)*yd9(i) + xd7(i)*yd8(i);
      sFb := sFb + xd8(i)*yd7(i) + xd9(i)*yd6(i) + xdA(i)*yd5(i)
                 + xdB(i)*yd4(i) + xdC(i)*yd3(i) + xdD(i)*yd2(i)
                 + xdE(i)*yd1(i) + xdF(i)*yd0(i);
    end loop;
    res := res
         + to_quad_double(s0,s1,s2,s3,s4,s5,s6,s7,
                          s8a,s9a,sAa,sBa,sCa,sDa,sEa,sFa,verbose=>true)
         + to_quad_double(s8b,s9b,sAb,sBb,sCb,sDb,sEb,sFb,verbose=>true);
    return res;
  end Product;

  function Product ( x,y : QuadDobl_Complex_Vectors.Vector;
                     verbose : boolean := true ) return Complex_Number is

    res : Complex_Number;
    resre,resim : quad_double;
    xre,xim : Quad_Double_Vectors.Vector(x'range);
    yre,yim : Quad_Double_Vectors.Vector(y'range);

  begin
    for i in x'range loop
      xre(i) := QuadDobl_Complex_Numbers.REAL_PART(x(i));
      xim(i) := QuadDobl_Complex_Numbers.IMAG_PART(x(i));
      yre(i) := QuadDobl_Complex_Numbers.REAL_PART(y(i));
      yim(i) := QuadDobl_Complex_Numbers.IMAG_PART(y(i));
    end loop;
    resre := Product(xre,yre,verbose) - Product(xim,yim,verbose);
    resim := Product(xre,yim,verbose) + Product(xim,yre,verbose);
    res := create(resre,resim);
    return res;
  end Product;

end Vectored_Quad_Doubles;
