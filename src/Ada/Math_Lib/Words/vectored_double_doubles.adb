with text_io;                            use text_io;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;
with Double_Double_Basics;
with Bits_of_Doubles;

package body Vectored_Double_Doubles is

  procedure Split ( v : in Double_Double_Vectors.Vector;
                    v0,v1,v2,v3 : out Standard_Floating_Vectors.Vector ) is

    nbr : double_double;
    flt : double_float;

  begin
    for i in v'range loop
      nbr := v(i);
      flt := hi_part(nbr);
      Double_Double_Basics.Split(flt,v0(i),v1(i));
      flt := lo_part(nbr);
      Double_Double_Basics.Split(flt,v2(i),v3(i));
    end loop;
  end Split;

  procedure Signed_Split
              ( v : in Double_Double_Vectors.Vector;
                p0,p1,p2,p3 : out Standard_Floating_Vectors.Vector;
                m0,m1,m2,m3 : out Standard_Floating_Vectors.Vector;
                np0,np1,np2,np3,nm0,nm1,nm2,nm3 : out integer32 ) is

    nbr : double_double;
    flt,hi,lo : double_float;

  begin
    np0 := 0; np1 := 0; np2 := 0; np3 := 0;
    nm0 := 0; nm1 := 0; nm2 := 0; nm3 := 0;
    for i in v'range loop
      nbr := v(i);
      flt := hi_part(nbr);
      Double_Double_Basics.Split(flt,hi,lo);
      if hi >= 0.0 then
        np0 := np0 + 1;
        p0(np0) := hi;
      else
        nm0 := nm0 + 1;
        m0(nm0) := hi;
      end if;
      if lo >= 0.0 then
        np1 := np1 + 1;
        p1(np1) := lo;
      else
        nm1 := nm1 + 1;
        m1(nm1) := lo;
      end if;
      flt := lo_part(nbr);
      Double_Double_Basics.Split(flt,hi,lo);
      if hi >= 0.0 then
        np2 := np2 + 1;
        p2(np2) := hi;
      else
        nm2 := nm2 + 1;
        m2(nm2) := hi;
      end if;
      if lo >= 0.0 then
        np3 := np3 + 1;
        p3(np3) := lo;
      else
        nm3 := nm3 + 1;
        m3(nm3) := lo;
      end if;
    end loop;
  end Signed_Split;

  procedure Split ( v : in DoblDobl_Complex_Vectors.Vector;
                    v0re,v1re : out Standard_Floating_Vectors.Vector;
                    v2re,v3re : out Standard_Floating_Vectors.Vector;
                    v0im,v1im : out Standard_Floating_Vectors.Vector;
                    v2im,v3im : out Standard_Floating_Vectors.Vector ) is

     nbr : double_double;
     flt : double_float;

  begin
    for i in v'range loop
      nbr := DoblDobl_Complex_Numbers.REAL_PART(v(i));
      flt := hi_part(nbr);
      Double_Double_Basics.Split(flt,v0re(i),v1re(i));
      flt := lo_part(nbr);
      Double_Double_Basics.Split(flt,v2re(i),v3re(i));
      nbr := DoblDobl_Complex_Numbers.IMAG_PART(v(i));
      flt := hi_part(nbr);
      Double_Double_Basics.Split(flt,v0im(i),v1im(i));
      flt := lo_part(nbr);
      Double_Double_Basics.Split(flt,v2im(i),v3im(i));
    end loop;
  end Split;

  procedure Quarter ( v : in Double_Double_Vectors.Vector;
                      v0,v1,v2,v3 : out Standard_Floating_Vectors.Vector;
                      v4,v5,v6,v7 : out Standard_Floating_Vectors.Vector ) is

    nbr : double_double;
    flt : double_float;

  begin
    for i in v'range loop
      nbr := v(i);
      flt := hi_part(nbr);
      Bits_of_Doubles.Split(flt,v0(i),v1(i),v2(i),v3(i));
      flt := lo_part(nbr);
      Bits_of_Doubles.Split(flt,v4(i),v5(i),v6(i),v7(i));
    end loop;
  end Quarter;

  procedure Quarter ( v : in DoblDobl_Complex_Vectors.Vector;
                      v0re,v1re : out Standard_Floating_Vectors.Vector;
                      v2re,v3re : out Standard_Floating_Vectors.Vector;
                      v4re,v5re : out Standard_Floating_Vectors.Vector;
                      v6re,v7re : out Standard_Floating_Vectors.Vector;
                      v0im,v1im : out Standard_Floating_Vectors.Vector;
                      v2im,v3im : out Standard_Floating_Vectors.Vector;
                      v4im,v5im : out Standard_Floating_Vectors.Vector;
                      v6im,v7im : out Standard_Floating_Vectors.Vector ) is

    nbr : double_double;
    flt : double_float;

  begin
    for i in v'range loop
      nbr := DoblDobl_Complex_Numbers.REAL_PART(v(i));
      flt := hi_part(nbr);
      Bits_of_Doubles.Split(flt,v0re(i),v1re(i),v2re(i),v3re(i));
      flt := lo_part(nbr);
      Bits_of_Doubles.Split(flt,v4re(i),v5re(i),v6re(i),v7re(i));
      nbr := DoblDobl_Complex_Numbers.IMAG_PART(v(i));
      flt := hi_part(nbr);
      Bits_of_Doubles.Split(flt,v0im(i),v1im(i),v2im(i),v3im(i));
      flt := lo_part(nbr);
      Bits_of_Doubles.Split(flt,v4im(i),v5im(i),v6im(i),v7im(i));
    end loop;
  end Quarter;

  function Signs ( f0,f1,f2,f3,f4,f5,f6,f7 : double_float ) return string is

    res : String(1..8);

  begin
    if f0 >= 0.0
     then res(1) := '+';
     else res(1) := '-';
    end if;
    if f1 >= 0.0
     then res(2) := '+';
     else res(2) := '-';
    end if;
    if f2 >= 0.0
     then res(3) := '+';
     else res(3) := '-';
    end if;
    if f3 >= 0.0
     then res(4) := '+';
     else res(4) := '-';
    end if;
    if f4 >= 0.0
     then res(5) := '+';
     else res(5) := '-';
    end if;
    if f5 >= 0.0
     then res(6) := '+';
     else res(6) := '-';
    end if;
    if f6 >= 0.0
     then res(7) := '+';
     else res(7) := '-';
    end if;
    if f7 >= 0.0
     then res(8) := '+';
     else res(8) := '-';
    end if;
    return res;
  end Signs;

  procedure Signed_Convolutions ( sx : in string ) is
  begin
    put("s0 := s0 + ("); put_line(sx(sx'first) & "," & sx(sx'first) & ")");
    put("s1 := s1 + (");
    put(sx(sx'first) & "," & sx(sx'first+1));
    put_line(") + (" & sx(sx'first+1) & "," & sx(sx'first) & ")");
    put("s2 := s2 + (");
    put(sx(sx'first) & "," & sx(sx'first+2));
    put(") + (" & sx(sx'first+1) & "," & sx(sx'first+1));
    put_line(") + (" & sx(sx'first+2) & "," & sx(sx'first) & ")");
    put("s3 := s3 + (");
    put(sx(sx'first) & "," & sx(sx'first+3));
    put(") + (" & sx(sx'first+1) & "," & sx(sx'first+2));
    put(") + (" & sx(sx'first+2) & "," & sx(sx'first+1));
    put_line(") + (" & sx(sx'first+3) & "," & sx(sx'first) & ")");
    put("s4 := s4 + (");
    put(sx(sx'first) & "," & sx(sx'first+4));
    put(") + (" & sx(sx'first+1) & "," & sx(sx'first+3));
    put(") + (" & sx(sx'first+2) & "," & sx(sx'first+2));
    put(") + (" & sx(sx'first+3) & "," & sx(sx'first+1));
    put_line(") + (" & sx(sx'first+4) & "," & sx(sx'first) & ")");
    put("s5 := s5 + (");
    put(sx(sx'first) & "," & sx(sx'first+5));
    put(") + (" & sx(sx'first+1) & "," & sx(sx'first+4));
    put(") + (" & sx(sx'first+2) & "," & sx(sx'first+3));
    put(") + (" & sx(sx'first+3) & "," & sx(sx'first+2));
    put(") + (" & sx(sx'first+4) & "," & sx(sx'first+1));
    put_line(") + (" & sx(sx'first+5) & "," & sx(sx'first) & ")");
    put("s6 := s6 + (");
    put(sx(sx'first) & "," & sx(sx'first+6));
    put(") + (" & sx(sx'first+1) & "," & sx(sx'first+5));
    put(") + (" & sx(sx'first+2) & "," & sx(sx'first+4));
    put(") + (" & sx(sx'first+3) & "," & sx(sx'first+3));
    put(") + (" & sx(sx'first+4) & "," & sx(sx'first+2));
    put(") + (" & sx(sx'first+5) & "," & sx(sx'first+1));
    put_line(") + (" & sx(sx'first+6) & "," & sx(sx'first) & ")");
    put("s7 := s7 + (");
    put(sx(sx'first) & "," & sx(sx'first+7));
    put(") + (" & sx(sx'first+1) & "," & sx(sx'first+6));
    put(") + (" & sx(sx'first+2) & "," & sx(sx'first+5));
    put(") + (" & sx(sx'first+3) & "," & sx(sx'first+4));
    put(") + (" & sx(sx'first+4) & "," & sx(sx'first+3));
    put(") + (" & sx(sx'first+5) & "," & sx(sx'first+2));
    put(") + (" & sx(sx'first+6) & "," & sx(sx'first+1));
    put_line(") + (" & sx(sx'first+7) & "," & sx(sx'first) & ")");
  end Signed_Convolutions;

  procedure Signed_Convolutions ( sx,sy : in string ) is
  begin
    put("s0 := s0 + ("); put_line(sx(sx'first) & "," & sy(sy'first) & ")");
    put("s1 := s1 + (");
    put(sx(sx'first) & "," & sy(sy'first+1));
    put_line(") + (" & sx(sx'first+1) & "," & sy(sy'first) & ")");
    put("s2 := s2 + (");
    put(sx(sx'first) & "," & sy(sy'first+2));
    put(") + (" & sx(sx'first+1) & "," & sy(sy'first+1));
    put_line(") + (" & sx(sx'first+2) & "," & sy(sy'first) & ")");
    put("s3 := s3 + (");
    put(sx(sx'first) & "," & sy(sy'first+3));
    put(") + (" & sx(sx'first+1) & "," & sy(sy'first+2));
    put(") + (" & sx(sx'first+2) & "," & sy(sy'first+1));
    put_line(") + (" & sx(sx'first+3) & "," & sy(sy'first) & ")");
    put("s4 := s4 + (");
    put(sx(sx'first) & "," & sy(sy'first+4));
    put(") + (" & sx(sx'first+1) & "," & sy(sy'first+3));
    put(") + (" & sx(sx'first+2) & "," & sy(sy'first+2));
    put(") + (" & sx(sx'first+3) & "," & sy(sy'first+1));
    put_line(") + (" & sx(sx'first+4) & "," & sy(sy'first) & ")");
    put("s5 := s5 + (");
    put(sx(sx'first) & "," & sy(sy'first+5));
    put(") + (" & sx(sx'first+1) & "," & sy(sy'first+4));
    put(") + (" & sx(sx'first+2) & "," & sy(sy'first+3));
    put(") + (" & sx(sx'first+3) & "," & sy(sy'first+2));
    put(") + (" & sx(sx'first+4) & "," & sy(sy'first+1));
    put_line(") + (" & sx(sx'first+5) & "," & sy(sy'first) & ")");
    put("s6 := s6 + (");
    put(sx(sx'first) & "," & sy(sy'first+6));
    put(") + (" & sx(sx'first+1) & "," & sy(sy'first+5));
    put(") + (" & sx(sx'first+2) & "," & sy(sy'first+4));
    put(") + (" & sx(sx'first+3) & "," & sy(sy'first+3));
    put(") + (" & sx(sx'first+4) & "," & sy(sy'first+2));
    put(") + (" & sx(sx'first+5) & "," & sy(sy'first+1));
    put_line(") + (" & sx(sx'first+6) & "," & sy(sy'first) & ")");
    put("s7 := s7 + (");
    put(sx(sx'first) & "," & sy(sy'first+7));
    put(") + (" & sx(sx'first+1) & "," & sy(sy'first+6));
    put(") + (" & sx(sx'first+2) & "," & sy(sy'first+5));
    put(") + (" & sx(sx'first+3) & "," & sy(sy'first+4));
    put(") + (" & sx(sx'first+4) & "," & sy(sy'first+3));
    put(") + (" & sx(sx'first+5) & "," & sy(sy'first+2));
    put(") + (" & sx(sx'first+6) & "," & sy(sy'first+1));
    put_line(") + (" & sx(sx'first+7) & "," & sy(sy'first) & ")");
  end Signed_Convolutions;

  procedure Signed_Quarter
              ( v : in Double_Double_Vectors.Vector;
                x0,x1,x2,x3 : out Standard_Floating_Vectors.Vector;
                x4,x5,x6,x7 : out Standard_Floating_Vectors.Vector;
                pm0,pm1,pm2,pm3 : out Standard_Floating_Vectors.Vector;
                pm4,pm5,pm6,pm7 : out Standard_Floating_Vectors.Vector;
                mp0,mp1,mp2,mp3 : out Standard_Floating_Vectors.Vector;
                mp4,mp5,mp6,mp7 : out Standard_Floating_Vectors.Vector;
                nbx,npm,nmp : out integer32;
                verbose : in boolean := true ) is

    nbr : double_double;
    flt,f0,f1,f2,f3,f4,f5,f6,f7 : double_float;
    sx : String(1..8);

  begin
    nbx := 0; npm := 0; nmp := 0;
    for i in v'range loop
      nbr := v(i);
      flt := hi_part(nbr);
      Bits_of_Doubles.Split(flt,f0,f1,f2,f3);
      flt := lo_part(nbr);
      Bits_of_Doubles.Split(flt,f4,f5,f6,f7);
      if f0 >= 0.0 then
        if f4 >= 0.0 then
          nbx := nbx + 1;
          x0(nbx) := f0; x1(nbx) := f1; x2(nbx) := f2; x3(nbx) := f3;
          x4(nbx) := f4; x5(nbx) := f5; x6(nbx) := f6; x7(nbx) := f7;
        else -- f4 < 0.0
          npm := npm + 1;
          pm0(npm) := f0; pm1(npm) := f1; pm2(npm) := f2; pm3(npm) := f3;
          pm4(npm) := f4; pm5(npm) := f5; pm6(npm) := f6; pm7(npm) := f7;
        end if;
      else -- f0 < 0.0 
        if f4 < 0.0 then
          nbx := nbx + 1;
          x0(nbx) := f0; x1(nbx) := f1; x2(nbx) := f2; x3(nbx) := f3;
          x4(nbx) := f4; x5(nbx) := f5; x6(nbx) := f6; x7(nbx) := f7;
        else -- f0 >= 0.0
          nmp := nmp + 1;
          mp0(nmp) := f0; mp1(nmp) := f1; mp2(nmp) := f2; mp3(nmp) := f3;
          mp4(nmp) := f4; mp5(nmp) := f5; mp6(nmp) := f6; mp7(nmp) := f7;
        end if;
      end if;
      if verbose then
        sx := Signs(f0,f1,f2,f3,f4,f5,f6,f7);
        put_line("sx : " & sx);
        Signed_Convolutions(sx);
      end if;
    end loop;
  end Signed_Quarter;

  procedure Signed_Quarter
              ( x,y : in Double_Double_Vectors.Vector;
                xss0,xss1,xss2,xss3 : out Standard_Floating_Vectors.Vector;
                xss4,xss5,xss6,xss7 : out Standard_Floating_Vectors.Vector;
                yss0,yss1,yss2,yss3 : out Standard_Floating_Vectors.Vector;
                yss4,yss5,yss6,yss7 : out Standard_Floating_Vectors.Vector;
                xsd0,xsd1,xsd2,xsd3 : out Standard_Floating_Vectors.Vector;
                xsd4,xsd5,xsd6,xsd7 : out Standard_Floating_Vectors.Vector;
                ysd0,ysd1,ysd2,ysd3 : out Standard_Floating_Vectors.Vector;
                ysd4,ysd5,ysd6,ysd7 : out Standard_Floating_Vectors.Vector;
                xds0,xds1,xds2,xds3 : out Standard_Floating_Vectors.Vector;
                xds4,xds5,xds6,xds7 : out Standard_Floating_Vectors.Vector;
                yds0,yds1,yds2,yds3 : out Standard_Floating_Vectors.Vector;
                yds4,yds5,yds6,yds7 : out Standard_Floating_Vectors.Vector;
                xdd0,xdd1,xdd2,xdd3 : out Standard_Floating_Vectors.Vector;
                xdd4,xdd5,xdd6,xdd7 : out Standard_Floating_Vectors.Vector;
                ydd0,ydd1,ydd2,ydd3 : out Standard_Floating_Vectors.Vector;
                ydd4,ydd5,ydd6,ydd7 : out Standard_Floating_Vectors.Vector;
                nss,nsd,nds,ndd : out integer32;
                verbose : in boolean := true ) is

    nbr : double_double;
    flt,x0,x1,x2,x3,x4,x5,x6,x7,y0,y1,y2,y3,y4,y5,y6,y7 : double_float;
    sx,sy : String(1..8);

  begin
    nss := 0; nsd := 0; nds := 0; ndd := 0;
    for i in x'range loop
      nbr := x(i);
      flt := hi_part(nbr);
      Bits_of_Doubles.Split(flt,x0,x1,x2,x3);
      flt := lo_part(nbr);
      Bits_of_Doubles.Split(flt,x4,x5,x6,x7);
      nbr := y(i);
      flt := hi_part(nbr);
      Bits_of_Doubles.Split(flt,y0,y1,y2,y3);
      flt := lo_part(nbr);
      Bits_of_Doubles.Split(flt,y4,y5,y6,y7);
      if x0 >= 0.0 then
        if y0 >= 0.0 then
          if x4 >= 0.0 then
            if y4 >= 0.0 then
              nss := nss + 1;
              xss0(nss) := x0; xss1(nss) := x1;
              xss2(nss) := x2; xss3(nss) := x3;
              xss4(nss) := x4; xss5(nss) := x5;
              xss6(nss) := x6; xss7(nss) := x7;
              yss0(nss) := y0; yss1(nss) := y1;
              yss2(nss) := y2; yss3(nss) := y3;
              yss4(nss) := y4; yss5(nss) := y5;
              yss6(nss) := y6; yss7(nss) := y7;
            else -- y4 < 0.0 then
              nsd := nsd + 1;
              xsd0(nsd) := x0; xsd1(nsd) := x1;
              xsd2(nsd) := x2; xsd3(nsd) := x3;
              xsd4(nsd) := x4; xsd5(nsd) := x5;
              xsd6(nsd) := x6; xsd7(nsd) := x7;
              ysd0(nsd) := y0; ysd1(nsd) := y1;
              ysd2(nsd) := y2; ysd3(nsd) := y3;
              ysd4(nsd) := y4; ysd5(nsd) := y5;
              ysd6(nsd) := y6; ysd7(nsd) := y7;
            end if;
          else -- x4 < 0.0
            if y4 >= 0.0 then
              nsd := nsd + 1;
              xsd0(nsd) := x0; xsd1(nsd) := x1;
              xsd2(nsd) := x2; xsd3(nsd) := x3;
              xsd4(nsd) := x4; xsd5(nsd) := x5;
              xsd6(nsd) := x6; xsd7(nsd) := x7;
              ysd0(nsd) := y0; ysd1(nsd) := y1;
              ysd2(nsd) := y2; ysd3(nsd) := y3;
              ysd4(nsd) := y4; ysd5(nsd) := y5;
              ysd6(nsd) := y6; ysd7(nsd) := y7;
            else -- y4 < 0.0
              nss := nss + 1;
              xss0(nss) := x0; xss1(nss) := x1;
              xss2(nss) := x2; xss3(nss) := x3;
              xss4(nss) := x4; xss5(nss) := x5;
              xss6(nss) := x6; xss7(nss) := x7;
              yss0(nss) := y0; yss1(nss) := y1;
              yss2(nss) := y2; yss3(nss) := y3;
              yss4(nss) := y4; yss5(nss) := y5;
              yss6(nss) := y6; yss7(nss) := y7;
            end if;
          end if;
        else -- x0 >= 0.0, y0 < 0.0
          if x4 >= 0.0 then
            if y4 >= 0.0 then
              nds := nds + 1;
              xds0(nds) := x0; xds1(nds) := x1;
              xds2(nds) := x2; xds3(nds) := x3;
              xds4(nds) := x4; xds5(nds) := x5;
              xds6(nds) := x6; xds7(nds) := x7;
              yds0(nds) := y0; yds1(nds) := y1;
              yds2(nds) := y2; yds3(nds) := y3;
              yds4(nds) := y4; yds5(nds) := y5;
              yds6(nds) := y6; yds7(nds) := y7;
            else -- x4 >= 0.0, y4 < 0.0
              ndd := ndd + 1;
              xdd0(ndd) := x0; xdd1(ndd) := x1;
              xdd2(ndd) := x2; xdd3(ndd) := x3;
              xdd4(ndd) := x4; xdd5(ndd) := x5;
              xdd6(ndd) := x6; xdd7(ndd) := x7;
              ydd0(ndd) := y0; ydd1(ndd) := y1;
              ydd2(ndd) := y2; ydd3(ndd) := y3;
              ydd4(ndd) := y4; ydd5(ndd) := y5;
              ydd6(ndd) := y6; ydd7(ndd) := y7;
            end if;
          else -- x0 >= 0.0, y0 < 0.0, x4 < 0.0 
            if y4 >= 0.0 then
              ndd := ndd + 1;
              xdd0(ndd) := x0; xdd1(ndd) := x1;
              xdd2(ndd) := x2; xdd3(ndd) := x3;
              xdd4(ndd) := x4; xdd5(ndd) := x5;
              xdd6(ndd) := x6; xdd7(ndd) := x7;
              ydd0(ndd) := y0; ydd1(ndd) := y1;
              ydd2(ndd) := y2; ydd3(ndd) := y3;
              ydd4(ndd) := y4; ydd5(ndd) := y5;
              ydd6(ndd) := y6; ydd7(ndd) := y7;
            else -- y4 < 0.0
              nds := nds + 1;
              xds0(nds) := x0; xds1(nds) := x1;
              xds2(nds) := x2; xds3(nds) := x3;
              xds4(nds) := x4; xds5(nds) := x5;
              xds6(nds) := x6; xds7(nds) := x7;
              yds0(nds) := y0; yds1(nds) := y1;
              yds2(nds) := y2; yds3(nds) := y3;
              yds4(nds) := y4; yds5(nds) := y5;
              yds6(nds) := y6; yds7(nds) := y7;
            end if;
          end if;
        end if;
      else -- x0 < 0.0
        if y0 >= 0.0 then
          if x4 >= 0.0 then
            if y4 >= 0.0 then
              nds := nds + 1;
              xds0(nds) := x0; xds1(nds) := x1;
              xds2(nds) := x2; xds3(nds) := x3;
              xds4(nds) := x4; xds5(nds) := x5;
              xds6(nds) := x6; xds7(nds) := x7;
              yds0(nds) := y0; yds1(nds) := y1;
              yds2(nds) := y2; yds3(nds) := y3;
              yds4(nds) := y4; yds5(nds) := y5;
              yds6(nds) := y6; yds7(nds) := y7;
            else -- y4 < 0.0 then
              ndd := ndd + 1;
              xdd0(ndd) := x0; xdd1(ndd) := x1;
              xdd2(ndd) := x2; xdd3(ndd) := x3;
              xdd4(ndd) := x4; xdd5(ndd) := x5;
              xdd6(ndd) := x6; xdd7(ndd) := x7;
              ydd0(ndd) := y0; ydd1(ndd) := y1;
              ydd2(ndd) := y2; ydd3(ndd) := y3;
              ydd4(ndd) := y4; ydd5(ndd) := y5;
              ydd6(ndd) := y6; ydd7(ndd) := y7;
            end if;
          else -- x0 < 0.0, y0 >= 0.0, x4 < 0.0
            if y4 >= 0.0 then
              ndd := ndd + 1;
              xdd0(ndd) := x0; xdd1(ndd) := x1;
              xdd2(ndd) := x2; xdd3(ndd) := x3;
              xdd4(ndd) := x4; xdd5(ndd) := x5;
              xdd6(ndd) := x6; xdd7(ndd) := x7;
              ydd0(ndd) := y0; ydd1(ndd) := y1;
              ydd2(ndd) := y2; ydd3(ndd) := y3;
              ydd4(ndd) := y4; ydd5(ndd) := y5;
              ydd6(ndd) := y6; ydd7(ndd) := y7;
            else -- y4 < 0.0
              nds := nds + 1;
              xds0(nds) := x0; xds1(nds) := x1;
              xds2(nds) := x2; xds3(nds) := x3;
              xds4(nds) := x4; xds5(nds) := x5;
              xds6(nds) := x6; xds7(nds) := x7;
              yds0(nds) := y0; yds1(nds) := y1;
              yds2(nds) := y2; yds3(nds) := y3;
              yds4(nds) := y4; yds5(nds) := y5;
              yds6(nds) := y6; yds7(nds) := y7;
            end if;
          end if;
        else -- x0 < 0.0, y0 < 0.0
          if x4 >= 0.0 then
            if y4 >= 0.0 then
              nss := nss + 1;
              xss0(nss) := x0; xss1(nss) := x1;
              xss2(nss) := x2; xss3(nss) := x3;
              xss4(nss) := x4; xss5(nss) := x5;
              xss6(nss) := x6; xss7(nss) := x7;
              yss0(nss) := y0; yss1(nss) := y1;
              yss2(nss) := y2; yss3(nss) := y3;
              yss4(nss) := y4; yss5(nss) := y5;
              yss6(nss) := y6; yss7(nss) := y7;
            else -- x4 >= 0.0, y4 < 0.0
              nsd := nsd + 1;
              xsd0(nsd) := x0; xsd1(nsd) := x1;
              xsd2(nsd) := x2; xsd3(nsd) := x3;
              xsd4(nsd) := x4; xsd5(nsd) := x5;
              xsd6(nsd) := x6; xsd7(nsd) := x7;
              ysd0(nsd) := y0; ysd1(nsd) := y1;
              ysd2(nsd) := y2; ysd3(nsd) := y3;
              ysd4(nsd) := y4; ysd5(nsd) := y5;
              ysd6(nsd) := y6; ysd7(nsd) := y7;
            end if;
          else -- x0 < 0.0, y0 < 0.0, x4 < 0.0 
            if y4 >= 0.0 then
              nsd := nsd + 1;
              xsd0(nsd) := x0; xsd1(nsd) := x1;
              xsd2(nsd) := x2; xsd3(nsd) := x3;
              xsd4(nsd) := x4; xsd5(nsd) := x5;
              xsd6(nsd) := x6; xsd7(nsd) := x7;
              ysd0(nsd) := y0; ysd1(nsd) := y1;
              ysd2(nsd) := y2; ysd3(nsd) := y3;
              ysd4(nsd) := y4; ysd5(nsd) := y5;
              ysd6(nsd) := y6; ysd7(nsd) := y7;
            else -- y4 < 0.0
              nss := nss + 1;
              xss0(nss) := x0; xss1(nss) := x1;
              xss2(nss) := x2; xss3(nss) := x3;
              xss4(nss) := x4; xss5(nss) := x5;
              xss6(nss) := x6; xss7(nss) := x7;
              yss0(nss) := y0; yss1(nss) := y1;
              yss2(nss) := y2; yss3(nss) := y3;
              yss4(nss) := y4; yss5(nss) := y5;
              yss6(nss) := y6; yss7(nss) := y7;
            end if;
          end if;
        end if;
      end if;
      if verbose then
        sx := Signs(x0,x1,x2,x3,x4,x5,x6,x7);
        sy := Signs(y0,y1,y2,y3,y4,y5,y6,y7);
        put_line("sx : " & sx & ", sy : " & sy);
        Signed_Convolutions(sx,sy);
      end if;
    end loop;
  end Signed_Quarter;

  procedure Sum ( v0,v1,v2,v3 : in Standard_Floating_Vectors.Vector;
                  s0,s1,s2,s3 : out double_float ) is
  begin
    s0 := 0.0;
    s1 := 0.0;
    s2 := 0.0;
    s3 := 0.0;
    for i in v0'range loop
      s0 := s0 + v0(i);
      s1 := s1 + v1(i);
      s2 := s2 + v2(i);
      s3 := s3 + v3(i);
    end loop;
  end Sum;

  procedure Sum ( v0re,v1re,v2re,v3re : in Standard_Floating_Vectors.Vector;
                  v0im,v1im,v2im,v3im : in Standard_Floating_Vectors.Vector;
                  s0re,s1re,s2re,s3re : out double_float;
                  s0im,s1im,s2im,s3im : out double_float ) is
  begin
    s0re := 0.0;
    s1re := 0.0;
    s2re := 0.0;
    s3re := 0.0;
    s0im := 0.0;
    s1im := 0.0;
    s2im := 0.0;
    s3im := 0.0;
    for i in v0re'range loop
      s0re := s0re + v0re(i);
      s1re := s1re + v1re(i);
      s2re := s2re + v2re(i);
      s3re := s3re + v3re(i);
      s0im := s0im + v0im(i);
      s1im := s1im + v1im(i);
      s2im := s2im + v2im(i);
      s3im := s3im + v3im(i);
    end loop;
  end Sum;

  procedure Product ( x0,x1,x2,x3 : in Standard_Floating_Vectors.Vector;
                      x4,x5,x6,x7 : in Standard_Floating_Vectors.Vector;
                      y0,y1,y2,y3 : in Standard_Floating_Vectors.Vector;
                      y4,y5,y6,y7 : in Standard_Floating_Vectors.Vector;
                      s0,s1,s2,s3,s4,s5,s6,s7 : out double_float ) is
  begin
    s0 := 0.0;
    s1 := 0.0;
    s2 := 0.0;
    s3 := 0.0;
    s4 := 0.0;
    s5 := 0.0;
    s6 := 0.0;
    s7 := 0.0;
    for i in x0'range loop
      s0 := s0 + x0(i)*y0(i);
      s1 := s1 + x0(i)*y1(i) + x1(i)*y0(i);
      s2 := s2 + x0(i)*y2(i) + x1(i)*y1(i) + x2(i)*y0(i);
      s3 := s3 + x0(i)*y3(i) + x1(i)*y2(i) + x2(i)*y1(i) + x3(i)*y0(i);
      s4 := s4 + x0(i)*y4(i) + x1(i)*y3(i) + x2(i)*y2(i) + x3(i)*y1(i)
               + x4(i)*y0(i);
      s5 := s5 + x0(i)*y5(i) + x1(i)*y4(i) + x2(i)*y3(i) + x3(i)*y2(i)
               + x4(i)*y1(i) + x5(i)*y0(i);
      s6 := s6 + x0(i)*y6(i) + x1(i)*y5(i) + x2(i)*y4(i) + x3(i)*y3(i)
               + x4(i)*y2(i) + x5(i)*y1(i) + x6(i)*y0(i);
      s7 := s7 + x0(i)*y7(i) + x1(i)*y6(i) + x2(i)*y5(i) + x3(i)*y4(i)
               + x4(i)*y3(i) + x5(i)*y2(i) + x6(i)*y1(i) + x7(i)*y0(i);
    end loop;
  end Product;

  procedure Product ( x0re,x1re : in Standard_Floating_Vectors.Vector;
                      x2re,x3re : in Standard_Floating_Vectors.Vector;
                      x4re,x5re : in Standard_Floating_Vectors.Vector;
                      x6re,x7re : in Standard_Floating_Vectors.Vector;
                      x0im,x1im : in Standard_Floating_Vectors.Vector;
                      x2im,x3im : in Standard_Floating_Vectors.Vector;
                      x4im,x5im : in Standard_Floating_Vectors.Vector;
                      x6im,x7im : in Standard_Floating_Vectors.Vector;
                      y0re,y1re : in Standard_Floating_Vectors.Vector;
                      y2re,y3re : in Standard_Floating_Vectors.Vector;
                      y4re,y5re : in Standard_Floating_Vectors.Vector;
                      y6re,y7re : in Standard_Floating_Vectors.Vector;
                      y0im,y1im : in Standard_Floating_Vectors.Vector;
                      y2im,y3im : in Standard_Floating_Vectors.Vector;
                      y4im,y5im : in Standard_Floating_Vectors.Vector;
                      y6im,y7im : in Standard_Floating_Vectors.Vector;
                      s0re,s1re,s2re,s3re : out double_float;
                      s4re,s5re,s6re,s7re : out double_float;
                      s0im,s1im,s2im,s3im : out double_float;
                      s4im,s5im,s6im,s7im : out double_float ) is
  begin
    s0re := 0.0; s1re := 0.0; s2re := 0.0; s3re := 0.0;
    s4re := 0.0; s5re := 0.0; s6re := 0.0; s7re := 0.0;
    s0im := 0.0; s1im := 0.0; s2im := 0.0; s3im := 0.0;
    s4im := 0.0; s5im := 0.0; s6im := 0.0; s7im := 0.0;
    for i in x0re'range loop
      s0re := s0re + x0re(i)*y0re(i) - x0im(i)*y0im(i);
      s1re := s1re + x0re(i)*y1re(i) - x0im(i)*y1im(i)
                   + x1re(i)*y0re(i) - x1im(i)*y0im(i);
      s2re := s2re + x0re(i)*y2re(i) - x0im(i)*y2im(i)
                   + x1re(i)*y1re(i) - x1im(i)*y1im(i)
                   + x2re(i)*y0re(i) - x2im(i)*y0im(i);
      s3re := s3re + x0re(i)*y3re(i) - x0im(i)*y3im(i)
                   + x1re(i)*y2re(i) - x1im(i)*y2im(i)
                   + x2re(i)*y1re(i) - x2im(i)*y1im(i)
                   + x3re(i)*y0re(i) - x3im(i)*y0im(i);
      s4re := s4re + x0re(i)*y4re(i) - x0im(i)*y4im(i)
                   + x1re(i)*y3re(i) - x1im(i)*y3im(i)
                   + x2re(i)*y2re(i) - x2im(i)*y2im(i)
                   + x3re(i)*y1re(i) - x3im(i)*y1im(i)
                   + x4re(i)*y0re(i) - x4im(i)*y0im(i);
      s5re := s5re + x0re(i)*y5re(i) - x0im(i)*y5im(i)
                   + x1re(i)*y4re(i) - x1im(i)*y4im(i)
                   + x2re(i)*y3re(i) - x2im(i)*y3im(i)
                   + x3re(i)*y2re(i) - x3im(i)*y2im(i)
                   + x4re(i)*y1re(i) - x4im(i)*y1im(i)
                   + x5re(i)*y0re(i) - x5im(i)*y0im(i);
      s6re := s6re + x0re(i)*y6re(i) - x0im(i)*y6im(i)
                   + x1re(i)*y5re(i) - x1im(i)*y5im(i)
                   + x2re(i)*y4re(i) - x2im(i)*y4im(i)
                   + x3re(i)*y3re(i) - x3im(i)*y3im(i)
                   + x4re(i)*y2re(i) - x4im(i)*y2im(i)
                   + x5re(i)*y1re(i) - x5im(i)*y1im(i)
                   + x6re(i)*y0re(i) - x6im(i)*y0im(i);
      s7re := s7re + x0re(i)*y7re(i) - x0im(i)*y7im(i)
                   + x1re(i)*y6re(i) - x1im(i)*y6im(i)
                   + x2re(i)*y5re(i) - x2im(i)*y5im(i)
                   + x3re(i)*y4re(i) - x3im(i)*y4im(i)
                   + x4re(i)*y3re(i) - x4im(i)*y3im(i)
                   + x5re(i)*y2re(i) - x5im(i)*y2im(i)
                   + x6re(i)*y1re(i) - x6im(i)*y1im(i)
                   + x7re(i)*y0re(i) - x7im(i)*y0im(i);
      s0im := s0im + x0re(i)*y0im(i) + x0im(i)*y0re(i);
      s1im := s1im + x0re(i)*y1im(i) + x1im(i)*y0re(i)
                   + x1re(i)*y0im(i) + x0im(i)*y1re(i);
      s2im := s2im + x0re(i)*y2im(i) + x0im(i)*y2re(i)
                   + x1re(i)*y1im(i) + x1im(i)*y1re(i)
                   + x2re(i)*y0im(i) + x2im(i)*y0re(i);
      s3im := s3im + x0re(i)*y3im(i) + x0im(i)*y3re(i)
                   + x1re(i)*y2im(i) + x1im(i)*y2re(i)
                   + x2re(i)*y1im(i) + x2im(i)*y1re(i)
                   + x3re(i)*y0im(i) + x3im(i)*y0re(i);
      s4im := s4im + x0re(i)*y4im(i) + x0im(i)*y4re(i)
                   + x1re(i)*y3im(i) + x1im(i)*y3re(i)
                   + x2re(i)*y2im(i) + x2im(i)*y2re(i)
                   + x3re(i)*y1im(i) + x3im(i)*y1re(i)
                   + x4re(i)*y0im(i) + x4im(i)*y0re(i);
      s5im := s5im + x0re(i)*y5im(i) + x0im(i)*y5re(i)
                   + x1re(i)*y4im(i) + x1im(i)*y4re(i)
                   + x2re(i)*y3im(i) + x2im(i)*y3re(i)
                   + x3re(i)*y2im(i) + x3im(i)*y2re(i)
                   + x4re(i)*y1im(i) + x4im(i)*y1re(i)
                   + x5re(i)*y0im(i) + x5im(i)*y0re(i);
      s6im := s6im + x0re(i)*y6im(i) + x0im(i)*y6re(i)
                   + x1re(i)*y5im(i) + x1im(i)*y5re(i)
                   + x2re(i)*y4im(i) + x2im(i)*y4re(i)
                   + x3re(i)*y3im(i) + x3im(i)*y3re(i)
                   + x4re(i)*y2im(i) + x4im(i)*y2re(i)
                   + x5re(i)*y1im(i) + x5im(i)*y1re(i)
                   + x6re(i)*y0im(i) + x6im(i)*y0re(i);
      s7im := s7im + x0re(i)*y7im(i) + x0im(i)*y7re(i)
                   + x1re(i)*y6im(i) + x1im(i)*y6re(i)
                   + x2re(i)*y5im(i) + x2im(i)*y5re(i)
                   + x3re(i)*y4im(i) + x3im(i)*y4re(i)
                   + x4re(i)*y3im(i) + x4im(i)*y3re(i)
                   + x5re(i)*y2im(i) + x5im(i)*y2re(i)
                   + x6re(i)*y1im(i) + x6im(i)*y1re(i)
                   + x7re(i)*y0im(i) + x7im(i)*y0re(i);
    end loop;
  end Product;

  function to_double_double
             ( s0,s1,s2,s3 : double_float;
               verbose : boolean := true ) return double_double is

    shi,slo,err : double_float;
    res : double_double;

  begin
    Double_Double_Basics.quick_two_sum(s0,s1,shi,err);
    res := create(shi,err);
    if verbose then
      put("shi : "); put(shi); new_line;
      put("err : "); put(err); new_line;
    end if;
    Double_Double_Basics.quick_two_sum(s2,s3,slo,err);
    res := res + create(slo,err);
    if verbose then
      put("slo : "); put(slo); new_line;
      put("err : "); put(err); new_line;
    end if;
    return res;
  end to_double_double;

  function to_Double_Double
             ( s0,s1,s2,s3,s4,s5,s6,s7 : double_float;
               verbose : boolean := true ) return double_double is

    res : double_double;
    z0,z1,z2,z3 : double_float;
    e1,e2,e3,e4 : double_float;
    dz0,dz1,dz2,dz3 : double_double;

  begin
    Double_Double_Basics.two_sum(s0,s1,z0,e1);
    if verbose then
      put(" z0 : "); put(z0); new_line;
      put("err : "); put(e1); new_line;
    end if;
    dz0 := create(z0,e1);
    Double_Double_Basics.two_sum(s2,s3,z1,e2);
    if verbose then
      put(" z1 : "); put(z1); new_line;
      put("err : "); put(e2); new_line;
    end if;
    dz1 := create(z1,e2);
    Double_Double_Basics.two_sum(s4,s5,z2,e3);
    if verbose then
      put(" z2 : "); put(z2); new_line;
      put("err : "); put(e3); new_line;
    end if;
    dz2 := create(z2,e3);
    Double_Double_Basics.two_sum(s6,s7,z3,e4);
    if verbose then
      put(" z3 : "); put(z3); new_line;
      put("err : "); put(e4); new_line;
    end if;
    dz3 := create(z3,e4);
    res := ((dz3 + dz2) + dz1) + dz0;
    return res;
  end to_Double_Double;

  function to_Complex_Double_Double
             ( s0re,s1re,s2re,s3re,s0im,s1im,s2im,s3im : double_float;
               verbose : boolean := true ) return Complex_Number is

    res : Complex_Number;
    sre : constant double_double
        := to_Double_Double(s0re,s1re,s2re,s3re,verbose);
    sim : constant double_double
        := to_Double_Double(s0im,s1im,s2im,s3im,verbose);

  begin
    res := Create(sre,sim);
    return res;
  end to_Complex_Double_Double;

  function to_Complex_Double_Double
             ( s0re,s1re,s2re,s3re,s4re,s5re,s6re,s7re : double_float;
               s0im,s1im,s2im,s3im,s4im,s5im,s6im,s7im : double_float;
               verbose : boolean := true ) return Complex_Number is

    res : Complex_Number;
    sre : constant double_double
        := to_Double_Double(s0re,s1re,s2re,s3re,s4re,s5re,s6re,s7re,verbose);
    sim : constant double_double
        := to_Double_Double(s0im,s1im,s2im,s3im,s4im,s5im,s6im,s7im,verbose);

  begin
    res := Create(sre,sim);
    return res;
  end to_Complex_Double_Double;

-- SIGN AWARE WRAPPERS :

  function Sum ( v : Double_Double_Vectors.Vector;
                 verbose : boolean := true ) return double_double is

    res : double_double := create(0.0);
    dim : constant integer32 := v'length;
    p0,p1,p2,p3 : Standard_Floating_Vectors.Vector(1..dim);
    m0,m1,m2,m3 : Standard_Floating_Vectors.Vector(1..dim);
    np0,np1,np2,np3,nm0,nm1,nm2,nm3 : integer32;
    s0,s1,s2,s3,s4,s5,s6,s7 : double_float;

  begin
    Signed_Split(v,p0,p1,p2,p3,m0,m1,m2,m3,np0,np1,np2,np3,nm0,nm1,nm2,nm3);
    if verbose then
      put("np0 : "); put(np0,1); 
      put(", np1 : "); put(np1,1); 
      put(", np2 : "); put(np2,1); 
      put(", np3 : "); put(np3,1); new_line;
      put("nm0 : "); put(nm0,1); 
      put(", nm1 : "); put(nm1,1); 
      put(", nm2 : "); put(nm2,1); 
      put(", nm3 : "); put(nm3,1); new_line;
    end if;
    s0 := 0.0; s1 := 0.0; s2 := 0.0; s3 := 0.0;
    for i in 1..np0 loop
      s0 := s0 + p0(i);
    end loop;
    for i in 1..np1 loop
      s1 := s1 + p1(i);
    end loop;
    for i in 1..np2 loop
      s2 := s2 + p2(i);
    end loop;
    for i in 1..np3 loop
      s3 := s3 + p3(i);
    end loop;
    res := to_double_double(s0,s1,s2,s3,verbose);
    s4 := 0.0; s5 := 0.0; s6 := 0.0; s7 := 0.0;
    for i in 1..nm0 loop
      s4 := s4 + m0(i);
    end loop;
    for i in 1..nm1 loop
      s5 := s5 + m1(i);
    end loop;
    for i in 1..nm2 loop
      s6 := s6 + m2(i);
    end loop;
    for i in 1..nm3 loop
      s7 := s7 + m3(i);
    end loop;
    res := res + to_double_double(s4,s5,s6,s7,verbose);
    return res;
  end Sum;

  function Sum ( v : DoblDobl_Complex_Vectors.Vector;
                 verbose : boolean := true ) return Complex_Number is

    res : Complex_Number;
    resre,resim : double_double;
    vre,vim : Double_Double_Vectors.Vector(v'range);

  begin
    for i in v'range loop
      vre(i) := REAL_PART(v(i));
      vim(i) := IMAG_PART(v(i));
    end loop;
    resre := Sum(vre,verbose);
    resim := Sum(vim,verbose);
    res := Create(resre,resim);
    return res;
  end Sum;

  function Squared_Norm
             ( x : Double_Double_Vectors.Vector;
               verbose : boolean := true ) return double_double is

    res : double_double;
    dim : constant integer32 := x'length;
    x0,x1,x2,x3 : Standard_Floating_Vectors.Vector(1..dim);
    x4,x5,x6,x7 : Standard_Floating_Vectors.Vector(1..dim);
    pm0,pm1,pm2,pm3 : Standard_Floating_Vectors.Vector(1..dim);
    pm4,pm5,pm6,pm7 : Standard_Floating_Vectors.Vector(1..dim);
    mp0,mp1,mp2,mp3 : Standard_Floating_Vectors.Vector(1..dim);
    mp4,mp5,mp6,mp7 : Standard_Floating_Vectors.Vector(1..dim);
    nbx,npm,nmp : integer32;
    s0,s1,s2,s3,s4,s5,s6,s7 : double_float;

  begin
    Signed_Quarter(x,x0,x1,x2,x3,x4,x5,x6,x7,
                   pm0,pm1,pm2,pm3,pm4,pm5,pm6,pm7,
                   mp0,mp1,mp2,mp3,mp4,mp5,mp6,mp7,nbx,npm,nmp,verbose);
    if verbose then
      put("#x : "); put(nbx,1); 
      put(", #pm : "); put(npm,1); 
      put(", #mp : "); put(nmp,1); new_line;
    end if;
    s0 := 0.0; s1 := 0.0; s2 := 0.0; s3 := 0.0;
    s4 := 0.0; s5 := 0.0; s6 := 0.0; s7 := 0.0;
    for i in 1..nbx loop
      s0 := s0 + x0(i)*x0(i);
      s1 := s1 + x0(i)*x1(i) + x1(i)*x0(i);
      s2 := s2 + x0(i)*x2(i) + x1(i)*x1(i) + x2(i)*x0(i);
      s3 := s3 + x0(i)*x3(i) + x1(i)*x2(i) + x2(i)*x1(i) + x3(i)*x0(i);
      s4 := s4 + x0(i)*x4(i) + x1(i)*x3(i) + x2(i)*x2(i) + x3(i)*x1(i)
               + x4(i)*x0(i);
      s5 := s5 + x0(i)*x5(i) + x1(i)*x4(i) + x2(i)*x3(i) + x3(i)*x2(i)
               + x4(i)*x1(i) + x5(i)*x0(i);
      s6 := s6 + x0(i)*x6(i) + x1(i)*x5(i) + x2(i)*x4(i) + x3(i)*x3(i)
               + x4(i)*x2(i) + x5(i)*x1(i) + x6(i)*x0(i);
      s7 := s7 + x0(i)*x7(i) + x1(i)*x6(i) + x2(i)*x5(i) + x3(i)*x4(i)
               + x4(i)*x3(i) + x5(i)*x2(i) + x6(i)*x1(i) + x7(i)*x0(i);
    end loop;
    res := to_double_double(s0,s1,s2,s3,s4,s5,s6,s7,verbose);
    s0 := 0.0; s1 := 0.0; s2 := 0.0; s3 := 0.0;
    s4 := 0.0; s5 := 0.0; s6 := 0.0; s7 := 0.0;
    for i in 1..npm loop
      s0 := s0 + pm0(i)*pm0(i);
      s1 := s1 + pm0(i)*pm1(i) + pm1(i)*pm0(i);
      s2 := s2 + pm0(i)*pm2(i) + pm1(i)*pm1(i) + pm2(i)*pm0(i);
      s3 := s3 + pm0(i)*pm3(i) + pm1(i)*pm2(i) + pm2(i)*pm1(i) + pm3(i)*pm0(i);
      s4 := s4 + pm1(i)*pm3(i) + pm2(i)*pm2(i) + pm3(i)*pm1(i);
      s5 := s5 + pm2(i)*pm3(i) + pm3(i)*pm2(i);
      s6 := s6 + pm3(i)*pm3(i);
    end loop;
    res := res + to_double_double(s0,s1,s2,s3,s4,s5,s6,s7,verbose);
    s0 := 0.0; s1 := 0.0; s2 := 0.0; s3 := 0.0;
    s4 := 0.0; s5 := 0.0; s6 := 0.0; s7 := 0.0;
    for i in 1..nmp loop
      s0 := s0 + mp0(i)*mp0(i);
      s1 := s1 + mp0(i)*mp1(i) + mp1(i)*mp0(i);
      s2 := s2 + mp0(i)*mp2(i) + mp1(i)*mp1(i) + mp2(i)*mp0(i);
      s3 := s3 + mp0(i)*mp3(i) + mp1(i)*mp2(i) + mp2(i)*mp1(i) + mp3(i)*mp0(i);
      s4 := s4 + mp1(i)*mp3(i) + mp2(i)*mp2(i) + mp3(i)*mp1(i);
      s5 := s5 + mp2(i)*mp3(i) + mp3(i)*mp2(i);
      s6 := s6 + mp3(i)*mp3(i);
    end loop;
    res := res + to_double_double(s0,s1,s2,s3,s4,s5,s6,s7,verbose);
    s0 := 0.0; s1 := 0.0; s2 := 0.0; s3 := 0.0;
    s4 := 0.0; s5 := 0.0; s6 := 0.0; s7 := 0.0;
    for i in 1..npm loop
      s4 := s4 + pm0(i)*pm4(i)
               + pm4(i)*pm0(i);
      s5 := s5 + pm0(i)*pm5(i) + pm1(i)*pm4(i)
               + pm4(i)*pm1(i) + pm5(i)*pm0(i);
      s6 := s6 + pm0(i)*pm6(i) + pm1(i)*pm5(i) + pm2(i)*pm4(i)
               + pm4(i)*pm2(i) + pm5(i)*pm1(i) + pm6(i)*pm0(i);
      s7 := s7 + pm0(i)*pm7(i) + pm1(i)*pm6(i) + pm2(i)*pm5(i) + pm3(i)*pm4(i)
               + pm4(i)*pm3(i) + pm5(i)*pm2(i) + pm6(i)*pm1(i) + pm7(i)*pm0(i);
    end loop;
    res := res + to_double_double(s0,s1,s2,s3,s4,s5,s6,s7,verbose);
    s0 := 0.0; s1 := 0.0; s2 := 0.0; s3 := 0.0;
    s4 := 0.0; s5 := 0.0; s6 := 0.0; s7 := 0.0;
    for i in 1..nmp loop
      s4 := s4 + mp0(i)*mp4(i)
               + mp4(i)*mp0(i);
      s5 := s5 + mp0(i)*mp5(i) + mp1(i)*mp4(i)
               + mp4(i)*mp1(i) + mp5(i)*mp0(i);
      s6 := s6 + mp0(i)*mp6(i) + mp1(i)*mp5(i) + mp2(i)*mp4(i)
               + mp4(i)*mp2(i) + mp5(i)*mp1(i) + mp6(i)*mp0(i);
      s7 := s7 + mp0(i)*mp7(i) + mp1(i)*mp6(i) + mp2(i)*mp5(i) + mp3(i)*mp4(i)
               + mp4(i)*mp3(i) + mp5(i)*mp2(i) + mp6(i)*mp1(i) + mp7(i)*mp0(i);
    end loop;
    res := res + to_double_double(s0,s1,s2,s3,s4,s5,s6,s7,verbose);
    return res;
  end Squared_Norm;

  function Squared_Norm
             ( x : DoblDobl_Complex_Vectors.Vector;
               verbose : boolean := true ) return double_double is

    res : double_double;
    xre,xim : Double_Double_Vectors.Vector(x'range);

  begin
    for i in x'range loop
      xre(i) := DoblDobl_Complex_Numbers.REAL_PART(x(i));
      xim(i) := DoblDobl_Complex_Numbers.IMAG_PART(x(i));
    end loop;
    res := Squared_Norm(xre,verbose) + Squared_Norm(xim,verbose);
    return res;
  end Squared_Norm;

  function Product ( x,y : Double_Double_Vectors.Vector;
                     verbose : boolean := true ) return double_double is

    res : double_double;
    dim : constant integer32 := x'length;
    xss0,xss1,xss2,xss3 : Standard_Floating_Vectors.Vector(1..dim);
    xss4,xss5,xss6,xss7 : Standard_Floating_Vectors.Vector(1..dim);
    yss0,yss1,yss2,yss3 : Standard_Floating_Vectors.Vector(1..dim);
    yss4,yss5,yss6,yss7 : Standard_Floating_Vectors.Vector(1..dim);
    xsd0,xsd1,xsd2,xsd3 : Standard_Floating_Vectors.Vector(1..dim);
    xsd4,xsd5,xsd6,xsd7 : Standard_Floating_Vectors.Vector(1..dim);
    ysd0,ysd1,ysd2,ysd3 : Standard_Floating_Vectors.Vector(1..dim);
    ysd4,ysd5,ysd6,ysd7 : Standard_Floating_Vectors.Vector(1..dim);
    xds0,xds1,xds2,xds3 : Standard_Floating_Vectors.Vector(1..dim);
    xds4,xds5,xds6,xds7 : Standard_Floating_Vectors.Vector(1..dim);
    yds0,yds1,yds2,yds3 : Standard_Floating_Vectors.Vector(1..dim);
    yds4,yds5,yds6,yds7 : Standard_Floating_Vectors.Vector(1..dim);
    xdd0,xdd1,xdd2,xdd3 : Standard_Floating_Vectors.Vector(1..dim);
    xdd4,xdd5,xdd6,xdd7 : Standard_Floating_Vectors.Vector(1..dim);
    ydd0,ydd1,ydd2,ydd3 : Standard_Floating_Vectors.Vector(1..dim);
    ydd4,ydd5,ydd6,ydd7 : Standard_Floating_Vectors.Vector(1..dim);
    nss,nsd,nds,ndd : integer32;
    s0,s1,s2,s3,s4,s5,s6,s7 : double_float;

  begin
    Signed_Quarter(x,y,
                   xss0,xss1,xss2,xss3,xss4,xss5,xss6,xss7,
                   yss0,yss1,yss2,yss3,yss4,yss5,yss6,yss7,
                   xsd0,xsd1,xsd2,xsd3,xsd4,xsd5,xsd6,xsd7,
                   ysd0,ysd1,ysd2,ysd3,ysd4,ysd5,ysd6,ysd7,
                   xds0,xds1,xds2,xds3,xds4,xds5,xds6,xds7,
                   yds0,yds1,yds2,yds3,yds4,yds5,yds6,yds7,
                   xdd0,xdd1,xdd2,xdd3,xdd4,xdd5,xdd6,xdd7,
                   ydd0,ydd1,ydd2,ydd3,ydd4,ydd5,ydd6,ydd7,
                   nss,nsd,nds,ndd,verbose);
    if verbose then
      put("#ss : "); put(nss,1); 
      put(", #sd : "); put(nsd,1); 
      put(", #ds : "); put(nds,1);
      put(", #dd : "); put(ndd,1); new_line;
    end if;
    s0 := 0.0; s1 := 0.0; s2 := 0.0; s3 := 0.0;
    s4 := 0.0; s5 := 0.0; s6 := 0.0; s7 := 0.0;
    for i in 1..nss loop
      s0 := s0 + xss0(i)*yss0(i);
      s1 := s1 + xss0(i)*yss1(i) + xss1(i)*yss0(i);
      s2 := s2 + xss0(i)*yss2(i) + xss1(i)*yss1(i) + xss2(i)*yss0(i);
      s3 := s3 + xss0(i)*yss3(i) + xss1(i)*yss2(i) + xss2(i)*yss1(i)
               + xss3(i)*yss0(i);
      s4 := s4 + xss0(i)*yss4(i) + xss1(i)*yss3(i) + xss2(i)*yss2(i)
               + xss3(i)*yss1(i) + xss4(i)*yss0(i);
      s5 := s5 + xss0(i)*yss5(i) + xss1(i)*yss4(i) + xss2(i)*yss3(i)
               + xss3(i)*yss2(i) + xss4(i)*yss1(i) + xss5(i)*yss0(i);
      s6 := s6 + xss0(i)*yss6(i) + xss1(i)*yss5(i) + xss2(i)*yss4(i)
               + xss3(i)*yss3(i) + xss4(i)*yss2(i) + xss5(i)*yss1(i)
               + xss6(i)*yss0(i);
      s7 := s7 + xss0(i)*yss7(i) + xss1(i)*yss6(i) + xss2(i)*yss5(i)
               + xss3(i)*yss4(i) + xss4(i)*yss3(i) + xss5(i)*yss2(i)
               + xss6(i)*yss1(i) + xss7(i)*yss0(i);
    end loop;
    res := to_double_double(s0,s1,s2,s3,s4,s5,s6,s7,verbose);
    s0 := 0.0; s1 := 0.0; s2 := 0.0; s3 := 0.0;
    s4 := 0.0; s5 := 0.0; s6 := 0.0; s7 := 0.0;
    for i in 1..nsd loop
      s0 := s0 + xsd0(i)*ysd0(i);
      s1 := s1 + xsd0(i)*ysd1(i) + xsd1(i)*ysd0(i);
      s2 := s2 + xsd0(i)*ysd2(i) + xsd1(i)*ysd1(i) + xsd2(i)*ysd0(i);
      s3 := s3 + xsd0(i)*ysd3(i) + xsd1(i)*ysd2(i) + xsd2(i)*ysd1(i)
               + xsd3(i)*ysd0(i);
      s4 := s4 + xsd1(i)*ysd3(i) + xsd2(i)*ysd2(i) + xsd3(i)*ysd1(i);
      s5 := s5 + xsd2(i)*ysd4(i) + xsd3(i)*ysd2(i);
      s6 := s6 + xds3(i)*ysd5(i);
    end loop;
    res := res + to_double_double(s0,s1,s2,s3,s4,s5,s6,s7,verbose);
    s0 := 0.0; s1 := 0.0; s2 := 0.0; s3 := 0.0;
    s4 := 0.0; s5 := 0.0; s6 := 0.0; s7 := 0.0;
    for i in 1..ndd loop
      s0 := s0 + xss0(i)*yss0(i);
      s1 := s1 + xss0(i)*yss1(i) + xss1(i)*yss0(i);
      s2 := s2 + xss0(i)*yss2(i) + xss1(i)*yss1(i) + xss2(i)*yss0(i);
      s3 := s3 + xss0(i)*yss3(i) + xss1(i)*yss2(i) + xss2(i)*yss1(i)
               + xss3(i)*yss0(i);
      s4 := s4 + xss0(i)*yss4(i) + xss1(i)*yss3(i) + xss2(i)*yss2(i)
               + xss3(i)*yss1(i) + xss4(i)*yss0(i);
      s5 := s5 + xss0(i)*yss5(i) + xss1(i)*yss4(i) + xss2(i)*yss3(i)
               + xss3(i)*yss2(i) + xss4(i)*yss1(i) + xss5(i)*yss0(i);
      s6 := s6 + xss0(i)*yss6(i) + xss1(i)*yss5(i) + xss2(i)*yss4(i)
               + xss3(i)*yss3(i) + xss4(i)*yss2(i) + xss5(i)*yss1(i)
               + xss6(i)*yss0(i);
      s7 := s7 + xss0(i)*yss7(i) + xss1(i)*yss6(i) + xss2(i)*yss5(i)
               + xss3(i)*yss4(i) + xss4(i)*yss3(i) + xss5(i)*yss2(i)
               + xss6(i)*yss1(i) + xss7(i)*yss0(i);
    end loop;
    res := res + to_double_double(s0,s1,s2,s3,s4,s5,s6,s7,verbose);
    return res;
  end Product;

end Vectored_Double_Doubles;
