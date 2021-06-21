with text_io;                           use text_io;
with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with Standard_Complex_Numbers;          use Standard_Complex_Numbers;
with Standard_Complex_Numbers_io;       use Standard_Complex_Numbers_io;

package body Double_Laurent_Series is

  procedure Write ( e : in integer32;
                    c : in Standard_Complex_Vectors.Vector ) is
  begin
    Write(standard_output,e,c);
  end Write;

  procedure Write ( file : in file_type; e : in integer32;
                    c : in Standard_Complex_Vectors.Vector ) is
  begin
    for i in c'range loop
      if i > c'first
       then put(file," + (");
       else put(file,"   (");
      end if;
      put(file,c(i)); put(file,")*t^"); put(file,e+i,1);
      new_line(file);
    end loop;
  end Write;

  procedure Multiply ( d,xe,ye : in integer32;
                       xc,yc : in Standard_Complex_Vectors.Vector;
                       ze : out integer32;
                       zc : out Standard_Complex_Vectors.Vector ) is
  begin
    ze := xe + ye;
    for i in 0..d loop
       zc(i) := xc(0)*yc(i);
       for j in 1..i loop
         zc(i) := zc(i) + xc(j)*yc(i-j);
       end loop;
    end loop;
  end Multiply;

  procedure Inverse ( d,xe : in integer32;
                      xc : in Standard_Complex_Vectors.Vector;
                      ye : out integer32;
                      yc : out Standard_Complex_Vectors.Vector ) is
  begin
    ye := -xe;
    yc(0) := 1.0/xc(0);
    for i in 1..d loop
      yc(i) := -xc(1)*yc(i-1);
      for j in 2..i loop
        yc(i) := yc(i) - xc(j)*yc(i-j);
      end loop;
      yc(i) := yc(i)/xc(0);
    end loop;
  end Inverse;

  procedure Divide ( d,xe,ye : in integer32;
                     xc,yc : in Standard_Complex_Vectors.Vector;
                     ze : out integer32;
                     zc : out Standard_Complex_Vectors.Vector;
                     iyc : out Standard_Complex_Vectors.Vector ) is

    iye : integer32;

  begin
    Inverse(d,ye,yc,iye,iyc);
    Multiply(d,xe,iye,xc,iyc,ze,zc);
  end Divide;

  function Is_Zero ( d : integer32;
                     c : Standard_Complex_Vectors.Vector;
                     tol : double_float := 1.0E-15 ) return boolean is

  begin
    for i in 0..d loop
      if (AbsVal(c(0)) > tol)
       then return false;
      end if;
    end loop;
    return true;
  end Is_Zero;

  procedure Normalize ( d : in integer32; e : in out integer32;
                        c : in out Standard_Complex_Vectors.Vector;
                        tol : in double_float := 1.0E-15 ) is

    allzero : boolean := true;

  begin
    for i in 0..d loop
      if (AbsVal(c(0)) > tol)
       then allzero := false; exit;
      end if;
      e := e + 1;
      for k in 1..d-i loop -- shift the coefficients
        c(k-1) := c(k);
      end loop;
      for k in (d-i)..d loop -- must insert zeros at the end
        c(k) := Standard_Complex_Numbers.Create(0.0);
      end loop;
    end loop;
    if allzero
     then e := 0;
    end if;
  end Normalize;

  function Exponent_Gap ( a,b : integer32 ) return integer32 is

    result : integer32 := 0;

  begin
    if a <= b then
      -- if a >= 0 then -- 0 <= a <= b
      --   result := b - a;
      -- else -- a < 0
      --   if b <= 0 then -- a <= b <= 0
      --     result := b - a; -- a is more negative than b, so result >= 0
      --   else -- a < 0 and b > 0
      --     result := b - a;
      --   end if;
      -- end if;
      result := b - a; -- in all the above case, the same formula
    else -- a > b
      result := a - b; -- by symmetry with above case analysis
    end if;
    return result;
  end Exponent_Gap;

  procedure Add ( d,xe,ye : in integer32;
                  xc,yc : in Standard_Complex_Vectors.Vector;
                  ze : out integer32;
                  zc : out Standard_Complex_Vectors.Vector;
                  tol : in double_float := 1.0E-15 ) is

    gap : integer32; -- gap between the leading coefficients

  begin
    if Is_Zero(d,xc,tol) then -- z is a copy of y
      ze := ye;
      for i in 0..d loop
        zc(i) := yc(i);
      end loop;
    elsif Is_Zero(d,yc,tol) then -- z is a copy of x
      ze := xe;
      for i in 0..d loop
        zc(i) := xc(i);
      end loop;
    else
      if xe < ye then
        ze := xe;
        gap := Exponent_Gap(xe,ye);
        for i in 0..gap-1 loop
          exit when (i > zc'last);
          zc(i) := xc(i);
        end loop;
        for i in gap..d loop
          zc(i) := xc(i) + yc(i-gap);
        end loop;
      elsif xe > ye then
        ze := ye;
        gap := Exponent_Gap(xe,ye);
        for i in 0..gap-1 loop
          exit when (i > zc'last);
          zc(i) := yc(i);
        end loop;
        for i in gap..d loop
          zc(i) := yc(i) + xc(i-gap);
        end loop;
      else -- xe = ye
        ze := xe;
        for i in 0..d loop
          zc(i) := xc(i) + yc(i);
        end loop;
        Normalize(d,ze,zc,tol);
      end if;
    end if;
  end Add;

  procedure Subtract ( d,xe,ye : in integer32;
                       xc,yc : in Standard_Complex_Vectors.Vector;
                       ze : out integer32;
                       zc : out Standard_Complex_Vectors.Vector;
                       tol : in double_float := 1.0E-15 ) is

    gap : integer32; -- gap between the leading coefficients

  begin
    if Is_Zero(d,xc,tol) then -- z equals -y
      ze := ye;
      for i in 0..d loop
        zc(i) := -yc(i);
      end loop;
    elsif Is_Zero(d,yc,tol) then -- z is a copy of x
      ze := xe;
      for i in 0..d loop
        zc(i) := xc(i);
      end loop;
    else
      if xe < ye then
        ze := xe;
        gap := Exponent_Gap(xe,ye);
        for i in 0..gap-1 loop
          exit when (i > zc'last);
          zc(i) := xc(i);
        end loop;
        for i in gap..d loop
          zc(i) := xc(i) - yc(i-gap);
        end loop;
      elsif xe > ye then
        ze := ye;
        gap := Exponent_Gap(xe,ye);
        for i in 0..gap-1 loop
          exit when (i > zc'last);
          zc(i) := -yc(i);
        end loop;
        for i in gap..d loop
          zc(i) := xc(i-gap) - yc(i);
        end loop;
      else -- xe = ye
        ze := xe;
        for i in 0..d loop
          zc(i) := xc(i) - yc(i);
        end loop;
        Normalize(d,ze,zc,tol);
      end if;
    end if;
  end Subtract;

end Double_Laurent_Series;
