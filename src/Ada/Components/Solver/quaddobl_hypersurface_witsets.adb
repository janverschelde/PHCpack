with Standard_Natural_Numbers_io;       use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with Quad_Double_Numbers_io;            use Quad_Double_Numbers_io;
with QuadDobl_Complex_Numbers_io;       use QuadDobl_Complex_Numbers_io;
with QuadDobl_Complex_Numbers_Polar;    use QuadDobl_Complex_Numbers_Polar;
with QuadDobl_Random_Numbers;           use QuadDobl_Random_Numbers;
with QuadDobl_Random_Vectors;           use QuadDobl_Random_Vectors;
with QuadDobl_Complex_Vector_Norms;     use QuadDobl_Complex_Vector_Norms;

package body QuadDobl_Hypersurface_Witsets is

-- PART I : the primitives in the Durand-Kerner method (aka Weierstrass)

  procedure Divided_Differences ( x : in Vector; f : in out Vector ) is
  begin
    for i in 1..f'last loop
      for j in 0..i-1 loop
        f(i) := (f(i)-f(j))/(x(i)-x(j));
      end loop;
    end loop;
  end Divided_Differences;

  function Roots_of_Unity ( d : natural32 ) return Vector is

    res : Vector(1..integer32(d));
    one : Complex_Number := Create(integer(1));

  begin
    for i in 1..d loop
      res(integer32(i)) := Root(one,d,i);
    end loop;
    return res;
  end Roots_of_Unity;

  function Compute_q ( i : integer32; a : Vector ) return Complex_Number is

    res : Complex_Number;

  begin
    res := Create(integer(1));
    for j in a'range loop
      if j /= i
       then res := res*(a(i)-a(j));
      end if;
    end loop;
    return res;
  end Compute_q;

  procedure Write ( file : in file_type; z,err,res : in Vector ) is

    f : quad_double;

  begin
    for i in z'range loop
      put(file,i,2); put(file," : ");
      put(file,z(i));
      put(file," : "); f := AbsVal(err(i)); put(file,f,3);
      put(file," : "); f := AbsVal(res(i)); put(file,f,3);
      new_line(file);
    end loop;
  end Write;

-- PART II : generic procedures for the root finders

  procedure Silent_Root_Finder0
              ( degree : in natural32; eps : in quad_double;
                max_it : in natural32; fail : out boolean;
                b,v : in Vector; t : out Vector; nrm : out quad_double ) is

    err,res : Vector(1..integer32(degree));

    procedure Find_Roots is new Silent_Root_Finder1(f);

  begin
    Find_Roots(degree,eps,max_it,fail,b,v,t,err,res,nrm);
  end Silent_Root_Finder0;

  procedure Silent_Root_Finder1
              ( degree : in natural32; eps : in quad_double;
                max_it : in natural32; fail : out boolean;
                b,v : in Vector; t,err,res : out Vector;
                nrm : out quad_double ) is

  -- IMPLEMENTATION NOTE :
  --   The iteration stops before the maximal allowed number
  --   if the desired accuracy is reached and the method is
  --   no longer converging.

    previous_nrm : quad_double;
    z : Vector(1..integer32(degree)) := Roots_of_Unity(degree);
    pts,dvd : Vector(0..integer32(degree));
    lc : Complex_Number;

    function Eval ( b,v : Vector; t : Complex_Number )
                  return Complex_Number is

    -- DESCRIPTION :
    --   Returns the function value of p at the point b + t*v.

      point : Vector(b'range) := b + t*v;

    begin
      return f(point);
    end Eval;

    procedure DK ( b,v : in Vector; lc : in Complex_Number;
                   z,err,res : in out Vector ) is

    -- DESCRIPTION :
    --   Computes one step in the Durand-Kerner iteration with
    --   the polynomial p/lc, where lc is the leading coefficient
    --   for the parameteric representation using the direction v.
    --   For high degrees, overflow exceptions often occur.
    --   To handle this, we just choose another random starting point.

    begin
      for i in z'range loop
        declare
        begin
          err(i) := - Eval(b,v,z(i))/(lc*Compute_q(i,z));
          z(i) := z(i) + err(i);
          res(i) := Eval(b,v,z(i))/lc;
        exception
          when others => --put("Exception occurred at component ");
                         --put(i,1); put_line(".");
                         z(i) := Random1;
        end;
      end loop;
    end DK;

  begin
    pts(0) := Create(integer(0));
    dvd(0) := Eval(b,v,pts(0));
    for i in 1..dvd'last loop
      pts(i) := z(i);
      dvd(i) := Eval(b,v,pts(i));
    end loop;
    Divided_Differences(pts,dvd);
    lc := dvd(dvd'last);
    previous_nrm := create(1.0E+8);
    fail := true;
    for k in 1..max_it loop
      DK(b,v,lc,z,err,res);
      nrm := Max_Norm(res);
      exit when (nrm <= eps) and (previous_nrm <= 10.0*nrm);
    end loop;
    fail := (nrm > eps);
    t := z;
  end Silent_Root_Finder1;

  procedure Reporting_Root_Finder0
              ( file : in file_type;
                degree : in natural32; eps : in quad_double;
                max_it : in natural32; fail : out boolean;
                b,v : in Vector; t : out Vector; nrm : out quad_double ) is

    err,res : Vector(1..integer32(degree));

    procedure Find_Roots is new Reporting_Root_Finder1(f);

  begin
    Find_Roots(file,degree,eps,max_it,fail,b,v,t,err,res,nrm);
  end Reporting_Root_Finder0;

  procedure Reporting_Root_Finder1
              ( file : in file_type;
                degree : in natural32; eps : in quad_double;
                max_it : in natural32; fail : out boolean;
                b,v : in Vector; t,err,res : out Vector;
                nrm : out quad_double ) is

  -- IMPLEMENTATION NOTE :
  --   The iteration stops before the maximal allowed number
  --   if the desired accuracy is reached and the method is
  --   no longer converging.

    previous_nrm : quad_double;
    z : Vector(1..integer32(degree)) := Roots_of_Unity(degree);
    pts,dvd : Vector(0..integer32(degree));
    lc : Complex_Number;

    function Eval ( b,v : Vector; t : Complex_Number )
                  return Complex_Number is

    -- DESCRIPTION :
    --   Returns the function value of p at the point b + t*v.

      point : Vector(b'range) := b + t*v;

    begin
      return f(point);
    end Eval;

    procedure DK ( b,v : in Vector; lc : in Complex_Number;
                   z,err,res : in out Vector ) is

    -- DESCRIPTION :
    --   Computes one step in the Durand-Kerner iteration with
    --   the polynomial p/lc, where lc is the leading coefficient
    --   for the parameteric representation using the direction v.
    --   For high degrees, overflow exceptions often occur.
    --   To handle this, we just choose another random starting point.

    begin
      for i in z'range loop
        declare
        begin
          err(i) := - Eval(b,v,z(i))/(lc*Compute_q(i,z));
          z(i) := z(i) + err(i);
          res(i) := Eval(b,v,z(i))/lc;
        exception
          when others => --put("Exception occurred at component ");
                         --put(i,1); put_line(".");
                         z(i) := Random1;
        end;
      end loop;
    end DK;

  begin
   -- put_line("Entered Reporting_Root_Finder1...");
    pts(0) := Create(integer(0));
    dvd(0) := Eval(b,v,pts(0));
    for i in z'range loop
      pts(i) := z(i);
      dvd(i) := Eval(b,v,pts(i));
    end loop;
   -- put_line("Calling Divided_Differences...");
    Divided_Differences(pts,dvd);
   -- put_line("... done with Divided_Differences.");
    lc := dvd(dvd'last);
    previous_nrm := create(1.0E+8);
    fail := true;
    for k in 1..max_it loop
      DK(b,v,lc,z,err,res);
      put(file,"The results after step ");
      put(file,k,1); put_line(file," : ");
      Write(file,z,err,res);
      nrm := Max_Norm(res);
      put(file,"The residual norm : ");
      put(file,nrm,3); new_line(file);
      exit when (nrm <= eps) and (previous_nrm <= 10.0*nrm);
      previous_nrm := nrm;
    end loop;
    fail := (nrm > eps);
    t := z;
  end Reporting_Root_Finder1;

end QuadDobl_Hypersurface_Witsets;
