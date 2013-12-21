with Standard_Natural_Numbers_io;       use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers_io;      use Standard_Floating_Numbers_io;
with Standard_Complex_Numbers_io;       use Standard_Complex_Numbers_io;
with Standard_Complex_Numbers_Polar;    use Standard_Complex_Numbers_Polar;
with Standard_Random_Numbers;           use Standard_Random_Numbers;
with Standard_Complex_Norms_Equals;     use Standard_Complex_Norms_Equals;

package body Hypersurface_Witness_Sets is

-- UTILITIES needed in method of Durand-Kerner

  procedure Divided_Differences ( x : in Vector; f : in out Vector ) is

  -- DESCRIPTION :
  --   Computes in f the divided differences using the points in x.

  begin
   -- put_line("Entered Divided_Differences...");
   -- put_line(" x = ");
   -- for i in x'range loop
   --   put(x(i)); new_line;
   -- end loop;
   -- put_line(" f = ");
   -- for i in x'range loop
   --   put(f(i)); new_line;
   -- end loop;
    for i in 1..f'last loop
      for j in 0..i-1 loop
       -- put("x("); put(i,1); put(")-x("); put(j,1); put(") :");
       -- put(x(i)-x(j)); new_line;
        f(i) := (f(i)-f(j))/(x(i)-x(j));
      end loop;
    end loop;
  end Divided_Differences;

  function Roots_of_Unity ( d : natural32 ) return Vector is

  -- DESCRIPTION :
  --   Returns a vector with the d complex roots of unity.
  --   We take roots as unity as start values for the univariate
  --   root finder.

    res : Vector(1..integer32(d));
    one : Complex_Number := Create(1.0);

  begin
    for i in 1..d loop
      res(integer32(i)) := Root(one,d,i);
    end loop;
    return res;
  end Roots_of_Unity;

  function Compute_q ( i : integer32; a : Vector ) return Complex_Number is

  -- DESCRIPTION :
  --   Computes the quotient needed in the Durand-Kerner step.

    res : Complex_Number;

  begin
    res := Create(1.0);
    for j in a'range loop
      if j /= i
       then res := res*(a(i)-a(j));
      end if;
    end loop;
    return res;
  end Compute_q;

  procedure Write ( file : in file_type; z,err,res : in Vector ) is

  -- DESCRIPTION :
  --   Writes the current approximations for the roots in z,
  --   along with their residuals in res on the file.

    f : double_float;

  begin
    for i in z'range loop
      put(file,i,2); put(file," : ");
      put(file,z(i));
      put(file," : "); f := AbsVal(err(i)); put(file,f,3);
      put(file," : "); f := AbsVal(res(i)); put(file,f,3);
      new_line(file);
    end loop;
  end Write;

-- TARGET PROCEDURES :

  procedure Silent_Root_Finder0
              ( degree : in natural32; eps : in double_float;
                max_it : in natural32; fail : out boolean;
                b,v : in Vector; t : out Vector; nrm : out double_float ) is

    err,res : Vector(1..integer32(degree));

    procedure Find_Roots is new Silent_Root_Finder1(f);

  begin
    Find_Roots(degree,eps,max_it,fail,b,v,t,err,res,nrm);
  end Silent_Root_Finder0;

  procedure Silent_Root_Finder1
              ( degree : in natural32; eps : in double_float;
                max_it : in natural32; fail : out boolean;
                b,v : in Vector; t,err,res : out Vector;
                nrm : out double_float ) is

  -- IMPLEMENTATION NOTE :
  --   The iteration stops before the maximal allowed number
  --   if the desired accuracy is reached and the method is
  --   no longer converging.

    previous_nrm : double_float;
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
    pts(0) := Create(0.0);
    dvd(0) := Eval(b,v,pts(0));
    for i in 1..dvd'last loop
      pts(i) := z(i);
      dvd(i) := Eval(b,v,pts(i));
    end loop;
    Divided_Differences(pts,dvd);
    lc := dvd(dvd'last);
    previous_nrm := 1.0E+8;
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
                degree : in natural32; eps : in double_float;
                max_it : in natural32; fail : out boolean;
                b,v : in Vector; t : out Vector; nrm : out double_float ) is

    err,res : Vector(1..integer32(degree));

    procedure Find_Roots is new Reporting_Root_Finder1(f);

  begin
    Find_Roots(file,degree,eps,max_it,fail,b,v,t,err,res,nrm);
  end Reporting_Root_Finder0;

  procedure Reporting_Root_Finder1
              ( file : in file_type;
                degree : in natural32; eps : in double_float;
                max_it : in natural32; fail : out boolean;
                b,v : in Vector; t,err,res : out Vector;
                nrm : out double_float ) is

  -- IMPLEMENTATION NOTE :
  --   The iteration stops before the maximal allowed number
  --   if the desired accuracy is reached and the method is
  --   no longer converging.

    previous_nrm : double_float;
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
    pts(0) := Create(0.0);
    dvd(0) := Eval(b,v,pts(0));
    for i in z'range loop
      pts(i) := z(i);
      dvd(i) := Eval(b,v,pts(i));
    end loop;
   -- put_line("Calling Divided_Differences...");
    Divided_Differences(pts,dvd);
   -- put_line("... done with Divided_Differences.");
    lc := dvd(dvd'last);
    previous_nrm := 1.0E+8;
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

end Hypersurface_Witness_Sets;
