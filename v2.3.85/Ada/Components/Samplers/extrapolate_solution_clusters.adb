with Standard_Integer_Numbers_io;      use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;        use Standard_Floating_Numbers;
with Standard_Mathematical_Functions;  use Standard_Mathematical_Functions;
with Multprec_Floating_Numbers;        use Multprec_Floating_Numbers;
with Multprec_Floating_Numbers_io;     use Multprec_Floating_Numbers_io;
with Multprec_Complex_Numbers;         use Multprec_Complex_Numbers;
with Multprec_Complex_Number_Tools;    use Multprec_Complex_Number_Tools;
with Multprec_Complex_Numbers_io;      use Multprec_Complex_Numbers_io;
with Multprec_Complex_Vectors;         use Multprec_Complex_Vectors;
with Multprec_Complex_Vectors_io;      use Multprec_Complex_Vectors_io;

package body Extrapolate_Solution_Clusters is

  function Extrapolate ( x1,x0,s1m,s0m,s10m : Complex_Number )
                       return Complex_Number is

  -- DESCRIPTION :
  --   Returns a first order approximation using x0 = x(s0), x1 = x(s0),
  --   and with s1m = s1^m, s0m = s0^m, and s10m = s1m - s0m.

    res,acc : Complex_Number;

  begin
    res := s1m*x0;
    acc := s0m*x1;
    Sub(res,acc);
    Clear(acc);
    Div(res,s10m);
    return res;
  end Extrapolate;

  function Extrapolate ( v1,v0 : Vector; s1m,s0m,s10m : Complex_Number )
                       return Vector is

  -- DESCRIPTION :
  --   Applies the scalar extrapolator above componentwise to the vectors.

    res : Vector(v1'range);

  begin
    for i in res'range loop
      res(i) := Extrapolate(v1(i),v0(i),s1m,s0m,s10m);
    end loop;
    return res;
  end Extrapolate;

  function Extrapolate ( x1,x0,r,onemr : Complex_Number )
                       return Complex_Number is

  -- DESCRIPTION :
  --   Returns (x0 - r*x1)/onemr.

    res : Complex_Number := r*x1;

  begin
    Min(res);
    Add(res,x0);
    Div(res,onemr);
    return res;
  end Extrapolate;

  function Extrapolate ( v1,v0 : Vector; r,onemr : Complex_Number )
                       return Vector is

    res : Vector(v1'range);

  begin
    for i in res'range loop
      res(i) := Extrapolate(v1(i),v0(1),r,onemr);
    end loop;
    return res;
  end Extrapolate;

  procedure Extrapolate
              ( file : in file_type; order,m,n : in integer32;
                size : in natural32;
                sols : in Solution_Array; extra : out Vector ) is

    s,sm : Multprec_Complex_Vectors.Vector(0..order);
    difs10,difs20,f10,f20,acc,ratio,onemratio : Complex_Number;
    extra10,extra20,extra210 : Vector(1..n);
    mpfs : Floating_Number;
    stfs,expo : double_float;

  begin
    for i in sols'range loop
      Copy(sols(i).v(n),s(i));
      put(file,"s("); put(file,i,1); put(file,") : ");
      put(file,s(i)); new_line(file);
      mpfs := REAL_PART(s(i));
      stfs := Round(mpfs);
      expo := 1.0/double_float(m);
      stfs := stfs**expo;
      Clear(mpfs);
      mpfs := Create(stfs);
      sm(i) := Create(mpfs);
      Set_Size(sm(i),size);
      put(file,"s("); put(file,i,1); put(file,")^1/");
      put(file,m,1); put(file," : "); put(file,sm(i));
      new_line(file);
    end loop;
    difs10 := sm(1) - sm(0);
    extra10 := Extrapolate(sols(1).v,sols(0).v,sm(1),sm(0),difs10);
    difs20 := sm(2) - sm(0);
    extra20 := Extrapolate(sols(2).v,sols(0).v,sm(2),sm(0),difs20);
    f10 := sm(1)*s(0);
    acc := sm(0)*s(1);
    Sub(f10,acc); Clear(acc);
    Div(f10,difs10);
    f20 := sm(2)*s(0);
    acc := sm(0)*s(2);
    Sub(f20,acc); Clear(acc);
    Div(f20,difs20);
    acc := f20 - f10;
    extra210 := Extrapolate(extra20,extra10,f20,f10,acc);
   -- ratio := sm(order)/sm(order-1);
   -- onemratio := Create(1.0) - ratio;
   -- extrasolv := Extrapolate(sols(order-1).v,sols(order).v,
   --                          ratio,onemratio);
    put_line(file,"The extrapolated solution vector : ");
    put_line(file,extra210);
    extra := extra210;
  end Extrapolate;

end Extrapolate_Solution_Clusters;
