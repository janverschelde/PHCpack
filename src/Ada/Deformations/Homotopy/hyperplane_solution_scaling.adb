with Standard_Complex_Numbers_Polar;
with DoblDobl_Complex_Numbers_Polar;
with QuadDobl_Complex_Numbers_Polar;
with Standard_Complex_Vector_Norms;
with DoblDobl_Complex_Vector_Norms;
with QuadDobl_Complex_Vector_Norms;

package body Hyperplane_Solution_Scaling is

  procedure Sub ( p : in out Standard_Complex_Polynomials.Poly;
                  c : in Standard_Complex_Numbers.Complex_Number ) is

    use Standard_Complex_Polynomials;

    t : Term;
    n : constant integer32 := integer32(Number_of_Unknowns(p));
    d : constant Standard_Natural_Vectors.Vector(1..n) := (1..n => 0);

  begin
    t.cf := c;
    t.dg := new Standard_Natural_Vectors.Vector'(d);
    Sub(p,t);
    Clear(t);
  end Sub;

  procedure Sub ( p : in out DoblDobl_Complex_Polynomials.Poly;
                  c : in DoblDobl_Complex_Numbers.Complex_Number ) is

    use DoblDobl_Complex_Polynomials;

    t : Term;
    n : constant integer32 := integer32(Number_of_Unknowns(p));
    d : constant Standard_Natural_Vectors.Vector(1..n) := (1..n => 0);

  begin
    t.cf := c;
    t.dg := new Standard_Natural_Vectors.Vector'(d);
    Sub(p,t);
    Clear(t);
  end Sub;

  procedure Sub ( p : in out QuadDobl_Complex_Polynomials.Poly;
                  c : in QuadDobl_Complex_Numbers.Complex_Number ) is

    use QuadDobl_Complex_Polynomials;

    t : Term;
    n : constant integer32 := integer32(Number_of_Unknowns(p));
    d : constant Standard_Natural_Vectors.Vector(1..n) := (1..n => 0);

  begin
    t.cf := c;
    t.dg := new Standard_Natural_Vectors.Vector'(d);
    Sub(p,t);
    Clear(t);
  end Sub;

  procedure Scale ( v : in out Standard_Complex_Vectors.Vector ) is

    mxv : constant double_float := Standard_Complex_Vector_Norms.Max_Norm(v);

    use Standard_Complex_Numbers;

  begin
    for i in v'range loop
      v(i) := v(i)/mxv;
    end loop;
  end Scale;

  procedure Scale ( v : in out DoblDobl_Complex_Vectors.Vector ) is

    mxv : constant double_double
        := DoblDobl_Complex_Vector_Norms.Max_Norm(v);

    use DoblDobl_Complex_Numbers;

  begin
    for i in v'range loop
      v(i) := v(i)/mxv;
    end loop;
  end Scale;

  procedure Scale ( v : in out QuadDobl_Complex_Vectors.Vector ) is

    mxv : constant quad_double
        := QuadDobl_Complex_Vector_Norms.Max_Norm(v);

    use QuadDobl_Complex_Numbers;

  begin
    for i in v'range loop
      v(i) := v(i)/mxv;
    end loop;
  end Scale;

  function Max_Norm ( v : Standard_Complex_Vectors.Vector;
                      k : natural32;
                      z : Standard_Natural_Vectors.Vector ) 
                    return double_float is

    res : double_float := -1.0;
    val : double_float;
    idx : constant integer32 := z'last + integer32(k);

  begin
    for i in z'range loop
      if z(i) = k then
        val := Standard_Complex_Numbers.AbsVal(v(i));
        if val > res
         then res := val;
        end if;
      end if;
    end loop;
    val := Standard_Complex_Numbers.AbsVal(v(idx));
    if val > res
     then res := val;
    end if;
    return res;
  end Max_Norm;

  function Max_Norm ( v : DoblDobl_Complex_Vectors.Vector;
                      k : natural32;
                      z : Standard_Natural_Vectors.Vector ) 
                    return double_double is

    res : double_double := create(-1.0);
    val : double_double;
    idx : constant integer32 := z'last + integer32(k);

  begin
    for i in z'range loop
      if z(i) = k then
        val := DoblDobl_Complex_Numbers.AbsVal(v(i));
        if val > res
         then res := val;
        end if;
      end if;
    end loop;
    val := DoblDobl_Complex_Numbers.AbsVal(v(idx));
    if val > res
     then res := val;
    end if;
    return res;
  end Max_Norm;

  function Max_Norm ( v : QuadDobl_Complex_Vectors.Vector;
                      k : natural32;
                      z : Standard_Natural_Vectors.Vector ) 
                    return quad_double is

    res : quad_double := create(-1.0);
    val : quad_double;
    idx : constant integer32 := z'last + integer32(k);

  begin
    for i in z'range loop
      if z(i) = k then
        val := QuadDobl_Complex_Numbers.AbsVal(v(i));
        if val > res
         then res := val;
        end if;
      end if;
    end loop;
    val := QuadDobl_Complex_Numbers.AbsVal(v(idx));
    if val > res
     then res := val;
    end if;
    return res;
  end Max_Norm;

  procedure Scale ( v : in out Standard_Complex_Vectors.Vector;
                    m : in natural32;
                    z : in Standard_Natural_Vectors.Vector ) is

    mxv : double_float;
    idx : integer32;

    use Standard_Complex_Numbers;

  begin
    for k in 1..m loop
      mxv := Max_Norm(v,k,z);
      for i in z'range loop
        if z(i) = k
         then v(i) := v(i)/mxv;
        end if;
      end loop;
      idx := z'last + integer32(k);
      v(idx) := v(idx)/mxv;
    end loop;
  end Scale;

  procedure Scale ( v : in out DoblDobl_Complex_Vectors.Vector;
                    m : in natural32;
                    z : in Standard_Natural_Vectors.Vector ) is

    mxv : double_double;
    idx : integer32;

    use DoblDobl_Complex_Numbers;

  begin
    for k in 1..m loop
      mxv := Max_Norm(v,k,z);
      for i in z'range loop
        if z(i) = k
         then v(i) := v(i)/mxv;
        end if;
      end loop;
      idx := z'last + integer32(k);
      v(idx) := v(idx)/mxv;
    end loop;
  end Scale;

  procedure Scale ( v : in out QuadDobl_Complex_Vectors.Vector;
                    m : in natural32;
                    z : in Standard_Natural_Vectors.Vector ) is

    mxv : quad_double;
    idx : integer32;

    use QuadDobl_Complex_Numbers;

  begin
    for k in 1..m loop
      mxv := Max_Norm(v,k,z);
      for i in z'range loop
        if z(i) = k
         then v(i) := v(i)/mxv;
        end if;
      end loop;
      idx := z'last + integer32(k);
      v(idx) := v(idx)/mxv;
    end loop;
  end Scale;

  procedure Adjust ( c : in Standard_Complex_Vectors.Link_to_Vector;
                     v : in Standard_Complex_Vectors.Vector ) is

    use Standard_Complex_Numbers;

    val : Complex_Number := c(c'last);

  begin
    for i in v'range loop
      val := val + c(i)*v(i);
    end loop;
    c(c'last) := c(c'last) - val;
  end Adjust;

  procedure Adjust ( c : in DoblDobl_Complex_Vectors.Link_to_Vector;
                     v : in DoblDobl_Complex_Vectors.Vector ) is

    use DoblDobl_Complex_Numbers;

    val : Complex_Number := c(c'last);

  begin
    for i in v'range loop
      val := val + c(i)*v(i);
    end loop;
    c(c'last) := c(c'last) - val;
  end Adjust;

  procedure Adjust ( c : in QuadDobl_Complex_Vectors.Link_to_Vector;
                     v : in QuadDobl_Complex_Vectors.Vector ) is

    use QuadDobl_Complex_Numbers;

    val : Complex_Number := c(c'last);

  begin
    for i in v'range loop
      val := val + c(i)*v(i);
    end loop;
    c(c'last) := c(c'last) - val;
  end Adjust;

-- MULTI-HOMOGENEOUS SOLUTION SCALING :

  procedure Scale ( v : in out Standard_Complex_Vectors.Vector;
                    idz : in Standard_Natural_Vectors.Link_to_Vector;
                    m,i : in integer32 ) is

    use Standard_Complex_Numbers;

    dim : constant integer32 := v'last - m;
    first : boolean := true;
    maxrad,radval : double_float;
    idx : constant integer32 := dim + i; -- index of i-th homogeneous variable

  begin
    for k in v'first..dim loop
      if integer32(idz(k)) = i then -- select variable from the i-th set
        if first then
          maxrad := Standard_Complex_Numbers_Polar.Radius(v(k));
          first := false;
        else
          radval := Standard_Complex_Numbers_Polar.Radius(v(k));
          if radval > maxrad
           then maxrad := radval;
          end if;
        end if;
      end if;
    end loop;
    radval := Standard_Complex_Numbers_Polar.Radius(v(idx));
    if radval > maxrad
     then maxrad := radval;
    end if;
    for k in v'first..dim loop
      if integer32(idz(k)) = i then -- select variable from the i-th set
        v(k) := v(k)/maxrad;
      end if;
    end loop;
    v(idx) := v(idx)/maxrad;
  end Scale;

  procedure Scale ( v : in out DoblDobl_Complex_Vectors.Vector;
                    idz : in Standard_Natural_Vectors.Link_to_Vector;
                    m,i : in integer32 ) is

    use DoblDobl_Complex_Numbers;

    dim : constant integer32 := v'last - m;
    first : boolean := true;
    maxrad,radval : double_double;
    idx : constant integer32 := dim + i; -- index of i-th homogeneous variable

  begin
    for k in v'first..dim loop
      if integer32(idz(k)) = i then -- select variable from the i-th set
        if first then
          maxrad := DoblDobl_Complex_Numbers_Polar.Radius(v(k));
          first := false;
        else
          radval := DoblDobl_Complex_Numbers_Polar.Radius(v(k));
          if radval > maxrad
           then maxrad := radval;
          end if;
        end if;
      end if;
    end loop;
    radval := DoblDobl_Complex_Numbers_Polar.Radius(v(idx));
    if radval > maxrad
     then maxrad := radval;
    end if;
    for k in v'first..dim loop
      if integer32(idz(k)) = i then -- select variable from the i-th set
        v(k) := v(k)/maxrad;
      end if;
    end loop;
    v(idx) := v(idx)/maxrad;
  end Scale;

  procedure Scale ( v : in out QuadDobl_Complex_Vectors.Vector;
                    idz : in Standard_Natural_Vectors.Link_to_Vector;
                    m,i : in integer32 ) is

    use QuadDobl_Complex_Numbers;

    dim : constant integer32 := v'last - m;
    first : boolean := true;
    maxrad,radval : quad_double;
    idx : constant integer32 := dim + i; -- index of i-th homogeneous variable

  begin
    for k in v'first..dim loop
      if integer32(idz(k)) = i then -- select variable from the i-th set
        if first then
          maxrad := QuadDobl_Complex_Numbers_Polar.Radius(v(k));
          first := false;
        else
          radval := QuadDobl_Complex_Numbers_Polar.Radius(v(k));
          if radval > maxrad
           then maxrad := radval;
          end if;
        end if;
      end if;
    end loop;
    radval := QuadDobl_Complex_Numbers_Polar.Radius(v(idx));
    if radval > maxrad
     then maxrad := radval;
    end if;
    for k in v'first..dim loop
      if integer32(idz(k)) = i then -- select variable from the i-th set
        v(k) := v(k)/maxrad;
      end if;
    end loop;
    v(idx) := v(idx)/maxrad;
  end Scale;

  procedure Scale ( v : in out Standard_Complex_Vectors.Vector;
                    idz : in Standard_Natural_Vectors.Link_to_Vector;
                    m : in integer32 ) is
  begin
    for i in 1..m loop
      Scale(v,idz,m,i);
    end loop;
  end Scale;

  procedure Scale ( v : in out DoblDobl_Complex_Vectors.Vector;
                    idz : in Standard_Natural_Vectors.Link_to_Vector;
                    m : in integer32 ) is
  begin
    for i in 1..m loop
      Scale(v,idz,m,i);
    end loop;
  end Scale;

  procedure Scale ( v : in out QuadDobl_Complex_Vectors.Vector;
                    idz : in Standard_Natural_Vectors.Link_to_Vector;
                    m : in integer32 ) is
  begin
    for i in 1..m loop
      Scale(v,idz,m,i);
    end loop;
  end Scale;

end Hyperplane_Solution_Scaling;
