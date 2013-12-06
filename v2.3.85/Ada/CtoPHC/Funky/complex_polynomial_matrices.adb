with unchecked_deallocation;

package body Complex_Polynomial_Matrices is

  function Degrees ( pm : Polynomial_Matrix )
                   return Standard_Integer_Vectors.Vector is

    res : Standard_Integer_Vectors.Vector(1..pm'length(1)*pm'length(2));
    cnt : integer32 := 0;

    use Standard_Complex_Vectors;

  begin
    for i in pm'range(1) loop
      for j in pm'range(2) loop
        cnt := cnt + 1;
        if pm(i,j) = null
         then res(cnt) := -1;
         else res(cnt) := pm(i,j)'last;
        end if;
      end loop;
    end loop;
    return res;
  end Degrees;

  function Degrees ( apm : Array_of_Polynomial_Matrices )
                   return Standard_Integer_Vectors.Vector is

    len : constant integer32 := apm'length*apm(1)'length(1)*apm(1)'length(2);
    res : Standard_Integer_Vectors.Vector(1..len);
    cnt : integer32 := 0;

  begin
    for i in apm'range loop
      declare
        d : constant Standard_Integer_Vectors.Vector := Degrees(apm(i).all);
      begin
        for j in d'range loop
          cnt := cnt + 1;
          res(cnt) := d(j);
        end loop;
      end;
    end loop;
    return res;
  end Degrees;

  function Coefficients ( k : integer32; pm : Polynomial_Matrix )
                        return Standard_Complex_Vectors.Vector is

    res : Standard_Complex_Vectors.Vector(1..k);
    cnt : integer32 := 0;

    use Standard_Complex_Vectors;

  begin
    for i1 in pm'range(1) loop
      for i2 in pm'range(2) loop
        if pm(i1,i2) /= null then
          for j in pm(i1,i2)'range loop
            cnt := cnt + 1;
            res(cnt) := pm(i1,i2)(j);
          end loop;
        end if;
      end loop;
    end loop;
    return res;
  end Coefficients;

  function Coefficients ( k : integer32; apm : Array_of_Polynomial_Matrices )
                        return Standard_Complex_Vectors.Vector is

    res : Standard_Complex_Vectors.Vector(1..k);
    cnt : integer32 := 0;

    use Standard_Complex_Vectors;

  begin
    for i in apm'range loop
      for i1 in apm(i)'range(1) loop
        for i2 in apm(i)'range(2) loop
          if apm(i)(i1,i2) /= null then
            for j in apm(i)(i1,i2)'range loop
              cnt := cnt + 1;
              res(cnt) := apm(i)(i1,i2)(j);
            end loop;
          end if;
        end loop;
      end loop;
    end loop;
    return res;
  end Coefficients;

  function Create ( n,m : integer32;
                    d : Standard_Integer_Vectors.Vector;
                    c : Standard_Complex_Vectors.Vector )
                  return Polynomial_Matrix is

    res : Polynomial_Matrix(1..n,1..m);
    ind_deg : integer32 := 0;
    ind_cff : integer32 := 0;

  begin
    for i1 in 1..n loop
      for i2 in 1..m loop
        ind_deg := ind_deg + 1;
        if d(ind_deg) >= 0 then
          res(i1,i2) := new Standard_Complex_Vectors.Vector(0..d(ind_deg));
          for j in 0..d(ind_deg) loop
            ind_cff := ind_cff + 1;
            res(i1,i2)(j) := c(ind_cff);
          end loop;
        end if;
      end loop;
    end loop;
    return res;
  end Create;

  function Left_Multiply
             ( a : Standard_Complex_Matrices.Matrix; b : Polynomial_Matrix ) 
             return Polynomial_Matrix is

    res : Polynomial_Matrix(a'range(1),b'range(2));
    use Standard_Complex_Vectors;

  begin
    for i in a'range(1) loop
      for j in b'range(2) loop
        for k in a'range(2) loop
          declare
            ab : Vector(b(k,j)'range);
          begin
            for kk in ab'range loop
              ab(kk) := a(i,k)*Conjugate(b(k,j)(kk));
            end loop;
            if res(i,j) = null then
              res(i,j) := new Vector'(ab);
            else
              if ab'length <= res(i,j)'length then
                for kk in ab'range loop
                  res(i,j)(kk) := res(i,j)(kk) + ab(kk);
                end loop;
              else
                for kk in res(i,j)'range loop
                  ab(kk) := ab(kk) + res(i,j)(kk);
                end loop;
                Clear(res(i,j));
                res(i,j) := new Vector'(ab);
              end if;
            end if;
          end;
        end loop;
      end loop;
    end loop;
    return res;
  end Left_Multiply;

  function Left_Multiply ( a : Standard_Complex_Matrices.Matrix;
                           b : Array_of_Polynomial_Matrices )
                         return Array_of_Polynomial_Matrices is

    res : Array_of_Polynomial_Matrices(b'range);
 
  begin
    for i in b'range loop
      res(i) := new Polynomial_Matrix'(Left_Multiply(a,b(i).all));
    end loop;
    return res;
  end Left_Multiply;

  function Left_Multiply
             ( a : Standard_Floating_Matrices.Matrix; b : Polynomial_Matrix ) 
             return Polynomial_Matrix is

    res : Polynomial_Matrix(a'range(1),b'range(2));
    use Standard_Complex_Vectors;

  begin
    for i in a'range(1) loop
      for j in b'range(2) loop
        for k in a'range(2) loop
          declare
            ab : Vector(b(k,j)'range);
          begin
            for kk in ab'range loop
              ab(kk) := Create(a(i,k))*b(k,j)(kk);
            end loop;
            if res(i,j) = null then
              res(i,j) := new Vector'(ab);
            else
              if ab'length <= res(i,j)'length then
                for kk in ab'range loop
                  res(i,j)(kk) := res(i,j)(kk) + ab(kk);
                end loop;
              else
                for kk in res(i,j)'range loop
                  ab(kk) := ab(kk) + res(i,j)(kk);
                end loop;
                Clear(res(i,j));
                res(i,j) := new Vector'(ab);
              end if;
            end if;
          end;
        end loop;
      end loop;
    end loop;
    return res;
  end Left_Multiply;

  function Left_Multiply ( a : Standard_Floating_Matrices.Matrix;
                           b : Array_of_Polynomial_Matrices )
                         return Array_of_Polynomial_Matrices is

    res : Array_of_Polynomial_Matrices(b'range);
 
  begin
    for i in b'range loop
      res(i) := new Polynomial_Matrix'(Left_Multiply(a,b(i).all));
    end loop;
    return res;
  end Left_Multiply;

  function Eval ( p : Standard_Complex_Vectors.Vector;
                  x : Complex_Number ) return Complex_Number is

    res : Complex_Number := p(p'last);

  begin
    for i in reverse 1..p'last loop
      res := res*x;
      res := res + p(i-1);
    end loop;
    return res;
  end Eval;

  function Eval ( pm : Polynomial_Matrix; x : Complex_Number )
                return Standard_Complex_Matrices.Matrix is

    res : Standard_Complex_Matrices.Matrix(pm'range(1),pm'range(2));

    use Standard_Complex_Vectors;

  begin
    for i in pm'range(1) loop
      for j in pm'range(2) loop
        if pm(i,j) = null
         then res(i,j) := Create(0.0);
         else res(i,j) := Eval(pm(i,j).all,x);
        end if;
      end loop;
    end loop;
    return res;
  end Eval;

  procedure Clear ( pm : in out Polynomial_Matrix ) is
  begin
    for i in pm'range(1) loop
      for j in pm'range(2) loop
        Standard_Complex_Vectors.Clear(pm(i,j));
      end loop;
    end loop;
  end Clear;

  procedure Clear ( lpm : in out Link_to_Polynomial_Matrix ) is

    procedure free is
      new unchecked_deallocation(Polynomial_Matrix,Link_to_Polynomial_Matrix);

  begin
    if lpm /= null
     then Clear(lpm.all);
          free(lpm);
    end if;
  end Clear;

end Complex_Polynomial_Matrices;
