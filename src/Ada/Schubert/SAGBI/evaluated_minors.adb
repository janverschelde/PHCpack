with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Vectors;           use Standard_Integer_Vectors;
with Standard_Floating_Linear_Solvers;   use Standard_Floating_Linear_Solvers;
with Standard_Complex_Linear_Solvers;    use Standard_Complex_Linear_Solvers;
with Double_Double_Linear_Solvers;       use Double_Double_Linear_Solvers;
with DoblDobl_Complex_Linear_Solvers;    use DoblDobl_Complex_Linear_Solvers;
with Quad_Double_Linear_Solvers;         use Quad_Double_Linear_Solvers;
with QuadDobl_Complex_Linear_Solvers;    use QuadDobl_Complex_Linear_Solvers;

package body Evaluated_Minors is

  function Determinant ( m : Standard_Floating_Matrices.Matrix )
                       return double_float is

    res : double_float;
    wrk : Standard_Floating_Matrices.Matrix(m'range(1),m'range(2));
    piv : Standard_Integer_Vectors.Vector(m'range(1));
    inf : integer32;

  begin
    for i in m'range(1) loop
      piv(i) := i;
      for j in m'range(2) loop
        wrk(i,j) := m(i,j);
      end loop;
    end loop;
    lufac(wrk,m'last(1),piv,inf);
    if inf /= 0 then
      res := 0.0;
    else
      res := 1.0;
      for i in m'range(1) loop
        res := res*wrk(i,i);
      end loop;
      for i in piv'range loop
        if piv(i) > i
         then res := -res;
        end if;
      end loop;
    end if;
    return res;
  end Determinant;

  function Determinant ( m : Double_Double_Matrices.Matrix )
                       return double_double is

    res : double_double;
    wrk : Double_Double_Matrices.Matrix(m'range(1),m'range(2));
    piv : Standard_Integer_Vectors.Vector(m'range(1));
    inf : integer32;

  begin
    for i in m'range(1) loop
      piv(i) := i;
      for j in m'range(2) loop
        wrk(i,j) := m(i,j);
      end loop;
    end loop;
    lufac(wrk,m'last(1),piv,inf);
    if inf /= 0 then
      res := create(0.0);
    else
      res := create(1.0);
      for i in m'range(1) loop
        res := res*wrk(i,i);
      end loop;
      for i in piv'range loop
        if piv(i) > i
         then res := -res;
        end if;
      end loop;
    end if;
    return res;
  end Determinant;

  function Determinant ( m : Quad_Double_Matrices.Matrix )
                       return quad_double is

    res : quad_double;
    wrk : Quad_Double_Matrices.Matrix(m'range(1),m'range(2));
    piv : Standard_Integer_Vectors.Vector(m'range(1));
    inf : integer32;

  begin
    for i in m'range(1) loop
      piv(i) := i;
      for j in m'range(2) loop
        wrk(i,j) := m(i,j);
      end loop;
    end loop;
    lufac(wrk,m'last(1),piv,inf);
    if inf /= 0 then
      res := create(0.0);
    else
      res := create(1.0);
      for i in m'range(1) loop
        res := res*wrk(i,i);
      end loop;
      for i in piv'range loop
        if piv(i) > i
         then res := -res;
        end if;
      end loop;
    end if;
    return res;
  end Determinant;

  function Determinant ( m : Standard_Floating_Matrices.Matrix; b : Bracket )
                       return double_float is

    res : double_float;
    sqm : Standard_Floating_Matrices.Matrix(b'range,b'range);
    piv : Standard_Integer_Vectors.Vector(b'range);
    inf : integer32;

  begin
    for i in b'range loop
      piv(i) := i;
      for j in b'range loop
        sqm(i,j) := m(integer32(b(i)),j);
      end loop;
    end loop;
    lufac(sqm,b'last,piv,inf);
    if inf /= 0 then
      res := 0.0;
    else
      res := 1.0;
      for i in b'range loop
        res := res*sqm(i,i);
      end loop;
      for i in piv'range loop
        if piv(i) > i
         then res := -res;
        end if;
      end loop;
    end if;
    return res;
  end Determinant;

  function Determinant ( m : Double_Double_Matrices.Matrix; b : Bracket )
                       return double_double is

    res : double_double;
    sqm : Double_Double_Matrices.Matrix(b'range,b'range);
    piv : Standard_Integer_Vectors.Vector(b'range);
    inf : integer32;

  begin
    for i in b'range loop
      piv(i) := i;
      for j in b'range loop
        sqm(i,j) := m(integer32(b(i)),j);
      end loop;
    end loop;
    lufac(sqm,b'last,piv,inf);
    if inf /= 0 then
      res := create(0.0);
    else
      res := create(1.0);
      for i in b'range loop
        res := res*sqm(i,i);
      end loop;
      for i in piv'range loop
        if piv(i) > i
         then res := -res;
        end if;
      end loop;
    end if;
    return res;
  end Determinant;

  function Determinant ( m : Quad_Double_Matrices.Matrix; b : Bracket )
                       return quad_double is

    res : quad_double;
    sqm : Quad_Double_Matrices.Matrix(b'range,b'range);
    piv : Standard_Integer_Vectors.Vector(b'range);
    inf : integer32;

  begin
    for i in b'range loop
      piv(i) := i;
      for j in b'range loop
        sqm(i,j) := m(integer32(b(i)),j);
      end loop;
    end loop;
    lufac(sqm,b'last,piv,inf);
    if inf /= 0 then
      res := create(0.0);
    else
      res := create(1.0);
      for i in b'range loop
        res := res*sqm(i,i);
      end loop;
      for i in piv'range loop
        if piv(i) > i
         then res := -res;
        end if;
      end loop;
    end if;
    return res;
  end Determinant;

  function Determinant ( m : Standard_Complex_Matrices.Matrix )
                       return Standard_Complex_Numbers.Complex_Number is

    use Standard_Complex_Numbers;

    res : Complex_Number;
    wrk : Standard_Complex_Matrices.Matrix(m'range(1),m'range(2));
    piv : Standard_Integer_Vectors.Vector(m'range(1));
    inf : integer32;

  begin
    for i in m'range(1) loop
      piv(i) := i;
      for j in m'range(2) loop
        wrk(i,j) := m(i,j);
      end loop;
    end loop;
    lufac(wrk,wrk'last(1),piv,inf);
    if inf /= 0 then
      res := Create(0.0);
    else
      res := Create(1.0);
      for i in wrk'range(1) loop
        res := res*wrk(i,i);
      end loop;
      for i in piv'range loop
        if piv(i) > i
         then res := -res;
        end if;
      end loop;
    end if;
    return res;
  end Determinant;

  function Determinant ( m : DoblDobl_Complex_Matrices.Matrix )
                       return DoblDobl_Complex_Numbers.Complex_Number is

    use DoblDobl_Complex_Numbers;

    res : Complex_Number;
    wrk : DoblDobl_Complex_Matrices.Matrix(m'range(1),m'range(2));
    piv : Standard_Integer_Vectors.Vector(m'range(1));
    inf : integer32;
    zero : constant double_double := create(0.0);
    one : constant double_double := create(1.0);

  begin
    for i in m'range(1) loop
      piv(i) := i;
      for j in m'range(2) loop
        wrk(i,j) := m(i,j);
      end loop;
    end loop;
    lufac(wrk,wrk'last(1),piv,inf);
    if inf /= 0 then
      res := Create(zero);
    else
      res := Create(one);
      for i in wrk'range(1) loop
        res := res*wrk(i,i);
      end loop;
      for i in piv'range loop
        if piv(i) > i
         then res := -res;
        end if;
      end loop;
    end if;
    return res;
  end Determinant;

  function Determinant ( m : QuadDobl_Complex_Matrices.Matrix )
                       return QuadDobl_Complex_Numbers.Complex_Number is

    use QuadDobl_Complex_Numbers;

    res : Complex_Number;
    wrk : QuadDobl_Complex_Matrices.Matrix(m'range(1),m'range(2));
    piv : Standard_Integer_Vectors.Vector(m'range(1));
    inf : integer32;
    zero : constant quad_double := create(0.0);
    one : constant quad_double := create(1.0);

  begin
    for i in m'range(1) loop
      piv(i) := i;
      for j in m'range(2) loop
        wrk(i,j) := m(i,j);
      end loop;
    end loop;
    lufac(wrk,wrk'last(1),piv,inf);
    if inf /= 0 then
      res := Create(zero);
    else
      res := Create(one);
      for i in wrk'range(1) loop
        res := res*wrk(i,i);
      end loop;
      for i in piv'range loop
        if piv(i) > i
         then res := -res;
        end if;
      end loop;
    end if;
    return res;
  end Determinant;

  function Determinant ( m : Standard_Complex_Matrices.Matrix; b : Bracket )
                       return Standard_Complex_Numbers.Complex_Number is

    use Standard_Complex_Numbers;

    res : Complex_Number;
    sqm : Standard_Complex_Matrices.Matrix(b'range,b'range);
    piv : Standard_Integer_Vectors.Vector(b'range);
    inf : integer32;

  begin
    for i in b'range loop
      piv(i) := i;
      for j in b'range loop
        sqm(i,j) := m(integer32(b(i)),j);
      end loop;
    end loop;
    lufac(sqm,b'last,piv,inf);
    if inf /= 0 then
      res := Create(0.0);
    else
      res := Create(1.0);
      for i in sqm'range(1) loop
        res := res*sqm(i,i);
      end loop;
      for i in piv'range loop
        if piv(i) > i
         then res := -res;
        end if;
      end loop;
    end if;
    return res;
  end Determinant;

  function Determinant ( m : DoblDobl_Complex_Matrices.Matrix; b : Bracket )
                       return DoblDobl_Complex_Numbers.Complex_Number is

    use DoblDobl_Complex_Numbers;

    res : Complex_Number;
    sqm : DoblDobl_Complex_Matrices.Matrix(b'range,b'range);
    piv : Standard_Integer_Vectors.Vector(b'range);
    inf : integer32;
    zero : constant double_double := create(0.0);
    one : constant double_double := create(1.0);

  begin
    for i in b'range loop
      piv(i) := i;
      for j in b'range loop
        sqm(i,j) := m(integer32(b(i)),j);
      end loop;
    end loop;
    lufac(sqm,b'last,piv,inf);
    if inf /= 0 then
      res := Create(zero);
    else
      res := Create(one);
      for i in sqm'range(1) loop
        res := res*sqm(i,i);
      end loop;
      for i in piv'range loop
        if piv(i) > i
         then res := -res;
        end if;
      end loop;
    end if;
    return res;
  end Determinant;

  function Determinant ( m : QuadDobl_Complex_Matrices.Matrix; b : Bracket )
                       return QuadDobl_Complex_Numbers.Complex_Number is

    use QuadDobl_Complex_Numbers;

    res : Complex_Number;
    sqm : QuadDobl_Complex_Matrices.Matrix(b'range,b'range);
    piv : Standard_Integer_Vectors.Vector(b'range);
    inf : integer32;
    zero : constant quad_double := create(0.0);
    one : constant quad_double := create(1.0);

  begin
    for i in b'range loop
      piv(i) := i;
      for j in b'range loop
        sqm(i,j) := m(integer32(b(i)),j);
      end loop;
    end loop;
    lufac(sqm,b'last,piv,inf);
    if inf /= 0 then
      res := Create(zero);
    else
      res := Create(one);
      for i in sqm'range(1) loop
        res := res*sqm(i,i);
      end loop;
      for i in piv'range loop
        if piv(i) > i
         then res := -res;
        end if;
      end loop;
    end if;
    return res;
  end Determinant;

  function Determinant ( m : Standard_Complex_Matrices.Matrix;
                         rows,cols : Bracket )
                       return Standard_Complex_Numbers.Complex_Number is

    use Standard_Complex_Numbers;

    res : Complex_Number;
    sqm : Standard_Complex_Matrices.Matrix(rows'range,cols'range);
    piv : Standard_Integer_Vectors.Vector(rows'range);
    inf : integer32;

  begin
    for i in rows'range loop
      piv(i) := i;
      for j in cols'range loop
        sqm(i,j) := m(integer32(rows(i)),integer32(cols(j)));
      end loop;
    end loop;
    lufac(sqm,rows'last,piv,inf);
    if inf /= 0 then
      res := Create(0.0);
    else
      res := Create(1.0);
      for i in sqm'range(1) loop
        res := res*sqm(i,i);
      end loop;
      for i in piv'range loop
        if piv(i) > i
         then res := -res;
        end if;
      end loop;
    end if;
    return res;
  end Determinant;

  function Determinant ( m : DoblDobl_Complex_Matrices.Matrix;
                         rows,cols : Bracket )
                       return DoblDobl_Complex_Numbers.Complex_Number is

    use DoblDobl_Complex_Numbers;

    res : Complex_Number;
    sqm : DoblDobl_Complex_Matrices.Matrix(rows'range,cols'range);
    piv : Standard_Integer_Vectors.Vector(rows'range);
    inf : integer32;
    zero : constant double_double := create(0.0);
    one : constant double_double := create(1.0);

  begin
    for i in rows'range loop
      piv(i) := i;
      for j in cols'range loop
        sqm(i,j) := m(integer32(rows(i)),integer32(cols(j)));
      end loop;
    end loop;
    lufac(sqm,rows'last,piv,inf);
    if inf /= 0 then
      res := Create(zero);
    else
      res := Create(one);
      for i in sqm'range(1) loop
        res := res*sqm(i,i);
      end loop;
      for i in piv'range loop
        if piv(i) > i
         then res := -res;
        end if;
      end loop;
    end if;
    return res;
  end Determinant;

  function Determinant ( m : QuadDobl_Complex_Matrices.Matrix;
                         rows,cols : Bracket )
                       return QuadDobl_Complex_Numbers.Complex_Number is

    use QuadDobl_Complex_Numbers;

    res : Complex_Number;
    sqm : QuadDobl_Complex_Matrices.Matrix(rows'range,cols'range);
    piv : Standard_Integer_Vectors.Vector(rows'range);
    inf : integer32;
    zero : constant quad_double := create(0.0);
    one : constant quad_double := create(1.0);

  begin
    for i in rows'range loop
      piv(i) := i;
      for j in cols'range loop
        sqm(i,j) := m(integer32(rows(i)),integer32(cols(j)));
      end loop;
    end loop;
    lufac(sqm,rows'last,piv,inf);
    if inf /= 0 then
      res := Create(zero);
    else
      res := Create(one);
      for i in sqm'range(1) loop
        res := res*sqm(i,i);
      end loop;
      for i in piv'range loop
        if piv(i) > i
         then res := -res;
        end if;
      end loop;
    end if;
    return res;
  end Determinant;

end Evaluated_Minors;
