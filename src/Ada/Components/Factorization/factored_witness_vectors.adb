with Double_Double_Numbers;              use Double_Double_Numbers;
with Quad_Double_Numbers;                use Quad_Double_Numbers;

package body Factored_Witness_Vectors is

  procedure Swap ( m : in out Standard_Natural_Vectors.Vector;
                   i,j : in integer32 ) is

    tmp : constant natural32 := m(i);

  begin
    m(i) := m(j);
    m(j) := tmp;
  end Swap;

  procedure Swap ( v : in out Standard_Complex_Vectors.Vector;
                   i,j : in integer32 ) is

    use Standard_Complex_Numbers;
    tmp : constant Complex_Number := v(i);

  begin
    v(i) := v(j);
    v(j) := tmp;
  end Swap;

  procedure Swap ( v : in out DoblDobl_Complex_Vectors.Vector;
                   i,j : in integer32 ) is

    use DoblDobl_Complex_Numbers;
    tmp : constant Complex_Number := v(i);

  begin
    v(i) := v(j);
    v(j) := tmp;
  end Swap;

  procedure Swap ( v : in out QuadDobl_Complex_Vectors.Vector;
                   i,j : in integer32 ) is

    use QuadDobl_Complex_Numbers;
    tmp : constant Complex_Number := v(i);

  begin
    v(i) := v(j);
    v(j) := tmp;
  end Swap;

  procedure Sort ( m : in out Standard_Natural_Vectors.Vector;
                   w : in out Standard_Complex_Vectors.Vector ) is
  begin
    for i in m'range loop
      for j in i+1..m'last loop
        if m(i) > m(j) then
          Swap(m,i,j);
          Swap(w,i,j);
        end if;
      end loop;
    end loop;
  end Sort;

  procedure Sort ( m : in out Standard_Natural_Vectors.Vector;
                   w : in out DoblDobl_Complex_Vectors.Vector ) is
  begin
    for i in m'range loop
      for j in i+1..m'last loop
        if m(i) > m(j) then
          Swap(m,i,j);
          Swap(w,i,j);
        end if;
      end loop;
    end loop;
  end Sort;

  procedure Sort ( m : in out Standard_Natural_Vectors.Vector;
                   w : in out QuadDobl_Complex_Vectors.Vector ) is
  begin
    for i in m'range loop
      for j in i+1..m'last loop
        if m(i) > m(j) then
          Swap(m,i,j);
          Swap(w,i,j);
        end if;
      end loop;
    end loop;
  end Sort;

  function Countmu ( m : Standard_Natural_Vectors.Vector; mu : natural32 )
                   return natural32 is

    res : natural32 := 0;

  begin
    for i in m'range loop
      if m(i) = mu 
       then res := res + 1;
      end if;
    end loop;
    return res;
  end Countmu;

  function Select_Multiple_Factors
               ( m : Standard_Natural_Vectors.Vector;
                 w : Standard_Complex_Vectors.Vector; mu : natural32 ) 
               return Standard_Complex_Vectors.Vector is

    res : Standard_Complex_Vectors.Vector(w'range);
    ind : integer32 := w'first-1;

  begin
    for i in m'range loop
      if m(i) = mu then
        ind := ind + 1;
        res(ind) := w(i);
      end if;
    end loop;
    return res(res'first..ind);
  end Select_Multiple_Factors;

  function Select_Multiple_Factors
               ( m : Standard_Natural_Vectors.Vector;
                 w : DoblDobl_Complex_Vectors.Vector; mu : natural32 ) 
               return DoblDobl_Complex_Vectors.Vector is

    res : DoblDobl_Complex_Vectors.Vector(w'range);
    ind : integer32 := w'first-1;

  begin
    for i in m'range loop
      if m(i) = mu then
        ind := ind + 1;
        res(ind) := w(i);
      end if;
    end loop;
    return res(res'first..ind);
  end Select_Multiple_Factors;

  function Select_Multiple_Factors
               ( m : Standard_Natural_Vectors.Vector;
                 w : QuadDobl_Complex_Vectors.Vector; mu : natural32 ) 
               return QuadDobl_Complex_Vectors.Vector is

    res : QuadDobl_Complex_Vectors.Vector(w'range);
    ind : integer32 := w'first-1;

  begin
    for i in m'range loop
      if m(i) = mu then
        ind := ind + 1;
        res(ind) := w(i);
      end if;
    end loop;
    return res(res'first..ind);
  end Select_Multiple_Factors;

  function Is_In ( v : Standard_Complex_Vectors.Vector;
                   x : Standard_Complex_Numbers.Complex_Number;
                   tol : double_float ) return boolean is

   use Standard_Complex_Numbers;

  begin
    for i in v'range loop
      if AbsVal(v(i)-x) <= tol
       then return true;
      end if;
    end loop;
    return false;
  end Is_In; 

  function Is_In ( v : DoblDobl_Complex_Vectors.Vector;
                   x : DoblDobl_Complex_Numbers.Complex_Number;
                   tol : double_float ) return boolean is

   use DoblDobl_Complex_Numbers;
   val : double_double;

  begin
    for i in v'range loop
      val := AbsVal(v(i)-x);
      if val <= tol
       then return true;
      end if;
    end loop;
    return false;
  end Is_In; 

  function Is_In ( v : QuadDobl_Complex_Vectors.Vector;
                   x : QuadDobl_Complex_Numbers.Complex_Number;
                   tol : double_float ) return boolean is

   use QuadDobl_Complex_Numbers;
   val : quad_double;

  begin
    for i in v'range loop
      val := AbsVal(v(i)-x);
      if val <= tol
       then return true;
      end if;
    end loop;
    return false;
  end Is_In; 

  function Position ( v : Standard_Complex_Vectors.Vector;
                      x : Standard_Complex_Numbers.Complex_Number;
                      tol : double_float ) return integer32 is

    use Standard_Complex_Numbers;

  begin
    for i in v'range loop
      if AbsVal(v(i)-x) <= tol
       then return i;
      end if;
    end loop;
    return v'first-1;
  end Position;

  function Position ( v : DoblDobl_Complex_Vectors.Vector;
                      x : DoblDobl_Complex_Numbers.Complex_Number;
                      tol : double_float ) return integer32 is

    use DoblDobl_Complex_Numbers;
    val : double_double;

  begin
    for i in v'range loop
      val := AbsVal(v(i)-x);
      if val <= tol
       then return i;
      end if;
    end loop;
    return v'first-1;
  end Position;

  function Position ( v : QuadDobl_Complex_Vectors.Vector;
                      x : QuadDobl_Complex_Numbers.Complex_Number;
                      tol : double_float ) return integer32 is

    use QuadDobl_Complex_Numbers;
    val : quad_double;

  begin
    for i in v'range loop
      val := AbsVal(v(i)-x);
      if val <= tol
       then return i;
      end if;
    end loop;
    return v'first-1;
  end Position;

  function Positions ( v,x : Standard_Complex_Vectors.Vector;
                       tol : double_float )
                     return Standard_Natural_Vectors.Vector is

    res : Standard_Natural_Vectors.Vector(x'range);

  begin
    for i in x'range loop
      res(i) := natural32(Position(v,x(i),tol));
    end loop;
    return res;
  end Positions;

  function Positions ( v,x : DoblDobl_Complex_Vectors.Vector;
                       tol : double_float )
                     return Standard_Natural_Vectors.Vector is

    res : Standard_Natural_Vectors.Vector(x'range);

  begin
    for i in x'range loop
      res(i) := natural32(Position(v,x(i),tol));
    end loop;
    return res;
  end Positions;

  function Positions ( v,x : QuadDobl_Complex_Vectors.Vector;
                       tol : double_float )
                     return Standard_Natural_Vectors.Vector is

    res : Standard_Natural_Vectors.Vector(x'range);

  begin
    for i in x'range loop
      res(i) := natural32(Position(v,x(i),tol));
    end loop;
    return res;
  end Positions;

  function Remove_Duplicates
               ( w : Standard_Complex_Vectors.Vector; tol : double_float )
               return Standard_Complex_Vectors.Vector is

    res : Standard_Complex_Vectors.Vector(w'range);
    ind : integer32 := w'first;

  begin
    res(ind) := w(ind);
    for i in w'first+1..w'last loop
      if not Is_In(res(res'first..ind),w(i),tol) then
        ind := ind + 1;
        res(ind) := w(i);
      end if;
    end loop;
    return res(res'first..ind);
  end Remove_Duplicates;

  function Remove_Duplicates
               ( w : DoblDobl_Complex_Vectors.Vector; tol : double_float )
               return DoblDobl_Complex_Vectors.Vector is

    res : DoblDobl_Complex_Vectors.Vector(w'range);
    ind : integer32 := w'first;

  begin
    res(ind) := w(ind);
    for i in w'first+1..w'last loop
      if not Is_In(res(res'first..ind),w(i),tol) then
        ind := ind + 1;
        res(ind) := w(i);
      end if;
    end loop;
    return res(res'first..ind);
  end Remove_Duplicates;

  function Remove_Duplicates
               ( w : QuadDobl_Complex_Vectors.Vector; tol : double_float )
               return QuadDobl_Complex_Vectors.Vector is

    res : QuadDobl_Complex_Vectors.Vector(w'range);
    ind : integer32 := w'first;

  begin
    res(ind) := w(ind);
    for i in w'first+1..w'last loop
      if not Is_In(res(res'first..ind),w(i),tol) then
        ind := ind + 1;
        res(ind) := w(i);
      end if;
    end loop;
    return res(res'first..ind);
  end Remove_Duplicates;

  function Remove_Duplicates
               ( w : Standard_Complex_Vectors.Vector; tol : double_float;
                 m : Standard_Natural_Vectors.Vector )
               return Standard_Natural_Vectors.Vector is

    res : Standard_Natural_Vectors.Vector(w'range) := m;
    ind : integer32 := m'first;

  begin
    for i in w'first+1..w'last loop
      if not Is_In(w(w'first..i-1),w(i),tol) then
        ind := ind + 1;
        res(ind) := m(i);
      end if;
    end loop;
    return res(res'first..ind);
  end Remove_Duplicates;

  function Remove_Duplicates
               ( w : DoblDobl_Complex_Vectors.Vector; tol : double_float;
                 m : Standard_Natural_Vectors.Vector )
               return Standard_Natural_Vectors.Vector is

    res : Standard_Natural_Vectors.Vector(w'range) := m;
    ind : integer32 := m'first;

  begin
    for i in w'first+1..w'last loop
      if not Is_In(w(w'first..i-1),w(i),tol) then
        ind := ind + 1;
        res(ind) := m(i);
      end if;
    end loop;
    return res(res'first..ind);
  end Remove_Duplicates;

  function Remove_Duplicates
               ( w : QuadDobl_Complex_Vectors.Vector; tol : double_float;
                 m : Standard_Natural_Vectors.Vector )
               return Standard_Natural_Vectors.Vector is

    res : Standard_Natural_Vectors.Vector(w'range) := m;
    ind : integer32 := m'first;

  begin
    for i in w'first+1..w'last loop
      if not Is_In(w(w'first..i-1),w(i),tol) then
        ind := ind + 1;
        res(ind) := m(i);
      end if;
    end loop;
    return res(res'first..ind);
  end Remove_Duplicates;

  function Select_Points ( w : in Standard_Complex_Vectors.Vector;
                           k : in Standard_Natural_Vectors.Vector )
                         return Standard_Complex_Vectors.Vector is

    res : Standard_Complex_Vectors.Vector(k'range);

  begin
    for i in k'range loop
      res(i) := w(integer32(k(i)));
    end loop;
    return res;
  end Select_Points;

  function Select_Points ( w : in DoblDobl_Complex_Vectors.Vector;
                           k : in Standard_Natural_Vectors.Vector )
                         return DoblDobl_Complex_Vectors.Vector is

    res : DoblDobl_Complex_Vectors.Vector(k'range);

  begin
    for i in k'range loop
      res(i) := w(integer32(k(i)));
    end loop;
    return res;
  end Select_Points;

  function Select_Points ( w : in QuadDobl_Complex_Vectors.Vector;
                           k : in Standard_Natural_Vectors.Vector )
                         return QuadDobl_Complex_Vectors.Vector is

    res : QuadDobl_Complex_Vectors.Vector(k'range);

  begin
    for i in k'range loop
      res(i) := w(integer32(k(i)));
    end loop;
    return res;
  end Select_Points;

  function Multiplicity_of_Factors
              ( factors : Standard_Natural_VecVecs.VecVec;
                m : in Standard_Natural_Vectors.Vector )
              return Standard_Natural_Vectors.Vector is

    use Standard_Natural_VecVecs,Standard_Natural_Vectors;

    res : Standard_Natural_Vectors.Vector(factors'range);
    ind : integer32 := res'first-1;

  begin
    for i in factors'range loop
      if factors(i) /= null then
        ind := ind + 1;
        res(ind) := m(integer32(factors(i)(1)));
      end if;
    end loop;
    return res(res'first..ind);
  end Multiplicity_of_Factors;

end Factored_Witness_Vectors;
