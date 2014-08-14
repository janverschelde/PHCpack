with text_io;                           use text_io;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with Double_Double_Numbers_io;          use Double_Double_Numbers_io;
with DoblDobl_Complex_Numbers_io;       use DoblDobl_Complex_Numbers_io;
with DoblDobl_Complex_Vector_Norms;     use DoblDobl_Complex_Vector_Norms;

package body DoblDobl_Rescaling_Coordinates is

  function Linear_Offset_Shift
             ( a,b : Vector; t : Complex_Number ) return Vector is

    res : Vector(b'range);
    dd_one : constant double_double := create(1.0);
    one : constant Complex_Number := Create(dd_one);
    s : constant Complex_Number := one - t;

  begin
    for i in b'range loop
      res(i) := s*a(i) + t*b(i);
    end loop;
    return res;
  end Linear_Offset_Shift;

  function Complement_of_Projection
             ( p : Matrix; v : Vector ) return Vector is

    res : Vector(v'range) := v;
    w : Vector(v'range);
    c : Complex_Number;

  begin
    for j in 1..p'last(2) loop
      for i in w'range loop
        w(i) := p(i,j);
      end loop;
      c := Conjugated_Inner_Product(w,res);
      for i in v'range loop
        res(i) := res(i) - c*w(i);
      end loop;
    end loop;
    return res;
  end Complement_of_Projection;

  function Distance ( p : Matrix; x : Vector ) return double_double is

    bx,c : Vector(p'range(1));

  begin
    for i in bx'range loop
      bx(i) := p(i,0) - x(i);
    end loop;
    c := Complement_of_Projection(p,bx);
    return Norm2(c);
  end Distance;

  procedure Check_Orthonormality
              ( p : in Matrix; tol : in double_double;
                fail : out boolean ) is

    k : constant integer32 := p'last(2);
    nrm : double_double;
    ip : Complex_Number;
    v,w : Vector(p'range(1));

  begin
    fail := false;
    put("Checking the orthonormality of a ");
    put(k,1); put_line("-plane ...");
    for i in 1..k loop
      for j in v'range loop
        v(j) := p(j,i);
      end loop;
      nrm := Norm2(v);
      put("vector "); put(i,1); put(" has norm ");
      put(nrm); new_line;
      fail := fail or (AbsVal(nrm-1.0) > tol);
      for j in i+1..k loop
        for ii in w'range loop 
          w(ii) := p(ii,j);
        end loop;
        ip := Conjugated_Inner_Product(v,w);
        put("product with vector "); put(j,1);
        put(" : "); put(ip); new_line;
        fail := fail or (AbsVal(ip) > tol);
      end loop;
    end loop;
  end Check_Orthonormality;

end DoblDobl_Rescaling_Coordinates;
