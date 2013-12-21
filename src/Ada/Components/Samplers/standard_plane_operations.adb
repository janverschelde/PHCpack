with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers_io;      use Standard_Floating_Numbers_io;
with Standard_Complex_Numbers_io;       use Standard_Complex_Numbers_io;
with Standard_Complex_Norms_Equals;     use Standard_Complex_Norms_Equals;
with Standard_Random_Vectors;           use Standard_Random_Vectors;
with Standard_Plane_Representations;    use Standard_Plane_Representations;
with Standard_Point_Coordinates;        use Standard_Point_Coordinates;

package body Standard_Plane_Operations is

  procedure Random_Affine_Plane ( n,k : in integer32;
                                  b : out Vector; v : out VecVec ) is
  begin
    b := Random_Vector(1,n);
    for i in v'range loop
      if v(i) = null
       then v(i) := new Standard_Complex_Vectors.Vector(1..n);
      end if;
      v(i).all := Random_Vector(1,n);
    end loop;
  end Random_Affine_Plane;

  function Random_Point ( b : Vector; v : VecVec ) return Vector is

    c : constant Vector(v'range) := Random_Vector(v'first,v'last);
    res : constant Vector(b'range) := Affine_Expand(c,b,v);

  begin
    return res;
  end Random_Point;

  function Evaluate ( h,x : Vector ) return Complex_Number is

    res : constant Complex_Number := h(0) + h(x'range)*x;

  begin
    return res;
  end Evaluate;

  function Evaluate ( h : VecVec; x : Vector ) return Vector is

    res : Vector(h'range);
    
  begin
    for i in h'range loop
      res(i) := Evaluate(h(i).all,x);
    end loop;
    return res;
  end Evaluate;

  function Orthogonalize ( v : Array_of_VecVecs ) return Array_of_VecVecs is

    res : Array_of_VecVecs(v'range);

  begin
    for i in v'range loop
      res(i) := new VecVec'(Orthogonalize(v(i).all));
    end loop;
    return res;
  end Orthogonalize;

  function In_Span ( v : VecVec; x : Standard_Complex_Vectors.Vector;
                     tol : double_float ) return boolean is

    p : Vector(x'range) := x;
    ip : Complex_Number;

  begin
    for i in v'range loop
      ip := Inner_Product(p,v(i).all);
      p := p - ip*v(i).all;
    end loop;
    return (Max_Norm(p) <= tol);
  end In_Span;

  procedure Affine_Orthonormal_Basis
              ( n,k : in integer32; slices : in VecVec;
                b : out Vector; v,w : out VecVec ) is

    eqs : VecVec(slices'range);

  begin
    for i in 1..k loop
      eqs(i) := new Standard_Complex_Vectors.Vector(0..n);
      for j in 0..n loop
        eqs(i)(j) := slices(i)(j);
      end loop;
    end loop;
    Generators(n,k,eqs,b,v);
    w := Orthogonalize(v);
  end Affine_Orthonormal_Basis;

  function Truncate ( v : in VecVec; n : in integer32 ) return VecVec is

    res : VecVec(v'range);

  begin
    for i in v'range loop
      if v(i) /= null
       then res(i) := new Standard_Complex_Vectors.Vector'(v(i)(1..n));
      end if;
    end loop;
    return res;
  end Truncate;

  function Truncate ( v : in Array_of_VecVecs; n : in integer32 )
                    return Array_of_VecVecs is

    res : Array_of_VecVecs(v'range);

  begin
    for i in v'range loop
      if v(i) /= null
       then res(i) := new VecVec'(Truncate(v(i).all,n));
      end if;
    end loop;
    return res;
  end Truncate;

  procedure Evaluate ( file : in file_type;
                       equ : in Matrix; v : in Vector;
                       res : out double_float ) is

    y : Complex_Number;

  begin
    res := 0.0;
    for i in equ'range loop
      put(file,"at equation ");
      put(file,i,1); put(file," : ");
      y := equ(i,0);
      for j in v'range loop
        y := y + equ(i,j)*v(j);
      end loop;
      put(file,y); new_line(file);
      res := res + AbsVal(y);
    end loop;
  end Evaluate;

  procedure Evaluate ( equ : in Matrix; v : in Vector;
                       res : out double_float ) is

    y : Complex_Number;

  begin
    res := 0.0;
    for i in equ'range loop
      y := equ(i,0);
      for j in v'range loop
        y := y + equ(i,j)*v(j);
      end loop;
      res := res + AbsVal(y);
    end loop;
  end Evaluate;

  procedure Evaluate ( file : in file_type;
                       p : in Matrix; g : in Matrix;
                       res : out double_float ) is

    y : Complex_Number;
    rgk : double_float;

  begin
    res := 0.0;
    for k in g'range(2) loop
      rgk := 0.0;
     -- put(file,"Verifying generator ");
     -- put(file,k,1); put_line(file,"...");
      for i in p'range(1) loop
       -- put(file,"  at equation "); put(file,i,1); 
        if k = 0
         then y := p(i,0);
         else y := Create(0.0);
        end if;
        for j in 1..p'last(2) loop
          y := y + p(i,j)*g(j,k);
        end loop;
       -- put(file," : "); put(file,y); new_line(file);
        res := res + AbsVal(y);
        rgk := rgk + AbsVal(y);
      end loop;
      put(file,"Value at generator "); put(file,k,1);
      put(file," : "); put(file,rgk,3); new_line(file);
    end loop;
  end Evaluate;

  procedure Evaluate ( p : in Matrix; g : in Matrix;
                       res : out double_float ) is

    y : Complex_Number;

  begin
    res := 0.0;
    for k in g'range(2) loop
      for i in p'range(1) loop
        if k = 0
         then y := p(i,0);
         else y := Create(0.0);
        end if;
        for j in 1..p'last(2) loop
          y := y + p(i,j)*g(j,k);
        end loop;
        res := res + AbsVal(y);
      end loop;
    end loop;
  end Evaluate;

  procedure Intersect ( e1,e2 : in Matrix; p1,p2 : in out Matrix ) is

    eva0,eva1,eva2,t : Complex_Number;
    v : Vector(p2'range(1));
   -- res : double_float;

  begin
    eva0 := Create(0.0);
    eva1 := Create(0.0);
    for i in 1..e1'last(2) loop
      eva0 := eva0 + e1(1,i)*p2(i,0);
      eva1 := eva1 + e1(1,i)*p2(i,1);
    end loop;
    for k in 2..p2'last(2) loop
      eva2 := Create(0.0);
      for i in 1..e1'last(2) loop
        eva2 := eva2 + e1(1,i)*p2(i,k);
      end loop;
      t := (e1(1,0) + eva0 + eva1)/eva2;
      for i in p2'range(1) loop
        v(i) := p2(i,0) + p2(i,1) - t*p2(i,k);
      end loop;
     -- put_line("Verify if vector satisfies 1st equations : ");
     -- Evaluate(Standard_Output,e1,v,res);
     -- put("  residual : "); put(res,3); new_line;
     -- Evaluate(e1,v,res);
     -- put("Value of vector at 1st plane : "); put(res,3); new_line;
     -- put_line("Verify if vector satisfies 2nd equations : ");
     -- Evaluate(Standard_Output,e2,v,res);
     -- put("  residual : "); put(res,3); new_line;
     -- Evaluate(e2,v,res);
     -- put("Value of vector at 1st plane : "); put(res,3); new_line;
      for i in v'range loop
        p1(i,k-1) := v(i) - p1(i,0);
      end loop;
    end loop;
    for k in 2..p2'last(2) loop
      for i in p2'range(1) loop
        p2(i,k-1) := p1(i,k-1);
      end loop;
    end loop;
    p1 := Orthogonalize(p1);
    p2 := Orthogonalize(p2);
  end Intersect;

  procedure Intersect ( file : in file_type;
                        e1,e2 : in Matrix; p1,p2 : in out Matrix ) is

    eva0,eva1,eva2,t : Complex_Number;
    v : Vector(p2'range(1));
    res : double_float;

  begin
    eva0 := Create(0.0);
    eva1 := Create(0.0);
    for i in 1..e1'last(2) loop
      eva0 := eva0 + e1(1,i)*p2(i,0);
      eva1 := eva1 + e1(1,i)*p2(i,1);
    end loop;
    for k in 2..p2'last(2) loop
      eva2 := Create(0.0);
      for i in 1..e1'last(2) loop
        eva2 := eva2 + e1(1,i)*p2(i,k);
      end loop;
      t := (e1(1,0) + eva0 + eva1)/eva2;
      for i in p2'range(1) loop
        v(i) := p2(i,0) + p2(i,1) - t*p2(i,k);
      end loop;
      put_line(file,"Verify if vector satisfies 1st equations : ");
      Evaluate(file,e1,v,res);
      put(file,"  residual : "); put(file,res,3); new_line(file);
      Evaluate(e1,v,res);
      put(file,"Value of vector at 1st plane : ");
      put(file,res,3); new_line(file);
      put_line(file,"Verify if vector satisfies 2nd equations : ");
      Evaluate(file,e2,v,res);
      put(file,"  residual : "); put(file,res,3); new_line(file);
      Evaluate(e2,v,res);
      put(file,"Value of vector at 1st plane : ");
      put(file,res,3); new_line(file);
      for i in v'range loop
        p1(i,k-1) := v(i) - p1(i,0);
      end loop;
    end loop;
    for k in 2..p2'last(2) loop
      for i in p2'range(1) loop
        p2(i,k-1) := p1(i,k-1);
      end loop;
    end loop;
    p1 := Orthogonalize(p1);
    p2 := Orthogonalize(p2);
  end Intersect;

end Standard_Plane_Operations;
