with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Double_Double_Numbers;             use Double_Double_Numbers;
with DoblDobl_Complex_Numbers;          use DoblDobl_Complex_Numbers;
with DoblDobl_Complex_Numbers_Polar;    use DoblDobl_Complex_Numbers_Polar;

package body DoblDobl_Radial_Solvers is

  function Radii ( c : DoblDobl_Complex_Vectors.Vector )
                 return Double_Double_Vectors.Vector is

    res : Double_Double_Vectors.Vector(c'range);

  begin
    for i in c'range loop
      res(i) := Radius(c(i));
    end loop;
    return res;
  end Radii;

  function Scale ( c : DoblDobl_Complex_Vectors.Vector;
                   r : Double_Double_Vectors.Vector )
                 return DoblDobl_Complex_Vectors.Vector is

    res : DoblDobl_Complex_Vectors.Vector(c'range);

  begin
    for i in c'range loop
      res(i) := c(i)/r(i);
    end loop;
    return res;
  end Scale;

  function Log10 ( r : Double_Double_Vectors.Vector ) 
                 return Double_Double_Vectors.Vector is

    res : Double_Double_Vectors.Vector(r'range);

  begin
    for i in r'range loop
      res(i) := LOG10(r(i));
    end loop;
    return res;
  end Log10;

  function Exp10 ( r : Double_Double_Vectors.Vector ) 
                 return Double_Double_Vectors.Vector is

    res : Double_Double_Vectors.Vector(r'range);
    ten : constant double_double := Double_Double_Numbers.create(10.0);
  
  begin
    for i in r'range loop
      res(i) := ten**r(i);
    end loop;
    return res;
  end Exp10;

  function Radial_Upper_Solve
              ( U : Standard_Integer64_Matrices.Matrix;
                logr : Double_Double_Vectors.Vector )
              return Double_Double_Vectors.Vector is

    res : Double_Double_Vectors.Vector(logr'range);
    eva : double_double;

  begin
    for i in res'range loop
      res(i) := create(0.0);
    end loop;
    for j in U'range loop
      eva := logr(j);
      for i in 1..j-1 loop
        eva := eva - create(integer(U(i,j)))*res(i);
      end loop;
      res(j) := eva/create(integer(U(j,j)));
    end loop;
    return res;
  end Radial_Upper_Solve;

  function Multiply ( A : in Standard_Integer64_Matrices.Matrix;
                      x : in Double_Double_Vectors.Vector )
                    return Double_Double_Vectors.Vector is

    res : Double_Double_Vectors.Vector(x'range);

  begin
    for i in x'range loop
      res(i) := create(0.0);
    end loop;
    for j in A'range(1) loop
      for i in A'range(2) loop
        res(j) := res(j) + create(integer(A(i,j)))*x(i);
      end loop;
    end loop;
    return res;
  end Multiply;

  function Eval ( A : in Standard_Integer64_Matrices.Matrix;
                  x : in Double_Double_Vectors.Vector )
                return Double_Double_Vectors.Vector is

    res : Double_Double_Vectors.Vector(x'range);

  begin
    for i in res'range loop
      res(i) := create(1.0);
    end loop;
    for j in A'range(2) loop
      for i in A'range(1) loop
         res(j) := res(j)*(x(i)**integer(A(i,j)));
      end loop;
    end loop;
    return res;
  end Eval;

  procedure Multiply ( s : in out DoblDobl_Complex_Vectors.Vector;
                       r : in Double_Double_Vectors.Vector ) is
  begin
    for i in s'range loop
      s(i) := s(i)*r(i);
    end loop;
  end Multiply;

  procedure Multiply ( s : in out Solution;
                       r : in Double_Double_Vectors.Vector )  is
  begin
    Multiply(s.v,r);
  end Multiply;

  procedure Multiply ( s : in out Solution_List;
                       r : in Double_Double_Vectors.Vector ) is

    tmp : Solution_List := s;
    ls : Link_to_Solution;

  begin
    while not Is_Null(tmp) loop
      ls := Head_Of(tmp);
      Multiply(ls.v,r);
      Set_Head(tmp,ls);
      tmp := Tail_Of(tmp);
    end loop;
  end Multiply;

end DoblDobl_Radial_Solvers;
