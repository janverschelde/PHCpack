with Standard_Integer_Vectors;
with Standard_Random_Vectors;
with Standard_Random_Matrices;

package body Standard_Random_Laurentials is

  function Create_Laurent_Polynomial
             ( c : Standard_Complex_Vectors.Vector;
               e : Standard_Integer_Matrices.Matrix ) return Poly is

    res : Poly := Null_Poly;
    t : Term;

  begin
    t.dg := new Standard_Integer_Vectors.Vector(e'range(1));
    for j in c'range loop
      t.cf := c(j);
      for i in e'range(1) loop
        t.dg(i) := e(i,j);
      end loop;
      Add(res,t);
    end loop;
    Clear(t);
    return res;
  end Create_Laurent_Polynomial;

  function Random_Laurent_Polynomial 
             ( n,m : natural32; lower,upper : integer32 ) return Poly is

    e : constant Standard_Integer_Matrices.Matrix
          (1..integer32(n),1..integer32(m))
      := Standard_Random_Matrices.Random_Matrix(n,m,lower,upper);
    c : constant Standard_Complex_Vectors.Vector
      := Standard_Random_Vectors.Random_Vector(1,integer32(m));
    res : constant Poly := Create_Laurent_Polynomial(c,e);

  begin
    return res;
  end Random_Laurent_Polynomial;

end Standard_Random_Laurentials;
