with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Quad_Double_Numbers;              use Quad_Double_Numbers;
with QuadDobl_Complex_Numbers;           use QuadDobl_Complex_Numbers;
with Standard_Integer_Vectors;
with Quad_Double_Vectors;
with QuadDobl_Complex_Vectors;
with Quad_Double_Linear_Solvers;       use Quad_Double_Linear_Solvers;
with QuadDobl_Complex_Linear_Solvers;    use QuadDobl_Complex_Linear_Solvers;

package body QuadDobl_Matrix_Inversion is

-- ALGORITHM :
--   Solves repeatedly A*x = b, with b the i-th standard basis vector.

  function Inverse ( m : Quad_Double_Matrices.Matrix )
                   return Quad_Double_Matrices.Matrix is

    n : constant integer32 := m'last(1);
    wrk : Quad_Double_Matrices.Matrix(m'range(1),m'range(2)) := m;
    res : Quad_Double_Matrices.Matrix(m'range(1),m'range(2));
    piv : Standard_Integer_Vectors.Vector(m'range(1));
    rhs : Quad_Double_Vectors.Vector(m'range(1));
    inf : integer32;
    one : constant quad_double := create(1.0);
    zero : constant quad_double := create(0.0);

  begin
    for i in piv'range loop
      piv(i) := i;
    end loop;
    lufac(wrk,n,piv,inf);
    if inf = 0 then
      for j in m'range(2) loop
        rhs := (rhs'range => zero);
        rhs(j) := one;
        lusolve(wrk,n,piv,rhs);
        for i in m'range(1) loop
          res(i,j) := rhs(i);
        end loop;
      end loop;
    end if;
    return res;
  end Inverse;

  function Inverse ( m : QuadDobl_Complex_Matrices.Matrix )
                   return QuadDobl_Complex_Matrices.Matrix is

    n : constant integer32 := m'last(1);
    wrk : QuadDobl_Complex_Matrices.Matrix(m'range(1),m'range(2)) := m;
    res : QuadDobl_Complex_Matrices.Matrix(m'range(1),m'range(2));
    piv : Standard_Integer_Vectors.Vector(m'range(1));
    rhs : QuadDobl_Complex_Vectors.Vector(m'range(1));
    inf : integer32;
    one : constant quad_double := create(1.0);
    zero : constant quad_double := create(0.0);

  begin
    for i in piv'range loop
      piv(i) := i;
    end loop;
    lufac(wrk,n,piv,inf);
    if inf = 0 then
      for j in m'range(2) loop
        rhs := (rhs'range => Create(zero));
        rhs(j) := Create(one);
        lusolve(wrk,n,piv,rhs);
        for i in m'range(1) loop
          res(i,j) := rhs(i);
        end loop;
      end loop;
    end if;
    return res;
  end Inverse;

end QuadDobl_Matrix_Inversion;
