with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Complex_Numbers;           use Standard_Complex_Numbers;
with Standard_Integer_Vectors;
with Standard_Floating_Vectors;
with Standard_Complex_Vectors;
with Standard_Floating_Linear_Solvers;   use Standard_Floating_Linear_Solvers;
with Standard_Complex_Linear_Solvers;    use Standard_Complex_Linear_Solvers;

package body Standard_Matrix_Inversion is

-- ALGORITHM :
--   Solves repeatedly A*x = b, with b the i-th standard basis vector.

  function Inverse ( m : Standard_Floating_Matrices.Matrix )
                   return Standard_Floating_Matrices.Matrix is

    n : constant integer32 := m'last(1);
    wrk : Standard_Floating_Matrices.Matrix(m'range(1),m'range(2)) := m;
    res : Standard_Floating_Matrices.Matrix(m'range(1),m'range(2));
    piv : Standard_Integer_Vectors.Vector(m'range(1));
    rhs : Standard_Floating_Vectors.Vector(m'range(1));
    inf : integer32;

  begin
    for i in piv'range loop
      piv(i) := i;
    end loop;
    lufac(wrk,n,piv,inf);
    if inf = 0 then
      for j in m'range(2) loop
        rhs := (rhs'range => 0.0);
        rhs(j) := 1.0;
        lusolve(wrk,n,piv,rhs);
        for i in m'range(1) loop
          res(i,j) := rhs(i);
        end loop;
      end loop;
    end if;
    return res;
  end Inverse;

  function Inverse ( m : Standard_Complex_Matrices.Matrix )
                   return Standard_Complex_Matrices.Matrix is

    n : constant integer32 := m'last(1);
    wrk : Standard_Complex_Matrices.Matrix(m'range(1),m'range(2)) := m;
    res : Standard_Complex_Matrices.Matrix(m'range(1),m'range(2));
    piv : Standard_Integer_Vectors.Vector(m'range(1));
    rhs : Standard_Complex_Vectors.Vector(m'range(1));
    inf : integer32;

  begin
    for i in piv'range loop
      piv(i) := i;
    end loop;
    lufac(wrk,n,piv,inf);
    if inf = 0 then
      for j in m'range(2) loop
        rhs := (rhs'range => Create(0.0));
        rhs(j) := Create(1.0);
        lusolve(wrk,n,piv,rhs);
        for i in m'range(1) loop
          res(i,j) := rhs(i);
        end loop;
      end loop;
    end if;
    return res;
  end Inverse;

end Standard_Matrix_Inversion;
