with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Double_Double_Numbers;              use Double_Double_Numbers;
with DoblDobl_Complex_Numbers;           use DoblDobl_Complex_Numbers;
with Standard_Integer_Vectors;
with Double_Double_Vectors;
with DoblDobl_Complex_Vectors;
with Double_Double_Linear_Solvers;       use Double_Double_Linear_Solvers;
with DoblDobl_Complex_Linear_Solvers;    use DoblDobl_Complex_Linear_Solvers;

package body DoblDobl_Matrix_Inversion is

-- ALGORITHM :
--   Solves repeatedly A*x = b, with b the i-th standard basis vector.

  function Inverse ( m : Double_Double_Matrices.Matrix )
                   return Double_Double_Matrices.Matrix is

    n : constant integer32 := m'last(1);
    wrk : Double_Double_Matrices.Matrix(m'range(1),m'range(2)) := m;
    res : Double_Double_Matrices.Matrix(m'range(1),m'range(2));
    piv : Standard_Integer_Vectors.Vector(m'range(1));
    rhs : Double_Double_Vectors.Vector(m'range(1));
    inf : integer32;
    one : constant double_double := create(1.0);
    zero : constant double_double := create(0.0);

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

  function Inverse ( m : DoblDobl_Complex_Matrices.Matrix )
                   return DoblDobl_Complex_Matrices.Matrix is

    n : constant integer32 := m'last(1);
    wrk : DoblDobl_Complex_Matrices.Matrix(m'range(1),m'range(2)) := m;
    res : DoblDobl_Complex_Matrices.Matrix(m'range(1),m'range(2));
    piv : Standard_Integer_Vectors.Vector(m'range(1));
    rhs : DoblDobl_Complex_Vectors.Vector(m'range(1));
    inf : integer32;
    one : constant double_double := create(1.0);
    zero : constant double_double := create(0.0);

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

end DoblDobl_Matrix_Inversion;
