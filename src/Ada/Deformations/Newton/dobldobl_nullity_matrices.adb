with Standard_Natural_Numbers_io;       use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with DoblDobl_Complex_Polynomials;      use DoblDobl_Complex_Polynomials;
with DoblDobl_Complex_Poly_Functions;   use DoblDobl_Complex_Poly_Functions;
with DoblDobl_Complex_Poly_SysFun;      use DoblDobl_Complex_Poly_SysFun;
with Monomial_Hashing;                  use Monomial_Hashing;
with DoblDobl_Nullity_Polynomials;      use DoblDobl_Nullity_Polynomials;
 
package body DoblDobl_Nullity_Matrices is

  procedure Dimensions_of_Nullity_Matrix
              ( nq,nv,k : in natural32; nr,nc : out natural32 ) is

    nb0,nb1,nr1,nc1 : natural32;

  begin
    if k = 1 then
      nr := nq;
      nc := nv+1;
    else
      nb0 := Monomial_Count(k,nv);
      nb1 := Monomial_Count(k-1,nv);
      Dimensions_of_Nullity_Matrix(nq,nv,k-1,nr1,nc1);
      nr := nr1 + nb1*nq;
      nc := nc1 + nb0;
    end if;
  end Dimensions_of_Nullity_Matrix;

  function Create_Nullity_Matrix
              ( nq,nv,nr,nc,k : natural32; f : Poly_Sys )
              return DoblDobl_Complex_Poly_Matrices.Matrix is

    res : DoblDobl_Complex_Poly_Matrices.Matrix
            (1..integer32(nr),1..integer32(nc));
    nr1,nc1 : natural32;

  begin
    if k = 1 then
      for i in res'range(1) loop
        Copy(f(i),res(i,1));
        for j in 1..integer32(nv) loop
          res(i,j+1) := Diff(f(i),j);
        end loop;
      end loop;
    else
      Dimensions_of_Nullity_Matrix(nq,nv,k-1,nr1,nc1);
      declare
        nm1 : constant DoblDobl_Complex_Poly_Matrices.Matrix
                (1..integer32(nr1),1..integer32(nc1))
            := Create_Nullity_Matrix(nq,nv,nr1,nc1,k-1,f);
      begin
        for i in nm1'range(1) loop
          for j in nm1'range(2) loop
            res(i,j) := nm1(i,j);
          end loop;
        end loop;
      end;
      Compute_Monomial_Multiples(res,nr1+1,1,nq,nv,k,nc1,f);
    end if;
    return res;
  end Create_Nullity_Matrix;

  function Create_Nullity_Matrix
              ( file : file_type;
                nq,nv,nr,nc,k : natural32; f : Poly_Sys )
              return DoblDobl_Complex_Poly_Matrices.Matrix is

    res : DoblDobl_Complex_Poly_Matrices.Matrix
            (1..integer32(nr),1..integer32(nc));
    nr1,nc1 : natural32;

  begin
    put(file,"creating a "); put(file,nr,1); put(file,"-by-");
    put(file,nc,1); put(file," nullity matrix for differential order ");
    put(file,k,1); put_line(file,"...");
    if k = 1 then
      for i in res'range(1) loop
        Copy(f(i),res(i,1));
        for j in 1..integer32(nv) loop
          res(i,j+1) := Diff(f(i),j);
        end loop;
      end loop;
    else
      Dimensions_of_Nullity_Matrix(nq,nv,k-1,nr1,nc1);
      declare
        nm1 : constant DoblDobl_Complex_Poly_Matrices.Matrix
                (1..integer32(nr1),1..integer32(nc1))
            := Create_Nullity_Matrix(file,nq,nv,nr1,nc1,k-1,f);
      begin
        for i in nm1'range(1) loop
          for j in nm1'range(2) loop
            res(i,j) := nm1(i,j);
          end loop;
        end loop;
      end;
      Compute_Monomial_Multiples(file,res,nr1+1,1,nq,nv,k,nc1,f);
    end if;
    return res;
  end Create_Nullity_Matrix;

  function Eval0 ( nm : DoblDobl_Complex_Poly_Matrices.Matrix;
                   z : DoblDobl_Complex_Vectors.Vector )
                 return DoblDobl_Complex_Matrices.Matrix is

    res : DoblDobl_Complex_Matrices.Matrix(nm'range(1),nm'range(2));

  begin
    for i in res'range(1) loop
      for j in res'range(2) loop
        res(i,j) := Eval(nm(i,j),z);
      end loop;
    end loop;
    return res;
  end Eval0;

  function Eval1 ( nm : DoblDobl_Complex_Poly_Matrices.Matrix;
                   z : DoblDobl_Complex_Vectors.Vector )
                 return DoblDobl_Complex_Matrices.Matrix is

    res : DoblDobl_Complex_Matrices.Matrix
            (nm'range(1),nm'first(2)..nm'last(2)-1);

  begin
    for i in res'range(1) loop
      for j in res'range(2) loop
        res(i,j) := Eval(nm(i,j+1),z);
      end loop;
    end loop;
    return res;
  end Eval1;

  function Evaluate_Nullity_Matrix
             ( nq,nv,nr,nc,k : natural32; f : Poly_Sys;
               z : DoblDobl_Complex_Vectors.Vector )
             return DoblDobl_Complex_Matrices.Matrix is

    res : DoblDobl_Complex_Matrices.Matrix
            (1..integer32(nr),1..integer32(nc));
    nr1,nc1 : natural32;
    y : DoblDobl_Complex_Vectors.Vector(1..integer32(nq));
    dp : Poly;

  begin
    if k = 1 then
      y := Eval(f,z);
      for i in y'range loop
        res(i,1) := y(i);
        for j in 1..integer32(nv) loop
          dp := Diff(f(i),j);
          res(i,j+1) := Eval(dp,z);
          Clear(dp);
        end loop;
      end loop;
    else
      Dimensions_of_Nullity_Matrix(nq,nv,k-1,nr1,nc1);
      declare
        eva : constant DoblDobl_Complex_Matrices.Matrix
                (1..integer32(nr1),1..integer32(nc1))
            := Evaluate_Nullity_Matrix(nq,nv,nr1,nc1,k-1,f,z);
      begin
        for i in eva'range(1) loop
          for j in eva'range(2) loop
            res(i,j) := eva(i,j);
          end loop;
        end loop;
      end;
      Evaluate_Monomial_Multiples(res,nr1+1,1,nq,nv,k,nc1,f,z);
    end if;
    return res;
  end Evaluate_Nullity_Matrix;

  function Evaluate_Nullity_Matrix
             ( nq,nv,nr,nc,k : natural32;
               a1 : DoblDobl_Complex_Matrices.Matrix; f : Poly_Sys;
               z : DoblDobl_Complex_Vectors.Vector )
             return DoblDobl_Complex_Matrices.Matrix is

    res : DoblDobl_Complex_Matrices.Matrix
            (1..integer32(nr),1..integer32(nc));
    nr1,nc1 : natural32;
    y : DoblDobl_Complex_Vectors.Vector(1..integer32(nq));

  begin
    if k = 1 then
      y := Eval(f,z);
      for i in y'range loop
        res(i,1) := y(i);
        for j in a1'range(2) loop
          res(i,j+1) := a1(i,j);
        end loop;
      end loop;
    else
      nr1 := natural32(a1'last(1));
      nc1 := natural32(a1'last(2));
      for i in 1..integer32(nr1) loop
        for j in 1..integer32(nc1) loop
          res(i,j) := a1(i,j);
        end loop;
      end loop;
      Evaluate_Monomial_Multiples(res,nr1+1,1,nq,nv,k,nc1,f,z);
    end if;
    return res;
  end Evaluate_Nullity_Matrix;

  function Evaluate_Nullity_Matrix
             ( file : file_type;
               nq,nv,nr,nc,k : natural32; f : Poly_Sys;
               z : DoblDobl_Complex_Vectors.Vector )
             return DoblDobl_Complex_Matrices.Matrix is

    res : DoblDobl_Complex_Matrices.Matrix
            (1..integer32(nr),1..integer32(nc));
    nr1,nc1 : natural32;
    y : DoblDobl_Complex_Vectors.Vector(1..integer32(nq));
    dp : Poly;

  begin
    put(file,"evaluating a "); put(file,nr,1); put(file,"-by-");
    put(file,nc,1); put(file," nullity matrix for differential order ");
    put(file,k,1); put_line(file,"...");
    if k = 1 then
      y := Eval(f,z);
      for i in y'range loop
        res(i,1) := y(i);
        for j in 1..integer32(nv) loop
          dp := Diff(f(i),j);
          res(i,j+1) := Eval(dp,z);
          Clear(dp);
        end loop;
      end loop;
    else
      Dimensions_of_Nullity_Matrix(nq,nv,k-1,nr1,nc1);
      declare
        eva : constant DoblDobl_Complex_Matrices.Matrix
                (1..integer32(nr1),1..integer32(nc1))
            := Evaluate_Nullity_Matrix(file,nq,nv,nr1,nc1,k-1,f,z);
      begin
        for i in eva'range(1) loop
          for j in eva'range(2) loop
            res(i,j) := eva(i,j);
          end loop;
        end loop;
      end;
      Evaluate_Monomial_Multiples(file,res,nr1+1,1,nq,nv,k,nc1,f,z);
    end if;
    return res;
  end Evaluate_Nullity_Matrix;

end DoblDobl_Nullity_Matrices;
