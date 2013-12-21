with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Chebychev_Polynomials;              use Chebychev_Polynomials;
with Standard_Integer_Vectors;
with Standard_Floating_Vectors;          use Standard_Floating_Vectors;
with Standard_Random_Numbers;            use Standard_Random_Numbers;
with Standard_Floating_Linear_Solvers;   use Standard_Floating_Linear_Solvers;
with Standard_Floating_QR_Least_Squares; use Standard_Floating_QR_Least_Squares;

package body Osculating_Planes is

  function Standard_Basis
                ( n,d : natural32; s : double_float ) return Matrix is

    res : Matrix(1..integer32(n),1..integer32(d));
    acc : double_float;
    fac,lim : integer32;

  begin
    for i in 1..integer32(d) loop         -- set the zeros and ones
      res(i,i) := 1.0;
      for j in (i+1)..integer32(d) loop
        res(i,j) := 0.0;
      end loop;
    end loop;
    for j in 1..integer32(d) loop         -- set the powers of the s-values
      acc := s;
      for i in (j+1)..integer32(n) loop
        res(i,j) := acc;
        acc := acc*s;
      end loop;
    end loop;
    for i in 3..integer32(n) loop      -- compute the factors from derivation
      fac := 1;
      if i-1 > integer32(d)
       then lim := integer32(d);
       else lim := i-1;
      end if;
      for j in 2..lim loop
        fac := fac*(i+1-j);
        res(i,j) := double_float(fac)*res(i,j);
      end loop;
      if i <= integer32(d)
       then res(i,i) := double_float(fac);
      end if;
    end loop;
    for j in 3..integer32(d) loop                     -- scale the columns
      for i in (j+1)..integer32(n) loop
        res(i,j) := res(i,j)/res(j,j);
      end loop;
      res(j,j) := 1.0;
    end loop;
    return res;
  end Standard_Basis;

  function Chebychev_Basis
                ( n,d : natural32; s : double_float ) return Matrix is

    res : Matrix(1..integer32(n),1..integer32(d));
    lim : integer32;

  begin
    for i in 1..integer32(d) loop              -- set the zeros and ones
      res(i,i) := 1.0;
      for j in (i+1)..integer32(d) loop
        res(i,j) := 0.0;
      end loop;
    end loop;
    for i in 2..integer32(n) loop              -- evaluate and differentiate
      declare
        p : constant Vector := Chebychev_Polynomials.Create(natural32(i-1));
      begin
        res(i,1) := Eval(p,s);
        if i > integer32(d)
         then lim := integer32(d);
         else lim := i;
        end if;
        for j in 2..lim loop
          declare
            dp : constant Vector
               := Chebychev_Polynomials.Diff(p,natural32(j-1));
          begin
            res(i,j) := Eval(dp,s);
          end;
        end loop;
      end;
    end loop;
    for j in 3..integer32(d) loop                     -- scale the columns
      for i in (j+1)..integer32(n) loop
        res(i,j) := res(i,j)/res(j,j);
      end loop;
      res(j,j) := 1.0;
    end loop;
    return res;
  end Chebychev_Basis;

  function Orthogonal_Basis
              ( n,d : natural32; s : double_float ) return Matrix is

    res : Matrix(1..integer32(n),1..integer32(d)) := Chebychev_Basis(n,d,s);
    wrk : Matrix(1..integer32(n),1..integer32(d));
    bas : Matrix(1..integer32(n),1..integer32(n));
    qra : Standard_Floating_Vectors.Vector(1..integer32(n))
        := (1..integer32(n) => 0.0);
    pvt : Standard_Integer_Vectors.Vector(1..integer32(n))
        := (1..integer32(n) => 0); 

  begin
    wrk := res;
    QRD(wrk,qra,pvt,false);
    for i in wrk'range(1) loop
      for j in wrk'range(2) loop
        bas(i,j) := wrk(i,j);
      end loop;
      for j in integer32(d)+1..integer32(n) loop
        bas(i,j) := 0.0;
      end loop;
    end loop;
    Basis(bas,res); 
    for i in res'range(1) loop
      for j in res'range(2) loop
        res(i,j) := bas(i,j);
      end loop;
    end loop;
    return res;
  end Orthogonal_Basis;

  function Determinant
              ( mat : Matrix; rows : Standard_Integer_Vectors.Vector )
                return double_float is

  -- DESCRIPTION :
  --   Computes the determinant of the matrix obtained by selecting rows.

    res : double_float := 1.0;
    sqm : Matrix(rows'range,rows'range);
    piv : Standard_Integer_Vectors.Vector(rows'range);
    inf : integer32;

  begin
    for i in rows'range loop
      piv(i) := i;
      for j in rows'range loop
        sqm(i,j) := mat(rows(i),j);
      end loop;
    end loop;
    lufac(sqm,rows'last,piv,inf);
    for i in rows'range loop
      res := res*sqm(i,i);
    end loop;
    for i in piv'range loop
      if piv(i) > i
       then res := -res;
      end if;
    end loop;
    return res;
  end Determinant;

  procedure Maximal_Minors ( n,d : in natural32; mat : in Matrix;
                             min,max : out double_float ) is

  -- DESCRIPTION :
  --   Computes all maximal minors of a (nxd)-matrix mat, d < n.

    rows : Standard_Integer_Vectors.Vector(1..integer32(d));
    first : boolean := true;
    mindet,maxdet : double_float;

    procedure Select_Rows ( k,start : in integer32 ) is

      det : double_float;

    begin
      if k > integer32(d) then
        det := Determinant(mat,rows);
        det := abs(det);
        if first then
          mindet := det; maxdet := det; first := false;
        else
          if det > maxdet then
            maxdet := det;
          elsif det < mindet then
            mindet := det;
          end if;
        end if;
      else 
        for j in start..integer32(n) loop
          rows(k) := j;
          Select_Rows(k+1,j+1);
        end loop;
      end if;
    end Select_Rows;

  begin
    Select_Rows(1,1);
    min := mindet; max := maxdet;
  end Maximal_Minors;

  procedure Sampled_Chebychev_Basis
                  ( n,d,m : in natural32; mat : out Matrix;
                    s,ratio : out double_float ) is                   

    rs,min,max,rat,bestrat,bests : double_float;
    first : boolean;
    themat : Matrix(1..integer32(n),1..integer32(d));

  begin
    first := true;
    for i in 1..m loop
      rs := Random;
      themat := Chebychev_Basis(n,d,rs); 
      Maximal_Minors(n,d,themat,min,max);
      rat := max/min; 
      if first then
        bestrat := rat; bests := rs; first := false;
      else
        if rat < bestrat then
          bestrat := rat; bests := rs;
          mat := themat;
        end if;
      end if;
    end loop;
    ratio := bestrat;
    s := bests;
  end Sampled_Chebychev_Basis;

end Osculating_Planes;
