with Standard_Complex_Numbers;
with DoblDobl_Complex_Numbers;
with QuadDobl_Complex_Numbers;

package body Checker_Localization_Patterns is

-- Part I

  function Moving_Flag ( p : Standard_Natural_Vectors.Vector )
                       return Standard_Natural_Matrices.Matrix is

    res : Standard_Natural_Matrices.Matrix(p'range,p'range);
    col : integer32;
    free : boolean;

  begin
    for i in res'range(1) loop
      for j in res'range(2) loop
        res(i,j) := 0;
      end loop;
    end loop;
    for i in p'range loop
      col := integer32(p(p'last+1-i));
      res(i,col) := 1;
    end loop;
    for i in p'range loop
      col := integer32(p(p'last+1-i));
      for k in 1..(i-1) loop
        free := true;
        for j in 1..(col-1) loop
          if res(k,j) = 1
           then free := false;
          end if;
          exit when not free;
        end loop;
        if free
         then res(k,col) := 2;
        end if;
      end loop;
    end loop;
    return res;
  end Moving_Flag;

  function Transformation ( n,p : integer32 )
                          return Standard_Natural_Matrices.Matrix is

    res : Standard_Natural_Matrices.Matrix(1..n,1..n);

  begin
    for i in 1..n loop
      for j in 1..n loop
        res(i,j) := 0;
      end loop;
      res(i,i) := 1;
    end loop;
    res(p,p) := 2;    res(p,p+1) := 1;
    res(p+1,p) := 1;  res(p+1,p+1) := 0;
    return res;
  end Transformation;

-- Part II

  procedure Sort_Rows ( r,c : in out Standard_Natural_Vectors.Vector ) is

    min : integer32;
    tmp : natural32;

  begin
    for i in r'first..r'last-1 loop
      min := i;
      for j in i+1..r'last loop
        if r(j) < r(min)
         then min := j;
        end if;
      end loop;
      if min /= i then
        tmp := r(min);
        r(min) := r(i);
        r(i) := tmp;
        tmp := c(c'last+1-min);
        c(c'last+1-min) := c(c'last+1-i);
        c(c'last+1-i) := tmp;
      end if;
    end loop;
  end Sort_Rows;

  function Column_Pattern
             ( n,k : integer32; p,r,c : Standard_Natural_Vectors.Vector )
             return Standard_Natural_Matrices.Matrix is

    res : Standard_Natural_Matrices.Matrix(1..n,1..k);
    sr : Standard_Natural_Vectors.Vector(r'range) := r;
    sc : Standard_Natural_Vectors.Vector(r'range) := c;
    row,col : integer32;
    white_on_black : boolean;

  begin
    for i in 1..n loop        -- initialize to zero
      for j in 1..k loop
        res(i,j) := 0;
      end loop;
    end loop;
    Sort_Rows(sr,sc);
    for i in 1..k loop        -- put one in rows of white checker
      res(integer32(sr(i)),i) := 1;
    end loop;
    for i in 1..n loop        -- put two in free places
      row := integer32(p(i)); -- row of i-th black checker
      col := n-i+1;           -- column of i-th black checker
      for j in 1..k loop      -- white checker at (r(j),c(k-j+1))
        if (integer32(sr(j)) >= row) and (integer32(sc(k-j+1)) >= col) then
          if res(row,j) = 0 then
            res(row,j) := 2;
          end if;
        end if;
      end loop;
    end loop;
    for j in 1..k loop       -- zeroes where white on black
      white_on_black := false;
      for i in 1..n loop
        if sr(j) = p(i) and integer32(sc(k-j+1)) = n-i+1
         then white_on_black := true;
        end if;
        exit when white_on_black;
      end loop;
      if white_on_black then
        for j2 in j+1..k loop
          res(integer32(sr(j)),j2) := 0;
        end loop;
      end if;
    end loop;
    return res;
  end Column_Pattern;

  function Row_Pattern
             ( n,k : integer32; p,r,c : Standard_Natural_Vectors.Vector ) 
             return Standard_Natural_Matrices.Matrix is

    res : Standard_Natural_Matrices.Matrix(1..k,1..n);
    clp : constant Standard_Natural_Matrices.Matrix(1..n,1..k)
        := Column_Pattern(n,k,p,r,c);

  begin
    for i in 1..k loop
      for j in 1..n loop
        res(i,j) := clp(j,i);
      end loop;
    end loop;
    return res;
  end Row_Pattern;

  function Degree_of_Freedom
             ( lp : Standard_Natural_Matrices.Matrix ) return natural32 is

    res : natural32 := 0;

  begin
    for i in lp'range(1) loop
      for j in lp'range(2) loop
        if lp(i,j) = 2
         then res := res + 1;
        end if;
      end loop;
    end loop;
    return res;
  end Degree_of_Freedom;

  function Degree_of_Freedom
             ( p,r,c : Standard_Natural_Vectors.Vector ) return natural32 is
  begin
    return Degree_of_Freedom(Column_Pattern(p'last,r'last,p,r,c));
  end Degree_of_Freedom;

  function Rank ( lp : Standard_Natural_Matrices.Matrix;
                  i,j : integer32 ) return integer32 is

    res : integer32 := 0;

  begin
    for i1 in lp'range(1) loop
      for j1 in lp'range(2) loop
        if lp(i1,j1) = 2
         then res := res + 1;
        end if;
        if ((i1 = i) and (j1 = j)) then
          if lp(i1,j1) /= 2
           then return 0;
           else return res;
          end if;
        end if;
      end loop;
    end loop;
    return res;
  end Rank;

  function Permute_Index
              ( p,q : Standard_Natural_Matrices.Matrix ) return integer32 is

    res : integer32 := 0;
    cnt : integer32 := 0;

  begin
    for i in p'range(1) loop
      for j in p'range(2) loop
         if p(i,j) = 2 then
           cnt := cnt + 1;
           if p(i,j) /= q(i,j) then
             for k in i+1..p'last(1) loop
               if p(k,j) = 2
                then res := cnt;
               end if;
             end loop;
           end if;
         end if;
         exit when (res > 0);
      end loop;
      exit when (res > 0);
    end loop;
    return res;
  end Permute_Index;

-- PART III : mapping localization patterns to a complex matrix

  function Map ( m : Standard_Natural_Matrices.Matrix;
                 x : Standard_Complex_Vectors.Vector )
               return Standard_Complex_Matrices.Matrix is

    use Standard_Complex_Numbers;

    res : Standard_Complex_Matrices.Matrix(m'range(1),m'range(2));
    ind : integer32 := x'first-1;

  begin
    for i in m'range(1) loop
      for j in m'range(2) loop
        if m(i,j) = 0 then
          res(i,j) := Create(0.0);
        elsif m(i,j) = 1 then
          res(i,j) := Create(1.0);
        else
          ind := ind + 1;
          res(i,j) := x(ind);
        end if;
      end loop;
    end loop;
    return res;
  end Map;

  function Map ( m : Standard_Natural_Matrices.Matrix;
                 x : DoblDobl_Complex_Vectors.Vector )
               return DoblDobl_Complex_Matrices.Matrix is

    use DoblDobl_Complex_Numbers;

    res : DoblDobl_Complex_Matrices.Matrix(m'range(1),m'range(2));
    ind : integer32 := x'first-1;

  begin
    for i in m'range(1) loop
      for j in m'range(2) loop
        if m(i,j) = 0 then
          res(i,j) := Create(integer(0));
        elsif m(i,j) = 1 then
          res(i,j) := Create(integer(1));
        else
          ind := ind + 1;
          res(i,j) := x(ind);
        end if;
      end loop;
    end loop;
    return res;
  end Map;

  function Map ( m : Standard_Natural_Matrices.Matrix;
                 x : QuadDobl_Complex_Vectors.Vector )
               return QuadDobl_Complex_Matrices.Matrix is

    use QuadDobl_Complex_Numbers;

    res : QuadDobl_Complex_Matrices.Matrix(m'range(1),m'range(2));
    ind : integer32 := x'first-1;

  begin
    for i in m'range(1) loop
      for j in m'range(2) loop
        if m(i,j) = 0 then
          res(i,j) := Create(integer(0));
        elsif m(i,j) = 1 then
          res(i,j) := Create(integer(1));
        else
          ind := ind + 1;
          res(i,j) := x(ind);
        end if;
      end loop;
    end loop;
    return res;
  end Map;

  function Map ( m : Standard_Natural_Matrices.Matrix;
                 x : Standard_Complex_Matrices.Matrix )
               return Standard_Complex_Vectors.Vector is

    deg : constant natural32 := Degree_of_Freedom(m);
    res : Standard_Complex_Vectors.Vector(1..integer32(deg));
    ind : integer32 := 0;

  begin
    for i in m'range(1) loop
      for j in m'range(2) loop
        if m(i,j) = 2 then
          ind := ind + 1;
          res(ind) := x(i,j);
        end if;
      end loop;
    end loop;
    return res;
  end Map;

  function Map ( m : Standard_Natural_Matrices.Matrix;
                 x : DoblDobl_Complex_Matrices.Matrix )
               return DoblDobl_Complex_Vectors.Vector is

    deg : constant natural32 := Degree_of_Freedom(m);
    res : DoblDobl_Complex_Vectors.Vector(1..integer32(deg));
    ind : integer32 := 0;

  begin
    for i in m'range(1) loop
      for j in m'range(2) loop
        if m(i,j) = 2 then
          ind := ind + 1;
          res(ind) := x(i,j);
        end if;
      end loop;
    end loop;
    return res;
  end Map;

  function Map ( m : Standard_Natural_Matrices.Matrix;
                 x : QuadDobl_Complex_Matrices.Matrix )
               return QuadDobl_Complex_Vectors.Vector is

    deg : constant natural32 := Degree_of_Freedom(m);
    res : QuadDobl_Complex_Vectors.Vector(1..integer32(deg));
    ind : integer32 := 0;

  begin
    for i in m'range(1) loop
      for j in m'range(2) loop
        if m(i,j) = 2 then
          ind := ind + 1;
          res(ind) := x(i,j);
        end if;
      end loop;
    end loop;
    return res;
  end Map;

end Checker_Localization_Patterns;
