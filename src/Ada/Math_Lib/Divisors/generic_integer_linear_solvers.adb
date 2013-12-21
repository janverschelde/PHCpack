with Standard_Natural_Numbers;           use Standard_Natural_Numbers;

package body Generic_Integer_Linear_Solvers is

  use Common_Divisors;

-- SCALERS :

  function Divisors ( a : Matrix ) return Vectors.Vector is

  -- DESCRIPTION :
  --   Returns a vector containing the gcd of the elements of each row.

    v : Vectors.Vector(a'range(1));
   
  begin
    for i in a'range(1) loop
      v(i) := a(i,a'first(2));
      if not Equal(v(i),one) then
        for j in (a'first(2)+1)..a'last(2) loop
          v(i) := gcd(v(i),a(i,j));
          exit when Equal(v(i),one);
        end loop;
      end if;
    end loop;
    return v;
  end Divisors;

  function Scale ( a : Matrix ) return Matrix is

    v : constant Vectors.Vector(a'range(1)) := Divisors(a);
    b : Matrix(a'range(1),a'range(2));

  begin
    for i in b'range(1) loop
      if (Equal(v(i),one) or Equal(v(i),zero)) then
        for j in b'range(2) loop
          b(i,j) := a(i,j);
        end loop;
      else
        for j in b'range(2) loop
          b(i,j) := a(i,j) / v(i);
        end loop;
      end if;
    end loop;
    return b;
  end Scale;

  procedure Scale ( a : in out Matrix; v : out Vectors.Vector ) is

    dv : constant Vectors.Vector(a'range(1)) := Divisors(a);

  begin
    for i in a'range(1) loop
      if (Equal(dv(i),one) or Equal(dv(i),zero)) then
        null;
      else
        for j in a'range(2) loop
          a(i,j) := a(i,j) / dv(i);
        end loop;
      end if;
    end loop;
    v := dv;
  end Scale;

  procedure Scale ( a : in out Matrix ) is

    v : constant Vectors.Vector(a'range(1)) := Divisors(a);

  begin
    for i in a'range(1) loop
      if (Equal(v(i),one) or Equal(v(i),zero)) then
        null;
      else
        for j in a'range(2) loop
          a(i,j) := a(i,j) / v(i);
        end loop;
      end if;
    end loop;
  end Scale;

  procedure Scale ( a : in out Matrix; row,col : in integer32 ) is

    g : number := a(row,col);

  begin
    if not Equal(g,one) then
      for l in (col+1)..a'last(2) loop
        g := gcd(g,a(row,l));
        exit when Equal(g,one);
      end loop;
    end if;
    if (not Equal(g,zero)) and (not Equal(g,one)) then
      for l in row..a'last(2) loop
        a(row,l) := a(row,l)/g;
      end loop;
    end if;
  end Scale;

-- STATIC TRIANGULATORS :

  procedure Upper_Triangulate ( l : out Matrix; a : in out Matrix ) is

    row,pivot : integer32;
    temp,aa,bb,ka,lb,d,a_rowk,a_ik : number;
    ll : Matrix(a'range(1),a'range(1));

  begin
    for i in ll'range(1) loop
      for j in ll'range(2) loop
        ll(i,j) := Create(0);
      end loop;
      ll(i,i) := Create(1);
    end loop;
    row := a'first(1);
    for j in a'first(2)..a'last(2) loop
      pivot := row-1;                                         -- find pivot
      for i in row..a'last(1) loop
        if not Equal(a(i,j),zero)
         then pivot := i; exit;
        end if;
      end loop;
      if pivot >= row then
        if pivot /= row then                    -- interchange if necessary
          for k in a'range(2) loop
            temp := Create(0);
            Copy(a(row,k),temp);
            Copy(a(pivot,k),a(row,k));
            Copy(temp,a(pivot,k));
          end loop;
          for k in ll'range(2) loop
            Copy(ll(row,k),temp);
            Copy(ll(pivot,k),ll(row,k));
            Copy(temp,ll(pivot,k));
          end loop;
        end if;
        for i in (row+1)..a'last(1) loop                  -- make zeroes
          if not Equal(a(i,j),zero) then
            Clear(ka); Clear(lb); Clear(d);
            gcd(a(row,j),a(i,j),ka,lb,d);
            aa := a(row,j)/d;     bb := a(i,j)/d;
            if Equal(aa,bb) and then Equal(ka,zero)
             then Copy(lb,ka); Copy(zero,lb);
            end if;
            if Equal(aa,-bb) and then Equal(ka,zero)
             then ka := -lb; Copy(zero,lb);
            end if;
            for k in j..a'last(2) loop
              a_rowk := Create(0);
              Copy(a(row,k),a_rowk); Clear(a(row,k));
              a_ik := Create(0);
              Copy(a(i,k),a_ik); Clear(a(i,k));
              a(row,k) :=  ka*a_rowk + lb*a_ik;
              a(i,k)   := (-bb)*a_rowk + aa*a_ik;
            end loop;
            for k in ll'range(2) loop
              a_rowk := ll(row,k);
              a_ik := ll(i,k);
              ll(row,k) :=  ka*a_rowk + lb*a_ik;
              ll(i,k)   := (-bb)*a_rowk + aa*a_ik;
            end loop;
          end if;
        end loop;
        row := row + 1;
      end if;
      exit when row > a'last(1);
    end loop;
    l := ll;
  end Upper_Triangulate;

  procedure Upper_Triangulate ( a : in out Matrix ) is

    row,pivot : integer32;
    temp,aa,bb,ka,lb,d,a_rowk,a_ik : number;

  begin
    row := a'first(1);
    for j in a'range(2) loop
     -- put("in column j = "); put(j,1); new_line;
      pivot := row-1;                                           -- find pivot
      for i in row..a'last(1) loop
        if not Equal(a(i,j),zero)
         then pivot := i; exit;
        end if;
      end loop;
      if pivot >= row then
        if pivot /= row then                      -- interchange if necessary 
          for k in a'range(2) loop
           -- temp := Create(0);
	    Copy(a(row,k),temp);
	    Copy(a(pivot,k),a(row,k));
	    Copy(temp,a(pivot,k)); Clear(temp);
          end loop;
        end if;
        for i in (row+1)..a'last(1) loop                    -- make zeroes
          if not Equal(a(i,j),zero) then
            Clear(ka); Clear(lb); Clear(d);
           -- put("making i = "); put(i,1); put_line(" zero...");
            gcd(a(row,j),a(i,j),ka,lb,d);
           -- if Equal(d,zero)
           --  then put_line("d is zero!");
           -- end if;
            aa := a(row,j)/d;     bb := a(i,j)/d;
            if Equal(aa,bb) and then Equal(ka,zero)
             then Copy(lb,ka); Copy(zero,lb);
            end if;
            if Equal(aa,-bb) and then Equal(ka,zero)
             then ka := -lb; Copy(zero,lb);
            end if;
           -- Copy(zero,a_rowk); -- a_rowk := Create(0);
            for k in j..a'last(2) loop
              Copy(a(row,k),a_rowk); Clear(a(row,k));
             -- Copy(zero,a_ik); -- a_ik := Create(0);
              Copy(a(i,k),a_ik); Clear(a(i,k));
              a(row,k) :=  ka*a_rowk + lb*a_ik;
              if Equal(bb,zero)
               then a(i,k) := aa*a_ik;
               else a(i,k) := (-bb)*a_rowk + aa*a_ik;
              end if;
            end loop;
          end if;
        end loop;
        row := row + 1;
      end if;
      exit when row > a'last(1);
    end loop;
 -- exception
 --   when others => put_line("exception in Upper_Triangulate");
 --                  raise;
  end Upper_Triangulate;

  procedure Upper_Triangulate 
                 ( a : in out Matrix;
                   ipvt : in out Standard_Integer_Vectors.Vector ) is

    row,pivot,tmppiv : integer32;
    temp,aa,bb,ka,lb,d,a_rowk,a_ik : number;

  begin
    row := a'first(1);
    for j in a'first(2)..a'last(2) loop
      pivot := row-1;                                           -- find pivot
      for i in row..a'last(1) loop
        if not Equal(a(i,j),zero)
         then pivot := i; exit;
        end if;
      end loop;
      if pivot >= row then
        if pivot /= row then                  -- interchange if necessary
          temp := Create(0);
          for k in a'range(2) loop
            Copy(a(row,k),temp);
            Copy(a(pivot,k),a(row,k));
            Copy(temp,a(pivot,k));
	  end loop;
          tmppiv := ipvt(row);
          ipvt(row) := ipvt(pivot);
          ipvt(pivot) := tmppiv;
        end if;
        for i in (row+1)..a'last(1) loop                   -- make zeroes
          if not Equal(a(i,j),zero) then
            gcd(a(row,j),a(i,j),ka,lb,d);
            aa := a(row,j)/d;     bb := a(i,j)/d;
            if Equal(aa,bb) and then Equal(ka,zero)
             then Copy(lb,ka); Copy(zero,lb);
            end if;
            if Equal(aa,-bb) and then Equal(ka,zero)
             then ka := -lb; Copy(zero,lb);
            end if;
            for k in j..a'last(2) loop
              a_rowk := Create(0);
              Copy(a(row,k),a_rowk); Clear(a(row,k));
              a_ik := Create(0);
              Copy(a(i,k),a_ik); Clear(a(i,k)); 
              a(row,k) :=  ka*a_rowk + lb*a_ik;
              a(i,k)   := (-bb)*a_rowk + aa*a_ik;
            end loop;
          end if;
        end loop;
        row := row + 1;
      end if;
      exit when row > a'last(1);
    end loop;
  end Upper_Triangulate;

  procedure Lower_Triangulate ( a : in out Matrix; u : out Matrix ) is

    column,pivot : integer32;
    temp,aa,bb,ka,lb,d,a_kcolumn,a_kj : number;
    uu : Matrix(a'range(2),a'range(2));

  begin
    for i in uu'range(1) loop
      for j in uu'range(2) loop
        uu(i,j) := Create(0);
      end loop;
      uu(i,i) := Create(1);
    end loop;
    column := a'first(2);
    for i in a'first(1)..a'last(1) loop
      pivot := column-1;                                      -- find pivot
      for j in column..a'last(2) loop
        if not Equal(a(i,j),zero)
         then pivot := j; exit;
        end if;
      end loop;
      if pivot >= column then
        if pivot /= column then                 -- interchange if necessary
          temp := Create(0);
          for k in a'range(1) loop
            Copy(a(k,column),temp);
            Copy(a(k,pivot),a(k,column));
            Copy(temp,a(k,pivot));
          end loop;
          for k in uu'range(1) loop
            Copy(uu(k,column),temp);
            Copy(uu(k,pivot),uu(k,column));
            Copy(temp,uu(k,pivot));
          end loop;
        end if;
        for j in (column+1)..a'last(2) loop               -- make zeroes
          if not Equal(a(i,j),zero) then
            gcd(a(i,column),a(i,j),ka,lb,d);
            aa := a(i,column)/d;        bb := a(i,j)/d;
            if Equal(aa,bb) and then Equal(ka,zero)
             then Copy(lb,ka); Copy(zero,lb);
            end if;
            if Equal(aa,-bb) and then Equal(ka,zero)
             then ka := -lb; Copy(zero,lb);
            end if;
            for k in i..a'last(1) loop
              a_kcolumn := Create(0);
              Copy(a(k,column),a_kcolumn); Clear(a(k,column));
              a_kj := Create(0);
              Copy(a(k,j),a_kj); Clear(a(k,j)); -- a_kj := a(k,j);
              a(k,column) := a_kcolumn*ka    + a_kj*lb;
              a(k,j)      := a_kcolumn*(-bb) + a_kj*aa;
            end loop;
            for k in uu'range(1) loop
              Copy(uu(k,column),a_kcolumn); Clear(uu(k,column));
              Copy(uu(k,j),a_kj); Clear(uu(k,j));
              uu(k,column) := a_kcolumn*ka    + a_kj*lb;
              uu(k,j)      := a_kcolumn*(-bb) + a_kj*aa;
            end loop;
          end if;
        end loop;
        column := column + 1;
      end if;
      exit when column > a'last(2);
    end loop;
    u := uu;
  end Lower_Triangulate;

  procedure Lower_Triangulate ( a : in out Matrix ) is

    column,pivot : integer32;
    temp,aa,bb,ka,lb,d,a_kcolumn,a_kj : number;

  begin
    column := a'first(2);
    for i in a'first(1)..a'last(1) loop
      pivot := column-1;                                      -- find pivot
      for j in column..a'last(2) loop
        if not Equal(a(i,j),zero)
         then pivot := j; exit;
        end if;
      end loop;
      if pivot >= column then
        if pivot /= column then                 -- interchange if necessary
          temp := Create(0);
          for k in a'range(1) loop
            Copy(a(k,column),temp);
            Copy(a(k,pivot),a(k,column));
            Copy(temp,a(k,pivot));
          end loop;
        end if;
        for j in (column+1)..a'last(2) loop               -- make zeroes
          if not Equal(a(i,j),zero) then
            gcd(a(i,column),a(i,j),ka,lb,d);
            aa := a(i,column)/d; bb := a(i,j)/d;
            if Equal(aa,bb) and then Equal(ka,zero)
             then Copy(lb,ka); Copy(zero,lb);
            end if;
            if Equal(aa,-bb) and then Equal(ka,zero)
             then ka := -lb; Copy(zero,lb);
            end if;
            for k in i..a'last(1) loop
              a_kcolumn := Create(0);
              Copy(a(k,column),a_kcolumn); Clear(a(k,column));
              a_kj := Create(0);
              Copy(a(k,j),a_kj); Clear(a(k,j));
              a(k,column) := a_kcolumn*ka    + a_kj*lb;
              a(k,j)      := a_kcolumn*(-bb) + a_kj*aa;
            end loop;
          end if;
        end loop;
        column := column + 1;
      end if;
      exit when column > a'last(2);
    end loop;
  end Lower_Triangulate;

  procedure Lower_Triangulate
               ( a : in out Matrix;
                 ipvt : in out Standard_Integer_Vectors.Vector ) is

    column,pivot,tmppiv : integer32;
    temp,aa,bb,ka,lb,d,a_kcolumn,a_kj : number;

  begin
    column := a'first(2);
    for i in a'first(1)..a'last(1) loop
      pivot := column-1;                                       -- find pivot
      for j in column..a'last(2) loop
        if not Equal(a(i,j),zero)
         then pivot := j; exit;
        end if;
      end loop;
      if pivot >= column then
        if pivot /= column then                  -- interchange if necessary
          temp := Create(0);
          for k in a'range(1) loop
            Copy(a(k,column),temp);
            Copy(a(k,pivot),a(k,column));
            Copy(temp,a(k,pivot));
          end loop;
          tmppiv := ipvt(column);
          ipvt(column) := ipvt(pivot);
          ipvt(pivot) := tmppiv;
        end if;
        for j in (column+1)..a'last(2) loop               -- make zeroes
          if not Equal(a(i,j),zero) then
            gcd(a(i,column),a(i,j),ka,lb,d);
            aa := a(i,column)/d; bb := a(i,j)/d;
            if Equal(aa,bb) and then Equal(ka,zero)
             then Copy(lb,ka); Copy(zero,lb);
            end if;
            if Equal(aa,-bb) and then Equal(ka,zero)
             then ka := -lb; Copy(zero,lb);
            end if;
            for k in i..a'last(1) loop
              a_kcolumn := Create(0);
              Copy(a(k,column),a_kcolumn); Clear(a(k,column));
              a_kj := Create(0);
              Copy(a(k,j),a_kj); Clear(a(k,j));
              a(k,column) := a_kcolumn*ka    + a_kj*lb;
              a(k,j)      := a_kcolumn*(-bb) + a_kj*aa;
            end loop;
          end if;
        end loop;
        column := column + 1;
      end if;
      exit when column > a'last(2);
    end loop;
  end Lower_Triangulate;

-- SELECTORS :

  function Det ( a : Matrix ) return number is

  -- NOTE :
  --   The triangulation is implemented independently to keep track
  --   of row interchanges.

    res : number := one;
    m : matrix(a'range(1),a'range(2));
    row,pivot : integer32;
    temp,aa,bb,ka,lb,d,m_rowk,m_ik : number;

  begin
    Copy(a,m);   -- triangulate m
    row := m'first(1);
    for j in m'first(2)..m'last(2) loop
      pivot := row-1;                                           -- find pivot
      for i in row..m'last(1) loop
        if not Equal(m(i,j),zero)
         then pivot := i; exit;
        end if;
      end loop;
      if pivot >= row then
        if pivot /= row then                      -- interchange if necessary
          temp := Create(0);
          for k in m'range(2) loop
            Copy(m(row,k),temp);
            Copy(m(pivot,k),m(row,k));
            Copy(temp,m(pivot,k));
          end loop;
          Min(res);
        end if;
        for i in (row+1)..m'last(1) loop                   -- make zeroes
          if not Equal(m(i,j),zero) then
            gcd(m(row,j),m(i,j),ka,lb,d);
            aa := m(row,j)/d; bb := m(i,j)/d;
            if Equal(aa,bb) and then Equal(ka,zero)
             then Copy(lb,ka); Copy(zero,lb);
            end if;
            if Equal(aa,-bb) and then Equal(ka,zero)
             then ka := -lb; Copy(zero,lb);
            end if;
            for k in j..m'last(2) loop
              m_rowk := Create(0);
              Copy(m(row,k),m_rowk); Clear(m(row,k));
              m_ik := Create(0);
              Copy(m(i,k),m_ik); Clear(m(i,k));
              m(row,k) :=  ka*m_rowk + lb*m_ik;
              m(i,k)   := (-bb)*m_rowk + aa*m_ik;
            end loop;
          end if;
        end loop;
        row := row + 1;
      end if;
      exit when row > m'last(1);
    end loop;
    for k in m'range(1) loop
      Mul(res,m(k,k));
    end loop;
    Clear(m);
    return res;
  end Det;

  function Per ( i : integer32; n : natural32; a : matrix;
                 kk : Standard_Integer_Vectors.Vector ) return number is
  begin
    if i = integer32(n+1) then
      return one;
    else
      declare
        res,acc : number;
        kkk : Standard_Integer_Vectors.Vector(kk'range) := kk;
      begin
        res := Create(0);
        for j in kk'range loop
          if not Equal(a(i,j),zero) and then (kk(j) /= 0) then
            kkk(j) := kkk(j) - 1;
            acc := a(i,j)*Per(i+1,n,a,kkk);
            Add(res,acc);
            Clear(acc);
            kkk(j) := kkk(j) + 1;
          end if;
        end loop;
        return res;
      end;
    end if;
  end Per;

  function Per ( i : integer32; n : natural32; a : matrix;
                 kk : Standard_Integer_Vectors.Vector; max : number )
               return number is
  begin
    if i = integer32(n+1) then
      return one;
    else
      declare
        res,acc : number;
        kkk : Standard_Integer_Vectors.Vector(kk'range) := kk;
      begin
        res := Create(0);
        for j in kk'range loop
          if not Equal(a(i,j),zero) and then (kk(j) /= 0) then
            kkk(j) := kkk(j) - 1;
            acc := a(i,j)*Per(i+1,n,a,kkk,max);
            Add(res,acc);
            Clear(acc);
            kkk(j) := kkk(j) + 1;
          end if;
          exit when (res > max) or Equal(res,max);
        end loop;
        return res;
      end;
    end if;
  end Per;

  function Per ( a : matrix; k : Standard_Integer_Vectors.Vector )
               return number is

  -- ALGORITHM :
  --   Row expansion without memory, as developed by C.W. Wampler,
  --   see `Bezout Number Calculations for Multi-Homogeneous Polynomial
  --        Systems', Appl. Math. Comput. 51:(2-3), 143-157, 1992.

  begin
    return Per(1,natural32(a'last(1)),a,k);
  end Per;

  function Per ( a : matrix; k : Standard_Integer_Vectors.Vector;
                 max : number ) return number is

  -- ALGORITHM :
  --   Row expansion without memory, as developed by C.W. Wampler,
  --   see `Bezout Number Calculations for Multi-Homogeneous Polynomial
  --        Systems', Appl. Math. Comput. 51:(2-3), 143-157, 1992.

  begin
    return Per(1,natural32(a'last(1)),a,k,max);
  end Per;

  function Rank ( a : Matrix ) return integer32 is

    res : integer32 := 0; 
    m : Matrix(a'range(1),a'range(2));
    column : integer32;

  begin
    Copy(a,m);
    Upper_Triangulate(m);
   -- compute the length of chain of nonzero elements in m :
   -- search first nonzero element in first row of m :
    column := m'first(2)-1;
    for k in m'range(2) loop
      if not Equal(m(m'first(1),k),zero)
       then column := k;
      end if;
      exit when (column = k);
    end loop;
    if column < m'first(2) then
      return 0;   -- all elements of m are zero
    else
      for k in m'range(1) loop
        exit when column > m'last(2);
        if not Equal(m(k,column),zero) then
          res := res + 1;
        else -- search for next nonzero element on row k :
          for l in column+1..m'last(2) loop
            if not Equal(m(k,l),zero)
             then column := l; res := res + 1;
            end if;
            exit when (column = l);
          end loop;
        end if;
        column := column + 1;
      end loop;
    end if;
    Clear(m);
    return res;
  end Rank;

-- DYNAMIC TRIANGULATOR :

  procedure Triangulate ( l : in integer32; m : in out matrix;
                          ipvt : in out Standard_Integer_Vectors.vector;
                          piv : out integer32 ) is

  -- DESCRIPTION :
  --   Updates lth row of m such that m remains upper triangular.

    tmp,a,b,lcmab,faca,facb : number;
    pivot,index,tmppiv : integer32;

  begin
    Switch(ipvt,l,m);                      -- pivoting for previous unknowns
    index := 1;                         -- update : make l-1 zeroes in row l
    for k in 1..(l-1) loop
      if (not Equal(m(l,index),zero))
         and then (not Equal(m(k,index),zero)) then   -- make m(l,index) zero
        a := m(k,index); b := m(l,index);
        lcmab := lcm(a,b);
        if lcmab < zero then lcmab := -lcmab; end if;
        facb := lcmab/b; faca := lcmab/a;
        if facb > zero then
          for i in index..m'last(2) loop
            m(l,i) :=  facb*m(l,i) - faca*m(k,i);
          end loop;
        else
          for i in index..m'last(2) loop
            m(l,i) := (-facb)*m(l,i) + faca*m(k,i);
          end loop;
        end if;
      end if;
      if not Equal(m(k,index),zero)
       then index := index + 1;
      end if;
    end loop;
    pivot := 0;                                              -- search pivot
    for k in l..m'last(2)-1 loop
      if not Equal(m(l,k),zero)
       then pivot := k;
      end if;
      exit when pivot /= 0;
    end loop;
    if pivot > l then
      tmp := Create(0);
      for k in 1..l loop              -- interchange columns l and pivot
        Copy(m(k,l),tmp);
        Copy(m(k,pivot),m(k,l));
        Copy(tmp,m(k,pivot));
      end loop;
      tmppiv := ipvt(l);
      ipvt(l) := ipvt(pivot);
      ipvt(pivot) := tmppiv;
    end if;
    piv := pivot;
  end Triangulate;

  procedure Switch ( ipvt : in Standard_Integer_Vectors.vector;
                     index : in integer32; m : in out matrix ) is

    tmpv : Vectors.Vector(m'range(2));

  begin
    for k in tmpv'range loop
      tmpv(k) := Create(0);
      Copy(m(index,k),tmpv(k));
    end loop;
    for k in tmpv'range loop
      Copy(tmpv(ipvt(k)),m(index,k));
    end loop;
    Vectors.Clear(tmpv);
  end Switch;

  procedure Switch ( ipvt : in Standard_Integer_Vectors.vector;
                     first,last : in integer32; m : in out matrix) is

    tmpv : Vectors.Vector(m'range(2));

  begin
    for index in first..last loop
      for k in tmpv'range loop
        tmpv(k) := Create(0);
        Copy(m(index,k),tmpv(k));
      end loop;
      for k in tmpv'range loop
        Copy(tmpv(ipvt(k)),m(index,k));
      end loop;
      Vectors.Clear(tmpv);
    end loop;
  end Switch;

  procedure Switch ( l,pivot,index : in integer32; m : in out matrix ) is

    tmp : number;

  begin
    if l /= pivot then
      tmp := Create(0); Copy(m(index,l),tmp);
      Copy(m(index,pivot),m(index,l));
      Copy(tmp,m(index,pivot));
    end if;
  end Switch;

  procedure Switch ( l,pivot : in integer32;
                     first,last : in integer32; m : in out matrix ) is

    tmp : number;

  begin
    if l /= pivot then
      for index in first..last loop
        tmp := Create(0); Copy(m(index,l),tmp);
        Copy(m(index,pivot),m(index,l));
        Copy(tmp,m(index,pivot));
      end loop;
    end if;
  end Switch;

-- SOLVERS :

--  function Check0 ( a : Matrix; x : Vectors.Vector ) return boolean is
--   
--   -- DESCRIPTION :
--   --   Returns true if x is a solution of the system a*x = 0.
--
--    tmp : Vectors.Vector(a'range(1));
--
--  begin
--    tmp := a*x;
--    for i in tmp'range loop
--      if not Equal(tmp(i),zero)
--       then Vectors.Clear(tmp); return false;
--      end if;
--    end loop;
--    Vectors.Clear(tmp);
--    return true;
--  end Check0;

  procedure Solve0 ( a : in Matrix; x : in out Vectors.Vector ) is

  -- ALGORITHM :
  --   An intermediate, generating matrix tmp will be constructed,
  --   such that 
  --      1) the solution x to tmp*x = 0 is the same of a*x = 0;
  --      2) tmp(i,i) and tmp(i,ind) are the only nonzero entries.
  --   Before this construction, it will be checked whether there
  --   exists a zero column and the index ind must be determined.
  --   After the definition of tmp, the back substitution process
  --   yields a solution.

    piv,ind : integer32;
    tmp : Matrix(a'first(1)..(a'last(1)+1),a'range(2));
    res : Vectors.Vector(x'range);
    pivots : Standard_Integer_Vectors.Vector(x'range);
    zero_column : boolean;

  begin
   -- initialization of tmp,ind and pivots :
    for i in tmp'range(1) loop
      for j in tmp'range(2) loop
	tmp(i,j) := Create(0);
      end loop;
    end loop;
    for i in pivots'range loop
      pivots(i) := i;
    end loop;
    ind := x'first(1)-1;
    for i in a'range(1) loop
      piv := pivots'first-1;
      for j in a'range(2) loop
	if not Equal(a(i,j),zero) then
          piv := pivots(j);
	  pivots(j) := pivots(i);
	  pivots(i) := piv;
	  exit;
        end if;
      end loop;
      zero_column := true;
      for j in a'first(1)..i loop
        Copy(a(j,pivots(i)),tmp(j,i));
	if zero_column and then not Equal(tmp(j,i),zero)
	 then zero_column := false;
        end if;
      end loop;
      if piv < pivots'first or else zero_column or else Equal(tmp(i,i),zero)
       then ind := i; exit;
      end if;
    end loop;
    if zero_column then
      for i in x'range loop
        Copy(zero,x(i));
      end loop;
      Copy(one,x(ind));
    else
      -- if ((ind < x'first(1)) and (a'last(1) >= a'last(2)));
      -- does no longer work for gnat15p under windows...
       if ((ind < x'first(1)) and (a'length(1) >= a'length(2))) then
         for i in x'range loop
           Copy(zero,x(i));
         end loop;
       else
         if ind < x'first(1) then
           ind := a'last(1)+1;
	   for j in tmp'range(2) loop
	     Copy(zero,tmp(ind,j));
           end loop;
           zero_column := true;
           for j in a'range(1) loop
             Copy(a(j,pivots(ind)),tmp(j,ind));
             if zero_column and then not Equal(tmp(j,ind),zero)
              then zero_column := false;
             end if;
             end loop;
           end if;
           if zero_column then
             for i in x'range loop
               Copy(zero,x(i));
             end loop;
	     Copy(one,x(ind));
           else -- construct generating matrix :
              for i in reverse (tmp'first(2)+1)..(ind-1) loop  -- i = column
                for k in tmp'first(1)..(i-1) loop
                  if not Equal(tmp(k,i),zero) then -- make tmp(k,i) zero
                    declare 
                      aa,bb,d : number;
                    begin
                      d := gcd(tmp(i,i),tmp(k,i));
                      aa := tmp(i,i)/d;     bb := tmp(k,i)/d;
                      for l in k..(i-1) loop
                        Mul(tmp(k,l),aa);
                      end loop;
                      Copy(zero,tmp(k,i)); --aa*tmp(k,i) - bb*tmp(i,i);
                      tmp(k,ind) := aa*tmp(k,ind) - bb*tmp(i,ind);
                    end;
                    Scale(tmp,k,k);  -- to avoid numeric_error
                  end if;
                end loop; -- upper half of ith colum consists of zero entries
              end loop;
             -- generate x by back substitution :
              x(ind) := tmp(x'first,x'first);
              for i in (x'first+1)..(ind-1) loop
	        if not Equal(tmp(i,i),zero)
                 then x(ind) := lcm(tmp(i,i),x(ind));
                end if;
              end loop;
              for i in x'first..(ind-1) loop
	        if Equal(tmp(i,i),zero)
	         then Copy(zero,x(i));
                 else x(i) := (-(tmp(i,ind)*x(ind)))/tmp(i,i);
                end if;
              end loop;
           end if;
      end if;
    end if;
    for i in res'range loop                      -- take pivots into account
      res(i) := Create(0);
    end loop;
    if ind = 0 then
      for i in x'range loop
        Copy(x(i),res(pivots(i)));
      end loop;
    else
      for i in x'first..ind loop
        Copy(x(i),res(pivots(i)));
      end loop;
    end if;
    Vectors.Copy(res,x);
  end Solve0;

  procedure Solve1 ( a : in Matrix; x : in out Vectors.Vector;
                     b : in Vectors.Vector; fail : out boolean ) is

    acc : number;

  begin
    fail := false;
    for i in reverse x'range loop
      Copy(b(i),x(i));
      for j in (i+1)..x'last loop
        acc := a(i,j)*x(j);
        Sub(x(i),acc);
        Clear(acc);
      end loop;
      if not Equal(x(i),zero) and then not Equal(a(i,i),zero) then
        acc := Rmd(x(i),a(i,i));
        if Equal(acc,zero)
         then Div(x(i),a(i,i));
         else fail := true;
        end if;
        Clear(acc);
        if fail
         then Vectors.Clear(x); return;
        end if;
      end if;
    end loop;
  end Solve1;

  procedure Solve1 ( a : in Matrix; b : in out Vectors.Vector;
                     fail : out boolean ) is

    acc : number;

  begin
    fail := false;
    for i in reverse b'range loop
      for j in (i+1)..b'last loop
        acc := a(i,j)*b(j);
        Sub(b(i),acc);
        Clear(acc);
      end loop;
      if not Equal(b(i),zero) and then not Equal(a(i,i),zero) then
        acc := Rmd(b(i),a(i,i));
        if Equal(acc,zero)
         then Div(b(i),a(i,i));
         else fail := true; 
        end if;
        Clear(acc);
        if fail
         then Vectors.Clear(b); return;
        end if;
      end if;
    end loop;
  end Solve1;

  function Solve ( m : Matrix; ipvt : Standard_Integer_Vectors.Vector )
                 return Vectors.Vector is

    x,res : Vectors.Vector(ipvt'range);
    a : matrix(m'first(1)..m'last(1)-1,m'range(2));
    ind : integer32 := a'first(1);       -- index for the current row number
    cnt0 : natural32 := 0;               -- counts the number of zero rows

  begin
    a(a'first(1),a'first(2)) := Create(0);
    for k in a'range(1) loop
      if not Equal(m(k,k),zero) then        -- otherwise : skip zero row !
        for l in a'range(2) loop
          Copy(m(k,l),a(ind,l));
        end loop;
        ind := ind + 1;
      else
        for l in a'range(2) loop
          Copy(m(k,l),a(a'last(1) - integer32(cnt0),l));
        end loop;
        cnt0 := cnt0 + 1;
      end if;
    end loop;
    for i in x'range loop
      x(i) := Create(0);
    end loop;
    Solve0(a,x);
    for k in res'range loop
      res(ipvt(k)) := x(k);
    end loop;
    if res(res'last) < zero
     then Vectors.Min(res);
    end if;
    Clear(a);
    return res;
  end Solve;

end Generic_Integer_Linear_Solvers;
