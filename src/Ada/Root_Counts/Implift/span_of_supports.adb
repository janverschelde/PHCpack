with text_io;                            use text_io;
with Standard_Random_Numbers;
with Standard_Common_Divisors;           use Standard_Common_Divisors;
with Standard64_Common_Divisors;         use Standard64_Common_Divisors;
with Standard_Integer_Vectors_io;        use Standard_Integer_Vectors_io;
with Standard_Integer64_Vectors;
with Standard_Integer_Matrices_io;       use Standard_Integer_Matrices_io;
with Standard_Integer64_Matrices_io;     use Standard_Integer64_Matrices_io;
with Standard_Integer_Norms;             use Standard_Integer_Norms;
with Standard_Integer_Linear_Solvers;    use Standard_Integer_Linear_Solvers;
with Standard_Integer64_Linear_Solvers;  use Standard_Integer64_Linear_Solvers;

package body Span_of_Supports is

  procedure Triangulate_Span
              ( s : in Lists_of_Integer_Vectors.List;
                n : in integer32; rank : out natural32;
                ipvt : out Standard_Integer_Vectors.Vector;
                span : out Standard_Integer_Matrices.Matrix ) is

    use Standard_Integer_Vectors;
    use Lists_of_Integer_Vectors;

    first : constant Link_to_Vector := Head_Of(s);
    tmp : List := Tail_Of(s);
    nextpt : Link_to_Vector;
    pivot,forswap,current,d,x,y : integer32;

  begin
    for i in 1..n loop
      ipvt(i) := i;
      for j in 1..n loop
        span(i,j) := 0;
      end loop;
    end loop;
    current := 1;
    while not Is_Null(tmp) loop
      nextpt := Head_Of(tmp);
      pivot := 0;
      for j in 1..n loop
        span(current,j) := nextpt(ipvt(j)) - first(ipvt(j));
        if pivot = 0 then
          if span(current,j) /= 0
           then pivot := j;         -- first nonzero column
          end if;
        end if;
      end loop;
      if pivot < current then       -- row reduction
        for j in 1..(current-1) loop
           if span(current,j) /= 0 then
             d := gcd(span(j,j),span(current,j));
             x := span(j,j)/d;
             y := span(current,j)/d;
             for k in j..n loop
               span(current,k) := x*span(current,k) - y*span(j,k);
             end loop;
           end if;
        end loop;
        pivot := 0;                 -- check if zero pivot
        for j in current..n loop
          if span(current,j) /= 0
           then pivot := j;
          end if;
          exit when (pivot /= 0);
        end loop;
      end if;
      if pivot > current then       -- swap for next row reduction
        forswap := ipvt(current);
        ipvt(current) := ipvt(pivot);
        ipvt(pivot) := forswap;
        span(current,current) := span(current,pivot);
        span(current,pivot) := 0;
        for i in 1..current-1 loop  -- also swap in previous rows
          forswap := span(i,pivot);
          span(i,pivot) := span(i,current);
          span(i,current) := forswap;
        end loop;
      end if;
      exit when (pivot /= 0 and current = n); -- full rank
      tmp := Tail_Of(tmp);
      if pivot > 0
       then current := current + 1;
      end if;
    end loop;
    if current = n and pivot /= 0
     then rank := natural32(current);     -- case of full rank
     else rank := natural32(current-1);   -- current row is zero
    end if;
  end Triangulate_Span;

  procedure Triangulate_Span
              ( s : in Lists_of_Integer_Vectors.List;
                n : in integer32; rank : out natural32;
                ipvt : out Standard_Integer_Vectors.Vector;
                span : out Standard_Integer64_Matrices.Matrix ) is

    use Standard_Integer_Vectors;
    use Lists_of_Integer_Vectors;

    first : constant Link_to_Vector := Head_Of(s);
    tmp : List := Tail_Of(s);
    nextpt : Link_to_Vector;
    pivot,current,forswap32 : integer32;
    forswap64,d,x,y : integer64;

  begin
    for i in 1..n loop
      ipvt(i) := i;
      for j in 1..n loop
        span(i,j) := 0;
      end loop;
    end loop;
    current := 1;
    while not Is_Null(tmp) loop
      nextpt := Head_Of(tmp);
      pivot := 0;
      for j in 1..n loop
        span(current,j) := integer64(nextpt(ipvt(j)) - first(ipvt(j)));
        if pivot = 0 then
          if span(current,j) /= 0
           then pivot := j;         -- first nonzero column
          end if;
        end if;
      end loop;
      if pivot < current then       -- row reduction
        for j in 1..(current-1) loop
           if span(current,j) /= 0 then
             d := gcd(span(j,j),span(current,j));
             x := span(j,j)/d;
             y := span(current,j)/d;
             for k in j..n loop
               span(current,k) := x*span(current,k) - y*span(j,k);
             end loop;
           end if;
        end loop;
        pivot := 0;                 -- check if zero pivot
        for j in current..n loop
          if span(current,j) /= 0
           then pivot := j;
          end if;
          exit when (pivot /= 0);
        end loop;
      end if;
      if pivot > current then       -- swap for next row reduction
        forswap32 := ipvt(current);
        ipvt(current) := ipvt(pivot);
        ipvt(pivot) := forswap32;
        span(current,current) := span(current,pivot);
        span(current,pivot) := 0;
        for i in 1..current-1 loop  -- also swap in previous rows
          forswap64 := span(i,pivot);
          span(i,pivot) := span(i,current);
          span(i,current) := forswap64;
        end loop;
      end if;
      exit when (pivot /= 0 and current = n); -- full rank
      tmp := Tail_Of(tmp);
      if pivot > 0
       then current := current + 1;
      end if;
    end loop;
    if current = n and pivot /= 0
     then rank := natural32(current);     -- case of full rank
     else rank := natural32(current-1);   -- current row is zero
    end if;
  end Triangulate_Span;

  function Normal_Vectors
              ( U : in Standard_Integer_Matrices.Matrix;
                n,r : in integer32 )
              return Standard_Integer_Matrices.Matrix is

    use Standard_Integer_Vectors;
    use Standard_Integer_Matrices;

    res : Matrix(1..n,1..n-r);
    B : Matrix(1..r,1..r+1);
    x : Vector(B'range(2));

  begin
    for i in 1..r loop
      for j in 1..r loop
         B(i,j) := U(i,j);
      end loop;
    end loop;
    for k in r+1..n loop
      for i in B'range(1) loop
         B(i,B'last(2)) := U(i,k);
      end loop;
      Solve0(B,x);
      Normalize(x);
      for i in 1..r loop
        res(i,k-r) := x(i);
      end loop;
      for i in r+1..n loop
        res(i,k-r) := 0;
      end loop;
      res(k,k-r) := x(r+1);
    end loop;
    return res;
  end Normal_Vectors;

  function Normal_Vectors
              ( U : in Standard_Integer64_Matrices.Matrix;
                n,r : in integer32 )
              return Standard_Integer64_Matrices.Matrix is

    use Standard_Integer64_Vectors;
    use Standard_Integer64_Matrices;

    res : Matrix(1..n,1..n-r);
    B : Matrix(1..r,1..r+1);
    x : Vector(B'range(2));

  begin
    for i in 1..r loop
      for j in 1..r loop
         B(i,j) := U(i,j);
      end loop;
    end loop;
    for k in r+1..n loop
      for i in B'range(1) loop
         B(i,B'last(2)) := U(i,k);
      end loop;
      Solve0(B,x);
      Normalize(x);
      for i in 1..r loop
        res(i,k-r) := x(i);
      end loop;
      for i in r+1..n loop
        res(i,k-r) := 0;
      end loop;
      res(k,k-r) := x(r+1);
    end loop;
    return res;
  end Normal_Vectors;

  function Apply_Pivots
              ( A : Standard_Integer_Matrices.Matrix;
                ipvt : Standard_Integer_Vectors.Vector )
              return Standard_Integer_Matrices.Matrix is

    use Standard_Integer_Matrices;

    res : Matrix(A'range(1),A'range(2));

  begin
    for i in A'range(1) loop
      for j in A'range(2) loop
        res(ipvt(i),j) := A(i,j);
      end loop;
    end loop;
    return res;
  end Apply_Pivots;

  function Apply_Pivots
              ( A : Standard_Integer64_Matrices.Matrix;
                ipvt : Standard_Integer_Vectors.Vector )
              return Standard_Integer64_Matrices.Matrix is

    use Standard_Integer64_Matrices;

    res : Matrix(A'range(1),A'range(2));

  begin
    for i in A'range(1) loop
      for j in A'range(2) loop
        res(ipvt(i),j) := A(i,j);
      end loop;
    end loop;
    return res;
  end Apply_Pivots;

  procedure Rank_of_Support
              ( s : in Lists_of_Integer_Vectors.List; n : in integer32;
                debug : in boolean; rank : out natural32;
                normals : out Standard_Integer_Matrices.Link_to_Matrix ) is

    use Standard_Integer_Vectors;
    use Standard_Integer_Matrices;
    use Lists_of_Integer_Vectors;

    r : integer32;
    vectors : Matrix(1..n,1..n);
    ipvt : Vector(1..n);

  begin
    Triangulate_Span(s,n,rank,ipvt,vectors);
    if debug then
      put("Pivoting vector : "); put(ipvt); new_line;
      put_line("Triangulated matrix of vectors : "); put(vectors);
    end if;
    r := integer32(rank);
    if r < n then
      declare
        nors : Matrix(1..n,1..n-r);
      begin
        nors := Normal_Vectors(vectors,n,r);
        if debug then
          put_line("The normal vectors : "); put(nors);
          put_line("span*normals :"); put(vectors*nors);
        end if;
        nors := Apply_Pivots(nors,ipvt);
        normals := new Matrix'(nors);
        if debug then
          declare
            supf : Matrix(1..integer32(Length_Of(s)),1..n-r);
          begin
            supf := Support_Function(s,nors);
            put_line("the support function :"); put(supf);
          end;
        end if;
      end;
    end if;
  end Rank_of_Support;

  procedure Rank_of_Support
              ( s : in Lists_of_Integer_Vectors.List; n : in integer32;
                debug : in boolean; rank : out natural32;
                normals : out Standard_Integer64_Matrices.Link_to_Matrix ) is

    use Standard_Integer_Vectors;
    use Standard_Integer64_Matrices;
    use Lists_of_Integer_Vectors;

    r : integer32;
    vectors : Matrix(1..n,1..n);
    ipvt : Vector(1..n);

  begin
    Triangulate_Span(s,n,rank,ipvt,vectors);
    if debug then
      put("Pivoting vector : "); put(ipvt); new_line;
      put_line("Triangulated matrix of vectors : "); put(vectors);
    end if;
    r := integer32(rank);
    if r < n then
      declare
        nors : Matrix(1..n,1..n-r);
      begin
        nors := Normal_Vectors(vectors,n,r);
        if debug then
          put_line("The normal vectors : "); put(nors);
          put_line("span*normals :"); put(vectors*nors);
        end if;
        nors := Apply_Pivots(nors,ipvt);
        normals := new Matrix'(nors);
        if debug then
          declare
            supf : Matrix(1..integer32(Length_Of(s)),1..n-r);
          begin
            supf := Support_Function(s,nors);
            put_line("the support function :"); put(supf);
          end;
        end if;
      end;
    end if;
  end Rank_of_Support;

  function Rank32_of_Support
              ( s : Lists_of_Integer_Vectors.List; n : integer32;
                debug : in boolean ) return natural32 is

    use Standard_Integer_Vectors;
    use Standard_Integer_Matrices;
    use Lists_of_Integer_Vectors;

    res : natural32;
    r : integer32;
    vectors : Matrix(1..n,1..n);
    ipvt : Vector(1..n);

  begin
    Triangulate_Span(s,n,res,ipvt,vectors);
    if debug then
      put("Pivoting vector : "); put(ipvt); new_line;
      put_line("Triangulated matrix of vectors : "); put(vectors);
      r := integer32(res);
      if r < n then
        declare
          nors : Matrix(1..n,1..n-r);
          supf : Matrix(1..integer32(Length_Of(s)),1..n-r);
        begin
          nors := Normal_Vectors(vectors,n,r);
          put_line("The normal vectors : "); put(nors);
          put_line("span*normals :"); put(vectors*nors);
          nors := Apply_Pivots(nors,ipvt);
          supf := Support_Function(s,nors);
          put_line("the support function :"); put(supf);
        end;
      end if;
    end if;
    return res;
  end Rank32_of_Support;

  function Rank64_of_Support
              ( s : Lists_of_Integer_Vectors.List; n : integer32;
                debug : in boolean ) return natural32 is

    use Standard_Integer_Vectors;
    use Standard_Integer64_Matrices;
    use Lists_of_Integer_Vectors;

    res : natural32;
    r : integer32;
    vectors : Matrix(1..n,1..n);
    ipvt : Vector(1..n);

  begin
    Triangulate_Span(s,n,res,ipvt,vectors);
    if debug then
      put("Pivoting vector : "); put(ipvt); new_line;
      put_line("Triangulated matrix of vectors : "); put(vectors);
      r := integer32(res);
      if r < n then
        declare
          nors : Matrix(1..n,1..n-r);
          supf : Matrix(1..integer32(Length_Of(s)),1..n-r);
        begin
          nors := Normal_Vectors(vectors,n,r);
          put_line("The normal vectors : "); put(nors);
          put_line("span*normals :"); put(vectors*nors);
          nors := Apply_Pivots(nors,ipvt);
          supf := Support_Function(s,nors);
          put_line("the support function :"); put(supf);
        end;
      end if;
    end if;
    return res;
  end Rank64_of_Support;

  function Cayley_Embedding
              ( point : Standard_Integer_Vectors.Vector; k,r : integer32 )
              return Standard_Integer_Vectors.Vector is

    res : Standard_Integer_Vectors.Vector(point'first..point'last+r);

  begin
    for i in point'range loop
      res(i) := point(i);
    end loop;
    for i in point'last+1..res'last loop
      res(i) := 0;
    end loop;
    if k > 0
     then res(point'last+k) := 1;
    end if;
    return res;
  end Cayley_Embedding;

  function Cayley_Embedding
              ( supports : Arrays_of_Integer_Vector_Lists.Array_of_Lists )
              return Lists_of_Integer_Vectors.List is

    use Standard_Integer_Vectors;
    use Lists_of_Integer_Vectors;

    res,res_last,tmp : List;
    dim : constant integer32 := supports'last - 1;
    lpt : Link_to_Vector;

  begin
    for k in supports'range loop
      tmp := supports(k);
      while not Is_Null(tmp) loop
        lpt := Head_Of(tmp);
        Append(res,res_last,Cayley_Embedding(lpt.all,k-1,dim));
        tmp := Tail_Of(tmp);
      end loop;
    end loop;
    return res;
  end Cayley_Embedding;

  function Remove_Cayley_Embedding
             ( transfo : Standard_Integer_Matrices.Matrix;
               dim : in integer32 )
             return Standard_Integer_Matrices.Matrix is

    res : Standard_Integer_Matrices.Matrix(1..dim,1..dim);

  begin
    for i in res'range(1) loop
      for j in res'range(2) loop
        res(i,j) := transfo(i,j);
      end loop;
    end loop;
    return res;
  end Remove_Cayley_Embedding;

  function Remove_Cayley_Rows
             ( mat : Standard_Integer_Matrices.Matrix;
               dim : in integer32 )
             return Standard_Integer_Matrices.Matrix is

    res : Standard_Integer_Matrices.Matrix(1..dim,mat'range(2));

  begin
    for i in res'range(1) loop
      for j in res'range(2) loop
        res(i,j) := mat(i,j);
      end loop;
    end loop;
    return res;
  end Remove_Cayley_Rows;

-- FOR TESTING PURPOSES :

  function Support_Function
              ( s : Lists_of_Integer_Vectors.List;
                normals : Standard_Integer_Matrices.Matrix )
              return Standard_Integer_Matrices.Matrix is

    use Standard_Integer_Vectors;
    use Standard_Integer_Matrices;
    use Lists_of_Integer_Vectors;

    len : constant integer32 := integer32(Length_Of(s));
    res : Matrix(1..len,normals'range(2));
    tmp : List := s;
    lpt : Link_to_Vector;

  begin
    for i in res'range(1) loop
      lpt := Head_Of(tmp);
      for j in normals'range(2) loop
        res(i,j) := 0;
        for k in normals'range(1) loop
          res(i,j) := res(i,j) + lpt(k)*normals(k,j);    
        end loop;
      end loop;
      tmp := Tail_Of(tmp);
    end loop;
    return res;
  end Support_Function;

  function Support_Function
              ( s : Lists_of_Integer_Vectors.List;
                normals : Standard_Integer64_Matrices.Matrix )
              return Standard_Integer64_Matrices.Matrix is

    use Standard_Integer_Vectors;
    use Standard_Integer64_Matrices;
    use Lists_of_Integer_Vectors;

    len : constant integer32 := integer32(Length_Of(s));
    res : Matrix(1..len,normals'range(2));
    tmp : List := s;
    lpt : Link_to_Vector;

  begin
    for i in res'range(1) loop
      lpt := Head_Of(tmp);
      for j in normals'range(2) loop
        res(i,j) := 0;
        for k in normals'range(1) loop
          res(i,j) := res(i,j) + integer64(lpt(k))*normals(k,j);    
        end loop;
      end loop;
      tmp := Tail_Of(tmp);
    end loop;
    return res;
  end Support_Function;

  function Random_Lower
              ( dim,low,upp : integer32 )
              return Standard_Integer_Matrices.Matrix is

    res : Standard_Integer_Matrices.Matrix(1..dim,1..dim);

  begin
    for i in 1..dim loop
      for j in 1..(i-1) loop
        res(i,j) := Standard_Random_Numbers.Random(low,upp);
      end loop;
      res(i,i) := 1;
      for j in (i+1)..dim loop
        res(i,j) := 0;
      end loop;
    end loop;
    return res;
  end Random_Lower;

  function Random_Support
              ( dim,nbp,low,upp,rnk : integer32 )
              return Lists_of_Integer_Vectors.List is

    use Standard_Integer_Vectors;
    use Standard_Integer_Matrices;
    use Lists_of_Integer_Vectors;

    res,res_last : List;
    point,transpoint : Vector(1..dim);
    trans : constant Matrix(1..dim,1..dim) := Random_Lower(dim,low,upp);

  begin
   -- put_line("The generated points : ");
    for i in 1..nbp loop
      for k in 1..rnk loop
        point(k) := Standard_Random_Numbers.Random(low,upp);
      end loop;
      for k in (rnk+1)..dim loop
        point(k) := 0;
      end loop;
     -- put(point); new_line;
      transpoint := trans*point;
      Append(res,res_last,transpoint);
    end loop;
   -- put_line("The coordinate transformation :"); put(trans);
    return res;
  end Random_Support;

  function Random_Tuple_of_Supports
             ( dim,mix,low,upp,rnk : integer32;
               nbp : Standard_Integer_Vectors.Vector )
             return Arrays_of_Integer_Vector_Lists.Array_of_Lists is

    use Standard_Integer_Vectors;
    use Standard_Integer_Matrices;
    use Lists_of_Integer_Vectors;
    use Arrays_of_Integer_Vector_Lists;

    res,res_last : Array_of_Lists(1..mix);
    point,transpoint : Vector(1..dim);
    trans : constant Matrix(1..dim,1..dim) := Random_Lower(dim,low,upp);

  begin
    for i in 1..mix loop
      for j in 1..nbp(i) loop
        for k in 1..rnk loop
          point(k) := Standard_Random_Numbers.Random(low,upp);
        end loop;
        for k in (rnk+1)..dim loop
          point(k) := 0;
        end loop;
        transpoint := trans*point;
        Append(res(i),res_last(i),transpoint);
      end loop;
    end loop;
    return res;
  end Random_Tuple_of_Supports;

end Span_of_Supports;
