with text_io;                            use text_io;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Integer_Vectors;           use Standard_Integer_Vectors;
with Standard_Integer_Vectors_io;        use Standard_Integer_Vectors_io;
with Standard_Integer_Matrices;          use Standard_Integer_Matrices;
with Standard_Integer_Matrices_io;       use Standard_Integer_Matrices_io;
with Standard_Integer_Linear_Solvers;    use Standard_Integer_Linear_Solvers;
with Lists_of_Integer_Vectors;           use Lists_of_Integer_Vectors;
with Arrays_of_Integer_Vector_Lists;     use Arrays_of_Integer_Vector_Lists;
with Arrays_of_Integer_Vector_Lists_io;  use Arrays_of_Integer_Vector_Lists_io;
with Standard_Complex_Laurentials;       use Standard_Complex_Laurentials;
with Standard_Complex_Laur_Systems;      use Standard_Complex_Laur_Systems;
with Standard_Complex_Laur_Systems_io;   use Standard_Complex_Laur_Systems_io;
with Supports_of_Polynomial_Systems;

procedure ts_supspan is

-- DESCRIPTION :
--   Computes the span of the supports of a polynomial system.
--   If the span is less than the dimension of the ambient space,
--   then the generators of the normals to the hyperplanes that
--   contain the supports are computed.

  function First_Vectors
             ( s : Array_of_Lists; n : integer32 ) return Matrix is

  -- DESCRIPTION :
  --   Returns a matrix with in its rows n-vectors defined by the
  --   differences of the first two vectors in s.

    res : Matrix(s'range,1..n);
    first,second : Link_to_Vector;

  begin
    for i in s'range loop
       if Length_Of(s(i)) < 2 then
         for j in 1..n loop
           res(i,j) := 0;
         end loop;
       else
         first := Head_Of(s(i));
         second := Head_Of(Tail_Of(s(i)));
         for j in 1..n loop
           res(i,j) := second(j) - first(j);
         end loop;
       end if;
    end loop;
    return res;
  end First_Vectors;

  procedure Search_Point
              ( s : in List; v : in Vector;
                k : in out integer32; w : out Link_to_Vector ) is

  -- DESCRIPTION :
  --   Searches the list s, starting at the k-th point k for the first point 
  --   for which the vector w defined as difference with the first point in s
  --   makes a nonzero inner product with v.
  --   If no such point exists, then k is zero on return,
  --   otherwise k is the position of the point in the list s.

  -- REQUIRED : Length_Of(s) >= k > 1.

    first : constant Link_to_Vector := Head_Of(s);
    next : Link_to_Vector;
    ptr : List := s;
    wrk : Vector(v'range);

  begin
    for i in 1..(k-1) loop
      ptr := Tail_Of(ptr);
    end loop;
    while not Is_Null(ptr) loop
      next := Head_Of(ptr);
      for i in v'range loop
        wrk(i) := next(i) - first(i);
      end loop;
      if v*wrk /= 0 
       then w := new Vector'(wrk); return;
      end if;
      ptr := Tail_Of(ptr); k := k + 1;
    end loop;
    k := 0;
  end Search_Point;

  procedure Span_of_System ( p : in Laur_Sys; n : in integer32 ) is

  -- DESCRIPTION :
  --   Computes the span of the supports of p.
  --   The ambient dimension is n.

    s : Array_of_Lists(p'range) := Supports_of_Polynomial_Systems.Create(p);
    F : Matrix(1..n,p'range) := First_Vectors(s,n);
    A : Matrix(1..n,p'range) := F;
    v : Vector(p'range);
    r,ind,pos : integer32;
    w : Link_to_Vector;

  begin
    new_line;
    put_line("The supports : "); put(s);
    put_line("The matrix of first vectors : "); put(A);
    Upper_Triangulate(A);
    put_line("After upper triangulation : "); put(A);
    r := Rank(A);
    put("The rank of A : "); put(r,1); put_line(".");
    if r = n then
      put_line("The supports span the ambient space.");
    else
      for k in 1..n loop
        Scale(A);
        Solve0(A,v);
        put("A vector in the kernel of A : "); put(v); new_line;
        ind := s'first;
        while ind <= s'last loop
          pos := 2;
          Search_Point(s(ind),v,pos,w);
          exit when (pos /= 0);
          ind := ind + 1;
        end loop;
        if pos = 0 then
          put_line("Supports do not span the ambient space ?"); exit;
        else
          put("Found nonparallel vector at position "); put(pos,1);
          put(" of support "); put(ind,1); 
          put(" : "); put(w.all); put_line(".");
          for i in A'range(1) loop
            if ind /= i then
              for j in 1..n loop
                A(i,j) := F(i,j);
              end loop;
            else
              for j in 1..n loop
                A(i,j) := w(j);
              end loop;
            end if;
          end loop;
          put_line("The new matrix A : "); put(A);
          F := A; -- back up for the next round
          Upper_Triangulate(A);
          put_line("After upper triangulation : "); put(A);
          r := Rank(A);
          put("The rank of A : "); put(r,1); put_line(".");
        end if;
        exit when (r >= n);
      end loop;
    end if;
  end Span_of_System;

  procedure Main is

  -- DESCRIPTION :
  --   Prompts the user for a Laurent polynomial system
  --   and computes the span of its support.

    lp : Link_to_Laur_Sys;
    nq,nv : integer32 := 0;

  begin
    new_line;
    put_line("Reading a Laurent polynomial system ...");
    get(lp);
    new_line;
    nq := lp'last;
    nv := integer32(Number_of_Unknowns(lp(lp'first)));
    put("-> read "); put(nq,1); put(" polynomials in ");
    put(nv,1); put_line(" variables");
    Span_of_System(lp.all,nv);
  end Main;

begin
  Main;
end ts_supspan;
