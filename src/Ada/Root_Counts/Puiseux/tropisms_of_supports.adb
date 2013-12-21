with text_io,integer_io;                 use text_io,integer_io;
with Standard_Natural_Vectors;
with Standard_Natural_VecVecs;
with Standard_Integer_Vectors;
with Standard_Integer_Vectors_io;        use Standard_Integer_Vectors_io;
with Standard_Integer_VecVecs;           use Standard_Integer_VecVecs;
with Standard_Integer_VecVecs_io;        use Standard_Integer_VecVecs_io;
with Standard_Floating_Vectors;
with Standard_Floating_Vectors_io;       use Standard_Floating_Vectors_io;
with Standard_Floating_VecVecs;
with Parallel_Directions;
with Parallel_Edges;

package body Tropisms_of_Supports is

  function Edges ( n : natural; s : Array_of_Lists ) return Array_of_Faces is

    res : Array_of_Faces(s'range);

  begin
    for i in s'range loop
      res(i) := Create(1,n,s(i));
    end loop;
    return res;
  end Edges;

  function Pivot ( v : Standard_Integer_Vectors.Vector ) return integer is

  -- DESCRIPTION :
  --   Returns either v'last+1 if v is zero, or the first nonzero entry in v.

  begin
    for i in v'range loop
      if v(i) /= 0
       then return i;
      end if;
    end loop;
    return v'last+1;
  end Pivot;

  function Is_Parallel ( v1,v2 : Standard_Integer_Vectors.Vector )
                       return boolean is

  -- DESCRIPTION :
  --   Returns true if the two vectors v1 and v2 are parallel.

    d : integer;
    p1 : constant integer := Pivot(v1);
    p2 : constant integer := Pivot(v2);

  begin
    if p1 /= p2 then
      return false;
    else
      for i in p1+1..v1'last loop
        d := v1(p1)*v2(i) - v2(p1)*v1(i);
        if d /= 0
         then return false;
        end if;
      end loop;
      return true;
    end if;
  end Is_Parallel;

  function Is_Parallel ( n : natural; f1,f2 : Face ) return boolean is

  -- DESCRIPTION :
  --   Returns true if the two edges are parallel.

    v1,v2 : Standard_Integer_Vectors.Vector(1..n);

  begin
   -- put_line("Checking faces "); 
   -- Standard_Integer_VecVecs_io.put(Standard_Integer_VecVecs.VecVec(f1.all));
   -- put_line(" and ");
   -- Standard_Integer_VecVecs_io.put(Standard_Integer_VecVecs.VecVec(f2.all));
    for i in 1..n loop
      v1(i) := f1(f1'first+1)(i) - f1(f1'first)(i);
      v2(i) := f2(f2'first+1)(i) - f2(f2'first)(i);
    end loop;
   -- put("  v1 : "); put(v1);
   -- put("  v2 : "); put(v2); new_line;
    return Is_Parallel(v1,v2);
  end Is_Parallel;

  procedure Show_Parallel_Edges
              ( n : in natural; s : in Array_of_Lists;
                e : in Array_of_Faces ) is

    e1,e2 : Faces;
    f1,f2 : Face;
    cnt : natural := 0;

  begin
    for i in e'first..e'last-1 loop
      put("#edges of support "); put(i,1);
      put(" : "); put(Length_Of(e(i)),1); new_line;
      e1 := e(i);
      for ii in 1..Length_Of(e1) loop
        f1 := Head_Of(e1);
        for j in i+1..e'last loop
          e2 := e(j);
          for jj in 1..Length_Of(e2) loop
            f2 := Head_Of(e2);
            if Is_Parallel(n,f1,f2) then
              put("edge "); put(ii,1);
              put(" of support ");  put(i,1);
              put(" and edge "); put(jj,1);
              put(" of support ");  put(j,1);
              put_line(" are parallel");
              cnt := cnt + 1;
            end if;
            e2 := Tail_Of(e2);
          end loop;
        end loop;
        e1 := Tail_Of(e1);
      end loop;
    end loop;
    put("found "); put(cnt,1); put_line(" parallel pairs ...");
  end Show_Parallel_Edges;

  procedure Sort_Parallel_Edges
              ( n : in natural; s : in Array_of_Lists;
                e : in Array_of_Faces ) is

    d : Standard_Floating_VecVecs.Array_of_VecVecs(e'range)
      := Parallel_Edges.Edge_Directions(n,e);
    sd : Standard_Natural_VecVecs.VecVec(d'range) := Parallel_Edges.Sort(d);
    ind : natural;
    sum : natural := Parallel_Edges.Sum_of_Lengths(d);
    snb,dst : Standard_Natural_Vectors.Vector(1..sum);
    v,w : Standard_Floating_Vectors.Link_to_Vector;

  begin
    for i in d'range loop
      put("sorted directions for support "); put(i,1); put_line(" :");
      for j in sd(i)'range loop
        ind := sd(i)(j);
        put(d(i)(ind).all,2); new_line;
      end loop;
    end loop;
    Parallel_Edges.Merge_Sort(d,sd,snb,dst);
    put_line("Sorted directions :");
    for i in 1..sum loop
      put(d(snb(i))(dst(i)).all,2); new_line;
    end loop;
    put_line("Parallel directions :");
    for i in 2..sum loop
      v := d(snb(i-1))(dst(i-1));
      w := d(snb(i))(dst(i));
      if Parallel_Directions.Equal(v,w) then
        put("edge "); put(dst(i-1),1);
        put(" of support "); put(snb(i-1),1);
        put(" is parallel to edge "); put(dst(i),1);
        put(" of support "); put(snb(i),1); new_line;
      end if;
    end loop;
  end Sort_Parallel_Edges;

end Tropisms_of_Supports;
