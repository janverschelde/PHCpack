with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Integer_Vectors;
with Parallel_Directions;

package body Parallel_Edges is

  function Edge_Direction
              ( n : natural; f : Face )
              return Standard_Floating_Vectors.Vector is

    res : Standard_Floating_Vectors.Vector(1..n);
    dir : Standard_Integer_Vectors.Vector(1..n);
    piv : integer;

  begin
    for i in 1..n loop
      dir(i) := f(f'first+1)(i) - f(f'first)(i);
    end loop;
    piv := Parallel_Directions.Pivot(dir);
    res := Parallel_Directions.Normalize(dir,piv);
    return res;
  end Edge_Direction;

  function Edge_Directions 
              ( n : natural; f : Faces ) 
              return Standard_Floating_VecVecs.VecVec is

    res : Standard_Floating_VecVecs.VecVec(1..Length_Of(f));
    p : Faces := f;
    e : Face;
    v : Standard_Floating_Vectors.Vector(1..n);

  begin
    for i in res'range loop
      e := Head_Of(p);
      v := Edge_Direction(n,e);
      res(i) := new Standard_Floating_Vectors.Vector'(v);
      p := Tail_Of(p);
    end loop;
    return res;
  end Edge_Directions;

  function Edge_Directions
              ( n : natural; f : Array_of_Faces )
              return Standard_Floating_VecVecs.Array_of_VecVecs is

    res : Standard_Floating_VecVecs.Array_of_VecVecs(f'range);

  begin
    for i in res'range loop
      declare
        L : constant natural := Length_Of(f(i));
        v : Standard_Floating_VecVecs.VecVec(1..L) := Edge_Directions(n,f(i));
      begin
        res(i) := new Standard_Floating_VecVecs.VecVec'(v);
      end;
    end loop;
    return res;
  end Edge_Directions;

  function Sort ( v : Standard_Floating_VecVecs.VecVec )
                return Standard_Natural_Vectors.Vector is

    res : Standard_Natural_Vectors.Vector(v'range);
    min : Standard_Floating_Vectors.Link_to_Vector;

    use Parallel_Directions;

  begin
    for i in v'range loop
      res(i) := i;
    end loop;
    for i in v'first..v'last-1 loop
      min := v(res(i));
      for j in i+1..v'last loop
        if min > v(res(j)) then
          Parallel_Directions.Swap(res,i,j);
          min := v(res(j));
        end if;
      end loop;
    end loop;
    return res;
  end Sort;

  function Sort ( v : Standard_Floating_VecVecs.Array_of_VecVecs )
                return Standard_Natural_VecVecs.VecVec is

    res : Standard_Natural_VecVecs.VecVec(v'range);

  begin
    for i in v'range loop
      declare
        ind : Standard_Natural_Vectors.Vector(v(i)'range) := Sort(v(i).all);
      begin
        res(i) := new Standard_Natural_Vectors.Vector'(ind);
      end;
    end loop;
    return res;
  end Sort;

  procedure Merge_Sort
              ( v,w : in Standard_Floating_VecVecs.VecVec;
                sv,sw : in Standard_Natural_Vectors.Vector;
                vnb,wnb : in natural;
                snb,dst : out Standard_Natural_Vectors.Vector ) is

    use Parallel_Directions;

    vk : integer := v'first;
    wk : integer := w'first;
    ind : integer := snb'first;

  begin
    while vk <= v'last and wk <= w'last loop
      if v(sv(vk)) > w(sw(wk)) then
        snb(ind) := wnb;
        dst(ind) := sw(wk);
        wk := wk + 1;
      else
        snb(ind) := vnb;
        dst(ind) := sv(vk);
        vk := vk + 1;
      end if;
      ind := ind + 1;
    end loop;
    if vk > v'last then
      while wk <= w'last loop
        snb(ind) := wnb;
        dst(ind) := sw(wk);
        wk := wk + 1;
        ind := ind + 1;
      end loop;
    elsif wk > w'last then
      while wk <= w'last loop
        snb(ind) := vnb;
        dst(ind) := sv(vk);
        vk := vk + 1;
        ind := ind + 1;
      end loop;
    end if;
  end Merge_Sort;

  function Sum_of_Lengths
              ( v : in Standard_Floating_VecVecs.Array_of_VecVecs ) 
              return natural is

    res : natural := 0;

  begin
    for i in v'range loop
      res := res + v(i)'last;
    end loop;
    return res;
  end Sum_of_Lengths;

  procedure Merge_Sort
              ( v : in Standard_Floating_VecVecs.Array_of_VecVecs;
                sv : in Standard_Natural_VecVecs.VecVec;
                snb,dst : out Standard_Natural_Vectors.Vector ) is

    p : integer := -1;
    ind : Standard_Natural_Vectors.Vector(v'range);
    pos : integer := snb'first;
    k : integer;

  begin
    for i in ind'range loop
      ind(i) := v(i)'first;
    end loop;
    while p < 0 loop
      k := Parallel_Directions.Min(v,sv,ind);
      snb(pos) := k;
      dst(pos) := sv(k)(ind(k));
      pos := pos + 1;
      ind(k) := ind(k) + 1;
      p := Parallel_Directions.Last_Index(v,ind);
    end loop;
    if p > 0 then
      while pos <= snb'last loop
        snb(pos) := p;
        dst(pos) := sv(p)(ind(p));
        pos := pos + 1;
        ind(p) := ind(p) + 1;
        exit when ind(p) > v(p)'last;
      end loop;
    end if;
  end Merge_Sort;

end Parallel_Edges;
