with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_VecVecs;           use Standard_Integer_VecVecs;
with Lists_of_Integer_Vectors;           use Lists_of_Integer_Vectors;
with Transforming_Integer32_Vector_Lists;
 use Transforming_Integer32_Vector_Lists;
with Arrays_of_Integer_Vector_Lists;     use Arrays_of_Integer_Vector_Lists;
with Standard_Dynamic32_Triangulations;  use Standard_Dynamic32_Triangulations;
with Unfolding_Subdivisions;             use Unfolding_Subdivisions;

package body Triangulations_and_Subdivisions is

-- REFINEMENT ROUTINES :

  procedure Refine ( n : in integer32; mic : in out Mixed_Cell ) is

  -- NOTE :
  --   Dynamic lifting will be applied with standard settings,
  --   under the assumption that there are only few points in the cell.

    support : constant List := Reduce(mic.pts(1),n+1);
    t : Triangulation;
    lifted,lifted_last : List;

  begin
    Dynamic_Lifting(support,false,true,0,lifted,lifted_last,t);
    mic.sub := new Mixed_Subdivision'(Deep_Create(n,t));
    Deep_Clear(lifted); Clear(t); 
      -- pity that Shallow_Clear(t) is not yet possible ...
  end Refine;

  procedure Refine ( n : in integer32; mixsub : in out Mixed_Subdivision ) is

  -- NOTE :
  --   Refines the mixed subdivision, under the safe assumption that 
  --   there is only one support set to deal with.

    res,res_last : Mixed_Subdivision;
    tmp : Mixed_Subdivision := mixsub;
    mic : Mixed_Cell;

  begin
    while not Is_Null(tmp) loop
      mic := Head_Of(tmp);
      if Length_Of(mic.pts(1)) > natural32(n)+1
       then Refine(n,mic);
      end if;
      Append(res,res_last,mic);
      tmp := Tail_Of(tmp);
    end loop;
    mixsub := res;
  end Refine;

-- TARGET PROCEDURES :

  function Deep_Create ( n : integer32; s : Simplex ) return Mixed_Cell is

    res : Mixed_Cell;
    ver : constant VecVec := Vertices(s);

  begin
    res.nor := new Standard_Integer_Vectors.Vector'(Normal(s));
    res.pts := new Array_of_Lists(1..1);
    res.pts(1) := Deep_Create(ver);
    return res;
  end Deep_Create;

  function Shallow_Create ( n : integer32; s : Simplex ) return Mixed_Cell is

    res : Mixed_Cell;
    ver : constant VecVec := Vertices(s);

  begin
    res.nor := new Standard_Integer_Vectors.Vector'(Normal(s));
    res.pts := new Array_of_Lists(1..1);
    res.pts(1) := Shallow_Create(ver);
    return res;
  end Shallow_Create;

  function Deep_Create ( n : integer32; t : Triangulation )
                       return Mixed_Subdivision is

    res,res_last : Mixed_Subdivision;
    tmp : Triangulation := t;

  begin
    while not Is_Null(tmp) loop
      Append(res,res_last,Deep_Create(n,Head_Of(tmp)));
      tmp := Tail_Of(tmp);
    end loop;
    return res;
  end Deep_Create;

  function Shallow_Create ( n : integer32; t : Triangulation )
                          return Mixed_Subdivision is

    res,res_last : Mixed_Subdivision;
    tmp : Triangulation := t;

  begin
    while not Is_Null(tmp) loop
      Append(res,res_last,Shallow_Create(n,Head_Of(tmp)));
      tmp := Tail_Of(tmp);
    end loop;
    return res;
  end Shallow_Create;

  function Deep_Create ( n : integer32; flatnor : Vector; t : Triangulation )
                       return Mixed_Subdivision is

    res,res_last : Mixed_Subdivision;
    tmp : Triangulation := t;
    s : Simplex;

  begin
    while not Is_Null(tmp) loop
      s := Head_Of(tmp);
      exit when (flatnor = Normal(s));
      Append(res,res_last,Deep_Create(n,s));
      tmp := Tail_Of(tmp);
    end loop;
    res := Merge(res);             -- merge cells with same inner normal
    Refine(n,res);                 -- refine the non-fine cells
    return res;
  end Deep_Create;

  function Shallow_Create
             ( n : integer32; flatnor : Vector; t : Triangulation )
             return Mixed_Subdivision is

    res,res_last : Mixed_Subdivision;
    tmp : Triangulation := t;
    s : Simplex;

  begin
    while not Is_Null(tmp) loop
      s := Head_Of(tmp);
      exit when (flatnor = Normal(s));
      Append(res,res_last,Shallow_Create(n,s));
      tmp := Tail_Of(tmp);
    end loop;
    res := Merge(res);             -- merge cells with same inner normal
    Refine(n,res);                 -- refine the non-fine cells
    return res;
  end Shallow_Create;

  function Non_Flat_Deep_Create ( n : integer32; t : Triangulation )
                                return Mixed_Subdivision is

    flatnor : Vector(1..n+1) := (1..n+1 => 0);

  begin
    flatnor(n+1) := 1;
    return Deep_Create(n,flatnor,t);
  end Non_Flat_Deep_Create;

  function Non_Flat_Shallow_Create ( n : integer32; t : Triangulation )
                                   return Mixed_Subdivision is

    flatnor : Vector(1..n+1) := (1..n+1 => 0);

  begin
    flatnor(n+1) := 1;
    return Shallow_Create(n,flatnor,t);
  end Non_Flat_Shallow_Create;

  function Deep_Create ( n : integer32; mixsub : Mixed_Subdivision )
                       return Triangulation is

    res : Triangulation;
    tmp : Mixed_Subdivision := mixsub;
    mic : Mixed_Cell;

  begin
    while not Is_Null(tmp) loop
      mic := Head_Of(tmp);
      declare
        v : VecVec(0..n);
        tmppts : List := mic.pts(mic.pts'first);
        s : Simplex;
      begin
        for i in v'range loop
          v(i) := new Standard_Integer_Vectors.Vector'(Head_Of(tmppts).all);
          tmppts := Tail_Of(tmppts);
          exit when Is_Null(tmppts);
        end loop;
        s := Create(v);
        Construct(s,res);
      end;
      tmp := Tail_Of(tmp);
    end loop;
    Connect(res);
    return res;
  end Deep_Create;

  function Shallow_Create ( n : integer32; mixsub : Mixed_Subdivision )
                          return Triangulation is

    res : Triangulation;
    tmp : Mixed_Subdivision := mixsub;
    mic : Mixed_Cell;

  begin
    while not Is_Null(tmp) loop
      mic := Head_Of(tmp);
      declare
        v : VecVec(0..n);
        tmppts : List := mic.pts(mic.pts'first);
        s : Simplex;
      begin
        for i in v'range loop
          v(i) := Head_Of(tmppts);
          tmppts := Tail_Of(tmppts);
          exit when Is_Null(tmppts);
        end loop;
        s := Create(v);
        Construct(s,res);
      end;
      tmp := Tail_Of(tmp);
    end loop;
    Connect(res);
    return res;
  end Shallow_Create;

end Triangulations_and_Subdivisions;
