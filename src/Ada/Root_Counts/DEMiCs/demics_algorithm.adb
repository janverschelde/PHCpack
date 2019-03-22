with text_io;                            use text_io;
with Interfaces.C;
with String_Splitters;
with Communications_with_User;           use Communications_with_User;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;
with Standard_Random_Numbers;
with Standard_Integer_Vectors_io;        use Standard_Integer_Vectors_io;
with Standard_Floating_VecVecs;
with Standard_Floating_VecVecs_io;       use Standard_Floating_VecVecs_io;
with Lists_of_Integer_Vectors;
with Arrays_of_Integer_Vector_Lists_io;  use Arrays_of_Integer_Vector_Lists_io;
with Supports_of_Polynomial_Systems;
with Floating_Lifting_Functions;
with Mixed_Volume_Computation;
with Floating_Mixed_Subdivisions_io;
with Lists_of_Strings;
with DEMiCs_Command_Line;
with DEMiCs_Output_Convertors;
with DEMiCs_Output_Data;

package body DEMiCs_Algorithm is

  procedure Extract_Supports 
               ( p : in Poly_Sys;
                 mix : out Standard_Integer_Vectors.Link_to_Vector;
                 supports : out Arrays_of_Integer_Vector_Lists.Array_of_Lists;
                 verbose : in boolean := true ) is

    perm : Standard_Integer_Vectors.Link_to_Vector;

  begin
    supports := Supports_of_Polynomial_Systems.Create(p);
    Mixed_Volume_Computation.Compute_Mixture(supports,mix,perm);
    if verbose then
      put_line("The supports : "); put(supports);
      put("The mixture type : "); put(mix.all); new_line;
    end if;
  end Extract_Supports;

  procedure Extract_Supports 
               ( p : in Poly_Sys;
                 mix,perm : out Standard_Integer_Vectors.Link_to_Vector;
                 supports : out Arrays_of_Integer_Vector_Lists.Array_of_Lists;
                 verbose : in boolean := true ) is
  begin
    supports := Supports_of_Polynomial_Systems.Create(p);
    Mixed_Volume_Computation.Compute_Mixture(supports,mix,perm);
    if verbose then
      put_line("The supports : "); put(supports);
      put("The mixture type : "); put(mix.all); new_line;
      put("The permutation : "); put(perm.all); new_line;
    end if;
  end Extract_Supports;

  procedure Extract_Supports 
              ( p : in Laur_Sys;
                mix : out Standard_Integer_Vectors.Link_to_Vector;
                supports : out Arrays_of_Integer_Vector_Lists.Array_of_Lists;
                verbose : in boolean := true ) is

    perm : Standard_Integer_Vectors.Link_to_Vector;

  begin
    supports := Supports_of_Polynomial_Systems.Create(p);
    Mixed_Volume_Computation.Compute_Mixture(supports,mix,perm);
    if verbose then
      put_line("The supports : "); put(supports);
      put("The mixture type : "); put(mix.all); new_line;
    end if;
  end Extract_Supports;

  procedure Extract_Supports 
              ( p : in Laur_Sys;
                mix,perm : out Standard_Integer_Vectors.Link_to_Vector;
                supports : out Arrays_of_Integer_Vector_Lists.Array_of_Lists;
                verbose : in boolean := true ) is
  begin
    supports := Supports_of_Polynomial_Systems.Create(p);
    Mixed_Volume_Computation.Compute_Mixture(supports,mix,perm);
    if verbose then
      put_line("The supports : "); put(supports);
      put("The mixture type : "); put(mix.all); new_line;
      put("The permutation : "); put(perm.all); new_line;
    end if;
  end Extract_Supports;

  procedure Add_Artificial_Origin
              ( dim : in integer32;
                sup : in out Lists_of_Integer_Vectors.List;
                added : out boolean ) is

    tmp : Lists_of_Integer_Vectors.List := sup;
    last : Lists_of_Integer_Vectors.List;
    lpt : Standard_Integer_Vectors.Link_to_Vector;
    found : boolean := false;
    origin : constant Standard_Integer_Vectors.Vector(1..dim)
           := (1..dim => 0);

  begin
    while not Lists_of_Integer_Vectors.Is_Null(tmp) loop
      lpt := Lists_of_Integer_Vectors.Head_Of(tmp);
      found := false;
      for k in lpt'range loop
        if lpt(k) /= 0
         then found := false; exit;
        end if;
      end loop;
      exit when found;
      last := tmp;
      tmp := Lists_of_Integer_Vectors.Tail_Of(tmp);
    end loop;
    if found then
      added := false;
    else
      Lists_of_Integer_Vectors.Append(sup,last,origin);
      added := true;
    end if;
  end Add_Artificial_Origin;

  procedure Add_Artificial_Origins
              ( dim : in integer32;
                sup : in out Arrays_of_Integer_Vector_Lists.Array_of_Lists;
                nbadd : out integer32;
                added : out Standard_Integer_Vectors.Vector ) is

    origin_added : boolean;

  begin
    nbadd := 0;
    for k in sup'range loop
      Add_Artificial_Origin(dim,sup(k),origin_added);
      if origin_added then
        nbadd := nbadd + 1;
        added(k) := 1;
      else
        added(k) := 0;
      end if;
    end loop;
  end Add_Artificial_Origins;

  function Mixture_Type
             ( mix : Standard_Integer_Vectors.Vector )
             return C_Integer_Array is

    lst : constant Interfaces.C.size_T := Interfaces.C.size_T(mix'last-1);
    res : C_Integer_Array(0..lst);

  begin
    for i in res'range loop
      res(i) := Interfaces.C.int(mix(integer32(i)+1));
    end loop;
    return res;
  end Mixture_Type;

  function Cardinalities
             ( mix : Standard_Integer_Vectors.Vector;
               sup : Arrays_of_Integer_Vector_Lists.Array_of_Lists )
             return C_Integer_Array is

    lst : constant Interfaces.C.size_T := Interfaces.C.size_T(mix'last-1);
    res : C_Integer_Array(0..lst);
    supidx : integer32 := sup'first;
    mixidx : integer32 := mix'first-1; -- bug fix 03/22/19, added -1

    use Lists_of_Integer_Vectors;

  begin
    for i in res'range loop
      res(i) := Interfaces.C.int(Length_Of(sup(supidx)));
      mixidx := mixidx + 1;
      supidx := supidx + mix(mixidx);
    end loop;
    return res;
  end Cardinalities;

  function Number_of_Points
             ( mix : Standard_Integer_Vectors.Vector;
               sup : Arrays_of_Integer_Vector_Lists.Array_of_Lists )
             return integer32 is

    res : integer32 := 0;
    idx : integer32 := sup'first;

  begin
    for i in mix'range loop
      res := res + integer32(Lists_of_Integer_Vectors.Length_Of(sup(idx)));
      idx := idx + mix(i);
    end loop;
    return res;
  end Number_of_Points;

  function Coordinates
             ( nbr : integer32;
               mix : Standard_Integer_Vectors.Vector;
               sup : Arrays_of_Integer_Vector_Lists.Array_of_Lists )
             return C_Integer_Array is

    use Interfaces.C;
    use Lists_of_Integer_Vectors;

    lst : constant Interfaces.C.size_T := Interfaces.C.size_T(nbr-1);
    res : C_Integer_Array(0..lst);
    idx : Interfaces.C.size_T := 0;
    tmp : List;
    lpt : Standard_Integer_Vectors.Link_to_Vector;
    idxsup : integer32 := sup'first;

  begin
    for i in mix'range loop
      tmp := sup(idxsup);
      while not Is_Null(tmp) loop
        lpt := Head_Of(tmp);
        for j in lpt'range loop
          res(idx) := Interfaces.C.int(lpt(j));
          idx := idx + 1;
        end loop;
        tmp := Tail_Of(tmp);
      end loop;
      idxsup := idxsup + mix(i);
    end loop;
    return res;
  end Coordinates;

  function Random_Lifting
             ( mix : Standard_Integer_Vectors.Link_to_Vector;
               sup : Arrays_of_Integer_Vector_Lists.Array_of_Lists )
             return Standard_Floating_VecVecs.Link_to_VecVec is

    res : Standard_Floating_VecVecs.Link_to_VecVec;
    resrep : Standard_Floating_VecVecs.VecVec(mix'range);
    idx : integer32 := 1;
    len : integer32;

  begin
    for i in resrep'range loop
      len := integer32(Lists_of_Integer_Vectors.Length_Of(sup(idx)));
      declare
        vals : Standard_Floating_Vectors.Vector(1..len);
      begin
        for j in 1..len loop
          vals(j) := Standard_Random_Numbers.Random;
        end loop;
        resrep(i) := new Standard_Floating_Vectors.Vector'(vals);
      end;
      idx := idx + mix(i);
    end loop;
    res := new Standard_Floating_VecVecs.VecVec'(resrep);
    return res;
  end Random_Lifting;

  function Random_Lifting
             ( mix : Standard_Integer_Vectors.Link_to_Vector;
               sup : Arrays_of_Integer_Vector_Lists.Array_of_Lists;
               stlb : double_float;
               added : Standard_Integer_Vectors.Vector )
             return Standard_Floating_VecVecs.Link_to_VecVec is

    res : Standard_Floating_VecVecs.Link_to_VecVec;
    resrep : Standard_Floating_VecVecs.VecVec(mix'range);
    idx : integer32 := 1;
    len : integer32;

  begin
    for i in resrep'range loop
      len := integer32(Lists_of_Integer_Vectors.Length_Of(sup(idx)));
      declare
        vals : Standard_Floating_Vectors.Vector(1..len);
      begin
        for j in 1..len loop
          vals(j) := Standard_Random_Numbers.Random;
        end loop;
        if added(idx) = 1          -- artificial origin is the last point
         then vals(len) := stlb;   -- and receives the lifting stlb
        end if;
        resrep(i) := new Standard_Floating_Vectors.Vector'(vals);
      end;
      idx := idx + mix(i);
    end loop;
    res := new Standard_Floating_VecVecs.VecVec'(resrep);
    return res;
  end Random_Lifting;

  function Random_Lifting ( nbr : integer32 ) return C_Double_Array is

    lst : constant Interfaces.C.size_T := Interfaces.C.size_T(nbr-1);
    res : C_Double_Array(0..lst);
    rnd : double_float;

  begin
    for i in res'range loop
      rnd := Floating_Lifting_Functions.Random_Lift(0.0,1.0);
      res(i) := Interfaces.C.double(rnd);
    end loop;
    return res;
  end Random_Lifting;

  function Size ( v : Standard_Floating_VecVecs.Link_to_VecVec )
                return integer32 is

    res : integer32 := 0;

  begin
    for i in v'range loop
      res := res + v(i)'last;
    end loop;
    return res;
  end Size;

  function Flatten ( v : Standard_Floating_VecVecs.Link_to_VecVec )
                   return Standard_Floating_Vectors.Vector is

    res : Standard_Floating_Vectors.Vector(1..Size(v));
    idx : integer32 := 0;

  begin
    for i in v'range loop
      for j in v(i)'range loop
        idx := idx + 1;
        res(idx) := v(i)(j);
      end loop;
    end loop;
    return res;
  end Flatten;

  function Copy_Lifting
             ( lif : Standard_Floating_Vectors.Vector )
             return C_Double_Array is

    nbr : constant integer32 := lif'last;
    lst : constant Interfaces.C.size_T := Interfaces.C.size_T(nbr-1);
    res : C_Double_Array(0..lst);

  begin
    for i in res'range loop
      res(i) := Interfaces.C.double(lif(integer32(i)+1));
    end loop;
    return res;
  end Copy_Lifting;

  procedure Call_DEMiCs
              ( mix : in Standard_Integer_Vectors.Link_to_Vector;
                supports : in Arrays_of_Integer_Vector_Lists.Array_of_Lists;
                nbrpts : in integer32; lif : in C_Double_Array;
                verbose : in boolean := true ) is

    dim : constant integer32 := supports'last;
    return_of_call : integer32;
    mixtoc : constant C_Integer_Array := Mixture_Type(mix.all);
    crdtoc : constant C_Integer_Array := Cardinalities(mix.all,supports);
    lenpts : constant integer32 := nbrpts*dim;
    ptstoc : constant C_Integer_Array
           := Coordinates(lenpts,mix.all,supports);
    crdsup : Standard_Integer_Vectors.Vector(mix'range);
    idx : Interfaces.C.size_T := lif'first;
    use Interfaces.C; -- needed to increment idx

  begin
    for i in crdtoc'range loop
      crdsup(integer32(i)+1) := integer32(crdtoc(i));
    end loop;
    if verbose then
      put("The dimension : "); put(dim,1); new_line;
      put("Cardinalities : "); put(crdsup); new_line;
      put("The total number of points : "); put(nbrpts,1); new_line;
    end if;
    DEMiCs_Output_Data.Initialize_Lifting(crdsup);
    for i in crdsup'range loop
      for j in 1..crdsup(i) loop
        DEMiCs_Output_Data.Assign_Lifting(i,j,double_float(lif(idx)));
        idx := idx + 1;
      end loop;
    end loop;
    if verbose then
     -- return_of_call := run_demics(1,dim,mix'last,mixtoc,crdtoc,ptstoc);
      put_line("The generated lifting values :");
      for i in lif'range loop
        put(" "); put(double_float(lif(i)));
      end loop;
      new_line;
      return_of_call := fly_demics(1,dim,mix'last,mixtoc,crdtoc,ptstoc,lif);
    else
     -- return_of_call := run_demics(0,dim,mix'last,mixtoc,crdtoc,ptstoc);
      return_of_call := fly_demics(0,dim,mix'last,mixtoc,crdtoc,ptstoc,lif);
    end if;
  end Call_DEMiCs;

  procedure Call_DEMiCs
              ( mix : in Standard_Integer_Vectors.Link_to_Vector;
                sup : in Arrays_of_Integer_Vector_Lists.Array_of_Lists;
                verbose : in boolean := true ) is

    nbrpts : constant integer32 := Number_of_Points(mix.all,sup);
    lif : C_Double_Array := Random_Lifting(nbrpts);

  begin
    Call_DEMiCs(mix,sup,nbrpts,lif,verbose);
  end Call_DEMiCs;

  procedure Call_DEMiCs
              ( mix : in Standard_Integer_Vectors.Link_to_Vector;
                sup : in out Arrays_of_Integer_Vector_Lists.Array_of_Lists;
                stable : in boolean; stlb : in double_float;
                verbose : in boolean := true ) is

    dim : constant integer32 := sup'last;
    nbrpts : constant integer32 := Number_of_Points(mix.all,sup);
    nbadded : integer32;
    added : Standard_Integer_Vectors.Vector(sup'range);
    lifting : Standard_Floating_VecVecs.Link_to_VecVec;

  begin
    if not stable then
      Call_DEMiCs(mix,sup,verbose);
    else
      Add_Artificial_Origins(dim,sup,nbadded,added);
      lifting := Random_Lifting(mix,sup,stlb,added);
      declare
        lifvals : constant Standard_Floating_Vectors.Vector
                := Flatten(lifting);
        lif : constant C_Double_Array := Copy_Lifting(lifvals);
      begin
        Call_DEMiCs(mix,sup,nbrpts,lif,verbose);
      end;
    end if;
  end Call_DEMiCs;

  procedure Call_DEMiCs
              ( mix : in Standard_Integer_Vectors.Link_to_Vector;
                supports : in Arrays_of_Integer_Vector_Lists.Array_of_Lists;
                lifting : in Standard_Floating_Vectors.Vector;
                verbose : in boolean := true ) is

    nbrpts : constant integer32 := Number_of_Points(mix.all,supports);
    lif : constant C_Double_Array := Copy_Lifting(lifting);

  begin
    if verbose then
      put("Total number of points : "); put(nbrpts,1); new_line;
      put("Number of lifting values : "); put(lifting'last,1);
      if nbrpts = lifting'last
       then put_line("  okay.");
       else put_line("  wrong!?");
      end if;
    end if;
    Call_DEMiCs(mix,supports,nbrpts,lif,verbose);
  end Call_DEMiCs;

  procedure Show_Output is

    lifting : constant Standard_Floating_VecVecs.Link_to_VecVec
            := DEMiCs_Output_Data.Lifting_Values;
    cells : constant Lists_of_Strings.List
          := DEMiCs_Output_Data.Retrieve_Cell_Indices;
    tmp : Lists_of_Strings.List := cells;
    ls : String_Splitters.Link_to_String;
    mv : constant integer32 := DEMiCs_Output_Data.mixed_volume;

  begin
    put_line("The lifting values :");
    put(lifting.all);
    put_line("The mixed cells :");
    while not Lists_of_Strings.Is_Null(tmp) loop
      ls := Lists_of_Strings.Head_Of(tmp);
      put_line(ls.all);
      tmp := Lists_of_Strings.Tail_Of(tmp);
    end loop;
    put("The mixed volume : "); put(mv,1); new_line;
  end Show_Output;

  function Extract_Indices ( s : string ) return string is

    idx : integer := s'first;

  begin
    while s(idx) /= ':' loop
      idx := idx + 1;
    end loop;
    return s(idx+1..s'last);
  end Extract_Indices;

  procedure Process_Output
              ( dim : in integer32;
                mix : in Standard_Integer_Vectors.Link_to_Vector;
                sup : in Arrays_of_Integer_Vector_Lists.Array_of_Lists;
                lif : out Arrays_of_Floating_Vector_Lists.Array_of_Lists;
                mcc : out Mixed_Subdivision;
                verbose : in boolean := true ) is

    use DEMiCs_Command_Line;
    use DEMiCs_Output_Convertors;

    lifting : constant Standard_Floating_VecVecs.Link_to_VecVec
            := DEMiCs_Output_Data.Lifting_Values;
    tmp : Lists_of_Strings.List := DEMiCs_Output_Data.Retrieve_Cell_Indices;
    ls : String_Splitters.Link_to_String;
    nbr : constant integer32 := Number_of_Points_in_Cell(mix.all);
    cells,cells_last : Lists_of_Integer_Vectors.List;

  begin
    lif := DEMiCs_Output_Convertors.Apply_Lifting(mix.all,sup,lifting.all);
    while not Lists_of_Strings.Is_Null(tmp) loop
      ls := Lists_of_Strings.Head_Of(tmp);
      declare
        strcell : constant string := Extract_Indices(ls.all);
        idxcell : constant Standard_Integer_Vectors.Vector
                := Extract_Cell_Indices(nbr,mix.all,strcell,verbose);
      begin
        if verbose then
          put("The cell indices : "); put(idxcell); new_line;
        end if;
        Lists_of_Integer_Vectors.Append(cells,cells_last,idxcell);
      end;
      tmp := Lists_of_Strings.Tail_Of(tmp);
    end loop;
    mcc := Make_Mixed_Cells(dim,mix.all,cells,lif,verbose);
  end Process_Output;

end DEMiCs_Algorithm;
