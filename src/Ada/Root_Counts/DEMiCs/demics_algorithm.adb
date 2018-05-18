with text_io;                            use text_io;
with Interfaces.C;
with String_Splitters;
with Communications_with_User;           use Communications_with_User;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;
with Standard_Integer_Vectors;
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
             ( sup : Arrays_of_Integer_Vector_Lists.Array_of_Lists )
             return C_Integer_Array is

    lst : constant Interfaces.C.size_T := Interfaces.C.size_T(sup'last-1);
    res : C_Integer_Array(0..lst);

    use Lists_of_Integer_Vectors;

  begin
    for i in res'range loop
      res(i) := Interfaces.C.int(Length_Of(sup(integer32(i)+1)));
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

  procedure Call_DEMiCs
              ( mix : in Standard_Integer_Vectors.Link_to_Vector;
                supports : in Arrays_of_Integer_Vector_Lists.Array_of_Lists;
                verbose : in boolean := true ) is

    dim : constant integer32 := supports'last;
    return_of_call : integer32;
    mixtoc : constant C_Integer_Array := Mixture_Type(mix.all);
    crdtoc : constant C_Integer_Array := Cardinalities(supports);
    nbrpts : constant integer32 := Number_of_Points(mix.all,supports);
    lenpts : constant integer32 := nbrpts*dim;
    ptstoc : constant C_Integer_Array
           := Coordinates(lenpts,mix.all,supports);
    lif : constant C_Double_Array := Random_Lifting(nbrpts);

  begin
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
