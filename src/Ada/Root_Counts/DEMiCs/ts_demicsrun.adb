with text_io;                            use text_io;
with Interfaces.C;
with String_Splitters;
with Communications_with_User;           use Communications_with_User;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Integer_Vectors;
with Standard_Integer_Vectors_io;        use Standard_Integer_Vectors_io;
with Standard_Floating_VecVecs;
with Standard_Floating_VecVecs_io;       use Standard_Floating_VecVecs_io;
with Lists_of_Integer_Vectors;
with Arrays_of_Integer_Vector_Lists;
with Arrays_of_Integer_Vector_Lists_io;  use Arrays_of_Integer_Vector_Lists_io;
with Arrays_of_Floating_Vector_Lists;
with Standard_Complex_Poly_Systems;      use Standard_Complex_Poly_Systems;
with Standard_Complex_Poly_Systems_io;   use Standard_Complex_Poly_Systems_io;
with Supports_of_Polynomial_Systems;
with Mixed_Volume_Computation;
with Floating_Mixed_Subdivisions;        use Floating_Mixed_Subdivisions;
with Floating_Mixed_Subdivisions_io;
with C_Integer_Arrays;                   use C_Integer_Arrays;
with Lists_of_Strings;
with DEMiCs_Command_Line;
with DEMiCs_Output_Convertors;
with DEMiCs_Output_Data;

procedure ts_demicsrun is

-- DESCRIPTION :
--   Development of a run with DEMiCs, invoked by an Ada main program.

  function Mixture_Type
             ( mix : Standard_Integer_Vectors.Vector )
             return C_Integer_Array is

  -- DESCRIPTION :
  --   Returns the type of mixtures as an array suitable to pass to C.

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

  -- DESCRIPTION :
  --   Returns the length of each support list in sup
  --   as an array suitable to pass to C.

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

  -- DESCRIPTION :
  --   Returns the total number of points in the supports sup.

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

  -- DESCRIPTION :
  --   Returns the length of each support list in sup
  --   as an array suitable to pass to C.
  --
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

  function run_demics
             ( v,d,r : integer32; mtp,crd,pts : C_Integer_Array )
             return integer32;
  pragma import(C, run_demics, "demicsrun");

  -- DESCRIPTION :
  --   Interface to the C++ function demicsrun.

  -- ON ENTRY :
  --   v       0 or 1 whether silent or verbose;
  --   d       dimension;
  --   r       number of distinct support sets;
  --   mtp     mixture type;
  --   crd     cardinalities of each support set;
  --   pts     coordinates of the points in the each support.

  procedure Extract_Supports 
               ( p : in Poly_Sys;
                 mix : out Standard_Integer_Vectors.Link_to_Vector;
                 supports : out Arrays_of_Integer_Vector_Lists.Array_of_Lists;
                 verbose : in boolean := true ) is

  -- DESCRIPTION :
  --   Extracts the supports and computes the type of mixtures.

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

  -- DESCRIPTION :
  --   Calls DEMiCs on the given supports.

    dim : constant integer32 := supports'last;
    return_of_call : integer32;
    mixtoc : constant C_Integer_Array := Mixture_Type(mix.all);
    crdtoc : constant C_Integer_Array := Cardinalities(supports);
    nbrpts : constant integer32 := Number_of_Points(mix.all,supports);
    lenpts : constant integer32 := nbrpts*dim;
    ptstoc : constant C_Integer_Array
           := Coordinates(lenpts,mix.all,supports);

  begin
    if verbose then
      return_of_call := run_demics(1,dim,mix'last,mixtoc,crdtoc,ptstoc);
    else
      return_of_call := run_demics(0,dim,mix'last,mixtoc,crdtoc,ptstoc);
    end if;
  end Call_DEMiCs;

  procedure Show_Output is

  -- DESCRIPITON :
  --   Shows the output stored in DEMiCs_Output_Data.

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

  -- DESCRIPTION :
  --   Scans s to the first colon and then returns the rest of s,
  --   all characters of s after the first colon.

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
                verbose : in boolean := true ) is

  -- DESCRIPTION :
  --   Converts the data computed by DEMiCs to lifted supports
  --   and a mixed cell configuration.

  -- ON ENTRY :
  --   mix      type of mixture;
  --   sup      supports of a polynomial system.

    use DEMiCs_Command_Line;
    use DEMiCs_Output_Convertors;

    lifting : constant Standard_Floating_VecVecs.Link_to_VecVec
            := DEMiCs_Output_Data.Lifting_Values;
    lifsup : Arrays_of_Floating_Vector_Lists.Array_of_lists(mix'range)
           := DEMiCs_Output_Convertors.Apply_Lifting(mix.all,sup,lifting.all);
    tmp : Lists_of_Strings.List := DEMiCs_Output_Data.Retrieve_Cell_Indices;
    ls : String_Splitters.Link_to_String;
    nbr : constant integer32 := Number_of_Points_in_Cell(mix.all);
    cells,cells_last : Lists_of_Integer_Vectors.List;
    mcc : Mixed_Subdivision;
    mv : natural32;

  begin
    put_line("The lifted supports :");
    Floating_Mixed_Subdivisions_io.put(lifsup);
    while not Lists_of_Strings.Is_Null(tmp) loop
      ls := Lists_of_Strings.Head_Of(tmp);
      declare
        strcell : constant string := Extract_Indices(ls.all);
        idxcell : constant Standard_Integer_Vectors.Vector
                := Extract_Cell_Indices(nbr,mix.all,strcell,verbose);
      begin
        put("The cell indices : "); put(idxcell); new_line;
        Lists_of_Integer_Vectors.Append(cells,cells_last,idxcell);
      end;
      tmp := Lists_of_Strings.Tail_Of(tmp);
    end loop;
    mcc := Make_Mixed_Cells(dim,mix.all,cells,lifsup,verbose);
    put_line("The mixed-cell configuration :");
    Floating_Mixed_Subdivisions_io.put(natural32(dim),mix.all,mcc,mv);
    put("The mixed volume : "); put(mv,1); new_line;
  end Process_Output;

  procedure Compute_Mixed_Volume ( p : in Poly_Sys ) is

  -- DESCRIPTION :
  --   Computes the mixed volume of p.

    dim : constant integer32 := p'last;
    ans : character;
    mix : Standard_Integer_Vectors.Link_to_Vector;
    sup : Arrays_of_Integer_Vector_Lists.Array_of_Lists(p'range);
    verbose : boolean;

  begin
    new_line;
    put("Verbose ? (y/n) ");
    Ask_Yes_or_No(ans);
    verbose := (ans = 'y');
    Extract_Supports(p,mix,sup,verbose);
    Call_DEMiCs(mix,sup,verbose);
    Show_Output;
    Process_Output(dim,mix,sup);
  end Compute_Mixed_Volume;

  procedure Main is

  -- DESCRIPTION :
  --   Prompts the user for a polynomial system
  --   and then prepares the input for DEMiCs.

    lp : Link_to_Poly_Sys;

  begin
    new_line;
    put_line("Reading a polynomial system ...");
    get(lp);
    Compute_Mixed_Volume(lp.all);
  end Main;

begin
  Main;
end ts_demicsrun;
