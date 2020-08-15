with text_io;                           use text_io;
with Interfaces.C;
with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Floating_Numbers;         use Standard_Floating_Numbers;
with Standard_Complex_Numbers;
with Standard_Floating_Vectors;
with Standard_Complex_Vectors;
with Standard_Complex_VecVecs;
with Standard_Complex_Laur_Systems;
with DoblDobl_Complex_Laur_Systems;
with QuadDobl_Complex_Laur_Systems;
with Standard_Complex_Solutions;
with DoblDobl_Complex_Solutions;
with QuadDobl_Complex_Solutions;
with Sampling_Laurent_Machine;
with DoblDobl_Sampling_Laurent_Machine;
with QuadDobl_Sampling_Laurent_Machine;
with Assignments_in_Ada_and_C;          use Assignments_in_Ada_and_C;
with Standard_LaurSys_Container;
with Standard_Solutions_Container;
with DoblDobl_LaurSys_Container;
with DoblDobl_Solutions_Container;
with QuadDobl_LaurSys_Container;
with QuadDobl_Solutions_Container;
with Standard_Sampling_Operations;
with DoblDobl_Sampling_Operations;
with QuadDobl_Sampling_Operations;
with Standard_Monodromy_Permutations;
with DoblDobl_Monodromy_Permutations;
with QuadDobl_Monodromy_Permutations;

with Witness_Interface;
with Monodromy_Interface;

function use_c2fac ( job : integer32;
                     a : C_intarrs.Pointer;
                     b : C_intarrs.Pointer;
                     c : C_dblarrs.Pointer;
                     vrblvl : integer32 := 0 ) return integer32 is

  procedure Write_Menu is
  begin
    new_line;
    put_line("General MENU to use factorization in PHCpack from C :");
    put_line("  0. display this menu;");
    put_line("  1. read an embedding of a polynomial system;");
    put_line("  2. initialize the sampling machine;");
    put_line("  3. assign coefficient of a slice;");
    put_line("  4. store a random gamma constant in PHCpack;");
    put_line("  5. compute a new witness set on the new slices;");
    put_line("  6. swaps slices and solution sets to turn back;");
    put_line("  7. copy embedded system from sampler to systems container;");
    put_line("  8. copy original solutions from sampler to container;");
    put_line("  9. copy solutions from monodromy grid in container;");
    put_line(" 10. initialize maximum number of monodromy loops;");
    put_line(" 11. copy solutions in container to Monodromy_Permutations;");
    put_line(" 12. compute permutation by last stored solution list;");
    put_line(" 13. update the monodromy breakup with a new permutation;");
    put_line(" 14. write the current monodromy breakup;");
    put_line(" 15. apply the linear trace to certify the monodromy breakup;");
    put_line(" 16. returns the diagnostics of the trace grid;");
    put_line(" 17. compute the difference in the trace sum for a factor;");
    put_line(" 18. find the index of a solution label in a slice;");
    put_line(" 19. initialize number of slices in Sampling_Operations;");
    put_line(" 20. adds a new slice to Sampling_Operations;");
    put_line(" 21. retrieves a slice from Sampling_Operations;");
    put_line(" 22. sets target slices to a previously stored set of slices;");
    put_line(" 23. completes one loop starting at one solution;");
    put_line(" 24. read a witness set from file, given by name;");
    put_line(" 25. writes system and solutions as a witness set to file;");
    put_line(" 26. returns current number of irreducible factors;");
    put_line(" 27. gets labels of witness points of an irreducible factor;");
    put_line(" 28. set the state of monodromy permutations to silent.");
  end Write_Menu;

  function Convert_to_Hyperplanes
             ( v : Standard_Floating_Vectors.Vector; k,n : integer32 )
             return Standard_Complex_VecVecs.VecVec is

  -- DESCRIPTION :
  --   Uses the numbers in v to create k hyperplanes in n-space.

    res : Standard_Complex_VecVecs.VecVec(1..k);
    ind : integer32 := v'first;

  begin
    for i in 1..k loop 
      declare
        hyp : Standard_Complex_Vectors.Vector(0..n);
      begin
        for j in 0..n loop
          hyp(j) := Standard_Complex_Numbers.Create(v(ind),v(ind+1));
          ind := ind+2;
        end loop;
        res(i) := new Standard_Complex_Vectors.Vector'(hyp);
      end;
    end loop;
    return res;
  end Convert_to_Hyperplanes;

  function Convert_to_Coefficients
             ( n : integer32; v : Standard_Complex_VecVecs.VecVec )
             return Standard_Floating_Vectors.Vector is

  -- DESCRIPTION :
  --   Returns a vector of range 1..n with the complex coefficients of v
  --   stored as sequences of real and imaginary parts.

    res : Standard_Floating_Vectors.Vector(1..n);
    lv : Standard_Complex_Vectors.Link_to_Vector;
    ind : integer32 := 0;

  begin
    for i in v'range loop
      lv := v(i);
      for j in lv'range loop
        ind := ind + 1; res(ind) := Standard_Complex_Numbers.REAL_PART(lv(j));
        ind := ind + 1; res(ind) := Standard_Complex_Numbers.IMAG_PART(lv(j));
      end loop;
    end loop;
    return res;
  end Convert_to_Coefficients;

  function Job18 return integer32 is -- index of standard solution label

    va : constant C_Integer_Array
       := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(2));
    label : constant integer32 := integer32(va(0));
    slice : constant integer32 := integer32(va(1));
    result : integer32;

  begin
    result := Standard_Monodromy_Permutations.In_Slice(label,slice);
    Assign(result,b);
    return 0;
  exception
    when others =>
      put_line("Exception when finding index of a solution label.");
      return 58;
  end Job18;

  function Job48 return integer32 is -- index of dobldobl solution label

    va : constant C_Integer_Array
       := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(2));
    label : constant integer32 := integer32(va(0));
    slice : constant integer32 := integer32(va(1));
    result : integer32;

  begin
    result := DoblDobl_Monodromy_Permutations.In_Slice(label,slice);
    Assign(result,b);
    return 0;
  exception
    when others =>
      put_line("Exception when finding index of a solution label.");
      return 648;
  end Job48;

  function Job78 return integer32 is -- index of quaddobl solution label

    va : constant C_Integer_Array
       := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(2));
    label : constant integer32 := integer32(va(0));
    slice : constant integer32 := integer32(va(1));
    result : integer32;

  begin
    result := QuadDobl_Monodromy_Permutations.In_Slice(label,slice);
    Assign(result,b);
    return 0;
  exception
    when others =>
      put_line("Exception when finding index of a solution label.");
      return 678;
  end Job78;

  function Job19 return integer32 is

    va : constant C_Integer_Array := C_intarrs.Value(a);
    nb : constant integer32 := integer32(va(va'first));

  begin
    Standard_Sampling_Operations.Initialize_Slices(nb);
    return 0;
  exception
    when others =>
      put_line("Exception at init slices in Standard_Sampling_Operations");
      return 59;
  end Job19;

  function Job49 return integer32 is

    va : constant C_Integer_Array := C_intarrs.Value(a);
    nb : constant integer32 := integer32(va(va'first));

  begin
    DoblDobl_Sampling_Operations.Initialize_Slices(nb);
    return 0;
  exception
    when others =>
      put_line("Exception at init slices in DoblDobl_Sampling_Operations");
      return 649;
  end Job49;

  function Job79 return integer32 is

    va : constant C_Integer_Array := C_intarrs.Value(a);
    nb : constant integer32 := integer32(va(va'first));

  begin
    QuadDobl_Sampling_Operations.Initialize_Slices(nb);
    return 0;
  exception
    when others =>
      put_line("Exception at init slices in QuadDobl_Sampling_Operations");
      return 679;
  end Job79;

  function Job20 return integer32 is -- add new slice to Sampling_Operations

    va : constant C_Integer_Array
       := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(3));
    nb_cff : constant integer32 := integer32(va(0));
    k : constant integer32 := integer32(va(1));
    n : constant integer32 := integer32(va(2));
    cff : Standard_Floating_Vectors.Vector(1..nb_cff);
    v : Standard_Complex_VecVecs.VecVec(1..k);

  begin
    Assign(natural32(nb_cff),c,cff);
    v := Convert_to_Hyperplanes(cff,k,n);
    Standard_Sampling_Operations.Add_Slices(v);
    return 0;
  exception
    when others =>
      put_line("Exception when adding new slice to Sampling_Operations");
      return 60;
  end Job20;

  function Job21 return integer32 is -- returning coefficients of slice

    vb : constant C_Integer_Array := C_intarrs.Value(b);
    i : constant integer32 := integer32(vb(vb'first));
    va : constant C_Integer_Array
       := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(3));
    nb_cff : constant integer32 := integer32(va(0));
   -- dim : constant natural32 := natural32(va(1));
   -- n : constant natural32 := natural32(va(2));
    v : Standard_Complex_VecVecs.Link_to_VecVec;
    cff : Standard_Floating_Vectors.Vector(1..nb_cff);
    use Standard_Complex_VecVecs;

  begin
    v := Standard_Sampling_Operations.Retrieve_Slices(i);
    if v /= null then
      cff := Convert_to_Coefficients(nb_cff,v.all);
      Assign(cff,c);
    end if;
    return 0;
  exception
    when others =>
      put_line("Exception when returning coefficients of a slice.");
      return 61;
  end Job21;

  function Job22 return integer32 is -- set target standard slices

    va : constant C_Integer_Array := C_intarrs.Value(a);
    i : constant integer32 := integer32(va(va'first));

  begin
    Standard_Sampling_Operations.Set_Target_Slices(i);
    return 0;
  exception
    when others =>
      put_line("Exception raised when setting target slices.");
      return 62;
  end Job22;

  function Job52 return integer32 is -- set target dobldobl slices

    va : constant C_Integer_Array := C_intarrs.Value(a);
    i : constant integer32 := integer32(va(va'first));

  begin
    DoblDobl_Sampling_Operations.Set_Target_Slices(i);
    return 0;
  exception
    when others =>
      put_line("Exception raised when setting dobldobl target slices.");
      return 652;
  end Job52;

  function Job82 return integer32 is -- set target quaddobl slices

    va : constant C_Integer_Array := C_intarrs.Value(a);
    i : constant integer32 := integer32(va(va'first));

  begin
    QuadDobl_Sampling_Operations.Set_Target_Slices(i);
    return 0;
  exception
    when others =>
      put_line("Exception raised when setting quaddobl target slices.");
      return 682;
  end Job82;

  function Job23 return integer32 is -- one standard sampling loop

    use Standard_Complex_Solutions;

    va : constant C_Integer_Array
       := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(2));
    vb : constant C_Integer_Array := C_intarrs.Value(b);
    start_slice : constant integer32 := integer32(va(0));
    target_slice : constant integer32 := integer32(va(1));
    start_label : constant integer32 := integer32(vb(vb'first));
    target_label : integer32;
    tol : constant double_float := 1.0E-8;
    sls,tls : Link_to_Solution;

  begin
    if start_slice = 0 then
      sls := Standard_Monodromy_Permutations.Retrieve
               (start_label,start_slice);
    else
      sls := Standard_Monodromy_Permutations.Retrieve
               (start_label,start_slice+2);
                         -- +2 for trace grid
    end if;
    tls := Standard_Sampling_Operations.Sample_Loop
             (start_slice,target_slice,sls);
    if target_slice = 0 then
      target_label := Standard_Monodromy_Permutations.Match
                        (tls,target_slice,tol);
    else
      target_label := Standard_Monodromy_Permutations.Match
                        (tls,target_slice+2,tol);
    end if;
    Assign(target_label,b);
    return 0;
  exception
    when others =>
      put_line("Exception raised when sampling one loop.");
      return 63;
  end Job23;

  function Job53 return integer32 is -- one dobldobl sampling loop

    use DoblDobl_Complex_Solutions;

    va : constant C_Integer_Array
       := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(2));
    vb : constant C_Integer_Array := C_intarrs.Value(b);
    start_slice : constant integer32 := integer32(va(0));
    target_slice : constant integer32 := integer32(va(1));
    start_label : constant integer32 := integer32(vb(vb'first));
    target_label : integer32;
    tol : constant double_float := 1.0E-8;
    sls,tls : Link_to_Solution;

  begin
    if start_slice = 0 then
      sls := DoblDobl_Monodromy_Permutations.Retrieve
               (start_label,start_slice);
    else
      sls := DoblDobl_Monodromy_Permutations.Retrieve
               (start_label,start_slice+2);
                         -- +2 for trace grid
    end if;
    tls := DoblDobl_Sampling_Operations.Sample_Loop
             (start_slice,target_slice,sls);
    if target_slice = 0 then
      target_label := DoblDobl_Monodromy_Permutations.Match
                        (tls,target_slice,tol);
    else
      target_label := DoblDobl_Monodromy_Permutations.Match
                        (tls,target_slice+2,tol);
    end if;
    Assign(target_label,b);
    return 0;
  exception
    when others =>
      put_line("Exception at sampling one loop in dobldobl precision.");
      return 653;
  end Job53;

  function Job83 return integer32 is -- one quaddobl sampling loop

    use QuadDobl_Complex_Solutions;

    va : constant C_Integer_Array
       := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(2));
    vb : constant C_Integer_Array := C_intarrs.Value(b);
    start_slice : constant integer32 := integer32(va(0));
    target_slice : constant integer32 := integer32(va(1));
    start_label : constant integer32 := integer32(vb(vb'first));
    target_label : integer32;
    tol : constant double_float := 1.0E-8;
    sls,tls : Link_to_Solution;

  begin
    if start_slice = 0 then
      sls := QuadDobl_Monodromy_Permutations.Retrieve
               (start_label,start_slice);
    else
      sls := QuadDobl_Monodromy_Permutations.Retrieve
               (start_label,start_slice+2);
                         -- +2 for trace grid
    end if;
    tls := QuadDobl_Sampling_Operations.Sample_Loop
             (start_slice,target_slice,sls);
    if target_slice = 0 then
      target_label := QuadDobl_Monodromy_Permutations.Match
                        (tls,target_slice,tol);
    else
      target_label := QuadDobl_Monodromy_Permutations.Match
                        (tls,target_slice+2,tol);
    end if;
    Assign(target_label,b);
    return 0;
  exception
    when others =>
      put_line("Exception at sampling one loop in quaddobl precision.");
      return 683;
  end Job83;

  function Job97 return integer32 is -- initializes standard Laurent sampler

    use Standard_Complex_Laur_Systems;
    use Standard_Complex_Solutions;

    lp : constant Link_to_Laur_Sys := Standard_LaurSys_Container.Retrieve;
    sols : constant Solution_List := Standard_Solutions_Container.Retrieve;
    va : constant C_Integer_Array := C_intarrs.Value(a);
    dim : constant integer32 := integer32(va(va'first));

  begin
    Standard_Sampling_Operations.Initialize(lp.all,sols,dim);
    return 0;
  exception
    when others =>
      put_line("Exception at initializing standard Laurent sampler.");
      return 804;
  end Job97;

  function Job98 return integer32 is -- initializes dobldobl Laurent sampler

    use DoblDobl_Complex_Laur_Systems;
    use DoblDobl_Complex_Solutions;

    lp : constant Link_to_Laur_Sys := DoblDobl_LaurSys_Container.Retrieve;
    sols : constant Solution_List := DoblDobl_Solutions_Container.Retrieve;
    va : constant C_Integer_Array := C_intarrs.Value(a);
    dim : constant integer32 := integer32(va(va'first));

  begin
    DoblDobl_Sampling_Operations.Initialize(lp.all,sols,dim);
    return 0;
  exception
    when others =>
      put_line("Exception at initializing dobldobl Laurent sampler.");
      return 805;
  end Job98;

  function Job99 return integer32 is -- initializes quaddobl Laurent sampler

    use QuadDobl_Complex_Laur_Systems;
    use QuadDobl_Complex_Solutions;

    lp : constant Link_to_Laur_Sys := QuadDobl_LaurSys_Container.Retrieve;
    sols : constant Solution_List := QuadDobl_Solutions_Container.Retrieve;
    va : constant C_Integer_Array := C_intarrs.Value(a);
    dim : constant integer32 := integer32(va(va'first));

  begin
    QuadDobl_Sampling_Operations.Initialize(lp.all,sols,dim);
    return 0;
  exception
    when others =>
      put_line("Exception at initializing quaddobl Laurent sampler.");
      return 806;
  end Job99;

  function Job100 return integer32 is -- standard sampler -> laur container

     use Standard_Complex_Laur_Systems;

     p : constant Laur_Sys := Sampling_Laurent_Machine.Embedded_System;

  begin
    Standard_LaurSys_Container.Initialize(p);
    return 0;
  exception
    when others =>
      put_line("Exception at copy standard sampler -> Laurent container.");
      return 807;
  end Job100;

  function Job101 return integer32 is -- dobldobl sampler -> laur container

     use DoblDobl_Complex_Laur_Systems;

     p : constant Laur_Sys
       := DoblDobl_Sampling_Laurent_Machine.Embedded_System;

  begin
    DoblDobl_LaurSys_Container.Initialize(p);
    return 0;
  exception
    when others =>
      put_line("Exception at copy dobldoblsampler -> Laurent container.");
      return 808;
  end Job101;

  function Job102 return integer32 is -- quaddobl sampler -> laur container

     use QuadDobl_Complex_Laur_Systems;

     p : constant Laur_Sys
       := QuadDobl_Sampling_Laurent_Machine.Embedded_System;

  begin
    QuadDobl_LaurSys_Container.Initialize(p);
    return 0;
  exception
    when others =>
      put_line("Exception at copy dobldoblsampler -> Laurent container.");
      return 808;
  end Job102;

  function Handle_Jobs return integer32 is

    use Witness_Interface;
    use Monodromy_Interface;

  begin
    case job is
      when   0 => Write_Menu; return 0;
      when   1 => return Witness_Standard_Polynomial_Prompt(a,b,vrblvl-1);
      when   2 => return Monodromy_Standard_Initialize_Sampler(a,vrblvl-1);
      when   3 => return Monodromy_Standard_Set_Coefficient(a,b,c,vrblvl-1);
      when   4 => return Monodromy_Standard_Store_Gamma(a,c,vrblvl-1);
      when   5 => return Monodromy_Standard_Sample(vrblvl-1);
      when   6 => return Monodromy_Standard_Swap_Slices(vrblvl-1);
      when   7 => return Monodromy_Standard_Copy_System(vrblvl-1);
      when   8 => return Monodromy_Standard_Copy_Solutions(vrblvl-1);
      when   9 => return Monodromy_Standard_Grid_Solutions(a,vrblvl-1);
      when  10 => return Monodromy_Standard_Init_Permutations(a,b,vrblvl-1);
      when  11 => return Monodromy_Standard_Perm_Solutions(vrblvl-1);
      when  12 => return Monodromy_Standard_Permutation(b,vrblvl-1);
      when  13 => return Monodromy_Standard_Update(a,b,vrblvl-1);
      when  14 => return Monodromy_Standard_Write(vrblvl-1);
      when  15 => return Monodromy_Standard_Trace_Test(a,vrblvl-1);
      when  16 => return Monodromy_Standard_Diagnostics(c,vrblvl-1);
      when  17 => return Monodromy_Standard_Trace_Sum(a,b,c,vrblvl-1);
      when  18 => return Job18; -- finding index of solution label
      when  19 => return Job19; -- init slices in Standard_Sampling_Operations
      when  20 => return Job20; -- adding new slice to Sampling_Operations
      when  21 => return Job21; -- returning coefficients of a slice
      when  22 => return Job22; -- setting standard target slices
      when  23 => return Job23; -- one sampling loop in double precision
      when  24 => return Witness_Standard_Polynomial_Read(a,b,vrblvl-1);
      when  25 => return Witness_Standard_Polynomial_Write(a,b,vrblvl-1);
      when  26 => return Monodromy_Standard_Factor_Count(a,vrblvl-1);
      when  27 => return Monodromy_Standard_Get_Factor(a,b,vrblvl-1);
      when  28 => return Monodromy_Standard_Set_Silent(vrblvl-1);
      when  29 => return Monodromy_Standard_Random(c,vrblvl-1);
      when  30 => return Monodromy_Standard_Set_Verbose(vrblvl-1);
      when  31 => return Witness_DoblDobl_Polynomial_Prompt(a,b,vrblvl-1);
      when  32 => return Monodromy_DoblDobl_Initialize_Sampler(a,vrblvl-1);
      when  33 => return Monodromy_DoblDobl_Set_Coefficient(a,b,c,vrblvl-1);
      when  34 => return Monodromy_DoblDobl_Store_Gamma(a,c,vrblvl-1);
      when  35 => return Monodromy_DoblDobl_Sample(vrblvl-1);
      when  36 => return Monodromy_DoblDobl_Swap_Slices(vrblvl-1);
      when  37 => return Monodromy_DoblDobl_Copy_System(vrblvl-1);
      when  38 => return Monodromy_DoblDobl_Copy_Solutions(vrblvl-1);
      when  39 => return Monodromy_DoblDobl_Grid_Solutions(a,vrblvl-1);
      when  40 => return Monodromy_DoblDobl_Init_Permutations(a,b,vrblvl-1);
      when  41 => return Monodromy_DoblDobl_Perm_Solutions(vrblvl-1);
      when  42 => return Monodromy_DoblDobl_Permutation(b,vrblvl-1);
      when  43 => return Monodromy_DoblDobl_Update(a,b,vrblvl-1);
      when  44 => return Monodromy_DoblDobl_Write(vrblvl-1);
      when  45 => return Monodromy_DoblDobl_Trace_Test(a,vrblvl-1);
      when  46 => return Monodromy_DoblDobl_Diagnostics(c,vrblvl-1);
      when  47 => return Monodromy_DoblDobl_Trace_Sum(a,b,c,vrblvl-1);
      when  48 => return Job48; -- index of dobldobl solution label
      when  49 => return Job49; -- init slices in DoblDobl_Sampling_Operations
      when  52 => return Job52; -- setting dobldobl target slices
      when  53 => return Job53; -- one sampling loop in dobldobl precision
      when  54 => return Witness_DoblDobl_Polynomial_Read(a,b,vrblvl-1);
      when  55 => return Witness_DoblDobl_Polynomial_Write(a,b,vrblvl-1);
      when  56 => return Monodromy_DoblDobl_Factor_Count(a,vrblvl-1);
      when  57 => return Monodromy_DoblDobl_Get_Factor(a,b,vrblvl-1);
      when  58 => return Monodromy_DoblDobl_Set_Silent(vrblvl-1);
      when  59 => return Monodromy_DoblDobl_Random(c,vrblvl-1);
      when  60 => return Monodromy_DoblDobl_Set_Verbose(vrblvl-1);
      when  61 => return Witness_QuadDobl_Polynomial_Prompt(a,b,vrblvl-1);
      when  62 => return Monodromy_QuadDobl_Initialize_Sampler(a,vrblvl-1);
      when  63 => return Monodromy_QuadDobl_Set_Coefficient(a,b,c,vrblvl-1);
      when  64 => return Monodromy_QuadDobl_Store_Gamma(a,c,vrblvl-1);
      when  65 => return Monodromy_QuadDobl_Sample(vrblvl-1);
      when  66 => return Monodromy_QuadDobl_Swap_Slices(vrblvl-1);
      when  67 => return Monodromy_QuadDobl_Copy_System(vrblvl-1);
      when  68 => return Monodromy_QuadDobl_Copy_Solutions(vrblvl-1);
      when  69 => return Monodromy_QuadDobl_Grid_Solutions(a,vrblvl-1);
      when  70 => return Monodromy_QuadDobl_Init_Permutations(a,b,vrblvl-1);
      when  71 => return Monodromy_QuadDobl_Perm_Solutions(vrblvl-1);
      when  72 => return Monodromy_QuadDobl_Permutation(b,vrblvl-1);
      when  73 => return Monodromy_QuadDobl_Update(a,b,vrblvl-1);
      when  74 => return Monodromy_QuadDobl_Write(vrblvl-1);
      when  75 => return Monodromy_QuadDobl_Trace_Test(a,vrblvl-1);
      when  76 => return Monodromy_QuadDobl_Diagnostics(c,vrblvl-1);
      when  77 => return Monodromy_QuadDobl_Trace_Sum(a,b,c,vrblvl-1);
      when  78 => return Job78; -- index of quaddobl solution label
      when  79 => return Job79; -- init slices in QuadDobl_Sampling_Operations
      when  82 => return Job82; -- setting quaddobl target slices
      when  83 => return Job83; -- one sampling loop in quaddobl precision
      when  84 => return Witness_QuadDobl_Polynomial_Read(a,b,vrblvl-1);
      when  85 => return Witness_QuadDobl_Polynomial_Write(a,b,vrblvl-1);
      when  86 => return Monodromy_QuadDobl_Factor_Count(a,vrblvl-1);
      when  87 => return Monodromy_QuadDobl_Get_Factor(a,b,vrblvl-1);
      when  88 => return Monodromy_QuadDobl_Set_Silent(vrblvl-1);
      when  89 => return Monodromy_QuadDobl_Random(c,vrblvl-1);
      when  90 => return Monodromy_QuadDobl_Set_Verbose(vrblvl-1);
      when  91 => return Witness_Standard_Laurent_Prompt(a,b,vrblvl-1);
      when  92 => return Witness_DoblDobl_Laurent_Prompt(a,b,vrblvl-1);
      when  93 => return Witness_QuadDobl_Laurent_Prompt(a,b,vrblvl-1);
      when  94 => return Witness_Standard_Laurent_Read(a,b,vrblvl-1);
      when  95 => return Witness_DoblDobl_Laurent_Read(a,b,vrblvl-1);
      when  96 => return Witness_QuadDobl_Laurent_Read(a,b,vrblvl-1);
      when  97 => return Job97; -- initialize standard Laurent sampler
      when  98 => return Job98; -- initialize dobldobl Laurent sampler
      when  99 => return Job99; -- initialize quaddobl Laurent sampler
      when 100 => return Job100; -- standard sampler -> Laurent container
      when 101 => return Job101; -- dobldobl sampler -> Laurent container
      when 102 => return Job102; -- quaddobl sampler -> Laurent container
      when others => put_line("  Sorry.  Invalid operation."); return -1;
    end case;
  end Handle_Jobs;

begin
  return Handle_Jobs;
end use_c2fac;
