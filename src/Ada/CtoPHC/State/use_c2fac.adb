with text_io;                           use text_io;
with Interfaces.C;
with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Floating_Numbers;         use Standard_Floating_Numbers;
with Double_Double_Numbers;             use Double_Double_Numbers;
with Quad_Double_Numbers;               use Quad_Double_Numbers;
with Standard_Complex_Numbers;
with DoblDobl_Complex_Numbers;
with QuadDobl_Complex_Numbers;
with Standard_Random_Numbers;
with DoblDobl_Random_Numbers;
with QuadDobl_Random_Numbers;
with Standard_Natural_Vectors;
with Standard_Natural_VecVecs;
with Standard_Floating_Vectors;
with Standard_Complex_Vectors;
with Standard_Complex_VecVecs;
with Standard_Complex_Poly_Systems;
with Standard_Complex_Laur_Systems;
with DoblDobl_Complex_Poly_Systems;
with DoblDobl_Complex_Laur_Systems;
with QuadDobl_Complex_Poly_Systems;
with QuadDobl_Complex_Laur_Systems;
with Standard_Complex_Solutions;
with DoblDobl_Complex_Solutions;
with QuadDobl_Complex_Solutions;
with Standard_System_and_Solutions_io;
with DoblDobl_System_and_Solutions_io;
with QuadDobl_System_and_Solutions_io;
with Sampling_Laurent_Machine;
with DoblDobl_Sampling_Laurent_Machine;
with QuadDobl_Sampling_Laurent_Machine;
with Monodromy_Partitions;
with Assignments_in_Ada_and_C;          use Assignments_in_Ada_and_C;
with Standard_PolySys_Container;
with Standard_LaurSys_Container;
with Standard_Solutions_Container;
with DoblDobl_PolySys_Container;
with DoblDobl_LaurSys_Container;
with DoblDobl_Solutions_Container;
with QuadDobl_PolySys_Container;
with QuadDobl_LaurSys_Container;
with QuadDobl_Solutions_Container;
with PHCpack_Operations;
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

  function Job10 return integer32 is

    va : constant C_Integer_Array := C_intarrs.Value(a);
    vb : constant C_Integer_Array
       := C_intarrs.Value(b,Interfaces.C.ptrdiff_t(2));
    n : constant integer32 := integer32(va(va'first));
    d : constant integer32 := integer32(vb(0));
    k : constant integer32 := integer32(vb(1));

  begin
    Standard_Monodromy_Permutations.Initialize(n,d,k);
    return 0;
  exception
    when others =>
      put_line("Exception at init of Monodromy_Permutations.");
      return 50;
  end Job10;

  function Job40 return integer32 is -- init dobldobl monodromy permutations

    va : constant C_Integer_Array := C_intarrs.Value(a);
    vb : constant C_Integer_Array
       := C_intarrs.Value(b,Interfaces.C.ptrdiff_t(2));
    n : constant integer32 := integer32(va(va'first));
    d : constant integer32 := integer32(vb(0));
    k : constant integer32 := integer32(vb(1));

  begin
    DoblDobl_Monodromy_Permutations.Initialize(n,d,k);
    return 0;
  exception
    when others =>
      put_line("Exception at init of DoblDobl_Monodromy_Permutations.");
      return 640;
  end Job40;

  function Job70 return integer32 is -- init quaddobl monodromy permutations

    va : constant C_Integer_Array := C_intarrs.Value(a);
    vb : constant C_Integer_Array
       := C_intarrs.Value(b,Interfaces.C.ptrdiff_t(2));
    n : constant integer32 := integer32(va(va'first));
    d : constant integer32 := integer32(vb(0));
    k : constant integer32 := integer32(vb(1));

  begin
    QuadDobl_Monodromy_Permutations.Initialize(n,d,k);
    return 0;
  exception
    when others =>
      put_line("Exception at init of QuadDobl_Monodromy_Permutations.");
      return 670;
  end Job70;

  function Job11 return integer32 is -- standard sols to Monodromy_Permutations

    use Standard_Complex_Solutions;

    sols : constant Solution_List := Standard_Solutions_Container.Retrieve;

  begin
    Standard_Monodromy_Permutations.Store(sols);
    return 0;
  exception 
    when others =>
      put_line("Exception at standard solutions to Monodromy_Permutations.");
      return 51;
  end Job11;

  function Job41 return integer32 is -- dobldobl sols to Monodromy_Permutations

    use DoblDobl_Complex_Solutions;

    sols : constant Solution_List := DoblDobl_Solutions_Container.Retrieve;

  begin
    DoblDobl_Monodromy_Permutations.Store(sols);
    return 0;
  exception 
    when others =>
      put_line("Exception at dobldobl solutions to Monodromy_Permutations.");
      return 641;
  end Job41;

  function Job71 return integer32 is -- quaddobl sols to Monodromy_Permutations

    use QuadDobl_Complex_Solutions;

    sols : constant Solution_List := QuadDobl_Solutions_Container.Retrieve;

  begin
    QuadDobl_Monodromy_Permutations.Store(sols);
    return 0;
  exception 
    when others =>
      put_line("Exception at quaddobl solutions to Monodromy_Permutations.");
      return 671;
  end Job71;

  function Job12 return integer32 is -- standard monodromy permutation

    perm : constant Standard_Natural_Vectors.Vector
         := Standard_Monodromy_Permutations.Permutation;

  begin
    Assign(perm,b);
    return 0;
  exception
    when others =>
      put_line("Exception at a standard double monodromy permutation.");
      return 52;
  end Job12;

  function Job42 return integer32 is -- dobldobl monodromy permutation

    perm : constant Standard_Natural_Vectors.Vector
         := DoblDobl_Monodromy_Permutations.Permutation;

  begin
    Assign(perm,b);
    return 0;
  exception
    when others =>
      put_line("Exception at a double double monodromy permutation.");
      return 642;
  end Job42;

  function Job72 return integer32 is -- quaddobl monodromy permutation

    perm : constant Standard_Natural_Vectors.Vector
         := QuadDobl_Monodromy_Permutations.Permutation;

  begin
    Assign(perm,b);
    return 0;
  exception
    when others =>
      put_line("Exception at a double double monodromy permutation.");
      return 672;
  end Job72;

  function Job13 return integer32 is -- update with standard permutation

    va : constant C_Integer_Array := C_intarrs.Value(a);
    n : constant integer32 := integer32(va(va'first));
    p : Standard_Natural_Vectors.Vector(1..n);
    nf : Standard_Natural_Vectors.Vector(1..2);

  begin
    Assign(natural32(n),b,p);
    Standard_Monodromy_Permutations.Update_Decomposition(p,nf(1),nf(2));
    Assign(nf,a);
    return 0;
  exception
    when others =>
      put_line("Exception at update with standard double permutation");
      return 53;
  end Job13;

  function Job43 return integer32 is -- update with dobldobl permutation

    va : constant C_Integer_Array := C_intarrs.Value(a);
    n : constant integer32 := integer32(va(va'first));
    p : Standard_Natural_Vectors.Vector(1..n);
    nf : Standard_Natural_Vectors.Vector(1..2);

  begin
    Assign(natural32(n),b,p);
    DoblDobl_Monodromy_Permutations.Update_Decomposition(p,nf(1),nf(2));
    Assign(nf,a);
    return 0;
  exception
    when others =>
      put_line("Exception at update with double double permutation");
      return 643;
  end Job43;

  function Job73 return integer32 is -- update with quaddobl permutation

    va : constant C_Integer_Array := C_intarrs.Value(a);
    n : constant integer32 := integer32(va(va'first));
    p : Standard_Natural_Vectors.Vector(1..n);
    nf : Standard_Natural_Vectors.Vector(1..2);

  begin
    Assign(natural32(n),b,p);
    QuadDobl_Monodromy_Permutations.Update_Decomposition(p,nf(1),nf(2));
    Assign(nf,a);
    return 0;
  exception
    when others =>
      put_line("Exception at update with quad double permutation");
      return 673;
  end Job73;

  function Job14 return integer32 is -- writes the standard decomposition

    deco : constant Standard_Natural_VecVecs.Link_to_VecVec
         := Standard_Monodromy_Permutations.Decomposition;
    use Standard_Natural_VecVecs;

  begin
    if deco /= null then
      if PHCpack_Operations.Is_File_Defined then
        Monodromy_Partitions.Write_Factors
          (PHCpack_Operations.output_file,deco.all);
      else
        Monodromy_Partitions.Write_Factors(standard_output,deco.all);
      end if;
    end if;
    return 0;
  exception
    when others =>
      put_line("Exception when writing the standard decomposition.");
      return 54;
  end Job14;

  function Job44 return integer32 is -- writes the dobldobl decomposition

    deco : constant Standard_Natural_VecVecs.Link_to_VecVec
         := DoblDobl_Monodromy_Permutations.Decomposition;
    use Standard_Natural_VecVecs;

  begin
    if deco /= null then
      if PHCpack_Operations.Is_File_Defined then
        Monodromy_Partitions.Write_Factors
          (PHCpack_Operations.output_file,deco.all);
      else
        Monodromy_Partitions.Write_Factors(standard_output,deco.all);
      end if;
    end if;
    return 0;
  exception
    when others =>
      put_line("Exception when writing the dobldobl decomposition.");
      return 644;
  end Job44;

  function Job74 return integer32 is -- writes the quaddobl decomposition

    deco : constant Standard_Natural_VecVecs.Link_to_VecVec
         := QuadDobl_Monodromy_Permutations.Decomposition;
    use Standard_Natural_VecVecs;

  begin
    if deco /= null then
      if PHCpack_Operations.Is_File_Defined then
        Monodromy_Partitions.Write_Factors
          (PHCpack_Operations.output_file,deco.all);
      else
        Monodromy_Partitions.Write_Factors(standard_output,deco.all);
      end if;
    end if;
    return 0;
  exception
    when others =>
      put_line("Exception when writing the quaddobl decomposition.");
      return 674;
  end Job74;

  function Monodromy_Standard_Trace_Test
             ( a : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

  -- DESCRIPTION :
  --   Applies the linear trace test to stop the monodromy loops
  --   in double precision.

  -- ON ENTRY :
  --   vrblvl  is the verbose level.

  -- ON RETURN :
  --   a       in a[0] is 1 if done, or 1 if not done.

    done : constant boolean
         := Standard_Monodromy_Permutations.Certify_with_Linear_Trace;

  begin
    if vrblvl > 0 then
      put("-> in monodromy_interface.");
      put_line("Monodromy_Standard_Trace_Test ...");
    end if;
    if done
     then Assign(1,a);
     else Assign(0,a);
    end if;
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in monodromy_interface.");
        put_line("Monodromy_Standard_Trace_Test.");
      end if;
      return 55;
  end Monodromy_Standard_Trace_Test;

  function Monodromy_DoblDobl_Trace_Test
             ( a : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

  -- DESCRIPTION :
  --   Applies the linear trace test to stop the monodromy loops
  --   in double double precision.

  -- ON ENTRY :
  --   vrblvl  is the verbose level.

  -- ON RETURN :
  --   a       in a[0] is 1 if done, or 1 if not done.

    done : constant boolean
         := DoblDobl_Monodromy_Permutations.Certify_with_Linear_Trace;

  begin
    if vrblvl > 0 then
      put("-> in monodromy_interface.");
      put_line("Monodromy_DoblDobl_Trace_Test ...");
    end if;
    if done
     then Assign(1,a);
     else Assign(0,a);
    end if;
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in monodromy_interface.");
        put_line("Monodromy_DoblDobl_Trace_Test.");
      end if;
      return 645;
  end Monodromy_DoblDobl_Trace_Test;

  function Monodromy_QuadDobl_Trace_Test
             ( a : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

  -- DESCRIPTION :
  --   Applies the linear trace test to stop the monodromy loops
  --   in quad double precision.

  -- ON ENTRY :
  --   vrblvl  is the verbose level.

  -- ON RETURN :
  --   a       in a[0] is 1 if done, or 1 if not done.

    done : constant boolean
         := QuadDobl_Monodromy_Permutations.Certify_with_Linear_Trace;

  begin
    if vrblvl > 0 then
      put("-> in monodromy_interface.");
      put_line("Monodromy_QuadDobl_Trace_Test ...");
    end if;
    if done
     then Assign(1,a);
     else Assign(0,a);
    end if;
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in monodromy_interface.");
        put_line("Monodromy_QuadDobl_Trace_Test.");
      end if;
      return 675;
  end Monodromy_QuadDobl_Trace_Test;

  function Job16 return integer32 is -- standard trace grid diagnostics

    use Standard_Complex_Numbers;

    err,dis : double_float;
    ada_c : Complex_Number;

  begin
    Standard_Monodromy_Permutations.Trace_Grid_Diagnostics(err,dis);
    ada_c := Create(err,dis);  -- a complex number is an array
    Assign(ada_c,c);
    return 0;
  exception
    when others =>
      put_line("Exception at diagnostics from standard trace grid.");
      return 56;
  end Job16;

  function Job46 return integer32 is -- dobldobl trace grid diagnostics

    use Standard_Complex_Numbers; -- return hi_parts of double doubles

    dd_err,dd_dis : double_double;
    st_err,st_dis : double_float;
    ada_c : Complex_Number;

  begin
    DoblDobl_Monodromy_Permutations.Trace_Grid_Diagnostics(dd_err,dd_dis);
    st_err := to_double(dd_err);
    st_dis := to_double(dd_dis);
    ada_c := Create(st_err,st_dis);  -- a complex number is an array
    Assign(ada_c,c);
    return 0;
  exception
    when others =>
      put_line("Exception at diagnostics from dobldobl trace grid.");
      return 646;
  end Job46;

  function Job76 return integer32 is -- quaddobl trace grid diagnostics

    use Standard_Complex_Numbers; -- return hi_parts of quad doubles

    qd_err,qd_dis : quad_double;
    st_err,st_dis : double_float;
    ada_c : Complex_Number;

  begin
    QuadDobl_Monodromy_Permutations.Trace_Grid_Diagnostics(qd_err,qd_dis);
    st_err := to_double(qd_err);
    st_dis := to_double(qd_dis);
    ada_c := Create(st_err,st_dis);  -- a complex number is an array
    Assign(ada_c,c);
    return 0;
  exception
    when others =>
      put_line("Exception at diagnostics from quaddobl trace grid.");
      return 676;
  end Job76;

  function Job17 return integer32 is -- standard trace sum differences

    va : constant C_Integer_Array := C_intarrs.Value(a);
    n : constant integer32 := integer32(va(va'first));
    f : Standard_Natural_Vectors.Vector(1..n);
    d : double_float;

  begin
    Assign(natural32(n),b,f);
    d := Standard_Monodromy_Permutations.Trace_Sum_Difference(f);
    Assign(d,c);
    return 0;
  exception
    when others =>
      put_line("Exception when comparing standard trace sum differences.");
      return 57;
  end Job17;

  function Job47 return integer32 is -- dobldobl trace sum differences

    va : constant C_Integer_Array := C_intarrs.Value(a);
    n : constant integer32 := integer32(va(va'first));
    f : Standard_Natural_Vectors.Vector(1..n);
    dd_d : double_double;
    st_d : double_float;

  begin
    Assign(natural32(n),b,f);
    dd_d := DoblDobl_Monodromy_Permutations.Trace_Sum_Difference(f);
    st_d := to_double(dd_d);
    Assign(st_d,c);
    return 0;
  exception
    when others =>
      put_line("Exception when comparing dobldobl trace sum differences.");
      return 647;
  end Job47;

  function Job77 return integer32 is -- quaddobl trace sum differences

    va : constant C_Integer_Array := C_intarrs.Value(a);
    n : constant integer32 := integer32(va(va'first));
    f : Standard_Natural_Vectors.Vector(1..n);
    qd_d : quad_double;
    st_d : double_float;

  begin
    Assign(natural32(n),b,f);
    qd_d := QuadDobl_Monodromy_Permutations.Trace_Sum_Difference(f);
    st_d := to_double(qd_d);
    Assign(st_d,c);
    return 0;
  exception
    when others =>
      put_line("Exception when comparing quaddobl trace sum differences.");
      return 677;
  end Job77;

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

  function Job25 return integer32 is -- write standard witness set to file

    use Standard_Complex_Poly_Systems;
    use Standard_Complex_Solutions;

    v_a : constant C_Integer_Array
        := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(1));
    n : constant natural32 := natural32(v_a(v_a'first));
    n1 : constant Interfaces.C.size_t := Interfaces.C.size_t(n-1);
    v_b : constant C_Integer_Array(0..n1)
        := C_intarrs.Value(b,Interfaces.C.ptrdiff_t(n));
    filename : constant String(1..integer(n))
             := C_Integer_Array_to_String(n,v_b);
    file : file_type;
    lp : constant Link_to_Poly_Sys := Standard_PolySys_Container.Retrieve;
    sols : constant Solution_List := Standard_Solutions_Container.Retrieve;

  begin
    Create(file,out_file,filename);
    Standard_System_and_Solutions_io.put(file,lp.all,sols);
    Close(file);
    return 0;
  exception 
    when others =>
      put_line("Exception when writing a standard witness set to file.");
      return 65;
  end Job25;

  function Job55 return integer32 is -- write dobldobl witness set to file

    use DoblDobl_Complex_Poly_Systems;
    use DoblDobl_Complex_Solutions;

    v_a : constant C_Integer_Array
        := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(1));
    n : constant natural32 := natural32(v_a(v_a'first));
    n1 : constant Interfaces.C.size_t := Interfaces.C.size_t(n-1);
    v_b : constant C_Integer_Array(0..n1)
        := C_intarrs.Value(b,Interfaces.C.ptrdiff_t(n));
    filename : constant String(1..integer(n))
             := C_Integer_Array_to_String(n,v_b);
    file : file_type;
    lp : constant Link_to_Poly_Sys := DoblDobl_PolySys_Container.Retrieve;
    sols : constant Solution_List := DoblDobl_Solutions_Container.Retrieve;

  begin
    Create(file,out_file,filename);
    DoblDobl_System_and_Solutions_io.put(file,lp.all,sols);
    Close(file);
    return 0;
  exception 
    when others =>
      put_line("Exception when writing a dobldobl witness set to file.");
      return 655;
  end Job55;

  function Job85 return integer32 is -- write quaddobl witness set to file

    use QuadDobl_Complex_Poly_Systems;
    use QuadDobl_Complex_Solutions;

    v_a : constant C_Integer_Array
        := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(1));
    n : constant natural32 := natural32(v_a(v_a'first));
    n1 : constant Interfaces.C.size_t := Interfaces.C.size_t(n-1);
    v_b : constant C_Integer_Array(0..n1)
        := C_intarrs.Value(b,Interfaces.C.ptrdiff_t(n));
    filename : constant String(1..integer(n))
             := C_Integer_Array_to_String(n,v_b);
    file : file_type;
    lp : constant Link_to_Poly_Sys := QuadDobl_PolySys_Container.Retrieve;
    sols : constant Solution_List := QuadDobl_Solutions_Container.Retrieve;

  begin
    Create(file,out_file,filename);
    QuadDobl_System_and_Solutions_io.put(file,lp.all,sols);
    Close(file);
    return 0;
  exception 
    when others =>
      put_line("Exception when writing a quaddobl witness set to file.");
      return 685;
  end Job85;

  function Job26 return integer32 is -- number of standard factors

    f : constant natural32
      := Standard_Monodromy_Permutations.Number_of_Irreducible_Factors;

  begin
    Assign(integer32(f),a);
    return 0;
  exception
    when others =>
      put_line("Exception at number of standard irreducible factors.");
      return 68;
  end Job26;

  function Job56 return integer32 is -- number of dobldobl factors

    f : constant natural32
      := DoblDobl_Monodromy_Permutations.Number_of_Irreducible_Factors;

  begin
    Assign(integer32(f),a);
    return 0;
  exception
    when others =>
      put_line("Exception at number of dobldobl irreducible factors.");
      return 656;
  end Job56;

  function Job86 return integer32 is -- number of quaddobl factors

    f : constant natural32
      := QuadDobl_Monodromy_Permutations.Number_of_Irreducible_Factors;

  begin
    Assign(integer32(f),a);
    return 0;
  exception
    when others =>
      put_line("Exception at number of quaddobl irreducible factors.");
      return 686;
  end Job86;

  function Job27 return integer32 is -- standard irreducible factor

    v_a : constant C_Integer_Array
        := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(1));
    k : constant integer32 := integer32(v_a(v_a'first));
    f : constant Standard_Natural_Vectors.Link_to_Vector
      := Standard_Monodromy_Permutations.Component(k);

  begin
    Assign(f'last,a);
    Assign(f.all,b);
    return 0;
  exception
    when others =>
      put_line("Exception when retrieving an irreducible factor.");
      return 69;
  end Job27;

  function Job57 return integer32 is -- dobldobl irreducible factor

    v_a : constant C_Integer_Array
        := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(1));
    k : constant integer32 := integer32(v_a(v_a'first));
    f : constant Standard_Natural_Vectors.Link_to_Vector
      := DoblDobl_Monodromy_Permutations.Component(k);

  begin
    Assign(f'last,a);
    Assign(f.all,b);
    return 0;
  exception
    when others =>
      put_line("Exception at retrieving a dobldobl irreducible factor.");
      return 657;
  end Job57;

  function Job87 return integer32 is -- quaddobl irreducible factor

    v_a : constant C_Integer_Array
        := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(1));
    k : constant integer32 := integer32(v_a(v_a'first));
    f : constant Standard_Natural_Vectors.Link_to_Vector
      := QuadDobl_Monodromy_Permutations.Component(k);

  begin
    Assign(f'last,a);
    Assign(f.all,b);
    return 0;
  exception
    when others =>
      put_line("Exception at retrieving a quaddobl irreducible factor.");
      return 687;
  end Job87;

  function Job28 return integer32 is -- set state of standard to silent
  begin
    Standard_Monodromy_Permutations.stay_silent := true;
    return 0;
  end Job28;

  function Job58 return integer32 is -- set state of dobldobl to silent
  begin
    DoblDobl_Monodromy_Permutations.stay_silent := true;
    return 0;
  end Job58;

  function Job88 return integer32 is -- set state of quaddobl to silent
  begin
    QuadDobl_Monodromy_Permutations.stay_silent := true;
    return 0;
  end Job88;

  function Job30 return integer32 is -- set state of standard to verbose
  begin
    Standard_Monodromy_Permutations.stay_silent := false;
    return 0;
  end Job30;

  function Job60 return integer32 is -- set state of dobldobl to verbose
  begin
    DoblDobl_Monodromy_Permutations.stay_silent := false;
    return 0;
  end Job60;

  function Job90 return integer32 is -- set state of quaddobl to verbose
  begin
    QuadDobl_Monodromy_Permutations.stay_silent := false;
    return 0;
  end Job90;

  function Job29 return integer32 is -- standard random complex number

    use Standard_Complex_Numbers;
    res : constant Complex_Number := Standard_Random_Numbers.Random1;

  begin
    Assign(res,c);
    return 0;
  exception
    when others =>
      put_line("Exception at generating standard random complex number.");
      return 280;
  end Job29;

  function Job59 return integer32 is -- dobldobl random complex number

    use DoblDobl_Complex_Numbers;
    res : constant Complex_Number := DoblDobl_Random_Numbers.Random1;

  begin
    Assign(res,c);
    return 0;
  exception
    when others =>
      put_line("Exception at generating dobldobl random complex number.");
      return 659;
  end Job59;

  function Job89 return integer32 is -- quaddobl random complex number

    use QuadDobl_Complex_Numbers;
    res : constant Complex_Number := QuadDobl_Random_Numbers.Random1;

  begin
    Assign(res,c);
    return 0;
  exception
    when others =>
      put_line("Exception at generating quaddobl random complex number.");
      return 689;
  end Job89;

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
      when  10 => return Job10; -- initializing Monodromy_Permutations
      when  11 => return Job11; -- standard solutions to Monodromy_Permutations
      when  12 => return Job12; -- compute standard monodromy permutation
      when  13 => return Job13; -- update with standard permutation
      when  14 => return Job14; -- writes the standard decomposition
      when  15 => return Monodromy_Standard_Trace_Test(a,vrblvl-1);
      when  16 => return Job16; -- return standard trace grid diagnostics
      when  17 => return Job17; -- comparing standard trace sum differences
      when  18 => return Job18; -- finding index of solution label
      when  19 => return Job19; -- init slices in Standard_Sampling_Operations
      when  20 => return Job20; -- adding new slice to Sampling_Operations
      when  21 => return Job21; -- returning coefficients of a slice
      when  22 => return Job22; -- setting standard target slices
      when  23 => return Job23; -- one sampling loop in double precision
      when  24 => return Witness_Standard_Polynomial_Read(a,b,vrblvl-1);
      when  25 => return Job25; -- writes standard witness set to file
      when  26 => return Job26; -- returns number of standard factors
      when  27 => return Job27; -- returns labels in standard component
      when  28 => return Job28; -- make standard monodromy permutations silent
      when  29 => return Job29; -- standard random complex number
      when  30 => return Job30; -- make standard monodromy permutations verbose 
      when  31 => return Witness_DoblDobl_Polynomial_Prompt(a,b,vrblvl-1);
      when  32 => return Monodromy_DoblDobl_Initialize_Sampler(a,vrblvl-1);
      when  33 => return Monodromy_DoblDobl_Set_Coefficient(a,b,c,vrblvl-1);
      when  34 => return Monodromy_DoblDobl_Store_Gamma(a,c,vrblvl-1);
      when  35 => return Monodromy_DoblDobl_Sample(vrblvl-1);
      when  36 => return Monodromy_DoblDobl_Swap_Slices(vrblvl-1);
      when  37 => return Monodromy_DoblDobl_Copy_System(vrblvl-1);
      when  38 => return Monodromy_DoblDobl_Copy_Solutions(vrblvl-1);
      when  39 => return Monodromy_DoblDobl_Grid_Solutions(a,vrblvl-1);
      when  40 => return Job40; -- initialize dobldobl monodromy permutations
      when  41 => return Job41; -- dobldobl solutions to Monodromy_Permutations
      when  42 => return Job42; -- compute dobldobl monodromy permutation
      when  43 => return Job43; -- update with dobldobl permutation
      when  44 => return Job44; -- writes the dobldobl decomposition
      when  45 => return Monodromy_DoblDobl_Trace_Test(a,vrblvl-1);
      when  46 => return Job46; -- return dobldobl trace grid diagnostics
      when  47 => return Job47; -- comparing dobldobl trace sum differences
      when  48 => return Job48; -- index of dobldobl solution label
      when  49 => return Job49; -- init slices in DoblDobl_Sampling_Operations
      when  52 => return Job52; -- setting dobldobl target slices
      when  53 => return Job53; -- one sampling loop in dobldobl precision
      when  54 => return Witness_DoblDobl_Polynomial_Read(a,b,vrblvl-1);
      when  55 => return Job55; -- writes dobldobl witness set to file
      when  56 => return Job56; -- returns number of dobldobl factors
      when  57 => return Job57; -- returns labels in dobldobl component
      when  58 => return Job58; -- make dobldobl monodromy permutations silent
      when  59 => return Job59; -- random dobldobl complex number
      when  60 => return Job60; -- make dobldobl monodromy permutations verbose 
      when  61 => return Witness_QuadDobl_Polynomial_Prompt(a,b,vrblvl-1);
      when  62 => return Monodromy_QuadDobl_Initialize_Sampler(a,vrblvl-1);
      when  63 => return Monodromy_QuadDobl_Set_Coefficient(a,b,c,vrblvl-1);
      when  64 => return Monodromy_QuadDobl_Store_Gamma(a,c,vrblvl-1);
      when  65 => return Monodromy_QuadDobl_Sample(vrblvl-1);
      when  66 => return Monodromy_QuadDobl_Swap_Slices(vrblvl-1);
      when  67 => return Monodromy_QuadDobl_Copy_System(vrblvl-1);
      when  68 => return Monodromy_QuadDobl_Copy_Solutions(vrblvl-1);
      when  69 => return Monodromy_QuadDobl_Grid_Solutions(a,vrblvl-1);
      when  70 => return Job70; -- initialize quaddobl monodromy permutations
      when  71 => return Job71; -- quaddobl solutions to Monodromy_Permutations
      when  72 => return Job72; -- compute quaddobl monodromy permutation
      when  73 => return Job73; -- update with quaddobl permutation
      when  74 => return Job74; -- writes the quaddobl decomposition
      when  75 => return Monodromy_QuadDobl_Trace_Test(a,vrblvl-1);
      when  76 => return Job76; -- return quaddobl trace grid diagnostics
      when  77 => return Job77; -- comparing quaddobl trace sum differences
      when  78 => return Job78; -- index of quaddobl solution label
      when  79 => return Job79; -- init slices in QuadDobl_Sampling_Operations
      when  82 => return Job82; -- setting quaddobl target slices
      when  83 => return Job83; -- one sampling loop in quaddobl precision
      when  84 => return Witness_QuadDobl_Polynomial_Read(a,b,vrblvl-1);
      when  85 => return Job85; -- writes quaddobl witness set to file
      when  86 => return Job86; -- returns number of quaddobl factors
      when  87 => return Job87; -- returns labels in quaddobl component
      when  88 => return Job88; -- make quaddobl monodromy permutations silent
      when  89 => return Job89; -- random quaddobl complex number
      when  90 => return Job90; -- make quaddobl monodromy permutations verbose 
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
