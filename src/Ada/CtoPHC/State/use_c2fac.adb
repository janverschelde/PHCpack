with text_io;                           use text_io;
with Interfaces.C;
with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Floating_Numbers;         use Standard_Floating_Numbers;
with Double_Double_Numbers;             use Double_Double_Numbers;
with Quad_Double_Numbers;               use Quad_Double_Numbers;
with Standard_Complex_Numbers;
with DoblDobl_Complex_Numbers;
with QuadDobl_Complex_Numbers;
with Standard_Natural_Vectors;
with Standard_Natural_VecVecs;
with Standard_Floating_Vectors;
with Standard_Complex_Vectors;
with Standard_Complex_VecVecs;
with Standard_Complex_Poly_Systems;
with DoblDobl_Complex_Poly_Systems;
with QuadDobl_Complex_Poly_Systems;
with Standard_Complex_Solutions;
with DoblDobl_Complex_Solutions;
with QuadDobl_Complex_Solutions;
with Standard_System_and_Solutions_io;
with Sampling_Machine;
with Witness_Sets_io;                   use Witness_Sets_io;
with Assignments_in_Ada_and_C;          use Assignments_in_Ada_and_C;
with Standard_PolySys_Container;
with Standard_Solutions_Container;
with DoblDobl_PolySys_Container;
with DoblDobl_Solutions_Container;
with QuadDobl_PolySys_Container;
with QuadDobl_Solutions_Container;
with PHCpack_Operations;
with Sampling_Operations;
with DoblDobl_Sampling_Operations;
with QuadDobl_Sampling_Operations;
with Monodromy_Partitions;
with Monodromy_Permutations;

function use_c2fac ( job : integer32;
                     a : C_intarrs.Pointer;
		     b : C_intarrs.Pointer;
                     c : C_dblarrs.Pointer ) return integer32 is

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

  function Job1 return integer32 is -- prompts for a standard witness set

    use Standard_Complex_Poly_Systems;
    use Standard_Complex_Solutions;

    p : Link_to_Poly_Sys;
    sols : Solution_List;
    dim : natural32;
    data : Standard_Natural_Vectors.Vector(1..2);

  begin
    Standard_Read_Embedding(p,sols,dim);
    Standard_PolySys_Container.Initialize(p.all);
    Standard_Solutions_Container.Initialize(sols);
   -- Sampling_Operations.Initialize(p.all,sols,dim);
    data(1) := dim;
    data(2) := Length_Of(sols);
   -- put("The dimension is "); put(data(1),1); new_line;
   -- put("The degree is "); put(data(2),1); new_line;
    Assign(p'last,a);
    Assign(data,b);
    return 0;
  exception 
    when others =>
      put_line("Exception raised when reading witness set.");
      return 41;
  end Job1;

  function Job31 return integer32 is -- prompts for a dobldobl witness set

    use DoblDobl_Complex_Poly_Systems;
    use DoblDobl_Complex_Solutions;

    p : Link_to_Poly_Sys;
    sols : Solution_List;
    dim : natural32;
    data : Standard_Natural_Vectors.Vector(1..2);

  begin
    DoblDobl_Read_Embedding(p,sols,dim);
    DoblDobl_PolySys_Container.Initialize(p.all);
    DoblDobl_Solutions_Container.Initialize(sols);
   -- Sampling_Operations.Initialize(p.all,sols,dim);
    data(1) := dim;
    data(2) := Length_Of(sols);
   -- put("The dimension is "); put(data(1),1); new_line;
   -- put("The degree is "); put(data(2),1); new_line;
    Assign(p'last,a);
    Assign(data,b);
    return 0;
  exception 
    when others =>
      put_line("Exception raised when reading witness set.");
      return 631;
  end Job31;

  function Job61 return integer32 is -- prompts for a quaddobl witness set

    use QuadDobl_Complex_Poly_Systems;
    use QuadDobl_Complex_Solutions;

    p : Link_to_Poly_Sys;
    sols : Solution_List;
    dim : natural32;
    data : Standard_Natural_Vectors.Vector(1..2);

  begin
    QuadDobl_Read_Embedding(p,sols,dim);
    QuadDobl_PolySys_Container.Initialize(p.all);
    QuadDobl_Solutions_Container.Initialize(sols);
   -- Sampling_Operations.Initialize(p.all,sols,dim);
    data(1) := dim;
    data(2) := Length_Of(sols);
   -- put("The dimension is "); put(data(1),1); new_line;
   -- put("The degree is "); put(data(2),1); new_line;
    Assign(p'last,a);
    Assign(data,b);
    return 0;
  exception 
    when others =>
      put_line("Exception raised when reading witness set.");
      return 661;
  end Job61;

  function Job2 return integer32 is -- initializes standard sampling machine

    use Standard_Complex_Poly_Systems;
    use Standard_Complex_Solutions;

    lp : constant Link_to_Poly_Sys := Standard_PolySys_Container.Retrieve;
    sols : constant Solution_List := Standard_Solutions_Container.Retrieve;
    va : constant C_Integer_Array := C_intarrs.Value(a);
    dim : constant integer32 := integer32(va(va'first));

  begin
   -- put("initializing sampler with #solutions = ");
   -- put(Length_Of(sols),1); new_line;
    Sampling_Operations.Initialize(lp.all,sols,dim);
    return 0;
  exception
    when others =>
      put_line("Exception raised when initializing sampling machine.");
      return 42;
  end Job2;

  function Job32 return integer32 is -- initializes dobldobl sampling machine

    use DoblDobl_Complex_Poly_Systems;
    use DoblDobl_Complex_Solutions;

    lp : constant Link_to_Poly_Sys := DoblDobl_PolySys_Container.Retrieve;
    sols : constant Solution_List := DoblDobl_Solutions_Container.Retrieve;
    va : constant C_Integer_Array := C_intarrs.Value(a);
    dim : constant integer32 := integer32(va(va'first));

  begin
   -- put("initializing sampler with #solutions = ");
   -- put(Length_Of(sols),1); new_line;
    DoblDobl_Sampling_Operations.Initialize(lp.all,sols,dim);
    return 0;
  exception
    when others =>
      put_line("Exception raised when initializing sampling machine.");
      return 632;
  end Job32;

  function Job62 return integer32 is -- initializes quaddobl sampling machine

    use QuadDobl_Complex_Poly_Systems;
    use QuadDobl_Complex_Solutions;

    lp : constant Link_to_Poly_Sys := QuadDobl_PolySys_Container.Retrieve;
    sols : constant Solution_List := QuadDobl_Solutions_Container.Retrieve;
    va : constant C_Integer_Array := C_intarrs.Value(a);
    dim : constant integer32 := integer32(va(va'first));

  begin
   -- put("initializing sampler with #solutions = ");
   -- put(Length_Of(sols),1); new_line;
    QuadDobl_Sampling_Operations.Initialize(lp.all,sols,dim);
    return 0;
  exception
    when others =>
      put_line("Exception raised when initializing sampling machine.");
      return 662;
  end Job62;

  function Job3 return integer32 is -- assign standard coefficient

    use Standard_Complex_Numbers;

    va : constant C_Integer_Array := C_intarrs.Value(a);
    vb : constant C_Integer_Array := C_intarrs.Value(b);
    i : constant integer32 := integer32(va(va'first));
    j : constant integer32 := integer32(vb(vb'first));
    vc : constant C_Double_Array
       := C_DblArrs.Value(c,Interfaces.C.ptrdiff_t(2));
    re : constant double_float := double_float(vc(0));
    im : constant double_float := double_float(vc(1));
    cf : constant Complex_Number := Create(re,im);

  begin
    Sampling_Operations.Assign_Slice(cf,i,j);
    return 0;
  exception
    when others =>
      put_line("Exception raised when assigning standard coefficient.");
      return 43;
  end Job3;

  function Job33 return integer32 is -- assign dobldobl coefficient

    use DoblDobl_Complex_Numbers;

    va : constant C_Integer_Array := C_intarrs.Value(a);
    vb : constant C_Integer_Array := C_intarrs.Value(b);
    i : constant integer32 := integer32(va(va'first));
    j : constant integer32 := integer32(vb(vb'first));
    vc : constant C_Double_Array
       := C_DblArrs.Value(c,Interfaces.C.ptrdiff_t(4));
    re_hi : constant double_float := double_float(vc(0));
    re_lo : constant double_float := double_float(vc(1));
    im_hi : constant double_float := double_float(vc(2));
    im_lo : constant double_float := double_float(vc(3));
    re : constant double_double := Create(re_hi,re_lo);
    im : constant double_double := Create(im_hi,im_lo);
    cf : constant Complex_Number := Create(re,im);

  begin
    DoblDobl_Sampling_Operations.Assign_Slice(cf,i,j);
    return 0;
  exception
    when others =>
      put_line("Exception raised when assigning dobldobl coefficient.");
      return 633;
  end Job33;

  function Job63 return integer32 is -- assign quaddobl coefficient

    use QuadDobl_Complex_Numbers;

    va : constant C_Integer_Array := C_intarrs.Value(a);
    vb : constant C_Integer_Array := C_intarrs.Value(b);
    i : constant integer32 := integer32(va(va'first));
    j : constant integer32 := integer32(vb(vb'first));
    vc : constant C_Double_Array
       := C_DblArrs.Value(c,Interfaces.C.ptrdiff_t(8));
    re_hihi : constant double_float := double_float(vc(0));
    re_lohi : constant double_float := double_float(vc(1));
    re_hilo : constant double_float := double_float(vc(2));
    re_lolo : constant double_float := double_float(vc(3));
    im_hihi : constant double_float := double_float(vc(4));
    im_lohi : constant double_float := double_float(vc(5));
    im_hilo : constant double_float := double_float(vc(6));
    im_lolo : constant double_float := double_float(vc(7));
    re : constant quad_double := Create(re_hihi,re_lohi,re_hilo,re_lolo);
    im : constant quad_double := Create(im_hihi,im_lohi,im_hilo,im_lolo);
    cf : constant Complex_Number := Create(re,im);

  begin
    QuadDobl_Sampling_Operations.Assign_Slice(cf,i,j);
    return 0;
  exception
    when others =>
      put_line("Exception raised when assigning quaddobl coefficient.");
      return 633;
  end Job63;

  function Job4 return integer32 is -- storing gamma constant

    use Standard_Complex_Numbers;

    va : constant C_Integer_Array := C_intarrs.Value(a);
    vc : constant C_Double_Array
       := C_DblArrs.Value(c,Interfaces.C.ptrdiff_t(2));
    re : constant double_float := double_float(vc(0));
    im : constant double_float := double_float(vc(1));
    gamma : constant Complex_Number := Create(re,im);
    i : constant integer32 := integer32(va(va'first));

  begin
    Sampling_Operations.Store_Gamma(gamma,i);
    return 0;
  exception
    when others =>
      put_line("Exception raised when storing gamma constant.");
      return 44;
  end Job4;

  function Job7 return integer32 is -- copy from sampler to container

     use Standard_Complex_Poly_Systems;

     p : constant Poly_Sys := Sampling_Machine.Embedded_System;

  begin
    Standard_PolySys_Container.Initialize(p);
    return 0;
  exception
    when others =>
      put_line("Exception raised when copying from sampler to container.");
      return 47;
  end Job7;

  function Job8 return integer32 is -- copy first solutions to container

    use Standard_Complex_Solutions;

    s : constant Solution_List := Sampling_Operations.Retrieve_First_Solutions;

  begin
    Standard_Solutions_Container.Clear;
    Standard_Solutions_Container.Initialize(s);
    return 0;
  exception
    when others =>
      put_line("Exception raised when copying first solutions to container.");
      return 48;
  end Job8;

  function Job9 return integer32 is -- solutions from grid to container

    use Standard_Complex_Solutions;

    va : constant C_Integer_Array := C_intarrs.Value(a);
    i : constant integer32 := integer32(va(va'first));
    s : constant Solution_List := Monodromy_Permutations.Retrieve(i);
    cp_s : Solution_List;  -- will be a copy of s

  begin -- since this will be traffic on the same node
    Copy(s,cp_s); -- we better make a copy
    Standard_Solutions_Container.Clear;
    Standard_Solutions_Container.Initialize(cp_s);
    return 0;
  exception
    when others =>
      put_line("Exception when copying from grid to container.");
      return 59;
  end Job9;

  function Job10 return integer32 is

    va : constant C_Integer_Array := C_intarrs.Value(a);
    vb : constant C_Integer_Array
       := C_intarrs.Value(b,Interfaces.C.ptrdiff_t(2));
    n : constant integer32 := integer32(va(va'first));
    d : constant integer32 := integer32(vb(0));
    k : constant integer32 := integer32(vb(1));

  begin
    -- put("initializing monodromy_permutations with ");
    -- put("  n = "); put(n,1);
    -- put("  d = "); put(d,1);
    -- put("  k = "); put(k,1); new_line;
    Monodromy_Permutations.Initialize(n,d,k);
    return 0;
  exception
    when others =>
      put_line("Exception when initializing Monodromy_Permutations.");
      return 50;
  end Job10;

  function Job11 return integer32 is -- solutions to Monodromy_Permutations

    use Standard_Complex_Solutions;

    sols : constant Solution_List := Standard_Solutions_Container.Retrieve;

  begin
   -- put("storing "); put(Length_Of(sols),1);
   -- put_line(" into Monodromy_Permutations ...");
    Monodromy_Permutations.Store(sols);
    return 0;
  exception 
    when others =>
      put_line("Exception when storing solutions to Monodromy_Permutations.");
      return 51;
  end Job11;

  function Job12 return integer32 is -- compute monodromy permutation

    perm : constant Standard_Natural_Vectors.Vector
         := Monodromy_Permutations.Permutation;

  begin
   -- put("Ada permutation : ");
   -- for i in perm'range loop
   --   put(" "); put(perm(i),1);
   -- end loop;
   -- new_line;
    Assign(perm,b);
    return 0;
  exception
    when others =>
      put_line("Exception when computing a monodromy permutation.");
      return 52;
  end Job12;

  function Job13 return integer32 is -- update with permutation

    va : constant C_Integer_Array := C_intarrs.Value(a);
    n : constant integer32 := integer32(va(va'first));
    p : Standard_Natural_Vectors.Vector(1..n);
    nf : Standard_Natural_Vectors.Vector(1..2);

  begin
    Assign(natural32(n),b,p);
   -- put("Processing");
   -- for i in p'range loop
   --   put(" "); put(p(i),1);
   -- end loop;
    Monodromy_Permutations.Update_Decomposition(p,nf(1),nf(2));
   -- put(" : "); put(nf(1),1); put(" -> "); put(nf(2),1);
   -- new_line;
    Assign(nf,a);
    return 0;
  exception
    when others =>
      put_line("Exception raised when updating with permutation");
      return 53;
  end Job13;

  function Job14 return integer32 is -- writes the decomposition

    deco : constant Standard_Natural_VecVecs.Link_to_VecVec
         := Monodromy_Permutations.Decomposition;
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
      put_line("Exception raised when writing the decomposition.");
      return 54;
  end Job14;

  function Job15 return integer32 is -- applies linear trace test

    done : constant boolean
         := Monodromy_Permutations.Certify_with_Linear_Trace;

  begin
    if done
     then Assign(1,a);
     else Assign(0,a);
    end if;
    return 0;
  exception
    when others =>
      put_line("Exception raised when applying linear trace test.");
      return 55;
  end Job15;

  function Job16 return integer32 is -- return trace grid diagnostics

    use Standard_Complex_Numbers;

    err,dis : double_float;
    ada_c : Complex_Number;

  begin
    Monodromy_Permutations.Trace_Grid_Diagnostics(err,dis);
    ada_c := Create(err,dis);  -- a complex number is an array
    Assign(ada_c,c);
    return 0;
  exception
    when others =>
      put_line("Exception when returning diagnostics from trace grid.");
      return 56;
  end Job16;

  function Job17 return integer32 is -- compare trace sum differences

    va : constant C_Integer_Array := C_intarrs.Value(a);
    n : constant integer32 := integer32(va(va'first));
    f : Standard_Natural_Vectors.Vector(1..n);
    d : double_float;

  begin
    Assign(natural32(n),b,f);
    d := Monodromy_Permutations.Trace_Sum_Difference(f);
    Assign(d,c);
    return 0;
  exception
    when others =>
      put_line("Exception when comparing trace sum differences.");
      return 57;
  end Job17;

  function Job18 return integer32 is -- finding index of solution label

    va : constant C_Integer_Array
       := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(2));
    label : constant integer32 := integer32(va(0));
    slice : constant integer32 := integer32(va(1));
    result : integer32;

  begin
    result := Monodromy_Permutations.In_Slice(label,slice);
    Assign(result,b);
    return 0;
  exception
    when others =>
      put_line("Exception when finding index of a solution label.");
      return 58;
  end Job18;

  function Job19 return integer32 is

    va : constant C_Integer_Array := C_intarrs.Value(a);
    nb : constant integer32 := integer32(va(va'first));

  begin
    Sampling_Operations.Initialize_Slices(nb);
    return 0;
  exception
    when others =>
      put_line("Exception when initializing slices in Sampling_Operations");
      return 59;
  end Job19;

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
    Sampling_Operations.Add_Slices(v);
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
    v := Sampling_Operations.Retrieve_Slices(i);
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

  function Job22 return integer32 is -- set target slices

    va : constant C_Integer_Array := C_intarrs.Value(a);
    i : constant integer32 := integer32(va(va'first));

  begin
    Sampling_Operations.Set_Target_Slices(i);
    return 0;
  exception
    when others =>
      put_line("Exception raised when setting target slices.");
      return 62;
  end Job22;

  function Job23 return integer32 is -- completes one sampling loop

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
      sls := Monodromy_Permutations.Retrieve(start_label,start_slice);
    else
      sls := Monodromy_Permutations.Retrieve(start_label,start_slice+2);
                         -- +2 for trace grid
    end if;
   -- put_line("The retrieved start solution :");
   -- put(sls.all);
    tls := Sampling_Operations.Sample_Loop(start_slice,target_slice,sls);
    if target_slice = 0 then
      target_label := Monodromy_Permutations.Match(tls,target_slice,tol);
    else
      target_label := Monodromy_Permutations.Match(tls,target_slice+2,tol);
    end if;
    Assign(target_label,b);
    return 0;
  exception
    when others =>
      put_line("Exception raised when sampling one loop.");
      return 63;
  end Job23;

  function Job24 return integer32 is -- reads witness set from file

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
    p : Link_to_Poly_Sys;
    sols : Solution_List;
    dim : natural32;
    data : Standard_Natural_Vectors.Vector(1..2);

  begin
    Open(file,in_file,filename);
    Standard_Read_Embedding(file,p,sols,dim);
    Standard_PolySys_Container.Initialize(p.all);
    Standard_Solutions_Container.Initialize(sols);
   -- Sampling_Operations.Initialize(p.all,sols,dim);
    data(1) := dim;
    data(2) := Length_Of(sols);
   -- put("The dimension is "); put(data(1),1); new_line;
   -- put("The degree is "); put(data(2),1); new_line;
    Assign(p'last,a);
    Assign(data,b);
    Close(file);
    return 0;
  exception 
    when others =>
      put_line("Exception when reading witness set from file.");
      return 64;
  end Job24;

  function Job25 return integer32 is -- write witness set to file

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
      put_line("Exception when writing witness set to file.");
      return 65;
  end Job25;

  function Job26 return integer32 is -- returns current #factors

    f : constant natural32
      := Monodromy_Permutations.Number_of_Irreducible_Factors;

  begin
    Assign(integer32(f),a);
    return 0;
  exception
    when others =>
      put_line("Exception when retrieving number of irreducible factors.");
      return 68;
  end Job26;

  function Job27 return integer32 is -- returns irreducible factor

    v_a : constant C_Integer_Array
        := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(1));
    k : constant integer32 := integer32(v_a(v_a'first));
    f : constant Standard_Natural_Vectors.Link_to_Vector
      := Monodromy_Permutations.Component(k);

  begin
    Assign(f'last,a);
    Assign(f.all,b);
    return 0;
  exception
    when others =>
      put_line("Exception when retrieving an irreducible factor.");
      return 69;
  end Job27;

  function Job28 return integer32 is -- set state to silent
  begin
    Monodromy_Permutations.stay_silent := true;
    return 0;
  end Job28;

  function Handle_Jobs return integer32 is
  begin
    case job is
      when 0 => Write_Menu; return 0;
      when 1 => return Job1; -- read witness set in double precision
      when 2 => return Job2; -- initialize standard sampling machine
      when 3 => return Job3; -- assigning standard coefficient of slice
      when 4 => return Job4; -- storing gamma constant
      when 5 => Sampling_Operations.Sample; return 0;
      when 6 => Sampling_Operations.Swap_Slices; return 0;
      when 7 => return Job7; -- copy from sampler to container
      when 8 => return Job8; -- copy first solutions to container
      when 9 => return Job9; -- solutions from grid to container
      when 10 => return Job10; -- initializing Monodromy_Permutations
      when 11 => return Job11; -- solutions to Monodromy_Permutations
      when 12 => return Job12; -- compute monodromy permutation
      when 13 => return Job13; -- update with permutation
      when 14 => return Job14; -- writes the decomposition
      when 15 => return Job15; -- apply linear trace test
      when 16 => return Job16; -- return trace grid diagnostics
      when 17 => return Job17; -- comparing trace sum differences
      when 18 => return Job18; -- finding index of solution label
      when 19 => return Job19; -- initializing slices in Sampling_Operations
      when 20 => return Job20; -- adding new slice to Sampling_Operations
      when 21 => return Job21; -- returning coefficients of a slice
      when 22 => return Job22; -- setting target slices
      when 23 => return Job23; -- completes one sampling loop
      when 24 => return Job24; -- reads witness set from file
      when 25 => return Job25; -- writes witness set to file
      when 26 => return Job26; -- returns current number of factors
      when 27 => return Job27; -- returns labels of points in component
      when 28 => return Job28; -- state of monodromy permutations to silent
      when 31 => return Job31; -- read witness set in double double precision
      when 32 => return Job32; -- initialize dobldobl sampling machine
      when 33 => return Job33; -- assign dobldobl coefficient of slice
      when 35 => DoblDobl_Sampling_Operations.Sample; return 0;
      when 36 => DoblDobl_Sampling_Operations.Swap_Slices; return 0;
      when 61 => return Job61; -- read witness set in quad double precision
      when 62 => return Job62; -- initialize quaddobl sampling machine
      when 63 => return Job63; -- assign quaddobl coefficient of slice
      when 65 => QuadDobl_Sampling_Operations.Sample; return 0;
      when 66 => QuadDobl_Sampling_Operations.Swap_Slices; return 0;
      when others => put_line("  Sorry.  Invalid operation."); return 1;
    end case;
  end Handle_Jobs;

begin
  return Handle_Jobs;
end use_c2fac;
