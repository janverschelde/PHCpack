with text_io;                           use text_io;
with Interfaces.C;                      use Interfaces.C;
with Communications_with_User;
with Timing_Package;                    use Timing_Package;
with Characters_and_Numbers;            use Characters_and_Numbers;
with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;       use Standard_Natural_Numbers_io;
with Multprec_Natural_Numbers;          use Multprec_Natural_Numbers;
with Multprec_Natural_Numbers_io;       use Multprec_Natural_Numbers_io;
with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;         use Standard_Floating_Numbers;
with Double_Double_Numbers;             use Double_Double_Numbers;
with Quad_Double_Numbers;               use Quad_Double_Numbers;
with Standard_Complex_Numbers;
with DoblDobl_Complex_Numbers;
with QuadDobl_Complex_Numbers;
with Standard_Natural_Vectors;
with Standard_Natural_Vectors_io;       use Standard_Natural_Vectors_io;
with Standard_Natural_VecVecs;
with Standard_Complex_VecMats;
with DoblDobl_Complex_VecMats;
with QuadDobl_Complex_VecMats;
with Standard_Complex_Poly_Systems;
with DoblDobl_Complex_Poly_Systems;
with QuadDobl_Complex_Poly_Systems;
with Standard_Complex_Solutions;
with DoblDobl_Complex_Solutions;
with QuadDobl_Complex_Solutions;
with Brackets;                          use Brackets;
with Brackets_io;                       use Brackets_io;
with Checker_Moves;
with Intersection_Posets;               use Intersection_Posets;
with Standard_Solution_Posets;
with DoblDobl_Solution_Posets;
with QuadDobl_Solution_Posets;
with Resolve_Schubert_Problems;         use Resolve_Schubert_Problems;
with Drivers_for_Schubert_Induction;    use Drivers_for_Schubert_Induction;
with Standard_PolySys_Container;
with Standard_Solutions_Container;
with DoblDobl_PolySys_Container;
with DoblDobl_Solutions_Container;
with QuadDobl_PolySys_Container;
with QuadDobl_Solutions_Container;
with Assignments_in_Ada_and_C;          use Assignments_in_Ada_and_C;

function use_c2lrhom ( job : integer32;
                       a : C_intarrs.Pointer;
		       b : C_intarrs.Pointer;
                       c : C_dblarrs.Pointer ) return integer32 is

  procedure Get_Dimensions
              ( a : C_intarrs.Pointer; n,k,c : out integer32;
                verbose : out boolean ) is

  -- DESCRIPTION :
  --   Extracts the dimensions (n,k,c) from the array a,
  --   where n is the ambient dimension, k the dimension of the solution
  --   planes, c the number of intersection conditions, and verbose flags
  --   whether addition output should be written during the resolution.

    v : constant C_Integer_Array(0..3)
      := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(4));

  begin
    n := integer32(v(0));
    k := integer32(v(1));
    c := integer32(v(2));
    if integer32(v(3)) = 1
     then verbose := true;
     else verbose := false;
    end if;
  end Get_Dimensions;

  procedure Get_Dimensions2
              ( a : C_intarrs.Pointer; n,k,c,nchar : out integer32;
                verbose,verify : out boolean ) is

  -- DESCRIPTION :
  --   Extracts the dimensions (n,k,c) from the array a,
  --   where n is the ambient dimension, k the dimension of the solution
  --   planes, c the number of intersection conditions, and verbose flag
  --   whether addition output should be written during the resolution,
  --   followed by the verify flag, whether diagnostic verification is needed.
  --   The last value in a is the number of characters in the string
  --   for the output file.

    v : constant C_Integer_Array(0..5)
      := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(6));

  begin
    n := integer32(v(0));
    k := integer32(v(1));
    c := integer32(v(2));
    if integer32(v(3)) = 1
     then verbose := true;
     else verbose := false;
    end if;
    if integer32(v(4)) = 1
     then verify := true;
     else verify := false;
    end if;
    nchar := integer32(v(5));
  end Get_Dimensions2;

  function Get_Conditions
             ( b : C_intarrs.Pointer; k,nb : integer32 )
             return Array_of_Brackets is

  -- DESCRIPTOIN :
  --   Extracts the brackets from the array b,
 
  -- REQUIRED :
  --   The array b must contain k times nb integers.

  -- ON ENTRY :
  --   b      an array of integers;
  --   k      dimension of the solution planes;
  --   nb     the number of conditions. 

    res : Array_of_Brackets(1..nb);
    dim : constant integer32 := k*nb;
    dmc : constant Interfaces.C.size_t := Interfaces.C.size_t(dim);
    dm1 : constant Interfaces.C.size_t := Interfaces.C.size_t(dim-1);
    v : constant C_Integer_Array(0..dm1)
      := C_intarrs.Value(b,Interfaces.C.ptrdiff_t(dmc));
    ind : Interfaces.C.size_t := 0;

  begin
    for i in res'range loop
      declare
        brk : Bracket(1..k);
      begin
        for j in 1..k loop
          brk(j) := natural32(v(ind));
          ind := ind + 1;
        end loop;
        res(i) := new Bracket'(brk);
      end;
    end loop;
    return res;
  end Get_Conditions;

  function Get_File_Name 
             ( c : C_dblarrs.Pointer; nbchar : integer32 ) return string is

  -- DESCRIPTION :
  --   Extracts the file name (as many characters as the value of nbchar)
  --   from the input parameter c.

    res : string(1..integer(nbchar)) := (1..integer(nbchar) => ' ');
    val : constant C_Double_Array(0..Interfaces.C.size_t(nbchar-1))
        := C_dblarrs.Value(c,Interfaces.C.ptrdiff_t(nbchar));
    ind : Interfaces.C.size_t := 0;

  begin
    for i in res'range loop
      res(i) := Integer_to_Character(integer32(val(ind)));
      ind := ind + 1;
    end loop;
    return res;
  end Get_File_Name;

  procedure Assign_Count_and_Flags
              ( n : in integer32; cnt : in double_float;
                f : in Standard_Complex_VecMats.VecMat;
                c_d : in C_DblArrs.Pointer ) is

  -- DESCRIPTION :
  --   Assigns the real and imaginary parts of the n-dimensional flags 
  --   in row wise fashion to the C pointer, starting at position 1,
  --   because at position 0, the root count cnt will be placed.

    use Standard_Complex_Numbers;

    m : constant integer32 := f'last;
    size : constant integer32 := 2*m*n*n+1;
    val : C_Double_Array(0..Interfaces.C.size_t(size-1))
        := C_dblarrs.Value(c_d);
    ind : Interfaces.C.size_t := 0;
    cff : Complex_Number;

  begin
    val(ind) := Interfaces.C.double(cnt); ind := ind + 1;
    for k in 1..m loop
      for i in 1..n loop
        for j in 1..n loop
          cff := f(k)(i,j);
          val(ind) := Interfaces.C.double(REAL_PART(cff)); ind := ind + 1;
          val(ind) := Interfaces.C.double(IMAG_PART(cff)); ind := ind + 1;
        end loop;
      end loop;
    end loop;
    C_dblarrs.Copy_Array(val(0)'unchecked_access,c_d,
                         Interfaces.C.ptrdiff_t(size));
  end Assign_Count_and_Flags;

  procedure Assign_Count_and_Flags
              ( n : in integer32; cnt : in double_float;
                f : in DoblDobl_Complex_VecMats.VecMat;
                c_d : in C_DblArrs.Pointer ) is

  -- DESCRIPTION :
  --   Assigns the real and imaginary parts of the n-dimensional flags 
  --   in row wise fashion to the C pointer, starting at position 1,
  --   because at position 0, the root count cnt will be placed.

    use DoblDobl_Complex_Numbers;

    m : constant integer32 := f'last;
    size : constant integer32 := 4*m*n*n+1;
    val : C_Double_Array(0..Interfaces.C.size_t(size-1))
        := C_dblarrs.Value(c_d);
    ind : Interfaces.C.size_t := 0;
    cff : Complex_Number;
    ddf : double_double;

  begin
    val(ind) := Interfaces.C.double(cnt); ind := ind + 1;
    for k in 1..m loop
      for i in 1..n loop
        for j in 1..n loop
          cff := f(k)(i,j);
          ddf := REAL_PART(cff);
          val(ind) := Interfaces.C.double(hi_part(ddf)); ind := ind + 1;
          val(ind) := Interfaces.C.double(lo_part(ddf)); ind := ind + 1;
          ddf := IMAG_PART(cff);
          val(ind) := Interfaces.C.double(hi_part(ddf)); ind := ind + 1;
          val(ind) := Interfaces.C.double(lo_part(ddf)); ind := ind + 1;
        end loop;
      end loop;
    end loop;
    C_dblarrs.Copy_Array(val(0)'unchecked_access,c_d,
                         Interfaces.C.ptrdiff_t(size));
  end Assign_Count_and_Flags;

  procedure Assign_Count_and_Flags
              ( n : in integer32; cnt : in double_float;
                f : in QuadDobl_Complex_VecMats.VecMat;
                c_d : in C_DblArrs.Pointer ) is

  -- DESCRIPTION :
  --   Assigns the real and imaginary parts of the n-dimensional flags 
  --   in row wise fashion to the C pointer, starting at position 1,
  --   because at position 0, the root count cnt will be placed.

    use QuadDobl_Complex_Numbers;

    m : constant integer32 := f'last;
    size : constant integer32 := 8*m*n*n+1;
    val : C_Double_Array(0..Interfaces.C.size_t(size-1))
        := C_dblarrs.Value(c_d);
    ind : Interfaces.C.size_t := 0;
    cff : Complex_Number;
    qdf : Quad_double;

  begin
    val(ind) := Interfaces.C.double(cnt); ind := ind + 1;
    for k in 1..m loop
      for i in 1..n loop
        for j in 1..n loop
          cff := f(k)(i,j);
          qdf := REAL_PART(cff);
          val(ind) := Interfaces.C.double(hihi_part(qdf)); ind := ind + 1;
          val(ind) := Interfaces.C.double(lohi_part(qdf)); ind := ind + 1;
          val(ind) := Interfaces.C.double(hilo_part(qdf)); ind := ind + 1;
          val(ind) := Interfaces.C.double(lolo_part(qdf)); ind := ind + 1;
          qdf := IMAG_PART(cff);
          val(ind) := Interfaces.C.double(hihi_part(qdf)); ind := ind + 1;
          val(ind) := Interfaces.C.double(lohi_part(qdf)); ind := ind + 1;
          val(ind) := Interfaces.C.double(hilo_part(qdf)); ind := ind + 1;
          val(ind) := Interfaces.C.double(lolo_part(qdf)); ind := ind + 1;
        end loop;
      end loop;
    end loop;
    C_dblarrs.Copy_Array(val(0)'unchecked_access,c_d,
                         Interfaces.C.ptrdiff_t(size));
  end Assign_Count_and_Flags;

  function Job0 return integer32 is -- resolve Schubert intersection conditions

    n,k,nbc : integer32;
    otp : boolean;
    rc : Natural_Number;
    nrc : natural32;

  begin
    Get_Dimensions(a,n,k,nbc,otp);
   -- new_line;
   -- put_line("The dimensions : ");
   -- put("  n = "); put(n,1);
   -- put("  k = "); put(k,1);
   -- put("  c = "); put(nbc,1);
   -- if otp
   --  then put_line("  output wanted");
   --  else put_line("  in silent mode");
   -- end if;
    declare
      cond : constant Array_of_Brackets(1..nbc) := Get_Conditions(b,k,nbc);
    begin
     -- put_line("The brackets : ");
     -- for i in cond'range loop
     --   put(cond(i).all);
     -- end loop;
     -- new_line;
      Create_Intersection_Poset(n,nbc,cond,not otp,rc);
    end;
    nrc := Multprec_Natural_Numbers.Create(rc);
   -- put("The formal root count : "); put(nrc,1); new_line;
    Assign(double_float(nrc),c);
    return 0;
  end Job0;

  function Job1 return integer32 is -- standard double L-R homotopies

    use Standard_Complex_Solutions;
    use Standard_Solution_Posets;

    n,k,nbc,nbchar : integer32;
    otp,verify : boolean;
    rc : Natural_Number;
    nrc : natural32;
    tol : constant double_float := 1.0E-5;
    minrep : boolean := true;

  begin
    Get_Dimensions2(a,n,k,nbc,nbchar,otp,verify);
    verify := false;
   -- new_line;
   -- put_line("The dimensions : ");
   -- put("  n = "); put(n,1);
   -- put("  k = "); put(k,1);
   -- put("  c = "); put(nbc,1);
   -- if otp
   --  then put_line("  output wanted");
   --  else put_line("  in silent mode");
   -- end if;
   -- put("Number of characters : "); put(nbchar,1); new_line;
    declare
      cond : constant Array_of_Brackets(1..nbc) := Get_Conditions(b,k,nbc);
      q : constant Standard_Natural_Vectors.Vector
        := Checker_Moves.Identity_Permutation(natural32(n));
      rows,cols : Standard_Natural_Vectors.Vector(1..k);
      cnds : Standard_Natural_VecVecs.Link_to_VecVec;
      ips : Intersection_Poset(nbc-1) := Process_Conditions(n,k,nbc,cond);
      sps : Solution_Poset(ips.m) := Create(ips);
      flgs : Standard_Complex_VecMats.VecMat(1..nbc-2)
           := Random_Flags(n,nbc-2);
      name : constant string := Get_File_Name(c,nbchar);
      monitor : constant boolean := true;
      timer : Timing_Widget;
      file : file_type;
      sols : Solution_List;
      fsys : Standard_Complex_Poly_Systems.Link_to_Poly_Sys;
    begin
      put_line("The brackets : ");
      for i in cond'range loop
        put(cond(i).all);
      end loop;
      new_line;
      for i in 1..k loop
        rows(i) := cond(1)(i);
        cols(i) := cond(2)(i);
      end loop;
      cnds := new Standard_Natural_VecVecs.VecVec(1..nbc-2);
      for j in 1..nbc-2 loop 
        cnds(j) := new Standard_Natural_Vectors.Vector(1..k);
        for i in 1..k loop
          cnds(j)(i) := cond(j+2)(i);
        end loop;
      end loop;
      Communications_with_User.Create_Output_File(file,name);
      Count_Roots(file,ips,rc);
      put("the root count : "); put(rc,1); new_line;
      tstart(timer);
      Resolve(file,monitor,otp,n,k,tol,ips,sps,verify,minrep,
              cnds.all,flgs,sols);
      tstop(timer);
      Write_Results(file,n,k,q,rows,cols,minrep,cnds,flgs,sols,fsys);
      Standard_PolySys_Container.Initialize(fsys.all);
      Standard_Solutions_Container.Initialize(sols);
      new_line(file);
      print_times(file,timer,"resolving a Schubert problem");
      Close(file);
      nrc := Multprec_Natural_Numbers.Create(rc);
     -- put("The formal root count : "); put(nrc,1); new_line;
      Assign_Count_and_Flags(n,double_float(nrc),flgs,c);
    end;
    return 0;
  end Job1;

  function Job2 return integer32 is -- double double L-R homotopies

    use DoblDobl_Complex_Solutions;
    use DoblDobl_Solution_Posets;

    n,k,nbc,nbchar : integer32;
    otp,verify : boolean;
    rc : Natural_Number;
    nrc : natural32;
    tol : constant double_float := 1.0E-5;
    minrep : boolean := true;

  begin
    Get_Dimensions2(a,n,k,nbc,nbchar,otp,verify);
    verify := false;
   -- new_line;
   -- put_line("The dimensions : ");
   -- put("  n = "); put(n,1);
   -- put("  k = "); put(k,1);
   -- put("  c = "); put(nbc,1);
   -- if otp
   --  then put_line("  output wanted");
   --  else put_line("  in silent mode");
   -- end if;
   -- put("Number of characters : "); put(nbchar,1); new_line;
    declare
      cond : constant Array_of_Brackets(1..nbc) := Get_Conditions(b,k,nbc);
      q : constant Standard_Natural_Vectors.Vector
        := Checker_Moves.Identity_Permutation(natural32(n));
      rows,cols : Standard_Natural_Vectors.Vector(1..k);
      cnds : Standard_Natural_VecVecs.Link_to_VecVec;
      ips : Intersection_Poset(nbc-1) := Process_Conditions(n,k,nbc,cond);
      sps : Solution_Poset(ips.m) := Create(ips);
      flgs : DoblDobl_Complex_VecMats.VecMat(1..nbc-2)
           := Random_Flags(n,nbc-2);
      name : constant string := Get_File_Name(c,nbchar);
      monitor : constant boolean := true;
      timer : Timing_Widget;
      file : file_type;
      sols : Solution_List;
      fsys : DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
    begin
      put_line("The brackets : ");
      for i in cond'range loop
        put(cond(i).all);
      end loop;
      new_line;
      for i in 1..k loop
        rows(i) := cond(1)(i);
        cols(i) := cond(2)(i);
      end loop;
      cnds := new Standard_Natural_VecVecs.VecVec(1..nbc-2);
      for j in 1..nbc-2 loop 
        cnds(j) := new Standard_Natural_Vectors.Vector(1..k);
        for i in 1..k loop
          cnds(j)(i) := cond(j+2)(i);
        end loop;
      end loop;
      Communications_with_User.Create_Output_File(file,name);
      Count_Roots(file,ips,rc);
      put("the root count : "); put(rc,1); new_line;
      tstart(timer);
      Resolve(file,monitor,otp,n,k,tol,ips,sps,verify,minrep,
              cnds.all,flgs,sols);
      tstop(timer);
      Write_Results(file,n,k,q,rows,cols,minrep,cnds,flgs,sols,fsys);
      DoblDobl_PolySys_Container.Initialize(fsys.all);
      DoblDobl_Solutions_Container.Initialize(sols);
      new_line(file);
      print_times(file,timer,"resolving a Schubert problem");
      Close(file);
      nrc := Multprec_Natural_Numbers.Create(rc);
     -- put("The formal root count : "); put(nrc,1); new_line;
      Assign_Count_and_Flags(n,double_float(nrc),flgs,c);
    end;
    return 0;
  end Job2;

  function Job3 return integer32 is -- double double L-R homotopies

    use QuadDobl_Complex_Solutions;
    use QuadDobl_Solution_Posets;

    n,k,nbc,nbchar : integer32;
    otp,verify : boolean;
    rc : Natural_Number;
    nrc : natural32;
    tol : constant double_float := 1.0E-5;
    minrep : boolean := true;

  begin
    Get_Dimensions2(a,n,k,nbc,nbchar,otp,verify);
    verify := false;
   -- new_line;
   -- put_line("The dimensions : ");
   -- put("  n = "); put(n,1);
   -- put("  k = "); put(k,1);
   -- put("  c = "); put(nbc,1);
   -- if otp
   --  then put_line("  output wanted");
   --  else put_line("  in silent mode");
   -- end if;
   -- put("Number of characters : "); put(nbchar,1); new_line;
    declare
      cond : constant Array_of_Brackets(1..nbc) := Get_Conditions(b,k,nbc);
      q : constant Standard_Natural_Vectors.Vector
        := Checker_Moves.Identity_Permutation(natural32(n));
      rows,cols : Standard_Natural_Vectors.Vector(1..k);
      cnds : Standard_Natural_VecVecs.Link_to_VecVec;
      ips : Intersection_Poset(nbc-1) := Process_Conditions(n,k,nbc,cond);
      sps : Solution_Poset(ips.m) := Create(ips);
      flgs : QuadDobl_Complex_VecMats.VecMat(1..nbc-2)
           := Random_Flags(n,nbc-2);
      name : constant string := Get_File_Name(c,nbchar);
      monitor : constant boolean := true;
      timer : Timing_Widget;
      file : file_type;
      sols : Solution_List;
      fsys : QuadDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
    begin
      put_line("The brackets : ");
      for i in cond'range loop
        put(cond(i).all);
      end loop;
      new_line;
      for i in 1..k loop
        rows(i) := cond(1)(i);
        cols(i) := cond(2)(i);
      end loop;
      cnds := new Standard_Natural_VecVecs.VecVec(1..nbc-2);
      for j in 1..nbc-2 loop 
        cnds(j) := new Standard_Natural_Vectors.Vector(1..k);
        for i in 1..k loop
          cnds(j)(i) := cond(j+2)(i);
        end loop;
      end loop;
      Communications_with_User.Create_Output_File(file,name);
      Count_Roots(file,ips,rc);
      put("the root count : "); put(rc,1); new_line;
      tstart(timer);
      Resolve(file,monitor,otp,n,k,tol,ips,sps,verify,minrep,
              cnds.all,flgs,sols);
      tstop(timer);
      Write_Results(file,n,k,q,rows,cols,minrep,cnds,flgs,sols,fsys);
      QuadDobl_PolySys_Container.Initialize(fsys.all);
      QuadDobl_Solutions_Container.Initialize(sols);
      new_line(file);
      print_times(file,timer,"resolving a Schubert problem");
      Close(file);
      nrc := Multprec_Natural_Numbers.Create(rc);
     -- put("The formal root count : "); put(nrc,1); new_line;
      Assign_Count_and_Flags(n,double_float(nrc),flgs,c);
    end;
    return 0;
  end Job3;

  function Handle_Jobs return integer32 is
  begin
    case job is
      when 0 => return Job0; -- resolve Schubert intersection conditions
      when 1 => return Job1; -- standard Littlewood-Richardson homotopies
      when 2 => return Job2; -- dobldobl Littlewood-Richardson homotopies
      when 3 => return Job3; -- quaddobl Littlewood-Richardson homotopies
      when others => put_line("  Sorry.  Invalid operation in use_c2lrhom.");
                     return 1;
    end case;
  end Handle_Jobs;

begin
  return Handle_Jobs;
end use_c2lrhom;
