with text_io;                           use text_io;
with Interfaces.C;                      use Interfaces.C;
with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;       use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;         use Standard_Floating_Numbers;
with Standard_Floating_Numbers_io;      use Standard_Floating_Numbers_io;
with Double_Double_Numbers;             use Double_Double_Numbers;
with Quad_Double_Numbers;               use Quad_Double_Numbers;
with Multprec_Floating_Numbers;
with Standard_Complex_Numbers;
with DoblDobl_Complex_Numbers;
with QuadDobl_Complex_Numbers;
with Multprec_Complex_Numbers;
with Standard_Natural_Vectors;
with Standard_Floating_Vectors;
with Symbol_Table;
with Standard_Complex_Polynomials;
with Standard_Complex_Poly_Strings;
with Standard_Complex_Poly_Systems;     use Standard_Complex_Poly_Systems;
with Standard_Complex_Poly_SysFun;      use Standard_Complex_Poly_SysFun;
with Standard_Complex_Jaco_Matrices;    use Standard_Complex_Jaco_Matrices;
with Standard_Complex_Solutions;
with DoblDobl_Complex_Solutions;
with QuadDobl_Complex_Solutions;
with Standard_Root_Refiners;            use Standard_Root_Refiners;
with Standard_Homotopy;
with Standard_Continuation_Data_io;
with DoblDobl_Homotopy;
with DoblDobl_Continuation_Data_io;
with QuadDobl_Homotopy;
with QuadDobl_Continuation_Data_io;
with Multprec_Homotopy;
with Witness_Sets,Witness_Sets_io;
with Extrinsic_Diagonal_Homotopies;
with Extrinsic_Diagonal_Solvers;
with Hypersurface_Witness_Sets_io;      use Hypersurface_Witness_Sets_io;
with Assignments_in_Ada_and_C;          use Assignments_in_Ada_and_C;
with Assignments_of_Solutions;          use Assignments_of_Solutions;
with Standard_PolySys_Container;
with Standard_Solutions_Container;
with PHCpack_Operations;
with PHCpack_Operations_io;

function use_track ( job : integer32;
                     a : C_intarrs.Pointer;
                     b : C_intarrs.Pointer;
                     c : C_dblarrs.Pointer ) return integer32 is
 
  function JobM1 return integer32 is -- refine solution with Newton

    ls : constant Link_to_Poly_Sys := Standard_PolySys_Container.Retrieve;
    ep : constant Link_to_Eval_Poly_Sys
       := Standard_PolySys_Container.Evaluator;
    jf : constant Link_to_Eval_Jaco_Mat
       := Standard_PolySys_Container.Jacobian_Evaluator;
    sol : Standard_Complex_Solutions.Solution(ls'last)
        := Convert_to_Solution(b,c);
    epsxa : constant double_float := 1.0E-14;
    epsfa : constant double_float := 1.0E-14;
    max : constant natural32 := 3;
    numit : natural32 := 0;
    fail : boolean;
   -- x : Standard_Complex_Vectors.Vector(ls'range)
   --   := (ls'range => Create(1.0));
   -- y : Standard_Complex_Vectors.Vector(ls'range);
   -- r : constant integer := Standard_Random_Numbers.Random(0,1000);

  begin
   -- put("starting evaluation with id = "); put(r,1); new_line;
   -- for i in 1..1000 loop
   --   y := Eval(ls.all,x);
   --   -- y := Eval(ls.all,sol.v); --y := Eval(ep.all,sol.v);
   -- end loop;
   -- put("ending evaluation with id = "); put(r,1); new_line;
    Silent_Newton(ep.all,jf.all,sol,epsxa,epsfa,numit,max,fail);
    Assign_Solution(sol,b,c);
    return 0;
  exception
    when others => put_line("exception occurred in root refiner...");
                   return 149;
  end JobM1;

  function Job3 return integer32 is -- standard homotopy with given gamma

    g : Standard_Floating_Vectors.Vector(1..2);
    gamma : Standard_Complex_Numbers.Complex_Number;

  begin
    Assign(2,c,g);
    gamma := Standard_Complex_Numbers.Create(g(1),g(2));
    PHCpack_Operations.Create_Standard_Homotopy(gamma);
    return 0;
  end Job3;

  function Job23 return integer32 is -- create homotopy with given gamma

    g : Standard_Floating_Vectors.Vector(1..2);
    g_re,g_im : double_double;
    gamma : DoblDobl_Complex_Numbers.Complex_Number;

  begin
    Assign(2,c,g);
    g_re := create(g(1));
    g_im := create(g(2));
    gamma := DoblDobl_Complex_Numbers.Create(g_re,g_im);
    PHCpack_Operations.Create_DoblDobl_Homotopy(gamma);
    return 0;
  end Job23;

  function Job33 return integer32 is -- quaddobl homotopy with given gamma

    g : Standard_Floating_Vectors.Vector(1..2);
    g_re,g_im : quad_double;
    gamma : QuadDobl_Complex_Numbers.Complex_Number;

  begin
    Assign(2,c,g);
    g_re := create(g(1));
    g_im := create(g(1));
    gamma := QuadDobl_Complex_Numbers.Create(g_re,g_im);
    PHCpack_Operations.Create_QuadDobl_Homotopy(gamma);
    return 0;
  end Job33;

  function Job53 return integer32 is -- multprec homotopy with given gamma

    g : Standard_Floating_Vectors.Vector(1..2);
    g_re,g_im : Multprec_Floating_Numbers.Floating_Number;
    gamma : Multprec_Complex_Numbers.Complex_Number;

  begin
    Assign(2,c,g);
    g_re := Multprec_Floating_Numbers.create(g(1));
    g_im := Multprec_Floating_Numbers.create(g(1));
    gamma := Multprec_Complex_Numbers.Create(g_re,g_im);
    Multprec_Floating_Numbers.Clear(g_re);
    Multprec_Floating_Numbers.Clear(g_im);
    PHCpack_Operations.Create_Multprec_Homotopy(gamma);
    return 0;
  end Job53;

  function Job5 return integer32 is -- track one path silently

    ls : Standard_Complex_Solutions.Link_to_Solution
       := Convert_to_Solution(b,c);
    length : double_float;
    nbstep,nbfail,nbiter,nbsyst : natural32;
    nbs : Standard_Natural_Vectors.Vector(1..4);
    crash : boolean;

  begin
    PHCpack_Operations.Silent_Path_Tracker
      (ls,length,nbstep,nbfail,nbiter,nbsyst,crash);
    nbs(1) := nbstep; nbs(2) := nbfail;
    nbs(3) := nbiter; nbs(4) := nbsyst;
    Assign(nbs,a);
    Assign_Solution(ls,b,c);
    Standard_Complex_Solutions.Clear(ls);
    if crash
     then return 5;
     else return 0;
    end if;
  end Job5;

  function Job25 return integer32 is -- track one path silently

    ls : DoblDobl_Complex_Solutions.Link_to_Solution
       := Convert_to_Solution(b,c);
    length : double_float;
    nbstep,nbfail,nbiter,nbsyst : natural32;
    nbs : Standard_Natural_Vectors.Vector(1..4);
    crash : boolean;

  begin
    PHCpack_Operations.Silent_Path_Tracker
      (ls,length,nbstep,nbfail,nbiter,nbsyst,crash);
    nbs(1) := nbstep; nbs(2) := nbfail;
    nbs(3) := nbiter; nbs(4) := nbsyst;
    Assign(nbs,a);
    Assign_Solution(ls,b,c);
    DoblDobl_Complex_Solutions.Clear(ls);
    if crash
     then return 25;
     else return 0;
    end if;
  end Job25;

  function Job35 return integer32 is -- track one path silently

    ls : QuadDobl_Complex_Solutions.Link_to_Solution
       := Convert_to_Solution(b,c);
    length : double_float;
    nbstep,nbfail,nbiter,nbsyst : natural32;
    nbs : Standard_Natural_Vectors.Vector(1..4);
    crash : boolean;

  begin
    PHCpack_Operations.Silent_Path_Tracker
      (ls,length,nbstep,nbfail,nbiter,nbsyst,crash);
    nbs(1) := nbstep; nbs(2) := nbfail;
    nbs(3) := nbiter; nbs(4) := nbsyst;
    Assign(nbs,a);
    Assign_Solution(ls,b,c);
    QuadDobl_Complex_Solutions.Clear(ls);
    if crash
     then return 25;
     else return 0;
    end if;
  end Job35;

  function Job6 return integer32 is -- track one path with reporting

    ls : Standard_Complex_Solutions.Link_to_Solution
       := Convert_to_Solution(b,c);
    length : double_float;
    nbstep,nbfail,nbiter,nbsyst : natural32;
    nbs : Standard_Natural_Vectors.Vector(1..4);
    crash : boolean;

  begin
    PHCpack_Operations.Reporting_Path_Tracker
      (ls,length,nbstep,nbfail,nbiter,nbsyst,crash);
    nbs(1) := nbstep; nbs(2) := nbfail;
    nbs(3) := nbiter; nbs(4) := nbsyst;
    Assign(nbs,a);
    Assign_Solution(ls,b,c);
    Standard_Complex_Solutions.Clear(ls);
    if crash
     then return 6;
     else return 0;
    end if;
  end Job6;

  function Job26 return integer32 is -- track one path with reporting

    ls : DoblDobl_Complex_Solutions.Link_to_Solution
       := Convert_to_Solution(b,c);
    length : double_float;
    nbstep,nbfail,nbiter,nbsyst : natural32;
    nbs : Standard_Natural_Vectors.Vector(1..4);
    crash : boolean;

  begin
    PHCpack_Operations.Reporting_Path_Tracker
      (ls,length,nbstep,nbfail,nbiter,nbsyst,crash);
    nbs(1) := nbstep; nbs(2) := nbfail;
    nbs(3) := nbiter; nbs(4) := nbsyst;
    Assign(nbs,a);
    Assign_Solution(ls,b,c);
    DoblDobl_Complex_Solutions.Clear(ls);
    if crash
     then return 26;
     else return 0;
    end if;
  end Job26;

  function Job36 return integer32 is -- track one path with reporting

    ls : QuadDobl_Complex_Solutions.Link_to_Solution
       := Convert_to_Solution(b,c);
    length : double_float;
    nbstep,nbfail,nbiter,nbsyst : natural32;
    nbs : Standard_Natural_Vectors.Vector(1..4);
    crash : boolean;

  begin
    PHCpack_Operations.Reporting_Path_Tracker
      (ls,length,nbstep,nbfail,nbiter,nbsyst,crash);
    nbs(1) := nbstep; nbs(2) := nbfail;
    nbs(3) := nbiter; nbs(4) := nbsyst;
    Assign(nbs,a);
    Assign_Solution(ls,b,c);
    QuadDobl_Complex_Solutions.Clear(ls);
    if crash
     then return 36;
     else return 0;
    end if;
  end Job36;

  function Job7 return integer32 is -- write solution with diagnostics

    use Standard_Complex_Numbers,Standard_Continuation_Data_io;

    ls : constant Standard_Complex_Solutions.Link_to_Solution
       := Convert_to_Solution(b,c);
    nb : Standard_Natural_Vectors.Vector(1..5);
    length_path : constant double_float := IMAG_PART(ls.t);
    nbfail,nbregu,nbsing,kind : natural32 := 0;
    tol_zero : constant double_float := 1.0E-12;
    tol_sing : constant double_float := 1.0E-8;


  begin
    Assign(5,a,nb);
    if PHCpack_Operations.Is_File_Defined then
      Write_Statistics(PHCpack_Operations.output_file,
                       integer32(nb(5)),nb(1),nb(2),nb(3),nb(4));
      Write_Diagnostics(PHCpack_Operations.output_file,
                        ls.all,tol_zero,tol_sing,nbfail,nbregu,nbsing,kind);
    else
      Write_Statistics(standard_output,
                       integer32(nb(5)),nb(1),nb(2),nb(3),nb(4));
      Write_Diagnostics(standard_output,
                        ls.all,tol_zero,tol_sing,nbfail,nbregu,nbsing,kind);
    end if;
    ls.t := Create(REAL_PART(ls.t),0.0);
    if PHCpack_Operations.Is_File_Defined
     then Write_Solution(PHCpack_Operations.output_file,ls.all,length_path);
     else Write_Solution(standard_output,ls.all,length_path);
    end if;
    return 0;
  end Job7;

  function Job27 return integer32 is -- write solution with diagnostics

    use DoblDobl_Complex_Numbers,DoblDobl_Continuation_Data_io;

    ls : constant DoblDobl_Complex_Solutions.Link_to_Solution
       := Convert_to_Solution(b,c);
    nb : Standard_Natural_Vectors.Vector(1..5);
    length_path : constant double_double := IMAG_PART(ls.t);
    nbfail,nbregu,nbsing,kind : natural32 := 0;
    tol_zero : constant double_float := 1.0E-12;
    tol_sing : constant double_float := 1.0E-8;


  begin
    Assign(5,a,nb);
    if PHCpack_Operations.Is_File_Defined then
      Write_Statistics(PHCpack_Operations.output_file,
                       nb(5),nb(1),nb(2),nb(3),nb(4));
      Write_Diagnostics(PHCpack_Operations.output_file,
                        ls.all,tol_zero,tol_sing,nbfail,nbregu,nbsing,kind);
    else
      Write_Statistics(standard_output,nb(5),nb(1),nb(2),nb(3),nb(4));
      Write_Diagnostics(standard_output,
                        ls.all,tol_zero,tol_sing,nbfail,nbregu,nbsing,kind);
    end if;
    ls.t := Create(REAL_PART(ls.t),Double_Double_Numbers.create(0.0));
    if PHCpack_Operations.Is_File_Defined
     then Write_Solution(PHCpack_Operations.output_file,ls.all,length_path);
     else Write_Solution(standard_output,ls.all,length_path);
    end if;
    return 0;
  end Job27;

  function Job37 return integer32 is -- write solution with diagnostics

    use QuadDobl_Complex_Numbers,QuadDobl_Continuation_Data_io;

    ls : constant QuadDobl_Complex_Solutions.Link_to_Solution
       := Convert_to_Solution(b,c);
    nb : Standard_Natural_Vectors.Vector(1..5);
    length_path : constant quad_double := IMAG_PART(ls.t);
    nbfail,nbregu,nbsing,kind : natural32 := 0;
    tol_zero : constant double_float := 1.0E-12;
    tol_sing : constant double_float := 1.0E-8;


  begin
    Assign(5,a,nb);
    if PHCpack_Operations.Is_File_Defined then
      Write_Statistics(PHCpack_Operations.output_file,
                       nb(5),nb(1),nb(2),nb(3),nb(4));
      Write_Diagnostics(PHCpack_Operations.output_file,
                        ls.all,tol_zero,tol_sing,nbfail,nbregu,nbsing,kind);
    else
      Write_Statistics(standard_output,nb(5),nb(1),nb(2),nb(3),nb(4));
      Write_Diagnostics(standard_output,
                        ls.all,tol_zero,tol_sing,nbfail,nbregu,nbsing,kind);
    end if;
    ls.t := Create(REAL_PART(ls.t),Quad_Double_Numbers.create(integer(0)));
    if PHCpack_Operations.Is_File_Defined
     then Write_Solution(PHCpack_Operations.output_file,ls.all,length_path);
     else Write_Solution(standard_output,ls.all,length_path);
    end if;
    return 0;
  end Job37;

  function Job8 return integer32 is -- write a string to defined output

    v_a : constant C_Integer_Array := C_intarrs.Value(a);
    n : constant integer := integer(v_a(v_a'first));
    n1 : constant Interfaces.C.size_t := Interfaces.C.size_t(n-1);
    v_b : constant C_Integer_Array(0..n1)
        := C_Intarrs.Value(b,Interfaces.C.ptrdiff_t(n));
    s : constant String(1..n) := C_Integer_Array_to_String(natural32(n),v_b);

  begin
    if PHCpack_Operations.Is_File_Defined
     then put(PHCpack_Operations.output_file,s);
     else put(standard_output,s);
    end if;
    return 0;
  end Job8;

  function Job9 return integer32 is -- writes integers to defined output

    v_a : constant C_Integer_Array := C_intarrs.Value(a);
    n : constant natural32 := natural32(v_a(v_a'first));
    n1 : constant Interfaces.C.size_t := Interfaces.C.size_t(n-1);
    v_b : constant C_Integer_Array(0..n1)
        := C_Intarrs.Value(b,Interfaces.C.ptrdiff_t(n));

   -- use Interfaces.C;

  begin
    if PHCpack_Operations.Is_File_Defined then
      put(PHCpack_Operations.output_file,integer32(v_b(v_b'first)),1);
      for i in v_b'first+1..v_b'last loop
        put(PHCpack_Operations.output_file," ");
        put(PHCpack_Operations.output_file,integer32(v_b(i)),1);
      end loop;
    else
      put(standard_output,integer32(v_b(v_b'first)),1);
      for i in v_b'first+1..v_b'last loop
        put(standard_output," ");
        put(standard_output,integer32(v_b(i)),1);
      end loop;
    end if;
    return 0;
  end Job9;

  function Job10 return integer32 is -- writes doubles to defined output

    v_a : constant C_Integer_Array := C_intarrs.Value(a);
    n : constant natural32 := natural32(v_a(v_a'first));
    n1 : constant Interfaces.C.size_t := Interfaces.C.size_t(n-1);
    v_c : C_Double_Array(0..n1)
        := C_dblarrs.Value(c,Interfaces.C.ptrdiff_t(n));
    d : double_float;

   -- use Interfaces.C;

  begin
    d := double_float(v_c(v_c'first));
    if PHCpack_Operations.Is_File_Defined then
      put(PHCpack_Operations.output_file,d);
      for i in v_c'first+1..v_c'last loop
        put(PHCpack_Operations.output_file," ");
        d := double_float(v_c(i));
        put(PHCpack_Operations.output_file,d);
      end loop;
    else
      put(standard_output,d);
      for i in v_c'first+1..v_c'last loop
        put(standard_output," ");
        d := double_float(v_c(i));
        put(standard_output,d);
      end loop;
    end if;
    return 0;
  end Job10;

  function Job11 return integer32 is -- file name to read target system

    v_a : constant C_Integer_Array := C_intarrs.Value(a);
    n : constant integer := integer(v_a(v_a'first));
    n1 : constant Interfaces.C.size_t := Interfaces.C.size_t(n-1);
    v_b : constant C_Integer_Array(0..n1)
        := C_intarrs.Value(b,Interfaces.C.ptrdiff_t(n));
    s : constant String(1..n) := C_Integer_Array_to_String(natural32(n),v_b);

  begin
   -- put_line("opening the file " & s & " for the target system ...");
    PHCpack_Operations_io.Read_Target_System_without_Solutions(s);
    return 0;
  exception
    when others =>
      put_line("Exception raised when opening " & s & " for target system.");
      return 161;
  end Job11;

  function Job12 return integer32 is -- file name to read start system

    v_a : constant C_Integer_Array := C_intarrs.Value(a);
    n : constant integer := integer(v_a(v_a'first));
    n1 : constant Interfaces.C.size_t := Interfaces.C.size_t(n-1);
    v_b : constant C_Integer_Array(0..n1)
        := C_intarrs.Value(b,Interfaces.C.ptrdiff_t(n));
    s : constant String(1..n) := C_Integer_Array_to_String(natural32(n),v_b);

  begin
   -- put_line("opening the file " & s & " for the start system ...");
    PHCpack_Operations_io.Read_Start_System_without_Solutions(s);
    return 0;
  exception
    when others =>
      put_line("Exception raised when opening " & s & " for start system.");
      return 12;
  end Job12;

  function Job13 return integer32 is -- name to read linear-product system

    v_a : constant C_Integer_Array := C_intarrs.Value(a);
    n : constant integer := integer(v_a(v_a'first));
    n1 : constant Interfaces.C.size_t := Interfaces.C.size_t(n-1);
    v_b : constant C_Integer_Array(0..n1)
        := C_intarrs.Value(b,Interfaces.C.ptrdiff_t(n));
    s : constant String(1..n) := C_Integer_Array_to_String(natural32(n),v_b);
    fail : boolean;

  begin
   -- put_line("opening the file " & s & " for the start system ...");
    PHCpack_Operations_io.Read_Linear_Product_Start_System(s,fail);
    if fail
     then return 163;
     else return 0;
    end if;
  exception
    when others =>
      put_line("Exception raised when opening " & s & " for start system.");
      return 13;
  end Job13;

  function Job15 return integer32 is -- create a diagonal homotopy

    v_a : constant C_Integer_Array := C_intarrs.Value(a);
    v_b : constant C_Integer_Array := C_intarrs.Value(b);
    a_dim : constant natural32 := natural32(v_a(v_a'first));
    b_dim : constant natural32 := natural32(v_b(v_b'first));

  begin
    PHCpack_Operations.Create_Diagonal_Homotopy(a_dim,b_dim);
    return 0;
  exception
    when others =>
      put_line("Exception raised when creating a diagonal homotopy.");
      return 15;
  end Job15;

  function Job16 return integer32 is -- read a witness set

    v_a : constant C_Integer_Array := C_intarrs.Value(a);
    k : constant natural32 := natural32(v_a(v_a'first));
    n,dim,deg : natural32;
    nbs : Standard_Natural_Vectors.Vector(1..2);
    fail : boolean;

  begin
    if k = 1 or k = 2 then
      PHCpack_Operations_io.Read_Witness_Set_for_Diagonal_Homotopy
        (k,n,dim,deg,fail);
      if fail then
        return 16;
      else
        Assign(integer32(n),a);
        nbs(1) := dim; nbs(2) := deg;
        Assign(nbs,b);
      end if;
    else
      put("Wrong value on input : "); put(k,1); new_line;
      return 16;
    end if;
    return 0;
  exception
    when others =>
      put_line("Exception raised when reading a witness set.");
      return 16;
  end Job16;

  function Job17 return integer32 is -- reset input file for witness set k

    v_a : constant C_Integer_Array := C_intarrs.Value(a);
    k : constant natural32 := natural32(v_a(v_a'first));
    deg,dim : natural32;
    nbs : Standard_Natural_Vectors.Vector(1..2);
    fail : boolean;

  begin
   -- put("resetting the input file for witness set "); put(k,1); new_line;
    PHCpack_Operations_io.Reset_Witness_Input_File(k,deg,dim,fail);
    if fail  then
      return 17;
    else
      -- put("  degree : "); put(deg,1);
      -- put("  n : "); put(dim,1); new_line;
      nbs(1) := deg; nbs(2) := dim;
      Assign(nbs,b);
    end if;
    return 0;
  exception
    when others =>
      put("Exception raised when resetting input file for witness set ");
      put(k,1); put_line(".");
      return 17;
  end Job17;

  function Job18 return integer32 is -- returns extrinsic cascade dimension

    v_a : constant C_Integer_Array
        := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(2));
    n1 : constant natural32 := natural32(v_a(v_a'first));
    n2 : constant natural32 := natural32(v_a(v_a'first+1));
    v_b : constant C_Integer_Array
        := C_intarrs.Value(b,Interfaces.C.ptrdiff_t(2));
    a_d : constant natural32 := natural32(v_b(v_b'first));
    b_d : constant natural32 := natural32(v_b(v_b'first+1));
    cd : natural32;

  begin
    if a_d >= b_d then
      cd := Extrinsic_Diagonal_Homotopies.Cascade_Dimension(n1,n2,a_d,b_d);
    else
      cd := Extrinsic_Diagonal_Homotopies.Cascade_Dimension(n2,n1,b_d,a_d);
    end if;
    Assign(integer32(cd),a);
    return 0;
  exception
    when others =>
      put_line("Exception raised when computing cascade dimension.");
      return 18;
  end Job18;

  function Job19 return integer32 is -- computes witness set for hypersurface

    v_a : constant C_Integer_Array
        := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(2));
    k : constant integer32 := integer32(v_a(v_a'first));
    n : constant integer := integer(v_a(v_a'first+1));
    lp : constant Link_to_Poly_Sys := Standard_PolySys_Container.Retrieve;
    n1 : constant Interfaces.C.size_t := Interfaces.C.size_t(n-1);
    v_b : constant C_Integer_Array(0..n1)
        := C_intarrs.Value(b,Interfaces.C.ptrdiff_t(n));
    filename : constant String(1..n)
             := C_Integer_Array_to_String(natural32(n),v_b);
    file : file_type;
    fail : boolean;

  begin
   -- put("The number of the equation : "); put(k,1); new_line;
   -- if lp = null then
   --   put_line("But the systems container is empty!");
   -- elsif lp'last < k then
   --   put("But there are only "); put(lp'last,1);
   --   put_line(" polynomials in container!");
   -- else
   --   put("Polynomial : "); put(lp(k)); new_line;
   --   put_line("Creating the output file with name " & filename);
      Create(file,out_file,filename);
      Call_Root_Finder(file,lp(k),false,1.0E-10,fail);
      Close(file);
   -- end if;
    return 0;
  exception
    when others =>
      put_line("Exception raised at hypersurface witness set computation.");
      return 19;
  end Job19;

  function Job20 return integer32 is -- collapses extrinsic diagonal

    v_a : constant C_Integer_Array
        := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(2));
    k : constant natural32 := natural32(v_a(v_a'first));
    d : constant natural32 := natural32(v_a(v_a'first+1));
    lp : constant Link_to_Poly_Sys := Standard_PolySys_Container.Retrieve;
    sols : constant Standard_Complex_Solutions.Solution_List
         := Standard_Solutions_Container.Retrieve;
    clps : Standard_Complex_Solutions.Solution_List;
    cp : Link_to_Poly_Sys;

  begin
   -- put("#equations in systems container : "); put(lp'last,1); new_line;
   -- put("#solutions in solutions container : ");
   -- put(Length_Of(sols),1); new_line;
    Standard_Complex_Solutions.Copy(sols,clps);
    Extrinsic_Diagonal_Solvers.Collapse_System(lp.all,clps,k,d,cp);
    Standard_PolySys_Container.Clear;
    Standard_PolySys_Container.Initialize(cp.all);
    Standard_Solutions_Container.Clear;
    Standard_Solutions_Container.Initialize(clps);
    return 0;
  exception 
     when others =>
       put_line("Exception raised at collapsing of extrinsic diagonal.");
       return 20;
  end Job20;

  function Job21 return integer32 is -- remove last slack variable

    v_a : constant C_Integer_Array
        := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(1));
    k : constant natural := natural(v_a(v_a'first));
    lp : constant Link_to_Poly_Sys := Standard_PolySys_Container.Retrieve;
    sols : Standard_Complex_Solutions.Solution_List
         := Standard_Solutions_Container.Retrieve;
    np : Poly_Sys(lp'first..lp'last-1);

  begin
   -- put("#equations in systems container : "); put(lp'last,1); new_line;
   -- put("#solutions in solutions container : ");
   -- put(Length_Of(sols),1); new_line;
    if k > 0 then
      Witness_Sets.Remove_Component(sols);
      np := Witness_Sets.Remove_Embedding1(lp.all,1);
      Standard_PolySys_Container.Clear;
      Standard_PolySys_Container.Initialize(np);
    end if;
    return 0;
  exception
    when others =>
      put_line("Exception raised when removing last slack variable.");
      return 21;
  end Job21;

  function Job40 return integer32 is -- witness set of hypersurface

    v_a : constant C_Integer_Array
        := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(2));
    nv : constant natural32 := natural32(v_a(v_a'first));
    nc : constant integer := integer(v_a(v_a'first+1));
    nc1 : constant Interfaces.C.size_t := Interfaces.C.size_t(nc-1);
    v_b : constant C_Integer_Array(0..nc1)
        := C_Intarrs.Value(b,Interfaces.C.ptrdiff_t(nc));
    s : constant String(1..nc) := C_Integer_Array_to_String(natural32(nc),v_b);
    p : Standard_Complex_Polynomials.Poly;
    eps : constant double_float := 1.0E-12;
    fail : boolean;
    e : Standard_Complex_Poly_Systems.Link_to_Poly_Sys;
    esols : Standard_Complex_Solutions.Solution_List;

  begin
    if Symbol_Table.Empty
     then Symbol_Table.Init(nv);
    end if;
    p := Standard_Complex_Poly_Strings.Parse(nv,s);
    Silent_Root_Finder(p,eps,fail,e,esols);
   -- if fail
   --  then put_line("a failure occurred");
   --  else put_line("no failure occurred");
   -- end if;
    Witness_Sets_io.Add_Embed_Symbols(nv-1);
    Standard_PolySys_Container.Clear;
    Standard_PolySys_Container.Initialize(e.all);
    Standard_Solutions_Container.Clear;
    Standard_Solutions_Container.Initialize(esols);
    return 0;
  exception
    when others =>
      put_line("Exception raised when making witness set for hypersurface.");
      return 270;
  end Job40;

  function Job41 return integer32 is -- solutions to start diagonal cascade

    v_a : constant C_Integer_Array := C_intarrs.Value(a);
    v_b : constant C_Integer_Array := C_intarrs.Value(b);
    a_dim : constant natural32 := natural32(v_a(v_a'first));
    b_dim : constant natural32 := natural32(v_b(v_b'first));

  begin
    PHCpack_Operations.Start_Diagonal_Cascade_Solutions(a_dim,b_dim);
    return 0;
  exception
    when others =>
      put_line("Exception raised when making solutions to start a cascade.");
      return 271;
  end Job41;

  function Handle_Jobs return integer32 is
  begin
    case job is
      when -1 => return JobM1; -- refine solution with Newton
      when 0 => PHCpack_Operations_io.Read_Target_System_without_Solutions;
                return 0;
      when 1 => PHCpack_Operations_io.Read_Start_System_without_Solutions;
                return 0;
      when 2 => PHCpack_Operations.Create_Standard_Homotopy; return 0;
      when 3 => return Job3;   -- create homotopy with given gamma
      when 4 => Standard_Homotopy.Clear; return 0;
      when 5 => return Job5;   -- track one path silently
      when 6 => return Job6;   -- track one path with reporting
      when 7 => return Job7;   -- write next solution with diagnostics
      when 8 => return Job8;   -- write string to defined output file
      when 9 => return Job9;   -- write integers to defined output file
      when 10 => return Job10; -- write doubles to defined output file
      when 11 => return Job11; -- file name to read target system
      when 12 => return Job12; -- file name to read start system
      when 13 => return Job13; -- file name to read linear-product system
      when 14 => PHCpack_Operations.Create_Cascade_Homotopy; return 0;
      when 15 => return Job15; -- create a diagonal homotopy
      when 16 => return Job16; -- read a witness set
      when 17 => return Job17; -- reset input file for witness set k
      when 18 => return Job18; -- returns the extrinsic cascade dimension
      when 19 => return Job19; -- witness set for one polynomial
      when 20 => return Job20; -- collapse extrinsic diagonal
      when 21 => return Job21; -- remove last slack variable
     -- tracking in double double precision :
      when 22 => PHCpack_Operations.Create_DoblDobl_Homotopy; return 0;
      when 23 => return Job23; -- double double homotopy with given gamma
      when 24 => DoblDobl_Homotopy.Clear; return 0;
      when 25 => return Job25; -- track one path silently
      when 26 => return Job26; -- track one path with reporting
      when 27 => return Job27; -- write solution with diagnostics
      when 28 => PHCpack_Operations.Create_DoblDobl_Cascade_Homotopy; return 0;
     -- tracking in quad double precision :
      when 32 => PHCpack_Operations.Create_QuadDobl_Homotopy; return 0;
      when 33 => return Job33; -- quad double homotopy with given gamma
      when 34 => QuadDobl_Homotopy.Clear; return 0;
      when 35 => return Job35; -- track one path silently
      when 36 => return Job36; -- track one path with reporting
      when 37 => return Job37; -- write solution with diagnostics
      when 38 => PHCpack_Operations.Create_QuadDobl_Cascade_Homotopy; return 0;
     -- redefining diagonal homotopies ...
      when 40 => return Job40; -- witness set of hypersurface
      when 41 => return Job41; -- solutions to start diagonal cascade
     -- multiprecision versions to create homotopy :
      when 52 => PHCpack_Operations.Create_Multprec_Homotopy; return 0;
      when 53 => return Job53; -- multiprecision homotopy with given gamma
      when 54 => Multprec_Homotopy.Clear; return 0;
      when others => put_line("  Sorry.  Invalid operation."); return 1;
    end case;
  end Handle_Jobs;

begin
  return Handle_Jobs;
end use_track;
