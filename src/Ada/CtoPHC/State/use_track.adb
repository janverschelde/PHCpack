with text_io;                           use text_io;
with Interfaces.C;                      use Interfaces.C;
with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;       use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;         use Standard_Floating_Numbers;
with Standard_Floating_Numbers_io;      use Standard_Floating_Numbers_io;
with Standard_Natural_Vectors;
with Standard_Complex_Poly_Systems;
with Standard_Complex_Poly_SysFun;      use Standard_Complex_Poly_SysFun;
with Standard_Complex_Jaco_Matrices;    use Standard_Complex_Jaco_Matrices;
with Standard_Complex_Solutions;
with Standard_Root_Refiners;            use Standard_Root_Refiners;
with Witness_Sets;
with Extrinsic_Diagonal_Homotopies;
with Standard_Hypersurface_Witdrivers;
with Assignments_in_Ada_and_C;          use Assignments_in_Ada_and_C;
with Assignments_of_Solutions;          use Assignments_of_Solutions;
with Standard_PolySys_Container;
with Standard_Solutions_Container;
with PHCpack_Operations;
with PHCpack_Operations_io;
with Crude_Path_Trackers;

with Diagonal_Homotopy_Interface;
with Path_Trackers_Interface;

function use_track ( job : integer32;
                     a : C_intarrs.Pointer;
                     b : C_intarrs.Pointer;
                     c : C_dblarrs.Pointer;
                     vrblvl : integer32 := 0 ) return integer32 is
 
  function JobM1 return integer32 is -- refine solution with Newton

    ls : constant Standard_Complex_Poly_Systems.Link_to_Poly_Sys
       := Standard_PolySys_Container.Retrieve;
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

  function Job8 return integer32 is -- write a string to defined output

    v_a : constant C_Integer_Array
        := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(1));
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

    v_a : constant C_Integer_Array
        := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(1));
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

    v_a : constant C_Integer_Array
        := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(1));
    n : constant natural32 := natural32(v_a(v_a'first));
    n1 : constant Interfaces.C.size_t := Interfaces.C.size_t(n-1);
    v_c : constant C_Double_Array(0..n1)
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

    v_a : constant C_Integer_Array
        := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(1));
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

    v_a : constant C_Integer_Array
        := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(1));
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

    v_a : constant C_Integer_Array
        := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(1));
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

  function Job16 return integer32 is -- read a witness set

    v_a : constant C_Integer_Array
        := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(1));
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

    v_a : constant C_Integer_Array
        := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(1));
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
   -- put("  n1 = "); put(n1,1);
   -- put("  n2 = "); put(n2,1);
   -- put("  a = "); put(a_d,1);
   -- put("  b = "); put(b_d,1); new_line;
    if a_d >= b_d then
      cd := Extrinsic_Diagonal_Homotopies.Cascade_Dimension(n1,n2,a_d,b_d);
    else
      cd := Extrinsic_Diagonal_Homotopies.Cascade_Dimension(n2,n1,b_d,a_d);
    end if;
   -- put("cascade dimension cd = "); put(cd,1); new_line;
    Assign(integer32(cd),a);
    return 0;
  exception
    when others =>
      put_line("Exception raised when computing cascade dimension.");
      return 18;
  end Job18;

  function Job19 return integer32 is -- computes witness set for hypersurface

    use Standard_Hypersurface_Witdrivers;

    v_a : constant C_Integer_Array
        := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(2));
    k : constant integer32 := integer32(v_a(v_a'first));
    n : constant integer := integer(v_a(v_a'first+1));
    lp : constant Standard_Complex_Poly_Systems.Link_to_Poly_Sys
       := Standard_PolySys_Container.Retrieve;
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

  function Job21 return integer32 is -- remove last slack variable

    v_a : constant C_Integer_Array
        := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(1));
    k : constant natural := natural(v_a(v_a'first));
    lp : constant Standard_Complex_Poly_Systems.Link_to_Poly_Sys
       := Standard_PolySys_Container.Retrieve;
    sols : Standard_Complex_Solutions.Solution_List
         := Standard_Solutions_Container.Retrieve;
    np : Standard_Complex_Poly_Systems.Poly_Sys(lp'first..lp'last-1);

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

 -- procedure Write_Symbols ( s : in Symbol_Table.Array_of_Symbols ) is

  -- DESCRIPTION :
  --   Writes the symbols in s to screen, useful for checking.

 -- begin
 --   for i in s'range loop
 --     put(" "); Symbol_Table_io.put(s(i));
 --   end loop;
 --   new_line;
 -- end Write_Symbols;

  function Job55 return integer32 is -- crude tracker in double precision

  -- DESCRIPTION :
  --   Gets the flag for the verbose mode from the first parameter a[0]
  --   and launches the crude path tracker in double precision.

    v_a : constant C_Integer_Array
        := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(1));
    val : constant natural32 := natural32(v_a(v_a'first));
    verbose : constant boolean := (val = 1);

  begin
    Crude_Path_Trackers.Standard_Track_Paths(verbose);
    return 0;
  exception
    when others =>
      put_line("Exception when launching crude tracker in double precision.");
      raise;
  end Job55;

  function Job56 return integer32 is -- crude tracker with double doubles

  -- DESCRIPTION :
  --   Gets the flag for the verbose mode from the first parameter a[0]
  --   and launches the crude path tracker in double double precision.

    v_a : constant C_Integer_Array
        := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(1));
    val : constant natural32 := natural32(v_a(v_a'first));
    verbose : constant boolean := (val = 1);

  begin
    Crude_Path_Trackers.DoblDobl_Track_Paths(verbose);
    return 0;
  exception
    when others =>
      put_line("Exception when launching crude tracker with double doubles.");
      raise;
  end Job56;

  function Job57 return integer32 is -- crude tracker with quad doubles

  -- DESCRIPTION :
  --   Gets the flag for the verbose mode from the first parameter a[0]
  --   and launches the crude path tracker in double double precision.

    v_a : constant C_Integer_Array
        := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(1));
    val : constant natural32 := natural32(v_a(v_a'first));
    verbose : constant boolean := (val = 1);

  begin
    Crude_Path_Trackers.QuadDobl_Track_Paths(verbose);
    return 0;
  exception
    when others =>
      put_line("Exception when launching crude tracker with double doubles.");
      raise;
  end Job57;

  function Handle_Jobs return integer32 is

    use Diagonal_Homotopy_Interface;
    use Path_Trackers_Interface;

  begin
    case job is
      when -1 => return JobM1; -- refine solution with Newton
      when 0 => PHCpack_Operations_io.Read_Target_System_without_Solutions;
                return 0;
      when 1 => PHCpack_Operations_io.Read_Start_System_without_Solutions;
                return 0;
      when 2 => return Path_Trackers_Standard_Homotopy_Random(vrblvl-1);
      when 3 => return Path_Trackers_Standard_Homotopy_Gamma(c,vrblvl-1);
      when 4 => return Path_Trackers_Standard_Homotopy_Clear(vrblvl-1);
      when 5 => return Path_Trackers_Standard_Silent_Track(a,b,c,vrblvl-1);
      when 6 => return Path_Trackers_Standard_Report_Track(a,b,c,vrblvl-1);
      when 7 => return Path_Trackers_Standard_Write_Solution(a,b,c,vrblvl-1);
      when 8 => return Job8;   -- write string to defined output file
      when 9 => return Job9;   -- write integers to defined output file
      when 10 => return Job10; -- write doubles to defined output file
      when 11 => return Job11; -- file name to read target system
      when 12 => return Job12; -- file name to read start system
      when 13 => return Job13; -- file name to read linear-product system
      when 14 => PHCpack_Operations.Standard_Cascade_Homotopy; return 0;
      when 15 => return Diagonal_Homotopy_Standard_Polynomial_Make(a,b,vrblvl);
      when 16 => return Job16; -- read a witness set
      when 17 => return Job17; -- reset input file for witness set k
      when 18 => return Job18; -- returns the extrinsic cascade dimension
      when 19 => return Job19; -- witness set for one polynomial
      when 20 => return Diagonal_Homotopy_Standard_Collapse(a,vrblvl);
      when 21 => return Job21; -- remove last slack variable
     -- tracking in double double precision :
      when 22 => return Path_Trackers_DoblDobl_Homotopy_Random(vrblvl-1);
      when 23 => return Path_Trackers_DoblDobl_Homotopy_Gamma(c,vrblvl-1);
      when 24 => return Path_Trackers_DoblDobl_Homotopy_Clear(vrblvl-1);
      when 25 => return Path_Trackers_DoblDobl_Silent_Track(a,b,c,vrblvl-1);
      when 26 => return Path_Trackers_DoblDobl_Report_Track(a,b,c,vrblvl-1);
      when 27 => return Path_Trackers_DoblDobl_Write_Solution(a,b,c,vrblvl-1);
      when 28 => PHCpack_Operations.DoblDobl_Cascade_Homotopy; return 0;
     -- tracking in quad double precision :
      when 32 => return Path_Trackers_QuadDobl_Homotopy_Random(vrblvl-1);
      when 33 => return Path_Trackers_QuadDobl_Homotopy_Gamma(c,vrblvl-1);
      when 34 => return Path_Trackers_QuadDobl_Homotopy_Clear(vrblvl-1);
      when 35 => return Path_Trackers_QuadDobl_Silent_Track(a,b,c,vrblvl-1);
      when 36 => return Path_Trackers_QuadDobl_Report_Track(a,b,c,vrblvl-1);
      when 37 => return Path_Trackers_QuadDobl_Write_Solution(a,b,c,vrblvl-1);
      when 38 => PHCpack_Operations.QuadDobl_Cascade_Homotopy; return 0;
     -- redefining diagonal homotopies ...
      when 40 => return Diagonal_Homotopy_Standard_Polynomial_Set(a,b,vrblvl);
      when 41 => return Diagonal_Homotopy_Standard_Start_Solutions(a,b,vrblvl);
      when 42 => return Diagonal_Homotopy_Symbols_Doubler(a,b,vrblvl);
     -- diagonal homotopy in double double and quad double precision
      when 43 => return Diagonal_Homotopy_DoblDobl_Polynomial_Make(a,b,vrblvl);
      when 44 => return Diagonal_Homotopy_QuadDobl_Polynomial_Make(a,b,vrblvl);
      when 45 => return Diagonal_Homotopy_DoblDobl_Start_Solutions(a,b,vrblvl);
      when 46 => return Diagonal_Homotopy_QuadDobl_Start_Solutions(a,b,vrblvl);
      when 47 => return Diagonal_Homotopy_DoblDobl_Collapse(a,vrblvl);
      when 48 => return Diagonal_Homotopy_QuadDobl_Collapse(a,vrblvl);
     -- double double and quad double witness sets for hypersurface
      when 49 => return Diagonal_Homotopy_DoblDobl_Polynomial_Set(a,b,vrblvl);
      when 50 => return Diagonal_Homotopy_QuadDobl_Polynomial_Set(a,b,vrblvl);
     -- multiprecision versions to create homotopy :
      when 52 => return Path_Trackers_Multprec_Homotopy_Random(vrblvl-1);
      when 53 => return Path_Trackers_Multprec_Homotopy_Gamma(c,vrblvl-1);
      when 54 => return Path_Trackers_Multprec_Homotopy_Clear(vrblvl-1);
     -- crude path trackers
      when 55 => return Job55; -- crude tracker in double precision
      when 56 => return Job56; -- crude tracker in double double precision
      when 57 => return Job57; -- crude tracker in quad double precision
     -- homotopies for Laurent systems
      when 58 =>
        PHCpack_Operations.Standard_Cascade_Laurent_Homotopy; return 0;
      when 59 =>
        PHCpack_Operations.DoblDobl_Cascade_Laurent_Homotopy; return 0;
      when 60 =>
        PHCpack_Operations.QuadDobl_Cascade_Laurent_Homotopy; return 0;
      when 61 => return Diagonal_Homotopy_Standard_Laurent_Make(a,b,vrblvl);
      when 62 => return Diagonal_Homotopy_DoblDobl_Laurent_Make(a,b,vrblvl);
      when 63 => return Diagonal_Homotopy_QuadDobl_Laurent_Make(a,b,vrblvl);
      when 64 => return Diagonal_Homotopy_Standard_Laurential_Set(a,b,vrblvl);
      when 65 => return Diagonal_Homotopy_DoblDobl_Laurential_Set(a,b,vrblvl);
      when 66 => return Diagonal_Homotopy_QuadDobl_Laurential_Set(a,b,vrblvl);
      when 67
        => PHCpack_Operations_io.Read_DoblDobl_Target_System_without_Solutions;
           return 0;
      when 68
        => PHCpack_Operations_io.Read_QuadDobl_Target_System_without_Solutions;
           return 0;
      when others => put_line("  Sorry.  Invalid operation."); return 1;
    end case;
  end Handle_Jobs;

begin
  return Handle_Jobs;
end use_track;
