with text_io;                           use text_io;
with Interfaces.C;                      use Interfaces.C;
with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;         use Standard_Floating_Numbers;
with Standard_Floating_Numbers_io;      use Standard_Floating_Numbers_io;
with Standard_Complex_Numbers;
with Double_Double_Numbers;             use Double_Double_Numbers;
with Quad_Double_Numbers;               use Quad_Double_Numbers;
with DoblDobl_Complex_Numbers;
with QuadDobl_Complex_Numbers;
with Standard_Complex_Vectors;
with Standard_Complex_Vectors_io;
with Standard_Complex_VecVecs;
with DoblDobl_Complex_Vectors;
with DoblDobl_Complex_Vectors_io;
with DoblDobl_Complex_VecVecs;
with QuadDobl_Complex_Vectors;
with QuadDobl_Complex_Vectors_io;
with QuadDobl_Complex_VecVecs;
with Standard_Complex_Poly_Systems;
with Standard_Complex_Laur_Systems;
with DoblDobl_Complex_Poly_Systems;
with DoblDobl_Complex_Laur_Systems;
with QuadDobl_Complex_Poly_Systems;
with QuadDobl_Complex_Laur_Systems;
with Standard_Complex_Solutions;
with Standard_Complex_Solutions_io;
with Standard_Solution_Strings;
with DoblDobl_Complex_Solutions;
with DoblDobl_Complex_Solutions_io;
with DoblDobl_Solution_Strings;
with QuadDobl_Complex_Solutions;
with QuadDobl_Complex_Solutions_io;
with QuadDobl_Solution_Strings;
with Sampling_Machine;
with Sampling_Laurent_Machine;
with DoblDobl_Sampling_Machine;
with DoblDobl_Sampling_Laurent_Machine;
with QuadDobl_Sampling_Machine;
with QuadDobl_Sampling_Laurent_Machine;
with Witness_Sets;
with Homotopy_Membership_Tests;         use Homotopy_Membership_Tests;
with Multitasking_Membership_Tests;     use Multitasking_Membership_Tests;
with Assignments_in_Ada_and_C;          use Assignments_in_Ada_and_C;
with Standard_PolySys_Container;
with Standard_LaurSys_Container;
with DoblDobl_PolySys_Container;
with DoblDobl_LaurSys_Container;
with QuadDobl_PolySys_Container;
with QuadDobl_LaurSys_Container;
with Standard_Solutions_Container;
with DoblDobl_Solutions_Container;
with QuadDobl_Solutions_Container;

package body Member_Interface is

  procedure Get_Input_Parameters
              ( a : in C_intarrs.Pointer; b : in C_intarrs.Pointer;
                verbose : out boolean; nbr,dim,nbt : out integer32 ) is

  -- DESCRIPTION :
  --   Extracts the input parameters from the a and b input.

  -- ON ENTRY :
  --   a        value for the verbose parameter;
  --   b        dimensions and number of tasks.

  -- ON RETURN :
  --   verbose  flag to determine the verbosity of the test
  --            the value for verbose is in the first a value;;
  --   nbr      number of coordinates in the test point,
  --            the value for nbr is in the first b value;;
  --   dim      dimension of the witness set,
  --            the value for dim is in the second b value;;
  --   nbt      number of tasks for the homotopy membership test,
  --            the value for nbt is in the third b value..

    va : constant C_Integer_Array := C_intarrs.Value(a);
    vrb : constant integer32 := integer32(va(va'first));
    vb : constant C_Integer_Array
       := C_intarrs.Value(b,Interfaces.C.ptrdiff_t(3));

  begin
    verbose := (vrb = 1);
    nbr := integer32(vb(vb'first));
   -- put("nbr = "); put(nbr,1); new_line;
    dim := integer32(vb(vb'first+1));
   -- put("dim = "); put(dim,1); new_line;
    nbt := integer32(vb(vb'first+2));
   -- put("nbt = "); put(nbt,1); new_line;
  end Get_Input_Parameters;

  procedure Get_Input_Integers
              ( a : in C_intarrs.Pointer; verbose : out boolean;
                nbr,dim,nbc,nbt : out integer32 ) is

  -- DESCRIPTION :
  --   Extracts the input parameters from the a input,
  --   expecting five integers on input,
  --   in the first values of a.

  -- ON ENTRY :
  --   a        pointer to the integer values.

  -- ON RETURN :
  --   verbose  flag to determine the verbosity of the test;
  --   nbr      number of coordinates in the test point;
  --   dim      dimension of the witness set.
  --   nbc      number of characters in the input string;
  --   nbt      number of tasks for the homotopy membership test.

    va : constant C_Integer_Array
       := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(5));
    vrb : constant integer32 := integer32(va(va'first));

  begin
    verbose := (vrb = 1);
    nbr := integer32(va(va'first+1));
    if verbose
     then put("in use_c2mbt, nbr = "); put(nbr,1); new_line;
    end if;
    dim := integer32(va(va'first+2));
    if verbose
     then put("in use_c2mbt, dim = "); put(dim,1); new_line;
    end if;
    nbc := integer32(va(va'first+3));
    if verbose
     then put("in use_c2mbt, nbc = "); put(nbc,1); new_line;
    end if;
    nbt := integer32(va(va'first+4));
    if verbose
     then put("in use_c2mbt, nbt = "); put(nbt,1); new_line;
    end if;
  end Get_Input_Integers;

  procedure Get_Input_Tolerances
              ( c : in C_dblarrs.Pointer;
                verbose : in boolean; restol,homtol : out double_float ) is

  -- DESCRIPTION :
  --   Extracts the values in the parameter c on input
  --   for the tolerances in the homotopy membership test.

  -- ON ENTRY :
  --   c        pointer to the values with the tolerances;
  --   verbose  if verbose, then the values of the tolerance are written.

  -- ON RETURN :
  --   restol   tolerance on the residual of the evaluation;
  --   homtol   tolerance on the membership for the new generic points.

    vc : constant C_Double_Array
       := C_DblArrs.Value(c,Interfaces.C.ptrdiff_t(2));

  begin
    restol := double_float(vc(0));
    if verbose
     then put("in use_c2mbt, restol = "); put(restol,3); new_line;
    end if;
    homtol := double_float(vc(1));
    if verbose
     then put("in use_c2mbt, homtol = "); put(homtol,3); new_line;
    end if;
  end Get_Input_Tolerances;

  procedure Get_Standard_Input_Values
              ( c : in C_dblarrs.Pointer;
                nbr : in integer32; restol,homtol : out double_float;
                pt : out Standard_Complex_Vectors.Vector ) is

  -- DESCRIPTION :
  --   Extracts the values in the parameter c on input,
  --   for a test point given in standard double precision.

  -- ON ENTRY :
  --   c        pointer to the values;
  --   nbr      the dimension of the test point.
 
  -- ON RETURN :
  --   restol   tolerance on the residual of the evaluation;
  --   homtol   tolerance on the membership for the new generic points;
  --   pt       coordinates of the test point.

    nv : constant integer32 := 2 + 2*nbr;
    vc : constant C_Double_Array
       := C_DblArrs.Value(c,Interfaces.C.ptrdiff_t(nv));
    ind : Interfaces.C.size_t := 2;
    re,im : double_float;

  begin
    restol := double_float(vc(0));
    homtol := double_float(vc(1));
    for k in 1..nbr loop
      re := double_float(vc(ind));
      im := double_float(vc(ind+1));
      pt(k) := Standard_Complex_Numbers.Create(re,im);
      ind := ind + 2;
    end loop;
   -- put_line("The coordinates of the test point : ");
   -- put_line(pt);
  end Get_Standard_Input_Values;

  procedure Get_DoblDobl_Input_Values
              ( c : in C_dblarrs.Pointer;
                nbr : in integer32; restol,homtol : out double_float;
                pt : out DoblDobl_Complex_Vectors.Vector ) is

  -- DESCRIPTION :
  --   Extracts the values in the parameter c on input,
  --   for a test point given in double double precision;

  -- ON ENTRY :
  --   c        points to the values;
  --   nbr      the dimension of the test point.
 
  -- ON RETURN :
  --   restol   tolerance on the residual of the evaluation;
  --   homtol   tolerance on the membership for the new generic points;
  --   pt       coordinates of the test point.

    nv : constant integer32 := 2 + 4*nbr;
    vc : constant C_Double_Array
       := C_DblArrs.Value(c,Interfaces.C.ptrdiff_t(nv));
    ind : Interfaces.C.size_t := 2;
    rehi,relo,imhi,imlo : double_float;
    re,im : double_double;

  begin
    restol := double_float(vc(0));
    homtol := double_float(vc(1));
    for k in 1..nbr loop
      rehi := double_float(vc(ind));
      relo := double_float(vc(ind+1));
      imhi := double_float(vc(ind+2));
      imlo := double_float(vc(ind+3));
      re := create(rehi,relo);
      im := create(imhi,imlo);
      pt(k) := DoblDobl_Complex_Numbers.Create(re,im);
      ind := ind + 4;
    end loop;
  end Get_DoblDobl_Input_Values;

  procedure Get_QuadDobl_Input_Values
              ( c : in C_dblarrs.Pointer;
                nbr : in integer32; restol,homtol : out double_float;
                pt : out QuadDobl_Complex_Vectors.Vector ) is

  -- DESCRIPTION :
  --   Extracts the values in the parameter c on input,
  --   for a test point given in quad double precision.

  -- ON ENTRY :
  --   c        points to the input values;
  --   nbr      the dimension of the test point.
 
  -- ON RETURN :
  --   restol   tolerance on the residual of the evaluation;
  --   homtol   tolerance on the membership for the new generic points;
  --   pt       coordinates of the test point.

    nv : constant integer32 := 2 + 8*nbr;
    vc : constant C_Double_Array
       := C_DblArrs.Value(c,Interfaces.C.ptrdiff_t(nv));
    ind : Interfaces.C.size_t := 2;
    rehihi,rehilo,relohi,relolo : double_float;
    imhihi,imhilo,imlohi,imlolo : double_float;
    re,im : quad_double;

  begin
    restol := double_float(vc(0));
    homtol := double_float(vc(1));
    for k in 1..nbr loop
      rehihi := double_float(vc(ind));
      rehilo := double_float(vc(ind+1));
      relohi := double_float(vc(ind+2));
      relolo := double_float(vc(ind+3));
      imhihi := double_float(vc(ind+4));
      imhilo := double_float(vc(ind+5));
      imlohi := double_float(vc(ind+6));
      imlolo := double_float(vc(ind+7));
      re := create(rehihi,rehilo,relohi,relolo);
      im := create(imhihi,imhilo,imlohi,imlolo);
      pt(k) := QuadDobl_Complex_Numbers.Create(re,im);
      ind := ind + 8;
    end loop;
  end Get_QuadDobl_Input_Values;

  procedure Standard_Parse_Test_Point
              ( b : in C_intarrs.Pointer;
                verbose : in boolean; nbr,nbc : in integer32;
                pt : out Standard_Complex_Vectors.Vector ) is

  -- DESCRIPTION :
  --   Extracts the values in the parameter b on input,
  --   for a test point given in standard double precision.

  -- ON ENTRY :
  --   b        values for the string;
  --   verbose  if verbose, then debugging info is written;
  --   nbr      the dimension of the test point;
  --   nbc      the number of characters in the string representation
  --            of the solution which containes the test point.
 
  -- ON RETURN :
  --   pt       coordinates of the test point.

    vb : constant C_Integer_Array
       := C_intarrs.Value(b,Interfaces.C.ptrdiff_t(nbc));
    str : constant string := C_Integer_Array_to_String(natural32(nbc),vb);
    sol : Standard_Complex_Solutions.Solution(nbr);
    idx : integer := str'first;
    fail : boolean;

  begin
    if verbose then
      put_line("in use_c2mbt, the string for the test point :");
      put_line(str);
    end if;
    Standard_Solution_Strings.Parse(str,idx,natural32(nbr),sol,fail);
    if verbose then
      put_line("in use_c2mbt, the parsed solution for the test point :");
      Standard_Complex_Solutions_io.put(sol); new_line;
    end if;
    if not fail then
      for k in sol.v'range loop
        pt(k) := sol.v(k);
      end loop;
    end if;
    if verbose then
      put_line("in use_c2mbt, the coordinates of the test point : ");
      Standard_Complex_Vectors_io.put_line(pt);
    end if;
  end Standard_Parse_Test_Point;

  procedure DoblDobl_Parse_Test_Point
              ( b : in C_intarrs.Pointer;
                verbose : in boolean; nbr,nbc : in integer32;
                pt : out DoblDobl_Complex_Vectors.Vector ) is

  -- DESCRIPTION :
  --   Extracts the values in the parameter b on input,
  --   for a test point given in double double precision.

  -- ON ENTRY :
  --   b        values for the string;
  --   verbose  if verbose, then debugging info is written;
  --   nbr      the dimension of the test point;
  --   nbc      the number of characters in the string representation
  --            of the solution which containes the test point.
 
  -- ON RETURN :
  --   pt       coordinates of the test point.

    vb : constant C_Integer_Array
       := C_intarrs.Value(b,Interfaces.C.ptrdiff_t(nbc));
    str : constant string := C_Integer_Array_to_String(natural32(nbc),vb);
    sol : DoblDobl_Complex_Solutions.Solution(nbr);
    idx : integer := str'first;
    fail : boolean;

  begin
    if verbose then
      put_line("in use_c2mbt, the string for the test point :");
      put_line(str);
    end if;
    DoblDobl_Solution_Strings.Parse(str,idx,natural32(nbr),sol,fail);
    if verbose then
      put_line("in use_c2mbt, the parsed solution for the test point :");
      DoblDobl_Complex_Solutions_io.put(sol); new_line;
    end if;
    if not fail then
      for k in sol.v'range loop
        pt(k) := sol.v(k);
      end loop;
    end if;
    if verbose then
      put_line("in use_c2mbt, the coordinates of the test point : ");
      DoblDobl_Complex_Vectors_io.put_line(pt);
    end if;
  end DoblDobl_Parse_Test_Point;

  procedure QuadDobl_Parse_Test_Point
              ( b : in C_intarrs.Pointer;
                verbose : in boolean; nbr,nbc : in integer32;
                pt : out QuadDobl_Complex_Vectors.Vector ) is

  -- DESCRIPTION :
  --   Extracts the values in the parameter b on input,
  --   for a test point given in double double precision.

  -- ON ENTRY :
  --   b        values for the string;
  --   verbose  if verbose, then debugging info is written;
  --   nbr      the dimension of the test point;
  --   nbc      the number of characters in the string representation
  --            of the solution which containes the test point.
 
  -- ON RETURN :
  --   pt       coordinates of the test point.

    vb : constant C_Integer_Array
       := C_intarrs.Value(b,Interfaces.C.ptrdiff_t(nbc));
    str : constant string := C_Integer_Array_to_String(natural32(nbc),vb);
    sol : QuadDobl_Complex_Solutions.Solution(nbr);
    idx : integer := str'first;
    fail : boolean;

  begin
    if verbose then
      put_line("in use_c2mbt, the string for the test point :");
      put_line(str);
    end if;
    QuadDobl_Solution_Strings.Parse(str,idx,natural32(nbr),sol,fail);
    if verbose then
      put_line("in use_c2mbt, the parsed solution for the test point :");
      QuadDobl_Complex_Solutions_io.put(sol); new_line;
    end if;
    if not fail then
      for k in sol.v'range loop
        pt(k) := sol.v(k);
      end loop;
    end if;
    if verbose then
      put_line("in use_c2mbt, the coordinates of the test point : ");
      QuadDobl_Complex_Vectors_io.put_line(pt);
    end if;
  end QuadDobl_Parse_Test_Point;

  procedure Assign_Results
              ( a : in C_intarrs.Pointer; b : in C_intarrs.Pointer;
                onpolsys,inwitset : in boolean ) is

  -- DESCRIPTION :
  --   Assigns the results of the test to a and b.

  -- ON ENTRY :
  --   onpolsys is true if the test point satisfied the polynomials,
  --            and is false otherwise;
  --   inwitset is true if the test point belongs to the witness set,
  --            and is false otherwise.

  -- ON RETURN :
  --   a        0 for false and 1 for true, if satisfies the equations;
  --   b        0 for false and 1 for true, if in the witness set.

  begin
    if not onpolsys then
      Assign(0,a);
      Assign(0,b);
    else
      Assign(1,a);
      if inwitset
       then Assign(1,b);
       else Assign(0,b);
      end if;
    end if;
  end Assign_Results;

  function Member_Standard_Polynomial_Numbers
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

  -- DESCRIPTION :
  --   Runs a membership test on polynomial systems in double precision.
  --   The test point is given as a sequence of numbers.

    use Standard_Complex_Poly_Systems;
    use Standard_Complex_Solutions;

    verbose,onp,inw : boolean;
    nbr,dim,nbtasks : integer32;
    restol,homtol : double_float;
    lp : constant Link_to_Poly_Sys := Standard_PolySys_Container.Retrieve;
    sols : constant Solution_List := Standard_Solutions_Container.Retrieve;

  begin
    if vrblvl > 0 then
      put("-> in member_interface.");
      put_line("Member_Standard_Polynomial_Numbers ...");
    end if;
    Sampling_Machine.Initialize(lp.all);
    Sampling_Machine.Default_Tune_Sampler(2);
    Sampling_Machine.Default_Tune_Refiner;
    Get_Input_Parameters(a,b,verbose,nbr,dim,nbtasks);
    declare
      tpt : Standard_Complex_Vectors.Vector(1..nbr);
      sli : Standard_Complex_VecVecs.VecVec(1..dim)
          := Witness_Sets.Slices(lp.all,natural32(dim));
    begin
      Get_Standard_Input_Values(c,nbr,restol,homtol,tpt);
      if nbtasks = 0 then
        Homotopy_Membership_Test
          (verbose,lp.all,natural32(dim),sli,sols,tpt,restol,homtol,onp,inw);
      else
        Homotopy_Membership_Test
          (verbose,natural32(nbtasks),lp.all,natural32(dim),sols,tpt,
           restol,homtol,onp,inw);
      end if;
      Standard_Complex_VecVecs.Clear(sli);
    end;
    Assign_Results(a,b,onp,inw);
    Sampling_Machine.Clear;
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in member_interface.");
        put_line("Member_Standard_Polynomial_Numbers");
      end if;
      return 537;
  end Member_Standard_Polynomial_Numbers;

  function Member_DoblDobl_Polynomial_Numbers
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

  -- DESCRIPTION :
  --   Runs a membership test on polynomial systems in double double precision.
  --   The test point is given as a sequence of numbers.

    use DoblDobl_Complex_Poly_Systems;
    use DoblDobl_Complex_Solutions;

    verbose,onp,inw : boolean;
    nbr,dim,nbtasks : integer32;
    restol,homtol : double_float;
    lp : constant Link_to_Poly_Sys := DoblDobl_PolySys_Container.Retrieve;
    sols : constant Solution_List := DoblDobl_Solutions_Container.Retrieve;

  begin
    if vrblvl > 0 then
      put("-> in member_interface.");
      put_line("Member_DoblDobl_Polynomial_Numbers ...");
    end if;
    DoblDobl_Sampling_Machine.Initialize(lp.all);
    DoblDobl_Sampling_Machine.Default_Tune_Sampler(0);
    DoblDobl_Sampling_Machine.Default_Tune_Refiner;
    Get_Input_Parameters(a,b,verbose,nbr,dim,nbtasks);
    declare
      tpt : DoblDobl_Complex_Vectors.Vector(1..nbr);
      sli : DoblDobl_Complex_VecVecs.VecVec(1..dim)
          := Witness_Sets.Slices(lp.all,natural32(dim));
    begin
      Get_DoblDobl_Input_Values(c,nbr,restol,homtol,tpt);
      if nbtasks = 0 then
        Homotopy_Membership_Test
          (verbose,lp.all,natural32(dim),sli,sols,tpt,restol,homtol,onp,inw);
      else
        Homotopy_Membership_Test
          (verbose,natural32(nbtasks),lp.all,natural32(dim),sols,tpt,
           restol,homtol,onp,inw);
      end if;
      DoblDobl_Complex_VecVecs.Clear(sli);
    end;
    Assign_Results(a,b,onp,inw);
    DoblDobl_Sampling_Machine.Clear;
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in member_interface.");
        put_line("Member_DoblDobl_Polynomial_Numbers");
      end if;
      return 538;
  end Member_DoblDobl_Polynomial_Numbers;

  function Member_QuadDobl_Polynomial_Numbers
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

  -- DESCRIPTION :
  --   Runs a membership test on polynomial systems in quad double precision.
  --   The test point is given as a sequence of numbers.

    use QuadDobl_Complex_Poly_Systems;
    use QuadDobl_Complex_Solutions;

    verbose,onp,inw : boolean;
    nbr,dim,nbtasks : integer32;
    restol,homtol : double_float;
    lp : constant Link_to_Poly_Sys := QuadDobl_PolySys_Container.Retrieve;
    sols : constant Solution_List := QuadDobl_Solutions_Container.Retrieve;

  begin
    if vrblvl > 0 then
      put("-> in member_interface.");
      put_line("Member_QuadDobl_Polynomial_Numbers ...");
    end if;
    QuadDobl_Sampling_Machine.Initialize(lp.all);
    QuadDobl_Sampling_Machine.Default_Tune_Sampler(0);
    QuadDobl_Sampling_Machine.Default_Tune_Refiner;
    Get_Input_Parameters(a,b,verbose,nbr,dim,nbtasks);
    declare
      tpt : QuadDobl_Complex_Vectors.Vector(1..nbr);
      sli : QuadDobl_Complex_VecVecs.VecVec(1..dim)
          := Witness_Sets.Slices(lp.all,natural32(dim));
    begin
      Get_QuadDobl_Input_Values(c,nbr,restol,homtol,tpt);
      if nbtasks = 0 then
        Homotopy_Membership_Test
          (verbose,lp.all,natural32(dim),sli,sols,tpt,restol,homtol,onp,inw);
      else
        Homotopy_Membership_Test
          (verbose,natural32(nbtasks),lp.all,natural32(dim),sols,tpt,
           restol,homtol,onp,inw);
      end if;
      QuadDobl_Complex_VecVecs.Clear(sli);
    end;
    Assign_Results(a,b,onp,inw);
    QuadDobl_Sampling_Machine.Clear;
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in member_interface.");
        put_line("Member_QuadDobl_Polynomial_Numbers");
      end if;
      return 539;
  end Member_QuadDobl_Polynomial_Numbers;

  function Member_Standard_Laurent_Numbers
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

  -- DESCRIPTION :
  --   Runs a membership test on Laurent systems in double precision.
  --   The test point is given as sequence of numbers.

    use Standard_Complex_Laur_Systems;
    use Standard_Complex_Solutions;

    verbose,onp,inw : boolean;
    nbr,dim,nbtasks : integer32;
    restol,homtol : double_float;
    lp : constant Link_to_Laur_Sys := Standard_LaurSys_Container.Retrieve;
    sols : constant Solution_List := Standard_Solutions_Container.Retrieve;

  begin
    if vrblvl > 0 then
      put("-> in member_interface.");
      put_line("Member_Standard_Laurent_Numbers ...");
    end if;
    Sampling_Laurent_Machine.Initialize(lp.all);
    Sampling_Laurent_Machine.Default_Tune_Sampler(2);
    Sampling_Laurent_Machine.Default_Tune_Refiner;
    Get_Input_Parameters(a,b,verbose,nbr,dim,nbtasks);
    declare
      tpt : Standard_Complex_Vectors.Vector(1..nbr);
      sli : Standard_Complex_VecVecs.VecVec(1..dim)
          := Witness_Sets.Slices(lp.all,natural32(dim));
    begin
      Get_Standard_Input_Values(c,nbr,restol,homtol,tpt);
      if nbtasks = 0 then
        Homotopy_Membership_Test
          (verbose,lp.all,natural32(dim),sli,sols,tpt,restol,homtol,onp,inw);
      else
        Homotopy_Membership_Test
          (verbose,natural32(nbtasks),lp.all,natural32(dim),sols,tpt,
           restol,homtol,onp,inw);
      end if;
      Standard_Complex_VecVecs.Clear(sli);
    end;
    Assign_Results(a,b,onp,inw);
    Sampling_Laurent_Machine.Clear;
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in member_interface.");
        put_line("Member_Standard_Laurent_Numbers");
      end if;
      return 795;
  end Member_Standard_Laurent_Numbers;

  function Member_DoblDobl_Laurent_Numbers
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

  -- DESCRIPTION :
  --   Runs a membership test on Laurent systems in double double precision.
  --   The test point is given as a sequence of numbers.

    use DoblDobl_Complex_Laur_Systems;
    use DoblDobl_Complex_Solutions;

    verbose,onp,inw : boolean;
    nbr,dim,nbtasks : integer32;
    restol,homtol : double_float;
    lp : constant Link_to_Laur_Sys := DoblDobl_LaurSys_Container.Retrieve;
    sols : constant Solution_List := DoblDobl_Solutions_Container.Retrieve;

  begin
    if vrblvl > 0 then
      put("-> in member_interface.");
      put_line("Member_DoblDobl_Laurent_Numbers ...");
    end if;
    DoblDobl_Sampling_Laurent_Machine.Initialize(lp.all);
    DoblDobl_Sampling_Laurent_Machine.Default_Tune_Sampler(0);
    DoblDobl_Sampling_Laurent_Machine.Default_Tune_Refiner;
    Get_Input_Parameters(a,b,verbose,nbr,dim,nbtasks);
    declare
      tpt : DoblDobl_Complex_Vectors.Vector(1..nbr);
      sli : DoblDobl_Complex_VecVecs.VecVec(1..dim)
          := Witness_Sets.Slices(lp.all,natural32(dim));
    begin
      Get_DoblDobl_Input_Values(c,nbr,restol,homtol,tpt);
      if nbtasks = 0 then
        Homotopy_Membership_Test
          (verbose,lp.all,natural32(dim),sli,sols,tpt,restol,homtol,onp,inw);
      else
        Homotopy_Membership_Test
          (verbose,natural32(nbtasks),lp.all,natural32(dim),sols,tpt,
           restol,homtol,onp,inw);
      end if;
      DoblDobl_Complex_VecVecs.Clear(sli);
    end;
    Assign_Results(a,b,onp,inw);
    DoblDobl_Sampling_Laurent_Machine.Clear;
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in member_interface.");
        put_line("Member_DoblDobl_Laurent_Numbers");
      end if;
      return 796;
  end Member_DoblDobl_Laurent_Numbers;

  function Member_QuadDobl_Laurent_Numbers
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

  -- DESCRIPTION :
  --   Runs a membership test on Laurent systems in quad double precision.
  --   The test point is given as a sequence of numbers.

    use QuadDobl_Complex_Laur_Systems;
    use QuadDobl_Complex_Solutions;

    verbose,onp,inw : boolean;
    nbr,dim,nbtasks : integer32;
    restol,homtol : double_float;
    lp : constant Link_to_Laur_Sys := QuadDobl_LaurSys_Container.Retrieve;
    sols : constant Solution_List := QuadDobl_Solutions_Container.Retrieve;

  begin
    if vrblvl > 0 then
      put("-> in member_interface.");
      put_line("Member_QuadDobl_Laurent_Numbers ...");
    end if;
    QuadDobl_Sampling_Laurent_Machine.Initialize(lp.all);
    QuadDobl_Sampling_Laurent_Machine.Default_Tune_Sampler(0);
    QuadDobl_Sampling_Laurent_Machine.Default_Tune_Refiner;
    Get_Input_Parameters(a,b,verbose,nbr,dim,nbtasks);
    declare
      tpt : QuadDobl_Complex_Vectors.Vector(1..nbr);
      sli : QuadDobl_Complex_VecVecs.VecVec(1..dim)
          := Witness_Sets.Slices(lp.all,natural32(dim));
    begin
      Get_QuadDobl_Input_Values(c,nbr,restol,homtol,tpt);
      if nbtasks = 0 then
        Homotopy_Membership_Test
          (verbose,lp.all,natural32(dim),sli,sols,tpt,restol,homtol,onp,inw);
      else
        Homotopy_Membership_Test
          (verbose,natural32(nbtasks),lp.all,natural32(dim),sols,tpt,
           restol,homtol,onp,inw);
      end if;
      QuadDobl_Complex_VecVecs.Clear(sli);
    end;
    Assign_Results(a,b,onp,inw);
    QuadDobl_Sampling_Laurent_Machine.Clear;
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in member_interface.");
        put_line("Member_QuadDobl_Laurent_Numbers");
      end if;
      return 797;
  end Member_QuadDobl_Laurent_Numbers;

  function Member_Standard_Polynomial_String
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

  -- DESCRIPTION :
  --   Runs a membership test on polynomials in double precision.
  --   The test point is given as a string.

    use Standard_Complex_Poly_Systems;
    use Standard_Complex_Solutions;

    verbose,onp,inw : boolean;
    nbr,dim,nbc,nbtasks : integer32;
    restol,homtol : double_float;
    lp : constant Link_to_Poly_Sys := Standard_PolySys_Container.Retrieve;
    sols : constant Solution_List := Standard_Solutions_Container.Retrieve;

  begin
    if vrblvl > 0 then
      put("-> in member_interface.");
      put_line("Member_Standard_Laurent_String ...");
    end if;
    Sampling_Machine.Initialize(lp.all);
    Sampling_Machine.Default_Tune_Sampler(2);
    Sampling_Machine.Default_Tune_Refiner;
    Get_Input_Integers(a,verbose,nbr,dim,nbc,nbtasks);
    Get_Input_Tolerances(c,verbose,restol,homtol);
    declare
      tpt : Standard_Complex_Vectors.Vector(1..nbr);
      sli : Standard_Complex_VecVecs.VecVec(1..dim)
          := Witness_Sets.Slices(lp.all,natural32(dim));
    begin
      Standard_Parse_Test_Point(b,verbose,nbr,nbc,tpt);
      if nbtasks = 0 then
        Homotopy_Membership_Test
          (verbose,lp.all,natural32(dim),sli,sols,tpt,restol,homtol,onp,inw);
      else
        Homotopy_Membership_Test
          (verbose,natural32(nbtasks),lp.all,natural32(dim),sols,tpt,
           restol,homtol,onp,inw);
      end if;
      Standard_Complex_VecVecs.Clear(sli);
    end;
    Assign_Results(a,b,onp,inw);
    Sampling_Machine.Clear;
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in member_interface.");
        put_line("Member_Standard_Polynomial_String");
      end if;
      return 822;
  end Member_Standard_Polynomial_String;

  function Member_DoblDobl_Polynomial_String
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

  -- DESCRIPTION :
  --   Runs a membership test on polynomial systems in double double precision.
  --   The test point is given as a string.

    use DoblDobl_Complex_Poly_Systems;
    use DoblDobl_Complex_Solutions;

    verbose,onp,inw : boolean;
    nbr,dim,nbc,nbtasks : integer32;
    restol,homtol : double_float;
    lp : constant Link_to_Poly_Sys := DoblDobl_PolySys_Container.Retrieve;
    sols : constant Solution_List := DoblDobl_Solutions_Container.Retrieve;

  begin
    if vrblvl > 0 then
      put("-> in member_interface.");
      put_line("Member_DoblDobl_Laurent_String ...");
    end if;
    DoblDobl_Sampling_Machine.Initialize(lp.all);
    DoblDobl_Sampling_Machine.Default_Tune_Sampler(2);
    DoblDobl_Sampling_Machine.Default_Tune_Refiner;
    Get_Input_Integers(a,verbose,nbr,dim,nbc,nbtasks);
    Get_Input_Tolerances(c,verbose,restol,homtol);
    declare
      tpt : DoblDobl_Complex_Vectors.Vector(1..nbr);
      sli : DoblDobl_Complex_VecVecs.VecVec(1..dim)
          := Witness_Sets.Slices(lp.all,natural32(dim));
    begin
      DoblDobl_Parse_Test_Point(b,verbose,nbr,nbc,tpt);
      if nbtasks = 0 then
        Homotopy_Membership_Test
          (verbose,lp.all,natural32(dim),sli,sols,tpt,restol,homtol,onp,inw);
      else
        Homotopy_Membership_Test
          (verbose,natural32(nbtasks),lp.all,natural32(dim),sols,tpt,
           restol,homtol,onp,inw);
      end if;
      DoblDobl_Complex_VecVecs.Clear(sli);
    end;
    Assign_Results(a,b,onp,inw);
    DoblDobl_Sampling_Machine.Clear;
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in member_interface.");
        put_line("Member_DoblDobl_Polynomial_String");
      end if;
      return 823;
  end Member_DoblDobl_Polynomial_String;

  function Member_QuadDobl_Polynomial_String
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

  -- DESCRIPTION :
  --   Runs a membership test on polynomial systems in quad double precision.
  --   The test point is given as a string.

    use QuadDobl_Complex_Poly_Systems;
    use QuadDobl_Complex_Solutions;

    verbose,onp,inw : boolean;
    nbr,dim,nbc,nbtasks : integer32;
    restol,homtol : double_float;
    lp : constant Link_to_Poly_Sys := QuadDobl_PolySys_Container.Retrieve;
    sols : constant Solution_List := QuadDobl_Solutions_Container.Retrieve;

  begin
    if vrblvl > 0 then
      put("-> in member_interface.");
      put_line("Member_QuadDobl_Polynomial_String ...");
    end if;
    QuadDobl_Sampling_Machine.Initialize(lp.all);
    QuadDobl_Sampling_Machine.Default_Tune_Sampler(2);
    QuadDobl_Sampling_Machine.Default_Tune_Refiner;
    Get_Input_Integers(a,verbose,nbr,dim,nbc,nbtasks);
    Get_Input_Tolerances(c,verbose,restol,homtol);
    declare
      tpt : QuadDobl_Complex_Vectors.Vector(1..nbr);
      sli : QuadDobl_Complex_VecVecs.VecVec(1..dim)
          := Witness_Sets.Slices(lp.all,natural32(dim));
    begin
      QuadDobl_Parse_Test_Point(b,verbose,nbr,nbc,tpt);
      if nbtasks = 0 then
        Homotopy_Membership_Test
          (verbose,lp.all,natural32(dim),sli,sols,tpt,restol,homtol,onp,inw);
      else
        Homotopy_Membership_Test
          (verbose,natural32(nbtasks),lp.all,natural32(dim),sols,tpt,
           restol,homtol,onp,inw);
      end if;
      QuadDobl_Complex_VecVecs.Clear(sli);
    end;
    Assign_Results(a,b,onp,inw);
    QuadDobl_Sampling_Machine.Clear;
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in member_interface.");
        put_line("Member_QuadDobl_Polynomial_String");
      end if;
      return 824;
  end Member_QuadDobl_Polynomial_String;

  function Member_Standard_Laurent_String
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

  -- DESCRIPTION :
  --   Runs a membership test on Laurent systems in double precision.
  --   The test point is given as a string.

    use Standard_Complex_Laur_Systems;
    use Standard_Complex_Solutions;

    verbose,onp,inw : boolean;
    nbr,dim,nbc,nbtasks : integer32;
    restol,homtol : double_float;
    lp : constant Link_to_Laur_Sys := Standard_LaurSys_Container.Retrieve;
    sols : constant Solution_List := Standard_Solutions_Container.Retrieve;

  begin
    if vrblvl > 0 then
      put("-> in member_interface.");
      put_line("Member_Standard_Laurent_String ...");
    end if;
    Sampling_Laurent_Machine.Initialize(lp.all);
    Sampling_Laurent_Machine.Default_Tune_Sampler(2);
    Sampling_Laurent_Machine.Default_Tune_Refiner;
    Get_Input_Integers(a,verbose,nbr,dim,nbc,nbtasks);
    Get_Input_Tolerances(c,verbose,restol,homtol);
    declare
      tpt : Standard_Complex_Vectors.Vector(1..nbr);
      sli : Standard_Complex_VecVecs.VecVec(1..dim)
          := Witness_Sets.Slices(lp.all,natural32(dim));
    begin
      Standard_Parse_Test_Point(b,verbose,nbr,nbc,tpt);
      if nbtasks = 0 then
        Homotopy_Membership_Test
          (verbose,lp.all,natural32(dim),sli,sols,tpt,restol,homtol,onp,inw);
      else
        Homotopy_Membership_Test
          (verbose,natural32(nbtasks),lp.all,natural32(dim),sols,tpt,
           restol,homtol,onp,inw);
      end if;
      Standard_Complex_VecVecs.Clear(sli);
    end;
    Assign_Results(a,b,onp,inw);
    Sampling_Laurent_Machine.Clear;
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in member_interface.");
        put_line("Member_Standard_Laurent_String");
      end if;
      return 825;
  end Member_Standard_Laurent_String;

  function Member_DoblDobl_Laurent_String
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

  -- DESCRIPTION :
  --   Runs a membership test on Laurent systems in double double precision.
  --   The test point is given as a string.

    use DoblDobl_Complex_Laur_Systems;
    use DoblDobl_Complex_Solutions;

    verbose,onp,inw : boolean;
    nbr,dim,nbc,nbtasks : integer32;
    restol,homtol : double_float;
    lp : constant Link_to_Laur_Sys := DoblDobl_LaurSys_Container.Retrieve;
    sols : constant Solution_List := DoblDobl_Solutions_Container.Retrieve;

  begin
    if vrblvl > 0 then
      put("-> in member_interface.");
      put_line("Member_DoblDobl_Laurent_String ...");
    end if;
    DoblDobl_Sampling_Laurent_Machine.Initialize(lp.all);
    DoblDobl_Sampling_Laurent_Machine.Default_Tune_Sampler(2);
    DoblDobl_Sampling_Laurent_Machine.Default_Tune_Refiner;
    Get_Input_Integers(a,verbose,nbr,dim,nbc,nbtasks);
    Get_Input_Tolerances(c,verbose,restol,homtol);
    declare
      tpt : DoblDobl_Complex_Vectors.Vector(1..nbr);
      sli : DoblDobl_Complex_VecVecs.VecVec(1..dim)
          := Witness_Sets.Slices(lp.all,natural32(dim));
    begin
      DoblDobl_Parse_Test_Point(b,verbose,nbr,nbc,tpt);
      if nbtasks = 0 then
        Homotopy_Membership_Test
          (verbose,lp.all,natural32(dim),sli,sols,tpt,restol,homtol,onp,inw);
      else
        Homotopy_Membership_Test
          (verbose,natural32(nbtasks),lp.all,natural32(dim),sols,tpt,
           restol,homtol,onp,inw);
      end if;
      DoblDobl_Complex_VecVecs.Clear(sli);
    end;
    Assign_Results(a,b,onp,inw);
    DoblDobl_Sampling_Laurent_Machine.Clear;
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in member_interface.");
        put_line("Member_DoblDobl_Laurent_String");
      end if;
      return 826;
  end Member_DoblDobl_Laurent_String;

  function Member_QuadDobl_Laurent_String
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

  -- DESCRIPTION :
  --   Runs a membership test on Laurent systems in quad double precision.
  --   The test point is given as a string.

    use QuadDobl_Complex_Laur_Systems;
    use QuadDobl_Complex_Solutions;

    verbose,onp,inw : boolean;
    nbr,dim,nbc,nbtasks : integer32;
    restol,homtol : double_float;
    lp : constant Link_to_Laur_Sys := QuadDobl_LaurSys_Container.Retrieve;
    sols : constant Solution_List := QuadDobl_Solutions_Container.Retrieve;

  begin
    if vrblvl > 0 then
      put("-> in member_interface.");
      put_line("Member_QuadDobl_Laurent_String ...");
    end if;
    QuadDobl_Sampling_Laurent_Machine.Initialize(lp.all);
    QuadDobl_Sampling_Laurent_Machine.Default_Tune_Sampler(2);
    QuadDobl_Sampling_Laurent_Machine.Default_Tune_Refiner;
    Get_Input_Integers(a,verbose,nbr,dim,nbc,nbtasks);
    Get_Input_Tolerances(c,verbose,restol,homtol);
    declare
      tpt : QuadDobl_Complex_Vectors.Vector(1..nbr);
      sli : QuadDobl_Complex_VecVecs.VecVec(1..dim)
          := Witness_Sets.Slices(lp.all,natural32(dim));
    begin
      QuadDobl_Parse_Test_Point(b,verbose,nbr,nbc,tpt);
      if nbtasks = 0 then
        Homotopy_Membership_Test
          (verbose,lp.all,natural32(dim),sli,sols,tpt,restol,homtol,onp,inw);
      else
        Homotopy_Membership_Test
          (verbose,natural32(nbtasks),lp.all,natural32(dim),sols,tpt,
           restol,homtol,onp,inw);
      end if;
      QuadDobl_Complex_VecVecs.Clear(sli);
    end;
    Assign_Results(a,b,onp,inw);
    QuadDobl_Sampling_Laurent_Machine.Clear;
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in member_interface.");
        put_line("Member_QuadDobl_Laurent_String");
      end if;
      return 827;
  end Member_QuadDobl_Laurent_String;

end Member_Interface;
