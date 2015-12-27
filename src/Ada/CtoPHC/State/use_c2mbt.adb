with text_io;                           use text_io;
with Interfaces.C;                      use Interfaces.C;
with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
--with Standard_Integer_Numbers_io;
-- use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;         use Standard_Floating_Numbers;
with Standard_Complex_Numbers;
with Double_Double_Numbers;             use Double_Double_Numbers;
with Quad_Double_Numbers;               use Quad_Double_Numbers;
with DoblDobl_Complex_Numbers;
with QuadDobl_Complex_Numbers;
with Standard_Complex_Vectors;
--with Standard_Complex_Vectors_io;
-- use Standard_Complex_Vectors_io;
with Standard_Complex_VecVecs;
with DoblDobl_Complex_Vectors;
with DoblDobl_Complex_VecVecs;
with QuadDobl_Complex_Vectors;
with QuadDobl_Complex_VecVecs;
with Standard_Complex_Poly_Systems;
with DoblDobl_Complex_Poly_Systems;
with QuadDobl_Complex_Poly_Systems;
with Standard_Complex_Solutions;
with DoblDobl_Complex_Solutions;
with QuadDobl_Complex_Solutions;
with Sampling_Machine;
with DoblDobl_Sampling_Machine;
with QuadDobl_Sampling_Machine;
with Witness_Sets;
with Homotopy_Membership_Tests;         use Homotopy_Membership_Tests;
with Assignments_in_Ada_and_C;          use Assignments_in_Ada_and_C;
with Standard_PolySys_Container;
with DoblDobl_PolySys_Container;
with QuadDobl_PolySys_Container;
with Standard_Solutions_Container;
with DoblDobl_Solutions_Container;
with QuadDobl_Solutions_Container;

function use_c2mbt ( job : integer32;
                     a : C_intarrs.Pointer;
		     b : C_intarrs.Pointer;
                     c : C_dblarrs.Pointer ) return integer32 is

  procedure Get_Input_Parameters
              ( verbose : out boolean; nbr,dim : out integer32 ) is

  -- DESCRIPTION :
  --   Extracts the input parameters from the a and b input.

  -- ON RETURN :
  --   verbose  flag to determine the verbosity of the test;
  --   nbr      number of coordinates in the test point;
  --   dim      dimension of the witness set.

    va : constant C_Integer_Array := C_intarrs.Value(a);
    vrb : constant integer32 := integer32(va(va'first));
    vb : constant C_Integer_Array
       := C_intarrs.Value(b,Interfaces.C.ptrdiff_t(2));

  begin
    verbose := (vrb = 1);
    nbr := integer32(vb(vb'first));
   -- put("nbr = "); put(nbr,1); new_line;
    dim := integer32(vb(vb'first+1));
   -- put("dim = "); put(dim,1); new_line;
  end Get_Input_Parameters;

  procedure Get_Standard_Input_Values
              ( nbr : in integer32; restol,homtol : out double_float;
                pt : out Standard_Complex_Vectors.Vector ) is

  -- DESCRIPTION :
  --   Extracts the values in the parameter c on input,
  --   for a test point given in standard double precision.

  -- ON ENTRY :
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
              ( nbr : in integer32; restol,homtol : out double_float;
                pt : out DoblDobl_Complex_Vectors.Vector ) is

  -- DESCRIPTION :
  --   Extracts the values in the parameter c on input,
  --   for a test point given in double double precision;

  -- ON ENTRY :
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
              ( nbr : in integer32; restol,homtol : out double_float;
                pt : out QuadDobl_Complex_Vectors.Vector ) is

  -- DESCRIPTION :
  --   Extracts the values in the parameter c on input,
  --   for a test point given in quad double precision.

  -- ON ENTRY :
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

  procedure Assign_Results
              ( onpolsys,inwitset : in boolean ) is

  -- DESCRIPTION :
  --   Assigns the results of the test to a and b.

  -- ON ENTRY :
  --   onpolsys is true if the test point satisfied the polynomials,
  --            and is false otherwise;
  --   inwitset is true if the test point belongs to the witness set,
  --            and is false otherwise.

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

  function Job0 return integer32 is -- standard double membership test

    use Standard_Complex_Poly_Systems;
    use Standard_Complex_Solutions;

    verbose,onp,inw : boolean;
    nbr,dim : integer32;
    restol,homtol : double_float;
    lp : constant Link_to_Poly_Sys := Standard_PolySys_Container.Retrieve;
    sols : constant Solution_List := Standard_Solutions_Container.Retrieve;

  begin
    Sampling_Machine.Initialize(lp.all);
    Sampling_Machine.Default_Tune_Sampler(2);
    Sampling_Machine.Default_Tune_Refiner;
    Get_Input_Parameters(verbose,nbr,dim);
    declare
      tpt : Standard_Complex_Vectors.Vector(1..nbr);
      sli : Standard_Complex_VecVecs.VecVec(1..dim)
          := Witness_Sets.Slices(lp.all,natural32(dim));
    begin
      Get_Standard_Input_Values(nbr,restol,homtol,tpt);
      Homotopy_Membership_Test
        (verbose,lp.all,natural32(dim),sli,sols,tpt,restol,homtol,onp,inw);
      Standard_Complex_VecVecs.Clear(sli);
    end;
    Assign_Results(onp,inw);
    Sampling_Machine.Clear;
    return 0;
  end Job0;

  function Job1 return integer32 is -- double double membership test

    use DoblDobl_Complex_Poly_Systems;
    use DoblDobl_Complex_Solutions;

    verbose,onp,inw : boolean;
    nbr,dim : integer32;
    restol,homtol : double_float;
    lp : constant Link_to_Poly_Sys := DoblDobl_PolySys_Container.Retrieve;
    sols : constant Solution_List := DoblDobl_Solutions_Container.Retrieve;

  begin
    DoblDobl_Sampling_Machine.Initialize(lp.all);
    DoblDobl_Sampling_Machine.Default_Tune_Sampler(0);
    DoblDobl_Sampling_Machine.Default_Tune_Refiner;
    Get_Input_Parameters(verbose,nbr,dim);
    declare
      tpt : DoblDobl_Complex_Vectors.Vector(1..nbr);
      sli : DoblDobl_Complex_VecVecs.VecVec(1..dim)
          := Witness_Sets.Slices(lp.all,natural32(dim));
    begin
      Get_DoblDobl_Input_Values(nbr,restol,homtol,tpt);
      Homotopy_Membership_Test
        (verbose,lp.all,natural32(dim),sli,sols,tpt,restol,homtol,onp,inw);
      DoblDobl_Complex_VecVecs.Clear(sli);
    end;
    Assign_Results(onp,inw);
    DoblDobl_Sampling_Machine.Clear;
    return 0;
  end Job1;

  function Job2 return integer32 is -- quad double membership test

    use QuadDobl_Complex_Poly_Systems;
    use QuadDobl_Complex_Solutions;

    verbose,onp,inw : boolean;
    nbr,dim : integer32;
    restol,homtol : double_float;
    lp : constant Link_to_Poly_Sys := QuadDobl_PolySys_Container.Retrieve;
    sols : constant Solution_List := QuadDobl_Solutions_Container.Retrieve;

  begin
    QuadDobl_Sampling_Machine.Initialize(lp.all);
    QuadDobl_Sampling_Machine.Default_Tune_Sampler(0);
    QuadDobl_Sampling_Machine.Default_Tune_Refiner;
    Get_Input_Parameters(verbose,nbr,dim);
    declare
      tpt : QuadDobl_Complex_Vectors.Vector(1..nbr);
      sli : QuadDobl_Complex_VecVecs.VecVec(1..dim)
          := Witness_Sets.Slices(lp.all,natural32(dim));
    begin
      Get_QuadDobl_Input_Values(nbr,restol,homtol,tpt);
      Homotopy_Membership_Test
        (verbose,lp.all,natural32(dim),sli,sols,tpt,restol,homtol,onp,inw);
      QuadDobl_Complex_VecVecs.Clear(sli);
    end;
    Assign_Results(onp,inw);
    QuadDobl_Sampling_Machine.Clear;
    return 0;
  end Job2;

  function Handle_Jobs return integer32 is
  begin
    case job is
      when 0 => return Job0; -- run membership test with standard doubles
      when 1 => return Job1; -- run membership test with double doubles
      when 2 => return Job2; -- run membership test with quad doubles
      when others => put_line("  Sorry.  Invalid operation."); return -1;
    end case;
  end Handle_Jobs;

begin
  return Handle_Jobs;
end use_c2mbt;
