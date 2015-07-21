with Interfaces.C;
with text_io;                            use text_io;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Double_Double_Numbers;              use Double_Double_Numbers;
with Quad_Double_Numbers;                use Quad_Double_Numbers;
with Multprec_Floating_Numbers;          use Multprec_Floating_Numbers;
with Standard_Complex_Numbers;
with Standard_Complex_Vectors;
with Standard_Complex_Poly_Systems;
with Standard_Complex_Solutions;
with Standard_PolySys_Container;
with Standard_Solutions_Container;
with Standard_Scaling;
with DoblDobl_Complex_Numbers;
with DoblDobl_Complex_Vectors;
with DoblDobl_Complex_Poly_Systems;
with DoblDobl_Complex_Solutions;
with DoblDobl_PolySys_Container;
with DoblDobl_Solutions_Container;
with DoblDobl_Scaling;
with QuadDobl_Complex_Numbers;
with QuadDobl_Complex_Vectors;
with QuadDobl_Complex_Vectors_cv;
with QuadDobl_Complex_Poly_Systems;
with QuadDobl_Complex_Solutions;
with QuadDobl_PolySys_Container;
with QuadDobl_Solutions_Container;
with QuadDobl_Scaling;
with Multprec_Complex_Numbers;
with Multprec_Complex_Vectors;
with Multprec_Complex_Poly_Systems;
with Multprec_Complex_Solutions;
with Multprec_PolySys_Container;
with Multprec_Solutions_Container;
with Multprec_Scaling;
with Assignments_in_Ada_and_C;           use Assignments_in_Ada_and_C;

function use_scaling ( job : integer32;
                       a : C_intarrs.Pointer;
                       b : C_intarrs.Pointer;
                       c : C_dblarrs.Pointer ) return integer32 is

  function Job1 return integer32 is -- scale standard system

  -- DESCRIPTION :
  --   Extracts from a the type of scaling, retrieves the
  --   system from the container and applies the scaling,
  --   in standard double precision.

    v_a : constant C_Integer_Array := C_intarrs.Value(a);
    stp : constant integer32 := integer32(v_a(v_a'first));
    rco : double_float;
    lp : Standard_Complex_Poly_Systems.Link_to_Poly_Sys
       := Standard_PolySys_Container.Retrieve;
    scf : Standard_Complex_Vectors.Vector(1..2*lp'last+1);

  begin
    if stp = 0 then
      Standard_Scaling.Scale(lp.all);
    else
      if stp = 1
       then Standard_Scaling.Scale(lp.all,10,false,rco,scf(1..2*lp'last));
       else Standard_Scaling.Scale(lp.all,10,true,rco,scf(1..2*lp'last));
      end if;
      scf(scf'last) := Standard_Complex_Numbers.Create(rco);
      Assign(scf,c);
    end if;
    return 0;
  exception
    when others => return 1;
  end Job1;

  function Job2 return integer32 is -- scale dobldobl system

  -- DESCRIPTION :
  --   Extracts from a the type of scaling, retrieves the
  --   system from the container and applies the scaling,
  --   in double double precision.

    v_a : constant C_Integer_Array := C_intarrs.Value(a);
    stp : constant integer32 := integer32(v_a(v_a'first));
    rco : double_double;
    lp : DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys
       := DoblDobl_PolySys_Container.Retrieve;
    scf : DoblDobl_Complex_Vectors.Vector(1..2*lp'last+1);

  begin
    if stp = 0 then
      DoblDobl_Scaling.Scale(lp.all);
    else
      if stp = 1
       then DoblDobl_Scaling.Scale(lp.all,10,false,rco,scf(1..2*lp'last));
       else DoblDobl_Scaling.Scale(lp.all,10,true,rco,scf(1..2*lp'last));
      end if;
      scf(scf'last) := DoblDobl_Complex_Numbers.Create(rco);
      Assign(scf,c);
    end if;
    return 0;
  exception
    when others => return 2;
  end Job2;

  function Job3 return integer32 is -- scale quaddobl system

  -- DESCRIPTION :
  --   Extracts from a the type of scaling, retrieves the
  --   system from the container and applies the scaling,
  --   in quad double precision.

    v_a : constant C_Integer_Array := C_intarrs.Value(a);
    stp : constant integer32 := integer32(v_a(v_a'first));
    rco : quad_double;
    lp : QuadDobl_Complex_Poly_Systems.Link_to_Poly_Sys
       := QuadDobl_PolySys_Container.Retrieve;
    scf : QuadDobl_Complex_Vectors.Vector(1..2*lp'last+1);

  begin
    if stp = 0 then
      QuadDobl_Scaling.Scale(lp.all);
    else
      if stp = 1
       then QuadDobl_Scaling.Scale(lp.all,10,false,rco,scf(1..2*lp'last));
       else QuadDobl_Scaling.Scale(lp.all,10,true,rco,scf(1..2*lp'last));
      end if;
      scf(scf'last) := QuadDobl_Complex_Numbers.Create(rco);
      Assign(scf,c);
    end if;
    return 0;
  exception
    when others => return 3;
  end Job3;

  function Job4 return integer32 is -- scale multprec system

  -- DESCRIPTION :
  --   Extracts from a the type of scaling, retrieves the
  --   system from the container and applies the scaling,
  --   in arbitrary multiprecision.

    use QuadDobl_Complex_Vectors_cv;

    v_a : constant C_Integer_Array := C_intarrs.Value(a);
    stp : constant integer32 := integer32(v_a(v_a'first));
    rco : Floating_Number;
    lp : Multprec_Complex_Poly_Systems.Link_to_Poly_Sys
       := Multprec_PolySys_Container.Retrieve;
    scf : Multprec_Complex_Vectors.Vector(1..2*lp'last+1);
    qd_scf : QuadDobl_Complex_Vectors.Vector(scf'range);

  begin
    if stp = 0 then
      Multprec_Scaling.Scale(lp.all);
    else
      if stp = 1
       then Multprec_Scaling.Scale(lp.all,10,false,rco,scf(1..2*lp'last));
       else Multprec_Scaling.Scale(lp.all,10,true,rco,scf(1..2*lp'last));
      end if;
      scf(scf'last) := Multprec_Complex_Numbers.Create(rco);
      qd_scf := Multprec_to_QuadDobl_Complex(scf);
      Assign(qd_scf,c);
      Multprec_Complex_Vectors.Clear(scf);
    end if;
    return 0;
  exception
    when others => return 4;
  end Job4;

  function Job5 return integer32 is -- scale standard solutions

  -- DESCRIPTION :
  --   Extracts the dimension, basis, and scaling coefficients
  --   of the parameters a, b, and c.  Then scales the solutions
  --   in the standard solutions container.

    use Standard_Complex_Solutions;

    v_a : constant C_Integer_Array := C_intarrs.Value(a);
    dim : constant integer32 := integer32(v_a(v_a'first));
    v_b : constant C_Integer_Array := C_intarrs.Value(a);
    bas : constant natural32 := natural32(v_a(v_a'first));
    sols : Solution_List := Standard_Solutions_Container.Retrieve;
    scf : Standard_Complex_Vectors.Vector(1..dim/2); 
 
  begin
    Assign(natural32(dim),c,scf);
    Standard_Scaling.Scale(bas,scf,sols);
    return 0;
  exception
    when others => return 5;
  end Job5;

  function Job6 return integer32 is -- scale dobldobl solutions

  -- DESCRIPTION :
  --   Extracts the dimension, basis, and scaling coefficients
  --   of the parameters a, b, and c.  Then scales the solutions
  --   in the dobldobl solutions container.

    use DoblDobl_Complex_Solutions;

    v_a : constant C_Integer_Array := C_intarrs.Value(a);
    dim : constant integer32 := integer32(v_a(v_a'first));
    v_b : constant C_Integer_Array := C_intarrs.Value(a);
    bas : constant natural32 := natural32(v_a(v_a'first));
    sols : Solution_List := DoblDobl_Solutions_Container.Retrieve;
    scf : DoblDobl_Complex_Vectors.Vector(1..dim/4); 
 
  begin
    Assign(natural32(dim),c,scf);
    DoblDobl_Scaling.Scale(bas,scf,sols);
    return 0;
  exception
    when others => return 6;
  end Job6;

  function Job7 return integer32 is -- scale quaddobl solutions

  -- DESCRIPTION :
  --   Extracts the dimension, basis, and scaling coefficients
  --   of the parameters a, b, and c.  Then scales the solutions
  --   in the quaddobl solutions container.

    use QuadDobl_Complex_Solutions;

    v_a : constant C_Integer_Array := C_intarrs.Value(a);
    dim : constant integer32 := integer32(v_a(v_a'first));
    v_b : constant C_Integer_Array := C_intarrs.Value(a);
    bas : constant natural32 := natural32(v_a(v_a'first));
    sols : Solution_List := QuadDobl_Solutions_Container.Retrieve;
    scf : QuadDobl_Complex_Vectors.Vector(1..dim/8); 
 
  begin
    Assign(natural32(dim),c,scf);
    QuadDobl_Scaling.Scale(bas,scf,sols);
    return 0;
  exception
    when others => return 7;
  end Job7;

  function Job8 return integer32 is -- scale multprec solutions

  -- DESCRIPTION :
  --   Extracts the dimension, basis, and scaling coefficients
  --   of the parameters a, b, and c.  Then scales the solutions
  --   in the multprec solutions container.

    use QuadDobl_Complex_Vectors_cv;
    use Multprec_Complex_Solutions;

    v_a : constant C_Integer_Array := C_intarrs.Value(a);
    dim : constant integer32 := integer32(v_a(v_a'first));
    v_b : constant C_Integer_Array := C_intarrs.Value(a);
    bas : constant natural32 := natural32(v_a(v_a'first));
    sols : Solution_List := Multprec_Solutions_Container.Retrieve;
    scf : QuadDobl_Complex_Vectors.Vector(1..dim/8); 
    mp_scf : Multprec_Complex_Vectors.Vector(scf'range); 
 
  begin
    Assign(natural32(dim),c,scf);
    mp_scf := QuadDobl_Complex_to_Multprec(scf);
    Multprec_Scaling.Scale(bas,mp_scf,sols);
    Multprec_Complex_Vectors.Clear(mp_scf);
    return 0;
  exception
    when others => return 8;
  end Job8;

  function Handle_Jobs return integer32 is
  begin
    case job is
      when 1 => return Job1;  -- scale system in the standard container
      when 2 => return Job2;  -- scale system in the dobldobl container
      when 3 => return Job3;  -- scale system in the quaddobl container
      when 4 => return Job4;  -- scale system in the multprec container
      when 5 => return Job5;  -- scale solutions in the standard container
      when 6 => return Job6;  -- scale solutions in the dobldobl container
      when 7 => return Job7;  -- scale solutions in the quaddobl container
      when 8 => return Job8;  -- scale solutions in the multprec container
      when others => put_line("  Sorry.  Invalid operation."); return 1;
    end case;
  exception
    when others => put("Exception raised in use_scaling handling job ");
                   put(job,1); put_line(".  Will not ignore."); raise;
  end Handle_Jobs;

begin
  return Handle_Jobs;
exception
  when others => put_line("Ignoring the exception, returning job number.");
                 return job;
end use_scaling;
