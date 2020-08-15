with text_io;                           use text_io;
with Interfaces.C;
with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;       use Standard_Natural_Numbers_io;
with Standard_Floating_Numbers;         use Standard_Floating_Numbers;
with Standard_Floating_Numbers_io;      use Standard_Floating_Numbers_io;
with Standard_Natural_Vectors;
with Standard_Complex_Poly_Systems;
with DoblDobl_Complex_Poly_Systems;
with QuadDobl_Complex_Poly_Systems;
with Standard_Complex_Solutions;
with DoblDobl_Complex_Solutions;
with QuadDobl_Complex_Solutions;
with Standard_Deflation_Methods;
with DoblDobl_Deflation_Methods;
with QuadDobl_Deflation_Methods;
with Standard_Multiplicity_Structure;
with DoblDobl_Multiplicity_Structure;
with QuadDobl_Multiplicity_Structure;
with Assignments_in_Ada_and_C;          use Assignments_in_Ada_and_C;
with Standard_PolySys_Container;
with DoblDobl_PolySys_Container;
with QuadDobl_PolySys_Container;
with Standard_Solutions_Container;
with DoblDobl_Solutions_Container;
with QuadDobl_Solutions_Container;

package body Deflation_Interface is

  function Deflation_Standard_Run
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    use Interfaces.C;
    use Standard_Complex_Poly_Systems;
    use Standard_Complex_Solutions;

    lp : constant Link_to_Poly_Sys := Standard_PolySys_Container.Retrieve;
    sols : constant Solution_List := Standard_Solutions_Container.Retrieve;
    work : Solution_List;
   -- symbolic,output : boolean;
   -- nbdgts : natural32;
    v_a : constant C_Integer_Array
        := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(1));
    maxitr : constant natural32 := natural32(v_a(v_a'first));
    v_b : constant C_Integer_Array
        := C_intarrs.Value(b,Interfaces.C.ptrdiff_t(1));
    maxdef : constant natural32 := natural32(v_b(v_b'first));
    v_c : constant C_Double_Array
        := C_dblarrs.Value(c,Interfaces.C.ptrdiff_t(3));
    tolerr : constant double_float := double_float(v_c(v_c'first));
    tolres : constant double_float := double_float(v_c(v_c'first+1));
    tolrnk : constant double_float := double_float(v_c(v_c'first+2));

  begin
    if vrblvl > 0 then
      put("-> in deflation_interface.");
      put_line("Deflation_Standard_Run ...");
    end if;
    Copy(sols,work);
   -- Deflate_Singularities(lp.all,work);
   -- Set_Default_Parameters
   --   (symbolic,output,maxitr,maxdef,nbdgts,tolerr,tolres,tolrnk);
    Standard_Deflation_Methods.Algorithmic_Deflation_and_Clustering
      (lp.all,work,maxitr,maxdef,tolerr,tolres,tolrnk);
    Standard_Solutions_Container.Clear;
    Standard_Solutions_Container.Initialize(work);
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in deflation_interface.");
        put_line("Deflation_Standard_Run.");
      end if;
      return 196;
  end Deflation_Standard_Run;

  function Deflation_DoblDobl_Run
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    use Interfaces.C;
    use DoblDobl_Complex_Poly_Systems;
    use DoblDobl_Complex_Solutions;

    lp : constant Link_to_Poly_Sys := DoblDobl_PolySys_Container.Retrieve;
    sols : constant Solution_List := DoblDobl_Solutions_Container.Retrieve;
    work : Solution_List;
   -- symbolic,output : boolean;
   -- nbdgts : natural32;
    v_a : constant C_Integer_Array
        := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(1));
    maxitr : constant natural32 := natural32(v_a(v_a'first));
    v_b : constant C_Integer_Array
        := C_intarrs.Value(b,Interfaces.C.ptrdiff_t(1));
    maxdef : constant natural32 := natural32(v_b(v_b'first));
    v_c : constant C_Double_Array
        := C_dblarrs.Value(c,Interfaces.C.ptrdiff_t(3));
    tolerr : constant double_float := double_float(v_c(v_c'first));
    tolres : constant double_float := double_float(v_c(v_c'first+1));
    tolrnk : constant double_float := double_float(v_c(v_c'first+2));

  begin
    if vrblvl > 0 then
      put("-> in deflation_interface.");
      put_line("Deflation_DoblDobl_Run ...");
    end if;
    Copy(sols,work);
   -- Deflate_Singularities(lp.all,work);
   -- Set_Default_Parameters
   --   (symbolic,output,maxitr,maxdef,nbdgts,tolerr,tolres,tolrnk);
    DoblDobl_Deflation_Methods.Algorithmic_Deflation_and_Clustering
      (lp.all,work,maxitr,maxdef,tolerr,tolres,tolrnk);
    DoblDobl_Solutions_Container.Clear;
    DoblDobl_Solutions_Container.Initialize(work);
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in deflation_interface.");
        put_line("Deflation_DoblDobl_Run.");
      end if;
      return 249;
  end Deflation_DoblDobl_Run;

  function Deflation_QuadDobl_Run
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    use Interfaces.C;
    use QuadDobl_Complex_Poly_Systems;
    use QuadDobl_Complex_Solutions;

    lp : constant Link_to_Poly_Sys := QuadDobl_PolySys_Container.Retrieve;
    sols : constant Solution_List := QuadDobl_Solutions_Container.Retrieve;
    work : Solution_List;
   -- symbolic,output : boolean;
   -- nbdgts : natural32;
    v_a : constant C_Integer_Array
        := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(1));
    maxitr : constant natural32 := natural32(v_a(v_a'first));
    v_b : constant C_Integer_Array
        := C_intarrs.Value(b,Interfaces.C.ptrdiff_t(1));
    maxdef : constant natural32 := natural32(v_b(v_b'first));
    v_c : constant C_Double_Array
        := C_dblarrs.Value(c,Interfaces.C.ptrdiff_t(3));
    tolerr : constant double_float := double_float(v_c(v_c'first));
    tolres : constant double_float := double_float(v_c(v_c'first+1));
    tolrnk : constant double_float := double_float(v_c(v_c'first+2));

  begin
    if vrblvl > 0 then
      put("-> in deflation_interface.");
      put_line("Deflation_QuadDobl_Run ...");
    end if;
    Copy(sols,work);
   -- Deflate_Singularities(lp.all,work);
   -- Set_Default_Parameters
   --   (symbolic,output,maxitr,maxdef,nbdgts,tolerr,tolres,tolrnk);
    QuadDobl_Deflation_Methods.Algorithmic_Deflation_and_Clustering
      (lp.all,work,maxitr,maxdef,tolerr,tolres,tolrnk);
    QuadDobl_Solutions_Container.Clear;
    QuadDobl_Solutions_Container.Initialize(work);
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in deflation_interface.");
        put_line("Deflation_QuadDobl_Run.");
      end if;
      return 250;
  end Deflation_QuadDobl_Run;

  procedure Extract_Parameters
              ( a : in C_intarrs.Pointer;
                b : in C_intarrs.Pointer;
                c : in C_dblarrs.Pointer;
                order : out natural32;
                verbose : out boolean;
                tol : out double_float ) is

    v_a : constant C_Integer_Array
        := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(1));
    v_b : constant C_Integer_Array
        := C_intarrs.Value(b,Interfaces.C.ptrdiff_t(1));
    v_c : constant C_Double_Array
        := C_dblarrs.Value(c,Interfaces.C.ptrdiff_t(1));

  begin
    order := natural32(v_a(v_a'first));
    verbose := (integer32(v_b(v_b'first)) = 1);
    tol := double_float(v_c(v_c'first));
    if verbose then
      put("The maximum differentiation order : "); put(order,1); new_line;
      put("Tolerance on the numerical rank : "); put(tol,3); new_line;
    end if;
  end Extract_Parameters;

  function Deflation_Standard_Multiplicity
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    order : natural32;
    verbose : boolean;
    tol : double_float;
    ls : constant Standard_Complex_Poly_Systems.Link_to_Poly_Sys
       := Standard_PolySys_Container.Retrieve;
    sols : constant Standard_Complex_Solutions.Solution_List
         := Standard_Solutions_Container.Retrieve;
    zero : constant Standard_Complex_Solutions.Link_to_Solution
         := Standard_Complex_Solutions.Head_Of(sols);

    use Standard_Multiplicity_Structure;

  begin
    if vrblvl > 0 then
      put("-> in deflation_interface.");
      put_line("Deflation_Standard_Multiplicity ...");
    end if;
    Extract_Parameters(a,b,c,order,verbose,tol);
    declare
      h : Standard_Natural_Vectors.Vector(0..integer32(order));
      sh : Standard_Natural_Vectors.Vector(1..h'last+1);
      m : natural32;
    begin
      if verbose then
        Multiplicity_Structure(standard_output,ls.all,zero.v,tol,order,h,m);
      else
        Multiplicity_Structure(ls.all,zero.v,tol,order,h,m);
      end if;
      Assign(integer32(m),a);
      for k in h'range loop
        sh(k+1) := h(k);
      end loop;
      Assign(sh,b);
    end;
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in deflation_interface.");
        put_line("Deflation_Standard_Multiplicity.");
      end if;
      return 732;
  end Deflation_Standard_Multiplicity;

  function Deflation_DoblDobl_Multiplicity
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    order : natural32;
    verbose : boolean;
    tol : double_float;
    ls : constant DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys
       := DoblDobl_PolySys_Container.Retrieve;
    sols : constant DoblDobl_Complex_Solutions.Solution_List
         := DoblDobl_Solutions_Container.Retrieve;
    zero : constant DoblDobl_Complex_Solutions.Link_to_Solution
         := DoblDobl_Complex_Solutions.Head_Of(sols);

    use DoblDobl_Multiplicity_Structure;

  begin
    if vrblvl > 0 then
      put("-> in deflation_interface.");
      put_line("Deflation_DoblDobl_Multiplicity ...");
    end if;
    Extract_Parameters(a,b,c,order,verbose,tol);
    declare
      h : Standard_Natural_Vectors.Vector(0..integer32(order));
      sh : Standard_Natural_Vectors.Vector(1..h'last+1);
      m : natural32;
    begin
      if verbose then
        Multiplicity_Structure(standard_output,ls.all,zero.v,tol,order,h,m);
      else
        Multiplicity_Structure(ls.all,zero.v,tol,order,h,m);
      end if;
      Assign(integer32(m),a);
      for k in h'range loop
        sh(k+1) := h(k);
      end loop;
      Assign(sh,b);
    end;
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in deflation_interface.");
        put_line("Deflation_DoblDobl_Multiplicity.");
      end if;
      return 733;
  end Deflation_DoblDobl_Multiplicity;

  function Deflation_QuadDobl_Multiplicity
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    order : natural32;
    verbose : boolean;
    tol : double_float;
    ls : constant QuadDobl_Complex_Poly_Systems.Link_to_Poly_Sys
       := QuadDobl_PolySys_Container.Retrieve;
    sols : constant QuadDobl_Complex_Solutions.Solution_List
         := QuadDobl_Solutions_Container.Retrieve;
    zero : constant QuadDobl_Complex_Solutions.Link_to_Solution
         := QuadDobl_Complex_Solutions.Head_Of(sols);

    use QuadDobl_Multiplicity_Structure;

  begin
    if vrblvl > 0 then
      put("-> in deflation_interface.");
      put_line("Deflation_QuadDobl_Multiplicity ...");
    end if;
    Extract_Parameters(a,b,c,order,verbose,tol);
    declare
      h : Standard_Natural_Vectors.Vector(0..integer32(order));
      sh : Standard_Natural_Vectors.Vector(1..h'last+1);
      m : natural32;
    begin
      if verbose then
        Multiplicity_Structure(standard_output,ls.all,zero.v,tol,order,h,m);
      else
        Multiplicity_Structure(ls.all,zero.v,tol,order,h,m);
      end if;
      Assign(integer32(m),a);
      for k in h'range loop
        sh(k+1) := h(k);
      end loop;
      Assign(sh,b);
    end;
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in deflation_interface.");
        put_line("Deflation_QuadDobl_Multiplicity.");
      end if;
      return 734;
  end Deflation_QuadDobl_Multiplicity;

end Deflation_Interface;
