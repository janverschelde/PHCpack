with Interfaces.C;
with text_io;                           use text_io;
with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;       use Standard_Natural_Numbers_io;
with Standard_Floating_Numbers;         use Standard_Floating_Numbers;
with Standard_Floating_Numbers_io;      use Standard_Floating_Numbers_io;
with Standard_Natural_Vectors;          use Standard_Natural_Vectors;
with Standard_Complex_Poly_Systems;
with DoblDobl_Complex_Poly_Systems;
with QuadDobl_Complex_Poly_Systems;
with Standard_Complex_Solutions;
with DoblDobl_Complex_Solutions;
with QuadDobl_Complex_Solutions;
with Standard_Multiplicity_Structure;
with DoblDobl_Multiplicity_Structure;
with QuadDobl_Multiplicity_Structure;
with Assignments_in_Ada_and_C;          use Assignments_in_Ada_and_C;
with Standard_PolySys_Container;
with Standard_Solutions_Container;
with DoblDobl_PolySys_Container;
with DoblDobl_Solutions_Container;
with QuadDobl_PolySys_Container;
with QuadDobl_Solutions_Container;

function use_multip ( job : integer32;
                      a : C_intarrs.Pointer;
                      b : C_intarrs.Pointer;
                      c : C_dblarrs.Pointer ) return integer32 is

  procedure Extract_Parameters
              ( order : out natural32;
                verbose : out boolean;
                tol : out double_float ) is

  -- DESCRIPTION :
  --   Extracts the parameters out of the a, b, c parameters.

  -- ON RETURN :
  --   order    maximum differentiation order;
  --   verbose  if true, then verbose, otherwise silent;
  --   tol      tolerance on the numerical rank.

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

  function Job0 return integer32 is -- in standard double precision

    order : natural32;
    verbose : boolean;
    tol : double_float;
    ls : Standard_Complex_Poly_Systems.Link_to_Poly_Sys
       := Standard_PolySys_Container.Retrieve;
    sols : Standard_Complex_Solutions.Solution_List
         := Standard_Solutions_Container.Retrieve;
    zero : Standard_Complex_Solutions.Link_to_Solution
         := Standard_Complex_Solutions.Head_Of(sols);

    use Standard_Multiplicity_Structure;

  begin
    Extract_Parameters(order,verbose,tol);
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
  end Job0;

  function Job1 return integer32 is -- in double double precision

    order : natural32;
    verbose : boolean;
    tol : double_float;
    ls : DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys
       := DoblDobl_PolySys_Container.Retrieve;
    sols : DoblDobl_Complex_Solutions.Solution_List
         := DoblDobl_Solutions_Container.Retrieve;
    zero : DoblDobl_Complex_Solutions.Link_to_Solution
         := DoblDobl_Complex_Solutions.Head_Of(sols);

    use DoblDobl_Multiplicity_Structure;

  begin
    Extract_Parameters(order,verbose,tol);
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
      Assign(h,b);
    end;
    return 0;
  end Job1;

  function Job2 return integer32 is -- in quad double precision

    order : natural32;
    verbose : boolean;
    tol : double_float;
    ls : QuadDobl_Complex_Poly_Systems.Link_to_Poly_Sys
       := QuadDobl_PolySys_Container.Retrieve;
    sols : QuadDobl_Complex_Solutions.Solution_List
         := QuadDobl_Solutions_Container.Retrieve;
    zero : QuadDobl_Complex_Solutions.Link_to_Solution
         := QuadDobl_Complex_Solutions.Head_Of(sols);

    use QuadDobl_Multiplicity_Structure;

  begin
    Extract_Parameters(order,verbose,tol);
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
  end Job2;

  function Handle_Jobs return integer32 is
  begin
    case job is
      when 0 => return Job0;   -- in standard double precision
      when 1 => return Job1;   -- in double double precision
      when 2 => return Job2;   -- in quad double precision
      when others => put_line("  Sorry.  Invalid operation."); return 1;
    end case;
  end Handle_Jobs;

begin
  return Handle_Jobs;
end use_multip;
