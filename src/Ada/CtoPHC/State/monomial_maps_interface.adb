with text_io;                           use text_io;
with Interfaces.C;
with Standard_Integer_Vectors;
--with Standard_Integer_Vectors_io;       use Standard_Integer_Vectors_io;
with Standard_Integer_VecVecs;
--with Standard_Integer_VecVecs_io;       use Standard_Integer_VecVecs_io;
with Standard_Complex_Vectors;
--with Standard_Complex_Vectors_io;       use Standard_Complex_Vectors_io;
with Standard_Complex_Laur_Systems;     use Standard_Complex_Laur_Systems;
with Standard_LaurSys_Container;
with Standard_Monomial_Maps;            use Standard_Monomial_Maps;
with Standard_Monomial_Maps_io;
with Black_Box_Binomial_Solvers;        use Black_Box_Binomial_Solvers;
with Monomial_Maps_Container;
with Assignments_in_Ada_and_C;          use Assignments_in_Ada_and_C;

package body Monomial_Maps_Interface is

  function Monomial_Maps_Solve
             ( a : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    lp : constant Link_to_Laur_Sys := Standard_LaurSys_Container.Retrieve;
    sols : Link_to_Array_of_Monomial_Map_Lists;
    fail : boolean;
    val : constant C_Integer_Array := C_intarrs.Value(a);
    puretopdim : constant integer32 := integer32(val(val'first));

  begin
    if vrblvl > 0
     then put_line("-> in monomial_maps_interface.Monomial_Maps_Solve ...");
    end if;
    if lp /= null then
      if puretopdim = 1
       then Black_Box_Binomial_Solver(lp.all,true,sols,fail);
       else Black_Box_Binomial_Solver(lp.all,false,sols,fail);
      end if;
      if not fail
       then Monomial_Maps_Container.Initialize(sols.all);
      end if;
    end if;
    return 0;
  exception
    when others =>
      if vrblvl > 0 then
        put("Exception raised in monomial_maps_interface.");
        put_line("Monomial_Maps_Solve.");
      end if;
    return 430;
  end Monomial_Maps_Solve;

  function Monomial_Maps_Write
             ( vrblvl : integer32 := 0 ) return integer32 is

    sols : constant Link_to_Array_of_Monomial_Map_Lists
         := Monomial_Maps_Container.Retrieve;

  begin
    if vrblvl > 0
     then put_line("-> in monomial_maps_interface.Monomial_Maps_Write ...");
    end if;
    if sols /= null
     then Standard_Monomial_Maps_io.put(sols.all);
    end if;
    return 0;
  exception
    when others =>
      if vrblvl > 0 then
        put("Exception raised in monomial_maps_interface.");
        put_line("Monomial_Maps_Write.");
      end if;
    return 431;
  end Monomial_Maps_Write;

  function Monomial_Maps_Top_Dimension
             ( a : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    topdim : integer32;

  begin
    if vrblvl > 0 then
      put("-> in monomial_maps_interface.");
      put_line("Monomial_Maps_Top_Dimension ...");
    end if;
    topdim := Monomial_Maps_Container.Top_Dimension;
   -- put("The top dimension : "); put(topdim,1); put_line(".");
    Assign(topdim,a);
    return 0;
  exception
    when others =>
      if vrblvl > 0 then
        put("Exception raised in monomial_maps_interface.");
        put_line("Monomial_Maps_Top_Dimension.");
      end if;
    return 433;
  end Monomial_Maps_Top_Dimension;

  function Monomial_Maps_Size
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    v : constant C_Integer_Array := C_intarrs.Value(a);
    d : constant integer32 := integer32(v(v'first));
    n : constant integer32 := Monomial_Maps_Container.Number_of_Maps(d);

  begin
    if vrblvl > 0
     then put_line("-> in monomial_maps_interface.Monomial_Maps_Size ...");
    end if;
   -- put("The number of maps of dimension "); put(d,1); put(" : ");
   -- put(n,1); put_line(".");
    Assign(n,b);
    return 0;
  exception
    when others =>
      if vrblvl > 0 then
        put("Exception raised in monomial_maps_interface.");
        put_line("Monomial_Maps_Size.");
      end if;
    return 434;
  end Monomial_Maps_Size;

  function Monomial_Maps_Degree
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    use Interfaces.C;

    v : constant C_Integer_Array := C_intarrs.Value(a);
    dim : constant integer32 := integer32(v(v'first));
    ind : constant integer32 := integer32(v(v'first+1));
    deg : constant integer32 := Monomial_Maps_Container.Degree(dim,ind);

  begin
    if vrblvl > 0
     then put_line("-> in monomial_maps_interface.Monomial_Maps_Degree ...");
    end if;
    Assign(deg,b);
    return 0;
  exception
    when others =>
      if vrblvl > 0 then
        put("Exception raised in monomial_maps_interface.");
        put_line("Monomial_Maps_Degree.");
      end if;
    return 435;
  end Monomial_Maps_Degree;

  function Monomial_Maps_Coefficients
             ( a : C_intarrs.Pointer;
               c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    use Interfaces.C;

    v : constant C_Integer_Array := C_intarrs.Value(a);
    dim : constant integer32 := integer32(v(v'first));
    ind : constant integer32 := integer32(v(v'first+1));
   -- nbv : constant integer32 := integer32(v(v'first+2));
    cff : constant Standard_Complex_Vectors.Vector
        := Monomial_Maps_Container.Coefficients(dim,ind);

  begin
    if vrblvl > 0 then
      put("-> in monomial_maps_interface.");
      put_line("Monomial_Maps_Coefficients ...");
    end if;
   -- put("the coefficients : "); put(cff); new_line;
    Assign(cff,c);
    return 0;
  exception
    when others =>
      if vrblvl > 0 then
        put("Exception raised in monomial_maps_interface.");
        put_line("Monomial_Maps_Clear.");
      end if;
    return 436;
  end Monomial_Maps_Coefficients;

  function Monomial_Maps_Exponents
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    use Interfaces.C;

    v : constant C_Integer_Array := C_intarrs.Value(a);
    dim : constant integer32 := integer32(v(v'first));
    ind : constant integer32 := integer32(v(v'first+1));
    nbv : constant integer32 := integer32(v(v'first+2));
    exp : constant Standard_Integer_VecVecs.VecVec
        := Monomial_Maps_Container.Exponents(dim,ind);
    flat : Standard_Integer_Vectors.Vector(1..dim*nbv);
    cnt : integer32 := 0;

  begin
    if vrblvl > 0 then
      put("-> in monomial_maps_interface.");
      put_line("Monomial_Maps_Exponents ...");
    end if;
   -- put_line("The exponents as vectors of vectors : "); put(exp);
    for i in exp'range loop
      for j in exp(i)'range loop
         cnt := cnt + 1;
         flat(cnt) := exp(i)(j);
      end loop;
    end loop;
   -- put("The flattened exponents : "); put(flat); new_line;
    Assign(flat,b);
    return 0;
  exception
    when others =>
      if vrblvl > 0 then
        put("Exception raised in monomial_maps_interface.");
        put_line("Monomial_Maps_Clear.");
      end if;
    return 437;
  end Monomial_Maps_Exponents;

  function Monomial_Maps_Data
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    use Interfaces.C;

    v : constant C_Integer_Array := C_intarrs.Value(a);
    dim : constant integer32 := integer32(v(v'first));
    ind : constant integer32 := integer32(v(v'first+1));
    nbv : constant integer32 := integer32(v(v'first+2));
    cff : Standard_Complex_Vectors.Vector(1..nbv);
    exp : Standard_Integer_VecVecs.VecVec(1..nbv);
    flat : Standard_Integer_Vectors.Vector(1..dim*nbv);
    cnt : integer32 := 0;

  begin
    if vrblvl > 0
     then put("-> in monomial_maps_interface.Monomial_Maps_Data ...");
    end if;
    Monomial_Maps_Container.Coefficients_and_Exponents(dim,ind,cff,exp);
    Assign(cff,c);
    for i in exp'range loop
      for j in exp(i)'range loop
         cnt := cnt + 1;
         flat(cnt) := exp(i)(j);
      end loop;
    end loop;
    Assign(flat,b);
    return 0;
  exception
    when others =>
      if vrblvl > 0 then
        put("Exception raised in monomial_maps_interface.");
        put_line("Monomial_Maps_Clear.");
      end if;
    return 438;
  end Monomial_Maps_Data;

  function Monomial_Maps_Clear
             ( vrblvl : integer32 := 0 ) return integer32 is
  begin
    if vrblvl > 0
     then put_line("-> in monomial_maps_interface.Monomial_Maps_Clear ...");
    end if;
    Monomial_Maps_Container.Clear;
    return 0;
  exception
    when others =>
      if vrblvl > 0 then
        put("Exception raised in monomial_maps_interface.");
        put_line("Monomial_Maps_Clear.");
      end if;
    return 432;
  end Monomial_Maps_Clear;

end Monomial_Maps_Interface;
