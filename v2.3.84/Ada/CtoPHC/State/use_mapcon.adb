with text_io;                           use text_io;
with Interfaces.C;
with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with Standard_Integer_Vectors;
--with Standard_Integer_Vectors_io;       use Standard_Integer_Vectors_io;
with Standard_Integer_VecVecs;
--with Standard_Integer_VecVecs_io;       use Standard_Integer_VecVecs_io;
with Standard_Complex_Vectors;
--with Standard_Complex_Vectors_io;       use Standard_Complex_Vectors_io;
with Standard_Complex_Laur_Systems;     use Standard_Complex_Laur_Systems;
with Laurent_Systems_Container;
with Standard_Monomial_Maps;            use Standard_Monomial_Maps;
with Standard_Monomial_Maps_io;
with Black_Box_Binomial_Solvers;        use Black_Box_Binomial_Solvers;
with Monomial_Maps_Container;
with Assignments_in_Ada_and_C;          use Assignments_in_Ada_and_C;

function use_mapcon ( job : integer32;
                      a : C_intarrs.Pointer;
		      b : C_intarrs.Pointer;
                      c : C_dblarrs.Pointer ) return integer32 is

  function Job0 return integer32 is -- solve binomial system

    lp : constant Link_to_Laur_Sys := Laurent_Systems_Container.Retrieve;
    sols : Link_to_Array_of_Monomial_Map_Lists;
    fail : boolean;
    val : constant C_Integer_Array := C_intarrs.Value(a);
    puretopdim : constant integer32 := integer32(val(val'first));

  begin
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
  end Job0;

  function Job1 return integer32 is -- write monomial maps

    sols : constant Link_to_Array_of_Monomial_Map_Lists
         := Monomial_Maps_Container.Retrieve;

  begin
    if sols /= null
     then Standard_Monomial_Maps_io.put(sols.all);
    end if;
    return 0;
  end Job1;

  function Job2 return integer32 is -- clear monomial maps
  begin
    Monomial_Maps_Container.Clear;
    return 0;
  end Job2;

  function Job3 return integer32 is -- return top dimension

    topdim : integer32;

  begin
    topdim := Monomial_Maps_Container.Top_Dimension;
   -- put("The top dimension : "); put(topdim,1); put_line(".");
    Assign(topdim,a);
    return 0;
  end Job3;

  function Job4 return integer32 is -- return #maps of given dimension

    v : constant C_Integer_Array := C_intarrs.Value(a);
    d : constant integer32 := integer32(v(v'first));
    n : constant integer32 := Monomial_Maps_Container.Number_of_Maps(d);

  begin
   -- put("The number of maps of dimension "); put(d,1); put(" : ");
   -- put(n,1); put_line(".");
    Assign(n,b);
    return 0;
  end Job4;

  function Job5 return integer32 is -- return degree of a map

    use Interfaces.C;

    v : constant C_Integer_Array := C_intarrs.Value(a);
    dim : constant integer32 := integer32(v(v'first));
    ind : constant integer32 := integer32(v(v'first+1));
    deg : constant integer32 := Monomial_Maps_Container.Degree(dim,ind);

  begin
    Assign(deg,b);
    return 0;
  end Job5;

  function Job6 return integer32 is -- return coefficients of a map

    use Interfaces.C;

    v : constant C_Integer_Array := C_intarrs.Value(a);
    dim : constant integer32 := integer32(v(v'first));
    ind : constant integer32 := integer32(v(v'first+1));
    nbv : constant integer32 := integer32(v(v'first+2));
    cff : constant Standard_Complex_Vectors.Vector
        := Monomial_Maps_Container.Coefficients(dim,ind);

  begin
   -- put("the coefficients : "); put(cff); new_line;
    Assign(cff,c);
    return 0;
  end Job6;

  function Job7 return integer32 is -- return exponents of a map

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
  end Job7;

  function Job8 return integer32 is -- return coefficients and exponents

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
  end Job8;

  function Handle_Jobs return integer32 is
  begin
    case job is
      when 0 => return Job0; -- solve binomial system
      when 1 => return Job1; -- write monomial maps
      when 2 => return Job2; -- clear monomial maps
      when 3 => return Job3; -- return top dimension
      when 4 => return Job4; -- return #maps of given dimension
      when 5 => return Job5; -- return degree of a map
      when 6 => return Job6; -- return coefficients of a map
      when 7 => return Job7; -- return exponents of a map
      when 8 => return Job8; -- return coefficients and exponents of a map
      when others => put_line("invalid operation"); return 1;
    end case;
  end Handle_Jobs;

begin
  return Handle_Jobs;
end use_mapcon;
