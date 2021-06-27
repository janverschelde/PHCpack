with text_io;                           use text_io;
with Communications_with_User;          use Communications_with_User;
with String_Splitters;                  use String_Splitters;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with Characters_and_Numbers;            use Characters_and_Numbers;
with Standard_Integer_Vectors;
with Standard_Integer_Vectors_io;       use Standard_Integer_Vectors_io;
with Standard_Integer_VecVecs;
with Standard_Integer_VecVecs_io;       use Standard_Integer_VecVecs_io;
with Standard_Complex_Poly_Systems;     use Standard_Complex_Poly_Systems;
with Standard_Complex_Poly_Systems_io;  use Standard_Complex_Poly_Systems_io;
with Supports_of_Polynomial_Systems;
with Lists_of_Integer_Vectors;          use Lists_of_Integer_Vectors;
with Arrays_of_Integer_Vector_Lists;    use Arrays_of_Integer_Vector_Lists;
with Arrays_of_Integer_Vector_Lists_io; use Arrays_of_Integer_Vector_Lists_io;
with Interfaces.C;
with C_Integer_Arrays,C_to_Ada_Arrays;  use C_Integer_Arrays,C_to_Ada_Arrays;

procedure ts_mv2c is

-- DESCRIPTION :
--   This test procedure is a prototype to interface with other
--   software to compute mixed volumes, written in C.

  function mv1 ( ns : integer32; fn : C_Integer_Array;
                 n,m : integer32; ind,cnt,sup : C_Integer_Array )
               return integer32;
  pragma Import(C,mv1,"mv1");

  function Mixed_Volume
             ( file : string; n,m : integer32;
               ind,cnt,sup : Standard_Integer_Vectors.Vector )
             return integer32 is

  -- DESCRIPTION :
  --   Returns the mixed volume of the polytopes spanned by given supports.

  -- ON ENTRY :
  --   file    name of output file for mixed-cell configuration;
  --   n       ambient dimension, length of the vectors in supports;
  --   m       total number of points in the supports;
  --   ind     ind(i) is the start of the i-th support;
  --   cnt     cnt(i) counts the length of the i-th support;
  --   sup     coordinates of the points in the supports.

    ptr : natural := 0;
    c_ind : constant C_Integer_Array := Convert(ind);
    c_cnt : constant C_Integer_Array := Convert(cnt);
    c_sup : constant C_Integer_Array := Convert(sup);
    fnv : Standard_Integer_Vectors.Vector
            (integer32(file'first)..integer32(file'last));
    c_fnv : C_Integer_Array(1..Interfaces.C.size_t(fnv'last));

  begin
    for i in file'range loop
      fnv(integer32(i)) := Character_to_Integer(file(i));
    end loop;
    c_fnv := Convert(fnv);
  --  for i in 1..n loop
  --    put("Points in Support "); put(i,1); put_line(" :");
  --    for j in 1..cnt(i) loop
  --      for k in 1..n loop
  --        ptr := ptr + 1;
  --        put(" "); put(sup(ptr),1);
  --      end loop;
  --      new_line;
  --    end loop;
  --  end loop;
  --  return mv(file,n,m,ind,cnt,sup); -- Ada
    return mv1(fnv'last,c_fnv,n,m,c_ind,c_cnt,c_sup); -- C
  end Mixed_Volume;

  procedure Flatten_Supports
             ( n,m : in integer32; s : in Array_of_Lists;
               ind : out Standard_Integer_Vectors.Vector;
               sup : out Standard_Integer_VecVecs.VecVec ) is

  -- DESCRIPTION :
  --   Flattens the array of lists s into one vector.

  -- ON ENTRY :
  --   n       number of supports, length of the vectors;
  --   m       total number of points in s;
  --   s       n supports of a polynomial system.

  -- ON RETURN :
  --   ind     indices to the start of each support: ind(i) points to
  --           the start of the i-th support s(i) in the vector sup;
  --   sup     contains points in the support.

    ptr : List;
    cnt : integer32 := 0;

  begin
    for i in s'range loop
      ptr := s(i); 
      ind(i) := cnt + 1;
      while not Is_Null(ptr) loop
        cnt := cnt + 1;
        sup(cnt) := Head_Of(ptr);
        ptr := Tail_Of(ptr);
      end loop;
    end loop;
  end Flatten_Supports;

  function Flatten ( n,m : integer32; v : Standard_Integer_VecVecs.VecVec )
                   return Standard_Integer_Vectors.Vector is

  -- DESCRIPTION :
  --   Returns an integer vector of size n*m with the coordinates
  --   of the vectors in v.

    res : Standard_Integer_Vectors.Vector(1..n*m);
    cnt : integer32 := 0;

  begin
    for i in v'range loop
      for j in v(i)'range loop
        cnt := cnt + 1;
        res(cnt) := v(i)(j);
      end loop;
    end loop;
    return res;
  end Flatten;

  procedure Extract_Supports
             ( n : in integer32; p : in Poly_Sys; m : out integer32;
               ind,cnt : out Standard_Integer_Vectors.Vector;
               sup : out Standard_Integer_Vectors.Link_to_Vector ) is

  -- DESCRIPTION :
  --   This procedure extracts the supports of the polynomial system
  --   and then computes the mixed volume of these supports.

  -- ON ENTRY :
  --   n       number of variables and equations in p;
  --   p       a polynomial system.

  -- ON RETURN :
  --   m       total number of points in the supports of p;
  --   ind     ind(i) marks the beginning of the i-th support;
  --   cnt     cnt(i) counts the number of elements in the i-th support;
  --   sup     vector of range 1..n*m with coordinates of all supports.

    s : Array_of_Lists(p'range) := Supports_of_Polynomial_Systems.Create(p);

  begin
   -- new_line;
   -- put_line("The supports of your system :"); put(s);
    m := 0;
    for i in s'range loop
      cnt(i) := integer32(Length_Of(s(i)));
      m := m + cnt(i);
    end loop;
   -- put("Cardinality of supports :"); put(cnt); new_line;
   -- put(" total number of points : "); put(m,1); put_line(".");
    declare
      spv : Standard_Integer_VecVecs.VecVec(1..m);
    begin
      Flatten_Supports(n,m,s,ind,spv);
     -- put_line("Flattened supports :"); put(spv);
     -- put("Index set :"); put(ind); new_line;
      sup := new Standard_Integer_Vectors.Vector'(Flatten(n,m,spv));
    end;
  end Extract_Supports;

  procedure Main is

  -- DESCRIPTION :
  --   Prompts the user for a polynomial system, computes its supports,
  --   and then calls the routine to compute the mixed volume.

    lp : Link_to_Poly_Sys;
    fn : Link_to_String;

  begin
    new_line;
    put_line("Reading a polynomial system ..."); get(lp);
    new_line;
    put_line("Reading the name of the output file...");
    fn := new String'(String_Splitters.Read_String);
    new_line;
   -- put_line("Your polynomial system is "); put(lp.all);
    declare
      n : constant integer32 := lp'last;
      m,mv : integer32;
      cnt,ind : Standard_Integer_Vectors.Vector(1..n);
      sup : Standard_Integer_Vectors.Link_to_Vector;
    begin
      Extract_Supports(n,lp.all,m,ind,cnt,sup);
     -- put("Cardinality of supports :"); put(cnt); new_line;
     -- put(" total number of points : "); put(m,1); new_line;
     -- put("The index set :"); put(ind); new_line;
     -- put_line("The support vector : "); put(sup); new_line;
      mv := Mixed_Volume(fn.all,n,m,ind,cnt,sup.all);
    end;
  end Main;

begin
  Main;
end ts_mv2c;
