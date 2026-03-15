with Interfaces.C;
with Ada.Text_IO;                       use Ada.Text_IO;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Integer_Numbers_IO;       use Standard_Integer_Numbers_IO;
with Standard_Floating_Numbers;         use Standard_Floating_Numbers;
with Standard_Complex_Numbers;          use Standard_Complex_Numbers;
with Standard_Complex_Vectors;
with Standard_Complex_Vectors_IO;       use Standard_Complex_Vectors_IO;
with Standard_Random_Vectors;
with C_Integer_Arrays;                  use C_Integer_Arrays;
with C_Double_Arrays;                   use C_Double_Arrays;
with DCMPLX_VecVecs_Interface;          use DCMPLX_VecVecs_Interface;

procedure ts_use_cvvcon is

-- DESCRIPTION :
--   Tests the interface to the arrays of vectors of complex vectors.

  function Get_Dimension ( vrb : integer32 ) return integer32 is

  -- DESCRIPTION :
  --   Returns the number of arrays in the container.
  --   The verbose level is given in vrb.

    ar : C_Integer_Array(0..Interfaces.C.size_t(1));
    a : constant C_IntArrs.Pointer := ar(0)'unchecked_access;
    r,dim : integer32;

  begin
    r := DCMPLX_VecVecs_Get_Dimension(a,vrb);
    if r /= 0 then
      put("DCMPLX_VecVecs_Get_Dimension returned ");
      put(r,1); new_line;
    end if;
    dim := integer32(ar(0));
    return dim;
  end Get_Dimension;

  function Get_Size ( idx,vrb : integer32 ) return integer32 is

  -- DESCRIPTION :
  --   Returns the size of the vector of array with index idx.
  --   The verbose level is given in vrb.

    ar : C_Integer_Array(0..Interfaces.C.size_t(1));
    br : C_Integer_Array(0..Interfaces.C.size_t(1));
    a : constant C_IntArrs.Pointer := ar(0)'unchecked_access;
    b : constant C_IntArrs.Pointer := br(0)'unchecked_access;
    r,size : integer32;

  begin
    ar(0) := Interfaces.C.int(idx);
    r := DCMPLX_VecVecs_Get_Size(a,b,vrb);
    if r /= 0 then
      put("DCMPLX_VecVecs_Get_Size returned ");
      put(r,1); new_line;
    end if;
    size := integer32(br(0));
    return size;
  end Get_Size;

  procedure Prompt_Dimensions ( vrb : in integer32 ) is

  -- DESCRIPTION :
  --   Prompts for the dimensions and initializes the container.

    n,r,m : integer32 := 0;

  begin
    new_line;
    put("Give the number of arrays : "); get(n);   
    declare
      ar : C_Integer_Array(0..Interfaces.C.size_t(1));
      br : C_Integer_Array(0..Interfaces.C.size_t(n));
      a : constant C_IntArrs.Pointer := ar(0)'unchecked_access;
      b : constant C_IntArrs.Pointer := br(0)'unchecked_access;
    begin
      ar(0) := Interfaces.C.int(n);
      for i in 1..n loop
        put("Give the size of array "); put(i,1); put(" : "); get(m);
        br(Interfaces.C.size_t(i-1)) := Interfaces.C.int(m);
      end loop;
      r := DCMPLX_VecVecs_Initialize(a,b,vrb);
      if r /= 0
       then put("DCMPLX_VecVecs_Initialize returned "); put(r,1); new_line;
      end if;
    end;
    put_line("-> retrieving the stored dimensions ...");
    declare
      dim : constant integer32 := Get_Dimension(vrb);
      size : integer32;
    begin
      put("  dim : "); put(dim,1); new_line;
      for i in 1..dim loop
        size := Get_Size(i,vrb);
        put("  size : "); put(size,1); new_line;
      end loop;
    end;
  end Prompt_Dimensions;

  procedure Add_Random_Vectors ( dim,vrb : in integer32 ) is

    data : Standard_Complex_Vectors.Vector(1..dim);
    anbr : constant integer32 := Get_Dimension(vrb);
    size : integer32;
    ar : C_Integer_Array(0..Interfaces.C.size_t(1));
    br : C_Integer_Array(0..Interfaces.C.size_t(2));
    cr : C_Double_Array(0..Interfaces.C.size_t(2*dim));
    a : constant C_IntArrs.Pointer := ar(0)'unchecked_access;
    b : constant C_IntArrs.Pointer := br(0)'unchecked_access;
    c : constant C_DblArrs.Pointer := cr(0)'unchecked_access;
    r : integer32;
    idx : Interfaces.C.size_t;
    use Interfaces.C;

  begin
    put("adding "); put(anbr,1); put_line(" arrays ...");
    for i in 1..anbr loop
      size := Get_Size(i,vrb);
      for j in 1..size loop
        data := Standard_Random_Vectors.Random_Vector(1,dim);
        put("-> adding "); put(j,1); put(" to array "); put(i,1); 
        put_line(" : "); put_line(data);
        ar(0) := Interfaces.C.int(dim);
        br(0) := Interfaces.C.int(i);
        br(1) := Interfaces.C.int(j);
        idx := 0;
        for k in 1..dim loop
          cr(idx) := Interfaces.C.double(REAL_PART(data(k))); idx := idx+1;
          cr(idx) := Interfaces.C.double(IMAG_PART(data(k))); idx := idx+1;
        end loop;
        r := DCMPLX_VecVecs_Set(a,b,c,vrb);
        if r /= 0
         then put("DCMPLX_VecVecs_Set returned "); put(r,1); new_line;
        end if;
      end loop;
    end loop;
  end Add_Random_Vectors;

  procedure Get_Vectors ( dim,vrb : in integer32 ) is

    data : Standard_Complex_Vectors.Vector(1..dim);
    anbr : constant integer32 := Get_Dimension(vrb);
    size : integer32;
    ar : C_Integer_Array(0..Interfaces.C.size_t(1));
    br : C_Integer_Array(0..Interfaces.C.size_t(2));
    cr : C_Double_Array(0..Interfaces.C.size_t(2*dim-1))
       := (0..Interfaces.C.size_t(2*dim-1) => Interfaces.C.double(0.0));
    a : constant C_IntArrs.Pointer := ar(0)'unchecked_access;
    b : constant C_IntArrs.Pointer := br(0)'unchecked_access;
    c : constant C_DblArrs.Pointer := cr(0)'unchecked_access;
    r : integer32;
    idx : Interfaces.C.size_t;

    use Interfaces.C;

  begin
    for i in 1..anbr loop
      size := Get_Size(i,vrb);
      for j in 1..size loop
        put("-> getting "); put(j,1); put(" from array "); put(i,1); 
        put_line(" ...");
        ar(0) := Interfaces.C.int(dim);
        br(0) := Interfaces.C.int(i);
        br(1) := Interfaces.C.int(j);
        r := DCMPLX_VecVecs_Get(a,b,c,vrb);
        if r /= 0
         then put("DCMPLX_VecVecs_Get returned "); put(r,1); new_line;
        end if;
        idx := 0;
        for k in 1..dim loop
          data(k) := Create(double_float(cr(idx)),double_float(cr(idx+1)));
          idx := idx + 2;
        end loop;
        put("vector "); put(j,1); put(" from array "); put(i,1);
        put_line(" :");
        put_line(data);
      end loop;
    end loop;
  end Get_Vectors;

  procedure Main is

    r : integer32 := 0;
    dim : integer32 := 0;

  begin
    Prompt_Dimensions(10);
    new_line;
    put("Give the dimension of the vectors : "); get(dim);
    Add_Random_Vectors(dim,10);
    new_line;
    put_line("Retrieving the vectors ...");
    Get_Vectors(dim,10);
    r := DCMPLX_VecVecs_Clear;
    if r /= 0
     then put("Double_VecVecs_Clear returned "); put(r,1); new_line;
    end if;
  end Main;

begin 
  Main;
end ts_use_cvvcon;
