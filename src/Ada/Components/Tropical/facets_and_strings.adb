with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Multprec_Integer_Numbers;           use Multprec_Integer_Numbers;
with Multprec_Integer_Numbers_io;
with Characters_and_Numbers;
with Multprec_Lattice_Supports;

package body Facets_and_Strings is

-- AUXILIARY ROUTINES :

  function Add_Coordinates
             ( v : Standard_Integer_Vectors.Vector; k : integer32;
               accu : string ) return string is

  -- DESCRIPTION :
  --   Appends the value of v(k) to the accumulator accu,
  --   or returns the accumulator if k exceeds the range of v.

  -- REQUIRED : k in v'range.

    strvk : constant string := Characters_and_Numbers.Convert(v(k));

  begin
    if k = v'last
     then return accu & strvk;
     else return Add_Coordinates(v,k+1,accu & strvk & ", ");
    end if;
  end Add_Coordinates;

  function Add_Coordinates
             ( v : Multprec_Integer_Vectors.Vector; k : integer32;
               accu : string ) return string is

  -- DESCRIPTION :
  --   Appends the value of v(k) to the accumulator accu,
  --   or returns the accumulator if k exceeds the range of v.

  -- REQUIRED : k in v'range.

    strvk : constant string
          := Multprec_Integer_Numbers_io.Convert_to_String(v(k));

  begin
    if k = v'last
     then return accu & strvk;
     else return Add_Coordinates(v,k+1,accu & strvk & ", ");
    end if;
  end Add_Coordinates;

  function Labels_of_Neighbors 
             ( f : Multprec_Lattice_3d_Facets.Facet_in_3d )
             return Standard_Integer_Vectors.Vector is

  -- DESCRIPTION :
  --   Returns the vector of the labels of those facets
  --   that are neighbors of f.

    res : Standard_Integer_Vectors.Vector(f.neighbors'range);
    lft : Multprec_Lattice_3d_Facets.Link_to_3d_Facet;

  begin
    for i in f.neighbors'range loop
      lft := f.neighbors(i);
      res(i) := lft.label;
    end loop;
    return res;
  end Labels_of_Neighbors;

  function Labels_of_Neighbors 
             ( f : Multprec_Lattice_4d_Facets.Facet_in_4d )
             return Standard_Integer_Vectors.Vector is

  -- DESCRIPTION :
  --   Returns the vector of the labels of those facets
  --   that are neighbors of f.

    res : Standard_Integer_Vectors.Vector(f.neighbors'range);
    lft : Multprec_Lattice_4d_Facets.Link_to_4d_Facet;

  begin
    for i in f.neighbors'range loop
      lft := f.neighbors(i);
      res(i) := lft.label;
    end loop;
    return res;
  end Labels_of_Neighbors;

-- TARGET FUNCTIONS :

  function write ( v : Standard_Integer_Vectors.Vector ) return string is

    res : constant string := Add_Coordinates(v,v'first,"[");
 
  begin
    return res & "]";
  end write;

  function write ( v : Multprec_Integer_Vectors.Vector ) return string is

    res : constant string := Add_Coordinates(v,v'first,"[");

  begin
    return res & "]";
  end write;

  function write ( A : Multprec_Integer_Matrices.Matrix;
                   f : Multprec_Lattice_3d_Facets.Facet_in_3d )
                 return string is

    val : Integer_Number
        := Multprec_Lattice_3d_Facets.Support_Value_of_Facet(A,f);
    res : constant string
        := "(" & Multprec_Integer_Numbers_io.Convert_to_String(val) & ", "
               & write(f.normal) & ", "
               & write(f.points) & ", "
               & write(Labels_of_Neighbors(f)) & ")";

  begin
    Clear(val);
    return res;
  end write;

  function write ( A : Multprec_Integer_Matrices.Matrix;
                   f : Multprec_Lattice_4d_Facets.Facet_in_4d )
                 return string is

    val : Integer_Number
        := Multprec_Lattice_4d_Facets.Support_Value_of_Facet(A,f);
    res : constant string
        := "(" & Multprec_Integer_Numbers_io.Convert_to_String(val) & ", "
               & write(f.normal) & ", " 
               & write(f.points) & ", "
               & write(Labels_of_Neighbors(f)) & ")";

  begin
    Clear(val);
    return res;
  end write;

end Facets_and_Strings;
