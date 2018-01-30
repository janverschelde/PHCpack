with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Integer_Vectors_io;        use Standard_Integer_Vectors_io;

package body Integer_Permanents is

  procedure Permanent
              ( row : in integer32;
                mat : in Standard_Integer_Matrices.Matrix;
                cols,cnts : in out Standard_Integer_Vectors.Vector;
                per : in out integer64 ) is

    acc : integer64;

  begin
    if row > mat'last(1) then
      acc := 1;
      for i in cols'range loop
        acc := acc*integer64(mat(i,cols(i)));
      end loop;
      per := per + acc;
    else
      for col in mat'range(2) loop
        if cnts(col) > 0 and mat(row,col) /= 0 then
          cnts(col) := 0;
          cols(row) := col;
          Permanent(row+1,mat,cols,cnts,per);
          cnts(col) := 1;
        end if;
      end loop;
    end if;
  end Permanent;

  procedure Permanent
              ( file : in file_type; row : in integer32;
                mat : in Standard_Integer_Matrices.Matrix;
                cols,cnts : in out Standard_Integer_Vectors.Vector;
                per : in out integer64 ) is

    acc : integer64;

  begin
    if row > mat'last(1) then
      acc := 1;
      for i in cols'range loop
        acc := acc*integer64(mat(i,cols(i)));
      end loop;
      per := per + acc;
      put(file,cols); put(file," : "); put(file,acc); new_line(file);
    else
      for col in mat'range(2) loop
        if cnts(col) > 0 and mat(row,col) /= 0 then
          cnts(col) := 0;
          cols(row) := col;
          Permanent(file,row+1,mat,cols,cnts,per);
          cnts(col) := 1;
        end if;
      end loop;
    end if;
  end Permanent;

  procedure Permanent
              ( row : in integer32;
                mat : in Boolean_Matrices.Matrix;
                cols,cnts : in out Standard_Integer_Vectors.Vector;
                per : in out integer64 ) is
  begin
    if row > mat'last(1) then
      per := per + 1;
    else
      for col in mat'range(2) loop
        if cnts(col) > 0 and mat(row,col) then
          cnts(col) := 0;
          cols(row) := col;
          Permanent(row+1,mat,cols,cnts,per);
          cnts(col) := 1;
        end if;
      end loop;
    end if;
  end Permanent;

  procedure Permanent
              ( file : in file_type; row : in integer32;
                mat : in Boolean_Matrices.Matrix;
                cols,cnts : in out Standard_Integer_Vectors.Vector;
                per : in out integer64 ) is
  begin
    if row > mat'last(1) then
      per := per + 1;
      put(file,cols); new_line(file);
    else
      for col in mat'range(2) loop
        if cnts(col) > 0 and mat(row,col) then
          cnts(col) := 0;
          cols(row) := col;
          Permanent(file,row+1,mat,cols,cnts,per);
          cnts(col) := 1;
        end if;
      end loop;
    end if;
  end Permanent;

  function Permanent ( mat : Standard_Integer_Matrices.Matrix )
                     return integer64 is

    res : integer64 := 0;
    cnts : Standard_Integer_Vectors.Vector(mat'range(2))
         := (mat'range(2) => 1);
    cols : Standard_Integer_Vectors.Vector(mat'range(2));

  begin
    Permanent(mat'first(1),mat,cols,cnts,res);
    return res;
  end Permanent;

  function Permanent ( file : file_type;
                       mat : Standard_Integer_Matrices.Matrix )
                     return integer64 is

    res : integer64 := 0;
    cnts : Standard_Integer_Vectors.Vector(mat'range(2))
         := (mat'range(2) => 1);
    cols : Standard_Integer_Vectors.Vector(mat'range(2));

  begin
    Permanent(file,mat'first(1),mat,cols,cnts,res);
    return res;
  end Permanent;

  function Permanent ( mat : Boolean_Matrices.Matrix )
                     return integer64 is

    res : integer64 := 0;
    cnts : Standard_Integer_Vectors.Vector(mat'range(2))
         := (mat'range(2) => 1);
    cols : Standard_Integer_Vectors.Vector(mat'range(2));

  begin
    Permanent(mat'first(1),mat,cols,cnts,res);
    return res;
  end Permanent;

  function Permanent ( file : file_type;
                       mat : Boolean_Matrices.Matrix )
                     return integer64 is

    res : integer64 := 0;
    cnts : Standard_Integer_Vectors.Vector(mat'range(2))
         := (mat'range(2) => 1);
    cols : Standard_Integer_Vectors.Vector(mat'range(2));

  begin
    Permanent(file,mat'first(1),mat,cols,cnts,res);
    return res;
  end Permanent;

end Integer_Permanents;
