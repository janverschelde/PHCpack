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

  function Number_of_Start_Columns
             ( dim,nbrows : integer32 ) return integer32 is

    res : integer32 := dim;
 
  begin
    for k in 1..nbrows-1 loop
      res := res*(dim-k);
    end loop;
    return res;
  end Number_of_Start_Columns;

  procedure Start_Columns
              ( row,nbrows : in integer32;
                mat : in Standard_Integer_Matrices.Matrix;
                cols,cnts : in out Standard_Integer_Vectors.Vector;
                idx : in out integer32;
                stc : in out Standard_Integer_VecVecs.VecVec ) is
  begin
    if row > nbrows then
      idx := idx + 1;
      stc(idx) := new Standard_Integer_Vectors.Vector'(cols);
    else
      for col in mat'range(2) loop
        if cnts(col) > 0 and mat(row,col) /= 0 then
          cnts(col) := 0;
          cols(row) := col;
          Start_Columns(row+1,nbrows,mat,cols,cnts,idx,stc);
          cnts(col) := 1;
        end if;
      end loop;
    end if;
  end Start_Columns;

  procedure Start_Columns
              ( row,nbrows : in integer32;
                mat : in Boolean_Matrices.Matrix;
                cols,cnts : in out Standard_Integer_Vectors.Vector;
                idx : in out integer32;
                stc : in out Standard_Integer_VecVecs.VecVec ) is
  begin
    if row > nbrows then
      idx := idx + 1;
      stc(idx) := new Standard_Integer_Vectors.Vector'(cols);
    else
      for col in mat'range(2) loop
        if cnts(col) > 0 and mat(row,col) then
          cnts(col) := 0;
          cols(row) := col;
          Start_Columns(row+1,nbrows,mat,cols,cnts,idx,stc);
          cnts(col) := 1;
        end if;
      end loop;
    end if;
  end Start_Columns;

  function Start_Permanent
             ( row : integer32;
               stc : Standard_Integer_VecVecs.VecVec;
               mat : Standard_Integer_Matrices.Matrix )
             return integer64 is

    res,acc : integer64 := 0;
    cnts : Standard_Integer_Vectors.Vector(mat'range(2));
    cols : Standard_Integer_Vectors.Vector(mat'range(2));

  begin
    for k in stc'range loop
      for i in cols'range loop
        cols(i) := stc(k)(i);
      end loop;
      cnts := (mat'range(2) => 1);
      for col in 1..(row-1) loop
        cnts(cols(col)) := 0;
      end loop;
      acc := 0;
      Permanent(row,mat,cols,cnts,acc);
      res := res + acc;
    end loop;
    return res;
  end Start_Permanent;

  function Start_Permanent
             ( file : file_type; row : integer32;
               stc : Standard_Integer_VecVecs.VecVec;
               mat : Standard_Integer_Matrices.Matrix )
             return integer64 is

    res,acc : integer64 := 0;
    cnts : Standard_Integer_Vectors.Vector(mat'range(2));
    cols : Standard_Integer_Vectors.Vector(mat'range(2));

  begin
    for k in stc'range loop
      for i in cols'range loop
        cols(i) := stc(k)(i);
      end loop;
      cnts := (mat'range(2) => 1);
      for col in 1..(row-1) loop
        cnts(cols(col)) := 0;
      end loop;
      acc := 0;
      Permanent(file,row,mat,cols,cnts,acc);
      res := res + acc;
    end loop;
    return res;
  end Start_Permanent;

  function Start_Permanent
             ( row : integer32;
               stc : Standard_Integer_VecVecs.VecVec;
               mat : Boolean_Matrices.Matrix )
             return integer64 is

    res,acc : integer64 := 0;
    cnts : Standard_Integer_Vectors.Vector(mat'range(2));
    cols : Standard_Integer_Vectors.Vector(mat'range(2));

  begin
    for k in stc'range loop
      for i in cols'range loop
        cols(i) := stc(k)(i);
      end loop;
      cnts := (mat'range(2) => 1);
      for col in 1..(row-1) loop
        cnts(cols(col)) := 0;
      end loop;
      acc := 0;
      Permanent(row,mat,cols,cnts,acc);
      res := res + acc;
    end loop;
    return res;
  end Start_Permanent;

  function Start_Permanent
             ( file : file_type; row : integer32;
               stc : Standard_Integer_VecVecs.VecVec;
               mat : Boolean_Matrices.Matrix )
             return integer64 is

    res,acc : integer64 := 0;
    cnts : Standard_Integer_Vectors.Vector(mat'range(2));
    cols : Standard_Integer_Vectors.Vector(mat'range(2));

  begin
    for k in stc'range loop
      for i in cols'range loop
        cols(i) := stc(k)(i);
      end loop;
      cnts := (mat'range(2) => 1);
      for col in 1..(row-1) loop
        cnts(cols(col)) := 0;
      end loop;
      acc := 0;
      Permanent(file,row,mat,cols,cnts,acc);
      res := res + acc;
    end loop;
    return res;
  end Start_Permanent;

end Integer_Permanents;
