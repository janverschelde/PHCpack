with text_io;                           use text_io;
with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Characters_and_Numbers;
with Multprec_Integer_Numbers;          use Multprec_Integer_Numbers;
with Multprec_Integer_Numbers_io;
with Standard_Integer_Vectors;

package body Point_Lists_and_Strings is

-- AUXILIARY ROUTINES :

  function Add_Coordinates
             ( A : Standard_Integer64_Matrices.Matrix; k,i : integer32;
               accu : string ) return string is

  -- DESCRIPTION :
  --   Returns accu added with the i-th coordinate of the point defined 
  --   by the k-th column of A.

  -- REQUIRED : k <= A'last(2) and i <= A'last(1).

    strAik : constant string := Characters_and_Numbers.Convert(A(i,k));

  begin
    if i = A'last(1) then
      return accu & strAik & ")";
    else
      return Add_Coordinates(A,k,i+1,accu & strAik & ", ");
    end if;
  end Add_Coordinates;

  function Add_Coordinates
             ( A : Multprec_Integer_Matrices.Matrix; k,i : integer32;
               accu : string ) return string is

  -- DESCRIPTION :
  --   Returns accu added with the i-th coordinate of the point defined 
  --   by the k-th column of A.

  -- REQUIRED : k <= A'last(2) and i <= A'last(1).

    strAik : constant string
           := Multprec_Integer_Numbers_io.Convert_to_String(A(i,k));

  begin
    if i = A'last(1) then
      return accu & strAik & ")";
    else
      return Add_Coordinates(A,k,i+1,accu & strAik & ", ");
    end if;
  end Add_Coordinates;

  function Add_Point
             ( A : Standard_Integer64_Matrices.Matrix; k : integer32;
               accu : string ) return string is

  -- DESCRIPTION :
  --   Returns accu added with the point defined by the k-th column of A.

    point : constant string := Add_Coordinates(A,k,A'first(1),"(");

  begin
    if k = A'last(2) then
      return accu & point;
    else
      return Add_Point(A,k+1,accu & point & ", ");
    end if;
  end Add_Point;

  function Add_Point
             ( A : Multprec_Integer_Matrices.Matrix; k : integer32;
               accu : string ) return string is

  -- DESCRIPTION :
  --   Returns accu added with the point defined by the k-th column of A.

    point : constant string := Add_Coordinates(A,k,A'first(1),"(");

  begin
    if k = A'last(2) then
      return accu & point;
    else
      return Add_Point(A,k+1,accu & point & ", ");
    end if;
  end Add_Point;

  procedure Scan ( s : in string; s_pos : in integer;
                   b : in out string; b_last : out integer ) is

  -- DESCRIPTION :
  --   Copies characters of s into b, starting at s_pos.
  --   The copying stops when a comma or closing bracket is encountered.
  --   The value of b_last on return points to the last copied character
  --   in b, so the result is b(b'first..b_last).

    ind : integer := b'first-1;

  begin
    for i in s_pos..s'last loop
      if s(i) /= ',' and s(i) /= ')' then
        ind := ind + 1;
        b(ind) := s(i);
      else
        b_last := ind; return;
      end if;
    end loop;
  end Scan;

  procedure Extract_Numbers
              ( s : in string; rows,cols : in integer32;
                A : out Standard_Integer64_Matrices.Matrix ) is

  -- DESCRIPTION :
  --   Extracts the numbers of the string representation s of the point
  --   configuration into the matrix A, of dimensions rows and cols.

    use Characters_and_Numbers;

    row,col : integer32 := 0;
    buffer : string(s'range);
    buffer_last : integer;

  begin
    for i in s'range loop
      if s(i) = '(' then
        col := col + 1;
        row := A'first(1);
        Scan(s,i+1,buffer,buffer_last);
       -- put("A("); put(row,1); put(","); put(col,1); put(") : ");
       -- put_line(buffer(buffer'first..buffer_last));
        A(row,col) := convert(buffer(buffer'first..buffer_last));
      elsif s(i) = ',' then
        row := row + 1;
        if row <= A'last(1) then
          Scan(s,i+1,buffer,buffer_last);
         -- put("A("); put(row,1); put(","); put(col,1); put(") : ");
         -- put_line(buffer(buffer'first..buffer_last));
          A(row,col) := convert(buffer(buffer'first..buffer_last));
        end if;
      end if;
    end loop;
  end Extract_Numbers;

  procedure Extract_Numbers
              ( s : in string; rows,cols : in integer32;
                A : out Multprec_Integer_Matrices.Matrix ) is

  -- DESCRIPTION :
  --   Extracts the numbers of the string representation s of the point
  --   configuration into the matrix A, of dimensions rows and cols.

    use Multprec_Integer_Numbers_io;

    row,col : integer32 := 0;
    buffer : string(s'range);
    buffer_last : integer;

  begin
    for i in s'range loop
      if s(i) = '(' then
        col := col + 1;
        row := A'first(1);
        Scan(s,i+1,buffer,buffer_last);
       -- put("A("); put(row,1); put(","); put(col,1); put(") : ");
       -- put_line(buffer(buffer'first..buffer_last));
        get(buffer(buffer'first..buffer_last),A(row,col));
      elsif s(i) = ',' then
        row := row + 1;
        if row <= A'last(1) then
          Scan(s,i+1,buffer,buffer_last);
         -- put("A("); put(row,1); put(","); put(col,1); put(") : ");
         -- put_line(buffer(buffer'first..buffer_last));
          get(buffer(buffer'first..buffer_last),A(row,col));
        end if;
      end if;
    end loop;
  end Extract_Numbers;

-- TARGET FUNCTIONS :

  function write ( A : Standard_Integer64_Matrices.Matrix ) return string is

    res : constant string := Add_Point(A,A'first(2),"[");

  begin
    return res & "]";
  end write;

  function write ( A : Multprec_Integer_Matrices.Matrix ) return string is

    res : constant string := Add_Point(A,A'first(2),"[");

  begin
    return res & "]";
  end write;

  function write ( A : Lists_of_Integer_Vectors.List ) return string is

    use Lists_of_Integer_Vectors;
    len : constant integer32 := integer32(Length_Of(A));
    dim : integer32;

  begin
    
    if len = 0 then
      return "";
    else
      dim := Head_Of(A)'last;
      declare
        mat : Standard_Integer64_Matrices.Matrix(1..dim,1..len);
        tmp : List := A;
        lpt : Standard_Integer_Vectors.Link_to_Vector;
        col : integer32 := 0;
      begin
        while not Is_Null(tmp) loop
          lpt := Head_Of(tmp);
          col := col + 1;
          for row in lpt'range loop
            mat(row,col) := integer64(lpt(row));
          end loop;
          tmp := Tail_Of(tmp);
        end loop;
        return write(mat);
      end;
    end if;
  end write;

  procedure Extract_Dimensions ( s : in string; rows,cols : out integer32 ) is
  begin
    rows := 0;
    cols := 0;
    for i in s'range loop
      if s(i) = '(' then
        cols := cols + 1;
      elsif s(i) = ',' then     -- as we count also the comma after )
        if cols = 1             -- the count will equal the dimension
         then rows := rows + 1;
        end if;
      end if;
    end loop;
  end Extract_Dimensions;

  function parse ( s : string; rows,cols : integer32 )
                 return Standard_Integer64_Matrices.Matrix is

    res : Standard_Integer64_Matrices.Matrix(1..rows,1..cols);

  begin
    Extract_Numbers(s,rows,cols,res);
    return res;
  end parse;

  function parse ( s : string; rows,cols : integer32 )
                 return Multprec_Integer_Matrices.Matrix is

    res : Multprec_Integer_Matrices.Matrix(1..rows,1..cols);

  begin
    Extract_Numbers(s,rows,cols,res);
    return res;
  end parse;

  function parse ( n : integer32; s : string )
                 return Standard_Integer_Vectors.Vector is

    use Characters_and_Numbers;

    res : Standard_Integer_Vectors.Vector(1..n);
    ind : integer32 := 0;
    buffer : string(s'range);
    buffer_last : integer;

  begin
    for i in s'range loop
      if s(i) = '(' then
        ind := ind + 1;
        Scan(s,i+1,buffer,buffer_last);
        res(ind) := convert(buffer(buffer'first..buffer_last));
      elsif s(i) = ',' then
        ind := ind + 1;
        if ind <= res'last then
          Scan(s,i+1,buffer,buffer_last);
          res(ind) := convert(buffer(buffer'first..buffer_last));
        end if;
      end if;
    end loop;
    return res;
  end parse;

end Point_Lists_and_Strings;
