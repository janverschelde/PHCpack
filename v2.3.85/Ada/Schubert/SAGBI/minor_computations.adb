with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;
with Standard_Integer_Vectors_io;        use Standard_Integer_Vectors_io;
with Standard_Floating_Linear_Solvers;   use Standard_Floating_Linear_Solvers;

package body Minor_Computations is

  function Determinant
              ( mat : Matrix; rows,cols : Standard_Integer_Vectors.Vector )
              return double_float is

  -- DESCRIPTION :
  --   Computes the determinant of the matrix selecting d rows and columns,
  --   where d = rows'length = cols'length.

    res : double_float := 1.0;
    sqm : Matrix(rows'range,rows'range);
    piv : Standard_Integer_Vectors.Vector(rows'range);
    inf : integer32;

  begin
    for i in rows'range loop
      piv(i) := i;
      for j in cols'range loop
        sqm(i,j) := mat(rows(i),cols(j));
      end loop;
    end loop;
    lufac(sqm,rows'last,piv,inf);
    for i in rows'range loop
      res := res*sqm(i,i);
    end loop;
    for i in piv'range loop
      if piv(i) > i
       then res := -res;
      end if;
    end loop;
    return res;
  end Determinant;

  procedure Minors ( file : in file_type;
                     n,m : in natural32; mat : in Matrix ) is

    procedure Minors ( d : in natural32 ) is
  
      rows,cols : Standard_Integer_Vectors.Vector(1..integer32(d));

      procedure Select_Rows ( k,start : in integer32 ) is

        det : double_float;

      begin
        if k > integer32(d) then
          det := Determinant(mat,rows,cols);
          put(file,"Minor ");
          put(file,rows); put(file," X"); put(file,cols);
          put(file," equals ");
          put(file,det); new_line(file);
        else 
          for j in start..integer32(n) loop
            rows(k) := j;
            Select_Rows(k+1,j+1);
          end loop;
        end if;
      end Select_Rows;

      procedure Select_Columns ( k,start : in integer32 ) is
      begin
        if k > integer32(d) then
          Select_Rows(1,1);
        else
          for j in start..integer32(m) loop
            cols(k) := j;
            Select_Columns(k+1,j+1);
          end loop;
        end if;
      end Select_Columns;

    begin
      Select_Columns(1,1);
    end Minors;

  begin
    for d in 2..m loop
      Minors(d);
    end loop;
  end Minors;

  function Number_of_Minors ( n,m : natural32 ) return natural32 is

    res : natural32 := 0;

    procedure Count_Minors ( d : in natural32 ) is
  
      procedure Select_Rows ( k,start : in natural32 ) is
      begin
        if k > d then
          res := res + 1;
        else
          for j in start..n loop
            Select_Rows(k+1,j+1);
          end loop;
        end if;
      end Select_Rows;

      procedure Select_Columns ( k,start : in natural32 ) is
      begin
        if k > d then
          Select_Rows(1,1);
        else
          for j in start..m loop
            Select_Columns(k+1,j+1);
          end loop;
        end if;
      end Select_Columns;

    begin
      Select_Columns(1,1);
    end Count_Minors;

  begin  
    for d in 2..m loop
      Count_Minors(d);
    end loop;
    return res;
  end Number_of_Minors;

  function Sign_of_Minors ( n,m : natural32; mat : Matrix ) return Vector is

    res : Vector(1..integer32(Number_of_Minors(n,m)));
    tol : constant double_float := 10.0**(-10);
    cnt : integer32 := 0;

    procedure Compute_Minors ( d : in natural32 ) is
  
      rows,cols : Standard_Integer_Vectors.Vector(1..integer32(d));

      procedure Select_Rows ( k,start : in integer32 ) is

        det : double_float;

      begin
        if k > integer32(d) then
          det := Determinant(mat,rows,cols);
          cnt := cnt + 1;
          if abs(det) < tol then
            res(cnt) := 0;
          elsif det > 0.0 then
            res(cnt) := +1;
          else
            res(cnt) := -1;
          end if;
        else
          for j in start..integer32(n) loop
            rows(k) := j;
            Select_Rows(k+1,j+1);
          end loop;
        end if;
      end Select_Rows;

      procedure Select_Columns ( k,start : in integer32 ) is
      begin
        if k > integer32(d) then
          Select_Rows(1,1);
        else
          for j in start..integer32(m) loop
            cols(k) := j;
            Select_Columns(k+1,j+1);
          end loop;
        end if;
      end Select_Columns;

    begin
      Select_Columns(1,1);
    end Compute_Minors;

  begin
    for d in 2..m loop
      Compute_Minors(d);
    end loop;
    return res;
  end Sign_of_Minors;

end Minor_Computations;
