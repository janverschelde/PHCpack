with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;
with Standard_Natural_Vectors;           use Standard_Natural_Vectors;
with Standard_Natural_Vectors_io;        use Standard_Natural_Vectors_io;
with Standard_Integer_Vectors;           use Standard_Integer_Vectors;
with Standard_Floating_Linear_Solvers;   use Standard_Floating_Linear_Solvers;

procedure Maximal_Minors ( file : in file_type;
                           n,d : in natural32; mat : in Matrix;
                           min,max : out double_float ) is

  function Determinant
              ( mat : Matrix; rows : Standard_Natural_Vectors.Vector )
              return double_float is

  -- DESCRIPTION :
  --   Computes the determinant of the matrix obtained by selecting rows.

    res : double_float := 1.0;
    sqm : Matrix(rows'range,rows'range);
    piv : Standard_Integer_Vectors.Vector(rows'range);
    inf : integer32;

  begin
    for i in rows'range loop
      piv(i) := i;
      for j in rows'range loop
        sqm(i,j) := mat(integer32(rows(i)),j);
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

  procedure Main is

    rows : Standard_Natural_Vectors.Vector(1..integer32(d));
    first : boolean := true;
    mindet,maxdet : double_float;

    procedure Select_Rows ( k,start : in integer32 ) is

      det : double_float;

    begin
      if k > integer32(d) then
        det := Determinant(mat,rows);
        put(file,"Minor "); put(file,rows); put(file," equals ");
        put(file,det); new_line(file);
        det := abs(det);
        if first then
          mindet := det; maxdet := det; first := false;
        else
          if det > maxdet then
            maxdet := det;
          elsif det < mindet then
            mindet := det;
          end if;
         end if;
       else
         for j in start..integer32(n) loop
           rows(k) := natural32(j);
           Select_Rows(k+1,j+1);
         end loop;
      end if;
    end Select_Rows;

  begin
    Select_Rows(1,1);
    put(file,"Min : ");       put(file,mindet,3,3,3);
    put(file,"  Max : ");     put(file,maxdet,3,3,3);
    put(file,"  Max/Min : "); put(file,maxdet/mindet,3,3,3); new_line(file);
    min := mindet; max := maxdet;
  end;

begin
  Main;
end Maximal_Minors;
