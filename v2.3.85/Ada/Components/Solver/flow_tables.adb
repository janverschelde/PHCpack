with Standard_Natural_Numbers;         use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;      use Standard_Natural_Numbers_io;
with Standard_Natural_Vectors;         use Standard_Natural_Vectors;

package body Flow_Tables is

  function Create ( nb_eqs,nb_var : integer32 ) return Flow_Table is

    res : Flow_Table;
    res_rep : Matrix(0..nb_eqs,0..nb_var);

  begin
    for i in res_rep'range(1) loop
      for j in res_rep'range(2) loop
        res_rep(i,j) := 0;
      end loop;
    end loop;
    res := new Matrix'(res_rep);
    return res;
  end Create;

  procedure Store_Hypersurface_Degree
              ( ft : in out Flow_Table; k,d : in integer32 ) is
  begin
    ft(k-1,0) := natural32(d);
  end Store_Hypersurface_Degree;

  procedure Update ( ft : in out Flow_Table; k : in integer32;
                     s : in Array_of_Solution_Lists ) is
  begin
    for j in s'range loop
      exit when j > ft'last(2);
      ft(k,j) := Length_Of(s(j));
    end loop;
  end Update;

  function Max_in_Column ( ft : Flow_Table; k : integer32 ) return natural32 is

  -- DESCRIPTION :
  --   Returns the maximal number in column k.

    res : natural32 := ft(ft'first(1),k);

  begin
    for i in ft'first(1)+1..ft'last(1) loop
      if ft(i,k) > res
       then res := ft(i,k);
      end if;
    end loop;
    return res;
  end Max_in_Column;

  function Number_of_Digits ( n : natural32 ) return natural32 is

  -- DESCRIPTION :
  --   Returns the number of decimal places occupied by the number n.

  begin
    if n < 10
     then return 1;
     else return 1 + Number_of_Digits(n/10);
    end if;
  end Number_of_Digits;

  function Digits_in_Columns ( ft : Flow_Table ) return Vector is

  -- DESCRIPTION :
  --   Returns the number of decimal places the largest number
  --   occupies in every column.

    res : Vector(ft'range(2));
    max : natural32;
  
  begin
    for k in ft'range(2) loop
      max := Max_in_Column(ft,k);
      res(k) := Number_of_Digits(max);
    end loop;
    return res;
  end Digits_in_Columns;

  procedure Write ( file : in file_type; ft : in Flow_Table ) is

    format : constant Vector(ft'range(2)) := Digits_in_Columns(ft);

  begin
    for i in ft'range(1) loop
      for j in 1..ft'last(2) loop
        put(file," "); put(file,ft(i,j),format(j));
      end loop;
      put(file," | ");
      put(file,ft(i,0),format(0));
      new_line(file);
    end loop;
  end Write;

  procedure Clear ( ft : in out Flow_Table ) is
  begin
    Standard_Natural_Matrices.Clear(Link_to_Matrix(ft));
  end Clear;

end Flow_Tables;
