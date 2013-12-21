with Standard_Integer_Numbers;           use Standard_Integer_Numbers;

package body Dictionaries is

-- INITIALIZERS :

  procedure Init_Basis
               ( in_bas,out_bas : in out Standard_Integer_Vectors.Vector ) is

    n : constant integer32 := out_bas'last;

  begin
    for i in in_bas'range loop
      in_bas(i) := n+i;              -- slack variable for ith constraint
    end loop;
    for i in out_bas'range loop
      out_bas(i) := i;
    end loop;
  end Init_Basis;

  function Init_Primal_Dictionary
               ( c : Standard_Floating_Vectors.Vector; a : Matrix )
               return Matrix is

    dic : Matrix(0..a'last(1),a'range(2));

  begin
    for i in c'range loop
      dic(0,i) := -c(i);
    end loop;
    for i in a'range(1) loop
      for j in a'range(2) loop
	dic(i,j) := a(i,j);
      end loop;
    end loop;
    return dic;
  end Init_Primal_Dictionary;

  function Init_Dual_Dictionary
              ( c : Standard_Floating_Vectors.Vector; a : Matrix )
              return Matrix is

    dic :  Matrix(0..a'last(1),a'range(2));

  begin
    for i in c'range loop
      dic(0,i) := -c(i);
    end loop;
    for i in a'range(1) loop
      for j in a'range(2) loop
        dic(i,j) := -a(i,j);
      end loop;
    end loop;
    return dic;
  end Init_Dual_Dictionary;

  procedure Primal_Init
               ( c : in Standard_Floating_Vectors.Vector;
                 a : in Matrix; dic : out Matrix;
                 in_bas,out_bas : in out Standard_Integer_Vectors.Vector ) is
  begin
    dic := Init_Primal_Dictionary(c,a);
    Init_Basis(in_bas,out_bas);
  end Primal_Init;

  procedure Dual_Init
               ( c : in Standard_Floating_Vectors.Vector;
                 a : in Matrix; dic : out Matrix;
                 in_bas,out_bas : in out Standard_Integer_Vectors.Vector ) is
  begin
    dic := Init_Dual_Dictionary(c,a);
    Init_Basis(in_bas,out_bas);
  end Dual_Init;

-- MODIFIERS :

  procedure Primal_Update
               ( dic : in out Matrix;
                 in_bas,out_bas : in out Standard_Integer_Vectors.Vector;
                 eps : in double_float; unbounded : out boolean ) is

    column_index : integer32 := 0;
    min : double_float := 0.0;
    row_index : integer32 := 0;
    temp : double_float;
    tmp : integer32;

  begin
    for i in (dic'first(2)+1)..dic'last(2) loop       -- which unknown enters?
      if dic(0,i) < min
       then min := dic(0,i); column_index := i;
      end if;
    end loop;
    if column_index /= 0 then               -- otherwise optimality is reached
      for i in (dic'first(1)+1)..dic'last(1) loop    -- which unknown leaves?
        temp := dic(i,column_index);
        if abs(temp) > eps then
          temp := dic(i,0)/temp;
          if temp > 0.0 then
            if row_index = 0 then
              row_index := i; min := temp;
            elsif temp < min then
              row_index := i; min := temp;
            end if;
          end if;
        end if;
      end loop;
      if row_index = 0 then                           -- solution is unbounded
        unbounded := true;
      else
        unbounded := false;
        temp := dic(row_index,column_index);                        -- pivot
        for j in dic'range(2) loop                       -- divide pivot row
          dic(row_index,j) := dic(row_index,j) / temp;
        end loop;
        for i in dic'range(1) loop -- update other rows, except pivot column
          if i /= row_index then
            for j in dic'range(2) loop
              if j /= column_index then
                dic(i,j) := dic(i,j)-dic(row_index,j)*dic(i,column_index);
              end if;
            end loop;
          end if;
        end loop;
        for i in dic'range(1) loop                    -- update pivot column
          if i = row_index
           then dic(i,column_index) := 1.0 / temp;
           else dic(i,column_index) := - dic(i,column_index) / temp;
          end if;
        end loop;
        tmp := in_bas(row_index);                -- update basis information
        in_bas(row_index) := out_bas(column_index);
        out_bas(column_index) := tmp;
      end if;
    end if;
  end Primal_Update;

  procedure Dual_Update
                ( dic : in out Matrix;
                  in_bas,out_bas : in out Standard_Integer_Vectors.Vector;
                  eps : in double_float; feasible : out boolean ) is

    column_index : integer32 := 0;
    min : double_float := 0.0;
    row_index : integer32 := 0;
    temp : double_float;
    tmp : integer32;

  begin
    for i in (dic'first(1)+1)..dic'last(1) loop       -- which unknown leaves?
      if dic(i,0) < min
       then min := dic(i,0); row_index := i;
      end if;
    end loop;
    if row_index /= 0 then                  -- otherwise optimality is reached
      for i in (dic'first(2)+1)..dic'last(2) loop    -- which unknown enters?
        temp := dic(row_index,i);
        if abs(temp) > eps and then (temp < 0.0) then
          temp := dic(0,i)/temp;
          if column_index = 0 then
            column_index := i; min := temp;
          elsif temp < min then
            column_index := i; min := temp;
          end if;
        end if;
      end loop;
      if column_index = 0 then                        -- problem is infeasible
        feasible := false;
      else
        feasible := true;
        temp := dic(row_index,column_index);                       -- pivot
        for j in dic'range(2) loop                      -- divide pivot row
          dic(row_index,j) := dic(row_index,j) / temp;
        end loop;
        for i in dic'range(1) loop -- update other rows, except pivot column
          if i /= row_index then
            for j in dic'range(2) loop
              if j /= column_index then
                dic(i,j) := dic(i,j)-dic(row_index,j)*dic(i,column_index);
              end if;
            end loop;
          end if;
        end loop;
        for i in dic'range(1) loop                -- update the pivot column
          if i = row_index
           then dic(i,column_index) := 1.0 / temp;
           else dic(i,column_index) := - dic(i,column_index) / temp;
          end if;
        end loop;
        tmp := in_bas(row_index);            -- update the basis information
        in_bas(row_index) := out_bas(column_index);
        out_bas(column_index) := tmp;
      end if;
    end if;
  end Dual_Update;

-- SELECTORS :

  function Primal_Optimal ( dic : Matrix; eps : double_float ) return boolean is
  begin
    for i in (dic'first(2)+1)..dic'last(2) loop
      if abs(dic(0,i)) > eps then
        if dic(0,i) < 0.0
         then return false;
        end if;
      end if;
    end loop;
    return true;
  end Primal_Optimal;

  function Dual_Optimal ( dic : Matrix; eps : double_float ) return boolean is
  begin
    for i in (dic'first(1)+1)..dic'last(1) loop
      if abs(dic(i,0)) > eps then
        if dic(i,0) < 0.0
         then return false;
        end if;
      end if;
    end loop;
    return true;
  end Dual_Optimal;

  function Optimum ( dic : Matrix ) return double_float is
  begin
    return dic(dic'first(1),dic'first(2));
  end Optimum;

  function Primal_Solution 
                  ( dic : Matrix;
                    in_bas,out_bas : Standard_Integer_Vectors.Vector)
                  return Standard_Floating_Vectors.Vector is

    res : Standard_Floating_Vectors.Vector((dic'first(2)+1)..dic'last(2));
    n : constant integer32 := dic'last(2);

  begin
    for i in in_bas'range loop
      if in_bas(i) <= n
       then res(in_bas(i)) := dic(i,0);
      end if;
    end loop;
    for i in out_bas'range loop
      if out_bas(i) <= n
       then res(out_bas(i)) := 0.0;
      end if;
    end loop;
    return res;
  end Primal_Solution;

  function Dual_Solution 
                ( dic : Matrix;
                  in_bas,out_bas : Standard_Integer_Vectors.Vector)
                return Standard_Floating_Vectors.Vector is

    res : Standard_Floating_Vectors.Vector((dic'first(1)+1)..dic'last(1));
    n : constant integer32 := dic'last(2);

  begin
    for i in in_bas'range loop
      if in_bas(i) > n
       then res(in_bas(i)-n) := dic(i,0);
      end if;
    end loop;
    for i in out_bas'range loop
      if out_bas(i) > n
       then res(out_bas(i)-n) := dic(0,i);
      end if;
    end loop;
    return res;
  end Dual_Solution;

end Dictionaries;
