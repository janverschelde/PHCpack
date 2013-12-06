with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with Standard_Complex_Numbers_io;       use Standard_Complex_Numbers_io;
with Standard_Complex_Vectors;          use Standard_Complex_Vectors;
with Standard_Complex_Vectors_io;       use Standard_Complex_Vectors_io;

package body Complex_Polynomial_Matrices_io is

  procedure Interactive_get ( lpm : out Link_to_Polynomial_Matrix ) is

    n,m : integer32 := 0;
  
  begin
    put("Give the number of rows : "); get(n);
    put("Give the number of columns : "); get(m);
    declare
      pm : Polynomial_Matrix(1..n,1..m);
    begin
      Interactive_get(pm);
      lpm := new Polynomial_Matrix'(pm);
    end;
  end Interactive_get;

  procedure Interactive_get ( pm : in out Polynomial_Matrix ) is

    d : integer32 := 0;
  
  begin
    for i in pm'range(1) loop
      for j in pm'range(2) loop
        put("Give degree of entry(");
        put(i,1); put(","); put(j,1); put(") : ");
        get(d);
        if d >= 0 then
          pm(i,j) := new Standard_Complex_Vectors.Vector(0..d);
          put("Reading "); put(d+1,1);
          put_line(" complex coefficients...");
          for k in 0..d loop
            put("  cff("); put(k,1); put(") : ");
            get(pm(i,j)(k));
          end loop;
        end if;
      end loop; 
    end loop;
  end Interactive_get;

  procedure get ( lpm : out Link_to_Polynomial_Matrix ) is
  begin
    get(Standard_Input,lpm);
  end get;

  procedure get ( pm : in out Polynomial_Matrix ) is
  begin
    get(Standard_Input,pm);
  end get;

  procedure get ( file : in file_type;
                  lpm : out Link_to_Polynomial_Matrix ) is

    n,m : integer32 := 0;
  
  begin
    get(n); get(m);
    declare
      pm : Polynomial_Matrix(1..n,1..m);
    begin
      get(file,pm);
      lpm := new Polynomial_Matrix'(pm);
    end;
  end get;

  procedure get ( file : in file_type; pm : in out Polynomial_Matrix ) is

    d : integer32 := 0;

  begin
    for i in pm'range(1) loop
      for j in pm'range(2) loop
        get(file,d);
        pm(i,j) := new Standard_Complex_Vectors.Vector(0..d);
        for k in 0..d loop
          get(file,pm(i,j)(k));
        end loop;
      end loop;
    end loop;
  end get;

  procedure put ( pm : in Polynomial_Matrix ) is
  begin
    put(Standard_Output,pm);
  end put;

  procedure put ( file : in file_type; pm : in Polynomial_Matrix ) is
  begin
    put(file,integer32(pm'length(1)),1);
    put(file," ");
    put(file,integer32(pm'length(2)),1);
    new_line(file);
    for i in pm'range(1) loop
      for j in pm'range(2) loop
        if pm(i,j) = null
         then put_line(file,"-1");
         else put(file,pm(i,j)'last,1); new_line(file);
              put_line(file,pm(i,j).all);
        end if;
      end loop;
    end loop;
  end put;

  procedure put ( apm : in Array_of_Polynomial_Matrices ) is
  begin
    put(Standard_Output,apm);
  end put;

  procedure put ( file : in file_type;
                  apm : in Array_of_Polynomial_Matrices ) is
  begin
    for i in apm'range loop
      if apm(i) /= null
       then put(file,apm(i).all);
      end if;
    end loop;
  end put;

end Complex_Polynomial_Matrices_io;
