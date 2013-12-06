with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Integer_Vectors;           use Standard_Integer_Vectors;
with Standard_Integer_VecVecs;           use Standard_Integer_VecVecs;

package body Standard_Integer_Transformations_io is

  procedure get ( n : in natural32; t : out Transfo ) is
  begin
    get(Standard_Input,n,t);
  end get;

  procedure get ( file : in file_type; n : in natural32; t : out Transfo ) is

    v : VecVec(1..integer32(n));

  begin
    for i in v'range loop
      v(i) := new Vector(1..integer32(n));
      for j in v(i)'range loop
	get(file,v(i)(j));
      end loop;
    end loop;
    t := Create(v);
    Clear(v);
  end get;

  procedure put ( t : in Transfo ) is
  begin
    put(Standard_Output,t);
  end put;

  procedure put ( file : in file_type; t : in Transfo ) is

    n : constant natural32 := Dimension(t);
    v : VecVec(1..integer32(n));

  begin
    for i in v'range loop
      v(i) := new Vector'(1..integer32(n) => 0);
      v(i)(i) := 1;
      Apply(t,v(i).all);
    end loop;
    for i in v'range loop
      for j in v(i)'range loop
	put(file,' '); put(file,v(j)(i),1);
      end loop;
      new_line(file);
    end loop;
    Clear(v);
  end put;

end Standard_Integer_Transformations_io;
