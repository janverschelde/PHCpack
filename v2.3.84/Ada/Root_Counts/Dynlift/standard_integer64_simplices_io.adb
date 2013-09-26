with Standard_Natural_Numbers_io;       use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with Standard_Integer64_Vectors_io;     use Standard_Integer64_Vectors_io;
with Standard_Integer64_VecVecs;        use Standard_Integer64_VecVecs;
with Standard_Integer64_VecVecs_io;     use Standard_Integer64_VecVecs_io;

package body Standard_Integer64_Simplices_io is

  procedure get ( s : in out Simplex ) is

    n : natural32 := 0;

  begin
    get(n);
    declare
      v : VecVec(1..integer32(n));
    begin
      get(n,v);
      s := Create(v);
    end;
  end get;

  procedure get ( n : in natural32; s : in out Simplex ) is
  
    v : VecVec(1..integer32(n));

  begin
    get(n,v);
    s := Create(v);
  end get;

  procedure get ( file : in file_type; s : in out Simplex ) is

    n : natural32;

  begin
    get(file,n);
    declare
      v : VecVec(1..integer32(n));
    begin
      get(file,n,v);
      s := Create(v);
    end;
  end get;

  procedure get ( file : in file_type; n : in natural32; s : in out Simplex ) is
  
    v : VecVec(1..integer32(n));

  begin
    get(file,n,v);
    s := Create(v);
  end get;

  procedure put ( s : in Simplex ) is
  begin
    put(Normal(s)); new_line;
    put(Normal(s)'last,1); new_line;
    put(Vertices(s));
  end put;

  procedure put ( file : in file_type; s : in Simplex ) is
  begin
    put(file,Normal(s)); new_line(file);
    put(file,Normal(s)'last,1); new_line(file);
    put(file,Vertices(s));
  end put;

end Standard_Integer64_Simplices_io;
