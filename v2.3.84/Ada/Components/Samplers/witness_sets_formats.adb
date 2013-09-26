with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io; 
with Standard_Complex_Numbers;          use Standard_Complex_Numbers;
with Standard_Random_Numbers;           use Standard_Random_Numbers;
with Standard_Natural_Vectors;
with Symbol_Table;
with Standard_Complex_Polynomials;      use Standard_Complex_Polynomials;
with Standard_Complex_Poly_Systems_io;  use Standard_Complex_Poly_Systems_io;
with Standard_Complex_Solutions_io;     use Standard_Complex_Solutions_io;
with Witness_Sets_io;                   use Witness_Sets_io;
with Standard_Embed_Polynomials;        use Standard_Embed_Polynomials;
with Planes_and_Polynomials;            use Planes_and_Polynomials;
--with Standard_Affine_Planes;            use Standard_Affine_Planes;
--with Standard_Affine_Solutions;         use Standard_Affine_Solutions;
with Standard_Plane_Representations;    use Standard_Plane_Representations;

package body Witness_Sets_Formats is

  function Embedded_System
              ( n : integer32; b,v : Vector; p : Poly_Sys ) return Poly_Sys is

    res : Poly_Sys(p'first..p'last+n-1);
    hyp : VecVec(1..n-1) := Equations1(b,v);
    et : Term;

  begin
    et.dg := new Standard_Natural_Vectors.Vector'(1..2*n-1 => 0);
    for i in p'range loop
      res(i) := Add_Variables(p(i),natural32(n-1));
      for j in 1..n-1 loop
        et.cf := Random1;
        et.dg(n+j) := 1;
        Add(res(i),et);
        et.dg(n+j) := 0;
      end loop;
    end loop;
    for i in hyp'range loop
      declare
        eh : Vector(0..2*n-1);
      begin
        for j in hyp(i)'range loop
          eh(j) := hyp(i)(j);
        end loop;
        for j in hyp(i)'last+1..eh'last loop
          eh(j) := Create(0.0);
        end loop;
        eh(n+i) := Create(1.0);
        res(p'last+i) := Hyperplane(eh);
      end;
    end loop;
    Clear(hyp); Clear(et);
    return res;
  end Embedded_System;

  function Embedded_System
              ( n,k : integer32; b : Vector; v : VecVec; p : Poly_Sys )
              return Poly_Sys is

    res : Poly_Sys(p'first..p'last+n-k);
    hyp : VecVec(1..n-k) := Equations(b,v);
    et : Term;

  begin
    et.dg := new Standard_Natural_Vectors.Vector'(1..2*n-k => 0);
    for i in p'range loop
      res(i) := Add_Variables(p(i),natural32(n-k));
      for j in 1..n-k loop
        et.cf := Random1;
        et.dg(n+j) := 1;
        Add(res(i),et);
        et.dg(n+j) := 0;
      end loop;
    end loop;
    for i in hyp'range loop
      declare
        eh : Vector(0..2*n-k);
      begin
        for j in hyp(i)'range loop
          eh(j) := hyp(i)(j);
        end loop;
        for j in hyp(i)'last+1..eh'last loop
          eh(j) := Create(0.0);
        end loop;
        eh(n+i) := Create(1.0);
        res(p'last+i) := Hyperplane(eh);
      end;
    end loop;
    Clear(hyp); Clear(et);
    return res;
  end Embedded_System;

  function Embedded_Extrinsic_Solutions
              ( n : integer32; b,v,sol : Vector ) return Solution_List is

    res,res_last : Solution_List;

  begin
    for i in sol'range loop
      declare
        s : Solution(2*n-1);
      begin
        s.t := Create(1.0);
        s.m := 1;
        for j in 1..n loop
          s.v(j) := b(j) + sol(i)*v(j);
        end loop;
        for j in n+1..s.v'last loop
          s.v(j) := Create(0.0);
        end loop;
        s.err := 0.0;
        s.rco := 1.0;
        s.res := 0.0;
        Append(res,res_last,s);
      end;
    end loop;
    return res;
  end Embedded_Extrinsic_Solutions;

  function Embedded_Extrinsic_Solutions
              ( n,k : integer32; b : Vector; v : VecVec;
                sols : Solution_List ) return Solution_List is

    res,res_last : Solution_List;
    tmp : Solution_List := sols;

  begin
    while not Is_Null(tmp) loop
      declare
        ls : constant Link_to_Solution := Head_Of(tmp);
        s : Solution(2*n-k);
      begin
        s.t := Create(1.0);
        s.m := 1;
        for i in 1..n loop
          s.v(i) := b(i);
          for j in 1..k loop
            s.v(i) := s.v(i) + ls.v(j)*v(j)(i);
          end loop;
        end loop;
        for j in n+1..s.v'last loop
          s.v(j) := Create(0.0);
        end loop;
        s.err := 0.0;
        s.rco := 1.0;
        s.res := 0.0;
        Append(res,res_last,s);
      end;
      tmp := Tail_Of(tmp);
    end loop;
    return res;
  end Embedded_Extrinsic_Solutions;

  procedure Write_Embedding
              ( file : in file_type; p : in Poly_Sys;
                n : in integer32; b,v : in Vector;
                sols : in Solution_List ) is

    ep : constant Poly_Sys(p'first..p'last+n-1) := Embedded_System(n,b,v,p);
    esols : constant Solution_List
          := Embedded_Extrinsic_Solutions(n,b,v,Head_Of(sols).v);

  begin
    if Symbol_Table.Number < natural32(2*n-1)
     then Add_Embed_Symbols(natural32(n-1));
    end if;
    new_line(file);
    put(file,"THE EMBEDDED SYSTEM for dimension ");
    put(file,n-1,1);  put_line(file," : ");
    new_line(file);
    put_line(file,ep);
    new_line(file);
    put_line(file,"THE SOLUTIONS : ");
    put(file,Length_Of(esols),natural32(2*n-1),esols);
    new_line(file);
  end Write_Embedding;

  procedure Write_Embedding
              ( file : in file_type; p : in Poly_Sys;
                n,k : in integer32; b : in Vector; v : in VecVec;
                sols : in Solution_List ) is

    ep : constant Poly_Sys(p'first..p'last+n-k) := Embedded_System(n,k,b,v,p);
    esols : constant Solution_List
          := Embedded_Extrinsic_Solutions(n,k,b,v,sols);

  begin
    if Symbol_Table.Number < natural32(n+k)
     then Add_Embed_Symbols(natural32(k));
    end if;
    new_line(file);
    put(file,"THE EMBEDDED SYSTEM for dimension ");
    put(file,n-k,1);  put_line(file," : ");
    new_line(file);
    put_line(file,ep);
    new_line(file);
    put_line(file,"THE SOLUTIONS : ");
    put(file,Length_Of(esols),natural32(2*n-k),esols);
    new_line(file);
  end Write_Embedding;

end Witness_Sets_Formats;
