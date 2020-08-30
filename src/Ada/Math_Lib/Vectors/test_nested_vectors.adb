with text_io;                            use text_io;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Complex_Vectors_io;        use Standard_Complex_Vectors_io;
with Standard_Complex_NesVecs_io;        use Standard_Complex_NesVecs_io;

package body Test_Nested_Vectors is
 
  procedure Interactive_Read ( nv : in out NesVec ) is
  begin
    case nv.n is
      when 1 => put("Give "); put(nv.b - nv.a + 1,1);
                put_line(" complex numbers : ");
                get(nv.v);
      when others => for i in nv.w'range loop
                       Interactive_Read(nv.n-1,nv.w(i));
                     end loop;
    end case;
  end Interactive_Read;

  procedure Interactive_Read ( nv : in out Link_to_NesVec ) is

    n : natural32 := 0;

  begin
    put("Give dimension : "); get(n);
    Interactive_Read(n,nv);
  end Interactive_Read;

  procedure Interactive_Read ( n : in natural32;
                               nv : in out Link_to_NesVec ) is

    a,b : integer32 := 0;

  begin
    put("At level "); put(n,1); put(" give start of range : "); get(a);
    put("At level "); put(n,1); put(" give end of range : "); get(b);
    declare
      nv_rep : NesVec(n,a,b);
    begin
      Interactive_Read(nv_rep);
      nv := new NesVec'(nv_rep);
    end;
  end Interactive_Read;

  procedure Main is

    nv,nw : Link_to_NesVec;
    name : string(1..80);
    cnt : natural := 0;
    ipt,opt : file_type;

  begin
    put_line("Reading a multi-dimensional matrix...");
    Interactive_Read(nv);
    put_line("Your multi-dimensional matrix : ");
    put(nv);
    new_line;
    put_line("Testing input/output with files.");
    new_line;
    skip_line;
    put("Give the name of the output file : ");
    get_line(name,cnt);
    create(opt,out_file,name(1..cnt));
    put(opt,nv);
    close(opt);
    open(ipt,in_file,name(1..cnt));
    get(ipt,nw);
    close(ipt);
    put_line("The multi-dimensional matrix read from file : ");
    put(nw);
  end Main;

end Test_Nested_Vectors;
