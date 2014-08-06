with Standard_Natural_Numbers_io;       use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with Monomial_Hashing;                  use Monomial_Hashing;
with QuadDobl_Deflation_Matrices;       use QuadDobl_Deflation_Matrices;

package body QuadDobl_Evaluate_Deflation_io is

  procedure Write_Derivative_Operator
              ( d : in Standard_Natural_Vectors.Vector ) is
  begin
    Write_Derivative_Operator(standard_output,d);
  end Write_Derivative_Operator;

  procedure Write_Derivative_Operator
              ( file : in file_type;
                d : in Standard_Natural_Vectors.Vector ) is
  begin
    put(file,"d(");
    for i in d'range loop
      put(file,d(i),1);
      if i < d'last
       then put(file,",");
       else put(file,")");
      end if;
    end loop;
  end Write_Derivative_Operator;

  procedure Write_Derivative_Operator
              ( d : in Standard_Natural_Vectors.Vector;
                k : in natural32 ) is
  begin
    Write_Derivative_Operator(standard_output,d,k);
  end Write_Derivative_Operator;

  procedure Write_Derivative_Operator
              ( file : in file_type;
                d : in Standard_Natural_Vectors.Vector;
                k : in natural32 ) is
  begin
    Write_Derivative_Operator(file,d);
    put(file,"A("); put(file,k,1); put(file,")");
  end Write_Derivative_Operator;

  procedure Write_Derivative_Operator
               ( d : in Standard_Natural_Vectors.Vector;
                 k,L : in natural32 ) is
  begin
    Write_Derivative_Operator(standard_output,d,k,L);
  end Write_Derivative_Operator;

  procedure Write_Derivative_Operator
               ( file : in file_type;
                 d : in Standard_Natural_Vectors.Vector;
                 k,L : in natural32 ) is
  begin
    for i in 1..L loop
      put(file,"  ");
    end loop;
    Write_Derivative_Operator(file,d,k);
  end Write_Derivative_Operator;

  procedure Write_Spaces ( L : in natural32 ) is
  begin
    Write_Spaces(standard_output,l);
  end Write_Spaces;

  procedure Write_Spaces ( file : in file_type; L : in natural32 ) is
  begin
    for i in 1..L loop
      put(file,"  ");
    end loop;
  end Write_Spaces;

  procedure Write_Zero ( L : in natural32 ) is
  begin
    Write_Zero(standard_output,L);
  end Write_Zero;

  procedure Write_Zero ( file : in file_type; L : in natural32 ) is
  begin
    Write_Spaces(file,L);
    put_line(file,"0");
  end Write_Zero;

  procedure Write ( evt : in Eval_Tree ) is
  begin
    Write(standard_output,evt);
  end Write;

  procedure Write ( file : in file_type; evt : in Eval_Tree ) is

    lt : constant Link_to_Eval_Tree := new Eval_Tree'(evt);

    procedure Write_Lookup ( key : in integer32 ) is

      ft : constant Link_to_Eval_Tree := Look_Up(lt,key);

    begin
      if ft /= null
       then Write_Derivative_Operator(file,ft.d,natural32(ft.m));
       else put(file,"  lookup failed!");
      end if;
    end Write_Lookup;

    procedure Write_Nodes ( t : in Eval_Tree; m : in natural32 ) is

      key : integer32;

    begin
      Write_Derivative_Operator(t.d,natural32(t.m),m);
      Write_Spaces(file,natural32(t.m));
      put(file,"  Node "); put(file,t.key,1);
      key := Key_In(evt,t.d,t.m,integer32(t.key)-1);
      if key = -1
       then put_line(file," does not occur yet.");
       else put(file," occurs as node "); put(file,key,1); put(file," : ");
            Write_Lookup(key); new_line(file);
      end if;
      if t.m > 0 and then not Is_Leaf(t) then
        for i in t.c'range loop
          if t.c(i) = null then
            if t.e(i) > -1 then
              Write_Spaces(file,m);
              put(file,"  occurs as node ");
              put(file,t.e(i),1); put(file," : ");
              Write_Lookup(t.e(i)); new_line(file);
            elsif i = 0 then
              Write_Zero(file,m+1);
            end if;
          else
            Write_Nodes(t.c(i).all,m+1);
          end if;
        end loop;
      end if;
    end Write_Nodes;

  begin
    Write_Nodes(evt,0);
  end Write;

  procedure Write ( evt : in Eval_Tree;
                    nv,nq,R1 : in Standard_Natural_Vectors.Vector ) is
  begin
    Write(standard_output,evt,nv,nq,R1);
  end Write;

  procedure Write ( file : in file_type; evt : in Eval_Tree;
                    nv,nq,R1 : in Standard_Natural_Vectors.Vector ) is

    lt : constant Link_to_Eval_Tree := new Eval_Tree'(evt);
    cnt_jacmat : natural32 := 0;
    cnt_defmat : natural32 := 0;

    procedure Write_Lookup ( key : in integer32 ) is

      ft : constant Link_to_Eval_Tree := Look_Up(lt,key);

    begin
      if ft /= null
       then Write_Derivative_Operator(file,ft.d,natural32(ft.m));
       else put(file,"  lookup failed!");
      end if;
    end Write_Lookup;

    procedure Display_Count ( d : in Standard_Natural_Vectors.Vector;
                              i : in integer32 ) is

      inc,moncnt,nbc : natural32;

    begin
      put("  ");
      if ((d(d'first) = 0) or (i > 0)) then
        nbc := Number_of_Columns(d,nv,R1,natural32(i));
        inc := nq(i)*nbc;
        put(file,nq(i),1); put(file,"-by-"); put(file,nbc,1);
        put(file," matrix");
        put(file," + "); put(file,inc,1); new_line(file);
        cnt_defmat := cnt_defmat + inc;
      else
        moncnt := Monomial_Count(d(d'first),nv(0));
        inc := moncnt*nq(i)*nv(i);
        put(file,moncnt,1); put(file," ");
        put(file,nq(i),1); put(file,"-by-"); put(file,nv(i),1);
        put(file," matrices");
        put(file," + "); put(file,inc,1); new_line(file);
        cnt_jacmat := cnt_jacmat + inc;
      end if;
    end Display_Count;

    procedure Write_Nodes ( t : in Eval_Tree; m : in natural32 ) is

      key : integer32;

    begin
      Write_Derivative_Operator(file,t.d,natural32(t.m),m);
      Write_Spaces(file,natural32(t.m));
      put(file,"  Node "); put(file,t.key,1);
      key := Key_In(evt,t.d,t.m,integer32(t.key)-1);
      if key = -1
       then Display_Count(t.d,t.m);
       else put(file," occurs as node "); put(file,key,1); put(file," : ");
            Write_Lookup(key); new_line(file);
      end if;
      if t.m > 0 and then not Is_Leaf(t) then
        for i in t.c'range loop
          if t.c(i) = null then
            if t.e(i) > -1 then
              Write_Spaces(file,m);
              put(file,"  occurs as node ");
              put(file,t.e(i),1); put(file," : ");
              Write_Lookup(t.e(i)); new_line(file);
            elsif i = 0 then
              Write_Zero(file,m+1);
            end if;
          else
            Write_Nodes(t.c(i).all,m+1);
          end if;
        end loop;
      end if;
    end Write_Nodes;

  begin
    Write_Nodes(evt,0);
    put(file,"#numbers in deflation matrices : ");
    put(file,cnt_defmat,1); new_line(file);
    put(file,"#numbers in Jacobian matrices : ");
    put(file,cnt_jacmat,1); new_line(file);
  end Write;

end QuadDobl_Evaluate_Deflation_io;
