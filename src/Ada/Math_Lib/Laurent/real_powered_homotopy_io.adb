with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Characters_and_Numbers;
with Standard_Complex_Laurentials_io;
with Symbol_Table;
with Standard_Complex_Poly_Strings;
with Real_Powered_Series_IO;

package body Real_Powered_Homotopy_IO is

-- WRITE OUTPUT :

  function to_string
             ( d : Standard_Complex_Laurentials.Degrees ) return string is

  -- DESCRIPTION :
  --   Returns the string representation of the terms defined
  --   by the degrees d.

    function write ( k : integer32; accu : string ) return string is

    -- DESCRIPTION :
    --   Recursive write where k is the position in the degrees d.

    begin
      if k > d'last then
        return accu;
      elsif d(k) = 0 then
        return write(k+1,accu);
      else
        declare
          sb : constant Symbol_Table.Symbol := Symbol_Table.Get(natural32(k));
        begin
          if d(k) = 1 then
            declare
              new_accu : constant string := accu & '*'
                & Standard_Complex_Poly_Strings.Concat_Symbol0("",sb);
            begin
              return write(k+1,new_accu);
            end;
          else
            declare
              new_accu : constant string := accu & '*'
                & Standard_Complex_Poly_Strings.Concat_Symbol0("",sb)
                & '^' & Characters_and_Numbers.Convert(d(k));
            begin
              return write(k+1,new_accu);
            end;
          end if;
        end;
      end if;
    end write;

  begin
    return write(d'first,"");
  end to_string;

  function to_string ( q : Standard_Complex_Laurentials.Poly;
                       c : Standard_Complex_VecVecs.VecVec;
                       p : Standard_Floating_VecVecs.VecVec;
                       t : character := 't'; vrblvl : integer32 := 0 )
                     return string is

    terms : array(c'range) of Standard_Complex_Laurentials.Term;
    idx : integer32 := 0;

    procedure Get_Term ( trm : in Standard_Complex_Laurentials.Term;
                         continue : out boolean ) is

    -- DESCRIPTION :
    --   Stores the current term trm in the array terms at position idx
    --   and continues if continue is set to true.

    begin
      idx := idx + 1;
      terms(idx) := trm;
      continue := true;
    end Get_Term;

    procedure Get_Terms is new 
      Standard_Complex_Laurentials.Visiting_Iterator(Get_Term);

    function Write ( k : integer32; accu : string ) return string is

    -- DESCRIPTION :
    --   Recursive function to write the terms, stops when k > c'last,
    --   returning the accumulator accu.

    begin
      if k > c'last then
        return accu;
      else
        declare
          sk : constant string
             := Real_Powered_Series_IO.to_string(c(k).all,p(k).all,t,vrblvl-1);
          str_trm : constant string := to_string(terms(k).dg);
        begin
          if k = 1 then -- do not write the first +
            declare
              new_accu : constant string := accu & "(" & sk & ")";
            begin
              return Write(k+1,new_accu & str_trm);
            end;
          else
            declare
              new_accu : constant string := accu & " + (" & sk & ")";
            begin
              return Write(k+1,new_accu & str_trm);
            end;
          end if;
        end;
      end if;
    end Write;

  begin
    if vrblvl > 0
     then put_line("-> in Real_Powered_Homotopy_IO.to_string ...");
    end if;
    Get_Terms(q);
    return Write(1,"");
  end to_string;

  procedure put ( q : in Standard_Complex_Laurentials.Poly;
                  c : in Standard_Complex_VecVecs.VecVec;
                  p : in Standard_Floating_VecVecs.VecVec;
                  t : in character := 't'; vrblvl : in integer32 := 0 ) is
  begin
    if vrblvl > 0
     then put_line("-> in Real_Powered_Homotopy_IO.put 0 ...");
    end if;
    put(standard_output,q,c,p,t,vrblvl);
  end put;

  procedure put ( file : in file_type;
                  q : in Standard_Complex_Laurentials.Poly;
                  c : in Standard_Complex_VecVecs.VecVec;
                  p : in Standard_Floating_VecVecs.VecVec;
                  t : in character := 't'; vrblvl : in integer32 := 0 ) is

    idx : integer32 := 0;

    procedure Write_Term ( trm : in Standard_Complex_Laurentials.Term;
                           continue : out boolean ) is

    -- DESCRIPTION :
    --   Writes the term trm, and continues if continue is set to true.

    begin
      if idx = 0
       then put(file,"(");
       else put(file," + (");
      end if;
      idx := idx + 1;
      Real_Powered_Series_io.put(file,c(idx).all,p(idx).all,t);
      put(file,")*");
      Standard_Complex_Laurentials_IO.put(file,trm.dg,std => true, pow => '^');
      continue := true;
    end Write_Term;

    procedure Write_Terms is
      new Standard_Complex_Laurentials.Visiting_Iterator(Write_Term);

  begin
    if vrblvl > 0
     then put_line("-> in Real_Powered_Homotopy_IO.put 1 ...");
    end if;
    Write_Terms(q); put(file,";");
  end put;

  procedure put_line ( q : in Standard_Complex_Laurentials.Poly;
                       c : in Standard_Complex_VecVecs.VecVec;
                       p : in Standard_Floating_VecVecs.VecVec;
                       t : in character := 't';
                       vrblvl : in integer32 := 0 ) is

  begin
    if vrblvl > 0
     then put_line("-> in Real_Powered_Homotopy_IO.put_line 0 ...");
    end if;
    put_line(standard_output,q,c,p,t,vrblvl);
  end put_line;

  procedure put_line ( file : in file_type;
                       q : in Standard_Complex_Laurentials.Poly;
                       c : in Standard_Complex_VecVecs.VecVec;
                       p : in Standard_Floating_VecVecs.VecVec;
                       t : in character := 't';
                       vrblvl : in integer32 := 0 ) is

    idx : integer32 := 0;

    procedure Write_Term ( trm : in Standard_Complex_Laurentials.Term;
                           continue : out boolean ) is

    -- DESCRIPTION :
    --   Writes the term trm, and continues if continue is set to true.

    begin
      if idx = 0
       then put(file,"(");
       else put(file," + (");
      end if;
      idx := idx + 1;
      Real_Powered_Series_io.put_line(file,c(idx).all,p(idx).all,t);
      put(file,")*");
      Standard_Complex_Laurentials_IO.put(file,trm.dg,std => true, pow => '^');
      continue := true;
    end Write_Term;

    procedure Write_Terms is
      new Standard_Complex_Laurentials.Visiting_Iterator(Write_Term);

  begin
    if vrblvl > 0
     then put_line("-> in Real_Powered_Homotopy_IO.put_line 1 ...");
    end if;
    Write_Terms(q);
    put_line(file,";");
  end put_line;

end Real_Powered_Homotopy_IO;
