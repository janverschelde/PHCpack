with Timing_Package;                    use Timing_Package;
with Standard_Natural_Numbers_io;       use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Floating_Numbers;         use Standard_Floating_Numbers;
with Standard_Complex_Numbers;
with Standard_Complex_Numbers_io;       use Standard_Complex_Numbers_io;
with Double_Double_Numbers;             use Double_Double_Numbers;
with Quad_Double_Numbers;               use Quad_Double_Numbers;
with DoblDobl_Complex_Numbers;
with DoblDobl_Complex_Numbers_io;       use DoblDobl_Complex_Numbers_io;
with QuadDobl_Complex_Numbers;
with QuadDobl_Complex_Numbers_io;       use QuadDobl_Complex_Numbers_io;
with Standard_Natural_Vectors_io;       use Standard_Natural_Vectors_io;
with Standard_Complex_Vectors;
with DoblDobl_Complex_Vectors;
with QuadDobl_Complex_Vectors;
with Standard_Trace_Interpolators;      use Standard_Trace_Interpolators;
with DoblDobl_Trace_Interpolators;      use DoblDobl_Trace_Interpolators;
with QuadDobl_Trace_Interpolators;      use QuadDobl_Trace_Interpolators;

package body Combinatorial_Factorization is

-- AUXILIARIES :

  procedure Write ( file : in file_type; f : in VecVec ) is

  -- DESCRIPTION :
  --   The vectors in f represent labels of witness points,
  --   grouped along their irreducible factors.

    first : boolean := true;

  begin
    for i in f'range loop
      if f(i) /= null then
        if first
         then first := false;
         else put(file," |");
        end if;
        put(file,f(i).all);
      end if;
    end loop;
  end Write;

  function Number_of_Factors ( f : VecVec ) return natural32 is

  -- DESCRIPTION :
  --   Returns the number of nonvoid pointers in the vector f.

    cnt : natural32 := 0;

  begin
    for i in f'range loop
      if f(i) /= null
       then cnt := cnt + 1;
      end if;
    end loop;
    return cnt;
  end Number_of_Factors;

  procedure Write_Labels ( file : in file_type; f : in VecVec ) is
  begin
    put(file,Number_of_Factors(f),1);
    new_line(file);
    for i in f'range loop
      if f(i) /= null then
        put(file,natural32(f(i)'length),1);
        put(file," :");
        put(file,f(i).all);
        new_line(file);
      end if;
    end loop;
  end Write_Labels;

  function Next_Factor ( n : natural32; f : Vector ) return Vector is

  -- DESCRIPTION :
  --   Returns the next factor in the sequence, starting first with the
  --   factors of least cardinality and building up lexicographically.

    res : Vector(f'first..f'last+1);
    ind : integer32;

  begin
    if f(f'last) < n then
      res(f'range) := f;
      res(f'last) := f(f'last)+1;
      return res(f'range);
    else
      ind := f'last-1;
      while ind >= f'first loop
        if f(ind) < f(ind+1)-1 then
          res(f'range) := f; 
          res(ind) := f(ind)+1;
          for j in ind+1..res'last loop
            res(j) := res(j-1)+1;
          end loop;
          return res(f'range);
        end if;
        ind := ind - 1;
      end loop;
      if f'length = integer32(n) then
        return f;
      else
        for i in res'range loop
          res(i) := natural32(i);
        end loop;
        return res;
      end if;
    end if;
  end Next_Factor;

  function Is_In ( v : Vector; k : natural32 ) return boolean is

  -- DESCRIPTION :
  --   Returns true if there is an index i in v such that v(i) = k.

  begin
    for i in v'range loop
      if v(i) = k
       then return true;
      end if;
    end loop;
    return false;
  end Is_In;

  function Sort ( f : Vector ) return Vector is

  -- DESCRIPTION :
  --   Sorts the vector in increasing order.

    res : Vector(f'range) := f;
    tmp : natural32;

  begin
    for i in res'range loop
      for j in i+1..res'last loop
        if res(i) > res(j) then
          tmp := res(i);
          res(i) := res(j);
          res(j) := tmp;
        end if;
      end loop;
    end loop;
    return res;
  end Sort;

  function Join ( f : Vector; k : natural32 ) return Vector is

  -- DESCRIPTION :
  --   Returns the vector f, with k added to it.
	  
  begin
    if Is_In(f,k) then
      return f;
    else
      declare
        fk : Vector(f'first..f'last+1);
      begin
        fk(fk'first) := k;
        for i in f'range loop
          fk(i+1) := f(i);
        end loop;
        return Sort(fk);
      end;
    end if;
  end Join;

  function Is_Disjoint ( f1,f2 : Vector ) return boolean is
 
  -- DESCRIPTION :
  --   Returns true if the two vectors have no element in common.

  begin
    for i in f2'range loop
      if Is_In(f1,f2(i))
       then return false;
      end if;
    end loop;
    return true;
  end Is_Disjoint;

  function Is_Disjoint ( factors : VecVec; f : Vector ) return boolean is

  -- DESCRIPTION :
  --   Returns true if the factor f contains no point in common with
  --   the factors in the vector of vectors.

  -- REQUIRED :
  --   factors contains no empty vectors and 
  --   the range of f contains no undefined entries.

  begin
    for i in factors'range loop
      if not Is_Disjoint(factors(i).all,f)
       then return false;
      end if;
    end loop;
    return true;
  end Is_Disjoint;

  procedure Next_Disjoint_Factor
               ( n : in natural32; factors : in VecVec; f : in Vector;
                 new_fac : out Vector; len : out integer32 ) is

  -- DESCRIPTION :
  --   Finds the next disjoint factor, following f.

  -- ON ENTRY :
  --   n         number of witness points;
  --   factors   the factors already stored;
  --   f         current factor;

  -- ON RETURN :
  --   new_fac   factor following f, disjoint from factors;
  --   len       valid range for new_fac,
  --             if 0, then at end of enumeration.

    fac : Vector(1..integer32(n));
    nb : integer32 := f'length;

  begin
    fac(1..nb) := f;
    loop
      declare
        nf : constant Vector := Next_Factor(n,fac(1..nb));
      begin
        if nf = fac(1..nb) then
          len := 0; return;
        elsif Is_Disjoint(factors,nf) then
          new_fac(nf'range) := nf; len := nf'length; return;   
        else
          nb := nf'length; fac(1..nb) := nf;
        end if;
      end;
    end loop;
  end Next_Disjoint_Factor;

  generic

    with procedure Output_Factors
                      ( fs : in VecVec; count,depth : in natural32 );

  procedure Enumerate ( n : in natural32; factors : in VecVec;
                        f : in Vector; cnt : in out natural32;
                        depth : in natural32 );

  -- DESCRIPTION :
  --   Enumerates all factorization of a degree n component.

  -- ON ENTRY :
  --   n         number of witness points;
  --   factors   current factors collected;
  --   f         candidate factor, to accept or to reject;
  --   cnt       counts number of factorizations;
  --   depth     depth of the enumeration tree.

  -- ON RETURN :
  --   cnt       updated counter for the number of factorizations.

  procedure Enumerate ( n : in natural32; factors : in VecVec;
                        f : in Vector; cnt : in out natural32;
                        depth : in natural32 ) is

    new_factors : VecVec(factors'first..factors'last+1);
    one : constant Vector(1..1) := (1..1 => 1);
    new_f : Vector(1..integer32(n));
    len : integer32;

  begin
    put("Accepting factor :"); put(f); new_line;
    new_factors(factors'range) := factors;
    new_factors(factors'last+1) := new Vector'(f);
    if Is_Disjoint(new_factors,one)
     then new_f(1) := 1; len := 1;
     else Next_Disjoint_Factor(n,new_factors,one,new_f,len);
    end if;
    if len = 0 then
      cnt := cnt + 1;
      Output_Factors(new_factors,cnt,depth);
    else
      Enumerate(n,new_factors,new_f(1..len),cnt,depth+1);
    end if;
    put("Rejecting factor :"); put(f); new_line;
    Next_Disjoint_Factor(n,factors,f,new_f,len);
    if len > 0 then
      if not Is_In(new_f(1..len),f(1)) then
        declare
          jf : constant Vector := Join(new_f(1..len),f(1));
        begin
          Enumerate(n,factors,jf,cnt,depth+1);
        end;
      else
        Enumerate(n,factors,new_f(1..len),cnt,depth+1);
      end if;
    end if;
  end Enumerate;

  generic

    with function Is_Factor ( f : Vector ) return boolean;

  function Search ( n : natural32; factors : VecVec;
                    f : Vector; depth : natural32 ) return VecVec;

  function Search ( n : natural32; factors : VecVec;
                    f : Vector; depth : natural32 ) return VecVec is

    new_factors : VecVec(factors'first..factors'last+1);
    one : constant Vector(1..1) := (1..1 => 1);
    new_f : Vector(1..integer32(n));
    len : integer32;

  begin
    if Is_Factor(f) then
      new_factors(factors'range) := factors;
      new_factors(factors'last+1) := new Vector'(f);
      if Is_Disjoint(new_factors,one)
       then new_f(1) := 1; len := 1;
       else Next_Disjoint_Factor(n,new_factors,one,new_f,len);
      end if;
      if len = 0
       then return new_factors;
       else return Search(n,new_factors,new_f(1..len),depth+1);
      end if;
    else
      Next_Disjoint_Factor(n,factors,f,new_f,len);
      if len > 0 then
        if not Is_In(new_f(1..len),f(1)) then
          declare
            jf : constant Vector := Join(new_f(1..len),f(1));
          begin
            return Search(n,factors,jf,depth+1);
          end;
        else
          return Search(n,factors,new_f(1..len),depth+1);
        end if;
      else
        return factors;
      end if;
    end if;
  end Search;

  generic

    with function Is_Factor ( f : Vector ) return boolean;

  function Search_with_Output
	          ( file : file_type; n : natural32; factors : VecVec;
                    f : Vector; depth : natural32 ) return VecVec;

  function Search_with_Output
	          ( file : file_type; n : natural32; factors : VecVec;
                    f : Vector; depth : natural32 ) return VecVec is

    new_factors : VecVec(factors'first..factors'last+1);
    one : constant Vector(1..1) := (1..1 => 1);
    new_f : Vector(1..integer32(n));
    len : integer32;

  begin
    if Is_Factor(f) then
      put(file,"  Accepted factor :"); put(file,f); new_line(file);
      new_factors(factors'range) := factors;
      new_factors(factors'last+1) := new Vector'(f);
      if Is_Disjoint(new_factors,one)
       then new_f(1) := 1; len := 1;
       else Next_Disjoint_Factor(n,new_factors,one,new_f,len);
      end if;
      if len = 0
       then return new_factors;
       else return Search_with_Output
                     (file,n,new_factors,new_f(1..len),depth+1);
      end if;
    else
      put(file,"  Rejected factor :"); put(file,f); new_line(file);
      Next_Disjoint_Factor(n,factors,f,new_f,len);
      if len > 0 then
        if not Is_In(new_f(1..len),f(1)) then
          declare
            jf : constant Vector := Join(new_f(1..len),f(1));
          begin
            return Search_with_Output(file,n,factors,jf,depth+1);
          end;
        else
          return Search_with_Output(file,n,factors,new_f(1..len),depth+1);
        end if;
      else
        return factors;
      end if;
    end if;
  end Search_with_Output;

-- TARGET ROUTINES :

  procedure Enumerate_Factors ( n : in natural32 ) is

    nb : integer32 := 1;
    f : Vector(1..integer32(n));

  begin
    f(1) := 1;
    put(f(1..nb)); new_line;
    loop
      declare
        nf : constant Vector := Next_Factor(n,f(1..nb));
      begin
        exit when nf = f(1..nb);
        nb := nf'length;
        f(1..nb) := nf;
        put(f(1..nb)); new_line;
      end;
    end loop;
  end Enumerate_Factors;

  procedure Enumerate_Factorizations ( n : in natural32 ) is

    factors : VecVec(1..0);
    f : Vector(1..1);
    cnt : natural32 := 0;

    procedure Enum is new Enumerate(Output_Factors);

  begin
    f(1) := 1;
    Enum(n,factors,f,cnt,0);
  end Enumerate_Factorizations;

  function Search_Factorization ( n : natural32 ) return VecVec is

    factors : VecVec(1..0);
    f : Vector(1..1);

    function Seek is new Search(Is_Factor);

  begin
    f(1) := 1;
    return Seek(n,factors,f,0);
  end Search_Factorization;

  function Search_Factorization_with_Output
             ( file : file_type; n : natural32 ) return VecVec is

    factors : VecVec(1..0);
    f : Vector(1..1);

    function Seek is new Search_with_Output(Is_Factor);

  begin
    f(1) := 1;
    return Seek(file,n,factors,f,0);
  end Search_Factorization_with_Output;

-- TARGET FACTORIZATION FUNCTIONS :

  function Factor ( n : natural32;
                    grid : Array_of_Standard_Sample_Lists ) return VecVec is

    tol : constant double_float := 1.0E-8;

    function Certify_Factor
               ( f : Standard_Natural_Vectors.Vector ) return boolean is

      use Standard_Complex_Numbers;

      sub_grid : constant Array_of_Standard_Sample_Lists(grid'range)
               := Select_Sub_Grid(grid,f);
      q : constant Standard_Complex_Vectors.Vector := Create(sub_grid,1);
      val,eva : Complex_Number;

    begin
      Eval_Trace(q,1,sub_grid(2),val,eva);
      return (Standard_Complex_Numbers.AbsVal(val-eva) < tol);
    end Certify_Factor;

    function CombFac is new Search_Factorization(Certify_Factor);

  begin
    declare
      factors : constant Standard_Natural_VecVecs.VecVec := CombFac(n);
    begin
      return factors;
    end;
  end Factor;

  function Factor ( file : file_type; n : natural32;
                    grid : Array_of_Standard_Sample_Lists ) return VecVec is

    timer : Timing_Widget;
    tol : constant double_float := 1.0E-8;

    function Certify_Factor
               ( f : Standard_Natural_Vectors.Vector ) return boolean is

      use Standard_Complex_Numbers;

      sub_grid : constant Array_of_Standard_Sample_Lists(grid'range)
               := Select_Sub_Grid(grid,f);
      q : constant Standard_Complex_Vectors.Vector := Create(sub_grid,1);
      val,eva : Complex_Number;

    begin
      Eval_Trace(q,1,sub_grid(2),val,eva);
      put(file,"calculated sum at samples : ");
      put(file,val); new_line(file);
      put(file,"value at the linear trace : ");
      put(file,eva); new_line(file);
      put(file,"The witness points"); put(file,f);
      if Standard_Complex_Numbers.AbsVal(val-eva) < tol
       then put_line(file," define a factor.");        return true;
       else put_line(file," do not define a factor."); return false;
      end if;
    end Certify_Factor;

    function CombFac is new Search_Factorization_with_Output(Certify_Factor);

  begin
    tstart(timer);
    declare
      factors : constant Standard_Natural_VecVecs.VecVec
              := CombFac(file,n);
    begin
      tstop(timer);
      put(file,"Factorization found : ");
      Write(file,factors); new_line(file);
      new_line(file);
      put_line(file,"The decomposition : ");
      Write_Labels(file,factors);
      new_line(file);
      print_times(file,timer,"combinatorial factorization");
      return factors;
    end;
  end Factor;

  function Factor ( n : natural32;
                    grid : Array_of_DoblDobl_Sample_Lists ) return VecVec is

    tol : constant double_float := 1.0E-8;

    function Certify_Factor
               ( f : Standard_Natural_Vectors.Vector ) return boolean is

      use DoblDobl_Complex_Numbers;

      sub_grid : constant Array_of_DoblDobl_Sample_Lists(grid'range)
               := Select_Sub_Grid(grid,f);
      q : constant DoblDobl_Complex_Vectors.Vector := Create(sub_grid,1);
      val,eva : Complex_Number;

    begin
      Eval_Trace(q,1,sub_grid(2),val,eva);
      return (AbsVal(val-eva) < tol);
    end Certify_Factor;

    function CombFac is new Search_Factorization(Certify_Factor);

  begin
    declare
      factors : constant Standard_Natural_VecVecs.VecVec := CombFac(n);
    begin
      return factors;
    end;
  end Factor;

  function Factor ( file : file_type; n : natural32;
                    grid : Array_of_DoblDobl_Sample_Lists ) return VecVec is

    timer : Timing_Widget;
    tol : constant double_float := 1.0E-8;

    function Certify_Factor
               ( f : Standard_Natural_Vectors.Vector ) return boolean is

      use DoblDobl_Complex_Numbers;

      sub_grid : constant Array_of_DoblDobl_Sample_Lists(grid'range)
               := Select_Sub_Grid(grid,f);
      q : constant DoblDobl_Complex_Vectors.Vector := Create(sub_grid,1);
      val,eva : Complex_Number;

    begin
      Eval_Trace(q,1,sub_grid(2),val,eva);
      put(file,"calculated sum at samples : ");
      put(file,val); new_line(file);
      put(file,"value at the linear trace : ");
      put(file,eva); new_line(file);
      put(file,"The witness points"); put(file,f);
      if AbsVal(val-eva) < tol
       then put_line(file," define a factor.");        return true;
       else put_line(file," do not define a factor."); return false;
      end if;
    end Certify_Factor;

    function CombFac is new Search_Factorization_with_Output(Certify_Factor);

  begin
    tstart(timer);
    declare
      factors : constant Standard_Natural_VecVecs.VecVec := CombFac(file,n);
    begin
      tstop(timer);
      put(file,"Factorization found : ");
      Write(file,factors); new_line(file);
      new_line(file);
      put_line(file,"The decomposition : ");
      Write_Labels(file,factors);
      new_line(file);
      print_times(file,timer,"combinatorial factorization");
      return factors;
    end;
  end Factor;

  function Factor ( n : natural32;
                    grid : Array_of_QuadDobl_Sample_Lists ) return VecVec is

    tol : constant double_float := 1.0E-8;

    function Certify_Factor
               ( f : Standard_Natural_Vectors.Vector ) return boolean is

      use QuadDobl_Complex_Numbers;

      sub_grid : constant Array_of_QuadDobl_Sample_Lists(grid'range)
               := Select_Sub_Grid(grid,f);
      q : constant QuadDobl_Complex_Vectors.Vector := Create(sub_grid,1);
      val,eva : Complex_Number;

    begin
      Eval_Trace(q,1,sub_grid(2),val,eva);
      return (AbsVal(val-eva) < tol);
    end Certify_Factor;

    function CombFac is new Search_Factorization(Certify_Factor);

  begin
    declare
      factors : constant Standard_Natural_VecVecs.VecVec := CombFac(n);
    begin
      return factors;
    end;
  end Factor;

  function Factor ( file : file_type; n : natural32;
                    grid : Array_of_QuadDobl_Sample_Lists ) return VecVec is

    timer : Timing_Widget;
    tol : constant double_float := 1.0E-8;

    function Certify_Factor
               ( f : Standard_Natural_Vectors.Vector ) return boolean is

      use QuadDobl_Complex_Numbers;

      sub_grid : constant Array_of_QuadDobl_Sample_Lists(grid'range)
               := Select_Sub_Grid(grid,f);
      q : constant QuadDobl_Complex_Vectors.Vector := Create(sub_grid,1);
      val,eva : Complex_Number;

    begin
      Eval_Trace(q,1,sub_grid(2),val,eva);
      put(file,"calculated sum at samples : ");
      put(file,val); new_line(file);
      put(file,"value at the linear trace : ");
      put(file,eva); new_line(file);
      put(file,"The witness points"); put(file,f);
      if AbsVal(val-eva) < tol
       then put_line(file," define a factor.");        return true;
       else put_line(file," do not define a factor."); return false;
      end if;
    end Certify_Factor;

    function CombFac is new Search_Factorization_with_Output(Certify_Factor);

  begin
    tstart(timer);
    declare
      factors : constant Standard_Natural_VecVecs.VecVec := CombFac(file,n);
    begin
      tstop(timer);
      put(file,"Factorization found : ");
      Write(file,factors); new_line(file);
      new_line(file);
      put_line(file,"The decomposition : ");
      Write_Labels(file,factors);
      new_line(file);
      print_times(file,timer,"combinatorial factorization");
      return factors;
    end;
  end Factor;

end Combinatorial_Factorization;
