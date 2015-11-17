with Communications_with_User;          use Communications_with_User;
with Numbers_io;                        use Numbers_io;
with Standard_Natural_Numbers_io;       use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Double_Double_Numbers;             use Double_Double_Numbers;
with Quad_Double_Numbers;               use Quad_Double_Numbers;
with Standard_Complex_Numbers;
with DoblDobl_Complex_Numbers;
with QuadDobl_Complex_Numbers;
with Standard_Integer_Vectors;          use Standard_Integer_Vectors;
with Standard_Integer_Vectors_io;       use Standard_Integer_Vectors_io;
with Standard_Complex_Vectors;
with Standard_Complex_Vectors_io;       use Standard_Complex_Vectors_io;
with Standard_Complex_VecVecs;          use Standard_Complex_VecVecs;
with DoblDobl_Complex_Vectors;
with DoblDobl_Complex_Vectors_io;       use DoblDobl_Complex_Vectors_io;
with QuadDobl_Complex_Vectors;
with QuadDobl_Complex_Vectors_io;       use QuadDobl_Complex_Vectors_io;
with Standard_Complex_Polynomials;
with Standard_Complex_Poly_SysFun;
with Standard_Complex_Poly_Systems_io;  use Standard_Complex_Poly_Systems_io;
with DoblDobl_Complex_Polynomials;
with DoblDobl_Complex_Poly_SysFun;
with DoblDobl_Complex_Poly_Systems_io;  use DoblDobl_Complex_Poly_Systems_io;
with QuadDobl_Complex_Polynomials;
with QuadDobl_Complex_Poly_SysFun;
with QuadDobl_Complex_Poly_Systems_io;  use QuadDobl_Complex_Poly_Systems_io;
with Standard_Complex_Solutions_io;     use Standard_Complex_Solutions_io;
with DoblDobl_Complex_Solutions_io;     use DoblDobl_Complex_Solutions_io;
with QuadDobl_Complex_Solutions_io;     use QuadDobl_Complex_Solutions_io;
with Permutations,Permute_Operations;   use Permutations,Permute_Operations;
with Witness_Sets,Witness_Sets_io;      use Witness_Sets,Witness_Sets_io;
with Standard_Diagonal_Polynomials;     use Standard_Diagonal_Polynomials;
with Standard_Diagonal_Solutions;       use Standard_Diagonal_Solutions;
with DoblDobl_Diagonal_Polynomials;     use DoblDobl_Diagonal_Polynomials;
with DoblDobl_Diagonal_Solutions;       use DoblDobl_Diagonal_Solutions;
with QuadDobl_Diagonal_Polynomials;     use QuadDobl_Diagonal_Polynomials;
with QuadDobl_Diagonal_Solutions;       use QuadDobl_Diagonal_Solutions;
with Extrinsic_Diagonal_Homotopies;     use Extrinsic_Diagonal_Homotopies;
with Extrinsic_Diagonal_Homotopies_io;  use Extrinsic_Diagonal_Homotopies_io;

package body Extrinsic_Diagonal_Solvers is

-- UTILITIES TO SET UP THE DIAGONAL HOMOTOPIES :

  function Is_Dummy ( p : Standard_Complex_Polynomials.Poly;
                      k : integer32 ) return boolean is

  -- DESCRIPTION :
  --   Returns true if the polynomial has only one term that is linear
  --   and where the occurring variables are among the k last ones.

    use Standard_Complex_Polynomials;
    n : constant integer32 := integer32(Number_of_Unknowns(p));

  begin
    if Number_of_Terms(p) > 1 then
      return false;
    elsif Degree(p) > 1 then
      return false;
    else
      for i in 1..n-k loop
        if Degree(p,i) = 1
         then return false;
        end if;
      end loop;
      return true;
    end if;
  end Is_Dummy;

  function Is_Dummy ( p : DoblDobl_Complex_Polynomials.Poly;
                      k : integer32 ) return boolean is

  -- DESCRIPTION :
  --   Returns true if the polynomial has only one term that is linear
  --   and where the occurring variables are among the k last ones.

    use DoblDobl_Complex_Polynomials;
    n : constant integer32 := integer32(Number_of_Unknowns(p));

  begin
    if Number_of_Terms(p) > 1 then
      return false;
    elsif Degree(p) > 1 then
      return false;
    else
      for i in 1..n-k loop
        if Degree(p,i) = 1
         then return false;
        end if;
      end loop;
      return true;
    end if;
  end Is_Dummy;

  function Is_Dummy ( p : QuadDobl_Complex_Polynomials.Poly;
                      k : integer32 ) return boolean is

  -- DESCRIPTION :
  --   Returns true if the polynomial has only one term that is linear
  --   and where the occurring variables are among the k last ones.

    use QuadDobl_Complex_Polynomials;
    n : constant integer32 := integer32(Number_of_Unknowns(p));

  begin
    if Number_of_Terms(p) > 1 then
      return false;
    elsif Degree(p) > 1 then
      return false;
    else
      for i in 1..n-k loop
        if Degree(p,i) = 1
         then return false;
        end if;
      end loop;
      return true;
    end if;
  end Is_Dummy;

  function Number_of_Dummies
             ( p : Standard_Complex_Poly_Systems.Poly_Sys;
               k : integer32 ) return natural32 is

  -- DESCRIPTION :
  --   Returns the number of equations that turn the slack variables
  --   into dummies.

    res : natural32 := 0;

  begin
    for i in p'range loop
      if Is_Dummy(p(i),k)
       then res := res + 1;
      end if;
    end loop;
    return res;
  end Number_of_Dummies;

  function Number_of_Dummies
             ( p : DoblDobl_Complex_Poly_Systems.Poly_Sys;
               k : integer32 ) return natural32 is

  -- DESCRIPTION :
  --   Returns the number of equations that turn the slack variables
  --   into dummies.

    res : natural32 := 0;

  begin
    for i in p'range loop
      if Is_Dummy(p(i),k)
       then res := res + 1;
      end if;
    end loop;
    return res;
  end Number_of_Dummies;

  function Number_of_Dummies
             ( p : QuadDobl_Complex_Poly_Systems.Poly_Sys;
               k : integer32 ) return natural32 is

  -- DESCRIPTION :
  --   Returns the number of equations that turn the slack variables
  --   into dummies.

    res : natural32 := 0;

  begin
    for i in p'range loop
      if Is_Dummy(p(i),k)
       then res := res + 1;
      end if;
    end loop;
    return res;
  end Number_of_Dummies;
 
-- UTILITIES TO SAVE THE RESULTS :

  procedure Test_Solutions
              ( file : in file_type;
                p : in Standard_Complex_Poly_Systems.Poly_Sys;
                s : in Standard_Complex_Solutions.Solution_List ) is

    use Standard_Complex_Solutions;

    y : Standard_Complex_Vectors.Vector(p'range);
    t : Solution_List := s;
    ls : Link_to_Solution;
    ind : natural32 := 0;

  begin
    while not Is_Null(t) loop
      ind := ind + 1;
      ls := Head_Of(t);
      y := Standard_Complex_Poly_SysFun.Eval(p,ls.v);
      put(file,"Value at solution : "); put(file,ind,1);
      put_line(file," :"); put_line(file,y);
      t := Tail_Of(t);
    end loop;
  end Test_Solutions;

  procedure Test_Solutions
              ( file : in file_type;
                p : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                s : in DoblDobl_Complex_Solutions.Solution_List ) is

    use DoblDobl_Complex_Solutions;

    y : DoblDobl_Complex_Vectors.Vector(p'range);
    t : Solution_List := s;
    ls : Link_to_Solution;
    ind : natural32 := 0;

  begin
    while not Is_Null(t) loop
      ind := ind + 1;
      ls := Head_Of(t);
      y := DoblDobl_Complex_Poly_SysFun.Eval(p,ls.v);
      put(file,"Value at solution : "); put(file,ind,1);
      put_line(file," :"); put_line(file,y);
      t := Tail_Of(t);
    end loop;
  end Test_Solutions;

  procedure Test_Solutions
              ( file : in file_type;
                p : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                s : in QuadDobl_Complex_Solutions.Solution_List ) is

    use QuadDobl_Complex_Solutions;

    y : QuadDobl_Complex_Vectors.Vector(p'range);
    t : Solution_List := s;
    ls : Link_to_Solution;
    ind : natural32 := 0;

  begin
    while not Is_Null(t) loop
      ind := ind + 1;
      ls := Head_Of(t);
      y := QuadDobl_Complex_Poly_SysFun.Eval(p,ls.v);
      put(file,"Value at solution : "); put(file,ind,1);
      put_line(file," :"); put_line(file,y);
      t := Tail_Of(t);
    end loop;
  end Test_Solutions;

  procedure Save_Start_System
              ( p : in Standard_Complex_Poly_Systems.Poly_Sys;
                s : in Standard_Complex_Solutions.Solution_List ) is

    use Standard_Complex_Solutions;

    file : file_type;

  begin
    new_line;
    put_line("Reading the name of the file to save the start system.");
    Read_Name_and_Create_File(file);
    put_line(file,p);
    new_line(file);
    put_line(file,"THE SOLUTIONS :");
    put(file,Length_Of(s),natural32(p'last),s);
  end Save_Start_System;

  procedure Save_Start_System
              ( p : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                s : in DoblDobl_Complex_Solutions.Solution_List ) is

    use DoblDobl_Complex_Solutions;

    file : file_type;

  begin
    new_line;
    put_line("Reading the name of the file to save the start system.");
    Read_Name_and_Create_File(file);
    put_line(file,p);
    new_line(file);
    put_line(file,"THE SOLUTIONS :");
    put(file,Length_Of(s),natural32(p'last),s);
  end Save_Start_System;

  procedure Save_Start_System
              ( p : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                s : in QuadDobl_Complex_Solutions.Solution_List ) is

    use QuadDobl_Complex_Solutions;

    file : file_type;

  begin
    new_line;
    put_line("Reading the name of the file to save the start system.");
    Read_Name_and_Create_File(file);
    put_line(file,p);
    new_line(file);
    put_line(file,"THE SOLUTIONS :");
    put(file,Length_Of(s),natural32(p'last),s);
  end Save_Start_System;

  procedure Save_Target_System
              ( p : in Standard_Complex_Poly_Systems.Poly_Sys ) is

    file : file_type;

  begin
    new_line;
    put_line("Reading the name of the file to save the target system.");
    Read_Name_and_Create_File(file);
    put_line(file,p);
  end Save_Target_System;

  procedure Save_Target_System
              ( p : in DoblDobl_Complex_Poly_Systems.Poly_Sys ) is

    file : file_type;

  begin
    new_line;
    put_line("Reading the name of the file to save the target system.");
    Read_Name_and_Create_File(file);
    put_line(file,p);
  end Save_Target_System;

  procedure Save_Target_System
              ( p : in QuadDobl_Complex_Poly_Systems.Poly_Sys ) is

    file : file_type;

  begin
    new_line;
    put_line("Reading the name of the file to save the target system.");
    Read_Name_and_Create_File(file);
    put_line(file,p);
  end Save_Target_System;

-- MAIN DRIVERS :

  procedure Standard_Randomize_System is

    use Standard_Complex_Polynomials;
    use Standard_Complex_Poly_Systems;

    lp : Link_to_Poly_Sys;
    k : natural32 := 0;
 
  begin
    new_line;
    put_line("Reading the polynomial system...");
    get(lp);
    new_line;
    put("Give the dimension of the solution component : ");
    get(k); skip_line;
    declare 
      n : constant natural32 := Number_of_Unknowns(lp(lp'first));
      rp : constant Poly_Sys := Complete(n,k,lp.all);
    begin
      Save_Target_System(rp);
    end;
  end Standard_Randomize_System;

  procedure DoblDobl_Randomize_System is

    use DoblDobl_Complex_Polynomials;
    use DoblDobl_Complex_Poly_Systems;

    lp : Link_to_Poly_Sys;
    k : natural32 := 0;
 
  begin
    new_line;
    put_line("Reading the polynomial system...");
    get(lp);
    new_line;
    put("Give the dimension of the solution component : ");
    get(k); skip_line;
    declare 
      n : constant natural32 := Number_of_Unknowns(lp(lp'first));
      rp : constant Poly_Sys := Complete(n,k,lp.all);
    begin
      Save_Target_System(rp);
    end;
  end DoblDobl_Randomize_System;

  procedure QuadDobl_Randomize_System is

    use QuadDobl_Complex_Polynomials;
    use QuadDobl_Complex_Poly_Systems;

    lp : Link_to_Poly_Sys;
    k : natural32 := 0;
 
  begin
    new_line;
    put_line("Reading the polynomial system...");
    get(lp);
    new_line;
    put("Give the dimension of the solution component : ");
    get(k); skip_line;
    declare 
      n : constant natural32 := Number_of_Unknowns(lp(lp'first));
      rp : constant Poly_Sys := Complete(n,k,lp.all);
    begin
      Save_Target_System(rp);
    end;
  end QuadDobl_Randomize_System;

  procedure Randomize_System is

    p : constant character := Prompt_for_Precision;

  begin
    case p is
      when '0' => Standard_Randomize_System;
      when '1' => DoblDobl_Randomize_System;
      when '2' => QuadDobl_Randomize_System;
      when others => null;
    end case;
  end Randomize_System;

  procedure Build_Cascade_Homotopy
              ( file : in file_type;
                p1e,p2e : in Standard_Complex_Poly_Systems.Poly_Sys;
                dim1,dim2 : in natural32;
                sols1e,sols2e : in Standard_Complex_Solutions.Solution_List;
                s1e,s2e : in Array_of_Symbols ) is

    use Standard_Complex_Numbers;
    use Standard_Complex_Polynomials;
    use Standard_Complex_Poly_Systems;
    use Standard_Complex_Solutions;

    dim : constant natural32 := Cascade_Dimension(p1e,p2e,dim1,dim2);
    k : constant natural32 := Number_of_Unknowns(p1e(p1e'first))-dim1;
    ch_start,ch_target : Poly_Sys(1..integer32(dim));
    sols1 : constant Solution_List := Remove_Embedding(sols1e,dim1);
    sols2 : constant Solution_List := Remove_Embedding(sols2e,dim2);
    sols : constant Solution_List := Product(sols1,sols2);
    embsols : Solution_List;
    s1 : constant Array_of_Symbols := Remove_Embed_Symbols(s1e);
    s11 : constant Array_of_Symbols(s1'range) := Add_Suffix(s1,'1');
    s2 : constant Array_of_Symbols := Remove_Embed_Symbols(s2e);
    s22 : constant Array_of_Symbols(s1'range) := Add_Suffix(s2,'2');
 
  begin
    new_line(file);
    put_line(file,"Building the homotopy to start the cascade...");
    put(file,"k : "); put(file,k,1);
    put(file,"  a : "); put(file,dim1,1);
    put(file,"  b : "); put(file,dim2,1); new_line(file);
    put(file,"cascade dimension : "); put(file,dim,1); new_line(file);
    if dim1+dim2 < k then
      put_line(file,"calling Extrinsic_Diagonal_Homotopies.Cascade1");
      Cascade1(p1e,p2e,dim1,dim2,ch_start,ch_target);
      embsols := Add_Embedding(sols,dim2);
    else
      put_line(file,"calling Extrinsic_Diagonal_Homotopies.Cascade2");
      Cascade2(p1e,p2e,dim1,dim2,ch_start,ch_target);
      embsols := Add_Embedding(sols,k-dim1);
    end if;
    Symbol_Table.Clear;
    Assign_Symbol_Table(s11,s22);
    Add_Embed_Symbols(dim2);
    Save_Target_System(ch_target);
    Set_Continuation_Parameter(embsols,Create(0.0));
    Save_Start_System(ch_start,embsols);
    Test_Solutions(file,ch_start,embsols);
  end Build_Cascade_Homotopy;

  procedure Build_Cascade_Homotopy
              ( file : in file_type;
                p1e,p2e : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                dim1,dim2 : in natural32;
                sols1e,sols2e : in DoblDobl_Complex_Solutions.Solution_List;
                s1e,s2e : in Array_of_Symbols ) is

    use DoblDobl_Complex_Numbers;
    use DoblDobl_Complex_Polynomials;
    use DoblDobl_Complex_Poly_Systems;
    use DoblDobl_Complex_Solutions;

    dim : constant natural32 := Cascade_Dimension(p1e,p2e,dim1,dim2);
    k : constant natural32 := Number_of_Unknowns(p1e(p1e'first))-dim1;
    ch_start,ch_target : Poly_Sys(1..integer32(dim));
    sols1 : constant Solution_List := Remove_Embedding(sols1e,dim1);
    sols2 : constant Solution_List := Remove_Embedding(sols2e,dim2);
    sols : constant Solution_List := Product(sols1,sols2);
    embsols : Solution_List;
    s1 : constant Array_of_Symbols := Remove_Embed_Symbols(s1e);
    s11 : constant Array_of_Symbols(s1'range) := Add_Suffix(s1,'1');
    s2 : constant Array_of_Symbols := Remove_Embed_Symbols(s2e);
    s22 : constant Array_of_Symbols(s1'range) := Add_Suffix(s2,'2');
    zero : constant double_double := create(0.0);
 
  begin
    new_line(file);
    put_line(file,"Building the homotopy to start the cascade...");
    put(file,"k : "); put(file,k,1);
    put(file,"  a : "); put(file,dim1,1);
    put(file,"  b : "); put(file,dim2,1); new_line(file);
    put(file,"cascade dimension : "); put(file,dim,1); new_line(file);
    if dim1+dim2 < k then
      put_line(file,"calling Extrinsic_Diagonal_Homotopies.Cascade1");
      Cascade1(p1e,p2e,dim1,dim2,ch_start,ch_target);
      embsols := Add_Embedding(sols,dim2);
    else
      put_line(file,"calling Extrinsic_Diagonal_Homotopies.Cascade2");
      Cascade2(p1e,p2e,dim1,dim2,ch_start,ch_target);
      embsols := Add_Embedding(sols,k-dim1);
    end if;
    Symbol_Table.Clear;
    Assign_Symbol_Table(s11,s22);
    Add_Embed_Symbols(dim2);
    Save_Target_System(ch_target);
    Set_Continuation_Parameter(embsols,Create(zero));
    Save_Start_System(ch_start,embsols);
    Test_Solutions(file,ch_start,embsols);
  end Build_Cascade_Homotopy;

  procedure Build_Cascade_Homotopy
              ( file : in file_type;
                p1e,p2e : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                dim1,dim2 : in natural32;
                sols1e,sols2e : in QuadDobl_Complex_Solutions.Solution_List;
                s1e,s2e : in Array_of_Symbols ) is

    use QuadDobl_Complex_Numbers;
    use QuadDobl_Complex_Polynomials;
    use QuadDobl_Complex_Poly_Systems;
    use QuadDobl_Complex_Solutions;

    dim : constant natural32 := Cascade_Dimension(p1e,p2e,dim1,dim2);
    k : constant natural32 := Number_of_Unknowns(p1e(p1e'first))-dim1;
    ch_start,ch_target : Poly_Sys(1..integer32(dim));
    sols1 : constant Solution_List := Remove_Embedding(sols1e,dim1);
    sols2 : constant Solution_List := Remove_Embedding(sols2e,dim2);
    sols : constant Solution_List := Product(sols1,sols2);
    embsols : Solution_List;
    s1 : constant Array_of_Symbols := Remove_Embed_Symbols(s1e);
    s11 : constant Array_of_Symbols(s1'range) := Add_Suffix(s1,'1');
    s2 : constant Array_of_Symbols := Remove_Embed_Symbols(s2e);
    s22 : constant Array_of_Symbols(s1'range) := Add_Suffix(s2,'2');
    zero : constant quad_double := create(0.0);
 
  begin
    new_line(file);
    put_line(file,"Building the homotopy to start the cascade...");
    put(file,"k : "); put(file,k,1);
    put(file,"  a : "); put(file,dim1,1);
    put(file,"  b : "); put(file,dim2,1); new_line(file);
    put(file,"cascade dimension : "); put(file,dim,1); new_line(file);
    if dim1+dim2 < k then
      put_line(file,"calling Extrinsic_Diagonal_Homotopies.Cascade1");
      Cascade1(p1e,p2e,dim1,dim2,ch_start,ch_target);
      embsols := Add_Embedding(sols,dim2);
    else
      put_line(file,"calling Extrinsic_Diagonal_Homotopies.Cascade2");
      Cascade2(p1e,p2e,dim1,dim2,ch_start,ch_target);
      embsols := Add_Embedding(sols,k-dim1);
    end if;
    Symbol_Table.Clear;
    Assign_Symbol_Table(s11,s22);
    Add_Embed_Symbols(dim2);
    Save_Target_System(ch_target);
    Set_Continuation_Parameter(embsols,Create(zero));
    Save_Start_System(ch_start,embsols);
    Test_Solutions(file,ch_start,embsols);
  end Build_Cascade_Homotopy;

  procedure Permute
              ( p : in Permutation;
                sols : in out Standard_Complex_Solutions.Solution_List ) is

  -- DESCRIPTION :
  --   Permutes the vectors in the solution list with the permutation p.

    use Standard_Complex_Solutions;

    tmp : Solution_List := sols;

  begin
    while not Is_Null(tmp) loop
      declare
        ls : constant Link_to_Solution := Head_Of(tmp);
      begin
        ls.v := p*ls.v;
        Set_Head(tmp,ls);
      end;
      tmp := Tail_Of(tmp);
    end loop;
  end Permute;

  procedure Permute
              ( p : in Permutation;
                sols : in out DoblDobl_Complex_Solutions.Solution_List ) is

  -- DESCRIPTION :
  --   Permutes the vectors in the solution list with the permutation p.

    use DoblDobl_Complex_Solutions;

    tmp : Solution_List := sols;

  begin
    while not Is_Null(tmp) loop
      declare
        ls : constant Link_to_Solution := Head_Of(tmp);
      begin
        ls.v := p*ls.v;
        Set_Head(tmp,ls);
      end;
      tmp := Tail_Of(tmp);
    end loop;
  end Permute;

  procedure Permute
              ( p : in Permutation;
                sols : in out QuadDobl_Complex_Solutions.Solution_List ) is

  -- DESCRIPTION :
  --   Permutes the vectors in the solution list with the permutation p.

    use QuadDobl_Complex_Solutions;

    tmp : Solution_List := sols;

  begin
    while not Is_Null(tmp) loop
      declare
        ls : constant Link_to_Solution := Head_Of(tmp);
      begin
        ls.v := p*ls.v;
        Set_Head(tmp,ls);
      end;
      tmp := Tail_Of(tmp);
    end loop;
  end Permute;

  procedure Permute_by_Symbols
              ( file : in file_type;
                p2e : in out Standard_Complex_Poly_Systems.Poly_Sys;
                sols2e : in out Standard_Complex_Solutions.Solution_List;
                s1,s2 : in out Array_of_Symbols;
                dim1,dim2 : in natural32 ) is

  -- DESCRIPTION :
  --   While the polynomial systems should have all their variables
  --   in common, they may have occurred at different places.
  --   This procedure permutes the variables in the system if needed.

    os1 : constant Array_of_Symbols := s1(s1'first..s1'last-integer32(dim1));
    os2 : Array_of_Symbols(s2'first..s2'last-integer32(dim2))
        := s2(s2'first..s2'last-integer32(dim2));
    prm : constant Permutation(os1'range) := Match_Symbols(os1,os2);
    full_prm,inv_prm : Permutation(s2'range);

  begin
    put(file,"The symbols from the 1st system : "); Write(file,s1);
    put(file,"The symbols from the 2nd system : "); Write(file,s2);
    put(file,"The original first symbols  : "); Write(file,os1);
    put(file,"The original second symbols : "); Write(file,os2);
    put(file,"The matching permutation : ");
    put(file,Standard_Integer_Vectors.Vector(prm)); new_line(file);
    Permute(prm,os2);
    s2(os2'range) := os2;
    put(file,"The permuted second symbols : "); Write(file,os2);
    full_prm(prm'range) := prm;
    for i in prm'last+1..full_prm'last loop
      full_prm(i) := i;
    end loop;
    put(file,"The full matching permutation : ");
    put(file,Standard_Integer_Vectors.Vector(full_prm)); new_line(file);
    inv_prm := inv(full_prm);
    put(file,"Inverse of full permutation   : ");
    put(file,Standard_Integer_Vectors.Vector(inv_prm)); new_line(file);
    p2e := p2e*inv_prm;
    Assign_Symbol_Table(s2);
    put(file,"The permuted second system : "); put(file,p2e);
    Permute(inv_prm,sols2e);
  end Permute_by_Symbols;

  procedure Permute_by_Symbols
              ( file : in file_type;
                p2e : in out DoblDobl_Complex_Poly_Systems.Poly_Sys;
                sols2e : in out DoblDobl_Complex_Solutions.Solution_List;
                s1,s2 : in out Array_of_Symbols;
                dim1,dim2 : in natural32 ) is

  -- DESCRIPTION :
  --   While the polynomial systems should have all their variables
  --   in common, they may have occurred at different places.
  --   This procedure permutes the variables in the system if needed.

    os1 : constant Array_of_Symbols := s1(s1'first..s1'last-integer32(dim1));
    os2 : Array_of_Symbols(s2'first..s2'last-integer32(dim2))
        := s2(s2'first..s2'last-integer32(dim2));
    prm : constant Permutation(os1'range) := Match_Symbols(os1,os2);
    full_prm,inv_prm : Permutation(s2'range);

  begin
    put(file,"The symbols from the 1st system : "); Write(file,s1);
    put(file,"The symbols from the 2nd system : "); Write(file,s2);
    put(file,"The original first symbols  : "); Write(file,os1);
    put(file,"The original second symbols : "); Write(file,os2);
    put(file,"The matching permutation : ");
    put(file,Standard_Integer_Vectors.Vector(prm)); new_line(file);
    Permute(prm,os2);
    s2(os2'range) := os2;
    put(file,"The permuted second symbols : "); Write(file,os2);
    full_prm(prm'range) := prm;
    for i in prm'last+1..full_prm'last loop
      full_prm(i) := i;
    end loop;
    put(file,"The full matching permutation : ");
    put(file,Standard_Integer_Vectors.Vector(full_prm)); new_line(file);
    inv_prm := inv(full_prm);
    put(file,"Inverse of full permutation   : ");
    put(file,Standard_Integer_Vectors.Vector(inv_prm)); new_line(file);
    p2e := p2e*inv_prm;
    Assign_Symbol_Table(s2);
    put(file,"The permuted second system : "); put(file,p2e);
    Permute(inv_prm,sols2e);
  end Permute_by_Symbols;

  procedure Permute_by_Symbols
              ( file : in file_type;
                p2e : in out QuadDobl_Complex_Poly_Systems.Poly_Sys;
                sols2e : in out QuadDobl_Complex_Solutions.Solution_List;
                s1,s2 : in out Array_of_Symbols;
                dim1,dim2 : in natural32 ) is

  -- DESCRIPTION :
  --   While the polynomial systems should have all their variables
  --   in common, they may have occurred at different places.
  --   This procedure permutes the variables in the system if needed.

    os1 : constant Array_of_Symbols := s1(s1'first..s1'last-integer32(dim1));
    os2 : Array_of_Symbols(s2'first..s2'last-integer32(dim2))
        := s2(s2'first..s2'last-integer32(dim2));
    prm : constant Permutation(os1'range) := Match_Symbols(os1,os2);
    full_prm,inv_prm : Permutation(s2'range);

  begin
    put(file,"The symbols from the 1st system : "); Write(file,s1);
    put(file,"The symbols from the 2nd system : "); Write(file,s2);
    put(file,"The original first symbols  : "); Write(file,os1);
    put(file,"The original second symbols : "); Write(file,os2);
    put(file,"The matching permutation : ");
    put(file,Standard_Integer_Vectors.Vector(prm)); new_line(file);
    Permute(prm,os2);
    s2(os2'range) := os2;
    put(file,"The permuted second symbols : "); Write(file,os2);
    full_prm(prm'range) := prm;
    for i in prm'last+1..full_prm'last loop
      full_prm(i) := i;
    end loop;
    put(file,"The full matching permutation : ");
    put(file,Standard_Integer_Vectors.Vector(full_prm)); new_line(file);
    inv_prm := inv(full_prm);
    put(file,"Inverse of full permutation   : ");
    put(file,Standard_Integer_Vectors.Vector(inv_prm)); new_line(file);
    p2e := p2e*inv_prm;
    Assign_Symbol_Table(s2);
    put(file,"The permuted second system : "); put(file,p2e);
    Permute(inv_prm,sols2e);
  end Permute_by_Symbols;

  procedure Standard_Diagonal_Cascade is

    use Standard_Complex_Poly_Systems;
    use Standard_Complex_Solutions;
 
    file : file_type;
    lp1,lp2 : Link_to_Poly_Sys;
    sols1e,sols2e : Solution_List;
    dim1,dum1,dim2,dum2 : natural32;
    lsym1,lsym2 : Link_to_Array_of_Symbols;

  begin
    new_line;
    put_line("Reading the first embedded polynomial system...");
    Standard_Read_Embedding(lp1,sols1e,dim1);
    dum1 := Number_of_Dummies(lp1.all,integer32(dim1));
    lsym1 := Get_Link_to_Symbols;
    new_line;
    put_line("Reading the second embedded polynomial system...");
    Standard_Read_Embedding(lp2,sols2e,dim2);
    dum2 := Number_of_Dummies(lp2.all,integer32(dim2));
    if dim1 < dim2 then
      new_line;
      put("WARNING:");
      put(" dim(component-1) = "); put(dim1,1);
      put(", dim(component-2) = "); put(dim2,1);
      put(", and "); put(dim1,1); put("<"); put(dim2,1); put_line(".");
      put_line("This will cause too many variables in the homotopies.");
      put_line("Be advised to stop program and switch order of systems.");
    end if;
    new_line;
    put_line("Reading the name of the output file...");
    Read_Name_and_Create_File(file);
    put_line(file,"The first polynomial system :");
    put(file,lp1.all,lsym1.all);
    put(file,"Dimension of component 1  : ");
    put(file,dim1,1); new_line(file);
    put(file,"Number of dummy slackvars : ");
    put(file,dum1,1); new_line(file);
    put_line(file,"The second polynomial system :"); put(file,lp2.all);
    put(file,"Dimension of component 2  : ");
    put(file,dim2,1); new_line(file);
    put(file,"Number of dummy slackvars : ");
    put(file,dum2,1); new_line(file);
    lsym2 := Get_Link_to_Symbols;
    put(file,"2nd Symbols :"); Write(file,lsym2.all);
    Permute_by_Symbols(file,lp2.all,sols2e,lsym1.all,lsym2.all,dim1,dim2);
    put(file,"1st Symbols :"); Write(file,lsym1.all);
    put(file,"2nd Symbols :"); Write(file,lsym2.all);
   -- Assign_Symbol_Table(lsym2.all);
    Build_Cascade_Homotopy
      (file,lp1.all,lp2.all,dim1,dim2,sols1e,sols2e,lsym1.all,lsym2.all);
  end Standard_Diagonal_Cascade;

  procedure DoblDobl_Diagonal_Cascade is

    use DoblDobl_Complex_Poly_Systems;
    use DoblDobl_Complex_Solutions;
 
    file : file_type;
    lp1,lp2 : Link_to_Poly_Sys;
    sols1e,sols2e : Solution_List;
    dim1,dum1,dim2,dum2 : natural32;
    lsym1,lsym2 : Link_to_Array_of_Symbols;

  begin
    new_line;
    put_line("Reading the first embedded polynomial system...");
    DoblDobl_Read_Embedding(lp1,sols1e,dim1);
    dum1 := Number_of_Dummies(lp1.all,integer32(dim1));
    lsym1 := Get_Link_to_Symbols;
    new_line;
    put_line("Reading the second embedded polynomial system...");
    DoblDobl_Read_Embedding(lp2,sols2e,dim2);
    dum2 := Number_of_Dummies(lp2.all,integer32(dim2));
    if dim1 < dim2 then
      new_line;
      put("WARNING:");
      put(" dim(component-1) = "); put(dim1,1);
      put(", dim(component-2) = "); put(dim2,1);
      put(", and "); put(dim1,1); put("<"); put(dim2,1); put_line(".");
      put_line("This will cause too many variables in the homotopies.");
      put_line("Be advised to stop program and switch order of systems.");
    end if;
    new_line;
    put_line("Reading the name of the output file...");
    Read_Name_and_Create_File(file);
    put_line(file,"The first polynomial system :");
    put(file,lp1.all,lsym1.all);
    put(file,"Dimension of component 1  : ");
    put(file,dim1,1); new_line(file);
    put(file,"Number of dummy slackvars : ");
    put(file,dum1,1); new_line(file);
    put_line(file,"The second polynomial system :"); put(file,lp2.all);
    put(file,"Dimension of component 2  : ");
    put(file,dim2,1); new_line(file);
    put(file,"Number of dummy slackvars : ");
    put(file,dum2,1); new_line(file);
    lsym2 := Get_Link_to_Symbols;
    put(file,"2nd Symbols :"); Write(file,lsym2.all);
    Permute_by_Symbols(file,lp2.all,sols2e,lsym1.all,lsym2.all,dim1,dim2);
    put(file,"1st Symbols :"); Write(file,lsym1.all);
    put(file,"2nd Symbols :"); Write(file,lsym2.all);
   -- Assign_Symbol_Table(lsym2.all);
    Build_Cascade_Homotopy
      (file,lp1.all,lp2.all,dim1,dim2,sols1e,sols2e,lsym1.all,lsym2.all);
  end DoblDobl_Diagonal_Cascade;

  procedure QuadDobl_Diagonal_Cascade is

    use QuadDobl_Complex_Poly_Systems;
    use QuadDobl_Complex_Solutions;
 
    file : file_type;
    lp1,lp2 : Link_to_Poly_Sys;
    sols1e,sols2e : Solution_List;
    dim1,dum1,dim2,dum2 : natural32;
    lsym1,lsym2 : Link_to_Array_of_Symbols;

  begin
    new_line;
    put_line("Reading the first embedded polynomial system...");
    QuadDobl_Read_Embedding(lp1,sols1e,dim1);
    dum1 := Number_of_Dummies(lp1.all,integer32(dim1));
    lsym1 := Get_Link_to_Symbols;
    new_line;
    put_line("Reading the second embedded polynomial system...");
    QuadDobl_Read_Embedding(lp2,sols2e,dim2);
    dum2 := Number_of_Dummies(lp2.all,integer32(dim2));
    if dim1 < dim2 then
      new_line;
      put("WARNING:");
      put(" dim(component-1) = "); put(dim1,1);
      put(", dim(component-2) = "); put(dim2,1);
      put(", and "); put(dim1,1); put("<"); put(dim2,1); put_line(".");
      put_line("This will cause too many variables in the homotopies.");
      put_line("Be advised to stop program and switch order of systems.");
    end if;
    new_line;
    put_line("Reading the name of the output file...");
    Read_Name_and_Create_File(file);
    put_line(file,"The first polynomial system :");
    put(file,lp1.all,lsym1.all);
    put(file,"Dimension of component 1  : ");
    put(file,dim1,1); new_line(file);
    put(file,"Number of dummy slackvars : ");
    put(file,dum1,1); new_line(file);
    put_line(file,"The second polynomial system :"); put(file,lp2.all);
    put(file,"Dimension of component 2  : ");
    put(file,dim2,1); new_line(file);
    put(file,"Number of dummy slackvars : ");
    put(file,dum2,1); new_line(file);
    lsym2 := Get_Link_to_Symbols;
    put(file,"2nd Symbols :"); Write(file,lsym2.all);
    Permute_by_Symbols(file,lp2.all,sols2e,lsym1.all,lsym2.all,dim1,dim2);
    put(file,"1st Symbols :"); Write(file,lsym1.all);
    put(file,"2nd Symbols :"); Write(file,lsym2.all);
   -- Assign_Symbol_Table(lsym2.all);
    Build_Cascade_Homotopy
      (file,lp1.all,lp2.all,dim1,dim2,sols1e,sols2e,lsym1.all,lsym2.all);
  end QuadDobl_Diagonal_Cascade;

  procedure Build_Diagonal_Cascade is

    p : constant character := Prompt_for_Precision;

  begin
    case p is
      when '0' => Standard_Diagonal_Cascade;
      when '1' => DoblDobl_Diagonal_Cascade;
      when '2' => QuadDobl_Diagonal_Cascade;
      when others => null;
    end case;
  end Build_Diagonal_Cascade;

--  procedure Replace_with_Dummies 
--               ( p : in out Poly_Sys; n,dim : in natural ) is
--
--  -- DESCRIPTION :
--  --   Equations which collapsed into constant = 0 equations are replaced
--  --   with equations setting the dummy embedding variables to zero.
--
--  -- REQUIRED : dim > 0.
--
--    nd : natural := 0;
--    dt : Term;
--
--  begin
--    dt.cf := Create(1.0);
--    dt.dg := new Standard_Natural_Vectors.Vector'(1..n+dim => 0);
--    for i in p'range loop
--      if Number_of_Terms(p(i)) <= 1 then
--        nd := nd + 1;
--        dt.dg(n+nd) := 1;
--        Clear(p(i));
--        p(i) := Create(dt);
--      end if;
--      exit when (nd = dim);
--    end loop;
--    Clear(dt);
--  end Replace_with_Dummies;

  procedure Add_Embedding
              ( p : in out Standard_Complex_Polynomials.Poly;
                n,k : in natural32; nbdumb : in out natural32 ) is

  -- DESCRIPTION :
  --   Adds an embedding with k extra variables to a polynomial p
  --   in n variables, where p may be an empty polynomial.
  --   The variable nbdumb keeps track of the empty polynomials.

    use Standard_Complex_Polynomials;

    newcp : Poly;

  begin
    if Number_of_Unknowns(p) < n+k then
      if p = Null_Poly then
        nbdumb := nbdumb + 1;
        if nbdumb <= k
         then p := Add_Dummy(n,k,nbdumb);
        end if;
      else
        newcp := Add_Embedding(p,k);
        Copy(newcp,p); Clear(newcp);
      end if;
    end if;
  end Add_Embedding;

  procedure Add_Embedding
              ( p : in out DoblDobl_Complex_Polynomials.Poly;
                n,k : in natural32; nbdumb : in out natural32 ) is

  -- DESCRIPTION :
  --   Adds an embedding with k extra variables to a polynomial p
  --   in n variables, where p may be an empty polynomial.
  --   The variable nbdumb keeps track of the empty polynomials.

    use DoblDobl_Complex_Polynomials;

    newcp : Poly;

  begin
    if Number_of_Unknowns(p) < n+k then
      if p = Null_Poly then
        nbdumb := nbdumb + 1;
        if nbdumb <= k
         then p := Add_Dummy(n,k,nbdumb);
        end if;
      else
        newcp := Add_Embedding(p,k);
        Copy(newcp,p); Clear(newcp);
      end if;
    end if;
  end Add_Embedding;

  procedure Add_Embedding
              ( p : in out QuadDobl_Complex_Polynomials.Poly;
                n,k : in natural32; nbdumb : in out natural32 ) is

  -- DESCRIPTION :
  --   Adds an embedding with k extra variables to a polynomial p
  --   in n variables, where p may be an empty polynomial.
  --   The variable nbdumb keeps track of the empty polynomials.

    use QuadDobl_Complex_Polynomials;

    newcp : Poly;

  begin
    if Number_of_Unknowns(p) < n+k then
      if p = Null_Poly then
        nbdumb := nbdumb + 1;
        if nbdumb <= k
         then p := Add_Dummy(n,k,nbdumb);
        end if;
      else
        newcp := Add_Embedding(p,k);
        Copy(newcp,p); Clear(newcp);
      end if;
    end if;
  end Add_Embedding;

  procedure Collapse_System
             ( p : in Standard_Complex_Poly_Systems.Poly_Sys;
               sols : in out Standard_Complex_Solutions.Solution_List;
               dim,add2dim : in natural32;
               r : out Standard_Complex_Poly_Systems.Link_to_Poly_Sys ) is

    use Standard_Complex_Polynomials;
    use Standard_Complex_Poly_Systems;
    use Standard_Complex_Solutions;

    nb : constant natural32 := Symbol_Table.Number;
    n : constant integer32 := integer32(nb-dim)/2;
    s1,s2 : Array_of_Symbols(1..n);
   -- ps2 : Array_of_Symbols(1..n);
    p1,p2,prm : Permutation(1..n);
    cbp : Permutation(1..integer32(nb));
    pp,cp : Poly_Sys(p'range);
    res : Poly_Sys(1..n+integer32(add2dim+dim));
    tsols : Solution_List;
    ind : integer32 := res'first;
    cnt : natural32 := 0;
   -- ans : character;

  begin
   -- put("The original ambient dimension (#variables) : ");
   -- put(n,1); new_line;
   -- Write_Symbol_Table;
    Retrieve_Suffixed_Symbols(n,'1',s1,p1);
    Retrieve_Suffixed_Symbols(n,'2',s2,p2);
   -- put("Symbols of 1st system : "); Write(s1);
   -- put("Symbols of 2nd system : "); Write(s2);
   -- put("Indices of 1st symbols : ");
   -- put(Standard_Integer_Vectors.Vector(p1)); new_line;
   -- put("Indices of 2nd symbols : ");
   -- put(Standard_Integer_Vectors.Vector(p2)); new_line;
    prm := Matching_Permutation(s1,s2);
   -- put("The permutation : ");
   -- put(Standard_Integer_Vectors.Vector(prm)); new_line;
   -- ps2 := s2;
   -- for i in s2'range loop
   --   ps2(i) := s2(i);
   -- end loop;
   -- Permute(prm,ps2);
   -- put("The original 1st symbols : "); Write(s1);
   -- put("The permuted 2nd symbols : "); Write(ps2);
    cbp := Combine_Permutations(n,integer32(dim),p1,p2,prm);
   -- put("The combined permutation : ");
   -- put(Standard_Integer_Vectors.Vector(cbp)); new_line;
   -- put_line("The symbols in permuted order : ");
   -- for i in cbp'range loop
   --   declare
   --     sb : Symbol := Symbol_Table.get(cbp(i));
   --   begin
   --     put(" "); Symbol_Table_io.put(sb);
   --     if i = n then new_line; end if;
   --   end;
   -- end loop;
   -- new_line;
   -- put("Continue ? (y/n) "); Ask_Yes_or_No(ans);
    pp := p*cbp;
    Permute(cbp,sols);
   -- put_line("The system before the collapse : "); put(pp);
    cp := Collapse(pp,n);
    for i in s1'range loop
      s1(i) := Remove_Suffix(s1(i));
    end loop;
    Assign_Symbol_Table(s1);
   -- put_line("The collapsed system : "); put(cp);
   -- put("Continue ? (y/n) "); Ask_Yes_or_No(ans);
    Add_Embed_Symbols(dim+add2dim);
    if dim + add2dim > 0 then
      for i in cp'range loop
        Add_Embedding(cp(i),natural32(n),dim+add2dim,cnt);
        if cp(i) /= Null_Poly then
          res(ind) := cp(i);
          ind := ind + 1;
        end if;
      end loop;
    else
      for i in res'range loop
        res(i) := cp(i);
      end loop;
    end if;
   -- put_line("After adding the embedding : "); put(res);
   -- put("Continue ? (y/n) "); Ask_Yes_or_No(ans);
    tsols := Truncate(sols,n);
    Clear(sols);
    if dim + add2dim > 0
     then sols := Add_Embedding(tsols,dim+add2dim); Clear(tsols);
     else sols := tsols;
    end if;
    Clear(pp);
   -- put_line("The collapsed system : "); put(res);
    r := new Poly_Sys'(res);
  end Collapse_System;

  procedure Collapse_System
             ( p : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
               sols : in out DoblDobl_Complex_Solutions.Solution_List;
               dim,add2dim : in natural32;
               r : out DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys ) is

    use DoblDobl_Complex_Polynomials;
    use DoblDobl_Complex_Poly_Systems;
    use DoblDobl_Complex_Solutions;

    nb : constant natural32 := Symbol_Table.Number;
    n : constant integer32 := integer32(nb-dim)/2;
    s1,s2 : Array_of_Symbols(1..n);
    p1,p2,prm : Permutation(1..n);
    cbp : Permutation(1..integer32(nb));
    pp,cp : Poly_Sys(p'range);
    res : Poly_Sys(1..n+integer32(add2dim+dim));
    tsols : Solution_List;
    ind : integer32 := res'first;
    cnt : natural32 := 0;

  begin
    Retrieve_Suffixed_Symbols(n,'1',s1,p1);
    Retrieve_Suffixed_Symbols(n,'2',s2,p2);
    prm := Matching_Permutation(s1,s2);
    cbp := Combine_Permutations(n,integer32(dim),p1,p2,prm);
    pp := p*cbp;
    Permute(cbp,sols);
    cp := Collapse(pp,n);
    for i in s1'range loop
      s1(i) := Remove_Suffix(s1(i));
    end loop;
    Assign_Symbol_Table(s1);
    Add_Embed_Symbols(dim+add2dim);
    if dim + add2dim > 0 then
      for i in cp'range loop
        Add_Embedding(cp(i),natural32(n),dim+add2dim,cnt);
        if cp(i) /= Null_Poly then
          res(ind) := cp(i);
          ind := ind + 1;
        end if;
      end loop;
    else
      for i in res'range loop
        res(i) := cp(i);
      end loop;
    end if;
    tsols := Truncate(sols,n);
    Clear(sols);
    if dim + add2dim > 0
     then sols := Add_Embedding(tsols,dim+add2dim); Clear(tsols);
     else sols := tsols;
    end if;
    Clear(pp);
    r := new Poly_Sys'(res);
  end Collapse_System;

  procedure Collapse_System
             ( p : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
               sols : in out QuadDobl_Complex_Solutions.Solution_List;
               dim,add2dim : in natural32;
               r : out QuadDobl_Complex_Poly_Systems.Link_to_Poly_Sys ) is

    use QuadDobl_Complex_Polynomials;
    use QuadDobl_Complex_Poly_Systems;
    use QuadDobl_Complex_Solutions;

    nb : constant natural32 := Symbol_Table.Number;
    n : constant integer32 := integer32(nb-dim)/2;
    s1,s2 : Array_of_Symbols(1..n);
    p1,p2,prm : Permutation(1..n);
    cbp : Permutation(1..integer32(nb));
    pp,cp : Poly_Sys(p'range);
    res : Poly_Sys(1..n+integer32(add2dim+dim));
    tsols : Solution_List;
    ind : integer32 := res'first;
    cnt : natural32 := 0;

  begin
    Retrieve_Suffixed_Symbols(n,'1',s1,p1);
    Retrieve_Suffixed_Symbols(n,'2',s2,p2);
    prm := Matching_Permutation(s1,s2);
    cbp := Combine_Permutations(n,integer32(dim),p1,p2,prm);
    pp := p*cbp;
    Permute(cbp,sols);
    cp := Collapse(pp,n);
    for i in s1'range loop
      s1(i) := Remove_Suffix(s1(i));
    end loop;
    Assign_Symbol_Table(s1);
    Add_Embed_Symbols(dim+add2dim);
    if dim + add2dim > 0 then
      for i in cp'range loop
        Add_Embedding(cp(i),natural32(n),dim+add2dim,cnt);
        if cp(i) /= Null_Poly then
          res(ind) := cp(i);
          ind := ind + 1;
        end if;
      end loop;
    else
      for i in res'range loop
        res(i) := cp(i);
      end loop;
    end if;
    tsols := Truncate(sols,n);
    Clear(sols);
    if dim + add2dim > 0
     then sols := Add_Embedding(tsols,dim+add2dim); Clear(tsols);
     else sols := tsols;
    end if;
    Clear(pp);
    r := new Poly_Sys'(res);
  end Collapse_System;

  procedure Standard_Collapse_Diagonal_System is

    use Standard_Complex_Poly_Systems;
    use Standard_Complex_Solutions;

    file : file_type;
    lp : Link_to_Poly_Sys;
    sols : Solution_List;
    dim,addtodim : natural32 := 0;

  begin
    new_line;
    put("Reading the diagonal system...");
    Standard_Read_Embedding(lp,sols,dim);
    new_line;
    put_line("Reading the name of the file to save the collapsed system.");
    Read_Name_and_Create_File(file);
    new_line;
    put("The dimension is "); put(dim,1); new_line;
    put("Give a natural number to add to the dimension : ");
    Read_Natural(addtodim);
   -- Write_Symbol_Table;
    declare
      cp : Link_to_Poly_Sys;
    begin
      Collapse_System(lp.all,sols,dim,addtodim,cp);
      put_line(file,cp.all);
      new_line(file);
      put_line(file,"THE SOLUTIONS :");
      put(file,Length_Of(sols),natural32(cp'last),sols);
    end;
  end Standard_Collapse_Diagonal_System;

  procedure DoblDobl_Collapse_Diagonal_System is

    use DoblDobl_Complex_Poly_Systems;
    use DoblDobl_Complex_Solutions;

    file : file_type;
    lp : Link_to_Poly_Sys;
    sols : Solution_List;
    dim,addtodim : natural32 := 0;

  begin
    new_line;
    put("Reading the diagonal system...");
    DoblDobl_Read_Embedding(lp,sols,dim);
    new_line;
    put_line("Reading the name of the file to save the collapsed system.");
    Read_Name_and_Create_File(file);
    new_line;
    put("The dimension is "); put(dim,1); new_line;
    put("Give a natural number to add to the dimension : ");
    Read_Natural(addtodim);
   -- Write_Symbol_Table;
    declare
      cp : Link_to_Poly_Sys;
    begin
      Collapse_System(lp.all,sols,dim,addtodim,cp);
      put_line(file,cp.all);
      new_line(file);
      put_line(file,"THE SOLUTIONS :");
      put(file,Length_Of(sols),natural32(cp'last),sols);
    end;
  end DoblDobl_Collapse_Diagonal_System;

  procedure QuadDobl_Collapse_Diagonal_System is

    use QuadDobl_Complex_Poly_Systems;
    use QuadDobl_Complex_Solutions;

    file : file_type;
    lp : Link_to_Poly_Sys;
    sols : Solution_List;
    dim,addtodim : natural32 := 0;

  begin
    new_line;
    put("Reading the diagonal system...");
    QuadDobl_Read_Embedding(lp,sols,dim);
    new_line;
    put_line("Reading the name of the file to save the collapsed system.");
    Read_Name_and_Create_File(file);
    new_line;
    put("The dimension is "); put(dim,1); new_line;
    put("Give a natural number to add to the dimension : ");
    Read_Natural(addtodim);
   -- Write_Symbol_Table;
    declare
      cp : Link_to_Poly_Sys;
    begin
      Collapse_System(lp.all,sols,dim,addtodim,cp);
      put_line(file,cp.all);
      new_line(file);
      put_line(file,"THE SOLUTIONS :");
      put(file,Length_Of(sols),natural32(cp'last),sols);
    end;
  end QuadDobl_Collapse_Diagonal_System;

  procedure Collapse_Diagonal_System is

    p : constant character := Prompt_for_Precision;

  begin
    case p is
      when '0' => Standard_Collapse_Diagonal_System; 
      when '1' => DoblDobl_Collapse_Diagonal_System; 
      when '2' => QuadDobl_Collapse_Diagonal_System; 
      when others => null;
    end case;
  end Collapse_Diagonal_System;

end Extrinsic_Diagonal_Solvers;
