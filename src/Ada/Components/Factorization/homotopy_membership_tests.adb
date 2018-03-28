with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;
with Standard_Complex_Numbers;           use Standard_Complex_Numbers;
with Double_Double_Numbers;              use Double_Double_Numbers;
with Double_Double_Numbers_io;           use Double_Double_Numbers_io;
with Quad_Double_Numbers;                use Quad_Double_Numbers;
with Quad_Double_Numbers_io;             use Quad_Double_Numbers_io;
with DoblDobl_Complex_Numbers;           use DoblDobl_Complex_Numbers;
with QuadDobl_Complex_Numbers;           use QuadDobl_Complex_Numbers;
with Standard_Complex_Norms_Equals;      use Standard_Complex_Norms_Equals;
with DoblDobl_Complex_Vector_Norms;      use DoblDobl_Complex_Vector_Norms;
with QuadDobl_Complex_Vector_Norms;      use QuadDobl_Complex_Vector_Norms;
with Standard_Complex_Poly_SysFun;
with Standard_Complex_Laur_SysFun;
with DoblDobl_Complex_Poly_SysFun;
with DoblDobl_Complex_Laur_SysFun;
with QuadDobl_Complex_Poly_SysFun;
with QuadDobl_Complex_Laur_SysFun;
with Standard_Complex_Solutions_io;      use Standard_Complex_Solutions_io;
with DoblDobl_Complex_Solutions_io;      use DoblDobl_Complex_Solutions_io;
with QuadDobl_Complex_Solutions_io;      use QuadDobl_Complex_Solutions_io;
with Witness_Sets;                       use Witness_Sets;
with Sampling_Machine;
with Sampling_Laurent_Machine;
with DoblDobl_Sampling_Machine;
with DoblDobl_Sampling_Laurent_Machine;
with QuadDobl_Sampling_Machine;
with QuadDobl_Sampling_Laurent_Machine;
with Homotopy_Membership_Target;         use Homotopy_Membership_Target;

package body Homotopy_Membership_Tests is

-- AUXILIARY OPERATIONS :

  function Embed_Solution
             ( s : Standard_Complex_Solutions.Solution;
               n : integer32 )
             return Standard_Complex_Solutions.Solution is

  -- DESCRIPTION :
  --   Adds as many slack variables to the solution s so that its
  --   number of variables equals the ambient dimension n.

  begin
    if s.n >= n
     then return s;
     else return Add_Embedding(s,natural32(n-s.n));
    end if;
  end Embed_Solution;

  function Embed_Solution
             ( s : DoblDobl_Complex_Solutions.Solution;
               n : integer32 )
             return DoblDobl_Complex_Solutions.Solution is

  -- DESCRIPTION :
  --   Adds as many slack variables to the solution s so that its
  --   number of variables equals the ambient dimension n.

  begin
    if s.n >= n
     then return s;
     else return Add_Embedding(s,natural32(n-s.n));
    end if;
  end Embed_Solution;

  function Embed_Solution
             ( s : QuadDobl_Complex_Solutions.Solution;
               n : integer32 )
             return QuadDobl_Complex_Solutions.Solution is

  -- DESCRIPTION :
  --   Adds as many slack variables to the solution s so that its
  --   number of variables equals the ambient dimension n.

  begin
    if s.n >= n
     then return s;
     else return Add_Embedding(s,natural32(n-s.n));
    end if;
  end Embed_Solution;

  function Is_Member
              ( verbose : boolean; idx : natural32;
                diff,homtol : double_float ) return boolean is

  -- DESCRIPTION:
  --   Returns true if the difference diff is within the tolerance,
  --   returns false otherwise.

  -- ON INPUT :
  --   verbose  if true, then results are written to standard output;
  --   idx      index of the test solution;
  --   diff     difference between test and new generic point;
  --   homtol   tolerance for equality.

    found : boolean := false;

  begin
    if diff <= homtol then
      found := true;
      if verbose then
        put("    match with generic point "); put(idx,1);
        put(", as difference is");
        put(diff,3); put(" <="); put(homtol,3); new_line;
      end if;
    else
      if verbose then
        put("    difference with generic point "); put(idx,1);
        put(" is "); put(diff,3); put(" >"); put(homtol,3); new_line;
      end if;
    end if;
    return found;
  end Is_Member;

  function Is_Member
              ( verbose : boolean; idx : natural32;
                diff : double_double; homtol : double_float )
              return boolean is

  -- DESCRIPTION:
  --   Returns true if the difference diff is within the tolerance,
  --   returns false otherwise.

  -- ON INPUT :
  --   verbose  if true, then results are written to standard output;
  --   idx      index of the test solution;
  --   diff     difference between test and new generic point;
  --   homtol   tolerance for equality.

    found : boolean := false;

  begin
    if diff <= homtol then
      found := true;
      if verbose then
        put("    match with generic point "); put(idx,1);
        put(", as difference is");
        put(diff,3); put(" <="); put(homtol,3); new_line;
      end if;
    else
      if verbose then
        put("    difference with generic point "); put(idx,1);
        put(" is "); put(diff,3); put(" >"); put(homtol,3); new_line;
      end if;
    end if;
    return found;
  end Is_Member;

  function Is_Member
              ( verbose : boolean; idx : natural32;
                diff : quad_double; homtol : double_float )
              return boolean is

  -- DESCRIPTION:
  --   Returns true if the difference diff is within the tolerance,
  --   returns false otherwise.

  -- ON INPUT :
  --   verbose  if true, then results are written to standard output;
  --   idx      index of the test solution;
  --   diff     difference between test and new generic point;
  --   homtol   tolerance for equality.

    found : boolean := false;

  begin
    if diff <= homtol then
      found := true;
      if verbose then
        put("    match with generic point "); put(idx,1);
        put(", as difference is");
        put(diff,3); put(" <="); put(homtol,3); new_line;
      end if;
    else
      if verbose then
        put("    difference with generic point "); put(idx,1);
        put(" is "); put(diff,3); put(" >"); put(homtol,3); new_line;
      end if;
    end if;
    return found;
  end Is_Member;

  function Is_Member
              ( file : file_type; idx : natural32;
                diff,homtol : double_float ) return boolean is

  -- DESCRIPTION:
  --   Returns true if the difference diff is within the tolerance,
  --   returns false otherwise.

  -- ON INPUT :
  --   verbose  if true, then results are written to standard output;
  --   idx      index of the test solution;
  --   diff     difference between test and new generic point;
  --   homtol   tolerance for equality.

    found : boolean := false;

  begin
    if diff <= homtol then
      put(file,"    match with generic point "); put(file,idx,1);
      put(file,", as difference is");
      put(file,diff,3); put(file," <="); put(file,homtol,3);
      new_line(file);
      found := true;
    else
      put(file,"    difference with generic point "); put(file,idx,1);
      put(file," is ");
      put(file,diff,3); put(file," >"); put(file,homtol,3);
      new_line(file);
    end if;
    return found;
  end Is_Member;

  function Is_Member
              ( file : file_type; idx : natural32;
                diff : double_double; homtol : double_float )
              return boolean is

  -- DESCRIPTION:
  --   Returns true if the difference diff is within the tolerance,
  --   returns false otherwise.

  -- ON INPUT :
  --   verbose  if true, then results are written to standard output;
  --   idx      index of the test solution;
  --   diff     difference between test and new generic point;
  --   homtol   tolerance for equality.

    found : boolean := false;

  begin
    if diff <= homtol then
      put(file,"    match with generic point "); put(file,idx,1);
      put(file,", as difference is");
      put(file,diff,3); put(file," <="); put(file,homtol,3);
      new_line(file);
      found := true;
    else
      put(file,"    difference with generic point "); put(file,idx,1);
      put(file," is ");
      put(file,diff,3); put(file," >"); put(file,homtol,3);
      new_line(file);
    end if;
    return found;
  end Is_Member;

  function Is_Member
              ( file : file_type; idx : natural32;
                diff : quad_double; homtol : double_float )
              return boolean is

  -- DESCRIPTION:
  --   Returns true if the difference diff is within the tolerance,
  --   returns false otherwise.

  -- ON INPUT :
  --   verbose  if true, then results are written to standard output;
  --   idx      index of the test solution;
  --   diff     difference between test and new generic point;
  --   homtol   tolerance for equality.

    found : boolean := false;

  begin
    if diff <= homtol then
      put(file,"    match with generic point "); put(file,idx,1);
      put(file,", as difference is");
      put(file,diff,3); put(file," <="); put(file,homtol,3);
      new_line(file);
      found := true;
    else
      put(file,"    difference with generic point "); put(file,idx,1);
      put(file," is ");
      put(file,diff,3); put(file," >"); put(file,homtol,3);
      new_line(file);
    end if;
    return found;
  end Is_Member;

  procedure Write_Outcome ( verbose,found : in boolean ) is

  -- DESCRIPTION :
  --   Writes the outcome of the homotopy membership test in found,
  --   if verbose.

  begin
    if verbose then
      if found
       then put_line("  Point lies on a solution component.");
       else put_line("  Point does not lie on a solution component.");
      end if;
    end if;
  end Write_Outcome;

  procedure Write_Outcome ( file : in file_type; found : in boolean ) is

  -- DESCRIPTION :
  --   Writes the outcome of the homotopy membership test in found,
  --   to file.

  begin
    if found
     then put_line(file,"  Point lies on a solution component.");
     else put_line(file,"  Point does not lie on a solution component.");
    end if;
  end Write_Outcome;

-- TARGET ROUTINES :

  procedure Homotopy_Membership_Test
              ( verbose : in boolean;
                genpts : in Standard_Complex_Solutions.Solution_List;
                x : in Standard_Complex_Vectors.Vector;
                adjsli : in Standard_Complex_VecVecs.VecVec;
                homtol : in double_float; found : out boolean ) is

    newsols,tmp : Standard_Complex_Solutions.Solution_List;
    test_sol : Standard_Complex_Solutions.Link_to_Solution;
    difference : Standard_Complex_Vectors.Vector(x'range);
    diff : double_float;

  begin
    if verbose then
      Sampling_Machine.Sample_with_Stop
        (standard_output,genpts,x,homtol,adjsli,newsols);
    else
      Sampling_Machine.Sample_with_Stop
        (genpts,x,homtol,adjsli,newsols);
    end if;
    tmp := newsols;
    for i in 1..Standard_Complex_Solutions.Length_Of(newsols) loop
      test_sol := Standard_Complex_Solutions.Head_Of(tmp);
      difference := test_sol.v(x'range) - x;
      diff := Max_Norm(difference);
      found := Is_Member(verbose,i,diff,homtol);
      exit when found;
      tmp := Standard_Complex_Solutions.Tail_Of(tmp);
    end loop;
    Write_Outcome(verbose,found);
    Standard_Complex_Solutions.Clear(newsols);
  end Homotopy_Membership_Test;

  procedure Laurent_Homotopy_Membership_Test
              ( verbose : in boolean;
                genpts : in Standard_Complex_Solutions.Solution_List;
                x : in Standard_Complex_Vectors.Vector;
                adjsli : in Standard_Complex_VecVecs.VecVec;
                homtol : in double_float; found : out boolean ) is

    newsols,tmp : Standard_Complex_Solutions.Solution_List;
    test_sol : Standard_Complex_Solutions.Link_to_Solution;
    difference : Standard_Complex_Vectors.Vector(x'range);
    diff : double_float;

  begin
    if verbose then
      Sampling_Laurent_Machine.Sample_with_Stop
        (standard_output,genpts,x,homtol,adjsli,newsols);
    else
      Sampling_Laurent_Machine.Sample_with_Stop
        (genpts,x,homtol,adjsli,newsols);
    end if;
    tmp := newsols;
    for i in 1..Standard_Complex_Solutions.Length_Of(newsols) loop
      test_sol := Standard_Complex_Solutions.Head_Of(tmp);
      difference := test_sol.v(x'range) - x;
      diff := Max_Norm(difference);
      found := Is_Member(verbose,i,diff,homtol);
      exit when found;
      tmp := Standard_Complex_Solutions.Tail_Of(tmp);
    end loop;
    Write_Outcome(verbose,found);
    Standard_Complex_Solutions.Clear(newsols);
  end Laurent_Homotopy_Membership_Test;

  procedure Homotopy_Membership_Test
              ( verbose : in boolean;
                genpts : in DoblDobl_Complex_Solutions.Solution_List;
                x : in DoblDobl_Complex_Vectors.Vector;
                adjsli : in DoblDobl_Complex_VecVecs.VecVec;
                homtol : in double_float; found : out boolean ) is

    newsols,tmp : DoblDobl_Complex_Solutions.Solution_List;
    test_sol : DoblDobl_Complex_Solutions.Link_to_Solution;
    difference : DoblDobl_Complex_Vectors.Vector(x'range);
    diff : double_double;

  begin
    if verbose then
      DoblDobl_Sampling_Machine.Sample_with_Stop
        (standard_output,genpts,x,homtol,adjsli,newsols);
    else
      DoblDobl_Sampling_Machine.Sample_with_Stop
        (genpts,x,homtol,adjsli,newsols);
    end if;
    tmp := newsols;
    for i in 1..DoblDobl_Complex_Solutions.Length_Of(newsols) loop
      test_sol := DoblDobl_Complex_Solutions.Head_Of(tmp);
      difference := test_sol.v(x'range) - x;
      diff := Max_Norm(difference);
      found := Is_Member(verbose,i,diff,homtol);
      exit when found;
      tmp := DoblDobl_Complex_Solutions.Tail_Of(tmp);
    end loop;
    Write_Outcome(verbose,found);
    DoblDobl_Complex_Solutions.Clear(newsols);
  end Homotopy_Membership_Test;

  procedure Laurent_Homotopy_Membership_Test
              ( verbose : in boolean;
                genpts : in DoblDobl_Complex_Solutions.Solution_List;
                x : in DoblDobl_Complex_Vectors.Vector;
                adjsli : in DoblDobl_Complex_VecVecs.VecVec;
                homtol : in double_float; found : out boolean ) is

    newsols,tmp : DoblDobl_Complex_Solutions.Solution_List;
    test_sol : DoblDobl_Complex_Solutions.Link_to_Solution;
    difference : DoblDobl_Complex_Vectors.Vector(x'range);
    diff : double_double;

  begin
    if verbose then
      DoblDobl_Sampling_Laurent_Machine.Sample_with_Stop
        (standard_output,genpts,x,homtol,adjsli,newsols);
    else
      DoblDobl_Sampling_Laurent_Machine.Sample_with_Stop
        (genpts,x,homtol,adjsli,newsols);
    end if;
    tmp := newsols;
    for i in 1..DoblDobl_Complex_Solutions.Length_Of(newsols) loop
      test_sol := DoblDobl_Complex_Solutions.Head_Of(tmp);
      difference := test_sol.v(x'range) - x;
      diff := Max_Norm(difference);
      found := Is_Member(verbose,i,diff,homtol);
      exit when found;
      tmp := DoblDobl_Complex_Solutions.Tail_Of(tmp);
    end loop;
    Write_Outcome(verbose,found);
    DoblDobl_Complex_Solutions.Clear(newsols);
  end Laurent_Homotopy_Membership_Test;

  procedure Homotopy_Membership_Test
              ( verbose : in boolean;
                genpts : in QuadDobl_Complex_Solutions.Solution_List;
                x : in QuadDobl_Complex_Vectors.Vector;
                adjsli : in QuadDobl_Complex_VecVecs.VecVec;
                homtol : in double_float; found : out boolean ) is

    newsols,tmp : QuadDobl_Complex_Solutions.Solution_List;
    test_sol : QuadDobl_Complex_Solutions.Link_to_Solution;
    difference : QuadDobl_Complex_Vectors.Vector(x'range);
    diff : quad_double;

  begin
    if verbose then
      QuadDobl_Sampling_Machine.Sample_with_Stop
        (standard_output,genpts,x,homtol,adjsli,newsols);
    else
      QuadDobl_Sampling_Machine.Sample_with_Stop
        (genpts,x,homtol,adjsli,newsols);
    end if;
    tmp := newsols;
    for i in 1..QuadDobl_Complex_Solutions.Length_Of(newsols) loop
      test_sol := QuadDobl_Complex_Solutions.Head_Of(tmp);
      difference := test_sol.v(x'range) - x;
      diff := Max_Norm(difference);
      found := Is_Member(verbose,i,diff,homtol);
      exit when found;
      tmp := QuadDobl_Complex_Solutions.Tail_Of(tmp);
    end loop;
    Write_Outcome(verbose,found);
    QuadDobl_Complex_Solutions.Clear(newsols);
  end Homotopy_Membership_Test;

  procedure Laurent_Homotopy_Membership_Test
              ( verbose : in boolean;
                genpts : in QuadDobl_Complex_Solutions.Solution_List;
                x : in QuadDobl_Complex_Vectors.Vector;
                adjsli : in QuadDobl_Complex_VecVecs.VecVec;
                homtol : in double_float; found : out boolean ) is

    newsols,tmp : QuadDobl_Complex_Solutions.Solution_List;
    test_sol : QuadDobl_Complex_Solutions.Link_to_Solution;
    difference : QuadDobl_Complex_Vectors.Vector(x'range);
    diff : quad_double;

  begin
    if verbose then
      QuadDobl_Sampling_Laurent_Machine.Sample_with_Stop
        (standard_output,genpts,x,homtol,adjsli,newsols);
    else
      QuadDobl_Sampling_Laurent_Machine.Sample_with_Stop
        (genpts,x,homtol,adjsli,newsols);
    end if;
    tmp := newsols;
    for i in 1..QuadDobl_Complex_Solutions.Length_Of(newsols) loop
      test_sol := QuadDobl_Complex_Solutions.Head_Of(tmp);
      difference := test_sol.v(x'range) - x;
      diff := Max_Norm(difference);
      found := Is_Member(verbose,i,diff,homtol);
      exit when found;
      tmp := QuadDobl_Complex_Solutions.Tail_Of(tmp);
    end loop;
    Write_Outcome(verbose,found);
    QuadDobl_Complex_Solutions.Clear(newsols);
  end Laurent_Homotopy_Membership_Test;

  procedure Homotopy_Membership_Test
              ( file : in file_type;
                genpts : in Standard_Complex_Solutions.Solution_List;
                x : in Standard_Complex_Vectors.Vector;
                adjsli : in Standard_Complex_VecVecs.VecVec;
                homtol : in double_float; found : out boolean ) is

    newsols,tmp : Standard_Complex_Solutions.Solution_List;
    test_sol : Standard_Complex_Solutions.Link_to_Solution;
    difference : Standard_Complex_Vectors.Vector(x'range);
    diff : double_float;

  begin
    Sampling_Machine.Sample_with_Stop
      (file,genpts,x,homtol,adjsli,newsols);
    put_line(file,"  Scanning the new list of generic points :");
    tmp := newsols;
    for i in 1..Standard_Complex_Solutions.Length_Of(newsols) loop
      test_sol := Standard_Complex_Solutions.Head_Of(tmp);
      difference := test_sol.v(x'range) - x;
      diff := Max_Norm(difference);
      found := Is_Member(file,i,diff,homtol);
      exit when found;
      tmp := Standard_Complex_Solutions.Tail_Of(tmp);
    end loop;
    Write_Outcome(file,found);
    Standard_Complex_Solutions.Clear(newsols);
  end Homotopy_Membership_Test;

  procedure Laurent_Homotopy_Membership_Test
              ( file : in file_type;
                genpts : in Standard_Complex_Solutions.Solution_List;
                x : in Standard_Complex_Vectors.Vector;
                adjsli : in Standard_Complex_VecVecs.VecVec;
                homtol : in double_float; found : out boolean ) is

    newsols,tmp : Standard_Complex_Solutions.Solution_List;
    test_sol : Standard_Complex_Solutions.Link_to_Solution;
    difference : Standard_Complex_Vectors.Vector(x'range);
    diff : double_float;

  begin
    Sampling_Laurent_Machine.Sample_with_Stop
      (file,genpts,x,homtol,adjsli,newsols);
    put_line(file,"  Scanning the new list of generic points :");
    tmp := newsols;
    for i in 1..Standard_Complex_Solutions.Length_Of(newsols) loop
      test_sol := Standard_Complex_Solutions.Head_Of(tmp);
      difference := test_sol.v(x'range) - x;
      diff := Max_Norm(difference);
      found := Is_Member(file,i,diff,homtol);
      exit when found;
      tmp := Standard_Complex_Solutions.Tail_Of(tmp);
    end loop;
    Write_Outcome(file,found);
    Standard_Complex_Solutions.Clear(newsols);
  end Laurent_Homotopy_Membership_Test;

  procedure Homotopy_Membership_Test
              ( file : in file_type;
                genpts : in DoblDobl_Complex_Solutions.Solution_List;
                x : in DoblDobl_Complex_Vectors.Vector;
                adjsli : in DoblDobl_Complex_VecVecs.VecVec;
                homtol : in double_float; found : out boolean ) is

    newsols,tmp : DoblDobl_Complex_Solutions.Solution_List;
    test_sol : DoblDobl_Complex_Solutions.Link_to_Solution;
    difference : DoblDobl_Complex_Vectors.Vector(x'range);
    diff : double_double;

  begin
    DoblDobl_Sampling_Machine.Sample_with_Stop
      (file,genpts,x,homtol,adjsli,newsols);
    put_line(file,"  Scanning the new list of generic points :");
    tmp := newsols;
    for i in 1..DoblDobl_Complex_Solutions.Length_Of(newsols) loop
      test_sol := DoblDobl_Complex_Solutions.Head_Of(tmp);
      difference := test_sol.v(x'range) - x;
      diff := Max_Norm(difference);
      found := Is_Member(file,i,diff,homtol);
      exit when found;
      tmp := DoblDobl_Complex_Solutions.Tail_Of(tmp);
    end loop;
    Write_Outcome(file,found);
    DoblDobl_Complex_Solutions.Clear(newsols);
  end Homotopy_Membership_Test;

  procedure Laurent_Homotopy_Membership_Test
              ( file : in file_type;
                genpts : in DoblDobl_Complex_Solutions.Solution_List;
                x : in DoblDobl_Complex_Vectors.Vector;
                adjsli : in DoblDobl_Complex_VecVecs.VecVec;
                homtol : in double_float; found : out boolean ) is

    newsols,tmp : DoblDobl_Complex_Solutions.Solution_List;
    test_sol : DoblDobl_Complex_Solutions.Link_to_Solution;
    difference : DoblDobl_Complex_Vectors.Vector(x'range);
    diff : double_double;

  begin
    DoblDobl_Sampling_Laurent_Machine.Sample_with_Stop
      (file,genpts,x,homtol,adjsli,newsols);
    put_line(file,"  Scanning the new list of generic points :");
    tmp := newsols;
    for i in 1..DoblDobl_Complex_Solutions.Length_Of(newsols) loop
      test_sol := DoblDobl_Complex_Solutions.Head_Of(tmp);
      difference := test_sol.v(x'range) - x;
      diff := Max_Norm(difference);
      found := Is_Member(file,i,diff,homtol);
      exit when found;
      tmp := DoblDobl_Complex_Solutions.Tail_Of(tmp);
    end loop;
    Write_Outcome(file,found);
    DoblDobl_Complex_Solutions.Clear(newsols);
  end Laurent_Homotopy_Membership_Test;

  procedure Homotopy_Membership_Test
              ( file : in file_type;
                genpts : in QuadDobl_Complex_Solutions.Solution_List;
                x : in QuadDobl_Complex_Vectors.Vector;
                adjsli : in QuadDobl_Complex_VecVecs.VecVec;
                homtol : in double_float; found : out boolean ) is

    newsols,tmp : QuadDobl_Complex_Solutions.Solution_List;
    test_sol : QuadDobl_Complex_Solutions.Link_to_Solution;
    difference : QuadDobl_Complex_Vectors.Vector(x'range);
    diff : quad_double;

  begin
    QuadDobl_Sampling_Machine.Sample_with_Stop
      (file,genpts,x,homtol,adjsli,newsols);
    put_line(file,"  Scanning the new list of generic points :");
    tmp := newsols;
    for i in 1..QuadDobl_Complex_Solutions.Length_Of(newsols) loop
      test_sol := QuadDobl_Complex_Solutions.Head_Of(tmp);
      difference := test_sol.v(x'range) - x;
      diff := Max_Norm(difference);
      found := Is_Member(file,i,diff,homtol);
      exit when found;
      tmp := QuadDobl_Complex_Solutions.Tail_Of(tmp);
    end loop;
    Write_Outcome(file,found);
    QuadDobl_Complex_Solutions.Clear(newsols);
  end Homotopy_Membership_Test;

  procedure Laurent_Homotopy_Membership_Test
              ( file : in file_type;
                genpts : in QuadDobl_Complex_Solutions.Solution_List;
                x : in QuadDobl_Complex_Vectors.Vector;
                adjsli : in QuadDobl_Complex_VecVecs.VecVec;
                homtol : in double_float; found : out boolean ) is

    newsols,tmp : QuadDobl_Complex_Solutions.Solution_List;
    test_sol : QuadDobl_Complex_Solutions.Link_to_Solution;
    difference : QuadDobl_Complex_Vectors.Vector(x'range);
    diff : quad_double;

  begin
    QuadDobl_Sampling_Laurent_Machine.Sample_with_Stop
      (file,genpts,x,homtol,adjsli,newsols);
    put_line(file,"  Scanning the new list of generic points :");
    tmp := newsols;
    for i in 1..QuadDobl_Complex_Solutions.Length_Of(newsols) loop
      test_sol := QuadDobl_Complex_Solutions.Head_Of(tmp);
      difference := test_sol.v(x'range) - x;
      diff := Max_Norm(difference);
      found := Is_Member(file,i,diff,homtol);
      exit when found;
      tmp := QuadDobl_Complex_Solutions.Tail_Of(tmp);
    end loop;
    Write_Outcome(file,found);
    QuadDobl_Complex_Solutions.Clear(newsols);
  end Laurent_Homotopy_Membership_Test;

  procedure Homotopy_Membership_Test 
              ( verbose : in boolean;
                ep : in Standard_Complex_Poly_Systems.Poly_Sys;
                dim : in natural32;
                sli : in Standard_Complex_VecVecs.VecVec;
                genpts : in Standard_Complex_Solutions.Solution_List;
                x : in Standard_Complex_Vectors.Vector;
                restol,homtol : in double_float;
                success,found : out boolean ) is

    ex : Standard_Complex_Vectors.Vector(ep'range);
    y : Standard_Complex_Vectors.Vector(ep'range);
    res : double_float;
    newsli : Standard_Complex_VecVecs.VecVec(sli'range);

  begin
    ex(x'range) := x;
    for k in x'last+1..ex'last loop
      ex(k) := Standard_Complex_Numbers.Create(0.0);
    end loop;
    y := Standard_Complex_Poly_SysFun.Eval(ep,ex);
    res := Max_Norm(y(y'first..y'last-integer32(dim)));
    if verbose
     then put("residual is "); put(res); new_line;
    end if;
    if res <= restol then
      success := true;
      if verbose then
        put("  point satisfies the equations, as residual <=");
        put(restol,3); new_line;
      end if;
      newsli := Adjusted_Slices(sli,ex);
      Homotopy_Membership_Test(verbose,genpts,ex,newsli,homtol,found);
    else
      success := false;
      if verbose then
        put("  point does not lie on the component, as residual >");
        put(restol,3); new_line;
      end if;
      found := false;
    end if;
  end Homotopy_Membership_Test;

  procedure Homotopy_Membership_Test 
              ( verbose : in boolean;
                ep : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                dim : in natural32;
                sli : in DoblDobl_Complex_VecVecs.VecVec;
                genpts : in DoblDobl_Complex_Solutions.Solution_List;
                x : in DoblDobl_Complex_Vectors.Vector;
                restol,homtol : in double_float;
                success,found : out boolean ) is

    ex : DoblDobl_Complex_Vectors.Vector(ep'range);
    y : DoblDobl_Complex_Vectors.Vector(ep'range);
    res : double_double;
    newsli : DoblDobl_Complex_VecVecs.VecVec(sli'range);

  begin
    ex(x'range) := x;
    for k in x'last+1..ex'last loop
      ex(k) := DoblDobl_Complex_Numbers.Create(integer(0));
    end loop;
    y := DoblDobl_Complex_Poly_SysFun.Eval(ep,ex);
    res := Max_Norm(y(y'first..y'last-integer32(dim)));
    if verbose
     then put("residual is "); put(res); new_line;
    end if;
    if res <= restol then
      success := true;
      if verbose then
        put("  point satisfies the equations, as residual <=");
        put(restol,3); new_line;
      end if;
      newsli := Adjusted_Slices(sli,ex);
      Homotopy_Membership_Test(verbose,genpts,ex,newsli,homtol,found);
    else
      success := false;
      if verbose then
        put("  point does not lie on the component, as residual >");
        put(restol,3); new_line;
      end if;
      found := false;
    end if;
  end Homotopy_Membership_Test;

  procedure Homotopy_Membership_Test 
              ( verbose : in boolean;
                ep : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                dim : in natural32;
                sli : in QuadDobl_Complex_VecVecs.VecVec;
                genpts : in QuadDobl_Complex_Solutions.Solution_List;
                x : in QuadDobl_Complex_Vectors.Vector;
                restol,homtol : in double_float;
                success,found : out boolean ) is

    ex : QuadDobl_Complex_Vectors.Vector(ep'range);
    y : QuadDobl_Complex_Vectors.Vector(ep'range);
    res : quad_double;
    newsli : QuadDobl_Complex_VecVecs.VecVec(sli'range);

  begin
    ex(x'range) := x;
    for k in x'last+1..ex'last loop
      ex(k) := QuadDobl_Complex_Numbers.Create(integer(0));
    end loop;
    y := QuadDobl_Complex_Poly_SysFun.Eval(ep,ex);
    res := Max_Norm(y(y'first..y'last-integer32(dim)));
    if verbose
     then put("residual is "); put(res); new_line;
    end if;
    if res <= restol then
      success := true;
      if verbose then
        put("  point satisfies the equations, as residual <=");
        put(restol,3); new_line;
      end if;
      newsli := Adjusted_Slices(sli,ex);
      Homotopy_Membership_Test(verbose,genpts,ex,newsli,homtol,found);
    else
      success := false;
      if verbose then
        put("  point does not lie on the component, as residual >");
        put(restol,3); new_line;
      end if;
      found := false;
    end if;
  end Homotopy_Membership_Test;

  procedure Homotopy_Membership_Test 
              ( verbose : in boolean;
                ep : in Standard_Complex_Laur_Systems.Laur_Sys;
                dim : in natural32;
                sli : in Standard_Complex_VecVecs.VecVec;
                genpts : in Standard_Complex_Solutions.Solution_List;
                x : in Standard_Complex_Vectors.Vector;
                restol,homtol : in double_float;
                success,found : out boolean ) is

    ex : Standard_Complex_Vectors.Vector(ep'range);
    y : Standard_Complex_Vectors.Vector(ep'range);
    res : double_float;
    newsli : Standard_Complex_VecVecs.VecVec(sli'range);

  begin
    ex(x'range) := x;
    for k in x'last+1..ex'last loop
      ex(k) := Standard_Complex_Numbers.Create(0.0);
    end loop;
    y := Standard_Complex_Laur_SysFun.Eval(ep,ex);
    res := Max_Norm(y(y'first..y'last-integer32(dim)));
    if verbose
     then put("residual is "); put(res); new_line;
    end if;
    if res <= restol then
      success := true;
      if verbose then
        put("  point satisfies the equations, as residual <=");
        put(restol,3); new_line;
      end if;
      newsli := Adjusted_Slices(sli,ex);
      Laurent_Homotopy_Membership_Test(verbose,genpts,ex,newsli,homtol,found);
    else
      success := false;
      if verbose then
        put("  point does not lie on the component, as residual >");
        put(restol,3); new_line;
      end if;
      found := false;
    end if;
  end Homotopy_Membership_Test;

  procedure Homotopy_Membership_Test 
              ( verbose : in boolean;
                ep : in DoblDobl_Complex_Laur_Systems.Laur_Sys;
                dim : in natural32;
                sli : in DoblDobl_Complex_VecVecs.VecVec;
                genpts : in DoblDobl_Complex_Solutions.Solution_List;
                x : in DoblDobl_Complex_Vectors.Vector;
                restol,homtol : in double_float;
                success,found : out boolean ) is

    ex : DoblDobl_Complex_Vectors.Vector(ep'range);
    y : DoblDobl_Complex_Vectors.Vector(ep'range);
    res : double_double;
    newsli : DoblDobl_Complex_VecVecs.VecVec(sli'range);

  begin
    ex(x'range) := x;
    for k in x'last+1..ex'last loop
      ex(k) := DoblDobl_Complex_Numbers.Create(integer(0));
    end loop;
    y := DoblDobl_Complex_Laur_SysFun.Eval(ep,ex);
    res := Max_Norm(y(y'first..y'last-integer32(dim)));
    if verbose
     then put("residual is "); put(res); new_line;
    end if;
    if res <= restol then
      success := true;
      if verbose then
        put("  point satisfies the equations, as residual <=");
        put(restol,3); new_line;
      end if;
      newsli := Adjusted_Slices(sli,ex);
      Laurent_Homotopy_Membership_Test(verbose,genpts,ex,newsli,homtol,found);
    else
      success := false;
      if verbose then
        put("  point does not lie on the component, as residual >");
        put(restol,3); new_line;
      end if;
      found := false;
    end if;
  end Homotopy_Membership_Test;

  procedure Homotopy_Membership_Test 
              ( verbose : in boolean;
                ep : in QuadDobl_Complex_Laur_Systems.Laur_Sys;
                dim : in natural32;
                sli : in QuadDobl_Complex_VecVecs.VecVec;
                genpts : in QuadDobl_Complex_Solutions.Solution_List;
                x : in QuadDobl_Complex_Vectors.Vector;
                restol,homtol : in double_float;
                success,found : out boolean ) is

    ex : QuadDobl_Complex_Vectors.Vector(ep'range);
    y : QuadDobl_Complex_Vectors.Vector(ep'range);
    res : quad_double;
    newsli : QuadDobl_Complex_VecVecs.VecVec(sli'range);

  begin
    ex(x'range) := x;
    for k in x'last+1..ex'last loop
      ex(k) := QuadDobl_Complex_Numbers.Create(integer(0));
    end loop;
    y := QuadDobl_Complex_Laur_SysFun.Eval(ep,ex);
    res := Max_Norm(y(y'first..y'last-integer32(dim)));
    if verbose
     then put("residual is "); put(res); new_line;
    end if;
    if res <= restol then
      success := true;
      if verbose then
        put("  point satisfies the equations, as residual <=");
        put(restol,3); new_line;
      end if;
      newsli := Adjusted_Slices(sli,ex);
      Laurent_Homotopy_Membership_Test(verbose,genpts,ex,newsli,homtol,found);
    else
      success := false;
      if verbose then
        put("  point does not lie on the component, as residual >");
        put(restol,3); new_line;
      end if;
      found := false;
    end if;
  end Homotopy_Membership_Test;

  procedure Homotopy_Membership_Test 
              ( file : in file_type;
                ep : in Standard_Complex_Poly_Systems.Poly_Sys;
                dim : in natural32;
                sli : in Standard_Complex_VecVecs.VecVec;
                genpts : in Standard_Complex_Solutions.Solution_List;
                s : in Standard_Complex_Solutions.Solution;
                restol,homtol : in double_float;
                success,found : out boolean ) is

   -- es : Solution(s.n + dim) := Embed(s,dim);
    es : constant Standard_Complex_Solutions.Solution(ep'last)
       := Embed_Solution(s,ep'last);
    y : constant Standard_Complex_Vectors.Vector(ep'range)
      := Standard_Complex_Poly_SysFun.Eval(ep,es.v);
    res : constant double_float := Max_Norm(y(y'first..y'last-integer32(dim)));
    newsli : Standard_Complex_VecVecs.VecVec(sli'range);
   -- eva : Complex_Number;

  begin
   -- put_line(file,"The solution vector : "); put_line(file,es.v);
   -- put_line(file,"The evaluation : "); put_line(file,y);
    put(file,"residual is "); put(file,res); new_line(file);
    if res <= restol then
      put(file,"  point satisfies the equations, as residual <=");
      put(file,restol,3); new_line(file);
      success := true;
      newsli := Adjusted_Slices(sli,es.v);
     -- for i in newsli'range loop
     --   put(file,"eval at adjusted slice : ");
     --   eva := newsli(i)(0) + newsli(i)(es.v'range)*es.v;
     --   put(file,eva); new_line(file);
     -- end loop;
      Homotopy_Membership_Test(file,genpts,es.v,newsli,homtol,found);
    else
      put(file,"  point does not lie on the component, as residual >");
      put(file,restol,3); new_line(file);
      success := false;
      found := false;
    end if;
  end Homotopy_Membership_Test;

  procedure Homotopy_Membership_Test 
              ( file : in file_type;
                ep : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                dim : in natural32;
                sli : in DoblDobl_Complex_VecVecs.VecVec;
                genpts : in DoblDobl_Complex_Solutions.Solution_List;
                s : in DoblDobl_Complex_Solutions.Solution;
                restol,homtol : in double_float;
                success,found : out boolean ) is

   -- es : Solution(s.n + dim) := Embed(s,dim);
    es : constant DoblDobl_Complex_Solutions.Solution(ep'last)
       := Embed_Solution(s,ep'last);
    y : constant DoblDobl_Complex_Vectors.Vector(ep'range)
      := DoblDobl_Complex_Poly_SysFun.Eval(ep,es.v);
    res : constant double_double
        := Max_Norm(y(y'first..y'last-integer32(dim)));
    newsli : DoblDobl_Complex_VecVecs.VecVec(sli'range);
   -- eva : Complex_Number;

  begin
   -- put_line(file,"The solution vector : "); put_line(file,es.v);
   -- put_line(file,"The evaluation : "); put_line(file,y);
    put(file,"residual is "); put(file,res); new_line(file);
    if res <= restol then
      put(file,"  point satisfies the equations, as residual <=");
      put(file,restol,3); new_line(file);
      success := true;
      newsli := Adjusted_Slices(sli,es.v);
     -- for i in newsli'range loop
     --   put(file,"eval at adjusted slice : ");
     --   eva := newsli(i)(0) + newsli(i)(es.v'range)*es.v;
     --   put(file,eva); new_line(file);
     -- end loop;
      Homotopy_Membership_Test(file,genpts,es.v,newsli,homtol,found);
    else
      put(file,"  point does not lie on the component, as residual >");
      put(file,restol,3); new_line(file);
      success := false;
      found := false;
    end if;
  end Homotopy_Membership_Test;

  procedure Homotopy_Membership_Test 
              ( file : in file_type;
                ep : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                dim : in natural32;
                sli : in QuadDobl_Complex_VecVecs.VecVec;
                genpts : in QuadDobl_Complex_Solutions.Solution_List;
                s : in QuadDobl_Complex_Solutions.Solution;
                restol,homtol : in double_float;
                success,found : out boolean ) is

   -- es : Solution(s.n + dim) := Embed(s,dim);
    es : constant QuadDobl_Complex_Solutions.Solution(ep'last)
       := Embed_Solution(s,ep'last);
    y : constant QuadDobl_Complex_Vectors.Vector(ep'range)
      := QuadDobl_Complex_Poly_SysFun.Eval(ep,es.v);
    res : constant quad_double
        := Max_Norm(y(y'first..y'last-integer32(dim)));
    newsli : QuadDobl_Complex_VecVecs.VecVec(sli'range);
   -- eva : Complex_Number;

  begin
   -- put_line(file,"The solution vector : "); put_line(file,es.v);
   -- put_line(file,"The evaluation : "); put_line(file,y);
    put(file,"residual is "); put(file,res); new_line(file);
    if res <= restol then
      put(file,"  point satisfies the equations, as residual <=");
      put(file,restol,3); new_line(file);
      success := true;
      newsli := Adjusted_Slices(sli,es.v);
     -- for i in newsli'range loop
     --   put(file,"eval at adjusted slice : ");
     --   eva := newsli(i)(0) + newsli(i)(es.v'range)*es.v;
     --   put(file,eva); new_line(file);
     -- end loop;
      Homotopy_Membership_Test(file,genpts,es.v,newsli,homtol,found);
    else
      put(file,"  point does not lie on the component, as residual >");
      put(file,restol,3); new_line(file);
      success := false;
      found := false;
    end if;
  end Homotopy_Membership_Test;

  procedure Homotopy_Membership_Test 
              ( file : in file_type;
                ep : in Standard_Complex_Laur_Systems.Laur_Sys;
                dim : in natural32;
                sli : in Standard_Complex_VecVecs.VecVec;
                genpts : in Standard_Complex_Solutions.Solution_List;
                s : in Standard_Complex_Solutions.Solution;
                restol,homtol : in double_float;
                success,found : out boolean ) is

   -- es : Solution(s.n + dim) := Embed(s,dim);
    es : constant Standard_Complex_Solutions.Solution(ep'last)
       := Embed_Solution(s,ep'last);
    y : constant Standard_Complex_Vectors.Vector(ep'range)
      := Standard_Complex_Laur_SysFun.Eval(ep,es.v);
    res : constant double_float := Max_Norm(y(y'first..y'last-integer32(dim)));
    newsli : Standard_Complex_VecVecs.VecVec(sli'range);
   -- eva : Complex_Number;

  begin
   -- put_line(file,"The solution vector : "); put_line(file,es.v);
   -- put_line(file,"The evaluation : "); put_line(file,y);
    put(file,"residual is "); put(file,res); new_line(file);
    if res <= restol then
      put(file,"  point satisfies the equations, as residual <=");
      put(file,restol,3); new_line(file);
      success := true;
      newsli := Adjusted_Slices(sli,es.v);
     -- for i in newsli'range loop
     --   put(file,"eval at adjusted slice : ");
     --   eva := newsli(i)(0) + newsli(i)(es.v'range)*es.v;
     --   put(file,eva); new_line(file);
     -- end loop;
      Laurent_Homotopy_Membership_Test(file,genpts,es.v,newsli,homtol,found);
    else
      put(file,"  point does not lie on the component, as residual >");
      put(file,restol,3); new_line(file);
      success := false;
      found := false;
    end if;
  end Homotopy_Membership_Test;

  procedure Homotopy_Membership_Test 
              ( file : in file_type;
                ep : in DoblDobl_Complex_Laur_Systems.Laur_Sys;
                dim : in natural32;
                sli : in DoblDobl_Complex_VecVecs.VecVec;
                genpts : in DoblDobl_Complex_Solutions.Solution_List;
                s : in DoblDobl_Complex_Solutions.Solution;
                restol,homtol : in double_float;
                success,found : out boolean ) is

   -- es : Solution(s.n + dim) := Embed(s,dim);
    es : constant DoblDobl_Complex_Solutions.Solution(ep'last)
       := Embed_Solution(s,ep'last);
    y : constant DoblDobl_Complex_Vectors.Vector(ep'range)
      := DoblDobl_Complex_Laur_SysFun.Eval(ep,es.v);
    res : constant double_double
        := Max_Norm(y(y'first..y'last-integer32(dim)));
    newsli : DoblDobl_Complex_VecVecs.VecVec(sli'range);
   -- eva : Complex_Number;

  begin
   -- put_line(file,"The solution vector : "); put_line(file,es.v);
   -- put_line(file,"The evaluation : "); put_line(file,y);
    put(file,"residual is "); put(file,res); new_line(file);
    if res <= restol then
      put(file,"  point satisfies the equations, as residual <=");
      put(file,restol,3); new_line(file);
      success := true;
      newsli := Adjusted_Slices(sli,es.v);
     -- for i in newsli'range loop
     --   put(file,"eval at adjusted slice : ");
     --   eva := newsli(i)(0) + newsli(i)(es.v'range)*es.v;
     --   put(file,eva); new_line(file);
     -- end loop;
      Laurent_Homotopy_Membership_Test(file,genpts,es.v,newsli,homtol,found);
    else
      put(file,"  point does not lie on the component, as residual >");
      put(file,restol,3); new_line(file);
      success := false;
      found := false;
    end if;
  end Homotopy_Membership_Test;

  procedure Homotopy_Membership_Test 
              ( file : in file_type;
                ep : in QuadDobl_Complex_Laur_Systems.Laur_Sys;
                dim : in natural32;
                sli : in QuadDobl_Complex_VecVecs.VecVec;
                genpts : in QuadDobl_Complex_Solutions.Solution_List;
                s : in QuadDobl_Complex_Solutions.Solution;
                restol,homtol : in double_float;
                success,found : out boolean ) is

   -- es : Solution(s.n + dim) := Embed(s,dim);
    es : constant QuadDobl_Complex_Solutions.Solution(ep'last)
       := Embed_Solution(s,ep'last);
    y : constant QuadDobl_Complex_Vectors.Vector(ep'range)
      := QuadDobl_Complex_Laur_SysFun.Eval(ep,es.v);
    res : constant quad_double
        := Max_Norm(y(y'first..y'last-integer32(dim)));
    newsli : QuadDobl_Complex_VecVecs.VecVec(sli'range);
   -- eva : Complex_Number;

  begin
   -- put_line(file,"The solution vector : "); put_line(file,es.v);
   -- put_line(file,"The evaluation : "); put_line(file,y);
    put(file,"residual is "); put(file,res); new_line(file);
    if res <= restol then
      put(file,"  point satisfies the equations, as residual <=");
      put(file,restol,3); new_line(file);
      success := true;
      newsli := Adjusted_Slices(sli,es.v);
     -- for i in newsli'range loop
     --   put(file,"eval at adjusted slice : ");
     --   eva := newsli(i)(0) + newsli(i)(es.v'range)*es.v;
     --   put(file,eva); new_line(file);
     -- end loop;
      Laurent_Homotopy_Membership_Test(file,genpts,es.v,newsli,homtol,found);
    else
      put(file,"  point does not lie on the component, as residual >");
      put(file,restol,3); new_line(file);
      success := false;
      found := false;
    end if;
  end Homotopy_Membership_Test;

  procedure Homotopy_Membership_Test 
              ( verbose : in boolean;
                ep : in Standard_Complex_Poly_Systems.Poly_Sys;
                dim : in natural32;
                genpts,sols : in Standard_Complex_Solutions.Solution_List;
                restol,homtol : in double_float ) is

    tmp : Standard_Complex_Solutions.Solution_List := sols;
    ls : Standard_Complex_Solutions.Link_to_Solution;
    sli : constant Standard_Complex_VecVecs.VecVec
        := Witness_Sets.Slices(ep,dim);
    success,found : boolean;
    cnt : natural32 := 0;

  begin
    Sampling_Machine.Initialize(ep);
    Sampling_Machine.Default_Tune_Sampler(2);
    Sampling_Machine.Default_Tune_Refiner;
    for i in 1..Standard_Complex_Solutions.Length_Of(sols) loop
      ls := Standard_Complex_Solutions.Head_Of(tmp);
      if verbose
       then put("Testing point "); put(i,1); put(" : ");
      end if;
      Homotopy_Membership_Test
        (verbose,ep,dim,sli,genpts,ls.v,restol,homtol,success,found);
      if success then
        if found
         then cnt := cnt + 1;
        end if;
      end if;
      tmp := Standard_Complex_Solutions.Tail_Of(tmp);
    end loop;
    if verbose then
      put("Tested ");
      put(Standard_Complex_Solutions.Length_Of(sols),1); 
      put_line(" candidates for solutions:");
      put("  Found "); put(cnt,1);
      put_line(" solutions on the components.");
    end if;
    Sampling_Machine.Clear;
  end Homotopy_Membership_Test;

  procedure Homotopy_Membership_Test 
              ( verbose : in boolean;
                ep : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                dim : in natural32;
                genpts,sols : in DoblDobl_Complex_Solutions.Solution_List;
                restol,homtol : in double_float ) is

    tmp : DoblDobl_Complex_Solutions.Solution_List := sols;
    ls : DoblDobl_Complex_Solutions.Link_to_Solution;
    sli : constant DoblDobl_Complex_VecVecs.VecVec
        := Witness_Sets.Slices(ep,dim);
    success,found : boolean;
    cnt : natural32 := 0;

  begin
    DoblDobl_Sampling_Machine.Initialize(ep);
    DoblDobl_Sampling_Machine.Default_Tune_Sampler(0);
    DoblDobl_Sampling_Machine.Default_Tune_Refiner;
    for i in 1..DoblDobl_Complex_Solutions.Length_Of(sols) loop
      ls := DoblDobl_Complex_Solutions.Head_Of(tmp);
      if verbose
       then put("Testing point "); put(i,1); put(" : ");
      end if;
      Homotopy_Membership_Test
        (verbose,ep,dim,sli,genpts,ls.v,restol,homtol,success,found);
      if success then
        if found
         then cnt := cnt + 1;
        end if;
      end if;
      tmp := DoblDobl_Complex_Solutions.Tail_Of(tmp);
    end loop;
    if verbose then
      put("Tested ");
      put(DoblDobl_Complex_Solutions.Length_Of(sols),1); 
      put_line(" candidates for solutions:");
      put("  Found "); put(cnt,1);
      put_line(" solutions on the components.");
    end if;
    DoblDobl_Sampling_Machine.Clear;
  end Homotopy_Membership_Test;

  procedure Homotopy_Membership_Test 
              ( verbose : in boolean;
                ep : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                dim : in natural32;
                genpts,sols : in QuadDobl_Complex_Solutions.Solution_List;
                restol,homtol : in double_float ) is

    tmp : QuadDobl_Complex_Solutions.Solution_List := sols;
    ls : QuadDobl_Complex_Solutions.Link_to_Solution;
    sli : constant QuadDobl_Complex_VecVecs.VecVec
        := Witness_Sets.Slices(ep,dim);
    success,found : boolean;
    cnt : natural32 := 0;

  begin
    QuadDobl_Sampling_Machine.Initialize(ep);
    QuadDobl_Sampling_Machine.Default_Tune_Sampler(0);
    QuadDobl_Sampling_Machine.Default_Tune_Refiner;
    for i in 1..QuadDobl_Complex_Solutions.Length_Of(sols) loop
      ls := QuadDobl_Complex_Solutions.Head_Of(tmp);
      if verbose
       then put("Testing point "); put(i,1); put(" : ");
      end if;
      Homotopy_Membership_Test
        (verbose,ep,dim,sli,genpts,ls.v,restol,homtol,success,found);
      if success then
        if found
         then cnt := cnt + 1;
        end if;
      end if;
      tmp := QuadDobl_Complex_Solutions.Tail_Of(tmp);
    end loop;
    if verbose then
      put("Tested ");
      put(QuadDobl_Complex_Solutions.Length_Of(sols),1); 
      put_line(" candidates for solutions:");
      put("  Found "); put(cnt,1);
      put_line(" solutions on the components.");
    end if;
    QuadDobl_Sampling_Machine.Clear;
  end Homotopy_Membership_Test;

  procedure Homotopy_Membership_Test 
              ( verbose : in boolean;
                ep : in Standard_Complex_Laur_Systems.Laur_Sys;
                dim : in natural32;
                genpts,sols : in Standard_Complex_Solutions.Solution_List;
                restol,homtol : in double_float ) is

    tmp : Standard_Complex_Solutions.Solution_List := sols;
    ls : Standard_Complex_Solutions.Link_to_Solution;
    sli : constant Standard_Complex_VecVecs.VecVec
        := Witness_Sets.Slices(ep,dim);
    success,found : boolean;
    cnt : natural32 := 0;

  begin
    Sampling_Laurent_Machine.Initialize(ep);
    Sampling_Laurent_Machine.Default_Tune_Sampler(2);
    Sampling_Laurent_Machine.Default_Tune_Refiner;
    for i in 1..Standard_Complex_Solutions.Length_Of(sols) loop
      ls := Standard_Complex_Solutions.Head_Of(tmp);
      if verbose
       then put("Testing point "); put(i,1); put(" : ");
      end if;
      Homotopy_Membership_Test
        (verbose,ep,dim,sli,genpts,ls.v,restol,homtol,success,found);
      if success then
        if found
         then cnt := cnt + 1;
        end if;
      end if;
      tmp := Standard_Complex_Solutions.Tail_Of(tmp);
    end loop;
    if verbose then
      put("Tested ");
      put(Standard_Complex_Solutions.Length_Of(sols),1); 
      put_line(" candidates for solutions:");
      put("  Found "); put(cnt,1);
      put_line(" solutions on the components.");
    end if;
    Sampling_Laurent_Machine.Clear;
  end Homotopy_Membership_Test;

  procedure Homotopy_Membership_Test 
              ( verbose : in boolean;
                ep : in DoblDobl_Complex_Laur_Systems.Laur_Sys;
                dim : in natural32;
                genpts,sols : in DoblDobl_Complex_Solutions.Solution_List;
                restol,homtol : in double_float ) is

    tmp : DoblDobl_Complex_Solutions.Solution_List := sols;
    ls : DoblDobl_Complex_Solutions.Link_to_Solution;
    sli : constant DoblDobl_Complex_VecVecs.VecVec
        := Witness_Sets.Slices(ep,dim);
    success,found : boolean;
    cnt : natural32 := 0;

  begin
    DoblDobl_Sampling_Laurent_Machine.Initialize(ep);
    DoblDobl_Sampling_Laurent_Machine.Default_Tune_Sampler(0);
    DoblDobl_Sampling_Laurent_Machine.Default_Tune_Refiner;
    for i in 1..DoblDobl_Complex_Solutions.Length_Of(sols) loop
      ls := DoblDobl_Complex_Solutions.Head_Of(tmp);
      if verbose
       then put("Testing point "); put(i,1); put(" : ");
      end if;
      Homotopy_Membership_Test
        (verbose,ep,dim,sli,genpts,ls.v,restol,homtol,success,found);
      if success then
        if found
         then cnt := cnt + 1;
        end if;
      end if;
      tmp := DoblDobl_Complex_Solutions.Tail_Of(tmp);
    end loop;
    if verbose then
      put("Tested ");
      put(DoblDobl_Complex_Solutions.Length_Of(sols),1); 
      put_line(" candidates for solutions:");
      put("  Found "); put(cnt,1);
      put_line(" solutions on the components.");
    end if;
    DoblDobl_Sampling_Laurent_Machine.Clear;
  end Homotopy_Membership_Test;

  procedure Homotopy_Membership_Test 
              ( verbose : in boolean;
                ep : in QuadDobl_Complex_Laur_Systems.Laur_Sys;
                dim : in natural32;
                genpts,sols : in QuadDobl_Complex_Solutions.Solution_List;
                restol,homtol : in double_float ) is

    tmp : QuadDobl_Complex_Solutions.Solution_List := sols;
    ls : QuadDobl_Complex_Solutions.Link_to_Solution;
    sli : constant QuadDobl_Complex_VecVecs.VecVec
        := Witness_Sets.Slices(ep,dim);
    success,found : boolean;
    cnt : natural32 := 0;

  begin
    QuadDobl_Sampling_Laurent_Machine.Initialize(ep);
    QuadDobl_Sampling_Laurent_Machine.Default_Tune_Sampler(0);
    QuadDobl_Sampling_Laurent_Machine.Default_Tune_Refiner;
    for i in 1..QuadDobl_Complex_Solutions.Length_Of(sols) loop
      ls := QuadDobl_Complex_Solutions.Head_Of(tmp);
      if verbose
       then put("Testing point "); put(i,1); put(" : ");
      end if;
      Homotopy_Membership_Test
        (verbose,ep,dim,sli,genpts,ls.v,restol,homtol,success,found);
      if success then
        if found
         then cnt := cnt + 1;
        end if;
      end if;
      tmp := QuadDobl_Complex_Solutions.Tail_Of(tmp);
    end loop;
    if verbose then
      put("Tested ");
      put(QuadDobl_Complex_Solutions.Length_Of(sols),1); 
      put_line(" candidates for solutions:");
      put("  Found "); put(cnt,1);
      put_line(" solutions on the components.");
    end if;
    QuadDobl_Sampling_Laurent_Machine.Clear;
  end Homotopy_Membership_Test;

  procedure Homotopy_Membership_Test 
              ( file : in file_type;
                ep : in Standard_Complex_Poly_Systems.Poly_Sys;
                dim : in natural32;
                genpts,sols : in Standard_Complex_Solutions.Solution_List;
                restol,homtol : in double_float ) is

    isols,isols_last : Standard_Complex_Solutions.Solution_List;

  begin
    Homotopy_Membership_Test
      (file,ep,dim,genpts,sols,restol,homtol,isols,isols_last);
    Standard_Complex_Solutions.Clear(isols);
  end Homotopy_Membership_Test;

  procedure Homotopy_Membership_Test 
              ( file : in file_type;
                ep : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                dim : in natural32;
                genpts,sols : in DoblDobl_Complex_Solutions.Solution_List;
                restol,homtol : in double_float ) is

    isols,isols_last : DoblDobl_Complex_Solutions.Solution_List;

  begin
    Homotopy_Membership_Test
      (file,ep,dim,genpts,sols,restol,homtol,isols,isols_last);
    DoblDobl_Complex_Solutions.Clear(isols);
  end Homotopy_Membership_Test;

  procedure Homotopy_Membership_Test 
              ( file : in file_type;
                ep : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                dim : in natural32;
                genpts,sols : in QuadDobl_Complex_Solutions.Solution_List;
                restol,homtol : in double_float ) is

    isols,isols_last : QuadDobl_Complex_Solutions.Solution_List;

  begin
    Homotopy_Membership_Test
      (file,ep,dim,genpts,sols,restol,homtol,isols,isols_last);
    QuadDobl_Complex_Solutions.Clear(isols);
  end Homotopy_Membership_Test;

  procedure Homotopy_Membership_Test
              ( file : in file_type;
                ep : in Standard_Complex_Laur_Systems.Laur_Sys;
                dim : in natural32;
                genpts,sols : in Standard_Complex_Solutions.Solution_List;
                restol,homtol : in double_float ) is

    isols,isols_last : Standard_Complex_Solutions.Solution_List;

  begin
    Homotopy_Membership_Test
      (file,ep,dim,genpts,sols,restol,homtol,isols,isols_last);
    Standard_Complex_Solutions.Clear(isols);
  end Homotopy_Membership_Test;

  procedure Homotopy_Membership_Test
              ( file : in file_type;
                ep : in DoblDobl_Complex_Laur_Systems.Laur_Sys;
                dim : in natural32;
                genpts,sols : in DoblDobl_Complex_Solutions.Solution_List;
                restol,homtol : in double_float ) is

    isols,isols_last : DoblDobl_Complex_Solutions.Solution_List;

  begin
    Homotopy_Membership_Test
      (file,ep,dim,genpts,sols,restol,homtol,isols,isols_last);
    DoblDobl_Complex_Solutions.Clear(isols);
  end Homotopy_Membership_Test;

  procedure Homotopy_Membership_Test
              ( file : in file_type;
                ep : in QuadDobl_Complex_Laur_Systems.Laur_Sys;
                dim : in natural32;
                genpts,sols : in QuadDobl_Complex_Solutions.Solution_List;
                restol,homtol : in double_float ) is

    isols,isols_last : QuadDobl_Complex_Solutions.Solution_List;

  begin
    Homotopy_Membership_Test
      (file,ep,dim,genpts,sols,restol,homtol,isols,isols_last);
    QuadDobl_Complex_Solutions.Clear(isols);
  end Homotopy_Membership_Test;

  procedure Homotopy_Membership_Test 
              ( file : in file_type;
                ep : in Standard_Complex_Poly_Systems.Poly_Sys;
                dim : in natural32;
                genpts,sols : in Standard_Complex_Solutions.Solution_List;
                restol,homtol : in double_float;
                isols,isols_last
                  : in out Standard_Complex_Solutions.Solution_List ) is

    tmp : Standard_Complex_Solutions.Solution_List := sols;
    ls : Standard_Complex_Solutions.Link_to_Solution;
    sli : constant Standard_Complex_VecVecs.VecVec
        := Witness_Sets.Slices(ep,dim);
    success,found : boolean;
    cnt : natural32 := 0;

  begin
    Sampling_Machine.Initialize(ep);
    Sampling_Machine.Default_Tune_Sampler(2);
    Sampling_Machine.Default_Tune_Refiner;
    for i in 1..Standard_Complex_Solutions.Length_Of(sols) loop
      ls := Standard_Complex_Solutions.Head_Of(tmp);
      put(file,"Testing point "); put(file,i,1); put(file," : ");
      Homotopy_Membership_Test
        (file,ep,dim,sli,genpts,ls.all,restol,homtol,success,found);
      if success then
        if found
         then cnt := cnt + 1;
         else Standard_Complex_Solutions.Append(isols,isols_last,ls.all);
        end if;
      end if;
      tmp := Standard_Complex_Solutions.Tail_Of(tmp);
    end loop;
    Sampling_Machine.Clear;
    put(file,"Tested ");
    put(file,Standard_Complex_Solutions.Length_Of(sols),1); 
    put_line(file," candidates for solutions:");
    put(file,"  Found "); put(file,cnt,1);
    put_line(file," solutions on the components.");
    if not Standard_Complex_Solutions.Is_Null(isols) then
      put(file,"  Found ");
      put(file,Standard_Complex_Solutions.Length_Of(isols),1);
      put_line(file," solutions not on the components.");
      new_line(file);
      put_line(file,"THE SOLUTIONS NOT ON THE COMPONENTS :");
      new_line(file);
      put(file,Standard_Complex_Solutions.Length_Of(isols),
          natural32(Standard_Complex_Solutions.Head_Of(isols).n),isols);
    else
      put_line(file,"  Found no other solutions.");
    end if;
  end Homotopy_Membership_Test;

  procedure Homotopy_Membership_Test 
              ( file : in file_type;
                ep : in Standard_Complex_Laur_Systems.Laur_Sys;
                dim : in natural32;
                genpts,sols : in Standard_Complex_Solutions.Solution_List;
                restol,homtol : in double_float;
                isols,isols_last
                  : in out Standard_Complex_Solutions.Solution_List ) is

    tmp : Standard_Complex_Solutions.Solution_List := sols;
    ls : Standard_Complex_Solutions.Link_to_Solution;
    sli : constant Standard_Complex_VecVecs.VecVec
        := Witness_Sets.Slices(ep,dim);
    success,found : boolean;
    cnt : natural32 := 0;

  begin
    Sampling_Laurent_Machine.Initialize(ep);
    Sampling_Laurent_Machine.Default_Tune_Sampler(2);
    Sampling_Laurent_Machine.Default_Tune_Refiner;
    for i in 1..Standard_Complex_Solutions.Length_Of(sols) loop
      ls := Standard_Complex_Solutions.Head_Of(tmp);
      put(file,"Testing point "); put(file,i,1); put(file," : ");
      Homotopy_Membership_Test
        (file,ep,dim,sli,genpts,ls.all,restol,homtol,success,found);
      if success then
        if found
         then cnt := cnt + 1;
         else Standard_Complex_Solutions.Append(isols,isols_last,ls.all);
        end if;
      end if;
      tmp := Standard_Complex_Solutions.Tail_Of(tmp);
    end loop;
    Sampling_Laurent_Machine.Clear;
    put(file,"Tested ");
    put(file,Standard_Complex_Solutions.Length_Of(sols),1); 
    put_line(file," candidates for solutions:");
    put(file,"  Found "); put(file,cnt,1);
    put_line(file," solutions on the components.");
    if not Standard_Complex_Solutions.Is_Null(isols) then
      put(file,"  Found ");
      put(file,Standard_Complex_Solutions.Length_Of(isols),1);
      put_line(file," solutions not on the components.");
      new_line(file);
      put_line(file,"THE SOLUTIONS NOT ON THE COMPONENTS :");
      new_line(file);
      put(file,Standard_Complex_Solutions.Length_Of(isols),
          natural32(Standard_Complex_Solutions.Head_Of(isols).n),isols);
    else
      put_line(file,"  Found no other solutions.");
    end if;
  end Homotopy_Membership_Test;

  procedure Homotopy_Membership_Test 
              ( file : in file_type;
                ep : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                dim : in natural32;
                genpts,sols : in DoblDobl_Complex_Solutions.Solution_List;
                restol,homtol : in double_float;
                isols,isols_last
                  : in out DoblDobl_Complex_Solutions.Solution_List ) is

    tmp : DoblDobl_Complex_Solutions.Solution_List := sols;
    ls : DoblDobl_Complex_Solutions.Link_to_Solution;
    sli : constant DoblDobl_Complex_VecVecs.VecVec
        := Witness_Sets.Slices(ep,dim);
    success,found : boolean;
    cnt : natural32 := 0;

  begin
    DoblDobl_Sampling_Machine.Initialize(ep);
    DoblDobl_Sampling_Machine.Default_Tune_Sampler(0);
    DoblDobl_Sampling_Machine.Default_Tune_Refiner;
    for i in 1..DoblDobl_Complex_Solutions.Length_Of(sols) loop
      ls := DoblDobl_Complex_Solutions.Head_Of(tmp);
      put(file,"Testing point "); put(file,i,1); put(file," : ");
      Homotopy_Membership_Test
        (file,ep,dim,sli,genpts,ls.all,restol,homtol,success,found);
      if success then
        if found
         then cnt := cnt + 1;
         else DoblDobl_Complex_Solutions.Append(isols,isols_last,ls.all);
        end if;
      end if;
      tmp := DoblDobl_Complex_Solutions.Tail_Of(tmp);
    end loop;
    DoblDobl_Sampling_Machine.Clear;
    put(file,"Tested ");
    put(file,DoblDobl_Complex_Solutions.Length_Of(sols),1); 
    put_line(file," candidates for solutions:");
    put(file,"  Found "); put(file,cnt,1);
    put_line(file," solutions on the components.");
    if not DoblDobl_Complex_Solutions.Is_Null(isols) then
      put(file,"  Found ");
      put(file,DoblDobl_Complex_Solutions.Length_Of(isols),1);
      put_line(file," solutions not on the components.");
      new_line(file);
      put_line(file,"THE SOLUTIONS NOT ON THE COMPONENTS :");
      new_line(file);
      put(file,DoblDobl_Complex_Solutions.Length_Of(isols),
          natural32(DoblDobl_Complex_Solutions.Head_Of(isols).n),isols);
    else
      put_line(file,"  Found no other solutions.");
    end if;
  end Homotopy_Membership_Test;

  procedure Homotopy_Membership_Test 
              ( file : in file_type;
                ep : in DoblDobl_Complex_Laur_Systems.Laur_Sys;
                dim : in natural32;
                genpts,sols : in DoblDobl_Complex_Solutions.Solution_List;
                restol,homtol : in double_float;
                isols,isols_last
                  : in out DoblDobl_Complex_Solutions.Solution_List ) is

    tmp : DoblDobl_Complex_Solutions.Solution_List := sols;
    ls : DoblDobl_Complex_Solutions.Link_to_Solution;
    sli : constant DoblDobl_Complex_VecVecs.VecVec
        := Witness_Sets.Slices(ep,dim);
    success,found : boolean;
    cnt : natural32 := 0;

  begin
    DoblDobl_Sampling_Laurent_Machine.Initialize(ep);
    DoblDobl_Sampling_Laurent_Machine.Default_Tune_Sampler(0);
    DoblDobl_Sampling_Laurent_Machine.Default_Tune_Refiner;
    for i in 1..DoblDobl_Complex_Solutions.Length_Of(sols) loop
      ls := DoblDobl_Complex_Solutions.Head_Of(tmp);
      put(file,"Testing point "); put(file,i,1); put(file," : ");
      Homotopy_Membership_Test
        (file,ep,dim,sli,genpts,ls.all,restol,homtol,success,found);
      if success then
        if found
         then cnt := cnt + 1;
         else DoblDobl_Complex_Solutions.Append(isols,isols_last,ls.all);
        end if;
      end if;
      tmp := DoblDobl_Complex_Solutions.Tail_Of(tmp);
    end loop;
    DoblDobl_Sampling_Laurent_Machine.Clear;
    put(file,"Tested ");
    put(file,DoblDobl_Complex_Solutions.Length_Of(sols),1); 
    put_line(file," candidates for solutions:");
    put(file,"  Found "); put(file,cnt,1);
    put_line(file," solutions on the components.");
    if not DoblDobl_Complex_Solutions.Is_Null(isols) then
      put(file,"  Found ");
      put(file,DoblDobl_Complex_Solutions.Length_Of(isols),1);
      put_line(file," solutions not on the components.");
      new_line(file);
      put_line(file,"THE SOLUTIONS NOT ON THE COMPONENTS :");
      new_line(file);
      put(file,DoblDobl_Complex_Solutions.Length_Of(isols),
          natural32(DoblDobl_Complex_Solutions.Head_Of(isols).n),isols);
    else
      put_line(file,"  Found no other solutions.");
    end if;
  end Homotopy_Membership_Test;

  procedure Homotopy_Membership_Test 
              ( file : in file_type;
                ep : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                dim : in natural32;
                genpts,sols : in QuadDobl_Complex_Solutions.Solution_List;
                restol,homtol : in double_float;
                isols,isols_last
                  : in out QuadDobl_Complex_Solutions.Solution_List ) is

    tmp : QuadDobl_Complex_Solutions.Solution_List := sols;
    ls : QuadDobl_Complex_Solutions.Link_to_Solution;
    sli : constant QuadDobl_Complex_VecVecs.VecVec
        := Witness_Sets.Slices(ep,dim);
    success,found : boolean;
    cnt : natural32 := 0;

  begin
    QuadDobl_Sampling_Machine.Initialize(ep);
    QuadDobl_Sampling_Machine.Default_Tune_Sampler(0);
    QuadDobl_Sampling_Machine.Default_Tune_Refiner;
    for i in 1..QuadDobl_Complex_Solutions.Length_Of(sols) loop
      ls := QuadDobl_Complex_Solutions.Head_Of(tmp);
      put(file,"Testing point "); put(file,i,1); put(file," : ");
      Homotopy_Membership_Test
        (file,ep,dim,sli,genpts,ls.all,restol,homtol,success,found);
      if success then
        if found
         then cnt := cnt + 1;
         else QuadDobl_Complex_Solutions.Append(isols,isols_last,ls.all);
        end if;
      end if;
      tmp := QuadDobl_Complex_Solutions.Tail_Of(tmp);
    end loop;
    QuadDobl_Sampling_Machine.Clear;
    put(file,"Tested ");
    put(file,QuadDobl_Complex_Solutions.Length_Of(sols),1); 
    put_line(file," candidates for solutions:");
    put(file,"  Found "); put(file,cnt,1);
    put_line(file," solutions on the components.");
    if not QuadDobl_Complex_Solutions.Is_Null(isols) then
      put(file,"  Found ");
      put(file,QuadDobl_Complex_Solutions.Length_Of(isols),1);
      put_line(file," solutions not on the components.");
      new_line(file);
      put_line(file,"THE SOLUTIONS NOT ON THE COMPONENTS :");
      new_line(file);
      put(file,QuadDobl_Complex_Solutions.Length_Of(isols),
          natural32(QuadDobl_Complex_Solutions.Head_Of(isols).n),isols);
    else
      put_line(file,"  Found no other solutions.");
    end if;
  end Homotopy_Membership_Test;

  procedure Homotopy_Membership_Test 
              ( file : in file_type;
                ep : in QuadDobl_Complex_Laur_Systems.Laur_Sys;
                dim : in natural32;
                genpts,sols : in QuadDobl_Complex_Solutions.Solution_List;
                restol,homtol : in double_float;
                isols,isols_last
                  : in out QuadDobl_Complex_Solutions.Solution_List ) is

    tmp : QuadDobl_Complex_Solutions.Solution_List := sols;
    ls : QuadDobl_Complex_Solutions.Link_to_Solution;
    sli : constant QuadDobl_Complex_VecVecs.VecVec
        := Witness_Sets.Slices(ep,dim);
    success,found : boolean;
    cnt : natural32 := 0;

  begin
    QuadDobl_Sampling_Laurent_Machine.Initialize(ep);
    QuadDobl_Sampling_Laurent_Machine.Default_Tune_Sampler(0);
    QuadDobl_Sampling_Laurent_Machine.Default_Tune_Refiner;
    for i in 1..QuadDobl_Complex_Solutions.Length_Of(sols) loop
      ls := QuadDobl_Complex_Solutions.Head_Of(tmp);
      put(file,"Testing point "); put(file,i,1); put(file," : ");
      Homotopy_Membership_Test
        (file,ep,dim,sli,genpts,ls.all,restol,homtol,success,found);
      if success then
        if found
         then cnt := cnt + 1;
         else QuadDobl_Complex_Solutions.Append(isols,isols_last,ls.all);
        end if;
      end if;
      tmp := QuadDobl_Complex_Solutions.Tail_Of(tmp);
    end loop;
    QuadDobl_Sampling_Laurent_Machine.Clear;
    put(file,"Tested ");
    put(file,QuadDobl_Complex_Solutions.Length_Of(sols),1); 
    put_line(file," candidates for solutions:");
    put(file,"  Found "); put(file,cnt,1);
    put_line(file," solutions on the components.");
    if not QuadDobl_Complex_Solutions.Is_Null(isols) then
      put(file,"  Found ");
      put(file,QuadDobl_Complex_Solutions.Length_Of(isols),1);
      put_line(file," solutions not on the components.");
      new_line(file);
      put_line(file,"THE SOLUTIONS NOT ON THE COMPONENTS :");
      new_line(file);
      put(file,QuadDobl_Complex_Solutions.Length_Of(isols),
          natural32(QuadDobl_Complex_Solutions.Head_Of(isols).n),isols);
    else
      put_line(file,"  Found no other solutions.");
    end if;
  end Homotopy_Membership_Test;

end Homotopy_Membership_Tests;
