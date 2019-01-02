with text_io;                            use text_io;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;
with Double_Double_Numbers;              use Double_Double_Numbers;
with Double_Double_Numbers_io;           use Double_Double_Numbers_io;
with Quad_Double_Numbers;                use Quad_Double_Numbers;
with Quad_Double_Numbers_io;             use Quad_Double_Numbers_io;
with Standard_Complex_Numbers_io;        use Standard_Complex_Numbers_io;
with Standard_Integer_Vectors_io;        use Standard_Integer_Vectors_io;
with Standard_Complex_Vectors_io;        use Standard_Complex_Vectors_io;
with Standard_Complex_Norms_Equals;
with DoblDobl_Complex_Vectors_io;        use DoblDobl_Complex_Vectors_io;
with DoblDobl_Complex_Vector_Norms;
with QuadDobl_Complex_Vectors_io;        use QuadDobl_Complex_Vectors_io;
with QuadDobl_Complex_Vector_Norms;
with Standard_Integer_Matrices_io;       use Standard_Integer_Matrices_io;
with Standard_Integer64_Matrices_io;     use Standard_Integer64_Matrices_io;
with Standard_Complex_Matrices_io;       use Standard_Complex_Matrices_io;
with DoblDobl_Complex_Matrices_io;       use DoblDobl_Complex_Matrices_io;
with QuadDobl_Complex_Matrices_io;       use QuadDobl_Complex_Matrices_io;
with Standard_Integer_Linear_Solvers;
with Standard_Integer64_Linear_Solvers;
with Standard_Complex_Linear_Solvers;
with DoblDobl_Complex_Linear_Solvers;
with QuadDobl_Complex_Linear_Solvers;
with Lists_of_Floating_Vectors;          use Lists_of_Floating_Vectors;
with Standard_Complex_Laur_Functions;
with Standard_Complex_Laur_Systems_io;   use Standard_Complex_Laur_Systems_io;
with Standard_Complex_Laur_SysFun;
with Standard_Tableau_Formats;
with DoblDobl_Complex_Laur_Functions;
with DoblDobl_Complex_Laur_Systems_io;   use DoblDobl_Complex_Laur_Systems_io;
with DoblDobl_Complex_Laur_SysFun;
with DoblDobl_Tableau_Formats;
with QuadDobl_Complex_Laur_Functions;
with QuadDobl_Complex_Laur_Systems_io;   use QuadDobl_Complex_Laur_Systems_io;
with QuadDobl_Complex_Laur_SysFun;
with QuadDobl_Tableau_Formats;
with Supports_of_Polynomial_Systems;
--with Mixed_Volume_Computation;
with Standard_Radial_Solvers;
with Standard_Binomial_Systems;
with Standard_Binomial_Solvers;
with DoblDobl_Radial_Solvers;
with DoblDobl_Binomial_Systems;
with DoblDobl_Binomial_Solvers;
with QuadDobl_Radial_Solvers;
with QuadDobl_Binomial_Systems;
with QuadDobl_Binomial_Solvers;

package body Polyhedral_Start_Systems is

-- (1) TABLEAU DATA STRUCTURES FOR START SYSTEM SELECTION :

  function Is_Equal
             ( x : Standard_Integer_Vectors.Link_to_Vector;
               y : Standard_Floating_Vectors.Link_to_Vector )
             return boolean is
  begin
    for i in x'range loop
      if x(i) /= integer32(y(i))
       then return false;
      end if;
    end loop;
    return true;
  end Is_Equal;

  function Coefficient
              ( cff : Standard_Complex_Vectors.Link_to_Vector;
                exp : Standard_Integer_VecVecs.Link_to_VecVec;
                pt : Standard_Floating_Vectors.Link_to_Vector )
              return Standard_Complex_Numbers.Complex_Number is

    use Standard_Complex_Numbers;

    res : Complex_Number := Create(0.0);

  begin
    for i in exp'range loop
      if Is_Equal(exp(i),pt)
       then res := cff(i); exit;
      end if;
    end loop;
    return res;
  end Coefficient;

  function Coefficient
              ( cff : DoblDobl_Complex_Vectors.Link_to_Vector;
                exp : Standard_Integer_VecVecs.Link_to_VecVec;
                pt : Standard_Floating_Vectors.Link_to_Vector )
              return DoblDobl_Complex_Numbers.Complex_Number is

    use DoblDobl_Complex_Numbers;

    zero : constant double_double := create(0.0);
    res : Complex_Number := Create(zero);

  begin
    for i in exp'range loop
      if Is_Equal(exp(i),pt)
       then res := cff(i); exit;
      end if;
    end loop;
    return res;
  end Coefficient;

  function Coefficient
              ( cff : QuadDobl_Complex_Vectors.Link_to_Vector;
                exp : Standard_Integer_VecVecs.Link_to_VecVec;
                pt : Standard_Floating_Vectors.Link_to_Vector )
              return QuadDobl_Complex_Numbers.Complex_Number is

    use QuadDobl_Complex_Numbers;

    zero : constant quad_double := create(0.0);
    res : Complex_Number := Create(zero);

  begin
    for i in exp'range loop
      if Is_Equal(exp(i),pt)
       then res := cff(i); exit;
      end if;
    end loop;
    return res;
  end Coefficient;

  procedure Select_Coefficients
              ( q_c : in Standard_Complex_VecVecs.VecVec;
                q_e : in Standard_Integer_VecVecs.Array_of_VecVecs;
                pts : in Arrays_of_Floating_Vector_Lists.Array_of_Lists;
                cff : out Standard_Complex_Vectors.Vector ) is

    ind : integer32 := cff'first-1;
    tmp : List;
    lpt : Standard_Floating_Vectors.Link_to_Vector;

  begin
    for i in pts'range loop
      tmp := pts(i);
      while not Is_Null(tmp) loop
        lpt := Head_Of(tmp);
        ind := ind + 1;
        cff(ind) := Coefficient(q_c(i),q_e(i),lpt);
        tmp := Tail_Of(tmp);
      end loop;
    end loop;
  end Select_Coefficients;

  procedure Select_Coefficients
              ( q_c : in Standard_Complex_VecVecs.VecVec;
                q_e : in Exponent_Vectors.Exponent_Vectors_Array;
                pts : in Arrays_of_Floating_Vector_Lists.Array_of_Lists;
                cff : out Standard_Complex_Vectors.Vector ) is

    ind : integer32 := cff'first-1;
    tmp : List;
    lpt : Standard_Floating_Vectors.Link_to_Vector;

  begin
    for i in pts'range loop
      tmp := pts(i);
      while not Is_Null(tmp) loop
        lpt := Head_Of(tmp);
        ind := ind + 1;
        cff(ind) := Coefficient(q_c(i),q_e(i),lpt);
        tmp := Tail_Of(tmp);
      end loop;
    end loop;
  end Select_Coefficients;

  procedure Select_Coefficients
              ( q_c : in DoblDobl_Complex_VecVecs.VecVec;
                q_e : in Standard_Integer_VecVecs.Array_of_VecVecs;
                pts : in Arrays_of_Floating_Vector_Lists.Array_of_Lists;
                cff : out DoblDobl_Complex_Vectors.Vector ) is

    ind : integer32 := cff'first-1;
    tmp : List;
    lpt : Standard_Floating_Vectors.Link_to_Vector;

  begin
    for i in pts'range loop
      tmp := pts(i);
      while not Is_Null(tmp) loop
        lpt := Head_Of(tmp);
        ind := ind + 1;
        cff(ind) := Coefficient(q_c(i),q_e(i),lpt);
        tmp := Tail_Of(tmp);
      end loop;
    end loop;
  end Select_Coefficients;

  procedure Select_Coefficients
              ( q_c : in DoblDobl_Complex_VecVecs.VecVec;
                q_e : in Exponent_Vectors.Exponent_Vectors_Array;
                pts : in Arrays_of_Floating_Vector_Lists.Array_of_Lists;
                cff : out DoblDobl_Complex_Vectors.Vector ) is

    ind : integer32 := cff'first-1;
    tmp : List;
    lpt : Standard_Floating_Vectors.Link_to_Vector;

  begin
    for i in pts'range loop
      tmp := pts(i);
      while not Is_Null(tmp) loop
        lpt := Head_Of(tmp);
        ind := ind + 1;
        cff(ind) := Coefficient(q_c(i),q_e(i),lpt);
        tmp := Tail_Of(tmp);
      end loop;
    end loop;
  end Select_Coefficients;

  procedure Select_Coefficients
              ( q_c : in QuadDobl_Complex_VecVecs.VecVec;
                q_e : in Standard_Integer_VecVecs.Array_of_VecVecs;
                pts : in Arrays_of_Floating_Vector_Lists.Array_of_Lists;
                cff : out QuadDobl_Complex_Vectors.Vector ) is

    ind : integer32 := cff'first-1;
    tmp : List;
    lpt : Standard_Floating_Vectors.Link_to_Vector;

  begin
    for i in pts'range loop
      tmp := pts(i);
      while not Is_Null(tmp) loop
        lpt := Head_Of(tmp);
        ind := ind + 1;
        cff(ind) := Coefficient(q_c(i),q_e(i),lpt);
        tmp := Tail_Of(tmp);
      end loop;
    end loop;
  end Select_Coefficients;

  procedure Select_Coefficients
              ( q_c : in QuadDobl_Complex_VecVecs.VecVec;
                q_e : in Exponent_Vectors.Exponent_Vectors_Array;
                pts : in Arrays_of_Floating_Vector_Lists.Array_of_Lists;
                cff : out QuadDobl_Complex_Vectors.Vector ) is

    ind : integer32 := cff'first-1;
    tmp : List;
    lpt : Standard_Floating_Vectors.Link_to_Vector;

  begin
    for i in pts'range loop
      tmp := pts(i);
      while not Is_Null(tmp) loop
        lpt := Head_Of(tmp);
        ind := ind + 1;
        cff(ind) := Coefficient(q_c(i),q_e(i),lpt);
        tmp := Tail_Of(tmp);
      end loop;
    end loop;
  end Select_Coefficients;

  procedure Write_Tableau
              ( c : in Standard_Complex_VecVecs.VecVec;
                e : in Standard_Integer_VecVecs.Array_of_VecVecs ) is
  begin
    put(c'last,1); new_line;
    for i in c'range loop
      put(c(i)'last,1); new_line;
      for j in c(i)'range loop
        put(c(i)(j)); put("  ");
        put(e(i)(j)); new_line;
      end loop;
    end loop;
  end Write_Tableau;

  procedure Write_Tableau
              ( c : in Standard_Complex_Vectors.Vector;
                e : in Arrays_of_Floating_Vector_Lists.Array_of_Lists ) is

    ind : integer32 := c'first - 1;
    tmp : List;
    lpt : Standard_Floating_Vectors.Link_to_Vector;

  begin
    put(e'last,1); new_line;
    for i in e'range loop
      put(Length_Of(e(i)),1); new_line;
      tmp := e(i);
      while not Is_Null(tmp) loop
        lpt := Head_Of(tmp);
        ind := ind + 1;
        put(c(ind)); put("  ");
        for k in lpt'first..lpt'last-1 loop -- drop lifting
          put(" "); put(integer32(lpt(k)),1);
        end loop;
        new_line;
        tmp := Tail_Of(tmp);
      end loop;
    end loop;
  end Write_Tableau;

  procedure Fully_Mixed_to_Binomial_Format
              ( c : in Standard_Complex_Vectors.Vector;
                e : in Arrays_of_Floating_Vector_Lists.Array_of_Lists;
                A : out Standard_Integer_Matrices.Matrix;
                b : out Standard_Complex_Vectors.Vector ) is

    use Standard_Complex_Numbers;

    first,second : Standard_Floating_Vectors.Link_to_Vector;
    ind : integer32 := c'first;

  begin
    for i in e'range loop
      first := Head_Of(e(i));
      second := Head_Of(Tail_Of(e(i)));
      for j in A'range(2) loop
        A(j,i) := integer32(first(j)) - integer32(second(j));
      end loop;
      b(i) := (-c(ind+1))/c(ind);
      ind := ind + 2;
    end loop;
  end Fully_Mixed_to_Binomial_Format;

  procedure Fully_Mixed_to_Binomial_Format
              ( c : in DoblDobl_Complex_Vectors.Vector;
                e : in Arrays_of_Floating_Vector_Lists.Array_of_Lists;
                A : out Standard_Integer64_Matrices.Matrix;
                b : out DoblDobl_Complex_Vectors.Vector ) is

    use DoblDobl_Complex_Numbers;

    first,second : Standard_Floating_Vectors.Link_to_Vector;
    ind : integer32 := c'first;

  begin
    for i in e'range loop
      first := Head_Of(e(i));
      second := Head_Of(Tail_Of(e(i)));
      for j in A'range(2) loop
        A(j,i) := integer64(first(j)) - integer64(second(j));
      end loop;
      b(i) := (-c(ind+1))/c(ind);
      ind := ind + 2;
    end loop;
  end Fully_Mixed_to_Binomial_Format;

  procedure Fully_Mixed_to_Binomial_Format
              ( c : in QuadDobl_Complex_Vectors.Vector;
                e : in Arrays_of_Floating_Vector_Lists.Array_of_Lists;
                A : out Standard_Integer64_Matrices.Matrix;
                b : out QuadDobl_Complex_Vectors.Vector ) is

    use QuadDobl_Complex_Numbers;

    first,second : Standard_Floating_Vectors.Link_to_Vector;
    ind : integer32 := c'first;

  begin
    for i in e'range loop
      first := Head_Of(e(i));
      second := Head_Of(Tail_Of(e(i)));
      for j in A'range(2) loop
        A(j,i) := integer64(first(j)) - integer64(second(j));
      end loop;
      b(i) := (-c(ind+1))/c(ind);
      ind := ind + 2;
    end loop;
  end Fully_Mixed_to_Binomial_Format;

  procedure Select_Subsystem_to_Matrix_Format 
              ( q_c : in Standard_Complex_VecVecs.VecVec;
                q_e : in Standard_Integer_VecVecs.Array_of_VecVecs;
                mix : in Standard_Integer_Vectors.Vector;
                pts : in Arrays_of_Floating_Vector_Lists.Array_of_Lists;
                A : out Standard_Integer_Matrices.Matrix;
                C : out Standard_Complex_Matrices.Matrix;
                b : out Standard_Complex_Vectors.Vector ) is

    use Standard_Complex_Numbers;

    ind : integer32 := q_c'first-1;
    first,second : Standard_Floating_Vectors.Link_to_Vector;
    tmp : List;
    col : integer32 := A'first(2) - 1;

  begin
    for i in C'range(1) loop
      for j in C'range(2) loop
        C(i,j) := Create(0.0);
      end loop;
    end loop;
    for i in mix'range loop
      first := Head_Of(pts(i));
      for k in 1..mix(i) loop
        b(ind+k) := -Coefficient(q_c(ind+k),q_e(ind+k),first);
      end loop;
      tmp := Tail_Of(pts(i));
      while not Is_Null(tmp) loop
        second := Head_Of(tmp);
        col := col + 1;
        for row in A'range(1) loop
          A(row,col) := integer32(second(row)) - integer32(first(row));
        end loop;
        for k in 1..mix(i) loop
          C(ind+k,col) := Coefficient(q_c(ind+k),q_e(ind+k),second);
        end loop;
        tmp := Tail_Of(tmp);
      end loop;
      ind := ind + mix(i);
    end loop;
  end Select_Subsystem_to_Matrix_Format;

  procedure Select_Subsystem_to_Matrix_Format 
              ( q_c : in Standard_Complex_VecVecs.VecVec;
                q_e : in Exponent_Vectors.Exponent_Vectors_Array;
                mix : in Standard_Integer_Vectors.Vector;
                pts : in Arrays_of_Floating_Vector_Lists.Array_of_Lists;
                A : out Standard_Integer_Matrices.Matrix;
                C : out Standard_Complex_Matrices.Matrix;
                b : out Standard_Complex_Vectors.Vector ) is

    use Standard_Complex_Numbers;

    ind : integer32 := q_c'first-1;
    first,second : Standard_Floating_Vectors.Link_to_Vector;
    tmp : List;
    col : integer32 := A'first(2) - 1;

  begin
    for i in C'range(1) loop
      for j in C'range(2) loop
        C(i,j) := Create(0.0);
      end loop;
    end loop;
    for i in mix'range loop
      first := Head_Of(pts(i));
      for k in 1..mix(i) loop
        b(ind+k) := -Coefficient(q_c(ind+k),q_e(ind+k),first);
      end loop;
      tmp := Tail_Of(pts(i));
      while not Is_Null(tmp) loop
        second := Head_Of(tmp);
        col := col + 1;
        for row in A'range(1) loop
          A(row,col) := integer32(second(row)) - integer32(first(row));
        end loop;
        for k in 1..mix(i) loop
          C(ind+k,col) := Coefficient(q_c(ind+k),q_e(ind+k),second);
        end loop;
        tmp := Tail_Of(tmp);
      end loop;
      ind := ind + mix(i);
    end loop;
  end Select_Subsystem_to_Matrix_Format;

  procedure Select_Subsystem_to_Matrix_Format 
              ( q_c : in DoblDobl_Complex_VecVecs.VecVec;
                q_e : in Standard_Integer_VecVecs.Array_of_VecVecs;
                mix : in Standard_Integer_Vectors.Vector;
                pts : in Arrays_of_Floating_Vector_Lists.Array_of_Lists;
                A : out Standard_Integer64_Matrices.Matrix;
                C : out DoblDobl_Complex_Matrices.Matrix;
                b : out DoblDobl_Complex_Vectors.Vector ) is

    use DoblDobl_Complex_Numbers;

    ind : integer32 := q_c'first-1;
    first,second : Standard_Floating_Vectors.Link_to_Vector;
    tmp : List;
    col : integer32 := A'first(2) - 1;
    zero : constant double_double := create(0.0);

  begin
    for i in C'range(1) loop
      for j in C'range(2) loop
        C(i,j) := Create(zero);
      end loop;
    end loop;
    for i in mix'range loop
      first := Head_Of(pts(i));
      for k in 1..mix(i) loop
        b(ind+k) := -Coefficient(q_c(ind+k),q_e(ind+k),first);
      end loop;
      tmp := Tail_Of(pts(i));
      while not Is_Null(tmp) loop
        second := Head_Of(tmp);
        col := col + 1;
        for row in A'range(1) loop
          A(row,col) := integer64(second(row)) - integer64(first(row));
        end loop;
        for k in 1..mix(i) loop
          C(ind+k,col) := Coefficient(q_c(ind+k),q_e(ind+k),second);
        end loop;
        tmp := Tail_Of(tmp);
      end loop;
      ind := ind + mix(i);
    end loop;
  end Select_Subsystem_to_Matrix_Format;

  procedure Select_Subsystem_to_Matrix_Format 
              ( q_c : in DoblDobl_Complex_VecVecs.VecVec;
                q_e : in Exponent_Vectors.Exponent_Vectors_Array;
                mix : in Standard_Integer_Vectors.Vector;
                pts : in Arrays_of_Floating_Vector_Lists.Array_of_Lists;
                A : out Standard_Integer64_Matrices.Matrix;
                C : out DoblDobl_Complex_Matrices.Matrix;
                b : out DoblDobl_Complex_Vectors.Vector ) is

    use DoblDobl_Complex_Numbers;

    ind : integer32 := q_c'first-1;
    first,second : Standard_Floating_Vectors.Link_to_Vector;
    tmp : List;
    col : integer32 := A'first(2) - 1;
    zero : constant double_double := create(0.0);

  begin
    for i in C'range(1) loop
      for j in C'range(2) loop
        C(i,j) := Create(zero);
      end loop;
    end loop;
    for i in mix'range loop
      first := Head_Of(pts(i));
      for k in 1..mix(i) loop
        b(ind+k) := -Coefficient(q_c(ind+k),q_e(ind+k),first);
      end loop;
      tmp := Tail_Of(pts(i));
      while not Is_Null(tmp) loop
        second := Head_Of(tmp);
        col := col + 1;
        for row in A'range(1) loop
          A(row,col) := integer64(second(row)) - integer64(first(row));
        end loop;
        for k in 1..mix(i) loop
          C(ind+k,col) := Coefficient(q_c(ind+k),q_e(ind+k),second);
        end loop;
        tmp := Tail_Of(tmp);
      end loop;
      ind := ind + mix(i);
    end loop;
  end Select_Subsystem_to_Matrix_Format;

  procedure Select_Subsystem_to_Matrix_Format 
              ( q_c : in QuadDobl_Complex_VecVecs.VecVec;
                q_e : in Standard_Integer_VecVecs.Array_of_VecVecs;
                mix : in Standard_Integer_Vectors.Vector;
                pts : in Arrays_of_Floating_Vector_Lists.Array_of_Lists;
                A : out Standard_Integer64_Matrices.Matrix;
                C : out QuadDobl_Complex_Matrices.Matrix;
                b : out QuadDobl_Complex_Vectors.Vector ) is

    use QuadDobl_Complex_Numbers;

    ind : integer32 := q_c'first-1;
    first,second : Standard_Floating_Vectors.Link_to_Vector;
    tmp : List;
    col : integer32 := A'first(2) - 1;
    zero : constant quad_double := create(0.0);

  begin
    for i in C'range(1) loop
      for j in C'range(2) loop
        C(i,j) := Create(zero);
      end loop;
    end loop;
    for i in mix'range loop
      first := Head_Of(pts(i));
      for k in 1..mix(i) loop
        b(ind+k) := -Coefficient(q_c(ind+k),q_e(ind+k),first);
      end loop;
      tmp := Tail_Of(pts(i));
      while not Is_Null(tmp) loop
        second := Head_Of(tmp);
        col := col + 1;
        for row in A'range(1) loop
          A(row,col) := integer64(second(row)) - integer64(first(row));
        end loop;
        for k in 1..mix(i) loop
          C(ind+k,col) := Coefficient(q_c(ind+k),q_e(ind+k),second);
        end loop;
        tmp := Tail_Of(tmp);
      end loop;
      ind := ind + mix(i);
    end loop;
  end Select_Subsystem_to_Matrix_Format;

  procedure Select_Subsystem_to_Matrix_Format 
              ( q_c : in QuadDobl_Complex_VecVecs.VecVec;
                q_e : in Exponent_Vectors.Exponent_Vectors_Array;
                mix : in Standard_Integer_Vectors.Vector;
                pts : in Arrays_of_Floating_Vector_Lists.Array_of_Lists;
                A : out Standard_Integer64_Matrices.Matrix;
                C : out QuadDobl_Complex_Matrices.Matrix;
                b : out QuadDobl_Complex_Vectors.Vector ) is

    use QuadDobl_Complex_Numbers;

    ind : integer32 := q_c'first-1;
    first,second : Standard_Floating_Vectors.Link_to_Vector;
    tmp : List;
    col : integer32 := A'first(2) - 1;
    zero : constant quad_double := create(0.0);

  begin
    for i in C'range(1) loop
      for j in C'range(2) loop
        C(i,j) := Create(zero);
      end loop;
    end loop;
    for i in mix'range loop
      first := Head_Of(pts(i));
      for k in 1..mix(i) loop
        b(ind+k) := -Coefficient(q_c(ind+k),q_e(ind+k),first);
      end loop;
      tmp := Tail_Of(pts(i));
      while not Is_Null(tmp) loop
        second := Head_Of(tmp);
        col := col + 1;
        for row in A'range(1) loop
          A(row,col) := integer64(second(row)) - integer64(first(row));
        end loop;
        for k in 1..mix(i) loop
          C(ind+k,col) := Coefficient(q_c(ind+k),q_e(ind+k),second);
        end loop;
        tmp := Tail_Of(tmp);
      end loop;
      ind := ind + mix(i);
    end loop;
  end Select_Subsystem_to_Matrix_Format;

-- (2) INPLACE BINOMIAL SYSTEM SOLVERS FOR START SOLUTIONS :

  function Create ( n : integer32 )
                  return Standard_Complex_Solutions.Solution is

    res : Standard_Complex_Solutions.Solution(n);

  begin
    res.t := Standard_Complex_Numbers.Create(0.0);
    res.m := 1;
    res.v := (1..n => Standard_Complex_Numbers.Create(0.0));
    res.err := 0.0;
    res.rco := 1.0;
    res.res := 0.0;
    return res;
  end Create;

  function Create ( n,m : integer32 )
                  return Standard_Complex_Solutions.Solution_List is

    use Standard_Complex_Solutions;
    res,res_last : Solution_List;

  begin
    for i in 1..m loop
      Append(res,res_last,Create(n));
    end loop;
    return res;
  end Create;

  function Create ( n : integer32 )
                  return DoblDobl_Complex_Solutions.Solution is

    res : DoblDobl_Complex_Solutions.Solution(n);
    zero : constant double_double := create(0.0);
    one : constant double_double := create(1.0);

  begin
    res.t := DoblDobl_Complex_Numbers.Create(zero);
    res.m := 1;
    res.v := (1..n => DoblDobl_Complex_Numbers.Create(zero));
    res.err := zero;
    res.rco := one;
    res.res := zero;
    return res;
  end Create;

  function Create ( n,m : integer32 )
                  return DoblDobl_Complex_Solutions.Solution_List is

    use DoblDobl_Complex_Solutions;
    res,res_last : Solution_List;

  begin
    for i in 1..m loop
      Append(res,res_last,Create(n));
    end loop;
    return res;
  end Create;

  function Create ( n : integer32 )
                  return QuadDobl_Complex_Solutions.Solution is

    res : QuadDobl_Complex_Solutions.Solution(n);
    zero : constant quad_double := create(0.0);
    one : constant quad_double := create(1.0);

  begin
    res.t := QuadDobl_Complex_Numbers.Create(zero);
    res.m := 1;
    res.v := (1..n => QuadDobl_Complex_Numbers.Create(zero));
    res.err := zero;
    res.rco := one;
    res.res := zero;
    return res;
  end Create;

  function Create ( n,m : integer32 )
                  return QuadDobl_Complex_Solutions.Solution_List is

    use QuadDobl_Complex_Solutions;
    res,res_last : Solution_List;

  begin
    for i in 1..m loop
      Append(res,res_last,Create(n));
    end loop;
    return res;
  end Create;

  procedure Allocate
              ( first,last : in out Standard_Complex_Solutions.Solution_List;
                n,m : in integer32 ) is
  begin
    for i in 1..m loop
      Standard_Complex_Solutions.Append(first,last,Create(n));
    end loop;
  end Allocate;

  procedure Allocate
              ( first,last : in out DoblDobl_Complex_Solutions.Solution_List;
                n,m : in integer32 ) is
  begin
    for i in 1..m loop
      DoblDobl_Complex_Solutions.Append(first,last,Create(n));
    end loop;
  end Allocate;

  procedure Allocate
              ( first,last : in out QuadDobl_Complex_Solutions.Solution_List;
                n,m : in integer32 ) is
  begin
    for i in 1..m loop
      QuadDobl_Complex_Solutions.Append(first,last,Create(n));
    end loop;
  end Allocate;

  function Product_of_Diagonal
             ( A : Standard_Integer_Matrices.Matrix ) return integer32 is

    res : integer32 := 1;

  begin
    for i in A'range loop
      res := res*A(i,i);
    end loop;
    return res;
  end Product_of_Diagonal;

  function Product_of_Diagonal
             ( A : Standard_Integer64_Matrices.Matrix ) return integer64 is

    res : integer64 := 1;

  begin
    for i in A'range loop
      res := res*A(i,i);
    end loop;
    return res;
  end Product_of_Diagonal;

  function Volume_of_Diagonal
             ( A : Standard_Integer_Matrices.Matrix ) return natural32 is

    res : constant integer32 := Product_of_Diagonal(A);

  begin
    if res >= 0
     then return natural32(res);
     else return natural32(-res);
    end if;
  end Volume_of_Diagonal;

  function Volume_of_Diagonal
             ( A : Standard_Integer64_Matrices.Matrix ) return natural64 is

    res : integer64 := 1; -- constant integer64 := Product_of_Diagonal(A);

  begin
    for i in A'range(1) loop
      res := res*A(i,i);
    end loop;
    if res >= 0
     then return natural64(res);
     else return natural64(-res);
    end if;
  end Volume_of_Diagonal;

  function Volume_of_Cell
             ( A : Standard_Integer_Matrices.Matrix ) return natural32 is

    res : constant integer32
        := Standard_Integer_Linear_Solvers.Det(A);

  begin
    if res >= 0
     then return natural32(res);
     else return natural32(-res);
    end if;
  end Volume_of_Cell;

  function Volume_of_Cell
             ( A : Standard_Integer64_Matrices.Matrix ) return natural64 is

    res : constant integer64
        := Standard_Integer64_Linear_Solvers.Det(A);

  begin
    if res >= 0
     then return natural64(res);
     else return natural64(-res);
    end if;
  end Volume_of_Cell;

  procedure Fully_Mixed_Start_Systems
              ( q : in Standard_Complex_Laur_Systems.Laur_Sys;
                mcc : in Mixed_Subdivision ) is

    use Standard_Complex_Solutions;

    n : constant integer32 := q'last;
    cff : Standard_Complex_VecVecs.VecVec(q'range);
    exp : Standard_Integer_VecVecs.Array_of_VecVecs(q'range);
    mic : Mixed_Cell;
    tmp : Mixed_Subdivision := mcc;
    s_c : Standard_Complex_Vectors.Vector(1..2*n); -- fully mixed
    A,M,U : Standard_Integer_Matrices.Matrix(q'range,q'range);
    prod : integer32;
    b,wrk : Standard_Complex_Vectors.Vector(q'range);
    brd,logbrd,logx,e10x : Standard_Floating_Vectors.Vector(b'range);
    bsc : Standard_Complex_Vectors.Vector(b'range);
    rnk : integer32 := 0;
    len,totlen : natural32 := 0;
    res,chksum : double_float := 0.0;
    sols,Asols : Solution_List;

  begin
    Standard_Tableau_Formats.Extract_Coefficients_and_Exponents(q,cff,exp);
   -- put_line("The tableau format : "); Write_Tableau(cff,exp);
    for i in 1..Length_Of(mcc) loop
      mic := Head_Of(tmp);
     -- setting up the binomial system
      Select_Coefficients(cff,exp,mic.pts.all,s_c);
      put("sys "); put(i,1); -- put_line(" :"); Write_Tableau(s_c,mic.pts.all);
      Fully_Mixed_To_Binomial_Format(s_c,mic.pts.all,A,b);
     -- solving in place :
      U := A;
      Standard_Integer_Linear_Solvers.Upper_Triangulate(M,U);
      prod := Product_of_Diagonal(U);
      if prod < 0
       then Asols := Create(n,-prod);
       else Asols := Create(n,prod);
      end if;
      put(" : det : "); put(prod,1);
      brd := Standard_Radial_Solvers.Radii(b);
      bsc := Standard_Radial_Solvers.Scale(b,brd);
     -- Asols := Standard_Binomial_Solvers.Solve_Upper_Square(U,bsc);
      Standard_Binomial_Solvers.Solve_Upper_Square(U,bsc,Asols);
      logbrd := Standard_Radial_Solvers.Log10(brd);
      logx := Standard_Radial_Solvers.Radial_Upper_Solve(U,logbrd);
      logx := Standard_Radial_Solvers.Multiply(M,logx);
      e10x := Standard_Radial_Solvers.Exp10(logx);
      Standard_Binomial_Systems.Eval(M,Asols,wrk);
      Standard_Radial_Solvers.Multiply(Asols,e10x);
     -- solving with one call :
      Standard_Binomial_Solvers.Solve(A,b,rnk,M,U,sols);
     -- checking the solutions :
      len := Length_Of(sols);
      res := Standard_Binomial_Solvers.Sum_Residuals(A,b,Asols);
      put(" det(A) = "); put(Standard_Integer_Linear_Solvers.Det(A),1);
      put(", found "); put(len,1);
      put(" solutions, residual = "); put(res,3); new_line;
      Clear(sols); Clear(Asols);
      chksum := chksum + res;
      totlen := totlen + len;
      tmp := Tail_Of(tmp);
    end loop;
    put("number of solutions found : "); put(totlen,1); new_line;
    put("total sum over all residuals : "); put(chksum,3); new_line;
  end Fully_Mixed_Start_Systems;

  procedure Fully_Mixed_Start_Systems
              ( q : in DoblDobl_Complex_Laur_Systems.Laur_Sys;
                mcc : in Mixed_Subdivision ) is

    use DoblDobl_Complex_Solutions;

    n : constant integer32 := q'last;
    cff : DoblDobl_Complex_VecVecs.VecVec(q'range);
    exp : Standard_Integer_VecVecs.Array_of_VecVecs(q'range);
    mic : Mixed_Cell;
    tmp : Mixed_Subdivision := mcc;
    s_c : DoblDobl_Complex_Vectors.Vector(1..2*n); -- fully mixed
    A,M,U : Standard_Integer64_Matrices.Matrix(q'range,q'range);
    prod : integer64;
    b,wrk : DoblDobl_Complex_Vectors.Vector(q'range);
    brd,logbrd,logx,e10x : Double_Double_Vectors.Vector(b'range);
    bsc : DoblDobl_Complex_Vectors.Vector(b'range);
    rnk : integer32 := 0;
    len,totlen : natural32 := 0;
    res,chksum : double_double := create(0.0);
    sols,Asols : Solution_List;

  begin
    DoblDobl_Tableau_Formats.Extract_Coefficients_and_Exponents(q,cff,exp);
   -- put_line("The tableau format : "); Write_Tableau(cff,exp);
    for i in 1..Length_Of(mcc) loop
      mic := Head_Of(tmp);
     -- setting up the binomial system
      Select_Coefficients(cff,exp,mic.pts.all,s_c);
      put("sys "); put(i,1); -- put_line(" :"); Write_Tableau(s_c,mic.pts.all);
      Fully_Mixed_To_Binomial_Format(s_c,mic.pts.all,A,b);
     -- solving in place :
      U := A;
      Standard_Integer64_Linear_Solvers.Upper_Triangulate(M,U);
      prod := Product_of_Diagonal(U);
      if prod < 0
       then Asols := Create(n,integer32(-prod));
       else Asols := Create(n,integer32(prod));
      end if;
      put(" : det : "); put(prod,1);
      brd := DoblDobl_Radial_Solvers.Radii(b);
      bsc := DoblDobl_Radial_Solvers.Scale(b,brd);
     -- Asols := DoblDobl_Binomial_Solvers.Solve_Upper_Square(U,bsc);
      DoblDobl_Binomial_Solvers.Solve_Upper_Square(U,bsc,Asols);
      logbrd := DoblDobl_Radial_Solvers.Log10(brd);
      logx := DoblDobl_Radial_Solvers.Radial_Upper_Solve(U,logbrd);
      logx := DoblDobl_Radial_Solvers.Multiply(M,logx);
      e10x := DoblDobl_Radial_Solvers.Exp10(logx);
      DoblDobl_Binomial_Systems.Eval(M,Asols,wrk);
      DoblDobl_Radial_Solvers.Multiply(Asols,e10x);
     -- solving with one call :
      DoblDobl_Binomial_Solvers.Solve(A,b,rnk,M,U,sols);
     -- checking the solutions :
      len := Length_Of(sols);
      res := DoblDobl_Binomial_Solvers.Sum_Residuals(A,b,Asols);
      put(" det(A) = "); put(Standard_Integer64_Linear_Solvers.Det(A),1);
      put(", found "); put(len,1);
      put(" solutions, residual = "); put(res,3); new_line;
      Clear(sols); Clear(Asols);
      chksum := chksum + res;
      totlen := totlen + len;
      tmp := Tail_Of(tmp);
    end loop;
    put("number of solutions found : "); put(totlen,1); new_line;
    put("total sum over all residuals : "); put(chksum,3); new_line;
  end Fully_Mixed_Start_Systems;

  procedure Fully_Mixed_Start_Systems
              ( q : in QuadDobl_Complex_Laur_Systems.Laur_Sys;
                mcc : in Mixed_Subdivision ) is

    use QuadDobl_Complex_Solutions;

    n : constant integer32 := q'last;
    cff : QuadDobl_Complex_VecVecs.VecVec(q'range);
    exp : Standard_Integer_VecVecs.Array_of_VecVecs(q'range);
    mic : Mixed_Cell;
    tmp : Mixed_Subdivision := mcc;
    s_c : QuadDobl_Complex_Vectors.Vector(1..2*n); -- fully mixed
    A,M,U : Standard_Integer64_Matrices.Matrix(q'range,q'range);
    prod : integer64;
    b,wrk : QuadDobl_Complex_Vectors.Vector(q'range);
    brd,logbrd,logx,e10x : Quad_Double_Vectors.Vector(b'range);
    bsc : QuadDobl_Complex_Vectors.Vector(b'range);
    rnk : integer32 := 0;
    len,totlen : natural32 := 0;
    res,chksum : quad_double := create(0.0);
    sols,Asols : Solution_List;

  begin
    QuadDobl_Tableau_Formats.Extract_Coefficients_and_Exponents(q,cff,exp);
   -- put_line("The tableau format : "); Write_Tableau(cff,exp);
    for i in 1..Length_Of(mcc) loop
      mic := Head_Of(tmp);
     -- setting up the binomial system
      Select_Coefficients(cff,exp,mic.pts.all,s_c);
      put("sys "); put(i,1); -- put_line(" :"); Write_Tableau(s_c,mic.pts.all);
      Fully_Mixed_To_Binomial_Format(s_c,mic.pts.all,A,b);
     -- solving in place :
      U := A;
      Standard_Integer64_Linear_Solvers.Upper_Triangulate(M,U);
      prod := Product_of_Diagonal(U);
      if prod < 0
       then Asols := Create(n,integer32(-prod));
       else Asols := Create(n,integer32(prod));
      end if;
      put(" : det : "); put(prod,1);
      brd := QuadDobl_Radial_Solvers.Radii(b);
      bsc := QuadDobl_Radial_Solvers.Scale(b,brd);
     -- Asols := QuadDobl_Binomial_Solvers.Solve_Upper_Square(U,bsc);
      QuadDobl_Binomial_Solvers.Solve_Upper_Square(U,bsc,Asols);
      logbrd := QuadDobl_Radial_Solvers.Log10(brd);
      logx := QuadDobl_Radial_Solvers.Radial_Upper_Solve(U,logbrd);
      logx := QuadDobl_Radial_Solvers.Multiply(M,logx);
      e10x := QuadDobl_Radial_Solvers.Exp10(logx);
      QuadDobl_Binomial_Systems.Eval(M,Asols,wrk);
      QuadDobl_Radial_Solvers.Multiply(Asols,e10x);
     -- solving with one call :
      QuadDobl_Binomial_Solvers.Solve(A,b,rnk,M,U,sols);
     -- checking the solutions :
      len := Length_Of(sols);
      res := QuadDobl_Binomial_Solvers.Sum_Residuals(A,b,Asols);
      put(" det(A) = "); put(Standard_Integer64_Linear_Solvers.Det(A),1);
      put(", found "); put(len,1);
      put(" solutions, residual = "); put(res,3); new_line;
      Clear(sols); Clear(Asols);
      chksum := chksum + res;
      totlen := totlen + len;
      tmp := Tail_Of(tmp);
    end loop;
    put("number of solutions found : "); put(totlen,1); new_line;
    put("total sum over all residuals : "); put(chksum,3); new_line;
  end Fully_Mixed_Start_Systems;

  procedure Semi_Mixed_Start_Systems
              ( q : in Standard_Complex_Laur_Systems.Laur_Sys;
                m : in natural32;
                mix : in Standard_Integer_Vectors.Vector;
                mcc : in Mixed_Subdivision ) is

    use Standard_Complex_Solutions;

    n : constant integer32 := q'last;
    cff : Standard_Complex_VecVecs.VecVec(q'range);
    exp : Standard_Integer_VecVecs.Array_of_VecVecs(q'range);
    mic : Mixed_Cell;
    tmp : Mixed_Subdivision := mcc;
    A,T,U : Standard_Integer_Matrices.Matrix(q'range,q'range);
    C : Standard_Complex_Matrices.Matrix(q'range,q'range);
    b : Standard_Complex_Vectors.Vector(q'range);
    piv : Standard_Integer_Vectors.Vector(q'range);
    info : integer32;
    vol,mv : natural32 := 0;
    sq : Standard_Complex_Laur_Systems.Laur_Sys(q'range);
    brd,logbrd,logx,e10x : Standard_Floating_Vectors.Vector(b'range);
    wrk,bsc : Standard_Complex_Vectors.Vector(b'range);
    res,chksum : double_float := 0.0;
    sols,sols_ptr : Solution_List;
    ls : Link_to_Solution;

  begin
    Standard_Tableau_Formats.Extract_Coefficients_and_Exponents(q,cff,exp);
    for i in 1..Length_Of(mcc) loop
      put("processing cell "); put(i,1); put_line(" :");
      mic := Head_Of(tmp);
      Select_Subsystem_to_Matrix_Format(cff,exp,mix,mic.pts.all,A,C,b);
      sq := Supports_of_Polynomial_Systems.Select_Terms(q,mix,mic.pts.all);
      put_line("The subsystem : "); put_line(sq);
      put_line("The matrix A :"); put(A);
      put_line("The matrix C :"); put(C);
      put_line("The right hand side vector b :"); put_line(b);
      Standard_Complex_Linear_Solvers.lufac(C,n,piv,info);
      Standard_Complex_Linear_Solvers.lusolve(C,n,piv,b);
      U := A;
      Standard_Integer_Linear_Solvers.Upper_Triangulate(T,U);
      vol := Volume_of_Diagonal(U); mv := mv + vol;
      sols := Create(n,integer32(vol));
      brd := Standard_Radial_Solvers.Radii(b);
      bsc := Standard_Radial_Solvers.Scale(b,brd);
      Standard_Binomial_Solvers.Solve_Upper_Square(U,bsc,sols);
      logbrd := Standard_Radial_Solvers.Log10(brd);
      logx := Standard_Radial_Solvers.Radial_Upper_Solve(U,logbrd);
      logx := Standard_Radial_Solvers.Multiply(T,logx);
      e10x := Standard_Radial_Solvers.Exp10(logx);
      Standard_Binomial_Systems.Eval(T,sols,wrk);
      Standard_Radial_Solvers.Multiply(sols,e10x);
      sols_ptr := sols;
      while not Is_Null(sols_ptr) loop
        ls := Head_Of(sols_ptr);
        wrk := Standard_Complex_Laur_SysFun.Eval(sq,ls.v);
        res := Standard_Complex_Norms_Equals.Max_Norm(wrk);
        put("residual : "); put(res); new_line;
        chksum := chksum + res;
        sols_ptr := Tail_Of(sols_ptr);
      end loop;
      Clear(sols); Standard_Complex_Laur_Systems.Clear(sq);
      tmp := Tail_Of(tmp);
    end loop;
    put("the mixed volume : "); put(mv,1); new_line;
    put("sum of residuals : "); put(chksum); new_line;
  end Semi_Mixed_Start_Systems;

  procedure Semi_Mixed_Start_Systems
              ( q : in DoblDobl_Complex_Laur_Systems.Laur_Sys;
                m : in natural32;
                mix : in Standard_Integer_Vectors.Vector;
                mcc : in Mixed_Subdivision ) is

    use DoblDobl_Complex_Solutions;

    n : constant integer32 := q'last;
    cff : DoblDobl_Complex_VecVecs.VecVec(q'range);
    exp : Standard_Integer_VecVecs.Array_of_VecVecs(q'range);
    mic : Mixed_Cell;
    tmp : Mixed_Subdivision := mcc;
    A,T,U : Standard_Integer64_Matrices.Matrix(q'range,q'range);
    C : DoblDobl_Complex_Matrices.Matrix(q'range,q'range);
    b : DoblDobl_Complex_Vectors.Vector(q'range);
    piv : Standard_Integer_Vectors.Vector(q'range);
    info : integer32;
    vol,mv : natural64 := 0;
    sq : DoblDobl_Complex_Laur_Systems.Laur_Sys(q'range);
    brd,logbrd,logx,e10x : Double_Double_Vectors.Vector(b'range);
    wrk,bsc : DoblDobl_Complex_Vectors.Vector(b'range);
    res,chksum : double_double := create(0.0);
    sols,sols_ptr : Solution_List;
    ls : Link_to_Solution;

  begin
    DoblDobl_Tableau_Formats.Extract_Coefficients_and_Exponents(q,cff,exp);
    for i in 1..Length_Of(mcc) loop
      put("processing cell "); put(i,1); put_line(" :");
      mic := Head_Of(tmp);
      Select_Subsystem_to_Matrix_Format(cff,exp,mix,mic.pts.all,A,C,b);
      sq := Supports_of_Polynomial_Systems.Select_Terms(q,mix,mic.pts.all);
      put_line("The subsystem : "); put_line(sq);
      put_line("The matrix A :"); put(A);
      put_line("The matrix C :"); put(C);
      put_line("The right hand side vector b :"); put_line(b);
      DoblDobl_Complex_Linear_Solvers.lufac(C,n,piv,info);
      DoblDobl_Complex_Linear_Solvers.lusolve(C,n,piv,b);
      U := A;
      Standard_Integer64_Linear_Solvers.Upper_Triangulate(T,U);
      vol := Volume_of_Diagonal(U); mv := mv + vol;
      sols := Create(n,integer32(vol));
      brd := DoblDobl_Radial_Solvers.Radii(b);
      bsc := DoblDobl_Radial_Solvers.Scale(b,brd);
      DoblDobl_Binomial_Solvers.Solve_Upper_Square(U,bsc,sols);
      logbrd := DoblDobl_Radial_Solvers.Log10(brd);
      logx := DoblDobl_Radial_Solvers.Radial_Upper_Solve(U,logbrd);
      logx := DoblDobl_Radial_Solvers.Multiply(T,logx);
      e10x := DoblDobl_Radial_Solvers.Exp10(logx);
      DoblDobl_Binomial_Systems.Eval(T,sols,wrk);
      DoblDobl_Radial_Solvers.Multiply(sols,e10x);
      sols_ptr := sols;
      while not Is_Null(sols_ptr) loop
        ls := Head_Of(sols_ptr);
        wrk := DoblDobl_Complex_Laur_SysFun.Eval(sq,ls.v);
        res := DoblDobl_Complex_Vector_Norms.Max_Norm(wrk);
        put("residual : "); put(res); new_line;
        chksum := chksum + res;
        sols_ptr := Tail_Of(sols_ptr);
      end loop;
      Clear(sols); DoblDobl_Complex_Laur_Systems.Clear(sq);
      tmp := Tail_Of(tmp);
    end loop;
    put("the mixed volume : "); put(mv,1); new_line;
    put("sum of residuals : "); put(chksum); new_line;
  end Semi_Mixed_Start_Systems;

  procedure Semi_Mixed_Start_Systems
              ( q : in QuadDobl_Complex_Laur_Systems.Laur_Sys;
                m : in natural32;
                mix : in Standard_Integer_Vectors.Vector;
                mcc : in Mixed_Subdivision ) is

    use QuadDobl_Complex_Solutions;

    n : constant integer32 := q'last;
    cff : QuadDobl_Complex_VecVecs.VecVec(q'range);
    exp : Standard_Integer_VecVecs.Array_of_VecVecs(q'range);
    mic : Mixed_Cell;
    tmp : Mixed_Subdivision := mcc;
    A,T,U : Standard_Integer64_Matrices.Matrix(q'range,q'range);
    C : QuadDobl_Complex_Matrices.Matrix(q'range,q'range);
    b : QuadDobl_Complex_Vectors.Vector(q'range);
    piv : Standard_Integer_Vectors.Vector(q'range);
    info : integer32;
    vol,mv : natural64 := 0;
    sq : QuadDobl_Complex_Laur_Systems.Laur_Sys(q'range);
    brd,logbrd,logx,e10x : Quad_Double_Vectors.Vector(b'range);
    wrk,bsc : QuadDobl_Complex_Vectors.Vector(b'range);
    res,chksum : quad_double := create(0.0);
    sols,sols_ptr : Solution_List;
    ls : Link_to_Solution;

  begin
    QuadDobl_Tableau_Formats.Extract_Coefficients_and_Exponents(q,cff,exp);
    for i in 1..Length_Of(mcc) loop
      put("processing cell "); put(i,1); put_line(" :");
      mic := Head_Of(tmp);
      Select_Subsystem_to_Matrix_Format(cff,exp,mix,mic.pts.all,A,C,b);
      sq := Supports_of_Polynomial_Systems.Select_Terms(q,mix,mic.pts.all);
      put_line("The subsystem : "); put_line(sq);
      put_line("The matrix A :"); put(A);
      put_line("The matrix C :"); put(C);
      put_line("The right hand side vector b :"); put_line(b);
      QuadDobl_Complex_Linear_Solvers.lufac(C,n,piv,info);
      QuadDobl_Complex_Linear_Solvers.lusolve(C,n,piv,b);
      U := A;
      Standard_Integer64_Linear_Solvers.Upper_Triangulate(T,U);
      vol := Volume_of_Diagonal(U); mv := mv + vol;
      sols := Create(n,integer32(vol));
      brd := QuadDobl_Radial_Solvers.Radii(b);
      bsc := QuadDobl_Radial_Solvers.Scale(b,brd);
      QuadDobl_Binomial_Solvers.Solve_Upper_Square(U,bsc,sols);
      logbrd := QuadDobl_Radial_Solvers.Log10(brd);
      logx := QuadDobl_Radial_Solvers.Radial_Upper_Solve(U,logbrd);
      logx := QuadDobl_Radial_Solvers.Multiply(T,logx);
      e10x := QuadDobl_Radial_Solvers.Exp10(logx);
      QuadDobl_Binomial_Systems.Eval(T,sols,wrk);
      QuadDobl_Radial_Solvers.Multiply(sols,e10x);
      sols_ptr := sols;
      while not Is_Null(sols_ptr) loop
        ls := Head_Of(sols_ptr);
        wrk := QuadDobl_Complex_Laur_SysFun.Eval(sq,ls.v);
        res := QuadDobl_Complex_Vector_Norms.Max_Norm(wrk);
        put("residual : "); put(res); new_line;
        chksum := chksum + res;
        sols_ptr := Tail_Of(sols_ptr);
      end loop;
      Clear(sols); QuadDobl_Complex_Laur_Systems.Clear(sq);
      tmp := Tail_Of(tmp);
    end loop;
    put("the mixed volume : "); put(mv,1); new_line;
    put("sum of residuals : "); put(chksum); new_line;
  end Semi_Mixed_Start_Systems;

  procedure Check_Solutions
              ( cff : in Standard_Complex_VecVecs.VecVec;
                exp : in Standard_Integer_VecVecs.Array_of_VecVecs;
                mcc : in Mixed_Subdivision;
                sols : in Standard_Complex_Solutions.Array_of_Solution_Lists;
                res : out Standard_Floating_Vectors.Vector ) is

    use Standard_Complex_Solutions;

    n : constant integer32 := cff'last;
    nt : constant integer32 := sols'last;
    tmp : Mixed_Subdivision := mcc;
    mic : Mixed_Cell;
    s_c : Standard_Complex_Vectors.Vector(1..2*n); -- fully mixed
    A : Standard_Integer_Matrices.Matrix(1..n,1..n);
    b : Standard_Complex_Vectors.Vector(1..n);
    pdetA : natural32;
    ptrs : Array_of_Solution_Lists(sols'range);
    Asols,Asols_last : Solution_List;
    ind : integer32;
    r : double_float;

  begin
    for i in sols'range loop
      ptrs(i) := sols(i);
      res(i) := 0.0;
    end loop;
    for k in 1..Length_Of(mcc) loop
      mic := Head_Of(tmp);
      Select_Coefficients(cff,exp,mic.pts.all,s_c);
      Fully_Mixed_To_Binomial_Format(s_c,mic.pts.all,A,b);
      pdetA := Volume_of_Cell(A);
      ind := (integer32(k) mod nt) + 1;     -- task index tells where computed
      for i in 1..pdetA loop   -- select next pdetA solutions
        Append(Asols,Asols_last,Head_Of(ptrs(ind)).all);
        ptrs(ind) := Tail_Of(ptrs(ind));
      end loop;
      r := Standard_Binomial_Solvers.Sum_Residuals(A,b,Asols);
      res(ind) := res(ind) + r;
      Clear(Asols);
      tmp := Tail_Of(tmp);
    end loop;
  end Check_Solutions;

  procedure Check_Solutions
              ( cff : in DoblDobl_Complex_VecVecs.VecVec;
                exp : in Standard_Integer_VecVecs.Array_of_VecVecs;
                mcc : in Mixed_Subdivision;
                sols : in DoblDobl_Complex_Solutions.Array_of_Solution_Lists;
                res : out Double_Double_Vectors.Vector ) is

    use DoblDobl_Complex_Solutions;

    n : constant integer32 := cff'last;
    nt : constant integer32 := sols'last;
    tmp : Mixed_Subdivision := mcc;
    mic : Mixed_Cell;
    s_c : DoblDobl_Complex_Vectors.Vector(1..2*n); -- fully mixed
    A : Standard_Integer64_Matrices.Matrix(1..n,1..n);
    b : DoblDobl_Complex_Vectors.Vector(1..n);
    pdetA : natural64;
    ptrs : Array_of_Solution_Lists(sols'range);
    Asols,Asols_last : Solution_List;
    ind : integer32;
    r : double_double;

  begin
    for i in sols'range loop
      ptrs(i) := sols(i);
      res(i) := create(0.0);
    end loop;
    for k in 1..Length_Of(mcc) loop
      mic := Head_Of(tmp);
      Select_Coefficients(cff,exp,mic.pts.all,s_c);
      Fully_Mixed_To_Binomial_Format(s_c,mic.pts.all,A,b);
      pdetA := Volume_of_Cell(A);
      ind := (integer32(k) mod nt) + 1;     -- task index tells where computed
      for i in 1..pdetA loop   -- select next pdetA solutions
        Append(Asols,Asols_last,Head_Of(ptrs(ind)).all);
        ptrs(ind) := Tail_Of(ptrs(ind));
      end loop;
      r := DoblDobl_Binomial_Solvers.Sum_Residuals(A,b,Asols);
      res(ind) := res(ind) + r;
      Clear(Asols);
      tmp := Tail_Of(tmp);
    end loop;
  end Check_Solutions;

  procedure Check_Solutions
              ( cff : in QuadDobl_Complex_VecVecs.VecVec;
                exp : in Standard_Integer_VecVecs.Array_of_VecVecs;
                mcc : in Mixed_Subdivision;
                sols : in QuadDobl_Complex_Solutions.Array_of_Solution_Lists;
                res : out Quad_Double_Vectors.Vector ) is

    use QuadDobl_Complex_Solutions;

    n : constant integer32 := cff'last;
    nt : constant integer32 := sols'last;
    tmp : Mixed_Subdivision := mcc;
    mic : Mixed_Cell;
    s_c : QuadDobl_Complex_Vectors.Vector(1..2*n); -- fully mixed
    A : Standard_Integer64_Matrices.Matrix(1..n,1..n);
    b : QuadDobl_Complex_Vectors.Vector(1..n);
    pdetA : natural64;
    ptrs : Array_of_Solution_Lists(sols'range);
    Asols,Asols_last : Solution_List;
    ind : integer32;
    r : quad_double;

  begin
    for i in sols'range loop
      ptrs(i) := sols(i);
      res(i) := create(0.0);
    end loop;
    for k in 1..Length_Of(mcc) loop
      mic := Head_Of(tmp);
      Select_Coefficients(cff,exp,mic.pts.all,s_c);
      Fully_Mixed_To_Binomial_Format(s_c,mic.pts.all,A,b);
      pdetA := Volume_of_Cell(A);
      ind := (integer32(k) mod nt) + 1;     -- task index tells where computed
      for i in 1..pdetA loop   -- select next pdetA solutions
        Append(Asols,Asols_last,Head_Of(ptrs(ind)).all);
        ptrs(ind) := Tail_Of(ptrs(ind));
      end loop;
      r := QuadDobl_Binomial_Solvers.Sum_Residuals(A,b,Asols);
      res(ind) := res(ind) + r;
      Clear(Asols);
      tmp := Tail_Of(tmp);
    end loop;
  end Check_Solutions;

  procedure Check_Solutions
              ( q : in Standard_Complex_Laur_Systems.Laur_Sys;
                mcc : in Mixed_Subdivision;
                sols : in Standard_Complex_Solutions.Array_of_Solution_Lists;
                res : out Standard_Floating_Vectors.Vector ) is

    cff : Standard_Complex_VecVecs.VecVec(q'range);
    exp : Standard_Integer_VecVecs.Array_of_VecVecs(q'range);

  begin
    Standard_Tableau_Formats.Extract_Coefficients_and_Exponents(q,cff,exp);
    Check_Solutions(cff,exp,mcc,sols,res);
  end Check_Solutions;

  procedure Check_Solutions
              ( q : in DoblDobl_Complex_Laur_Systems.Laur_Sys;
                mcc : in Mixed_Subdivision;
                sols : in DoblDobl_Complex_Solutions.Array_of_Solution_Lists;
                res : out Double_Double_Vectors.Vector ) is

    cff : DoblDobl_Complex_VecVecs.VecVec(q'range);
    exp : Standard_Integer_VecVecs.Array_of_VecVecs(q'range);

  begin
    DoblDobl_Tableau_Formats.Extract_Coefficients_and_Exponents(q,cff,exp);
    Check_Solutions(cff,exp,mcc,sols,res);
  end Check_Solutions;

  procedure Check_Solutions
              ( q : in QuadDobl_Complex_Laur_Systems.Laur_Sys;
                mcc : in Mixed_Subdivision;
                sols : in QuadDobl_Complex_Solutions.Array_of_Solution_Lists;
                res : out Quad_Double_Vectors.Vector ) is

    cff : QuadDobl_Complex_VecVecs.VecVec(q'range);
    exp : Standard_Integer_VecVecs.Array_of_VecVecs(q'range);

  begin
    QuadDobl_Tableau_Formats.Extract_Coefficients_and_Exponents(q,cff,exp);
    Check_Solutions(cff,exp,mcc,sols,res);
  end Check_Solutions;

  procedure Check_Solutions
              ( q : in Standard_Complex_Laur_Systems.Laur_Sys;
                mix : in Standard_Integer_Vectors.Vector;
                mcc : in Mixed_Subdivision;
                sols : in Standard_Complex_Solutions.Array_of_Solution_Lists;
                res : out Standard_Floating_Vectors.Vector ) is

    use Standard_Complex_Solutions;

    nt : constant integer32 := sols'last;
    ptrs : Array_of_Solution_Lists(sols'range);
    ls : Link_to_Solution;
    cff : Standard_Complex_VecVecs.VecVec(q'range);
    exp : Standard_Integer_VecVecs.Array_of_VecVecs(q'range);
    tmp : Mixed_Subdivision := mcc;
    mic : Mixed_Cell;
    A : Standard_Integer_Matrices.Matrix(q'range,q'range);
    C : Standard_Complex_Matrices.Matrix(q'range,q'range);
    b : Standard_Complex_Vectors.Vector(q'range);
    vol : natural32;
    ind : integer32;
    r : double_float;
    sq : Standard_Complex_Laur_Systems.Laur_Sys(q'range);

  begin
    Standard_Tableau_Formats.Extract_Coefficients_and_Exponents(q,cff,exp);
    for i in sols'range loop
      ptrs(i) := sols(i);
      res(i) := 0.0;
    end loop;
    for k in 1..Length_Of(mcc) loop
      mic := Head_Of(tmp);
      ind := (integer32(k) mod nt) + 1;
      Select_Subsystem_to_Matrix_Format(cff,exp,mix,mic.pts.all,A,C,b);
      vol := Volume_of_Cell(A);
      sq := Supports_of_Polynomial_Systems.Select_Terms(q,mix,mic.pts.all);
      for i in 1..vol loop
        ls := Head_Of(ptrs(ind));
        b := Standard_Complex_Laur_SysFun.Eval(sq,ls.v);
        r := Standard_Complex_Norms_Equals.Max_Norm(b);
        res(ind) := res(ind) + r;
        ptrs(ind) := Tail_Of(ptrs(ind));
      end loop;
      tmp := Tail_Of(tmp);
    end loop;
  end Check_Solutions;

  procedure Check_Solutions
              ( q : in DoblDobl_Complex_Laur_Systems.Laur_Sys;
                mix : in Standard_Integer_Vectors.Vector;
                mcc : in Mixed_Subdivision;
                sols : in DoblDobl_Complex_Solutions.Array_of_Solution_Lists;
                res : out Double_Double_Vectors.Vector ) is

    use DoblDobl_Complex_Solutions;

    nt : constant integer32 := sols'last;
    ptrs : Array_of_Solution_Lists(sols'range);
    ls : Link_to_Solution;
    cff : DoblDobl_Complex_VecVecs.VecVec(q'range);
    exp : Standard_Integer_VecVecs.Array_of_VecVecs(q'range);
    tmp : Mixed_Subdivision := mcc;
    mic : Mixed_Cell;
    A : Standard_Integer64_Matrices.Matrix(q'range,q'range);
    C : DoblDobl_Complex_Matrices.Matrix(q'range,q'range);
    b : DoblDobl_Complex_Vectors.Vector(q'range);
    vol : natural64;
    ind : integer32;
    r : double_double;
    sq : DoblDobl_Complex_Laur_Systems.Laur_Sys(q'range);

  begin
    DoblDobl_Tableau_Formats.Extract_Coefficients_and_Exponents(q,cff,exp);
    for i in sols'range loop
      ptrs(i) := sols(i);
      res(i) := create(0.0);
    end loop;
    for k in 1..Length_Of(mcc) loop
      mic := Head_Of(tmp);
      ind := (integer32(k) mod nt) + 1;
      Select_Subsystem_to_Matrix_Format(cff,exp,mix,mic.pts.all,A,C,b);
      vol := Volume_of_Cell(A);
      sq := Supports_of_Polynomial_Systems.Select_Terms(q,mix,mic.pts.all);
      for i in 1..vol loop
        ls := Head_Of(ptrs(ind));
        b := DoblDobl_Complex_Laur_SysFun.Eval(sq,ls.v);
        r := DoblDobl_Complex_Vector_Norms.Max_Norm(b);
        res(ind) := res(ind) + r;
        ptrs(ind) := Tail_Of(ptrs(ind));
      end loop;
      tmp := Tail_Of(tmp);
    end loop;
  end Check_Solutions;

  procedure Check_Solutions
              ( q : in QuadDobl_Complex_Laur_Systems.Laur_Sys;
                mix : in Standard_Integer_Vectors.Vector;
                mcc : in Mixed_Subdivision;
                sols : in QuadDobl_Complex_Solutions.Array_of_Solution_Lists;
                res : out Quad_Double_Vectors.Vector ) is

    use QuadDobl_Complex_Solutions;

    nt : constant integer32 := sols'last;
    ptrs : Array_of_Solution_Lists(sols'range);
    ls : Link_to_Solution;
    cff : QuadDobl_Complex_VecVecs.VecVec(q'range);
    exp : Standard_Integer_VecVecs.Array_of_VecVecs(q'range);
    tmp : Mixed_Subdivision := mcc;
    mic : Mixed_Cell;
    A : Standard_Integer64_Matrices.Matrix(q'range,q'range);
    C : QuadDobl_Complex_Matrices.Matrix(q'range,q'range);
    b : QuadDobl_Complex_Vectors.Vector(q'range);
    vol : natural64;
    ind : integer32;
    r : quad_double;
    sq : QuadDobl_Complex_Laur_Systems.Laur_Sys(q'range);

  begin
    QuadDobl_Tableau_Formats.Extract_Coefficients_and_Exponents(q,cff,exp);
    for i in sols'range loop
      ptrs(i) := sols(i);
      res(i) := create(0.0);
    end loop;
    for k in 1..Length_Of(mcc) loop
      mic := Head_Of(tmp);
      ind := (integer32(k) mod nt) + 1;
      Select_Subsystem_to_Matrix_Format(cff,exp,mix,mic.pts.all,A,C,b);
      vol := Volume_of_Cell(A);
      sq := Supports_of_Polynomial_Systems.Select_Terms(q,mix,mic.pts.all);
      for i in 1..vol loop
        ls := Head_Of(ptrs(ind));
        b := QuadDobl_Complex_Laur_SysFun.Eval(sq,ls.v);
        r := QuadDobl_Complex_Vector_Norms.Max_Norm(b);
        res(ind) := res(ind) + r;
        ptrs(ind) := Tail_Of(ptrs(ind));
      end loop;
      tmp := Tail_Of(tmp);
    end loop;
  end Check_Solutions;

-- (3) ALLOCATING WORK SPACE FOR EXPONENTS AND COEFFICIENTS :

  procedure Allocate_Workspace_for_Exponents
              ( epv : in Exponent_Vectors.Exponent_Vectors_Array;
                dpw : in out Standard_Floating_VecVecs.Array_of_VecVecs ) is
  begin
    for i in dpw'range loop
      dpw(i) := new Standard_Floating_VecVecs.VecVec(epv'range);
      for k in dpw(i)'range loop
        dpw(i)(k) := new Standard_Floating_Vectors.Vector(epv(k)'range);
      end loop;
    end loop;
  end Allocate_Workspace_for_Exponents;

  procedure Allocate_Workspace_for_Coefficients
              ( cff : in Standard_Complex_VecVecs.VecVec;
                cft : in out Standard_Complex_VecVecs.Array_of_VecVecs ) is
  begin
    for i in cft'range loop
      cft(i) := new Standard_Complex_VecVecs.VecVec(cff'range);
      for k in cft(i)'range loop
        cft(i)(k) := new Standard_Complex_Vectors.Vector(cff(k)'range);
      end loop;
    end loop;
  end Allocate_Workspace_for_Coefficients;

  procedure Allocate_Workspace_for_Coefficients
              ( cff : in DoblDobl_Complex_VecVecs.VecVec;
                cft : in out DoblDobl_Complex_VecVecs.Array_of_VecVecs ) is
  begin
    for i in cft'range loop
      cft(i) := new DoblDobl_Complex_VecVecs.VecVec(cff'range);
      for k in cft(i)'range loop
        cft(i)(k) := new DoblDobl_Complex_Vectors.Vector(cff(k)'range);
      end loop;
    end loop;
  end Allocate_Workspace_for_Coefficients;

  procedure Allocate_Workspace_for_Coefficients
              ( cff : in QuadDobl_Complex_VecVecs.VecVec;
                cft : in out QuadDobl_Complex_VecVecs.Array_of_VecVecs ) is
  begin
    for i in cft'range loop
      cft(i) := new QuadDobl_Complex_VecVecs.VecVec(cff'range);
      for k in cft(i)'range loop
        cft(i)(k) := new QuadDobl_Complex_Vectors.Vector(cff(k)'range);
      end loop;
    end loop;
  end Allocate_Workspace_for_Coefficients;

end Polyhedral_Start_Systems;
