with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;       use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;         use Standard_Floating_Numbers;
with Double_Double_Numbers;             use Double_Double_Numbers;
with Quad_Double_Numbers;               use Quad_Double_Numbers;
with Standard_Complex_Numbers;
with DoblDobl_Complex_Numbers;
with QuadDobl_Complex_Numbers;
with Standard_Natural_Vectors_io;       use Standard_Natural_Vectors_io;
with Standard_Natural_Matrices_io;      use Standard_Natural_Matrices_io;
with Standard_Complex_Vectors_io;       use Standard_Complex_Vectors_io;
with Standard_Complex_Matrices_io;      use Standard_Complex_Matrices_io;
with DoblDobl_Complex_Vectors_io;       use DoblDobl_Complex_Vectors_io;
with DoblDobl_Complex_Matrices_io;      use DoblDobl_Complex_Matrices_io;
with QuadDobl_Complex_Vectors_io;       use QuadDobl_Complex_Vectors_io;
with QuadDobl_Complex_Matrices_io;      use QuadDobl_Complex_Matrices_io;
with Standard_Complex_Polynomials;
with Standard_Complex_Poly_Functions;   use Standard_Complex_Poly_Functions;
with DoblDobl_Complex_Polynomials;
with DoblDobl_Complex_Poly_Functions;   use DoblDobl_Complex_Poly_Functions;
with QuadDobl_Complex_Polynomials;
with QuadDobl_Complex_Poly_Functions;   use QuadDobl_Complex_Poly_Functions;
with Standard_Complex_Poly_Matrices_io; use Standard_Complex_Poly_Matrices_io;
with DoblDobl_Complex_Poly_Matrices_io; use DoblDobl_Complex_Poly_Matrices_io;
with QuadDobl_Complex_Poly_Matrices_io; use QuadDobl_Complex_Poly_Matrices_io;
with Standard_Complex_Solutions_io;     use Standard_Complex_Solutions_io;
with DoblDobl_Complex_Solutions_io;     use DoblDobl_Complex_Solutions_io;
with QuadDobl_Complex_Solutions_io;     use QuadDobl_Complex_Solutions_io;
with Standard_Matrix_Inversion;
with DoblDobl_Matrix_Inversion;
with QuadDobl_Matrix_Inversion;
with Checker_Moves;                     use Checker_Moves;
with Checker_Localization_Patterns;     use Checker_Localization_Patterns;

package body Checker_Homotopies is

-- PART I : classification of Littlewood-Richardson homotopies

  procedure Define_Specializing_Homotopy
              ( n : in integer32;
                p,row,col : in Standard_Natural_Vectors.Vector;
                homtp,ctr : out integer32 ) is
  begin
    Define_Specializing_Homotopy(standard_output,n,p,row,col,homtp,ctr);
  end Define_Specializing_Homotopy;

  procedure Define_Specializing_Homotopy
              ( file : in file_type; n : in integer32;
                p,row,col : in Standard_Natural_Vectors.Vector;
                homtp,ctr : out integer32 ) is

    d : constant integer32 := Descending_Checker(p);
    r,cr,cd : integer32;
    np : Standard_Natural_Vectors.Vector(p'range) := p;
   -- b : Board(1..n,1..n);
    stay,swap : boolean;

  begin
    if d = 0 then
      put_line(file,"there is no homotopy to define...");
    else 
      put(file,"The descending checker is at (");
      put(file,p(d),1); put(file,","); put(file,p'last-d+1,1);
      put_line(file,")");
      r := Rising_Checker(p,d);
      put(file,"The rising checker is at (");
      put(file,p(r),1); put(file,","); put(file,p'last-r+1,1);
      put_line(file,")");
      Move(np,d,r);
     -- b := Configuration(np);
     -- put_line(file,"The board after moving black checkers : ");
     -- Write(file,b);
      cr := Critical_Row(integer32(p(d)),n-d+1,row,col);
      put(file,"The critical row r = "); put(file,p(d),1); new_line(file);
      ctr := cr;
      if cr = 0 then
        put_line(file,"There are no white checkers in critical row.");
        put_line(file,"Trivial stay case, no homotopy,"
                    & " only change of variables."); homtp := 0;
        swap := false; stay := true;
      else
        put_line(file,"There are white checkers in critical row.");
        cd := Top_White_Checker(integer32(p(r)),n-r+1,n,row,col);
        put(file,"The critical diagonal : "); put(file,cd,1); new_line(file);
        if cd = 0 then
          swap := false; stay := true;
        else
          case Central_Choice(file,p,d,r,row,col,cr,cd,true) is
            when 0 => swap := false; stay := true;
            when 1 => swap := true;  stay := false;
            when others => swap := true; stay := true;
          end case;
        end if;
        if stay
         then put_line(file,"homotopy needed for stay case"); homtp := 1;
        end if;
        if swap
         then put_line(file,"homotopy needed for swap case"); homtp := 2;
        end if;
      end if;
    end if;
  end Define_Specializing_Homotopy;

  procedure Define_Generalizing_Homotopy 
              ( file : in file_type; n : in integer32;
                p,row,col : in Standard_Natural_Vectors.Vector;
                stay_child : in boolean; homtp,ctr : out integer32 ) is

    d : constant integer32 := Descending_Checker(p);
    r,cr,cd,little_r,big_R : integer32;
    np : Standard_Natural_Vectors.Vector(p'range) := p;
   -- b : Board(1..n,1..n);
   -- stay : boolean;
    swap : boolean;

  begin
    if d = 0 then
      put_line(file,"there is no homotopy to define...");
    else 
      put(file,"The descending checker is at (");
      put(file,p(d),1); put(file,","); put(file,d,1); put_line(file,")");
      r := Rising_Checker(p,d);
      put(file,"The rising checker is at (");
      put(file,p(r),1); put(file,","); put(file,p'last-r+1,1);
      put_line(file,")");
      Move(np,d,r);
     -- b := Configuration(np);
     -- put_line(file,"The board after moving black checkers : ");
     -- Write(file,b);
      cr := Critical_Row(integer32(p(d)),n-d+1,row,col);
      little_r := integer32(p(d)); ctr := little_r;
      put(file,"The critical row r = "); put(file,little_r,1); new_line(file);
      if cr = 0 then
        put_line(file,"There are no white checkers in critical row.");
        put_line(file,"Trivial stay case, no homotopy,"
                    & " only change of variables."); homtp := 0;
        swap := false; -- stay := true;
      else
        put_line(file,"There are white checkers in critical row.");
        cd := Top_White_Checker(integer32(p(r)),n-r+1,n,row,col);
        put(file,"The critical diagonal : "); put(file,cd,1); new_line(file);
        if cd = 0 then
          swap := false; -- stay := true;
        else
          case Central_Choice(file,p,d,r,row,col,cr,cd,true) is
            when 0 => swap := false; -- stay := true;
            when 1 => swap := true; -- stay := false;
            when others => swap := true; -- stay := true;
          end case;
          if swap then
            put(file,"row of white checker in critical diagonal R = ");
            big_R := integer32(row(cd));
            put(file,big_R,1); new_line(file);
          end if;
        end if;
        if stay_child then
          put_line(file,"homotopy needed for stay case"); homtp := 1;
        else
          put(file,"homotopy needed for swap case "); homtp := 2;
          if big_R = little_r + 1 then
            put(file,"R(="); put(file,big_R,1); put(file,") = ");
            put(file,"r(="); put(file,little_r,1); put_line(file,") + 1");
          else
            put(file,"R(="); put(file,big_R,1); put(file,") > ");
            put(file,"r(="); put(file,little_r,1); put_line(file,") + 1");
          end if;
        end if;
      end if;
    end if;
  end Define_Generalizing_Homotopy;

  procedure Define_Generalizing_Homotopy 
              ( n : in integer32;
                p,row,col : in Standard_Natural_Vectors.Vector;
                stay_child : in boolean; homtp,ctr : out integer32 ) is

    d : constant integer32 := Descending_Checker(p);
    r,cr,cd,little_r,big_R : integer32;
    np : Standard_Natural_Vectors.Vector(p'range) := p;
    swap : boolean;

  begin
    if d /= 0 then
      r := Rising_Checker(p,d);
      Move(np,d,r);
      cr := Critical_Row(integer32(p(d)),n-d+1,row,col);
      little_r := integer32(p(d)); ctr := little_r;
      if cr = 0 then
        homtp := 0;
        swap := false;
      else
        cd := Top_White_Checker(integer32(p(r)),n-r+1,n,row,col);
        if cd = 0 then
          swap := false;
        else
          case Central_Choice(p,d,r,row,col,cr,cd,true) is
            when 0 => swap := false; 
            when 1 => swap := true;
            when others => swap := true;
          end case;
          if swap
           then big_R := integer32(row(cd));
          end if;
        end if;
        if stay_child
         then homtp := 1;
         else homtp := 2;
        end if;
      end if;
    end if;
  end Define_Generalizing_Homotopy;

-- PART II : coordinate transformations on input flags and solution plane

  procedure Inverse_Row_Transformation
              ( r : in integer32;
                x : in out Standard_Complex_Matrices.Matrix ) is

    use Standard_Complex_Numbers;
    tmp : Complex_Number;

  begin
    for k in x'range(2) loop
      tmp := x(r,k);
      x(r,k) := -x(r+1,k);
      x(r+1,k) := tmp + x(r+1,k);
    end loop;
  end Inverse_Row_Transformation;

  procedure Inverse_Row_Transformation
              ( r : in integer32;
                x : in out DoblDobl_Complex_Matrices.Matrix ) is

    use DoblDobl_Complex_Numbers;
    tmp : Complex_Number;

  begin
    for k in x'range(2) loop
      tmp := x(r,k);
      x(r,k) := -x(r+1,k);
      x(r+1,k) := tmp + x(r+1,k);
    end loop;
  end Inverse_Row_Transformation;

  procedure Inverse_Row_Transformation
              ( r : in integer32;
                x : in out QuadDobl_Complex_Matrices.Matrix ) is

    use QuadDobl_Complex_Numbers;
    tmp : Complex_Number;

  begin
    for k in x'range(2) loop
      tmp := x(r,k);
      x(r,k) := -x(r+1,k);
      x(r+1,k) := tmp + x(r+1,k);
    end loop;
  end Inverse_Row_Transformation;

  procedure Inverse_Row_Transformation
              ( mf : in Standard_Complex_Matrices.Matrix;
                x : in out Standard_Complex_Matrices.Matrix ) is

    invm : constant Standard_Complex_Matrices.Matrix(mf'range(1),mf'range(2))
         := Standard_Matrix_Inversion.Inverse(mf);

    use Standard_Complex_Matrices;

  begin
    x := invm*x;
  end Inverse_Row_Transformation;

  procedure Inverse_Row_Transformation
              ( mf : in DoblDobl_Complex_Matrices.Matrix;
                x : in out DoblDobl_Complex_Matrices.Matrix ) is

    invm : constant DoblDobl_Complex_Matrices.Matrix(mf'range(1),mf'range(2))
         := DoblDobl_Matrix_Inversion.Inverse(mf);

    use DoblDobl_Complex_Matrices;

  begin
    x := invm*x;
  end Inverse_Row_Transformation;

  procedure Inverse_Row_Transformation
              ( mf : in QuadDobl_Complex_Matrices.Matrix;
                x : in out QuadDobl_Complex_Matrices.Matrix ) is

    invm : constant QuadDobl_Complex_Matrices.Matrix(mf'range(1),mf'range(2))
         := QuadDobl_Matrix_Inversion.Inverse(mf);

    use QuadDobl_Complex_Matrices;

  begin
    x := invm*x;
  end Inverse_Row_Transformation;

  procedure Normalize_to_Fit
              ( pattern : in Standard_Natural_Matrices.Matrix;
                x : in out Standard_Complex_Matrices.Matrix ) is

    use Standard_Complex_Numbers;
    pivot : integer32;

  begin
    for k in pattern'range(2) loop
      pivot := 0;
      for i in pattern'range(1) loop
        if pattern(i,k) = 1
         then pivot := i;
        end if;
        exit when (pivot > 0);
      end loop;
      for i in x'first(1)..pivot-1 loop
        Div(x(i,k),x(pivot,k));
      end loop;
      for i in pivot+1..x'last(1) loop
        Div(x(i,k),x(pivot,k));
      end loop;
      x(pivot,k) := Create(1.0);
    end loop;
  end Normalize_to_Fit;

  procedure Normalize_to_Fit
              ( pattern : in Standard_Natural_Matrices.Matrix;
                x : in out DoblDobl_Complex_Matrices.Matrix ) is

    use DoblDobl_Complex_Numbers;
    pivot : integer32;

  begin
    for k in pattern'range(2) loop
      pivot := 0;
      for i in pattern'range(1) loop
        if pattern(i,k) = 1
         then pivot := i;
        end if;
        exit when (pivot > 0);
      end loop;
      for i in x'first(1)..pivot-1 loop
        Div(x(i,k),x(pivot,k));
      end loop;
      for i in pivot+1..x'last(1) loop
        Div(x(i,k),x(pivot,k));
      end loop;
      x(pivot,k) := Create(integer(1));
    end loop;
  end Normalize_to_Fit;

  procedure Normalize_to_Fit
              ( pattern : in Standard_Natural_Matrices.Matrix;
                x : in out QuadDobl_Complex_Matrices.Matrix ) is

    use QuadDobl_Complex_Numbers;
    pivot : integer32;

  begin
    for k in pattern'range(2) loop
      pivot := 0;
      for i in pattern'range(1) loop
        if pattern(i,k) = 1
         then pivot := i;
        end if;
        exit when (pivot > 0);
      end loop;
      for i in x'first(1)..pivot-1 loop
        Div(x(i,k),x(pivot,k));
      end loop;
      for i in pivot+1..x'last(1) loop
        Div(x(i,k),x(pivot,k));
      end loop;
      x(pivot,k) := Create(integer(1));
    end loop;
  end Normalize_to_Fit;

  procedure Reduce_to_Fit
              ( pattern : in Standard_Natural_Matrices.Matrix;
                x : in out Standard_Complex_Matrices.Matrix ) is

    use Standard_Complex_Numbers;
    tol : constant double_float := 1.0e-8;
    avx : double_float;
    pivot : integer32;

  begin
    for k in pattern'range(2) loop
      for i in pattern'range(1) loop
        if pattern(i,k) = 0 then
          avx := AbsVal(x(i,k));
          if avx > tol then
            pivot := 0;
            for j in 1..k-1 loop
              if pattern(i,j) = 1
               then pivot := j;
              end if;
              exit when (pivot > 0);
            end loop;
            if pivot > 0 then
              for ii in x'first(1)..i-1 loop
                x(ii,k) := x(ii,k) - x(i,k)*x(ii,pivot);
              end loop;
              for ii in i+1..x'last(1) loop
                x(ii,k) := x(ii,k) - x(i,k)*x(ii,pivot);
              end loop;
              x(i,k) := Create(0.0);
            end if;
          end if;
        end if;
      end loop;
    end loop;
  end Reduce_to_Fit;

  procedure Reduce_to_Fit
              ( pattern : in Standard_Natural_Matrices.Matrix;
                x : in out DoblDobl_Complex_Matrices.Matrix ) is

    use DoblDobl_Complex_Numbers;
    tol : constant double_float := 1.0e-8;
    avx : double_double;
    pivot : integer32;

  begin
    for k in pattern'range(2) loop
      for i in pattern'range(1) loop
        if pattern(i,k) = 0 then
          avx := AbsVal(x(i,k));
          if avx > tol then
            pivot := 0;
            for j in 1..k-1 loop
              if pattern(i,j) = 1
               then pivot := j;
              end if;
              exit when (pivot > 0);
            end loop;
            if pivot > 0 then
              for ii in x'first(1)..i-1 loop
                x(ii,k) := x(ii,k) - x(i,k)*x(ii,pivot);
              end loop;
              for ii in i+1..x'last(1) loop
                x(ii,k) := x(ii,k) - x(i,k)*x(ii,pivot);
              end loop;
              x(i,k) := Create(integer(0));
            end if;
          end if;
        end if;
      end loop;
    end loop;
  end Reduce_to_Fit;

  procedure Reduce_to_Fit
              ( pattern : in Standard_Natural_Matrices.Matrix;
                x : in out QuadDobl_Complex_Matrices.Matrix ) is

    use QuadDobl_Complex_Numbers;
    tol : constant double_float := 1.0e-8;
    avx : quad_double;
    pivot : integer32;

  begin
    for k in pattern'range(2) loop
      for i in pattern'range(1) loop
        if pattern(i,k) = 0 then
          avx := AbsVal(x(i,k));
          if avx > tol then
            pivot := 0;
            for j in 1..k-1 loop
              if pattern(i,j) = 1
               then pivot := j;
              end if;
              exit when (pivot > 0);
            end loop;
            if pivot > 0 then
              for ii in x'first(1)..i-1 loop
                x(ii,k) := x(ii,k) - x(i,k)*x(ii,pivot);
              end loop;
              for ii in i+1..x'last(1) loop
                x(ii,k) := x(ii,k) - x(i,k)*x(ii,pivot);
              end loop;
              x(i,k) := Create(integer(0));
            end if;
          end if;
        end if;
      end loop;
    end loop;
  end Reduce_to_Fit;

  procedure Normalize_and_Reduce_to_Fit
              ( pattern : in Standard_Natural_Matrices.Matrix;
                x : in out Standard_Complex_Matrices.Matrix ) is

  begin
    Normalize_to_Fit(pattern,x);
    Reduce_to_Fit(pattern,x);
  end Normalize_and_Reduce_to_Fit;

  procedure Normalize_and_Reduce_to_Fit
              ( pattern : in Standard_Natural_Matrices.Matrix;
                x : in out DoblDobl_Complex_Matrices.Matrix ) is

  begin
    Normalize_to_Fit(pattern,x);
    Reduce_to_Fit(pattern,x);
  end Normalize_and_Reduce_to_Fit;

  procedure Normalize_and_Reduce_to_Fit
              ( pattern : in Standard_Natural_Matrices.Matrix;
                x : in out QuadDobl_Complex_Matrices.Matrix ) is

  begin
    Normalize_to_Fit(pattern,x);
    Reduce_to_Fit(pattern,x);
  end Normalize_and_Reduce_to_Fit;

  procedure Trivial_Stay_Coordinates
              ( n,k,r : in integer32;
                q,p,qr,qc,pr,pc : in Standard_Natural_Vectors.Vector;
                x : in out Standard_Complex_Vectors.Vector ) is

    plocmap : constant Standard_Natural_Matrices.Matrix(1..n,1..k)
            := Checker_Localization_Patterns.Column_Pattern(n,k,p,pr,pc);
    qlocmap : constant Standard_Natural_Matrices.Matrix(1..n,1..k)
            := Checker_Localization_Patterns.Column_Pattern(n,k,q,qr,qc);
    y : Standard_Complex_Matrices.Matrix(1..n,1..k);

  begin
    y := Map(plocmap,x);
    Inverse_Row_Transformation(r,y);
    Normalize_and_Reduce_to_Fit(qlocmap,y);
    x := Map(qlocmap,y);
  end Trivial_Stay_Coordinates;

  procedure Trivial_Stay_Coordinates
              ( n,k,r : in integer32;
                q,p,qr,qc,pr,pc : in Standard_Natural_Vectors.Vector;
                x : in out DoblDobl_Complex_Vectors.Vector ) is

    plocmap : constant Standard_Natural_Matrices.Matrix(1..n,1..k)
            := Checker_Localization_Patterns.Column_Pattern(n,k,p,pr,pc);
    qlocmap : constant Standard_Natural_Matrices.Matrix(1..n,1..k)
            := Checker_Localization_Patterns.Column_Pattern(n,k,q,qr,qc);
    y : DoblDobl_Complex_Matrices.Matrix(1..n,1..k);

  begin
    y := Map(plocmap,x);
    Inverse_Row_Transformation(r,y);
    Normalize_and_Reduce_to_Fit(qlocmap,y);
    x := Map(qlocmap,y);
  end Trivial_Stay_Coordinates;

  procedure Trivial_Stay_Coordinates
              ( n,k,r : in integer32;
                q,p,qr,qc,pr,pc : in Standard_Natural_Vectors.Vector;
                x : in out QuadDobl_Complex_Vectors.Vector ) is

    plocmap : constant Standard_Natural_Matrices.Matrix(1..n,1..k)
            := Checker_Localization_Patterns.Column_Pattern(n,k,p,pr,pc);
    qlocmap : constant Standard_Natural_Matrices.Matrix(1..n,1..k)
            := Checker_Localization_Patterns.Column_Pattern(n,k,q,qr,qc);
    y : QuadDobl_Complex_Matrices.Matrix(1..n,1..k);

  begin
    y := Map(plocmap,x);
    Inverse_Row_Transformation(r,y);
    Normalize_and_Reduce_to_Fit(qlocmap,y);
    x := Map(qlocmap,y);
  end Trivial_Stay_Coordinates;

  procedure Trivial_Stay_Coordinates
              ( n,k,r : in integer32;
                q,p,qr,qc,pr,pc : in Standard_Natural_Vectors.Vector;
                sols : in out Standard_Complex_Solutions.Solution_List) is

    use Standard_Complex_Solutions;

    plocmap : constant Standard_Natural_Matrices.Matrix(1..n,1..k)
            := Checker_Localization_Patterns.Column_Pattern(n,k,p,pr,pc);
    qlocmap : constant Standard_Natural_Matrices.Matrix(1..n,1..k)
            := Checker_Localization_Patterns.Column_Pattern(n,k,q,qr,qc);
    y : Standard_Complex_Matrices.Matrix(1..n,1..k);
    tmp : Solution_List := sols;
    ls : Link_to_Solution;

  begin
    while not Is_Null(tmp) loop
      ls := Head_Of(tmp);
      y := Map(plocmap,ls.v);
      Inverse_Row_Transformation(r,y);
      Normalize_and_Reduce_to_Fit(qlocmap,y);
      ls.v := Map(qlocmap,y);
      Set_Head(tmp,ls);
      tmp := Tail_Of(tmp);
    end loop;
  end Trivial_Stay_Coordinates;

  procedure Trivial_Stay_Coordinates
              ( n,k,r : in integer32;
                q,p,qr,qc,pr,pc : in Standard_Natural_Vectors.Vector;
                sols : in out DoblDobl_Complex_Solutions.Solution_List) is

    use DoblDobl_Complex_Solutions;

    plocmap : constant Standard_Natural_Matrices.Matrix(1..n,1..k)
            := Checker_Localization_Patterns.Column_Pattern(n,k,p,pr,pc);
    qlocmap : constant Standard_Natural_Matrices.Matrix(1..n,1..k)
            := Checker_Localization_Patterns.Column_Pattern(n,k,q,qr,qc);
    y : DoblDobl_Complex_Matrices.Matrix(1..n,1..k);
    tmp : Solution_List := sols;
    ls : Link_to_Solution;

  begin
    while not Is_Null(tmp) loop
      ls := Head_Of(tmp);
      y := Map(plocmap,ls.v);
      Inverse_Row_Transformation(r,y);
      Normalize_and_Reduce_to_Fit(qlocmap,y);
      ls.v := Map(qlocmap,y);
      Set_Head(tmp,ls);
      tmp := Tail_Of(tmp);
    end loop;
  end Trivial_Stay_Coordinates;

  procedure Trivial_Stay_Coordinates
              ( n,k,r : in integer32;
                q,p,qr,qc,pr,pc : in Standard_Natural_Vectors.Vector;
                sols : in out QuadDobl_Complex_Solutions.Solution_List) is

    use QuadDobl_Complex_Solutions;

    plocmap : constant Standard_Natural_Matrices.Matrix(1..n,1..k)
            := Checker_Localization_Patterns.Column_Pattern(n,k,p,pr,pc);
    qlocmap : constant Standard_Natural_Matrices.Matrix(1..n,1..k)
            := Checker_Localization_Patterns.Column_Pattern(n,k,q,qr,qc);
    y : QuadDobl_Complex_Matrices.Matrix(1..n,1..k);
    tmp : Solution_List := sols;
    ls : Link_to_Solution;

  begin
    while not Is_Null(tmp) loop
      ls := Head_Of(tmp);
      y := Map(plocmap,ls.v);
      Inverse_Row_Transformation(r,y);
      Normalize_and_Reduce_to_Fit(qlocmap,y);
      ls.v := Map(qlocmap,y);
      Set_Head(tmp,ls);
      tmp := Tail_Of(tmp);
    end loop;
  end Trivial_Stay_Coordinates;

  procedure Trivial_Stay_Coordinates
              ( file : in file_type; n,k,r : in integer32;
                q,p,qr,qc,pr,pc : in Standard_Natural_Vectors.Vector;
                x : in out Standard_Complex_Vectors.Vector ) is

    plocmap : constant Standard_Natural_Matrices.Matrix(1..n,1..k)
            := Checker_Localization_Patterns.Column_Pattern(n,k,p,pr,pc);
    qlocmap : constant Standard_Natural_Matrices.Matrix(1..n,1..k)
            := Checker_Localization_Patterns.Column_Pattern(n,k,q,qr,qc);
    inrowrp1 : integer32 := 0;
    y : Standard_Complex_Matrices.Matrix(1..n,1..k);

  begin
    put(file,"Trivial Stay with critical row = "); put(file,r,1);
    put_line(file,".");
   -- put_line(file,"The previous localization map : "); put(file,plocmap);
   -- put_line(file,"The current localization map : "); put(file,qlocmap);
    put(file,"rows of white checkers : "); put(file,qr); new_line(file);
    for j in qlocmap'range(2) loop
      if qlocmap(r+1,j) = 1 then
        inrowrp1 := j;
      end if;
    end loop;
    if inrowrp1 = 0 then
      put(file,"no white checker in row ");
    else
      put(file,"white checker "); 
      put(file,inrowrp1,1); put(file," in row ");
    end if;
    put(file,r+1,1); new_line(file);
    y := Map(plocmap,x);
   -- put_line(file,"the given solution plane :"); put(file,y,3);
    Inverse_Row_Transformation(r,y);
   -- put_line(file,"after the inverse transformation :"); put(file,y,3);
    Normalize_and_Reduce_to_Fit(qlocmap,y);
   -- put_line(file,"The transformed plane :"); put(file,y,3);
    x := Map(qlocmap,y);
  end Trivial_Stay_Coordinates;

  procedure Trivial_Stay_Coordinates
              ( file : in file_type; n,k,r : in integer32;
                q,p,qr,qc,pr,pc : in Standard_Natural_Vectors.Vector;
                x : in out DoblDobl_Complex_Vectors.Vector ) is

    plocmap : constant Standard_Natural_Matrices.Matrix(1..n,1..k)
            := Checker_Localization_Patterns.Column_Pattern(n,k,p,pr,pc);
    qlocmap : constant Standard_Natural_Matrices.Matrix(1..n,1..k)
            := Checker_Localization_Patterns.Column_Pattern(n,k,q,qr,qc);
    inrowrp1 : integer32 := 0;
    y : DoblDobl_Complex_Matrices.Matrix(1..n,1..k);

  begin
    put(file,"Trivial Stay with critical row = "); put(file,r,1);
    put_line(file,".");
   -- put_line(file,"The previous localization map : "); put(file,plocmap);
   -- put_line(file,"The current localization map : "); put(file,qlocmap);
    put(file,"rows of white checkers : "); put(file,qr); new_line(file);
    for j in qlocmap'range(2) loop
      if qlocmap(r+1,j) = 1 then
        inrowrp1 := j;
      end if;
    end loop;
    if inrowrp1 = 0 then
      put(file,"no white checker in row ");
    else
      put(file,"white checker "); 
      put(file,inrowrp1,1); put(file," in row ");
    end if;
    put(file,r+1,1); new_line(file);
    y := Map(plocmap,x);
   -- put_line(file,"the given solution plane :"); put(file,y,3);
    Inverse_Row_Transformation(r,y);
   -- put_line(file,"after the inverse transformation :"); put(file,y,3);
    Normalize_and_Reduce_to_Fit(qlocmap,y);
   -- put_line(file,"The transformed plane :"); put(file,y,3);
    x := Map(qlocmap,y);
  end Trivial_Stay_Coordinates;

  procedure Trivial_Stay_Coordinates
              ( file : in file_type; n,k,r : in integer32;
                q,p,qr,qc,pr,pc : in Standard_Natural_Vectors.Vector;
                x : in out QuadDobl_Complex_Vectors.Vector ) is

    plocmap : constant Standard_Natural_Matrices.Matrix(1..n,1..k)
            := Checker_Localization_Patterns.Column_Pattern(n,k,p,pr,pc);
    qlocmap : constant Standard_Natural_Matrices.Matrix(1..n,1..k)
            := Checker_Localization_Patterns.Column_Pattern(n,k,q,qr,qc);
    inrowrp1 : integer32 := 0;
    y : QuadDobl_Complex_Matrices.Matrix(1..n,1..k);

  begin
    put(file,"Trivial Stay with critical row = "); put(file,r,1);
    put_line(file,".");
   -- put_line(file,"The previous localization map : "); put(file,plocmap);
   -- put_line(file,"The current localization map : "); put(file,qlocmap);
    put(file,"rows of white checkers : "); put(file,qr); new_line(file);
    for j in qlocmap'range(2) loop
      if qlocmap(r+1,j) = 1 then
        inrowrp1 := j;
      end if;
    end loop;
    if inrowrp1 = 0 then
      put(file,"no white checker in row ");
    else
      put(file,"white checker "); 
      put(file,inrowrp1,1); put(file," in row ");
    end if;
    put(file,r+1,1); new_line(file);
    y := Map(plocmap,x);
   -- put_line(file,"the given solution plane :"); put(file,y,3);
    Inverse_Row_Transformation(r,y);
   -- put_line(file,"after the inverse transformation :"); put(file,y,3);
    Normalize_and_Reduce_to_Fit(qlocmap,y);
   -- put_line(file,"The transformed plane :"); put(file,y,3);
    x := Map(qlocmap,y);
  end Trivial_Stay_Coordinates;

  procedure Trivial_Stay_Coordinates
              ( file : in file_type; n,k,r : in integer32;
                q,p,qr,qc,pr,pc : in Standard_Natural_Vectors.Vector;
                sols : in out Standard_Complex_Solutions.Solution_List ) is

    use Standard_Complex_Solutions;

    plocmap : constant Standard_Natural_Matrices.Matrix(1..n,1..k)
            := Checker_Localization_Patterns.Column_Pattern(n,k,p,pr,pc);
    qlocmap : constant Standard_Natural_Matrices.Matrix(1..n,1..k)
            := Checker_Localization_Patterns.Column_Pattern(n,k,q,qr,qc);
    inrowrp1 : integer32 := 0;
    y : Standard_Complex_Matrices.Matrix(1..n,1..k);
    tmp : Solution_List := sols;
    ls : Link_to_Solution;

  begin
    put(file,"Trivial Stay with critical row = "); put(file,r,1);
    put_line(file,".");
   -- put_line(file,"The previous localization map : "); put(file,plocmap);
   -- put_line(file,"The current localization map : "); put(file,qlocmap);
    put(file,"rows of white checkers : "); put(file,qr); new_line(file);
    for j in qlocmap'range(2) loop
      if qlocmap(r+1,j) = 1 then
        inrowrp1 := j;
      end if;
    end loop;
    if inrowrp1 = 0 then
      put(file,"no white checker in row ");
    else
      put(file,"white checker "); 
      put(file,inrowrp1,1); put(file," in row ");
    end if;
    put(file,r+1,1); new_line(file);
    while not Is_Null(tmp) loop
      ls := Head_Of(tmp);
      y := Map(plocmap,ls.v);
     -- put_line(file,"the given solution plane :"); put(file,y,3);
      Inverse_Row_Transformation(r,y);
     -- put_line(file,"after the inverse transformation :"); put(file,y,3);
      Normalize_and_Reduce_to_Fit(qlocmap,y);
     -- put_line(file,"The transformed plane :"); put(file,y,3);
      ls.v := Map(qlocmap,y);
      Set_Head(tmp,ls);
      tmp := Tail_Of(tmp);
    end loop;
  end Trivial_Stay_Coordinates;

  procedure Trivial_Stay_Coordinates
              ( file : in file_type; n,k,r : in integer32;
                q,p,qr,qc,pr,pc : in Standard_Natural_Vectors.Vector;
                sols : in out DoblDobl_Complex_Solutions.Solution_List ) is

    use DoblDobl_Complex_Solutions;

    plocmap : constant Standard_Natural_Matrices.Matrix(1..n,1..k)
            := Checker_Localization_Patterns.Column_Pattern(n,k,p,pr,pc);
    qlocmap : constant Standard_Natural_Matrices.Matrix(1..n,1..k)
            := Checker_Localization_Patterns.Column_Pattern(n,k,q,qr,qc);
    inrowrp1 : integer32 := 0;
    y : DoblDobl_Complex_Matrices.Matrix(1..n,1..k);
    tmp : Solution_List := sols;
    ls : Link_to_Solution;

  begin
    put(file,"Trivial Stay with critical row = "); put(file,r,1);
    put_line(file,".");
   -- put_line(file,"The previous localization map : "); put(file,plocmap);
   -- put_line(file,"The current localization map : "); put(file,qlocmap);
    put(file,"rows of white checkers : "); put(file,qr); new_line(file);
    for j in qlocmap'range(2) loop
      if qlocmap(r+1,j) = 1 then
        inrowrp1 := j;
      end if;
    end loop;
    if inrowrp1 = 0 then
      put(file,"no white checker in row ");
    else
      put(file,"white checker "); 
      put(file,inrowrp1,1); put(file," in row ");
    end if;
    put(file,r+1,1); new_line(file);
    while not Is_Null(tmp) loop
      ls := Head_Of(tmp);
      y := Map(plocmap,ls.v);
     -- put_line(file,"the given solution plane :"); put(file,y,3);
      Inverse_Row_Transformation(r,y);
     -- put_line(file,"after the inverse transformation :"); put(file,y,3);
      Normalize_and_Reduce_to_Fit(qlocmap,y);
     -- put_line(file,"The transformed plane :"); put(file,y,3);
      ls.v := Map(qlocmap,y);
      Set_Head(tmp,ls);
      tmp := Tail_Of(tmp);
    end loop;
  end Trivial_Stay_Coordinates;

  procedure Trivial_Stay_Coordinates
              ( file : in file_type; n,k,r : in integer32;
                q,p,qr,qc,pr,pc : in Standard_Natural_Vectors.Vector;
                sols : in out QuadDobl_Complex_Solutions.Solution_List ) is

    use QuadDobl_Complex_Solutions;

    plocmap : constant Standard_Natural_Matrices.Matrix(1..n,1..k)
            := Checker_Localization_Patterns.Column_Pattern(n,k,p,pr,pc);
    qlocmap : constant Standard_Natural_Matrices.Matrix(1..n,1..k)
            := Checker_Localization_Patterns.Column_Pattern(n,k,q,qr,qc);
    inrowrp1 : integer32 := 0;
    y : QuadDobl_Complex_Matrices.Matrix(1..n,1..k);
    tmp : Solution_List := sols;
    ls : Link_to_Solution;

  begin
    put(file,"Trivial Stay with critical row = "); put(file,r,1);
    put_line(file,".");
   -- put_line(file,"The previous localization map : "); put(file,plocmap);
   -- put_line(file,"The current localization map : "); put(file,qlocmap);
    put(file,"rows of white checkers : "); put(file,qr); new_line(file);
    for j in qlocmap'range(2) loop
      if qlocmap(r+1,j) = 1 then
        inrowrp1 := j;
      end if;
    end loop;
    if inrowrp1 = 0 then
      put(file,"no white checker in row ");
    else
      put(file,"white checker "); 
      put(file,inrowrp1,1); put(file," in row ");
    end if;
    put(file,r+1,1); new_line(file);
    while not Is_Null(tmp) loop
      ls := Head_Of(tmp);
      y := Map(plocmap,ls.v);
     -- put_line(file,"the given solution plane :"); put(file,y,3);
      Inverse_Row_Transformation(r,y);
     -- put_line(file,"after the inverse transformation :"); put(file,y,3);
      Normalize_and_Reduce_to_Fit(qlocmap,y);
     -- put_line(file,"The transformed plane :"); put(file,y,3);
      ls.v := Map(qlocmap,y);
      Set_Head(tmp,ls);
      tmp := Tail_Of(tmp);
    end loop;
  end Trivial_Stay_Coordinates;

  function Eval ( m : Standard_Complex_Poly_Matrices.Matrix;
                  x : Standard_Complex_Vectors.Vector )
                return Standard_Complex_Matrices.Matrix is

  -- DESCRIPTION :
  --   Returns the value of m evaluated at  x.

    res : Standard_Complex_Matrices.Matrix(m'range(1),m'range(2));

  begin
    for i in res'range(1) loop
      for j in res'range(2) loop
        res(i,j) := Eval(m(i,j),x);
      end loop;
    end loop;
    return res;
  end Eval;

  function Eval ( m : DoblDobl_Complex_Poly_Matrices.Matrix;
                  x : DoblDobl_Complex_Vectors.Vector )
                return DoblDobl_Complex_Matrices.Matrix is

  -- DESCRIPTION :
  --   Returns the value of m evaluated at  x.

    res : DoblDobl_Complex_Matrices.Matrix(m'range(1),m'range(2));

  begin
    for i in res'range(1) loop
      for j in res'range(2) loop
        res(i,j) := Eval(m(i,j),x);
      end loop;
    end loop;
    return res;
  end Eval;

  function Eval ( m : QuadDobl_Complex_Poly_Matrices.Matrix;
                  x : QuadDobl_Complex_Vectors.Vector )
                return QuadDobl_Complex_Matrices.Matrix is

  -- DESCRIPTION :
  --   Returns the value of m evaluated at  x.

    res : QuadDobl_Complex_Matrices.Matrix(m'range(1),m'range(2));

  begin
    for i in res'range(1) loop
      for j in res'range(2) loop
        res(i,j) := Eval(m(i,j),x);
      end loop;
    end loop;
    return res;
  end Eval;

  procedure Homotopy_Stay_Coordinates
              ( n,k,r : in integer32;
                p,rows,cols : in Standard_Natural_Vectors.Vector;
                mf : in Standard_Complex_Matrices.Matrix;
                xtm : in Standard_Complex_Poly_Matrices.Matrix;
                x : in out Standard_Complex_Vectors.Vector ) is

    use Standard_Complex_Numbers;

    locmap : constant Standard_Natural_Matrices.Matrix(1..n,1..k)
           := Checker_Localization_Patterns.Column_Pattern(n,k,p,rows,cols);
    y : Standard_Complex_Matrices.Matrix(1..n,1..k);
    xt : Standard_Complex_Vectors.Vector(x'first..x'last+1);

  begin
    xt(x'range) := x;
    xt(xt'last) := Create(1.0);
    y := Eval(xtm,xt);
    Inverse_Row_Transformation(mf,y);
    Normalize_and_Reduce_to_Fit(locmap,y);
    x := Map(locmap,y);
  end Homotopy_Stay_Coordinates;

  procedure Homotopy_Stay_Coordinates
              ( n,k,r : in integer32;
                p,rows,cols : in Standard_Natural_Vectors.Vector;
                mf : in DoblDobl_Complex_Matrices.Matrix;
                xtm : in DoblDobl_Complex_Poly_Matrices.Matrix;
                x : in out DoblDobl_Complex_Vectors.Vector ) is

    use DoblDobl_Complex_Numbers;

    locmap : constant Standard_Natural_Matrices.Matrix(1..n,1..k)
           := Checker_Localization_Patterns.Column_Pattern(n,k,p,rows,cols);
    y : DoblDobl_Complex_Matrices.Matrix(1..n,1..k);
    xt : DoblDobl_Complex_Vectors.Vector(x'first..x'last+1);

  begin
    xt(x'range) := x;
    xt(xt'last) := Create(integer(1));
    y := Eval(xtm,xt);
    Inverse_Row_Transformation(mf,y);
    Normalize_and_Reduce_to_Fit(locmap,y);
    x := Map(locmap,y);
  end Homotopy_Stay_Coordinates;

  procedure Homotopy_Stay_Coordinates
              ( n,k,r : in integer32;
                p,rows,cols : in Standard_Natural_Vectors.Vector;
                mf : in QuadDobl_Complex_Matrices.Matrix;
                xtm : in QuadDobl_Complex_Poly_Matrices.Matrix;
                x : in out QuadDobl_Complex_Vectors.Vector ) is

    use QuadDobl_Complex_Numbers;

    locmap : constant Standard_Natural_Matrices.Matrix(1..n,1..k)
           := Checker_Localization_Patterns.Column_Pattern(n,k,p,rows,cols);
    y : QuadDobl_Complex_Matrices.Matrix(1..n,1..k);
    xt : QuadDobl_Complex_Vectors.Vector(x'first..x'last+1);

  begin
    xt(x'range) := x;
    xt(xt'last) := Create(integer(1));
    y := Eval(xtm,xt);
    Inverse_Row_Transformation(mf,y);
    Normalize_and_Reduce_to_Fit(locmap,y);
    x := Map(locmap,y);
  end Homotopy_Stay_Coordinates;

  procedure Homotopy_Stay_Coordinates
              ( n,k,r : in integer32;
                p,rows,cols : in Standard_Natural_Vectors.Vector;
                mf : in Standard_Complex_Matrices.Matrix;
                xtm : in Standard_Complex_Poly_Matrices.Matrix;
                sols : in out Standard_Complex_Solutions.Solution_List) is

    use Standard_Complex_Numbers;
    use Standard_Complex_Solutions;

    locmap : constant Standard_Natural_Matrices.Matrix(1..n,1..k)
           := Checker_Localization_Patterns.Column_Pattern(n,k,p,rows,cols);
    y : Standard_Complex_Matrices.Matrix(1..n,1..k);
    tmp : Solution_List := sols;
    ls : Link_to_Solution;

  begin
    while not Is_Null(tmp) loop
      ls := Head_Of(tmp);
      declare
        xt : Standard_Complex_Vectors.Vector(ls.v'first..ls.v'last+1);
      begin
        xt(ls.v'range) := ls.v;
        xt(xt'last) := Create(1.0);
        y := Eval(xtm,xt);
      end;
      Inverse_Row_Transformation(mf,y);
      Normalize_and_Reduce_to_Fit(locmap,y);
      ls.v := Map(locmap,y);
      Set_Head(tmp,ls);
      tmp := Tail_Of(tmp);
    end loop;
  end Homotopy_Stay_Coordinates;

  procedure Homotopy_Stay_Coordinates
              ( n,k,r : in integer32;
                p,rows,cols : in Standard_Natural_Vectors.Vector;
                mf : in DoblDobl_Complex_Matrices.Matrix;
                xtm : in DoblDobl_Complex_Poly_Matrices.Matrix;
                sols : in out DoblDobl_Complex_Solutions.Solution_List) is

    use DoblDobl_Complex_Numbers;
    use DoblDobl_Complex_Solutions;

    locmap : constant Standard_Natural_Matrices.Matrix(1..n,1..k)
           := Checker_Localization_Patterns.Column_Pattern(n,k,p,rows,cols);
    y : DoblDobl_Complex_Matrices.Matrix(1..n,1..k);
    tmp : Solution_List := sols;
    ls : Link_to_Solution;

  begin
    while not Is_Null(tmp) loop
      ls := Head_Of(tmp);
      declare
        xt : DoblDobl_Complex_Vectors.Vector(ls.v'first..ls.v'last+1);
      begin
        xt(ls.v'range) := ls.v;
        xt(xt'last) := Create(integer(1));
        y := Eval(xtm,xt);
      end;
      Inverse_Row_Transformation(mf,y);
      Normalize_and_Reduce_to_Fit(locmap,y);
      ls.v := Map(locmap,y);
      Set_Head(tmp,ls);
      tmp := Tail_Of(tmp);
    end loop;
  end Homotopy_Stay_Coordinates;

  procedure Homotopy_Stay_Coordinates
              ( n,k,r : in integer32;
                p,rows,cols : in Standard_Natural_Vectors.Vector;
                mf : in QuadDobl_Complex_Matrices.Matrix;
                xtm : in QuadDobl_Complex_Poly_Matrices.Matrix;
                sols : in out QuadDobl_Complex_Solutions.Solution_List) is

    use QuadDobl_Complex_Numbers;
    use QuadDobl_Complex_Solutions;

    locmap : constant Standard_Natural_Matrices.Matrix(1..n,1..k)
           := Checker_Localization_Patterns.Column_Pattern(n,k,p,rows,cols);
    y : QuadDobl_Complex_Matrices.Matrix(1..n,1..k);
    tmp : Solution_List := sols;
    ls : Link_to_Solution;

  begin
    while not Is_Null(tmp) loop
      ls := Head_Of(tmp);
      declare
        xt : QuadDobl_Complex_Vectors.Vector(ls.v'first..ls.v'last+1);
      begin
        xt(ls.v'range) := ls.v;
        xt(xt'last) := Create(integer(1));
        y := Eval(xtm,xt);
      end;
      Inverse_Row_Transformation(mf,y);
      Normalize_and_Reduce_to_Fit(locmap,y);
      ls.v := Map(locmap,y);
      Set_Head(tmp,ls);
      tmp := Tail_Of(tmp);
    end loop;
  end Homotopy_Stay_Coordinates;

  procedure Homotopy_Stay_Coordinates
              ( file : in file_type; n,k,r : in integer32;
                p,rows,cols : in Standard_Natural_Vectors.Vector;
                mf : in Standard_Complex_Matrices.Matrix;
                xtm : in Standard_Complex_Poly_Matrices.Matrix;
                x : in out Standard_Complex_Vectors.Vector ) is

    use Standard_Complex_Numbers;

    locmap : constant Standard_Natural_Matrices.Matrix(1..n,1..k)
           := Checker_Localization_Patterns.Column_Pattern(n,k,p,rows,cols);
    y : Standard_Complex_Matrices.Matrix(1..n,1..k);
    xt : Standard_Complex_Vectors.Vector(x'first..x'last+1);

  begin
    put(file,"Homotopy Stay with critical row = "); put(file,r,1);
    put_line(file,".");
   -- put_line(file,"The localization map : "); put(file,locmap);
   -- put_line(file,"The matrix xtm : "); put(file,xtm);
    xt(x'range) := x;
    xt(xt'last) := Create(1.0);
    put_line(file,"The vector xt : "); put_line(file,xt);
    y := Eval(xtm,xt);
   -- put_line(file,"The matrix xtm evaluated at the solution : ");
   -- put(file,y,2);
    Inverse_Row_Transformation(mf,y);
   -- put_line(file,"after the inverse transformation :"); put(file,y,3);
    Normalize_and_Reduce_to_Fit(locmap,y);
   -- put_line(file,"The transformed plane :"); put(file,y,3);
    x := Map(locmap,y);
  end Homotopy_Stay_Coordinates;

  procedure Homotopy_Stay_Coordinates
              ( file : in file_type; n,k,r : in integer32;
                p,rows,cols : in Standard_Natural_Vectors.Vector;
                mf : in DoblDobl_Complex_Matrices.Matrix;
                xtm : in DoblDobl_Complex_Poly_Matrices.Matrix;
                x : in out DoblDobl_Complex_Vectors.Vector ) is

    use DoblDobl_Complex_Numbers;

    locmap : constant Standard_Natural_Matrices.Matrix(1..n,1..k)
           := Checker_Localization_Patterns.Column_Pattern(n,k,p,rows,cols);
    y : DoblDobl_Complex_Matrices.Matrix(1..n,1..k);
    xt : DoblDobl_Complex_Vectors.Vector(x'first..x'last+1);

  begin
    put(file,"Homotopy Stay with critical row = "); put(file,r,1);
    put_line(file,".");
   -- put_line(file,"The localization map : "); put(file,locmap);
   -- put_line(file,"The matrix xtm : "); put(file,xtm);
    xt(x'range) := x;
    xt(xt'last) := Create(integer(1));
    put_line(file,"The vector xt : "); put_line(file,xt);
    y := Eval(xtm,xt);
   -- put_line(file,"The matrix xtm evaluated at the solution : ");
   -- put(file,y,2);
    Inverse_Row_Transformation(mf,y);
   -- put_line(file,"after the inverse transformation :"); put(file,y,3);
    Normalize_and_Reduce_to_Fit(locmap,y);
   -- put_line(file,"The transformed plane :"); put(file,y,3);
    x := Map(locmap,y);
  end Homotopy_Stay_Coordinates;

  procedure Homotopy_Stay_Coordinates
              ( file : in file_type; n,k,r : in integer32;
                p,rows,cols : in Standard_Natural_Vectors.Vector;
                mf : in QuadDobl_Complex_Matrices.Matrix;
                xtm : in QuadDobl_Complex_Poly_Matrices.Matrix;
                x : in out QuadDobl_Complex_Vectors.Vector ) is

    use QuadDobl_Complex_Numbers;

    locmap : constant Standard_Natural_Matrices.Matrix(1..n,1..k)
           := Checker_Localization_Patterns.Column_Pattern(n,k,p,rows,cols);
    y : QuadDobl_Complex_Matrices.Matrix(1..n,1..k);
    xt : QuadDobl_Complex_Vectors.Vector(x'first..x'last+1);

  begin
    put(file,"Homotopy Stay with critical row = "); put(file,r,1);
    put_line(file,".");
   -- put_line(file,"The localization map : "); put(file,locmap);
   -- put_line(file,"The matrix xtm : "); put(file,xtm);
    xt(x'range) := x;
    xt(xt'last) := Create(integer(1));
    put_line(file,"The vector xt : "); put_line(file,xt);
    y := Eval(xtm,xt);
   -- put_line(file,"The matrix xtm evaluated at the solution : ");
   -- put(file,y,2);
    Inverse_Row_Transformation(mf,y);
   -- put_line(file,"after the inverse transformation :"); put(file,y,3);
    Normalize_and_Reduce_to_Fit(locmap,y);
   -- put_line(file,"The transformed plane :"); put(file,y,3);
    x := Map(locmap,y);
  end Homotopy_Stay_Coordinates;

  procedure Homotopy_Stay_Coordinates
              ( file : in file_type; n,k,r : in integer32;
                p,rows,cols : in Standard_Natural_Vectors.Vector;
                mf : in Standard_Complex_Matrices.Matrix;
                xtm : in Standard_Complex_Poly_Matrices.Matrix;
                sols : in out Standard_Complex_Solutions.Solution_List ) is

    use Standard_Complex_Numbers;
    use Standard_Complex_Solutions;

    locmap : constant Standard_Natural_Matrices.Matrix(1..n,1..k)
           := Checker_Localization_Patterns.Column_Pattern(n,k,p,rows,cols);
    y : Standard_Complex_Matrices.Matrix(1..n,1..k);
    tmp : Solution_List := sols;
    ls : Link_to_Solution;

  begin
    put(file,"Homotopy Stay with critical row = "); put(file,r,1);
    put_line(file,".");
   -- put_line(file,"The localization map : "); put(file,locmap);
   -- put_line(file,"The matrix xtm : "); put(file,xtm);
    while not Is_Null(tmp) loop
      ls := Head_Of(tmp);
      declare
        xt : Standard_Complex_Vectors.Vector(ls.v'first..ls.v'last+1);
      begin
        xt(ls.v'range) := ls.v;
        xt(xt'last) := Create(1.0);
        put_line(file,"The vector xt : "); put_line(file,xt);
        y := Eval(xtm,xt);
      end;
     -- put_line(file,"The matrix xtm evaluated at the solution : ");
     -- put(file,y,2);
      Inverse_Row_Transformation(mf,y);
     -- put_line(file,"after the inverse transformation :"); put(file,y,3);
      Normalize_and_Reduce_to_Fit(locmap,y);
     -- put_line(file,"The transformed plane :"); put(file,y,3);
      ls.v := Map(locmap,y);
     -- Set_Head(tmp,ls);
      tmp := Tail_Of(tmp);
    end loop;
    put_line(file,"The transformed solution list :");
    put(file,Length_Of(sols),natural32(Head_Of(sols).n),sols);
  end Homotopy_Stay_Coordinates;

  procedure Homotopy_Stay_Coordinates
              ( file : in file_type; n,k,r : in integer32;
                p,rows,cols : in Standard_Natural_Vectors.Vector;
                mf : in DoblDobl_Complex_Matrices.Matrix;
                xtm : in DoblDobl_Complex_Poly_Matrices.Matrix;
                sols : in out DoblDobl_Complex_Solutions.Solution_List ) is

    use DoblDobl_Complex_Numbers;
    use DoblDobl_Complex_Solutions;

    locmap : constant Standard_Natural_Matrices.Matrix(1..n,1..k)
           := Checker_Localization_Patterns.Column_Pattern(n,k,p,rows,cols);
    y : DoblDobl_Complex_Matrices.Matrix(1..n,1..k);
    tmp : Solution_List := sols;
    ls : Link_to_Solution;

  begin
    put(file,"Homotopy Stay with critical row = "); put(file,r,1);
    put_line(file,".");
   -- put_line(file,"The localization map : "); put(file,locmap);
   -- put_line(file,"The matrix xtm : "); put(file,xtm);
    while not Is_Null(tmp) loop
      ls := Head_Of(tmp);
      declare
        xt : DoblDobl_Complex_Vectors.Vector(ls.v'first..ls.v'last+1);
      begin
        xt(ls.v'range) := ls.v;
        xt(xt'last) := Create(integer(1));
        put_line(file,"The vector xt : "); put_line(file,xt);
        y := Eval(xtm,xt);
      end;
     -- put_line(file,"The matrix xtm evaluated at the solution : ");
     -- put(file,y,2);
      Inverse_Row_Transformation(mf,y);
     -- put_line(file,"after the inverse transformation :"); put(file,y,3);
      Normalize_and_Reduce_to_Fit(locmap,y);
     -- put_line(file,"The transformed plane :"); put(file,y,3);
      ls.v := Map(locmap,y);
     -- Set_Head(tmp,ls);
      tmp := Tail_Of(tmp);
    end loop;
    put_line(file,"The transformed solution list :");
    put(file,Length_Of(sols),natural32(Head_Of(sols).n),sols);
  end Homotopy_Stay_Coordinates;

  procedure Homotopy_Stay_Coordinates
              ( file : in file_type; n,k,r : in integer32;
                p,rows,cols : in Standard_Natural_Vectors.Vector;
                mf : in QuadDobl_Complex_Matrices.Matrix;
                xtm : in QuadDobl_Complex_Poly_Matrices.Matrix;
                sols : in out QuadDobl_Complex_Solutions.Solution_List ) is

    use QuadDobl_Complex_Numbers;
    use QuadDobl_Complex_Solutions;

    locmap : constant Standard_Natural_Matrices.Matrix(1..n,1..k)
           := Checker_Localization_Patterns.Column_Pattern(n,k,p,rows,cols);
    y : QuadDobl_Complex_Matrices.Matrix(1..n,1..k);
    tmp : Solution_List := sols;
    ls : Link_to_Solution;

  begin
    put(file,"Homotopy Stay with critical row = "); put(file,r,1);
    put_line(file,".");
   -- put_line(file,"The localization map : "); put(file,locmap);
   -- put_line(file,"The matrix xtm : "); put(file,xtm);
    while not Is_Null(tmp) loop
      ls := Head_Of(tmp);
      declare
        xt : QuadDobl_Complex_Vectors.Vector(ls.v'first..ls.v'last+1);
      begin
        xt(ls.v'range) := ls.v;
        xt(xt'last) := Create(integer(1));
        put_line(file,"The vector xt : "); put_line(file,xt);
        y := Eval(xtm,xt);
      end;
     -- put_line(file,"The matrix xtm evaluated at the solution : ");
     -- put(file,y,2);
      Inverse_Row_Transformation(mf,y);
     -- put_line(file,"after the inverse transformation :"); put(file,y,3);
      Normalize_and_Reduce_to_Fit(locmap,y);
     -- put_line(file,"The transformed plane :"); put(file,y,3);
      ls.v := Map(locmap,y);
     -- Set_Head(tmp,ls);
      tmp := Tail_Of(tmp);
    end loop;
    put_line(file,"The transformed solution list :");
    put(file,Length_Of(sols),natural32(Head_Of(sols).n),sols);
  end Homotopy_Stay_Coordinates;

  procedure Update_Swap_Column
              ( x : in out Standard_Complex_Matrices.Matrix;
                s : in integer32 ) is

  -- DESCRIPTION :
  --   Replaces column s+1 by the sum of columns s and s+1.

    use Standard_Complex_Numbers;

  begin
    for i in x'range(1) loop
      x(i,s+1) := x(i,s+1) + x(i,s);
    end loop;
  end Update_Swap_Column;

  procedure Update_Swap_Column
              ( x : in out DoblDobl_Complex_Matrices.Matrix;
                s : in integer32 ) is

  -- DESCRIPTION :
  --   Replaces column s+1 by the sum of columns s and s+1.

    use DoblDobl_Complex_Numbers;

  begin
    for i in x'range(1) loop
      x(i,s+1) := x(i,s+1) + x(i,s);
    end loop;
  end Update_Swap_Column;

  procedure Update_Swap_Column
              ( x : in out QuadDobl_Complex_Matrices.Matrix;
                s : in integer32 ) is

  -- DESCRIPTION :
  --   Replaces column s+1 by the sum of columns s and s+1.

    use QuadDobl_Complex_Numbers;

  begin
    for i in x'range(1) loop
      x(i,s+1) := x(i,s+1) + x(i,s);
    end loop;
  end Update_Swap_Column;

  procedure First_Swap_Coordinates
              ( n,k,r,big_r,dc,s : in integer32;
                q,p,qr,qc,pr,pc : in Standard_Natural_Vectors.Vector;
                mf : in Standard_Complex_Matrices.Matrix;
                xtm : in Standard_Complex_Poly_Matrices.Matrix;
                x : in out Standard_Complex_Vectors.Vector ) is

    use Standard_Complex_Numbers;

    qlocmap : constant Standard_Natural_Matrices.Matrix(1..n,1..k)
            := Checker_Localization_Patterns.Column_Pattern(n,k,q,qr,qc);
    y : Standard_Complex_Matrices.Matrix(1..n,1..k);
    xt : Standard_Complex_Vectors.Vector(x'first..x'last+1);

  begin
    xt(x'range) := x;
    xt(xt'last) := Create(1.0);
    y := Eval(xtm,xt);
    Inverse_Row_Transformation(mf,y);
    Update_Swap_Column(y,s);
    Normalize_and_Reduce_to_Fit(qlocmap,y);
    x := Map(qlocmap,y);
  end First_Swap_Coordinates;

  procedure First_Swap_Coordinates
              ( n,k,r,big_r,dc,s : in integer32;
                q,p,qr,qc,pr,pc : in Standard_Natural_Vectors.Vector;
                mf : in DoblDobl_Complex_Matrices.Matrix;
                xtm : in DoblDobl_Complex_Poly_Matrices.Matrix;
                x : in out DoblDobl_Complex_Vectors.Vector ) is

    use DoblDobl_Complex_Numbers;

    qlocmap : constant Standard_Natural_Matrices.Matrix(1..n,1..k)
            := Checker_Localization_Patterns.Column_Pattern(n,k,q,qr,qc);
    y : DoblDobl_Complex_Matrices.Matrix(1..n,1..k);
    xt : DoblDobl_Complex_Vectors.Vector(x'first..x'last+1);

  begin
    xt(x'range) := x;
    xt(xt'last) := Create(integer(1));
    y := Eval(xtm,xt);
    Inverse_Row_Transformation(mf,y);
    Update_Swap_Column(y,s);
    Normalize_and_Reduce_to_Fit(qlocmap,y);
    x := Map(qlocmap,y);
  end First_Swap_Coordinates;

  procedure First_Swap_Coordinates
              ( n,k,r,big_r,dc,s : in integer32;
                q,p,qr,qc,pr,pc : in Standard_Natural_Vectors.Vector;
                mf : in QuadDobl_Complex_Matrices.Matrix;
                xtm : in QuadDobl_Complex_Poly_Matrices.Matrix;
                x : in out QuadDobl_Complex_Vectors.Vector ) is

    use QuadDobl_Complex_Numbers;

    qlocmap : constant Standard_Natural_Matrices.Matrix(1..n,1..k)
            := Checker_Localization_Patterns.Column_Pattern(n,k,q,qr,qc);
    y : QuadDobl_Complex_Matrices.Matrix(1..n,1..k);
    xt : QuadDobl_Complex_Vectors.Vector(x'first..x'last+1);

  begin
    xt(x'range) := x;
    xt(xt'last) := Create(integer(1));
    y := Eval(xtm,xt);
    Inverse_Row_Transformation(mf,y);
    Update_Swap_Column(y,s);
    Normalize_and_Reduce_to_Fit(qlocmap,y);
    x := Map(qlocmap,y);
  end First_Swap_Coordinates;

  procedure First_Swap_Coordinates
              ( n,k,r,big_r,dc,s : in integer32;
                q,p,qr,qc,pr,pc : in Standard_Natural_Vectors.Vector;
                mf : in Standard_Complex_Matrices.Matrix;
                xtm : in Standard_Complex_Poly_Matrices.Matrix;
                sols : in out Standard_Complex_Solutions.Solution_List ) is

    use Standard_Complex_Numbers;
    use Standard_Complex_Solutions;

    qlocmap : constant Standard_Natural_Matrices.Matrix(1..n,1..k)
            := Checker_Localization_Patterns.Column_Pattern(n,k,q,qr,qc);
    y : Standard_Complex_Matrices.Matrix(1..n,1..k);
    tmp : Solution_List := sols;
    ls : Link_to_Solution;

  begin
    while not Is_Null(tmp) loop
      ls := Head_Of(tmp);
      declare
        xt : Standard_Complex_Vectors.Vector(ls.v'first..ls.v'last+1);
      begin
        xt(ls.v'range) := ls.v;
        xt(xt'last) := Create(1.0);
        y := Eval(xtm,xt);
      end;
      Inverse_Row_Transformation(mf,y);
      Update_Swap_Column(y,s);
      Normalize_and_Reduce_to_Fit(qlocmap,y);
      ls.v := Map(qlocmap,y);
      Set_Head(tmp,ls);
      tmp := Tail_Of(tmp);
    end loop;
  end First_Swap_Coordinates;

  procedure First_Swap_Coordinates
              ( n,k,r,big_r,dc,s : in integer32;
                q,p,qr,qc,pr,pc : in Standard_Natural_Vectors.Vector;
                mf : in DoblDobl_Complex_Matrices.Matrix;
                xtm : in DoblDobl_Complex_Poly_Matrices.Matrix;
                sols : in out DoblDobl_Complex_Solutions.Solution_List ) is

    use DoblDobl_Complex_Numbers;
    use DoblDobl_Complex_Solutions;

    qlocmap : constant Standard_Natural_Matrices.Matrix(1..n,1..k)
            := Checker_Localization_Patterns.Column_Pattern(n,k,q,qr,qc);
    y : DoblDobl_Complex_Matrices.Matrix(1..n,1..k);
    tmp : Solution_List := sols;
    ls : Link_to_Solution;

  begin
    while not Is_Null(tmp) loop
      ls := Head_Of(tmp);
      declare
        xt : DoblDobl_Complex_Vectors.Vector(ls.v'first..ls.v'last+1);
      begin
        xt(ls.v'range) := ls.v;
        xt(xt'last) := Create(integer(1));
        y := Eval(xtm,xt);
      end;
      Inverse_Row_Transformation(mf,y);
      Update_Swap_Column(y,s);
      Normalize_and_Reduce_to_Fit(qlocmap,y);
      ls.v := Map(qlocmap,y);
      Set_Head(tmp,ls);
      tmp := Tail_Of(tmp);
    end loop;
  end First_Swap_Coordinates;

  procedure First_Swap_Coordinates
              ( n,k,r,big_r,dc,s : in integer32;
                q,p,qr,qc,pr,pc : in Standard_Natural_Vectors.Vector;
                mf : in QuadDobl_Complex_Matrices.Matrix;
                xtm : in QuadDobl_Complex_Poly_Matrices.Matrix;
                sols : in out QuadDobl_Complex_Solutions.Solution_List ) is

    use QuadDobl_Complex_Numbers;
    use QuadDobl_Complex_Solutions;

    qlocmap : constant Standard_Natural_Matrices.Matrix(1..n,1..k)
            := Checker_Localization_Patterns.Column_Pattern(n,k,q,qr,qc);
    y : QuadDobl_Complex_Matrices.Matrix(1..n,1..k);
    tmp : Solution_List := sols;
    ls : Link_to_Solution;

  begin
    while not Is_Null(tmp) loop
      ls := Head_Of(tmp);
      declare
        xt : QuadDobl_Complex_Vectors.Vector(ls.v'first..ls.v'last+1);
      begin
        xt(ls.v'range) := ls.v;
        xt(xt'last) := Create(integer(1));
        y := Eval(xtm,xt);
      end;
      Inverse_Row_Transformation(mf,y);
      Update_Swap_Column(y,s);
      Normalize_and_Reduce_to_Fit(qlocmap,y);
      ls.v := Map(qlocmap,y);
      Set_Head(tmp,ls);
      tmp := Tail_Of(tmp);
    end loop;
  end First_Swap_Coordinates;

  procedure First_Swap_Coordinates
              ( file : in file_type; n,k,r,big_r,dc,s : in integer32;
                q,p,qr,qc,pr,pc : in Standard_Natural_Vectors.Vector;
                mf : in Standard_Complex_Matrices.Matrix;
                xtm : in Standard_Complex_Poly_Matrices.Matrix;
                x : in out Standard_Complex_Vectors.Vector ) is

    use Standard_Complex_Numbers;

    plocmap : constant Standard_Natural_Matrices.Matrix(1..n,1..k)
            := Checker_Localization_Patterns.Column_Pattern(n,k,p,pr,pc);
    qlocmap : constant Standard_Natural_Matrices.Matrix(1..n,1..k)
            := Checker_Localization_Patterns.Column_Pattern(n,k,q,qr,qc);
    y : Standard_Complex_Matrices.Matrix(1..n,1..k);
    xt : Standard_Complex_Vectors.Vector(x'first..x'last+1);

  begin
    put(file,"Swap type I with critical row = "); put(file,r,1);
    put(file,", R = "); put(file,big_r,1);
    put(file," and s = "); put(file,s,1); put_line(file,".");
    put_line(file,"The start localization map : "); put(file,plocmap);
    put_line(file,"The target localization map : "); put(file,qlocmap);
   -- put_line(file,"The matrix xtm : "); put(file,xtm);
    xt(x'range) := x;
    xt(xt'last) := Create(1.0);
    put_line(file,"The vector xt : "); put_line(file,xt);
    y := Eval(xtm,xt);
   -- put_line(file,"The matrix xtm evaluated at the solution : ");
   -- put(file,y,2);
    Inverse_Row_Transformation(mf,y);
   -- put_line(file,"after the inverse transformation :"); put(file,y,3);
    Update_Swap_Column(y,s);
   -- put_line(file,"After updating the swap column :"); put(file,y,3);
    Normalize_and_Reduce_to_Fit(qlocmap,y);
   -- put_line(file,"The transformed plane :"); put(file,y,3);
    x := Map(qlocmap,y);
  end First_Swap_Coordinates;

  procedure First_Swap_Coordinates
              ( file : in file_type; n,k,r,big_r,dc,s : in integer32;
                q,p,qr,qc,pr,pc : in Standard_Natural_Vectors.Vector;
                mf : in DoblDobl_Complex_Matrices.Matrix;
                xtm : in DoblDobl_Complex_Poly_Matrices.Matrix;
                x : in out DoblDobl_Complex_Vectors.Vector ) is

    use DoblDobl_Complex_Numbers;

    plocmap : constant Standard_Natural_Matrices.Matrix(1..n,1..k)
            := Checker_Localization_Patterns.Column_Pattern(n,k,p,pr,pc);
    qlocmap : constant Standard_Natural_Matrices.Matrix(1..n,1..k)
            := Checker_Localization_Patterns.Column_Pattern(n,k,q,qr,qc);
    y : DoblDobl_Complex_Matrices.Matrix(1..n,1..k);
    xt : DoblDobl_Complex_Vectors.Vector(x'first..x'last+1);

  begin
    put(file,"Swap type I with critical row = "); put(file,r,1);
    put(file,", R = "); put(file,big_r,1);
    put(file," and s = "); put(file,s,1); put_line(file,".");
    put_line(file,"The start localization map : "); put(file,plocmap);
    put_line(file,"The target localization map : "); put(file,qlocmap);
   -- put_line(file,"The matrix xtm : "); put(file,xtm);
    xt(x'range) := x;
    xt(xt'last) := Create(integer(1));
    put_line(file,"The vector xt : "); put_line(file,xt);
    y := Eval(xtm,xt);
   -- put_line(file,"The matrix xtm evaluated at the solution : ");
   -- put(file,y,2);
    Inverse_Row_Transformation(mf,y);
   -- put_line(file,"after the inverse transformation :"); put(file,y,3);
    Update_Swap_Column(y,s);
    put_line(file,"After updating the swap column :"); put(file,y,3);
    Normalize_and_Reduce_to_Fit(qlocmap,y);
    put_line(file,"The transformed plane :"); put(file,y,3);
    x := Map(qlocmap,y);
  end First_Swap_Coordinates;

  procedure First_Swap_Coordinates
              ( file : in file_type; n,k,r,big_r,dc,s : in integer32;
                q,p,qr,qc,pr,pc : in Standard_Natural_Vectors.Vector;
                mf : in QuadDobl_Complex_Matrices.Matrix;
                xtm : in QuadDobl_Complex_Poly_Matrices.Matrix;
                x : in out QuadDobl_Complex_Vectors.Vector ) is

    use QuadDobl_Complex_Numbers;

    plocmap : constant Standard_Natural_Matrices.Matrix(1..n,1..k)
            := Checker_Localization_Patterns.Column_Pattern(n,k,p,pr,pc);
    qlocmap : constant Standard_Natural_Matrices.Matrix(1..n,1..k)
            := Checker_Localization_Patterns.Column_Pattern(n,k,q,qr,qc);
    y : QuadDobl_Complex_Matrices.Matrix(1..n,1..k);
    xt : QuadDobl_Complex_Vectors.Vector(x'first..x'last+1);

  begin
    put(file,"Swap type I with critical row = "); put(file,r,1);
    put(file,", R = "); put(file,big_r,1);
    put(file," and s = "); put(file,s,1); put_line(file,".");
    put_line(file,"The start localization map : "); put(file,plocmap);
    put_line(file,"The target localization map : "); put(file,qlocmap);
   -- put_line(file,"The matrix xtm : "); put(file,xtm);
    xt(x'range) := x;
    xt(xt'last) := Create(integer(1));
    put_line(file,"The vector xt : "); put_line(file,xt);
    y := Eval(xtm,xt);
   -- put_line(file,"The matrix xtm evaluated at the solution : ");
   -- put(file,y,2);
    Inverse_Row_Transformation(mf,y);
   -- put_line(file,"after the inverse transformation :"); put(file,y,3);
    Update_Swap_Column(y,s);
   -- put_line(file,"After updating the swap column :"); put(file,y,3);
    Normalize_and_Reduce_to_Fit(qlocmap,y);
   -- put_line(file,"The transformed plane :"); put(file,y,3);
    x := Map(qlocmap,y);
  end First_Swap_Coordinates;

  procedure First_Swap_Coordinates
              ( file : in file_type; n,k,r,big_r,dc,s : in integer32;
                q,p,qr,qc,pr,pc : in Standard_Natural_Vectors.Vector;
                mf : in Standard_Complex_Matrices.Matrix;
                xtm : in Standard_Complex_Poly_Matrices.Matrix;
                sols : in out Standard_Complex_Solutions.Solution_List ) is

    use Standard_Complex_Numbers;
    use Standard_Complex_Solutions;

    plocmap : constant Standard_Natural_Matrices.Matrix(1..n,1..k)
            := Checker_Localization_Patterns.Column_Pattern(n,k,p,pr,pc);
    qlocmap : constant Standard_Natural_Matrices.Matrix(1..n,1..k)
            := Checker_Localization_Patterns.Column_Pattern(n,k,q,qr,qc);
    y : Standard_Complex_Matrices.Matrix(1..n,1..k);
    tmp : Solution_List := sols;
    ls : Link_to_Solution;

  begin
    put(file,"Swap type I with critical row = "); put(file,r,1);
    put(file,", R = "); put(file,big_r,1);
    put(file," and s = "); put(file,s,1); put_line(file,".");
    put_line(file,"The start localization map : "); put(file,plocmap);
    put_line(file,"The target localization map : "); put(file,qlocmap);
   -- put_line(file,"The matrix xtm : "); put(file,xtm);
    while not Is_Null(tmp) loop
      ls := Head_Of(tmp);
      declare
        xt : Standard_Complex_Vectors.Vector(ls.v'first..ls.v'last+1);
      begin
        xt(ls.v'range) := ls.v;
        xt(xt'last) := Create(1.0);
        put_line(file,"The vector xt : "); put_line(file,xt);
        y := Eval(xtm,xt);
      end;
     -- put_line(file,"The matrix xtm evaluated at the solution : ");
     -- put(file,y,2);
      Inverse_Row_Transformation(mf,y);
     -- put_line(file,"after the inverse transformation :"); put(file,y,3);
      Update_Swap_Column(y,s);
     -- put_line(file,"After updating the swap column :"); put(file,y,3);
      Normalize_and_Reduce_to_Fit(qlocmap,y);
     -- put_line(file,"The transformed plane :"); put(file,y,3);
      ls.v := Map(qlocmap,y);
      Set_Head(tmp,ls);
      tmp := Tail_Of(tmp);
    end loop;
  end First_Swap_Coordinates;

  procedure First_Swap_Coordinates
              ( file : in file_type; n,k,r,big_r,dc,s : in integer32;
                q,p,qr,qc,pr,pc : in Standard_Natural_Vectors.Vector;
                mf : in DoblDobl_Complex_Matrices.Matrix;
                xtm : in DoblDobl_Complex_Poly_Matrices.Matrix;
                sols : in out DoblDobl_Complex_Solutions.Solution_List ) is

    use DoblDobl_Complex_Numbers;
    use DoblDobl_Complex_Solutions;

    plocmap : constant Standard_Natural_Matrices.Matrix(1..n,1..k)
            := Checker_Localization_Patterns.Column_Pattern(n,k,p,pr,pc);
    qlocmap : constant Standard_Natural_Matrices.Matrix(1..n,1..k)
            := Checker_Localization_Patterns.Column_Pattern(n,k,q,qr,qc);
    y : DoblDobl_Complex_Matrices.Matrix(1..n,1..k);
    tmp : Solution_List := sols;
    ls : Link_to_Solution;

  begin
    put(file,"Swap type I with critical row = "); put(file,r,1);
    put(file,", R = "); put(file,big_r,1);
    put(file," and s = "); put(file,s,1); put_line(file,".");
    put_line(file,"The start localization map : "); put(file,plocmap);
    put_line(file,"The target localization map : "); put(file,qlocmap);
   -- put_line(file,"The matrix xtm : "); put(file,xtm);
    while not Is_Null(tmp) loop
      ls := Head_Of(tmp);
      declare
        xt : DoblDobl_Complex_Vectors.Vector(ls.v'first..ls.v'last+1);
      begin
        xt(ls.v'range) := ls.v;
        xt(xt'last) := Create(integer(1));
        put_line(file,"The vector xt : "); put_line(file,xt);
        y := Eval(xtm,xt);
      end;
     -- put_line(file,"The matrix xtm evaluated at the solution : ");
     -- put(file,y,2);
      Inverse_Row_Transformation(mf,y);
     -- put_line(file,"after the inverse transformation :"); put(file,y,3);
      Update_Swap_Column(y,s);
     -- put_line(file,"After updating the swap column :"); put(file,y,3);
      Normalize_and_Reduce_to_Fit(qlocmap,y);
     -- put_line(file,"The transformed plane :"); put(file,y,3);
      ls.v := Map(qlocmap,y);
      Set_Head(tmp,ls);
      tmp := Tail_Of(tmp);
    end loop;
  end First_Swap_Coordinates;

  procedure First_Swap_Coordinates
              ( file : in file_type; n,k,r,big_r,dc,s : in integer32;
                q,p,qr,qc,pr,pc : in Standard_Natural_Vectors.Vector;
                mf : in QuadDobl_Complex_Matrices.Matrix;
                xtm : in QuadDobl_Complex_Poly_Matrices.Matrix;
                sols : in out QuadDobl_Complex_Solutions.Solution_List ) is

    use QuadDobl_Complex_Numbers;
    use QuadDobl_Complex_Solutions;

    plocmap : constant Standard_Natural_Matrices.Matrix(1..n,1..k)
            := Checker_Localization_Patterns.Column_Pattern(n,k,p,pr,pc);
    qlocmap : constant Standard_Natural_Matrices.Matrix(1..n,1..k)
            := Checker_Localization_Patterns.Column_Pattern(n,k,q,qr,qc);
    y : QuadDobl_Complex_Matrices.Matrix(1..n,1..k);
    tmp : Solution_List := sols;
    ls : Link_to_Solution;

  begin
    put(file,"Swap type I with critical row = "); put(file,r,1);
    put(file,", R = "); put(file,big_r,1);
    put(file," and s = "); put(file,s,1); put_line(file,".");
    put_line(file,"The start localization map : "); put(file,plocmap);
    put_line(file,"The target localization map : "); put(file,qlocmap);
   -- put_line(file,"The matrix xtm : "); put(file,xtm);
    while not Is_Null(tmp) loop
      ls := Head_Of(tmp);
      declare
        xt : QuadDobl_Complex_Vectors.Vector(ls.v'first..ls.v'last+1);
      begin
        xt(ls.v'range) := ls.v;
        xt(xt'last) := Create(integer(1));
        put_line(file,"The vector xt : "); put_line(file,xt);
        y := Eval(xtm,xt);
      end;
     -- put_line(file,"The matrix xtm evaluated at the solution : ");
     -- put(file,y,2);
      Inverse_Row_Transformation(mf,y);
     -- put_line(file,"after the inverse transformation :"); put(file,y,3);
      Update_Swap_Column(y,s);
     -- put_line(file,"After updating the swap column :"); put(file,y,3);
      Normalize_and_Reduce_to_Fit(qlocmap,y);
     -- put_line(file,"The transformed plane :"); put(file,y,3);
      ls.v := Map(qlocmap,y);
      Set_Head(tmp,ls);
      tmp := Tail_Of(tmp);
    end loop;
  end First_Swap_Coordinates;

  procedure Second_Swap_Coordinates
              ( n,k,r,s : in integer32;
                p,rows,cols : in Standard_Natural_Vectors.Vector;
                mf : in Standard_Complex_Matrices.Matrix;
                xtm : in Standard_Complex_Poly_Matrices.Matrix;
                x : in out Standard_Complex_Vectors.Vector ) is

    use Standard_Complex_Numbers;

    locmap : constant Standard_Natural_Matrices.Matrix(1..n,1..k)
           := Checker_Localization_Patterns.Column_Pattern(n,k,p,rows,cols);
    y : Standard_Complex_Matrices.Matrix(1..n,1..k); -- := Map(locmap,x);
    xt : Standard_Complex_Vectors.Vector(x'first..x'last+1);

  begin
    xt(x'range) := x;
    xt(xt'last) := Create(1.0);
    y := Eval(xtm,xt);
    Inverse_Row_Transformation(mf,y);
    Update_Swap_Column(y,s);
    Normalize_and_Reduce_to_Fit(locmap,y);
    x := Map(locmap,y);
  end Second_Swap_Coordinates;

  procedure Second_Swap_Coordinates
              ( n,k,r,s : in integer32;
                p,rows,cols : in Standard_Natural_Vectors.Vector;
                mf : in DoblDobl_Complex_Matrices.Matrix;
                xtm : in DoblDobl_Complex_Poly_Matrices.Matrix;
                x : in out DoblDobl_Complex_Vectors.Vector ) is

    use DoblDobl_Complex_Numbers;

    locmap : constant Standard_Natural_Matrices.Matrix(1..n,1..k)
           := Checker_Localization_Patterns.Column_Pattern(n,k,p,rows,cols);
    y : DoblDobl_Complex_Matrices.Matrix(1..n,1..k); -- := Map(locmap,x);
    xt : DoblDobl_Complex_Vectors.Vector(x'first..x'last+1);

  begin
    xt(x'range) := x;
    xt(xt'last) := Create(integer(1));
    y := Eval(xtm,xt);
    Inverse_Row_Transformation(mf,y);
    Update_Swap_Column(y,s);
    Normalize_and_Reduce_to_Fit(locmap,y);
    x := Map(locmap,y);
  end Second_Swap_Coordinates;

  procedure Second_Swap_Coordinates
              ( n,k,r,s : in integer32;
                p,rows,cols : in Standard_Natural_Vectors.Vector;
                mf : in QuadDobl_Complex_Matrices.Matrix;
                xtm : in QuadDobl_Complex_Poly_Matrices.Matrix;
                x : in out QuadDobl_Complex_Vectors.Vector ) is

    use QuadDobl_Complex_Numbers;

    locmap : constant Standard_Natural_Matrices.Matrix(1..n,1..k)
           := Checker_Localization_Patterns.Column_Pattern(n,k,p,rows,cols);
    y : QuadDobl_Complex_Matrices.Matrix(1..n,1..k); -- := Map(locmap,x);
    xt : QuadDobl_Complex_Vectors.Vector(x'first..x'last+1);

  begin
    xt(x'range) := x;
    xt(xt'last) := Create(integer(1));
    y := Eval(xtm,xt);
    Inverse_Row_Transformation(mf,y);
    Update_Swap_Column(y,s);
    Normalize_and_Reduce_to_Fit(locmap,y);
    x := Map(locmap,y);
  end Second_Swap_Coordinates;

  procedure Second_Swap_Coordinates
              ( n,k,r,s : in integer32;
                p,rows,cols : in Standard_Natural_Vectors.Vector;
                mf : in Standard_Complex_Matrices.Matrix;
                xtm : in Standard_Complex_Poly_Matrices.Matrix;
                sols : in out Standard_Complex_Solutions.Solution_List ) is

    use Standard_Complex_Numbers;
    use Standard_Complex_Solutions;

    locmap : constant Standard_Natural_Matrices.Matrix(1..n,1..k)
           := Checker_Localization_Patterns.Column_Pattern(n,k,p,rows,cols);
    y : Standard_Complex_Matrices.Matrix(1..n,1..k); -- := Map(locmap,x);
    tmp : Solution_List := sols;
    ls : Link_to_Solution;

  begin
    while not Is_Null(tmp) loop
      ls := Head_Of(tmp);
      declare
        xt : Standard_Complex_Vectors.Vector(ls.v'first..ls.v'last+1);
      begin
        xt(ls.v'range) := ls.v;
        xt(xt'last) := Create(1.0);
        y := Eval(xtm,xt);
      end;
      Inverse_Row_Transformation(mf,y);
      Update_Swap_Column(y,s);
      Normalize_and_Reduce_to_Fit(locmap,y);
      ls.v := Map(locmap,y);
      Set_Head(tmp,ls);
      tmp := Tail_Of(tmp);
    end loop;
  end Second_Swap_Coordinates;

  procedure Second_Swap_Coordinates
              ( n,k,r,s : in integer32;
                p,rows,cols : in Standard_Natural_Vectors.Vector;
                mf : in DoblDobl_Complex_Matrices.Matrix;
                xtm : in DoblDobl_Complex_Poly_Matrices.Matrix;
                sols : in out DoblDobl_Complex_Solutions.Solution_List ) is

    use DoblDobl_Complex_Numbers;
    use DoblDobl_Complex_Solutions;

    locmap : constant Standard_Natural_Matrices.Matrix(1..n,1..k)
           := Checker_Localization_Patterns.Column_Pattern(n,k,p,rows,cols);
    y : DoblDobl_Complex_Matrices.Matrix(1..n,1..k); -- := Map(locmap,x);
    tmp : Solution_List := sols;
    ls : Link_to_Solution;

  begin
    while not Is_Null(tmp) loop
      ls := Head_Of(tmp);
      declare
        xt : DoblDobl_Complex_Vectors.Vector(ls.v'first..ls.v'last+1);
      begin
        xt(ls.v'range) := ls.v;
        xt(xt'last) := Create(integer(1));
        y := Eval(xtm,xt);
      end;
      Inverse_Row_Transformation(mf,y);
      Update_Swap_Column(y,s);
      Normalize_and_Reduce_to_Fit(locmap,y);
      ls.v := Map(locmap,y);
      Set_Head(tmp,ls);
      tmp := Tail_Of(tmp);
    end loop;
  end Second_Swap_Coordinates;

  procedure Second_Swap_Coordinates
              ( n,k,r,s : in integer32;
                p,rows,cols : in Standard_Natural_Vectors.Vector;
                mf : in QuadDobl_Complex_Matrices.Matrix;
                xtm : in QuadDobl_Complex_Poly_Matrices.Matrix;
                sols : in out QuadDobl_Complex_Solutions.Solution_List ) is

    use QuadDobl_Complex_Numbers;
    use QuadDobl_Complex_Solutions;

    locmap : constant Standard_Natural_Matrices.Matrix(1..n,1..k)
           := Checker_Localization_Patterns.Column_Pattern(n,k,p,rows,cols);
    y : QuadDobl_Complex_Matrices.Matrix(1..n,1..k); -- := Map(locmap,x);
    tmp : Solution_List := sols;
    ls : Link_to_Solution;

  begin
    while not Is_Null(tmp) loop
      ls := Head_Of(tmp);
      declare
        xt : QuadDobl_Complex_Vectors.Vector(ls.v'first..ls.v'last+1);
      begin
        xt(ls.v'range) := ls.v;
        xt(xt'last) := Create(integer(1));
        y := Eval(xtm,xt);
      end;
      Inverse_Row_Transformation(mf,y);
      Update_Swap_Column(y,s);
      Normalize_and_Reduce_to_Fit(locmap,y);
      ls.v := Map(locmap,y);
      Set_Head(tmp,ls);
      tmp := Tail_Of(tmp);
    end loop;
  end Second_Swap_Coordinates;

  procedure Second_Swap_Coordinates
              ( file : in file_type; n,k,r,s : in integer32;
                p,rows,cols : in Standard_Natural_Vectors.Vector;
                mf : in Standard_Complex_Matrices.Matrix;
                xtm : in Standard_Complex_Poly_Matrices.Matrix;
                x : in out Standard_Complex_Vectors.Vector ) is

    use Standard_Complex_Numbers;

    locmap : constant Standard_Natural_Matrices.Matrix(1..n,1..k)
           := Checker_Localization_Patterns.Column_Pattern(n,k,p,rows,cols);
    y : Standard_Complex_Matrices.Matrix(1..n,1..k); -- := Map(locmap,x);
    xt : Standard_Complex_Vectors.Vector(x'first..x'last+1);

  begin
    put(file,"Swap type II with critical row = "); put(file,r,1);
    put(file," and s = "); put(file,s,1); put_line(file,".");
   -- put_line(file,"The localization map : "); put(file,locmap);
   -- put_line(file,"The matrix xtm : "); put(file,xtm);
    xt(x'range) := x;
    xt(xt'last) := Create(1.0);
    put_line(file,"The vector xt : "); put_line(file,xt);
    y := Eval(xtm,xt);
   -- put_line(file,"The matrix xtm evaluated at the solution : ");
   -- put(file,y,2);
    Inverse_Row_Transformation(mf,y);
   -- put_line(file,"after the inverse transformation :"); put(file,y,3);
    Update_Swap_Column(y,s);
   -- put_line(file,"After updating the swap column :"); put(file,y,3);
    Normalize_and_Reduce_to_Fit(locmap,y);
   -- put_line(file,"The transformed plane :"); put(file,y,3);
    x := Map(locmap,y);
  end Second_Swap_Coordinates;

  procedure Second_Swap_Coordinates
              ( file : in file_type; n,k,r,s : in integer32;
                p,rows,cols : in Standard_Natural_Vectors.Vector;
                mf : in DoblDobl_Complex_Matrices.Matrix;
                xtm : in DoblDobl_Complex_Poly_Matrices.Matrix;
                x : in out DoblDobl_Complex_Vectors.Vector ) is

    use DoblDobl_Complex_Numbers;

    locmap : constant Standard_Natural_Matrices.Matrix(1..n,1..k)
           := Checker_Localization_Patterns.Column_Pattern(n,k,p,rows,cols);
    y : DoblDobl_Complex_Matrices.Matrix(1..n,1..k); -- := Map(locmap,x);
    xt : DoblDobl_Complex_Vectors.Vector(x'first..x'last+1);

  begin
    put(file,"Swap type II with critical row = "); put(file,r,1);
    put(file," and s = "); put(file,s,1); put_line(file,".");
   -- put_line(file,"The localization map : "); put(file,locmap);
   -- put_line(file,"The matrix xtm : "); put(file,xtm);
    xt(x'range) := x;
    xt(xt'last) := Create(integer(1));
    put_line(file,"The vector xt : "); put_line(file,xt);
    y := Eval(xtm,xt);
   -- put_line(file,"The matrix xtm evaluated at the solution : ");
   -- put(file,y,2);
    Inverse_Row_Transformation(mf,y);
   -- put_line(file,"after the inverse transformation :"); put(file,y,3);
    Update_Swap_Column(y,s);
   -- put_line(file,"After updating the swap column :"); put(file,y,3);
    Normalize_and_Reduce_to_Fit(locmap,y);
   -- put_line(file,"The transformed plane :"); put(file,y,3);
    x := Map(locmap,y);
  end Second_Swap_Coordinates;

  procedure Second_Swap_Coordinates
              ( file : in file_type; n,k,r,s : in integer32;
                p,rows,cols : in Standard_Natural_Vectors.Vector;
                mf : in QuadDobl_Complex_Matrices.Matrix;
                xtm : in QuadDobl_Complex_Poly_Matrices.Matrix;
                x : in out QuadDobl_Complex_Vectors.Vector ) is

    use QuadDobl_Complex_Numbers;

    locmap : constant Standard_Natural_Matrices.Matrix(1..n,1..k)
           := Checker_Localization_Patterns.Column_Pattern(n,k,p,rows,cols);
    y : QuadDobl_Complex_Matrices.Matrix(1..n,1..k); -- := Map(locmap,x);
    xt : QuadDobl_Complex_Vectors.Vector(x'first..x'last+1);

  begin
    put(file,"Swap type II with critical row = "); put(file,r,1);
    put(file," and s = "); put(file,s,1); put_line(file,".");
   -- put_line(file,"The localization map : "); put(file,locmap);
   -- put_line(file,"The matrix xtm : "); put(file,xtm);
    xt(x'range) := x;
    xt(xt'last) := Create(integer(1));
    put_line(file,"The vector xt : "); put_line(file,xt);
    y := Eval(xtm,xt);
   -- put_line(file,"The matrix xtm evaluated at the solution : ");
   -- put(file,y,2);
    Inverse_Row_Transformation(mf,y);
   -- put_line(file,"after the inverse transformation :"); put(file,y,3);
    Update_Swap_Column(y,s);
   -- put_line(file,"After updating the swap column :"); put(file,y,3);
    Normalize_and_Reduce_to_Fit(locmap,y);
   -- put_line(file,"The transformed plane :"); put(file,y,3);
    x := Map(locmap,y);
  end Second_Swap_Coordinates;

  procedure Second_Swap_Coordinates
              ( file : in file_type; n,k,r,s : in integer32;
                p,rows,cols : in Standard_Natural_Vectors.Vector;
                mf : in Standard_Complex_Matrices.Matrix;
                xtm : in Standard_Complex_Poly_Matrices.Matrix;
                sols : in out Standard_Complex_Solutions.Solution_List ) is

    use Standard_Complex_Numbers;
    use Standard_Complex_Solutions;

    locmap : constant Standard_Natural_Matrices.Matrix(1..n,1..k)
           := Checker_Localization_Patterns.Column_Pattern(n,k,p,rows,cols);
    y : Standard_Complex_Matrices.Matrix(1..n,1..k); -- := Map(locmap,x);
    tmp : Solution_List := sols;
    ls : Link_to_Solution;

  begin
    put(file,"Swap type II with critical row = "); put(file,r,1);
    put(file," and s = "); put(file,s,1); put_line(file,".");
   -- put_line(file,"The localization map : "); put(file,locmap);
   -- put_line(file,"The matrix xtm : "); put(file,xtm);
    while not Is_Null(tmp) loop
      ls := Head_Of(tmp);
      declare
        xt : Standard_Complex_Vectors.Vector(ls.v'first..ls.v'last+1);
      begin
        xt(ls.v'range) := ls.v;
        xt(xt'last) := Create(1.0);
        put_line(file,"The vector xt : "); put_line(file,xt);
        y := Eval(xtm,xt);
      end;
     -- put_line(file,"The matrix xtm evaluated at the solution : ");
     -- put(file,y,2);
      Inverse_Row_Transformation(mf,y);
     -- put_line(file,"after the inverse transformation :"); put(file,y,3);
      Update_Swap_Column(y,s);
     -- put_line(file,"After updating the swap column :"); put(file,y,3);
      Normalize_and_Reduce_to_Fit(locmap,y);
     -- put_line(file,"The transformed plane :"); put(file,y,3);
      ls.v := Map(locmap,y);
      Set_Head(tmp,ls);
      tmp := Tail_Of(tmp);
    end loop;
  end Second_Swap_Coordinates;

  procedure Second_Swap_Coordinates
              ( file : in file_type; n,k,r,s : in integer32;
                p,rows,cols : in Standard_Natural_Vectors.Vector;
                mf : in DoblDobl_Complex_Matrices.Matrix;
                xtm : in DoblDobl_Complex_Poly_Matrices.Matrix;
                sols : in out DoblDobl_Complex_Solutions.Solution_List ) is

    use DoblDobl_Complex_Numbers;
    use DoblDobl_Complex_Solutions;

    locmap : constant Standard_Natural_Matrices.Matrix(1..n,1..k)
           := Checker_Localization_Patterns.Column_Pattern(n,k,p,rows,cols);
    y : DoblDobl_Complex_Matrices.Matrix(1..n,1..k); -- := Map(locmap,x);
    tmp : Solution_List := sols;
    ls : Link_to_Solution;

  begin
    put(file,"Swap type II with critical row = "); put(file,r,1);
    put(file," and s = "); put(file,s,1); put_line(file,".");
   -- put_line(file,"The localization map : "); put(file,locmap);
   -- put_line(file,"The matrix xtm : "); put(file,xtm);
    while not Is_Null(tmp) loop
      ls := Head_Of(tmp);
      declare
        xt : DoblDobl_Complex_Vectors.Vector(ls.v'first..ls.v'last+1);
      begin
        xt(ls.v'range) := ls.v;
        xt(xt'last) := Create(integer(1));
        put_line(file,"The vector xt : "); put_line(file,xt);
        y := Eval(xtm,xt);
      end;
     -- put_line(file,"The matrix xtm evaluated at the solution : ");
     -- put(file,y,2);
      Inverse_Row_Transformation(mf,y);
     -- put_line(file,"after the inverse transformation :"); put(file,y,3);
      Update_Swap_Column(y,s);
     -- put_line(file,"After updating the swap column :"); put(file,y,3);
      Normalize_and_Reduce_to_Fit(locmap,y);
     -- put_line(file,"The transformed plane :"); put(file,y,3);
      ls.v := Map(locmap,y);
      Set_Head(tmp,ls);
      tmp := Tail_Of(tmp);
    end loop;
  end Second_Swap_Coordinates;

  procedure Second_Swap_Coordinates
              ( file : in file_type; n,k,r,s : in integer32;
                p,rows,cols : in Standard_Natural_Vectors.Vector;
                mf : in QuadDobl_Complex_Matrices.Matrix;
                xtm : in QuadDobl_Complex_Poly_Matrices.Matrix;
                sols : in out QuadDobl_Complex_Solutions.Solution_List ) is

    use QuadDobl_Complex_Numbers;
    use QuadDobl_Complex_Solutions;

    locmap : constant Standard_Natural_Matrices.Matrix(1..n,1..k)
           := Checker_Localization_Patterns.Column_Pattern(n,k,p,rows,cols);
    y : QuadDobl_Complex_Matrices.Matrix(1..n,1..k); -- := Map(locmap,x);
    tmp : Solution_List := sols;
    ls : Link_to_Solution;

  begin
    put(file,"Swap type II with critical row = "); put(file,r,1);
    put(file," and s = "); put(file,s,1); put_line(file,".");
   -- put_line(file,"The localization map : "); put(file,locmap);
   -- put_line(file,"The matrix xtm : "); put(file,xtm);
    while not Is_Null(tmp) loop
      ls := Head_Of(tmp);
      declare
        xt : QuadDobl_Complex_Vectors.Vector(ls.v'first..ls.v'last+1);
      begin
        xt(ls.v'range) := ls.v;
        xt(xt'last) := Create(integer(1));
        put_line(file,"The vector xt : "); put_line(file,xt);
        y := Eval(xtm,xt);
      end;
     -- put_line(file,"The matrix xtm evaluated at the solution : ");
     -- put(file,y,2);
      Inverse_Row_Transformation(mf,y);
     -- put_line(file,"after the inverse transformation :"); put(file,y,3);
      Update_Swap_Column(y,s);
     -- put_line(file,"After updating the swap column :"); put(file,y,3);
      Normalize_and_Reduce_to_Fit(locmap,y);
     -- put_line(file,"The transformed plane :"); put(file,y,3);
      ls.v := Map(locmap,y);
      Set_Head(tmp,ls);
      tmp := Tail_Of(tmp);
    end loop;
  end Second_Swap_Coordinates;

-- PART III : coordinate definitions for the stay and swap homotopies

  procedure Initialize_Moving_Plane
              ( x : in out Standard_Complex_Poly_Matrices.Matrix;
                m : in Standard_Natural_Matrices.Matrix ) is

    use Standard_Complex_Numbers;
    use Standard_Complex_Polynomials;

    n : constant integer32 := integer32(Degree_of_Freedom(m));
    ind : integer32 := 0;
    t : Term;

  begin
    t.dg := new Standard_Natural_Vectors.Vector'(1..n+1 => 0);
    t.cf := Create(1.0);
    for i in x'range(1) loop
      for j in x'range(2) loop
        if m(i,j) = 2           -- need to count order of variables
         then ind := ind + 1;   -- no matter whether used or not
        end if;
        if m(i,j) = 0 then
          x(i,j) := Null_Poly;
        elsif m(i,j) = 1 then
          x(i,j) := Create(t);
        else
          t.dg(ind) := 1;         
          x(i,j) := Create(t);
          t.dg(ind) := 0;
        end if;
      end loop;
    end loop;
  end Initialize_Moving_Plane;

  procedure Initialize_Moving_Plane
              ( x : in out DoblDobl_Complex_Poly_Matrices.Matrix;
                m : in Standard_Natural_Matrices.Matrix ) is

    use DoblDobl_Complex_Numbers;
    use DoblDobl_Complex_Polynomials;

    n : constant integer32 := integer32(Degree_of_Freedom(m));
    ind : integer32 := 0;
    t : Term;

  begin
    t.dg := new Standard_Natural_Vectors.Vector'(1..n+1 => 0);
    t.cf := Create(integer(1));
    for i in x'range(1) loop
      for j in x'range(2) loop
        if m(i,j) = 2           -- need to count order of variables
         then ind := ind + 1;   -- no matter whether used or not
        end if;
        if m(i,j) = 0 then
          x(i,j) := Null_Poly;
        elsif m(i,j) = 1 then
          x(i,j) := Create(t);
        else
          t.dg(ind) := 1;         
          x(i,j) := Create(t);
          t.dg(ind) := 0;
        end if;
      end loop;
    end loop;
  end Initialize_Moving_Plane;

  procedure Initialize_Moving_Plane
              ( x : in out QuadDobl_Complex_Poly_Matrices.Matrix;
                m : in Standard_Natural_Matrices.Matrix ) is

    use QuadDobl_Complex_Numbers;
    use QuadDobl_Complex_Polynomials;

    n : constant integer32 := integer32(Degree_of_Freedom(m));
    ind : integer32 := 0;
    t : Term;

  begin
    t.dg := new Standard_Natural_Vectors.Vector'(1..n+1 => 0);
    t.cf := Create(integer(1));
    for i in x'range(1) loop
      for j in x'range(2) loop
        if m(i,j) = 2           -- need to count order of variables
         then ind := ind + 1;   -- no matter whether used or not
        end if;
        if m(i,j) = 0 then
          x(i,j) := Null_Poly;
        elsif m(i,j) = 1 then
          x(i,j) := Create(t);
        else
          t.dg(ind) := 1;         
          x(i,j) := Create(t);
          t.dg(ind) := 0;
        end if;
      end loop;
    end loop;
  end Initialize_Moving_Plane;

  procedure Initialize_Moving_Plane
              ( x : in out Standard_Complex_Poly_Matrices.Matrix;
                m : in Standard_Natural_Matrices.Matrix; s : in integer32 ) is

    use Standard_Complex_Numbers;
    use Standard_Complex_Polynomials;

    n : constant integer32 := integer32(Degree_of_Freedom(m));
    ind : integer32 := 0;
    t : Term;

  begin
    t.dg := new Standard_Natural_Vectors.Vector'(1..n+1 => 0);
    t.cf := Create(1.0);
    for i in x'range(1) loop
      for j in x'range(2) loop
        if m(i,j) = 2           -- need to count order of variables
         then ind := ind + 1;   -- no matter whether used or not
        end if;
        if j = s or j = s+1 or m(i,j) = 0 then
          x(i,j) := Null_Poly;
        elsif m(i,j) = 1 then
          x(i,j) := Create(t);
        else
          t.dg(ind) := 1;         
          x(i,j) := Create(t);
          t.dg(ind) := 0;
        end if;
      end loop;
    end loop;
  end Initialize_Moving_Plane;

  procedure Initialize_Moving_Plane
              ( x : in out DoblDobl_Complex_Poly_Matrices.Matrix;
                m : in Standard_Natural_Matrices.Matrix; s : in integer32 ) is

    use DoblDobl_Complex_Numbers;
    use DoblDobl_Complex_Polynomials;

    n : constant integer32 := integer32(Degree_of_Freedom(m));
    ind : integer32 := 0;
    t : Term;

  begin
    t.dg := new Standard_Natural_Vectors.Vector'(1..n+1 => 0);
    t.cf := Create(integer(1));
    for i in x'range(1) loop
      for j in x'range(2) loop
        if m(i,j) = 2           -- need to count order of variables
         then ind := ind + 1;   -- no matter whether used or not
        end if;
        if j = s or j = s+1 or m(i,j) = 0 then
          x(i,j) := Null_Poly;
        elsif m(i,j) = 1 then
          x(i,j) := Create(t);
        else
          t.dg(ind) := 1;         
          x(i,j) := Create(t);
          t.dg(ind) := 0;
        end if;
      end loop;
    end loop;
  end Initialize_Moving_Plane;

  procedure Initialize_Moving_Plane
              ( x : in out QuadDobl_Complex_Poly_Matrices.Matrix;
                m : in Standard_Natural_Matrices.Matrix; s : in integer32 ) is

    use QuadDobl_Complex_Numbers;
    use QuadDobl_Complex_Polynomials;

    n : constant integer32 := integer32(Degree_of_Freedom(m));
    ind : integer32 := 0;
    t : Term;

  begin
    t.dg := new Standard_Natural_Vectors.Vector'(1..n+1 => 0);
    t.cf := Create(integer(1));
    for i in x'range(1) loop
      for j in x'range(2) loop
        if m(i,j) = 2           -- need to count order of variables
         then ind := ind + 1;   -- no matter whether used or not
        end if;
        if j = s or j = s+1 or m(i,j) = 0 then
          x(i,j) := Null_Poly;
        elsif m(i,j) = 1 then
          x(i,j) := Create(t);
        else
          t.dg(ind) := 1;         
          x(i,j) := Create(t);
          t.dg(ind) := 0;
        end if;
      end loop;
    end loop;
  end Initialize_Moving_Plane;

  function Is_Zone_A_Empty
             ( locmap : Standard_Natural_Matrices.Matrix;
               p : Standard_Natural_Vectors.Vector;
               r,s,dc : integer32 ) return boolean is

    res : boolean := true;

  begin
    for i in p'range loop
      if integer32(p(i)) < r then
       -- if p'last+1-i > p'last-dc+1
        if p'last+1-i > p'last-dc+1 and locmap(integer32(p(i)),s+1) = 2
         then res := false; exit;
        end if;
      end if;
    end loop;
    return res;
  end Is_Zone_A_Empty;

  procedure First_Swap_Plane
              ( x : in out Standard_Complex_Poly_Matrices.Matrix;
                r,big_r,dc,s : in integer32;
                p : in Standard_Natural_Vectors.Vector;
                locmap : in Standard_Natural_Matrices.Matrix ) is

    use Standard_Complex_Numbers;
    use Standard_Complex_Polynomials;

    dim : constant natural32 := Degree_of_Freedom(locmap);
    np1 : constant integer32 := integer32(dim)+1;
    ind : integer32;
    piv : constant integer32
        := Checker_Localization_Patterns.Rank(locmap,r+1,s+1);
    t : Term;
    empty_zone_A : constant boolean := Is_Zone_A_Empty(locmap,p,r,s,dc);

  begin
    t.dg := new Standard_Natural_Vectors.Vector'(1..np1 => 0);
    t.cf := Create(1.0);
    if not empty_zone_A
     then t.dg(piv) := 1;   -- multiply by y(r+1,s+1)
    end if;
    x(r,s) := Create(t);    -- m(r) of 1st swapped column s
    t.dg(np1) := 1; t.cf := Create(-1.0);
    x(r+1,s) := Create(t);  -- t*m(r+1) of 1st swapped column s
    t.dg(np1) := 0; t.dg(piv) := 0; t.cf := Create(1.0);
    for i in p'range loop
      if integer32(p(i)) < r then      -- in zones A and B
        if p'last+1-i > p'last-dc+1 then    -- in zone A
          if locmap(integer32(p(i)),s+1) = 2 then
            ind := Checker_Localization_Patterns.Rank
                     (locmap,integer32(p(i)),s+1);
            t.dg(ind) := 1; t.dg(np1) := 1; t.cf := Create(-1.0);
            x(integer32(p(i)),s) := Create(t);
            t.dg(ind) := 0; t.dg(np1) := 0; t.cf := Create(1.0);
          end if;
        elsif p'last+1-i < p'last-r+1 then -- in zone B
          if locmap(integer32(p(i)),s) = 2 then
            ind := Checker_Localization_Patterns.Rank
                     (locmap,integer32(p(i)),s);
            t.dg(ind) := 1;
            if not empty_zone_A
             then t.dg(piv) := 1;
            end if;
            x(integer32(p(i)),s) := Create(t);
            t.dg(ind) := 0; t.dg(piv) := 0;
          end if;
        elsif locmap(integer32(p(i)),s) = 2 then -- do not forget variables!
          ind := Checker_Localization_Patterns.Rank
                   (locmap,integer32(p(i)),s);
          if ind in t.dg'range then
            t.dg(ind) := 1;
            if not empty_zone_A
             then t.dg(piv) := 1;
            end if;
            x(integer32(p(i)),s) := Create(t);
            t.dg(ind) := 0; t.dg(piv) := 0;
          end if;
        end if;
      end if;
    end loop;
    t.dg(piv) := 1;
    x(r+1,s+1) := Create(t);    -- x(r+1,s+1)*m(r+1) for column s+1
    t.dg(piv) := 0;
    x(big_r,s+1) := Create(t);  -- m(R) in column s+1 
    for i in p'range loop
      if integer32(p(i)) /= r and integer32(p(i)) /= r+1
                              and integer32(p(i)) /= big_r then
       -- if locmap(p(i),s+1) = 2 then
          ind := Checker_Localization_Patterns.Rank
                   (locmap,integer32(p(i)),s+1);
        if ind in t.dg'range then
          t.dg(ind) := 1;
          x(integer32(p(i)),s+1) := Create(t);
          t.dg(ind) := 0;
        end if;
       -- end if;
      end if;
    end loop;
    Clear(t);
  end First_Swap_Plane;

  procedure First_Swap_Plane
              ( x : in out DoblDobl_Complex_Poly_Matrices.Matrix;
                r,big_r,dc,s : in integer32;
                p : in Standard_Natural_Vectors.Vector;
                locmap : in Standard_Natural_Matrices.Matrix ) is

    use DoblDobl_Complex_Numbers;
    use DoblDobl_Complex_Polynomials;

    dim : constant natural32 := Degree_of_Freedom(locmap);
    np1 : constant integer32 := integer32(dim)+1;
    ind : integer32;
    piv : constant integer32
        := Checker_Localization_Patterns.Rank(locmap,r+1,s+1);
    t : Term;
    empty_zone_A : constant boolean := Is_Zone_A_Empty(locmap,p,r,s,dc);

  begin
    t.dg := new Standard_Natural_Vectors.Vector'(1..np1 => 0);
    t.cf := Create(integer(1));
    if not empty_zone_A
     then t.dg(piv) := 1;   -- multiply by y(r+1,s+1)
    end if;
    x(r,s) := Create(t);    -- m(r) of 1st swapped column s
    t.dg(np1) := 1; t.cf := Create(integer(-1));
    x(r+1,s) := Create(t);  -- t*m(r+1) of 1st swapped column s
    t.dg(np1) := 0; t.dg(piv) := 0; t.cf := Create(integer(1));
    for i in p'range loop
      if integer32(p(i)) < r then      -- in zones A and B
        if p'last+1-i > p'last-dc+1 then    -- in zone A
          if locmap(integer32(p(i)),s+1) = 2 then
            ind := Checker_Localization_Patterns.Rank
                     (locmap,integer32(p(i)),s+1);
            t.dg(ind) := 1; t.dg(np1) := 1; t.cf := Create(integer(-1));
            x(integer32(p(i)),s) := Create(t);
            t.dg(ind) := 0; t.dg(np1) := 0; t.cf := Create(integer(1));
          end if;
        elsif p'last+1-i < p'last-r+1 then -- in zone B
          if locmap(integer32(p(i)),s) = 2 then
            ind := Checker_Localization_Patterns.Rank
                     (locmap,integer32(p(i)),s);
            t.dg(ind) := 1;
            if not empty_zone_A
             then t.dg(piv) := 1;
            end if;
            x(integer32(p(i)),s) := Create(t);
            t.dg(ind) := 0; t.dg(piv) := 0;
          end if;
        elsif locmap(integer32(p(i)),s) = 2 then -- do not forget variables!
          ind := Checker_Localization_Patterns.Rank
                   (locmap,integer32(p(i)),s);
          if ind in t.dg'range then
            t.dg(ind) := 1;
            if not empty_zone_A
             then t.dg(piv) := 1;
            end if;
            x(integer32(p(i)),s) := Create(t);
            t.dg(ind) := 0; t.dg(piv) := 0;
          end if;
        end if;
      end if;
    end loop;
    t.dg(piv) := 1;
    x(r+1,s+1) := Create(t);    -- x(r+1,s+1)*m(r+1) for column s+1
    t.dg(piv) := 0;
    x(big_r,s+1) := Create(t);  -- m(R) in column s+1 
    for i in p'range loop
      if integer32(p(i)) /= r and integer32(p(i)) /= r+1
                              and integer32(p(i)) /= big_r then
       -- if locmap(p(i),s+1) = 2 then
          ind := Checker_Localization_Patterns.Rank
                   (locmap,integer32(p(i)),s+1);
        if ind in t.dg'range then
          t.dg(ind) := 1;
          x(integer32(p(i)),s+1) := Create(t);
          t.dg(ind) := 0;
        end if;
       -- end if;
      end if;
    end loop;
    Clear(t);
  end First_Swap_Plane;

  procedure First_Swap_Plane
              ( x : in out QuadDobl_Complex_Poly_Matrices.Matrix;
                r,big_r,dc,s : in integer32;
                p : in Standard_Natural_Vectors.Vector;
                locmap : in Standard_Natural_Matrices.Matrix ) is

    use QuadDobl_Complex_Numbers;
    use QuadDobl_Complex_Polynomials;

    dim : constant natural32 := Degree_of_Freedom(locmap);
    np1 : constant integer32 := integer32(dim)+1;
    ind : integer32;
    piv : constant integer32
        := Checker_Localization_Patterns.Rank(locmap,r+1,s+1);
    t : Term;
    empty_zone_A : constant boolean := Is_Zone_A_Empty(locmap,p,r,s,dc);

  begin
    t.dg := new Standard_Natural_Vectors.Vector'(1..np1 => 0);
    t.cf := Create(integer(1));
    if not empty_zone_A
     then t.dg(piv) := 1;   -- multiply by y(r+1,s+1)
    end if;
    x(r,s) := Create(t);    -- m(r) of 1st swapped column s
    t.dg(np1) := 1; t.cf := Create(integer(-1));
    x(r+1,s) := Create(t);  -- t*m(r+1) of 1st swapped column s
    t.dg(np1) := 0; t.dg(piv) := 0; t.cf := Create(integer(1));
    for i in p'range loop
      if integer32(p(i)) < r then      -- in zones A and B
        if p'last+1-i > p'last-dc+1 then    -- in zone A
          if locmap(integer32(p(i)),s+1) = 2 then
            ind := Checker_Localization_Patterns.Rank
                     (locmap,integer32(p(i)),s+1);
            t.dg(ind) := 1; t.dg(np1) := 1; t.cf := Create(integer(-1));
            x(integer32(p(i)),s) := Create(t);
            t.dg(ind) := 0; t.dg(np1) := 0; t.cf := Create(integer(1));
          end if;
        elsif p'last+1-i < p'last-r+1 then -- in zone B
          if locmap(integer32(p(i)),s) = 2 then
            ind := Checker_Localization_Patterns.Rank
                     (locmap,integer32(p(i)),s);
            t.dg(ind) := 1;
            if not empty_zone_A
             then t.dg(piv) := 1;
            end if;
            x(integer32(p(i)),s) := Create(t);
            t.dg(ind) := 0; t.dg(piv) := 0;
          end if;
        elsif locmap(integer32(p(i)),s) = 2 then -- do not forget variables!
          ind := Checker_Localization_Patterns.Rank
                   (locmap,integer32(p(i)),s);
          if ind in t.dg'range then
            t.dg(ind) := 1;
            if not empty_zone_A
             then t.dg(piv) := 1;
            end if;
            x(integer32(p(i)),s) := Create(t);
            t.dg(ind) := 0; t.dg(piv) := 0;
          end if;
        end if;
      end if;
    end loop;
    t.dg(piv) := 1;
    x(r+1,s+1) := Create(t);    -- x(r+1,s+1)*m(r+1) for column s+1
    t.dg(piv) := 0;
    x(big_r,s+1) := Create(t);  -- m(R) in column s+1 
    for i in p'range loop
      if integer32(p(i)) /= r and integer32(p(i)) /= r+1
                              and integer32(p(i)) /= big_r then
       -- if locmap(p(i),s+1) = 2 then
          ind := Checker_Localization_Patterns.Rank
                   (locmap,integer32(p(i)),s+1);
        if ind in t.dg'range then
          t.dg(ind) := 1;
          x(integer32(p(i)),s+1) := Create(t);
          t.dg(ind) := 0;
        end if;
       -- end if;
      end if;
    end loop;
    Clear(t);
  end First_Swap_Plane;

  procedure First_Swap_Plane
              ( file : in file_type;
                x : in out Standard_Complex_Poly_Matrices.Matrix;
                r,big_r,dc,s : in integer32;
                p : in Standard_Natural_Vectors.Vector;
                locmap : in Standard_Natural_Matrices.Matrix ) is

    use Standard_Complex_Numbers;
    use Standard_Complex_Polynomials;

    dim : constant natural32 := Degree_of_Freedom(locmap);
    np1 : constant integer32 := integer32(dim)+1;
    ind : integer32;
    piv : constant integer32
        := Checker_Localization_Patterns.Rank(locmap,r+1,s+1);
    t : Term;
    empty_zone_A : constant boolean := Is_Zone_A_Empty(locmap,p,r,s,dc);
    failed : boolean := false;

  begin
   -- put_line(file,"the localization map : "); put(file,locmap);
    put(file,"dim = "); put(file,dim,1); new_line(file);
    put(file,"np1 = "); put(file,np1,1); new_line(file);
    t.dg := new Standard_Natural_Vectors.Vector'(1..np1 => 0);
    t.cf := Create(1.0);
    if not empty_zone_A then
      put_line(file,"empty_zone_A is false");
      t.dg(piv) := 1;   -- multiply by y(r+1,s+1)
      put(file,"multiply by x(");
      put(file,r+1,1); put(file,",");
      put(file,s+1,1); put(file,")");
      put(file," piv = "); put(file,piv,1); new_line(file);
    else
      put_line(file,"empty_zone_A is true");
    end if;
    x(r,s) := Create(t);    -- m(r) of 1st swapped column s
    t.dg(np1) := 1; t.cf := Create(-1.0);
    x(r+1,s) := Create(t);  -- t*m(r+1) of 1st swapped column s
    t.dg(np1) := 0; t.dg(piv) := 0; t.cf := Create(1.0);
    for i in p'range loop
      if integer32(p(i)) < r then      -- in zones A and B
        if p'last+1-i > p'last-dc+1 then    -- in zone A
          if locmap(integer32(p(i)),s+1) = 2 then -- why was this s+1?
         -- if locmap(integer32(p(i)),s) = 2 then -- it must be s+1!
            ind := Checker_Localization_Patterns.Rank -- because we take
                     (locmap,integer32(p(i)),s+1);    -- the variable at s+1
            t.dg(ind) := 1; t.dg(np1) := 1; t.cf := Create(-1.0);
            x(integer32(p(i)),s) := Create(t); -- is this a bug ???
           -- x(integer32(p(i)),s+1) := Create(t); -- assign to column s?
            t.dg(ind) := 0; t.dg(np1) := 0; t.cf := Create(1.0);
          end if;
        elsif p'last+1-i < p'last-r+1 then -- in zone B
          put(file,"checker p("); put(file,i,1); put(file,") = ");
          put(file,p(i),1); put_line(file," is in zone B");
          if locmap(integer32(p(i)),s) = 2 then
            put(file,"checking locmap("); put(file,p(i),1);
            put(file,","); put(file,s,1); put(file,") = 2");
            ind := Checker_Localization_Patterns.Rank
                     (locmap,integer32(p(i)),s);
            t.dg(ind) := 1;
            put(file,"  ind = "); put(file,ind,1); new_line(file);
            if not empty_zone_A
             then t.dg(piv) := 1;
            end if;
            x(integer32(p(i)),s) := Create(t);
            t.dg(ind) := 0; t.dg(piv) := 0;
          end if;
        elsif locmap(integer32(p(i)),s) = 2 then -- do not forget variables!
          ind := Checker_Localization_Patterns.Rank
                   (locmap,integer32(p(i)),s);
          if ind in t.dg'range then
            t.dg(ind) := 1;
            if not empty_zone_A
             then t.dg(piv) := 1;
            end if;
            put(file,"assigning forgotten variable x(");
            put(file,integer32(p(i)),1); put(file,",");
            put(file,s,1); put_line(file,")");
            x(integer32(p(i)),s) := Create(t);
            t.dg(ind) := 0; t.dg(piv) := 0;
          end if;
        end if;
      end if;
    end loop;
   -- if locmap(r+1,s+1) = 2 then   -- only assign to free position!
      t.dg(piv) := 1;
      x(r+1,s+1) := Create(t);    -- x(r+1,s+1)*m(r+1) for column s+1
      t.dg(piv) := 0;
   -- end if;
    x(big_r,s+1) := Create(t);  -- m(R) in column s+1 
    for i in p'range loop
      if integer32(p(i)) /= r and integer32(p(i)) /= r+1
                              and integer32(p(i)) /= big_r then
        --if locmap(integer32(p(i)),s+1) = 2 then -- was guarded before!
          ind := Checker_Localization_Patterns.Rank
                   (locmap,integer32(p(i)),s+1);
        if ind in t.dg'range then
          t.dg(ind) := 1;
          x(integer32(p(i)),s+1) := Create(t);
          t.dg(ind) := 0;
          if locmap(integer32(p(i)),s+1) /= 2 then
            put(file,"found index = "); put(file,ind,1);
            put(file," for p("); put(file,i,1); put(file,") = ");
            put(file,p(i),1); put(file," and s+1 = "); put(file,s+1,1);
            new_line(file);
          end if;
        else
          put(file,"failed index = "); put(file,ind,1);
          put(file," for p("); put(file,i,1); put(file,") = ");
          put(file,p(i),1); put(file," and s+1 = "); put(file,s+1,1);
          new_line(file);
          failed := true;
        end if;
        --end if; -- placed if again
      end if;
    end loop;
   -- if failed then
   --   put_line(file,"the localization map for failed indices :");
   --   put(file,locmap);
   -- end if;
    put_line(file,"the localization map : "); put(file,locmap);
    put_line(file,"the polynomial matrix of indeterminates :"); put(file,x);
    Clear(t);
  end First_Swap_Plane;

  procedure First_Swap_Plane
              ( file : in file_type;
                x : in out DoblDobl_Complex_Poly_Matrices.Matrix;
                r,big_r,dc,s : in integer32;
                p : in Standard_Natural_Vectors.Vector;
                locmap : in Standard_Natural_Matrices.Matrix ) is

    use DoblDobl_Complex_Numbers;
    use DoblDobl_Complex_Polynomials;

    dim : constant natural32 := Degree_of_Freedom(locmap);
    np1 : constant integer32 := integer32(dim)+1;
    ind : integer32;
    piv : constant integer32
        := Checker_Localization_Patterns.Rank(locmap,r+1,s+1);
    t : Term;
    empty_zone_A : constant boolean := Is_Zone_A_Empty(locmap,p,r,s,dc);

  begin
   -- put_line(file,"the localization map : "); put(file,locmap);
    put(file,"dim = "); put(file,dim,1); new_line(file);
    put(file,"np1 = "); put(file,np1,1); new_line(file);
    t.dg := new Standard_Natural_Vectors.Vector'(1..np1 => 0);
    t.cf := Create(integer(1));
    if not empty_zone_A then
      put_line(file,"empty_zone_A is false");
      t.dg(piv) := 1;   -- multiply by y(r+1,s+1)
      put(file,"multiply by x(");
      put(file,r+1,1); put(file,",");
      put(file,s+1,1); put(file,")");
      put(file," piv = "); put(file,piv,1); new_line(file);
    else
      put_line(file,"empty_zone_A is true");
    end if;
    x(r,s) := Create(t);    -- m(r) of 1st swapped column s
    t.dg(np1) := 1; t.cf := Create(integer(-1));
    x(r+1,s) := Create(t);  -- t*m(r+1) of 1st swapped column s
    t.dg(np1) := 0; t.dg(piv) := 0; t.cf := Create(integer(1));
    for i in p'range loop
      if integer32(p(i)) < r then      -- in zones A and B
        if p'last+1-i > p'last-dc+1 then    -- in zone A
          if locmap(integer32(p(i)),s+1) = 2 then
            ind := Checker_Localization_Patterns.Rank
                     (locmap,integer32(p(i)),s+1);
            t.dg(ind) := 1; t.dg(np1) := 1; t.cf := Create(integer(-1));
            x(integer32(p(i)),s) := Create(t);
            t.dg(ind) := 0; t.dg(np1) := 0; t.cf := Create(integer(1));
          end if;
        elsif p'last+1-i < p'last-r+1 then -- in zone B
          put(file,"checker p("); put(file,i,1); put(file,") = ");
          put(file,p(i),1); put_line(file," is in zone B");
          if locmap(integer32(p(i)),s) = 2 then
            put(file,"checking locmap("); put(file,p(i),1);
            put(file,","); put(file,s,1); put(file,") = 2");
            ind := Checker_Localization_Patterns.Rank
                     (locmap,integer32(p(i)),s);
            t.dg(ind) := 1;
            put(file,"  ind = "); put(file,ind,1); new_line(file);
            if not empty_zone_A
             then t.dg(piv) := 1;
            end if;
            x(integer32(p(i)),s) := Create(t);
            t.dg(ind) := 0; t.dg(piv) := 0;
          end if;
        elsif locmap(integer32(p(i)),s) = 2 then -- do not forget variables!
          ind := Checker_Localization_Patterns.Rank
                   (locmap,integer32(p(i)),s);
          if ind in t.dg'range then
            t.dg(ind) := 1;
            if not empty_zone_A
             then t.dg(piv) := 1;
            end if;
            put(file,"assigning forgotten variable x(");
            put(file,integer32(p(i)),1); put(file,",");
            put(file,s,1); put_line(file,")");
            x(integer32(p(i)),s) := Create(t);
            t.dg(ind) := 0; t.dg(piv) := 0;
          end if;
        end if;
      end if;
    end loop;
    t.dg(piv) := 1;
    x(r+1,s+1) := Create(t);    -- x(r+1,s+1)*m(r+1) for column s+1
    t.dg(piv) := 0;
    x(big_r,s+1) := Create(t);  -- m(R) in column s+1 
    for i in p'range loop
      if integer32(p(i)) /= r and integer32(p(i)) /= r+1
                              and integer32(p(i)) /= big_r then
       -- if locmap(p(i),s+1) = 2 then
          ind := Checker_Localization_Patterns.Rank
                   (locmap,integer32(p(i)),s+1);
        if ind in t.dg'range then
          t.dg(ind) := 1;
          x(integer32(p(i)),s+1) := Create(t);
          t.dg(ind) := 0;
          if locmap(integer32(p(i)),s+1) /= 2 then
            put(file,"found index = "); put(file,ind,1);
            put(file," for p("); put(file,i,1); put(file,") = ");
            put(file,p(i),1); put(file," and s+1 = "); put(file,s+1,1);
            new_line(file);
          end if;
        else
          put(file,"failed index = "); put(file,ind,1);
          put(file," for p("); put(file,i,1); put(file,") = ");
          put(file,p(i),1); put(file," and s+1 = "); put(file,s+1,1);
          new_line(file);
        end if;
       -- end if;
      end if;
    end loop;
    Clear(t);
  end First_Swap_Plane;

  procedure First_Swap_Plane
              ( file : in file_type;
                x : in out QuadDobl_Complex_Poly_Matrices.Matrix;
                r,big_r,dc,s : in integer32;
                p : in Standard_Natural_Vectors.Vector;
                locmap : in Standard_Natural_Matrices.Matrix ) is

    use QuadDobl_Complex_Numbers;
    use QuadDobl_Complex_Polynomials;

    dim : constant natural32 := Degree_of_Freedom(locmap);
    np1 : constant integer32 := integer32(dim)+1;
    ind : integer32;
    piv : constant integer32
        := Checker_Localization_Patterns.Rank(locmap,r+1,s+1);
    t : Term;
    empty_zone_A : constant boolean := Is_Zone_A_Empty(locmap,p,r,s,dc);

  begin
   -- put_line(file,"the localization map : "); put(file,locmap);
    put(file,"dim = "); put(file,dim,1); new_line(file);
    put(file,"np1 = "); put(file,np1,1); new_line(file);
    t.dg := new Standard_Natural_Vectors.Vector'(1..np1 => 0);
    t.cf := Create(integer(1));
    if not empty_zone_A then
      put_line(file,"empty_zone_A is false");
      t.dg(piv) := 1;   -- multiply by y(r+1,s+1)
      put(file,"multiply by x(");
      put(file,r+1,1); put(file,",");
      put(file,s+1,1); put(file,")");
      put(file," piv = "); put(file,piv,1); new_line(file);
    else
      put_line(file,"empty_zone_A is true");
    end if;
    x(r,s) := Create(t);    -- m(r) of 1st swapped column s
    t.dg(np1) := 1; t.cf := Create(integer(-1));
    x(r+1,s) := Create(t);  -- t*m(r+1) of 1st swapped column s
    t.dg(np1) := 0; t.dg(piv) := 0; t.cf := Create(integer(1));
    for i in p'range loop
      if integer32(p(i)) < r then      -- in zones A and B
        if p'last+1-i > p'last-dc+1 then    -- in zone A
          if locmap(integer32(p(i)),s+1) = 2 then
            ind := Checker_Localization_Patterns.Rank
                     (locmap,integer32(p(i)),s+1);
            t.dg(ind) := 1; t.dg(np1) := 1; t.cf := Create(integer(-1));
            x(integer32(p(i)),s) := Create(t);
            t.dg(ind) := 0; t.dg(np1) := 0; t.cf := Create(integer(1));
          end if;
        elsif p'last+1-i < p'last-r+1 then -- in zone B
          put(file,"checker p("); put(file,i,1); put(file,") = ");
          put(file,p(i),1); put_line(file," is in zone B");
          if locmap(integer32(p(i)),s) = 2 then
            put(file,"checking locmap("); put(file,p(i),1);
            put(file,","); put(file,s,1); put(file,") = 2");
            ind := Checker_Localization_Patterns.Rank
                     (locmap,integer32(p(i)),s);
            t.dg(ind) := 1;
            put(file,"  ind = "); put(file,ind,1); new_line(file);
            if not empty_zone_A
             then t.dg(piv) := 1;
            end if;
            x(integer32(p(i)),s) := Create(t);
            t.dg(ind) := 0; t.dg(piv) := 0;
          end if;
        elsif locmap(integer32(p(i)),s) = 2 then -- do not forget variables!
          ind := Checker_Localization_Patterns.Rank
                   (locmap,integer32(p(i)),s);
          if ind in t.dg'range then
            t.dg(ind) := 1;
            if not empty_zone_A
             then t.dg(piv) := 1;
            end if;
            put(file,"assigning forgotten variable x(");
            put(file,integer32(p(i)),1); put(file,",");
            put(file,s,1); put_line(file,")");
            x(integer32(p(i)),s) := Create(t);
            t.dg(ind) := 0; t.dg(piv) := 0;
          end if;
        end if;
      end if;
    end loop;
    t.dg(piv) := 1;
    x(r+1,s+1) := Create(t);    -- x(r+1,s+1)*m(r+1) for column s+1
    t.dg(piv) := 0;
    x(big_r,s+1) := Create(t);  -- m(R) in column s+1 
    for i in p'range loop
      if integer32(p(i)) /= r and integer32(p(i)) /= r+1
                              and integer32(p(i)) /= big_r then
       -- if locmap(p(i),s+1) = 2 then
          ind := Checker_Localization_Patterns.Rank
                   (locmap,integer32(p(i)),s+1);
        if ind in t.dg'range then
          t.dg(ind) := 1;
          x(integer32(p(i)),s+1) := Create(t);
          t.dg(ind) := 0;
          if locmap(integer32(p(i)),s+1) /= 2 then
            put(file,"found index = "); put(file,ind,1);
            put(file," for p("); put(file,i,1); put(file,") = ");
            put(file,p(i),1); put(file," and s+1 = "); put(file,s+1,1);
            new_line(file);
          end if;
        else
          put(file,"failed index = "); put(file,ind,1);
          put(file," for p("); put(file,i,1); put(file,") = ");
          put(file,p(i),1); put(file," and s+1 = "); put(file,s+1,1);
          new_line(file);
        end if;
       -- end if;
      end if;
    end loop;
    Clear(t);
  end First_Swap_Plane;

  procedure Second_Swap_Plane
              ( x : in out Standard_Complex_Poly_Matrices.Matrix;
                r,dc,s : in integer32;
                p : in Standard_Natural_Vectors.Vector;
                locmap : in Standard_Natural_Matrices.Matrix ) is

    use Standard_Complex_Numbers;
    use Standard_Complex_Polynomials;

    dim : constant natural32 := Degree_of_Freedom(locmap);
    np1 : constant integer32 := integer32(dim)+1;
    ind : integer32;
    t : Term;

  begin
    t.dg := new Standard_Natural_Vectors.Vector'(1..np1 => 0);
    t.cf := Create(1.0);
    x(r,s) := Create(t);    -- m(r) of 1st swapped column s
    t.dg(np1) := 1; t.cf := Create(-1.0);
    x(r+1,s) := Create(t);  -- t*m(r+1) of 1st swapped column s
    t.dg(np1) := 0; t.cf := Create(1.0);
    for i in p'range loop
      if integer32(p(i)) < r then      -- in zones A and B
        if p'last+1-i > p'last-dc+1 then    -- in zone A
          if locmap(integer32(p(i)),s+1) = 2 then
            ind := Checker_Localization_Patterns.Rank
                     (locmap,integer32(p(i)),s+1);
            t.dg(ind) := 1; t.dg(np1) := 1; t.cf := Create(-1.0);
            x(integer32(p(i)),s) := Create(t);
            t.dg(ind) := 0; t.dg(np1) := 0; t.cf := Create(1.0);
          end if;
        elsif p'last+1-i < p'last-r+1 then -- in zone B
          if locmap(integer32(p(i)),s) = 2 then
            ind := Checker_Localization_Patterns.Rank
                     (locmap,integer32(p(i)),s);
            t.dg(ind) := 1;
            x(integer32(p(i)),s) := Create(t);
            t.dg(ind) := 0;
          end if;
        elsif locmap(integer32(p(i)),s) = 2 then -- do not forget variables!
          ind := Checker_Localization_Patterns.Rank
                   (locmap,integer32(p(i)),s);
          if ind in t.dg'range then
            t.dg(ind) := 1;
            x(integer32(p(i)),s) := Create(t);
            t.dg(ind) := 0;
          end if;
        end if;
      end if;
    end loop;
    x(r+1,s+1) := Create(t);  -- m(r+1) in 2nd swapped column s+1
    for i in p'range loop
      if integer32(p(i)) < r then      -- in zones A and B
        if locmap(integer32(p(i)),s+1) = 2 then
          ind := Checker_Localization_Patterns.Rank
                   (locmap,integer32(p(i)),s+1);
          t.dg(ind) := 1;
          x(integer32(p(i)),s+1) := Create(t);
          t.dg(ind) := 0;
        end if;
      end if;
    end loop;
    Clear(t);
  end Second_Swap_Plane;

  procedure Second_Swap_Plane
              ( x : in out DoblDobl_Complex_Poly_Matrices.Matrix;
                r,dc,s : in integer32;
                p : in Standard_Natural_Vectors.Vector;
                locmap : in Standard_Natural_Matrices.Matrix ) is

    use DoblDobl_Complex_Numbers;
    use DoblDobl_Complex_Polynomials;

    dim : constant natural32 := Degree_of_Freedom(locmap);
    np1 : constant integer32 := integer32(dim)+1;
    ind : integer32;
    t : Term;

  begin
    t.dg := new Standard_Natural_Vectors.Vector'(1..np1 => 0);
    t.cf := Create(integer(1));
    x(r,s) := Create(t);    -- m(r) of 1st swapped column s
    t.dg(np1) := 1; t.cf := Create(integer(-1));
    x(r+1,s) := Create(t);  -- t*m(r+1) of 1st swapped column s
    t.dg(np1) := 0; t.cf := Create(integer(1));
    for i in p'range loop
      if integer32(p(i)) < r then      -- in zones A and B
        if p'last+1-i > p'last-dc+1 then    -- in zone A
          if locmap(integer32(p(i)),s+1) = 2 then
            ind := Checker_Localization_Patterns.Rank
                     (locmap,integer32(p(i)),s+1);
            t.dg(ind) := 1; t.dg(np1) := 1; t.cf := Create(integer(-1));
            x(integer32(p(i)),s) := Create(t);
            t.dg(ind) := 0; t.dg(np1) := 0; t.cf := Create(integer(1));
          end if;
        elsif p'last+1-i < p'last-r+1 then -- in zone B
          if locmap(integer32(p(i)),s) = 2 then
            ind := Checker_Localization_Patterns.Rank
                     (locmap,integer32(p(i)),s);
            t.dg(ind) := 1;
            x(integer32(p(i)),s) := Create(t);
            t.dg(ind) := 0;
          end if;
        elsif locmap(integer32(p(i)),s) = 2 then -- do not forget variables!
          ind := Checker_Localization_Patterns.Rank
                   (locmap,integer32(p(i)),s);
          if ind in t.dg'range then
            t.dg(ind) := 1;
            x(integer32(p(i)),s) := Create(t);
            t.dg(ind) := 0;
          end if;
        end if;
      end if;
    end loop;
    x(r+1,s+1) := Create(t);  -- m(r+1) in 2nd swapped column s+1
    for i in p'range loop
      if integer32(p(i)) < r then      -- in zones A and B
        if locmap(integer32(p(i)),s+1) = 2 then
          ind := Checker_Localization_Patterns.Rank
                   (locmap,integer32(p(i)),s+1);
          t.dg(ind) := 1;
          x(integer32(p(i)),s+1) := Create(t);
          t.dg(ind) := 0;
        end if;
      end if;
    end loop;
    Clear(t);
  end Second_Swap_Plane;

  procedure Second_Swap_Plane
              ( x : in out QuadDobl_Complex_Poly_Matrices.Matrix;
                r,dc,s : in integer32;
                p : in Standard_Natural_Vectors.Vector;
                locmap : in Standard_Natural_Matrices.Matrix ) is

    use QuadDobl_Complex_Numbers;
    use QuadDobl_Complex_Polynomials;

    dim : constant natural32 := Degree_of_Freedom(locmap);
    np1 : constant integer32 := integer32(dim)+1;
    ind : integer32;
    t : Term;

  begin
    t.dg := new Standard_Natural_Vectors.Vector'(1..np1 => 0);
    t.cf := Create(integer(1));
    x(r,s) := Create(t);    -- m(r) of 1st swapped column s
    t.dg(np1) := 1; t.cf := Create(integer(-1));
    x(r+1,s) := Create(t);  -- t*m(r+1) of 1st swapped column s
    t.dg(np1) := 0; t.cf := Create(integer(1));
    for i in p'range loop
      if integer32(p(i)) < r then      -- in zones A and B
        if p'last+1-i > p'last-dc+1 then    -- in zone A
          if locmap(integer32(p(i)),s+1) = 2 then
            ind := Checker_Localization_Patterns.Rank
                     (locmap,integer32(p(i)),s+1);
            t.dg(ind) := 1; t.dg(np1) := 1; t.cf := Create(integer(-1));
            x(integer32(p(i)),s) := Create(t);
            t.dg(ind) := 0; t.dg(np1) := 0; t.cf := Create(integer(1));
          end if;
        elsif p'last+1-i < p'last-r+1 then -- in zone B
          if locmap(integer32(p(i)),s) = 2 then
            ind := Checker_Localization_Patterns.Rank
                     (locmap,integer32(p(i)),s);
            t.dg(ind) := 1;
            x(integer32(p(i)),s) := Create(t);
            t.dg(ind) := 0;
          end if;
        elsif locmap(integer32(p(i)),s) = 2 then -- do not forget variables!
          ind := Checker_Localization_Patterns.Rank
                   (locmap,integer32(p(i)),s);
          if ind in t.dg'range then
            t.dg(ind) := 1;
            x(integer32(p(i)),s) := Create(t);
            t.dg(ind) := 0;
          end if;
        end if;
      end if;
    end loop;
    x(r+1,s+1) := Create(t);  -- m(r+1) in 2nd swapped column s+1
    for i in p'range loop
      if integer32(p(i)) < r then      -- in zones A and B
        if locmap(integer32(p(i)),s+1) = 2 then
          ind := Checker_Localization_Patterns.Rank
                   (locmap,integer32(p(i)),s+1);
          t.dg(ind) := 1;
          x(integer32(p(i)),s+1) := Create(t);
          t.dg(ind) := 0;
        end if;
      end if;
    end loop;
    Clear(t);
  end Second_Swap_Plane;

  procedure Second_Swap_Plane
              ( file : in file_type;
                x : in out Standard_Complex_Poly_Matrices.Matrix;
                r,dc,s : in integer32;
                p : in Standard_Natural_Vectors.Vector;
                locmap : in Standard_Natural_Matrices.Matrix ) is

    use Standard_Complex_Numbers;
    use Standard_Complex_Polynomials;

    dim : constant natural32 := Degree_of_Freedom(locmap);
    np1 : constant integer32 := integer32(dim)+1;
    ind : integer32;
    t : Term;

  begin
   -- put_line(file,"the localization map : "); put(file,locmap);
    put(file,"dim = "); put(file,dim,1);
    put(file,"  np1 = "); put(file,np1,1); new_line(file);
    t.dg := new Standard_Natural_Vectors.Vector'(1..np1 => 0);
    t.cf := Create(1.0);
    x(r,s) := Create(t);    -- m(r) of 1st swapped column s
    t.dg(np1) := 1; t.cf := Create(-1.0);
    x(r+1,s) := Create(t);  -- t*m(r+1) of 1st swapped column s
    t.dg(np1) := 0; t.cf := Create(1.0);
    for i in p'range loop
      if integer32(p(i)) < r then      -- in zones A and B
        if p'last+1-i > p'last-dc+1 then    -- in zone A
          put(file,"checker "); put(file,p(i),1); put_line(file," in zone A");
          if locmap(integer32(p(i)),s+1) = 2 then
            ind := Checker_Localization_Patterns.Rank
                     (locmap,integer32(p(i)),s+1);
            put(file," -> assigning to x("); put(file,p(i),1); put(file,",");
            put(file,s,1); put_line(file,")...");
            t.dg(ind) := 1; t.dg(np1) := 1; t.cf := Create(-1.0);
            x(integer32(p(i)),s) := Create(t);
            t.dg(ind) := 0; t.dg(np1) := 0; t.cf := Create(1.0);
          end if;
        elsif p'last+1-i < p'last-r+1 then -- in zone B
          if locmap(integer32(p(i)),s) = 2 then
            ind := Checker_Localization_Patterns.Rank
                     (locmap,integer32(p(i)),s);
            t.dg(ind) := 1;
            x(integer32(p(i)),s) := Create(t);
            t.dg(ind) := 0;
          end if;
        elsif locmap(integer32(p(i)),s) = 2 then -- do not forget variables!
          ind := Checker_Localization_Patterns.Rank
                   (locmap,integer32(p(i)),s);
          if ind in t.dg'range then
            t.dg(ind) := 1;
            x(integer32(p(i)),s) := Create(t);
            t.dg(ind) := 0;
          end if;
        end if;
      end if;
    end loop;
    x(r+1,s+1) := Create(t);  -- m(r+1) in 2nd swapped column s+1
    for i in p'range loop
      if integer32(p(i)) < r then      -- in zones A and B
        if locmap(integer32(p(i)),s+1) = 2 then
          ind := Checker_Localization_Patterns.Rank
                   (locmap,integer32(p(i)),s+1);
          t.dg(ind) := 1;
          x(integer32(p(i)),s+1) := Create(t);
          t.dg(ind) := 0;
        end if;
      end if;
    end loop;
    Clear(t);
  end Second_Swap_Plane;

  procedure Second_Swap_Plane
              ( file : in file_type;
                x : in out DoblDobl_Complex_Poly_Matrices.Matrix;
                r,dc,s : in integer32;
                p : in Standard_Natural_Vectors.Vector;
                locmap : in Standard_Natural_Matrices.Matrix ) is

    use DoblDobl_Complex_Numbers;
    use DoblDobl_Complex_Polynomials;

    dim : constant natural32 := Degree_of_Freedom(locmap);
    np1 : constant integer32 := integer32(dim)+1;
    ind : integer32;
    t : Term;

  begin
   -- put_line(file,"the localization map : "); put(file,locmap);
    put(file,"dim = "); put(file,dim,1);
    put(file,"  np1 = "); put(file,np1,1); new_line(file);
    t.dg := new Standard_Natural_Vectors.Vector'(1..np1 => 0);
    t.cf := Create(integer(1));
    x(r,s) := Create(t);    -- m(r) of 1st swapped column s
    t.dg(np1) := 1; t.cf := Create(integer(-1));
    x(r+1,s) := Create(t);  -- t*m(r+1) of 1st swapped column s
    t.dg(np1) := 0; t.cf := Create(integer(1));
    for i in p'range loop
      if integer32(p(i)) < r then      -- in zones A and B
        if p'last+1-i > p'last-dc+1 then    -- in zone A
          put(file,"checker "); put(file,p(i),1); put_line(file," in zone A");
          if locmap(integer32(p(i)),s+1) = 2 then
            ind := Checker_Localization_Patterns.Rank
                     (locmap,integer32(p(i)),s+1);
            put(file," -> assigning to x("); put(file,p(i),1); put(file,",");
            put(file,s,1); put_line(file,")...");
            t.dg(ind) := 1; t.dg(np1) := 1; t.cf := Create(integer(-1));
            x(integer32(p(i)),s) := Create(t);
            t.dg(ind) := 0; t.dg(np1) := 0; t.cf := Create(integer(1));
          end if;
        elsif p'last+1-i < p'last-r+1 then -- in zone B
          if locmap(integer32(p(i)),s) = 2 then
            ind := Checker_Localization_Patterns.Rank
                     (locmap,integer32(p(i)),s);
            t.dg(ind) := 1;
            x(integer32(p(i)),s) := Create(t);
            t.dg(ind) := 0;
          end if;
        elsif locmap(integer32(p(i)),s) = 2 then -- do not forget variables!
          ind := Checker_Localization_Patterns.Rank
                   (locmap,integer32(p(i)),s);
          if ind in t.dg'range then
            t.dg(ind) := 1;
            x(integer32(p(i)),s) := Create(t);
            t.dg(ind) := 0;
          end if;
        end if;
      end if;
    end loop;
    x(r+1,s+1) := Create(t);  -- m(r+1) in 2nd swapped column s+1
    for i in p'range loop
      if integer32(p(i)) < r then      -- in zones A and B
        if locmap(integer32(p(i)),s+1) = 2 then
          ind := Checker_Localization_Patterns.Rank
                   (locmap,integer32(p(i)),s+1);
          t.dg(ind) := 1;
          x(integer32(p(i)),s+1) := Create(t);
          t.dg(ind) := 0;
        end if;
      end if;
    end loop;
    Clear(t);
  end Second_Swap_Plane;

  procedure Second_Swap_Plane
              ( file : in file_type;
                x : in out QuadDobl_Complex_Poly_Matrices.Matrix;
                r,dc,s : in integer32;
                p : in Standard_Natural_Vectors.Vector;
                locmap : in Standard_Natural_Matrices.Matrix ) is

    use QuadDobl_Complex_Numbers;
    use QuadDobl_Complex_Polynomials;

    dim : constant natural32 := Degree_of_Freedom(locmap);
    np1 : constant integer32 := integer32(dim)+1;
    ind : integer32;
    t : Term;

  begin
   -- put_line(file,"the localization map : "); put(file,locmap);
    put(file,"dim = "); put(file,dim,1);
    put(file,"  np1 = "); put(file,np1,1); new_line(file);
    t.dg := new Standard_Natural_Vectors.Vector'(1..np1 => 0);
    t.cf := Create(integer(1));
    x(r,s) := Create(t);    -- m(r) of 1st swapped column s
    t.dg(np1) := 1; t.cf := Create(integer(-1));
    x(r+1,s) := Create(t);  -- t*m(r+1) of 1st swapped column s
    t.dg(np1) := 0; t.cf := Create(integer(1));
    for i in p'range loop
      if integer32(p(i)) < r then      -- in zones A and B
        if p'last+1-i > p'last-dc+1 then    -- in zone A
          put(file,"checker "); put(file,p(i),1); put_line(file," in zone A");
          if locmap(integer32(p(i)),s+1) = 2 then
            ind := Checker_Localization_Patterns.Rank
                     (locmap,integer32(p(i)),s+1);
            put(file," -> assigning to x("); put(file,p(i),1); put(file,",");
            put(file,s,1); put_line(file,")...");
            t.dg(ind) := 1; t.dg(np1) := 1; t.cf := Create(integer(-1));
            x(integer32(p(i)),s) := Create(t);
            t.dg(ind) := 0; t.dg(np1) := 0; t.cf := Create(integer(1));
          end if;
        elsif p'last+1-i < p'last-r+1 then -- in zone B
          if locmap(integer32(p(i)),s) = 2 then
            ind := Checker_Localization_Patterns.Rank
                     (locmap,integer32(p(i)),s);
            t.dg(ind) := 1;
            x(integer32(p(i)),s) := Create(t);
            t.dg(ind) := 0;
          end if;
        elsif locmap(integer32(p(i)),s) = 2 then -- do not forget variables!
          ind := Checker_Localization_Patterns.Rank
                   (locmap,integer32(p(i)),s);
          if ind in t.dg'range then
            t.dg(ind) := 1;
            x(integer32(p(i)),s) := Create(t);
            t.dg(ind) := 0;
          end if;
        end if;
      end if;
    end loop;
    x(r+1,s+1) := Create(t);  -- m(r+1) in 2nd swapped column s+1
    for i in p'range loop
      if integer32(p(i)) < r then      -- in zones A and B
        if locmap(integer32(p(i)),s+1) = 2 then
          ind := Checker_Localization_Patterns.Rank
                   (locmap,integer32(p(i)),s+1);
          t.dg(ind) := 1;
          x(integer32(p(i)),s+1) := Create(t);
          t.dg(ind) := 0;
        end if;
      end if;
    end loop;
    Clear(t);
  end Second_Swap_Plane;

  function Swap_Column
              ( r : integer32; m : Standard_Natural_Matrices.Matrix )
              return integer32 is
  begin
    for j in m'first(2)..m'last(2)-1 loop
      if m(r,j) = 1
       then return j;
      end if;
    end loop;
    return 0;
  end Swap_Column;

  function Swap_Column
              ( r : integer32; rows : Standard_Natural_Vectors.Vector )
              return integer32 is

  -- NOTE : rows of checkers are labeled in reverse order,
  --   i.e.: the row of the 1st checker is in row(row'last).

    ind : integer32 := 0;

  begin
    for i in reverse rows'range loop
      ind := ind + 1;
      if rows(i) = natural32(r)
       then return ind-1; -- we find checker that has been swapped
      end if;
    end loop;
    return 0;
  end Swap_Column;

  function Swap_Checker ( q,rows,cols : Standard_Natural_Vectors.Vector )
                        return integer32 is

    n : constant integer32 := q'last;
    dc : constant integer32 := Descending_Checker(q);
    rc : constant integer32 := Rising_Checker(q,dc);
    cd : constant integer32
       := Top_White_Checker(integer32(q(rc)),n-rc+1,n,rows,cols);
    res : constant natural32 := rows(cd);

  begin
    return integer32(res);
  end Swap_Checker;

  function Stay_Moving_Plane
              ( n,k,r : in integer32;
                p,rows,cols : in Standard_Natural_Vectors.Vector ) 
              return Standard_Complex_Poly_Matrices.Matrix is

    use Standard_Complex_Numbers;
    use Standard_Complex_Polynomials;

    res : Standard_Complex_Poly_Matrices.Matrix(1..n,1..k);
    locmap : constant Standard_Natural_Matrices.Matrix(1..n,1..k)
           := Checker_Localization_Patterns.Column_Pattern(n,k,p,rows,cols);
    dim : constant integer32
        := integer32(Checker_Localization_Patterns.Degree_of_Freedom(locmap));
    t : Term;

  begin
    Initialize_Moving_Plane(res,locmap);
    t.dg := new Standard_Natural_Vectors.Vector'(1..dim+1 => 0);
    t.dg(dim+1) := 1;
    t.cf := Create(-1.0);
    for j in locmap'range(2) loop  -- find spot to place t in row r+1
      if locmap(r,j) = 1
       then res(r+1,j) := Create(t); exit;
      end if;
    end loop;
    Clear(t);
    return res;
  end Stay_Moving_Plane;

  function Stay_Moving_Plane
              ( n,k,r : in integer32;
                p,rows,cols : in Standard_Natural_Vectors.Vector ) 
              return DoblDobl_Complex_Poly_Matrices.Matrix is

    use DoblDobl_Complex_Numbers;
    use DoblDobl_Complex_Polynomials;

    res : DoblDobl_Complex_Poly_Matrices.Matrix(1..n,1..k);
    locmap : constant Standard_Natural_Matrices.Matrix(1..n,1..k)
           := Checker_Localization_Patterns.Column_Pattern(n,k,p,rows,cols);
    dim : constant integer32
        := integer32(Checker_Localization_Patterns.Degree_of_Freedom(locmap));
    t : Term;

  begin
    Initialize_Moving_Plane(res,locmap);
    t.dg := new Standard_Natural_Vectors.Vector'(1..dim+1 => 0);
    t.dg(dim+1) := 1;
    t.cf := Create(integer(-1));
    for j in locmap'range(2) loop  -- find spot to place t in row r+1
      if locmap(r,j) = 1
       then res(r+1,j) := Create(t); exit;
      end if;
    end loop;
    Clear(t);
    return res;
  end Stay_Moving_Plane;

  function Stay_Moving_Plane
              ( n,k,r : in integer32;
                p,rows,cols : in Standard_Natural_Vectors.Vector ) 
              return QuadDobl_Complex_Poly_Matrices.Matrix is

    use QuadDobl_Complex_Numbers;
    use QuadDobl_Complex_Polynomials;

    res : QuadDobl_Complex_Poly_Matrices.Matrix(1..n,1..k);
    locmap : constant Standard_Natural_Matrices.Matrix(1..n,1..k)
           := Checker_Localization_Patterns.Column_Pattern(n,k,p,rows,cols);
    dim : constant integer32
        := integer32(Checker_Localization_Patterns.Degree_of_Freedom(locmap));
    t : Term;

  begin
    Initialize_Moving_Plane(res,locmap);
    t.dg := new Standard_Natural_Vectors.Vector'(1..dim+1 => 0);
    t.dg(dim+1) := 1;
    t.cf := Create(integer(-1));
    for j in locmap'range(2) loop  -- find spot to place t in row r+1
      if locmap(r,j) = 1
       then res(r+1,j) := Create(t); exit;
      end if;
    end loop;
    Clear(t);
    return res;
  end Stay_Moving_Plane;

  function Swap_Moving_Plane
              ( n,k,r,big_r,s : in integer32;
                q,p,rows,cols : in Standard_Natural_Vectors.Vector ) 
              return Standard_Complex_Poly_Matrices.Matrix is

    res : Standard_Complex_Poly_Matrices.Matrix(1..n,1..k);
    dc : constant integer32 := Descending_Checker(q);
    locmap : constant Standard_Natural_Matrices.Matrix(1..n,1..k)
           := Checker_Localization_Patterns.Column_Pattern(n,k,p,rows,cols);

  begin
    Initialize_Moving_Plane(res,locmap,s);
    if big_r = r + 1
     then Second_Swap_Plane(res,r,dc,s,p,locmap);
     else First_Swap_Plane(res,r,big_r,dc,s,p,locmap);
    end if;
    return res;
  end Swap_Moving_Plane;

  function Swap_Moving_Plane
              ( n,k,r,big_r,s : in integer32;
                q,p,rows,cols : in Standard_Natural_Vectors.Vector ) 
              return DoblDobl_Complex_Poly_Matrices.Matrix is

    res : DoblDobl_Complex_Poly_Matrices.Matrix(1..n,1..k);
    dc : constant integer32 := Descending_Checker(q);
    locmap : constant Standard_Natural_Matrices.Matrix(1..n,1..k)
           := Checker_Localization_Patterns.Column_Pattern(n,k,p,rows,cols);

  begin
    Initialize_Moving_Plane(res,locmap,s);
    if big_r = r + 1
     then Second_Swap_Plane(res,r,dc,s,p,locmap);
     else First_Swap_Plane(res,r,big_r,dc,s,p,locmap);
    end if;
    return res;
  end Swap_Moving_Plane;

  function Swap_Moving_Plane
              ( n,k,r,big_r,s : in integer32;
                q,p,rows,cols : in Standard_Natural_Vectors.Vector ) 
              return QuadDobl_Complex_Poly_Matrices.Matrix is

    res : QuadDobl_Complex_Poly_Matrices.Matrix(1..n,1..k);
    dc : constant integer32 := Descending_Checker(q);
    locmap : constant Standard_Natural_Matrices.Matrix(1..n,1..k)
           := Checker_Localization_Patterns.Column_Pattern(n,k,p,rows,cols);

  begin
    Initialize_Moving_Plane(res,locmap,s);
    if big_r = r + 1
     then Second_Swap_Plane(res,r,dc,s,p,locmap);
     else First_Swap_Plane(res,r,big_r,dc,s,p,locmap);
    end if;
    return res;
  end Swap_Moving_Plane;

  function Swap_Moving_Plane
              ( file : in file_type; n,k,r,big_r,s : in integer32;
                q,p,rows,cols : in Standard_Natural_Vectors.Vector ) 
              return Standard_Complex_Poly_Matrices.Matrix is

    res : Standard_Complex_Poly_Matrices.Matrix(1..n,1..k);
    dc : constant integer32 := Descending_Checker(q);
    rc : constant integer32 := Rising_Checker(q,dc);
    cd : constant integer32
       := Top_White_Checker(integer32(q(rc)),n-rc+1,n,rows,cols);
    locmap : constant Standard_Natural_Matrices.Matrix(1..n,1..k)
           := Checker_Localization_Patterns.Column_Pattern(n,k,p,rows,cols);

  begin
    put_line(file,"defining coordinates of moving plane for swap homotopy");
   -- put_line(file,"the localization map :"); put(file,locmap);
    put(file,"critical row : "); put(file,r,1);
    put(file,"  s = "); put(file,s,1);
    put(file,"  top white checker : "); put(file,cd,1);
    put(file,"  R : "); put(file,big_r,1); new_line(file);
    put(file,"ascending checker index : "); put(file,dc,1); new_line(file);
    -- descending in specializing, but ascending in generalizing poset
    put(file,"  p = "); put(file,p);
    put(file,"  rows = "); put(file,rows);
    put(file,"  cols = "); put(file,cols); new_line(file);
    Initialize_Moving_Plane(res,locmap,s);
   -- put_line(file,"After initialization of the moving plane :");
   -- put(file,res);
    put(file,"Black checkers in upper right (zone A) : ");
    for i in p'range loop
      if integer32(p(i)) < r and p'last+1-i > p'last-dc+1
       then put(file," "); put(file,p(i),1);
      end if;
    end loop;
    new_line(file);
    put(file,"Black checkers in upper left (zone B) : ");
    for i in p'range loop
      if integer32(p(i)) < r and p'last+1-i < p'last-r+1
       then put(file," "); put(file,p(i),1);
      end if;
    end loop;
    new_line(file);
    put(file,"Black checkers in critical diagonal (zone E) : ");
    for i in p'range loop
      if integer32(p(i)) > r + 1
       then put(file," "); put(file,p(i),1);
      end if;
    end loop;
    new_line(file);
    if big_r = r + 1 then
      put(file,"R(="); put(file,big_r,1); put(file,") = ");
      put(file,"r(="); put(file,r,1); put_line(file,") + 1");
      put_line(file,"second type of swap homotopy");
      Second_Swap_Plane(file,res,r,dc,s,p,locmap);
    else
      put(file,"R(="); put(file,big_r,1); put(file,") > ");
      put(file,"r(="); put(file,r,1); put_line(file,") + 1");
      put_line(file,"first type of swap homotopy");
      First_Swap_Plane(file,res,r,big_r,dc,s,p,locmap);
    end if;
    return res;
  end Swap_Moving_Plane;
 
  function Swap_Moving_Plane
              ( file : in file_type; n,k,r,big_r,s : in integer32;
                q,p,rows,cols : in Standard_Natural_Vectors.Vector ) 
              return DoblDobl_Complex_Poly_Matrices.Matrix is

    res : DoblDobl_Complex_Poly_Matrices.Matrix(1..n,1..k);
    dc : constant integer32 := Descending_Checker(q);
    rc : constant integer32 := Rising_Checker(q,dc);
    cd : constant integer32
       := Top_White_Checker(integer32(q(rc)),n-rc+1,n,rows,cols);
    locmap : constant Standard_Natural_Matrices.Matrix(1..n,1..k)
           := Checker_Localization_Patterns.Column_Pattern(n,k,p,rows,cols);

  begin
    put_line(file,"defining coordinates of moving plane for swap homotopy");
   -- put_line(file,"the localization map :"); put(file,locmap);
    put(file,"critical row : "); put(file,r,1);
    put(file,"  s = "); put(file,s,1);
    put(file,"  top white checker : "); put(file,cd,1);
    put(file,"  R : "); put(file,big_r,1); new_line(file);
    put(file,"ascending checker index : "); put(file,dc,1); new_line(file);
    -- descending in specializing, but ascending in generalizing poset
    put(file,"  p = "); put(file,p);
    put(file,"  rows = "); put(file,rows);
    put(file,"  cols = "); put(file,cols); new_line(file);
    Initialize_Moving_Plane(res,locmap,s);
   -- put_line(file,"After initialization of the moving plane :");
   -- put(file,res);
    put(file,"Black checkers in upper right (zone A) : ");
    for i in p'range loop
      if integer32(p(i)) < r and p'last+1-i > p'last-dc+1
       then put(file," "); put(file,p(i),1);
      end if;
    end loop;
    new_line(file);
    put(file,"Black checkers in upper left (zone B) : ");
    for i in p'range loop
      if integer32(p(i)) < r and p'last+1-i < p'last-r+1
       then put(file," "); put(file,p(i),1);
      end if;
    end loop;
    new_line(file);
    put(file,"Black checkers in critical diagonal (zone E) : ");
    for i in p'range loop
      if integer32(p(i)) > r + 1
       then put(file," "); put(file,p(i),1);
      end if;
    end loop;
    new_line(file);
    if big_r = r + 1 then
      put(file,"R(="); put(file,big_r,1); put(file,") = ");
      put(file,"r(="); put(file,r,1); put_line(file,") + 1");
      put_line(file,"second type of swap homotopy");
      Second_Swap_Plane(file,res,r,dc,s,p,locmap);
    else
      put(file,"R(="); put(file,big_r,1); put(file,") > ");
      put(file,"r(="); put(file,r,1); put_line(file,") + 1");
      put_line(file,"first type of swap homotopy");
      First_Swap_Plane(file,res,r,big_r,dc,s,p,locmap);
    end if;
    return res;
  end Swap_Moving_Plane;
 
  function Swap_Moving_Plane
              ( file : in file_type; n,k,r,big_r,s : in integer32;
                q,p,rows,cols : in Standard_Natural_Vectors.Vector ) 
              return QuadDobl_Complex_Poly_Matrices.Matrix is

    res : QuadDobl_Complex_Poly_Matrices.Matrix(1..n,1..k);
    dc : constant integer32 := Descending_Checker(q);
    rc : constant integer32 := Rising_Checker(q,dc);
    cd : constant integer32
       := Top_White_Checker(integer32(q(rc)),n-rc+1,n,rows,cols);
    locmap : constant Standard_Natural_Matrices.Matrix(1..n,1..k)
           := Checker_Localization_Patterns.Column_Pattern(n,k,p,rows,cols);

  begin
    put_line(file,"defining coordinates of moving plane for swap homotopy");
   -- put_line(file,"the localization map :"); put(file,locmap);
    put(file,"critical row : "); put(file,r,1);
    put(file,"  s = "); put(file,s,1);
    put(file,"  top white checker : "); put(file,cd,1);
    put(file,"  R : "); put(file,big_r,1); new_line(file);
    put(file,"ascending checker index : "); put(file,dc,1); new_line(file);
    -- descending in specializing, but ascending in generalizing poset
    put(file,"  p = "); put(file,p);
    put(file,"  rows = "); put(file,rows);
    put(file,"  cols = "); put(file,cols); new_line(file);
    Initialize_Moving_Plane(res,locmap,s);
   -- put_line(file,"After initialization of the moving plane :");
   -- put(file,res);
    put(file,"Black checkers in upper right (zone A) : ");
    for i in p'range loop
      if integer32(p(i)) < r and p'last+1-i > p'last-dc+1
       then put(file," "); put(file,p(i),1);
      end if;
    end loop;
    new_line(file);
    put(file,"Black checkers in upper left (zone B) : ");
    for i in p'range loop
      if integer32(p(i)) < r and p'last+1-i < p'last-r+1
       then put(file," "); put(file,p(i),1);
      end if;
    end loop;
    new_line(file);
    put(file,"Black checkers in critical diagonal (zone E) : ");
    for i in p'range loop
      if integer32(p(i)) > r + 1
       then put(file," "); put(file,p(i),1);
      end if;
    end loop;
    new_line(file);
    if big_r = r + 1 then
      put(file,"R(="); put(file,big_r,1); put(file,") = ");
      put(file,"r(="); put(file,r,1); put_line(file,") + 1");
      put_line(file,"second type of swap homotopy");
      Second_Swap_Plane(file,res,r,dc,s,p,locmap);
    else
      put(file,"R(="); put(file,big_r,1); put(file,") > ");
      put(file,"r(="); put(file,r,1); put_line(file,") + 1");
      put_line(file,"first type of swap homotopy");
      First_Swap_Plane(file,res,r,big_r,dc,s,p,locmap);
    end if;
    return res;
  end Swap_Moving_Plane;
 
end Checker_Homotopies;
