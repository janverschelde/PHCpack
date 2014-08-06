with Standard_Natural_Numbers_io;      use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;         use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;      use Standard_Integer_Numbers_io;
with Double_Double_Numbers;            use Double_Double_Numbers;
with DoblDobl_Complex_Numbers;         use DoblDobl_Complex_Numbers;
with Standard_Natural_Vectors_io;      use Standard_Natural_Vectors_io;
with Monomial_Hashing;                 use Monomial_Hashing;

package body DoblDobl_Deflation_Matrices is

  function Number_of_Columns
              ( d,nv,R1 : Standard_Natural_Vectors.Vector;
                m : natural32 ) return natural32 is

    res : natural32 := nv(integer32(m));

  begin
    for j in 1..d(0) loop
      res := res*nv(0);
    end loop;
    for i in 1..d'last loop
      for j in 1..d(i) loop
        res := res*R1(i);
      end loop;
    end loop;
    return res;
  end Number_of_Columns;

  function Zero_Matrix ( rows,cols : natural32 ) return Matrix is

    zero : constant double_double := create(0.0);

  begin
    if rows*cols > storage_threshold then
      raise STORAGE_ERROR;
    else
      declare
        res : Matrix(1..integer32(rows),1..integer32(cols));
      begin
        for i in res'range(1) loop
          for j in res'range(2) loop
            res(i,j) := Create(zero);
          end loop;
        end loop;
        return res;
      end;
    end if;
  end Zero_Matrix;

  procedure Multi_Loop
              ( s : in natural32; d,n : in Standard_Natural_Vectors.Vector ) is

    index : Standard_Natural_Vectors.Vector(1..integer32(s))
          := (1..integer32(s) => 0);
    accum : Standard_Natural_Vectors.Vector(d'range) := d;
    p : integer32 := 0;

    procedure Run_Multi_Loop ( k : integer32 ) is
    begin
      if k > accum'last then
        Loop_Body(index);
      elsif accum(k) = 0 then
        Run_Multi_Loop(k+1);
      else
        p := p+1;
        accum(k) := accum(k)-1;
        for i in 1..n(k) loop
          index(p) := i;
          Run_Multi_Loop(k);
        end loop;
        p := p-1;
        accum(k) := accum(k)+1;
      end if;
    end Run_Multi_Loop;

  begin
    Run_Multi_Loop(d'first);
  end Multi_Loop;

-- MULTIPLIERS :

  function Multiply ( B : VecMat; x : DoblDobl_Complex_VecVecs.VecVec )
                    return DoblDobl_Complex_VecVecs.VecVec is

    res : DoblDobl_Complex_VecVecs.VecVec(B'range);

  begin
    for i in res'range loop
      res(i) := new DoblDobl_Complex_Vectors.Vector'(B(i).all*x(i).all);
    end loop;
    return res;
  end Multiply;

  procedure Multiply ( A : in out Link_to_Matrix; r,c : in integer32;
                       JM,B : in Link_to_Matrix ) is

    acc : Complex_Number;

  begin
    for i in JM'range(1) loop
      for j in B'range(2) loop
        acc := JM(i,JM'first(2))*B(B'first(1),j);
        for k in JM'first(2)+1..B'last(1) loop
          acc := acc + JM(i,k)*B(k,j);
        end loop;
        A(r+i-1,c+j-1) := acc;
      end loop;
    end loop;
  end Multiply;

  procedure Multiply ( A : in out Link_to_Matrix; r,c : in integer32;
                       JM : in Link_to_Matrix;
                       Bl : in DoblDobl_Complex_Vectors.Link_to_Vector ) is

    acc : Complex_Number;
    m : constant integer32 := JM'last(2)/Bl'last;
    offset : integer32 := 0;

  begin
   -- put(" JM'last(1) = "); put(JM'last(1),1);
   -- put(" JM'last(2) = "); put(JM'last(2),1);
   -- put(" Bl'last = "); put(Bl'last,1); new_line;
    for k in 0..m-1 loop
      for i in JM'range(1) loop
        acc := JM(i,JM'first(2)+offset)*Bl(Bl'first);
        for j in Bl'first+1..Bl'last loop
          acc := acc + JM(i,j+offset)*Bl(j);
        end loop;
        A(r+i-1,c+k) := acc;
      end loop;
      offset := offset + Bl'last;
    end loop;
  end Multiply;

  procedure Multiply ( A : in out Link_to_Matrix; r,c,m : in integer32;
                       JM : in Link_to_Matrix;
                       Bl : in DoblDobl_Complex_Vectors.Link_to_Vector ) is

    acc : Complex_Number;
    offset : integer32 := 0;
    index : integer32;
    zero : constant double_double := create(0.0);

  begin
    put(" JM'last(1) = "); put(JM'last(1),1);
    put(" JM'last(2) = "); put(JM'last(2),1);
    put(" Bl'last = "); put(Bl'last,1);
    put(" m = "); put(m,1); new_line;
    for k in 0..m-1 loop
      put("doing column "); put(k+1,1); new_line;
      for i in JM'range(1) loop
        offset := k*m;
        index := Bl'first;
        acc := Create(zero);
        while index <= Bl'last loop
          for j in 1..m loop
            put("  index = "); put(index,1);
            put("  offset = "); put(offset,1);
            put("  j = "); put(j,1); new_line;
            acc := acc + JM(i,j+offset)*Bl(index);
            index := index + 1;
          end loop;
          offset := offset + Bl'last;
        end loop;
        A(r+i-1,c+k) := acc;
      end loop;
    end loop;
  end Multiply;

  procedure Alternating_Permute
              ( A : in out Link_to_Matrix;
                row,col,nbrows,nbcols,nbvars : in integer32 ) is

    buffer : Matrix(1..nbrows,1..nbcols);
    offset,index : integer32;

  begin
    if nbvars > 1 then
      index := 0;
      for k in 1..nbvars loop           -- select from A into buffer
        offset := k;
        while offset <= nbcols loop
          index := index + 1;
          for i in 1..nbrows loop
            if row+i-1 > A'last(1) or col-1+offset > A'last(2)
             then return;  -- PATCH
             else buffer(i,index) := A(row+i-1,col-1+offset);
            end if;
          end loop;
          offset := offset + nbvars;
        end loop;
      end loop;
      for i in 1..nbrows loop            -- copy buffer into A
        for j in 1..nbcols loop
          A(row+i-1,col+j-1) := buffer(i,j);
        end loop;
      end loop;
    end if;
  end Alternating_Permute;

  procedure Alternating_Permute
              ( file : in file_type; A : in out Link_to_Matrix;
                row,col,nbrows,nbcols,nbvars : in integer32 ) is

    buffer : Matrix(1..nbrows,1..nbcols);
    offset,index : integer32;

  begin
    put(file,"-> Alternating_Permute on A,");
    put(file," A'last(1) = "); put(file,A'last(1),1);
    put(file," A'last(2) = "); put(file,A'last(2),1); new_line(file);
    put(file," row = "); put(file,row,1);
    put(file," col = "); put(file,col,1);
    put(file," nbrows = "); put(file,nbrows,1);
    put(file," nbcols = "); put(file,nbcols,1);
    put(file," nbvars = "); put(file,nbvars,1); new_line(file);
    if nbvars > 1 then
      index := 0;
      for k in 1..nbvars loop           -- select from A into buffer
        offset := k;
        while offset <= nbcols loop
          index := index + 1;
          for i in 1..nbrows loop
            if row+i-1 > A'last(1) or col-1+offset > A'last(2)
             then return; -- PATCH
             else buffer(i,index) := A(row+i-1,col-1+offset);
            end if;
          end loop;
          offset := offset + nbvars;
        end loop;
      end loop;
      for i in 1..nbrows loop            -- copy buffer into A
        for j in 1..nbcols loop
          A(row+i-1,col+j-1) := buffer(i,j);
        end loop;
      end loop;
    end if;
  end Alternating_Permute;

  procedure One_Right_Multiply_Deflation
              ( A : in out Link_to_Matrix; m,r,c,k : in integer32;
                R0 : in Standard_Natural_Vectors.Vector;
                JM : in Link_to_Matrix;
                Bl : in DoblDobl_Complex_Vectors.Link_to_Vector ) is

    JM_offset,Bl_offset : integer32;
    acc : Complex_Number;
    zero : constant double_double := create(0.0);

  begin
    for i in JM'range(1) loop       -- run over all rows of JM
      JM_offset := 0;
      Bl_offset := 0;
      for j0 in 0..m-1 loop         -- m stages: m-column block structure
        for j1 in 1..integer32(R0(k-1)) loop
          acc := Create(zero);
          for j2 in 1..integer32(R0(j0)) loop
            acc := acc + JM(i,JM_offset+j2)*Bl(Bl_offset+j2);
          end loop;
          A(r+i-1,c+j1-1) := A(r+i-1,c+j1-1) + acc;
          JM_offset := JM_offset + integer32(R0(j0));
        end loop;
        Bl_offset := Bl_offset + integer32(R0(j0));
      end loop;
    end loop;
  end One_Right_Multiply_Deflation;

  procedure One_Right_Multiply_Deflation
              ( file : in file_type;
                A : in out Link_to_Matrix; m,r,c,k : in integer32;
                R0 : in Standard_Natural_Vectors.Vector;
                JM : in Link_to_Matrix;
                Bl : in DoblDobl_Complex_Vectors.Link_to_Vector ) is

    JM_offset,Bl_offset : integer32;
    acc : Complex_Number;
    zero : constant double_double := create(0.0);

  begin
    put(file,"-> executing One_Right_Multiply_Deflation on Bl, k = ");
    put(file,k,1); put(file,", r = ");
    put(file,r,1); put(file," and c = "); put(file,c,1); new_line(file);
    for i in JM'range(1) loop       -- run over all rows of JM
      JM_offset := 0;
      Bl_offset := 0;
      for j0 in 0..m-1 loop         -- m stages: m-column block structure
        for j1 in 1..integer32(R0(k-1)) loop
          acc := Create(zero);
          for j2 in 1..integer32(R0(j0)) loop
            acc := acc + JM(i,JM_offset+j2)*Bl(Bl_offset+j2);
          end loop;
          A(r+i-1,c+j1-1) := A(r+i-1,c+j1-1) + acc;
          JM_offset := JM_offset + integer32(R0(j0));
        end loop;
        Bl_offset := Bl_offset + integer32(R0(j0));
      end loop;
    end loop;
  end One_Right_Multiply_Deflation;

  function Is_Mixed ( d : Standard_Natural_Vectors.Vector ) return boolean is

  -- DESCRIPTION :
  --   The derivative operator d is mixed if more than one entry is nonzero.
  --   If d is mixed, then true is returned, otherwise false is returned.

    cnt : integer32 := 0;  -- counts nonzero elements in d

  begin
    for i in d'range loop
      if d(i) /= 0
       then cnt := cnt + 1;
      end if;
    end loop;
    return (cnt > 1);
  end Is_Mixed;

  function Permute ( A : Link_to_Matrix; nq,k : integer32;
                     R0 : Standard_Natural_Vectors.Vector )
                   return Link_to_Matrix is

  -- DESCRIPTION :
  --   Depending on the value of k, the Alternating_Permute is
  --   called in once, or iteratively.

  -- ON ENTRY :
  --   A         second-order derivative of A(1) of mixed type;
  --   nq        number of equations in the original system;
  --   k         label of the child in A(2);
  --   R0        number of variables added in each stage in the deflation.

  -- ON RETURN :
  --   an alternatingly permuted second-order derivative of A(1).

    res : Link_to_Matrix;
    nc,row,col : integer32;

  begin
    res := new Matrix'(A.all);
    row := A'first(1) + nq;
    if k = 1 then
      nc := integer32(R0(0)*R0(0)*R0(1));
      Alternating_Permute(res,row,A'first(2),nq,nc,integer32(R0(1)));
    else
      nc := integer32(R0(0)*R0(1));
      col := A'first(2);
      for i in 1..R0(0) loop
        Alternating_Permute(res,row,col,nq,nc,integer32(R0(1)));
        col := col + nc;
      end loop;
    end if;
    return res;
  end Permute;

  function Permute ( file : file_type;
                     A : Link_to_Matrix; nq,k : integer32;
                     R0 : Standard_Natural_Vectors.Vector )
                   return Link_to_Matrix is

  -- DESCRIPTION :
  --   Depending on the value of k, the Alternating_Permute is
  --   called in once, or iteratively.

  -- ON ENTRY :
  --   file      to write diagnostics;
  --   A         second-order derivative of A(1) of mixed type;
  --   nq        number of equations in the original system;
  --   k         label of the child in A(2);
  --   R0        number of variables added in each stage in the deflation.

  -- ON RETURN :
  --   an alternatingly permuted second-order derivative of A(1).

    res : Link_to_Matrix;
    nc,row,col : integer32;

  begin
    res := new Matrix'(A.all);
    row := A'first(1) + nq;
    if k = 1 then
      nc := integer32(R0(0)*R0(0)*R0(1));
      Alternating_Permute(file,res,row,A'first(2),nq,nc,integer32(R0(1)));
    else
      nc := integer32(R0(0)*R0(1));
      col := A'first(2);
      for i in 1..R0(0) loop
        Alternating_Permute(file,res,row,col,nq,nc,integer32(R0(1)));
        col := col + nc;
      end loop;
    end if;
    return res;
  end Permute;

  procedure Multi_Right_Multiply_Deflation
              ( A : in out Link_to_Matrix; m,r,c,k,s : in integer32;
                d,R0 : in Standard_Natural_Vectors.Vector;
                JM : in Link_to_Matrix;
                Bl : in DoblDobl_Complex_Vectors.Link_to_Vector ) is

    JM_offset,Bl_offset,block,row,pivot : integer32;

    procedure Multiply ( index : in Standard_Natural_Vectors.Vector ) is

      zero : constant double_double := create(0.0);
      acc : Complex_Number := Create(zero);
      j1 : constant integer32 := integer32(index(1));
      j2 : constant integer32 := integer32(index(2));
      i : constant integer32 := row;
      A_offset : integer32;

    begin
      for j3 in 1..integer32(R0(block)) loop
        acc := acc + JM(i,JM_offset+j3)*Bl(Bl_offset+j3);
      end loop;
      --if k = 1
      -- then A_offset := (j1-1)*R0(pivot-1) + j2-1;
      -- else A_offset := (j2-1)*R0(pivot-1) + j1-1;
      --end if;
      A_offset := (j1-1)*integer32(R0(pivot-1)) + j2-1;
      if r+i-1 > A'last(1) or c+A_offset > A'last(2)
       then return;
       else A(r+i-1,c+A_offset) := A(r+i-1,c+A_offset) + acc;
      end if;
      JM_offset := JM_offset + integer32(R0(block));
    end Multiply;
    procedure Run_Multiply is new Multi_Loop(Multiply);

  begin
    if d(k) = 0
     then pivot := k;
     else pivot := k+1;
    end if;
    for i in JM'range(1) loop       -- run over all rows of JM
      JM_offset := 0;
      Bl_offset := 0;
      row := i;
      for j0 in 0..m-1 loop         -- m stages: m-column block structure
        block := j0;
        Run_Multiply(natural32(s),d,R0);
        Bl_offset := Bl_offset + integer32(R0(j0));
      end loop;
    end loop;
  end Multi_Right_Multiply_Deflation;

  procedure Multi_Right_Multiply_Deflation
              ( file : in file_type;
                A : in out Link_to_Matrix; m,r,c,k,s : in integer32;
                d,R0 : in Standard_Natural_Vectors.Vector;
                JM : in Link_to_Matrix;
                Bl : in DoblDobl_Complex_Vectors.Link_to_Vector ) is

    JM_offset,Bl_offset,block,row,pivot : integer32;

    procedure Multiply ( index : in Standard_Natural_Vectors.Vector ) is

      zero : constant double_double := create(0.0);
      acc : Complex_Number := Create(zero);
      j1 : constant integer32 := integer32(index(1));
      j2 : constant integer32 := integer32(index(2));
      i : constant integer32 := row;
      A_offset : integer32;

    begin
      for j3 in 1..integer32(R0(block)) loop
        if row = 1 then
          put(file," j1 = "); put(file,j1,1);
          put(file," j2 = "); put(file,j2,1);
          put(file," j3 = "); put(file,j3,1);
          put(file," Bl_off = "); put(file,Bl_offset,1);
          put(file," JM_off = "); put(file,JM_offset,1);
          put(file," Bl_off + j3 = "); put(file,Bl_offset+j3,1);
          put(file," JM_off + j3 = "); put(file,JM_offset+j3,1);
          new_line(file);
        end if;
        acc := acc + JM(i,JM_offset+j3)*Bl(Bl_offset+j3);
      end loop;
      --if k = 1
      -- then A_offset := (j1-1)*R0(pivot-1) + j2-1;
      -- else A_offset := (j2-1)*R0(pivot-1) + j1-1;
      --end if;
      A_offset := (j1-1)*integer32(R0(pivot-1)) + j2-1;
      if row = 1 then
        put(file,"block = "); put(file,block,1);
        put(file,"  A_offset = "); put(file,A_offset,1); new_line(file);
      end if;
      if r+i-1 > A'last(1) or c+A_offset > A'last(2)
       then return;  -- PATCH !!!
       else A(r+i-1,c+A_offset) := A(r+i-1,c+A_offset) + acc;
      end if;
      JM_offset := JM_offset + integer32(R0(block));
    end Multiply;
    procedure Run_Multiply is new Multi_Loop(Multiply);

  begin
    put(file,"-> executing Multi_Right_Multiply_Deflation on Bl, k = ");
    put(file,k,1); put(file,", r = ");
    put(file,r,1); put(file," and c = "); put(file,c,1); new_line(file);
    put(file,"   d = "); put(file,d);
    if d(k) = 0
     then pivot := k;
     else pivot := k+1;
    end if;
    put(file,", pivot = "); put(file,pivot,1);
    put(file,"  JM'last(2) = "); put(file,JM'last(2),1); new_line(file);
    for i in JM'range(1) loop       -- run over all rows of JM
      JM_offset := 0;
      Bl_offset := 0;
      row := i;
      for j0 in 0..m-1 loop         -- m stages: m-column block structure
        block := j0;
        Run_Multiply(natural32(s),d,R0);
        Bl_offset := Bl_offset + integer32(R0(j0));
      end loop;
    end loop;
  end Multi_Right_Multiply_Deflation;

  procedure One_Right_Multiply_Deflation
              ( A : in out Link_to_Matrix; m,r,c,k,p : in integer32;
                R0 : in Standard_Natural_Vectors.Vector;
                JM,B : in Link_to_Matrix ) is

  -- DESCRIPTION :
  --   Applies the code as in Right_Multiply_Deflation from above
  --   to the p-th column of B.

    JM_offset,B_offset : integer32;
    acc : Complex_Number;
    zero : constant double_double := create(0.0);

  begin
    for i in JM'range(1) loop       -- run over all rows of JM
      JM_offset := 0;
      B_offset := 0;
      for j0 in 0..m-1 loop         -- m stages: m-column block structure
        for j1 in 1..integer32(R0(k-1)) loop
          acc := Create(zero);
          for j2 in 1..integer32(R0(j0)) loop
            if JM_offset+j2 > JM'last(2) or B_offset+j2 > B'last(1)
             then return; -- PATCH!!!
             else acc := acc + JM(i,JM_offset+j2)*B(B_offset+j2,p);
            end if;
          end loop;
          if r+i-1 > A'last(1) or c+j1-1 > A'last(2)
           then return; -- PATCH !!!
           else A(r+i-1,c+j1-1) := A(r+i-1,c+j1-1) + acc;
          end if;
          JM_offset := JM_offset + integer32(R0(j0));
        end loop;
        B_offset := B_offset + integer32(R0(j0));
      end loop;
    end loop;
  end One_Right_Multiply_Deflation;

  procedure One_Right_Multiply_Deflation
              ( file : in file_type;
                A : in out Link_to_Matrix; m,r,c,k,p : in integer32;
                R0 : in Standard_Natural_Vectors.Vector;
                JM,B : in Link_to_Matrix ) is

  -- DESCRIPTION :
  --   Applies the code as in Right_Multiply_Deflation from above
  --   to the p-th column of B.

    JM_offset,B_offset : integer32;
    acc : Complex_Number;
    zero : constant double_double := create(0.0);

  begin
    put(file,"-> executing One_Right_Multiply_Deflation on B, p = ");
    put(file,p,1); put(file,", k = "); put(file,k,1);
    put(file,", r = "); put(file,r,1); 
    put(file," and c = "); put(file,c,1); new_line(file);
    for i in JM'range(1) loop       -- run over all rows of JM
      JM_offset := 0;
      B_offset := 0;
      for j0 in 0..m-1 loop         -- m stages: m-column block structure
        for j1 in 1..integer32(R0(k-1)) loop
          acc := Create(zero);
          for j2 in 1..integer32(R0(j0)) loop
            if JM_offset+j2 > JM'last(2) or B_offset+j2 > B'last(1)
             then return; -- PATCH !!!
             else acc := acc + JM(i,JM_offset+j2)*B(B_offset+j2,p);
            end if;
          end loop;
          if r+i-1 > A'last(1) or c+j1-1 > A'last(2)
           then return; -- PATCH !!!
           else A(r+i-1,c+j1-1) := A(r+i-1,c+j1-1) + acc;
          end if;
          JM_offset := JM_offset + integer32(R0(j0));
        end loop;
        B_offset := B_offset + integer32(R0(j0));
      end loop;
    end loop;
  end One_Right_Multiply_Deflation;

  procedure One_Right_Multiply_Deflation
              ( A : in out Link_to_Matrix; m,r,c,k : in integer32;
                R0 : in Standard_Natural_Vectors.Vector;
                JM,B : in Link_to_Matrix ) is

    offset : integer32 := 0;
    nc : constant integer32 := B'last(2)*integer32(R0(k-1));

  begin
    for j in B'range(2) loop
      One_Right_Multiply_Deflation(A,m,r,c+offset,k,j,R0,JM,B);
      offset := offset + integer32(R0(k-1));
    end loop;
    Alternating_Permute(A,r,c,JM'last(1),nc,integer32(R0(k-1)));
  end One_Right_Multiply_Deflation;

  procedure One_Right_Multiply_Deflation
              ( file : in file_type;
                A : in out Link_to_Matrix; m,r,c,k : in integer32;
                R0 : in Standard_Natural_Vectors.Vector;
                JM,B : in Link_to_Matrix ) is

    offset : integer32 := 0;
    nc : constant integer32 := B'last(2)*integer32(R0(k-1));

  begin
    for j in B'range(2) loop
      One_Right_Multiply_Deflation(file,A,m,r,c+offset,k,j,R0,JM,B);
      offset := offset + integer32(R0(k-1));
    end loop;
    Alternating_Permute(file,A,r,c,JM'last(1),nc,integer32(R0(k-1)));
  end One_Right_Multiply_Deflation;

  function Offset_for_One_Right_Multiply
              ( R0 : Standard_Natural_Vectors.Vector; k : integer32 )
              return integer32 is

  -- DESCRIPTION :
  --   Returns the offset to add to the column number in case k > 1
  --   to assign children when right multiplying with a vector.
  --   To compute the offset for a matrix with N columns,
  --   multiply the outcome of this function with N.
  --   Assumes the first order derivatives are used.

    res : integer32;

  begin
    if k = 1 then
      res := 0;
    else
      res := integer32(R0(0));
      for i in 1..k-2 loop
        res := res + integer32(R0(i));
      end loop;
    end if;
    return res;
  end Offset_for_One_Right_Multiply;

  procedure One_Right_Multiply
              ( A : in out Link_to_Matrix; m,r,c,k : in integer32;
                R0 : in Standard_Natural_Vectors.Vector;
                JM : in Link_to_Matrix;
                Bl : in DoblDobl_Complex_Vectors.Link_to_Vector ) is

    A_offset : integer32;

  begin
    if k = 1
     then One_Right_Multiply_Deflation(A,m,r,c,k,R0,JM,Bl);
     else A_offset := Offset_for_One_Right_Multiply(R0,k);
	  One_Right_Multiply_Deflation(A,m,r,c+A_offset,k,R0,JM,Bl);
    end if;
  end One_Right_Multiply;

  procedure One_Right_Multiply
              ( file : in file_type;
                A : in out Link_to_Matrix; m,r,c,k : in integer32;
                R0 : in Standard_Natural_Vectors.Vector;
                JM : in Link_to_Matrix;
                Bl : in DoblDobl_Complex_Vectors.Link_to_Vector ) is

    A_offset : integer32;

  begin
    if k = 1
     then One_Right_Multiply_Deflation(file,A,m,r,c,k,R0,JM,Bl);
     else A_offset := Offset_for_One_Right_Multiply(R0,k);
	  One_Right_Multiply_Deflation(file,A,m,r,c+A_offset,k,R0,JM,Bl);
    end if;
  end One_Right_Multiply;
 
  procedure One_Right_Multiply
              ( A : in out Link_to_Matrix; m,r,c,k : in integer32;
                R0 : in Standard_Natural_Vectors.Vector;
                JM,B : in Link_to_Matrix ) is

    A_offset : integer32;

  begin
    if k = 1
     then One_Right_Multiply_Deflation(A,m,r,c,k,R0,JM,B);
     else A_offset := Offset_for_One_Right_Multiply(R0,k)*B'last(2);
          One_Right_Multiply_Deflation(A,m,r,c+A_offset,k,R0,JM,B);
    end if;
  end One_Right_Multiply;
 
  procedure One_Right_Multiply
              ( file : in file_type;
                A : in out Link_to_Matrix; m,r,c,k : in integer32;
                R0 : in Standard_Natural_Vectors.Vector;
                JM,B : in Link_to_Matrix ) is

    A_offset : integer32;

  begin
    if k = 1
     then One_Right_Multiply_Deflation(file,A,m,r,c,k,R0,JM,B);
     else A_offset := Offset_for_One_Right_Multiply(R0,k)*B'last(2);
          One_Right_Multiply_Deflation(file,A,m,r,c+A_offset,k,R0,JM,B);
    end if;
  end One_Right_Multiply;

-- ASSIGNMENTS :

  procedure Assign_Scaling_Coefficients
               ( A : in out Link_to_Matrix;
                 h : in DoblDobl_Complex_Vectors.Link_to_Vector ) is

    len : constant integer32 := h'length;
    ind : integer32 := A'last(2)-len;

  begin
    for i in h'range loop
      ind := ind + 1;
      A(A'last(1),ind) := h(i);
    end loop;
  end Assign_Scaling_Coefficients;

  procedure Assign_Lower_Jacobian_Matrices
               ( A : in out Link_to_Matrix; row,col : in integer32;
                 jms : in Link_to_VecMat;
                 Bl : in DoblDobl_Complex_Vectors.Link_to_Vector ) is

    wrkcol : integer32 := col;
    JM : Link_to_Matrix;

  begin
    for k in jms'range loop
      if jms(k) /= null then
        JM := jms(k);        
        Multiply(A,row,wrkcol,JM,Bl);
      end if;
      wrkcol := wrkcol + 1;
    end loop;
  end Assign_Lower_Jacobian_Matrices;

  procedure Assign_Lower_Jacobian_Matrices
               ( A : in out Link_to_Matrix; row,col : in integer32;
                 jms : in Link_to_VecMat; B : in Link_to_Matrix ) is

    wrkcol : integer32 := col;
    JM : Link_to_Matrix;

  begin
    for k in jms'range loop
      if jms(k) /= null then
        JM := jms(k);
        Multiply(A,row,wrkcol,JM,B);
      end if;
      wrkcol := wrkcol + B'last(2);
    end loop;
  end Assign_Lower_Jacobian_Matrices;

  procedure Assign_from_Jacobian_Matrices
               ( A : in out Link_to_Matrix; col : in out integer32;
                 jms : in Link_to_VecMat; nbcols : in natural32 ) is

    JM : Link_to_Matrix;

  begin
    for k in jms'range loop
      if jms(k) /= null then
        JM := jms(k);
        for i in JM'range(1) loop
          for j in JM'range(2) loop
            A(i,j+col-1) := JM(i,j);
          end loop;
        end loop;
      end if;
      col := col + integer32(nbcols);       -- we skip void matrices
    end loop;
  end Assign_from_Jacobian_Matrices;

  procedure Assign_from_Jacobian_Matrices
               ( A : in out Link_to_Matrix; row : in integer32;
                 col : in out integer32; jms : in Link_to_VecMat;
                 B : in Link_to_Matrix ) is

    JM : Link_to_Matrix;

  begin
    for k in jms'range loop
      if jms(k) /= null then
        JM := jms(k);
        Multiply(A,row,col,JM,B);
      end if;
      col := col + B'last(2);
    end loop;
  end Assign_from_Jacobian_Matrices;

  procedure Assign_Higher_Jacobian_Matrices
               ( A : in out Link_to_Matrix; row : in integer32;
                 col : in out integer32; jms : in Link_to_VecMat;
                 monkeys : in Standard_Natural64_VecVecs.VecVec;
                 k,n : in integer32 ) is

    procedure Assign ( m : in Standard_Natural_Vectors.Vector ) is

      ind : constant integer32 := integer32(Search(monkeys,m));
      JM : constant Link_to_Matrix := jms(ind);

    begin
      if JM /= null then
        for i in JM'range(1) loop
          for j in JM'range(2) loop
            A(row+i-1,col+j-1) := JM(i,j);
          end loop;
        end loop;
      end if;
      col := col + n;
    end Assign;
    procedure Enum is new Enumerate_Leaves_of_Monomial_Tree(Assign);

  begin
    Enum(natural32(k),natural32(n));
  end Assign_Higher_Jacobian_Matrices;

  procedure Assign_Higher_Jacobian_Matrices
               ( file : in file_type;
                 A : in out Link_to_Matrix; row : in integer32;
                 col : in out integer32; jms : in Link_to_VecMat;
                 monkeys : in Standard_Natural64_VecVecs.VecVec;
                 k,n : in integer32 ) is

    procedure Assign ( m : in Standard_Natural_Vectors.Vector ) is

      ind : constant integer32 := integer32(Search(monkeys,m));
      JM : constant Link_to_Matrix := jms(ind);

    begin
      if JM /= null then
        put(file,"assigning Jacobian "); put(file,m);
        put(file," into row "); put(file,row,1);
        put(file," and column "); put(file,col,1); new_line(file);
        for i in JM'range(1) loop
          for j in JM'range(2) loop
            A(row+i-1,col+j-1) := JM(i,j);
          end loop;
        end loop;
      end if;
      col := col + n;
    end Assign;
    procedure Enum is new Enumerate_Leaves_of_Monomial_Tree(Assign);

  begin
    Enum(natural32(k),natural32(n));
  end Assign_Higher_Jacobian_Matrices;

  procedure Assign_Higher_Jacobian_Matrices
               ( A : in out Link_to_Matrix; row : in integer32;
                 col : in out integer32; jms : in Link_to_VecMat;
                 monkeys : in Standard_Natural64_VecVecs.VecVec;
                 Bl : in DoblDobl_Complex_Vectors.Link_to_Vector;
                 k,n : in integer32 ) is

    procedure Assign ( m : in Standard_Natural_Vectors.Vector ) is

      ind : constant integer32 := integer32(Search(monkeys,m));
      JM : constant Link_to_Matrix := jms(ind);

    begin
      if JM /= null then
        Multiply(A,row,col,JM,Bl);
      end if;
      col := col+1;
    end Assign;
    procedure Enum is new Enumerate_Leaves_of_Monomial_Tree(Assign);

  begin
    Enum(natural32(k),natural32(n));
  end Assign_Higher_Jacobian_Matrices;

  procedure Assign_Higher_Jacobian_Matrices
               ( file : in file_type;
                 A : in out Link_to_Matrix; row : in integer32;
                 col : in out integer32; jms : in Link_to_VecMat;
                 monkeys : in Standard_Natural64_VecVecs.VecVec;
                 Bl : in DoblDobl_Complex_Vectors.Link_to_Vector;
                 k,n : in integer32 ) is

    procedure Assign ( m : in Standard_Natural_Vectors.Vector ) is

      ind : constant integer32 := integer32(Search(monkeys,m));
      JM : constant Link_to_Matrix := jms(ind);

    begin
      put(file,"multiplying Jacobian matrix "); put(file,m);
      put(file," with Bl into row = "); put(file,row,1);
      put(file," and column = "); put(file,col,1); put(file," ... ");
      if JM /= null then
        Multiply(A,row,col,JM,Bl);
        put_line(file," done.");
      else
        put_line(file," empty Jacobian.");
      end if;
      col := col + 1;
    end Assign;
    procedure Enum is new Enumerate_Leaves_of_Monomial_Tree(Assign);

  begin
    Enum(natural32(k),natural32(n));
  end Assign_Higher_Jacobian_Matrices;

  procedure Assign_Higher_Jacobian_Matrices
               ( A : in out Link_to_Matrix; row : in integer32;
                 col : in out integer32; jms : in Link_to_VecMat;
                 monkeys : in Standard_Natural64_VecVecs.VecVec;
                 B : in Link_to_Matrix; k,n : in integer32 ) is

    procedure Assign ( m : in Standard_Natural_Vectors.Vector ) is

      ind : constant integer32 := integer32(Search(monkeys,m));
      JM : constant Link_to_Matrix := jms(ind);

    begin
      if JM /= null then
        Multiply(A,row,col,JM,B);
      end if;
      col := col + B'last(2);
    end Assign;
    procedure Enum is new Enumerate_Leaves_of_Monomial_Tree(Assign);

  begin
    Enum(natural32(k),natural32(n));
  end Assign_Higher_Jacobian_Matrices;

  procedure Assign_Higher_Jacobian_Matrices
               ( file : in file_type;
                 A : in out Link_to_Matrix; row : in integer32;
                 col : in out integer32; jms : in Link_to_VecMat;
                 monkeys : in Standard_Natural64_VecVecs.VecVec;
                 B : in Link_to_Matrix; k,n : in integer32 ) is

    procedure Assign ( m : in Standard_Natural_Vectors.Vector ) is

      ind : constant integer32 := integer32(Search(monkeys,m));
      JM : constant Link_to_Matrix := jms(ind);

    begin
      if JM /= null then
        put(file,"multiplying JM "); put(file,m);
        put(file,"(ind = "); put(file,ind,1);
        put(file,") with B in row = "); put(file,row,1); 
        put(file," and column = "); put(file,col,1); new_line(file);
        Multiply(A,row,col,JM,B);
      end if;
      col := col + B'last(2);
    end Assign;
    procedure Enum is new Enumerate_Leaves_of_Monomial_Tree(Assign);

  begin
    Enum(natural32(k),natural32(n));
  end Assign_Higher_Jacobian_Matrices;

  procedure Assign_Root_Child
               ( A : in out Link_to_Matrix; d0,m : in natural32;
                 R0 : in Standard_Natural_Vectors.Vector;
                 v0 : in Link_to_Matrix;
                 monkeys : in Standard_Natural64_VecVecs.VecVec;
                 jrt : in Jacobian_Remember_Table; B : in Link_to_Matrix ) is

  -- NOTE : There are two cases to consider, depending on the value of v0
  --   (1) v0 is an index to the Jacobian remember table;
  --   (2) v0 is an evaluated deflation matrix.

    row,col,ind,nbrows,nbcols : integer32;

  begin
    if ((v0'first(1) = v0'last(1)) and
        (v0'first(2) = v0'last(2))) then  -- case (1)
      ind := integer32(to_double(REAL_PART(v0(1,1))));
      col := A'first(2);
      Dimensions(jrt(ind),nbrows,nbcols);
      if nbrows > 0 then  -- otherwise everything is zero
        if ind <= 1 then
          Assign_from_Jacobian_Matrices(A,col,jrt(ind),natural32(nbcols));
          row := A'first(1) + nbrows;       -- start row to place jrt(ind)*B
          Assign_from_Jacobian_Matrices(A,row,col,jrt(ind),B);
        else
          row := A'first(1);
          Assign_Higher_Jacobian_Matrices
            (A,row,col,jrt(ind),monkeys,ind,nbcols);
          row := A'first(1) + nbrows;
          Assign_Higher_Jacobian_Matrices
            (A,row,col,jrt(ind),monkeys,B,ind,nbcols);
        end if;
      end if;
    else                                  -- case(2)
      for i in v0'range(1) loop           -- copy v0 into A
        for j in v0'range(2) loop
          A(i,j) := v0(i,j);
        end loop;
      end loop;
      row := A'first(1)+v0'last(1);       -- start row to place v0*B
      col := A'first(2)+v0'last(2);       -- start column for v0*B
      if v0'last(2) = B'last(1) then      -- assign v0*B into A
        Multiply(A,row,col,v0,B);         -- regular matrix multiplication
      elsif d0 > 0 then
        One_Right_Multiply_Deflation(A,integer32(m),row,col,1,R0,v0,B);
      else
        One_Right_Multiply_Deflation(A,integer32(m),row,col,
                                     integer32(m),R0,v0,B);
      end if;
    end if;
  end Assign_Root_Child;

  procedure Assign_Root_Child
               ( file : in file_type;
                 A : in out Link_to_Matrix; d0,m : in natural32;
                 R0 : in Standard_Natural_Vectors.Vector;
                 v0 : in Link_to_Matrix;
                 monkeys : in Standard_Natural64_VecVecs.VecVec;
                 jrt : in Jacobian_Remember_Table; B : in Link_to_Matrix ) is

  -- NOTE : There are two cases to consider, depending on the value of v0
  --   (1) v0 is an index to the Jacobian remember table;
  --   (2) v0 is an evaluated deflation matrix.

    row,col,ind,nbrows,nbcols : integer32;

  begin
    if ((v0'first(1) = v0'last(1)) and
        (v0'first(2) = v0'last(2))) then  -- case (1)
      ind := integer32(to_double(REAL_PART(v0(1,1))));
      col := A'first(2);
      Dimensions(jrt(ind),nbrows,nbcols);
      put(file,"Dimensions ind = "); put(file,ind,1);
      put(file,"  nbcols = "); put(file,nbcols,1); new_line(file);
      if nbrows > 0 then  -- otherwise everything is zero
        if ind <= 1 then
          Assign_from_Jacobian_Matrices(A,col,jrt(ind),natural32(nbcols));
          row := A'first(1) + nbrows;       -- start row to place jrt(ind)*B
          Assign_from_Jacobian_Matrices(A,row,col,jrt(ind),B);
        else
          row := A'first(1);
          put(file,"Higher 1 with row = "); put(file,row,1); 
          put(file," and col = "); put(file,col,1); new_line(file);
          Assign_Higher_Jacobian_Matrices
            (file,A,row,col,jrt(ind),monkeys,ind,nbcols);
          row := A'first(1) + nbrows;
          put(file,"Higher 2 with row = "); put(file,row,1);
          put(file," and col = "); put(file,col,1); new_line(file);
          Assign_Higher_Jacobian_Matrices
            (file,A,row,col,jrt(ind),monkeys,B,ind,nbcols);
        end if;
      end if;
    else                                  -- case(2)
      put(file,"Assign_Root_Child case(2), m = "); put(file,m,1);
      put(file," v0'last(2) = "); put(file,v0'last(2),1);
      put(file," B'last(1) = ");  put(file,B'last(1),1);
      put(file," d0 = "); put(file,d0,1); new_line(file);
      for i in v0'range(1) loop           -- copy v0 into A
        for j in v0'range(2) loop
          A(i,j) := v0(i,j);
        end loop;
      end loop;
      row := A'first(1)+v0'last(1);       -- start row to place v0*B
      col := A'first(2)+v0'last(2);       -- start column for v0*B
      if v0'last(2) = B'last(1) then      -- assign v0*B into A
        Multiply(A,row,col,v0,B);         -- regular matrix multiplication
      elsif d0 > 0 then
        One_Right_Multiply_Deflation(file,A,integer32(m),row,col,1,R0,v0,B);
      else
        One_Right_Multiply_Deflation(file,A,integer32(m),row,col,
                                     integer32(m),R0,v0,B);
      end if;
    end if;
  end Assign_Root_Child;

  procedure Assign_Jacobian_Child
               ( A : in out Link_to_Matrix; m : in integer32;
                 R0 : in Standard_Natural_Vectors.Vector;
                 v0 : in Link_to_Matrix;
                 monkeys : in Standard_Natural64_VecVecs.VecVec;
                 jrt : in Jacobian_Remember_Table; B : in Link_to_Matrix;
                 Bl : in DoblDobl_Complex_Vectors.Link_to_Vector;
                 ind : in integer32; col : in out integer32 ) is

  -- DESCRIPTION :
  --   Assigns into A a Jacobian matrix of the original system.

    row,nbrows,nbcols,nc : integer32;

  begin
    Dimensions(jrt(ind),nbrows,nbcols);
    if nbrows > 0 then  -- otherwise there is nothing to do
      if ind <= 1 then
        row := A'first(1) + nbrows;
        if v0 = null then
          Assign_Lower_Jacobian_Matrices(A,row,col,jrt(ind),B);
          nc := jrt(ind)'length*B'last(2);
          Alternating_Permute(A,row,col,nbrows,nc,integer32(R0(m)));
        else
          Assign_Lower_Jacobian_Matrices(A,row,col,jrt(ind),Bl);
          nc := jrt(ind)'length;
        end if;
        col := col + nc;
      else
        row := nbrows+1;
        if v0 = null then
          Assign_Higher_Jacobian_Matrices
            (A,row,col,jrt(ind),monkeys,B,ind,nbcols);
        else
          Assign_Higher_Jacobian_Matrices
            (A,row,col,jrt(ind),monkeys,Bl,ind,nbcols);
          nc := integer32(R0(0))**integer(ind);
        end if;
      end if;
    end if;
  end Assign_Jacobian_Child;

  procedure Assign_Jacobian_Child
               ( file : in file_type;
                 A : in out Link_to_Matrix; m : in integer32; 
                 R0 : in Standard_Natural_Vectors.Vector;
                 v0 : in Link_to_Matrix;
                 monkeys : in Standard_Natural64_VecVecs.VecVec;
                 jrt : in Jacobian_Remember_Table; B : in Link_to_Matrix;
                 Bl : in DoblDobl_Complex_Vectors.Link_to_Vector;
                 ind : in integer32; col : in out integer32 ) is

  -- DESCRIPTION :
  --   Assigns into A a Jacobian matrix of the original system.

    row,nbrows,nbcols,nc : integer32;

  begin
    Dimensions(jrt(ind),nbrows,nbcols);
    put(file," nbrows = "); put(file,nbrows,1);
    put(file," and nbcols = "); put(file,nbcols,1); new_line(file);
    if nbrows > 0 then  -- otherwise there is nothing to do
      if ind <= 1 then
        row := A'first(1) + nbrows;
        if v0 = null then
          put_line(file,"Calling first lower ...");
          Assign_Lower_Jacobian_Matrices(A,row,col,jrt(ind),B);
          nc := jrt(ind)'length*B'last(2);
          Alternating_Permute(file,A,row,col,nbrows,nc,integer32(R0(m)));
        else
          put_line(file,"Calling second lower ...");
          Assign_Lower_Jacobian_Matrices(A,row,col,jrt(ind),Bl);
          nc := jrt(ind)'length;
        end if;
        col := col + nc;
      else
        row := nbrows+1;
        if v0 = null then
          put(file," higher 1 with row = "); put(file,row,1);
          put(file," and col = "); put(file,col,1); new_line(file);
          Assign_Higher_Jacobian_Matrices
            (file,A,row,col,jrt(ind),monkeys,B,ind,nbcols);
         -- nc := (R0(0)**ind)*B'last(2);
         -- Alternating_Permute(A,row,backup_col,nbrows,nc,R0(m));
         -- code below only works for ind = 2, need multi-loop in general
         -- nc := R0(0)*B'last(2);
         -- col := backup_col;     
         -- for i in 1..R0(0) loop
         --   Alternating_Permute(A,row,col,nbrows,nc,R0(m));
         --   col := col + nc;
         -- end loop;
        else
          put(file," higher 2 with row = "); put(file,row,1);
          put(file," and col = "); put(file,col,1); new_line(file);
          Assign_Higher_Jacobian_Matrices
            (file,A,row,col,jrt(ind),monkeys,Bl,ind,nbcols);
          nc := integer32(R0(0))**integer(ind);
        end if;
      end if;
    end if;
  end Assign_Jacobian_Child;

  procedure Assign_Children
               ( A : in out Link_to_Matrix; m : in integer32;
                 R0 : in Standard_Natural_Vectors.Vector;
                 v : in VecMat; B : in Link_to_Matrix ) is

    row,col,nc,A_offset : integer32;
    JM : Link_to_Matrix;

  begin
    col := A'first(2);
    for k in 1..v'last loop
      if v(k) /= null then
        JM := v(k);
        row := A'first(1) + JM'last(1);
        One_Right_Multiply(A,m,row,col,k,R0,JM,B);
        A_offset := Offset_for_One_Right_Multiply(R0,k)*B'last(2);
        nc := integer32(R0(k-1))*B'last(2);
        Alternating_Permute(A,row,col+A_offset,JM'last(1),nc,integer32(R0(m)));
      end if;
    end loop;
  end Assign_Children;

  procedure Assign_Children
               ( file : in file_type;
                 A : in out Link_to_Matrix; m : in integer32;
                 R0 : in Standard_Natural_Vectors.Vector;
                 v : in VecMat; B : in Link_to_Matrix ) is

    row,col,nc,A_offset : integer32;
    JM : Link_to_Matrix;

  begin
    col := A'first(2);
    for k in 1..v'last loop
      if v(k) /= null then
        JM := v(k);
        row := A'first(1) + JM'last(1);
        One_Right_Multiply(file,A,m,row,col,k,R0,JM,B);
        A_offset := Offset_for_One_Right_Multiply(R0,k)*B'last(2);
        nc := integer32(R0(k-1))*B'last(2);
        Alternating_Permute
          (file,A,row,col+A_offset,JM'last(1),nc,integer32(R0(m)));
      end if;
    end loop;
  end Assign_Children;

  procedure Assign_Children
               ( A : in out Link_to_Matrix; nq,d0,m,order : in natural32; 
                 chtp : in Standard_Natural_VecVecs.VecVec;
                 R0 : in Standard_Natural_Vectors.Vector; v : in VecMat;
                 monkeys : in Standard_Natural64_VecVecs.VecVec;
                 jrt : in Jacobian_Remember_Table; B : in Link_to_Matrix;
                 Bl : in DoblDobl_Complex_Vectors.Link_to_Vector ) is

    AA : Link_to_Matrix;
    row,col,ind : integer32;

    use Standard_Natural_Vectors;

  begin
    if v(0) /= null then             -- assign v(0) into A
      Assign_Root_Child(A,d0,m,R0,v(0),monkeys,jrt,B);
    end if;
    col := A'first(2);               -- start column for v(1..m)*Bl
    for k in 1..v'last loop          -- compute v(k)*Bl
      if v(k) /= null then
        AA := v(k);
        if AA'first(1) = AA'last(1) and AA'first(2) = AA'last(2) then
          ind := integer32(to_double(REAL_PART(AA(AA'first(1),AA'first(2)))));
          Assign_Jacobian_Child(A,integer32(m),R0,v(0),monkeys,
                                jrt,B,Bl,ind,col);
        else
          row := A'first(1) + AA'last(1);    -- start row to place v(k)*Bl
          if order = 1 then
            One_Right_Multiply(A,integer32(m),row,A'first(2),k,R0,AA,Bl);
          else
            if chtp(k) /= null then
              if ((R0(1) > 1) and then Is_Mixed(chtp(k).all)) then
                declare
                  pAA : Link_to_Matrix := Permute(AA,integer32(nq),k,R0);
                begin
                  Multi_Right_Multiply_Deflation
                    (A,integer32(m),row,col,k,integer32(order),
                     chtp(k).all,R0,pAA,Bl);
                  Clear(pAA);
                end;
              else
                Multi_Right_Multiply_Deflation
                  (A,integer32(m),row,col,k,integer32(order),
                   chtp(k).all,R0,AA,Bl);
              end if;
            end if;
            col := col + integer32(R0(k-1)*R0(k-1));
          end if;
        end if;
      end if;
    end loop;
  end Assign_Children;

  procedure Assign_Children
               ( file : in file_type;
                 A : in out Link_to_Matrix; nq,d0,m,order : in natural32; 
                 chtp : in Standard_Natural_VecVecs.VecVec;
                 R0 : in Standard_Natural_Vectors.Vector; v : in VecMat;
                 monkeys : in Standard_Natural64_VecVecs.VecVec;
                 jrt : in Jacobian_Remember_Table; B : in Link_to_Matrix;
                 Bl : in DoblDobl_Complex_Vectors.Link_to_Vector ) is

    AA : Link_to_Matrix;
    row,col,ind : integer32;

    use Standard_Natural_Vectors;

  begin
    if v(0) /= null then             -- assign v(0) into A
      Assign_Root_Child(file,A,d0,m,R0,v(0),monkeys,jrt,B);
    end if;
    col := A'first(2);               -- start column for v(1..m)*Bl
    for k in 1..v'last loop          -- compute v(k)*Bl
      if v(k) /= null then
        AA := v(k);
        if AA'first(1) = AA'last(1) and AA'first(2) = AA'last(2) then
          ind := integer32(to_double(REAL_PART(AA(AA'first(1),AA'first(2)))));
          put(file,"  ind = "); put(file,ind,1);
          put(file,"  jrt(ind)'length = ");
          put(file,natural32(jrt(ind)'length),1);
          Assign_Jacobian_Child(file,A,integer32(m),R0,v(0),monkeys,
                                jrt,B,Bl,ind,col);
        else
          put(file," child "); put(file,k,1); put(file," is AA :");
         -- Write_Matrix(AA.all);
          row := A'first(1) + AA'last(1);    -- start row to place v(k)*Bl
          put(file," AA'last(1) = "); put(file,AA'last(1),1);
          put(file," AA'last(2) = "); put(file,AA'last(2),1);
          put(file," row = "); put(file,row,1);
          put(file," col = "); put(file,col,1); new_line(file);
          if order = 1 then
            put(file,"One_Right_Multiply with Bl in Assign_Children");
            put(file," with m = "); put(file,m,1); put_line(file,"...");
            One_Right_Multiply(file,A,integer32(m),row,A'first(2),k,R0,AA,Bl);
          else
            put(file,"Multi_Right_Multiply with Bl in Assign_Children");
            put(file," with m = ");
            put(file,m,1); put(file," type"); put(file,chtp(k));
            if chtp(k) /= null then
              if ((R0(1) > 1) and then Is_Mixed(chtp(k).all)) then
                put_line(file," mixed");
                declare
                  pAA : Link_to_Matrix := Permute(file,AA,integer32(nq),k,R0);
                begin
                  Multi_Right_Multiply_Deflation
                    (file,A,integer32(m),row,col,k,integer32(order),
                     chtp(k).all,R0,pAA,Bl);
                  Clear(pAA);
                end;
              else
                put_line(file," pure");
                Multi_Right_Multiply_Deflation
                  (file,A,integer32(m),row,col,k,integer32(order),
                   chtp(k).all,R0,AA,Bl);
              end if;
            end if;
            col := col + integer32(R0(k-1)*R0(k-1));
          end if;
        end if;
      end if;
    end loop;
  end Assign_Children;

end DoblDobl_Deflation_Matrices;
