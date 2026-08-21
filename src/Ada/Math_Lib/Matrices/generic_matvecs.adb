with unchecked_deallocation;

package body Generic_MatVecs is

  procedure Copy ( A : in MatVec; B : in out MatVec ) is
  begin
    Clear(B);
    for i in A'range(1) loop
      for j in B'range(1) loop
        declare
          vec : constant Vectors.Vector := A(i,j).all;
        begin
          B(i,j) := new Vectors.Vector'(vec);
        end;
      end loop; 
    end loop; 
  end Copy;

  procedure Clear ( A : in out MatVec ) is
  begin
    for i in A'range(1) loop
      for j in A'range(2) loop
        Vectors.Clear(A(i,j));
      end loop;
    end loop;
  end Clear;

  procedure Shallow_Clear ( A : in out Link_to_MatVec ) is

    procedure free is new unchecked_deallocation(MatVec,Link_to_MatVec);

  begin
    free(A);
  end Shallow_Clear;

  procedure Deep_Clear ( A : in out Link_to_MatVec ) is
  begin
    if A /= null
     then Clear(A.all); Shallow_Clear(A);
    end if;
  end Deep_Clear;

  procedure Clear ( A : in out MatVec_Array ) is
  begin
    for k in A'range loop
      Deep_Clear(A(k));
    end loop;
  end Clear;

end Generic_MatVecs;
