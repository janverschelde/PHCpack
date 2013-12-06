with unchecked_deallocation;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;

package body Matrix_Homotopies is

-- INTERNAL DATA STRUCTURES :

  type Matrix_Homotopy ( n,m : integer32 ) is record
    start,target : Matrix(1..n,1..m);
  end record;
  type Link_to_Matrix_Homotopy is access Matrix_Homotopy;

  type Matrix_Homotopy_Array is
    array ( integer32 range <> ) of Link_to_Matrix_Homotopy;
  type Link_to_Matrix_Homotopy_Array is access Matrix_Homotopy_Array;

-- INTERNAL DATA :

  mathom : Link_to_Matrix_Homotopy_Array;
  curmat : integer32 := 0;

-- CREATORS :

  procedure Init ( n : in natural32 ) is
  begin
    mathom := new Matrix_Homotopy_Array(1..integer32(n));
    curmat := 0;
  end Init;

  procedure Add ( start,target : in Matrix ) is
  begin
    curmat := curmat+1;
    mathom(curmat) := new Matrix_Homotopy(start'last(1),start'last(2));
    mathom(curmat).start := start;
    mathom(curmat).target := target;
  end Add;

  procedure Add_Start ( mapno : in natural32; start : in Matrix ) is
  begin
    if mathom(integer32(mapno)) = null then
      mathom(integer32(mapno))
        := new Matrix_Homotopy(start'last(1),start'last(2));
      curmat := integer32(mapno);
    end if;
    mathom(integer32(mapno)).start := start;
  end Add_Start;

  procedure Add_Target ( mapno : in natural32; target : in Matrix ) is
  begin
    if mathom(integer32(mapno)) = null then
      mathom(integer32(mapno))
        := new Matrix_Homotopy(target'last(1),target'last(2));
      curmat := integer32(mapno);
    end if;
    mathom(integer32(mapno)).target := target;
  end Add_Target;

-- SELECTORS :

  function Empty ( mapno : natural32 ) return boolean is
  begin
    return (mathom(integer32(mapno)) = null);
  end Empty;

  function Cardinality return natural32 is
  begin
    return natural32(curmat);
  end Cardinality;

-- EVALUATOR :

  function Eval ( mapno : natural32; t : Complex_Number ) return Matrix is

    mho : constant Link_to_Matrix_Homotopy := mathom(integer32(mapno));
    res : Matrix(1..mho.n,1..mho.m);
    m1t : constant Complex_Number := Create(1.0) - t;

  begin
    if t = Create(0.0) then
      res := mho.start;
    elsif t = Create(1.0) then
      res := mho.target;
    else
      for i in res'range(1) loop
        for j in res'range(2) loop
          res(i,j) := m1t*mho.start(i,j) + t*mho.target(i,j);
        end loop;
      end loop;
    end if;
    return res;
  end Eval;

-- DESTRUCTOR :

  procedure Clear ( mh : in out Link_to_Matrix_Homotopy ) is

    procedure free is
      new unchecked_deallocation(Matrix_Homotopy,Link_to_Matrix_Homotopy); 

  begin
    free(mh);
  end Clear;

  procedure Clear is

    procedure free is
      new unchecked_deallocation(Matrix_Homotopy_Array,
                                 Link_to_Matrix_Homotopy_Array); 

  begin
    for i in 1..curmat loop
      Clear(mathom(i));
    end loop;
    free(mathom);
  end Clear;

end Matrix_Homotopies;
