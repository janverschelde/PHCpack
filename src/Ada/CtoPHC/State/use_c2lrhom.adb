with text_io;                           use text_io;
with Interfaces.C;                      use Interfaces.C;
with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
-- with Standard_Natural_Numbers_io;       use Standard_Natural_Numbers_io;
with Multprec_Natural_Numbers;          use Multprec_Natural_Numbers;
-- with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;         use Standard_Floating_Numbers;
with Brackets;                          use Brackets;
-- with Brackets_io;                       use Brackets_io;
with Drivers_for_Schubert_Induction;    use Drivers_for_Schubert_Induction;
with Assignments_in_Ada_and_C;          use Assignments_in_Ada_and_C;

function use_c2lrhom ( job : integer32;
                       a : C_intarrs.Pointer;
		       b : C_intarrs.Pointer;
                       c : C_dblarrs.Pointer ) return integer32 is

  procedure Get_Dimensions
              ( a : C_intarrs.Pointer; n,k,c : out integer32;
                verbose : out boolean ) is

  -- DESCRIPTION :
  --   Extracts the dimension (n,k,c) from the array a,
  --   where n is the ambient dimension, k the dimension of the solution
  --   planes, c the number of intersection conditions, and verbose flags
  --   whether addition output should be written during the resolution.

    v : constant C_Integer_Array(0..3)
      := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(4));

  begin
    n := integer32(v(0));
    k := integer32(v(1));
    c := integer32(v(2));
    if integer32(v(3)) = 1
     then verbose := true;
     else verbose := false;
    end if;
  end Get_Dimensions;

  function Get_Conditions
             ( b : C_intarrs.Pointer; k,nb : integer32 )
             return Array_of_Brackets is

  -- DESCRIPTOIN :
  --   Extracts the brackets from the array b,
 
  -- REQUIRED :
  --   The array b must contain k times nb integers.

  -- ON ENTRY :
  --   b      an array of integers;
  --   k      dimension of the solution planes;
  --   nb     the number of conditions. 

    res : Array_of_Brackets(1..nb);
    dim : constant integer32 := k*nb;
    dmc : constant Interfaces.C.size_t := Interfaces.C.size_t(dim);
    dm1 : constant Interfaces.C.size_t := Interfaces.C.size_t(dim-1);
    v : constant C_Integer_Array(0..dm1)
      := C_intarrs.Value(b,Interfaces.C.ptrdiff_t(dmc));
    ind : Interfaces.C.size_t := 0;

  begin
    for i in res'range loop
      declare
        brk : Bracket(1..k);
      begin
        for j in 1..k loop
          brk(j) := natural32(v(ind));
          ind := ind + 1;
        end loop;
        res(i) := new Bracket'(brk);
      end;
    end loop;
    return res;
  end Get_Conditions;

  function Job0 return integer32 is -- resolve Schubert intersection conditions

    n,k,nbc : integer32;
    otp : boolean;
    rc : Natural_Number;
    nrc : natural32;

  begin
    Get_Dimensions(a,n,k,nbc,otp);
   -- new_line;
   -- put_line("The dimensions : ");
   -- put("  n = "); put(n,1);
   -- put("  k = "); put(k,1);
   -- put("  c = "); put(nbc,1);
   -- if otp
   --  then put_line("  output wanted");
   --  else put_line("  in silent mode");
   -- end if;
    declare
      cond : constant Array_of_Brackets(1..nbc) := Get_Conditions(b,k,nbc);
    begin
     -- put_line("The brackets : ");
     -- for i in cond'range loop
     --   put(cond(i).all);
     -- end loop;
     -- new_line;
      Create_Intersection_Poset(n,nbc,cond,not otp,rc);
    end;
    nrc := Multprec_Natural_Numbers.Create(rc);
   -- put("The formal root count : "); put(nrc,1); new_line;
    Assign(double_float(nrc),c);
    return 0;
  end Job0;

  function Handle_Jobs return integer32 is
  begin
    case job is
      when 0 => return Job0; -- resolve Schubert intersection conditions
      when others => put_line("  Sorry.  Invalid operation in use_c2lrhom.");
                     return 1;
    end case;
  end Handle_Jobs;

begin
  return Handle_Jobs;
end use_c2lrhom;
