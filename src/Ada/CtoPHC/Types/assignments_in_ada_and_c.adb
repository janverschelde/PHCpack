with Interfaces.C;                      use Interfaces.C;
with Characters_and_Numbers;            use Characters_and_Numbers;

package body Assignments_in_Ada_and_C is

  procedure Assign ( this : in integer32; to_that : in C_intarrs.Pointer ) is

    val : C_Integer_Array(0..0); -- := C_intarrs.Value(to_that);

  begin
    val(0) := Interfaces.C.int(this);
    C_intarrs.Copy_Array(val(0)'unchecked_access,to_that,1);
  end Assign;

  procedure Assign ( this : in C_intarrs.Pointer; to_that : out integer32 ) is

    val : constant C_Integer_Array := C_intarrs.Value(this);
 
  begin
    to_that := integer32(val(val'first));
  end Assign;

  procedure Assign ( this : in double_float;
                     to_that : in C_dblarrs.Pointer ) is

    val : C_Double_Array(0..0); -- := C_dblarrs.Value(to_that);

  begin
    val(0) := Interfaces.C.double(this);
    C_dblarrs.Copy_Array(val(0)'unchecked_access,to_that,1);
  end Assign;

  procedure Assign ( this : in C_dblarrs.Pointer;
                     to_that : out double_float ) is

    val : constant C_Double_Array := C_dblarrs.Value(this);
 
  begin
    to_that := double_float(val(val'first));
  end Assign;

  procedure Assign ( this : in double_double;
                     to_that : in C_dblarrs.Pointer ) is

    val : C_Double_Array(0..1);

  begin
    val(0) := Interfaces.C.double(hi_part(this));
    val(1) := Interfaces.C.double(lo_part(this));
    C_dblarrs.Copy_Array(val(0)'unchecked_access,to_that,2);
  end Assign;

  procedure Assign ( this : in C_dblarrs.Pointer;
                     to_that : out double_double ) is

    val : constant C_Double_Array := C_dblarrs.Value(this,2);

  begin
    to_that := Create(double_float(val(0)),double_float(val(1)));
  end Assign;

  procedure Assign ( this : in quad_double;
                     to_that : in C_dblarrs.Pointer ) is

    val : C_Double_Array(0..3);

  begin
    val(0) := Interfaces.C.double(hihi_part(this));
    val(1) := Interfaces.C.double(lohi_part(this));
    val(2) := Interfaces.C.double(hilo_part(this));
    val(3) := Interfaces.C.double(lolo_part(this));
    C_dblarrs.Copy_Array(val(0)'unchecked_access,to_that,4);
  end Assign;

  procedure Assign ( this : in C_dblarrs.Pointer;
                     to_that : out quad_double ) is

    val : constant C_Double_Array := C_dblarrs.Value(this,4);
    hihi : constant double_float := double_float(val(0));
    lohi : constant double_float := double_float(val(1));
    hilo : constant double_float := double_float(val(2));
    lolo : constant double_float := double_float(val(3));

  begin
    to_that := create(hihi,lohi,hilo,lolo);
  end Assign;

  procedure Assign ( ada_cf : in Standard_Complex_Numbers.Complex_Number;
                     c_cf : in C_dblarrs.Pointer ) is

    use Standard_Complex_Numbers;

    val : C_Double_Array(0..1);

  begin
    val(0) := Interfaces.C.double(REAL_PART(ada_cf));
    val(1) := Interfaces.C.double(IMAG_PART(ada_cf));
    C_dblarrs.Copy_Array(val(0)'unchecked_access,c_cf,2);
  end Assign;

  procedure Assign ( ada_cf : in DoblDobl_Complex_Numbers.Complex_Number;
                     c_cf : in C_dblarrs.Pointer ) is

    use DoblDobl_Complex_Numbers;

    r : constant double_double := REAL_PART(ada_cf);
    i : constant double_double := IMAG_PART(ada_cf);

    val : C_Double_Array(0..3);

  begin
    val(0) := Interfaces.C.double(hi_part(r));
    val(1) := Interfaces.C.double(lo_part(r));
    val(2) := Interfaces.C.double(hi_part(i));
    val(3) := Interfaces.C.double(lo_part(i));
    C_dblarrs.Copy_Array(val(0)'unchecked_access,c_cf,4);
  end Assign;

  procedure Assign ( ada_cf : in QuadDobl_Complex_Numbers.Complex_Number;
                     c_cf : in C_dblarrs.Pointer ) is

    use QuadDobl_Complex_Numbers;

    r : constant quad_double := REAL_PART(ada_cf);
    i : constant quad_double := IMAG_PART(ada_cf);

    val : C_Double_Array(0..7);

  begin
    val(0) := Interfaces.C.double(hihi_part(r));
    val(1) := Interfaces.C.double(lohi_part(r));
    val(2) := Interfaces.C.double(hilo_part(r));
    val(3) := Interfaces.C.double(lolo_part(r));
    val(4) := Interfaces.C.double(hihi_part(i));
    val(5) := Interfaces.C.double(lohi_part(i));
    val(6) := Interfaces.C.double(hilo_part(i));
    val(7) := Interfaces.C.double(lolo_part(i));
    C_dblarrs.Copy_Array(val(0)'unchecked_access,c_cf,8);
  end Assign;

  procedure Assign ( c_cf : in C_dblarrs.Pointer;
                     ada_cf : out Standard_Complex_Numbers.Complex_Number ) is

    use Standard_Complex_Numbers;

    val : constant C_Double_Array := C_dblarrs.Value(c_cf,2);

  begin
    ada_cf := Create(double_float(val(0)),double_float(val(1)));
  end Assign;

  procedure Assign ( c_cf : in C_dblarrs.Pointer;
                     ada_cf : out DoblDobl_Complex_Numbers.Complex_Number ) is

    val : constant C_Double_Array := C_dblarrs.Value(c_cf,4);

    r : constant double_double 
      := create(double_float(val(0)),double_float(val(1)));
    i : constant double_double
      := create(double_float(val(0)),double_float(val(1)));

  begin
    ada_cf := DoblDobl_Complex_Numbers.Create(r,i);
  end Assign;

  procedure Assign ( c_cf : in C_dblarrs.Pointer;
                     ada_cf : out QuadDobl_Complex_Numbers.Complex_Number ) is

    val : constant C_Double_Array := C_dblarrs.Value(c_cf,8);

    r : constant quad_double
      := create(double_float(val(0)),double_float(val(1)),
                double_float(val(2)),double_float(val(3)));
    i : constant quad_double
      := create(double_float(val(4)),double_float(val(5)),
                double_float(val(6)),double_float(val(7)));

  begin
    ada_cf := QuadDobl_Complex_Numbers.Create(r,i);
  end Assign;

  procedure Assign ( ada_d : in Standard_Natural_Vectors.Vector;
                     c_d : in C_intarrs.Pointer ) is

    nd : constant natural32 := natural32(ada_d'last);
    val : C_Integer_Array(0..Interfaces.C.size_t(nd-1));
    ind : Interfaces.C.size_t := 0;

  begin
    for i in ada_d'range loop
      val(ind) := Interfaces.C.int(ada_d(i));
      ind := ind + 1;
    end loop;
    C_intarrs.Copy_Array(val(0)'unchecked_access,c_d,
                         Interfaces.C.ptrdiff_t(nd));
  end Assign;

  procedure Assign ( ada_d : in Standard_Integer_Vectors.Vector;
                     c_d : in C_intarrs.Pointer ) is

    nd : constant natural32 := natural32(ada_d'last);
    val : C_Integer_Array(0..Interfaces.C.size_t(nd-1));
    ind : Interfaces.C.size_t := 0;

  begin
    for i in ada_d'range loop
      val(ind) := Interfaces.C.int(ada_d(i));
      ind := ind + 1;
    end loop;
    C_intarrs.Copy_Array(val(0)'unchecked_access,c_d,
                         Interfaces.C.ptrdiff_t(nd));
  end Assign;

  procedure Assign ( ada_d : in Standard_Floating_Vectors.Vector;
                     c_d : in C_dblarrs.Pointer ) is

    nd : constant natural32 := natural32(ada_d'last);
    val : C_Double_Array(0..Interfaces.C.size_t(nd-1));
    ind : Interfaces.C.size_t := 0;

  begin
    for i in ada_d'range loop
      val(ind) := Interfaces.C.double(ada_d(i));
      ind := ind + 1;
    end loop;
    C_dblarrs.Copy_Array(val(0)'unchecked_access,c_d,
                         Interfaces.C.ptrdiff_t(nd));
  end Assign;

  procedure Assign ( ada_d : in Double_Double_Vectors.Vector;
                     c_d : in C_dblarrs.Pointer ) is

    nd : constant natural32 := natural32(ada_d'last);
    val : C_Double_Array(0..Interfaces.C.size_t(2*nd-1));
    ind : Interfaces.C.size_t := 0;
    nbr : double_double;

  begin
    for i in ada_d'range loop
      nbr := ada_d(i);
      val(ind) := Interfaces.C.double(hi_part(nbr)); ind := ind + 1;
      val(ind) := Interfaces.C.double(lo_part(nbr)); ind := ind + 1;
    end loop;
    C_dblarrs.Copy_Array(val(0)'unchecked_access,c_d,
                         Interfaces.C.ptrdiff_t(2*nd));
  end Assign;

  procedure Assign ( ada_d : in Quad_Double_Vectors.Vector;
                     c_d : in C_dblarrs.Pointer ) is

    nd : constant natural32 := natural32(ada_d'last);
    val : C_Double_Array(0..Interfaces.C.size_t(4*nd-1));
    ind : Interfaces.C.size_t := 0;
    nbr : quad_double;

  begin
    for i in ada_d'range loop
      nbr := ada_d(i);
      val(ind) := Interfaces.C.double(hihi_part(nbr)); ind := ind + 1;
      val(ind) := Interfaces.C.double(lohi_part(nbr)); ind := ind + 1;
      val(ind) := Interfaces.C.double(hilo_part(nbr)); ind := ind + 1;
      val(ind) := Interfaces.C.double(lolo_part(nbr)); ind := ind + 1;
    end loop;
    C_dblarrs.Copy_Array(val(0)'unchecked_access,c_d,
                         Interfaces.C.ptrdiff_t(4*nd));
  end Assign;

  procedure Assign ( ada_d : in Standard_Complex_Vectors.Vector;
                     c_d : in C_dblarrs.Pointer ) is

    use Standard_Complex_Numbers;

    v_n : constant natural32 := natural32(ada_d'last);
    val : C_Double_Array(0..Interfaces.C.size_t(2*v_n-1))
        := C_dblarrs.Value(c_d);
    ind : Interfaces.C.size_t := 0;

  begin
    for i in 1..integer32(v_n) loop
      val(ind) := Interfaces.C.double(REAL_PART(ada_d(i)));
      ind := ind + 1;
      val(ind) := Interfaces.C.double(IMAG_PART(ada_d(i)));
      ind := ind + 1;
    end loop;    
    C_dblarrs.Copy_Array(val(0)'unchecked_access,c_d,
                         Interfaces.C.ptrdiff_t(2*v_n));
  end Assign;

  procedure Assign ( ada_d : in DoblDobl_Complex_Vectors.Vector;
                     c_d : in C_dblarrs.Pointer ) is

    use DoblDobl_Complex_Numbers;

    v_n : constant natural32 := natural32(ada_d'last);
    val : C_Double_Array(0..Interfaces.C.size_t(4*v_n-1))
        := C_dblarrs.Value(c_d);
    ind : Interfaces.C.size_t := 0;
    nbr : double_double;

  begin
    for i in 1..integer32(v_n) loop
      nbr := REAL_PART(ada_d(i));
      val(ind) := Interfaces.C.double(hi_part(nbr)); ind := ind + 1;
      val(ind) := Interfaces.C.double(lo_part(nbr)); ind := ind + 1;
      nbr := IMAG_PART(ada_d(i));
      val(ind) := Interfaces.C.double(hi_part(nbr)); ind := ind + 1;
      val(ind) := Interfaces.C.double(lo_part(nbr)); ind := ind + 1;
    end loop;    
    C_dblarrs.Copy_Array(val(0)'unchecked_access,c_d,
                         Interfaces.C.ptrdiff_t(4*v_n));
  end Assign;

  procedure Assign ( ada_d : in QuadDobl_Complex_Vectors.Vector;
                     c_d : in C_dblarrs.Pointer ) is

    use QuadDobl_Complex_Numbers;

    v_n : constant natural32 := natural32(ada_d'last);
    val : C_Double_Array(0..Interfaces.C.size_t(8*v_n-1))
        := C_dblarrs.Value(c_d);
    ind : Interfaces.C.size_t := 0;
    nbr : quad_double;

  begin
    for i in 1..integer32(v_n) loop
      nbr := REAL_PART(ada_d(i));
      val(ind) := Interfaces.C.double(hihi_part(nbr)); ind := ind + 1;
      val(ind) := Interfaces.C.double(lohi_part(nbr)); ind := ind + 1;
      val(ind) := Interfaces.C.double(hilo_part(nbr)); ind := ind + 1;
      val(ind) := Interfaces.C.double(lolo_part(nbr)); ind := ind + 1;
      nbr := IMAG_PART(ada_d(i));
      val(ind) := Interfaces.C.double(hihi_part(nbr)); ind := ind + 1;
      val(ind) := Interfaces.C.double(lohi_part(nbr)); ind := ind + 1;
      val(ind) := Interfaces.C.double(hilo_part(nbr)); ind := ind + 1;
      val(ind) := Interfaces.C.double(lolo_part(nbr)); ind := ind + 1;
    end loop;    
    C_dblarrs.Copy_Array(val(0)'unchecked_access,c_d,
                         Interfaces.C.ptrdiff_t(8*v_n));
  end Assign;

  procedure Assign ( v_n : in natural32; c_d : in C_intarrs.Pointer;
                     ada_d : out Standard_Natural_Vectors.Vector ) is

    val : constant C_Integer_Array(0..Interfaces.C.size_t(v_n-1))
        := C_intarrs.Value(c_d,Interfaces.C.ptrdiff_t(v_n));
    ind : Interfaces.C.size_t := 0;

  begin
    for i in 1..integer32(v_n) loop
      ada_d(i) := natural32(val(ind));
      ind := ind + 1;
    end loop;
  end Assign;

  procedure Assign ( v_n : in natural32; c_d : in C_intarrs.Pointer;
                     ada_d : out Standard_Integer_Vectors.Vector ) is

    val : constant C_Integer_Array(0..Interfaces.C.size_t(v_n-1))
        := C_intarrs.Value(c_d,Interfaces.C.ptrdiff_t(v_n));
    ind : Interfaces.C.size_t := 0;

  begin
    for i in 1..integer32(v_n) loop
      ada_d(i) := integer32(val(ind));
      ind := ind + 1;
    end loop;
  end Assign;

  procedure Assign ( v_n : in natural32; c_d : in C_dblarrs.Pointer;
                     ada_d : out Standard_Floating_Vectors.Vector ) is

    val : constant C_Double_Array(0..Interfaces.C.size_t(v_n-1))
        := C_dblarrs.Value(c_d,Interfaces.C.ptrdiff_t(v_n));
    ind : Interfaces.C.size_t := 0;

  begin
    for i in 1..integer32(v_n) loop
      ada_d(i) := double_float(val(ind));
      ind := ind + 1;
    end loop;
  end Assign;

  procedure Assign ( v_n : in natural32; c_d : in C_dblarrs.Pointer;
                     ada_d : out Double_Double_Vectors.Vector ) is

    val : constant C_Double_Array(0..Interfaces.C.size_t(v_n-1))
        := C_dblarrs.Value(c_d,Interfaces.C.ptrdiff_t(v_n));
    ind : Interfaces.C.size_t := 0;
    lo,hi : double_float;
    dim : constant integer32 := integer32(v_n)/2;

  begin
    for i in 1..dim loop
      hi := double_float(val(ind)); ind := ind + 1;
      lo := double_float(val(ind)); ind := ind + 1;
      ada_d(i) := Double_Double_Numbers.Create(hi,lo);
    end loop;
  end Assign;

  procedure Assign ( v_n : in natural32; c_d : in C_dblarrs.Pointer;
                     ada_d : out Quad_Double_Vectors.Vector ) is

    val : constant C_Double_Array(0..Interfaces.C.size_t(v_n-1))
        := C_dblarrs.Value(c_d,Interfaces.C.ptrdiff_t(v_n));
    ind : Interfaces.C.size_t := 0;
    hihi,lohi,hilo,lolo : double_float;
    dim : constant integer32 := integer32(v_n)/4;

  begin
    for i in 1..dim loop
      hihi := double_float(val(ind)); ind := ind + 1;
      lohi := double_float(val(ind)); ind := ind + 1;
      hilo := double_float(val(ind)); ind := ind + 1;
      lolo := double_float(val(ind)); ind := ind + 1;
      ada_d(i) := Quad_Double_Numbers.Create(hihi,lohi,hilo,lolo);
    end loop;
  end Assign;

  procedure Assign ( v_n : in natural32; c_d : in C_dblarrs.Pointer;
                     ada_d : out Standard_Complex_Vectors.Vector ) is

    use Standard_Complex_Numbers;

    val : constant C_Double_Array(0..Interfaces.C.size_t(v_n-1))
        := C_dblarrs.Value(c_d,Interfaces.C.ptrdiff_t(v_n));
    ind : Interfaces.C.size_t := 0;

  begin
    for i in 1..integer32(v_n)/2 loop
      ada_d(i) := Create(double_float(val(ind)),double_float(val(ind+1)));
      ind := ind + 2;
    end loop;
  end Assign;

  procedure Assign ( v_n : in natural32; c_d : in C_dblarrs.Pointer;
                     ada_d : out DoblDobl_Complex_Vectors.Vector ) is

    use DoblDobl_Complex_Numbers;

    val : constant C_Double_Array(0..Interfaces.C.size_t(v_n-1))
        := C_dblarrs.Value(c_d,Interfaces.C.ptrdiff_t(v_n));
    ind : Interfaces.C.size_t := 0;
    lo,hi : double_float;
    re,im : double_double;
    dim : constant integer32 := integer32(v_n)/4;

  begin
    for i in 1..dim loop
      hi := double_float(val(ind)); ind := ind + 1;
      lo := double_float(val(ind)); ind := ind + 1;
      re := Double_Double_Numbers.Create(hi,lo);
      hi := double_float(val(ind)); ind := ind + 1;
      lo := double_float(val(ind)); ind := ind + 1;
      im := Double_Double_Numbers.Create(hi,lo);
      ada_d(i) := Create(re,im);
    end loop;
  end Assign;

  procedure Assign ( v_n : in natural32; c_d : in C_dblarrs.Pointer;
                     ada_d : out QuadDobl_Complex_Vectors.Vector ) is

    use QuadDobl_Complex_Numbers;

    val : constant C_Double_Array(0..Interfaces.C.size_t(v_n-1))
        := C_dblarrs.Value(c_d,Interfaces.C.ptrdiff_t(v_n));
    ind : Interfaces.C.size_t := 0;
    lolo,hilo,lohi,hihi : double_float;
    re,im : quad_double;
    dim : constant integer32 := integer32(v_n)/8;

  begin
    for i in 1..dim loop
      hihi := double_float(val(ind)); ind := ind + 1;
      lohi := double_float(val(ind)); ind := ind + 1;
      hilo := double_float(val(ind)); ind := ind + 1;
      lolo := double_float(val(ind)); ind := ind + 1;
      re := Quad_Double_Numbers.Create(hihi,lohi,hilo,lolo);
      hihi := double_float(val(ind)); ind := ind + 1;
      lohi := double_float(val(ind)); ind := ind + 1;
      hilo := double_float(val(ind)); ind := ind + 1;
      lolo := double_float(val(ind)); ind := ind + 1;
      im := Quad_Double_Numbers.Create(hihi,lohi,hilo,lolo);
      ada_d(i) := Create(re,im);
    end loop;
  end Assign;

  function C_Integer_Array_to_String
             ( n : natural32; v : C_Integer_Array ) return String is

    res : String(1..integer(n));
    ch : character;

  begin
    for i in v'range loop
      exit when (integer(i)+1 > res'last);
      ch := Integer_to_Character(integer32(v(i)));
      res(integer(i)+1) := ch;
    end loop;
    return res;
  end C_Integer_Array_to_String;

  function Pad_with_Spaces ( n : natural32; s : string ) return string is
  begin
    if s'last >= integer(n) then
      return s;
    else
      declare
        res : string(1..integer(n));
      begin
        res(s'range) := s;
        for i in s'last+1..integer(n) loop
          res(i) := ' ';
        end loop;
        return res;
      end;
    end if;
  end Pad_with_Spaces;

  function String_to_C_Integer_Array
             ( n : natural32; s : string ) return C_Integer_Array is


    res : C_Integer_Array(0..Interfaces.C.size_t(n-1));
    ind : Interfaces.C.size_t := 0;

  begin
    for i in s'range loop
      res(ind) := Interfaces.C.int(Character_to_Integer(s(i)));
      ind := ind + 1;
    end loop;
    return res;
  end String_to_C_Integer_Array;

  function String_to_Integer_Vector
             ( s : string ) return Standard_Integer_Vectors.Vector is

    res : Standard_Integer_Vectors.Vector
            (integer32(s'first)..integer32(s'last));

  begin
    for i in s'range loop
      res(integer32(i)) := Character_to_Integer(s(i));
    end loop;
    return res;
  end String_to_Integer_Vector;

end Assignments_in_Ada_and_C;
