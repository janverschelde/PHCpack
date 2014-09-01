with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Characters_and_Numbers;
with Multprec_Integer_Numbers;          use Multprec_Integer_Numbers;
with Multprec_Integer_Numbers_io;

package body Point_Lists_and_Strings is

-- AUXILIARY ROUTINES :

  function Add_Coordinates
             ( A : Standard_Integer64_Matrices.Matrix; k,i : integer32;
               accu : string ) return string is

  -- DESCRIPTION :
  --   Returns accu added with the i-th coordinate of the point defined 
  --   by the k-th column of A.

  -- REQUIRED : k <= A'last(2) and i <= A'last(1).

    strAik : constant string := Characters_and_Numbers.Convert(A(i,k));

  begin
    if i = A'last(1) then
      return accu & strAik & ")";
    else
      return Add_Coordinates(A,k,i+1,accu & strAik & ", ");
    end if;
  end Add_Coordinates;

  function One_if_Negative ( i : Integer_Number ) return natural32 is

  -- DESCRIPTION :
  --   Returns one if i < 0, returns 0 otherwise.

  begin
    if i < 0
     then return 1;
     else return 0;
    end if;
  end One_if_Negative;

  function Convert_to_String ( i : Integer_Number ) return string is

    dp : constant natural32 := Multprec_Integer_Numbers.Decimal_Places(i);
    dp1 : constant natural32 := dp + One_if_Negative(i);
    res : string(1..integer(dp1));

  begin
    Multprec_Integer_Numbers_io.put(res,i);
    return res;
  end Convert_to_String;

  function Add_Coordinates
             ( A : Multprec_Integer_Matrices.Matrix; k,i : integer32;
               accu : string ) return string is

  -- DESCRIPTION :
  --   Returns accu added with the i-th coordinate of the point defined 
  --   by the k-th column of A.

  -- REQUIRED : k <= A'last(2) and i <= A'last(1).

    strAik : constant string := Convert_to_String(A(i,k));

  begin
    if i = A'last(1) then
      return accu & strAik & ")";
    else
      return Add_Coordinates(A,k,i+1,accu & strAik & ", ");
    end if;
  end Add_Coordinates;

  function Add_Point
             ( A : Standard_Integer64_Matrices.Matrix; k : integer32;
               accu : string ) return string is

  -- DESCRIPTION :
  --   Returns accu added with the point defined by the k-th column of A.

    point : constant string := Add_Coordinates(A,k,A'first(1),"(");

  begin
    if k = A'last(2) then
      return accu & point;
    else
      return Add_Point(A,k+1,accu & point & ", ");
    end if;
  end Add_Point;

  function Add_Point
             ( A : Multprec_Integer_Matrices.Matrix; k : integer32;
               accu : string ) return string is

  -- DESCRIPTION :
  --   Returns accu added with the point defined by the k-th column of A.

    point : constant string := Add_Coordinates(A,k,A'first(1),"(");

  begin
    if k = A'last(2) then
      return accu & point;
    else
      return Add_Point(A,k+1,accu & point & ", ");
    end if;
  end Add_Point;

-- TARGET FUNCTIONS :

  function convert ( A : Standard_Integer64_Matrices.Matrix ) return string is

    res : constant string := Add_Point(A,A'first(2),"[");

  begin
    return res & "]";
  end convert;

  function convert ( A : Multprec_Integer_Matrices.Matrix ) return string is

    res : constant string := Add_Point(A,A'first(2),"[");

  begin
    return res & "]";
  end Convert;

end Point_Lists_and_Strings;
