with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Integer_Matrices;
with Standard_Integer_VecMats;
with Standard_Integer_VecMats_io;
with Standard_Integer_Matrices_io;
with Standard_Integer64_Matrices;
with Standard_Integer64_Matrices_io;
with Standard_Integer64_VecMats;
with Standard_Integer64_VecMats_io;
with Multprec_Integer_Vectors;
with Multprec_Integer_Vectors_io;
with Multprec_Integer_Matrices;
with Multprec_Integer_Matrices_io;

package body Test_Integer_Matrices is

  procedure Test_Standard_io is

    use Standard_Integer_Matrices;
    use Standard_Integer_Matrices_io;

    n,m : integer32 := 0;

  begin
    put("Give the number of rows : "); get(n);
    put("Give the number of columns : "); get(m);
    declare
      mat : Matrix(1..n,1..m);
    begin
      put("Give "); put(n,1); put("x"); put(m,1);
      put_line(" integer matrix : "); get(mat);
      put_line("Your matrix : "); put(mat);
    end;
  end Test_Standard_io;

  procedure Test_Standard64_io is

    use Standard_Integer64_Matrices;
    use Standard_Integer64_Matrices_io;

    n,m : integer32 := 0;

  begin
    put("Give the number of rows : "); get(n);
    put("Give the number of columns : "); get(m);
    declare
      mat : Matrix(1..n,1..m);
    begin
      put("Give "); put(n,1); put("x"); put(m,1);
      put_line(" integer matrix : "); get(mat);
      put_line("Your matrix : "); put(mat);
    end;
  end Test_Standard64_io;

  procedure Test_Standard_VecMat_io is

    use Standard_Integer_Matrices;
    use Standard_Integer_Matrices_io;
    use Standard_Integer_VecMats;
    use Standard_Integer_VecMats_io;

    n,n1,n2 : natural32 := 0;
    lv : Link_to_VecMat;

  begin
    put("Give the number of matrices : "); get(n);
    put("Give #rows : "); get(n1);
    put("Give #columns : "); get(n2);
    put("Give "); put(n,1);
    put(" "); put(n1,1);
    put("-by-"); put(n2,1);
    put_line(" integer matrices : ");
    get(n,n1,n2,lv);
    put_line("The vector of matrices :"); put(lv);
  end Test_Standard_VecMat_io;

  procedure Test_Standard64_VecMat_io is

    use Standard_Integer64_Matrices;
    use Standard_Integer64_Matrices_io;
    use Standard_Integer64_VecMats;
    use Standard_Integer64_VecMats_io;

    n,n1,n2 : natural32 := 0;
    lv : Link_to_VecMat;

  begin
    put("Give the number of matrices : "); get(n);
    put("Give #rows : "); get(n1);
    put("Give #columns : "); get(n2);
    put("Give "); put(n,1);
    put(" "); put(n1,1);
    put("-by-"); put(n2,1);
    put_line(" integer matrices : ");
    get(n,n1,n2,lv);
    put_line("The vector of matrices :"); put(lv);
  end Test_Standard64_VecMat_io;

  procedure Test_Multprec_io is

    use Multprec_Integer_Matrices;
    use Multprec_Integer_Matrices_io;

    n,m : integer32 := 0;

  begin
    put("Give the number of rows : "); get(n);
    put("Give the number of columns : "); get(m);
    declare
      mat : Matrix(1..n,1..m);
    begin
      put("Give "); put(n,1); put("x"); put(m,1);
      put_line(" integer matrix : "); get(mat);
      put_line("Your matrix : "); put(mat); new_line;
    end;
  end Test_Multprec_io;

  procedure Test_Multprec_Matrix_Vector_Product is

    use Multprec_Integer_Vectors;
    use Multprec_Integer_Vectors_io;
    use Multprec_Integer_Matrices;
    use Multprec_Integer_Matrices_io; 

    n,m : integer32 := 0;

  begin
    put("Give the number of rows : "); get(n);
    put("Give the number of columns : "); get(m);
    declare
      mat : Matrix(1..n,1..m);
      vec : Vector(1..m);
      prod : Vector(1..n);
    begin
      put("Give "); put(n,1); put("x"); put(m,1);
      put_line(" integer matrix : "); get(mat);
      put_line("Your matrix : "); put(mat);   
      put("Give an " ); put(m,1); put("-vector : "); get(vec);
      put("Your vector : "); put(vec); new_line;
      prod := mat*vec;
      put_line("The matrix-vector product : "); put(prod); new_line;
    end;
  end Test_Multprec_Matrix_Vector_Product;

  procedure Main is

    ans,lng : character;

  begin
    new_line;
    put_line("Interactive testing of matrices of integer numbers");
    new_line;
    loop
      put_line("Choose one of the following : ");
      put_line("  0. exit this program.");
      put_line("  1. io of matrices of standard numbers.");
      put_line("  2. io of vectors of matrices of standard numbers.");
      put_line("  3. io of matrices of multi-precision numbers.");
      put_line("  4. test multi-precision matrix-vector product.");
      put("Type 0, 1, 2, 3, or 4 to choose: ");
      Ask_Alternative(ans,"01234");
      exit when (ans = '0');
      if ans = '1' or ans = '2' 
       then put("Use 64-bit integers ? (y/n) "); Ask_Yes_or_No(lng);
      end if;
      new_line;
      case ans is
        when '1' => if lng = 'y'
                     then Test_Standard64_io;
                     else Test_Standard_io;
                    end if;
        when '2' => if lng = 'y'
                     then Test_Standard64_VecMat_io;
                     else Test_Standard_VecMat_io;
                    end if;
        when '3' => Test_Multprec_io;
        when '4' => Test_Multprec_Matrix_Vector_Product;
        when others => null;
      end case;
    end loop;
  end Main;

end Test_Integer_Matrices;
