with text_io;                           use text_io;
with Communications_with_User;          use Communications_with_User;
with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;       use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;         use Standard_Floating_Numbers;
with Standard_Floating_Numbers_io;      use Standard_Floating_Numbers_io;
with Multprec_Floating_Numbers;         use Multprec_Floating_Numbers;
with Multprec_Floating_Numbers_io;      use Multprec_Floating_Numbers_io;
with Standard_Complex_Vectors;
with Standard_Random_Vectors;           use Standard_Random_Vectors;
with Standard_Complex_Matrices;
with Multprec_Complex_Vectors;
with Multprec_Random_Vectors;           use Multprec_Random_Vectors;
with Multprec_Complex_Matrices;
with Standard_Complex_Polynomials;      use Standard_Complex_Polynomials;
with Standard_Complex_Poly_Systems;
with Standard_Complex_Poly_Systems_io;  use Standard_Complex_Poly_Systems_io;
with Standard_Complex_Poly_Matrices;
with Standard_to_Multprec_Convertors;
with Multprec_Complex_Polynomials;      use Multprec_Complex_Polynomials;
with Multprec_Complex_Poly_Systems;
with Multprec_Complex_Poly_Systems_io;  use Multprec_Complex_Poly_Systems_io;
with Multprec_Complex_Poly_Matrices;
with Standard_Complex_Solutions;        use Standard_Complex_Solutions;
with Standard_System_and_Solutions_io;
with Multprec_Complex_Solutions;        use Multprec_Complex_Solutions;
with Multprec_System_and_Solutions_io;
with Standard_Deflation_Trees_io;       use Standard_Deflation_Trees_io;
with Standard_Query_Matrices;
with Standard_Nullity_Matrices;
with Standard_Multiple_Deflation;
with Multprec_Query_Matrices;
with Multprec_Nullity_Matrices;
with Multprec_Multiple_Deflation;

procedure ts_muldef is

-- DESCRIPTION :
--   Implementation of higher-order deflation, based on multiplicity
--   structure, whence the contraction of the name "multiple deflation"
--   in the name of the procedure.

  procedure Standard_Evaluation_of_Nullity_Matrix_at_Random_Points
              ( p : in Standard_Complex_Poly_Systems.Poly_Sys;
                d,nq,nv,nr,nc : in natural32 ) is

    use Standard_Query_Matrices;
    use Standard_Nullity_Matrices;

    ans : character;
    nulmat : Standard_Complex_Poly_Matrices.Matrix
               (1..integer32(nr),1..integer32(nc));
    x : Standard_Complex_Vectors.Vector(1..integer32(nv));
    y1,y2 : Standard_Complex_Matrices.Matrix
              (1..integer32(nr),1..integer32(nc));
    diff : double_float;

  begin
    put_line("Creating of the symbolic nullity matrix ...");
    nulmat := Create_Nullity_Matrix(standard_output,nq,nv,nr,nc,d,p);
    put_line("Evaluation at random points ...");
    loop
      x := Random_Vector(1,integer32(nv));
      y1 := Evaluate_Nullity_Matrix(nq,nv,nr,nc,d,p,x);
      y2 := Eval0(nulmat,x);
      diff := Difference_of_Matrices(y1,y2);
      put("Difference between the two matrices : "); put(diff,3); new_line;
      Show_Differences_in_Matrices(y1,y2);
      put("Do you want more tests ? (y/n) ");
      Ask_Yes_or_No(ans);
      exit when (ans /= 'y');
    end loop;
  end Standard_Evaluation_of_Nullity_Matrix_at_Random_Points;

  procedure Multprec_Evaluation_of_Nullity_Matrix_at_Random_Points
              ( p : in Multprec_Complex_Poly_Systems.Poly_Sys;
                sz,d,nq,nv,nr,nc : in natural32 ) is

    use Multprec_Query_Matrices;
    use Multprec_Nullity_Matrices;

    ans : character;
    nulmat : Multprec_Complex_Poly_Matrices.Matrix
               (1..integer32(nr),1..integer32(nc));
    x : Multprec_Complex_Vectors.Vector(1..integer32(nv));
    y1,y2 : Multprec_Complex_Matrices.Matrix
              (1..integer32(nr),1..integer32(nc));
    diff : Floating_Number;

  begin
    put_line("Creating of the symbolic nullity matrix ...");
    nulmat := Create_Nullity_Matrix(standard_output,nq,nv,nr,nc,d,p);
    put_line("Evaluation at random points ...");
    loop
      x := Random_Vector(1,integer32(nv),sz);
      y1 := Evaluate_Nullity_Matrix(nq,nv,nr,nc,d,p,x);
      y2 := Eval0(nulmat,x);
      diff := Difference_of_Matrices(y1,y2);
      put("Difference between the two matrices : "); put(diff,3); new_line;
      Show_Differences_in_Matrices(y1,y2);
      put("Do you want more tests ? (y/n) ");
      Ask_Yes_or_No(ans);
      exit when (ans /= 'y');
    end loop;
  end Multprec_Evaluation_of_Nullity_Matrix_at_Random_Points;

  procedure Standard_Evaluation_of_Nullity_Matrix
              ( p : in Standard_Complex_Poly_Systems.Poly_Sys ) is
 
  -- DESCRIPTION :
  --   Test evaluation of nullity matrices of order d on the system p.

    use Standard_Nullity_Matrices;

    nq : constant natural32 := natural32(p'last);
    nv : constant natural32 := Number_of_Unknowns(p(p'first));
    d,nr,nc : natural32 := 0;
 
  begin
    put("Give the order of the deflation : "); get(d);
    Dimensions_of_Nullity_Matrix(nq,nv,d,nr,nc);
    put("Nullity matrix of order "); put(d,1); put(" is ");
    put(nr,1); put("-by-"); put(nc,1); put_line("-matrix.");
    put_line("Testing the creation and evaluation ...");
    Standard_Evaluation_of_Nullity_Matrix_at_Random_Points(p,d,nq,nv,nr,nc);
  end Standard_Evaluation_of_Nullity_Matrix;

  procedure Multprec_Evaluation_of_Nullity_Matrix
              ( p : in Multprec_Complex_Poly_Systems.Poly_Sys;
                s : natural32 ) is
 
  -- DESCRIPTION :
  --   Test evaluation of nullity matrices of order d on the system p,
  --   using multiprecision arithmetic with numbers of size s.

    use Multprec_Nullity_Matrices;

    nq : constant natural32 := natural32(p'last);
    nv : constant natural32 := Number_of_Unknowns(p(p'first));
    d,nr,nc : natural32 := 0;
 
  begin
    put("Give the order of the deflation : "); get(d);
    Dimensions_of_Nullity_Matrix(nq,nv,d,nr,nc);
    put("Nullity matrix of order "); put(d,1); put(" is ");
    put(nr,1); put("-by-"); put(nc,1); put_line("-matrix.");
    put_line("Testing the creation and evaluation ...");
    Multprec_Evaluation_of_Nullity_Matrix_at_Random_Points(p,s,d,nq,nv,nr,nc);
  end Multprec_Evaluation_of_Nullity_Matrix;

  procedure Test_Symbolic_Creation_of_Deflation
              ( p : in Standard_Complex_Poly_Systems.Poly_Sys ) is

  -- DESCRIPTION :
  --   Creates a symbolic deflation of the given system.

    use Standard_Complex_Poly_Systems;
    use Standard_Multiple_Deflation;

    file : file_type;
    d,r : natural32 := 0;

  begin
    new_line;
    put_line("Reading the name of the output file ...");
    Read_Name_and_Create_File(file);
    new_line;
    put("Give the order of the deflation : "); get(d);
    put("Give the corank of the nullity matrix : "); get(r);
    put_line("The system on input : "); put(p);
    declare
      hdp : constant Poly_Sys := Symbolic_Deflate(p,d,r);
      nq : constant natural32 := natural32(hdp'last);
      nv : constant natural32 := Number_of_Unknowns(hdp(hdp'first));
      nm : constant natural32 := nv - Number_of_Unknowns(p(p'first));
    begin
      Add_Multiplier_Symbols(1,nm);
      put_line(file,hdp);
    end;
    close(file);
  end Test_Symbolic_Creation_of_Deflation;

  procedure Adjust_Tolerance ( tol : in out double_float ) is

    ans : character;

  begin
    loop
      put("The current value of the tolerance is ");
      put(tol,3); put_line(".");
      put("Do you want to change this value ? (y/n) ");
      Ask_Yes_or_No(ans);
      exit when (ans /= 'y');
      put("Give a new value for the tolerance : "); get(tol);
    end loop;
  end Adjust_Tolerance;

  procedure Standard_Test_Order_Prediction
              ( p : in Standard_Complex_Poly_Systems.Poly_Sys;
                sols : in out Standard_Complex_Solutions.Solution_List ) is

    use Standard_Multiple_Deflation;

    tol : double_float := 1.0E-5;
    d,c : natural32 := 0;

  begin
    Adjust_Tolerance(tol);
    Predict_Order(standard_output,p,Head_Of(sols).v,tol,c,d);
    put("The predicted order of deflation is "); put(d,1); put_line(".");
  end Standard_Test_Order_Prediction;

  procedure Multprec_Test_Order_Prediction
              ( p : in Multprec_Complex_Poly_Systems.Poly_Sys;
                sols : in out Multprec_Complex_Solutions.Solution_List;
                size : in natural32 ) is

    use Multprec_Multiple_Deflation;

    tol : double_float := 1.0E-5;
    d,c : natural32 := 0;

  begin
    Adjust_Tolerance(tol);
    Predict_Order(standard_output,p,Head_Of(sols).v,size,tol,c,d);
    put("The predicted order of deflation is "); put(d,1); put_line(".");
  end Multprec_Test_Order_Prediction;

  procedure Standard_Execute_Multiple_Deflation
              ( p : in Standard_Complex_Poly_Systems.Poly_Sys;
                sols : in out Standard_Complex_Solutions.Solution_List ) is

    use Standard_Multiple_Deflation;

    file : file_type;
    tol : double_float := 1.0E-5;

  begin
    new_line;
    put_line("Reading the name of the output file...");
    Read_Name_and_Create_File(file);
    Adjust_Tolerance(tol);
    Interactive_Symbolic_Deflation(file,p,sols,tol);
  end Standard_Execute_Multiple_Deflation;

  procedure Multprec_Execute_Multiple_Deflation
              ( p : in Multprec_Complex_Poly_Systems.Poly_Sys;
                sols : in out Multprec_Complex_Solutions.Solution_List;
                size : in natural32 ) is

    use Multprec_Multiple_Deflation;

    file : file_type;
    tol : double_float := 1.0E-5;

  begin
    new_line;
    put_line("Reading the name of the output file...");
    Read_Name_and_Create_File(file);
    Adjust_Tolerance(tol);
    Interactive_Symbolic_Deflation(file,p,sols,size,tol);
  end Multprec_Execute_Multiple_Deflation;

  procedure Main is

    stlp : Standard_Complex_Poly_Systems.Link_to_Poly_Sys;
    mplp : Multprec_Complex_Poly_Systems.Link_to_Poly_Sys;
    stsols : Standard_Complex_Solutions.Solution_List;
    mpsols : Multprec_Complex_Solutions.Solution_List;
    deci,size : natural32 := 0;
    ans : character;

  begin
    new_line;
    put_line("Testing higher-order deflation methods...");
    new_line;
    put_line("Choose one of the following : ");
    put_line("  1. Evaluation of nullity matrix in standard arithmetic;");
    put_line("  2.                              in multiprecision arithmetic;");
    put_line("  3. Creation of a symbolic higher order deflation;");
    put_line("  4. Predict order of deflation in standard arithmetic;");
    put_line("  5.                            in multiprecision arithmetic;");
    put_line("  6. Higher order deflation in standard arithmetic;");
    put_line("  7.                        in multiprecision arithmetic.");
    put("Type 1, 2, 3, 4, 5, 6, or 7 to make your choice : ");
    Ask_Alternative(ans,"1234567");
    case ans is
      when '1' => new_line; get(stlp);
                  Standard_Evaluation_of_Nullity_Matrix(stlp.all);
      when '2' => new_line; get(mplp);
                  put("Give the number of decimal places : "); get(deci);
                  size := Multprec_Floating_Numbers.Decimal_to_Size(deci);
                  Multprec_Evaluation_of_Nullity_Matrix(mplp.all,size);
      when '3' => new_line; get(stlp);
                  Test_Symbolic_Creation_of_Deflation(stlp.all);
      when '4' => Standard_System_and_Solutions_io.get(stlp,stsols);
                  Standard_Test_Order_Prediction(stlp.all,stsols);
      when '5' => Multprec_System_and_Solutions_io.get(mplp,mpsols);
                  put("Give the number of decimal places : "); get(deci);
                  skip_line;
                  size := Multprec_Floating_Numbers.Decimal_to_Size(deci);
                  Standard_to_Multprec_Convertors.Set_Size(mplp.all,size);
                  Multprec_Complex_Solutions.Set_Size(mpsols,size);
                  Multprec_Test_Order_Prediction(mplp.all,mpsols,size);
      when '6' => Standard_System_and_Solutions_io.get(stlp,stsols);
                  Standard_Execute_Multiple_Deflation(stlp.all,stsols);
      when '7' => Multprec_System_and_Solutions_io.get(mplp,mpsols);
                  put("Give the number of decimal places : "); get(deci);
                  skip_line;
                  size := Multprec_Floating_Numbers.Decimal_to_Size(deci);
                  Standard_to_Multprec_Convertors.Set_Size(mplp.all,size);
                  Multprec_Complex_Solutions.Set_Size(mpsols,size);
                  Multprec_Execute_Multiple_Deflation(mplp.all,mpsols,size);
      when others => null;
    end case;
  end Main;

begin
  Main;
end ts_muldef;
