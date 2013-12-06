with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Double_Double_Numbers;              use Double_Double_Numbers;
with Quad_Double_Numbers;                use Quad_Double_Numbers;
with DoblDobl_Complex_Numbers;
with DoblDobl_Complex_Numbers_io;        use DoblDobl_Complex_Numbers_io;
with DoblDobl_Complex_Vectors;
with DoblDobl_Complex_Vectors_io;        use DoblDobl_Complex_Vectors_io;
with DoblDobl_Complex_Matrices;
with DoblDobl_Random_Vectors;
with QuadDobl_Complex_Numbers;
with QuadDobl_Complex_Numbers_io;        use QuadDobl_Complex_Numbers_io;
with QuadDobl_Complex_Vectors;
with QuadDobl_Complex_Vectors_io;        use QuadDobl_Complex_Vectors_io;
with QuadDobl_Complex_Matrices;
with DoblDobl_Complex_Poly_Systems;
with DoblDobl_Complex_Poly_Systems_io;   use DoblDobl_Complex_Poly_Systems_io;
with QuadDobl_Complex_Poly_Systems;
with QuadDobl_Complex_Poly_Systems_io;   use QuadDobl_Complex_Poly_Systems_io;
with DoblDobl_Homotopy;
with QuadDobl_Homotopy;

procedure ts_qdhom is

  procedure Read_DoblDobl_Homotopy ( n : out integer32 ) is

    use DoblDobl_Complex_Numbers;
    use DoblDobl_Complex_Poly_Systems;

    lp,lq : Link_to_Poly_Sys;
    targetfile,startfile : file_type;
    one : constant double_double := create(1.0);
    gamma : constant Complex_Number := create(one);

  begin
    new_line;
    put_line("Reading the name of the file for the target system.");
    Read_Name_and_Open_File(targetfile);
    get(targetfile,lp);
    Close(targetfile);
    new_line;
    put_line("Reading the name of the file for the start system.");
    Read_Name_and_Open_File(startfile);
    get(startfile,lq);
    Close(startfile);
    DoblDobl_Homotopy.Create(lp.all,lq.all,1,gamma);
    n := lp'last;
  end Read_DoblDobl_Homotopy;

  procedure Read_QuadDobl_Homotopy ( n : out integer32 ) is

    use QuadDobl_Complex_Numbers;
    use QuadDobl_Complex_Poly_Systems;

    lp,lq : Link_to_Poly_Sys;
    targetfile,startfile : file_type;
    one : constant quad_double := create(1.0);
    gamma : constant Complex_Number := create(one);

  begin
    new_line;
    put_line("Reading the name of the file for the target system.");
    Read_Name_and_Open_File(targetfile);
    get(targetfile,lp);
    Close(targetfile);
    new_line;
    put_line("Reading the name of the file for the start system.");
    Read_Name_and_Open_File(startfile);
    get(startfile,lq);
    Close(startfile);
    QuadDobl_Homotopy.Create(lp.all,lq.all,1,gamma);
    n := lp'last;
  end Read_QuadDobl_Homotopy;

  procedure Compare_Evaluation ( n : in integer32 ) is

    x_dd,y_dd : DoblDobl_Complex_Vectors.Vector(1..n);
    x_qd,y_qd : QuadDobl_Complex_Vectors.Vector(1..n);
    a_dd : constant double_double := create(0.123123);
    t_dd : DoblDobl_Complex_Numbers.Complex_Number 
         := DoblDobl_Complex_Numbers.create(a_dd);
    a_qd : constant quad_double := create(0.123123);
    t_qd : QuadDobl_Complex_Numbers.Complex_Number
         := QuadDobl_Complex_Numbers.create(a_qd);
    m_dd : DoblDobl_Complex_Matrices.Matrix(1..n,1..n);
    m_qd : QuadDobl_Complex_Matrices.Matrix(1..n,1..n);

  begin
    x_dd := DoblDobl_Random_Vectors.Random_Vector(1,n);
    for i in x_qd'range loop
      declare
        re : quad_double 
           := create(DoblDobl_Complex_Numbers.REAL_PART(x_dd(i)));
        im : quad_double 
           := create(DoblDobl_Complex_Numbers.IMAG_PART(x_dd(i)));
      begin
        x_qd(i) := QuadDobl_Complex_Numbers.create(re,im);
      end;
    end loop;
    y_dd := DoblDobl_Homotopy.Eval(x_dd,t_dd);
    put_line(y_dd);
    y_qd := QuadDobl_Homotopy.Eval(x_qd,t_qd);
    put_line(y_qd);
    m_dd := DoblDobl_Homotopy.Diff(x_dd,t_dd);
    m_qd := QuadDobl_Homotopy.Diff(x_qd,t_qd);
    for i in 1..n loop
      for j in 1..n loop
        put("m_dd("); put(i,1); put(","); put(j,1); put(") : ");
        put(m_dd(i,j)); new_line;
        put("m_qd("); put(i,1); put(","); put(j,1); put(") : ");
        put(m_qd(i,j),32); new_line;
      end loop;
    end loop;
  end Compare_Evaluation;

  procedure Main is

    n : integer32;
    ans : character;

  begin
    new_line;
    put_line("Comparing homotopies with double double and quad doubles...");
    new_line;
    put_line("reading data for homotopy of double doubles...");
    Read_DoblDobl_Homotopy(n);
    new_line;
    put_line("reading data for homotopy of quad doubles...");
    Read_QuadDobl_Homotopy(n);
    loop
      Compare_Evaluation(n);
      put("more tests ? (y/n) "); Ask_Yes_or_No(ans);
      exit when ans /= 'y';
    end loop;
  end Main;

begin
  Main;
end ts_qdhom;
