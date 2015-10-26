with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Double_Double_Numbers_io;           use Double_Double_Numbers_io;
with DoblDobl_Complex_Numbers;           use DoblDobl_Complex_Numbers;
with DoblDobl_Complex_Numbers_io;        use DoblDobl_Complex_Numbers_io;
with Symbol_Table,Symbol_Table_io;

package body DoblDobl_Quad_Turn_Points_io is

  procedure Write_Vector ( x : in Double_Double_Vectors.Vector ) is
  begin
    Write_Vector(standard_output,x);
  end Write_Vector;

  procedure Write_Vector ( file : in file_type;
                           x : in Double_Double_Vectors.Vector ) is
  begin
    for i in x'range loop
      put(file,"  ");
      Symbol_Table_io.put(file,Symbol_Table.get(natural32(i)));
      put(file," : ");
      put(file,x(i)); new_line(file);
    end loop;
  end Write_Vector;

  procedure Write_Vector ( x : in DoblDobl_Complex_Vectors.Vector ) is
  begin
    Write_Vector(standard_output,x);
  end Write_Vector;

  procedure Write_Vector ( file : in file_type;
                           x : in DoblDobl_Complex_Vectors.Vector ) is
  begin
    for i in x'range loop
      put(file,"  ");
      Symbol_Table_io.put(file,Symbol_Table.get(natural32(i)));
      put(file," : ");
      put(file,x(i)); new_line(file);
    end loop;
  end Write_Vector;

  procedure Write_Tangent ( x : in Double_Double_Vectors.Vector ) is
  begin
    Write_Tangent(standard_output,x);
  end Write_Tangent;

  procedure Write_Tangent ( file : in file_type;
                            x : in Double_Double_Vectors.Vector ) is
  begin
    for i in x'range loop
      put(file,"  d");
      Symbol_Table_io.put(file,Symbol_Table.get(natural32(i)));
      put(file," : ");
      put(file,x(i)); new_line(file);
    end loop;
  end Write_Tangent;

  procedure Write_Tangent ( x : in DoblDobl_Complex_Vectors.Vector ) is
  begin
    Write_Tangent(standard_output,x);
  end Write_Tangent;

  procedure Write_Tangent ( file : in file_type;
                            x : in DoblDobl_Complex_Vectors.Vector ) is
  begin
    for i in x'range loop
      put(file,"  d");
      Symbol_Table_io.put(file,Symbol_Table.get(natural32(i)));
      put(file," : ");
      put(file,x(i)); new_line(file);
    end loop;
  end Write_Tangent;

  procedure Write_Vector_and_its_Tangent
              ( x,t : in Double_Double_Vectors.Vector ) is
  begin
    Write_Vector_and_its_Tangent(standard_output,x,t);
  end Write_Vector_and_its_Tangent;

  procedure Write_Vector_and_its_Tangent
              ( file : in file_type;
                x,t : in Double_Double_Vectors.Vector ) is
  begin
    for i in x'range loop
      put(file,"  ");
      Symbol_Table_io.put(file,Symbol_Table.get(natural32(i)));
      put(file," : ");
      put(file,x(i));
      put(file,"  |  d");
      Symbol_Table_io.put(file,Symbol_Table.get(natural32(i)));
      put(file," : ");
      put(file,t(i)); new_line(file);
    end loop;
  end Write_Vector_and_its_Tangent;

  procedure Write_Vector_Tangent_and_Determinants
              ( x,t,det : in Double_Double_Vectors.Vector ) is
  begin
    Write_Vector_Tangent_and_Determinants(standard_output,x,t,det);
  end Write_Vector_Tangent_and_Determinants;

  procedure Write_Vector_Tangent_and_Determinants
              ( file : in file_type;
                x,t,det : in Double_Double_Vectors.Vector ) is
  begin
    Write_Vector_and_its_Tangent(file,x,t);
    put(file,"Determinants :");
    for i in det'range loop
      put(file," "); put(file,det(i),1);
    end loop;
    new_line(file);
  end Write_Vector_Tangent_and_Determinants;

  procedure Write_Vector_and_its_Tangent
              ( x,t : in DoblDobl_Complex_Vectors.Vector ) is
  begin
    Write_Vector_and_its_Tangent(standard_output,x,t);
  end Write_Vector_and_its_Tangent;

  procedure Write_Vector_and_its_Tangent
              ( file : in file_type;
                x,t : in DoblDobl_Complex_Vectors.Vector ) is
  begin
    Write_Vector(file,x);
    Write_Tangent(file,t);
  end Write_Vector_and_its_Tangent;

  procedure Read_Initial_Vector
              ( x : in out Double_Double_Vectors.Vector ) is
  begin
    put("Reading "); put(x'last,1); put_line(" real initial values : ");
    for i in x'range loop
      put("Give value for ");
      Symbol_Table_io.put(Symbol_Table.get(natural32(i)));
      put(" : ");
      get(x(i));
    end loop;
    put_line("The initial vector is "); Write_Vector(x);
  end Read_Initial_Vector;

  procedure Read_Initial_Vector
              ( x : in out DoblDobl_Complex_Vectors.Vector ) is
  begin
    put("Reading "); put(x'last,1); put_line(" complex initial values : ");
    for i in x'range loop
      put("Give value for ");
      Symbol_Table_io.put(Symbol_Table.get(natural32(i)));
      put(" : ");
      get(x(i));
    end loop;
    put_line("The initial vector is "); Write_Vector(x);
  end Read_Initial_Vector;

  procedure Write_Corrector_Information
              ( x,y : in Double_Double_Vectors.Vector;
                err,res : in double_double ) is
  begin
    Write_Corrector_Information(standard_output,x,y,err,res);
  end Write_Corrector_Information;

  procedure Write_Corrector_Information
              ( x,y : in Double_Double_Vectors.Vector;
                err,res,det : in double_double ) is
  begin
    Write_Corrector_Information(standard_output,x,y,err,res,det);
  end Write_Corrector_Information;

  procedure Write_Corrector_Information_Data
              ( file : in file_type;
                x,y : in Double_Double_Vectors.Vector;
                err,res : in double_double ) is

  -- DESCRIPTION :
  --   Internal auxiliary procedure to implement the procedures
  --   Write_Corrector_Information with optional determinant.

  begin
    put_line(file,"Corrector information : ");
    for i in x'range loop
      put(file,"  ");
      Symbol_Table_io.put(file,Symbol_Table.get(natural32(i)));
      put(file," : ");
      put(file,x(i));
      if i < x'last then
        put(file," | y("); put(file,i,1); put(file,") = ");
        put(file,y(i)); new_line(file);
      else
        put(file," | err ="); put(file,err,2);
        put(file," res ="); put(file,res,2);
      end if;
    end loop;
  end Write_Corrector_Information_Data;

  procedure Write_Corrector_Information
              ( file : in file_type;
                x,y : in Double_Double_Vectors.Vector;
                err,res : in double_double ) is
  begin
    Write_Corrector_Information_Data(file,x,y,err,res);
    new_line(file);
  end Write_Corrector_Information;

  procedure Write_Corrector_Information
              ( file : in file_type;
                x,y : in Double_Double_Vectors.Vector;
                err,res,det : in double_double ) is
  begin
    Write_Corrector_Information_Data(file,x,y,err,res);
    put(file," det ="); put(file,det,2); new_line(file);
  end Write_Corrector_Information;

  procedure Write_Corrector_Information
              ( x : in DoblDobl_Complex_Vectors.Vector;
                err,res : in double_double ) is
  begin
    Write_Corrector_Information(standard_output,x,err,res);
  end Write_Corrector_Information;

  procedure Write_Corrector_Information_Data
              ( file : in file_type;
                x : in DoblDobl_Complex_Vectors.Vector;
                err,res : in double_double ) is

  -- DESCRIPTION :
  --   Internal auxiliary procedure for Write_Corrector_Information
  --   to implement for an optional determinant.
  
  begin
    put_line(file,"Corrector information : ");
    for i in x'range loop
      put(file,"  ");
      Symbol_Table_io.put(file,Symbol_Table.get(natural32(i)));
      put(file," : ");
      put(file,x(i)); new_line(file);
    end loop;
    put(file,"  err = "); put(file,err,3);
    put(file,"  res = "); put(file,res,3);
  end Write_Corrector_Information_Data;

  procedure Write_Corrector_Information
              ( file : in file_type;
                x : in DoblDobl_Complex_Vectors.Vector;
                err,res : in double_double ) is
  begin
    Write_Corrector_Information_Data(file,x,err,res);
    new_line(file);
  end Write_Corrector_Information;

  procedure Write_Corrector_Information
              ( x : in DoblDobl_Complex_Vectors.Vector;
                err,res,det : in double_double ) is
  begin
    Write_Corrector_Information(standard_output,x,err,res,det);
  end Write_Corrector_Information;

  procedure Write_Corrector_Information
              ( file : in file_type;
                x : in DoblDobl_Complex_Vectors.Vector;
                err,res,det : in double_double ) is
  begin
    Write_Corrector_Information_Data(file,x,err,res);
    put(file,"  det = "); put(file,det,3); new_line(file);
  end Write_Corrector_Information;

  function Maximal_Imaginary_Part
              ( v : DoblDobl_Complex_Vectors.Vector ) return double_double is

  -- DESCRIPTION :
  --   Returns the maximal imaginary part of the vector v.

    res : double_double := abs(IMAG_PART(v(v'first)));

  begin
    for i in v'first+1..v'last loop
      if abs(IMAG_PART(v(i))) > res
       then res := abs(IMAG_PART(v(i)));
      end if;
    end loop;
    return res;
  end Maximal_Imaginary_Part;

  procedure Write_Sweep_Summary
              ( file : in file_type; sols : in Solution_List;
                tol : in double_double;
                mint : out double_double; nbreal : out natural32 ) is

    tmp : Solution_List := sols;
    nbs : constant natural32 := Length_Of(sols);
    ls : Link_to_Solution;
    ima : double_double;

  begin
    mint := create(1.0E+8); nbreal := 0;
    new_line(file);
    put(file,"SWEEP SUMMARY for ");
    put(file,nbs,1); put_line(file," solutions :");
    put_line(file,"    : end value t : max imag part : real ?");
    for i in 1..nbs loop
      ls := Head_Of(tmp);
      put(file,i,3); put(file," : ");
      put(file,REAL_PART(ls.t),3); 
      if REAL_PART(ls.t) < mint
       then mint := REAL_PART(ls.t);
      end if;
      ima := Maximal_Imaginary_Part(ls.v);
      put(file," : "); put(file,ima,3);
      if ima > tol
       then put_line(file,"   : imaginary");
       else put_line(file,"   : real");
            nbreal := nbreal + 1;
      end if;
      tmp := Tail_Of(tmp);
    end loop;
    put(file,"Minimal t value : "); put(file,mint); new_line(file);
    put(file,"Number of real solutions : "); put(file,nbreal,1);
    new_line(file);
  end Write_Sweep_Summary;

  procedure Write_Sweep_Summary
              ( file : in file_type; sols : in Solution_List ) is

    tol : constant double_double := create(1.0E-8);
    min : double_double;
    nbreal : natural32;

  begin
    Write_Sweep_Summary(file,sols,tol,min,nbreal);
  end Write_Sweep_Summary;

end DoblDobl_Quad_Turn_Points_io;
