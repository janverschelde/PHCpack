with text_io;                           use text_io;
with Communications_with_User;          use Communications_with_User;
with Standard_Natural_Numbers_io;       use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;         use Standard_Floating_Numbers;
with Standard_Floating_Numbers_io;      use Standard_Floating_Numbers_io;
with Standard_Complex_Norms_Equals;     use Standard_Complex_Norms_Equals;
with Standard_Complex_VecVecs;
with Standard_Complex_Polynomials;      use Standard_Complex_Polynomials;
with Multprec_Complex_Polynomials;      use Multprec_Complex_Polynomials;
with Standard_Complex_Poly_SysFun;      use Standard_Complex_Poly_SysFun;
with Multprec_Floating_Numbers;         use Multprec_Floating_Numbers;
with Multprec_Complex_Vectors;
with Multprec_Complex_Norms_Equals;     use Multprec_Complex_Norms_Equals;
with Multprec_Complex_Poly_Systems;
with Multprec_Complex_Poly_SysFun;      use Multprec_Complex_Poly_SysFun;
with Witness_Sets,Witness_Sets_io;      use Witness_Sets,Witness_Sets_io;
with Sampling_Machine;
with Sample_Points;                     use Sample_Points;
with Sample_Point_Lists;                use Sample_Point_Lists;
with Span_of_Component;                 use Span_of_Component;
with Span_of_Component_io;              use Span_of_Component_io;
with Span_of_Component_Creators;        use Span_of_Component_Creators;

package body Span_of_Witness_Set is

  procedure Test_Standard_Eval
              ( sp : in Standard_Span; tol : in double_float;
                samples : in Standard_Sample_List ) is

  -- DESCRIPTION :
  --   Tests whether all samples belong to the span or not.

    tmp : Standard_Sample_List := samples;
    spt : Standard_Sample;
    equ : constant Standard_Complex_Poly_Systems.Poly_Sys := Equations(sp);
    eva : double_float;

  begin
    for i in equ'range loop
      put("equ("); put(i,1); put(") has ");
      put(Number_of_Unknowns(equ(i)),1); put_line(" unknowns.");
    end loop;
    while not Is_Null(tmp) loop
      spt := Head_Of(tmp);
      eva := Max_Norm(Eval(equ,Sample_Point(spt).v));
      put("Residual : "); put(eva);
      if In_Span(sp,tol,Sample_Point(spt).v)
       then put_line(" belongs to span.");
       else put_line(" does not belong to span.");
      end if;
      tmp := Tail_Of(tmp);
    end loop;
  end Test_Standard_Eval;

  procedure Test_Multprec_Eval
              ( sp : in Multprec_Span; tol : in double_float;
                samples : in Multprec_Sample_List ) is

  -- DESCRIPTION :
  --   Tests whether all samples belong to the span or not.

    tmp : Multprec_Sample_List := samples;
    spt : Multprec_Sample;
    equ : constant Multprec_Complex_Poly_Systems.Poly_Sys := Equations(sp);

  begin
    while not Is_Null(tmp) loop
      spt := Head_Of(tmp);
      declare
        eva : Multprec_Complex_Vectors.Vector(equ'range)
            := Eval(equ,Sample_Point(spt).v); 
        val : Floating_Number := Max_Norm(eva);
        fltval : constant double_float := Trunc(val);
      begin
        put("Residual : "); put(fltval);
        if In_Span(sp,tol,Sample_Point(spt).v)
         then put_line(" belongs to span.");
         else put_line(" does not belong to span.");
        end if;
        Multprec_Complex_Vectors.Clear(eva);
        Clear(val);
      end;
      tmp := Tail_Of(tmp);
    end loop;
  end Test_Multprec_Eval;

  procedure Standard_Enumerate_Linear_Spans
              ( p : in Standard_Complex_Poly_Systems.Poly_Sys;
                sols : in Standard_Complex_Solutions.Solution_List;
                dim : in natural32 )  is

    file : file_type;
    sp : Standard_Span;
    spt : Standard_Sample;
    hyp : constant Standard_Complex_VecVecs.VecVec(1..integer32(dim))
        := Slices(p,dim);
    tol : constant double_float := 1.0E-8;
    sps,sps_last : Standard_Sample_List;
    use Standard_Complex_Solutions;
    tmp : Solution_List := sols;

  begin
    new_line;
    put_line("Reading the name of the output file ...");
    Read_Name_and_Create_File(file);
    Sampling_Machine.Initialize(p);
    Sampling_Machine.Default_Tune_Sampler(0);
    Sampling_Machine.Default_Tune_Refiner;
    new_line;
    for i in 1..Length_Of(sols) loop
      spt := Create(Head_Of(tmp).all,hyp);
      Create_Span(spt,tol,sps,sps_last,sp);
      put("The span at point "); put(i,1);
      put(file,"The span at point "); put(file,i,1);
      if Empty(sp) then
        put_line(" is empty.");
        put_line(file," is empty.");
      else
        put_line(" is : "); put(sp);
        put_line(file," is : "); put(file,sp);
        Test_Standard_Eval(sp,tol,sps);
      end if;
      tmp := Tail_Of(tmp);
      Shallow_Clear(spt);
      Deep_Clear(sps);
      Clear(sp);
    end loop;
    Sampling_Machine.Clear;
  end Standard_Enumerate_Linear_Spans;

  procedure Multprec_Enumerate_Linear_Spans
              ( p : in Standard_Complex_Poly_Systems.Poly_Sys;
                sols : in Standard_Complex_Solutions.Solution_List;
                dim : in natural32 )  is

    file : file_type;
    sp : Multprec_Span;
    mp : Multprec_Complex_Poly_Systems.Link_to_Poly_Sys;
    hyp : constant Standard_Complex_VecVecs.VecVec(1..integer32(dim))
        := Slices(p,dim);
    spt : Standard_Sample;
    tol : constant double_float := 1.0E-8;
    sps,sps_last : Multprec_Sample_List;
    size,deci : natural32 := 0;
    use Standard_Complex_Solutions;
    tmp : Solution_List := sols;

  begin
    new_line;
    put("Give the number of decimal places : "); get(deci);
    size := Multprec_Floating_Numbers.Decimal_to_Size(deci);
    Get_Multprec_System(p,mp,size,dim);
    new_line;
    put_line("Reading the name of the output file ...");
    Read_Name_and_Create_File(file);
    Sampling_Machine.Initialize(p,mp.all,integer32(dim),size);
    Sampling_Machine.Default_Tune_Sampler(0);
    Sampling_Machine.Default_Tune_Refiner(size);
    Sampling_Machine.Default_Tune_Refiner;
    new_line;
    for i in 1..Length_Of(sols) loop
      spt := Create(Head_Of(tmp).all,hyp);
      Create_Span(spt,tol,size,sps,sps_last,sp);
      put("The span at point "); put(i,1);
      put(file,"The span at point "); put(file,i,1);
      if Empty(sp) then
        put_line(" is empty.");
        put_line(file," is empty.");
      else
        put_line(" is : "); put(sp);
        put_line(file," is : "); put(file,sp);
        Test_Multprec_Eval(sp,tol,sps);
      end if;
      tmp := Tail_Of(tmp);
      Shallow_Clear(spt);
      Deep_Clear(sps);
      Clear(sp);
    end loop;
    Sampling_Machine.Clear;
  end Multprec_Enumerate_Linear_Spans;

end Span_of_Witness_Set;
