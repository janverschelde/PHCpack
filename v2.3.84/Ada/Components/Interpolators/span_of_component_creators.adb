package body Span_of_Component_Creators is

  procedure Create_Span
               ( spt : in Standard_Sample; tol : in double_float;
                 sps,sps_last : in out Standard_Sample_List;
                 sp : out Standard_Span ) is

    n : constant natural32 := natural32(Number_of_Variables(spt));

  begin
    Sample(spt,3,sps,sps_last);
    sp := Create(sps,tol);
    while Empty(sp) and (Length_Of(sps) <= n) loop
      Sample(spt,1,sps,sps_last);
      sp := Create(sps,tol);
    end loop;
  end Create_Span;

  procedure Create_Span
               ( spt : in Standard_Sample; tol : in double_float;
                 size : in natural32;
                 sps,sps_last : in out Multprec_Sample_List;
                 sp : out Multprec_Span ) is

    n : constant natural32 := natural32(Number_of_Variables(spt));

  begin
    Sample(spt,3,sps,sps_last);
    sp := Create(sps,size,tol);
    while Empty(sp) and (Length_Of(sps) <= n) loop
      Sample(spt,1,sps,sps_last);
      sp := Create(sps,size,tol);
    end loop;
  end Create_Span;

end Span_of_Component_Creators;
