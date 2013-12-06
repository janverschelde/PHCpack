with text_io;                            use text_io;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Multprec_Integer_Numbers;           use Multprec_Integer_Numbers;
with Multprec_Integer_Numbers_io;        use Multprec_Integer_Numbers_io;
with Standard_Integer_Vectors;
with Standard_Integer_Vectors_io;        use Standard_Integer_Vectors_io;
with Multprec_Integer_Vectors_io;        use Multprec_Integer_Vectors_io;
with Multprec_Lattice_Supports;

package body Multprec_Lattice_3d_Facets_io is

  procedure Write_Coordinates ( A : in Matrix; k : in integer32 ) is
  begin
    for i in A'range(1) loop
      put(" "); put(A(i,k));
    end loop;
  end Write_Coordinates;

  procedure Write_Supported_Points
               ( A : in Matrix; v : in Multprec_Integer_Vectors.Vector ) is

    s : constant Standard_Integer_Vectors.Vector
      := Multprec_Lattice_Supports.Support(A,v);
    m : Integer_Number := Multprec_Lattice_Supports.Minimum(A,v);

  begin
    put("Vector"); put(v);
    put(" supports "); put(s'last,1); put(" points"); put(s);
    put(" with minimum "); put(m,1); new_line;
    Multprec_Integer_Numbers.Clear(m);
  end Write_Supported_Points;

  procedure Write_Facet ( A : in Matrix; f : in Facet_in_3d ) is

    m : Integer_Number := Multprec_Lattice_Supports.Minimum(A,f.normal);

  begin
    put("facet "); put(f.label,1);
    put(" spanned by"); put(f.points);
    put(" has normal"); put(f.normal);
    put(" and value "); put(m); new_line;
    Clear(m);
    put("  IP :");
    put(Multprec_Lattice_Supports.Inner_Products(A,f.normal));
    put("  support :");
    put(Multprec_Lattice_Supports.Support(A,f.normal)); new_line;
    put("  neighboring facets :");
    for i in f.neighbors'range loop
      if f.neighbors(i) /= null
       then put(" "); put(f.neighbors(i).label,1);
      end if;
    end loop;
    new_line;
  end Write_Facet;

  procedure Write_Facets ( A : in Matrix; f : in Facet_3d_List ) is

    tmp : Facet_3d_List := f;
    lft : Link_to_3d_Facet;

  begin
    while not Is_Null(tmp) loop
      lft := Head_Of(tmp);
      Write_Facet(A,lft.all);
      tmp := Tail_Of(tmp);
    end loop;
  end Write_Facets;

end Multprec_Lattice_3d_Facets_io;
