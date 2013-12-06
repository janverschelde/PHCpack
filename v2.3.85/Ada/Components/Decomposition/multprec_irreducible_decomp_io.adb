with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Multprec_Complex_VecVecs_io;        use Multprec_Complex_VecVecs_io;
with Standard_Complex_Solutions;         use Standard_Complex_Solutions;
with Standard_Complex_Solutions_io;      use Standard_Complex_Solutions_io;
with Multprec_Complex_Polynomials_io;    use Multprec_Complex_Polynomials_io;
with Multprec_Complex_Poly_Systems;      use Multprec_Complex_Poly_Systems;
with Multprec_Complex_Poly_Systems_io;   use Multprec_Complex_Poly_Systems_io;
with Planes_and_Polynomials;             use Planes_and_Polynomials;

package body Multprec_Irreducible_Decomp_io is

  procedure put ( h : in Hypersurfaces ) is
  begin
    put(Standard_Output,h);
  end put;

  procedure put ( file : in file_type; h : in Hypersurfaces ) is

    lenpts : constant natural32 := Length_Of(h.pts);

  begin
    put_line(file,"DIMENSIONS :");
    put(file,"  ");
    put(file,h.n,1);
    put(file,"D work space, ");
    if h.d = 0 then
      if lenpts = 0 then
        put_line(file,"no isolated solutions.");
      else
        put(file,lenpts,1);
        put(file," isolated solution");
        if lenpts = 1
         then put_line(file,".");
         else put_line(file,"s.");
        end if;
      end if;
      if h.equ /= null then
        put_line(file,"POLYNOMIAL SYSTEM : ");
        put(file,h.equ.all);
      end if;
      if lenpts > 0 then
        put_line(file,"THE ISOLATED SOLUTIONS : ");
        put(file,lenpts,natural32(Head_Of(h.pts).n),h.pts);
      end if;
    else
      if h.equ = null then
        put(file,"no components of dimension ");
        put(file,h.d,1);
        put_line(file,".");
      else
        put(file,natural32(h.equ'length),1);
        put(file," component");
        if h.equ'length = 1
         then put(file," of dimension ");
         else put(file,"s of dimension ");
        end if;
        put(file,h.d,1);
        put(file,", with sum of degrees = ");
        put(file,lenpts,1); put_line(file,".");
        put_line(file,"HYPERSURFACE EQUATIONS : ");
        put_line(file,h.equ.all);
        put_line(file,"PROJECTION OPERATORS AS POLYNOMIALS :");
        for i in h.hyp'range loop
          put_line(file,Hyperplane(h.hyp(i).all));
        end loop;                                                   
      end if;
      if not Is_Null(h.pts) then
        put_line(file,"GENERIC POINTS :");
        put(file,lenpts,natural32(Head_Of(h.pts).n),h.pts);
      end if;
      put_line(file,"HYPERPLANES FOR THE PROJECTION OPERATORS :");
      put_line(file,h.hyp);
    end if;
  end put;

  procedure put ( s : in Solution_Components ) is
  begin
    put(Standard_Output,s);
  end put;
 
  procedure put ( file : in file_type; s : in Solution_Components ) is
  begin
    for i in s'range loop
      if s(i) /= null
       then put(file,s(i).all);
      end if;
    end loop;
  end put;

end Multprec_Irreducible_Decomp_io;
