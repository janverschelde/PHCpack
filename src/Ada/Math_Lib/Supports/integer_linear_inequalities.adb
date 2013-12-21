with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Floating_Vectors;          use Standard_Floating_Vectors;
with Standard_Floating_Matrices;
with Floating_Linear_Inequalities;       use Floating_Linear_Inequalities;

package body Integer_Linear_Inequalities is

  procedure Integer_Complementary_Slackness
                  ( tableau : in out Matrix; feasible : out boolean ) is
  begin
    Integer_Complementary_Slackness(tableau,tableau'last(2)-1,feasible);
  end Integer_Complementary_Slackness;

  procedure Integer_Complementary_Slackness
                  ( tableau : in out Matrix; lastcol : in integer32;
                    feasible : out boolean ) is

    tab : Standard_Floating_Matrices.Matrix
                 (tableau'range(1),tableau'first(2)..lastcol);
    rhs,sol : Standard_Floating_Vectors.Vector(tab'range(1));
    tol : constant double_float := 10.0**(-8); --10.0**(-12);
    columns : Standard_Integer_Vectors.Vector(sol'range);

  begin
    for i in tab'range(1) loop
      for j in tab'range(2) loop
        tab(i,j) := double_float(tableau(i,j));
      end loop;
    end loop;
    for i in rhs'range loop
      rhs(i) := double_float(tableau(i,tableau'last(2)));
    end loop;
    Complementary_Slackness(tab,lastcol,rhs,tol,sol,columns,feasible);
  end Integer_Complementary_Slackness;

end Integer_Linear_Inequalities;
