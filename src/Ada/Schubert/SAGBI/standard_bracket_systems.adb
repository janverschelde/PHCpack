with Straightening_Syzygies;             use Straightening_Syzygies;

package body Standard_Bracket_Systems is

  function Straightening_Syzygies ( n,d : natural32 ) return Bracket_System is

    nonstd : constant Bracket_Polynomial := nonStandard_Monomials(n,d);
    res : Bracket_System(1..integer32(Number_of_Monomials(nonstd)));
    cnt : integer32 := 0;

    procedure Store_Syzygy ( t : in Bracket_Term; continue : out boolean ) is
    begin
      cnt := cnt + 1;
      res(cnt) := Straightening_Syzygy(t.monom);
      continue := true;
    end Store_Syzygy;
    procedure Store_Syzygies is new Enumerate_Terms(Store_Syzygy);

  begin
    Store_Syzygies(nonstd);
    return res;
  end Straightening_Syzygies;

  procedure Clear ( s : in out Bracket_System ) is
  begin
    for i in s'range loop
      Clear(s(i));
    end loop;
  end Clear;

end Standard_Bracket_Systems;
