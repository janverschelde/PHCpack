with Main_Pade_Trackers;

procedure mainseries ( nt : in natural32; precision : in character;
                       infilename,outfilename : in string;
                       verbose : in integer32 := 0 ) is
begin
  Main_Pade_Trackers.Main(infilename,outfilename,nt,precision,verbose);
end Mainseries;
