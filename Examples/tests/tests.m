load "../Gm-Reduce/Examples/tests/test_flight.m";

//filenames := BelyiDBFilenames(1);
//filenames cat:= BelyiDBFilenames(2);
filenames := BelyiDBFilenames(3);
filenames cat:= BelyiDBFilenames(4);
filenames cat:= BelyiDBFilenames(5);
filenames cat:= BelyiDBFilenames(6);
filenames cat:= BelyiDBFilenames(7);
filenames cat:= BelyiDBFilenames(8);
filenames cat:= BelyiDBFilenames(9);



for filename in filenames do
  X := BelyiDBRead(filename)`BelyiDBBelyiCurves[1];
  if Genus(X) in [0,1] then
    ReductionTest(filename: degree_extra:=1, effort:=10);
  end if;
end for;
