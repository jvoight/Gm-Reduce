load "../../Gm-Reduce/Examples/tests/test_flight.m";

/*filenames := BelyiDBFilenames(1);
filenames cat:= BelyiDBFilenames(2);
filenames cat:= BelyiDBFilenames(3);
filenames cat:= BelyiDBFilenames(4);
filenames cat:= BelyiDBFilenames(5);
filenames cat:= BelyiDBFilenames(6);
filenames cat:= BelyiDBFilenames(7);
filenames cat:= BelyiDBFilenames(8);
filenames cat:= BelyiDBFilenames(9);
*/


filenames:=BelyiDBFilenames(4);

for filename in filenames do
  ReductionTest(filename: degree_extra:=2, effort:=20);
end for;


filename:="4T1-[4,4,1]-4-4-1111-g0.m";
