AttachSpec("code/spec_database");
AttachSpec("code/spec_triangle");
AttachSpec("../Gm-Reduce/spec");
//AttachSpec("../Belyi/Code/spec");
SetProfile(true);

s := BelyiDBRead(filename);
Ks := [el[1] : el in s`BelyiDBBaseFieldData];
d := Max([Degree(el) : el in Ks]);
if BelyiMapIsComputed(s) and (d le 20) then
  try
    print BelyiDBToRows(s);
  catch e
    print "================================";
    printf "%o failed\n", s`BelyiDBFilename;
    print e`Object;
    print "================================";
  end try;
end if;
G := ProfileGraph();
ProfilePrintByTotalTime(G : Max := 15);
quit;
