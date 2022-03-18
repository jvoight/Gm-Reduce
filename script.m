load "config.m";
AttachSpec("../Gm-Reduce/spec");
//AttachSpec("../Belyi/Code/spec");
s := BelyiDBRead(filename);
if BelyiMapIsComputed(s) then
  print BelyiDBToRows(s);
end if;
quit;
