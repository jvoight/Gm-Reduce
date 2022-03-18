load "config.m";
AttachSpec("../Gm-Reduce/spec");
AttachSpec("../Belyi/spec");
s := BelyiDBRead(filename);
if BelyiMapIsComputed(s) then
  print BelyiDBToRows(s);
end if;
