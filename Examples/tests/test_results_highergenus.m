Magma V2.26-10    Fri Mar 18 2022 13:45:45 on LEGENDRE [Seed = 458069464]

+-------------------------------------------------------------------+
|       This copy of Magma has been made available through a        |
|                   generous initiative of the                      |
|                                                                   |
|                         Simons Foundation                         |
|                                                                   |
| covering U.S. Colleges, Universities, Nonprofit Research entities,|
|               and their students, faculty, and staff              |
+-------------------------------------------------------------------+

Type ? for help.  Type <Ctrl>-D to quit.
Loading "../Gm-Reduce/Examples/tests/test_flight.m"

In file "../Gm-Reduce/Examples/tests/test_flight.m", line 30, column 10:
>>     X := BelyiDBRead(filename)`BelyiDBBelyiCurves[1];
            ^
User error: Identifier 'BelyiDBRead' has not been declared or assigned

>> filenames := BelyiDBFilenames(3);
                ^
User error: Identifier 'BelyiDBFilenames' has not been declared or assigned

>> filenames cat:= BelyiDBFilenames(4);
                   ^
User error: Identifier 'BelyiDBFilenames' has not been declared or assigned

>> filenames cat:= BelyiDBFilenames(5);
                   ^
User error: Identifier 'BelyiDBFilenames' has not been declared or assigned

>> filenames cat:= BelyiDBFilenames(6);
                   ^
User error: Identifier 'BelyiDBFilenames' has not been declared or assigned

>> filenames cat:= BelyiDBFilenames(7);
                   ^
User error: Identifier 'BelyiDBFilenames' has not been declared or assigned

>> filenames cat:= BelyiDBFilenames(8);
                   ^
User error: Identifier 'BelyiDBFilenames' has not been declared or assigned

>> filenames cat:= BelyiDBFilenames(9);
                   ^
User error: Identifier 'BelyiDBFilenames' has not been declared or assigned

>>   if Genus(X) in gt 1 then
                    ^
User error: bad syntax

>>     ReductionTest(filename: degree_extra:=0, effort:=5);
                     ^
User error: Identifier 'filename' has not been declared or assigned

>>   end if;
     ^
User error: bad syntax

>> end for;
   ^
User error: bad syntax

Total time: 0.030 seconds, Total memory usage: 32.09MB
