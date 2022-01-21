//must be in BelyiDB directory
AttachSpec("code/spec_database");
load "../../Gm Reduce/Code/reducecurve.m";

/*
filenames:=[
"6T15-[5,5,4]-51-51-42-g1.m", unit group in magma not same as lmfdb??
"7T7-[5,12,12]-511-43-43-g1.m",
"5T5-[5,4,4]-5-41-41-g1.m",
"5T4-[5,5,2]-5-5-221-g1.m",
"4T5-[4,4,3]-4-4-31-g1.m",
"6T16-[6,4,6]-6-42-321-g1.m",
"7T6-[7,4,4]-7-421-421-g1.m",
"9T32-[9,3,7]-9-33111-711-g1.m",
"7T7-[6,10,4]-61-52-421-g1.m",
"7T7-[5,6,12]-511-3211-43-g0.m"
];
*/

filenames := BelyiDBFilenames(1);
filenames cat:= BelyiDBFilenames(2);
filenames cat:= BelyiDBFilenames(3);
filenames cat:= BelyiDBFilenames(4);
filenames cat:= BelyiDBFilenames(5);
filenames cat:= BelyiDBFilenames(6);
filenames cat:= BelyiDBFilenames(7);
filenames cat:= BelyiDBFilenames(8);
filenames cat:= BelyiDBFilenames(9);

filenames := [ file : file in filenames |
BelyiMapIsComputed(BelyiDBRead(file)) and BelyiDBRead(file)`BelyiDBGenus eq 1
  and Degree(BelyiDBRead(file)`BelyiDBBaseFieldData[1,1]) eq 2 ];

curves_n_maps:= < <BelyiDBRead(filename)`BelyiDBBelyiCurves[1],
                    BelyiDBRead(filename)`BelyiDBBelyiMaps[1]> : filename in filenames >;

//Index(filenames, "9T9-[4,4,2]-441-441-22221-g0.m");
for item in curves_n_maps do

  X := item[1];
  phi := item[2];

  //print "====================================";
  //filename;   X;  phi;

  RsandPs := Support(Divisor(phi));
  RsandQs := Support(Divisor(phi-1));
  PsQsRs := SetToSequence(SequenceToSet(RsandPs cat RsandQs));
  //xs := SmallFunctions(PsQsRs, 2*Genus(X)+1);
  //"The number of small functions is"; #xs;

  //SetProfile(true);
  ffs:=[];
  multiplicities:=[];
  support:=[];
  for xx in SmallFunctions(PsQsRs, 2*Genus(X)+1) do
    // this is too brutal, surely we want x to have zeros/poles at oo as well;
    // is there any pattern in what works best here?
    //try 1/phi etc
    S3orbit:=[ phi, 1/phi, phi-1, 1/(phi-1), 1/phi -1 ];
    for belyimap in [S3orbit[1]] do
      f := model(belyimap, xx);

      Append(~ffs,<#Sprint(f),f, Index(S3orbit, belyimap), Degree(xx) >);
      sup,mult:=Support(Divisor(xx));
      Append(~support, mult);
    end for;
  end for;
  ParallelSort(~ffs,~support);

  shortest_ffs:=[ ffs[i] : i in [1..Min(#ffs,5)] ];
  shortest_ffs;
  //[ support[i] : i in [1..5] ];

  freds := [];
  for ff in shortest_ffs do
    Append(~freds, reducemodel(ff[2]));
  end for;



      // sorry that these are in terms of variables u,v; this helped me to get
      // unconfused about something (fixed a stupid bug),
      // but ultimately it's more meaningful to write this in terms of x,t
      sup,mult:=Support(Divisor(xx));

      //print "f is"; print f; print "f reduced is";  print fred; print "x is"; print x; print Degree(x); print "==++==++==";
      //Append(~ffs, < 1,#Sprint(f),f,x,Degree(x) >);
      Append(~freds, <#Sprint(fred), fred, f,xx,Degree(xx)>);
      //if Max([ Abs(c) : c in mult ]) eq 1 then
        //Append(~data_mult1, [#Sprint(f), #Sprint(fred)]);
      //else
        //Append(~data_mult2, [#Sprint(f), #Sprint(fred)]);
      //end if;

      //Append(~multiplicities,<#Sprint(f), mult, x>);
      Append(~support,sup cat mult);
  end for;

ParallelSort(~freds, ~support);
//freds;

filename; X; phi;
for i in [1..1] do
  print "f is"; print freds[i,3]; print "f reduced is";
  print freds[i,2]; print "x is"; print freds[i,4]; print "Degree of x is";
  print freds[i,5]; print "support is"; print support[i];
  print "==++==++==";
end for;

/*
Sprintf("pts1 = points( %o, color='darkgreen', pointsize=50)",  [ Sprintf("(%o, %o)", data[1], data[2]) : data in data_mult1 ]);
Sprintf("pts2 = points( %o, color='magenta', pointsize=50)",  [ Sprintf("(%o, %o)", data[1], data[2]) : data in data_mult2 ]);
print "show(pts1+pts2,aspect_ratio=1)";
" ";
*/


//ParallelSort(~multiplicities,~support);
//[ <multiplicities[i], support[i]> : i in [1..#multiplicities] ];

//ProfilePrintByTotalTime(:Max:=40);

//print "f is"; print freds[1,3]; print "f reduced is";  print freds[1,2]; print "x is"; print freds[1,4];


end for;
