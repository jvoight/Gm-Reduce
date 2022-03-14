//must be in Belyi/BelyiDB directory
AttachSpec("../Code/spec");
AttachSpec("../../Gm-Reduce/spec");
SetDebugOnError(true);

/*Things to test:
-If IgnorePrimes() is removing the right primes
-The timings
-The index of the best small function
-The best degree
*/

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



WriteTestToFile:=function(filename,dextra,small_functions_size)
  //dextra is how far past d_init we go. small_functions_size is number of small functions to loop over
  writeto:="../../Gm-Reduce/test_results.m";
  try
    t:=Realtime();
    X := BelyiDBRead(filename)`BelyiDBBelyiCurves[1];
    phi := BelyiDBRead(filename)`BelyiDBBelyiMaps[1];

    RsandPs := Support(Divisor(phi));
    RsandQs := Support(Divisor(phi-1));
    PsQsRs := SetToSequence(SequenceToSet(RsandPs cat RsandQs));
    d_init:=Floor((Genus(X)+3)/2);

    fred_best:=1;
    ffs:=[];
    no_models:=small_functions_size;
    for deg in [0..d_init+dextra] do
      deg;
      xs := SmallFunctions(PsQsRs, deg);
      xs_sorted := SortSmallFunctions(phi,xs);
      xs_sorted := [ xs_sorted[i] : i in [1..Min(#xs_sorted,no_models)] ];

      ff_ds:=[];
      for x_op in xs_sorted do
        fred:=ReducedModelS3Orbit(phi,x_op);
        Append(~ff_ds,<#Sprint(fred),fred>);
      end for;

      Append(~ffs,ff_ds);
    end for;

    ff_all:=&cat(ffs);
    ffs_sorted:=Sort(ff_all);
    deg_best:= [ i : i in [1..#ffs] | ffs_sorted[1] in ffs[i] ][1]-1;

    Write(writeto, "====================================");
    Write(writeto,filename);
    Write(writeto, Sprintf("Original curve: \n %o \n", X));
    Write(writeto, Sprintf("Original map: \n %o \n",  phi));
    Write(writeto, Sprintf("Genus(X) = %o", Genus(X)));
    Write(writeto, Sprintf("[ Degree, Index of best model out of %o (or less), best model]", no_models));
    for i in [1..#ffs] do
      if ffs[i] ne [] then
        fdsort:=Sort(ffs[i]);
        fdsort_best:=fdsort[1];
        Write(writeto, [i-1, Index(ffs[i],fdsort_best), fdsort_best[2]]);
      end if;
    end for;

    Write(writeto, Sprintf("Ceiling((Genus(X)+3)/2) = %o \n", Ceiling((Genus(X)+3)/2)));
    Write(writeto, Sprintf("The degree of the best model is: %o \n", deg_best));

    //printf "The degree of the best function is \n %o = %o + %o = Ceiling((Genus(X)+3)/2) + %o", dextra_best+d_init, d_init, dextra_best, dextra_best;
    //printf "The best small function in SortSmallFunctions() has index: %o", Index(ffs[dextra_best+1],ffs_sorted[1]);
    Write(writeto, Sprintf("The number of models that were reduced is %o \n", #ff_all));
    Write(writeto, Sprintf("The total time this took was %o \n", Realtime(t)));

  catch e
    Write(writeto,e);
  end try;
  return "";
end function;

filenames:=BelyiDBFilenames(6);

for filename in filenames do
  WriteTestToFile(filename,3,20);
end for;


filename:="4T1-[4,4,1]-4-4-1111-g0.m";


/*

    dextra; ffs_sorted[1,2];
    if ffs_sorted[1,2] eq fred_best then
      break;
    else
      fred_best:=ffs_sorted[1,2];
    end if;
  end for;*/

/*    //fp,solns:=reducemodel_padic_old(fmod : Integral:=true, ClearDenominators:=false);
      fmod:=model(phi,x_op);
      ss_init:=CoefficientSupport(fmod);
      for i in [1..#ss_init] do
        cvals := [ Valuation(c,ss_init[i]) : c in Coefficients(fmod)  ];
        ignore:=(Set(cvals) in [{0,1},{0}]) and (#[ a : a in cvals | a eq 1 ] in [0,1]);
        padic_solns:=padic_LPsolutions(fmod,ss_init[i]);
        //i; Norm(ss_init[i]); cvals; ignore; padic_solns;
        if ignore then
          assert [0,0,0] in padic_solns;
        end if;
      end for;
      print "basic ignore prime worked";
*/
