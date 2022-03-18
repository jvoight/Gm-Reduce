//must be in BelyiDB directory
AttachSpec("../Gm-Reduce/spec");
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



ReductionTest := function(filename : degree_extra:=1, effort:=10)
  //{Take the Belyi map and output where in the sorted functions the best reduced map is. Also output the time. degree_extra is how far past Floor((Genus(X)+3)/2) we go. effort is how many small functions to look at for each degree}
  //writeto:="../../Gm-Reduce/test_results.m";

    t:=Realtime();
    X := BelyiDBRead(filename)`BelyiDBBelyiCurves[1];
    phi := BelyiDBRead(filename)`BelyiDBBelyiMaps[1];

    RsandPs := Support(Divisor(phi));
    RsandQs := Support(Divisor(phi-1));
    PsQsRs := SetToSequence(SequenceToSet(RsandPs cat RsandQs));
    degree_init:=Floor((Genus(X)+3)/2);

    fred_best:=1;
    ffs:=[];
    for deg in [0..degree_init+degree_extra] do
      xs := SmallFunctionsExactDegree(PsQsRs, deg);
      xs_sorted := SortSmallFunctions(phi,xs : effort := effort);
      //xs_sorted := [ xs_sorted[i] : i in [1..Min(#xs_sorted,effort)] ];

      ff_ds:=[];
      for tx in xs_sorted do
        fred:=ReducedModel(tx[1],tx[2]);
        Append(~ff_ds,<#Sprint(fred),fred>);
      end for;

      Append(~ffs,ff_ds);
    end for;

    ff_all:=&cat(ffs);
    ffs_sorted:=Sort(ff_all);
    deg_best:= [ i : i in [1..#ffs] | ffs_sorted[1] in ffs[i] ][1]-1;

    "====================================";
    filename;
    Sprintf("Original curve: \n %o \n", X);
    Sprintf("Original map: \n %o \n",  phi);
    Sprintf("Genus(X) = %o", Genus(X));
    Sprintf("[ Degree, Index of best model out of %o (or less), best model]", effort);
    for i in [1..#ffs] do
      if ffs[i] ne [] then
        fdsort:=Sort(ffs[i]);
        fdsort_best:=fdsort[1];
        [i-1, Index(ffs[i],fdsort_best), fdsort_best[2]];
      end if;
    end for;

    Sprintf("Floor((Genus(X)+3)/2) = %o \n", Floor((Genus(X)+3)/2));
    Sprintf("The degree of the best model is: %o \n", deg_best);

    //printf "The degree of the best function is \n %o = %o + %o = Ceiling((Genus(X)+3)/2) + %o", dextra_best+d_init, d_init, dextra_best, dextra_best;
    //printf "The best small function in SortSmallFunctions() has index: %o", Index(ffs[dextra_best+1],ffs_sorted[1]);
    Sprintf("The number of models that were reduced is %o \n", #ff_all);
    Sprintf("The total time this took was %o \n", Realtime(t));


  return "";
end function;


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
