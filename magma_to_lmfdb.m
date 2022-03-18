intrinsic BelyiDBToRow(s::BelyiDB) -> MonStgElt
  {}
  row := "";
  // galmaps dictionaries
  gal_orbits_before_sorting := GaloisOrbits(s); // we will sort by size (increasing)
  gal_orbits := gal_orbits_before_sorting;
  gal_orbits_sizes := [#orbit : orbit in gal_orbits_before_sorting];
  ParallelSort(~gal_orbits_sizes, ~gal_orbits);
  pass := PointedPassport(s);
  for i := 1 to #gal_orbits do
    if i lt #gal_orbits then assert #gal_orbits[i] le #gal_orbits[i+1]; end if;
    gal_orbit := gal_orbits[i];
    inds := [Index(pass, triple) : triple in gal_orbit];
    belyi_label := BelyiDB_label(s, inds, i);
    lmfdb_label := GalmapLabel(s, inds, i);
    phi := BelyiMap(s, inds, i);
    f := BestModel(phi);
    row *:= Sprintf("%o|%o|%o\n", lmfdb_label, belyi_label, f);
  end for;
  return row;
end intrinsic;
