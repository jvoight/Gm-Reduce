intrinsic model(phi::FldFunFracSchElt, x_op::FldFunFracSchElt) -> RngMPolElt
  {Given a Belyi map phi and a rational function x_op, find a plane model for the curve with phi and x_op as coordinates}
  fu := MinimalPolynomial(phi);
  fv := MinimalPolynomial(x_op);
  K := BaseRing(BaseRing(Parent(phi)));
  Kuvz<vz,uz,z> := PolynomialRing(K, 3);
  _<vf,uf,zf> := FieldOfFractions(Kuvz);

  fues := Eltseq(fu);
  fves := Eltseq(fv);
  fu := Numerator(&+[Evaluate(fues[i],zf)*uf^(i-1) : i in [1..#fues]]);
  fv := Numerator(&+[Evaluate(fves[i],zf)*vf^(i-1) : i in [1..#fves]]);
  fuv := Resultant(fu,fv,z);

  Kphi := Parent(phi);
  Kx_op := Parent(x_op);
  X := Curve(Kphi);
  // JV see monomials.m
  // if Genus(X) eq 0 then
  try
    iso_x := hom< Kx_op -> Kphi | [Kphi.1, Kphi.2]>;
  catch e;
    iso_x := hom< Kx_op -> Kphi | [Kphi.1]>;
  end try;

  fuvFact := Factorization(fuv);
  if #fuvFact gt 1 then
    for j := 1 to #fuvFact do
      if Evaluate(fuvFact[j][1], [iso_x(x_op),phi,0]) eq 0 then
        fuv := fuvFact[j][1];
        assert &and[Evaluate(fuvFact[k][1], [iso_x(x_op),phi,0]) ne 0 : k in [j+1..#fuvFact]];
      end if;
    end for;
  else
    fuv:=fuvFact[1][1];
  end if;

  _<t,x> := PolynomialRing(K,2);
  f_plane:=Evaluate(fuv,[x,t,0]);
  return f_plane;
  //x=v, t=u
end intrinsic;

intrinsic ReducedEquation(f::RngMPolElt) -> RngMPolElt
  {Given a multivariate polynomial return its reduction}
  t0:=Cputime();
  print "Starting p-adic reduction";
  f_padic, scalars1  := reducemodel_padic(f);
  t1:=Cputime();
  printf "Done with p-adic, it took %o seconds\n", t1-t0;

  t0:=Cputime();
  print "Starting unit reduction";
  f_unit, scalars2 := reducemodel_units(f_padic);
  t1:=Cputime();
  printf "Done with units, it took %o seconds\n", t1-t0;
  return f_unit, [ scalars1[i]*scalars2[i] : i in [1..#scalars1] ];
end intrinsic;

intrinsic ReducedModel(phi::FldFunFracSchElt, x_op::FldFunFracSchElt) -> RngMPolElt
  {}
  f_plane := model(phi,x_op);
  f_reduced, scalars := ReducedEquation(f_plane);
  return f_reduced, scalars;
end intrinsic;

intrinsic ReducedModelsS3Orbit(phi::FldFunFracSchElt, x_op::FldFunFracSchElt) -> RngMPolElt
  {}
  /*
    phis := S3Orbit(phi);
    return [* ReducedModel(el,x_op) : el in phis *];
  */
  f := ReducedModel(phi, x_op);
  return S3Orbit(f);
end intrinsic;

intrinsic ReducedModelS3Orbit(phi::FldFunFracSchElt, x_op::FldFunFracSchElt) -> RngMPolElt
  {find a reduced model of phi w.r.t. to x_op, compute the S3-orbit and return the smallest one}
  s3_orbit:=ReducedModelsS3Orbit(phi,x_op);
  s3_size:= [ < #Sprint(g), g > : g in s3_orbit ];
  Sort(~s3_size);
  return s3_size[1,2];
end intrinsic;

intrinsic AllReducedModels(phi::FldFunFracSchElt : effort := 0, degree := 0) -> SeqEnum
  {}

  Kinit:=BaseRing(BaseRing(Parent(phi)));
  if effort eq 0 then
    //wild effort hack
    effort:=Max(Floor(-2*(Degree(Kinit))/3+11),1);
  end if;
  if degree eq 0 then
    degree:=Floor((Genus(Curve(Parent(phi)))+3)/2);
  end if;
  RsandPs := Support(Divisor(phi));
  RsandQs := Support(Divisor(phi-1));
  PsQsRs := SetToSequence(SequenceToSet(RsandPs cat RsandQs));

  t0:=Cputime();
  print "Starting to compute SmallFunctions()";
  xs := SmallFunctions(PsQsRs, degree);
  t1:=Cputime();
  printf "Done with SmallFunctions(), it took %o seconds\n", t1-t0;

  t0:=Cputime();
  print "Starting to compute SortSmallFunctions()";
  ts_xs_Fs_sorted := SortSmallFunctions(phi, xs : effort := effort);

  while #ts_xs_Fs_sorted eq 0 do
    degree +:= 1;
    printf "degree is now %o", degree;
    xs := SmallFunctions(PsQsRs, degree);
    ts_xs_Fs_sorted := SortSmallFunctions(phi, xs : effort := effort);
  end while;
  t1:=Cputime();
  printf "Done with SortSmallFunctions(), it took %o seconds\n", t1-t0;

  printf "Computing reduced models...";
  t0 := Cputime();
  reduced_models := [];
  for tup in ts_xs_Fs_sorted do
    t, x, F := Explode(tup);
    fred,scalars := ReducedModel(t, x);
    // printf "t = %o,\nx = %o,\nreduced model = %o\n\n", t, x, fred;
    Append(~reduced_models, <#Sprint(fred), t, x, fred, scalars>);
  end for;
  t1 := Cputime();
  printf "done. That took %o seconds\n", t1 - t0;
  Sort(~reduced_models);
  // return reduced_models;
  return [ <reddat[4], reddat[5]> : reddat in reduced_models];
end intrinsic;

intrinsic BestModel(phi::FldFunFracSchElt : effort := 10, degree := 0) -> RngMPolElt
  {return then best model with some search parameters}
  if degree eq 0 then
    degree:=Floor((Genus(Curve(Parent(phi)))+3)/2);
  end if;
  list:=AllReducedModels(phi : effort:=effort, degree:=degree);
  f := list[1][1];
  return f, BaseRing(Parent(f))!1/list[1][2][1];
end intrinsic;

intrinsic PlaneModel(phi::FldFunFracSchElt, x_op::FldFunFracSchElt) -> RngMPolElt
  {}
  return model(phi,x_op);
end intrinsic;
