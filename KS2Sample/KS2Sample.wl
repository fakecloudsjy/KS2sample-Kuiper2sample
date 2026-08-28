(* KS2Sample.wl
   Exact weighted two-sample Kolmogorov-Smirnov test.
   Dimitrova, D. S., Jia, Y., Kaishev, V. K. (2026).
   Efficient Exact Calculation of p-values of the Two-sample
   Kolmogorov-Smirnov and Kuiper Tests.
   To appear in Journal of Statistical Computation and Simulation. *)

Options[KS2Sample] = {"Alternative" -> "TwoSided", "Weight" -> 0};

(* Dispatch: no-ties form -> ties form *)
KS2Sample[m_Integer, n_Integer, q_?NumericQ, opts:OptionsPattern[]] :=
  KS2Sample[m, n, ConstantArray[1, m + n], q, opts]

(* Core: Nikiforov (1994) recurrence, extended to ties and weight functions *)
KS2Sample[m_Integer, n_Integer, multiplicity_List, q_?NumericQ, opts:OptionsPattern[]] :=
  Module[
    {newrct, p, r, ni, slope, deviat, nties, ile, icat, ileft,
     iright, l2, jupp, jlow, iri, iles, iris, iri2, ile2, dl,
     delta, alternative, altCode, weight, weightFun, wvec},

    alternative = OptionValue["Alternative"];
    weight      = OptionValue["Weight"];

    If[!MemberQ[{"TwoSided", "Greater", "Less"}, alternative], Return[$Failed]];
    altCode = Switch[alternative, "TwoSided", 1, "Greater", 2, "Less", 3];

    weightFun = Which[
      NumericQ[weight] && weight == 0,
        Function[t, 1],
      NumericQ[weight] && 0 < weight <= 1,
        Function[t, (t * (1 - t))^(-weight)],
      Head[weight] === Function || Head[weight] === Symbol,
        weight,
      True,
        $Failed
    ];
    If[weightFun === $Failed, Return[$Failed]];

    ni   = m + n;
    wvec = Table[weightFun[i/ni], {i, 1, ni - 1}];

    newrct = True;  icat = 2;
    p = Table[0, m + 2];  p[[1]] = 1;
    slope = m/ni;
    delta = slope * q * n;
    ile = iri = ileft = iright = ile2 = iri2 = l2 = jupp = jlow = dl = 0;
    nties = multiplicity[[1]];
    r = Prepend[p, 0];

    Do[
      r = Prepend[p, 0];

      If[nties == 1,
        (* New distinct value: recompute active column bounds *)
        dl     = L * slope;
        deviat = delta / wvec[[L]];
        iri    = Min[Ceiling[dl + deviat] - 1, L, m];
        ile    = Max[Floor[dl - deviat + 1],   L - n, 0];
        iri    = If[altCode == 3, Min[m, L], iri];
        ile    = If[altCode == 2, Max[0, L - n], ile];
        nties  = multiplicity[[icat]];
        icat   = icat + 1;
        newrct = True
        ,
        (* Inside a tied block: carry bounds from block entry to block exit *)
        nties = nties - 1;
        If[newrct,
          newrct = False;
          l2     = L + nties;
          dl     = l2 * slope;
          If[l2 == ni,
            iri2 = m;  ile2 = m
            ,
            deviat = delta / wvec[[l2]];
            iri2   = Min[Ceiling[dl + deviat] - 1, l2, m];
            ile2   = Max[Floor[dl - deviat + 1],   l2 - n, 0]
          ];
          ileft  = ile;
          iright = iri2;
          jupp   = l2 - ile2;
          jlow   = L  - iri - 1
        ];
        ile = Max[ileft, L - jupp];
        iri = Min[iright, L - jlow]
      ];

      iri = If[altCode == 3, Min[m, L], iri];
      ile = If[altCode == 2, Max[0, L - n], ile];

      (* Recurrence B_S: count trajectories staying inside region S *)
      iles = Max[0, ile];   iris = Min[L, iri];
      p[[(iles + 1) ;; (iris + 1)]] =
        r[[(iles + 2) ;; (iris + 2)]] + r[[(iles + 1) ;; (iris + 1)]];

      (* Boundary sentinels *)
      iles = Max[1, ile];   iris = Min[L - 1, iri];
      p[[iles]]     = If[ile == 0, 1, 0];
      p[[iris + 2]] = If[iri == L, 1, 0],
      {L, 1, ni - 1}
    ];

    1 - (p[[m + 1]] + p[[m]]) / Binomial[ni, m]
  ]

(* Data interface: computes statistic from raw samples, dispatches to core *)
KS2Sample[data1_List, data2_List, OptionsPattern[]] :=
  Module[
    {m, n, ni, joint, multiplicity, wvec, edf1, edf2, z,
     dstat, pval, weightFun, alternative, weight},

    alternative = OptionValue["Alternative"];
    weight      = OptionValue["Weight"];

    m = Length[data1];  n = Length[data2];
    ni = m + n;

    If[!MemberQ[{"TwoSided", "Greater", "Less"}, alternative], Return[$Failed]];

    weightFun = Which[
      NumericQ[weight] && weight == 0,
        Function[t, 1],
      NumericQ[weight] && 0 < weight <= 1,
        Function[t, (t * (1 - t))^(-weight)],
      Head[weight] === Function || Head[weight] === Symbol,
        weight,
      True,
        $Failed
    ];
    If[weightFun === $Failed, Return[$Failed]];

    wvec         = Table[weightFun[i/ni], {i, 1, ni - 1}];
    joint        = Sort[Join[data1, data2]];
    multiplicity = Values[Counts[joint]];
    edf1         = EmpiricalDistribution[data1];
    edf2         = EmpiricalDistribution[data2];

    joint = DeleteDuplicates[joint];
    z = CDF[edf1, joint] - CDF[edf2, joint];
    z = z[[1 ;; (Length[joint] - 1)]] *
        wvec[[Accumulate[multiplicity][[1 ;; (Length[multiplicity] - 1)]]]];

    dstat = Switch[alternative,
      "TwoSided", Max[Abs[z]],
      "Greater",  Max[z],
      "Less",     Max[-z]
    ];

    pval = KS2Sample[m, n, multiplicity, dstat,
             "Alternative" -> alternative, "Weight" -> weight];
    {dstat, pval}
  ]
