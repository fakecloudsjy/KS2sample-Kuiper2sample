(* Kuiper2Sample.wl
   Exact two-sample Kuiper test for data on the real line or circle,
   with support for ties.
   Dimitrova, D. S., Jia, Y., Kaishev, V. K. (2026).
   Efficient Exact Calculation of p-values of the Two-sample
   Kolmogorov-Smirnov and Kuiper Tests.
   To appear in Journal of Statistical Computation and Simulation. *)

(* Helper: Nikiforov (1994) recurrence for one pair of one-sided bounds.
   The 10^-9 offset enforces strict inequality (D+ < qup, D- < qdown),
   consistent with the Nikiforov Fortran convention. *)
kuiperks[m1_Integer, n1_Integer, qup_, qdown_, multiplicity_List] :=
  Module[
    {p, r, ni, slope, deviatUp, deviatDown, icat, nties, ile, iri,
     iles, iris, iri2, ile2, dl, newrct, l2, ileft, iright, jupp, jlow},

    newrct = True;  ni = m1 + n1;  icat = 2;
    p = Table[0, m1 + 2];  p[[1]] = 1;
    slope      = m1/ni;
    deviatUp   = slope * (qup   + 10^-9) * n1;
    deviatDown = slope * (qdown + 10^-9) * n1;
    ile = iri = ile2 = iri2 = l2 = dl = 0;
    nties = multiplicity[[1]];
    r = Prepend[p, 0];

    Do[
      r = Prepend[p, 0];
      If[nties == 1,
        dl    = L * slope;
        iri   = Min[Floor[dl + deviatDown], L, m1];
        ile   = Max[Floor[dl - deviatUp + 1], L - n1, 0];
        nties = multiplicity[[icat]];
        icat  = icat + 1;
        newrct = True
        ,
        nties = nties - 1;
        If[newrct,
          newrct = False;
          l2     = L + nties;
          dl     = l2 * slope;
          iri2   = Min[Floor[dl + deviatDown], l2, m1];
          ile2   = Max[Floor[dl - deviatUp + 1], l2 - n1, 0];
          ileft  = ile;
          iright = iri2;
          jupp   = l2 - ile2;
          jlow   = L  - iri - 1
        ];
        ile = Max[ileft, L - jupp];
        iri = Min[iright, L - jlow]
      ];

      iles = Max[0, ile];   iris = Min[L, iri];
      p[[(iles + 1) ;; (iris + 1)]] =
        r[[(iles + 2) ;; (iris + 2)]] + r[[(iles + 1) ;; (iris + 1)]];
      iles = Max[1, ile];   iris = Min[L - 1, iri];
      p[[iles]]     = If[ile == 0, 1, 0];
      p[[iris + 2]] = If[iri == L, 1, 0],
      {L, 1, ni - 1}
    ];

    p[[m1 + 1]] + p[[m1]]
  ]

(* Dispatch: no-ties form -> ties form *)
Kuiper2Sample[m_Integer, n_Integer, q_?NumericQ] :=
  Kuiper2Sample[m, n, ConstantArray[1, m + n], q]

(* Core: exact p-value via Proposition 2 of Dimitrova, Jia, Kaishev (2026) *)
Kuiper2Sample[m_Integer, n_Integer, multiplicity_List, q_?NumericQ] :=
  Module[
    {nxy, ci, pos, neg, p},

    nxy = LCM[m, n];
    ci  = Ceiling[q * nxy - 1];

    Which[
      ci <= 0, 1,
      q >= 2,  0,
      ci <= nxy - 1,
        pos = Join[
          {{m, n, ci/nxy, 0, multiplicity}},
          Table[{m, n, i/nxy, (ci - i)/nxy,     multiplicity}, {i, 0, ci - 1}]
        ];
        neg = Table[{m, n, i/nxy, (ci - i - 1)/nxy, multiplicity}, {i, 0, ci - 1}];
        1 - (Total[MapApply[kuiperks, pos]] -
             Total[MapApply[kuiperks, neg]]) / Binomial[m + n, m]
        ,
        True,
        p = 2 * nxy - ci;
        pos = Join[
          {{m, n, 1, ci/nxy - 1, multiplicity}},
          Table[{m, n, (ci - nxy + i)/nxy, (nxy - i)/nxy,     multiplicity}, {i, 0, p - 1}]
        ];
        neg = Table[{m, n, (ci - nxy + i)/nxy, (nxy - i - 1)/nxy, multiplicity}, {i, 0, p - 1}];
        1 - (Total[MapApply[kuiperks, pos]] -
             Total[MapApply[kuiperks, neg]]) / Binomial[m + n, m]
    ]
  ]

(* Data interface: computes statistic from raw samples, dispatches to core *)
Kuiper2Sample[data1_List, data2_List] :=
  Module[
    {m, n, joint, multiplicity, edf1, edf2, z, dstat, pval},

    m = Length[data1];  n = Length[data2];
    joint        = Sort[Join[data1, data2]];
    multiplicity = Values[Counts[joint]];
    edf1         = EmpiricalDistribution[data1];
    edf2         = EmpiricalDistribution[data2];

    (* V = Delta+ + Delta- at all distinct pooled-sample points *)
    joint = DeleteDuplicates[joint];
    z     = CDF[edf1, joint] - CDF[edf2, joint];
    dstat = Max[z] + Max[-z];
    pval  = Kuiper2Sample[m, n, multiplicity, dstat];

    {dstat, pval}
  ]
