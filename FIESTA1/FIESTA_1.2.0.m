(*
    Copyright (C) Alexander Smirnov and Mikhail Tentyukov.
    The program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License version 2 as
    published by the Free Software Foundation.

    The program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.
*)

<<Combinatorica`;
Unprotect[Combinatorica`EmptyQ];
Remove[Combinatorica`EmptyQ];
pm;
If[Not[ValueQ[NumberOfSubkernels]],
    NumberOfSubkernels=0;
];

If[Not[ValueQ[CacheSize]],
    CacheSize=10000000
];
If[Not[ValueQ[ReturnErrorWithBrackets]],
    ReturnErrorWithBrackets=False
];

SDFile;
CurrentIntegratorSettings;
If[Not[ValueQ[d0]],
    d0=4;
];
If[Not[ValueQ[UsingC]],
    UsingC=True;
];
If[Not[ValueQ[CIntegratePath]],
    CIntegratePath="CIntegrate";
];
If[Not[ValueQ[NumberOfLinks]],
    NumberOfLinks=1;
];
If[Not[ValueQ[UsingQLink]],
    UsingQLink=True;
];
If[Not[ValueQ[QLinkPath]],
    QLinkPath="QLink64";
];
If[Not[ValueQ[DataPath]],
    DataPath="temp";
];
If[Not[ValueQ[PrepareStringsWhileIntegrating]],
    PrepareStringsWhileIntegrating=True;
];
If[Not[ValueQ[IntegrationCut]],
    IntegrationCut=0;
];
If[Not[ValueQ[IfCut]],
    IfCut=0.00;
];
If[Not[ValueQ[VarExpansionDegree]],
    VarExpansionDegree=1;
];
If[Not[ValueQ[ForceMixingOfSectors]],
    ForceMixingOfSectors=False;
];
If[Not[ValueQ[CurrentIntegrator]],
    CurrentIntegrator="vegasf";
]


STRATEGY_A;
STRATEGY_B;
STRATEGY_S;
STRATEGY_X;
STRATEGY_0;
STRATEGY_SS;
PrimarySectorCoefficients;

If[Not[ValueQ[STRATEGY]],
    STRATEGY=STRATEGY_S;
]; (*might be STRATEGY_0, STRATEGY_A, STRATEGY_B, STRATEGY_S, STRATEGY_X*)
If[Not[ValueQ[ResolveNegativeTerms]],
    ResolveNegativeTerms=True;
];
If[Not[ValueQ[UsingComplexNumbers]],
    UsingComplexNumbers=False;
];
If[Not[ValueQ[ComplexShift]],
    ComplexShift=- (I*(10^(-8)));
];
If[Not[ValueQ[AbortOnNaN]],
    AbortOnNaN=False;
];
If[Not[ValueQ[IntegrationDebug]],
    IntegrationDebug=False;
];
If[Not[ValueQ[NoIntegration]],
    NoIntegration=False;
];
If[Not[ValueQ[ExportData]],
    ExportData=False;
];
If[Not[ValueQ[ExportDataCounter]],
    ExportDataCounter=0;
];
If[Not[ValueQ[PMCounter]],
    PMCounter=1;
];
If[Not[ValueQ[VegasSettings]],
    VegasSettings={10000,5,100000,15};
]


WhyIndeterminate[]:=Module[{temp},
    RawPrintLn["Possible reasons for an Indeterminate result:"];
    RawPrintLn["1) Negative F provided."];
    RawPrintLn["This may happen in case the propagators of the form 1/(p^2-m^2)."];
    RawPrintLn["Solution: change propagators sign."];
    RawPrintLn["2) Complex number as an answer."];
    RawPrintLn["Currently the C integration cannot handle complex numbers."];
    RawPrintLn["Solution: try UsingC=False or wait for new versions."];
    RawPrintLn["3) Special singularities."];
    RawPrintLn["The standard sector decomposition approach works only for singularities for small values of variables."];
    RawPrintLn["Our code can resolve some extra singularities, but surely all types cannot be covered."];
    RawPrintLn["Solution: Try different values of the ResolveNegativeTerms option."];
    RawPrintLn["Contact the authors in order to try to make the resolution of singularities of this type automatic."];
    RawPrintLn["4) Numerical instability."];
    RawPrintLn["If all the singularities are for small values of integration variables, they are resolved in the code."];
    RawPrintLn["However, the numerical integration may experience problems with small values of variables."];
    RawPrintLn["Solution: use the IntegrationCut or IfCut together with VarExpansionDegree options."];
    RawPrintLn["IntegrationCut is faster but results in poor precision (integrating from IntegrationCut instead of 0)."];
    RawPrintLn["IfCut replaces the integrand with its expansion for small values of problematic variables."];
    RawPrintLn["The integration might be slow especially if VarExpansionDegree>1 (but even VarExpansionDegree=0 might be good enough)."];
]


Begin["FIESTA`"];

SDEvaluateG::Usage="SDEvaluate[{graph,external},{U,F,loops},indices,order] evaluates the integral"
SDEvaluate::Usage="SDEvaluate[{U,F,loops},indices,order] evaluates the integral"
ClearResults::Usage="ClearResults[] clears the results from memory"
UF::Usage="UF[LoopMomenta,Propagators,subst] generates the functions U and F"






UF[xx_, yy_, z_] := Module[{degree, coeff, i, t2, t1, t0},
  degree =
   -Sum[yy[[i]]*x[i], {i, 1, Length[yy]}];
  coeff = 1;
  For[i = 1, i <= Length[xx], i++,
   t2 = Coefficient[degree, xx[[i]], 2];
   t1 = Coefficient[degree, xx[[i]], 1];
   t0 = Coefficient[degree, xx[[i]], 0];
   coeff = coeff*t2;
   degree = Together[t0 - ((t1^2)/(4 t2))];
   ];
  degree = Together[-coeff*degree] //. z;
  coeff = Together[coeff] //. z;
  {coeff, Expand[degree], Length[xx]}
]

ClearResults[]:=Clear[SDIntegral, SDResolve];

VersionString:="FIESTA, version 1.2.0";
CodeInfo[]:=Module[{temp},
    RawPrintLn["UsingC: ",UsingC];
    RawPrintLn["NumberOfLinks: ",NumberOfLinks];
    RawPrintLn["UsingQLink: ",UsingQLink];
    RawPrintLn["IntegrationCut: ",IntegrationCut];
    RawPrintLn["IfCut: ",ToString[IfCut,InputForm]];
    If[IfCut>0,
        RawPrintLn["VarExpansionDegree: ", VarExpansionDegree]
    ];
    RawPrintLn["Strategy: ",STRATEGY];
    If[And[Not[STRATEGY===STRATEGY_SS],GraphUsed===True],
        RawPrintLn["WARNING: a graph has been specified -> STRATEGY_SS recomended"];
    ];
]
Print[VersionString];



CutExtraDigits[xx_]:=Chop[N[Round[xx,10^-6]],10^-6];
DoubleCutExtraDigits[xx_]:=Module[{temp},
    If[xx===Indeterminate,Return[Indeterminate]];
    temp=Variables[xx];
    If[Length[temp]>1,Return[xx]];
    If[Length[temp]==0,Return[CutExtraDigits[xx]]];
    CutExtraDigits[xx/.(temp[[1]]->0)]+CutExtraDigits[(xx-(xx/.(temp[[1]]->0)))/.(temp[[1]]->1)]*temp[[1]]
]

PMSymbol:=If[$FrontEnd===Null,"+-","±"]

PMSeriesForm[xx_] :=
 Module[{i},
 If[xx===Indeterminate,Return["INDETERMINATE"]];
 If[NumberQ[xx],Return[ToString[xx]]];
  StringReplace[StringReplace[(StringJoin @@ (Table[
        "(" <> PMForm[xx[[3]][[i]]] <> ")" <>
         If[xx[[4]] + (i - 1)*xx[[6]] === 0, "",
          If[xx[[4]] + (i - 1)*xx[[6]] === 1, "ep",
           If[xx[[4]] + (i - 1)*xx[[6]] > 0,
            "ep^" <> ToString[xx[[4]] + (i - 1)*xx[[6]]],
            "/ep^" <> ToString[-xx[[4]] - (i - 1)*xx[[6]]] <> ""]]] <>
          " + ", {i, Length[xx[[3]]]}])) <> "+", {"+ +" -> "",
    "^1 +" -> " +"}],{"+ +" -> ""}]]

PMForm[xx_]:=Module[{temp},
    If[xx===Indeterminate,Return["INDETERMINATE"]];
    temp=Variables[xx];
    If[Length[temp]>1,Return[ToString[xx,InputForm]]];
    If[Length[temp]==0,Return[ToString[xx,InputForm]]];
    ToString[xx/.(temp[[1]]->0),InputForm]<>" "<>PMSymbol<>" "<>ToString[(xx-(xx/.(temp[[1]]->0)))/.(temp[[1]]->1),InputForm]
]







RawPrint[x__]:=WriteString[$Output,x];
RawPrintLn[x__]:=WriteString[$Output,x,"\n"];
RawPrintLn[]:=WriteString[$Output,"\n"];
PrepareDots[n_]:=Module[{temp},gn=n;gl={};For[gi=1,gi<=10,gi++,AppendTo[gl,Quotient[gi*gn,10]]];gi=0;]
Dots[]:=(gi++;While[And[Length[gl]>0,gi>=gl[[1]]],gl=Drop[gl,1];RawPrint["."]];)

PrepareDots2[n_]:=Module[{temp},gn2=n;gl2={};For[gi2=1,gi2<=10,gi2++,AppendTo[gl2,Quotient[gi2*gn2,10]]];gi2=0;]
Dots2[]:=(gi2++;While[And[Length[gl2]>0,gi2>=gl2[[1]]],gl2=Drop[gl2,1];RawPrint["."]];)


CC[m_]:=(If[Mod[gi,m]===0,MyClearCache[]])
MyTimingForm[xx_]:=ToString[Chop[N[Round[xx,10^-4]],10^-4],InputForm]


If[Not[TrueQ[$VersionNumber>=6.0]],
  <<"Developer`";
  MyClearCache[] := Developer`ClearCache[];
  Print["WARNING! Full compatibility with Mathematica versions lower than 6.0 not guaranteed!"];
,
  MyClearCache[] := ClearSystemCache[];
];




FactorMonom[xxx_]:=Module[{xx,vars,i,monom,temp},
    xx=xxx;
    If[xx===0,Return[{0,0}]];
    vars=Sort[Variables[xx]];
    temp=-Exponent[xx,##^-1]&/@vars;
    monom=Times@@(vars^temp);
    {monom,Expand[xx/(monom)]}
]


ConstructTerm[xx_]:={Flatten[Position[xx[[1]],2]],Inner[(#1^#2)&,Array[x,Length[xx[[1]]]],xx[[2]],Times]*xx[[3]]};

EpPosExpand[xx_,order_]:= Module[{temp,shift,i,der,result},
    shift=0;temp=xx;
    While[TrueQ[MemberQ[{ComplexInfinity,Indeterminate},Quiet[(temp) /. ep -> 0]]],
        temp=Expand[temp*ep];
        shift++;
    ];
    result={{-shift,temp/.ep->0}};
    der=temp;
    For[i=1,i<=order+shift,i++,
        der=D[der,ep];
        temp=((der/Factorial[i])/.ep->0);
        If[Not[temp===0],AppendTo[result,{-shift+i,temp}]];
    ];
    result
]






MyExponent[xxx_, yy_] := Module[{temp},
    xx=Expand[xxx/(xxx/.x[aaa_]->1)];
  temp = Exponent[xx /. ep -> 0, yy];
  Expand[(Exponent[xx /. ep -> 1, yy] - temp)ep + temp]
  ]



ExpandZ[xx_] :=
 Module[{temp, degs,i, n, j,baddegs,gooddegs,monom},
  temp = {xx[[1]], xx[[2]], ##, xx[[4]]} & /@ xx[[3]];
  n = Length[temp[[1]][[1]]];

  temp = (degs=##[[2]] + Table[MyExponent[##[[3]], x[j]], {j, n}];
            {baddegs,gooddegs}=Transpose[Table[If[(degs[[j]]/.ep->0)>=0,{0,degs[[j]]},{degs[[j]],0}], {j, n}]];
            monom=Inner[Power,x/@Range[n],gooddegs,Times];
        {##[[1]], baddegs, (monom)*(##[[3]] /. x[aa_]->1),##[[4]]}

         )& /@ temp;
  temp]

GroupTerms[xx_] := Module[{temp},
   temp =
    Reap[Sow[##[[3]], {{##[[1]], ##[[2]], ##[[4]]}}] & /@ xx, _,
      List][[2]];
   temp = {##[[1]][[1]], ##[[1]][[2]], (Plus @@ ##[[2]])*##[[1]][[3]]} & /@ temp;
   temp
   ];



FirstVars[xx__]:=Module[{vars},
  vars = Union[Cases[xx,x[aaa_],Infinity]];
  {xx/.Apply[Rule,Transpose[{vars,Array[y,Length[vars]]}],{1}],Length[vars]}
]

CountVars[xx__]:=Module[{vars},
  vars = Union[Cases[xx,y[aaa_],Infinity]];
  {xx,Max@@((##[[1]])&/@vars),(*Length[vars]*)10}
]


AdvancedFirstVars[{yyy_,xx__}]:=Module[{yy=x/@yyy,vars,rules},
   vars = Union[Cases[xx,x[aaa_],Infinity]];
  vars=Sort[vars,(If[TrueQ[MemberQ[yy,#1]],If[TrueQ[MemberQ[yy,#2]],#1[[1]]<#2[[1]],True],False])&];
  rules=Apply[Rule,Transpose[{vars,Array[y,Length[vars]]}],{1}];
  {Sort[(##[[1]])&/@(yy/.rules)],xx/.rules}
]









UnsortedUnion[xx_] := Reap[Sow[1, xx], _, #1 &][[2]]

Format[if[x1_, x2_, x3_], InputForm] :=
  "if(" <> ToString[x1, InputForm] <> ")>(f[1])(" <>
   ToString[x2, InputForm] <> ")(" <> ToString[x3, InputForm] <> ")";


MyString[xx_]:=StringReplace[ToString[xx/.{Power->p,Log->l,y->x},InputForm],{"\\"->"","\""->""}];

DoIntegrate[xx__] := Module[{temp2,vars, temp,i,rules,met,res},saved=xx;
vars =
Union[Cases[xx,y[aaa_],Infinity],Cases[xx,x[aaa_],{0,Infinity}]];
If[Length[vars]===0,Return[Plus@@xx]];
        temp= Join[{Plus@@(xx(*If[UsingComplexNumbers,xx /. Log[aaa__]-> Log[aaa + ComplexShift],xx]*))}, {##, IntegrationCut,1} & /@ vars,
        {Method -> AdaptiveQuasiMonteCarlo, Compiled->True,MinRecursion -> 100, MaxRecursion -> 10000, MaxPoints -> 100000, PrecisionGoal -> 4}];
        Return[Apply[NIntegrate, temp, {0}]]
]





EpNegExpand[xx_]:=Module[{i,j,deg,ders,der},
    For[i=1,i<=Length[xx[[1]]],i++,
        deg=(xx[[2]][[i]]/.ep->0);
        If[And[xx[[1]][[i]]===1,deg<0],
            Return[Flatten[EpNegExpand/@Reap[
                If[xx[[2]][[i]]-deg==0,RawPrintLn["Impossible degree"];Abort[]];
                ders={};
                For[j=0,j<Abs[deg],j++,
                    If[j===0,der=xx[[3]],der=D[der,x[i]]/j];
                    AppendTo[ders,der/.x[i]->0];
                    Sow[{ReplacePart[xx[[1]],0,i],ReplacePart[xx[[2]],0,i],ders[[j+1]]/(xx[[2]][[i]]+j+1)}];
                ];
                Sow[{ReplacePart[xx[[1]],2,i],xx[[2]],xx[[3]]-(ders.Table[x[i]^(j-1),{j,Length[ders]}])}];
            ][[2]][[1]],1]];
        ];



    ];
    {xx}

];




KillE[xx_] := Module[{deg, val},
  If[MemberQ[Variables[xx], e],
   deg = xx /. e -> 0;
   val = Cancel[(xx - deg)/e];
   If[deg < -3, 0, val*(10^deg)]
   ,
   xx
   ]
  ]

PMSimplify[xx_]:=PMSimplify[xx,False]
PMSimplify[xx_,MakeBrackets_] := Module[{temp},
  temp = Sort[Variables[Expand[xx]]];
  If[MemberQ[temp,Indeterminate],Return[Indeterminate]];
  (xx /. (Rule[##, 0] & /@
       temp)) + ((Plus @@ (Abs[Coefficient[xx, ##]]^2 & /@ temp))^(1/2))*
    If[MakeBrackets,pm[PMCounter++],ToExpression["pm" <> ToString[PMCounter++]]]
  ]



VarDiffReplace[xx_, yy_] := Module[{temp},
  temp = xx /. VarDiffRule[yy];
  Expand[temp /. y -> x]
]

VarDiffRule[xx_] := Module[{temp, i, j},
  temp =
   Table[If[i === j, If[i < Length[xx], 1, 2],
     If[j === Length[xx], -1, 0]], {i, 1, Length[xx]}, {j, 1,
     Length[xx]}];
  temp = Inverse[temp];
  Apply[Rule,
   Transpose[{(x[##] & /@ xx), temp.(y[##] & /@ xx)}], {1}]
]

NegTerms[0] := {}

NegTerms[xx_] :=
 Select[List @@
   Expand[xx], ((## /. {x[ii_] -> 1, y[ii_] -> 1}) < 0) &]

KillNegTerms[ZZ_]:=Module[{temp,vars,subsets,U,F,UF,Z,i},
(*Return[{ZZ}];*)

    {Z,MM,UF}=ZZ;
    temp=NegTerms[UF[[2]]];
    If[temp==={},Return[{ZZ}]];
    RawPrintLn["Negative terms encountered: ",InputForm[temp]];
    vars=#[[1]]&/@Variables[temp];
    subsets=Subsets[vars,{2}];
    For[i=1,i<=Length[subsets],i++,
        If[And[Length[NegTerms[VarDiffReplace[UF[[2]],subsets[[i]]]]]<Length[temp],
               Length[NegTerms[VarDiffReplace[UF[[2]],{subsets[[i]][[2]],subsets[[i]][[1]]}]]]<Length[temp]]
                ,
            RawPrintLn["Decomposing into two sectors corresponding to ",subsets[[i]]];
            Z={#[[1]]/2,#[[2]]}&/@Z;
            Return[{VarDiffReplace[{Z,MM,UF},subsets[[i]]],VarDiffReplace[{Z,MM,UF},{subsets[[i]][[2]],subsets[[i]][[1]]}]}]
        ];
    ];
    Return[{ZZ}]
]

UseEpMonom[ZZ_]:=Module[{temp,vars,subsets,U,F,Z,MM,i,UF},
    {Z,MM,UF}=ZZ;
    MM=Expand[MM];
    If[And[Head[MM]===Power,Head[Expand[MM[[1]]]]===Plus],
        {{#[[1]],Prepend[#[[2]],MM[[2]]]}&/@Z,Prepend[UF,MM[[1]]]}
    ,
        {{#[[1]]*MM,#[[2]]}&/@Z,UF}
    ]
]

MyDot[xx_,yy_]:=(xx.##)&/@yy

BlockMatrix[MatrixList_, pos_, n_] :=
 Normal[SparseArray[
        Join[Rule[{##,##},1]&/@Complement[Range[n],pos],
     Flatten[Map[Rule[##[[1]][[1]], ##[[2]][[1]]] &, Map[UU, Outer[List, pos, pos], {2}] + Map[VV,##,{2}], {2}], 1]], {n,n}
        ]] &/@MatrixList

MatrixPart[MatrixList_, pos_] :=
 Transpose[Transpose[##[[pos]]][[pos]]] & /@ MatrixList

Reasonable[xx_,uuu_] := If[Abs[xx] > 100uuu, Sign[xx]*100uuu, xx]

PolDegree[xxx_] := Module[{temp,xx},
xx=xxx/.ep->E;
  temp = Variables[xx];
  temp = (## -> t) & /@ temp;
  Exponent[xx /. temp, t]
  ]

PolDegree[xx_,yy_] := Module[{temp},
  temp = (## -> t) & /@ yy;
  Exponent[xx /. temp, t]
  ]


VariableSeries[xx_,var_]:=Module[{temp,shift,s,result,i},

Check[
    temp=xx/.(Log[var]->CCCC);
    temp=temp/.(Log[aaa_]->Normal[Series[Log[aaa],{var,0,VarExpansionDegree}]]);
    temp=Expand[temp];
    shift=0;
    While[TrueQ[MemberQ[{ComplexInfinity,Indeterminate},Quiet[(temp) /. var -> 0]]],
        temp=Expand[temp*var];
        shift++;
    ];
    s=shift;
    While[shift>0,
        temp=D[temp,var];
        shift--;
    ];
    temp=temp/(s!);
    result=0;
    For[i=0,i<=VarExpansionDegree,i++,
        result=result+ (var^i)*(temp/.var->0)*(s!/((s+i)!));
        temp=D[temp,var];
    ];
    result/.(CCCC->Log[var])
   (* ((temp/.var->0)+var((D[temp,var]/.var->0)-s(temp/.var->0)))/.(CCCC->Log[var])*)
,Print[xx];Abort[]]
]

MakeTripple[var_,xx_]:=
    If[UsingC,
        If[Head[xx]===if,
            if[xx[[1]],MakeTripple[var,xx[[2]]],MakeTripple[var,xx[[3]]]],
            if[var,xx,VariableSeries[xx,var]]
        ]
        ,
        If[Head[xx]===If,
            If[xx[[1]],MakeTripple[var,xx[[2]]],MakeTripple[var,xx[[3]]]],
            If[var>IfCut,xx,Evaluate[VariableSeries[xx,var]]]
        ]
    ]


MakeIfConditions[xx_]:=Module[{i,result,llist,rules,var},
 (*   Print[xx];*)
    result=xx[[2]];
    For[i=1,i<=Length[Sort[xx[[1]]]],i++,
        var=y[xx[[1]][[i]]];
        result=MakeTripple[var,result];
    ];
(*    Print[{Sort[xx[[1]]],result}];*)
    {Sort[xx[[1]]],result}
]

ClearDatabase[index_]:=Module[{name},
    name=RemoveFIESTAName[DataPath<>ToString[index]<>"/"];
    QClose[name];
    Pause[1];
    QRemoveDatabase[name];
    QOpen[name];
]

RemoveFIESTAName[xx_]:=StringReplace[xx,{"FIESTA`"->""}];


FreeLink[]:=Max@@(Append[Flatten[Position[LinkFreeQ,True]],0])
ReadyLink[]:=Max@@(Append[Flatten[Position[LinkReadyQ/@CIntegrateLinks,True]],0])

ReadResult[]:=Module[{j,i,result},
    i=0;
    While[i==0,
        If[And[CurrentOrder<MaxOrder,PrepareStringsWhileIntegrating,StringPrepareCounter==TermNumber[CurrentOrder+1]+1],
            If[IntegrationDebug,
                RawPrintLn["     Expression of next order ready for integration!"]
                ,
                RawPrint["!"]
            ];
            StringPrepareCounter=TermNumber[CurrentOrder+1]+2;
        ];
        If[And[CurrentOrder<MaxOrder,PrepareStringsWhileIntegrating,StringPrepareCounter<=TermNumber[CurrentOrder+1]],
            For[j=Min[TermNumber[CurrentOrder+1],StringPrepareCounter+4],StringPrepareCounter<=j,StringPrepareCounter++,
                PrepareString[CurrentOrder+1];
            ];
        ,
            Pause[0.001];
        ];
        i=ReadyLink[]
    ];
(*    RawPrintLn["Received ",i];*)

    LinkFreeQ[[i]]=True;
    result=LinkRead[CIntegrateLinks[[i]]];
    If[And[result[[1]]===Indeterminate,AbortOnNaN],
        RawPrintLn["INDETERMINATE"];
        RawPrintLn[CurrentString[i]];
        Abort[];
    ];
(*    result=CutExtraDigits/@result;*)
   (* RawPrintLn[result];*)
    If[And[NoNan,result==={Indeterminate,Indeterminate}],
        RawPrintLn[""];
        RawPrintLn["Integration returned INDETERMINATE as an answer."];
        RawPrintLn["Run WhyIndeterminate[] for possible reasons;"];
        NoNan=False;
    ];
    result=result[[1]]+result[[2]]*ToExpression["pm"<>ToString[PMCounter++]];
    (*RawPrintLn[result];*)
    If[IntegrationDebug,
        RawPrintLn["     Link ",i," time: ",CutExtraDigits[AbsoluteTime[]-LinkStart[[i]]]," seconds; returned answer: ",PMForm[result]]
        ,
        Dots2[];
    ];
    result3=result3+result;
]

FindFreeLink[]:=Module[{i},
    i=FreeLink[];
    If[i==0,
        ReadResult[];
    ];
    i=FreeLink[];
    i
]

SubmitIntegration[xx_]:=Module[{i},
    i=FindFreeLink[];
    LinkFreeQ[[i]]=False;
    LinkStart[[i]]=AbsoluteTime[];
    If[AbortOnNaN,
        CurrentString[i]=xx;
    ];
    LinkWrite[CIntegrateLinks[[i]],CallPacket[0,{RemoveFIESTAName[xx]}]];
]

BeforePrepareStrings[i_,n_]:=Module[{temp},
        StringPrepareCounter=1;
        PrepareDots[TermNumber[i]];
        ClearDatabase[TargetDatabase];
        NumberOfVariables[TargetDatabase]=0;
        (PSector[TargetDatabase,##]=0)&/@Range[n];
]

PrepareString[i_]:=Module[{temp,num,DName,NOV},
        BlockSize=0;
        ttable=Reap[
        For[SKcounter=1,Or[SKcounter<=NumberOfTerms,BlockSize<CacheSize],SKcounter++,
            If[i===CurrentOrder,Dots[],gi++];CC[100];
            If[StringPrepareCounter<=TermNumber[i],
                BlockSize=BlockSize+StringLength[Sow[QGet[DataPath<>"2/",ToString[i]<>"-"<>ToString[StringPrepareCounter]]]];
                StringPrepareCounter++;
                ,
                Break[]
            ];
        ];
        ][[2]][[1]];
        StringPrepareCounter--;
(*        RawPrint["-",Length[ttable],":",BlockSize,"-"];*)
     (*   NOV=NumberOfVariables[TargetDatabase];
        SetSharedVariable[NOV];*)
         If[NumberOfSubkernels>0,
        DistributeDefinitions[IfCut,VarExpansionDegree,UsingC,if];
        ParallelEvaluate[        Format[if[x1_, x2_, x3_], InputForm] :=
  "if(" <> ToString[x1, InputForm] <> ")>(f[1])(" <>
   ToString[x2, InputForm] <> ")(" <> ToString[x3, InputForm] <> ")"];
   ];
        ttable=MyParallelize[(
            temp=ToExpression[##];
            num=temp[[1]];
            temp=temp[[2]];
            temp[[2]]=Expand[temp[[2]]];
            If[temp[[2]]===0,
                0
                ,
                If[IfCut>0,
                    temp=MakeIfConditions[temp];
                ];
                temp=temp[[2]];
                vars=CountVars[temp];
(*                NOV=Max[NOV,vars[[2]]];             *)
                {num,MyString[temp]<>";\n",vars[[2]]}
            ]
        )&/@ttable];
        ttable=DeleteCases[ttable,0];
   (*     NumberOfVariables[TargetDatabase]=NOV;*)
        DName=RemoveFIESTAName[DataPath<>ToString[TargetDatabase]<>"/"];
        (NumberOfVariables[TargetDatabase]=Max[NumberOfVariables[TargetDatabase],##[[3]]];
            QPut[DName,ToString[##[[1]]]<>"-"<>ToString[++PSector[TargetDatabase,##[[1]]]],##[[2]]];
        )&/@ttable;
]

MyParallelize[x_]:=If[NumberOfSubkernels>0,Parallelize[x(*,Method->"FinestGrained"*)],x];
SetAttributes[MyParallelize,HoldFirst];

SectorDecomposition[zs_,U_,intvars_,n_,j_]:=Module[{res,pol,active,SDDone},
  rules=Rule[x[##],1]&/@(zs);
            If[STRATEGY===STRATEGY_0,
                res={IdentityMatrix[n]};
                Goto[SDDone];
            ];
            If[STRATEGY===STRATEGY_SS,
                pol=Expand[Times@@U]/.rules;
                active=x[##]&/@Range[Length[intvars]];
                res=FindSD[{MyDegrees[pol,active]},MyDeleteEdges[CurrentGraph,zs]];
                active=##[[1]]&/@active;
                res=BlockMatrix[res,active,n];
                Goto[SDDone];
            ];
            If[STRATEGY===STRATEGY_X,
                pol=Expand[##/.rules]&/@U;
                active=(Union@@(Sort[Variables[##]]&/@pol));
                (*Print[MyDegrees[##,active]&/@pol];*)
                res=FindSD[MyDegrees[##,active]&/@pol];
                active=##[[1]]&/@active;
                res=BlockMatrix[res,active,n];
                Goto[SDDone];
            ];
            (* other strategies *)
                pol=Expand[Times@@U]/.rules;
                active=##[[1]]&/@Sort[Variables[pol]];
                res=FindSD[{MyDegrees[pol]}];
                res=BlockMatrix[res,active,n];
            Label[SDDone];
            RawPrintLn["Primary sector ",j," resulted in ",Length[res]," sectors."];
            res
]

SDIntegrate[intvars_,ZZ_,order_Integer,deltas_]:=
Module[{Z,U,F,SD,f,forsd,ii,i,vars,n,m,j,l,rule,Z3,Z2,U2,F2,coeff,k,Jac,res,result,md,result2,pol,active,zsets,SDCoeffs,timecounter},
    NoNan=True;
    timecounter=AbsoluteTime[];
    If[And@@((Not[Head[SDIntegral[intvars,ZZ,deltas,order+1-##]]===SDIntegral])&/@Range[40]),
        Return[];
    ];
fff={};
If[UsingQLink,
    ClearDatabase[1];
    run1=1;
];
    SCounter=0;
For[l=1,l<=Length[ZZ],l++,
 {Z,U}=ZZ[[l]];
 If[Times@@U===0,Continue[]];
 If[Or[Head[SDResolve[{intvars,Z,U}]]===SDResolve,UsingQLink],
    n=Length[intvars];
    vars=Apply[x, Position[intvars, 1], {1}];
    zsets=Tuples[deltas];
    If[Head[PrimarySectorCoefficients]===Symbol,SDCoeffs=Table[1,{Length[zsets]}],SDCoeffs=PrimarySectorCoefficients];
    zsets=Delete[zsets,Position[SDCoeffs,0]];
    SD=Transpose[{zsets,Table[U,{Length[zsets]}],Table[intvars,{Length[zsets]}],Table[n,{Length[zsets]}],Range[Length[zsets]]}];

    RawPrintLn["Sector decomposition - ",Length[SD]," sectors"];
    RawPrintLn["Totally: ",MyTimingForm[AbsoluteTiming[
        If[NumberOfSubkernels>0,
            DistributeDefinitions[STRATEGY];
            SetSharedFunction[RawPrintLn];
        ];
       SD=MyParallelize[(SectorDecomposition@@##)&/@SD];
        If[NumberOfSubkernels>0,
            UnsetShared[RawPrintLn];
        ];

    ][[1]]]," seconds; ",Total[Length/@SD]," sectors."];

    If[ValueQ[SDFile],Put[Sort[SD],SDFile];Abort[]];

    RawPrint["Variable substitution"];
    PrepareDots[Length[SD]];
    RawPrintLn[MyTimingForm[AbsoluteTiming[
    ff=Reap[
        For[j=1,j<=Length[SD],j++,
            If[ForceMixingOfSectors,SCounter++];
            Dots[];
            If[SDCoeffs[[j]]===0,Continue[]];
            rules=Rule[x[##],1]&/@zsets[[j]];

            (*If[ForceMixingOfSectors,
                SD[[j]]={j,##}&/@SD[[j]]
                ,
                SD[[j]]=Transpose[{Range[SCounter+1,SCounter+Length[SD[[j]]]],SD[[j]]}];
                SCounter=SCounter+Length[SD[[j]]];
            ];*)
            If[NumberOfSubkernels>0,
                DistributeDefinitions[n,x,y,rules,rule,j,zsets,SDCoeffs,Z,U,UsingComplexNumbers,ComplexShift];
            ];
             rrr=
            Reap[(
                If[Not[ForceMixingOfSectors],SCounter++];
                reps=Inner[Power[#2, #1] &, Transpose[##], Array[y,n],Times];
                rule=Apply[Rule,Transpose[{Array[x,n],reps}],{1}];
                Jac=Det[Outer[D,reps,Array[y,n]]];
                Z2=(Z/.rules)//.rule;
                U2=(U/.rules)//.rule;
                For[ii=1,ii<=Length[Z2],ii++,

                    NewDegrees=Table[0,{n}];
                    Z3=Z2[[ii]];
                    For[m=1,m<=Length[U2],m++,
                        temp=FactorMonom[U2[[m]] /. {0.->0}];
                        U3[m]=temp[[2]];
                        NewDegrees=NewDegrees+Z3[[2]][[m]]*Table[Exponent[temp[[1]],y[i]],{i,n}];
                    ];
                    Z3[[1]]=Expand[Z3[[1]]];
                    If[Head[Z3[[1]]]===Plus,Z3[[1]]=List@@Z3[[1]],Z3[[1]]={Z3[[1]]}];
                    NewDegrees=NewDegrees+Table[Exponent[Jac,y[i]],{i,n}];
                    f={ReplacePart[intvars,0,List/@zsets[[j]]],Expand[NewDegrees],Z3[[1]],SDCoeffs[[j]]*(Jac/.y[aaa_]->1)*Times@@Table[((If[UsingComplexNumbers,U3[m]+ComplexShift,U3[m]])^(Z3[[2]][[m]])),{m,Length[U2]}]};
                    f=f//.y->x;                     
                    Sow[{SCounter,f}]
                ];
            )&/@SD[[j]]][[2]][[1]];
            Sow[##]&/@rrr;

        ];
    ][[2]][[1]];
    ][[1]]]," seconds. "];
    
   
    
    If[UsingQLink,
        RawPrint["Preparing database: "];
        RawPrintLn[AbsoluteTiming[
            ClearDatabase[3];
            run3=1;
            If[##[[2]][[3]]=!=0,QPut[DataPath<>"3/",ToString[run3,InputForm],ToString[##,InputForm]];run3++]&/@ff;
            Clear[ff];
        ][[1]]," seconds. "]
    ];
    If[UsingQLink,
        DLength[3]=run3-1;
        ClearDatabase[2];
        RawPrint["Decomposing ep-independent term"];
        PrepareDots[DLength[3]];
        RawPrintLn[MyTimingForm[AbsoluteTiming[
            For[run3=1;run2=1,run3<=DLength[3],Null,
                BlockSize=0;
                ttable=Reap[
                For[SKcounter=1,Or[SKcounter<=NumberOfTerms,BlockSize<CacheSize],SKcounter++,
                     Dots[];CC[100];
                     If[run3<=DLength[3],
                         BlockSize=BlockSize+StringLength[Sow[QGet[DataPath<>"3/",ToString[run3]]]];
                         run3++;
                         ,
                         Break[]
                     ];
                ];
                ][[2]][[1]];
(*                RawPrint["-",Length[ttable],":",BlockSize,"-"];*)
                ttable=MyParallelize[(
                    temp=ToExpression[##];
                    num=temp[[1]];
                    temp=ExpandZ[temp[[2]]];
                    ToString[{num,##},InputForm]&/@temp
                )&/@ttable];
                ttable=Flatten[ttable,1];
                (QPut[DataPath<>"2/",ToString[run2,InputForm],##];run2++)&/@ttable;
            ];
            DLength[2]=run2-1;
        ][[1]]]," seconds; ",run2-1," terms."];
        RawPrint["Pole resolution"];
        PrepareDots[DLength[2]];
        RawPrintLn[MyTimingForm[AbsoluteTiming[
            For[run2=1,run2<=DLength[2],Null,
                BlockSize=0;
                ttable=Reap[
                For[SKcounter=1,Or[SKcounter<=NumberOfTerms,BlockSize<CacheSize],SKcounter++,
                    Dots[];CC[100];
                    If[run2<=DLength[2],
                        BlockSize=BlockSize+StringLength[Sow[QGet[DataPath<>"2/",ToString[run2]]]];
                        run2++;
                         ,
                         Break[]
                    ];
                ];
                ][[2]][[1]];
(*                RawPrint["-",Length[ttable],":",BlockSize,"-"];*)
                ttable=MyParallelize[(
                    temp=ToExpression[##];
                    num=temp[[1]];
                    temp=temp[[2]];
                    temp={temp[[1]],temp[[2]],temp[[3]]*temp[[4]]};
                    temp=EpNegExpand[temp];
                    ToString[{num,##},InputForm]&/@temp
                )&/@ttable];
                ttable=Flatten[ttable,1];
                (QPut[DataPath<>"1/",ToString[run1,InputForm],##];run1++)&/@ttable;
          ];
        ][[1]]]," seconds; ",run1-1," terms."];
    ,
        ff=DeleteCases[ff,{aaa_,{a1_,a2_,{0},a4_}}];
        RawPrint["Decomposing ep-independent term"];
        If[NumberOfSubkernels>0,
            SetSharedVariable[gi,gl]
        ];
        PrepareDots[Length[ff]];
        RawPrintLn[MyTimingForm[AbsoluteTiming[
            For[gi=Length[gl],gi>0,gi--,gl[[gi]]={If[gi===1,0,gl[[gi-1]]]+1,gl[[gi]]}];
            ff=Take[ff,##]&/@gl;
            ff=Flatten[   Reap[(
                        Sow[Flatten[MyParallelize[(gi++;CC[100];
                        num=##[[1]];
                        temp=ExpandZ[##[[2]]];
                        {num,##}&/@temp)&/@##],1]];RawPrint["."];
                        )&/@ff
                      ][[2]][[1]],1]
        ][[1]]]," seconds"];
(*        ff=GroupTerms[ff];*)
        RawPrint["Pole resolution"];
        PrepareDots[Length[ff]];
        RawPrintLn[MyTimingForm[AbsoluteTiming[
            For[gi=Length[gl],gi>0,gi--,gl[[gi]]={If[gi===1,0,gl[[gi-1]]]+1,gl[[gi]]}];
            ff=Take[ff,##]&/@gl;
            ff=Flatten[   Reap[(
                        Sow[Flatten[MyParallelize[(gi++;CC[100];
                        num=##[[1]];
                        temp=##[[2]];
                        temp={temp[[1]],temp[[2]],temp[[3]]*temp[[4]]};
                        temp=EpNegExpand[temp];
                        {num,##}&/@temp)&/@##],1]];RawPrint["."];
                        )&/@ff
                      ][[2]][[1]],1]
        ][[1]]]," seconds; ",Length[ff]," terms."];
        SDResolve[{intvars,Z,U}]=ff;
    ]
  ,
  RawPrintLn["Using stored expression with pole resolution."];
  ff=SDResolve[{intvars,Z,U}];
  ];
  If[Not[UsingQLink],
    AppendTo[fff,ff];
  ]
];

    If[UsingQLink,
        DLength[1]=run1-1;
        RawPrint["Expression construction"];
        PrepareDots[DLength[1]];
        RawPrintLn[MyTimingForm[AbsoluteTiming[
            For[run1=1,run1<=DLength[1],Null,
                BlockSize=0;
                ttable=Reap[
                For[SKcounter=1,Or[SKcounter<=NumberOfTerms,BlockSize<CacheSize],SKcounter++,
                    Dots[];CC[100];
                    If[run1<=DLength[1],
                        temp=QGet[DataPath<>"1/",ToString[run1]];
                        Sow[{run1,temp}];
                        BlockSize=BlockSize+StringLength[temp];
                        run1++;
                         ,
                         Break[]
                    ];
                ];
                ][[2]][[1]];
(*                RawPrint["-",Length[ttable],":",BlockSize,"-"];*)
                ttable=MyParallelize[(
                    num=##[[1]];
                    temp=ToExpression[##[[2]]];
                    {ToString[num],ToString[{temp[[1]],ConstructTerm[temp[[2]]]},InputForm]}
                )&/@ttable];
                QPut[DataPath<>"1/",##[[1]],##[[2]]]&/@ttable;
            ];
        ][[1]]]," seconds."];

        RawPrint["Replacing variables"];
        PrepareDots[DLength[1]];
        RawPrintLn[MyTimingForm[AbsoluteTiming[
            For[run1=1,run1<=DLength[1],Null,
                BlockSize=0;
                ttable=Reap[
                For[SKcounter=1,Or[SKcounter<=NumberOfTerms,BlockSize<CacheSize],SKcounter++,
                    Dots[];CC[100];
                    If[run1<=DLength[1],
                        temp=QGet[DataPath<>"1/",ToString[run1]];
                        Sow[{run1,temp}];
                        BlockSize=BlockSize+StringLength[temp];
                        run1++;
                         ,
                         Break[]
                    ];
                ];
                ][[2]][[1]];
(*                RawPrint["-",Length[ttable],":",BlockSize,"-"];*)
                ttable=MyParallelize[(
                    num=##[[1]];
                    temp=ToExpression[##[[2]]];
                    {ToString[num],ToString[{temp[[1]],AdvancedFirstVars[temp[[2]]]},InputForm]}
                )&/@ttable];
                QPut[DataPath<>"1/",##[[1]],##[[2]]]&/@ttable;
            ];
        ][[1]]]," seconds."];
        RawPrint["Epsilon expansion"];
        ClearDatabase[2];
        PrepareDots[DLength[1]];
        min=order+1;
        (TermNumber[##]=0)&/@Range[-20,20];
        RawPrintLn[MyTimingForm[AbsoluteTiming[

            For[run1=1,run1<=DLength[1],Null,
                BlockSize=0;
                ttable=Reap[
                For[SKcounter=1,Or[SKcounter<=NumberOfTerms,BlockSize<CacheSize],SKcounter++,
                    Dots[];CC[100];
                    If[run1<=DLength[1],
                        BlockSize=BlockSize+StringLength[Sow[QGet[DataPath<>"1/",ToString[run1]]]];
                        run1++;
                         ,
                         Break[]
                    ];
                ];
                ][[2]][[1]];
(*                RawPrint["-",Length[ttable],":",BlockSize,"-"];        *)
                ttable=MyParallelize[(
                    temp=ToExpression[##];
                    num=temp[[1]];
                    temp=temp[[2]];
                    temp1=temp[[1]];
                    temp2=EpPosExpand[temp[[2]],order];
                    temp={num,temp1,##}&/@temp2;
                    temp={##[[3]][[1]],ToString[{##[[1]],{##[[2]],##[[3]][[2]]}},InputForm]}&/@temp
                )&/@ttable];
                ttable=Flatten[ttable,1];
                (TermNumber[##[[1]]]=TermNumber[##[[1]]]+1;
                    min=Min[min,##[[1]]];
                    QPut[DataPath<>"2/",ToString[##[[1]]]<>"-"<>ToString[TermNumber[##[[1]]]],##[[2]]]
                )&/@ttable;
            ];
        ][[1]]]," seconds."];
    ,(*not UsingQLink*)
        ff=Flatten[fff,1];
        RawPrint["Expression construction"];
        PrepareDots[Length[ff]];
        For[gi=Length[gl],gi>0,gi--,gl[[gi]]={If[gi===1,0,gl[[gi-1]]]+1,gl[[gi]]}];
        ff=Take[ff,##]&/@gl;
        RawPrintLn[MyTimingForm[AbsoluteTiming[
            ff=Flatten[   Reap[(
                Sow[MyParallelize[(gi++;CC[100];{##[[1]],ConstructTerm[##[[2]]]})&/@##]];RawPrint["."];
            )&/@ff][[2]][[1]],1]
        ][[1]]]," seconds."];
        RawPrint["Replacing variables"];
        PrepareDots[Length[ff]];
        For[gi=Length[gl],gi>0,gi--,gl[[gi]]={If[gi===1,0,gl[[gi-1]]]+1,gl[[gi]]}];
        ff=Take[ff,##]&/@gl;
        RawPrintLn[MyTimingForm[AbsoluteTiming[
            ff=Flatten[   Reap[(
                Sow[MyParallelize[(gi++;CC[100];{##[[1]],AdvancedFirstVars[##[[2]]]})&/@##]];RawPrint["."];
            )&/@ff][[2]][[1]],1]
        ][[1]]]," seconds."];
        RawPrint["Epsilon expansion"];
        PrepareDots[Length[ff]];
        RawPrintLn[MyTimingForm[AbsoluteTiming[
            For[gi=Length[gl],gi>0,gi--,gl[[gi]]={If[gi===1,0,gl[[gi-1]]]+1,gl[[gi]]}];
            ff=Take[ff,##]&/@gl;
            fff=Flatten[   Reap[(
                        Sow[Flatten[MyParallelize[(gi++;CC[100];
                        num=##[[1]];
                        temp=##[[2]];
                        temp1=temp[[1]];
                        temp2=EpPosExpand[Expand[temp[[2]]],order];
                        {num,{##[[1]],{temp1,##[[2]]}}}&/@temp2)&/@##],1]];RawPrint["."];
                        )&/@ff
                      ][[2]][[1]],1]
        ][[1]]]," seconds."];

        fff=DeleteCases[fff,{aaaaaaa_,{aaa_,0}}];
        If[fff==={},
            min=order+1,
            min=Min@@(##[[2]][[1]]&/@fff);
        ];
    ];



    For[i=min-1,i>=min-50,i--,SDIntegral[intvars,ZZ,deltas,i]=0];
    If[min===order+1,Return[]];


    If[UsingQLink,
        TargetDatabase=1;
        BeforePrepareStrings[min,SCounter];
        MaxOrder=order;
    ];

    For[CurrentOrder=min,CurrentOrder<=order,CurrentOrder++,
        If[Not[Head[SDIntegral[intvars,ZZ,deltas,CurrentOrder]]===SDIntegral],Continue[]];
        u=1;
        If[Not[UsingQLink],
            ff={##[[1]],##[[2]][[2]]}&/@Cases[fff,{aaaaa_,{CurrentOrder,aaa__}}];
        ];
        If[UsingC,
            If[UsingQLink,
             If[Not[StringPrepareCounter>TermNumber[CurrentOrder]],
                If[IfCut>0,
                    RawPrint["Expanding, making if conditions and integration string"];
                ,
                    RawPrint["Expanding, making integration string"]
                ];

                If[StringPrepareCounter>1,RawPrint[" (",Round[100*(StringPrepareCounter-1)/TermNumber[CurrentOrder]],"% done)"];Dots[]];
                RawPrintLn[MyTimingForm[AbsoluteTiming[
                    For[Null,StringPrepareCounter<=TermNumber[CurrentOrder],StringPrepareCounter++,
                        PrepareString[CurrentOrder]
                    ];
                ][[1]]]," seconds."];
             ];

                ReadDatabase=TargetDatabase;
                TargetDatabase=4-TargetDatabase;
                If[CurrentOrder<MaxOrder,
                    BeforePrepareStrings[CurrentOrder+1,SCounter];
                ];



                RawPrintLn["Terms of order ",CurrentOrder,": ",Plus@@(PSector[ReadDatabase,##]&/@Range[SCounter])," (",NumberOfVariables[ReadDatabase],"-fold integrals)."];
                RawPrint["Numerical integration: "];
                RawPrintLn[SCounter-Length[Position[PSector[ReadDatabase,##]&/@Range[SCounter],0]]," parts; ",NumberOfLinks," links;"];
                PrepareDots2[SCounter];
                If[Not[IntegrationDebug],RawPrint["Integrating"]];
                RawPrintLn[If[IntegrationDebug,"Total: ",""],AbsoluteTiming[
                    result3=0;
                    For[j=1,j<=SCounter,j++,
                        If[PSector[ReadDatabase,j]==0,If[Not[IntegrationDebug],Dots2[]];Continue[]];
                        If[NumberOfVariables[ReadDatabase]===0,
                            result3=result3+Plus@@((ToExpression[StringDrop[QGet[RemoveFIESTAName[DataPath<>ToString[ReadDatabase]<>"/"],ToString[j]<>"-"<>ToString[##]],-2]])&/@Range[PSector[ReadDatabase,j]]);
                        ,
                            iii=FindFreeLink[];
                            LinkFreeQ[[iii]]=False;
                            LinkStart[[iii]]=AbsoluteTime[];
                            LinkWrite[CIntegrateLinks[[iii]],CallPacket[3,{RemoveFIESTAName[ToString[Round[QSize[DataPath<>ToString[ReadDatabase]<>"/"]]]<>";\n"]}]];
                            LinkRead[CIntegrateLinks[[iii]]];
                            LinkWrite[CIntegrateLinks[[iii]],CallPacket[3,{RemoveFIESTAName[ToString[NumberOfVariables[ReadDatabase]]<>";\n"<>ToString[PSector[ReadDatabase,j]+1]<>";\n"<>StringReplace[ToString[IfCut,InputForm],"*^" -> "e"]<>";\n"]}]];
                            LinkRead[CIntegrateLinks[[iii]]];
(*                          (Print[{RemoveFIESTAName[QGet[RemoveFIESTAName[DataPath<>ToString[ReadDatabase]<>"/"],ToString[j]<>"-"<>ToString[##]]]}];)&/@Range[PSector[ReadDatabase,j]];                            *)
                            (LinkWrite[CIntegrateLinks[[iii]],CallPacket[3,{RemoveFIESTAName[QGet[RemoveFIESTAName[DataPath<>ToString[ReadDatabase]<>"/"],ToString[j]<>"-"<>ToString[##]]]}]];LinkRead[CIntegrateLinks[[iii]]];)&/@Range[PSector[ReadDatabase,j]];
                            LinkWrite[CIntegrateLinks[[iii]],CallPacket[0,{RemoveFIESTAName[(StringJoin@@(("f["<>ToString[##]<>"]+")&/@Range[2,PSector[ReadDatabase,j]+1]))<>"0;\n"]}]];
                        ];
                    ];
                While[Not[And@@LinkFreeQ],ReadResult[]];
                ][[1]]," seconds; returned answer: ",PMForm[DoubleCutExtraDigits[PMSimplify[result3]]]];
            ,   (*not UsingQLink*)
                RawPrint["Expanding"];
                PrepareDots[Length[ff]];
                For[gi=Length[gl],gi>0,gi--,gl[[gi]]={If[gi===1,0,gl[[gi-1]]]+1,gl[[gi]]}];
                ff=Take[ff,##]&/@gl;
                RawPrintLn[MyTimingForm[AbsoluteTiming[
                    ff=Flatten[   Reap[(
                        Sow[MyParallelize[(gi++;CC[100];{##[[1]],Expand[##[[2]]]})&/@##]];RawPrint["."];
                    )&/@ff][[2]][[1]],1]
                ][[1]]]," seconds."];
                If[IfCut>0,
                    RawPrint["Making if conditions"];
                    PrepareDots[Length[ff]];
                    For[gi=Length[gl],gi>0,gi--,gl[[gi]]={If[gi===1,0,gl[[gi-1]]]+1,gl[[gi]]}];
                    ff=Take[ff,##]&/@gl;
                    RawPrintLn[MyTimingForm[AbsoluteTiming[
                        ff=Flatten[   Reap[(
                            Sow[MyParallelize[(gi++;CC[100];{##[[1]],MakeIfConditions[##[[2]]]})&/@##]];RawPrint["."];
                        )&/@ff][[2]][[1]],1]
                    ][[1]]]," seconds."];
                ];
                ff={##[[1]],##[[2]][[2]]}&/@ff;
                RawPrint["Counting variables: "];
                PrepareDots[Length[ff]];
                RawPrintLn[MyTimingForm[AbsoluteTiming[
                    NumberOfVariables=Max@@MyParallelize[((CountVars[##[[2]]][[2]])&/@ff)];
                ][[1]]]," seconds."];
                If[NumberOfVariables===-Infinity,
                    RawPrintLn["Terms of order ",CurrentOrder,": ",Length[ff]," (0-fold integrals)."];
                    RawPrint["Numerical integration: "];
                    result3=Plus@@((##[[2]])&/@ff);
                    RawPrintLn[result3];
                ,

                    RawPrint["Preparing integration string"];
                    PrepareDots[Length[ff]];
                    RawPrintLn[MyTimingForm[AbsoluteTiming[
                        ff=(Dots[];CC[100];{##[[1]],MyString[##[[2]]]<>";\n"})&/@ff;
                    ][[1]]]," seconds."];
                    result3=0;
                    RawPrintLn["Terms of order ",CurrentOrder,": ",Length[ff]," (",NumberOfVariables,"-fold integrals)."];
                    RawPrint["Numerical integration: "];
                    RawPrintLn[Length[Union[(##[[1]])&/@ff]]," parts; ",NumberOfLinks," links;"];
                    If[Not[IntegrationDebug],RawPrint["Integrating"]];
                    PrepareDots2[SCounter];
                    RawPrintLn[If[IntegrationDebug,"Total: ",""],AbsoluteTiming[
                        result3=0;
                        For[j=1,j<=SCounter,j++,
                            uuu=Cases[ff,{j,aaa_}];
                            If[Length[uuu]==0,If[Not[IntegrationDebug],Dots2[]];Continue[]];
                            uuu=(##[[2]])&/@uuu;
                            OutputString=ToString[NumberOfVariables]<>";\n"<>ToString[Length[uuu]+1]<>";\n"<>StringReplace[ToString[IfCut,InputForm],"*^" -> "e"]<>";\n"<>(StringJoin@@uuu);
                            OutputString=OutputString<>(StringJoin@@(("f["<>ToString[##]<>"]+")&/@Range[2,Length[uuu]+1]))<>"0;\n";
                            OutputString=ToString[StringLength[OutputString]]<>";\n"<>OutputString;
                            result=SubmitIntegration[OutputString];
                        ];
                        While[Not[And@@LinkFreeQ],ReadResult[]];
                    ][[1]]," seconds; returned answer: ",PMForm[DoubleCutExtraDigits[PMSimplify[result3]]]];
                ];
            ];
        , (* not UsingC  *)
            ff=(##[[2]])&/@ff;


                RawPrint["Expanding"];
                PrepareDots[Length[ff]];
                For[gi=Length[gl],gi>0,gi--,gl[[gi]]={If[gi===1,0,gl[[gi-1]]]+1,gl[[gi]]}];
                ff=Take[ff,##]&/@gl;
                RawPrintLn[MyTimingForm[AbsoluteTiming[
                    ff=Flatten[   Reap[(
                        Sow[MyParallelize[(gi++;CC[100];Expand[##])&/@##]];RawPrint["."];
                    )&/@ff][[2]][[1]],1]
                ][[1]]]," seconds."];
                If[IfCut>0,
                    RawPrint["Making if conditions"];
                    PrepareDots[Length[ff]];
                    For[gi=Length[gl],gi>0,gi--,gl[[gi]]={If[gi===1,0,gl[[gi-1]]]+1,gl[[gi]]}];
                    ff=Take[ff,##]&/@gl;
                    RawPrintLn[MyTimingForm[AbsoluteTiming[
                        ff=Flatten[   Reap[(
                            Sow[MyParallelize[(gi++;CC[100];MakeIfConditions[##])&/@##]];RawPrint["."];
                        )&/@ff][[2]][[1]],1]
                    ][[1]]]," seconds."];
                ];
            ff=(##[[2]])&/@ff;

            result3=0;
            RawPrintLn["Terms of order ",CurrentOrder,": ",Length[ff],"."];


                RawPrint["Numerical integration: "];
                NC=If[NumberOfSubkernels>0,NumberOfSubkernels*10,1];
                ff=Sort[ff,(Random[Integer]===1)&];
                gn=Length[ff];gl={};For[gi=1,gi<=NC,gi++,AppendTo[gl,Quotient[gi*gn,NC]]];gi=0;
                For[gi=Length[gl],gi>0,gi--,gl[[gi]]={If[gi===1,0,gl[[gi-1]]]+1,gl[[gi]]}];
                ff=Take[ff,##]&/@gl;
                DistributeDefinitions[IntegrationCut,ComplexShift,UsingComplexNumbers];
                RawPrintLn[AbsoluteTiming[
                result3=Plus@@MyParallelize[
                    (
                    result2=Reap[Block[{$MessagePrePrint = Sow, $Messages = {}},
                        Quiet[
                            result=DoIntegrate[##]
                        ,{Power::infy, Power::indet, Infinity::indet,General::stop}
                    ];
                    ]][[2]];
                    result2=DeleteCases[##, HoldForm[$MessageList]] & /@result2;
                    Good=False;
                    If[And[Length[$MessageList]>=1,Last[$MessageList]===HoldForm[NIntegrate::"maxp"],Length[result2]===1,
                        ReleaseHold[result2[[1]][[2]]]===result],
                        (*result=ReleaseHold[result2[[1]][[2]]]+ReleaseHold[result2[[1]][[3]]]*ToExpression["pm"<>ToString[PMCounter++]];*)
                        Good=True;
                        result3=ReleaseHold/@Drop[result2[[1]],1];
                        (*result3=result3[[1]]+result3[[2]]*ToExpression["pm"<>ToString[PMCounter++]];*)
                    ];
                    If[Length[result2]===0,
                        Good=True;
                        result3={result,0};
                    ];
                    If[Not[Good],
                        RawPrintLn[PMForm[result]];
                        RawPrintLn["Something went wrong"];
                        RawPrintLn[$MessageList];
                        RawPrintLn[result2];
                        Abort[];
                    ];
                    result3
                    )&/@ff];
                 result3=result3[[1]]+result3[[2]]*ToExpression["pm"<>ToString[PMCounter++]];
                ][[1]]," seconds."];


        ];
        result3=PMSimplify[result3];
        result3=DoubleCutExtraDigits[result3];
        RawPrintLn["Integration of order ",CurrentOrder,": ",PMForm[result3]];
        SDIntegral[intvars,ZZ,deltas,CurrentOrder]=result3;
        RawPrintLn[PMSeriesForm[GenerateAnswer[intvars,ZZ,CurrentOrder,CurrentOrder-(order-ORDER),deltas,OUTSIDE,DENOM]]];
    ];


    RawPrintLn["Total time used: ",CutExtraDigits[AbsoluteTime[]-timecounter]," seconds."];
    Return[];


];


MyDelta[xx_,yy_]:=Table[If[iii===xx,1,0],{iii,yy}]


TestLink[]:=Module[{temp,i,j},
    Print[VersionString];
    If[UsingC==False,UsingQLink=False];
    If[UsingQLink,
        If[Head[QLink]===Symbol,
            QLink=Install[QLinkPath];
            If[Not[Head[QLink]===LinkObject],
                Print["Could not start QLink"];
                Abort[];
            ];
            If[ToString[FileType[DataPath<>"1"]] === "Directory",
                QRemoveDatabase[DataPath<>"1/"];
            ];
            If[Not[QOpen[DataPath<>"1/"]],
                Print["Could not open database"];
                Abort[];
            ];
            If[ToString[FileType[DataPath<>"2"]] === "Directory",
                QRemoveDatabase[DataPath<>"2/"];
            ];
            If[Not[QOpen[DataPath<>"2/"]],
                Print["Could not open database"];
                Abort[];
            ];
            If[ToString[FileType[DataPath<>"3"]] === "Directory",
                QRemoveDatabase[DataPath<>"3/"];
            ];
            If[Not[QOpen[DataPath<>"3/"]],
                Print["Could not open database"];
                Abort[];
            ];
        ];
        Print["Database ready"];
    ];
    If[UsingC,
        If[Head[CIntegrateLinks]===Symbol,CIntegrateLinks=Table[False,{NumberOfLinks}];
                                        LinkFreeQ=Table[False,{NumberOfLinks}];
                                        LinkStart=Table[0,{NumberOfLinks}]
                                        ];
        For[i=1,i<=NumberOfLinks,i++,
            Quiet[While[LinkReadyQ[CIntegrateLinks[[i]]],
                LinkRead[CIntegrateLinks[[i]]]
            ]];
            If[Not[MemberQ[Links[],CIntegrateLinks[[i]]]],
                CIntegrateLinks[[i]] = Install[CIntegratePath];
            ];
            LinkFreeQ[[i]]=True;
            LinkWrite[CIntegrateLinks[[i]],CallPacket[1,VegasSettings]];
            LinkRead[CIntegrateLinks[[i]]];

            LinkWrite[CIntegrateLinks[[i]],CallPacket[5,{CurrentIntegrator}]];
            NewBinary=(LinkRead[CIntegrateLinks[[i]]]=!=$Failed);
            If[NewBinary,
                If[ValueQ[CurrentIntegratorSettings],
                   For[j=1,j<=Length[CurrentIntegratorSettings],j++,
                       LinkWrite[CIntegrateLinks[[i]],CallPacket[6,CurrentIntegratorSettings[[j]]]];
                       If[i == 1,
                          RawPrintLn[CurrentIntegratorSettings[[j]],
                            ": ",LinkRead[CIntegrateLinks[[i]]]];,
                          LinkRead[CIntegrateLinks[[i]]]
                       ];
                   ];
                ];
            ];
            If[IntegrationCut>0,
                LinkWrite[CIntegrateLinks[[i]],CallPacket[2,{0,IntegrationCut}]];
                LinkRead[CIntegrateLinks[[i]]];
            ];

            If[Not[Head[CIntegrateLinks[[i]]]===LinkObject],
                Print["Could not start external intgration link"];
                Abort[];
            ];
        ];
        If[NewBinary,
            RawPrintLn["CurrentIntegrator=",CurrentIntegrator,";"];
            LinkWrite[CIntegrateLinks[[1]],CallPacket[7,{}]];
            RawPrintLn["CurrentIntegratorSettings=\n",LinkRead[CIntegrateLinks[[1]]],";"];
            If[ValueQ[CurrentIntegratorSettings],
               Clear[CurrentIntegratorSettings];
            ];
        ];
    ];

    If[Or[Not[IntegerQ[NumberOfSubkernels]],NumberOfSubkernels<0],
        RawPrintLn["Incorrect number of Subkernels"];
        RawPrintLn["NumberOfSubkernels set to 0"];
        NumberOfSubkernels=0;
    ];
    If[And[NumberOfSubkernels>0,Not[TrueQ[$VersionNumber>=7.0]]],
        RawPrintLn["Subkernels can be used only starting from Mathematica 7.0"];
        RawPrintLn["NumberOfSubkernels set to 0"];
        NumberOfSubkernels=0;
    ];
(*    If[NumberOfSubkernels===1,
        RawPrintLn["The number of Subkernels should be 0 or greater than 1"];
        RawPrintLn["NumberOfSubkernels set to 0"];
        NumberOfSubkernels=0;
    ];*)
    If[NumberOfSubkernels>0,
        RawPrintLn["Starting ",NumberOfSubkernels," Subkernels"];
        PrepareStringsWhileIntegrating=False;
        RawPrintLn["PrepareStringsWhileIntegrating set to False"];
        CloseKernels[];
        LaunchKernels[NumberOfSubkernels];
        (*Quiet[LaunchKernels[NumberOfSubkernels], {Parallel`Queue`Interface`EmptyQ::shdw}];*)
        DistributeDefinitions["FIESTA`"];
    ];
    NumberOfTerms=If[NumberOfSubkernels===0,CacheSize=0;1,NumberOfSubkernels*100];
    CodeInfo[];
]





(*
SDEvaluate[{U_,F1_,F2_,h_},degrees1_,order_]:=Module[{part1,ZZ,n,degrees,epdegrees,runorder,outside,a1,a2},
    GraphUsed=False;
    TestLink[];
    n=Length[degrees1];
    degrees=degrees1/.{ep->0};
    epdegrees=Expand[degrees1-(degrees1/.{ep->0})]/ep;
    part1=If[Or[MemberQ[Variables[F1],x[##]],And[degrees[[##]]<=0,epdegrees[[##]]==0]],1,0]&/@Range[n];
    RawPrintLn[part1];
    SDEvaluate[{U,F1,F2,h},degrees1,part1,order]
]
SDEvaluate[{U_,F1_,F2_,h_},degrees1_,part1_,order_]:=Module[{ZZ,n,degrees,epdegrees,runorder,outside,a1,a2},
    n=Length[degrees1];
    degrees=degrees1/.{ep->0};
    epdegrees=Expand[degrees1-(degrees1/.{ep->0})]/ep;
    a1=Plus@@(degrees1*part1);
    a2=Plus@@(degrees1-degrees1*part1);
    ZZ={{{1,{a1+a2-(h+1)(2-ep),-(a1+a2/2-h(2-ep)),-a2/2}}},{U,F1,F2}};
    runorder=order;
    If[((a1+a2/2-h(2-ep))/.ep->0)<=0,runorder++];
    For[j=1,j<=Length[degrees],j++,
		If[And[degrees[[j]]<=0,Not[epdegrees[[j]]==0]],runorder--]
    ];
    outside=(E^(h EulerGamma ep))*Gamma[a1+a2/2-h(2-ep)]*Gamma[a2/2]/2;
    deltas={Flatten[Position[part1,1]],Flatten[Position[part1,0]]};
    {ZZ,deltas,degrees,epdegrees}=KillNegativeIndices[ZZ,deltas,degrees,epdegrees];
    SDEvaluate1[ZZ,degrees,epdegrees,runorder,order,outside,deltas]
]*)



SDEvaluateG[{graph_,external_},{U_,F_,h_},degrees1_,order_]:=Module[{temp},
    GraphInput=True;
    InfinityVertex=Apply[Max,graph,{0,Infinity}]+1;
    CurrentGraph=MyGraph[Join[graph,{##,InfinityVertex}&/@external]];
    Exte=Length[external];
    External=external;
    SDEvaluate[{U,F,h},degrees1,order]
]


SDEvaluate[{U_,F_,h_},degrees1_,order_]:=Module[{ZZ,n,degrees,epdegrees,runorder,outside,a1,deltas},
    If[GraphInput,GraphUsed=True,GraphUsed=False,GraphUsed=False];
    TestLink[];
    GraphInput=False;
    n=Length[degrees1];
    If[GraphUsed,
        If[Not[M[CurrentGraph]-Exte===Length[degrees1]],
            RawPrintLn["ERROR: length of indices is different from the number of lines in the graph"];
            Abort[];
        ];
    ];
    degrees=degrees1/.{ep->0};
    epdegrees=Expand[degrees1-(degrees1/.{ep->0})]/ep;
    a1=Plus@@degrees1;
    ZZ={{{1,{a1-(h+1)(d0/2-ep),-(a1-h(d0/2-ep))}}},{U,F}};
    runorder=order;
    If[((a1-(h)(d0/2-ep))/.ep->0)<=0,runorder++];
    For[j=1,j<=Length[degrees],j++,
		If[And[degrees[[j]]<=0,Not[epdegrees[[j]]==0]],runorder--]
    ];
    outside=(E^(h EulerGamma ep))*Gamma[a1-h(d0/2-ep)];
    deltas={Range[n]};
    {ZZ,deltas,degrees,epdegrees}=KillNegativeIndices[ZZ,deltas,degrees,epdegrees];

    SDEvaluate1[ZZ,degrees,epdegrees,runorder,order,outside,deltas]
]



KillNegativeIndices[ZZ1_,deltas1_,degrees1_,epdegrees1_]:=Module[{temp,i,j,var,ZZ,degrees,epdegrees,result},
    ZZ=ZZ1;
    deltas=deltas1;
    degrees=degrees1;
    epdegrees=epdegrees1;
	For[i=1,i<=Length[degrees],i++,var=x[i];
        If[epdegrees[[i]]===0,
            While[degrees[[i]]<0,
                ZZ[[1]]=Flatten[
                    Append[
                        Table[
                            {-##[[1]]*D[ZZ[[2]][[j]],var]*##[[2]][[j]],##[[2]]-MyDelta[j,Length[ZZ[[2]]]]},
                            {j,1,Length[ZZ[[2]]]}
                        ],
                        {-D[##[[1]],var],##[[2]]}
                   ]&/@ZZ[[1]],1];

                degrees[[i]]=degrees[[i]]+1;
            ];
            If[degrees[[i]]===0,
                ZZ=Expand[ZZ/.var->0];
                deltas=DeleteCases[##,i]&/@deltas;
                If[GraphUsed,
                    CurrentGraph=ContractEdge[CurrentGraph,i];
                    If[CurrentGraph===False,
                        RawPrintLn["WARNING: loop contracted, possible zero in the answer"];
                        Abort[];
                    ];
                ];
            ];
        ]
	];
    ZZ[[1]]=Expand[ZZ[[1]]];
    ZZ[[1]]=DeleteCases[ZZ[[1]],{0,{aaa1_,aaa2_}}];
    Return[{ZZ,deltas,degrees,epdegrees}];
];


SDEvaluate1[ZZ1_,degrees_,epdegrees_,runorder_,order_,outside_,deltas_]:=Module[{temp,i,j,var,ZZ},
    ZZ=ZZ1;
    If[Length[ZZ[[1]]]===0,Return[0]];
    n=Length[degrees];
    intvars=Table[(If[Or[degrees[[i]]>0,Not[epdegrees[[i]]===0]],1,0]),{i,1,n}];
    ZZ[[1]]={#[[1]]*(Times@@Table[x[i]^(If[degrees[[i]]>0,degrees[[i]]-1,0]),{i,1,n}]),#[[2]]}&/@ZZ[[1]];
    ZZ={{ZZ[[1]],(Times@@Table[x[i]^(If[intvars[[i]]===1,epdegrees[[i]]*ep+If[degrees[[i]]<=0,degrees[[i]]-1,0],0]),{i,1,n}]),ZZ[[2]]}};
    While[True,
        ZZOld=ZZ;
        If[ResolveNegativeTerms,
            ZZ=Flatten[KillNegTerms/@ZZ,1];
        ];
        If[ZZ===ZZOld,Break[]];
    ];
    ZZ=UseEpMonom/@ZZ;
    denom=Times@@Table[If[And[degrees[[i]]<=0,epdegrees[[i]]==0],1,Gamma[degrees[[i]]+ep*epdegrees[[i]]]],{i,Length[degrees]}];
    RawPrintLn["Integration has to be performed up to order ",runorder];
    OUTSIDE=outside;
    DENOM=denom;
    ORDER=order;
    SDIntegrate[intvars,ZZ,runorder,deltas];
    result=GenerateAnswer[intvars,ZZ,runorder,ORDER,deltas,OUTSIDE,DENOM];
(*    RawPrintLn[PMSeriesForm[result]];*)
    Normal[result]
]

GenerateAnswer[intvars_,ZZ_,runorder_,order_,deltas_,outside_,denom_]:=Module[{temp},
    temp=Series[outside*(Plus@@Table[(SDIntegral[intvars,ZZ,deltas,j])*ep^j,{j,runorder-30,runorder}])/denom,{ep,0,order}];
    If[temp===0,Return[0]];
    If[temp===Indeterminate,Return[Indeterminate]];
    temp[[3]]=PMSimplify[##,ReturnErrorWithBrackets]&/@temp[[3]];
    temp[[3]]=DoubleCutExtraDigits/@temp[[3]];
    temp
]


	
	



















MyRank[xx_] := MatrixRank[(## - xx[[1]]) & /@ xx]
VNorm[xx_] := Times @@ ((## + 1) & /@ xx);
PointOver[xx_, yy_] := And @@ ((## >= 0) & /@ (xx - yy));
MinVector[v_] := v/(GCD @@ v)
OnlyLowPoints[xx_] := Module[{temp, i, j, Good},
   temp = {};
   For[i = 1, i <= Length[xx], i++,
    j = 1;
    Good = True;
    While[j <= Length[temp],
     If[PointOver[xx[[i]], temp[[j]]], Good = False; Break[]];
     If[PointOver[temp[[j]], xx[[i]]], temp = Delete[temp, j];
      Continue[]];
     j++
     ];
    If[Good, AppendTo[temp, xx[[i]]]];
    ];
   temp
   ];

MyDegrees[pol_]:=MyDegrees[pol,Sort[Variables[pol]]];

MyDegrees[pol_,var_] := Module[{rule, temp, degrees, i},
(*    var=Array[x,n];*)
  rule = Rule[##, 1] & /@ var;
  temp = List @@ Expand[pol/.(0.->0)];
  temp = Union[(##/(## /. rule)) & /@ temp];
  degrees =
   Table[Exponent[##, var[[i]]], {i, Length[var]}] & /@ temp;
  degrees = OnlyLowPoints[degrees];
  degrees
]

SubSpace[xxx_] :=
  Module[{av, m, n, k, i, j, ind, facet, vector, temp, scalar, pr,
    newpoint, symplex, xx, sbs, km},
    If[Length[xxx[[1]]]===1,Return[{1}]];
   xx = Sort[xxx, (VNorm[#1] < VNorm[#2]) &];
   n = Length[xx[[1]]];
   k = n;
   m = MyRank[OnlyLowPoints[xx]];
   sbs = Subsets[Range[n], {2, Max[m + 1, n]}];
   sbs = Sort[sbs, (Plus @@ #1 > Plus @@ #2) &];
   For[i = 1, i <= Length[sbs], i++,
    av =
     Normal[SparseArray[Apply[Rule, (({##, 1}) & /@ sbs[[i]]), {1}],n]];
    If[MyRank[OnlyLowPoints[(##*av) & /@ xx]] >= Length[sbs[[i]]] - 1,
     Return[av];
     ];
    ];
   RawPrintLn["No subspace found"];
   Return[Table[0, {n}]];
   ];
Facet[xxx_] :=
 Module[{av, aa, bb, m, n, k, i, j, ind, pos, facet, vector, temp,
   scalar, pr, newpoint, symplex, xx, sbs, km, yy},
  xx = xxx;
  n = Length[xx[[1]]];
  If[n===1,Return[{{1},{1}}]];
  av = SubSpace[xx];
  Label[avlabel];
  xx = Sort[xx, (VNorm[#1*av] < VNorm[#2*av]) &];
  yy = {};
  For[i = 1, i <= Length[xx], i++,
   If[And[Length[yy] > 0,
     MyRank[(##*av) & /@ yy] ==
      MyRank[(##*av) & /@ Append[yy, xx[[i]]]]], Continue[]];
   AppendTo[yy, xx[[i]]];
   yy = OnlyLowPoints[yy];
   If[MyRank[(##*av) & /@ yy] < Total[av], Continue[]];
   aa = Table[0, {Length[yy]}];
   bb = Table[0, {Length[yy]}];
   If[Length[yy] == 1, Continue[]];
   For[j = 1, j <= Length[yy], j++,
    facet = Delete[##, Position[av, 0]] & /@ yy;
    newpoint = facet[[j]];
    facet = Delete[facet, j];
    pr = MyNullSpace[(## - facet[[1]]) & /@ facet];
    If[Length[pr] === 1,
     vector = pr[[1]];
     If[vector.facet[[1]] < 0, vector = -vector];
  (*   If[Or @@ ((## < 0) & /@ vector), Continue[]];*)
     If[vector.newpoint >= vector.facet[[1]], aa[[j]] = 1(*,aa[[j]]=2*)];
      bb[[j]] = vector;
     ];
    ];
   pos = Flatten[Position[aa, 1]];
   If[Length[pos] == 0,
        pos = Flatten[Position[aa, 2]];
        If[Length[pos] == 0,
           (* RawPrintLn["No facet found"]; RawPrintLn[yy]; Abort[]*)
           RawPrintLn[xxx];
           Return[False];
        ];
   ];

   j=Last[pos];
  (*j = Sort[pos, (Total[bb[[#1]]] > Total[bb[[#2]]]) &][[1]];*)
  (* j = Sort[pos, (VP[bb[[#1]]] < VP[bb[[#2]]]) &][[1]];*)
   yy = Delete[yy, j];
   ];
   If[Length[yy] == 1, Return[False]];
  facet = Delete[##, Position[av, 0]] & /@ yy;
  pr = MyNullSpace[(## - facet[[1]]) & /@ facet];
  pr = Normal[
      SparseArray[
       Apply[Rule, Transpose[{Flatten[Position[av, 1]], ##}], {1}],
       n]] & /@ pr;
  pr = pr[[1]];
  pr = If[## < 0, 0, ##] & /@ pr;
  If[Total[pr] == 1, av[[Flatten[Position[pr, 1]]]] = 0;
   Goto[avlabel]];
  pr
  ]

VP[xx_]:=(Length[Select[xx,(## > 0) &]])*10000+Total[xx]

 NormShift[xx_] := Module[{temp},
  temp = Apply[Min, Transpose[xx], {1}];
  (## - temp) & /@ xx
  ]

(*WD[aaa_,{}]:=Infinity;
WTilda[{}]:={};*)
  WMinVector[xx_] := Apply[Min, Transpose[xx], {1}];
WTilda[xx_] := Module[{temp = WMinVector[xx]}, (## - temp) & /@ xx];
WD[ind_, xx_] := Min @@ (ind.## & /@ xx);
WBSequence[xx_] :=
 Module[{n = Length[xx[[1]]], result = {}, ind, delta, d, H, temp,
   indN},
  ind = Table[1, {n}];
  delta = xx;
  result = {};
  While[True,
   AppendTo[result, {ind, delta}];
   If[delta == {}, Break[]];
   temp = WTilda[delta];
   d = WD[Table[1, {n}], temp];
   If[d == 0, Break[]];
   temp = Select[temp, ((ind.##) === d) &];
   H = If[##>0,1,0] & /@ Apply[Max, Transpose[temp], {1}];
   ind = Max[##, 0] & /@ (ind - H);
   temp = Union[(##/d) & /@ WTilda[delta], delta];
   temp = Select[temp, (H.## < 1) &];
   delta = Union[((##*ind)/(1 - ##.H)) & /@ temp];
   ];
  result
  ]

  WBNSequence[xx_]:= Module[{n = Length[xx[[1]]],temp=WBSequence[xx]},
    Append[Flatten[{Length[##[[1]]],WD[Table[1, {n}], WTilda[##[[2]]]]}&/@temp],WD[Table[1, {n}],Last[temp][[2]]]]
  ]

  WFindA[ind_, delta_] := Module[{temp, i, n, restart, temp2},
  n = Length[ind];
  temp = ind;
  Label[restart];
  For[i = 1, i <= n, i++,
   If[temp[[i]] == 1,
    temp2 = temp;
    temp2[[i]] = 0;
    If[And @@ ((temp2.## >= 1) & /@ delta),
     temp = temp2;
     Goto[restart];
     ]
    ]
   ];
  temp
  ]
WBSet[xx_] := Module[{n = Length[xx[[1]]], temp, ind, delta},
  temp = WBSequence[xx];
  {ind, delta} = Last[temp];
  If[delta === {}, Return[Table[1, {n}] - ind]];
  Return[Table[1, {n}] - ind + WFindA[ind,delta]]
  ]



MakeOneStep[xx_,facet_]:=MakeOneStep[xx,facet,Table[xx[[3]],{Length[facet]-Length[Position[facet,0]]}]]

MakeOneStep[xx_,facet_,graphs_]:=Module[{res,n},

n=Length[xx[[1]][[1]]];
  active = Complement[Range[n], Flatten[Position[facet, 0]]];
  res = {};
  For[i = 1, i <= Length[active], i++,
        Matrix=ReplacePart[IdentityMatrix[n],facet,active[[i]]];
        AppendTo[res, {(Matrix.##) & /@ xx[[1]],Matrix.xx[[2]],graphs[[i]]}];
  ];
  res = {NormShift[#[[1]]],#[[2]],#[[3]]}&/@ res;
  res = {OnlyLowPoints[#[[1]]],#[[2]],#[[3]]}&/@ res;
    res
];


MLN[xx_] := Module[{n, vectors, i, j, ll, nn,l,v,m},
  n = Length[xx];
  If[n===1,Return[{0,0,0}]];
  vectors =
   DeleteCases[
    Flatten[Table[xx[[i]] - xx[[j]], {i, 1, n}, {j, 1, n}], 1],
    Evaluate[Table[0, {Length[xx[[1]]]}]]];
  ll = (Max @@ ## - Min @@ ##) & /@ vectors;
  nn = (Length[Cases[##, Max @@ ##]] +
       Length[Cases[##, Min @@ ##]]) & /@ vectors;
  l=Min @@ ll;
  m=Min @@ ((##[[2]])&/@Cases[Transpose[{ll,nn}],{l,a_}]);

(*  m=Min@@Table[If[ll[[i]]===l,nn[[i]],Infinity],{i,Length[vectors]}];*)
  v=vectors[[Position[Transpose[{ll,nn}],{l,m}][[1]][[1]]]];
  {n, l, m, Position[v,Min@@v][[1]][[1]], Position[v,Max@@v][[1]][[1]]}
]

TriLower[xx_, yy_] :=
 If[xx[[1]] < yy[[1]], True,
  If[xx[[2]] < yy[[2]], True, If[xx[[3]] < yy[[3]], True, False],
   False], False]

 MakeOneStep[xxx_] := Module[{xx,temp, active, n,m, sbs,i,j, res,Matrix,newfacet,restarted=False},
 xx={xxx[[1]][[1]],xxx[[2]],xxx[[3]]};
 If[STRATEGY==STRATEGY_SS,
    n=Length[xx[[1]][[1]]];
    temp=GraphStep[xxx[[3]]];
    If[temp===True,
        Print["ERROR: strategy failed"];
        Print[{{{{xx[[1]][[1]]}},xx[[2]],xx[[3]]}}];
        Print[GetEdgeWeights[xxx[[3]]]];
        Abort[];
        Return[{{{{xx[[1]][[1]]}},xx[[2]],xx[[3]]}}];
    ];
    res=MakeOneStep[xx,Normal[SparseArray[Apply[Rule, (({##, 1}) & /@ temp[[1]]), {1}],n]],temp[[2]]];
  (*  Print[temp[[1]]];
    If[temp[[1]]==={2,5,3,4},res={res[[4]]}];*)
    Return[({{##[[1]]},##[[2]],##[[3]]}&/@res)]
 ];
 If[STRATEGY==STRATEGY_B,
    res=MakeOneStep[xx,WBSet[xx[[1]]]];
    Return[{{##[[1]]},##[[2]],##[[3]]}&/@res]
 ];
 If[STRATEGY==STRATEGY_S,
  n=Length[xx[[1]][[1]]];
  facet = TimeConstrained[Facet[xx[[1]]],1000,False];
  If[Not[facet===False],res=MakeOneStep[xx,facet]];
  If[Or[facet===False,Not[And@@((Length[##[[1]]]<Length[xx[[1]]])&/@res)]],
      {m0,l0,n0,vmin,vmax}=MLN[xx[[1]]];
      If[Or[facet===False,Not[And@@(TriLower[MLN[##[[1]]],{m0,l0,n0}]&/@res)]],
        facet=Normal[SparseArray[(##->1)&/@{vmin,vmax},n]];
        res=MakeOneStep[xx,facet];
        If[And@@(TriLower[MLN[##[[1]]],{m0,l0,n0}]&/@res),
            If[NumberQ[StopCounter],StopCounter++];
            If[StopCounter<0,Abort[]];
            Return[{{##[[1]]},##[[2]],##[[3]]}&/@res]
        ];
        RawPrintLn["Failed to resolve"];
        RawPrintLn[xx[[1]]];
        Abort[];
      ];
  ];
  Return[{{##[[1]]},##[[2]],##[[3]]}&/@res]
  ];



  If[STRATEGY==STRATEGY_A,
    n=Length[xx[[1]][[1]]];
    {m0,l0,n0,vmin,vmax}=MLN[xx[[1]]];
    facet=Normal[SparseArray[(##->1)&/@{vmin,vmax},n]];
    res=MakeOneStep[xx,facet];
    If[And@@(TriLower[MLN[##[[1]]],{m0,l0,n0}]&/@res),
         Return[{{##[[1]]},##[[2]],##[[3]]}&/@res]
    ];
    RawPrintLn["Failed to resolve"];
    RawPrintLn[xx[[1]]];
    Abort[];
  ];
  If[STRATEGY==STRATEGY_X,
    n=Length[xx[[1]][[1]]];
    For[j=1,j<=Length[xxx[[1]]],j++,

       For[m=2,m<=n,m++,
            sbs=Subsets[Range[n],{m}];
            For[i=1,i<=Length[sbs],i++,
                facet=Normal[SparseArray[(##->1)&/@sbs[[i]],n]];
                If[And@@((Not[(facet.##)===0])&/@(xxx[[1]][[j]])),
                    res=MakeOneStep[##,facet]&/@({##,xxx[[2]],xxx[[3]]}&/@(xxx[[1]]));
                    res=Transpose[res];
                    res=Transpose/@res;
                    res={##[[1]],##[[2]][[1]],##[[3]][[1]]}&/@res;
                    Return[res]
                ];
            ]
        ];


     ];


  ];




]











MyNullSpace[xx_]:=Module[{temp},
    If[xx==={},Return[{{1}}]];
    NullSpace[xx]
]

FindSD[xxx_]:=FindSD[xxx,False];
FindSD[xxx_,Extra_]:=Module[{temp,i}, yy = {};
        
  xx = {{OnlyLowPoints/@xxx,IdentityMatrix[Length[xxx[[1]][[1]]]],Extra}};
   temp=xx;
   yy = Join[yy, Select[temp, (Times@@(Length/@(##[[1]])) === 1) &]];
   xx = Select[temp, Not[Times@@(Length/@(##[[1]])) === 1] &];
                             
  While[True,
   temp = MakeOneStep /@ xx;
   temp = Flatten[temp, 1];
   yy = Join[yy, Select[temp, (Times@@(Length/@(##[[1]])) === 1) &]];
   xx = Select[temp, Not[Times@@(Length/@(##[[1]])) === 1] &];

   If[Length[xx] === 0,
       yy=#[[2]]&/@yy;
       Return[yy]
   ];
   ];

 ]

(* Graph compinatorics *)

MyGraph[edges_] := Module[{temp},
  temp = FromUnorderedPairs[edges];
  Graph[Transpose[{Edges[temp], (EdgeWeight -> ##) & /@
      Range[Length[Edges[temp]]]}], {##} & /@ Vertices[temp]
   ]
  ]
MyEdges[graph_] := {##[[1]], ##[[2]][[2]]} & /@ graph[[1]];
MyDeleteEdge[graph_, weight_] :=
  Graph[DeleteCases [graph[[1]], {aa_, EdgeWeight -> weight}],
   graph[[2]]];
MyDeleteEdges[graph_, {}] := graph;
MyDeleteEdges[graph_, weights_] :=
  Graph[Intersection @@ (DeleteCases [
        graph[[1]], {aa_, EdgeWeight -> ##}] & /@ weights),
   graph[[2]]];
MyBiconnectedComponents[graph_] :=
 InduceSubgraph[graph, DeleteCases[##,InfinityVertex]] & /@ BiconnectedComponents[graph]
CleanGraph[graph_] := Module[{g, temp},
  g = RemoveSelfLoops[graph];
  temp = MyBiconnectedComponents[g];
  temp = Select[temp, (M[##] <= 1) &];
  temp = Union @@ (GetEdgeWeights /@ temp);
  MyDeleteEdges[g, temp]
  ]

ExternalNotConnected[graph_, external_] := Module[{temp},
  temp = ConnectedComponents[graph];
  temp = Union[(##[[1]]) & /@ Position[temp, ##] & /@ external];
  Or[temp === {{}}, Length[temp] =!= 1]
]

GraphStep[graph_] := Module[{g, temp},
  g = CleanGraph[graph];
  If[M[g] === Exte, Return[True]];
  If[Not[ExternalNotConnected[DeleteVertex[graph,InfinityVertex],External]],
      temp = Sort[DeleteCases[##, InfinityVertex] & /@ BiconnectedComponents[g], (Length[#1] > Length[#2]) &];
  ,
      temp = Sort[BiconnectedComponents[DeleteVertex[graph,InfinityVertex]], (Length[#1] > Length[#2]) &];
  ];
  temp = GetEdgeWeights[InduceSubgraph[g,temp[[1]]]];
  temp=Sort[temp];
  Return[{temp, (MyDeleteEdge[g, ##] & /@ temp)}];
]
ContractEdge[graph_, weight_] := Module[{temp},
  temp =
   Rule @@ Cases [graph[[1]], {aa_, EdgeWeight -> weight}][[1]][[1]];
  If[temp[[1]] === temp[[2]], Return[False]];
  MyDeleteEdge[
   Graph[ReplacePart[##, ##[[1]] /. temp, 1] & /@ graph[[1]],
    graph[[2]]], weight]
  ]
ContractEdges[graph_, weights_] := Module[{temp, i},
  temp = graph;
  For[i = 1, i <= Length[weights], i++,
   If[temp === False, Return[False]];
   temp = ContractEdge[temp, weights[[i]]];
   ];
  temp
  ]


(* Graph compinatorics end *)




ep:=Global`ep;
QOpen[xx__]:=Global`QOpen[xx];
QRead[xx__]:=Global`QRead[xx];
QRemoveDatabase[xx__]:=Global`QRemoveDatabase[xx];
QRepair[xx__]:=Global`QRepair[xx];
QClose[xx__]:=Global`QClose[xx];
QList[xx__]:=Global`QList[xx];
QSize[xx__]:=Global`QSize[xx];
QPut[xx__]:=Global`QPut[xx];
QGet[xx__]:=Global`QGet[xx];
QSafeGet[xx__]:=Global`QSafeGet[xx];
QCheck[xx__]:=Global`QCheck[xx];
QRemove[xx__]:=Global`QRemove[xx];
CIntegrate[xx__]:=Global`CIntegrate[xx];
SetPoints[xx__]:=Global`SetPoints[xx];
SetCut[xx__]:=Global`SetCut[xx];
CAddString[xx__]:=Global`CAddString[xx];
CClearString[xx__]:=Global`CClearString[xx];



End[];

SDEvaluate[x__]:=FIESTA`SDEvaluate[x]
SDEvaluateG[x__]:=FIESTA`SDEvaluateG[x]
UF[x__]:=FIESTA`UF[x]
ClearResults[]:=FIESTA`ClearResults[]
