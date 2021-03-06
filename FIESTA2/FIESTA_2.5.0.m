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
Group;
<<Combinatorica`;
Unprotect[Combinatorica`EmptyQ];
Remove[Combinatorica`EmptyQ];
Unprotect[SeriesCoefficient];
SeriesCoefficient[x_, y_Integer] := 0 /; And[NumberQ[x], y =!= 0]
pm;
ep;
delta;
RemoteLinkName;
RemoteLinksFile;
StartingStage;
AbortStage;
KnownIntegrals;


SetSystemOptions[CacheOptions -> Symbolic -> False];
SetSystemOptions[CacheOptions -> Numeric -> False];
SetSystemOptions[CacheOptions -> Constants -> False];
SetSystemOptions[CacheOptions -> CacheKeyMaxBytes -> 1];
SetSystemOptions[CacheOptions -> CacheResultMaxBytes -> 1];


If[Not[ValueQ[RestartLimit]],
    RestartLimit=8;
];


If[Not[ValueQ[BisectionVariables]],
    BisectionVariables={};
];

If[Not[ValueQ[NegativeTermsHandling]],
    NegativeTermsHandling="Squares";
]; (* "AdvancedSquares", "AdvancedSquares2", "AdvancedSquares3", "Squares", "None", "MB" *)

If[Not[ValueQ[ExactIntegrationOrder]],
    ExactIntegrationOrder=-Infinity;
];

If[Not[ValueQ[ExactIntegrationTimeout]],
    ExactIntegrationTimeout=10;
];

If[Not[ValueQ[PrintInsteadOfIntegrating]],
    PrintInsteadOfIntegrating=False;
];

If[Not[ValueQ[ResolutionMode]],
    ResolutionMode="Taylor";
]; (* "Taylor", "IBP0", "IBP1" *)




If[Not[ValueQ[RemoteLinkTimeout]],
    RemoteLinkTimeout=500;
];
If[Not[ValueQ[RemoteLinkInstallTimeout]],
    RemoteLinkInstallTimeout=20;
];

If[Not[ValueQ[NumberOfSubkernels]],
    NumberOfSubkernels=0;
];

If[Not[ValueQ[ReturnErrorWithBrackets]],
    ReturnErrorWithBrackets=False
];

If[Not[ValueQ[MemoryDebug]],
    MemoryDebug=10000000
];
If[Not[ValueQ[IntegrationDebug]],
    IntegrationDebug=False;
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
    CIntegratePath="CIntegrateMP";
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
If[Not[ValueQ[IntegrationCut]],
    IntegrationCut=0;
];
If[Not[ValueQ[ZIntegrationCut]],
    ZIntegrationCut=0;
];
If[Not[ValueQ[MixSectors]],
    MixSectors=0;
];
If[Not[ValueQ[CurrentIntegrator]],
    CurrentIntegrator="vegasCuba";
]




DefaultPrecision;
If[Not[ValueQ[SmallX]],
    SmallX=0.00001;   (* evaluation should be possible at {SmallX,...,SmallX} *)
];


If[ValueQ[MPThreshold],MPThreshhold=MPThreshold;];
If[Not[ValueQ[MPThreshhold]],
     MPThreshhold=10^-9;    (*  going to high precision if the Monom is less than MPThreshhold *)
];

MPMin;
If[Not[ValueQ[PrecisionShift]],
    PrecisionShift=38;   (* additional bits in multiprecision *)
];


STRATEGY_A;
STRATEGY_B;
STRATEGY_S;
STRATEGY_X;
STRATEGY_0;
STRATEGY_SS;
STRATEGY_KU0;
STRATEGY_KU;
STRATEGY_KU2;
PrimarySectorCoefficients;
PolesMultiplicity;

If[Not[ValueQ[STRATEGY]],
    STRATEGY=STRATEGY_S;
]; (*might be STRATEGY_0, STRATEGY_A, STRATEGY_B, STRATEGY_S, STRATEGY_X, STRATEGY_KU0, STRATEGY_KU, STRATEGY_KU2*)

If[Not[ValueQ[AbortOnNaN]],
    AbortOnNaN=False;
];
If[Not[ValueQ[PMCounter]],
    PMCounter=1;
];

If[Not[ValueQ[RemoveDatabases]],
    RemoveDatabases=True;
];

If[Not[ValueQ[QHullPath]],
    QHullPath="./qhull";
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
    RawPrintLn["Ensure that you are using the multi-precision CIntegrate binary."];
    RawPrintLn["-----------------------------------------"];
    RawPrintLn["If none of the above points help, please contact the authors and provide your example."];
]


Begin["FIESTA`"];

SDEvaluateG::Usage="SDEvaluate[{graph,external},{U,F,loops},indices,order] evaluates the integral"
SDEvaluate::Usage="SDEvaluate[{U,F,loops},indices,order] evaluates the integral"
SDEvaluateDirect::Usage="SDEvaluate[degrees,functions,order] evaluates the integral (no delta function)"
SDExpandG::Usage="SDExpandG[{graph,external},{U,F,loops},indices,order,expand_var,expand_degree] expands the integral"
SDExpand::Usage="SDExpand[{U,F,loops},indices,order,expand_var,expand_degree] evaluates the integral"
ClearResults::Usage="ClearResults[] clears the results from memory"
UF::Usage="UF[LoopMomenta,Propagators,subst] generates the functions U and F"

ZREPLACEMENT:={z->-0.2,z2->-0.3};

MyWorstPower[x_Power, z_] := 
  If[Head[x[[1]]] === y, Min[##, 0] & /@ Exponent[x, z], 
   Abs[x[[2]]]*MyWorstPower[x[[1]], z]];
MyWorstPower[x_Plus, z_] := 
  Apply[Min, Transpose[MyWorstPower[##, z] & /@ List @@ x], {1}];
MyWorstPower[x_Times, z_] := 
  Apply[Plus, Transpose[MyWorstPower[##, z] & /@ List @@ x], {1}];
MyWorstPower[x_Integer, z_] := 0 & /@ z;
MyWorstPower[x_Rational, z_] := 0 & /@ z;
MyWorstPower[x_Log, z_] := 0 & /@ z;
MyWorstPower[x_y, z_] := 0 & /@ z;
MyWorstPower[x_, {}] := {};
MyWorstPower[EulerGamma, z_] := 0 & /@ z;
MyWorstPower[Pi, z_] := 0 & /@ z;
MyWorstPower[PolyGamma[_,_], z_] := 0 & /@ z;
MyWorstPower[Gamma[_], z_] := 0 & /@ z;



ZeroCheck[x_] := Module[{vars, i, temp},
  vars = Variables[x];
  For[i = 1, i <= 10, i++,
   temp = 
    Expand[x /. (Rule[##, 1/RandomInteger[{1, 10}]] & /@ vars)];
   If[temp =!= 0, Return[False]];
   ];
  True
  ]




MyToString[x__]:=If[UsingQLink,ToString[x,InputForm],x]
MyToExpression[x__]:=If[UsingQLink,ToExpression[x],x]
MyStringLength[x__]:=If[UsingQLink,StringLength[x],ByteCount[x]]

MyMemoryInUse[]:=If[NumberOfSubkernels>0,MaxMemoryUsed[]+Plus@@(ParallelEvaluate[MaxMemoryUsed[]]),MaxMemoryUsed[]]

BuryLink[i_]:=Module[{temp},
    If[MyLinkType[i]=="LocalCIntegrate",
        RawPrintLn["Dead local link!"];	
    ];
    If[MyLinkType[i]=="RemoteCIntegrate",
        RawPrintLn["Link at ",MyLinkPath[i]," is dead!"];
    ];    
    MyLinkType[i]="DeadLink";
]

BanLink[i_]:=Module[{temp},
    If[MyLinkType[i]=="LocalCIntegrate",
        RawPrintLn["Frozen local link!"];
    ];
    If[MyLinkType[i]=="RemoteCIntegrate",
        RawPrintLn["Link at ",MyLinkPath[i]," is frozen!"];
    ];
    MyLinkBannedUntill[i]=AbsoluteTime[]+(60*(2^MyLinkBanCounter));
    MyLinkBanCounter[i]++;
]


MyLinkRead[i_]:=Module[{temp},
    If[MyLinkBannedUntill[i]>0,Return[$Failed]];
    If[MyLinkType[i]=="DeadLink",Return[$Failed]];
    If[MyLinkType[i]=="LocalCIntegrate",
        temp=TimeConstrained[Quiet[LinkRead[MyLinks[[i]]]],RemoteLinkTimeout,$Timeout];
    ];
    If[MyLinkType[i]=="RemoteCIntegrate",
        temp=TimeConstrained[Quiet[LinkRead[MyLinks[[i]]]],RemoteLinkTimeout,$Timeout];
    ];
    If[temp===$Failed,BuryLink[i]];
    If[temp===$Timeout,BanLink[i];temp=$Failed];
    Return[temp];
]

MyLinkWrite[i_,string_]:=MyLinkWrite[i,string,False]
MyLinkSubmit[i_,string_]:=MyLinkWrite[i,string,True]

MyLinkWrite[i_,string_,submit_]:=Module[{temp},
    If[MyLinkBannedUntill[i]>0,Return[$Failed]];
    If[MyLinkType[i]=="DeadLink",Return[$Failed]];
    If[MyLinkType[i]=="MainKernel",
        If[submit,
            temp=IntegrateHere[{MyLinkBuffer[i]},CurrentIntegralExact];
            If[temp[[3]]>0,CurrentIntegralExact=False];
            temp=(temp[[1]]+temp[[2]]*ToExpression["pm"<>ToString[PMCounter++]]);
            result3=result3+temp;
            MyLinkBuffer[i]={},
            AppendTo[MyLinkBuffer[i],MyToExpression[string]];
        ];
    ];
    If[MyLinkType[i]=="SubKernel",
        If[submit,
            temp={MyLinkBuffer[i]};
            MyLinkBuffer[i]=ParallelSubmit[{temp,CurrentIntegralExact},IntegrateHere[temp,CurrentIntegralExact]],
            AppendTo[MyLinkBuffer[i],MyToExpression[string]];
        ];
    ];
    If[MyLinkType[i]=="LocalCIntegrate",
        temp=TimeConstrained[Quiet[LinkWrite[MyLinks[[i]],CallPacket[If[submit,0,2],{RemoveFIESTAName[string]}]]],RemoteLinkTimeout,$Timeout];
        If[IntegrationDebug,
            AppendTo[MyLinkBuffer[i],RemoveFIESTAName[string]];
        ];
    ];
    If[MyLinkType[i]=="RemoteCIntegrate",
        temp=TimeConstrained[Quiet[LinkWrite[MyLinks[[i]],CallPacket[If[submit,0,2],{RemoveFIESTAName[string]}]]],RemoteLinkTimeout,$Timeout];
    ];
    If[temp===$Failed,BuryLink[i]];
    If[temp===$Timeout,BanLink[i];temp=$Failed];
    Return[temp];
]


MyLinkReadyQ[i_]:=Module[{temp},
    If[MyLinkType[i]=="MainKernel",Return[False]];
    If[MyLinkType[i]=="SubKernel",Return[False]];
    If[MyLinkBannedUntill[i]>0,Return[False]];
    If[MyLinkType[i]=="DeadLink",Return[False]];
    If[MyLinkType[i]=="LocalCIntegrate",
        temp=TimeConstrained[Quiet[LinkReadyQ[MyLinks[[i]]]],RemoteLinkTimeout,$Timeout];
    ];
    If[MyLinkType[i]=="RemoteCIntegrate",
        temp=TimeConstrained[Quiet[LinkReadyQ[MyLinks[[i]]]],RemoteLinkTimeout,$Timeout];
    ];
    If[temp===$Failed,BuryLink[i]];
    If[temp===$Timeout,BanLink[i];temp=$Failed];
    Return[temp];
]

MyLinkWriteOption[i_,option_,value_]:=Module[{temp},
    If[MyLinkType[i]=="DeadLink",Return[$Failed]];
    If[MyLinkType[i]=="LocalCIntegrate",
        temp=Quiet[LinkWrite[MyLinks[[i]],CallPacket[option,value]]];
    ];
    If[MyLinkType[i]=="RemoteCIntegrate",
        temp=TimeConstrained[Quiet[LinkWrite[MyLinks[[i]],CallPacket[option,value]]],RemoteLinkTimeout,$Failed];
    ];
    If[temp===$Failed,BuryLink[i]];
    Return[temp];
]

(*
option 1: IntegrationCut {0,0.1}
option 4: Integrator {"vegasf"}
option 5: Integrator option {"key","value"}
option 6: request to send back integrator settings {}
option 7: MPPrecision
option 8: SmallX
option 9: MPThreshhold
option 10: MPMin
option 11: MPPrecisionShift
*)

FreeLink[]:=Module[{temp},
    temp=If[Or[MyLinkType[##]==="DeadLink",MyLinkBannedUntill[##]>0],Infinity,MyLinkTask[##]]
                    &/@Range[Length[MyLinks]];
    Max@@(Append[Flatten[Position[temp,0]],0])
]
(* ready to accept a new task *)

ReadyLink[Final_]:=Module[{temp,i},
    For[i=1,i<=Length[MyLinks],i++,
        If[MyLinkType[i]==="DeadLink",Continue[]];
        If[MyLinkBannedUntill[i]>0,Continue[]];
        If[MyLinkReadyQ[i],Return[i]];
        If[And[Final,MyLinkType[i]=!="LocalCIntegrate",AbsoluteTime[]-MyLinkStart[i]>10*MaxEvaluationTime],
            BanLink[i];
            Return[i];
        ];
    ];
    Return[0];
        Max@@(Append[Flatten[Position[
                If[MyLinkType[##]==="DeadLink",Infinity,MyLinkReadyQ[##]]
                    &/@Range[Length[MyLinks]],
            True]],0])
]
(* has something that we should read from it *)

MyLinkInitialize[i_]:=Module[{temp},

    If[Not[MemberQ[Links[],MyLinks[[i]]]],  (* installing *)
        If[MyLinkType[i]=="LocalCIntegrate",
            MyLinks[[i]] = Install[CIntegratePath];
        ];
        If[MyLinkType[i]=="RemoteCIntegrate",
            MyLinks[[i]] = TimeConstrained[Quiet[Install[LinkConnect[MyLinkPath[i],LinkProtocol->"TCPIP"]]],RemoteLinkInstallTimeout,$Failed];
            If[MyLinks[[i]] === $Failed,
                Print["Link at ",MyLinkPath[i]," could not be started"];
                MyLinkType[i]="DeadLink"
            ];
        ];
    ];
    If[And[MyLinkType[i]=!="DeadLink",Not[Head[MyLinks[[i]]]===LinkObject]],
        RawPrintLn["Could not start external integration link"];
        Abort[];
    ];
    Quiet[While[MyLinkReadyQ[i],
        MyLinkRead[i]
    ]];
    MyLinkStart[i]=0;
    MyLinkTask[i]=0;
    MyLinkBanCounter[i]=0;
    MyLinkBannedUntill[i]=0;
    MyLinkBuffer[i]={};
]

AddLocalLink[]:=Module[{temp},
    AppendTo[MyLinks,Null];
    MyLinkType[Length[MyLinks]]="LocalCIntegrate";
]


AddRemoteLink[name_]:=Module[{temp},
    AppendTo[MyLinks,Null];
    MyLinkType[Length[MyLinks]]="RemoteCIntegrate";
    MyLinkPath[Length[MyLinks]]=name;
]

SetLinkOptions[]:=Module[{i},
  For[i=1,i<=Length[MyLinks],i++,
    SetOneLinkOptions[i];
  ]
]

SetOneLinkOptions[i_]:=Module[{temp},
            If[MyLinkType[i]=="DeadLink",Return[]];
            If[IntegrationCut>0,
                MyLinkWriteOption[i,1,{0,IntegrationCut}];
                MyLinkRead[i];
            ];
            MyLinkWriteOption[i,4,{CurrentIntegrator}];
            MyLinkRead[i];
            If[ValueQ[CurrentIntegratorSettings],
               For[j=1,j<=Length[CurrentIntegratorSettings],j++,
                    MyLinkWriteOption[i,5,CurrentIntegratorSettings[[j]]];
                        If[i == 1,
                            RawPrintLn[CurrentIntegratorSettings[[j]],
                            ": ",MyLinkRead[i]];,
                          MyLinkRead[i]
                       ];
               ];
            ];
           If[ValueQ[MPPrecision],
                MyLinkWriteOption[i,7,{MPPrecision}];
                MyLinkRead[i];
           ];
           MyLinkWriteOption[i,8,{SmallX}];
           MyLinkRead[i];
           MyLinkWriteOption[i,9,{N[MPThreshhold]}];
           MyLinkRead[i];
           If[ValueQ[MPMin],
                MyLinkWriteOption[i,10,{MPMin}];
                MyLinkRead[i];
           ];
           MyLinkWriteOption[i,11,{PrecisionShift}];
           MyLinkRead[i];        
]


ReadLinkOptions[i_]:=Module[{temp},
            RawPrintLn["CurrentIntegrator=",CurrentIntegrator,";"];
            MyLinkWriteOption[i,6,{}];
            RawPrintLn["CurrentIntegratorSettings=\n",MyLinkRead[i],";"];
            If[ValueQ[CurrentIntegratorSettings],
               Clear[CurrentIntegratorSettings];
            ];
]


InitializeLinks[]:=Module[{temp,i,j},
    Print[VersionString];
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
        ];
        Print["Database ready"];
    ];
    If[UsingC,
        If[Head[MyLinks]===Symbol,
            MyLinks={};
            Do[AddLocalLink[],{NumberOfLinks}];
            If[Head[RemoteLinkName]===String,AddRemoteLink[RemoteLinkName]];
            If[Head[RemoteLinksFile]===String,
                Quiet[
                    temp=ReadList[RemoteLinksFile, Record, RecordSeparators -> "||"];
                    DeleteFile[RemoteLinksFile];
                ];
                If[Head[temp]===Symbol,Print["Cannot read from ",RemoteLinksFile]];
                If[Length[temp]>0,
                    Print[Length[temp]," remote links found."];
                    AddRemoteLink/@temp;
                ];
            ];
        ];
        For[i=1,i<=Length[MyLinks],i++,
            MyLinkInitialize[i]
        ];

        SetLinkOptions[];
        ReadLinkOptions[1];

    ,
        MyLinks={};
        If[NumberOfSubkernels>1,
            For[i=1,i<=NumberOfSubkernels,i++,
                AppendTo[MyLinks,LinkObject[]];
                MyLinkType[Length[MyLinks]]="SubKernel";
            ];
        ,
            AppendTo[MyLinks,LinkObject[]];
            MyLinkType[Length[MyLinks]]="MainKernel";
        ];
        For[i=1,i<=Length[MyLinks],i++,
            MyLinkInitialize[i]
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
    If[NumberOfSubkernels===1,
        RawPrintLn["The number of Subkernels should be 0 or greater than 1"];
        RawPrintLn["NumberOfSubkernels set to 0"];
        NumberOfSubkernels=0;
    ];
    If[NumberOfSubkernels>0,
        RawPrintLn["Starting ",NumberOfSubkernels," Subkernels"];	
	If[Length[Kernels[]]>0,
(*	  Print[1];*)
	  ParallelEvaluate[Abort[]];
(*Print[2];*)
Parallel`Developer`ClearKernels[];
(*Print[3];*)
	  Parallel`Developer`ResetQueues[];
(*Print[4];*)
	
	  CloseKernels[];	
	];
        LaunchKernels[NumberOfSubkernels];
        DistributeDefinitions["FIESTA`"];
        ParallelEvaluate[
            SetSystemOptions[CacheOptions -> Symbolic -> False];
            SetSystemOptions[CacheOptions -> Numeric -> False];
            SetSystemOptions[CacheOptions -> Constants -> False];
            SetSystemOptions[CacheOptions -> CacheKeyMaxBytes -> 1];
            SetSystemOptions[CacheOptions -> CacheResultMaxBytes -> 1];
        ];
    ];
    CodeInfo[];
]




UF[xx_, yy_, z_] := Module[{degree, coeff, i, t2, t1, t0, zz},
    zz = Map[Rationalize[##,0]&, z, {0, Infinity}];
    degree = -Sum[yy[[i]]*x[i], {i, 1, Length[yy]}];
    coeff = 1;
    For[i = 1, i <= Length[xx], i++,
        t2 = Coefficient[degree, xx[[i]], 2];
        t1 = Coefficient[degree, xx[[i]], 1];
        t0 = Coefficient[degree, xx[[i]], 0];
        coeff = coeff*t2;
        degree = Together[t0 - ((t1^2)/(4 t2))];
    ];
    degree = Together[-coeff*degree] //. zz;
    coeff = Together[coeff] //. zz;
    {coeff, Expand[degree], Length[xx]}
]

ClearResults[]:=Clear[SDIntegral,SDExact];

VersionString:="FIESTA 2.5.0";
CodeInfo[]:=Module[{temp},
    RawPrintLn["UsingC: ",UsingC];
    RawPrintLn["NumberOfLinks: ",NumberOfLinks];
    RawPrintLn["UsingQLink: ",UsingQLink];
    RawPrintLn["IntegrationCut: ",IntegrationCut];
    RawPrintLn["Strategy: ",STRATEGY];
    If[And[Not[STRATEGY===STRATEGY_SS],GraphUsed===True],
        RawPrintLn["WARNING: a graph has been specified -> STRATEGY_SS recomended"];
    ];
]
Print[VersionString];



CutExtraDigits[xx_]:=Map[If[And[NumberQ[##],Not[IntegerQ[##]]], Chop[N[Round[##,10^-6]],10^-6], ##] &,xx,{0,Infinity}];
DoubleCutExtraDigits[xx_]:=Module[{temp},
    If[xx===Indeterminate,Return[Indeterminate]];
    temp=DeleteCases[DeleteCases[Variables[xx],ExV],Log[ExV]];
    If[Length[temp]>1,Return[xx]];
    If[Length[temp]==0,Return[CutExtraDigits[xx]]];
    CutExtraDigits[xx/.(temp[[1]]->0)]+CutExtraDigits[(xx-(xx/.(temp[[1]]->0)))/.(temp[[1]]->1)]*temp[[1]]
]

PMSymbol:=If[$FrontEnd===Null,"+-","±"]

PMForm[xx_]:=Module[{temp},
    If[xx===Indeterminate,Return["INDETERMINATE"]];
    temp=DeleteCases[DeleteCases[Variables[xx],ExV],Log[ExV]];
    If[Length[temp]>1,Return[ToString[xx,InputForm]]];
    If[Length[temp]==0,Return[ToString[xx,InputForm]]];
    ToString[xx/.(temp[[1]]->0),InputForm]<>" "<>PMSymbol<>" "<>ToString[(xx-(xx/.(temp[[1]]->0)))/.(temp[[1]]->1),InputForm]
]







RawPrint[x__]:=WriteString[$Output,x];
RawPrintLn[x__]:=WriteString[$Output,x,"\n"];
RawPrintLn[]:=WriteString[$Output,"\n"];
PrepareDots[n_]:=Module[{temp},gn=n;gl={};For[gi=1,gi<=10,gi++,AppendTo[gl,Quotient[gi*gn,10]]];gi=0;]
Dots[]:=(gi++;While[And[Length[gl]>0,gi>=gl[[1]]],gl=Drop[gl,1];RawPrint["."]];)


CC[m_]:=(If[Mod[gi,m]===0,ClearSystemCache[]])
MyTimingForm[xx_]:=ToString[Chop[N[Round[xx,10^-4]],10^-4],InputForm]


If[Not[TrueQ[$VersionNumber>=6.0]],
  Print["WARNING! Full compatibility with Mathematica versions lower than 6.0 not guaranteed!"];
];





FactorMonom[xxx_]:=Module[{xx,vars,i,monom,temp},
    xx=xxx;
    If[xx===0,Return[{0,0}]];
    vars=Sort[Cases[Variables[xx],y[_]]];
    temp=-Exponent[xx,##^-1]&/@vars;
    monom=Times@@(vars^temp);
    {monom,Expand[xx/(monom)]}
]

FactorLogMonom[xxx_]:=Module[{xx,vars,i,monom,temp},
    xx=xxx;
    If[xx===0,Return[{0,0}]];
    vars=Sort[Cases[Variables[xx],Log[y[_]]]];
    temp=-Exponent[xx,##^-1]&/@vars;
    monom=Times@@(vars^temp);
    {monom,Expand[xx/(monom)]}
]


ConstructTerm[xx_]:=Module[{temp,i},
                temp=xx;
                If[ExpandMode,temp[[1]]=Drop[temp[[1]],-1];temp[[2]]=Drop[temp[[2]],-1];temp[[3]]=Drop[temp[[3]],-1]];
               temp=  {Flatten[Position[Table[If[FreeQ[temp[[2]][[i]],z],temp[[1]][[i]],0],{i,Length[temp[[1]]]}],2]],
                        Inner[(#1^#2)&,Array[x,Length[temp[[1]]]],temp[[2]],Times]*
                        Inner[(Log[#1]^#2)&,Array[x,Length[temp[[1]]]],temp[[3]],Times]*
                        xx[[4]],If[ExpandMode,{Last[xx[[2]]],Last[xx[[3]]]},{0,0}],
                        xx[[5]]
                        };
               temp
]



MyExpandDenominator[x_] := Module[{part1, part2,temp},
  If[Head[x] === Power,
      temp={x};
  ];
  If[Head[x] === Times,
      temp=List@@x;
  ];
  If[Head[temp]=!=List,
    Return[x];
  ];
  part1 = Select[temp, (And[Head[##] === Power, ##[[2]] < 0]) &];
  part2 = Select[temp, Not[TrueQ[(And[Head[##] === Power, ##[[2]] < 0])]] &];
  part2=Times@@part2;
  part1=ExpandAll[##,ep]&/@part1;
  part1=Times@@part1;
  Return[part1*part2];
]


EpPosExpand[xx_, extra_, inshift_, order_] := Module[{temp, shift, i, der, result},
    Check[shift = 0; temp = MyExpandDenominator[xx] /. MBshiftRules[ep];
        While[Quiet[Check[ep=0;Evaluate[temp],False]===False],
            Clear[ep];
            temp = Expand[temp*ep,ep];
            shift++;
            If[shift > 20,Print["big shift"]; Print[xx];
                Print[extra];
		Print[inshift];
                Print[order];
                Abort[];
            ];
        ];
        Clear[ep];
        result = Reap[
            Sow[{-shift + inshift, {temp, extra}}];
            der = {{temp, extra}};
            For[i = 1, i <= order -inshift + shift, i++,
                der = Flatten[If[
                    D[##[[2]][[1]], ep] =!= 0,
                    {{D[##[[1]], ep], ##[[2]]}, {##[[1]]*D[##[[2]][[1]], ep], {##[[2]][[1]], ##[[2]][[2]] + 1}}},
                    {{D[##[[1]], ep], ##[[2]]}}
                ] & /@ der, 1];
                Sow[{-shift + i + inshift, {##[[1]]/Factorial[i], ##[[2]]}}] & /@ der
            ]
        ][[2]][[1]];
        ep = 0;
        result = result;
        result=DeleteCases[result,{uu_,{0,vv_}}];
        Clear[ep];
   ,
   Print["ep pos expand"]; Print[xx];Print[extra];Print[inshift]; Print[order]; Abort[]];
   result
]




MyExponent[xxx_, yy_] := Module[{temp},
    xx=Expand[xxx/(xxx/.x[aaa_]->1)];
  temp = Exponent[xx /. ep -> 0, yy];
  Expand[(Exponent[xx /. ep -> 1, yy] - temp)ep + temp]
  ]



ExpandZ[xx_] :=
 Module[{temp, degs,i, n, j,baddegs,gooddegs,monom},
temp=xx;

If[temp[[1]]==={},Return[{xx}]];
   If[ExpandMode,temp={Append[temp[[1]],2],AppendTo[temp[[2]],(*If[Power2,z+2z2,z]*)Exponent[temp[[3]][[1]],ExV]],temp[[3]]/.ExV->1,temp[[4]],temp[[5]],Append[temp[[6]],0]}];

  temp = {temp[[1]], temp[[2]], ##, temp[[4]], {##,{},0}, temp[[5]],temp[[6]]} & /@ temp[[3]];
  n = Length[temp[[1]][[1]]];
  temp = (degs=##[[2]] + Table[MyExponent[PowerExpand[##[[3]]], x[j]], {j, n}];
            {baddegs,gooddegs}=Transpose[Table[If[(degs[[j]]/.Flatten[{ep->0,delta[_]->0,ZREPLACEMENT}])>=0,{0,degs[[j]]},{degs[[j]],0}], {j, n}]];
            monom=Inner[Power,x/@Range[n],gooddegs,Times];
     (*   {##[[1]], baddegs, (monom)*(##[[3]] /. x[aa_]->1),##[[4]],##[[5]]}*)
   {##[[1]], degs, (##[[3]] /. x[aa_]->1),##[[4]],##[[5]],##[[6]],##[[7]]}
         )& /@ temp;
   temp={##[[1]],##[[2]],##[[3]],##[[4]],{##[[5]][[1]],##[[5]][[2]],
        ExpandOrder(*-(Plus@@(If[Coefficient[##,z]>=0,0,(##+1)/.{ep->0,z->0}]&/@Union[##[[2]]]))*)
    },##[[6]],##[[7]]}&/@temp;
    If[Not[HASZ],temp={##[[1]],##[[2]],##[[7]],##[[3]]*##[[4]],##[[6]],##[[5]]}&/@temp];

  temp]

GroupTerms[xx_] := Module[{temp},
   temp =
    Reap[Sow[##[[3]], {{##[[1]], ##[[2]], ##[[4]]}}] & /@ xx, _,
      List][[2]];
   temp = {##[[1]][[1]], ##[[1]][[2]], (Plus @@ ##[[2]])*##[[1]][[3]]} & /@ temp;
   temp
   ];



FirstVars[xx__]:=Module[{vars},
  vars = Union[Cases[xx,x[aaa_],{0,Infinity}]];
  {xx/.Apply[Rule,Transpose[{vars,Array[y,Length[vars]]}],{1}],Length[vars]}
]

CountVars[xx__]:=Module[{vars},
  vars = Union[Cases[xx,y[aaa_],{0,Infinity}]];
  {xx,Max@@((##[[1]])&/@vars),(*Length[vars]*)10}
]


AdvancedFirstVars[{yyy_,xx_,zz_,uu_}]:=Module[{yy=x/@yyy,vars,rules},
   vars = Union[Cases[xx,x[aaa_],{0,Infinity}]];
  vars=Sort[vars,(If[TrueQ[MemberQ[yy,#1]],If[TrueQ[MemberQ[yy,#2]],#1[[1]]<#2[[1]],True],False])&];
  rules=Apply[Rule,Transpose[{vars,Array[y,Length[vars]]}],{1}];
  {Sort[(##[[1]])&/@(yy/.rules)],xx/.rules,zz,uu}
]









UnsortedUnion[xx_] := Reap[Sow[1, xx], _, #1 &][[2]]

Format[if[x1_, x2_, x3_], InputForm] :=
  "if(" <> ToString[x1, InputForm] <> ")>(f[1])(" <>
   ToString[x2, InputForm] <> ")(" <> ToString[x3, InputForm] <> ")";


MyString[xx_]:=StringReplace[ToString[xx/.{Power->p,Log->l,y->x,Pi->P,EulerGamma->G},InputForm],{"\\"->"","\""->""}];
(*MyString[xx_]:=StringReplace[ToString[xx//.{Power[a_,b_]:>p[a,If[IntegerQ[b],b,N[b]]]}/.{Log->l,y->x,Pi->P,EulerGamma->G},InputForm],{"\\"->"","\""->""}];*)

(*MyString[xx_] := 
  StringReplace[
   ToString[xx, InputForm], {"\\" -> "", "\"" -> "", "y" -> "x", 
    "Pi" -> "P"}];
*)

DoIntegrate[{hasz_,xx_},tryexact_] := Module[{temp2,vars, temp,i,rules,met,res},saved=xx;
vars =
Union[Cases[xx,y[aaa_],{0,Infinity}],Cases[xx,x[aaa_],{0,Infinity}]];
If[And[Length[vars]===0,Not[hasz]],Return[{True,Plus@@xx}]];
If[Plus@@xx===0,Return[{True,0}]];


        temp= Join[{Plus@@(xx)},
            If[hasz,
                Append[{##, IntegrationCut,1} & /@ vars,{im,0+ZIntegrationCut,1-ZIntegrationCut}],
(*                Append[{##, IntegrationCut,1} & /@ vars,{im,-1+ZIntegrationCut,1-ZIntegrationCut}],*)
(*                Append[{##, IntegrationCut,1} & /@ vars,{im,-Pi/2+ZIntegrationCut,Pi/2-ZIntegrationCut}],*)
(*                Append[{##, IntegrationCut,1} & /@ vars,{im,-1000,1000}],*)
                {##, IntegrationCut,1} & /@ vars
            ]
            ];
        If[hasz,
            temp[[1]]=temp[[1]]/(2Pi);

(* log *)
            temp=temp/.{z->-0.5 - I Log[im / (1-im)]};
            temp[[1]]=temp[[1]]/(im (1-im))


(* squares *)
(*            temp=temp/.{z->-0.5+(I im / (1-im^2))};
            temp[[1]]=temp[[1]]*(1+im^2)/((-1+im^2)^2)*)

(* tangent *)
(*            temp=temp/.{z->-0.5+I Tan[im]};
            temp[[1]]=temp[[1]]/(Cos[im]^2)*)

(* direct *)
(*            temp=temp/.{z->-0.5+I im};*)
        ];
        If[tryexact,
            res=Quiet[TimeConstrained[Integrate@@temp,ExactIntegrationTimeout,$Timeout]];
            If[And[Head[res]=!=Integrate,res=!=$Timeout],
                Return[{True,res}];
            ]
        ];
        temp=Join[temp,
            {Method -> AdaptiveQuasiMonteCarlo,
             Compiled->True,
             MinRecursion -> 100,
             MaxRecursion -> 10000,
             MaxPoints -> 100000,
             PrecisionGoal -> 4}];
        res=Apply[NIntegrate, temp, {0}];
(*        Print[{hasz,Length[temp[[1]]]},Drop[temp,1]];
        Print[res];
        *)

        {False,res}
]



ZNegExpand[xx_]:=Module[{i,j,deg,ders,der,result},
Check[
    i=1;
    If[ExpandMode,
        For[j=1,j<=Length[xx[[1]]],j++,
            If[And[xx[[1]][[j]]===1,
		Coefficient[xx[[2]][[j]],z]<0,
		Coefficient[xx[[2]][[j]],z2]<=0,
		(xx[[2]][[j]]/.{ep->0,z->0,z2->0})<=0,-xx[[4]][[3]]+(xx[[2]][[j]]/.{ep->0,z->0,z2->0} )<0],
                i=j
            ]
        ]
    ];
    For[Null,i<=Length[xx[[1]]],i++,
        If[Or[xx[[1]][[i]]=!=1,Coefficient[xx[[2]][[i]],z]===0],Continue[]];
        If[Or[xx[[1]][[i]]=!=1,And[ExpandMode,Coefficient[xx[[2]][[i]],z]>0]],Continue[]];
        If[ExpandMode,
            deg=-xx[[4]][[3]]+(xx[[2]][[i]]/.{ep->0,z->0,z2->0} ),
            deg=1+(xx[[2]][[i]]//.Flatten[{ep->0,ZREPLACEMENT}]);
        ];
        If[deg<0,
            Return[Flatten[ZNegExpand/@Reap[
                If[xx[[2]][[i]]-deg==0,Print["Impossible degree (z)"];Abort[]];
                ders={};(*vvv={};*)
                For[j=0,j+deg<0,j++,
                    If[j===0,der=xx[[3]],der=D[der,x[i]]/j];
                    (*AppendTo[vvv,(xx[[2]][[i]]+j+1)];*)
                    AppendTo[ders,der/.x[i]->0];
	      If[Expand[ders[[j+1]]]=!=0,
                    Sow[{ReplacePart[xx[[1]],0,i],
                         ReplacePart[xx[[2]],0,i],
                         ders[[j+1]]/(xx[[2]][[i]]+j+1),
                          {xx[[4]][[1]],Append[xx[[4]][[2]],(xx[[2]][[i]]+j+1)],xx[[4]][[3]]-If[Length[Position[xx[[2]],xx[[2]][[i]]]]>0,0,j]},
                          xx[[5]]
                    }];
	      ];
                ];
	      If[Expand[xx[[3]]-(ders.Table[x[i]^(j-1),{j,Length[ders]}])]=!=0,
                Sow[{ReplacePart[xx[[1]],2,i],
                    If[TrueQ[ExpandDirect],ReplacePart[xx[[2]],delta[j]+xx[[2]][[i]],i],xx[[2]]],
                    If[TrueQ[ExpandDirect],xx[[3]],xx[[3]]-(ders.Table[x[i]^(j-1),{j,Length[ders]}])],
                    xx[[4]],
                    xx[[5]]
                }];
	      ]
            ][[2]][[1]],1]];
        ];



    ];
,Print[xx];Abort[]];
{xx}

];


EpNegExpand[xx2_]:=Module[{i,j,deg,ders,der,HasDelta,temp,xx},
Check[
    xx=xx2;
    For[i=1,i<=Length[xx[[1]]],i++,
        If[xx[[2]][[i]]===0,Continue[]];
	HasDelta=MemberQ[Variables[xx[[2]][[i]]],delta[_]];
	If[HasDelta,DeltaPower=Cases[Variables[xx[[2]][[i]]],delta[_]][[1]][[1]]];
	temp=xx[[2]];
	temp[[i]]=temp[[i]]/.delta[_]->0;
	xx[[2]]=temp;
        deg=(xx[[2]][[i]]/.Flatten[{ep->0,ZREPLACEMENT}]);
        If[And[xx[[1]][[i]]===1,deg<=0],
            Return[Flatten[EpNegExpand/@Reap[
            If[And[xx[[3]][[i]]===0,ResolutionMode==="IBP0"],
		If[HasDelta,Print["Resolution mode does not work with delta"];Abort[]];
                der=xx[[4]];
                For[j=0,j+deg<0,j++,
                    der=der/(j+1+xx[[2]][[i]]);
                    Sow[{ReplacePart[xx[[1]],0,i],
                         ReplacePart[xx[[2]],0,i],
                         ReplacePart[xx[[3]],0,i],
                         der/.(x[i]->0),
                         xx[[5]]
                    }];
                    der=-D[der,x[i]];
                    Sow[{ReplacePart[xx[[1]],2,i],
                         ReplacePart[xx[[2]],0,i],
                         ReplacePart[xx[[3]],0,i],
                         -der,
                         xx[[5]]
                    }];
                ];

                Sow[{ReplacePart[xx[[1]],2,i],
                    ReplacePart[xx[[2]],xx[[2]][[i]]+j,i],
                    xx[[3]],
                    der,
                    xx[[5]]
                }];
          ];
          If[And[xx[[3]][[i]]===0,ResolutionMode==="IBP1"],
		If[HasDelta,Print["Resolution mode does not work with delta"];Abort[]];
                ders={};
                der=xx[[4]];
                For[j=0,j+deg<0,j++,
                    der=der/(j+1+xx[[2]][[i]]);
                    AppendTo[ders,der/.x[i]->1];
                    Sow[{ReplacePart[xx[[1]],0,i],
                         ReplacePart[xx[[2]],0,i],
                         ReplacePart[xx[[3]],0,i],
                         Last[ders],
                         xx[[5]]
                    }];
                    der=D[der,x[i]]/i;
                    der=-der*i;
                ];
                Sow[{ReplacePart[xx[[1]],2,i],
                    ReplacePart[xx[[2]],xx[[2]][[i]]+j,i],
                    xx[[3]],
                    der,
                    xx[[5]]
                }];
          ];
          If[Or[xx[[3]][[i]]=!=0,ResolutionMode==="Taylor"],
                           ders={};
                For[j=0,If[HasDelta,j<DeltaPower,j+deg<0],j++,
                    If[j===0,der=xx[[4]],der=D[der,x[i]]/j];
                    AppendTo[ders,der/.x[i]->0];
		  If[Not[HasDelta],
                    Sow[{ReplacePart[xx[[1]],0,i],
                         ReplacePart[xx[[2]],0,i],
                         ReplacePart[xx[[3]],0,i],
                    If[((xx[[2]][[i]]+j+1)/.ep->0)===0,
                        ders[[j+1]] * ((-1)^xx[[3]][[i]]) * (xx[[3]][[i]]!)/(
                                                ((Coefficient[xx[[2]][[i]]+j+1,ep])^(xx[[3]][[i]]+1))
                                                )
                    ,
                            ders[[j+1]] * ((-1)^xx[[3]][[i]]) * (xx[[3]][[i]]!)/(
                                                ((xx[[2]][[i]]+j+1)^(xx[[3]][[i]]+1))
                                                )
                    ],
                    If[((xx[[2]][[i]]+j+1)/.ep->0)===0,
                            xx[[5]]-xx[[3]][[i]]-1,
                         xx[[5]]
                    ]

                    }];
		  ];
                ];
                Sow[{ReplacePart[xx[[1]],2,i],
                    xx[[2]],
                    xx[[3]],
                    xx[[4]]-(ders.Table[x[i]^(j-1),{j,Length[ders]}]),
                    xx[[5]]
                }];

          ];



            ][[2]][[1]],1]];
        ];



    ];

,Print[xx2];Abort[]];    

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
  temp = Sort[DeleteCases[DeleteCases[Variables[Expand[xx]],ExV],Log[ExV]]];
  If[MemberQ[temp,Indeterminate],Return[Indeterminate]];
  (xx /. (Rule[##, 0] & /@
       temp)) + ((Plus @@ (Abs[(Coefficient[xx, ##]/.Log[ExV]->0 )/.ExV->0]^2 & /@ temp))^(1/2))*
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

AdvancedVarDiffReplace[xx_, yy_] := Module[{temp, table,rule},
Check[
  temp = If[Head[xx]===List,xx[[3]][[2]],xx];
  temp = NegTerms[temp];
  temp = Variables[temp];
  temp = If[Head[xx]===List,xx[[3]][[2]],xx] /. {x[i_] :> If[MemberQ[temp, x[i]], x[i], 0]};
  temp = {Coefficient[Coefficient[temp, x[yy[[1]]]], x[yy[[2]]]],
    Coefficient[temp, x[yy[[1]]], 2] /. {x[yy[[1]]] -> 0, 
      x[yy[[2]]] -> 0},
    Coefficient[temp, x[yy[[2]]], 2] /. {x[yy[[1]]] -> 0, 
      x[yy[[2]]] -> 0}
    };
  temp=Expand/@temp;
  If[And[temp[[3]]=!=0,temp[[2]]=!=0],
  If[NumberQ[Together[temp[[1]]^2/(temp[[2]]*temp[[3]])]],
   temp = Sqrt[Together[temp[[3]]/temp[[2]]]];
   If[Or[Head[temp] === Rational, Head[temp] === Integer],
    temp = {{1, temp/2}, {0, 1/2}};
    rule=Apply[Rule, 
      Transpose[{(x[##] & /@ yy), temp.(y[##] & /@ yy)}], {1}];
    temp = xx /. rule;
    Return[Expand[temp /. y -> x]];
    ,
    RawPrintLn["STRANGE SQUARE ROOT!"];
    ]
   ]];
  Return[VarDiffReplace[xx, yy]];
,
Print[xx];
Print[yy];
Abort[];
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



VariableSeries[xx_,var_]:=Module[{temp,shift,s,result,i},

Check[
    temp=xx/.(Log[var]->CCCC);
    temp=temp/.(Log[aaa_]->Normal[Series[Log[aaa],{var,0,VarExpansionDegree}]]);
    temp=Expand[temp];
    shift=0;
    While[TrueQ[MemberQ[{ComplexInfinity,Indeterminate},Quiet[(temp) /. var -> 0]]],
        temp=Expand[temp*var,var];
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



ClearDatabase[index_,reopen_]:=Module[{name},
    name=RemoveFIESTAName[DataPath<>ToString[index]<>"/"];
    Quiet[
        QClose[name];
        Pause[1];
        QRemoveDatabase[name];
    ];
    If[reopen,QOpen[name]];
]

RemoveFIESTAName[xx_]:=StringReplace[xx,{"FIESTA`"->""}];


ReadResult[Final_]:=Module[{j,i,result,request},
    Label[start];
    If[UsingC,
        i=0;
        While[i==0,
            Pause[0.001];
            i=ReadyLink[Final]
        ];
        result=MyLinkRead[i];
        If[Head[result]===EvaluatePacket,	    
            result=Evaluate@@result;
	    request=StringDrop[result,2];
	    request=StringDrop[request,-6];
	    If[Head[REQUESTS[request]]===REQUESTS,
  	      Print["Requested by CIntegrate: ",request];
	      REQUESTS[request]=True;
	    ];
	    result=Evaluate[ToExpression[result]];
            LinkWrite[MyLinks[[i]],ToString[result]];
            Goto[start]
        ];
    ,
        temp={##,MyLinkBuffer[##]}&/@Select[Range[Length[MyLinks]],(MyLinkTask[##]>0)&];
        result=WaitNext[(##[[2]])&/@temp];
        i=temp[[Position[(##[[2]])&/@temp,result[[2]]][[1]][[1]]]][[1]];
        result=result[[1]];
        If[result[[3]]>0,CurrentIntegralExact=False];
    ];
    If[result===False,
        RawPrintLn["Internal error"];
        RawPrintLn[MyLinkBuffer[i]];
        Abort[];
    ];
    If[Or[result===$Failed,result===$Timeout,And[2*Abs[result[[2]]]>Abs[result[[1]]],Abs[result[[2]]]>1],result[[1]]===Indeterminate],
	If[IntegrationDebug,
	  Print[MyLinkBuffer[i]];
	  Print[result];
	  Abort[];
	];
        AppendTo[BadParts,MyLinkTask[i]];
        MyLinkTask[i]=0;
        MyLinkBuffer[i]={};
        Return[]
    ];
    If[And[result[[1]]===Indeterminate,AbortOnNaN],
        RawPrintLn["INDETERMINATE"];
        Print[MyLinkBuffer[i]];
        Abort[];
    ];
    MyLinkTask[i]=0;
    MyLinkBuffer[i]={};
    If[And[NoNan,result==={Indeterminate,Indeterminate}],
        RawPrintLn[""];
        RawPrintLn["Integration returned INDETERMINATE as an answer."];
        RawPrintLn["Run WhyIndeterminate[] for possible reasons;"];
        NoNan=False;
    ];
    result=result[[1]]+result[[2]]*ToExpression["pm"<>ToString[PMCounter++]];
    MaxEvaluationTime=Max[MaxEvaluationTime,AbsoluteTime[]-MyLinkStart[i]];
    (*RawPrintLn[result];*)
    Dots[];
    result3=result3+result;
]

FindFreeLink[]:=Module[{i},
    For[i=1,i<=Length[MyLinks],i++,
        If[And[MyLinkType[i]=!="DeadLink",MyLinkBannedUntill[i]>0,MyLinkBannedUntill[i]<AbsoluteTime[]],
            MyLinkBannedUntill[i]=0;
            While[MyLinkReadyQ[i],MyLinkRead[i]];
            RawPrintLn["Unfreezing link ",MyLinkPath[i]];
        ];
	If[And[MyLinkType[i]==="DeadLink",MyLinkTask[i]==0,Restarts<=RestartLimit],	  
	      RestartLimit--;
	      If[MemberQ[Links[],MyLinks[[i]]],LinkClose[MyLinks[[i]]]];
	      MyLinkType[i]="LocalCIntegrate";
	      RawPrintLn["Restarting link"];
	      MyLinkInitialize[i];
	      SetOneLinkOptions[i];
	]
    ];
    i=FreeLink[];
    If[i==0,
        ReadResult[False];
    ];
    i=FreeLink[];
    i
]






MyParallelize[x_]:=If[NumberOfSubkernels>0,Parallelize[x(*,Method->"FinestGrained"*)],x];
SetAttributes[MyParallelize,HoldFirst];

SectorDecomposition[zs_,U_,intvars_,n_,j_]:=Module[{res,pol,active,SDDone,i,temp,degs},
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
                active=(Union@@(Sort[Cases[Variables[##],x[_]]]&/@pol));
                res=FindSD[MyDegrees[##,active]&/@pol];
                active=##[[1]]&/@active;
                res=BlockMatrix[res,active,n];
                Goto[SDDone];
            ];
	    If[Or[STRATEGY===STRATEGY_KU0,STRATEGY===STRATEGY_KU,STRATEGY===STRATEGY_KU2],
		pol=Expand[Times@@U]/.rules;
		active=##[[1]]&/@Sort[Cases[Variables[pol],x[_]]];
		degs=MyDegrees[pol];
		res=Reap[
		  For[i=1,i<=Length[degs],i++,
		    temp=(##-degs[[i]])&/@degs;
		    temp=Join[temp,IdentityMatrix[Length[active]]];
		    temp=SimplexCones[temp];
		    Sow/@temp;
		  ];
		][[2]][[1]];
		res=Map[(##/(GCD@@##))&,res,{2}];
		res=(If[Det[##]>0,##,Join[{##[[2]],##[[1]]},Drop[##,2]]])&/@res;
                res=BlockMatrix[res,active,n];
		Goto[SDDone];
	    ];
            (* other strategies *)
                pol=Expand[Times@@U]/.rules;
                active=##[[1]]&/@Sort[Cases[Variables[pol],x[_]]];
                res=FindSD[{MyDegrees[pol]}];
                res=BlockMatrix[res,active,n];
            Label[SDDone];

            RawPrintLn["Primary sector ",j," resulted in ",Length[res]," sectors."];
            res
]


RequiredResidues[gamm_,var_]:=Module[{gammas,temp,temp1,temp2,gam},
(*Print[gamm];Print[var];Abort[];*)
Check[
gam=gamm//. (Power[aa_, bb_] :> Timess @@ (Table[aa, {Abs[bb]}]));
    gammas = Join[Cases[gam,  Gamma[aa__], Infinity],Cases[gam,  Ggamma[aa__], Infinity]];
    gammas = Select[gammas, (Length[Cases[##[[1]], z, Infinity]] > 0) &];
(*Print[gammas];*)

    temp1=Flatten[(res = {};
            i=0;
	If[ExpandMode,
	    If[(Coefficient[##[[1]], z] /. {ep->0})==0,i=-Infinity];
	    If[(Coefficient[##[[1]], z2] /. {ep->0})>0,i=-Infinity];
	    If[(Coefficient[##[[1]], z] /. {ep->0})>0,
		If[Length[##]===2,
		    While[i>-##[[2]],
			AppendTo[res,{((i-(##[[1]]/.{z->0}))/Coefficient[##[[1]],z]),-Sign[Coefficient[##[[1]],z]]}];i--;
		    ]
		]
	    ,
		ExpandOrder2=ExpandOrder;
		If[(##[[1]] /.{z->0,z2->0,ep->0})>0,
		  ExpandOrder2=ExpandOrder-(##[[1]] /.{z->0,z2->0,ep->0})
		];
		While[If[Length[##]===1,i>=-ExpandOrder2,i>=-ExpandOrder2+##[[2]]],                               
		    AppendTo[res,{((i-(##[[1]]/.{z->0}))/Coefficient[##[[1]],z]),Sign[Coefficient[##[[1]],z]]}];i--;
		];
	    ]		
	,
	    If[(Coefficient[##[[1]], z] /. {ep->0})==0,i=-Infinity];
	    While[((i-(##[[1]]/.{z->0,ep->0}))/Coefficient[##[[1]],z]+(1/2))*Sign[Coefficient[##[[1]],z]]>=0,                
		AppendTo[res,{((i-(##[[1]]/.{z->0}))/Coefficient[##[[1]],z]),Sign[Coefficient[##[[1]],z]]}];i--;
	    ];
	];
	res
    ) & /@ gammas, 1];


    temp2=Flatten[(res = {};
        If[(Coefficient[##[[1]], z] /. {ep->0})==0,Continue[]];
        i=0;
(*        If[((i-(##[[1]]/.{z->0,ep->0}))/Coefficient[##[[1]],z]+(1/2))*Sign[Coefficient[##[[1]],z]]>=0,*)
            AppendTo[res,{((i-(##[[1]]/.{z->0}))/Coefficient[##[[1]],z]),Sign[Coefficient[##[[1]],z]]}];i--;
(*        ];*)
        res
    ) & /@ ({##}&/@var), 1];



    temp=Join[temp1,temp2];



temp=ExpandAll[temp];

If[Not[And@@(Reap[Sow[HHH[##[[2]]],##[[1]]]&/@temp,_,(Length[Union[#2]]===1)&][[2]])],
  Print["something wrong in RequiredResidues"];
  Print[gam];Print[var];Abort[];
];
    (*temp=Append[##[[1]], Length[##]] & /@ Reap[Sow[##, ##[[1]]] & /@ temp][[2]];*)
temp=##[[1]] & /@ Reap[Sow[##, ##[[1]]] & /@ temp][[2]];

temp=Append[##,tt=##[[1]];Length[Join[Select[gammas,And[IntegerQ[Expand[##[[1]]/.z->tt]],
			    Expand[##[[1]]/.z->tt]<=0
			      ]&],
			Select[var,(Expand[##/.z->tt]===0)&]
		    ]]
		  ]&/@temp;
,
  Print[gam];Print[var];Abort[];
];
    temp
]

(*
MyExpand[f_,var_]:=f//. RuleDelayed[
  var^(aa_: 1)*(bb_:1)*expr_Plus, bb*((var^aa*##) & /@ expr)]*)

MyExpand[f_,var_]:=f//. RuleDelayed[
  var^(aa_: 1)*expr_Plus, ((var^aa*##) & /@ expr)]


MyResidue[vardegrees_,xx_,x0_,maxorder_]:=Module[{temp,i,j},
(*TimeConstrained[*)
        
    Check[
	If[Head[vardegrees]===Pair,
	  terms={{vardegrees[[1]],vardegrees[[2]],xx}}
	,
	  terms={{vardegrees,Table[0,{Length[vardegrees]}],xx}}
	];
(*Print[Timing[*)
        terms=terms //. ((aa_:1) z :> Expand[aa x0] + aa dx);
        terms=terms/. MBshiftRules[dx];

        terms={##[[1]],##[[2]],Distribute[(dx^maxorder)*##[[3]]]}&/@terms;

        terms = terms /. MBexpansionRules[dx, -1+maxorder];
      
(*]];*)


        For[i=1,i<=-1+maxorder,i++,
            terms=Append[Table[{##[[1]],ReplacePart[##[[2]],#[[2]][[j]]+1,j],##[[3]]*D[##[[1]][[j]],dx]},{j,Length[##[[1]]]}],
                        {##[[1]],##[[2]],D[##[[3]],dx]}]&/@terms;
            terms=Flatten[terms,1];
            terms=DeleteCases[terms,{aa_,bb_,0}];
        ];

        terms=({##[[1]],##[[2]],Expand[##[[3]],dx]/(-1+maxorder)!}&/@terms) /. (dx -> 0);

	terms=DeleteCases[terms,{aa_,bb_,0}];
        temp=terms;
        ,
        Print[vardegrees];Print[xx];Print[x0];Print[maxorder];Abort[]
    ];
(*Print[temp];*)
(*
,300,    Print[timeout];Print[vardegrees];Print[xx];Print[x0];Print[maxorder];Abort[]
];
*)
    Return[temp];
]


DefinedFor[x_] := (ReleaseHold[Apply[List, ##[[1]], 1]]) & /@ DownValues[x]
QPutWrapper[dbase_,key_,value_]:=If[UsingQLink,QPut[RemoveFIESTAName[DataPath<>ToString[dbase]<>"/"],key,value],Evaluate[ToExpression["db"<>ToString[dbase]]][key]=value]
QGetWrapper[dbase_,key_]:=If[UsingQLink,QGet[RemoveFIESTAName[DataPath<>ToString[dbase]<>"/"],key],Evaluate[ToExpression["db"<>ToString[dbase]]][key]]
QSizeWrapper[dbase_]:=If[UsingQLink,QSize[RemoveFIESTAName[DataPath<>ToString[dbase]<>"/"]],1000]
QListWrapper[dbase_]:=If[UsingQLink,QList[RemoveFIESTAName[DataPath<>ToString[dbase]<>"/"]],(##[[1]])&/@DefinedFor[Evaluate[ToExpression["db"<>ToString[dbase]]]]]
ClearDatabaseWrapper[dbase_]:=ClearDatabaseWrapper[dbase,True];
ClearDatabaseWrapper[dbase_,reopen_]:=If[UsingQLink,ClearDatabase[dbase,reopen],Clear[Evaluate[ToExpression["db"<>ToString[dbase]]]]]


    PerformStage[s_,termss_,d_,prefix_,flatten_]:=Module[{temp,result,ev,terms},tasks={};
            RawPrintLn[MyTimingForm[AbsoluteTiming[
                terms=termss;
                If[terms===0,(*If[UsingQLink,terms=0;terms=Length[DefinedFor[ToExpression["db"<>ToString[s]]]]]*)
                    terms=Length[QListWrapper[s]]
                ];
                ClearDatabaseWrapper[d];
                PrepareDots[terms];
                For[runS=1;runD=1,runS<=terms,runS++,
                    Dots[];CC[100];
                    temp=QGetWrapper[s,prefix<>ToString[runS]];
                    If[NumberOfSubkernels===0,
                        result=EvalFunction[temp];
                        If[flatten,result=result,result={result}];
                        WriteFunction[d,##]&/@result;
                    ,
                        If[Length[tasks]===NumberOfSubkernels,
                            finished={};
                            Parallel`Developer`QueueRun[];
                            finished=Select[tasks,(Head[##[[4]]]===Parallel`Developer`finished)&];
                            working=Select[tasks,(Head[##[[4]]]=!=Parallel`Developer`finished)&];
                            If[Length[finished]===0,
                                {result,ev,tasks}=WaitNext[tasks];
                                result={result};
                            ,
                                result=WaitAll[finished];
                                tasks=working;
                            ];
                            AppendTo[tasks,ParallelSubmit[{temp},EvalFunction[temp]]];
                            Parallel`Developer`QueueRun[];
                            If[flatten,result=Flatten[result,1]];
                            WriteFunction[d,##]&/@result;
                        ,
                            AppendTo[tasks,ParallelSubmit[{temp},EvalFunction[temp]]];
                            Parallel`Developer`QueueRun[];
                        ];
                    ]
                ];
                If[And[NumberOfSubkernels>0,Length[tasks]>0],
                    Parallel`Developer`QueueRun[];
                    result=WaitAll[tasks];
                    If[flatten,result=Flatten[result,1]];
                    WriteFunction[d,##]&/@result;
                ];
            ][[1]]]," seconds; ",runD-1," terms."];
            If[MemoryDebug,RawPrintLn[MyMemoryInUse[]]];
            Return[runD-1];
    ]



PrintAllEntries[n_]:=Print[MyToExpression[FIESTA`QGetWrapper[n,##]]&/@Sort[FIESTA`QListWrapper[n],(MyToExpression[#1]<MyToExpression[#2])&]];


FindPoint[gamm_, var_] := 
 Module[{gammas, temp, temp1, temp2, gam}, 
  Check[gam = 
    gamm //. (Power[aa_, bb_] :> Timess @@ (Table[aa, {Abs[bb]}]));
   gammas = 
    Join[Cases[gam, Gamma[aa__], Infinity], 
     Cases[gam, FIESTA`Ggamma[aa__], Infinity]];
   gammas = 
    Union[Select[gammas, (Length[Cases[##[[1]], z, Infinity]] > 0) &],
      Select[gammas, (Length[Cases[##[[1]], z2, Infinity]] > 0) &]];
   temp1 = (If[
        Or[Length[##] === 1, ##[[2]] === 0], {##[[1]] > 
          1/111}, {##[[1]] > -##[[2]]+1/111, ##[[1]] < -##[[2]] + 1-1/111}]) & /@ 
     gammas;
   temp1 = Flatten[temp1, 1];
   temp2 = (## > 1/111) & /@ var;
(*Print[temp1];
Print[temp2];*)
   temp = Union[temp1, temp2] /. ep -> 1/1111;
   temp = FindInstance[temp, {z, z2}];
   If[Length[temp] == 0, Print["FindPoint error"]; Print[gamm]; 
    Print[var]; Abort[]];
  (z2 /. temp[[1]])
   , Print["FindPoint error"]; Print[gamm]; Print[var]; Abort[];]
  
  ]


SDIntegrate[intvars_,ZZ_,order_Integer,deltas_]:=
Module[{Z,U,F,SD,f,forsd,ii,i,vars,n,m,j,l,rule,Z3,Z2,U2,F2,coeff,k,Jac,res,result,md,result2,pol,active,zsets,SDCoeffs,timecounter,HasExtra},
    Restarts=0;
    Clear[TermNumber];
(*    temp=DeleteCases[DeleteCases[DeleteCases[Variables[ZZ],x[_]],ep],ExpandVariable];
    Print[temp];*)
   (* Print[ExpandVariable];
    HasExtra=(Length[temp]>0);
    If[HasExtra,RawPrintLn["WARNING: VARIABLES WITHOUT VALUES, NO INTEGRATION WILL BE PERFORMED"]];*)


    MixSectorCounter=-1;
    NoNan=True;
    HadIntegrationTimeout=False;
    timecounter=AbsoluteTime[];
    n=Length[intvars];
        If[NumberOfSubkernels>0,
            DistributeDefinitions[DataPath,QHullPath,ResolutionMode,UsingQLink,HASZ,ExpandOrder,ExV,ZREPLACEMENT,n,x,y,Power2,ExpandDirect,ExpandMode];
        ];
    Terms=0;
    (SDIntegral[intvars,ZZ,deltas,##[[1]],{##[[2]],##[[3]]}]=(##[[4]]+##[[5]]*ToExpression["pm"<>ToString[PMCounter++]]))
                &/@KnownIntegrals;
    (SDExact[intvars,ZZ,deltas,##[[1]],{##[[2]],##[[3]]}]=False)
                &/@KnownIntegrals;
    If[StartingStage>1,
        If[UsingQLink,Get[DataPath<>"3"]];
        Goto[StartingStage]
    ];

    ClearDatabaseWrapper[1];
    run1=1;
For[l=1,l<=Length[ZZ],l++,
 {Z,U}=ZZ[[l]];
 If[Times@@U===0,Continue[]];

    vars=Apply[x, Position[intvars, 1], {1}];

If[Length[vars]>0,

    zsets=Tuples[deltas];
    If[Head[PrimarySectorCoefficients]===Symbol,SDCoeffs=Table[1,{Length[zsets]}],SDCoeffs=PrimarySectorCoefficients];
    zsets=Delete[zsets,Position[SDCoeffs,0]];
    SDCoeffs=Delete[SDCoeffs,Position[SDCoeffs,0]];
    SD=Transpose[{zsets,Table[U,{Length[zsets]}],Table[intvars,{Length[zsets]}],Table[n,{Length[zsets]}],Range[Length[zsets]]}]/.Log[_]->1;
(*Print[SD];*)
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

(*    RawPrintLn[n];*)

        RawPrint["Preparing database: "];
        RawPrintLn[MyTimingForm[AbsoluteTiming[
            For[i=1,i<=Length[SD],i++,
                For[j=1,j<=Length[SD[[i]]],j++,
                    QPutWrapper[1,ToString[run1,InputForm],MyToString[{SD[[i]][[j]],SDCoeffs[[i]],zsets[[i]],Z,U,intvars,n}]];run1++;
                ];
            ];
        ][[1]]]," seconds. "];
,
    QPutWrapper[1,ToString[run1,InputForm],MyToString[{{},{},{},Times@@U,0}]];run1++;
];

];Terms=run1-1;


If[AllEntriesDebug,PrintAllEntries[1]];

    Label[2];
    If[UsingQLink,
        Put[Null,DataPath<>"3"];
        Save[DataPath<>"3", SCounter]
    ];

    If[AbortStage===2,Abort[]];


     (*  PrintAllEntries[1];Abort[];*)
    RawPrint["Variable substitution"];
        SCounter=1;



        EvalFunction[value_]:=Module[{temp,ThisSD,rules,reps,rule,Jac,U2,Z2,Z3,U3,NewDegrees},(ClearSystemCache[];
                    ThisSD=MyToExpression[value];
		    If[Length[ThisSD]===5,Return[{ThisSD}]];
											       n=ThisSD[[7]];
                    If[ThisSD[[2]]===0,{},
                        rules=Rule[x[##],1]&/@ThisSD[[3]];
                        reps=Inner[Power[#2, #1] &, Transpose[ThisSD[[1]]], Array[y,n],Times];
                        rule=Apply[Rule,Transpose[{Array[x,n],reps}],{1}];
                        Jac=Det[Outer[D,reps,Array[y,n]]];
                        Z2=(ThisSD[[4]]/.rules)//.rule;
                        U2=(ThisSD[[5]]/.rules)//.rule;
                        (
                            NewDegrees=Table[0,{n}];
			    NewLogDegrees=Table[0,{n}];
                            Z3=##;
                            For[m=1,m<=Length[U2],m++,
                                temp=FactorMonom[U2[[m]] /. {0.->0}];
                                U3[m]=temp[[2]];
                                NewDegrees=NewDegrees+Z3[[2]][[m]]*Table[Exponent[temp[[1]],y[i]],{i,n}];
				temp=FactorLogMonom[U3[m] /. {0.->0}];
                                U3[m]=temp[[2]];
                                NewLogDegrees=NewLogDegrees+Z3[[2]][[m]]*Table[Exponent[temp[[1]],Log[y[i]]],{i,n}];
                            ];
                            Z3[[1]]=Expand[Z3[[1]]];
                            If[Head[Z3[[1]]]===Plus,Z3[[1]]=List@@Z3[[1]],Z3[[1]]={Z3[[1]]}];
                            NewDegrees=NewDegrees+Table[Exponent[Jac,y[i]],{i,n}];
                            f={ReplacePart[ThisSD[[6]],0,List/@ThisSD[[3]]],Expand[NewDegrees],Z3[[1]],ThisSD[[2]]*(Jac/.y[aaa_]->1)*Times@@Table[((U3[m])^(Z3[[2]][[m]])),{m,Length[U2]}],0,NewLogDegrees};
                            f=f//.y->x;
                            f
                        )&/@Z2
                    ])];
        DistributeDefinitions[EvalFunction];
        WriteFunction[d_,value_]:=(MixSectorCounter++;
                 If[MixSectorCounter>MixSectors,SCounter++;MixSectorCounter=0];
                 QPutWrapper[d,ToString[runD,InputForm],MyToString[{SCounter,value}]];runD++);


        Terms=PerformStage[1,Terms,2,"",True];
 (*  PrintAllEntries[1];Abort[];*)

        (* {vars, indices, List of coeffs, expression, ep degree, log degree} *)

    Label[3];
    If[UsingQLink,
        Put[Null,DataPath<>"3"];
        Save[DataPath<>"3", SCounter]
    ];
    If[AbortStage===3,Abort[]];

If[AllEntriesDebug,PrintAllEntries[2]];

    RawPrint["Decomposing ep-independent term"];

        EvalFunction[value_]:=(
                    ClearSystemCache[];
                    temp=MyToExpression[value];
                    num=temp[[1]];
(*If[num===123,Print[temp[[2]]]];*)
                    temp=ExpandZ[temp[[2]]];
                    MyToString[{num,##}]&/@temp
                );
                DistributeDefinitions[EvalFunction];
        WriteFunction[d_,value_]:=(QPutWrapper[d,ToString[runD,InputForm],value];runD++);


        Terms=PerformStage[2,Terms,1,"",True];
        (* {vars, indices, log degrees, expression, {gammas,{}}, ep degree}   *)



If[DMode,

      result=Flatten[MyToExpression[QGetWrapper[1,ToString[##]]][[2]][[2]]&/@Range[Terms],1];
      result=Union[result];
      result=result/.{ep->(d0-d)/2};
      result=Expand/@result;
      result=Select[result,(Variables[##]=!={})&];
      result={##,Sort[{##/.d->DMIN,##/.d->DMAX}]}&/@result;
      result={##[[1]],{Ceiling[##[[2]][[1]]],Min[Floor[##[[2]][[2]]],-1]}}&/@result;
      result=Select[result,(##[[2]][[1]]<=##[[2]][[2]])&];
      result = {##[[1]], Range @@ (##[[2]])} & /@ result;
      result={Table[##[[1]],{Length[##[[2]]]}],##[[2]]}&/@result;
      result=Transpose/@result;
      result=Flatten[result,1];
      result=((##[[2]]-(##[[1]]/.d->0))/Coefficient[##[[1]],d])&/@result;
      If[TrueQ[PolesMultiplicity],
	  result=Reap[Sow[##, ##] & /@ result, _, {#1, Length[#2]} &][[2]]
      ,
	  result=Union[result];
      ];
      Return[result];
    ];

If[AllEntriesDebug,PrintAllEntries[1]];
    If[AbortStage===4,Abort[]];



If[HASZ,

    Label[4];
    If[UsingQLink,
        Put[Null,DataPath<>"3"];
        Save[DataPath<>"3", SCounter]
    ];


    RawPrint["Pole resolution"];

        EvalFunction[value_]:=(
                    ClearSystemCache[];
                    temp=MyToExpression[value];
                    num=temp[[1]];
                    temp=temp[[2]];
		    templogs=temp[[7]];
                    temp={temp[[1]],temp[[2]],temp[[3]]*temp[[4]],temp[[5]],temp[[6]]};
                    temp=ZNegExpand[temp];
		    temp=Append[##,templogs]&/@temp;
                    MyToString[{num,##}]&/@temp
                );
                DistributeDefinitions[EvalFunction];
        WriteFunction[d_,value_]:=(QPutWrapper[d,ToString[runD,InputForm],value];runD++);

        Terms=PerformStage[1,Terms,2,"",True];

        (* {vars, indices, expression, {gammas,variables divided}, ep degree} *)


     Label[5];
    If[UsingQLink,
        Put[Null,DataPath<>"3"];
        Save[DataPath<>"3", SCounter]
    ];
        If[AbortStage===5,Abort[]];


If[AllEntriesDebug,PrintAllEntries[2]];


(*PrintAllEntries[2];Abort[];    *)
(*PrintAllEntries[2];*)
     RawPrint["Taking residues"];

(*NUM=1;*)
        EvalFunction[value_]:=(

                        ClearSystemCache[];
                        temp=MyToExpression[value];
                        num=temp[[1]];
(*num=NUM;NUM++;*)
                        temp=temp[[2]];
                        temp1=temp[[1]];
                        temp0=Part[temp,-2];
			temp[[4]]=temp[[4]]/.Gamma[arg_]->Ggamma[arg,0];
			ReqR=RequiredResidues@@(Take[temp[[4]],2]);
			temp2={};
			
(*pass=0;*)
(*We need to find a starting point*)
(*Print[Take[temp[[4]],2]];Abort[];*)

SmallPoint=1/57;
If[Power2,
  SmallPoint=FindPoint@@(Take[temp[[4]],2]);
];

(*Print[SmallPoint];*)
StrangeReplacement={z2->SmallPoint,ep->1/1111};
(*Print[ReqR];
Print[StrangeReplacement];*)
ReqR=Sort[ReqR, ((#1[[1]]/.StrangeReplacement)<(#2[[1]]/.StrangeReplacement))&];
(*Print[ReqR];*)
(*Abort[];*)

(*If[num==4,*)
(*ReqR=Select[ReqR,(##[[1]]==-2-4 ep-2 z2)&];*)
(*ReqR={}];*)  (*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*)
(*Print[ReqR]; *)



			While[Length[ReqR]>0,	
(*Print[pass];		    
Print[temp];
Print[ReqR];*)(*Print[ReqR];
Abort[];*)
			    res1=First[ReqR];
(*Print[uuuu];*)
			    uuuu={z->res1[[1]]};

(*Print["Part ",num,": taking residue in z1 = ",res1[[1]]," with coefficient ",res1[[2]]];*)
			    res1taken=MyResidue[Pair[temp[[2]],Last[temp]],res1[[2]]*temp[[3]],res1[[1]],res1[[3]]];


			    temp2=Join[temp2,
				Append[Append[##,
				      Flatten[Table[Select[(temp[[4]][[2]]/.uuuu),MemberQ[Variables[Expand[##]],z2]&],{res1[[3]]}],1]
				      ],
				      Table[(ExpandAll[(temp[[4]][[1]]/.Gamma->Ggamma)/.uuuu]/.(Ggamma[_Integer,_]->1)),{res1[[3]]}]
				]&/@res1taken
			    ];
(*Print["!"];
Print[temp[[4]]];*)
			    temp[[4]]={(temp[[4]][[1]]/.(Ggamma[arg_,shift_]:>If[Expand[arg/.uuuu]===0,
						      Ggamma[arg,shift+1],Ggamma[arg,shift]
							]))*
				      Times@@(Ggamma[##,1]&/@(Select[temp[[4]][[2]],(Expand[##/.uuuu]===0)&])),
						      Select[temp[[4]][[2]],(Expand[##/.uuuu]=!=0)&],temp[[4]][[3]]};			    
(*Print["!"];
Print[temp[[4]]];*)

			    ReqR=RequiredResidues@@(Take[temp[[4]],2]);
ReqR=Sort[ReqR, ((#1[[1]]/.StrangeReplacement)<(#2[[1]]/.StrangeReplacement))&];
(*ReqR=Select[ReqR,(##[[1]]==-2-4 ep-2 z2)&]; *)(*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*)

(*pass++;
If[pass>=5,Abort[]];*)
			];
(*Print[temp];

Print[temp2];*)

(*Print[temp2];Abort[];*)

                        temp2=Select[temp2,(##[[3]]=!=0)&];
			temp2=If[TrueQ[Not[Coefficient[Last[##[[1]]],z2]===0]],
			      jac=Coefficient[Last[##[[1]]],z2];
				  ReplacePart[##,##[[3]]/Abs[jac],3]/.(z2->(z2-(Last[##[[1]]]/.z2->0))/jac)  
			      ,##]&/@temp2;

(*			temp2=If[TrueQ[Coefficient[Last[##[[1]]],z2]<0],##/.(z2->-z2),##]&/@temp2;*)
                        If[Not[ExpandMode],AppendTo[temp2,{temp[[2]],Table[0,{Length[temp[[2]]]}],temp[[3]],{}}]];
			(*temp2=temp2/.{Gamma[garg__]:>Gamma@@(Expand/@{garg})};*)
			temp2=temp2/. Gamma[garg_] :> Gamma[Expand[garg]];
			temp2=temp2/. Ggamma[garg1_,garg2_] :> Ggamma[Expand[garg1],Expand[garg2]];
			temp2=temp2/. Power[garg1_,garg2_] :> Power[Expand[garg1],Expand[garg2]];
                        MyToString[{num, {temp1,ExpandAll[##[[1]]],ExpandAll[##[[2]]],##[[3]],
			temp0,##[[4]],##[[5]]}}]& /@ temp2
                );
                DistributeDefinitions[EvalFunction];
        WriteFunction[d_,value_]:=(QPutWrapper[d,ToString[runD,InputForm],value];runD++);

        Terms=PerformStage[2,Terms,1,"",True];
  If[AllEntriesDebug,PrintAllEntries[1]];
  
(* {vars, indices, log indices, expression, ep degree} ??????????????????????????????*)

(*Print["aga"];*)

	  If[Power2,


		Power2=False;

		RawPrint["Pole resolution"];

		EvalFunction[value_]:=(
                    ClearSystemCache[];
                    temp=MyToExpression[value];
                    num=temp[[1]];
                    temp=temp[[2]];
		    templogs=temp[[3]];
                    temp={temp[[1]],temp[[2]],temp[[4]],{temp[[7]],temp[[6]],0},temp[[5]]}/.(z2->z);      
                    temp=ZNegExpand[temp];
		    temp=Append[##,templogs]&/@temp;
                    MyToString[{num,##}]&/@temp
                );
                DistributeDefinitions[EvalFunction];
		WriteFunction[d_,value_]:=(QPutWrapper[d,ToString[runD,InputForm],value];runD++);

		Terms=PerformStage[1,Terms,2,"",True];

		(* {vars, indices, expression, {gammas,variables divided}, ep degree} *)

		If[AllEntriesDebug,PrintAllEntries[2]];
		





		RawPrint["Taking residues"];

		EvalFunction[value_]:=(
                        ClearSystemCache[];
                        temp=MyToExpression[value];
                        num=temp[[1]];
                        temp=temp[[2]];
                        temp1=temp[[1]];
                        temp0=temp[[5]];

ReqR=(RequiredResidues@@(Take[temp[[4]],2]));

(*ReqR=Select[ReqR,(##[[1]]===-2-4ep)&];*)
(*Print[ReqR];*)


                        temp2=Flatten[(
  (*Print["Part ",num,": taking residue in z2 = ",##[[1]]," with coefficient ",##[[2]]];*)
  MyResidue[Pair[temp[[2]],temp[[6]]],##[[2]]*temp[[3]],##[[1]],##[[3]]])
                                                        &/@ReqR,1];

                        If[Not[ExpandMode],AppendTo[temp2,{temp[[2]],Table[0,{Length[temp[[2]]]}],temp[[3]]}]];
                        MyToString[{num, Append[Prepend[##,temp1],temp0]}]& /@ temp2
                );
                DistributeDefinitions[EvalFunction];
		WriteFunction[d_,value_]:=(QPutWrapper[d,ToString[runD,InputForm],value];runD++);

		Terms=PerformStage[2,Terms,1,"",True];



		If[AllEntriesDebug,PrintAllEntries[1]];
		(*Abort[];*)



	  ];




];  (*HASZ*)



        (* {vars, indices, log indices, expression, ep degree} *)


    Label[6];
    If[UsingQLink,
        Put[Null,DataPath<>"3"];
        Save[DataPath<>"3", SCounter]
    ];
    If[AbortStage===6,Abort[]];

If[ExpandDirect,

(*PrintAllEntries[1];*)
      powers=g[(Last[MyToExpression[FIESTA`QGetWrapper[1,##]][[2]][[2]]]/.{ep->0,RegVar->0}),
	      (Last[MyToExpression[FIESTA`QGetWrapper[1,##]][[2]][[3]]]/.{ep->0,RegVar->0})]&/@FIESTA`QListWrapper[1];
      powers=Union[powers];
      powers=powers/.g->List;
      powers=Sort[powers,(If[#1[[1]]<#2[[1]],True,If[#1[[1]]>#2[[1]],False,#1[[2]]<#2[[2]]]])&];

      While[True,
	  If[Length[powers]===0,Print[0];Abort[]];
	  MinYPower=powers[[1]][[1]];
	  MinYLPower=powers[[1]][[2]];

	  result=Select[MyToExpression[FIESTA`QGetWrapper[1,##]]&/@FIESTA`QListWrapper[1],And[Last[##[[2]][[2]]/.{ep->0,RegVar->0}]===MinYPower,Last[##[[2]][[3]]/.{ep->0,RegVar->0}]===MinYLPower]&];
	  result={y^MinYPower Log[y]^MinYLPower,y^(Last[##[[2]][[2]]])*Log[y]^(Last[##[[2]][[3]]])*
		  Inner[(#1^#2)&,Array[x,Length[##[[2]][[2]]]-1],Drop[##[[2]][[2]],-1],Times]*
		  Inner[(Log[#1]^#2)&,Array[x,Length[##[[2]][[3]]]-1],Drop[##[[2]][[3]],-1],Times]*
		  ##[[2]][[4]]*If[Length[##[[2]][[5]]]===1,(ep^##[[2]][[5]]),(ep^##[[2]][[5]][[1]])*(Log[ep]^##[[2]][[5]][[2]])]
		   &/@result};

	  If[ExpandAll[Plus@@(result[[2]])]===0,
		(*Print[result[[1]]];
		Print[0];*)
		powers=Drop[powers,1];
	  ,

		result[[2]]=MySimplify/@(result[[2]]);
		
		If[Or@@((Length[##]>1)&/@(result[[2]])),Print["WARNING!!!!!!!!!!!!!!!!!"];  
		    Print[Select[result[[2]],(Length[##]>1)&]];
		];

		result[[2]]=Flatten[result[[2]],1];
(*		result[[2]]=Apply[Group,result[[2]],{1}];*)
		result[[2]]=Map[MyFactor,(result[[2]]),{1}];
		
		RawPrint["Degree is ",ToString[result[[1]],InputForm]];
		Return[Plus@@result[[2]]];		
	  ]
      ]
];

If[AllEntriesDebug,PrintAllEntries[1]];


   RawPrint["Pole resolution"];

        EvalFunction[value_]:=(
                    ClearSystemCache[];
                    temp=MyToExpression[value];
(*         If[Length[Cases[temp,ComplexInfinity,{0,Infinity}]]>0,Print[temp];Print[RequiredResidues@@(Take[temp[[4]],2])];Abort[]];
*)
                    num=temp[[1]];
                    temp=temp[[2]];
		   (* templogs=temp[[6]];*)
                    temp=EpNegExpand[temp];
		   (* temp=Append[##,templogs]&/@temp;*)
                    MyToString[{num,##}]&/@temp
                );

                DistributeDefinitions[EvalFunction];
        WriteFunction[d_,value_]:=(QPutWrapper[d,ToString[runD,InputForm],value];runD++);

        Terms=PerformStage[1,Terms,2,"",True];

        (* {vars, indices, log indices, expression, ep degree} *)

    Label[7];
    If[UsingQLink,
        Put[Null,DataPath<>"3"];
        Save[DataPath<>"3", SCounter]
    ];
    If[AbortStage===7,Abort[]];


If[AllEntriesDebug,PrintAllEntries[2]];

    RawPrint["Expression preparation"];

        EvalFunction[value_]:=(
                    ClearSystemCache[];
                    temp=MyToExpression[value];
                    MyToString[{temp[[1]],AdvancedFirstVars[ConstructTerm[temp[[2]]]]}]
                );
                DistributeDefinitions[EvalFunction];
        WriteFunction[d_,value_]:=(QPutWrapper[d,ToString[runD,InputForm],value];runD++);

        Terms=PerformStage[2,Terms,1,"",False];

        (* {variables, expression, ExpandVarDegrees, ep degree} *)
  (*
    RawPrint["Replacing variables"];

        EvalFunction[value_]:=(
                    ClearSystemCache[];
                    temp=MyToExpression[value];
                    MyToString[{temp[[1]],AdvancedFirstVars[temp[[2]]]}]
                );
                DistributeDefinitions[EvalFunction];
        WriteFunction[d_,value_]:=(QPutWrapper[d,ToString[runD,InputForm],value];runD++);

        Terms=PerformStage[2,Terms,1,"",False];
*)
        (* {variables, expression, ExpandVar degrees, ep degree} *)

If[AllEntriesDebug,PrintAllEntries[1]];

    Label[8];
    If[UsingQLink,
        Put[Null,DataPath<>"3"];
        Save[DataPath<>"3", SCounter]
    ];
    If[AbortStage===8,Abort[]];

        RawPrint["Epsilon expansion"];
        min=order+1;

        EvalFunction[value_]:=(
                    ClearSystemCache[];
                    temp=MyToExpression[value];
                    num=temp[[1]];
                    temp=temp[[2]];
                    temp1=temp[[1]];
                    temp2=EpPosExpand[temp[[2]],temp[[3]],temp[[4]],order];
                    temp={num,temp1,##}&/@temp2;
                    temp={##[[3]][[1]],MyToString[{##[[1]],{##[[2]],##[[3]][[2]][[1]]}}],##[[3]][[2]][[2]]}&/@temp
                );
                DistributeDefinitions[EvalFunction];
        WriteFunction[d_,value_]:=(
                    If[Head[TermNumber[{value[[1]],value[[3]]}]]===TermNumber,TermNumber[{value[[1]],value[[3]]}]=0];
                    TermNumber[{value[[1]],value[[3]]}]=TermNumber[{value[[1]],value[[3]]}]+1;
                    min=Min[min,value[[1]]];
                    QPutWrapper[d,ToString[value[[1]]]<>"-"<>ToString[value[[3]],InputForm]<>"-"<>ToString[TermNumber[{value[[1]],value[[3]]}]],value[[2]]];
                    runD++
                );

        Terms=PerformStage[1,Terms,2,"",True];

If[AllEntriesDebug,PrintAllEntries[2]];


        Label[9];
    If[UsingQLink,
        Put[Null,DataPath<>"3"];
        Save[DataPath<>"3", SCounter];
        Save[DataPath<>"3", min];
        Save[DataPath<>"3", TermNumber];
    ];
        If[AbortStage===9,Abort[]];



        If[MemoryDebug,RawPrintLn[MyMemoryInUse[]]];
        ForEvaluation=(##[[1]])&/@DefinedFor[TermNumber];
        ForEvaluation=Select[ForEvaluation,(##[[1]]<=order)&];
        If[min===order+1,Return[]];
        If[NumberOfSubkernels>0,
            DistributeDefinitions[UsingC];
        ];

        If[PrintInsteadOfIntegrating,
            For[EvC=1,EvC<=Length[ForEvaluation],EvC++,
                CurrentOrder=ForEvaluation[[EvC]][[1]];
                CurrentVarDegree=ForEvaluation[[EvC]][[2]];
                If[ExpandMode,
                    Print["Integral of ep order ",CurrentOrder," and ",ExpandVariable," order ",ToString[CurrentVarDegree,InputForm],": "],
                    Print["Integral of order ",CurrentOrder,": "]
                ];
                Print[Expand[Plus@@(
                QGetWrapper[2,ToString[CurrentOrder]<>"-"<>ToString[CurrentVarDegree,InputForm]<>"-"<>ToString[##]][[2]][[2]]&/@
                    Range[TermNumber[{CurrentOrder,CurrentVarDegree}]])]
                ];

            ];
            Abort[];
        ];

        For[EvC=1,EvC<=Length[ForEvaluation],EvC++,
            CurrentOrder=ForEvaluation[[EvC]][[1]];
            CurrentVarDegree=ForEvaluation[[EvC]][[2]];
            If[Not[Head[SDIntegral[intvars,ZZ,deltas,CurrentOrder,CurrentVarDegree]]===SDIntegral],
                Continue[]
            ];
            NumberOfVariables=0;
            (PSector[##]=0)&/@Range[SCounter];

            RawPrint["Preparing integration string"];


            EvalFunction[value_]:=(
                            ClearSystemCache[];
                            temp=MyToExpression[value];
                            num=temp[[1]];
                            temp=temp[[2]][[2]];
			    If[ZeroCheck[temp],
				{}
			    ,
			    vars=CountVars[Variables[temp]];		
			    temp2=MyWorstPower[temp,vars[[1]]];
			    temp2=Ceiling[Abs[##]]&/@temp2;
                            {{num,If[UsingC,"("<>MyString[temp]<>")*("<>MyString[Inner[Power,vars[[1]],temp2,Times]]<>")*("<>MyString[Inner[Power,vars[[1]],-temp2,Times]]<>");\n",MyToString[temp]],vars[[2]]}}
                           ]
                        );
            DistributeDefinitions[ZeroCheck,EvalFunction];
	(*   If[UsingC],
	    ParallelEvaluate[
	      Unprotect[Power, Log];
	      Power /: Format[Power[a_, b_], InputForm] := p[a, b];
	      Log /: Format[Log[a_], InputForm] := Global`l[a];
	      Protect[Power, Log];
	    ];
	    Unprotect[Power, Log];
	    Power /: Format[Power[a_, b_], InputForm] := p[a, b];
	    Log /: Format[Log[a_], InputForm] := Global`l[a];
	    Protect[Power, Log];
           ]; *)
	    WriteFunction[d_,value_]:=(NumberOfVariables=Max[NumberOfVariables,value[[3]]];
                            QPutWrapper[d,ToString[value[[1]]]<>"-"<>ToString[++PSector[value[[1]]]],value[[2]]];
                            runD++
                        );

            PerformStage[2,TermNumber[{CurrentOrder,CurrentVarDegree}],1,ToString[CurrentOrder]<>"-"<>ToString[CurrentVarDegree,InputForm]<>"-",True];
      (*     If[UsingC,
	    ParallelEvaluate[
	      Unprotect[Power, Log];
	      Power /: Format[Power[a_, b_], InputForm] =.;
	      Log /: Format[Log[a_], InputForm] =.;
	      Protect[Power, Log];
	    ];
            Unprotect[Power, Log];
	    Power /: Format[Power[a_, b_], InputForm] =.;
	    Log /: Format[Log[a_], InputForm] =.;
	    Protect[Power, Log];
	   ]; *)

            If[ExpandMode,
                RawPrint["Terms of ep order ",CurrentOrder," and ",ExpandVariable," order ",ToString[CurrentVarDegree,InputForm],": "],
                RawPrint["Terms of order ",CurrentOrder,": "]
            ];

            RawPrintLn[Plus@@(PSector[##]&/@Range[SCounter])," (",NumberOfVariables,"-fold integrals)."];
            RawPrint["Numerical integration: "];
            RawPrintLn[SCounter-Length[Position[PSector[##]&/@Range[SCounter],0]]," parts; ",Length[MyLinks]," links;"];
            PrepareDots[SCounter];
            RawPrint["Integrating"];
            CurrentIntegralExact=And[Not[UsingC],CurrentOrder<=ExactIntegrationOrder];
            CurrentIntegralExact=And[Not[UsingC],CurrentOrder<=ExactIntegrationOrder];
            RawPrintLn[MyTimingForm[AbsoluteTiming[
                If[NumberOfSubkernels>0,
                    DistributeDefinitions[ExactIntegrationTimeout,ZIntegrationCut,IntegrationCut];
                ];
                result3=0;
                Parts=Range[SCounter];
                MaxEvaluationTime=0;
                While[Length[Parts]>0,
                    BadParts={};
                    For[j=1,j<=Length[Parts],j++,
                        SubmitIntegration[Parts[[j]]];
                    ];
                    While[(Plus@@(MyLinkTask[##]&/@Range[Length[MyLinks]]))>0,ReadResult[True]];
                    If[Length[BadParts]>0,
                        RawPrintLn["BAD PARTS: ",BadParts];
                        Parts=BadParts;
                    ,
                        Parts={}
                    ]
                ];

                If[Not[CurrentIntegralExact],
                    result3=PMSimplify[result3];
                    result3=DoubleCutExtraDigits[result3];
                ,
                    result3=FullSimplify[result3];
                ];
            ][[1]]]," seconds; returned answer: ",ToString[result3,InputForm]];
            SDIntegral[intvars,ZZ,deltas,CurrentOrder,CurrentVarDegree]=result3;
            SDExact[intvars,ZZ,deltas,CurrentOrder,CurrentVarDegree]=CurrentIntegralExact;
            result3=GenerateAnswer[intvars,ZZ,CurrentOrder,CurrentOrder-(order-ORDER),deltas,OUTSIDE,DENOM];
            For[j=1,j<=Length[result3],j++,
                RawPrint["(",ToString[result3[[j]][[1]],InputForm],")*",ToString[result3[[j]][[2]],InputForm]];
                If[j<Length[result3],RawPrint["+"],RawPrintLn[""]];
            ];
        ]; (*for EvC*)

        RawPrintLn["Total time used: ",CutExtraDigits[AbsoluteTime[]-timecounter]," seconds."];
        Return[];
];



SubmitIntegration[j_]:=Module[{temp,iii},
        If[PSector[j]==0,Dots[];Continue[]];
        If[NumberOfVariables===0,
             Dots[];result3=result3+Plus@@(
                        If[UsingC,
                          (*  Print[StringDrop[QGetWrapper[1,ToString[j]<>"-"<>ToString[##]],-2]];*)
                          Chop[N[Round[ToExpression[StringDrop[QGetWrapper[1,ToString[j]<>"-"<>ToString[##]],-2]]/.{P->Pi,l->Log,p->Power},10^-6]],10^-6],
                            MyToExpression[QGetWrapper[1,ToString[j]<>"-"<>ToString[##]]]
                        ]&/@Range[PSector[j]]);
        ,
            ress=$Failed;
            While[ress===$Failed,
                 iii=FindFreeLink[];
                 If[UsingC,
                    MyLinkWrite[iii,ToString[Round[QSizeWrapper[1]]]<>";\n"];
                    MyLinkRead[iii];
                    MyLinkWrite[iii,ToString[NumberOfVariables]<>";\n"<>ToString[PSector[j]+1]<>";\n"<>"0."<>";\n"];
                    MyLinkRead[iii];
                 ];
                 (MyLinkWrite[iii,QGetWrapper[1,ToString[j]<>"-"<>ToString[##]]];MyLinkRead[iii];)&/@Range[PSector[j]];
                 ress=MyLinkSubmit[iii,(StringJoin@@(("f["<>ToString[##]<>"]+")&/@Range[2,PSector[j]+1]))<>"0;\n"];
            ];
            If[And[MyLinkType[iii]==="MainKernel"],
                Dots[]
            ,
                MyLinkStart[iii]=AbsoluteTime[];
                MyLinkTask[iii]=j;
            ];
        ];
]



IntegrateHere[expression_,tryexact_]:=Module[{ff,result,result2},
                ff =Flatten[{{True,Select[##, (Length[Cases[##, z, Infinity]] > 0) &]},
                           {False,Select[##, (Length[Cases[##, z, Infinity]] == 0) &]}}&/@expression,1];
                ff=DeleteCases[ff,{aa_,{}}];
                Localresult3=
                    (
                    result2=Reap[Block[{$MessagePrePrint = Sow,$Messages = {}},
                        Quiet[
                            result=DoIntegrate[##,tryexact];
                        ,{Power::infy, Power::indet, Infinity::indet,General::stop}
                    ];
                    ]][[2]];
                    If[result[[1]],
                        Localresult3={result[[2]],0,0}
                    ,
                        result=result[[2]];
                        result2=DeleteCases[##, HoldForm[$MessageList]] & /@result2;
			result2=DeleteCases[##, HoldForm[y]] & /@result2;
			result2=DeleteCases[##, HoldForm[y[_]]] & /@result2;
			result2=DeleteCases[result2,{}];
                        Good=False;
			If[And[Length[$MessageList]>=1,Last[$MessageList]===HoldForm[NIntegrate::"maxp"],Length[result2]===1,
                            ReleaseHold[result2[[1]][[2]]]===result],
                            Good=True;
                            Localresult3=Append[ReleaseHold/@Drop[result2[[1]],1],1];
                        ];
			If[Length[result2]===0,
                            Good=True;
                            Localresult3={result,0,1};
                        ];
                        If[Not[Good],
                            RawPrintLn[PMForm[result]];
                            RawPrintLn["Something went wrong"];
                            RawPrintLn[$MessageList];
                            RawPrintLn[result2];
                            Abort[];
                        ];
                    ];
                    Localresult3
                    )&/@ff;
                 Plus@@Localresult3
]


MyDelta[xx_,yy_]:=Table[If[iii===xx,1,0],{iii,yy}]

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
    vars=#[[1]]&/@Cases[Variables[temp],x[_]];
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

AdvancedKillNegTerms[ZZ_,level_]:=Module[{temp,vars,subsets,U,F,UF,Z,i,mins,rs,r1,r2},
(*Return[{ZZ}];*)    
    If[level===1,RawPrintLn["Advanced netative terms resolution (might take some time)"]];
    {Z,MM,UF}=ZZ;
    temp=NegTerms[UF[[2]]];
    If[temp==={},Return[{{ZZ},0}]];
  (*  RawPrintLn["Negative terms encountered: ",InputForm[temp]];*)
    vars=#[[1]]&/@Cases[Variables[temp],x[_]];
    subsets=Subsets[vars,{2}];
    mins=Table[Infinity,{Length[subsets]}];
    rs=Table[{},{Length[subsets]}];
    For[i=1,i<=Length[subsets],i++,
        If[And[Length[NegTerms[VarDiffReplace[UF[[2]],subsets[[i]]]]]<Length[temp],
               Length[NegTerms[VarDiffReplace[UF[[2]],{subsets[[i]][[2]],subsets[[i]][[1]]}]]]<Length[temp]]
                ,
	  (*Print[subsets[[i]]];*)
            r1=AdvancedKillNegTerms[VarDiffReplace[{{#[[1]]/2,#[[2]]}&/@Z,MM,UF},subsets[[i]]],level+1];
            r2=AdvancedKillNegTerms[VarDiffReplace[{{#[[1]]/2,#[[2]]}&/@Z,MM,UF},{subsets[[i]][[2]],subsets[[i]][[1]]}],level+1];
            rs[[i]]=Join[r1[[1]],r2[[1]]];
            mins[[i]]=r1[[2]]+r2[[2]];
        ];
    ];
    

    i=Position[mins,Min[mins],{1}][[1]][[1]];
    If[mins[[i]]===Infinity,
       (*If[Length[temp]>0,Print[ZZ];Print[temp];Abort[]];*)
      (*  Print["Recursion level ",level,": ",Length[temp]," negative terms remaining"];*)
        Return[{{ZZ},Length[temp]}],
(*        Print["Recursion level ",level,": total number of negative terms below =",mins[[i]]];*)
        If[level===1,RawPrintLn["Total number of negative terms remaining in subexpressions: ",mins[[i]]]];
        If[level===1,RawPrintLn["Total number of subexpressions: ",Length[rs[[i]]]]];
        Return[{rs[[i]],mins[[i]]}];  
    ];
]

AdvancedKillNegTerms2[ZZ_,level_]:=Module[{temp,vars,subsets,U,F,UF,Z,i,mins,rs,r1,r2,pairs,variants},
(*Return[{ZZ}];*)    
    If[level===1,RawPrintLn["Advanced netative terms resolution, version 2 (might take some time)"]];
    {Z,MM,UF}=ZZ;
    temp=NegTerms[UF[[2]]];
    If[temp==={},Return[{{{ZZ},0}}]];
  (*  RawPrintLn["Negative terms encountered: ",InputForm[temp]];*)
    vars=#[[1]]&/@Cases[Variables[temp],x[_]];
    subsets=Subsets[vars,{2}];
    variants={};
    (*mins=Table[Infinity,{Length[subsets]}];
    rs=Table[{},{Length[subsets]}];*)
    For[i=1,i<=Length[subsets],i++,
        If[And[Length[NegTerms[VarDiffReplace[UF[[2]],subsets[[i]]]]]<Length[temp],
               Length[NegTerms[VarDiffReplace[UF[[2]],{subsets[[i]][[2]],subsets[[i]][[1]]}]]]<Length[temp]]
                ,
(*Print[subsets[[i]]];
Print["in"];*)
            r1=AdvancedKillNegTerms2[VarDiffReplace[{{#[[1]]/2,#[[2]]}&/@Z,MM,UF},subsets[[i]]],level+1];
            r2=AdvancedKillNegTerms2[VarDiffReplace[{{#[[1]]/2,#[[2]]}&/@Z,MM,UF},{subsets[[i]][[2]],subsets[[i]][[1]]}],level+1];
(*If[level==1,AAAr1=r1;
AAAr2=r2;Abort[]];*)



(*Print["out"];*)
	    pairs=Tuples[{r1,r2}];
	    variants=Join[variants,{Join[##[[1]][[1]],##[[2]][[1]]],##[[1]][[2]]+##[[2]][[2]]}&/@pairs];
(*If[level==1,Print[Length[variants[[1]][[1]]]];Abort[];];*)

        ];
    ];
      If[Length[variants]==0,Return[{{{ZZ},Length[temp]}}]];
	If[level===1,
	  min=Min@@((##[[2]])&/@variants);
	  i=Position[((##[[2]])&/@variants),min][[1]][[1]];	
          RawPrintLn["Total number of negative terms remaining in subexpressions: ",variants[[i]][[2]]];
	  RawPrintLn["Total number of subexpressions: ",Length[variants[[i]][[1]]]];
	  Return[{variants[[i]]}];	
	,
	  Return[variants];
	];  
   
]


AdvancedKillNegTerms3[ZZ_,level_]:=Module[{temp,vars,subsets,U,F,UF,Z,i,mins,rs,r1,r2,pairs,variants},
(*Return[{ZZ}];*)    
    If[level===1,RawPrintLn["Advanced netative terms resolution, version 3 (might take some time)"]];
    {Z,MM,UF}=ZZ;
    temp=NegTerms[UF[[2]]];
    If[temp==={},Return[{{{ZZ},0}}]];
  (*  RawPrintLn["Negative terms encountered: ",InputForm[temp]];*)
    vars=#[[1]]&/@Cases[Variables[temp],x[_]];
    subsets=Subsets[vars,{2}];
    variants={};
    (*mins=Table[Infinity,{Length[subsets]}];
    rs=Table[{},{Length[subsets]}];*)
    For[i=1,i<=Length[subsets],i++,
        If[And[Length[NegTerms[AdvancedVarDiffReplace[UF[[2]],subsets[[i]]]]]<Length[temp],
               Length[NegTerms[AdvancedVarDiffReplace[UF[[2]],{subsets[[i]][[2]],subsets[[i]][[1]]}]]]<Length[temp]]
                ,
(*Print[subsets[[i]]];
Print["in"];*)
            r1=AdvancedKillNegTerms3[AdvancedVarDiffReplace[{{#[[1]]/2,#[[2]]}&/@Z,MM,UF},subsets[[i]]],level+1];
            r2=AdvancedKillNegTerms3[AdvancedVarDiffReplace[{{#[[1]]/2,#[[2]]}&/@Z,MM,UF},{subsets[[i]][[2]],subsets[[i]][[1]]}],level+1];
(*If[level==1,AAAr1=r1;
AAAr2=r2;Abort[]];*)



(*Print["out"];*)
	    pairs=Tuples[{r1,r2}];
	    variants=Join[variants,{Join[##[[1]][[1]],##[[2]][[1]]],##[[1]][[2]]+##[[2]][[2]]}&/@pairs];
(*If[level==1,Print[Length[variants[[1]][[1]]]];Abort[];];*)

        ];
    ];
      If[Length[variants]==0,Return[{{{ZZ},Length[temp]}}]];
	If[level===1,
	  min=Min@@((##[[2]])&/@variants);
	  i=Position[((##[[2]])&/@variants),min][[1]][[1]];	
          RawPrintLn["Total number of negative terms remaining in subexpressions: ",variants[[i]][[2]]];
	  RawPrintLn["Total number of subexpressions: ",Length[variants[[i]][[1]]]];
	  Return[{variants[[i]]}];	
	,
	  Return[variants];
	];  
   
]


SeparateNegTerms[ZZ_]:=Module[{temp,vars,subsets,U,F,UF,Z,i},
    {Z,MM,UF}=ZZ;
    temp=Expand[UF[[2]]];
    temp=If[Head[temp]===Plus,List@@temp,{temp}];
    temp1=Select[temp, ((## /. {x[ii_] -> 1, y[ii_] -> 1}) < 0) &];
    If[And[Length[temp1]>0,Length[temp]>Length[temp1]],
        RawPrintLn["Terms of different signs encountered"];
        UF={UF[[1]],Expand[UF[[2]]-Plus@@temp1],-Plus@@temp1};
        MM=MM*(Exp[-I Pi z]);
        Z={##[[1]]*Gamma[-##[[2]][[2]]+z]*Gamma[-z]/Gamma[-##[[2]][[2]]], {##[[2]][[1]], ##[[2]][[2]]-z,z}} & /@Z;
        HASZ=True;
    ];
(*        Z={##[[1]]*Gamma[-##[[2]][[2]]], ##[[2]]} & /@Z;*)
  (*      OUTSIDE=OUTSIDE*Gamma[-Z[[1]][[2]][[2]]];
        If[(-Z[[1]][[2]][[2]]/.ep->0)<=0,RUNORDER++];*)

    {Z,MM,UF}
]

DecomposeF[ZZ_]:=Module[{temp,vars,subsets,U,F,UF,Z,i,shift},

    {Z,MM,UF}=ZZ;
    temp=Expand[UF[[2]]];
    temp1=Coefficient[temp,ExV,1];
    temp1=Expand[temp1];
    temp2=Coefficient[temp,ExV,2];
    temp=If[Head[temp]===Plus,List@@temp,{temp}];
    temp1=If[Head[temp1]===Plus,List@@temp1,{temp1}];
    temp2=If[Head[temp2]===Plus,List@@temp2,{temp2}];


    If[And[Length[temp1]>0,Length[temp]>Length[temp1]],
        RawPrintLn["Expansion variable found"];
	shift=0;
	While[Expand[UF[[1]]/.(ExV->0)]===0,
	  shift+=1;
	  UF[[1]]=Expand[UF[[1]]/ExV];
	];
(*Print[shift];*)
MM*=ExV^shift;
        
        MM=MM*(ExV^(z));
	If[temp2==={0},
	  UF=Join[{UF[[1]]/.(ExV->0),Expand[UF[[2]]-(ExV* Plus@@temp1)]//.(ExV->0),Plus@@temp1},Drop[UF,2]];  
	  Z={##[[1]]*Gamma[-##[[2]][[2]]+z]*Gamma[-z]/Gamma[-##[[2]][[2]]], Join[{##[[2]][[1]], ##[[2]][[2]]-z,z},Drop[##[[2]],2]]} & /@Z;
	  Power2=False;
	,
	  MM=MM*(ExV^(2 z2));
	  UF=Join[{UF[[1]]/.(ExV->0),Expand[UF[[2]]-(ExV* Plus@@temp1)-(ExV^2* Plus@@temp2)]//.(ExV->0),Plus@@temp1,Plus@@temp2},Drop[UF,2]];  
	  Z={##[[1]]*Gamma[-##[[2]][[2]]+z+z2]*Gamma[-z]*Gamma[-z2]/Gamma[-##[[2]][[2]]], Join[{##[[2]][[1]], ##[[2]][[2]]-z-z2,z,z2},Drop[##[[2]],2]]} & /@Z;
	  Power2=True;
	];

        HASZ=True;
    ];
(*Print[{Z,MM,UF}];*)
    If[Power2,
      {Z,MM,UF}={Z,MM,UF}//.{z->2u2,z2->(u/2-u2)};
(*{Z,MM,UF}={Z,MM,UF}//.{z->(u-2u2)/3,z2->(u+u2)/3};*)
      {Z,MM,UF}=ExpandAll[{Z,MM,UF}//.{u2->z2,u->z}];
     ];



    If[False,
      {Z,MM,UF}={Z,MM,UF}//.{z->z3};
      {Z,MM,UF}={Z,MM,UF}//.{z2->z};
      {Z,MM,UF}={Z,MM,UF}//.{z3->z2};
     ];
(*Print[{Z,MM,UF}];*)

    {Z,MM,UF}
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



SDExpandG[{graph_,external_},{U_,F_,h_},degrees1_,order_,exv_,exo_]:=Module[{temp},
    ExpandDirect=False;
    ExpandMode=True;
    DMode=False;
    ExpandVariable=exv;
    ExpandOrder=exo;
    SDEvaluateG1[{graph,external},{U,F,h},degrees1,order]
]

SDExpand[{U_,F_,h_},degrees1_,order_,exv_,exo_]:=Module[{temp},
    ExpandDirect=False;
    ExpandMode=True;
    DMode=False;
    ExpandVariable=exv;
    ExpandOrder=exo;
    SDEvaluate1[{U,F,h},degrees1,order]
]



SDEvaluateG[{graph_,external_},{U_,F_,h_},degrees1_,order_]:=Module[{temp},
    ExpandDirect=False;
    ExpandMode=False;
    DMode=False;
    SDEvaluateG1[{graph,external},{U,F,h},degrees1,order]
]

SDEvaluate[{U_,F_,h_},degrees1_,order_]:=Module[{temp},
    ExpandDirect=False;
    ExpandMode=False;
    DMode=False;
    SDEvaluate1[{U,F,h},degrees1,order]
]

SDAnalyze[{U_,F_,h_},degrees1_,order_,dmin_,dmax_]:=Module[{temp},
    ExpandDirect=False;
    ExpandMode=False;
    DMode=True;
    DMIN=dmin;
    DMAX=dmax;
    SDEvaluate1[{U,F,h},degrees1,order]
]



SDEvaluateG1[{graph_,external_},{U_,F_,h_},degrees1_,order_]:=Module[{temp},
    ExpandDirect=False;
    GraphInput=True;
    InfinityVertex=Apply[Max,graph,{0,Infinity}]+1;
    CurrentGraph=MyGraph[Join[graph,{##,InfinityVertex}&/@external]];
    Exte=Length[external];
    External=external;
    SDEvaluate1[{U,F,h},degrees1,order]
]


SDEvaluate1[{U_,F_,h_},degrees1_,order_]:=Module[{ZZ,n,degrees,epdegrees,runorder,outside,a1,deltas},
    ExpandDirect=False;
    If[GraphInput,GraphUsed=True,GraphUsed=False,GraphUsed=False];
    InitializeLinks[];
    GraphInput=False;
    n=Length[degrees1];
    If[GraphUsed,
        If[Not[M[CurrentGraph]-Exte===Length[degrees1]],
            RawPrintLn["ERROR: length of indices is different from the number of lines in the graph"];
            Abort[];
        ];
    ];
    If[Length[NegTerms[F]]==Length[Expand[F]],
      RawPrintLn["!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"];
      RawPrintLn["All terms in F are negative. Please consider to change the sign of all propagators"];
      ,
      If[Length[NegTerms[F]]>Length[Expand[F]]/2,
	RawPrintLn["There are more negative terms in F than positive ones. Please consider to change the sign of all propagators"]
      ]
    ];

    degrees=degrees1/.{ep->0};
    epdegrees=Expand[degrees1-(degrees1/.{ep->0})]/ep;
    a1=Plus@@degrees1;
    ZZ={{{1,{a1-(h+1)(d0/2-ep),-(a1-h(d0/2-ep))}}},{U,F}};
    runorder=order;
    For[j=1,j<=Length[degrees],j++,
        If[And[degrees[[j]]<=0,Not[epdegrees[[j]]==0]],runorder--]
    ];
    outside=(E^(h EulerGamma ep));

    If[((a1-(h)(d0/2-ep))/.ep->0)<=0,runorder++];
    outside=outside*Gamma[a1-h(d0/2-ep)];

    deltas={Range[n]};
    {ZZ,deltas,degrees,epdegrees}=KillNegativeIndices[ZZ,deltas,degrees,epdegrees];

    SDEvaluate2[ZZ,degrees,epdegrees,runorder,order,outside,deltas]
]
SDEvaluateDirect[var_, functions_, degrees_, order_] := SDEvaluateDirect[var, functions, degrees, order,{}];

SDEvaluateDirect[var_, functions_, degrees_, order_,deltas_] :=  Module[{temp, n, ZZ, result, intvars,i},
        If[GraphInput,GraphUsed=True,GraphUsed=False,GraphUsed=False];
            GraphInput=False;
    ExpandMode=False;
    DMode=False;
ExpandDirect=False;
  VariablesHere = Sort[Union[Cases[functions, var[_], {0, Infinity}]]];
  If[Length[VariablesHere]===0,
    n=0,
    n = Last[VariablesHere][[1]]
  ];
  temp = Transpose[{functions, degrees}];
  ZZ = {{{{Times @@ ((##[[1]]) & /@ Select[temp, (##[[2]] === 1) &]), DeleteCases[degrees, 1]}}, 
                 Rationalize[(##[[1]]) & /@ Select[temp, (##[[2]] =!= 1) &]]/.var->x}};
  If[Not[MemberQ[Variables[degrees],delta[_]]],
     (* Print[ZZ];*)
  For[i=1,i<=Length[BisectionVariables],i++,    
    If[MemberQ[Variables[ZZ],x[BisectionVariables[[i]]]],
      ZZ=Flatten[{{{{##[[1]][[1]][[1]]/2,##[[1]][[1]][[2]]}},Expand[##[[2]]/.(x[BisectionVariables[[i]]]->x[BisectionVariables[[i]]]/2)]},
                {{{##[[1]][[1]][[1]]/2,##[[1]][[1]][[2]]}},Expand[##[[2]]/.(x[BisectionVariables[[i]]]->1-x[BisectionVariables[[i]]]/2)]}}
            &/@ZZ,1];
    ]
  ]
];
(*  Print[ZZ];*)
 (* Abort[];*)
  OUTSIDE = 1;
  DENOM = 1;
  ORDER = order;
  RUNORDER = order;
  HASZ=False;
  (*deltas = {};*)
      InitializeLinks[];
  intvars = Table[1, {n}];
  SDIntegrate[intvars, ZZ, RUNORDER, deltas];
  result =
   GenerateAnswer[intvars, ZZ, RUNORDER, ORDER, deltas, OUTSIDE, DENOM];
   If[RemoveDatabases, ClearDatabaseWrapper[1, False]; ClearDatabaseWrapper[2, False];If[UsingQLink,Quiet[DeleteFile[DataPath<>"3"]]]];
  Plus @@ (((##[[1]])*##[[2]]) & /@ result)
 ]

SDExpandDirect[var_, functions_, degrees_, exv_]:=SDExpandDirect[var, functions, degrees, exv,{}]


SDExpandDirect[var_, functions_, degrees_, exv_,deltas_] :=  Module[{temp, VVariables, n, ZZ, result, intvars,i},
        If[GraphInput,GraphUsed=True,GraphUsed=False,GraphUsed=False];
            GraphInput=False;
    ExpandMode=True;
    DMode=False;
    ExpandVariable=exv;
    ExpandOrder=0;
    ExpandDirect=True;


  VVariables = Sort[Union[Cases[functions, var[_], {0, Infinity}]]];
  If[Length[VVariables]===0,Abort[]];
  n = Last[VVariables][[1]];
  ZZ = {{{{1, degrees}}, 1, Rationalize[functions]/.var->x}};



  OUTSIDE = 1;
  DENOM = 1;
  ORDER = Infinity;
  RUNORDER = 0;
  
      InitializeLinks[];
  intvars = Table[1, {n}];

    ZZ=ZZ/.ExpandVariable->ExV;

        ZZ=DecomposeF/@ZZ;
ZZ=UseEpMonom/@ZZ;



(*  Print[ZZ];*)
  For[i=1,i<=Length[BisectionVariables],i++,
If[MemberQ[Variables[ZZ],x[BisectionVariables[[i]]]],
    
    ZZ=Flatten[{{{{##[[1]][[1]][[1]]/2,##[[1]][[1]][[2]]}},Expand[##[[2]]/.(x[BisectionVariables[[i]]]->x[BisectionVariables[[i]]]/2)]},
                {{{##[[1]][[1]][[1]]/2,##[[1]][[1]][[2]]}},Expand[##[[2]]/.(x[BisectionVariables[[i]]]->1-x[BisectionVariables[[i]]]/2)]}}
            &/@ZZ,1];
]
  ];

  Return[SDIntegrate[intvars, ZZ, RUNORDER, deltas]];



  (*result =
   GenerateAnswer[intvars, ZZ, RUNORDER, ORDER, deltas, OUTSIDE, DENOM];
   If[RemoveDatabases, ClearDatabaseWrapper[1, False]; ClearDatabaseWrapper[2, False];If[UsingQLink,Quiet[DeleteFile[DataPath<>"3"]]]];
  Plus @@ (((##[[1]])*##[[2]]) & /@ result)*)
 ]




SDEvaluate2[ZZ1_,degrees_,epdegrees_,runorder_,order_,outside_,deltas_]:=Module[{temp,i,j,var,ZZ,aaa},

    Clear[TermNumbersForExpandOrder];
    ZZ=ZZ1;
    If[Length[ZZ[[1]]]===0,Return[0]];
    n=Length[degrees];
    intvars=Table[(If[Or[degrees[[i]]>0,Not[epdegrees[[i]]===0]],1,0]),{i,1,n}];
    ZZ[[1]]={#[[1]]*(Times@@Table[x[i]^(If[degrees[[i]]>0,degrees[[i]]-1,0]),{i,1,n}]),#[[2]]}&/@ZZ[[1]];
    ZZ={{ZZ[[1]],(Times@@Table[x[i]^(If[intvars[[i]]===1,epdegrees[[i]]*ep+If[degrees[[i]]<=0,degrees[[i]]-1,0],0]),{i,1,n}]),ZZ[[2]]}};

    OUTSIDE=outside;
    ORDER=order;
    RUNORDER=runorder;


    If[And[ExpandMode,NegativeTermsHandling==="MB"],
        RaPrintLn["Cannot use MB negative terms handling in expand mode, set to Squares"];
        NegativeTermsHandling="Squares"
    ];

    HASZ=False;

    If[NegativeTermsHandling==="Squares",
        While[True,
            ZZOld=ZZ;
            ZZ=Flatten[KillNegTerms/@ZZ,1];
            If[ZZ===ZZOld,Break[]];
        ];
    ];
    
    If[NegativeTermsHandling==="AdvancedSquares",
        ZZ=Flatten[AdvancedKillNegTerms[##,1][[1]]&/@ZZ,1]
    ];

    If[NegativeTermsHandling==="AdvancedSquares2",
        ZZ=Flatten[AdvancedKillNegTerms2[##,1][[1]][[1]]&/@ZZ,1]
    ];    

    If[NegativeTermsHandling==="AdvancedSquares3",
        ZZ=Flatten[AdvancedKillNegTerms3[##,1][[1]][[1]]&/@ZZ,1]
    ];    
        

        ZZ=ZZ/.ExpandVariable->ExV;
    If[ExpandMode,
        ZZ=DecomposeF/@ZZ
    ];

    If[NegativeTermsHandling==="MB",
        ZZ=SeparateNegTerms/@ZZ;
    ];

    ZZ=UseEpMonom/@ZZ;

   (* Print[ZZ];*)
   (* Abort[];*)
    DENOM=Times@@Table[If[And[degrees[[i]]<=0,epdegrees[[i]]==0],1,Gamma[degrees[[i]]+ep*epdegrees[[i]]]],{i,Length[degrees]}];
    RawPrintLn["Integration has to be performed up to order ",runorder];
   (* ZZZ145=ZZ;
    Abort[];*)
    aaa=SDIntegrate[intvars,ZZ,RUNORDER,deltas];
    If[DMode,Return[aaa]];
    result=GenerateAnswer[intvars,ZZ,RUNORDER,ORDER,deltas,OUTSIDE,DENOM];
    If[RemoveDatabases,ClearDatabaseWrapper[1,False];ClearDatabaseWrapper[2,False]];
    Plus@@(((##[[1]])*##[[2]])&/@result)

]

GenerateAnswer[intvars_,ZZ_,runorder_,order_,deltas_,outside_,denom_]:=Module[{temp,i,j,min,degs,result},


    degs=Union[##[[2]]&/@ForEvaluation];
    result=Flatten[(
        For[min=runorder-30,min<=runorder,min++,If[And[Head[SDIntegral[intvars,ZZ,deltas,min,##]]=!=SDIntegral,SDIntegral[intvars,ZZ,deltas,min,##]=!=0],Break[]]];
        If[min>runorder,{},
            temp=Series[outside/denom,{ep,0,order-min}];
            temp=Table[{
                    And@@Table[SDExact[intvars,ZZ,deltas,i,##],{i,min,j+runorder-order}],
                    Plus@@Table[SDIntegral[intvars,ZZ,deltas,i,##]*SeriesCoefficient[temp,j-i],{i,min,j+runorder-order}],
                ep^j*(ExpandVariable^##[[1]])*(Log[ExpandVariable]^##[[2]])},{j,min+order-runorder,order}];
            temp={If[##[[1]],FullSimplify[##[[2]]],DoubleCutExtraDigits[PMSimplify[##[[2]]]]],##[[3]]}&/@temp;
            temp
        ]
   )&/@degs,1];
    Return[result];

(*        Return[Plus@@(((##[[1]])*##[[2]])&/@temp)];*)

 (*   temp=Series[outside*(Plus@@Table[(SDIntegral[intvars,ZZ,deltas,j])*ep^j,{j,runorder-30,runorder}])/denom,{ep,0,order}];
    If[temp===0,Return[0]];
    If[temp===Indeterminate,Return[Indeterminate]];
    temp[[3]]=PMSimplify[##,ReturnErrorWithBrackets]&/@temp[[3]];
    temp[[3]]=DoubleCutExtraDigits/@temp[[3]];
    temp*)
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

MyDegrees[pol_]:=MyDegrees[pol,Sort[Cases[Variables[pol],x[_]]]];

MyDegrees[pol_,var_] := Module[{rule, temp, degrees, i},
(*    var=Array[x,n];*)
  rule = Rule[##, 1] & /@ var;
  temp = Expand[pol/.(0.->0)];
  If[Head[temp]===Plus,temp = List @@ temp,temp={temp}];
  
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


MLN[xx_] := Module[{n, vectors, i, j, ll, nn,l,v},
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
  If[xx[[1]] > yy[[1]], False, 
   If[xx[[2]] < yy[[2]], True, 
    If[xx[[2]] > yy[[2]], False, 
     If[xx[[3]] < yy[[3]], True, False]]]]]


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



MBshiftRules[dx_] := {
  Gamma[m_.+a_.*dx] :> Gamma[1+a*dx]/Product[a*dx-i, {i,0,-m}] /;
  IntegerQ[m] && m <= 0,

  PolyGamma[n_, m_.+a_.*dx] :> (-a*dx)^(-n-1)*
  ((-a*dx)^(n+1)*PolyGamma[n, 1+a*dx] + n!*Sum[(a*dx/(a*dx-i))^(n+1),
  {i,0,-m}]) /; IntegerQ[m] && m <= 0};

MBexpansionRules::series = "exhausted precomputed expansion of Gamma's (`1`)";

MBexpansionRules[dx_, order_] := {
  Gamma[m_+a_.*dx] :> If[order <= 20,
  Gamma[m]*Sum[(a*dx)^i/i!*MBexpGam[m, i], {i,0,order}],
  Message[MBexpansionRules::series, order];
  Normal[Series[Gamma[m+a*dx], {dx,0,order}]]] /; !IntegerQ[m] || m > 0,

  PolyGamma[n_, m_+a_.*dx] :> Sum[(a*dx)^i/i!*PolyGamma[n+i, m],
  {i,0,order}] /; !IntegerQ[m] || m > 0};





(* generated automatically with MATHEMATICA *)

MBexpGam[a_, 0] = 1;

MBexpGam[a_, 1] = PolyGamma[0, a];

MBexpGam[a_, 2] = PolyGamma[0, a]^2 + PolyGamma[1, a];

MBexpGam[a_, 3] = PolyGamma[0, a]^3 + 3*PolyGamma[0, a]*PolyGamma[1, a] +
     PolyGamma[2, a];

MBexpGam[a_, 4] = PolyGamma[0, a]^4 + 6*PolyGamma[0, a]^2*PolyGamma[1, a] +
     3*PolyGamma[1, a]^2 + 4*PolyGamma[0, a]*PolyGamma[2, a] + PolyGamma[3, a];

MBexpGam[a_, 5] = PolyGamma[0, a]^5 + 10*PolyGamma[0, a]^3*PolyGamma[1, a] +
     15*PolyGamma[0, a]*PolyGamma[1, a]^2 + 10*PolyGamma[0, a]^2*
      PolyGamma[2, a] + 10*PolyGamma[1, a]*PolyGamma[2, a] +
     5*PolyGamma[0, a]*PolyGamma[3, a] + PolyGamma[4, a];

MBexpGam[a_, 6] = PolyGamma[0, a]^6 + 15*PolyGamma[0, a]^4*PolyGamma[1, a] +
     45*PolyGamma[0, a]^2*PolyGamma[1, a]^2 + 15*PolyGamma[1, a]^3 +
     20*PolyGamma[0, a]^3*PolyGamma[2, a] + 60*PolyGamma[0, a]*
      PolyGamma[1, a]*PolyGamma[2, a] + 10*PolyGamma[2, a]^2 +
     15*PolyGamma[0, a]^2*PolyGamma[3, a] + 15*PolyGamma[1, a]*
      PolyGamma[3, a] + 6*PolyGamma[0, a]*PolyGamma[4, a] + PolyGamma[5, a];

MBexpGam[a_, 7] = PolyGamma[0, a]^7 + 21*PolyGamma[0, a]^5*PolyGamma[1, a] +
     105*PolyGamma[0, a]^3*PolyGamma[1, a]^2 + 105*PolyGamma[0, a]*
      PolyGamma[1, a]^3 + 35*PolyGamma[0, a]^4*PolyGamma[2, a] +
     210*PolyGamma[0, a]^2*PolyGamma[1, a]*PolyGamma[2, a] +
     105*PolyGamma[1, a]^2*PolyGamma[2, a] + 70*PolyGamma[0, a]*
      PolyGamma[2, a]^2 + 35*PolyGamma[0, a]^3*PolyGamma[3, a] +
     105*PolyGamma[0, a]*PolyGamma[1, a]*PolyGamma[3, a] +
     35*PolyGamma[2, a]*PolyGamma[3, a] + 21*PolyGamma[0, a]^2*
      PolyGamma[4, a] + 21*PolyGamma[1, a]*PolyGamma[4, a] +
     7*PolyGamma[0, a]*PolyGamma[5, a] + PolyGamma[6, a];

MBexpGam[a_, 8] = PolyGamma[0, a]^8 + 28*PolyGamma[0, a]^6*PolyGamma[1, a] +
     210*PolyGamma[0, a]^4*PolyGamma[1, a]^2 + 420*PolyGamma[0, a]^2*
      PolyGamma[1, a]^3 + 105*PolyGamma[1, a]^4 + 56*PolyGamma[0, a]^5*
      PolyGamma[2, a] + 560*PolyGamma[0, a]^3*PolyGamma[1, a]*
      PolyGamma[2, a] + 840*PolyGamma[0, a]*PolyGamma[1, a]^2*
      PolyGamma[2, a] + 280*PolyGamma[0, a]^2*PolyGamma[2, a]^2 +
     280*PolyGamma[1, a]*PolyGamma[2, a]^2 + 70*PolyGamma[0, a]^4*
      PolyGamma[3, a] + 420*PolyGamma[0, a]^2*PolyGamma[1, a]*
      PolyGamma[3, a] + 210*PolyGamma[1, a]^2*PolyGamma[3, a] +
     280*PolyGamma[0, a]*PolyGamma[2, a]*PolyGamma[3, a] +
     35*PolyGamma[3, a]^2 + 56*PolyGamma[0, a]^3*PolyGamma[4, a] +
     168*PolyGamma[0, a]*PolyGamma[1, a]*PolyGamma[4, a] +
     56*PolyGamma[2, a]*PolyGamma[4, a] + 28*PolyGamma[0, a]^2*
      PolyGamma[5, a] + 28*PolyGamma[1, a]*PolyGamma[5, a] +
     8*PolyGamma[0, a]*PolyGamma[6, a] + PolyGamma[7, a];

MBexpGam[a_, 9] = PolyGamma[0, a]^9 + 36*PolyGamma[0, a]^7*PolyGamma[1, a] +
     378*PolyGamma[0, a]^5*PolyGamma[1, a]^2 + 1260*PolyGamma[0, a]^3*
      PolyGamma[1, a]^3 + 945*PolyGamma[0, a]*PolyGamma[1, a]^4 +
     84*PolyGamma[0, a]^6*PolyGamma[2, a] + 1260*PolyGamma[0, a]^4*
      PolyGamma[1, a]*PolyGamma[2, a] + 3780*PolyGamma[0, a]^2*
      PolyGamma[1, a]^2*PolyGamma[2, a] + 1260*PolyGamma[1, a]^3*
      PolyGamma[2, a] + 840*PolyGamma[0, a]^3*PolyGamma[2, a]^2 +
     2520*PolyGamma[0, a]*PolyGamma[1, a]*PolyGamma[2, a]^2 +
     280*PolyGamma[2, a]^3 + 126*PolyGamma[0, a]^5*PolyGamma[3, a] +
     1260*PolyGamma[0, a]^3*PolyGamma[1, a]*PolyGamma[3, a] +
     1890*PolyGamma[0, a]*PolyGamma[1, a]^2*PolyGamma[3, a] +
     1260*PolyGamma[0, a]^2*PolyGamma[2, a]*PolyGamma[3, a] +
     1260*PolyGamma[1, a]*PolyGamma[2, a]*PolyGamma[3, a] +
     315*PolyGamma[0, a]*PolyGamma[3, a]^2 + 126*PolyGamma[0, a]^4*
      PolyGamma[4, a] + 756*PolyGamma[0, a]^2*PolyGamma[1, a]*
      PolyGamma[4, a] + 378*PolyGamma[1, a]^2*PolyGamma[4, a] +
     504*PolyGamma[0, a]*PolyGamma[2, a]*PolyGamma[4, a] +
     126*PolyGamma[3, a]*PolyGamma[4, a] + 84*PolyGamma[0, a]^3*
      PolyGamma[5, a] + 252*PolyGamma[0, a]*PolyGamma[1, a]*PolyGamma[5, a] +
     84*PolyGamma[2, a]*PolyGamma[5, a] + 36*PolyGamma[0, a]^2*
      PolyGamma[6, a] + 36*PolyGamma[1, a]*PolyGamma[6, a] +
     9*PolyGamma[0, a]*PolyGamma[7, a] + PolyGamma[8, a];

MBexpGam[a_, 10] = PolyGamma[0, a]^10 + 45*PolyGamma[0, a]^8*
      PolyGamma[1, a] + 630*PolyGamma[0, a]^6*PolyGamma[1, a]^2 +
     3150*PolyGamma[0, a]^4*PolyGamma[1, a]^3 + 4725*PolyGamma[0, a]^2*
      PolyGamma[1, a]^4 + 945*PolyGamma[1, a]^5 + 120*PolyGamma[0, a]^7*
      PolyGamma[2, a] + 2520*PolyGamma[0, a]^5*PolyGamma[1, a]*
      PolyGamma[2, a] + 12600*PolyGamma[0, a]^3*PolyGamma[1, a]^2*
      PolyGamma[2, a] + 12600*PolyGamma[0, a]*PolyGamma[1, a]^3*
      PolyGamma[2, a] + 2100*PolyGamma[0, a]^4*PolyGamma[2, a]^2 +
     12600*PolyGamma[0, a]^2*PolyGamma[1, a]*PolyGamma[2, a]^2 +
     6300*PolyGamma[1, a]^2*PolyGamma[2, a]^2 + 2800*PolyGamma[0, a]*
      PolyGamma[2, a]^3 + 210*PolyGamma[0, a]^6*PolyGamma[3, a] +
     3150*PolyGamma[0, a]^4*PolyGamma[1, a]*PolyGamma[3, a] +
     9450*PolyGamma[0, a]^2*PolyGamma[1, a]^2*PolyGamma[3, a] +
     3150*PolyGamma[1, a]^3*PolyGamma[3, a] + 4200*PolyGamma[0, a]^3*
      PolyGamma[2, a]*PolyGamma[3, a] + 12600*PolyGamma[0, a]*PolyGamma[1, a]*
      PolyGamma[2, a]*PolyGamma[3, a] + 2100*PolyGamma[2, a]^2*
      PolyGamma[3, a] + 1575*PolyGamma[0, a]^2*PolyGamma[3, a]^2 +
     1575*PolyGamma[1, a]*PolyGamma[3, a]^2 + 252*PolyGamma[0, a]^5*
      PolyGamma[4, a] + 2520*PolyGamma[0, a]^3*PolyGamma[1, a]*
      PolyGamma[4, a] + 3780*PolyGamma[0, a]*PolyGamma[1, a]^2*
      PolyGamma[4, a] + 2520*PolyGamma[0, a]^2*PolyGamma[2, a]*
      PolyGamma[4, a] + 2520*PolyGamma[1, a]*PolyGamma[2, a]*
      PolyGamma[4, a] + 1260*PolyGamma[0, a]*PolyGamma[3, a]*
      PolyGamma[4, a] + 126*PolyGamma[4, a]^2 + 210*PolyGamma[0, a]^4*
      PolyGamma[5, a] + 1260*PolyGamma[0, a]^2*PolyGamma[1, a]*
      PolyGamma[5, a] + 630*PolyGamma[1, a]^2*PolyGamma[5, a] +
     840*PolyGamma[0, a]*PolyGamma[2, a]*PolyGamma[5, a] +
     210*PolyGamma[3, a]*PolyGamma[5, a] + 120*PolyGamma[0, a]^3*
      PolyGamma[6, a] + 360*PolyGamma[0, a]*PolyGamma[1, a]*PolyGamma[6, a] +
     120*PolyGamma[2, a]*PolyGamma[6, a] + 45*PolyGamma[0, a]^2*
      PolyGamma[7, a] + 45*PolyGamma[1, a]*PolyGamma[7, a] +
     10*PolyGamma[0, a]*PolyGamma[8, a] + PolyGamma[9, a];

MBexpGam[a_, 11] = PolyGamma[0, a]^11 + 55*PolyGamma[0, a]^9*
      PolyGamma[1, a] + 990*PolyGamma[0, a]^7*PolyGamma[1, a]^2 +
     6930*PolyGamma[0, a]^5*PolyGamma[1, a]^3 + 17325*PolyGamma[0, a]^3*
      PolyGamma[1, a]^4 + 10395*PolyGamma[0, a]*PolyGamma[1, a]^5 +
     165*PolyGamma[0, a]^8*PolyGamma[2, a] + 4620*PolyGamma[0, a]^6*
      PolyGamma[1, a]*PolyGamma[2, a] + 34650*PolyGamma[0, a]^4*
      PolyGamma[1, a]^2*PolyGamma[2, a] + 69300*PolyGamma[0, a]^2*
      PolyGamma[1, a]^3*PolyGamma[2, a] + 17325*PolyGamma[1, a]^4*
      PolyGamma[2, a] + 4620*PolyGamma[0, a]^5*PolyGamma[2, a]^2 +
     46200*PolyGamma[0, a]^3*PolyGamma[1, a]*PolyGamma[2, a]^2 +
     69300*PolyGamma[0, a]*PolyGamma[1, a]^2*PolyGamma[2, a]^2 +
     15400*PolyGamma[0, a]^2*PolyGamma[2, a]^3 + 15400*PolyGamma[1, a]*
      PolyGamma[2, a]^3 + 330*PolyGamma[0, a]^7*PolyGamma[3, a] +
     6930*PolyGamma[0, a]^5*PolyGamma[1, a]*PolyGamma[3, a] +
     34650*PolyGamma[0, a]^3*PolyGamma[1, a]^2*PolyGamma[3, a] +
     34650*PolyGamma[0, a]*PolyGamma[1, a]^3*PolyGamma[3, a] +
     11550*PolyGamma[0, a]^4*PolyGamma[2, a]*PolyGamma[3, a] +
     69300*PolyGamma[0, a]^2*PolyGamma[1, a]*PolyGamma[2, a]*
      PolyGamma[3, a] + 34650*PolyGamma[1, a]^2*PolyGamma[2, a]*
      PolyGamma[3, a] + 23100*PolyGamma[0, a]*PolyGamma[2, a]^2*
      PolyGamma[3, a] + 5775*PolyGamma[0, a]^3*PolyGamma[3, a]^2 +
     17325*PolyGamma[0, a]*PolyGamma[1, a]*PolyGamma[3, a]^2 +
     5775*PolyGamma[2, a]*PolyGamma[3, a]^2 + 462*PolyGamma[0, a]^6*
      PolyGamma[4, a] + 6930*PolyGamma[0, a]^4*PolyGamma[1, a]*
      PolyGamma[4, a] + 20790*PolyGamma[0, a]^2*PolyGamma[1, a]^2*
      PolyGamma[4, a] + 6930*PolyGamma[1, a]^3*PolyGamma[4, a] +
     9240*PolyGamma[0, a]^3*PolyGamma[2, a]*PolyGamma[4, a] +
     27720*PolyGamma[0, a]*PolyGamma[1, a]*PolyGamma[2, a]*PolyGamma[4, a] +
     4620*PolyGamma[2, a]^2*PolyGamma[4, a] + 6930*PolyGamma[0, a]^2*
      PolyGamma[3, a]*PolyGamma[4, a] + 6930*PolyGamma[1, a]*PolyGamma[3, a]*
      PolyGamma[4, a] + 1386*PolyGamma[0, a]*PolyGamma[4, a]^2 +
     462*PolyGamma[0, a]^5*PolyGamma[5, a] + 4620*PolyGamma[0, a]^3*
      PolyGamma[1, a]*PolyGamma[5, a] + 6930*PolyGamma[0, a]*
      PolyGamma[1, a]^2*PolyGamma[5, a] + 4620*PolyGamma[0, a]^2*
      PolyGamma[2, a]*PolyGamma[5, a] + 4620*PolyGamma[1, a]*PolyGamma[2, a]*
      PolyGamma[5, a] + 2310*PolyGamma[0, a]*PolyGamma[3, a]*
      PolyGamma[5, a] + 462*PolyGamma[4, a]*PolyGamma[5, a] +
     330*PolyGamma[0, a]^4*PolyGamma[6, a] + 1980*PolyGamma[0, a]^2*
      PolyGamma[1, a]*PolyGamma[6, a] + 990*PolyGamma[1, a]^2*
      PolyGamma[6, a] + 1320*PolyGamma[0, a]*PolyGamma[2, a]*
      PolyGamma[6, a] + 330*PolyGamma[3, a]*PolyGamma[6, a] +
     165*PolyGamma[0, a]^3*PolyGamma[7, a] + 495*PolyGamma[0, a]*
      PolyGamma[1, a]*PolyGamma[7, a] + 165*PolyGamma[2, a]*PolyGamma[7, a] +
     55*PolyGamma[0, a]^2*PolyGamma[8, a] + 55*PolyGamma[1, a]*
      PolyGamma[8, a] + 11*PolyGamma[0, a]*PolyGamma[9, a] + PolyGamma[10, a]

MBexpGam[a_, 12] = PolyGamma[0, a]^12 + 66*PolyGamma[0, a]^10*
      PolyGamma[1, a] + 1485*PolyGamma[0, a]^8*PolyGamma[1, a]^2 +
     13860*PolyGamma[0, a]^6*PolyGamma[1, a]^3 + 51975*PolyGamma[0, a]^4*
      PolyGamma[1, a]^4 + 62370*PolyGamma[0, a]^2*PolyGamma[1, a]^5 +
     10395*PolyGamma[1, a]^6 + 220*PolyGamma[0, a]^9*PolyGamma[2, a] +
     7920*PolyGamma[0, a]^7*PolyGamma[1, a]*PolyGamma[2, a] +
     83160*PolyGamma[0, a]^5*PolyGamma[1, a]^2*PolyGamma[2, a] +
     277200*PolyGamma[0, a]^3*PolyGamma[1, a]^3*PolyGamma[2, a] +
     207900*PolyGamma[0, a]*PolyGamma[1, a]^4*PolyGamma[2, a] +
     9240*PolyGamma[0, a]^6*PolyGamma[2, a]^2 + 138600*PolyGamma[0, a]^4*
      PolyGamma[1, a]*PolyGamma[2, a]^2 + 415800*PolyGamma[0, a]^2*
      PolyGamma[1, a]^2*PolyGamma[2, a]^2 + 138600*PolyGamma[1, a]^3*
      PolyGamma[2, a]^2 + 61600*PolyGamma[0, a]^3*PolyGamma[2, a]^3 +
     184800*PolyGamma[0, a]*PolyGamma[1, a]*PolyGamma[2, a]^3 +
     15400*PolyGamma[2, a]^4 + 495*PolyGamma[0, a]^8*PolyGamma[3, a] +
     13860*PolyGamma[0, a]^6*PolyGamma[1, a]*PolyGamma[3, a] +
     103950*PolyGamma[0, a]^4*PolyGamma[1, a]^2*PolyGamma[3, a] +
     207900*PolyGamma[0, a]^2*PolyGamma[1, a]^3*PolyGamma[3, a] +
     51975*PolyGamma[1, a]^4*PolyGamma[3, a] + 27720*PolyGamma[0, a]^5*
      PolyGamma[2, a]*PolyGamma[3, a] + 277200*PolyGamma[0, a]^3*
      PolyGamma[1, a]*PolyGamma[2, a]*PolyGamma[3, a] +
     415800*PolyGamma[0, a]*PolyGamma[1, a]^2*PolyGamma[2, a]*
      PolyGamma[3, a] + 138600*PolyGamma[0, a]^2*PolyGamma[2, a]^2*
      PolyGamma[3, a] + 138600*PolyGamma[1, a]*PolyGamma[2, a]^2*
      PolyGamma[3, a] + 17325*PolyGamma[0, a]^4*PolyGamma[3, a]^2 +
     103950*PolyGamma[0, a]^2*PolyGamma[1, a]*PolyGamma[3, a]^2 +
     51975*PolyGamma[1, a]^2*PolyGamma[3, a]^2 + 69300*PolyGamma[0, a]*
      PolyGamma[2, a]*PolyGamma[3, a]^2 + 5775*PolyGamma[3, a]^3 +
     792*PolyGamma[0, a]^7*PolyGamma[4, a] + 16632*PolyGamma[0, a]^5*
      PolyGamma[1, a]*PolyGamma[4, a] + 83160*PolyGamma[0, a]^3*
      PolyGamma[1, a]^2*PolyGamma[4, a] + 83160*PolyGamma[0, a]*
      PolyGamma[1, a]^3*PolyGamma[4, a] + 27720*PolyGamma[0, a]^4*
      PolyGamma[2, a]*PolyGamma[4, a] + 166320*PolyGamma[0, a]^2*
      PolyGamma[1, a]*PolyGamma[2, a]*PolyGamma[4, a] +
     83160*PolyGamma[1, a]^2*PolyGamma[2, a]*PolyGamma[4, a] +
     55440*PolyGamma[0, a]*PolyGamma[2, a]^2*PolyGamma[4, a] +
     27720*PolyGamma[0, a]^3*PolyGamma[3, a]*PolyGamma[4, a] +
     83160*PolyGamma[0, a]*PolyGamma[1, a]*PolyGamma[3, a]*PolyGamma[4, a] +
     27720*PolyGamma[2, a]*PolyGamma[3, a]*PolyGamma[4, a] +
     8316*PolyGamma[0, a]^2*PolyGamma[4, a]^2 + 8316*PolyGamma[1, a]*
      PolyGamma[4, a]^2 + 924*PolyGamma[0, a]^6*PolyGamma[5, a] +
     13860*PolyGamma[0, a]^4*PolyGamma[1, a]*PolyGamma[5, a] +
     41580*PolyGamma[0, a]^2*PolyGamma[1, a]^2*PolyGamma[5, a] +
     13860*PolyGamma[1, a]^3*PolyGamma[5, a] + 18480*PolyGamma[0, a]^3*
      PolyGamma[2, a]*PolyGamma[5, a] + 55440*PolyGamma[0, a]*PolyGamma[1, a]*
      PolyGamma[2, a]*PolyGamma[5, a] + 9240*PolyGamma[2, a]^2*
      PolyGamma[5, a] + 13860*PolyGamma[0, a]^2*PolyGamma[3, a]*
      PolyGamma[5, a] + 13860*PolyGamma[1, a]*PolyGamma[3, a]*
      PolyGamma[5, a] + 5544*PolyGamma[0, a]*PolyGamma[4, a]*
      PolyGamma[5, a] + 462*PolyGamma[5, a]^2 + 792*PolyGamma[0, a]^5*
      PolyGamma[6, a] + 7920*PolyGamma[0, a]^3*PolyGamma[1, a]*
      PolyGamma[6, a] + 11880*PolyGamma[0, a]*PolyGamma[1, a]^2*
      PolyGamma[6, a] + 7920*PolyGamma[0, a]^2*PolyGamma[2, a]*
      PolyGamma[6, a] + 7920*PolyGamma[1, a]*PolyGamma[2, a]*
      PolyGamma[6, a] + 3960*PolyGamma[0, a]*PolyGamma[3, a]*
      PolyGamma[6, a] + 792*PolyGamma[4, a]*PolyGamma[6, a] +
     495*PolyGamma[0, a]^4*PolyGamma[7, a] + 2970*PolyGamma[0, a]^2*
      PolyGamma[1, a]*PolyGamma[7, a] + 1485*PolyGamma[1, a]^2*
      PolyGamma[7, a] + 1980*PolyGamma[0, a]*PolyGamma[2, a]*
      PolyGamma[7, a] + 495*PolyGamma[3, a]*PolyGamma[7, a] +
     220*PolyGamma[0, a]^3*PolyGamma[8, a] + 660*PolyGamma[0, a]*
      PolyGamma[1, a]*PolyGamma[8, a] + 220*PolyGamma[2, a]*PolyGamma[8, a] +
     66*PolyGamma[0, a]^2*PolyGamma[9, a] + 66*PolyGamma[1, a]*
      PolyGamma[9, a] + 12*PolyGamma[0, a]*PolyGamma[10, a] + PolyGamma[11, a]

MBexpGam[a_, 13] = PolyGamma[0, a]^13 + 78*PolyGamma[0, a]^11*
      PolyGamma[1, a] + 2145*PolyGamma[0, a]^9*PolyGamma[1, a]^2 +
     25740*PolyGamma[0, a]^7*PolyGamma[1, a]^3 + 135135*PolyGamma[0, a]^5*
      PolyGamma[1, a]^4 + 270270*PolyGamma[0, a]^3*PolyGamma[1, a]^5 +
     135135*PolyGamma[0, a]*PolyGamma[1, a]^6 + 286*PolyGamma[0, a]^10*
      PolyGamma[2, a] + 12870*PolyGamma[0, a]^8*PolyGamma[1, a]*
      PolyGamma[2, a] + 180180*PolyGamma[0, a]^6*PolyGamma[1, a]^2*
      PolyGamma[2, a] + 900900*PolyGamma[0, a]^4*PolyGamma[1, a]^3*
      PolyGamma[2, a] + 1351350*PolyGamma[0, a]^2*PolyGamma[1, a]^4*
      PolyGamma[2, a] + 270270*PolyGamma[1, a]^5*PolyGamma[2, a] +
     17160*PolyGamma[0, a]^7*PolyGamma[2, a]^2 + 360360*PolyGamma[0, a]^5*
      PolyGamma[1, a]*PolyGamma[2, a]^2 + 1801800*PolyGamma[0, a]^3*
      PolyGamma[1, a]^2*PolyGamma[2, a]^2 + 1801800*PolyGamma[0, a]*
      PolyGamma[1, a]^3*PolyGamma[2, a]^2 + 200200*PolyGamma[0, a]^4*
      PolyGamma[2, a]^3 + 1201200*PolyGamma[0, a]^2*PolyGamma[1, a]*
      PolyGamma[2, a]^3 + 600600*PolyGamma[1, a]^2*PolyGamma[2, a]^3 +
     200200*PolyGamma[0, a]*PolyGamma[2, a]^4 + 715*PolyGamma[0, a]^9*
      PolyGamma[3, a] + 25740*PolyGamma[0, a]^7*PolyGamma[1, a]*
      PolyGamma[3, a] + 270270*PolyGamma[0, a]^5*PolyGamma[1, a]^2*
      PolyGamma[3, a] + 900900*PolyGamma[0, a]^3*PolyGamma[1, a]^3*
      PolyGamma[3, a] + 675675*PolyGamma[0, a]*PolyGamma[1, a]^4*
      PolyGamma[3, a] + 60060*PolyGamma[0, a]^6*PolyGamma[2, a]*
      PolyGamma[3, a] + 900900*PolyGamma[0, a]^4*PolyGamma[1, a]*
      PolyGamma[2, a]*PolyGamma[3, a] + 2702700*PolyGamma[0, a]^2*
      PolyGamma[1, a]^2*PolyGamma[2, a]*PolyGamma[3, a] +
     900900*PolyGamma[1, a]^3*PolyGamma[2, a]*PolyGamma[3, a] +
     600600*PolyGamma[0, a]^3*PolyGamma[2, a]^2*PolyGamma[3, a] +
     1801800*PolyGamma[0, a]*PolyGamma[1, a]*PolyGamma[2, a]^2*
      PolyGamma[3, a] + 200200*PolyGamma[2, a]^3*PolyGamma[3, a] +
     45045*PolyGamma[0, a]^5*PolyGamma[3, a]^2 + 450450*PolyGamma[0, a]^3*
      PolyGamma[1, a]*PolyGamma[3, a]^2 + 675675*PolyGamma[0, a]*
      PolyGamma[1, a]^2*PolyGamma[3, a]^2 + 450450*PolyGamma[0, a]^2*
      PolyGamma[2, a]*PolyGamma[3, a]^2 + 450450*PolyGamma[1, a]*
      PolyGamma[2, a]*PolyGamma[3, a]^2 + 75075*PolyGamma[0, a]*
      PolyGamma[3, a]^3 + 1287*PolyGamma[0, a]^8*PolyGamma[4, a] +
     36036*PolyGamma[0, a]^6*PolyGamma[1, a]*PolyGamma[4, a] +
     270270*PolyGamma[0, a]^4*PolyGamma[1, a]^2*PolyGamma[4, a] +
     540540*PolyGamma[0, a]^2*PolyGamma[1, a]^3*PolyGamma[4, a] +
     135135*PolyGamma[1, a]^4*PolyGamma[4, a] + 72072*PolyGamma[0, a]^5*
      PolyGamma[2, a]*PolyGamma[4, a] + 720720*PolyGamma[0, a]^3*
      PolyGamma[1, a]*PolyGamma[2, a]*PolyGamma[4, a] +
     1081080*PolyGamma[0, a]*PolyGamma[1, a]^2*PolyGamma[2, a]*
      PolyGamma[4, a] + 360360*PolyGamma[0, a]^2*PolyGamma[2, a]^2*
      PolyGamma[4, a] + 360360*PolyGamma[1, a]*PolyGamma[2, a]^2*
      PolyGamma[4, a] + 90090*PolyGamma[0, a]^4*PolyGamma[3, a]*
      PolyGamma[4, a] + 540540*PolyGamma[0, a]^2*PolyGamma[1, a]*
      PolyGamma[3, a]*PolyGamma[4, a] + 270270*PolyGamma[1, a]^2*
      PolyGamma[3, a]*PolyGamma[4, a] + 360360*PolyGamma[0, a]*
      PolyGamma[2, a]*PolyGamma[3, a]*PolyGamma[4, a] +
     45045*PolyGamma[3, a]^2*PolyGamma[4, a] + 36036*PolyGamma[0, a]^3*
      PolyGamma[4, a]^2 + 108108*PolyGamma[0, a]*PolyGamma[1, a]*
      PolyGamma[4, a]^2 + 36036*PolyGamma[2, a]*PolyGamma[4, a]^2 +
     1716*PolyGamma[0, a]^7*PolyGamma[5, a] + 36036*PolyGamma[0, a]^5*
      PolyGamma[1, a]*PolyGamma[5, a] + 180180*PolyGamma[0, a]^3*
      PolyGamma[1, a]^2*PolyGamma[5, a] + 180180*PolyGamma[0, a]*
      PolyGamma[1, a]^3*PolyGamma[5, a] + 60060*PolyGamma[0, a]^4*
      PolyGamma[2, a]*PolyGamma[5, a] + 360360*PolyGamma[0, a]^2*
      PolyGamma[1, a]*PolyGamma[2, a]*PolyGamma[5, a] +
     180180*PolyGamma[1, a]^2*PolyGamma[2, a]*PolyGamma[5, a] +
     120120*PolyGamma[0, a]*PolyGamma[2, a]^2*PolyGamma[5, a] +
     60060*PolyGamma[0, a]^3*PolyGamma[3, a]*PolyGamma[5, a] +
     180180*PolyGamma[0, a]*PolyGamma[1, a]*PolyGamma[3, a]*PolyGamma[5, a] +
     60060*PolyGamma[2, a]*PolyGamma[3, a]*PolyGamma[5, a] +
     36036*PolyGamma[0, a]^2*PolyGamma[4, a]*PolyGamma[5, a] +
     36036*PolyGamma[1, a]*PolyGamma[4, a]*PolyGamma[5, a] +
     6006*PolyGamma[0, a]*PolyGamma[5, a]^2 + 1716*PolyGamma[0, a]^6*
      PolyGamma[6, a] + 25740*PolyGamma[0, a]^4*PolyGamma[1, a]*
      PolyGamma[6, a] + 77220*PolyGamma[0, a]^2*PolyGamma[1, a]^2*
      PolyGamma[6, a] + 25740*PolyGamma[1, a]^3*PolyGamma[6, a] +
     34320*PolyGamma[0, a]^3*PolyGamma[2, a]*PolyGamma[6, a] +
     102960*PolyGamma[0, a]*PolyGamma[1, a]*PolyGamma[2, a]*PolyGamma[6, a] +
     17160*PolyGamma[2, a]^2*PolyGamma[6, a] + 25740*PolyGamma[0, a]^2*
      PolyGamma[3, a]*PolyGamma[6, a] + 25740*PolyGamma[1, a]*PolyGamma[3, a]*
      PolyGamma[6, a] + 10296*PolyGamma[0, a]*PolyGamma[4, a]*
      PolyGamma[6, a] + 1716*PolyGamma[5, a]*PolyGamma[6, a] +
     1287*PolyGamma[0, a]^5*PolyGamma[7, a] + 12870*PolyGamma[0, a]^3*
      PolyGamma[1, a]*PolyGamma[7, a] + 19305*PolyGamma[0, a]*
      PolyGamma[1, a]^2*PolyGamma[7, a] + 12870*PolyGamma[0, a]^2*
      PolyGamma[2, a]*PolyGamma[7, a] + 12870*PolyGamma[1, a]*PolyGamma[2, a]*
      PolyGamma[7, a] + 6435*PolyGamma[0, a]*PolyGamma[3, a]*
      PolyGamma[7, a] + 1287*PolyGamma[4, a]*PolyGamma[7, a] +
     715*PolyGamma[0, a]^4*PolyGamma[8, a] + 4290*PolyGamma[0, a]^2*
      PolyGamma[1, a]*PolyGamma[8, a] + 2145*PolyGamma[1, a]^2*
      PolyGamma[8, a] + 2860*PolyGamma[0, a]*PolyGamma[2, a]*
      PolyGamma[8, a] + 715*PolyGamma[3, a]*PolyGamma[8, a] +
     286*PolyGamma[0, a]^3*PolyGamma[9, a] + 858*PolyGamma[0, a]*
      PolyGamma[1, a]*PolyGamma[9, a] + 286*PolyGamma[2, a]*PolyGamma[9, a] +
     78*PolyGamma[0, a]^2*PolyGamma[10, a] + 78*PolyGamma[1, a]*
      PolyGamma[10, a] + 13*PolyGamma[0, a]*PolyGamma[11, a] +
     PolyGamma[12, a]

MBexpGam[a_, 14] = PolyGamma[0, a]^14 + 91*PolyGamma[0, a]^12*
      PolyGamma[1, a] + 3003*PolyGamma[0, a]^10*PolyGamma[1, a]^2 +
     45045*PolyGamma[0, a]^8*PolyGamma[1, a]^3 + 315315*PolyGamma[0, a]^6*
      PolyGamma[1, a]^4 + 945945*PolyGamma[0, a]^4*PolyGamma[1, a]^5 +
     945945*PolyGamma[0, a]^2*PolyGamma[1, a]^6 + 135135*PolyGamma[1, a]^7 +
     364*PolyGamma[0, a]^11*PolyGamma[2, a] + 20020*PolyGamma[0, a]^9*
      PolyGamma[1, a]*PolyGamma[2, a] + 360360*PolyGamma[0, a]^7*
      PolyGamma[1, a]^2*PolyGamma[2, a] + 2522520*PolyGamma[0, a]^5*
      PolyGamma[1, a]^3*PolyGamma[2, a] + 6306300*PolyGamma[0, a]^3*
      PolyGamma[1, a]^4*PolyGamma[2, a] + 3783780*PolyGamma[0, a]*
      PolyGamma[1, a]^5*PolyGamma[2, a] + 30030*PolyGamma[0, a]^8*
      PolyGamma[2, a]^2 + 840840*PolyGamma[0, a]^6*PolyGamma[1, a]*
      PolyGamma[2, a]^2 + 6306300*PolyGamma[0, a]^4*PolyGamma[1, a]^2*
      PolyGamma[2, a]^2 + 12612600*PolyGamma[0, a]^2*PolyGamma[1, a]^3*
      PolyGamma[2, a]^2 + 3153150*PolyGamma[1, a]^4*PolyGamma[2, a]^2 +
     560560*PolyGamma[0, a]^5*PolyGamma[2, a]^3 + 5605600*PolyGamma[0, a]^3*
      PolyGamma[1, a]*PolyGamma[2, a]^3 + 8408400*PolyGamma[0, a]*
      PolyGamma[1, a]^2*PolyGamma[2, a]^3 + 1401400*PolyGamma[0, a]^2*
      PolyGamma[2, a]^4 + 1401400*PolyGamma[1, a]*PolyGamma[2, a]^4 +
     1001*PolyGamma[0, a]^10*PolyGamma[3, a] + 45045*PolyGamma[0, a]^8*
      PolyGamma[1, a]*PolyGamma[3, a] + 630630*PolyGamma[0, a]^6*
      PolyGamma[1, a]^2*PolyGamma[3, a] + 3153150*PolyGamma[0, a]^4*
      PolyGamma[1, a]^3*PolyGamma[3, a] + 4729725*PolyGamma[0, a]^2*
      PolyGamma[1, a]^4*PolyGamma[3, a] + 945945*PolyGamma[1, a]^5*
      PolyGamma[3, a] + 120120*PolyGamma[0, a]^7*PolyGamma[2, a]*
      PolyGamma[3, a] + 2522520*PolyGamma[0, a]^5*PolyGamma[1, a]*
      PolyGamma[2, a]*PolyGamma[3, a] + 12612600*PolyGamma[0, a]^3*
      PolyGamma[1, a]^2*PolyGamma[2, a]*PolyGamma[3, a] +
     12612600*PolyGamma[0, a]*PolyGamma[1, a]^3*PolyGamma[2, a]*
      PolyGamma[3, a] + 2102100*PolyGamma[0, a]^4*PolyGamma[2, a]^2*
      PolyGamma[3, a] + 12612600*PolyGamma[0, a]^2*PolyGamma[1, a]*
      PolyGamma[2, a]^2*PolyGamma[3, a] + 6306300*PolyGamma[1, a]^2*
      PolyGamma[2, a]^2*PolyGamma[3, a] + 2802800*PolyGamma[0, a]*
      PolyGamma[2, a]^3*PolyGamma[3, a] + 105105*PolyGamma[0, a]^6*
      PolyGamma[3, a]^2 + 1576575*PolyGamma[0, a]^4*PolyGamma[1, a]*
      PolyGamma[3, a]^2 + 4729725*PolyGamma[0, a]^2*PolyGamma[1, a]^2*
      PolyGamma[3, a]^2 + 1576575*PolyGamma[1, a]^3*PolyGamma[3, a]^2 +
     2102100*PolyGamma[0, a]^3*PolyGamma[2, a]*PolyGamma[3, a]^2 +
     6306300*PolyGamma[0, a]*PolyGamma[1, a]*PolyGamma[2, a]*
      PolyGamma[3, a]^2 + 1051050*PolyGamma[2, a]^2*PolyGamma[3, a]^2 +
     525525*PolyGamma[0, a]^2*PolyGamma[3, a]^3 + 525525*PolyGamma[1, a]*
      PolyGamma[3, a]^3 + 2002*PolyGamma[0, a]^9*PolyGamma[4, a] +
     72072*PolyGamma[0, a]^7*PolyGamma[1, a]*PolyGamma[4, a] +
     756756*PolyGamma[0, a]^5*PolyGamma[1, a]^2*PolyGamma[4, a] +
     2522520*PolyGamma[0, a]^3*PolyGamma[1, a]^3*PolyGamma[4, a] +
     1891890*PolyGamma[0, a]*PolyGamma[1, a]^4*PolyGamma[4, a] +
     168168*PolyGamma[0, a]^6*PolyGamma[2, a]*PolyGamma[4, a] +
     2522520*PolyGamma[0, a]^4*PolyGamma[1, a]*PolyGamma[2, a]*
      PolyGamma[4, a] + 7567560*PolyGamma[0, a]^2*PolyGamma[1, a]^2*
      PolyGamma[2, a]*PolyGamma[4, a] + 2522520*PolyGamma[1, a]^3*
      PolyGamma[2, a]*PolyGamma[4, a] + 1681680*PolyGamma[0, a]^3*
      PolyGamma[2, a]^2*PolyGamma[4, a] + 5045040*PolyGamma[0, a]*
      PolyGamma[1, a]*PolyGamma[2, a]^2*PolyGamma[4, a] +
     560560*PolyGamma[2, a]^3*PolyGamma[4, a] + 252252*PolyGamma[0, a]^5*
      PolyGamma[3, a]*PolyGamma[4, a] + 2522520*PolyGamma[0, a]^3*
      PolyGamma[1, a]*PolyGamma[3, a]*PolyGamma[4, a] +
     3783780*PolyGamma[0, a]*PolyGamma[1, a]^2*PolyGamma[3, a]*
      PolyGamma[4, a] + 2522520*PolyGamma[0, a]^2*PolyGamma[2, a]*
      PolyGamma[3, a]*PolyGamma[4, a] + 2522520*PolyGamma[1, a]*
      PolyGamma[2, a]*PolyGamma[3, a]*PolyGamma[4, a] +
     630630*PolyGamma[0, a]*PolyGamma[3, a]^2*PolyGamma[4, a] +
     126126*PolyGamma[0, a]^4*PolyGamma[4, a]^2 + 756756*PolyGamma[0, a]^2*
      PolyGamma[1, a]*PolyGamma[4, a]^2 + 378378*PolyGamma[1, a]^2*
      PolyGamma[4, a]^2 + 504504*PolyGamma[0, a]*PolyGamma[2, a]*
      PolyGamma[4, a]^2 + 126126*PolyGamma[3, a]*PolyGamma[4, a]^2 +
     3003*PolyGamma[0, a]^8*PolyGamma[5, a] + 84084*PolyGamma[0, a]^6*
      PolyGamma[1, a]*PolyGamma[5, a] + 630630*PolyGamma[0, a]^4*
      PolyGamma[1, a]^2*PolyGamma[5, a] + 1261260*PolyGamma[0, a]^2*
      PolyGamma[1, a]^3*PolyGamma[5, a] + 315315*PolyGamma[1, a]^4*
      PolyGamma[5, a] + 168168*PolyGamma[0, a]^5*PolyGamma[2, a]*
      PolyGamma[5, a] + 1681680*PolyGamma[0, a]^3*PolyGamma[1, a]*
      PolyGamma[2, a]*PolyGamma[5, a] + 2522520*PolyGamma[0, a]*
      PolyGamma[1, a]^2*PolyGamma[2, a]*PolyGamma[5, a] +
     840840*PolyGamma[0, a]^2*PolyGamma[2, a]^2*PolyGamma[5, a] +
     840840*PolyGamma[1, a]*PolyGamma[2, a]^2*PolyGamma[5, a] +
     210210*PolyGamma[0, a]^4*PolyGamma[3, a]*PolyGamma[5, a] +
     1261260*PolyGamma[0, a]^2*PolyGamma[1, a]*PolyGamma[3, a]*
      PolyGamma[5, a] + 630630*PolyGamma[1, a]^2*PolyGamma[3, a]*
      PolyGamma[5, a] + 840840*PolyGamma[0, a]*PolyGamma[2, a]*
      PolyGamma[3, a]*PolyGamma[5, a] + 105105*PolyGamma[3, a]^2*
      PolyGamma[5, a] + 168168*PolyGamma[0, a]^3*PolyGamma[4, a]*
      PolyGamma[5, a] + 504504*PolyGamma[0, a]*PolyGamma[1, a]*
      PolyGamma[4, a]*PolyGamma[5, a] + 168168*PolyGamma[2, a]*
      PolyGamma[4, a]*PolyGamma[5, a] + 42042*PolyGamma[0, a]^2*
      PolyGamma[5, a]^2 + 42042*PolyGamma[1, a]*PolyGamma[5, a]^2 +
     3432*PolyGamma[0, a]^7*PolyGamma[6, a] + 72072*PolyGamma[0, a]^5*
      PolyGamma[1, a]*PolyGamma[6, a] + 360360*PolyGamma[0, a]^3*
      PolyGamma[1, a]^2*PolyGamma[6, a] + 360360*PolyGamma[0, a]*
      PolyGamma[1, a]^3*PolyGamma[6, a] + 120120*PolyGamma[0, a]^4*
      PolyGamma[2, a]*PolyGamma[6, a] + 720720*PolyGamma[0, a]^2*
      PolyGamma[1, a]*PolyGamma[2, a]*PolyGamma[6, a] +
     360360*PolyGamma[1, a]^2*PolyGamma[2, a]*PolyGamma[6, a] +
     240240*PolyGamma[0, a]*PolyGamma[2, a]^2*PolyGamma[6, a] +
     120120*PolyGamma[0, a]^3*PolyGamma[3, a]*PolyGamma[6, a] +
     360360*PolyGamma[0, a]*PolyGamma[1, a]*PolyGamma[3, a]*PolyGamma[6, a] +
     120120*PolyGamma[2, a]*PolyGamma[3, a]*PolyGamma[6, a] +
     72072*PolyGamma[0, a]^2*PolyGamma[4, a]*PolyGamma[6, a] +
     72072*PolyGamma[1, a]*PolyGamma[4, a]*PolyGamma[6, a] +
     24024*PolyGamma[0, a]*PolyGamma[5, a]*PolyGamma[6, a] +
     1716*PolyGamma[6, a]^2 + 3003*PolyGamma[0, a]^6*PolyGamma[7, a] +
     45045*PolyGamma[0, a]^4*PolyGamma[1, a]*PolyGamma[7, a] +
     135135*PolyGamma[0, a]^2*PolyGamma[1, a]^2*PolyGamma[7, a] +
     45045*PolyGamma[1, a]^3*PolyGamma[7, a] + 60060*PolyGamma[0, a]^3*
      PolyGamma[2, a]*PolyGamma[7, a] + 180180*PolyGamma[0, a]*
      PolyGamma[1, a]*PolyGamma[2, a]*PolyGamma[7, a] +
     30030*PolyGamma[2, a]^2*PolyGamma[7, a] + 45045*PolyGamma[0, a]^2*
      PolyGamma[3, a]*PolyGamma[7, a] + 45045*PolyGamma[1, a]*PolyGamma[3, a]*
      PolyGamma[7, a] + 18018*PolyGamma[0, a]*PolyGamma[4, a]*
      PolyGamma[7, a] + 3003*PolyGamma[5, a]*PolyGamma[7, a] +
     2002*PolyGamma[0, a]^5*PolyGamma[8, a] + 20020*PolyGamma[0, a]^3*
      PolyGamma[1, a]*PolyGamma[8, a] + 30030*PolyGamma[0, a]*
      PolyGamma[1, a]^2*PolyGamma[8, a] + 20020*PolyGamma[0, a]^2*
      PolyGamma[2, a]*PolyGamma[8, a] + 20020*PolyGamma[1, a]*PolyGamma[2, a]*
      PolyGamma[8, a] + 10010*PolyGamma[0, a]*PolyGamma[3, a]*
      PolyGamma[8, a] + 2002*PolyGamma[4, a]*PolyGamma[8, a] +
     1001*PolyGamma[0, a]^4*PolyGamma[9, a] + 6006*PolyGamma[0, a]^2*
      PolyGamma[1, a]*PolyGamma[9, a] + 3003*PolyGamma[1, a]^2*
      PolyGamma[9, a] + 4004*PolyGamma[0, a]*PolyGamma[2, a]*
      PolyGamma[9, a] + 1001*PolyGamma[3, a]*PolyGamma[9, a] +
     364*PolyGamma[0, a]^3*PolyGamma[10, a] + 1092*PolyGamma[0, a]*
      PolyGamma[1, a]*PolyGamma[10, a] + 364*PolyGamma[2, a]*
      PolyGamma[10, a] + 91*PolyGamma[0, a]^2*PolyGamma[11, a] +
     91*PolyGamma[1, a]*PolyGamma[11, a] + 14*PolyGamma[0, a]*
      PolyGamma[12, a] + PolyGamma[13, a]

MBexpGam[a_, 15] = PolyGamma[0, a]^15 + 105*PolyGamma[0, a]^13*
      PolyGamma[1, a] + 4095*PolyGamma[0, a]^11*PolyGamma[1, a]^2 +
     75075*PolyGamma[0, a]^9*PolyGamma[1, a]^3 + 675675*PolyGamma[0, a]^7*
      PolyGamma[1, a]^4 + 2837835*PolyGamma[0, a]^5*PolyGamma[1, a]^5 +
     4729725*PolyGamma[0, a]^3*PolyGamma[1, a]^6 + 2027025*PolyGamma[0, a]*
      PolyGamma[1, a]^7 + 455*PolyGamma[0, a]^12*PolyGamma[2, a] +
     30030*PolyGamma[0, a]^10*PolyGamma[1, a]*PolyGamma[2, a] +
     675675*PolyGamma[0, a]^8*PolyGamma[1, a]^2*PolyGamma[2, a] +
     6306300*PolyGamma[0, a]^6*PolyGamma[1, a]^3*PolyGamma[2, a] +
     23648625*PolyGamma[0, a]^4*PolyGamma[1, a]^4*PolyGamma[2, a] +
     28378350*PolyGamma[0, a]^2*PolyGamma[1, a]^5*PolyGamma[2, a] +
     4729725*PolyGamma[1, a]^6*PolyGamma[2, a] + 50050*PolyGamma[0, a]^9*
      PolyGamma[2, a]^2 + 1801800*PolyGamma[0, a]^7*PolyGamma[1, a]*
      PolyGamma[2, a]^2 + 18918900*PolyGamma[0, a]^5*PolyGamma[1, a]^2*
      PolyGamma[2, a]^2 + 63063000*PolyGamma[0, a]^3*PolyGamma[1, a]^3*
      PolyGamma[2, a]^2 + 47297250*PolyGamma[0, a]*PolyGamma[1, a]^4*
      PolyGamma[2, a]^2 + 1401400*PolyGamma[0, a]^6*PolyGamma[2, a]^3 +
     21021000*PolyGamma[0, a]^4*PolyGamma[1, a]*PolyGamma[2, a]^3 +
     63063000*PolyGamma[0, a]^2*PolyGamma[1, a]^2*PolyGamma[2, a]^3 +
     21021000*PolyGamma[1, a]^3*PolyGamma[2, a]^3 +
     7007000*PolyGamma[0, a]^3*PolyGamma[2, a]^4 + 21021000*PolyGamma[0, a]*
      PolyGamma[1, a]*PolyGamma[2, a]^4 + 1401400*PolyGamma[2, a]^5 +
     1365*PolyGamma[0, a]^11*PolyGamma[3, a] + 75075*PolyGamma[0, a]^9*
      PolyGamma[1, a]*PolyGamma[3, a] + 1351350*PolyGamma[0, a]^7*
      PolyGamma[1, a]^2*PolyGamma[3, a] + 9459450*PolyGamma[0, a]^5*
      PolyGamma[1, a]^3*PolyGamma[3, a] + 23648625*PolyGamma[0, a]^3*
      PolyGamma[1, a]^4*PolyGamma[3, a] + 14189175*PolyGamma[0, a]*
      PolyGamma[1, a]^5*PolyGamma[3, a] + 225225*PolyGamma[0, a]^8*
      PolyGamma[2, a]*PolyGamma[3, a] + 6306300*PolyGamma[0, a]^6*
      PolyGamma[1, a]*PolyGamma[2, a]*PolyGamma[3, a] +
     47297250*PolyGamma[0, a]^4*PolyGamma[1, a]^2*PolyGamma[2, a]*
      PolyGamma[3, a] + 94594500*PolyGamma[0, a]^2*PolyGamma[1, a]^3*
      PolyGamma[2, a]*PolyGamma[3, a] + 23648625*PolyGamma[1, a]^4*
      PolyGamma[2, a]*PolyGamma[3, a] + 6306300*PolyGamma[0, a]^5*
      PolyGamma[2, a]^2*PolyGamma[3, a] + 63063000*PolyGamma[0, a]^3*
      PolyGamma[1, a]*PolyGamma[2, a]^2*PolyGamma[3, a] +
     94594500*PolyGamma[0, a]*PolyGamma[1, a]^2*PolyGamma[2, a]^2*
      PolyGamma[3, a] + 21021000*PolyGamma[0, a]^2*PolyGamma[2, a]^3*
      PolyGamma[3, a] + 21021000*PolyGamma[1, a]*PolyGamma[2, a]^3*
      PolyGamma[3, a] + 225225*PolyGamma[0, a]^7*PolyGamma[3, a]^2 +
     4729725*PolyGamma[0, a]^5*PolyGamma[1, a]*PolyGamma[3, a]^2 +
     23648625*PolyGamma[0, a]^3*PolyGamma[1, a]^2*PolyGamma[3, a]^2 +
     23648625*PolyGamma[0, a]*PolyGamma[1, a]^3*PolyGamma[3, a]^2 +
     7882875*PolyGamma[0, a]^4*PolyGamma[2, a]*PolyGamma[3, a]^2 +
     47297250*PolyGamma[0, a]^2*PolyGamma[1, a]*PolyGamma[2, a]*
      PolyGamma[3, a]^2 + 23648625*PolyGamma[1, a]^2*PolyGamma[2, a]*
      PolyGamma[3, a]^2 + 15765750*PolyGamma[0, a]*PolyGamma[2, a]^2*
      PolyGamma[3, a]^2 + 2627625*PolyGamma[0, a]^3*PolyGamma[3, a]^3 +
     7882875*PolyGamma[0, a]*PolyGamma[1, a]*PolyGamma[3, a]^3 +
     2627625*PolyGamma[2, a]*PolyGamma[3, a]^3 + 3003*PolyGamma[0, a]^10*
      PolyGamma[4, a] + 135135*PolyGamma[0, a]^8*PolyGamma[1, a]*
      PolyGamma[4, a] + 1891890*PolyGamma[0, a]^6*PolyGamma[1, a]^2*
      PolyGamma[4, a] + 9459450*PolyGamma[0, a]^4*PolyGamma[1, a]^3*
      PolyGamma[4, a] + 14189175*PolyGamma[0, a]^2*PolyGamma[1, a]^4*
      PolyGamma[4, a] + 2837835*PolyGamma[1, a]^5*PolyGamma[4, a] +
     360360*PolyGamma[0, a]^7*PolyGamma[2, a]*PolyGamma[4, a] +
     7567560*PolyGamma[0, a]^5*PolyGamma[1, a]*PolyGamma[2, a]*
      PolyGamma[4, a] + 37837800*PolyGamma[0, a]^3*PolyGamma[1, a]^2*
      PolyGamma[2, a]*PolyGamma[4, a] + 37837800*PolyGamma[0, a]*
      PolyGamma[1, a]^3*PolyGamma[2, a]*PolyGamma[4, a] +
     6306300*PolyGamma[0, a]^4*PolyGamma[2, a]^2*PolyGamma[4, a] +
     37837800*PolyGamma[0, a]^2*PolyGamma[1, a]*PolyGamma[2, a]^2*
      PolyGamma[4, a] + 18918900*PolyGamma[1, a]^2*PolyGamma[2, a]^2*
      PolyGamma[4, a] + 8408400*PolyGamma[0, a]*PolyGamma[2, a]^3*
      PolyGamma[4, a] + 630630*PolyGamma[0, a]^6*PolyGamma[3, a]*
      PolyGamma[4, a] + 9459450*PolyGamma[0, a]^4*PolyGamma[1, a]*
      PolyGamma[3, a]*PolyGamma[4, a] + 28378350*PolyGamma[0, a]^2*
      PolyGamma[1, a]^2*PolyGamma[3, a]*PolyGamma[4, a] +
     9459450*PolyGamma[1, a]^3*PolyGamma[3, a]*PolyGamma[4, a] +
     12612600*PolyGamma[0, a]^3*PolyGamma[2, a]*PolyGamma[3, a]*
      PolyGamma[4, a] + 37837800*PolyGamma[0, a]*PolyGamma[1, a]*
      PolyGamma[2, a]*PolyGamma[3, a]*PolyGamma[4, a] +
     6306300*PolyGamma[2, a]^2*PolyGamma[3, a]*PolyGamma[4, a] +
     4729725*PolyGamma[0, a]^2*PolyGamma[3, a]^2*PolyGamma[4, a] +
     4729725*PolyGamma[1, a]*PolyGamma[3, a]^2*PolyGamma[4, a] +
     378378*PolyGamma[0, a]^5*PolyGamma[4, a]^2 + 3783780*PolyGamma[0, a]^3*
      PolyGamma[1, a]*PolyGamma[4, a]^2 + 5675670*PolyGamma[0, a]*
      PolyGamma[1, a]^2*PolyGamma[4, a]^2 + 3783780*PolyGamma[0, a]^2*
      PolyGamma[2, a]*PolyGamma[4, a]^2 + 3783780*PolyGamma[1, a]*
      PolyGamma[2, a]*PolyGamma[4, a]^2 + 1891890*PolyGamma[0, a]*
      PolyGamma[3, a]*PolyGamma[4, a]^2 + 126126*PolyGamma[4, a]^3 +
     5005*PolyGamma[0, a]^9*PolyGamma[5, a] + 180180*PolyGamma[0, a]^7*
      PolyGamma[1, a]*PolyGamma[5, a] + 1891890*PolyGamma[0, a]^5*
      PolyGamma[1, a]^2*PolyGamma[5, a] + 6306300*PolyGamma[0, a]^3*
      PolyGamma[1, a]^3*PolyGamma[5, a] + 4729725*PolyGamma[0, a]*
      PolyGamma[1, a]^4*PolyGamma[5, a] + 420420*PolyGamma[0, a]^6*
      PolyGamma[2, a]*PolyGamma[5, a] + 6306300*PolyGamma[0, a]^4*
      PolyGamma[1, a]*PolyGamma[2, a]*PolyGamma[5, a] +
     18918900*PolyGamma[0, a]^2*PolyGamma[1, a]^2*PolyGamma[2, a]*
      PolyGamma[5, a] + 6306300*PolyGamma[1, a]^3*PolyGamma[2, a]*
      PolyGamma[5, a] + 4204200*PolyGamma[0, a]^3*PolyGamma[2, a]^2*
      PolyGamma[5, a] + 12612600*PolyGamma[0, a]*PolyGamma[1, a]*
      PolyGamma[2, a]^2*PolyGamma[5, a] + 1401400*PolyGamma[2, a]^3*
      PolyGamma[5, a] + 630630*PolyGamma[0, a]^5*PolyGamma[3, a]*
      PolyGamma[5, a] + 6306300*PolyGamma[0, a]^3*PolyGamma[1, a]*
      PolyGamma[3, a]*PolyGamma[5, a] + 9459450*PolyGamma[0, a]*
      PolyGamma[1, a]^2*PolyGamma[3, a]*PolyGamma[5, a] +
     6306300*PolyGamma[0, a]^2*PolyGamma[2, a]*PolyGamma[3, a]*
      PolyGamma[5, a] + 6306300*PolyGamma[1, a]*PolyGamma[2, a]*
      PolyGamma[3, a]*PolyGamma[5, a] + 1576575*PolyGamma[0, a]*
      PolyGamma[3, a]^2*PolyGamma[5, a] + 630630*PolyGamma[0, a]^4*
      PolyGamma[4, a]*PolyGamma[5, a] + 3783780*PolyGamma[0, a]^2*
      PolyGamma[1, a]*PolyGamma[4, a]*PolyGamma[5, a] +
     1891890*PolyGamma[1, a]^2*PolyGamma[4, a]*PolyGamma[5, a] +
     2522520*PolyGamma[0, a]*PolyGamma[2, a]*PolyGamma[4, a]*
      PolyGamma[5, a] + 630630*PolyGamma[3, a]*PolyGamma[4, a]*
      PolyGamma[5, a] + 210210*PolyGamma[0, a]^3*PolyGamma[5, a]^2 +
     630630*PolyGamma[0, a]*PolyGamma[1, a]*PolyGamma[5, a]^2 +
     210210*PolyGamma[2, a]*PolyGamma[5, a]^2 + 6435*PolyGamma[0, a]^8*
      PolyGamma[6, a] + 180180*PolyGamma[0, a]^6*PolyGamma[1, a]*
      PolyGamma[6, a] + 1351350*PolyGamma[0, a]^4*PolyGamma[1, a]^2*
      PolyGamma[6, a] + 2702700*PolyGamma[0, a]^2*PolyGamma[1, a]^3*
      PolyGamma[6, a] + 675675*PolyGamma[1, a]^4*PolyGamma[6, a] +
     360360*PolyGamma[0, a]^5*PolyGamma[2, a]*PolyGamma[6, a] +
     3603600*PolyGamma[0, a]^3*PolyGamma[1, a]*PolyGamma[2, a]*
      PolyGamma[6, a] + 5405400*PolyGamma[0, a]*PolyGamma[1, a]^2*
      PolyGamma[2, a]*PolyGamma[6, a] + 1801800*PolyGamma[0, a]^2*
      PolyGamma[2, a]^2*PolyGamma[6, a] + 1801800*PolyGamma[1, a]*
      PolyGamma[2, a]^2*PolyGamma[6, a] + 450450*PolyGamma[0, a]^4*
      PolyGamma[3, a]*PolyGamma[6, a] + 2702700*PolyGamma[0, a]^2*
      PolyGamma[1, a]*PolyGamma[3, a]*PolyGamma[6, a] +
     1351350*PolyGamma[1, a]^2*PolyGamma[3, a]*PolyGamma[6, a] +
     1801800*PolyGamma[0, a]*PolyGamma[2, a]*PolyGamma[3, a]*
      PolyGamma[6, a] + 225225*PolyGamma[3, a]^2*PolyGamma[6, a] +
     360360*PolyGamma[0, a]^3*PolyGamma[4, a]*PolyGamma[6, a] +
     1081080*PolyGamma[0, a]*PolyGamma[1, a]*PolyGamma[4, a]*
      PolyGamma[6, a] + 360360*PolyGamma[2, a]*PolyGamma[4, a]*
      PolyGamma[6, a] + 180180*PolyGamma[0, a]^2*PolyGamma[5, a]*
      PolyGamma[6, a] + 180180*PolyGamma[1, a]*PolyGamma[5, a]*
      PolyGamma[6, a] + 25740*PolyGamma[0, a]*PolyGamma[6, a]^2 +
     6435*PolyGamma[0, a]^7*PolyGamma[7, a] + 135135*PolyGamma[0, a]^5*
      PolyGamma[1, a]*PolyGamma[7, a] + 675675*PolyGamma[0, a]^3*
      PolyGamma[1, a]^2*PolyGamma[7, a] + 675675*PolyGamma[0, a]*
      PolyGamma[1, a]^3*PolyGamma[7, a] + 225225*PolyGamma[0, a]^4*
      PolyGamma[2, a]*PolyGamma[7, a] + 1351350*PolyGamma[0, a]^2*
      PolyGamma[1, a]*PolyGamma[2, a]*PolyGamma[7, a] +
     675675*PolyGamma[1, a]^2*PolyGamma[2, a]*PolyGamma[7, a] +
     450450*PolyGamma[0, a]*PolyGamma[2, a]^2*PolyGamma[7, a] +
     225225*PolyGamma[0, a]^3*PolyGamma[3, a]*PolyGamma[7, a] +
     675675*PolyGamma[0, a]*PolyGamma[1, a]*PolyGamma[3, a]*PolyGamma[7, a] +
     225225*PolyGamma[2, a]*PolyGamma[3, a]*PolyGamma[7, a] +
     135135*PolyGamma[0, a]^2*PolyGamma[4, a]*PolyGamma[7, a] +
     135135*PolyGamma[1, a]*PolyGamma[4, a]*PolyGamma[7, a] +
     45045*PolyGamma[0, a]*PolyGamma[5, a]*PolyGamma[7, a] +
     6435*PolyGamma[6, a]*PolyGamma[7, a] + 5005*PolyGamma[0, a]^6*
      PolyGamma[8, a] + 75075*PolyGamma[0, a]^4*PolyGamma[1, a]*
      PolyGamma[8, a] + 225225*PolyGamma[0, a]^2*PolyGamma[1, a]^2*
      PolyGamma[8, a] + 75075*PolyGamma[1, a]^3*PolyGamma[8, a] +
     100100*PolyGamma[0, a]^3*PolyGamma[2, a]*PolyGamma[8, a] +
     300300*PolyGamma[0, a]*PolyGamma[1, a]*PolyGamma[2, a]*PolyGamma[8, a] +
     50050*PolyGamma[2, a]^2*PolyGamma[8, a] + 75075*PolyGamma[0, a]^2*
      PolyGamma[3, a]*PolyGamma[8, a] + 75075*PolyGamma[1, a]*PolyGamma[3, a]*
      PolyGamma[8, a] + 30030*PolyGamma[0, a]*PolyGamma[4, a]*
      PolyGamma[8, a] + 5005*PolyGamma[5, a]*PolyGamma[8, a] +
     3003*PolyGamma[0, a]^5*PolyGamma[9, a] + 30030*PolyGamma[0, a]^3*
      PolyGamma[1, a]*PolyGamma[9, a] + 45045*PolyGamma[0, a]*
      PolyGamma[1, a]^2*PolyGamma[9, a] + 30030*PolyGamma[0, a]^2*
      PolyGamma[2, a]*PolyGamma[9, a] + 30030*PolyGamma[1, a]*PolyGamma[2, a]*
      PolyGamma[9, a] + 15015*PolyGamma[0, a]*PolyGamma[3, a]*
      PolyGamma[9, a] + 3003*PolyGamma[4, a]*PolyGamma[9, a] +
     1365*PolyGamma[0, a]^4*PolyGamma[10, a] + 8190*PolyGamma[0, a]^2*
      PolyGamma[1, a]*PolyGamma[10, a] + 4095*PolyGamma[1, a]^2*
      PolyGamma[10, a] + 5460*PolyGamma[0, a]*PolyGamma[2, a]*
      PolyGamma[10, a] + 1365*PolyGamma[3, a]*PolyGamma[10, a] +
     455*PolyGamma[0, a]^3*PolyGamma[11, a] + 1365*PolyGamma[0, a]*
      PolyGamma[1, a]*PolyGamma[11, a] + 455*PolyGamma[2, a]*
      PolyGamma[11, a] + 105*PolyGamma[0, a]^2*PolyGamma[12, a] +
     105*PolyGamma[1, a]*PolyGamma[12, a] + 15*PolyGamma[0, a]*
      PolyGamma[13, a] + PolyGamma[14, a]

MBexpGam[a_, 16] = PolyGamma[0, a]^16 + 120*PolyGamma[0, a]^14*
      PolyGamma[1, a] + 5460*PolyGamma[0, a]^12*PolyGamma[1, a]^2 +
     120120*PolyGamma[0, a]^10*PolyGamma[1, a]^3 + 1351350*PolyGamma[0, a]^8*
      PolyGamma[1, a]^4 + 7567560*PolyGamma[0, a]^6*PolyGamma[1, a]^5 +
     18918900*PolyGamma[0, a]^4*PolyGamma[1, a]^6 +
     16216200*PolyGamma[0, a]^2*PolyGamma[1, a]^7 +
     2027025*PolyGamma[1, a]^8 + 560*PolyGamma[0, a]^13*PolyGamma[2, a] +
     43680*PolyGamma[0, a]^11*PolyGamma[1, a]*PolyGamma[2, a] +
     1201200*PolyGamma[0, a]^9*PolyGamma[1, a]^2*PolyGamma[2, a] +
     14414400*PolyGamma[0, a]^7*PolyGamma[1, a]^3*PolyGamma[2, a] +
     75675600*PolyGamma[0, a]^5*PolyGamma[1, a]^4*PolyGamma[2, a] +
     151351200*PolyGamma[0, a]^3*PolyGamma[1, a]^5*PolyGamma[2, a] +
     75675600*PolyGamma[0, a]*PolyGamma[1, a]^6*PolyGamma[2, a] +
     80080*PolyGamma[0, a]^10*PolyGamma[2, a]^2 + 3603600*PolyGamma[0, a]^8*
      PolyGamma[1, a]*PolyGamma[2, a]^2 + 50450400*PolyGamma[0, a]^6*
      PolyGamma[1, a]^2*PolyGamma[2, a]^2 + 252252000*PolyGamma[0, a]^4*
      PolyGamma[1, a]^3*PolyGamma[2, a]^2 + 378378000*PolyGamma[0, a]^2*
      PolyGamma[1, a]^4*PolyGamma[2, a]^2 + 75675600*PolyGamma[1, a]^5*
      PolyGamma[2, a]^2 + 3203200*PolyGamma[0, a]^7*PolyGamma[2, a]^3 +
     67267200*PolyGamma[0, a]^5*PolyGamma[1, a]*PolyGamma[2, a]^3 +
     336336000*PolyGamma[0, a]^3*PolyGamma[1, a]^2*PolyGamma[2, a]^3 +
     336336000*PolyGamma[0, a]*PolyGamma[1, a]^3*PolyGamma[2, a]^3 +
     28028000*PolyGamma[0, a]^4*PolyGamma[2, a]^4 +
     168168000*PolyGamma[0, a]^2*PolyGamma[1, a]*PolyGamma[2, a]^4 +
     84084000*PolyGamma[1, a]^2*PolyGamma[2, a]^4 +
     22422400*PolyGamma[0, a]*PolyGamma[2, a]^5 + 1820*PolyGamma[0, a]^12*
      PolyGamma[3, a] + 120120*PolyGamma[0, a]^10*PolyGamma[1, a]*
      PolyGamma[3, a] + 2702700*PolyGamma[0, a]^8*PolyGamma[1, a]^2*
      PolyGamma[3, a] + 25225200*PolyGamma[0, a]^6*PolyGamma[1, a]^3*
      PolyGamma[3, a] + 94594500*PolyGamma[0, a]^4*PolyGamma[1, a]^4*
      PolyGamma[3, a] + 113513400*PolyGamma[0, a]^2*PolyGamma[1, a]^5*
      PolyGamma[3, a] + 18918900*PolyGamma[1, a]^6*PolyGamma[3, a] +
     400400*PolyGamma[0, a]^9*PolyGamma[2, a]*PolyGamma[3, a] +
     14414400*PolyGamma[0, a]^7*PolyGamma[1, a]*PolyGamma[2, a]*
      PolyGamma[3, a] + 151351200*PolyGamma[0, a]^5*PolyGamma[1, a]^2*
      PolyGamma[2, a]*PolyGamma[3, a] + 504504000*PolyGamma[0, a]^3*
      PolyGamma[1, a]^3*PolyGamma[2, a]*PolyGamma[3, a] +
     378378000*PolyGamma[0, a]*PolyGamma[1, a]^4*PolyGamma[2, a]*
      PolyGamma[3, a] + 16816800*PolyGamma[0, a]^6*PolyGamma[2, a]^2*
      PolyGamma[3, a] + 252252000*PolyGamma[0, a]^4*PolyGamma[1, a]*
      PolyGamma[2, a]^2*PolyGamma[3, a] + 756756000*PolyGamma[0, a]^2*
      PolyGamma[1, a]^2*PolyGamma[2, a]^2*PolyGamma[3, a] +
     252252000*PolyGamma[1, a]^3*PolyGamma[2, a]^2*PolyGamma[3, a] +
     112112000*PolyGamma[0, a]^3*PolyGamma[2, a]^3*PolyGamma[3, a] +
     336336000*PolyGamma[0, a]*PolyGamma[1, a]*PolyGamma[2, a]^3*
      PolyGamma[3, a] + 28028000*PolyGamma[2, a]^4*PolyGamma[3, a] +
     450450*PolyGamma[0, a]^8*PolyGamma[3, a]^2 + 12612600*PolyGamma[0, a]^6*
      PolyGamma[1, a]*PolyGamma[3, a]^2 + 94594500*PolyGamma[0, a]^4*
      PolyGamma[1, a]^2*PolyGamma[3, a]^2 + 189189000*PolyGamma[0, a]^2*
      PolyGamma[1, a]^3*PolyGamma[3, a]^2 + 47297250*PolyGamma[1, a]^4*
      PolyGamma[3, a]^2 + 25225200*PolyGamma[0, a]^5*PolyGamma[2, a]*
      PolyGamma[3, a]^2 + 252252000*PolyGamma[0, a]^3*PolyGamma[1, a]*
      PolyGamma[2, a]*PolyGamma[3, a]^2 + 378378000*PolyGamma[0, a]*
      PolyGamma[1, a]^2*PolyGamma[2, a]*PolyGamma[3, a]^2 +
     126126000*PolyGamma[0, a]^2*PolyGamma[2, a]^2*PolyGamma[3, a]^2 +
     126126000*PolyGamma[1, a]*PolyGamma[2, a]^2*PolyGamma[3, a]^2 +
     10510500*PolyGamma[0, a]^4*PolyGamma[3, a]^3 +
     63063000*PolyGamma[0, a]^2*PolyGamma[1, a]*PolyGamma[3, a]^3 +
     31531500*PolyGamma[1, a]^2*PolyGamma[3, a]^3 +
     42042000*PolyGamma[0, a]*PolyGamma[2, a]*PolyGamma[3, a]^3 +
     2627625*PolyGamma[3, a]^4 + 4368*PolyGamma[0, a]^11*PolyGamma[4, a] +
     240240*PolyGamma[0, a]^9*PolyGamma[1, a]*PolyGamma[4, a] +
     4324320*PolyGamma[0, a]^7*PolyGamma[1, a]^2*PolyGamma[4, a] +
     30270240*PolyGamma[0, a]^5*PolyGamma[1, a]^3*PolyGamma[4, a] +
     75675600*PolyGamma[0, a]^3*PolyGamma[1, a]^4*PolyGamma[4, a] +
     45405360*PolyGamma[0, a]*PolyGamma[1, a]^5*PolyGamma[4, a] +
     720720*PolyGamma[0, a]^8*PolyGamma[2, a]*PolyGamma[4, a] +
     20180160*PolyGamma[0, a]^6*PolyGamma[1, a]*PolyGamma[2, a]*
      PolyGamma[4, a] + 151351200*PolyGamma[0, a]^4*PolyGamma[1, a]^2*
      PolyGamma[2, a]*PolyGamma[4, a] + 302702400*PolyGamma[0, a]^2*
      PolyGamma[1, a]^3*PolyGamma[2, a]*PolyGamma[4, a] +
     75675600*PolyGamma[1, a]^4*PolyGamma[2, a]*PolyGamma[4, a] +
     20180160*PolyGamma[0, a]^5*PolyGamma[2, a]^2*PolyGamma[4, a] +
     201801600*PolyGamma[0, a]^3*PolyGamma[1, a]*PolyGamma[2, a]^2*
      PolyGamma[4, a] + 302702400*PolyGamma[0, a]*PolyGamma[1, a]^2*
      PolyGamma[2, a]^2*PolyGamma[4, a] + 67267200*PolyGamma[0, a]^2*
      PolyGamma[2, a]^3*PolyGamma[4, a] + 67267200*PolyGamma[1, a]*
      PolyGamma[2, a]^3*PolyGamma[4, a] + 1441440*PolyGamma[0, a]^7*
      PolyGamma[3, a]*PolyGamma[4, a] + 30270240*PolyGamma[0, a]^5*
      PolyGamma[1, a]*PolyGamma[3, a]*PolyGamma[4, a] +
     151351200*PolyGamma[0, a]^3*PolyGamma[1, a]^2*PolyGamma[3, a]*
      PolyGamma[4, a] + 151351200*PolyGamma[0, a]*PolyGamma[1, a]^3*
      PolyGamma[3, a]*PolyGamma[4, a] + 50450400*PolyGamma[0, a]^4*
      PolyGamma[2, a]*PolyGamma[3, a]*PolyGamma[4, a] +
     302702400*PolyGamma[0, a]^2*PolyGamma[1, a]*PolyGamma[2, a]*
      PolyGamma[3, a]*PolyGamma[4, a] + 151351200*PolyGamma[1, a]^2*
      PolyGamma[2, a]*PolyGamma[3, a]*PolyGamma[4, a] +
     100900800*PolyGamma[0, a]*PolyGamma[2, a]^2*PolyGamma[3, a]*
      PolyGamma[4, a] + 25225200*PolyGamma[0, a]^3*PolyGamma[3, a]^2*
      PolyGamma[4, a] + 75675600*PolyGamma[0, a]*PolyGamma[1, a]*
      PolyGamma[3, a]^2*PolyGamma[4, a] + 25225200*PolyGamma[2, a]*
      PolyGamma[3, a]^2*PolyGamma[4, a] + 1009008*PolyGamma[0, a]^6*
      PolyGamma[4, a]^2 + 15135120*PolyGamma[0, a]^4*PolyGamma[1, a]*
      PolyGamma[4, a]^2 + 45405360*PolyGamma[0, a]^2*PolyGamma[1, a]^2*
      PolyGamma[4, a]^2 + 15135120*PolyGamma[1, a]^3*PolyGamma[4, a]^2 +
     20180160*PolyGamma[0, a]^3*PolyGamma[2, a]*PolyGamma[4, a]^2 +
     60540480*PolyGamma[0, a]*PolyGamma[1, a]*PolyGamma[2, a]*
      PolyGamma[4, a]^2 + 10090080*PolyGamma[2, a]^2*PolyGamma[4, a]^2 +
     15135120*PolyGamma[0, a]^2*PolyGamma[3, a]*PolyGamma[4, a]^2 +
     15135120*PolyGamma[1, a]*PolyGamma[3, a]*PolyGamma[4, a]^2 +
     2018016*PolyGamma[0, a]*PolyGamma[4, a]^3 + 8008*PolyGamma[0, a]^10*
      PolyGamma[5, a] + 360360*PolyGamma[0, a]^8*PolyGamma[1, a]*
      PolyGamma[5, a] + 5045040*PolyGamma[0, a]^6*PolyGamma[1, a]^2*
      PolyGamma[5, a] + 25225200*PolyGamma[0, a]^4*PolyGamma[1, a]^3*
      PolyGamma[5, a] + 37837800*PolyGamma[0, a]^2*PolyGamma[1, a]^4*
      PolyGamma[5, a] + 7567560*PolyGamma[1, a]^5*PolyGamma[5, a] +
     960960*PolyGamma[0, a]^7*PolyGamma[2, a]*PolyGamma[5, a] +
     20180160*PolyGamma[0, a]^5*PolyGamma[1, a]*PolyGamma[2, a]*
      PolyGamma[5, a] + 100900800*PolyGamma[0, a]^3*PolyGamma[1, a]^2*
      PolyGamma[2, a]*PolyGamma[5, a] + 100900800*PolyGamma[0, a]*
      PolyGamma[1, a]^3*PolyGamma[2, a]*PolyGamma[5, a] +
     16816800*PolyGamma[0, a]^4*PolyGamma[2, a]^2*PolyGamma[5, a] +
     100900800*PolyGamma[0, a]^2*PolyGamma[1, a]*PolyGamma[2, a]^2*
      PolyGamma[5, a] + 50450400*PolyGamma[1, a]^2*PolyGamma[2, a]^2*
      PolyGamma[5, a] + 22422400*PolyGamma[0, a]*PolyGamma[2, a]^3*
      PolyGamma[5, a] + 1681680*PolyGamma[0, a]^6*PolyGamma[3, a]*
      PolyGamma[5, a] + 25225200*PolyGamma[0, a]^4*PolyGamma[1, a]*
      PolyGamma[3, a]*PolyGamma[5, a] + 75675600*PolyGamma[0, a]^2*
      PolyGamma[1, a]^2*PolyGamma[3, a]*PolyGamma[5, a] +
     25225200*PolyGamma[1, a]^3*PolyGamma[3, a]*PolyGamma[5, a] +
     33633600*PolyGamma[0, a]^3*PolyGamma[2, a]*PolyGamma[3, a]*
      PolyGamma[5, a] + 100900800*PolyGamma[0, a]*PolyGamma[1, a]*
      PolyGamma[2, a]*PolyGamma[3, a]*PolyGamma[5, a] +
     16816800*PolyGamma[2, a]^2*PolyGamma[3, a]*PolyGamma[5, a] +
     12612600*PolyGamma[0, a]^2*PolyGamma[3, a]^2*PolyGamma[5, a] +
     12612600*PolyGamma[1, a]*PolyGamma[3, a]^2*PolyGamma[5, a] +
     2018016*PolyGamma[0, a]^5*PolyGamma[4, a]*PolyGamma[5, a] +
     20180160*PolyGamma[0, a]^3*PolyGamma[1, a]*PolyGamma[4, a]*
      PolyGamma[5, a] + 30270240*PolyGamma[0, a]*PolyGamma[1, a]^2*
      PolyGamma[4, a]*PolyGamma[5, a] + 20180160*PolyGamma[0, a]^2*
      PolyGamma[2, a]*PolyGamma[4, a]*PolyGamma[5, a] +
     20180160*PolyGamma[1, a]*PolyGamma[2, a]*PolyGamma[4, a]*
      PolyGamma[5, a] + 10090080*PolyGamma[0, a]*PolyGamma[3, a]*
      PolyGamma[4, a]*PolyGamma[5, a] + 1009008*PolyGamma[4, a]^2*
      PolyGamma[5, a] + 840840*PolyGamma[0, a]^4*PolyGamma[5, a]^2 +
     5045040*PolyGamma[0, a]^2*PolyGamma[1, a]*PolyGamma[5, a]^2 +
     2522520*PolyGamma[1, a]^2*PolyGamma[5, a]^2 + 3363360*PolyGamma[0, a]*
      PolyGamma[2, a]*PolyGamma[5, a]^2 + 840840*PolyGamma[3, a]*
      PolyGamma[5, a]^2 + 11440*PolyGamma[0, a]^9*PolyGamma[6, a] +
     411840*PolyGamma[0, a]^7*PolyGamma[1, a]*PolyGamma[6, a] +
     4324320*PolyGamma[0, a]^5*PolyGamma[1, a]^2*PolyGamma[6, a] +
     14414400*PolyGamma[0, a]^3*PolyGamma[1, a]^3*PolyGamma[6, a] +
     10810800*PolyGamma[0, a]*PolyGamma[1, a]^4*PolyGamma[6, a] +
     960960*PolyGamma[0, a]^6*PolyGamma[2, a]*PolyGamma[6, a] +
     14414400*PolyGamma[0, a]^4*PolyGamma[1, a]*PolyGamma[2, a]*
      PolyGamma[6, a] + 43243200*PolyGamma[0, a]^2*PolyGamma[1, a]^2*
      PolyGamma[2, a]*PolyGamma[6, a] + 14414400*PolyGamma[1, a]^3*
      PolyGamma[2, a]*PolyGamma[6, a] + 9609600*PolyGamma[0, a]^3*
      PolyGamma[2, a]^2*PolyGamma[6, a] + 28828800*PolyGamma[0, a]*
      PolyGamma[1, a]*PolyGamma[2, a]^2*PolyGamma[6, a] +
     3203200*PolyGamma[2, a]^3*PolyGamma[6, a] + 1441440*PolyGamma[0, a]^5*
      PolyGamma[3, a]*PolyGamma[6, a] + 14414400*PolyGamma[0, a]^3*
      PolyGamma[1, a]*PolyGamma[3, a]*PolyGamma[6, a] +
     21621600*PolyGamma[0, a]*PolyGamma[1, a]^2*PolyGamma[3, a]*
      PolyGamma[6, a] + 14414400*PolyGamma[0, a]^2*PolyGamma[2, a]*
      PolyGamma[3, a]*PolyGamma[6, a] + 14414400*PolyGamma[1, a]*
      PolyGamma[2, a]*PolyGamma[3, a]*PolyGamma[6, a] +
     3603600*PolyGamma[0, a]*PolyGamma[3, a]^2*PolyGamma[6, a] +
     1441440*PolyGamma[0, a]^4*PolyGamma[4, a]*PolyGamma[6, a] +
     8648640*PolyGamma[0, a]^2*PolyGamma[1, a]*PolyGamma[4, a]*
      PolyGamma[6, a] + 4324320*PolyGamma[1, a]^2*PolyGamma[4, a]*
      PolyGamma[6, a] + 5765760*PolyGamma[0, a]*PolyGamma[2, a]*
      PolyGamma[4, a]*PolyGamma[6, a] + 1441440*PolyGamma[3, a]*
      PolyGamma[4, a]*PolyGamma[6, a] + 960960*PolyGamma[0, a]^3*
      PolyGamma[5, a]*PolyGamma[6, a] + 2882880*PolyGamma[0, a]*
      PolyGamma[1, a]*PolyGamma[5, a]*PolyGamma[6, a] +
     960960*PolyGamma[2, a]*PolyGamma[5, a]*PolyGamma[6, a] +
     205920*PolyGamma[0, a]^2*PolyGamma[6, a]^2 + 205920*PolyGamma[1, a]*
      PolyGamma[6, a]^2 + 12870*PolyGamma[0, a]^8*PolyGamma[7, a] +
     360360*PolyGamma[0, a]^6*PolyGamma[1, a]*PolyGamma[7, a] +
     2702700*PolyGamma[0, a]^4*PolyGamma[1, a]^2*PolyGamma[7, a] +
     5405400*PolyGamma[0, a]^2*PolyGamma[1, a]^3*PolyGamma[7, a] +
     1351350*PolyGamma[1, a]^4*PolyGamma[7, a] + 720720*PolyGamma[0, a]^5*
      PolyGamma[2, a]*PolyGamma[7, a] + 7207200*PolyGamma[0, a]^3*
      PolyGamma[1, a]*PolyGamma[2, a]*PolyGamma[7, a] +
     10810800*PolyGamma[0, a]*PolyGamma[1, a]^2*PolyGamma[2, a]*
      PolyGamma[7, a] + 3603600*PolyGamma[0, a]^2*PolyGamma[2, a]^2*
      PolyGamma[7, a] + 3603600*PolyGamma[1, a]*PolyGamma[2, a]^2*
      PolyGamma[7, a] + 900900*PolyGamma[0, a]^4*PolyGamma[3, a]*
      PolyGamma[7, a] + 5405400*PolyGamma[0, a]^2*PolyGamma[1, a]*
      PolyGamma[3, a]*PolyGamma[7, a] + 2702700*PolyGamma[1, a]^2*
      PolyGamma[3, a]*PolyGamma[7, a] + 3603600*PolyGamma[0, a]*
      PolyGamma[2, a]*PolyGamma[3, a]*PolyGamma[7, a] +
     450450*PolyGamma[3, a]^2*PolyGamma[7, a] + 720720*PolyGamma[0, a]^3*
      PolyGamma[4, a]*PolyGamma[7, a] + 2162160*PolyGamma[0, a]*
      PolyGamma[1, a]*PolyGamma[4, a]*PolyGamma[7, a] +
     720720*PolyGamma[2, a]*PolyGamma[4, a]*PolyGamma[7, a] +
     360360*PolyGamma[0, a]^2*PolyGamma[5, a]*PolyGamma[7, a] +
     360360*PolyGamma[1, a]*PolyGamma[5, a]*PolyGamma[7, a] +
     102960*PolyGamma[0, a]*PolyGamma[6, a]*PolyGamma[7, a] +
     6435*PolyGamma[7, a]^2 + 11440*PolyGamma[0, a]^7*PolyGamma[8, a] +
     240240*PolyGamma[0, a]^5*PolyGamma[1, a]*PolyGamma[8, a] +
     1201200*PolyGamma[0, a]^3*PolyGamma[1, a]^2*PolyGamma[8, a] +
     1201200*PolyGamma[0, a]*PolyGamma[1, a]^3*PolyGamma[8, a] +
     400400*PolyGamma[0, a]^4*PolyGamma[2, a]*PolyGamma[8, a] +
     2402400*PolyGamma[0, a]^2*PolyGamma[1, a]*PolyGamma[2, a]*
      PolyGamma[8, a] + 1201200*PolyGamma[1, a]^2*PolyGamma[2, a]*
      PolyGamma[8, a] + 800800*PolyGamma[0, a]*PolyGamma[2, a]^2*
      PolyGamma[8, a] + 400400*PolyGamma[0, a]^3*PolyGamma[3, a]*
      PolyGamma[8, a] + 1201200*PolyGamma[0, a]*PolyGamma[1, a]*
      PolyGamma[3, a]*PolyGamma[8, a] + 400400*PolyGamma[2, a]*
      PolyGamma[3, a]*PolyGamma[8, a] + 240240*PolyGamma[0, a]^2*
      PolyGamma[4, a]*PolyGamma[8, a] + 240240*PolyGamma[1, a]*
      PolyGamma[4, a]*PolyGamma[8, a] + 80080*PolyGamma[0, a]*PolyGamma[5, a]*
      PolyGamma[8, a] + 11440*PolyGamma[6, a]*PolyGamma[8, a] +
     8008*PolyGamma[0, a]^6*PolyGamma[9, a] + 120120*PolyGamma[0, a]^4*
      PolyGamma[1, a]*PolyGamma[9, a] + 360360*PolyGamma[0, a]^2*
      PolyGamma[1, a]^2*PolyGamma[9, a] + 120120*PolyGamma[1, a]^3*
      PolyGamma[9, a] + 160160*PolyGamma[0, a]^3*PolyGamma[2, a]*
      PolyGamma[9, a] + 480480*PolyGamma[0, a]*PolyGamma[1, a]*
      PolyGamma[2, a]*PolyGamma[9, a] + 80080*PolyGamma[2, a]^2*
      PolyGamma[9, a] + 120120*PolyGamma[0, a]^2*PolyGamma[3, a]*
      PolyGamma[9, a] + 120120*PolyGamma[1, a]*PolyGamma[3, a]*
      PolyGamma[9, a] + 48048*PolyGamma[0, a]*PolyGamma[4, a]*
      PolyGamma[9, a] + 8008*PolyGamma[5, a]*PolyGamma[9, a] +
     4368*PolyGamma[0, a]^5*PolyGamma[10, a] + 43680*PolyGamma[0, a]^3*
      PolyGamma[1, a]*PolyGamma[10, a] + 65520*PolyGamma[0, a]*
      PolyGamma[1, a]^2*PolyGamma[10, a] + 43680*PolyGamma[0, a]^2*
      PolyGamma[2, a]*PolyGamma[10, a] + 43680*PolyGamma[1, a]*
      PolyGamma[2, a]*PolyGamma[10, a] + 21840*PolyGamma[0, a]*
      PolyGamma[3, a]*PolyGamma[10, a] + 4368*PolyGamma[4, a]*
      PolyGamma[10, a] + 1820*PolyGamma[0, a]^4*PolyGamma[11, a] +
     10920*PolyGamma[0, a]^2*PolyGamma[1, a]*PolyGamma[11, a] +
     5460*PolyGamma[1, a]^2*PolyGamma[11, a] + 7280*PolyGamma[0, a]*
      PolyGamma[2, a]*PolyGamma[11, a] + 1820*PolyGamma[3, a]*
      PolyGamma[11, a] + 560*PolyGamma[0, a]^3*PolyGamma[12, a] +
     1680*PolyGamma[0, a]*PolyGamma[1, a]*PolyGamma[12, a] +
     560*PolyGamma[2, a]*PolyGamma[12, a] + 120*PolyGamma[0, a]^2*
      PolyGamma[13, a] + 120*PolyGamma[1, a]*PolyGamma[13, a] +
     16*PolyGamma[0, a]*PolyGamma[14, a] + PolyGamma[15, a]

MBexpGam[a_, 17] = PolyGamma[0, a]^17 + 136*PolyGamma[0, a]^15*
      PolyGamma[1, a] + 7140*PolyGamma[0, a]^13*PolyGamma[1, a]^2 +
     185640*PolyGamma[0, a]^11*PolyGamma[1, a]^3 + 2552550*PolyGamma[0, a]^9*
      PolyGamma[1, a]^4 + 18378360*PolyGamma[0, a]^7*PolyGamma[1, a]^5 +
     64324260*PolyGamma[0, a]^5*PolyGamma[1, a]^6 +
     91891800*PolyGamma[0, a]^3*PolyGamma[1, a]^7 +
     34459425*PolyGamma[0, a]*PolyGamma[1, a]^8 + 680*PolyGamma[0, a]^14*
      PolyGamma[2, a] + 61880*PolyGamma[0, a]^12*PolyGamma[1, a]*
      PolyGamma[2, a] + 2042040*PolyGamma[0, a]^10*PolyGamma[1, a]^2*
      PolyGamma[2, a] + 30630600*PolyGamma[0, a]^8*PolyGamma[1, a]^3*
      PolyGamma[2, a] + 214414200*PolyGamma[0, a]^6*PolyGamma[1, a]^4*
      PolyGamma[2, a] + 643242600*PolyGamma[0, a]^4*PolyGamma[1, a]^5*
      PolyGamma[2, a] + 643242600*PolyGamma[0, a]^2*PolyGamma[1, a]^6*
      PolyGamma[2, a] + 91891800*PolyGamma[1, a]^7*PolyGamma[2, a] +
     123760*PolyGamma[0, a]^11*PolyGamma[2, a]^2 + 6806800*PolyGamma[0, a]^9*
      PolyGamma[1, a]*PolyGamma[2, a]^2 + 122522400*PolyGamma[0, a]^7*
      PolyGamma[1, a]^2*PolyGamma[2, a]^2 + 857656800*PolyGamma[0, a]^5*
      PolyGamma[1, a]^3*PolyGamma[2, a]^2 + 2144142000*PolyGamma[0, a]^3*
      PolyGamma[1, a]^4*PolyGamma[2, a]^2 + 1286485200*PolyGamma[0, a]*
      PolyGamma[1, a]^5*PolyGamma[2, a]^2 + 6806800*PolyGamma[0, a]^8*
      PolyGamma[2, a]^3 + 190590400*PolyGamma[0, a]^6*PolyGamma[1, a]*
      PolyGamma[2, a]^3 + 1429428000*PolyGamma[0, a]^4*PolyGamma[1, a]^2*
      PolyGamma[2, a]^3 + 2858856000*PolyGamma[0, a]^2*PolyGamma[1, a]^3*
      PolyGamma[2, a]^3 + 714714000*PolyGamma[1, a]^4*PolyGamma[2, a]^3 +
     95295200*PolyGamma[0, a]^5*PolyGamma[2, a]^4 +
     952952000*PolyGamma[0, a]^3*PolyGamma[1, a]*PolyGamma[2, a]^4 +
     1429428000*PolyGamma[0, a]*PolyGamma[1, a]^2*PolyGamma[2, a]^4 +
     190590400*PolyGamma[0, a]^2*PolyGamma[2, a]^5 +
     190590400*PolyGamma[1, a]*PolyGamma[2, a]^5 + 2380*PolyGamma[0, a]^13*
      PolyGamma[3, a] + 185640*PolyGamma[0, a]^11*PolyGamma[1, a]*
      PolyGamma[3, a] + 5105100*PolyGamma[0, a]^9*PolyGamma[1, a]^2*
      PolyGamma[3, a] + 61261200*PolyGamma[0, a]^7*PolyGamma[1, a]^3*
      PolyGamma[3, a] + 321621300*PolyGamma[0, a]^5*PolyGamma[1, a]^4*
      PolyGamma[3, a] + 643242600*PolyGamma[0, a]^3*PolyGamma[1, a]^5*
      PolyGamma[3, a] + 321621300*PolyGamma[0, a]*PolyGamma[1, a]^6*
      PolyGamma[3, a] + 680680*PolyGamma[0, a]^10*PolyGamma[2, a]*
      PolyGamma[3, a] + 30630600*PolyGamma[0, a]^8*PolyGamma[1, a]*
      PolyGamma[2, a]*PolyGamma[3, a] + 428828400*PolyGamma[0, a]^6*
      PolyGamma[1, a]^2*PolyGamma[2, a]*PolyGamma[3, a] +
     2144142000*PolyGamma[0, a]^4*PolyGamma[1, a]^3*PolyGamma[2, a]*
      PolyGamma[3, a] + 3216213000*PolyGamma[0, a]^2*PolyGamma[1, a]^4*
      PolyGamma[2, a]*PolyGamma[3, a] + 643242600*PolyGamma[1, a]^5*
      PolyGamma[2, a]*PolyGamma[3, a] + 40840800*PolyGamma[0, a]^7*
      PolyGamma[2, a]^2*PolyGamma[3, a] + 857656800*PolyGamma[0, a]^5*
      PolyGamma[1, a]*PolyGamma[2, a]^2*PolyGamma[3, a] +
     4288284000*PolyGamma[0, a]^3*PolyGamma[1, a]^2*PolyGamma[2, a]^2*
      PolyGamma[3, a] + 4288284000*PolyGamma[0, a]*PolyGamma[1, a]^3*
      PolyGamma[2, a]^2*PolyGamma[3, a] + 476476000*PolyGamma[0, a]^4*
      PolyGamma[2, a]^3*PolyGamma[3, a] + 2858856000*PolyGamma[0, a]^2*
      PolyGamma[1, a]*PolyGamma[2, a]^3*PolyGamma[3, a] +
     1429428000*PolyGamma[1, a]^2*PolyGamma[2, a]^3*PolyGamma[3, a] +
     476476000*PolyGamma[0, a]*PolyGamma[2, a]^4*PolyGamma[3, a] +
     850850*PolyGamma[0, a]^9*PolyGamma[3, a]^2 + 30630600*PolyGamma[0, a]^7*
      PolyGamma[1, a]*PolyGamma[3, a]^2 + 321621300*PolyGamma[0, a]^5*
      PolyGamma[1, a]^2*PolyGamma[3, a]^2 + 1072071000*PolyGamma[0, a]^3*
      PolyGamma[1, a]^3*PolyGamma[3, a]^2 + 804053250*PolyGamma[0, a]*
      PolyGamma[1, a]^4*PolyGamma[3, a]^2 + 71471400*PolyGamma[0, a]^6*
      PolyGamma[2, a]*PolyGamma[3, a]^2 + 1072071000*PolyGamma[0, a]^4*
      PolyGamma[1, a]*PolyGamma[2, a]*PolyGamma[3, a]^2 +
     3216213000*PolyGamma[0, a]^2*PolyGamma[1, a]^2*PolyGamma[2, a]*
      PolyGamma[3, a]^2 + 1072071000*PolyGamma[1, a]^3*PolyGamma[2, a]*
      PolyGamma[3, a]^2 + 714714000*PolyGamma[0, a]^3*PolyGamma[2, a]^2*
      PolyGamma[3, a]^2 + 2144142000*PolyGamma[0, a]*PolyGamma[1, a]*
      PolyGamma[2, a]^2*PolyGamma[3, a]^2 + 238238000*PolyGamma[2, a]^3*
      PolyGamma[3, a]^2 + 35735700*PolyGamma[0, a]^5*PolyGamma[3, a]^3 +
     357357000*PolyGamma[0, a]^3*PolyGamma[1, a]*PolyGamma[3, a]^3 +
     536035500*PolyGamma[0, a]*PolyGamma[1, a]^2*PolyGamma[3, a]^3 +
     357357000*PolyGamma[0, a]^2*PolyGamma[2, a]*PolyGamma[3, a]^3 +
     357357000*PolyGamma[1, a]*PolyGamma[2, a]*PolyGamma[3, a]^3 +
     44669625*PolyGamma[0, a]*PolyGamma[3, a]^4 + 6188*PolyGamma[0, a]^12*
      PolyGamma[4, a] + 408408*PolyGamma[0, a]^10*PolyGamma[1, a]*
      PolyGamma[4, a] + 9189180*PolyGamma[0, a]^8*PolyGamma[1, a]^2*
      PolyGamma[4, a] + 85765680*PolyGamma[0, a]^6*PolyGamma[1, a]^3*
      PolyGamma[4, a] + 321621300*PolyGamma[0, a]^4*PolyGamma[1, a]^4*
      PolyGamma[4, a] + 385945560*PolyGamma[0, a]^2*PolyGamma[1, a]^5*
      PolyGamma[4, a] + 64324260*PolyGamma[1, a]^6*PolyGamma[4, a] +
     1361360*PolyGamma[0, a]^9*PolyGamma[2, a]*PolyGamma[4, a] +
     49008960*PolyGamma[0, a]^7*PolyGamma[1, a]*PolyGamma[2, a]*
      PolyGamma[4, a] + 514594080*PolyGamma[0, a]^5*PolyGamma[1, a]^2*
      PolyGamma[2, a]*PolyGamma[4, a] + 1715313600*PolyGamma[0, a]^3*
      PolyGamma[1, a]^3*PolyGamma[2, a]*PolyGamma[4, a] +
     1286485200*PolyGamma[0, a]*PolyGamma[1, a]^4*PolyGamma[2, a]*
      PolyGamma[4, a] + 57177120*PolyGamma[0, a]^6*PolyGamma[2, a]^2*
      PolyGamma[4, a] + 857656800*PolyGamma[0, a]^4*PolyGamma[1, a]*
      PolyGamma[2, a]^2*PolyGamma[4, a] + 2572970400*PolyGamma[0, a]^2*
      PolyGamma[1, a]^2*PolyGamma[2, a]^2*PolyGamma[4, a] +
     857656800*PolyGamma[1, a]^3*PolyGamma[2, a]^2*PolyGamma[4, a] +
     381180800*PolyGamma[0, a]^3*PolyGamma[2, a]^3*PolyGamma[4, a] +
     1143542400*PolyGamma[0, a]*PolyGamma[1, a]*PolyGamma[2, a]^3*
      PolyGamma[4, a] + 95295200*PolyGamma[2, a]^4*PolyGamma[4, a] +
     3063060*PolyGamma[0, a]^8*PolyGamma[3, a]*PolyGamma[4, a] +
     85765680*PolyGamma[0, a]^6*PolyGamma[1, a]*PolyGamma[3, a]*
      PolyGamma[4, a] + 643242600*PolyGamma[0, a]^4*PolyGamma[1, a]^2*
      PolyGamma[3, a]*PolyGamma[4, a] + 1286485200*PolyGamma[0, a]^2*
      PolyGamma[1, a]^3*PolyGamma[3, a]*PolyGamma[4, a] +
     321621300*PolyGamma[1, a]^4*PolyGamma[3, a]*PolyGamma[4, a] +
     171531360*PolyGamma[0, a]^5*PolyGamma[2, a]*PolyGamma[3, a]*
      PolyGamma[4, a] + 1715313600*PolyGamma[0, a]^3*PolyGamma[1, a]*
      PolyGamma[2, a]*PolyGamma[3, a]*PolyGamma[4, a] +
     2572970400*PolyGamma[0, a]*PolyGamma[1, a]^2*PolyGamma[2, a]*
      PolyGamma[3, a]*PolyGamma[4, a] + 857656800*PolyGamma[0, a]^2*
      PolyGamma[2, a]^2*PolyGamma[3, a]*PolyGamma[4, a] +
     857656800*PolyGamma[1, a]*PolyGamma[2, a]^2*PolyGamma[3, a]*
      PolyGamma[4, a] + 107207100*PolyGamma[0, a]^4*PolyGamma[3, a]^2*
      PolyGamma[4, a] + 643242600*PolyGamma[0, a]^2*PolyGamma[1, a]*
      PolyGamma[3, a]^2*PolyGamma[4, a] + 321621300*PolyGamma[1, a]^2*
      PolyGamma[3, a]^2*PolyGamma[4, a] + 428828400*PolyGamma[0, a]*
      PolyGamma[2, a]*PolyGamma[3, a]^2*PolyGamma[4, a] +
     35735700*PolyGamma[3, a]^3*PolyGamma[4, a] + 2450448*PolyGamma[0, a]^7*
      PolyGamma[4, a]^2 + 51459408*PolyGamma[0, a]^5*PolyGamma[1, a]*
      PolyGamma[4, a]^2 + 257297040*PolyGamma[0, a]^3*PolyGamma[1, a]^2*
      PolyGamma[4, a]^2 + 257297040*PolyGamma[0, a]*PolyGamma[1, a]^3*
      PolyGamma[4, a]^2 + 85765680*PolyGamma[0, a]^4*PolyGamma[2, a]*
      PolyGamma[4, a]^2 + 514594080*PolyGamma[0, a]^2*PolyGamma[1, a]*
      PolyGamma[2, a]*PolyGamma[4, a]^2 + 257297040*PolyGamma[1, a]^2*
      PolyGamma[2, a]*PolyGamma[4, a]^2 + 171531360*PolyGamma[0, a]*
      PolyGamma[2, a]^2*PolyGamma[4, a]^2 + 85765680*PolyGamma[0, a]^3*
      PolyGamma[3, a]*PolyGamma[4, a]^2 + 257297040*PolyGamma[0, a]*
      PolyGamma[1, a]*PolyGamma[3, a]*PolyGamma[4, a]^2 +
     85765680*PolyGamma[2, a]*PolyGamma[3, a]*PolyGamma[4, a]^2 +
     17153136*PolyGamma[0, a]^2*PolyGamma[4, a]^3 +
     17153136*PolyGamma[1, a]*PolyGamma[4, a]^3 + 12376*PolyGamma[0, a]^11*
      PolyGamma[5, a] + 680680*PolyGamma[0, a]^9*PolyGamma[1, a]*
      PolyGamma[5, a] + 12252240*PolyGamma[0, a]^7*PolyGamma[1, a]^2*
      PolyGamma[5, a] + 85765680*PolyGamma[0, a]^5*PolyGamma[1, a]^3*
      PolyGamma[5, a] + 214414200*PolyGamma[0, a]^3*PolyGamma[1, a]^4*
      PolyGamma[5, a] + 128648520*PolyGamma[0, a]*PolyGamma[1, a]^5*
      PolyGamma[5, a] + 2042040*PolyGamma[0, a]^8*PolyGamma[2, a]*
      PolyGamma[5, a] + 57177120*PolyGamma[0, a]^6*PolyGamma[1, a]*
      PolyGamma[2, a]*PolyGamma[5, a] + 428828400*PolyGamma[0, a]^4*
      PolyGamma[1, a]^2*PolyGamma[2, a]*PolyGamma[5, a] +
     857656800*PolyGamma[0, a]^2*PolyGamma[1, a]^3*PolyGamma[2, a]*
      PolyGamma[5, a] + 214414200*PolyGamma[1, a]^4*PolyGamma[2, a]*
      PolyGamma[5, a] + 57177120*PolyGamma[0, a]^5*PolyGamma[2, a]^2*
      PolyGamma[5, a] + 571771200*PolyGamma[0, a]^3*PolyGamma[1, a]*
      PolyGamma[2, a]^2*PolyGamma[5, a] + 857656800*PolyGamma[0, a]*
      PolyGamma[1, a]^2*PolyGamma[2, a]^2*PolyGamma[5, a] +
     190590400*PolyGamma[0, a]^2*PolyGamma[2, a]^3*PolyGamma[5, a] +
     190590400*PolyGamma[1, a]*PolyGamma[2, a]^3*PolyGamma[5, a] +
     4084080*PolyGamma[0, a]^7*PolyGamma[3, a]*PolyGamma[5, a] +
     85765680*PolyGamma[0, a]^5*PolyGamma[1, a]*PolyGamma[3, a]*
      PolyGamma[5, a] + 428828400*PolyGamma[0, a]^3*PolyGamma[1, a]^2*
      PolyGamma[3, a]*PolyGamma[5, a] + 428828400*PolyGamma[0, a]*
      PolyGamma[1, a]^3*PolyGamma[3, a]*PolyGamma[5, a] +
     142942800*PolyGamma[0, a]^4*PolyGamma[2, a]*PolyGamma[3, a]*
      PolyGamma[5, a] + 857656800*PolyGamma[0, a]^2*PolyGamma[1, a]*
      PolyGamma[2, a]*PolyGamma[3, a]*PolyGamma[5, a] +
     428828400*PolyGamma[1, a]^2*PolyGamma[2, a]*PolyGamma[3, a]*
      PolyGamma[5, a] + 285885600*PolyGamma[0, a]*PolyGamma[2, a]^2*
      PolyGamma[3, a]*PolyGamma[5, a] + 71471400*PolyGamma[0, a]^3*
      PolyGamma[3, a]^2*PolyGamma[5, a] + 214414200*PolyGamma[0, a]*
      PolyGamma[1, a]*PolyGamma[3, a]^2*PolyGamma[5, a] +
     71471400*PolyGamma[2, a]*PolyGamma[3, a]^2*PolyGamma[5, a] +
     5717712*PolyGamma[0, a]^6*PolyGamma[4, a]*PolyGamma[5, a] +
     85765680*PolyGamma[0, a]^4*PolyGamma[1, a]*PolyGamma[4, a]*
      PolyGamma[5, a] + 257297040*PolyGamma[0, a]^2*PolyGamma[1, a]^2*
      PolyGamma[4, a]*PolyGamma[5, a] + 85765680*PolyGamma[1, a]^3*
      PolyGamma[4, a]*PolyGamma[5, a] + 114354240*PolyGamma[0, a]^3*
      PolyGamma[2, a]*PolyGamma[4, a]*PolyGamma[5, a] +
     343062720*PolyGamma[0, a]*PolyGamma[1, a]*PolyGamma[2, a]*
      PolyGamma[4, a]*PolyGamma[5, a] + 57177120*PolyGamma[2, a]^2*
      PolyGamma[4, a]*PolyGamma[5, a] + 85765680*PolyGamma[0, a]^2*
      PolyGamma[3, a]*PolyGamma[4, a]*PolyGamma[5, a] +
     85765680*PolyGamma[1, a]*PolyGamma[3, a]*PolyGamma[4, a]*
      PolyGamma[5, a] + 17153136*PolyGamma[0, a]*PolyGamma[4, a]^2*
      PolyGamma[5, a] + 2858856*PolyGamma[0, a]^5*PolyGamma[5, a]^2 +
     28588560*PolyGamma[0, a]^3*PolyGamma[1, a]*PolyGamma[5, a]^2 +
     42882840*PolyGamma[0, a]*PolyGamma[1, a]^2*PolyGamma[5, a]^2 +
     28588560*PolyGamma[0, a]^2*PolyGamma[2, a]*PolyGamma[5, a]^2 +
     28588560*PolyGamma[1, a]*PolyGamma[2, a]*PolyGamma[5, a]^2 +
     14294280*PolyGamma[0, a]*PolyGamma[3, a]*PolyGamma[5, a]^2 +
     2858856*PolyGamma[4, a]*PolyGamma[5, a]^2 + 19448*PolyGamma[0, a]^10*
      PolyGamma[6, a] + 875160*PolyGamma[0, a]^8*PolyGamma[1, a]*
      PolyGamma[6, a] + 12252240*PolyGamma[0, a]^6*PolyGamma[1, a]^2*
      PolyGamma[6, a] + 61261200*PolyGamma[0, a]^4*PolyGamma[1, a]^3*
      PolyGamma[6, a] + 91891800*PolyGamma[0, a]^2*PolyGamma[1, a]^4*
      PolyGamma[6, a] + 18378360*PolyGamma[1, a]^5*PolyGamma[6, a] +
     2333760*PolyGamma[0, a]^7*PolyGamma[2, a]*PolyGamma[6, a] +
     49008960*PolyGamma[0, a]^5*PolyGamma[1, a]*PolyGamma[2, a]*
      PolyGamma[6, a] + 245044800*PolyGamma[0, a]^3*PolyGamma[1, a]^2*
      PolyGamma[2, a]*PolyGamma[6, a] + 245044800*PolyGamma[0, a]*
      PolyGamma[1, a]^3*PolyGamma[2, a]*PolyGamma[6, a] +
     40840800*PolyGamma[0, a]^4*PolyGamma[2, a]^2*PolyGamma[6, a] +
     245044800*PolyGamma[0, a]^2*PolyGamma[1, a]*PolyGamma[2, a]^2*
      PolyGamma[6, a] + 122522400*PolyGamma[1, a]^2*PolyGamma[2, a]^2*
      PolyGamma[6, a] + 54454400*PolyGamma[0, a]*PolyGamma[2, a]^3*
      PolyGamma[6, a] + 4084080*PolyGamma[0, a]^6*PolyGamma[3, a]*
      PolyGamma[6, a] + 61261200*PolyGamma[0, a]^4*PolyGamma[1, a]*
      PolyGamma[3, a]*PolyGamma[6, a] + 183783600*PolyGamma[0, a]^2*
      PolyGamma[1, a]^2*PolyGamma[3, a]*PolyGamma[6, a] +
     61261200*PolyGamma[1, a]^3*PolyGamma[3, a]*PolyGamma[6, a] +
     81681600*PolyGamma[0, a]^3*PolyGamma[2, a]*PolyGamma[3, a]*
      PolyGamma[6, a] + 245044800*PolyGamma[0, a]*PolyGamma[1, a]*
      PolyGamma[2, a]*PolyGamma[3, a]*PolyGamma[6, a] +
     40840800*PolyGamma[2, a]^2*PolyGamma[3, a]*PolyGamma[6, a] +
     30630600*PolyGamma[0, a]^2*PolyGamma[3, a]^2*PolyGamma[6, a] +
     30630600*PolyGamma[1, a]*PolyGamma[3, a]^2*PolyGamma[6, a] +
     4900896*PolyGamma[0, a]^5*PolyGamma[4, a]*PolyGamma[6, a] +
     49008960*PolyGamma[0, a]^3*PolyGamma[1, a]*PolyGamma[4, a]*
      PolyGamma[6, a] + 73513440*PolyGamma[0, a]*PolyGamma[1, a]^2*
      PolyGamma[4, a]*PolyGamma[6, a] + 49008960*PolyGamma[0, a]^2*
      PolyGamma[2, a]*PolyGamma[4, a]*PolyGamma[6, a] +
     49008960*PolyGamma[1, a]*PolyGamma[2, a]*PolyGamma[4, a]*
      PolyGamma[6, a] + 24504480*PolyGamma[0, a]*PolyGamma[3, a]*
      PolyGamma[4, a]*PolyGamma[6, a] + 2450448*PolyGamma[4, a]^2*
      PolyGamma[6, a] + 4084080*PolyGamma[0, a]^4*PolyGamma[5, a]*
      PolyGamma[6, a] + 24504480*PolyGamma[0, a]^2*PolyGamma[1, a]*
      PolyGamma[5, a]*PolyGamma[6, a] + 12252240*PolyGamma[1, a]^2*
      PolyGamma[5, a]*PolyGamma[6, a] + 16336320*PolyGamma[0, a]*
      PolyGamma[2, a]*PolyGamma[5, a]*PolyGamma[6, a] +
     4084080*PolyGamma[3, a]*PolyGamma[5, a]*PolyGamma[6, a] +
     1166880*PolyGamma[0, a]^3*PolyGamma[6, a]^2 + 3500640*PolyGamma[0, a]*
      PolyGamma[1, a]*PolyGamma[6, a]^2 + 1166880*PolyGamma[2, a]*
      PolyGamma[6, a]^2 + 24310*PolyGamma[0, a]^9*PolyGamma[7, a] +
     875160*PolyGamma[0, a]^7*PolyGamma[1, a]*PolyGamma[7, a] +
     9189180*PolyGamma[0, a]^5*PolyGamma[1, a]^2*PolyGamma[7, a] +
     30630600*PolyGamma[0, a]^3*PolyGamma[1, a]^3*PolyGamma[7, a] +
     22972950*PolyGamma[0, a]*PolyGamma[1, a]^4*PolyGamma[7, a] +
     2042040*PolyGamma[0, a]^6*PolyGamma[2, a]*PolyGamma[7, a] +
     30630600*PolyGamma[0, a]^4*PolyGamma[1, a]*PolyGamma[2, a]*
      PolyGamma[7, a] + 91891800*PolyGamma[0, a]^2*PolyGamma[1, a]^2*
      PolyGamma[2, a]*PolyGamma[7, a] + 30630600*PolyGamma[1, a]^3*
      PolyGamma[2, a]*PolyGamma[7, a] + 20420400*PolyGamma[0, a]^3*
      PolyGamma[2, a]^2*PolyGamma[7, a] + 61261200*PolyGamma[0, a]*
      PolyGamma[1, a]*PolyGamma[2, a]^2*PolyGamma[7, a] +
     6806800*PolyGamma[2, a]^3*PolyGamma[7, a] + 3063060*PolyGamma[0, a]^5*
      PolyGamma[3, a]*PolyGamma[7, a] + 30630600*PolyGamma[0, a]^3*
      PolyGamma[1, a]*PolyGamma[3, a]*PolyGamma[7, a] +
     45945900*PolyGamma[0, a]*PolyGamma[1, a]^2*PolyGamma[3, a]*
      PolyGamma[7, a] + 30630600*PolyGamma[0, a]^2*PolyGamma[2, a]*
      PolyGamma[3, a]*PolyGamma[7, a] + 30630600*PolyGamma[1, a]*
      PolyGamma[2, a]*PolyGamma[3, a]*PolyGamma[7, a] +
     7657650*PolyGamma[0, a]*PolyGamma[3, a]^2*PolyGamma[7, a] +
     3063060*PolyGamma[0, a]^4*PolyGamma[4, a]*PolyGamma[7, a] +
     18378360*PolyGamma[0, a]^2*PolyGamma[1, a]*PolyGamma[4, a]*
      PolyGamma[7, a] + 9189180*PolyGamma[1, a]^2*PolyGamma[4, a]*
      PolyGamma[7, a] + 12252240*PolyGamma[0, a]*PolyGamma[2, a]*
      PolyGamma[4, a]*PolyGamma[7, a] + 3063060*PolyGamma[3, a]*
      PolyGamma[4, a]*PolyGamma[7, a] + 2042040*PolyGamma[0, a]^3*
      PolyGamma[5, a]*PolyGamma[7, a] + 6126120*PolyGamma[0, a]*
      PolyGamma[1, a]*PolyGamma[5, a]*PolyGamma[7, a] +
     2042040*PolyGamma[2, a]*PolyGamma[5, a]*PolyGamma[7, a] +
     875160*PolyGamma[0, a]^2*PolyGamma[6, a]*PolyGamma[7, a] +
     875160*PolyGamma[1, a]*PolyGamma[6, a]*PolyGamma[7, a] +
     109395*PolyGamma[0, a]*PolyGamma[7, a]^2 + 24310*PolyGamma[0, a]^8*
      PolyGamma[8, a] + 680680*PolyGamma[0, a]^6*PolyGamma[1, a]*
      PolyGamma[8, a] + 5105100*PolyGamma[0, a]^4*PolyGamma[1, a]^2*
      PolyGamma[8, a] + 10210200*PolyGamma[0, a]^2*PolyGamma[1, a]^3*
      PolyGamma[8, a] + 2552550*PolyGamma[1, a]^4*PolyGamma[8, a] +
     1361360*PolyGamma[0, a]^5*PolyGamma[2, a]*PolyGamma[8, a] +
     13613600*PolyGamma[0, a]^3*PolyGamma[1, a]*PolyGamma[2, a]*
      PolyGamma[8, a] + 20420400*PolyGamma[0, a]*PolyGamma[1, a]^2*
      PolyGamma[2, a]*PolyGamma[8, a] + 6806800*PolyGamma[0, a]^2*
      PolyGamma[2, a]^2*PolyGamma[8, a] + 6806800*PolyGamma[1, a]*
      PolyGamma[2, a]^2*PolyGamma[8, a] + 1701700*PolyGamma[0, a]^4*
      PolyGamma[3, a]*PolyGamma[8, a] + 10210200*PolyGamma[0, a]^2*
      PolyGamma[1, a]*PolyGamma[3, a]*PolyGamma[8, a] +
     5105100*PolyGamma[1, a]^2*PolyGamma[3, a]*PolyGamma[8, a] +
     6806800*PolyGamma[0, a]*PolyGamma[2, a]*PolyGamma[3, a]*
      PolyGamma[8, a] + 850850*PolyGamma[3, a]^2*PolyGamma[8, a] +
     1361360*PolyGamma[0, a]^3*PolyGamma[4, a]*PolyGamma[8, a] +
     4084080*PolyGamma[0, a]*PolyGamma[1, a]*PolyGamma[4, a]*
      PolyGamma[8, a] + 1361360*PolyGamma[2, a]*PolyGamma[4, a]*
      PolyGamma[8, a] + 680680*PolyGamma[0, a]^2*PolyGamma[5, a]*
      PolyGamma[8, a] + 680680*PolyGamma[1, a]*PolyGamma[5, a]*
      PolyGamma[8, a] + 194480*PolyGamma[0, a]*PolyGamma[6, a]*
      PolyGamma[8, a] + 24310*PolyGamma[7, a]*PolyGamma[8, a] +
     19448*PolyGamma[0, a]^7*PolyGamma[9, a] + 408408*PolyGamma[0, a]^5*
      PolyGamma[1, a]*PolyGamma[9, a] + 2042040*PolyGamma[0, a]^3*
      PolyGamma[1, a]^2*PolyGamma[9, a] + 2042040*PolyGamma[0, a]*
      PolyGamma[1, a]^3*PolyGamma[9, a] + 680680*PolyGamma[0, a]^4*
      PolyGamma[2, a]*PolyGamma[9, a] + 4084080*PolyGamma[0, a]^2*
      PolyGamma[1, a]*PolyGamma[2, a]*PolyGamma[9, a] +
     2042040*PolyGamma[1, a]^2*PolyGamma[2, a]*PolyGamma[9, a] +
     1361360*PolyGamma[0, a]*PolyGamma[2, a]^2*PolyGamma[9, a] +
     680680*PolyGamma[0, a]^3*PolyGamma[3, a]*PolyGamma[9, a] +
     2042040*PolyGamma[0, a]*PolyGamma[1, a]*PolyGamma[3, a]*
      PolyGamma[9, a] + 680680*PolyGamma[2, a]*PolyGamma[3, a]*
      PolyGamma[9, a] + 408408*PolyGamma[0, a]^2*PolyGamma[4, a]*
      PolyGamma[9, a] + 408408*PolyGamma[1, a]*PolyGamma[4, a]*
      PolyGamma[9, a] + 136136*PolyGamma[0, a]*PolyGamma[5, a]*
      PolyGamma[9, a] + 19448*PolyGamma[6, a]*PolyGamma[9, a] +
     12376*PolyGamma[0, a]^6*PolyGamma[10, a] + 185640*PolyGamma[0, a]^4*
      PolyGamma[1, a]*PolyGamma[10, a] + 556920*PolyGamma[0, a]^2*
      PolyGamma[1, a]^2*PolyGamma[10, a] + 185640*PolyGamma[1, a]^3*
      PolyGamma[10, a] + 247520*PolyGamma[0, a]^3*PolyGamma[2, a]*
      PolyGamma[10, a] + 742560*PolyGamma[0, a]*PolyGamma[1, a]*
      PolyGamma[2, a]*PolyGamma[10, a] + 123760*PolyGamma[2, a]^2*
      PolyGamma[10, a] + 185640*PolyGamma[0, a]^2*PolyGamma[3, a]*
      PolyGamma[10, a] + 185640*PolyGamma[1, a]*PolyGamma[3, a]*
      PolyGamma[10, a] + 74256*PolyGamma[0, a]*PolyGamma[4, a]*
      PolyGamma[10, a] + 12376*PolyGamma[5, a]*PolyGamma[10, a] +
     6188*PolyGamma[0, a]^5*PolyGamma[11, a] + 61880*PolyGamma[0, a]^3*
      PolyGamma[1, a]*PolyGamma[11, a] + 92820*PolyGamma[0, a]*
      PolyGamma[1, a]^2*PolyGamma[11, a] + 61880*PolyGamma[0, a]^2*
      PolyGamma[2, a]*PolyGamma[11, a] + 61880*PolyGamma[1, a]*
      PolyGamma[2, a]*PolyGamma[11, a] + 30940*PolyGamma[0, a]*
      PolyGamma[3, a]*PolyGamma[11, a] + 6188*PolyGamma[4, a]*
      PolyGamma[11, a] + 2380*PolyGamma[0, a]^4*PolyGamma[12, a] +
     14280*PolyGamma[0, a]^2*PolyGamma[1, a]*PolyGamma[12, a] +
     7140*PolyGamma[1, a]^2*PolyGamma[12, a] + 9520*PolyGamma[0, a]*
      PolyGamma[2, a]*PolyGamma[12, a] + 2380*PolyGamma[3, a]*
      PolyGamma[12, a] + 680*PolyGamma[0, a]^3*PolyGamma[13, a] +
     2040*PolyGamma[0, a]*PolyGamma[1, a]*PolyGamma[13, a] +
     680*PolyGamma[2, a]*PolyGamma[13, a] + 136*PolyGamma[0, a]^2*
      PolyGamma[14, a] + 136*PolyGamma[1, a]*PolyGamma[14, a] +
     17*PolyGamma[0, a]*PolyGamma[15, a] + PolyGamma[16, a]

MBexpGam[a_, 18] = PolyGamma[0, a]^18 + 153*PolyGamma[0, a]^16*
      PolyGamma[1, a] + 9180*PolyGamma[0, a]^14*PolyGamma[1, a]^2 +
     278460*PolyGamma[0, a]^12*PolyGamma[1, a]^3 + 4594590*PolyGamma[0, a]^10*
      PolyGamma[1, a]^4 + 41351310*PolyGamma[0, a]^8*PolyGamma[1, a]^5 +
     192972780*PolyGamma[0, a]^6*PolyGamma[1, a]^6 +
     413513100*PolyGamma[0, a]^4*PolyGamma[1, a]^7 +
     310134825*PolyGamma[0, a]^2*PolyGamma[1, a]^8 +
     34459425*PolyGamma[1, a]^9 + 816*PolyGamma[0, a]^15*PolyGamma[2, a] +
     85680*PolyGamma[0, a]^13*PolyGamma[1, a]*PolyGamma[2, a] +
     3341520*PolyGamma[0, a]^11*PolyGamma[1, a]^2*PolyGamma[2, a] +
     61261200*PolyGamma[0, a]^9*PolyGamma[1, a]^3*PolyGamma[2, a] +
     551350800*PolyGamma[0, a]^7*PolyGamma[1, a]^4*PolyGamma[2, a] +
     2315673360*PolyGamma[0, a]^5*PolyGamma[1, a]^5*PolyGamma[2, a] +
     3859455600*PolyGamma[0, a]^3*PolyGamma[1, a]^6*PolyGamma[2, a] +
     1654052400*PolyGamma[0, a]*PolyGamma[1, a]^7*PolyGamma[2, a] +
     185640*PolyGamma[0, a]^12*PolyGamma[2, a]^2 +
     12252240*PolyGamma[0, a]^10*PolyGamma[1, a]*PolyGamma[2, a]^2 +
     275675400*PolyGamma[0, a]^8*PolyGamma[1, a]^2*PolyGamma[2, a]^2 +
     2572970400*PolyGamma[0, a]^6*PolyGamma[1, a]^3*PolyGamma[2, a]^2 +
     9648639000*PolyGamma[0, a]^4*PolyGamma[1, a]^4*PolyGamma[2, a]^2 +
     11578366800*PolyGamma[0, a]^2*PolyGamma[1, a]^5*PolyGamma[2, a]^2 +
     1929727800*PolyGamma[1, a]^6*PolyGamma[2, a]^2 +
     13613600*PolyGamma[0, a]^9*PolyGamma[2, a]^3 +
     490089600*PolyGamma[0, a]^7*PolyGamma[1, a]*PolyGamma[2, a]^3 +
     5145940800*PolyGamma[0, a]^5*PolyGamma[1, a]^2*PolyGamma[2, a]^3 +
     17153136000*PolyGamma[0, a]^3*PolyGamma[1, a]^3*PolyGamma[2, a]^3 +
     12864852000*PolyGamma[0, a]*PolyGamma[1, a]^4*PolyGamma[2, a]^3 +
     285885600*PolyGamma[0, a]^6*PolyGamma[2, a]^4 +
     4288284000*PolyGamma[0, a]^4*PolyGamma[1, a]*PolyGamma[2, a]^4 +
     12864852000*PolyGamma[0, a]^2*PolyGamma[1, a]^2*PolyGamma[2, a]^4 +
     4288284000*PolyGamma[1, a]^3*PolyGamma[2, a]^4 +
     1143542400*PolyGamma[0, a]^3*PolyGamma[2, a]^5 +
     3430627200*PolyGamma[0, a]*PolyGamma[1, a]*PolyGamma[2, a]^5 +
     190590400*PolyGamma[2, a]^6 + 3060*PolyGamma[0, a]^14*PolyGamma[3, a] +
     278460*PolyGamma[0, a]^12*PolyGamma[1, a]*PolyGamma[3, a] +
     9189180*PolyGamma[0, a]^10*PolyGamma[1, a]^2*PolyGamma[3, a] +
     137837700*PolyGamma[0, a]^8*PolyGamma[1, a]^3*PolyGamma[3, a] +
     964863900*PolyGamma[0, a]^6*PolyGamma[1, a]^4*PolyGamma[3, a] +
     2894591700*PolyGamma[0, a]^4*PolyGamma[1, a]^5*PolyGamma[3, a] +
     2894591700*PolyGamma[0, a]^2*PolyGamma[1, a]^6*PolyGamma[3, a] +
     413513100*PolyGamma[1, a]^7*PolyGamma[3, a] + 1113840*PolyGamma[0, a]^11*
      PolyGamma[2, a]*PolyGamma[3, a] + 61261200*PolyGamma[0, a]^9*
      PolyGamma[1, a]*PolyGamma[2, a]*PolyGamma[3, a] +
     1102701600*PolyGamma[0, a]^7*PolyGamma[1, a]^2*PolyGamma[2, a]*
      PolyGamma[3, a] + 7718911200*PolyGamma[0, a]^5*PolyGamma[1, a]^3*
      PolyGamma[2, a]*PolyGamma[3, a] + 19297278000*PolyGamma[0, a]^3*
      PolyGamma[1, a]^4*PolyGamma[2, a]*PolyGamma[3, a] +
     11578366800*PolyGamma[0, a]*PolyGamma[1, a]^5*PolyGamma[2, a]*
      PolyGamma[3, a] + 91891800*PolyGamma[0, a]^8*PolyGamma[2, a]^2*
      PolyGamma[3, a] + 2572970400*PolyGamma[0, a]^6*PolyGamma[1, a]*
      PolyGamma[2, a]^2*PolyGamma[3, a] + 19297278000*PolyGamma[0, a]^4*
      PolyGamma[1, a]^2*PolyGamma[2, a]^2*PolyGamma[3, a] +
     38594556000*PolyGamma[0, a]^2*PolyGamma[1, a]^3*PolyGamma[2, a]^2*
      PolyGamma[3, a] + 9648639000*PolyGamma[1, a]^4*PolyGamma[2, a]^2*
      PolyGamma[3, a] + 1715313600*PolyGamma[0, a]^5*PolyGamma[2, a]^3*
      PolyGamma[3, a] + 17153136000*PolyGamma[0, a]^3*PolyGamma[1, a]*
      PolyGamma[2, a]^3*PolyGamma[3, a] + 25729704000*PolyGamma[0, a]*
      PolyGamma[1, a]^2*PolyGamma[2, a]^3*PolyGamma[3, a] +
     4288284000*PolyGamma[0, a]^2*PolyGamma[2, a]^4*PolyGamma[3, a] +
     4288284000*PolyGamma[1, a]*PolyGamma[2, a]^4*PolyGamma[3, a] +
     1531530*PolyGamma[0, a]^10*PolyGamma[3, a]^2 +
     68918850*PolyGamma[0, a]^8*PolyGamma[1, a]*PolyGamma[3, a]^2 +
     964863900*PolyGamma[0, a]^6*PolyGamma[1, a]^2*PolyGamma[3, a]^2 +
     4824319500*PolyGamma[0, a]^4*PolyGamma[1, a]^3*PolyGamma[3, a]^2 +
     7236479250*PolyGamma[0, a]^2*PolyGamma[1, a]^4*PolyGamma[3, a]^2 +
     1447295850*PolyGamma[1, a]^5*PolyGamma[3, a]^2 +
     183783600*PolyGamma[0, a]^7*PolyGamma[2, a]*PolyGamma[3, a]^2 +
     3859455600*PolyGamma[0, a]^5*PolyGamma[1, a]*PolyGamma[2, a]*
      PolyGamma[3, a]^2 + 19297278000*PolyGamma[0, a]^3*PolyGamma[1, a]^2*
      PolyGamma[2, a]*PolyGamma[3, a]^2 + 19297278000*PolyGamma[0, a]*
      PolyGamma[1, a]^3*PolyGamma[2, a]*PolyGamma[3, a]^2 +
     3216213000*PolyGamma[0, a]^4*PolyGamma[2, a]^2*PolyGamma[3, a]^2 +
     19297278000*PolyGamma[0, a]^2*PolyGamma[1, a]*PolyGamma[2, a]^2*
      PolyGamma[3, a]^2 + 9648639000*PolyGamma[1, a]^2*PolyGamma[2, a]^2*
      PolyGamma[3, a]^2 + 4288284000*PolyGamma[0, a]*PolyGamma[2, a]^3*
      PolyGamma[3, a]^2 + 107207100*PolyGamma[0, a]^6*PolyGamma[3, a]^3 +
     1608106500*PolyGamma[0, a]^4*PolyGamma[1, a]*PolyGamma[3, a]^3 +
     4824319500*PolyGamma[0, a]^2*PolyGamma[1, a]^2*PolyGamma[3, a]^3 +
     1608106500*PolyGamma[1, a]^3*PolyGamma[3, a]^3 +
     2144142000*PolyGamma[0, a]^3*PolyGamma[2, a]*PolyGamma[3, a]^3 +
     6432426000*PolyGamma[0, a]*PolyGamma[1, a]*PolyGamma[2, a]*
      PolyGamma[3, a]^3 + 1072071000*PolyGamma[2, a]^2*PolyGamma[3, a]^3 +
     402026625*PolyGamma[0, a]^2*PolyGamma[3, a]^4 +
     402026625*PolyGamma[1, a]*PolyGamma[3, a]^4 + 8568*PolyGamma[0, a]^13*
      PolyGamma[4, a] + 668304*PolyGamma[0, a]^11*PolyGamma[1, a]*
      PolyGamma[4, a] + 18378360*PolyGamma[0, a]^9*PolyGamma[1, a]^2*
      PolyGamma[4, a] + 220540320*PolyGamma[0, a]^7*PolyGamma[1, a]^3*
      PolyGamma[4, a] + 1157836680*PolyGamma[0, a]^5*PolyGamma[1, a]^4*
      PolyGamma[4, a] + 2315673360*PolyGamma[0, a]^3*PolyGamma[1, a]^5*
      PolyGamma[4, a] + 1157836680*PolyGamma[0, a]*PolyGamma[1, a]^6*
      PolyGamma[4, a] + 2450448*PolyGamma[0, a]^10*PolyGamma[2, a]*
      PolyGamma[4, a] + 110270160*PolyGamma[0, a]^8*PolyGamma[1, a]*
      PolyGamma[2, a]*PolyGamma[4, a] + 1543782240*PolyGamma[0, a]^6*
      PolyGamma[1, a]^2*PolyGamma[2, a]*PolyGamma[4, a] +
     7718911200*PolyGamma[0, a]^4*PolyGamma[1, a]^3*PolyGamma[2, a]*
      PolyGamma[4, a] + 11578366800*PolyGamma[0, a]^2*PolyGamma[1, a]^4*
      PolyGamma[2, a]*PolyGamma[4, a] + 2315673360*PolyGamma[1, a]^5*
      PolyGamma[2, a]*PolyGamma[4, a] + 147026880*PolyGamma[0, a]^7*
      PolyGamma[2, a]^2*PolyGamma[4, a] + 3087564480*PolyGamma[0, a]^5*
      PolyGamma[1, a]*PolyGamma[2, a]^2*PolyGamma[4, a] +
     15437822400*PolyGamma[0, a]^3*PolyGamma[1, a]^2*PolyGamma[2, a]^2*
      PolyGamma[4, a] + 15437822400*PolyGamma[0, a]*PolyGamma[1, a]^3*
      PolyGamma[2, a]^2*PolyGamma[4, a] + 1715313600*PolyGamma[0, a]^4*
      PolyGamma[2, a]^3*PolyGamma[4, a] + 10291881600*PolyGamma[0, a]^2*
      PolyGamma[1, a]*PolyGamma[2, a]^3*PolyGamma[4, a] +
     5145940800*PolyGamma[1, a]^2*PolyGamma[2, a]^3*PolyGamma[4, a] +
     1715313600*PolyGamma[0, a]*PolyGamma[2, a]^4*PolyGamma[4, a] +
     6126120*PolyGamma[0, a]^9*PolyGamma[3, a]*PolyGamma[4, a] +
     220540320*PolyGamma[0, a]^7*PolyGamma[1, a]*PolyGamma[3, a]*
      PolyGamma[4, a] + 2315673360*PolyGamma[0, a]^5*PolyGamma[1, a]^2*
      PolyGamma[3, a]*PolyGamma[4, a] + 7718911200*PolyGamma[0, a]^3*
      PolyGamma[1, a]^3*PolyGamma[3, a]*PolyGamma[4, a] +
     5789183400*PolyGamma[0, a]*PolyGamma[1, a]^4*PolyGamma[3, a]*
      PolyGamma[4, a] + 514594080*PolyGamma[0, a]^6*PolyGamma[2, a]*
      PolyGamma[3, a]*PolyGamma[4, a] + 7718911200*PolyGamma[0, a]^4*
      PolyGamma[1, a]*PolyGamma[2, a]*PolyGamma[3, a]*PolyGamma[4, a] +
     23156733600*PolyGamma[0, a]^2*PolyGamma[1, a]^2*PolyGamma[2, a]*
      PolyGamma[3, a]*PolyGamma[4, a] + 7718911200*PolyGamma[1, a]^3*
      PolyGamma[2, a]*PolyGamma[3, a]*PolyGamma[4, a] +
     5145940800*PolyGamma[0, a]^3*PolyGamma[2, a]^2*PolyGamma[3, a]*
      PolyGamma[4, a] + 15437822400*PolyGamma[0, a]*PolyGamma[1, a]*
      PolyGamma[2, a]^2*PolyGamma[3, a]*PolyGamma[4, a] +
     1715313600*PolyGamma[2, a]^3*PolyGamma[3, a]*PolyGamma[4, a] +
     385945560*PolyGamma[0, a]^5*PolyGamma[3, a]^2*PolyGamma[4, a] +
     3859455600*PolyGamma[0, a]^3*PolyGamma[1, a]*PolyGamma[3, a]^2*
      PolyGamma[4, a] + 5789183400*PolyGamma[0, a]*PolyGamma[1, a]^2*
      PolyGamma[3, a]^2*PolyGamma[4, a] + 3859455600*PolyGamma[0, a]^2*
      PolyGamma[2, a]*PolyGamma[3, a]^2*PolyGamma[4, a] +
     3859455600*PolyGamma[1, a]*PolyGamma[2, a]*PolyGamma[3, a]^2*
      PolyGamma[4, a] + 643242600*PolyGamma[0, a]*PolyGamma[3, a]^3*
      PolyGamma[4, a] + 5513508*PolyGamma[0, a]^8*PolyGamma[4, a]^2 +
     154378224*PolyGamma[0, a]^6*PolyGamma[1, a]*PolyGamma[4, a]^2 +
     1157836680*PolyGamma[0, a]^4*PolyGamma[1, a]^2*PolyGamma[4, a]^2 +
     2315673360*PolyGamma[0, a]^2*PolyGamma[1, a]^3*PolyGamma[4, a]^2 +
     578918340*PolyGamma[1, a]^4*PolyGamma[4, a]^2 +
     308756448*PolyGamma[0, a]^5*PolyGamma[2, a]*PolyGamma[4, a]^2 +
     3087564480*PolyGamma[0, a]^3*PolyGamma[1, a]*PolyGamma[2, a]*
      PolyGamma[4, a]^2 + 4631346720*PolyGamma[0, a]*PolyGamma[1, a]^2*
      PolyGamma[2, a]*PolyGamma[4, a]^2 + 1543782240*PolyGamma[0, a]^2*
      PolyGamma[2, a]^2*PolyGamma[4, a]^2 + 1543782240*PolyGamma[1, a]*
      PolyGamma[2, a]^2*PolyGamma[4, a]^2 + 385945560*PolyGamma[0, a]^4*
      PolyGamma[3, a]*PolyGamma[4, a]^2 + 2315673360*PolyGamma[0, a]^2*
      PolyGamma[1, a]*PolyGamma[3, a]*PolyGamma[4, a]^2 +
     1157836680*PolyGamma[1, a]^2*PolyGamma[3, a]*PolyGamma[4, a]^2 +
     1543782240*PolyGamma[0, a]*PolyGamma[2, a]*PolyGamma[3, a]*
      PolyGamma[4, a]^2 + 192972780*PolyGamma[3, a]^2*PolyGamma[4, a]^2 +
     102918816*PolyGamma[0, a]^3*PolyGamma[4, a]^3 +
     308756448*PolyGamma[0, a]*PolyGamma[1, a]*PolyGamma[4, a]^3 +
     102918816*PolyGamma[2, a]*PolyGamma[4, a]^3 + 18564*PolyGamma[0, a]^12*
      PolyGamma[5, a] + 1225224*PolyGamma[0, a]^10*PolyGamma[1, a]*
      PolyGamma[5, a] + 27567540*PolyGamma[0, a]^8*PolyGamma[1, a]^2*
      PolyGamma[5, a] + 257297040*PolyGamma[0, a]^6*PolyGamma[1, a]^3*
      PolyGamma[5, a] + 964863900*PolyGamma[0, a]^4*PolyGamma[1, a]^4*
      PolyGamma[5, a] + 1157836680*PolyGamma[0, a]^2*PolyGamma[1, a]^5*
      PolyGamma[5, a] + 192972780*PolyGamma[1, a]^6*PolyGamma[5, a] +
     4084080*PolyGamma[0, a]^9*PolyGamma[2, a]*PolyGamma[5, a] +
     147026880*PolyGamma[0, a]^7*PolyGamma[1, a]*PolyGamma[2, a]*
      PolyGamma[5, a] + 1543782240*PolyGamma[0, a]^5*PolyGamma[1, a]^2*
      PolyGamma[2, a]*PolyGamma[5, a] + 5145940800*PolyGamma[0, a]^3*
      PolyGamma[1, a]^3*PolyGamma[2, a]*PolyGamma[5, a] +
     3859455600*PolyGamma[0, a]*PolyGamma[1, a]^4*PolyGamma[2, a]*
      PolyGamma[5, a] + 171531360*PolyGamma[0, a]^6*PolyGamma[2, a]^2*
      PolyGamma[5, a] + 2572970400*PolyGamma[0, a]^4*PolyGamma[1, a]*
      PolyGamma[2, a]^2*PolyGamma[5, a] + 7718911200*PolyGamma[0, a]^2*
      PolyGamma[1, a]^2*PolyGamma[2, a]^2*PolyGamma[5, a] +
     2572970400*PolyGamma[1, a]^3*PolyGamma[2, a]^2*PolyGamma[5, a] +
     1143542400*PolyGamma[0, a]^3*PolyGamma[2, a]^3*PolyGamma[5, a] +
     3430627200*PolyGamma[0, a]*PolyGamma[1, a]*PolyGamma[2, a]^3*
      PolyGamma[5, a] + 285885600*PolyGamma[2, a]^4*PolyGamma[5, a] +
     9189180*PolyGamma[0, a]^8*PolyGamma[3, a]*PolyGamma[5, a] +
     257297040*PolyGamma[0, a]^6*PolyGamma[1, a]*PolyGamma[3, a]*
      PolyGamma[5, a] + 1929727800*PolyGamma[0, a]^4*PolyGamma[1, a]^2*
      PolyGamma[3, a]*PolyGamma[5, a] + 3859455600*PolyGamma[0, a]^2*
      PolyGamma[1, a]^3*PolyGamma[3, a]*PolyGamma[5, a] +
     964863900*PolyGamma[1, a]^4*PolyGamma[3, a]*PolyGamma[5, a] +
     514594080*PolyGamma[0, a]^5*PolyGamma[2, a]*PolyGamma[3, a]*
      PolyGamma[5, a] + 5145940800*PolyGamma[0, a]^3*PolyGamma[1, a]*
      PolyGamma[2, a]*PolyGamma[3, a]*PolyGamma[5, a] +
     7718911200*PolyGamma[0, a]*PolyGamma[1, a]^2*PolyGamma[2, a]*
      PolyGamma[3, a]*PolyGamma[5, a] + 2572970400*PolyGamma[0, a]^2*
      PolyGamma[2, a]^2*PolyGamma[3, a]*PolyGamma[5, a] +
     2572970400*PolyGamma[1, a]*PolyGamma[2, a]^2*PolyGamma[3, a]*
      PolyGamma[5, a] + 321621300*PolyGamma[0, a]^4*PolyGamma[3, a]^2*
      PolyGamma[5, a] + 1929727800*PolyGamma[0, a]^2*PolyGamma[1, a]*
      PolyGamma[3, a]^2*PolyGamma[5, a] + 964863900*PolyGamma[1, a]^2*
      PolyGamma[3, a]^2*PolyGamma[5, a] + 1286485200*PolyGamma[0, a]*
      PolyGamma[2, a]*PolyGamma[3, a]^2*PolyGamma[5, a] +
     107207100*PolyGamma[3, a]^3*PolyGamma[5, a] + 14702688*PolyGamma[0, a]^7*
      PolyGamma[4, a]*PolyGamma[5, a] + 308756448*PolyGamma[0, a]^5*
      PolyGamma[1, a]*PolyGamma[4, a]*PolyGamma[5, a] +
     1543782240*PolyGamma[0, a]^3*PolyGamma[1, a]^2*PolyGamma[4, a]*
      PolyGamma[5, a] + 1543782240*PolyGamma[0, a]*PolyGamma[1, a]^3*
      PolyGamma[4, a]*PolyGamma[5, a] + 514594080*PolyGamma[0, a]^4*
      PolyGamma[2, a]*PolyGamma[4, a]*PolyGamma[5, a] +
     3087564480*PolyGamma[0, a]^2*PolyGamma[1, a]*PolyGamma[2, a]*
      PolyGamma[4, a]*PolyGamma[5, a] + 1543782240*PolyGamma[1, a]^2*
      PolyGamma[2, a]*PolyGamma[4, a]*PolyGamma[5, a] +
     1029188160*PolyGamma[0, a]*PolyGamma[2, a]^2*PolyGamma[4, a]*
      PolyGamma[5, a] + 514594080*PolyGamma[0, a]^3*PolyGamma[3, a]*
      PolyGamma[4, a]*PolyGamma[5, a] + 1543782240*PolyGamma[0, a]*
      PolyGamma[1, a]*PolyGamma[3, a]*PolyGamma[4, a]*PolyGamma[5, a] +
     514594080*PolyGamma[2, a]*PolyGamma[3, a]*PolyGamma[4, a]*
      PolyGamma[5, a] + 154378224*PolyGamma[0, a]^2*PolyGamma[4, a]^2*
      PolyGamma[5, a] + 154378224*PolyGamma[1, a]*PolyGamma[4, a]^2*
      PolyGamma[5, a] + 8576568*PolyGamma[0, a]^6*PolyGamma[5, a]^2 +
     128648520*PolyGamma[0, a]^4*PolyGamma[1, a]*PolyGamma[5, a]^2 +
     385945560*PolyGamma[0, a]^2*PolyGamma[1, a]^2*PolyGamma[5, a]^2 +
     128648520*PolyGamma[1, a]^3*PolyGamma[5, a]^2 +
     171531360*PolyGamma[0, a]^3*PolyGamma[2, a]*PolyGamma[5, a]^2 +
     514594080*PolyGamma[0, a]*PolyGamma[1, a]*PolyGamma[2, a]*
      PolyGamma[5, a]^2 + 85765680*PolyGamma[2, a]^2*PolyGamma[5, a]^2 +
     128648520*PolyGamma[0, a]^2*PolyGamma[3, a]*PolyGamma[5, a]^2 +
     128648520*PolyGamma[1, a]*PolyGamma[3, a]*PolyGamma[5, a]^2 +
     51459408*PolyGamma[0, a]*PolyGamma[4, a]*PolyGamma[5, a]^2 +
     2858856*PolyGamma[5, a]^3 + 31824*PolyGamma[0, a]^11*PolyGamma[6, a] +
     1750320*PolyGamma[0, a]^9*PolyGamma[1, a]*PolyGamma[6, a] +
     31505760*PolyGamma[0, a]^7*PolyGamma[1, a]^2*PolyGamma[6, a] +
     220540320*PolyGamma[0, a]^5*PolyGamma[1, a]^3*PolyGamma[6, a] +
     551350800*PolyGamma[0, a]^3*PolyGamma[1, a]^4*PolyGamma[6, a] +
     330810480*PolyGamma[0, a]*PolyGamma[1, a]^5*PolyGamma[6, a] +
     5250960*PolyGamma[0, a]^8*PolyGamma[2, a]*PolyGamma[6, a] +
     147026880*PolyGamma[0, a]^6*PolyGamma[1, a]*PolyGamma[2, a]*
      PolyGamma[6, a] + 1102701600*PolyGamma[0, a]^4*PolyGamma[1, a]^2*
      PolyGamma[2, a]*PolyGamma[6, a] + 2205403200*PolyGamma[0, a]^2*
      PolyGamma[1, a]^3*PolyGamma[2, a]*PolyGamma[6, a] +
     551350800*PolyGamma[1, a]^4*PolyGamma[2, a]*PolyGamma[6, a] +
     147026880*PolyGamma[0, a]^5*PolyGamma[2, a]^2*PolyGamma[6, a] +
     1470268800*PolyGamma[0, a]^3*PolyGamma[1, a]*PolyGamma[2, a]^2*
      PolyGamma[6, a] + 2205403200*PolyGamma[0, a]*PolyGamma[1, a]^2*
      PolyGamma[2, a]^2*PolyGamma[6, a] + 490089600*PolyGamma[0, a]^2*
      PolyGamma[2, a]^3*PolyGamma[6, a] + 490089600*PolyGamma[1, a]*
      PolyGamma[2, a]^3*PolyGamma[6, a] + 10501920*PolyGamma[0, a]^7*
      PolyGamma[3, a]*PolyGamma[6, a] + 220540320*PolyGamma[0, a]^5*
      PolyGamma[1, a]*PolyGamma[3, a]*PolyGamma[6, a] +
     1102701600*PolyGamma[0, a]^3*PolyGamma[1, a]^2*PolyGamma[3, a]*
      PolyGamma[6, a] + 1102701600*PolyGamma[0, a]*PolyGamma[1, a]^3*
      PolyGamma[3, a]*PolyGamma[6, a] + 367567200*PolyGamma[0, a]^4*
      PolyGamma[2, a]*PolyGamma[3, a]*PolyGamma[6, a] +
     2205403200*PolyGamma[0, a]^2*PolyGamma[1, a]*PolyGamma[2, a]*
      PolyGamma[3, a]*PolyGamma[6, a] + 1102701600*PolyGamma[1, a]^2*
      PolyGamma[2, a]*PolyGamma[3, a]*PolyGamma[6, a] +
     735134400*PolyGamma[0, a]*PolyGamma[2, a]^2*PolyGamma[3, a]*
      PolyGamma[6, a] + 183783600*PolyGamma[0, a]^3*PolyGamma[3, a]^2*
      PolyGamma[6, a] + 551350800*PolyGamma[0, a]*PolyGamma[1, a]*
      PolyGamma[3, a]^2*PolyGamma[6, a] + 183783600*PolyGamma[2, a]*
      PolyGamma[3, a]^2*PolyGamma[6, a] + 14702688*PolyGamma[0, a]^6*
      PolyGamma[4, a]*PolyGamma[6, a] + 220540320*PolyGamma[0, a]^4*
      PolyGamma[1, a]*PolyGamma[4, a]*PolyGamma[6, a] +
     661620960*PolyGamma[0, a]^2*PolyGamma[1, a]^2*PolyGamma[4, a]*
      PolyGamma[6, a] + 220540320*PolyGamma[1, a]^3*PolyGamma[4, a]*
      PolyGamma[6, a] + 294053760*PolyGamma[0, a]^3*PolyGamma[2, a]*
      PolyGamma[4, a]*PolyGamma[6, a] + 882161280*PolyGamma[0, a]*
      PolyGamma[1, a]*PolyGamma[2, a]*PolyGamma[4, a]*PolyGamma[6, a] +
     147026880*PolyGamma[2, a]^2*PolyGamma[4, a]*PolyGamma[6, a] +
     220540320*PolyGamma[0, a]^2*PolyGamma[3, a]*PolyGamma[4, a]*
      PolyGamma[6, a] + 220540320*PolyGamma[1, a]*PolyGamma[3, a]*
      PolyGamma[4, a]*PolyGamma[6, a] + 44108064*PolyGamma[0, a]*
      PolyGamma[4, a]^2*PolyGamma[6, a] + 14702688*PolyGamma[0, a]^5*
      PolyGamma[5, a]*PolyGamma[6, a] + 147026880*PolyGamma[0, a]^3*
      PolyGamma[1, a]*PolyGamma[5, a]*PolyGamma[6, a] +
     220540320*PolyGamma[0, a]*PolyGamma[1, a]^2*PolyGamma[5, a]*
      PolyGamma[6, a] + 147026880*PolyGamma[0, a]^2*PolyGamma[2, a]*
      PolyGamma[5, a]*PolyGamma[6, a] + 147026880*PolyGamma[1, a]*
      PolyGamma[2, a]*PolyGamma[5, a]*PolyGamma[6, a] +
     73513440*PolyGamma[0, a]*PolyGamma[3, a]*PolyGamma[5, a]*
      PolyGamma[6, a] + 14702688*PolyGamma[4, a]*PolyGamma[5, a]*
      PolyGamma[6, a] + 5250960*PolyGamma[0, a]^4*PolyGamma[6, a]^2 +
     31505760*PolyGamma[0, a]^2*PolyGamma[1, a]*PolyGamma[6, a]^2 +
     15752880*PolyGamma[1, a]^2*PolyGamma[6, a]^2 +
     21003840*PolyGamma[0, a]*PolyGamma[2, a]*PolyGamma[6, a]^2 +
     5250960*PolyGamma[3, a]*PolyGamma[6, a]^2 + 43758*PolyGamma[0, a]^10*
      PolyGamma[7, a] + 1969110*PolyGamma[0, a]^8*PolyGamma[1, a]*
      PolyGamma[7, a] + 27567540*PolyGamma[0, a]^6*PolyGamma[1, a]^2*
      PolyGamma[7, a] + 137837700*PolyGamma[0, a]^4*PolyGamma[1, a]^3*
      PolyGamma[7, a] + 206756550*PolyGamma[0, a]^2*PolyGamma[1, a]^4*
      PolyGamma[7, a] + 41351310*PolyGamma[1, a]^5*PolyGamma[7, a] +
     5250960*PolyGamma[0, a]^7*PolyGamma[2, a]*PolyGamma[7, a] +
     110270160*PolyGamma[0, a]^5*PolyGamma[1, a]*PolyGamma[2, a]*
      PolyGamma[7, a] + 551350800*PolyGamma[0, a]^3*PolyGamma[1, a]^2*
      PolyGamma[2, a]*PolyGamma[7, a] + 551350800*PolyGamma[0, a]*
      PolyGamma[1, a]^3*PolyGamma[2, a]*PolyGamma[7, a] +
     91891800*PolyGamma[0, a]^4*PolyGamma[2, a]^2*PolyGamma[7, a] +
     551350800*PolyGamma[0, a]^2*PolyGamma[1, a]*PolyGamma[2, a]^2*
      PolyGamma[7, a] + 275675400*PolyGamma[1, a]^2*PolyGamma[2, a]^2*
      PolyGamma[7, a] + 122522400*PolyGamma[0, a]*PolyGamma[2, a]^3*
      PolyGamma[7, a] + 9189180*PolyGamma[0, a]^6*PolyGamma[3, a]*
      PolyGamma[7, a] + 137837700*PolyGamma[0, a]^4*PolyGamma[1, a]*
      PolyGamma[3, a]*PolyGamma[7, a] + 413513100*PolyGamma[0, a]^2*
      PolyGamma[1, a]^2*PolyGamma[3, a]*PolyGamma[7, a] +
     137837700*PolyGamma[1, a]^3*PolyGamma[3, a]*PolyGamma[7, a] +
     183783600*PolyGamma[0, a]^3*PolyGamma[2, a]*PolyGamma[3, a]*
      PolyGamma[7, a] + 551350800*PolyGamma[0, a]*PolyGamma[1, a]*
      PolyGamma[2, a]*PolyGamma[3, a]*PolyGamma[7, a] +
     91891800*PolyGamma[2, a]^2*PolyGamma[3, a]*PolyGamma[7, a] +
     68918850*PolyGamma[0, a]^2*PolyGamma[3, a]^2*PolyGamma[7, a] +
     68918850*PolyGamma[1, a]*PolyGamma[3, a]^2*PolyGamma[7, a] +
     11027016*PolyGamma[0, a]^5*PolyGamma[4, a]*PolyGamma[7, a] +
     110270160*PolyGamma[0, a]^3*PolyGamma[1, a]*PolyGamma[4, a]*
      PolyGamma[7, a] + 165405240*PolyGamma[0, a]*PolyGamma[1, a]^2*
      PolyGamma[4, a]*PolyGamma[7, a] + 110270160*PolyGamma[0, a]^2*
      PolyGamma[2, a]*PolyGamma[4, a]*PolyGamma[7, a] +
     110270160*PolyGamma[1, a]*PolyGamma[2, a]*PolyGamma[4, a]*
      PolyGamma[7, a] + 55135080*PolyGamma[0, a]*PolyGamma[3, a]*
      PolyGamma[4, a]*PolyGamma[7, a] + 5513508*PolyGamma[4, a]^2*
      PolyGamma[7, a] + 9189180*PolyGamma[0, a]^4*PolyGamma[5, a]*
      PolyGamma[7, a] + 55135080*PolyGamma[0, a]^2*PolyGamma[1, a]*
      PolyGamma[5, a]*PolyGamma[7, a] + 27567540*PolyGamma[1, a]^2*
      PolyGamma[5, a]*PolyGamma[7, a] + 36756720*PolyGamma[0, a]*
      PolyGamma[2, a]*PolyGamma[5, a]*PolyGamma[7, a] +
     9189180*PolyGamma[3, a]*PolyGamma[5, a]*PolyGamma[7, a] +
     5250960*PolyGamma[0, a]^3*PolyGamma[6, a]*PolyGamma[7, a] +
     15752880*PolyGamma[0, a]*PolyGamma[1, a]*PolyGamma[6, a]*
      PolyGamma[7, a] + 5250960*PolyGamma[2, a]*PolyGamma[6, a]*
      PolyGamma[7, a] + 984555*PolyGamma[0, a]^2*PolyGamma[7, a]^2 +
     984555*PolyGamma[1, a]*PolyGamma[7, a]^2 + 48620*PolyGamma[0, a]^9*
      PolyGamma[8, a] + 1750320*PolyGamma[0, a]^7*PolyGamma[1, a]*
      PolyGamma[8, a] + 18378360*PolyGamma[0, a]^5*PolyGamma[1, a]^2*
      PolyGamma[8, a] + 61261200*PolyGamma[0, a]^3*PolyGamma[1, a]^3*
      PolyGamma[8, a] + 45945900*PolyGamma[0, a]*PolyGamma[1, a]^4*
      PolyGamma[8, a] + 4084080*PolyGamma[0, a]^6*PolyGamma[2, a]*
      PolyGamma[8, a] + 61261200*PolyGamma[0, a]^4*PolyGamma[1, a]*
      PolyGamma[2, a]*PolyGamma[8, a] + 183783600*PolyGamma[0, a]^2*
      PolyGamma[1, a]^2*PolyGamma[2, a]*PolyGamma[8, a] +
     61261200*PolyGamma[1, a]^3*PolyGamma[2, a]*PolyGamma[8, a] +
     40840800*PolyGamma[0, a]^3*PolyGamma[2, a]^2*PolyGamma[8, a] +
     122522400*PolyGamma[0, a]*PolyGamma[1, a]*PolyGamma[2, a]^2*
      PolyGamma[8, a] + 13613600*PolyGamma[2, a]^3*PolyGamma[8, a] +
     6126120*PolyGamma[0, a]^5*PolyGamma[3, a]*PolyGamma[8, a] +
     61261200*PolyGamma[0, a]^3*PolyGamma[1, a]*PolyGamma[3, a]*
      PolyGamma[8, a] + 91891800*PolyGamma[0, a]*PolyGamma[1, a]^2*
      PolyGamma[3, a]*PolyGamma[8, a] + 61261200*PolyGamma[0, a]^2*
      PolyGamma[2, a]*PolyGamma[3, a]*PolyGamma[8, a] +
     61261200*PolyGamma[1, a]*PolyGamma[2, a]*PolyGamma[3, a]*
      PolyGamma[8, a] + 15315300*PolyGamma[0, a]*PolyGamma[3, a]^2*
      PolyGamma[8, a] + 6126120*PolyGamma[0, a]^4*PolyGamma[4, a]*
      PolyGamma[8, a] + 36756720*PolyGamma[0, a]^2*PolyGamma[1, a]*
      PolyGamma[4, a]*PolyGamma[8, a] + 18378360*PolyGamma[1, a]^2*
      PolyGamma[4, a]*PolyGamma[8, a] + 24504480*PolyGamma[0, a]*
      PolyGamma[2, a]*PolyGamma[4, a]*PolyGamma[8, a] +
     6126120*PolyGamma[3, a]*PolyGamma[4, a]*PolyGamma[8, a] +
     4084080*PolyGamma[0, a]^3*PolyGamma[5, a]*PolyGamma[8, a] +
     12252240*PolyGamma[0, a]*PolyGamma[1, a]*PolyGamma[5, a]*
      PolyGamma[8, a] + 4084080*PolyGamma[2, a]*PolyGamma[5, a]*
      PolyGamma[8, a] + 1750320*PolyGamma[0, a]^2*PolyGamma[6, a]*
      PolyGamma[8, a] + 1750320*PolyGamma[1, a]*PolyGamma[6, a]*
      PolyGamma[8, a] + 437580*PolyGamma[0, a]*PolyGamma[7, a]*
      PolyGamma[8, a] + 24310*PolyGamma[8, a]^2 + 43758*PolyGamma[0, a]^8*
      PolyGamma[9, a] + 1225224*PolyGamma[0, a]^6*PolyGamma[1, a]*
      PolyGamma[9, a] + 9189180*PolyGamma[0, a]^4*PolyGamma[1, a]^2*
      PolyGamma[9, a] + 18378360*PolyGamma[0, a]^2*PolyGamma[1, a]^3*
      PolyGamma[9, a] + 4594590*PolyGamma[1, a]^4*PolyGamma[9, a] +
     2450448*PolyGamma[0, a]^5*PolyGamma[2, a]*PolyGamma[9, a] +
     24504480*PolyGamma[0, a]^3*PolyGamma[1, a]*PolyGamma[2, a]*
      PolyGamma[9, a] + 36756720*PolyGamma[0, a]*PolyGamma[1, a]^2*
      PolyGamma[2, a]*PolyGamma[9, a] + 12252240*PolyGamma[0, a]^2*
      PolyGamma[2, a]^2*PolyGamma[9, a] + 12252240*PolyGamma[1, a]*
      PolyGamma[2, a]^2*PolyGamma[9, a] + 3063060*PolyGamma[0, a]^4*
      PolyGamma[3, a]*PolyGamma[9, a] + 18378360*PolyGamma[0, a]^2*
      PolyGamma[1, a]*PolyGamma[3, a]*PolyGamma[9, a] +
     9189180*PolyGamma[1, a]^2*PolyGamma[3, a]*PolyGamma[9, a] +
     12252240*PolyGamma[0, a]*PolyGamma[2, a]*PolyGamma[3, a]*
      PolyGamma[9, a] + 1531530*PolyGamma[3, a]^2*PolyGamma[9, a] +
     2450448*PolyGamma[0, a]^3*PolyGamma[4, a]*PolyGamma[9, a] +
     7351344*PolyGamma[0, a]*PolyGamma[1, a]*PolyGamma[4, a]*
      PolyGamma[9, a] + 2450448*PolyGamma[2, a]*PolyGamma[4, a]*
      PolyGamma[9, a] + 1225224*PolyGamma[0, a]^2*PolyGamma[5, a]*
      PolyGamma[9, a] + 1225224*PolyGamma[1, a]*PolyGamma[5, a]*
      PolyGamma[9, a] + 350064*PolyGamma[0, a]*PolyGamma[6, a]*
      PolyGamma[9, a] + 43758*PolyGamma[7, a]*PolyGamma[9, a] +
     31824*PolyGamma[0, a]^7*PolyGamma[10, a] + 668304*PolyGamma[0, a]^5*
      PolyGamma[1, a]*PolyGamma[10, a] + 3341520*PolyGamma[0, a]^3*
      PolyGamma[1, a]^2*PolyGamma[10, a] + 3341520*PolyGamma[0, a]*
      PolyGamma[1, a]^3*PolyGamma[10, a] + 1113840*PolyGamma[0, a]^4*
      PolyGamma[2, a]*PolyGamma[10, a] + 6683040*PolyGamma[0, a]^2*
      PolyGamma[1, a]*PolyGamma[2, a]*PolyGamma[10, a] +
     3341520*PolyGamma[1, a]^2*PolyGamma[2, a]*PolyGamma[10, a] +
     2227680*PolyGamma[0, a]*PolyGamma[2, a]^2*PolyGamma[10, a] +
     1113840*PolyGamma[0, a]^3*PolyGamma[3, a]*PolyGamma[10, a] +
     3341520*PolyGamma[0, a]*PolyGamma[1, a]*PolyGamma[3, a]*
      PolyGamma[10, a] + 1113840*PolyGamma[2, a]*PolyGamma[3, a]*
      PolyGamma[10, a] + 668304*PolyGamma[0, a]^2*PolyGamma[4, a]*
      PolyGamma[10, a] + 668304*PolyGamma[1, a]*PolyGamma[4, a]*
      PolyGamma[10, a] + 222768*PolyGamma[0, a]*PolyGamma[5, a]*
      PolyGamma[10, a] + 31824*PolyGamma[6, a]*PolyGamma[10, a] +
     18564*PolyGamma[0, a]^6*PolyGamma[11, a] + 278460*PolyGamma[0, a]^4*
      PolyGamma[1, a]*PolyGamma[11, a] + 835380*PolyGamma[0, a]^2*
      PolyGamma[1, a]^2*PolyGamma[11, a] + 278460*PolyGamma[1, a]^3*
      PolyGamma[11, a] + 371280*PolyGamma[0, a]^3*PolyGamma[2, a]*
      PolyGamma[11, a] + 1113840*PolyGamma[0, a]*PolyGamma[1, a]*
      PolyGamma[2, a]*PolyGamma[11, a] + 185640*PolyGamma[2, a]^2*
      PolyGamma[11, a] + 278460*PolyGamma[0, a]^2*PolyGamma[3, a]*
      PolyGamma[11, a] + 278460*PolyGamma[1, a]*PolyGamma[3, a]*
      PolyGamma[11, a] + 111384*PolyGamma[0, a]*PolyGamma[4, a]*
      PolyGamma[11, a] + 18564*PolyGamma[5, a]*PolyGamma[11, a] +
     8568*PolyGamma[0, a]^5*PolyGamma[12, a] + 85680*PolyGamma[0, a]^3*
      PolyGamma[1, a]*PolyGamma[12, a] + 128520*PolyGamma[0, a]*
      PolyGamma[1, a]^2*PolyGamma[12, a] + 85680*PolyGamma[0, a]^2*
      PolyGamma[2, a]*PolyGamma[12, a] + 85680*PolyGamma[1, a]*
      PolyGamma[2, a]*PolyGamma[12, a] + 42840*PolyGamma[0, a]*
      PolyGamma[3, a]*PolyGamma[12, a] + 8568*PolyGamma[4, a]*
      PolyGamma[12, a] + 3060*PolyGamma[0, a]^4*PolyGamma[13, a] +
     18360*PolyGamma[0, a]^2*PolyGamma[1, a]*PolyGamma[13, a] +
     9180*PolyGamma[1, a]^2*PolyGamma[13, a] + 12240*PolyGamma[0, a]*
      PolyGamma[2, a]*PolyGamma[13, a] + 3060*PolyGamma[3, a]*
      PolyGamma[13, a] + 816*PolyGamma[0, a]^3*PolyGamma[14, a] +
     2448*PolyGamma[0, a]*PolyGamma[1, a]*PolyGamma[14, a] +
     816*PolyGamma[2, a]*PolyGamma[14, a] + 153*PolyGamma[0, a]^2*
      PolyGamma[15, a] + 153*PolyGamma[1, a]*PolyGamma[15, a] +
     18*PolyGamma[0, a]*PolyGamma[16, a] + PolyGamma[17, a]

MBexpGam[a_, 19] = PolyGamma[0, a]^19 + 171*PolyGamma[0, a]^17*
      PolyGamma[1, a] + 11628*PolyGamma[0, a]^15*PolyGamma[1, a]^2 +
     406980*PolyGamma[0, a]^13*PolyGamma[1, a]^3 + 7936110*PolyGamma[0, a]^11*
      PolyGamma[1, a]^4 + 87297210*PolyGamma[0, a]^9*PolyGamma[1, a]^5 +
     523783260*PolyGamma[0, a]^7*PolyGamma[1, a]^6 +
     1571349780*PolyGamma[0, a]^5*PolyGamma[1, a]^7 +
     1964187225*PolyGamma[0, a]^3*PolyGamma[1, a]^8 +
     654729075*PolyGamma[0, a]*PolyGamma[1, a]^9 + 969*PolyGamma[0, a]^16*
      PolyGamma[2, a] + 116280*PolyGamma[0, a]^14*PolyGamma[1, a]*
      PolyGamma[2, a] + 5290740*PolyGamma[0, a]^12*PolyGamma[1, a]^2*
      PolyGamma[2, a] + 116396280*PolyGamma[0, a]^10*PolyGamma[1, a]^3*
      PolyGamma[2, a] + 1309458150*PolyGamma[0, a]^8*PolyGamma[1, a]^4*
      PolyGamma[2, a] + 7332965640*PolyGamma[0, a]^6*PolyGamma[1, a]^5*
      PolyGamma[2, a] + 18332414100*PolyGamma[0, a]^4*PolyGamma[1, a]^6*
      PolyGamma[2, a] + 15713497800*PolyGamma[0, a]^2*PolyGamma[1, a]^7*
      PolyGamma[2, a] + 1964187225*PolyGamma[1, a]^8*PolyGamma[2, a] +
     271320*PolyGamma[0, a]^13*PolyGamma[2, a]^2 +
     21162960*PolyGamma[0, a]^11*PolyGamma[1, a]*PolyGamma[2, a]^2 +
     581981400*PolyGamma[0, a]^9*PolyGamma[1, a]^2*PolyGamma[2, a]^2 +
     6983776800*PolyGamma[0, a]^7*PolyGamma[1, a]^3*PolyGamma[2, a]^2 +
     36664828200*PolyGamma[0, a]^5*PolyGamma[1, a]^4*PolyGamma[2, a]^2 +
     73329656400*PolyGamma[0, a]^3*PolyGamma[1, a]^5*PolyGamma[2, a]^2 +
     36664828200*PolyGamma[0, a]*PolyGamma[1, a]^6*PolyGamma[2, a]^2 +
     25865840*PolyGamma[0, a]^10*PolyGamma[2, a]^3 +
     1163962800*PolyGamma[0, a]^8*PolyGamma[1, a]*PolyGamma[2, a]^3 +
     16295479200*PolyGamma[0, a]^6*PolyGamma[1, a]^2*PolyGamma[2, a]^3 +
     81477396000*PolyGamma[0, a]^4*PolyGamma[1, a]^3*PolyGamma[2, a]^3 +
     122216094000*PolyGamma[0, a]^2*PolyGamma[1, a]^4*PolyGamma[2, a]^3 +
     24443218800*PolyGamma[1, a]^5*PolyGamma[2, a]^3 +
     775975200*PolyGamma[0, a]^7*PolyGamma[2, a]^4 +
     16295479200*PolyGamma[0, a]^5*PolyGamma[1, a]*PolyGamma[2, a]^4 +
     81477396000*PolyGamma[0, a]^3*PolyGamma[1, a]^2*PolyGamma[2, a]^4 +
     81477396000*PolyGamma[0, a]*PolyGamma[1, a]^3*PolyGamma[2, a]^4 +
     5431826400*PolyGamma[0, a]^4*PolyGamma[2, a]^5 +
     32590958400*PolyGamma[0, a]^2*PolyGamma[1, a]*PolyGamma[2, a]^5 +
     16295479200*PolyGamma[1, a]^2*PolyGamma[2, a]^5 +
     3621217600*PolyGamma[0, a]*PolyGamma[2, a]^6 +
     3876*PolyGamma[0, a]^15*PolyGamma[3, a] + 406980*PolyGamma[0, a]^13*
      PolyGamma[1, a]*PolyGamma[3, a] + 15872220*PolyGamma[0, a]^11*
      PolyGamma[1, a]^2*PolyGamma[3, a] + 290990700*PolyGamma[0, a]^9*
      PolyGamma[1, a]^3*PolyGamma[3, a] + 2618916300*PolyGamma[0, a]^7*
      PolyGamma[1, a]^4*PolyGamma[3, a] + 10999448460*PolyGamma[0, a]^5*
      PolyGamma[1, a]^5*PolyGamma[3, a] + 18332414100*PolyGamma[0, a]^3*
      PolyGamma[1, a]^6*PolyGamma[3, a] + 7856748900*PolyGamma[0, a]*
      PolyGamma[1, a]^7*PolyGamma[3, a] + 1763580*PolyGamma[0, a]^12*
      PolyGamma[2, a]*PolyGamma[3, a] + 116396280*PolyGamma[0, a]^10*
      PolyGamma[1, a]*PolyGamma[2, a]*PolyGamma[3, a] +
     2618916300*PolyGamma[0, a]^8*PolyGamma[1, a]^2*PolyGamma[2, a]*
      PolyGamma[3, a] + 24443218800*PolyGamma[0, a]^6*PolyGamma[1, a]^3*
      PolyGamma[2, a]*PolyGamma[3, a] + 91662070500*PolyGamma[0, a]^4*
      PolyGamma[1, a]^4*PolyGamma[2, a]*PolyGamma[3, a] +
     109994484600*PolyGamma[0, a]^2*PolyGamma[1, a]^5*PolyGamma[2, a]*
      PolyGamma[3, a] + 18332414100*PolyGamma[1, a]^6*PolyGamma[2, a]*
      PolyGamma[3, a] + 193993800*PolyGamma[0, a]^9*PolyGamma[2, a]^2*
      PolyGamma[3, a] + 6983776800*PolyGamma[0, a]^7*PolyGamma[1, a]*
      PolyGamma[2, a]^2*PolyGamma[3, a] + 73329656400*PolyGamma[0, a]^5*
      PolyGamma[1, a]^2*PolyGamma[2, a]^2*PolyGamma[3, a] +
     244432188000*PolyGamma[0, a]^3*PolyGamma[1, a]^3*PolyGamma[2, a]^2*
      PolyGamma[3, a] + 183324141000*PolyGamma[0, a]*PolyGamma[1, a]^4*
      PolyGamma[2, a]^2*PolyGamma[3, a] + 5431826400*PolyGamma[0, a]^6*
      PolyGamma[2, a]^3*PolyGamma[3, a] + 81477396000*PolyGamma[0, a]^4*
      PolyGamma[1, a]*PolyGamma[2, a]^3*PolyGamma[3, a] +
     244432188000*PolyGamma[0, a]^2*PolyGamma[1, a]^2*PolyGamma[2, a]^3*
      PolyGamma[3, a] + 81477396000*PolyGamma[1, a]^3*PolyGamma[2, a]^3*
      PolyGamma[3, a] + 27159132000*PolyGamma[0, a]^3*PolyGamma[2, a]^4*
      PolyGamma[3, a] + 81477396000*PolyGamma[0, a]*PolyGamma[1, a]*
      PolyGamma[2, a]^4*PolyGamma[3, a] + 5431826400*PolyGamma[2, a]^5*
      PolyGamma[3, a] + 2645370*PolyGamma[0, a]^11*PolyGamma[3, a]^2 +
     145495350*PolyGamma[0, a]^9*PolyGamma[1, a]*PolyGamma[3, a]^2 +
     2618916300*PolyGamma[0, a]^7*PolyGamma[1, a]^2*PolyGamma[3, a]^2 +
     18332414100*PolyGamma[0, a]^5*PolyGamma[1, a]^3*PolyGamma[3, a]^2 +
     45831035250*PolyGamma[0, a]^3*PolyGamma[1, a]^4*PolyGamma[3, a]^2 +
     27498621150*PolyGamma[0, a]*PolyGamma[1, a]^5*PolyGamma[3, a]^2 +
     436486050*PolyGamma[0, a]^8*PolyGamma[2, a]*PolyGamma[3, a]^2 +
     12221609400*PolyGamma[0, a]^6*PolyGamma[1, a]*PolyGamma[2, a]*
      PolyGamma[3, a]^2 + 91662070500*PolyGamma[0, a]^4*PolyGamma[1, a]^2*
      PolyGamma[2, a]*PolyGamma[3, a]^2 + 183324141000*PolyGamma[0, a]^2*
      PolyGamma[1, a]^3*PolyGamma[2, a]*PolyGamma[3, a]^2 +
     45831035250*PolyGamma[1, a]^4*PolyGamma[2, a]*PolyGamma[3, a]^2 +
     12221609400*PolyGamma[0, a]^5*PolyGamma[2, a]^2*PolyGamma[3, a]^2 +
     122216094000*PolyGamma[0, a]^3*PolyGamma[1, a]*PolyGamma[2, a]^2*
      PolyGamma[3, a]^2 + 183324141000*PolyGamma[0, a]*PolyGamma[1, a]^2*
      PolyGamma[2, a]^2*PolyGamma[3, a]^2 + 40738698000*PolyGamma[0, a]^2*
      PolyGamma[2, a]^3*PolyGamma[3, a]^2 + 40738698000*PolyGamma[1, a]*
      PolyGamma[2, a]^3*PolyGamma[3, a]^2 + 290990700*PolyGamma[0, a]^7*
      PolyGamma[3, a]^3 + 6110804700*PolyGamma[0, a]^5*PolyGamma[1, a]*
      PolyGamma[3, a]^3 + 30554023500*PolyGamma[0, a]^3*PolyGamma[1, a]^2*
      PolyGamma[3, a]^3 + 30554023500*PolyGamma[0, a]*PolyGamma[1, a]^3*
      PolyGamma[3, a]^3 + 10184674500*PolyGamma[0, a]^4*PolyGamma[2, a]*
      PolyGamma[3, a]^3 + 61108047000*PolyGamma[0, a]^2*PolyGamma[1, a]*
      PolyGamma[2, a]*PolyGamma[3, a]^3 + 30554023500*PolyGamma[1, a]^2*
      PolyGamma[2, a]*PolyGamma[3, a]^3 + 20369349000*PolyGamma[0, a]*
      PolyGamma[2, a]^2*PolyGamma[3, a]^3 + 2546168625*PolyGamma[0, a]^3*
      PolyGamma[3, a]^4 + 7638505875*PolyGamma[0, a]*PolyGamma[1, a]*
      PolyGamma[3, a]^4 + 2546168625*PolyGamma[2, a]*PolyGamma[3, a]^4 +
     11628*PolyGamma[0, a]^14*PolyGamma[4, a] + 1058148*PolyGamma[0, a]^12*
      PolyGamma[1, a]*PolyGamma[4, a] + 34918884*PolyGamma[0, a]^10*
      PolyGamma[1, a]^2*PolyGamma[4, a] + 523783260*PolyGamma[0, a]^8*
      PolyGamma[1, a]^3*PolyGamma[4, a] + 3666482820*PolyGamma[0, a]^6*
      PolyGamma[1, a]^4*PolyGamma[4, a] + 10999448460*PolyGamma[0, a]^4*
      PolyGamma[1, a]^5*PolyGamma[4, a] + 10999448460*PolyGamma[0, a]^2*
      PolyGamma[1, a]^6*PolyGamma[4, a] + 1571349780*PolyGamma[1, a]^7*
      PolyGamma[4, a] + 4232592*PolyGamma[0, a]^11*PolyGamma[2, a]*
      PolyGamma[4, a] + 232792560*PolyGamma[0, a]^9*PolyGamma[1, a]*
      PolyGamma[2, a]*PolyGamma[4, a] + 4190266080*PolyGamma[0, a]^7*
      PolyGamma[1, a]^2*PolyGamma[2, a]*PolyGamma[4, a] +
     29331862560*PolyGamma[0, a]^5*PolyGamma[1, a]^3*PolyGamma[2, a]*
      PolyGamma[4, a] + 73329656400*PolyGamma[0, a]^3*PolyGamma[1, a]^4*
      PolyGamma[2, a]*PolyGamma[4, a] + 43997793840*PolyGamma[0, a]*
      PolyGamma[1, a]^5*PolyGamma[2, a]*PolyGamma[4, a] +
     349188840*PolyGamma[0, a]^8*PolyGamma[2, a]^2*PolyGamma[4, a] +
     9777287520*PolyGamma[0, a]^6*PolyGamma[1, a]*PolyGamma[2, a]^2*
      PolyGamma[4, a] + 73329656400*PolyGamma[0, a]^4*PolyGamma[1, a]^2*
      PolyGamma[2, a]^2*PolyGamma[4, a] + 146659312800*PolyGamma[0, a]^2*
      PolyGamma[1, a]^3*PolyGamma[2, a]^2*PolyGamma[4, a] +
     36664828200*PolyGamma[1, a]^4*PolyGamma[2, a]^2*PolyGamma[4, a] +
     6518191680*PolyGamma[0, a]^5*PolyGamma[2, a]^3*PolyGamma[4, a] +
     65181916800*PolyGamma[0, a]^3*PolyGamma[1, a]*PolyGamma[2, a]^3*
      PolyGamma[4, a] + 97772875200*PolyGamma[0, a]*PolyGamma[1, a]^2*
      PolyGamma[2, a]^3*PolyGamma[4, a] + 16295479200*PolyGamma[0, a]^2*
      PolyGamma[2, a]^4*PolyGamma[4, a] + 16295479200*PolyGamma[1, a]*
      PolyGamma[2, a]^4*PolyGamma[4, a] + 11639628*PolyGamma[0, a]^10*
      PolyGamma[3, a]*PolyGamma[4, a] + 523783260*PolyGamma[0, a]^8*
      PolyGamma[1, a]*PolyGamma[3, a]*PolyGamma[4, a] +
     7332965640*PolyGamma[0, a]^6*PolyGamma[1, a]^2*PolyGamma[3, a]*
      PolyGamma[4, a] + 36664828200*PolyGamma[0, a]^4*PolyGamma[1, a]^3*
      PolyGamma[3, a]*PolyGamma[4, a] + 54997242300*PolyGamma[0, a]^2*
      PolyGamma[1, a]^4*PolyGamma[3, a]*PolyGamma[4, a] +
     10999448460*PolyGamma[1, a]^5*PolyGamma[3, a]*PolyGamma[4, a] +
     1396755360*PolyGamma[0, a]^7*PolyGamma[2, a]*PolyGamma[3, a]*
      PolyGamma[4, a] + 29331862560*PolyGamma[0, a]^5*PolyGamma[1, a]*
      PolyGamma[2, a]*PolyGamma[3, a]*PolyGamma[4, a] +
     146659312800*PolyGamma[0, a]^3*PolyGamma[1, a]^2*PolyGamma[2, a]*
      PolyGamma[3, a]*PolyGamma[4, a] + 146659312800*PolyGamma[0, a]*
      PolyGamma[1, a]^3*PolyGamma[2, a]*PolyGamma[3, a]*PolyGamma[4, a] +
     24443218800*PolyGamma[0, a]^4*PolyGamma[2, a]^2*PolyGamma[3, a]*
      PolyGamma[4, a] + 146659312800*PolyGamma[0, a]^2*PolyGamma[1, a]*
      PolyGamma[2, a]^2*PolyGamma[3, a]*PolyGamma[4, a] +
     73329656400*PolyGamma[1, a]^2*PolyGamma[2, a]^2*PolyGamma[3, a]*
      PolyGamma[4, a] + 32590958400*PolyGamma[0, a]*PolyGamma[2, a]^3*
      PolyGamma[3, a]*PolyGamma[4, a] + 1222160940*PolyGamma[0, a]^6*
      PolyGamma[3, a]^2*PolyGamma[4, a] + 18332414100*PolyGamma[0, a]^4*
      PolyGamma[1, a]*PolyGamma[3, a]^2*PolyGamma[4, a] +
     54997242300*PolyGamma[0, a]^2*PolyGamma[1, a]^2*PolyGamma[3, a]^2*
      PolyGamma[4, a] + 18332414100*PolyGamma[1, a]^3*PolyGamma[3, a]^2*
      PolyGamma[4, a] + 24443218800*PolyGamma[0, a]^3*PolyGamma[2, a]*
      PolyGamma[3, a]^2*PolyGamma[4, a] + 73329656400*PolyGamma[0, a]*
      PolyGamma[1, a]*PolyGamma[2, a]*PolyGamma[3, a]^2*PolyGamma[4, a] +
     12221609400*PolyGamma[2, a]^2*PolyGamma[3, a]^2*PolyGamma[4, a] +
     6110804700*PolyGamma[0, a]^2*PolyGamma[3, a]^3*PolyGamma[4, a] +
     6110804700*PolyGamma[1, a]*PolyGamma[3, a]^3*PolyGamma[4, a] +
     11639628*PolyGamma[0, a]^9*PolyGamma[4, a]^2 +
     419026608*PolyGamma[0, a]^7*PolyGamma[1, a]*PolyGamma[4, a]^2 +
     4399779384*PolyGamma[0, a]^5*PolyGamma[1, a]^2*PolyGamma[4, a]^2 +
     14665931280*PolyGamma[0, a]^3*PolyGamma[1, a]^3*PolyGamma[4, a]^2 +
     10999448460*PolyGamma[0, a]*PolyGamma[1, a]^4*PolyGamma[4, a]^2 +
     977728752*PolyGamma[0, a]^6*PolyGamma[2, a]*PolyGamma[4, a]^2 +
     14665931280*PolyGamma[0, a]^4*PolyGamma[1, a]*PolyGamma[2, a]*
      PolyGamma[4, a]^2 + 43997793840*PolyGamma[0, a]^2*PolyGamma[1, a]^2*
      PolyGamma[2, a]*PolyGamma[4, a]^2 + 14665931280*PolyGamma[1, a]^3*
      PolyGamma[2, a]*PolyGamma[4, a]^2 + 9777287520*PolyGamma[0, a]^3*
      PolyGamma[2, a]^2*PolyGamma[4, a]^2 + 29331862560*PolyGamma[0, a]*
      PolyGamma[1, a]*PolyGamma[2, a]^2*PolyGamma[4, a]^2 +
     3259095840*PolyGamma[2, a]^3*PolyGamma[4, a]^2 +
     1466593128*PolyGamma[0, a]^5*PolyGamma[3, a]*PolyGamma[4, a]^2 +
     14665931280*PolyGamma[0, a]^3*PolyGamma[1, a]*PolyGamma[3, a]*
      PolyGamma[4, a]^2 + 21998896920*PolyGamma[0, a]*PolyGamma[1, a]^2*
      PolyGamma[3, a]*PolyGamma[4, a]^2 + 14665931280*PolyGamma[0, a]^2*
      PolyGamma[2, a]*PolyGamma[3, a]*PolyGamma[4, a]^2 +
     14665931280*PolyGamma[1, a]*PolyGamma[2, a]*PolyGamma[3, a]*
      PolyGamma[4, a]^2 + 3666482820*PolyGamma[0, a]*PolyGamma[3, a]^2*
      PolyGamma[4, a]^2 + 488864376*PolyGamma[0, a]^4*PolyGamma[4, a]^3 +
     2933186256*PolyGamma[0, a]^2*PolyGamma[1, a]*PolyGamma[4, a]^3 +
     1466593128*PolyGamma[1, a]^2*PolyGamma[4, a]^3 +
     1955457504*PolyGamma[0, a]*PolyGamma[2, a]*PolyGamma[4, a]^3 +
     488864376*PolyGamma[3, a]*PolyGamma[4, a]^3 + 27132*PolyGamma[0, a]^13*
      PolyGamma[5, a] + 2116296*PolyGamma[0, a]^11*PolyGamma[1, a]*
      PolyGamma[5, a] + 58198140*PolyGamma[0, a]^9*PolyGamma[1, a]^2*
      PolyGamma[5, a] + 698377680*PolyGamma[0, a]^7*PolyGamma[1, a]^3*
      PolyGamma[5, a] + 3666482820*PolyGamma[0, a]^5*PolyGamma[1, a]^4*
      PolyGamma[5, a] + 7332965640*PolyGamma[0, a]^3*PolyGamma[1, a]^5*
      PolyGamma[5, a] + 3666482820*PolyGamma[0, a]*PolyGamma[1, a]^6*
      PolyGamma[5, a] + 7759752*PolyGamma[0, a]^10*PolyGamma[2, a]*
      PolyGamma[5, a] + 349188840*PolyGamma[0, a]^8*PolyGamma[1, a]*
      PolyGamma[2, a]*PolyGamma[5, a] + 4888643760*PolyGamma[0, a]^6*
      PolyGamma[1, a]^2*PolyGamma[2, a]*PolyGamma[5, a] +
     24443218800*PolyGamma[0, a]^4*PolyGamma[1, a]^3*PolyGamma[2, a]*
      PolyGamma[5, a] + 36664828200*PolyGamma[0, a]^2*PolyGamma[1, a]^4*
      PolyGamma[2, a]*PolyGamma[5, a] + 7332965640*PolyGamma[1, a]^5*
      PolyGamma[2, a]*PolyGamma[5, a] + 465585120*PolyGamma[0, a]^7*
      PolyGamma[2, a]^2*PolyGamma[5, a] + 9777287520*PolyGamma[0, a]^5*
      PolyGamma[1, a]*PolyGamma[2, a]^2*PolyGamma[5, a] +
     48886437600*PolyGamma[0, a]^3*PolyGamma[1, a]^2*PolyGamma[2, a]^2*
      PolyGamma[5, a] + 48886437600*PolyGamma[0, a]*PolyGamma[1, a]^3*
      PolyGamma[2, a]^2*PolyGamma[5, a] + 5431826400*PolyGamma[0, a]^4*
      PolyGamma[2, a]^3*PolyGamma[5, a] + 32590958400*PolyGamma[0, a]^2*
      PolyGamma[1, a]*PolyGamma[2, a]^3*PolyGamma[5, a] +
     16295479200*PolyGamma[1, a]^2*PolyGamma[2, a]^3*PolyGamma[5, a] +
     5431826400*PolyGamma[0, a]*PolyGamma[2, a]^4*PolyGamma[5, a] +
     19399380*PolyGamma[0, a]^9*PolyGamma[3, a]*PolyGamma[5, a] +
     698377680*PolyGamma[0, a]^7*PolyGamma[1, a]*PolyGamma[3, a]*
      PolyGamma[5, a] + 7332965640*PolyGamma[0, a]^5*PolyGamma[1, a]^2*
      PolyGamma[3, a]*PolyGamma[5, a] + 24443218800*PolyGamma[0, a]^3*
      PolyGamma[1, a]^3*PolyGamma[3, a]*PolyGamma[5, a] +
     18332414100*PolyGamma[0, a]*PolyGamma[1, a]^4*PolyGamma[3, a]*
      PolyGamma[5, a] + 1629547920*PolyGamma[0, a]^6*PolyGamma[2, a]*
      PolyGamma[3, a]*PolyGamma[5, a] + 24443218800*PolyGamma[0, a]^4*
      PolyGamma[1, a]*PolyGamma[2, a]*PolyGamma[3, a]*PolyGamma[5, a] +
     73329656400*PolyGamma[0, a]^2*PolyGamma[1, a]^2*PolyGamma[2, a]*
      PolyGamma[3, a]*PolyGamma[5, a] + 24443218800*PolyGamma[1, a]^3*
      PolyGamma[2, a]*PolyGamma[3, a]*PolyGamma[5, a] +
     16295479200*PolyGamma[0, a]^3*PolyGamma[2, a]^2*PolyGamma[3, a]*
      PolyGamma[5, a] + 48886437600*PolyGamma[0, a]*PolyGamma[1, a]*
      PolyGamma[2, a]^2*PolyGamma[3, a]*PolyGamma[5, a] +
     5431826400*PolyGamma[2, a]^3*PolyGamma[3, a]*PolyGamma[5, a] +
     1222160940*PolyGamma[0, a]^5*PolyGamma[3, a]^2*PolyGamma[5, a] +
     12221609400*PolyGamma[0, a]^3*PolyGamma[1, a]*PolyGamma[3, a]^2*
      PolyGamma[5, a] + 18332414100*PolyGamma[0, a]*PolyGamma[1, a]^2*
      PolyGamma[3, a]^2*PolyGamma[5, a] + 12221609400*PolyGamma[0, a]^2*
      PolyGamma[2, a]*PolyGamma[3, a]^2*PolyGamma[5, a] +
     12221609400*PolyGamma[1, a]*PolyGamma[2, a]*PolyGamma[3, a]^2*
      PolyGamma[5, a] + 2036934900*PolyGamma[0, a]*PolyGamma[3, a]^3*
      PolyGamma[5, a] + 34918884*PolyGamma[0, a]^8*PolyGamma[4, a]*
      PolyGamma[5, a] + 977728752*PolyGamma[0, a]^6*PolyGamma[1, a]*
      PolyGamma[4, a]*PolyGamma[5, a] + 7332965640*PolyGamma[0, a]^4*
      PolyGamma[1, a]^2*PolyGamma[4, a]*PolyGamma[5, a] +
     14665931280*PolyGamma[0, a]^2*PolyGamma[1, a]^3*PolyGamma[4, a]*
      PolyGamma[5, a] + 3666482820*PolyGamma[1, a]^4*PolyGamma[4, a]*
      PolyGamma[5, a] + 1955457504*PolyGamma[0, a]^5*PolyGamma[2, a]*
      PolyGamma[4, a]*PolyGamma[5, a] + 19554575040*PolyGamma[0, a]^3*
      PolyGamma[1, a]*PolyGamma[2, a]*PolyGamma[4, a]*PolyGamma[5, a] +
     29331862560*PolyGamma[0, a]*PolyGamma[1, a]^2*PolyGamma[2, a]*
      PolyGamma[4, a]*PolyGamma[5, a] + 9777287520*PolyGamma[0, a]^2*
      PolyGamma[2, a]^2*PolyGamma[4, a]*PolyGamma[5, a] +
     9777287520*PolyGamma[1, a]*PolyGamma[2, a]^2*PolyGamma[4, a]*
      PolyGamma[5, a] + 2444321880*PolyGamma[0, a]^4*PolyGamma[3, a]*
      PolyGamma[4, a]*PolyGamma[5, a] + 14665931280*PolyGamma[0, a]^2*
      PolyGamma[1, a]*PolyGamma[3, a]*PolyGamma[4, a]*PolyGamma[5, a] +
     7332965640*PolyGamma[1, a]^2*PolyGamma[3, a]*PolyGamma[4, a]*
      PolyGamma[5, a] + 9777287520*PolyGamma[0, a]*PolyGamma[2, a]*
      PolyGamma[3, a]*PolyGamma[4, a]*PolyGamma[5, a] +
     1222160940*PolyGamma[3, a]^2*PolyGamma[4, a]*PolyGamma[5, a] +
     977728752*PolyGamma[0, a]^3*PolyGamma[4, a]^2*PolyGamma[5, a] +
     2933186256*PolyGamma[0, a]*PolyGamma[1, a]*PolyGamma[4, a]^2*
      PolyGamma[5, a] + 977728752*PolyGamma[2, a]*PolyGamma[4, a]^2*
      PolyGamma[5, a] + 23279256*PolyGamma[0, a]^7*PolyGamma[5, a]^2 +
     488864376*PolyGamma[0, a]^5*PolyGamma[1, a]*PolyGamma[5, a]^2 +
     2444321880*PolyGamma[0, a]^3*PolyGamma[1, a]^2*PolyGamma[5, a]^2 +
     2444321880*PolyGamma[0, a]*PolyGamma[1, a]^3*PolyGamma[5, a]^2 +
     814773960*PolyGamma[0, a]^4*PolyGamma[2, a]*PolyGamma[5, a]^2 +
     4888643760*PolyGamma[0, a]^2*PolyGamma[1, a]*PolyGamma[2, a]*
      PolyGamma[5, a]^2 + 2444321880*PolyGamma[1, a]^2*PolyGamma[2, a]*
      PolyGamma[5, a]^2 + 1629547920*PolyGamma[0, a]*PolyGamma[2, a]^2*
      PolyGamma[5, a]^2 + 814773960*PolyGamma[0, a]^3*PolyGamma[3, a]*
      PolyGamma[5, a]^2 + 2444321880*PolyGamma[0, a]*PolyGamma[1, a]*
      PolyGamma[3, a]*PolyGamma[5, a]^2 + 814773960*PolyGamma[2, a]*
      PolyGamma[3, a]*PolyGamma[5, a]^2 + 488864376*PolyGamma[0, a]^2*
      PolyGamma[4, a]*PolyGamma[5, a]^2 + 488864376*PolyGamma[1, a]*
      PolyGamma[4, a]*PolyGamma[5, a]^2 + 54318264*PolyGamma[0, a]*
      PolyGamma[5, a]^3 + 50388*PolyGamma[0, a]^12*PolyGamma[6, a] +
     3325608*PolyGamma[0, a]^10*PolyGamma[1, a]*PolyGamma[6, a] +
     74826180*PolyGamma[0, a]^8*PolyGamma[1, a]^2*PolyGamma[6, a] +
     698377680*PolyGamma[0, a]^6*PolyGamma[1, a]^3*PolyGamma[6, a] +
     2618916300*PolyGamma[0, a]^4*PolyGamma[1, a]^4*PolyGamma[6, a] +
     3142699560*PolyGamma[0, a]^2*PolyGamma[1, a]^5*PolyGamma[6, a] +
     523783260*PolyGamma[1, a]^6*PolyGamma[6, a] + 11085360*PolyGamma[0, a]^9*
      PolyGamma[2, a]*PolyGamma[6, a] + 399072960*PolyGamma[0, a]^7*
      PolyGamma[1, a]*PolyGamma[2, a]*PolyGamma[6, a] +
     4190266080*PolyGamma[0, a]^5*PolyGamma[1, a]^2*PolyGamma[2, a]*
      PolyGamma[6, a] + 13967553600*PolyGamma[0, a]^3*PolyGamma[1, a]^3*
      PolyGamma[2, a]*PolyGamma[6, a] + 10475665200*PolyGamma[0, a]*
      PolyGamma[1, a]^4*PolyGamma[2, a]*PolyGamma[6, a] +
     465585120*PolyGamma[0, a]^6*PolyGamma[2, a]^2*PolyGamma[6, a] +
     6983776800*PolyGamma[0, a]^4*PolyGamma[1, a]*PolyGamma[2, a]^2*
      PolyGamma[6, a] + 20951330400*PolyGamma[0, a]^2*PolyGamma[1, a]^2*
      PolyGamma[2, a]^2*PolyGamma[6, a] + 6983776800*PolyGamma[1, a]^3*
      PolyGamma[2, a]^2*PolyGamma[6, a] + 3103900800*PolyGamma[0, a]^3*
      PolyGamma[2, a]^3*PolyGamma[6, a] + 9311702400*PolyGamma[0, a]*
      PolyGamma[1, a]*PolyGamma[2, a]^3*PolyGamma[6, a] +
     775975200*PolyGamma[2, a]^4*PolyGamma[6, a] + 24942060*PolyGamma[0, a]^8*
      PolyGamma[3, a]*PolyGamma[6, a] + 698377680*PolyGamma[0, a]^6*
      PolyGamma[1, a]*PolyGamma[3, a]*PolyGamma[6, a] +
     5237832600*PolyGamma[0, a]^4*PolyGamma[1, a]^2*PolyGamma[3, a]*
      PolyGamma[6, a] + 10475665200*PolyGamma[0, a]^2*PolyGamma[1, a]^3*
      PolyGamma[3, a]*PolyGamma[6, a] + 2618916300*PolyGamma[1, a]^4*
      PolyGamma[3, a]*PolyGamma[6, a] + 1396755360*PolyGamma[0, a]^5*
      PolyGamma[2, a]*PolyGamma[3, a]*PolyGamma[6, a] +
     13967553600*PolyGamma[0, a]^3*PolyGamma[1, a]*PolyGamma[2, a]*
      PolyGamma[3, a]*PolyGamma[6, a] + 20951330400*PolyGamma[0, a]*
      PolyGamma[1, a]^2*PolyGamma[2, a]*PolyGamma[3, a]*PolyGamma[6, a] +
     6983776800*PolyGamma[0, a]^2*PolyGamma[2, a]^2*PolyGamma[3, a]*
      PolyGamma[6, a] + 6983776800*PolyGamma[1, a]*PolyGamma[2, a]^2*
      PolyGamma[3, a]*PolyGamma[6, a] + 872972100*PolyGamma[0, a]^4*
      PolyGamma[3, a]^2*PolyGamma[6, a] + 5237832600*PolyGamma[0, a]^2*
      PolyGamma[1, a]*PolyGamma[3, a]^2*PolyGamma[6, a] +
     2618916300*PolyGamma[1, a]^2*PolyGamma[3, a]^2*PolyGamma[6, a] +
     3491888400*PolyGamma[0, a]*PolyGamma[2, a]*PolyGamma[3, a]^2*
      PolyGamma[6, a] + 290990700*PolyGamma[3, a]^3*PolyGamma[6, a] +
     39907296*PolyGamma[0, a]^7*PolyGamma[4, a]*PolyGamma[6, a] +
     838053216*PolyGamma[0, a]^5*PolyGamma[1, a]*PolyGamma[4, a]*
      PolyGamma[6, a] + 4190266080*PolyGamma[0, a]^3*PolyGamma[1, a]^2*
      PolyGamma[4, a]*PolyGamma[6, a] + 4190266080*PolyGamma[0, a]*
      PolyGamma[1, a]^3*PolyGamma[4, a]*PolyGamma[6, a] +
     1396755360*PolyGamma[0, a]^4*PolyGamma[2, a]*PolyGamma[4, a]*
      PolyGamma[6, a] + 8380532160*PolyGamma[0, a]^2*PolyGamma[1, a]*
      PolyGamma[2, a]*PolyGamma[4, a]*PolyGamma[6, a] +
     4190266080*PolyGamma[1, a]^2*PolyGamma[2, a]*PolyGamma[4, a]*
      PolyGamma[6, a] + 2793510720*PolyGamma[0, a]*PolyGamma[2, a]^2*
      PolyGamma[4, a]*PolyGamma[6, a] + 1396755360*PolyGamma[0, a]^3*
      PolyGamma[3, a]*PolyGamma[4, a]*PolyGamma[6, a] +
     4190266080*PolyGamma[0, a]*PolyGamma[1, a]*PolyGamma[3, a]*
      PolyGamma[4, a]*PolyGamma[6, a] + 1396755360*PolyGamma[2, a]*
      PolyGamma[3, a]*PolyGamma[4, a]*PolyGamma[6, a] +
     419026608*PolyGamma[0, a]^2*PolyGamma[4, a]^2*PolyGamma[6, a] +
     419026608*PolyGamma[1, a]*PolyGamma[4, a]^2*PolyGamma[6, a] +
     46558512*PolyGamma[0, a]^6*PolyGamma[5, a]*PolyGamma[6, a] +
     698377680*PolyGamma[0, a]^4*PolyGamma[1, a]*PolyGamma[5, a]*
      PolyGamma[6, a] + 2095133040*PolyGamma[0, a]^2*PolyGamma[1, a]^2*
      PolyGamma[5, a]*PolyGamma[6, a] + 698377680*PolyGamma[1, a]^3*
      PolyGamma[5, a]*PolyGamma[6, a] + 931170240*PolyGamma[0, a]^3*
      PolyGamma[2, a]*PolyGamma[5, a]*PolyGamma[6, a] +
     2793510720*PolyGamma[0, a]*PolyGamma[1, a]*PolyGamma[2, a]*
      PolyGamma[5, a]*PolyGamma[6, a] + 465585120*PolyGamma[2, a]^2*
      PolyGamma[5, a]*PolyGamma[6, a] + 698377680*PolyGamma[0, a]^2*
      PolyGamma[3, a]*PolyGamma[5, a]*PolyGamma[6, a] +
     698377680*PolyGamma[1, a]*PolyGamma[3, a]*PolyGamma[5, a]*
      PolyGamma[6, a] + 279351072*PolyGamma[0, a]*PolyGamma[4, a]*
      PolyGamma[5, a]*PolyGamma[6, a] + 23279256*PolyGamma[5, a]^2*
      PolyGamma[6, a] + 19953648*PolyGamma[0, a]^5*PolyGamma[6, a]^2 +
     199536480*PolyGamma[0, a]^3*PolyGamma[1, a]*PolyGamma[6, a]^2 +
     299304720*PolyGamma[0, a]*PolyGamma[1, a]^2*PolyGamma[6, a]^2 +
     199536480*PolyGamma[0, a]^2*PolyGamma[2, a]*PolyGamma[6, a]^2 +
     199536480*PolyGamma[1, a]*PolyGamma[2, a]*PolyGamma[6, a]^2 +
     99768240*PolyGamma[0, a]*PolyGamma[3, a]*PolyGamma[6, a]^2 +
     19953648*PolyGamma[4, a]*PolyGamma[6, a]^2 + 75582*PolyGamma[0, a]^11*
      PolyGamma[7, a] + 4157010*PolyGamma[0, a]^9*PolyGamma[1, a]*
      PolyGamma[7, a] + 74826180*PolyGamma[0, a]^7*PolyGamma[1, a]^2*
      PolyGamma[7, a] + 523783260*PolyGamma[0, a]^5*PolyGamma[1, a]^3*
      PolyGamma[7, a] + 1309458150*PolyGamma[0, a]^3*PolyGamma[1, a]^4*
      PolyGamma[7, a] + 785674890*PolyGamma[0, a]*PolyGamma[1, a]^5*
      PolyGamma[7, a] + 12471030*PolyGamma[0, a]^8*PolyGamma[2, a]*
      PolyGamma[7, a] + 349188840*PolyGamma[0, a]^6*PolyGamma[1, a]*
      PolyGamma[2, a]*PolyGamma[7, a] + 2618916300*PolyGamma[0, a]^4*
      PolyGamma[1, a]^2*PolyGamma[2, a]*PolyGamma[7, a] +
     5237832600*PolyGamma[0, a]^2*PolyGamma[1, a]^3*PolyGamma[2, a]*
      PolyGamma[7, a] + 1309458150*PolyGamma[1, a]^4*PolyGamma[2, a]*
      PolyGamma[7, a] + 349188840*PolyGamma[0, a]^5*PolyGamma[2, a]^2*
      PolyGamma[7, a] + 3491888400*PolyGamma[0, a]^3*PolyGamma[1, a]*
      PolyGamma[2, a]^2*PolyGamma[7, a] + 5237832600*PolyGamma[0, a]*
      PolyGamma[1, a]^2*PolyGamma[2, a]^2*PolyGamma[7, a] +
     1163962800*PolyGamma[0, a]^2*PolyGamma[2, a]^3*PolyGamma[7, a] +
     1163962800*PolyGamma[1, a]*PolyGamma[2, a]^3*PolyGamma[7, a] +
     24942060*PolyGamma[0, a]^7*PolyGamma[3, a]*PolyGamma[7, a] +
     523783260*PolyGamma[0, a]^5*PolyGamma[1, a]*PolyGamma[3, a]*
      PolyGamma[7, a] + 2618916300*PolyGamma[0, a]^3*PolyGamma[1, a]^2*
      PolyGamma[3, a]*PolyGamma[7, a] + 2618916300*PolyGamma[0, a]*
      PolyGamma[1, a]^3*PolyGamma[3, a]*PolyGamma[7, a] +
     872972100*PolyGamma[0, a]^4*PolyGamma[2, a]*PolyGamma[3, a]*
      PolyGamma[7, a] + 5237832600*PolyGamma[0, a]^2*PolyGamma[1, a]*
      PolyGamma[2, a]*PolyGamma[3, a]*PolyGamma[7, a] +
     2618916300*PolyGamma[1, a]^2*PolyGamma[2, a]*PolyGamma[3, a]*
      PolyGamma[7, a] + 1745944200*PolyGamma[0, a]*PolyGamma[2, a]^2*
      PolyGamma[3, a]*PolyGamma[7, a] + 436486050*PolyGamma[0, a]^3*
      PolyGamma[3, a]^2*PolyGamma[7, a] + 1309458150*PolyGamma[0, a]*
      PolyGamma[1, a]*PolyGamma[3, a]^2*PolyGamma[7, a] +
     436486050*PolyGamma[2, a]*PolyGamma[3, a]^2*PolyGamma[7, a] +
     34918884*PolyGamma[0, a]^6*PolyGamma[4, a]*PolyGamma[7, a] +
     523783260*PolyGamma[0, a]^4*PolyGamma[1, a]*PolyGamma[4, a]*
      PolyGamma[7, a] + 1571349780*PolyGamma[0, a]^2*PolyGamma[1, a]^2*
      PolyGamma[4, a]*PolyGamma[7, a] + 523783260*PolyGamma[1, a]^3*
      PolyGamma[4, a]*PolyGamma[7, a] + 698377680*PolyGamma[0, a]^3*
      PolyGamma[2, a]*PolyGamma[4, a]*PolyGamma[7, a] +
     2095133040*PolyGamma[0, a]*PolyGamma[1, a]*PolyGamma[2, a]*
      PolyGamma[4, a]*PolyGamma[7, a] + 349188840*PolyGamma[2, a]^2*
      PolyGamma[4, a]*PolyGamma[7, a] + 523783260*PolyGamma[0, a]^2*
      PolyGamma[3, a]*PolyGamma[4, a]*PolyGamma[7, a] +
     523783260*PolyGamma[1, a]*PolyGamma[3, a]*PolyGamma[4, a]*
      PolyGamma[7, a] + 104756652*PolyGamma[0, a]*PolyGamma[4, a]^2*
      PolyGamma[7, a] + 34918884*PolyGamma[0, a]^5*PolyGamma[5, a]*
      PolyGamma[7, a] + 349188840*PolyGamma[0, a]^3*PolyGamma[1, a]*
      PolyGamma[5, a]*PolyGamma[7, a] + 523783260*PolyGamma[0, a]*
      PolyGamma[1, a]^2*PolyGamma[5, a]*PolyGamma[7, a] +
     349188840*PolyGamma[0, a]^2*PolyGamma[2, a]*PolyGamma[5, a]*
      PolyGamma[7, a] + 349188840*PolyGamma[1, a]*PolyGamma[2, a]*
      PolyGamma[5, a]*PolyGamma[7, a] + 174594420*PolyGamma[0, a]*
      PolyGamma[3, a]*PolyGamma[5, a]*PolyGamma[7, a] +
     34918884*PolyGamma[4, a]*PolyGamma[5, a]*PolyGamma[7, a] +
     24942060*PolyGamma[0, a]^4*PolyGamma[6, a]*PolyGamma[7, a] +
     149652360*PolyGamma[0, a]^2*PolyGamma[1, a]*PolyGamma[6, a]*
      PolyGamma[7, a] + 74826180*PolyGamma[1, a]^2*PolyGamma[6, a]*
      PolyGamma[7, a] + 99768240*PolyGamma[0, a]*PolyGamma[2, a]*
      PolyGamma[6, a]*PolyGamma[7, a] + 24942060*PolyGamma[3, a]*
      PolyGamma[6, a]*PolyGamma[7, a] + 6235515*PolyGamma[0, a]^3*
      PolyGamma[7, a]^2 + 18706545*PolyGamma[0, a]*PolyGamma[1, a]*
      PolyGamma[7, a]^2 + 6235515*PolyGamma[2, a]*PolyGamma[7, a]^2 +
     92378*PolyGamma[0, a]^10*PolyGamma[8, a] + 4157010*PolyGamma[0, a]^8*
      PolyGamma[1, a]*PolyGamma[8, a] + 58198140*PolyGamma[0, a]^6*
      PolyGamma[1, a]^2*PolyGamma[8, a] + 290990700*PolyGamma[0, a]^4*
      PolyGamma[1, a]^3*PolyGamma[8, a] + 436486050*PolyGamma[0, a]^2*
      PolyGamma[1, a]^4*PolyGamma[8, a] + 87297210*PolyGamma[1, a]^5*
      PolyGamma[8, a] + 11085360*PolyGamma[0, a]^7*PolyGamma[2, a]*
      PolyGamma[8, a] + 232792560*PolyGamma[0, a]^5*PolyGamma[1, a]*
      PolyGamma[2, a]*PolyGamma[8, a] + 1163962800*PolyGamma[0, a]^3*
      PolyGamma[1, a]^2*PolyGamma[2, a]*PolyGamma[8, a] +
     1163962800*PolyGamma[0, a]*PolyGamma[1, a]^3*PolyGamma[2, a]*
      PolyGamma[8, a] + 193993800*PolyGamma[0, a]^4*PolyGamma[2, a]^2*
      PolyGamma[8, a] + 1163962800*PolyGamma[0, a]^2*PolyGamma[1, a]*
      PolyGamma[2, a]^2*PolyGamma[8, a] + 581981400*PolyGamma[1, a]^2*
      PolyGamma[2, a]^2*PolyGamma[8, a] + 258658400*PolyGamma[0, a]*
      PolyGamma[2, a]^3*PolyGamma[8, a] + 19399380*PolyGamma[0, a]^6*
      PolyGamma[3, a]*PolyGamma[8, a] + 290990700*PolyGamma[0, a]^4*
      PolyGamma[1, a]*PolyGamma[3, a]*PolyGamma[8, a] +
     872972100*PolyGamma[0, a]^2*PolyGamma[1, a]^2*PolyGamma[3, a]*
      PolyGamma[8, a] + 290990700*PolyGamma[1, a]^3*PolyGamma[3, a]*
      PolyGamma[8, a] + 387987600*PolyGamma[0, a]^3*PolyGamma[2, a]*
      PolyGamma[3, a]*PolyGamma[8, a] + 1163962800*PolyGamma[0, a]*
      PolyGamma[1, a]*PolyGamma[2, a]*PolyGamma[3, a]*PolyGamma[8, a] +
     193993800*PolyGamma[2, a]^2*PolyGamma[3, a]*PolyGamma[8, a] +
     145495350*PolyGamma[0, a]^2*PolyGamma[3, a]^2*PolyGamma[8, a] +
     145495350*PolyGamma[1, a]*PolyGamma[3, a]^2*PolyGamma[8, a] +
     23279256*PolyGamma[0, a]^5*PolyGamma[4, a]*PolyGamma[8, a] +
     232792560*PolyGamma[0, a]^3*PolyGamma[1, a]*PolyGamma[4, a]*
      PolyGamma[8, a] + 349188840*PolyGamma[0, a]*PolyGamma[1, a]^2*
      PolyGamma[4, a]*PolyGamma[8, a] + 232792560*PolyGamma[0, a]^2*
      PolyGamma[2, a]*PolyGamma[4, a]*PolyGamma[8, a] +
     232792560*PolyGamma[1, a]*PolyGamma[2, a]*PolyGamma[4, a]*
      PolyGamma[8, a] + 116396280*PolyGamma[0, a]*PolyGamma[3, a]*
      PolyGamma[4, a]*PolyGamma[8, a] + 11639628*PolyGamma[4, a]^2*
      PolyGamma[8, a] + 19399380*PolyGamma[0, a]^4*PolyGamma[5, a]*
      PolyGamma[8, a] + 116396280*PolyGamma[0, a]^2*PolyGamma[1, a]*
      PolyGamma[5, a]*PolyGamma[8, a] + 58198140*PolyGamma[1, a]^2*
      PolyGamma[5, a]*PolyGamma[8, a] + 77597520*PolyGamma[0, a]*
      PolyGamma[2, a]*PolyGamma[5, a]*PolyGamma[8, a] +
     19399380*PolyGamma[3, a]*PolyGamma[5, a]*PolyGamma[8, a] +
     11085360*PolyGamma[0, a]^3*PolyGamma[6, a]*PolyGamma[8, a] +
     33256080*PolyGamma[0, a]*PolyGamma[1, a]*PolyGamma[6, a]*
      PolyGamma[8, a] + 11085360*PolyGamma[2, a]*PolyGamma[6, a]*
      PolyGamma[8, a] + 4157010*PolyGamma[0, a]^2*PolyGamma[7, a]*
      PolyGamma[8, a] + 4157010*PolyGamma[1, a]*PolyGamma[7, a]*
      PolyGamma[8, a] + 461890*PolyGamma[0, a]*PolyGamma[8, a]^2 +
     92378*PolyGamma[0, a]^9*PolyGamma[9, a] + 3325608*PolyGamma[0, a]^7*
      PolyGamma[1, a]*PolyGamma[9, a] + 34918884*PolyGamma[0, a]^5*
      PolyGamma[1, a]^2*PolyGamma[9, a] + 116396280*PolyGamma[0, a]^3*
      PolyGamma[1, a]^3*PolyGamma[9, a] + 87297210*PolyGamma[0, a]*
      PolyGamma[1, a]^4*PolyGamma[9, a] + 7759752*PolyGamma[0, a]^6*
      PolyGamma[2, a]*PolyGamma[9, a] + 116396280*PolyGamma[0, a]^4*
      PolyGamma[1, a]*PolyGamma[2, a]*PolyGamma[9, a] +
     349188840*PolyGamma[0, a]^2*PolyGamma[1, a]^2*PolyGamma[2, a]*
      PolyGamma[9, a] + 116396280*PolyGamma[1, a]^3*PolyGamma[2, a]*
      PolyGamma[9, a] + 77597520*PolyGamma[0, a]^3*PolyGamma[2, a]^2*
      PolyGamma[9, a] + 232792560*PolyGamma[0, a]*PolyGamma[1, a]*
      PolyGamma[2, a]^2*PolyGamma[9, a] + 25865840*PolyGamma[2, a]^3*
      PolyGamma[9, a] + 11639628*PolyGamma[0, a]^5*PolyGamma[3, a]*
      PolyGamma[9, a] + 116396280*PolyGamma[0, a]^3*PolyGamma[1, a]*
      PolyGamma[3, a]*PolyGamma[9, a] + 174594420*PolyGamma[0, a]*
      PolyGamma[1, a]^2*PolyGamma[3, a]*PolyGamma[9, a] +
     116396280*PolyGamma[0, a]^2*PolyGamma[2, a]*PolyGamma[3, a]*
      PolyGamma[9, a] + 116396280*PolyGamma[1, a]*PolyGamma[2, a]*
      PolyGamma[3, a]*PolyGamma[9, a] + 29099070*PolyGamma[0, a]*
      PolyGamma[3, a]^2*PolyGamma[9, a] + 11639628*PolyGamma[0, a]^4*
      PolyGamma[4, a]*PolyGamma[9, a] + 69837768*PolyGamma[0, a]^2*
      PolyGamma[1, a]*PolyGamma[4, a]*PolyGamma[9, a] +
     34918884*PolyGamma[1, a]^2*PolyGamma[4, a]*PolyGamma[9, a] +
     46558512*PolyGamma[0, a]*PolyGamma[2, a]*PolyGamma[4, a]*
      PolyGamma[9, a] + 11639628*PolyGamma[3, a]*PolyGamma[4, a]*
      PolyGamma[9, a] + 7759752*PolyGamma[0, a]^3*PolyGamma[5, a]*
      PolyGamma[9, a] + 23279256*PolyGamma[0, a]*PolyGamma[1, a]*
      PolyGamma[5, a]*PolyGamma[9, a] + 7759752*PolyGamma[2, a]*
      PolyGamma[5, a]*PolyGamma[9, a] + 3325608*PolyGamma[0, a]^2*
      PolyGamma[6, a]*PolyGamma[9, a] + 3325608*PolyGamma[1, a]*
      PolyGamma[6, a]*PolyGamma[9, a] + 831402*PolyGamma[0, a]*
      PolyGamma[7, a]*PolyGamma[9, a] + 92378*PolyGamma[8, a]*
      PolyGamma[9, a] + 75582*PolyGamma[0, a]^8*PolyGamma[10, a] +
     2116296*PolyGamma[0, a]^6*PolyGamma[1, a]*PolyGamma[10, a] +
     15872220*PolyGamma[0, a]^4*PolyGamma[1, a]^2*PolyGamma[10, a] +
     31744440*PolyGamma[0, a]^2*PolyGamma[1, a]^3*PolyGamma[10, a] +
     7936110*PolyGamma[1, a]^4*PolyGamma[10, a] + 4232592*PolyGamma[0, a]^5*
      PolyGamma[2, a]*PolyGamma[10, a] + 42325920*PolyGamma[0, a]^3*
      PolyGamma[1, a]*PolyGamma[2, a]*PolyGamma[10, a] +
     63488880*PolyGamma[0, a]*PolyGamma[1, a]^2*PolyGamma[2, a]*
      PolyGamma[10, a] + 21162960*PolyGamma[0, a]^2*PolyGamma[2, a]^2*
      PolyGamma[10, a] + 21162960*PolyGamma[1, a]*PolyGamma[2, a]^2*
      PolyGamma[10, a] + 5290740*PolyGamma[0, a]^4*PolyGamma[3, a]*
      PolyGamma[10, a] + 31744440*PolyGamma[0, a]^2*PolyGamma[1, a]*
      PolyGamma[3, a]*PolyGamma[10, a] + 15872220*PolyGamma[1, a]^2*
      PolyGamma[3, a]*PolyGamma[10, a] + 21162960*PolyGamma[0, a]*
      PolyGamma[2, a]*PolyGamma[3, a]*PolyGamma[10, a] +
     2645370*PolyGamma[3, a]^2*PolyGamma[10, a] + 4232592*PolyGamma[0, a]^3*
      PolyGamma[4, a]*PolyGamma[10, a] + 12697776*PolyGamma[0, a]*
      PolyGamma[1, a]*PolyGamma[4, a]*PolyGamma[10, a] +
     4232592*PolyGamma[2, a]*PolyGamma[4, a]*PolyGamma[10, a] +
     2116296*PolyGamma[0, a]^2*PolyGamma[5, a]*PolyGamma[10, a] +
     2116296*PolyGamma[1, a]*PolyGamma[5, a]*PolyGamma[10, a] +
     604656*PolyGamma[0, a]*PolyGamma[6, a]*PolyGamma[10, a] +
     75582*PolyGamma[7, a]*PolyGamma[10, a] + 50388*PolyGamma[0, a]^7*
      PolyGamma[11, a] + 1058148*PolyGamma[0, a]^5*PolyGamma[1, a]*
      PolyGamma[11, a] + 5290740*PolyGamma[0, a]^3*PolyGamma[1, a]^2*
      PolyGamma[11, a] + 5290740*PolyGamma[0, a]*PolyGamma[1, a]^3*
      PolyGamma[11, a] + 1763580*PolyGamma[0, a]^4*PolyGamma[2, a]*
      PolyGamma[11, a] + 10581480*PolyGamma[0, a]^2*PolyGamma[1, a]*
      PolyGamma[2, a]*PolyGamma[11, a] + 5290740*PolyGamma[1, a]^2*
      PolyGamma[2, a]*PolyGamma[11, a] + 3527160*PolyGamma[0, a]*
      PolyGamma[2, a]^2*PolyGamma[11, a] + 1763580*PolyGamma[0, a]^3*
      PolyGamma[3, a]*PolyGamma[11, a] + 5290740*PolyGamma[0, a]*
      PolyGamma[1, a]*PolyGamma[3, a]*PolyGamma[11, a] +
     1763580*PolyGamma[2, a]*PolyGamma[3, a]*PolyGamma[11, a] +
     1058148*PolyGamma[0, a]^2*PolyGamma[4, a]*PolyGamma[11, a] +
     1058148*PolyGamma[1, a]*PolyGamma[4, a]*PolyGamma[11, a] +
     352716*PolyGamma[0, a]*PolyGamma[5, a]*PolyGamma[11, a] +
     50388*PolyGamma[6, a]*PolyGamma[11, a] + 27132*PolyGamma[0, a]^6*
      PolyGamma[12, a] + 406980*PolyGamma[0, a]^4*PolyGamma[1, a]*
      PolyGamma[12, a] + 1220940*PolyGamma[0, a]^2*PolyGamma[1, a]^2*
      PolyGamma[12, a] + 406980*PolyGamma[1, a]^3*PolyGamma[12, a] +
     542640*PolyGamma[0, a]^3*PolyGamma[2, a]*PolyGamma[12, a] +
     1627920*PolyGamma[0, a]*PolyGamma[1, a]*PolyGamma[2, a]*
      PolyGamma[12, a] + 271320*PolyGamma[2, a]^2*PolyGamma[12, a] +
     406980*PolyGamma[0, a]^2*PolyGamma[3, a]*PolyGamma[12, a] +
     406980*PolyGamma[1, a]*PolyGamma[3, a]*PolyGamma[12, a] +
     162792*PolyGamma[0, a]*PolyGamma[4, a]*PolyGamma[12, a] +
     27132*PolyGamma[5, a]*PolyGamma[12, a] + 11628*PolyGamma[0, a]^5*
      PolyGamma[13, a] + 116280*PolyGamma[0, a]^3*PolyGamma[1, a]*
      PolyGamma[13, a] + 174420*PolyGamma[0, a]*PolyGamma[1, a]^2*
      PolyGamma[13, a] + 116280*PolyGamma[0, a]^2*PolyGamma[2, a]*
      PolyGamma[13, a] + 116280*PolyGamma[1, a]*PolyGamma[2, a]*
      PolyGamma[13, a] + 58140*PolyGamma[0, a]*PolyGamma[3, a]*
      PolyGamma[13, a] + 11628*PolyGamma[4, a]*PolyGamma[13, a] +
     3876*PolyGamma[0, a]^4*PolyGamma[14, a] + 23256*PolyGamma[0, a]^2*
      PolyGamma[1, a]*PolyGamma[14, a] + 11628*PolyGamma[1, a]^2*
      PolyGamma[14, a] + 15504*PolyGamma[0, a]*PolyGamma[2, a]*
      PolyGamma[14, a] + 3876*PolyGamma[3, a]*PolyGamma[14, a] +
     969*PolyGamma[0, a]^3*PolyGamma[15, a] + 2907*PolyGamma[0, a]*
      PolyGamma[1, a]*PolyGamma[15, a] + 969*PolyGamma[2, a]*
      PolyGamma[15, a] + 171*PolyGamma[0, a]^2*PolyGamma[16, a] +
     171*PolyGamma[1, a]*PolyGamma[16, a] + 19*PolyGamma[0, a]*
      PolyGamma[17, a] + PolyGamma[18, a]

MBexpGam[a_, 20] = PolyGamma[0, a]^20 + 190*PolyGamma[0, a]^18*
      PolyGamma[1, a] + 14535*PolyGamma[0, a]^16*PolyGamma[1, a]^2 +
     581400*PolyGamma[0, a]^14*PolyGamma[1, a]^3 +
     13226850*PolyGamma[0, a]^12*PolyGamma[1, a]^4 +
     174594420*PolyGamma[0, a]^10*PolyGamma[1, a]^5 +
     1309458150*PolyGamma[0, a]^8*PolyGamma[1, a]^6 +
     5237832600*PolyGamma[0, a]^6*PolyGamma[1, a]^7 +
     9820936125*PolyGamma[0, a]^4*PolyGamma[1, a]^8 +
     6547290750*PolyGamma[0, a]^2*PolyGamma[1, a]^9 +
     654729075*PolyGamma[1, a]^10 + 1140*PolyGamma[0, a]^17*PolyGamma[2, a] +
     155040*PolyGamma[0, a]^15*PolyGamma[1, a]*PolyGamma[2, a] +
     8139600*PolyGamma[0, a]^13*PolyGamma[1, a]^2*PolyGamma[2, a] +
     211629600*PolyGamma[0, a]^11*PolyGamma[1, a]^3*PolyGamma[2, a] +
     2909907000*PolyGamma[0, a]^9*PolyGamma[1, a]^4*PolyGamma[2, a] +
     20951330400*PolyGamma[0, a]^7*PolyGamma[1, a]^5*PolyGamma[2, a] +
     73329656400*PolyGamma[0, a]^5*PolyGamma[1, a]^6*PolyGamma[2, a] +
     104756652000*PolyGamma[0, a]^3*PolyGamma[1, a]^7*PolyGamma[2, a] +
     39283744500*PolyGamma[0, a]*PolyGamma[1, a]^8*PolyGamma[2, a] +
     387600*PolyGamma[0, a]^14*PolyGamma[2, a]^2 +
     35271600*PolyGamma[0, a]^12*PolyGamma[1, a]*PolyGamma[2, a]^2 +
     1163962800*PolyGamma[0, a]^10*PolyGamma[1, a]^2*PolyGamma[2, a]^2 +
     17459442000*PolyGamma[0, a]^8*PolyGamma[1, a]^3*PolyGamma[2, a]^2 +
     122216094000*PolyGamma[0, a]^6*PolyGamma[1, a]^4*PolyGamma[2, a]^2 +
     366648282000*PolyGamma[0, a]^4*PolyGamma[1, a]^5*PolyGamma[2, a]^2 +
     366648282000*PolyGamma[0, a]^2*PolyGamma[1, a]^6*PolyGamma[2, a]^2 +
     52378326000*PolyGamma[1, a]^7*PolyGamma[2, a]^2 +
     47028800*PolyGamma[0, a]^11*PolyGamma[2, a]^3 +
     2586584000*PolyGamma[0, a]^9*PolyGamma[1, a]*PolyGamma[2, a]^3 +
     46558512000*PolyGamma[0, a]^7*PolyGamma[1, a]^2*PolyGamma[2, a]^3 +
     325909584000*PolyGamma[0, a]^5*PolyGamma[1, a]^3*PolyGamma[2, a]^3 +
     814773960000*PolyGamma[0, a]^3*PolyGamma[1, a]^4*PolyGamma[2, a]^3 +
     488864376000*PolyGamma[0, a]*PolyGamma[1, a]^5*PolyGamma[2, a]^3 +
     1939938000*PolyGamma[0, a]^8*PolyGamma[2, a]^4 +
     54318264000*PolyGamma[0, a]^6*PolyGamma[1, a]*PolyGamma[2, a]^4 +
     407386980000*PolyGamma[0, a]^4*PolyGamma[1, a]^2*PolyGamma[2, a]^4 +
     814773960000*PolyGamma[0, a]^2*PolyGamma[1, a]^3*PolyGamma[2, a]^4 +
     203693490000*PolyGamma[1, a]^4*PolyGamma[2, a]^4 +
     21727305600*PolyGamma[0, a]^5*PolyGamma[2, a]^5 +
     217273056000*PolyGamma[0, a]^3*PolyGamma[1, a]*PolyGamma[2, a]^5 +
     325909584000*PolyGamma[0, a]*PolyGamma[1, a]^2*PolyGamma[2, a]^5 +
     36212176000*PolyGamma[0, a]^2*PolyGamma[2, a]^6 +
     36212176000*PolyGamma[1, a]*PolyGamma[2, a]^6 +
     4845*PolyGamma[0, a]^16*PolyGamma[3, a] + 581400*PolyGamma[0, a]^14*
      PolyGamma[1, a]*PolyGamma[3, a] + 26453700*PolyGamma[0, a]^12*
      PolyGamma[1, a]^2*PolyGamma[3, a] + 581981400*PolyGamma[0, a]^10*
      PolyGamma[1, a]^3*PolyGamma[3, a] + 6547290750*PolyGamma[0, a]^8*
      PolyGamma[1, a]^4*PolyGamma[3, a] + 36664828200*PolyGamma[0, a]^6*
      PolyGamma[1, a]^5*PolyGamma[3, a] + 91662070500*PolyGamma[0, a]^4*
      PolyGamma[1, a]^6*PolyGamma[3, a] + 78567489000*PolyGamma[0, a]^2*
      PolyGamma[1, a]^7*PolyGamma[3, a] + 9820936125*PolyGamma[1, a]^8*
      PolyGamma[3, a] + 2713200*PolyGamma[0, a]^13*PolyGamma[2, a]*
      PolyGamma[3, a] + 211629600*PolyGamma[0, a]^11*PolyGamma[1, a]*
      PolyGamma[2, a]*PolyGamma[3, a] + 5819814000*PolyGamma[0, a]^9*
      PolyGamma[1, a]^2*PolyGamma[2, a]*PolyGamma[3, a] +
     69837768000*PolyGamma[0, a]^7*PolyGamma[1, a]^3*PolyGamma[2, a]*
      PolyGamma[3, a] + 366648282000*PolyGamma[0, a]^5*PolyGamma[1, a]^4*
      PolyGamma[2, a]*PolyGamma[3, a] + 733296564000*PolyGamma[0, a]^3*
      PolyGamma[1, a]^5*PolyGamma[2, a]*PolyGamma[3, a] +
     366648282000*PolyGamma[0, a]*PolyGamma[1, a]^6*PolyGamma[2, a]*
      PolyGamma[3, a] + 387987600*PolyGamma[0, a]^10*PolyGamma[2, a]^2*
      PolyGamma[3, a] + 17459442000*PolyGamma[0, a]^8*PolyGamma[1, a]*
      PolyGamma[2, a]^2*PolyGamma[3, a] + 244432188000*PolyGamma[0, a]^6*
      PolyGamma[1, a]^2*PolyGamma[2, a]^2*PolyGamma[3, a] +
     1222160940000*PolyGamma[0, a]^4*PolyGamma[1, a]^3*PolyGamma[2, a]^2*
      PolyGamma[3, a] + 1833241410000*PolyGamma[0, a]^2*PolyGamma[1, a]^4*
      PolyGamma[2, a]^2*PolyGamma[3, a] + 366648282000*PolyGamma[1, a]^5*
      PolyGamma[2, a]^2*PolyGamma[3, a] + 15519504000*PolyGamma[0, a]^7*
      PolyGamma[2, a]^3*PolyGamma[3, a] + 325909584000*PolyGamma[0, a]^5*
      PolyGamma[1, a]*PolyGamma[2, a]^3*PolyGamma[3, a] +
     1629547920000*PolyGamma[0, a]^3*PolyGamma[1, a]^2*PolyGamma[2, a]^3*
      PolyGamma[3, a] + 1629547920000*PolyGamma[0, a]*PolyGamma[1, a]^3*
      PolyGamma[2, a]^3*PolyGamma[3, a] + 135795660000*PolyGamma[0, a]^4*
      PolyGamma[2, a]^4*PolyGamma[3, a] + 814773960000*PolyGamma[0, a]^2*
      PolyGamma[1, a]*PolyGamma[2, a]^4*PolyGamma[3, a] +
     407386980000*PolyGamma[1, a]^2*PolyGamma[2, a]^4*PolyGamma[3, a] +
     108636528000*PolyGamma[0, a]*PolyGamma[2, a]^5*PolyGamma[3, a] +
     4408950*PolyGamma[0, a]^12*PolyGamma[3, a]^2 +
     290990700*PolyGamma[0, a]^10*PolyGamma[1, a]*PolyGamma[3, a]^2 +
     6547290750*PolyGamma[0, a]^8*PolyGamma[1, a]^2*PolyGamma[3, a]^2 +
     61108047000*PolyGamma[0, a]^6*PolyGamma[1, a]^3*PolyGamma[3, a]^2 +
     229155176250*PolyGamma[0, a]^4*PolyGamma[1, a]^4*PolyGamma[3, a]^2 +
     274986211500*PolyGamma[0, a]^2*PolyGamma[1, a]^5*PolyGamma[3, a]^2 +
     45831035250*PolyGamma[1, a]^6*PolyGamma[3, a]^2 +
     969969000*PolyGamma[0, a]^9*PolyGamma[2, a]*PolyGamma[3, a]^2 +
     34918884000*PolyGamma[0, a]^7*PolyGamma[1, a]*PolyGamma[2, a]*
      PolyGamma[3, a]^2 + 366648282000*PolyGamma[0, a]^5*PolyGamma[1, a]^2*
      PolyGamma[2, a]*PolyGamma[3, a]^2 + 1222160940000*PolyGamma[0, a]^3*
      PolyGamma[1, a]^3*PolyGamma[2, a]*PolyGamma[3, a]^2 +
     916620705000*PolyGamma[0, a]*PolyGamma[1, a]^4*PolyGamma[2, a]*
      PolyGamma[3, a]^2 + 40738698000*PolyGamma[0, a]^6*PolyGamma[2, a]^2*
      PolyGamma[3, a]^2 + 611080470000*PolyGamma[0, a]^4*PolyGamma[1, a]*
      PolyGamma[2, a]^2*PolyGamma[3, a]^2 + 1833241410000*PolyGamma[0, a]^2*
      PolyGamma[1, a]^2*PolyGamma[2, a]^2*PolyGamma[3, a]^2 +
     611080470000*PolyGamma[1, a]^3*PolyGamma[2, a]^2*PolyGamma[3, a]^2 +
     271591320000*PolyGamma[0, a]^3*PolyGamma[2, a]^3*PolyGamma[3, a]^2 +
     814773960000*PolyGamma[0, a]*PolyGamma[1, a]*PolyGamma[2, a]^3*
      PolyGamma[3, a]^2 + 67897830000*PolyGamma[2, a]^4*PolyGamma[3, a]^2 +
     727476750*PolyGamma[0, a]^8*PolyGamma[3, a]^3 +
     20369349000*PolyGamma[0, a]^6*PolyGamma[1, a]*PolyGamma[3, a]^3 +
     152770117500*PolyGamma[0, a]^4*PolyGamma[1, a]^2*PolyGamma[3, a]^3 +
     305540235000*PolyGamma[0, a]^2*PolyGamma[1, a]^3*PolyGamma[3, a]^3 +
     76385058750*PolyGamma[1, a]^4*PolyGamma[3, a]^3 +
     40738698000*PolyGamma[0, a]^5*PolyGamma[2, a]*PolyGamma[3, a]^3 +
     407386980000*PolyGamma[0, a]^3*PolyGamma[1, a]*PolyGamma[2, a]*
      PolyGamma[3, a]^3 + 611080470000*PolyGamma[0, a]*PolyGamma[1, a]^2*
      PolyGamma[2, a]*PolyGamma[3, a]^3 + 203693490000*PolyGamma[0, a]^2*
      PolyGamma[2, a]^2*PolyGamma[3, a]^3 + 203693490000*PolyGamma[1, a]*
      PolyGamma[2, a]^2*PolyGamma[3, a]^3 + 12730843125*PolyGamma[0, a]^4*
      PolyGamma[3, a]^4 + 76385058750*PolyGamma[0, a]^2*PolyGamma[1, a]*
      PolyGamma[3, a]^4 + 38192529375*PolyGamma[1, a]^2*PolyGamma[3, a]^4 +
     50923372500*PolyGamma[0, a]*PolyGamma[2, a]*PolyGamma[3, a]^4 +
     2546168625*PolyGamma[3, a]^5 + 15504*PolyGamma[0, a]^15*
      PolyGamma[4, a] + 1627920*PolyGamma[0, a]^13*PolyGamma[1, a]*
      PolyGamma[4, a] + 63488880*PolyGamma[0, a]^11*PolyGamma[1, a]^2*
      PolyGamma[4, a] + 1163962800*PolyGamma[0, a]^9*PolyGamma[1, a]^3*
      PolyGamma[4, a] + 10475665200*PolyGamma[0, a]^7*PolyGamma[1, a]^4*
      PolyGamma[4, a] + 43997793840*PolyGamma[0, a]^5*PolyGamma[1, a]^5*
      PolyGamma[4, a] + 73329656400*PolyGamma[0, a]^3*PolyGamma[1, a]^6*
      PolyGamma[4, a] + 31426995600*PolyGamma[0, a]*PolyGamma[1, a]^7*
      PolyGamma[4, a] + 7054320*PolyGamma[0, a]^12*PolyGamma[2, a]*
      PolyGamma[4, a] + 465585120*PolyGamma[0, a]^10*PolyGamma[1, a]*
      PolyGamma[2, a]*PolyGamma[4, a] + 10475665200*PolyGamma[0, a]^8*
      PolyGamma[1, a]^2*PolyGamma[2, a]*PolyGamma[4, a] +
     97772875200*PolyGamma[0, a]^6*PolyGamma[1, a]^3*PolyGamma[2, a]*
      PolyGamma[4, a] + 366648282000*PolyGamma[0, a]^4*PolyGamma[1, a]^4*
      PolyGamma[2, a]*PolyGamma[4, a] + 439977938400*PolyGamma[0, a]^2*
      PolyGamma[1, a]^5*PolyGamma[2, a]*PolyGamma[4, a] +
     73329656400*PolyGamma[1, a]^6*PolyGamma[2, a]*PolyGamma[4, a] +
     775975200*PolyGamma[0, a]^9*PolyGamma[2, a]^2*PolyGamma[4, a] +
     27935107200*PolyGamma[0, a]^7*PolyGamma[1, a]*PolyGamma[2, a]^2*
      PolyGamma[4, a] + 293318625600*PolyGamma[0, a]^5*PolyGamma[1, a]^2*
      PolyGamma[2, a]^2*PolyGamma[4, a] + 977728752000*PolyGamma[0, a]^3*
      PolyGamma[1, a]^3*PolyGamma[2, a]^2*PolyGamma[4, a] +
     733296564000*PolyGamma[0, a]*PolyGamma[1, a]^4*PolyGamma[2, a]^2*
      PolyGamma[4, a] + 21727305600*PolyGamma[0, a]^6*PolyGamma[2, a]^3*
      PolyGamma[4, a] + 325909584000*PolyGamma[0, a]^4*PolyGamma[1, a]*
      PolyGamma[2, a]^3*PolyGamma[4, a] + 977728752000*PolyGamma[0, a]^2*
      PolyGamma[1, a]^2*PolyGamma[2, a]^3*PolyGamma[4, a] +
     325909584000*PolyGamma[1, a]^3*PolyGamma[2, a]^3*PolyGamma[4, a] +
     108636528000*PolyGamma[0, a]^3*PolyGamma[2, a]^4*PolyGamma[4, a] +
     325909584000*PolyGamma[0, a]*PolyGamma[1, a]*PolyGamma[2, a]^4*
      PolyGamma[4, a] + 21727305600*PolyGamma[2, a]^5*PolyGamma[4, a] +
     21162960*PolyGamma[0, a]^11*PolyGamma[3, a]*PolyGamma[4, a] +
     1163962800*PolyGamma[0, a]^9*PolyGamma[1, a]*PolyGamma[3, a]*
      PolyGamma[4, a] + 20951330400*PolyGamma[0, a]^7*PolyGamma[1, a]^2*
      PolyGamma[3, a]*PolyGamma[4, a] + 146659312800*PolyGamma[0, a]^5*
      PolyGamma[1, a]^3*PolyGamma[3, a]*PolyGamma[4, a] +
     366648282000*PolyGamma[0, a]^3*PolyGamma[1, a]^4*PolyGamma[3, a]*
      PolyGamma[4, a] + 219988969200*PolyGamma[0, a]*PolyGamma[1, a]^5*
      PolyGamma[3, a]*PolyGamma[4, a] + 3491888400*PolyGamma[0, a]^8*
      PolyGamma[2, a]*PolyGamma[3, a]*PolyGamma[4, a] +
     97772875200*PolyGamma[0, a]^6*PolyGamma[1, a]*PolyGamma[2, a]*
      PolyGamma[3, a]*PolyGamma[4, a] + 733296564000*PolyGamma[0, a]^4*
      PolyGamma[1, a]^2*PolyGamma[2, a]*PolyGamma[3, a]*PolyGamma[4, a] +
     1466593128000*PolyGamma[0, a]^2*PolyGamma[1, a]^3*PolyGamma[2, a]*
      PolyGamma[3, a]*PolyGamma[4, a] + 366648282000*PolyGamma[1, a]^4*
      PolyGamma[2, a]*PolyGamma[3, a]*PolyGamma[4, a] +
     97772875200*PolyGamma[0, a]^5*PolyGamma[2, a]^2*PolyGamma[3, a]*
      PolyGamma[4, a] + 977728752000*PolyGamma[0, a]^3*PolyGamma[1, a]*
      PolyGamma[2, a]^2*PolyGamma[3, a]*PolyGamma[4, a] +
     1466593128000*PolyGamma[0, a]*PolyGamma[1, a]^2*PolyGamma[2, a]^2*
      PolyGamma[3, a]*PolyGamma[4, a] + 325909584000*PolyGamma[0, a]^2*
      PolyGamma[2, a]^3*PolyGamma[3, a]*PolyGamma[4, a] +
     325909584000*PolyGamma[1, a]*PolyGamma[2, a]^3*PolyGamma[3, a]*
      PolyGamma[4, a] + 3491888400*PolyGamma[0, a]^7*PolyGamma[3, a]^2*
      PolyGamma[4, a] + 73329656400*PolyGamma[0, a]^5*PolyGamma[1, a]*
      PolyGamma[3, a]^2*PolyGamma[4, a] + 366648282000*PolyGamma[0, a]^3*
      PolyGamma[1, a]^2*PolyGamma[3, a]^2*PolyGamma[4, a] +
     366648282000*PolyGamma[0, a]*PolyGamma[1, a]^3*PolyGamma[3, a]^2*
      PolyGamma[4, a] + 122216094000*PolyGamma[0, a]^4*PolyGamma[2, a]*
      PolyGamma[3, a]^2*PolyGamma[4, a] + 733296564000*PolyGamma[0, a]^2*
      PolyGamma[1, a]*PolyGamma[2, a]*PolyGamma[3, a]^2*PolyGamma[4, a] +
     366648282000*PolyGamma[1, a]^2*PolyGamma[2, a]*PolyGamma[3, a]^2*
      PolyGamma[4, a] + 244432188000*PolyGamma[0, a]*PolyGamma[2, a]^2*
      PolyGamma[3, a]^2*PolyGamma[4, a] + 40738698000*PolyGamma[0, a]^3*
      PolyGamma[3, a]^3*PolyGamma[4, a] + 122216094000*PolyGamma[0, a]*
      PolyGamma[1, a]*PolyGamma[3, a]^3*PolyGamma[4, a] +
     40738698000*PolyGamma[2, a]*PolyGamma[3, a]^3*PolyGamma[4, a] +
     23279256*PolyGamma[0, a]^10*PolyGamma[4, a]^2 +
     1047566520*PolyGamma[0, a]^8*PolyGamma[1, a]*PolyGamma[4, a]^2 +
     14665931280*PolyGamma[0, a]^6*PolyGamma[1, a]^2*PolyGamma[4, a]^2 +
     73329656400*PolyGamma[0, a]^4*PolyGamma[1, a]^3*PolyGamma[4, a]^2 +
     109994484600*PolyGamma[0, a]^2*PolyGamma[1, a]^4*PolyGamma[4, a]^2 +
     21998896920*PolyGamma[1, a]^5*PolyGamma[4, a]^2 +
     2793510720*PolyGamma[0, a]^7*PolyGamma[2, a]*PolyGamma[4, a]^2 +
     58663725120*PolyGamma[0, a]^5*PolyGamma[1, a]*PolyGamma[2, a]*
      PolyGamma[4, a]^2 + 293318625600*PolyGamma[0, a]^3*PolyGamma[1, a]^2*
      PolyGamma[2, a]*PolyGamma[4, a]^2 + 293318625600*PolyGamma[0, a]*
      PolyGamma[1, a]^3*PolyGamma[2, a]*PolyGamma[4, a]^2 +
     48886437600*PolyGamma[0, a]^4*PolyGamma[2, a]^2*PolyGamma[4, a]^2 +
     293318625600*PolyGamma[0, a]^2*PolyGamma[1, a]*PolyGamma[2, a]^2*
      PolyGamma[4, a]^2 + 146659312800*PolyGamma[1, a]^2*PolyGamma[2, a]^2*
      PolyGamma[4, a]^2 + 65181916800*PolyGamma[0, a]*PolyGamma[2, a]^3*
      PolyGamma[4, a]^2 + 4888643760*PolyGamma[0, a]^6*PolyGamma[3, a]*
      PolyGamma[4, a]^2 + 73329656400*PolyGamma[0, a]^4*PolyGamma[1, a]*
      PolyGamma[3, a]*PolyGamma[4, a]^2 + 219988969200*PolyGamma[0, a]^2*
      PolyGamma[1, a]^2*PolyGamma[3, a]*PolyGamma[4, a]^2 +
     73329656400*PolyGamma[1, a]^3*PolyGamma[3, a]*PolyGamma[4, a]^2 +
     97772875200*PolyGamma[0, a]^3*PolyGamma[2, a]*PolyGamma[3, a]*
      PolyGamma[4, a]^2 + 293318625600*PolyGamma[0, a]*PolyGamma[1, a]*
      PolyGamma[2, a]*PolyGamma[3, a]*PolyGamma[4, a]^2 +
     48886437600*PolyGamma[2, a]^2*PolyGamma[3, a]*PolyGamma[4, a]^2 +
     36664828200*PolyGamma[0, a]^2*PolyGamma[3, a]^2*PolyGamma[4, a]^2 +
     36664828200*PolyGamma[1, a]*PolyGamma[3, a]^2*PolyGamma[4, a]^2 +
     1955457504*PolyGamma[0, a]^5*PolyGamma[4, a]^3 +
     19554575040*PolyGamma[0, a]^3*PolyGamma[1, a]*PolyGamma[4, a]^3 +
     29331862560*PolyGamma[0, a]*PolyGamma[1, a]^2*PolyGamma[4, a]^3 +
     19554575040*PolyGamma[0, a]^2*PolyGamma[2, a]*PolyGamma[4, a]^3 +
     19554575040*PolyGamma[1, a]*PolyGamma[2, a]*PolyGamma[4, a]^3 +
     9777287520*PolyGamma[0, a]*PolyGamma[3, a]*PolyGamma[4, a]^3 +
     488864376*PolyGamma[4, a]^4 + 38760*PolyGamma[0, a]^14*PolyGamma[5, a] +
     3527160*PolyGamma[0, a]^12*PolyGamma[1, a]*PolyGamma[5, a] +
     116396280*PolyGamma[0, a]^10*PolyGamma[1, a]^2*PolyGamma[5, a] +
     1745944200*PolyGamma[0, a]^8*PolyGamma[1, a]^3*PolyGamma[5, a] +
     12221609400*PolyGamma[0, a]^6*PolyGamma[1, a]^4*PolyGamma[5, a] +
     36664828200*PolyGamma[0, a]^4*PolyGamma[1, a]^5*PolyGamma[5, a] +
     36664828200*PolyGamma[0, a]^2*PolyGamma[1, a]^6*PolyGamma[5, a] +
     5237832600*PolyGamma[1, a]^7*PolyGamma[5, a] +
     14108640*PolyGamma[0, a]^11*PolyGamma[2, a]*PolyGamma[5, a] +
     775975200*PolyGamma[0, a]^9*PolyGamma[1, a]*PolyGamma[2, a]*
      PolyGamma[5, a] + 13967553600*PolyGamma[0, a]^7*PolyGamma[1, a]^2*
      PolyGamma[2, a]*PolyGamma[5, a] + 97772875200*PolyGamma[0, a]^5*
      PolyGamma[1, a]^3*PolyGamma[2, a]*PolyGamma[5, a] +
     244432188000*PolyGamma[0, a]^3*PolyGamma[1, a]^4*PolyGamma[2, a]*
      PolyGamma[5, a] + 146659312800*PolyGamma[0, a]*PolyGamma[1, a]^5*
      PolyGamma[2, a]*PolyGamma[5, a] + 1163962800*PolyGamma[0, a]^8*
      PolyGamma[2, a]^2*PolyGamma[5, a] + 32590958400*PolyGamma[0, a]^6*
      PolyGamma[1, a]*PolyGamma[2, a]^2*PolyGamma[5, a] +
     244432188000*PolyGamma[0, a]^4*PolyGamma[1, a]^2*PolyGamma[2, a]^2*
      PolyGamma[5, a] + 488864376000*PolyGamma[0, a]^2*PolyGamma[1, a]^3*
      PolyGamma[2, a]^2*PolyGamma[5, a] + 122216094000*PolyGamma[1, a]^4*
      PolyGamma[2, a]^2*PolyGamma[5, a] + 21727305600*PolyGamma[0, a]^5*
      PolyGamma[2, a]^3*PolyGamma[5, a] + 217273056000*PolyGamma[0, a]^3*
      PolyGamma[1, a]*PolyGamma[2, a]^3*PolyGamma[5, a] +
     325909584000*PolyGamma[0, a]*PolyGamma[1, a]^2*PolyGamma[2, a]^3*
      PolyGamma[5, a] + 54318264000*PolyGamma[0, a]^2*PolyGamma[2, a]^4*
      PolyGamma[5, a] + 54318264000*PolyGamma[1, a]*PolyGamma[2, a]^4*
      PolyGamma[5, a] + 38798760*PolyGamma[0, a]^10*PolyGamma[3, a]*
      PolyGamma[5, a] + 1745944200*PolyGamma[0, a]^8*PolyGamma[1, a]*
      PolyGamma[3, a]*PolyGamma[5, a] + 24443218800*PolyGamma[0, a]^6*
      PolyGamma[1, a]^2*PolyGamma[3, a]*PolyGamma[5, a] +
     122216094000*PolyGamma[0, a]^4*PolyGamma[1, a]^3*PolyGamma[3, a]*
      PolyGamma[5, a] + 183324141000*PolyGamma[0, a]^2*PolyGamma[1, a]^4*
      PolyGamma[3, a]*PolyGamma[5, a] + 36664828200*PolyGamma[1, a]^5*
      PolyGamma[3, a]*PolyGamma[5, a] + 4655851200*PolyGamma[0, a]^7*
      PolyGamma[2, a]*PolyGamma[3, a]*PolyGamma[5, a] +
     97772875200*PolyGamma[0, a]^5*PolyGamma[1, a]*PolyGamma[2, a]*
      PolyGamma[3, a]*PolyGamma[5, a] + 488864376000*PolyGamma[0, a]^3*
      PolyGamma[1, a]^2*PolyGamma[2, a]*PolyGamma[3, a]*PolyGamma[5, a] +
     488864376000*PolyGamma[0, a]*PolyGamma[1, a]^3*PolyGamma[2, a]*
      PolyGamma[3, a]*PolyGamma[5, a] + 81477396000*PolyGamma[0, a]^4*
      PolyGamma[2, a]^2*PolyGamma[3, a]*PolyGamma[5, a] +
     488864376000*PolyGamma[0, a]^2*PolyGamma[1, a]*PolyGamma[2, a]^2*
      PolyGamma[3, a]*PolyGamma[5, a] + 244432188000*PolyGamma[1, a]^2*
      PolyGamma[2, a]^2*PolyGamma[3, a]*PolyGamma[5, a] +
     108636528000*PolyGamma[0, a]*PolyGamma[2, a]^3*PolyGamma[3, a]*
      PolyGamma[5, a] + 4073869800*PolyGamma[0, a]^6*PolyGamma[3, a]^2*
      PolyGamma[5, a] + 61108047000*PolyGamma[0, a]^4*PolyGamma[1, a]*
      PolyGamma[3, a]^2*PolyGamma[5, a] + 183324141000*PolyGamma[0, a]^2*
      PolyGamma[1, a]^2*PolyGamma[3, a]^2*PolyGamma[5, a] +
     61108047000*PolyGamma[1, a]^3*PolyGamma[3, a]^2*PolyGamma[5, a] +
     81477396000*PolyGamma[0, a]^3*PolyGamma[2, a]*PolyGamma[3, a]^2*
      PolyGamma[5, a] + 244432188000*PolyGamma[0, a]*PolyGamma[1, a]*
      PolyGamma[2, a]*PolyGamma[3, a]^2*PolyGamma[5, a] +
     40738698000*PolyGamma[2, a]^2*PolyGamma[3, a]^2*PolyGamma[5, a] +
     20369349000*PolyGamma[0, a]^2*PolyGamma[3, a]^3*PolyGamma[5, a] +
     20369349000*PolyGamma[1, a]*PolyGamma[3, a]^3*PolyGamma[5, a] +
     77597520*PolyGamma[0, a]^9*PolyGamma[4, a]*PolyGamma[5, a] +
     2793510720*PolyGamma[0, a]^7*PolyGamma[1, a]*PolyGamma[4, a]*
      PolyGamma[5, a] + 29331862560*PolyGamma[0, a]^5*PolyGamma[1, a]^2*
      PolyGamma[4, a]*PolyGamma[5, a] + 97772875200*PolyGamma[0, a]^3*
      PolyGamma[1, a]^3*PolyGamma[4, a]*PolyGamma[5, a] +
     73329656400*PolyGamma[0, a]*PolyGamma[1, a]^4*PolyGamma[4, a]*
      PolyGamma[5, a] + 6518191680*PolyGamma[0, a]^6*PolyGamma[2, a]*
      PolyGamma[4, a]*PolyGamma[5, a] + 97772875200*PolyGamma[0, a]^4*
      PolyGamma[1, a]*PolyGamma[2, a]*PolyGamma[4, a]*PolyGamma[5, a] +
     293318625600*PolyGamma[0, a]^2*PolyGamma[1, a]^2*PolyGamma[2, a]*
      PolyGamma[4, a]*PolyGamma[5, a] + 97772875200*PolyGamma[1, a]^3*
      PolyGamma[2, a]*PolyGamma[4, a]*PolyGamma[5, a] +
     65181916800*PolyGamma[0, a]^3*PolyGamma[2, a]^2*PolyGamma[4, a]*
      PolyGamma[5, a] + 195545750400*PolyGamma[0, a]*PolyGamma[1, a]*
      PolyGamma[2, a]^2*PolyGamma[4, a]*PolyGamma[5, a] +
     21727305600*PolyGamma[2, a]^3*PolyGamma[4, a]*PolyGamma[5, a] +
     9777287520*PolyGamma[0, a]^5*PolyGamma[3, a]*PolyGamma[4, a]*
      PolyGamma[5, a] + 97772875200*PolyGamma[0, a]^3*PolyGamma[1, a]*
      PolyGamma[3, a]*PolyGamma[4, a]*PolyGamma[5, a] +
     146659312800*PolyGamma[0, a]*PolyGamma[1, a]^2*PolyGamma[3, a]*
      PolyGamma[4, a]*PolyGamma[5, a] + 97772875200*PolyGamma[0, a]^2*
      PolyGamma[2, a]*PolyGamma[3, a]*PolyGamma[4, a]*PolyGamma[5, a] +
     97772875200*PolyGamma[1, a]*PolyGamma[2, a]*PolyGamma[3, a]*
      PolyGamma[4, a]*PolyGamma[5, a] + 24443218800*PolyGamma[0, a]*
      PolyGamma[3, a]^2*PolyGamma[4, a]*PolyGamma[5, a] +
     4888643760*PolyGamma[0, a]^4*PolyGamma[4, a]^2*PolyGamma[5, a] +
     29331862560*PolyGamma[0, a]^2*PolyGamma[1, a]*PolyGamma[4, a]^2*
      PolyGamma[5, a] + 14665931280*PolyGamma[1, a]^2*PolyGamma[4, a]^2*
      PolyGamma[5, a] + 19554575040*PolyGamma[0, a]*PolyGamma[2, a]*
      PolyGamma[4, a]^2*PolyGamma[5, a] + 4888643760*PolyGamma[3, a]*
      PolyGamma[4, a]^2*PolyGamma[5, a] + 58198140*PolyGamma[0, a]^8*
      PolyGamma[5, a]^2 + 1629547920*PolyGamma[0, a]^6*PolyGamma[1, a]*
      PolyGamma[5, a]^2 + 12221609400*PolyGamma[0, a]^4*PolyGamma[1, a]^2*
      PolyGamma[5, a]^2 + 24443218800*PolyGamma[0, a]^2*PolyGamma[1, a]^3*
      PolyGamma[5, a]^2 + 6110804700*PolyGamma[1, a]^4*PolyGamma[5, a]^2 +
     3259095840*PolyGamma[0, a]^5*PolyGamma[2, a]*PolyGamma[5, a]^2 +
     32590958400*PolyGamma[0, a]^3*PolyGamma[1, a]*PolyGamma[2, a]*
      PolyGamma[5, a]^2 + 48886437600*PolyGamma[0, a]*PolyGamma[1, a]^2*
      PolyGamma[2, a]*PolyGamma[5, a]^2 + 16295479200*PolyGamma[0, a]^2*
      PolyGamma[2, a]^2*PolyGamma[5, a]^2 + 16295479200*PolyGamma[1, a]*
      PolyGamma[2, a]^2*PolyGamma[5, a]^2 + 4073869800*PolyGamma[0, a]^4*
      PolyGamma[3, a]*PolyGamma[5, a]^2 + 24443218800*PolyGamma[0, a]^2*
      PolyGamma[1, a]*PolyGamma[3, a]*PolyGamma[5, a]^2 +
     12221609400*PolyGamma[1, a]^2*PolyGamma[3, a]*PolyGamma[5, a]^2 +
     16295479200*PolyGamma[0, a]*PolyGamma[2, a]*PolyGamma[3, a]*
      PolyGamma[5, a]^2 + 2036934900*PolyGamma[3, a]^2*PolyGamma[5, a]^2 +
     3259095840*PolyGamma[0, a]^3*PolyGamma[4, a]*PolyGamma[5, a]^2 +
     9777287520*PolyGamma[0, a]*PolyGamma[1, a]*PolyGamma[4, a]*
      PolyGamma[5, a]^2 + 3259095840*PolyGamma[2, a]*PolyGamma[4, a]*
      PolyGamma[5, a]^2 + 543182640*PolyGamma[0, a]^2*PolyGamma[5, a]^3 +
     543182640*PolyGamma[1, a]*PolyGamma[5, a]^3 + 77520*PolyGamma[0, a]^13*
      PolyGamma[6, a] + 6046560*PolyGamma[0, a]^11*PolyGamma[1, a]*
      PolyGamma[6, a] + 166280400*PolyGamma[0, a]^9*PolyGamma[1, a]^2*
      PolyGamma[6, a] + 1995364800*PolyGamma[0, a]^7*PolyGamma[1, a]^3*
      PolyGamma[6, a] + 10475665200*PolyGamma[0, a]^5*PolyGamma[1, a]^4*
      PolyGamma[6, a] + 20951330400*PolyGamma[0, a]^3*PolyGamma[1, a]^5*
      PolyGamma[6, a] + 10475665200*PolyGamma[0, a]*PolyGamma[1, a]^6*
      PolyGamma[6, a] + 22170720*PolyGamma[0, a]^10*PolyGamma[2, a]*
      PolyGamma[6, a] + 997682400*PolyGamma[0, a]^8*PolyGamma[1, a]*
      PolyGamma[2, a]*PolyGamma[6, a] + 13967553600*PolyGamma[0, a]^6*
      PolyGamma[1, a]^2*PolyGamma[2, a]*PolyGamma[6, a] +
     69837768000*PolyGamma[0, a]^4*PolyGamma[1, a]^3*PolyGamma[2, a]*
      PolyGamma[6, a] + 104756652000*PolyGamma[0, a]^2*PolyGamma[1, a]^4*
      PolyGamma[2, a]*PolyGamma[6, a] + 20951330400*PolyGamma[1, a]^5*
      PolyGamma[2, a]*PolyGamma[6, a] + 1330243200*PolyGamma[0, a]^7*
      PolyGamma[2, a]^2*PolyGamma[6, a] + 27935107200*PolyGamma[0, a]^5*
      PolyGamma[1, a]*PolyGamma[2, a]^2*PolyGamma[6, a] +
     139675536000*PolyGamma[0, a]^3*PolyGamma[1, a]^2*PolyGamma[2, a]^2*
      PolyGamma[6, a] + 139675536000*PolyGamma[0, a]*PolyGamma[1, a]^3*
      PolyGamma[2, a]^2*PolyGamma[6, a] + 15519504000*PolyGamma[0, a]^4*
      PolyGamma[2, a]^3*PolyGamma[6, a] + 93117024000*PolyGamma[0, a]^2*
      PolyGamma[1, a]*PolyGamma[2, a]^3*PolyGamma[6, a] +
     46558512000*PolyGamma[1, a]^2*PolyGamma[2, a]^3*PolyGamma[6, a] +
     15519504000*PolyGamma[0, a]*PolyGamma[2, a]^4*PolyGamma[6, a] +
     55426800*PolyGamma[0, a]^9*PolyGamma[3, a]*PolyGamma[6, a] +
     1995364800*PolyGamma[0, a]^7*PolyGamma[1, a]*PolyGamma[3, a]*
      PolyGamma[6, a] + 20951330400*PolyGamma[0, a]^5*PolyGamma[1, a]^2*
      PolyGamma[3, a]*PolyGamma[6, a] + 69837768000*PolyGamma[0, a]^3*
      PolyGamma[1, a]^3*PolyGamma[3, a]*PolyGamma[6, a] +
     52378326000*PolyGamma[0, a]*PolyGamma[1, a]^4*PolyGamma[3, a]*
      PolyGamma[6, a] + 4655851200*PolyGamma[0, a]^6*PolyGamma[2, a]*
      PolyGamma[3, a]*PolyGamma[6, a] + 69837768000*PolyGamma[0, a]^4*
      PolyGamma[1, a]*PolyGamma[2, a]*PolyGamma[3, a]*PolyGamma[6, a] +
     209513304000*PolyGamma[0, a]^2*PolyGamma[1, a]^2*PolyGamma[2, a]*
      PolyGamma[3, a]*PolyGamma[6, a] + 69837768000*PolyGamma[1, a]^3*
      PolyGamma[2, a]*PolyGamma[3, a]*PolyGamma[6, a] +
     46558512000*PolyGamma[0, a]^3*PolyGamma[2, a]^2*PolyGamma[3, a]*
      PolyGamma[6, a] + 139675536000*PolyGamma[0, a]*PolyGamma[1, a]*
      PolyGamma[2, a]^2*PolyGamma[3, a]*PolyGamma[6, a] +
     15519504000*PolyGamma[2, a]^3*PolyGamma[3, a]*PolyGamma[6, a] +
     3491888400*PolyGamma[0, a]^5*PolyGamma[3, a]^2*PolyGamma[6, a] +
     34918884000*PolyGamma[0, a]^3*PolyGamma[1, a]*PolyGamma[3, a]^2*
      PolyGamma[6, a] + 52378326000*PolyGamma[0, a]*PolyGamma[1, a]^2*
      PolyGamma[3, a]^2*PolyGamma[6, a] + 34918884000*PolyGamma[0, a]^2*
      PolyGamma[2, a]*PolyGamma[3, a]^2*PolyGamma[6, a] +
     34918884000*PolyGamma[1, a]*PolyGamma[2, a]*PolyGamma[3, a]^2*
      PolyGamma[6, a] + 5819814000*PolyGamma[0, a]*PolyGamma[3, a]^3*
      PolyGamma[6, a] + 99768240*PolyGamma[0, a]^8*PolyGamma[4, a]*
      PolyGamma[6, a] + 2793510720*PolyGamma[0, a]^6*PolyGamma[1, a]*
      PolyGamma[4, a]*PolyGamma[6, a] + 20951330400*PolyGamma[0, a]^4*
      PolyGamma[1, a]^2*PolyGamma[4, a]*PolyGamma[6, a] +
     41902660800*PolyGamma[0, a]^2*PolyGamma[1, a]^3*PolyGamma[4, a]*
      PolyGamma[6, a] + 10475665200*PolyGamma[1, a]^4*PolyGamma[4, a]*
      PolyGamma[6, a] + 5587021440*PolyGamma[0, a]^5*PolyGamma[2, a]*
      PolyGamma[4, a]*PolyGamma[6, a] + 55870214400*PolyGamma[0, a]^3*
      PolyGamma[1, a]*PolyGamma[2, a]*PolyGamma[4, a]*PolyGamma[6, a] +
     83805321600*PolyGamma[0, a]*PolyGamma[1, a]^2*PolyGamma[2, a]*
      PolyGamma[4, a]*PolyGamma[6, a] + 27935107200*PolyGamma[0, a]^2*
      PolyGamma[2, a]^2*PolyGamma[4, a]*PolyGamma[6, a] +
     27935107200*PolyGamma[1, a]*PolyGamma[2, a]^2*PolyGamma[4, a]*
      PolyGamma[6, a] + 6983776800*PolyGamma[0, a]^4*PolyGamma[3, a]*
      PolyGamma[4, a]*PolyGamma[6, a] + 41902660800*PolyGamma[0, a]^2*
      PolyGamma[1, a]*PolyGamma[3, a]*PolyGamma[4, a]*PolyGamma[6, a] +
     20951330400*PolyGamma[1, a]^2*PolyGamma[3, a]*PolyGamma[4, a]*
      PolyGamma[6, a] + 27935107200*PolyGamma[0, a]*PolyGamma[2, a]*
      PolyGamma[3, a]*PolyGamma[4, a]*PolyGamma[6, a] +
     3491888400*PolyGamma[3, a]^2*PolyGamma[4, a]*PolyGamma[6, a] +
     2793510720*PolyGamma[0, a]^3*PolyGamma[4, a]^2*PolyGamma[6, a] +
     8380532160*PolyGamma[0, a]*PolyGamma[1, a]*PolyGamma[4, a]^2*
      PolyGamma[6, a] + 2793510720*PolyGamma[2, a]*PolyGamma[4, a]^2*
      PolyGamma[6, a] + 133024320*PolyGamma[0, a]^7*PolyGamma[5, a]*
      PolyGamma[6, a] + 2793510720*PolyGamma[0, a]^5*PolyGamma[1, a]*
      PolyGamma[5, a]*PolyGamma[6, a] + 13967553600*PolyGamma[0, a]^3*
      PolyGamma[1, a]^2*PolyGamma[5, a]*PolyGamma[6, a] +
     13967553600*PolyGamma[0, a]*PolyGamma[1, a]^3*PolyGamma[5, a]*
      PolyGamma[6, a] + 4655851200*PolyGamma[0, a]^4*PolyGamma[2, a]*
      PolyGamma[5, a]*PolyGamma[6, a] + 27935107200*PolyGamma[0, a]^2*
      PolyGamma[1, a]*PolyGamma[2, a]*PolyGamma[5, a]*PolyGamma[6, a] +
     13967553600*PolyGamma[1, a]^2*PolyGamma[2, a]*PolyGamma[5, a]*
      PolyGamma[6, a] + 9311702400*PolyGamma[0, a]*PolyGamma[2, a]^2*
      PolyGamma[5, a]*PolyGamma[6, a] + 4655851200*PolyGamma[0, a]^3*
      PolyGamma[3, a]*PolyGamma[5, a]*PolyGamma[6, a] +
     13967553600*PolyGamma[0, a]*PolyGamma[1, a]*PolyGamma[3, a]*
      PolyGamma[5, a]*PolyGamma[6, a] + 4655851200*PolyGamma[2, a]*
      PolyGamma[3, a]*PolyGamma[5, a]*PolyGamma[6, a] +
     2793510720*PolyGamma[0, a]^2*PolyGamma[4, a]*PolyGamma[5, a]*
      PolyGamma[6, a] + 2793510720*PolyGamma[1, a]*PolyGamma[4, a]*
      PolyGamma[5, a]*PolyGamma[6, a] + 465585120*PolyGamma[0, a]*
      PolyGamma[5, a]^2*PolyGamma[6, a] + 66512160*PolyGamma[0, a]^6*
      PolyGamma[6, a]^2 + 997682400*PolyGamma[0, a]^4*PolyGamma[1, a]*
      PolyGamma[6, a]^2 + 2993047200*PolyGamma[0, a]^2*PolyGamma[1, a]^2*
      PolyGamma[6, a]^2 + 997682400*PolyGamma[1, a]^3*PolyGamma[6, a]^2 +
     1330243200*PolyGamma[0, a]^3*PolyGamma[2, a]*PolyGamma[6, a]^2 +
     3990729600*PolyGamma[0, a]*PolyGamma[1, a]*PolyGamma[2, a]*
      PolyGamma[6, a]^2 + 665121600*PolyGamma[2, a]^2*PolyGamma[6, a]^2 +
     997682400*PolyGamma[0, a]^2*PolyGamma[3, a]*PolyGamma[6, a]^2 +
     997682400*PolyGamma[1, a]*PolyGamma[3, a]*PolyGamma[6, a]^2 +
     399072960*PolyGamma[0, a]*PolyGamma[4, a]*PolyGamma[6, a]^2 +
     66512160*PolyGamma[5, a]*PolyGamma[6, a]^2 + 125970*PolyGamma[0, a]^12*
      PolyGamma[7, a] + 8314020*PolyGamma[0, a]^10*PolyGamma[1, a]*
      PolyGamma[7, a] + 187065450*PolyGamma[0, a]^8*PolyGamma[1, a]^2*
      PolyGamma[7, a] + 1745944200*PolyGamma[0, a]^6*PolyGamma[1, a]^3*
      PolyGamma[7, a] + 6547290750*PolyGamma[0, a]^4*PolyGamma[1, a]^4*
      PolyGamma[7, a] + 7856748900*PolyGamma[0, a]^2*PolyGamma[1, a]^5*
      PolyGamma[7, a] + 1309458150*PolyGamma[1, a]^6*PolyGamma[7, a] +
     27713400*PolyGamma[0, a]^9*PolyGamma[2, a]*PolyGamma[7, a] +
     997682400*PolyGamma[0, a]^7*PolyGamma[1, a]*PolyGamma[2, a]*
      PolyGamma[7, a] + 10475665200*PolyGamma[0, a]^5*PolyGamma[1, a]^2*
      PolyGamma[2, a]*PolyGamma[7, a] + 34918884000*PolyGamma[0, a]^3*
      PolyGamma[1, a]^3*PolyGamma[2, a]*PolyGamma[7, a] +
     26189163000*PolyGamma[0, a]*PolyGamma[1, a]^4*PolyGamma[2, a]*
      PolyGamma[7, a] + 1163962800*PolyGamma[0, a]^6*PolyGamma[2, a]^2*
      PolyGamma[7, a] + 17459442000*PolyGamma[0, a]^4*PolyGamma[1, a]*
      PolyGamma[2, a]^2*PolyGamma[7, a] + 52378326000*PolyGamma[0, a]^2*
      PolyGamma[1, a]^2*PolyGamma[2, a]^2*PolyGamma[7, a] +
     17459442000*PolyGamma[1, a]^3*PolyGamma[2, a]^2*PolyGamma[7, a] +
     7759752000*PolyGamma[0, a]^3*PolyGamma[2, a]^3*PolyGamma[7, a] +
     23279256000*PolyGamma[0, a]*PolyGamma[1, a]*PolyGamma[2, a]^3*
      PolyGamma[7, a] + 1939938000*PolyGamma[2, a]^4*PolyGamma[7, a] +
     62355150*PolyGamma[0, a]^8*PolyGamma[3, a]*PolyGamma[7, a] +
     1745944200*PolyGamma[0, a]^6*PolyGamma[1, a]*PolyGamma[3, a]*
      PolyGamma[7, a] + 13094581500*PolyGamma[0, a]^4*PolyGamma[1, a]^2*
      PolyGamma[3, a]*PolyGamma[7, a] + 26189163000*PolyGamma[0, a]^2*
      PolyGamma[1, a]^3*PolyGamma[3, a]*PolyGamma[7, a] +
     6547290750*PolyGamma[1, a]^4*PolyGamma[3, a]*PolyGamma[7, a] +
     3491888400*PolyGamma[0, a]^5*PolyGamma[2, a]*PolyGamma[3, a]*
      PolyGamma[7, a] + 34918884000*PolyGamma[0, a]^3*PolyGamma[1, a]*
      PolyGamma[2, a]*PolyGamma[3, a]*PolyGamma[7, a] +
     52378326000*PolyGamma[0, a]*PolyGamma[1, a]^2*PolyGamma[2, a]*
      PolyGamma[3, a]*PolyGamma[7, a] + 17459442000*PolyGamma[0, a]^2*
      PolyGamma[2, a]^2*PolyGamma[3, a]*PolyGamma[7, a] +
     17459442000*PolyGamma[1, a]*PolyGamma[2, a]^2*PolyGamma[3, a]*
      PolyGamma[7, a] + 2182430250*PolyGamma[0, a]^4*PolyGamma[3, a]^2*
      PolyGamma[7, a] + 13094581500*PolyGamma[0, a]^2*PolyGamma[1, a]*
      PolyGamma[3, a]^2*PolyGamma[7, a] + 6547290750*PolyGamma[1, a]^2*
      PolyGamma[3, a]^2*PolyGamma[7, a] + 8729721000*PolyGamma[0, a]*
      PolyGamma[2, a]*PolyGamma[3, a]^2*PolyGamma[7, a] +
     727476750*PolyGamma[3, a]^3*PolyGamma[7, a] + 99768240*PolyGamma[0, a]^7*
      PolyGamma[4, a]*PolyGamma[7, a] + 2095133040*PolyGamma[0, a]^5*
      PolyGamma[1, a]*PolyGamma[4, a]*PolyGamma[7, a] +
     10475665200*PolyGamma[0, a]^3*PolyGamma[1, a]^2*PolyGamma[4, a]*
      PolyGamma[7, a] + 10475665200*PolyGamma[0, a]*PolyGamma[1, a]^3*
      PolyGamma[4, a]*PolyGamma[7, a] + 3491888400*PolyGamma[0, a]^4*
      PolyGamma[2, a]*PolyGamma[4, a]*PolyGamma[7, a] +
     20951330400*PolyGamma[0, a]^2*PolyGamma[1, a]*PolyGamma[2, a]*
      PolyGamma[4, a]*PolyGamma[7, a] + 10475665200*PolyGamma[1, a]^2*
      PolyGamma[2, a]*PolyGamma[4, a]*PolyGamma[7, a] +
     6983776800*PolyGamma[0, a]*PolyGamma[2, a]^2*PolyGamma[4, a]*
      PolyGamma[7, a] + 3491888400*PolyGamma[0, a]^3*PolyGamma[3, a]*
      PolyGamma[4, a]*PolyGamma[7, a] + 10475665200*PolyGamma[0, a]*
      PolyGamma[1, a]*PolyGamma[3, a]*PolyGamma[4, a]*PolyGamma[7, a] +
     3491888400*PolyGamma[2, a]*PolyGamma[3, a]*PolyGamma[4, a]*
      PolyGamma[7, a] + 1047566520*PolyGamma[0, a]^2*PolyGamma[4, a]^2*
      PolyGamma[7, a] + 1047566520*PolyGamma[1, a]*PolyGamma[4, a]^2*
      PolyGamma[7, a] + 116396280*PolyGamma[0, a]^6*PolyGamma[5, a]*
      PolyGamma[7, a] + 1745944200*PolyGamma[0, a]^4*PolyGamma[1, a]*
      PolyGamma[5, a]*PolyGamma[7, a] + 5237832600*PolyGamma[0, a]^2*
      PolyGamma[1, a]^2*PolyGamma[5, a]*PolyGamma[7, a] +
     1745944200*PolyGamma[1, a]^3*PolyGamma[5, a]*PolyGamma[7, a] +
     2327925600*PolyGamma[0, a]^3*PolyGamma[2, a]*PolyGamma[5, a]*
      PolyGamma[7, a] + 6983776800*PolyGamma[0, a]*PolyGamma[1, a]*
      PolyGamma[2, a]*PolyGamma[5, a]*PolyGamma[7, a] +
     1163962800*PolyGamma[2, a]^2*PolyGamma[5, a]*PolyGamma[7, a] +
     1745944200*PolyGamma[0, a]^2*PolyGamma[3, a]*PolyGamma[5, a]*
      PolyGamma[7, a] + 1745944200*PolyGamma[1, a]*PolyGamma[3, a]*
      PolyGamma[5, a]*PolyGamma[7, a] + 698377680*PolyGamma[0, a]*
      PolyGamma[4, a]*PolyGamma[5, a]*PolyGamma[7, a] +
     58198140*PolyGamma[5, a]^2*PolyGamma[7, a] + 99768240*PolyGamma[0, a]^5*
      PolyGamma[6, a]*PolyGamma[7, a] + 997682400*PolyGamma[0, a]^3*
      PolyGamma[1, a]*PolyGamma[6, a]*PolyGamma[7, a] +
     1496523600*PolyGamma[0, a]*PolyGamma[1, a]^2*PolyGamma[6, a]*
      PolyGamma[7, a] + 997682400*PolyGamma[0, a]^2*PolyGamma[2, a]*
      PolyGamma[6, a]*PolyGamma[7, a] + 997682400*PolyGamma[1, a]*
      PolyGamma[2, a]*PolyGamma[6, a]*PolyGamma[7, a] +
     498841200*PolyGamma[0, a]*PolyGamma[3, a]*PolyGamma[6, a]*
      PolyGamma[7, a] + 99768240*PolyGamma[4, a]*PolyGamma[6, a]*
      PolyGamma[7, a] + 31177575*PolyGamma[0, a]^4*PolyGamma[7, a]^2 +
     187065450*PolyGamma[0, a]^2*PolyGamma[1, a]*PolyGamma[7, a]^2 +
     93532725*PolyGamma[1, a]^2*PolyGamma[7, a]^2 +
     124710300*PolyGamma[0, a]*PolyGamma[2, a]*PolyGamma[7, a]^2 +
     31177575*PolyGamma[3, a]*PolyGamma[7, a]^2 + 167960*PolyGamma[0, a]^11*
      PolyGamma[8, a] + 9237800*PolyGamma[0, a]^9*PolyGamma[1, a]*
      PolyGamma[8, a] + 166280400*PolyGamma[0, a]^7*PolyGamma[1, a]^2*
      PolyGamma[8, a] + 1163962800*PolyGamma[0, a]^5*PolyGamma[1, a]^3*
      PolyGamma[8, a] + 2909907000*PolyGamma[0, a]^3*PolyGamma[1, a]^4*
      PolyGamma[8, a] + 1745944200*PolyGamma[0, a]*PolyGamma[1, a]^5*
      PolyGamma[8, a] + 27713400*PolyGamma[0, a]^8*PolyGamma[2, a]*
      PolyGamma[8, a] + 775975200*PolyGamma[0, a]^6*PolyGamma[1, a]*
      PolyGamma[2, a]*PolyGamma[8, a] + 5819814000*PolyGamma[0, a]^4*
      PolyGamma[1, a]^2*PolyGamma[2, a]*PolyGamma[8, a] +
     11639628000*PolyGamma[0, a]^2*PolyGamma[1, a]^3*PolyGamma[2, a]*
      PolyGamma[8, a] + 2909907000*PolyGamma[1, a]^4*PolyGamma[2, a]*
      PolyGamma[8, a] + 775975200*PolyGamma[0, a]^5*PolyGamma[2, a]^2*
      PolyGamma[8, a] + 7759752000*PolyGamma[0, a]^3*PolyGamma[1, a]*
      PolyGamma[2, a]^2*PolyGamma[8, a] + 11639628000*PolyGamma[0, a]*
      PolyGamma[1, a]^2*PolyGamma[2, a]^2*PolyGamma[8, a] +
     2586584000*PolyGamma[0, a]^2*PolyGamma[2, a]^3*PolyGamma[8, a] +
     2586584000*PolyGamma[1, a]*PolyGamma[2, a]^3*PolyGamma[8, a] +
     55426800*PolyGamma[0, a]^7*PolyGamma[3, a]*PolyGamma[8, a] +
     1163962800*PolyGamma[0, a]^5*PolyGamma[1, a]*PolyGamma[3, a]*
      PolyGamma[8, a] + 5819814000*PolyGamma[0, a]^3*PolyGamma[1, a]^2*
      PolyGamma[3, a]*PolyGamma[8, a] + 5819814000*PolyGamma[0, a]*
      PolyGamma[1, a]^3*PolyGamma[3, a]*PolyGamma[8, a] +
     1939938000*PolyGamma[0, a]^4*PolyGamma[2, a]*PolyGamma[3, a]*
      PolyGamma[8, a] + 11639628000*PolyGamma[0, a]^2*PolyGamma[1, a]*
      PolyGamma[2, a]*PolyGamma[3, a]*PolyGamma[8, a] +
     5819814000*PolyGamma[1, a]^2*PolyGamma[2, a]*PolyGamma[3, a]*
      PolyGamma[8, a] + 3879876000*PolyGamma[0, a]*PolyGamma[2, a]^2*
      PolyGamma[3, a]*PolyGamma[8, a] + 969969000*PolyGamma[0, a]^3*
      PolyGamma[3, a]^2*PolyGamma[8, a] + 2909907000*PolyGamma[0, a]*
      PolyGamma[1, a]*PolyGamma[3, a]^2*PolyGamma[8, a] +
     969969000*PolyGamma[2, a]*PolyGamma[3, a]^2*PolyGamma[8, a] +
     77597520*PolyGamma[0, a]^6*PolyGamma[4, a]*PolyGamma[8, a] +
     1163962800*PolyGamma[0, a]^4*PolyGamma[1, a]*PolyGamma[4, a]*
      PolyGamma[8, a] + 3491888400*PolyGamma[0, a]^2*PolyGamma[1, a]^2*
      PolyGamma[4, a]*PolyGamma[8, a] + 1163962800*PolyGamma[1, a]^3*
      PolyGamma[4, a]*PolyGamma[8, a] + 1551950400*PolyGamma[0, a]^3*
      PolyGamma[2, a]*PolyGamma[4, a]*PolyGamma[8, a] +
     4655851200*PolyGamma[0, a]*PolyGamma[1, a]*PolyGamma[2, a]*
      PolyGamma[4, a]*PolyGamma[8, a] + 775975200*PolyGamma[2, a]^2*
      PolyGamma[4, a]*PolyGamma[8, a] + 1163962800*PolyGamma[0, a]^2*
      PolyGamma[3, a]*PolyGamma[4, a]*PolyGamma[8, a] +
     1163962800*PolyGamma[1, a]*PolyGamma[3, a]*PolyGamma[4, a]*
      PolyGamma[8, a] + 232792560*PolyGamma[0, a]*PolyGamma[4, a]^2*
      PolyGamma[8, a] + 77597520*PolyGamma[0, a]^5*PolyGamma[5, a]*
      PolyGamma[8, a] + 775975200*PolyGamma[0, a]^3*PolyGamma[1, a]*
      PolyGamma[5, a]*PolyGamma[8, a] + 1163962800*PolyGamma[0, a]*
      PolyGamma[1, a]^2*PolyGamma[5, a]*PolyGamma[8, a] +
     775975200*PolyGamma[0, a]^2*PolyGamma[2, a]*PolyGamma[5, a]*
      PolyGamma[8, a] + 775975200*PolyGamma[1, a]*PolyGamma[2, a]*
      PolyGamma[5, a]*PolyGamma[8, a] + 387987600*PolyGamma[0, a]*
      PolyGamma[3, a]*PolyGamma[5, a]*PolyGamma[8, a] +
     77597520*PolyGamma[4, a]*PolyGamma[5, a]*PolyGamma[8, a] +
     55426800*PolyGamma[0, a]^4*PolyGamma[6, a]*PolyGamma[8, a] +
     332560800*PolyGamma[0, a]^2*PolyGamma[1, a]*PolyGamma[6, a]*
      PolyGamma[8, a] + 166280400*PolyGamma[1, a]^2*PolyGamma[6, a]*
      PolyGamma[8, a] + 221707200*PolyGamma[0, a]*PolyGamma[2, a]*
      PolyGamma[6, a]*PolyGamma[8, a] + 55426800*PolyGamma[3, a]*
      PolyGamma[6, a]*PolyGamma[8, a] + 27713400*PolyGamma[0, a]^3*
      PolyGamma[7, a]*PolyGamma[8, a] + 83140200*PolyGamma[0, a]*
      PolyGamma[1, a]*PolyGamma[7, a]*PolyGamma[8, a] +
     27713400*PolyGamma[2, a]*PolyGamma[7, a]*PolyGamma[8, a] +
     4618900*PolyGamma[0, a]^2*PolyGamma[8, a]^2 + 4618900*PolyGamma[1, a]*
      PolyGamma[8, a]^2 + 184756*PolyGamma[0, a]^10*PolyGamma[9, a] +
     8314020*PolyGamma[0, a]^8*PolyGamma[1, a]*PolyGamma[9, a] +
     116396280*PolyGamma[0, a]^6*PolyGamma[1, a]^2*PolyGamma[9, a] +
     581981400*PolyGamma[0, a]^4*PolyGamma[1, a]^3*PolyGamma[9, a] +
     872972100*PolyGamma[0, a]^2*PolyGamma[1, a]^4*PolyGamma[9, a] +
     174594420*PolyGamma[1, a]^5*PolyGamma[9, a] + 22170720*PolyGamma[0, a]^7*
      PolyGamma[2, a]*PolyGamma[9, a] + 465585120*PolyGamma[0, a]^5*
      PolyGamma[1, a]*PolyGamma[2, a]*PolyGamma[9, a] +
     2327925600*PolyGamma[0, a]^3*PolyGamma[1, a]^2*PolyGamma[2, a]*
      PolyGamma[9, a] + 2327925600*PolyGamma[0, a]*PolyGamma[1, a]^3*
      PolyGamma[2, a]*PolyGamma[9, a] + 387987600*PolyGamma[0, a]^4*
      PolyGamma[2, a]^2*PolyGamma[9, a] + 2327925600*PolyGamma[0, a]^2*
      PolyGamma[1, a]*PolyGamma[2, a]^2*PolyGamma[9, a] +
     1163962800*PolyGamma[1, a]^2*PolyGamma[2, a]^2*PolyGamma[9, a] +
     517316800*PolyGamma[0, a]*PolyGamma[2, a]^3*PolyGamma[9, a] +
     38798760*PolyGamma[0, a]^6*PolyGamma[3, a]*PolyGamma[9, a] +
     581981400*PolyGamma[0, a]^4*PolyGamma[1, a]*PolyGamma[3, a]*
      PolyGamma[9, a] + 1745944200*PolyGamma[0, a]^2*PolyGamma[1, a]^2*
      PolyGamma[3, a]*PolyGamma[9, a] + 581981400*PolyGamma[1, a]^3*
      PolyGamma[3, a]*PolyGamma[9, a] + 775975200*PolyGamma[0, a]^3*
      PolyGamma[2, a]*PolyGamma[3, a]*PolyGamma[9, a] +
     2327925600*PolyGamma[0, a]*PolyGamma[1, a]*PolyGamma[2, a]*
      PolyGamma[3, a]*PolyGamma[9, a] + 387987600*PolyGamma[2, a]^2*
      PolyGamma[3, a]*PolyGamma[9, a] + 290990700*PolyGamma[0, a]^2*
      PolyGamma[3, a]^2*PolyGamma[9, a] + 290990700*PolyGamma[1, a]*
      PolyGamma[3, a]^2*PolyGamma[9, a] + 46558512*PolyGamma[0, a]^5*
      PolyGamma[4, a]*PolyGamma[9, a] + 465585120*PolyGamma[0, a]^3*
      PolyGamma[1, a]*PolyGamma[4, a]*PolyGamma[9, a] +
     698377680*PolyGamma[0, a]*PolyGamma[1, a]^2*PolyGamma[4, a]*
      PolyGamma[9, a] + 465585120*PolyGamma[0, a]^2*PolyGamma[2, a]*
      PolyGamma[4, a]*PolyGamma[9, a] + 465585120*PolyGamma[1, a]*
      PolyGamma[2, a]*PolyGamma[4, a]*PolyGamma[9, a] +
     232792560*PolyGamma[0, a]*PolyGamma[3, a]*PolyGamma[4, a]*
      PolyGamma[9, a] + 23279256*PolyGamma[4, a]^2*PolyGamma[9, a] +
     38798760*PolyGamma[0, a]^4*PolyGamma[5, a]*PolyGamma[9, a] +
     232792560*PolyGamma[0, a]^2*PolyGamma[1, a]*PolyGamma[5, a]*
      PolyGamma[9, a] + 116396280*PolyGamma[1, a]^2*PolyGamma[5, a]*
      PolyGamma[9, a] + 155195040*PolyGamma[0, a]*PolyGamma[2, a]*
      PolyGamma[5, a]*PolyGamma[9, a] + 38798760*PolyGamma[3, a]*
      PolyGamma[5, a]*PolyGamma[9, a] + 22170720*PolyGamma[0, a]^3*
      PolyGamma[6, a]*PolyGamma[9, a] + 66512160*PolyGamma[0, a]*
      PolyGamma[1, a]*PolyGamma[6, a]*PolyGamma[9, a] +
     22170720*PolyGamma[2, a]*PolyGamma[6, a]*PolyGamma[9, a] +
     8314020*PolyGamma[0, a]^2*PolyGamma[7, a]*PolyGamma[9, a] +
     8314020*PolyGamma[1, a]*PolyGamma[7, a]*PolyGamma[9, a] +
     1847560*PolyGamma[0, a]*PolyGamma[8, a]*PolyGamma[9, a] +
     92378*PolyGamma[9, a]^2 + 167960*PolyGamma[0, a]^9*PolyGamma[10, a] +
     6046560*PolyGamma[0, a]^7*PolyGamma[1, a]*PolyGamma[10, a] +
     63488880*PolyGamma[0, a]^5*PolyGamma[1, a]^2*PolyGamma[10, a] +
     211629600*PolyGamma[0, a]^3*PolyGamma[1, a]^3*PolyGamma[10, a] +
     158722200*PolyGamma[0, a]*PolyGamma[1, a]^4*PolyGamma[10, a] +
     14108640*PolyGamma[0, a]^6*PolyGamma[2, a]*PolyGamma[10, a] +
     211629600*PolyGamma[0, a]^4*PolyGamma[1, a]*PolyGamma[2, a]*
      PolyGamma[10, a] + 634888800*PolyGamma[0, a]^2*PolyGamma[1, a]^2*
      PolyGamma[2, a]*PolyGamma[10, a] + 211629600*PolyGamma[1, a]^3*
      PolyGamma[2, a]*PolyGamma[10, a] + 141086400*PolyGamma[0, a]^3*
      PolyGamma[2, a]^2*PolyGamma[10, a] + 423259200*PolyGamma[0, a]*
      PolyGamma[1, a]*PolyGamma[2, a]^2*PolyGamma[10, a] +
     47028800*PolyGamma[2, a]^3*PolyGamma[10, a] + 21162960*PolyGamma[0, a]^5*
      PolyGamma[3, a]*PolyGamma[10, a] + 211629600*PolyGamma[0, a]^3*
      PolyGamma[1, a]*PolyGamma[3, a]*PolyGamma[10, a] +
     317444400*PolyGamma[0, a]*PolyGamma[1, a]^2*PolyGamma[3, a]*
      PolyGamma[10, a] + 211629600*PolyGamma[0, a]^2*PolyGamma[2, a]*
      PolyGamma[3, a]*PolyGamma[10, a] + 211629600*PolyGamma[1, a]*
      PolyGamma[2, a]*PolyGamma[3, a]*PolyGamma[10, a] +
     52907400*PolyGamma[0, a]*PolyGamma[3, a]^2*PolyGamma[10, a] +
     21162960*PolyGamma[0, a]^4*PolyGamma[4, a]*PolyGamma[10, a] +
     126977760*PolyGamma[0, a]^2*PolyGamma[1, a]*PolyGamma[4, a]*
      PolyGamma[10, a] + 63488880*PolyGamma[1, a]^2*PolyGamma[4, a]*
      PolyGamma[10, a] + 84651840*PolyGamma[0, a]*PolyGamma[2, a]*
      PolyGamma[4, a]*PolyGamma[10, a] + 21162960*PolyGamma[3, a]*
      PolyGamma[4, a]*PolyGamma[10, a] + 14108640*PolyGamma[0, a]^3*
      PolyGamma[5, a]*PolyGamma[10, a] + 42325920*PolyGamma[0, a]*
      PolyGamma[1, a]*PolyGamma[5, a]*PolyGamma[10, a] +
     14108640*PolyGamma[2, a]*PolyGamma[5, a]*PolyGamma[10, a] +
     6046560*PolyGamma[0, a]^2*PolyGamma[6, a]*PolyGamma[10, a] +
     6046560*PolyGamma[1, a]*PolyGamma[6, a]*PolyGamma[10, a] +
     1511640*PolyGamma[0, a]*PolyGamma[7, a]*PolyGamma[10, a] +
     167960*PolyGamma[8, a]*PolyGamma[10, a] + 125970*PolyGamma[0, a]^8*
      PolyGamma[11, a] + 3527160*PolyGamma[0, a]^6*PolyGamma[1, a]*
      PolyGamma[11, a] + 26453700*PolyGamma[0, a]^4*PolyGamma[1, a]^2*
      PolyGamma[11, a] + 52907400*PolyGamma[0, a]^2*PolyGamma[1, a]^3*
      PolyGamma[11, a] + 13226850*PolyGamma[1, a]^4*PolyGamma[11, a] +
     7054320*PolyGamma[0, a]^5*PolyGamma[2, a]*PolyGamma[11, a] +
     70543200*PolyGamma[0, a]^3*PolyGamma[1, a]*PolyGamma[2, a]*
      PolyGamma[11, a] + 105814800*PolyGamma[0, a]*PolyGamma[1, a]^2*
      PolyGamma[2, a]*PolyGamma[11, a] + 35271600*PolyGamma[0, a]^2*
      PolyGamma[2, a]^2*PolyGamma[11, a] + 35271600*PolyGamma[1, a]*
      PolyGamma[2, a]^2*PolyGamma[11, a] + 8817900*PolyGamma[0, a]^4*
      PolyGamma[3, a]*PolyGamma[11, a] + 52907400*PolyGamma[0, a]^2*
      PolyGamma[1, a]*PolyGamma[3, a]*PolyGamma[11, a] +
     26453700*PolyGamma[1, a]^2*PolyGamma[3, a]*PolyGamma[11, a] +
     35271600*PolyGamma[0, a]*PolyGamma[2, a]*PolyGamma[3, a]*
      PolyGamma[11, a] + 4408950*PolyGamma[3, a]^2*PolyGamma[11, a] +
     7054320*PolyGamma[0, a]^3*PolyGamma[4, a]*PolyGamma[11, a] +
     21162960*PolyGamma[0, a]*PolyGamma[1, a]*PolyGamma[4, a]*
      PolyGamma[11, a] + 7054320*PolyGamma[2, a]*PolyGamma[4, a]*
      PolyGamma[11, a] + 3527160*PolyGamma[0, a]^2*PolyGamma[5, a]*
      PolyGamma[11, a] + 3527160*PolyGamma[1, a]*PolyGamma[5, a]*
      PolyGamma[11, a] + 1007760*PolyGamma[0, a]*PolyGamma[6, a]*
      PolyGamma[11, a] + 125970*PolyGamma[7, a]*PolyGamma[11, a] +
     77520*PolyGamma[0, a]^7*PolyGamma[12, a] + 1627920*PolyGamma[0, a]^5*
      PolyGamma[1, a]*PolyGamma[12, a] + 8139600*PolyGamma[0, a]^3*
      PolyGamma[1, a]^2*PolyGamma[12, a] + 8139600*PolyGamma[0, a]*
      PolyGamma[1, a]^3*PolyGamma[12, a] + 2713200*PolyGamma[0, a]^4*
      PolyGamma[2, a]*PolyGamma[12, a] + 16279200*PolyGamma[0, a]^2*
      PolyGamma[1, a]*PolyGamma[2, a]*PolyGamma[12, a] +
     8139600*PolyGamma[1, a]^2*PolyGamma[2, a]*PolyGamma[12, a] +
     5426400*PolyGamma[0, a]*PolyGamma[2, a]^2*PolyGamma[12, a] +
     2713200*PolyGamma[0, a]^3*PolyGamma[3, a]*PolyGamma[12, a] +
     8139600*PolyGamma[0, a]*PolyGamma[1, a]*PolyGamma[3, a]*
      PolyGamma[12, a] + 2713200*PolyGamma[2, a]*PolyGamma[3, a]*
      PolyGamma[12, a] + 1627920*PolyGamma[0, a]^2*PolyGamma[4, a]*
      PolyGamma[12, a] + 1627920*PolyGamma[1, a]*PolyGamma[4, a]*
      PolyGamma[12, a] + 542640*PolyGamma[0, a]*PolyGamma[5, a]*
      PolyGamma[12, a] + 77520*PolyGamma[6, a]*PolyGamma[12, a] +
     38760*PolyGamma[0, a]^6*PolyGamma[13, a] + 581400*PolyGamma[0, a]^4*
      PolyGamma[1, a]*PolyGamma[13, a] + 1744200*PolyGamma[0, a]^2*
      PolyGamma[1, a]^2*PolyGamma[13, a] + 581400*PolyGamma[1, a]^3*
      PolyGamma[13, a] + 775200*PolyGamma[0, a]^3*PolyGamma[2, a]*
      PolyGamma[13, a] + 2325600*PolyGamma[0, a]*PolyGamma[1, a]*
      PolyGamma[2, a]*PolyGamma[13, a] + 387600*PolyGamma[2, a]^2*
      PolyGamma[13, a] + 581400*PolyGamma[0, a]^2*PolyGamma[3, a]*
      PolyGamma[13, a] + 581400*PolyGamma[1, a]*PolyGamma[3, a]*
      PolyGamma[13, a] + 232560*PolyGamma[0, a]*PolyGamma[4, a]*
      PolyGamma[13, a] + 38760*PolyGamma[5, a]*PolyGamma[13, a] +
     15504*PolyGamma[0, a]^5*PolyGamma[14, a] + 155040*PolyGamma[0, a]^3*
      PolyGamma[1, a]*PolyGamma[14, a] + 232560*PolyGamma[0, a]*
      PolyGamma[1, a]^2*PolyGamma[14, a] + 155040*PolyGamma[0, a]^2*
      PolyGamma[2, a]*PolyGamma[14, a] + 155040*PolyGamma[1, a]*
      PolyGamma[2, a]*PolyGamma[14, a] + 77520*PolyGamma[0, a]*
      PolyGamma[3, a]*PolyGamma[14, a] + 15504*PolyGamma[4, a]*
      PolyGamma[14, a] + 4845*PolyGamma[0, a]^4*PolyGamma[15, a] +
     29070*PolyGamma[0, a]^2*PolyGamma[1, a]*PolyGamma[15, a] +
     14535*PolyGamma[1, a]^2*PolyGamma[15, a] + 19380*PolyGamma[0, a]*
      PolyGamma[2, a]*PolyGamma[15, a] + 4845*PolyGamma[3, a]*
      PolyGamma[15, a] + 1140*PolyGamma[0, a]^3*PolyGamma[16, a] +
     3420*PolyGamma[0, a]*PolyGamma[1, a]*PolyGamma[16, a] +
     1140*PolyGamma[2, a]*PolyGamma[16, a] + 190*PolyGamma[0, a]^2*
      PolyGamma[17, a] + 190*PolyGamma[1, a]*PolyGamma[17, a] +
     20*PolyGamma[0, a]*PolyGamma[18, a] + PolyGamma[19, a]

(******************************************************************************)



MySimplify[expr_] := Module[{temp, parts},
  temp = Together[expr] /. 
    Power[xxx_, yyy_] :> 
     If[(*And[MemberQ[Variables[yyy],ep],Head[Expand[xxx]]===Plus]*)
      And[IntegerQ[yyy], yyy > 0] =!= True, p[xxx, yyy], 
      Power[xxx, yyy]];
  temp = {temp};
  parts = Cases[temp, p[__], {0, Infinity}];
  While[Length[parts] > 0,
   temp = temp /. parts[[1]] -> exp1;
   
   temp = (min = -Exponent[##, 1/exp1];
       max = Exponent[##, exp1];
       (*Print[min];
       Print[max];*)
       
       Table[Coefficient[##, exp1, iii]*(exp1)^iii, {iii, min, 
         max}]) & /@ temp;
   temp = Flatten[temp, 1];
   temp = temp /. exp1 -> (parts[[1]] /. p -> pp);
   parts = Cases[temp, p[__], {0, Infinity}];
   ];
  (*temp = Together[temp];*)
temp=Together/@temp;
(*If[Length[temp]>1,Print[expr]];*)
  temp
  ]
MyFactor[xxx_] := Module[{temp,bad},
  temp=Together[xxx];
  If[Head[temp===Plus],Return[SDEvaluateDirect2[x,{temp},{1},0]]];
  temp = Select[
    Select[Cases[temp, pp[_, _]], MemberQ[Variables[##[[1]]], x[_]] &],
     Or[MemberQ[Variables[##[[2]]], ep],MemberQ[Variables[##[[2]]], delta[_]]] &];
  temp = Transpose[
     Prepend[Apply[List, temp, {1}], {Together[xxx/Times @@ temp], 
       1}]] /. pp -> Power;
  (*Print[temp];*)
  
  temp = 
    Reap[Apply[Sow, 
       Transpose[
        temp], {1}], _, {Together[Times @@ (#2)], #1} &][[2]];
  bad=Select[temp,And[MemberQ[Variables[##[[2]]],ep],Or[Head[Together[##[[1]]]]===Plus,Depth[##[[1]]]>3]]&];
  If[Length[bad]>1,Print[temp]];
  If[Length[bad]==1,
    temp=DeleteCases[temp,bad[[1]]];
    temp=Join[{temp[[1]],bad[[1]]},Drop[temp,1]];
  ];
  temp=Transpose[temp];
  (*temp*)
  Global`ParInt@@temp
(*  SDEvaluateDirect2 @@ Append[Prepend[temp, x], 0]*)
  ]







(*functions required for the new strategy*)
(*******************************************************************************)
(* QHull[pts,dim]: simplistic interface to QHull package by Alexey Pak

   Parameters:
     dim -- dimension of space (e.g., "2"),
     pts -- list of dim-dimensional points
       (e.g., "{{1,0},{0,1},{1,1},{0,0}}").
   Returns:
     list of convex hull faces, each represented by a
     list of point indices that are incident to the face
     (e.g., "{{3,2},{1,3},{2,4},{4,1}}").
   Note: necessary options from qhull manual:
   >
   > Fv - print vertices for each facet
   >
   > The first line is the number of facets. Then each facet is printed,
   > one per line. Each line is the number of vertices followed by the
   > corresponding point ids. Vertices are listed in the order they were
   > added to the hull (the last one added is the first listed).
*)

ClearAll[QHull];

Options[QHull] = {
  Verbose    -> False,
  Executable -> "./qhull",
  Mode       -> "Fv",
  InFile     -> "./_qhi.tmp",
  OutFile    -> "./_qho.tmp" };

QHull[pts_List, dim_Integer, ops___Rule] :=
  Module[{vb,ex,md,fi,fo,str,np,nf,fs,tm,i,s},

  vb = Verbose    /. {ops} /. Options[QHull];
  ex = Executable /. {ops} /. Options[QHull];
  md = Mode       /. {ops} /. Options[QHull];
  fi = InFile     /. {ops} /. Options[QHull];
  fo = OutFile    /. {ops} /. Options[QHull];

  np = Length[pts]; (* number of points *)
(*
  (* write input file *)
  str = OpenWrite[fi];
  WriteString[str,ToString[dim] <> " # dimension\n"];
  WriteString[str,ToString[np] <> " # number of points\n"];
  For[i = 1, i <= np, i++,
    s = StringJoin[Riffle[ToString /@ pts[[i]],{" "}]];
    WriteString[str, s <> "\n"]];
  Close[str];

  (* run qhull program *)
  If[vb, Print["running qhull..."]];
  tm = Timing[Run[ex <> " " <> md <> " < " <> fi <> " > " <> fo]];
  If[vb, Print["done: ",tm[[1]]," seconds"]];

  (* read qhull output *)
  str = OpenRead[fo];
  nf = Read[str,Number]; (* number of facets *)
  fs = ReadList[str,Number,RecordLists->True];
  Close[str];

  DeleteFile[fi];
  DeleteFile[fo];
*)

instr=ToString[dim] <> " # dimension\n" <> ToString[np] <> " # number of points\n";
For[i = 1, i <= np, i++,
    s = StringJoin[Riffle[ToString /@ pts[[i]],{" "}]];
    instr= instr <> s <> "\n"];

    outstr=ReadList["!echo \"" <> instr <> "\" | "<> ex <>  " "<> md, Number, RecordLists -> True];
    fs=Drop[outstr,1];
(*
Print[instr];
Print[outstr];/home/sander/FIESTA2/qhull
Abort[];
*)

  (* adjust facet notations for Mathematica: *)
  (* - remove vertex counter (first element) from every entry *)
  (* - increment vertex indices (qhull assumes first index 0) *)
  Return[Drop[# + 1,1]& /@ fs]]

(*******************************************************************************)


QHullRun[pts_List, dim_Integer]:=QHull[pts,dim,Executable->QHullPath,InFile->DataPath<>"q"<>ToString[$KernelID]<>"i",OutFile->DataPath<>"q"<>ToString[$KernelID]<>"o",Mode->"Fv 2>/dev/null",Verbose->False];

(*faces of a ploytope containing the origin*)
ZeroFaces[list_List]:=Module[{n,pos,temp},
    If[Length[list]==0,Print["Zero list"];Abort[]];
    n=Length[list[[1]]];
    pos=Position[list,Table[0,{n}]];
    If[pos=={},Print["List does not contain zero"];Abort[]];
    pos=pos[[1]][[1]];
    temp=QHullRun[list,n];
    If[temp=={},Print["Something went wrong with QHull"];Print[list];Abort[]];
    temp=Select[temp,MemberQ[##,pos]&];
    If[temp=={},Return[{}]];
    temp=Map[list[[##]]&,temp,{2}];
    Return[temp];
]

(* normal vector to a facet containing zero*)
NormVector[list_List,direction_List]:=Module[{n,pos,temp,det,vec,rank,i,j,temp2},
    If[Length[list]==0,Print["Empty facet"];Abort[]];
    n=Length[list[[1]]];
    pos=Position[list,Table[0,{n}]];
    If[pos=={},Print["Facet does not contain zero"];Abort[]];
    temp=DeleteCases[list,Table[0,{n}]];

(*    rank=0;i=1;temp2=Table[Table[0,{n}],{n}];
    While[rank<n-1,
      If[MatrixRank[ReplacePart[temp2,temp[[i]],rank+1]]>MatrixRank[temp2],
	temp2[[rank+1]]=temp[[i]];
	rank++;
      ];
      i++;
    ];
    temp=temp2;
    temp[[n]]=Table[1,{n}];
*)

      For[i=1,i<=Length[temp],i++,
	  If[temp[[i]]==Table[0,{Length[temp[[i]]]}],Continue[]];
	  pos=First[Complement[Range[Length[temp[[i]]]],Flatten[Position[temp[[i]],0]]]];
	  For[j=i+1,j<=Length[temp],j++,
	      If[temp[[j]][[pos]]!=0,temp[[j]]=(temp[[i]][[pos]])*temp[[j]]-(((temp[[j]][[pos]]))*temp[[i]])];
	  ];
      ];
      temp=DeleteCases[temp,Table[0,{Length[temp[[1]]]}]];
      AppendTo[temp,Table[1,{n}]];

    
    det=Det[temp];
    vec=Table[0,{n}];
    vec[[n]]=det;
    vec=Inverse[temp].vec;
    det=GCD@@vec;
    vec=(##/det)&/@vec;
    If[vec.direction<0,vec=-vec];
    Return[vec];
]

(*dual cone to the original polytope*)
DualCone[list_List]:=Module[{n,pos,temp,c},
    temp=ZeroFaces[list];
    If[temp=={},Return[{}]];
    n=Length[list[[1]]];
    If[Length[temp]<n,Return[{}]]; (*zero is not a facet in fact*)
    (*If[Length[temp]<n,Print["The number of facets is too small"];Print[list];Abort[];];*)
    temp=NormVector[##,DeleteCases[DeleteCases[temp,##][[1]],Table[0,{n}]][[1]]]&/@temp;
    temp={Plus@@##,##}&/@temp;
    c=LCM@@((##[[1]])&/@temp);
    temp=((c/##[[1]])*(##[[2]]))&/@temp;
    Return[temp];
]

(*qhull of a polytope of a smaller dimension*)
QHullWrapper[list_List,dim_Integer]:=Module[{n,temp,rank,i,ort,found,minor,positions,pos},
    temp=(##-list[[1]])&/@list;
    temp=Drop[temp,1];
    n=Length[temp[[1]]];
    rank=dim;
    If[rank<n,
      temp=Transpose[temp];
      For[i=1,i<=n,i++,
	  If[temp[[i]]==Table[0,{Length[temp[[i]]]}],Continue[]];
	  pos=First[Complement[Range[Length[temp[[i]]]],Flatten[Position[temp[[i]],0]]]];
	  For[j=i+1,j<=n,j++,
	      If[temp[[j]][[pos]]!=0,temp[[j]]=temp[[j]]-(((temp[[j]][[pos]])/(temp[[i]][[pos]]))*temp[[i]])];
	  ];
      ];
      temp=DeleteCases[temp,Table[0,{Length[temp[[1]]]}]];
      newtemp=Transpose[temp];
      newtemp=Prepend[newtemp,Table[0,{rank}]];


      newtemp=(LCM@@(Denominator/@Flatten[newtemp]))*newtemp;


      Return[QHullRun[newtemp,rank]];
    ,
      Return[QHullRun[list,n]];
    ]
]

(*recursive simplexification of a cone*)
Simplexify[list_List,dim_Integer]:=Module[{temp},
    If[Length[list]<=dim,Print["Simplefify error"];Abort[]];
    If[Length[list]==dim+1,Return[{list}]];
    temp=QHullWrapper[list,dim];
    temp=Select[temp,Not[MemberQ[##,1]]&];
    temp=Map[list[[##]]&,temp,{2}];
    If[STRATEGY===STRATEGY_KU2,
      temp=(SimplexifyWrapper[##,dim-1])&/@temp;
    ,
      temp=(Simplexify[##,dim-1])&/@temp;
    ];
    temp=Flatten[temp,1];
    temp=Prepend[##,list[[1]]]&/@temp;
    Return[temp];
]

SimplexifyWrapper[list_List,dim_Integer]:=Module[{temp,min,i,tempres,res},
    min=Infinity; 
    i=1;
    temp=list;
    While[And[i<=Length[temp],min>2],
	tempres=Simplexify[temp,dim];
	If[Length[tempres]<min,
	  min=Length[tempres];
	  res=tempres;
	];
	i++;
	temp=RotateRight[temp];
    ];
    Return[res];
]

(*simplex cones corresponding to the original polytope*)
SimplexCones[list_List]:=Module[{n,pos,temp,facets,i,min,res,tempres},
    temp=DualCone[list];
    If[temp=={},Return[{}]];
    n=Length[temp[[1]]];
    If[Length[temp]<n,Print["Not enough edges in the dual cone"];Abort[]];
    If[MatrixRank[temp]<n,Return[{}]]; (*zero was not the vertex*)
    If[STRATEGY===STRATEGY_KU0,
      Return[Simplexify[temp,n-1]];
    ,
      Return[SimplexifyWrapper[temp,n-1]];
    ];
(*
    min=Infinity;
    i=1;
    While[And[i<=Length[temp],min>2],
	tempres=Simplexify[temp,n-1];
	If[Length[tempres]<min,
	  min=Length[tempres];
	  res=tempres;
	];
	i++;
	temp=RotateRight[temp];
    ];
    Return[res];*)
]
    
    












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


SDEvaluateDirect[x__]:=FIESTA`SDEvaluateDirect[x]
SDExpandDirect[x__]:=FIESTA`SDExpandDirect[x]
SDEvaluate[x__]:=FIESTA`SDEvaluate[x]
SDEvaluateG[x__]:=FIESTA`SDEvaluateG[x]
SDExpand[x__]:=FIESTA`SDExpand[x]
SDExpandG[x__]:=FIESTA`SDExpandG[x]
SDAnalyze[x__]:=FIESTA`SDAnalyze[x]
UF[x__]:=FIESTA`UF[x]
ClearResults[]:=FIESTA`ClearResults[]
