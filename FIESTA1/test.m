<< FIESTA_1.0.0.m;
UsingQLink=False;
SDEvaluate[UF[{k},{-k^2,-(k+p1)^2,-(k+p1+p2)^2,-(k+p1+p2+p4)^2},
{p1^2->0,p2^2->0,p4^2->0,p1 p2->-S/2,p2 p4->-T/2,p1 p4->(S+T)/2,
S->3, T->1}],{1,1,1,1},0]
Quit;
