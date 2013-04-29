NumberOfLinks=8;
IfCut=0.1; 
VarExpansionDegree=0;
STRATEGY=STRATEGY_S;
Get["/users/guest/asmirnov/FIESTA/FIESTA_1.2.0.m"];
VegasSettings={10000,5,100000,15};
NumberOfSubkernels=8;
CIntegratePath="/users/guest/asmirnov/FIESTA/CIntegrate64";
QLinkPath="/users/guest/asmirnov/FIESTA/QLink64";
SDEvaluate[
UF[{k,l,r},{-(k + q)^2, -(r + q)^2, -r^2, -l^2, -k^2, -(k - r)^2, -(l - r)^2, 
-(k - l)^2, -v k, -v  r, -v (k - l), -(l + q)^2},{q^2->-QQ,v^2->vv,q*v->0,
QQ->1,vv->1}],
{1, 0, 0, 1, -1, 1, 1, 0, 0, 1, 1, 0},2]
Exit[]

