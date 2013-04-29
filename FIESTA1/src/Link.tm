
:Begin:
:Function:      Integrate
:Pattern:       CIntegrate[s_String]
:Arguments:     {s}
:ArgumentTypes: {String}
:ReturnType:    Manual
:End:

:Begin:
:Function:      setpoints
:Pattern:       SetPoints[i1_Integer,i2_Integer,i11_Integer,i22_Integer]
:Arguments:     {i1,i2,i11,i22}
:ArgumentTypes: {Integer,Integer,Integer,Integer}
:ReturnType:    Integer
:End:

:Begin:
:Function:      setcut
:Pattern:       SetCut[x_Integer,ep_Real]
:Arguments:     {x,ep}
:ArgumentTypes: {Integer,Real}
:ReturnType:    Integer
:End:

:Begin:
:Function:      AddString
:Pattern:       CAddString[s_String]
:Arguments:     {s}
:ArgumentTypes: {String}
:ReturnType:    Manual
:End:

:Begin:
:Function:      ClearString
:Pattern:       CClearString[]
:Arguments:     {}
:ArgumentTypes: {}
:ReturnType:    Manual
:End:


:Evaluate:      CIntegrate::usage = "CIntegrate[string] performs an external integration"
:Evaluate:      CIntegrate::failed = "`1`"
:Evaluate:		Print["External integration ready! Use CIntegrate to perform calls"];

