
:Begin:
:Function:      Integrate
:Pattern:       CIntegrate[s_String]
:Arguments:     {s}
:ArgumentTypes: {String}
:ReturnType:    Manual
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

:Begin:
:Function:      setIntegrator
:Pattern:       SetIntegrator[s_String]
:Arguments:     {s}
:ArgumentTypes: {String}
:ReturnType:    Manual
:End:

:Begin:
:Function:      setCurrentIntegratorParameter
:Pattern:       SetIntegratorParameter[name_String,value_String]
:Arguments:     {name,value}
:ArgumentTypes: {String,String}
:ReturnType:    Manual
:End:

:Begin:
:Function:      getCurrentIntegratorParameters
:Pattern:       GetIntegratorParameters[]
:Arguments:     {}
:ArgumentTypes: {}
:ReturnType:    Manual
:End:

:Begin:
:Function:      setDefaultPrecision
:Pattern:       SetMPPrecision[x_Integer]
:Arguments:     {x}
:ArgumentTypes: {Integer}
:ReturnType:    Integer
:End:

:Begin:
:Function:      setMPSmallX
:Pattern:       SetSmallX[ep_Real]
:Arguments:     {ep}
:ArgumentTypes: {Real}
:ReturnType:    Integer
:End:

:Begin:
:Function:      setMPthreshold
:Pattern:       SetMPThreshold[ep_Real]
:Arguments:     {ep}
:ArgumentTypes: {Real}
:ReturnType:    Integer
:End:

:Begin:
:Function:      setMPmin
:Pattern:       SetMPMin[ep_Real]
:Arguments:     {ep}
:ArgumentTypes: {Real}
:ReturnType:    Integer
:End:

:Begin:
:Function:      setMPPrecisionShift
:Pattern:       SetMPPrecisionShift[x_Integer]
:Arguments:     {x}
:ArgumentTypes: {Integer}
:ReturnType:    Integer
:End:


:Evaluate:      CIntegrate::usage = "CIntegrate[string] performs an external integration"
:Evaluate:      CIntegrate::failed = "`1`"
:Evaluate:      SetMPPrecision::failed = "`1`"
:Evaluate:      SetCut::failed = "`1`"
:Evaluate:      SetSmallX::failed = "`1`"
:Evaluate:      SetMPThreshold::failed = "`1`"
:Evaluate:		Print["External integration ready! Use CIntegrate to perform calls"];

