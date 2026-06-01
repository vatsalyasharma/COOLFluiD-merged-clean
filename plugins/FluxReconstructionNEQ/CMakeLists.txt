LIST( APPEND FluxReconstructionNEQ_files
BCFarField.cxx
BCFarField.hh
BCFarFieldVS.cxx
BCFarFieldVS.hh
BCRankineHugoniotVS.cxx
BCRankineHugoniotVS.hh
BCSubInletSFInternalVS.cxx
BCSubInletSFInternalVS.hh
ComputeShockVelocityVS.cxx
ComputeShockVelocityVS.hh
FluxReconstructionNEQ.hh
CNEQSourceTerm.cxx
CNEQSourceTerm.hh
CNEQSourceTermVS.cxx
CNEQSourceTermVS.hh
TNEQSourceTerm.cxx
TNEQSourceTerm.hh
TNEQSourceTermVS.cxx
TNEQSourceTermVS.hh
)

LIST ( APPEND FluxReconstructionNEQ_cflibs FluxReconstructionMethod NEQ Framework)
CF_ADD_PLUGIN_LIBRARY ( FluxReconstructionNEQ )

CF_WARN_ORPHAN_FILES()
