IF ( CF_ENABLE_PLATO )
   LIST ( APPEND PlatoI_files
   Plato.hh
   PlatoLibrary.hh
   PlatoLibrary.cxx
   )
   LIST ( APPEND PlatoI_includedirs ${PLATO_INCLUDE_DIR} )
   LIST ( APPEND PlatoI_libs ${PLATO_LIBRARY} gfortran ) #Vatsalya: libplato.so (built with gfortran) did not declare libgfortran dependency; masked when PETSc loaded transitively, crashed for explicit (no PETSc) runs
   LIST ( APPEND PlatoI_cflibs Framework )		
   CF_ADD_PLUGIN_LIBRARY ( PlatoI )

   CF_WARN_ORPHAN_FILES()
ENDIF()
