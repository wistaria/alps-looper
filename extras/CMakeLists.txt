option(LOOPER_BUILD_EXTRAS "Build in extras subdirectory" OFF)

option(LOOPER_BUILD_EXTRAS_CORRLEN "Build in extras/corrlen subdirectory" OFF)
if(LOOPER_BUILD_EXTRAS OR LOOPER_BUILD_EXTRAS_CORRLEN)
  add_subdirectory(corrlen)
endif(LOOPER_BUILD_EXTRAS OR LOOPER_BUILD_EXTRAS_CORRLEN)

option(LOOPER_BUILD_EXTRAS_GAP "Build in extras/gap subdirectory" OFF)
if(LOOPER_BUILD_EXTRAS OR LOOPER_BUILD_EXTRAS_GAP)
  add_subdirectory(gap)
endif(LOOPER_BUILD_EXTRAS OR LOOPER_BUILD_EXTRAS_GAP)

option(LOOPER_BUILD_EXTRAS_GAP4 "Build in extras/gap4 subdirectory" OFF)
if(LOOPER_BUILD_EXTRAS OR LOOPER_BUILD_EXTRAS_GAP4)
  add_subdirectory(gap4)
endif(LOOPER_BUILD_EXTRAS OR LOOPER_BUILD_EXTRAS_GAP4)

# option(LOOPER_BUILD_EXTRAS_LOCALSUS "Build in extras/localsus subdirectory" OFF)
# if(LOOPER_BUILD_EXTRAS OR LOOPER_BUILD_EXTRAS_LOCALSUS)
#   add_subdirectory(localsus)
# endif(LOOPER_BUILD_EXTRAS OR LOOPER_BUILD_EXTRAS_LOCALSUS)

option(LOOPER_BUILD_EXTRAS_TOP "Build in extras/top subdirectory" OFF)
if(LOOPER_BUILD_EXTRAS OR LOOPER_BUILD_EXTRAS_TOP)
  add_subdirectory(top)
endif(LOOPER_BUILD_EXTRAS OR LOOPER_BUILD_EXTRAS_TOP)

option(LOOPER_BUILD_EXTRAS_TRANSMAG "Build in extras/transmag subdirectory" OFF)
if(LOOPER_BUILD_EXTRAS OR LOOPER_BUILD_EXTRAS_TRANSMAG)
  add_subdirectory(transmag)
endif(LOOPER_BUILD_EXTRAS OR LOOPER_BUILD_EXTRAS_TRANSMAG)

option(LOOPER_BUILD_EXTRAS_XYZ "Build in extras/xyz subdirectory" OFF)
if(LOOPER_BUILD_EXTRAS OR LOOPER_BUILD_EXTRAS_XYZ)
  add_subdirectory(xyz)
endif(LOOPER_BUILD_EXTRAS OR LOOPER_BUILD_EXTRAS_XYZ)
