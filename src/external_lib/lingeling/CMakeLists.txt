# TODO set correct CFLAGS
add_library(lgl STATIC
    "${CMAKE_SOURCE_DIR}/external/lingeling-ala/lglib.c"
    "${CMAKE_SOURCE_DIR}/external/lingeling-ala/lglbnr.c"
)
