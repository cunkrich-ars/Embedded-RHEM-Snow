MODULE snow

  USE, INTRINSIC::ISO_C_BINDING

  IMPLICIT NONE

  INTERFACE
    INTEGER (C_LONG) FUNCTION snow_run(id, cligen_file, soil_class, slope, aspect) BIND(C)
      USE, INTRINSIC::ISO_C_BINDING
      CHARACTER, DIMENSION(*) :: id
      CHARACTER, DIMENSION(*) :: cligen_file
      CHARACTER, DIMENSION(*) :: soil_class
      REAL (C_DOUBLE), VALUE :: slope
      REAL (C_DOUBLE), VALUE :: aspect
    END FUNCTION snow_run
  END INTERFACE

  INTERFACE
    INTEGER (C_LONG) FUNCTION snow_next_event( date ) BIND(C)
      USE, INTRINSIC::ISO_C_BINDING
      TYPE (C_PTR), VALUE :: date ! INTEGER (C_LONG)
    END FUNCTION snow_next_event
  END INTERFACE

  INTERFACE
    INTEGER (C_LONG) FUNCTION snow_npoints() BIND(C)
      USE, INTRINSIC::ISO_C_BINDING
    END FUNCTION snow_npoints
  END INTERFACE

  INTERFACE
    INTEGER (C_LONG) FUNCTION snow_times( t ) BIND(C)
      USE, INTRINSIC::ISO_C_BINDING
      TYPE (C_PTR), VALUE :: t ! REAL (C_DOUBLE)
    END FUNCTION snow_times
  END INTERFACE

  INTERFACE
    INTEGER (C_LONG) FUNCTION snow_depths( d ) BIND(C)
      USE, INTRINSIC::ISO_C_BINDING
      TYPE (C_PTR), VALUE :: d ! REAL (C_DOUBLE)
    END FUNCTION snow_depths
  END INTERFACE

  INTERFACE
    REAL (C_DOUBLE) FUNCTION snow_sat() BIND(C)
      USE, INTRINSIC::ISO_C_BINDING
    END FUNCTION snow_sat
  END INTERFACE

  INTERFACE
    REAL (C_DOUBLE) FUNCTION snow_ice() BIND(C)
      USE, INTRINSIC::ISO_C_BINDING
    END FUNCTION snow_ice
  END INTERFACE

  INTERFACE
    INTEGER (C_LONG) FUNCTION snow_finale() BIND(C)
      USE, INTRINSIC::ISO_C_BINDING
    END FUNCTION snow_finale
  END INTERFACE

END MODULE snow
