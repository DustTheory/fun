      INTEGER FUNCTION ISMIN(N,SX,INCX)
      INTEGER INCX,N
      REAL SX(*)
      REAL SMIN
      INTEGER I,IX
      ISMIN = 0
      IF (N.LT.1 .OR. INCX.LE.0) RETURN
      ISMIN = 1
      IF (N.EQ.1) RETURN
      IF (INCX.EQ.1) THEN
        SMIN = SX(1)
        DO I = 2,N
          IF (SX(I).LT.SMIN) THEN
            ISMIN = I
            SMIN = SX(I)
          END IF
        END DO
      ELSE
        IX = 1
        SMIN = SX(1)
        IX = IX + INCX
        DO I = 2,N
          IF (SX(IX).LT.SMIN) THEN
            ISMIN = I
            SMIN = SX(IX)
          END IF
          IX = IX + INCX
        END DO
      END IF
      RETURN
      END
