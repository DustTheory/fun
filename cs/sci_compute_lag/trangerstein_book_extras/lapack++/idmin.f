c adapted from LAPACK function IDAMAX
      INTEGER FUNCTION IDMIN(N,DX,INCX)
      INTEGER INCX,N
      DOUBLE PRECISION DX(*)
      DOUBLE PRECISION DMIN
      INTEGER I,IX
      IDMIN = 0
      IF (N.LT.1 .OR. INCX.LE.0) RETURN
      IDMIN = 1
      IF (N.EQ.1) RETURN
      IF (INCX.EQ.1) THEN
        DMIN = DX(1)
        DO I = 2,N
          IF (DX(I).LT.DMIN) THEN
            IDMIN = I
            DMIN = DX(I)
          END IF
        END DO
      ELSE
        IX = 1
        DMIN = DX(1)
        IX = IX + INCX
        DO I = 2,N
          IF (DX(IX).LT.DMIN) THEN
            IDMIN = I
            DMIN = DX(IX)
          END IF
          IX = IX + INCX
        END DO
      END IF
      RETURN
      END
