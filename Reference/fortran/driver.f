C     driver.f: print the partonic NLO coefficient pieces of the reference
C     code (kang2.f: unpolarized, Aversa et al.; ourh_p.f: longitudinally
C     polarized, Jaeger et al.) for the channel J0 at the point (V, W),
C     partonic S, factorization scales Q2FAC (initial), Q2FRAG (final),
C     renormalization scale Q2MU and NF flavours, with Nc = 3.
C     Input lines: J0 V W S Q2FAC Q2MU Q2FRAG NF
C     Output: F0(J0), AVDELMU, AVDEL, AVWPL(W), AVWPL(1), AVLO(W), AVLO(1),
C             AVGO(W,V)+STRUV(W,V,X3,S), CC(J0)
      PROGRAM DRIVER
      IMPLICIT REAL*8 (A-H,L-Z)
      INTEGER J0, NFI, IOS
      COMMON /PARAMS/ CF,NC,NF,PI,HC2
      COMMON /CONS/ PPI,S,GV,GW,NNC,GTR,CCF,VC
      COMMON /ORDE/ AL,CQ,V1,V2,V3,V4
      COMMON /EDF/ J0
      COMMON /SCALES/ Q2FAC,Q2MU,Q2FRAG
      DIMENSION F0(16)
      PI=3.141592653589793238D0
      NC=3.D0
      CF=4.D0/3.D0
      HC2=1.D0
 10   READ(*,*,IOSTAT=IOS) J0,V,W,S,Q2FAC,Q2MU,Q2FRAG,NFI
      IF (IOS.NE.0) STOP
      NF=DBLE(NFI)
      PPI=PI
      NNC=NC
      CCF=CF
      AL=1.D0
      CQ=0.D0
      VC=NNC**2-1.D0
      V1=VC**2/NNC
      V2=VC/NNC
      V3=(NNC**4-1.D0)/2.D0/NNC**2
      V4=VC**2/2.D0/NNC**2
      GTR=NF/2.D0
      GV=V
      GW=W
      UN=1.D0
      IF (J0.EQ.16.OR.J0.EQ.15) THEN
        CC=8.D0*VC**2
      ELSE IF (J0.EQ.14.OR.J0.EQ.13.OR.J0.EQ.10.OR.J0.EQ.9.OR.J0.EQ.8)
     #  THEN
        CC=8.D0*VC*NNC
      ELSE
        CC=8.D0*NNC**2
      ENDIF
      CALL FBOR(V,S,F0)
      X3=1.D0
      WRITE(*,'(I3,9(1X,1PE23.15))') J0, F0(J0), AVDELMU(V,S,J0),
     #  AVDEL(V,S), AVWPL(W,V,S), AVWPL(UN,V,S), AVLO(W,V,S),
     #  AVLO(UN,V,S), AVGO(W,V)+STRUV(W,V,X3,S), CC
      GOTO 10
      END
