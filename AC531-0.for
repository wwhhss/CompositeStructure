      SUBROUTINE UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,
     1 RPL,DDSDDT,DRPLDE,DRPLDT,STRAN,DSTRAN,
     2 TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,MATERL,NDI,NSHR,NTENS,
     3 NSTATV,PROPS,NPROPS,COORDS,DROT,PNEWDT,CELENT,
     4 DFGRD0,DFGRD1,NOEL,NPT,KSLAY,KSPT,KSTEP,KINC)
C
      INCLUDE 'ABA_PARAM.INC'
C
      CHARACTER*80 MATERL
      DIMENSION STRESS(NTENS),STATEV(NSTATV),
     1 DDSDDE(NTENS,NTENS),DDSDDT(NTENS),DRPLDE(NTENS),
     2 STRAN(NTENS),DSTRAN(NTENS),TIME(2),PREDEF(1),DPRED(1),
     3 PROPS(NPROPS),COORDS(3),DROT(3,3),
     4 DFGRD0(3,3),DFGRD1(3,3)
C
      !real*8:: TTdeta,alpha22,alpha11,alpha33
      !real*8:: epsilongT(6)
C
      PARAMETER (ONE=1.0D0,TWO=2.0D0,THREE=3.0D0,SIX=6.0D0)
      DATA NEWTON,TOLER/10,1.D-6/
C
C -----------------------------------------------------------
C
      IF (NDI.NE.3) THEN
         WRITE(6,1)
1       FORMAT(//,30X,'***ERROR - THIS UMAT MAY ONLY BE USED FOR ',
     1          'ELEMENTS WITH THREE DIRECT STRESS COMPONENTS')
      ENDIF
C
C     ELASTIC STIFFNESS
C
      DO K1=1,6
        DO K2=1,6
           DDSDDE(K2,K1)=0.0
        ENDDO
      ENDDO
C
      Er = 5908
      Gr = 1.233
      Vf = 0.57
      E1f = 294000
      G12f = 15000
      G23f = 5500
      v12f = 0.2
      vr = 0.38
!
      kf = E1f / (2 * (1 - v12f - 2 * v12f**2))
      kr = Er / (2 * (1 - vr - 2 * vr**2))
      kt = ((kf + Gr) * kr + (kf - kr) * Gr * Vf) / ((kf + Gr) - (kf - 
     & kr) * Vf)
      G12 = Gr * ((G12f + Gr) + (G12f - Gr) * Vf) / ((G12f + Gr) - (G12f
     & - Gr) * Vf)
      G13 = G12
      G23 = Gr * (kr * (Gr + G23f) + 2 * G23f * Gr + kr * (G23f - Gr) * 
     & Vf) / (kr * (Gr + G23f) + 2 * G23f * Gr - (kr + 2 * Gr) * (G23f -
     & Gr) * Vf)
      v12 = v12f * Vf + vr * (1 - Vf) + ((vr - v12f) * (kr - kf) * Gr * 
     & (1 - Vf) * Vf) / ((kf + Gr) * kr + (kf - kr) * Gr * Vf)
      v13 = v12
      E11 = E1f * Vf + Er * (1 - Vf) + (4 * (vr - v12f**2) * kf * kr * 
     & Gr * (1 - Vf) * Vf) / ((kf + Gr) * kr + (kf - kr) * Gr * Vf)
      E22 = 1 / ((4 * kt)**(-1) + (4 * G23)**(-1) + (v12**2 / E11))
      E33 = E22
      v23 = (2 * E11 * kt - E11 * E22 - 4 * v12**2 * kt * E22) / (2 * 
     & E11 * kt)
      v21 = v12 * E22 / E11
      v31 = v13 * E33 / E11
      v32 = v23 * E33 / E22
C
C     Insert your calculations for the material behavior here
C
      DETA = (1 - V12 * V21 - V23 * V32 - V13 * V31 - 2 * V21 * V32 * 
     & V13) / (E11 * E22 * E33)
      DDSDDE(1,1) = (1 - V23 * V32) / (E22 * E33 * DETA)
      DDSDDE(2,2) = (1 - V13 * V31) / (E11 * E33 * DETA)
      DDSDDE(3,3) = (1 - V12 * V21) / (E11 * E22 * DETA)
      DDSDDE(2,1) = (V21 + V23 * V31) / (E22 * E33 * DETA)
      DDSDDE(3,1) = (V31 + V21 * V32) / (E22 * E33 * DETA)
      DDSDDE(3,2) = (V32 + V12 * V31) / (E11 * E33 * DETA)
      DDSDDE(4,4) = G12
      DDSDDE(5,5) = G13
      DDSDDE(6,6) = G23
      DDSDDE(1,2) = DDSDDE(2,1)
      DDSDDE(1,3) = DDSDDE(3,1)
      DDSDDE(2,3) = DDSDDE(3,2)
      
       
       if(1.LE.TIME(2)) then 
       DDSDDE(1,1)= 175252.70779471454
       DDSDDE(2,1)= 9559.064071797018
       DDSDDE(3,1)= 9559.064071797018 
       DDSDDE(1,2)= DDSDDE(2,1)
       DDSDDE(2,2)= 20903.1459946621
       DDSDDE(3,2)= 14708.015078008937
       DDSDDE(1,3)= DDSDDE(3,1)
       DDSDDE(2,3)= DDSDDE(3,2)
       DDSDDE(3,3)= 20903.1459946621 
       DDSDDE(4,4)= 4565.714285714285
       DDSDDE(5,5)= 4565.714285714285
       DDSDDE(6,6)= 3097.5654583265814
       Endif

C
C    CALCULATE STRESS FROM ELASTIC STRAINS
C
      DO  K1=1,6
      DO  K2=1,6
       STRESS(K2)=STRESS(K2)+DDSDDE(K2,K1)*(DSTRAN(K1))
      ENDDO
      ENDDO
C
C
      !if(TIME(2) .LE.2) then
      !print *,DDSDDE(1,1)
      !endif
      
      RETURN
      END