C     ABAQUS Subroutine for Isotropic material - Plane Stress

C---------------------------------------------------------------------------------------

C Start of Base code, DO NOT change

C---------------------------------------------------------------------------------------
      SUBROUTINE UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,
     1 RPL,DDSDDT,DRPLDE,DRPLDT,
     2 STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,
     3 NDI,NSHR,NTENS,NSTATV,PROPS,NPROPS,COORDS,DROT,PNEWDT,
     4 CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,JSTEP,KINC)
C
      INCLUDE 'ABA_PARAM.INC'
C
      CHARACTER*80 CMNAME
      DIMENSION STRESS(NTENS),STATEV(NSTATV),
     1 DDSDDE(NTENS,NTENS),DDSDDT(NTENS),DRPLDE(NTENS),
     2 STRAN(NTENS),DSTRAN(NTENS),TIME(2),PREDEF(1),DPRED(1),
     3 PROPS(NPROPS),COORDS(3),DROT(3,3),DFGRD0(3,3),DFGRD1(3,3),
     4 JSTEP(4)

C---------------------------------------------------------------------------------------

C End of Base code

C---------------------------------------------------------------------------------------

        PARAMETER (ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0)

C---------------------------------------------------------------------------------------

C Start of USER code

C---------------------------------------------------------------------------------------
C           PROPS(NPROPS)
C           User-specified array of material constants associated with this user material.
C
C           NPROPS
C           User-defined number of material constants associated with this user material.
C
C           PROPS(1) = Elasticity Modulus
C           PROPS(2) = Poisson's Ratio     
     
		   REAL*8 I1, I2, Jc, E11, E22, E12, devE11, devE22, devE12
		   REAL*8 Ci, C0, C1, C2, C3, C4, C5, eta
		   REAL*8 dh, ddh, T, T1, RPL, II

				II = STATEV(1) !Loading I1 value from last increment, it's only used for the T calc not stress
				!Strains and Invariants
				E11 = STRAN(1) + DSTRAN(1)
				E22 = STRAN(2) + DSTRAN(2)				
				E12 = 0.5*(STRAN(3) + DSTRAN(3))
				T=TEMP+DTEMP

				I1 = E11+E22

				!Deviatoric Strains
				devE11 = 0.5*(E11-E22)
				devE22 = 0.5*(E22-E11) 
				devE12 = E12
					

				!Jacobian (detF)				
				!Jc = DFGRD1(1,1)*DFGRD1(2,2) - DFGRD1(1,2)*DFGRD1(2,1)
				Jc = DFGRD1(1,1)*DFGRD1(2,2)*DFGRD1(3,3) - DFGRD1(1,2)*DFGRD1(2,1)*DFGRD1(3,3)


				!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
				! Stored Energy Density
				! Psi= Ci*I_2^2 + h(I_1)
				! i.e., stress = del Psi/del I_1 * delta + del Psi/del I_2 * (dev strain / I_2) = dh * delta + 2 * C * dev strain  (+ eta * del strain / del time)
				!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

				!Parameters
				Ci = PROPS(1)
				!
				C0 = 0.0				
				C1 = PROPS(3)
				C2 = PROPS(4)
				C3 = PROPS(5)
				C4 = PROPS(6)
				C5 = PROPS(7)

				!
				eta = PROPS(2) ! eta/DTIME

				!del Psi /del I1		
				dh = C0 + C1*T + C2*T**2.0 + C3*T**3.0 + C4*T**4.0 + C5*T**5.0

				!del^2 Psi / del I1^2	
				ddh = C1 + 2.0*C2*T + 3.0*C3*T**2.0 + 4.0*C4*T**3.0 + 5.0*C5*T**4.0

				!Stress
				STRESS(1)=(1.0/Jc) * (dh + 2.0*Ci*devE11 + eta*DSTRAN(1)/DTIME)
				STRESS(2)=(1.0/Jc) * (dh + 2.0*Ci*devE22 + eta*DSTRAN(2)/DTIME)
				STRESS(3)=(1.0/Jc) * (2.0*Ci*devE12 + eta*DSTRAN(3)/DTIME)
						   
				!Consistent Jacobian						
				DDSDDE(1,1) = (1.0/Jc) * (Ci + eta/DTIME)
				DDSDDE(1,2) = (1.0/Jc) * (- Ci) 
				DDSDDE(1,3) = 0.0
				DDSDDE(2,1) = (1.0/Jc) * (- Ci)
				DDSDDE(2,2) = (1.0/Jc) * (Ci + eta/DTIME)
				DDSDDE(2,3) = 0.0
				DDSDDE(3,1) = 0.0
				DDSDDE(3,2) = 0.0
				DDSDDE(3,3) = (1.0/Jc) * (Ci + eta/DTIME) !with respect to engineering shear strain, not shear strain
				
				!DelS/DelT
				DDSDDT(1) = 1.0/Jc*ddh
				DDSDDT(2) = 1.0/Jc*ddh
				DDSDDT(3) = 0
				
				!Heat source
				RPL = I1 - T 	! If using II "Staggered", if I1 "Monolithic"
				!delr/delT 
				DRPLDT = -1 
				DRPLDE(1) = 1  
				DRPLDE(2) = 1
				DRPLDE(3) = 0
				
				STATEV(1) = I1
				STATEV(2) = dh
				STATEV(3) = T

C---------------------------------------------------------------------------------------

C End of USER code

C---------------------------------------------------------------------------------------

      RETURN
      END	