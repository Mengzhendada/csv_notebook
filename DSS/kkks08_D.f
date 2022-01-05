c ... kkks08_D.f (modified by TRK, origninal version by IJS bkk05_d3.f)
c ... 26.03.2007 
c ... packed by H. Spiesberger, 12.07.2010

c ... Thresholds at Heavy Quark Masses
c ... NLO FFs for D^*, D^0, D^+ with and without mass effects
c ... extracted from OPAL, ALEPH, BELLE and CLEO

********************************************************************
*   Returns Fragmentation Functions D_i(x,Q^2) (NOT x*FFs !)       * 
*   for D^*,D^0,D^+                                                *
*                                                                  *
*   NOTE: D^* = D^*+ + D^*- ; D^+ + D^-                            *
*                                                                  *
*   INPUT:                                                         *
* ... NO = 1 (inlcudes mass effects), 2 (without mass effects)     *
* ... NH = Global: 1 KKK06 D^*, 2 KKK06 D^0, 3 KKK06 D^+           * 
*          BELLE:  4 KKK06 D^*, 5 KKK06 D^0, 6 KKK06 D^+           *
*          OPAL: 7 KKK06 D^*, 8 KKK06 D^0, 9 KKK06 D^+             *
*          ALEPH: 10 KKK06 D^*                                     *
*          CLEO: 11 KKK06 D^*, 12 KKK06 D^0, 13 KKK06 D^+          *
*     Note: CLEO FF only availabe with mass effects!               *
*                                                                  *
*   OUTPUT: DH(0:10) <-> D_i(x,Q^2) (i=0,...,10)                   *
* ... DH: 0    1    2    3    4    5    6    7    8    9    10     *
* ...     g    u  ubar   d  dbar   s  sbar   c  cbar   b   bbar    *
*                                                                  *
* ... Note that: D_u=D_d=D_s; D_i = D_ibar                         *
*                                                                  * 
*            x  = frag. variable  (between  1.E-4  and  1)         *
*            Q2 = scale in GeV**2 (between  2.25 and  1.E6)        *
*                                                                  *
*   last change: 26.03.2007                                        *
********************************************************************

      subroutine kkks08_D(NO,NH,x,Q2,DH)
      implicit none
      integer NO,NH   ! Input
      double precision X,Q2     ! Input
      double precision DH(0:10) ! Output

      integer NOold,NHold  
      SAVE NOold, NHold
c      
      integer m,n,i
      integer IQ,IX
      integer npart,nx,nq,narg
      parameter(NPART=4, NX=131, NQ=29, NARG=2)
      double precision XB0,XB1
      double precision PARTON(0:NPART-1,NQ,NX-1), QS(NQ), XB(NX)
      integer NA(NARG)
      double precision XT(NARG), ARRF(NX+NQ) 
      double precision XDG(NX,NQ),XDU(NX,NQ)!,XDD(NX,NQ),XDS(NX,NQ)
      double precision XDC(NX,NQ),XDB(NX,NQ)
      character*3 chno
      character*9 chnh
      integer fragini
      SAVE fragini
c      COMMON / inifrag / fragini
      SAVE XDU, XDG, XDC, XDB, NA, ARRF
      double precision FF_FINT
      double precision XBG,XBU,XBC,XBB
      double precision xg,xu,xc,xbot
      double precision rg,ru,rc,rb
      double precision rg1,ru1,rc1,rb1
      SAVE rg, ru, rc, rb, rg1, ru1, rc1, rb1

      double precision XB1G,XB1U,XB1C,XB1B
      double precision xg1,xu1,xc1,xbot1

c ... NQ=29
*...X AND Q**2 VALUES OF THE GRID :
       DATA QS / 2.25d0, 3d0, 5d0, 7d0, 10D0, 16D0, 25D0, 40D0, 
     #           64D0, 81d0,
     #           1.0D2, 1.5D2, 2.0D2, 2.5D2, 5.0d2,
     #           1.0D3, 1.8D3, 3.2D3, 5.7D3, 8315.18d0, 
     #           1.0D4, 1.8D4, 3.2D4, 5.7D4, 1.0d5,
     #           1.8D5, 3.2D5, 5.7D5, 1.0d6/

c ... NX=131
       DATA XB / 1.0D-4,1.4D-4,2.0D-4,3.0D-4,4.5D-4,6.7D-4,8.0d-4,
     #           1.0D-3,1.4D-3,2.0D-3,3.0D-3,4.5D-3,6.7D-3,8.0d-3,
     #           1.0D-2, 1.4D-2, 2.0D-2, 3.0D-2, 4.0d-2,
     #0.046d0,0.050d0,0.054d0,0.059d0,0.064d0,0.068d0,
     #0.072d0,0.076d0,0.080d0,0.085d0,0.09d0,0.095d0,0.1d0,0.105d0,
     #0.110d0,0.117d0,0.125d0,0.130d0,0.135d0,0.140d0,0.15d0,0.16d0,
     #0.17d0,0.18d0,0.19d0,0.2d0,0.21d0,0.22d0,0.23d0,0.24d0,
     #0.246d0,0.252d0,0.260d0,0.2660,0.272d0,0.280d0,0.29d0,0.3d0,
     #0.31d0,0.32d0,0.33d0,0.34d0,0.35d0,0.36d0,0.37d0,
     #0.38d0,0.39d0,0.40d0,0.41d0,0.42d0,0.43d0,
     #0.44d0,0.45d0,0.46d0,0.47d0,0.48d0,0.49d0,0.5d0,
     #0.51d0,0.52d0,0.53d0,0.54d0,0.55d0,0.56d0,
     #0.57d0,0.58d0,0.59d0,0.6d0,0.61d0,0.62d0,0.63d0,
     #0.64d0,0.65d0,0.66d0,0.67d0,0.68d0,0.69d0,
     #0.70d0,0.71d0,0.72d0,0.73d0,0.74d0,0.75d0,0.76d0,
     #0.77d0,0.78d0,0.79d0,0.80d0,0.81d0,0.82d0,
     #0.83d0,0.84d0,0.85d0,0.86d0,0.87d0,0.888d0,0.89d0,
     #0.90d0,0.91d0,0.92d0,0.93d0,0.94d0,0.95d0,0.96d0,
     #0.965d0,0.97d0,0.975d0,0.98d0,0.985d0,0.99d0,0.995d0,0.9999d0 /

      integer iroutx1,iroutx2,iroutq1,iroutq2
      data iroutx1/0/,iroutx2/0/
      data iroutq1/0/,iroutq2/0/

*...CHECK OF X AND Q2 VALUES : 
       IF ( (X.LT.0.99E-4) ) THEN
         iroutx1=iroutx1+1
         if (iroutx1.le.20) then
           WRITE(6,91) 
  91     FORMAT (2X,'PARTON INTERPOLATIO: in kkks08_D x TOO SMALL')
         endif
         if (iroutx1.eq.20) then
           write(6,92)
  92     format (2x,'x TOO SMALL in kkks08_D has occured 20 times',
     .              /,'no further warnings printed')
         endif
         X=0.99E-4
       ENDIF
       IF ( (X.GT.1.0) ) THEN
         iroutx2=iroutx2+1
         if (iroutx2.le.20) then
           WRITE(6,93) 
  93     FORMAT (2X,'PARTON INTERPOLATION in kkks08_D: x TOO LARGE')
         endif
         if (iroutx2.eq.20) then
           write(6,94)
  94     format (2x,'x TOO LARGE in kkks08_D has occured 20 times',
     .              /,'no further warnings printed')
         endif
         X=1.0
       ENDIF
       IF ( (Q2.LT.2.25d0) ) THEN
         iroutq1=iroutq1+1
         if (iroutq1.le.20) then
           WRITE(6,95) 
  95     FORMAT (2X,'PARTON INTERPOLATION in kkks08_D: Q2 TOO SMALL')
         endif
         if (iroutq1.eq.20) then
           write(6,94)
  96     format (2x,'Q2 TOO SMALL in kkks08_D has occured 20 times',
     .              /,'no further warnings printed')
         endif 
         Q2=2.25d0
       ENDIF
       IF ( (Q2.GT.1d4) ) THEN
         iroutq2=iroutq2+1
         if (iroutq2.le.20) then
           WRITE(6,97) 
  97     FORMAT (2X,'PARTON INTERPOLATION in kkks08_D: Q2 TOO LARGE')
         endif
         if (iroutq2.eq.20) then
           write(6,98)
  98     format (2x,'Q2 TOO LARGE in kkks08_D has occured 20 times',
     .              /,'no further warnings printed')
         endif 
       ENDIF

*
*...INITIALIZATION :
*    SELECTION AND READING OF THE GRID : 
c check whether parameters NO or NH changed
c or this is a the first call (after the firtst call fragini is set to
c the arbitrary number 12345)
c      IF (fragini.NE.0) GOTO 16
      IF (fragini.EQ.12345 .AND. NO.EQ.NOold .AND. NH.EQ.NHold) GOTO 16

      if ((NO.EQ.1) .and. (NH.EQ.1)) then ! Global with mass D^*
               chno='mas'
               chnh='glob_d+st'

               rg = -1d0
               ru = -2d0
               rc = 1d0
               rb = 1d0

               rg1 = 0d0
               ru1 = 2d0
               rc1 = 2d0
               rb1 = 3d0

      elseif ((NO.EQ.1) .and. (NH.EQ.2)) then ! Global with mass D^0         
               chno='mas'
               chnh='global_d0'

               rg = -2d0
               ru = -2d0
               rc = 0d0
               rb = 1d0

               rg1 = 1d0
               ru1 = 2d0
               rc1 = 2d0
               rb1 = 3d0

      elseif ((NO.EQ.1) .and. (NH.EQ.3)) then ! Global with mass D^+
               chno='mas'
               chnh='global_d+'

               rg = -1d0
               ru = -2d0
               rc = 0d0
               rb = 1d0

               rg1 = 0d0
               ru1 = 0d0
               rc1 = 2d0
               rb1 = 3d0


      elseif ((NO.EQ.1) .and. (NH.EQ.4)) then ! BELLE fit with mass D^*
               chno='mas'
               chnh='bell_d+st'

               rg = -1d0
               ru = -2d0
               rc = 1d0
               rb = 1d0

               rg1 = 0d0
               ru1 = 2d0
               rc1 = 2d0
               rb1 = 3d0

      elseif ((NO.EQ.1) .and. (NH.EQ.5)) then ! BELLE fit with mass D^0
               chno='mas'
               chnh='belle_d0_'

               rg = -2d0
               ru = -2d0
               rc = 0d0
               rb = 1d0

               rg1 = 1d0
               ru1 = 2d0
               rc1 = 2d0
               rb1 = 3d0

      elseif ((NO.EQ.1) .and. (NH.EQ.6)) then ! BELLE fit with mass D^+
               chno='mas'
               chnh='belle_d+_'

               rg = -1d0
               ru = -2d0
               rc = 0d0
               rb = 1d0

               rg1 = 0d0
               ru1 = 0d0
               rc1 = 2d0
               rb1 = 3d0

      elseif ((NO.EQ.1) .and. (NH.EQ.7)) then ! OPAL fit with mass D^*
               chno='mas'
               chnh='opal_d+st'

               rg = -1d0
               ru = -2d0
               rc = 1d0
               rb = 1d0

               rg1 = 0d0
               ru1 = 2d0
               rc1 = 2d0
               rb1 = 3d0

      elseif ((NO.EQ.1) .and. (NH.EQ.8)) then ! OPAL fit with mass D^0
               chno='mas'
               chnh='opal_d0__'

               rg = -2d0
               ru = -2d0
               rc = 0d0
               rb = 1d0

               rg1 = 1d0
               ru1 = 2d0
               rc1 = 2d0
               rb1 = 3d0

      elseif ((NO.EQ.1) .and. (NH.EQ.9)) then ! OPAL fit with mass D^+
               chno='mas'
               chnh='opal_d+__'

               rg = -1d0
               ru = -2d0
               rc = 0d0
               rb = 1d0

               rg1 = 0d0
               ru1 = 0d0
               rc1 = 2d0
               rb1 = 3d0

      elseif ((NO.EQ.1) .and. (NH.EQ.10)) then ! ALEPH fit with mass D^*
               chno='mas'
               chnh='alep_d+st'

               rg = -1d0
               ru = -2d0
               rc = 1d0
               rb = 1d0

               rg1 = 0d0
               ru1 = 2d0
               rc1 = 2d0
               rb1 = 3d0

      elseif ((NO.EQ.1) .and. (NH.EQ.11)) then ! CLEO fit with mass D^*
               chno='mas'
               chnh='cleo_d+st'
               stop

               rg = -1d0
               ru = -2d0
               rc = 1d0
               rb = 1d0

               rg1 = 0d0
               ru1 = 2d0
               rc1 = 2d0
               rb1 = 3d0

      elseif ((NO.EQ.1) .and. (NH.EQ.12)) then ! CLEO fit with mass D^0
               chno='mas'
               chnh='cleo_d0__'

               rg = -2d0
               ru = -2d0
               rc = 0d0
               rb = 1d0

               rg1 = 1d0
               ru1 = 2d0
               rc1 = 2d0
               rb1 = 3d0

      elseif ((NO.EQ.1) .and. (NH.EQ.13)) then ! CLEO fit with mass D^+
               chno='mas'
               chnh='cleo_d+__'
               stop

               rg = -1d0
               ru = -2d0
               rc = 0d0
               rb = 1d0

               rg1 = 0d0
               ru1 = 0d0
               rc1 = 2d0
               rb1 = 3d0

      elseif ((NO.EQ.2) .and. (NH.EQ.1)) then ! Global without mass D^*
               chno='m00'
               chnh='glob_d+st'

               rg = -1d0
               ru = -2d0
               rc = 1d0
               rb = 1d0

               rg1 = 0d0
               ru1 = 2d0
               rc1 = 2d0
               rb1 = 3d0

      elseif ((NO.EQ.2) .and. (NH.EQ.2)) then ! Global without mass D^0         
               chno='m00'
               chnh='global_d0'

               rg = -2d0
               ru = -2d0
               rc = 0d0
               rb = 1d0

               rg1 = 1d0
               ru1 = 2d0
               rc1 = 2d0
               rb1 = 3d0

      elseif ((NO.EQ.2) .and. (NH.EQ.3)) then ! Global without mass D^+
               chno='m00'
               chnh='global_d+'

               rg = -1d0
               ru = -2d0
               rc = 0d0
               rb = 1d0

               rg1 = 0d0
               ru1 = 0d0
               rc1 = 2d0
               rb1 = 3d0


      elseif ((NO.EQ.2) .and. (NH.EQ.4)) then ! BELLE fit without mass D^*
               chno='m00'
               chnh='bell_d+st'

               rg = -1d0
               ru = -2d0
               rc = 1d0
               rb = 1d0

               rg1 = 0d0
               ru1 = 2d0
               rc1 = 2d0
               rb1 = 3d0

      elseif ((NO.EQ.2) .and. (NH.EQ.5)) then ! BELLE fit without mass D^0
               chno='m00'
               chnh='belle_d0_'

               rg = -2d0
               ru = -2d0
               rc = 0d0
               rb = 1d0

               rg1 = 1d0
               ru1 = 2d0
               rc1 = 2d0
               rb1 = 3d0

      elseif ((NO.EQ.2) .and. (NH.EQ.6)) then ! BELLE fit without mass D^+
               chno='m00'
               chnh='belle_d+_'

               rg = -1d0
               ru = -2d0
               rc = 0d0
               rb = 1d0

               rg1 = 0d0
               ru1 = 0d0
               rc1 = 2d0
               rb1 = 3d0

      elseif ((NO.EQ.2) .and. (NH.EQ.7)) then ! OPAL fit without mass D^*
               chno='m00'
               chnh='opal_d+st'

               rg = -1d0
               ru = -2d0
               rc = 1d0
               rb = 1d0

               rg1 = 0d0
               ru1 = 2d0
               rc1 = 2d0
               rb1 = 3d0

      elseif ((NO.EQ.2) .and. (NH.EQ.8)) then ! OPAL fit without mass D^0
               chno='m00'
               chnh='opal_d0__'

               rg = -2d0
               ru = -2d0
               rc = 0d0
               rb = 1d0

               rg1 = 1d0
               ru1 = 2d0
               rc1 = 2d0
               rb1 = 3d0

      elseif ((NO.EQ.2) .and. (NH.EQ.9)) then ! OPAL fit without mass D^+
               chno='m00'
               chnh='opal_d+__'

               rg = -1d0
               ru = -2d0
               rc = 0d0
               rb = 1d0

               rg1 = 0d0
               ru1 = 0d0
               rc1 = 2d0
               rb1 = 3d0

      elseif ((NO.EQ.2) .and. (NH.EQ.10)) then ! ALEPH fit without mass D^*
               chno='m00'
               chnh='alep_d+st'

               rg = -1d0
               ru = -2d0
               rc = 1d0
               rb = 1d0

               rg1 = 0d0
               ru1 = 2d0
               rc1 = 2d0
               rb1 = 3d0

      elseif ((NO.EQ.2) .and. (NH.EQ.11)) then ! CLEO fit without mass D^*
               chno='m00'
               chnh='cleo_d+st'
               stop

               rg = -1d0
               ru = -2d0
               rc = 1d0
               rb = 1d0

               rg1 = 0d0
               ru1 = 2d0
               rc1 = 2d0
               rb1 = 3d0

      elseif ((NO.EQ.2) .and. (NH.EQ.12)) then ! CLEO fit without mass D^0
               chno='m00'
               chnh='cleo_d0__'
               stop

               rg = -2d0
               ru = -2d0
               rc = 0d0
               rb = 1d0

               rg1 = 1d0
               ru1 = 2d0
               rc1 = 2d0
               rb1 = 3d0

      elseif ((NO.EQ.2) .and. (NH.EQ.13)) then ! CLEO fit without mass D^+
               chno='m00'
               chnh='cleo_d+__'
               stop

               rg = -1d0
               ru = -2d0
               rc = 0d0
               rb = 1d0

               rg1 = 0d0
               ru1 = 0d0
               rc1 = 2d0
               rb1 = 3d0


      else
         print*,'NO =',NO,' and ','NH =',NH,' not supported'
         stop
      end if    

      fragini = 12345
      NOold = NO
      NHold = NH
c
      OPEN(UNIT=8,FILE=chnh//'_'//chno//'.grid',STATUS='old')

            DO M = 1, NX-1  ! Loop over x
               DO N = 1, NQ ! Loop over Q^2

               READ(8,90) PARTON(0,N,M), PARTON(1,N,M), PARTON(2,N,M), 
     #                    PARTON(3,N,M)

               enddo ! N-loop
            enddo ! M-loop
  90  FORMAT (6(1PE13.5,1X))

      CLOSE(8)
*
*....ARRAYS FOR THE INTERPOLATION SUBROUTINE :
      DO 10 IQ = 1, NQ
      DO 20 IX = 1, NX-1
c      DO 20 IX = 1, NX
        XB0 = XB(IX) 
        XB1 = 1.-XB0

        XBG = XB0**rg
        XBU = XB0**ru
        XBC = XB0**rc
        XBB = XB0**rb

        XB1G = XB1**rg1
        XB1U = XB1**ru1
        XB1C = XB1**rc1
        XB1B = XB1**rb1

        XDG(IX,IQ) = PARTON(0,IQ,IX) / (XBG * XB1G)
        XDU(IX,IQ) = PARTON(1,IQ,IX) / (XBU * XB1U)
        XDC(IX,IQ) = PARTON(2,IQ,IX) / (XBC * XB1C)
        XDB(IX,IQ) = PARTON(3,IQ,IX) / (XBB * XB1B)
  20  CONTINUE
        XDG(NX,IQ) = 0.E0
        XDU(NX,IQ) = 0.E0
        XDC(NX,IQ) = 0.E0
        XDB(NX,IQ) = 0.E0
  10  CONTINUE  
      NA(1) = NX
      NA(2) = NQ
      DO 30 IX = 1, NX
        ARRF(IX) = DLOG(XB(IX))
  30  CONTINUE
      DO 40 IQ = 1, NQ
        ARRF(NX+IQ) = DLOG(QS(IQ))
  40  CONTINUE
  16  CONTINUE
*
*...INTERPOLATION :
      XT(1) = DLOG(X)
      XT(2) = DLOG(Q2)

      xg = x**rg
      xu = x**ru
      xc = x**rc
      xbot = x**rb

      xg1 = (1.-x)**rg1
      xu1 = (1.-x)**ru1
      xc1 = (1.-x)**rc1
      xbot1 = (1.-x)**rb1

c ... flavour symmetric: D_u=D_d=D_s
c ... charge symmetry: D_i = D_anti-i

      DH(0) = FF_FINT(NARG,XT,NA,ARRF,XDG) * xg * xg1      ! D_g
      DH(1) = FF_FINT(NARG,XT,NA,ARRF,XDU) * xu * xu1      ! D_u
      DH(2) = DH(1)                                     ! D_ubar
      DH(3) = DH(1)                                     ! D_d
      DH(4) = DH(3)                                     ! D_dbar
      DH(5) = DH(1)                                     ! D_s
      DH(6) = DH(5)                                     ! D_sbar
      DH(7) = FF_FINT(NARG,XT,NA,ARRF,XDC) * xc * xc1      ! D_c
      DH(8) = DH(7)                                     ! D_cbar
      DH(9) = FF_FINT(NARG,XT,NA,ARRF,XDB) * xbot * xbot1  ! D_b
      DH(10)= DH(9)                                     ! D_bbar

c ... to avoid crash for scales smaller than q2 data/         
      if (Q2.lt.2.25d0) then
c         write(6,*) 'q**2 must be > 2.25 GeV^2',qp2
         do i=0,10
            DH(i)=0d0
         enddo   
      return
      endif

      if (x.eq.1.d0) then
         do i=0,10
            DH(i)=0d0
         enddo   
      return
      endif

      RETURN
      END
