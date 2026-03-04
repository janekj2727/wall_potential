c*******************************************************************
c*                                                                 *
c* Module wall_sph                                                 *
c* - external potential of a sphere generating constant density    *
c*   profile                                                       *
c* - classical dft with original rosenfeld's fmt functional        *
c* - provides analytical expression for the individual parts of    *
c*   the direct correlation function c(1) for hard spheres         *
c*   assuming constant density                                     *
c* - attractive interactions in the fluid due to cut lj potential  *
c*   can be assumed by adding analytical expression for the approp.*
c*   contribution to the direct correlation function c(1)_att      *
c* - polyroots module by Jacob Williams is used to calculate roots *
c*   of the quartic polynoms                                       *
c*   see: https://github.com/jacobwilliams/polyroots-fortran       *
c*                                                                 *
c* author: Jiří Janek                                              *
c* date: January 2026                                              *
c*                                                                 *
c* compilation:                                                    *
c*   > gfortran -Wall -Wno-tabs -c -O3 wall_sph.for                *
c* the module must then be linked against polyroots_module.o and   *
c* with flags: -llapack -lblas -lm                                 *
c*                                                                 *
c*******************************************************************
      module wall_sph

        use polyroots_module

        implicit none

        real*8, parameter :: m_pi=3.141592653589793d0
        real*8, parameter :: m_inf=1.0d6

      contains

c*******************************************************************
c*                                                                 *
c* Roots of the quartic equation p_4^a (needed for summation)      *
c*     3*etab*rho^2*(rho+1)^2                                      *
c*       -2*[1+etab*rho*(4*rho^2 + 6*rho + 3)]t                    *
c*       +3*etab(2*rho^2 + 2*rho + 1)*t^2 - etab*t^4 = 0           *
c*                                                                 *
c* parameters:                                                     *
c*     etab ... bulk packing fraction                              *
c*     rho  ... radius of the sphere (spherical wall)              *
c*                                                                 *
c* output:                                                         *
c*     t(4) ... four complex roots of the equation                 *
c*                                                                 *
c*******************************************************************
        subroutine roots4a(etab, rho, t)
c
          real*8, intent(in) :: etab, rho
          complex*16, intent(out) :: t(4)
          real*8 :: coeffs(5)
          real*8 :: treal(4), timag(4)
          integer :: i

          coeffs(1) = 3.0d0*etab*rho**2*(rho+1.0d0)**2
          coeffs(2) = -2.0d0 * (1.0d0
     &          + etab*rho*(4.0d0*rho**2 + 6.0d0*rho + 3.0d0))
          coeffs(3) = 3.0d0*etab*(2.0d0*rho**2 + 2.0d0*rho + 1.0d0)
          coeffs(4) = 0.0d0
          coeffs(5) = -etab

          call dqtcrt(coeffs, treal, timag)
          do 101 i = 1, 4
            t(i) = cmplx(treal(i), timag(i), kind=8)
  101     continue

          return
        end

c*******************************************************************
c*                                                                 *
c* roots of the quartic equation p_4^b (needed for summation)      *
c*     2*rho + 2*t - 3*etab*(2*rho+1)^2*t^2                        *
c*       + 4*rho*etab*t^3 + etab*t^4 = 0                           *
c*                                                                 *
c* parameters:                                                     *
c*     etab ... bulk packing fraction                              *
c*     rho  ... radius of the sphere (spherical wall)              *
c*                                                                 *
c* output:                                                         *
c*     t(4) ... four complex roots of the equation                 *
c*                                                                 *
c*******************************************************************
        subroutine roots4b(etab, rho, t)
c
          real*8, intent(in) :: etab, rho
          complex*16, intent(out) :: t(4)
          real*8 :: coeffs(5)
          real*8 :: treal(4), timag(4)
          integer :: i

          coeffs(1) = 2.0d0 * rho
          coeffs(2) = 2.0d0
          coeffs(3) = -3.0d0*etab*(2.0d0*rho + 1.0d0)
          coeffs(4) = 4.0d0*etab*rho
          coeffs(5) = etab

          call dqtcrt(coeffs, treal, timag)
          do 101 i = 1, 4
            t(i) = cmplx(treal(i), timag(i), kind=8)
  101     continue

          return
        end

        real*8 function C0(r, rho, etab)
c====================================================================
c   Direct correlation function (part) c^(1)_0 in spherical geometry.
c
c   input:
c        r    ... distance from the sphere center
c        rho  ... radius of the sphere (spherical wall)
c        etab ... bulk packing fraction
c
c   output:
c        Direct correlation function c^(1)_0(r) for given rho and etab
c====================================================================

          real*8 :: r, rho, etab
          integer :: k

          real*8 :: logp1, nons
          complex*16 :: t(4), sumnum, sumden, logarg, sum

          call roots4b(etab, rho, t)


c====================================================================
c   region 1: r < rho
c====================================================================
          if (r .lt. rho) then
            c0 = 0.0d0
            return
c====================================================================
c   region 2: rho <= r < rho + 1/2
c====================================================================
          else if (r .ge. rho .and. r .lt. rho + 0.5d0) then

c----- prefactor part: (2r+2rho+1)(2r-2rho+1)/(16r) * [2 ln(..) - 3]
            logp1 =
     &      ((4d0*r*r + 8d0*r*rho - 12d0*rho*rho - 20d0*rho
     &        + 4d0*r - 11d0)
     &        * etab * (2d0*r - 2d0*rho + 1d0)**2 )
     &        / (16d0*(2d0*r + 1d0)) + 1d0

            nons = ((2d0*r + 2d0*rho + 1d0)*(2d0*r - 2d0*rho + 1d0))
     &           /(16d0*r)
     &           * (2d0*dlog(logp1) - 3d0)

c----- rho^2/(2r) * ln(2 rho / (2r+1))
            nons = nons 
     &           + (rho*rho)/(2d0*r) * dlog(2d0*rho/(2d0*r + 1d0))

c----- summation: - sum(...)
            sum = (0d0, 0d0)

            do 102 k=1, 4
              sumnum =
     &        ((4d0*rho*rho + 6d0*rho + 3d0)*etab)*t(k)**3
     &        - 3d0*t(k)**2
     &        - 6d0*t(k)*rho
     &        - 4d0*rho*rho

              sumden = 2d0*((2d0*t(k)**2
     &               + 6d0*(t(k)-1d0)*rho
     &               - 3d0)*etab*r*t(k)) + 2d0*r

              logarg = (2d0*r - 2d0*t(k) - 2d0*rho + 1d0)/(-2d0*t(k))

              sum = sum + sumnum/sumden * zlog(logarg)
  102       continue


            c0 = nons - dble(sum)
            return
c====================================================================
c   region 3: rho + 1/2 <= r < rho + 3/2
c====================================================================
          else if (r .ge. rho + 0.5d0 .and. r .lt. rho + 1.5d0) then

c----- term 1
            logp1 =
     &      ((4d0*r*r - 12d0*rho*rho + 8d0*r*rho
     &        - 4d0*r - 28d0*rho - 11d0)
     &        * etab*(2d0*r - 2d0*rho - 1d0)**2 )
     &        /(16d0*(2d0*r - 1d0)) + 1d0

            nons = -((2d0*r + 2d0*rho - 1d0)*(2d0*r - 2d0*rho - 1d0))
     &           /(8d0*r) * dlog(logp1)

c----- term 2
            nons =  nons + (3d0*(2d0*r + 2d0*rho + 1d0)*
     &           (2d0*r - 2d0*rho - 3d0))/(16d0*r)

c----- term 3
            nons =  nons + (rho*rho)/(2d0*r)
     &           * dlog((2d0*r - 1d0)/(2d0*rho + 2d0))

c----- term 4
            nons = nons + ((2d0*r + 2d0*rho + 1d0)*
     &             (2d0*r - 2d0*rho + 1d0))/(8d0*r) * dlog(1d0 - etab)

c----- summation
            sum = (0d0, 0d0)

            do 202 k=1, 4
              sumnum =
     &        ((4d0*rho*rho + 6d0*rho + 3d0)*etab)*t(k)**3
     &        - 3d0*t(k)**2
     &        - 6d0*t(k)*rho
     &        - 4d0*rho*rho

              sumden = 2d0*((2d0*t(k)**2
     &               + 6d0*(t(k)-1d0)*rho
     &               - 3d0)*etab*r*t(k)) + 2d0*r

              logarg =
     &        (2d0*r - 2d0*t(k) - 2d0*rho - 1d0)/(2d0 - 2d0*t(k))

              sum = sum + sumnum/sumden * zlog(logarg)
  202       continue

            C0 = nons + dble(sum)
            return
c====================================================================
c   region 4: r >= rho + 3/2
c====================================================================
          else if (r .ge. rho + 1.5d0) then
            c0 = dlog(1d0 - etab)
            return
c====================================================================
c   fallback (should not be reached)
c====================================================================
          else
            C0 = 0d0
            return
          endif
        end

        real*8 function C1(r, rho, etab)
c====================================================================
c   Direct correlation function (part) c^(1)_1 in spherical geometry.
c
c   input:
c        r    ... distance from the sphere center
c        rho  ... radius of the sphere (spherical wall)
c        etab ... bulk packing fraction
c
c   output: value of c_1(r) for given rho and etab
c
c====================================================================

          real*8 r, rho, etab
          complex*16 t(4)
          complex*16 num, den, st, logarg, sum
          integer*4 k
          sum = (0.0d0, 0.0d0)

          call roots4b(etab, rho, t)

c====================================================================
c   region 1: r < rho
c====================================================================
          if (r .lt. rho) then
            C1 = 0.0d0
            return
cc====================================================================
c   region 2: rho <= r < rho + 1/2
c====================================================================
          else if (r .ge. rho .and. r .lt. rho + 0.5d0) then
            sum = (0.0d0, 0.0d0)
            do k = 1, 4
              num = (t(k) + 2.0d0*rho + 1.0d0)*(t(k) + rho)*t(k)
              den = ((2.0d0*t(k)*t(k) + 6.0d0*(t(k) - 1.0d0)*rho
     &                  - 3.0d0)*etab*r*t(k) + r)
              st  = num / den
              logarg = (2.0d0*r - 2.0d0*t(k) - 2.0d0*rho + 1.0d0) /
     &           (-2.0d0*t(k))
              sum = sum + st * zlog(logarg)
            end do
            C1 = -1.5d0 * etab * dble(sum)
            return
c====================================================================
c   region 3: rho + 1/2 <= r < rho + 3/2
c====================================================================
          else if (r .ge. rho + 0.5d0 .and. r .lt. rho + 1.5d0) then
            sum = (0.0d0, 0.0d0)
            do k = 1, 4
              num = (t(k) + 2.0d0*rho + 1.0d0)*(t(k) + rho)*t(k)
              den = ((2.0d0*t(k)*t(k) + 6.0d0*(t(k) - 1.0d0)*rho
     &                  - 3.0d0)*etab*r*t(k) + r)
              st  = num / den
              logarg = (2.0d0*r - 2.0d0*t(k) - 2.0d0*rho - 1.0d0) /
     &                   (2.0d0 - 2.0d0*t(k))
              sum = sum + st * zlog(logarg)
            end do
            C1 = 3.0d0*(2.0d0*r + 2.0d0*rho + 3.0d0)*
     &           (2.0d0*r - 2.0d0*rho - 1.0d0)*etab /
     &           (8.0d0*(etab - 1.0d0)*r) + 1.5d0 * etab * dble(sum)
            return
c====================================================================
c   region 4: r >= rho + 3/2
c====================================================================
          else if (r .ge. rho + 1.5d0) then
            C1 = 3.0d0 * etab / (etab - 1.0d0)
            return
c====================================================================
c   fallback (should not be reached)
c====================================================================
            C1 = 0.0d0
            return
          else
            C1 = 0.0d0
            return
          end if

          return
        end

c========= c_2 calculation ===========================================

        complex*16 function Q2(t, rho, etab)
c=====================================================================
c Auxiliary function – numerator for the sum part in c^(1)_2 calculation
c
c input: t ... one of the roots of the quartic polynomial over which the sum is done
c      rho ... radius of the sphere (spherical wall)
c     etab ... bulk packing fraction
c
c output: value of Q_2(t, rho, etab)
c
c======================================================================

          complex*16 t
          real*8 rho, etab, rho2, rho3, rho4, rho5
          complex*16 teta1, teta2, teta3

          rho2 = rho*rho
          rho3 = rho2*rho
          rho4 = rho2*rho2
          rho5 = rho2*rho3

          ! teta1: etab^3 piece
          teta1 = 2.d0 * ( 8.d0*(4.d0*rho2 + 6.d0*rho + 3.d0)*t**2
     &             + (10.d0*rho2 + 15.d0*rho + 6.d0)*(2.d0*rho + 3.d0)*t
     &             - 5.2d1*rho4 - 1.56d2*rho3 - 1.77d2*rho2
     &               - 9.0d1*rho - 1.8d1)
     &             * (rho - t) * rho * etab**3

          ! teta2: etab^2 piece
          teta2 = (8.d0*(8.d0*rho3 + 12.d0*rho2 + 6.d0*rho - 1.d0)*t**3
     &        + 2.d0*(48.d0*rho3 + 72.d0*rho2 + 36.d0*rho - 3.d0)*t**2
     &        - 6.d0*(16.d0*rho5 + 40.d0*rho4 + 28.d0*rho3
     &          - 4.d0*rho2 - 10.d0*rho - 1.d0)*t
     &            + (8.d0*rho3 + 12.d0*rho2 + 6.d0*rho - 3.d0)*
     &                (4.d0*rho2 + 6.d0*rho + 3.d0)*rho ) * etab**2

          ! teta3: etab piece
          teta3 = (8.d0*t**3 + 1.2d1*t**2 
     &            - 6.d0*(2.d0*rho2+2.d0*rho-1.d0)*t
     &            + 3.d0*(4.d0*rho2 + 6.d0*rho + 3.d0)*rho 
     &            - 1.d0 ) * etab

          ! final result
          Q2 = teta1 + teta2 + teta3 + (1.d0, 0.d0)

          return
        end


        real*8 function S2(r, rho, etab)
c=====================================================================
c     Auxiliary function – numerator for the non-sum part in c^(1)_2 calculation
c         for rho < r < rho + 1/2
c
c     input: r ... position (from the sphere center)
c          rho ... radius of the sphere (spherical wall)
c         etab ... bulk packing fraction
c
c     output: value of S_2(r, rho, etab)
c
c======================================================================

          real*8 rho, r, etab
          real*8 teta3, teta2, teta1
          real*8 r2, r3, rho2, rho3, rho4, rho5, rho6

c     precompute powers of r and rho for clarity
          r2 = r*r
          r3 = r2*r
          rho2 = rho*rho
          rho3 = rho2*rho
          rho4 = rho2*rho2
          rho5 = rho3*rho2
          rho6 = rho3*rho3

c     term multiplying etab**3
          teta3 = 8.0d0 * (12.0d0*(2.0d0*rho2 + 3.0d0*rho + 2.0d0)
     &        *(rho + 1.0d0)*r2
     &        + 6.0d0*(2.0d0*rho2 + 3.0d0*rho + 2.0d0)
     &         *(4.0d0*rho2 + 8.0d0*rho + 5.0d0)*r
     &        - 3.0d0*(12.0d0*rho3 + 40.0d0*rho2 + 49.0d0*rho + 22.0d0)
     &         *(2.0d0*rho + 1.0d0)*rho + 6.0d0)
     &        *(2.0d0*r - 2.0d0*rho + 1.0d0)**2 * rho * etab**3

c     term multiplying etab**2
          teta2 = (24.0d0*(8.0d0*rho3 + 12.0d0*rho2 
     &          + 8.0d0*rho + 3.0d0)*r3
     &        + 12.0d0*(16.0d0*rho4 + 48.0d0*rho3 + 52.0d0*rho2
     &          + 30.0d0*rho + 17.0d0)*r2
     &        - 6.0d0*(32.0d0*rho5 + 80.0d0*rho4 + 248.0d0*rho3
     &          + 328.0d0*rho2 + 164.0d0*rho - 13.0d0)*r
     &        - 3.0d0*(64.0d0*rho6 + 256.0d0*rho5 + 352.0d0*rho4
     &          + 376.0d0*rho3 + 352.0d0*rho2 + 178.0d0*rho + 1.0d0))
     &        * (2.0d0*r - 2.0d0*rho + 1.0d0) * etab**2

c     term multiplying etab
          teta1 = 3.0d0*(8.0d0*r3 + 4.0d0*(2.0d0*rho + 3.0d0)*r2
     &        - (8.0d0*rho2 + 8.0d0*rho + 34.0d0)*r
     &        - 8.0d0*rho3 - 20.0d0*rho2 - 14.0d0*rho - 19.0d0)
     &        * (2.0d0*r - 2.0d0*rho + 1.0d0) * etab

          S2 = teta3 + teta2 + teta1

          return
        end

        real*8 function T2(r, rho, etab)
c=====================================================================
c     Auxiliary function – numerator for the non-sum part in c^(1)_2 calculation
c         for rho + 1/2 < r < rho + 3/2
c
c     input: r ... position (from the sphere center)
c          rho ... radius of the sphere (spherical wall)
c         etab ... bulk packing fraction
c
c     output: value of T_2(r, rho, etab)
c
c======================================================================

          real*8 rho,r,etab
          real*8 teta1, teta2, teta3, teta4
          real*8 r2,r3,r4,r5
          real*8 rho2,rho3,rho4,rho5,rho6,rho7,rho8

          r2 = r*r
          r3 = r2*r
          r4 = r3*r
          r5 = r4*r

          rho2 = rho*rho
          rho3 = rho2*rho
          rho4 = rho2*rho2
          rho5 = rho3*rho2
          rho6 = rho3*rho3
          rho7 = rho4*rho3
          rho8 = rho4*rho4

c--- etab^4 term -------------------------------------------------------
          teta4 = -(384d0*(4d0*rho2+6d0*rho+3d0)*r5*(r-1d0)
     &      -96d0*(56d0*rho3+112d0*rho2+88d0*rho+21d0)*(2d0*rho+1d0)*r4
     &      +192d0*(64d0*rho5+208d0*rho4+288d0*rho3
     &           +192d0*rho2+48d0*rho-3d0)*r3
     &      +24d0*(48d0*rho4+184d0*rho3+184d0*rho2
     &          +6d0*rho-39d0)*(2d0*rho+1d0)**2*r2
     &      -24d0*(128d0*rho5+528d0*rho4+728d0*rho3
     &          +280d0*rho2-146d0*rho-105d0)*(2d0*rho+1d0)**2*r
     &      +6d0*(48d0*rho4+136d0*rho3+80d0*rho2
     &         -46d0*rho-45d0)*(2d0*rho+3d0)*(2d0*rho+1d0)**3)
     &      *rho*etab**4

c--- etab^3 term -------------------------------------------------------
          teta3 = -(192d0*(16d0*rho3+24d0*rho2+12d0*rho+1d0)*r5*(r-1d0)
     &      -48d0*(448d0*rho5+1216d0*rho4+1472d0*rho3
     &          +932d0*rho2+280d0*rho+7d0)*r4
     &      +96d0*(256d0*rho6+832d0*rho5+1200d0*rho4
     &          +992d0*rho3+496d0*rho2+124d0*rho-1d0)*r3
     &      +12d0*(768d0*rho7+6016d0*rho6+16704d0*rho5
     &          +22896d0*rho4+16784d0*rho3
     &          +6144d0*rho2+828d0*rho-13d0)*r2
     &      -12d0*(2048d0*rho8+13568d0*rho7+37888d0*rho6
     &          +58048d0*rho5+52880d0*rho4+28176d0*rho3
     &          +7688d0*rho2+636d0*rho-35d0)*r
     &      +3d0*(1536d0*rho8+8960d0*rho7+21888d0*rho6
     &         +29408d0*rho5+23696d0*rho4+11184d0*rho3
     &         +2600d0*rho2+90d0*rho-15d0)*(2d0*rho+3d0))*etab**3

c--- etab^2 term -------------------------------------------------------
          teta2 = -(384d0*(r-1d0)*r5
     &      -96d0*(8d0*rho3+40d0*rho2+40d0*rho+21d0)*r4
     &      +192d0*(56d0*rho3+88d0*rho2+52d0*rho+9d0)*r3
     &      +24d0*(64d0*rho5+208d0*rho4+688d0*rho3
     &          +976d0*rho2+544d0*rho+71d0)*r2
     &      -24d0*(448d0*rho5+1840d0*rho4+3184d0*rho3
     &          +2640d0*rho2+968d0*rho+65d0)*r
     &      -6d0*(64d0*rho6+32d0*rho5-592d0*rho4
     &         -1296d0*rho3-1108d0*rho2-398d0*rho-15d0)
     &         *(2d0*rho+3d0))*etab**2

c--- etab term ---------------------------------------------------------
          teta1 = (96d0*(r-10d0)*r3
     &      -48d0*(4d0*rho2+4d0*rho+15d0)*r2
     &      +48d0*(20d0*rho2+44d0*rho+37d0)*r
     &      +6d0*(4d0*rho2-4d0*rho-11d0)*(2d0*rho+3d0)**2)*etab

          T2 = teta4 + teta3 + teta2 + teta1
          return
        end

        real*8 function C2(r, rho, etab)
c=======================================================================
c  Direct correlation function (part) c^(1)_2 in spherical geometry
c
c  inputs: r    ... distance from the sphere center
c          rho  ... sphere (spherical wall) radius
c          etab ... bulk packing fraction
c
c  output: value of c_2(r) for given rho and etab
c
c=======================================================================


          real*8 r, rho, etab
          complex*16 t(4), sum, lga, cden
          integer k
          real*8 den1,den2,logterm

          call roots4a(etab, rho, t)

          den1 = 2.0d0 * (4.0d0*rho**2 + 6.0d0*rho + 3.0d0)*etab*rho
          den1 = (den1+1.0d0)*(etab-1.0d0)

c====================================================================
c   region 1: r < rho
c====================================================================
          if (r.lt.rho) then
            C2 = 0.0d0
c====================================================================
c   region 2: rho <= r < rho + 1/2
c====================================================================
          else if ((r .ge. rho) .and. (r .lt. rho + 0.5d0)) then
            den2 = (4.0d0*(r**2 + r + 2.0d0*r*rho - 3.0d0*rho**2
     &        - 5.0d0*rho) - 1.1d1)*etab
     &        *((2.0d0*r - 2.0d0*rho + 1.0d0)**2)
     &        + 1.6d1*(2.0d0*r + 1.0d0)
            logterm = log((2.0d0*r + 1.0d0)/(2.0d0*rho))/(8.0d0*r)
            sum = (0.0d0, 0.0d0)
            do k = 1,4
              cden = (2.0d0*(t(k)+2.0d0*rho)*(t(k)-rho)-6.0d0*rho-3.0d0)
     &             *etab*(t(k)-rho)+1.0d0
              lga = (2.0d0*r-2.0d0*t(k)+1.0d0)/(2.0d0*rho-2.0d0*t(k))
              sum = sum 
     &            + Q2(t(k), rho, etab)/(8.0d0*r*cden*den1)*zlog(lga)
            end do
            C2 = S2(r, rho, etab)/(8.0*r*den2*den1) 
     &         + logterm + dble(sum)
c====================================================================
c   region 3: rho + 1/2 <= r < rho + 3/2
c====================================================================
          else if ((r .ge. rho+0.5d0).and.(r .lt. rho+1.5d0)) then
            den2 = (4.0d0*(r**2 - r + 2.0d0*r*rho - 3.0d0*rho**2
     &        - 7.0d0*rho) - 1.1d1)*etab
     &        *((2.0d0*r - 2.0d0*rho - 1.0d0)**2)
     &        + 1.6d1*(2.0d0*r - 1.0d0)
            logterm = -log((2.0d0*r - 1.0d0)/(2.0d0*rho + 2.0d0))
     &              /(8.0d0*r)
            sum = (0.0d0, 0.0d0)
            do k = 1,4
              cden = (2.0d0*(t(k)+2.0d0*rho)*(t(k)-rho)-6.0d0*rho-3.0d0)
     &             *etab*(t(k)-rho)+1.0d0
              lga = (2.0d0*r-2.0d0*t(k)-1.0d0)/2.0d0/(rho-t(k)+1.0d0)
              sum = sum 
     &            + Q2(t(k), rho, etab)/(8.0*r*cden*den1)*zlog(lga)
            end do
            C2 = t2(r, rho, etab)/(1.6d1*r*den2*den1*(etab-1.0d0))
     &           + logterm - dble(sum)
c====================================================================
c   region 4: r >= rho + 3/2
c====================================================================
          else if (r .ge. rho+1.5d0) then
            C2 = -3.0d0*etab*(2.0d0+etab)/(2.0d0*(1.0d0-etab)**2)
            return
c====================================================================
c   fallback (should not be reached)
c====================================================================
          else
            C2 = 0.0d0
            return
          end if

          return
        end

c===== end of c_2 calculation ==========================================

c===== c_3 calculation =================================================

        complex*16 function Q3(t,r,rho,etab)
c=====================================================================
c     Auxiliary function – numerator for the sum part in c^(1)_3 calculation
c
c     input: t ... one of the roots of the quartic polynomial over which the sum is done
c            r ... distance from the sphere center
c          rho ... radius of the sphere (spherical wall)
c         etab ... bulk packing fraction
c
c     output: complex*16 value of Q_3(t, rho, etab)
c
c======================================================================

          complex*16 t
          real*8 r,rho,etab
          complex*16 termt3,teta5,teta4,teta3,teta2,teta1,teta0
          real*8 rho2, rho3, rho4, rho5, rho6
          real*8 rho7, rho8, rho9, rho10, rho11

          rho2 = rho*rho
          rho3 = rho2*rho
          rho4 = rho2*rho2
          rho5 = rho3*rho2
          rho6 = rho3*rho3
          rho7 = rho3*rho4
          rho8 = rho4*rho4
          rho9 = rho4*rho5
          rho10 = rho5*rho5
          rho11 = rho5*rho6

          termt3 =2.d0*(
     &        4.d0*(4.d0*rho2+6.d0*rho+3.d0)**2*rho2*etab**4
     &       -4.d0*(4.d0*rho2+6.d0*rho+3.d0)
     &            *(8.d0*rho3+12.d0*rho2+6.d0*rho-1.d0)
     &            *rho*etab**3
     &       +(64.d0*rho6+192.d0*rho5+240.d0*rho4
     &         +112.d0*rho3-12.d0*rho2-24.d0*rho+1.d0)
     &            *etab**2
     &       +2.d0*(8.d0*rho3+12.d0*rho2+6.d0*rho-1.d0)
     &            *etab
     &       +1.d0
     &      )*(
     &        4.d0*(2.d0*rho+1.d0)**2*r**2
     &       -96.d0*(rho+1.d0)*rho*r
     &       -48.d0*rho4-96.d0*rho3
     &       -16.d0*rho2+32.d0*rho-1.d0
     &      )*etab*t**3

          teta5 = (
     &        4.d0*(4.d0*(328.d0*rho4+948.d0*rho3
     &        +1122.d0*rho2+639.d0*rho+153.d0)*r**2
     &        +8.d0*(88.d0*rho5-52.d0*rho4
     &        -590.d0*rho3-879.d0*rho2
     &        -546.d0*rho-135.d0)*r
     &        +(3072.d0*rho7+22560.d0*rho6
     &        +61584.d0*rho5+94112.d0*rho4
     &        +90600.d0*rho3+56442.d0*rho2
     &        +21285.d0*rho+3843.d0))
     &        *(rho+1.d0)*rho3*t**2
     &       -8.d0*(12.d0*(64.d0*rho8+344.d0*rho7
     &        +820.d0*rho6+1122.d0*rho5
     &        +979.d0*rho4+576.d0*rho3
     &        +234.d0*rho2+63.d0*rho+9.d0)*r**2
     &        -8.d0*(104.d0*rho6+596.d0*rho5
     &        +1198.d0*rho4+1143.d0*rho3
     &        +522.d0*rho2+81.d0*rho-9.d0)
     &        *(rho+1.d0)*rho*r
     &        +768.d0*rho10+11296.d0*rho9
     &        +50032.d0*rho8+115248.d0*rho7
     &        +161352.d0*rho6+144018.d0*rho5
     &        +80763.d0*rho4+25920.d0*rho3
     &        +3402.d0*rho2-189.d0*rho-27.d0)
     &        *rho2*t
     &       +4.d0*(4.d0*(256.d0*rho8+1480.d0*rho7
     &        +3964.d0*rho6+6390.d0*rho5
     &        +6849.d0*rho4+5040.d0*rho3
     &        +2493.d0*rho2+756.d0*rho+108.d0)*r**2
     &        -24.d0*(56.d0*rho4+140.d0*rho3
     &        +106.d0*rho2+21.d0*rho-6.d0)
     &        *(rho+1.d0)**3*rho*r
     &        +5664.d0*rho9+34736.d0*rho8
     &        +93488.d0*rho7+142856.d0*rho6
     &        +133218.d0*rho5+74799.d0*rho4
     &        +22176.d0*rho3+1503.d0*rho2
     &        -756.d0*rho-108.d0)*rho3
     &     )*etab**5

          teta4 =
     &     -(
     &       -(
     &         4.d0*(576.d0*rho6-256.d0*rho5-3504.d0*rho4
     &         -5280.d0*rho3-3012.d0*rho2-396.d0*rho
     &         +207.d0)*r**2
     &        -8.d0*(576.d0*rho7+3520.d0*rho6
     &         +5360.d0*rho5+1600.d0*rho4
     &         -3308.d0*rho3-3052.d0*rho2
     &         -611.d0*rho+180.d0)*r
     &        -(22272.d0*rho8+167424.d0*rho7
     &         +462336.d0*rho6+710528.d0*rho5
     &         +676320.d0*rho4+405264.d0*rho3
     &         +135072.d0*rho2+12564.d0*rho
     &         -5121.d0)
     &       )*(rho+1.d0)*rho2*t**2
     &       -2.d0*(12.d0*(704.d0*rho9+3904.d0*rho8
     &         +9584.d0*rho7+13424.d0*rho6
     &         +11724.d0*rho5+6548.d0*rho4
     &         +2289.d0*rho3+437.d0*rho2
     &         +12.d0*rho-12.d0)*r**2
     &        -8.d0*(576.d0*rho8+4160.d0*rho7
     &         +13040.d0*rho6+20992.d0*rho5
     &         +18348.d0*rho4+8088.d0*rho3
     &         +1173.d0*rho2-168.d0*rho
     &         +12.d0)*(rho+1.d0)*rho*r
     &        +8448.d0*rho11+107776.d0*rho10
     &         +457216.d0*rho9+1025664.d0*rho8
     &         +1388864.d0*rho7+1163552.d0*rho6
     &         +562344.d0*rho5+108240.d0*rho4
     &         -23427.d0*rho3-11823.d0*rho2
     &         -36.d0*rho+36.d0
     &       )*rho*t
     &       +(4.d0*(3776.d0*rho9+21056.d0*rho8
     &         +53168.d0*rho7+78832.d0*rho6
     &         +75172.d0*rho5+46864.d0*rho4
     &         +17865.d0*rho3+2859.d0*rho2
     &         -585.d0*rho-279.d0)*r**2
     &        -24.d0*(576.d0*rho6+2240.d0*rho5
     &         +3440.d0*rho4+2368.d0*rho3
     &         +596.d0*rho2-52.d0*rho
     &         +1.d0)*(rho+1.d0)**3*rho*r
     &        +6912.d0*rho11+92928.d0*rho10
     &         +419584.d0*rho9+985216.d0*rho8
     &         +1379680.d0*rho7+1191248.d0*rho6
     &         +603104.d0*rho5+139292.d0*rho4
     &         -8973.d0*rho3-8007.d0*rho2
     &         +585.d0*rho+279.d0
     &       )*rho2)*etab**4

          teta3 =((4.d0*(1536.d0*rho6+4608.d0*rho5
     &         +5832.d0*rho4+3132.d0*rho3
     &         +306.d0*rho2-279.d0*rho
     &         +17.d0)*r**2
     &        -8.d0*(1920.d0*rho6+5976.d0*rho5
     &         +7892.d0*rho4+4546.d0*rho3
     &         +647.d0*rho2-324.d0*rho
     &         +15.d0)*r
     &        +12288.d0*rho9+92160.d0*rho8
     &         +258048.d0*rho7+395040.d0*rho6
     &         +353520.d0*rho5+177408.d0*rho4
     &         +26712.d0*rho3-18942.d0*rho2
     &         -9573.d0*rho+427.d0
     &       )*(rho+1.d0)*rho*t**2
     &       -((6144.d0*rho10+30720.d0*rho9
     &         +66048.d0*rho8+73536.d0*rho7
     &         +37728.d0*rho6-3504.d0*rho5
     &         -15672.d0*rho4-9000.d0*rho3
     &         -2712.d0*rho2-504.d0*rho
     &         +24.d0)*r**2
     &        +16.d0*(1224.d0*rho6+4260.d0*rho5
     &         +6258.d0*rho4+4335.d0*rho3
     &         +1199.d0*rho2-25.d0*rho
     &         +1.d0)*(rho+1.d0)*rho*r
     &        +6144.d0*rho**12+86016.d0*rho11
     &         +365568.d0*rho10+795456.d0*rho9
     &         +968416.d0*rho8+573184.d0*rho7
     &         -45744.d0*rho6-307484.d0*rho5
     &         -182498.d0*rho4-30390.d0*rho3
     &         +4950.d0*rho2+126.d0*rho
     &         -6.d0
     &       )*t
     &       +(4.d0*(1024.d0*rho10+5632.d0*rho9
     &         +14848.d0*rho8+22952.d0*rho7
     &         +21428.d0*rho6+10262.d0*rho5
     &         -577.d0*rho4-4018.d0*rho3
     &         -2233.d0*rho2-378.d0*rho+60.d0)*r**2
     &        +24.d0*(216.d0*rho3+388.d0*rho2
     &         +226.d0*rho-1.d0)*(rho+1.d0)**3*rho2*r
     &        +18432.d0*rho11+109568.d0*rho10
     &         +283808.d0*rho9+400592.d0*rho8
     &         +304432.d0*rho7+79480.d0*rho6
     &         -52286.d0*rho5-43487.d0*rho4
     &         -5918.d0*rho3+2653.d0*rho2+378.d0*rho-60.d0
     &       )*rho
     &     )*etab**3

          teta2 = ((4.d0*(384.d0*rho3+576.d0*rho2+288.d0*rho-23.d0)*r**2
     &        -24.d0*(160.d0*rho3+246.d0*rho2+126.d0*rho-9.d0)*r
     &        +6144.d0*rho6+27648.d0*rho5+46080.d0*rho4
     &         +41088.d0*rho3+19212.d0*rho2+3948.d0*rho-853.d0
     &       )*(rho+1.d0)*rho*t**2
     &       +(
     &         (-1536.d0*rho7-5376.d0*rho6-7296.d0*rho5
     &         -4272.d0*rho4-864.d0*rho3-24.d0*rho2
     &         -72.d0*rho+48.d0)*r**2
     &        -16.d0*(288.d0*rho3+481.d0*rho2+265.d0*rho-1.d0)
     &         *(rho+1.d0)*rho*r
     &        -1536.d0*rho9-31488.d0*rho8-110592.d0*rho7
     &         -176112.d0*rho6-138992.d0*rho5-38948.d0*rho4
     &         +13896.d0*rho3+9038.d0*rho2-262.d0*rho-12.d0
     &       )*t
     &
     &       +(
     &         (2048.d0*rho8+8192.d0*rho7+15872.d0*rho6
     &         +17792.d0*rho5+11804.d0*rho4+3608.d0*rho3
     &         -308.d0*rho2-464.d0*rho+16.d0)*r**2
     &        +432.d0*(rho+1.d0)**3*rho3*r
     &        + (1152.d0*rho7+3904.d0*rho6+4800.d0*rho5
     &         +1929.d0*rho4-614.d0*rho3-435.d0*rho2
     &         +132.d0*rho-4.d0)*(2.d0*rho+1.d0)**2
     &       )
     &     )*etab**2

          teta1 =(48.d0*(2.d0*r**2-5.d0*r+20.d0*rho3+42.d0*rho2
     &            +27.d0*rho+9.d0)*(rho+1.d0)*rho*t**2
     &       -2.d0*(12.d0*(4.d0*rho4+8.d0*rho3+4.d0*rho2+1.d0)*r**2
     &          +144.d0*(rho+1.d0)*rho*r
     &          +48.d0*rho6+2064.d0*rho5+4932.d0*rho4
     &          +4344.d0*rho3+1120.d0*rho2-308.d0*rho-3.d0)*t
     &       +(16.d0*(20.d0*rho3+30.d0*rho2+15.d0*rho-2.d0)*
     &            (rho2+rho+1.d0)*r**2
     &          +288.d0*rho6+784.d0*rho5+664.d0*rho4
     &          +28.d0*rho3-172.d0*rho2-52.d0*rho+8.d0)
     &     )*etab


          teta0 = 48.d0*(rho+1.d0)*rho*t*(t-4.d0)
     &     +4.d0*(rho2+rho+1.d0)*(2.d0*r+1.d0)*(2.d0*r-1.d0)

          Q3 = termt3 + teta5 + teta4 + teta3 + teta2 + teta1 + teta0

          return
        end

        real*8 function S3(r,rho,etab)
c=====================================================================
c     Auxiliary function – numerator for the non-sum part in c^(1)_3 calculation
c         for rho < r < rho + 1/2
c
c     input: r ... position (from the sphere center)
c          rho ... radius of the sphere (spherical wall)
c         etab ... bulk packing fraction
c
c     output: value of S_3(r, rho, etab)
c
c======================================================================

          real*8 r,rho,etab
          real*8 pre,teta6,teta5,teta4,teta3,teta2,teta1,teta0
          real*8 rho2, rho3, rho4, rho5, rho6
          real*8 rho7, rho8, rho9, rho10, rho11

          rho2 = rho*rho
          rho3 = rho2*rho
          rho4 = rho2*rho2
          rho5 = rho3*rho2
          rho6 = rho3*rho3
          rho7 = rho3*rho4
          rho8 = rho4*rho4
          rho9 = rho4*rho5
          rho10 = rho5*rho5
          rho11 = rho5*rho6

          pre = (2.d0*r-2.d0*rho+1.d0)

          teta6 = -16.d0*(4.d0*r*(r+2.d0*rho+1.d0)
     &       -12.d0*rho**2-20.d0*rho-11.d0)
     &     *(4.d0*rho**2+6.d0*rho+3.d0)**2
     &     *(2.d0*r-2.d0*rho+1.d0)**2
     &     *(2.d0*r-2.d0*rho-1.d0)
     &     *(2.d0*rho+1.d0)
     &     *rho**2*etab**6

          teta5 = 4.d0*(
     &      32.d0*(320.d0*rho6+976.d0*rho5+1200.d0*rho4
     &            +636.d0*rho3+72.d0*rho**2-54.d0*rho-9.d0)
     &            *r**5
     &     -16.d0*(640.d0*rho7+1376.d0*rho6-1136.d0*rho5
     &            -6440.d0*rho4-8236.d0*rho3-4728.d0*rho**2
     &            -1056.d0*rho+27.d0)*r**4
     &     -16.d0*(3840.d0*rho8+16832.d0*rho7+32384.d0*rho6
     &            +30592.d0*rho5+9168.d0*rho4-8748.d0*rho3
     &            -9072.d0*rho**2-2718.d0*rho-45.d0)*r**3
     &     +8.d0*(17920.d0*rho9+78464.d0*rho8+128832.d0*rho7
     &           +66112.d0*rho6-76896.d0*rho5-146856.d0*rho4
     &           -101124.d0*rho3-33336.d0*rho**2-4176.d0*rho
     &           +63.d0)*r**2
     &     -2.d0*(56320.d0*rho10+257792.d0*rho9+385792.d0*rho8
     &           +46016.d0*rho7-546624.d0*rho6-686224.d0*rho5
     &           -298304.d0*rho4+48980.d0*rho3+91896.d0*rho**2
     &           +26214.d0*rho+81.d0)*r
     &     +(30720.d0*rho11+147968.d0*rho10+202496.d0*rho9
     &      -70784.d0*rho8-448064.d0*rho7-374048.d0*rho6
     &      +47440.d0*rho5+230808.d0*rho4+99516.d0*rho3
     &      -12024.d0*rho**2-13536.d0*rho-99.d0)
     &     )*rho*etab**5

          teta4 = (
     &      32.d0*(384.d0*rho7+3456.d0*rho6+8832.d0*rho5
     &            +10992.d0*rho4+6864.d0*rho3+1812.d0*rho**2
     &            -28.d0*rho-5.d0)*r**5
     &     -16.d0*(768.d0*rho8+4992.d0*rho7+20864.d0*rho6
     &            +41888.d0*rho5+43600.d0*rho4+21352.d0*rho3
     &            +2220.d0*rho**2-1482.d0*rho+15.d0)*r**4
     &     -16.d0*(4608.d0*rho9+47616.d0*rho8+161664.d0*rho7
     &            +317248.d0*rho6+393792.d0*rho5+310176.d0*rho4
     &            +138336.d0*rho3+23792.d0*rho**2-3220.d0*rho
     &            -25.d0)*r**3
     &     +8.d0*(21504.d0*rho10+207360.d0*rho9+862464.d0*rho8
     &           +1933824.d0*rho7+2565824.d0*rho6+2048256.d0*rho5
     &           +922336.d0*rho4+177152.d0*rho3-10920.d0*rho**2
     &           -5474.d0*rho+35.d0)*r**2
     &     -2.d0*(67584.d0*rho11+632832.d0*rho10+3058688.d0*rho9
     &           +7361792.d0*rho8+9771392.d0*rho7+6720960.d0*rho6
     &           +1059264.d0*rho5-1829952.d0*rho4-1292528.d0*rho3
     &           -243324.d0*rho**2+33516.d0*rho+45.d0)*r
     &     +(36864.d0*rho**12+337920.d0*rho11+1855488.d0*rho10
     &      +4647424.d0*rho9+5665792.d0*rho8+2344192.d0*rho7
     &      -1824448.d0*rho6-2218048.d0*rho5-203296.d0*rho4
     &      +609032.d0*rho3+221516.d0*rho**2-17638.d0*rho-55.d0)
     &     )*etab**4

          teta3 = -2.d0*(
     &      32.d0*(192.d0*rho6+576.d0*rho5+672.d0*rho4
     &            +36.d0*rho3-444.d0*rho**2-281.d0*rho-1.d0)*r**5
     &     +16.d0*(384.d0*rho7-960.d0*rho6-4800.d0*rho5
     &            -6648.d0*rho4-2852.d0*rho3+834.d0*rho**2
     &            +885.d0*rho-67.d0)*r**4
     &     -16.d0*(2304.d0*rho8+8448.d0*rho7+19008.d0*rho6
     &            +24432.d0*rho5+15792.d0*rho4-2544.d0*rho3
     &            -9836.d0*rho**2-4817.d0*rho+119.d0)*r**3
     &     +8.d0*(1536.d0*rho9+39168.d0*rho8+138624.d0*rho7
     &           +203232.d0*rho6+112560.d0*rho5-40720.d0*rho4
     &           -83096.d0*rho3-34818.d0*rho**2-2509.d0*rho
     &           +243.d0)*r**2
     &     +2.d0*(15360.d0*rho10-125952.d0*rho9-577536.d0*rho8
     &           -904896.d0*rho7-322816.d0*rho6+722512.d0*rho5
     &           +1019792.d0*rho4+459692.d0*rho3+2836.d0*rho**2
     &           -49481.d0*rho+1359.d0)*r
     &     -(18432.d0*rho11-70656.d0*rho10-405504.d0*rho9
     &      -649344.d0*rho8-275520.d0*rho7+141664.d0*rho6
     &      -195376.d0*rho5-712936.d0*rho4-530084.d0*rho3
     &      -94034.d0*rho**2+35995.d0*rho-701.d0)
     &     )*etab**3

          teta2 = -(
     &      96.d0*(32.d0*rho3+48.d0*rho**2+22.d0*rho-15.d0)
     &            *r**5
     &     +16.d0*(192.d0*rho4-768.d0*rho3-1428.d0*rho**2
     &            -768.d0*rho+149.d0)*r**4
     &     -16.d0*(1152.d0*rho5+2496.d0*rho4+5112.d0*rho3
     &            +4692.d0*rho**2+1722.d0*rho-785.d0)*r**3
     &     +8.d0*(3840.d0*rho6+27648.d0*rho5+53616.d0*rho4
     &           +39264.d0*rho3+4248.d0*rho**2-6312.d0*rho
     &           -577.d0)*r**2
     &     +2.d0*(19968.d0*rho7-123648.d0*rho6-397920.d0*rho5
     &           -439056.d0*rho4-100112.d0*rho3+118600.d0*rho**2
     &           +72082.d0*rho-8117.d0)*r
     &     -(9216.d0*rho8-61440.d0*rho7-82368.d0*rho6
     &      +136704.d0*rho5+360912.d0*rho4+171200.d0*rho3
     &      -84340.d0*rho**2-82096.d0*rho+5587.d0)
     &     )*etab**2

          teta1 = -2.d0*(
     &       96.d0*r**5
     &     + 48.d0*(2.d0*rho-11.d0)*r**4
     &     - 48.d0*(12.d0*rho**2+8.d0*rho+33.d0)*r**3
     &     + 24.d0*(136.d0*rho3+372.d0*rho**2+262.d0*rho-29.d0)*r**2
     &     + 2.d0*(1776.d0*rho4-11136.d0*rho3
     &             -16968.d0*rho**2-7824.d0*rho+3043.d0)*r
     &     - (288.d0*rho5-3504.d0*rho4+7632.d0*rho3
     &        +16968.d0*rho**2+9450.d0*rho-3055.d0)
     &     )*etab


          teta0  = -96.d0*(2.d0*r+2.d0*rho-15.d0)*(2.d0*r+1.d0)

          S3 = pre*(teta6+teta5+teta4+teta3+teta2+teta1+teta0)

          return
        end

c======================================================================

        real*8 function T3(r,rho,etab)
c=====================================================================
c     Auxiliary function – numerator for the non-sum part in c^(1)_3 calculation
c         for rho + 1/2 < r < rho + 3/2
c
c     input: r ... position (from the sphere center)
c          rho ... radius of the sphere (spherical wall)
c         etab ... bulk packing fraction
c
c     output: value of T_3(r, rho, etab)
c
c======================================================================

          real*8 r,rho,etab
          real*8 teta6,teta5,teta4,teta3,teta2,teta1,teta0
          real*8 rho2, rho3, rho4, rho5, rho6
          real*8 rho7, rho8, rho9, rho10, rho11

          rho2 = rho*rho
          rho3 = rho2*rho
          rho4 = rho2*rho2
          rho5 = rho3*rho2
          rho6 = rho3*rho3
          rho7 = rho3*rho4
          rho8 = rho4*rho4
          rho9 = rho4*rho5
          rho10 = rho5*rho5
          rho11 = rho5*rho6

          teta6 = 16.d0*(
     &      64.d0*(2.d0-r)*(4.d0*rho**2+6.d0*rho+3.d0)**2*r**7
     &     +16.d0*(768.d0*rho6+3200.d0*rho5+6256.d0*rho4
     &             +7248.d0*rho3+5232.d0*rho**2
     &             +2196.d0*rho+423.d0)*r**6
     &     -32.d0*(512.d0*rho7+3008.d0*rho6+7760.d0*rho5
     &             +12536.d0*rho4+13924.d0*rho3
     &             +10446.d0*rho**2+4710.d0*rho
     &             +981.d0)*r**5
     &     -4.d0*(7680.d0*rho8+35840.d0*rho7+83264.d0*rho6
     &            +117888.d0*rho5+88080.d0*rho4
     &            +1632.d0*rho3-55692.d0*rho**2
     &            -42768.d0*rho-11241.d0)*r**4
     &     +8.d0*(12288.d0*rho9+76544.d0*rho8
     &            +226688.d0*rho7+434176.d0*rho6
     &            +587328.d0*rho5+559456.d0*rho4
     &            +356376.d0*rho3+137616.d0*rho**2
     &            +25848.d0*rho+855.d0)*r**3
     &     -(102400.d0*rho10+739328.d0*rho9
     &       +2481920.d0*rho8+5415680.d0*rho7
     &       +8795008.d0*rho6+10944192.d0*rho5
     &       +10115488.d0*rho4+6578384.d0*rho3
     &       +2805000.d0*rho**2+696876.d0*rho
     &       +75843.d0)*r**2
     &     +2.d0*(6144.d0*rho9+43776.d0*rho8
     &            +140224.d0*rho7+295328.d0*rho6
     &            +494784.d0*rho5+667552.d0*rho4
     &            +660492.d0*rho3+432246.d0*rho**2
     &            +165654.d0*rho+28341.d0)
     &            *(2.d0*rho+1.d0)**2*r
     &     -(1152.d0*rho9+7488.d0*rho8+21376.d0*rho7
     &       +42816.d0*rho6+77096.d0*rho5
     &       +113852.d0*rho4+116712.d0*rho3
     &       +74748.d0*rho**2+26793.d0*rho
     &       +4059.d0)*(2.d0*rho+3.d0)*(2.d0*rho+1.d0)**2
     &     )*rho**2*etab**6
c======================================================================
          teta5 = 2.d0*(
     &     512.d0*(2.d0-r)*(4.d0*rho**2+6.d0*rho+3.d0)
     &            *(4.d0*rho**2+2.d0*rho+1.d0)*(rho+1.d0)*r**7
     &     +64.d0*(1536.d0*rho7+7232.d0*rho6
     &            +12832.d0*rho5+10704.d0*rho4
     &            +3552.d0*rho3+72.d0*rho**2
     &            +204.d0*rho+285.d0)*r**6
     &     -64.d0*(2048.d0*rho8+13696.d0*rho7
     &            +36352.d0*rho6+40384.d0*rho5
     &            +7456.d0*rho4-24504.d0*rho3
     &            -20208.d0*rho**2-3514.d0*rho
     &            +1317.d0)*r**5
     &     -16.d0*(15360.d0*rho9+88320.d0*rho8
     &            +167552.d0*rho7+141056.d0*rho6
     &            +149152.d0*rho5+336624.d0*rho4
     &            +443824.d0*rho3+266444.d0*rho**2
     &            +53604.d0*rho-7491.d0)*r**4
     &     +32.d0*(24576.d0*rho10+186368.d0*rho9
     &            +548864.d0*rho8+743680.d0*rho7
     &            +332288.d0*rho6-312928.d0*rho5
     &            -459552.d0*rho4-172392.d0*rho3
     &            +21564.d0*rho**2+22860.d0*rho
     &            +603.d0)*r**3
     &     -4.d0*(204800.d0*rho11+1811456.d0*rho10
     &           +6457856.d0*rho9+10749696.d0*rho8
     &           +4965632.d0*rho7-12269760.d0*rho6
     &           -25457312.d0*rho5-22228640.d0*rho4
     &           -10114944.d0*rho3-2004320.d0*rho**2
     &           +48548.d0*rho+50661.d0)*r**2
     &     +4.d0*(98304.d0*rho**12+985088.d0*rho11
     &           +4075520.d0*rho10+8044544.d0*rho9
     &           +4541952.d0*rho8-13309440.d0*rho7
     &           -35129088.d0*rho6-40819168.d0*rho5
     &           -27170736.d0*rho4-10160248.d0*rho3
     &           -1636536.d0*rho**2+81714.d0*rho
     &           +37731.d0)*r
     &     -(36864.d0*rho**12+356352.d0*rho11
     &       +1398784.d0*rho10+2332672.d0*rho9
     &       -331264.d0*rho8-8488960.d0*rho7
     &       -16720576.d0*rho6-17258112.d0*rho5
     &       -10497696.d0*rho4-3555024.d0*rho3
     &       -473100.d0*rho**2+46152.d0*rho
     &       +10791.d0)*(2.d0*rho+3.d0)
     &     )*rho*etab**5

          teta4 = 4.d0*(
     &      64.d0*(2.d0-r)*(2.d0*rho+1.d0)**6*r**7
     &     +16.d0*(3072.d0*rho8+14208.d0*rho7+32640.d0*rho6
     &            +45600.d0*rho5+40128.d0*rho4+19944.d0*rho3
     &            +3936.d0*rho**2-416.d0*rho+47.d0)*r**6
     &     -32.d0*(2048.d0*rho9+13440.d0*rho8+40896.d0*rho7
     &            +83232.d0*rho6+114352.d0*rho5+100080.d0*rho4
     &            +47040.d0*rho3+6522.d0*rho**2
     &            -2315.d0*rho+109.d0)*r**5
     &     -4.d0*(30720.d0*rho10+171520.d0*rho9
     &            +489984.d0*rho8+782336.d0*rho7
     &            +523264.d0*rho6-278016.d0*rho5
     &            -750304.d0*rho4-448248.d0*rho3
     &            -34052.d0*rho**2+44872.d0*rho
     &            -1249.d0)*r**4
     &     +8.d0*(49152.d0*rho11+362496.d0*rho10
     &            +1310208.d0*rho9+3083520.d0*rho8
     &            +4941056.d0*rho7+5278912.d0*rho6
     &            +3522528.d0*rho5+1257904.d0*rho4
     &            +124792.d0*rho3-26752.d0*rho**2
     &            +6858.d0*rho+95.d0)*r**3
     &     -(409600.d0*rho**12+3520512.d0*rho11
     &       +14518272.d0*rho10+39896576.d0*rho9
     &       +79227904.d0*rho8+112583424.d0*rho7
     &       +109728640.d0*rho6+68306784.d0*rho5
     &       +23047216.d0*rho4+1666408.d0*rho3
     &       -1144440.d0*rho**2-166128.d0*rho
     &       +8427.d0)*r**2
     &     +2.d0*(98304.d0*rho**13+956416.d0*rho**12
     &            +4396032.d0*rho11+13556224.d0*rho10
     &            +31287040.d0*rho9+53866368.d0*rho8
     &            +66479040.d0*rho7+55691456.d0*rho6
     &            +28814144.d0*rho5+6900176.d0*rho4
     &            -715544.d0*rho3-686514.d0*rho**2
     &            -52929.d0*rho+3149.d0)*r
     &     -(18432.d0*rho**13+172032.d0*rho**12
     &       +751104.d0*rho11+2294016.d0*rho10
     &       +5449472.d0*rho9+9633856.d0*rho8
     &       +11960512.d0*rho7+9817680.d0*rho6
     &       +4782872.d0*rho5+907196.d0*rho4
     &       -257442.d0*rho3-137050.d0*rho**2
     &       -4813.d0*rho+451.d0)*(2.d0*rho+3.d0)
     &     )*etab**4

          teta3 = (256.d0*(2.d0-r)*(16.d0*rho3+24.d0*rho**2
     &         +12.d0*rho+1.d0)*r**7
     &     -128.d0*(192.d0*rho6+192.d0*rho5
     &         -480.d0*rho4-1596.d0*rho3
     &         -1638.d0*rho**2-676.d0*rho
     &         +19.d0)*r**6
     &     +128.d0*(1344.d0*rho6+2976.d0*rho5
     &         +1008.d0*rho4-5784.d0*rho3
     &         -7834.d0*rho**2-3437.d0*rho
     &         +198.d0)*r**5
     &     +64.d0*(2688.d0*rho8+6144.d0*rho7
     &         -8160.d0*rho6-43080.d0*rho5
     &         -54644.d0*rho4-16706.d0*rho3
     &         +15085.d0*rho**2+10263.d0*rho
     &         -929.d0)*r**4
     &     -32.d0*(6144.d0*rho9+56832.d0*rho8
     &         +130560.d0*rho7+48384.d0*rho6
     &         -267264.d0*rho5-486560.d0*rho4
     &         -345992.d0*rho3-92044.d0*rho**2
     &         +2650.d0*rho-377.d0)*r**3
     &     -8.d0*(9216.d0*rho10-232448.d0*rho9
     &         -1256448.d0*rho8-2428224.d0*rho7
     &         -1331808.d0*rho6+2314912.d0*rho5
     &         +4628128.d0*rho4+3289308.d0*rho3
     &         +941094.d0*rho**2+17232.d0*rho
     &         -8719.d0)*r**2
     &     +8.d0*(24576.d0*rho11-64512.d0*rho10
     &         -970240.d0*rho9-3068160.d0*rho8
     &         -4269952.d0*rho7-1520992.d0*rho6
     &         +3514576.d0*rho5+5456064.d0*rho4
     &         +3202096.d0*rho3+682854.d0*rho**2
     &         -49157.d0*rho-5952.d0)*r
     &     -(36864.d0*rho11-24576.d0*rho10
     &       -946176.d0*rho9-3121920.d0*rho8
     &       -4634496.d0*rho7-2574464.d0*rho6
     &       +1765216.d0*rho5+3792160.d0*rho4
     &       +2305160.d0*rho3+449576.d0*rho**2
     &       -69262.d0*rho-2607.d0)
     &       *(2.d0*rho+3.d0)
     &     )*etab**3

          teta2 = -2.d0 * (128.d0*(r-2.d0)*r**7
     &     +96.d0*(32.d0*rho3+32.d0*rho**2
     &         -2.d0*rho-35.d0)*r**6
     &     -64.d0*(432.d0*rho3+582.d0*rho**2
     &         +171.d0*rho-277.d0)*r**5
     &     -8.d0*(2688.d0*rho5+3552.d0*rho4
     &         -9624.d0*rho3-17988.d0*rho**2
     &         -7706.d0*rho+3411.d0)*r**4
     &     +16.d0*(3072.d0*rho6+18816.d0*rho5
     &         +29472.d0*rho4+8568.d0*rho3
     &         -16428.d0*rho**2-11974.d0*rho
     &         -71.d0)*r**3
     &     +2.d0*(4608.d0*rho7-259072.d0*rho6
     &         -904608.d0*rho5-1174896.d0*rho4
     &         -475504.d0*rho3+231848.d0*rho**2
     &         +222510.d0*rho+11293.d0)*r**2
     &     -4.d0*(12288.d0*rho8-48384.d0*rho7
     &         -424480.d0*rho6-945360.d0*rho5
     &         -923200.d0*rho4-270072.d0*rho3
     &         +173214.d0*rho**2+121135.d0*rho
     &         -51.d0)*r
     &     + (4608.d0*rho8-5568.d0*rho7
     &       -145344.d0*rho6-368496.d0*rho5
     &       -401040.d0*rho4-143908.d0*rho3
     &       +56548.d0*rho**2+51463.d0*rho
     &       -1467.d0)*(2.d0*rho+3.d0)
     &     )*etab**2

          teta1 = 2.d0 * (64.d0*(29.d0-3.d0*r)*r**5
     &     +16.d0*(84.d0*rho**2-325.d0)*r**4
     &     -96.d0*(80.d0*rho3+228.d0*rho**2
     &         +120.d0*rho-77.d0)*r**3
     &     -4.d0*(144.d0*rho4-19456.d0*rho3
     &         -33144.d0*rho**2-16416.d0*rho
     &         +4253.d0)*r**2
     &     +4.d0*(1920.d0*rho5-9456.d0*rho4
     &         -48992.d0*rho3-58296.d0*rho**2
     &         -21104.d0*rho+6229.d0)*r
     &     - (288.d0*rho5+336.d0*rho4
     &       -15216.d0*rho3-22680.d0*rho**2
     &       -11046.d0*rho+2953.d0)
     &       *(2.d0*rho+3.d0)
     &     )*etab

          teta0 = -96.d0*(2.d0*r+2.d0*rho-15.d0)
     &                    *(2.d0*r-2.d0*rho-3.d0)*(2.d0*r-1.d0)

          T3 = teta6+teta5+teta4+teta3+teta2+teta1+teta0
          return
        end

        real*8 function C3(r, rho, etab)
c=======================================================================
c  Direct correlation function (part) c^(1)_3 in spherical geometry
c
c  inputs: r    ... distance from the sphere center
c          rho  ... sphere radius
c          etab ... packing fraction
c
c  output: value of c_3(r) for given rho and etab
c
c=======================================================================


          real*8 r, rho, etab
          complex*16 t(4), sum, lga, cden
          integer k

          real*8 den1,den2,logterm

          call roots4a(etab, rho, t)

          den1 = 2.0d0 * (4.0d0*rho**2 + 6.0d0*rho + 3.0d0)*etab*rho
     &         +1.0d0

c====================================================================
c   region 1: r < rho
c====================================================================
          if (r.lt.rho) then
            C3 = 0.0d0
            return
c====================================================================
c   region 2: rho <= r < rho + 1/2
c====================================================================
          else if ((r .ge. rho) .and. (r .lt. rho + 0.5d0)) then
            den2 = (4.0d0*(r**2 + r + 2.0d0*r*rho - 3.0d0*rho**2
     &        - 5.0d0*rho) - 1.1d1)*etab
     &          *((2.0d0*r - 2.0d0*rho + 1.0d0)**2)
     &        + 1.6d1*(2.0d0*r + 1.0d0)
            logterm = log((2.0d0*r + 1.0d0)/(2.0d0*rho))
            logterm = logterm*(2.0d0*r+1.0d0)*(2.0d0*r-1.0d0)
     &         *(rho**2 + rho + 1.0d0)
            logterm = logterm/(8.0d0*r*rho*(rho+1.0d0))
            sum = (0.0d0, 0.0d0)
            do k = 1,4
              cden = (2.0d0*(t(k)+2.0d0*rho)*(t(k)-rho)-6.0d0*rho-3.0d0)
     &             *etab*(t(k)-rho)+1.0d0
              lga = (2.0d0*r-2.0d0*t(k)+1.0d0)/(2.0d0*rho-2.0d0*t(k))
              sum = sum + Q3(t(k), r, rho, etab)*zlog(lga)/
     &               (3.2d1*r*rho*(rho+1.0d0)*cden*den1**2
     &                 *(etab-1.0d0)**2)
            end do
            C3 = S3(r, rho, etab)/(3.2d1*r*den2*den1**2
     &             *(etab-1.0d0)**2)
            C3 = C3 - logterm + dble(sum)
            return
c====================================================================
c   region 3: rho + 1/2 <= r < rho + 3/2
c====================================================================
          else if ((r .ge. rho+0.5d0).and.(r .lt. rho+1.5d0)) then
            den2 = (4.0d0*(r**2 - r + 2.0d0*r*rho - 3.0d0*rho**2
     &        - 7.0d0*rho) - 1.1d1)*etab
     &          *((2.0d0*r - 2.0d0*rho - 1.0d0)**2)
     &        + 1.6d1*(2.0d0*r - 1.0d0)
            logterm = log((2.0d0*r - 1.0d0)/(2.0d0*rho + 2.0d0))
            logterm = logterm*(2.0d0*r+1.0d0)*(2.0d0*r-1.0d0)
            logterm = logterm*(rho**2 + rho + 1.0d0)
            logterm = logterm/(8.0d0*r*rho*(rho+1.0d0))
            sum = (0.0d0, 0.0d0)
            do k = 1,4
              cden = (2.0d0*(t(k)+2.0d0*rho)*(t(k)-rho)-6.0d0*rho-3.0d0)
     &             *etab*(t(k)-rho)+1.0d0
              lga = (2.0d0*r-2.0d0*t(k)-1.0d0)
     &            /(2.0d0*rho-2.0d0*t(k)+2.0d0)
              sum = sum + Q3(t(k), r, rho, etab)*zlog(lga)/
     &               (3.2d1*r*rho*(rho+1.0d0)*cden*den1**2
     &                 *(etab-1.0d0)**2)
            end do
            C3 = T3(r, rho, etab)/(3.2d1*r*den2*den1**2*(etab-1.0d0)**3)
            C3 = C3 + logterm - dble(sum)
            return
c====================================================================
c   region 4: r >= rho + 3/2
c====================================================================
          else if (r .ge. rho + 1.5d0) then
            C3 = -etab*(1.0d0+etab+etab**2)/((1.0d0-etab)**3)
            return
c====================================================================
c   fallback (should not be reached)
c====================================================================
          else
            C3 = 0.0d0
            return
          end if

          return
        end

c===== end of c_3 calculation ==========================================

        real*8 function Cv1(r, rho, etab)
c=======================================================================
c  Direct correlation function (part) \vec{c}^{(1)}_1 in spherical geometry
c     radial part of the vector function (normal to the wall)
c
c  inputs: r    ... distance from the sphere center
c          rho  ... sphere radius
c          etab ... packing fraction
c
c  output: value of c_N1(r) for given rho and etab
c
c=======================================================================

          integer k
          real*8  r, rho, etab
          complex*16 t(4)

          real*8  r2
          real*8  nsum
          complex*16 sumval, denom, term, tval
          complex*16 logterm

          call roots4b(etab, rho, t)

          r2 = r*r

c====================================================================
c   region 1: r < rho
c====================================================================
          if (r .lt. rho) then
            Cv1 = 0.0d0
            return
c====================================================================
c   region 2: rho <= r < rho + 1/2
c====================================================================
          else if ( (r .ge. rho) .and. (r .lt. rho + 0.5d0) ) then
c  nsum polynomial factor:
            nsum = (4.0d0*r2 - 2.0d0*r - 2.0d0*r*rho
     &              - 26.0d0*rho*rho - 25.0d0*rho - 11.0d0)
     &              * (2.0d0*r - 2.0d0*rho + 1.0d0)
            nsum = - nsum / (8.0d0*r2)

c  sum over roots t(k)
            sumval = (0.0d0, 0.0d0)
            do k = 1, 4
              tval = t(k)

c---------- denominator d = (2 t^2 + 6(t-1)rho - 3)*etab*t + 1 --------
              denom = (2.0d0*tval*tval +
     &                 6.0d0*(tval-(1.0d0,0.0d0))*rho -
     &                 3.0d0) * etab * tval + (1.0d0,0.0d0)

c---------- log term: ln( (2r - 2t - 2rho +1)/(-2t) ) ----------
              logterm = zlog( (2.0d0*r - 2.0d0*tval - 2.0d0*rho + 1.0d0)
     &                         / (-2.0d0*tval) )

c---------- numerator polynomial ---------------------------------------
c  t^3 term
              term = (4.0d0*((4.0d0*rho*rho + 6.0d0*rho + 3.0d0)
     &                *rho*etab) + 4.0d0) * tval*tval*tval

c  t^2 term
              term = term + (
     &           ( (8.0d0*rho*rho + 8.0d0*rho + 4.0d0)*etab*r2
     &           - (8.0d0*rho**4 + 40.0d0*rho**3 + 66.0d0*rho*rho
     &             + 46.0d0*rho + 11.0d0)*etab
     &           + 12.0d0*rho ) * tval*tval )

c  t^1 term
              term = term + (
     &          ( (2.0d0*rho+1.0d0)**2 * (2.0d0*rho-1.0d0)*rho*etab
     &          - 4.0d0*(2.0d0*rho+1.0d0)*rho*etab*r2
     &          - 4.0d0*r2 + 28.0d0*rho*rho + 16.0d0*rho + 7.0d0 )
     &          * tval )

c  t^0 term
              term = term - ( (4.0d0*r2 - 20.0d0*rho*rho
     &            - 16.0d0*rho - 7.0d0)*rho )

c---------- accumulate
              sumval = sumval 
     &               - (3.0d0/(8.0d0*r2)) * term/denom * logterm
            end do

            Cv1 = nsum + dble(sumval)
            return
c====================================================================
c   region 3: rho + 1/2 <= r < rho + 3/2
c====================================================================
          else if ((r .ge. rho+0.5d0) .and. (r .lt. rho+1.5d0)) then
            nsum = (4.0d0*r2 - 2.0d0*r*rho - 26.0d0*rho*rho
     &              - 27.0d0*rho - 12.0d0)
     &              * (2.0d0*r - 2.0d0*rho - 3.0d0)
            nsum = nsum / (8.0d0*r2)

            sumval = (0.0d0,0.0d0)
            do k = 1, 4
              tval = t(k)

c---------- denominator ----------
              denom = (2.0d0*tval*tval +
     &                 6.0d0*(tval-(1.0d0,0.0d0))*rho -
     &                 3.0d0) * etab * tval + (1.0d0,0.0d0)

c---------- log term: ln( (2r -2t -2rho -1)/(2-2t) ) ----------
              logterm = zlog( (2.0d0*r - 2.0d0*tval
     &                      - 2.0d0*rho - 1.0d0)
     &                   / (2.0d0 - 2.0d0*tval) )

c---------- numerator polynomial (same as corrected region-2) ---------

              term = (4.0d0*((4.0d0*rho*rho + 6.0d0*rho + 3.0d0)
     &                *rho*etab) + 4.0d0) * tval*tval*tval

c  t^2 term
              term = term + (
     &          ( (8.0d0*rho*rho + 8.0d0*rho + 4.0d0)*etab*r2
     &          - (8.0d0*rho**4 + 40.0d0*rho**3 + 66.0d0*rho*rho
     &            + 46.0d0*rho + 11.0d0)*etab
     &          + 12.0d0*rho ) * tval*tval )

c  t^1 term
              term = term + (
     &          ( (2.0d0*rho+1.0d0)**2*(2.0d0*rho-1.0d0)*rho*etab
     &          - 4.0d0*(2.0d0*rho+1.0d0)*rho*etab*r2
     &          - 4.0d0*r2 + 28.0d0*rho*rho + 16.0d0*rho + 7.0d0)
     &          * tval )

c  t^0 term
              term = term - ( (4.0d0*r2 - 20.0d0*rho*rho
     &            - 16.0d0*rho - 7.0d0)*rho )

c---------- accumulate
              sumval = sumval 
     &               + (3.0d0/(8.0d0*r2)) * term/denom * logterm
            end do

            Cv1 = nsum + dble(sumval)
            return
c====================================================================
c   otherwise (r > rho + 3/2)
c====================================================================
          else
            Cv1 = 0.0d0
            return
          end if

          return
        end

        real*8 function Cv2(r, rho, etab)
c=======================================================================
c  Direct correlation function (part) \vec{c}^{(1)}_2 in spherical geometry
c     radial part of the vector function (normal to the wall)
c
c  inputs: r    ... distance from the sphere center
c          rho  ... sphere radius
c          etab ... packing fraction
c
c  output: value of c_N2(r) for given rho and etab
c
c=======================================================================


          real*8 r, rho, etab
          real*8 rr, rh2, rh3, rh4
          real*8 term1, term2, term3
          integer k

          complex*16 t(4)
          complex*16 sumval, logarg, logterm
          complex*16 termt3, termt2, termt1, termt0, tval
          complex*16 den_t, den_t4

          call roots4b(etab, rho, t)

          rr = r*r
          rh2 = rho*rho
          rh3 = rh2*rho
          rh4 = rh2*rh2
          sumval = (0.0d0, 0.0d0)

c====================================================================
c   region 1: r < rho
c====================================================================
          if (r .lt. rho) then
            Cv2 = 0.0d0
            return
c====================================================================
c   region 2: rho <= r < rho + 1/2
c====================================================================
          else if ( r .ge. rho  .and. r .lt. rho + 0.5d0 ) then

c----- term 1 --------------------------------------------------
            term1 =
     &       -3d0*(2d0*r + 2d0*rho + 3d0)*(2d0*r - 2d0*rho + 1d0)
     &          *(2d0*r + 1d0)*etab
     &       / ( ((4d0*rr + 4d0*r + 8d0*r*rho - 12d0*rh2
     &           -20d0*rho - 11d0)*etab*(2d0*r - 2d0*rho + 1d0)**2
     &           + 16d0*(2d0*r + 1d0))*r )

c----- term 2 --------------------------------------------------
            term2 =
     &      -(4d0*rr - 2d0*r*(rho + 1d0) - 26d0*rh2 - 25d0*rho - 35d0)
     &        *(2d0*r - 2d0*rho + 1d0) / (8d0*rr)

c----- summation ----------------------------------------------------
            sumval = (0d0, 0d0)

            do k = 1, 4
              tval = t(k)

c----- den_t = (2 t^2 + 6 (t-1) rho - 3)*etab*t + 1
              den_t = ((2d0*tval*tval + 6d0*(tval - 1d0)*rho - 3d0)
     &          *etab*tval) + 1d0
              den_t4 = 4d0*den_t

c----- term t^3
              termt3 = ((8d0*rh3 + 12d0*rh2 + 6d0*rho - 3d0)*etab + 2d0)
     &           /den_t

c----- term t^2
              termt2 = (4d0*(rh2 + rho + 1d0)*etab*rr
     &            - (4d0*rh4 + 20d0*rh3 + 55d0*rh2 
     &              + 54d0*rho + 17d0)*etab
     &            + 6d0*rho) / den_t

c----- term t^1
              termt1 = -(
     &           4d0*(4d0*rh2 - 2d0*rho - 1d0)*etab*rr
     &         - (8d0*rh3 - 24d0*rh2 + 1d0)*(2d0*rho + 1d0)*etab
     &         + 8d0*rr - 56d0*rh2 - 32d0*rho - 46d0 ) / den_t4

c----- term t^0
              termt0 =
     &        ( ((4d0*rr - 4d0*rh2 + 1d0)*etab*(2d0*rho + 1d0)
     &          - 8d0*rr + 40d0*rh2 + 32d0*rho + 46d0)*rho ) / den_t4

c----- log argument
              logarg = (2d0*r - 2d0*tval - 2d0*rho + 1d0)/(-2d0*tval)
              logterm = zlog(logarg)

c----- add contributions
              sumval = sumval +
     &          (termt3*tval**3 + termt2*tval**2 
     &           + termt1*tval + termt0)*logterm
            end do

            Cv2 = term1 + term2 - 3d0*dble(sumval)/(4d0*rr)
            return
c====================================================================
c   region 3: rho + 1/2 <= r < rho + 3/2
c====================================================================
          else if ( r .ge. rho+0.5d0 .and. r .lt. rho+1.5d0 ) then

c----- term 1 --------------------------------------------------
            term1 =
     &       (4d0*rr - 2d0*r*(rho - 1d0) - 26d0*rh2 - 23d0*rho - 29d0)
     &       *(2d0*r - 2d0*rho - 1d0)/(8d0*rr)

c----- term 2 --------------------------------------------------
            term2 =
     &      -((12d0*rr*(rho+2d0) - 12d0*rh3 
     &        - 96d0*rh2 - 105d0*rho - 94d0)
     &          *etab - 12d0*rr + 60d0*rh2 + 72d0*rho + 85d0)
     &        / (8d0*(etab - 1d0)*rr)

c----- term 3 --------------------------------------------------
            term3 =
     &      -3d0*(
     &       (16d0*r**4 + 16d0*r**3*rho
     &        - 8d0*(10d0*rh2 + 11d0*rho + 6d0)*rr
     &        + 12d0*(4d0*rh2 + 8d0*rho + 3d0)*rho*r
     &        + (6d0*rho + 11d0)*(2d0*rho + 1d0)**2 )
     &        *etab*(2d0*r - 2d0*rho - 1d0)
     &        + 16d0*(2d0*r + 1d0)*(2d0*r - 1d0)
     &      )
     &      / (4d0*((4d0*rr - 4d0*r + 8d0*r*rho
     &              - 12d0*rh2 - 28d0*rho - 11d0)
     &              *etab*(2d0*r - 2d0*rho - 1d0)**2
     &              + 16d0*(2d0*r - 1d0))*rr)

c----- summation ----------------------------------------------------
            sumval = (0d0, 0d0)

            do k = 1, 4
              tval = t(k)

              den_t = ((2d0*tval*tval + 6d0*(tval - 1d0)*rho - 3d0)
     &          *etab*tval) + 1d0
              den_t4 = 4d0*den_t

              termt3 = ((8d0*rh3 + 12d0*rh2 + 6d0*rho - 3d0)*etab + 2d0)
              termt3 = termt3/den_t

              termt2 = (4d0*(rh2 + rho + 1d0)*etab*rr
     &            - (4d0*rh4 + 20d0*rh3 + 55d0*rh2 
     &              + 54d0*rho + 17d0)*etab
     &            + 6d0*rho ) / den_t

              termt1 = -(
     &         4d0*(4d0*rh2 - 2d0*rho - 1d0)*etab*rr
     &        - (8d0*rh3 - 24d0*rh2 + 1d0)*(2d0*rho + 1d0)*etab
     &        + 8d0*rr - 56d0*rh2 - 32d0*rho - 46d0 ) / den_t4

              termt0 =
     &        (((4d0*rr - 4d0*rh2 + 1d0)*etab*(2d0*rho + 1d0)
     &         - 8d0*rr + 40d0*rh2 + 32d0*rho + 46d0)*rho) / den_t4

c----- region 2 log argument
              logarg = (2d0*r - 2d0*tval - 2d0*rho - 1d0)
     &               /(2d0 - 2d0*tval)
              logterm = zlog(logarg)

              sumval = sumval +
     &         (termt3*tval**3 + termt2*tval**2 
     &         + termt1*tval + termt0)*logterm
            end do

            Cv2 =
     &         term1 + term2 + term3 + 3d0*dble(sumval)/(4d0*rr)
            return
c====================================================================
c   otherwise (r > rho + 3/2)
c====================================================================
          else
            Cv2 = 0.0d0
            return
          endif
        end

c===================================================================
c*
c* Direct correlation function c^(1)_HS
c*     Hard-sphere contribution sum
c*
c* Parameters:
c*     r       ... distance from the sphere center
c*     rho     ... radius of the sphere (spherical wall)
c*     etab    ... bulk density fraction
c*
c* Output:
c*      CHS    ... sum of C3, C2, C1, C0, Cv1 and Cv2
c*
c===================================================================
        real*8 function CHS(r, rho, etab)

          real*8 r, rho, etab

          CHS = C3(r, rho, etab) + C2(r, rho, etab) + C1(r, rho, etab)
     &      + C0(r, rho, etab) + Cv1(r, rho, etab) + Cv2(r, rho, etab)

          return
        end

c*******************************************************************
c*                                                                 *
c* correlation function c^(1)_att (piecewise defined)              *
c*     attractive contribution to direct correlation function due  *
c*     to cut lj potential                                         *
c*     spherical geometry                                          *
c*                                                                 *
c* parameters:                                                     *
c*     r       ... distance from the sphere center                 *
c*     rho     ... radius of the sphere (spherical wall)           *
c*     etab    ... bulk packing fraction                           *
c*     beta    ... inverse temperature 1/kt                        *
c*     epsilon ... depth of the lj potential well                  *
c*     rcut    ... cut-off distance of the lj potential            *
c*                                                                 *
c* output:                                                         *
c*     catt1   ... contribution of attractive interparticle poten- *
c*                 tial to direct correlation function             *
c*                                                                 *
c*******************************************************************
        real*8 function Catt(r, rho, etab, beta, epsilon, rcut)
c
          real*8 r, rho, etab, beta, epsilon, rcut
c
          if (r .ge. (0.5d0 + rcut + rho)) then
            Catt = -3.2d1*etab*beta*epsilon*
     &           (1.0d0/(rcut**3)-1.0d0)
            return
          else if (r .ge. (rho + 1.5d0)) then
            Catt = 6.0d0*etab*beta*epsilon*
     &           (8.0d0*(2.0d0*rcut**3-1.0d0)/(3.0*rcut**3)
     &           -r/(rcut**4)+((2.0d0*rho+1.0d0)**2-8.0d0*rcut**2)/
     &           (4.0d0*r*rcut**4)
     &           +(4.0d0*(6.0d0*rho-2.0d0*r+3.0d0))/
     &           (3.0d0*r*(2.0d0*rho-2.0d0*r+1.0d0)**3))
            return
          else if (r .ge. rho) then
            Catt = 6.0d0*etab*beta*epsilon*
     &           (8.0d0*(rcut**3-1.0d0)/(3.0*rcut**3)
     &           +r*(rcut**4-1.0d0)/(rcut**4)
     &           +((2.0d0*rho+1.0d0)**2-8.0d0*rcut**2)/
     &           (4.0d0*r*rcut**4)
     &           -(4.0d0*rho**2+4.0d0*rho-7.0d0)/(4.0d0*r))
            return
          else
            Catt = 0.0d0
            return
          end if
          return
        end
c


c*******************************************************************
c*                                                                 *
c* Bulk (reduced) chemical potential as a function of density      *
c*    for hard-sphere fluid                                        *
c*                                                                 *
c* parameters:                                                     *
c*     etab ... bulk packing fraction                              *
c*                                                                 *
c* output:                                                         *
c*     BetabMuHS   ... chemical potential \betab\mu                *
c*                                                                 *
c*******************************************************************
        real*8 function BetaMuHS(etab)
c
          real*8 :: etab
c
          BetaMuHS=(5.0d0*etab**3-1.3d1*etab**2+1.4d1*etab)/2.0d0/
     &           (1.0d0-etab)**3-log(1.0d0-etab)+log(etab*6.0d0/m_pi)
c
          return
        end


c*******************************************************************
c*                                                                 *
c* Bulk (reduced) chemical potential as a function of density      *
c*    for hard-sphere fluid with cut lj attraction                 *
c*                                                                 *
c* parameters:                                                     *
c*     etab    ... bulk packing fraction                           *
c*     betab   ... strength of attraction                          *
c*     epsilon ... depth of the lj potential well                  *
c*     rcut    ... cut-off distance of the lj potential            *
c*                                                                 *
c* output:                                                         *
c*     BetaMuLJ   ... chemical potential \betab\mu                *
c*                                                                 *
c*******************************************************************
        real*8 function BetaMuLJ(etab, betab, epsilon, rcut)
c
          real*8 :: etab, betab, epsilon, rcut
c
          BetaMuLJ=(5.0d0*etab**3-1.3d1*etab**2+1.4d1*etab)/2.0d0/
     &           (1.0d0-etab)**3-log(1.0d0-etab)+log(etab*6.0d0/m_pi)
     &           +3.2d1*etab*betab*epsilon*(1.0d0/(rcut**3)-1.0d0)
c
          return
        end


c*******************************************************************
c*                                                                 *
c* External potential (reduced) inducing constant density profile  *
c*    for hard-sphere fluid in a spherical geometry                   *
c*                                                                 *
c* parameters:                                                     *
c*     r    ... distance from the sphere center                    *
c*     rho  ... sphere (spherical wall) radius                     *
c*     etab ... bulk packing fraction                              *
c*                                                                 *
c* output:                                                         *
c*     BetaVHS ... reduced external potential \betab v_ext        *
c*                                                                 *
c*******************************************************************
        real*8 function BetaVHS(r, rho, etab)
c
          real*8 :: r, rho, etab
c
          if (r .lt. (rho + 0.5d0)) then
            BetaVHS = m_inf ! effectively infinite potential
          else
            BetaVHS=BetaMuHS(etab) + CHS(r, rho, etab)
     &      - log(etab*6.0d0/m_pi)
          end if
c
          return
        end


c*******************************************************************
c*                                                                 *
c* external potential (reduced) inducing constant density profile  *
c*    for hard-sphere fluid with cut lj attraction                 *
c*    spherical geometry                                              *
c*                                                                 *
c* parameters:                                                     *
c*     r       ... distance from the sphere center                 *
c*     rho     ... radius of the sphere (spherical wall)           *
c*     etab    ... bulk packing fraction                           *
c*     betab    ... inverse temperature 1/kt                       *
c*     epsilon ... depth of the lj potential well                  *
c*     rcut    ... cut-off distance of the lj potential            *
c*                                                                 *
c* output:                                                         *
c*     BetaVLJ ... reduced external potential \betab v_ext          *
c*                                                                 *
c*******************************************************************
        real*8 function BetaVLJ(r, rho, etab, betab, epsilon, rcut)
c
          real*8 :: r, rho, etab, betab, epsilon, rcut
c
          if (r .lt. (rho + 0.5d0)) then
            BetaVLJ = m_inf ! effectively infinite potential
          else
            BetaVLJ=BetaMuLJ(etab, betab, epsilon, rcut)
     &              + CHS(r, rho, etab) - log(etab*6.0d0/m_pi)
     &              + Catt(r, rho, etab, betab, epsilon, rcut)
          end if
c
          return
        end


c*******************************************************************
c*                                                                 *
c* Boltzmann factor inducing constant density profile              *
c*    boltzmann factor exp(-betab v_ext)                           *
c*    for hard-sphere fluid in a spherical geometry                *
c*                                                                 *
c* parameters:                                                     *
c*     r    ... distance from the sphere center                    *
c*     rho  ... sphere (spherical wall) radius                     *
c*     etab ... bulk packing fraction                              *
c*     etab ... bulk packing fraction                              *
c*                                                                 *
c* output:                                                         *
c*     BoltzFHS ... boltzmann factor exp(-betab v_ext)              *
c*                                                                 *
c*******************************************************************
        real*8 function BoltzFHS(r, rho, etab)
c
          real*8 :: r, rho, etab
c
          if (r .lt. (rho + 0.5d0)) then
            BoltzFHS = 0.0d0 ! effectively infinite potential
          else
            BoltzFHS=BetaMuHS(etab) + CHS(r, rho, etab)
     &               - log(etab*6.0d0/m_pi)
            BoltzFHS=exp(-BoltzFHS)
          end if
c
          return
        end


c*******************************************************************
c*                                                                 *
c* Boltzmann factor inducing constant density profile              *
c*    boltzmann factor exp(-betab v_ext)                           *
c*    for hard-sphere fluid with cut lj attraction                 *
c*    spherical geometry                                           *
c*                                                                 *
c* parameters:                                                     *
c*     r       ... distance from the sphere center                 *
c*     rho     ... sphere (spherical wall) radius                  *
c*     etab    ... bulk packing fraction                           *
c*     etab    ... bulk packing fraction                           *
c*     betab   ... inverse temperature 1/kt                        *
c*     epsilon ... depth of the lj potential well                  *
c*     rcut    ... cut-off distance of the lj potential            *
c*                                                                 *
c* output:                                                         *
c*     BoltzFLJ ... boltzmann factor exp(-betab v_ext)             *
c*                                                                 *
c*******************************************************************
        real*8 function BoltzFLJ(r, rho, etab, beta, epsilon, rcut)
c
          real*8 :: r, rho, etab, beta, epsilon, rcut
c
          if (r .lt. (rho + 0.5d0)) then
            BoltzFLJ = 0.0d0 ! effectively infinite potential
          else
            BoltzFLJ=BetaMuLJ(etab, beta, epsilon, rcut)
     &              + CHS(r, rho, etab)- log(etab*6.0d0/m_pi)
     &              + Catt(r, rho, etab, beta, epsilon, rcut)
            BoltzFLJ=exp(-BoltzFLJ)
          end if
c
          return
        end

      end module wall_sph

