C*******************************************************************
C*                                                                 *
C* Module wall_pl                                                  *
C* - external potential of a planar wall generating constant den-  *
C*   sity profile                                                  *
C* - classical DFT with original Rosenfeld's FMT functional        *
C* - provides analytical expression for the direct correlation     *
C*   function c(1) for hard spheres assuming constant density      *
C* - attractive interactions in the fluid due to cut LJ potential  *
C*   can be assumed by adding analytical expression for the approp.*
C*   contribution to the direct correlation function c(1)_att      *
C*                                                                 *
C* Author: Jiří Janek                                              *
C* Date: January 2026                                              *
C*                                                                 *
C* Compilation:                                                    *
C*   > gfortran -Wall -Wno-tabs -c -O3 wall_pl.for                 *
C*                                                                 *
C*******************************************************************
      module wall_pl

        implicit none

        real*8, parameter :: m_pi=3.141592653589793D0
        real*8, parameter :: m_inf=1.0d6

      contains

C*******************************************************************
C*                                                                 *
C* Roots of the cubic equation (needed for summation)              *
C*     2*eta*t**3-3*eta*t**2+1 = 0                                 *
C*                                                                 *
C* Parameters:                                                     *
C*     etab ... bulk packing fraction                              *
C*                                                                 *
C* Output:                                                         *
C*     t(3) ... three complex roots ot the equation                *
C*                                                                 *
C*******************************************************************
        subroutine Roots(etab, t)
C
          real*8, INTENT(IN) :: etab
          complex*16, INTENT(OUT) :: t(3)
          complex*16 :: aux1, aux2, aux3
C
C   auxiliary variables
C     aux1 = ((eta-2+2*sqrt(1-eta))*eta**2)**(1/3)
C     aux2 = -aux1/4/eta-eta/4/aux1+1/2
C     aux3 = sqrt(3)/2*(aux1/2/eta-eta/2/aux1)
C
          aux1=(((1.0D0, 0.0D0)*etab-2.0D0
     &    +2.0D0*sqrt(1.0D0-etab))*etab**2)**(1.0D0/3.0D0)
          aux2 = -aux1/(4.0D0*etab)-(1.0D0,0.0D0)*etab/(4.0D0*aux1)
     &     +0.5D0
          aux3 = 0.5D0*sqrt(3.0D0)*(aux1/(2.0D0*etab)
     &    -(1.0D0,0.0D0)*etab/(2.0D0*aux1))
C
C   roots
C     t1 = aux1/2/etab+etab/2/aux1+1/2
C     t2 = aux2 + I*aux3
C     t3 = aux2 - I*aux3
C
          t(1) = aux1/(2.0D0*etab)+etab/(2.0D0*aux1)+0.5D0
          t(2) = aux2 + (0.0D0,1.0D0)*aux3
          t(3) = aux2 - (0.0D0,1.0D0)*aux3
C
          return
        end


C*******************************************************************
C*                                                                 *
C* Correlation function c(1) for z > 3/2                           *
C*                                                                 *
C* Parameters:                                                     *
C*     z    ... distance from the wall                             *
C*     etab ... bulk packing fraction                              *
C*                                                                 *
C* Output:                                                         *
C*     c1   ... direct correlation function (for z > 3/2)          *
C*                                                                 *
C*******************************************************************
        real*8 function CHSa(etab)
C
          real*8 :: etab
C
C     c1 for z > 3R (z > 3/2) (independent on z)
          CHSa =-(5.0D0 * etab**3 - 1.3D1 * etab**2 + 1.4D1 * etab)
     &    /(2.0D0 * (1.0D0 - etab)**3) + log(1.0D0 - etab)
C
          return
        end


C*******************************************************************
C*                                                                 *
C* Correlation function c(1) for 1/2 < z < 3/2                     *
C*                                                                 *
C* Parameters:                                                     *
C*     z    ... distance from the wall                             *
C*     etab ... bulk packing fraction                              *
C*                                                                 *
C* Output:                                                         *
C*     c1   ... direct correlation function (for z in [1/2,3/2])   *
C*                                                                 *
C*******************************************************************
        real*8 function CHSb(z, etab)
C
          real*8 :: z, etab, mcbc
          complex*16 :: tt(3), sum, t
          real*8 :: aux, aux1, z2, z3, z4, z5, z6
          integer*4 :: i
          call Roots(etab, tt)

          aux1 = (2.0D0 + (2.0D0*z - 1.0D0)**2 *(z - 2.0D0)*etab)
          aux = (etab-1.0D0)**3 * aux1
          z2 = z*z
          z3 = z2*z
          z4 = z2*z2
          z5 = z3*z2
          z6 = z3*z3
          sum = (0.0D0,0.0D0)
C
C   sum over roots
C
          do 980 i = 1, 3
            t = tt(i)
            sum = sum + (zlog((1.0D0,0.0D0)-t) - zlog((2.0D0,0.0D0)*z
     &      -2.0D0*t-(1.0D0,0.0D0)))*(2.0D0*z - 1.0D0)
     &      /((9.6D1,0.0D0)*(t-(1.0D0,0.0D0))*t*(etab - 1.0D0)**2)
     &      * ((4.8D1*t*t - 1.6D1*t -(1.6D1,0.0D0)*z - 8.0D0)*etab**2
     &      -(9.6D1*t*t + 3.6D1*t*z - 7.0D1*t - 2.2D1*z - 2.3D1)*etab
     &      +4.8D1*t*t - (2.4D1,0.0D0)*z + (1.2D1,0.0D0))
  980     continue
C
C   non-sum part
C
          mcbc =
     &    (6.4D1*z6-3.2D2*z5+5.92D2*z4-5.44D2*z3+3.24D2*z2
     &    -1.08D2*z+1.3D1)*etab**4/8.0D0/aux
     &    +(1.28D2*z6-8.48D2*z5+2.176D3*z4-2.488D3*z3+1.032D3*z2
     &    -1.41D2*z-3.6D1)*etab**3/1.6D1/aux
     &    +(1.28D2*z6-8.16D2*z5+1.824D3*z4-1.864D3*z3+1.032D3*z2
     &    -1.23D2*z+9.0D1)*etab**2/1.6D1/aux
     &    +(4.0D0*z3-1.2D1*z2+3.0D0*z-5.0D0)*etab/aux
     &    -(z+0.5D0)*log(1.0D0-etab)
     &    +(z-0.5D0)*(log(aux1)+2.0D0*log(2.0D0))
C
          CHSb = -(mcbc + dble((1.0D0, 0.0D0)*sum))
          return
        end


C*******************************************************************
C*                                                                 *
C* Correlation function c(1) for 0 < z < 1/2                       *
C*                                                                 *
C* Parameters:                                                     *
C*     z    ... distance from the wall                             *
C*     etab ... bulk packing fraction                              *
C*                                                                 *
C* Output:                                                         *
C*     c1   ... direct correlation function (for z in [0,1/2])     *
C*                                                                 *
C*******************************************************************
        real*8 function CHSc(z, etab)
C
          real*8 :: z, etab
          complex*16 :: tt(3), sum, t
          real*8 :: mccc, z2, z3, z4, z5, z6
          integer*4 :: i
          call Roots(etab, tt)
          z2 = z*z
          z3 = z2*z
          z4 = z2*z2
          z5 = z3*z2
          z6 = z3*z3
          sum = (0.0D0,0.0D0)
C
C   sum over roots
C
          do 981 i = 1, 3
            t = tt(i)
            sum = sum+
     &      (zlog(-t)-zlog((2.0D0,0.0D0)*z-2.0D0*t+(1.0D0,0.0D0)))
     &      *(2.0D0*z-1.0D0)/(9.6D1*(t-1.0D0)*(etab-1.0D0)**2*t)
     &      *((-4.8D1*t*t+1.6D1*t+1.6D1*z+8.0D0)*etab**2
     &      +(3.6D1*t*z-2.2D1*z+9.6D1*t*t-7.0D1*t-2.3D1)*etab
     &      -4.8D1*t*t+2.4D1*z-1.2D1)
  981     continue
C
C   non-sum part
C
          mccc = 2.0D0*etab*(2.0D0*z+1.0D0)
     &    /((etab-1.0D0)**2*(2.0D0+(2.0D0*z+1.0D0)**2*(z-1.0D0)*etab))
     &    *((z4-0.5D0*z3-7.5D-1*z2+1.25D-1*z+1.25D-1)*etab**3
     &    -(2.5D0*z4-1.25D0*z3-1.125D0*z2-4.375D-1*z-2.5D-1)*etab**2
     &    -(7.5D-1*z4-2.625D0*z3-1.3125D0*z2+3.65625D0*z
     &    +1.03125D0)*etab + 1.5D0)
     &    -(z+0.5D0)*log(2.0D0+(2.0D0*z+1)**2*(z-1.0D0)*etab)
     &    -2.0D0*(z-1.0D0)*log(2.0D0)
C      c1 for z>0 and z<1/2
          CHSc = -(mccc + dble((1.0D0,0.0D0)*sum))
          return
        end


C*******************************************************************
C*                                                                 *
C* Direct correlation function c^(1)_hs (piecewise defined)        *
C*     Call CHSa, CHSb or CHSc depending on the value of z         *
C*     Hard-sphere contribution to direct correlation function     *
C*                                                                 *
C* Parameters:                                                     *
C*     z    ... distance from the wall                             *
C*     etab ... bulk packing fraction                              *
C*                                                                 *
C* Output:                                                         *
C*     CHS  ... direct correlation function                        *
C*                                                                 *
C*******************************************************************
        real*8 function CHS(z, etab)
C
          real*8 :: z, etab

C      WRITE(*,*) r, z, x
          if (z .gt. 1.5D0) then
            CHS = CHSa(etab)
          else if (z .gt. 0.5000D0) then
            CHS = CHSb(z, etab)
          else if (z .ge. 0.0D0) then
            CHS = CHSc(z, etab)
          else
            CHS =0.0D0
          end if
          return
        end


C*******************************************************************
C*                                                                 *
C* Padé approximation of direct correlation function c^(1)_hs      *
C*     Call CHSa, CHSb or CHSc depending on the value of z         *
C*     Hard-sphere contribution to direct correlation function     *
C*                                                                 *
C* Parameters:                                                     *
C*     z    ... distance from the wall                             *
C*     etab ... bulk packing fraction                              *
C*                                                                 *
C* Output:                                                         *
C*     CHSPade  ... direct correlation function                    *
C*                                                                 *
C*******************************************************************
        real*8 function CHSPade(z, etab)
C
          real*8 :: z, etab
          real*8 :: A0, A1, A2, B1, B2
          real*8 :: aux1, aux2, aux3, ln1me, den1, den2, den3
          real*8 :: e2, e3, e4, e5, e6
          real*8 :: zarg

          e2 = etab*etab
          e3 = e2*etab
          e4 = e2*e2
          e5 = e2*e3
          e6 = e3*e3

          den1 = 2.0d0*(1.0d0 - etab)**3 ! denominator for A0, A1 and A2
          den2 = 3.0d0*(1.0d0 - etab)*(etab + 2.0d0) ! denom. for A1, A2 and B1
          den3 = e2 - 2.0d0*etab + 4.0d0 ! denom for A2 and B2

          aux1 = etab*(5.0d0*e2 - 1.3d1*etab + 1.4d1)/den1
          aux2 = -2.5d1*e4 + 4.0d0*e3 - 9.0d0*e2 - 2.0d1*etab - 4.0d0
          aux3 = 1.07d2*e6 - 3.45d2*e5 + 5.01d2*e4 - 2.53d2*e3
     &           - 1.14d2*e2 + 6.6d2*etab - 2.32d2
          ln1me = log(1.0d0 - etab)

          B1 = (1.3d1*e2 + 7.0d0*etab - 2.0d0)/den2
          A0 = aux1 - ln1me
          A1 = -B1*ln1me + aux1*B1
          A2 = (aux2*ln1me + aux3*etab/den1)/den2/den3
          B2 = -aux2/den2/den3

          zarg = z - 1.5d0


C      WRITE(*,*) r, z, x
          if (z .gt. 1.5D0) then
            CHSPade = CHSa(etab) ! bulk limit - simple expression
          else if (z .gt. 0.5000D0) then
            CHSPade = -(A0 + A1*zarg + A2*zarg**2)
     &              /(1.0d0 + B1*zarg + B2*zarg**2)
          else if (z .ge. 0.0D0) then
            CHSPade = -(A0 + A1*zarg + A2*zarg**2)
     &              /(1.0d0 + B1*zarg + B2*zarg**2) ! irrelevant for external potential
          else
            CHSPade =0.0D0
          end if
          return
        end          


C*******************************************************************
C*                                                                 *
C* Correlation function c^(1)_att (piecewise defined)              *
C*     Attractive contribution to direct correlation function due  *
C*     to cut LJ potential                                         *
C*                                                                 *
C* Parameters:                                                     *
C*     z       ... distance from the wall                          *
C*     etab    ... bulk packing fraction                           *
C*     beta    ... inverse temperature 1/kT                        *
C*     epsilon ... depth of the LJ potential well                  *
C*     rcut    ... cut-off distance of the LJ potential            *
C*                                                                 *
C* Output:                                                         *
C*     catt1   ... contribution of attractive interparticle poten- *
C*                 tial to direct correlation function             *
C*                                                                 *
C*******************************************************************
        real*8 function Catt(z, etab, beta, epsilon, rcut)
C
          real*8 :: z, etab, beta, epsilon, rcut
C
          if (z .gt. (0.5d0+rcut)) then
            Catt = -3.2d1*etab*beta*epsilon*
     &           (1.0d0/(rcut**3)-1.0d0)
          else if (z .gt. 1.5d0) then
            Catt = -9.6d1*etab*beta*epsilon*
     &           (1.0d0/(6.0*rcut**3)-1.0d0/3.0d0
     &           +(2.0d0*z-1.0d0)/(1.6d1*rcut**4)
     &           -1.0d0/(3.0d0*(1.0d0-2.0d0*z)**3))
          else if (z .ge. 0.0D0) then
            Catt = -4.8d1*etab*beta*epsilon*
     &           (1.0d0/(3.0d0*rcut**3)-5.0d0/2.4d1
     &           +(2.0d0*z-1.0d0)/(8.0d0*rcut**4) - z/4.0d0)
          else
            Catt = 0.0D0
          end if
          return
        end
C


C*******************************************************************
C*                                                                 *
C* Bulk (reduced) chemical potential as a function of density      *
C*    For hard-sphere fluid                                        *
C*                                                                 *
C* Parameters:                                                     *
C*     etab ... bulk packing fraction                              *
C*                                                                 *
C* Output:                                                         *
C*     BetaMuHS   ... chemical potential \beta\mu                  *
C*                                                                 *
C*******************************************************************
        real*8 function BetaMuHS(etab)
C
          real*8 :: etab
C
          BetaMuHS=(5.0D0*etab**3-1.3D1*etab**2+1.4D1*etab)/2.0D0/
     &           (1.0D0-etab)**3-log(1.0D0-etab)+log(etab*6.0D0/m_pi)
C
          return
        end


C*******************************************************************
C*                                                                 *
C* Bulk (reduced) chemical potential as a function of density      *
C*    For hard-sphere fluid with cut LJ attraction                 *
C*                                                                 *
C* Parameters:                                                     *
C*     etab    ... bulk packing fraction                           *
C*     beta    ... strength of attraction                          *
C*     epsilon ... depth of the LJ potential well                  *
C*     rcut    ... cut-off distance of the LJ potential            *
C*                                                                 *
C* Output:                                                         *
C*     BetaMuLJ   ... chemical potential \beta\mu                  *
C*                                                                 *
C*******************************************************************
        real*8 function BetaMuLJ(etab, beta, epsilon, rcut)
C
          real*8 :: etab, beta, epsilon, rcut
C
          BetaMuLJ=(5.0D0*etab**3-1.3D1*etab**2+1.4D1*etab)/2.0D0/
     &           (1.0D0-etab)**3-log(1.0D0-etab)+log(etab*6.0D0/m_pi)
     &           +3.2d1*etab*beta*epsilon*(1.0d0/(rcut**3)-1.0d0)
C
          return
        end


C*******************************************************************
C*                                                                 *
C* External potential (reduced) inducing constant density profile  *
C*    For hard-sphere fluid in a planar geometry                   *
C*                                                                 *
C* Parameters:                                                     *
C*     z    ... distance from the wall                             *
C*     etab ... bulk packing fraction                              *
C*                                                                 *
C* Output:                                                         *
C*     BetaVHS ... reduced external potential \beta V_ext          *
C*                                                                 *
C*******************************************************************
        real*8 function BetaVHS(z, etab)
C
          real*8 :: z, etab
C
          if (z .lt. 0.5d0) then
            BetaVHS = m_inf ! effectively infinite potential
          else
            BetaVHS=BetaMuHS(etab) + CHS(z, etab)
     &              - log(etab*6.0D0/m_pi)
          end if
C
          return
        end


C*******************************************************************
C*                                                                 *
C* External potential (reduced) inducing constant density profile  *
C*    For hard-sphere fluid in a planar geometry                   *
C*    Padé approximation for c_HS is used                          *
C*                                                                 *
C* Parameters:                                                     *
C*     z    ... distance from the wall                             *
C*     etab ... bulk packing fraction                              *
C*                                                                 *
C* Output:                                                         *
C*     BetaVHSPade ... reduced external potential \beta V_ext      *
C*                                                                 *
C*******************************************************************
        real*8 function BetaVHSPade(z, etab)
C
          real*8 :: z, etab
C
          if (z .lt. 0.5d0) then
            BetaVHSPade = m_inf ! effectively infinite potential
          else
            BetaVHSPade=BetaMuHS(etab) + CHSPade(z, etab)
     &              - log(etab*6.0D0/m_pi)
          end if
C
          return
        end


C*******************************************************************
C*                                                                 *
C* External potential (reduced) inducing constant density profile  *
C*    For hard-sphere fluid with cut LJ attraction                 *
C*    Planar geometry                                              *
C*                                                                 *
C* Parameters:                                                     *
C*     z       ... distance from the wall                          *
C*     etab    ... bulk packing fraction                           *
C*     beta    ... inverse temperature 1/kT                        *
C*     epsilon ... depth of the LJ potential well                  *
C*     rcut    ... cut-off distance of the LJ potential            *
C*                                                                 *
C* Output:                                                         *
C*     BetaVLJ ... reduced external potential \beta V_ext          *
C*                                                                 *
C*******************************************************************
        real*8 function BetaVLJ(z, etab, beta, epsilon, rcut)
C
          real*8 :: z, etab, beta, epsilon, rcut
C
          if (z .lt. 0.5d0) then
            BetaVLJ = m_inf ! effectively infinite potential
          else
            BetaVLJ=BetaMuLJ(etab, beta, epsilon, rcut)
     &              + CHS(z, etab)- log(etab*6.0D0/m_pi)
     &              + Catt(z, etab, beta, epsilon, rcut)
          end if
C
          return
        end


C*******************************************************************
C*                                                                 *
C* External potential (reduced) inducing constant density profile  *
C*    For hard-sphere fluid with cut LJ attraction                 *
C*    Planar geometry                                              *
C*    Padé approximation for c_HS is used                          *
C*                                                                 *
C* Parameters:                                                     *
C*     z       ... distance from the wall                          *
C*     etab    ... bulk packing fraction                           *
C*     beta    ... inverse temperature 1/kT                        *
C*     epsilon ... depth of the LJ potential well                  *
C*     rcut    ... cut-off distance of the LJ potential            *
C*                                                                 *
C* Output:                                                         *
C*     BetaVLJPade ... reduced external potential \beta V_ext      *
C*                                                                 *
C*******************************************************************
        real*8 function BetaVLJPade(z, etab, beta, epsilon, rcut)
C
          real*8 :: z, etab, beta, epsilon, rcut
C
          if (z .lt. 0.5d0) then
            BetaVLJPade = m_inf ! effectively infinite potential
          else
            BetaVLJPade=BetaMuLJ(etab, beta, epsilon, rcut)
     &              + CHSPade(z, etab)- log(etab*6.0D0/m_pi)
     &              + Catt(z, etab, beta, epsilon, rcut)
          end if
C
          return
        end


C*******************************************************************
C*                                                                 *
C* Boltzmann factor inducing constant density profile              *
C*    Boltzmann factor exp(-beta V_ext)                            *
C*    For hard-sphere fluid in a planar geometry                   *
C*                                                                 *
C* Parameters:                                                     *
C*     z    ... distance from the wall                             *
C*     etab ... bulk packing fraction                              *
C*                                                                 *
C* Output:                                                         *
C*     BoltzFHS ... Boltzmann factor exp(-beta V_ext)              *
C*                                                                 *
C*******************************************************************
        real*8 function BoltzFHS(z, etab)
C
          real*8 :: z, etab
C
          if (z .lt. 0.5d0) then
            BoltzFHS = 0.0d0 ! effectively infinite potential
          else
            BoltzFHS = BetaMuHS(etab) + CHS(z, etab) 
     &        - log(etab*6.0D0/m_pi)
            BoltzFHS = exp(-BoltzFHS)
          end if
C
          return
        end


C*******************************************************************
C*                                                                 *
C* Boltzmann factor inducing constant density profile              *
C*    Boltzmann factor exp(-beta V_ext)                            *
C*    For hard-sphere fluid in a planar geometry                   *
C*    Padé approximation for c_HS is used                          *
C*                                                                 *
C* Parameters:                                                     *
C*     z    ... distance from the wall                             *
C*     etab ... bulk packing fraction                              *
C*                                                                 *
C* Output:                                                         *
C*     BoltzFHSPade ... Boltzmann factor exp(-beta V_ext)          *
C*                                                                 *
C*******************************************************************
        real*8 function BoltzFHSPade(z, etab)
C
          real*8 :: z, etab
C
          if (z .lt. 0.5d0) then
            BoltzFHSPade = 0.0d0 ! effectively infinite potential
          else
            BoltzFHSPade = BetaMuHS(etab) + CHSPade(z, etab) 
     &        - log(etab*6.0D0/m_pi)
            BoltzFHSPade = exp(-BoltzFHSPade)
          end if
C
          return
        end



C*******************************************************************
C*                                                                 *
C* Boltzmann factor inducing constant density profile              *
C*    Boltzmann factor exp(-beta V_ext)                            *
C*    For hard-sphere fluid with cut LJ attraction                 *
C*    Planar geometry                                              *
C*                                                                 *
C* Parameters:                                                     *
C*     z       ... distance from the wall                          *
C*     etab    ... bulk packing fraction                           *
C*     beta    ... inverse temperature 1/kT                        *
C*     epsilon ... depth of the LJ potential well                  *
C*     rcut    ... cut-off distance of the LJ potential            *
C*                                                                 *
C* Output:                                                         *
C*     BoltzFLJ ... Boltzmann factor exp(-beta V_ext)              *
C*                                                                 *
C*******************************************************************
        real*8 function BoltzFLJ(z, etab, beta, epsilon, rcut)
C
          real*8 :: z, etab, beta, epsilon, rcut
C
          if (z .lt. 0.5d0) then
            BoltzFLJ = 0.0d0 ! effectively infinite potential
          else
            BoltzFLJ=BetaMuLJ(etab, beta, epsilon, rcut) 
     &              + CHS(z, etab) - log(etab*6.0D0/m_pi)
     &              + Catt(z, etab, beta, epsilon, rcut)
            BoltzFLJ=exp(-BoltzFLJ)
          end if
C
          return
        end

    
C*******************************************************************
C*                                                                 *
C* Boltzmann factor inducing constant density profile              *
C*    Boltzmann factor exp(-beta V_ext)                            *
C*    For hard-sphere fluid with cut LJ attraction                 *
C*    Planar geometry                                              *
C*    Padé approximation for c_HS is used                          *
C*                                                                 *
C* Parameters:                                                     *
C*     z       ... distance from the wall                          *
C*     etab    ... bulk packing fraction                           *
C*     beta    ... inverse temperature 1/kT                        *
C*     epsilon ... depth of the LJ potential well                  *
C*     rcut    ... cut-off distance of the LJ potential            *
C*                                                                 *
C* Output:                                                         *
C*     BoltzFLJPade ... Boltzmann factor exp(-beta V_ext)          *
C*                                                                 *
C*******************************************************************
        real*8 function BoltzFLJPade(z, etab, beta, epsilon, rcut)
C
          real*8 :: z, etab, beta, epsilon, rcut
C
          if (z .lt. 0.5d0) then
            BoltzFLJPade = 0.0d0 ! effectively infinite potential
          else
            BoltzFLJPade=BetaMuLJ(etab, beta, epsilon, rcut) 
     &              + CHSPade(z, etab) - log(etab*6.0D0/m_pi)   
     &              + Catt(z, etab, beta, epsilon, rcut)
            BoltzFLJPade=exp(-BoltzFLJPade)
          end if
C
          return
        end

      end module wall_pl
