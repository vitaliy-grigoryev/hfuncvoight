! This module describes all Global constants, such as pi, m_e etc.
module Fundamental_data
implicit none

integer, parameter :: k_p = 8
real(k_p) :: pi, me, mp, c, e, h, hp, g, kb
real(k_p) :: kev, rel, finestr, mee, thom, r0, ks, stef, lam_e, B_q
real(k_p) :: m_sun

      parameter ( pi = 3.14159265358979d0 )	!the Pi [1]
      parameter ( me = 9.10938189d-28 )		!the electron mass [g]
      parameter ( mp = 1.67262158d-24 )		!the proton mass [g] 
      parameter ( c  = 2.99792458d10 )		!the speed of light [cm s^{-1}] 
      parameter ( e  = 4.8032068d-10 )		!the electron charge [CGS] 
      parameter ( h  = 6.6260755d-27 )		!the Plank constant [erg s]
      parameter ( hp = h/(2.0d0*pi) )		!the modified Plank constant [erg s rad^{-1}]	
      parameter ( G  = 6.67259d-8 )			!the gravitational constant [dyn cm^2 g^{-2}]
      parameter ( kb = 1.380658d-16 )		!the Boltzmann constant [erg K^{-1}]   

      parameter ( kev = 1.0d3*1.602192d-12 )   					!kev [erg/keV]
      parameter ( rel=e*e/(me*c*c) )          					!the classical electron radius [cm] 
      parameter ( finestr=e*e/(hp*c) )          				!the fine structure constant alpha [1]
      parameter ( mee=me*c*c )               					!the electron mass in energy units [erg]
      parameter ( thom=(8.0d0*pi/3.0d0)*rel*rel )  				!the Thomson crossection [cm^2]
      parameter ( r0=2.0d0*pi*rel*rel*mee*c )    				!a typical energy-loss scale [see LR'82] [cm]
      parameter ( ks=me/mp )                  					!the electron-to-proton mass ratio [1]
      parameter ( stef=(pi*kb*kb*kb*kb)/(60.0d0*hp*hp*hp*c*c) )	!the Stefan-Boltzman constant [erg s^{-1} cm^{-2} K^{-4}]
      parameter ( lam_e=hp/(me*c) )           					!the Compton wavelength of an electron [cm]
      parameter ( B_q = me*me*c*c*c/(e*hp))						!the critical magnetic field (electron energy=cyclotrone energy) [G]

      parameter ( m_sun=1.98892d33 )          !the solar mass [g]    
 
end module Fundamental_data
