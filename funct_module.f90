module Funct_module

use Fundamental_data
use Integrating_module
use Global_module

implicit none

integer, parameter, private :: knots_l = 8 ! count of Legendre`s roots for integrating (Gauss` integral) (Nagirner`s count = 8)

real(k_p), parameter, private :: acc = dble(1.0e-30) ! accuracy
real(k_p), parameter, private :: x_as = 10.0_k_p*2.718281828_k_p**t_up ! x-asymptotical in formula (84) = x(knots)
real(k_p), parameter, private :: step = 0.001_k_p ! step on sum of integrals (empirically: 0.001 enough for pages 54-55)

contains


!*************************************
! Calculating function U(a,x) --- see page 21 (Труды астрономической обсерватории XXXI) 
! a is Voight`s parameter
! x is dimensionless frequency [-inf; +inf]
!*************************************
elemental function calc_u(a, x) result(u) ! FIXME FIXME FIXME FIXME FIXME FIXME FIXME FIXME FIXME and intu2 too

implicit none 
real(k_p), intent(in) :: a, x
real(k_p) :: u, uu
real(k_p) :: x_m2, a_2, s1, s2, up, down, term
integer :: n, k

	uu = 0.0_k_p
	if (a .eq. 0.0_k_p) then
		uu = exp(-x**2)/sqrt(pi)
	else
		if (x .gt. 10000.0_k_p*x_middle) then ! asymptotics - formula (27) !x_middle 50000.0_k_p 
			x_m2 = x**(-2)
			a_2 = a**2
!			uu = (a/pi) * x_m2 * ( 1.0_k_p + (1.5_k_p - a_2)*x_m2 + (3.75_k_p - 5.0_k_p*a_2 + a_2*a_2)*x_m2*x_m2 ) ! formula (83)
!write(*,*)'formula x = ',x,' uu=',uu	! for x ~ 10 n = 27; for x ~3e6   n=3		
			
			n = 0
			k = 0
			s1 = 1.0_k_p
			s2 = a
			uu = s1*s2 ! 0-th term
			do while ( (abs(s1*s2) .gt. acc) .or. (n .eq. 0) )
				n = n + 1
				s1 = s1*x_m2*real(n+1)*real(n+2)
				term = a*(0.5_k_p**(2*n))/fact(n) ! k = 0
				s2 = term
				do k = 1, n-1
					term = -term*a_2*4.0_k_p*(n-k)/( (2.0_k_p*k+1) )
					s2 = s2 + term 
				end do
				term = -term*a_2*4.0_k_p/( (2.0_k_p*n+1) ) ! k = n
				s2 = s2 + term
				uu = uu + s1*s2
			end do ! while
			uu = uu/(pi * x**2)
		else
			s1 = 10.0_k_p + acc
			s2 = 10.0_k_p + acc
			down = 0.0_k_p
			up = down + step

			do while ((abs(s1) .gt. acc) .and. (abs(s2) .gt. acc))
				
				if (s1 .gt. acc) then
					s1 = int_ey2(down, up, a, x)
				end if
				if (s2 .gt. acc) then
					s2 = int_ey2(-up, -down, a, x)
				end if
				uu = uu + s1 + s2
				
				down = down + step
				up = up + step
			end do
			uu = uu * a/(pi*sqrt(pi))
		end if
	end if

	u = uu

end function calc_u

!*************************************
! Calculating function integral in function U(a,x) (part of formula 16) from down to up
! a is Voigt`s parameter
! x is dimensionless frequency [-inf; +inf]
!*************************************
pure function int_ey2(down, up, a, x) result(ii)

implicit none 
real(k_p), intent(in) :: down, up, a, x
real(k_p) :: ii 
real(k_p), dimension(1:knots_l) :: y, w

	call knots_and_weights(knots_l, down, up, y, w)
	ii = sum( w(1:knots_l) * exp(-y(1:knots_l)**2)/( (y(1:knots_l) - x)**2 + a**2 ) )
	
end function int_ey2


!*************************************
! Calculating derivative of odd function df/dx ( used for see page 21 (Труды астрономической обсерватории XXXI) )
! n - count of knots
! x - vector of arguments (from 0 to x_max)
! f0 - input values of f in knots. They are continued to negative x in even order
!*************************************
function calc_df_dx(n, x0, f0) result(df)
implicit none 
integer, parameter :: order = 10 ! count of coefficients-1 => accuracy is O(h^order)
integer, parameter :: k = 1 ! first derivative

integer, intent(in) :: n
real(k_p), dimension(0:n), intent(in) :: x0, f0
real(k_p), dimension(0:n) :: df

real(k_p), dimension(-n/10:n) :: f, x
real(k_p), dimension(0:order, 0:n) :: c ! coefficients for calculating ( стр. 234, п.4 [Методы Вычислений (1962)] )
real(k_p), dimension(0:order, 0:order, 0:n) :: matr_c
real(k_p), dimension(0:order, 0:n) :: vec_c
integer :: i, j

	f(0:n) = f0(0:n)
	x(0:n) = x0(0:n)
	forall (i = 1:n/10)
		f(-i) = f0(i)
		x(-i) = -x0(i)
	end forall ! making function odd
	matr_c(0:order, 0, 0:n) = 1.0_k_p
	vec_c(0, 0:n) = 1.0_k_p
	do j = 0, n-order/2
		do i = 1, order
			matr_c(0:order, i, j) = matr_c(0:order, i-1, j)*x(j-order/2:j+order/2)
			if (i .gt. k) then
				vec_c(i,j) = vec_c(i-1,j)*real(i)*x(j)/real(i-k)
			else
				vec_c(i,j) = vec_c(i-1,j)*real(i)
			end if
		end do
	end do
	do j = n-order/2 + 1, n
		do i = 1, order
			matr_c(0:order, i, j) = matr_c(0:order, i-1,j)*x(n-order:n)
			if (i .gt. k) then
				vec_c(i,j) = vec_c(i-1,j)*real(i)*x(j)/real(i-k)
			else
				vec_c(i,j) = vec_c(i-1,j)*real(i)
			end if
		end do
	end do
	vec_c(0:k-1, 0:n) = 0.0_k_p

	forall (i = 0:n)
		c(0:order, i) = CHOOSE(matr_c(0:order,0:order,i), vec_c(0:order,i),order+1)
	end forall

	forall (i = 0:n-order/2) 
		df(i) = sum(c(0:order,i)*f(i-order/2:i+order/2))
	end forall
	
	forall (i = n-order/2+1:n) 
		df(i) = sum(c(0:order,i)*f(n-order:n))
	end forall

	return 
end function calc_df_dx


!*************************************
! Calculating function \int^{\inf}_{x} U(a,x) dx --- see page 21 (Труды астрономической обсерватории XXXI) 
! Since \int_{0}^{\inf} U(a,x)dx = 0.5, \int^{\inf}_{x} = \int_{0}^{\inf} - \int_{0}^{x} = 0.5 - \int_{0}^{x}
! a is Voigt`s parameter
! x is dimensionless frequency [-inf; +inf]
!*************************************
elemental function calc_int_u(a, xin) result(int_u)
implicit none 
real(k_p), intent(in) :: a, xin
real(k_p) :: int_u, int1, int2
integer, parameter :: knots_y = 10!knots/20
real(k_p), dimension(1:knots_l, 1:knots_y) :: w, xx, tt
real(k_p), dimension(0:knots_y) :: y, fy
integer :: knot_brdr, i

	int1 = 0.0_k_p
	fy = 0.0_k_p
	if (xin .le. x_middle) then
		
		forall (i=0:knots_y)
			y(i) = (xin)*dble(i)/dble(knots_y)
		end forall
		do i=1, knots_y
			call knots_and_weights(knots_l, y(i-1), y(i), xx(:,i), w(:,i))
		end do
		forall (i = 1:knots_y)
			fy(i) = sum(w(1:knots_l,i)*calc_u(a,xx(1:knots_l,i)))  
		end forall
		int1 = sum(fy)

	else !(xin .gt. x_middle)
	
		knot_brdr = knot_x_last
		do while ((x(knot_brdr) .lt. xin) .and. (knot_brdr .lt. knots))
			knot_brdr = knot_brdr + 1
		end do
		
		if (xin .ge. 100.0_k_p*x_middle) then !x(knot_x_last)
! 			fy = 0.0_k_p
! 			forall (i=0:knots_y)
! 				y(i) = (log(x(knots)/x_middle) - log(xin/x_middle))*dble(i)/dble(knots_y)
! 			end forall
! 			do i=1, knots_y
! 				call knots_and_weights(knots_l, y(i-1), y(i), tt(:,i), w(:,i))
! 			end do
! 			forall (i = 1:knots_y)
! 				fy(i) = x_middle*sum(w(1:knots_l,i)*exp(tt(1:knots_l,i))*calc_u(a,x_middle*exp(tt(1:knots_l,i))))  
! 			end forall
! 			int2 = sum(fy) ! посчитано то, что и нужно
			int1 = simpson_int_knots_fx(knot_brdr, knot_x_last, x(0:knot_brdr), x_middle, uax(0:knot_brdr))
		else
			forall (i=0:knots_y)
				y(i) = x_middle*dble(i)/dble(knots_y)
			end forall
			do i=1, knots_y
				call knots_and_weights(knots_l, y(i-1), y(i), xx(:,i), w(:,i))
			end do
			forall (i = 1:knots_y)
				fy(i) = sum(w(1:knots_l,i)*calc_u(a,xx(1:knots_l,i)))  
			end forall
			int1 = sum(fy)
			
			fy = 0.0_k_p
			forall (i=0:knots_y)
				y(i) = log(xin/x_middle)*dble(i)/dble(knots_y)
			end forall
			do i=1, knots_y
				call knots_and_weights(knots_l, y(i-1), y(i), tt(:,i), w(:,i))
			end do
			forall (i = 1:knots_y)
				fy(i) = sum(w(1:knots_l,i)*x_middle*exp(tt(1:knots_l,i))*calc_u(a,x_middle*exp(tt(1:knots_l,i))))  
			end forall
			!call knots_and_weights(knots_l, x_middle, x(knot_brdr), xx(1:knots_l,1), w(1:knots_l,1))!xin
			int2 = sum(fy)!sum(w(1:knots_l,1)*calc_u(a,xx(1:knots_l,1))**2)    
   
			int1 = int1 + int2
		end if
		
	end if
	
	int_u = max(0.5_k_p - int1, 0.0_k_p) ! = \int^{\inf}_{0} - \int^{x}_{0} ! FIXME

	return 
end function calc_int_u


!*************************************
! Calculating function \int^{\inf}_{x} U(a,x)^2 dx --- see page 21 (Труды астрономической обсерватории XXXI) 
! Calculating it like previous integral ( see calc_int_u ): int^{\inf}_{x} = int^{\inf}_{0} - int^{x}_{0}
! a is Vooight`s parameter
! xin = x is dimensionless frequency [-inf; +inf]
!*************************************
elemental function calc_int_u2(a, xin) result(int_u)
implicit none 
real(k_p), intent(in) :: a, xin
real(k_p) :: int_u, int1, int2, int_ux0
integer, parameter :: knots_y = 10!knots/20
real(k_p), dimension(1:knots_l, 1:knots_y) :: w, xx, tt
real(k_p), dimension(0:knots_y) :: y, fy
integer :: knot_brdr, i

	int_ux0 = uasqrt20/sqrt(8.0_k_p) ! exact value of integral, if x = 0 -- formula (26), page 10.
	int1 = 0.0_k_p
	fy = 0.0_k_p
	if (xin .le. x_middle) then
		
		forall (i=0:knots_y)
			y(i) = (xin)*dble(i)/dble(knots_y)
		end forall
		do i=1, knots_y
			call knots_and_weights(knots_l, y(i-1), y(i), xx(:,i), w(:,i))
		end do
		forall (i = 1:knots_y)
			fy(i) = sum(w(1:knots_l,i)*calc_u(a,xx(1:knots_l,i))**2)  
		end forall
		int1 = sum(fy)

	else !(xin .gt. x_middle)
		
		if (a .lt. 1.0_k_p) then
			if (xin .ge. 1000.0_k_p*x_middle) then !x(knot_x_last)
				fy = 0.0_k_p
				forall (i=0:knots_y)
					y(i) = (log(x(knots)/x_middle) - log(xin/x_middle))*dble(i)/dble(knots_y)
				end forall
				do i=1, knots_y
					call knots_and_weights(knots_l, y(i-1), y(i), tt(:,i), w(:,i))
				end do
				forall (i = 1:knots_y)
					fy(i) = x_middle*sum(w(1:knots_l,i)*exp(tt(1:knots_l,i))*calc_u(a,x_middle*exp(tt(1:knots_l,i)))**2)  
				end forall
				int2 = sum(fy) ! посчитано то, что и нужно
				int1 = int_ux0 - int2 !
			else
				knot_brdr = knot_x_last
				do while ((x(knot_brdr) .lt. xin) .and. (knot_brdr .lt. knots))
					knot_brdr = knot_brdr + 1
				end do
				int1 = simpson_int_knots_fx(knot_brdr, knot_x_last, x(0:knot_brdr), x_middle, (uax(0:knot_brdr)**2))
			end if
 		else
			forall (i=0:knots_y)
				y(i) = x_middle*dble(i)/dble(knots_y)
			end forall
			do i=1, knots_y
				call knots_and_weights(knots_l, y(i-1), y(i), xx(:,i), w(:,i))
			end do
			forall (i = 1:knots_y)
				fy(i) = sum(w(1:knots_l,i)*calc_u(a,xx(1:knots_l,i))**2)  
			end forall
			int1 = sum(fy)
			
			fy = 0.0_k_p
			forall (i=0:knots_y)
				y(i) = log(xin/x_middle)*dble(i)/dble(knots_y)
			end forall
			do i=1, knots_y
				call knots_and_weights(knots_l, y(i-1), y(i), tt(:,i), w(:,i))
			end do
			forall (i = 1:knots_y)
				fy(i) = sum(w(1:knots_l,i)*x_middle*exp(tt(1:knots_l,i))*calc_u(a,x_middle*exp(tt(1:knots_l,i)))**2)  
			end forall
			!call knots_and_weights(knots_l, x_middle, x(knot_brdr), xx(1:knots_l,1), w(1:knots_l,1))!xin
			int2 = sum(fy)!sum(w(1:knots_l,1)*calc_u(a,xx(1:knots_l,1))**2)    
   
			int1 = int1 + int2
		end if
		
	end if
	
	int_u = max(int_ux0 - int1, 0.0_k_p) ! ! = \int^{\inf}_{0} - \int^{x}_{0} ! FIXME

	return 
end function calc_int_u2

!*************************************
! Calculating function \int^{\inf}_{x} U(a,x)^l dx --- see page 91 (Труды астрономической обсерватории XXXI) 
! Calculating it like previous integral ( see calc_int_u ): int^{\inf}_{x} = int^{\inf}_{0} - int^{x}_{0}
! a is Vooight`s parameter
! xin = x is dimensionless frequency [-inf; +inf]
! l -- power
! i -- number if x-knot
!*************************************
elemental function calc_int_ul(a, xin, l, i) result(int_u)
implicit none 
real(k_p), intent(in) :: a, xin
integer, intent(in) :: l, i
real(k_p) :: int_u

	if (l .eq. 1) then
		if (x(i) .ne. xin) then
			int_u = calc_int_u(a, xin)
		else
			int_u = intu1(i)
		end if
	else
		if (l .eq. 2) then
			if (x(i) .ne. xin) then
				int_u = calc_int_u2(a, xin)
			else
				int_u = intu2(i)
			end if
		else
			int_u = 0.0_k_p ! TODO!!!!!!==========================================================FIXME here!!!!!!
		end if
	end if

	return 
end function calc_int_ul

!*************************************
! Calculating function \delta(\beta) --- see page 22, formula 84 (Труды астрономической обсерватории XXXI) 
! b = \beta
! a is Vooight`s parameter
! formula 84 used without out-of-integral term
!*************************************
pure function calc_delta(b, a) result(d)
implicit none 
real(k_p), intent(in) :: b, a
real(k_p) :: d, dd
real(k_p) :: buf, int1
real(k_p), dimension(0:knots) :: frac

	if (b .ne. 0.0_k_p) then
		frac(0:knots) = uax(0:knots)/(b*ua0 + uax(0:knots))
		int1 = simpson_int_knots_fx( knots, knot_x_last, x(0:knots), x_middle, frac(0:knots) ) 
		buf = sqrt(a/(pi*b*ua0))
		dd = int1 + buf*atan(buf)/x(knots)
		d = 2.0_k_p * b * ua0 * dd	
	else ! b = 0
		d = 0.0_k_p			
	end if
 
end function calc_delta


!*************************************
! Calculating function V(u, \beta) --- see page 22, formulas 85-87 (Труды астрономической обсерватории XXXI) 
! u = 1/v
! b = \beta
! 
! formula (95) doesn`t work because it is for Lorentz profile only
! 
!*************************************
elemental subroutine calc_v(u, b, vv, dv, vv_mod) !result(vv)

implicit none 
real(k_p), intent(in)  :: u, b
real(k_p), intent(out) :: vv, dv, vv_mod
real(k_p) :: v

real(k_p), dimension(0:knots) :: frac, y, brackets, y2, ayy, ayy1, arctg
real(k_p) :: buf, int1
integer :: i, knot_brdr

	if (u .ge. 1.0_k_p) then  ! (formula 85) ! series calculates if abs(v)<1
!write(*,*)' (formula 85) 1'		
		v = 1.0_k_p/u
		buf = 0.5_k_p*pi*uasqrt20/sqrt(2.0_k_p)

		y(0:knots) = v*(b*ua0 + uax(0:knots))/ ua0 
		knot_brdr = 0!knots
!write(*,*)' all y=',y(0:knots)		
		do while ((2.0_k_p*y(knot_brdr) .lt. 1.0_k_p) .and. (knot_brdr .lt. knots) ) 
			knot_brdr = knot_brdr + 1
		end do 
			
		y2 = y**2
		ayy = 1.0_k_p!y
		ayy1 = 1.0_k_p
		arctg = 0.0_k_p
		i = 0
!write(*,*)' (formula 85) 2, knot_brdr=',knot_brdr, ' y=',y(knot_brdr), ' y(-1)=',y(knot_brdr-1)
		if (knot_brdr .gt. 0) then		
		do while (maxval(abs(ayy1(0:knot_brdr))) .gt. acc) ! calculating (atan(v*[ ]/ua0))
			arctg(0:knot_brdr) = arctg(0:knot_brdr) + ayy(0:knot_brdr)/(2.0_k_p*i + 1.0_k_p)
			ayy1(0:knot_brdr) = ayy(0:knot_brdr)
			ayy(0:knot_brdr) = -ayy(0:knot_brdr)*y2(0:knot_brdr)
			i = i + 1
!write(*,*)'ayy=',ayy(0),' ',ayy(knot_brdr)
		end do
		arctg = arctg * y
		end if
!write(*,*)' (85). 1'		
		arctg(knot_brdr:knots) = atan2(v*(b*ua0 + uax(knot_brdr:knots)), ua0)	
!write(*,*)' (85). 2'
		frac(0:knots) = arctg(0:knots)*(uax(0:knots))**2	
		int1 = simpson_int_knots_fx( knots, knot_x_last, x(0:knots), x_middle, frac(0:knots) ) 
		vv = (buf - 2.0_k_p*int1 )*v/ua0
		dv = (1.0_k_p - delta_b - vv)!/(v*ua0)
		vv_mod = dv/sqrt(api29*ua0*u) ! в этом случае оно совсем не нужно
!write(*,*)' (85). 3'	
	else ! (formula 86) 		

		y(0:knots) = u*ua0/(b*ua0 + uax(0:knots)) ! (formula 87)
		knot_brdr = 0 ! y(knot_brdr) < border_value (in article border_value = 1/3)
		do while ((2.0_k_p*y(knot_brdr) .lt. 1.0_k_p) .and. (knot_brdr .lt. knots) ) 
			knot_brdr = knot_brdr + 1
		end do

		y2 = y**2
		ayy = 1.0_k_p!y2*y
		ayy1 = 1.0_k_p
		brackets = 0.0_k_p
		i = 1
!write(*,*)'calc_v. knot_brdr=',knot_brdr		
		if (knot_brdr .ne. 0) then
			do while (maxval(abs(ayy1(0:knot_brdr))) .ge. acc) ! calculating (y - atan(y))
				brackets(0:knot_brdr) = brackets(0:knot_brdr) + ayy(0:knot_brdr)/(2.0_k_p*i + 1.0_k_p)
				ayy1(0:knot_brdr) = ayy(0:knot_brdr)
				ayy(0:knot_brdr) = -ayy(0:knot_brdr)*y2(0:knot_brdr)
				i = i + 1
!write(*,*)'calc_v. i=',i,' ayy=',ayy(0),' ',ayy(knot_brdr)
			end do
			brackets(0:knot_brdr) = brackets(0:knot_brdr)*y2(0:knot_brdr)*y(0:knot_brdr)
		end if
		
		brackets(knot_brdr:knots) = (1.0_k_p - atan(y(knot_brdr:knots))/y(knot_brdr:knots))*y(knot_brdr:knots)
		!y(knot_brdr:knots) - atan( y(knot_brdr:knots) )
		frac(0:knots) = brackets(0:knots) * (uax(0:knots)/ua0)**2 
		int1 = simpson_int_knots_fx( knots, knot_x_last, x(0:knots), x_middle, frac(0:knots) ) 	
		
		dv = 2.0_k_p * int1 * ua0/u ! по тетрадке, а не по формуле. Отличие -- в расположении дроби.
		vv = 1.0_k_p - delta_b - dv!*ua0/u		
		vv_mod = dv/sqrt(api29*ua0*u)
	end if !end of formulas (85) and (86)

end subroutine calc_v


!*************************************
! Calculating function lnH(z, l, b) --- see page 23, formulas 88-90 (Труды астрономической обсерватории XXXI) 
! z - argument
! l = \lambda
! b = \beta
! a - argument of U = U(a,x)
!*************************************
elemental function calc_lnH(z, l, b) result(lnH)

implicit none 
real(k_p), intent(in) :: z, l, b
real(k_p) :: lnH 

integer, parameter :: cntv = v_cnt/2

real(k_p) :: s1, int1, int2, ss, q, q4, lnhq, lerchphi, buf1, buf2, buf3, buf4
real(k_p), dimension(0:v_cnt) :: func, log1lvub, lvu, lvu1
real(k_p), dimension(0:knots) :: func2, log1lvxb, lvx, lvx1
integer :: i
integer :: knot_brdr, knot_brdr1
	
if (z .eq. 0.0_k_p) then
	lnH = -0.5_k_p*log(1.0_k_p - l + l*delta_b)
else
	
	if (z .lt. 1.0_k_p) then ! formula (88) why 10 ?
!write(*,*) 'formula (88) z=',z

			s1 = l*uasqrt20*z*log( (1.0_k_p + z**2)/(z**2) )/( 4.0_k_p*sqrt(2.0_k_p)*ua0 )   
			
			lvu = l*vub
			lvu1 = 1.0_k_p!lvu
			log1lvub = lvu1
			i = 1
!write(*,*) 'calc lnH 1'			
			do while (maxval(abs(lvu1)) .ge. acc) ! calculating ln( 1 - l*V(1/v,b) )
				i = i + 1
				lvu1 = lvu1*lvu
				log1lvub = log1lvub + lvu1/(1.0_k_p*i)
			end do
			log1lvub = -log1lvub*lvu
			
			func(0:v_cnt) = ( log1lvub(0:v_cnt) + 0.5_k_p*l*pi*uasqrt20*v(0:v_cnt)/(ua0*sqrt(2.0_k_p)) )/&
																										( v(0:v_cnt)**2 + z**2 )
			int1 = simpson_int(v_cnt, v(0), v(v_cnt), func(0:v_cnt))
			
			if (l .eq. 1.0_k_p) then
				
			!	log1lvub(0:v_cnt) = log( dvv(0:v_cnt) + delta_b)
				log1lvxb(0:knots) = log( dv(0:knots) + delta_b )	

			else		
		
				knot_brdr = 0 ! vxb(knot_brdr) < border_value (in article no border_value)
				do while ((vxb(knot_brdr) .le. 0.5_k_p) .and. (knot_brdr .lt. knots) ) ! FIXME slooooow
					knot_brdr = knot_brdr + 1
				end do
		
				lvx = 0.0_k_p
				lvx(0:knot_brdr)= l*vxb(0:knot_brdr)
				lvx1 = -lvx
				log1lvxb = lvx1
				i = 1
!write(*,*) 'calc lnH 4 knot_brdr=',knot_brdr, ' vxb=',vxb(knot_brdr:knots)		
				do while (maxval(abs(lvx1)) .ge. acc) ! calculating ln( 1 - l*V )
					i = i + 1
					lvx1 = lvx1*lvx
					log1lvxb = log1lvxb + lvx1/(1.0_k_p*i)
				end do
!write(*,*) 'calc lnH 5 knot_brdr=',knot_brdr!, ' vxb=',vxb(knot_brdr:knots)
				knot_brdr1 = knots ! FIXME ============================================== FIXME     костыыыыыыыыыль   !!!!
!				do while ((vxb(knot_brdr1) .eq. 1.0_k_p) .and. (knot_brdr1 .gt. knot_brdr) ) 
!					knot_brdr1 = knot_brdr1 - 1
!				end do
				log1lvxb(knot_brdr:knot_brdr1) = log(1.0_k_p - l*vxb(knot_brdr:knot_brdr1))	
				log1lvxb(knot_brdr1:knots) = log1lvxb(knot_brdr1)
			end if
!write(*,*) 'calc lnH 5 knot_brdr1=',knot_brdr1			
!write(*,*)log1lvxb(knot_brdr:knots)
!read(*,*)
!write(*,*) 'calc lnH 5.1'
			func2(0:knots) = -log1lvxb(0:knots)*dudx(0:knots)/(ua0**2 + (z*uax(0:knots))**2 )
!write(*,*) 'calc lnH 5.2'
			int2 = simpson_int_knots_fx( knots, knot_x_last, x(0:knots), x_middle, func2(0:knots) ) 
!write(*,*) 'calc lnH 6'			
			lnH = s1 - (int1 + ua0*int2)*z/pi
!		end if
	
	else ! z .gt. 1
!write(*,*) 'z>1'

		if ( (l .eq. 1.0_k_p) .and. (b .eq. 0.0_k_p) ) then ! formula (90) ! 
			
!write(*,*) 'formula (90)'!, 'l*vub=',l*vub
			!log1lvub = log(1.0_k_p - l*vub)
			lvu = l*vub
			lvu1 = 1.0_k_p!lvu
			log1lvub = lvu1
			i = 1
!write(*,*) 'calc lnH 1'			
			do while (maxval(abs(lvu1)) .ge. acc) ! calculating ln( 1 - l*V(1/v,b) )
				i = i + 1
				lvu1 = lvu1*lvu
				log1lvub = log1lvub + lvu1/(1.0_k_p*i)
			end do
			log1lvub = -log1lvub*lvu
			func(0:v_cnt) = ( log1lvub(0:v_cnt) )/( 1.0_k_p + (v(0:v_cnt)/z)**2 )	
			int1 = simpson_int(v_cnt, v(0), v(v_cnt), func(0:v_cnt))
!write(*,*) 'calc lnH 2', vxb_mod(0:knots)

!			lvx = 0.0_k_p
!			lvx(0:knot_brdr)= l*vxb(0:knot_brdr)
!			lvx1 = -lvx
!			log1lvxb = 0.0_k_p
!			i = 1	
!			do while (maxval(abs(lvx1)) .ge. acc) ! calculating ln( 1 - l*V )
!				log1lvxb = log1lvxb + lvx1/(1.0_k_p*i)
!				lvx1 = lvx1*lvx
!				i = i + 1
!			end do
!			knot_brdr1 = knots ! FIXME ============================================== FIXME     костыыыыыыыыыль   !!!!
!			do while ((vxb(knot_brdr1) .eq. 1.0_k_p) .and. (knot_brdr1 .gt. knot_brdr) ) 
!				knot_brdr1 = knot_brdr1 - 1
!			end do
			log1lvxb(0:knots) = log(vxb_mod(0:knots))	
!			log1lvxb(knot_brdr:knot_brdr1) = log(vxb_mod(knot_brdr:knot_brdr1))	
!			log1lvxb(knot_brdr1:knots) = log1lvxb(knot_brdr1)

!write(*,*) 'calc lnH 3'
			func2(0:knots) = ( log1lvxb(0:knots) ) * ( -dudx(0:knots)/((z*uax(0:knots))**2 + (ua0)**2) )
			!func2(0:knots) = ( log1lvxb(0:knots) - log(sqrt(api29*uax(0:knots)) ) )& 
			!								* ( -dudx(0:knots)/((z*uax(0:knots))**2 + (ua0)**2) )  ! *uax(0:knots)/uax
!write(*,*) 'calc lnH 4'
			int2 = 	simpson_int_knots_fx( knots, knot_x_last, x(0:knots), x_middle, func2(0:knots) ) 			

!write(*,*) 'calc series. z=', z,' int2=',int2
			! calculating the series in formula (90)
			if (abs(z-1.0_k_p) .lt. acc) then ! честно посчитан ряд
				s1 = 0.915965594177219015054603514932384110774149374281672134266_k_p
			else
				i = 0 
				ss = 1.0_k_p/z
				s1 = ss
				do while (abs(ss/(2.0_k_p*(i+1) + 1.0_k_p)**2) .ge. acc)
					i = i + 1 
					ss = -ss/z**2
					s1 = s1 + ss/(2.0_k_p*i + 1.0_k_p)**2	
!write(*,*)'s1=',s1
				end do
			end if
		
			lnH = -int1/(pi*z) - z*ua0*int2/(pi) - (atan(z)*log(ua0*api29)  - s1)*0.5_k_p/pi + 0.25_k_p*log(z)
		
		else ! formula (89): z>"1"
!write(*,*)'formula (89)'
			s1 = -atan(z)*log(1.0_k_p - l + l*delta_b)/pi
			
			lvu = l*vub
			lvu1 = 1.0_k_p!lvu
			log1lvub = lvu1
			i = 1
!write(*,*) 'calc lnH 1'			
			do while (maxval(abs(lvu1)) .ge. acc) ! calculating ln( 1 - l*V(1/v,b) )
				i = i + 1
				lvu1 = lvu1*lvu
				log1lvub = log1lvub + lvu1/(1.0_k_p*i)
			end do
			log1lvub = -log1lvub*lvu
			
			!func(0:v_cnt) = log(1.0_k_p - l*vub(0:v_cnt))/(1.0_k_p + (v(0:v_cnt)/z)**2 ) !
			func(0:v_cnt) = log1lvub(0:v_cnt)/(1.0_k_p + (v(0:v_cnt)/z)**2 )
			int1 = simpson_int(v_cnt, v(0), v(v_cnt), func(0:v_cnt))
			
			knot_brdr = 0 ! 1 + l/(1-l + l*delta) * (1 - delta_b - vxb(knot_brdr)) > border_value (in article no border_value)
			do while ((l*dv(knot_brdr)/(1.0_k_p - l + l*delta_b) .ge. 0.5_k_p) .and. (knot_brdr .lt. knots) ) 
				knot_brdr = knot_brdr + 1
			end do
	
			lvx = 0.0_k_p
			!lvx(knot_brdr:knots)= l*(1.0_k_p - delta_b - vxb(knot_brdr:knots))/(1.0_k_p - l + l*delta_b)
			lvx(knot_brdr:knots)= l*dv(knot_brdr:knots)/(1.0_k_p - l + l*delta_b) !!!!uax(knot_brdr:knots)*
			lvx1 = lvx ! uax -?
			log1lvxb = lvx1
			i = 1
			
			do while (maxval(abs(lvx1)) .ge. acc) ! calculating ln( 1 + l/(1-l + l*delta) * (1 - delta - vxb) )
				i = i + 1
				lvx1 = -lvx1*lvx
				log1lvxb = log1lvxb + lvx1/(1.0_k_p*i)
			end do
		
			!log1lvxb(0:knot_brdr) = log(1.0_k_p + l*(1.0_k_p - delta_b - vxb(0:knot_brdr))/(1.0_k_p - l + l*delta_b))
			log1lvxb(0:knot_brdr) = log(1.0_k_p + l*dv(0:knot_brdr)/(1.0_k_p - l + l*delta_b)) !!!! uax(0:knot_brdr)*
			func2(0:knots) = log1lvxb(0:knots) * ( -dudx(0:knots)/((uax(0:knots))**2 + (ua0/z)**2) )  ! *uax(0:knots)/uax
			int2 = 	simpson_int_knots_fx( knots, knot_x_last, x(0:knots), x_middle, func2(0:knots) ) 

			lnH = s1 - (int1 + ua0*int2)/(z*pi) !/(pi*z)
		end if
	end if
end if ! z = 0

end function calc_lnH



!*************************************
!Calculating function N_{knl}(m) --- see page 24.
! k,n,l -- some powers
! m - argument
! lambda -- lambda
! b -- beta
! a -- Voight paraneter
!*************************************
elemental function calc_Nknlm(k, n, l, m, lambda, b, a) result(NN)

implicit none
integer, intent(in) :: k, n, l
real(k_p), intent(in) :: m
real(k_p), intent(in) :: lambda, b, a
real(k_p) :: NN

integer, parameter :: knots_l2 = 2*knots_l ! is nesessary? difference between knots_l and 2*knots_2 is ~10^(-7)

real(k_p) :: int1, int2, int3
real(k_p), dimension(0:knots) :: func3, hfunc3, int4 
real(k_p), dimension(1:knots_l2) :: z, w, func2, hfunc2

integer :: i

	int1 = calc_int_ul(a, 0.0_k_p, l+1, 0)/ua0**l
	
	call knots_and_weights(knots_l2, 0.0_k_p, 1.0_k_p, z, w)
!	do i = 1, knots_l
!		hfunc2(i) = exp(calc_lnH(z(i)/(1.0_k_p + b*z(i)), lambda, b))
!	end do
	hfunc2(1:knots_l2) = exp(calc_lnH(z(1:knots_l2)/(1.0_k_p + b*z(1:knots_l2)), lambda, b))
	func2 = hfunc2/((1.0_k_p + b*z + m*z)**k)
	func2 = func2*(z**(k+n-2))/((1.0_k_p + b*z)**n)
	int2 = sum( w * func2 )
	
	!do i = 0, knots
	!	hfunc3(i) = exp(calc_lnH(ua0/(uax(i) + b*ua0), lambda, b))
	!end do
	hfunc3(0:knots) = exp(calc_lnH(ua0/(uax(0:knots) + b*ua0), lambda, b))
!	do i = 0, knots
!		int4(i) = calc_int_ul(a, x(i), l+1, i)
!	end do
	forall (i = 0:knots)
		int4(i) = calc_int_ul(a, x(i), l+1, i)
	end forall

	
	func3(0:knots) = hfunc3(0:knots)/((uax(0:knots) + b*ua0 + m*ua0)**k)
	func3(0:knots) = func3(0:knots)*(-dudx(0:knots)*int4(0:knots))/((uax(0:knots) + b*ua0)**n)
	int3 = simpson_int_knots_fx(knots, knot_x_last, x(0:knots), x_middle, func3(0:knots))
	
	NN = 2.0_k_p*int1*int2 + 2.0_k_p*int3*(ua0**(k+n-l-1))

end function calc_Nknlm

!*************************************
!Calculating moment H_1 of function H
! z -- argument
! l,b -- lambda and beta
! a -- Voight parameter
!*************************************
elemental function calc_H1(z, l, b, a) result(H1)
implicit none
real(k_p), intent(in) :: z, l, b, a
real(k_p) :: H1

	H1 = calc_Nknlm(2, 1, 1, z, l, b, a) !1.0_k_p/

end function calc_H1

!*************************************
!Calculating moment H_0 of function H
! z -- argument
! l,b -- lambda and beta
! a -- Voight parameter
!*************************************
!elemental function calc_H0(z, l, b, a) result(H0)
!implicit none
!real(k_p), intent(in) :: z, l, b, a
!real(k_p) :: H0

!	H0 = calc_Nknlm(1, 1, 1, 0.0_k_p, l, b, a)

!end function calc_H0

!*************************************
!Calculating moment M_0 of function H
! z -- argument
! l,b -- lambda and beta
! a -- Voight parameter
!*************************************
elemental function calc_M0(z, l, b, a) result(M0)
implicit none
real(k_p), intent(in) :: z, l, b, a
real(k_p) :: M0

	M0 = calc_Nknlm(1, 1, 0, 1.0_k_p/z, l, b, a)

end function calc_M0

!*************************************
!Calculating function M(1/m)  --- formula (92)
! q = 1/m -- argument
! l,b -- lambda and beta
! a -- Voight parameter
!*************************************
elemental function calc_M(q, l, b, a) result(MM)
implicit none
integer, parameter :: knots_l2 = knots_l
real(k_p), intent(in) :: q, l, b, a
real(k_p) :: MM, m, int1, int2, hfuncq, h1func
real(k_p), dimension(1:knots_l2) :: z, w, func1, hfunc1, sum1
real(k_p), dimension(0:knots) :: func2, hfunc2, sum2 

integer :: i

	m = 1.0_k_p/q
	hfuncq = exp(calc_lnH(q, l, b))
	h1func = calc_H1(q, l, b, a)
	call knots_and_weights(knots_l2, 0.0_k_p, 1.0_k_p, z, w)
	hfunc1 = exp(calc_lnH(z/(1.0_k_p + b*z), l, b))
!write(*,*)'2', h1func
	where (z(1:knots_l2)*(m - b) .ne. 1.0_k_p)
		sum1 = hfuncq/(1.0_k_p + b*z + m*z) + (hfuncq - hfunc1)/(1.0_k_p + b*z - m*z)
	elsewhere
		sum1 = hfuncq/(1.0_k_p + b*z + m*z) + l*(hfuncq**2)*h1func/(2.0_k_p*z)
	end where	
	func1 = sum1/(1.0_k_p + b*z)
	int1 = sum(w*func1)
	
!write(*,*)'3'	
	hfunc2(0:knots) = hfunc_glob(0:knots)!exp(calc_lnH(ua0/(uax(0:knots) + b*ua0), l, b))
!write(*,*)'3.5',minloc(abs(uax(0:knots)-ua0*(m-b)))
!do i = 0, knots
!	if (abs(uax(i) - ua0*(m - b)) .ge. 1.0e-15 ) then
!		write(*,*)i,')', hfuncq/(uax(i) + b*ua0 + m*ua0),' + ',(hfuncq - hfunc2(i)),'/',(uax(i) - (m - b)*ua0)
!	else
!		write(*,*)i,')',hfuncq/(uax(i) + b*ua0 + m*ua0),' + ',l*(hfuncq**2)*h1func/(2.0_k_p*ua0)
!	end if
!end do
	sum2 = 0.0_k_p
	where (abs(uax(0:knots) - ua0*(m - b)) .ge. 1.0e-15)
		sum2(0:knots) = hfuncq/(uax(0:knots) + b*ua0 + m*ua0) + (hfuncq - hfunc2(0:knots))/(uax(0:knots) + b*ua0 - m*ua0)
	elsewhere
		sum2(0:knots) = hfuncq/(uax(0:knots) + b*ua0 + m*ua0) + l*(hfuncq**2)*h1func/(2.0_k_p*ua0)
	end where
!write(*,*)'4',sum2(0),'!'!,sum2(4:10)!,intu1(0:knots)
	func2(0:knots) = -dudx(0:knots)*intu1(0:knots)/(uax(0:knots) + b*ua0) *sum2(0:knots)
!write(*,*)'5'!, func2(0:knots)
	int2 = simpson_int_knots_fx(knots, knot_x_last, x(0:knots), x_middle, func2(0:knots))
!write(*,*)'q=',q,' int1=',int1,' int2=',int2	

	MM = int1 + 2.0_k_p*ua0*int2

end function calc_M


!*************************************
!Calculating function K_{knl}(m) --- see page 17.
! Calculating is the same as N_{knl}, but H is equal to 1
! k,n,l -- some powers
! m - argument
! b -- beta
! a -- Voight paraneter
!*************************************
!elemental function calc_Kknlm(k, n, l, m, b, a) result(KK)

!implicit none
!integer, intent(in) :: k, n, l
!real(k_p), intent(in) :: m
!real(k_p), intent(in) :: b, a
!real(k_p) :: KK

!integer, parameter :: knots_l2 = 2*knots_l ! is nesessary? difference between knots_l and 2*knots_2 is ~10^(-7)

!real(k_p) :: int1, int2, int3
!real(k_p), dimension(0:knots) :: func3, hfunc3, int4 
!real(k_p), dimension(1:knots_l2) :: z, w, func2, hfunc2
!integer :: i

!	int1 = calc_int_ul(a, 0.0_k_p, l+1, 0)/ua0**l	
!	call knots_and_weights(knots_l2, 0.0_k_p, 1.0_k_p, z, w)

!	hfunc2(1:knots_l2) = 1.0_k_p!exp(calc_lnH(z(1:knots_l2)/(1.0_k_p + b*z(1:knots_l2)), lambda, b))
!	func2 = hfunc2/((1.0_k_p + b*z + m*z)**k)
!	func2 = func2*(z**(k+n-2))/((1.0_k_p + b*z)**n)
!	int2 = sum( w * func2 )	
!	
!	hfunc3(0:knots) = 1.0_k_p!exp(calc_lnH(ua0/(uax(0:knots) + b*ua0), lambda, b))

!	forall (i = 0:knots)
!		int4(i) = calc_int_ul(a, x(i), l+1, i)
!	end forall


!	
!	func3(0:knots) = hfunc3(0:knots)/((uax(0:knots) + b*ua0 + m*ua0)**k)
!	func3(0:knots) = func3(0:knots)*(-dudx(0:knots)*int4(0:knots))/((uax(0:knots) + b*ua0)**n)

!	int3 = simpson_int_knots_fx(knots, knot_x_last, x(0:knots), x_middle, func3(0:knots))
!	
!	KK = 2.0_k_p*int1*int2 + 2.0_k_p*int3*(ua0**(k+n-l-1))

!end function calc_Kknlm


!*************************************
!Calculating function M_*(inf) --- see page 20.
! lambda -- lambda
! b -- beta
! a -- Voight paraneter
!*************************************
!function calc_Mstarinf(lambda, b, a) result(MM)

!implicit none
!real(k_p), intent(in) :: lambda, b, a
!real(k_p) :: MM, hfunc

!	hfunc = exp(calc_lnH(infinity,lambda,b))
!	MM = b*( lambda*calc_Kknlm(1,1,0, 0.0_k_p,b,a)*(hfunc**2)*calc_Nknlm(2,1,1, 0.0_k_p,lambda,b,a) + &
!																					calc_Nknlm(2,1,0, 0.0_k_p,lambda,b,a))

!end function calc_Mstarinf


!*************************************
!Calculating function I(eta1, eta2, XX, l, b) --- formula ?? of paper.
! [eta1,eta2] -- range of angles
! XX -- the edge frequency of F
! l -- lambda
! b - beta
!*************************************
pure function calc_reflected_intensity_line(eta1, eta2, XX, l, b) result(II)

implicit none
integer, parameter :: knots_i = 5
real(k_p), intent(in) :: eta1, eta2, XX, l, b
real(k_p), dimension(0:nmu, 0:knots_all) :: II
integer :: knot_brdr
integer :: i, j, k
real(k_p), dimension(1:knots_i) :: y0, wy0, hfuncy
real(k_p), dimension(0:nmu, 0:knots_all, 0:knots_all) :: intfunc

	knot_brdr = 0
	do while ((x(knot_brdr) .lt. XX) .and. (knot_brdr .le. knots))
		knot_brdr = knot_brdr + 1
	end do

	intfunc = 0.0_k_p
	do k = 0, knot_brdr
		call knots_and_weights(knots_i, (alpha(k) + b)/eta2, (alpha(k) + b)/eta1 , y0, wy0) 
		hfuncy(1:knots_i) = exp( calc_lnH(1.0_k_p/y0(1:knots_i),l,b) )
		forall (j = 0:nmu, i = 0:knots) ! integrating on y0
			intfunc(j, i, k) = sum(wy0(:)*hfuncy(:)/((p(j, i) + y0(:))*y0(:)**2))
		end forall
	end do !k
	
	forall (k = 0:knot_brdr)
		intfunc(:, :, k) = alpha(k)*(alpha(k) + b)*intfunc(:, :, k)
	end forall
	
	! TODO: from x(knot_brdr) to x0
	
	forall (i = 0:knots, j = 0:nmu) !integrating on x0
		II(j,i) = 0.25_k_p*l*ua0*alpha(i)*hfunc(j, i)/mu(j) * &
				  simpson_int_knots_fx(knot_brdr, knot_x_last, x(0:knot_brdr), x_middle, intfunc(j,i, 0:knot_brdr))
	end forall
	
	return

end function calc_reflected_intensity_line



!*************************************
!Calculating function I(eta1, eta2, XX, l, b) --- formula ?? of paper.
! [eta1,eta2] -- range of angles
! XX -- the edge frequency of F
! l -- lambda
! b -- beta
! a -- Voight parameter 
!*************************************
function calc_reflected_intensity_continuum(eta1, eta2, XX, l, b, a) result(II)

implicit none
integer, parameter :: knots_i = 5
real(k_p), intent(in) :: eta1, eta2, XX, l, b, a
real(k_p), dimension(0:nmu, 0:knots_all) :: II
integer :: knot_brdr
integer :: i, j, k
real(k_p), dimension(0:knots_all) :: y1, y2
real(k_p), dimension(1:knots_i) :: y0, wy0, hfuncy, mfuncy
real(k_p), dimension(0:nmu, 0:knots_all, 0:knots_all) :: func, func1, func2, intfunc

	knot_brdr = 0
	do while ((x(knot_brdr) .lt. XX) .and. (knot_brdr .le. knots))
		knot_brdr = knot_brdr + 1
	end do
	
	y2(:) = (alpha(:) + b)/eta2
	y1(:) = (alpha(:) + b)/eta1

	intfunc = 0.0_k_p
	do k = 0, knot_brdr
		call knots_and_weights(knots_i, y2(k), y1(k), y0, wy0) 
		mfuncy(1:knots_i) =  calc_M(1.0_k_p/y0(1:knots_i),l,b,a) 
! ================ TODO: можно вынести M0 за интеграл, и тогда получим уже проинтегрированный (аналитически) результат
		forall(j = 0:nmu, i = 0:knots)! integrating on y0  
			intfunc(j,i,k) = sum( wy0(:)*( &
				(m0func(j,i) + mfuncy(:)) ) &
				/( (p(j,i) + y0(:))*y0(:)**2 ) & ! y0 = y0(k)
				) 
		end forall
	end do !k
	
	func = 0.0_k_p
	forall (k = 0:knot_brdr)
		func(:, :, k) = 0.5_k_p*(1.0_k_p/y2(k)**2 - 1.0_k_p/y1(k)**2)
	end forall
	forall (k = 0:knot_brdr, i = 0:knots, j = 0:nmu, p(j,i) .ne. 0.0_k_p)
		func(j, i, k) = ( p(j, i)*(y1(k) - y2(k)) - y1(k)*y2(k)* &
						log( ( y1(k)*( p(j, i) + y2(k) ) )/( y2(k)*( p(j, i) + y1(k) ) ) ) )/&
						( y1(k)*y2(k)*p(j, i)**2 )  
	end forall
	
	func1 = 0.0_k_p
	func2 = 0.0_k_p
	forall (k = 0:knot_brdr)
		func1(:,:,k) = alpha(k)*(alpha(k) + b)*func(:,:,k)
		func2(:,:,k) = alpha(k)*(alpha(k) + b)*intfunc(:,:,k)
	end forall
	
	! TODO: from x(knot_brdr) to x0
	
	forall (i = 0:knots, j = 0:nmu) ! integrating on x0
		II(j, i) = &
		0.25_k_p*l*ua0/mu(j) * ( &
		simpson_int_knots_fx(knot_brdr, knot_x_last, x(0:knot_brdr), x_middle, func1(j, i, 0:knot_brdr)) + &
		0.125_k_p*(l**2)*ua0*alpha(i)*hfunc(j, i) * &
		simpson_int_knots_fx(knot_brdr, knot_x_last, x(0:knot_brdr), x_middle, func2(j, i, 0:knot_brdr))   &
		)
	end forall
	
	return

end function calc_reflected_intensity_continuum

!*************************************
!*************************************
!*************************************

!subroutine bounds(u_mu)!, time

!implicit none 
!real(k_p), dimension(0:m_mu+1, 0:m_chi+1), intent(inout) :: u_mu

!end subroutine bounds 

end module Funct_module

