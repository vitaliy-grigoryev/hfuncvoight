module Global_module

use Fundamental_data

implicit none
integer, parameter :: knots = 200!180!360!180 ! 1800 count for knots on x
integer, parameter :: nmu = 3 ! knots on angles (mu = cos(angle))
integer, parameter :: knots_all = knots + knots/100 ! real count for knots on x
integer, parameter :: knot_x_last = knots/2 ! last knot of x. After this goes t
integer, parameter :: v_cnt = 50 ! count for knots on u = 1/v (integrating from 0 to 1)
real(k_p), parameter :: x_middle = 10.0_k_p ! value of x_* in paper, x in x_last
real(k_p), parameter :: v_up = 1.0_k_p
integer, parameter :: t_up = 20
integer, parameter :: knot_x_1 = int(knots/x_middle) ! knot x = 1
integer, parameter :: nfact = 170 ! factorials

real(k_p), parameter :: infinity = dble(1.0d100)!1.0_k_p/0.0_k_p

real(k_p), dimension(0:knots_all) :: x, uax, dudx, vxb, intu1, intu2, dv, vxb_mod, hfunc_glob, alpha
real(k_p), dimension(0:nmu) :: mu
real(k_p), dimension(0:nmu, 0:knots_all) :: p, hfunc, m0func
real(k_p) :: uasqrt20, ua0, delta_b, api29
real(k_p), dimension(0:v_cnt) :: vub, v, dvv, vub_mod
integer :: knots_xas
real(k_p), dimension(0:nfact) :: fact

 contains

subroutine initialization_global()
implicit none

integer :: i	

	fact(0) = 1.0_k_p
	do i = 1, nfact
		fact(i) = fact(i-1)*real(i)
	end do
	
	forall (i = 0:nmu)
		mu(i) = cos(0.5*pi*dble(i)/dble(nmu))
	end forall

	forall (i = 0:knot_x_last)
		x(i) = x_middle*dble(i)/dble(knot_x_last)
	end forall
	
	forall (i = 0:knots-knot_x_last)
		x(knot_x_last+i) = x_middle*exp( dble(t_up)*dble(i)/dble(knots - knot_x_last) )
	end forall
	
	forall (i = knots+1:knots_all)
		x(i) = x_middle*exp( dble(t_up)*dble(i-knot_x_last)/dble(knots - knot_x_last) )
	end forall
	
	knots_xas = knots*(1 + 7/t_up)/2 + 1
	
	forall (i = 0:v_cnt)
		v(i) = v_up*dble(i)/dble(v_cnt)
	end forall
	
	uax = 0.0_k_p
	dudx = 0.0_k_p
	intu1 = 0.0_k_p
	intu2 = 0.0_k_p

end subroutine initialization_global

end module Global_module
