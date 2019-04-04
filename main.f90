program Profiles
use Global_module
use Fundamental_data
use Funct_module
use Integrating_module

implicit none

real(k_p) :: a, u, b, d, z, l, lnh, int1, H1, m0, m, h1inf, m0inf, n210b0, mstarinf
real(k_p) :: eta1, eta2, AA, hfunc0

integer, parameter :: nx0 = 1

real(k_p), dimension(1:nx0) :: x0, alpha0

real(k_p), dimension(0:knots_all) :: test_f, test_df
real(k_p) :: test_f1
integer :: i,j,i1,j1, nn, n1, buf_knots,  ix0
integer, dimension(1:nx0) :: knot_brdr


integer, parameter :: knots_i = 5
real(k_p), dimension(0:nmu, 0:knots_all) ::  mfunc, yy!,  m0func ! hfunc,
real(k_p), dimension(0:nmu, 0:knots_all, 0:knots_all) :: func, func1, func2, intfunc
real(k_p), dimension(0:nmu, 0:knots_all, 1:nx0) :: intensity
real(k_p), dimension(1:nx0) :: z0, Mz0, H0, y0
real(k_p), dimension(0:knots_all) :: zz, bufz, y1, y2

real(k_p), dimension(0:nmu,0:knots_all) :: zzmu

real(k_p), dimension(1:knots_i) :: hfuncy, mfuncy

character(2) :: calc_flag
	

	open(unit = 100, file='files_created.txt')
	call initialization_global()
!=== TESTS =====
!	do i = 0, knots_all
!		test_f(i) = test_func(x(i))
!	end do
!	
!	test_df = calc_df_dx(knots_all, x, test_f)
!	
!	open(unit = 10, file = "temp.txt")
!	do i = 0, knots_all
!		write(10,*) i, test_df(i) , (-2.0_k_p*x(i)/(x(i)**2 + 1.0_k_p)**2)
!	end do
!	close(10)
!	
!	nn = knot_x_last
!	n1 = 0
!	write(*,*) simpson_int(nn-n1, x(n1), x(nn), test_f(n1:nn)), log(x(nn) + 1.0_k_p) , log(x(n1) + 1.0_k_p)
!	
!	
!	write(*,*) simpson_int_knots_fx(knots, knot_x_last, x(0:knots), x_middle, test_f(0:knots)) - &
!		log(x(knots) + 1.0_k_p) + log(x(0) + 1.0_k_p)
!	write(*,*)"Test finished. Press Ctrl+C"
!	read(*,*)	
!=== TESTS =====	

	
	write(*,*)"Our goal is to calculate reflected intensity (formula 42) for line or continuum"
	write(*,*)"It depends on x1 and eta --- the frequency and the cos(falling angle)"
	write(*,*)"Also we need a and b --- parameters of Voight function"
	write(*,*)"Input all parameters in order: a, b, x0, l, eta0"
	open(unit = 11, file='in.txt')
	read(11,*) a
	read(11,*) b
	read(11,*) l
	read(11,*) x0
	read(11,*) eta1
	read(11,*) eta2
	read(11,*) calc_flag
	close(11)
	
	l = 1.0_k_p - l
	write(*,*)"So, a = ", a
	write(*,*)"    b = ", b
	write(*,*)"    l = ", l
	write(*,*)"   x1 = ", x0
	write(*,*)"  eta = [", eta2,":",eta1,"] degree"
	write(*,*)"Calculating: [",calc_flag,"]"
	eta1 = cos(eta1*pi/180.0_k_p)
	eta2 = cos(eta2*pi/180.0_k_p)
	
	write(*,*)"cos eta = [", eta2,":",eta1,"] degree"
	

	write(*,*)"Start initialization "
	call initialization() ! initialize knots, weights in Integrating_module
	call initialization_global()
	ua0 = calc_u(a, 0.0_k_p)
	uasqrt20 = calc_u(a*sqrt(2.0_k_p), 0.0_k_p)
	api29 = 2.0_k_p*a*pi/9.0_k_p
	write(*,*)"Finish initialization "

	uax(0:knots_all) = calc_u(a, x(0:knots_all))
	write(*,*)"U(a,x) calculated"
	dudx(0:knots_all) = calc_df_dx(knots_all, x(0:knots_all), uax(0:knots_all))
	write(*,*)"dU(a,x)/dx calculated"

	
	AA = ua0!/intu1(0)
	alpha = uax/ua0

!	write(*,*)"Calculating Delta(b,a)"
	delta_b = calc_delta(b, a)
	write(*,*)"Delta(b,a) = ", delta_b

	write(*,*)"Calculating V(u, b) for 1/v"
	vub(0) = 0.0_k_p ! formula (85) -- calculating V(inf, b)
	dvv(0) = 1.0_k_p - delta_b - vub(0)
	vub_mod(0) = 0.0_k_p
	call calc_v(1.0_k_p/v(1:v_cnt), b, vub(1:v_cnt), dvv(1:v_cnt), vub_mod(1:v_cnt)) ! case u>1

	write(*,*)"Calculating V(u, b) for u"
	call calc_v(uax(0:knots)/ua0, b, vxb(0:knots), dv(0:knots), vxb_mod(0:knots)) ! case u<1
	
	write(*,*)"V(u, b) calculated"

	write(*,*)"Calculating lnH(z). (=>   z = mu*U(a,0)/(U(a,x) + b*U(a,0)))"
	zz = ua0/(uax + b*ua0) ! zz  == 1/p
	forall (i = 0:knots_all, j = 0:nmu)
		zzmu(j,i) = zz(i)*mu(j)
	end forall
	p = 1.0_k_p/zzmu
	hfunc(0:nmu,0:knots_all) = exp(calc_lnH(zzmu(0:nmu,0:knots_all), l, b))
	hfunc_glob(0:knots_all) = hfunc(0,0:knots_all)

	if (calc_flag .eq. 'mb' .or. calc_flag .eq. 'ml' .or. calc_flag .eq. 'bl') then
		write(*,*)'Calculating reflected intensity for monochrome line, F = 1, cos(angle of incidence)=',eta2
		intensity = 0.0_k_p
		forall (ix0 = 1:nx0)
			y0(ix0) = ( calc_u(a,x0(ix0))/ua0 + b)/eta2
		end forall
		H0(1:nx0) = exp(calc_lnH(1.0_k_p/y0(1:nx0), l, b))
		
		forall (i = 0:knots, j = 0:nmu, ix0 = 1:nx0)
			intensity(j,i,ix0) = 0.25_k_p*l*ua0*alpha(i)*(hfunc(j, i)/mu(j)) * H0(ix0)/(p(j,i) + y0(ix0))
		end forall
		do ix0 = 1, nx0
			call data_output('ml',intensity(0:nmu,0:knots,ix0), x0(ix0), nint(180*acos(eta2)/pi))
		end do
		write(*,*)'>>>> Data for reflected intensity for monochrome line has been written to ',nx0, ' files <<<<'
	end if
	
	if (calc_flag .eq. 'la' ) then
		write(*,*)'Calculating reflected intensity for line, monoangle, Q = 1'
		intensity = 0.0_k_p
		do ix0 = 1, nx0
			write(*,*) 'Integrating to x=',x0(ix0)!,' knot_brdr=',knot_brdr(ix0)
			intensity(:, :, ix0) = calc_reflected_intensity_line_angle(eta1, x0(ix0), l, b)
		end do ! ix0
	
	
		do ix0 = 1, nx0
			call data_output('la',intensity(0:nmu,0:knots,ix0), x0(ix0), nint(180*(acos(eta1) - acos(eta2))/pi))
		enddo ! ix0
		write(*,*)'>>>> Data for reflected intensity in line has been written to ',nx0, ' files <<<<'
	end if
	
	if (calc_flag .eq. 'bb' .or. calc_flag .eq. 'll' .or. calc_flag .eq. 'bl') then
		write(*,*)'Calculating reflected intensity for line, F = 1'
		intensity = 0.0_k_p
		do ix0 = 1, nx0
			write(*,*) 'Integrating to x=',x0(ix0)!,' knot_brdr=',knot_brdr(ix0)
			intensity(:, :, ix0) = calc_reflected_intensity_line_full(eta1, eta2, x0(ix0), l, b)
		end do ! ix0
		
		do ix0 = 1, nx0
			call data_output('ll',intensity(0:nmu,0:knots,ix0), x0(ix0), nint(180*(acos(eta1) - acos(eta2))/pi))
		enddo ! ix0
		write(*,*)'>>>> Data for reflected intensity in line has been written to ',nx0, ' files <<<<'
		
	end if ! if line
	
	if (calc_flag .eq. 'oc' ) then
		write(*,*)'Calculating own radiation, continuum sources '
		intensity = 0.0_k_p
		
		call get_u_integrals()
		
		write(*,*)'Calculating M0(z=1/m)'
		m0func(0:nmu, 0:knots_all) = calc_M0(zzmu(0:nmu,0:knots_all), l, b, a)
		call store_data_h()
		
		do ix0 = 1, nx0
			write(*,*)'q = alpha(x0) + beta = ',calc_u(a, x0(ix0))/ua0 + b, " x0=",x0(ix0)
			intensity(:,:,ix0) = calc_own_intensity_continuum(calc_u(a, x0(ix0))/ua0 + b, l, b, a)
		end do ! ix0
		
		do ix0 = 1, nx0
			call data_output("oc", intensity(0:nmu,0:knots,ix0), x0(ix0), nint(180*acos(eta2)/pi))
		end do
		write(*,*)'>>>> Data for reflected intensity for monoangle continuum has been written to ',nx0, ' files <<<<'
	end if
	
	if (calc_flag .eq. 'mb' .or. calc_flag .eq. 'mc' .or. calc_flag .eq. 'bc') then
		write(*,*)'Calculating reflected intensity for monoangle continuum, Q = 1, cos(angle of incidence)=',eta2
		intensity = 0.0_k_p
		
		call get_u_integrals()
		
		write(*,*)'Calculating M0(z=1/m)'
		m0func(0:nmu, 0:knots_all) = calc_M0(zzmu(0:nmu,0:knots_all), l, b, a)
		call store_data_h()
		
		do ix0 = 1, nx0
			write(*,*)'Integrating from 0 to x=',x0(ix0)!,' knot_brdr=',knot_brdr(ix0)
			intensity(:,:,ix0) = calc_reflected_intensity_continuum(eta2, x0(ix0), l, b, a)
		end do ! ix0
		
		do ix0 = 1, nx0
			call data_output('mc', intensity(0:nmu,0:knots,ix0), x0(ix0), nint(180*acos(eta2)/pi))
		end do
		write(*,*)'>>>> Data for reflected intensity for monoangle continuum has been written to ',nx0, ' files <<<<'
	end if
	
	if (calc_flag .eq. 'bb' .or. calc_flag .eq. 'cc' .or. calc_flag .eq. 'bc') then
		write(*,*)'Calculating intensity for continuum sources, F=1'! (see page 28, formula I.2 + modification)'
		
		call get_u_integrals()
		
		write(*,*)'Calculating M0(z=1/m)' ! zzmu = ua0*mu/(uax + b*ua0)
		m0func(0:nmu, 0:knots_all) = calc_M0(zzmu(0:nmu,0:knots_all), l, b, a)
		call store_data_h()
		
		do ix0 = 1, nx0
			write(*,*)'Integrating from 0 to x=',x0(ix0)!,' knot_brdr=',knot_brdr(ix0)
			intensity(:,:,ix0) = calc_reflected_intensity_continuum_full(eta1, eta2, x0(ix0), l, b, a)
		end do ! ix0
	
		do ix0 = 1, nx0
			call data_output('cc',intensity(0:nmu,0:knots,ix0), x0(ix0), nint(180*acos(eta1)/pi - 180*acos(eta2)/pi))
		enddo ! ix0
		write(*,*)'>>>> Data for continuum sources intensity has been written to ',nx0, ' files <<<<'
	end if ! continuum

	call store_data_h()

!	===== For Alexander Mushtukov ==============================================================
! 	write(buf_str,'("out/hp_",es6.1e1,"_",es8.1e2,"_",es6.1e1,".txt")') a, 1.0_k_p-l, b
! 	open(unit = 100, file = buf_str)
!	write(100,'(" # lambda = ",e21.14,"    beta = ",e21.14)') l, b
!	write(100,*) ' ' 
!	do i = knots, 0, -1
! 		write(100,"(i12,2e27.16)")knots-i+1, 1.0_k_p/zz(i), hfunc(0, i)
!	end do
! 	close(100)
!	write(*,*)'All Data of hfunc has been written to  ', buf_str
!	===== For Alexander Mushtukov ==============================================================	

	write(*,*)'Complete!' !See file "int<i>.txt"
	close(100) ! files_created.txt
	

contains


subroutine data_output(flag, II, xx, eta)
implicit none
character(2), intent(in) :: flag 						! calc_flag
real(k_p), dimension(0:nmu, 0:knots), intent(in) :: II 	! array of intensity
real(k_p), intent(in) :: xx 							! border frequency
integer, intent(in) :: eta 								! incident angle [deg]
character(45) :: buf_str

	write(buf_str,'("out/",a,"_",es6.1e1,"_",es6.1e1,"_",es6.1e1,"_",es6.1e1,"_",i0".txt")') flag, a, b, 1.0_k_p-l, xx, eta
	open(file = buf_str, unit = 15) 
	write(15,'("# a = ",e13.7,"; beta = b = ",e13.7,"; 1-lambda = ",e10.3,"; delta=",e13.7)') a, b, 1.0_k_p-l, delta_b
	write(15,'("# x0 =",e15.7)') xx
	do i = 0, knots 	! mu: 0 = 0 deg, 1 = 15 deg, 2 = 30 deg, 3 = 45 deg, 4 = 60 deg, 5 = 75 deg, 6 = 90 deg
		write(15,'(50e14.6)')x(i), II(0:nmu, i)/II(0:nmu,0)
	end do
	close(15)
	
	write(100,'(a)')buf_str ! files created
	
	open(file = "out.txt", unit=16) ! for comfortable using
	write(16,'("# a = ",e13.7,"; beta = b = ",e13.7,"; 1-lambda = ",e10.3,"; delta=",e13.7)') a, b, 1.0_k_p-l, delta_b
	write(16,'("# x0 =",e15.7)') xx
	do i = 0, knots 	! mu: 0 = 0 deg, 1 = 15 deg, 2 = 30 deg, 3 = 45 deg, 4 = 60 deg, 5 = 75 deg, 6 = 90 deg
		write(16,'(50e14.6)')x(i), II(0:nmu, i)/II(0:nmu,0)
	end do
	close(16)
	

end subroutine data_output

subroutine get_u_integrals()
implicit none
integer :: ioflag, rknots, i
real(k_p) :: ra, rb, rl, rd
real(k_p), dimension(1:10) :: buf
logical :: needcalc ! Do I need to calculate integrals?

	needcalc = .true.
!	open(unit=99, file="data_h.txt", status="old", iostat=ioflag)
	intu1 = 0.0_k_p
	intu2 = 0.0_k_p
! 	if (ioflag .eq. 0) then ! file exists
! 		write(*,*)'Try reading integrals of U(a,x)'
! 		read(99,'(12x,i10)') rknots
! 		read(99,'(6x,g13.7,13x,g13.7,13x,g10.3,8x,g13.7)') ra, rb, rl, rd
! 		read(99,*) ! x, uax...
! 		if (rknots .eq. knots_all .and. ra .eq. a ) then !.and. rb .eq. b .and. 1.0_k_p-rl .eq. l
! 			read(99,"(9e27.16)")buf(1), buf(2), buf(3), buf(4), buf(5), buf(6), intu1(0), intu2(0), buf(7)
! 			if (intu1(0) .ne. 0.0_k_p) then
! 				write(*,*)'Parameter a is correct, integrals are not empty'
! 				do i = 1, knots
! 					read(99,"(9e27.16)")buf(1), buf(2), buf(3), buf(4), buf(5), buf(6), intu1(i), intu2(i), buf(7)
! 				end do
! 				write(*,*)'Integrals have been succesfully read'
! 				needcalc = .false.
! 			end if
! 			close(99)
! 		end if
! 	end if
	
	if (needcalc) then
		write(*,*)'Calculating integrals of U(a,x)'
! 		do i = 1, knots
! 		write(*,*) i
! 			intu1(i) = table_int(i, knots, x(0:knots), uax(0:knots))
! 			intu2(i) = table_int(i, knots, x(0:knots), uax(0:knots)**2)
! 		end do
		intu1(0:knots) = calc_int_u(a, x(0:knots))
		write(*,*)"Int U(a,x) dx calculated"
		intu2(0:knots) = calc_int_u2(a, x(0:knots))
		write(*,*)"Int U(a,x)^2 dx calculated"

	end if

end subroutine get_u_integrals

subroutine store_data_h()
implicit none
integer :: i
	! H-function data output
	open(unit = 10, file = "data_h.txt", status="replace")
	write(10,'("# knots_all=",i0)') knots_all
	write(10,'("# a = ",e13.7,"; beta = b = ",e13.7,"; 1-lambda = ",e10.3,"; delta=",e13.7)') a, b, 1.0_k_p-l, delta_b
	write(10,'("#",13x,"x",25x,"U(a,x)",19x,"dU(a,x)/dx",16x,"V(x, b)",11x,"z=U(a,0)/(U(a,x)+b*U(a,0))",9x,"H(z)")') !"\int U(a,x)dx",16x,
	do i = 0, knots_all
 		write(10,"(10e27.16)")x(i), uax(i), dudx(i), vxb(i), zz(i), hfunc(0, i), intu1(i), intu2(i), m0func(0, i)
	end do
 	close(10)
 	write(*,*)'>>>> All Data of H-function has been written to file "data_h.txt" <<<<'

end subroutine store_data_h

end program Profiles
