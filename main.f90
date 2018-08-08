program Profiles
use Global_module
use Fundamental_data
use Funct_module
use Integrating_module

implicit none

real(k_p) :: a, u, b, d, z, l, lnh, int1, H1, m0, m, h1inf, m0inf, n210b0, mstarinf
real(k_p) :: eta1, eta2, AA, hfunc0

integer, parameter :: nx0 = 4

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


logical :: dataread ! read data from file?
 character(45) :: buf_str
 character(2) :: calc_flag
	

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
	
	

	dataread = .false.
	
	write(*,*)"Our goal is to calculate reflected intensity (formula 42) for line or continuum"
	write(*,*)"It depends on x1 and eta --- the frequency and the cos(falling angle)"
	write(*,*)"Also we need a and b --- parameters of Voight function"
	write(*,*)"Input all parameters in order: a, b, x0, l, eta0"
	open(unit = 11, file='in.txt')
	read(11,*) a
	read(11,*) b
	read(11,*) x0
	read(11,*) l
	read(11,*) eta1
	read(11,*) eta2
	read(11,*) calc_flag
	close(11)
	
	write(*,*)"So, a = ", a
	write(*,*)"    b = ", b
	write(*,*)"   x1 = ", x0
	write(*,*)"    l = ", l
	write(*,*)"  eta = [", eta2,":",eta1,"] degree"
	write(*,*)"Calculating: [",calc_flag,"] (b - both, l - line, c - continuum)"
	eta1 = cos(eta1*pi/180.0_k_p)
	eta2 = cos(eta2*pi/180.0_k_p)
	

	write(*,*)"Start initialization "
	call initialization() ! initialize knots, weights in Integrating_module
	call initialization_global()
	ua0 = calc_u(a, 0.0_k_p)
	uasqrt20 = calc_u(a*sqrt(2.0_k_p), 0.0_k_p)
	api29 = 2.0_k_p*a*pi/9.0_k_p
	write(*,*)"Finish initialization "

	if (dataread) then
	
		open(unit = 12, file = "data_u.txt")
		read(12,'(a2,i5, 4e27.16)') buf_str, buf_knots, a, b, x0, l
		if (buf_knots .ne. knots_all) then
 			write(*,*) "Counts of knots aren`t equal! Aborting..."
 			call exit(0)
 		end if
		do i = 0, knots_all
 			read(12,"(7e27.16)")x(i), uax(i), dudx(i), intu1(i), intu2(i), hfunc(0, i), mfunc(0, i)
		end do
 		close(12)
 		write(*,*)"From file values had been read:"
 		write(*,*)"    a = ", a
		write(*,*)"    b = ", b
		write(*,*)"    x = ", x0
		write(*,*)"    l = ", l, ' <- we use them.'
		write(*,*)'All x(i), uax(i), dudx(i) intu1(i), intu2(i)'		
	
	else

! 		if (a .gt. 1.0_k_p) then
! 			write(*,*)'Attention! a > 1, so x = x/a (turn to Lorentz`s halfwidth)'
! 			x = x/a ! переход к частоте, считаемой в лоренцовских ширинах
! 		end if
		
!		write(*,*)"Start calculating U(a,x)"
		uax(0:knots_all) = calc_u(a, x(0:knots_all))
!		do i = 0, knots_all
!write(*,*)'U(a,x) i=',i, ' x=',x(i)
!			uax(i) = calc_u(a, x(i))
!		end do
!read(*,*)
!=== TESTS =====
		!uax = test_f
!=== TESTS =====
		write(*,*)"U(a,x) calculated"
		dudx(0:knots_all) = calc_df_dx(knots_all, x(0:knots_all), uax(0:knots_all))
		write(*,*)"dU(a,x)/dx calculated"
				
	end if
	
	AA = ua0!/intu1(0)
	alpha = uax/ua0

!	write(*,*)"Calculating Delta(b,a)"
	delta_b = calc_delta(b, a)
	write(*,*)"Delta(b,a) = ", delta_b

	write(*,*)"Calculating V(u, b) for 1/v"
	!read(*,*) b
	vub(0) = 0.0_k_p ! formula (85) -- calculating V(inf, b)
	dvv(0) = 1.0_k_p - delta_b - vub(0)
	vub_mod(0) = 0.0_k_p
	call calc_v(1.0_k_p/v(1:v_cnt), b, vub(1:v_cnt), dvv(1:v_cnt), vub_mod(1:v_cnt)) ! case u>1
	
! 	do i = 1, v_cnt
! 		call calc_v(1.0_k_p/v(i), b, vub(i), dvv(i),vub_mod(i))
! write(*,*)"V.",i," calc_v(1/v,b) = ",vub(i)!, " dv=",dvv(i)
! 	end do

	write(*,*)"Calculating V(u, b) for u"
	call calc_v(uax(0:knots)/ua0, b, vxb(0:knots), dv(0:knots), vxb_mod(0:knots)) ! case u<1
	
! 	do i = 0, knots
! 		call calc_v(uax(i)/ua0, b, vxb(i), dv(i),vxb_mod(i))
! write(*,*)"V.",i," calc_v(u,b) = ",vxb(i),' x(i)=',x(i)!, " dv=",dv(i)
!	end do
	
	write(*,*)"V(u,b) calculated"
!write(*,*)"V(1/v,b) = ", vub
!write(*,*)"V(u,b)   =", vxb(0:10)," !!! ", vxb(knots-10:knots)


! open(unit = 10, file = "t.txt")
! do i = 0, v_cnt
! 	write(10,"(4e27.16)")v(i), 0.0, 0.0, vub(i)
! end do
! do i = 0, knots
! 	write(10,"(5e27.16)")x(i), uax(i), uax(i)/ua0 , dudx(i), vxb(i)
! end do
! close(10)

!----------------------------	
	write(*,*)"Calculating lnH(z). (=>   z = mu*U(a,0)/(U(a,x) + b*U(a,0)))"
	mu(0) = 1.0_k_p
	mu(1) = sqrt(3.0_k_p)/2.0_k_p
	mu(2) = 0.75_k_p
	mu(3) = sqrt(2.0_k_p)/2.0_k_p
	mu(4) = 0.5_k_p
	mu(5) = 0.258819_k_p
	zz = ua0/(uax + b*ua0) ! zz  == 1/p
!!	do i = 0, knots_all
!!		do j = 0, nmu
!!write(*,*)i,j, ' zz(i)=',zz(i),' mu(j)=',mu(j), ' x(i)=',x(i)
	forall (i = 0:knots_all, j = 0:nmu)
		zzmu(j,i) = zz(i)*mu(j)
	end forall
	p = 1.0_k_p/zzmu
!		end do
!	end do
!write(*,*)zz
	hfunc(0:nmu,0:knots_all) = exp(calc_lnH(zzmu(0:nmu,0:knots_all), l, b))
	hfunc_glob(0:knots_all) = hfunc(0,0:knots_all)


!	do i = 0,knots_all
!		do j = 0, nmu
!!write(*,*)i,j, ' zz(i)=',zz(i),' mu(j)=',mu(j), ' x(i)=',x(i)
!			hfunc(j,i) = exp(calc_lnH(zzmu(j,i), l, b))
!		end do
!	end do

		
	forall (i = 0:knots, j = 0:nmu)
		yy(j,i) = 1.0_k_p/zzmu(j,i)!(alpha(i) + b)/mu(j) ! = 1/zzmu
	end forall
	
	if (calc_flag .eq. 'mb' .or. calc_flag .eq. 'ml') then
		write(*,*)'Calculating reflected intensity for monochrome line, F = 1, cos(angle of incidence)=',eta2
		intensity = 0.0_k_p
		forall (ix0 = 1:nx0)
			y0(ix0) = ( calc_u(a,x0(ix0))/ua0 + b)/eta2
		end forall
		H0(1:nx0) = exp(calc_lnH(1.0_k_p/y0(1:nx0), l, b))
		
		forall (i = 0:knots, j = 0:nmu, ix0 = 1:nx0)
			intensity(j,i,ix0) = 0.25_k_p*l*ua0*alpha(i)*(hfunc(j, i)/mu(j)) * H0(ix0)/(yy(j,i) + y0(ix0))
		end forall
		do ix0 = 1, nx0
			call data_output('ml',intensity(0:nmu,0:knots,ix0), ix0)
		end do
		write(*,*)'>>>> Data for reflected intensity for monochrome line has been written to ',nx0, ' files <<<<'
	end if
	
		
		
	if (calc_flag .eq. 'bb' .or. calc_flag .eq. 'll') then
		write(*,*)'Calculating reflected intensity for line, F = 1'
		intensity = 0.0_k_p
		do ix0 = 1, nx0
			write(*,*) 'Integrating to x=',x0(ix0)!,' knot_brdr=',knot_brdr(ix0)
			intensity(:, :, ix0) = calc_reflected_intensity_line(eta1, eta2, x0(ix0), l, b)
		end do ! ix0
		
		! data output
		do ix0 = 1, nx0
			!write(buf_str,'("ri",i0,"-",f4.1,".txt")')ix0,a	
			!write(buf_str,'("out/l",i0,"_",es6.1e1,"_",es8.1e2,"_",es6.1e1,".txt")') ix0,a, 1.0_k_p-l, b
			write(buf_str,'("out/l",i0,".txt")') ix0
			open(file = buf_str, unit = 15+ix0) 
			write(15+ix0,'("# a = ",e13.7,"; beta = b = ",e13.7,"; 1-lambda = ",e10.3,"; delta=",e13.7)') a,b,1.0_k_p-l,delta_b
			write(15+ix0,'("# x_border=",e15.7)') x0(ix0)
			do i = 0, knots 					! mu: 0 = 0 deg, 1 = 30 deg, 2 = 41.5 deg, 3 = 45 deg, 4 = 60 deg, 5 = 75 deg
				write(15+ix0,'(50e14.6)')x(i), intensity(0:nmu, i, ix0)
			end do
			close(15+ix0)
		enddo ! ix0
		write(*,*)'>>>> Data for reflected intensity in line has been written to ',nx0, ' files <<<<'
		
	end if ! if line
	
	
	
	
	if (calc_flag .eq. 'bb' .or. calc_flag .eq. 'cc') then
		write(*,*)'Calculating intensity for continuum sources, F=1'! (see page 28, formula I.2 + modification)'
		
		write(*,*)'Calculating integrals of U(a,x)'
		intu1(0:knots) = calc_int_u(a, x(0:knots))
		write(*,*)"Int U(a,x) dx calculated"
		intu2(0:knots) = calc_int_u2(a, x(0:knots))
		write(*,*)"Int U(a,x)^2 dx calculated"
		
		write(*,*)'Calculating M0(z=1/m)' ! zzmu = ua0*mu/(uax + b*ua0)
		m0func(0:nmu, 0:knots_all) = calc_M0(zzmu(0:nmu,0:knots_all), l, b, a)
		
! 		intensity = 0.0_k_p
! 		! common multiplier
! 		forall (i = 0:knots, j = 0:nmu, ix0 = 1:nx0) 
! 			intensity(j,i,ix0) = 0.25_k_p*l*AA/mu(j)
! 		end forall
		
! 		write(*,*)'Calculating knot borders'
! 		knot_brdr = 0
! 		do ix0 = 1, nx0
! 			do while ((x(knot_brdr(ix0)) .lt. x0(ix0)) .and. (knot_brdr(ix0) .le. knots))
! 				knot_brdr(ix0) = knot_brdr(ix0) + 1
! 			end do
! 		end do
		
! 		y2(:) = (alpha(:) + b)/eta2
! 		y1(:) = (alpha(:) + b)/eta1
		
		do ix0 = 1, nx0
			write(*,*)'Integrating from 0 to x=',x0(ix0)!,' knot_brdr=',knot_brdr(ix0)
			intensity(:,:,ix0) = calc_reflected_intensity_continuum(eta1, eta2, x0(ix0), l, b, a)
			!write(*,*)'Integrating on y0, part 1'
! 			do i = 0, knot_brdr(ix0)
! 				call knots_and_weights(knots_i, y2(i), y1(i) , y0, wy0) 
! 				mfuncy(1:knots_i) =  calc_M(1.0_k_p/y0(1:knots_i),l,b,a) 
! ! ================ TODO: можно вынести M0 за интеграл, и тогда получим уже проинтегрированный (аналитически) результат
! 				forall(j1 = 0:nmu, i1 = 0:knots )! integrating on y0  
! 					intfunc(j1,i1,i) = sum( wy0(:)*( &
! 						(m0func(j1,i1) + mfuncy(:)) ) &
! 						/( (yy(j1,i1) + y0(:))*y0(:)**2 ) &
! 						) 
! 				end forall
! 			end do

			!write(*,*)'Integrating on y0, part 2'
! 			func = 0.0_k_p
! 			!lim y->0
! 			forall (i = 0:knot_brdr(ix0), j1=0:nmu, i1 = 0:knots)
! 				func(j1,i1,i) = 0.5_k_p*(1.0_k_p/y2(i)**2 - 1.0_k_p/y1(i)**2 )
! 			end forall
! 			forall (i = 0:knot_brdr(ix0), i1 = 0:knots, j1 = 0:nmu, yy(j1,i1) .ne. 0.0_k_p)
! 				func(j1,i1,i) = ( yy(j1,i1)*(y1(i) - y2(i)) - y1(i)*y2(i)* &
! 								log( ( y1(i)*( yy(j1,i1) + y2(i) ) )/( y2(i)*( yy(j1,i1) + y1(i) ) ) ) )/&
! 								( y1(i)*y2(i)*yy(j1,i1)**2 )  
! 			end forall
			
! 			func1 = 0.0_k_p
! 			func2 = 0.0_k_p
! 			forall (i = 0:knot_brdr(ix0))
! 				func1(:,:,i) = alpha(i)*(alpha(i) + b)*func(:,:,i)
! 				func2(:,:,i) = alpha(i)*(alpha(i) + b)*intfunc(:,:,i)
! 			end forall
		
! 			forall (i = 0:knots, j = 0:nmu) ! integrating on x0
! 				intensity(j,i, ix0) = &
! 					intensity(j,i, ix0) * ( &
! 					simpson_int_knots_fx(knot_brdr(ix0), knot_x_last, x(0:knot_brdr(ix0)), x_middle, func1(j,i, 0:knot_brdr(ix0))) +&
! 					0.125_k_p*(l**2)*AA*alpha(i)*hfunc(j,i)* &
! 					simpson_int_knots_fx(knot_brdr(ix0), knot_x_last, x(0:knot_brdr(ix0)), x_middle, func2(j,i, 0:knot_brdr(ix0))) &
! 					)
! 			end forall
		end do ! ix0
	
		! data output
		do ix0 = 1, nx0
			!write(buf_str,'("c",i0,"-",f4.1,".txt")')ix0,a
			!write(buf_str,'("out/c",i0,"_",es6.1e1,"_",es8.1e2,"_",es6.1e1,".txt")') ix0, a, 1.0_k_p-l, b
			write(buf_str,'("out/c",i0,".txt")') ix0
		
			open(file = buf_str, unit = 15+ix0) ! 0 = 0 deg, 1 = 30 deg, 2 = 41.5 deg, 3 = 45 deg, 4 = 60 deg, 5 = 75 deg
			write(15+ix0,'("# a = ",e13.7,"; beta = b = ",e13.7,"; 1-lambda = ",e10.3,"; delta=",e13.7)') a,b,1.0_k_p-l,delta_b
			write(15+ix0,'("# ")')
			do i = 0, knots
				write(15+ix0,'(50e14.6)')x(i), intensity(0:nmu, i, ix0)
			end do
			close(15+ix0)
	
		enddo ! ix0
		write(*,*)'>>>> Data for continuum sources intensity has been written to ',nx0, ' files <<<<'
	end if ! continuum

	
	
! H-function data output
	open(unit = 10, file = "data_h.txt")
	write(10,'("# knots_all=",i0)') knots_all
	write(10,'("# a = ",e13.7,"; beta = b = ",e13.7,"; 1-lambda = ",e10.3,"; delta=",e13.7)') a, b, 1.0_k_p-l, delta_b
	write(10,'("#",13x,"x",25x,"U(a,x)",19x,"dU(a,x)/dx",16x,"V(x, b)",11x,"z=U(a,0)/(U(a,x)+b*U(a,0))",9x,"H(z)")') !"\int U(a,x)dx",16x,
	do i = 0, knots_all
 		write(10,"(10e27.16)")x(i), uax(i), dudx(i), vxb(i), zz(i), hfunc(0, i), intu1(i), intu2(i), m0func(0, i)
	end do
 	close(10)
 	write(*,*)'>>>> All Data of H-function has been written to file "data_h.txt" <<<<'

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
	

contains

elemental function test_func(x) result(f)
real(k_p), intent(in) :: x
real(k_p) :: f
	
	f = 1.0_k_p/(x**2 + 1.0_k_p)

end function test_func

subroutine data_output(flag, II, par)
implicit none
character(2), intent(in) :: flag
real(k_p), dimension(0:nmu, 0:knots), intent(in) :: II
integer, intent(in) :: par

	write(buf_str,'("out/",a,i0,"_",es6.1e1,"_",es6.1e1,"_",es6.1e1,".txt")') flag, par, a, b, 1.0_k_p-l
	open(file = buf_str, unit = 15) 
	write(15,'("# a = ",e13.7,"; beta = b = ",e13.7,"; 1-lambda = ",e10.3,"; delta=",e13.7)') a,b,1.0_k_p-l,delta_b
	write(15,'("# x_border=",e15.7)') x0(par)
	do i = 0, knots 					! mu: 0 = 0 deg, 1 = 30 deg, 2 = 41.5 deg, 3 = 45 deg, 4 = 60 deg, 5 = 75 deg
		write(15,'(50e14.6)')x(i), II(0:nmu, i)
	end do
	close(15)

end subroutine data_output

end program Profiles
