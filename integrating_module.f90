module Integrating_module

use Fundamental_data

implicit none
integer, parameter, private :: max_dim = 28
real(k_p), parameter, private :: acc = dble(1.0e-150) !?????

real(k_p), dimension(1:max_dim, 1:max_dim), private :: roots, weights, table_c ! saved  Legendre`s roots and their weights; weights of table indefinite integration

contains

subroutine initialization()
implicit none
integer :: i

	roots = 0.0_k_p
	weights = 0.0_k_p
	
	do i = 1, max_dim
		 call knots_and_weights_init(i, roots(i,1:i), weights(i, 1:i))
	end do
	
	table_c = 0.0_k_p
	do i = 1, max_dim
		 call table_weights_init(i, table_c(i, 1:i))
	end do

end subroutine initialization

!*************************************
! Returns values of Legendre`s roots and weights for integrating (init - for initialization, normal - for computation)
! n - length of l
! a, b - [a,b]
! l - vector of roots
! aa - vector of weights
! 
! For true integrating find the sum of multipliers:  
!    do i = 1, n {ii=ii+aa(i)*f(x(i))} end do  ....or.... sum(aa(1:n)*f(x(1:n)))
!*************************************	
pure subroutine knots_and_weights_init(n, l, aa) 
implicit none

integer, intent(in) :: n
integer :: i
real(k_p), dimension(1:n), intent(out) :: l, aa !l - roots of legendre`s polynomial, k - result vector in system, aa - A
real(k_p), dimension(1:n) :: k
real(k_p), dimension(1:n, 1:n) :: t ! SIGMA(A[i]*t[i]^k)=2/(k+1) or 0

	!maybe some optimization: save and read files with roots and weights?
	l = legendre(n)
	
	t(1:n, 1) = 1.0_k_p
	k(1) = 2.0_k_p
	do i = 2, n
		t(1:n, i) = t(1:n, i-1)*l(1:n)
		k(i) = dble(mod(i, 2))*2.0_k_p/dble(i)
	enddo
	aa = CHOOSE(t, k, n)
	
end subroutine knots_and_weights_init

pure subroutine knots_and_weights(n, a, b, l, aa) 
implicit none

integer, intent(in) :: n
real(k_p), intent(in) :: a, b 
real(k_p), dimension(1:n), intent(out) :: l, aa !l - roots of legendre`s polynomial, k - result vector in system, aa - A
	
	l(1:n) = roots(n, 1:n)
	aa(1:n) = weights(n, 1:n)
	
	l(1:n) = a + (1.0_k_p + l(1:n))*(b - a)/2.0_k_p ! scaling roots to [a,b]
	aa = aa * (b - a)/2.0_k_p ! scale weights
	
end subroutine knots_and_weights


!*************************************
! Returns values of C_m weights for table indefinite integrating (init - for initialization, normal - for computation)
! n - length of c
! c - vector of weights
! 
! For true integrating find the sum of multipliers:  
!    do i = 1, n {ii=ii+cs(i)*f(x(i))} end do  ....or.... sum(c(1:n)*f(x(1:n)))
!*************************************	
pure subroutine table_weights_init(n, c) 
implicit none

integer, intent(in) :: n
integer :: i, j
real(k_p), dimension(1:n), intent(out) :: c ! result vector
real(k_p), dimension(1:n) :: b
real(k_p), dimension(1:n, 1:n) :: a 

	
	a(1:n, 1) = 1.0_k_p ! first row
	b(1) = 0.5_k_p
	forall (i = 2:n)
		b(i) = 1.0_k_p/dble(2.0_k_p*i - 1)
	end forall
	
	forall (i = 2:n, j = 1:n)
		a(j,i) = (j-1)**(2*(i-1)) + j**(2*(i-1))
	end forall
	
	! Solve A*c = b
	c = CHOOSE(a, b, n)
	
end subroutine table_weights_init



!*************************************
! Returns integral f(x)dx on some segment, calculated by Simpson`s method
! n - count of parts of segments of [a, b]
! a, b - segment [a,b]
! f - vector of f(x)
!*************************************
pure function simpson_int(n, a, b, f) result(ii)
implicit none
integer, intent(in) :: n
real(k_p), intent(in) :: a, b
real(k_p), dimension(0:n), intent(in) :: f
real(k_p) :: ii
	
	ii = 0.0_k_p
	if (n .ne. 0) then
		ii = (b - a)/dble(3*n)  
		ii = ii*( f(0) + 2.0_k_p*sum(f(2:n-2:2)) + 4.0_k_p*sum(f(1:n-1:2)) + f(n) ) 
	end if
	
end function simpson_int

!*************************************
! Returns integral f(x)dx on some segment, calculated by Simpson`s method, but using values of argument in all knots
! n - count of parts of segments of [a, b]
! m - position of x_middle
! x - vector of x in segment [a,b]
! f - vector of f(x)
!*************************************
pure function simpson_int_knots_fx(n, m, x, f) result(ii)
implicit none
integer, intent(in) :: n, m ! full count of elements in vector x, start of t
real(k_p), dimension(0:n), intent(in) :: x, f
!real(k_p), intent(in) :: x_m ! x = x_m*exp(t)
real(k_p) :: ii, int1, int2
	
	int1 = 0.0_k_p
	int2 = 0.0_k_p
	if (n .gt. m) then
		int1 = simpson_int(m, x(0), x(m), f(0:m))
		int2 = simpson_int(n-m, 0.0_k_p, log(x(n)/x(m)), x(m:n)*f(m:n))
	else
		int1 = simpson_int(n, x(0), x(n), f(0:n))
	end if
	ii = int1 + int2
	
end function simpson_int_knots_fx

!*************************************
! Returns integral f(x)dx on some segment, calculated by table indefinite integrating method
! n - length of vectors x and f
! i - number of element in given vector x, where integration is needed, 0 < i <=n
! x - vector of x
! f - vector of f = f(x)
! Calculation integral of \int^{x_i}_{x_{i-1}} f(x) dx
!*************************************
pure function table_knots_int(n, i, x, f) result(ii)
implicit none
integer, intent(in) :: n, i
real(k_p), dimension(0:n-1), intent(in) :: x, f
real(k_p) :: ii
integer :: k
real(k_p), dimension(0:n-1) :: c ! result vector
real(k_p), dimension(0:n-1) :: b
real(k_p), dimension(0:n-1, 0:n-1) :: a 

	forall (k = 0:n-1)
		b(k) = (x(i)**(k+1) - x(i-1)**(k+1))/dble(k+1)
	end forall
	
	forall (k = 0:n-1)
		a(0:n-1,k) = x(0:n-1)**k
	end forall
	
	! Solve A*c = b
	c = CHOOSE(a, b, n)
	
	ii = sum(c(:)*f(:))
	
end function table_knots_int

!*************************************
! Returns integral f(x)dx on some segment, calculated by table indefinite integrating method
! n - number of knot to integrate: \int_0^x[n]
! knots - legth of x and f
! x - vector of x in segment [a,b]
! f - vector of f(x)
! Yes, this method could be optimized very well, but it is not nesessary
!*************************************
pure function table_int(n, knots, x, f) result(integral)
implicit none
integer, intent(in) :: n, knots
real(k_p), dimension(0:knots), intent(in) :: x, f
real(k_p) :: integral

integer, parameter :: order = 18
real(k_p), dimension(0:knots) :: ii
integer :: i
	
	ii = 0.0_k_p
	!if (n .gt. 0) then !.and. n .le. order
		forall (i = 1:order/2)
			ii(i) = table_knots_int(order+1, i, x(0:order), f(0:order))
		end forall
	!end if
	if (n+order/2-1 .ge. knots-order/2) then
		forall (i = order/2+1:knots-order/2-1)
			ii(i) = table_knots_int(order+1, order/2, x(i-order/2:i+order/2), f(i-order/2:i+order/2))
		end forall
		forall (i = knots-order/2:knots)
			ii(i) = table_knots_int(order+1, i, x(knots-order:knots), f(knots-order:knots))
		end forall
	else
 		forall (i = order/2+1:n+order/2+1)
			ii(i) = table_knots_int(order+1, order/2, x(i-order/2:i+order/2), f(i-order/2:i+order/2))
		end forall
	end if

	integral = sum(ii(0:n))
	
end function table_int

!*************************************
! Returns integral f(x)dx on some segment, calculated by table indefinite integrating method, but using values of argument in all knots
! n - count of parts of segments of [a, b]
! order - order of integrating
! m - position of x_middle
! x - vector of x (-order+1:order)
! f - vector of f(x)
!*************************************
! pure function table_int_knots_fx(n, order, m,  x, f) result(ii)
! implicit none
! integer, intent(in) :: n, order, m ! full count of elements in vector x, order if integrating
! real(k_p), dimension(-order+1:order), intent(in) :: x, f
! !real(k_p), intent(in) :: x_m ! x = x_m*exp(t)
! real(k_p) :: ii
! 	
! 	ii = 0.0_k_p
! 	if (n .lt. order) then
! 		! TODO: extrapolating <-
! 		ii = 0.0_k_p
! 	else
! 		if (n .gt. m-order) then
! 			if (n .le. m) then
! 				! TODO: extrapolating ->
! 				ii = 0.0_k_p!table_int(m, order, x(0:m), f(0:m)) 
! 			else
! 				if (n .lt. m+order) then
! 				! TODO: extrapolating <-
! 				ii = 0.0_k_p!table_int(n-m, order, log(x(m:n+order)/x(m)), x(m:n+order)*f(m:n+order))
! 				else ! n >= m+order
! 					ii = table_int(order, x(n-order+1:n+order)*f(n-order+1:n+order))
! 				end if
! 			end if
! 		else
! 			ii = table_int(order, f(m-order+1:n+order))
! 		end if
! 	end if
! 	
! end function table_int_knots_fx

!*************************************
! Returns integral f(x)dx on some segment, calculated by trapezoid method
! n - count of parts of segments of [a, b]
! x - vector of x
! f - vector of f(x)
!*************************************
pure function trapezoid_int(n, x, f) result(ii)
implicit none
integer, intent(in) :: n
real(k_p), dimension(0:n), intent(in) :: x, f
real(k_p) :: ii

	ii = 0.5_k_p*sum( (f(1:n) + f(0:n-1)) * (x(1:n) - x(0:n-1)) )

end function trapezoid_int


!*************************************
! Finds max(abs(x0)), where x0 - root of polynomial with a - koefficients
!*************************************
pure function bernoulli(n, a) result (x) 
implicit none
integer, intent(in) :: n
real(k_p), dimension(0:n), intent(in) :: a
real(k_p) :: x
integer :: i
real(k_p), dimension(0:n) :: y
real(k_p) :: x1, buf
	
		y = 1.0_k_p
		x = 0.0_k_p
		x1 = 1.0_k_p
		do while (abs(x-x1) .gt. acc) 
			y(1:n)=y(0:n-1)
			x1=x
			buf=0.0_k_p
			do i=1,n 
				buf=buf+a(i)*y(i)
			enddo
			y(0)=-buf/a(0)
			x=y(0)/y(1)
		enddo

end function bernoulli


!*************************************
! Realizes Gorner`s scheme: divide one polynomial on another one
! a0x^n + a1x^(n-1) + ... divide (x-x0)
!*************************************	
pure function gorner(n,a,x0) result (aa)
implicit none
integer, intent(in) :: n
real(k_p), dimension(0:n),intent(in) :: a
real(k_p), intent(in) :: x0
integer :: i
real(k_p), dimension(0:n) :: aa
	
	aa(0)=a(0)
	do i=1,n
		aa(i)=a(i)+x0*aa(i-1)
	enddo
	
end function gorner

!*************************************
! Returns all roots of given polynom
! a0x^n + a1x^(n-1) + ... + an
!*************************************	
pure function allroots(n, aa) result (x)
implicit none
integer, intent(in) :: n
integer :: i, nn
real(k_p), dimension(0:n), intent(in) :: aa
real(k_p), dimension(0:n) :: a
real(k_p), dimension(1:n) :: x
	
	nn = n
	a = aa
	do i = 1, n
		x(i) = bernoulli(nn, a)
		nn=nn-1
		a(0:nn) = gorner(nn+1, a, x(i))              
	enddo
	
end function allroots


!*************************************
! Returns roots of Legendre`s polynom degree = n
! a0x^n + a1x^(n-1) + ... + an
!*************************************	
pure function legendre(n) result (l) 
implicit none
integer, intent(in) :: n
integer :: i,j,nn
real(k_p), dimension (0:n) :: p0,p1,p2
real(k_p), dimension (1:n) :: l
	
	l = 0.0_k_p
	p0 = 0.0_k_p
	p1 = 0.0_k_p
	p2 = 0.0_k_p !calculate Legendre`s polynomial without coefficients equivalent 0
	p0(0) = 1.0_k_p
	p1(0) = 1.0_k_p
	nn = n
	do j = 2, nn
		p2 = p1*(2*j-1)/j
		p2(1:nn) = p2(1:nn) - p0(0:nn-1)*(j-1)/j
		p0 = p1
		p1 = p2
	enddo
				
	l(1:n/2) = allroots(n/2,p2) ! now we have half of count of roots, but squared
	do i = 1, n/2
		l(i) = -sqrt(l(i))    ! now we have half of roots
	enddo
		
	do i = 1 + mod(n,2) + n/2, n
		l(i) = -l(n-i+1)
	enddo !reverse roots. now we have all roots, increased
		
	! полином Лежандра имеет сииметричные корни. Коэффициенты - либо при четных степенях, либо
	! при нечетных. Т.о. надо найти квадраты корней, решая полином относительно x^2
	
end function legendre

!*************************************
! Solves system of OLE Ax = b
! aa - matrix A [n*n]
! b - vector b
! n - length of b
!*************************************	
pure function CHOOSE(aa, b, n) result(x) ! copied from module4 and changed a few
implicit none
integer :: i, j, k
integer, intent(in) :: n
real(k_p) :: q
real(k_p), dimension(1:n,1:n), intent(in) :: aa
real(k_p), dimension(1:n+1,1:n) :: a
real(k_p), dimension(1:n), intent(in) :: b
real(k_p), dimension(1:n) :: x, stl
real(k_p), dimension(1:n+1) :: str
integer, dimension(1:n) :: chstl
integer, dimension(1:n) :: chstr
integer, dimension(1:2) :: xy
integer :: buf
	
	
	a(1:n,1:n) = aa(1:n,1:n)
	a(n+1,1:n) = b(1:n)	
	
	forall (i = 1:n)
		chstr(i) = i
		chstl(i) = i
	end forall
	
	do i = 1, n-1
		xy = maxloc(abs(a(i:n,i:n)))
		xy(1) = xy(1)+i-1
		xy(2) = xy(2)+i-1
		str = a(:,xy(2))
		a(:,xy(2)) = a(:,i)
		a(:,i) = str
		buf = chstr(i)
		chstr(i) = chstr(xy(2))
		chstr(xy(2)) = buf
			
		stl = a(xy(1),:)
		a(xy(1),:) = a(i,:)
		a(i,:) = stl
		buf = chstl(i)
		chstl(i) = chstl(xy(1))
		chstl(xy(1)) = buf
		q = a(i,i)
		a(:,i) = a(:,i)/q
		
		!dividing string by element
		do k = i+1, n
			q = a(i,k)
			a(:,k) = a(:,k) - a(:,i)*q
		enddo
		
	enddo
	
	a(n+1,n) = a(n+1,n)/a(n,n)
	a(n,n)= 1.0_k_p
	
	do i = n-1, 1, -1
   		do j = n, i+1, -1
			a(n+1,i) = a(n+1,i) - a(j,i)*a(n+1,j)
		enddo
	enddo
	
	do i = 1, n
		x(chstl(i)) = a(n+1,i)
	enddo! change to another way
	
end function CHOOSE
	

end module Integrating_module
