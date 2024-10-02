module inverse_mod

  implicit none

  private

  public :: matinv

contains
  
  FUNCTION matinv(mmat)

    double precision, DIMENSION(:,:) :: mmat
    double precision, DIMENSION(size(mmat,1),size(mmat,2)) :: matinv
    double precision, dimension(size(mmat,1),size(mmat,2)) :: aa, cc

    !
    aa = mmat
    cc = 0
    call inverse(aa , cc,size(mmat,1))
    matinv = cc
    !

  contains

    subroutine inverse(a,c,nnn)
      !============================================================
      ! Inverse matrix
      ! Method: Based on Doolittle LU factorization for Ax=b
      ! Alex G. December 2009
      !-----------------------------------------------------------
      ! input ...
      ! a(n,n) - array of coefficients for matrix A
      ! n      - dimension
      ! output ...
      ! c(n,n) - inverse matrix of A
      ! comments ...
      ! the original matrix a(n,n) will be destroyed
      ! during the calculation
      !===========================================================

      implicit none
      integer :: nnn
      double precision :: a(nnn,nnn), c(nnn,nnn)
      double precision :: L(nnn,nnn), U(nnn,nnn), b(nnn), d(nnn), x(nnn)
      double precision :: coeff
      integer :: i, j, k

      ! step 0: initialization for matrices L and U and b
      ! Fortran 90/95 aloows such operations on matrices

      L=0.0
      U=0.0
      b=0.0

      ! step 1: forward elimination
      do k=1, nnn-1
         do i=k+1,nnn
            coeff=a(i,k)/a(k,k)
            L(i,k) = coeff
            do j=k+1,nnn
               a(i,j) = a(i,j)-coeff*a(k,j)
            end do
         end do
      end do

      ! Step 2: prepare L and U matrices
      ! L matrix is a matrix of the elimination coefficient
      ! + the diagonal elements are 1.0
      do i=1,nnn
         L(i,i) = 1.0
      end do

      ! U matrix is the upper triangular part of A
      do j=1,nnn
         do i=1,j
            U(i,j) = a(i,j)
         end do
      end do

      ! Step 3: compute columns of the inverse matrix C
      do k=1,nnn
         b(k)=1.0
         d(1) = b(1)
         ! Step 3a: Solve Ld=b using the forward substitution
         do i=2,nnn
            d(i)=b(i)
            do j=1,i-1
               d(i) = d(i) - L(i,j)*d(j)
            end do
         end do
         ! Step 3b: Solve Ux=d using the back substitution
         x(nnn)=d(nnn)/U(nnn,nnn)
         do i = nnn-1,1,-1
            x(i) = d(i)
            do j=nnn,i+1,-1
               x(i)=x(i)-U(i,j)*x(j)
            end do
            x(i) = x(i)/u(i,i)
         end do
         ! Step 3c: fill the solutions x(n) into column k of C
         do i=1,nnn
            c(i,k) = x(i)
         end do
         b(k)=0.0
      end do
    end subroutine inverse

  END FUNCTION matinv

end module inverse_mod

!==================================================================================

module bspline

  implicit none
  private
  public :: bspln

contains

  recursive function bspln(t,n,i,k,x) result(b)

    implicit none

    integer :: n,i,k
    double precision :: t(:),x,b
    double precision :: c1, c2

    if(k==1)then
       if(x/=t(n+k))then
          if((x>=t(i).and.x<t(i+1)))then
             b=1.d0
          else
             b=0.d0
          endif
       else
          if((x>t(i).and.x<=t(i+1)))then
             b=1.d0
          else
             b=0.d0
          endif
       endif
    else
       if((x-t(i))==0.d0.or.(t(i+k-1)-t(i))==0.d0)then
          c1=0.d0
       else
          c1 = (x-t(i))/(t(i+k-1)-t(i))
       endif
       if((t(i+k)-x)==0.d0.or.(t(i+k)-t(i+1))==0.d0)then
          c2 = 0.d0
       else
          c2 = (t(i+k)-x)/(t(i+k)-t(i+1))
       endif
       b=c1*bspln(t,n,i,k-1,x)+c2*bspln(t,n,i+1,k-1,x)
    endif

  end function bspln
end module bspline

!==================================================================================

module bspline_compression

  use inverse_mod
  use bspline
  
  implicit none

  private

  public :: get_coeffs_spline
  ! public :: t, bb

contains

  ! subroutine get_coeffs_spline(t_length, dt, n_spline, poly_order, t, bb)
  subroutine get_coeffs_spline(t_length, dt, inc_basis, n_spline, poly_order, t, bb)

    implicit none
    
    integer, intent(in)             :: poly_order
    integer, intent(out)            :: n_spline
    double precision, intent(in)    :: dt, inc_basis
    double precision, intent(inout) :: t_length ! t_length in seconds
    double precision, allocatable   :: t(:), bb(:) ! Allocatable ?? Like, really?

    ! local
    ! double precision, allocatable :: dir(:)
    double precision, allocatable :: A(:,:), At(:,:), invAtA(:,:)
    double precision, allocatable :: time(:), signal(:)
    integer          :: nt
    integer          :: i, j

    ! double precision, external :: bspln
    ! double precision, external :: matinv

    ! Determining the number of time iterations
    nt = floor(dble(t_length) / dt)

    if (mod(nt, 2) .ne. 1)then
       nt = nt + 1;
    endif

    !    nt = int(nt)
    ! Add something if nt isn't an odd integer

    allocate(time(nt))
    allocate(signal(nt))
    
    time(nt/2 + 1) = 0

    do i = nt/2 + 2,nt
       time(i) = time(i-1) + dt
    enddo

    do i = nt/2,1,-1
       time(i) = time(i+1) - dt
    enddo

    t_length = time(nt) - time(1)

    ! Add a print statement
    n_spline = ceiling((t_length) / inc_basis) + poly_order

    signal(:) = 0
    signal(nt/2+1) = 1

 
    allocate (t(n_spline + poly_order + 1))
    call Get_Knot_Vector(n_spline,poly_order+1, t, 2)

    t = t*(t_length) + time(1)

    allocate(A(nt,n_spline))

    do i = 1,n_spline
       do j = 1,nt
          A(j,i) = bspln(t, n_spline, i, poly_order+1, time(j))
       enddo
    enddo

    allocate(At(n_spline, nt))

    At = transpose(A)

    ! take the inverse of the product At*A
    
    allocate(invAtA(n_spline, n_spline))

    invAtA = matinv(matmul(At, A))

    ! now inverse the system A*b = d (in the least square sense since A isn't a square matrix)
    ! here d is the original signal and we search for the b-spline coefficients b
    
    allocate(bb(n_spline))

    bb = matmul(invAtA, matmul(At, signal))

    deallocate(A)
    deallocate(At)
    deallocate(invAtA)
    deallocate(time)
    deallocate(signal)

    
  contains
  
    subroutine Get_Knot_Vector(nn, k, t, opt)
      !-----------------------------------    
      !                                                                                                                    
      implicit none
      !                                                                                                                     
      integer :: nn,k
      double precision :: t(nn+k)
      !
      integer :: i,opt
      double precision :: dx
      !
      dx = 1.d0 / (nn + 1 - k)

      !
      do i = k, nn+1
         t(i) = (i-k)*dx
      enddo
      !
      select case(opt)
      case(1)
         !
         ! not-a-knot at both ends
         !
         t(1  :k-1) = t(k  )
         t(nn+2:nn+k) = t(nn+1)
         !
      case(2)
         !
         ! uniform knot
         !
         do i = k-1,1,-1
            t(i) = t(i+1)-dx
         enddo
         !
         do i = nn+2,nn+k
            t(i) = t(i-1)+dx
         enddo
         !
      case default
         !
      end select
      !
      return
      !
      !-----------------------------
    end subroutine Get_Knot_Vector
    !-----------------------------   

  end subroutine get_coeffs_spline
  
end module bspline_compression
