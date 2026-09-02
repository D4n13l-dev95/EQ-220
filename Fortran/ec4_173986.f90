!> @file      ec4_173986.f90
!> @author    Daniel Campiotti
!> @date      2026-09-01
!> @brief     Translation of a python code to Fortran for the ec4_173986 project
!
!  Record of revisions:
!  Date           Programmer              Description of change
!  ====           ==========              =====================
!  2026-09-01 -- Daniel Campiotti      Original code

program ec4_173986
    implicit none

    ! --- Constants ---
    integer, parameter :: dp = selected_real_kind(15, 307)

    ! --- Variables ---
    real(dp) :: error, x_new, x_old, y, dy
    real(dp), parameter :: tolerance = 1.0e-7 
    integer :: i ! counter for iterations
    integer :: j, n ! counter for the interval in loop
    real(dp) :: k ! counter for the roots found in the interval
    real(dp), dimension(:), allocatable :: interval ! interval from -1 to 5 with step of 0.1
    real(dp), dimension(:), allocatable :: ls_root ! putting the roots found in the interval

    ! --- Initialization / Input ---

    error = 1.0e0
    x_new = 0.0e0
    y = f(x_new)
    dy = df(x_new)
    i = 0

    do while (error > tolerance)
        x_new = x_new - (y / dy)
        y = f(x_new)
        dy = df(x_new)
        x_old = x_new
        error = abs(x_new - x_old)
        i = i + 1
    end do
    
    print *, "Root found: ", x_new
    print *, "Number of iterations: ", i

    ! running through the interval from -1 to 5 with a step of 0.1
    
    n = 0
    allocate(interval(61))
    do j = 0, 60
        interval(j + 1) = -1.0d0 + j * 0.1d0
    end do

    print *, "Interval: ", interval

    ! second iteration to calculate the roots in the interval
    do k = 1, size(interval) - 1
        if (f(interval(k)) * f(interval(k + 1)) < 0.0d0) then
            call add_root(ls_root, n, interval(k))
            call add_root(ls_root, n, interval(k + 1))
        end if
    end do

    print *, "Roots found in the interval: ", ls_root

    ! finding the roots in the intervals found
    
    do k = 1, size(ls_root)
        x_new = ls_root(k)
        y = f(x_new)    
        dy = df(x_new)
        error = 1.0e0
        i = 0
        do while (error > tolerance)
            x_new = x_new - (y / dy)
            y = f(x_new)    
            dy = df(x_new)
            x_old = x_new
            error = abs(x_new - x_old)
            i = i + 1
        end do

        print *, "Root found: ", x_new
        print *, "Number of iterations: ", i
    
    
    ! functions definitions

    
    contains
      function f(x)
        implicit none
          real(dp), intent(in) :: x
          real(dp), parameter :: a = 5.0d0, b = 6.0d0, c = 1.0d0
          real(dp) :: f
          
          f = x**3 - a * (x*x) + b * x - c
      end function f

      function df(x)  
        implicit none 
          real(dp), intent(in) :: x
          real(dp), parameter :: a = 5.0d0, b = 6.0d0
          real(dp) :: df

          df = 3 * (x*x) - 2 * a * x + b
      end function df

      ! subroutine to add a root to the list of roots found in the interval
      subroutine add_root(ls_root, n, value)
        real(dp), allocatable, intent(inout) :: ls_root(:)
        integer, intent(inout) :: n
        real(dp), intent(in) :: value
        real(dp), allocatable :: tmp(:)

        n = n + 1
        if (.not. allocated(ls_root)) then
            allocate(ls_root(1))
        else if (n > size(ls_root)) then
            allocate(tmp(size(ls_root) * 2))
            tmp(1:size(ls_root)) = ls_root
            call move_alloc(tmp, ls_root)
        end if
        ls_root(n) = value
      end subroutine

end program ec4_173986

! import numpy as np
! import matplotlib.pyplot as plt
! import sympy as sp

! # Definição inicial da função

! def f(x, a = 5, b = 6, c = 1):

!   y = x**3 - a * (x*x) + b * x - c  
  
!   return  y


! def df(x, a = 5, b = 6):
!   dy = 3 * (x**2) - 2 * a * x + b
!   return dy  

! tolerance = 1.0e-7   # definição da tolerância
! i = 0                # contador de iterações

! # initialization of the first guess

! x_new = 0.0e0       
! y = f(x_new)
! dy = df(x_new)
! erro = 1.0e0

! while erro > tolerance:
!     x_new = x_new - y / dy
!     y = f(x_new)
!     dy = df(x_new)
!     x_old = x_new
!     erro = abs(x_new - x_old)
!     i += 1

! print("Root finded: ", x_new)

! # interval from -1 to 5 with step of 0.1

! interval = np.arange(-1, 5, 0.1)
! j = 0  
! ls_root = []  # list to store the roots found in the interval

! # second iteration to calculate the roots in the interval

! for j in range(len(interval) - 1): # evaluating intervals with roots
!     x1 = interval[j]

!     x2 = interval[j + 1] 

!     f1 = f(x1)

!     f2 = f(x2)

!     if f1 * f2 < 0:  # check if the function changes sign in the interval
!         print("Root found in the interval: [", x1, ",", x2, "]")

!         ls_root.append(x1)

!         ls_root.append(x2)
!     else:           # exit the loop if no root is found in the interval
!       continue 

! for k in range(0, len(ls_root), 1):  # finding the roots in the intervals found
!     x_new = ls_root[k]

!     y = f(x_new)

!     dy = df(x_new)

!     erro = 1.0e0

!     while erro > tolerance:
!         x_new = x_new - y/dy

!         y = f(x_new)

!         dy = df(x_new)

!         x_old = x_new

!         erro = abs(x_new - x_old)

!         i += 1

!     print("Root finded: ", x_new)
!     print("Number of iterations: ", i)