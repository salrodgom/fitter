! Michael F. Hutt
! http://www.mikehutt.com
! Mar. 31, 1998
! $Id: frosen.f90,v 1.4 2007/07/10 12:45:32 mike Exp $
!
! This program will attempt to minimize Rosenbrock's function using the 
! Nelder-Mead simplex method. The program was originally developed in C. 
! To be consistent with the way arrays are handled in C, all arrays will 
! start from 0.
!
! to compile this program with g77 use:
! g77 -ffree-form -o frosen.exe frosen.f

! * Copyright (c) 1998-2004 <Michael F. Hutt>
! *
! * Permission is hereby granted, free of charge, to any person obtaining
! * a copy of this software and associated documentation files (the
! * "Software"), to deal in the Software without restriction, including
! * without limitation the rights to use, copy, modify, merge, publish,
! * distribute, sublicense, and/or sell copies of the Software, and to
! * permit persons to whom the Software is furnished to do so, subject to
! * the following conditions:
! *
! * The above copyright notice and this permission notice shall be
! * included in all copies or substantial portions of the Software.
! *
! * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
! * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
! * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
! * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
! * LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
! * OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
! * WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
!
!
! Converted program from Fortran 77 to Fortran 90
! This program will attempt to minimize Rosenbrock's function using the 
! Nelder-Mead simplex method. The program was originally developed in C. 
! To be consistent with the way arrays are handled in C, all arrays will
! start from 0.
! compiles with ELF90
! ======================================================================
! Start of main program
module mod_simplex
 !implicit none
 !real :: sp(0:1) 
 !integer :: d=2, iprint
 !real :: e=1.0e-8, scale=1.0
 !sp(0) = -1.2
 !sp(1) = 1.0
 !iprint = 0 ! set to 0 to get step by step print out
 private
 public  :: simplex
 ! call simplex(sp,d,e,scale,iprint)
CONTAINS
! ======================================================================
! This is the function to be minimized

real function func(x) result(rosen) 
  implicit none
  real, intent (in) :: x(:)
  rosen = (100*((x(2)-x(1)**2)**2)+(1.0-x(1))**2)

  return 
end function func

! ======================================================================
! This is the simplex routine

subroutine simplex(start, n, EPSILON, scale, iprint)
  implicit none

  integer, intent (in) :: n, iprint
  real, intent (in), dimension(0:n-1) :: start
  real, intent (in) :: EPSILON, scale

! Define Constants
  integer, parameter :: MAX_IT = 1000
  real, parameter :: ALPHA=1.0
  real, parameter :: BETA=0.5
  real, parameter :: GAMMA=2.0

! ======================================================================
! Variable Definitions
! 
! Integer vs = vertex with the smallest value
! Integer vh = vertex with next smallest value 
! Integer vg = vertex with largest value 
! Integer i,j,m,row
! Integer k = track the number of function evaluations 
! Integer itr = track the number of iterations
! real v = holds vertices of simplex 
! real pn,qn = values used to create initial simplex 
! real f = value of function at each vertex 
! real fr = value of function at reflection point 
! real fe = value of function at expansion point 
! real fc = value of function at contraction point 
! real vr = reflection - coordinates 
! real ve = expansion - coordinates 
! real vc = contraction - coordinates 
! real vm = centroid - coordinates 
! real min
! real fsum,favg,s,cent
! real vtmp = temporary array passed to FUNC
! ======================================================================

  Integer :: vs,vh,vg
  Integer :: i,j,k,itr,m,row
  real, dimension(:,:), allocatable :: v
  real, dimension(:), allocatable  :: f
  real, dimension(:), allocatable :: vr
  real, dimension(:), allocatable :: ve
  real, dimension(:), allocatable :: vc
  real, dimension(:), allocatable :: vm
  real, dimension(:), allocatable :: vtmp
  real :: pn,qn
  real :: fr,fe,fc
  real :: min,fsum,favg,cent,s

  allocate (v(0:n,0:n-1))
  allocate (f(0:n))
  allocate (vr(0:n-1))
  allocate (ve(0:n-1))
  allocate (vc(0:n-1))
  allocate (vm(0:n-1))
  allocate (vtmp(0:n-1))

! create the initial simplex
! assume one of the vertices is 0.0

  pn = scale*(sqrt(n+1.)-1.+n)/(n*sqrt(2.))
  qn = scale*(sqrt(n+1.)-1.)/(n*sqrt(2.))

  DO i=0,n-1
    v(0,i) = start(i)
  END DO

  DO i=1,n
    DO j=0,n-1
      IF (i-1 == j) THEN
        v(i,j) = pn + start(j)
      ELSE
        v(i,j) = qn + start(j)
      END IF
    END DO
  END DO


! find the initial function values

  DO j=0,n
! put coordinates into single dimension array
! to pass it to FUNC
    DO m=0,n-1
      vtmp(m) = v(j,m)
    END DO
    f(j) = FUNC(vtmp)
  END DO

! Print out the initial simplex
! Print out the initial function values

  IF (iprint == 0) THEN
    Write(*,*) "Initial Values"
    Write(*,300) ((v(i,j),j=0,n-1),f(i),i=0,n)
  END IF

  k = n+1

! begin main loop of the minimization

DO itr=1,MAX_IT
! find the index of the largest value
  vg = 0
  DO j=0,n
    IF (f(j) .GT. f(vg)) THEN
      vg = j
    END IF
  END DO

! find the index of the smallest value
  vs = 0
  DO j=0,n
    If (f(j) .LT. f(vs)) Then
      vs = j
    END IF
  END DO

! find the index of the second largest value
  vh = vs
  Do j=0,n
    If ((f(j) .GT. f(vh)) .AND. (f(j) .LT. f(vg))) Then
      vh = j
    END IF
  END DO

! calculate the centroid
  DO j=0,n-1
  cent = 0.0
    DO m=0,n
      If (m .NE. vg) Then
        cent = cent + v(m,j)
      END IF
    END DO
    vm(j) = cent/n
  END DO

! reflect vg to new vertex vr
  DO j=0,n-1
    vr(j) = (1+ALPHA)*vm(j) - ALPHA*v(vg,j)
  END DO
  fr = FUNC(vr)
  k = k+1

  If ((fr .LE. f(vh)) .AND. (fr .GT. f(vs))) Then
    DO j=0,n-1
      v(vg,j) = vr(j)
    END DO
    f(vg) = fr
  END IF

! investigate a step further in this direction
  If (fr .LE. f(vs)) Then
    DO j=0,n-1
      ve(j) = GAMMA*vr(j) + (1-GAMMA)*vm(j)
    END DO
    fe = FUNC(ve)
    k = k+1

! by making fe < fr as opposed to fe < f(vs), Rosenbrocks function
! takes 62 iterations as opposed to 64. 

    If (fe .LT. fr) Then
      DO j=0,n-1
        v(vg,j) = ve(j)
      END DO
      f(vg) = fe
    Else
      DO j=0,n-1
        v(vg,j) = vr(j)
      END DO
      f(vg) = fr
    END IF
  END IF

! check to see if a contraction is necessary
  If (fr .GT. f(vh)) Then
    DO j=0,n-1
      vc(j) = BETA*v(vg,j) + (1-BETA)*vm(j)
    END DO
    fc = FUNC(vc)
    k = k+1
    If (fc .LT. f(vg)) Then
      DO j=0,n-1
        v(vg,j) = vc(j)
      END DO
    f(vg) = fc

! at this point the contraction is not successful,
! we must halve the distance from vs to all the
! vertices of the simplex and then continue.
! 10/31/97 - modified C program to account for 
! all vertices.

  Else
    DO row=0,n
      If (row .NE. vs) Then
        DO j=0,n-1
          v(row,j) = v(vs,j)+(v(row,j)-v(vs,j))/2.0
        END DO
      END IF
    END DO
    DO m=0,n-1
      vtmp(m) = v(vg,m)
    END DO
    f(vg) = FUNC(vtmp)
    k = k+1

    DO m=0,n-1
      vtmp(m) = v(vh,m)
    END DO
    f(vh) = FUNC(vtmp)
    k = k+1
    END IF
  END IF

! print out the value at each iteration 
  IF (iprint == 0) THEN
    Write(*,*) "Iteration ",itr
    Write(*,300) ((v(i,j),j=0,n-1),f(i),i=0,n)
  END IF

! test for convergence
  fsum = 0.0
  DO j=0,n
    fsum = fsum + f(j)
  END DO
  favg = fsum/(n+1.)
  s = 0.0
  DO j=0,n
    s = s + ((f(j)-favg)**2.)/n
  END DO
  s = sqrt(s)
  If (s .LT. EPSILON) Then
    EXIT ! Nelder Mead has converged - exit main loop
  END IF
END DO
! end main loop of the minimization
! find the index of the smallest value

  vs = 0
  DO j=0,n
    If (f(j) .LT. f(vs)) Then
      vs = j
    END IF
  END DO

!  print out the minimum

  DO m=0,n-1
    vtmp(m) = v(vs,m)
  END DO

  min = FUNC(vtmp)
  k = k+1
  write(*,*)'The minimum was found at ',v(vs,0),v(vs,1)
  write(*,250)'The value at the minimum is ',min
  write(*,*)'The number of function evaluations was',k
  write(*,*)'The number of iterations was',itr
250  FORMAT(A29,F7.4)
300  FORMAT(F11.6,F11.6,F11.6)

  return
  end subroutine simplex
! ======================================================================
end module mod_simplex
