!******************************************************************************
!*                    Code generated with sympy 0.7.5-git                     *
!*                                                                            *
!*              See http://www.sympy.org/ for more information.               *
!*                                                                            *
!*                      This file is part of 'libslater'                      *
!******************************************************************************


interface
REAL*8 function cou_0_0(a, b, R)
implicit none
REAL*8, intent(in) :: a
REAL*8, intent(in) :: b
REAL*8, intent(in) :: R
end function
end interface
interface
REAL*8 function cou_taylor_0_0(a, b, R)
implicit none
REAL*8, intent(in) :: a
REAL*8, intent(in) :: b
REAL*8, intent(in) :: R
end function
end interface
interface
REAL*8 function cou_1_0(a, b, R)
implicit none
REAL*8, intent(in) :: a
REAL*8, intent(in) :: b
REAL*8, intent(in) :: R
end function
end interface
interface
REAL*8 function cou_taylor_1_0(a, b, R)
implicit none
REAL*8, intent(in) :: a
REAL*8, intent(in) :: b
REAL*8, intent(in) :: R
end function
end interface
interface
REAL*8 function cou_1_1(a, b, R)
implicit none
REAL*8, intent(in) :: a
REAL*8, intent(in) :: b
REAL*8, intent(in) :: R
end function
end interface
interface
REAL*8 function cou_taylor_1_1(a, b, R)
implicit none
REAL*8, intent(in) :: a
REAL*8, intent(in) :: b
REAL*8, intent(in) :: R
end function
end interface
interface
REAL*8 function cou_2_0(a, b, R)
implicit none
REAL*8, intent(in) :: a
REAL*8, intent(in) :: b
REAL*8, intent(in) :: R
end function
end interface
interface
REAL*8 function cou_taylor_2_0(a, b, R)
implicit none
REAL*8, intent(in) :: a
REAL*8, intent(in) :: b
REAL*8, intent(in) :: R
end function
end interface
interface
REAL*8 function cou_2_1(a, b, R)
implicit none
REAL*8, intent(in) :: a
REAL*8, intent(in) :: b
REAL*8, intent(in) :: R
end function
end interface
interface
REAL*8 function cou_taylor_2_1(a, b, R)
implicit none
REAL*8, intent(in) :: a
REAL*8, intent(in) :: b
REAL*8, intent(in) :: R
end function
end interface
interface
REAL*8 function cou_2_2(a, b, R)
implicit none
REAL*8, intent(in) :: a
REAL*8, intent(in) :: b
REAL*8, intent(in) :: R
end function
end interface
interface
REAL*8 function cou_taylor_2_2(a, b, R)
implicit none
REAL*8, intent(in) :: a
REAL*8, intent(in) :: b
REAL*8, intent(in) :: R
end function
end interface
interface
REAL*8 function cou_3_0(a, b, R)
implicit none
REAL*8, intent(in) :: a
REAL*8, intent(in) :: b
REAL*8, intent(in) :: R
end function
end interface
interface
REAL*8 function cou_taylor_3_0(a, b, R)
implicit none
REAL*8, intent(in) :: a
REAL*8, intent(in) :: b
REAL*8, intent(in) :: R
end function
end interface
interface
REAL*8 function cou_3_1(a, b, R)
implicit none
REAL*8, intent(in) :: a
REAL*8, intent(in) :: b
REAL*8, intent(in) :: R
end function
end interface
interface
REAL*8 function cou_taylor_3_1(a, b, R)
implicit none
REAL*8, intent(in) :: a
REAL*8, intent(in) :: b
REAL*8, intent(in) :: R
end function
end interface
interface
REAL*8 function cou_3_2(a, b, R)
implicit none
REAL*8, intent(in) :: a
REAL*8, intent(in) :: b
REAL*8, intent(in) :: R
end function
end interface
interface
REAL*8 function cou_taylor_3_2(a, b, R)
implicit none
REAL*8, intent(in) :: a
REAL*8, intent(in) :: b
REAL*8, intent(in) :: R
end function
end interface
interface
REAL*8 function cou_3_3(a, b, R)
implicit none
REAL*8, intent(in) :: a
REAL*8, intent(in) :: b
REAL*8, intent(in) :: R
end function
end interface
interface
REAL*8 function cou_taylor_3_3(a, b, R)
implicit none
REAL*8, intent(in) :: a
REAL*8, intent(in) :: b
REAL*8, intent(in) :: R
end function
end interface
interface
REAL*8 function olp_0_0(a, b, R)
implicit none
REAL*8, intent(in) :: a
REAL*8, intent(in) :: b
REAL*8, intent(in) :: R
end function
end interface
interface
REAL*8 function olp_taylor_0_0(a, b, R)
implicit none
REAL*8, intent(in) :: a
REAL*8, intent(in) :: b
REAL*8, intent(in) :: R
end function
end interface
interface
REAL*8 function olp_1_0(a, b, R)
implicit none
REAL*8, intent(in) :: a
REAL*8, intent(in) :: b
REAL*8, intent(in) :: R
end function
end interface
interface
REAL*8 function olp_taylor_1_0(a, b, R)
implicit none
REAL*8, intent(in) :: a
REAL*8, intent(in) :: b
REAL*8, intent(in) :: R
end function
end interface
interface
REAL*8 function olp_1_1(a, b, R)
implicit none
REAL*8, intent(in) :: a
REAL*8, intent(in) :: b
REAL*8, intent(in) :: R
end function
end interface
interface
REAL*8 function olp_taylor_1_1(a, b, R)
implicit none
REAL*8, intent(in) :: a
REAL*8, intent(in) :: b
REAL*8, intent(in) :: R
end function
end interface
interface
REAL*8 function olp_2_0(a, b, R)
implicit none
REAL*8, intent(in) :: a
REAL*8, intent(in) :: b
REAL*8, intent(in) :: R
end function
end interface
interface
REAL*8 function olp_taylor_2_0(a, b, R)
implicit none
REAL*8, intent(in) :: a
REAL*8, intent(in) :: b
REAL*8, intent(in) :: R
end function
end interface
interface
REAL*8 function olp_2_1(a, b, R)
implicit none
REAL*8, intent(in) :: a
REAL*8, intent(in) :: b
REAL*8, intent(in) :: R
end function
end interface
interface
REAL*8 function olp_taylor_2_1(a, b, R)
implicit none
REAL*8, intent(in) :: a
REAL*8, intent(in) :: b
REAL*8, intent(in) :: R
end function
end interface
interface
REAL*8 function olp_2_2(a, b, R)
implicit none
REAL*8, intent(in) :: a
REAL*8, intent(in) :: b
REAL*8, intent(in) :: R
end function
end interface
interface
REAL*8 function olp_taylor_2_2(a, b, R)
implicit none
REAL*8, intent(in) :: a
REAL*8, intent(in) :: b
REAL*8, intent(in) :: R
end function
end interface
interface
REAL*8 function olp_3_0(a, b, R)
implicit none
REAL*8, intent(in) :: a
REAL*8, intent(in) :: b
REAL*8, intent(in) :: R
end function
end interface
interface
REAL*8 function olp_taylor_3_0(a, b, R)
implicit none
REAL*8, intent(in) :: a
REAL*8, intent(in) :: b
REAL*8, intent(in) :: R
end function
end interface
interface
REAL*8 function olp_3_1(a, b, R)
implicit none
REAL*8, intent(in) :: a
REAL*8, intent(in) :: b
REAL*8, intent(in) :: R
end function
end interface
interface
REAL*8 function olp_taylor_3_1(a, b, R)
implicit none
REAL*8, intent(in) :: a
REAL*8, intent(in) :: b
REAL*8, intent(in) :: R
end function
end interface
interface
REAL*8 function olp_3_2(a, b, R)
implicit none
REAL*8, intent(in) :: a
REAL*8, intent(in) :: b
REAL*8, intent(in) :: R
end function
end interface
interface
REAL*8 function olp_taylor_3_2(a, b, R)
implicit none
REAL*8, intent(in) :: a
REAL*8, intent(in) :: b
REAL*8, intent(in) :: R
end function
end interface
interface
REAL*8 function olp_3_3(a, b, R)
implicit none
REAL*8, intent(in) :: a
REAL*8, intent(in) :: b
REAL*8, intent(in) :: R
end function
end interface
interface
REAL*8 function olp_taylor_3_3(a, b, R)
implicit none
REAL*8, intent(in) :: a
REAL*8, intent(in) :: b
REAL*8, intent(in) :: R
end function
end interface

