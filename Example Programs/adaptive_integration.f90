!-----------------------------------------------------------------------------!
! MAE 5093 DOCUMENTATION | Engineering Numerical Analysis
!-----------------------------------------------------------------------------!
! >>> Adaptive Numerical Integration 
!     Simpson base
!-----------------------------------------------------------------------------!
! References: 
! * Fundamentals of Engineering Numerical Analysis by P. Moin (2012) -- Ch.1
! * Numerical Recipes: The Art of Scientific Computing, 2nd Edition (1992) 
!-----------------------------------------------------------------------------!
! Written by Omer San
!            CFDLab, Oklahoma State University, cfdlab.osu@gmail.com
!            www.cfdlab.org
! 
! Last updated: Sep. 8, 2015
!-----------------------------------------------------------------------------!


program adaptive_integration
implicit none
real*8 ::S,eI,a,b,geps

!exact solution for the given function:
eI = -0.56681975015d0

geps = 1.0d-5  !error criterion

a = -1.0d0 !lower bound
b =  1.0d0 !upper bound

!Adaptive Simpson Integration
call adaptsimpson(a,b,geps,S)

write(*,19)"exact: ", eI, 100.0d0*dabs((eI-eI)/eI)
write(*,19)"numerical: ", S, 100.0d0*dabs((S-eI)/eI)

19 format(A20,2F20.10)



end


!-----------------------------------------------------------------------------!
!Given function to integrate
!-----------------------------------------------------------------------------!
real*8 function f(x)
implicit none
real*8 :: x
f = 10.0d0*dexp(-50.0d0*dabs(x)) &
  - 0.01d0/((x-0.5d0)*(x-0.5d0) + 0.001d0) &
  + 5.0d0*dsin(5.0d0*x)
end function f

!-----------------------------------------------------------------------------!
!Adaptive integration with Simpson base
!Compute the integral of f(x) within the domain [a,b]
!uses external function f(x)
!a: lower bound of the function
!b: upper bound of the function
!geps: gloabal error criterion (user defined)
!S: integral of f(x) within the domain
!-----------------------------------------------------------------------------!
subroutine adaptsimpson(a,b,geps,S)
implicit none
real*8 ::a,b,geps,S
real*8 ::f
real*8 ::S1,S2,h,e1,e2,x,eps

!Adaptive integration:
S = 0.0d0
x = a
eps= geps/(b-a)

!span along x direction from left to right
200 continue
h  = (b-x)

if (x.ge.b) goto 300

!trial solution 
100 continue
S1 = h/6.0d0*(f(x)+4.0d0*f(x+0.5d0*h)+f(x+h))
S2 = h/12.0d0*(f(x)+4.0d0*f(x+0.25d0*h)+2.0d0*f(x+0.5d0*h)+4.0d0*f(x+0.75d0*h)+f(x+h))

e1 = 1.0d0/15.0d0*dabs(S2-S1)
e2 = h*eps

!error check
if (e1.le.e2) then !accept step
	S = S + (16.0d0*S2-S1)/15.0d0
	x = x + h
	goto 200    
else !reduce step size and perform trial solution
	h = 0.5d0*h
	goto 100
end if
300 continue !done

return
end 








