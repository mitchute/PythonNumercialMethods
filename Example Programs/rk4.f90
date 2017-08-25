!-----------------------------------------------------------------------------!
! MAE 5093 DOCUMENTATION | Engineering Numerical Analysis
!-----------------------------------------------------------------------------!
! >>> Basic ODE solvers
!     Euler forward, Euler backward, Trapezoidal (Crank-Nicolson)
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


program rk4
implicit none
real*8 ::h,tmax,t,w
integer::j,k,n,np,ne
real*8,allocatable ::y(:),k1(:),k2(:),k3(:),k4(:),ye(:)

!Solve y'' + w*w*y = 0

w = 4.0d0
h = 0.01d0
tmax = 6.0d0

n = nint(tmax/h)

ne = 2

allocate(y(ne))

allocate(ye(ne))

allocate(k1(ne))
allocate(k2(ne))
allocate(k3(ne))
allocate(k4(ne))

!Initial condition
y(1) = 1.0d0
y(2) = 0.0d0

ye(1) = 1.0d0
ye(2) = 0.0d0

open(12, file="numerical.plt")
write(12,*)'variables ="t","y","dy"'
t = 0.0d0
write(12,*)t,y(1),y(2)

open(15, file="euler.plt")
write(15,*)'variables ="t","y","dy"'
write(15,*)t,ye(1),ye(2)



!RK4
do j=1,n

   t = dfloat(j)*h

   call RHS(ne,y,t,w,h,k1)
   
   call RHS(ne,y+k1/2.0d0,t+h/2.0d0,w,h,k2)

   call RHS(ne,y+k2/2.0d0,t+h/2.0d0,w,h,k3)

   call RHS(ne,y+k3,t+h,w,h,k4)

   do k=1,ne
   y(k) = y(k) + (k1(k)+2.0d0*(k2(k)+k3(k))+k4(k))/6.0d0
   end do

   do k=1,ne
   ye(k) = ye(k) + k1(k)
   end do
   
write(12,*)t,y(1),y(2)
write(15,*)t,ye(1),ye(2)

end do

!Plot
close(12)
close(15)
   
! Writing exact solution using np points
np = 2000 !use 2000 points to plot curve
open(12, file="exact_curve.plt")
write(12,*)'variables ="t","y","dy"'
	do j=0,np
		t = dfloat(j)*(tmax)/(dfloat(np))
		write(12,*) t,dcos(w*t),-w*dsin(w*t)
	end do
close(12)



end


subroutine RHS(ne,y,t,w,h,f)
implicit none
integer::ne,k
real*8 ::t,w,h,y(ne),f(ne)

f(1) = y(2)
f(2) =-w*w*y(1)


do k=1,ne
f(k) = h*f(k) 
end do
 
return
end




