program Pendulum


real::dt,y,x,f,t0,tf,x0,y0,t,k11,k21,k31,k41,l11,l21,l31,l41,theta0,thetarad,thetarad0,omega0,v0,theta,omega,g , &
      k1,k22,k32,k42,l12,l22,l32,l42,xdot,g2,s,r,s0
      
real,parameter::pi=acos(-1.0),l=1.0,xor=0.0,yor=0.0,gc=9.81,k=5.0,m=1.0
integer::n,i




print*,"give the initial time and final time"
read*,t0,tf
print*,"give the number of interval"
read*,n
!print*,"give the length of the pendulum "
!read*,R
print*,"give the initial value of stretch x ,theta(less than 40 deg) and omega"
read*,x0,theta0,omega0



thetarad0=(pi*theta0/180)

r=(l+x0)

open(30,file="elasticpen.dat")
!open(35,file="modeuler")
dt=(tf-t0)/n
print*,"dt=",dt
omega=omega0
thetarad=thetarad0
t=t0


y0=yor-(l+x0)*cos(thetarad)

s0=xor+(l+x0)*sin(thetarad)

y=yor+y0
s=xor+s0

do i=0,n
!==================================================================================================
k11=dt*f(thetarad,x,omega,xdot,t)
l11=dt*g(thetarad,x,omega,xdot,t)
k21=dt*f(thetarad+k11*0.5,x+k11*0.5,omega+l11*0.5,xdot+l11*0.5,t+dt*0.5)
l21=dt*g(thetarad+k11*0.5,x+k11*0.5,omega+l11*0.5,xdot+l11*0.5,t+dt*0.5)
k31=dt*f(thetarad+k21*0.5,x+k21*0.5,omega+l21*0.5,xdot+l21*0.5,t+dt*0.5)
l31=dt*g(thetarad+k21*0.5,x+k21*0.5,omega+l21*0.5,xdot+l21*0.5,t+dt*0.5)
k41=dt*f(thetarad+k31*0.5,x+k31*0.5,omega+l31*0.5,xdot+l31*0.5,t+dt*0.5)
l41=dt*g(thetarad+k31*0.5,x+k31*0.5,omega+l31*0.5,xdot+l31*0.5,t+dt*0.5)
!===================================================================================================
k12=dt*f2(thetarad,x,omega,xdot,t)
l12=dt*g2(thetarad,x,omega,xdot,t)
k22=dt*f2(thetarad+k12*0.5,x+k12*0.5,omega+l12*0.5,xdot+l12*0.5,t+dt*0.5)
l22=dt*g2(thetarad+k12*0.5,x+k12*0.5,omega+l12*0.5,xdot+l12*0.5,t+dt*0.5)
k32=dt*f2(thetarad+k22*0.5,x+k22*0.5,omega+l22*0.5,xdot+l22*0.5,t+dt*0.5)
l32=dt*g2(thetarad+k22*0.5,x+k22*0.5,omega+l22*0.5,xdot+l22*0.5,t+dt*0.5)
k42=dt*f2(thetarad+k32*0.5,x+k32*0.5,omega+l32*0.5,xdot+l32*0.5,t+dt*0.5)
l42=dt*g2(thetarad+k32*0.5,x+k32*0.5,omega+l32*0.5,xdot+l32*0.5,t+dt*0.5)
!=============================================================================================

theta=180*thetarad/pi

if(mod(i,5)==0)then
write(30,*)0,0,s,y,t,(l+x),thetarad,omega
end if

!===================================================================================================
t=i*dt+t0
thetarad=thetarad+(k11+2*k21+2*k31+k41)/6.0
omega=omega+(l11+2*l11+2*l31+l41)/6.0

x=x+(k12+2*k22+2*k32+k42)/6.0
xdot=xdot+(l12+2*l12+2*l32+l42)/6.0

y=yor-(l+x)*cos(thetarad)

s=xor+(l+x)*sin(thetarad)


!print*,x,y
!write(30,*)x,y,t

end do

!close(35)
!call system ('gnuplot -p euler_plot.plt')
!call system ('gnuplot -p modeuler_plot.plt')
end program
!---------------------------------------------------------------------
real function f(thetarad,x,omega,xdot,t)


real,intent(in)::thetarad,omega,t
real,parameter::pi=acos(-1.0),l=1.0,xor=0.0,yor=0.0,gc=9.81,k=5.0,m=1.0

f=omega
end function f
!======================================================================

real function f2(thetarad,x,omega,xdot,t)


real,intent(in)::thetarad,omega,t,xdot
real,parameter::pi=acos(-1.0),l=1.0,xor=0.0,yor=0.0,gc=9.81,k=5.0,m=1.0

f=xdot
end function f2

!----------------------------------------------------------------------
real function g(thetarad,x,omega,xdot,t)


real,intent(in)::thetarad,omega,t,xdot
real,parameter::pi=acos(-1.0),l=1.0,xor=0.0,yor=0.0,gc=9.81,k=5.0,m=1.0

g=-gc*sin(thetarad)/(l+x)-2*xdot*omega/(l+x)

end function g
!========================================================================
real function g2(thetarad,x,omega,xdot,t)


real,intent(in)::thetarad,omega,t,xdot,x
real,parameter::pi=acos(-1.0),l=1.0,xor=0.0,yor=0.0,gc=9.81,k=5.0,m=1.0

g2=(l+x)*(omega)**2-(k*x/m)+gc*cos(thetarad)

end function g2
!====================================================================



