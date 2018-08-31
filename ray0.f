    program rays  

    implicit double precision(a-h,o-z)
    double prevision kx0,ky0,kxray(5000),kyray(5000),ktotal
    dimension xray(5000),yray(5000)
    dimension hplot(100,200)
    common/depth/alpha,eps,beta
    common/xor/xor,yor
    data pi/3.1415926535d0/
    open(8,file="rays.ps")
    call opengr(8)

c gravity constant:

    grav=9.8d0

c x-coordinate of the beginning point of each ray (1000 meters offshore):

    x0=-1000.d0

c deep water values:

    period=12.d0
    omega=2.d0*pi/period
    ktotal=omega**2/grav
    write(6,'("wavelength=",f6.1)')2.d0*pi/ktotal
    kx0=ktotal/dsqrt(2.d0)
    ky0=-ktotal/dsqrt(2.d0)

c length of the domain along the coast:

    dist=2000.d0

c parameters that define the bathymetry:

    alpha=1.d0/20.d0
    beta=2.d0*pi/1000.d0
    eps=.3d0

c draws a contour plot of the bathymetry
    xor=0.d0
    yor=0.d0
    width=.8d0
    height=1.6d0
    do j=1,200
    y=dble(j-1)/dble(199)*dist
    do i=1,100
    x=dble(i-100)/dble(99)*.5d0*dist
    hplot(i,j)=h(x,y)
    enddo
    enddo
    call plot2d(20,hplot,100,200,8,width,height,cint)

c number of rays to be followed:

    nrays=20

c maximum group velocity:

    cgmax=dsqrt(grav*h(x0,0.d0))

c time step:

    dt=10.d0/cgmax

c loop through the rays:

    do 2000 j=1,nrays

c for each ray set the starting values of location and wavenumber 

    xray(1)=x0
    yray(1)=dble(j-1)/dble(nrays)*dist
    kxray(1)=kx0
    kyray(1)=ky0

c follow each ray toward shore:

    do 1000 i=1,1000

c check to see if the ray is still within bounds:

    if(xray(i).lt.x0)go to 1001
    if(xray(i).ge.0.d0)go to 1001
    if(yray(i).lt.0.d0)go to 1001
    if(yray(i).gt.dist)go to 1001

c temporarily store the depth and its derivatives:

    htem=h(xray(i),yray(i))
    hxtem=hx(xray(i),yray(i))
    hytem=hy(xray(i),yray(i))

c compute the total wavenumber, etc:

    ktotal=dsqrt(kxray(i)**2+kyray(i)**2)
    arg=ktotal*htem
    omega=dsqrt(grav)*ktotal*dtanh(arg))
    groupspeed=.50/omega*grav*
    *  (dtanh(arg)+arg/dcosh(arg)**2)
     dwdh=.5d0/omega*grav*
     *  ktotal**2/dcosh(arg)**2

c advance ray one time stop:

    if(i.eq.1)then

c the first step is a forward step:

    xray(2)=xray(1)+dt*dxdt
    yray(2)=yray(1)+dt*dydt
    kxray(2)=kxray(1)+dt*dkxdt
    kyray(2)=kyray(1)+dt*dkydt
    else

c subsequent steps are leapfrog steps:

    xray(i+1)=xray(i-1)+2.d0*dt*dxdt
    yray(i+1)=yray(i-1)+2.d0*dt*dydt
    kxray(i+1)=kxray(i-1)+2.d0*dt*dkdxdt
    kyray(i+1)=kyray(i-1)+2.d0*dt*dkdydt
    endif

c draw the current ray segment:

    x1=xor+(xray(i)-x0)/dabs(x0)*width
    y1=yor+(yray(i)/dist)*height
    x2=xor+(xray(i+1)-x0)/dabs(x0)*width
    y2=yor+(yray(i+1)/dist)*height
    call line(8,x1,y1,x2,y2,4)

1000 continue
1001 continue

c continue with the next ray 

2000 continue

c exit after the last ray

    call closegr(8)
    stop
    end

    function h(x,y)
    implicit double precision(a-h,o-z)
    common/depth/alpha/eps/beta
    h=1.d0-alpha*x*(1.d0+eps*dcos(beta*y))
    return
    end

    function hx(x,y)
    implicit double precision(a-h,o-z)
    common/depth/alpha/eps/beta
    h=1.d0-alpha*x*(1.d0+eps*dcos(beta*y))
    return
    end

    function hy(x,y)
    implicit double precision(a-h,o-z)
    common/depth/alpha/eps/beta
    hy=alpha**eps*beta*x*dsin(beta*y)
    return
    end