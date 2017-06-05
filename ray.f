    program rays  
c
c finds the ray paths beginning at x=x0 and various y
c in an ocean on x<0 with depth h(x,y)
c the x-derivative of h is hx(x,y)
c the y-derivative of h is hy(x,y)
c the initial number wavenumber is (kx0, ky0)
c the location along the ray path is (xray,yray)
c the wavenumber along the ray path is (kxray,kyray)
c mks units for everything
c
    implicit double precision(a-h,o-z)
    double prevision kx0,ky0,kxray(5000),kyray(5000),ktotal
    dimension xray(5000),yray(5000)
    dimension hplot(100,200)
    common/depth/alpha,eps,beta
    common/xor/xor,yor
    data pi/3.1415926535d0/
    open(8,file="rays.ps")
    call opengr(8)
c
c set the starting parameters
c gravity constant:
    grav=9.8d0
c x-coordinate of the beginning point of each ray
c (1000 meters offshore):
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
c
c draw a contour plot of the bathymetry
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
c
c number of rays to be followed:
    nrays=20
c maximum group velocity:
    cgmax=dsqrt(grav*h(x0,0.d0))
c time step:
    dt=10.d0/cgmax
c
c loop through the rays:
    do 2000 j=1,nrays
c for each ray set the starting values of location and wavenumber 
    xray(1)=x0
    yray(1)=dble(j-1)/dble(nrays)*dist
    kxray(1)=kx0
    kyray(1)=ky0
c
c follow each ray toward shore:
    do 1000 i=1,1000
c
c check to see if the ray is still within bounds:
    if(xray(i).lt.x0)go to 1001
    if(xray(i).ge.0.d0)go to 1001
    if(yray(i).lt.0.d0)go to 1001
    if(yray(i).gt.dist)go to 1001
c
c temporarily store the depth and its derivatives:
    htem=h(xray(i),yray(i))
    hxtem=hx(xray(i),yray(i))
    hytem=hy(xray(i),yray(i))
c
c compute the total wavenumber, etc:
    
