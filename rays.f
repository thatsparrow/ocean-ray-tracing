c ============================================================
c  RAYS.F  –  Ocean wave ray-tracing program
c  author: jann benjamin (la jolla, 2016)
c  last edited: 4/1/2026
c
c  What this program does:
c    Traces the paths of ocean wave "rays" as they travel
c    from deep water toward shore, bending as the water
c    gets shallower (refraction).
c
c  Fortran basics for non-Fortran readers:
c    - Lines starting with "c" are comments (ignored by compiler)
c    - "d0" suffix on numbers means double-precision float, e.g. 1.0d0
c    - "**" means "raise to the power of", e.g. x**2 = x squared
c    - Arrays are declared with "dimension" and indexed from 1
c    - "common" blocks are how Fortran shares variables between
c      functions without passing them as arguments
c    - "implicit double precision(a-h,o-z)" means any variable
c      whose name starts with a-h or o-z is automatically treated
c      as a double-precision float (variables starting i-n default
c      to integers, which is a Fortran convention)
c ============================================================

      program rays

      implicit double precision(a-h,o-z)

c     -- Declare kx0,ky0 etc. as double precision explicitly.
c        (BUG 1 FIX: was "double prevision", a typo)
      double precision kx0,ky0,kxray(5000),kyray(5000),ktotal

c     -- 1-D arrays to hold the (x,y) position along each ray
      dimension xray(5000),yray(5000)

c     -- 2-D array to hold depth values for the contour plot
      dimension hplot(100,200)

c     -- Shared constants for the bathymetry functions h, hx, hy.
c        "common" lets these variables be seen by the subroutines
c        at the bottom of the file without passing them explicitly.
      common/depth/alpha,eps,beta

c     -- Shared origin for the graphics coordinate system
      common/xor/xor,yor

c     -- Define pi as a double-precision constant
      data pi/3.1415926535d0/

c     -- Open a PostScript file for graphical output
      open(8,file="rays.ps")
      call opengr(8)

c ---- Physical and numerical setup -------------------------

c     Gravitational acceleration (m/s^2)
      grav=9.8d0

c     x-coordinate where each ray starts (1000 m offshore).
c     x<0 is the ocean; x=0 is the shoreline.
      x0=-1000.d0

c     Wave period (seconds) and angular frequency omega = 2*pi/T
      period=12.d0
      omega=2.d0*pi/period

c     Deep-water wavenumber from the dispersion relation
c     for deep water:  omega^2 = g * k
      ktotal=omega**2/grav
      write(6,'("wavelength=",f6.1)')2.d0*pi/ktotal

c     Initial wavenumber components: waves arrive at 45 degrees
c     to the shore, so kx = ky = ktotal / sqrt(2)
      kx0=ktotal/dsqrt(2.d0)
      ky0=-ktotal/dsqrt(2.d0)

c     Length of the coastal domain in the y-direction (m)
      dist=2000.d0

c     Bathymetry parameters (see function h below):
c       alpha  - controls how steeply the seafloor slopes up
c       beta   - alongshore wavenumber of the depth variations
c       eps    - amplitude of the alongshore depth variation
      alpha=1.d0/20.d0
      beta=2.d0*pi/1000.d0
      eps=.3d0

c ---- Draw a contour plot of the bathymetry ----------------

      xor=0.d0
      yor=0.d0
      width=.8d0
      height=1.6d0

c     Fill hplot(i,j) with depth values on a 100x200 grid
      do j=1,200
        y=dble(j-1)/dble(199)*dist
        do i=1,100
          x=dble(i-100)/dble(99)*.5d0*dist
          hplot(i,j)=h(x,y)
        enddo
      enddo
      call plot2d(20,hplot,100,200,8,width,height,cint)

c ---- Ray-tracing setup ------------------------------------

c     Number of rays to trace (evenly spaced along the coast)
      nrays=20

c     Maximum group velocity (shallow-water limit: cg = sqrt(g*h))
c     used to set a stable time step
      cgmax=dsqrt(grav*h(x0,0.d0))

c     Time step in seconds
      dt=10.d0/cgmax

c ---- Main loop: trace each ray ----------------------------

      do 2000 j=1,nrays

c       Starting position: x=x0 (offshore), y spread along coast
        xray(1)=x0
        yray(1)=dble(j-1)/dble(nrays)*dist

c       Starting wavenumber (same for all rays: same wave direction)
        kxray(1)=kx0
        kyray(1)=ky0

c       Follow this ray for up to 1000 time steps
        do 1000 i=1,1000

c         Stop if the ray has left the domain
          if(xray(i).lt.x0)go to 1001
          if(xray(i).ge.0.d0)go to 1001
          if(yray(i).lt.0.d0)go to 1001
          if(yray(i).gt.dist)go to 1001

c         Cache depth and its gradients at the current position
          htem=h(xray(i),yray(i))
          hxtem=hx(xray(i),yray(i))
          hytem=hy(xray(i),yray(i))

c         Total wavenumber magnitude at current position
          ktotal=dsqrt(kxray(i)**2+kyray(i)**2)

c         Non-dimensional depth argument for tanh dispersion relation
          arg=ktotal*htem

c         Angular frequency from the full linear dispersion relation:
c           omega = sqrt(g * k * tanh(k*h))
c         (BUG 2 FIX: removed extra closing parenthesis after dtanh(arg))
          omega=dsqrt(grav)*ktotal*dtanh(arg)

c         Group velocity cg = d(omega)/dk  (speed at which wave energy travels).
c         (BUG 3 FIX: was ".50" (single precision); changed to ".5d0")
          groupspeed=.5d0/omega*grav*
     *      (dtanh(arg)+arg/dcosh(arg)**2)

c         Partial derivative of omega with respect to depth h.
c         Needed to compute how the wavenumber changes as depth changes.
          dwdh=.5d0/omega*grav*
     *      ktotal**2/dcosh(arg)**2

c         Ray equations (Hamilton's equations for wave rays):
c           dx/dt  = +d(omega)/d(kx)  =  cg * kx/|k|
c           dy/dt  = +d(omega)/d(ky)  =  cg * ky/|k|
c           dkx/dt = -d(omega)/dx     = -(d(omega)/dh) * dh/dx
c           dky/dt = -d(omega)/dy     = -(d(omega)/dh) * dh/dy
c         (BUG 4 FIX: these four RHS terms were never computed)
          dxdt  = groupspeed*kxray(i)/ktotal
          dydt  = groupspeed*kyray(i)/ktotal
          dkxdt = -dwdh*hxtem
          dkydt = -dwdh*hytem

c         Time integration using leapfrog (2nd-order accurate):
c           first step must be a simple forward Euler step because
c           there is no "previous" point yet.
          if(i.eq.1)then

c           Forward Euler step for the first point only
            xray(2) =xray(1) +dt*dxdt
            yray(2) =yray(1) +dt*dydt
            kxray(2)=kxray(1)+dt*dkxdt
            kyray(2)=kyray(1)+dt*dkydt

          else

c           Leapfrog: new value = value two steps ago + 2*dt*(current rate)
c           (BUG 5 FIX: was dkdxdt/dkdydt, which don't exist; correct names
c            are dkxdt/dkydt as assigned above)
            xray(i+1) =xray(i-1) +2.d0*dt*dxdt
            yray(i+1) =yray(i-1) +2.d0*dt*dydt
            kxray(i+1)=kxray(i-1)+2.d0*dt*dkxdt
            kyray(i+1)=kyray(i-1)+2.d0*dt*dkydt

          endif

c         Convert ray positions to graphics page coordinates and draw segment
          x1=xor+(xray(i)  -x0)/dabs(x0)*width
          y1=yor+(yray(i)      /dist)     *height
          x2=xor+(xray(i+1)-x0)/dabs(x0)*width
          y2=yor+(yray(i+1)    /dist)     *height
          call line(8,x1,y1,x2,y2,4)

1000    continue
1001    continue

2000  continue

      call closegr(8)
      stop
      end

c ============================================================
c  BATHYMETRY FUNCTIONS
c
c  The seafloor depth is modelled as a planar slope (deepening
c  offshore in x) with a gentle sinusoidal ripple along the
c  coast (in y):
c
c    h(x,y) = 1 - alpha * x * (1 + eps * cos(beta * y))
c
c  Because x <= 0 in the ocean, alpha*x <= 0, so h >= 1 > 0.
c
c  Note: all three functions share alpha, eps, beta via the
c  common block /depth/.  The common block separator is a
c  comma, not a slash.
c  (BUG 8 FIX, applies to all three functions below:
c   was "common/depth/alpha/eps/beta" which accidentally placed
c   eps and beta into separate unnamed common blocks instead of
c   /depth/, so they were never shared with the main program)
c ============================================================

      function h(x,y)
      implicit double precision(a-h,o-z)
      common/depth/alpha,eps,beta
      h=1.d0-alpha*x*(1.d0+eps*dcos(beta*y))
      return
      end

c ------------------------------------------------------------
c  HX: partial derivative of h with respect to x
c    d/dx [ 1 - alpha*x*(1 + eps*cos(beta*y)) ]
c        = -alpha * (1 + eps*cos(beta*y))
c
c  (BUG 6 FIX: the original function body was a copy-paste of h(x,y)
c   and stored the result in "h" instead of "hx", so it returned
c   the depth rather than its x-derivative)
c ------------------------------------------------------------

      function hx(x,y)
      implicit double precision(a-h,o-z)
      common/depth/alpha,eps,beta
      hx=-alpha*(1.d0+eps*dcos(beta*y))
      return
      end

c ------------------------------------------------------------
c  HY: partial derivative of h with respect to y
c    d/dy [ 1 - alpha*x*(1 + eps*cos(beta*y)) ]
c        = alpha * eps * beta * x * sin(beta*y)
c
c  (BUG 7 FIX: was "alpha**eps" meaning alpha raised to the power
c   eps, which is dimensionally wrong and numerically very different;
c   the correct operation is multiplication: alpha*eps)
c ------------------------------------------------------------

      function hy(x,y)
      implicit double precision(a-h,o-z)
      common/depth/alpha,eps,beta
      hy=alpha*eps*beta*x*dsin(beta*y)
      return
      end
