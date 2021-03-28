c
c	program DMC_Basic.f
c
c	A teaching tool of the class "Spectroscopy and Computer Modeling
c	of Molecular Systems", M.S. Physics, Bari University, Italy.
c
c	Details of the algorithm are found in the reference paper Longo, 
c	GM, et al. "The unbiased Diffusion Monte Carlo: a versatile tool
c	for two-electron systems confined in different geometries.", EPJ
c	D 2021. Please make reference to the paper if this code is used
c	as a research tool.
c
c	It allows to perform many experiments, provided the two-electron
c	system is in its singlet ground state. A list of quantum systems
c	to study, with the related coulomb energy functions, is provided
c	inside the code. Examples are H2+, HeH+, H2. A simple empty-core
c	pseudo-potential is proposed to make experiments with e.g. LiH+,
c	LiH, Li2+.
c
c	For the one-electron case, excited states are readily selected 
c	using symmetry-based nodal planes. For details: Longo GM et al,
c	Plasma Sources Science and Technology 24.6 (2015): 065019.

c
c	This program uses atomic units
c
c	This version by Savino Longo, Gaia Micca Longo, 
c			    University of Bari, 2021
c

c	ndim is 3 for 1-e systems
c	ndim is 6 for 2-e systems
c	note "activate the following 3 lines", later in the code
c

	parameter (n=100000, nist=100, ndim = 3)            

c
c
c	x(i,j) j-component of the 3N-D position of the i-th walker


        dimension x(n,ndim), jfl(n)      ! position vector and flags
	integer isto(nist)               ! histogram
	integer sto(nist)
	rmax=10
	dr = rmax/float(nist)
	xmin=-1.5; xmax=1.5              ! range istogramma


c	with these parameters we have a reasonably good performance on
c	most systems:
c 	eau=-.5				 ! initial energy estimate
c	tau=.01                  	 ! imaginary time step
c	nstep = 50000			 ! no. of time steps for each PES point
c	np0=1000                         ! initial number of walkers

c	H2

 	eau=-.5				 ! initial energy estimate
	tau=.01                  	 ! imaginary time step
	nstep = 100000			 ! no. of time steps for each PES point
	np0=1000                         ! initial number of walkers

	np=np0

	iplot=0
        do k=1,np 
c
c	"3" here becomes "6" for 2-electrons calculations
c
	 do j=1,ndim                      	 ! initialization
          x(k,j)=xmin+(xmax-xmin)*rand()
	 end do
         jfl(k)=1
	end do

	do m=1,nist                      ! histogram initialization
	 isto(m)=0
	 sto(m)=0
	end do


c	print *," "

c
c	cycle of the internuclear distances
c
c	The aim of this cycle is to provide a scan of values of the parameter
c	d, the internuclear distance. In this way a Potential Energy Surface
c	is automatically produced. In case d is not used do not erase this 
c	line but replace with "do iter = 1,1"
c

c	for H2
c
c	do iter = 0,15
c	d = 0.75 + float(iter)/4
		

	do iter = 0,20
	 d = 1 + float(iter)/4.

	 emedia=0
	 pmedio=0
	 istart=000

c
c	accelerated diffusion propagator (see the reference paper)
c
	 amp=2*sqrt(3.)*sqrt(tau)


	 do i=1,nstep		! (imaginary) time cycle

	  npn=np             ! new number of walkers (to be updated later)

          do k=1,np		! walker cycle

c
c	accelerated diffusion propagator
c	a uniform distribution is used instead of gaussian
c	see the reference paper for explanation of how it works
c
	   do j=1,ndim
            x(k,j)=x(k,j)+amp*rand()-amp/2 ! random displacement
	   end do                      
                 
c   	Green function

c
c	local energy U(x1,y1, ... )
c	
c	in this calculations the vector variables are stored in shorter
c	symbols a,b,c... for a more readable expression of the energy
c
c	a = x1, b = y1, c=z1, u=x2, v=y2, c=z2
c

	a=x(k,1)
	b=x(k,2)
	c=x(k,3)

c	activate the following line for 2-electrons calculations

c	u=x(k,4);v=x(k,5);w=x(k,6)

c	Coulomb potential energy:

c ------------------------------ begin of the "example database" -------------
c
c	activate an example hamiltonian to use it.
c	"g" is the expression of the coulomb potential g(a,b,...,d)
c

c
c	H2+ here, d is the internuclear distance. Nuclei along the z-axis
c
	g=-1/sqrt(a**2+b**2+(c-d/2)**2)-1/sqrt(a**2+b**2+(c+d/2)**2)+1/d

c	nodal plane for excited sigma_u (see explanation above)

c	if(c<0)g=1e6

c	if(a<0)g=1e6

c
c	Single H atom

c	r = sqrt(a**2+b**2+c**2)
c	g=-1/r

c	for several excited states of 1-e systems the nodal surface can be
c	found by using symmetry considerations, look for explanations: 
c	Micca Longo, G., S. Longo, and D. Giordano. "Spherically confined 
c	H2+:^{2}\Sigma_{g}^{+} and^{2}\Sigma_{u}^{+} states." Physica 
c	Scripta 90.2 (2015): 025403.
c


c	use of symmetry for 2p

c	if(a+b+c<0)g=1e6  ! a=0 --> x=0 is a node

c	if(r.lt.d)g=1e6

c
c	H2 (2 electrons)
c
c	g=-1/sqrt(a**2+b**2+(c-d/2)**2)-1/sqrt(a**2+b**2+(c+d/2)**2)
c     ^-1/sqrt(u**2+v**2+(w-d/2)**2)-1/sqrt(u**2+v**2+(w+d/2)**2)	
c     ^+1/sqrt((a-u)**2+(b-v)**2+(c-w)**2)
c     ^+1/d


c
c	He (2 electrons)
c
c	g=-2/sqrt(a**2+b**2+c**2)
c     ^ -2/sqrt(u**2+v**2+w**2)	
c     ^ +1/sqrt((a-u)**2+(b-v)**2+(c-w)**2)

c
c	Example of pseudo-potential: Li atom
c
c	Ashcroft empty-core PP with U(r) =0 for r<Rc, U=-Z/r for r> Rc.
c	Rc = 1.26 a.u. in Bastea, Marina, and Sorin Bastea. "Electrical 
c	conductivity of lithium at mega-bar pressures." Physical Review B 
c	65.19 (2002): 193104.
c

c
c	Single Li atom
c

c	r = sqrt(a**2+b**2+c**2)
c	g = 0		! zero for r < Rc
c	if(r.gt.1.75) g=-1/r     ! Li+    Z = 3+ -2e (core electrons, 1s2)

c	if(c.lt.0)g=1e6 		! to select the 2pz orbital


c
c	Application to Li_2 + molecule
c

c	g = 1/d
c	rc = 1.26
c	r1 = sqrt(a**2+b**2+(c-d/2)**2)
c	if(r1.gt.rc) g=g-1/r1
c	r2 = sqrt(a**2+b**2+(c+d/2)**2)
c	if(r2.gt.rc) g=g-1/r2

c ------------------------------ end of the "example database" ---------------

c ---------------avoid changing anything below this point: core algorithm ----

          cc=rand()  
	  pw=exp(-tau*(g-eau)) 
	  if (cc>pw) jfl(k)=0           !walker survives with probability px
	  if (cc<(pw-1)) then  		!create new walker with probability px-1
	   npn=npn+1
c
c	"3" here becomes "6" for 2-electrons calculations
c
	     do j=1,ndim
	      x(npn,j)=x(k,j)         	!download the old vector to the new one 
	     end do

c	   x(npn,1)=x(k,1); x(npn,2)=x(k,2); x(npn,3)=x(k,3); 
	   jfl(npn)=1
	  end if
	 end do
	 np=npn                     	!new walker number
	 iconta=0
	 do k=1,np                      !compact the list
	  if (jfl(k)==1) then
           iconta=iconta+1              !count active walkers
c
c	"3" here becomes "6" for 2-electrons calculations
c
	     do j=1,ndim
	     x(iconta,j)=x(k,j)         !download the old vector to the new one
	     end do
	   jfl(iconta)=1                !activate new walkers
          end if                     
	 end do
	 np=iconta
	 npold=np
	 np=iconta
	 if(i.gt.1)then
  	  eau=eau+tau*alog(float(np0)/float(np)); !energy self-regulation
	  if(i.gt.1000)eau=eau+.001*alog(pmedio/float(np)); 
	 end if
	 em=float(i-istart)
	 emedia=(em*emedia+eau)/(em+1)           !time average
	 pmedio=(em*pmedio+float(np))/(em+1)     
	 iplot=iplot+1
	 if(iplot.ge.nstep)then
	  iplot=0
	  print *, d, emedia
	 end if

	 if(i.gt.istart+1000)then
	  err=abs((float(np)-float(np0))/float(np0))
	  if(err.lt..01)eau=emedia
	 end if	

c ------------------------------ avoid changing anything above this point ----


c
c	histogram production and print (if required)
c	preliminar version: walkers vs. r
c	dr = r step
c

	 jsto=jsto+1
	 if(i.gt.10000.and.jsto.ge.10)then
	  jsto=0
	  do k=1,np
 	   ii = (sqrt(x(k,1)**2+x(k,2)**2+x(k,3)**2))/dr
	   if(ii>0.and.ii<nist)	   isto(ii)=isto(ii)+1
	  end do
	 end if
	
	 if(i.eq.nstep)then		!print at the end
c	 print *,d,"histogram"
	  do j=1,nist
c	   print *,j*dr,isto(j)
	  end do
	 end if

 
	 end do		!end time cycle

	end do		!end parameter cycle

	end	
	
