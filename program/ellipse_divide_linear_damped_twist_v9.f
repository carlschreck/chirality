      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!
      !!  Completely overdamped dynamcics with horizontal PBC's
      !!
      !!  Cells modeled as ellipsoidal obejects
      !!
      !!  Layer of cells near boundary is propagated, all other
      !!  particles are removed each growth step
      !!
      !!  Constant torque is applied to all growing cells
      !!
      !!  F = b*m*dr/dt (m=1 implicit)   
      !!  T = b*I*dth/dt (I=inertia, th=orientation angle)   
      !!
      !!  Carl Schreck
      !!  8/15/2019
      !!  Wolfram Moebius & Anton Souslov
      !!  from 11/2020
      !!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      program ellipse

      implicit none
      integer Ntot
      parameter(Ntot=2**18)
      double precision pi
      parameter(pi=3.1415926535897932d0)
      double precision x(Ntot),y(Ntot),th(Ntot),D(Ntot),D1,exp,ran2
      double precision ftol1,fret,alpha0,width,rate0(Ntot),alpha(Ntot)
      double precision dt,rate(Ntot),scale(Ntot),layerwidth,layerdepth
      double precision b,s,phi,tdiv,Lx,desync,propdepth,bounddepth,xi
      double precision xij,ymin,ymax,depth(Ntot),rate00,frontdepth
      double precision traildepth,maxdis,ftol,P,V,xp(Ntot),yp(Ntot)
      double precision vx(Ntot),vy(Ntot),vth(Ntot),ax(Ntot),ay(Ntot)
      double precision ath(Ntot),bx(Ntot),by(Ntot),bth(Ntot),fx(Ntot)
      double precision fy(Ntot),fth(Ntot),inert(Ntot),typ_inert,torque
      double precision twist,twist_radians,t_double
      double precision bani
      integer N,seed,k,c(Ntot),i,Nforce,bound(Ntot),ntUL,ntDL
      integer forcelist(Ntot),keeppart(Ntot),prodskip,numsteps
      integer traillist(Ntot),countn(Ntot),nl(100,Ntot),dataskip
      integer Ntrail,Nsum,Nkeep,layerskip,Nprev
      character file1*80,file2*80,file3*80
      logical movie
      integer id(Ntot),idmax

      common /f2com/ width
      common /f3com/ alpha,Lx
      COMMON /f4com/ exp
      common /f8com/ forcelist,Nforce
      common /f10com/ bound
      common /f12com/ nl,countn
      
      ! read geometric parameters
      read(*,*) alpha0        ! aspect ratio at birth
      read(*,*) Lx            ! periodic box width
      read(*,*) D1            ! minor axis of cell

      ! read twist rate parameter
      read(*,*) twist         ! twist rate (degrees / cell cycle)

      ! read growth layer parameters
      read(*,*) layerwidth    ! width of bin that tracks furthest
                              !    forward cells for depth calc
      read(*,*) layerdepth    ! growth layer depth
      read(*,*) frontdepth    ! depth where cells considered to be in front
      
      ! read layer parameters for force calc
      read(*,*) propdepth     ! propagation layer depth - cells in force
                              !    calc, move when pushed
      read(*,*) bounddepth    ! boundary layer depth - cells in force
                              !    calc, do not move when pushed
      read(*,*) traildepth    ! trail layer depth - cells not in force calc

      ! read run parameters
      read(*,*) rate00        ! fundemental growth rate
      read(*,*) desync        ! desynchronization parameter
      read(*,*) seed          ! random seed
            
      ! read overdamped MD parameters
      read(*,*) numsteps      ! number of growth steps in simulation 
      read(*,*) dataskip      ! number of steps between outputting to screen
      read(*,*) prodskip      ! number of steps between outputting cell info
      read(*,*) layerskip     ! number of steps between calulating depth 
      read(*,*) dt            ! time-step
      read(*,*) b             ! damping coefficient
      read(*,*) bani          ! damping anisotropy (between -b and b)
      
      ! read output files
      read(*,*) movie         ! logical parameter for writing movie
      read(*,*) file1         ! file name

      ! additional parameters
      exp=2d0    ! 2 =  linear spring, 2.5 = Hertzian, >2.9 = RLJ
      width=D1   ! width of neighborlist 

      ! calculate # steps until division
      tdiv=dlog10(2d0)/dlog10(1d0+dt*rate00)

      ! calculate constant value of torque
      typ_inert=D1**2/32d0*dlog(4d0)/datanh(3d0/(8d0*alpha0**2+5d0))
      twist_radians=twist*pi/180d0
      t_double=dt*tdiv
      torque=(twist_radians*typ_inert)/(t_double*b)
      
      ! open movie file
      if(movie) then
         open(unit=1,file=TRIM(file1)//'.dat')
         open(unit=2,file=TRIM(file1)//'_adv.dat')
         open(unit=3,file=TRIM(file1)//'_lin.dat')
      endif

      ! initialize positions, velocities, etc
      call initialize(N,c,d,x,y,th,inert,rate0,rate,vx,vy,vth,ax,ay,ath,
     +     bx,by,bth,alpha0,rate00,desync,D1,torque,xp,yp,b,bani,seed,
     +     id,idmax)

      ! output initial configuration
      write(1,*) N
      do i=1,N                        
         write(1,'(6E20.12,I8)') x(i),y(i),th(i),d(i),
     +        d(i)*alpha(i),depth(i),c(i)
         write(2,'(6E20.12,3I16,3E20.12)') x(i),y(i),th(i),d(i),
     +        d(i)*alpha(i),depth(i),c(i),id(i),0,
     +        vx(i),vy(i),vth(i)
      enddo
         
      ! loop over time-steps 
      Nsum=N
      do k=1,numsteps
         
         ! calculate depth from front
         if(mod(k,layerskip).eq.0) then
            call calc_depth(N,Lx,d,x,y,layerwidth,depth)
            call remove_cells(N,Nforce,Ntrail,d,x,y,th,alpha,inert,
     +           c,rate,rate0,depth,vx,vy,vth,ax,ay,ath,bx,by,bth,
     +           forcelist,traillist,xp,yp,nl,countn,bound,id)
         endif

         ! grow and divide cells
         Nprev=N
         call grow_cells(N,Nsum,dt,x,y,th,D,c,rate,alpha0,rate0,
     +        seed,desync,traillist,depth,layerdepth,propdepth,
     +        bounddepth,vx,vy,vth,ax,ay,ath,bx,by,bth,id,idmax)
         if(N.gt.Nprev) then
            call makelist(N,x,y,D,xp,yp,countn,nl,alpha0)
         endif

         ! calculate inertia
         do i=1,N
            inert(i)=(1d0+alpha(i)**2)/16d0*d(i)**2
         enddo

         ! update neighborlist
         call checklist(N,x,y,xp,yp,maxdis)
         if(maxdis.gt.width) then ! fix - update if cell divides
            call makelist(N,x,y,D,xp,yp,countn,nl,alpha0)
         end if

         ! Gear precictor corrector
         call predict(dt,N,x,y,th,vx,vy,vth,ax,ay,ath,bx,by,bth)
         call force(N,x,y,th,D,V,P,fx,fy,fth,torque)
         call correct(dt,N,x,y,th,vx,vy,vth,ax,ay,ath,bx,by,bth,
     +        fx,fy,fth,inert,b,bani)

         ! output configuration
         if(mod(k,prodskip).eq.0) then
            write(1,*) N
            do i=1,N
               write(1,'(6E20.12,I8)') x(i),y(i),th(i),d(i),
     +              d(i)*alpha(i),depth(i),c(i)
               write(2,'(6E20.12,3I16,3E20.12)') x(i),y(i),th(i),d(i),
     +              d(i)*alpha(i),depth(i),c(i),id(i),k/prodskip,
     +              vx(i),vy(i),vth(i)
            enddo
         endif

         ! count # cells in growth layer
         ntUL=0  ! # in top layer
         ntDL=0  ! # on bottom layer
         do i=1,N
            if(depth(i).lt.layerdepth) then
               if(y(i).gt.0) then
                  ntUL=ntUL+1
               else
                  ntDL=ntDL+1
               endif
            endif
         enddo

         ! print summary to screen
         if(mod(k,dataskip).eq.0) then
            write(*,'(6I8,E20.10)')k,Nsum,N,Nforce,ntUL,ntDL,V/dble(N)
         endif
         
      enddo
      
      end ! end program


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine initialize(N,c,d,x,y,th,inert,rate0,rate,vx,vy,vth,ax,
     +     ay,ath,bx,by,bth,alpha0,rate00,desync,D1,torque,
     +     xp,yp,b,bani,seed,id,idmax)
      
      integer Ntot
      parameter(Ntot=2**18)      
      double precision pi
      parameter(pi=3.1415926535897932d0)
      double precision Lx,d(Ntot),x(Ntot),y(Ntot),th(Ntot),inert(Ntot)
      double precision rate0(Ntot),rate(Ntot),vx(Ntot),vy(Ntot),ax(Ntot)
      double precision ay(Ntot),bx(Ntot),by(Ntot),vth(Ntot),ath(Ntot)
      double precision bth(Ntot),fx(Ntot),fy(Ntot),fth(Ntot),alpha(Ntot)
      double precision ran2,alpha0,rate00,desync,xp(Ntot),yp(Ntot)
      double precision V,P,D1,torque,b,bani
      integer i,N,c(Ntot),bound(Ntot),Nforce,forcelist(Ntot)
      integer countn(Ntot),nl(100,Ntot),seed
      integer id(Ntot),idmax
      common /f3com/ alpha,Lx
      common /f8com/ forcelist,Nforce
      common /f10com/ bound
      common /f12com/ nl,countn
      
      ! random initial config for line of cells
      N=floor(Lx/2d0/D1/alpha0)
      do i=1,N     
         c(i)=i
         d(i)=D1
         x(i)=(dble(i)-dble(N+1)/2d0)*D1*2d0*alpha0
         y(i)=0d0
         th(i)=(ran2(seed)-0.5d0)*2d0*pi 
         alpha(i)=alpha0*(1d0+ran2(seed))
         inert(i)=(1d0+alpha(i)**2)/16d0*d(i)**2
         rate0(i)=rate00
         rate(i)=(1d0+(ran2(seed)-0.5d0)*desync)*rate0(i)
         id(i)=i
      enddo
      idmax=N
      
       ! set initial velocities     
      Nforce=N
      do i=1,N
         forcelist(i)=i
         bound(i)=0
      enddo
      call makelist(N,x,y,D,xp,yp,countn,nl,alpha0)
      call force(N,x,y,th,D,V,P,fx,fy,fth,torque)
      do i=1,N
         vx(i)=b*fx(i) + bani*(dcos(2*th(i))*fx(i)+dsin(2*th(i))*fy(i))
         vy(i)=b*fy(i) - bani*(-dsin(2*th(i))*fx(i)+dcos(2*th(i))*fy(i))
         vth(i)=b*fth(i)/inert(i)
         ax(i)=0d0
         ay(i)=0d0
         ath(i)=0d0
         bx(i)=0d0
         by(i)=0d0
         bth(i)=0d0         
      enddo
      
      end


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine calc_depth(N,Lx,d,x,y,layerwidth,depth)

      integer Ntot
      parameter(Ntot=2**18)      
      double precision Lx,x(Ntot),y(Ntot),xij,xi,ymax,ymin,layerwidth
      double precision depth(Ntot),depthU(Ntot),depthD(Ntot)
      double precision D(Ntot),dx,drsq
      integer N,i,j,bin,numbins,bini(Ntot),npart(100),ibin(100,Ntot)
      integer dbin,front(Ntot),numfrontD,numfrontU
      
      ! number of bins for front calc
      numbins=nint(Lx/layerwidth)

      ! calc distance to front
      do bin=1,numbins
         npart(bin)=0
      enddo
      
      ! calc distance to front
      offset=floor(Lx/2d0/layerwidth)+1
      do i=1,N
         xi=x(i)-dnint(x(i)/Lx)*Lx
         bin=floor(xi/layerwidth)+offset
         bini(i)=bin             ! x bin of cell i
         npart(bin)=npart(bin)+1 ! # cells in bin
         ibin(bin,npart(bin))=i  ! ibin = cell index in bin
      enddo

      ! calc distance to front      
      do i=1,N
         ymin=0d0
         ymax=0d0
         do dbin=-1,1
            bin=mod(numbins+bini(i)+dbin-1,numbins)+1
            do nj=1,npart(bin)                  
               j=ibin(bin,nj)
               xij=x(i)-x(j)
               xij=xij-dnint(xij/Lx)*Lx
               if(dabs(xij)<layerwidth) then
                  if(y(j).gt.ymax) then
                     ymax=y(j)
                  elseif(y(j).lt.ymin) then
                     ymin=y(j)
                  endif
               endif
            enddo
         enddo
         depth(i)=min(ymax-y(i),y(i)-ymin)       
      enddo

      ! assign cells near front to be at front
      numfront=0
      do i=1,N
         if(depth(i).lt.D(i)) then 
            numfront=numfront+1
            front(numfront)=i
            depth(i)=0d0
         endif
      enddo
      
      ! calc distance to nearest cell at front 
      do i=1,N
         do jj=1,numfront
            j=front(jj)
            dx=x(i)-x(j)
            if(dabs(dx).lt.depth(i)) then
               dy=y(i)-y(j)
               drsq=dx*dx+dy*dy
               if(drsq.lt.depth(i)**2) then
                  depth(i)=dsqrt(drsq)
               endif
            endif
         enddo
      enddo
      
      end


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine grow_cells(N,Nsum,dt,x,y,th,D,c,rate,alpha0,rate0,
     +     seed,desync,traillist,depth,layerdepth,propdepth,
     +     bounddepth,vx,vy,vth,ax,ay,ath,bx,by,bth,id,idmax)
      
      integer Ntot
      parameter(Ntot=2**18)      
      double precision alpha(Ntot),x(Ntot),y(Ntot),th(Ntot),d(Ntot)
      double precision rate0(Ntot),rate(Ntot),alpha0,dt,desync,ran2
      double precision depth(Ntot),layerdepth,propdepth,bounddepth
      double precision vx(Ntot),vy(Ntot),vth(Ntot),ax(Ntot),ay(Ntot)
      double precision ath(Ntot),bx(Ntot),by(Ntot),bth(Ntot),alphadiv
      integer i,N,Nsum,c(Ntot),bound(Ntot),Nforce,seed
      integer forcelist(Ntot),traillist(Ntot)
      integer id(Ntot),idmax
      common /f3com/ alpha,Lx
      common /f8com/ forcelist,Nforce
      common /f10com/ bound

      ! grow & divide cells
      Nforce=0
      Ntrail=0
      do i=1,N
            
         ! grow cells
         if(depth(i).lt.layerdepth) then
            alpha(i)=alpha(i)*(1d0+dt*rate(i))
         endif
            
         if(alpha(i).gt.2d0*alpha0) then
            alphadiv=alpha(i)
            
            ! divide into 2 - 1st assigned index N+1
            id(N)=idmax+1
            idmax=idmax+1
            write(*,'(I8)')idmax
            N=N+1
            Nsum=Nsum+1
            c(N)=c(i)
            D(N)=D(i)
            x(N)=x(i)+D(i)*alphadiv/4d0*dcos(th(i))
            y(N)=y(i)+D(i)*alphadiv/4d0*dsin(th(i))
            th(N)=th(i)
            rate(N)=(1d0+(ran2(seed)-0.5d0)*desync)*rate0(i) 
            alpha(N)=alpha0               
            rate0(N)=rate0(i)
            vx(N)=vx(i)
            vy(N)=vy(i)
            vth(N)=vth(i)               
            ax(N)=ax(i)
            ay(N)=ay(i)
            ath(N)=ath(i)               
            bx(N)=bx(i)
            by(N)=by(i)
            bth(N)=bth(i)

            write(3,'(2I8,2E20.12)')N,id(i),x(N),y(N)
            
            ! divide into 2 - 1st assigned index i
            x(i)=x(i)-D(i)*alphadiv/4d0*dcos(th(i))
            y(i)=y(i)-D(i)*alphadiv/4d0*dsin(th(i))
            rate(i)=(1d0+(ran2(seed)-0.5d0)*desync)*rate0(i)
            alpha(i)=alpha0

            ! add particle N to propagation &/or force calc lists
            forcelist(Nforce+1)=i
            forcelist(Nforce+2)=N
            Nforce=Nforce+2
            bound(i)=0
            bound(N)=0
         else
            if(depth(i).lt.layerdepth+propdepth+bounddepth) then
               forcelist(Nforce+1)=i
               Nforce=Nforce+1
               if (depth(i).gt.layerdepth+propdepth) then
                  bound(i)=1
               else
                  bound(i)=0
               endif                  
            else if(depth(i).lt.layerdepth+propdepth+
     +              bounddepth+traildepth) then
               traillist(Ntrail+1)=i
               Ntrail=Ntrail+1
            endif                  
         endif            
      enddo

      end

      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine remove_cells(N,Nforce,Ntrail,d,x,y,th,alpha,
     +     inert,c,rate,rate0,depth,vx,vy,vth,ax,ay,ath,
     +     bx,by,bth,forcelist,traillist,xp,yp,nl,countn,bound,id)
 
      integer Ntot
      parameter(Ntot=2**18)      
      double precision x(Ntot),y(Ntot),th(Ntot),alpha(Ntot),depth(Ntot)
      Double precision vx(Ntot),vy(Ntot),vth(Ntot),rate(Ntot),d(Ntot)
      double precision ax(Ntot),ay(Ntot),ath(Ntot),rate0(Ntot)
      double precision bx(Ntot),by(Ntot),bth(Ntot),inert(Ntot)
      double precision xp(Ntot),yp(Ntot)
      integer id(Ntot)
      integer N,Nforce,Ntrail,Nkeep,c(Ntot),keeppart(Ntot)
      integer forcelist(Ntot),traillist(Ntot),renumber(Ntot)
      integer countn(Ntot),nl(100,Ntot),ikeep,jkeep,bound(Ntot)

      ! make list of cells to keep 
      do i=1,N         
         keeppart(i)=0
      enddo
      do i=1,Nforce            
         keeppart(forcelist(i))=1
      enddo
      do i=1,Ntrail
         keeppart(traillist(i))=1
      enddo

      ! remove cells & shift values to lower indices
      ikeep=0
      do i=1,N
         if(keeppart(i).eq.1) then
            ikeep=ikeep+1
            id(ikeep)=id(i)
            x(ikeep)=x(i)
            y(ikeep)=y(i)
            th(ikeep)=th(i)
            d(ikeep)=d(i)
            alpha(ikeep)=alpha(i)
            inert(ikeep)=inert(i)
            c(ikeep)=c(i)
            rate(ikeep)=rate(i)
            rate0(ikeep)=rate0(i)   
            depth(ikeep)=depth(i)
            vx(ikeep)=vx(i)
            ax(ikeep)=ax(i)
            bx(ikeep)=bx(i)
            vy(ikeep)=vy(i)
            ay(ikeep)=ay(i)
            by(ikeep)=by(i)
            vth(ikeep)=vth(i)
            ath(ikeep)=ath(i)
            bth(ikeep)=bth(i)
            xp(ikeep)=xp(i)
            yp(ikeep)=yp(i)
            bound(ikeep)=bound(i)

            ! index map for subsequent lists
            renumber(i)=ikeep
         endif
      enddo
      Nkeep=ikeep

      ! replace lists with new indices
      do i=1,Nforce
         forcelist(i)=renumber(forcelist(i))
      enddo
      do i=1,Ntrail
         traillist(i)=renumber(traillist(i)) 
      enddo

      ! remove cells from neighborlist
      ikeep=0
      do i=1,N	
         if(keeppart(i).eq.1) then
            ikeep=ikeep+1
            jkeep=0
            do j=1,countn(i)
               if(keeppart(nl(j,i)).eq.1) then
                  jkeep=jkeep+1
                  nl(jkeep,ikeep)=renumber(nl(j,i))
               endif
            enddo
            countn(ikeep)=jkeep
         endif
      enddo

      N=Nkeep         
 
      end
      
         
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      double precision function sigmasq(rijsq,xij,yij,thi,thj,
     +     di,dj,alphai,alphaj)
      
      implicit none

      double precision xij,yij,thi,thj,di,dj
      double precision rmui,rmuj,muij
      double precision rij, li, lj, alphai, alphaj
      double precision chi, alp, sig_0sq, rijsq, sigsq
      double precision n1, n2, d
      common /f11com/ rmui,rmuj,muij,chi,alp,sig_0sq,n1,n2,d,sigsq

      if(alphai.eq.1d0.and.alphaj.eq.1d0) then
         sigmasq = (di+dj)*(di+dj)/4d0
      else
         li = alphai*di
         lj = alphaj*dj
         
         rmui=xij*dcos(thi)+yij*dsin(thi)
         rmuj=xij*dcos(thj)+yij*dsin(thj)
         muij=dcos(thi-thj)

         chi = dsqrt((li*li-di*di)*(lj*lj-dj*dj)/
     +      ((lj*lj+di*di)*(li*li+dj*dj)))
         alp = dsqrt(dsqrt((li*li-di*di)*(lj*lj+di*di)/
     +        ((lj*lj-dj*dj)*(li*li+dj*dj))))
         sig_0sq = (di*di+dj*dj)/2d0

         n1=alp*rmui+rmuj/alp
         n2=alp*rmui-rmuj/alp
         d=chi*muij

         sigmasq = sig_0sq
     +        /(1d0-0.5d0*chi*(n1*n1/(1d0+d)+n2*n2/(1d0-d))/rijsq)
         sigsq=sigmasq
      end if

      return
      end

      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine sigma(rijsq,xij,yij,thi,thj,di,dj,
     +     sig,dlnsig,fthi,fthj,alphai,alphaj)
      
      double precision xij,yij,thi,thj,di,dj
      double precision rmui,rmuj,muij
      double precision rij, li, lj, alphai, alphaj
      double precision chi, alp, sig_0
      double precision d,n1,n2,n3,n4,n5,n6,n7,n8
      double precision sig,dlnsig,pre
      double precision s2,rmui2,rmuj2,muij2,fthi,fthj
      double precision sig_0sq,sigsq,rijsq
      common /f11com/ rmui,rmuj,muij,chi,alp,sig_0sq,n1,n3,d,sigsq

      if(alphai.eq.1d0.and.alphaj.eq.1d0) then
         sig = 0.5d0*(di+dj)
         dlnsig = 0d0
         fthi = 0d0
         fthj = 0d0
      else
         sig=dsqrt(sigsq)

         n2=n1/(1d0+d)
         n4=n3/(1d0-d)

         s2 = sigsq/sig_0sq
         rmui2 = (yij*dcos(thi)-xij*dsin(thi))
         rmuj2 = (yij*dcos(thj)-xij*dsin(thj))
         n5 = alp*rmui2+rmuj2/alp
         n6 = alp*rmui2-rmuj2/alp
         pre=0.5d0*chi*s2
         dlnsig = -pre*(n2*n5+n4*n6) /rijsq

         muij2 = dsin(thi-thj)
         n8 = chi*muij2*(n2*n2-n4*n4)
         fthi = 0.5d0*pre*(2d0*alp*rmui2*(n2+n4)+n8)/rijsq
         fthj = 0.5d0*pre*(2d0/alp*rmuj2*(n2-n4)-n8)/rijsq   
      end if

      return
      end

      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!!!!!!!!!!!!!    predicts new positions and velocities    !!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
      subroutine predict(dt,N,x,y,th,vx,vy,vth,ax,ay,ath,bx,by,bth)     
      parameter(Ntot=2**18)
      integer N,i,ensemble
      double precision dt,x(Ntot),y(Ntot),vx(Ntot),vy(Ntot)
      double precision ax(Ntot),ay(Ntot),bx(Ntot),by(Ntot),v
      double precision th(Ntot),vth(Ntot),ath(Ntot),bth(Ntot),c1,c2,c3,L

      c1 = dt
      c2 = c1*dt/2d0
      c3 = c2*dt/3d0

      do i=1,N
         x(i) = x(i) + c1*vx(i) + c2*ax(i) + c3*bx(i)
         y(i) = y(i) + c1*vy(i) + c2*ay(i) + c3*by(i)
         th(i) = th(i) + c1*vth(i) + c2*ath(i) + c3*bth(i)         
         vx(i) = vx(i) + c1*ax(i) + c2*bx(i)
         vy(i) = vy(i) + c1*ay(i) + c2*by(i)     
         vth(i) = vth(i) + c1*ath(i) + c2*bth(i)     
         ax(i) = ax(i) + c1*bx(i)
         ay(i) = ay(i) + c1*by(i)
         ath(i) = ath(i) + c1*bth(i)
      enddo

      end


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!!!!!!!!!!!!!!!!!!!!!!!   corrects prediction   !!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
      subroutine correct(dt,N,x,y,th,vx,vy,vth,ax,ay,ath,bx,by,bth,
     +     fx,fy,fth,inert,b,bani)
      parameter(Ntot=2**18)
      integer i,N,ensemble
      double precision dt,x(Ntot),y(Ntot),vx(Ntot),vy(Ntot),ax(Ntot)
      double precision ay(Ntot),bx(Ntot),by(Ntot),th(Ntot),vth(Ntot)
      double precision ath(Ntot),bth(Ntot),fx(Ntot),fy(Ntot),fth(Ntot)
      double precision Lx,c1,c2,c3,gear0,gear2,gear3,cg0,cg2,cg3
      double precision vxi,vyi,vthi,pvxi,pvyi,pvthi,corrx,corry,b,bani
      double precision corrpx,corrpy,corr,corrth,corrpth,inert(Ntot)

      gear0 = 3d0/8d0
      gear2 = 3d0/4d0
      gear3 = 1d0/6d0

      c1 = dt
      c2 = c1*dt/2d0
      c3 = c2*dt/3d0

      cg0 = gear0*c1
      cg2 = gear2*c1/c2
      cg3 = gear3*c1/c3

      do i=1,N
         vxi = b*fx(i) + bani*(dcos(2*th(i))*fx(i)+dsin(2*th(i))*fy(i))
         vyi = b*fy(i) - bani*(-dsin(2*th(i))*fx(i)+dcos(2*th(i))*fy(i))
         vthi = b*fth(i)/inert(i)
         corrx = vxi - vx(i)
         corry = vyi - vy(i)
         corrth = vthi - vth(i)
         x(i) = x(i) + cg0*corrx
         y(i) = y(i) + cg0*corry
         th(i) = th(i) + cg0*corrth        
         vx(i) = vxi
         vy(i) = vyi
         vth(i) = vthi
         ax(i) = ax(i) + cg2*corrx
         ay(i) = ay(i) + cg2*corry
         ath(i) = ath(i) + cg2*corrth
         bx(i) = bx(i) + cg3*corrx
         by(i) = by(i) + cg3*corry
         bth(i) = bth(i) + cg3*corrth         
      enddo

      end

            
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!           force            !!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine force(N,x,y,th,D,V,P,fx,fy,fth,torque) ! dimer

      parameter(Ntot=2**18)
      double precision pi
      parameter(pi=3.1415926535897932d0)
      double precision x(Ntot),y(Ntot),th(Ntot),ep,D(Ntot),dij
      double precision fx(Ntot),fy(Ntot),fth(Ntot),rij,xij,yij,fr,exp
      double precision dij_up,alpha(Ntot),LJ,fc,ft,f_x,f_y,scale(Ntot)
      double precision fthi,fthj,fth_c,rijsq,dijsq_up,c(Ntot),s(Ntot),dd
      double precision dr(Ntot,2),xa(Ntot,2),ya(Ntot,2),dk(Ntot,2),V
      double precision di_up(Ntot),di1j1,Vij,Lx,ap,Pij,P,dijsq
      double precision outang,num,dem1,dem2,yang,soutang,coutang
      double precision soutang2,coutang2,xdiagu,xdiagb,FwallR,PwallR
      double precision FwallT,PwallT,FwallB,PwallB,FwallL,PwallL,torque
      double precision sigmasq
      integer countn(Ntot),nl(100,Ntot),N,ki,kj,jj,up,down
      integer forcelist(Ntot),Nforce,bound(Ntot)

      common /f8com/ forcelist,Nforce
      common /f10com/ bound
      common /f4com/ exp
      common /f12com/ nl,countn
      common /f3com/ alpha,Lx ! aspect ratio
      
      do i=1,N
         fx(i)=0d0
         fy(i)=0d0
         if(bound(i).eq.1) then
            fth(i)=0d0
         else
            fth(i)=torque
         endif
      enddo
      V=0d0
      P=0d0
      
      ! inter-particle interactions
      do ii=1,Nforce
         i=forcelist(ii)
         if(countn(i).ge.1) then
            do j=1,countn(i)
               if(bound(i).eq.0.or.bound(nl(j,i)).eq.0) then
                  dijsq_up=(alpha(i)**2*D(i)**2+
     +                 alpha(nl(j,i))**2*D(nl(j,i))**2)/2d0
                  xij=x(i)-x(nl(j,i))               
                  xij=xij-dnint(xij/Lx)*Lx
                  if(xij*xij.lt.dijsq_up) then 
                     yij=y(i)-y(nl(j,i))                  
                     rijsq=xij**2+yij**2
                     if(rijsq.lt.dijsq_up) then 
                        dijsq=sigmasq(rijsq,xij,yij,th(i), th(nl(j,i)),
     +                       D(i),D(nl(j,i)),alpha(i),alpha(nl(j,i)))   
                        if(rijsq.lt.dijsq) then
                           call sigma(rijsq,xij,yij,th(i),
     +                          th(nl(j,i)),D(i),
     +                          D(nl(j,i)),dij,ft,fthi,fthj,
     +                          alpha(i),alpha(nl(j,i)))
                           rij=dsqrt(rijsq)
                           if(rij.lt.dij) then
                              if(exp .gt. 2.9) then
                                 LJ = (dij/rij)*(dij/rij)
                                 LJ = LJ*LJ*LJ
                                 fc = 1d0/rij*LJ*(LJ-1d0)
                                 V=V+(LJ-1d0)**2
                              else
                                 fc=(1d0-rij/dij)**(exp-1d0)/dij 
                                 V=V+(1-rij/dij)**exp/exp
                              endif
                              fr=-fc/rij
                              fth_c=rij*fc
                              f_x=fr*(xij+yij*ft)
                              f_y=fr*(yij-xij*ft)                   
                              if(bound(i).eq.0) then
                                 fx(i)=fx(i)-f_x
                                 fy(i)=fy(i)-f_y
                                 fth(i)=fth(i)-fthi*fth_c
                              endif                           
                              if(bound(nl(j,i)).eq.0) then
                                 fx(nl(j,i))=fx(nl(j,i))+f_x
                                 fy(nl(j,i))=fy(nl(j,i))+f_y
                                 fth(nl(j,i))=fth(nl(j,i))-fthj*fth_c
                              endif
                           end if
                        endif
                     end if
                  end if
               endif
            enddo
         end if
      enddo
      
      if(exp .gt. 2.9) then
         do i=1, N
            fx(i)=fx(i)/6d0
            fy(i)=fy(i)/6d0
            fth(i)=fth(i)/6d0 
         enddo
      endif
            
      return							
      end


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine checklist(N,x,y,xp,yp,maxdis)
      parameter(Ntot = 2**18)
      double precision maxdis,x(Ntot),y(Ntot),xp(Ntot),yp(Ntot)
      integer N

      maxdis=0d0
      do i=1,N
         maxdis=max(dabs(x(i)-xp(i)),maxdis)
         maxdis=max(dabs(y(i)-yp(i)),maxdis)
      enddo
      maxdis=2d0*dsqrt(2d0*maxdis*maxdis)

      return
      end


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine makelist(N,x,y,D,xp,yp,countn,nl,alpha0)

      parameter(Ntot = 2**18)
      double precision x(Ntot),y(Ntot),xp(Ntot),yp(Ntot),D(Ntot),Lx
      double precision xij,yij,rij,dij,rijsq,alpha(Ntot),alpha0,width
      integer countn(Ntot),nl(100,Ntot),N,ii,jj,i,j
      integer forcelist(Ntot),Nforce
      common /f2com/ width
      common /f3com/ alpha,Lx ! aspect ratio
      common /f8com/ forcelist,Nforce

      do i=1,N
         countn(i)=0
      enddo

      do ii=1,Nforce-1
         i=forcelist(ii)
         do jj=ii+1,Nforce
            j=forcelist(jj)
 
            xij=x(i)-x(j)
            xij=xij-dnint(xij/Lx)*Lx
            yij=y(i)-y(j)
            rijsq=xij*xij+yij*yij
            dij=2d0*alpha0*dsqrt((D(i)**2+D(j)**2)/2d0) ! dij of major axes
            if(rijsq.lt.(dij+width)**2) then
               countn(i)=countn(i)+1
               nl(countn(i),i)=j
            end if
         enddo
      enddo
      
      do i=1,N
         xp(i)=x(i)
         yp(i)=y(i)
      enddo

      return
      end

      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      FUNCTION ran2(idum)
      INTEGER idum,IM1,IM2,IMM1,IA1,IA2,IQ1,IQ2,IR1,IR2,NTAB,NDIV
      double precision ran2,AM,EPS,RNMX
      PARAMETER (IM1=2147483563,IM2=2147483399,AM=1./IM1,IMM1=IM1-1,
     *IA1=40014,IA2=40692,IQ1=53668,IQ2=52774,IR1=12211,IR2=3791,
     *NTAB=32,NDIV=1+IMM1/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
      INTEGER idum2,j,k,iv(NTAB),iy
      SAVE iv,iy,idum2
      DATA idum2/123456789/, iv/NTAB*0/, iy/0/
      if (idum.le.0) then
        idum=max(-idum,1)
        idum2=idum
        do 11 j=NTAB+8,1,-1
          k=idum/IQ1
          idum=IA1*(idum-k*IQ1)-k*IR1
          if (idum.lt.0) idum=idum+IM1
          if (j.le.NTAB) iv(j)=idum
11      continue
        iy=iv(1)
      endif
      k=idum/IQ1
      idum=IA1*(idum-k*IQ1)-k*IR1
      if (idum.lt.0) idum=idum+IM1
      k=idum2/IQ2
      idum2=IA2*(idum2-k*IQ2)-k*IR2
      if (idum2.lt.0) idum2=idum2+IM2
      j=1+iy/NDIV
      iy=iv(j)-idum2
      iv(j)=idum
      if(iy.lt.1)iy=iy+IMM1
      ran2=min(AM*iy,RNMX)
      return
      END
