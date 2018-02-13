      subroutine stokesSLP(s1,s2,xs,ys,ns,xt,yt,nt,riequal,riprec,
     $      u1,u2)
c     Evaluates the single-layer potential for Stokes' equation
c     (s1,s2) is the source strength of length ns
c     xs and ys are the source points and xt and yt are the target 
c     locations of length ns
c     iequal is 1 if sources .eq. targest
c     iequal is 0 if sources .ne. targets
c     iprec controls the accuracy of the FMM
c         -2 => tolerance =.5d0
c         -1 => tolerance =.5d-1
c          0 => tolerance =.5d-2
c          1 => tolerance =.5d-3
c          2 => tolerance =.5d-6
c          3 => tolerance =.5d-9
c          4 => tolerance =.5d-12
c          5 => tolerance =.5d-15
c     (u1,u2) are the two components of the velocity field
      implicit real*8 (a-h,o-z)

      integer error
      real *8 s1(ns),s2(ns)
c     strength of Stokes SLP
      real *8 xs(ns),ys(ns)
c     location of the source points
      real *8 xt(nt),yt(nt)
c     location of the target points
      real *8 u1(nt),u2(nt)
c     x and y components of the velocity field
      
      real *8, allocatable :: sources(:,:)
c     location of sources
      real *8, allocatable :: targets(:,:)
c     location of targets
      real *8, allocatable :: charges(:)
c     charge strength of single-layer potential term
      real *8, allocatable :: dipstr(:),dipvec(:,:)
c     charge strength and direction of the 
c     double-layer potential term
      real *8, allocatable :: pot1(:),pot2(:)
c     room for two components that have to be summed
c     to form the velocity field

      allocate(targets(2,nt),stat=error)
c      if (error .ne. 0) then
c        print*,'ERROR WITH ALLOCATING TARGETS'
c      endif
      allocate(sources(2,ns),stat=error)
c      if (error .ne. 0) then
c        print*,'ERROR WITH ALLOCATING SOURCES'
c      endif
      allocate(charges(ns),stat=error)
c      if (error .ne. 0) then
c        print*,'ERROR WITH ALLOCATING CHARGES'
c      endif
      allocate(dipstr(ns),stat=error)
c      if (error .ne. 0) then
c        print*,'ERROR WITH ALLOCATING DIPSTR'
c      endif
      allocate(dipvec(2,ns),stat=error)
c      if (error .ne. 0) then
c        print*,'ERROR WITH ALLOCATING DIPVEC'
c      endif
      allocate(pot1(nt),stat=error)
c      if (error .ne. 0) then
c        print*,'ERROR WITH ALLOCATING POT1'
c      endif
      allocate(pot2(nt),stat=error)
c      if (error .ne. 0) then
c        print*,'ERROR WITH ALLOCATING POT2'
c      endif
c     allocate memory for temporary variables

      twopi = 8.d0*datan(1.d0)

c      iprec = 4 ! enough for 12 digits
      iequal = int(riequal)
      iprec = int(riprec)

      if (iequal .eq. 1) then
        ifpot = 1 ! need the potential since sources .eq. targets
      else
        ifpot = 0 ! don't need the potential since sources ~= targets
      endif
      ifgrad = 0 ! don't need the gradient
      ifhess = 0 ! don't need the Hessian
      ifpottarg = 1 ! need the potential at the targets
      ifgradtarg = 0 ! don't need the gradient at the targets
      ifhesstarg = 0 ! don't need the Hessian at the targets
c     set flags for what output values we require

      do i=1,ns
        sources(1,i) = xs(i)
        sources(2,i) = ys(i)
      enddo
c     set source locations
      do i=1,nt
        targets(1,i) = xt(i)
        targets(2,i) = yt(i)
      enddo
c     set target locations

c     START OF FORMING FIRST COMPONENET OF VELCOTIY
      ifcharge = 1 ! need a charge component
      ifdipole = 1 ! need a dipole componenet

      do i=1,ns
        charges(i) = -s1(i)
        dipstr(i) = xs(i)
        dipvec(1,i) = s1(i)
        dipvec(2,i) = s2(i)
      enddo

      if (iequal .eq. 1) then
        call rfmm2dpart(ierr,iprec,ns,sources,
     $     ifcharge,charges,ifdipole,dipstr,dipvec,
     $     ifpot,pot1,ifgrad,grad,ifhess,hess)
      else
        call rfmm2dparttarg(ierr,iprec,ns,sources,
     $     ifcharge,charges,ifdipole,dipstr,dipvec,
     $     ifpot,pot,ifgrad,grad,ifhess,hess,
     $     nt,targets,ifpottarg,pot1,ifgradtarg,gradtarg,
     $     ifhesstarg,hesstarg)
      endif
c     compute the first component in the Stokes SLP

      ifcharge = 0 ! don't need a charge component
      ifdipole = 1 ! need a dipole component

      do i=1,ns
        dipstr(i) = -1.d0
      enddo

      if (iequal .eq. 1) then
        call rfmm2dpart(ierr,iprec,ns,sources,
     $     ifcharge,charges,ifdipole,dipstr,dipvec,
     $     ifpot,pot2,ifgrad,grad,ifhess,hess)
      else
        call rfmm2dparttarg(ierr,iprec,ns,sources,
     $     ifcharge,charges,ifdipole,dipstr,dipvec,
     $     ifpot,temp1,ifgrad,temp1,ifhess,temp1,
     $     nt,targets,ifpottarg,pot2,ifgradtarg,temp1,
     $     ifhesstarg,temp1)
      endif
c     compute the first component in the Stokes SLP

      do i=1,nt
        u1(i) = pot1(i) + xt(i)*pot2(i)
        u1(i) = 5.d-1*u1(i)/twopi
      enddo
c     END OF FORMING FIRST COMPONENT OF VELOCITY


c     START OF FORMING SECOND COMPONENET OF VELCOTIY
      ifcharge = 1 ! need a charge component
      ifdipole = 1 ! need a dipole componenet

      do i=1,ns
        charges(i) = -s2(i)
        dipstr(i) = ys(i)
        dipvec(1,i) = s1(i)
        dipvec(2,i) = s2(i)
      enddo

      if (iequal .eq. 1) then
        call rfmm2dpart(ierr,iprec,ns,sources,
     $     ifcharge,charges,ifdipole,dipstr,dipvec,
     $     ifpot,pot1,ifgrad,grad,ifhess,hess)
      else
        call rfmm2dparttarg(ierr,iprec,ns,sources,
     $     ifcharge,charges,ifdipole,dipstr,dipvec,
     $     ifpot,pot,ifgrad,grad,ifhess,hess,
     $     nt,targets,ifpottarg,pot1,ifgradtarg,gradtarg,
     $     ifhesstarg,hesstarg)
      endif
c     compute the first component in the Stokes SLP

      ifcharge = 0 ! don't need a charge component
      ifdipole = 1 ! need a dipole componenet

      do i=1,ns
        dipstr(i) = -1.d0
      enddo

      if (iequal .eq. 1) then
        call rfmm2dpart(ierr,iprec,ns,sources,
     $     ifcharge,charges,ifdipole,dipstr,dipvec,
     $     ifpot,pot2,ifgrad,grad,ifhess,hess)
      else
        call rfmm2dparttarg(ierr,iprec,ns,sources,
     $     ifcharge,charges,ifdipole,dipstr,dipvec,
     $     ifpot,temp1,ifgrad,temp1,ifhess,temp1,
     $     nt,targets,ifpottarg,pot2,ifgradtarg,temp1,
     $     ifhesstarg,temp1)
      endif
c     compute the second component in the Stokes SLP

      do i=1,nt
        u2(i) = pot1(i) + yt(i)*pot2(i)
        u2(i) = 5.d-1*u2(i)/twopi
      enddo
c     END OF FORMING SECOND COMPONENET OF VELCOTIY
          


      end
 
