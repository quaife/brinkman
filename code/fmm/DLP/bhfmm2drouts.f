cc Copyright (C) 2010-2011: Leslie Greengard and Zydrunas Gimbutas
cc Contact: greengard@cims.nyu.edu
cc 
cc This program is free software; you can redistribute it and/or modify 
cc it under the terms of the GNU General Public License as published by 
cc the Free Software Foundation; either version 2 of the License, or 
cc (at your option) any later version.  This program is distributed in 
cc the hope that it will be useful, but WITHOUT ANY WARRANTY; without 
cc even the implied warranty of MERCHANTABILITY or FITNESS FOR A 
cc PARTICULAR PURPOSE.  See the GNU General Public License for more 
cc details. You should have received a copy of the GNU General Public 
cc License along with this program; 
cc if not, see <http://www.gnu.org/licenses/>.
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c    $Date$
c    $Revision$
c
c       
c
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c        this is the end of the debugging code and the beginning 
c        of the routines for Helmholtz FMM in R^2
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c
        subroutine bhfmm2dparttree(ier,iprec,
     $     nsource,source,ntarget,target,
     $     nbox,epsfmm,iisource,iitarget,iwlists,lwlists,
     $     nboxes,laddr,nlev,center,size,
     $     w,lw,lused7)
        implicit real *8 (a-h,o-z)
c       
c       Helmholtz FMM in R^2: build the quad-tree
c
c     INPUT PARAMETERS:
c
c       nsource: integer:  number of sources
c       source: real *8 (2,n):  source locations
c       w: real *8 (lw): workspace
c       lw:  length of workspace
c
c     OUTPUT PARAMETERS:
c
c       ier   =  error return code
c       lused7 = the amount of workspace used
c
c
        real *8 source(2,*),target(2,*)
c       
        real *8 center(2)
c       
        integer laddr(2,200)
        integer box(15)
        real *8 center0(2),corners0(2,4)
c       
        integer box1(15)
        real *8 center1(2),corners1(2,4)
c
        real *8 w(*)
c       
        ier=0
c       
        done=1
        pi=4*atan(done)
c       
        lused7=0
        ifprint=0
c       
        iisource=1
        lused7=lused7+nsource
        if (ifprint.eq.1) call prinf('lused7=*',lused7,1)
        if (lused7 .ge. lw) ier=128
        if( ier .ne. 0 ) return
c
        iitarget=iisource+nsource
        lused7=lused7+ntarget
        if (ifprint.eq.1) call prinf('lused7=*',lused7,1)
        if (lused7 .ge. lw) ier=128
        if( ier .ne. 0 ) return
c
        iwlists=iisource+lused7+10
c
c       ... construct the adaptive FMM quad-tree structure
c
c        call d2tstrcr(ier,source,nsource,nbox,
c     $     nboxes,w(iisource),laddr,nlev,center,size,
c     $     target,ntarget,w(iitarget),w(iwlists),lw-lused7,lused)
        ifempty=0
        minlevel=0
        maxlevel=30
        call d2tstrcrem(ier,source,nsource,nbox,
     $     nboxes,w(iisource),laddr,nlev,center,size,
     $     target,ntarget,w(iitarget),w(iwlists),lw-lused7,lused,
     $     ifempty,minlevel,maxlevel)
c
        if( ier .ne. 0 ) return
c
        lwlists=lused
        lused7=lused7+lwlists
        if (lused7 .ge. lw) ier=128
        if( ier .ne. 0 ) return
c       
        if (ifprint.eq.1) 
     $       call prin2('after d2tstrcr, center=*',center,2)
        if (ifprint.eq.1) 
     $       call prin2('after d2tstrcr, size=*',size,1)
        if (ifprint.eq.1) 
     $       call prinf('after d2tstrcr, nlev=*',nlev,1)
        if (ifprint.eq.1) 
     $       call prinf('after d2tstrcr, nbox=*',nbox,1)
        if (ifprint.eq.1) 
     $       call prinf('after d2tstrcr, laddr=*',laddr,2*(nlev+1))
c
ccc        call prinf('after d2tstrcr, isource=*',w(iisource),nsource)
ccc        call prinf('after d2tstrcr, itarget=*',w(iitarget),ntarget)
c
ccc        call d2tprint(w(iwlists),lwlists)
c
c       ... optional, plot the oct-tree in gnuplot compatible format
c
        ifplot = 0
        if (ifplot .eq. 1 .and. nsource .lt. 10000 ) then
c
c       ... plot the boxes
c
        iw=51
        call plot_box2d(iw,center,size)
c       
        itag=1
        iw=52
        call plot_label2d(iw,center,size,itag,itag)
c
        iw=60
        call plot_points2d(iw,source,nsource)
c       
        iw=63
        call plot_points2d(iw,target,ntarget)
c       
        do ibox=1,nboxes
           call d2tgetb(ier,ibox,box,center0,corners0,w(iwlists))
           level=box(1)
           size0=size/2**level
c
           iw=61
           call plot_box2d(iw,center0,size0)
c       
           itag=ibox
           iw=62
           call plot_label2d(iw,center0,size0,itag,itag)
        enddo  
c      
        endif
c
        return
        end
c
c
c
c
c
        subroutine bhfmm2d_list2
     $     (bsize,nlev,laddr,scale,nterms,rmlexp,iaddr,epsfmm,
     $     timeinfo,wlists,mptemp1,lmptemp1,mptemp2,lmptemp2,
     $     mptemp3,lmptemp3,mptemp4,lmptemp4,mptemp5,lmptemp5,
     $     ifprune_list2)
   
        implicit real *8 (a-h,o-z)
c
        integer iaddr(10,*),laddr(2,*),nterms(0:*)
        real *8 rmlexp(*),scale(0:*)
        integer itable(-3:3,-3:3)
c
        integer list(10 000)
c
        integer box(15)
        real *8 bsize(0:200)
        real *8 center0(2),corners0(2,4)
c       
        integer box1(15)
        real *8 center1(2),corners1(2,4)
c
        real *8 wlists(*)
        complex *16 mptemp1(lmptemp1)
        complex *16 mptemp2(lmptemp2)
        complex *16 mptemp3(lmptemp3)
        complex *16 mptemp4(lmptemp4)
        complex *16 mptemp5(lmptemp5)
        complex *16 vtemp,gatemp,gaatemp
c
        real *8 timeinfo(10)
c
        real *8, allocatable :: carray(:,:)
c               
c
        ldc = 100
        allocate(carray(0:ldc,0:ldc))
        call bh2d_init_carray(carray,ldc)
c
c     ifprint is an internal information printing flag. 
c     Suppressed if ifprint=0.
c     Prints timing breakdown and other things if ifprint=1.
c     Prints timing breakdown, list information, and other things if ifprint=2.
c       
        ifprint=1
c
        if (ifprint .ge. 1) 
     $     call prinf('=== STEP 3 (merge mp) ===*',i,0)
c        t1=second()
c        t1=omp_get_wtime() 
c       ... step 3, merge all multipole expansions
c       
        do 2300 ilev=nlev,3,-1
!$OMP PARALLEL DO DEFAULT(SHARED)
!$OMP$PRIVATE(ibox,box,center0,corners0,level0,level,npts,nkids,radius)
!$OMP$PRIVATE(jbox,box1,center1,corners1,level1)
!$OMP$PRIVATE(mptemp1,lused,ier,i,j,vtemp,gatemp,cd)
!$OMP$PRIVATE(mptemp2,mptemp3,mptemp4,mptemp5,gaatemp)
ccc!$OMP$SCHEDULE(DYNAMIC)
ccc!$OMP$NUM_THREADS(4) 
           do 2200 ibox=laddr(1,ilev),laddr(1,ilev)+laddr(2,ilev)-1
              call d2tgetb(ier,ibox,box,center0,corners0,wlists)
              call d2tnkids(box,nkids)
c       ... prune all sourceless boxes
              if( box(10) .eq. 0 ) goto 2200
              if (nkids .ne. 0) then
                 level0=box(1)
                 if( level0 .ge. 2 ) then
                    radius = (corners0(1,1) - center0(1))**2
                    radius = radius + (corners0(2,1) - center0(2))**2
                    radius = sqrt(radius)
c       
                    if( ifprint .ge. 2 ) then
                       call prin2('radius=*',radius,1)
                       call prinf('ibox=*',ibox,1)
                       call prinf('box=*',box,15)
                       call prinf('nkids=*',nkids,1)
                    endif
c
c       ... merge multipole expansions of the kids 
c
                    call bh2dzero(rmlexp(iaddr(1,ibox)),nterms(level0))
                    call bh2dzero(rmlexp(iaddr(2,ibox)),nterms(level0))
                    call bh2dzero(rmlexp(iaddr(3,ibox)),nterms(level0))
                    call bh2dzero(rmlexp(iaddr(4,ibox)),nterms(level0))
                    call bh2dzero(rmlexp(iaddr(5,ibox)),nterms(level0))
                    if (ifprint .ge. 2) then
                       call prin2('center0=*',center0,2)
                    endif
c
                    do 2100 i = 1,4
                       jbox = box(4+i)
                       if (jbox.eq.0) goto 2100
                       call d2tgetb(ier,jbox,box1,center1,corners1,
     $                 wlists)
                       if (ifprint .ge. 2) then
                          call prinf('jbox=*',jbox,1)
                          call prin2('center1=*',center1,2)
                       endif
                       level1=box1(1)
                       if(nterms(level0)+nterms(level1).gt.95) then
                          call bh2dmpmp(scale(level1),center1,
     1                    rmlexp(iaddr(1,jbox)),
     2                    rmlexp(iaddr(2,jbox)),
     3                    rmlexp(iaddr(3,jbox)),
     4                    rmlexp(iaddr(4,jbox)),
     5                    rmlexp(iaddr(5,jbox)),
     6                    nterms(level1),scale(level0),center0,
     7                    mptemp1,mptemp2,mptemp3,mptemp4,mptemp5,
     8                    nterms(level0))
                          
                       else
                          call bh2dmpmp_carray(scale(level1),center1,
     1                    rmlexp(iaddr(1,jbox)),
     2                    rmlexp(iaddr(2,jbox)),
     3                    rmlexp(iaddr(3,jbox)),
     4                    rmlexp(iaddr(4,jbox)),
     5                    rmlexp(iaddr(5,jbox)),
     6                    nterms(level1),scale(level0),center0,
     7                    mptemp1,mptemp2,mptemp3,mptemp4,mptemp5,
     8                    nterms(level0),carray,ldc)

                       endif
                       call bh2dadd(mptemp1,rmlexp(iaddr(1,ibox)),
     1                 nterms(level0))
                       call bh2dadd(mptemp2,rmlexp(iaddr(2,ibox)),
     1                 nterms(level0))
                       call bh2dadd(mptemp3,rmlexp(iaddr(3,ibox)),
     1                 nterms(level0))
                       call bh2dadd(mptemp4,rmlexp(iaddr(4,ibox)),
     1                 nterms(level0))
                       call bh2dadd(mptemp5,rmlexp(iaddr(5,ibox)),
     1                 nterms(level0))
 2100               continue
                    if (ifprint .ge. 2) then
                       call prinf('=============*',x,0)
                    endif
                 endif
              endif
 2200      continue
!$OMP END PARALLEL DO
 2300   continue
c
c------------------------------------------------------------
c      DEBUGGING SEGMENT - once all multipole expansions are merged
c      to top level, one can compare it to the direct formation of the
c      expansion at the top level from the source locations.
c
ccc        call prinm(rmlexp(iaddr(1,1)),nterms(0))
c
ccc        call h2dformmp(ier,scale(0),source,charge,n,
ccc     1  	center,nterms(0),mptemp)
c
ccc        call prinm(mptemp,nterms(0))
c
ccc        call h2dmperr(rmlexp(iaddr(1,1)),mptemp,nterms(0),d)
ccc        call prin2('error in upward pass=*',d,1)
c
ccc        pause
ccc        stop
c      END DEBUGGING SEGMENT
c------------------------------------------------------------
c
c        t2=second()
c        t2=omp_get_wtime() 
c        timeinfo(3)=t2-t1
c
        if (ifprint .ge. 1) 
     $     call prinf('=== STEP 4 (mp to lo) ===*',i,0)
c        t1=second()
c        t1=omp_get_wtime() 
c       ... step 4, convert multipole expansions into the local ones
c
        call l2dterms_list2(epsfmm, itable, ier)
        do 4300 ilev=3,nlev+1
!$OMP PARALLEL DO DEFAULT(SHARED)
!$OMP$PRIVATE(ibox,box,center0,corners0,level0,list,nlists,nlist,itype)
!$OMP$PRIVATE(jbox,box1,center1,corners1,level1,ifdirect2,radius)
!$OMP$PRIVATE(mptemp1,lused,ier,i,j,vtemp,gatemp,gaatemp,cd,ilist)
!$OMP$PRIVATE(mptemp2,mptemp3,mptemp4,mptemp5)
!$OMP$PRIVATE(if_use_trunc,nterms_trunc,ii,jj) 
!$OMP$SCHEDULE(DYNAMIC)
ccc!$OMP$NUM_THREADS(1)
           do 4200 ibox=laddr(1,ilev),laddr(1,ilev)+laddr(2,ilev)-1
              call d2tgetb(ier,ibox,box,center0,corners0,wlists)
              if (ifprint .ge. 2) then
                 call prinf('ibox=*',ibox,1)
                 call prinf('box=*',box,15)
              endif
              level0=box(1)
              if (level0 .ge. 2) then
                 itype=2
                 call d2tgetl(ier,ibox,itype,list,nlist,wlists)
                 if (ifprint .ge. 2) then
                    call prinf('list2=*',list,nlist)
                 endif
                 
                 do 4150 ilist=1,nlist
                    jbox=list(ilist)
                    call d2tgetb(ier,jbox,box1,center1,corners1,wlists)
                    if( box1(10) .eq. 0 ) goto 4150
                    if (jbox.eq.0) goto 4150
                    if ((box(12).eq.0).and.(ifprune_list2.eq.1))
     $                 goto 4150
c
c       ... convert multipole expansions for all boxes in list 2 in local exp
c       ... if source is childless, evaluate directly (if cheaper)
c

                    level1=box1(1)
                    ifdirect2 = 0
                    if_use_trunc = 0
                    if (ifdirect2 .eq. 0) then
                       if(if_use_trunc .eq. 0) then
                          call bh2dzero(mptemp1,nterms(level1))
                          call bh2dzero(mptemp2,nterms(level1))
                          call bh2dzero(mptemp3,nterms(level1))
                          call bh2dzero(mptemp4,nterms(level1))
                          call bh2dzero(mptemp5,nterms(level1))
                          if(nterms(level0)+nterms(level1).gt.95) then
                             call bh2dmploc(scale(level1),center1,
     $                       rmlexp(iaddr(1,jbox)),
     1                       rmlexp(iaddr(2,jbox)),
     2                       rmlexp(iaddr(3,jbox)),
     3                       rmlexp(iaddr(4,jbox)),
     4                       rmlexp(iaddr(5,jbox)),
     5                       nterms(level1),scale(level0),center0,
     6                       mptemp1,mptemp2,mptemp3,mptemp4,mptemp5,
     7                       nterms(level0))

                          else
                             call bh2dmploc_carray(scale(level1),
     $                       center1,rmlexp(iaddr(1,jbox)),
     1                       rmlexp(iaddr(2,jbox)),
     2                       rmlexp(iaddr(3,jbox)),
     3                       rmlexp(iaddr(4,jbox)),
     4                       rmlexp(iaddr(5,jbox)),
     5                       nterms(level1),scale(level0),center0,
     6                       mptemp1,mptemp2,mptemp3,mptemp4,mptemp5,
     7                       nterms(level0),carray,ldc)

                          endif
                          call bh2dadd(mptemp1,rmlexp(iaddr(6,ibox)),
     1                    nterms(level0))
c
                          call bh2dadd(mptemp2,rmlexp(iaddr(7,ibox)),
     1                    nterms(level0))
c
                          call bh2dadd(mptemp3,rmlexp(iaddr(8,ibox)),
     1                    nterms(level0))
c
                          call bh2dadd(mptemp4,rmlexp(iaddr(9,ibox)),
     1                    nterms(level0))
c
                          call bh2dadd(mptemp5,rmlexp(iaddr(10,ibox)),
     1                    nterms(level0))
c
                       else
                          ii=box1(2)-box(2)
                          jj=box1(3)-box(3)
                          nterms_trunc=itable(ii,jj)
                          nterms_trunc=min(nterms(level0),nterms_trunc)
                          nterms_trunc=min(nterms(level1),nterms_trunc)

                          call bh2dzero(mptemp1,nterms_trunc)
                          call bh2dzero(mptemp2,nterms_trunc)
                          call bh2dzero(mptemp3,nterms_trunc)
                          call bh2dzero(mptemp4,nterms_trunc)
                          call bh2dzero(mptemp5,nterms_trunc)
                          if(nterms_trunc+nterms_trunc .gt. 95) then
                             call bh2dmploc(scale(level1),center1,
     $                       rmlexp(iaddr(1,jbox)),
     1                       rmlexp(iaddr(2,jbox)),
     2                       rmlexp(iaddr(3,jbox)),
     3                       rmlexp(iaddr(4,jbox)),
     4                       rmlexp(iaddr(5,jbox)),
     5                       nterms_trunc,scale(level0),center0,
     6                       mptemp1,mptemp2,mptemp3,mptemp4,mptemp5,
     7                       nterms_trunc)

                          else
                             call bh2dmploc_carray(scale(level1),
     $                       center1,rmlexp(iaddr(1,jbox)),
     1                       rmlexp(iaddr(2,jbox)),
     2                       rmlexp(iaddr(3,jbox)),
     3                       rmlexp(iaddr(4,jbox)),
     4                       rmlexp(iaddr(5,jbox)),
     5                       nterms_trunc,scale(level0),center0,
     6                       mptemp1,mptemp2,mptemp3,mptemp4,mptemp5,
     7                       nterms_trunc,carray,ldc)

                          endif
                          call bh2dadd(mptemp1,rmlexp(iaddr(6,ibox)),
     1                    nterms_trunc)

                          call bh2dadd(mptemp2,rmlexp(iaddr(7,ibox)),
     1                    nterms_trunc)

                          call bh2dadd(mptemp3,rmlexp(iaddr(8,ibox)),
     1                    nterms_trunc)

                          call bh2dadd(mptemp4,rmlexp(iaddr(9,ibox)),
     1                    nterms_trunc)

                          call bh2dadd(mptemp5,rmlexp(iaddr(10,ibox)),
     1                    nterms_trunc)

                       endif
                    endif
 4150            continue
              endif
 4200      continue
!$OMP END PARALLEL DO
 4300   continue
c
c        t2=second()
c        t2=omp_get_wtime() 
c        timeinfo(4)=t2-t1
c       
        if (ifprint .ge. 1) 
     $    call prinf('=== STEP 5 (split lo) ===*',i,0)
c        t1=second()
c        t1=omp_get_wtime() 
c       ... step 5, split all local expansions
        
        do 5300 ilev=3,nlev
!$OMP PARALLEL DO DEFAULT(SHARED)
!$OMP$PRIVATE(ibox,box,center0,corners0,level0,level,npts,nkids,radius)
!$OMP$PRIVATE(jbox,box1,center1,corners1,level1)
!$OMP$PRIVATE(mptemp1,lused,ier,i,j,vtemp,gatemp,cd)
!$OMP$PRIVATE(mptemp2,mptemp3,mptemp4,mptemp5,gaatemp) 
ccc!$OMP$SCHEDULE(DYNAMIC)
ccc!$OMP$NUM_THREADS(4) 
           do 5200 ibox=laddr(1,ilev),laddr(1,ilev)+laddr(2,ilev)-1
              call d2tgetb(ier,ibox,box,center0,corners0,wlists)
              call d2tnkids(box,nkids)
              if (nkids .ne. 0) then
                 level0=box(1)
                 if (level0 .ge. 2) then
                    if (ifprint .ge. 2) then
                       call prinf('ibox=*',ibox,1)
                       call prinf('box=*',box,15)
                       call prinf('nkids=*',nkids,1)
                       call prin2('center0=*',center0,2)
                    endif
c
                    do 5100 i = 1,4
                       jbox = box(4+i)
                       if (jbox.eq.0) goto 5100
                       call d2tgetb(ier,jbox,box1,center1,corners1,
     1                 wlists)
                       radius = (corners1(1,1) - center1(1))**2
                       radius = radius + (corners1(2,1)-center1(2))**2
                       radius = sqrt(radius)
                       if (ifprint .ge. 2) then
                          call prinf('jbox=*',jbox,1)
                          call prin2('radius=*',radius,1)
                          call prin2('center1=*',center1,2)
                       endif
                       level1=box1(1)
                       if(nterms(level0)+nterms(level1).gt.95) then
                          call bh2dlocloc(scale(level0),center0,
     1                    rmlexp(iaddr(6,ibox)),
     2                    rmlexp(iaddr(7,ibox)),
     3                    rmlexp(iaddr(8,ibox)),
     4                    rmlexp(iaddr(9,ibox)),
     5                    rmlexp(iaddr(10,ibox)),
     6                    nterms(level0),scale(level1),center1,
     7                    mptemp1,mptemp2,mptemp3,mptemp4,mptemp5,
     8                    nterms(level1))

                       else
                          call bh2dlocloc_carray(scale(level0),center0,
     1                    rmlexp(iaddr(6,ibox)),
     2                    rmlexp(iaddr(7,ibox)),
     3                    rmlexp(iaddr(8,ibox)),
     4                    rmlexp(iaddr(9,ibox)),
     5                    rmlexp(iaddr(10,ibox)),
     6                    nterms(level0),scale(level1),center1,
     7                    mptemp1,mptemp2,mptemp3,mptemp4,mptemp5,
     8                    nterms(level1),carray,ldc)

                       endif
                       call bh2dadd(mptemp1,rmlexp(iaddr(6,jbox)),
     1                 nterms(level1))

                       call bh2dadd(mptemp2,rmlexp(iaddr(7,jbox)),
     1                 nterms(level1))

                       call bh2dadd(mptemp3,rmlexp(iaddr(8,jbox)),
     1                 nterms(level1))

                       call bh2dadd(mptemp4,rmlexp(iaddr(9,jbox)),
     1                 nterms(level1))

                       call bh2dadd(mptemp5,rmlexp(iaddr(10,jbox)),
     1                 nterms(level1))
 5100               continue
                    if (ifprint .ge. 2) call prinf('=============*',x,0)
                 endif
              endif
              if (nkids .ne. 0) then
                 level=box(1)
                 if (level .ge. 2) then
                    if( ifprint .ge. 2 ) then
                       call prinf('ibox=*',ibox,1)
                       call prinf('box=*',box,15)
                       call prinf('nkids=*',nkids,1)
                    endif
                 endif
              endif
 5200      continue
!$OMP END PARALLEL DO
 5300   continue
c        t2=second()
c        t2=omp_get_wtime() 
c        timeinfo(5)=t2-t1
        return
        end
