      subroutine stokesSLPdirect(s1,s2,xs,ys,ns,xt,yt,nt,iequal,u1,u2)

      implicit real*8 (a-h,o-z)

      real *8 s1(ns),s2(ns)
      real *8 xs(ns),ys(ns)
      real *8 xt(nt),yt(nt)
      real *8 u1(nt),u2(nt)

      twopi = 8.d0*datan(1.d0)

      do i=1,nt
        u1(i) = 0.d0
        u2(i) = 0.d0
        do j=1,ns
          if ((i .ne. j .and. iequal .eq. 1) .or. iequal .eq. 0) then
            dis2 = (xt(i) - xs(j))**2.d0 + (yt(i) - ys(j))**2.d0
            rdots = ((xt(i) - xs(j))*s1(j) + 
     $          (yt(i) - ys(j))*s2(j))/dis2
            u1(i) = u1(i) - 5.d-1*dlog(dis2)*s1(j) + 
     $          rdots*(xt(i) - xs(j))
            u2(i) = u2(i) - 5.d-1*dlog(dis2)*s2(j) + 
     $          rdots*(yt(i) - ys(j))
          endif
        enddo
        u1(i) = u1(i)*5.d-1/twopi
        u2(i) = u2(i)*5.d-1/twopi
      enddo


      end
