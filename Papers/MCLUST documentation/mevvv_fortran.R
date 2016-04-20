
subroutine mevvv ( EQPRO, x, n, p, G, Vinv, z, maxi, tol, eps, 
     *                   mu, U, pro, w, S)

c This function is part of the MCLUST software described at
c       http://www.stat.washington.edu/mclust
c Copyright information and conditions for use of MCLUST are given at
c        http://www.stat.washington.edu/mclust/license.txt

      implicit NONE

      logical            EQPRO

      integer            n, p, G, maxi

      double precision   Vinv, eps, tol

c     double precision   x(n,p), z(n,G), w(p)
      double precision   x(n,*), z(n,*), w(*)

c     double precision   mu(p,G), U(p,p,G), pro(G), S(p,p)
      double precision   mu(p,*), U(p,p,*), pro(*), S(p,*)

      integer                 nz, p1, iter, i, j, k, l, j1

      double precision        piterm, hold, rcmin, rteps
      double precision        temp, term, cs, sn, umin, umax
      double precision        sumz, sum, detlog, const, hood, err
      double precision        prok, tmin, tmax, ViLog, zsum

      double precision        zero, one, two
      parameter              (zero = 0.d0, one = 1.d0, two = 2.d0)

      double precision        pi2log
      parameter              (pi2log = 1.837877066409345d0)

      double precision        FLMAX
      parameter              (FLMAX = 1.7976931348623157d308)

      double precision        SMALOG
      parameter              (SMALOG = -708.d0)

      external                ddot
      double precision        ddot

c------------------------------------------------------------------------------

      if (maxi .le. 0) return

      if (Vinv .gt. zero) then
        nz = G + 1
        ViLog = log(Vinv)
      else
        nz = G
        if (EQPRO) call dcopy( G, one/dble(G), 0, pro, 1)
      end if

      piterm = dble(p)*pi2log/two

      p1     = p + 1

      eps    = max(eps,zero)
      rteps  = sqrt(eps)

      tol    = max(tol,zero)

c     FLMAX  = d1mach(2)
      hold   = FLMAX/two
      hood   = FLMAX
      err    = FLMAX

c zero out the lower triangle of array U
      do k = 1, G

        do j = 1, p
          do l = 1, p
            S(l,j) = U(l,j,k)
          end do
        end do

        i = 1
        do j = 2, p
          call dcopy( p-i, zero, 0, S(j,i), 1)
          i = j
        end do

        do j = 1, p
          do l = 1, p
            U(l,j,k) = S(l,j)
          end do
        end do

      end do




      iter = 0
100   continue

      iter = iter + 1

      zsum = one
      do k = 1, G

        do j = 1, p
          do l = 1, p
            S(l,j) = U(l,j,k)
          end do
        end do

        do j = 1, p
          call dcopy( j, zero, 0, S(1,j), 1) # zero out entire S
        end do

        call dcopy( p, zero, 0, mu(1,k), 1) # zero out k-th column of mu

        sumz = zero
        do i = 1, n
          temp = z(i,k)
          sumz = sumz + temp # sumz = sum_i z_ik
          call daxpy( p, temp, x(i,1), n, mu(1,k), 1) # replace k-th column of mu by z_ik * y_i
        end do

        if (.not. EQPRO) pro(k) = sumz / dble(n) # obtain proportion estimates tau_k

        zsum = min(sumz,zsum) 
        if (sumz .gt. rteps) then
          call dscal( p, (one/sumz), mu(1,k), 1) # obtain mean estimates mu_k
          do i = 1, n
            call dcopy( p, x(i,1), n, w, 1)
            call daxpy( p, (-one), mu(1,k), 1, w, 1)
            call dscal( p, sqrt(z(i,k)), w, 1)
            j = 1
            do j1 = 2, p
              call drotg( S(j,j), w(j), cs, sn)
              call drot( p-j, S(j,j1), p, w(j1), 1, cs, sn)
              j = j1
            end do
            call drotg( S(p,p), w(p), cs, sn)
          end do
          do j = 1, p
            call dscal( j, one/sqrt(sumz), S(1,j), 1) # obtain covariance estimate?
          end do
        else
          call dcopy( p, FLMAX, 0,  z(1,k), 1) 
        end if

        do j = 1, p
          do l = 1, p
            U(l,j,k) = S(l,j)
          end do
        end do

      end do


      if (zsum .le. rteps) then
        tol  = zsum
        eps  = -FLMAX
        maxi = iter
        return
      end if


      if (Vinv .gt. zero) then

        term = zero
        do i = 1, n
          term = term + z(i,nz)
        end do
        temp    = term / dble(n)
        pro(nz) = temp

        call dcopy( n, ViLog, 0, z(1,nz), 1)

        if (EQPRO) then
          temp = (one - pro(nz))/dble(G)
          call dcopy( G, temp, 0, pro, 1)
        end if

      end if

      rcmin = FLMAX
      do k = 1, G
        do j = 1, p
          do l = 1, p
            S(l,j) = U(l,j,K) # S is the upper triangular cholesky decomposition of sigma_k
          end do
        end do
        call absrng( p, S, p1, umin, umax) # self-defined function?
        rcmin = min(umin/(one+umax),rcmin)
      end do

      if (rcmin .le. rteps) then
        tol  = rcmin
        eps  = FLMAX
        maxi = iter
        return
      end if


      ## E-step
      do k = 1, G

        do j = 1, p
          do l = 1, p
            S(l,j) = U(l,j,K)
          end do
        end do

c       temp = pro(k)

        detlog = zero
        do j = 1, p
          detlog = detlog + log(abs(S(j,j)))
        end do

        const = piterm+detlog

        do i = 1, n
          call dcopy( p, x(i,1), n, w, 1)
          call daxpy( p, (-one), mu(1,k), 1, w, 1)
          call dtrsv( 'U', 'T', 'N', p, S, p, w, 1)
          sum    = ddot( p, w, 1, w, 1)/two
c         z(i,k) = temp*exp(-(const+sum))
          z(i,k) = -(const+sum)
        end do

      end do









      hood = zero
      do i = 1, n
        tmin =  FLMAX
        tmax = -FLMAX

        do k = 1, nz
          prok = pro(k)
          if (prok .eq. zero) then
            z(i,k) = zero
          else
            temp   = log(prok) + z(i,k)
            tmin   = min(tmin,temp)
            tmax   = max(tmax,temp)
            z(i,k) = temp
          end if
        end do

        sum   = zero
        do k = 1, nz
          if (pro(k) .ne. zero) then
            temp =  z(i,k) - tmax
            if (temp .ge. SMALOG) then
              z(i,k) = exp(temp)
              sum    = sum + z(i,k)
            else
              z(i,k) = zero
            end if
          end if
        end do

        hood  = hood + (log(sum)+tmax)
        call dscal( nz, (one/sum), z(i,1), n) # obtain estimates for z

      end do

      err  = abs(hold-hood)/(one+abs(hood))
      hold = hood

      if (err .gt. tol .and. iter .lt. maxi) goto 100

c     w(1) = rcmin

      tol  = err
      eps  = hood
      maxi = iter

      return
      end

