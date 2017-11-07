
      subroutine init_grid ( nx, ny, nz, xx, yy, zz )
!
! Initializes parameters 
!
      implicit none
      integer*4 ndim
      parameter(ndim=3)
      integer*4 nd(ndim), m, nx, ny, nz, n, i, j, nseed
      real*8 xx(nx), yy(ny), zz(nz), x, vol, d, big
      real*8 dmin(3), dmax(3), volume, energy, 
     & dist, dist_max, dist_min, move_max, move_min
      common/geo/ dmin, dmax, 
     & dist, dist_max, dist_min, move_max, move_min
      real*8 sigma, emax, emin
      common/potential/sigma, emin, emax
      parameter(nseed=86456, big=1.e20)
      call srand(nseed)
!     Find the bounding box of the domain

      nd(1) = nx
      nd(2) = ny
      nd(3) = nz
      dmin(1) = xx(2)
      dmax(1) = xx(nx-1)
      dmin(2) = yy(2)
      dmax(2) = yy(ny-1)
      dmin(3) = zz(2)
      dmax(3) = zz(nz-1)
!     Estimate volume per node
      write(*,*)'Bounding box:',nx-2,ny-2,nz-2,n
      volume = 1.0
      dist = 0.0
      j = 0
      m = 1
      do i=1,ndim
         d = dmax(i) - dmin(i)
         write(*,*)i,dmin(i),dmax(i),d
         volume = volume * d
         if (dist < d) then
           dist = d
           j = i
         endif
         m = m * (nd(i)-2)
      enddo

      write(*,*)'Volume: ',volume
!     Compute distance based on volume:
!      dist = 4*(volume/(real(m)))**(1./3.)

!     Compute distance based on grid cell size:
      dist = 3*dist/real(nd(j)-2)
      write(*,*)'Average distance between nodes: ',dist

      dist_max = 0.3 * dist ! max displacement after which
      ! triangulation is repeated - used only in Persson's algorithm
      dist_min = 1.e-3 * dist

      sigma = 1.5 * dist ! interaction radius of potential function
      emax = energy( dist )
      emin = energy( sigma )
      move_max = 0.3 * sigma
      move_min = 0.01 * move_max
      write(*,*)'Inter-node distance: ',dist
      write(*,*)'Min Displacement: ',dist_min
      write(*,*)'Max Displacement: ',dist_max
      return
      end

!
! Read grid dimensions from file
!
      subroutine read_grid_dimensions (file_name, nx, ny, nz)
      implicit none

      integer*4 nx,ny,nz
      logical got_one
      character ( len = * ) file_name
      integer ( kind = 4 ) ierror
      integer ( kind = 4 ) num
      integer ( kind = 4 ) input_unit
      integer ( kind = 4 ) ios
      character ( len = 255 ) line
      ierror = 0
      
      call get_unit ( input_unit )
      
      open(unit=input_unit,file=file_name,status='old',
     & form='formatted', access='sequential', iostat=ios)
      if ( ios /= 0 ) then
        num = -1
        write (*,'(a)') ' '
        write (*,'(a)') 'FILE_COLUMN_COUNT - Fatal error!'
        write (*,'(a)') '  Could not open the file:'
        write (*,'(a)') trim ( file_name )
        return
      end if

      write(*,*)'READ GRID DIMENSIONS:'
      read ( input_unit, *)nx,ny,nz
      write(*,*)'Dimensions: ',nx,ny,nz

      close ( unit = input_unit )
      return
      end

!
! Read grid coordinates
!
      subroutine load_grid_coordinates (file_name, x, y, z)

      implicit none

      integer*4 i,j,k,nx,ny,nz
      real*8 x(*),y(*),z(*)
      character ( len = * ) file_name
      integer ( kind = 4 ) ierror
      integer ( kind = 4 ) num
      integer ( kind = 4 ) input_unit
      integer ( kind = 4 ) output
      integer ( kind = 4 ) ios
      character ( len = 255 ) line
      ierror = 0
      
      call get_unit ( input_unit )
      
      open(unit=input_unit,file=file_name,status='old',
     & form='formatted', access='sequential', iostat=ios)
      if ( ios /= 0 ) then
        num = -1
        write (*,'(a)') ' '
        write (*,'(a)') 'FILE_COLUMN_COUNT - Fatal error!'
        write (*,'(a)') '  Could not open the file:'
        write (*,'(a)') trim ( file_name )
        return
      end if
      read ( input_unit, *)nx,ny,nz
      write(*,*)'Grid Dimensions: ',nx,ny,nz
      read ( input_unit, *) 
     & (x(i),i=1,nx),
     & (y(j),j=1,ny),
     & (z(k),k=1,nz)
      close ( unit = input_unit )
      return
      end

!
! Read a scalar variable defined on a rectilinear grid
!
      subroutine load_grid_scalar (file_name, 
     & nx, ny, nz, phi)

      implicit none

      integer*4 nx,ny,nz,ni,nj,nk,i,j,k
      real*8 phi(nx,ny,nz)
      character ( len = * ) file_name
      integer ( kind = 4 ) ierror
      integer ( kind = 4 ) input_unit
      integer ( kind = 4 ) output
      integer ( kind = 4 ) ios
      character ( len = 255 ) line
      ierror = 0

      write(*,*)'Load grid scalar:',nx,ny,nz
      
      call get_unit ( input_unit )
      
      open(unit=input_unit,file=file_name,status='old',
     & form='formatted', access='sequential', iostat=ios)
      if ( ios /= 0 ) then
        write (*,'(a)') ' '
        write (*,'(a)') 'FILE_COLUMN_COUNT - Fatal error!'
        write (*,'(a)') '  Could not open the file:'
        write (*,'(a)') trim ( file_name )
        return
      end if
      write(*,*)'READ GRID:'
      read ( input_unit, *)ni,nj,nk
      write(*,*)'Dimensions: ',ni,nj,nk
      if (ni .ne. nx .or. nj .ne. ny .or.
     &    nk .ne. nz) then
         write(*,*)'Dimensions mismatch'
         return
      endif
      write(*,*)'Reading ',file_name
      read ( input_unit, *)
     & (((phi(i,j,k),i=1,nx),j=1,ny),k=1,nz)

      close ( unit = input_unit )
      write (*,*) nx*ny*nz,'scalars read'

      return
      end
!
! Compute gradient of a scalar variable
!
      subroutine gradient ( nx, ny, nz, xx, yy, zz, phi, grad )
      parameter (ndim=3)
      integer*4 nx,ny,nz,ix,iy,iz
      real*8 xx(nx), yy(ny), zz(nz), 
     &       phi(nx,ny,nz), grad(ndim, nx, ny, nz),
     &       dx, dy, dz, dpx, dpy, dpz
      do iz=1,nz
      do iy=1,ny
      do ix=1,nx
      do i=1,3
      grad(i,ix,iy,iz) = 0.0
      enddo
      enddo
      enddo
      enddo
! Compute at the internal nodes
      do iz=2,nz-1
      dz = zz(iz+1)-zz(iz-1)
      do iy=2,ny-1
      dy = yy(iy+1)-yy(iy-1)
      do ix=2,nx-1
      dx = xx(ix+1)-xx(ix-1)
      dpx = phi(ix+1,iy  ,iz  )-phi(ix-1,iy  ,iz  )
      dpy = phi(ix  ,iy+1,iz  )-phi(ix  ,iy-1,iz  )
      dpz = phi(ix  ,iy  ,iz+1)-phi(ix  ,iy  ,iz-1)
      grad(1,ix,iy,iz) = dpx/dx
      grad(2,ix,iy,iz) = dpy/dy
      grad(3,ix,iy,iz) = dpz/dz
c      write(*,*)grad(1,ix,iy,iz),grad(2,ix,iy,iz),grad(3,ix,iy,iz)
      enddo
      enddo
      enddo

! Compute at the boundary excluding the corners:
      iz=1
      dz = zz(iz+1)-zz(iz)
      do iy=2,ny-1
      dy = yy(iy+1)-yy(iy-1)
      do ix=2,nx-1
      dx = xx(ix+1)-xx(ix-1)
      dpx = phi(ix+1,iy  ,iz  )-phi(ix-1,iy  ,iz  )
      dpy = phi(ix  ,iy+1,iz  )-phi(ix  ,iy-1,iz  )
      dpz = phi(ix  ,iy  ,iz+1)-phi(ix  ,iy  ,iz  )
      grad(1,ix,iy,iz) = dpx/dx  
      grad(2,ix,iy,iz) = dpy/dy   
      grad(3,ix,iy,iz) = dpz/dz   
      enddo
      enddo

      iz=nz
      dz = zz(iz)-zz(iz-1)
      do iy=2,ny-1
      dy = yy(iy+1)-yy(iy-1)
      do ix=2,nx-1
      dx = xx(ix+1)-xx(ix-1)
      dpx = phi(ix+1,iy  ,iz  )-phi(ix-1,iy  ,iz  )
      dpy = phi(ix  ,iy+1,iz  )-phi(ix  ,iy-1,iz  )
      dpz = phi(ix  ,iy  ,iz  )-phi(ix  ,iy  ,iz-1)
      grad(1,ix,iy,iz) = dpx/dx  
      grad(2,ix,iy,iz) = dpy/dy   
      grad(3,ix,iy,iz) = dpz/dz   
      enddo
      enddo

      do iz=2,nz-1
      dz = zz(iz+1)-zz(iz-1)
      iy=1
      dy = yy(iy+1)-yy(iy  )
      do ix=2,nx-1
      dx = xx(ix+1)-xx(ix-1)
      dpx = phi(ix+1,iy  ,iz  )-phi(ix-1,iy  ,iz  )
      dpy = phi(ix  ,iy+1,iz  )-phi(ix  ,iy  ,iz  )
      dpz = phi(ix  ,iy  ,iz+1)-phi(ix  ,iy  ,iz-1)
      grad(1,ix,iy,iz) = dpx/dx  
      grad(2,ix,iy,iz) = dpy/dy   
      grad(3,ix,iy,iz) = dpz/dz   
      enddo
      enddo

      do iz=2,nz-1
      dz = zz(iz+1)-zz(iz-1)
      iy=ny
      dy = yy(iy  )-yy(iy-1)
      do ix=2,nx-1
      dx = xx(ix+1)-xx(ix-1)
      dpx = phi(ix+1,iy  ,iz  )-phi(ix-1,iy  ,iz  )
      dpy = phi(ix  ,iy  ,iz  )-phi(ix  ,iy-1,iz  )
      dpz = phi(ix  ,iy  ,iz+1)-phi(ix  ,iy  ,iz-1)
      grad(1,ix,iy,iz) = dpx/dx  
      grad(2,ix,iy,iz) = dpy/dy   
      grad(3,ix,iy,iz) = dpz/dz   
      enddo
      enddo

      do iz=2,nz-1
      dz = zz(iz+1)-zz(iz-1)
      do iy=2,ny-1
      dy = yy(iy+1)-yy(iy-1)
      ix=1
      dx = xx(ix+1)-xx(ix  )
      dpx = phi(ix+1,iy  ,iz  )-phi(ix  ,iy  ,iz  )
      dpy = phi(ix  ,iy+1,iz  )-phi(ix  ,iy-1,iz  )
      dpz = phi(ix  ,iy  ,iz+1)-phi(ix  ,iy  ,iz-1)
      grad(1,ix,iy,iz) = dpx/dx  
      grad(2,ix,iy,iz) = dpy/dy   
      grad(3,ix,iy,iz) = dpz/dz   
      enddo
      enddo

      do iz=2,nz-1
      dz = zz(iz+1)-zz(iz-1)
      do iy=2,ny-1
      dy = yy(iy+1)-yy(iy-1)
      ix=nx
      dx = xx(ix  )-xx(ix-1)
      dpx = phi(ix+1,iy  ,iz  )-phi(ix  ,iy  ,iz  )
      dpy = phi(ix  ,iy+1,iz  )-phi(ix  ,iy-1,iz  )
      dpz = phi(ix  ,iy  ,iz+1)-phi(ix  ,iy  ,iz-1)
      grad(1,ix,iy,iz) = dpx/dx  
      grad(2,ix,iy,iz) = dpy/dy   
      grad(3,ix,iy,iz) = dpz/dz   
      enddo                       
      enddo                       
      return
      end
!
! Tri-Linear Interpolation
!https://en.wikipedia.org/wiki/Trilinear_interpolation
!
      subroutine interp3 (x, xmin, xmax, v, u)
! Interpolate variable v given at eight corners of 
! the cube defined by xmin,xmax to the point x
      real*8 
     & x(3), ! point of interpolation
     & xmin(3),xmax(3), ! two limitig corners
     & v(2,2,2), ! variable values at 8 corners
     & u, ! the result of interpolation
     & xd1,yd1,zd1,xd2,yd2,zd2,
     & c11,c12,c21,c22,c1,c2
      xd2 = (x(1)-xmin(1))/(xmax(1)-xmin(1))
      yd2 = (x(2)-xmin(2))/(xmax(2)-xmin(2))
      zd2 = (x(3)-xmin(3))/(xmax(3)-xmin(3))
      xd1 = 1.0 - xd2
      yd1 = 1.0 - yd2
      zd1 = 1.0 - zd2
      c11 = v(1,1,1)*xd1 + v(2,1,1)*xd2
      c12 = v(1,1,2)*xd1 + v(2,1,2)*xd2
      c21 = v(1,2,1)*xd1 + v(2,2,1)*xd2
      c22 = v(1,2,2)*xd1 + v(2,2,2)*xd2
      c1 = c11*yd1 + c21*yd2
      c2 = c12*yd1 + c22*yd2
      u = c1*zd1 + c2*zd2
      end
!
! Find index of a variabe in a sorged array
! (linear search)
!
      function ind0(x, n, arr)
      integer*4 ind, n, i
      real*8 x, arr(n)
      if (arr(1) > x) then
         ind = 0
         return
      endif
      if (arr(n) < x) then
         ind = n+1
         return
      endif
      do i=2,n
         if (arr(i) >= x) then
            exit
         endif
      enddo
      ind = i
      return
      end
!
! Find point's index in a sorted array
! (binary search)
!
      function ind(x, n, arr)
      parameter(max_iter = 999)
      integer*4 ind, n, i0, i1, i, j
      real*8 x, arr(n)
      if (arr(1) > x) then
         ind = 0
         return
      endif
      if (arr(n) < x) then
         ind = n+1
         return
      endif
      i0 = 1
      i1 = n-1
      do j=1,max_iter
         i=(i0+i1)/2
         a0 = arr(i)
         a1 = arr(i+1)
         if (a0 < x .and. x <= a1) then
            ind = i+1
            return
         endif
         if (x <= a0) then
            i1 = i
            cycle
         endif
         if (x > a1) then
            i0 = i+1
            cycle
         endif
      enddo
      write(*,*)'WARNING: index search exceeded max interations of ',j
      ind = 0
      return
      end
!
! Creating nodes inside the iso-surface of a variable
!
      subroutine make_nodes (nx, ny, nz, xx, yy, zz,
     & q, nn, nodes, nodeq)
! Nodes are inserted in a uniform rectangular grid
! with grid cell size d
      integer*4 nx, ny, nz, 
     & nn, ! number of nodes given
     & n ! number of nodes calculated
     & mx, my, mz, ! grid dimensions
     & ix, iy, iz, jx, jy, jz, i, j, k
c     & tet(3,4) ! initial tetrahedron
      real*8 xx(nx), yy(ny), zz(nz), q(nx,ny,nz), 
     & nodes(3,nn), nodeq(nn), pmin(3), pmax(3),
     & xmin, xmax, ymin, ymax, zmin, zmax, 
     & d, dx, dy, dz, x, y, z, v(2,2,2), p(3), u
      real*8 dmin(3), dmax(3), 
     & dist, dist_max, dist_min, move_max, move_min
      common/geo/ dmin, dmax, 
     & dist, dist_max, dist_min, move_max, move_min

! dist: node-insertion interval
!
      n = 0 ! number of nodes 
      xmin = 0.5*(xx(1)+xx(2))
      ymin = 0.5*(yy(1)+yy(2))
      zmin = 0.5*(zz(1)+zz(2))
      xmax = xx(nx-1)
      ymax = yy(ny-1)
      zmax = zz(nz-1)
      d = 0.5 * dist
      write(*,*)'X:',xmin,xmax
      write(*,*)'Y:',ymin,ymax
      write(*,*)'Z:',zmin,zmax
      write(*,*)'dist=',dist

      dx = xmax-xmin
      dy = ymax-ymin
      dz = zmax-zmin
      mx = int(dx/dist)
      my = int(dy/dist)
      mz = int(dz/dist)
      write(*,*)'Max number of points:',mx,my,mz !-
!     new z:
      z0 = zmin
      do jz=1,mz ! LOOP Z
      z = z0 + (jz-1)*dist
      iz1 = ind(z, nz, zz)
      if (iz1 > nz) then
         exit
      else 
         if (iz1 == 0) then
            cycle
         endif
      endif
      iz = iz1-1
!     new y:
      y0 = ymin
      if (mod(jz,2)==0) then
      ! shift every other layer
         y0 = ymin + 0.5*dist
      endif
      do jy=1,my !LOOP Y
      y = y0 + (jy-1)*dist
      iy1 = ind(y, ny, yy)
      if (iy1 > ny) then
         exit
      else 
         if (iy1 == 0) then
            cycle
         endif
      endif
      iy = iy1-1
!     new x:
      x0 = xmin
      if (mod(jy,2)==0) then
         x0 = xmin + 0.5*dist
      endif
      do jx=1,mx ! LOOP X
      x = x0 + (jx-1)*dist 
      ix1 = ind(x, nx, xx)
      if (ix1 > nx) then
         exit
      else 
         if (ix1 == 0) then
            cycle
         endif
      endif
      ix = ix1-1
      p(1) = x
      p(2) = y
      p(3) = z
      pmin(1) = xx(ix)
      pmax(1) = xx(ix1)
      pmin(2) = yy(iy)
      pmax(2) = yy(iy1)
      pmin(3) = zz(iz)
      pmax(3) = zz(iz1)
      do i=0,1
      do j=0,1
      do k=0,1
      v(i+1,j+1,k+1) = q(ix+i,iy+j,iz+k)
      enddo
      enddo
      enddo
      call interp3 (p, pmin, pmax, v, u)
      if (u > 0.0) then
         n = n + 1
         if (nn > 0) then
            do i=1,3
               nodes(i,n) = p(i)
            enddo
            nodeq(n) = u
         endif
      endif
      enddo ! next x
      enddo ! next y
      enddo ! next z
      if (nn == 0) then
         nn = n
      else 
         if (nn .ne. n) then
            write(*,*)'Given node number',nn,
     &                ' is not equal to computed',n
            nn = 0
            return
         endif
      endif
      return
      end
!==================================================
!
! Node moving routines below use boxing technique 
! to achieve a quasi-linear scale-up in the number 
! of nodes.
!
! The boxing technique is as follows:
!
! - Domain is split into equally sized boxes 
! - Each point is assigned to a box where it is currently located
! - Interaction between the points is restricted to within a single
!   box and its neighboring boxes.
!--------------------------------------------------
!
! Assign a point to a box
!
      subroutine assign_box(ip, np, 
     & nx, ny, nz, mpb,
     & ix, iy, iz, ipx, ipy, ipz,
     & ipb, ipp, ierr)
      implicit none
      integer*4 ip, np, nx, ny, nz, mpb, mb, 
     & ix, iy, iz, ierr,
     & ipx(np), ipy(np), ipz(np), ipb(np),
     & ipp(3,mpb+2,nx,ny,nz),
     & i, j, k, n, pind, prev, next, head, free,
     & head_next, new_slot, next_to_new, free_next,
     & inext, iprev, jx, jy, jz, ib
      parameter(pind=1, prev=2, next=3, head=1, free=2)
!     ip = current particle index
!     np = total number of particles
!     nx,ny,nz = number of boxes in three directions
!     ix,iy,iz = current particle location indices
!     mpb = max number of particles in each box
!     ipp(1,i,x,y,z) = index of particle in a box 'i'
!     at x,y,z
!     ipb(i) = index of a box
!     ipx(i),ipy(i),ipz(i) = indexes of a box on a grid
!                            for the i-th particle
      if (ix == ipx(ip) .and. iy == ipy(ip) .and. iz == ipz(ip)) return

      ierr = 0

!     Link to the new box:

      n = ipp(pind,free,ix,iy,iz)
      if (n == 0) then
        write(*,*)'ERROR: box full at ',ix,iy,iz
        ierr = 1
        return
      endif
! Use next of free as a new slot
      new_slot = ipp(next,free,ix,iy,iz)
! Unlink new slot:
      next_to_new = ipp(next,new_slot,ix,iy,iz)
      ipp(prev,next_to_new,ix,iy,iz) = free
      ipp(next,free,ix,iy,iz) = next_to_new
! Decrement free slots count:
      ipp(pind,free,ix,iy,iz) = n - 1
! Store new point to new_slot:
      ipp(pind,new_slot,ix,iy,iz) = ip
! Link new point to head list:
      head_next = ipp(next,head,ix,iy,iz)
      ipp(prev,head_next,ix,iy,iz) = new_slot
      ipp(next,head,ix,iy,iz) = new_slot
      ipp(next,new_slot,ix,iy,iz) = head_next
      ipp(prev,new_slot,ix,iy,iz) = head
! Increment new slots count
      n = ipp(pind,head,ix,iy,iz)
      ipp(pind,head,ix,iy,iz) = n + 1
      
!     Unlink from the old box:

      ib = ipb(ip)

      if (ib > 0) then ! old particle
!     i.e. particle was assigned to a box                  

      jx = ipx(ip)
      jy = ipy(ip)
      jz = ipz(ip)
      
      ipp(pind,ib,ix,iy,iz)=0
! Unlink
      inext = ipp(next,ib,ix,iy,iz)
      iprev = ipp(prev,ib,ix,iy,iz)
      ipp(prev,inext,ix,iy,iz) = iprev
      ipp(next,iprev,ix,iy,iz) = inext
      n = ipp(pind,head,ix,iy,iz)
      ipp(pind,head,ix,iy,iz) = n - 1
! Reassign ib pointer to the free list:
! Put it ahead of free
      free_next = ipp(next,free,jx,jy,jz)
      ipp(prev,free_next,jx,jy,jz) = ib
      ipp(next,free,jx,jy,jz) = ib
      ipp(pind,ib,jx,jy,jz) = 0
      n = ipp(pind,free,jx,jy,jz) 
      ipp(pind,free,jx,jy,jz) = n + 1

      endif

      ipb(ip) = new_slot
      ipx(ip) = ix
      ipy(ip) = iy
      ipz(ip) = iz

! Debugging output if needed:
!      write(*,*)'ASSIGNED P=',ip
!      call raw_box_info (
!     &     ix,iy,iz,
!     &     nx,ny,nz,
!     &     mpb,ipp
!     & )
!      call box_info (
!     &     ix,iy,iz,
!     &     nx,ny,nz,
!     &     mpb,ipp
!     & )

      end
!
!-------------------------------------------
!
! MOVING NODES
!
! Node positions are adjusted to a possibly
! uniform distribution.
! Here a Monte-Carlo technique is used.
! A faster implementation based on a forcing
! function is possible, but is less stable.
! 
! This routine uses boxing for quasi-linear speedup
! 
      subroutine move_nodes ( 
     & nx, ny, nz, xx, yy, zz,
     & field, nn, nodes, nodeq,
     & nbx, nby, nbz, mpb, 
     & ipb, ipx, ipy, ipz, 
     & ipp, toteng, nmoves, disp)

      implicit none

      integer*4 ndim ! number of nodes given
      parameter(ndim=3)

      integer*4 
     & nn, mn, n, m,   ! number of nodes and loop counters
     & nx, ny, nz, ! grid dimensions
     & ix, iy, iz, jx, jy, jz, i, j, k, l,
     & itry, ntry, ierr, nmoves
      parameter(ntry=99)
      real*8 xx(nx), yy(ny), zz(nz), field(nx,ny,nz), 
     & nodes(ndim,nn), nodeq(nn), xmin(ndim), xmax(ndim), toteng,
     & f, x, e0, e1, pos1(ndim), pos0(ndim), 
     & dmin(ndim), dmax(ndim), qold, qnew,
     & dist, dist_max, dist_min, move_max, move_min,
     & d, disp

! Domain dimensions and scales:
      common/geo/ dmin, dmax, 
     & dist, dist_max, dist_min, move_max, move_min
      integer*4 ipb(nn), ipx(nn), ipy(nn), ipz(nn),
     & ipp(ndim, mpb+2, nbx, nby, nbz)
      integer*4 nbx, nby, nbz, mpb 

      toteng = 0.0 ! total energy
      mn = 0
      nmoves = 0
      disp = 0 ! max displacement
!
! Looping over all nodes
!
      do n=1,nn

        qold = nodeq(n) !-
        qnew = qold !-

        do i=1,ndim ! remember old position
          pos0(i)=nodes(i,n)
        enddo
        do i=1,ndim
          x = pos0(i)
          if (x <= dmin(i) .or. x >= dmax(i)) exit
        enddo
        if (i < ndim) cycle ! frozen

        do itry=1,ntry ! Make trial move

          do i=1,ndim
            x = pos0(i) + move_max*(2*rand()-1.0)
            if (x > dmin(i) .and. x < dmax(i)) then
              pos1(i) = x
            else 
              pos1(i) = pos0(i)
            endif
          enddo 

          call getvar ( pos1, nx, ny, nz, xx, yy, zz,
     &    field, qnew, ierr )

          if (ierr > 0) then
            cycle
          endif

          if ( qnew > 0.0 ) exit

        enddo ! END TRIAL MOVE

        if (itry > ntry) cycle

        call energies ( 
     &       pos0, pos1, n, nn, nodes,
     &       nbx, nby, nbz, mpb, 
     &       ipp, e0, e1
     & )

        if (e1 < e0) then
! If potential decreased accept the new position
          d = 0.0
          do i=1,ndim
            x = pos1(i)
            d = d + (x - pos0(i))**2
            nodes(i,n) = x
            ipx(n)=ix
            ipy(n)=iy
            ipz(n)=iz
          enddo
          nodeq(n) = qnew 
          toteng = toteng + e1
          nmoves = nmoves + 1
          if (disp < d) disp = d
        else ! keep the old energy and position
          toteng = toteng + e0
        endif
        mn = mn + 1
      enddo ! next n
      toteng = toteng / mn
      disp = sqrt(disp)
      return
      end
!---------------------------------------------
!
! Normal projection of surface nodes
!
! Moving all near surface nodes to the surface
! along the gradient lines of a field.
! 
      subroutine project_surface_nodes (
     & nx, ny, nz, xx, yy, zz,
     & field, grad, nn, nodes, nodeq,
     & nbx, nby, nbz, mpb, 
     & ipb, ipx, ipy, ipz, 
     & ipp, toteng)

      implicit none

      integer*4 ndim ! number of nodes given
      parameter(ndim=3)

      integer*4 
     & nn, n,    ! number of nodes and loop counters
     & nx, ny, nz, ! grid dimensions
     & ix, iy, iz, jx, jy, jz, i, j, k, l,
     & ierr
      real*8 xx(nx), yy(ny), zz(nz), field(nx,ny,nz), 
     & nodes(ndim,nn), nodeq(nn), xmin(ndim), xmax(ndim), toteng,
     & r, f, x, e, e0, e1, pos(ndim), 
     & dmin(ndim), dmax(ndim), 
     & dist, dist_max, dist_min, move_max, move_min

! Domain dimensions and scales:
      common/geo/ dmin, dmax, 
     & dist, dist_max, dist_min, move_max, move_min
      integer*4 ipb(nn), ipx(nn), ipy(nn), ipz(nn),
     & ipp(ndim, mpb+2, nbx, nby, nbz)
      integer*4 nbx, nby, nbz, mpb, 
     &          ibx, iby, ibz, ibn

! Surface projection variables:
      real*8 g, q, u, d(3), qold, qnew, 
     & vert(2,2,2), grad(ndim,nx,ny,nz)
      real*8 surface_proximity
      common/surface/surface_proximity
!
! Looping over all nodes
!
      do n=1,nn

        qold = nodeq(n) 

        do i=1,ndim ! remember old position
          pos(i)=nodes(i,n)
        enddo
        do i=1,ndim
          x = pos(i)
          if (x <= dmin(i) .or. x >= dmax(i)) exit
        enddo
        if (i < ndim) cycle ! outside domain -> freeze

        call getixyz (pos, 
     &  nx, ny, nz, 
     &  xx, yy, zz, 
     &  ix, iy, iz, 
     &  xmin, xmax, 
     &  ierr)

        if (ierr > 0) then
          cycle
        endif

        g = 0.0
        do l=1,3

        do k=0,1
        do j=0,1
        do i=0,1
        vert(i+1,j+1,k+1) = grad(l,ix+i,iy+j,iz+k)
        enddo
        enddo
        enddo

        call interp3 ( pos, xmin, xmax, vert, u )

        d(l) = u
        g = g + u*u

        enddo !l

        ! Project onto the zero surface:
        do i=1,ndim
          pos(i) = pos(i) - qold * d(i) / g
        enddo

        call getvar ( pos, nx, ny, nz, xx, yy, zz,
     &  field, qnew, ierr )

        if (ierr > 0) then
          cycle
        endif

        if ( abs(qnew) > abs(qold) ) cycle 

        do i=1,ndim
          nodes(i,n) = pos(i)
          ipx(n)=ix
          ipy(n)=iy
          ipz(n)=iz
        enddo
        nodeq(n) = qnew 
      enddo ! next n
      return
      end
!
! Potential energy to make MC move
!
      function energy ( r )
      implicit none
      real*8 energy, r, eta, eps, x
      parameter(eta=1.0, eps=1.0e-20)
      real*8 dmin(3), dmax(3), 
     & dist, dist_max, dist_min, move_max, move_min
      common/geo/ dmin, dmax, 
     & dist, dist_max, dist_min, move_max, move_min
      real*8 sigma, emin, emax
      common/potential/sigma, emin, emax
      x = sigma/max(r, eps)
      energy = x**2
!?      energy = x**12 - 2*x**6 !LJ potential
      return
      end
!
! Compute energies for two positions
!
      subroutine energies ( 
     &  pos0, pos1,
     &  n, np, p,
     &  nbx, nby, nbz, mpb, 
     &  ipp, 
     &  e0, e1
     & )
      implicit none
      integer*4 ndim
      parameter(ndim=3)
      integer*4 m, n, np, nn, nbx, nby, nbz, mpb, 
     &          i, j, ix, iy, iz, jx, jy, jz,
     &          pind, prev, next, head,
     & ipp(ndim, mpb+2, nbx, nby, nbz)
! Boxing variables for quasi-linear scale-up:
      real*8 p(ndim, np), pos0(ndim), pos1(ndim),
     & db, bmin(ndim), bmax(ndim), 
     & x, d0, d1, e, e0, e1, r0, r1, energy
      common/pb/ db, bmin, bmax
      parameter(pind=1,prev=2,next=3,head=1)

        e0 = 0.0 ! initial potential
        e1 = 0.0 ! potential after a trial move
        nn = 0   ! number of neighbors 

        ix = max(min(int((pos1(1)-bmin(1))/db),nbx-1),2)
        iy = max(min(int((pos1(2)-bmin(2))/db),nby-1),2)
        iz = max(min(int((pos1(3)-bmin(3))/db),nbz-1),2)

        do jz=iz-1,iz+1 ! loop over the current
        do jy=iy-1,iy+1 ! and neighboring boxes
        do jx=ix-1,ix+1

        if (ipp(pind,head,jx,jy,jz) > 0) then
        ! the above checks the nuber of nodes in the box         
        ! ipp(pind,head,*) = number of nodes in the box

        i = head ! start with the head of the node list
        do ! loop over all points in the box
          i = ipp(next,i,jx,jy,jz) ! next point index
          m = ipp(pind,i,jx,jy,jz) ! point index in 
                                   ! p(3,np) array
          if (i == head) exit 
          if (m == n) cycle

          r0=0.0
          r1=0.0
          do j=1,ndim
            x = p(j, m)
            d0 = x - pos0(j)
            d1 = x - pos1(j)
            r0 = r0 + d0*d0
            r1 = r1 + d1*d1
          enddo
          e0 = e0 + energy(sqrt(r0))
          e1 = e1 + energy(sqrt(r1))
          nn = nn + 1

        enddo
        endif

        enddo ! next jx
        enddo ! next jy
        enddo ! next jz
        e0 = e0/real(nn)
        e1 = e1/real(nn)
      end

!------------------------------------------------------------
!
! SURFACE NODES PROJECTION ROUTINES
!
      subroutine count_surface_nodes (
     &        np, pv, ! coordinates of all points in a volume
     &        nx, ny, nz, ! grid dimensions
     &        xx, yy, zz, ! grid coordinates
     &        fv,         ! scalar field in a volume
     &        ns  ! returned number of surface nodes
     & )

      implicit none
      integer*4 ndim
      parameter(ndim=3)

      integer*4 np, ns, nx, ny, nz,
     & ip,i,j,k,n, ierr
      real*8 pv(ndim, np), 
     & xx(nx), yy(ny), zz(nz), fv(nx,ny,nz),
     & pos(ndim), f

      real*8 surface_proximity
      common/surface/surface_proximity
      ! Surface proximity is measure in terms of 
      ! the scalar field. 

      surface_proximity = 0.01 
      ! layer thickness around the surface to define 
      ! the nodes at the surface

      ns = 0
      do ip=1,np
        do i=1,ndim
          pos(i) = pv(i,ip)
        enddo
        call getvar(pos, nx, ny, nz, xx, yy, zz, fv, f, ierr)
        if (abs(f) <= surface_proximity) ns = ns + 1
      enddo

      end

      subroutine extract_surface_nodes (
     &        np, pv, ! points in a volume
     &        nx, ny, nz, xx, yy, zz, fv, ! scalar field
     &        ns, ps, ss ! points and scalars on the surface
     & )

      implicit none

      integer*4 ndim, np, ns, nx, ny, nz, ip, i, ierr
      parameter(ndim=3)

      real*8 pv(ndim,np), ps(ndim,ns), fv(np), ss(ns),
     & xx(nx), yy(ny), zz(nz), field(nx,ny,nz),
     & pos(ndim), f

      real*8 surface_proximity
      common/surface/surface_proximity

      ns = 0
      do ip=1,np
        do i=1,ndim
          pos(i) = pv(i,ip)
        enddo
        call getvar(pos, nx, ny, nz, xx, yy, zz, fv, f, ierr)
        if (abs(f) <= surface_proximity) then
          ns = ns + 1
          do i=1,ndim
            ps(i,ns) = pv(i,ip)
          enddo
          ss(ns) = f
        endif
      enddo

      end

!----------------------------------------
!
! Utilities for Box Debugging 
!
      subroutine box_stat (
     & nx,ny,nz,
     & np,ipp
     & )
      implicit none
      integer*4 ix,iy,iz,nx,ny,nz,np,ip,m,n,
     & ipp(3,np+2,nx,ny,nz), 
     & pind, prev, next, head, free,
     & inext, iprev
      
      parameter(pind=1, prev=2, next=3, head=1, free=2)

      write(*,*)'>>>>> BOX FILLUP:',nx,ny,nz,np

      m=0
      do iz=1,nz
      do iy=1,ny
      do ix=1,nx
      n = ipp(pind,head,ix,iy,iz)
      if (n>0) write(*,*)ix,iy,iz,n
      m = m + n
      enddo
      enddo
      enddo
      n = nx*ny*nz
      write(*,*)'Total particles: ',m,' for ',n,'boxes, utilization ',
     & real(m)/real(n)*100,'%'
      end
!
! Show box information
! (for debugging)
!
      subroutine box_info (
     & ix,iy,iz,
     & nx,ny,nz,
     & mpb,ipp
     & )
      implicit none
      integer*4 ix,iy,iz,nx,ny,nz,mpb,ib,jb,i,ip,n,
     & ipp(3,mpb+2,nx,ny,nz), 
     & pind, prev, next, head, free,
     & inext, iprev
      
      parameter(pind=1, prev=2, next=3, head=1, free=2)

      write(*,*)'>>>>> BOX ',ix,iy,iz
      do i=1,2
      if (i==head) then
        write(*,*)'PARTICLES '
      else
        cycle !-
        write(*,*)'FREE      '
      endif
      ib = i
      n = 1
      do 
        ip = ipp(pind,ib,ix,iy,iz)
        inext = ipp(next,ib,ix,iy,iz) 
        iprev = ipp(prev,ib,ix,iy,iz)
        write(*,*)n,ip,ib,iprev,inext
        ib = ipp(next,ib,ix,iy,iz)
        if (ib == i) exit
        n = n+1
      enddo
      enddo
      end
!
! Show raw box contents
! (for debugging)
!
      subroutine raw_box_info (
     & ix,iy,iz,
     & nx,ny,nz,
     & mpb,ipp
     & )
      implicit none
      integer*4 ix,iy,iz,nx,ny,nz,mpb,ib,jb,i,ip,n,
     & ipp(3,mpb+2,nx,ny,nz), inext, iprev, np, nf,
     & pind, prev, next, head, free
      
      parameter(pind=1, prev=2, next=3, head=1, free=2)

      write(*,*)'>>>>> RAW BOX ',ix,iy,iz
      np = ipp(pind,head,ix,iy,iz)
      iprev = ipp(prev,head,ix,iy,iz)
      inext = ipp(next,head,ix,iy,iz)
      write(*,*)'np=',np,iprev,inext
      nf = ipp(pind,free,ix,iy,iz)
      iprev = ipp(prev,free,ix,iy,iz)
      inext = ipp(next,free,ix,iy,iz)
      write(*,*)'nf=',nf,iprev,inext

      do i=3,mpb+2
        ip    = ipp(pind,i,ix,iy,iz)
        iprev = ipp(prev,i,ix,iy,iz)
        inext = ipp(next,i,ix,iy,iz)
        write(*,*)i,ip,iprev,inext
      enddo
      end

!________________________________________________________
!
! Compute normals to the surface
!________________________________________________________

      subroutine surface_normals (
     & nx, ny, nz, xx, yy, zz,
     & grad, nn, nodes, 
     & nbx, nby, nbz, mpb, 
     & ipb, ipx, ipy, ipz, 
     & ipp, snorm)

      implicit none

      integer*4 ndim ! number of nodes given
      parameter(ndim=3)

      integer*4 
     & nn, n,    ! number of nodes and loop counters
     & nx, ny, nz, ! grid dimensions
     & ix, iy, iz, jx, jy, jz, i, j, k, l,
     & ierr
      real*8 xx(nx), yy(ny), zz(nz), 
     & nodes(ndim,nn), xmin(ndim), xmax(ndim), 
     & r, f, x, e, e0, e1, pos(ndim), 
     & dmin(ndim), dmax(ndim), snorm(ndim,nn),
     & dist, dist_max, dist_min, move_max, move_min

! Domain dimensions and scales:
      common/geo/ dmin, dmax, 
     & dist, dist_max, dist_min, move_max, move_min
      integer*4 ipb(nn), ipx(nn), ipy(nn), ipz(nn),
     & ipp(ndim, mpb+2, nbx, nby, nbz)
      integer*4 nbx, nby, nbz, mpb, 
     &          ibx, iby, ibz, ibn

! Surface projection variables:
      real*8 g, q, u, d(3), qold, qnew, 
     & vert(2,2,2), grad(ndim,nx,ny,nz)
      real*8 surface_proximity
      common/surface/surface_proximity
!
! Looping over all nodes
!
      do n=1,nn

        do i=1,ndim ! remember old position
          pos(i)=nodes(i,n)
        enddo
        do i=1,ndim
          x = pos(i)
          if (x <= dmin(i) .or. x >= dmax(i)) exit
        enddo
        if (i < ndim) cycle ! outside domain -> freeze

        call getixyz (pos, 
     &  nx, ny, nz, 
     &  xx, yy, zz, 
     &  ix, iy, iz, 
     &  xmin, xmax, 
     &  ierr)

        if (ierr > 0) then
          cycle
        endif

        g = 0.0
        do l=1,3

        do k=0,1
        do j=0,1
        do i=0,1
        vert(i+1,j+1,k+1) = grad(l,ix+i,iy+j,iz+k)
        enddo
        enddo
        enddo

        call interp3 ( pos, xmin, xmax, vert, u )

        d(l) = u
        g = g + u*u

        enddo !l

        g = sqrt(g)

        do i=1,ndim
          snorm(i,n) = d(i)/g
        enddo

      enddo ! next n
      return
      end

!=================================================
! Implementation of Persson's algorithm:
! persson.berkeley.edu/distmesh
! with triangualtion by Burkardt:
! people.sc.fsu.edu/~jburkardt/f_src/table_tet_mesh/
!
! EDGE EXTRACTION ROUTINES
!
! Appends a new edge to the list of extracted edges
!
      subroutine append_edge ( n1, n2, nt, tree, nl, list )
      integer*4 n1, n2, nt, tree(3,nt), nl, list(3,nl)
      integer*4 none_ind,edge_ind,next_ind,left_ind,right_ind
      common/indices/none_ind,edge_ind,next_ind,left_ind,right_ind
!     Insert the first edge with the first index n1:
      nt = nt+1
      nl = nl+1
      tree( edge_ind , nt ) = nl
      tree( left_ind , nt ) = none_ind
      tree( right_ind, nt ) = none_ind
      list( 1        , nl ) = n1 ! index of node 1
      list( 2        , nl ) = n2 ! index of node 2
      list( next_ind , nl ) = none_ind
      end
! 
! Inserts the new edge into a sorted binary tree.
!
      subroutine insert_edge ( n1, n2, mt, nt, tree, nl, list )
      integer*4 n1, n2, mt, nt, tree(3,nt), nl, list(3,nl),
     & i, n, edge, next
      integer*4 none_ind,edge_ind,next_ind,left_ind,right_ind
      common/indices/none_ind,edge_ind,next_ind,left_ind,right_ind
!     mt: maximum limit on the size of the tree
!     nt: actual size of the tree
!      if ( n1 == n2 ) then
!        return
!      endif
!      if ( n1 > n2 ) then
!        write(*,*)'WARNING: SWAP INDICES',n1,n2
!        i = n1
!        n1 = n2
!        n2 = i
!      endif
!      if (n1 < 1 .or. n2 < 1) then
!         write(*,*) 'ERROR: Insert Edge: ',n1,',',n2
!         return
!      endif
!      write(*,*)'LIST of nl=',nl,':'
      if ( nt == 0 ) then
        call append_edge ( n1, n2, nt, tree, nl, list )
        return
      endif
      i = 1
      do n = 1, mt
        edge = tree ( edge_ind, i )
        i1 = list ( 1, edge )
        if ( n1 < i1 ) then 
          next = tree ( left_ind, i )
          if ( next == none_ind ) then
            call append_edge ( n1, n2, nt, tree, nl, list )
            tree ( left_ind, i ) = nt
            return
          else
            i = next
            cycle
          endif
        endif
        if ( n1 > i1 ) then 
          next = tree ( right_ind, i )
          if ( next == none_ind ) then
            call append_edge ( n1, n2, nt, tree, nl, list )
            tree ( right_ind, i ) = nt
            return
          else
            i = next
            cycle
          endif
        endif
        exit
      enddo
      if ( n == mt ) then
        write(*,*)'ERROR: Exceeded tree capacity of ', mt
        return
      endif
!     Now n1 == i1: first edge node is found
!     Check the second edge node:
!      write(*,*)'DO Traverse list n2=',n2
      do n=0,nl
        i2   = list ( 2, edge )
        next = list ( next_ind, edge )
        if ( next .ne. none_ind .and. i2 .ne. n2 ) then
          edge = next
        else
          exit
        endif
      enddo
      if ( i2 == n2 ) then
        return
      endif
      if ( n == nl ) then
        return
      endif
!     Add edge to the list:
      nl = nl + 1
      list ( next_ind, edge ) = nl
      list ( 1       , nl   ) = n1
      list ( 2       , nl   ) = n2
      list ( next_ind, nl   ) = none_ind
      end

!
! Extracts all edges from the 3D mesh obtained from triangulation routines.
! The edges are extracted into the sorted binary tree to increase speed and 
! drop the dupplicates.
!
! The triangulation routine returns the list of all triangular faces of the mesh. 
! The list of edges should be extracted from the list of faces. Since the same 
! edge can be referenced from multiple faces such dupplicates should be dropped. 
! A time efficient way of achieving this is to successively extract each edge into 
! a sorted binary tree. This enables checking of dupplicates on-the-fly in log-time.
!

      subroutine extract_edges ( node_num,
     & max_edge_num, edge_num, edges, 
     & tree, fc_num, fc, vm ) 
      integer*4 node_num, 
     & max_edge_num, edge_num, edges(3, max_edge_num),
     & tree(3, max_edge_num), fc_num, fc(7, fc_num), vm(node_num)
      integer*4 tree_num, i, j, k, n1, n2, edge, next, error
      integer*4 none_ind,edge_ind,next_ind,left_ind,right_ind
      common/indices/none_ind,edge_ind,next_ind,left_ind,right_ind
      none_ind = 0
      edge_ind = 1
      left_ind = 2
      right_ind = 3
      next_ind = 3
!
! Initialize edges and the tree
!
      do i=1,max_edge_num
      do j=1,3
      edges(j,i) = none_ind
      tree(j,i) = none_ind
      enddo
      enddo
      write(*,*)'Extracting edges:'
      write(*,*)'node_num=',node_num
      write(*,*)'max_edge_num=',max_edge_num
      write(*,*)'edge_num=',edge_num
      write(*,*)'fc_num=',fc_num
      edge_num = 0
      tree_num = 0
      error = 0
      do i = 1, fc_num
!        write(*,*)'i=',i
        if (fc(1, i) <= 0 .and. fc(2, i) <= 0) then
          cycle
        endif
        do j = 1, 3
          k = mod(j, 3) + 1
          if (k > 1) then
             n1 = vm(fc( j, i ))
             n2 = vm(fc( k, i ))
!?             n1 = fc( j, i )
!?             n2 = fc( k, i )
          else
             n1 = vm(fc( k, i ))
             n2 = vm(fc( j, i ))
!?             n1 = fc( k, i )
!?             n2 = fc( j, i )
          endif
          if (n1 < 1 .or. n2 < 1
     & .or. n1 > node_num .or. n2 > node_num ) then
             write(*,*)'WARNING: Wrong indexes ',n1,', ',n2
!?             error = 2
             cycle             
          endif
          if ( edge_num == max_edge_num ) then
            write(*,*)'ERROR: Reached the maximum edge limit of ',
     &                 edge_num
            error = 1
            exit
          endif
          call insert_edge 
     &    ( 
     &      n1, n2, max_edge_num, 
     &      tree_num, tree,
     &      edge_num, edges
     &    )
        enddo
        if (error > 0) then
           exit
        endif 
      end do
      end

      subroutine print_list ( n, list )
      integer*4 n, list(3,*), i
      integer*4 none_ind,edge_ind,next_ind,left_ind,right_ind
      common/indices/none_ind,edge_ind,next_ind,left_ind,right_ind
      do i=1,n
        write(*,*)list( 1, i ), list( 2, i )
      enddo
      end
 
      subroutine print_sorted ( tree, list )
      integer*4 tree(3,*), list(3,*), next, edge
      integer*4 none_ind,edge_ind,next_ind,left_ind,right_ind
      common/indices/none_ind,edge_ind,next_ind,left_ind,right_ind
      parameter(nstack=1000)
      integer*4 stack(nstack), tail
      next = 1
      tail = 1
      do
        if ( tree(left_ind, next) .ne. none_ind ) then
          if (tail == nstack) then
            write(*,*)'ERROR: Stack size of ',nstack,' exceeded'
            return
          endif
          stack(tail) = next
          tail = tail + 1
          next = tree(left_ind, next)
          cycle
        endif
        edge = tree(edge_ind, next)
        do 
          n1 = list(1, edge)
          n2 = list(2, edge)
          edge = list( next_ind, edge )
          if (edge == none_ind) then
            exit
          endif
        enddo
        if ( tree(right_ind, next) .ne. none_ind ) then
          if (tail == nstack) then
            write(*,*)'ERROR: Stack size of ',nstack,' exceeded'
            return
          endif
          stack(tail) = next
          tail = tail + 1
          next = tree(right_ind, next)
          cycle
        endif
        if (tail > 1) then
          tail = tail - 1
          next = stack(tail)
        endif
      enddo
      end 

!
! Routine to compute the force between the two particles.
!
      subroutine force ( r, f )
      parameter(force_coeff=0.1)
      real*8 r, f, p
      real*8 dmin(3), dmax(3), 
     & dist, dist_max, dist_min, move_max, move_min
      common/geo/ dmin, dmax, 
     & dist, dist_max, dist_min, move_max, move_min
      p = max(r, dist_min)
      f = force_coeff * min(dist, 1./p)
      end

      subroutine save_vector (nvec, vec, buf)
      parameter(ndim=3)
      integer*4 nvec, i, j
      real*8 vec( ndim, nvec ), buf( ndim, nvec )
      do i=1,nvec
         do j=1,ndim
            buf(j, i) = vec(j, i)
         enddo
      enddo
      end
!
! Get grid indexes for a given position
!
      subroutine getixyz (pos, 
     & nx, ny, nz, 
     & xx, yy, zz,
     & ix, iy, iz, 
     & xmin, xmax,
     & ierr)
      integer*4 nx, ny, nz, ix, iy, iz, ierr
      real*8 pos(3), xx(nx), yy(ny), zz(nz),
     & xmin(3), xmax(3)
      ierr = 0
      ix = ind(pos(1), nx, xx) - 1
      if ( ix <= 0 .or. ix >= nx) then
         write(*,*)'WARNING: particle moved out of range in X'
         ierr = 1
         return
      endif
      iy = ind(pos(2), ny, yy) - 1
      if ( iy <= 0 .or. iy >= ny) then
         write(*,*)'WARNING: particle moved out of range in Y'
         ierr = 2
         return
      endif
      iz = ind(pos(3), nz, zz) - 1
      if ( iz <= 0 .or. iz >= nz) then
         write(*,*)'WARNING: particle moved out of range in Z'
         ierr = 3
         return
      endif
      if (ix*iy == 1 .or. ix*iz == 1 .or. iy*iz == 1) then
         write(*,*)'WARNING: particle at the domain corner:'
         ierr = 4
         return
      endif
      xmin(1) = xx(ix)
      xmax(1) = xx(ix+1)
      xmin(2) = yy(iy)
      xmax(2) = yy(iy+1)
      xmin(3) = zz(iz)
      xmax(3) = zz(iz+1)
      return
      end
!
! Interpolate variable 
!
      subroutine getvar ( pos,
     & nx, ny, nz, xx, yy, zz,
     & field, u, ierr)

      integer*4 nx, ny, nz, ix, iy, iz, i, j, k, ierr
      real*8 pos(3), field(nx,ny,nz), vert(2,2,2), u, 
     & xx(nx), yy(ny), zz(nz), xmin(3), xmax(3)

      call getixyz (pos, 
     & nx, ny, nz, 
     & xx, yy, zz,
     & ix, iy, iz, 
     & xmin, xmax, 
     & ierr)
      if (ierr > 0) then
        write(*,*)'WARNING: getvar indexes out of range' !-
        return
      endif
      do i=0,1
      do j=0,1
      do k=0,1
      vert(i+1,j+1,k+1) = field(ix+i,iy+j,iz+k)
      enddo
      enddo
      enddo
  
      call interp3 ( pos, xmin, xmax, vert, u )
  
      return
      end

!
! MOVES PARTICLES BY FORCES
! Algorithm of Persson:
! persson.berkeley.edu/distmesh/
!
      subroutine push (
     & nx, ny, nz,
     & edge_num, node_num, 
     & edges, nodes, forces, 
     & nodeq, !-
     & xx, yy, zz, field, grad
     & )
!
! Moves particles by applying the forces 
! The forces are computed between the particles 
! connected by a common edge. 
! In the first part all the particle forces are 
! assembled by looping through all edges. 
! Then each particle is moved by the applied force.
!
      parameter(ndim=3,small=1.e-30)
      integer*4 nx, ny, nz, ix, iy, iz,
     &          edge_num, node_num,
     &          edges(3, edge_num)
      real*8 nodes(ndim, node_num), forces(ndim, node_num),
     & nodeq(node_num), !-
     & xx(nx), yy(ny), zz(nz), xmin(3), xmax(3),
     & field(nx,ny,nz), grad(ndim,nx,ny,nz), 
     & d(3), pos(3), vert(2,2,2), r, x, x1, x2, f, 
     & g, q, u
      real*8 dmin(3), dmax(3), 
     & dist, dist_max, dist_min, move_max, move_min
      common/geo/ dmin, dmax, 
     & dist, dist_max, dist_min, move_max, move_min
      integer*4 i,j,k,l,m,n,i1,i2,ierr
      do i=1,node_num
         do j=1,ndim
            forces(j, i) = 0.0
         enddo
      enddo
!     Assemble forces
      do i=1,edge_num
        i1 = edges( 1, i ) ! index of node 1
        i2 = edges( 2, i ) ! index of node 2
!#ifdef DEBUG
        if (i1 > node_num .or. i2 > node_num
     &      .or. i1 < 1 .or. i2 < 1) then
           write(*,*)'ERROR: Wrong indexing in move:'
           write(*,*)'edge ',i,' of ',edge_num
           write(*,*)'i1=',i1,', i2=',i2,' of ',node_num
           return
        endif
!#endif
        r = 0.0
        do j=1,ndim
          x1 = nodes( j, i1 )
          x2 = nodes( j, i2 )
          x = x2 - x1
          r = r + x*x
          d(j) = x
        enddo
        r = sqrt(r)
        do j=1,ndim
          target_cell_size = 1.5*dist
          x = max(target_cell_size - r, 0.0)
          f = x * d(j) / r
          forces(j, i1) = forces(j, i1) - f
          forces(j, i2) = forces(j, i2) + f
        enddo
      enddo

!     Move nodes
      nmov=0 !- moved by force
      nprj=0 !- normal projected
      mprj=0 !- max projected
      nfrz=0 !- number frozen
      nerr=0 !-

      write(*,*)'Looping over ',node_num,' nodes'

      do n=1,node_num

        do i=1,ndim
           pos( i ) = nodes( i, n )
        enddo

        q = nodeq(n)

        if (q > 0.0) then ! in the positive field: moving

          f = 0.0 !-
          do i=1,ndim
             r = forces( i, n ) !-
             pos( i ) = pos( i ) + r
             f = f + r*r
          enddo

          call getvar (pos, 
     &    nx, ny, nz, xx, yy, zz,
     &    field, u, ierr)

          if (ierr > 0) then
            nfrz = nfrz + 1
            write(*,*)'Particle ',n,' stuck' !-
            cycle ! particle stuck at the grid boundary
          endif    
  
          nodeq(n) = u
          nmov=nmov+1 !-

        else ! in the negative field: projecting

! Move particle to zero-filed level along the direction
! of the gradient

         iprj = 0 !-
         jprj = 0 !-


           call getixyz (pos, 
     &     nx, ny, nz, 
     &     xx, yy, zz, 
     &     ix, iy, iz, 
     &     xmin, xmax, 
     &     ierr)

         g = 0.0
         do l=1,3

           do k=0,1
           do j=0,1
           do i=0,1
           vert(i+1,j+1,k+1) = grad(l,ix+i,iy+j,iz+k)
           enddo
           enddo
           enddo

           call interp3 ( pos, xmin, xmax, vert, f )

           d(l) = f
           g = g + f*f

         enddo !l

         ! Project onto the zero surface:
         do i=1,ndim
           pos(i) = pos(i) + 0.5 * f * d(i) / g
         enddo

         call getvar (pos, 
     &   nx, ny, nz, xx, yy, zz,
     &   field, u, ierr)

         if (ierr > 0) then
           nfrz = nfrz + 1 !-
           nprj = nprj - iprj
           mprj = mprj - jprj
           nerr = nerr + 1
       write (*,*)'Particle ',n, ' stuck at the boundary' !-
           cycle ! particle stuck at the grid boundary
         endif    

         nodeq(n) = u
       write(*,*)n,'Particle u=',u,' nprj=',nprj,' mprj=',mprj !-

        endif
! Update node positions
        do i=1,ndim
           nodes( i, n ) = pos(i)
        enddo
      enddo ! next n
      write(*,*)'MOVE: moved=',nmov,' nprj=',nprj,
     &' mprj=',mprj,' nfrz=',nfrz,' nerr=',nerr !-
      return
      end

!
! Compute the number of displaced particles and the maximum displacement. 
! The displacement is measured by comparing with the position of the 
! particle before the particles were moved. This requires to store a copy 
! of the particles initial positions.
! Computing the displacement is needed to decide of whether to continue the 
! particle moving loop or run re-triangulation. When at least one particle
! is moved over the pre-defined limit the mesh is re-generated.
!

      subroutine displacement(num, old, new, nd, del, dmx)
      parameter(ndim=3)
      real*8 dmin(3), dmax(3), 
     & dist, dist_max, dist_min, move_max, move_min
      common/geo/ dmin, dmax, 
     & dist, dist_max, dist_min, move_max, move_min
      integer*4 num, i, j, nd
      real*8 old( ndim, num), new( ndim, num),
     & d, del
      dmx = dist_max
      nd = 0
      del = 0.0
      do i=1,num
         d = 0.0
         do j=1,ndim
            d = d + (new(j,i) - old(j,i))**2
         enddo
         d = sqrt(d)
         if (d > dist_max) then
            nd = nd + 1
         endif
         del = del + d
      enddo
      if (nd > 0) then
         del = del / real(num)
      endif
      end

      subroutine get_unit ( iunit )
      implicit none
      integer ( kind = 4 ) i
      integer ( kind = 4 ) ios
      integer ( kind = 4 ) iunit
      logical lopen
      iunit = 0
      do i = 1, 99
      if ( i /= 5 .and. i /= 6 .and. i /= 9 ) then
        inquire ( unit = i, opened = lopen, iostat = ios )
        if ( ios == 0 ) then
          if ( .not. lopen ) then
            iunit = i
            return
          end if
        end if
      end if
      end do
      return
      end

