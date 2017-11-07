!  This program is designed to generate surface mesh based on
!  the level-set result
!
!  Written by Andrei from 07/11/2017
!
!  This is D&P's intellectual property. Do not distribute it outside of
!  the company!
!______________________________________________________________________
!
PROGRAM MAIN

      IMPLICIT NONE
      INTEGER  IMAX,JMAX,KMAX,i
      REAL*8,ALLOCATABLE:: X(:),Y(:),Z(:),PHI(:,:,:)
!
!..Read in the background mesh and level sets
!
      WRITE(*,*)'READING DATA'
      CALL INPUT(X,Y,Z,PHI,IMAX,JMAX,KMAX)

!      write(*,*)'After INPUT: X:'
!      do i=1,IMAX
!         write(*,*)i,X(i)
!      enddo
!      write(*,*)'Y:'
!      do i=1,JMAX
!         write(*,*)i,Y(i)
!      enddo
!      write(*,*)'Z:'
!      do i=1,KMAX
!         write(*,*)i,Z(i)
!      enddo
!
!..Obtain surface points by normal projection
!
      CALL SURFER(IMAX,JMAX,KMAX,X,Y,Z,PHI)
!
!..Triangulate surface points with Delauney algorithm
!
!      CALL TRIANG()
!
      DEALLOCATE(X,Y,Z,PHI)
!
      STOP

contains

!______________________________________________________________________
!
      SUBROUTINE INPUT(X,Y,Z,PHI,IMAX,JMAX,KMAX)
!
!     Read in the background mesh and level sets
!______________________________________________________________________
!
      IMPLICIT NONE
      INTEGER*4 I,IMAX,J,JMAX,K,KMAX,N
      REAL*8,ALLOCATABLE:: X(:),Y(:),Z(:),PHI(:,:,:)
!
      READ(8,*) IMAX,JMAX,KMAX

      WRITE(*,*)'GRID DIMENSIONS:',IMAX,JMAX,KMAX

      ALLOCATE(X(IMAX),Y(JMAX),Z(KMAX))
!      READ(8,*)(((X(I),I=1,IMAX),J=1,JMAX),K=1,KMAX), &
!               (((Y(J),I=1,IMAX),J=1,JMAX),K=1,KMAX), &
!               (((Z(K),I=1,IMAX),J=1,JMAX),K=1,KMAX)

! The input was restructued as a rectilinear grid to save space:
      READ(8,*)(X(I),I=1,IMAX), &
               (Y(J),J=1,JMAX), &
               (Z(K),K=1,KMAX)
      
      CLOSE(8)

      READ(9,*) IMAX,JMAX,KMAX
      ALLOCATE(PHI(IMAX,JMAX,KMAX))
      READ(9,*)(((PHI(I,J,K),I=1,IMAX),J=1,JMAX),K=1,KMAX)
      CLOSE(9)

      RETURN
      END SUBROUTINE

!______________________________________________________________________
!
! CREATE SURFACE MESH
!______________________________________________________________________
!
      SUBROUTINE SURFER(NX,NY,NZ,X,Y,Z,Q)

      IMPLICIT NONE

      INTEGER*4 NDIM, NX,NY,NZ,NP,IERR
      PARAMETER(NDIM=3)
      REAL*8 X(NX),Y(NY),Z(NZ),Q(NX*NY*NZ)
      REAL*8,ALLOCATABLE:: & 
             P(:,:),       & ! Points coordinates
             F(:,:),       & ! Forces at points
             S(:),         & ! Scalar at the points
             G(:,:,:,:)      ! Gradient of Q


! NODE PLACEMENT VARIABLES:
      INTEGER*4 MAX_ITER, ITER, & 
      NDISP, & ! number of consecutive iterations with 
               ! displacement below minimum to satisfy convergence
      NMOVES, & ! number of consecutive iteration satisfying
                ! max displaement criterion for convergence
      NFAIL, IFAIL ! number of iterations with displacements
                   ! above the minimum
      REAL*8 ENG, EMIN, EMAX, ER, SIGMA, DISP_MIN, &
      CONVER, & ! convergence factor based on maximum displacment
      DISP      ! max displacement computed in one iteration
      COMMON/POTENTIAL/SIGMA,EMIN,EMAX

      PARAMETER(MAX_ITER=999, NDISP=1000, CONVER=0.1)

      REAL*8 DMIN(3), DMAX(3), &
      DIST, DIST_MAX, DIST_MIN, MOVE_MAX, MOVE_MIN
      COMMON/GEO/ DMIN, DMAX, &
      DIST, DIST_MAX, DIST_MIN, MOVE_MAX, MOVE_MIN


! ACCELERATION VARIABLES:
      REAL*8 DB, BMIN(3), BMAX(3)
      INTEGER*4 NBX, NBY, NBZ, MPB
      INTEGER*4,ALLOCATABLE:: IPB(:), IPX(:), IPY(:), IPZ(:), &
                              IPP(:, :, :, :, :)  
      COMMON/PB/ DB, BMIN, BMAX

! SURFACE VARIABLES
      INTEGER*4          & 
        SURFACE_PROJECTION, & ! use surface projection
        NS                    ! Number of surface points
      REAL*8,ALLOCATABLE:: &
        PS(:,:),           &  ! Surface points coordinates 
        SN(:,:),           &  ! Surface normal vectors 
        SS(:)                 ! Surface scalar coordinates 

      WRITE(*,*)'GRID:',NX,NY,NZ

 CALL INIT_GRID ( NX, NY, NZ, X, Y, Z )

      NP = 0

 CALL MAKE_NODES (   & ! If NP=0 then
        NX, NY, NZ,  & ! this will return NP
        X, Y, Z,     & ! so we can allocate here
        Q, NP, P, S  & ! This is because F77
        )              ! cannot allocate

       ALLOCATE(P(NDIM,NP))
       ALLOCATE(S(NP))

 CALL MAKE_NODES (   & ! now create points
        NX, NY, NZ,  &
        X, Y, Z,     &
        Q, NP, P, S  &
      )

 
 CALL WRITE_VTP ('initial_points.vtp', NP, P, S)


! Initialize particle boxing for quasi-linear scale-up in NP


 CALL ALLOC_PB ( NP, NBX, NBY, NBZ, MPB,           &
                     IPB, IPX, IPY, IPZ, IPP, IERR )

 CALL INIT_PB  ( NP, P,                         &
                 NBX, NBY, NBZ, MPB,            &
                 IPB, IPX, IPY, IPZ, IPP, IERR )

      IF (IERR > 0) THEN
        WRITE(*,*)'ERROR INITILIZING PARTICLE BOXES'
        RETURN
      ENDIF

!
! MOVE NODES TO OBTAIN UNIFORM DISTRIBUTION
!

WRITE(*,*)'ADJUSTING NODES TO MINIMIZE VOLUME ENERGY:'

      DISP_MIN = CONVER * DIST 
WRITE(*,*)'MINIMUM DISPLACEMENT FOR CONVERGENCE: '
WRITE(*,*)'   DISP_MIN=',CONVER,'*',DIST,'=',DISP_MIN
      NFAIL = 0

      DO ITER=1,MAX_ITER

 CALL   MOVE_NODES (  & ! moving volume nodes
           NX,  NY,  NZ, X, Y, Z, &
            Q,  NP,   P, S,       &
          NBX, NBY, NBZ, MPB,     &
          IPB, IPX, IPY, IPZ,     &
          IPP, ENG, NMOVES, DISP )

        WRITE(*,*)ITER,'MOVES:',NMOVES,'DISP:',DISP/DISP_MIN,'FAIL:',NFAIL,NDISP-IFAIL,'ENERGY: ',ENG

        IF (DISP > DISP_MIN) NFAIL = NFAIL + 1
        IFAIL = MOD(ITER,NDISP)
        IF (IFAIL == 0) THEN
          IF(NFAIL == 0) THEN
            WRITE(*,*)'CONVERGED AFTER',NDISP,' ITEREATIONS WITH MAX DISPLACEMENT BELOW: ',DISP_MIN
            EXIT
          ENDIF
          NFAIL = 0
        ENDIF

      ENDDO ! NEXT ITER


! Output the nodes to vtp file for Paraview:
 CALL WRITE_VTP ( 'points.vtp', NP, P, S )

! Output the nodes to dat file for Qhull:
 CALL WRITE_DAT ( 'points.dat', NP, P )

! Output the nodes to file for Tetgen:
 CALL WRITE_NODE ( 'surf.node', NP, P )

!
! COMPUTE SURFACE NORMAL VECTORS
! This is required for frontmesh

!
! PROJECT SURFACE NODES
! 
! Apply normal projection procedure to place nodes 
! on the surface.
! This step is not necessary, since the previous step
! already produces a good mesh. However, with a more
! refined level-set function it may be useful.

WRITE(*,*)'EXTRACTING SURFACE NODES'
!
! First count all surface particles:
!
 CALL COUNT_SURFACE_NODES (     &
        NP,  P,                 & ! all nodes
        NX, NY, NZ, X, Y, Z, Q, & ! scalar field
        NS                      & ! returns number of 
                                  ! nodes on the surface
      )

! 
! Allocate surface nodes:
! PS: node coordinates
!
      ALLOCATE(PS(NDIM,NS))
      ALLOCATE(SS(NS))

!
! Extract surface particles from all particles
!
 CALL EXTRACT_SURFACE_NODES (   &
        NP,  P,                 & ! all nodes
        NX, NY, NZ, X, Y, Z, Q, & ! scalar field
        NS, PS, SS              & ! extracted surface nodes
      )

WRITE(*,*)'EXTRACTED ',NS,' SURFACE NODES'
! Output the nodes to file for Paraview:
 CALL WRITE_VTP ('extracted_points.vtp', NS, PS, SS)
!
! Compute gradient of Q:
!
    ALLOCATE (G(NDIM,NX,NY,NZ))

 CALL GRADIENT( NX, NY, NZ, X, Y, Z, Q, G )

WRITE(*,*) 'INITIALIZING PARTICLE BOXES'

 CALL INIT_PB (  NS,  PS,                      &
                NBX, NBY, NBZ, MPB,            &
                IPB, IPX, IPY, IPZ, IPP, IERR )

WRITE(*,*)'PROJECTING SURFACE NODES'

 CALL PROJECT_SURFACE_NODES (       &
         NX,  NY,  NZ,   X,  Y, Z,  &
          Q,   G,  NS,  PS, SS,     &
        NBX, NBY, NBZ, MPB,         &
        IPB, IPX, IPY, IPZ,         &
        IPP, ENG )

 CALL WRITE_VTP ('projected_points.vtp', NS, PS, SS)
! Output nodes to file for Tetgen and Frontmesh:
 CALL WRITE_NODE ( 'surf.node', NS, PS )
! Output nodes to file for Qhull and Frontmesh:
 CALL WRITE_DAT ( 'surf.dat', NS, PS )

WRITE(*,*)'COMPUTING SURFACE NORMALS'

      ALLOCATE(SN(NDIM,NS))

 CALL SURFACE_NORMALS (       &
         NX,  NY,  NZ,   X,  Y, Z,  &
          G,  NS,  PS,              &
        NBX, NBY, NBZ, MPB,         &
        IPB, IPX, IPY, IPZ,         &
        IPP, SN )

! Output surface normals for frontmesh in node-format:
 CALL WRITE_NODE ( 'snorm.node', NS, SN )
! Output normals for Frontmesh in dat format:
 CALL WRITE_DAT ( 'snorm.dat', NS, SN )
! Output surface nomrals for Paraview
 CALL WRITE_VECTOR_VTP ('snorm.vtp', NS, NDIM, PS, SN)

      DEALLOCATE(SN)
      DEALLOCATE(G)
      DEALLOCATE(PS, SS)
      DEALLOCATE(IPB, IPX, IPY, IPZ, IPP)
      DEALLOCATE(S, P)

      END
!______________________________________________________________________
!
! WRITE POINTS IN TETGEN FORMAT  
!______________________________________________________________________
!
      SUBROUTINE WRITE_NODE (FILENAME, NP, P )

      IMPLICIT NONE
      CHARACTER*(*) FILENAME
      INTEGER*4 NP, I
      REAL*8 P(3,NP)
      WRITE(*,*)'OUTPUT TETGEN DATA TO ',FILENAME
      OPEN(1,FILE=FILENAME)
      WRITE(1,*) NP,3,0,0
      DO I=1,NP
        WRITE(1,*)I,P(1,I),P(2,I),P(3,I)
      ENDDO
      CLOSE(1)

      END



!______________________________________________________________________
!
! MESH THE DOMAIN BY 3D DELAUNAY TRIANGULATION
!______________________________________________________________________
!
      SUBROUTINE DELAUNAY3D (FILENAME, NP, P )

      IMPLICIT NONE

      CHARACTER*(*) FILENAME

      INTEGER*4 NP, I, ERR,                     &
                MAX_NODE_EDGES, MAX_NODE_FACES, &
                FC_MAX, BF_MAX,                 &
                FC_NUM, BF_NUM, HT_NUM,         &
                FACE_NUM, CELL_NUM
      INTEGER*4,ALLOCATABLE:: HT(:), VM(:), BF(:,:), FC(:,:)
      REAL*8 P(3,NP)

      WRITE (*,*)'GENERATE MESH FROM ',NP,' NODES'

      MAX_NODE_EDGES = 21 ! MAX EDGES CONNECTED TO A NODE

      MAX_NODE_FACES = MAX_NODE_EDGES
      FC_MAX = NP * MAX_NODE_FACES
      BF_MAX = 4*INT(REAL(FC_MAX)**(2./3.))
      HT_NUM = 8*NP 

      ALLOCATE(HT(HT_NUM))
      ALLOCATE(VM(NP))
      ALLOCATE(BF(3,BF_MAX))
      ALLOCATE(FC(7,FC_MAX))

      DO I=1,NP
        VM(I) = I
      ENDDO

 CALL DTRIS3( NP, HT_NUM, BF_MAX, FC_MAX,                &
              P, VM, BF_NUM, FC_NUM, FACE_NUM, CELL_NUM, &
              BF, FC, HT, ERR )

      IF (ERR /= 0) THEN
        WRITE(*,*)'ERROR ',ERR,': TRIANGULATION FAILED'
        RETURN
      ENDIF

  IF ( ERR /= 0 ) THEN
    WRITE ( *, '(A)' ) ' '
    WRITE ( *, '(A)' ) 'TABLE_TET_MESH - FATAL ERROR!'
    WRITE ( *, '(A,I8)' ) '  DTRIS3 RETURNED ERR = ', ERR
    WRITE ( *, '(A)' ) ' '
    WRITE ( *, '(A)' ) 'TABLE_TET_MESH'
    WRITE ( *, '(A)' ) '  ABNORMAL END OF EXECUTION!'
    STOP
  END IF

  WRITE ( *, '(A)' ) 'AFTER: '
  WRITE ( *, '(A,I8)' ) '  BF_MAX = ', BF_MAX
  WRITE ( *, '(A,I8)' ) '  BF_NUM = ', BF_NUM
  WRITE ( *, '(A)' ) ' '
  WRITE ( *, '(A,I8)' ) '  FC_MAX = ', FC_MAX
  WRITE ( *, '(A,I8)' ) '  FC_NUM = ', FC_NUM
  WRITE ( *, '(A)' ) ' '
  WRITE ( *, '(A,I8)' ) '  HT_NUM = ', HT_NUM
  WRITE ( *, '(A)' ) ' '
  WRITE ( *, '(A,I8)' ) '  TETRA_NUM = ', CELL_NUM

 CALL OUTPUT_BOUNDARY ( FILENAME, NP, FC_NUM, P, FC, VM )

      DEALLOCATE(HT)
      DEALLOCATE(VM)
      DEALLOCATE(BF)
      DEALLOCATE(FC)

      END SUBROUTINE

!______________________________________________________________________
!
! PARTICLE BOXING ROUTINES
!______________________________________________________________________
!
! Domain is split into boxes to provide a quasi-linear speedup in NP
!
! ALLOCATE PARTICLE BOXES
!
      SUBROUTINE ALLOC_PB ( NP, & ! Number of particles
                 NBX, NBY, NBZ, MPB, IPB, &
                 IPX, IPY, IPZ, IPP, IERR)

      REAL*8 DMIN(3), DMAX(3), &
      DIST, DIST_MAX, DIST_MIN, MOVE_MAX, MOVE_MIN
      COMMON/GEO/ DMIN, DMAX, &
      DIST, DIST_MAX, DIST_MIN, MOVE_MAX, MOVE_MIN
      INTEGER*4  NP, NBX, NBY, NBZ, MPB, MB, IERR
      REAL*8 DB, BMIN(3), BMAX(3)
      COMMON/PB/ DB, BMIN, BMAX
      INTEGER*4,ALLOCATABLE:: IPB(:), IPX(:), IPY(:), IPZ(:), &
                              IPP(:, :, :, :, :)   
      INTEGER*4 IX,IY,IZ,IB,IP,IBX,IBY,IBZ,       &
                PIND,PREV,NEXT,HEAD,FREE
      REAL*8 X
      PARAMETER(PIND=1,PREV=2,NEXT=3,HEAD=1,FREE=2)
!     NP: number of particles
!     MPB: number of particles in a box
!     NBX,NBY,NBZ: numbex of boxes in respective directions
!     MPB: max particles in each box
!     PTS(i,j): i-th coordinate of j-th particle
!     IPB(i): index of particle in ipc array
!     IPP: list of indexs to particles in boxes
!     i=IPP(ix,iy,iz,1,j)=index to the particle in 
!     PTS array, i.e. PTS(*,i)
!     j=location of that index in the grid box at
!     ix,iy,iz
!
!     PIND: (particle index) is the element of IPP array 
!     holding the particle index, i.e.:
!     if i=IPP(PIND,:,:,:,:) then
!     PTS(j,i) is the j-th coordinate of i-th particle
!     PREV,NEXT: pointers to next/prev particles in IPP array, ie:
!     if n1 is the index of a particle in the box at x,y,z then
!     n2=IPP(NEXT,n1,x,y,z) is the next particle to n1 
!        in the box at x,y,z. Then
!     i1=IPP(PIND,n1,x,y,z) and 
!     i2=IPP(PIND,n2,x,y,z) are the indexes to the 
!     global array of particles for these two particles, i.e.:
!     PTS(j,i1) and PTS(j,i2) are j-th coordinates of the 
!     n1-th and n2-th particles in the box at x,y,z.
!
!     PTS(j,i) is the j-th coordinate of i-th particle 
!     in a box at x,y,z. 
!     IPB(i) is the index of the i-th particle in the 
!     box at x,y,z
!     IPX/Y/Z(i) the coordinates of the box for the i-th
!     particle. For example:
!     If i is a global particle index in PTS(:,i) array, then
!     if n=IPB(i), x=IPX(i), y=IPY(i), z=IPZ(i) then
!     PTS(j,IPP(PIND,n,x,y,z)) is the j-th coordinate of 
!     the i-th particle. Thus this is an identity:
!     
!     PTS(j,IPP(PIND,IPB(i),IPX(i),IPY(i),IPZ(i))) 
!   = PTS(j,i)
!
      MPB = 40 ! MAX NUMBER OF PARTICLES IN A BOX
      MB = MPB+2 ! BOX SIZE INCLUDING SLOTS FOR 
!                 'HEAD' AND 'FREE' COUNTS

      DB = 2*DIST
      WRITE(*,*)'INIT_PB: DB=',DB
      WRITE(*,*)'BMIN/MAX:'
      DO I=1,3
! domain is extended to avoid processing 
! boundary exceptions
        BMIN(I) = DMIN(I) - DB
        BMAX(I) = DMAX(I) + DB
        WRITE(*,*)I,'D:',DMIN(I),DMAX(I)
        WRITE(*,*)I,'C:',BMIN(I),BMAX(I)
      ENDDO
      NBX = (BMAX(1)-BMIN(1))/DB + 2
      NBY = (BMAX(2)-BMIN(2))/DB + 2
      NBZ = (BMAX(3)-BMIN(3))/DB + 2

      WRITE(*,*)'NB:',NBX,NBY,NBZ

      ALLOCATE(IPB(NP))
      ALLOCATE(IPX(NP))
      ALLOCATE(IPY(NP))
      ALLOCATE(IPZ(NP))

      ALLOCATE(IPP(3,MB,NBX,NBY,NBZ))
      END

!______________________________________________________________________
!
! INITIALIZE PARTICLE BOXES
!______________________________________________________________________
!
      SUBROUTINE INIT_PB (NP, PTS,                 &
                          NBX, NBY, NBZ, MPB, IPB, &
                          IPX, IPY, IPZ, IPP, IERR)

      REAL*8 PTS(3,NP), DB, BMIN(3), BMAX(3)
      INTEGER*4 NP, NBX, NBY, NBZ, MPB, MB, IERR
      COMMON/PB/ DB, BMIN, BMAX
      INTEGER*4 IPB(NP), IPX(NP), IPY(NP), IPZ(NP), &
                IPP(3, MPB+2, NBX, NBY, NBZ)   
      INTEGER*4 IX,IY,IZ,IB,IP,IBX,IBY,IBZ,       &
                PIND,PREV,NEXT,HEAD,FREE
      REAL*8 X
      PARAMETER(PIND=1,PREV=2,NEXT=3,HEAD=1,FREE=2)
!     NP: number of particles
!     MPB: number of particles in a box
!     NBX,NBY,NBZ: numbex of boxes in respective directions
!     MPB: max particles in each box
!     PTS(i,j): i-th coordinate of j-th particle
!     IPB(i): index of particle in ipc array
!     IPP: list of indexs to particles in boxes
!     i=IPP(ix,iy,iz,1,j)=index to the particle in 
!     PTS array, i.e. PTS(*,i)
!     j=location of that index in the grid box at
!     ix,iy,iz
!
!     PIND: (particle index) is the element of IPP array 
!     holding the particle index, i.e.:
!     if i=IPP(PIND,:,:,:,:) then
!     PTS(j,i) is the j-th coordinate of i-th particle
!     PREV,NEXT: pointers to next/prev particles in IPP array, ie:
!     if n1 is the index of a particle in the box at x,y,z then
!     n2=IPP(NEXT,n1,x,y,z) is the next particle to n1 
!        in the box at x,y,z. Then
!     i1=IPP(PIND,n1,x,y,z) and 
!     i2=IPP(PIND,n2,x,y,z) are the indexes to the 
!     global array of particles for these two particles, i.e.:
!     PTS(j,i1) and PTS(j,i2) are j-th coordinates of the 
!     n1-th and n2-th particles in the box at x,y,z.
!
!     PTS(j,i) is the j-th coordinate of i-th particle 
!     in a box at x,y,z. 
!     IPB(i) is the index of the i-th particle in the 
!     box at x,y,z
!     IPX/Y/Z(i) the coordinates of the box for the i-th
!     particle. For example:
!     If i is a global particle index in PTS(:,i) array, then
!     if n=IPB(i), x=IPX(i), y=IPY(i), z=IPZ(i) then
!     PTS(j,IPP(PIND,n,x,y,z)) is the j-th coordinate of 
!     the i-th particle. Thus this is an identity:
!     
!     PTS(j,IPP(PIND,IPB(i),IPX(i),IPY(i),IPZ(i))) 
!   = PTS(j,i)
!
      MPB = 40 ! MAX NUMBER OF PARTICLES IN A BOX
      MB = MPB+2 ! BOX SIZE INCLUDING SLOTS FOR 
!                 'HEAD' AND 'FREE' COUNTS
!      NB = NP/MPB ! NUMBER OF BOXES = NBX*NBY*NBZ

!
! INITIALIZE ARRAYS
!
      DO IZ=1,NBZ
      DO IY=1,NBY
      DO IX=1,NBX
      DO IB=FREE+1,MB-1
      IPP(PIND,IB,IX,IY,IZ) = 0
      IPP(PREV,IB,IX,IY,IZ) = IB-1
      IPP(NEXT,IB,IX,IY,IZ) = IB+1
      ENDDO 
      IPP(PIND,FREE,IX,IY,IZ) = MPB ! stores 
!     the number of free particle slots
      IPP(PREV,FREE,IX,IY,IZ) = MB
      IPP(NEXT,FREE,IX,IY,IZ) = FREE+1
      IPP(PIND,MB,IX,IY,IZ) = 0
      IPP(PREV,MB,IX,IY,IZ) = MB-1
      IPP(NEXT,MB,IX,IY,IZ) = FREE
!     Assign first particle at element 1
      IPP(PIND,HEAD,IX,IY,IZ) = 0 ! number of particles
      IPP(PREV,HEAD,IX,IY,IZ) = HEAD
      IPP(NEXT,HEAD,IX,IY,IZ) = HEAD
      ENDDO
      ENDDO
      ENDDO
         
      DO IP=1,NP
        IPB(IP) = 0
        IPX(IP) = 0
        IPY(IP) = 0
        IPZ(IP) = 0
      ENDDO
! ASSIGN PARTICLES TO BOXES:
! (to boost the efficiency to quasi-linear in NP)
      DO IP=1,NP
! Determin box x,y,z indexes:
        IBX = MAX(MIN(INT((PTS(1,IP)-BMIN(1))/DB),NBX-1),2)
        IBY = MAX(MIN(INT((PTS(2,IP)-BMIN(2))/DB),NBY-1),2)
        IBZ = MAX(MIN(INT((PTS(3,IP)-BMIN(3))/DB),NBZ-1),2)

   CALL ASSIGN_BOX( IP, NP,                            &
                    NBX, NBY, NBZ, MPB,                &
                    IBX, IBY, IBZ, IPX, IPY, IPZ,      &
                    IPB, IPP, IERR                     &
        )
      IF (IERR > 0) THEN
        WRITE(*,*)'FAILED TO ASSIGN BOX'
        RETURN
      ENDIF
      ENDDO
      WRITE(*,*)'INITIAL NODES CREATED'
      END SUBROUTINE

END PROGRAM MAIN

