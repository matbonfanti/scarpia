MODULE PES
   implicit none

real,parameter :: pi = 4.0*atan(1.0)
real,parameter :: degtorad = pi/180.0

integer :: frag,syssize,ndofs

real,dimension(:),allocatable :: paramguess
real :: torsval

! real :: torsval,firstev

! real,dimension(:),allocatable :: paramarray, eigenvalues
! real,dimension(:,:),allocatable :: matrix

!
! call read_inputfiles
! call build_matrix
! call diagonalize_matrix
! call write_outputfiles
!



contains

!##################################################################################################################################

subroutine read_inputfiles

implicit none

integer :: inputfile,iter

real :: torsvalread

open(newunit=inputfile,file='input.inp',status='old',action='read')

read(inputfile,*) frag
if (frag < 2) then
    write(*,*) 'frag must be at least 2!'
    stop
end if

syssize = 2*frag
ndofs = 2*syssize-1

read(inputfile,*) torsvalread

torsval = torsvalread*degtorad



! allocate(paramarray(ndofs),paramguess(ndofs))
allocate(paramguess(ndofs))

!first (syssize)-many values for the ring-breathing modes
!second (syssize-1)-many values for the bond-stretch modes
do iter=1,ndofs
    read(inputfile,*) paramguess(iter)
end do
! paramarray(:) = 0.0

close(inputfile)



! allocate(matrix(syssize,syssize))
! allocate(eigenvalues(syssize))

! matrix = 0.0
! eigenvalues = 0.0

end subroutine read_inputfiles

!#################################################################################################################################

subroutine build_matrix( paramarray, matrix )

implicit none
real,dimension(ndofs),intent(in) :: paramarray
real,dimension(syssize,syssize),intent(out) :: matrix

integer :: iter

!real,parameter :: w=-0.0458471238

matrix=0.0

do iter=1,syssize
    if ((iter == frag) .or. (iter == frag+1)) then
        matrix(iter,iter) = matrix(iter,iter)+tor(torsval,0)+tor(torsval,1)
    else
        matrix(iter,iter) = matrix(iter,iter)+2.0*tor(torsval,0)
    end if

    if (iter == syssize) exit

    matrix(iter,iter+1) = tor(torsval,7)
    matrix(iter+1,iter) = tor(torsval,7)
!    matrix(iter,iter+1) = w
!    matrix(iter+1,iter) = w
end do



call buildbs( paramarray, matrix )
call buildint( paramarray, matrix )

end subroutine build_matrix

!#################################################################################################################################

subroutine buildbs( paramarray, matrix )

implicit none
real,dimension(ndofs),intent(in) :: paramarray
real,dimension(syssize,syssize),intent(inout) :: matrix

integer :: siteiter,dofiter
real,dimension(:),allocatable :: potbs

allocate(potbs(syssize))

potbs=0.0

do siteiter=1,syssize
    do dofiter=1,syssize-1
            potbs(siteiter)=potbs(siteiter)+2.0*bs(paramarray(dofiter+syssize),.true.)
    end do

    if (siteiter == 1) then
            potbs(siteiter)=potbs(siteiter)+bs(paramarray(siteiter+syssize),.false.)-bs(paramarray(siteiter+syssize),.true.)
    elseif (siteiter == syssize) then
            potbs(siteiter)=potbs(siteiter)+bs(paramarray(siteiter-1+syssize),.false.)-bs(paramarray(siteiter-1+syssize),.true.)
    else
            potbs(siteiter)=potbs(siteiter)+bs(paramarray(siteiter-1+syssize),.false.)-bs(paramarray(siteiter-1+syssize),.true.) &
                                           +bs(paramarray(siteiter+syssize),.false.)-bs(paramarray(siteiter+syssize),.true.)
    end if
end do

do siteiter=1,syssize
    matrix(siteiter,siteiter)=matrix(siteiter,siteiter)+potbs(siteiter)
end do

deallocate(potbs)

end subroutine buildbs

!#################################################################################################################################

subroutine buildint( paramarray, matrix )

implicit none
real,dimension(ndofs),intent(in) :: paramarray
real,dimension(syssize,syssize),intent(inout) :: matrix

integer :: siteiter,dofiter
real,dimension(:),allocatable :: potint

allocate(potint(syssize))

potint=0.0

do siteiter=1,syssize
    do dofiter=1,syssize
        potint(siteiter)=potint(siteiter)+intf(paramarray(dofiter),.true.)
    end do

    potint(siteiter)=potint(siteiter)+intf(paramarray(siteiter),.false.)-intf(paramarray(siteiter),.true.)
end do

do siteiter=1,syssize
    matrix(siteiter,siteiter)=matrix(siteiter,siteiter)+potint(siteiter)
end do

deallocate(potint)

end subroutine buildint

!##################################################################################################################################

subroutine diagonalize_matrix( matrix, eigenvalues, firstev )

implicit none
real,dimension(syssize,syssize),intent(in) :: matrix
real,dimension(syssize),intent(out) :: eigenvalues
real,intent(out) :: firstev
integer :: lwork,info
real,dimension(:),allocatable :: eigv,work
real,dimension(:,:),allocatable :: tmpmat

allocate(tmpmat(syssize,syssize))
allocate(eigv(syssize))

tmpmat = matrix
eigv = 0.0

allocate(work(1))

call DSYEV('V','U',syssize,tmpmat,syssize,eigv,work,-1,info)

lwork = ceiling(work(1))

deallocate(work)

tmpmat = matrix
eigv = 0.0

allocate(work(lwork))

call DSYEV('V','U',syssize,tmpmat,syssize,eigv,work,lwork,info)
if (info /= 0) then
    write(*,*) 'DSYEV error!'
    stop
end if

deallocate(work)



eigenvalues = eigv
firstev = eigv(1)



deallocate(tmpmat)
deallocate(eigv)

end subroutine diagonalize_matrix

!##################################################################################################################################

subroutine write_outputfiles( paramarray, eigenvalues )

implicit none
real,dimension(ndofs),intent(in) :: paramarray
real,dimension(syssize),intent(in) :: eigenvalues
integer :: bsfile,intfile,eigvfile,iter

open(newunit=bsfile,file='bs.out',status='replace',action='write')
open(newunit=intfile,file='int.out',status='replace',action='write')
open(newunit=eigvfile,file='eigenvalues.out',status='replace',action='write')

write(intfile,*) (paramarray(iter),iter=1,syssize)
write(bsfile,*) (paramarray(iter),iter=syssize+1,ndofs)
write(eigvfile,*) (eigenvalues(iter),iter=1,syssize)

close(bsfile)
close(intfile)
close(eigvfile)

end subroutine write_outputfiles

!##################################################################################################################################



!#################################################################################################################################

real function tor(val,stat)

implicit none

real,intent(in) :: val
integer,intent(in) :: stat

integer :: iter
real,dimension(7) :: coeff

if (stat == 0) then !GS
    coeff(1) = -0.000140462
    coeff(2) = -0.00111639
    coeff(3) = -0.000106193
    coeff(4) = 0.000351812
    coeff(5) = -0.00000883105
    coeff(6) = 0.0000387425
    coeff(7) = 0.00100871
elseif (stat == 1) then !ES
    coeff(1) = 0.00450247
    coeff(2) = -0.0314782
    coeff(3) = -0.00467622
    coeff(4) = 0.00375236
    coeff(5) = 0.00057716
    coeff(6) = -0.0000485252
    coeff(7) = 0.0213311
elseif (stat == 7) then !coupling
    coeff(1) = 0.000371217
    coeff(2) = 0.00271413
    coeff(3) = -0.000518811
    coeff(4) = 0.0000837664
    coeff(5) = 0.000812405
    coeff(6) = -0.000513085
    coeff(7) = -0.0483783
else
    write(*,*) 'second argument to tor must be 0, 1 or 7!'
    stop
end if

tor = coeff(7)
do iter=1,6
    tor=tor+coeff(iter)*cos(dble(iter)*val)
end do

end function tor

!#################################################################################################################################

real function bs(val,stat)

implicit none

real,intent(in) :: val
logical,intent(in) :: stat

real,dimension(4) :: coeff

if (stat) then !GS
    coeff(1) = 0.103532
    coeff(2) = 0.953379
    coeff(3) = 0.0
    coeff(4) = 0.0
else !ES
    coeff(1) = 0.173653
    coeff(2) = 1.05246
    coeff(3) = -0.265807
    coeff(4) = -0.0102905
end if

bs = coeff(1)*(exp(-coeff(2)*(val-coeff(3)))-1.0)*(exp(-coeff(2)*(val-coeff(3)))-1.0)+coeff(4)

end function bs

!#################################################################################################################################

real function intf(val,stat)

implicit none

real,intent(in) :: val
logical,intent(in) :: stat

real,dimension(4) :: coeff

if (stat) then !GS
    coeff(1) = 2.15672
    coeff(2) = 0.485731
    coeff(3) = 0.0
    coeff(4) = 0.0
else !ES
    coeff(1) = 0.100369
    coeff(2) = 3.07391
    coeff(3) = 0.138426
    coeff(4) = -0.0353169
end if

intf = coeff(1)*(exp(-coeff(2)*(val-coeff(3)))-1.0)*(exp(-coeff(2)*(val-coeff(3)))-1.0)+coeff(4)

end function intf

!#################################################################################################################################

END MODULE PES
