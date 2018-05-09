module types
 implicit none
 integer,parameter   :: max_atom_number = 1000
 integer,parameter   :: maxnp = 10, maxcompounds = 1, maxdata = 1000
 integer,parameter   :: maxlinelength=maxnp*32
 type   CIFfile
  character(len=100) :: filename
  integer            :: n_atoms
  real               :: rv(1:3,1:3)
  real               :: vr(1:3,1:3)
  real               :: cell_0(1:6)
  real               :: atom_xcrystal(1:3,1:max_atom_number)
  real               :: atom_charge(1:max_atom_number)
  character(len=4)   :: atom_label(1:max_atom_number)
  character(len=2)   :: type_symbol(1:max_atom_number)
  real               :: obs_energy
  real               :: cal_energy
 end type
 type                          :: typ_ga
  character(len=maxlinelength) :: genotype
  real                         :: phenotype(1:maxnp)
  real                         :: fitness
 end type  
end module types
!
module mod_random
! module for pseudo random numbers
 implicit none
 private
 public init_random_seed, randint, r4_uniform
contains
 subroutine init_random_seed(seed)
  implicit none
  integer, intent(out) :: seed
  integer   day,hour,i4_huge,milli,minute,month,second,year
  parameter (i4_huge=2147483647)
  double precision temp
  character*(10) time
  character*(8) date
  call date_and_time (date,time)
  read (date,'(i4,i2,i2)')year,month,day
  read (time,'(i2,i2,i2,1x,i3)')hour,minute,second,milli
  temp=0.0D+00
  temp=temp+dble(month-1)/11.0D+00
  temp=temp+dble(day-1)/30.0D+00
  temp=temp+dble(hour)/23.0D+00
  temp=temp+dble(minute)/59.0D+00
  temp=temp+dble(second)/59.0D+00
  temp=temp+dble(milli)/999.0D+00
  temp=temp/6.0D+00
  doext: do
    if(temp<=0.0D+00 )then
       temp=temp+1.0D+00
       cycle doext
    else
       exit doext
    end if
  enddo doext
  doext2: do
    if (1.0D+00<temp) then
       temp=temp-1.0D+00
       cycle doext2
    else
       exit doext2
    end if
  end do doext2
  seed=int(dble(i4_huge)*temp)
  if(seed == 0)       seed = 1
  if(seed == i4_huge) seed = seed-1
  return
 end subroutine init_random_seed
 integer function randint(i,j,seed)
  real               ::  a,b
  integer,intent(in) ::  i,j,seed
  a = real(i)
  b = real(j)
  randint=int(r4_uniform(a,b+1.0,seed))
 end function randint
 real function r4_uniform(b1,b2,seed)
  implicit none
  real b1,b2
  integer i4_huge,k,seed
  parameter (i4_huge=2147483647)
  if(seed == 0) then
   write(*,'(b1)')'R4_UNIFORM - Fatal error!'
   write(*,'(b1)')'Input value of SEED = 0.'
   stop '[ERROR] Chiquitan chiquitintatan'
  end if
  k=seed/127773
  seed=16807*(seed-k*17773)-k*2836
  if(seed<0) then
    seed=seed+i4_huge
  endif
  r4_uniform=b1+(b2-b1)*real(dble(seed)* 4.656612875D-10)
  return
 end function r4_uniform
end module
!
module GeometricProperties
 implicit none
 private
 public cell,uncell,output_gulp,output_gulp_fit,Clen,Clen_trim
 contains
 PURE INTEGER FUNCTION Clen(s)      ! returns same result as LEN unless:
 CHARACTER(*),INTENT(IN) :: s       ! last non-blank char is null
 INTEGER :: i
 Clen = LEN(s)
 i = LEN_TRIM(s)
 IF (s(i:i) == CHAR(0)) Clen = i-1  ! len of C string
 END FUNCTION Clen
 PURE INTEGER FUNCTION Clen_trim(s) ! returns same result as LEN_TRIM unless:
 CHARACTER(*),INTENT(IN) :: s       ! last char non-blank is null, if true:
 INTEGER :: i                       ! then len of C string is returned, note:
                                    ! Ctrim is only user of this function
 i = LEN_TRIM(s) ; Clen_trim = i
 IF (s(i:i) == CHAR(0)) Clen_trim = Clen(s)   ! len of C string
 END FUNCTION Clen_trim
 SUBROUTINE cell(rv,vr,cell_0)
 implicit none
 real, intent(in)  :: cell_0(6)
 real, intent(out) :: rv(3,3),vr(3,3)
 real, parameter   :: pi = ACOS(-1.0)
 real :: alp,bet
 real :: cosa,cosb,cosg
 real :: gam,sing
 real :: DEGTORAD
 DEGTORAD=pi/180.0
 IF(cell_0(4) == 90.0) THEN
   cosa = 0.0
 ELSE
   ALP=cell_0(4)*degtorad
   COSA=cos(ALP)
 ENDIF
 IF(cell_0(5) == 90.0) THEN
   cosb = 0.0
 ELSE
   bet = cell_0(5)*degtorad
   cosb = cos(bet)
 ENDIF
 IF(cell_0(6) == 90.0) then
   sing = 1.0
   cosg = 0.0
 ELSE
   gam = cell_0(6)*degtorad
   sing = sin(gam)
   cosg = cos(gam)
 ENDIF
 rv(1,1) = cell_0(1)
 rv(1,2) = cell_0(2)*cosg
 rv(1,3) = cell_0(3)*cosb
 rv(2,1) = 0.0
 rv(2,2) = cell_0(2)*sing
 rv(2,3) = cell_0(3)*(cosa - cosb*cosg)/sing
 rv(3,1) = 0.0
 rv(3,2) = 0.0
 rv(3,3) = sqrt( cell_0(3)*cell_0(3) - rv(1,3)*rv(1,3) - rv(2,3)*rv(2,3))
 call inverse(rv,vr,3)
 RETURN
 END SUBROUTINE cell
!
 SUBROUTINE uncell(rv,cell_0)
  implicit none
  real,intent(out)   :: cell_0(6)
  real,intent(in)    :: rv(3,3)
  integer            :: i,j
  real               :: temp(6)
  REAL               :: radtodeg
  REAL, PARAMETER    :: pi=ACOS(-1.0)
  radtodeg=180.0/PI
  do i = 1,3
    temp(i) = 0.0
    do j = 1,3
      temp(i) = temp(i) + rv(j,i)**2
    enddo
    temp(i) = sqrt(temp(i))
  enddo
  cell_0(1) = abs(temp(1))
  cell_0(2) = abs(temp(2))
  cell_0(3) = abs(temp(3))
  do i = 1,3
    temp(3+i) = 0.0
  enddo
  do j = 1,3
    temp(4) = temp(4) + rv(j,2)*rv(j,3)
    temp(5) = temp(5) + rv(j,1)*rv(j,3)
    temp(6) = temp(6) + rv(j,1)*rv(j,2)
  enddo
  temp(4) = temp(4)/(temp(2)*temp(3))
  temp(5) = temp(5)/(temp(1)*temp(3))
  temp(6) = temp(6)/(temp(1)*temp(2))
  cell_0(4) = radtodeg*acos(temp(4))
  cell_0(5) = radtodeg*acos(temp(5))
  cell_0(6) = radtodeg*acos(temp(6))
  DO i=4,6
     if (abs(cell_0(i) - 90.0 ).lt.0.00001) cell_0(i) = 90.0
     if (abs(cell_0(i) - 120.0).lt.0.00001) cell_0(i) = 120.0
  ENDDO
  RETURN
 END SUBROUTINE uncell
 SUBROUTINE inverse(a,c,n)
 implicit none
 integer n
 real a(n,n), c(n,n)
 real L(n,n), U(n,n), b(n), d(n), x(n)
 real coeff
 integer i, j, k
 L=0.0
 U=0.0
 b=0.0
 do k=1, n-1
   do i=k+1,n
      coeff=a(i,k)/a(k,k)
      L(i,k) = coeff
      do j=k+1,n
         a(i,j) = a(i,j)-coeff*a(k,j)
      end do
   end do
 end do
 do i=1,n
  L(i,i) = 1.0
 end do
 do j=1,n
  do i=1,j
    U(i,j) = a(i,j)
  end do
 end do
 do k=1,n
  b(k)=1.0
  d(1) = b(1)
  do i=2,n
    d(i)=b(i)
    do j=1,i-1
      d(i) = d(i) - L(i,j)*d(j)
    end do
  end do
  x(n)=d(n)/U(n,n)
  do i = n-1,1,-1
    x(i) = d(i)
    do j=n,i+1,-1
      x(i)=x(i)-U(i,j)*x(j)
    end do
    x(i) = x(i)/u(i,i)
  end do
  do i=1,n
    c(i,k) = x(i)
  end do
  b(k)=0.0
 end do
 RETURN
 END SUBROUTINE inverse
 subroutine output_gulp(u,CIFFiles,GULPFilename)
  use types
  implicit none
  type(CIFfile),intent(inout)       :: CIFFiles
  character(len=100),intent(out)    :: GULPFilename
  integer,intent(in)                :: u
  integer                           :: i,k
  real               :: mmm,rrr
  integer            :: zzz
  character(len=2)   :: zlz
  character(len=4)   :: extension=".gin"
  GULPFilename=CIFFiles%filename(1:Clen_trim(CIFFiles%filename))//extension
  GULPFilename=adjustl(GULPfilename)
  open(u,file=GULPFilename)
  write(u,'(a)')'single conv mol'
  write(u,'(A)')'cell'
  write(u,'(6(f9.5,1x))') (CIFFiles%cell_0(i) , i=1,6)
  write(u,'(A)')'fractional'
  do i=1,CIFFiles%n_atoms
   write(u,'(a4,1x,3(f14.7,1x),1x,f14.7)')CIFFiles%atom_label(i),&
    (CIFFiles%atom_xcrystal(k,i),k=1,3),CIFFiles%atom_charge(i)
  end do
  write(u,'(a)')'supercell 2 1 1'
  write(u,'(A)')'library peros'
  !write(u,'(a,1x,a)')'output lammps ','test'
  !write(u,'(a,1x,a)')'output cif ','test'
  close(u)
 end subroutine output_gulp
 subroutine output_gulp_fit(u,n_files,CIFFiles)
  use types
  implicit none
  integer,intent(in)                :: u,n_files
  type(CIFfile),intent(in   )       :: CIFFiles(n_files)
  character(len=100)                :: GULPFilename="fit.gin"
  integer                           :: i,j,k
  real               :: mmm,rrr,obs_energy_min
  integer            :: zzz
  character(len=2)   :: zlz
  obs_energy_min=minval(CIFFiles%obs_energy)
  open(u,file=GULPFilename)
  write(u,'(a)')'fit conv molecule'
  do i=1,n_files
   write(u,'(A)')'cell'
   write(u,'(6(f9.5,1x))') (CIFFiles(i)%cell_0(k) , k=1,6)
   write(u,'(A)')'fractional'
   do j=1,CIFFiles(i)%n_atoms
    write(u,'(a4,1x,3(f14.7,1x),1x,f14.7)')CIFFiles(i)%atom_label(j),&
    (CIFFiles(i)%atom_xcrystal(k,j),k=1,3),CIFFiles(i)%atom_charge(j)
   end do
   write(u,'(a)')'supercell 2 1 1'
   write(u,'(a)')'observable'
   write(u,*)'energy eV'
   write(u,*)CIFFiles(i)%obs_energy - obs_energy_min, 100.0
   write(u,*)'end'
  end do
  !# Tell the program to fit the overall shift
  write(u,'(A)')'vary'
  write(u,'(A)')' shift'
  write(u,'(A)')'end'
  write(u,'(A)')'library peros'
  close(u)
 end subroutine output_gulp_fit
end module GeometricProperties
!
module get_structures
 implicit none
 private
 public GenerateCIFFileList, ReadListOfCIFFiles, ReadCIFFiles, WriteEnergies
 contains
  subroutine GenerateCIFFileList()
   implicit none
   character(len=100)  :: string=" "
   call system("if [ -f list ] ; then rm list ; touch list ; fi")
   write(6,'(a)') "if [ -f list ] ; then rm list ; touch list ; fi"
   string="ls struc/*.cif > list"
   call system(string)
   write(6,'(a)')string
  end subroutine GenerateCIFFileList
!
  subroutine ReadListOfCIFFiles(n_files)
   implicit none
   character(len=100)  :: line = " "
   integer             :: ierr = 0
   integer,intent(out) :: n_files
   n_files = 0
   open(111,file="list",iostat=ierr)
   read_list: do
    read(111,'(a)',iostat=ierr) line
    if(ierr/=0) exit read_list
    n_files=n_files+1
   end do read_list
   rewind(111)
   close(111)
   return
  end subroutine
!
  subroutine ReadCIFFiles(n_files,CIFFiles)
   use types
   use GeometricProperties
   implicit none
   real                              :: energy = 0.0
   integer                           :: i,j,k
   integer                           :: ierr = 0
   integer,intent(in)                :: n_files
   type(CIFfile),intent(inout)       :: CIFFiles(n_files)
   character(len=100)                :: line = " "
   character(len=100)                :: filename = " "
   character(len=20)                 :: spam
   character(len=80)                 :: string_stop_head= "_atom_site_charge"
   real                              :: infinite = 3.4028e38
   open(111,file="list",iostat=ierr)
   if(ierr/=0)stop
   do i=1,n_files
    read(111,'(a)')line
    read(line(1:39),'(a)') CIFFiles(i)%filename
    read(line(40:),*)      CIFFiles(i)%obs_energy
    write(6,'(a)')trim(CIFFiles(i)%filename)
    open(100,file=trim(CIFFiles(i)%filename),status='old',iostat=ierr)
    if(ierr/=0) stop 'CIFFile does not found'
    read_cif: do
     read(100,'(a)',iostat=ierr) line
     if(ierr/=0) exit read_cif
     if(line(1:14)=="_cell_length_a")then
      read(line,*)spam,CIFFiles(i)%cell_0(1)
      cycle read_cif
     end if
     if(line(1:14)=="_cell_length_b")then
      read(line,*)spam,CIFFiles(i)%cell_0(2)
      cycle read_cif
     end if
     if(line(1:14)=="_cell_length_c")then
      read(line,*)spam,CIFFiles(i)%cell_0(3)
      cycle read_cif
     end if
     if(line(1:17)=="_cell_angle_alpha")then
      read(line,*)spam,CIFFiles(i)%cell_0(4)
      cycle read_cif
     end if
     if(line(1:16)=="_cell_angle_beta")then
      read(line,*)spam,CIFFiles(i)%cell_0(5)
      cycle read_cif
     end if
     if(line(1:17)=="_cell_angle_gamma")then
      read(line,*)spam,CIFFiles(i)%cell_0(6)
      cycle read_cif
     end if
     if(line(1:)==string_stop_head) exit read_cif
    end do read_cif
    call cell(CIFFiles(i)%rv,CIFFiles(i)%vr,CIFFiles(i)%cell_0)
    CIFFiles(i)%n_atoms=0
    read_natoms: do
     read(100,'(a)',iostat=ierr) line
     if(ierr/=0) exit read_natoms
     CIFFiles(i)%n_atoms=CIFFiles(i)%n_atoms+1
    end do read_natoms
    rewind(100)
    write(6,'(80a)')('=',k=1,80)
    write(6,*)'Observable, energy:',CIFFiles(i)%obs_energy
    write(6,*)'Atoms:', CIFFiles(i)%n_atoms, CIFFiles(i)%filename
    write(6,'(3(f14.7,1x))') ( CIFFiles(i)%rv(1,j), j=1,3 )
    write(6,'(3(f14.7,1x))') ( CIFFiles(i)%rv(2,j), j=1,3 )
    write(6,'(3(f14.7,1x))') ( CIFFiles(i)%rv(3,j), j=1,3 )
    do
     read(100,'(a)') line
     if(line(1:)==string_stop_head) exit
    end do
    do j=1,CIFFiles(i)%n_atoms
     read(100,'(a)') line
     read(line,*) CIFFiles(i)%atom_label(j),CIFFiles(i)%type_symbol(j),&
                 (CIFFiles(i)%atom_xcrystal(k,j),k=1,3),CIFFiles(i)%atom_charge(j)
     write(6,'(a2,1x,3(f14.7,1x))')CIFFiles(i)%type_symbol(j),(CIFFiles(i)%atom_xcrystal(k,j),k=1,3)
    end do
    close(100)
    call output_gulp(444,CIFFiles(i),filename )
    write(6,*)'GULP file:',filename
    !line="~/GULP-4.2.0/Src/gulp-4.2.0 < "//filename(1:50)//" > tmp "
    line="~/bin/gulp < "//filename(1:50)//" > tmp "
    call system(line)
    line="grep 'Total lattice energy       =' tmp | grep 'eV' | awk '{print $5}' > c"
    call system(line)
    open(456,file="c")
    read(456,'(a)')line
    if(line(1:20)=="********************")then
     CIFFiles(i)%cal_energy=infinite
    else
     read(line,*) CIFFiles(i)%cal_energy
    end if
    close(456)
    line="rm c tmp "//filename//" "
    call system(line)
    write(6,*)'Calculated, energy', CIFFiles(i)%cal_energy
   end do
   close(111)
   call output_gulp_fit(444,n_files,CIFFiles)
   return
  end subroutine ReadCIFFiles
  subroutine WriteEnergies(n_files,CIFFiles,add)
   use types
   implicit none
   integer, intent(in)         :: n_files
   character(len=3),intent(in) :: add
   type(CIFfile),intent(in)    :: CIFFiles(n_files)
   integer                     :: i,u=123,ierr=0
   real                        :: obs_energy_min,cal_energy_min
   character(len=20)           :: filename = " "
   filename="fitness_"//add(1:3)//".txt"
   obs_energy_min=minval(CIFFiles%obs_energy)
   cal_energy_min=minval(CIFFiles%cal_energy)
   open(u,file=filename,iostat=ierr)
   if(ierr/=0) stop "fitness.txt can not be open"
   write(u,'(a)')"# struc.;  cell_size/A ;  Rela. Energy Obs / eV; DIFF ; Rela. Energy Cal. / eV ; Obs.Energy ; Cal.Energy"
   do i=1,n_files
    write(u,*)i,CIFFiles(i)%cell_0(1),CIFFiles(i)%obs_energy-obs_energy_min,&
                CIFFiles(i)%cal_energy-cal_energy_min,abs((CIFFiles(i)%obs_energy-obs_energy_min)-&
               (CIFFiles(i)%cal_energy-cal_energy_min)),&
                CIFFiles(i)%obs_energy,CIFFiles(i)%cal_energy
   end do
   close(u)
  end subroutine WriteEnergies
!
end module get_structures
!
module qsort_c_module
! Recursive Fortran 95 quicksort routine sorts real numbers into ascending numerical order
! Author: Juli Rew, SCD Consulting (juliana@ucar.edu), 9/03
! Based on algorithm from Cormen et al., Introduction to Algorithms, 1997 printing
! Made F conformant by Walt Brainerd
 implicit none
 public  :: QsortC
 private :: Partition
 contains
 recursive subroutine QsortC(A)
  real(16), intent(in out), dimension(:) :: A
  integer                            :: iq
  if(size(A) > 1) then
     call Partition(A, iq)
     call QsortC(A(:iq-1))
     call QsortC(A(iq:))
  endif
 end subroutine QsortC
 subroutine Partition(A, marker)
  real(16), intent(in out), dimension(:) :: A
  integer, intent(out)               :: marker
  integer                            :: i, j
  real(16)                           :: temp
  real(16)                           :: x
  x = A(1)
  i= 0
  j= size(A) + 1
  do
     j = j-1
     do
        if (A(j) <= x) exit
        j = j-1
     end do
     i = i+1
     do
        if (A(i) >= x) exit
        i = i+1
     end do
     if (i < j) then
        temp = A(i)
        A(i) = A(j)
        A(j) = temp
     elseif (i == j) then
        marker = i+1
        return
     else
        marker = i
        return
     endif
  end do
 end subroutine Partition
end module qsort_c_module
!
module fitter_globals
 use types
 use mod_random
 use get_structures
 implicit none
 integer                    :: npar,i,seed = 0
 integer,allocatable        :: np(:)
 integer                    :: err_apertura,ii,intervalos,j
 !integer,parameter          :: integration_points = 10000
 !real,parameter             :: precision_Newton = 1e-5
 !real                       :: tol = 0.001, tolfire = 0.25
 real,allocatable           :: param(:,:)
 !,concy(:),pi(:,:),iso(:,:)
 !integer,allocatable        :: npress(:)
 real,target                :: datas(2,maxcompounds,maxdata),f(maxdata)
 real,pointer               :: x(:),y(:),alldat(:,:)
 character(15),allocatable  :: ajuste(:,:)
 character(32*maxnp),allocatable :: string_IEEE(:)
 character(100)             :: line,string
 character(5)               :: inpt
 logical                    :: flag = .true., FlagFire = .false.,seed_flag=.true.
 logical                    :: physical_constrains = .false., range_flag =.true.
 logical                    :: refit_flag = .false., flag_shift=.false.
 real,parameter             :: R = 0.008314472 ! kJ / mol / K
 real                       :: T = 298.0
 !real                       :: inferior
 !real                       :: superior
 contains
  subroutine read_input()
  implicit none
  character(len=32)  :: chain
  allocate(np(1))
  npar = 0
  read_input_do: do
   read(5,'(A)',iostat=err_apertura)line
   if ( err_apertura /= 0 ) exit read_input_do
   if(line(1:1)=='#') cycle read_input_do
   if(line(1:7)=="shifted") then
    write(6,*)'[WARN] Shifting ( n_par -> n_par + 1 )'
    read(line,*)inpt,flag_shift
    if( flag_shift ) then
     npar=npar+1
    end if
    cycle read_input_do
   end if
   if(line(1:5)=='n_par')then
    read(line,*)inpt,np(1)
    npar=npar+np(1)
    allocate(ajuste(1:2,1:npar))
    allocate(param(1,0:npar-1))
    allocate(string_IEEE(1))
    np(1) = npar
    param(1,:) = 0.0
    do i=1,npar
     read(5,'(a)') line
     read(line( 1:10),'(a)') ajuste(1,i)
     read(line(12:),'(a)')   ajuste(2,i)
     write(6,'(a10,1x,a10)') ajuste(1,i),ajuste(2,i)
    end do
    cycle read_input_do
   end if
   if(line(1:5)=='ffit?')then
    read(line,*)inpt,flag
    cycle read_input_do
   end if
   if(line(1:5)=='refit') then
    write(6,*)'[WARN] Refitting parameters'
    read(line,*)inpt,refit_flag
    if(refit_flag)then
     do j=1,npar
      read(5,'(a32)')chain(1:32)
      write(6,'(a32)')chain(1:32)
      string_IEEE(1)(32*(j-1)+1:32*j)=chain(1:32)
     end do
     write(6,'(a)')string_IEEE(1)(1:32*npar)
    else
     string_IEEE(1) = ' '
    end if
   end if
   if(line(1:5)=='RSeed') then
    seed_flag=.false.
    read(line,*)inpt, seed
   end if
   if(line(1:22)=='physically_constrained') then
    physical_constrains=.true.
    write(6,'(a)') '[WARN] The fits are physically constrained'
   end if
   if(err_apertura/=0) exit read_input_do
  end do read_input_do
 end subroutine read_input
!
 subroutine ReadObservables(n_files,CIFFiles)
  implicit none
  integer                        :: i
  integer,intent(in)             :: n_files
  type(CIFfile),intent(in)       :: CIFFiles(n_files)
  do i=1,n_files
   datas(1,1,i)=real(i)
   datas(2,1,i)=CIFFiles(i)%obs_energy
   write(6,'(2f20.5)') datas(1,1,i),datas(2,1,i)
  end do
 end subroutine ReadObservables
end module fitter_globals
!
module mod_genetic
 use types
 use mod_random
 use fitter_globals
 use qsort_c_module
 use get_structures
 use GeometricProperties
 implicit none
 private
 public fit
 integer,parameter             :: ga_size         = 32 !2**10 ! numero de cromosomas
 real,parameter                :: ga_mutationrate = 0.3333 !2000/real(ga_size) ! ga_mutationrate=0.333
 real,parameter                :: ga_eliterate= 0.25, GA_DisasterRate = 0.0000001
 integer,parameter             :: ga_elitists = int( ga_size * ga_eliterate)                   
 type(typ_ga), pointer         :: parents(:)
 type(typ_ga), pointer         :: children(:)
 type(typ_ga), target          :: pop_alpha( ga_size )
 type(typ_ga), target          :: pop_beta( ga_size )
 contains
!
 type(typ_ga) function new_citizen(compound,seed,n_files,CIFFiles)
  implicit none
  integer                     :: i,j,k
  integer,intent(in)          :: n_files,compound,seed
  type(CIFFile),intent(inout) :: CIFFiles(n_files)
  real                        :: infinite = 3.4028e38
  new_citizen%fitness = infinite
  new_citizen%genotype = ' '
  !do while ( new_citizen%fitness == infinite )
   do i = 1,32*np(compound)
    new_citizen%genotype(i:i) = achar(randint(48,49,seed))
   end do
   do i = 1,np(compound)
    read(new_citizen%genotype(32*(i-1)+1:32*i),'(b32.32)') new_citizen%phenotype(i)
   end do
   new_citizen%fitness = fitness( new_citizen%phenotype,compound,n_files,CIFFiles)
  !end do
  return
 end function new_citizen
!
 subroutine UpdateCitizen( axolotl ,compound , n_files, CIFFiles)
  implicit none
  integer,intent(in)          :: compound,n_files
  type(CIFFile),intent(inout)    :: CIFFiles(n_files)
  integer                     :: i,GA_ELITISTS
  real                        :: infinite = 0.0
  type(typ_ga), intent(inout) :: axolotl
  do i = 1,np(compound)
   read(axolotl%genotype(32*(i-1)+1:32*i),'(b32.32)') axolotl%phenotype(i)
  end do
  axolotl%fitness = fitness( axolotl%phenotype,compound,n_files,CIFFiles)
  return
 end subroutine UpdateCitizen
 integer*4 function get_file_unit (lu_max)
  integer*4 lu_max,  lu, m, iostat
  logical   opened
  m = lu_max  ;  if (m < 1) m = 97
  do lu = m,1,-1
   inquire (unit=lu, opened=opened, iostat=iostat)
   if (iostat.ne.0) cycle
   if (.not.opened) exit
  end do
  get_file_unit = lu
  return
 end function get_file_unit
!
 real function Fitness(phenotype,compound,n_files,CIFFiles)
  use omp_lib
  implicit none
  real, intent(in)    :: phenotype(maxnp)
  integer,intent(in)  :: n_files
  type(CIFfile),intent(inout) :: CIFFIles(n_files)
  integer             :: i,compound,k = 0, u,ii,np_real
  real                :: a(0:np(compound)-1),xx,yy,penalty,obs_energy_min,cal_energy_min,obs_energy_max
  character(len=100)  :: funk,filename(n_files)
  character(len=200)  :: line
  logical             :: flagzero = .false.
  real                :: infinite = 3.4028e38
  call system("cp peros_input.lib peros.lib")
  penalty=0.0
  !if( flag_shift ) then 
  ! np_real = np(compound) - 1
  ! fitness = phenotype(np(compound))
  !else
  ! np_real = np(compound) 
   fitness = 0.0
  !end if
  refer: do i = 1,np(compound)
   write(line,*)"sed -i 's/",ajuste(1,i)(1:Clen_trim(ajuste(1,i))),"/",phenotype(i),"/g' peros.lib"
   call system(line)
   phys_constrains: if ( physical_constrains ) then
    funk=ajuste(2,i)
    select case(funk)
     case("A_buck")
      if (phenotype(i)<=1e-2.or.isnan(phenotype(i)).or.phenotype(i)>=1e15)then
       penalty=infinite
       exit refer
      end if
     case("A_lj")
      if (phenotype(i)<=1e-5.or.isnan(phenotype(i)).or.phenotype(i)>=1e10)then
       penalty=infinite
       exit refer
      end if
     case("rho_buck")
      if (phenotype(i)<1.0e-8.or.phenotype(i)>1.0.or.isnan(phenotype(i)))then
       penalty=infinite
       exit refer
      end if
     case("C_buck")
      if(phenotype(i)<0.0.or.phenotype(i)>1e7.or.isnan(phenotype(i)))then 
       penalty=infinite
       exit refer
      end if
     case("E_shift")
      if(phenotype(i)>=maxval(CIFFiles%obs_energy))then
       penalty=infinite
       exit refer
      end if
    end select
   end if phys_constrains
  end do refer
  calgulp: if ( penalty < 1.0 ) then
  !$omp parallel default(private) shared(n_files, CIFFiles, datas, filename, u,ii, line)
  ii=0
  scan_: do
   !$omp critical
   ii=ii+1
   i=ii
   if(i .le. n_files)then
    u=get_file_unit(444)
    call output_gulp(u,CIFFiles(i),filename(i))
    write(line,*)"~/bin/gulp < ",filename(i)(1:Clen_trim(filename(i)))," > ",&
     filename(i)(1:Clen_trim(filename(i))),".gout "
    call system(line)
    write(line,*)"grep 'Total lattice energy       =' ",filename(i)(1:Clen_trim(filename(i))),&
     ".gout | grep 'eV' | awk '{print $5}' > ",filename(i)(1:Clen_trim(filename(i))),".tmp "
    call system(line)
    write(line,*)"grep 'Total lattice energy       =' ",filename(i)(1:Clen_trim(filename(i))),&
     ".gout | grep 'eV' | awk '{print $5}' > ",filename(i)(1:Clen_trim(filename(i))),".tmp "
    call system(line)
    u=get_file_unit(444)
    open(u,file=filename(i)(1:Clen_trim(filename(i)))//".tmp")
    read(u,'(a)')line
    if(line(1:20)=="********************")then
     CIFFiles(i)%cal_energy=infinite
    else
     read(line,*) CIFFiles(i)%cal_energy
    end if
    close(u)
    end if
   !$omp end critical
   if(i.gt.n_files) exit
  end do scan_
  !$omp end parallel
  obs_energy_min=minval(CIFFiles%obs_energy)
  cal_energy_min=minval(CIFFiles%cal_energy)
  fitness = fitness + sum(abs( CIFFiles(1:n_files)%obs_energy-obs_energy_min -&
   (CIFFiles(1:n_files)%cal_energy-cal_energy_min))**2)/real(2*n_files)
  else
   fitness = fitness + penalty
  end if calgulp
  !call system("rm if [ $(echo '$(ls *.tmp  | wc -l) > 0' | bc -l) == 1 ] ; then rm -rf *.tmp  ; fi")
  !call system("rm if [ $(echo '$(ls *.gout | wc -l) > 0' | bc -l) == 1 ] ; then rm -rf *.gout ; fi")
  return
 end function Fitness

 subroutine WriteCitizen(k,kk,kkk,compound,lod,vgh)
  implicit none
  integer,intent(in)             :: k,kk,compound,lod,vgh,kkk
  integer                        :: i
  character(len=100)             :: fmt_
  character(len=32*np(compound)) :: wnowaste
  real                           :: wnowasteparam(1:32*np(compound)),wfitness
  do i=1,32*np(compound)
   wnowaste(i:i)=' '
  end do
  ! ...
  wnowaste = parents(k)%genotype
  do i=1,np(compound)
   wnowasteparam(i) = parents(k)%phenotype(i)
  end do
  wfitness = parents(k)%fitness
  write(6,'(i2,1x,i5,1x,a32,1x,e25.12,1x,f14.7,1x,a,1x,a,i8,a,i8,a,1x,a)')compound,kk,wnowaste(1:32),&
       wnowasteparam(1),wfitness,'[Fitness]','(',kkk,'/',vgh,')','[Similarity]' !,k
  do i=2,np(compound)
   if(lod>0.and.i==3)then
    write(6,'(9x,a32,1x,e25.12,10x,a,1x,i2,a)')wnowaste(1+32*(i-1):32*i),wnowasteparam(i),'Finishing:',lod,'/10'
   else
    write(6,'(9x,a32,1x,e25.12)')wnowaste(1+32*(i-1):32*i),wnowasteparam(i)
   end if
  end do
 end subroutine WriteCitizen

 subroutine piksrt(n,arr)
 ! Sort real array() of n elements
 implicit none
 integer :: n,j,i = 0
 REAL    :: a,arr(n)
 do j=2,n
   a = arr(j)
   do i=j-1,1,-1
     if (arr(i)<=a) goto 10
     arr(i+1)=arr(i)
   end do
   i=0
10  arr(i+1)=a
 end do
 return
 end subroutine piksrt

 subroutine SortByFitness()
  type(typ_ga)             ::  sorted(1:ga_size)
  integer                  ::  k,i
  real(16)                 ::  ftnss(1:ga_size)
  do k=1,ga_size
   ftnss(k)=dble(parents(k)%fitness)
   if(isnan(parents(k)%fitness)) ftnss(k) = 9999999999.d99
  end do
  call QsortC( ftnss )
  exter:do k=1,ga_size ! <- ordered
   inter:do i=1,ga_size
   if( dble(parents(i)%fitness) == ftnss(k))then
     sorted(k) = parents(i)
     cycle inter
   end if
   end do inter
  end do exter
  parents=sorted
  return
 end subroutine SortByFitness

 integer function Biodiversity( compound, animalito)
  implicit none
  integer,intent(in)             :: Compound
  type(typ_ga), intent(in)       :: animalito(1:ga_size)
  integer                        :: suma
  integer                        :: i,j,k,cont
  character(len=20)              :: mode = 'Normal'
  logical                        :: flag = .false.
  real                           :: error_d = 1e-3
  select case (mode)
   case('None')
    Biodiversity = 0
   case('Normal')
    suma=0
    Biodiversity = 0
    suma=0.5*ga_size*ga_size-ga_size
    do k =1,ga_size
     do j=k+1,ga_size 
      !suma = suma + 1
      if( animalito(k)%genotype(1:32*np(Compound)) == animalito(j)%genotype(1:32*np(Compound)) )then
       Biodiversity = Biodiversity + 1
      end if
     end do
    end do
   case('Superficial')
    Biodiversity = 0.0
    suma=0
    do k = 1, ga_size
     do j = k+1,ga_size
      cont=0
      suma=suma+1
      dbio: do i = 1,np(compound)
       bio: if( abs(animalito(k)%phenotype(i) - animalito(j)%phenotype(i)) <= error_d )then
        cont=cont+1
       end if bio
      end do dbio
      if( cont == np(compound) ) Biodiversity = Biodiversity + 1
     end do
    end do
  end select
  return
 end function Biodiversity

 subroutine Mutate( macrophage , compound )
  implicit none
  type(typ_ga), intent(inout) :: macrophage
  integer                     :: ipos,compound
  do i = 1,np(compound)
   ipos = randint(32*(i-1)+1,32*i,seed)
   macrophage%genotype(ipos:ipos) = achar(randint(48,49,seed))
  end do
  return
 end subroutine Mutate

 subroutine NuclearDisaster(Compound,n_files,CIFFiles)
  implicit none
  integer,intent(in)      ::  Compound,n_files
  type(CIFFile),intent(inout):: CIFFiles(n_files)
  integer                 :: k = 0, i, j
  real                    :: rrr
  do i = GA_ELITISTS + 1, GA_Size
   do j=1,32*np(compound)
    Children%genotype(j:j) = achar(randint(48,49,seed))
   end do
   call UpdateCitizen(Children(i),Compound,n_files,CIFFiles)
  end do
  return
 end subroutine NuclearDisaster

 subroutine Swap()
  if (associated(parents, target=pop_alpha)) then
      parents => pop_beta
      children => pop_alpha
  else
      parents => pop_alpha
      children => pop_beta
  end if
  return
 end subroutine Swap

 subroutine Elitism()
  children(:GA_ELITISTS) = parents(:GA_ELITISTS)
  return
 end subroutine

 subroutine Mate(compound,n_files,CIFFiles)
  integer             :: i, i1, i2, spos
  integer, intent(in)      :: compound,n_files
  type(CIFFile),intent(inout) :: CIFFiles(n_files)
  real                :: rrr
  call Elitism()
  do i = GA_ELITISTS + 1, ga_size
   ! Crossover:
   ! {{ eleccion random del primer 50% de la tabla
   call choose_randomly(i1,i2)
   ! }}
   ! {{ eleccion proporcionalmente a su fitness
   !call choose_propto_fitness(i1,i2)
   !write(6,*)i1,i2
   ! }}
   spos = randint(0, 32*np(compound), seed )
   children(i)%genotype = parents(i1)%genotype(:spos) // parents(i2)%genotype(spos+1:)
   ! Mutate and NuclearDisaster:
   rrr = r4_uniform(0.0,1.0,seed)
   if ( rrr < GA_MUTATIONRATE) then
    call Mutate(children(i),compound)
   else if ( rrr >= GA_MutationRate .and. rrr <= GA_MutationRate + GA_DisasterRate ) then
    call NuclearDisaster(Compound,n_files,CIFFiles)
   end if
  end do
  do i = 1, ga_size
   call UpdateCitizen(children(i),compound,n_files,CIFFiles)
  end do
  return
 end subroutine Mate

 subroutine choose_randomly(j1,j2)
  implicit none
  integer,intent(out) :: j1,j2
  j1  = randint(1, int(ga_size/2),seed)
  j2  = randint(1, int(ga_size/2),seed)
  do while ( j1 == j2 )
   j2 = randint(1, int(ga_size/2),seed)
  end do
  return
 end subroutine choose_randomly

 subroutine choose_propto_fitness(j1,j2)
  implicit none
  integer,intent(out) :: j1,j2
  integer             :: i
  real                :: ftnss(ga_size),prop(0:ga_size)=0.0,rrr1,rrr2
  real                :: infinity = HUGE(2147483647)
  rrr1 = 0.0
  do i = 1, ga_size
   ftnss(i) = 1.0/parents(i)%fitness
   if ( isnan( parents(i)%fitness ) ) ftnss(i) = 0.0
   if ( parents(i)%fitness > infinity ) ftnss(i) = 0.0
   prop(i) = ftnss(i)
   if( ftnss(i) >= infinity ) then
     rrr1 = rrr1 + infinity
   else
     rrr1 = rrr1 + ftnss(i)
   end if
  end do
  prop = prop / rrr1
  ! select 1:
   rrr1 = r4_uniform(0.0,1.0,seed)
   slct1: do i=1,ga_size
    if(rrr1<=prop(i-1).and.rrr1>prop(i))then
     j1 = i
    end if
   end do slct1
   ! select 2:
   rrr2 = r4_uniform(0.0,1.0,seed)
   do while ( rrr1 == rrr2 )
    rrr2 = r4_uniform(0.0,1.0,seed)
   end do
   slct2: do i=1,ga_size
    if(rrr2<=prop(i-1).and.rrr2>prop(i))then
     j2 = i
    end if
   end do slct2
  return
 end subroutine choose_propto_fitness

 subroutine Fit(Compound,Seed,n_files,CIFFiles)
  implicit none
  integer,intent(in)       :: Compound, Seed,n_files
  type(CIFFile),intent(inout) :: CIFFiles(n_files)
  integer,parameter  :: maxstep = 1000, minstep = 10
  integer            :: kk, ii, i, k,vgh
  real               :: diff = 0.0, fit0 = 0.0
  integer            :: eps 
  kk = 0
  ii = 0
  pop_alpha = [(new_citizen(compound,seed,n_files,CIFFiles), i = 1,ga_size)]
  parents =>  pop_alpha
  children => pop_beta
  if ( refit_flag ) then
   do i = 1, 2
    parents(i)%genotype=string_IEEE(compound)
    children(i)%genotype=string_IEEE(compound)
    call UpdateCitizen(parents(i),compound,n_files,CIFFiles)
    call UpdateCitizen(children(i),compound,n_files,CIFFiles)
   end do
   call Mate(compound,n_files,CIFFiles)
   call Swap()
   call SortByFitness()
  end if
  call WriteCitizen(1,ii,1,compound,kk, int(0.5*ga_size*ga_size-ga_size) )
  converge: do while ( .true. )
   ii=ii+1
   call SortByFitness()
   eps = Biodiversity( compound, children)
   diff = eps
   call WriteCitizen(1,ii,eps,compound,kk, int(0.5*ga_size*ga_size-ga_size) )
   !fire: if ( FlagFire ) then
   ! if ( ii >= minstep .and. parents(1)%fitness <= TolFire ) exit converge
   !else
    if( ii>=minstep .and. parents(1)%fitness <= 0.1 .and. abs(parents(1)%fitness-fit0) <= 1e-4)then
    !if( abs(diff - eps) <= 1e-2 .and. ii >= minstep .and. &
    ! parents(1)%fitness - fit0 == 0 ) then
     kk = kk + 1
    else
     kk = 0
    end if
    if ( ii >= maxstep .or. kk >= 10 ) exit converge
   !end if fire
   call Mate(compound,n_files,CIFFiles)
   call Swap()
   fit0 = parents(1)%fitness
  end do converge
  do i = 0, np( compound )-1
   param( compound,i ) = children(1)%phenotype(i+1)
  end do
  write(6,*)'#',(param(compound,i),i=0,np(compound )-1)
  write(6,*)'#','Fitness:',fit0,'Similarity:',eps,'Rseed',seed
  return
 end subroutine fit
end module mod_genetic
!
program fitter
 use iso_fortran_env
 use types
 use mod_random
 use get_structures
 use fitter_globals
 use mod_genetic
 implicit none
 type(CIFfile),allocatable       :: CIFFiles(:)
 integer                         :: n_files=0
 print '(4a)', 'This file was compiled by ', &
       compiler_version(), ' using the options ', &
       compiler_options()
 call read_input()
 if (seed_flag) then
  seed = 73709
  call init_random_seed(seed)
  seed_flag=.false.
 end if
 write(6,'("Random Seed:",1x,i10)') seed
 !call GenerateCIFFileList()
 call ReadListOfCIFFiles(n_files)
 allocate( CIFFiles(1:n_files) )
 call ReadCIFFiles(n_files,CIFFiles)
 call WriteEnergies(n_files,CIFFiles,"ini")
 call ReadObservables(n_files,CIFFiles)
 if (flag) call fit(1,seed,n_files,CIFFiles)
 call WriteEnergies(n_files,CIFFiles,"end") 
end program fitter
