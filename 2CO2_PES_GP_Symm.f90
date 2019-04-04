!!CO2-Ar PES
!!Distances in Angstrom and energies in Hartrees
module PES_2CO2_details
  double precision :: gpRMax = 9.0 
  double precision :: gpRMin = 1.5  
  double precision :: lCO = 1.1632
  double precision :: AngToBohr =  1.8897259885789
  interface PES_2CO2_GP 
     function PES_2CO2_GP(xStar) 
       implicit none 
       double precision:: PES_2CO2_GP
       double precision, dimension(:) ::  xStar
     end function PES_2CO2_GP
  end interface PES_2CO2_GP
end module PES_2CO2_details


module GP_2CO2_variables
  double precision, allocatable :: alpha (:), lScale(:), xTraining(:,:), xTrainingPerm(:,:,:)
  double precision expVar,NuggVar, gpEmax
  integer :: nDim=9
  integer :: nTraining=146
  integer :: nPerms=8
end module GP_2CO2_variables


module perm_array_2CO2
  integer, allocatable:: perm(:,:)
end module perm_array_2CO2


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!Gaussian Process Code!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!pom


  
subroutine load_GP_2CO2_Data
  use GP_2CO2_variables
  use PES_2CO2_details
  use perm_array_2CO2
  implicit none
  
  double precision, allocatable::  xStar(:)
  integer i,j,k
  double precision :: dum, expVar1, expVar2,expVar3
  character (len=90) :: filename
  !integer, allocatable:: perm(:,:)
  CHARACTER(len=255) :: homedir,codedir

  allocate (alpha(nTraining), lScale(nDim), xTraining(nDim,nTraining),xTrainingPerm(nDim,nPerms,nTraining),  &
       xStar(nDim), perm(nDim,nPerms))


  CALL getenv("HOME", homedir)
  codedir=TRIM(homedir) // '/source/2CO2_PES'
  !!call chdir(codedir)

  !====Load hyperparameters====
  write (filename, '( "/TrainingData/HyperParams_Symm", I3.3, ".dat" )' )  nTraining
  filename = trim(codedir) // filename
  !write (filename, '( "TrainingData/myTest.dat" )' )  
  !open (unit = 7, file = filename)
  !! Older symmetric way
  !Only need to read some as others are tied.
  !!read (7,*) lScale(1),lScale(2), lScale(5), expVar3, expVar2, expVar1,NuggVar, gpEmax
  
  open (unit = 7, file = filename)
  do i=1,nDim
     read (7,*) dum
     j=int(dum)
     !print *,i,j
     read (7,*) lScale(j)
  end do
  
  read (7,*) expVar
  read (7,*) NuggVar
  read (7,*) gpEMax
 
  !!print *,"HyperParams, lScale=",lScale(1), lScale(2),lScale(3),lScale(4), lScale(5), &
    !!   lScale(6),lScale(7),lScale(8),lScale(9)
  
  !!print *,"HyperParams",expVar,NuggVar, gpEmax
  close(7)
  
  !====Load alpha coefficients====
  write (filename, '( "/TrainingData/alpha_Symm", I3.3, ".dat" )' )  nTraining
  filename = trim(codedir) // filename
  open (unit = 7, file = filename)
  do i=1,nTraining
     read (7,*) alpha(i)
     !!print *,"alpha ",i, alpha(i)
  end do
  close(7)


  !====Load training data x values ====
  write (filename, '( "/TrainingData/xTraining", I3.3, ".dat" )' )  nTraining
  filename = trim(codedir) // filename
  open (unit = 7, file = filename)
    
  do i=1,nTraining
     read (7,*) xTraining(1,i), xTraining(2,i), xTraining(3,i), xTraining(4,i), xTraining(5,i),&
          xTraining(6,i), xTraining(7,i), xTraining(8,i), xTraining(9,i)
     !print *,xTraining(1,i), xTraining(2,i), xTraining(3,i), xTraining(4,i), xTraining(5,i), &
     !     xTraining(6,i), xTraining(7,i), xTraining(8,i), xTraining(9,i)
  end do
  close(7)

  write (filename, '( "/2CO2.sym" )' )
  filename = trim(codedir) // filename

  open (unit = 7, file = filename,STATUS='OLD')
  read(7,*) perm

  !do i=1,nPerms
   !  print *, perm(1,i), perm(2,i), perm(3,i), perm(4,i), perm(5,i), perm(6,i), perm(7,i),perm(8,i), perm(9,i)
  !end do
  
  !! Permute the training vectors
  
  do i=1,nDim
     do j=1,nPerms
        do k=1,nTraining
           xTrainingPerm(i,j,k)=xTraining(perm(i,j),k)
        end do
     end do
  end do

end subroutine load_GP_2CO2_Data
  
function PES_2CO2_GP(xStar)
  use GP_2CO2_variables
  implicit none
  double precision, dimension(:) :: xStar
  double precision:: PES_2CO2_GP
  integer i,j,k
  double precision kSqExpAllPerms, kSqExpJthPerm, kKernTotal


  !! Non-symmetric way
  !kKernTotal=0
  !do i=1,nTraining
  !   kSqExpAllPerms=1.0;
  !   do k=1,nDim
  !      kSqExpAllPerms=kSqExpAllPerms * &
  !           ( exp( - (xStar(k)-xTraining(k,i))**2 /2.0/lScale(k)**2) )
  !   end do
  !   kKernTotal = kKernTotal + alpha(i)*kSqExpAllPerms
  !end do

  !Symmetric way
  kKernTotal=0
  do i=1,nTraining
     kSqExpAllPerms=0
     do j=1,nPerms
        kSqExpJthPerm=1        
        do k=1,nDim
           kSqExpJthPerm  =  kSqExpJthPerm * &
                ( exp( - (xStar(k)-xTrainingPerm(k,j,i))**2 /2.0/lScale(k)**2) )
        end do !Dimensions (k)
        kSqExpAllPerms = kSqExpAllPerms + kSqExpJthPerm
     end do !Permuations (
     kKernTotal = kKernTotal + alpha(i) * kSqExpAllPerms
  end do !Training points (i)
  
  PES_2CO2_GP=kKernTotal * expVar
end function PES_2CO2_GP



function PES_2CO2( rab )
  !! Takes in rab in Angstrom
  use PES_2CO2_details
  use GP_2CO2_variables
  implicit none
  double precision rab(9), xStar(9), asymp_2CO2
  double precision  PES_2CO2, sum
  double precision repFactor
  integer i

  repFactor=1.0
  
  if( minval(rab) > gpRMax ) then !!Use asymptotic function
     PES_2CO2 = asymp_2CO2(rab, lCO)
     
  else if (minval(rab)< gpRMin/repFactor ) then !! Use repulsive approximation function
     sum = 0.0
     do i=1,nDim
        sum = sum + 1.0/rab(i)**12
     enddo
     PES_2CO2=gpEmax* sum *gpRMin **12/(1.0*nDim)
     
  else !! Use the Guassian Process function
     xStar(:) = 1/rab(:)
     PES_2CO2 = PES_2CO2_GP( xStar)
  end if

  !PES_2CO2=gpEmax/3.0* (1.0/rab(1)**12+1.0/rab(2)**12+1.0/rab(3)**12) *gpRMin **12

  
end function PES_2CO2

double precision function asymp_2CO2(rab,rco2)
! Work out the asymptotic energy for CO2-CO2
! rab are the interatomic distances OO, OC, OO', CO, CC, CO', O'O, O'C', O'O'
! rco2 is the CO distances in the molecules
! Energy is in Atomic Units (Hartree)
! Distance is in Angstrom (Not atomic units)
implicit none
double precision rab(9),rco2,emult
double precision qc,qo,qcp,qop,uc,uo,ucp,uop,tcp,top
double precision cozz,cczz,cozx,cczx,coxz,ccxz,coxx,ccxx
double precision oozz,oczz,oozx,oczx,ooxz,ocxz,ooxx,ocxx
double precision asymp
! The energy consists of nibe pair contributions, each containing
! seven terms (six electrostatic and one dispersion).
! The function emult is used to get the seven terms.
! This function calls emult nine times, and adds the results.
asymp=0d0
! The arguments of emult are the relevant coordinates and parameters.
! Coordinates are rab, rab', ra'b, ra'b', raa', rbb' where rab is the distance
! between the atoms of interest, a' is bonded to a, b' is bonded to b.
! Parameters are charge of a, dipole of a (in the direction away from
! a'), charge of b, dipole of b (away from b'),
! and four ab dispersion energy coefficients, ||, |-, -| and --
! where | means parallel to the bond and - means perpendicular.
!
! Parameters
qc=0.8670d0 ! charges
qo=-0.4335d0
uc=0d0 ! dipoles; C doesn't have one
uo=0.1347d0 ! away from C
cczz=2.04d0
cozz=3.79d0 ! dispersion, z means parallel
oczz=cozz
oozz=7.09d0
cczx=1.42d0
ccxz=cczx
cozx=2.04d0 ! x means perpendicular
ocxz=cozx
coxz=2.69d0
oczx=coxz
oozx=3.87d0
ooxz=oozx
ccxx=1.03d0
coxx=1.49d0
ocxx=coxx
ooxx=2.16d0
! OO (so a=O, b=O, a'=C, b'=C)
asymp=asymp+emult(rab(1),rab(2),rab(4),rab(5),rco2,rco2,qo,uo,qo,uo,oozz,oozx,ooxz,ooxx)
! OC (so a=O, b=C, a'=C, b'=O)
asymp=asymp+emult(rab(2),rab(1),rab(5),rab(4),rco2,rco2,qo,uo,qc,uc,oczz,oczx,ocxz,ocxx)
! OO' (so a=O, b=O', a'=C, b'=C)
asymp=asymp+emult(rab(3),rab(2),rab(6),rab(5),rco2,rco2,qo,uo,qo,uo,oozz,oozx,ooxz,ooxx)
! CO (so a=C, b=O, a'=O, b'=C)
asymp=asymp+emult(rab(4),rab(5),rab(1),rab(2),rco2,rco2,qc,uc,qo,uo,cozz,cozx,coxz,coxx)
! CC (so a=C, b=C, a'=O, b'=O)
asymp=asymp+emult(rab(5),rab(4),rab(2),rab(1),rco2,rco2,qc,uc,qc,uc,cczz,cczx,ccxz,ccxx)
! CO' (so a=C, b=O', a'=O, b'=C)
asymp=asymp+emult(rab(6),rab(5),rab(3),rab(2),rco2,rco2,qc,uc,qo,uo,cozz,cozx,coxz,coxx)
! O'O (so a=O', b=O, a'=C, b'=C)
asymp=asymp+emult(rab(7),rab(8),rab(4),rab(5),rco2,rco2,qo,uo,qo,uo,oozz,oozx,ooxz,ooxx)
! O'C (so a=O', b=C, a'=C, b'=O)
asymp=asymp+emult(rab(8),rab(7),rab(5),rab(4),rco2,rco2,qo,uo,qc,uc,oczz,oczx,ocxz,ocxx)
! O'O' (so a=O', b=O', a'=C, b'=C)
asymp=asymp+emult(rab(9),rab(8),rab(6),rab(5),rco2,rco2,qo,uo,qo,uo,oozz,oozx,ooxz,ooxx)
asymp_2CO2=asymp
end
!
double precision function emult(rabang,rabpang,rapbang,rapbpang,raapang,rbbpang,qa,ua,qb,ub,d1,d2,d3,d4)
implicit none
double precision, parameter :: angstrom=1.88973d0
double precision rabang,rabpang,rapbang,rapbpang,raapang,rbbpang,qa,ua,qb,ub,d1,d2,d3,d4
double precision rab,rabp,rapb,rapbp,raap,rbbp
double precision eadd, costa, costb, cosphi, dc1, dc2, dc3, dc4
! Convert distances to Bohr
rab=rabang*angstrom
rabp=rabpang*angstrom
rapb=rapbang*angstrom
rapbp=rapbpang*angstrom
raap=raapang*angstrom
rbbp=rbbpang*angstrom
emult=0d0
! Calculate the multipolar energy (Coulomb+dispersion)
!write(6,*)'emult called with distances ',rab,rabp,rapb,rapbp,raap,rbbp
!write(6,*)'and parameters ',qa,ua,qb,ub,d1,d2,d3,d4
! Charges
eadd=qa*qb/rab
!write(6,*)'Charge-charge ',eadd
emult=emult+eadd
! Dipole of A with charge of B, for which we need the angle ta between
! the ap->a axis and the a->b vector
costa=(rapb**2-raap**2-rab**2)/(2d0*raap*rab)
!write(6,*)'cos(ta)=',costa
eadd=ua*qb*costa/rab**2
!write(6,*)'Dipole-charge ',eadd
emult=emult+eadd
! Charge of A with dipole of B, for which we need the angle tb between
! the bp->b axis and the b->a vector
costb=(rabp**2-rbbp**2-rab**2)/(2d0*rbbp*rab)
!write(6,*)'cos(tb)=',costb
eadd=qa*ub*costb/rab**2
!write(6,*)'Charge-dipole ',eadd
emult=emult+eadd
! Dipole of A with dipole of B, for which we need the angle phi between
! the ap->a axis and the bp->b axis (NB, not usual definition of phi)
cosphi=(rabp**2+rapb**2-rab**2-rapbp**2)/(2d0*raap*rbbp)
!write(6,*)'cos(phi)=',cosphi
eadd=ua*ub*(cosphi+3d0*costa*costb)/rab**3
!write(6,*)'Dipole-dipole ',eadd
emult=emult+eadd
!write(6,*)'Elec ',emult
! Dispersion par-par
dc1=d1*(cosphi+3*costa*costb)**2
! Dispersion par-perp
dc2=d2*((1+3*costa**2)-(cosphi+3*costa*costb)**2)
! Dispersion perp-par
dc3=d3*((1+3*costb**2)-(cosphi+3*costa*costb)**2)
! Dispersion perp-perp
dc4=d4*((4-3*costa**2-3*costb**2)+(cosphi+3*costa*costb)**2)
!write(6,*)'DCs ',dc1,dc2,dc3,dc4
eadd=-(dc1+dc2+dc3+dc4)/rab**6
!write(6,*)'Dispersion ',eadd
emult=emult+eadd
!write(6,*)'Total ',emult
end
