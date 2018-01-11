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


double precision function asymp_2CO2(rab)
  use PES_2CO2_details
  implicit none
  double precision rab(3)
  double precision c1,c2,c3, cAr, o1Ar, o2Ar

  c1= (  rab(1)**2 + lCO**2 - rab(2)**2)/2.0/rab(1)/lCO
  
  c2= (  rab(2)**2 + lCO**2 - rab(1)**2)/2.0/rab(2)/lCO
  c3= (  rab(3)**2 + lCO**2 - rab(1)**2)/2.0/rab(3)/lCO


  cAr  = - ( 4.64 * (1 + 3*c1**2)  +  3.30 * (5 - 3*c1**2) ) / (rab(1)*AngToBohr)**6
  o1Ar = - ( 8.69 * (1 + 3*c2**2)  +  4.76 * (5 - 3*c2**2) ) / (rab(2)*AngToBohr)**6
  o2Ar = - ( 8.69 * (1 + 3*c3**2)  +  4.76 * (5 - 3*c3**2) ) / (rab(3)*AngToBohr)**6
  
  asymp_2CO2= cAr + o1Ar + o2Ar
end


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!Gaussian Process Code!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!pom


  
subroutine load_GP_2CO2_Data
  use GP_2CO2_variables
  use PES_2CO2_details
  implicit none
  
  double precision, allocatable::  xStar(:)
  integer i,j,k
  double precision :: dum, expVar1, expVar2,expVar3
  character (len=90) :: filename
  integer, allocatable:: perm(:,:)
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

  write (filename, '( "2CO2.sym" )' )
  open (unit = 7, file = filename)
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
  double precision  PES_2CO2
  double precision repFactor

  repFactor=1.0
  
  if( rab(1) > gpRMax  .AND.  rab(2) > gpRMax .AND.  rab(3) > gpRMax &
       ) then !!Use asymptotic function
     PES_2CO2 = 0
     
  else if (rab(1) < gpRMin/repFactor  .OR.  rab(2) < gpRMin/repFactor  .OR.  rab(3) < gpRMin/repFactor &
       ) then !! Use repulsive approximation function
     PES_2CO2=gpEmax* (1.0/rab(1)**12+1.0/rab(2)**12+1.0/rab(3)**12) *gpRMin **12
     
  else !! Use the Guassian Process function
     xStar(:) = 1/rab(:)
     PES_2CO2 = PES_2CO2_GP( xStar)
  end if

  !PES_2CO2=gpEmax/3.0* (1.0/rab(1)**12+1.0/rab(2)**12+1.0/rab(3)**12) *gpRMin **12

  
end function PES_2CO2
