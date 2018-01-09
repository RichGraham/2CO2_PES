!!CO2-Ar PES
!!Distances in Angstrom and energies in Hartrees
module PES_details
  double precision :: gpRMax = 9.0  
  double precision :: gpRMin = 1.5  
  double precision :: lCO = 1.1632
  double precision :: AngToBohr =  1.8897259885789
  interface PES_GP 
     function PES_GP(xStar) 
       implicit none 
       double precision:: PES_GP
       double precision, dimension(:) ::  xStar
     end function PES_GP
  end interface PES_GP
end module PES_details


module GP_variables
  double precision, allocatable :: alpha (:), lScale(:), xTraining(:,:), xTrainingPerm(:,:,:)
  double precision expVar,NuggVar, gpEmax
  integer :: nDim=9
  integer :: nTraining=146
  integer :: nPerms=8
end module GP_variables


! Test program
use PES_details
implicit none
 
integer k,i, choice

double precision rab(9),e, PES, xStar(9)

xStar(1:9)=(/0.20971571,  0.22605512,  0.22877215,  0.20158676,  0.2298365 , &
     0.2503203 ,  0.18510657,  0.2182286 ,  0.25385964/)
!xStar(1:9)=(/ 0.3018619 ,  0.33131717,  0.31699997,  0.30609471,  0.34708315, &
!     0.34030303,  0.27652502,  0.3131383 ,  0.31583372 /)

!xStar(1:9)=(/ 0.16918005,  0.20963065,  0.27441895,  0.15582496,  0.18869334, &
!        0.23773219,  0.14123042,  0.16638506,  0.20074532 /)

!xStar(1:9)=(/ 0.37918623,  0.43618246,  0.39957494,  0.32684715,  0.42808256, &
!     0.48455317,  0.26286132,  0.34577986,  0.44897025 /)

!xStar(1:9)=(/0.3791862300000000129962530, 0.4361824599999999940713735, &
!     0.3995749399999999895705116, 0.3268471499999999752006374, 0.4280825600000000008549250,&
!     0.4845531699999999775130277, 0.2628613200000000094114228, 0.3457798599999999944465401, &
!     0.4489702500000000151558766/)

!xStar(1:9)=(/ 0.24927635,  0.29087296,  0.3123213 ,  0.277394  ,  0.30196993, &
!        0.29306223,  0.28166435,  0.27929086,  0.25204695 /)

call load_GP_Data
call fixedAngleSlice

!e=PES_GP( xStar)
!e=PES( rab)
!write(6,*)e

end
!

subroutine fixedAngleSlice()
  use PES_details
  implicit none
  double precision rab(9)
  integer i, itot
  double precision  r, beta1,beta2, alpha2, e, e_GP, asymp, PES, PI

  PI=4.D0*DATAN(1.D0)
    
  itot=500

  !T shape
  !beta1 = acos(0.0)
  !beta2 = acos(1.0)
  !alpha2 = 0.0

  !I shape
  !beta1 = acos(1.0)
  !beta2 = acos(1.0)
  !alpha2 = 0.0
  
  !Para shape
  !beta1 = acos(0.0)
  !beta2 = acos(0.0)
  !alpha2 = 0.0
  
  !X shape
  beta1 = acos(0.0)
  beta2 = acos(0.0)
  alpha2 = 0.5 * PI

  
  open (unit=15, file="PES_Out.dat ", status='replace')
  
  do i=0, itot

     ! specify centre-to-centre separation
     r = (  1.5 + 8.5*i/(1.0*itot) ) 

     call computeDistances(r,alpha2,beta1,beta2,rab)
     
     
     e=PES( rab)
     !e_GP = PES_GP( xStar)
     write(15,*) r , e 
     
  enddo
  write(6,*)'Written to file: PES_Out.dat '
  close(15)

end subroutine fixedAngleSlice
  
  
  
subroutine computeDistances(r,alpha2,beta1, beta2, rDist)
  use PES_details
  implicit none
  double precision  rDist(9), r, beta1, beta2, alpha2
  double precision O1(3),O2(3),O3(3),O4(3),C1(3),C2(3)
  double precision sqLength

  O1(:)=0
  O2(:)=0
  O3(:)=0
  O4(:)=0

  C1(:)=0
  C2(:)=0
  
  !CO2 1
  !C is already at the origin
  O1(1)=lCO * sin( beta1)
  O1(3)=lCO *cos(beta1)
  
  O2(1)=-lCO * sin( beta1)
  O2(3)=-lCO * cos(beta1)
  
 
  !CO2 2
  C2(3)= r
  
  O3(1) =  lCO * sin(beta2) * cos(alpha2)
  O3(2) =  lCO * sin(beta2) * sin(alpha2)
  O3(3) =  r  +  lCO * cos(beta2)
  
  O4(1) =  -lCO * sin(beta2) * cos(alpha2)
  O4(2) =  -lCO * sin(beta2) * sin(alpha2)
  O4(3) =  r  -  lCO * cos(beta2)
  
  
  rDist(1)= SQRT( sqLength(O1,O3) )
  rDist(2)= SQRT( sqLength(O1,C2) )
  rDist(3)= SQRT( sqLength(O1,O4) )


  rDist(4)= SQRT( sqLength(C1,O3) )
  rDist(5)= SQRT( sqLength(C1,C2) )
  rDist(6)= SQRT( sqLength(C1,O4) )

  rDist(7)= SQRT( sqLength(O2,O3) )
  rDist(8)= SQRT( sqLength(O2,C2) )
  rDist(9)= SQRT( sqLength(O2,O4) )

  
   
end subroutine computeDistances


double precision function sqLength( u,v)
  !returns the square of the distance between two vectors
  implicit none
  double precision sum
  double precision u(3), v(3)
  integer i

  sum=0.0
  do i = 1,3 
    sum = sum  +  (  v(i)-u(i)  )**2
  enddo
  
  sqLength= sum
end function sqLength



double precision function asymp(rab)
  use PES_details
  implicit none
  double precision rab(3)
  double precision c1,c2,c3, cAr, o1Ar, o2Ar

  c1= (  rab(1)**2 + lCO**2 - rab(2)**2)/2.0/rab(1)/lCO
  
  c2= (  rab(2)**2 + lCO**2 - rab(1)**2)/2.0/rab(2)/lCO
  c3= (  rab(3)**2 + lCO**2 - rab(1)**2)/2.0/rab(3)/lCO


  cAr  = - ( 4.64 * (1 + 3*c1**2)  +  3.30 * (5 - 3*c1**2) ) / (rab(1)*AngToBohr)**6
  o1Ar = - ( 8.69 * (1 + 3*c2**2)  +  4.76 * (5 - 3*c2**2) ) / (rab(2)*AngToBohr)**6
  o2Ar = - ( 8.69 * (1 + 3*c3**2)  +  4.76 * (5 - 3*c3**2) ) / (rab(3)*AngToBohr)**6
  
  asymp= cAr + o1Ar + o2Ar
end


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!Gaussian Process Code!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!pom


  
subroutine load_GP_Data
  use GP_variables
  use PES_details
  implicit none
  
  double precision, allocatable::  xStar(:)
  integer i,j,k
  double precision :: dum, expVar1, expVar2,expVar3
  character (len=90) :: filename
  integer, allocatable:: perm(:,:)

  allocate (alpha(nTraining), lScale(nDim), xTraining(nDim,nTraining),xTrainingPerm(nDim,nPerms,nTraining), xStar(nDim), &
       perm(nDim,nPerms))

  !====Load hyperparameters====
  write (filename, '( "TrainingData/HyperParams_Symm", I3.3, ".dat" )' )  nTraining
  !write (filename, '( "TrainingData/myTest.dat" )' )  
  !open (unit = 7, file = filename)
  !! Older symmetric way
  !Only need to read some as others are tied.
  !!read (7,*) lScale(1),lScale(2), lScale(5), expVar3, expVar2, expVar1,NuggVar, gpEmax
  
  open (unit = 7, file = filename)
  do i=1,nDim
     read (7,*) dum
     j=int(dum)
     print *,i,j
     read (7,*) lScale(j)
  end do
  
  read (7,*) expVar
  read (7,*) NuggVar
  read (7,*) gpEMax
 
  print *,"HyperParams, lScale=",lScale(1), lScale(2),lScale(3),lScale(4), lScale(5), &
       lScale(6),lScale(7),lScale(8),lScale(9)
  
  print *,"HyperParams",expVar,NuggVar, gpEmax
  close(7)
  
  !====Load alpha coefficients====
  write (filename, '( "TrainingData/alpha_Symm", I3.3, ".dat" )' )  nTraining
  open (unit = 7, file = filename)
  do i=1,nTraining
     read (7,*) alpha(i)
     !!print *,"alpha ",i, alpha(i)
  end do
  close(7)


  !====Load training data x values ====
  write (filename, '( "TrainingData/xTraining", I3.3, ".dat" )' )  nTraining
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

end subroutine load_GP_Data
  
function PES_GP(xStar)
  use GP_variables
  implicit none
  double precision, dimension(:) :: xStar
  double precision:: PES_GP
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
  
  PES_GP=kKernTotal * expVar
end function PES_GP



function PES( rab )
  !! Takes in rab in Angstrom
  use PES_details
  use GP_variables
  implicit none
  double precision rab(9), xStar(9), asymp
  double precision  PES
  double precision repFactor

  repFactor=1.0
  
  if( rab(1) > gpRMax  .AND.  rab(2) > gpRMax .AND.  rab(3) > gpRMax &
       ) then !!Use asymptotic function
     PES = 0
     
  else if (rab(1) < gpRMin/repFactor  .OR.  rab(2) < gpRMin/repFactor  .OR.  rab(3) < gpRMin/repFactor &
       ) then !! Use repulsive approximation function
     PES=gpEmax* (1.0/rab(1)**12+1.0/rab(2)**12+1.0/rab(3)**12) *gpRMin **12
     
  else !! Use the Guassian Process function
     xStar(:) = 1/rab(:)
     PES = PES_GP( xStar)
  end if

  !PES=gpEmax/3.0* (1.0/rab(1)**12+1.0/rab(2)**12+1.0/rab(3)**12) *gpRMin **12

  
end function PES
