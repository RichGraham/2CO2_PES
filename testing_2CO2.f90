module testing_2CO2
  implicit none

    public :: testPerms
  private


  
contains


  

subroutine testPerms(rab)
  use GP_2CO2_variables
  use PES_2CO2_details
  use perm_array_2CO2!!! Uses the permutation group that was loaded at start of program
  implicit none
  double precision, allocatable:: rabPerm(:)
  integer perm_count
  double precision, dimension(:) :: rab
  double precision, allocatable:: EPerm(:)
  double precision :: toll = 1e-7, toll2=1e-8, error2, PES_2CO2
  integer :: i
  
  allocate( rabPerm(nDim))
  allocate( EPerm(nPerms))
  !return
  !print *,'*********Perm 1 *********************'
  EPerm(1) = PES_2CO2( rab)
  do i=2,nPerms
     EPerm(i) = EPerm(1)
  enddo
  
  do perm_count=2,nperms
   !  print *,'*********Perm ',perm_count,' *********************'
     !! Carry out the permutation on the interatomic distances
     do i=1,nDim
        rabPerm(i)  =  rab( perm( i, perm_count ) )
     enddo
     EPerm( perm_count) =  PES_2CO2( rabPerm)
     !print *,'Perm ',perm_count,EPerm(perm_count),EPerm(1)

     !! Catch large absolute or % differences
     if( abs(EPerm(1))<1e-17) then
        error2 = 0.0
     else
        error2  =  (EPerm(1) - EPerm(perm_count))  /  EPerm(1)
     endif

     !print *, EPerm(1), error2,abs(EPerm(1)-EPerm(perm_count))
     
     if( abs(EPerm(1)-EPerm(perm_count))>toll .and. abs( error2 )>toll2) then
        print *,'****************** PERMUTATION ERROR (Perm ',perm_count,')*************************'
        print *,' rab =',rab
        print *, 'E(Perm1)',EPerm(1), 'E(other perm)',EPerm(perm_count), 'Abs err:',abs(EPerm(1)-EPerm(perm_count)), &
             '  % err:', abs(error2)*100
        ERROR STOP
     endif
  enddo

  return;

end subroutine testPerms
  



!   initialize a random seed from the system clock at every run (fortran 95 code)
subroutine init_random_seed()
  INTEGER :: i, n, clock
  INTEGER, DIMENSION(:), ALLOCATABLE :: seed
  CALL RANDOM_SEED(size = n)
  ALLOCATE(seed(n))
  CALL SYSTEM_CLOCK(COUNT=clock)
  seed = clock + 37 * (/ (i - 1, i = 1, n) /)
  CALL RANDOM_SEED(PUT = seed)
    DEALLOCATE(seed)
  end subroutine init_random_seed

end module testing_2CO2
