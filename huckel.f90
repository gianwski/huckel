! Author: Gianluca Regni
! Copyright (c) 2024 Gianluca Regni
! License: GPL-3.0 (https://www.gnu.org/licenses/)
! Credits: Please cite the author if you use or reference this program.
!
! Description: 
! The program constructs and diagonalizes the Hamiltonian of a 1D polyene (Hückel Matrix),
! whether linear or cyclic, based on an input keyword. Additionally, it calculates the
! HOMO, LUMO, and the band gap.
! The program considers two types of dimerization as defined in the input:
! - Bond alternation (e.g., different bond lengths)
! - Atom alternation (e.g., different atom types)
!
! Before running the program, ensure the LAPACK software package (LAPACK Users' Guide (Third) (1999);
! DOI: 10.1137/1.9780898719604) is installed in your machine.
! It is available at https://www.netlib.org/lapack/.
! 
! Input file example:
! 20            ! number of atoms (integer)
! linear        ! set if the system is linear or cyclic (character)
! .TRUE.        ! dimerization: Bond alternation (logical)
! .FALSE.       ! dimerization: Atom alternation (logical)
! 1             ! number of atom in the unit cell (integer)
! 1             ! number of electron per unit cell (integer)
! 0             ! alpha/alpha1 alpha2 ... (real)
! -1.0 -0.8     ! beta/beta+ beta- (real)
! 4             ! number of eigenvector to print (integer)
! 1 2 3 4    ! eigenvector you want to print (integer)
!
! Usage:
! gfortran huckel.f90 -o huckel -llapack
! ./huckel

  program main

  implicit none

! Definition of variables

  integer          :: i, n, j, k, nevect, ntype, index, nelec, electot
  double precision :: beta_plus, beta_minus, band_gap
  character(len=6) :: structure
  logical          :: bond_dim, atom_dim

  double precision, dimension(:, :), allocatable :: H, Eigenvectors
  double precision, dimension(:), allocatable    :: Eigenvalues, eval_norm, alpha
  integer, dimension(:), allocatable             :: vector

! Open input and output files

  open(1,file="huckel.in",status="old")
  open(2,file="evalu.out")
  open(3,file="evect.out")
  open(4,file="huckel.out")

! Read input variables

  read(1,*) n
  read(1,*) structure

  if(structure /= "linear" .and. structure /= "cyclic" ) then
     write(6,*) "error: structure keyword is not set correctly"
     stop
  endif

  read(1,*) bond_dim
  read(1,*) atom_dim
  read(1,*) ntype

  if(atom_dim .and. ntype < 2) then
     write(6,*) "error: atom dimerization is set to true, but there is less than two atom per unit cell."
     stop
  endif
  read(1,*) nelec

  allocate(alpha(ntype))

  if(atom_dim) then
     read(1,*) (alpha(i), i=1,ntype)
  else
     read(1,*) alpha(1)
  endif

  if(bond_dim) then
     read(1,*) beta_plus, beta_minus
  else
     read(1,*) beta_plus
  endif

  if(structure == 'cyclic' .and. bond_dim .and. mod(n,2) /= 0) then
     write(6,*) "error: number of atom should be even"
     stop
  endif

  read(1,*) nevect

  if(nevect .gt. n) then 
     write(6,*) "error: number of eigenvector should be less or equal than number of atoms"
     stop
  endif

  allocate(vector(nevect))

  read(1,*) (vector(i), i=1, nevect)

! Print input variables

  write(6,*) " "
  write(6,"(A)") "--------------------- Input Variables ---------------------"
  write(6,*) " "
  write(6,"(A,I3)") "number of atoms:   ", n
  write(6,"(A,A)")  "structure:         ", structure
  write(6,"(A,L1)") "bond dimerization: ", bond_dim
  write(6,"(A,L1)") "atom dimerization: ", atom_dim
  write(6,"(A,I2)") "number of atom in the unit cell:     ", ntype
  write(6,"(A,I2)") "number of electron in the unit cell: ", nelec
  write(6,"(A,90(F6.3,1X))")  "alpha:          ", alpha(:)
  if(bond_dim) then
     write(6,"(A,2(F6.3,1X))") "beta (+ and -): ", beta_plus, beta_minus
  else
     write(6,"(A,F6.3)") "beta: ", beta_plus
  endif
  write(6,"(A,I4)") "number of eigenvector to print: ", nevect
  write(6,"(A, 90(I2, 1X))")  "index of eigenvector to print: ",  (vector(i), i=1, nevect)
  write(6, *) " "

! Write to the summary output
  
  write(4,"(A)") "--------------------- Input Variables ---------------------"
  write(4,*) " "
  write(4,"(A,I3)") "number of atoms:   ", n
  write(4,"(A,A)")  "structure:         ", structure
  write(4,"(A,L1)") "bond dimerization: ", bond_dim
  write(4,"(A,L1)") "atom dimerization: ", atom_dim
  write(4,"(A,I2)") "number of atom in the unit cell:     ", ntype
  write(4,"(A,I2)") "number of electron in the unit cell: ", nelec
  write(4,"(A,90(F6.3,1X))")  "alpha:          ", alpha(:)
  if(bond_dim) then
     write(4,"(A,2(F6.3,1X))") "beta (+ and -): ", beta_plus, beta_minus
  else
     write(4,"(A,F6.3)") "beta: ", beta_plus
  endif
  write(4,"(A,I4)") "number of eigenvector to print: ", nevect
  write(4,"(A, 90(I2, 1X))")  "index of eigenvector to print: ",  (vector(i), i=1, nevect)
  write(4, *) " "

! Allocate memory

  allocate(H(n,n), Eigenvectors(n,n), Eigenvalues(n), eval_norm(n))

! Initialize Hamiltonian matrix

  H = 0.d0

! Construction of the Hamiltonian

  do i = 1, n
     index = mod(i-1, ntype) + 1
     H(i, i) = alpha(index)
  enddo

  if(bond_dim) then
     do i = 1, n-1
        if(mod(i,2) == 0) then
           H(i, i+1) = beta_minus
           H(i+1, i) = beta_minus
        else
           H(i, i+1) = beta_plus
           H(i+1, i) = beta_plus
        endif
     end do
   else
     do i = 1, n-1
        H(i, i+1) = beta_plus
        H(i+1, i) = beta_plus
     enddo
   endif   

! Addition of the matrix elements to the corners if cyclic

  if(structure == 'cyclic') then
     if(bond_dim) then
        H(1, n) = beta_minus
        H(n, 1) = beta_minus
     else
        H(1, n) = beta_plus
        H(n, 1) = beta_plus 
     endif 
  endif

! Print Huckel Matrix and write to the summary output

  write(6,*) " "
  write(4,*) " "
  if(n < 21) then
     write(6,"(A)") "---------------------- Hückel Matrix ----------------------"
     write(4,"(A)") "---------------------- Hückel Matrix ----------------------"
     write(6,*) " "
     write(4,*) " "
     write(6,"(5X,30(I2,4X))") (k, k=1,n)
     do i=1, n
        write(6,"(I3,1X,30(F5.2,1X))") i, H(i, :)
     end do
     write(4,"(5X,30(I2,4X))") (k, k=1,n)
     do i=1, n
        write(4,"(I3,1X,30(F5.2,1X))") i, H(i, :)
     end do
  else
     write(6,"(A)") "Hückel Matrix is too big to be printed here"
     write(4,"(A)") "Hückel Matrix is too big to be printed here"
  endif

! Diagonalization of Huckel matrix

  call diagonalization(H,Eigenvalues,Eigenvectors,N)

! Calculation of norm for eigenvalues for plot purpose

  do i = 1, n
     eval_norm(i) = (real(i)-1)/(real(n)-1) 
  end do

! Write eigenvalues to the output

  write(2,*) "# Eigenvalues"
  do i=1, n
        write(2,"(I3,1X,30(F9.6,1X))") i, eval_norm(i), Eigenvalues(i)
  end do

! Write summary to the output

  write(4,*) " "
  write(4,"(A)") "----------------------- Eigenvalues -----------------------"
  write(4,*) " "
  do i=1, n
        write(4,"(I3,1X,30(F11.6,1X))") i, eval_norm(i), Eigenvalues(i)
  end do

! Calculation of total π-electron

  electot = n/ntype * nelec
  
! Calculation of band gap

  if (mod(electot,2) == 0) then
      band_gap = Eigenvalues(electot/2+1)-Eigenvalues(electot/2)
  else
      band_gap = Eigenvalues(int(electot/2)+2)-Eigenvalues(int(electot/2)+1)
  endif

! Print HOMO, LUMO and the band gap

  write(6,*) " "
  write(6,"(A)") "------------------------ Energies -------------------------"
  write(6,*) " "
  write(6,"(A,I4)") "electrons: ", electot
  if (mod(electot,2) == 0) then
      write(6,"(A, 1X, F9.6)") "HOMO:     ", Eigenvalues(electot/2)
      write(6,"(A, 1X, F9.6)") "LUMO:     ", Eigenvalues(electot/2+1)
  else
      write(6,"(A, 1X, F9.6)") "HOMO:     ", Eigenvalues(int(electot/2)+1)
      write(6,"(A, 1X, F9.6)") "LUMO:     ", Eigenvalues(int(electot/2)+2)
  endif
  write(6,"(A, 1X, F9.6)") "band gap: ", band_gap
  write(6, *) " "

! Write summary to the output

  write(4,*) " "
  write(4,"(A)") "------------------------ Energies -------------------------"
  write(4,*) " "
  write(4,"(A,I4)") "electrons: ", electot
  if (mod(electot,2) == 0) then
      write(4,"(A, 1X, F9.6)") "HOMO:     ", Eigenvalues(electot/2)
      write(4,"(A, 1X, F9.6)") "LUMO:     ", Eigenvalues(electot/2+1)
  else
      write(4,"(A, 1X, F9.6)") "HOMO:     ", Eigenvalues(int(electot/2)+1)
      write(4,"(A, 1X, F9.6)") "LUMO:     ", Eigenvalues(int(electot/2)+2)
  endif
  write(4,"(A, 1X, F9.6)") "band gap: ", band_gap
  write(4, *) " "

! Write eigenvector defined to the output

  write(3,"(A,2X,A,4X,30(I2,8X))") "#", "norm", (vector(i), i=1, nevect)
  do i = 1, n
     write(3,"(I3,1X,30(F9.6,1X))") i, eval_norm(i), (Eigenvectors(i, vector(j)), j=1, nevect)
  end do

! Write to the summary output

  write(4,*) " "
  write(4,"(A)") "---------------------- Eigenvectors -----------------------"
  write(4,*) " "
  write(4,"(7X,A,7X,30(I2,8X))") "norm", (vector(i), i=1, nevect)
  do i = 1, n
     write(4,"(I3,1X,30(F9.6,1X))") i, eval_norm(i), (Eigenvectors(i, vector(j)), j=1, nevect)
  end do

  close(1)
  close(2)
  close(3)

  end program main

  subroutine diagonalization(matrix, eigenvalues, eigenvectors, N)

! The subroutine performs the diagonalization of a real symmetric matrix

  implicit none
                                                   
  double precision, dimension(N,N), intent(in)    :: matrix
  double precision, dimension(N,N), intent(inout) :: eigenvectors
  double precision, dimension(N),   intent(inout) :: eigenvalues
  double precision, dimension(:), allocatable     :: WORK
  integer, intent(in)                             :: N
  integer                                         :: INFO, LWORK

  allocate(WORK(1))

  eigenvectors = matrix

  call DSYEV('V', 'U', N, eigenvectors, N, eigenvalues, WORK, -1, INFO)

! "V" states the subroutine computes eigenvectors (and eigenvalues)
! "U" states the input matrix is memorized in the upper triangular part 
! N defines the dimension of the square matrix
! matrix_temp is the matrix that contains input matrix
! N defines the dimension of eigenvalues array
! eigenvalues is the array where the subroutine will store the eigenvalues
! WORK is an array used for temporary data
! -1 indicates that the subroutine had not to computes eigenvalues or/and eigenvector,
! but it have to return back the optimal dimension of WORK as first element of the array
! INFO is an output parameter to communicate the output state of the subroutine

  if (INFO /= 0) stop ! INFO equal 0 means that the calculation is successful

! Optimal lenght of WORK array is stored in the first element of WORK

  LWORK = WORK(1)
  deallocate(WORK)
  allocate(WORK(LWORK))

  call DSYEV('V', 'U', N, eigenvectors, N, eigenvalues, WORK, LWORK, INFO)

! LWORK replace -1 and now indicates the lenght of WORK array

  if (INFO < 0) write(6,*) "diagonalization failure: wrong argument"
  if (INFO > 0) write(6,*) "diagonalization failure: convergence not reached"
  if (INFO /= 0) stop

  deallocate(WORK)

  endsubroutine diagonalization
