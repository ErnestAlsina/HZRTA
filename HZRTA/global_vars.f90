MODULE global_vars
!----------------------------------------------------------------------
IMPLICIT NONE
!------------------------------------------------------------------------
INTEGER :: nwl, nz
DOUBLE PRECISION :: nu0
DOUBLE PRECISION :: Aul, Hu
DOUBLE PRECISION :: Bint, them, chim
DOUBLE PRECISION :: gu, gl, nular,atwe
INTEGER :: ju2,jl2
INTEGER :: df
DOUBLE PRECISION, DIMENSION(:), allocatable :: Hup2
DOUBLE PRECISION, DIMENSION(:), allocatable :: nu, ht, dnd, ad, fr3,epsprm
DOUBLE PRECISION, DIMENSION(:), allocatable :: uj
DOUBLE PRECISION, DIMENSION(:), allocatable :: dnul !Zeeman-shifted change in freq
DOUBLE PRECISION, DIMENSION(:,:),allocatable :: Hup
DOUBLE PRECISION, DIMENSION(:,:,:,:), allocatable :: prof_vo, prof_fa
DOUBLE PRECISION, DIMENSION(:,:,:,:,:),allocatable :: gprof_re, gprof_im
!----------------------
DOUBLE PRECISION, DIMENSION(:,:,:,:,:),allocatable :: eta_re, eta_im, rho_re, rho_im
DOUBLE PRECISION, DIMENSION(:,:,:,:,:), allocatable :: gen_prof
DOUBLE PRECISION, DIMENSION(:,:,:,:,:,:,:),allocatable :: cfac
DOUBLE PRECISION, DIMENSION(:), allocatable :: tau !keep for now
DOUBLE PRECISION, DIMENSION(:,:,:), allocatable :: Cij
DOUBLE PRECISION, DIMENSION(:,:,:,:,:), allocatable :: jrad
DOUBLE PRECISION, DIMENSION(:,:,:,:,:), allocatable :: jrout
DOUBLE PRECISION, DIMENSION(:,:,:), allocatable :: stokh 

DOUBLE PRECISION, DIMENSION(0:6,-6:6,-6:6) :: DR,DI
DOUBLE PRECISION, DIMENSION(:,:,:,:,:), allocatable :: TTR,TTI !Global variable for all the T^K_Q tensors
DOUBLE PRECISION, DIMENSION(0:3,0:3) :: iden
END MODULE global_vars
