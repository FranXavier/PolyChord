!> This allows for a simple C interface... 


module InitSampler

    ! ~~~~~~~ Loaded Modules ~~~~~~~
    use ini_module,               only: read_params,initialise_program
    use params_module,            only: add_parameter,param_type
    use priors_module
    use settings_module,          only: program_settings,initialise_settings
    use random_module,            only: initialise_random
    use nested_sampling_module,   only: NestedSampling
    use utils_module,             only: STR_LENGTH
    use abort_module,             only: halt_program
#ifdef MPI
    use mpi_module,               only: initialise_mpi, finalise_mpi
    use mpi,                      only: MPI_COMM_WORLD
#endif


    contains

subroutine DoSamplingFromC(Lfunc, nDims, inifile, context)

    ! ~~~~~~~ Local Variable Declaration ~~~~~~~
    implicit none

    !Lfunc is the likelihood to be evaluated by PolyChord
    !inifile is the configuration file
    !context is a pointer to a c struct
    
    !Input parameters for initialisation
    character(STR_LENGTH) :: inifile
    integer nDims, context
        
    double precision, dimension(5) :: output_info

    type(program_settings)    :: settings  ! The program settings 
    type(prior), dimension(:),allocatable     :: priors    ! The details of the priors

    character(len=STR_LENGTH)                 :: input_file     ! input file
    type(param_type),dimension(:),allocatable :: params         ! Parameter array
    type(param_type),dimension(:),allocatable :: derived_params ! Derived parameter array

    ! Temporary variables for initialising loglikelihoods
    double precision :: loglike

    !Had to change interface to likelihood to include the size of the arrays theta and phi, otherwise it seg faulted in C, alternatives welcome.
    interface
        function Lfunc(nDims, theta, nDerived, phi)
            integer,          intent(in)                 :: nDims
            integer,          intent(in)                 :: nDerived
            double precision, intent(in), dimension(nDims) ::theta
            double precision, intent(out),  dimension(nDerived) :: phi
            double precision :: Lfunc
        end function
    end interface

#ifdef MPI
    call initialise_mpi
#endif
    call initialise_random()
    settings%nDims = nDims
    settings%nDerived = nDims
    settings%precision_criterion = 0.5
    call read_params(trim(inifile),settings,params,derived_params)
    
    call initialise_program(settings,priors,params,derived_params)

#ifdef MPI
    output_info = NestedSampling(Lfunc,priors,settings,MPI_COMM_WORLD) 
#else
    output_info = NestedSampling(Lfunc,priors,settings,0) 
#endif
#ifdef MPI
    call finalise_mpi
#endif

end subroutine

end module InitSampler
