
#define M_PI_2 (1.57079632679489661923_rk)

module Global
    integer, parameter  ::  rk   = 8
    
    ! this should stay a parameter so that the
    ! compiler is able to optimize loops involving dims
    integer, parameter  ::  dims = 2
    
    ! we need PI/2 for the potential
    real(rk), parameter :: pi_2 = M_PI_2
end module
