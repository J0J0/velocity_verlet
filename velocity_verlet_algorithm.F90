
#ifndef M_PI_2
#define M_PI_2 (1.57079632679489661923_rk)
#endif

program velocity_verlet_algorithm
    
    integer, parameter          ::  rk      = 8,    &
                                    foffset = 20,   &
                                    dims    = 2
    real(rk), parameter         ::  tmax    = 10.0
    
    integer                     ::  i, j, n
    real(rk)                    ::  t, dt
    
    real(rk), dimension(:,:)    ::  r, v, a, a_old
    allocatable                 ::  r, v, a, a_old
    
    real(rk)                    ::  d, dcut
    real(rk), dimension(dims)   ::  rij, force_ij, lbounds, ubounds

    integer                     ::  stat, sysclock, seed_size
    integer, allocatable        ::  seed(:)
    character(4)                ::  arg
    character(32)               ::  tmp_file_str
    
    call get_command_argument(1, arg, status=stat)
    if(stat /= 0) stop 1
    read(arg, *) n
    allocate(r(dims,n), v(dims,n), a(dims,n), a_old(dims,n))
    
    do i = 1, n
        write(tmp_file_str, "(a,i0.4,a)") "output/particle_", i, ".txt"
        open(foffset+i, file=trim(tmp_file_str), status='replace')
    end do
    
    call system_clock(count=sysclock)
    call random_seed(size=seed_size)
    allocate(seed(seed_size))
    seed = (/ (sysclock+i, i=1,seed_size) /)
    call random_seed(put=seed)
    deallocate(seed)
    
    lbounds = 0.0
    ubounds = 1.0

    call random_number(r)
    v = 0.0
    a = 0.0
    do i = 1, n
        do j = i+1, n
            rij  = r(:,j) - r(:,i)
            d    = sqrt(sum(rij * rij))
            dcut = min(d,M_PI_2)
            force_ij = sin(2.0*dcut) * rij/d
            
            a(:,i) = a(:,i) + force_ij
            a(:,j) = a(:,j) - force_ij
        end do
    end do
    
    t  = 0.0
    dt = 0.01
    
    do i = 1, n
        write(foffset+i, *) t, r(:,i)
    end do
    
    do
        t = t + dt
        if( t > tmax ) exit
        
        do i = 1, n
            r(:,i) = r(:,i) + v(:,i)*dt + 0.5*a(:,i)*dt*dt
            
            do j = 1, dims
                if( r(j,i) <= lbounds(j) ) then
                    r(j,i) = 2.0*lbounds(j) - r(j,i)
                    v(j,i) = -v(j,i)
                elseif( r(j,i) >= ubounds(j) ) then
                    r(j,i) = 2.0*ubounds(j) - r(j,i)
                    v(j,i) = -v(j,i)
                end if
            end do

            a_old(:,i) = a(:,i)
            a(:,i) = 0
            
            write(foffset+i, *) t, r(:,i)
        end do
        
        do i = 1, n
            do j = i+1, n
                rij  = r(:,j) - r(:,i)
                d    = sqrt(sum(rij * rij))
                dcut = min(d,M_PI_2)
                force_ij = sin(2.0*dcut) * rij/d
                
                a(:,i) = a(:,i) + force_ij
                a(:,j) = a(:,j) - force_ij
            end do
        end do
        
        do i = 1, n
            v(:,i) = v(:,i) + 0.5*(a_old(:,i)+a(:,i))*dt
        end do
    end do
    
    
end program


