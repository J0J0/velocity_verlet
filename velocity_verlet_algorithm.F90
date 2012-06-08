
#ifndef M_PI_2
#define M_PI_2 (1.57079632679489661923_rk)
#endif

program velocity_verlet_algorithm
    
    integer, parameter          ::  rk      = 8,    &
                                    foffset = 20,   &
                                    dims    = 2
    real(rk), parameter         ::  tmax    = 2.0
    integer                     ::  i, j, n
    real(rk), dimension(:,:)    ::  r, v, a
    allocatable                 ::  r, v, a
    real(rk)                    ::  t, dt
    real(rk)                    ::  d, dcut
    real(rk), dimension(dims)   ::  rij, force, lbounds, ubounds
    integer                     ::  stat, sysclock, seed_size
    integer, allocatable        ::  seed(:)
    character(4)                ::  arg
    character(32)               ::  tmp_file_str
    
    call get_command_argument(1, arg, status=stat)
    if(stat /= 0) stop 1
    read(arg, *) n
    allocate(r(dims,n), v(dims,n), a(dims,n))
    
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
    
    t  = 0.0
    dt = 0.01
    
    do i = 1, n
        write(foffset+i, *) t, r(:,i)
    end do
    
    do
        t = t + dt
        
        do i = 1, n
            r(:,i) = r(:,i) + v(:,i)*dt + 0.5*a(:,i)*dt*dt
            force = 0
            do j = 1, n
                if(i == j) cycle
                rij  = r(:,j) - r(:,i)
                d    = sqrt(sum(rij * rij))
                dcut = min(d,M_PI_2)
                force = force + sin(2.0*dcut) * rij/d
            end do
            v(:,i) = v(:,i) + 0.5*(a(:,i)+force)
            a(:,i) = force
            
            do j = 1, dims
                if( r(j,i) <= lbounds(j) ) then
                    if( r(j,i) < lbounds(j) ) &
                        r(j,i) = lbounds(j)
                    v(j,i) = -v(j,i)
                elseif( r(j,i) >= ubounds(j) ) then
                    if( r(j,i) > ubounds(j) ) &
                        r(j,i) = ubounds(j)
                    v(j,i) = -v(j,i)
                end if
            end do
            
            write(foffset+i, *) t, r(:,i)
        end do
        
        if( t > tmax ) exit
    end do
    
    
end program


