
#ifndef M_PI_2
#define M_PI_2 (1.57079632679489661923_rk)
#endif

program velocity_verlet_algorithm
    
    integer, parameter          ::  rk      = 8,    &
                                    foffset = 20,   &
                                    dims    = 2
    real(rk), parameter         ::  tmax    = 0.10
    
    integer                     ::  i, j, n
    integer                     ::  loop_count
    real(rk)                    ::  t, dt, write_t, write_dt 
    
    real(rk), dimension(:,:)    ::  r, v, a, a_old
    allocatable                 ::  r, v, a, a_old
    
    real(rk)                    ::  d, dcut
    real(rk), dimension(dims)   ::  rij, force_ij, lbounds, ubounds
    real(rk)                    ::  kin, pot, energy, energy_old
    real(rK)                    ::  energy_init, energy_err_acum, energy_err_max

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
    
    open(foffset, file="output/energy.txt", status='replace')
    
    call system_clock(count=sysclock)
    call random_seed(size=seed_size)
    allocate(seed(seed_size))
    seed = (/ (sysclock+i, i=1,seed_size) /)
    call random_seed(put=seed)
    deallocate(seed)
    
    lbounds = 0.0
    ubounds = 1.0

    call random_number(r)
    call random_number(v)
    !v = 0.0
    a = 0.0
    
    pot = 0
    kin = 0
    do i = 1, n
        do j = i+1, n
            rij  = r(:,j) - r(:,i)
            d    = sqrt(sum(rij * rij))
            dcut = min(d,M_PI_2)
            pot  = pot - sin(dcut)*sin(dcut)
            force_ij = sin(2.0*dcut) * rij/d
            
            a(:,i) = a(:,i) + force_ij
            a(:,j) = a(:,j) - force_ij
        end do
        kin = kin + 0.5*sum(v(:,i)*v(:,i))
    end do
    energy = pot + kin
    energy_init = energy
    
    t  = 0.0
    dt = 1e-9
    write_t  = 0
    write_dt = 0.005
    
    do i = 1, n
        write(foffset+i, *) t, r(:,i)
    end do
    write(foffset, *) t, energy, 0.0, 0.0
    
    energy_err_acum = 0
    energy_err_max  = 0
    loop_count = 0
    do
        t = t + dt
        if( t > tmax ) exit
        loop_count = loop_count + 1
        
        pot = 0
        kin = 0
        
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
            
            if( t > write_t ) &
                write(foffset+i, *) t, r(:,i)
        end do
        
        do i = 1, n
            do j = i+1, n
                rij  = r(:,j) - r(:,i)
                d    = sqrt(sum(rij * rij))
                dcut = min(d,M_PI_2)
                pot  = pot - sin(dcut)*sin(dcut)
                force_ij = sin(2.0*dcut) * rij/d
                
                a(:,i) = a(:,i) + force_ij
                a(:,j) = a(:,j) - force_ij
            end do
        end do
        
        do i = 1, n
            v(:,i) = v(:,i) + 0.5*(a_old(:,i)+a(:,i))*dt
            kin = kin + 0.5*sum(v(:,i)*v(:,i))
        end do
        
        energy_old = energy
        energy = pot + kin
        energy_err_acum = energy_err_acum + abs(energy-energy_init)
        energy_err_max  = max(abs(energy-energy_init), energy_err_max)
        
        if( t > write_t) then
            write(foffset, *) t, energy, (energy-energy_init), abs(energy-energy_old)
            
            write_t = write_t + write_dt
        end if
        
    end do
    
    print *, "energy errors:"
    print *, " av. abs. err: ", energy_err_acum/loop_count
    print *, "max. abs. err: ", energy_err_max
    
end program


