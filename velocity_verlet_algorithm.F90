
program velocity_verlet_algorithm
    use Global
    use Helpers
    
    integer                     ::  i, j, n
    integer(8)                  ::  loop_count
    real(rk)                    ::  t, dt, write_t, write_dt 
    real(rk)                    ::  tmax
    
    real(rk), dimension(:,:)    ::  r, v, a, a_old
    allocatable                 ::  r, v, a, a_old
    
    real(rk)                    ::  d, dcut
    real(rk), dimension(dims)   ::  rij, force_ij, lbounds, ubounds
    real(rk)                    ::  kin, pot, energy, energy_old
    real(rK)                    ::  energy_init, energy_err_acum, energy_err_max

    ! set the simulation parameters
    ! (hardcoded within the source for now)
    t    = 0
    dt   = 1e-07
    tmax =  0.01
    write_t  = t
    write_dt = 0.01
    
    lbounds = 0.0
    ubounds = 5.0
    
    ! get particle count from first cmd line argument
    ! and allocate all arrays
    call get_particle_count(n)
    allocate(r(dims,n), v(dims,n), a(dims,n), a_old(dims,n))
    
    ! initialize the intern RNG
    call init_rand()

    ! open the data (output) files
    call open_files(n)
    
    ! set initial values for r, v and a
    ! (see Helpers module)
    call init_rva(r, v, a, lbounds, ubounds)
   
    ! calculate the initial energy
    pot = 0
    kin = 0
    do i = 1, n
        do j = i+1, n
            rij  = r(:,j) - r(:,i)
            d    = sqrt(sum(rij * rij))
            dcut = min(d,pi_2)
            pot  = pot - sin(dcut)*sin(dcut)
            force_ij = sin(2.0*dcut) * rij/d
            
            a(:,i) = a(:,i) + force_ij
            a(:,j) = a(:,j) - force_ij
        end do
        kin = kin + 0.5*sum(v(:,i)*v(:,i))
    end do
    energy      = pot + kin
    energy_init = energy
    
    ! and write it to the file
    call write_energy(t, energy, 0.0, 0.0)
    
    ! initally we've got no errors at all
    energy_err_acum = 0
    energy_err_max  = 0
    
    ! reset the loop counter and start the main loop
    loop_count = 0
    
    !$omp parallel default(none)                                &
    !$omp shared(loop_count,t,dt,tmax,write_t,write_dt)         &
    !$omp shared(r,v,a,a_old)                                   &
    !$omp shared(n,lbounds,ubounds)                             &
    !$omp shared(pot,kin,energy,energy_init)                    &
    !$omp shared(energy_old,energy_err_acum,energy_err_max)
    
    do while( t < tmax )
        !$omp single
        
        ! increment the time
        t = t + dt
        loop_count = loop_count + 1
       
        ! reset energy accumulators
        pot = 0  
        kin = 0 

        !$omp end single
                
        ! update positions (r)
        ! (and clear the acceleration)
        ! also write r to file if neccessary
        !$omp do private(i,j)
        do i = 1, n
            r(:,i) = r(:,i) + v(:,i)*dt + 0.5*a(:,i)*dt*dt
            
            ! boundary checks
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
            
            if( t >= write_t ) &
                call write_single_particle(t, i, r(:,i))
        end do
        !$omp end do
        
        ! update accelerations (a) and pot energy
        !$omp do private(i,j,rij,d,dcut,force_ij) reduction(+:a,pot)
        do i = 1, n
            do j = i+1, n
                rij  = r(:,j) - r(:,i)
                d    = sqrt(sum(rij * rij))
                dcut = min(d,pi_2)
                pot  = pot - sin(dcut)*sin(dcut)
                force_ij = sin(2.0*dcut) * rij/d
                
                a(:,i) = a(:,i) + force_ij
                a(:,j) = a(:,j) - force_ij
            end do
        end do
        !$omp end do
        
        ! update velocities (v) and kin energy
        !$omp do private(i) reduction(+:kin)
        do i = 1, n
            v(:,i) = v(:,i) + 0.5*(a_old(:,i)+a(:,i))*dt
            kin = kin + 0.5*sum(v(:,i)*v(:,i))
        end do
        !$omp end do
        
        !$omp flush
        !$omp single
        
        ! recompute energy (errors)
        energy_old = energy
        energy = pot + kin
        energy_err_acum = energy_err_acum + abs(energy-energy_init)
        energy_err_max  = max(abs(energy-energy_init), energy_err_max)
        
        ! write energy values to file (and increment write_t) if neccessary
        if( t >= write_t) then
            call write_energy(t, energy, (energy-energy_init), abs(energy-energy_old))
            
            write_t = write_t + write_dt
        end if
        
        !$omp end single
        
    end do

    !$omp end parallel

    ! end of main loop
    
    print *, "energy errors:"
    print *, "E(tmax) - E(0): ", energy - energy_init
    print *, "  av. abs. err: ", energy_err_acum/loop_count
    print *, " max. abs. err: ", energy_err_max
    
end program


