
module Helpers
    use iso_fortran_env, only: ERROR_UNIT
    use Global
    
    integer,      parameter     ::  foffset = 20
    character(*), parameter     ::  outpath = "output/"
    
contains

    subroutine init_rand()
        integer                 ::  i, sysclock, seed_size
        integer, allocatable    ::  seed(:)
        
        call system_clock(count=sysclock)
        call random_seed(size=seed_size)
        allocate(seed(seed_size))
        seed = (/ (sysclock+i, i=1,seed_size) /)
        call random_seed(put=seed)
        deallocate(seed)
    end subroutine
    
    subroutine get_particle_count(n)
        integer, intent(out)     ::  n
        integer                  ::  stat
        character(4)             ::  arg
        
        call get_command_argument(1, arg, status=stat)
        call assert("couldn't receive particle count from first argument", &
                    stat == 0)
        read(arg, *) n
        if(n <= 0) n = 1
    end subroutine
    
    subroutine open_files(n)
        integer, intent(in)  ::  n
        integer              ::  i
        character(32)        ::  tmp_file_str
        
        do i = 1, n
            write(tmp_file_str, "(a,a,i0.4,a)") outpath, "particle_", i, ".dat"
            open(foffset+i, file=trim(tmp_file_str), status='replace')
        end do
        
        open(foffset, file=outpath//"energy.dat", status='replace')
    end subroutine
    
    subroutine init_rva(r, v, a, lb, ub)
        real(rk), dimension(:,:),  intent(out)    ::  r, v, a
        real(rk), dimension(dims), intent(in)     ::  lb, ub

        integer  ::  i
        
        ! place particles randomly within the volume
        call random_number(r)
        do i = 1, dims
            r(i,:) = lb(i)+(ub(i)-lb(i))*r(i,:)
        end do

        ! just give the particles a small initial velocity for now
        call random_number(v)

        ! set initial acceleration to zero
        a = 0
    end subroutine
    
    subroutine write_single_particle(t, i, ri)
        real(rk),               intent(in)    ::  t
        integer,                intent(in)    ::  i
        real(rk), dimension(:), intent(in)    ::  ri 
        
        write(foffset+i, *) t, ri
    end subroutine
    
    subroutine write_particles(t, r)
        real(rk),                 intent(in)    ::  t
        real(rk), dimension(:,:), intent(in)    ::  r 

        integer  ::  i

        do i = 1, size(r,2)
            call write_single_particle(t, i, r(:,i))
        end do
    end subroutine
    
    subroutine write_energy(t, energy, abserr, changeerr)
        real(rk), intent(in)    ::  t, energy, abserr, changeerr

        write(foffset, *) t, energy, abserr, changeerr
    end subroutine
        
    subroutine assert(str, bool)
        character(*), intent(in)    ::  str
        logical,      intent(in)    ::  bool
        
        if( .not. bool ) then
            write(ERROR_UNIT, "(a,a)") "Error: ", str
            stop 1
        end if
    end subroutine

    
end module







