module output
use machina_basic, only: i4, f8

implicit none 

contains

subroutine out_title()
    implicit none 

    write(99,'(1x,a,/)')'======================================================================================='
    write(99,'(1x,20x,a,20x,/)') 'LuAIMD: An ab initio molecular dynamics program'
    write(99,'(1x,32x,a,26x,/)') 'Author: Yun-ting Bu'
    write(99,'(1x,a,/)')'======================================================================================='
end subroutine

subroutine out_mole_info(n,basis,charge,multiplicity,atom,x,y,z)
    implicit none
    character(len=*),intent(in) :: basis 
    integer,intent(in) :: charge
    integer,intent(in) :: n
    integer,intent(in) :: multiplicity
    character(len=*),allocatable,dimension(:),intent(in) :: atom
    real(kind=f8),allocatable,dimension(:),intent(in) :: x
    real(kind=f8),allocatable,dimension(:),intent(in) :: y
    real(kind=f8),allocatable,dimension(:),intent(in) :: z 
    integer :: i

    write(99,'(1x,a,a)') 'Basis: ',trim(basis) 
    write(99,'(1x,a,i3,3x,a,i2)') 'Charge =',charge,'Spin multiplicity =',multiplicity
    write(99,'(/,1x,a)') 'Cartesian coordinates (in Angstroms):'
    write(99,'(/,2x,a)') ' Atom           X                 Y                 Z       '
    write(99,'(2x,a)')   '------  ----------------  ----------------  ----------------'
    do i=1,n 
        write(99,'(4x,a2,5x,f12.8,6x,f12.8,7x,f12.8)') atom(i),x(i),y(i),z(i)
    end do 
end subroutine

subroutine out_scf_set(guess,damp)
    implicit none
    character(len=*),intent(in) :: guess
    real(kind=f8),intent(in) :: damp 

    write(99,'(/,1x,a)') 'SCF information:'
    write(99,'(1x,a,a)') 'The init guess: ',trim(guess)
    write(99,'(1x,a)') 'The first two iterations will use DAMP to accelerate the SCF procedure.'
    write(99,'(1x,a,f5.2)') 'The number of damp is: ',damp
    write(99,'(1x,a)') 'If not converged in two cycles, we will start DIIS.'
    end subroutine

subroutine out_md(nstep,time,md_job,init_temp,temp_job,bath_temp,con_time)
    implicit none
    character(len=*),intent(in) :: md_job
    character(len=*),intent(in) :: temp_job
    integer,intent(in) :: nstep
    real(kind=f8),intent(in) :: time 
    real(kind=f8),intent(in) :: bath_temp
    real(kind=f8),intent(in) :: init_temp
    real(kind=f8),intent(in) :: con_time

    write(99,'(/,a)') '======================================='
    write(99,'(10x,a)') 'Molecular Dynamics'
    write(99,'(a,/)') '======================================='
    write(99,'(1x,a,i5)') 'Nstep= ',nstep
    write(99,'(1x,a,f8.2,1x,a)') 'Time= ',time,'fs' 
    write(99,'(1x,a,f8.2,a)') 'Velocity initialization: Maxwell-Boltzmann distribution, at ',init_temp,' K.'
    write(99,'(1x,a,a)') 'MD intergral method: ',trim(md_job)
    write(99,'(1x,a,a)') 'Thermostat: ', trim(temp_job)
    write(99,'(1x,a,f8.2,1x,a)') 'Bath temperature: ',bath_temp,'K'
    write(99,'(1x,a,f6.2,1x,a)') 'Relaxation time: ',con_time,'fs'
    write(99,'(1x,a)') 'Cartesian coordinates will be given in: pos.xyz.'
    write(99,'(1x,a,/)') 'Bond length and bond angle you have chosen will be given in: mole_geom.dat.' 
    write(99,'(a,/)') '------MD step------'

end subroutine

subroutine energy_out(Ent,iter)
    implicit none
    integer,intent(in) :: iter 
    real(kind=f8),intent(in) :: Ent 
    write(*,'(/,1x,a,f16.10,2x,a,i4,a,/)') '[ Total energy = ',Ent,'Hartree, after',iter,' cycles. ]'
end subroutine

subroutine out_mp2(Emp2,Ent)
    implicit none
    real(kind=f8),intent(in) :: Ent
    real(kind=f8),intent(in) :: Emp2
    write(*,'(/,1x,a,f16.10,2x,a,f16.10,a,/)') '[ Total energy = ',Ent+Emp2,&
                                           'Hartree, correlation energy = ',Emp2,' Hartree. ]'
end subroutine

subroutine out_cputime(hour,min,sec,cputime)
    implicit none
    integer,intent(in) :: hour
    integer,intent(in) :: min 
    integer,intent(in) :: sec
    real(kind=f8),intent(in) :: cputime
    write(99,'(/,a)')'=========================================================='
    if (sec /= 0) then 
        write(99,'(/,1x,a,i3,1x,a,1x,i3,1x,a,i3,1x,a)') 'CPU time: ',hour,'h',min,'min',sec,'sec.'
    else 
        write(99,'(/,1x,a,i3,1x,a,1x,i3,1x,a,f10.4,1x,a)') 'CPU time: ',hour,'h',min,'min',cputime,'sec.'
    end if
    write(99,'(1x,a)')'Normal termination, have a nice day! (^ _ ^)'
end subroutine



end module