module scf 
use machina_basic, only: i4, f8
use read_input
use basis 
use cint
use output 
implicit none 

type,public :: diis
      logical :: use_diis = .true.
      integer :: n_errmat = 6
      integer :: n_active = 0
      integer :: iter = 0
      real(kind=f8), allocatable :: e(:,:,:)
      real(kind=f8), allocatable :: F(:,:,:)
      real(kind=f8), allocatable :: B(:,:)
      real(kind=f8), allocatable :: c(:)
      real(kind=f8), allocatable :: rhs(:)

      contains 
        procedure,public :: init_diis => init_diis_sub
        procedure,public :: diis_pro => diis_pro_sub
end type diis

type(mole_input_info) :: mole_in_scf
type(cint_info) :: cint_in_scf 

type,public :: scf_info
    real(kind=f8),allocatable,dimension(:,:) :: Hcore
    real(kind=f8),allocatable,dimension(:,:) :: H_GWH
    real(kind=f8),allocatable,dimension(:,:) :: P 
    real(kind=f8),allocatable,dimension(:,:) :: P_read 
    real(kind=f8),allocatable,dimension(:,:) :: F
    real(kind=f8),allocatable,dimension(:,:) :: C 
    real(kind=f8),allocatable,dimension(:,:) :: C_mp2
    real(kind=f8),allocatable,dimension(:,:) :: G
    real(kind=f8),allocatable,dimension(:) :: E 
    real(kind=f8),allocatable,dimension(:,:) :: X 
    real(kind=f8),allocatable,dimension(:,:) :: Bond_R
    real(kind=f8),allocatable,dimension(:,:,:,:) :: eri_mp2
    real(kind=f8),allocatable,dimension(:,:,:,:) :: eri_MO
    real(kind=f8) :: damp
    real(kind=f8) :: Ee 
    real(kind=f8) :: Enuc
    real(kind=f8) :: Ent 
    real(kind=f8) :: Emp2
    character(len=200) :: job
    integer :: nocc
    integer :: iter_cong

    contains
        procedure,public :: guess => guess_sub
        procedure,public :: form_H_GWH => form_H_GWH_sub
        procedure,public :: form_X => form_X_sub
        procedure,public :: form_F => form_F_sub
        procedure,public :: ex_Huckel => ex_Huckel_sub
        procedure,public :: scf_pro => scf_pro_sub
        procedure,public :: bond_lengths => bond_lengths_sub
        procedure,public :: init_mp2 => init_mp2_sub
        procedure,public :: cal_mp2 => cal_mp2_sub

end type scf_info 

contains

subroutine guess_sub(this,cint_in_scf,guess)
implicit none
    class(scf_info) :: this
    type(cint_info),intent(in) :: cint_in_scf
    character(len=*),intent(in) :: guess 
    integer :: status,i,j

    allocate(this%Hcore(cint_in_scf%nbas,cint_in_scf%nbas),stat=status)
    allocate(this%H_GWH(cint_in_scf%nbas,cint_in_scf%nbas),stat=status)
    allocate(this%P(cint_in_scf%nbas,cint_in_scf%nbas),stat=status)
    allocate(this%P_read(cint_in_scf%nbas,cint_in_scf%nbas),stat=status)
    allocate(this%C(cint_in_scf%nbas,cint_in_scf%nbas),stat=status)
    allocate(this%F(cint_in_scf%nbas,cint_in_scf%nbas),stat=status)
    allocate(this%G(cint_in_scf%nbas,cint_in_scf%nbas),stat=status)
    allocate(this%E(cint_in_scf%nbas),stat=status)
    do i=1,cint_in_scf%nbas 
        do j=1,cint_in_scf%nbas
            this%Hcore(i,j)=cint_in_scf%T(i,j)+cint_in_scf%V(i,j)
            this%P(i,j)=0.0_f8
        end do 
    end do 
    call this%form_H_GWH(cint_in_scf)
    if (guess=='Huckel') then
        call this%ex_Huckel(cint_in_scf)
    else if (guess=='read') then 
        this%P=this%P_read
    else 
        this%P=0.0_f8
    end if
end subroutine

subroutine form_X_sub(this,cint_in_scf,mole_in_scf)
implicit none
    class(scf_info) :: this
    type(cint_info),intent(in) :: cint_in_scf 
    type(mole_input_info),intent(in) :: mole_in_scf
    real(kind=f8),allocatable,dimension(:,:) :: oldS
    real(kind=f8),allocatable,dimension(:,:) :: tempS
    real(kind=f8),allocatable,dimension(:) :: eivuS
    real(kind=f8),allocatable,dimension(:,:) :: eivtS
    real(kind=f8),dimension(3*cint_in_scf%nbas-1) :: w
    integer :: status,i,j,k,info

    allocate(oldS(cint_in_scf%nbas,cint_in_scf%nbas),stat=status)
    allocate(tempS(cint_in_scf%nbas,cint_in_scf%nbas),stat=status)
    allocate(eivuS(cint_in_scf%nbas),stat=status)
    allocate(eivtS(cint_in_scf%nbas,cint_in_scf%nbas),stat=status)
    allocate(this%X(cint_in_scf%nbas,cint_in_scf%nbas),stat=status)
    oldS=cint_in_scf%S 
    call dsyev('V', 'L', cint_in_scf%nbas, oldS, cint_in_scf%nbas, eivuS, w, 3*cint_in_scf%nbas-1, info)
    eivtS=oldS
    do i=1,cint_in_scf%nbas
        do j=1,cint_in_scf%nbas
            tempS(i,j)=eivtS(i,j)/sqrt(eivuS(j))
        end do
    end do
    do i=1,cint_in_scf%nbas
        do j=1,cint_in_scf%nbas
            this%X(i,j)=0.0_f8
            do k=1,cint_in_scf%nbas
                this%X(i, j) = this%X(i, j) + eivtS(i, k)*tempS(j, k)
            end do
        end do
    end do
    deallocate(oldS)
    deallocate(tempS)
    deallocate(eivuS)
    deallocate(eivtS)
    end subroutine

subroutine form_F_sub(this,cint_in_scf,mole_in_scf,iter)
implicit none
    class(scf_info) :: this
    type(cint_info),intent(in) :: cint_in_scf 
    type(mole_input_info),intent(in) :: mole_in_scf
    integer,intent(in) :: iter
    integer :: status,i,j,k,l,info

        do i=1,cint_in_scf%nbas
            do j=1,cint_in_scf%nbas
                this%G(i,j)=0.0_f8
                do k=1,cint_in_scf%nbas
                    do l=1,cint_in_scf%nbas
                        this%G(i,j)=this%G(i,j)+this%P(k,l)*(cint_in_scf%eri(i,j,k,l)-0.5_f8*cint_in_scf%eri(i,l,k,j))
                    end do
                end do
            end do
        end do
        do i=1,cint_in_scf%nbas
            do j=1,cint_in_scf%nbas
                if (trim(mole_in_scf%guess)=='GWH' .and. iter==0) then
                    this%F(i,j)=this%H_GWH(i,j)+this%G(i,j)
                else
                    this%F(i,j)=this%Hcore(i,j)+this%G(i,j)
                end if
            end do 
        end do
end subroutine


subroutine scf_pro_sub(this,cint_in_scf,mole_in_scf)
implicit none
    class(scf_info) :: this
    type(diis) :: diis0
    type(cint_info),intent(in) :: cint_in_scf 
    type(mole_input_info),intent(in) :: mole_in_scf
    real(kind=f8),allocatable,dimension(:,:) :: oldP 
    real(kind=f8),allocatable,dimension(:,:) :: Fprime
    real(kind=f8),allocatable,dimension(:,:) :: Cprime
    real(kind=f8),allocatable,dimension(:,:) :: XT
    real(kind=f8),dimension(3*cint_in_scf%nbas-1) :: w
    real(kind=f8) :: oldEnt
    real(kind=f8) :: delta
    integer :: inocc
    integer :: iter 
    integer :: maxcycle = 128
    integer :: status,i,j,k,l,info

    this%job = mole_in_scf%QC_method
    inocc = 0
    do i=1,cint_in_scf%natom
        inocc = cint_in_scf%charge(i)+inocc
    end do 
    this%nocc = (inocc-mole_in_scf%mole_charge)/2
    this%damp = mole_in_scf%damp

    call this%form_X(cint_in_scf,mole_in_scf)
    call this%guess(cint_in_scf,trim(mole_in_scf%guess))
    allocate(oldP(cint_in_scf%nbas,cint_in_scf%nbas),stat=status)
    allocate(Fprime(cint_in_scf%nbas,cint_in_scf%nbas),stat=status)
    allocate(Cprime(cint_in_scf%nbas,cint_in_scf%nbas),stat=status)
    allocate(XT(cint_in_scf%nbas,cint_in_scf%nbas),stat=status)
    XT=transpose(this%X)
    oldEnt=0.0_f8

    iter = 0
    call this%form_F(cint_in_scf,mole_in_scf,iter)
    call diis0%init_diis(cint_in_scf)

    do iter=0,maxcycle
        this%Ee=0.0_f8
        do i=1,cint_in_scf%nbas
            do j=1,cint_in_scf%nbas
                this%Ee=this%Ee+0.5_f8*this%P(i,j)*(this%Hcore(i,j)+this%F(i,j))
            end do
        end do
        Fprime=matmul(matmul(XT,this%F),this%X)
        call dsyev('V', 'L', cint_in_scf%nbas, Fprime, cint_in_scf%nbas, this%E, w, 3*cint_in_scf%nbas-1,info)
        Cprime=Fprime
        this%C=matmul(this%X,Cprime)
        oldP=this%P
        do i=1,cint_in_scf%nbas
            do j=1,cint_in_scf%nbas
                this%P(i,j)=0.0_f8
            end do
        end do 
        this%P = 2.0_f8*matmul(this%C(:,1:this%nocc), transpose(this%C(:,1:this%nocc)))
        if (iter <= 2) then 
            diis0%use_diis = .false.
            if(iter==2) this%P = this%damp*oldP + (1.0-this%damp)*this%P
        else 
            diis0%use_diis = .true.
        end if
        this%Enuc=0.0_f8
        do i=1,cint_in_scf%natom
            do j=i+1,cint_in_scf%natom
                call this%bond_lengths(cint_in_scf,mole_in_scf)
                this%Enuc=this%Enuc+((cint_in_scf%charge(i)*cint_in_scf%charge(j))/(this%Bond_R(i,j)/0.529177249_f8))
            end do
        end do
        this%Ent=this%Ee+this%Enuc
        if (iter==0) then
            delta = 0.0_f8
        else
            delta=abs(this%Ent-oldEnt)
        end if
        if (delta<1.0D-8 .and. iter /= 0) exit
        !Write (*, '(1x,a7,f16.10,2x,a,2x,a9,f16.8,2x,a,i3)') 'E(HF)= ', this%Ent,'a.u.', 'delta= ', delta,'Cycle=',iter
        oldEnt=this%Ent
        call this%form_F(cint_in_scf,mole_in_scf,iter)
        call diis0%diis_pro(this%F,this%P,cint_in_scf%S)
    end do
    !write(*,'(/,1x,a,/)') '==> Converged Energies <=='
    !write(*,'(1x,a,f16.10,2x,a)') 'Electronic energy = ',this%Ee,'a.u.'
    !write(*,'(1x,a,f16.10,2x,a)') 'Nuclear repulsive energy = ',this%Enuc,'a.u.'
    this%iter_cong = iter
    this%P_read = this%P
    deallocate(oldP)
    deallocate(Fprime)
    deallocate(Cprime)
    deallocate(XT)
    end subroutine

    subroutine form_H_GWH_sub(this,cint_in_scf)
    implicit none
        class(scf_info) :: this
        type(cint_info),intent(in) :: cint_in_scf
        integer :: i,j,status

        do i=1,cint_in_scf%nbas
            do j=1,cint_in_scf%nbas
                this%H_GWH(i,j) = 0.5_f8*1.75_f8 &
                                  *(this%Hcore(i,i)+this%Hcore(j,j))*cint_in_scf%S(i,j)
            end do
        end do 
    end subroutine

    subroutine ex_Huckel_sub(this,cint_in_scf)
    implicit none
        class(scf_info) :: this
        type(cint_info),intent(in) :: cint_in_scf
        real(kind=f8),allocatable,dimension(:,:) :: Hprime
        real(kind=f8),allocatable,dimension(:,:) :: Cprime
        real(kind=f8),allocatable,dimension(:) :: E
        real(kind=f8),dimension(3*cint_in_scf%nbas-1) :: w
        integer :: i,j,status,info

        allocate(Hprime(cint_in_scf%nbas,cint_in_scf%nbas),stat=status)
        allocate(Cprime(cint_in_scf%nbas,cint_in_scf%nbas),stat=status)
        allocate(E(cint_in_scf%nbas),stat=status)
        Hprime=matmul(matmul(this%X,this%H_GWH),this%X)
        Call dsyev('V', 'L', cint_in_scf%nbas, Hprime, cint_in_scf%nbas, E, w, 3*cint_in_scf%nbas-1,info)
        Cprime=Hprime
        this%C=matmul(this%X,Cprime)
        do i=1,cint_in_scf%nbas
            do j=1,cint_in_scf%nbas
                this%P(i,j)=0.0_f8
            end do
        end do 
        this%P = 2.0_f8*matmul(this%C(:,1:this%nocc), transpose(this%C(:,1:this%nocc)))
        deallocate(Hprime)
        deallocate(Cprime)
        deallocate(E)
    end subroutine

    subroutine bond_lengths_sub(this,cint_in_scf,mole_in_scf)
    implicit none
        class(scf_info) :: this
        type(cint_info),intent(in) :: cint_in_scf 
        type(mole_input_info),intent(in) :: mole_in_scf
        integer :: i,j,status

        allocate (this%Bond_R(cint_in_scf%natom,cint_in_scf%natom),stat=status)
        do i=1,cint_in_scf%natom
            do j=1,cint_in_scf%natom
                this%Bond_R(i,j)=sqrt(((mole_in_scf%x_list(i)-mole_in_scf%x_list(j))**2.0) &
                +((mole_in_scf%y_list(i)-mole_in_scf%y_list(j))**2.0)+((mole_in_scf%z_list(i)-mole_in_scf%z_list(j))**2.0))
            end do 
        end do 
    end subroutine    

    subroutine init_diis_sub(this, cint_in_scf)
    implicit none
        class(diis) :: this
        type(cint_info),intent(in) :: cint_in_scf 

        if (this%n_errmat <2 ) then
        this%use_diis = .false.
        else
            associate(n=>cint_in_scf%nbas)
               allocate(this%e(n,n,this%n_errmat), source=0.0_f8)
               allocate(this%F(n,n,this%n_errmat), source=0.0_f8)
            end associate
        end if
    end subroutine init_diis_sub

    subroutine diis_pro_sub(this, F, P, S)
    implicit none
        class(diis) :: this
        real(kind=f8), intent(inout) :: F(:,:)
        real(kind=f8), intent(in) :: S(:,:), P(:,:)
        real(kind=f8), allocatable :: work(:)
        real(kind=f8), allocatable :: ipiv(:)
        integer :: i, j, ierr, lwork, il

        lwork = -1
        if (this%use_diis) then
        this%iter = this%iter+1
        if (this%iter > this%n_errmat) this%iter = this%iter - this%n_errmat
        if (this%n_active < this%n_errmat) this%n_active = this%n_active+1
        this%F(:,:,this%iter) = F(:,:)
        this%e(:,:,this%iter) = matmul(F, matmul(P,S)) &
                                - matmul(S, matmul(P, F))
        associate(n=>this%n_active, nerr=>this%n_errmat)
        if (n > 1) then
        ! construct the B matrix
            if (n <= nerr) then
               if (allocated(this%B)) deallocate(this%B, this%c, this%rhs)
               allocate(this%B(n+1,n+1), this%c(n+1), this%rhs(n+1), source=0.0_f8)
            end if
            this%B(n+1,:) = -1.0_f8
            this%B(n+1,n+1) = 0.0_f8
            this%rhs(n+1) = -1.0_f8
            this%c = this%rhs
            do i = 1, n
               do j = 1, i
                  this%B(i,j) = sum(this%e(:,:,i)*this%e(:,:,j))
               end do
            end do
            do i = 1, n+1
            end do
            allocate(ipiv(size(this%c)))
            do il = 1, 2
                allocate(work(abs(lwork)))
                call dsysv('L', size(this%c), 1, this%B, size(this%c), ipiv, this%c, size(this%c), work, lwork, ierr)
                lwork = nint(work(1))
                deallocate(work)
            end do
            F = 0.0_f8
            do i = 1, n
               F = F + this%c(i) * this%F(:,:,i)
            end do
         end if
         end associate
         end if
      end subroutine

    subroutine init_mp2_sub(this,cint_in_scf)
    implicit none
    class(scf_info) :: this
    type(cint_info),intent(in) :: cint_in_scf

    allocate(this%C_mp2(cint_in_scf%nbas,cint_in_scf%nbas))
    allocate(this%eri_mp2(cint_in_scf%nbas,cint_in_scf%nbas,cint_in_scf%nbas,cint_in_scf%nbas))
    allocate(this%eri_MO(cint_in_scf%nbas,cint_in_scf%nbas,cint_in_scf%nbas,cint_in_scf%nbas))
    this%C_mp2 = this%C 
    this%eri_mp2 = cint_in_scf%eri
    
    end subroutine

    subroutine cal_mp2_sub(this,cint_in_scf)
    implicit none
    class(scf_info) :: this
    type(cint_info),intent(in) :: cint_in_scf 
    real(kind=f8),allocatable,dimension(:,:,:,:) :: a_m
    real(kind=f8),allocatable,dimension(:,:,:,:) :: b_m
    real(kind=f8),allocatable,dimension(:,:,:,:) :: c_m
    real(kind=f8) :: temp
    integer :: i,j,k,l,r,p,q,s,a,b,n 

    allocate(a_m(cint_in_scf%nbas,cint_in_scf%nbas,cint_in_scf%nbas,cint_in_scf%nbas))
    allocate(b_m(cint_in_scf%nbas,cint_in_scf%nbas,cint_in_scf%nbas,cint_in_scf%nbas))
    allocate(c_m(cint_in_scf%nbas,cint_in_scf%nbas,cint_in_scf%nbas,cint_in_scf%nbas))

    do l = 1, cint_in_scf%nbas 
        do r = 1, cint_in_scf%nbas 
            do q = 1, cint_in_scf%nbas 
                do p = 1, cint_in_scf%nbas 
                    a_m(p, q, r, l) = 0.0_f8
                    do s = 1, cint_in_scf%nbas 
                        a_m(p, q, r, l) = a_m(p, q, r, l) + this%C_mp2(s, l)*this%eri_mp2(p, q, r, s)
                    end do
                end do
            end do
        end do
    end do
    do l = 1, cint_in_scf%nbas 
        do k = 1, cint_in_scf%nbas 
            do q = 1, cint_in_scf%nbas 
                do p = 1, cint_in_scf%nbas 
                    b_m(p, q, k, l) = 0.0_f8
                    do r = 1, cint_in_scf%nbas
                        b_m(p, q, k, l) = b_m(p, q, k, l) + this%C_mp2(r, k)*a_m(p, q, r, l)
                    end do
                end do
            end do
        end do
    end do
    do l = 1, cint_in_scf%nbas 
        do k = 1, cint_in_scf%nbas 
            do j = 1, cint_in_scf%nbas 
                do p = 1, cint_in_scf%nbas
                    c_m(p, j, k, l) = 0.0_f8
                    do q = 1, cint_in_scf%nbas
                        c_m(p, j, k, l) = c_m(p, j, k, l) + this%C_mp2(q, j)*b_m(p, q, k, l)
                    end do
                end do
            end do
        end do
    end do
    do l = 1, cint_in_scf%nbas 
        do k = 1, cint_in_scf%nbas 
            do j = 1, cint_in_scf%nbas 
                do i = 1, cint_in_scf%nbas
                    this%eri_MO(i, j, k, l) = 0.0_f8
                    do p = 1, cint_in_scf%nbas
                        this%eri_MO(i, j, k, l) = this%eri_MO(i, j, k, l) + this%C_mp2(p, i)*c_m(p, j, k, l)
                    end do
                end do
            end do
        end do
    end do
    n = 0
    if (trim(this%job)=='mp2(fc)') then 
        do i = 1,cint_in_scf%natom
            n = n + cint_in_scf%fz_core(i)
        end do 
    else 
        n = 1
    end if 
    temp = 0.0_f8
    do k = n, this%nocc
        do l = n, this%nocc
            do a = this%nocc + 1, cint_in_scf%nbas
                do b = this%nocc + 1, cint_in_scf%nbas
                    temp = temp + this%eri_MO(k, a, l, b) &
                                    *(2.0_f8*this%eri_MO(k,a,l,b)-this%eri_MO(k,b,l,a)) &
                                    /(-this%E(a)-this%E(b)+this%E(l)+this%E(k))
                End Do
            End Do
        End Do
    End Do
    this%Emp2 = temp
    deallocate(a_m)
    deallocate(b_m)
    deallocate(c_m)
    deallocate(this%C_mp2)
    deallocate(this%eri_mp2)
    deallocate(this%eri_MO)
    end subroutine
    

end module


