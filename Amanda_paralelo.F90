
program Argonio
   
   use omp_lib
   implicit none
   
!!! vetores   
   real(8), allocatable ::   x(:),   y(:),   z(:)
   real(8), allocatable ::  vx(:),  vy(:),  vz(:)
   real(8), allocatable :: vx0(:), vy0(:), vz0(:)
   real(8), allocatable ::  ax(:),  ay(:),  az(:)
   real(8), allocatable :: Etot(:), Ecin(:), Temp(:), Pres(:)

   
!!!variáveis tiradas do arquivo
   real(8)		:: dens, temp0, dt
   integer(4)           :: ncell, nstep, Nm
!!!
   integer(4)           :: natm, t, stepRm, r,r2
   real(8)		:: Epot, dcell, Rcut, Rm, dx, dy, dz
!!! para os loops
   integer(4)           :: iatm, istep, istepRm, aux, i, j, auxstep 
!!!para contar os tempos
   integer              :: start,start1, finish,finish1, clock_rate,clock_rate1, clock_max, clock_max1
   real                 :: sumtime


   call system_clock(start1, clock_rate1, clock_max1)
   open(unit=12, file="tempos.txt")  
  
!!!pegamos os dados do arquivo de entrada
   open(10,file="entrada.dat")
   read(10,*) dens      ! densidade
   read(10,*) Temp0 		! temperatura
   read(10,*) dt 	        ! delta tempo
   read(10,*) ncell 		! número de células
   read(10,*) nstep 		! quantidade de passos
   read(10,*) rcut              ! raio de corte
   close(10)
   
   
   natm = ncell**3 ! cálculo da quantidade de partículas
   dcell = ncell / dens**(1.d0/3.d0) ! 
   nm=100
   auxstep = nstep/Nm
   allocate(x(natm), y(natm), z(natm))
   allocate(vx(natm), vy(natm), vz(natm))
   allocate(vx0(natm), vy0(natm), vz0(natm))
   allocate(ax(natm), ay(natm), az(natm))
   allocate(Etot(auxstep), Ecin(auxstep),Temp(auxstep),Pres(auxstep))
   
   call system_clock(start, clock_rate, clock_max)
  
   call ConfigInicial(x,y,z,vx,vy,vz,ax,ay,az,natm,dens,ncell,Temp0)
   
   vx0 = vx
   vy0 = vy
   vz0 = vz
   
   call system_clock(finish, clock_rate, clock_max)
   
   write(12,*) "Configuração Inicial", real(finish-start)/real(clock_rate)
   
   call system_clock(start, clock_rate, clock_max)
   sumtime=0.0d0 
   open(unit=14, file="saidagraf.dat")
   t = 0.d0
   do istep = 1,auxstep
     
      do isteprm=1, nm
         t = t + dt                                                  
         call Verlet(x,y,z,vx,vy,vz,ax,ay,az,Rcut,natm,dt,dens,dcell,Ecin(istep),&
      							Epot,Pres(istep),sumtime)
      end do
      Etot(istep) = Ecin(istep) + Epot
      Temp(istep) = natm * Ecin(istep) / (0.5d0*(3.d0*natm-3.d0))
      write(14, *) istep, Temp(istep), pres(istep), Ecin(istep), Etot(istep)
   end do
   close(14)
   
   call system_clock(finish, clock_rate, clock_max)
   write(12,*) "Verlet", real(finish-start)/real(clock_rate)
   write(12,*) "Força", sumtime
   call system_clock(finish1, clock_rate1, clock_max1)

   write(12,*) "tempo total", real(finish1-start1)/real(clock_rate1) 
   close(14)
   deallocate(x, y, z)
   deallocate(vx, vy, vz)
   deallocate(vx0, vy0, vz0)
   deallocate(ax, ay, az)
   deallocate(Etot, Ecin,Temp,Pres)
   
end program Argonio

!===============================================================================
subroutine ConfigInicial(x,y,z,vx,vy,vz,ax,ay,az,natm,dens,ncell,T)
!-------------------------------------------------------------------------------

   implicit none

   
   real(8)	       ::  x(natm),  y(natm),  z(natm)
   real(8)         :: vx(natm), vy(natm), vz(natm)
   real(8)         :: ax(natm), ay(natm), az(natm)
   
   integer(4)      :: natm, ncell
   real(8)	       :: dens, T

   real(8)         :: h, xi, yj, zk, Ecin
   integer(8)      :: iatm, i, j, k
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   h = 1.d0 / dens**(1.d0/3.d0)
   iatm = 0
   do i = 1,ncell
      xi = (i - 0.5) * h
      do j = 1,ncell
         yj = (j - 0.5) * h
         do k = 1,ncell
            zk = (k - 0.5) * h
            iatm = iatm + 1
            x(iatm) = xi
            y(iatm) = yj
            z(iatm) = zk
         end do
      end do
   end do

   Ecin = 0.5d0 * (3*natm - 6) * T
   call Heat(1.d0,vx,vy,vz,natm,Ecin)

   do iatm = 1,natm
      ax(iatm) = 0.d0; ay(iatm) = 0.d0; az(iatm) = 0.d0
   end do
end subroutine
!----------------------------------------------------------------------------

!============================================================================
subroutine MoveCM(x,y,z,natm)
!----------------------------------------------------------------------------
   implicit none
   
    
   real(8)     :: x(natm), y(natm), z(natm)
   real(8)     :: xCM, yCM, zCM
   integer(4)  :: iatm, natm
   xCM = 0.d0
   yCM = 0.d0
   zCM = 0.d0
   do iatm = 1,natm
      xCM = xCM + x(iatm)
      yCM = yCM + y(iatm)
      zCM = zCM + z(iatm)
   end do
   xCM = xCM / natm
   yCM = yCM / natm
   zCM = zCM / natm

   do iatm = 1,natm
      x(iatm) = x(iatm) - xCM
      y(iatm) = y(iatm) - yCM
      z(iatm) = z(iatm) - zCM
   end do
end

!============================================================================
subroutine Vrescal(amass,vx,vy,vz,natm,Ecin)
!----------------------------------------------------------------------------

   implicit none
	
   real(8)      :: vx(natm), vy(natm), vz(natm)
   integer(4)   :: natm, iatm
   real(8)      :: Ecin0, Ecin, amass, f

   Ecin0 = 0.d0
   do iatm = 1,natm
      Ecin0 = Ecin0 + vx(iatm)*vx(iatm) + vy(iatm)*vy(iatm) + vz(iatm)*vz(iatm)
   end do
   Ecin0 = 0.5d0 * amass *Ecin0

   f = sqrt(Ecin/Ecin0)                                     
   do iatm = 1,natm
      vx(iatm) = f * vx(iatm)
      vy(iatm) = f * vy(iatm)
      vz(iatm) = f * vz(iatm)
   end do
end

!============================================================================
subroutine Heat(amass,vx,vy,vz,natm,Ecin)

   implicit none

   real(8)       :: vx(natm), vy(natm), vz(natm)
   real(8)       :: ex, ey, ez
   integer(4)    :: iatm, natm
   real(8)       :: Ecin, amass

   do iatm = 1,natm
      call RandVec(ex,ey,ez)
      vx(iatm) = ex
      vy(iatm) = ey
      vz(iatm) = ez
   end do

   call MoveCM(vx,vy,vz,natm)
                                                                
   call Vrescal(amass,vx,vy,vz,natm,Ecin)
end


!===============================================================================
subroutine RandVec(e1,e2,e3)
!-------------------------------------------------------------------------------

   implicit none
   
   real(8)     :: f, rnd
   real(8)     :: e1, e2, e3
   f = 2.d0
   do
      if (f <= 1.d0) exit
      call random_number(rnd); e1 = 2.d0* rnd - 1.d0
      call random_number(rnd); e2 = 2.d0* rnd - 1.d0
      f = e1**2 + e2**2
   end do
   e3 = 1.d0 - 2.d0 * f; f = 2.d0 * sqrt(1.d0 - f)
   e1 = f * e1
   e2 = f * e2
end

!============================================================================
subroutine Verlet(x,y,z,vx,vy,vz,ax,ay,az,Rcut,natm,dt,dens,dcell,Ecin,Epot,pres,sumtime)
!----------------------------------------------------------------------------

   implicit none

   real(8)     ::  x(natm),  y(natm),  z(natm)
   real(8)     :: vx(natm), vy(natm), vz(natm)
   real(8)     :: ax(natm), ay(natm), az(natm)
   
   integer(4)  :: natm, iatm
   real(8)     :: dt, dt2, dens, dcell, rcut
   real(8)     :: Ecin, Epot, pres, virial
   
   integer     :: start, finish, clock_rate, clock_max
   real        :: sumtime

   dt2 = 0.5d0 * dt
   vx = vx + ax * dt2                        ! v(t+dt/2)
   vy = vy + ay * dt2
   vz = vz + az * dt2
   x = x + vx * dt                             ! r(t+dt)
   y = y + vy * dt
   z = z + vz * dt
  
   x = x - dcell*floor(x/dcell)
   y = y - dcell*floor(y/dcell)
   z = z - dcell*floor(z/dcell)
 
   call system_clock(start, clock_rate, clock_max)

   call Forcas(x,y,z,ax,ay,az,Rcut,natm,dcell,Epot,virial)   ! f(t+dt)

   call system_clock(finish, clock_rate, clock_max)
   sumtime = sumtime + real(finish-start)/real(clock_rate)

   vx = vx + ax * dt2                        ! v(t+dt)
   vy = vy + ay * dt2
   vz = vz + az * dt2
   
   Ecin = sum(vx(:)*vx(:)+vy(:)*vy(:)+vz(:)*vz(:))

   Ecin = 0.5d0 * Ecin / natm

   pres = dens * (2.d0*Ecin + virial) / 3.d0
end
!---------------------------------------------------------------------------

!============================================================================
subroutine Forcas(x,y,z,fx,fy,fz,Rcut,natm,dcell,Epot,virial)
!----------------------------------------------------------------------------
!  ENERGIA: eV
!  FORÇA eV/Angstroem
!  DISTÂNCIA Angstroem

    implicit none
    
!!!eps = epsilon   
   real, parameter  :: eps = 1.d0, sigma = 1.d0, sigma2 = sigma * sigma	

   real(8)          ::  x(natm),  y(natm),  z(natm)
   real(8)          :: fx(natm), fy(natm), fz(natm)
   real(8)          :: fr2, fr6, fpr, r2
   real(8)          :: dcell, rcut
    
   integer(4)       :: natm, iatm, jatm, j
   real(8)          :: dx, dy, dz
   real(8)          :: Epot, virial
   real             :: Rcut2

   
   fx=0.0d0; fy=0.0d0; fz=0.d0
   Rcut2 = Rcut * Rcut
   Epot   = 0.d0
   virial = 0.d0

!$omp parallel do private(iatm,jatm,dx,dy,dz,r2,fr2,fr6,fpr), &
!$omp & shared(natm,x,y,z,rcut2,dcell),reduction(+:epot,virial,fx,fy,fz) &
!$omp & schedule(runtime) 
   do iatm = 1,natm-1
      do jatm = iatm+1,natm
         dx = x(iatm) - x(jatm)
         dy = y(iatm) - y(jatm)
         dz = z(iatm) - z(jatm)
         ! sign(x,y) == fornece valor positivo de x se y>= 0 e negativo se y < 0
         if (abs(dx) > 0.5d0*dcell) dx = dx - sign(dcell,dx)
         if (abs(dy) > 0.5d0*dcell) dy = dy - sign(dcell,dy)
         if (abs(dz) > 0.5d0*dcell) dz = dz - sign(dcell,dz)

         r2 = dx*dx + dy*dy + dz*dz

         if (r2 < Rcut2) then
            fr2 = sigma2 / r2
            fr6 = fr2 * fr2 * fr2
            fpr = 48.d0 * eps * fr6 * (fr6 - 0.5d0) / r2            ! f/r

            fx(iatm) = fx(iatm) + fpr * dx; fx(jatm) = fx(jatm) - fpr * dx
            fy(iatm) = fy(iatm) + fpr * dy; fy(jatm) = fy(jatm) - fpr * dy
            fz(iatm) = fz(iatm) + fpr * dz; fz(jatm) = fz(jatm) - fpr * dz

            Epot = Epot + 4.d0 * eps * fr6 * (fr6 - 1.d0)
            virial = virial + fpr * r2

         end if
      end do
   end do
!$omp end parallel do
   Epot   = Epot   / natm
   virial = virial / natm
end
!=================================================================================
