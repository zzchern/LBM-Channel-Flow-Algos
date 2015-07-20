!===================================================================    
      subroutine collision_MRT1
      use mpi
      use var_inc
      implicit none 

      real, dimension(0:npop-1) :: Fbar
      real rho9, ux9, uy9, uz9, ux9s, uy9s, uz9s
      real t1, tl1, tl2, tl3, tl4, tl5, tl6, tl7, tl8, tl9, tl10, tl11,&
           tl12, tl13, tl14, tl15, tl16, tl17, tl18, tl19, tl20, tl21
      real eqm1, eqm2, eqm3, eqm4, eqm5, eqm6, eqm7, eqm8, eqm9, eqm10,&
           eqm11, eqm12, eqm13, eqm14, eqm15 
      real sum1, sum2, sum3, sum4, sum5, sum6, sum7, sum8, sum9, sum10,&
           sum11
      real evlm1, evlm2, evlm3, evlm4, evlm5, evlm6, evlm7, evlm8,     &
           evlm9, evlm10, evlm11, evlm12, evlm13, evlm14, evlm15
      real eqmc1, eqmc2, eqmc3, eqmc4, eqmc5, eqmc6, eqmc7, eqmc8,     &
           eqmc9, eqmc10, eqmc11, eqmc12, eqmc13, eqmc14, eqmc15
      real suma, sumb, sumc, sumd, sume, sumf, sumg, sumh, sumi, sumk, &
           sump, sum67, sum89, sum1011
      integer ix, iy, iz, imove, jmove, kmove, sp
      integer ip, ilen
      integer status_array(MPI_STATUS_SIZE,4), req(4)
      real fx9,fy9,fz9,G1,G2,G3

      do iz = 1,lz
      do iy = 1,ly
      do ix = 1,lx

        !Determining macroscopiv variables
        rho9 = rho(ix,iy,iz)
        ux9 = ux(ix,iy,iz)
        uy9 = uy(ix,iy,iz)
        uz9 = uz(ix,iy,iz)
        ux9s = ux9*ux9
        uy9s = uy9*uy9
        uz9s = uz9*uz9

        !I'm merging the forcing term with the collision operator
        fx9 = force_realx(ix,iy,iz)
        fy9 = force_realy(ix,iy,iz)
        fz9 = force_realz(ix,iy,iz)
        G3 = ux9*fx9 + uy9*fy9 + uz9*fz9

        Fbar(0) = -G3

        do ip=1,6
          G1 = cix(ip)*fx9 + ciy(ip)*fy9 + ciz(ip)*fz9
          G2 = cix(ip)*ux9 + ciy(ip)*uy9 + ciz(ip)*uz9
          Fbar(ip) = ww1*(3.*G1 + 9.*G1*G2 - 3.*G3)
        enddo

        do ip=7,(npop-1)
          G1 = cix(ip)*fx9 + ciy(ip)*fy9 + ciz(ip)*fz9
          G2 = cix(ip)*ux9 + ciy(ip)*uy9 + ciz(ip)*uz9
          Fbar(ip) = ww2*(3.*G1 + 9.*G1*G2 - 3.*G3)
        enddo 

        f(:,ix,iy,iz) = f(:,ix,iy,iz) + Fbar(:)/2.

        t1 = ux9s + uy9s + uz9s
        eqm1 = -11.0*rho9 + 19.0*t1
        eqm2 = omegepsl*rho9 + omegepslj*t1
        eqm3 = coef1*ux9
        eqm4 = coef1*uy9
        eqm5 = coef1*uz9
        eqm6 = 2.0*ux9s - uy9s - uz9s
        eqm7 = omegxx*eqm6
        eqm8 = uy9s - uz9s
        eqm9 = omegxx*eqm8
        eqm10 = ux9*uy9
        eqm11 = uy9*uz9
        eqm12 = ux9*uz9
        eqm13 = 0.0
        eqm14 = 0.0
        eqm15 = 0.0

        sum1 = f(1,ix,iy,iz) + f(2,ix,iy,iz) + f(3,ix,iy,iz) &
               + f(4,ix,iy,iz) + f(5,ix,iy,iz) + f(6,ix,iy,iz)        
        sum2 = f(7,ix,iy,iz) + f(8,ix,iy,iz) + f(9,ix,iy,iz)&
               + f(10,ix,iy,iz) + f(11,ix,iy,iz) + f(12,ix,iy,iz) &
               + f(13,ix,iy,iz) + f(14,ix,iy,iz) + f(15,ix,iy,iz) &
               + f(16,ix,iy,iz) + f(17,ix,iy,iz) + f(18,ix,iy,iz)      
        sum3 = f(7,ix,iy,iz) - f(8,ix,iy,iz) + f(9,ix,iy,iz) &
               - f(10,ix,iy,iz) + f(11,ix,iy,iz) - f(12,ix,iy,iz) &
               + f(13,ix,iy,iz) - f(14,ix,iy,iz)
        sum4 = f(7,ix,iy,iz) + f(8,ix,iy,iz) - f(9,ix,iy,iz) &
               - f(10,ix,iy,iz) + f(15,ix,iy,iz) - f(16,ix,iy,iz) &
               + f(17,ix,iy,iz) - f(18,ix,iy,iz)
        sum5 = f(11,ix,iy,iz) + f(12,ix,iy,iz) - f(13,ix,iy,iz) &
               - f(14,ix,iy,iz) + f(15,ix,iy,iz) + f(16,ix,iy,iz) &
               - f(17,ix,iy,iz) - f(18,ix,iy,iz)
        sum6 = f(1,ix,iy,iz) + f(2,ix,iy,iz)
        sum7 = f(3,ix,iy,iz) + f(4,ix,iy,iz) + f(5,ix,iy,iz) &
               + f(6,ix,iy,iz)
        sum8 = f(7,ix,iy,iz) + f(8,ix,iy,iz) + f(9,ix,iy,iz) &
               + f(10,ix,iy,iz) + f(11,ix,iy,iz) + f(12,ix,iy,iz) &
               + f(13,ix,iy,iz) + f(14,ix,iy,iz)
        sum9 = f(15,ix,iy,iz) + f(16,ix,iy,iz) + f(17,ix,iy,iz) + f(18,ix,iy,iz)
        sum10 = f(3,ix,iy,iz) + f(4,ix,iy,iz) - f(5,ix,iy,iz) - f(6,ix,iy,iz)
        sum11 = f(7,ix,iy,iz) + f(8,ix,iy,iz) + f(9,ix,iy,iz) &
                + f(10,ix,iy,iz) - f(11,ix,iy,iz) - f(12,ix,iy,iz) &
              - f(13,ix,iy,iz) - f(14,ix,iy,iz)
    
        evlm1 = -30.0*f(0,ix,iy,iz) + coef2*sum1 + coef3*sum2
        evlm2 = 12.0*f(0,ix,iy,iz) + coef4*sum1 + sum2
        evlm3 = coef4*(f(1,ix,iy,iz) - f(2,ix,iy,iz)) + sum3
        evlm4 = coef4*(f(3,ix,iy,iz) - f(4,ix,iy,iz)) + sum4
        evlm5 = coef4*(f(5,ix,iy,iz) - f(6,ix,iy,iz)) + sum5
        evlm6 = coef5*sum6 - sum7 + sum8 - coef5*sum9
        evlm7 = coef4*sum6 + coef5*sum7 + sum8 - coef5*sum9
        evlm8 = sum10 + sum11
        evlm9 =-coef5*sum10 + sum11
        evlm10 = f(7,ix,iy,iz) - f(8,ix,iy,iz) - f(9,ix,iy,iz) + f(10,ix,iy,iz)
        evlm11 = f(15,ix,iy,iz) - f(16,ix,iy,iz) - f(17,ix,iy,iz) + f(18,ix,iy,iz)
        evlm12 = f(11,ix,iy,iz) - f(12,ix,iy,iz) - f(13,ix,iy,iz) + f(14,ix,iy,iz)
        evlm13 = f(7,ix,iy,iz) - f(8,ix,iy,iz) + f(9,ix,iy,iz) &
               - f(10,ix,iy,iz) - f(11,ix,iy,iz) + f(12,ix,iy,iz) &
               - f(13,ix,iy,iz) + f(14,ix,iy,iz)
        evlm14 =-f(7,ix,iy,iz) - f(8,ix,iy,iz) + f(9,ix,iy,iz) &
                + f(10,ix,iy,iz) + f(15,ix,iy,iz) - f(16,ix,iy,iz) &
               + f(17,ix,iy,iz) - f(18,ix,iy,iz)
        evlm15 = f(11,ix,iy,iz) + f(12,ix,iy,iz) - f(13,ix,iy,iz) &
                - f(14,ix,iy,iz) - f(15,ix,iy,iz) - f(16,ix,iy,iz) &
                + f(17,ix,iy,iz) + f(18,ix,iy,iz)
   
        eqmc1 = evlm1 - s1*(evlm1 - eqm1)
        eqmc2 = evlm2 - s2*(evlm2 - eqm2)
        eqmc3 = evlm3 - s4*(evlm3 - eqm3)
        eqmc4 = evlm4 - s4*(evlm4 - eqm4)
        eqmc5 = evlm5 - s4*(evlm5 - eqm5)
        eqmc6 = evlm6 - s9*(evlm6 - eqm6)
        eqmc7 = evlm7 - s10*(evlm7 - eqm7)
        eqmc8 = evlm8 - s9*(evlm8 - eqm8)
        eqmc9 = evlm9 - s10*(evlm9 - eqm9)
        eqmc10 = evlm10 - s13*(evlm10 - eqm10)
        eqmc11 = evlm11 - s13*(evlm11 - eqm11)
        eqmc12 = evlm12 - s13*(evlm12 - eqm12)
        eqmc13 = evlm13 - s16*(evlm13 - eqm13)
        eqmc14 = evlm14 - s16*(evlm14 - eqm14)
        eqmc15 = evlm15 - s16*(evlm15 - eqm15)


        tl1 = val1i*rho9
        tl2 = coef2*val2i*eqmc1
        tl3 = coef3*val2i*eqmc1
        tl4 = coef4*val3i*eqmc2
        tl5 = val3i*eqmc2
        tl6 = val4i*ux9
        tl7 = val5i*eqmc3
        tl8 = val4i*uy9
        tl9 = val5i*eqmc4
        tl10 = val4i*uz9
        tl11 = val5i*eqmc5
        tl12 = val6i*eqmc6
        tl13 = val7i*eqmc7
        tl14 = val8i*eqmc8
        tl15 = val9i*eqmc9
        tl16 = -coef4i*eqmc10
        tl17 = -coef4i*eqmc11
        tl18 = -coef4i*eqmc12
        tl19 = coef3i*eqmc13
        tl20 = coef3i*eqmc14
        tl21 = coef3i*eqmc15


        f(0,ix,iy,iz) = tl1 - 30.0*val2i*eqmc1 + val8*val3i*eqmc2 + 0.5*Fbar(0)

        suma = tl1 + tl2 + tl4
        sumb = tl1 + tl3 + tl5
        sumc = tl6 + coef4*tl7
        sumd = coef5*tl12 + coef4*tl13
        sume = tl8 + coef4*tl9
        sumf = -tl12 + coef5*tl13 + tl14 - coef5*tl15
        sumg = tl10 + coef4*tl11
        sumh = -tl12 + coef5*tl13 - tl14 + coef5*tl15
 
        sumi = tl12 + tl13 + tl14 + tl15
        sumk = tl12 + tl13 - tl14 - tl15
 
        sump = -coef5*tl12 - coef5*tl13

        sum67 = tl6 + tl7
        sum89 = tl8 + tl9
        sum1011 = tl10 + tl11

        f(2,ix,iy,iz) = suma + sumc + sumd + 0.5*Fbar(1)
        f(1,ix,iy,iz) = suma - sumc + sumd + 0.5*Fbar(2)
        f(4,ix,iy,iz) = suma + sume + sumf + 0.5*Fbar(3)
        f(3,ix,iy,iz) = suma - sume + sumf + 0.5*Fbar(4)
        f(6,ix,iy,iz) = suma + sumg + sumh + 0.5*Fbar(5)
        f(5,ix,iy,iz) = suma - sumg + sumh + 0.5*Fbar(6)

        f(10,ix,iy,iz) = sumb + sum67 + sum89 + sumi + tl16 + tl19 - tl20 + 0.5*Fbar(7)
        f(9,ix,iy,iz) = sumb - sum67 + sum89 + sumi - tl16 - tl19 - tl20 + 0.5*Fbar(8)
        f(8,ix,iy,iz) = sumb + sum67 - sum89 + sumi - tl16 + tl19 + tl20 + 0.5*Fbar(9)
        f(7,ix,iy,iz) = sumb - sum67 - sum89 + sumi + tl16 - tl19 + tl20 + 0.5*Fbar(10)
 
        f(14,ix,iy,iz) = sumb + sum67 + sum1011 + sumk + tl18 - tl19 + tl21 + 0.5*Fbar(11)
        f(13,ix,iy,iz) = sumb - sum67 + sum1011 + sumk - tl18 + tl19 + tl21 + 0.5*Fbar(12)
        f(12,ix,iy,iz) = sumb + sum67 - sum1011 + sumk - tl18 - tl19 - tl21 + 0.5*Fbar(13)
        f(11,ix,iy,iz) = sumb - sum67 - sum1011 + sumk + tl18 + tl19 - tl21 + 0.5*Fbar(14)
 
        f(18,ix,iy,iz) = sumb + sum89 + sum1011 + sump + tl17 + tl20 - tl21 + 0.5*Fbar(15)
        f(17,ix,iy,iz) = sumb - sum89 + sum1011 + sump - tl17 - tl20 - tl21 + 0.5*Fbar(16)
        f(16,ix,iy,iz) = sumb + sum89 - sum1011 + sump - tl17 + tl20 + tl21 + 0.5*Fbar(17)
        f(15,ix,iy,iz) = sumb - sum89 - sum1011 + sump + tl17 - tl20 + tl21 + 0.5*Fbar(18)
          
      end do
      end do
      end do

    end subroutine collision_MRT1

!===========================================================================
! This subroutine is based on element-wise operations to increase efficiency.
! Tests show that this can reduce the running time compared
! with array operations, due to the elimination of unnecessary operations
! involving the "zero" elements of the sparce matrix [M], which transfers 
! velocity space into moment space.

      subroutine collision_MRT2     
      use mpi
      use var_inc
      implicit none 

      real, dimension(0:npop-1) :: f9
      real, dimension(0:npop-1) :: Fbar
      real, dimension(0:npop-1,lx,ly)     :: tmpzmR,tmpzpR, tmpzpS, tmpzmS
      real, dimension(0:npop-1,lx,0:lz+1) :: tmpypS,tmpypR,tmpymS,tmpymR
      real rho9, ux9, uy9, uz9, ux9s, uy9s, uz9s
      real t1, tl1, tl2, tl3, tl4, tl5, tl6, tl7, tl8, tl9, tl10, tl11,&
           tl12, tl13, tl14, tl15, tl16, tl17, tl18, tl19, tl20, tl21
      real eqm1, eqm2, eqm3, eqm4, eqm5, eqm6, eqm7, eqm8, eqm9, eqm10,&
           eqm11, eqm12, eqm13, eqm14, eqm15 
      real sum1, sum2, sum3, sum4, sum5, sum6, sum7, sum8, sum9, sum10,&
           sum11
      real evlm1, evlm2, evlm3, evlm4, evlm5, evlm6, evlm7, evlm8,     &
           evlm9, evlm10, evlm11, evlm12, evlm13, evlm14, evlm15
      real eqmc1, eqmc2, eqmc3, eqmc4, eqmc5, eqmc6, eqmc7, eqmc8,     &
           eqmc9, eqmc10, eqmc11, eqmc12, eqmc13, eqmc14, eqmc15
      real suma, sumb, sumc, sumd, sume, sumf, sumg, sumh, sumi, sumk, &
           sump, sum67, sum89, sum1011
      integer ix, iy, iz, imove, jmove, kmove, sp
      integer ip, ilen
      integer status_array(MPI_STATUS_SIZE,4), req(4)
      real fx9,fy9,fz9,G1,G2,G3

      tmpzmS(:,:,:) = f(:,:,1:ly,1)
      tmpzpS(:,:,:) = f(:,:,1:ly,lz)

      ilen = npop*lx*ly
      call MPI_IRECV(tmpzmR,ilen,MPI_REAL8,mzm,0,MPI_COMM_WORLD,req(1),ierr)
      call MPI_IRECV(tmpzpR,ilen,MPI_REAL8,mzp,1,MPI_COMM_WORLD,req(2),ierr)

      call MPI_ISEND(tmpzmS,ilen,MPI_REAL8,mzm,1,MPI_COMM_WORLD,req(3),ierr)
      call MPI_ISEND(tmpzpS,ilen,MPI_REAL8,mzp,0,MPI_COMM_WORLD,req(4),ierr)
      call MPI_WAITALL(4,req,status_array,ierr)

      tmpypS(:,:,0) = tmpzmR(:,:,ly)
      tmpypS(:,:,1:lz) = f(:,:,ly,1:lz)
      tmpypS(:,:,lz+1) = tmpzpR(:,:,ly)

      tmpymS(:,:,0) = tmpzmR(:,:,1)
      tmpymS(:,:,1:lz) = f(:,:,1,1:lz)
      tmpymS(:,:,lz+1) = tmpzpR(:,:,1)

      ilen = npop * lx * (lz + 2)
      call MPI_IRECV(tmpymR,ilen,MPI_REAL8,mym,0,MPI_COMM_WORLD,req(1),ierr)
      call MPI_IRECV(tmpypR,ilen,MPI_REAL8,myp,1,MPI_COMM_WORLD,req(2),ierr)

      call MPI_ISEND(tmpymS,ilen,MPI_REAL8,mym,1,MPI_COMM_WORLD,req(3),ierr)
      call MPI_ISEND(tmpypS,ilen,MPI_REAL8,myp,0,MPI_COMM_WORLD,req(4),ierr)
      call MPI_WAITALL(4,req,status_array,ierr)

      f(:,:,0,:) = tmpymR(:,:,:)
      f(:,:,ly+1,:) = tmpypR(:,:,:)

      f(:,:,1:ly,0) = tmpzmR(:,:,:)
      f(:,:,1:ly,lz+1) = tmpzpR(:,:,:)

      do iz = 1,lz
      do iy = 1,ly
      do ix = 1,lx

        do ip = 0,npop-1
          imove = ix - cix(ip) 
          jmove = iy - ciy(ip)
          kmove = iz - ciz(ip)
          if(imove < 1 .or. imove > lx)then
            f9(ip) = f(ip,ix,iy,iz)
          else
            f9(ip) = f(ipopp(ip),imove,jmove,kmove)
          endif
        enddo
        !Determining macroscopiv variables
        sum1 = f9(7) - f9(10)
        sum2 = f9(9) - f9(8)

        sum3 = f9(11) - f9(14)
        sum4 = f9(13) - f9(12)

        sum5 = f9(15) - f9(18)
        sum6 = f9(17) - f9(16)

        ux9 = f9(1) - f9(2) + sum1 + sum2 + sum3 + sum4 + force_realx(ix,iy,iz)/2.
        uy9 = f9(3) - f9(4) + sum1 - sum2 + sum5 + sum6 + force_realy(ix,iy,iz)/2.
        uz9 = f9(5) - f9(6) + sum3 - sum4 + sum5 - sum6 + force_realz(ix,iy,iz)/2.
        ux9s = ux9*ux9
        uy9s = uy9*uy9
        uz9s = uz9*uz9

        rho9 = f9(0)+f9(1)+f9(2)+f9(3)+f9(4)+f9(5)+f9(6)&
           +f9(7)+f9(8)+f9(9)+f9(10)+f9(11)+f9(12)&
          + f9(13)+f9(14)+f9(15)+f9(16)+f9(17)+f9(18)

        !I'm merging the forcing term with the collision operator
        fx9 = force_realx(ix,iy,iz)
        fy9 = force_realy(ix,iy,iz)
        fz9 = force_realz(ix,iy,iz)
        G3 = ux9*fx9 + uy9*fy9 + uz9*fz9

        Fbar(0) = -G3

        do ip=1,6
          G1 = cix(ip)*fx9 + ciy(ip)*fy9 + ciz(ip)*fz9
          G2 = cix(ip)*ux9 + ciy(ip)*uy9 + ciz(ip)*uz9
          Fbar(ip) = ww1*(3.*G1 + 9.*G1*G2 - 3.*G3)
        enddo

        do ip=7,(npop-1)
          G1 = cix(ip)*fx9 + ciy(ip)*fy9 + ciz(ip)*fz9
          G2 = cix(ip)*ux9 + ciy(ip)*uy9 + ciz(ip)*uz9
          Fbar(ip) = ww2*(3.*G1 + 9.*G1*G2 - 3.*G3)
        enddo 

        f9(:) = f9(:) + Fbar(:)/2.

        t1 = ux9s + uy9s + uz9s
        eqm1 = -11.0*rho9 + 19.0*t1
        eqm2 = omegepsl*rho9 + omegepslj*t1
        eqm3 = coef1*ux9
        eqm4 = coef1*uy9
        eqm5 = coef1*uz9
        eqm6 = 2.0*ux9s - uy9s - uz9s
        eqm7 = omegxx*eqm6
        eqm8 = uy9s - uz9s
        eqm9 = omegxx*eqm8
        eqm10 = ux9*uy9
        eqm11 = uy9*uz9
        eqm12 = ux9*uz9
        eqm13 = 0.0
        eqm14 = 0.0
        eqm15 = 0.0

        sum1 = f9(1) + f9(2) + f9(3) + f9(4) + f9(5) + f9(6)        
        sum2 = f9(7) + f9(8) + f9(9) + f9(10) + f9(11) + f9(12)        &
             + f9(13) + f9(14) + f9(15) + f9(16) + f9(17) + f9(18)      
        sum3 = f9(7) - f9(8) + f9(9) - f9(10) + f9(11) - f9(12)        &
             + f9(13) - f9(14)
        sum4 = f9(7) + f9(8) - f9(9) - f9(10) + f9(15) - f9(16)        &
             + f9(17) - f9(18)
        sum5 = f9(11) + f9(12) - f9(13) - f9(14) + f9(15) + f9(16)     &
             - f9(17) - f9(18)
        sum6 = f9(1) + f9(2)
        sum7 = f9(3) + f9(4) + f9(5) + f9(6)
        sum8 = f9(7) + f9(8) + f9(9) + f9(10) + f9(11) + f9(12)        &
             + f9(13) + f9(14)
        sum9 = f9(15) + f9(16) + f9(17) + f9(18)
        sum10 = f9(3) + f9(4) - f9(5) - f9(6)
        sum11 = f9(7) + f9(8) + f9(9) + f9(10) - f9(11) - f9(12)       &
              - f9(13) - f9(14)
    
        evlm1 = -30.0*f9(0) + coef2*sum1 + coef3*sum2
        evlm2 = 12.0*f9(0) + coef4*sum1 + sum2
        evlm3 = coef4*(f9(1) - f9(2)) + sum3
        evlm4 = coef4*(f9(3) - f9(4)) + sum4
        evlm5 = coef4*(f9(5) - f9(6)) + sum5
        evlm6 = coef5*sum6 - sum7 + sum8 - coef5*sum9
        evlm7 = coef4*sum6 + coef5*sum7 + sum8 - coef5*sum9
        evlm8 = sum10 + sum11
        evlm9 =-coef5*sum10 + sum11
        evlm10 = f9(7) - f9(8) - f9(9) + f9(10)
        evlm11 = f9(15) - f9(16) - f9(17) + f9(18)
        evlm12 = f9(11) - f9(12) - f9(13) + f9(14)
        evlm13 = f9(7) - f9(8) + f9(9) - f9(10) - f9(11) + f9(12)      &
               - f9(13) + f9(14)
        evlm14 =-f9(7) - f9(8) + f9(9) + f9(10) + f9(15) - f9(16)      &
               + f9(17) - f9(18)
        evlm15 = f9(11) + f9(12) - f9(13) - f9(14) - f9(15) - f9(16)   &
               + f9(17) + f9(18)
   
        eqmc1 = evlm1 - s1*(evlm1 - eqm1)
        eqmc2 = evlm2 - s2*(evlm2 - eqm2)
        eqmc3 = evlm3 - s4*(evlm3 - eqm3)
        eqmc4 = evlm4 - s4*(evlm4 - eqm4)
        eqmc5 = evlm5 - s4*(evlm5 - eqm5)
        eqmc6 = evlm6 - s9*(evlm6 - eqm6)
        eqmc7 = evlm7 - s10*(evlm7 - eqm7)
        eqmc8 = evlm8 - s9*(evlm8 - eqm8)
        eqmc9 = evlm9 - s10*(evlm9 - eqm9)
        eqmc10 = evlm10 - s13*(evlm10 - eqm10)
        eqmc11 = evlm11 - s13*(evlm11 - eqm11)
        eqmc12 = evlm12 - s13*(evlm12 - eqm12)
        eqmc13 = evlm13 - s16*(evlm13 - eqm13)
        eqmc14 = evlm14 - s16*(evlm14 - eqm14)
        eqmc15 = evlm15 - s16*(evlm15 - eqm15)


        tl1 = val1i*rho9
        tl2 = coef2*val2i*eqmc1
        tl3 = coef3*val2i*eqmc1
        tl4 = coef4*val3i*eqmc2
        tl5 = val3i*eqmc2
        tl6 = val4i*ux9
        tl7 = val5i*eqmc3
        tl8 = val4i*uy9
        tl9 = val5i*eqmc4
        tl10 = val4i*uz9
        tl11 = val5i*eqmc5
        tl12 = val6i*eqmc6
        tl13 = val7i*eqmc7
        tl14 = val8i*eqmc8
        tl15 = val9i*eqmc9
        tl16 = -coef4i*eqmc10
        tl17 = -coef4i*eqmc11
        tl18 = -coef4i*eqmc12
        tl19 = coef3i*eqmc13
        tl20 = coef3i*eqmc14
        tl21 = coef3i*eqmc15


        f9(0) = tl1 - 30.0*val2i*eqmc1 + val8*val3i*eqmc2

        suma = tl1 + tl2 + tl4
        sumb = tl1 + tl3 + tl5
        sumc = tl6 + coef4*tl7
        sumd = coef5*tl12 + coef4*tl13
        sume = tl8 + coef4*tl9
        sumf = -tl12 + coef5*tl13 + tl14 - coef5*tl15
        sumg = tl10 + coef4*tl11
        sumh = -tl12 + coef5*tl13 - tl14 + coef5*tl15
 
        sumi = tl12 + tl13 + tl14 + tl15
        sumk = tl12 + tl13 - tl14 - tl15
 
        sump = -coef5*tl12 - coef5*tl13

        sum67 = tl6 + tl7
        sum89 = tl8 + tl9
        sum1011 = tl10 + tl11

        f9(1) = suma + sumc + sumd
        f9(2) = suma - sumc + sumd
        f9(3) = suma + sume + sumf
        f9(4) = suma - sume + sumf
        f9(5) = suma + sumg + sumh
        f9(6) = suma - sumg + sumh

        f9(7) = sumb + sum67 + sum89 + sumi + tl16 + tl19 - tl20
        f9(8) = sumb - sum67 + sum89 + sumi - tl16 - tl19 - tl20
        f9(9) = sumb + sum67 - sum89 + sumi - tl16 + tl19 + tl20
        f9(10) = sumb - sum67 - sum89 + sumi + tl16 - tl19 + tl20
 
        f9(11) = sumb + sum67 + sum1011 + sumk + tl18 - tl19 + tl21
        f9(12) = sumb - sum67 + sum1011 + sumk - tl18 + tl19 + tl21
        f9(13) = sumb + sum67 - sum1011 + sumk - tl18 - tl19 - tl21
        f9(14) = sumb - sum67 - sum1011 + sumk + tl18 + tl19 - tl21
 
        f9(15) = sumb + sum89 + sum1011 + sump + tl17 + tl20 - tl21
        f9(16) = sumb - sum89 + sum1011 + sump - tl17 - tl20 - tl21
        f9(17) = sumb + sum89 - sum1011 + sump - tl17 + tl20 + tl21
        f9(18) = sumb - sum89 - sum1011 + sump + tl17 - tl20 + tl21

        do ip = 0,npop-1
          imove = ix + cix(ip) 
          jmove = iy + ciy(ip)
          kmove = iz + ciz(ip)
          if(imove < 1 .or. imove > lx)then
            f(ipopp(ip),ix,iy,iz) = f9(ip) + 0.5*Fbar(ip)
          else
            f(ip,imove,jmove,kmove) = f9(ip) + 0.5*Fbar(ip)
          endif
        enddo
          
      end do
      end do
      end do

      tmpymS(:,:,:) = f(:,:,0,:)
      tmpypS(:,:,:) =  f(:,:,ly+1,:)

      tmpzmS(:,:,:) = f(:,:,1:ly,0)
      tmpzpS(:,:,:) = f(:,:,1:ly,lz+1)

      !The data is sent in reverse order, y first then z
      ilen = npop * lx * (lz + 2)
      call MPI_IRECV(tmpymR,ilen,MPI_REAL8,mym,0,MPI_COMM_WORLD,req(1),ierr)
      call MPI_IRECV(tmpypR,ilen,MPI_REAL8,myp,1,MPI_COMM_WORLD,req(2),ierr)

      call MPI_ISEND(tmpymS,ilen,MPI_REAL8,mym,1,MPI_COMM_WORLD,req(3),ierr)
      call MPI_ISEND(tmpypS,ilen,MPI_REAL8,myp,0,MPI_COMM_WORLD,req(4),ierr)
      call MPI_WAITALL(4,req,status_array,ierr)

      !Attach diagonal moving data to the z shifts
      tmpzpS(15,:,1) = tmpymR(15,:,lz+1)
      tmpzpS(16,:,ly) = tmpypR(16,:,lz+1)

      tmpzmS(17,:,1) = tmpymR(17,:,0)
      tmpzmS(18,:,ly) = tmpypR(18,:,0)

      ilen = npop * lx * ly
      call MPI_IRECV(tmpzmR,ilen,MPI_REAL8,mzm,0,MPI_COMM_WORLD,req(1),ierr)
      call MPI_IRECV(tmpzpR,ilen,MPI_REAL8,mzp,1,MPI_COMM_WORLD,req(2),ierr)

      call MPI_ISEND(tmpzmS,ilen,MPI_REAL8,mzm,1,MPI_COMM_WORLD,req(3),ierr)
      call MPI_ISEND(tmpzpS,ilen,MPI_REAL8,mzp,0,MPI_COMM_WORLD,req(4),ierr)

      call MPI_WAITALL(4,req,status_array,ierr)

      f(4,:,ly,1:lz) = tmpypR(4,:,1:lz)
      f(9,2:lx,ly,1:lz) = tmpypR(9,2:lx,1:lz)
      f(10,1:lx-1,ly,1:lz) = tmpypR(10,1:lx-1,1:lz)
      f(16,:,ly,1:lz) = tmpypR(16,:,1:lz) !Note this is incorrect for y/z edge
      f(18,:,ly,1:lz) = tmpypR(18,:,1:lz) !Note this is incorrect for y/z edge

      f(3,:,1,1:lz) = tmpymR(3,:,1:lz)
      f(7,2:lx,1,1:lz) = tmpymR(7,2:lx,1:lz)
      f(8,1:lx-1,1,1:lz) = tmpymR(8,1:lx-1,1:lz)
      f(15,:,1,1:lz) = tmpymR(15,:,1:lz) !Note this is incorrect for y/z edge
      f(17,:,1,1:lz) = tmpymR(17,:,1:lz) !Note this is incorrect for y/z edge
      
      f(6,:,1:ly,lz) = tmpzpR(6,:,:)
      f(13,2:lx,1:ly,lz) = tmpzpR(13,2:lx,:)
      f(14,1:lx-1,1:ly,lz) = tmpzpR(14,1:lx-1,:)
      f(17,:,1:ly,lz) = tmpzpR(17,:,:)
      f(18,:,1:ly,lz) = tmpzpR(18,:,:)

      f(5,:,1:ly,1) = tmpzmR(5,:,:)
      f(11,2:lx,1:ly,1) = tmpzmR(11,2:lx,:)
      f(12,1:lx-1,1:ly,1) = tmpzmR(12,1:lx-1,:)
      f(15,:,1:ly,1) = tmpzmR(15,:,:)
      f(16,:,1:ly,1) = tmpzmR(16,:,:)

      end subroutine collision_MRT2
!==========================================================================

      subroutine rhoupdat 
      use var_inc
      implicit none 

!      rho = f(0,:,:,:)+f(1,:,:,:)+f(2,:,:,:)+f(3,:,:,:)+f(4,:,:,:)     &
!          +f(5,:,:,:)+f(6,:,:,:)+f(7,:,:,:)+f(8,:,:,:)+f(9,:,:,:)     &
!         +f(10,:,:,:)+f(11,:,:,:)+f(12,:,:,:)+f(13,:,:,:)+f(14,:,:,:)&
!        +f(15,:,:,:)+f(16,:,:,:)+f(17,:,:,:)+f(18,:,:,:)

      integer ip

      rho = f(0,:,:,:)
      do ip = 1,npop-1
        rho = rho + f(ip,:,:,:)
      end do
      end subroutine rhoupdat 
!===================================================================

      subroutine macrovar 
      use var_inc
      implicit none 

      integer ip, i, j, k, id             
      real xc, yc, zc, xpnt, ypnt, zpnt, xx0, yy0, zz0     
      real w1, w2, w3, omg1, omg2, omg3 


      integer  iz,iy,ix
      real  sum1,sum2,sum3,sum4,sum5,sum6,ux9,uy9,uz9
      real  rho9
      real, dimension(0:npop-1) :: f9

      do iz = 1,lz
      do iy = 1,ly
      do ix = 1,lx
      f9 = f(:,ix,iy,iz)

      sum1 = f9(7) - f9(10)
      sum2 = f9(9) - f9(8)

      sum3 = f9(11) - f9(14)
      sum4 = f9(13) - f9(12)

      sum5 = f9(15) - f9(18)
      sum6 = f9(17) - f9(16)

      ux9 = f9(1) - f9(2) + sum1 + sum2 + sum3 + sum4
      uy9 = f9(3) - f9(4) + sum1 - sum2 + sum5 + sum6
      uz9 = f9(5) - f9(6) + sum3 - sum4 + sum5 - sum6

      rho9 = f9(0)+f9(1)+f9(2)+f9(3)+f9(4)+f9(5)+f9(6)&
           +f9(7)+f9(8)+f9(9)+f9(10)+f9(11)+f9(12)&
          + f9(13)+f9(14)+f9(15)+f9(16)+f9(17)+f9(18)

      ux(ix,iy,iz) = ux9 + force_realx(ix,iy,iz)/2.
      uy(ix,iy,iz) = uy9 + force_realy(ix,iy,iz)/2.
      uz(ix,iy,iz) = uz9 + force_realz(ix,iy,iz)/2.
      rho(ix,iy,iz) = rho9

      enddo
      enddo
      enddo
      
! for solid particle nodes
      if(ipart)then

      do k = 1,lz
      do j = 1,ly
      do i = 1,lx
      if(ibnodes(i,j,k) > 0)then
        id = isnodes(i,j,k) 

        xpnt = real(i) - 0.5
        ypnt = real(j) - 0.5 + real(indy*ly)
        zpnt = real(k) - 0.5 + real(indz*lz)

        xc = ypglb(1,id)
        yc = ypglb(2,id)
        zc = ypglb(3,id)

! use the nearest particle center instead of the real center
!       if((xc - xpnt) > real(nxh)) xc = xc - real(nx)
!       if((xc - xpnt) < -real(nxh)) xc = xc + real(nx)
!        if(abs(xc -xpnt) <= real(nxh)) xc = xc

        if((yc - ypnt) > real(nyh)) yc = yc - real(ny)
        if((yc - ypnt) < -real(nyh)) yc = yc + real(ny)
!        if(abs(yc -ypnt) <= real(nyh)) yc = yc

        if((zc - zpnt) > real(nzh)) zc = zc - real(nz)
        if((zc - zpnt) < -real(nzh)) zc = zc + real(nz)
!        if(abs(zc -zpnt) <= real(nzh)) zc = zc

        xx0 = xpnt - xc
        yy0 = ypnt - yc
        zz0 = zpnt - zc

        w1 = wp(1,id)
        w2 = wp(2,id)
        w3 = wp(3,id)

        omg1 = omgp(1,id)
        omg2 = omgp(2,id)
        omg3 = omgp(3,id)

        ux(i,j,k) = w1 + (omg2*zz0 - omg3*yy0)
        uy(i,j,k) = w2 + (omg3*xx0 - omg1*zz0)
        uz(i,j,k) = w3 + (omg1*yy0 - omg2*xx0)

        rho(i,j,k) = rhopart
      end if
      end do
      end do
      end do

      end if

     
      end subroutine macrovar
!===================================================================
!==================================================================
      SUBROUTINE FORCING
      use mpi
      use var_inc
      implicit none
      integer ixs,ihh,i,j,k,jj,kk
      real x9,y9,z9

      force_realx(:,:,:) = 0.0
      force_realy(:,:,:) = force_in_y*force_mag
      force_realz(:,:,:) = 0.0

      RETURN
      END SUBROUTINE FORCING
!=================================
      SUBROUTINE FORCINGP
      use mpi
      use var_inc
      implicit none
      integer ixs,ihh,i,j,k,jj,kk
      real x9,y9,z9,Amp0,beta9,gamma9,Tpd, phase9,Tpdp,ixs0

      Tpd = 2000.
      Tpdp = 1500.
      beta9 = 3.0
      gamma9=2.0
      ixs0 = 2
      Amp0 = 40.00*beta9/real(ny)*sin(pi2*real(istep)/Tpd)
!     phase9 = sin(pi2*real(istep)/Tpdp)
      phase9 = 0.25

      force_realx(:,:,:) = 0.0
      force_realy(:,:,:) = force_in_y
      force_realz(:,:,:) = 0.0
! add some perturbation
      ixs = ixs0
      ihh = lxh/2
      do k=1,lz
      kk = k + indz*lz
      z9 = pi2*(real(kk)-0.5)/real(nz)
      do j=1,ly
      jj = j + indy*ly
      y9 = pi2*(real(jj)-0.5)/real(ny)
      do i=1,ihh
      x9 = pi2*(real(i)-0.5)/real(ihh)
      force_realx(i+ixs,j,k) = force_in_y*0.5*Amp0*real(ihh)*(1.-cos(x9))  &
                                 *cos(beta9*y9)*cos(gamma9*z9)
      force_realy(i+ixs,j,k) = force_in_y*(1.0-Amp0*real(ny)/beta9    &
                       *sin(x9)*sin(beta9*y9)*cos(gamma9*z9))
      force_realz(i+ixs,j,k) = force_in_y*0.5*Amp0*real(nz)/gamma9*   &
                        sin(x9)*cos(beta9*y9)*sin(gamma9*z9)
      end do
      end do
      end do

      ihh = lxh/2
      ixs = nx - ixs0 - ihh
      do k=1,lz
      kk = k + indz*lz
      z9 = pi2*(real(kk)-0.5)/real(nz)
      do j=1,ly
      jj = j + indy*ly
      y9 = pi2*( (real(jj)-0.5)/real(ny) + phase9 )
      do i=1,ihh
      x9 = pi2*(real(i)-0.5)/real(ihh)
      force_realx(i+ixs,j,k) = -force_in_y*0.5*Amp0*real(ihh)*(1.-cos(x9))  &
                                 *cos(beta9*y9)*cos(gamma9*z9)
      force_realy(i+ixs,j,k) = force_in_y*(1.0+Amp0*real(ny)/beta9    &
                       *sin(x9)*sin(beta9*y9)*cos(gamma9*z9))
      force_realz(i+ixs,j,k) = -force_in_y*0.5*Amp0*real(nz)/gamma9*   &
                        sin(x9)*cos(beta9*y9)*sin(gamma9*z9)
      end do
      end do
      end do

      RETURN
      END SUBROUTINE FORCINGP

