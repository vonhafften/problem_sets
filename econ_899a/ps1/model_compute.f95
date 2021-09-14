! Modified by Alex von Hafften on Sept 14, 2021
! Based on script written by Pablo D'Erasmo   6 Sep 2007    5:56 pm
! Compile with 'gfortran -o model_compute model_compute.f95'
! Run with './model_compute'

module parameters
   implicit none

   REAL,    PARAMETER :: b = 0.99, d = 0.025, a = 0.36
   REAL,    PARAMETER :: z_low = 0.2, z_high = 1.25
   REAL,    PARAMETER :: pi_h_h = 0.977, pi_h_l = 0.023, pi_l_h = 0.074, pi_l_l = 0.926
   REAL,    PARAMETER :: klb = 0.01, inc = 0.025, kub = 45.0
   INTEGER, PARAMETER :: length_grid_k = (kub-klb)/inc + 1
   REAL ,   PARAMETER :: toler   = 1.e-4  ! Numerical tolerance

end module

! ============================================================================================

module global
   USE parameters
   implicit none
   
   REAL  :: Kgrid(length_grid_k)
   REAL  :: value_low(length_grid_k), g_k_low(length_grid_k)
   REAL  :: vtmp_low(length_grid_k,length_grid_k), value_new_low(length_grid_k)
   REAL  :: value_high(length_grid_k), g_k_high(length_grid_k)
   REAL  :: vtmp_high(length_grid_k,length_grid_k), value_new_high(length_grid_k)

end module

! ============================================================================================

PROGRAM  compute_model
   REAL               :: total, etime, dist
   REAL, DIMENSION(2) :: elapsed

   call solution

   total=etime(elapsed)

   PRINT*, '--------------------------------------------------'
   PRINT*, 'total time elpased =', total
   PRINT*, '--------------------------------------------------'

END PROGRAM compute_model

! ============================================================================================

subroutine solution

   USE parameters
   USE global

   IMPLICIT  NONE

   INTEGER :: iter, i_k, i_kp
   REAL    :: diff, diff_low, diff_high, k, kp, c_low, c_high
   INTEGER :: i = 1

   do while (i <= length_grid_k)   !do loop for assigning capital grid K
     Kgrid(i) = klb + (i - 1) * inc
     i = i + 1
   end do


   iter = 1
   diff = 1000.d0
   value_low = 0.0 * Kgrid   ! Initial Value guess
   value_high = 0.0 * Kgrid   ! Initial Value guess
   
	do while (diff >= toler)

		do i_k = 1, length_grid_k     ! Capital grid

         k = Kgrid(i_k)
         vtmp_low(i_k,:) = -1.0e-16
         vtmp_high(i_k,:) = -1.0e-16

			do i_kp = 1, length_grid_k
            kp = Kgrid(i_kp)

            c_low = z_low * k ** a + (1.0 - d) * k - kp
            c_high = z_high * k ** a + (1.0 - d) * k - kp

				if (c_low > 0.0) then
               vtmp_low(i_k, i_kp) = log(c_low) + b * value_low(i_kp) * pi_l_l + b * value_high(i_kp) * pi_l_h
               vtmp_high(i_k, i_kp) = log(c_high) + b * value_low(i_kp) * pi_h_l + b * value_high(i_kp) * pi_h_h
	         endif

			enddo

         value_new_low(i_k) = MAXVAL(vtmp_low(i_k,:))
         value_new_high(i_k) = MAXVAL(vtmp_high(i_k,:))

         g_k_low(i_k) = Kgrid(MAXLOC(vtmp_low(i_k,:),1))
         g_k_high(i_k) = Kgrid(MAXLOC(vtmp_high(i_k,:),1))

      enddo

      diff_low  = maxval(abs(value_new_low-value_low))/ABS(value_new_low(length_grid_k))
      diff_high  = maxval(abs(value_new_high-value_high))/ABS(value_new_high(length_grid_k))

      diff = MAX(diff_low, diff_high)

      value_low = value_new_low
      value_high = value_new_high

      print*, 'Iteration =',iter,'sup_norm =', diff
      iter = iter+1

	enddo

   print *, ' '
   print *, 'Successfully converged with sup_norm ', diff

   open (UNIT=1,FILE='fortran_output.csv',STATUS='replace')

   do i_k = 1, length_grid_k
      WRITE(UNIT=1,FMT=*) Kgrid(i_k), ",", value_low(i_k), ",", value_high(i_k), ",", g_k_low(i_k), ",", g_k_high(i_k)
   end do

   close (UNIT=1)

end subroutine
