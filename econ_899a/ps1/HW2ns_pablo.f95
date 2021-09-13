!     Last change:  PND   6 Sep 2007    5:56 pm

module parameters
implicit none
   REAL, PARAMETER  	:: b = 0.99, d = 0.025, a = 0.36
   REAL, PARAMETER  	:: klb = 0.01, inc = 0.025, kub = 45.0
   INTEGER, PARAMETER 	:: length_grid_k = (kub-klb)/inc + 1
   REAL , PARAMETER 	:: toler   = 1.e-4						! Numerical tolerance
end module

! ============================================================================================
module global
USE parameters
implicit none
   REAL 		:: Kgrid(length_grid_k), value(length_grid_k), g_k(length_grid_k)
   REAL                 :: vtmp(length_grid_k,length_grid_k), value_new(length_grid_k)
end module

! ============================================================================================

PROGRAM  HW2NonStochastic
   REAL 			:: total, etime, dist
   REAL, DIMENSION(2)  		:: elapsed

   call solution

   total=etime(elapsed)


	PRINT*,'--------------------------------------------------'
	PRINT*,'total time elpased =',total
	PRINT*,'--------------------------------------------------'


END PROGRAM HW2NonStochastic

! ============================================================================================
subroutine solution
USE parameters
USE global

   IMPLICIT  NONE

   INTEGER :: iter, index_k, index_kp
   REAL :: diff, k, kp, c
   
   INTEGER :: i = 1


   do while (i<=length_grid_k)   !do loop for assigning capital grid K
     Kgrid(i) = klb + (i-1)*inc
     !write(*,*) i, Kgrid(i)
     i = i + 1
   end do


   iter = 1
   diff = 1000.d0
   value = 0.*Kgrid		!Initial Value guess
   
	do while (diff>= toler)

		do index_k = 1, length_grid_k				! Capital grid
			k = Kgrid(index_k)
                        vtmp(index_k,:) = -1.0e-16

			do index_kp = 1, length_grid_k
				kp = Kgrid(index_kp)
				c = k**a+(1.-d)*k-kp

				if (c>0.) then
	                        	vtmp(index_k,index_kp) = log(c)+b*value(index_kp)
	          		endif

			enddo

			value_new(index_k) = MAXVAL(vtmp(index_k,:))
                        g_k(index_k) 	   = Kgrid(MAXLOC(vtmp(index_k,:),1))
                enddo

		diff  = maxval(abs(value_new-value))/ABS(value_new(length_grid_k))
		value = value_new


		print*, 'Iteration =',iter,'sup_norm =',diff
                iter = iter+1

	enddo

	print *, ' '
	print *, 'Successfully converged with sup_norm ', diff
    	!print *, g_k
    
    !CALL vcDrawCurve@(d, Kgrid, g_k, length_grid_k)


        open (UNIT=1,FILE='valuefun',STATUS='replace')
	do index_k = 1, length_grid_k
        	WRITE(UNIT=1,FMT=*) value(index_k)
        end do
	close (UNIT=1)

end subroutine
