module Signal
    use JobConstants
    use ShockSignal, only: applyShock
    implicit none
    
    public:: nonDimensionalizeSignal
    public:: amplifySignal
    public:: shearSignal
    public:: spaceSignal
	public:: fineSignal
	public:: coarseSignal
    private:: calculateGamma
    private:: applyShock
    
    contains
!*****************************************************************************
subroutine fineSignal(nt, time_input, pressure_amplified, nt_fine, time_input_fine, pressure_amplified_fine)
	integer,intent(in):: nt,nt_fine
	real,dimension(nt),intent(in):: time_input, pressure_amplified
	real,dimension(nt_fine),intent(out):: time_input_fine, pressure_amplified_fine
	
	integer:: i,j
	
	j = 1
	do i = 1,nt_fine
		if (mod(i,2) .ne. 0) then
			time_input_fine(i) = time_input(j)
			pressure_amplified_fine(i) = pressure_amplified(j)
			j = j + 1
		else if (mod(i,2) .eq. 0) then
			time_input_fine(i) = (time_input(j) + time_input(j - 1)) / 2
			pressure_amplified_fine(i) = (pressure_amplified(j) + pressure_amplified(j - 1)) / 2
		end if
	end do
	
end subroutine fineSignal
!*****************************************************************************
subroutine nonDimensionalizeSignal(nt, pressure_amplified, pressure_amplified_nondim)
    integer,intent(in):: nt
    real,dimension(nt),intent(in):: pressure_amplified
    real,dimension(nt),intent(out):: pressure_amplified_nondim
    
    real:: minimum
    
    minimum = abs(minval(pressure_amplified))
    
    pressure_amplified_nondim = pressure_amplified / minimum
end subroutine nonDimensionalizeSignal
!*****************************************************************************
subroutine shearSignal(nt, mach, chord, pressure_amplified, pressure_amplified_nondim, time_input, time_sheared)
    integer,intent(in):: nt
    real,intent(in):: mach,chord
    real,dimension(nt),intent(in):: time_input, pressure_amplified_nondim
    real,dimension(nt),intent(out):: time_sheared 
    real,dimension(nt),intent(inout):: pressure_amplified
    
    integer:: i, j, time_maximum_index, time_minimum_index
	integer:: time_shock_top_index, time_shock_bottom_index
	integer:: max_shock_index, min_shock_index, mid_line_index
    integer,dimension(1):: min_line_index, max_line_index
    real:: gamma_angle, tan_gamma
    real:: delta_time, time_maximum, time_minimum, time_average
	real:: time_shock_top, time_current, time_shock_bottom
	real:: difference, previous_difference, time_shift
    real,dimension(nt):: time_widened
    logical:: time_max, time_min, time_half
    
    call calculateGamma(nt,mach,chord,gamma_angle,time_input,pressure_amplified_nondim)
    
    time_max = .false.
    time_min = .false.
    time_half = .true.
    delta_time = pressure_amplified_nondim(1) * tan(gamma_angle)
    tan_gamma = tan(gamma_angle)
    time_sheared(1) = time_input(1) - delta_time
    time_average = 1
    
    ! This loop shears and widens the signal
    do i = 2,nt
    
        ! Shears the signal
        delta_time = pressure_amplified_nondim(i) * tan_gamma
        time_sheared(i) = time_input(i) - delta_time
        
        ! If the signal begins curving back on itself
        if(time_sheared(i) .lt. time_sheared(i - 1) .and. time_max .eq. .false.) then
        
            ! Point where signal curves on itself
            time_max = .true.
            time_maximum_index = i - 1
            time_maximum = time_sheared(i - 1)
            
        end if
        
        ! Check if time value is time_minimum if signal begins to curve back on itself
        if(time_sheared(i) .gt. time_sheared(i - 1) .and. time_max .eq. .true.) then
        
            time_max = .false.
            time_minimum_index = i - 1
            time_minimum = time_sheared(i - 1)
            
            time_min = .true.
            
!            time_half = .true.
!            time_average = (time_minimum + time_maximum) / 2

        end if
        
        if(time_sheared(i) - time_maximum .ge. 0 .and. time_min .eq. .true.) then
            
            max_shock_index = i - 1
            
            time_min = .false.
            
            j = time_maximum_index 
            
            do while(time_sheared(j) - time_minimum .gt. 0)
            
                j = j - 1
            
            end do
            
            min_shock_index = j + 1
            
            max_line_index = minloc(pressure_amplified(min_shock_index:time_maximum_index))
            max_line_index = max_line_index + min_shock_index
            min_line_index = maxloc(pressure_amplified(time_minimum_index:max_shock_index))
            min_line_index = min_line_index + time_minimum_index
            mid_line_index = int((max_line_index(1) + min_line_index(1)) / 2)
            
            call applyShock(min_line_index(1), max_line_index(1), min_shock_index, max_shock_index, nt, time_sheared, pressure_amplified)
            
        end if
        
!        ! Find the time where the shock ends
!        difference = time_average - time_sheared(i)
!        previous_difference = time_average - time_sheared(i - 1)
!        if(abs(previous_difference) .lt. abs(difference) .and. time_half .eqv. .true.) then
!        
!            time_half = .false.
!            time_shock_top_index = i - 1
!            time_shock_top = time_sheared(i - 1)
!            time_current = time_sheared(time_maximum_index)
!            j = time_maximum_index
!            
!            ! Find the time where the shock starts
!            difference = time_average - time_current
!            previous_difference = time_average - time_sheared(j - 1)
!            do while(abs(previous_difference) .lt. abs(difference))
!            
!                j = j - 1
!                time_current = time_sheared(j)
!                difference = time_average - time_current
!                previous_difference = time_average - time_sheared(j - 1)
!                
!            end do
!            
!            time_shock_bottom_index = j + 1
!            time_shock_bottom = time_sheared(j + 1)
!            ! Call subroutine to insert shocks where signal loops back on itself
!            ! ERROR: IF signal starts or ends in middle of signal, then issue with inserting shock
!            call applyShock(nt,time_sheared,time_shock_bottom_index,time_shock_top_index,pressure_amplified)
!            
!        end if
        
        
        
    end do
    
end subroutine shearSignal
!*****************************************************************************
subroutine amplifySignal(nt,mach,pressure_amplified,pressure_amplified_amplified)
    integer,intent(in):: nt
    real,intent(in):: mach
    real,dimension(nt),intent(in):: pressure_amplified
    real,dimension(nt),intent(out):: pressure_amplified_amplified
    
    real:: mach_factor,HSI_factor
	
	
	! ERROR: Insert failsafe if Mach = 1 or greater
    mach_factor = 1 / (1 - mach)

    if (mach_factor .lt. 2.5) then
        HSI_factor = 1
    else if (mach_factor .ge. 2.5 .and. mach_factor .lt. 3.33) then
        HSI_factor = 0.1287 * mach + 1.1316
    else if (mach_factor .ge. 3.33 .and. mach_factor .lt. 5) then
        HSI_factor = 0.7321 * mach + 0.7092
    ! Baeder's
    else if (mach_factor .ge. 5 .and. mach_factor .lt. 10) then
        HSI_factor = 87.225 * (mach ** 2) - 140.52 * mach + 57.89
    ! Plots
!    elseif (mach_factor .ge. 5 .and. mach_factor .lt. 10) then
!        HSI_factor = 172.47 * (mach ** 2) - 287.49 * mach + 121.35
    else if (mach_factor .ge. 10 .and. mach_factor .lt. 13.33) then
        HSI_factor = 2.7613 * mach - 0.3968
    else if (mach_factor .ge. 13.33) then
        HSI_factor = 2646.1 * (mach ** 3) - 7593.1 * (mach ** 2) + 7238.2 * mach - 2290.6
        ! HSI_factor = 1.085336538 * HSI_factor
    else
        HSI_factor = 1
    end if
    
    pressure_amplified_amplified = pressure_amplified * HSI_factor
end subroutine amplifySignal
!*****************************************************************************
subroutine coarseSignal(nt, time, pressure, nt_fine, time_sheared, pressure_amplified_amplified)
	integer,intent(in):: nt, nt_fine
	real,dimension(nt_fine),intent(in):: time_sheared, pressure_amplified_amplified
	real,dimension(nt),intent(out):: time, pressure
	
	integer:: i,j
	
	j = 1
	do i = 1,nt_fine
		if (mod(i,2) .ne. 0) then	
			time(j) = time_sheared(i)
			pressure(j) = pressure_amplified_amplified(i)
			j = j + 1
		end if
	end do
	
end subroutine coarseSignal
!*****************************************************************************
subroutine calculateGamma(nt,mach,chord,gamma_angle,time,pressure)
    integer,intent(in):: nt
    real,intent(in):: mach,chord
    real,intent(out):: gamma_angle
    real,dimension(nt),intent(in):: time,pressure
    
    integer:: i,j,number_blades
    integer,dimension(1):: minimum_pressure_index,maximum_pressure_index
    real:: percent_shift,check,time_minimum,time_maximum,delta_time,time_shift
    
    if (mach .lt. 0.6) then
        percent_shift = 0.0
    elseif (mach .ge. 0.6 .and. mach .lt. 0.8) then
        percent_shift = -0.001 * (mach ** 2) + 0.0018 * mach - 0.008 
    elseif (mach .ge. 0.8 .and. mach .lt. 0.85) then
        percent_shift = 0.0004 * mach - 0.0003
    elseif (mach .ge. 0.85 .and. mach .lt. 0.88) then
        percent_shift = 0.0013 * mach - 0.0011
    elseif (mach .ge. 0.88 .and. mach .lt. 0.9) then
        percent_shift = 0.0035 * mach - 0.003
    elseif (mach .ge. 0.9 .and. mach .lt. 1.0) then
        percent_shift = 0.0048 * mach - 0.0042
    elseif (mach .gt. 1.0) then
        percent_shift = 0
    endif
    
!    i = 0
!    j = 0
!    
!    minimum_pressure_index = minloc(pressure(1:int(nt/blades)))
!    maximum_pressure_index = maxloc(pressure(1:int(nt/blades)))
!    
!    time_maximum = time(minimum_pressure_index(1))
!    time_minimum = time(maximum_pressure_index(1))
!    delta_time = time_maximum - time_minimum
    
!    time_shift = percent_shift * delta_time
    
	! 0.25ft is the chord length of Baeder's UH-1H blade 
    gamma_angle = percent_shift * chord / 0.25 
end subroutine
!*****************************************************************************
!subroutine applyShock(nt,time_sheared,time_shock_bottom_index,time_shock_top_index,pressure_thickness_amplified)
!    integer,intent(in):: nt,time_shock_bottom_index,time_shock_top_index
!    real,dimension(nt),intent(inout):: time_sheared,pressure_thickness_amplified
!    
!    integer:: i
!    real:: slope_time,y_int_time,slope_pressure,y_int_pressure
!    real:: pressure_top,pressure_bottom
!    
!    slope_time = (time_sheared(time_shock_top_index) - time_sheared(time_shock_bottom_index)) / (time_shock_top_index - time_shock_bottom_index) 
!    y_int_time = time_sheared(time_shock_top_index) - slope_time * time_shock_top_index
!    pressure_top = pressure_thickness_amplified(time_shock_top_index)
!    pressure_bottom = pressure_thickness_amplified(time_shock_bottom_index)
!    slope_pressure = (pressure_top - pressure_bottom) / (time_shock_top_index - time_shock_bottom_index)
!    y_int_pressure = pressure_top - slope_pressure * time_shock_top_index
!    
!    do i = time_shock_bottom_index,time_shock_top_index
!    
!        time_sheared(i) = slope_time * i + y_int_time
!        pressure_thickness_amplified(i) = slope_pressure * i + y_int_pressure
!        
!    end do
!    
!end subroutine applyShock
!*****************************************************************************
subroutine spaceSignal(n,y1,y2,x2,x1) 
    integer,intent(in):: n
    real,dimension(n),intent(in):: y1,x2,x1
    real,dimension(n),intent(out):: y2
    
    integer:: i,j
    real:: diff,prev_diff
    logical:: array_done
    
    j = 2
    array_done = .false.
    
    do i = 1,n
        
        diff = abs(x2(i) - x1(j))
        prev_diff = abs(x2(i) - x1(j - 1))
        
        do while(diff .lt. prev_diff .and. array_done .eq. .false.)
            
            j = j + 1
            diff = abs(x2(i) - x1(j))
            prev_diff = abs(x2(i) - x1(j - 1))
            
            if(j .eq. n) then
                array_done = .true.
            end if
            
        end do
        
        if(array_done .eq. .false. .and. j .ne. 2) then
            j = j - 1
        else if(array_done .eq. .true.) then
            j = j
        end if
        
        y2(i) = y1(j) - (((y1(j) - y1(j - 1)) * (x1(j) - x2(i))) / (x1(j) - x1(j - 1)))
         
    end do
    
end subroutine spaceSignal
!*****************************************************************************
end module
