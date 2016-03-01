module ShockSignal
    implicit none
    public:: applyShock
    private:: locateLineIndex
    private:: integrate
    private:: shockLine
    contains
!*****************************************************************************
subroutine applyShock(max_press_index,min_press_index,min_index,max_index,nt,time_sheared,pressure_thickness)
    integer,intent(in):: max_press_index,min_press_index,min_index,max_index,nt
    real,dimension(nt),intent(inout):: time_sheared,pressure_thickness
    
    integer:: lower_shock_line_index,upper_shock_line_index,mid_shock_line_index
    real:: area_lower,area_lower_1,area_lower_2,area_upper,area_upper_1,area_upper_2,net_area,net_area_previous
    logical:: stop_loop
    
    ! Draw first shock line adjacent to maximum pressure index
    upper_shock_line_index = int((max_index + max_press_index) / 2)
    stop_loop = .false.
    net_area_previous = 0
    
    do while(stop_loop .eq. .false.)
    
        call locateLineIndex(min_index, min_press_index, lower_shock_line_index, time_sheared(upper_shock_line_index), nt, time_sheared)
        call locateLineIndex(min_press_index, max_press_index, mid_shock_line_index, time_sheared(upper_shock_line_index), nt, time_sheared)
        
        call integrate(area_lower_1, lower_shock_line_index, min_press_index, time_sheared(upper_shock_line_index), nt, time_sheared, pressure_thickness)
        call integrate(area_lower_2, min_press_index, mid_shock_line_index, time_sheared(upper_shock_line_index), nt, time_sheared, pressure_thickness)
        
        area_lower = abs(area_lower_2 - area_lower_1)
        
        call integrate(area_upper_1, mid_shock_line_index, max_press_index, time_sheared(upper_shock_line_index), nt, time_sheared, pressure_thickness)
        call integrate(area_upper_2, max_press_index, upper_shock_line_index, time_sheared(upper_shock_line_index), nt, time_sheared, pressure_thickness)
        
        area_upper = area_upper_2 - area_upper_1 
        net_area = area_upper + area_lower
        
        if(net_area * net_area_previous .lt. 0) then
            stop_loop = .true.
        else if(net_area .lt. 0) then
            upper_shock_line_index = upper_shock_line_index - 1
        else if(net_area .gt. 0) then
            upper_shock_line_index = upper_shock_line_index + 1
        else if(net_area .eq. 0) then
            stop_loop = .true.
        end if
        write(*,*) net_area
        net_area_previous = net_area
        
    end do
    
    upper_shock_line_index = upper_shock_line_index 
    lower_shock_line_index = lower_shock_line_index
    
    call shockLine(nt, time_sheared, lower_shock_line_index, upper_shock_line_index, pressure_thickness)
    
end subroutine applyShock
!*****************************************************************************
subroutine locateLineIndex(start_index, end_index, locate_index, target_value, n, array)
    integer,intent(in):: start_index,end_index,n
    integer,intent(out):: locate_index
    real,intent(in):: target_value
    real,dimension(n),intent(in):: array
    
    real:: test, start_value, current_value, start_diff, current_diff
    integer:: i
    
    test = 1.0
    i = start_index
    
    do while(test .gt. 0)
    
        i = i + 1
        start_value = array(start_index)
        start_diff = start_value - target_value
        current_value = array(i)
        current_diff = current_value - target_value
        test = start_diff * current_diff
        if(i .gt. end_index) then
            write(*,*) "ERROR: Index exceeds end_index"
            exit
        end if
        
    end do
    
    locate_index = i - 1
    
end subroutine locateLineIndex
!*****************************************************************************
subroutine integrate(area,start_index,end_index,time_initial,n,time,pressure)
    integer,intent(in):: start_index,end_index,n
    real,intent(in):: time_initial
    real,intent(out):: area
    real,dimension(n),intent(in):: time,pressure
    
    integer:: i
    
    area = 0
    
    do i = start_index + 1,end_index
        area = area + (abs(time(i) + time(i - 1) - 2 * time_initial) * abs(pressure(i) - pressure(i-1)) / 2)
    end do
    
end subroutine integrate
!*****************************************************************************
subroutine shockLine(nt,time_sheared,time_shock_bottom_index,time_shock_top_index,pressure_thickness)
    integer,intent(in):: nt,time_shock_bottom_index,time_shock_top_index
    real,dimension(nt),intent(inout):: time_sheared,pressure_thickness
    
    integer:: i
    real:: slope_time,y_int_time,slope_pressure,y_int_pressure
    real:: pressure_top,pressure_bottom
    
    slope_time = (time_sheared(time_shock_top_index) - time_sheared(time_shock_bottom_index)) / (time_shock_top_index - time_shock_bottom_index) 
    y_int_time = time_sheared(time_shock_top_index) - slope_time * time_shock_top_index
    pressure_top = pressure_thickness(time_shock_top_index)
    pressure_bottom = pressure_thickness(time_shock_bottom_index)
    slope_pressure = (pressure_top - pressure_bottom) / (time_shock_top_index - time_shock_bottom_index)
    y_int_pressure = pressure_top - slope_pressure * time_shock_top_index
    
    do i = time_shock_bottom_index,time_shock_top_index
    
        time_sheared(i) = slope_time * i + y_int_time
        pressure_thickness(i) = slope_pressure * i + y_int_pressure
        
    end do
    
end subroutine shockLine
!*****************************************************************************
end module
