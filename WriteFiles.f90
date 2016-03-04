module WriteOut
    use JobConstants
    implicit none
    
    public:: writeTecplot
    
    contains
!*****************************************************************************
subroutine writeTecplot(nt,time_sheared,pressure_thickness,pressure_loading,pressure_total)
    integer,intent(in):: nt
    real,dimension(nt),intent(in):: time_sheared, pressure_thickness, pressure_loading, pressure_total
    
    integer :: i
	
	pressure_HSI = pressure_total - (pressure_loading + pressure_thickness)
	
	! Tecplot file to be displayed
    open(unit = tecplot_unit_out, file = tecplot_file_out)
    write(tecplot_unit_out,'(a19)') 'TITLE="PSU-WOPWOP"'
    write(tecplot_unit_out,'(a162)') 'VARIABLES= "Observer Time (s)", "Thickness Acoustic Pressure (Pa)", "Loading Acoustic Pressure (Pa)", "HSI Acoustic Pressure (Pa)", "Total Acoustic Pressure (Pa)"'
    write(tecplot_unit_out,'(a185)') 'TEXT X=100 Y=95 AN=HEADRIGHT T="Observer Time (s) vs \\nThickness Acoustic Pressure (Pa)\\nLoading Acoustic Pressure (Pa)\\nHSI Acoustic Pressure (Pa)\\nTotal Acoustic Pressure (Pa)\\ "'
    write(tecplot_unit_out,*) 
    write(tecplot_unit_out,'(ES15.8,ES16.8,ES16.8,ES16.8,ES16.8)') (time_sheared(i), pressure_thickness(i), pressure_loading(i), pressure_HSI(i), pressure_total(i), i=1,nt)
    Close(tecplot_unit_out)
	
	open(unit = tecplot_unit_out_wopwop, file = tecplot_file_out_wopwop)
    write(tecplot_unit_out_wopwop,'(a19)') 'TITLE="PSU-WOPWOP"'
    write(tecplot_unit_out_wopwop,'(a162)') 'VARIABLES= "Observer Time (s)", "Thickness Acoustic Pressure (Pa)", "Loading Acoustic Pressure (Pa)", "HSI Acoustic Pressure (Pa)", "Total Acoustic Pressure (Pa)"'
    write(tecplot_unit_out_wopwop,'(a185)') 'TEXT X=100 Y=95 AN=HEADRIGHT T="Observer Time (s) vs \\nThickness Acoustic Pressure (Pa)\\nLoading Acoustic Pressure (Pa)\\nHSI Acoustic Pressure (Pa)\\nTotal Acoustic Pressure (Pa)\\ "'
    write(tecplot_unit_out_wopwop,*) 
    write(tecplot_unit_out_wopwop,'(ES15.8,ES16.8,ES16.8,ES16.8,ES16.8)') (time_sheared(i), pressure_thickness(i), pressure_loading(i), pressure_HSI(i), pressure_total(i), i=1,nt)
    Close(tecplot_unit_out_wopwop)
end subroutine writeTecplot
!*****************************************************************************
end module