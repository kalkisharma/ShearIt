module ReadIn
    use JobConstants
    implicit none
    
    public:: ReadParameters
    public:: readPressureTimeHistory
    
    contains
!*****************************************************************************
subroutine readParameters(chord,mach,nt)
    integer,intent(out):: nt
    real,intent(out):: mach,chord
    
    open(unit = HSI_unit, file = HSI_file)
    read(HSI_unit,*) nt
    read(HSI_unit,*) mach
    read(HSI_unit,*) chord
    close(HSI_unit)
end subroutine ReadParameters
!*****************************************************************************
subroutine readPressureTimeHistory(nt,time_input,pressure_thickness,pressure_loading,pressure_total)
    integer,intent(in):: nt
    real,dimension(nt),intent(out):: time_input,pressure_thickness,pressure_loading,pressure_total
    
    integer:: i
    
    open(unit = tecplot_unit_in, file = tecplot_file_in)
    read(tecplot_unit_in,*) 
    read(tecplot_unit_in,*) 
    read(tecplot_unit_in,*) 
    ! ERROR: Four lines in text file, but only need to read three lines before data. WHY?
!    read(tecplot_unit,*) 
    read(tecplot_unit_in,*) (time_input(i),pressure_thickness(i),pressure_loading(i),pressure_total(i),i=1,nt)
    close(tecplot_unit_in)
end subroutine readPRessureTimeHistory
!*****************************************************************************
end module