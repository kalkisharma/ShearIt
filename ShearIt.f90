!*****************************************************************************
!Author: Kalki Sharma
!Date: 1/27/2016
!Purpose: This program takes in a pressure time history tecplot file and then 
!   amplifies, shears, and widens the thickness pressure time history using 
!   ad hoc methods to mimic high speed impulsive (HSI) effects
!Input Files: .tec file & HSIParameters.txt
!*****************************************************************************
program ShearIt
    use ReadIn
    use Signal
    use WriteOut
    implicit none
    
    integer:: nt !Time Steps 
    real:: mach,chord !Mach number and Chord
    real,dimension(:),allocatable:: time_input,pressure_thickness,pressure_thickness_nondim,pressure_thickness_amplified,time_sheared,pressure_loading,pressure_total,pressure_thickness_spaced
    
    ! Reads in the mach number and time steps from the HSIParameters.txt file
    call readParameters(chord,mach,nt)
    
    allocate(time_input(nt),pressure_thickness(nt),pressure_thickness_nondim(nt),pressure_thickness_amplified(nt),time_sheared(nt),pressure_loading(nt),pressure_total(nt),pressure_thickness_spaced(nt))
    
    ! Reads the time and pressure values from the tecplot file
    call readPressureTimeHistory(nt,time_input,pressure_thickness,pressure_loading,pressure_total)
    ! Nondimensionalize pressure array with respect to the negative peak pressure
    call nonDimensionalizeSignal(nt,pressure_thickness,pressure_thickness_nondim)
    ! Shear the amplified signal by gamma angle value
    call shearSignal(nt,mach,chord,pressure_thickness,pressure_thickness_nondim,time_input,time_sheared)
    ! Multiplies thickness pressure by amplification factor dependent on Mach #
    call amplifySignal(nt,mach,pressure_thickness,pressure_thickness_amplified)
    call spaceSignal(nt,pressure_thickness_amplified,pressure_thickness_spaced,time_input,time_sheared)
    ! Writes out tecplot file
    call writeTecplot(nt,time_input,pressure_thickness_spaced,pressure_loading,pressure_total)
    
    deallocate(time_input,pressure_thickness,pressure_thickness_nondim,pressure_thickness_amplified,time_sheared,pressure_loading,pressure_total,pressure_thickness_spaced)
    
end program ShearIt