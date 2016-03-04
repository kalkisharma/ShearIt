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
    
    integer:: nt,nt_fine !Time Steps 
    real:: mach,chord !Mach number and Chord
    real,dimension(:),allocatable:: time_input, pressure_thickness, pressure_adjustable_nondim
    real,dimension(:),allocatable:: pressure_adjustable_amplified, time_sheared, pressure_loading
	real,dimension(:),allocatable:: pressure_total, pressure_adjustable_spaced, pressure_adjustable
	real,dimension(:),allocatable:: time_input_fine, pressure_adjustable_fine
    real,dimension(:),allocatable:: time, pressure
	
    ! Reads in the mach number and time steps from the HSIParameters.txt file
    call readParameters(chord,mach,nt,nt_fine)
    
    allocate(time_input(nt),pressure_thickness(nt), pressure_loading(nt), pressure_total(nt), pressure_adjustable(nt))
	allocate(pressure_adjustable_nondim(nt_fine), pressure_adjustable_amplified(nt_fine)) 
	allocate(time_sheared(nt_fine), pressure_adjustable_spaced(nt_fine))
    allocate(time_input_fine(nt_fine), pressure_adjustable_fine(nt_fine))
	allocate(time(nt), pressure(nt))
	
    ! Reads the time and pressure values from the tecplot file
    call readPressureTimeHistory(nt, time_input, pressure_thickness, pressure_loading, pressure_total)
	
	! Acoustic pressure history to be altered
	pressure_adjustable = pressure_total
	
	! Increases resolution of signal 
	call fineSignal(nt, time_input, pressure_adjustable, nt_fine, time_input_fine, pressure_adjustable_fine)
	
    ! Nondimensionalize pressure array with respect to the negative peak pressure
    call nonDimensionalizeSignal(nt_fine, pressure_adjustable_fine, pressure_adjustable_nondim)
    
	! Shear the amplified signal by gamma angle value
    call shearSignal(nt_fine, mach, chord, pressure_adjustable_fine, pressure_adjustable_nondim, &
						time_input_fine, time_sheared)

    ! Multiplies thickness pressure by amplification factor dependent on Mach #
    call amplifySignal(nt_fine, mach, pressure_adjustable_fine, pressure_adjustable_amplified)
    
	! Changes fine signal back to original resolution
	call coarseSignal(nt, time, pressure, nt_fine, time_sheared, pressure_adjustable_amplified)
	
	! Spaces out adjusted signal to same spacing as inputed signal
	call spaceSignal(nt, pressure, pressure_adjustable_spaced, time_input, time)
	
    ! Writes out tecplot file
    call writeTecplot(nt, time_input, pressure_thickness, pressure_loading, pressure_adjustable_spaced)
    
    deallocate(time_input, pressure_thickness)
	deallocate(pressure_adjustable_nondim, pressure_adjustable_amplified, time_sheared) 
	deallocate(pressure_loading, pressure_total, pressure_adjustable_spaced)
    deallocate(time_input_fine, pressure_adjustable_fine)
	deallocate(time, pressure)
end program ShearIt