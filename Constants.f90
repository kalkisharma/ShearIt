module JobConstants
    implicit none
    
    ! file unit numbers
    integer,parameter:: HSI_unit = 10
    integer,parameter:: tecplot_unit_in = 11
    integer,parameter:: tecplot_unit_out = 12
    ! math constants
    real,parameter:: deg2rad = 0.01745329251
    ! filenames
    character(len=17),parameter:: HSI_file = "HSIParameters.txt"
    character(len=12),parameter:: tecplot_file_in = "pressure.tec"
    character(len=20),parameter:: tecplot_file_out = "pressure_sheared.plt"
end module