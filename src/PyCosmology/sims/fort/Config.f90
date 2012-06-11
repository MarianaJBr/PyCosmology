module Config
  implicit none
  
  
  !Input File names
    character(len=80), parameter :: DIRECTORY = "/Users/smurray/Documents/PhD/"&
                                               &//"Simulations/data_for_steven/"
    character(len=50), parameter :: sim_file = "snapshot_130"
    character(len=50), parameter :: fof_file = "groups/groups_130.fofcat.for"
    character(len=50), parameter :: ids_file = "groups/groups_130.pos.ids"
    
  !Output File names
    character(len=50), parameter :: pos_file = 's130_pos.dat'
    character(len=50), parameter :: slice_file = "s130_pos_zslice.dat"
    character(len=50), parameter :: vel_file = 's130_vel.dat'
    character(len=50), parameter :: pos_file_c = 's130_pos_centred.dat'
    character(len=50), parameter :: pos_file_p = 's130_pos_principle_axes.dat'
    character(len=50), parameter :: pos_file_r = 's130_pos_rotated.dat'
    character(len=50), parameter :: dens_file = 's130_densities.dat'
    character(len=50), parameter :: vel_file_r = 's130_vel_rotated.dat'
    
  !Position filtering options for filtered output
    logical, parameter          :: pos_filter   = .FALSE.
    integer, parameter          :: filter_coord = 3
    real, parameter             :: min_filt     = 9.8
    real, parameter             :: max_filt     = 10.2
    
  !Halo groups comparison options
    logical, parameter          :: find_groups = .TRUE.
    logical, parameter          :: center_group = .TRUE.
    integer, parameter          :: number_of_groups = 1 !In order of size
    character(len=50), parameter:: group_dir = "groups/GroupData/Group_"
    logical, parameter          :: rotate   = .TRUE.
    real, parameter             :: u(3) = (/0.,0.,1.0/)
    real, parameter             :: theta = 3.14
    real, parameter             :: min_bin = 0.002 !Mpc
    integer, parameter          :: bins = 10
end module
