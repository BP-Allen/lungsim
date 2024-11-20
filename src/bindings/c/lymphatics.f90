module lymphatics_c
  implicit none

  private

contains
!
!###################################################################################
!
!*alveolar_capillary_flux:* 
  subroutine alveolar_capillary_flux_c(num_nodes) bind(C, name="alveolar_capillary_flux_c")
    use lymphatics,only: alveolar_capillary_flux
    implicit none

    integer,intent(in) :: num_nodes

    call alveolar_capillary_flux(num_nodes)
    
  end subroutine alveolar_capillary_flux_c
!
!###################################################################################
!
!*lymphatic_transport* 
  subroutine lymphatic_transport_c(filename, filename_len) bind(C, name="lymphatic_transport_c")
    use iso_c_binding, only: c_ptr
    use utils_c, only: strncpy
    use other_consts, only: MAX_FILENAME_LEN
    use lymphatics,only: lymphatic_transport
    implicit none
    
    integer,intent(in) :: filename_len
    type(c_ptr), value, intent(in) :: filename
    character(len=MAX_FILENAME_LEN) :: filename_f

    call strncpy(filename_f, filename, filename_len)
    call lymphatic_transport(filename_f)
    
  end subroutine lymphatic_transport_c

!
!###################################################################################
!

end module lymphatics_c
