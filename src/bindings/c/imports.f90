module imports_c
  implicit none
  private

contains
!
!###################################################################################
!
  subroutine import_capillary_c(FLOWFILE, filename_len) bind(C, name="import_capillary_c")
  
    use iso_c_binding, only: c_ptr
    use utils_c, only: strncpy
    use imports, only: import_capillary
    use other_consts, only: MAX_STRING_LEN, MAX_FILENAME_LEN
    implicit none
    integer,intent(in) :: filename_len
    type(c_ptr), value, intent(in) :: FLOWFILE
    character(len=MAX_FILENAME_LEN) :: filename_f

    call strncpy(filename_f, FLOWFILE, filename_len)
    call import_capillary(filename_f)

  end subroutine import_capillary_c
!
!###################################################################################
!
  subroutine import_terminal_c(FLOWFILE, filename_len) bind(C, name="import_terminal_c")
  
    use iso_c_binding, only: c_ptr
    use utils_c, only: strncpy
    use imports, only: import_terminal
    use other_consts, only: MAX_STRING_LEN, MAX_FILENAME_LEN
    implicit none
    integer,intent(in) :: filename_len
    type(c_ptr), value, intent(in) :: FLOWFILE
    character(len=MAX_FILENAME_LEN) :: filename_f

    call strncpy(filename_f, FLOWFILE, filename_len)
    call import_terminal(filename_f)

  end subroutine import_terminal_c
!
!###################################################################################
!
  subroutine import_ventilation_c(FLOWFILE, filename_len) bind(C, name="import_ventilation_c")
  
    use iso_c_binding, only: c_ptr
    use utils_c, only: strncpy
    use imports, only: import_ventilation
    use other_consts, only: MAX_STRING_LEN, MAX_FILENAME_LEN
    implicit none
    integer,intent(in) :: filename_len
    type(c_ptr), value, intent(in) :: FLOWFILE
    character(len=MAX_FILENAME_LEN) :: filename_f

    call strncpy(filename_f, FLOWFILE, filename_len)

    call import_ventilation(filename_f)

  end subroutine import_ventilation_c
  
!
!###################################################################################
!
  subroutine import_perfusion_c(FLOWFILE, filename_len) bind(C, name="import_perfusion_c")

    use iso_c_binding, only: c_ptr
    use utils_c, only: strncpy
    use imports, only: import_perfusion
    use other_consts, only: MAX_STRING_LEN, MAX_FILENAME_LEN
    implicit none
    integer,intent(in) :: filename_len
    type(c_ptr), value, intent(in) :: FLOWFILE
    character(len=MAX_FILENAME_LEN) :: filename_f

    call strncpy(filename_f, FLOWFILE, filename_len)

    call import_perfusion(filename_f)

  end subroutine import_perfusion_c

!
!###################################################################################
!
 subroutine import_exnodefield_c(NODEFILE,filename_len) bind(C, name="import_exnodefield_c")
    use iso_c_binding, only: c_ptr
    use utils_c, only: strncpy
    use arrays, only: dp    
    use imports, only: import_exnodefield
    use other_consts, only: MAX_STRING_LEN, MAX_FILENAME_LEN
    implicit none
    integer,intent(in) :: filename_len
    type(c_ptr), value, intent(in) :: NODEFILE
    character(len=MAX_FILENAME_LEN) :: filename_f

    call strncpy(filename_f, NODEFILE, filename_len)

    call import_exnodefield(filename_f)

  end subroutine import_exnodefield_c
      
    
 

end module imports_c

