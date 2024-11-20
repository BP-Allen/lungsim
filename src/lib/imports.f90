module imports
!*Brief Description:* This module contains all the subroutines required to
!import fields, previous model results, etc.
!*LICENSE:*
!
!
!
!*Full Description:*
!
  !
  use arrays
  use diagnostics
  use geometry
  use indices
  use other_consts
  use ventilation
  
  implicit none

  !Module parameters

  !Module types

  !Module variables

  !Interfaces
  private 
  public import_capillary
  public import_ventilation
  public import_perfusion
  public import_terminal
  public import_exnodefield

contains
 !
  !###########################################################################################
  !
  !>*import_capillary:* This subroutine reads in the results for the micro-circulatory
  ! components (capillary bed within 'units') of a perfusion model .
  subroutine import_capillary(micro_unit_file)
    
    character(len=MAX_FILENAME_LEN),intent(in) :: micro_unit_file
    !local variables
    integer :: count_units,ios,ne,nunit
    real(dp) :: Pin,Pout,Qtot,temp(7),TOTAL_CAP_VOL,TOTAL_SHEET_SA,TT_TOTAL
    
    character(len=60) :: sub_name
    
    sub_name = 'import_capillary'
    call enter_exit(sub_name,1)
    
    open(10, file = micro_unit_file, status='old')
    
    ! ios is negative if an end of record condition is encountered or if
    ! an endfile condition was detected.  It is positive if an error was
    ! detected.  ios is zero otherwise.
    ios = 0
    count_units = 0
    do while (ios == 0) ! temp = (x,y,z,Q01_mthrees*1.d9,Rtot/1000.d0**3.d0,TOTAL_SHEET_H,Ppl)
       read(10, '(I6,X,5(F9.2,X),2(F8.5,X),F10.2,X,F8.4,X,F10.4,X,F10.3,X,F8.4,X,F9.4,X)', iostat=ios)  &
            ne,temp(1:3),Pin,Pout,temp(4),Qtot,temp(5),TOTAL_CAP_VOL,TOTAL_SHEET_SA,TT_TOTAL,temp(6:7)
       if(ios == 0)then
          ! record the unit values for mean pressure, transit time, surface area.
          ! ne is the 'linker' element, so nunit is for its parent element
          nunit = int(elem_field(ne_unit,elem_cnct(-1,1,ne))) 
          unit_field(nu_blood_press,nunit) = (Pin+Pout)/2.0_dp
          unit_field(nu_tt,nunit) = TT_TOTAL
          unit_field(nu_vol,nunit) = TOTAL_CAP_VOL
          unit_field(nu_sheet_area,nunit) = TOTAL_SHEET_SA
          count_units = count_units + 1
       endif
    enddo

    if(count_units.ne.num_units)then
       write(*,'(''WARNING: the number of capillary units ('',i6,'') does not match the geometric model ('',i6,'')'')') &
            count_units,num_units
    endif

    close(10)
    
   call enter_exit(sub_name,2)
   
 end subroutine import_capillary
!
!##############################################################################
!
!>*import_ventilation:* This subroutine reads in the results of a ventilation model that
! has been saved in an exelem format as a single flow field (elements listed with
! ventilation as field values).
 subroutine import_ventilation(FLOWFILE)

   character(len=MAX_FILENAME_LEN),intent(in) :: FLOWFILE
   !local variables
   integer :: ierror,ne,nunit
   character(LEN=132) :: ctemp1,exfile
   real(dp) :: flow,flow_unit,maxflow

   character(len=60) :: sub_name

   sub_name = 'import_ventilation'
   call enter_exit(sub_name,1)

   print *, 'Reading in ventilation results'
   call import_exelemfield(FLOWFILE,ne_Vdot)
   do nunit = 1,num_units
     ne = units(nunit)
     if(elem_field(ne_Vdot,ne).lt.0.0_dp) elem_field(ne_Vdot,ne) = zero_tol
     unit_field(nu_Vdot0,nunit) = elem_field(ne_Vdot,ne)
   enddo

!!! sum the fields up the tree
   call sum_elem_field_from_periphery(ne_Vdot) !sum the air flows recursively UP the tree
   maxflow = elem_field(ne_Vdot,1)


   call enter_exit(sub_name,2)
 end subroutine import_ventilation

!
!###########################################################################################
!
!>*import_perfusion:* This subroutine reads in the results of a ventilation model that
! has been saved in an exelem format as a single flow field (elements listed with
! ventilation as field values).
 subroutine import_perfusion(FLOWFILE)

   character(len=MAX_FILENAME_LEN),intent(in) :: FLOWFILE
   !local variables
   integer :: ierror,ne,nunit
   character(LEN=132) :: ctemp1,exfile
   real(dp) :: flow,flow_unit,maxflow

   character(len=60) :: sub_name

   sub_name = 'import_perfusion'
   call enter_exit(sub_name,1)

   print *, 'Reading in perfusion results'
   call import_exelemfield(FLOWFILE,ne_Qdot)
   do nunit = 1,num_units
     ne = units(nunit)
     if(elem_field(ne_Qdot,ne).lt.0.0_dp) elem_field(ne_Qdot,ne) = zero_tol
     unit_field(nu_perf,nunit) = elem_field(ne_Qdot,ne)
   enddo

!!! sum the fields up the tree
   call sum_elem_field_from_periphery(ne_Qdot) !sum the air flows recursively UP the tree
   maxflow = elem_field(ne_Qdot,1)

   call enter_exit(sub_name,2)
 end subroutine import_perfusion

!
!##############################################################################
!
!>*import_exelemfield:* This subroutine reads in the content of an exelem field file (1 field)
 subroutine import_exelemfield(FLOWFILE,field_no)

   character(len=MAX_FILENAME_LEN),intent(in) :: FLOWFILE
   integer, intent(in) :: field_no
   !local variables
   integer :: ierror,ne,nunit
   character(LEN=132) :: ctemp1,exfile
   real(dp) :: flow,flow_unit,maxflow

   character(len=60) :: sub_name

   sub_name = 'import_exelemfield'
   call enter_exit(sub_name,1)

   open(10, file=FLOWFILE, status='old')
   ne = 0
   read_elem_flow : do !define a do loop name
     !.......read element flow
     read(unit=10, fmt="(a)", iostat=ierror) ctemp1
     if(index(ctemp1, "Values:")> 0) then
       ne = ne+1
       read(unit=10, fmt="(a)", iostat=ierror) ctemp1
       flow = get_final_real(ctemp1)
       if(flow.lt.0.0_dp) flow = zero_tol
         elem_field(field_no,ne) = flow! read it in
       end if
       if(ne.ge.num_elems) exit read_elem_flow
     end do read_elem_flow

   close(10)

   call enter_exit(sub_name,2)
   
 end subroutine import_exelemfield
 

!
!##############################################################################
!
!>*import_terminal:* This subroutine reads in the content of an exnode field file
 subroutine import_terminal(EXFILE)

   character(len=MAX_FILENAME_LEN),intent(in) :: EXFILE
   !local variables
   integer :: count_units,i,ibeg,iend,ierror,ne,nunit,n_fields
   integer :: np,np_global,field_label(10) !!!! SHOULD IN FUTURE HAVE FIELD_LABEL AS ALLOCATABLE
   real(dp) :: rtemp
   character(LEN=132) :: ctemp,label
   character(len=300) :: readfile

   character(len=60) :: sub_name

   sub_name = 'import_terminal'
   call enter_exit(sub_name,1)

   if(index(EXFILE, ".exnode")> 0) then !full filename is given
      readfile = EXFILE
   else ! need to append the correct filename extension
      readfile = trim(EXFILE)//'.exnode'
   endif
   
   open(10, file=readfile, status='old')

   n_fields = 0
   ierror = 0

   read_field_labels : do
      read(10, fmt="(a)", iostat=ierror) ctemp
      if(index(ctemp, "Node:")> 0) exit read_field_labels
      if(index(ctemp, ") ")> 0) then
         ibeg = index(ctemp, ") ")+1 ! beginning of label
         iend = index(ctemp, ",")-1  ! end of label
         label = adjustl(ctemp(ibeg:iend))
         n_fields = n_fields + 1
         field_label(n_fields) = 0
         
         select case(label)
         case('terminal_element')
         case('pleural_pressure')
            field_label(n_fields) = nu_pe
         case('tidal_volume')
            field_label(n_fields) = nu_vt ! THIS ISN'T INITIALISED IN LYMPHATIC INDICES????????
         case('max_Pe')
            field_label(n_fields) = nu_Pe_max
         case('min_Pe')
            field_label(n_fields) = nu_Pe_min
         end select
      endif
         
   end do read_field_labels

   ierror = 0
   
   do while (ierror == 0)
      if(index(ctemp, "Node:")> 0) then !! apparently no integer at the end here??
         do i = 1,3  ! read coordinates; not used 
            read(10, *, iostat=ierror) ctemp
         enddo
         
         ! THIS SHOULD BE USED TO ENSURE THE RIGHT ORDERING IS USED
         read(10, *, iostat=ierror) rtemp ! read element: used for setting the correct unit number

         nunit = elem_field(ne_unit,int(rtemp))
         
         !elem_field(ne_unit,int(rtemp)) = nunit !!! tobe used if counting the nunit values instead of the append_units() settings

         
         do i = 3,n_fields
            read(10, *, iostat=ierror) rtemp
            if(field_label(i).gt.0)then
               unit_field(field_label(i),nunit) = rtemp
            endif
         enddo
         
         read(10, *, iostat=ierror) ctemp
      endif
   enddo
      
   close(10)
   
   call enter_exit(sub_name,2)
   
 end subroutine import_terminal
!
!##############################################################################
!
!>*import_exnodefield:* This subroutine reads in the content of an exnode field file (N fields)
!! BE CAREFUL: This subroutine must either have the nodefile account for the entire network (arteries,veins,capillaries) or be done before the matching stage
 subroutine import_exnodefield(NODEFILE) 

   character(len=MAX_FILENAME_LEN),intent(in) :: NODEFILE
   !local variables
   integer :: ierror,np_global,ne,num_fields,counter,temp_index,np
   character(LEN=132) :: ctemp1
   character(LEN=10) :: temp_str
   character(len=60) :: sub_name
   ! temporary local variables
   type(character(len=10)), allocatable, dimension(:) :: field_categories
   real(dp),allocatable :: temp_values(:)
   real(dp) :: np1,np2
   
   sub_name = 'import_exnodefield'
   call enter_exit(sub_name,1)

   open(10, file=NODEFILE, status='old')
   read(unit=10, fmt="(a)", iostat=ierror) ctemp1   
   temp_index = 1
   read_node_file : do !define a do loop name
     !.......read node file
     read(unit=10, fmt="(a)", iostat=ierror) ctemp1
     if (ierror.eq.5001.or.ierror.eq.-1) exit read_node_file
     ! allocating the temporary arrays based on number of fields in the node file
     if(index(ctemp1, "Fields=")> 0) then
       read(ctemp1(10:11), *) num_fields
       num_fields = num_fields - 1
       if(allocated(temp_values)) deallocate(temp_values)
       if(allocated(field_categories)) deallocate(field_categories)
       allocate(temp_values(num_fields + 3))
       allocate(field_categories(num_fields))
       temp_index = 1
     elseif(index(ctemp1, "field")> 0) then
     ! storing the field categories (may be erroneously ordered)
       read(unit=10, fmt="(a)", iostat=ierror) ctemp1
       field_categories(temp_index) = ctemp1(1:index(ctemp1, "."))
       temp_index = temp_index + 1       
     ! storing the values for node fields 
     elseif(index(ctemp1, "Node:")> 0) then
       np_global = get_final_integer(ctemp1)
       call get_local_node(np_global,np) ! get local node np for global node
       read(unit=10, fmt="(a)", iostat=ierror) ctemp1
       read(ctemp1,*) temp_values 
       read_fields: do counter=1,num_fields
         temp_str = field_categories(counter)
         if(index(temp_str,'tumour')>0) then
           temp_index = nj_gtv
         elseif(index(temp_str,'dose')>0) then
           temp_index = nj_dose
         elseif(index(temp_str,'emph')>0) then
           temp_index = nj_emph
         endif 
         call set_node_field_value(temp_index,np,temp_values(counter + 3))
       enddo read_fields
     endif
     if(np.ge.num_nodes) exit read_node_file
   enddo read_node_file

   close(10)
   deallocate(temp_values)
   deallocate(field_categories)
   
   do ne=1,num_elems
      np1 = node_field(nj_emph,elem_nodes(1,ne))
      np2 = node_field(nj_emph,elem_nodes(2,ne))
      if (np1.le.-950_dp .or.np2.le.-950_dp)then
         !elem_field(ne_emph,ne) = 0.9_dp -- BPOA masters versions
         !elem_field(ne_emph_c,ne) = 0.92_dp 
         elem_field(ne_emph,ne) = -0.001_dp*(np1 + np2)/2.0_dp - 0.05_dp ! Averaging may not be the best course of action here
         elem_field(ne_emph_c,ne) = -0.0012_dp*(np1 + np2)/2.0_dp - 0.22_dp         
      endif
   enddo
   
   call enter_exit(sub_name,2)
 end subroutine import_exnodefield


end module imports
