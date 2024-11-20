module field_utilities
  !*Brief Description:* This module contains all the subroutines that perform general operations
  !on fields
  !*LICENSE:*
  !
  !
  !
  !*Full Description:*
  !More info on what the module does if necessary
  !
  use arrays
  use diagnostics
  use geometry
  use indices
  use other_consts
  use precision
  
  implicit none
  
  !Module parameters
  
  !Module types
  
  !Module variables
  
  !Interfaces
  private
  public scale_flow_to_inlet,calculate_blood_volume
  
contains

!!!#############################################################################
  
  subroutine scale_flow_to_inlet(inlet_flow,VorQ)
    !*scale_flow_field:* Scales a flow field to an 'inlet flow' value (real units).

    real(dp),intent(in) :: inlet_flow
    character(len=1), intent(in) :: VorQ
    ! Local variables
    real(dp) :: ratio
    character(len=60) :: sub_name

    ! --------------------------------------------------------------------------

    sub_name = 'scale_flow_to_inlet'
    call enter_exit(sub_name,1)

    if(VorQ.eq.'V')then
      if(abs(elem_field(ne_Vdot,1)).gt.zero_tol)then
         ratio = inlet_flow/elem_field(ne_Vdot,1)
         unit_field(nu_Vdot0,1:num_units) = unit_field(nu_Vdot0,1:num_units)*ratio
         elem_field(ne_Vdot,1:num_elems) = elem_field(ne_Vdot,1:num_elems)*ratio
         elem_field(ne_dvdt,1:num_elems) = elem_field(ne_dvdt,1:num_elems)*ratio
      else
         write(*,'('' Cannot scale to zero flow'')')
      endif
    elseif(VorQ.eq.'Q')then
      if(abs(elem_field(ne_Qdot,1)).gt.zero_tol)then
        ratio = inlet_flow/elem_field(ne_Qdot,1)
        unit_field(nu_perf,1:num_units) = unit_field(nu_perf,1:num_units)*ratio
        elem_field(ne_Qdot,1:num_elems) = elem_field(ne_Qdot,1:num_elems)*ratio
      else
        write(*,'('' Cannot scale to zero flow'')')
      endif
    else
     print*, ' WARNING: Not a ventilation or perfusion scaling, no calculations done'
    endif

    call enter_exit(sub_name,2)
    
  end subroutine scale_flow_to_inlet

!!!#############################################################################

  subroutine calculate_blood_volume(out_format) 
    !*calculate_blood_volume:* Calculates blood volumes for a range of arterial vessel sizes.

    character(len=4),intent(in) :: out_format ! whether to output results as continuous ('cont') or discrete ('disc')
    ! Local variables
    character(len=60) :: sub_name
    integer :: ne,i,bn,vessel_type
    real(dp) :: L, R, area,vol,cbv
    real(dp) :: bins(18), bin_storage(19,2)
    character(len=7):: string(19)   
    bins = [5.0_dp, 10.0_dp, 15.0_dp, 20.0_dp, 25.0_dp, 30.0_dp, 35.0_dp, 40.0_dp, 45.0_dp, 50.0_dp, 55.0_dp,&
    &60.0_dp, 65.0_dp, 70.0_dp, 75.0_dp, 80.0_dp, 85.0_dp, 90.0_dp]
    bin_storage = 0.0_dp
    cbv = 0.0_dp
    string = ['BV00_05','BV05_10','BV10_15','BV15_20','BV20_25','BV25_30','BV30_35','BV35_40','BV40_45','BV45_50',&
    &'BV50_55','BV55_60','BV60_65','BV65_70','BV70_75','BV75_80','BV80_85','BV85_90','BV90_up']
    ! --------------------------------------------------------------------------

    sub_name = 'calculate_blood_volume'

    call enter_exit(sub_name,1)
    
    if (out_format.eq.'disc') then
      do ne=1,num_elems
        L = elem_field(ne_length,ne)
        !! Calculating the area and determining vessel type from current element
        R = elem_field(ne_radius,ne)  
        if (elem_field(ne_group,ne).eq.0) then
           area = PI*R**2.0_dp
           vol = area*L
           vessel_type = 1
        elseif (elem_field(ne_group,ne).eq.2) then
           area = PI*R**2.0_dp
           vol = area*L
           vessel_type = 2
        else
           vessel_type = 0
        endif
        !! Categorising volume based on the area of the vessel
        if (vessel_type.ne.0) then
          bin_find: do bn=1,18 
             if (area.lt.bins(bn)) then !! store volumes/areas for EACH element
                bin_storage(bn,vessel_type) = bin_storage(bn,vessel_type) + vol
                exit bin_find
             elseif (area.gt.bins(18)) then
                bin_storage(19,vessel_type) = bin_storage(19,vessel_type) + vol
                exit bin_find               
             endif
          enddo bin_find
        else
          cbv = cbv + elem_field(nu_vol,ne)
        endif
      enddo 
      101 format(1x, *(g0, ","))        
      open(10, file='TBV_disc.out', status='replace') 
      write(10, 101) string
      write(10, *) "Arterial Volumes"
      write(10, 101) bin_storage(:,1)
      write(10, *) "Venous Volumes"       
      write(10, 101) bin_storage(:,2)
      write(10, *) "Capillary Volume"
      write(10, *) cbv
      close(10)
    elseif (out_format.eq.'cont') then  
       open(10, file='TBV_cont.exelem', status='replace')
       write(10,'( '' Group name: Blood_volumes'' )')
       write(10,'( '' Shape.  Dimension=1'' )')
       write(10,'( '' #Scale factor sets=0'' )')
       write(10,'( '' #Nodes= 0'' )')
       write(10,'( '' #Fields= 1'' )')
       write(10,'( '' 1) volume, field, rectangular cartesian, #Components=1'')')
       write(10,'(''   volume.  l.Lagrange, no modify, grid based.'')')                 
       write(10,'( ''     #xi1=1'')')
       do ne=1,num_elems
          if (elem_field(ne_group,ne).ne.1) then
            L = elem_field(ne_length,ne)
            !! Calculating the area and determining vessel type from current element
            R = elem_field(ne_radius,ne)    
            area = PI*R**2.0_dp
            vol = area*L
            !**               write the element
            write(10,'(1X,''Element: '',I12,'' 0 0'' )') ne
            !**                 write the scale factors
            write(10,'(3X,''Values:'' )')
            write(10,'(4X,2(1X,E12.5))') vol,vol
          else !! Capillary Volume          
            vol = elem_field(nu_vol,ne)
            !**               write the element
            write(10,'(1X,''Element: '',I12,'' 0 0'' )') ne
            !**                 write the scale factors
            write(10,'(3X,''Values:'' )')
            write(10,'(4X,2(1X,E12.5))') vol,vol
          endif
       enddo 
       close(10)
    endif
    
  end subroutine calculate_blood_volume
   



end module field_utilities
