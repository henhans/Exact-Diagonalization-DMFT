!+------------------------------------------------------+
!PROGRAM  : 
!TYPE     : subroutine
!PURPOSE  :On OUTPUT \Delta = g
!+------------------------------------------------------+
subroutine getdelta_bethe_analytic()
  if(trim(order)=="lro")then
     if(allocated(deltaup))deallocate(deltaup)
     if(allocated(deltadw))deallocate(deltadw)
     allocate(deltaup(NL),deltadw(NL))
     deltaup=(d**2/4.d0)*Giwup
     deltaup=(d**2/4.d0)*Giwdw
  else
     if(allocated(delta))deallocate(delta)
     allocate(delta(NL))
     delta=(d**2/4.d0)*Giw
  endif
  return    
end subroutine getdelta_bethe_analytic
!********************************************************
!********************************************************
!********************************************************
