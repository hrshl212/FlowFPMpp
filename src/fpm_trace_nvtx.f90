module fpm_trace
  ! NVTX range markers for Nsight Systems timelines (NVHPC, link -cudalib=nvtx).
  use nvtx,only:nvtxStartRange,nvtxEndRange
  implicit none
  private
  public::trace_push,trace_pop
contains
  subroutine trace_push(name)
    character(*),intent(in)::name
    call nvtxStartRange(name)
  end subroutine trace_push

  subroutine trace_pop()
    call nvtxEndRange()
  end subroutine trace_pop
end module fpm_trace
