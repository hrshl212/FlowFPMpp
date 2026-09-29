module fpm_trace
  ! No-op range markers for CPU builds; OpenACC builds link fpm_trace_nvtx.f90.
  implicit none
  private
  public::trace_push,trace_pop
contains
  subroutine trace_push(name)
    character(*),intent(in)::name
    if(len(name)<0)return
  end subroutine trace_push

  subroutine trace_pop()
  end subroutine trace_pop
end module fpm_trace
