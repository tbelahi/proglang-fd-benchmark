program test

      implicit none

      character(len=256) :: s
      integer :: i,l,n,m
      integer, allocatable :: b(:,:), c(:,:)

      s='salut les copains'

      print *, s
      
      n = 5
      m = 10
      allocate(b(m,n), c(m,n))
      b = 1
      c = 2*b

      do i=1,m
      write(*,'(100i2)')  (c(i,l), l=1,n) 
      enddo
      print *, b

end program
