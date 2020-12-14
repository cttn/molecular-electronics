module plotmod
 implicit none
 integer,parameter :: jtotal=100 !! Nº de archivos
end module plotmod

program gnuplotconf
 use plotmod
 implicit none
  integer :: i,j,fi
  real(kind(0d0)) :: x,y,ang,coseno,seno,T
  character(len=20) :: entrada,salida,ext*10

 do j=1,jtotal
   write(ext,'(I10)') j*10
   entrada='fort.'//adjustl(ext)
   salida ='sal.'//adjustl(ext)

   open(j*10,file=entrada,status='old')
   open(j*1000,file=salida)
     write(j*1000,*) 'set xrange[-1:34]'
     write(j*1000,*) 'set yrange[-1:34]'
     write(j*1000,*) 'unset xtics'
     write(j*1000,*) 'unset ytics'
     write(j*1000,*) 'unset key'

     do i=1,1024
       read(j*10,1130) x,y,ang,coseno,seno,T
       write(j*1000,1131)'set arrow from ',x,',',y,' to ',x+coseno,',',y+seno
     enddo

     write(j*1000,'(A23,F8.3)') 'set title "Temperatura:',T
     write(j*1000,*) 'plot 500'
   close(j*10)
   close(j*1000)
 enddo

 1130 format(6E16.8)
 1131 format(A16,F4.1,A,F4.1,A5,F5.2,A,F5.2)
end program gnuplotconf
