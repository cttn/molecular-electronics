program gnuplotter
use plotmod
implicit none
 integer :: j
 character :: archivo*9,ext*5

 open(44,file='gnuplot.me')
 do j=1,jtotal
    write(ext,'(I5)') j*10

   archivo ='sal.'//adjustl(ext)

   if (Len_Trim(archivo) .eq. 6) then
     write(44,'(A8,A6,A)')"load './",archivo(1:len_trim(archivo)),"'"
   elseif (Len_Trim(archivo) .eq. 7) then
     write(44,'(A8,A7,A)')"load './",archivo(1:len_trim(archivo)),"'"
   elseif (Len_Trim(archivo) .eq. 8) then
     write(44,'(A8,A8,A)')"load './",archivo(1:len_trim(archivo)),"'"
   elseif (Len_Trim(archivo) .eq. 9) then
     write(44,'(A8,A9,A)')"load './",archivo(1:len_trim(archivo)),"'"
   elseif (Len_Trim(archivo) .eq. 10) then
     write(44,'(A8,A10,A)')"load './",archivo(1:len_trim(archivo)),"'"
   else
     stop 'situacion no contemplada'
   endif

   write(44,'(A26)')"set term jpeg size 700,700"

   if (len_trim(archivo) .eq. 6) then
     write(44,'(A12,A6,A5)')"set output '",archivo(1:len_trim(archivo)),".jpg'"
   elseif (len_trim(archivo) .eq. 7) then
     write(44,'(A12,A7,A5)')"set output '",archivo(1:len_trim(archivo)),".jpg'"
   elseif (len_trim(archivo) .eq. 8) then
     write(44,'(A12,A8,A5)')"set output '",archivo(1:len_trim(archivo)),".jpg'"
   elseif (len_trim(archivo) .eq. 9) then
     write(44,'(A12,A9,A5)')"set output '",archivo(1:len_trim(archivo)),".jpg'"
   elseif (len_trim(archivo) .eq. 10) then
     write(44,'(A12,A10,A5)')"set output '",archivo(1:len_trim(archivo)),".jpg'"
   else
     stop 'situacion no contemplada'
   endif

   write(44,'(A6)' )"replot"
   write(44,'(A12)')"set term wxt"
   write(44,'(A6)' )"replot"
   write(44,'(A11)')"unset arrow"
   write(44,'(A)'  )"#"
   write(44,'(A)'  )"#"
 enddo
 close(44)
end program gnuplotter
