subroutine logsum(log_terms, n, lsum)
!   Compute the log of a sum of terms whose logarithms are provided.

!   REQUIRED ARGUMENTS  
!     log_terms is the array (possibly multidimensional) containing the logs of the terms to be summed.

!   RETURN VALUES
!     log_sum is the log of the sum of the terms.
   implicit none
   integer, intent(IN) :: n
   double precision, intent(IN) :: log_terms(n)
   double precision, intent(OUT) :: lsum
   integer j
   double precision lt
   lsum = log_terms(1)
   do j = 2,n
      lt = log_terms(j)
      if (lsum > lt) then
         lsum = lsum + log(1.D0 + exp(-lsum + lt) )
      else
         lsum = lt + log(1.D0 + exp(lsum - lt) )
      endif
   enddo
end subroutine logsum
