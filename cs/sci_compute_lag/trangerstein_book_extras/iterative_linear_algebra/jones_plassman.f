*****************************************************************
*     perform standard incomplete choleski: column oriented
*****************************************************************
      integer function istdic(n,diag,a,ia,ja,ta,ifirst,list)

*     if the factorization was p.d. then 0 is returned
*     otherwise a negative value is returned that indicates
*     the column number where a negative diagonal was encountered

*     the order of the matrix
*     input only
      integer n
*     the diagonals of a
*     input/output
      double precision diag(*)
*     the off-diagonals of a
*     input/output
      double precision a(*)
*     pointers to the columns of a
*     ia(k) is the index in a() and ja() where column k starts
*     only the strictly lower triangle of a is stored
*     ia is length n+1 (position n+1 indicates where column n+1
*     would start if it existed)
*     input
      integer ia(*)
*     the row numbers of the off-diagonals of a
*     input
      integer ja(*)
*     a temporary work vector of length n to keep the current column
*     contents destroyed
      double precision ta(*)
*     ifirst(j) points to the next value in column j to use (length n)
*     ifirst also has a dual use.  at step k, only the first k-1 
*     elements are used for the above purpose.  for the last n-k 
*     elements, ifirst(j) indicates if if a nonzero value exists in 
*     position j of column k.
*     contents destroyed
      integer ifirst(*)
*     list(j) points to a linked list of columns that will update
*     column j (length n)
*     contents destroyed
      integer list(*)

*     variables used
      integer isk, iek, isj, iej
      integer i, j, k
      integer row
      double precision lval, t
      integer iptr
 
*****************************************
*      start of executable statements
*****************************************

      do 100 j = 1, n
        ifirst(j) = 0
        list(j) = 0
100   continue

*     loop over all columns
      do 700 k = 1,n
*        load column k into ta
         isk = ia(k)
         iek = ia(k+1)-1
         do 200 j = isk, iek
           row = ja(j)
           ta(row) = a(j)
           ifirst(row) = 1
200      continue

*        make sure the diagonal of k is okay and then take the sqrt
         if (diag(k).le.0.0d0) then
           diag(k) = -1.0d0
           goto 800
         endif
         diag(k) = sqrt(diag(k))

*        update column k using the previous columns
         j = list(k)
300      continue
         if (j.eq.0) goto 500
           isj = ifirst(j)
           iej = ia(j+1)-1
           lval = a(isj)
           isj = isj + 1
           if (isj.lt.iej) then
             ifirst(j) = isj
             iptr = j
             j = list(j)
             list(iptr) = list(ja(isj))
             list(ja(isj)) = iptr
           else
             j = list(j)
           endif
           do 400 i = isj, iej
             row = ja(i)
             if (ifirst(row).ne.0) then
               ta(row) = ta(row) - lval*a(i)
             endif
400        continue
         goto 300
500      continue

*        ifirst and list keep track of where in column k we are
         if (isk.lt.iek) then
           iptr = ja(isk)
           list(k) = list(iptr)
           list(iptr) = k
           ifirst(k) = isk
         endif

*        update remaining diagonals using column k
         lval = diag(k)
         do 600 j = isk, iek
           row = ja(j)
           t = ta(row)
           ifirst(row) = 0
           t = t/lval
           diag(row) = diag(row) - t*t
           a(j) = t
600      continue

700   continue

*     return a non-negative value
      istdic = 0
      return

*     if an error occured, return a negative value
800   continue
      istdic = -k
      return
*****************************************************************
*     end of istdic
*****************************************************************
      end
*****************************************************************
*     perform jones/plassmann incomplete choleski:column oriented
*     description of column storage format
*      only the lower triangular portion of the matrix is stored.
*      (1) the n diagonals are stored in diag
*      (2) the nz offdiagonals are stored in a
*      (3) a pointer to the beginning of each column
*          is given in ia.  ia(i) points to the
*          beginning of column i in a.  the last element
*          in this column is at position ia(i+1)-1.
*          this vector is length n+1 because the end of column
*          n must be given.
*      (4) the row numbers for the nz offdiagonals are given
*          in ja.  ja(i) is the column number for the
*          offdiagonal in a(i).
*****************************************************************
      integer function jpicc(n,diag,a,ia,ja,ta,itcol,ifirst,list)

*     if the factorization was p.d. then 0 is returned
*     otherwise a negative value is returned that indicates
*     the column number where a negative diagonal was encountered

*     the order of the matrix
*     input only
      integer n
*     the diagonals of a
*     input/output
      double precision diag(*)
*     the off-diagonals of a
*     input/output
      double precision a(*)
*     pointers to the columns of a
*     ia(k) is the index in a() and ja() where column k starts
*     only the strictly lower triangle of a is stored
*     ia is length n+1 (position n+1 indicates where column n+1
*     would start if it existed)
*     input
      integer ia(*)
*     the row numbers of the off-diagonals of a
*     input/output
      integer ja(*)
*     a temporary work vector of length n to keep the current column
*     contents destroyed
      double precision ta(*)
*     a temporary work vector of length n to keep track of the row
*     values in the current column
*     contents destroyed
      integer itcol(*)
*     ifirst(j) points to the next value in column j to use (length n)
*     ifirst also has a dual use.  at step k, only the first k-1 
*     elements are used for the above purpose.  for the last n-k 
*     elements, ifirst(j) indicates if if a nonzero value exists in 
*     position j of column k.
*     contents destroyed
      integer ifirst(*)
*     list(j) points to a linked list of columns that will update
*     column j (length n)
*     contents destroyed
      integer list(*)

*     subroutines used 
      external ibsort, dbsort, dhsort

*     variables used
      integer isj, iej, isk, iek
      integer i, j, k
      integer talen
      integer row, count
      double precision lval
      integer iptr
 
*****************************************
*      start of executable statements
*****************************************

      do 100 j = 1, n
        ifirst(j) = 0
        list(j) = 0
100    continue

*     loop over all columns
      do 900 k = 1,n
*        load column k into ta
         talen = 0
         isk = ia(k)
         iek = ia(k+1)-1
         do 200 j = isk, iek
           row = ja(j)
           ta(row) = a(j)
           talen = talen + 1
           itcol(talen) = row
           ifirst(row) = 1
200       continue

*        make sure the diagonal of k is okay and then take the sqrt
         if (diag(k).le.0.0d0) then
           diag(k) = -1.0d0
           goto 1000
         endif
         diag(k) = sqrt(diag(k))

*        update column k using the previous columns
         j = list(k)
300      continue
         if (j.eq.0) goto 500
           isj = ifirst(j)
           iej = ia(j+1)-1
           lval = a(isj)
           isj = isj + 1
           if (isj.lt.iej) then
             ifirst(j) = isj
             iptr = j
             j = list(j)
             list(iptr) = list(ja(isj))
             list(ja(isj)) = iptr
           else
             j = list(j)
           endif
           do 400 i = isj, iej
             row = ja(i)
             if (ifirst(row).ne.0) then
               ta(row) = ta(row) - lval*a(i)
             else
               ifirst(row) = 1
               talen = talen + 1
               itcol(talen) = row
               ta(row) = - lval*a(i)
             endif
400        continue
         goto 300
500      continue

*        update remaining diagonals using column k
         do 600 j = 1, talen
           row = itcol(j)
           ta(row) = ta(row)/diag(k)
           diag(row) = diag(row) - ta(row)*ta(row)
600      continue

*        find the largest elements in column k now
         count = min(iek-isk+1,talen)
         call ibsort(talen,count,ta,itcol)
         if (count.lt.20) then
           call dbsort(count,itcol)
         else
           call dhsort(count,itcol)
         endif

*        put the largest elements back into the sparse data structure
         count = 1
         do 700 j = isk, iek
           a(j) = ta(itcol(count))
           ja(j) = itcol(count)
           count = count + 1
700      continue

*        ifirst and list keep track of where in column k we are
         if (isk.lt.iek) then
           iptr = ja(isk)
           list(k) = list(iptr)
           list(iptr) = k
           ifirst(k) = isk
         endif

         do 800 j = 1, talen
           ifirst(itcol(j)) = 0
800      continue

900   continue

      jpicc = 0
      return

*     if an error occured, return a negative value
1000  continue
      jpicc = -k
      return
*****************************************************************
*     end of jpicc
*****************************************************************
      end
*****************************************************************
*     perform jones/plassmann incomplete choleski: row oriented
*     description of row storage format
*      only the lower triangular portion of the matrix is stored.
*      (1) the n diagonals are stored in rdiag
*      (2) the nz offdiagonals are stored in rofdag
*      (3) a pointer to the beginning of each row
*          is given in rrwptr.  rrwptr(i) points to the
*          beginning of row i in rofdag.  the last element
*          in this row is at position rrwptr(i+1)-1.
*          this vector is length n+1 because the end of row
*          n must be given.
*      (4) the column numbers for the nz offdiagonals are given
*          in rclnum.  rclnum(i) is the row number for the
*          offdiagonal in rofdag(i).
*****************************************************************
      integer function jpicr(n,diag,a,ia,ja,ta,itcol,iclptr,
     c iclist,irvals,ioccpy)

*     if the factorization was p.d. then 0 is returned
*     otherwise a negative value is returned that indicates
*     the row number where a negative diagonal was encountered
*     the order of the matrix
*     input only
      integer n
*     the diagonals of a
*     input/output
      double precision diag(*)
*     the off-diagonals of a
*     input/output
      double precision a(*)
*     pointers to the rows of a
*     ia(k) is the index in a() and ja() where row k starts
*     only the strictly lower triangle of a is stored
*     ia is length n+1 (position n+1 indicates where row n+1
*     would start if it existed)
*     ia is length n+1
*     input
      integer ia(*)
*     the column numbers of the off-diagonals of a
*     input/output
      integer ja(*)
*     a temporary work vector of length n to keep the current row
*     contents destroyed
      double precision ta(*)
*     a temporary work vector of length n to keep the column values in
*     the current row
*     contents destroyed
      integer itcol(*)
*     iclptr(j) points to the head of the list containing column j
*     length n
*     contents destroyed
      integer iclptr(*)
*     a temporary work vector of the same length as ja
*     iclist(i) points to the next value in the column list that 
*     element i belongs to
*     contents destroyed
      integer iclist(*)
*     a temporary work vector of the same length as ja
*     irvals(i) is the row number for element i
*     contents destroyed
      integer irvals(*)
*     a temporary work vector of length n
*     ioccpy(i) indicates whether or not an element is present in
*     column i in the current row
*     contents destroyed
      integer ioccpy(*)

*     subroutines used 
      external ibsort, dbsort, dhsort

*     variables used
      integer isk, iek
      integer i, j, k
      integer talen
      integer col, count, row
      double precision lval
      integer iptr, rptr
      integer left, mid, right
 
*****************************************
*      start of executable statements
*****************************************

      do 100 j = 1, n
        ioccpy(j) = 0
        iclptr(j) = 0
100   continue

*     loop over all columns
      do 600 k = 1,n
*        load row k into ta
         talen = 0
         isk = ia(k)
         iek = ia(k+1)-1
         do 150 j = isk, iek
           col = ja(j)
           ta(col) = a(j)
           talen = talen + 1
           itcol(talen) = col
           ioccpy(col) = 1
150      continue

         iptr = 1
200      continue
         if (iptr.le.talen) then
           col = itcol(iptr)
           lval = ta(col)/diag(col)
           ta(col) = lval
           diag(k) = diag(k) - lval*lval
           rptr = iclptr(col)
250        continue
           if (rptr.ne.0) then
             row = irvals(rptr)
             if (ioccpy(row).ne.0) then
               ta(row) = ta(row) - lval*a(rptr)
             else
               ioccpy(row) = 1
               ta(row) = - lval*a(rptr)
               if (iptr.lt.talen) then
                 if (row.lt.itcol(iptr+1)) then
                   do 300 i = talen, iptr+1, -1
                     itcol(i+1) = itcol(i)
300                continue
                   talen = talen + 1
                   itcol(iptr+1) = row
                   goto 450
                 else
                   left = iptr+1
                   if (row.gt.itcol(talen)) then
                     talen = talen + 1
                     itcol(talen) = row
                     goto 450
                   else
                     right = talen
                   endif
                 endif
               else
                 talen = talen + 1
                 itcol(talen) = row
                 goto 450
               endif
350            continue
               if (right.gt.left+1) then
                 mid = (right + left)/2
                 if (row.gt.itcol(mid)) then
                   left = mid
                 else
                   right = mid
                 endif
                 goto 350
               else
                 do 400 i = talen, right, -1
                   itcol(i+1) = itcol(i)
400              continue
                 talen = talen + 1
                 itcol(right) = row
               endif
450            continue
             endif
             rptr = iclist(rptr)
             goto 250
           endif
           iptr = iptr + 1
           goto 200
         endif

*        make sure the diagonal of k is okay and then take the sqrt
         if (diag(k).le.0.0d0) then
           diag(k) = -1.0d0
           goto 650
         endif
         diag(k) = sqrt(diag(k))

*        find the largest elements in row k now
         count = min(iek-isk+1,talen)
         call ibsort(talen,count,ta,itcol)
         if (count.lt.20) then
           call dbsort(count,itcol)
         else
           call dhsort(count,itcol)
         endif

*        put the largest elements back into the sparse data structure
         count = 1
         do 500 j = isk, iek
           a(j) = ta(itcol(count))
           ja(j) = itcol(count)
           iclist(j) = iclptr(itcol(count))
           iclptr(itcol(count)) = j
           irvals(j) = k
           count = count + 1
500      continue

         do 550 j = 1, talen
           ioccpy(itcol(j)) = 0
550      continue

600   continue

*     return a non-negative value
      jpicr = 0
      return

*     if an error occured, return a negative value
650   continue
      jpicr = -k
      return
*****************************************************************
*     end of jpicr
*****************************************************************
      end
*****************************************************************
*      get the k largest nonzeroes in akeys indirectly addressed 
*      by indvec; upon exit the first k elements in indvec will
*      contain the indices of the k largest elements in akeys
*****************************************************************
       subroutine ibsort(n,k,akeys,indvec)

*      the length of the integer vector
       integer n
*      the number wanted
       integer k
*      the double precision vector to be sorted
       double precision akeys(*)
*      the integer vector associated with akeys
*      indvec(i) gives the position in akeys of the ith element
       integer indvec(*)
 
*      the rest are internal variables
       integer i,j
       integer itemp, curptr, right, left
       double precision curmin, newval, curval, lval

       external ihsort
 
*****************************************
*      start of executable statements
*****************************************

*      if the list is small or the number required is 0 then 
*      return
       if ((n.le.1).or.(k.le.0)) return

*      heap sort the first k elements of the vector
       call ihsort(k,indvec,akeys)

*      loop through the rest of the vector and find any elements
*      that are larger than any of the first k elements
       curmin = abs(akeys(indvec(k)))
       do 400 i = k+1, n
          itemp = indvec(i)
          newval = abs(akeys(itemp))
          if (newval.gt.curmin) then
*           find position for new value
            left = 1
            lval = abs(akeys(indvec(1)))
            if (newval.gt.lval) then
              curptr = 1
              goto 200
            endif
            right = k
            curptr = (k+1)/2
100         continue            
            if (right.gt.left+1) then
              curval = abs(akeys(indvec(curptr)))
              if (curval.lt.newval) then
                right = curptr
              else
                left = curptr
                lval = curval
              endif
              curptr = (right+left)/2
              goto 100
            endif
            curptr = right

*           shift sorted values and insert new value
200         continue
            indvec(i) = indvec(k)
            do 300 j = k, curptr+1, -1
              indvec(j) = indvec(j-1)
300         continue
            indvec(curptr) = itemp
            curmin = abs(akeys(indvec(k)))
          endif
400    continue
 
       return
*****************************************************************
*      end of ibsort
*****************************************************************
       end
*****************************************************************
*      sorts an integer vector (uses bubble sort)
*      ascending order
*****************************************************************
       subroutine dbsort(n,keyvec)

*      the length of the vector
       integer n
*      the integer vector to be sorted
       integer keyvec(*)
 
*      the rest are internal variables
       integer i,j
       integer temp
 
*****************************************
*      start of executable statements
*****************************************
       do 200 i = 1, n-1
          do 100 j = i+1, n
             if (keyvec(i).gt.keyvec(j)) then
                temp = keyvec(i)
                keyvec(i) = keyvec(j)
                keyvec(j) = temp
             endif
100       continue
200    continue
 
       return
*****************************************************************
*      end of dbsort
*****************************************************************
       end
*****************************************************************
*      sorts an integer vector (uses heap sort)
*      ascending order
*****************************************************************
       subroutine dhsort(len,keys)
*      the length of the array
       integer len
*      the array to be sorted
       integer keys(*)

*      the rest are internal variables
       integer k, m, lheap, rheap, mid
       integer x

*****************************************
*      start of executable statements
*****************************************
       if (len.le.1) return

*      build the heap
       mid = len/2
       do 300 k = mid, 1, -1
          x = keys(k)
          lheap = k
          rheap = len
          m = lheap*2
100       continue
             if (m.gt.rheap) then
                keys(lheap) = x
                goto 200
             endif
             if (m.lt.rheap) then
                if (keys(m) .lt. keys(m+1)) m = m+1
             endif
             if (x.ge.keys(m)) then
                m = rheap + 1
             else
                keys(lheap) = keys(m)
                lheap = m
                m = 2*lheap
             endif
             goto 100
200       continue
300    continue

*      sort the heap
       do 600 k = len, 2, -1
          x = keys(k)
          keys(k) = keys(1)
          lheap = 1
          rheap = k-1
          m = 2
400       continue
             if (m.gt.rheap) then
                keys(lheap) = x
                goto 500
             endif
             if (m.lt.rheap) then
                if (keys(m) .lt. keys(m+1)) m = m+1
             endif
             if (x.ge.keys(m)) then
                m = rheap + 1
             else
                keys(lheap) = keys(m)
                lheap = m
                m = 2*lheap
             endif
             goto 400
500       continue
600    continue

       return
*****************************************************************
*      end of dhsort
*****************************************************************
       end
*****************************************************************
*      sorts a double precision vector indirectly addressed by
*      an integer vector(uses heap sort)
*      the indvec is rearranged such that indvec(1) addresses
*      the largest element in akeys, indvec(2) addresses the
*      next largest ....
*****************************************************************
       subroutine ihsort(len,indvec,akeys)
*      the length of the integer array
       integer len
*      the integer array that indirectly addresses the d.p. array
       integer indvec(*)
*      the array to be sorted
       double precision akeys(*)

*      the rest are internal variables
       integer k, m, lheap, rheap, mid
       integer x

*****************************************
*      start of executable statements
*****************************************
       if (len.le.1) return

*      build the heap
       mid = len/2
       do 300 k = mid, 1, -1
          x = indvec(k)
          lheap = k
          rheap = len
          m = lheap*2
100       continue
             if (m.gt.rheap) then
                indvec(lheap) = x
                goto 200
             endif
             if (m.lt.rheap) then
                if (abs(akeys(indvec(m))).gt.abs(akeys(indvec(m+1))))
     c              m = m + 1
             endif
             if (abs(akeys(x)).le.abs(akeys(indvec(m)))) then
                m = rheap + 1
             else
                indvec(lheap) = indvec(m)
                lheap = m
                m = 2*lheap
             endif
             goto 100
200       continue
300    continue

*      sort the heap
       do 600 k = len, 2, -1
          x = indvec(k)
          indvec(k) = indvec(1)
          lheap = 1
          rheap = k-1
          m = 2
400       continue
             if (m.gt.rheap) then
                indvec(lheap) = x
                goto 500
             endif
             if (m.lt.rheap) then
                if (abs(akeys(indvec(m))).gt.abs(akeys(indvec(m+1))))
     c              m = m + 1
             endif
             if (abs(akeys(x)).le.abs(akeys(indvec(m)))) then
                m = rheap + 1
             else
                indvec(lheap) = indvec(m)
                lheap = m
                m = 2*lheap
             endif
             goto 400
500       continue
600    continue

       return
*****************************************************************
*      end of ihsort
*****************************************************************
       end
