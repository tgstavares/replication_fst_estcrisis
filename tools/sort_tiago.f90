module sort_tiago
  use prec
contains

  subroutine qsort_dble(n,v,idx)
    implicit none
    integer, intent(in)::n
    real(dp), intent(in)::v(n)
    integer, intent(out)::idx(n)
    integer iii
    real(dp) a(n)
    a = v
    !write(*,'(5f12.4)')a
    idx = (/(dble(iii),iii=1,n)/)
    call qsort(a,idx)
    !write(*,'(5f12.4)')a
    !write(*,'(5i12)')idx    
  end subroutine qsort_dble
  
  subroutine sort_dble(n,v,idx)
    implicit none
    integer, intent(in)::n
    real(dp), intent(in)::v(n)
    integer, intent(out)::idx(n)

    real(dp) a(n)

    integer ( kind = 4 ) i
    integer ( kind = 4 ) i_save
    integer ( kind = 4 ) indx
    integer ( kind = 4 ) isgn
    integer ( kind = 4 ) j
    integer ( kind = 4 ) j_save
    integer ( kind = 4 ) k
    integer ( kind = 4 ) kk
    integer ( kind = 4 ) k_save
    integer ( kind = 4 ) l_save
    integer ( kind = 4 ) n_save

    integer iii
    
    a = v
    idx = (/(dble(iii),iii=1,n)/)

    write(*,'(5f12.4)')a
    
    do

       call sort_safe_rc ( n, indx, i, j, isgn, &
            i_save, j_save, k_save, l_save, n_save )

       if ( indx < 0 ) then
          isgn = 1
          if ( a(i) <= a(j) ) then
             isgn = -1
          end if
       else if ( 0 < indx ) then
          k    = a(i)
          a(i) = a(j)
          a(j) = k

          kk     = idx(i)
          idx(i) = idx(j)
          idx(j) = kk
       else
          exit
       end if

    end do

    write(*,'(5f12.4)')a
    write(*,'(5i12)')idx

  end subroutine sort_dble

  subroutine sort_rc ( n, indx, i, j, isgn )

    !*****************************************************************************80
    !
    !! SORT_RC externally sorts a list of items into ascending order.
    !
    !  Discussion:
    !
    !    The actual list of data is not passed to the routine.  Hence this
    !    routine may be used to sort integers, reals, numbers, names,
    !    dates, shoe sizes, and so on.  After each call, the routine asks
    !    the user to compare or interchange two items, until a special
    !    return value signals that the sorting is completed.
    !
    !    Note that this function uses internal persistent memory during the sort.
    !
    !  Example:
    !
    !    n = 100
    !    indx = 0
    !    i = 0
    !    j = 0
    !    isgn = 0
    !
    !    do
    !
    !      call sort_rc ( n, indx, i, j, isgn )
    !
    !      if ( indx < 0 ) then
    !
    !        isgn = 1
    !        if ( a(i) <= a(j) ) then
    !          isgn = -1
    !        end if
    !
    !      else if ( 0 < indx ) then
    !
    !        k    = a(i)
    !        a(i) = a(j)
    !        a(j) = k
    !
    !      else
    !
    !        exit
    !
    !      end if
    !
    !    end do
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    05 February 2004
    !
    !  Author:
    !
    !    Original FORTRAN77 version by Albert Nijenhuis, Herbert Wilf.
    !    FORTRAN90 version by John Burkardt.
    !
    !  Reference:
    !
    !    Albert Nijenhuis, Herbert Wilf,
    !    Combinatorial Algorithms for Computers and Calculators,
    !    Academic Press, 1978,
    !    ISBN: 0-12-519260-6,
    !    LC: QA164.N54.
    !
    !  Parameters:
    !
    !    Input, integer ( kind = 4 ) N, the number of items to be sorted.
    !
    !    Input/output, integer ( kind = 4 ) INDX, the main communication signal.
    !    The user must set INDX to 0 before the first call.
    !    Thereafter, the user should not change the value of INDX until
    !    the sorting is done.
    !    On return, if INDX is
    !    * greater than 0,
    !      > interchange items I and J;
    !      > call again.
    !    * less than 0,
    !      > compare items I and J;
    !      > set ISGN = -1 if I < J, ISGN = +1 if J < I;
    !      > call again.
    !    * equal to 0, the sorting is done.
    !
    !    Output, integer ( kind = 4 ) I, J, the indices of two items.
    !    On return with INDX positive, elements I and J should be interchanged.
    !    On return with INDX negative, elements I and J should be compared, and
    !    the result reported in ISGN on the next call.
    !
    !    Input, integer ( kind = 4 ) ISGN, results of comparison of elements
    !    I and J. (Used only when the previous call returned INDX less than 0).
    !    ISGN <= 0 means I is less than or equal to J;
    !    0 <= ISGN means I is greater than or equal to J.
    !
    implicit none

    integer ( kind = 4 ) i
    integer ( kind = 4 ), save :: i_save = 0
    integer ( kind = 4 ) indx
    integer ( kind = 4 ) isgn
    integer ( kind = 4 ) j
    integer ( kind = 4 ), save :: j_save = 0
    integer ( kind = 4 ), save :: k_save = 0
    integer ( kind = 4 ), save :: l_save = 0
    integer ( kind = 4 ) n
    integer ( kind = 4 ), save :: n_save = 0
    !
    !  INDX = 0: This is the first call.
    !
    if ( indx == 0 ) then

       i_save = 0
       j_save = 0
       k_save = n / 2
       l_save = n / 2
       n_save = n
       !
       !  INDX < 0: The user is returning the results of a comparison.
       !
    else if ( indx < 0 ) then

       if ( indx == -2 ) then

          if ( isgn < 0 ) then
             i_save = i_save + 1
          end if

          j_save = l_save
          l_save = i_save
          indx = -1
          i = i_save
          j = j_save
          return

       end if

       if ( 0 < isgn ) then
          indx = 2
          i = i_save
          j = j_save
          return
       end if

       if ( k_save <= 1 ) then

          if ( n_save == 1 ) then
             i_save = 0
             j_save = 0
             indx = 0
          else
             i_save = n_save
             n_save = n_save - 1
             j_save = 1
             indx = 1
          end if

          i = i_save
          j = j_save
          return

       end if

       k_save = k_save - 1
       l_save = k_save
       !
       !  0 < INDX, the user was asked to make an interchange.
       !
    else if ( indx == 1 ) then

       l_save = k_save

    end if

    do

       i_save = 2 * l_save

       if ( i_save == n_save ) then
          j_save = l_save
          l_save = i_save
          indx = -1
          i = i_save
          j = j_save
          return
       else if ( i_save <= n_save ) then
          j_save = i_save + 1
          indx = -2
          i = i_save
          j = j_save
          return
       end if

       if ( k_save <= 1 ) then
          exit
       end if

       k_save = k_save - 1
       l_save = k_save

    end do

    if ( n_save == 1 ) then
       i_save = 0
       j_save = 0
       indx = 0
       i = i_save
       j = j_save
    else
       i_save = n_save
       n_save = n_save - 1
       j_save = 1
       indx = 1
       i = i_save
       j = j_save
    end if

    return
  end subroutine sort_rc

  subroutine sort_safe_rc ( n, indx, i, j, isgn, i_save, j_save, &
       k_save, l_save, n_save )

    !*****************************************************************************80
    !
    !! SORT_SAFE_RC externally ascending sorts a list of items.
    !
    !  Discussion:
    !
    !    This is a version of SORT_RC which does not rely on
    !    storing certain work variables internally to the function.  This makes
    !    the function somewhat more awkward to call, but easier to program
    !    in a variety of languages, and safe to use in a parallel programming
    !    environment, or in cases where the sorting of several vectors is to
    !    be carried out at more or less the same time.
    !
    !    The actual list of data is not passed to the routine.  Hence this
    !    routine may be used to sort integers, reals, numbers, names,
    !    dates, shoe sizes, and so on.  After each call, the routine asks
    !    the user to compare or interchange two items, until a special
    !    return value signals that the sorting is completed.
    !
    !  Example:
    !
    !    n = 100
    !    indx = 0
    !
    !    do
    !
    !      call sort_safe_rc ( n, indx, i, j, isgn, i_save, j_save,
    !        k_save, l_save, n_save )
    !
    !      if ( indx < 0 ) then
    !
    !        isgn = 1
    !        if ( a(i) <= a(j) ) then
    !          isgn = -1
    !        end if
    !
    !      else if ( 0 < indx ) then
    !
    !        k    = a(i)
    !        a(i) = a(j)
    !        a(j) = k
    !
    !      else
    !
    !        exit
    !
    !      end if
    !
    !    end do
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    09 March 2015
    !
    !  Author:
    !
    !    John Burkardt.
    !
    !  Reference:
    !
    !    Albert Nijenhuis, Herbert Wilf,
    !    Combinatorial Algorithms for Computers and Calculators,
    !    Academic Press, 1978,
    !    ISBN: 0-12-519260-6,
    !    LC: QA164.N54.
    !
    !  Parameters:
    !
    !    Input, integer ( kind = 4 ) N, the number of items to be sorted.
    !
    !    Input/output, integer ( kind = 4 ) INDX, the main communication signal.
    !    The user must set INDX to 0 before the first call.
    !    Thereafter, the user should not change the value of INDX until
    !    the sorting is done.
    !    On return, if INDX is
    !    * greater than 0,
    !      > interchange items I and J;
    !      > call again.
    !    * less than 0,
    !      > compare items I and J;
    !      > set ISGN = -1 if I < J, ISGN = +1 if J < I;
    !      > call again.
    !    * equal to 0, the sorting is done.
    !
    !    Output, integer ( kind = 4 ) I, J, the indices of two items.
    !    On return with INDX positive, elements I and J should be interchanged.
    !    On return with INDX negative, elements I and J should be compared, and
    !    the result reported in ISGN on the next call.
    !
    !    Input, integer ( kind = 4 ) ISGN, results of comparison of elements
    !    I and J. (Used only when the previous call returned INDX less than 0).
    !    ISGN <= 0 means I is less than or equal to J;
    !    0 <= ISGN means I is greater than or equal to J.
    !
    !    Input/output, integer ( kind = 4 ) I_SAVE, J_SAVE, K_SAVE, L_SAVE,
    !    N_SAVE, workspace needed by the routine.  Before calling the function,
    !    the user should declare variables to hold these values, but should
    !    not change them, and need not ever examine them.
    !
    implicit none

    integer ( kind = 4 ) i
    integer ( kind = 4 ) i_save
    integer ( kind = 4 ) indx
    integer ( kind = 4 ) isgn
    integer ( kind = 4 ) j
    integer ( kind = 4 ) j_save
    integer ( kind = 4 ) k_save
    integer ( kind = 4 ) l_save
    integer ( kind = 4 ) n
    integer ( kind = 4 ) n_save
    !
    !  INDX = 0: This is the first call.
    !
    if ( indx == 0 ) then

       i_save = 0
       j_save = 0
       k_save = n / 2
       l_save = n / 2
       n_save = n
       !
       !  INDX < 0: The user is returning the results of a comparison.
       !
    else if ( indx < 0 ) then

       if ( indx == -2 ) then

          if ( isgn < 0 ) then
             i_save = i_save + 1
          end if

          j_save = l_save
          l_save = i_save
          indx = -1
          i = i_save
          j = j_save
          return

       end if

       if ( 0 < isgn ) then
          indx = 2
          i = i_save
          j = j_save
          return
       end if

       if ( k_save <= 1 ) then

          if ( n_save == 1 ) then
             i_save = 0
             j_save = 0
             indx = 0
          else
             i_save = n_save
             n_save = n_save - 1
             j_save = 1
             indx = 1
          end if

          i = i_save
          j = j_save
          return

       end if

       k_save = k_save - 1
       l_save = k_save
       !
       !  0 < INDX, the user was asked to make an interchange.
       !
    else if ( indx == 1 ) then

       l_save = k_save

    end if

    do

       i_save = 2 * l_save

       if ( i_save == n_save ) then
          j_save = l_save
          l_save = i_save
          indx = -1
          i = i_save
          j = j_save
          return
       else if ( i_save <= n_save ) then
          j_save = i_save + 1
          indx = -2
          i = i_save
          j = j_save
          return
       end if

       if ( k_save <= 1 ) then
          exit
       end if

       k_save = k_save - 1
       l_save = k_save

    end do

    if ( n_save == 1 ) then
       i_save = 0
       j_save = 0
       indx = 0
       i = i_save
       j = j_save
    else
       i_save = n_save
       n_save = n_save - 1
       j_save = 1
       indx = 1
       i = i_save
       j = j_save
    end if

    return
  end subroutine sort_safe_rc

  ! ***********************************
  ! *
  Subroutine Qsort(X, Ipt)
    ! *
    ! ***********************************
    ! * Sort Array X(:) in ascendent order 
    ! * If present Ipt, a pointer with the 
    ! * changes is returned in Ipt.
    ! ***********************************

    Type Limits
       Integer :: Ileft, Iright
    End Type Limits

    ! For a list with Isw number of elements or
    ! less use Insrt
    Integer, Parameter :: Isw = 10

    Real (kind=8), Intent (inout) :: X(:)
    Integer, Intent (out), Optional :: Ipt(:)

    Integer :: I, Ipvn, Ileft, Iright, ISpos, ISmax
    Integer, Allocatable :: IIpt(:)
    Type (Limits), Allocatable :: Stack(:)


    Allocate(Stack(Size(X)))

    Stack(:)%Ileft = 0
    If (Present(Ipt)) Then
       Forall (I=1:Size(Ipt)) Ipt(I) = I

       ! Iniitialize the stack
       Ispos = 1
       Ismax = 1
       Stack(ISpos)%Ileft  = 1
       Stack(ISpos)%Iright = Size(X)

       Do While (Stack(ISpos)%Ileft /= 0)

          Ileft = Stack(ISPos)%Ileft
          Iright = Stack(ISPos)%Iright
          If (Iright-Ileft <= Isw) Then
             CALL InsrtLC(X, Ipt, Ileft,Iright)
             ISpos = ISPos + 1
          Else
             Ipvn = ChoosePiv(X, Ileft, Iright)
             Ipvn = Partition(X, Ileft, Iright, Ipvn, Ipt)

             Stack(ISmax+1)%Ileft = Ileft
             Stack(ISmax+1) %Iright = Ipvn-1
             Stack(ISmax+2)%Ileft = Ipvn + 1
             Stack(ISmax+2)%Iright = Iright
             ISpos = ISpos + 1
             ISmax = ISmax + 2
          End If
       End Do

    Else

       ! Iniitialize the stack
       Ispos = 1
       Ismax = 1
       Stack(ISpos)%Ileft  = 1
       Stack(ISpos)%Iright = Size(X)

       Allocate(IIpt(10))
       Do While (Stack(ISpos)%Ileft /= 0)
          !          Write(*,*)Ispos, ISmax

          Ileft = Stack(ISPos)%Ileft
          Iright = Stack(ISPos)%Iright
          If (Iright-Ileft <= Isw) Then
             CALL InsrtLC(X, IIpt, Ileft, Iright)
             ISpos = ISPos + 1
          Else
             Ipvn = ChoosePiv(X, Ileft, Iright)
             Ipvn = Partition(X, Ileft, Iright, Ipvn)

             Stack(ISmax+1)%Ileft = Ileft
             Stack(ISmax+1) %Iright = Ipvn-1
             Stack(ISmax+2)%Ileft = Ipvn + 1
             Stack(ISmax+2)%Iright = Iright
             ISpos = ISpos + 1
             ISmax = ISmax + 2
          End If
       End Do
       Deallocate(IIpt)

    End If

    Deallocate(Stack)

    Return

  CONTAINS

    ! ***********************************
    Integer Function ChoosePiv(XX, IIleft, IIright) Result (IIpv)
      ! ***********************************
      ! * Choose a Pivot element from XX(Ileft:Iright)
      ! * for Qsort. This routine chooses the median
      ! * of the first, last and mid element of the 
      ! * list.
      ! ***********************************

      Real (kind=8), Intent (in) :: XX(:)
      Integer, Intent (in) :: IIleft, IIright

      Real (kind=8) :: XXcp(3)
      Integer :: IIpt(3), IImd

      IImd = Int((IIleft+IIright)/2)
      XXcp(1) = XX(IIleft)
      XXcp(2) = XX(IImd)
      XXcp(3) = XX(IIright)
      IIpt = (/1,2,3/)

      CALL InsrtLC(XXcp, IIpt, 1, 3)

      Select Case (IIpt(2))
      Case (1)
         IIpv = IIleft
      Case (2)
         IIpv = IImd
      Case (3)
         IIpv = IIright
      End Select

      Return
    End Function ChoosePiv

    ! ***********************************
    Subroutine InsrtLC(XX, IIpt, IIl, IIr)
      ! ***********************************
      ! * Perform an insertion sort of the list 
      ! * XX(:) between index values IIl and IIr.
      ! * IIpt(:) returns the permutations
      ! * made to sort.
      ! ***********************************

      Real (kind=8), Intent (inout) :: XX(:)
      Integer, Intent (inout) :: IIpt(:)
      Integer, Intent (in) :: IIl, IIr

      Real (kind=8) :: RRtmp
      Integer :: II, JJ, IItmp

      Do II = IIl+1, IIr
         RRtmp = XX(II)
         Do JJ = II-1, 1, -1
            If (RRtmp < XX(JJ)) Then
               XX(JJ+1) = XX(JJ)
               CALL Swap_IN(IIpt, JJ, JJ+1)
            Else
               Exit
            End If
         End Do
         XX(JJ+1) = RRtmp
      End Do

      Return
    End Subroutine InsrtLC

  End Subroutine Qsort

  ! ***********************************
  ! *
  Integer Function Partition(X, Ileft, Iright, Ipv, Ipt) Result (Ipvfn)
    ! *
    ! ***********************************
    ! * This routine arranges the array X
    ! * between the index values Ileft and Iright
    ! * positioning elements smallers than
    ! * X(Ipv) at the left and the others 
    ! * at the right.
    ! * Internal routine used by Qsort.
    ! ***********************************

    Real (kind=8), Intent (inout) :: X(:)
    Integer, Intent (in) :: Ileft, Iright, Ipv
    Integer, Intent (inout), Optional :: Ipt(:)

    Real (kind=8) :: Rpv
    Integer :: I

    Rpv = X(Ipv)
    CALL Swap(X, Ipv, Iright)
    If (Present(Ipt)) CALL Swap_IN(Ipt, Ipv, Iright)
    Ipvfn = Ileft

    If (Present(Ipt))  Then
       Do I = Ileft, Iright-1
          If (X(I) <= Rpv) Then
             CALL Swap(X, I, Ipvfn)
             CALL Swap_IN(Ipt, I, Ipvfn)
             Ipvfn = Ipvfn + 1
          End If
       End Do
    Else
       Do I = Ileft, Iright-1
          If (X(I) <= Rpv) Then
             CALL Swap(X, I, Ipvfn)
             Ipvfn = Ipvfn + 1
          End If
       End Do
    End If

    CALL Swap(X, Ipvfn, Iright)
    If (Present(Ipt)) CALL Swap_IN(Ipt, Ipvfn, Iright)

    Return
  End Function Partition

  ! ***********************************
  ! *
  Subroutine Swap(X, I, J)
    ! *
    ! ***********************************
    ! * Swaps elements I and J of array X(:). 
    ! ***********************************

    Real (kind=8), Intent (inout) :: X(:)
    Integer, Intent (in) :: I, J

    Real (kind=8) :: Itmp

    Itmp = X(I)
    X(I) = X(J)
    X(J) = Itmp

    Return
  End Subroutine Swap

  ! ***********************************
  ! *
  Subroutine Swap_IN(X, I, J)
    ! *
    ! ***********************************
    ! * Swaps elements I and J of array X(:). 
    ! ***********************************

    Integer, Intent (inout) :: X(:)
    Integer, Intent (in) :: I, J

    Integer :: Itmp

    Itmp = X(I)
    X(I) = X(J)
    X(J) = Itmp

    Return
  End Subroutine Swap_IN

end module sort_tiago
