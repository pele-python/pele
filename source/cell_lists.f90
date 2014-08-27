!   GMIN: A program for finding global minima
!   Copyright (C) 1999-2006 David J. Wales
!   This file is part of GMIN.
!
!   GMIN is free software; you can redistribute it and/or modify
!   it under the terms of the GNU General Public License as published by
!   the Free Software Foundation; either version 2 of the License, or
!   (at your option) any later version.
!
!   GMIN is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!   GNU General Public License for more details.
!
!   You should have received a copy of the GNU General Public License
!   along with this program; if not, write to the Free Software
!   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
!
!

module cell_lists_mod
   !a module for creating cell lists using linked-lists
   ! http://en.wikipedia.org/wiki/Cell_lists
   !currently only works for cubic periodic systems, but could be easily
   !modified
   implicit none
   private

   double precision rcut  !potential cutoff
   double precision rcell !lenth of a cell
   double precision boxl  !box size (assume cubic)
   integer ncellx !number of cells in the x direction (assume cubic box)
   integer ncellx2
   integer ncells !total number of cells
   integer natoms 

   integer, allocatable :: hoc(:) !head of cell:  hoc(icell) is the first atom in cell icell
   integer, allocatable :: ll(:)  !linked list:    ll(i) is the next atom after i in the cell
                                  !             if ll(i)==0 then there are no more atoms in the cell

   integer, allocatable :: cell_neib(:,:)  ! a list containing neighboring cells
   integer :: ncell_neib                   ! the length of cell_neib

   integer, allocatable :: cell_lists_neib_list(:,:)  ! a list containing neighboring atoms
   integer cell_lists_neib_nlist


   public :: cell_lists_setup, make_neighbor_list, cell_lists_neib_list, cell_lists_neib_nlist

   contains

   subroutine cell_lists_setup(natoms_i, boxl_i, rcut_i)
      !allocate space, initialize values
      implicit none
      integer, intent(in) :: natoms_i
      double precision, intent(in) :: boxl_i, rcut_i

      rcut = rcut_i
      natoms = natoms_i
      boxl = boxl_i

      !determine the number of cells in the linear direction.  This is arbitrary,
      !and should be chosen based on what is fastest.  rcell should not be larger
      !than rcut, but could be smaller.
      ncellx = int( boxl / rcut ) * 2
      if (ncellx .lt. 8) ncellx = 8
      ncellx2 = ncellx**2

      rcell = boxl / ncellx !rcell must be greater than rcut
      ncells = ncellx**3

      write(*,*) "cell_lists_mod> rcut   ", rcut
      write(*,*) "cell_lists_mod> rcell  ", rcell
      write(*,*) "cell_lists_mod> ncellx ", ncellx
      write(*,*) "cell_lists_mod> ncells ", ncells

      allocate(hoc( ncells ) )
      allocate(ll( natoms ) )
      allocate(cell_neib( 2, ncells*(ncells+1)/2 ) ) !this allocates way too much memory
      allocate(cell_lists_neib_list( 2, natoms*(natoms-1)/2 ) ) !this allocates way too much memory

      call setup_cell_neibs()

   end subroutine cell_lists_setup

   subroutine xyz2icell( xyz_i, icell )
      !determine the cell index from the xyz coordinates
      implicit none
      double precision, intent(in) :: xyz_i(3)
      integer, intent(out) :: icell
      double precision :: xyz(3)
      integer ijk(3), i
      xyz(:) = xyz_i(:)
      !put x, y, z in box [0,boxl)
      do i=1,3
         xyz(i) = xyz(i) - boxl * nint( xyz(i) / boxl )
         if (xyz(i) .lt. 0) xyz(i) = xyz(i) + boxl
         ijk(i) = floor(xyz(i) / rcell)
      enddo
      icell = ijk(1) + ijk(2)*ncellx + ijk(3) * ncellx2 +1 !+1 for fortran indexing
   end subroutine xyz2icell

   subroutine icell2xyz( icell_i, xyz )
      !return the xyz coordinates of the corner of the cell
      implicit none
      double precision, intent(out) :: xyz(3)
      integer, intent(in) :: icell_i
      integer icell, ijk(3)
      icell = icell_i - 1

      ijk(3) = int( (icell) / ncellx2 ) 
      ijk(2) = int( (icell - ijk(3)*ncellx2 ) / ncellx ) 
      ijk(1) = int( (icell - ijk(3)*ncellx2 - ijk(2)*ncellx ) ) 

      xyz(3) = rcell * ijk(3)
      xyz(2) = rcell * ijk(2)
      xyz(1) = rcell * ijk(1)
      !write(*,*) "icell2xyz", icell_i, xyz, ijk, ncellx2, icell/ncellx2
   end subroutine icell2xyz

   function are_cell_neibs(icell, jcell) result( b )
      !determine if icell and jcell are close enough together
      !find the minimum distance between the two cells
      implicit none
      double precision :: xyz1(3), xyz2(3), r2, x, dxmin
      logical b
      integer i,k, icell, jcell
      call icell2xyz( icell, xyz1 )
      call icell2xyz( jcell, xyz2 )
      xyz2 = xyz2(:) - xyz1(:)
      !write(*,*) "xyz2", xyz2(:)
      do i=1,3 !loop over x,y,z
         dxmin = 1000.d0
         do k=-1,1 !determine minimum distance in this direction
            x = xyz2(i) + dble(k) * rcell
            x = x - boxl*nint(x/boxl)
            if (abs(x) .lt. abs(dxmin)) then
               dxmin = x
               !write(*,*) "dxmin k", dxmin, k
            endif
         enddo
         xyz2(i) = dxmin
         !write(*,*) "xyz2(i), i", xyz2(i), i
      enddo
      !write(*,*) "xyz2", xyz2(:)
      r2 = sum( xyz2(:)**2 )
      !write(*,*) r2, icell, jcell
      b = (r2 .lt. rcut**2)
      !if ( b) write(*,*) sqrt(r2), rcut, icell, jcell
   end function are_cell_neibs


   subroutine setup_cell_neibs()
      !determine which cells are neighbors.
      !a cell is it's own neighbor
      implicit none
      integer n, icell, jcell
      n = 0
      do icell = 1,ncells
         n = n + 1
         !write(*,*) "cell neibs", icell, icell, n
         cell_neib(1,n) = icell
         cell_neib(2,n) = icell
      enddo
      do icell = 1,ncells
         do jcell = icell+1,ncells
            if ( are_cell_neibs(icell, jcell) ) then
               n = n + 1
               !write(*,*) "cell neibs", icell, jcell, n
               cell_neib(1,n) = icell
               cell_neib(2,n) = jcell
            !else
               !write(*,*) "not cell neibs", icell, jcell, n
            endif
         enddo
      enddo
      ncell_neib = n
      write(*,*) "number of cell-cell neighbors", ncell_neib, ncells*14
   end subroutine setup_cell_neibs

   subroutine generate_cell_lists(coords)
      !determine which cell each atom is in
      !populate the arrays hoc and ll
      implicit none
      double precision, intent(in) :: coords(3*natoms)
      integer j, icell
      hoc(:) = 0
      ll(:) = 0
      do j=1,natoms
         call xyz2icell( coords(3*(j-1) + 1 : 3*(j-1) + 3), icell )
         ll(j) = hoc(icell)
         hoc(icell) = j
      enddo
   end subroutine generate_cell_lists

   subroutine make_neighbor_list(coords)
      !convert cell lists to neighbor lists
      !this is a waste of time, but it puts the output in a form that's easy to
      !use
      implicit none
      double precision, intent(in) :: coords(3*natoms)
      integer :: i,j, icell, jcell, n, nlist
      !double precision :: time0, time1, time2, time01 = 0.d0, time12 = 0.d0

      !CALL MYCPU_TIME(TIME0)
      call generate_cell_lists(coords)
      !CALL MYCPU_TIME(TIME1)

      nlist = 0
      do n = 1,ncell_neib
         icell = cell_neib(1,n)
         jcell = cell_neib(2,n)
         if (icell .ne. jcell) then
            !loop over all atoms in each cell
            i = hoc( icell )
            do while (i .ne. 0)
               j = hoc( jcell )
               do while (j .ne. 0)
                  nlist = nlist + 1
                  cell_lists_neib_list(1,nlist) = i
                  cell_lists_neib_list(2,nlist) = j
                  j = ll(j)
               enddo
               i = ll(i)
            enddo
         else
            !loop over all atoms in each cell, avoiding duplicates
            i = hoc( icell )
            do while (i .ne. 0)
               j = hoc( jcell )
               do while (j .ne. i)
                  nlist = nlist + 1
                  cell_lists_neib_list(1,nlist) = i
                  cell_lists_neib_list(2,nlist) = j
                  j = ll(j)
               enddo
               i = ll(i)
            enddo
         endif
      enddo
      cell_lists_neib_nlist = nlist

      !CALL MYCPU_TIME(TIME2)
      !time01 = time01 + time1 - time0
      !time12 = time12 + time2 - time1
   end subroutine make_neighbor_list



end module cell_lists_mod
