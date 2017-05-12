! Function to get C6 coefficients
!
! Copyright (C) 2016, Bálint Aradi
! Copyright (C) 2017, Anders Steen Christensen
!
! This program is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 1, or (at your option)
! any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
! GNU General Public License for more details.
!
! For the GNU General Public License, see <http://www.gnu.org/licenses/>
!

function fgetc6_periodic(coords, nAtoms, atnum, latvecs) result(c6)

  use dftd3_api
  implicit none

  ! Working precision

  ! Same conversion factor as in dftd3
  double precision, parameter :: AA__Bohr = 1.0d0 / 0.52917726d0

  integer, intent(in) :: nAtoms
  
  double precision, intent(in) :: coords(:, :)
  
  ! integer, intent(in) :: species(nAtoms)
  integer, intent(in) :: atnum(:)

  double precision, dimension(nAtoms, nAtoms) :: c6

  double precision, dimension(3, 3):: latVecs

  ! ! Lattice vectors in Angstrom as found in dna.xyz/dna.poscar
  ! ! They must be converted to Bohr before passed to dftd3
  ! real(wp), parameter :: latVecs(3, 3) = reshape([&
  !      &  8.0000000000E+00,   0.0000000000E+00,   0.0000000000E+00, &
  !      &  0.0000000000E+00,   8.0000000000E+00,   0.0000000000E+00, &
  !      &  0.0000000000E+00,   0.0000000000E+00,   1.5000000000E+01  &
  !      & ] * AA__Bohr, [3, 3])

  ! integer, parameter :: nSpecies = 4
  ! character(2), parameter :: speciesNames(nSpecies) = [ 'N ', 'C ', 'O ', 'H ']
  

  type(dftd3_input) :: input
  type(dftd3_calc) :: dftd3
  double precision :: edisp
  double precision :: grads(3, nAtoms), stress(3, 3)

  ! Initialize dftd3
  call dftd3_init(dftd3, input)
  call dftd3_set_functional(dftd3, func='dftb3', version=4, tz=.false.)

  call get_c6_pbc(dftd3, coords*AA__Bohr, atnum, latvecs*AA__Bohr, edisp, grads, stress, c6)

end function fgetc6_periodic

function fgetc6(coords, nAtoms, atnum) result(c6)

  use dftd3_api
  implicit none

  ! Working precision

  ! Same conversion factor as in dftd3
  double precision, parameter :: AA__Bohr = 1.0d0 / 0.52917726d0

  integer, intent(in) :: nAtoms
  
  double precision, intent(in) :: coords(:, :)
  
  ! integer, intent(in) :: species(nAtoms)
  integer, intent(in) :: atnum(:)

  double precision, dimension(nAtoms, nAtoms) :: c6

  ! ! Lattice vectors in Angstrom as found in dna.xyz/dna.poscar
  ! ! They must be converted to Bohr before passed to dftd3
  ! real(wp), parameter :: latVecs(3, 3) = reshape([&
  !      &  8.0000000000E+00,   0.0000000000E+00,   0.0000000000E+00, &
  !      &  0.0000000000E+00,   8.0000000000E+00,   0.0000000000E+00, &
  !      &  0.0000000000E+00,   0.0000000000E+00,   1.5000000000E+01  &
  !      & ] * AA__Bohr, [3, 3])

  ! integer, parameter :: nSpecies = 4
  ! character(2), parameter :: speciesNames(nSpecies) = [ 'N ', 'C ', 'O ', 'H ']
  

  type(dftd3_input) :: input
  type(dftd3_calc) :: dftd3
  double precision :: edisp
  double precision :: grads(3, nAtoms)

  ! Initialize dftd3
  call dftd3_init(dftd3, input)
  call dftd3_set_functional(dftd3, func='dftb3', version=4, tz=.false.)

  call get_c6(dftd3, coords*AA__Bohr, atnum, edisp, grads, c6)

end function fgetc6


