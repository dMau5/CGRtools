# -*- coding: utf-8 -*-
#
#  Copyright 2017-2019 Ramil Nugmanov <stsouko@live.ru>
#  This file is part of CGRtools.
#
#  CGRtools is free software; you can redistribute it and/or modify
#  it under the terms of the GNU Lesser General Public License as published by
#  the Free Software Foundation; either version 3 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
#  GNU Lesser General Public License for more details.
#
#  You should have received a copy of the GNU Lesser General Public License
#  along with this program; if not, see <https://www.gnu.org/licenses/>.
#


class EmptyMolecule(ValueError):
    """
    bad files parsing
    """


class MappingError(ValueError):
    """
    bad files parsing
    """


class InvalidFileType(TypeError):
    """
    bad files parsing
    """


class InvalidAromaticRing(ValueError):
    """
    impossible kekule structure
    """


class InvalidAtomNumber(KeyError):
    """
    atom number not found
    """


class InvalidWedgeMark(ValueError):
    """
    wedge mark error
    """


class InvalidStereoCenter(KeyError):
    """
    impossible stereo
    """
