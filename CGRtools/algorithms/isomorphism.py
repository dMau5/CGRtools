# -*- coding: utf-8 -*-
#
#  Copyright 2018, 2019 Ramil Nugmanov <stsouko@live.ru>
#  This file is part of CGRtools.
#
#  CGRtools is free software; you can redistribute it and/or modify
#  it under the terms of the GNU Affero General Public License as published by
#  the Free Software Foundation; either version 3 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
#  GNU Lesser General Public License for more details.
#
#  You should have received a copy of the GNU Affero General Public License
#  along with this program; if not, see <https://www.gnu.org/licenses/>.
#
from abc import abstractmethod
from CachedMethods import cached_property
from collections import defaultdict
from itertools import permutations, product
from typing import Dict, Iterator


class Isomorphism:
    __slots__ = ()

    def __lt__(self, other):
        if len(self) >= len(other):
            return False
        return self.is_substructure(other)

    def __le__(self, other):
        return self.is_substructure(other)

    def __gt__(self, other):
        if len(self) <= len(other):
            return False
        return other.is_substructure(self)

    def __ge__(self, other):
        return other.is_substructure(self)

    def is_substructure(self, other) -> bool:
        """
        test self is substructure of other
        """
        try:
            next(self.get_mapping(other))
        except StopIteration:
            return False
        return True

    def is_equal(self, other) -> bool:
        """
        test self is structure of other
        """
        if len(self) != len(other):
            return False
        try:
            next(self.get_mapping(other))
        except StopIteration:
            return False
        return True

    @abstractmethod
    def get_mapping(self, other, *, automorphism_filter: bool = True) -> Iterator[Dict[int, int]]:
        """
        get self to other substructure mapping generator
        """
        seen = set()
        components, closures = self.__compiled_query
        for candidates in permutations(other.connected_components, len(components)):
            for match in product(*(self.__get_mapping(order, closures, other, component)
                                   for order, component in zip(components, candidates))):
                mapping = match[0]
                for m in match[1:]:
                    mapping.update(m)
                if automorphism_filter:
                    atoms = frozenset(mapping.values())
                    if atoms in seen:
                        continue
                    seen.add(atoms)
                yield mapping

    @staticmethod
    def __get_mapping(order, closures, other, atoms):
        size = len(order) - 1
        order_depth = {v[0]: k for k, v in enumerate(order)}

        o_atoms = other._atoms
        o_bonds = other._bonds

        stack = []
        path = []
        mapping = {}
        reversed_mapping = {}

        s_atom = order[0][1]
        for n, o_atom in o_atoms.items():
            if n in atoms and s_atom == o_atom:
                stack.append((n, 0))

        while stack:
            o_atom, depth = stack.pop()
            s_atom = order[depth][0]
            if depth == size:
                yield {s_atom: o_atom, **mapping}
            else:
                if len(path) != depth:
                    for x in path[depth:]:
                        del mapping[reversed_mapping[x]]
                    path = path[:depth]

                back = order[depth + 1][1]
                if back != s_atom:
                    fork = path[order_depth[back]]
                else:
                    fork = o_atom

                path.append(o_atom)
                mapping[s_atom] = o_atom
                reversed_mapping[o_atom] = s_atom

                lp = len(path)
                for o_n, o_bond in o_bonds[fork].items():
                    if o_n not in atoms:
                        continue
                    s_n, _, s_atom, s_bond = order[lp]
                    if o_n not in path and s_bond == o_bond and s_atom == o_atoms[o_n] \
                            and all(bond == o_bonds[mapping[m]].get(o_n) for m, bond in closures[s_n]):
                        stack.append((o_n, lp))

    @cached_property
    def __compiled_query(self):
        atoms = self._atoms
        bonds = self._bonds
        atoms_order = self.atoms_order

        closures = defaultdict(list)
        components = []
        seen = set()
        while len(seen) < len(atoms):
            start = max(atoms.keys() - seen, key=lambda x: atoms_order[x])
            seen.add(start)
            stack = [(n, start, atoms[n], bond) for n, bond in sorted(bonds[start].items(),
                                                                      key=lambda x: atoms_order[x[0]])]
            order = [(start, atoms[start])]
            components.append(order)

            while stack:
                front, back, *_ = atom = stack.pop()
                if front not in seen:
                    order.append(atom)
                    for n, bond in sorted(bonds[front].items(), key=lambda x: atoms_order[x[0]]):
                        if n != back:
                            if n not in seen:
                                stack.append((n, front, atoms[n], bond))
                            else:
                                closures[front].append((n, bond))
                    seen.add(front)
        return components, closures


__all__ = ['Isomorphism']
