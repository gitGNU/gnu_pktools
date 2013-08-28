/**********************************************************************
IndexValue.h: class to work with structs containing an index and a value
Copyright (C) 2008-2012 Pieter Kempeneers

This file is part of pktools

pktools is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

pktools is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with pktools.  If not, see <http://www.gnu.org/licenses/>.
***********************************************************************/
#ifndef _INDEXVALUE_H_
#define _INDEXVALUE_H_

struct IndexValue{
  int position;
  double value;
};
class Compare_IndexValue{
public:
  int operator() (const IndexValue& pv1, const IndexValue& pv2) const{
    return pv1.value>pv2.value;//for decreasing order
  }
};
class Decrease_IndexValue{
public:
  int operator() (const IndexValue& pv1, const IndexValue& pv2) const{
    return pv1.value>pv2.value;//for decreasing order
  }
};
class Increase_IndexValue{
public:
  int operator() (const IndexValue& pv1, const IndexValue& pv2) const{
    return pv1.value<pv2.value;//for increasing order
  }
};
#endif /* _INDEXVALUE_H_ */
