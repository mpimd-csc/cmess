//
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program; if not, see <http://www.gnu.org/licenses/>.
//
// Copyright (C) Peter Benner, Martin Koehler, Jens Saak and others
//               2009-2018
//

#ifndef MESS_OPTIONS_OCTAVE_H_
#define MESS_OPTIONS_OCTAVE_H_

#include <octave/oct.h>
#include "mess/mess.h"


class MESS_options : public octave_base_value {

    public:

        /*-----------------------------------------------------------------------------
         *  constructors
         *-----------------------------------------------------------------------------*/
        MESS_options(void);
        MESS_options(mess_options opt);
        MESS_options(const MESS_options & opt);


        /*-----------------------------------------------------------------------------
         *  getter
         *-----------------------------------------------------------------------------*/
        mess_options get_ptr(void) const;


        /*-----------------------------------------------------------------------------
         *  deep copy
         *-----------------------------------------------------------------------------*/
        MESS_options deepcopy(void) const;


        /*-----------------------------------------------------------------------------
         *  print function for octave
         *-----------------------------------------------------------------------------*/
        void print(std::ostream & os, bool pr_as_read_syntax = false) const;
        void print(std::ostream & os, bool pr_as_read_syntax = false);
        void print_raw(std::ostream & os, bool pr_as_read_syntax) const;
        bool print_name_tag(std::ostream & os, const std::string & name) const;


        /*-----------------------------------------------------------------------------
         * overrides for octave
         *-----------------------------------------------------------------------------*/
        bool is_constant(void) const;
        bool is_defined(void) const;
        bool print_as_scalar(void) const;
        bool is_empty(void) const;
        bool is_map (void) const;
        size_t byte_size(void) const;

        string_vector map_keys (void) const;
        octave_value_list subsref(const std::string& type, const std::list<octave_value_list>& idx, int k);
        octave_value subsref (const std::string&, const std::list<octave_value_list>&);
        octave_value subsasgn (const std::string& type, const std::list<octave_value_list>& idx, const octave_value& rhs);


        /*-----------------------------------------------------------------------------
         *  conversion operator
         *-----------------------------------------------------------------------------*/
        operator mess_options() const;
        operator bool() const;


        /*-----------------------------------------------------------------------------
         *  destructor
         *-----------------------------------------------------------------------------*/
        ~MESS_options(void);


        /*-----------------------------------------------------------------------------
         * installation routines for interpreter
         *-----------------------------------------------------------------------------*/
        static bool type_installed;
        static void install(void);


    private:

        mess_options ptr;

        DECLARE_OV_TYPEID_FUNCTIONS_AND_DATA;

};


#endif







