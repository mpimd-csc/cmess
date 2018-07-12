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

#include <octave/oct.h>
#include <octave/ops.h>
#include <octave/ov-complex.h>
#include <octave/ov-re-sparse.h>
#include <math.h>
#include "mess/mess.h"
#include "MESS_options.h"
#include "MESS_error.h"


DEFINE_OV_TYPEID_FUNCTIONS_AND_DATA (MESS_options, "mess_options", "mess_options");


static const char * const  fieldnames [] =
{
    /* type of equation */
    "type",

    /* shift parameter related fields */
    "adi_shifts_p", "adi_shifts_paratype",
    "adi_shifts_arp_p", "adi_shifts_arp_m", "adi_shifts_l0", "adi_shifts_b0",
    "proj_space",

    /* ADI related fields */
    "adi_maxit", "adi_res2_tol", "adi_res2c_tol", "adi_rel_change_tol",
    "adi_output", "residual_method",

    /* NM related fields */
    "nm_maxit", "nm_res2_tol", "nm_output", "nm_linesearch",
    "nm_singleshifts", "nm_gpStep",

    /* others and internal */
    "memory_usage", "K0", "W"
};

static const size_t len_fieldnames = 23;


/*-----------------------------------------------------------------------------
 *  constructors
 *-----------------------------------------------------------------------------*/
MESS_options::MESS_options(void):ptr(NULL){ }
MESS_options::MESS_options(mess_options opt):ptr(opt) { }
MESS_options::MESS_options(const MESS_options & opt):ptr(opt.ptr){ };


/*-----------------------------------------------------------------------------
 *  getter
 *-----------------------------------------------------------------------------*/
mess_options MESS_options::get_ptr(void) const{ return this->ptr; }


/*-----------------------------------------------------------------------------
 *  deep copy
 *-----------------------------------------------------------------------------*/
MESS_options MESS_options::deepcopy(void) const{
    MESS_options opt;
    if(*this){
        mess_options ptr = opt.get_ptr();
        mess_options_init(&ptr);
        mess_options_copy(this->ptr,ptr);
    }
    return opt;
}


/*-----------------------------------------------------------------------------
 *  print functions octave
 *-----------------------------------------------------------------------------*/
void MESS_options::print(std::ostream & os, bool pr_as_read_syntax) const {
    print_raw(os, pr_as_read_syntax);
    newline(os);
}


void MESS_options::print(std::ostream & os, bool pr_as_read_syntax) {
    print_raw(os, pr_as_read_syntax);
    newline(os);
}


void MESS_options::print_raw(std::ostream & os, bool pr_as_read_syntax) const {
    os << MESS_options::static_type_name();
#ifdef MESS_DEBUG
    os << " ( " << this->ptr << " ) ";
#endif
}


bool MESS_options::print_name_tag (std::ostream& os, const std::string& name) const {
    os << name << " = ";
    return true;
}


/*-----------------------------------------------------------------------------
 * overrides for octave
 *-----------------------------------------------------------------------------*/
bool MESS_options::is_constant(void) const { return true; }
bool MESS_options::is_defined(void) const { return true; }
bool MESS_options::print_as_scalar(void) const { return true; }
bool MESS_options::is_empty(void) const { return this->ptr?false:true; }
bool MESS_options::is_map (void) const { return true; }
size_t MESS_options::byte_size(void) const { return this->ptr? mess_options_memsize(*this):0; }
string_vector MESS_options::map_keys (void) const { return string_vector(fieldnames,len_fieldnames);}

octave_value MESS_options::subsref (const std::string& type, const std::list<octave_value_list>& idx){
    return subsref(type,idx,0);
}


octave_value_list MESS_options::subsref (const std::string& type, const std::list<octave_value_list>& idx, int k){


    /*-----------------------------------------------------------------------------
     * check index operation
     *-----------------------------------------------------------------------------*/
    if (idx.size()!=1){
        error ("invalid index operation for %s", MESS_options::static_type_name().c_str());
        return octave_value();
    }

    octave_value_list idxlist = idx.front();

    if(idxlist.length()!=1){
        error ("invalid index operation for %s", MESS_options::static_type_name().c_str());
        return octave_value();
    }

    std::string indexname = idxlist(0).string_value();

    if(!(*this)){
        error("%s instance points to NULL.",MESS_options::static_type_name().c_str());
        return octave_value();
    }

    /*-----------------------------------------------------------------------------
     *  perform index operation
     *-----------------------------------------------------------------------------*/
    switch (type[0]){
        case '(':
        case '{':
            error ("%s cannot be indexed with %c", MESS_options::static_type_name().c_str(), type[0]);
            break;
        case '.':

            /*-----------------------------------------------------------------------------
             * type of equation
             *-----------------------------------------------------------------------------*/
            if(indexname == "type"){
                return octave_value(this->ptr->type);

                /*-----------------------------------------------------------------------------
                 * Shift Parameter related fields
                 *-----------------------------------------------------------------------------*/
            }else if(indexname == "adi_shifts_p"){
                error("not implemented");
                return octave_value();
            }else if(indexname == "adi_shifts_paratype"){
                return octave_value(this->ptr->adi_shifts_paratype);
            }else if(indexname == "adi_shifts_arp_p"){
                return octave_value(this->ptr->adi_shifts_arp_p);
            }else if(indexname == "adi_shifts_arp_m"){
                return octave_value(this->ptr->adi_shifts_arp_m);
            }else if(indexname == "adi_shifts_l0"){
                return octave_value(this->ptr->adi_shifts_l0);
            }else if(indexname == "adi_shifts_b0"){
                error("not implemented");
                return octave_value();
            }else if(indexname == "proj_space"){
                error("not implemented");
                return octave_value();

                /*-----------------------------------------------------------------------------
                 * ADI related fields
                 *-----------------------------------------------------------------------------*/
            }else if(indexname == "adi_maxit"){
                return octave_value(this->ptr->adi_maxit);
            }else if(indexname == "adi_res2_tol"){
                return octave_value(this->ptr->adi_res2_tol);
            }else if(indexname == "adi_res2c_tol"){
                return octave_value(this->ptr->adi_res2c_tol);
            }else if(indexname == "adi_rel_change_tol"){
                return octave_value(this->ptr->adi_rel_change_tol);
            }else if(indexname == "adi_output"){
                return octave_value(this->ptr->adi_output);
            }else if(indexname == "residual_method"){
                return octave_value(this->ptr->residual_method);

                /*-----------------------------------------------------------------------------
                 * NM related fields
                 *-----------------------------------------------------------------------------*/
            }else if(indexname == "nm_maxit"){
                return octave_value(this->ptr->nm_maxit);
            }else if(indexname == "nm_res2_tol"){
                return octave_value(this->ptr->nm_res2_tol);
            }else if(indexname == "nm_output"){
                return octave_value(this->ptr->nm_output);
            }else if(indexname == "nm_linesearch"){
                return octave_value(this->ptr->nm_linesearch);
            }else if(indexname == "nm_singleshifts"){
                return octave_value(this->ptr->nm_singleshifts);
            }else if(indexname == "nm_gpStep"){
                return octave_value(this->ptr->nm_gpStep);

                /*-----------------------------------------------------------------------------
                 *  others and internal
                 *-----------------------------------------------------------------------------*/
            }else if(indexname == "memory_usage"){
                return octave_value(this->ptr->memory_usage);
            }else if(indexname == "K0"){
                return octave_value(this->ptr->K0);
            }else if(indexname == "W"){
                return octave_value(this->ptr->W);

                /*-----------------------------------------------------------------------------
                 * not implemented
                 *-----------------------------------------------------------------------------*/
            }else{
                error ("%s cannot be indexed with %s", MESS_options::static_type_name().c_str(), indexname.c_str());
                return octave_value();
            }
            break;
        default:
            panic_impossible ();
    }

    return octave_value();

}


octave_value MESS_options::subsasgn (const std::string& type, const std::list<octave_value_list>& idx, const octave_value& rhs){

#if 0
    /*-----------------------------------------------------------------------------
     * check index operation
     *-----------------------------------------------------------------------------*/
    if (idx.size()!=1){
        error ("invalid index operation for %s", MESS_options::static_type_name().c_str());
        return octave_value();
    }

    octave_value_list idxlist = idx.front();

    if(idxlist.length()!=1){
        error ("invalid index operation for %s", MESS_options::static_type_name().c_str());
        return octave_value();
    }

    std::string indexname = idxlist(0).string_value();



    /*-----------------------------------------------------------------------------
     *  perform index operation
     *-----------------------------------------------------------------------------*/
    switch (type[0]){
        case '(':
        case '{':
            error ("%s cannot be indexed with %c", MESS_options::static_type_name().c_str(), type[0]);
            break;
        case '.':

            /*-----------------------------------------------------------------------------
             * type of equation
             *-----------------------------------------------------------------------------*/
            if(indexname == "type"){
                return octave_value(this->ptr->type);

                /*-----------------------------------------------------------------------------
                 * Shift Parameter related fields
                 *-----------------------------------------------------------------------------*/
            }else if(indexname == "adi_shifts_p"){
                error("not implemented");
                return octave_value();
            }else if(indexname == "adi_shifts_paratype"){
                return octave_value(this->ptr->adi_shifts_paratype);
            }else if(indexname == "adi_shifts_arp_p"){
                return octave_value(this->ptr->adi_shifts_arp_p);
            }else if(indexname == "adi_shifts_arp_m"){
                return octave_value(this->ptr->adi_shifts_arp_m);
            }else if(indexname == "adi_shifts_l0"){
                return octave_value(this->ptr->adi_shifts_l0);
            }else if(indexname == "adi_shifts_b0"){
                error("not implemented");
                return octave_value();
            }else if(indexname == "proj_space"){
                error("not implemented");
                return octave_value();

                /*-----------------------------------------------------------------------------
                 * ADI related fields
                 *-----------------------------------------------------------------------------*/
            }else if(indexname == "adi_maxit"){
                return octave_value(this->ptr->adi_maxit);
            }else if(indexname == "adi_res2_tol"){
                return octave_value(this->ptr->adi_res2_tol);
            }else if(indexname == "adi_res2c_tol"){
                return octave_value(this->ptr->adi_res2c_tol);
            }else if(indexname == "adi_rel_change_tol"){
                return octave_value(this->ptr->adi_rel_change_tol);
            }else if(indexname == "adi_output"){
                return octave_value(this->ptr->adi_output);
            }else if(indexname == "residual_method"){
                return octave_value(this->ptr->residual_method);

                /*-----------------------------------------------------------------------------
                 * NM related fields
                 *-----------------------------------------------------------------------------*/
            }else if(indexname == "nm_maxit"){
                return octave_value(this->ptr->nm_maxit);
            }else if(indexname == "nm_res2_tol"){
                return octave_value(this->ptr->nm_res2_tol);
            }else if(indexname == "nm_output"){
                return octave_value(this->ptr->nm_output);
            }else if(indexname == "nm_linesearch"){
                return octave_value(this->ptr->nm_linesearch);
            }else if(indexname == "nm_singleshifts"){
                return octave_value(this->ptr->nm_singleshifts);
            }else if(indexname == "nm_gpStep"){
                return octave_value(this->ptr->nm_gpStep);

                /*-----------------------------------------------------------------------------
                 *  others and internal
                 *-----------------------------------------------------------------------------*/
            }else if(indexname == "memory_usage"){
                return octave_value(this->ptr->memory_usage);
            }else if(indexname == "K0"){
                return octave_value(this->ptr->K0);
            }else if(indexname == "W"){
                return octave_value(this->ptr->W);

                /*-----------------------------------------------------------------------------
                 * not implemented
                 *-----------------------------------------------------------------------------*/
            }else{
                error ("%s cannot be indexed with %s", MESS_options::static_type_name().c_str(), indexname.c_str());
                return octave_value();
            }
            break;
        default:
            panic_impossible ();
    }
#endif
    return octave_value();


}


/*-----------------------------------------------------------------------------
 *  conversion operator
 *-----------------------------------------------------------------------------*/
MESS_options::operator mess_options(void) const { return this->ptr; }
MESS_options::operator bool(void) const { return this->ptr != NULL; }


/*-----------------------------------------------------------------------------
 *  destructor
 *-----------------------------------------------------------------------------*/
MESS_options::~MESS_options(void){
    if(this->ptr){
        mess_options_clear(&(this->ptr));
    }
    this->ptr = NULL;
}




