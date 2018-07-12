#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, see <http://www.gnu.org/licenses/>.
#
# Copyright (C) Peter Benner, Martin Koehler, Jens Saak and others
#               2009-2018
#

r"""Interface to the
    mess_eps, mess_version, mess_version_verbose,
    mess_is_debug,
    mess_have_zlib, mess_have_bzip2, mess_have_umfpack,
    mess_have_amd, mess_have_colamd, mess_have_cholmod,
    mess_have_csparse, mess_have_superlu, mess_have_mklpardiso,
    mess_have_arpack, mess_have_matio, mess_have_openmp, mess_have_mess64,
    mess_git_id, mess_git_branch,
    mess_set_errorlevel functions and definitions of Exceptions."""

from pymess._c_interface import _eps, _mess_version, _mess_version_verbose,                     \
                                _mess_version_major, _mess_version_minor, _mess_version_patch,  \
                                _mess_is_debug,                                                 \
                                _mess_have_bzip2, _mess_have_zlib,                              \
                                _mess_have_umfpack, _mess_have_amd, _mess_have_colamd,          \
                                _mess_have_cholmod, _mess_have_csparse,                         \
                                _mess_have_superlu, _mess_have_mklpardiso,                      \
                                _mess_have_arpack, _mess_have_matio, _mess_have_openmp,         \
                                _mess_have_mess64, _mess_git_branch, _mess_git_id,              \
                                _set_errorlevel


def eps():
    r"""This function represents the interface to the mess_eps function.

    Return the machine epsilon.


    Parameters
    ----------

    Returns
    -------
    eps: double, machine epsilon.

    Raises
    ------

    See Also
    --------

    Notes
    -----

    References
    ----------

    Examples
    --------


    .. doctest::

        >>> from pymess import *
        >>> print(eps())   # doctest: +SKIP

    """
    return _eps()




def mess_version():
    r"""This function represents the interface to the mess_version function.

    Print information about C-M.E.S.S..


    Parameters
    ----------

    Returns
    -------

    Raises
    ------

    See Also
    --------

    Notes
    -----

    References
    ----------

    Examples
    --------


    .. doctest::

        >>> from pymess import *
        >>> mess_version()   # doctest: +SKIP

    """
    return _mess_version()

def mess_version_verbose():
    r"""This function represents the interface to the mess_version_verbose function.

    Print information about C-M.E.S.S..


    Parameters
    ----------

    Returns
    -------

    Raises
    ------

    See Also
    --------

    Notes
    -----

    References
    ----------

    Examples
    --------


    .. doctest::

        >>> from pymess import *
        >>> mess_version_verbose()   # doctest: +SKIP

    """
    return _mess_version_verbose()

def mess_version_major():
    r"""This function represents the interface to the mess_version_major function.

    Return major Version of C-M.E.S.S.


    Parameters
    ----------

    Returns
    -------
    ver: int, major version of C-M.E.S.S.


    Raises
    ------

    See Also
    --------

    Notes
    -----

    References
    ----------

    Examples
    --------


    .. doctest::

        >>> from pymess import *
        >>> mess_version_major()   # doctest: +SKIP

    """
    return _mess_version_major()

def mess_version_minor():
    r"""This function represents the interface to the mess_version_minor function.

    Return minor Version of C-M.E.S.S.


    Parameters
    ----------

    Returns
    -------
    ver: int, minor version of C-M.E.S.S.


    Raises
    ------

    See Also
    --------

    Notes
    -----

    References
    ----------

    Examples
    --------


    .. doctest::

        >>> from pymess import *
        >>> mess_version_minor()   # doctest: +SKIP

    """
    return _mess_version_minor()


def mess_version_patch():
    r"""This function represents the interface to the mess_version_patch function.

    Return patch Version of C-M.E.S.S.


    Parameters
    ----------

    Returns
    -------
    ver: int, patch version of C-M.E.S.S.


    Raises
    ------

    See Also
    --------

    Notes
    -----

    References
    ----------

    Examples
    --------

    .. doctest::

        >>> from pymess import *
        >>> mess_version_patch()   # doctest: +SKIP

    """
    return _mess_version_patch()


def mess_git_id():
    r"""This function represents the interface to the mess_git_id function.

    Return the git id which was used to configure C-M.E.S.S. build.


    Parameters
    ----------

    Returns
    -------
    gitid: string, git id


    Raises
    ------

    See Also
    --------

    Notes
    -----

    References
    ----------

    Examples
    --------

    .. doctest::

        >>> from pymess import *
        >>> mess_git_id()   # doctest: +SKIP

    """
    return _mess_git_id()


def mess_git_branch():
    r"""This function represents the interface to the mess_git_branch function.

    Return the git branch which was used to configure C-M.E.S.S. build.


    Parameters
    ----------

    Returns
    -------
    gitbranch: string, git branch


    Raises
    ------

    See Also
    --------

    Notes
    -----

    References
    ----------

    Examples
    --------

    .. doctest::

        >>> from pymess import *
        >>> mess_git_branch()   # doctest: +SKIP

    """
    return _mess_git_branch()


def mess_is_debug():
    r"""This function represents the interface to the mess_is_debug function.

    Return a True if debug modus is active otherwise False.


    Parameters
    ----------

    Returns
    -------

    Raises
    ------

    See Also
    --------

    Notes
    -----

    References
    ----------

    Examples
    --------

    .. doctest::

        >>> from pymess import *
        >>> mess_is_debug()   # doctest: +SKIP

    """
    return _mess_is_debug()


def mess_have_zlib():
    r"""This function represents the interface to the mess_have_zlib function.

    Return a True if zlib is available otherwise False.


    Parameters
    ----------

    Returns
    -------

    Raises
    ------

    See Also
    --------

    Notes
    -----

    References
    ----------

    Examples
    --------


    .. doctest::

        >>> from pymess import *
        >>> mess_have_zlib()   # doctest: +SKIP

    """
    return _mess_have_zlib()


def mess_have_bzip2():
    r"""This function represents the interface to the mess_have_bzip2 function.

    Return a True if bzip2 is available otherwise False.


    Parameters
    ----------

    Returns
    -------

    Raises
    ------

    See Also
    --------

    Notes
    -----

    References
    ----------

    Examples
    --------


    .. doctest::

        >>> from pymess import *
        >>> mess_have_bzip2()   # doctest: +SKIP

    """
    return _mess_have_bzip2()



def mess_have_umfpack():
    r"""This function represents the interface to the mess_have_umfpack function.

    Return a True if umfpack is available otherwise False.


    Parameters
    ----------

    Returns
    -------

    Raises
    ------

    See Also
    --------

    Notes
    -----

    References
    ----------

    Examples
    --------


    .. doctest::

        >>> from pymess import *
        >>> mess_have_umfpack()   # doctest: +SKIP

    """
    return _mess_have_umfpack()


def mess_have_amd():
    r"""This function represents the interface to the mess_have_amd function.

    Return a True if amd is available otherwise False.


    Parameters
    ----------

    Returns
    -------

    Raises
    ------

    See Also
    --------

    Notes
    -----

    References
    ----------

    Examples
    --------


    .. doctest::

        >>> from pymess import *
        >>> mess_have_amd()   # doctest: +SKIP

    """
    return _mess_have_amd()


def mess_have_colamd():
    r"""This function represents the interface to the mess_have_colamd function.

    Return a True if colamd is available otherwise False.


    Parameters
    ----------

    Returns
    -------

    Raises
    ------

    See Also
    --------

    Notes
    -----

    References
    ----------

    Examples
    --------


    .. doctest::

        >>> from pymess import *
        >>> mess_have_colamd()   # doctest: +SKIP

    """
    return _mess_have_colamd()



def mess_have_cholmod():
    r"""This function represents the interface to the mess_have_cholmod function.

    Return a True if cholmod is available otherwise False.


    Parameters
    ----------

    Returns
    -------

    Raises
    ------

    See Also
    --------

    Notes
    -----

    References
    ----------

    Examples
    --------


    .. doctest::

        >>> from pymess import *
        >>> mess_have_cholmod()   # doctest: +SKIP

    """
    return _mess_have_cholmod()



def mess_have_csparse():
    r"""This function represents the interface to the mess_have_csparse function.

    Return a True if csparse is available otherwise False.


    Parameters
    ----------

    Returns
    -------

    Raises
    ------

    See Also
    --------

    Notes
    -----

    References
    ----------

    Examples
    --------


    .. doctest::

        >>> from pymess import *
        >>> mess_have_csparse()   # doctest: +SKIP

    """
    return _mess_have_csparse()


def mess_have_superlu():
    r"""This function represents the interface to the mess_have_superlu function.

    Return a True if superlu is available otherwise False.


    Parameters
    ----------

    Returns
    -------

    Raises
    ------

    See Also
    --------

    Notes
    -----

    References
    ----------

    Examples
    --------


    .. doctest::

        >>> from pymess import *
        >>> mess_have_superlu()   # doctest: +SKIP

    """
    return _mess_have_superlu()


def mess_have_mklpardiso():
    r"""This function represents the interface to the mess_have_mklpardiso function.

    Return a True if MKL Pardiso is available otherwise False.


    Parameters
    ----------

    Returns
    -------

    Raises
    ------

    See Also
    --------

    Notes
    -----

    References
    ----------

    Examples
    --------


    .. doctest::

        >>> from pymess import *
        >>> mess_have_mklpardiso()   # doctest: +SKIP

    """
    return _mess_have_mklpardiso()


def mess_have_arpack():
    r"""This function represents the interface to the mess_have_arpack function.

    Return a True if arpack is available otherwise False.


    Parameters
    ----------

    Returns
    -------

    Raises
    ------

    See Also
    --------

    Notes
    -----

    References
    ----------

    Examples
    --------


    .. doctest::

        >>> from pymess import *
        >>> mess_have_arpack()   # doctest: +SKIP

    """
    return _mess_have_arpack()


def mess_have_matio():
    r"""This function represents the interface to the mess_have_matio function.

    Return a True if matio is available otherwise False.


    Parameters
    ----------

    Returns
    -------

    Raises
    ------

    See Also
    --------

    Notes
    -----

    References
    ----------

    Examples
    --------


    .. doctest::

        >>> from pymess import *
        >>> mess_have_matio()   # doctest: +SKIP

    """
    return _mess_have_matio()


def mess_have_openmp():
    r"""This function represents the interface to the mess_have_openmp function.

    Return a True if openmp is available otherwise False.


    Parameters
    ----------

    Returns
    -------

    Raises
    ------

    See Also
    --------

    Notes
    -----

    References
    ----------

    Examples
    --------


    .. doctest::

        >>> from pymess import *
        >>> mess_have_openmp()   # doctest: +SKIP

    """
    return _mess_have_openmp()


def mess_have_mess64():
    r"""This function represents the interface to the mess_have_mess64 function.

    Return a True if 64 bit integers are available otherwise False.


    Parameters
    ----------

    Returns
    -------

    Raises
    ------

    See Also
    --------

    Notes
    -----

    References
    ----------

    Examples
    --------


    .. doctest::

        >>> from pymess import *
        >>> mess_have_mess64()   # doctest: +SKIP

    """
    return _mess_have_mess64()




def set_errorlevel(lvl):
    r"""This function represents the interface to the mess_set_errorlevel function.
    Set the verbosity of the error output.

    Returns None.


    Parameters
    ----------
    lvl, int

        * 0, nothing is displayed
        * 1, only errors are displayed
        * 2, errors and warnings are displayed
        * 3, errors, warnings and debug messages are displayed

    Returns
    -------

    Raises
    -------
    MESSErrorArguments

    See Also
    --------

    Notes
    -----

    References
    ----------

    Examples
    --------


    .. doctest::

        >>> from pymess import *
        >>> set_errorlevel(0)   # doctest: +SKIP
        >>> set_errorlevel(1)   # doctest: +SKIP
        >>> set_errorlevel(2)   # doctest: +SKIP
        >>> set_errorlevel(3)   # doctest: +SKIP

    """
    if (not isinstance(lvl, int)) or lvl < 0 or lvl > 4:
        raise MESSErrorArguments("lvl is not a valid int.")
    return _set_errorlevel(lvl)



class MESSErrorArguments(Exception):
    r"""
    Exception raised for errors in the input data.
    """
    pass


class MESSErrorData(Exception):
    r"""
    Exception raised for errors in the input data.
    """
    pass


class MESSErrorNotSupported(Exception):
    r"""
    Exception raised if function is not supported.
    """
    pass


class MESSErrorEquationType(Exception):
    r"""
    Exception raised if wrong equation type is passed.
    """
    pass
