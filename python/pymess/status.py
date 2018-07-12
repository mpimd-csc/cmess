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

r"""This module contains :class:`Status` class, which represents the status structure mess_status."""


class Status(object):
    r"""This class represents the status structure mess_status.

    Parameters
    ----------

    res2_norms: vector
        Absolute 2-norm in every step.

    rel_changes: vector
        Relative change of the Iterate in every step in
        Frobenius norm for :func:`~pymess.lradi.lradi`.
        Relative change of the 2-Norm residual for :func:`~pymess.lrnm.lrnm`.

    res2_norm: real
        Final absolute 2-norm residual.

    it: integer
        Number of iteration the algorithm took.

    res2_change: real
        Final relative 2-norm residual.

    res2_0: real
        2-norm of the right hand side.

    rel_change: real
        Final relative change of :math:`Z`.

    stop_res2: integer
        Number of iteration the algorithm stops with respect to relative 2-norm residual.

    stop_res2c: integer
        Number of iteration the algorithm stops with respect to  2-norm relative change.

    stop_rel: integer
        Number of iteration the algorithm stops with respect to relative change in :math:`Z`.

    stop_user: integer
        Flag if iteration was cancelled.

    time_all: real
        Overall time.

    time_adi: real
        Time for ADI.

    unstable: integer
        True if instable ritz values occur.

    n_internal_status: integer
        Number of internal :class:`.Status` instances in internal_status.

    internal_status: list of status
        List of :class:`.Status` instances.

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

    """

    _name = "pymess.status"

    def __init__(self):
        self.res2_norms = None
        self.rel_changes = None
        self.res2_norm = None
        self.it = None
        self.res2_change = None
        self.res2_0 = None
        self.rel_change = None
        self.stop_res2 = None
        self.stop_res2c = None
        self.stop_rel = None
        self.stop_user = None
        self.time_all = None
        self.time_adi = None
        self.unstable = None
        self.n_internal_status = 0
        self.internal_status = None

    def __repr__(self):
        ret = ""

        if self.n_internal_status == 0:
            # lradi
            ret += "iterations:    \t {0:d}\n".format(self.it)
            ret += "res2_norm:     \t {0:8e}\n".format(self.res2_norm)
            ret += "res2_change:   \t {0:8e}\n".format(self.res2_change)
            ret += "rel_change:    \t {0:8e}\n".format(self.rel_change)
            ret += "res2_0:        \t {0:8e}\n".format(self.res2_0)
            ret += "stop_res2:     \t {0:d}\n".format(self.stop_res2)
            ret += "stop_res2c:    \t {0:d}\n".format(self.stop_res2c)
            ret += "stop_rel:      \t {0:d}\n".format(self.stop_rel)
            ret += "stop_user:     \t {0:d}\n".format(self.stop_user)
            ret += "time - all:    \t {0:8e}\n".format(self.time_all)
            ret += "time - adi:    \t {0:8e}\n".format(self.time_adi)
            ret += "unstable:      \t {0:d}\n\n".format(self.unstable)
            ret += " it   |  ||V||_F/||Z||_F  |  ||L(ZZ')||_2/||RHS||_2  |  ||L(ZZ')||_2  \n"
            ret += "------+-------------------+--------------------------+----------------\n"
            for j in range(self.it):
                relc = self.rel_changes[j]
                res2 = self.res2_norms[j]
                res2_0 = self.res2_0
                ret += "{0:5d} |    {1:9e}   |       {2:9e}       |   {3:9e}\n".format(
                    j, relc, res2 / res2_0, res2)
        else:
            # lrnm
            for stat in [s for s in self.internal_status if s is not None]:
                ret += "iterations:    \t {0:d}\n".format(stat.it)
                ret += "res2_norm:     \t {0:8e}\n".format(stat.res2_norm)
                ret += "res2_change:   \t {0:8e}\n".format(stat.res2_change)
                ret += "rel_change:    \t {0:8e}\n".format(stat.rel_change)
                ret += "res2_0:        \t {0:8e}\n".format(stat.res2_0)
                ret += "stop_res2:     \t {0:d}\n".format(stat.stop_res2)
                ret += "stop_res2c:    \t {0:d}\n".format(stat.stop_res2c)
                ret += "stop_rel:      \t {0:d}\n".format(stat.stop_rel)
                ret += "stop_user:     \t {0:d}\n".format(stat.stop_user)
                ret += "time - all:    \t {0:8e}\n".format(stat.time_all)
                ret += "time - adi:    \t {0:8e}\n".format(stat.time_adi)
                ret += "unstable:      \t {0:d}\n\n".format(stat.unstable)
                ret += " it   |  ||V||_F/||Z||_F  |  ||L(ZZ')||_2/||RHS||_2  |  ||L(ZZ')||_2  \n"
                ret += "------+-------------------+--------------------------+----------------\n"
                for j in range(stat.it):
                    relc = stat.rel_changes[j]
                    res2 = stat.res2_norms[j]
                    res2_0 = stat.res2_0
                    ret += "{0:5d} |    {1:9e}   |       {2:9e}       |   {3:9e}\n".format(
                        j, relc, res2 / res2_0, res2)
                ret += "\n"

            ret += "iterations:    \t {0:d}\n".format(self.it)
            ret += "res2_norm:     \t {0:8e}\n".format(self.res2_norm)
            ret += "res2_change:   \t {0:8e}\n".format(self.res2_change)
            ret += "rel_change:    \t {0:8e}\n".format(self.rel_change)
            ret += "res2_0:        \t {0:8e}\n".format(self.res2_0)
            ret += "stop_res2:     \t {0:d}\n".format(self.stop_res2)
            ret += "stop_res2c:    \t {0:d}\n".format(self.stop_res2c)
            ret += "stop_rel:      \t {0:d}\n".format(self.stop_rel)
            ret += "stop_user:     \t {0:d}\n".format(self.stop_user)
            ret += "time - all:    \t {0:8e}\n".format(self.time_all)
            ret += "time - adi:    \t {0:8e}\n".format(self.time_adi)
            ret += "unstable:      \t {0:d}\n\n".format(self.unstable)
            ret += " it   |   rel_change   |  ||R(X)||_2/||RHS||_2  |     ||R(X)||_2     | Inner ADI Steps\n"
            ret += "------+----------------+------------------------+--------------------+----------------\n"
            for j in range(self.it):
                relc = self.rel_changes[j]
                res2 = self.res2_norms[j]
                res2_0 = self.res2_0
                it = self.internal_status[j].it
                ret += "{0:5d} |  {1:8e}  |      {2:8e}      |    {3:8e}    |      {4:d}\n".format(
                    j, relc, res2 / res2_0, res2, it)

        #return ret.encode("utf-8")
        return ret

    def lradi_stat(self):
        """Return a string which represents the status instance, after :func:`lradi` was called."""
        ret = ""
        # lradi
        ret += "iterations:    \t {0:d}\n".format(self.it)
        ret += "res2_norm:     \t {0:8e}\n".format(self.res2_norm)
        ret += "res2_change:   \t {0:8e}\n".format(self.res2_change)
        ret += "rel_change:    \t {0:8e}\n".format(self.rel_change)
        ret += "res2_0:        \t {0:8e}\n".format(self.res2_0)
        ret += "stop_res2:     \t {0:d}\n".format(self.stop_res2)
        ret += "stop_res2c:    \t {0:d}\n".format(self.stop_res2c)
        ret += "stop_rel:      \t {0:d}\n".format(self.stop_rel)
        ret += "stop_user:     \t {0:d}\n".format(self.stop_user)
        ret += "time - all:    \t {0:8e}\n".format(self.time_all)
        ret += "time - adi:    \t {0:8e}\n".format(self.time_adi)
        ret += "unstable:      \t {0:d}\n\n".format(self.unstable)
        ret += " it   |  ||V||_F/||Z||_F  |  ||L(ZZ')||_2/||RHS||_2  |   ||L(ZZ')||_2  \n"
        ret += "------+-------------------+--------------------------+------------------\n"
        for j in range(self.it):
            relc = self.rel_changes[j]
            res2 = self.res2_norms[j]
            res2_0 = self.res2_0
            ret += "{0:5d} |    {1:9e}   |       {2:9e}       |   {3:9e}\n".format(
                j, relc, res2 / res2_0, res2)

        #return ret.encode("utf-8")
        return ret

    def lrnm_stat(self):
        """Return a string which represents the status instance, after :func:`lrnm` was called."""
        ret = ""
        ret += "iterations:    \t {0:d}\n".format(self.it)
        ret += "res2_norm:     \t {0:8e}\n".format(self.res2_norm)
        ret += "res2_change:   \t {0:8e}\n".format(self.res2_change)
        ret += "rel_change:    \t {0:8e}\n".format(self.rel_change)
        ret += "res2_0:        \t {0:8e}\n".format(self.res2_0)
        ret += "stop_res2:     \t {0:d}\n".format(self.stop_res2)
        ret += "stop_res2:     \t {0:d}\n".format(self.stop_res2)
        ret += "stop_res2c:    \t {0:d}\n".format(self.stop_res2c)
        ret += "stop_rel:      \t {0:d}\n".format(self.stop_rel)
        ret += "stop_user:     \t {0:d}\n".format(self.stop_user)
        ret += "time - all:    \t {0:8e}\n".format(self.time_all)
        ret += "time - adi:    \t {0:8e}\n".format(self.time_adi)
        ret += "unstable:      \t {0:d}\n\n".format(self.unstable)
        ret += " it   |  ||R(X)||_2/||RHS||_2  |     ||R(X)||_2     | (rel.) Res-2 Change   | Inner ADI Steps \n"
        ret += "------+------------------------+--------------------+-----------------------+-----------------\n"
        for j in range(self.it):
            relc = self.rel_changes[j]
            res2 = self.res2_norms[j]
            res2_0 = self.res2_0
            it = self.internal_status[j].it
            ret += "{0:5d} |      {1:8e}      |    {2:8e}    |     {3:8e}      |      {4:d}\n".format(
                j, res2 / res2_0, res2, relc, it)

        #return ret.encode("utf-8")
        return ret

    def __delattr__(self, key):
        raise TypeError(
            "Trying to delete attributes of %r is not allowed" % self)
