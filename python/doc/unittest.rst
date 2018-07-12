..
   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.
   You should have received a copy of the GNU General Public License
   along with this program; if not, see <http://www.gnu.org/licenses/>.
   Copyright (C) Peter Benner, Martin Koehler, Jens Saak and others
                 2009-2018

.. include:: global.rst

.. _unittest:



Unit Testing Framework
======================

After successfull :ref:`Installation <installation>`
you should change directory to :code:`<MY BUILD>/python`.
You can run all unittests using

.. code-block:: bash

    python unittests/run.py


The unittest framework of |pymess| is divided into several units.

You can run specific units by calling :code:`run.py` with an input argument.


======================================================== ===================================================
Description of Tests                                      Command
======================================================== ===================================================
Datatype Conversions                                      :code:`python unittests/run.py CONVERSION`
:func:`~pymess.easy.lyap` and :func:`~pymess.easy.care`   :code:`python unittests/run.py EASY`
:func:`~pymess.lradi.lradi`                               :code:`python unittests/run.py LRADI`
:func:`~pymess.lrnm.lrnm`                                 :code:`python unittests/run.py LRNM`
Callback Functionality                                    :code:`python unittests/run.py CALLBACK`
======================================================== ===================================================

