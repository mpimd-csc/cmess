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

/**
 * @file tutorials/plot/tutorial_plot.c
 * @brief Tutorial for the plotting interface and plot the sine and cosine function.
 * @author @mbehr
 * @author @koehlerm
 *
 * # Tutorial Plotting Interface
 * In this tutorial we use the plotting interface to plot the sine and cosine function.\n
 * More precisely we plot  \f$ 20\sin(x),\ 20\cos(x),\ 20\sin(x)+20\cos(x)\f$.
 * We export the graphics to a @gnuplot script as well as a @tikz script.
 *
 * ## 1. Include Header Files
 * We include standard @c C header files and in addtion the @c mess.h header file
 * to make the @mess library available.
 * @snippet "tutorials/plot/tutorial_plot.c" HEADER
 * \n
 *
 * ## 2. Declare mess_plotExport and additional variables.
 * We use the \ref mess_plotExport type to generate a plot interface object
 * and some additional variables to generate the data.
 * @snippet "tutorials/plot/tutorial_plot.c" DECLARE
 * \n
 *
 * ## 3. Init mess_plotExport Instance
 * Before we are able to use the @ref mess_plotExport instance we have to init the fields of the structure.
 * Therefore we use @ref mess_plotExport_init. We name our plot @c TestPlot choose @ref MESS_PLOT_LIN for linear scaling the @c x-axis
 * as well as the @c y-axis. For logarithmic scaling use @ref MESS_PLOT_LOG. The label of the @c x-axis
 * is @c x and the lable fo the @c y-axis is @c y. \n
 * The last argument indicates that we want to store @c 3 datasets.
 * @snippet "tutorials/plot/tutorial_plot.c" INIT
 * \n
 *
 * ## 4. Adjust the Plot
 * Now we are ready to use the @c plot instance of @ref mess_plotExport.
 * We use the functions @ref mess_plotExport_setLegendPos, @ref mess_plotExport_setLabel,
 * @ref mess_plotExport_setType to specify the position of the legend and for each
 * dataset the label, type and color of the plot. The second argument is used to reference to the dataset. \n
 * You can find more options by reading the corresponding documentation of the functions.
 * @snippet "tutorials/plot/tutorial_plot.c" ADJUST
 * \n
 *
 * ## 5. Fill data
 * To fill the @c plot instance with data we use a simple @c for loop compute the
 * function values of the functions and use @ref mess_plotExport_addData to add the
 * arguments as well as the function values to the corresponding data set.
 * @snippet "tutorials/plot/tutorial_plot.c" DATA
 * \n
 *
 * ## 6. Export
 * Exporting the plot is very easy using @ref mess_plotExport_createGnuScript and
 * @ref mess_plotExport_createTikzScript to export the data.
 * The function create a @gnuplot script and a datafile and a
 * @tikz file with corresponding data file. \n
 *
 * The name of the generated files
 * starts with @c TestPlot*.
 * You can also use
 * @ref  mess_plotExport_createTikzScript_filename or
 * @ref mess_plotExport_createGnuScript_filename to specify the name of the output.
 *
 * @snippet "tutorials/plot/tutorial_plot.c" EXPORT
 * \n
 *
 * ## 7. Clear mess_plotExport instance.
 * Do not forget to release the memory allocated by @ref mess_plotExport_init by calling
 * @ref mess_plotExport_clear.
 * @snippet "tutorials/plot/tutorial_plot.c" CLEAR
 * @sa @ref tutorials
 */

//cond -> code is ignored by doxygen documentation, but snippet references are working

///@cond
///[HEADER]
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <math.h>
#include "mess/mess.h"
///[HEADER]

int main ( int argc, char **argv){

    /// [DECLARE]
    mess_int_t i, N=100;
    double x, y1, y2, h = 2.0*3.14159/N;
    mess_plotExport plot;
    /// [DECLARE]

    /// [INIT]
    mess_plotExport_init(&plot, "TestPlot", MESS_PLOT_LIN, MESS_PLOT_LIN, "x", "y", 3);
    /// [INIT]

    /// [ADJUST]
    mess_plotExport_setLegendPos(plot, "belowCenter");
    mess_plotExport_setLabel(plot,  0, "20\\sin(x)");
    mess_plotExport_setColor(plot,  0, "blue");
    mess_plotExport_setType(plot,   0, "lines");
    mess_plotExport_setLabel(plot,  1, "20\\cos(x)");
    mess_plotExport_setColor(plot,  1, "yellow");
    mess_plotExport_setType(plot,   1, "lines");
    mess_plotExport_setLabel(plot,  2, "20\\sin(x)+20\\cos(x)");
    mess_plotExport_setColor(plot,  2, "green");
    mess_plotExport_setType(plot,   2, "lines");
    printf("\nNew Plot initialised: \"TestPlot\" with 3 datasets. Keybox orientation: belowCenter\n");
    /// [ADJUST]

    /// [DATA]
    for ( i = 0; i <= N; i++){
        x = h*i;
        y1 = sin(x)*20;
        y2 = cos(x)*20;
        mess_plotExport_addData(plot,0, x,y1 );
        mess_plotExport_addData(plot,1, x,y2 );
        mess_plotExport_addData(plot,2, x,y1+y2 );
    }
    /// [DATA]

    /// [EXPORT]
    printf("\ngnu_plot_script runs\n");
    mess_plotExport_createGnuScript(plot);

    printf("Tikz_plot_script runs\n");
    mess_plotExport_createTikzScript(plot);
    /// [EXPORT]

    /// [CLEAR]
    mess_plotExport_clear(&plot);
    /// [CLEAR]

    return 0;
}
///@endcond



