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
 * @file lib/plot/plot_scriptExporter.c
 * @brief A small plot script creator.
 * @author @weiss
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <complex.h>
#include <math.h>
#include <pthread.h>
#include <unistd.h>
#include <string.h>
#include "mess/config.h"
#include "mess/mess.h"
#include "mess/error_macro.h"

/**
 * @brief Initialize all data structures and variables.
 * @param[in,out] p plot
 * @param[in] title input plot title
 * @param[in] xscale input scale of the x axis (not used)
 * @param[in] yscale  input scale of the y axis (not used)
 * @param[in] xLabel  input label of the x axis
 * @param[in] yLabel  input label of the y axis
 * @param[in] data  input number of data elements
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref mess_plotExport_init function initializes all variables of a plot.
 */
int mess_plotExport_init(mess_plotExport *p, char *title, mess_plot_axis_scale_t xscale, mess_plot_axis_scale_t  yscale, char *xLabel, char *yLabel, int data){
    MSG_FNAME(__func__);
    int i = 0;

    mess_check_positive(data);

    mess_try_alloc(*p , mess_plotExport, sizeof (struct __mess_plotterScript) );
    (*p)->plot_len = 0;
    mess_try_calloc((*p)->title, char * , (strlen(title)+1) , sizeof ( char ));

    (*p)->plot_len = data;
    mess_try_alloc((*p)->plot, mess_plotExportData *, data*sizeof(mess_plotExportData));
    for (i = 0; i < data; i++) {
        (*p)->plot[i].x = NULL;
        (*p)->plot[i].y = NULL;
        (*p)->plot[i].len = 0;
        strncpy((*p)->plot[i].color , "",40);
        strncpy((*p)->plot[i].type, "", 20);
        strncpy((*p)->plot[i].label, "", 40);
    }
    strncpy((*p)->title, title, strlen(title));
    (*p)->legendGnuPos= "top right";
    (*p)->legendTikzPos= ", legend pos=north east";
    (*p)->xLabel = xLabel;
    (*p)->yLabel = yLabel;
    (*p)->xscale = xscale;
    (*p)->yscale = yscale;

    return 0;
}

/**
 * @brief Clear all data in a plot.
 * @param[in,out] pl plot with data
 * @return always zero
 *
 * The @ref mess_plotExport_clear function clears all information in the plot \f$ pl \f$.
 */
int mess_plotExport_clear(mess_plotExport *pl) {
    if ( pl == NULL ) return 0;
    if ( *pl == NULL ) return 0;
    int data = (*pl)->plot_len;
    int i ;
    for (i = 0; i < data; i++) {
        if ( (*pl)->plot[i].x !=  NULL)mess_free((*pl)->plot[i].x );
        if ( (*pl)->plot[i].y !=  NULL)mess_free((*pl)->plot[i].y );
    }
    if ((*pl)->plot != NULL)mess_free((*pl)->plot);
    if ((*pl)->title != NULL)mess_free((*pl)->title);
    // if ((*pl)->legendGnuPos != NULL)mess_free((*pl)->legendGnuPos);
    // if ((*pl)->legendTikzPos != NULL)mess_free((*pl)->legendTikzPos);
    mess_free(*pl);
    *pl = NULL;
    return 0;
}

/**
 * @brief Clear all data of a data element in a plot.
 * @param[in,out] p      plot with data
 * @param[in] data input number of data element
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref mess_plotExport_clearData function clears all data in a data element with number \f$ data \f$ in a plot \f$ p \f$.
 */
int mess_plotExport_clearData (mess_plotExport p, int data){
    MSG_FNAME(__func__);
    if ( p == NULL) {
        MSG_ERROR("p points to NULL\n");
        return MESS_ERROR_NULLPOINTER;
    }
    if ( data < 0 || data >= p->plot_len){
        MSG_ERROR("data is out of range\n");
        return MESS_ERROR_ARGUMENTS;
    }

    p->plot[data].len = 0;
    mess_free(  p->plot[data].x ) ;  p->plot[data].x = NULL;
    mess_free(  p->plot[data].y ) ;  p->plot[data].y = NULL;

    return 0;
}

/**
 * @brief Add new data values to a plot.
 * @param[in] p      input plot to add data
 * @param[in] data  input number of data
 * @param[in] x     input  \f$ x \f$-coordinate of data
 * @param[in] y     input  \f$ y \f$-coordinate of data
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref mess_plotExport_addData function adds a new value \f$ (x,y) \f$ to the data element with
 * the number \f$ data \f$ and updates the graphics.
 */
int mess_plotExport_addData(mess_plotExport p, int data, double x, double y){
    MSG_FNAME(__func__);
    mess_int_t i;
    if ( p == NULL){
        MSG_ERROR("p points to NULL\n");
        return MESS_ERROR_NULLPOINTER;
    }
    if ( data < 0 || data >= p->plot_len){
        MSG_ERROR("data is out of range\n");
        return MESS_ERROR_ARGUMENTS;
    }

    i = p->plot[data].len;

    p->plot[data].len++;
    mess_try_realloc(p->plot[data].x, double *, sizeof(double) * p->plot[data].len);
    mess_try_realloc(p->plot[data].y, double *, sizeof(double) * p->plot[data].len);
    p->plot[data].x[i] =x;
    p->plot[data].y[i] =y;

    return 0;
}

/**
 * @brief Set position of the legend.
 * @param[in,out] p      plot to add legend
 * @param[in] orientation input orientation of the key-box used in the diagramm
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref mess_plotExport_setLegendPos function sets the position or orientation of the legend used in the
 * diagramm. \n
 * Possible orientation values are
 * <ul>

 * <li> "belowCenter",
 * <li> "bottomRight",
 * <li> "bottomLeft",
 * <li> "bottom",
 * <li> "left",
 * <li> "right",
 * <li> "topRight",
 * <li> "topLeft",
 * <li> "top".
 * </ul>
 */
int mess_plotExport_setLegendPos(mess_plotExport p, char *orientation){
    MSG_FNAME(__func__);
    if ( p == NULL){
        MSG_ERROR("p points to NULL\n");
        return MESS_ERROR_NULLPOINTER;
    }

    if(strcmp(orientation,"topRight")==0 || strcmp(orientation, "tr")==0){
        p->legendTikzPos = ", legend pos=north east";
        p->legendGnuPos = "top right";
    }
    else if(strcmp(orientation,"topLeft")==0 || strcmp(orientation, "tl")==0){
        p->legendTikzPos = ", legend pos=north west";
        p->legendGnuPos = "top left";
    }
    else if(strcmp(orientation, "top")==0 || strcmp(orientation,"t")==0){
        p->legendTikzPos = ", legend pos=north";
        p->legendGnuPos = "top";
    }
    else if(strcmp(orientation, "bottomRight")==0 || strcmp(orientation,"br")==0){
        p->legendTikzPos = ", legend pos=south east";
        p->legendGnuPos = "bottom right";
    }
    else if(strcmp(orientation, "bottomLeft")==0 || strcmp(orientation,"bl")==0){
        p->legendTikzPos = ", legend pos=south west";
        p->legendGnuPos = "bottom left";
    }
    else if(strcmp(orientation,"bottom")==0 || strcmp(orientation, "b")==0){
        p->legendTikzPos = ", legend pos=south";
        p->legendGnuPos = "bottom";
    }
    else if(strcmp(orientation,"left")==0 || strcmp(orientation,"l")==0){
        p->legendTikzPos = ", legend pos=west";
        p->legendGnuPos = "left";
    }
    else if(strcmp(orientation,"right")==0 || strcmp(orientation,"r")==0){
        p->legendTikzPos = ", legend pos=east";
        p->legendGnuPos = "top left";
    }
    else if(strcmp(orientation,"belowCenter")==0){
        p->legendTikzPos = ", legend pos=outer south";
        p->legendGnuPos = "bmargin center horizontal";
    }
    else
        MSG_PRINT("unknown orientation. Legend Position is set to top right");
    return 0;
}

/**
 * @brief Set label of a data element.
 * @param[in,out] p plot
 * @param[in] data  input number of the data element
 * @param[in] label input the label
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref mess_plotExport_setLabel function sets the label of a data element with the number \f$ data \f$.\n
 * The label will be cut off after \f$ 40 \f$ characters.
 *
 */
int mess_plotExport_setLabel(mess_plotExport p, int data, char *label){
    MSG_FNAME(__func__);
    if ( p == NULL){
        MSG_ERROR("p points to NULL\n");
        return MESS_ERROR_NULLPOINTER;
    }
    if ( data < 0 || data >= p->plot_len){
        MSG_ERROR("data is out of range\n");
        return MESS_ERROR_ARGUMENTS;
    }

    strncpy(p->plot[data].label, label, 40);

    return 0;
}
/**
 * @brief Set the color of a plot.
 * @param[in,out] p   plot
 * @param[in] data     input number of the data element
 * @param[in] color   input colorstyle to be used for this dataset
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref mess_plotExport_setColor function sets the color of the plot. \n
 * Color need to be given as a character colourname. Possible characters are:
 * <ul>
 * <li> "blue"
 * <li> "green"
 * <li> "orange"
 * <li> "red"
 * <li> etc.
 * </ul>
 * For more colornames look up 'show palette colornames'.
 *
 */
int mess_plotExport_setColor(mess_plotExport p, int data, char *color){
    MSG_FNAME(__func__);
    if ( p == NULL){
        MSG_ERROR("p points to NULL\n");
        return MESS_ERROR_NULLPOINTER;
    }
    if ( data < 0 || data >= p->plot_len){
        MSG_ERROR("data is out of range\n");
        return MESS_ERROR_ARGUMENTS;
    }
    strncpy(p->plot[data].color, color, 40);

    return 0;
}
/**
 * @brief Set the style of a plot.
 * @param[in,out] p       plot
 * @param[in] data    input number of the data element
 * @param[in] type  input type of style to use for this dataset
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref mess_plotExport_setType function sets the type of the dataset. \n
 * Type need to be given as a character typename.Possible types are
 * <ul>
 * <li> "boxes",
 * <li> "dots",
 * <li> "impules",
 * <li> "lines" ,
 * <li> "points",
 * <li> "steps",
 * <li>etc.
 * </ul>
 *
 *
 */
int mess_plotExport_setType(mess_plotExport p, int data, char *type){
    MSG_FNAME(__func__);
    if ( p == NULL){
        MSG_ERROR("p points to NULL\n");
        return MESS_ERROR_NULLPOINTER;
    }
    if ( data < 0 || data >= p->plot_len){
        MSG_ERROR("data is out of range\n");
        return MESS_ERROR_ARGUMENTS;
    }
    strncpy(p->plot[data].type, type, 40);

    return 0;
}
/**
 * @brief Save all data to a .dat file.
 * @param[in] p     input plot
 * @param[in] path input destination folder
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref mess_plotExport_save function saves all data added to a .dat file
 *
 */
int mess_plotExport_save(mess_plotExport p, const char *path){

    MSG_FNAME(__func__);
    if ( p == NULL){
        MSG_ERROR("p points to NULL\n");
        return MESS_ERROR_NULLPOINTER;
    }
    if ( p->plot_len <= 0){
        MSG_ERROR("data is out of range\n");
        return MESS_ERROR_ARGUMENTS;
    }

    int i, j, k=0, temp=0;

    /* char filename[80] ="";
       strcat(filename, path);
       strcat(filename, "_data.dat");        */

    FILE *data;
    data = fopen(path, "w");
    int numData = p->plot_len;
    for(temp =0 ; temp <numData; temp++){
        if(p->plot[temp].len>k)
            k = p->plot[temp].len;
    }
    fprintf(data,"########################\n#       Data\n########################\n\n");
    for(j=0; j<k ;j++){
        for(i=0; i<numData ;i++){
            if(j<p->plot[i].len)
                fprintf(data,"%20lg %20lg \t", p->plot[i].x[j], p->plot[i].y[j]);
            else
                fprintf(data, "%20lg %20lg \t", p->plot[i].x[p->plot[i].len-1], p->plot[i].y[p->plot[i].len-1]);

        }
        fprintf(data,"\n");
    }
    fclose(data);

    return 0;
}

/**
 * @brief Create a GNU plot script file.
 * @param[in] p    input plot
 * @param[in] path input destination folder to store the GNU plot script
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref mess_plotExport_createGnuScript_int function will use all known data from the plot to create a
 * script file. \n
 * This script file will be saved as a _script.gnu file with the plots title as a prefix.
 *
 */
static int mess_plotExport_createGnuScript_int(mess_plotExport p, const char *path ){

    MSG_FNAME(__func__);
    if ( p == NULL){
        MSG_ERROR("p points to NULL\n");
        return MESS_ERROR_NULLPOINTER;
    }
    if ( p->plot_len <= 0){
        MSG_ERROR("data is out of range\n");
        return MESS_ERROR_ARGUMENTS;
    }

    int i = 0, j=1;
    char filename[2048];
    strncpy(filename, path, 2048);
    strncat(filename, "_script.gnu",2047);
    char dataname[2048];
    strncpy(dataname, path, 2048);
    strncat(dataname, "_data.dat",2047);
    // MSG_PRINT("Filename = %s Dataname = %s\n", filename, dataname);
    mess_plotExport_save(p, dataname);

    FILE *data;
    data = fopen(filename, "w");
    fprintf(data, "#------- GnuPlot_Script --------\n");
    fprintf(data,"#Titel \nset title '%s'\n", p->title);

    /**fprintf(data,"\n#Diagramrange\nset xrange [%f:%f]\n", p->minX, p->maxX);
      fprintf(data,"set yrange [%f:%f]\n", p->minY, p->maxY); **/
    fprintf(data,"\n#Diagramlabels\nset xlabel '%s'\n", p->xLabel);
    fprintf(data,"set ylabel '%s'\n", p->yLabel);

    fprintf(data, "set key %s",p->legendGnuPos);

    fprintf(data,"\n#Plotbefehle \nplot ");
    for(i=0 ; i < p->plot_len ; i++){
        fprintf(data, "'%s' using %i:%i", dataname, j, j+1);
        if(strcmp("",p->plot[i].label) != 0)
            fprintf(data, " title '%s'", p->plot[i].label);
        if(strcmp("",p->plot[i].type) != 0)
            fprintf(data, " with %s", p->plot[i].type);
        if(p->plot[i].color >= 0)
            fprintf(data, " lc rgb \"%s\"", p->plot[i].color);

        j=j+2;

        if(i+1 < p->plot_len)
            fprintf(data,", ");
    }
    //pause command
    fprintf(data, "\npause -1\n");

    fclose(data);

    return 0;
}

/**
 * @brief Create a GNU plot script file.
 * @param[in] p    input plot
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref mess_plotExport_createGnuScript function will use all known data from the plot to create a
 * script file. \n
 * This script file will be saved as a _script.gnu file with the plots title as a prefix in a folder with the same title.
 *
 *
 */
int mess_plotExport_createGnuScript(mess_plotExport p){
    return mess_plotExport_createGnuScript_int(p, p->title);
}

/**
 * @brief Create a GNU plot script file.
 * @param[in] p    input plot
 * @param[in] path input destination folder to store the GNU plot script
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref mess_plotExport_createGnuScript_filename function will use all known data from the plot to create a
 * script file. \n
 * This script file will be saved as a _script.gnu file with the plots title as a prefix.
 *
 */
int mess_plotExport_createGnuScript_filename(mess_plotExport p, const char *path){
    return mess_plotExport_createGnuScript_int(p, path);
}

/**
 * @brief Generate a TikZ-diagram and store it in a .tex-file.
 * @param[in] p     input plot with data
 * @param[in] path  input destination folder to store the TikZ Script
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref mess_plotExport_createTikzScript_int evaluates information stored in a plot and creates a
 * TikZ-picture. \n
 * This will be saved as a .tex-file and can be included to a TeX-document.
 */
static int mess_plotExport_createTikzScript_int(mess_plotExport p, const char *path){

    MSG_FNAME(__func__);
    if ( p == NULL){
        MSG_ERROR("p points to NULL\n");
        return MESS_ERROR_NULLPOINTER;
    }
    if ( p->plot_len <= 0){
        MSG_ERROR("data is out of range\n");
        return MESS_ERROR_ARGUMENTS;
    }

    int i=0, j=0;

    char filename[80] ="";
    strcat(filename, path);
    strcat(filename, "_tikz.tex");

    char dataname[80] ="";
    strcat(dataname, path);
    strcat(dataname, "_data.dat");

    mess_plotExport_save(p, path);

    FILE *data;
    data = fopen(filename, "w");
    fprintf(data, "%% --------- TikZ Plot --------\n%%\n%%----------------------------\n");
    fprintf(data, "\\begin{figure}\n\t");
    fprintf(data, "\\begin{tikzpicture}\n\t\t");

    if(p->xscale == 1 && p->yscale == 1)
        fprintf (data, "\\begin{loglogaxis}[xlabel=%s,ylabel=%s%s]\n\t\t\t", p->xLabel, p->yLabel, p->legendTikzPos);
    else if(p->xscale == 0 && p->yscale == 1)
        fprintf (data, "\\begin{semilogyaxis}[xlabel=%s,ylabel=%s%s]\n\t\t\t", p->xLabel, p->yLabel, p->legendTikzPos);
    else if(p->xscale == 1 && p->yscale == 0)
        fprintf (data, "\\begin{semilogxaxis}[xlabel=%s,ylabel=%s%s]\n\t\t\t", p->xLabel, p->yLabel, p->legendTikzPos);
    else
        fprintf (data, "\\begin{axis}[xlabel=%s,ylabel=%s%s]\n\t\t\t", p->xLabel, p->yLabel, p->legendTikzPos);


    for( i=0; i<p->plot_len; i++){
        if(p->plot[i].color>=0)
            fprintf(data, "\\addplot[%s] coordinates{", p->plot[i].color);
        else
            fprintf(data, "\\addplot coordinates {");
        for( j=0; j<p->plot[i].len; j++){
            fprintf(data, "\n\t\t\t\t(%f, \t%f)", p->plot[i].x[j], p->plot[i].y[j]);
        }
        fprintf(data, "\n\t\t\t};\n\n\t\t");
    }

    fprintf(data, "\\legend{");
    for( i=0;i<p->plot_len-1; i++){
        fprintf(data, "%s, ", p->plot[i].label);
    }

    fprintf(data, "%s}\n\n\t\t", p->plot[p->plot_len-1].label);

    if(p->xscale == 1 && p->yscale == 1)
        fprintf (data, "\\end{loglogaxis}\n\t");
    else if(p->xscale == 0 && p->yscale == 1)
        fprintf (data, "\\end{semilogyaxis}\n\t");
    else if(p->xscale == 1 && p->yscale == 0)
        fprintf (data, "\\end{semilogxaxis}\n\t");
    else
        fprintf (data, "\\end{axis}\n\t");

    fprintf(data, "\\end{tikzpicture}\n");
    fprintf(data, "\\end{figure}\n");
    fclose(data);

    return 0;
}

/**
 * @brief Generate a TikZ-diagram and store it in a .tex-file.
 * @param[in] p input plot with data
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref mess_plotExport_createTikzScript function evaluates the information stored in a plot and creates a
 * TikZ-picture. \n
 * This will be saved as a .tex-file and can be included to a TeX-document.
 */
int mess_plotExport_createTikzScript(mess_plotExport p){
    return mess_plotExport_createTikzScript_int(p, p->title);
}

/**
 * @brief Generate a TikZ-diagram and store it in a .tex-file.
 * @param[in] p input plot with data
 * @param[in] path input destination folder to store the TikZ Script
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref mess_plotExport_createTikzScript_filename function evaluates the information stored in a plot and creates a
 * TikZ-picture.\n
 * This will be saved as a .tex-file in the destination folder \f$ path \f$ and can be included to a TeX-document.
 */
int mess_plotExport_createTikzScript_filename(mess_plotExport p, const char *path){
    return mess_plotExport_createTikzScript_int(p, path);
}






