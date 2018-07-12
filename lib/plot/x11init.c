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
 * @file lib/plot/x11init.c
 * @brief A small plotter interface.
 * @author @koehlerm
 *
 * This file provides a small and easy plotting interface on top of the X11 library.
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

#ifdef MESS_HAVE_X11
#include <X11/Xlib.h>
#include <X11/Xutil.h>
#ifdef MESS_HAVE_XPM
#include <X11/xpm.h>
#endif

/** Displaypoints for the X and Y axis */
typedef struct _scalepoints{
    double val;
    char label[11];
} ScalePoints;

/**
 * @brief Map a point \f$ x \f$ from an interval \f$ [c,d] \f$ to a point in \f$ [a,b] \f$
 * @param[in] a input lower bound of destination interval
 * @param[in] b input upper bound of destination interval
 * @param[in] c input lower bound of source interval
 * @param[in] d input upper bound of source interval
 * @param[in] x input point to map in \f$ [a,b] \f$
 * @return mapped point in \f$ {a,b} \f$
 *
 * The @ref map_to_win function computes the linear mapping of a point \f$ x \f$ from an interval \f$ [c,d] \f$ to a
 * point in the interval \f$ [a,b] \f$ and prints out the mapped point.
 */
static double map_to_win( double a, double b, double c, double d, double x) {
    return ((b-a)/(d-c)*(x-c)+a);
}


/**
 * @brief Draw \f$ x \f$ and \f$ y \f$ axis on a plot.
 * @param[in,out] win    window to draw on
 * @param[in] xmin  input minimal value on \f$ x \f$ axis
 * @param[in] xmax  input maximal value on \f$ x \f$ axis
 * @param[in] ymin  input minimal value on \f$ y \f$ axis
 * @param[in] ymax  input maximal value on \f$ y \f$ axis
 * @param[in] XScale input points and labels on \f$ x \f$ axis
 * @param[in] XP    input number of labels on \f$ x \f$ axis
 * @param[in] YScale input points and labels on \f$ y \f$ axis
 * @param[in] YP    input number of labels on \f$ y \f$ axis
 * @param[in] axis_pos input position of axis (0 if they should go through (0,0) or 1 if they should be on the left and bottom )
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref drawaxis function draws the axis on a plot.
 *
 */
static int drawaxis(mess_plot win, double xmin, double xmax, double ymin, double ymax, ScalePoints *XScale,int XP,  ScalePoints *YScale, int YP , int axis_pos){
    int x0_pos=0, y0_pos=0;
    int x, y;
    int i;
    switch (axis_pos) {
        case 0:
            y0_pos = map_to_win(win->height-10, 10, ymin, ymax, 0);
            x0_pos = map_to_win( 10, win->width-10, xmin, xmax, 0);
            break;
        case 1:
            y0_pos = win->height-10;
            x0_pos = 10;
            break;
        default:
            return 1;
    }

    XDrawLine(win->display, win->win, win->gc, 10, y0_pos, win->width-10, y0_pos);
    XDrawLine(win->display, win->win, win->gc, x0_pos, 10, x0_pos, win->height-10);
    XDrawString(win->display, win->win, win->gc, x0_pos-8, y0_pos-8, "0", 1);
    for ( i = 0; i < XP; i++){
        x = map_to_win(10, win->width-10, xmin, xmax, XScale[i].val);
        XDrawLine(win->display, win->win, win->gc, x, y0_pos, x, y0_pos-10);
        XDrawString(win->display, win->win, win->gc, x-5, y0_pos+10, XScale[i].label, strlen(XScale[i].label));
    }
    for ( i = 0; i < YP; i++){
        y = map_to_win(win->height-10,10, ymin, ymax, YScale[i].val);
        XDrawLine(win->display, win->win, win->gc, x0_pos, y, x0_pos+10, y);
        XDrawString(win->display, win->win, win->gc, x0_pos+2, y-5, YScale[i].label, strlen(YScale[i].label));
    }

    XFlush(win->display);
    return 0;
}


/**
 * @brief Set drawing color.
 * @param[in,out] win   plotting window
 * @param[in] colorname input name of the color ( X11 color names)
 *
 * The @ref SetColor function sets the foreground color in a plot.\n
 * Allowed colornames are all that are used by X11.
 *
 */
static void SetColor(mess_plot win, char *colorname){
    XColor system_color;
    XColor exact_color;
    Colormap cmap;
    int rc;
    cmap = DefaultColormap (win->display,win->screen_num);
    /* allocate a "red" color map entry. */
    rc = XAllocNamedColor(win->display, cmap, colorname, &system_color, &exact_color);
    /* make sure the allocation succeeded. */
    if (rc == 0) {
        fprintf(stderr, "XAllocNamedColor - allocation of '%s' color failed. Using Black\n", colorname);
        XSetForeground(win->display, win->gc, BlackPixel(win->display, win->screen_num));
    }
    else {
        XSetForeground(win->display, win->gc, system_color.pixel);
        //  XFreeColors(win->display, cmap, &system_color.pixel, 1, 0);
    }

}



/**
 * @brief Draw a complete plot.
 * @param[in,out] p plot to draw
 *
 * The @ref plot function is the internal main function of the plotter.
 *
 */
static void plot(mess_plot p){
    MSG_FNAME(__func__);
    double xmax = 1, xmin = 1;
    double ymax = 1, ymin = 1;
    int x0maped = 0, x1maped = 0;
    int y0maped = 0, y1maped = 0;
    int plots;
    int i;
    int nXSP = 0;
    int nYSP = 0;
    ScalePoints *XSP = NULL;
    ScalePoints *YSP = NULL;
    if ( p->plot_len == 0) {
        MSG_INFO("nothing to plot\n");
        return;
    }
    XClearWindow(p->display, p->win);
    pthread_mutex_lock(&(p->lock));
    while (p->editing == 1)
        pthread_cond_wait(&(p->sig),&(p->lock));

    for ( plots = 0; plots < p->plot_len; plots++ ) {
        if ( plots == 0 && p->plot[0].len > 0) {
            xmin = xmax = p->plot[0].x[0];
            ymin = ymax = p->plot[0].y[0];
        }
        for ( i = 0; i < p->plot[plots].len; i++) {
            if ( xmax < p->plot[plots].x[i] ) xmax=p->plot[plots].x[i];
            if ( xmin > p->plot[plots].x[i] ) xmin=p->plot[plots].x[i];
            if ( ymax < p->plot[plots].y[i] ) ymax=p->plot[plots].y[i];
            if ( ymin > p->plot[plots].y[i] ) ymin=p->plot[plots].y[i];
        }
    }

    if ( p->xscale == MESS_PLOT_LOG) {
        if ( xmin != 0.0)
            xmin = floor(log ( xmin  )/ log(10));
        else
            xmin = -16;

        xmax = ceil( log ( xmax ) / log (10));
        nXSP = (int)(fabs(xmax - xmin)) +1 ;
        mess_try_alloc2(    XSP,  ScalePoints *,  sizeof(ScalePoints) * (nXSP));
        for ( i = 0; i < nXSP; i++){
            XSP[i].val=xmin+i;
            sprintf(XSP[i].label, "%.1e", pow(10,xmin+i));
        }
    } else if ( p->xscale == MESS_PLOT_LIN) {
        double ten = floor(log(fabs(xmax-xmin))/log(10));
        double step = pow ( 10.0, ten);
        double pos;
        step = (step > 0) ? step : 1;
        nXSP = 0;
        for ( pos = xmin; pos <=xmax; pos+=step) nXSP++;
        nXSP++;
        mess_try_alloc2(XSP, ScalePoints *, sizeof(ScalePoints) * (nXSP));
        for ( i = 0; i < nXSP; i++){
            XSP[i].val=xmin+i*step;
            sprintf(XSP[i].label, "%.1e", xmin+i*step);
        }

    }
    if ( p->yscale == MESS_PLOT_LOG) {
        if ( ymin != 0.0)
            ymin = floor(log ( ymin  )/ log(10));
        else
            ymin = -16;

        ymax = ceil( log ( ymax ) / log (10));
        nYSP = (int)(fabs(ymax - ymin)) +1 ;
        mess_try_alloc2(YSP, ScalePoints *, sizeof(ScalePoints) * (nYSP));
        for ( i = 0; i < nYSP; i++){
            YSP[i].val=ymin+i;
            sprintf(YSP[i].label, "%.1e", pow(10,ymin+i));
        }
    } else if ( p->yscale == MESS_PLOT_LIN) {
        double ten = floor(log(fabs(ymax-ymin))/log(10));
        double step = pow ( 10.0, ten);
        double pos;
        step = (step > 0) ? step : 1;
        nYSP = 0;
        for ( pos = ymin; pos <=ymax; pos+=step) nYSP++;
        nYSP++;
        mess_try_alloc2(YSP, ScalePoints *, sizeof(ScalePoints) * (nYSP));
        for ( i = 0; i < nYSP; i++){
            YSP[i].val=ymin+i*step;
            sprintf(YSP[i].label, "%.1e", ymin+i*step);
        }
    }

    XSetForeground(p->display, p->gc, BlackPixel(p->display, p->screen_num));
    drawaxis(p, xmin, xmax, ymin, ymax, XSP, nXSP, YSP, nYSP, 1);
    for ( plots = 0; plots < p->plot_len; plots++ ) {
        if ( strlen(p->plot[plots].color)!=0) {
            SetColor(p, p->plot[plots].color);
        }
        XDrawString(p->display, p->win, p->gc, 60, 20+plots*15, p->plot[plots].label, strlen(p->plot[plots].label));
        if ( p->plot[plots].len == 0) continue;
        if ( p->xscale == MESS_PLOT_LOG) {
            x0maped = map_to_win(10, p->width-10, xmin, xmax, log(p->plot[plots].x[0])/log(10));
        } else if (p->xscale == MESS_PLOT_LIN) {
            x0maped = map_to_win(10, p->width-10, xmin, xmax, p->plot[plots].x[0]);
        }
        if ( p->yscale == MESS_PLOT_LOG) {
            y0maped = map_to_win(p->height-10, 10, ymin, ymax, log(p->plot[plots].y[0])/log(10));
        } else if (p->yscale == MESS_PLOT_LIN) {
            y0maped = map_to_win(p->height-10, 10, ymin, ymax,p->plot[plots].y[0]);
        }

        for ( i = 0; i < p->plot[plots].len-1; i++) {
            if ( p->xscale == MESS_PLOT_LOG) {
                if (p->plot[plots].x[i+1] != 0.0)
                    x1maped = map_to_win(10, p->width-10, xmin, xmax, log(p->plot[plots].x[i+1])/log(10));
            } else if (p->xscale == MESS_PLOT_LIN) {
                x1maped = map_to_win(10, p->width-10, xmin, xmax, p->plot[plots].x[i+1]);
            }
            if ( p->yscale == MESS_PLOT_LOG) {
                if (p->plot[plots].y[i+1] != 0)
                    y1maped = map_to_win(p->height-10, 10, ymin, ymax, log(p->plot[plots].y[i+1])/log(10));
            } else if (p->yscale == MESS_PLOT_LIN) {
                y1maped = map_to_win(p->height-10, 10, ymin, ymax,p->plot[plots].y[i+1]);
            }
            XDrawLine(p->display, p->win, p->gc, x0maped, y0maped, x1maped, y1maped);
            x0maped=x1maped;
            y0maped=y1maped;
        }
    }

    mess_free(XSP);
    mess_free(YSP);
    XFlush(p->display);
    p->editing = 0;
    pthread_mutex_unlock(&(p->lock));
    pthread_cond_broadcast(&(p->sig));
}


/**
 * @brief Mainloop of the plot-window thread.
 * @param[in] arg input plot to work on
 *
 * The @ref __plotter function is the mainloop and the event handler of
 * the plotting window.
 *
 */
static void * __plotter(void *arg) {
    //  MSG_FNAME(__func__);
    mess_plot p = (mess_plot) arg;
    XEvent an_event;
    while (!(p->done)) {
        if ( XPending(p->display) == 0) {
            struct timespec req;
            req.tv_sec = 0;
            req.tv_nsec = 20000*1000;
            nanosleep(&req, &req);
            continue;
        }
        XNextEvent(p->display, &an_event);
        switch (an_event.type) {
            case ClientMessage:
                //  MSG_INFO("client message event\n");
                //   if (an_event.xclient.data.l[0] == p->delete_atom)
                //      break;
                break;
            case ConfigureNotify:
                //MSG_INFO("Configure Notify\n");
                p->width = an_event.xconfigure.width;
                p->height = an_event.xconfigure.height;
            case Expose:
                plot(p);
                break;
            case ButtonPress:
                break;
            case KeyPress: /* exit the application by braking out of the events loop. */
                //p->done = 1;
                break;
            case DestroyNotify:
                // printf("destroy event\n");
                p->done = 1;
                break;
            default: /* ignore any other event types. */
                break;
        } /* end switch on event type */
    } /* end while events handling */
    XFlush(p->display);
    XFreeGC(p->display, p->gc);
    XCloseDisplay(p->display);
    return NULL;
}


/**
 * @brief Create a plotting window.
 * @param[in,out] pl pointer to plot structure
 * @param[in] width input width of plot
 * @param[in] height input height of plot
 * @param[in] title input title of plotting window
 * @param[in] xscale input type of x scaling
 * @param[in] yscale input type of y scaling
 * @param[in] data input number of different data
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref mess_plot_create function creates a new plotting window of the size \f$ height \times width \f$ with \f$ title \f$
 * as title.\n
 * The \f$ xscale \f$ and \f$ yscale \f$ arguments can be
 * <ul>
 * <li> MESS_PLOT_LIN,
 * <li> MESS_PLOT_LOG
 * </ul>
 * to specify if the scaling is linear or logarithmic.\n
 * The \f$ data \f$ arguments specify the number of different data in the plot.
 *
 */
int mess_plot_create(mess_plot *pl, mess_int_t width, mess_int_t height, char *title, int xscale, int yscale, int data){
    MSG_FNAME(__func__);
    int i;
    mess_plot p;
    char *display_name = getenv("DISPLAY");  /* address of the X display.      */
    int win_border_width = 2;
    int screen_width;
    int screen_height;
    unsigned long valuemask = 0;        /* which values in 'values' to  */
    XGCValues values;           /* initial values for the GC.   */
    unsigned int line_width = 1;        /* line width for the GC.       */
    int line_style = LineSolid;     /* style for lines drawing and  */
    int cap_style = CapButt;        /* style of the line's edje and */
    int join_style = JoinBevel;     /*  joined lines.       */

    mess_try_alloc(*pl, mess_plot, sizeof(struct __mess_plotter));
    (*pl)->plot_len = 0;
    (*pl)->done = 0;
    (*pl)->height = height;
    (*pl)->width = width;
    (*pl)->xscale = xscale;
    (*pl)->yscale = yscale;
    mess_try_alloc((*pl)->title, char *, (strlen(title)+1) * sizeof ( char ));
    (*pl)->plot_len = data;
    mess_try_alloc((*pl)->plot, mess_plotdata *, data*sizeof(mess_plotdata));
    (*pl)->display = NULL;
    for ( i = 0; i < data; i++) {
        (*pl)->plot[i].x = NULL;
        (*pl)->plot[i].y = NULL;
        (*pl)->plot[i].len = 0;
        strncpy((*pl)->plot[i].color, "", 40);
        strncpy((*pl)->plot[i].label, "", 40);

    }
    strncpy((*pl)->title, title, strlen(title));
    pthread_mutex_init (&((*pl)->lock),NULL);
    pthread_cond_init (&((*pl)->sig),NULL);
    p = *pl;
    p->editing = 0;


    p->display = XOpenDisplay(display_name);
    if (p->display == NULL) {
        MSG_ERROR("cannot connect to X server '%s'\n", display_name);
        return 1;
    }

    p->screen_num =  DefaultScreen(p->display);
    screen_width = DisplayWidth(p->display, p->screen_num);
    screen_height = DisplayHeight(p->display, p->screen_num);
    if ( p->width == 0 || p->height == 0) {
        p->width = screen_width / 3;
        p->height = screen_height / 3;
    }
    p->win = XCreateSimpleWindow(p->display, RootWindow(p->display, p->screen_num),
            0, 0, p->width, p->height, win_border_width,
            BlackPixel(p->display, p->screen_num),
            WhitePixel(p->display, p->screen_num));
    p->delete_atom = XInternAtom(p->display, "WM_DELETE_WINDOW", True);
    if (p->delete_atom)
        XSetWMProtocols(p->display, p->win, &(p->delete_atom), 1);

    XStoreName(p->display, p->win, title);
    XMapWindow(p->display, p->win);
    XFlush(p->display);

    // generate GC
    p->gc = XCreateGC(p->display, p->win, valuemask, &values);

    if (p->gc < 0) {
        fprintf(stderr, "XCreateGC: \n");
        return 1;
    }
    XSetForeground(p->display, p->gc, BlackPixel(p->display, p->screen_num));
    XSetBackground(p->display, p->gc, WhitePixel(p->display, p->screen_num));
    XSetLineAttributes(p->display, p->gc,line_width, line_style, cap_style, join_style);
    XSetFillStyle(p->display, p->gc, FillSolid);
    XSync(p->display, False);
    XFlush(p->display);
    XSelectInput(p->display, p->win, ExposureMask | KeyPressMask | ButtonPressMask | StructureNotifyMask);

    pthread_create(&(p->thread),NULL, __plotter, (void *)p);
    return 0;
}



/**
 * @brief Add new data values to a plot.
 * @param[in,out] p plot to add the data
 * @param[in] data input number of data
 * @param[in] x  input \f$ x \f$ coordinate of data
 * @param[in] y  input \f$ y \f$ coordinate of data
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref mess_plot_addData function adds a new value \f$ (x,y) \f$ to the data element with
 * number \f$ data \f$ and updates the graphics.
 */
int mess_plot_addData(mess_plot p, int data, double x, double y){
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

    pthread_mutex_lock(&(p->lock));
    while (p->editing == 1)
        pthread_cond_wait(&(p->sig),&(p->lock));
    i = p->plot[data].len;
    p->plot[data].len++;
    mess_try_realloc(p->plot[data].x, double *, sizeof(double) * p->plot[data].len);
    mess_try_realloc(p->plot[data].y, double *, sizeof(double) * p->plot[data].len);
    p->plot[data].x[i] =x;
    p->plot[data].y[i] =y;
    p->editing = 0;
    pthread_mutex_unlock(&(p->lock));
    pthread_cond_broadcast(&(p->sig));
    return 0;
}


/**
 * @brief Clear all data of a data element in a plot.
 * @param[in,out] p plot with data
 * @param[in] data input number of the data element
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref mess_plot_clearData function clears all data in the data element with the number \f$ data \f$.
 *
 */
int mess_plot_clearData (mess_plot p, int data) {
    MSG_FNAME(__func__);
    if ( p == NULL) {
        MSG_ERROR("p points to NULL\n");
        return MESS_ERROR_NULLPOINTER;
    }
    if ( data < 0 || data >= p->plot_len){
        MSG_ERROR("data is out of range\n");
        return MESS_ERROR_ARGUMENTS;
    }

    pthread_mutex_lock(&(p->lock));
    while (p->editing == 1)
        pthread_cond_wait(&(p->sig),&(p->lock));
    p->plot[data].len = 0;
    mess_free(  p->plot[data].x ) ;  p->plot[data].x = NULL;
    mess_free(  p->plot[data].y ) ;  p->plot[data].y = NULL;
    p->editing = 0;
    pthread_mutex_unlock(&(p->lock));
    pthread_cond_broadcast(&(p->sig));
    return 0;


}


/**
 * @brief Set label of a data element.
 * @param[in,out] p  plot
 * @param[in] data   input number of the data element
 * @param[in] label  input label
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref mess_plot_setLabel function sets the label of the data element with the number \f$ data \f$. \n
 * The label will be cut off after \f$ 40 \f$ characters.
 *
 */
int mess_plot_setLabel(mess_plot p, int data, char *label){
    MSG_FNAME(__func__);
    if ( p == NULL){
        MSG_ERROR("p points to NULL\n");
        return MESS_ERROR_NULLPOINTER;
    }
    if ( data < 0 || data >= p->plot_len){
        MSG_ERROR("data is out of range\n");
        return MESS_ERROR_ARGUMENTS;
    }

    pthread_mutex_lock(&(p->lock));
    while (p->editing == 1)
        pthread_cond_wait(&(p->sig),&(p->lock));
    strncpy(p->plot[data].label, label, 40);
    p->editing = 0;
    pthread_mutex_unlock(&(p->lock));
    pthread_cond_broadcast(&(p->sig));
    return 0;
}


/**
 * @brief Set color of a data element.
 * @param[in,out] p plot
 * @param[in] data  input number of the data element
 * @param[in] color input colorname
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref mess_plot_setColor function sets the color of the data element
 * with the number \f$ data \f$. \n
 * Colorname can be one of all names your X11 implementation allow.
 *
 */
int mess_plot_setColor(mess_plot p, int data, char *color){
    MSG_FNAME(__func__);
    if ( p == NULL){
        MSG_ERROR("p points to NULL\n");
        return MESS_ERROR_NULLPOINTER;
    }
    if ( data < 0 || data >= p->plot_len){
        MSG_ERROR("data is out of range\n");
        return MESS_ERROR_ARGUMENTS;
    }

    pthread_mutex_lock(&(p->lock));
    while (p->editing == 1)
        pthread_cond_wait(&(p->sig),&(p->lock));
    strncpy(p->plot[data].color,color, 40);
    p->editing = 0;
    pthread_mutex_unlock(&(p->lock));
    pthread_cond_broadcast(&(p->sig));
    return 0;
}


/**
 * @brief Update / redraw a plot.
 * @param[in,out] p plot
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref mess_plot_update function updates/redraws a given plot.
 *
 */
int mess_plot_update(mess_plot p) {
    MSG_FNAME(__func__);
    XExposeEvent event;
    if ( p == NULL){
        MSG_ERROR("p points to NULL\n");
        return MESS_ERROR_NULLPOINTER;
    }
    event.type = Expose;
    event.display = p->display;
    event.window = p->win;
    event.send_event = True;
    XSendEvent(p->display,p->win, True, 0,(XEvent *)&event);
    return 0;
}


/**
 * @brief Internal helper function to save a plot as xpm.
 * @param[in] display input pointer to the display
 * @param[in] screen          input screennumber (unused)
 * @param[in] win_info  input Window information
 * @param[in] pixmap    input pixmap to save
 * @param[in] width     input width of image
 * @param[in] height            input height of image
 * @param[in] depth     input color depth of image (unused)
 * @param[in] filename  input filename of the output image
 * @return zero on success input or a non-zero error value otherwise
 *
 * The @ref save_as_xpm_file function saves a given pixmap to a file. \n
 * It is only a internal helper function.
 *
 */
static int save_as_xpm_file (Display *display, int screen, XWindowAttributes *win_info, Pixmap pixmap, unsigned int width, unsigned int height, unsigned int depth, char *filename)
{
#ifdef MESS_HAVE_XPM
    XpmAttributes attributes;
    //    attributes.colormap = private_cmap;  /* win_info->colormap; */
    attributes.colormap = win_info->colormap;
    attributes.width = width;
    attributes.height = height;
    attributes.valuemask = XpmColormap | XpmSize;
    XpmWriteFileFromPixmap(display, filename,pixmap,(Pixmap)0,&attributes);
#endif
    return 0;
}


/**
 * @brief Save a plot to a xpm image.
 * @param[in] p input plot to save
 * @param[in] filename input filename of image
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref mess_plot_save function saves a given plot \f$ p \f$ to a file called \f$ filename \f$ in a XPM image format.
 *
 */
int mess_plot_save(mess_plot p, char *filename) {
    MSG_FNAME(__func__);
    Pixmap tosave;
    int junk_left, junk_top;
    unsigned int junk_width, junk_height, junk_border_width;
    Window junk_root;
    XWindowAttributes win_info;
    unsigned int depth;

    mess_check_nullpointer(p);
    mess_check_nullpointer(filename);

    XGetWindowAttributes (p->display, p->win, &win_info);
    /* Use the depth of `a_window' for the depth of the pixmap */
    XGetGeometry (p->display, p->win, &junk_root, &junk_left, &junk_top, &junk_width, &junk_height, &junk_border_width, &depth);

    tosave = XCreatePixmap(p->display, DefaultRootWindow(p->display), p->width, p->height, depth);

    /*  now copy the area we specified  */

    XCopyArea(p->display, p->win, tosave, p->gc, 0, 0,   p->width, p->height, 0, 0);

    save_as_xpm_file(p->display, p->screen_num, &win_info, tosave, p->width, p->height, depth, filename);
    XFreePixmap(p->display, tosave);
    return 0;
}


/**
 * @brief Close a plotting window and clean up all internal data.
 * @param[in,out] p pointer to plot
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref mess_plot_close function closes a plot and destroys all internal data.
 *
 */
int mess_plot_close(mess_plot *p) {
    MSG_FNAME(__func__);
    XExposeEvent event;
    int i;
    if ( p == NULL ){
        MSG_ERROR("p points to NULL\n");
        return MESS_ERROR_NULLPOINTER;
    }

    (*p)->done = 1;
    event.type = Expose;
    XSendEvent((*p)->display, (*p)->win, True,0L, (XEvent *)&event);

    pthread_join( (*p)->thread, NULL);
    mess_free((*p)->title);
    for (i = 0; i < (*p)->plot_len; i++){
        if ( (*p)->plot[i].x != NULL)mess_free((*p)->plot[i].x);
        if ( (*p)->plot[i].y != NULL)mess_free((*p)->plot[i].y);
    }
    mess_free((*p)->plot);
    mess_free((*p));
    *p = NULL;
    return 0;
}


// Dummy functions if we don't have a useable X11 implementation.
#else

/**
 * @brief Create a plotting window.
 * @param[in,out] pl pointer to plot structure
 * @param[in] width input width of plot
 * @param[in] height input height of plot
 * @param[in] title input title of plotting window
 * @param[in] xscale input type of x scaling
 * @param[in] yscale input type of y scaling
 * @param[in] data input number of different data
 * @return always zero
 *
 * The @ref mess_plot_create function creates a new plotting window of the size \f$ height \times width \f$ with \f$ title \f$
 * as title.\n
 * The \f$ xscale \f$ and \f$ yscale \f$ arguments can be
 * <ul>
 * <li> MESS_PLOT_LIN,
 * <li> MESS_PLOT_LOG
 * </ul>
 * to specify if the scaling is linear or logarithmic.\n
 * The \f$ data \f$ arguments specify the number of different data in the plot.
 *
 * @attention Only a function for a non available X11 implementation.
 *
 */
int mess_plot_create(mess_plot *pl, mess_int_t width, mess_int_t height, char *title, int xscale, int yscale, int data){
    return 0;
}

/**
 * @brief Add new data values to a plot.
 * @param[in,out] p plot to add the data
 * @param[in] data input number of data
 * @param[in] x input \f$ x \f$ coordinate of data
 * @param[in] y input \f$ y \f$ coordinate of data
 * @return always zero
 *
 * The @ref mess_plot_addData function adds a new value \f$ (x,y) \f$ to the data element with
 * number \f$ data \f$ and updates the graphics.
 *
 * @attention Only a function for a non available X11 implementation.
 *
 */
int mess_plot_addData(mess_plot p, int data, double x, double y){
    return 0;
}

/**
 * @brief Clear all data of a data element in a plot.
 * @param[in,out] p plot with data
 * @param[in] data input number of the data element
 * @return always zero
 *
 * The @ref mess_plot_clearData function clears all data in the data element with the number \f$ data \f$.
 *
 * @attention Only a function for a non available X11 implementation.
 *
 */
int mess_plot_clearData (mess_plot p, int data) {
    return 0;
}

/**
 * @brief Set label of a data element.
 * @param[in,out] p  plot
 * @param[in] data   input number of the data element
 * @param[in] label  input label
 * @return always zero
 *
 * The @ref mess_plot_setLabel function sets the label of the data element with the number \f$ data \f$. \n
 * The label will be cut off after \f$ 40 \f$ characters.
 *
 * @attention Only a function for a non available X11 implementation.
 *
 */
int mess_plot_setLabel(mess_plot p, int data, char *label){
    return 0;
}

/**
 * @brief Set color of a data element.
 * @param[in,out] p plot
 * @param[in] data  input number of the data element
 * @param[in] color input colorname
 * @return always zero
 *
 * The @ref mess_plot_setColor function sets the color of the data element
 * with the number \f$ data \f$. \n
 * Colorname can be one of all names your X11 implementation allow.
 *
 * @attention Only a function for a non available X11 implementation.
 *
 */
int mess_plot_setColor(mess_plot p, int data, char *color){
    return 0;
}

/**
 * @brief Update / redraw a plot.
 * @param[in,out] p plot
 * @return always zero
 *
 * The @ref mess_plot_update function updates/redraws a given plot.
 *
 * @attention Only a function for a non available X11 implementation.
 *
 */
int mess_plot_update(mess_plot p) {
    return 0;
}

/**
 * @brief Save a plot to a xpm image.
 * @param[in] p input plot to save
 * @param[in] filename input filename of image
 * @return always zero
 *
 * The @ref mess_plot_save function saves a given plot \f$ p \f$ to a file called \f$ filename \f$ in a XPM image format.
 *
 * @attention Only a function for a non available X11 implementation.
 *
 */
int mess_plot_save(mess_plot p, char *filename) {
    return 0;
}

/**
 * @brief Close a plotting window and clean up all internal data.
 * @param[in,out] p pointer to plot
 * @return always zero
 *
 * The @ref mess_plot_close function closes a plot and destroys all internal data.
 *
 * @attention Only a function for a non available X11 implementation.
 *
 */
int mess_plot_close(mess_plot *p){
    return 0;
}

#endif


