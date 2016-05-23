#ifndef LP_MP_ASCII_PLOT_VISITOR_HXX
#define LP_MP_ASCII_PLOT_VISITOR_HXX

#include "standard_visitor.hxx"
#include <ncurses.h>
//#include <form.h>

namespace LP_MP {

// visitor deriving from StandardVisitor and printing convergence plots via curses and using ascii art
// also templatize for visitor class to make tightening visitor possible
template<class PROBLEM_DECOMPOSITION>
class AsciiPlotVisitor : public StandardVisitor<PROBLEM_DECOMPOSITION> {
public:
   using BaseVisitor = StandardVisitor<PROBLEM_DECOMPOSITION>;
   AsciiPlotVisitor(TCLAP::CmdLine& cmd, PROBLEM_DECOMPOSITION& pd) 
      : BaseVisitor(cmd,pd)
   {}

   LPVisitorReturnType begin(const LP* lp)
   {
      auto ret = BaseVisitor::begin(lp);
      //const REAL lowerBound = StandardVisitor::GetDualBound();
      //const REAL upperBound = lp->BestPrimalBound(); // note: LP does not compute initial primal bound right now. Implement this
      std::cout << "note: LP does not compute initial primal bound right now. Implement this\n";
      return ret;
   }

   template<LPVisitorReturnType LP_STATE>
   LPVisitorReturnType visit(LP* lp)
   {
      auto ret_state = BaseVisitor::template visit<LP_STATE>(lp);
      //if(LP_STATE == LPVisitorReturnType::ReparametrizeAndComputePrimal) {
         //const REAL lowerBound = StandardVisitor::GetDualBound();
         //const REAL upperBound = lp->BestPrimalBound();

      //}
      return ret_state; 
   }

private:
   typedef double (*yfunction)(double x);
   struct _viewwin {
      double xmin, xmax;
      double ymin, ymax;
      double xscl, yscl;
   };
   typedef struct _viewwin viewwin;

   double estimateSlope(yfunction func, double x, double accuracy) const
   {
      double y1 = func(x - accuracy);
      double y2 = func(x + accuracy);
      return (y2 - y1) / (2 * accuracy);
   }
   double scale(double value, double omin, double omax, double nmin, double nmax)
   {
      /* Useful function to scale a value in one range to a different range.
         omin/omax - old range
         nmin/nmax - new range
         */
      double x = (value - omin) / (omax - omin);
      return x * (nmax - nmin) + nmin;
   }
   char slopeChar(double slope)
   {
      // Gets the character to display at a point in the graph with a given slope.
      double a = std::abs(slope);
      if (a < 0.5)        return '=';
      else if (a < 1.5)   return slope>0 ? '/' : '\\';
      else                return '|';
   }

   void plotPoint(WINDOW *win, const viewwin *view, double x, double y, char ch, int *scrY, int *scrX)
   {
      /* Displays a point on the screen at a location determined by graph coordinates.
         win       - ncurses window for drawing (can be NULL to only set scrY and scrX w/o drawing)
         view      - view parameters structure
         x/y       - graph coordinates for point
         ch        - character to display
         scrY/scrX - screen coordinates where point was drawn are saved here if not NULL
         */
      int xm, ym; getmaxyx(win, ym, xm);
      int xp = scale(x, view->xmin, view->xmax, 0, xm);
      int yp = scale(y, view->ymin, view->ymax, ym, 0);

      if (scrX) *scrX = xp;
      if (scrY) *scrY = yp;

      if (win) mvwaddch(win, yp, xp, ch);
   }


   void drawAxes(WINDOW *win, const viewwin *view)
   {
      // This function is what draws the axes on the screen.

      int xm, ym; getmaxyx(win, ym, xm);
      double x0 = scale(0, view->xmin, view->xmax, 0, xm);
      double y0 = scale(0, view->ymin, view->ymax, ym, 0);

      double xstep, ystep; getViewStep(win, view, &xstep, &ystep);

      int i; for (i=0; i<=xm; i++) {
         double plotx = view->xmin + xstep * i;
         int tick = fabs(fmod(plotx, view->xscl)) < xstep;
         mvwaddch(win, y0, i, tick ? '+':'-');
      }
      for (i=0; i<=ym; i++) {
         double ploty = view->ymin + ystep * i;
         int tick = fabs(fmod(ploty, view->yscl)) < ystep;
         mvwaddch(win, i, x0, tick ? '+':'|');
      }

      mvwaddch(win, y0, x0, '+');
   }

   void drawGraph(WINDOW *win, const viewwin *view, yfunction yfunc, int enableSlopeChars)
   {
      /* Draws a graph on the screen without axes.
         win              - ncurses window for drawing
         view             - view parameters structure
         yfunc            - function to graph (function pointer)
         enableSlopeChars - whether or not to call slopeChar to determine characters
         */
      int xm, ym; getmaxyx(win, ym, xm);
      double step; getViewStep(win, view, &step, NULL);
      double x; for (x = view->xmin; x <= view->xmax; x += step)
      {
         double y = yfunc(x);
         double d = estimateSlope(yfunc, x, step/2);
         plotPoint(win, view, x, y, enableSlopeChars ? slopeChar(d):'#', NULL, NULL);
      }
   }






   double xMin, xMax;
   double yMin, yMax;
};

} // end namespace LP_MP

#endif // LP_MP_ASCII_PLOT_VISITOR_HXX
